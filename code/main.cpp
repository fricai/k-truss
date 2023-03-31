#include <fcntl.h>
#include <mpi.h>
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <execution>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

#include "argparse.hpp"
#include "debug.hpp"
#include "dsu.hpp"
#include "graph.hpp"
#include "types.hpp"
#include "utils.hpp"

int mpi_world_size, mpi_rank;

graph_t::graph_t() : actors{2 * omp_get_max_threads()} {}

void graph_t::init_io() {
    create_mmap();
    read_header();
    assign_owners();
    read_owned_graph();
}

void graph_t::init() {
    k1 = args.startk;
    k2 = args.endk;
    init_triangles();
    local_init();
}

void graph_t::local_init() {
    bucket.resize(k2 + 1);
    bucket_iter.resize(n);
    dead_triangle.resize(n);

    for (auto u : owned_vertices) {
        bucket_iter[u].resize(dodg[u].size());
        dead_triangle[u].resize(dodg[u].size());

        rep(idx_v, 0u, dodg[u].size()) {
            const auto tau = std::min(supp[u][idx_v], k2);
            bucket[tau].push_front({u, idx_v});
            bucket_iter[u][idx_v] = bucket[tau].begin();

            dead_triangle[u][idx_v].assign(inc_tri[u][idx_v].size(), false);
        }
    }
}

void graph_t::compute_truss() {
    fin_bucket.resize(k2 + 1);

    std::list<edge_idx_t> cur;

    std::vector<std::vector<triangle_t>> disjoint_queries(mpi_world_size);
    std::vector<triangle_t> send_buffer, recv_buffer;
    std::vector<int> send_offsets, send_cnts;
    std::vector<int> recv_cnts(mpi_world_size),
        recv_offsets(mpi_world_size + 1);

    std::vector<int8_t> responses, answers;

    for (int k = k1; k <= k2; ++k) {
        while (true) {
            {
                bool local_stop = bucket[k - 1].empty(), global_stop;
                MPI_Allreduce(&local_stop, &global_stop, 1, MPI_CXX_BOOL,
                              MPI_LAND, MPI_COMM_WORLD);
                if (global_stop) break;
            }

            cur.splice(cur.end(), bucket[k - 1]);

            // ensures the triangle deletions are mutually exclusive

#define F(x) (supp[w][get_edge_idx(w, x)] < k)
            // we mark dead triangles
            for (auto [u, idx_v] : cur) {
                const auto v = dodg[u][idx_v];
                // u -> v edge

                rep(tri_idx_w, 0u, inc_tri[u][idx_v].size()) {
                    if (dead_triangle[u][idx_v][tri_idx_w]) continue;
                    const auto w = inc_tri[u][idx_v][tri_idx_w];
                    if (rnk[w] < rnk[u]) {
                        // w -> u -> v
                        if (owner[w] == mpi_rank)  // local
                            dead_triangle[u][idx_v][tri_idx_w] = F(u) or F(v);
                        else
                            disjoint_queries[owner[w]].push_back({{w, u, v}});
                    } else if (rnk[w] < rnk[v])
                        // u -> w -> v
                        dead_triangle[u][idx_v][tri_idx_w] =
                            supp[u][get_edge_idx(u, w)] < k;
                }
            }

            flatten_vector(disjoint_queries, send_buffer, send_cnts,
                           send_offsets);
            for (auto& v : disjoint_queries) v.clear();

            MPI_Alltoall(send_cnts.data(), 1, MPI_INT, recv_cnts.data(), 1,
                         MPI_INT, MPI_COMM_WORLD);
            rep(i, 0, mpi_world_size) {
                recv_offsets[i + 1] = recv_offsets[i] + recv_cnts[i];
            }

            recv_buffer.resize(recv_offsets[mpi_world_size]);
            MPI_Alltoallv(send_buffer.data(), send_cnts.data(),
                          send_offsets.data(), mpi_array_int_3,
                          recv_buffer.data(), recv_cnts.data(),
                          recv_offsets.data(), mpi_array_int_3, MPI_COMM_WORLD);

            responses.resize(recv_offsets[mpi_world_size]);
            rep(i, 0, recv_offsets[mpi_world_size]) {
                const auto [w, u, v] = recv_buffer[i];
                responses[i] = F(u) or F(v);
            }

            answers.resize(send_offsets[mpi_world_size]);
            MPI_Alltoallv(responses.data(), recv_cnts.data(),
                          recv_offsets.data(), MPI_INT8_T, answers.data(),
                          send_cnts.data(), send_offsets.data(), MPI_INT8_T,
                          MPI_COMM_WORLD);

            std::vector<int> ctr(mpi_world_size);
            for (auto [u, idx_v] : cur) {
                // u -> v edge
                rep(tri_idx_w, 0u, inc_tri[u][idx_v].size()) {
                    if (dead_triangle[u][idx_v][tri_idx_w]) continue;
                    const auto w = inc_tri[u][idx_v][tri_idx_w];
                    if (rnk[w] < rnk[u] and owner[w] != mpi_rank) {
                        dead_triangle[u][idx_v][tri_idx_w] =
                            answers[send_offsets[owner[w]] + ctr[owner[w]]++];
                    }
                }
            }
#undef F

            const auto update_triangle = [&](int p, int q, int r) {
                assert(owner[p] == mpi_rank);

                const auto idx_q = get_edge_idx(p, q);
                const auto idx_r = get_triangle_idx(p, idx_q, r);

                if (dead_triangle[p][idx_q][idx_r]) return;
                dead_triangle[p][idx_q][idx_r] = true;

                auto& supp_val = supp[p][idx_q];
                if (supp_val >= k) {
                    --supp_val;
                    if (supp_val < k2) {
                        bucket[supp_val].splice(bucket[supp_val].end(),
                                                bucket[supp_val + 1],
                                                bucket_iter[p][idx_q]);
                    }
                }
            };

            for (auto [u, idx_v] : cur) {
                const auto v = dodg[u][idx_v];
                rep(tri_idx_w, 0u, inc_tri[u][idx_v].size()) {
                    if (dead_triangle[u][idx_v][tri_idx_w]) continue;

                    // decrement edge {u, w} and {v, w}
                    const auto w = inc_tri[u][idx_v][tri_idx_w];
                    auto tmp = [&](int p, int q, int r) {
                        if (rnk[p] > rnk[q]) std::swap(p, q);  // p -> q
                        if (owner[p] == mpi_rank)
                            update_triangle(p, q, r);
                        else
                            disjoint_queries[owner[p]].push_back({{p, q, r}});
                    };
                    tmp(w, u, v);
                    tmp(w, v, u);
                }
            }

            flatten_vector(disjoint_queries, send_buffer, send_cnts,
                           send_offsets);
            for (auto& v : disjoint_queries) v.clear();

            MPI_Alltoall(send_cnts.data(), 1, MPI_INT, recv_cnts.data(), 1,
                         MPI_INT, MPI_COMM_WORLD);

            rep(i, 0, mpi_world_size) {
                recv_offsets[i + 1] = recv_offsets[i] + recv_cnts[i];
            }

            recv_buffer.resize(recv_offsets[mpi_world_size]);
            MPI_Alltoallv(send_buffer.data(), send_cnts.data(),
                          send_offsets.data(), mpi_array_int_3,
                          recv_buffer.data(), recv_cnts.data(),
                          recv_offsets.data(), mpi_array_int_3, MPI_COMM_WORLD);

            for (auto [p, q, r] : recv_buffer) update_triangle(p, q, r);

            // finalize these edges
            fin_bucket[k - 1].splice(fin_bucket[k - 1].end(), cur);
        }
    }

    fin_bucket[k2].splice(fin_bucket[k2].end(), bucket[k2]);
}

void graph_t::create_mmap() {
    int fd = open(args.inputpath.data(), O_RDONLY);

    struct stat statbuf;
    fstat(fd, &statbuf);
    file_len = statbuf.st_size;

    file_map =
        (uint32_t*)mmap(nullptr, file_len, PROT_READ, MAP_PRIVATE, fd, 0);
}

bool graph_t::edge_oracle(int u, int v) {
    if (owner[u] == mpi_rank)
        return local_edge_oracle(u, v);
    else
        return file_edge_oracle(u, v);
}
bool graph_t::local_edge_oracle(int u, int v) {
    return std::binary_search(dodg[u].begin(), dodg[u].end(), v);
}
bool graph_t::file_edge_oracle(int u, int v) {
    return std::binary_search(file_map + offset[u] + 2,
                              file_map + offset[u + 1], (uint32_t)v);
}

idx_t graph_t::get_edge_idx(vertex_t x, vertex_t y) {
    return std::lower_bound(dodg[x].begin(), dodg[x].end(), y) -
           dodg[x].begin();
}

idx_t graph_t::get_triangle_idx(vertex_t x, idx_t idx_y, vertex_t z) {
    return std::lower_bound(inc_tri[x][idx_y].begin(), inc_tri[x][idx_y].end(),
                            z) -
           inc_tri[x][idx_y].begin();
}

void graph_t::assign_owners() {
    using I = int64_t;
    std::vector<I> load(n);
    rep(u, 0, n) load[u] = (rnk[u] >> 32) * (rnk[u] >> 32);

    // do static partioning instead
    // gives you locality

    I tot = 0;
    rep(u, 0, n) tot += load[u];

    owner.resize(n);

    int l = 0, r = 0;
    rep(i, 0, mpi_world_size) {
        const auto target = tot / (mpi_world_size - i);
        I cur_sum = 0;
        while (cur_sum < target) {
            assert(r < n);
            cur_sum += load[r];
            ++r;
        }
        rep(j, l, r) owner[j] = i;

        if (i == mpi_rank) {
            owned_vertices.resize(r - l);
            std::iota(owned_vertices.begin(), owned_vertices.end(), l);
        }

        l = r;
        tot -= cur_sum;
    }

    owner_actor.assign(n, -1);
    actor_owned_vertices.resize(actors);
    for (auto& v : actor_owned_vertices)
        v.reserve(owned_vertices.size() / actors + 1);
    for (auto u : owned_vertices) {
        owner_actor[u] = u % actors;
        actor_owned_vertices[owner_actor[u]].push_back(u);
    }
}

void graph_t::read_header() {
    n = file_map[0];
    m = file_map[1];

    std::ifstream hf(args.headerpath, std::ios::binary);
    offset.resize(n + 1);
    hf.read((char*)offset.data(), 4 * n);
    hf.close();

    rep(i, 0, n) offset[i] /= 4;

    int sum = 0;
    rnk.resize(n);
    rep(i, 0, n - 1) {
        rnk[i] = offset[i + 1] - offset[i] - 2;
        sum += rnk[i];
    }
    rnk[n - 1] = 2 * m - sum;
    offset[n] = offset[n - 1] + rnk[n - 1] + 2;

    rep(i, 0, n) rnk[i] = rnk[i] << 32 | i;
}

void graph_t::read_owned_graph() {
    dodg.resize(n);

    for (auto u : owned_vertices) {
        const auto nbh = file_map + offset[u] + 2;

        const auto id = nbh[-2], deg = nbh[-1];
        assert(id == (uint32_t)u);
        assert(deg == (rnk[u] >> 32));

        rep(i, 0u, deg) {
            const auto v = nbh[i];
            if (rnk[u] < rnk[v]) dodg[u].push_back(v);
        }

        assert(std::is_sorted(dodg[u].begin(), dodg[u].end()));
    }
}

void graph_t::init_triangles() {
    /*
     * inc_tri[p][idx_q] stores the triangle (p, q, r)
     * if and only if in the triangle
     *    q
     *   /  \
     *  p -> r
     *
     *  No notion of triangle support needed in mintruss
     */

    supp.resize(n);
    inc_tri.resize(n);

#pragma omp parallel for schedule(dynamic)
    rep(actor, 0, actors) {
        for (auto p : actor_owned_vertices[actor]) {
            supp[p].assign(dodg[p].size(), 0);
            inc_tri[p].assign(dodg[p].size(), {});
        }
    }

    std::vector queued_triangles(mpi_world_size,
                                 std::vector<std::vector<triangle_t>>(actors));
    // std::vector<std::vector<triangle_t>> queued_triangles(mpi_world_size);

#pragma omp parallel for schedule(dynamic)
    rep(actor, 0, actors) {
        for (auto p : actor_owned_vertices[actor]) {
            rep(idx_q_, 0u, dodg[p].size()) {
                rep(idx_r_, idx_q_ + 1, dodg[p].size()) {
                    vertex_t q = dodg[p][idx_q_], r = dodg[p][idx_r_];
                    idx_t idx_q = idx_q_, idx_r = idx_r_;

                    if (rnk[q] > rnk[r]) {
                        std::swap(idx_q, idx_r);
                        std::swap(q, r);
                    }

                    // rnk[p] < rnk[q] < rnk[r]

                    // p -> q edge exists
                    // p -> r edge exists
                    // if q -> r edge exists

                    if (edge_oracle(q, r)) {
                        inc_tri[p][idx_q].push_back(r);
                        inc_tri[p][idx_r].push_back(q);

                        queued_triangles[owner[q]][actor].push_back(
                            {{q, r, p}});
                    }
                }
            }
        }
    }

    std::vector<int> send_cnt, send_offsets;
    std::vector<triangle_t> send_triangles;

    flatten_3d_vector(queued_triangles, send_triangles, send_cnt, send_offsets);
    clear_vector(queued_triangles);

    std::vector<int> recv_cnt(mpi_world_size), recv_offsets(mpi_world_size + 1);
    MPI_Alltoall(send_cnt.data(), 1, MPI_INT, recv_cnt.data(), 1, MPI_INT,
                 MPI_COMM_WORLD);
    rep(i, 0, mpi_world_size) {
        recv_offsets[i + 1] = recv_offsets[i] + recv_cnt[i];
    }

    std::vector<triangle_t> recv_triangles(recv_offsets[mpi_world_size]);
    MPI_Alltoallv(send_triangles.data(), send_cnt.data(), send_offsets.data(),
                  mpi_array_int_3, recv_triangles.data(), recv_cnt.data(),
                  recv_offsets.data(), mpi_array_int_3, MPI_COMM_WORLD);
    clear_vector(send_triangles);

    // Would want to parallelize this
    std::sort(recv_triangles.begin(), recv_triangles.end(),
              [&](const auto& a, const auto& b) {
                  return owner_actor[a[0]] < owner_actor[b[0]];
              });

    // Note the best way of doing things :/

    std::vector<int> actor_offsets(actors + 1, 0);
#pragma omp parallel for schedule(dynamic)
    rep(i, 0, actors) {
        // I want first element with owner > i
        actor_offsets[i + 1] =
            std::lower_bound(
                recv_triangles.begin(), recv_triangles.end(), i,
                [&](const auto& a, int x) { return owner_actor[a[0]] <= x; }) -
            recv_triangles.begin();
    }

#pragma omp parallel for schedule(dynamic)
    rep(i, 0, actors) {
        rep(j, actor_offsets[i], actor_offsets[i + 1]) {
            const auto [q, r, p] = recv_triangles[j];
            const auto idx_r = get_edge_idx(q, r);
            inc_tri[q][idx_r].push_back(p);
        }
    }

#pragma omp parallel for schedule(dynamic)
    rep(i, 0, actors) {
        for (auto p : actor_owned_vertices[i]) {
            rep(idx, 0u, inc_tri[p].size()) {
                std::sort(inc_tri[p][idx].begin(), inc_tri[p][idx].end());
                // inc_tri[p][idx].shrink_to_fit();
                /*
                 * Should I keep this or not?
                 */

                supp[p][idx] = (int)inc_tri[p][idx].size();
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char** argv) {
    int required = MPI_THREAD_FUNNELED, provided;
    MPI_Init_thread(&argc, &argv, required, &provided);

    if (required != provided) {
        std::cout << "Threads not properly supported\n";
        MPI_Finalize();
        return 0;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    init_mpi_datatypes();

    args.taskid = 1;
    args.verbose = false;

    argp_parse(&argp, argc, argv, 0, 0, &args);

    const int k1 = args.startk;
    const int k2 = args.endk;
    assert(k1 <= k2);

    graph_t g;

    g.init_io();
    g.init();

    g.compute_truss();

    {
        std::vector<std::array<int, 3>> res;
        rep(k, 0, k2 + 1) {
            for (auto [u, idx_v] : g.fin_bucket[k]) {
                auto v = g.dodg[u][idx_v];
                if (u > v) std::swap(u, v);
                res.push_back({{u, v, k}});
            }
        }

        const int local_cnt = (int)res.size();

        std::vector<int> cnts(mpi_world_size);
        MPI_Gather(&local_cnt, 1, MPI_INT, cnts.data(), 1, MPI_INT, 0,
                   MPI_COMM_WORLD);

        std::vector<int> offsets(mpi_world_size + 1);
        rep(i, 0, mpi_world_size) offsets[i + 1] = offsets[i] + cnts[i];

        std::vector<std::array<int, 3>> collate(offsets[mpi_world_size]);
        MPI_Gatherv(res.data(), local_cnt, mpi_array_int_3, collate.data(),
                    cnts.data(), offsets.data(), mpi_array_int_3, 0,
                    MPI_COMM_WORLD);
        std::sort(collate.begin(), collate.end());

        std::cin.tie(nullptr)->sync_with_stdio(false);
        for (auto [u, v, k] : collate)
            std::cout << u << ' ' << v << ' ' << k << '\n';
    }

    MPI_Finalize();
}
