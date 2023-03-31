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
    bucket.resize(actors);
    rep(i, 0, actors) bucket[i].resize(k2 + 1);

    bucket_iter.resize(n);
    dead_triangle.resize(n);

#pragma omp parallel for schedule(dynamic, 1)
    rep(i, 0, actors) {
        for (auto u : actor_owned_vertices[i]) {
            bucket_iter[u].resize(dodg[u].size());
            dead_triangle[u].resize(dodg[u].size());

            rep(idx_v, 0u, dodg[u].size()) {
                const auto tau = std::min(supp[u][idx_v], k2);
                bucket[i][tau].push_front({u, idx_v});
                bucket_iter[u][idx_v] = bucket[i][tau].begin();

                dead_triangle[u][idx_v].assign(inc_tri[u][idx_v].size(), false);
            }
        }
    }
}

void graph_t::compute_truss() {
    fin_bucket.resize(k2 + 1);

    std::vector<std::list<edge_idx_t>> cur(actors);

    std::vector<triangle_t> send_buffer, recv_buffer;
    std::vector<int> send_offsets, send_cnts;
    std::vector<int> recv_cnts(mpi_world_size),
        recv_offsets(mpi_world_size + 1);

    std::vector disjoint_queries(mpi_world_size,
                                 std::vector<std::vector<triangle_t>>(actors));

    std::vector<std::vector<int>> cnt2, off2;

    std::vector<int8_t> responses, answers;

    for (int k = k1; k <= k2; ++k) {
        while (true) {
            {
                bool local_stop = true, global_stop;
                rep(i, 0, actors) local_stop &= bucket[i][k - 1].empty();
                MPI_Allreduce(&local_stop, &global_stop, 1, MPI_CXX_BOOL,
                              MPI_LAND, MPI_COMM_WORLD);
                if (global_stop) break;
            }

            rep(i, 0, actors) cur[i].splice(cur[i].end(), bucket[i][k - 1]);

            // ensures the triangle deletions are mutually exclusive

#define F(x) (supp[w][get_edge_idx(w, x)] < k)
            // we mark dead triangles
#pragma omp parallel for schedule(dynamic, 1)
            rep(i, 0, actors) {
                for (auto [u, idx_v] : cur[i]) {
                    const auto v = dodg[u][idx_v];
                    // u -> v edge

                    rep(tri_idx_w, 0u, inc_tri[u][idx_v].size()) {
                        if (dead_triangle[u][idx_v][tri_idx_w]) continue;
                        const auto w = inc_tri[u][idx_v][tri_idx_w];
                        if (rnk[w] < rnk[u]) {
                            // w -> u -> v
                            if (owner[w] == mpi_rank)  // local
                                dead_triangle[u][idx_v][tri_idx_w] =
                                    F(u) or F(v);
                            else
                                disjoint_queries[owner[w]][i].push_back(
                                    {{w, u, v}});
                        } else if (rnk[w] < rnk[v])
                            // u -> w -> v
                            dead_triangle[u][idx_v][tri_idx_w] =
                                supp[u][get_edge_idx(u, w)] < k;
                    }
                }
            }

            flatten_3d_vector(disjoint_queries, send_buffer, send_cnts,
                              send_offsets, cnt2, off2);
            for (auto& w : disjoint_queries)
                for (auto& v : w) v.clear();

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
#pragma omp parallel for
            rep(i, 0, recv_offsets[mpi_world_size]) {
                const auto [w, u, v] = recv_buffer[i];
                responses[i] = F(u) or F(v);
            }

            answers.resize(send_offsets[mpi_world_size]);
            MPI_Alltoallv(responses.data(), recv_cnts.data(),
                          recv_offsets.data(), MPI_INT8_T, answers.data(),
                          send_cnts.data(), send_offsets.data(), MPI_INT8_T,
                          MPI_COMM_WORLD);

#pragma omp parallel for schedule(dynamic, 1)
            rep(i, 0, actors) {
                std::vector ctr(mpi_world_size, 0);
                for (auto [u, idx_v] : cur[i]) {
                    // u -> v edge
                    rep(tri_idx_w, 0u, inc_tri[u][idx_v].size()) {
                        if (dead_triangle[u][idx_v][tri_idx_w]) continue;
                        const auto w = inc_tri[u][idx_v][tri_idx_w];
                        if (rnk[w] < rnk[u] and owner[w] != mpi_rank) {
                            dead_triangle[u][idx_v][tri_idx_w] =
                                answers[send_offsets[owner[w]] +
                                        off2[owner[w]][i] + ctr[owner[w]]++];
                        }
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
                        auto& b = bucket[owner_actor[p]];
                        b[supp_val].splice(b[supp_val].end(), b[supp_val + 1],
                                           bucket_iter[p][idx_q]);
                    }
                }
            };

#pragma omp parallel for schedule(dynamic, 1)
            rep(i, 0, actors) {
                for (auto [u, idx_v] : cur[i]) {
                    const auto v = dodg[u][idx_v];
                    rep(tri_idx_w, 0u, inc_tri[u][idx_v].size()) {
                        if (dead_triangle[u][idx_v][tri_idx_w]) continue;

                        // decrement edge {u, w} and {v, w}
                        const auto w = inc_tri[u][idx_v][tri_idx_w];
                        auto tmp = [&](int p, int q, int r) {
                            if (rnk[p] > rnk[q]) std::swap(p, q);  // p -> q
                            if (owner[p] == mpi_rank and owner_actor[p] == i)
                                update_triangle(p, q, r);
                            else
                                disjoint_queries[owner[p]][i].push_back(
                                    {{p, q, r}});
                        };
                        tmp(w, u, v);
                        tmp(w, v, u);
                    }
                }
            }

            flatten_3d_vector(disjoint_queries, send_buffer, send_cnts,
                              send_offsets);
            for (auto& w : disjoint_queries)
                for (auto& v : w) v.clear();

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

            // Can I parallelize this?
            std::sort(recv_buffer.begin(), recv_buffer.end(),
                      [&](const auto& a, const auto& b) {
                          return owner_actor[a[0]] < owner_actor[b[0]];
                      });

            std::vector<int> actor_offsets(actors + 1, 0);
#pragma omp parallel for schedule(dynamic, 1)
            rep(i, 0, actors) {
                // I want first element with owner > i
                actor_offsets[i + 1] =
                    std::lower_bound(recv_buffer.begin(), recv_buffer.end(), i,
                                     [&](const auto& a, int x) {
                                         return owner_actor[a[0]] <= x;
                                     }) -
                    recv_buffer.begin();
            }

#pragma omp parallel for schedule(dynamic, 1)
            rep(i, 0, actors) {
                rep(j, actor_offsets[i], actor_offsets[i + 1]) {
                    const auto [p, q, r] = recv_buffer[j];
                    update_triangle(p, q, r);
                }
            }

            // finalize these edges
            // DONT parallelize
            rep(i, 0, actors) fin_bucket[k - 1].splice(fin_bucket[k - 1].end(),
                                                       cur[i]);
        }
    }

    // DONT parallelize
    rep(i, 0, actors) fin_bucket[k2].splice(fin_bucket[k2].end(),
                                            bucket[i][k2]);
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

#pragma omp parallel for schedule(dynamic, 1)
    rep(actor, 0, actors) {
        for (auto p : actor_owned_vertices[actor]) {
            supp[p].assign(dodg[p].size(), 0);
            inc_tri[p].assign(dodg[p].size(), {});
        }
    }

    std::vector queued_triangles(mpi_world_size,
                                 std::vector<std::vector<triangle_t>>(actors));
    // std::vector<std::vector<triangle_t>> queued_triangles(mpi_world_size);

#pragma omp parallel for schedule(dynamic, 1)
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
#pragma omp parallel for schedule(dynamic, 1)
    rep(i, 0, actors) {
        // I want first element with owner > i
        actor_offsets[i + 1] =
            std::lower_bound(
                recv_triangles.begin(), recv_triangles.end(), i,
                [&](const auto& a, int x) { return owner_actor[a[0]] <= x; }) -
            recv_triangles.begin();
    }

#pragma omp parallel for schedule(dynamic, 1)
    rep(i, 0, actors) {
        rep(j, actor_offsets[i], actor_offsets[i + 1]) {
            const auto [q, r, p] = recv_triangles[j];
            const auto idx_r = get_edge_idx(q, r);
            inc_tri[q][idx_r].push_back(p);
        }
    }

#pragma omp parallel for schedule(dynamic, 1)
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

void output_2(const graph_t& g, int k1, int k2) {
    const int n = g.n;

    dsu F(n);

    std::vector<edge_t> queued_edges;

    std::vector<edge_t> spanning_forest;

    for (auto [u, idx_v] : g.fin_bucket[k2]) {
        auto v = g.dodg[u][idx_v];

        if (F.merge(u, v)) {
            // send to 0 to add to 0's DSU
            (mpi_rank == 0 ? spanning_forest : queued_edges)
                .push_back(edge_t({u, v}));
        }
    }

    std::vector<edge_t> edge_buffer;

    std::vector<int> queue_lengths(mpi_rank == 0 ? mpi_world_size : 0);

    int queued_edges_size = queued_edges.size();
    MPI_Gather(&queued_edges_size, 1, MPI_INT, queue_lengths.data(), 1, MPI_INT,
               0, MPI_COMM_WORLD);

    std::vector<int> displacements(mpi_rank == 0 ? mpi_world_size + 1 : 0);

    if (mpi_rank == 0) {
        displacements[0] = 0;
        rep(i, 0, mpi_world_size) displacements[i + 1] =
            displacements[i] + queue_lengths[i];
        edge_buffer.resize(displacements[mpi_world_size]);
    }

    MPI_Gatherv(queued_edges.data(), queued_edges_size, mpi_edge_t,
                edge_buffer.data(), queue_lengths.data(), displacements.data(),
                mpi_edge_t, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        for (auto [u, v] : edge_buffer) {
            if (F.merge(u, v)) {
                spanning_forest.push_back(edge_t({u, v}));
            }
        }
    }

    {
        int spanning_forest_len = (int)spanning_forest.size();
        MPI_Bcast(&spanning_forest_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        spanning_forest.resize(spanning_forest_len);
        MPI_Bcast(spanning_forest.data(), spanning_forest_len, MPI_INT, 0,
                  MPI_COMM_WORLD);
        if (mpi_rank != 0)
            for (auto [u, v] : spanning_forest) F.merge(u, v);
    }

    std::vector<std::string> grps(n);
    rep(u, 0, n) {
        if (F.size(u) > 1) grps[F.head(u)] += std::to_string(u) + " ";
    }

    std::vector<bool> vis(n, false);

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, args.outputpath.data(),
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    std::string buf;

    int fst = (n + mpi_world_size - 1) / mpi_world_size * mpi_rank,
        lst = std::min((n + mpi_world_size - 1) / mpi_world_size * (mpi_rank + 1), n);
    std::cout << "rank = " << mpi_rank << ", fst = " << fst << ", lst = " << lst << "\n";
    assert(fst <= n);

    std::vector<std::vector<int>> heads(lst-fst);
    auto proxy = [&](int v, int u) {
        u = F.head(u);
        if (grps[u].empty() or vis[u]) return;
        vis[u] = true;
        heads[v].push_back(u);
    };
    
    int n_infls = 0;
    rep(u, fst, lst) {
        int i = u-fst;
        uint32_t const* const neighbourhood = g.file_map + g.offset[u] + 2;
        const int deg = g.offset[u + 1] - g.offset[u] - 2;

        rep(j, 0, deg) proxy(i, neighbourhood[j]);
        
        for (auto v : heads[i]) vis[v] = false;

        if ((int)heads[i].size() < args.p) continue;

        n_infls++;
    }

    int global_n_infls;
    MPI_Allreduce(&n_infls, &global_n_infls, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        std::string global_n_infls_str = std::to_string(global_n_infls) + "\n";
        MPI_File_write_shared(fh, global_n_infls_str.data(), 
                global_n_infls_str.size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    rep(u, fst, lst) {
        int i = u-fst;
        if (heads[i].size() < args.p) continue;
        buf.clear();
        buf += std::to_string(u) + " ";
        if (args.verbose) {
            buf += "\n";
            for (auto r : heads[i]) buf += grps[r];
            buf += "\n";
        }
        MPI_File_write_shared(fh, buf.data(), buf.size(), MPI_CHAR,
                              MPI_STATUS_IGNORE);
    }

    MPI_File_close(&fh);
}

void output_1(const graph_t& g, int k1, int k2, int largest_truss) {
    const int n = g.n;
    std::ofstream outf(args.outputpath);
    if (args.verbose) {
        dsu F(n);

        std::vector<std::vector<std::vector<vertex_t>>> trusses(
            mpi_rank == 0 ? std::max(0, std::min(k2, largest_truss) - k1 + 1)
                          : 0);

        std::vector<std::vector<vertex_t>> scratch(n);
        std::vector<edge_t> queued_edges;

        std::vector<edge_t> edge_buffer;
        for (int k = largest_truss; k >= k1; --k) {
            for (auto [u, idx_v] : g.fin_bucket[k]) {
                auto v = g.dodg[u][idx_v];
                bool result = F.merge(u, v);
                if (result and mpi_rank != 0) {
                    // send to 0 to add to 0's DSU
                    queued_edges.push_back(edge_t({u, v}));
                }
            }

            std::vector<int> queue_lengths(mpi_rank == 0 ? mpi_world_size : 0);

            int queued_edges_size = queued_edges.size();
            MPI_Gather(&queued_edges_size, 1, MPI_INT, queue_lengths.data(), 1,
                       MPI_INT, 0, MPI_COMM_WORLD);

            std::vector<int> displacements(mpi_rank == 0 ? mpi_world_size + 1
                                                         : 0);

            if (mpi_rank == 0) {
                displacements[0] = 0;
                rep(i, 0, mpi_world_size) displacements[i + 1] =
                    displacements[i] + queue_lengths[i];
                edge_buffer.resize(displacements[mpi_world_size]);
            }

            MPI_Gatherv(queued_edges.data(), queued_edges_size, mpi_edge_t,
                        edge_buffer.data(), queue_lengths.data(),
                        displacements.data(), mpi_edge_t, 0, MPI_COMM_WORLD);

            if (mpi_rank == 0) {
                for (auto [u, v] : edge_buffer) F.merge(u, v);

                if (k <= k2) {
                    // k - k1
                    rep(u, 0, n) scratch[F.head(u)].push_back(u);
                    rep(u, 0, n) {
                        if (scratch[u].size() > 1)
                            trusses[k - k1].push_back(scratch[u]);
                        scratch[u].clear();
                    }
                }
            }
        }

        if (mpi_rank == 0) {
            for (int k = k1; k <= k2; ++k) {
                if (k <= largest_truss) {
                    outf << "1\n";
                    outf << trusses[k - k1].size() << '\n';
                    for (const auto& vec : trusses[k - k1]) {
                        for (auto v : vec) outf << v << ' ';
                        outf << '\n';
                    }
                } else {
                    outf << "0\n";
                }
            }
        }
    } else {
        if (mpi_rank == 0) {
            for (int k = k1; k < k2; ++k) outf << (k <= largest_truss) << ' ';
            outf << (k2 <= largest_truss);
        }
    }
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

    if (mpi_rank == 0) {
        std::cout << "inputpath = " << args.inputpath << "\n";
        std::cout << "headerpath = " << args.headerpath << "\n";
        std::cout << "outputpath = " << args.outputpath << "\n";
        std::cout << "taskid = " << args.taskid << "\n";
        std::cout << "verbose = " << args.verbose << "\n";
        std::cout << "startk = " << args.startk << "\n";
        std::cout << "endk = " << args.endk << "\n";
    }


    const int k1 = args.startk;
    const int k2 = args.endk;
    assert(k1 <= k2);

    graph_t g;

    g.init_io();
    g.init();

    g.compute_truss();
   
    // this has to be synced

    int largest_truss;
    {
        int locally_largest_truss = 0;
        for (int k = k2; k > 0; --k)
            if (!g.fin_bucket[k].empty()) {
                locally_largest_truss = k;
                break;
            }
        MPI_Allreduce(&locally_largest_truss, &largest_truss, 1, MPI_INT,
                      MPI_MAX, MPI_COMM_WORLD);
    }

    if (args.taskid == 1)
        output_1(g, k1, k2, largest_truss);
    else
        output_2(g, k1, k2);

    MPI_Finalize();
}
