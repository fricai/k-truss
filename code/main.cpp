#include <fcntl.h>
#include <mpi.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <cassert>
#include <chrono>
#include <iostream>
#include <queue>
#include <thread>
#include <vector>

#include "argparse.hpp"
#include "debug_output.hpp"
#include "dsu.hpp"
#include "graph.hpp"
#include "macros.hpp"
#include "types.hpp"

int mpi_world_size, mpi_rank;

void graph_t::create_mmap() {
    int fd = open(args.inputpath.data(), O_RDONLY);

    struct stat statbuf;
    fstat(fd, &statbuf);
    file_len = statbuf.st_size;

    file_map =
        (uint32_t*)mmap(nullptr, file_len, PROT_READ, MAP_PRIVATE, fd, 0);
}

bool graph_t::file_edge_oracle(int u, int v) {
    return std::binary_search(file_map + offset[u] + 2,
                              file_map + offset[u + 1], (uint32_t)v);
}

void graph_t::init() {
    create_mmap();

    init_triangles();

    compute_truss_range();
    local_init();
}

void graph_t::assign_owners() {
    using I = int64_t;
    std::vector<I> load(n);
    rep(u, 0, n) load[u] = (rnk[u] >> 32) * (rnk[u] >> 32);

    std::priority_queue<std::pair<I, int>> worlds;
    rep(i, 0, mpi_world_size) worlds.push({0, i});

    owner.resize(n);
    rep(u, 0, n) {
        auto [ld, id] = worlds.top();
        worlds.pop();

        ld -= load[u];
        owner[u] = id;

        if (id == mpi_rank) owned_vertices.push_back(u);

        worlds.push({ld, id});
    }
}

constexpr supp_t TRUSS_INF = 1e6 + 20;

void graph_t::read_header(std::ifstream& gf, std::ifstream& hf) {
    gf.read((char*)&n, 4);
    gf.read((char*)&m, 4);

    offset.resize(n + 1);
    hf.read((char*)offset.data(), 4 * n);

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

void graph_t::read_owned_graph(std::ifstream& gf) {
    dodg.resize(n);

    std::vector<uint32_t> scratch;
    scratch.reserve(n);
    for (auto u : owned_vertices) {
        gf.seekg(4 * offset[u]);

        uint32_t id, deg;
        gf.read((char*)&id, 4);
        gf.read((char*)&deg, 4);

        assert(id == (uint32_t)u);
        assert(deg == (rnk[u] >> 32));

        scratch.resize(deg);
        gf.read((char*)scratch.data(), 4 * deg);

        for (auto v : scratch)
            if (rnk[u] < rnk[v]) dodg[u].push_back(v);

        assert(std::is_sorted(dodg[u].begin(), dodg[u].end()));
    }
}

// must be called after read_header
void graph_t::read_complete_graph(std::ifstream& gf) {
    gf.seekg(4 + 4);

    dodg.resize(n);

    std::vector<uint32_t> scratch;
    scratch.reserve(n);

    rep(u, 0, n) {
        uint32_t id, deg;
        gf.read((char*)&id, 4);
        gf.read((char*)&deg, 4);
        assert(id == (uint32_t)u);
        assert(deg == (rnk[u] >> 32));

        scratch.resize(deg);
        gf.read((char*)scratch.data(), 4 * deg);

        for (auto v : scratch)
            if (rnk[u] < rnk[v]) dodg[u].push_back(v);
        assert(std::is_sorted(dodg[u].begin(), dodg[u].end()));
    }
}

void graph_t::read_from_stdin() {
    std::cin >> n >> m;

    undirected_graph.resize(n);
    rep(i, 0, n) {
        int idx, deg;
        std::cin >> idx >> deg;

        undirected_graph[idx].resize(deg);
        for (auto& v : undirected_graph[idx]) std::cin >> v;
    }

    construct_dodg();
}

void graph_t::construct_dodg() {
    dodg.resize(n);
    rnk.resize(n);

    rep(i, 0, n) rnk[i] = rnk_t(undirected_graph[i].size()) << 32 | i;
    rep(i, 0, n) {
        for (auto j : undirected_graph[i])
            if (rnk[i] < rnk[j]) dodg[i].push_back(j);
        undirected_graph[i].clear();
    }

    for (auto& v : dodg) std::sort(v.begin(), v.end());
}

bool graph_t::nonlocal_triangle_edge_query(vertex_t q, vertex_t r, vertex_t p) {
    assert(owner[q] == mpi_rank);

    auto idx_qr = get_edge_idx(q, r);
    if ((size_t)idx_qr == dodg[q].size() || dodg[q][idx_qr] != r) {
        assert(false);
        return false;
    }
    // now, (q, r) edge exists
    ++supp[q][idx_qr];
    inc_tri[q][idx_qr].push_back(p);
    return true;
}

bool graph_t::triangle_edge_query(vertex_t q, vertex_t r, vertex_t p) {
    assert(owner[q] == mpi_rank);

    auto idx_qr = get_edge_idx(q, r);
    if ((size_t)idx_qr == dodg[q].size() || dodg[q][idx_qr] != r) return false;
    // now, (q, r) edge exists
    ++supp[q][idx_qr];
    inc_tri[q][idx_qr].push_back(p);
    return true;
}

void graph_t::local_triangle_update(vertex_t p, idx_t idx_q, idx_t idx_r) {
    auto q = dodg[p][idx_q];
    auto r = dodg[p][idx_r];

    ++supp[p][idx_q];
    ++supp[p][idx_r];

    inc_tri[p][idx_q].push_back(r);
    inc_tri[p][idx_r].push_back(q);
}

void graph_t::init_triangles() {
    supp.resize(n);
    inc_tri.resize(n);
    tri_supp.resize(n);

    for (auto p : owned_vertices) {
        supp[p].assign(dodg[p].size(), 0);
        inc_tri[p].assign(dodg[p].size(), {});
        tri_supp[p].assign(dodg[p].size(), {});
    }

    // compute support

    // p -> q

    std::vector<std::vector<std::array<int, 3>>> queries(mpi_world_size);
    std::vector<std::vector<std::array<int, 2>>> query_idxs(mpi_world_size);

    for (auto p : owned_vertices) {
        rep(idx_q_, 0u, dodg[p].size()) {
            rep(idx_r_, idx_q_ + 1, dodg[p].size()) {
                vertex_t q = dodg[p][idx_q_], r = dodg[p][idx_r_];
                idx_t idx_q = idx_q_, idx_r = idx_r_;

                if (rnk[q] > rnk[r]) {
                    std::swap(idx_q, idx_r);
                    std::swap(q, r);
                }

                // rnk[p] < rnk[q] < rnk[r]

                if (owner[q] == mpi_rank) {
                    if (triangle_edge_query(q, r, p))
                        local_triangle_update(p, idx_q, idx_r);
                } else if (file_edge_oracle(q, r)) {
                    queries[owner[q]].push_back(std::array<int, 3>({q, r, p}));
                    query_idxs[owner[q]].push_back(
                        std::array<int, 2>({idx_q, idx_r}));
                }

                /*
                if (triangle_edge_query(q, r, p)) {
                    local_triangle_update(p, idx_q, idx_r);
                    // it updates the values stored at q too
                    // so the function only needs to think about local
                    // updates

                    // in async mode, the result would be received much
                    // later
                }
                */
            }
        }
    }

    std::vector<MPI_Request> query_requests(mpi_world_size);
    rep(i, 0, mpi_world_size) {
        if (i == mpi_rank) continue;
        MPI_Isend(queries[i].data(), queries[i].size(), mpi_array_int_3, i, 0,
                  MPI_COMM_WORLD, &query_requests[i]);
    }

    // sent out requests to everyone

    std::vector<std::vector<uint8_t>> replies_sent(mpi_world_size);
    std::vector<MPI_Request> reply_requests(mpi_world_size);

    {
        std::vector<std::array<int, 3>> scratch;
        rep(i, 0, mpi_world_size) {
            if (i == mpi_rank) continue;

            MPI_Status status;
            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);

            int buf_size;
            MPI_Get_count(&status, mpi_array_int_3, &buf_size);

            scratch.resize(buf_size);
            MPI_Recv(scratch.data(), buf_size, mpi_array_int_3, i, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            replies_sent[i].resize(buf_size);

            rep(j, 0, buf_size) {
                auto [q, r, p] = scratch[j];
                replies_sent[i][j] = nonlocal_triangle_edge_query(q, r, p);
            }

            MPI_Isend(replies_sent[i].data(), buf_size, MPI_BYTE, i, 0,
                      MPI_COMM_WORLD, &reply_requests[i]);
        }
    }

    {
        std::vector<uint8_t> replies_recvd;
        rep(i, 0, mpi_world_size) {
            if (i == mpi_rank) continue;

            const int len = queries[i].size();

            replies_recvd.resize(len);

            MPI_Recv(replies_recvd.data(), len, MPI_BYTE, i, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            rep(j, 0, len) {
                if (replies_recvd[j]) {
                    auto p = queries[i][j][2];
                    auto [idx_q, idx_r] = query_idxs[i][j];
                    local_triangle_update(p, idx_q, idx_r);
                }
            }
        }
    }

    rep(i, 0, mpi_world_size) {
        if (i == mpi_rank) continue;
        MPI_Wait(&query_requests[i], MPI_STATUS_IGNORE);
        MPI_Wait(&reply_requests[i], MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // I don't really need this barrier I think

    // tri_supp do be empty tho
    for (auto p : owned_vertices) {
        rep(idx, 0u, inc_tri[p].size()) {
            std::sort(inc_tri[p][idx].begin(), inc_tri[p][idx].end());
            tri_supp[p][idx].assign(inc_tri[p][idx].size(), TRUSS_INF);
        }
    }
}

void graph_t::local_init() {
    g = supp;
    h.resize(n);

    was_changed.resize(n);

    bucket.assign(k_max + 1, {});
    bucket_iter.resize(n);
    for (auto p : owned_vertices) {
        const auto out_deg = dodg[p].size();

        was_changed[p].assign(out_deg, false);

        h[p].resize(out_deg, {});
        rep(idx_q, 0u, out_deg) h[p][idx_q].assign(supp[p][idx_q], 0);

        bucket_iter[p].resize(out_deg);
        rep(idx_q, 0u, out_deg) {
            auto& cur_bucket = bucket[supp[p][idx_q]];
            cur_bucket.emplace_front(p, idx_q);
            bucket_iter[p][idx_q] = cur_bucket.begin();
        }
    }
}

void graph_t::compute_truss_range() {
    k_min = TRUSS_INF;
    k_max = 0;

    for (auto p : owned_vertices) {
        if (!supp[p].empty()) {
            auto [mn, mx] = minmax_element(supp[p].begin(), supp[p].end());
            k_min = std::min(k_min, *mn);
            k_max = std::max(k_max, *mx);
        }
    }

    int k_min_recv, k_max_recv;

    MPI_Allreduce(&k_min, &k_min_recv, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&k_max, &k_max_recv, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    k_min = k_min_recv;
    k_max = k_max_recv;
}

void graph_t::compute_truss(supp_t k1, supp_t k2) {
    // if k1 == 0, then k1 - 1 = 0, but window_l >= 0
    supp_t window_l = k_min, window_r = window_l;

    // active set is set of edges

    static constexpr float delta = 0.1;

    std::vector<edge_idx_t> active;

    using I = int64_t;
    I gamma_max = 0, gamma_act = 0;

    // gamma[e] = supp[e]
    //
    // and gamma[set] = sum(gamma[e] | e in set)

    auto reserve_active = [&](size_t target) {
        auto sz = std::max((size_t)1, active.capacity());
        while (sz < target) sz = 2 * sz;
        if (sz > active.capacity()) active.reserve(sz);
    };

    auto add_to_active = [&](const auto& to_add) {
        reserve_active(active.size() + to_add.size());

        for (auto [u, idx_v] : to_add) {
            // I could have stored (u, idx_v) instead of e
            gamma_act += supp[u][idx_v];
            active.push_back({u, idx_v});
            // gamma[e] = supp[e[0]][get_edge_idx(e[0], e[1])]
        }
    };

    auto clear_active = [&]() {
        gamma_act = 0;
        active.clear();
    };

    auto window_expansion = [&]() {
        auto& k = window_r;
        while (k < k_max and gamma_act <= delta * gamma_max) {
            // add bucket[k + 1] to active
            ++k;
            add_to_active(bucket[k]);
        }
    };

    rep(k, window_l, window_r + 1) add_to_active(bucket[k]);

    // assume that each machine is responsible for precisely one edge

    buffered_edge_updates.resize(mpi_world_size);

    std::vector<edge_update_t> _edge_update_scratch;

    // globalize this check across all nodes
    bool global_not_running;
    while (true) {
        bool not_running = active.empty() and window_r == k_max;
        MPI_Allreduce(&not_running, &global_not_running, 1, MPI_CXX_BOOL,
                      MPI_LAND, MPI_COMM_WORLD);
        if (global_not_running) break;

        window_expansion();

        for (auto [u, idx_v] : active) {
            // u -> v

            auto v = dodg[u][idx_v];

            rep(idx_w, 0u, inc_tri[u][idx_v].size()) {
                vertex_t w = inc_tri[u][idx_v][idx_w];
                auto& val = tri_supp[u][idx_v][idx_w];

                if (supp[u][idx_v] < val) {
                    val = supp[u][idx_v];

                    // send out updates
                    send_update_edge(u, w, v, val);
                    send_update_edge(v, w, u, val);
                }
            }
        }

        std::vector<MPI_Request> send_requests(mpi_world_size);
        rep(i, 0, mpi_world_size) {
            if (i == mpi_rank) continue;
            MPI_Isend(buffered_edge_updates[i].data(),
                      buffered_edge_updates[i].size(), mpi_edge_update_t, i, 0,
                      MPI_COMM_WORLD, &send_requests[i]);
        }

        rep(i, 0, mpi_world_size) {
            if (i == mpi_rank) continue;

            MPI_Status status;
            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);

            int buf_size;
            MPI_Get_count(&status, mpi_edge_update_t, &buf_size);

            _edge_update_scratch.resize(buf_size);
            MPI_Recv(_edge_update_scratch.data(), buf_size, mpi_edge_update_t,
                     i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (auto [u, v, w, val_new] : _edge_update_scratch)
                update_edge(u, v, w, val_new);
        }
        rep(i, 0, mpi_world_size) if (i != mpi_rank)
            MPI_Wait(&send_requests[i], MPI_STATUS_IGNORE);

        /*
        // send out all buffered updates
        for (auto& upd : buffered_edge_updates) {
            MPI_Isend(upd.data(), 1, mpi_edge_update_t, owner[upd[0]], 0,
                      MPI_COMM_WORLD, &__update_req_sent);
            MPI_Request_free(&__update_req_sent);
        }
        // broadcast we're done with sending queries
        edge_update_t update_done({-1, -1, -1, -1});
        MPI_Ibcast(update_done.data(), 1, mpi_edge_update_t, mpi_rank,
                   MPI_COMM_WORLD, &__update_req_sent);
        MPI_Request_free(&__update_req_sent);

        int n_nodes_done = 0;
        edge_update_t rcvd_update;
        MPI_Status recv_status;
        // receive updates and update the edges
        while (n_nodes_done != mpi_world_size - 1) {
            MPI_Recv(rcvd_update.data(), 1, mpi_edge_update_t, MPI_ANY_SOURCE,
                     0, MPI_COMM_WORLD, &recv_status);
            if (rcvd_update[0] == -1) {
                n_nodes_done += 1;
            } else {
                update_edge(rcvd_update[0], rcvd_update[1], rcvd_update[2],
                            rcvd_update[3]);
            }
        }
        */

        // barrier here
        MPI_Barrier(MPI_COMM_WORLD);
        // may be don't need this^

        clear_active();

        reserve_active(active.size() + changed_tau.size());
        for (auto [u, idx_v] : changed_tau) {
            was_changed[u][idx_v] = false;  // unset was_changed

            if (supp[u][idx_v] <= window_r) {
                gamma_act += supp[u][idx_v];
                active.push_back({u, idx_v});
            }
            // gamma[e] = supp(e)
        }

        changed_tau.clear();

        gamma_max = std::max(gamma_max, gamma_act);
    }
}

void graph_t::send_update_edge(vertex_t u, vertex_t v, vertex_t w,
                               supp_t val_new) {
    // not guaranteed that it's a u->v edge;
    if (rnk[u] > rnk[v]) std::swap(u, v);

    // tell owner of u to do update_edge
    // when owner receives update, it will call update_edge

    if (owner[u] == mpi_rank) {
        // we own it so do a local update
        update_edge(u, v, w, val_new);
    } else {
        // ~~send out the update~~ enqueue the update for sending
        edge_update_t update({u, v, w, val_new});
        buffered_edge_updates[owner[u]].push_back(update);
        // MPI_Isend(update, 1, mpi_edge_update_t, owner[u], 0, MPI_COMM_WORLD,
        //           &__update_req_sent);
        // MPI_Request_free(&__update_req_sent);
    }
}

void graph_t::update_edge(vertex_t u, vertex_t v, vertex_t w, supp_t val_new) {
    assert(owner[u] == mpi_rank);

    auto idx_uv = get_edge_idx(u, v);
    auto idx_uvw = get_triangle_idx(u, idx_uv, w);
    auto& val_old = tri_supp[u][idx_uv][idx_uvw];

    auto& tau = supp[u][idx_uv];

    if (val_new < val_old) {
        if (val_new < tau) {
            if (val_old < tau)
                --h[u][idx_uv][val_old];
            else
                --g[u][idx_uv];
            ++h[u][idx_uv][val_new];
        }
        val_old = val_new;
    }

    if (g[u][idx_uv] < tau) {
        bucket[tau].erase(bucket_iter[u][idx_uv]);
        --tau;
        bucket[tau].emplace_front(u, idx_uv);
        bucket_iter[u][idx_uv] = bucket[tau].begin();

        g[u][idx_uv] += h[u][idx_uv][tau];

        // changed tau[e]
        if (!was_changed[u][idx_uv]) {
            changed_tau.emplace_back(u, idx_uv);
            was_changed[u][idx_uv] = true;
        }
    }
}

void output_2(const graph_t& g, int k1, int k2, int largest_truss) {
    const int n = g.n;

    dsu F(n);

    std::vector<edge_t> queued_edges;

    std::vector<edge_t> spanning_forest;

    for (int k = largest_truss; k >= k2; --k) {
        for (auto [u, idx_v] : g.bucket[k]) {
            auto v = g.dodg[u][idx_v];

            if (F.merge(u, v)) {
                // send to 0 to add to 0's DSU
                (mpi_rank == 0 ? spanning_forest : queued_edges)
                    .push_back(edge_t({u, v}));
            }
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
    rep(u, 0, n) if (!grps[u].empty()) grps[u].back() = '\n';

    std::vector<bool> vis(n);
    std::vector<int> heads;
    auto proxy = [&](int u) {
        u = F.head(u);
        if (grps[u].empty() or vis[u]) return;
        vis[u] = true;
        heads.push_back(u);
    };

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, args.outputpath.data(),
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    std::string buf;

    bool found_any = false;

    int fst = (n + mpi_world_size - 1) / mpi_world_size * mpi_rank,
        lst = std::min(fst / mpi_rank * (mpi_rank + 1), n);
    assert(fst <= n);

    rep(u, fst, lst) {
        uint32_t const* const neighbourhood = g.file_map + g.offset[u] + 2;
        const int deg = g.offset[u + 1] - g.offset[u] - 2;

        for (auto v : heads) vis[v] = false;
        heads.clear();

        rep(i, 0, deg) proxy(neighbourhood[i]);

        if ((int)heads.size() < args.p) continue;

        found_any = true;

        buf.clear();
        buf += std::to_string(u) + "\n";
        if (args.verbose) {
            buf += std::to_string(heads.size()) + "\n";
            for (auto r : heads) buf += grps[r];
        }
        MPI_File_write_shared(fh, buf.data(), buf.size(), MPI_CHAR,
                              MPI_STATUS_IGNORE);
    }

    bool global_found_any;
    MPI_Allreduce(&found_any, &global_found_any, 1, MPI_CXX_BOOL, MPI_LOR,
                  MPI_COMM_WORLD);

    if (!global_found_any and mpi_rank == 0) {
        buf = "-1\n";
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
            for (auto [u, idx_v] : g.bucket[k]) {
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
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    init_mpi_datatypes();

    args.taskid = 1;
    args.verbose = false;

    argp_parse(&argp, argc, argv, 0, 0, &args);

    std::ifstream gf(args.inputpath, std::ios::binary);
    std::ifstream hf(args.headerpath, std::ios::binary);

    const int k1 = args.startk;
    const int k2 = args.endk;
    assert(k1 <= k2);

    graph_t g;

    g.read_header(gf, hf);

    g.assign_owners();
    g.read_owned_graph(gf);

    g.init();
    g.compute_truss(k1, k2);

    // this has to be synced

    int largest_truss;
    {
        int locally_largest_truss = 0;
        for (int k = g.k_max; k > 0; --k)
            if (!g.bucket[k].empty()) {
                locally_largest_truss = k;
                break;
            }
        MPI_Allreduce(&locally_largest_truss, &largest_truss, 1, MPI_INT,
                      MPI_MAX, MPI_COMM_WORLD);
    }

    if (args.taskid == 1)
        output_1(g, k1, k2, largest_truss);
    else
        output_2(g, k1, k2, largest_truss);

    MPI_Finalize();
}
