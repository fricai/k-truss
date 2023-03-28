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
#include "debug.hpp"
#include "dsu.hpp"
#include "graph.hpp"
#include "macros.hpp"
#include "types.hpp"

int mpi_world_size, mpi_rank;

void graph_t::init() {
    create_mmap();
    init_triangles();
    local_init();
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

    for (auto p : owned_vertices) {
        supp[p].assign(dodg[p].size(), 0);
        inc_tri[p].assign(dodg[p].size(), {});
    }

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

                // p -> q edge exists
                // p -> r edge exists
                // if q -> r edge exists

                if (edge_oracle(q, r)) {
                    ++supp[p][idx_q];
                    ++supp[p][idx_r];

                    inc_tri[p][idx_q].push_back(r);
                    inc_tri[p][idx_r].push_back(q);
                }
            }
        }
    }

    for (auto p : owned_vertices) {
        rep(idx, 0u, inc_tri[p].size()) {
            std::sort(inc_tri[p][idx].begin(), inc_tri[p][idx].end());
            // inc_tri[p][idx].shrink_to_fit();
            /*
             * Should I keep this or not?
             */
        }
    }
}

void compute_truss(int k1, int k2) {

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

    MPI_Finalize();
}
