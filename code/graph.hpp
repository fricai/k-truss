#pragma once

#include <algorithm>
#include <fstream>
#include <list>
#include <vector>

#include "types.hpp"

struct graph_t {
    std::vector<uint32_t> offset;
    uint32_t const* file_map;  // underlying data is const
    uint64_t file_len;
    void create_mmap();

    ~graph_t() { munmap((void*)file_map, file_len); }

    bool file_edge_oracle(int u, int v);

    void read_header(std::ifstream& gf, std::ifstream& hf);
    void read_complete_graph(std::ifstream& gf);
    void read_owned_graph(std::ifstream& gf);
    void read_from_stdin();

    std::vector<int> owner;
    std::vector<vertex_t> owned_vertices;
    void assign_owners();

    bool triangle_edge_query(vertex_t x, vertex_t y, vertex_t z);
    bool nonlocal_triangle_edge_query(vertex_t x, vertex_t y, vertex_t z);
    void local_triangle_update(vertex_t p, idx_t idx_q, idx_t idx_r);

    void init();
    void compute_truss_range();
    void construct_dodg();
    void init_triangles();
    void local_init();

    void compute_truss(supp_t k1, supp_t k2);
    void update_edge(vertex_t u, vertex_t v, vertex_t w, supp_t val_new);

    idx_t get_edge_idx(vertex_t x, vertex_t y) {
        return std::lower_bound(dodg[x].begin(), dodg[x].end(), y) -
               dodg[x].begin();
    }

    idx_t get_triangle_idx(vertex_t x, idx_t idx_y, vertex_t z) {
        return std::lower_bound(inc_tri[x][idx_y].begin(),
                                inc_tri[x][idx_y].end(), z) -
               inc_tri[x][idx_y].begin();
    }

    std::vector<std::vector<edge_update_t>> buffered_edge_updates;
    void send_update_edge(vertex_t u, vertex_t v, vertex_t w, supp_t val_new);

    int n, m;
    adj_t undirected_graph;

    adj_t dodg;
    std::vector<rnk_t> rnk;

    std::vector<std::vector<supp_t>> supp;

    std::vector<std::vector<supp_t>> g;
    std::vector<std::vector<std::vector<supp_t>>> h;

    std::vector<std::vector<std::vector<vertex_t>>> inc_tri;
    std::vector<std::vector<std::vector<supp_t>>> tri_supp;

    std::vector<std::list<edge_idx_t>> bucket;
    std::vector<std::vector<std::list<edge_idx_t>::iterator>> bucket_iter;

    std::vector<edge_idx_t> changed_tau;
    std::vector<std::vector<bool>> was_changed;

    supp_t k_min, k_max;
};

