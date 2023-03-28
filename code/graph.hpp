#pragma once

#include <sys/mman.h>

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

    bool edge_oracle(int u, int v);
    bool local_edge_oracle(int u, int v);
    bool file_edge_oracle(int u, int v);

    void read_header(std::ifstream& gf, std::ifstream& hf);
    void read_owned_graph(std::ifstream& gf);

    std::vector<int> owner;
    std::vector<vertex_t> owned_vertices;
    void assign_owners();

    void init();

    void init_triangles();

    void compute_truss(supp_t k1, supp_t k2);

    idx_t get_edge_idx(vertex_t x, vertex_t y);
    idx_t get_triangle_idx(vertex_t x, idx_t idx_y, vertex_t z);

    int n, m;

    adj_t dodg;
    std::vector<rnk_t> rnk;

    std::vector<std::vector<supp_t>> supp;
    std::vector<std::vector<std::vector<vertex_t>>> inc_tri;
};

