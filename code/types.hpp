#pragma once

#include <mpi.h>

#include <array>
#include <vector>

using vertex_t = int;
using edge_t = std::array<vertex_t, 2>;
using triangle_t = std::array<vertex_t, 3>;

using idx_t = int;
using supp_t = int;

using edge_idx_t = std::pair<vertex_t, idx_t>;

using adj_t = std::vector<std::vector<vertex_t>>;
using rnk_t = uint64_t;

using edge_update_t = std::array<vertex_t, 4>;

// index, p, q, r, idx_p, idx_q
MPI_Datatype mpi_array_int_3;
MPI_Datatype mpi_edge_update_t;
MPI_Datatype mpi_edge_t;

void init_mpi_datatypes() {
    MPI_Type_contiguous(3, MPI_INT, &mpi_array_int_3);
    MPI_Type_commit(&mpi_array_int_3);

    MPI_Type_contiguous(4, MPI_INT, &mpi_edge_update_t);
    MPI_Type_commit(&mpi_edge_update_t);

    MPI_Type_contiguous(2, MPI_INT, &mpi_edge_t);
    MPI_Type_commit(&mpi_edge_t);
}
