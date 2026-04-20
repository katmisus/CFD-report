#ifndef MPI_MESH_HPP
#define MPI_MESH_HPP

#include "mesh/mesh.hpp"
#include "mpi/decomposition.hpp"
#include <mpi.h>
#include <vector>


struct HaloExchange {
    int remote_rank;              // ранг соседа

    std::vector<int> send_ids;    // локальные индексы МОИХ ячеек, которые сосед
                                  // использует как ghost

    std::vector<int> recv_ids;    // локальные индексы GHOST-ячеек в mesh.cells[],
                                  // куда складываем принятые данные от соседа
};


struct LocalMesh {
    Mesh mesh;                        // локальная сетка (свои + ghost ячейки)
    int  n_own = 0;                   // число «своих» ячеек (mesh.cells[0..n_own-1])
    std::vector<HaloExchange> halos;  // по одному элементу на каждого соседа
};


LocalMesh distribute_mesh(const Mesh& global_mesh, MPI_Comm comm);

void exchange_halo(LocalMesh& lm, MPI_Comm comm);


void gather_and_write_vtk(const LocalMesh& lm,
                           const Mesh& global_mesh,
                           const std::vector<int>& part,
                           const std::string& filename,
                           double gamma, MPI_Comm comm);

#endif 
