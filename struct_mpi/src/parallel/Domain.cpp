#include "Domain.h"
#include "Parallel.h"
#include <mpi.h>
#include <fstream>
#include <vector>
#include <filesystem>

extern double Lx, Ly;

Domain BuildDomain(int rank,
                   int px,
                   int py,
                   int Nx_glob, 
                   int Ny_glob) {
    Domain d;

    d.rank = rank;

    d.px = px;
    d.py = py;

    d.rx = rank % px;
    d.ry = rank / px;

    int Nx_glob_cells = Nx_glob - 1;
    int Ny_glob_cells = Ny_glob - 1;

    d.Nx_cells = LocalSize(Nx_glob_cells, d.rx, px);
    d.Ny_cells = LocalSize(Ny_glob_cells, d.ry, py);

    d.offset_x = Offset(Nx_glob_cells, d.rx, px);
    d.offset_y = Offset(Ny_glob_cells, d.ry, py);

    d.Nx = d.Nx_cells + 1;
    d.Ny = d.Ny_cells + 1;

    return d;
}


void BuildGrid(Domain& dom, int fict, int Nx_glob, int Ny_glob) {

    double dx = Lx / (Nx_glob - 1);
    double dy = Ly / (Ny_glob - 1);

    dom.x.resize(dom.Nx + 2*fict);
    dom.y.resize(dom.Ny + 2*fict);

    for (int i = 0; i < dom.Nx + 2*fict; i++) {
        int global_i = dom.offset_x + i - fict;
        dom.x[i] = global_i * dx;
    }

    for (int j = 0; j < dom.Ny + 2*fict; j++) {
        int global_j = dom.offset_y + j - fict;
        dom.y[j] = global_j * dy;
    }
}


void WriteDomainSizes(const Domain& dom,
                      int p,
                      const std::filesystem::path& base_output) {

    const int NDATA = 8;

    int local_data[NDATA] = {
        dom.rank,
        dom.rx,
        dom.ry,
        dom.Nx,
        dom.Ny,
        dom.Nx_cells,
        dom.offset_x,
        dom.offset_y
    };

    std::vector<int> all_data;

    if (dom.rank == 0) {
        all_data.resize(NDATA * p);
    }

    MPI_Gather(local_data, NDATA, MPI_INT,
               all_data.data(), NDATA, MPI_INT,
               0, MPI_COMM_WORLD);

    if (dom.rank == 0) {

        std::ofstream file(base_output / "domain_sizes.txt");

        file << "rank rx ry Nx Ny Nx_cells offset_x offset_y\n";

        for (int i = 0; i < p; i++) {

            for (int j = 0; j < NDATA; j++) {
                file << all_data[i * NDATA + j] << " ";
            }

            file << "\n";
        }
    }
}

void WriteDomains(const Domain& dom,
                  const std::vector<double>& x,
                  const std::vector<double>& y,
                  int fict,
                  int p,
                  const std::filesystem::path& base_output) {

    double xmin = x[fict];
    double xmax = x[fict + dom.Nx];

    double ymin = y[fict];
    double ymax = y[fict + dom.Ny];

    double xmin_g = x[0];
    double xmax_g = x[dom.Nx + 2*fict - 1];

    double ymin_g = y[0];
    double ymax_g = y[dom.Ny + 2*fict - 1];

    const int NDATA = 11;

    double local_data[NDATA] = {
        (double)dom.rank,
        (double)dom.rx,
        (double)dom.ry,
        xmin,
        xmax,
        ymin,
        ymax,
        xmin_g,
        xmax_g,
        ymin_g,
        ymax_g
    };

    std::vector<double> all_data;

    if (dom.rank == 0) {
        all_data.resize(NDATA * p);
    }

    MPI_Gather(local_data, NDATA, MPI_DOUBLE,
               all_data.data(), NDATA, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    if (dom.rank == 0) {

        std::ofstream file(base_output / "domains.txt");

        file << "rank rx ry xmin xmax ymin ymax xmin_g xmax_g ymin_g ymax_g\n";

        for (int i = 0; i < p; i++) {

            for (int j = 0; j < NDATA; j++) {
                file << all_data[i * NDATA + j] << " ";
            }

            file << "\n";
        }
    }
}