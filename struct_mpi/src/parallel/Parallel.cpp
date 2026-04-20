#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

// extern int Nx_glob_cells, Ny_glob_cells;
// extern double Lx, Ly;

void DecomposeProcesses(int p,
                        int dim,
                        const std::vector<int>& global_sizes,
                        std::vector<int>& proc_dims) {
    proc_dims.assign(dim, 1);

    if (dim == 1) {
        proc_dims[0] = p;
        return;
    }

    double best_score = std::numeric_limits<double>::max();

    if (dim == 2) {
        double Nx = global_sizes[0];
        double Ny = global_sizes[1];

        double r_grid = std::min(Nx, Ny) / std::max(Nx, Ny);

        double best_score = std::numeric_limits<double>::max();

        for (int px = 1; px <= p; px++) {
            if (p % px != 0) continue;

            int py = p / px;

            double r_proc = (double)std::min(px, py) / std::max(px, py);

            double score = std::abs(r_grid - r_proc);

            if (score < best_score) {
                best_score = score;
                proc_dims[0] = px;
                proc_dims[1] = py;
            }
        }

        if (Nx > Ny && proc_dims[0] < proc_dims[1])
        std::swap(proc_dims[0], proc_dims[1]);

        if (Ny > Nx && proc_dims[1] < proc_dims[0])
        std::swap(proc_dims[0], proc_dims[1]);
    }

    if (dim == 3) {
        for (int px = 1; px <= p; px++) {
            if (p % px) continue;

            int rest = p / px;

            for (int py = 1; py <= rest; py++) {
                if (rest % py) continue;

                int pz = rest / py;

                double score = 0;

                double rg_xy = (double)global_sizes[0] / global_sizes[1];
                double rp_xy = (double)px / py;

                double rg_xz = (double)global_sizes[0] / global_sizes[2];
                double rp_xz = (double)px / pz;

                double rg_yz = (double)global_sizes[1] / global_sizes[2];
                double rp_yz = (double)py / pz;

                score += std::abs(rg_xy - rp_xy);
                score += std::abs(rg_xz - rp_xz);
                score += std::abs(rg_yz - rp_yz);

                if (score < best_score) {
                    best_score = score;

                    proc_dims[0] = px;
                    proc_dims[1] = py;
                    proc_dims[2] = pz;
                }
            }
        }
    }
}

int LocalSize(int N, int p_coord, int p_dim) {
    int base = N / p_dim;
    int rest = N % p_dim;

    if (p_coord < rest)
        return base + 1;
    else
        return base;
}

int Offset(int N_cells, int coord, int p_dim) {

    int base = N_cells / p_dim;
    int rest = N_cells % p_dim;

    if (coord < rest) {
        return coord * (base + 1);
    } else {
        return rest * (base + 1) + (coord - rest) * base;
    }
}

