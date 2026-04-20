#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <map>
#include <filesystem>
#include <numeric>
#include <mpi.h>
#include "FileProcessing.h"
#include "Init.h"
#include "BoundCond.h"
#include "RiemannSolver.h"
#include "Parallel.h"
#include "Domain.h"
#include "GeneralFunctions.h"
#include "Types.h"

int Nx_glob, Ny_glob;
int Nx, Ny;

Domain dom;

int step_fo, step_max, bound_case;

double Lx, Ly, t_max, time_fo, x0, gamm, CFL, Q, C1, C2,
     T_init, R_gas, M, P_min, E_act, Z_freq, VISC, MINWT, GASW, MINGRHO,
     Omega_L, Omega_R, Q_chem;

std::string x_left_bound, x_right_bound,
            y_up_bound, y_down_bound;

std::string high_order_method, TVD_solver, TVD_limiter;
std::string method, solver, time_method, rec_limiter;
bool Diffusion_flag, Viscous_flag, TVD_flag;
int fict = 1;

Field mass_fraction;

std::vector<std::vector<bool>> reacted;

// CHECK: CFL_2D
void GetDt(const Field& W,
           const std::vector<double>& x,
           const std::vector<double>& y,
           double& dt) {

    double dx = Lx / (Nx_glob - 1);
    double dy = Ly / (Ny_glob - 1);

    double max_lambda_x = 0.0;
    double max_lambda_y = 0.0;

    int Nx_cells = Nx - 1;
    int Ny_cells = Ny - 1;

    for (int i = fict; i < Nx_cells + fict; i++) {
        for (int j = fict; j < Ny_cells + fict; ++j) {
            double rho = W[i][j][0];
            double u   = W[i][j][1];
            double v   = W[i][j][2];
            double P   = W[i][j][NEQ - 1];

            double c = std::sqrt(gamm * P / rho);

            max_lambda_x = std::max(max_lambda_x, std::abs(u) + c);
            max_lambda_y = std::max(max_lambda_y, std::abs(v) + c);
        }
    }

    double dt_x = dx / max_lambda_x;
    double dt_y = dy / max_lambda_y;
    dt = CFL * std::min(dt_x, dt_y);
}


int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int rank, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    double t_start, t_end;

    if (argc < 2) {
        if (rank == 0)
            std::cerr << "Используйте: ./test <config_set_name>\n";
        MPI_Finalize();
        return 1;
    }

    fs::path config_path   = argv[1];
    std::string config_name = config_path.stem().string();

    fs::path base_output = fs::path("output") / config_name;
    fs::path csv_folder  = base_output / "CSV";
    fs::path pics_folder = base_output / "pics";

    if (rank == 0) {
        if (fs::exists(csv_folder)) {
            std::cout << "Очищаем папку: " << csv_folder << std::endl;
            fs::remove_all(csv_folder);
        }
        fs::create_directories(csv_folder);
        fs::create_directories(pics_folder);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    readConfig(config_path.string());

    std::vector<int> global_sizes = {Nx_glob, Ny_glob};
    std::vector<int> proc_dims;
    DecomposeProcesses(p, 2, global_sizes, proc_dims);

    int px = proc_dims[0];
    int py = proc_dims[1];

    dom = BuildDomain(rank, px, py, Nx_glob, Ny_glob);
    Nx  = dom.Nx;
    Ny  = dom.Ny;
    BuildGrid(dom, fict, Nx_glob, Ny_glob);

    int Nx_tot = Nx + 2*fict - 1;
    int Ny_tot = Ny + 2*fict - 1;

    mass_fraction = Field(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    reacted = std::vector<std::vector<bool>>(Nx_tot, std::vector<bool>(Ny_tot, false));

    Field W_0(Nx_tot, std::vector<State>(Ny_tot));
    InitValues(W_0, dom.x, dom.y, config_path);

    for (int i = 0; i < Nx_tot; i++)
        for (int j = 0; j < Ny_tot; j++)
            if (mass_fraction[i][j][0] < GASW)
                reacted[i][j] = true;

    BoundCond(W_0, dom);
    ExchangeGhostCells(W_0, dom);

    Field W     = W_0;
    Field W_new = W_0;

    double t = 0.0, dt = 1.0;
    int step = 0;

    std::ofstream outFile("values.txt", std::ios::app);
    if (method == "FLIC") {
        if (outFile.is_open()) {
            outFile << "CFL = " << CFL << std::endl;
        } else {
            std::cerr << "Не удалось открыть файл!" << std::endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t_start = MPI_Wtime();

    std::vector<double> dt_for_FLIC;
    while (step <= step_max && t <= t_max) {

        GetDt(W, dom.x, dom.y, dt);

        MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        if (t + dt > t_max)
            dt = t_max - t;

        t += dt;

        UpdateArrays(W, W_new, dom.x, dom.y, dt);
        // CHECK: FLIC_CFL
        if (method == "FLIC") dt_for_FLIC.push_back(dt);

        BoundCond(W, dom);
        ExchangeGhostCells(W, dom);

        if (step % step_fo == 0) {
            std::string filename = (csv_folder /
                                   (std::to_string(step) + "_step.csv")).string();
            for (int r = 0; r < p; r++) {
                MPI_Barrier(MPI_COMM_WORLD);
                if (dom.rank == r) {
                    bool append = (r != 0);
                    SaveFieldToCSV(W, dom.x, dom.y, t, filename, append);
                }
            }
        }

        if (t >= t_max) {
            std::string filename_step  = (csv_folder /
                                         (std::to_string(step) + "_step.csv")).string();
            std::string filename_final = (csv_folder / "Final.csv").string();

            for (int r = 0; r < p; r++) {
                MPI_Barrier(MPI_COMM_WORLD);
                if (dom.rank == r) {
                    bool append = (r != 0);
                    SaveFieldToCSV(W, dom.x, dom.y, t, filename_step,  append);
                    SaveFieldToCSV(W, dom.x, dom.y, t, filename_final, append);
                }
            }
            break;
        }

        step++;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    t_end = MPI_Wtime();
    
    if (method == "FLIC") {
        double sum = std::accumulate(dt_for_FLIC.begin(), dt_for_FLIC.end(), 0.0);
        double mean = sum / dt_for_FLIC.size(); 
		if (outFile.is_open()) {
			outFile << "dt_mean = " << mean << std::endl;;
			outFile.close();
		} else {
			std::cerr << "Не удалось открыть файл!" << std::endl;
		}
    }

    if (rank == 0)
        std::cout << "Время расчёта: " << (t_end - t_start) << " секунд\n";
    if (rank == 0)
        std::cout << "Завершено успешно." << std::endl;

    MPI_Finalize();
    return 0;
}