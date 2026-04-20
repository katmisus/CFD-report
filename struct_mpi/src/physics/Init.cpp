#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

#include "ParseTOML.h"
#include "Types.h"

extern int Nx, Ny;
extern int Nx_glob, Ny_glob;
extern int step_fo, step_max, bound_case;

extern double Lx, Ly, t_max, time_fo, x0, gamm, CFL, Q, C1, C2,
              T_init, R_gas, M, P_min, E_act, Z_freq, VISC, MINWT, GASW, MINGRHO,
              Omega_L, Omega_R, D, Q_chem;

extern std::string x_left_bound, x_right_bound,
                   y_up_bound, y_down_bound;

extern std::string high_order_method, TVD_solver, TVD_limiter;
extern std::string method, solver, time_method, rec_limiter;
extern bool Diffusion_flag, Viscous_flag, TVD_flag;
extern int fict;
extern Field mass_fraction;

void readConfig(const std::string& config_path) {

    SimpleToml toml;

    if (!toml.load(config_path)) {
        std::cerr << "Нет config-файла" << std::endl;
        return;
    }

    gamm = toml.root["simulation"].table["gamma"].number;
    auto& scheme = toml.root["scheme"].table;

    method      = scheme["method"].str;
    rec_limiter = scheme["rec_limiter"].str;
    solver      = scheme["solver"].str;
    time_method = scheme["time_integration_method"].str;

    Diffusion_flag = (scheme["Diffusion"].str == "On") ? true : false;
    Q = (Diffusion_flag == true) ? scheme["Q"].number : 0.0;

    Viscous_flag = (scheme["Viscous"].str == "On") ? true : false;
    C1 = (Viscous_flag == true) ? scheme["C1"].number : 0.0;
    C2 = (Viscous_flag == true) ? scheme["C2"].number : 0.0;

    TVD_flag          = (scheme["TVD"].str == "On") ? true : false;
    high_order_method = scheme["High_order_method"].str;
    TVD_solver        = scheme["High_order_method"].str;
    TVD_limiter       = scheme["High_order_method"].str;

    Nx_glob = scheme["N_x"].number;
    Ny_glob = scheme["N_y"].number;
    Lx      = scheme["L_x"].number;
    Ly      = scheme["L_y"].number;

    CFL           = scheme["CFL"].number;
    x_left_bound  = scheme["x_left_bound"].str;
    x_right_bound = scheme["x_right_bound"].str;
    y_up_bound    = scheme["y_up_bound"].str;
    y_down_bound  = scheme["y_down_bound"].str;

    T_init  = scheme["T_init"].number;
    R_gas   = scheme["R_gas"].number;
    M       = scheme["M"].number;
    P_min   = scheme["P_min"].number;
    E_act   = scheme["E_act"].number;
    Z_freq  = scheme["Z_freq"].number;
    VISC    = scheme["VISC"].number;
    MINWT   = scheme["MINWT"].number;
    GASW    = scheme["GASW"].number;
    MINGRHO = scheme["MINGRHO"].number;
    Q_chem  = scheme["Q_chem"].number;

    step_fo  = toml.root["recording"].table["step_fo"].number;
    time_fo  = toml.root["recording"].table["time_fo"].number;
    step_max = toml.root["recording"].table["step_max"].number;
}


void Grid(std::vector<double>& x, std::vector<double>& y,
          int offset_x, int offset_y) {

    double dx = Lx / (Nx_glob - 1);
    double dy = Ly / (Ny_glob - 1);

    for (int i = 0; i < Nx + 2*fict; i++) {
        int global_i = offset_x + i - fict;
        x[i] = global_i * dx;
    }
    for (int j = 0; j < Ny + 2*fict; j++) {
        int global_j = offset_y + j - fict;
        y[j] = global_j * dy;
    }
}


void InitValues(Field& W,
                const std::vector<double>& x,
                const std::vector<double>& y,
                const std::string& config_path) {

    SimpleToml config, test;

    config.load(config_path);
    std::string Test      = config.root["simulation"].table["Test"].str;
    std::string direction = config.root["simulation"].table["direction"].str;

    if (Test == "custom") {
        std::cerr << "Кастомный тест пока не настроен!" << std::endl;
        return;
    }

    // D: x < x_det, y > y_corner
    // Стенка: x < x_corner, y < y_corner
    // ВВ канал: x < x_corner, y > y_corner, x >= x_det
    // ВВ широкая: x >= x_corner 
    if (Test == "Mader3") {

        if (!test.load("tests.toml")) {
            std::cerr << "Нет файла с тестами" << std::endl;
            return;
        }
        t_max = test.root[Test].table["max_t"].number;

        // Геометрические границы
        const double x_corner = 2.0;   // конец канала
        const double y_corner = 1.0;   // высота стенки
        const double x_det    = 0.15;  // ширина детонатора

        // Параметры ВВ
        const double rho_bb   = 1.84;
        const double P_bb     = 0.0;
        const double mf_bb    = 1.0;   // непрореагировавшее

        // Параметры детонатора
		const double u_det    = 0.22;
        const double rho_det  = 1.84;
        const double P_det    = 0.3562; // P_CJ
        const double mf_det   = 0.0;   // уже прореагировало

        // Параметры стенки
        const double rho_wall = 1000.0;
        const double P_wall   = 0.0;
        const double mf_wall  = 0.0;

        size_t Nx_tot = Nx + 2*fict - 1;
        size_t Ny_tot = Ny + 2*fict - 1;

        auto xc = [&](size_t i) { return 0.5 * (x[i] + x[i+1]); };
        auto yc = [&](size_t j) { return 0.5 * (y[j] + y[j+1]); };

        for (size_t i = fict; i < Nx_tot - fict; i++) {
            for (size_t j = fict; j < Ny_tot - fict; j++) {

                double xi = xc(i);
                double yj = yc(j);

                // Стенка
                bool is_wall = (xi < x_corner) && (yj >= y_corner);

                // Детонатор
                bool is_det  = (xi < x_det) && (yj < y_corner);

                if (is_wall) {
                    W[i][j] = {rho_wall, 0.0, 0.0, P_wall};
                    mass_fraction[i][j][0] = mf_wall;
                }
                else if (is_det) {
                    W[i][j] = {rho_det, u_det, 0.0, P_det};
                    mass_fraction[i][j][0] = mf_det;
                }
                else {
                    // ВВ 
                    W[i][j] = {rho_bb, 0.0, 0.0, P_bb};
                    mass_fraction[i][j][0] = mf_bb;
                }
            }
        }

        return;
    }

    if (!test.load("tests.toml")) {
        std::cerr << "Нет файла с тестами" << std::endl;
        return;
    }

    auto& values = test.root[Test].table;
    double rho_L, u_L, P_L;
    double rho_R, u_R, P_R;

    rho_L   = values["rho_L"].number;
    u_L     = values["u_L"].number;
    P_L     = values["P_L"].number;
    Omega_L = values["Omega_L"].number;

    rho_R   = values["rho_R"].number;
    u_R     = values["u_R"].number;
    P_R     = values["P_R"].number;
    Omega_R = values["Omega_R"].number;

    x0    = values["x_gap"].number;
    t_max = values["max_t"].number;

    size_t Nx_tot = Nx + 2*fict - 1;
    size_t Ny_tot = Ny + 2*fict - 1;

    auto cell_center = [&](size_t k) {
        return 0.5 * (x[k] + x[k+1]);
    };

    if (direction == "x") {
        for (size_t i = fict; i < Nx_tot - fict; i++) {
            for (size_t j = fict; j < Ny_tot - fict; j++) {
                double xc = cell_center(i);
                if (xc < x0) {
                    W[i][j] = {rho_L, u_L, 0.0, P_L};
                    mass_fraction[i][j][0] = Omega_L;
                } else {
                    W[i][j] = {rho_R, u_R, 0.0, P_R};
                    mass_fraction[i][j][0] = Omega_R;
                }
            }
        }
    }
    else if (direction == "y") {
        for (size_t i = fict; i < Nx_tot - fict; i++) {
            for (size_t j = fict; j < Ny_tot - fict; j++) {
                double xc = cell_center(i);
                if (xc < x0)
                    W[j][i] = {rho_L, 0.0, u_L, P_L};
                else
                    W[j][i] = {rho_R, 0.0, u_R, P_R};
            }
        }
    }
}