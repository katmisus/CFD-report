#include "Predictor.h"
#include <algorithm>
#include <cmath>

// Upwind-схема первого порядка: выбор значения по знаку потока
static inline double upwind(double phi_L, double phi_R, double F) {
    return (F >= 0.0) ? phi_L : phi_R;
}

void NonStatMomentumPredictor(Fields& f,
                       const Fields& f_old,
                       const SimParams& par,
                       double alpha_u) {
    const int    Nx   = par.Nx;
    const int    Ny   = par.Ny;
    const double dx   = par.Lx / Nx;
    const double dy   = par.Ly / Ny;
    const double nu   = par.nu;
    const double rho  = par.rho;
    const double dt   = par.dt;
    const double U    = par.U_lid;
    const bool   periodic = (par.task == "taylor_green");

    const double dxdy = dx * dy;
    const double time_coeff = rho * dxdy / dt;  // член по времени в диагонали

    // Вспомогательные функции для периодических/непериодических индексов
    auto u_nb = [&](int i, int j) -> double {
        // доступ к f_old.u с периодикой по j (для TG)
        if (periodic) {
            int jp = (j + Ny) % Ny;
            return f_old.u[i][jp];
        }
        return f_old.u[i][j];
    };
    auto v_nb = [&](int i, int j) -> double {
        if (periodic) {
            int ip = (i + Nx) % Nx;
            int jp = (j + Ny + 1) % (Ny + 1);
            return f_old.v[ip][jp];
        }
        return f_old.v[i][j];
    };

    // u-предиктор
    for (int i = 1; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            const double u0 = f_old.u[i][j];

            double ue = 0.5 * (f_old.u[i][j] + f_old.u[i+1][j]);
            double uw = 0.5 * (f_old.u[i][j] + f_old.u[i-1][j]);
            double vn = 0.5 * (v_nb(i-1, j+1) + v_nb(i, j+1));
            double vs = 0.5 * (v_nb(i-1, j)   + v_nb(i, j));

            double Fe = ue * dy;
            double Fw = uw * dy;
            double Fn = vn * dx;
            double Fs = vs * dx;

            double aE_c = std::max(-Fe, 0.0);
            double aW_c = std::max( Fw, 0.0);
            double aN_c = std::max(-Fn, 0.0);
            double aS_c = std::max( Fs, 0.0);

            double aE_d = nu * dy / dx;
            double aW_d = nu * dy / dx;
            double aN_d = nu * dx / dy;
            double aS_d = nu * dx / dy;

            double aE = aE_c + aE_d;
            double aW = aW_c + aW_d;
            double aN = aN_c + aN_d;
            double aS = aS_c + aS_d;

            double u_S_val, u_N_val;
            double aS_wall = 0.0, aN_wall = 0.0;

            if (periodic) {
                // Периодические ГУ: соседи через границу
                u_S_val = u_nb(i, (j - 1 + Ny) % Ny);
                u_N_val = u_nb(i, (j + 1) % Ny);
            } else {
                if (j == 0) {
                    u_S_val  = 0.0;
                    aS_wall  = 2.0 * nu * dx / dy;
                    aS       = 0.0;
                } else {
                    u_S_val = f_old.u[i][j-1];
                }

                if (j == Ny - 1) {
                    u_N_val  = U;
                    aN_wall  = 2.0 * nu * dx / dy;
                    aN       = 0.0;
                } else {
                    u_N_val = f_old.u[i][j+1];
                }
            }

            // Диагональ с временным членом
            double aP = (aE + aW + aN + aS + aN_wall + aS_wall
                        + time_coeff) / alpha_u;

            // Оператор H
            double H = aE * f_old.u[i+1][j]
                     + aW * f_old.u[i-1][j]
                     + aN * u_N_val
                     + aS * u_S_val
                     + aN_wall * U
                     + aS_wall * 0.0
                     + time_coeff * u0;

            // Градиент давления
            double dp = (f.p[i][j] - f.p[i-1][j]) * dy;

            f.u[i][j] = (H - dp) / aP
                      + (1.0 - alpha_u) * u0 / alpha_u;

            f.au[i][j] = aP;
        }
    }

    // v-предиктор
    int j_v_start = periodic ? 0 : 1;
    int j_v_end   = periodic ? Ny : Ny;

    for (int i = 0; i < Nx; ++i) {
        for (int j = j_v_start; j < j_v_end; ++j) {
            const double v0 = f_old.v[i][j];

            double ue, uw, vn, vs;
            if (periodic) {
                int jm1 = (j - 1 + Ny) % Ny;
                ue = 0.5 * (u_nb(i+1, jm1) + u_nb(i+1, j));
                uw = 0.5 * (u_nb(i,   jm1) + u_nb(i,   j));
                vn = 0.5 * (f_old.v[i][j] + f_old.v[i][(j+1) % (Ny+1)]);
                vs = 0.5 * (f_old.v[i][j] + f_old.v[i][(j-1+Ny+1) % (Ny+1)]);
            } else {
                ue = 0.5 * (f_old.u[i+1][j-1] + f_old.u[i+1][j]);
                uw = 0.5 * (f_old.u[i][j-1]   + f_old.u[i][j]);
                vn = 0.5 * (f_old.v[i][j]   + f_old.v[i][j+1]);
                vs = 0.5 * (f_old.v[i][j]   + f_old.v[i][j-1]);
            }

            double Fe = ue * dy;
            double Fw = uw * dy;
            double Fn = vn * dx;
            double Fs = vs * dx;

            double aE_c = std::max(-Fe, 0.0);
            double aW_c = std::max( Fw, 0.0);
            double aN_c = std::max(-Fn, 0.0);
            double aS_c = std::max( Fs, 0.0);

            double aE_d = nu * dy / dx;
            double aW_d = nu * dy / dx;
            double aN_d = nu * dx / dy;
            double aS_d = nu * dx / dy;

            double aE = aE_c + aE_d;
            double aW = aW_c + aW_d;
            double aN = aN_c + aN_d;
            double aS = aS_c + aS_d;

            double v_W_val, v_E_val, v_N_val, v_S_val;
            double aW_wall = 0.0, aE_wall = 0.0;

            if (periodic) {
                v_W_val = f_old.v[(i - 1 + Nx) % Nx][j];
                v_E_val = f_old.v[(i + 1) % Nx][j];
                v_N_val = f_old.v[i][(j + 1) % (Ny + 1)];
                v_S_val = f_old.v[i][(j - 1 + Ny + 1) % (Ny + 1)];
            } else {
                v_N_val = f_old.v[i][j+1];
                v_S_val = f_old.v[i][j-1];
                if (i == 0) {
                    v_W_val = 0.0;
                    aW_wall = 2.0 * nu * dy / dx;
                    aW      = 0.0;
                } else {
                    v_W_val = f_old.v[i-1][j];
                }
                if (i == Nx - 1) {
                    v_E_val = 0.0;
                    aE_wall = 2.0 * nu * dy / dx;
                    aE      = 0.0;
                } else {
                    v_E_val = f_old.v[i+1][j];
                }
            }

            double aP = (aE + aW + aN + aS + aE_wall + aW_wall
                        + time_coeff) / alpha_u;

            if (aP < 1e-30) continue;

            double dp_j  = periodic ? (j % Ny) : j;
            double dp_jm = periodic ? ((j - 1 + Ny) % Ny) : (j - 1);
            double dp = (f.p[i][dp_j] - f.p[i][dp_jm]) * dx;

            double H = aE * v_E_val
                     + aW * v_W_val
                     + aN * v_N_val
                     + aS * v_S_val
                     + aE_wall * 0.0
                     + aW_wall * 0.0
                     + time_coeff * v0;

            f.v[i][j] = (H - dp) / aP
                      + (1.0 - alpha_u) * v0 / alpha_u;

            f.av[i][j] = aP;
        }
    }
}

// Стационарный предиктор импульса
void StatMomentumPredictor(Fields& f, const SimParams& par) {

    const int    Nx  = par.Nx;
    const int    Ny  = par.Ny;
    const double dx  = par.Lx / Nx;
    const double dy  = par.Ly / Ny;
    const double nu  = par.nu;
    const double au  = par.alpha_u;
    const double U   = par.U_lid;

    // u-предиктор
    // CHECK: PREDICTOR_U 
    for (int i = 1; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            const double u_old = f.u[i][j];

            double ue = 0.5 * (f.u[i][j] + f.u[i+1][j]);
            double uw = 0.5 * (f.u[i][j] + f.u[i-1][j]);
            double vn = 0.5 * (f.v[i-1][j+1] + f.v[i][j+1]);
            double vs = 0.5 * (f.v[i-1][j]   + f.v[i][j]);

            double Fe = ue * dy;
            double Fw = uw * dy;
            double Fn = vn * dx;
            double Fs = vs * dx;

            // Конвективные + диффузионные коэффициенты
            double aE = std::max(-Fe, 0.0) + nu * dy / dx;
            double aW = std::max( Fw, 0.0) + nu * dy / dx;
            double aN = std::max(-Fn, 0.0) + nu * dx / dy;
            double aS = std::max( Fs, 0.0) + nu * dx / dy;

            // Вклад от границ
            double u_S_val = 0.0, u_N_val = 0.0;
            double aS_wall = 0.0, aN_wall = 0.0;

            if (j == 0) {
                aS_wall = 2.0 * nu * dx / dy;
                u_S_val = 0.0;
                aS      = 0.0;
            } else {
                u_S_val = f.u[i][j-1];
            }

            if (j == Ny - 1) {
                // Верхняя крышка: u_wall = U_lid
                aN_wall = 2.0 * nu * dx / dy;
                u_N_val = U;
                aN      = 0.0;
            } else {
                u_N_val = f.u[i][j+1];
            }

            // Даигональный член
            double aP = aE + aW + aN + aS + aN_wall + aS_wall;

            if (aP < 1e-30) continue;

            // Оператор H
            double H = aE * f.u[i+1][j]
                     + aW * f.u[i-1][j]
                     + aN * u_N_val
                     + aS * u_S_val
                     + aN_wall * U
                     + aS_wall * 0.0;

            // Градиент давления
            double dp = (f.p[i][j] - f.p[i-1][j]) * dy;

            // Скорость с под-релаксацией
            f.u[i][j] = (1.0 - au) * u_old
                      + au / aP * (H - dp);

            // Эффективная диагональ a_tilde = aP / alpha_u
            f.au[i][j] = aP / au;
        }
    }

    // v-предиктор
    // CHECK: PREDICTOR_V
    for (int i = 0; i < Nx; ++i) {
        for (int j = 1; j < Ny; ++j) {
            const double v_old = f.v[i][j];

            double ue = 0.5 * (f.u[i+1][j-1] + f.u[i+1][j]);
            double uw = 0.5 * (f.u[i][j-1]   + f.u[i][j]);
            double vn = 0.5 * (f.v[i][j]   + f.v[i][j+1]);
            double vs = 0.5 * (f.v[i][j]   + f.v[i][j-1]);

            double Fe = ue * dy;
            double Fw = uw * dy;
            double Fn = vn * dx;
            double Fs = vs * dx;

            double aE = std::max(-Fe, 0.0) + nu * dy / dx;
            double aW = std::max(Fw, 0.0) + nu * dy / dx;
            double aN = std::max(-Fn, 0.0) + nu * dx / dy;
            double aS = std::max(Fs, 0.0) + nu * dx / dy;

            double v_W_val = 0.0, v_E_val = 0.0;
            double aW_wall = 0.0, aE_wall = 0.0;

            if (i == 0) {
                aW_wall = 2.0 * nu * dy / dx;
                v_W_val = 0.0;
                aW      = 0.0;
            } else {
                v_W_val = f.v[i-1][j];
            }

            if (i == Nx - 1) {
                aE_wall = 2.0 * nu * dy / dx;
                v_E_val = 0.0;
                aE      = 0.0;
            } else {
                v_E_val = f.v[i+1][j];
            }

            double aP = aE + aW + aN + aS + aE_wall + aW_wall;

            if (aP < 1e-30) continue;

            double H = aE * v_E_val
                     + aW * v_W_val
                     + aN * f.v[i][j+1]
                     + aS * f.v[i][j-1];

            double dp = (f.p[i][j] - f.p[i][j-1]) * dx;

            f.v[i][j] = (1.0 - au) * v_old
                      + au / aP * (H - dp);

            f.av[i][j] = aP / au;
        }
    }
}