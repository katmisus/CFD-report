#include "Poisson.h"
#include <cmath>
#include <algorithm>

// CHECK: POISSON 
void solvePoissonSOR(Grid2D& pp,
                     const Grid2D& u_w,
                     const Grid2D& v_w,
                     const Grid2D& au_w,
                     const Grid2D& av_w,
                     const SimParams& par,
                     int    max_iter,
                     double tol,
                     double omega) {
                        
    const int    Nx = par.Nx;
    const int    Ny = par.Ny;
    const double dx = par.Lx / Nx;
    const double dy = par.Ly / Ny;
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    // Инициализируем pp нулями
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            pp[i][j] = 0.0;

    for (int iter = 0; iter < max_iter; iter++) {
        double res = 0.0;

        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {

                if (i == 0 && j == 0) {
                    pp[0][0] = 0.0;
                    continue;
                }

                // Неоднородность
                double b = -((u_w[i+1][j] - u_w[i][j])*dy
                           + (v_w[i][j+1] - v_w[i][j])*dx);

                // Коэффициенты Пуассона
                double aE = (i < Nx - 1) ? dy2 / au_w[i+1][j] : 0.0;
                double aW = (i > 0)      ? dy2 / au_w[i][j]   : 0.0;
                double aN = (j < Ny - 1) ? dx2 / av_w[i][j+1] : 0.0;
                double aS = (j > 0)      ? dx2 / av_w[i][j]   : 0.0;
                double aP = aE + aW + aN + aS;

                if (aP < 1e-30) continue; 

                // Соседние значения p штрих
                double ppE = (i < Nx - 1) ? pp[i+1][j] : pp[i][j];
                double ppW = (i > 0)      ? pp[i-1][j] : pp[i][j];
                double ppN = (j < Ny - 1) ? pp[i][j+1] : pp[i][j];
                double ppS = (j > 0)      ? pp[i][j-1] : pp[i][j];

                double pp_new = (aE*ppE + aW*ppW + aN*ppN + aS*ppS + b) / aP;
                double delta  = pp_new - pp[i][j];
                pp[i][j]     += omega * delta;

                res = std::max(res, std::abs(delta));
            }
        }

        if (res < tol) break;
    }
}

double massResidual(const Grid2D& u_w,
                    const Grid2D& v_w,
                    const SimParams& par) {
    const int    Nx = par.Nx;
    const int    Ny = par.Ny;
    const double dx = par.Lx / Nx;
    const double dy = par.Ly / Ny;

    double maxR = 0.0;
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j) {
            double r = std::abs((u_w[i+1][j] - u_w[i][j]) * dy
                                + (v_w[i][j+1] - v_w[i][j]) * dx);
            maxR = std::max(maxR, r);
        }
    return maxR;
}
