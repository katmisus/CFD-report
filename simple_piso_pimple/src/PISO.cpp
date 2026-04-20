#include "PISO.h"
#include "BoundCond.h"
#include "Predictor.h"
#include "Poisson.h"
#include "IO.h"
#include <iostream>
#include <cmath>
#include <algorithm>

void stepPISO(Fields& f, const SimParams& par) {
    const int    Nx  = par.Nx;
    const int    Ny  = par.Ny;
    const double dx  = par.Lx / Nx;
    const double dy  = par.Ly / Ny;


    Fields f_old = f;  
    NonStatMomentumPredictor(f, f_old, par, /*alpha_u=*/1.0);
    applyBoundaryConditions(f, par);

    Grid2D pp(Nx, std::vector<double>(Ny, 0.0));   

    // Цикл коррекций  
    for (int m = 0; m < par.nCorr; ++m) {

        // Решаем уравнение Пуассона для pp
        solvePoissonSOR(pp, f.u, f.v, f.au, f.av, par);

        // Корректируем давление 
        for (int i = 0; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
                f.p[i][j] += pp[i][j];

        // Корректируем u 
        for (int i = 1; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j) {
                double pp_E = pp[i][j];
                double pp_W = pp[i-1][j];
                f.u[i][j] -= (dy / f.au[i][j]) * (pp_E - pp_W);
            }

        // Корректируем v 
        for (int i = 0; i < Nx; ++i)
            for (int j = 1; j < Ny; ++j) {
                double pp_N = pp[i][j];
                double pp_S = pp[i][j-1];
                f.v[i][j] -= (dx / f.av[i][j]) * (pp_N - pp_S);
            }
        // Применяем ГУ после каждой коррекции
        applyBoundaryConditions(f, par);
    }
}

// CHECK: PISO_LOOP 
void runPISO(Fields& f, const SimParams& par,
             const std::string& output_dir) {
    const int    Nx   = par.Nx;
    const int    Ny   = par.Ny;
    const double dt   = par.dt;
    const double T    = par.T_end;
    const int    Nsteps = static_cast<int>(T / dt);

    std::cout << "Use PISO solver\n";

    saveFieldsVTK(f, par, output_dir, 0, 0.0);

    for (int step = 1; step <= Nsteps; ++step) {

        stepPISO(f, par);
        double t = step * dt;

        // Запись
        if ((step % par.step_fo == 0) or (step == Nsteps)) {
            double mr = massResidual(f.u, f.v, par);
            std::cout << "  step=" << step
                      << "  t=" << t
                      << "  mass_res=" << mr
                      << "\n";
            saveFieldsVTK(f, par, output_dir, step, t);
        }
    }

    // CSV для верификации
    std::string csv_dir = "CSV/PISO";
    saveFieldsCSV(f, par, csv_dir, Nsteps, par.T_end);
    std::cout << "CSV saved: " << csv_dir << "\n";

    std::cout << "PISO done. Output: " << output_dir << "\n";
}
