#include "PIMPLE.h"
#include "BoundCond.h"
#include "Predictor.h"
#include "Poisson.h"
#include "IO.h"
#include <iostream>
#include <cmath>


void stepPIMPLE(Fields& f, const SimParams& par) {
    const int    Nx  = par.Nx;
    const int    Ny  = par.Ny;
    const double dx  = par.Lx / Nx;
    const double dy  = par.Ly / Ny;


    const Fields f_n = f;

    Grid2D pp(Nx, std::vector<double>(Ny, 0.0));


    // Внешний цикл k = 1..nOuter

    for (int k = 1; k <= par.nOuter; ++k) {

        const bool last_iter = (k == par.nOuter);

        // Релаксация: на последней итерации снимаем (alpha=1)
        const double au = last_iter ? 1.0 : par.alpha_u;
        const double ap = last_iter ? 1.0 : par.alpha_p;


        // Предиктор импульса
        NonStatMomentumPredictor(f, f_n, par, au);
        applyBoundaryConditions(f, par);


        // Цикл PISO-коррекций давления
        for (int m = 0; m < par.nCorr; ++m) {

            solvePoissonSOR(pp, f.u, f.v, f.au, f.av, par);

            // Коррекция давления с релаксацией ap
            for (int i = 0; i < Nx; ++i)
                for (int j = 0; j < Ny; ++j)
                    f.p[i][j] += ap * pp[i][j];

            // Коррекция u: релаксация давления уже в pp*ap,
            for (int i = 1; i < Nx; ++i)
                for (int j = 0; j < Ny; ++j)
                    f.u[i][j] -= (dy / f.au[i][j]) * (pp[i][j] - pp[i-1][j]);

            // Коррекция v
            for (int i = 0; i < Nx; ++i)
                for (int j = 1; j < Ny; ++j)
                    f.v[i][j] -= (dx / f.av[i][j]) * (pp[i][j] - pp[i][j-1]);

            applyBoundaryConditions(f, par);
        }
        // После nCorr коррекций f содержит u^(k), v^(k), p^(k)
    }
    // После nOuter итераций f = u^{n+1}, v^{n+1}, p^{n+1}
}


// CHECK: PIMPLE_LOOP 
void runPIMPLE(Fields& f, const SimParams& par,
               const std::string& output_dir) {
    const double dt     = par.dt;
    const double T      = par.T_end;
    const int    Nsteps = static_cast<int>(T / dt);

    std::cout << "Use PIMPLE solver\n";

    saveFieldsVTK(f, par, output_dir, 0, 0.0);

    for (int step = 1; step <= Nsteps; step++) {

        stepPIMPLE(f, par);

        double t = step * dt;

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
    std::string csv_dir = "CSV/PIMPLE";
    saveFieldsCSV(f, par, csv_dir, Nsteps, par.T_end);
    std::cout << "PIMPLE done. CSV: " << csv_dir << "\n";
}
