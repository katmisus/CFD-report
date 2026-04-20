#include "SIMPLE.h"
#include "Predictor.h"
#include "BoundCond.h"
#include "Poisson.h"
#include "IO.h"
#include <iostream>
#include <cmath>
#include <algorithm>


// CHECK: SIMPLE_LOOP 
void runSIMPLE(Fields& f, const SimParams& par,
               const std::string& output_dir /*, const SimpleParams& spar*/) {

    const int    Nx  = par.Nx;
    const int    Ny  = par.Ny;
    const double dx  = par.Lx / Nx;
    const double dy  = par.Ly / Ny;
    const double ap  = par.alpha_p;

    std::cout << "Use SIMPLE solver\n";

    // Начальное состояние
    saveFieldsVTK(f, par, output_dir, 0, 0.0);

    // p со штрихом
    Grid2D pp(Nx, std::vector<double>(Ny, 0.0));

    int iter = 0;
    double delta = 1.0;

    while (iter < par.max_iter && delta > par.tol) {
        
        // Сохранить старые скорости для проверки сходимости
        Grid2D u_old = f.u;

        // Предиктор импульса
        StatMomentumPredictor(f, par);
        applyBoundaryConditions(f, par);

        // Пуассон для штрихованных p
        solvePoissonSOR(pp, f.u, f.v, f.au, f.av, par);

        // Коррекция давления (с alpha_p)
        // CHECK: CORRECTION
        for (int i = 0; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
                f.p[i][j] += ap * pp[i][j];

        // Коррекция u
        for (int i = 1; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
                f.u[i][j] -= (dy / f.au[i][j]) * (pp[i][j] - pp[i-1][j]);

        // Коррекция v
        for (int i = 0; i < Nx; ++i)
            for (int j = 1; j < Ny; ++j)
                f.v[i][j] -= (dx / f.av[i][j]) * (pp[i][j] - pp[i][j-1]);


        applyBoundaryConditions(f, par);

        // Сходимость
        delta = 0.0;
        for (int i = 1; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
                delta = std::max(delta, std::abs(f.u[i][j] - u_old[i][j]));

        // Запись
        if ((iter % par.step_fo == 0) or (delta <= par.tol)) {
            std::cout << "  iter=" << iter
                      << "  delta=" << delta
                      << "\n";
            saveFieldsVTK(f, par, output_dir, iter, static_cast<double>(iter));
        }
        iter++;
    }

    // CSV для верификации
    std::string csv_dir = "CSV/SIMPLE";
    saveFieldsCSV(f, par, csv_dir, iter, static_cast<double>(iter));
    std::cout << "SIMPLE done. CSV: " << csv_dir << "\n";
}
