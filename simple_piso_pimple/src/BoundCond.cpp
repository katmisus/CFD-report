#include "BoundCond.h"

void applyBoundaryConditions(Fields& f, const SimParams& par) {
    if (par.task == "taylor_green")
        applyBCTaylorGreen(f, par);
    else if (par.task == "step")
        applyBCStep(f, par);
    else
        applyBCCavity(f, par);
}

void applyBCCavity(Fields& f, const SimParams& par) {
    const int Nx   = par.Nx;
    const int Ny   = par.Ny;

    // u-компонента
    // Левая стенка
    for (int j = 0; j < Ny; ++j)
        f.u[0][j] = 0.0;
    // Правая стенка: i=Nx, все j
    for (int j = 0; j < Ny; ++j)
        f.u[Nx][j] = 0.0;

    // v-компонента
    // Нижняя стенка
    for (int i = 0; i < Nx; ++i)
        f.v[i][0] = 0.0;
    // Верхняя стенка
    for (int i = 0; i < Nx; ++i)
        f.v[i][Ny] = 0.0;

}


void applyBCTaylorGreen(Fields& f, const SimParams& par) {
    const int Nx = par.Nx;
    const int Ny = par.Ny;

    // периодика - u по x
    for (int j = 0; j < Ny; ++j) {
        f.u[0][j]  = f.u[Nx][j];    // левая граница = правая
    }

    // периодика v по y
    for (int i = 0; i < Nx; ++i) {
        f.v[i][0]  = f.v[i][Ny];    // нижняя граница = верхняя
    }
    
}


void applyBCStep(Fields& f, const SimParams& par) {
    const int    Nx      = par.Nx;
    const int    Ny      = par.Ny;
    const int    Nys     = par.Ny_step;         // высота ступеньки в ячейках
    const double dy      = par.Ly / Ny;
    const double Ly_ch   = par.Ly - Nys * dy;  // высота канала выше ступеньки
    const double U_in    = par.U_in;

    // --- Вход: x=0, j >= Nys (выше ступеньки) ---
    // Параболический профиль Poiseuille: u = 6*U_in * (y-y_s)*(y_top-y)/(Ly_ch^2)
    // где y_s = Nys*dy, y_top = Ny*dy
    double y_s   = Nys * dy;
    for (int j = Nys; j < Ny; ++j) {
        double y_face = (j + 0.5) * dy;  // центр ячейки по y
        double xi = (y_face - y_s) / Ly_ch;
        double u_p = 6.0 * U_in * xi * (1.0 - xi);  // пуазейлевский профиль
        f.u[0][j] = u_p;
    }
    // Ниже ступеньки: u=0 на входе (стенка ступеньки)
    for (int j = 0; j < Nys; ++j)
        f.u[0][j] = 0.0;

    // --- Выход: x=Lx — Neumann (du/dx=0) => u[Nx][j] = u[Nx-1][j] ---
    for (int j = 0; j < Ny; ++j)
        f.u[Nx][j] = f.u[Nx-1][j];

    // --- Верхняя стенка: v[i][Ny]=0, u no-slip ---
    for (int i = 0; i < Nx; ++i)
        f.v[i][Ny] = 0.0;

    // --- Нижняя стенка (за ступенькой, i >= Nx_step): v[i][0]=0 ---
    for (int i = 0; i < Nx; ++i)
        f.v[i][0] = 0.0;

    // --- Ступенька: no-slip на поверхности (j = Nys, i < Nx_step) ---
    // Горизонтальная поверхность ступеньки (v=0 на j=Nys для i < Nx_step)
    // В данной геометрии Nx_step=0, ступенька начинается у входа —
    // no-slip реализуется через u=0 на грани i=0 при j < Nys (выше сделано)
    // и v=0 на нижних гранях:
    for (int i = 0; i < Nx; ++i)
        if (i == 0)
            for (int j = 0; j < Nys; ++j)
                f.v[i][j] = 0.0;

    // v на выходе: Neumann
    for (int j = 0; j <= Ny; ++j)
        f.v[Nx-1][j] = f.v[Nx-2][j];
}
