#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <string>

using Grid2D = std::vector<std::vector<double>>;

// Создать двумерный массив
inline Grid2D makeGrid(int rows, int cols, double val = 0.0) {
    return Grid2D(rows, std::vector<double>(cols, val));
}

// Параметры задачи
struct SimParams {

    std::string method = "PISO"; // "SIMPLE", "PISO" или "PIMPLE"
    
    int    Nx      = 41;     // число ячеек по x
    int    Ny      = 41;     // число ячеек по y
    double Lx      = 1.0;    // длина области по x
    double Ly      = 1.0;    // длина области по y
    double Re      = 100.0;  // число Рейнольдса
    double rho     = 1.0;    // плотность (const для несжимаемого)
    double nu      = 0.0;    // кинематическая вязкость = U*Lx/Re
    double dt      = 0.01;   // шаг по времени
    double T_end   = 15.0;   // конечное время
    int    step_fo = 100;    // частота записи

    // SIMPLE - параметры
    double tol = 1e-6;     // порог сходимости
    int max_iter = 5000;   // максимум итераций

    // PISO - параметры
    int nCorr = 2;  // число PISO-коррекций давления внутри шага

    // PIMPLE - параметры
    int nOuter = 3;           // число внешних SIMPLE-итераций (PIMPLE)
    double alpha_u  = 0.7;    // релаксация скорости
    double alpha_p  = 0.3;    // релаксация давления

    // Тип задачи
    std::string task = "cavity";  // "cavity", "taylor_green" или "step"

    // cavity
    double U_lid   = 1.0;    // скорость крышки

    // step
    int    Nx_step  = 0;     // ячеек ступеньки по x
    int    Ny_step  = 0;     // ячеек ступеньки по y 
    double U_in     = 1.0;   // максимальная скорость на входе 
};

// Поля скоростей и давления
// CHECK: STAGGERED_GRID
struct Fields {
    Grid2D u;   // (Nx+1) x Ny
    Grid2D v;   // Nx x (Ny+1)
    Grid2D p;   // Nx x Ny

    // диагональные коэффициенты
    Grid2D au;  // (Nx+1) x Ny
    Grid2D av;  // Nx x (Ny+1)

    Fields() = default;
    Fields(const SimParams& par) {
        int Nx = par.Nx, Ny = par.Ny;
        u  = makeGrid(Nx + 1, Ny, 0.0);
        v  = makeGrid(Nx, Ny + 1, 0.0);
        p  = makeGrid(Nx, Ny, 0.0);
        au = makeGrid(Nx + 1, Ny, 1.0);
        av = makeGrid(Nx, Ny + 1, 1.0);
    }
};

#endif
