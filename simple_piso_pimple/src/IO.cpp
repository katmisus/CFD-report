#include "IO.h"
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <algorithm>

namespace fs = std::filesystem;


static std::string stepFilename(int step) {
    std::ostringstream ss;
    ss << "step_" << std::setw(6) << std::setfill('0') << step << ".vtk";
    return ss.str();
}

void saveFieldsVTK(const Fields& f,
                   const SimParams& par,
                   const std::string& dir,
                   int step,
                   double t) {
    fs::create_directories(dir);

    const int    Nx  = par.Nx;
    const int    Ny  = par.Ny;
    const double dx  = par.Lx / Nx;
    const double dy  = par.Ly / Ny;
    const int    Nxn = Nx + 1;
    const int    Nyn = Ny + 1;
    const int    Npts = Nxn * Nyn;

    // Интерполяция в узлы

    auto u_node = [&](int i, int j) -> double {
        if      (j == 0)  return f.u[i][0];
        else if (j == Ny) return f.u[i][Ny-1];
        else              return 0.5 * (f.u[i][j-1] + f.u[i][j]);
    };

    auto v_node = [&](int i, int j) -> double {
        if      (i == 0)  return f.v[0][j];
        else if (i == Nx) return f.v[Nx-1][j];
        else              return 0.5 * (f.v[i-1][j] + f.v[i][j]);
    };

    // p_node(i,j): билинейная интерполяция из четырёх соседних
    // центров ячеек; на границе — нейман (зеркало)
    auto p_node = [&](int i, int j) -> double {
        int i0 = std::max(0, i-1), i1 = std::min(Nx-1, i);
        int j0 = std::max(0, j-1), j1 = std::min(Ny-1, j);
        return 0.25 * (f.p[i0][j0] + f.p[i1][j0]
                     + f.p[i0][j1] + f.p[i1][j1]);
    };

    // div в центре ячейки (i,j), затем интерполяция в узел
    auto div_cell = [&](int i, int j) -> double {
        return (f.u[i+1][j] - f.u[i][j]) / dx
             + (f.v[i][j+1] - f.v[i][j]) / dy;
    };
    auto div_node = [&](int i, int j) -> double {
        int i0 = std::max(0, i-1), i1 = std::min(Nx-1, i);
        int j0 = std::max(0, j-1), j1 = std::min(Ny-1, j);
        return 0.25 * (div_cell(i0,j0) + div_cell(i1,j0)
                     + div_cell(i0,j1) + div_cell(i1,j1));
    };

    // Запись VTK Legacy ASCII
    std::string fpath = dir + "/" + stepFilename(step);
    std::ofstream out(fpath);
    if (!out.is_open())
        throw std::runtime_error("Cannot write VTK: " + fpath);

    out << std::setprecision(8) << std::scientific;

    // Заголовок
    out << "# vtk DataFile Version 3.0\n";
    out << "PISO  t=" << t << "  step=" << step << "\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_GRID\n";
    out << "DIMENSIONS " << Nxn << " " << Nyn << " 1\n";
    out << "POINTS " << Npts << " double\n";

    // Координаты узлов
    for (int j = 0; j < Nyn; ++j)
        for (int i = 0; i < Nxn; ++i)
            out << (i * dx) << " " << (j * dy) << " 0.0\n";

    // Данные в узлах
    out << "\nPOINT_DATA " << Npts << "\n";

    // Вектор скорости
    out << "VECTORS Velocity double\n";
    for (int j = 0; j < Nyn; ++j)
        for (int i = 0; i < Nxn; ++i)
            out << u_node(i,j) << " " << v_node(i,j) << " 0.0\n";

    // Давление
    out << "\nSCALARS Pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int j = 0; j < Nyn; ++j)
        for (int i = 0; i < Nxn; ++i)
            out << p_node(i,j) << "\n";

    // // Дивергенция
    // out << "\nSCALARS Divergence double 1\n";
    // out << "LOOKUP_TABLE default\n";
    // for (int j = 0; j < Nyn; ++j)
    //     for (int i = 0; i < Nxn; ++i)
    //         out << div_node(i,j) << "\n";

    out.close();
}

void saveFieldsCSV(const Fields& f,
                   const SimParams& par,
                   const std::string& dir,
                   int step,
                   double t) {
    fs::create_directories(dir);

    const int    Nx = par.Nx;
    const int    Ny = par.Ny;
    const double dx = par.Lx / Nx;
    const double dy = par.Ly / Ny;

    std::ostringstream fname;
    fname << dir << "/" << step << "_step.csv";

    std::ofstream out(fname.str());
    if (!out.is_open())
        throw std::runtime_error("Cannot write CSV: " + fname.str());

    // Заголовок
    out << "# t=" << t << "  step=" << step << "\n";
    out << "x,y,u,v,p,div\n";
    out << std::setprecision(10) << std::scientific;

    for (int j = 0; j < Ny; ++j) {
        double y_c = (j + 0.5) * dy;
        for (int i = 0; i < Nx; ++i) {
            double x_c = (i + 0.5) * dx;

            double u_c  = 0.5 * (f.u[i][j]   + f.u[i+1][j]);
            double v_c  = 0.5 * (f.v[i][j]   + f.v[i][j+1]);
            double p_c  = f.p[i][j];
            double div  = (f.u[i+1][j] - f.u[i][j]) / dx
                        + (f.v[i][j+1] - f.v[i][j]) / dy;

            out << x_c << "," << y_c << ","
                << u_c << "," << v_c << ","
                << p_c << "," << div << "\n";
        }
    }
    out.close();
}

