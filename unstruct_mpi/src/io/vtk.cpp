#include "io/vtk.hpp"
#include <fstream>
#include <cmath>
#include <stdexcept>

// VTK cell type по числу узлов
static int vtk_cell_type(int n_nodes) {
    switch (n_nodes) {
        case 3: return 5;   // VTK_TRIANGLE
        case 4: return 9;   // VTK_QUAD
        default: return 7;  // VTK_POLYGON 
    }
}

void write_vtk(const Mesh& mesh, const std::string& filename, double gamma) {
    std::ofstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("vtk: cannot open '" + filename + "'");

    int nc = mesh.nc();
    int nn = mesh.nn();

    // ── Header ────────────────────────────────────────────────────────────────
    f << "# vtk DataFile Version 3.0\n";
    f << "CFD flow\n";
    f << "ASCII\n";
    f << "DATASET UNSTRUCTURED_GRID\n";

    // ── Координаты узлов ──────────────────────────────────────────────────────
    f << "POINTS " << nn << " float\n";
    for (const auto& n : mesh.nodes)
        f << n.x << " " << n.y << " 0.0\n";

    // ── Топология ячеек ───────────────────────────────────────────────────────
    // Подсчитываем суммарный размер connectivity (для каждой ячейки: 1 + n_nodes)
    int connectivity_size = 0;
    for (const auto& c : mesh.cells)
        connectivity_size += 1 + (int)c.node_ids.size();

    f << "CELLS " << nc << " " << connectivity_size << "\n";
    for (const auto& c : mesh.cells) {
        f << c.node_ids.size();
        for (int id : c.node_ids) f << " " << id;
        f << "\n";
    }

    // ── Типы ячеек ────────────────────────────────────────────────────────────
    f << "CELL_TYPES " << nc << "\n";
    for (const auto& c : mesh.cells)
        f << vtk_cell_type((int)c.node_ids.size()) << "\n";

    // ── Физические поля ───────────────────────────────────────────────────────
    f << "CELL_DATA " << nc << "\n";

    // Плотность
    f << "SCALARS density float 1\n";
    f << "LOOKUP_TABLE default\n";
    for (const auto& c : mesh.cells)
        f << c.U[0] << "\n";

    // Давление: p = (γ-1)*(E - (ρu²+ρv²)/(2ρ))
    f << "SCALARS pressure float 1\n";
    f << "LOOKUP_TABLE default\n";
    for (const auto& c : mesh.cells) {
        double rho = c.U[0];
        double p = (rho > 1e-10)
            ? (gamma - 1.0) * (c.U[3] - 0.5*(c.U[1]*c.U[1] + c.U[2]*c.U[2])/rho)
            : 0.0;
        f << p << "\n";
    }

    // Число Маха
    f << "SCALARS mach float 1\n";
    f << "LOOKUP_TABLE default\n";
    for (const auto& c : mesh.cells) {
        double rho = c.U[0];
        double mach_val = 0.0;
        if (rho > 1e-10) {
            double u = c.U[1]/rho, v = c.U[2]/rho;
            double p = (gamma-1.0)*(c.U[3] - 0.5*rho*(u*u+v*v));
            double a = std::sqrt(gamma * std::abs(p) / rho);
            mach_val = std::sqrt(u*u + v*v) / (a > 1e-10 ? a : 1e-10);
        }
        f << mach_val << "\n";
    }

    // Вектор скорости
    f << "VECTORS velocity float\n";
    for (const auto& c : mesh.cells) {
        double rho = c.U[0];
        if (rho > 1e-10)
            f << c.U[1]/rho << " " << c.U[2]/rho << " 0.0\n";
        else
            f << "0.0 0.0 0.0\n";
    }
}