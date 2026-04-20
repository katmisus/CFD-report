#include "physics/init.hpp"
#include "mesh/mesh.hpp"
#include <cmath>

void init_flow(Mesh& mesh, const FlowParams& fp) {
    for (auto& cell : mesh.cells) {
        cell.U[0] = fp.rho_inf;             // ρ
        cell.U[1] = fp.rho_inf * fp.u_inf;  // ρu
        cell.U[2] = fp.rho_inf * fp.v_inf;  // ρv
        cell.U[3] = fp.E_inf;               // E
    }
}