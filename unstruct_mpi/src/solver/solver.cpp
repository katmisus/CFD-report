#include "mesh/mesh.hpp"
#include "physics/init.hpp" 
#include "physics/bc.hpp" 
#include "solver/hllc.hpp" 
#include "solver/exact_riemann.hpp" 
#include <cmath>

double compute_dt(const Mesh& mesh, double cfl, const FlowParams& fp) {
    double dt = 1e8;

    for (const Face& f : mesh.faces) {

        for (int s = 0; s < 2; s++) {
            int ci = (s == 0) ? f.left : f.right;
            if (ci < 0) continue;  // граничное ребро

            const Vec4& U = mesh.cells[ci].U;
            double rho = U[0];
            double u = U[1] / rho;
            double v = U[2] / rho;
            double p = std::max((fp.gamma - 1.0)
                        * (U[3] - 0.5*rho*(u*u + v*v)), 1e-12);
            double a = std::sqrt(fp.gamma * p / rho);

            double h = mesh.cells[ci].vol / f.length;

            double lam = std::abs(u*f.nx + v*f.ny) + a;

            dt = std::min(dt, cfl * h / (lam + 1e-300));
        }
    }
    return dt;
}

// CHECK: UNSTRUCT_SCHEMES 
void residuals(Mesh& mesh, const FlowParams& fp) {
    mesh.zero_res();

    for (const Face& f : mesh.faces) {
        const Vec4& UL = mesh.cells[f.left].U;

        Vec4 UR = f.is_boundary()
                    ? ghost_state(UL, f.nx, f.ny, f.bc, fp)
                    : mesh.cells[f.right].U;

        Vec4 flux = exact_riemann(UL, UR, f.nx, f.ny, fp.gamma);

        mesh.cells[f.left].res -= flux * f.length;

        if (f.right >= 0)
            mesh.cells[f.right].res += flux * f.length;
    }
}

void euler_step(Mesh& mesh, double dt) {
    for (Cell& c : mesh.cells)
        c.U += c.res * (dt / c.vol);
}