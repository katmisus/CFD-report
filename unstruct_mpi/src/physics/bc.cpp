#include "mesh/mesh.hpp"
#include "physics/init.hpp"   
#include <cmath>
#include <algorithm>

double pressure(const Vec4& U, double gamma) {
    return (gamma - 1.0) * (U[3] - 0.5*(U[1]*U[1] + U[2]*U[2]) / U[0]);
}

double sound_speed(const Vec4& U, double gamma) {
    double p = pressure(U, gamma);
    return std::sqrt(gamma * std::max(p, 1e-12) / std::max(U[0], 1e-12));
}

Vec4 bc_wall(const Vec4& UL, double nx, double ny) {
    double rho = UL[0];
    double u   = UL[1] / rho;
    double v   = UL[2] / rho;
    double un  = u*nx + v*ny;   

    return { rho,
             rho * (u - 2.0*un*nx),   
             rho * (v - 2.0*un*ny),   
             UL[3] };                 
}

Vec4 bc_inflow(const Vec4& /*UL*/, double /*nx*/, double /*ny*/,
                      const FlowParams& fp) {
    double u    = fp.u_inf;
    double v    = fp.v_inf;
    double E    = fp.p_inf / (fp.gamma - 1.0)
                + 0.5 * fp.rho_inf * (u*u + v*v);

    return { fp.rho_inf,
             fp.rho_inf * u,
             fp.rho_inf * v,
             E };
}

Vec4 bc_outflow(const Vec4& UL, double /*nx*/, double /*ny*/,
                       const FlowParams& /*fp*/) {
    return UL;
}

Vec4 bc_symmetry(const Vec4& UL, double nx, double ny,
                        const FlowParams& /*fp*/) {
    return bc_wall(UL, nx, ny);
}

// CHECK: GHOST_STATE 
Vec4 ghost_state(const Vec4& UL, double nx, double ny,
                        Face::BC bc, const FlowParams& fp) {
    switch (bc) {
        case Face::BC::Wall:     return bc_wall    (UL, nx, ny);
        case Face::BC::Inflow:   return bc_inflow  (UL, nx, ny, fp);
        case Face::BC::Outflow:  return bc_outflow (UL, nx, ny, fp);
        case Face::BC::Symmetry: return bc_symmetry(UL, nx, ny, fp);
        default:                 return UL;
    }
}
