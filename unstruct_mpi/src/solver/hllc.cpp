#include "mesh/mesh.hpp"
#include <cmath>
#include <algorithm>

Vec4 flux_normal(const Vec4& U, double nx, double ny, double gamma) {
    double rho = U[0];
    double u   = U[1] / rho;
    double v   = U[2] / rho;
    double E   = U[3];
    double p   = (gamma - 1.0) * (E - 0.5*rho*(u*u + v*v));
    double un  = u*nx + v*ny;   

    return { rho*un,
             rho*u*un + p*nx,
             rho*v*un + p*ny,
             (E + p)*un };
}

Vec4 U_star(const Vec4& U, double nx, double ny,
                   double SK, double Ss, double gamma) {
    double rho = U[0];
    double u   = U[1] / rho;
    double v   = U[2] / rho;
    double E   = U[3];
    double p   = (gamma - 1.0) * (E - 0.5*rho*(u*u + v*v));
    double un  = u*nx + v*ny;

    double coeff = rho * (SK - un) / (SK - Ss);

    return { coeff,
             coeff * (u  + (Ss - un)*nx),
             coeff * (v  + (Ss - un)*ny),
             coeff * (E/rho + (Ss - un)*(Ss + p/(rho*(SK - un)))) };
}

Vec4 hllc(const Vec4& UL, const Vec4& UR,
                 double nx, double ny, double gamma) {

    double rL  = UL[0];
    double uL  = UL[1] / rL;
    double vL  = UL[2] / rL;
    double EL  = UL[3];
    double pL  = (gamma - 1.0) * (EL - 0.5*rL*(uL*uL + vL*vL));
    double aL  = std::sqrt(gamma * std::max(pL, 1e-12) / rL);
    double unL = uL*nx + vL*ny;
    // Полная энтальпия H = (E + p) / rho — нужна для Roe-среднего
    double HL  = (EL + pL) / rL;

    double rR  = UR[0];
    double uR  = UR[1] / rR;
    double vR  = UR[2] / rR;
    double ER  = UR[3];
    double pR  = (gamma - 1.0) * (ER - 0.5*rR*(uR*uR + vR*vR));
    double aR  = std::sqrt(gamma * std::max(pR, 1e-12) / rR);
    double unR = uR*nx + vR*ny;
    double HR  = (ER + pR) / rR;


    double sqL = std::sqrt(rL);
    double sqR = std::sqrt(rR);
    double inv = 1.0 / (sqL + sqR);

    double ur  = (sqL*uL + sqR*uR) * inv;
    double vr  = (sqL*vL + sqR*vR) * inv;
    double Hr  = (sqL*HL + sqR*HR) * inv;
    double unr = ur*nx + vr*ny;
    double ar  = std::sqrt(std::max(1e-12, (gamma - 1.0)*(Hr - 0.5*(ur*ur + vr*vr))));


    double SL = std::min(unL - aL, unr - ar);
    double SR = std::max(unR + aR, unr + ar);

    double denom = rL*(SL - unL) - rR*(SR - unR);
    double Ss    = (pR - pL + rL*unL*(SL - unL) - rR*unR*(SR - unR))
                 / (denom + 1e-300);  // 1e-300 защита от деления на ноль


    if (SL >= 0.0) {
        return flux_normal(UL, nx, ny, gamma);
    }
    if (SR <= 0.0) {
        return flux_normal(UR, nx, ny, gamma);
    }
    if (Ss >= 0.0) {
        Vec4 FL  = flux_normal(UL, nx, ny, gamma);
        Vec4 ULs = U_star(UL, nx, ny, SL, Ss, gamma);
        return FL + (ULs - UL) * SL;
    } else {
        Vec4 FR  = flux_normal(UR, nx, ny, gamma);
        Vec4 URs = U_star(UR, nx, ny, SR, Ss, gamma);
        return FR + (URs - UR) * SR;
    }
}