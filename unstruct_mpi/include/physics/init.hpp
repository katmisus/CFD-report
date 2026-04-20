#ifndef INIT_HPP
#define INIT_HPP

#include "mesh/mesh.hpp"
#include <cmath>

struct FlowParams {
    // Входные параметры
    double mach    = 0.3;
    double alpha   = 0.0;
    double gamma   = 1.4;
    double p_inf   = 1.0 / 1.4;
    double rho_inf = 1.0;

    // Производные величины (заполняются через update())
    double a_inf = 0.0;
    double u_inf = 0.0;
    double v_inf = 0.0;
    double E_inf = 0.0;

    void update() {
        a_inf = std::sqrt(gamma * p_inf / rho_inf);
        u_inf = mach * a_inf * std::cos(alpha);
        v_inf = mach * a_inf * std::sin(alpha);
        E_inf = p_inf / (gamma - 1.0)
              + 0.5 * rho_inf * (u_inf*u_inf + v_inf*v_inf);
    }
};

void init_flow(Mesh& mesh, const FlowParams& fp);

#endif

