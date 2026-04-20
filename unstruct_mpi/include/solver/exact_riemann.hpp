#ifndef EXACT_RIEMANN_HPP
#define EXACT_RIEMANN_HPP

#include "mesh/mesh.hpp"

Vec4 flux_n(const Vec4& U, double nx, double ny, double g);

double f_wave(double p, double pK, double rhoK, double aK, double g,
                             double& df);

double rho_star(double ps, double pK, double rhoK, double g);

Vec4 exact_riemann(const Vec4& UL, const Vec4& UR,
                           double nx, double ny, double g);                                  


#endif 