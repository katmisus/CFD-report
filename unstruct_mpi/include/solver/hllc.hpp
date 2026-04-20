#ifndef HLLC_HPP
#define HLLC_HPP

#include "mesh/mesh.hpp"

Vec4 flux_normal(const Vec4& U, double nx, double ny, double gamma);

Vec4 U_star(const Vec4& U, double nx, double ny,
                double SK, double Ss, double gamma);

Vec4 hllc(const Vec4& UL, const Vec4& UR,
                double nx, double ny, double gamma);

#endif

