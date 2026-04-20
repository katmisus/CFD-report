#ifndef BC_HPP
#define BC_HPP

#include "mesh/mesh.hpp"
#include "physics/init.hpp"   

double pressure(const Vec4& U, double gamma);

double sound_speed(const Vec4& U, double gamma);

Vec4 bc_wall(const Vec4& UL, double nx, double ny);

Vec4 bc_inflow(const Vec4& /*UL*/, double /*nx*/, double /*ny*/,
                      const FlowParams& fp);

Vec4 bc_outflow(const Vec4& UL, double /*nx*/, double /*ny*/,
                       const FlowParams& /*fp*/);

Vec4 bc_symmetry(const Vec4& UL, double nx, double ny,
                        const FlowParams& /*fp*/);

Vec4 ghost_state(const Vec4& UL, double nx, double ny,
                        Face::BC bc, const FlowParams& fp);                                              

#endif 
