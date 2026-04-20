#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "mesh/mesh.hpp"
#include "physics/init.hpp" 

double compute_dt(const Mesh& mesh, double cfl, const FlowParams& fp);

void residuals(Mesh& mesh, const FlowParams& fp);

void euler_step(Mesh& mesh, double dt);

#endif 
