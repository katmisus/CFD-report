#ifndef BUILD_FACES_HPP
#define BUILD_FACES_HPP

#include "mesh/mesh.hpp"
#include <map>

void build_faces(Mesh& mesh,
                 const std::map<std::pair<int,int>, Face::BC>& boundary_bc);

#endif