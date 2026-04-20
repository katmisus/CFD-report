#ifndef MESH_IO_HPP
#define MESH_IO_HPP

#include "mesh.hpp"
#include <string>
#include <map>

Mesh load_gmsh(const std::string& filename);

#endif
