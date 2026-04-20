#ifndef VTK_HPP
#define VTK_HPP

#include "mesh/mesh.hpp"
#include <string>

void write_vtk(const Mesh& mesh, const std::string& filename,
               double gamma = 1.4);

#endif
