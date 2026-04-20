#ifndef DECOMPOSITION_HPP
#define DECOMPOSITION_HPP

#include "mesh/mesh.hpp"
#include <vector>

std::vector<int> rcb_partition(const Mesh& mesh, int n_parts);

void print_partition_stats(const Mesh& mesh,
                           const std::vector<int>& part,
                           int n_parts);

#endif 
