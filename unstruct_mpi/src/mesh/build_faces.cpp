#include "mesh/mesh.hpp"
#include <map>
#include <cmath>
#include <stdexcept>
#include <iostream>

// CHECK: BUILD_FACES 
void build_faces(Mesh& mesh,
                 const std::map<std::pair<int,int>, Face::BC>& boundary_bc) {

    // Словарик ребро -> индексы ячеек-владельцев
    using EdgeKey = std::pair<int,int>;
    std::map<EdgeKey, std::vector<int>> edge_map;

    for (int ci = 0; ci < (int)mesh.cells.size(); ci++) {
        const auto& nd = mesh.cells[ci].node_ids;
        int nv = (int)nd.size();
        for (int k = 0; k < nv; k++) {
            int a = nd[k];
            int b = nd[(k+1) % nv];
            EdgeKey key = {std::min(a,b), std::max(a,b)};
            edge_map[key].push_back(ci);
        }
    }

    // Создаем ребра
    mesh.faces.clear();
    mesh.faces.reserve(edge_map.size());

    for (auto& [key, cells_list] : edge_map) {
        Face f;
        f.node0 = key.first;
        f.node1 = key.second;

        f.left  = cells_list[0];
        f.right = (cells_list.size() == 2) ? cells_list[1] : -1;


        // Геометрия ребра
        const Node& N0 = mesh.nodes[f.node0];
        const Node& N1 = mesh.nodes[f.node1];

        double dx = N1.x - N0.x;
        double dy = N1.y - N0.y;

        // Длина
        f.length = std::sqrt(dx*dx + dy*dy);

        // Центр
        f.mx = 0.5*(N0.x + N1.x);
        f.my = 0.5*(N0.y + N1.y);

        // Нормаль
        double nx_c = -dy / f.length;
        double ny_c =  dx / f.length;

        // Ориентируем ее L -> R
        const Cell& lc = mesh.cells[f.left];
        double ref_x, ref_y;
        if (f.right >= 0) {
            const Cell& rc = mesh.cells[f.right];
            ref_x = rc.cx - lc.cx;
            ref_y = rc.cy - lc.cy;
        } else {
            ref_x = f.mx - lc.cx;
            ref_y = f.my - lc.cy;
        }

        if (nx_c*ref_x + ny_c*ref_y < 0.0) {
            nx_c = -nx_c;
            ny_c = -ny_c;
        }

        f.nx = nx_c;
        f.ny = ny_c;

        // Граничное условие
        if (f.is_boundary()) {
            auto it = boundary_bc.find(key);
            if (it != boundary_bc.end()) {
                f.bc = it->second;
            } else {
                std::cerr << "build_faces: WARNING — boundary edge ("
                          << f.node0 << ", " << f.node1
                          << ") not found in boundary_bc dict!"
                          << " Defaulting to Wall.\n";
                f.bc = Face::BC::Wall;
            }
        } else {
            f.bc = Face::BC::Interior;
        }

        int fi = (int)mesh.faces.size();
        mesh.faces.push_back(f);

        mesh.cells[f.left].face_ids.push_back(fi);
        if (f.right >= 0)
            mesh.cells[f.right].face_ids.push_back(fi);
    }
}
