#include "mesh/mesh_io.hpp"
#include "mesh/mesh.hpp"
#include "mesh/build_faces.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <map>
#include <cmath>

void seek_section(std::ifstream& f, const std::string& name) {
    std::string line;
    while (std::getline(f, line)) {
        // Убираем \r на случай файла с Windows line endings
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line == name) return;
    }
    throw std::runtime_error("mesh_io: section '" + name + "' not found");
}

Face::BC tag_to_bc(int phys_tag) {
    switch (phys_tag) {
        case 10: return Face::BC::Inflow;
        case 11: return Face::BC::Outflow;
        case 12: return Face::BC::Wall;
        default:
            std::cerr << "mesh_io: unknown phys_tag=" << phys_tag
                      << ", using Wall\n";
            return Face::BC::Wall;
    }
}


Mesh load_gmsh(const std::string& filename) {

    std::ifstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("mesh_io: cannot open '" + filename + "'");

    Mesh mesh;

    seek_section(f, "$Nodes");

    int n_nodes;
    f >> n_nodes;
    mesh.nodes.resize(n_nodes);

    for (int i = 0; i < n_nodes; i++) {
        int id; double x, y, z;
        f >> id >> x >> y >> z;   // id 1-based → индекс id-1
        mesh.nodes[id-1] = {x, y};
    }


    seek_section(f, "$Elements");

    int n_elem;
    f >> n_elem;

    // Слоаврик ребро -> граничное условие
    std::map<std::pair<int,int>, Face::BC> boundary_bc;

    for (int i = 0; i < n_elem; i++) {
        int id, type, ntags;
        f >> id >> type >> ntags;

        std::vector<int> tags(ntags);
        for (int j = 0; j < ntags; j++) f >> tags[j];
        int phys_tag = ntags > 0 ? tags[0] : 0;

        if (type == 1) {
            int a, b;
            f >> a >> b;
            auto key = std::make_pair(std::min(a-1,b-1), std::max(a-1,b-1));
            boundary_bc[key] = tag_to_bc(phys_tag);

        } else if (type == 2) {

            int a, b, c;
            f >> a >> b >> c;

            Cell cell;
            cell.node_ids = {a-1, b-1, c-1};

            // Вычисляем площадь
            const Node& A = mesh.nodes[a-1];
            const Node& B = mesh.nodes[b-1];
            const Node& C = mesh.nodes[c-1];
            double cross = (B.x-A.x)*(C.y-A.y) - (B.y-A.y)*(C.x-A.x);
            cell.vol = 0.5 * std::abs(cross);
            
            // if (cell.vol < 1e-15)
            //     std::cerr << "mesh_io: degenerate cell " << id
            //                 << " (vol=" << cell.vol << ")\n";
            
            // Центроид
            cell.cx = (A.x + B.x + C.x) / 3.0;
            cell.cy = (A.y + B.y + C.y) / 3.0;
            
            cell.U.fill(0.0);
            cell.res.fill(0.0);

            mesh.cells.push_back(std::move(cell));

        } else if (type == 3) {

            int a, b, c, d;
            f >> a >> b >> c >> d;
 
            Cell cell;
            cell.node_ids = {a-1, b-1, c-1, d-1};
 
            
            // Площадь = сумма двух треугольников
            const Node& A = mesh.nodes[a-1];
            const Node& B = mesh.nodes[b-1];
            const Node& C = mesh.nodes[c-1];
            const Node& D = mesh.nodes[d-1];
            double cross1 = (B.x-A.x)*(C.y-A.y) - (B.y-A.y)*(C.x-A.x);
            double cross2 = (C.x-A.x)*(D.y-A.y) - (C.y-A.y)*(D.x-A.x);
            cell.vol = 0.5 * (std::abs(cross1) + std::abs(cross2));
            
            // if (cell.vol < 1e-15)
            //     std::cerr << "mesh_io: degenerate cell " << id
            //                 << " (vol=" << cell.vol << ")\n";

            // Центроид
            cell.cx = (A.x + B.x + C.x + D.x) / 4.0;
            cell.cy = (A.y + B.y + C.y + D.y) / 4.0;
            
            cell.U.fill(0.0);
            cell.res.fill(0.0);
 
            mesh.cells.push_back(std::move(cell));
        }

        else {

            std::string rest;
            std::getline(f, rest);
        }
    }

    build_faces(mesh, boundary_bc);

    // Статистика
    int n_int=0, n_wall=0, n_in=0, n_out=0;
    for (const auto& face : mesh.faces) {
        switch(face.bc) {
            case Face::BC::Interior: n_int++;  break;
            case Face::BC::Wall:     n_wall++; break;
            case Face::BC::Inflow:   n_in++;   break;
            case Face::BC::Outflow:  n_out++;  break;
            default: break;
        }
    }
    std::cout << "Mesh loaded: "
              << mesh.nodes.size() << " nodes, "
              << mesh.cells.size() << " cells, "
              << mesh.faces.size() << " faces\n"
              << "  Interior=" << n_int
              << "  Wall="     << n_wall
              << "  Inflow="   << n_in
              << "  Outflow="  << n_out << "\n";
    
    return mesh;
}
