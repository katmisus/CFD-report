#include "mpi/mpi_mesh.hpp"
#include "mpi/decomposition.hpp"
#include "mesh/build_faces.hpp"
#include "io/vtk.hpp"   

#include <fstream>
#include <cmath>
#include <mpi.h>
#include <algorithm>
#include <map>
#include <set>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <cassert>


// Передаём геометрию ячейки
struct CellPOD {
    double cx, cy, vol;
    int    n_nodes;          // 3 или 4
    int    node_ids[4];      // глобальные индексы узлов
};

// Передаём грань
struct FacePOD {
    int    left_global;      // глобальный индекс левой ячейки
    int    right_global;     // глобальный индекс правой ячейки
    int    node0, node1;     // глобальные индексы узлов
    double nx, ny, length, mx, my;
    int    bc_int;          
};


static int bc_to_int(Face::BC bc) { return static_cast<int>(bc); }
static Face::BC int_to_bc(int i)  { return static_cast<Face::BC>(i); }


LocalMesh distribute_mesh(const Mesh& global_mesh, MPI_Comm comm) {
    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    // 1. разбиение RCB
    int global_nc = 0;
    std::vector<int> part;

    if (rank == 0) {
        global_nc = global_mesh.nc();
        part = rcb_partition(global_mesh, nprocs);
        if (rank == 0) print_partition_stats(global_mesh, part, nprocs);
    }

    MPI_Bcast(&global_nc, 1, MPI_INT, 0, comm);
    part.resize(global_nc);
    MPI_Bcast(part.data(), global_nc, MPI_INT, 0, comm);

    int global_nn = 0;
    std::vector<double> node_buf; 

    if (rank == 0) {
        global_nn = global_mesh.nn();
        node_buf.resize(2 * global_nn);
        for (int i = 0; i < global_nn; i++) {
            node_buf[2*i]   = global_mesh.nodes[i].x;
            node_buf[2*i+1] = global_mesh.nodes[i].y;
        }
    }

    MPI_Bcast(&global_nn, 1, MPI_INT, 0, comm);
    node_buf.resize(2 * global_nn);
    MPI_Bcast(node_buf.data(), 2*global_nn, MPI_DOUBLE, 0, comm);


    // 2. Собираем массивы ячеек для каждого процесса
    // Определяем, какие ячейки идут какому процессу
    std::vector<std::vector<int>> cells_for(nprocs); 
    if (rank == 0) {
        for (int ci = 0; ci < global_nc; ci++)
            cells_for[part[ci]].push_back(ci);
    }

    // Broadcast cells_for (каждый узнаёт свои глобальные индексы)
    std::vector<int> my_global_cells; // глобальные индексы моих ячеек

    if (rank == 0) {
        for (int r = 0; r < nprocs; r++) {
            int cnt = (int)cells_for[r].size();
            if (r == 0) {
                my_global_cells = cells_for[0];
            } else {
                MPI_Send(&cnt, 1, MPI_INT, r, 0, comm);
                MPI_Send(cells_for[r].data(), cnt, MPI_INT, r, 1, comm);
            }
        }
    } else {
        int cnt;
        MPI_Recv(&cnt, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
        my_global_cells.resize(cnt);
        MPI_Recv(my_global_cells.data(), cnt, MPI_INT, 0, 1, comm, MPI_STATUS_IGNORE);
    }

    const int n_own = (int)my_global_cells.size();

    // Строим обратный словарь: глобальный индекс - локальный
    std::map<int,int> global_to_local; 


    // 3. Необходимые ghost-ячейки
    std::set<int> my_set(my_global_cells.begin(), my_global_cells.end());

    // ghost_cells[r] = список глобальных индексов ghost-ячеек от ранга r
    std::map<int, std::vector<int>> need_from;

    // Аналогично: send_to[remote_rank] = мои глоб. индексы, которые нужны соседу
    // Смотрим все грани глобальной сетки
    if (rank == 0) {
        std::vector<std::vector<std::pair<int,int>>> ghost_needed(nprocs);

        for (const Face& f : global_mesh.faces) {
            if (f.right < 0) continue; // физ. граница
            int pr_l = part[f.left];
            int pr_r = part[f.right];
            if (pr_l == pr_r) continue; // внутри одного процесса

            // Левому процессу нужна правая ячейка как ghost
            ghost_needed[pr_l].push_back({pr_r, f.right});
            // Правому процессу нужна левая ячейка как ghost
            ghost_needed[pr_r].push_back({pr_l, f.left});
        }

        for (int r = 0; r < nprocs; r++) {
            auto& gn = ghost_needed[r];
            std::sort(gn.begin(), gn.end());
            gn.erase(std::unique(gn.begin(), gn.end()), gn.end());

            int sz = (int)gn.size();
            std::vector<int> buf(2*sz);
            for (int i = 0; i < sz; i++) {
                buf[2*i]   = gn[i].first;  
                buf[2*i+1] = gn[i].second; 
            }

            if (r == 0) {
                for (auto& p : gn)
                    need_from[p.first].push_back(p.second);
            } else {
                MPI_Send(&sz, 1, MPI_INT, r, 10, comm);
                if (sz > 0)
                    MPI_Send(buf.data(), 2*sz, MPI_INT, r, 11, comm);
            }
        }
    } else {
        int sz;
        MPI_Recv(&sz, 1, MPI_INT, 0, 10, comm, MPI_STATUS_IGNORE);
        if (sz > 0) {
            std::vector<int> buf(2*sz);
            MPI_Recv(buf.data(), 2*sz, MPI_INT, 0, 11, comm, MPI_STATUS_IGNORE);
            for (int i = 0; i < sz; i++) {
                int owner = buf[2*i];
                int gid   = buf[2*i+1];
                need_from[owner].push_back(gid);
            }
        }
    }


    // 4. Обмен ghost-ячейками

    std::vector<MPI_Request> reqs;
    std::map<int, std::vector<int>> send_request_buf; // что отправляем

    for (auto& [remote, ids] : need_from) {
        int sz = (int)ids.size();
        send_request_buf[remote] = std::vector<int>(ids.begin(), ids.end());

        MPI_Request rq;
        int cnt_msg = sz + 1;
        std::vector<int>& buf = send_request_buf[remote];
        // Prepend size
        buf.insert(buf.begin(), sz);
        MPI_Isend(buf.data(), (int)buf.size(), MPI_INT, remote, 20, comm, &rq);
        reqs.push_back(rq);
    }


    std::map<int, std::vector<int>> send_to; 

    std::set<int> neighbor_set;
    for (auto& p : need_from)
        neighbor_set.insert(p.first);

    int n_neighbors = (int)neighbor_set.size();
    for (int i = 0; i < n_neighbors; i++) {
        MPI_Status st;
        MPI_Probe(MPI_ANY_SOURCE, 20, comm, &st);
        int remote = st.MPI_SOURCE;
        int count;
        MPI_Get_count(&st, MPI_INT, &count);
        std::vector<int> buf(count);
        MPI_Recv(buf.data(), count, MPI_INT, remote, 20, comm, MPI_STATUS_IGNORE);
        int sz = buf[0];
        send_to[remote].assign(buf.begin()+1, buf.begin()+1+sz);
    }

    if (!reqs.empty())
        MPI_Waitall((int)reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    reqs.clear();


    // 5. Отправляем геометрию ghost-ячеек соседям 
    std::vector<CellPOD> all_cell_pods(global_nc);

    if (rank == 0) {
        for (int ci = 0; ci < global_nc; ci++) {
            const Cell& c = global_mesh.cells[ci];
            all_cell_pods[ci].cx  = c.cx;
            all_cell_pods[ci].cy  = c.cy;
            all_cell_pods[ci].vol = c.vol;
            all_cell_pods[ci].n_nodes = (int)c.node_ids.size();
            for (int k = 0; k < 4; k++)
                all_cell_pods[ci].node_ids[k] = (k < (int)c.node_ids.size()) ? c.node_ids[k] : -1;
        }
    }

    MPI_Bcast(all_cell_pods.data(), global_nc * sizeof(CellPOD), MPI_BYTE, 0, comm);


    // 6. Собираем локальную сетку
    LocalMesh lm;
    lm.n_own = n_own;

    // Сначала заполняем global_to_local для своих ячеек
    for (int li = 0; li < n_own; li++)
        global_to_local[my_global_cells[li]] = li;

    // Добавляем узлы (все, как broadcast выше)
    lm.mesh.nodes.resize(global_nn);
    for (int i = 0; i < global_nn; i++) {
        lm.mesh.nodes[i].x = node_buf[2*i];
        lm.mesh.nodes[i].y = node_buf[2*i+1];
    }

    // Добавляем свои ячейки
    lm.mesh.cells.resize(n_own);
    for (int li = 0; li < n_own; li++) {
        int gi = my_global_cells[li];
        const CellPOD& pod = all_cell_pods[gi];
        Cell& c = lm.mesh.cells[li];
        c.cx  = pod.cx;
        c.cy  = pod.cy;
        c.vol = pod.vol;
        c.node_ids.resize(pod.n_nodes);
        for (int k = 0; k < pod.n_nodes; k++)
            c.node_ids[k] = pod.node_ids[k];
        c.U.fill(0.0);
        c.res.fill(0.0);
    }

    // Добавляем ghost-ячейки
    std::map<int,int> ghost_global_to_local; // глоб - локальный
    int ghost_offset = n_own;

    for (auto& [remote, gids] : need_from) {
    
        std::sort(gids.begin(), gids.end());
        gids.erase(std::unique(gids.begin(), gids.end()), gids.end());

        HaloExchange halo;
        halo.remote_rank = remote;

        for (int gi : gids) {
            int li = ghost_offset++;
            global_to_local[gi] = li;
            ghost_global_to_local[gi] = li;

            const CellPOD& pod = all_cell_pods[gi];
            Cell ghost;
            ghost.cx  = pod.cx;
            ghost.cy  = pod.cy;
            ghost.vol = pod.vol;
            ghost.node_ids.resize(pod.n_nodes);
            for (int k = 0; k < pod.n_nodes; k++)
                ghost.node_ids[k] = pod.node_ids[k];
            ghost.U.fill(0.0);
            ghost.res.fill(0.0);
            lm.mesh.cells.push_back(std::move(ghost));

            halo.recv_ids.push_back(li); // куда принимать
        }

        // send_ids: какие мои ячейки отправлять этому соседу
        auto it = send_to.find(remote);
        if (it != send_to.end()) {
            for (int gi : it->second) {
                auto jt = global_to_local.find(gi);
                if (jt != global_to_local.end())
                    halo.send_ids.push_back(jt->second);
            }
        }

        lm.halos.push_back(std::move(halo));
    }


    // 7. Строим локальные грани

    int global_nf = 0;
    std::vector<FacePOD> all_face_pods;

    if (rank == 0) {
        global_nf = global_mesh.nf();
        all_face_pods.resize(global_nf);
        for (int fi = 0; fi < global_nf; fi++) {
            const Face& f = global_mesh.faces[fi];
            all_face_pods[fi] = {
                f.left, f.right,
                f.node0, f.node1,
                f.nx, f.ny, f.length, f.mx, f.my,
                bc_to_int(f.bc)
            };
        }
    }

    MPI_Bcast(&global_nf, 1, MPI_INT, 0, comm);
    all_face_pods.resize(global_nf);
    MPI_Bcast(all_face_pods.data(), global_nf * sizeof(FacePOD), MPI_BYTE, 0, comm);

    lm.mesh.faces.clear();
    for (const FacePOD& fp : all_face_pods) {
        bool left_mine  = my_set.count(fp.left_global) > 0;
        bool right_mine = (fp.right_global >= 0) && my_set.count(fp.right_global) > 0;

        if (!left_mine && !right_mine) continue;

        int gl = fp.left_global;
        int gr = fp.right_global;
        double nx = fp.nx, ny = fp.ny;

        bool need_flip = false;
        if (!left_mine && right_mine) {
   
            std::swap(gl, gr);
            nx = -nx; ny = -ny;
            need_flip = true;
        }

        Face face;
        face.node0  = fp.node0;
        face.node1  = fp.node1;
        face.nx     = nx;
        face.ny     = ny;
        face.length = fp.length;
        face.mx     = fp.mx;
        face.my     = fp.my;

        // left — всегда наша ячейка
        face.left = global_to_local.at(gl);

        if (gr < 0) {
            // Физическая граница
            face.right = -1;
            face.bc    = int_to_bc(fp.bc_int);
        } else if (my_set.count(gr) > 0) {
            face.right = global_to_local.at(gr);
            face.bc    = Face::BC::Interior;
        } else {
            // Правая — ghost
            auto jt = global_to_local.find(gr);
            if (jt == global_to_local.end()) continue;
            face.right       = jt->second;
            face.bc          = Face::BC::MPIBound;
            face.remote_rank = part[gr];
        }

        int fi = (int)lm.mesh.faces.size();
        lm.mesh.faces.push_back(face);

        // Обновляем face_ids ячеек
        lm.mesh.cells[face.left].face_ids.push_back(fi);
        if (face.right >= 0 && face.right < n_own)
            lm.mesh.cells[face.right].face_ids.push_back(fi);
    }

    return lm;
}


// CHECK: HALO_EXCHANGE
void exchange_halo(LocalMesh& lm, MPI_Comm comm) {
    // Пересылаем Vec4 U (4 double) для каждой граничной ячейки.

    const int FIELDS = 4; 

    std::vector<std::vector<double>> send_bufs(lm.halos.size());
    std::vector<std::vector<double>> recv_bufs(lm.halos.size());
    std::vector<MPI_Request> reqs(2 * lm.halos.size());

    // 1. Post Irecv
    for (int h = 0; h < (int)lm.halos.size(); h++) {
        const HaloExchange& halo = lm.halos[h];
        int n_recv = (int)halo.recv_ids.size();
        recv_bufs[h].resize(n_recv * FIELDS);
        MPI_Irecv(recv_bufs[h].data(), n_recv * FIELDS, MPI_DOUBLE,
                  halo.remote_rank, 100, comm, &reqs[2*h]);
    }

    // 2. Pack и Isend
    for (int h = 0; h < (int)lm.halos.size(); h++) {
        const HaloExchange& halo = lm.halos[h];
        int n_send = (int)halo.send_ids.size();
        send_bufs[h].resize(n_send * FIELDS);

        for (int j = 0; j < n_send; j++) {
            const Vec4& U = lm.mesh.cells[halo.send_ids[j]].U;
            for (int k = 0; k < FIELDS; k++)
                send_bufs[h][j*FIELDS + k] = U[k];
        }

        MPI_Isend(send_bufs[h].data(), n_send * FIELDS, MPI_DOUBLE,
                  halo.remote_rank, 100, comm, &reqs[2*h + 1]);
    }

    // 3. Waitall
    if (!reqs.empty())
        MPI_Waitall((int)reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

    // 4. Unpack в ghost-ячейки
    for (int h = 0; h < (int)lm.halos.size(); h++) {
        const HaloExchange& halo = lm.halos[h];
        int n_recv = (int)halo.recv_ids.size();
        for (int j = 0; j < n_recv; j++) {
            Vec4& U = lm.mesh.cells[halo.recv_ids[j]].U;
            for (int k = 0; k < FIELDS; k++)
                U[k] = recv_bufs[h][j*FIELDS + k];
        }
    }
}

// Запись в втк.....
static int vtk_type(int n) {
    if (n == 3) return 5;  // VTK_TRIANGLE
    if (n == 4) return 9;  // VTK_QUAD
    return 7;
}

void gather_and_write_vtk(const LocalMesh& lm,
                           const Mesh& global_mesh,
                           const std::vector<int>& part,
                           const std::string& filename,
                           double gamma, MPI_Comm comm) {
    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    const int FIELDS = 4;  // только U[4] на ячейку

    // ── 1. Пакуем свои n_own ячеек в плоский send-буфер ──────────────────
    const int n_own = lm.n_own;
    std::vector<double> send_buf(n_own * FIELDS);
    for (int i = 0; i < n_own; i++) {
        const Vec4& U = lm.mesh.cells[i].U;
        for (int k = 0; k < FIELDS; k++)
            send_buf[i * FIELDS + k] = U[k];
    }

    // ── 2. Rank 0 собирает количества от каждого ранга ────────────────────
    std::vector<int> counts(nprocs, 0);   // counts[r] = n_own[r] * FIELDS
    int my_count = n_own * FIELDS;
    MPI_Gather(&my_count, 1, MPI_INT,
               counts.data(), 1, MPI_INT,
               0, comm);

    // ── 3. MPI_Gatherv — собираем все данные на rank 0 ───────────────────
    std::vector<int>    displs(nprocs, 0);
    std::vector<double> recv_buf;

    if (rank == 0) {
        for (int r = 1; r < nprocs; r++)
            displs[r] = displs[r-1] + counts[r-1];
        int total = displs[nprocs-1] + counts[nprocs-1];
        recv_buf.resize(total);
    }

    MPI_Gatherv(send_buf.data(), my_count, MPI_DOUBLE,
                recv_buf.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, comm);

    // ── 4. Rank 0 пишет VTK ──────────────────────────────────────────────
    if (rank != 0) return;

    // Раскладываем принятые данные по глобальным индексам.
    // part[gi] = ранг-владелец ячейки gi.
    // Для каждого ранга r: его ячейки идут в recv_buf в том же порядке,
    // в каком они записаны в send_buf (локальные 0..n_own-1).
    // Нам нужен обратный словарь: (ранг r, локальный индекс) → глобальный индекс.
    // Строим его из part[]: cells_of[r][local_idx] = global_idx.

    const int global_nc = global_mesh.nc();
    std::vector<std::vector<int>> cells_of(nprocs);
    for (int gi = 0; gi < global_nc; gi++)
        cells_of[part[gi]].push_back(gi);

    // Теперь заполняем U_global[gi] из recv_buf
    // U_global хранит Vec4 для каждой глобальной ячейки
    std::vector<Vec4> U_global(global_nc);
    for (int r = 0; r < nprocs; r++) {
        int base = displs[r];  // начало данных ранга r в recv_buf
        const auto& gids = cells_of[r];
        for (int li = 0; li < (int)gids.size(); li++) {
            int gi = gids[li];
            for (int k = 0; k < FIELDS; k++)
                U_global[gi][k] = recv_buf[base + li * FIELDS + k];
        }
    }

    // Пишем VTK
    std::ofstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("gather_and_write_vtk: cannot open '" + filename + "'");

    const int nn = global_mesh.nn();

    f << "# vtk DataFile Version 3.0\n";
    f << "CFD MPI flow\n";
    f << "ASCII\n";
    f << "DATASET UNSTRUCTURED_GRID\n";

    // Узлы
    f << "POINTS " << nn << " float\n";
    for (const auto& n : global_mesh.nodes)
        f << n.x << " " << n.y << " 0.0\n";

    // Топология
    int conn_size = 0;
    for (const auto& c : global_mesh.cells)
        conn_size += 1 + (int)c.node_ids.size();

    f << "CELLS " << global_nc << " " << conn_size << "\n";
    for (const auto& c : global_mesh.cells) {
        f << c.node_ids.size();
        for (int id : c.node_ids) f << " " << id;
        f << "\n";
    }

    f << "CELL_TYPES " << global_nc << "\n";
    for (const auto& c : global_mesh.cells)
        f << vtk_type((int)c.node_ids.size()) << "\n";

    // Данные ячеек
    f << "CELL_DATA " << global_nc << "\n";

    // Ранг-владелец (для визуализации разбиения)
    f << "SCALARS partition int 1\nLOOKUP_TABLE default\n";
    for (int gi = 0; gi < global_nc; gi++)
        f << part[gi] << "\n";

    // Плотность
    f << "SCALARS density float 1\nLOOKUP_TABLE default\n";
    for (int gi = 0; gi < global_nc; gi++)
        f << U_global[gi][0] << "\n";

    // Давление
    f << "SCALARS pressure float 1\nLOOKUP_TABLE default\n";
    for (int gi = 0; gi < global_nc; gi++) {
        const Vec4& U = U_global[gi];
        double rho = U[0];
        double p = (rho > 1e-10)
            ? (gamma - 1.0) * (U[3] - 0.5*(U[1]*U[1] + U[2]*U[2])/rho)
            : 0.0;
        f << p << "\n";
    }

    // Число Маха
    f << "SCALARS mach float 1\nLOOKUP_TABLE default\n";
    for (int gi = 0; gi < global_nc; gi++) {
        const Vec4& U = U_global[gi];
        double rho = U[0];
        double mach_val = 0.0;
        if (rho > 1e-10) {
            double u = U[1]/rho, v = U[2]/rho;
            double p = (gamma-1.0)*(U[3] - 0.5*rho*(u*u+v*v));
            double a = std::sqrt(gamma * std::max(p, 1e-12) / rho);
            mach_val = std::sqrt(u*u + v*v) / a;
        }
        f << mach_val << "\n";
    }

    // Вектор скорости
    f << "VECTORS velocity float\n";
    for (int gi = 0; gi < global_nc; gi++) {
        const Vec4& U = U_global[gi];
        double rho = U[0];
        if (rho > 1e-10)
            f << U[1]/rho << " " << U[2]/rho << " 0.0\n";
        else
            f << "0.0 0.0 0.0\n";
    }
}
