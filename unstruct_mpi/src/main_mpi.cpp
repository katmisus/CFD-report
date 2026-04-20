#include "mesh/mesh_io.hpp"
#include "physics/init.hpp"
#include "solver/solver.hpp"
#include "io/vtk.hpp"
#include "io/config.hpp"
#include "mpi/mpi_mesh.hpp"
#include "mpi/decomposition.hpp"

#include <mpi.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <filesystem>

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Конфиг и параметры
    std::string cfg_file = (argc > 1) ? argv[1] : "config.txt";
    Config cfg(cfg_file);

    FlowParams fp;
    fp.mach    = cfg.get_double("mach",    0.3);
    fp.alpha   = cfg.get_double("alpha",   0.0) * M_PI / 180.0;
    fp.gamma   = cfg.get_double("gamma",   1.4);
    fp.p_inf   = cfg.get_double("p_inf",   1.0 / 1.4);
    fp.rho_inf = cfg.get_double("rho_inf", 1.0);
    fp.update();

    std::string mesh_file = cfg.get_string("mesh_file", "Meshes/cylinder.msh");
    const double CFL      = cfg.get_double("CFL",    0.5);
    const double T_end    = cfg.get_double("T_end",  10.0);
    const int    io_step  = cfg.get_int   ("io_step", 100);

    // Загрузка глобальной сетки
    Mesh global_mesh;
    if (rank == 0)
        global_mesh = load_gmsh(mesh_file);

    // Разбиение и раздача
    LocalMesh lm = distribute_mesh(global_mesh, MPI_COMM_WORLD);

    int global_nc = 0;
    if (rank == 0) global_nc = global_mesh.nc();
    MPI_Bcast(&global_nc, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> part(global_nc);
    if (rank == 0)
        part = rcb_partition(global_mesh, nprocs);
    MPI_Bcast(part.data(), global_nc, MPI_INT, 0, MPI_COMM_WORLD);

    // Подготовка директории вывода
    namespace fs = std::filesystem;
    std::string mesh_name = fs::path(mesh_file).stem().string();
    std::string out_dir   = cfg.get_string("out_dir", "output") + "/" + mesh_name;

    if (rank == 0) {
        if (fs::exists(out_dir)) fs::remove_all(out_dir);
        fs::create_directories(out_dir);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    auto make_filename = [&](int fr) {
        std::ostringstream ss;
        ss << out_dir << "/flow_" << std::setw(4) << std::setfill('0') << fr << ".vtk";
        return ss.str();
    };

    // Инициализация
    double t = 0.0;
    int step  = 0;
    int frame = 0;

    init_flow(lm.mesh, fp);
    gather_and_write_vtk(lm, global_mesh, part,
                         make_filename(frame++), fp.gamma, MPI_COMM_WORLD);

    // Главный цикл 
    while (t < T_end) {

        // обновляем ghost-ячейки перед вычислением потоков
        exchange_halo(lm, MPI_COMM_WORLD);

        // Локальный dt
        double local_dt = compute_dt(lm.mesh, CFL, fp);
        double dt;
        MPI_Allreduce(&local_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        if (t + dt > T_end) dt = T_end - t;

        // Вычисляем невязки
        residuals(lm.mesh, fp);

        // Шаг Эйлера
        for (int i = 0; i < lm.n_own; i++)
            lm.mesh.cells[i].U += lm.mesh.cells[i].res * (dt / lm.mesh.cells[i].vol);

        // Обнуляем res ghost-ячеек
        for (int i = lm.n_own; i < lm.mesh.nc(); i++)
            lm.mesh.cells[i].res.fill(0.0);

        t += dt;
        step++;

        if (step % io_step == 0 or t == T_end) {
            gather_and_write_vtk(lm, global_mesh, part,
                                 make_filename(frame++), fp.gamma, MPI_COMM_WORLD);
            if (rank == 0)
                std::cout << "t = " << std::fixed << std::setprecision(4) << t
                          << "  step = " << step
                          << "  frame = " << frame - 1 << "\n";
        }
    }

    if (rank == 0)
        std::cout << "Done. " << frame << " frames written to " << out_dir << "/\n";

    MPI_Finalize();
    return 0;
}
