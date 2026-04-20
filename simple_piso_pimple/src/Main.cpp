#include <iostream>
#include <string>
#include <stdexcept>
#include <filesystem>

#include "Types.h"
#include "Init.h"
#include "PISO.h"
#include "PIMPLE.h"
#include "SIMPLE.h"

namespace fs = std::filesystem;

void runSolver(Fields& f, const SimParams& par,
               const std::string& output_dir) {
    if (par.method == "SIMPLE") 
        runSIMPLE(f, par, output_dir);
    else if (par.method == "PIMPLE") 
        runPIMPLE(f, par, output_dir);
    else 
        runPISO(f, par, output_dir);
}

int main(int argc, char* argv[]) {
    SimParams par;
    std::string output_dir;

    if (argc >= 2) {
        std::string config_path = argv[1];
        par = readParamsTOML(config_path);
        
        fs::path cp(config_path);
        output_dir = "/mnt/d/Study/cfd/output/" + cp.stem().string();
    } else {
        // Задача по умолчанию — каверна Re=100
        output_dir = "/mnt/d/Study/cfd/output/cavity_Re100";
    }

    if (fs::exists(output_dir))
        fs::remove_all(output_dir);
    fs::create_directories(output_dir);

    Fields f;

    // Инициализация полей
    if (par.task == "taylor_green") {
        std::cout << "Task: Taylor-Green Vortices\n";
        f = initTaylorGreen(par);
    } else if (par.task == "step") {
        std::cout << "Task: Backward-Facing Step\n";
        f = initStep(par);
    } else {
        std::cout << "Task: Lid-Driven Cavity\n";
        f = initCavity(par);
    }

    try {
        runSolver(f, par, output_dir);
    } catch (const std::exception& e) {
        std::cerr << "Solver error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
