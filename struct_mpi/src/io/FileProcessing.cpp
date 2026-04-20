#include <filesystem>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "Types.h"

extern int Nx, Ny, fict;
extern double gamm;
extern Field mass_fraction;

void SaveFieldToCSV(const Field& W,
                    const std::vector<double>& x,
                    const std::vector<double>& y,
                    const double& time,
                    const std::string& filename,
                    bool append) {

    size_t Nx_cells = Nx - 1;
    size_t Ny_cells = Ny - 1;

    std::ofstream file;

    if (append)
        file.open(filename, std::ios::app);
    else {
        file.open(filename);
        file << "t,x,y,rho,u,v,P,e,mf\n";  // ← добавлен столбец mf
    }

    if (!file.is_open()) {
        std::cerr << "Невозможно открыть файл" << filename << std::endl;
        return;
    }

    for (size_t i = fict; i < Nx_cells + fict; i++) {

        double xc = 0.5 * (x[i] + x[i + 1]);

        for (size_t j = fict; j < Ny_cells + fict; j++) {
            double yc = 0.5 * (y[j] + y[j + 1]);

            double rho = W[i][j][0];
            double P   = W[i][j][NEQ - 1];
            double e   = P / (rho * (gamm - 1.0));
            double mf  = mass_fraction[i][j][0];  // ← массовая доля

            file << time << ","
                 << xc   << ","
                 << yc   << ","
                 << W[i][j][0] << ","   /* rho */
                 << W[i][j][1] << ","   /* u   */
                 << W[i][j][2] << ","   /* v   */
                 << W[i][j][3] << ","   /* P   */
                 << e           << ","  /* e   */
                 << mf          << "\n";/* mf  */
        }
    }

    file.close();
}

void SaveFluxToCSV(const Field& Flux,
                   const std::vector<double>& x,
                   const std::vector<double>& y,
                   const double& time,
                   const std::string& filename,
                   int dir)
{
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Невозможно открыть файл " << filename << std::endl;
        return;
    }

    file << "t,x,y,rho_flux,momx_flux,momy_flux,E_flux\n";

    size_t Nx_tot = Nx + 2*fict - 1;
    size_t Ny_tot = Ny + 2*fict - 1;

    if (dir == 0) {
        for (size_t i = fict; i < Nx_tot + 1 - fict; i++) {
            double xf = x[i];
            for (size_t j = fict; j < Ny_tot - fict; j++) {
                double yc = 0.5 * (y[j] + y[j+1]);
                file << time << "," << xf << "," << yc << ","
                     << Flux[i][j][0] << "," << Flux[i][j][1] << ","
                     << Flux[i][j][2] << "," << Flux[i][j][3] << "\n";
            }
        }
    }
    else if (dir == 1) {
        for (size_t i = fict; i < Nx_tot - fict; i++) {
            double xc = 0.5 * (x[i] + x[i+1]);
            for (size_t j = fict; j < Ny_tot + 1 - fict; j++) {
                double yg = y[j];
                file << time << "," << xc << "," << yg << ","
                     << Flux[i][j][0] << "," << Flux[i][j][1] << ","
                     << Flux[i][j][2] << "," << Flux[i][j][3] << "\n";
            }
        }
    }

    file.close();
}

fs::path CreateDirFromPath(const std::string& file_path) {
    try {
        fs::path input_path(file_path);
        fs::path dir_path = input_path.parent_path() / input_path.stem();
        if (!fs::exists(dir_path))
            fs::create_directories(dir_path);
        return dir_path;
    }
    catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        return fs::path();
    }
}
