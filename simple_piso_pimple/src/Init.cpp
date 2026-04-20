#include "Init.h"
#include "BoundCond.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <cmath>

SimParams readParamsTOML(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Cannot open config: " + path);

    SimParams par;
    double Co_target = 0.4;
    // int Nx = 41;
    // int Ny = 41;

    std::string line;
    while (std::getline(file, line)) {
        // Убрать пробелы
        auto trim = [](std::string& s) {
            size_t a = s.find_first_not_of(" \t\r\n");
            size_t b = s.find_last_not_of(" \t\r\n");
            s = (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
        };
        trim(line);

        // Пропустить пустые строки и комментарии
        if (line.empty() || line[0] == '#' || line[0] == '[')
            continue;

        // Разбить на key = value
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = line.substr(0, eq);
        std::string val = line.substr(eq + 1);
        trim(key); trim(val);

        // Удалить комментарий
        auto hash = val.find('#');
        if (hash != std::string::npos) val = val.substr(0, hash);
        trim(val);
        // Снять кавычки со строковых значений: "PIMPLE" -> PIMPLE
        if (val.size() >= 2 && val.front() == '"' && val.back() == '"')
            val = val.substr(1, val.size() - 2);

        try {
            if      (key == "Nx")       par.Nx            = std::stoi(val);
            else if (key == "Ny")       par.Ny            = std::stoi(val);
            else if (key == "Lx")       par.Lx        = std::stod(val);
            else if (key == "Ly")       par.Ly        = std::stod(val);
            else if (key == "Re")       par.Re        = std::stod(val);
            else if (key == "U_lid")    par.U_lid     = std::stod(val);
            else if (key == "T_end")    par.T_end     = std::stod(val);
            else if (key == "Co")       Co_target     = std::stod(val);
            else if (key == "nCorr")    par.nCorr     = std::stoi(val);
            else if (key == "step_fo")  par.step_fo   = std::stoi(val);
            else if (key == "method")   par.method    = val;
            else if (key == "nOuter")   par.nOuter    = std::stoi(val);
            else if (key == "alpha_u")  par.alpha_u   = std::stod(val);
            else if (key == "alpha_p")  par.alpha_p   = std::stod(val);
            else if (key == "tol")       par.tol       = std::stod(val);
            else if (key == "max_iter")  par.max_iter  = std::stoi(val);
            else if (key == "task")      par.task      = val;
            else if (key == "U_in")      par.U_in      = std::stod(val);
            else if (key == "Nx_step")   par.Nx_step   = std::stoi(val);
            else if (key == "Ny_step")   par.Ny_step   = std::stoi(val);
        } catch (...) {
            std::cerr << "[Init] Warning: cannot parse: " << line << "\n";
        }
    }

    double dx = par.Lx / par.Nx;
        double dy = par.Ly / par.Ny;

    
    if (par.task == "taylor_green") {
        par.nu      = par.U_in * par.Lx / par.Re;
        par.dt = Co_target * std::min(dx, dy) / par.U_in;
    }
    else {
        par.nu = par.U_lid * par.Lx / par.Re;
        par.dt = Co_target * std::min(dx, dy) / par.U_lid;
    }


    // Достроить зависимые параметры (nu, dt из Co)
    // SimParams result = makeCavityParams(Nx, Ny,
    //                                     par.Re,
    //                                     par.U_lid,
    //                                     par.Lx,
    //                                     par.Ly,
    //                                     par.T_end,
    //                                     Co_target,
    //                                     par.nCorr,
    //                                     par.step_fo);

    // PIMPLE-параметры
    // result.method  = par.method;
    // result.nOuter  = par.nOuter;
    // result.alpha_u = par.alpha_u;
    // result.alpha_p = par.alpha_p;
    // result.tol      = par.tol;
    // result.max_iter = par.max_iter;
    // // Расширенные параметры задачи
    // result.task    = par.task;
    // result.U_in    = par.U_in;
    // result.Nx_step = par.Nx_step;
    // result.Ny_step = par.Ny_step;
    return par;
}

// ============================================================
//  Каверна
// ============================================================
SimParams makeCavityParams(int Nx, int Ny,
                            double Re,
                            double U_lid,
                            double Lx,
                            double Ly,
                            double T_end,
                            double Co_target,
                            int    nCorr,
                            int    step_fo) {
    SimParams par;
    par.Nx      = Nx;
    par.Ny      = Ny;
    par.Lx      = Lx;
    par.Ly      = Ly;
    par.Re      = Re;
    par.U_lid   = U_lid;
    par.rho     = 1.0;
    par.nu      = U_lid * Lx / Re;
    par.T_end   = T_end;
    par.nCorr   = nCorr;
    par.step_fo = step_fo;

    // CFL
    double dx = Lx / Nx;
    double dy = Ly / Ny;
    par.dt = Co_target * std::min(dx, dy) / U_lid;

    // std::cout << "[Init] Cavity params:\n"
    //           << "  Grid  : " << Nx << " x " << Ny << "\n"
    //           << "  Re    : " << Re << "\n"
    //           << "  nu    : " << par.nu << "\n"
    //           << "  dx    : " << dx << "  dy: " << dy << "\n"
    //           << "  dt    : " << par.dt << "  (Co=" << Co_target << ")\n"
    //           << "  T_end : " << T_end << "\n";

    return par;
}

// Все поля нулевые и задано граничное условие для стенок и крышки
Fields initCavity(const SimParams& par) {

    Fields f(par);
    applyBCCavity(f, par);

    return f;
}


// ============================================================
//  Вихри Тейлора–Грина
// ============================================================


// Распределение в соответствии с функциями
// CHECK: TG_INIT 
Fields initTaylorGreen(const SimParams& par) {
    Fields f(par);

    const int    Nx = par.Nx;
    const int    Ny = par.Ny;
    const double dx = par.Lx / Nx;
    const double dy = par.Ly / Ny;
    const double pi = M_PI;

    // u[i][j] — на грани x=i*dx, y=(j+0.5)*dy
    for (int i = 0; i <= Nx; ++i)
        for (int j = 0; j < Ny; ++j) {
            double x = i * dx;
            double y = (j + 0.5) * dy;
            f.u[i][j] = -std::cos(pi * x) * std::sin(pi * y);
        }

    // v[i][j] — на грани x=(i+0.5)*dx, y=j*dy
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j <= Ny; ++j) {
            double x = (i + 0.5) * dx;
            double y = j * dy;
            f.v[i][j] =  std::sin(pi * x) * std::cos(pi * y);
        }

    // p = 0 (аналитически p = -(cos(2πx)+cos(2πy))/4, но старт с 0 — норма)
    return f;
}

// Аналитическое решение в момент t — в центрах ячеек
// CHECK: TG_BC  
void taylorGreenExact(const SimParams& par, double t,
                      Grid2D& u_ex, Grid2D& v_ex) {
    const int    Nx  = par.Nx;
    const int    Ny  = par.Ny;
    const double dx  = par.Lx / Nx;
    const double dy  = par.Ly / Ny;
    const double pi  = M_PI;
    const double nu  = par.nu;
    double decay = std::exp(-2.0 * pi * pi * nu * t);

    u_ex = makeGrid(Nx, Ny, 0.0);
    v_ex = makeGrid(Nx, Ny, 0.0);

    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j) {
            double x = (i + 0.5) * dx;
            double y = (j + 0.5) * dy;
            u_ex[i][j] = -std::cos(pi * x) * std::sin(pi * y) * decay;
            v_ex[i][j] =  std::sin(pi * x) * std::cos(pi * y) * decay;
        }
}


// ============================================================
//  Обратная ступенька
// ============================================================

// Поля нулевые — ГУ применятся снаружи
Fields initStep(const SimParams& par) {
    Fields f(par);
    return f;
}
