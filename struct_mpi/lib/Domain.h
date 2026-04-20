#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include <vector>
#include <filesystem>

struct Domain {

    int rank;               // Ранг домена

    int px, py;             // Число процессов по x и y
    int rx, ry;             // "Номер" домена в строке (rx) или столбце (ry)

    int Nx_cells;           // ЛОКАЛЬНОЕ число ячеек домена
    int Ny_cells;

    int Nx;                 // ЛОКАЛЬНОЕ число границ домена
    int Ny;

    int offset_x;           // Сдвиг домена по x и y
    int offset_y;

    std::vector<double> x;  // ЛОКАЛЬНЫЕ координаты x, y (вместе с фиктивными ячейками)
    std::vector<double> y;
};

Domain BuildDomain(int rank,
                   int px, int py,
                   int Nx_glob, 
                   int Ny_glob);
void BuildGrid(Domain& dom, int fict, int Nx_glob, int Ny_glob);
void WriteDomainSizes(const Domain& dom,
                      int p,
                      const std::filesystem::path& base_output);

void WriteDomains(const Domain& dom,
                  const std::vector<double>& x,
                  const std::vector<double>& y,
                  int fict,
                  int p,
                  const std::filesystem::path& base_output);


#endif