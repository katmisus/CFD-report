#include "mpi/decomposition.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>

// ────────────────────────────────────────────────────────────────────────────
// Внутренняя рекурсивная функция RCB
//
//   cell_ids  — индексы ячеек, которые надо разбить на данном уровне
//   part      — выходной массив
//   cells     — ссылка на вектор ячеек сетки (для центроидов)
//   rank_base — первый номер процесса в текущей группе
//   n_procs   — сколько процессов надо покрыть
// ────────────────────────────────────────────────────────────────────────────
// CHECK: DECOMPOSITION 
void rcb_recursive(std::vector<int>& cell_ids,
                   std::vector<int>& part,
                   const std::vector<Cell>& cells,
                   int rank_base, int n_procs) {
    // Базовый случай: все ячейки принадлежат одному процессу
    if (n_procs == 1) {
        for (int ci : cell_ids)
            part[ci] = rank_base;
        return;
    }
    if (cell_ids.empty()) return;

    const int N = static_cast<int>(cell_ids.size());

    // 1. размах по X и Y 
    double xmin =  1e300, xmax = -1e300;
    double ymin =  1e300, ymax = -1e300;

    for (int ci : cell_ids) {
        double cx = cells[ci].cx;
        double cy = cells[ci].cy;
        if (cx < xmin) xmin = cx;
        if (cx > xmax) xmax = cx;
        if (cy < ymin) ymin = cy;
        if (cy > ymax) ymax = cy;
    }

    const double span_x = xmax - xmin;
    const double span_y = ymax - ymin;

    // 2. выбираем ось
    const bool split_by_x = (span_x >= span_y);

    // 3. сортируем по выбранной оси
    if (split_by_x) {
        std::sort(cell_ids.begin(), cell_ids.end(),
                  [&](int a, int b){ return cells[a].cx < cells[b].cx; });
    } else {
        std::sort(cell_ids.begin(), cell_ids.end(),
                  [&](int a, int b){ return cells[a].cy < cells[b].cy; });
    }

    // 4. точка разреза пропорциональна числу процессов
    const int K_left  = n_procs / 2;
    const int K_right = n_procs - K_left;

    const int cut = static_cast<int>(std::round(static_cast<double>(N) * K_left / n_procs));

    // Гарантируем, что обе части не пустые
    const int cut_safe = std::max(1, std::min(cut, N - 1));

    // 5. рекурсия
    std::vector<int> left_ids (cell_ids.begin(), cell_ids.begin() + cut_safe);
    std::vector<int> right_ids(cell_ids.begin() + cut_safe, cell_ids.end());

    rcb_recursive(left_ids,  part, cells, rank_base, K_left);
    rcb_recursive(right_ids, part, cells, rank_base + K_left,  K_right);
}


// Публичный интерфейс
std::vector<int> rcb_partition(const Mesh& mesh, int n_parts) {

    const int N = mesh.nc();
    std::vector<int> part(N, 0);

    if (n_parts == 1) 
        return part;

    // Начальный список — все ячейки
    std::vector<int> all_ids(N);
    std::iota(all_ids.begin(), all_ids.end(), 0);

    rcb_recursive(all_ids, part, mesh.cells, /*rank_base=*/0, n_parts);

    return part;
}


// Статистика разбиения
void print_partition_stats(const Mesh& mesh,
                            const std::vector<int>& part,
                            int n_parts) {
    const int N = mesh.nc();

    // Число ячеек на каждый процесс
    std::vector<int> count(n_parts, 0);
    for (int ci = 0; ci < N; ci++)
        count[part[ci]]++;

    // Число разрезанных граней (обе ячейки существуют, но в разных частях)
    int cut_faces = 0;
    for (const Face& f : mesh.faces) {
        if (f.right >= 0 && part[f.left] != part[f.right])
            cut_faces++;
    }

    // Вывод
    std::cout << "\n=== RCB Partition Statistics ===\n";
    std::cout << "  Processes    : " << n_parts << "\n";
    std::cout << "  Total cells  : " << N        << "\n";
    std::cout << "  Cut faces    : " << cut_faces << "\n";
    std::cout << "  Cells per rank:\n";

    int ideal = N / n_parts;
    for (int r = 0; r < n_parts; r++) {
        double imbalance = 100.0 * (count[r] - ideal) / std::max(1, ideal);
        std::cout << "    rank " << std::setw(3) << r
                  << " : " << std::setw(7) << count[r] << " cells"
                  << "  (" << std::showpos << std::fixed
                  << std::setprecision(1) << imbalance << "%)\n";
    }
    std::cout << std::noshowpos;
    std::cout << "================================\n\n";
}
