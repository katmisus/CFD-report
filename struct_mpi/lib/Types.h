#ifndef _TYPES_H_
#define _TYPES_H_

#include <vector>
#include <array>
#include <filesystem>
#include <functional>

// CHECK: STATE_VECTOR_2D
constexpr int NEQ = 4; // Количество переменных в векторе состояния (1D - 3; 2D - 4; 3D - 5)
using State = std::array<double, NEQ>; // Состояние (включает в себя NEQ переменных)
using Field = std::vector<std::vector<State>>; // Двумерное поле состояний (X x Y -> State)

namespace fs = std::filesystem;

using LimiterFunction = std::function<double(double)>;
using RecLimiterFunction = std::function<double(double, double)>;

#endif