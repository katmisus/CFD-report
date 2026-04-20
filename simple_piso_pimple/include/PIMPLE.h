#pragma once
#include "Types.h"
#include <string>

// Выполнить один временной шаг PIMPLE
void stepPIMPLE(Fields& f, const SimParams& par);

// Запустить полный расчёт PIMPLE (временной цикл)
void runPIMPLE(Fields& f, const SimParams& par,
               const std::string& output_dir);
