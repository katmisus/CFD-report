#pragma once
#include "Types.h"

// Выполнить один временной шаг PISO
void stepPISO(Fields& f, const SimParams& par);

// Запустить полный расчёт PISO (временной цикл)
void runPISO(Fields& f, const SimParams& par,
             const std::string& output_dir);
