#pragma once
#include "Types.h"

void applyBoundaryConditions(Fields& f, const SimParams& par);

void applyBCCavity(Fields& f, const SimParams& par);

// Периодические ГУ для вихрей Тейлора–Грина
void applyBCTaylorGreen(Fields& f, const SimParams& par);

// ГУ для обратной ступеньки:
void applyBCStep(Fields& f, const SimParams& par);
