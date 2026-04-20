#pragma once
#include "Types.h"

void NonStatMomentumPredictor(Fields& f,
                       const Fields& f_old,
                       const SimParams& par,
                       double alpha_u = 1.0);

void StatMomentumPredictor(Fields& f, const SimParams& par);