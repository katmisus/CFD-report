#pragma once
#include "Types.h"


void solvePoissonSOR(Grid2D& pp,
                     const Grid2D& u_w,
                     const Grid2D& v_w,
                     const Grid2D& au_w,
                     const Grid2D& av_w,
                     const SimParams& par,
                     int max_iter = 5000,
                     double tol   = 1e-7,
                     double omega = 1.6);

// Максимальная невязка массы
double massResidual(const Grid2D& u_w,
                    const Grid2D& v_w,
                    const SimParams& par);
