#pragma once
#include "Types.h"
#include <string>


SimParams readParamsTOML(const std::string& path);

SimParams makeCavityParams(int Nx, int Ny,
                            double Re,
                            double U_lid = 1.0,
                            double Lx   = 1.0,
                            double Ly   = 1.0,
                            double T_end = 15.0,
                            double Co_target = 0.4,
                            int nCorr = 2,
                            int step_fo = 100);

Fields initCavity(const SimParams& par);


// Заполнить поля аналитическим начальным условием t=0
Fields initTaylorGreen(const SimParams& par);

// Аналитическое решение в момент времени t (для верификации)
void taylorGreenExact(const SimParams& par, double t,
                      Grid2D& u_ex, Grid2D& v_ex);

Fields initStep(const SimParams& par);
