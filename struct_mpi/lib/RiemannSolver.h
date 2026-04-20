#ifndef _RIEMANN_SOLVER_H_
#define _RIEMANN_SOLVER_H_
#include "Types.h"

void GetComponents(State W, double& rho, double& u_n, double& u_t_1, double& P, int dir);

void PressureInitialGuess(State W_L, State W_R, double& P_prev, int dir);

void FindValuesOfFunctions(double P, double P_k, double rho_k, double a_k, double& func_k, double& der_func_k);

void NewtonForPressure(State W_L, State W_R, State& W_star, double eps, int dir);

State GetParamsFromChoosingWave(State W_L, State W_R, State W_star, double x, double t, int dir);

#endif
