#ifndef _MADER_H_
#define _MADER_H

#include "Types.h"

void EOS(Field rho, Field I, Field& P, Field& T);
void Arrenius(Field& mass_fraction, Field T, double dt, int Nx_tot, int Ny_tot);
void Viscosity(Field& q1, Field& q2, Field& q3, Field& q4, Field rho, Field u, Field v);
void VelocityTilde(Field& u_tilde, Field& v_tilde, Field P, Field rho, const std::vector<double>& x, const std::vector<double>& y, Field q1, Field q2, Field q3, Field q4, Field u, Field v, double dt);
void ZIPEnergy(Field I, Field& I_tilde, Field P, double dt, Field rho, Field u, Field u_tilde, Field v, Field v_tilde, Field q1, Field q2, Field q3, Field q4, const std::vector<double>& x, const std::vector<double>& y);
void ChangingFluxes(Field& alpha, Field& beta, Field& DM, Field& DE, Field& DW, Field& DPU, Field& DPV, Field u_tilde, Field v_tilde, const std::vector<double>& x, const std::vector<double>& y, double dt, Field rho, Field I_tilde, Field u, Field v, Field mass_fraction);
void Repartition(Field& rho, Field& u, Field& v, Field& I, Field& mass_fraction, Field I_tilde, Field u_tilde, Field v_tilde, Field DM, Field DE, Field DW, Field DPU, Field DPV);

void Mader(Field& W_new, const Field& W, const std::vector<double>& x, const std::vector<double>& y, double dt);

#endif