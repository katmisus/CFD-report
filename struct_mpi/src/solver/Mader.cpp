#include <cmath>
#include <vector>
#include <iostream>
#include "Mader.h"
#include "FileProcessing.h"
#include "BoundCond.h"
#include "Init.h"
#include "Types.h"

extern double gamm, Lx, Ly, T_init, R_gas, M, P_min, E_act, Z_freq,
              VISC, MINWT, GASW, MINGRHO, Q_chem;
extern int Nx, Ny, fict;
extern Field mass_fraction;
extern std::vector<std::vector<bool>> reacted;

static constexpr double WALL_RHO = 100.0;

static inline bool IsWall(double rho) { return rho > WALL_RHO; }

// CHECK: MADER_EOS
static void EOS(const Field& W,
                std::vector<std::vector<double>>& P,
                std::vector<std::vector<double>>& T,
                std::vector<std::vector<double>>& I,
                int Nx_tot, int Ny_tot)
{
    for (int i = 0; i < Nx_tot; i++) {
        for (int j = 0; j < Ny_tot; j++) {
            double rho = std::max(W[i][j][0], 1e-14);
            double p   = W[i][j][NEQ - 1];
            I[i][j]    = p / ((gamm - 1.0) * rho);
            P[i][j]    = std::max((gamm - 1.0) * rho * I[i][j], P_min);
            T[i][j]    = P[i][j] / (rho * R_gas * 1e-5 / M);
        }
    }
}

// CHECK: MADER_ARRHENIUS
static void ChemicalKinetics(const std::vector<std::vector<double>>& P,
                             const std::vector<std::vector<double>>& rho,
                             const std::vector<std::vector<double>>& T,
                             std::vector<std::vector<double>>& I,
                             int Nx_tot, int Ny_tot)
{
    const double P_threshold = 0.05;

    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            if (IsWall(rho[i][j])) continue;
            if (reacted[i][j]) continue;
            double mf = mass_fraction[i][j][0];
            // double t = T[i][j];
            if (mf <= GASW) { reacted[i][j] = true; continue; }
            //if (t <= MINWT) { reacted[i][j] = true; continue; }

            if (P[i][j] > P_threshold) {
                I[i][j]               += Q_chem * mf;
                mass_fraction[i][j][0] = 0.0;
                reacted[i][j]          = true;
            }
        }
    }
}

// CHECK: MADER_VISC
static void ComputeViscosity(const std::vector<std::vector<double>>& U,
                              const std::vector<std::vector<double>>& V,
                              const std::vector<std::vector<double>>& rho,
                              std::vector<std::vector<double>>& Q1,
                              std::vector<std::vector<double>>& Q2,
                              std::vector<std::vector<double>>& Q3,
                              std::vector<std::vector<double>>& Q4,
                              int Nx_tot, int Ny_tot)
{
    for (int i = fict; i < Nx + fict - 1; i++)
        for (int j = fict; j < Ny + fict - 1; j++) {
            if (IsWall(rho[i][j])) { Q1[i][j]=Q2[i][j]=0.0; continue; }

            if (!IsWall(rho[i][j+1]) && V[i][j] >= V[i][j+1])
                Q1[i][j] = VISC * rho[i][j] * (V[i][j] - V[i][j+1]);
            else
                Q1[i][j] = 0.0;

            if (!IsWall(rho[i+1][j]) && U[i][j] >= U[i+1][j])
                Q2[i][j] = VISC * rho[i][j] * (U[i][j] - U[i+1][j]);
            else
                Q2[i][j] = 0.0;
        }
    for (int i = fict; i < Nx + fict - 1; i++)
        for (int j = fict; j < Ny + fict - 1; j++) {
            Q3[i][j] = (j > fict) ? Q1[i][j-1] : 0.0;
            Q4[i][j] = (i > fict) ? Q2[i-1][j] : 0.0;
        }
}

// CHECK: MADER_VELOCITY
static void VelocityTilde(const std::vector<std::vector<double>>& U,
                           const std::vector<std::vector<double>>& V,
                           const std::vector<std::vector<double>>& P,
                           const std::vector<std::vector<double>>& rho,
                           const std::vector<std::vector<double>>& Q1,
                           const std::vector<std::vector<double>>& Q2,
                           const std::vector<std::vector<double>>& Q3,
                           const std::vector<std::vector<double>>& Q4,
                           std::vector<std::vector<double>>& U_tilde,
                           std::vector<std::vector<double>>& V_tilde,
                           const std::vector<double>& x,
                           const std::vector<double>& y,
                           double dt, int Nx_tot, int Ny_tot)
{
    for (int i = fict; i < Nx + fict - 1; i++)
        for (int j = fict; j < Ny + fict - 1; j++) {
            if (IsWall(rho[i][j])) {
                U_tilde[i][j] = 0.0;
                V_tilde[i][j] = 0.0;
                continue;
            }
            double dR = x[i] - x[i-1];
            double dZ = y[j] - y[j-1];
            double r  = std::max(rho[i][j], 1e-14);
            V_tilde[i][j] = V[i][j] - dt / (r * 2.0 * dZ) * ((P[i][j+1]-P[i][j-1]) + (Q3[i][j]-Q1[i][j]));
            U_tilde[i][j] = U[i][j] - dt / (r * 2.0 * dR) * ((P[i+1][j]-P[i-1][j]) + (Q4[i][j]-Q2[i][j]));
        }
}

// CHECK: MADER_ZIP
static void ZIPEnergy(const std::vector<std::vector<double>>& I,
                       std::vector<std::vector<double>>& I_tilde,
                       const std::vector<std::vector<double>>& P,
                       const std::vector<std::vector<double>>& rho,
                       const std::vector<std::vector<double>>& U,
                       const std::vector<std::vector<double>>& V,
                       const std::vector<std::vector<double>>& U_tilde,
                       const std::vector<std::vector<double>>& V_tilde,
                       const std::vector<std::vector<double>>& Q1,
                       const std::vector<std::vector<double>>& Q2,
                       const std::vector<std::vector<double>>& Q3,
                       const std::vector<std::vector<double>>& Q4,
                       const std::vector<double>& x,
                       const std::vector<double>& y,
                       double dt, int Nx_tot, int Ny_tot)
{
    for (int i = fict; i < Nx + fict - 1; i++)
        for (int j = fict; j < Ny + fict - 1; j++) {
            if (IsWall(rho[i][j])) { I_tilde[i][j] = I[i][j]; continue; }
            double dR = x[i] - x[i-1];
            double dZ = y[j] - y[j-1];
            double r  = std::max(rho[i][j], 1e-14);

            double U1 = U[i-1][j] + U_tilde[i-1][j];
            double U2 = U[i+1][j] + U_tilde[i+1][j];
            double V1 = V[i][j-1] + V_tilde[i][j-1];
            double V2 = V[i][j+1] + V_tilde[i][j+1];
            double T3 = U[i][j]   + U_tilde[i][j];
            double T1 = V[i][j]   + V_tilde[i][j];

            double zip = (P[i][j]/dR)*(U2-U1)
                       + (Q4[i][j]/dR)*(U2-T3)
                       + (Q2[i][j]/dR)*(T3-U1)
                       + (P[i][j]/dZ)*(V2-V1)
                       + (Q3[i][j]/dZ)*(V2-T1)
                       + (Q1[i][j]/dZ)*(T1-V1);

            I_tilde[i][j] = std::max(I[i][j] - dt/(4.0*r)*zip, 0.0);
        }
}




static void TransportAndRepartition(
    std::vector<std::vector<double>>& rho,
    std::vector<std::vector<double>>& U,
    std::vector<std::vector<double>>& V,
    std::vector<std::vector<double>>& I,
    const std::vector<std::vector<double>>& U_tilde,
    const std::vector<std::vector<double>>& V_tilde,
    const std::vector<std::vector<double>>& I_tilde,
    const std::vector<double>& x,
    const std::vector<double>& y,
    double dt, int Nx_tot, int Ny_tot)
{
    std::vector<std::vector<double>> DM (Nx_tot, std::vector<double>(Ny_tot, 0.0));
    std::vector<std::vector<double>> DE (Nx_tot, std::vector<double>(Ny_tot, 0.0));
    std::vector<std::vector<double>> DW (Nx_tot, std::vector<double>(Ny_tot, 0.0));
    std::vector<std::vector<double>> DPU(Nx_tot, std::vector<double>(Ny_tot, 0.0));
    std::vector<std::vector<double>> DPV(Nx_tot, std::vector<double>(Ny_tot, 0.0));

    // CHECK: MADER_DONOR
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {

            if (IsWall(rho[i][j])) continue;

            double dR = x[i] - x[i-1];
            double dZ = y[j] - y[j-1];

            {
                bool skip_r = (i > fict && IsWall(rho[i-1][j]));
                if (!skip_r) 
                {
                    double denom = 1.0 + (U_tilde[i-1][j] - U_tilde[i][j]) * dt / dR;
                    if (std::abs(denom) < 1e-14) denom = (denom >= 0) ? 1e-14 : -1e-14;
                    double alpha = 0.5*(U_tilde[i-1][j]+U_tilde[i][j])*dt/dR / denom;

                    int di = (alpha >= 0.0) ? i-1 : i;
                    int ai = (alpha >= 0.0) ? i   : i-1;

                    if (!IsWall(rho[di][j])) 
                    {
                        double rho_d = std::max(rho[di][j], 1e-14);
                        double DMASS = rho_d * std::abs(alpha);
                        double E_d   = I_tilde[di][j] + 0.5*(U_tilde[di][j]*U_tilde[di][j] + V_tilde[di][j]*V_tilde[di][j]);
                        double mf_d  = mass_fraction[di][j][0];

                        // CHECK: SHARGATOV 
                        double mf_transfer = (reacted[ai][j] || reacted[di][j]) ? 0.0 : mf_d;

                        DM [ai][j] += DMASS;       
                        DM [di][j] -= DMASS;
                        DE [ai][j] += E_d*DMASS;   
                        DE [di][j] -= E_d*DMASS;
                        DPU[ai][j] += U_tilde[di][j]*DMASS;
                        DPU[di][j] -= U_tilde[di][j]*DMASS;
                        DPV[ai][j] += V_tilde[di][j]*DMASS;
                        DPV[di][j] -= V_tilde[di][j]*DMASS;
                        DW [ai][j] += mf_transfer*DMASS;
                        DW [di][j] -= mf_transfer*DMASS;
                    }
                }
            }

            {
                bool skip_z = (j > fict && IsWall(rho[i][j-1]));
                if (!skip_z)
                {
                    double denom = 1.0 + (V_tilde[i][j-1] - V_tilde[i][j]) * dt / dZ;
                    if (std::abs(denom) < 1e-14) denom = (denom >= 0) ? 1e-14 : -1e-14;
                    double beta = 0.5*(V_tilde[i][j-1]+V_tilde[i][j])*dt/dZ / denom;

                    int dj = (beta >= 0.0) ? j-1 : j;
                    int aj = (beta >= 0.0) ? j   : j-1;

                    if (!IsWall(rho[i][dj]))
                    {
                        double rho_d = std::max(rho[i][dj], 1e-14);
                        double DMASS = rho_d * std::abs(beta);
                        double E_d   = I_tilde[i][dj]
                                    + 0.5*(U_tilde[i][dj]*U_tilde[i][dj]
                                    + V_tilde[i][dj]*V_tilde[i][dj]);
                        double mf_d  = mass_fraction[i][dj][0];

                        // CHECK: SHARGATOV 
                        double mf_transfer = (reacted[i][aj] || reacted[i][dj]) ? 0.0 : mf_d;

                        DM [i][aj] += DMASS;
                        DM [i][dj] -= DMASS;
                        DE [i][aj] += E_d*DMASS;   
                        DE [i][dj] -= E_d*DMASS;
                        DPU[i][aj] += U_tilde[i][dj]*DMASS;
                        DPU[i][dj] -= U_tilde[i][dj]*DMASS;
                        DPV[i][aj] += V_tilde[i][dj]*DMASS;
                        DPV[i][dj] -= V_tilde[i][dj]*DMASS;
                        DW [i][aj] += mf_transfer*DMASS;
                        DW [i][dj] -= mf_transfer*DMASS;
                    }
                }
            }
        }
    }

    // CHECK: MADER_REPARTITION

    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {

            if (IsWall(rho[i][j])) continue;

            double rho_new = rho[i][j] + DM[i][j];
            if (rho_new <= 1e-12) {
                rho[i][j] = MINGRHO; 
                U[i][j] = 0.0; 
                V[i][j] = 0.0;
                I[i][j]   = 0.0; 
                mass_fraction[i][j][0] = 0.0;
                continue;
            }

            double U_new = (rho[i][j]*U_tilde[i][j] + DPU[i][j]) / rho_new;
            double V_new = (rho[i][j]*V_tilde[i][j] + DPV[i][j]) / rho_new;
            double E_old = I_tilde[i][j] + 0.5*(U_tilde[i][j]*U_tilde[i][j]+V_tilde[i][j]*V_tilde[i][j]);
            double I_new = std::max((rho[i][j]*E_old + DE[i][j])/rho_new - 0.5*(U_new*U_new + V_new*V_new), 0.0);
            double mf_new = (rho[i][j]*mass_fraction[i][j][0] + DW[i][j]) / rho_new;

            rho[i][j] = rho_new;
            U[i][j]   = U_new;
            V[i][j]   = V_new;
            I[i][j]   = I_new;
            mass_fraction[i][j][0] = reacted[i][j] ? 0.0 : std::max(mf_new, 0.0);
        }
    }
}


static void FillGhost(std::vector<std::vector<double>>& F, int Nx_tot, int Ny_tot)
{
    int il = fict, ir = Nx + fict - 2;
    int jb = fict, jt = Ny + fict - 2;
    for (int j = 0; j < Ny_tot; j++)
        for (int g = 0; g < fict; g++) {
            F[g][j]          = F[il][j];
            F[Nx_tot-1-g][j] = F[ir][j];
        }
    for (int i = 0; i < Nx_tot; i++)
        for (int g = 0; g < fict; g++) {
            F[i][g]          = F[i][jb];
            F[i][Ny_tot-1-g] = F[i][jt];
        }
}

static void FillGhostMF(int Nx_tot, int Ny_tot)
{
    int il = fict, ir = Nx + fict - 2;
    int jb = fict, jt = Ny + fict - 2;
    for (int j = 0; j < Ny_tot; j++)
        for (int g = 0; g < fict; g++) {
            mass_fraction[g][j]          = mass_fraction[il][j];
            mass_fraction[Nx_tot-1-g][j] = mass_fraction[ir][j];
        }
    for (int i = 0; i < Nx_tot; i++)
        for (int g = 0; g < fict; g++) {
            mass_fraction[i][g]          = mass_fraction[i][jb];
            mass_fraction[i][Ny_tot-1-g] = mass_fraction[i][jt];
        }
}


void Mader(Field& W_new, const Field& W,
           const std::vector<double>& x,
           const std::vector<double>& y,
           double dt)
{
    int Nx_tot = Nx + 2*fict - 1;
    int Ny_tot = Ny + 2*fict - 1;

    std::vector<std::vector<double>> rho(Nx_tot, std::vector<double>(Ny_tot));
    std::vector<std::vector<double>> U  (Nx_tot, std::vector<double>(Ny_tot));
    std::vector<std::vector<double>> V  (Nx_tot, std::vector<double>(Ny_tot));
    std::vector<std::vector<double>> P  (Nx_tot, std::vector<double>(Ny_tot));
    std::vector<std::vector<double>> T  (Nx_tot, std::vector<double>(Ny_tot));
    std::vector<std::vector<double>> I  (Nx_tot, std::vector<double>(Ny_tot));

    for (int i = 0; i < Nx_tot; i++)
        for (int j = 0; j < Ny_tot; j++) {
            rho[i][j] = W[i][j][0];
            U[i][j]   = W[i][j][1];
            V[i][j]   = W[i][j][2];
            P[i][j]   = W[i][j][NEQ-1];
        }

    // ЕОS
    EOS(W, P, T, I, Nx_tot, Ny_tot);
    FillGhost(P, Nx_tot, Ny_tot);
    FillGhost(T, Nx_tot, Ny_tot);
    FillGhost(I, Nx_tot, Ny_tot);

    // Реакция
    ChemicalKinetics(P, rho, T, I, Nx_tot, Ny_tot);
    FillGhost(I, Nx_tot, Ny_tot);
    FillGhostMF(Nx_tot, Ny_tot);

    // Вязкость
    std::vector<std::vector<double>> Q1(Nx_tot, std::vector<double>(Ny_tot, 0.0));
    std::vector<std::vector<double>> Q2(Nx_tot, std::vector<double>(Ny_tot, 0.0));
    std::vector<std::vector<double>> Q3(Nx_tot, std::vector<double>(Ny_tot, 0.0));
    std::vector<std::vector<double>> Q4(Nx_tot, std::vector<double>(Ny_tot, 0.0));

    ComputeViscosity(U, V, rho, Q1, Q2, Q3, Q4, Nx_tot, Ny_tot);
    FillGhost(Q1, Nx_tot, Ny_tot);
    FillGhost(Q2, Nx_tot, Ny_tot);
    FillGhost(Q3, Nx_tot, Ny_tot);
    FillGhost(Q4, Nx_tot, Ny_tot);

    // Скорости
    std::vector<std::vector<double>> U_tilde(Nx_tot, std::vector<double>(Ny_tot));
    std::vector<std::vector<double>> V_tilde(Nx_tot, std::vector<double>(Ny_tot));
    for (int i = 0; i < Nx_tot; i++)
        for (int j = 0; j < Ny_tot; j++) {
            U_tilde[i][j] = U[i][j];
            V_tilde[i][j] = V[i][j];
        }

    VelocityTilde(U, V, P, rho, Q1, Q2, Q3, Q4,
                  U_tilde, V_tilde, x, y, dt, Nx_tot, Ny_tot);
    FillGhost(U_tilde, Nx_tot, Ny_tot);
    FillGhost(V_tilde, Nx_tot, Ny_tot);

    // ZIP Energy
    std::vector<std::vector<double>> I_tilde(Nx_tot, std::vector<double>(Ny_tot));
    for (int i = 0; i < Nx_tot; i++)
        for (int j = 0; j < Ny_tot; j++)
            I_tilde[i][j] = I[i][j];

    ZIPEnergy(I, I_tilde, P, rho, U, V, U_tilde, V_tilde,
              Q1, Q2, Q3, Q4, x, y, dt, Nx_tot, Ny_tot);
    FillGhost(I_tilde, Nx_tot, Ny_tot);

    // Перераспределение
    TransportAndRepartition(rho, U, V, I,
                             U_tilde, V_tilde, I_tilde,
                             x, y, dt, Nx_tot, Ny_tot);

    // Финал
    for (int i = fict; i < Nx + fict - 1; i++)
        for (int j = fict; j < Ny + fict - 1; j++) {
            if (IsWall(rho[i][j])) {
                W_new[i][j] = W[i][j];
                continue;
            }
            W_new[i][j][0] = rho[i][j];
            W_new[i][j][1] = U[i][j];
            W_new[i][j][2] = V[i][j];

            double kinetic = 0.5*rho[i][j]*(U[i][j]*U[i][j]+V[i][j]*V[i][j]);
            double E       = rho[i][j]*I[i][j] + kinetic;
            W_new[i][j][NEQ-1] = std::max(P_min, (gamm-1.0)*(E-kinetic));
        }
}