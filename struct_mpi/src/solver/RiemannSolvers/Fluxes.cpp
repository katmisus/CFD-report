#include <vector>
#include <cmath>
#include "TransformValues.h"
#include "RiemannSolver.h"
#include "Types.h"

extern double gamm;

State PhysicalFlux(const State& W, int dir) {
    State F;

    double rho = W[0];
    double u   = W[1];
    double v   = W[2];
    double P   = W[NEQ - 1];

    double E = P/(gamm-1.0)
               + 0.5*rho*(u*u + v*v);

    if (dir == 0) {
        F[0] = rho*u;
        F[1] = rho*u*u + P;
        F[2] = rho*u*v;
        F[NEQ - 1] = u*(E + P);
    }
    if (dir == 1) {
        F[0] = rho*v;
        F[1] = rho*u*v;
        F[2] = rho*v*v + P;
        F[NEQ - 1] = v*(E + P);
    }

    return F;
}

State HLLFlux(const State& WL,
              const State& WR,
              int dir) {
    State F;

    double rhoL = WL[0];
    double rhoR = WR[0];

    double unL = (dir == 0) ? WL[1] : WL[2];
    double unR = (dir == 0) ? WR[1] : WR[2];

    double PL = WL[3];
    double PR = WR[3];

    double aL = std::sqrt(gamm * PL / rhoL);
    double aR = std::sqrt(gamm * PR / rhoR);

    double SL = std::min(unL - aL, unR - aR);
    double SR = std::max(unL + aL, unR + aR);

    State FL = PhysicalFlux(WL, dir);
    State FR = PhysicalFlux(WR, dir);

    State UL, UR;
    ConvertWtoU(WL, UL);
    ConvertWtoU(WR, UR);

    if (SL >= 0.0)
        return FL;

    if (SR <= 0.0)
        return FR;

    for (int k = 0; k < NEQ; k++)
        F[k] = (SR*FL[k] - SL*FR[k]
               + SL*SR*(UR[k] - UL[k]))
               / (SR - SL);

    return F;
}

// CHECK: HLLC_SOLVER
State HLLCFlux(const State& WL,
               const State& WR,
               int dir) {
    State F;

    double rhoL = WL[0];
    double rhoR = WR[0];

    double uL = WL[1];
    double vL = WL[2];
    double PL = WL[3];

    double uR = WR[1];
    double vR = WR[2];
    double PR = WR[3];

    // --- Нормальная и тангенциальная скорости ---
    double unL = (dir == 0) ? uL : vL;
    double unR = (dir == 0) ? uR : vR;

    double ut1L = (dir == 0) ? vL : uL;
    double ut1R = (dir == 0) ? vR : uR;

    // --- Полная энергия ---
    double EL = PL/(gamm-1.0) + 0.5*rhoL*(uL*uL + vL*vL);
    double ER = PR/(gamm-1.0) + 0.5*rhoR*(uR*uR + vR*vR);

    // --- Скорости звука ---
    double aL = std::sqrt(gamm * PL / rhoL);
    double aR = std::sqrt(gamm * PR / rhoR);

    // --- Волновые скорости ---
    double SL = std::min(unL - aL, unR - aR);
    double SR = std::max(unL + aL, unR + aR);

    // --- Контактная скорость ---
    double SM =
        (PR - PL
         + rhoL*unL*(SL - unL)
         - rhoR*unR*(SR - unR))
        /
        (rhoL*(SL - unL)
         - rhoR*(SR - unR));

    // --- Физические потоки ---
    State FL = PhysicalFlux(WL, dir);
    State FR = PhysicalFlux(WR, dir);

    State UL, UR;
    ConvertWtoU(WL, UL);
    ConvertWtoU(WR, UR);

    // --- Левая область ---
    if (SL >= 0.0)
        return FL;

    // --- Правая область ---
    if (SR <= 0.0)
        return FR;

    // --- Левая звёздная область ---
    if (SL <= 0.0 && SM >= 0.0) {
        double rho_star =
            rhoL * (SL - unL) / (SL - SM);

        double un_star = SM;
        double ut1_star = ut1L;

        double P_star =
            PL + rhoL*(SL - unL)*(SM - unL);

        double E_star =
            ( (SL - unL)*EL
              - PL*unL
              + P_star*SM )
            / (SL - SM);

        State U_star;

        U_star[0] = rho_star;

        if (dir == 0) {
            U_star[1] = rho_star * un_star;
            U_star[2] = rho_star * ut1_star;
        } else {
            U_star[1] = rho_star * ut1_star;
            U_star[2] = rho_star * un_star;
        }

        U_star[3] = E_star;

        State F_star;

        for (int k = 0; k < NEQ; k++)
            F_star[k] =
                FL[k] + SL*(U_star[k] - UL[k]);

        return F_star;
    }

    // --- Правая звёздная область ---
    if (SM <= 0.0 && SR >= 0.0) {
        double rho_star =
            rhoR * (SR - unR) / (SR - SM);

        double un_star = SM;
        double ut1_star = ut1R;

        double P_star =
            PR + rhoR*(SR - unR)*(SM - unR);

        double E_star =
            ( (SR - unR)*ER
              - PR*unR
              + P_star*SM )
            / (SR - SM);

        State U_star;

        U_star[0] = rho_star;

        if (dir == 0) {
            U_star[1] = rho_star * un_star;
            U_star[2] = rho_star * ut1_star;
        } else {
            U_star[1] = rho_star * ut1_star;
            U_star[2] = rho_star * un_star;
        }

        U_star[3] = E_star;

        State F_star;

        for (int k = 0; k < NEQ; k++)
            F_star[k] =
                FR[k] + SR*(U_star[k] - UR[k]);

        return F_star;
    }

    return F; // сюда обычно не попадает
}

State RusanovFlux(const State& WL,
                  const State& WR,
                  int dir)
{
    const double gamma = gamm;

    // ===== Primitive =====
    double rhoL = WL[0];
    double uL   = WL[1];
    double vL   = WL[2];
    double pL   = WL[3];

    double rhoR = WR[0];
    double uR   = WR[1];
    double vR   = WR[2];
    double pR   = WR[3];

    // ===== Normal velocity =====
    double unL = (dir==0)? uL : vL;
    double unR = (dir==0)? uR : vR;

    // ===== Sound speed =====
    double aL = std::sqrt(gamma*pL/rhoL);
    double aR = std::sqrt(gamma*pR/rhoR);

    // ===== Max wave speed (Toro) =====
    double Splus =
        std::max(std::abs(unL)+aL,
                 std::abs(unR)+aR);

    // ===== Conservative states =====
    State UL, UR;
    ConvertWtoU(WL, UL);
    ConvertWtoU(WR, UR);

    // ===== Physical fluxes =====
    State FL = PhysicalFlux(WL, dir);
    State FR = PhysicalFlux(WR, dir);

    // ===== Rusanov flux =====
    State F;
    for(int k=0; k<NEQ; k++)
        F[k] = 0.5*(FL[k] + FR[k])
               - 0.5*Splus*(UR[k] - UL[k]);

    return F;
}


// double phi_lambda(double lambda, double delta = 1e-6)
// {
//     if (std::abs(lambda) >= delta)
//         return std::abs(lambda);
//     else
//         return (lambda*lambda + delta*delta)/(2.0*delta);
// }

State RoeFlux(const State& WL,
              const State& WR,
              int dir) {
    State F;

    // --- Primitive ---
    double rhoL = WL[0];
    double uL   = WL[1];
    double vL   = WL[2];
    double PL   = WL[3];

    double rhoR = WR[0];
    double uR   = WR[1];
    double vR   = WR[2];
    double PR   = WR[3];

    // --- Normal & tangential velocities ---
    double unL = (dir == 0) ? uL : vL;
    double unR = (dir == 0) ? uR : vR;

    double utL = (dir == 0) ? vL : uL;
    double utR = (dir == 0) ? vR : uR;

    // --- Total energy ---
    double EL = PL/(gamm-1.0)
              + 0.5*rhoL*(uL*uL + vL*vL);

    double ER = PR/(gamm-1.0)
              + 0.5*rhoR*(uR*uR + vR*vR);

    double HL = (EL + PL)/rhoL;
    double HR = (ER + PR)/rhoR;

    // --- Roe averages ---
    double sqL = std::sqrt(rhoL);
    double sqR = std::sqrt(rhoR);
    double denom = sqL + sqR;

    double un_tilde = (sqL*unL + sqR*unR)/denom;
    double ut_tilde = (sqL*utL + sqR*utR)/denom;
    double H_tilde  = (sqL*HL  + sqR*HR )/denom;

    double a_tilde =
        std::sqrt((gamm-1.0) *
        (H_tilde - 0.5*(un_tilde*un_tilde
                      + ut_tilde*ut_tilde)));

    double rho_tilde = sqL*sqR;

    // --- Eigenvalues ---
    double lambda1 = un_tilde - a_tilde;
    double lambda2 = un_tilde;
    double lambda3 = un_tilde + a_tilde;

    auto EntrophyFix = [&](double& lambda) {
        double delta = 1e-6;
            if (std::abs(lambda) >= delta)
            return std::abs(lambda);
        else
            return (lambda*lambda + delta*delta)/(2.0*delta);
    };

    lambda1 = EntrophyFix(lambda1);
    lambda2 = EntrophyFix(lambda2);
    lambda3 = EntrophyFix(lambda3);

    // --- Conservative vars ---
    State UL, UR;
    ConvertWtoU(WL, UL);
    ConvertWtoU(WR, UR);

    State dU;
    for (int k = 0; k < NEQ; k++)
        dU[k] = UR[k] - UL[k];

    double dP   = PR - PL;
    double drho = rhoR - rhoL;
    double dun  = unR - unL;

    // --- Wave strengths ---
    double alpha2 =
        drho - dP/(a_tilde*a_tilde);

    double alpha1 =
        0.5*(dP - rho_tilde*a_tilde*dun)
        /(a_tilde*a_tilde);

    double alpha3 =
        0.5*(dP + rho_tilde*a_tilde*dun)
        /(a_tilde*a_tilde);

    // --- Eigenvectors ---
    State r1, r2, r3;

    if (dir == 0)
    {
        r1 = {1.0,
              un_tilde - a_tilde,
              ut_tilde,
              H_tilde - un_tilde*a_tilde};

        r2 = {1.0,
              un_tilde,
              ut_tilde,
              0.5*(un_tilde*un_tilde
                 + ut_tilde*ut_tilde)};

        r3 = {1.0,
              un_tilde + a_tilde,
              ut_tilde,
              H_tilde + un_tilde*a_tilde};
    }
    else
    {
        r1 = {1.0,
              ut_tilde,
              un_tilde - a_tilde,
              H_tilde - un_tilde*a_tilde};

        r2 = {1.0,
              ut_tilde,
              un_tilde,
              0.5*(un_tilde*un_tilde
                 + ut_tilde*ut_tilde)};

        r3 = {1.0,
              ut_tilde,
              un_tilde + a_tilde,
              H_tilde + un_tilde*a_tilde};
    }

    // --- Physical flux ---
    State FL = PhysicalFlux(WL, dir);
    State FR = PhysicalFlux(WR, dir);

    for (int k = 0; k < NEQ; k++)
    {
        double diss =
            lambda1*alpha1*r1[k]
          + lambda2*alpha2*r2[k]
          + lambda3*alpha3*r3[k];

        F[k] = 0.5*(FL[k] + FR[k])
             - 0.5*diss;
    }

    return F;
}

State OsherFlux(const State& WL,
                const State& WR,
                int dir) {
    State F_osher;

    // --- Primitive ---
    double rhoL = WL[0];
    double uL   = WL[1];
    double vL   = WL[2];
    double PL   = WL[3];

    double rhoR = WR[0];
    double uR   = WR[1];
    double vR   = WR[2];
    double PR   = WR[3];

    // --- Normal & tangential ---
    double unL = (dir==0)? uL : vL;
    double unR = (dir==0)? uR : vR;

    //double utL = (dir==0)? vL : uL;
    //double utR = (dir==0)? vR : uR;

    double cL = std::sqrt(gamm * PL / rhoL);
    double cR = std::sqrt(gamm * PR / rhoR);

    // --- Star pressure ---
    double z = (gamm - 1.0)/(2.0*gamm);

    double P_star = std::pow(
                    (cL + cR
                        - 0.5*(unR-unL)*(gamm-1.0))
                        /
                    (cL/std::pow(PL,z)
                        + cR/std::pow(PR,z)),
                    1.0/z);

    // --- Star velocity ---
    double H = std::pow(PL/PR, z);

    double un_star =
        (H*unL/cL + unR/cR
         + 2.0*(H-1.0)*(gamm-1.0))
        /
        (H/cL + 1.0/cR);

    // --- Star densities ---
    double rho_star_L =
        rhoL*std::pow(P_star/PL, 1.0/gamm);

    double rho_star_R =
        rhoR*std::pow(P_star/PR, 1.0/gamm);

    double c_star_L =
        std::sqrt(gamm*P_star/rho_star_L);

    double c_star_R =
        std::sqrt(gamm*P_star/rho_star_R);

    // --- Build star states ---
    State W_star_L = WL;
    State W_star_R = WR;

    W_star_L[0] = rho_star_L;
    W_star_L[NEQ - 1] = P_star;
    W_star_R[0] = rho_star_R;
    W_star_R[NEQ - 1] = P_star;

    if (dir==0) {
        W_star_L[1] = un_star;
        W_star_R[1] = un_star;
    }
    else {
        W_star_L[2] = un_star;
        W_star_R[2] = un_star;
    }

    // --- Fluxes ---
    State F_L = PhysicalFlux(WL, dir);
    State F_R = PhysicalFlux(WR, dir);

    State F_star_L = PhysicalFlux(W_star_L, dir);
    State F_star_R = PhysicalFlux(W_star_R, dir);

    // --- Sonic states (как у тебя) ---
    double un_S0 =
        (gamm-1)/(gamm+1)*unL
        + 2*cL/(gamm+1);

    double rho_S0 =
        rhoL*std::pow(un_S0/cL,
        2.0/(gamm-1));

    double P_S0 =
        PL*std::pow(rho_S0/rhoL,
        gamm);

    double un_S1 =
        (gamm-1)/(gamm+1)*unR
        - 2*cR/(gamm+1);

    double rho_S1 =
        rhoR*std::pow((-un_S1)/cR,
        2.0/(gamm-1));

    double P_S1 =
        PR*std::pow(rho_S1/rhoR,
        gamm);

    State W_S0 = WL;
    State W_S1 = WR;

    W_S0[0] = rho_S0;
    W_S0[NEQ - 1] = P_S0;
    W_S1[0] = rho_S1;
    W_S1[NEQ - 1] = P_S1;

    if (dir==0) {
        W_S0[1] = un_S0;
        W_S1[1] = un_S1;
    }
    else {
        W_S0[2] = un_S0;
        W_S1[2] = un_S1;
    }

    State F_S0 = PhysicalFlux(W_S0, dir);
    State F_S1 = PhysicalFlux(W_S1, dir);

    // --- Integrals ---
    State I1 = {0,0,0,0};
    State I2 = {0,0,0,0};
    State I3 = {0,0,0,0};

    double lam_L  = unL - cL;
    double lam_13 = un_star - c_star_L;

    if (lam_L <= 0.0 && lam_13 <= 0.0)
        for(int k=0;k<NEQ;k++)
            I1[k] = F_star_L[k] - F_L[k];

    else if (lam_L >= 0.0 && lam_13 <= 0.0)
        for(int k=0;k<NEQ;k++)
            I1[k] = F_star_L[k] - F_S0[k];

    else if (lam_L <= 0.0 && lam_13 >= 0.0)
        for(int k=0;k<NEQ;k++)
            I1[k] = F_S0[k] - F_L[k];

    if (un_star < 0.0)
        for(int k=0;k<NEQ;k++)
            I2[k] = F_star_R[k] - F_star_L[k];

    double lam_23 = un_star + c_star_R;
    double lam_R = unR + cR;

    if (lam_23 <= 0.0 && lam_R <= 0.0)
        for(int k=0; k<NEQ; k++)
            I3[k] = F_R[k] - F_star_R[k];
    else if (lam_23 >= 0.0 && lam_R <= 0.0)
        for(int k=0; k<NEQ; k++)
            I3[k] = F_R[k] - F_S1[k];
    else if (lam_23 <= 0.0 && lam_R >= 0.0)
        for(int k=0; k<NEQ; k++)
            I3[k] = F_S1[k] - F_star_R[k];

    for(int k=0; k<NEQ; k++)
        F_osher[k] = F_L[k] + I1[k] + I2[k] + I3[k];

    return F_osher;
}

State ExactFlux(const State& WL,
                const State& WR,
                int dir) {
    State W_star;

    NewtonForPressure(WL, WR, W_star, 1e-6, dir);

    W_star = GetParamsFromChoosingWave(WL, WR, W_star,
            						   0.0, 1.0, dir);

    return PhysicalFlux(W_star, dir);
}