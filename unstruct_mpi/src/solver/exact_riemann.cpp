#include "mesh/mesh.hpp"
#include <cmath>
#include <algorithm>

Vec4 flux_n(const Vec4& U, double nx, double ny, double g) {
    double rho = U[0];
    double u   = U[1] / rho;
    double v   = U[2] / rho;
    double E   = U[3];
    double p   = (g - 1.0) * (E - 0.5*rho*(u*u + v*v));
    double un  = u*nx + v*ny;
    return { rho*un,
             rho*u*un + p*nx,
             rho*v*un + p*ny,
             (E + p)*un };
}

double f_wave(double p, double pK, double rhoK, double aK, double g,
                             double& df) {
    // Ударная волна
    if (p > pK) {
        double A = 2.0 / ((g + 1.0) * rhoK);
        double B = (g - 1.0) / (g + 1.0) * pK;
        double s = std::sqrt(A / (p + B));
        double f = (p - pK) * s;
        df = s * (1.0 - (p - pK) / (2.0*(p + B)));
        return f;
    // Волна разрежения
    } else {
        double exp = (g - 1.0) / (2.0 * g);
        double ratio = std::pow(p / pK, exp);
        double f = 2.0 * aK / (g - 1.0) * (ratio - 1.0);
        df = (1.0 / (rhoK * aK)) * std::pow(p / pK, -(g + 1.0) / (2.0 * g));
        return f;
    }
}

double rho_star(double ps, double pK, double rhoK, double g) {
    if (ps > pK) {
        double mu = (g - 1.0) / (g + 1.0);
        return rhoK * (ps/pK + mu) / (mu * ps/pK + 1.0);
    } else {
        return rhoK * std::pow(ps / pK, 1.0 / g);
    }
}

Vec4 exact_riemann(const Vec4& UL, const Vec4& UR,
                   double nx, double ny, double g) {

    double rL  = UL[0];         
    double uL  = UL[1] / rL;    
    double vL  = UL[2] / rL;    
    double EL  = UL[3];         
    double pL  = (g - 1.0) * (EL - 0.5*rL*(uL*uL + vL*vL));
    double aL  = std::sqrt(g * std::max(pL, 1e-12) / rL);
    double unL = uL*nx + vL*ny;   // нормальная скорость
    double utL_x = uL - unL*nx;   // тангенциальная скорость
    double utL_y = vL - unL*ny;

    double rR  = UR[0];
    double uR  = UR[1] / rR;
    double vR  = UR[2] / rR;
    double ER  = UR[3];
    double pR  = (g - 1.0) * (ER - 0.5*rR*(uR*uR + vR*vR));
    double aR  = std::sqrt(g * std::max(pR, 1e-12) / rR);
    double unR = uR*nx + vR*ny;
    double utR_x = uR - unR*nx;
    double utR_y = vR - unR*ny;

   
    // Метод Ньютона для давления
    // приближение двух волн разрежения
    double exp1 = (g - 1.0) / (2.0 * g);
    double pTR  = std::pow(
        (aL + aR - 0.5*(g - 1.0)*(unR - unL)) /
        (aL / std::pow(pL, exp1) + aR / std::pow(pR, exp1)),
        1.0 / exp1
    );
    double ps = std::max(pTR, 1e-12);

    const int MAX_ITER = 50;
    const double TOL = 1e-8;  

    for (int iter = 0; iter < MAX_ITER; iter++) {
        double dfL, dfR;
        double fL = f_wave(ps, pL, rL, aL, g, dfL);
        double fR = f_wave(ps, pR, rR, aR, g, dfR);

        double f  = fL + fR + (unR - unL);
        double df = dfL + dfR;

        double ps_new = ps - f / (df + 1e-300);
        ps_new = std::max(ps_new, 1e-12); 

        if (std::abs(ps_new - ps) / (0.5*(ps_new + ps)) < TOL) {
            ps = ps_new;
            break;
        }
        ps = ps_new;
    }

    // Определение оставшихся звездных пар-мов
    double dfL_dummy, dfR_dummy;
    double fL = f_wave(ps, pL, rL, aL, g, dfL_dummy);
    double fR = f_wave(ps, pR, rR, aR, g, dfR_dummy);

    double us  = 0.5*(unL + unR) + 0.5*(fR - fL);
    double rLs = rho_star(ps, pL, rL, g);
    double rRs = rho_star(ps, pR, rR, g);

 
    double SL, SHL, STL;  // скорость левой УВ или голова/хвост ЦВР
    double SR, SHR, STR;

    bool left_shock = (ps > pL);
    if (left_shock) {
        SL = unL - aL * std::sqrt((g + 1.0)/(2.0*g) * ps/pL + (g - 1.0)/(2.0*g));
        SHL = STL = SL;
    } else {
        double aLs = aL * std::pow(ps/pL, (g - 1.0)/(2.0*g));
        SHL = unL - aL;   
        STL = us  - aLs;  
        SL  = SHL;        
    }

    bool right_shock = (ps > pR);
    if (right_shock) {
        SR = unR + aR * std::sqrt((g+1.0)/(2.0*g) * ps/pR + (g-1.0)/(2.0*g));
        SHR = STR = SR;
    } else {
        double aRs = aR * std::pow(ps/pR, (g-1.0)/(2.0*g));
        SHR = unR + aR;
        STR = us  + aRs; 
        SR  = SHR;
    }

    // Регион и поток
    auto make_U = [&](double rho, double un, double ut_x, double ut_y, double p) -> Vec4 {
        double u = un*nx + ut_x;
        double v = un*ny + ut_y;
        double E = p/(g - 1.0) + 0.5*rho*(u*u + v*v);
        return {rho, rho*u, rho*v, E};
    };

    // Левее всех волн
    if (SL >= 0.0)
        return flux_n(UL, nx, ny, g);

    // Правее всех волн
    if (SR <= 0.0)
        return flux_n(UR, nx, ny, g);

    // Внутри левой ЦВР
    if (!left_shock && SHL < 0.0 && STL >= 0.0) {
        double un_fan = 2.0/(g+1.0) * (aL + (g-1.0)/2.0 * unL);
        double a_fan  = un_fan - unL + aL;  
        a_fan = aL - (g-1.0)/2.0*(un_fan - unL);
        double r_fan  = rL * std::pow(a_fan/aL, 2.0/(g-1.0));
        double p_fan  = pL * std::pow(a_fan/aL, 2.0*g/(g-1.0));
        Vec4 U_fan = make_U(r_fan, un_fan, utL_x, utL_y, p_fan);
        return flux_n(U_fan, nx, ny, g);
    }

    // Внутри правой ЦВР
    if (!right_shock && STR < 0.0 && SHR >= 0.0) {
        double un_fan = 2.0/(g+1.0) * (-aR + (g-1.0)/2.0 * unR);
        double a_fan  = aR + (g-1.0)/2.0*(unR - un_fan);
        double r_fan  = rR * std::pow(a_fan/aR, 2.0/(g-1.0));
        double p_fan  = pR * std::pow(a_fan/aR, 2.0*g/(g-1.0));
        Vec4 U_fan = make_U(r_fan, un_fan, utR_x, utR_y, p_fan);
        return flux_n(U_fan, nx, ny, g);
    }

    // Звёздная область L*
    if (us >= 0.0) {
        Vec4 ULs = make_U(rLs, us, utL_x, utL_y, ps);
        return flux_n(ULs, nx, ny, g);
    }

    // Звёздная область R*
    Vec4 URs = make_U(rRs, us, utR_x, utR_y, ps);
    return flux_n(URs, nx, ny, g);
}