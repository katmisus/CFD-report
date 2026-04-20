#include <cmath>
#include <vector>
#include <iostream>
#include "FLIC.h"
#include "FileProcessing.h"
#include "BoundCond.h"
#include "Init.h"
#include "Types.h"

extern double gamm, Lx, Ly, C1, C2;
extern int Nx, Ny, fict;

// CHECK: FLIC_LAGRANGE
void FLIC_L(const Field& W,
            Field& W_tilde,
            const std::vector<double>& x,
            const std::vector<double>& y,
            double dt,
            int dir) {

    if (dir == 0) {
        for (int i = fict; i < Nx + fict - 1; i++) {
            for (int j = fict; j < Ny + fict; j++) {
                double dx = x[i + 1] - x[i];

                double rho = W[i][j][0];
                double u = W[i][j][1];
                double v = W[i][j][2];
                double P = W[i][j][NEQ - 1];

                double P_ip = 0.5 * (W[i][j][NEQ - 1] + W[i + 1][j][NEQ - 1]);
                double P_im = 0.5 * (W[i - 1][j][NEQ - 1] + W[i][j][NEQ - 1]);

                W_tilde[i][j][0] = rho;
    
                
                double u_tilde = u - dt / rho * (P_ip - P_im) / dx;
                W_tilde[i][j][1] = u_tilde;

                W_tilde[i][j][2] = v; 

                double E = P / (gamm - 1.0) + 0.5 * rho * (u * u + v * v);
                double u_ip = 0.5 * (W[i][j][1] + W[i + 1][j][1]);
                double u_im = 0.5 * (W[i - 1][j][1] + W[i][j][1]);
                double E_tilde = E - dt * (P_ip * u_ip - P_im * u_im) / dx;
                double kinetic_tilde = 0.5 * rho * (u_tilde * u_tilde + v * v);
                double P_tilde = (gamm - 1.0) * (E_tilde - kinetic_tilde);
        
                W_tilde[i][j][NEQ - 1] = std::max(1e-8, P_tilde);
            }
        }
    }

    if (dir == 1) {
        for (int i = fict; i < Nx + fict; i++) {
            for (int j = fict; j < Ny + fict - 1; j++) {
                double dy = y[j + 1] - y[j];

                double rho = W[i][j][0];
                double u = W[i][j][1];
                double v = W[i][j][2];
                double P = W[i][j][NEQ - 1];

                double P_jp = 0.5 * (W[i][j][NEQ - 1] + W[i][j + 1][NEQ - 1]);
                double P_jm = 0.5 * (W[i][j - 1][NEQ - 1] + W[i][j][NEQ - 1]);

                W_tilde[i][j][0] = rho;
                W_tilde[i][j][1] = u;

                double dPdy = (P_jp - P_jm) / dy;
                if (std::abs(dPdy) < 1e-14) dPdy = 0.0;
                double v_tilde = v - dt / rho * dPdy;

                W_tilde[i][j][2] = v_tilde;

                double E = P / (gamm - 1.0) + 0.5 * rho * (u * u + v * v);
                double v_jp = 0.5 * (W[i][j][2] + W[i][j + 1][2]);
                double v_jm = 0.5 * (W[i][j - 1][2] + W[i][j][2]);

                double dPVdy = (P_jp * v_jp - P_jm * v_jm) / dy;
                if (std::abs(dPVdy) < 1e-14) dPVdy = 0.0;
                
                double E_tilde = E - dt  * dPVdy;
                double kinetic_tilde = 0.5 * rho * (u * u + v_tilde * v_tilde);
                double P_tilde = (gamm - 1.0) * (E_tilde - kinetic_tilde);

                W_tilde[i][j][NEQ - 1] = std::max(1e-8, P_tilde);
            }
        }
    }
}


// CHECK: FLIC_EULER
void FLIC_E(const Field& W_tilde,
            Field& W_new,
            const std::vector<double>& x,
            const std::vector<double>& y,
            double dt,
            int dir) {

    const double vel_eps = 1e-12;

    if (dir == 0) {  
        for (int i = fict; i < Nx + fict - 1; i++) {
            for (int j = fict; j < Ny + fict; j++) {
                double dx = x[i + 1] - x[i];

                double rho_flux  = 0.0;
                double momx_flux = 0.0;
                double momy_flux = 0.0;
                double E_flux    = 0.0;

                // upwind на границе i+1/2 
                double u_face = 0.5 * (W_tilde[i][j][1] + W_tilde[i+1][j][1]);
                if (std::abs(u_face) > vel_eps) {
                    int up = (u_face > 0.0) ? i : i+1;

                    rho_flux  = W_tilde[up][j][0] * W_tilde[up][j][1];
                    momx_flux = rho_flux * W_tilde[up][j][1];
                    momy_flux = rho_flux * W_tilde[up][j][2];

                    double E = W_tilde[up][j][NEQ - 1] / (gamm - 1.0) + 
                                0.5 * W_tilde[up][j][0] * 
                                (W_tilde[up][j][1]*W_tilde[up][j][1] + W_tilde[up][j][2] * W_tilde[up][j][2]);

                    E_flux = W_tilde[up][j][1] * E;
                }


                double rho_flux_m  = 0.0;
                double momx_flux_m = 0.0;
                double momy_flux_m = 0.0;
                double E_flux_m    = 0.0;
                // аналогично для i-1/2 
                double u_face_m = 0.5 * (W_tilde[i - 1][j][1] + W_tilde[i][j][1]);
                if (std::abs(u_face_m) > vel_eps) {
                    int up_m = (u_face_m > 0.0) ? i-1 : i;

                    rho_flux_m  = W_tilde[up_m][j][0] * W_tilde[up_m][j][1];
                    momx_flux_m = rho_flux_m * W_tilde[up_m][j][1];
                    momy_flux_m = rho_flux_m * W_tilde[up_m][j][2];

                    double E_m = W_tilde[up_m][j][NEQ - 1] / (gamm - 1.0) + 0.5 * W_tilde[up_m][j][0] * (W_tilde[up_m][j][1]*W_tilde[up_m][j][1] + W_tilde[up_m][j][2]*W_tilde[up_m][j][2]);

                    E_flux_m = W_tilde[up_m][j][1] * E_m;
                }
                // CHECK: FLIC_CONSERV
                double rho_new = W_tilde[i][j][0] - dt/dx * (rho_flux - rho_flux_m);

                double momx_new = W_tilde[i][j][0] * W_tilde[i][j][1] - dt/dx * (momx_flux - momx_flux_m);

                double momy_new = W_tilde[i][j][0] * W_tilde[i][j][2] - dt/dx * (momy_flux - momy_flux_m);

                double E_new = (W_tilde[i][j][NEQ - 1] / (gamm - 1.0) + 0.5 * W_tilde[i][j][0] * (W_tilde[i][j][1]*W_tilde[i][j][1] + W_tilde[i][j][2]*W_tilde[i][j][2])) - dt/dx * (E_flux - E_flux_m);
                
                rho_new = std::max(1e-8, rho_new);

                double u_new = momx_new / rho_new;
                double v_new = momy_new / rho_new;

                double kinetic = 0.5 * rho_new * (u_new * u_new + v_new * v_new);
                double P_new = (gamm - 1.0) * (E_new - kinetic);

                if (std::abs(u_new) < 1e-9)
                    u_new = 0.0;
                if (std::abs(v_new) < 1e-9)
                    v_new = 0.0;

                W_new[i][j][0] = rho_new;
                W_new[i][j][1] = u_new;
                W_new[i][j][2] = v_new;
                W_new[i][j][3] = std::max(1e-8, P_new);
                
            }
        }
    }

    if (dir == 1) {   
        for (int i = fict; i < Nx + fict; i++) {
            for (int j = fict; j < Ny + fict - 1; j++) {
                double dy = y[j + 1] - y[j];

                double rho_flux = 0.0;
                double momx_flux = 0.0;
                double momy_flux = 0.0;
                double E_flux = 0.0;

                double v_face = 0.5 * (W_tilde[i][j][2] + W_tilde[i][j+1][2]);
                if (std::abs(v_face) >  vel_eps) {
                    int up = (v_face > 0.0) ? j : j+1;

                    rho_flux  = W_tilde[i][up][0] * W_tilde[i][up][2];
                    momx_flux = rho_flux * W_tilde[i][up][1];
                    momy_flux = rho_flux * W_tilde[i][up][2];

                    double E = W_tilde[i][up][3] / (gamm - 1.0) + 0.5 * W_tilde[i][up][0] * (W_tilde[i][up][1]*W_tilde[i][up][1] + W_tilde[i][up][2]*W_tilde[i][up][2]);

                    E_flux = W_tilde[i][up][2] * E;
                }

                double rho_flux_m = 0.0;
                double momx_flux_m = 0.0;
                double momy_flux_m = 0.0;
                double E_flux_m = 0.0;

                double v_face_m = 0.5 * (W_tilde[i][j-1][2] + W_tilde[i][j][2]);
                if (std::abs(v_face_m) >  vel_eps) {
                    int up_m = (v_face_m > 0.0) ? j-1 : j;

                    rho_flux_m  = W_tilde[i][up_m][0] * W_tilde[i][up_m][2];
                    momx_flux_m = rho_flux_m * W_tilde[i][up_m][1];
                    momy_flux_m = rho_flux_m * W_tilde[i][up_m][2];

                    double E_m = W_tilde[i][up_m][NEQ - 1] / (gamm - 1.0) + 0.5 * W_tilde[i][up_m][0] * (W_tilde[i][up_m][1]*W_tilde[i][up_m][1] + W_tilde[i][up_m][2]*W_tilde[i][up_m][2]);

                    E_flux_m = W_tilde[i][up_m][2] * E_m;
                }
                // CHECK: FLIC_CONSERV
                double rho_new = W_tilde[i][j][0] - dt/dy * (rho_flux - rho_flux_m);

                double momx_new = W_tilde[i][j][0] * W_tilde[i][j][1] - dt/dy * (momx_flux - momx_flux_m);

                double momy_new = W_tilde[i][j][0] * W_tilde[i][j][2] - dt/dy * (momy_flux - momy_flux_m);

                double E_new = (W_tilde[i][j][3] / (gamm - 1.0) + 0.5 * W_tilde[i][j][0] * (W_tilde[i][j][1]*W_tilde[i][j][1] + W_tilde[i][j][2]*W_tilde[i][j][2])) - dt/dy * (E_flux - E_flux_m);
                rho_new = std::max(1e-8, rho_new);

                double u_new = momx_new / rho_new;
                double v_new = momy_new / rho_new;

                double kinetic = 0.5 * rho_new * (u_new * u_new + v_new * v_new);
                double P_new = (gamm - 1.0) * (E_new - kinetic);

                W_new[i][j][0] = rho_new;
                if (std::abs(u_new) < 1e-9)
                    u_new = 0.0;
                if (std::abs(v_new) < 1e-9)
                    v_new = 0.0;
                
                W_new[i][j][1] = u_new;
                W_new[i][j][2] = v_new;
                W_new[i][j][3] = std::max(1e-8, P_new);
            }
        }
    }
}


void FLIC(Field& W_new, 
		  const Field& W,
          const std::vector<double>& x, 
		  const std::vector<double>& y, 
		  double dt) {

    static bool swap = false;   // чередование направлений
     
    Field W_tilde = W;
    Field W_tmp = W;
    Field W_tmp_2 = W;

    double dt_half = 0.5 * dt;

    // По поводу swap - я не заметила, чтобы он что-то менял.
    // FLIC_L(W, W_tilde, x, y, dt_half, 0);
    // BoundCond(W_tilde);
    // FLIC_E(W_tilde, W_tmp, x, y, dt_half, 0);
    // BoundCond(W_tmp);

    // FLIC_L(W_tmp, W_tilde, x, y, dt, 1);
    // BoundCond(W_tilde);
    // FLIC_E(W_tilde, W_tmp_2, x, y, dt, 1);
    // BoundCond(W_tmp_2);

    // FLIC_L(W_tmp_2, W_tilde, x, y, dt_half, 0);
    // BoundCond(W_tilde);
    // FLIC_E(W_tilde, W_new, x, y, dt_half, 0);
    // BoundCond(W_tmp_2);

    if (!swap) {
        FLIC_L(W, W_tilde, x, y, dt_half, 0);
        BoundCond(W_tilde);
        // SaveFieldToCSV(W_tilde, x, y, dt,"data_track/W_tilde0_x.csv");
        FLIC_E(W_tilde, W_tmp, x, y, dt_half, 0);
        BoundCond(W_tmp);
        // SaveFieldToCSV(W_tmp, x, y, dt, "data_track/W_tmp0_x.csv");

        FLIC_L(W_tmp, W_tilde, x, y, dt, 1);
        BoundCond(W_tilde);
        // SaveFieldToCSV(W_tilde, x, y, dt,"data_track/W_tilde1_x.csv");
        FLIC_E(W_tilde, W_tmp_2, x, y, dt, 1);
        BoundCond(W_tmp_2);
        // SaveFieldToCSV(W_tmp_2, x, y, dt, "data_track/W_tmp1_x.csv");

        FLIC_L(W_tmp_2, W_tilde, x, y, dt_half, 0);
        BoundCond(W_tilde);
        // SaveFieldToCSV(W_tilde, x, y, dt,"data_track/W_tilde2_x.csv");
        FLIC_E(W_tilde, W_new, x, y, dt_half, 0);
        BoundCond(W_new);

    } else {
        FLIC_L(W, W_tilde, x, y, dt_half, 1);
        BoundCond(W_tilde);
        // SaveFieldToCSV(W_tilde, x, y, dt,"data_track/W_tilde0_y.csv");
        FLIC_E(W_tilde, W_tmp, x, y, dt_half, 1);
        BoundCond(W_tmp);
        // SaveFieldToCSV(W_tmp, x, y, dt, "data_track/W_tmp0_y.csv");

        FLIC_L(W_tmp, W_tilde, x, y, dt, 0);
        BoundCond(W_tilde);
        // SaveFieldToCSV(W_tilde, x, y, dt,"data_track/W_tilde1_y.csv");
        FLIC_E(W_tilde, W_tmp_2, x, y, dt, 0);
        BoundCond(W_tmp_2);
        // SaveFieldToCSV(W_tmp_2, x, y, dt, "data_track/W_tmp1_y.csv");

        FLIC_L(W_tmp_2, W_tilde, x, y, dt_half, 1);
        BoundCond(W_tilde);
        // SaveFieldToCSV(W_tilde, x, y, dt,"data_track/W_tilde2_y.csv");
        FLIC_E(W_tilde, W_new, x, y, dt_half, 1);
        BoundCond(W_new);            
    }

    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            if (std::abs(W_new[i][j][1]) < 1e-9) W_new[i][j][1] = 0.0;
            if (std::abs(W_new[i][j][2]) < 1e-9) W_new[i][j][2] = 0.0;
        }
    }

    swap = !swap;

    

    return;
}

