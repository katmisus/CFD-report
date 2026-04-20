#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
//#include "BoundCond.h"
#include "Fluxes.h"
#include "Types.h"
#include "FileProcessing.h"
#include "Reconstruction.h"
#include "FLIC.h"
#include "Mader.h"


extern double CFL, gamm, Lx, Ly, C1, C2, T_init, R_gas, M, P_min, E_act, Z_freq, VISC, MINWT, GASW, MINGRHO;
extern int Nx, Ny, fict;

extern std::string method, solver, time_method;
//extern std::string high_order_method, TVD_solver;
//extern bool Viscous_flag, TVD_flag;


State ComputeNumericalFlux(const State& WL,
                           const State& WR,
                           int dir) {
    if (solver == "Exact")
        return ExactFlux(WL, WR, dir);

    else if (solver == "HLL")
        return HLLFlux(WL, WR, dir);

    else if (solver == "HLLC")
        return HLLCFlux(WL, WR, dir);

	else if (solver == "Rusanov")
        return RusanovFlux(WL, WR, dir);

	else if (solver == "Osher")
        return OsherFlux(WL, WR, dir);

	else
        return RoeFlux(WL, WR, dir);
}

void ComputeFluxes(const Field& W,
					Field& F,
					const std::vector<double>& x, 
					const std::vector<double>& y,
					double dt,
					int dir) {

	size_t Nx_tot = Nx + 2*fict - 1;
	size_t Ny_tot = Ny + 2*fict - 1;

    if (dir == 0) {
		Field W_L(Nx_tot + 1,
                  std::vector<State>(Ny_tot));
        Field W_R(Nx_tot + 1,
                  std::vector<State>(Ny_tot));
												//добавил x, y, dt для родионова
        Reconstruct(W, W_L, W_R, x, y, dt, 0); // Сюда вшивать реконструкции Колгана/Родионова/ENO
        for (int i = fict; i < Nx + fict; i++) {
            for (int j = fict; j < Ny - 1 + fict; j++) {
            	F[i][j] = ComputeNumericalFlux(W_L[i][j], W_R[i][j], dir);
			}
		}
	}

	if (dir == 1) {
		Field W_L(Nx_tot,
                  std::vector<State>(Ny_tot + 1));
        Field W_R(Nx_tot,
                  std::vector<State>(Ny_tot + 1));

        Reconstruct(W, W_L, W_R, x, y, dt, 1);
        for (int i = fict; i < Nx - 1 + fict; i++) {
            for (int j = fict; j < Ny + fict; j++) {
            	F[i][j] = ComputeNumericalFlux(W_L[i][j], W_R[i][j], dir);
			}
		}
	}

}


void Euler(const Field& W,
		   Field& W_new,
		   const std::vector<double>& x, 
		   const std::vector<double>& y,
		   double dt) {

	size_t Nx_tot = Nx + 2*fict - 1;
	size_t Ny_tot = Ny + 2*fict - 1;

	Field F(Nx_tot + 1, std::vector<State>(Ny_tot));
	Field G(Nx_tot,     std::vector<State>(Ny_tot + 1));
	
	ComputeFluxes(W, F, x, y, dt, 0);
	ComputeFluxes(W, G, x, y, dt, 1);
	

	for (int i = fict; i < Nx + fict - 1; i++) {
		double dx = x[i + 1] - x[i];

		for (int j = fict; j < Ny + fict - 1; j++) {
			double dy = y[j + 1] - y[j];
		
			W_new[i][j][0] = std::max(1e-7,
							W[i][j][0]
							- dt/dx * (F[i + 1][j][0] - F[i][j][0])
							- dt/dy * (G[i][j + 1][0] - G[i][j][0]));


			double mom_x_new = W[i][j][0] * W[i][j][1]
							- dt/dx * (F[i + 1][j][1] - F[i][j][1])
							- dt/dy * (G[i][j + 1][1] - G[i][j][1]);  

			W_new[i][j][1] = mom_x_new / W_new[i][j][0];
			if (std::abs(W_new[i][j][1]) < 1e-9) W_new[i][j][1] = 0.0;

			double mom_y_new = W[i][j][0] * W[i][j][2]
				- dt/dx * (F[i + 1][j][2] - F[i][j][2])
				- dt/dy * (G[i][j + 1][2] - G[i][j][2]);

			W_new[i][j][2] = mom_y_new / W_new[i][j][0];
			if (std::abs(W_new[i][j][2]) < 1e-9) W_new[i][j][2] = 0.0;
			
			double E_old =
					W[i][j][3] / (gamm - 1.0)
					+ 0.5 * W[i][j][0] *
					(W[i][j][1]*W[i][j][1] + W[i][j][2]*W[i][j][2]);

			double E_new =
					E_old
					- dt/dx * (F[i + 1][j][3] - F[i][j][3])
					- dt/dy * (G[i][j + 1][3] - G[i][j][3]);

			double kinetic_new =
					0.5 * W_new[i][j][0] *
					(W_new[i][j][1]*W_new[i][j][1] +
					W_new[i][j][2]*W_new[i][j][2]);

			double P_new =
					(gamm - 1.0) * (E_new - kinetic_new);

			W_new[i][j][NEQ - 1] = std::max(1e-6, P_new);

		}
    }
}


void UpdateArrays(Field& W, 
				  Field W_new,
				  std::vector<double> x, 
				  std::vector<double> y,
				  double dt) {
	
	if (method == "FLIC") {
        FLIC(W_new, W, x, y, dt);
	} else if (method == "Mader") {
        Mader(W_new, W, x, y, dt);
	} else if (time_method == "Euler") {
        Euler(W, W_new, x, y, dt);
	}
	// Для остальных методов
	// else if (time_method == "RK3") {
	// 	RK3(W_new, W, method, solver, func, fict, N + fict - 1, x, dt, Viscous_flag);
	// } 

	std::swap(W, W_new);
    // for (int i = fict; i < Nx + fict - 1; i++) {
    //     for (int j = fict; j < Ny + fict - 1; j++) {
	// 		W[i][j] = W_new[i][j];
    //     }
    // }
}