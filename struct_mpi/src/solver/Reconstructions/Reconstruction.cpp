#include <vector>
#include "Limiters.h"
#include "Types.h"
#include "BoundCond.h"

//extern double gamm, Lx, Ly, C1, C2;
extern double gamm;
extern int Nx, Ny, fict;
extern std::string method, rec_limiter;

void ReconstructGodunov(const Field& W,
                        Field& W_L,
                        Field& W_R,
                        int dir) {
    if (dir == 0) {
        for (int j = fict; j < Ny - 1 + fict; j++) {
            for (int i = fict; i < Nx + fict; i++) {
                W_L[i][j] = W[i - 1][j]; // левая ячейка
                W_R[i][j] = W[i][j];     // правая ячейка
            }
        }
    }

    else if (dir == 1) {
        for (int i = fict; i < Nx - 1 + fict; i++) {
            for (int j = fict; j < Ny + fict; j++) {
                W_L[i][j] = W[i][j - 1]; // нижняя ячейка
                W_R[i][j] = W[i][j];     // верхняя ячейка
            }
        }
    }
}


void ReconstructKolgan(const Field& W,
                       Field& W_L,
                       Field& W_R,
                       int dir) {
    if (dir == 0) {
        for (int j = fict; j < Ny - 1 + fict; j++) {
            for (int i = fict; i < Nx + fict; i++) {
                State dWm, dWp;

                for (int k = 0; k < NEQ; k++) {
                    dWm[k] = W[i-1][j][k] - W[i-2][j][k];
                    dWp[k] = W[i][j][k]   - W[i-1][j][k];
                }

                // CHECK: RECONSTRUCTION
                State slope;
                if (rec_limiter == "minmod")
                    slope = Minmod(dWm, dWp);
                else if (rec_limiter == "superbee")
                    slope = Superbee(dWm, dWp);
                else if (rec_limiter == "vanleer")
                    slope = Vanleer(dWm, dWp);

                for (int k = 0; k < NEQ; k++) {
                    W_L[i][j][k] =
                        W[i-1][j][k] + 0.5 * slope[k];

                    W_R[i][j][k] =
                        W[i][j][k]   - 0.5 * slope[k];
                }
            }
        }
    }

    else if (dir == 1) {
        for (int i = fict; i < Nx - 1 + fict; i++) {
            for (int j = fict; j < Ny + fict; j++) {
                State dWm, dWp;

                for (int k = 0; k < NEQ; k++) {
                    dWm[k] = W[i][j-1][k] - W[i][j-2][k];
                    dWp[k] = W[i][j][k]   - W[i][j-1][k];
                }

                State slope;
                if (rec_limiter == "minmod")
                    slope = Minmod(dWm, dWp);
                else if (rec_limiter == "superbee")
                    slope = Superbee(dWm, dWp);
                else if (rec_limiter == "vanleer")
                    slope = Vanleer(dWm, dWp);

                for (int k = 0; k < NEQ; k++) {
                    W_L[i][j][k] =
                        W[i][j-1][k] + 0.5 * slope[k];

                    W_R[i][j][k] =
                        W[i][j][k]   - 0.5 * slope[k];
                }
            }
        }
    }
}


void PredictorRodionov(Field& W_new, 
                    const Field& W, 
                    const Field& F_R, 
                    const Field& F_L,
                    const Field& G_R,
                    const Field& G_L, 		   
                    const std::vector<double>& x, 
                    const std::vector<double>& y,
                    double dt){

	for (int i = fict; i < Nx + fict - 1; i++) {
		double dx = x[i + 1] - x[i];

		for (int j = fict; j < Ny + fict - 1; j++) {
			double dy = y[j + 1] - y[j];
		
			W_new[i][j][0] = std::max(1e-7,
							W[i][j][0]
							- dt/dx * (F_L[i + 1][j][0] - F_R[i][j][0])
							- dt/dy * (G_L[i][j + 1][0] - G_R[i][j][0]));


			double mom_x_new = W[i][j][0] * W[i][j][1]
							- dt/dx * (F_L[i + 1][j][1] - F_R[i][j][1])
							- dt/dy * (G_L[i][j + 1][1] - G_R[i][j][1]);  

			W_new[i][j][1] = mom_x_new / W_new[i][j][0];
			if (std::abs(W_new[i][j][1]) < 1e-9) W_new[i][j][1] = 0.0;

			double mom_y_new = W[i][j][0] * W[i][j][2]
				- dt/dx * (F_L[i + 1][j][2] - F_R[i][j][2])
				- dt/dy * (G_L[i][j + 1][2] - G_R[i][j][2]);

			W_new[i][j][2] = mom_y_new / W_new[i][j][0];
			if (std::abs(W_new[i][j][2]) < 1e-9) W_new[i][j][2] = 0.0;
			
			double E_old =
					W[i][j][3] / (gamm - 1.0)
					+ 0.5 * W[i][j][0] *
					(W[i][j][1]*W[i][j][1] + W[i][j][2]*W[i][j][2]);

			double E_new =
					E_old
					- dt/dx * (F_L[i + 1][j][3] - F_R[i][j][3])
					- dt/dy * (G_L[i][j + 1][3] - G_R[i][j][3]);

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

State FluxRodionov(const State& W, int dir) {
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



void ReconstructRodionov(const Field& W,
                       Field& W_L,
                       Field& W_R,
                       const std::vector<double>& x, 
		               const std::vector<double>& y,
                       double dt,
                       int dir) {

    size_t Nx_tot = Nx + 2*fict - 1;
	size_t Ny_tot = Ny + 2*fict - 1; 
    Field W_tilde(Nx_tot, std::vector<State>(Ny_tot));   
    Field  Slope_X(Nx_tot + 1, std::vector<State>(Ny_tot)); 
    Field  Slope_Y(Nx_tot,  std::vector<State>(Ny_tot + 1));      
    Field  F_L(Nx_tot + 1, std::vector<State>(Ny_tot));
    Field  F_R(Nx_tot + 1, std::vector<State>(Ny_tot));
    Field  G_L(Nx_tot,  std::vector<State>(Ny_tot + 1));
    Field  G_R(Nx_tot,  std::vector<State>(Ny_tot + 1));
    State dWmx, dWpx;
    State dWmy, dWpy;
    for (int j = fict; j < Ny + fict; j++) {                            //сделал одной прогонкой, а то и так долго
        for (int i = fict; i < Nx + fict; i++) {

            if(j < Ny - 1 + fict){

                for (int k = 0; k < NEQ; k++) {
                    dWmx[k] = W[i-1][j][k] - W[i-2][j][k];
                    dWpx[k] = W[i][j][k]   - W[i-1][j][k];
                }

                if (rec_limiter == "minmod")
                    Slope_X[i][j] = Minmod(dWmx, dWpx);
                else if (rec_limiter == "superbee")
                    Slope_X[i][j] = Superbee(dWmx, dWpx);
                else if (rec_limiter == "vanleer")
                    Slope_X[i][j] = Vanleer(dWmx, dWpx);

                for (int k = 0; k < NEQ; k++) {
                    W_L[i][j][k] =
                        W[i-1][j][k] + 0.5 * Slope_X[i][j][k];

                    W_R[i][j][k] =
                        W[i][j][k]   - 0.5 * Slope_X[i][j][k];
                }
                F_L[i][j] = FluxRodionov(W_L[i][j], 0);
                F_R[i][j] = FluxRodionov(W_R[i][j], 0);
            }

            if(i < Nx - 1 + fict){

                for (int k = 0; k < NEQ; k++) {
                    dWmy[k] = W[i][j-1][k] - W[i][j-2][k];
                    dWpy[k] = W[i][j][k]   - W[i][j-1][k];
                }

                if (rec_limiter == "minmod")
                    Slope_Y[i][j] = Minmod(dWmy, dWpy);
                else if (rec_limiter == "superbee")
                    Slope_Y[i][j] = Superbee(dWmy, dWpy);
                else if (rec_limiter == "vanleer")
                    Slope_Y[i][j] = Vanleer(dWmy, dWpy);

                for (int k = 0; k < NEQ; k++) {
                    W_L[i][j][k] =
                        W[i][j-1][k] + 0.5 * Slope_Y[i][j][k];

                    W_R[i][j][k] =
                        W[i][j][k]   - 0.5 * Slope_Y[i][j][k];
                }

                G_L[i][j] = FluxRodionov(W_L[i][j], 1);
                G_R[i][j] = FluxRodionov(W_R[i][j], 1);
            }
        }
    }

    PredictorRodionov(W_tilde, W, F_R, F_L, G_R, G_L, x, y, dt);    //предиктор 

    for (int i = fict; i < Nx - 1 + fict; i++) {                        //полушаг по времени
        for (int j = fict; j < Ny + fict; j++) {
            for (int k = 0; k < NEQ; k++) {
               W_tilde[i][j][k] = 0.5 * (W_tilde[i][j][k] + W[i][j][k]);
            }
        }
    }

    BoundCond(W_tilde);                                                 //восстанавливаем граничные условия

    if (dir == 0) {                                                     //реконструкция
        for (int j = fict; j < Ny - 1 + fict; j++) {
            for (int i = fict; i < Nx + fict; i++) {
                for (int k = 0; k < NEQ; k++) {
                    W_L[i][j][k] =
                        W_tilde[i-1][j][k] + 0.5 * Slope_X[i][j][k];

                    W_R[i][j][k] =
                        W_tilde[i][j][k]   - 0.5 * Slope_X[i][j][k];
                }
            }
        }
    }

    else if (dir == 1) {
        for (int i = fict; i < Nx - 1 + fict; i++) {
            for (int j = fict; j < Ny + fict; j++) {
                for (int k = 0; k < NEQ; k++) {
                    W_L[i][j][k] =
                        W_tilde[i][j-1][k] + 0.5 * Slope_Y[i][j][k];

                    W_R[i][j][k] =
                        W_tilde[i][j][k]   - 0.5 * Slope_Y[i][j][k];
                }
            }
        }
    }
}



void Reconstruct(const Field& W,
                 Field& W_L,
                 Field& W_R,
                 const std::vector<double>& x, 
		         const std::vector<double>& y,
                 double dt,
                 int dir) {
	
	if (method == "Godunov")
        ReconstructGodunov(W, W_L, W_R, dir);

    else if (method == "Kolgan")
        ReconstructKolgan(W, W_L, W_R, dir);

    
    else if (method == "Rodionov")
        ReconstructRodionov(W, W_L, W_R, x, y, dt, dir);

    
	return;
}