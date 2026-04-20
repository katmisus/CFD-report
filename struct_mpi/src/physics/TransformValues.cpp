#include <vector>
#include <cmath>
#include <iostream>
#include "RiemannSolver.h"
#include "Types.h"

extern int Nx, Ny, fict;
extern double gamm;

void ConvertWtoU(Field W, Field& U, int dir) {
	int ny = W.size();
	U.resize(ny);
    
	for (int j = 0; j < ny; j++) {
		int nx = W[j].size();
		U[j].resize(nx);
		for (int i = 0; i < nx; i++) {
            		double rho = W[j][i][0];
					double u   = W[j][i][1];
					double v   = W[j][i][2];
					double P   = W[j][i][3];

					double E = P / (rho * (gamm - 1.0)) + 0.5 * (u * u + v * v);

					U[j][i][0] = rho;
					U[j][i][1] = rho * u;
					U[j][i][2] = rho * v;
					U[j][i][3] = rho * E;
			
		}
	}
}

void ConvertWtoU(State W, State& U) {
    
	double rho = W[0];
	double u   = W[1];
	double v   = W[2];
	double P   = W[NEQ - 1];

	double E = P / (rho * (gamm - 1.0)) 
				+ 0.5 * (u * u + v * v);

	U[0] = rho;
	U[1] = rho * u;
	U[2] = rho * v;
	U[3] = rho * E;

	return;	
}

void ConvertUtoW(Field& W, Field U, int dir) {
	int ny = U.size();
    	W.resize(ny);

	for (int j = 0; j < ny; j++) {
        	int nx = U[j].size();
        	W[j].resize(nx);
        
        	for (int i = 0; i < nx; i++) {
            		double rho = std::max(1e-6, U[j][i][0]);           
            		double u = U[j][i][1] / rho;
					double v = U[j][i][2] / rho;
            		double E = U[j][i][3] / rho;                      
            
            
            		double P = std::max(1e-6, (gamm - 1.0) * rho * (E - 0.5 * (u * u + v * v)));
            
            		W[j][i][0] = rho;
            		W[j][i][1] = u;
            		W[j][i][2] = v;
            		W[j][i][3] = P;
        	}
    	}
}

