#ifndef _RECONSTRUCTION_H_
#define _RECONSTRUCTION_H_

#include "Types.h"

void ReconstructGodunov(const Field& W,
                        Field& W_L,
                        Field& W_R,
                        int dir);

void ReconstructKolgan(const Field& W,
                       Field& W_L,
                       Field& W_R,
                       int dir);
void PredictorRodionov( Field& W_new, 
                        const Field& W, 
                        const Field& F_R, 
                        const Field& F_L,
                        const Field& G_R,
                        const Field& G_L, 		   
                        const std::vector<double>& x, 
                        const std::vector<double>& y,
                        double dt);
                        
State FluxRodionov(const State& W, int dir);

void ReconstructRodionov(const Field& W,
                         Field& W_L,
                         Field& W_R,
                         const std::vector<double>& x, 
                         const std::vector<double>& y,
                         double dt,
                         int dir);

void Reconstruct(const Field& W,
                 Field& W_L,
                 Field& W_R,
                 const std::vector<double>& x, 
		         const std::vector<double>& y,
                 double dt,
                 int dir);


#endif