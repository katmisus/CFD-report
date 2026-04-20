#ifndef _FLUXES_H_
#define _FLUXES_H_

#include "Types.h"

State PhysicalFlux(const State& W, int dir);

State HLLFlux(const State& WL,
              const State& WR,
              int dir);

State HLLCFlux(const State& WL,
               const State& WR,
               int dir);

State RusanovFlux(const State& WL,
                  const State& WR,
                  int dir);
                  
// double phi_lambda(double lambda, double delta = 1e-6);

State RoeFlux(const State& WL,
              const State& WR,
              int dir);

State OsherFlux(const State& WL,
                const State& WR,
                int dir);

State ExactFlux(const State& WL,
                const State& WR,
                int dir);

#endif