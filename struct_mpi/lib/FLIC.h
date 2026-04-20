#ifndef _FLIC_H_
#define _FLIC_H

#include <functional>
#include "Types.h"

void FLIC_L(const Field& W,
            Field& W_tilde,
            const std::vector<double>& x,
            const std::vector<double>& y,
            double dt,
            int dir);

void FLIC_E(const Field& W_tilde,
            Field& W_new,
            const std::vector<double>& x,
            const std::vector<double>& y,
            double dt,
            int dir);

void FLIC(Field& W_new, 
		  const Field& W,
          const std::vector<double>& x, 
		  const std::vector<double>& y, 
		  double dt);

#endif