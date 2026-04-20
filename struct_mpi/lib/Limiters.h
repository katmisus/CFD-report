#ifndef _LIMITERS_H_
#define _LIMITERS_H_

#include "Types.h"

State Minmod(const State& a,
             const State& b);

State Superbee(const State& a,
               const State& b);

State Vanleer(const State& a,
              const State& b);

#endif