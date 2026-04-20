#include <vector>
#include "Types.h"

State Minmod(const State& a,
             const State& b) {
    State s;

    for (int k = 0; k < NEQ; k++) {
        if (a[k]*b[k] <= 0.0)
            s[k] = 0.0;
        else
            s[k] = (std::abs(a[k]) < std::abs(b[k])) ? a[k] : b[k];
    }
    return s;
}

State Superbee(const State& a,
               const State& b) {
    State s;

    for (int k = 0; k < NEQ; k++) {
        if (a[k] * b[k] <= 0.0) 
            s[k] = 0.0;
        else {
            double sign = (a[k] > 0.0) ? 1.0 : -1.0;

            double abs_a = std::abs(a[k]);
            double abs_b = std::abs(b[k]);

            double val1 = std::min(2.0 * abs_a, abs_b);
            double val2 = std::min(abs_a, 2.0 * abs_b);

            s[k] = sign * std::max(val1, val2);
        }
    }
    return s;
}

State Vanleer(const State& a,
              const State& b) {
    State s;

    for (int k = 0; k < NEQ; k++) {
        if (a[k] * b[k] <= 0.0)
            s[k] = 0.0;
        else
            s[k] = 2.0 * a[k] * b[k] / (a[k] + b[k]);
    }
    return s;
}