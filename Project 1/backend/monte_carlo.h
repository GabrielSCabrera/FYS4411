#include "../wavefunctions/Psi.h"

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

double* monte_carlo(Psi* PDF, int cycles, int N, double x_max, int equi_steps, double dt, double D);

#endif
