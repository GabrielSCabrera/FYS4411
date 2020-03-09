#include "../wavefunctions/Psi.h"

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

double* monte_carlo(Psi* PDF, int N, int dim, double x_max, int cycles, int equi_steps);

#endif
