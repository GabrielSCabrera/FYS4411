#include "../wavefunctions/Psi.h"

#ifndef IMPORTANCE_SAMPLING_H
#define IMPORTANCE_SAMPLING_H

double* monte_carlo(Psi* PDF, int N, int dim, double x_max, int cycles, int equi_steps);

#endif
