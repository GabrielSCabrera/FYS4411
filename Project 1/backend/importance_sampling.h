#include "../wavefunctions/Psi.h"

#ifndef IMPORTANCE_SAMPLING_H
#define IMPORTANCE_SAMPLING_H

double* importance_sampling(Psi* PDF, int cycles, int N, double x_max, int equi_steps,
                    double dt, double D);

#endif
