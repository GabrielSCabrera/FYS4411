#include "../backend/Psi.h"

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

double* monte_carlo(Psi PDF, int steps, int cycles, int N, double x_max, int equi_steps, double dt, double D);

void run_monte_carlo();

#endif
