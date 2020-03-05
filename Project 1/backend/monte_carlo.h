#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

double* monte_carlo(int steps, int cycles, int N, double x_max, double alpha,
                    double beta, double a, double omega, double omega_z,
                    int equi_steps, double dt, double D);

void run_monte_carlo();

#endif
