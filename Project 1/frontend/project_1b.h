#include "../matpak/Mat.h"

#ifndef PROJECT_1B_H
#define PROJECT_1B_H

double random(double N_max);

double normal_distribution(double N_max, double mu = 1, double sigma = 1);

Mat random_particles(int N, double x_max);

void project_1b_main();

#endif
