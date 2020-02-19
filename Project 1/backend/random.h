#include "../matpak/Mat.h"

#ifndef RANDOM_H
#define RANDOM_H

double rand_double(double N_max);

double normal_distribution(double N_max, double mu = 1, double sigma = 1);

Mat random_particles(int N, double x_max);

Mat random_walk(Mat P, double step_size);

#endif
