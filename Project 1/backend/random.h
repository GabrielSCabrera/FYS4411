#include "../matpak/Mat.h"
#include "Psi.h"

#ifndef RANDOM_H
#define RANDOM_H

double rand_double(double min, double max);

double normal_distribution(double N_max, double mu = 1, double sigma = 1);

Mat random_particles(int N, double min, double max);

Mat random_walk(Psi PDF, Mat P, double step_size, long idx, double Ddt);

#endif
