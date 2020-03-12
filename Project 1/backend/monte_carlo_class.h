#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"

#ifndef MONTE_CARLO_CLASS_H
#define MONTE_CARLO_CLASS_H

class Monte_Carlo {
protected:
	int N, dim;
	double step_length, x_max;
	double E, EE;
	double grad_alpha, grad_beta;
	double accepted_moves_ratio;
	double variance;
public:
	Psi* PDF;
	Monte_Carlo(Psi* trial_wave_function, int N_particles, int dimensions);
	void set_to_zero();
	double rand_double(double min, double max);
	virtual double acceptance_ratio(double Psi_new, double Psi_old, Mat R_new, Mat R_old, int index) = 0;
	virtual Mat random_walk(Mat R) = 0;
	Mat get_initial_R();
	Mat equilibriation(Mat R, int cycles);
	Mat sample_energy(Mat R, int cycles);
	Mat sample_variational_derivatives(Mat R, int cycles);
};

#endif

