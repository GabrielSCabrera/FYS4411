#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

class Monte_Carlo {
protected:
	int N, dim;
	double L;
	double E, EE;
	double grad_alpha, grad_beta;
	double accepted_moves_ratio;
	double variance;
	double *E_cycles;
	void set_to_zero();
public:
	Psi* PDF;
	Monte_Carlo(Psi* trial_wave_function, int N_particles, int dimensions);
	~Monte_Carlo();
	double get_energy();
	double get_energy_mean();
	double get_grad_alpha();
	double get_grad_beta();
	double get_variance();
	double get_accepted_moves_ratio();
	void print_info();
	double rand_double(double min, double max);
	virtual double acceptance_ratio(Mat R_new, Mat R_old, int index) = 0;
	virtual Mat random_walk(Mat R, int index) = 0;
	Mat get_initial_R();
	Mat get_initial_R_no_overlap();
	Mat equilibriation(Mat R, int cycles);
	Mat sample_energy(Mat R, int cycles);
	Mat sample_variational_derivatives(Mat R, int cycles);
};

#endif

