#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include <fstream>

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

class Monte_Carlo {
protected:
	int N, dim;
	int N_rho;
	double L;
	double E, EE;
	double E_alpha;		// dE/dalpha
	double accepted_moves_ratio;
	double variance;
	double *E_cycles;
	int** rho;
	int MC_cycles = 1;
	void set_to_zero();
public:
	double E_2alpha;	// d^2E/dalpha^2
	Psi* bose;
	Monte_Carlo(Psi* trial_wave_function, int N_particles, int dimensions);
	~Monte_Carlo();
	double get_energy();
	double get_energy_mean();
	double get_grad_alpha();
	double get_variance();
	double get_accepted_moves_ratio();
	double* get_E_cycles();
	int** get_rho(int *length_rho); 
	void print_info();
	void copy_step(Mat* from, Mat* to, int index);
	double random_normal_distribution();
	virtual double acceptance_ratio(Mat* R_new, Mat* R_old, int k) = 0;
	virtual void random_walk(Mat* R, int index) = 0;
	Mat get_initial_R();
	Mat get_initial_R_no_overlap();
	Mat equilibriation(Mat R, int cycles);
	Mat sample_energy(Mat R, int cycles);
	Mat sample_variational_derivatives(Mat R, int cycles);
	Mat one_body_density(Mat R, int cycles, int anticipated_max);
};

#endif

