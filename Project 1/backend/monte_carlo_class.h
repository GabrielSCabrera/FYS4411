#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"

#ifndef MONTE_CARLO_CLASS_H
#define MONTE_CARLO_CLASS_H

class Monte_Carlo {
protected:
	/*
	We would need getters for the expectation values
	and setters for step length and max value of x
	if we have a default
	*/
	int N, dim;
	double step_length, x_max;
	double E, EE;
	double grad_alpha, grad_beta;
	double accepted_moves_ratio;
	double variance;
public:
	Monte_Carlo(int, int, int, Psi*);
	Psi* PDF;
	void set_to_zero();
	double rand_double(double min, double max);
	virtual Mat random_walk(Mat R) = 0;
	virtual double acceptance_ratio(double Psi_new, double Psi_old, Mat R_new, Mat R_old, int index) = 0;
	Mat get_initial_R();
	Mat equilibriation(Mat R, int cycles);
	Mat sample_energy(Mat R, int cycles);
	Mat sample_variational_derivatives(Mat R, int cycles);
};

#endif


/*
Metropolis and importance sampling can be subclasses
They only need those two virtual methods + constructors
*/
