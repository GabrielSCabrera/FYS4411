#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "../wavefunctions/Psi_T.h"
#include "../wavefunctions/Psi_OB.h"
#include "monte_carlo_class.h"

Monte_Carlo::Monte_Carlo(Psi* trial_wave_function, int N_particles, int dimensions) {
	N = N_particles;
	dim = dimensions;
	PDF = trial_wave_function;
  double hard_radius = trial_wave_function->get_a();
	step_length = 100*hard_radius;
	// x_max = 100*N*dim*hard_radius;
  x_max = 1;
}

// GETTERS

double Monte_Carlo::get_energy() {
  return E;
}

double Monte_Carlo::get_energy_mean() {
  return E/(N*dim);
}

double Monte_Carlo::get_grad_alpha() {
  return grad_alpha;
}

double Monte_Carlo::get_grad_beta() {
  return grad_beta;
}

double Monte_Carlo::get_variance() {
  return variance;
}

double Monte_Carlo::get_accepted_moves_ratio() {
  return accepted_moves_ratio;
}

void Monte_Carlo::print_info() {
  printf("E: %.6lf, var: %.6lf , acceptance: %.6lf\n", E/(N*dim), variance, accepted_moves_ratio);
}

double Monte_Carlo::rand_double(double min, double max) {
  /*
    Returns a randomly generated <double> in range [-N_max, N_max]
  */
  return (max-min)*(static_cast<double>(rand()) / RAND_MAX) + min;
}

void Monte_Carlo::set_to_zero() {
	E = 0.0;
	EE = 0.0;
	grad_alpha = 0.0;
	grad_beta = 0.0;
	variance = 0.0;
}

Mat Monte_Carlo::get_initial_R() {
	/*
	Returns a <Mat> instance of shape (N, dim) full of randomly generated values
	in the range [-x_max, x_max].
	*/
	Mat R(N, dim);
	for (int i = 0; i < N*dim; i++) {
      R.set_raw(rand_double(-x_max, x_max), i);
	}
	return R;
}

Mat Monte_Carlo::equilibriation(Mat R, int cycles) {
	Mat R_new(N, dim);
	double psi = PDF->operator()(R);
	double psi_new;
	double A; 			// acceptance ratio
	for (int i = 0; i < cycles; i++) {
		for (int j = 0; j < N; j++) {

			R_new = random_walk(R, j);

			// Get New Probability
			psi_new = PDF->operator()(R_new);

			// Determine whether or not to accept movement
			A = acceptance_ratio(psi_new, psi, R_new, R, j);
			if (A > rand_double(0, 1)) {
				R = R_new;
				psi = psi_new;
			}
		}
	}
	return R;
}

Mat Monte_Carlo::sample_energy(Mat R, int cycles) {
  set_to_zero();
  // Initialize Secondary particle Array
  Mat R_new(N, dim);
  double psi_new, psi;
  double* E_cycle = new double[cycles];         // Cycle-Wise Energy
  double accepted_moves = 0;
  double A; 			// acceptance ratio

  // #pragma omp parallel for
  // Monte-Carlo Cycles
  psi = PDF->operator()(R);
  for (int i = 0; i < cycles; i++) {
      for (int j = 0; j < N; j++) {

				R_new = random_walk(R, j);

				// Get New Probability
        psi_new = PDF->operator()(R_new);

				// Determine whether or not to accept movement
        printf("old :%lf new: %lf\n", psi, psi_new);
        A = acceptance_ratio(psi_new, psi, R_new, R, j);
        if (A > rand_double(0, 1)) {
          // printf("A: %lf\n", A);
          // printf("old :%lf new: %lf\n", psi, psi_new);
          R = R_new;
          psi = psi_new;
          accepted_moves++;
        }
    } // End cycle
    E_cycle[i] = PDF->energy(R);
    // Save Total Energy
    E += E_cycle[i];
    EE += E_cycle[i]*E_cycle[i];

  } // End Monte Carlo
  E /= cycles;
  EE /= cycles;
  accepted_moves_ratio = accepted_moves/(cycles*N);
  variance = E*E - EE;
  delete [] E_cycle;
  return R;
}

Mat Monte_Carlo::sample_variational_derivatives(Mat R, int cycles) {
  set_to_zero();
  Mat R_new(N, dim);
  double psi_new, psi;
  double A; 			// acceptance ratio

  double psi_alpha = 0.0;			//  derivative of Psi with respect to alpha
  double E_psi_alpha = 0.0;		// (derivative of Psi with respect to alpha)*E
  double psi_beta = 0.0;			//  derivative of Psi with respect to beta
  double E_psi_beta = 0.0;		// (derivative of Psi with respect to beta)*E
  double temp, E_cycle;

  double accepted_moves = 0;
  // Monte-Carlo Cycles
  psi = PDF->operator()(R);
  for (int i = 0; i < cycles; i++) {
      for (int j = 0; j < N; j++) {

				R_new = random_walk(R, j);

				// Get New Probability
        psi_new = PDF->operator()(R_new);

				// Determine whether or not to accept movement
        A = acceptance_ratio(psi_new, psi, R_new, R, j);
        if (A > rand_double(0, 1)) {
          R = R_new;
          psi = psi_new;
          accepted_moves++;
        }
    } // End cycle

    E_cycle = PDF->energy(R);
    E += E_cycle;

    temp = PDF->grad_alpha(R);
    psi_alpha += temp;
    E_psi_alpha += temp*E_cycle;

    temp = PDF->grad_beta(R);
    psi_beta += temp;
    E_psi_beta += temp*E_cycle;
  } // End Monte Carlo

  E /= cycles;
  psi_alpha /= cycles;
  psi_beta /= cycles;
  E_psi_alpha /= cycles;
  E_psi_beta /= cycles;

  accepted_moves_ratio = accepted_moves/(cycles*N);

  grad_alpha = 2*(E_psi_alpha - psi_alpha*E);
  grad_beta = 2*(E_psi_beta  - psi_beta*E);
  return R;
}
