#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "../wavefunctions/Psi_T.h"
#include "../wavefunctions/Psi_OB.h"



Monte_Carlo::Monte_Carlo(int N_particles, int dimensions, int x_limit, Psi* PDF_) {
	N = N_particles;
	dim = dimensions;
	x_max = x_limit;
	PDF = PDF_;
	step_length = 0.1;
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
	for (int i = 0; i < N; i++) {
			for (int j = 0; j < dim; j++) {
				R.set_raw(rand_double(-x_max, x_max), i, j);
			}
	}
	return R;
}


Mat Monte_Carlo::equilibriation(Mat R, int cycles) {
	Mat R_new(N, dim);
	double Psi = PDF->operator()(R);
	double Psi_new;
	double A; 			// acceptance ratio
	for (int i = 0; i < cycles; i++) {
		for (int j = 0; j < N; j++) {
			R_new = random_walk(R);
			// Get New Probability
			Psi_new = PDF->operator()(R_new);
			// Determine whether or not to accept movement
			A = acceptance_ratio(Psi_new, Psi, R_new, R, j);
			if (A > rand_double(0, 1)) {
				R = R_new;
				Psi = Psi_new;
			}
		}
	}
	return R;
}


/*
I think I would rather store the expectation values as variables in the 
class, and then return the latest set of coordiantes, R

*/

Mat Monte_Carlo::sample_energy(Mat R, int cycles) {
  set_to_zero();
  // Initialize Secondary particle Array
  Mat R_new(N, dim);
  double Psi_new, Psi;
  double* E_cycle = new double[cycles];         // Cycle-Wise Energy
  int accepted_moves = 0;
  double A; 			// acceptance ratio

  // #pragma omp parallel for
  // Monte-Carlo Cycles
  Psi = PDF->operator()(R);
  for (int i = 0; i < cycles; i++) {
      for (int j = 0; j < N; j++) {
  		R_new = random_walk(R);
        // Get New Probability
        Psi_new = PDF->operator()(R_new);
        // Determine whether or not to accept movement
        A = acceptance_ratio(Psi_new, Psi, R_new, R, j);
        if (A > rand_double(0, 1)) {
          R = R_new;
          Psi = Psi_new;
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
  accepted_moves_ratio = (double) accepted_moves/(cycles*N);

  /* Calculating Variance Iteratively

  !!!!!! GABRIEL
  I think the variance is just 
  	E*E - EE 
  or the other way around
  */
  for (int i = 0; i < cycles; i++) {
    variance += (E_cycle[i]*E_cycle[i] - output[0]*output[0]) *
                 (E_cycle[i] - output[0]) * (E_cycle[i] - output[0]);
  }
  // Cleanup
  delete [] E_cycle;
  return R;
}


Mat Monte_Carlo::sample_variational_derivatives(Mat R, int cycles) {
  set_to_zero();
  Mat R_new(N, dim);
  double Psi_new, Psi;
  double A; 			// acceptance ratio

  double psi_alpha = 0.0;		// derivative of Psi with respect to alpha
  double E_psi_alpha = 0.0;		// (derivative of Psi with respect to alpha)*E
  double psi_beta = 0.0;		// derivative of Psi with respect to beta
  double E_psi_beta = 0.0;		// (derivative of Psi with respect to beta)*E
  double temp, E_cycle;

  // #pragma omp parallel for
  // Monte-Carlo Cycles
  Psi = PDF->operator()(R);
  for (int i = 0; i < cycles; i++) {
      for (int j = 0; j < N; j++) {
		R_new = random_walk(R);
        // Get New Probability
        Psi_new = PDF->operator()(R_new);
        // Determine whether or not to accept movement
        A = acceptance_ratio(Psi_new, Psi, R_new, R, j);
        if (A > rand_double(0, 1)) {
          R = R_new;
          Psi = Psi_new;
        }
    } // End cycle

    E_cycle = PDF->energy(R);
    E += E_cycle;

    temp = PDF->grad_alpha(R);
    psi_alpha += temp;
    E_Psi_alpha += temp*E_cycle;

    temp = PDF->grad_beta(R);
    psi_beta += temp;
    E_psi_beta += temp*E_cycle;
  } // End Monte Carlo

  E /= cycles;
  psi_alpha /= cycles;
  psi_beta /= cycles;
  E_psi_alpha /= cycles;
  E_psi_beta /= cycles;

  grad_alpha = 2*(E_psi_alpha - psi_alpha*E);
  grad_beta = 2*(E_psi_beta  - psi_beta*E);
  return R;
}




