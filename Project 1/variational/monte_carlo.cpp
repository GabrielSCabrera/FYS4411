#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include <fstream>

Monte_Carlo::Monte_Carlo(Psi* trial_wave_function, int N_particles, int dimensions) {
	N = N_particles;
	dim = dimensions;
	PDF = trial_wave_function;
	L = 10;
  E_cycles = new double [1];
  PDF->L = 0.5*L;
  if (dim < 3) {
    PDF->update_beta(1.0);
    PDF->update_gamma(1.0);
  }
}

// DESTRUCTOR
Monte_Carlo::~Monte_Carlo() {
  delete[] E_cycles;
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
      R.set_raw(rand_double(0.01*L, 0.99*L), i);
	}

	return R;
}

Mat Monte_Carlo::get_initial_R_no_overlap() {
  Mat R(N, dim);
  double r_ij, dx;
  bool overlap;
  double hard_radius = PDF->get_a();
  hard_radius *= hard_radius;
  for (int i = 0; i < N; i++) {
    overlap = true;
    while (overlap) {
      for (int l = 0; l < dim; l++) {
        R.set(rand_double(0.01*L, 0.99*L), i, l);
      }
      overlap = false;
      for (int j = 0; j < i; j++) {
        r_ij = 0.0;
        for (int l = 0; l < dim; l++) {
          dx = R.get(j, l) - R.get(i, l);
          r_ij += dx*dx;
        }
        if (r_ij <= hard_radius) {
          overlap = true;
          continue;
        }
      }
    }
  }
  return R;
}

Mat Monte_Carlo::equilibriation(Mat R, int cycles) {
	Mat R_new(N, dim);
	double A; 			// acceptance ratio
	for (int i = 0; i < cycles; i++) {
		for (int j = 0; j < N; j++) {
			R_new = random_walk(R, j);

			// Determine whether or not to accept movement
			A = acceptance_ratio(R_new, R, j);
			if (A > rand_double(0.0, 1.0)) {
				R = R_new;
			}
		}
	}
	return R;
}

Mat Monte_Carlo::sample_energy(Mat R, int cycles) {
  set_to_zero();
  delete[] E_cycles;
  // Initialize Secondary particle Array
  Mat R_new(N, dim);
  E_cycles = new double[cycles];         // Cycle-Wise Energy
  MC_cycles = cycles;
  double accepted_moves = 0;
  double A;															// Acceptance Ratio

  double r;

  // Monte-Carlo Cycles
  for (int i = 0; i < cycles; i++) {
      for (int j = 0; j < N; j++) {
				R_new = random_walk(R, j);

				// Determine whether or not to accept movement
        A = acceptance_ratio(R_new, R, j);
        r = rand_double(0.0, 1.0);
        printf("(%5.2lf,%5.2lf)", A, r);
        if (A > r) {
          R = R_new;
          accepted_moves++;
          printf("+ ");
        } else {
          printf("- ");
        }
    } // End cycle
    E_cycles[i] = PDF->energy(R);
    printf("\t-> %8.3lf : %6.3e \n", E_cycles[i], PDF->operator()(R));
 
    // Save Total Energy
    E += E_cycles[i];
    EE += (E_cycles[i]*E_cycles[i])/cycles;

  } // End Monte Carlo
  E /= cycles;
  //EE /= cycles;
  accepted_moves_ratio = accepted_moves/(cycles*N);
  variance = EE - E*E;
  return R;
}

Mat Monte_Carlo::sample_variational_derivatives(Mat R, int cycles) {
  set_to_zero();
  Mat R_new(N, dim);
  double A; 									// Acceptance Ratio

  double psi_alpha = 0.0;			//  derivative of Psi with respect to alpha
  double E_psi_alpha = 0.0;		// (derivative of Psi with respect to alpha)*E
  double psi_beta = 0.0;			//  derivative of Psi with respect to beta
  double E_psi_beta = 0.0;		// (derivative of Psi with respect to beta)*E
  double temp, E_cycle;

  double accepted_moves = 0;

	// Monte-Carlo Cycles
  for (int i = 0; i < cycles; i++) {
      for (int j = 0; j < N; j++) {
				R_new = random_walk(R, j);
				// Determine whether or not to accept movement
        A = acceptance_ratio(R_new, R, j);
        if (A > rand_double(0, 1)) {
          R = R_new;
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



void Monte_Carlo::write_E_to_file(std::ofstream& outfile) {
  outfile << E_cycles[0];
  for (int i = 1; i < MC_cycles; i++) {
    outfile << "\n" << E_cycles[i];
  }
}

void Monte_Carlo::write_val_to_file(std::ofstream& outfile) {
  outfile << "alpha:  " << PDF->get_alpha() << "\n";
  outfile << "beta:   " << PDF->get_beta() << "\n";
  outfile << "<E>:    " << E << "\n";
  outfile << "accept: " << accepted_moves_ratio;
}






