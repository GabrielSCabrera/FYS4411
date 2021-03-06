#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include <fstream>
#include <random>

// initializing random number generators
std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_real_distribution<double> UniformNumberGenerator(0.0, 1.0);
std::normal_distribution<double> Normaldistribution(0.0, 1.0);

// constructor
Monte_Carlo::Monte_Carlo(Psi* bosonic_system, int N_particles, int dimensions) {
	N = N_particles;
	dim = dimensions;
	bose = bosonic_system;
  E_cycles = new double [1];
  rho = new int*[3];
  rho[0] = new int[1];
  rho[1] = new int[1];
  rho[2] = new int[1];
  L = 0.5*N;
}

// destructor
Monte_Carlo::~Monte_Carlo() {
  delete[] E_cycles;
  delete[] rho[0];
  delete[] rho[1];
  delete[] rho[2];
  delete[] rho;
}

// getters
double Monte_Carlo::get_energy() {return E;}
double Monte_Carlo::get_energy_mean() {return E/(N*dim);}
double Monte_Carlo::get_grad_alpha() {return E_alpha;}
double Monte_Carlo::get_variance() {return EE - E*E;}
double Monte_Carlo::get_accepted_moves_ratio() {return accepted_moves_ratio;}
double* Monte_Carlo::get_E_cycles() {return E_cycles;}
int** Monte_Carlo::get_rho(int *length_rho) {
  (*length_rho) = N_rho;
  return rho;
}

// info (for debugging)
void Monte_Carlo::print_info() {
  printf("E: %.6lf +/- %8.2e , accept: %.4lf\n", E/(N*dim), std::sqrt((EE - E*E)/(N*dim)), accepted_moves_ratio);
}

// wrapper method for imported normal_distribution object
double Monte_Carlo::random_normal_distribution() {
  return Normaldistribution(gen);
}

// resetter method for energy and squared energy values
void Monte_Carlo::set_to_zero() {
	E = 0.0;
  EE = 0.0;
}

// generates random particle positions
Mat Monte_Carlo::get_initial_R() {
  Mat R(N, dim);
	for (int i = 0; i < N*dim; i++) {
    R.set_raw(UniformNumberGenerator(gen)*L, i);
	}
	// confirms that the particles do not overlap each other if we are looking
	// at the non-interacting case
  if (bose->interaction()) {
    double r_ij, dx;
    bool overlap;
    double hard_radius = bose->get_a();
    hard_radius *= hard_radius;
    for (int i = 0; i < N; i++) {
      overlap = true;
      while (overlap) {
        overlap = false;
        for (int j = 0; j < i; j++) {
          r_ij = 0.0;
          for (int l = 0; l < dim; l++) {
            dx = R.get(j, l) - R.get(i, l);
            r_ij += dx*dx;
          }
          if (r_ij <= hard_radius) {
            overlap = true;
            for (int l = 0; l < dim; l++) {
              R.set(UniformNumberGenerator(gen)*L, i, l);
            }
            continue;
          }
        }
      }
    }
  }
  return R;
}

// to prevent issues with pointers and overwriting of array values
void Monte_Carlo::copy_step(Mat* from, Mat* to, int index) {
  for (int i = 0; i < dim; i++) {
    to->set(from->get(index, i), index, i);
  }
}

// implements systemwise equilibriation to prepare for the true MC-cycles
Mat Monte_Carlo::equilibriation(Mat R, int cycles) {
	Mat R_new = R;
	for (int i = 0; i < cycles; i++) {
		for (int j = 0; j < N; j++) {
			random_walk(&R_new, j);
			if (acceptance_ratio(&R_new, &R, j) > UniformNumberGenerator(gen)) {
				copy_step(&R_new, &R, j);
      } else {
        copy_step(&R, &R_new, j);
			}
		}
	}
  return R;
}

// the true MC cycle from which energies are returned (then later written to file)
Mat Monte_Carlo::sample_energy(Mat R, int cycles) {
  set_to_zero();
  Mat R_new = R;
  delete[] E_cycles;
  E_cycles = new double[cycles];         // Cycle-Wise Energy
  MC_cycles = cycles;
  int accepted_moves = 0;
  double E_L;
  // Monte-Carlo Cycles
  for (int i = 0; i < cycles; i++) {
      for (int j = 0; j < N; j++) {
				random_walk(&R_new, j);
				// Determine whether or not to accept movement
        if (acceptance_ratio(&R_new, &R, j) > UniformNumberGenerator(gen)) {  // !!!!! FOR DEBUGGING
          copy_step(&R_new, &R, j);
          accepted_moves++;
        } else {
          copy_step(&R, &R_new, j);
        }
    } // End cycle
    E_L = bose->energy(&R);
    E += E_L;
    EE += E_L*E_L;
    E_cycles[i] = E_L;
  } // End Monte Carlo
  E /= cycles;
  EE /= cycles;
  accepted_moves_ratio = (double) accepted_moves/(cycles*N);
  return R;
}

// Calculate the particle density using Monte Carlo
Mat Monte_Carlo::sample_variational_derivatives(Mat R, int cycles) {
  set_to_zero();
  Mat R_new = R;
  double psi_alpha = 0.0;			      //  derivative of Psi with respect to alpha
  double psi_2alpha = 0.0;
  double E_x_psi_alpha = 0.0;		    // (derivative of Psi with respect to alpha)*E
  double E_x_psi_2alpha = 0.0;
  double psi_alpha_cycle, E_L;
	// Monte-Carlo Cycles
  for (int i = 0; i < cycles; i++) {
    for (int j = 0; j < N; j++) {
    	random_walk(&R_new, j);
      if (acceptance_ratio(&R_new, &R, j) > UniformNumberGenerator(gen)){
        copy_step(&R_new, &R, j);
      } else {
        copy_step(&R, &R_new, j);
      }
    }
    E_L = bose->energy(&R);
    E += E_L;

    psi_alpha_cycle = bose->grad_alpha(&R);
    psi_alpha  += psi_alpha_cycle;
    psi_2alpha += psi_alpha_cycle*psi_alpha_cycle;

    E_x_psi_alpha  += psi_alpha_cycle*E_L;
    E_x_psi_2alpha += psi_alpha_cycle*psi_alpha_cycle*E_L;
  } // End Monte Carlo
  E /= cycles;
  psi_alpha /= cycles;
  psi_2alpha /= cycles;
  E_x_psi_alpha /= cycles;
  E_x_psi_2alpha /= cycles;
  // E' (with respect to alpha)
  E_alpha = 2*(E_x_psi_alpha - psi_alpha*E);
  // E'' (with respect to alpha)
  E_2alpha  = 4*(E_x_psi_2alpha - E_x_psi_alpha*psi_alpha);
  E_2alpha -= 2*(psi_alpha*E_alpha) + 4*(psi_2alpha - psi_alpha*psi_alpha);
  return R;
}

/* Calculate the particle density using Monte Carlo.
        FOR THREE DIM ONLY!!
*/
Mat Monte_Carlo::one_body_density(Mat R, int cycles, int anticipated_max) {
  Mat R_new = R;
  delete[] rho[0];
  delete[] rho[1];
  delete[] rho[2];
  delete[] rho;
  int fac = 1e2;
  N_rho = 2*anticipated_max*fac + 1;
  rho = new int* [3];
  rho[0] = new int [N_rho]; // x-direction
  rho[1] = new int [N_rho]; // y-direction
  rho[2] = new int [N_rho]; // z-direction
	for (int i = 0; i < N_rho; i++) {
		rho[0][i] = 0;
    rho[1][i] = 0;
    rho[2][i] = 0;
	}
  int index;
  for (int i = 0; i < cycles; i++) {
    for (int j = 0; j < N; j++) {
      random_walk(&R_new, j);
      if (acceptance_ratio(&R_new, &R, j) > UniformNumberGenerator(gen)) {
        copy_step(&R_new, &R, j);
      } else {
        copy_step(&R, &R_new, j);
      }
      index = (R.get(j, 0) + anticipated_max)*fac + 0.5;
      rho[0][index]++;
      index = (R.get(j, 1) + anticipated_max)*fac + 0.5;
      rho[1][index]++;
      index = (R.get(j, 2) + anticipated_max)*fac + 0.5;
      rho[2][index]++;
    }
  }
  return R;
}
