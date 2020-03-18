#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include <fstream>
#include <random>

std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_real_distribution<double> UniformNumberGenerator(0.0, 1.0);
std::normal_distribution<double> Normaldistribution(0.0, 1.0);

Monte_Carlo::Monte_Carlo(Psi* bosonic_system, int N_particles, int dimensions) {
	N = N_particles;
	dim = dimensions;
	bose = bosonic_system;
  E_cycles = new double [1];
  L = 0.5*N;
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
  return E_alpha;
}

double Monte_Carlo::get_variance() {
  return EE - E*E;
}

double Monte_Carlo::get_accepted_moves_ratio() {
  return accepted_moves_ratio;
}

double* Monte_Carlo::get_E_cycles() {
  return E_cycles;
}

void Monte_Carlo::print_info() {
  printf("E: %.6lf, var: %.6lf , acceptance: %.6lf\n", E/(N*dim), EE - E*E, accepted_moves_ratio);
}

double Monte_Carlo::random_normal_distribution() {
  return Normaldistribution(gen);
}

void Monte_Carlo::set_to_zero() {
	E = 0.0;
	EE = 0.0;
  E_alpha = 0.0;
  E_2alpha = 0.0;
}


Mat Monte_Carlo::get_initial_R() {
  Mat R(N, dim);
	for (int i = 0; i < N*dim; i++) {
    R.set_raw(UniformNumberGenerator(gen)*L, i);
	}
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


void Monte_Carlo::copy_step(Mat* from, Mat* to, int index) {
  for (int i = 0; i < dim; i++) {
    to->set(from->get(index, i), index, i);
  }
}


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


Mat Monte_Carlo::sample_energy(Mat R, int cycles) {
  set_to_zero();
  Mat R_new = R;
  delete[] E_cycles;
  E_cycles = new double[cycles];         // Cycle-Wise Energy
  MC_cycles = cycles;
  int accepted_moves = 0;
  double A;															// Acceptance Ratio
  double r, E_L;
  //double psi_R;
  // Monte-Carlo Cycles
  for (int i = 0; i < cycles; i++) {
      for (int j = 0; j < N; j++) {
				random_walk(&R_new, j);

				// Determine whether or not to accept movement
        A = acceptance_ratio(&R_new, &R, j);
        r = UniformNumberGenerator(gen);
        //printf("(%5.2lf,%5.2lf)", A, r);
        if (A > r) {  // !!!!! FOR DEBUGGING
          copy_step(&R_new, &R, j);
          accepted_moves++;
          //printf("+ ");
        } else {
          copy_step(&R, &R_new, j);
          //printf("- ");
        }
    } // End cycle
    E_L = bose->energy(&R);
    //psi_R = bose->psi(&R);
    // debugging:
    //printf("\t-> %8.3lf : %6.3e \n", E_L/(N*dim), psi_R*psi_R);
    // Save Total Energy
    E += E_L;
    EE += E_L*E_L;
    E_cycles[i] = E_L;

  } // End Monte Carlo
  E /= cycles;
  EE /= cycles;
  accepted_moves_ratio = (double) accepted_moves/(cycles*N);
  return R;
}

Mat Monte_Carlo::sample_variational_derivatives(Mat R, int cycles) {
  set_to_zero();
  Mat R_new = R;
  double psi_alpha = 0.0;			      //  derivative of Psi with respect to alpha
  double E_x_psi_alpha = 0.0;		    // (derivative of Psi with respect to alpha)*E
  double psi_2alpha = 0.0;          // double derivative of Psi with respect to alpha
  double E_x_psi_2alpha = 0.0;
  double psi_alpha_cycle, psi_2alpha_cycle, E_L;

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
    psi_2alpha_cycle = psi_alpha_cycle*psi_alpha_cycle;
    psi_alpha += psi_alpha_cycle;
    psi_2alpha += psi_2alpha_cycle;

    E_x_psi_alpha += psi_alpha_cycle*E_L;
    E_x_psi_2alpha += psi_2alpha_cycle*E_L;
  } // End Monte Carlo
  E /= cycles;
  psi_alpha /= cycles;
  psi_2alpha /= cycles;
  E_x_psi_alpha /= cycles;
  E_x_psi_2alpha /= cycles;

  E_alpha = 2*(E_x_psi_alpha - psi_alpha*E);
  E_2alpha = 2*(E_x_psi_2alpha - E_x_psi_alpha*psi_alpha);
  //E_2alpha = 2*(E_x_psi_2alpha - E*psi_2alpha);
  E_2alpha -= (psi_alpha*E_alpha);// + E*psi_2alpha);
  return R;
}


void Monte_Carlo::write_E_to_file(std::ofstream& outfile) {
  outfile << E_cycles[0];
  for (int i = 1; i < MC_cycles; i++) {
    outfile << "\n" << E_cycles[i];
  }
}

void Monte_Carlo::write_val_to_file(std::ofstream& outfile) {
  outfile << "alpha " << bose->get_alpha() << "\n";
  outfile << "beta " << bose->get_beta() << "\n";
  outfile << "E " << E << "\n";
	outfile << "accept " << accepted_moves_ratio << "\n";
  outfile << "cycles " << MC_cycles;
}
