#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include "hastings.h"

// Biased walk (by the quantum force) of k'th particle
void Hastings::random_walk(Mat* R, int k) {
	double r_new;
	double* F = bose->drift_force(R, k);
	for (int l = 0; l < dim; l++) {
    r_new = R->get(k, l) + D*dt*F[l] + dt_sqrt*random_normal_distribution();
    R->set(r_new, k, l);
  }
  delete[] F;
}

// Acceptance ratio for new position of k'th particle
double Hastings::acceptance_ratio(Mat* R_new, Mat* R_old, int k) {
	double P_ratio = bose->probability_density_ratio(R_new, R_old, k);
	double G_ratio = bose->greens_ratio(R_new, R_old, dt, k);
  return G_ratio*P_ratio;
}

// update time step
void Hastings::set_dt(double new_dt) {
  dt = new_dt;
  dt_sqrt = std::sqrt(dt); 
}
