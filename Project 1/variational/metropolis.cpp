#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include "metropolis.h"


void Metropolis::random_walk(Mat* R, int k) {
  double r_new;
	for (int l = 0; l < dim; l++) {
      r_new = R->get(k, l) + step_length*random_normal_distribution();
      R->set(r_new, k, l);
  }
}

double Metropolis::acceptance_ratio(Mat* R_new, Mat* R_old, int k) {
  return bose->probability_density_ratio(R_new, R_old, k);
}
