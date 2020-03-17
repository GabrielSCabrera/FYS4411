#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include "metropolis.h"


void Metropolis::random_walk(Mat* R, int index) {
  double r_new;
	for (int k = 0; k < dim; k++) {
      r_new = R->get(index, k) + step_length*random_normal_distribution();
      R->set(r_new, index, k);
  }
}

double Metropolis::acceptance_ratio(Mat R_new, Mat R_old, int index) {
  return PDF->probability_density_ratio(R_new, R_old, index);
}

std::string Metropolis::filename_E() {
  std::string path = "results/";
  path.append(PDF->name());
  path.append("/metropolis/dim");
  path.append(std::to_string(dim));
  path.append("/E_");
  path.append(std::to_string(N));
  path.append(".dat");
  return path;
}

std::string Metropolis::filename_val() {
  std::string path = "results/";
  path.append(PDF->name());
  path.append("/metropolis/dim");
  path.append(std::to_string(dim));
  path.append("/val_");
  path.append(std::to_string(N));
  path.append(".dat");
  return path;
}
