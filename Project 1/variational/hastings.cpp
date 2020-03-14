#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include "hastings.h"


Mat Hastings::random_walk(Mat R, int index) {
	double r_new;
	double* F = PDF->drift_force(R, index);
	for (int k = 0; k < dim; k++) {
      r_new = R.get(index, k) + D*dt*F[k] + rand_double(-dt_sqrt , dt_sqrt);
      if (r_new < 0.0) {
        r_new += L; 
      } else if (r_new > L) {
        r_new -= L;
      }
      if (r_new < 0.0 || r_new > L) {
        printf("Outside bounds... %lf\n", F[k]);
      }
      R.set(r_new, index, k);
    }
    delete[] F;
    return R;
}


double Hastings::acceptance_ratio(Mat R_new, Mat R_old, int index) {
	double P_ratio = PDF->probability_density_ratio(R_new, R_old, index);
	double G_ratio = PDF->greens_ratio(R_old, R_new, dt, index);
  return G_ratio*P_ratio;
}

std::string Hastings::filename_E() {
  std::string path = "results/hastings/dim";
  path.append(std::to_string(dim));
  path.append("/E_");
  path.append(std::to_string(N));
  path.append(".dat");
  return path;
}

std::string Hastings::filename_val() {
  std::string path = "results/hastings/dim";
  path.append(std::to_string(dim));
  path.append("/val_");
  path.append(std::to_string(N));
  path.append(".dat");
  return path;
}
