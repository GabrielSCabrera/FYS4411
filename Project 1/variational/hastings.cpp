#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include "hastings.h"


Mat Hastings::random_walk(Mat R, int index) {
	double nu;
	double* F = PDF->drift_force(R, index);
	for (int k = 0; k < dim; k++) {
      nu = rand_double(-dt_sqrt , dt_sqrt);
      R.set(R.get(index, k) + D*dt*F[k] + nu, index, k);
    }
    delete[] F;
    return R;
}


double Hastings::acceptance_ratio(double Psi_new, double Psi_old, Mat R_new, Mat R_old, int index) {
	double W = Psi_new/Psi_old;
	double G_ratio = PDF->greens_ratio(R_old, R_new, dt, index);
  return G_ratio*W*W;
}
