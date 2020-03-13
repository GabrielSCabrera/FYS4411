#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "../wavefunctions/Psi_T.h"
#include "../wavefunctions/Psi_OB.h"
#include "monte_carlo_class.h"
#include "Importance_Sampling.h"


Mat Importance_Sampling::random_walk(Mat R, int index) {
	double r;
	double* F = PDF->drift_force(R, index);
	for (int k = 0; k < dim; k++) {
      r = dt_sqrt*rand_double(-step_length, step_length);
      R.set(R.get(index, k) + D*dt*F[k] + r, index, k);
    }
    delete[] F;
    return R;
}


double Importance_Sampling::acceptance_ratio(double Psi_new, double Psi_old, Mat R_new, Mat R_old, int index) {
	double W = Psi_new/Psi_old;
	double G_ratio = PDF->greens_ratio(R_old, R_new, dt, index);
  return G_ratio*W*W;
}