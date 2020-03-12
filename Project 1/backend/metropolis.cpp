#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "../wavefunctions/Psi_T.h"
#include "../wavefunctions/Psi_OB.h"
#include "monte_carlo_class.h"
#include "metropolis.h"


Mat Metropolis::random_walk(Mat R, int index) {
  double r;
	for (int k = 0; k < dim; k++) {
      r = rand_double(-step_length, step_length);
      R.set(R.get(index, k) + r, index, k);
  }
  return R;
}

double Metropolis::acceptance_ratio(double psi_new, double psi_old, Mat R_new, Mat R_old, int index) {
	double ratio = psi_new/psi_old;
  return ratio*ratio;
}
