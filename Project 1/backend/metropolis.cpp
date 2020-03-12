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

/*
Mat Importance_Sampling::random_walk(Mat R) {
	for (int k = 0; k < dim; k++) {
      r = rand_double(-step_size, step_size);
      R.set(R.get(j, k) + r, j, k);
    }
    return R;
}

double Importance_Sampling::acceptance_ratio(double Psi_new, double Psi_old, Mat R_new, Mat R_old, int index) {
	double W = Psi_new/Psi_old;
	W *= W;
	double G_ratio = PDF->greens_ratio(R_old, R_new, G_K, index);
    return G_ratio*W*W;
}*/
