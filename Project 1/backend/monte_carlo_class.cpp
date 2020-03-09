#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "../wavefunctions/Psi_T.h"
#include "../wavefunctions/Psi_OB.h"



Monte_Carlo::Monte_Carlo(int N_particles, int dimensions, int x_limit, Psi* PDF_) {
	N = N_particles;
	dim = dimensions;
	x_max = x_limit;
	PDF = PDF_;
}

double Monte_Carlo::rand_double(double min, double max) {
  /*
    Returns a randomly generated <double> in range [-N_max, N_max]
  */
  return (max-min)*(static_cast<double>(rand()) / RAND_MAX) + min;
}


Mat Monte_Carlo::get_initial_R() {
	/*
	Returns a <Mat> instance of shape (N, dim) full of randomly generated values
	in the range [-x_max, x_max].
	*/
	Mat R(N, dim);
	for (int i = 0; i < N; i++) {
			for (int j = 0; j < dim; j++) {
				R.set_raw(rand_double(-x_max, x_max), i, j);
			}
	}
	return R;
}


Mat Monte_Carlo::equilibriation(Mat R, cycles) {
	Mat R_new(N, dim);
	Psi_old = PDF->operator()(R);
	for (int i = 0; i < cycles; i++) {
		for (int j = 0; j < N; j++) {
			R_new = random_walk();
			// Get New Probability
			Psi_new = PDF->operator()(P_new);
			// Calculate Ratio of Probabilities
			W = Psi_new/Psi_old;
			W *= W;
			// Determine whether or not to accept movement
			if (W > rand_double(0, 1)) {
				R = R_new;
			}
		}
	}
	return R;
}