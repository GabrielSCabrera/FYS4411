

Mat Metropolis::random_walk(Mat R) {
	for (int k = 0; k < dim; k++) {
      r = rand_double(-step_size, step_size);
      R.set(R.get(j, k) + r, j, k);
    }
    return R;
}

double Metropolis::acceptance_ratio(double Psi_new, double Psi_old, Mat R_new, Mat R_old, int index) {
	double ratio = Psi_new/Psi_old;
    return ratio*ratio;
}



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
	/*
	G_ratio must be open for 1, 2 and 3 dim.
	Can send in more parameters if needed
	*/
	double G_ratio = PDF->greens_ratio(R_old, R_new, G_K, index);
    return G_ratio*W*W;
}