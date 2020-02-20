#include <cmath>

#include "../matpak/Mat.h"

double rand_double(double N_max) {
  /*
    Returns a randomly generated <double> in range [-N_max, N_max]
  */
  return N_max*(2*static_cast<double>(rand()) / RAND_MAX - 1);
}

double normal_distribution(double N_max, double mu = 1, double sigma = 0.5) {
  /*
    Returns a randomly generated <double> in range [0, N_max] according to a
    normal distribution.

    A random number in the range [-4*sigma, 4*sigma] is generated, then is
    subsequently passed through a Gaussian function with parameters 'mu' and
    'sigma' given as arguments.
  */
  return N_max * std::exp(-0.5*std::pow((rand_double(sigma*4) - mu)/sigma, 2)) /
                         (sigma*std::sqrt(2*3.141592653589793238));
}

Mat random_particles(int N, double x_max) {
  /*
    Returns a <Mat> instance of shape (N,3) full of randomly generated values
    in the range [-x_max, x_max].
    Note: Loop is parallelized.
  */
  Mat P(N,3);
  #pragma omp parallel for
  for (int i = 0; i < P.length(); i++) {
      P.set_raw(rand_double(x_max), i);
  }
  return P;
}

Mat random_walk(Mat P, double step_size) {
  /*
    Given a <Mat> instance 'P' of shape (N, 3), selects a random index in
    the range (0, N-1) and moves it by a random direction with a magnitude
    given by the 'step_size' parameter.

    Note: New <Mat> instance is created

  */
  Mat P_new(P);

  int N = P_new.shape().get(0);
  int idx = rand() % N;

  double x = rand_double(1);
  double y = rand_double(1);
  double z = rand_double(1);

  double mag = std::sqrt(x*x + y*y + z*z) / step_size;

  P_new.set(P_new.get(idx,0) + x/mag, idx, 0);
  P_new.set(P_new.get(idx,1) + y/mag, idx, 1);
  P_new.set(P_new.get(idx,2) + z/mag, idx, 2);

  return P_new;
}
