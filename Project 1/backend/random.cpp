#include <cmath>

#include "../matpak/Mat.h"
#include "Psi.h"

double rand_double(double min, double max) {
  /*
    Returns a randomly generated <double> in range [-N_max, N_max]
  */
  return (max-min)*(static_cast<double>(rand()) / RAND_MAX) + min;
}

double normal_distribution(double N_max, double mu = 1, double sigma = 0.5) {
  /*
    Returns a randomly generated <double> in range [0, N_max] according to a
    normal distribution.

    A random number in the range [-4*sigma, 4*sigma] is generated, then is
    subsequently passed through a Gaussian function with parameters 'mu' and
    'sigma' given as arguments.
  */
  return N_max * std::exp(-0.5*std::pow((rand_double(-sigma*4, sigma*4) - mu)/sigma, 2)) /
                         (sigma*std::sqrt(2*3.141592653589793238));
}

Mat random_particles(int N, double min, double max) {
  /*
    Returns a <Mat> instance of shape (N,3) full of randomly generated values
    in the range [-x_max, x_max].
    Note: Loop is parallelized.
  */
  Mat P(N,3);
  #pragma omp parallel for
  for (int i = 0; i < P.length(); i++) {
      P.set_raw(rand_double(min, max), i);
  }
  return P;
}

Mat random_walk(Psi PDF, Mat P, double step_size) {
  /*
    Given a <Mat> instance 'P' of shape (N, 3), selects a random index in
    the range (0, N-1) and moves it by a random direction with a magnitude
    given by the 'step_size' parameter.

    Note: New <Mat> instance is created

  */
  Mat P_new(P);

  int N = P_new.shape().get(0);
  int idx = rand() % N;

  double x = rand_double(-1, 1);
  double y = rand_double(-1, 1);
  double z = rand_double(-1, 1);

  double mag = std::sqrt(x*x + y*y + z*z) / step_size;
  double* force = PDF.drift(x, y, z);

  P_new.set(P_new.get(idx,0) + x/mag + force[0], idx, 0);
  P_new.set(P_new.get(idx,1) + y/mag + force[1], idx, 1);
  P_new.set(P_new.get(idx,2) + z/mag + force[2], idx, 2);

  delete [] force;
  return P_new;
}
