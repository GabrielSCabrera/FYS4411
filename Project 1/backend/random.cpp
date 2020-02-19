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
    Iterates through a <Mat> instance 'P' of shape (N,3) and moves each vector
    (or row in 'P') by a fixed 'step_size', but in a random direction.

    Note: New <Mat> instance is created

    Note: Loop is parallelized.
  */
  Mat P_new(P);
  int N = P_new.shape().get(0);
  double x; double y; double z; double mag;
  #pragma omp parallel for
  for (int i = 0; i < N; i++) {
    x = rand_double(1);
    y = rand_double(1);
    z = rand_double(1);
    mag = std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2)) / step_size;
    P_new.set(P_new.get(i,0) + x/mag, i, 0);
    P_new.set(P_new.get(i,1) + y/mag, i, 1);
    P_new.set(P_new.get(i,2) + z/mag, i, 2);
  }

  return P_new;
}
