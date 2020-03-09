#include <cmath>

#include "../matpak/Mat.h"
#include "Psi.h"
#include "../wavefunctions/Psi.h"
#include "../wavefunctions/Psi_T.h"
#include "../wavefunctions/Psi_OB.h"

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

    2.5066282746310002 is sqrt(2*pi)
  */
  double randval = (rand_double(-sigma*4, sigma*4)- mu)/sigma;
  return N_max * std::exp(-0.5*randval*randval / (sigma*2.5066282746310002));
}

Mat random_particles(int N, double min, double max) {
  /*
    Returns a <Mat> instance of shape (N,3) full of randomly generated values
    in the range [-x_max, x_max].
    Note: Loop is parallelized.
  */
  Mat P(N,3);
  for (int i = 0; i < P.length(); i++) {
      P.set_raw(rand_double(min, max), i);
  }
  return P;
}

Mat random_walk(Psi PDF, Mat P, double step_size, long idx, double Ddt) {
  /*
    Given a <Mat> instance 'P' of shape (N, 3), selects a random index in
    the range (0, N-1) and moves it by a random direction with a magnitude
    given by the 'step_size' parameter.

    Note: New <Mat> instance is created
  */
  Mat P_new(P);

  double* force = PDF.drift(P.get(idx, 0), P.get(idx, 1), P.get(idx, 2));

  double x = rand_double(-1, 1)*step_size;
  double y = rand_double(-1, 1)*step_size;
  double z = rand_double(-1, 1)*step_size;

  P_new.set(P_new.get(idx,0) + x + force[0]*Ddt, idx, 0);
  P_new.set(P_new.get(idx,1) + y + force[1]*Ddt, idx, 1);
  P_new.set(P_new.get(idx,2) + z + force[2]*Ddt, idx, 2);

  delete [] force;

  return P_new;
}
