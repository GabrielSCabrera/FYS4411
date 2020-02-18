#include <cmath>

#include "../backend/Psi.h"
#include "../matpak/Mat.h"

double step = 1;

double random(double N_max) {
  return N_max*(2*static_cast<double>(rand()) / RAND_MAX - 1);
}

double normal_distribution(double N_max, double mu = 1, double sigma = 0.5) {
  return N_max * std::exp(-0.5*std::pow((random(sigma*5) - mu)/sigma, 2)) /
                         (sigma*std::sqrt(2*3.141592653589793238));
}

Mat random_particles(int N, double x_max) {
  Mat P(N,3);
  for (int i = 0; i < P.length(); i++) {
      P.set_raw(random(x_max), i);
  }
  return P;
}

void project_1b_main() {
  random_particles(5, 10);
}
