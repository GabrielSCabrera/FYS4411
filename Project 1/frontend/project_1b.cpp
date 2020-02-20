#include <cmath>

#include "../matpak/Mat.h"
#include "../backend/Psi.h"
#include "../backend/random.h"



double monte_carlo_1b(double step_size, int steps, int cycles, int N,
                      double x_max, double alpha, double beta, double a,
                      double omega, double omega_z, double mass) {
  // Initialize Wave Function
  Psi PDF(alpha, beta, a, omega, omega_z, mass);
  // Initialize Random Particle Array
  Mat P = random_particles(N, x_max);
  // Initialize Secondary Random Particle Array
  Mat P_new(N,3);
  // Old Probability
  double Psi_old = PDF(P);
  // New Probability
  double Psi_new;
  // Energy
  double E = PDF.energy(P);
  // Acceptance Ratio
  double W;
  // Monte-Carlo Cycles
  for (int i = 0; i < cycles; i++) {
    // Samples per Cycle
    for (int j = 0; j < steps; j++) {

      P_new = random_walk(P, step_size);
      Psi_new = PDF(P_new);

      W = std::pow(std::abs(Psi_new/Psi_old), 2);
      if (W >= 0.5) {
        P = P_new;
      }

    }
    E = PDF.energy(P);
    std::cout << E << std::endl;
  }
  return E;
}
