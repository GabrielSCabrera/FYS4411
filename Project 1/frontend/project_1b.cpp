#include <cmath>

#include "../matpak/Mat.h"
#include "../backend/Psi.h"
#include "../backend/random.h"

double monte_carlo_1b(double step_size, int steps, int cycles, int N,
                      double x_max, double alpha, double beta, double a,
                      double omega, double omega_z, int equi_steps) {
  // Initialize Wave Function
  Psi PDF(alpha, beta, a, omega, omega_z);
  // Initialize Secondary Random Particle Array
  Mat P_new(N,3);
  // New Probability
  double Psi_new;
  // New Energy
  double E_new;
  // Total Energy
  double E_tot = 0;
  // Averaging Factor
  int accum_cycles = 0;
  // Acceptance Ratio
  double W;
  // Random Number between 0 and 1
  double r;
  // Loading Percentage
  double perc = 0;
  // Loading New Percentage
  double New_perc = 0;
  // Loading Bar Output
  double perc_disp = 0;
  // Monte-Carlo Cycles
  for (int i = 0; i < cycles; i++) {
    // Initialize Random Particle Array
    Mat P = random_particles(N, -x_max, x_max);
    // Old Probability
    double Psi_old = PDF(P);
    // Old Energy
    double E_old = PDF.energy(P);
    // Samples per Cycle
    for (int j = 0; j < steps; j++) {

      P_new = random_walk(PDF, P, step_size);
      Psi_new = PDF(P_new);

      W = std::pow(std::abs(Psi_new/Psi_old), 2);
      r = rand_double(0, 1);
      E_new = PDF.energy(P_new);

      if ((j < equi_steps && W > r) || (j >= equi_steps && (E_new < E_old || W > r))) {
        P = P_new;
        E_old = E_new;
        Psi_old = Psi_new;
      }
    }
    E_tot += E_old;
    accum_cycles += 1;
    std::cout << "\r" << ""
    // std::cout << E_tot/accum_cycles << " " << std::endl;
  }
  return E_tot/accum_cycles;
}

void run_1b() {

  double step_size = 1E-5;    // Step size during random walk
  int steps = 1E5;            // Number of Monte-Carle steps per cycle
  int cycles = 1E4;           // Number of Monte-Carlo cycles
  int N = 1;                  // Number of Particles
  int x_max = 1;              // Maximum Initial Distance From Origin
  double a = 1E-8;            // Atomic Radius
  double omega = 1;           // Harmonic Oscillator Frequency
  double omega_z = 1;         // Harmonic Oscillator Z-Frequency
  int equi_steps = 1E3;       // Number of Steps Dedicated to Equilibriation

  double* alphas = new double[5] {6E-1, 5E-1, 4E-1, 3E-1};  // Hyperparameter
  double* betas = new double[5] {6E-1, 5E-1, 4E-1, 3E-1};   // Hyperparameter
  double* energies = new double[5];                         // Final Energies

  double alpha; double beta;

  for (int i = 0; i < 5; i++) {
    alpha = alphas[i];
    beta = betas[i];
    energies[i] = monte_carlo_1b(step_size, steps, cycles, N, x_max, alpha,
                                 beta, a, omega, omega_z, equi_steps);
  }
}
