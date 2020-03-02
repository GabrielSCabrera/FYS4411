#include <cmath>

#include "../matpak/Mat.h"
#include "../backend/Psi.h"
#include "../backend/random.h"

double* monte_carlo_1b(double step_size, int steps, int cycles, int N,
                       double x_max, double alpha, double beta, double a,
                       double omega, double omega_z, int equi_steps) {
  // Initialize Wave Function
  Psi PDF(alpha, beta, a, omega, omega_z);
  // Initialize Secondary Random Particle Array
  Mat P_new(N,3);
  // New Probability
  double Psi_new = 0;
  // New Energy
  double E_new;
  // Total Energy
  double E_tot = 0;
  // Total Squared Energy
  double E_squared = 0;
  // Averaging Factor
  int accum_cycles = 0;
  // Acceptance Ratio
  double W;
  // Loading Percentage
  double perc = 0;
  // Loading New Percentage
  double perc_new;
  // Old Probability
  double Psi_old = 0;
  // Initialize Random Particle Array
  Mat P = random_particles(N, -x_max, x_max);
  // Monte-Carlo Cycles
  for (int i = 0; i < cycles; i++) {
    // Initialize Random Particle Array
    Mat P = random_particles(N, -x_max, x_max);
    // Old Probability
    Psi_old = PDF(P);
    // Equilibriation
    for (int j = 0; j < equi_steps; j++) {
      // Generate New Movement
      P_new = random_walk(PDF, P, step_size);
      // Get New Probability
      Psi_new = PDF(P_new);
      // Calculate Ratio of Probabilities
      W = Psi_new/Psi_old;
      W *= W;
      // Determine whether or not to accept movement
      if (W > rand_double(0, 1)) {
        P = P_new;
      }
    }

    // Samples per Cycle
    for (int j = 0; j < steps; j++) {
      // Generate New Movement
      P_new = random_walk(PDF, P, step_size);
      // Get New Probability
      Psi_new = PDF(P_new);
      // Calculate Ratio of Probabilities
      W = Psi_new/Psi_old;
      W *= W;

      // Determine whether or not to accept movement
      if (W > rand_double(0, 1)) {
        P = P_new;
        Psi_old = Psi_new;
      }
    }

    // Set New Energy
    E_new = PDF.energy(P_new);

    E_tot += E_new;
    E_squared += E_new*E_new;
    accum_cycles += 1;
    perc_new = std::trunc(100*i/cycles);
    if (perc_new > perc) {
      perc = perc_new;
      if (perc < 10) {
        std::cout << "\r   " << perc << "%" << std::flush;
      } else if (10 <= perc && perc < 100) {
        std::cout << "\r  " << perc << "%" << std::flush;
      } else {
        std::cout << "\r " << perc << "%" << std::flush;
      }
    }
  }

  // Preparing Function Output
  double* output = new double[3];
  output[0] = E_tot/accum_cycles;
  output[1] = E_squared/accum_cycles;
  output[2] = Psi_old;

  return output;
}

void run_1b() {

  double step_size = 1E-5;    // Step size during random walk
  int steps = 1E4;            // Number of Monte-Carle steps per cycle
  int cycles = 1E3;           // Number of Monte-Carlo cycles
  int N = 1;                  // Number of Particles
  int x_max = 1;              // Maximum Initial Distance From Origin
  double a = 1E-8;            // Atomic Radius
  double omega = 1;           // Harmonic Oscillator Frequency
  double omega_z = 1;         // Harmonic Oscillator Z-Frequency
  int equi_steps = 1E3;       // Number of Steps Dedicated to Equilibriation

  // Final Energies
  double* energies = new double[5];
  double* energies_squared = new double[5];
  double* probabilities = new double[5];
  double alpha; double beta;

  for (int i = 0; i < 5; i++) {
    alpha = 0.2;
    beta = 1.5;

    // Function Outputs
    double* output = new double [3];
    output = monte_carlo_1b(step_size, steps, cycles, N, x_max, alpha,
                            beta, a, omega, omega_z, equi_steps);
    delete[] output;
  }
  delete[] energies;
}
