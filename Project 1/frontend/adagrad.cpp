#include "adagrad.h"
#include "../backend/monte_carlo.h"
#include "../matpak/Mat.h"

double* adagrad() {

  // CONSTANT TERMS

  int steps = 1E4;            // Number of Monte-Carle steps per cycle
  int cycles = 1E4;           // Number of Monte-Carlo cycles
  int N = 25;                // Number of Particles
  int x_max = 1;              // Maximum Initial Distance From Origin
  double a = 1E-8;            // Atomic Radius
  double omega = 1;           // Harmonic Oscillator Frequency
  double omega_z = 1;         // Harmonic Oscillator Z-Frequency
  int equi_steps = 1E3;       // Number of Steps Dedicated to Equilibriation
  double dt = 1E-3;           // Time Step
  double D = 0.5;             // Diffusion Constant

  // Final Energies
  double* energies = new double[5];
  double* variances = new double[5];
  double* probabilities = new double[5];
  double* acceptance_rate = new double[5];
  double alpha; double beta;

  for (int i = 0; i < 1; i++) {
    alpha = 0.5;
    beta = 1;

    // Function Outputs
    double* output = monte_carlo(steps, cycles, N, x_max, alpha, beta, a, omega,
                         omega_z, equi_steps, dt, D);
    energies[i] = output[0];
    variances[i] = output[1];
    probabilities[i] = output[2];
    acceptance_rate[i] = output[3];
    std::cout << "E: " << output[0] << ", var: " << output[1] << ", P: "
              << output[2] << ", Acc. Rate: " << output[3] << std::endl;
    delete[] output;
  }
  delete[] energies;
  delete[] variances;
  delete[] probabilities;
  delete[] acceptance_rate;
  return 0;
}
