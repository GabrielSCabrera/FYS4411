#include "adagrad.h"
#include "../backend/monte_carlo.h"
#include "../matpak/Mat.h"

double get_energy(double alpha, double beta) {

  // HARDCODED CONSTANTS

  int steps = 1E4;            // Number of Monte-Carle steps per cycle
  int cycles = 1E4;           // Number of Monte-Carlo cycles
  int N = 25;                 // Number of Particles
  int x_max = 1;              // Maximum Initial Distance From Origin
  double a = 1E-8;            // Atomic Radius
  double omega = 1;           // Harmonic Oscillator Frequency
  double omega_z = 1;         // Harmonic Oscillator Z-Frequency
  int equi_steps = 1E3;       // Number of Steps Dedicated to Equilibriation
  double dt = 1E-3;           // Time Step
  double D = 0.5;             // Diffusion Constant

  // Function Outputs
  double* output = monte_carlo(steps, cycles, N, x_max, alpha, beta, a, omega,
                       omega_z, equi_steps, dt, D);

  double energy = output[0];
  delete[] output;
  return energy;

}

double* adagrad() {

  // Gradient Descent Parameters

  int N_steps = 1E2;      // Number of AdaGrad steps

  double alpha_0 = 0.5;   // Initial value of alpha
  double alpha_1 = 0.51;  // Second value of alpha

  double beta_0 = 1;      // Initial value of beta
  double beta_1 = 1.1;    // Second value of beta

  // Arrays

  double* energies = new double[N_steps];
  double* alphas = new double[N_steps+2];
  double* betas = new double[N_steps+2];

  alphas[0] = alpha_0;
  alphas[1] = alpha_1;

  betas[0] = beta_0;
  betas[1] = beta_1;

  // Allocation of Temporary Variables

  double learn_rate;

  for (int i = 0; i < N_steps; i++) {
    
  }

  delete [] energies;
  delete [] alphas;
  delete [] betas;
  return 0;
}
