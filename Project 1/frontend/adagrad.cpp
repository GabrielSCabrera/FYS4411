#include <cmath>

#include "../backend/monte_carlo.h"
#include "../wavefunctions/Psi.h"
#include "../wavefunctions/Psi_OB.h"
#include "../wavefunctions/Psi_T.h"
#include "../matpak/Mat.h"
#include "adagrad.h"


double* adagrad(Psi* PDF) {
  // Monte Carlo Parameters
  int steps = 1E4;            // Number of Monte-Carle steps per cycle
  int cycles = 2E2;           // Number of Monte-Carlo cycles
  int N = 1;                  // Number of Particles
  int x_max = 1;              // Maximum Initial Distance From Origin
  int equi_steps = 1E3;       // Number of Steps Dedicated to Equilibriation
  double dt = 1E-2;           // Time Step
  double D = 0.5;             // Diffusion Constant
  double eta = 1E-4;          // Learning Rate
  // double eps = 1E-8;          // Small Number Correction
  // double gamma = 0.9;         // Momentum Term

  // Wavefunction and Potential Constants
  double a = 1E-10;           // Atomic Radius
  double gamma = 1;           // Potential Elongation Factor
  double mass = 1;

  // Gradient Descent Parameters
  int N_steps = 1E2;          // Number of AdaGrad steps
  double alpha_0 = 1;         // Second value of alpha
  double alpha_1 = 0.9;       // First  value of alpha
  double beta_0 = 1;          // Second value of beta
  double beta_1 = 1;          // First  value of beta

  // Arrays
  double* energies = new double[N_steps];
  double* alphas = new double[N_steps+2];
  double* betas = new double[N_steps+2];

  // Setting Initial Conditions
  alphas[0] = alpha_0;
  alphas[1] = alpha_1;
  betas[0] = beta_0;
  betas[1] = beta_1;

  // Initialize Wavefunctions
  PDF->update_alpha(alpha_0);
  PDF->update_beta(beta_0);
  PDF->update_a(a);
  PDF->update_gamma(gamma);
  PDF->update_mass(mass);

  // Allocation of Temporary Variables for the Learning Rate
  // in Both the Directions of Alpha and Beta
  double* learn_rates = new double[2];

  // Allocation of Temporary Variables for the Derivatives of
  // <E> with Respect to Both Alpha and Beta
  double* derivatives = new double[2];

  // Setting Initial Values
  double* output = monte_carlo(PDF, steps, cycles, N, x_max, equi_steps, dt, D);
  // learn_rates[0] = (1-gamma)*output[4]*output[4];
  // learn_rates[1] = (1-gamma)*output[5]*output[5];

  // Clearing Allocated Memory
  delete [] output;

  for (int i = 0; i < N_steps; i++) {

    // Updating Wavefunctions
    PDF->update_alpha(alphas[i+1]);
    // PDF->update_beta(betas[i+1]);

    // Running Monte-Carlo
    double* output = monte_carlo(PDF, steps, cycles, N, x_max, equi_steps, dt, D);

    // Displaying Stats
    std::cout << "MC-Cycle – alpha = " << alphas[i+1] << ", beta = " << betas[i]
              << ", E = " << output[0] << std::endl;

    // Updating Derivatives
    derivatives[0] = output[4];
    derivatives[1] = output[5];

    // Updating Learning Rates
    // learn_rates[0] = gamma*learn_rates[0] + (1-gamma)*output[4]*output[4];
    // learn_rates[1] = gamma*learn_rates[1] + (1-gamma)*output[5]*output[5];

    // Updating Parameters
    // alphas[i+2] = alphas[i+1] - eta*derivatives[0]/std::sqrt(learn_rates[0] + eps);
    // betas[i+2]  = betas[i+1]  - eta*derivatives[1]/std::sqrt(learn_rates[1] + eps);
    alphas[i+2] = alphas[i+1] - eta*derivatives[0];
    // betas[i+2]  = betas[i+1]  - eta*derivatives[1];
    betas[i+2]  = betas[i+1];

    // Clearing Allocated Memory
    delete [] output;
  }

  // Clearing Allocated Memory
  delete [] energies;
  delete [] alphas;
  delete [] betas;
  delete [] learn_rates;
  delete [] derivatives;
  return 0;
}
