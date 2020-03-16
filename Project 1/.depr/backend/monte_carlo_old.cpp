#include <cmath>
#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "../wavefunctions/Psi_T.h"
#include "../wavefunctions/Psi_OB.h"
#include "random.h"
#include "monte_carlo.h"

/*
This is now the metropolis algorithm. 
All the original code is still in importance_sampling

*/

double* monte_carlo(Psi* PDF, int N, int dim, double x_max, int cycles, int equi_steps) {
  // Initialize Secondary Random Particle Array
  Mat P = random_particles(N, -x_max, x_max);
  Mat P_new(N, dim);
  double Psi_new, Psi_old;
  double* E_cycle = new double[cycles];         // Cycle-Wise Energy
  double W;                                     // Acceptance Ratio
  double* output = new double[6] {0,0,0,0,0,0}; // Preparing Function Output
  double step_size = 0.1;

  // Equilibriation
  double r;     // random number
  Psi_old = PDF->operator()(P);
  for (int i = 0; i < equi_steps; i++) {
    for (int j = 0; j < N; j++) {
        for (int k = 0; k < dim; k++) {
          r = rand_double(-step_size, step_size);
          P_new.set(P.get(j, k) + r, j, k);
        }
      // Get New Probability
      Psi_new = PDF->operator()(P_new);
      // Calculate Ratio of Probabilities
      W = Psi_new/Psi_old;
      W *= W;
      // Determine whether or not to accept movement
      if (W > rand_double(0, 1)) {
        P = P_new;
        Psi_old = Psi_new;
      }
    } // End Equilibriation
  }


  double Psi_alpha = 0.0;
  double E_Psi_alpha = 0.0;
  double Psi_beta = 0.0;
  double E_Psi_beta = 0.0;
  double temp;

  // #pragma omp parallel for
  // Monte-Carlo Cycles
  for (int i = 0; i < cycles; i++) {
      // Samples per Cycle
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < dim; k++) {
          r = rand_double(-step_size, step_size);
          P_new.set(P.get(j, k) + r, j, k);
        }
        // Get New Probability
        Psi_new = PDF->operator()(P_new);
        // Calculate Ratio of Probabilities
        W = Psi_new/Psi_old;
        W *= W;

        // Determine whether or not to accept movement
        if (W > rand_double(0, 1)) {
          P = P_new;
          Psi_old = Psi_new;
          output[3] += 1;
        }
    } // End cycle

    E_cycle[i] = PDF->energy(P);
    temp = PDF->grad_alpha(P);
    Psi_alpha += temp;
    E_Psi_alpha += temp*E_cycle[i];

    temp = PDF->grad_beta(P);
    Psi_beta += temp;
    E_Psi_beta += temp*E_cycle[i];
   
    // Save Total Energy
    output[0] += E_cycle[i];

  } // End Monte Carlo

  output[0] /= cycles;

  // Calculating Variance Iteratively
  for (int i = 0; i < cycles; i++) {
    output[1] += (E_cycle[i]*E_cycle[i] - output[0]*output[0]) *
                 (E_cycle[i] - output[0]) * (E_cycle[i] - output[0]);
  }
  E_Psi_alpha /= cycles;
  Psi_alpha /= cycles;
  E_Psi_beta /= cycles;
  Psi_beta /= cycles;

  output[2] = Psi_old;
  output[3] /= cycles*N;
  output[4] = 2*(E_Psi_alpha - Psi_alpha*output[0]);
  output[5] = 2*(E_Psi_beta  - Psi_beta*output[0]);
  // Cleanup
  delete [] E_cycle;

  return output;
}
