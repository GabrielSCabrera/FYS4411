
/*

NOT DONE

*/

#include <cmath>

#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "../wavefunctions/Psi_T.h"
#include "../wavefunctions/Psi_OB.h"
#include "random.h"
#include "monte_carlo.h"

double* monte_carlo(Psi* PDF, int cycles, int N, double x_max, int equi_steps,
                    double dt, double D) {
  // Step Size
  double step_size = std::sqrt(dt);
  // Initialize Secondary Random Particle Array
  Mat P_new(N,3);
  // New, Old Probability
  double Psi_new = 0; double Psi_old = 0;
  // Cycle-Wise Energy
  double* E_cycle = new double[cycles];
  // Averaging Factor
  int accum_cycles = 0;
  // Acceptance Ratio
  double W;
  // Loading Old, New Percentage
  double perc = 0; double perc_new;
  // Preparing Function Output
  double* output = new double[6] {0,0,0,0,0,0};
  // Initialize Random Particle Array
  Mat P = random_particles(N, -x_max, x_max);
  // Green's Function Ratio, constant
  double G_ratio; double G_K = D*dt;
  // Old Probability
  Psi_old = PDF->operator()(P);

  for (int i = 0; i < equi_steps; i++) {
    // Equilibriation
    for (int j = 0; j < N; j++) {
      P_new = random_walk(PDF, P, step_size, j, G_K);
      // Get New Probability
      Psi_new = PDF->operator()(P_new);
      // Calculate Ratio of Probabilities
      W = Psi_new/Psi_old;
      W *= W;
      // Include Drift Force
      G_ratio = PDF->greens_ratio(P.get(j, 0), P.get(j, 1),
                                  P.get(j, 2), P_new.get(j, 0),
                                  P_new.get(j, 1), P_new.get(j, 2), G_K);

      W *= G_ratio;
      // Determine whether or not to accept movement
      if (W > rand_double(0, 1)) {
        P = P_new;
      }
    } // End Equilibriation
  }

  double Psi_alpha = 0.0;
  double E_Psi_alpha = 0.0;
  double Psi_beta = 0.0;
  double E_Psi_beta = 0.0;
  double temp;
  double E_beta = 0.0;

  // #pragma omp parallel for
  // Monte-Carlo Cycles
  for (int i = 0; i < cycles; i++) {
      // Samples per Cycle
      for (int j = 0; j < N; j++) {
        P_new = random_walk(PDF, P, step_size, j, G_K);
        // Get New Probability
        Psi_new = PDF->operator()(P_new);
        // Calculate Ratio of Probabilities
        W = Psi_new/Psi_old;
        W *= W;
        // Include Drift Force
        G_ratio = PDF->greens_ratio(P.get(j, 0), P.get(j, 1), P.get(j, 2),
                                    P_new.get(j, 0), P_new.get(j, 1),
                                    P_new.get(j, 2), G_K);
        W *= G_ratio;
        // Determine whether or not to accept movement
        if (W > rand_double(0, 1)) {
          P = P_new;
          Psi_old = Psi_new;
          output[3] += 1;
        }
    } // End cycle
    E_cycle[i] = PDF->energy(P_new);
    temp = PDF->grad_alpha(P);
    Psi_alpha += temp;
    E_Psi_alpha += temp*E_cycle[i];

    temp = PDF->grad_beta(P);
    Psi_beta += temp;
    E_Psi_beta += temp*E_cycle[i];
    // Set New Energy

    // Save Total Energy
    output[0] += E_cycle[i];
    accum_cycles += 1;
    perc_new = std::trunc(100*accum_cycles/cycles);

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

  } // End Monte Carlo
  std::cout << "\r 100%" << std::endl;

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
