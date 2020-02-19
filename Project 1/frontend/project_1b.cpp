#include "../matpak/Mat.h"
#include "../backend/Psi.h"
#include "../backend/random.h"

double step_size = 1;     // Step size during random walk
int samples = 1E2;        // Number of Monte-Carle samples per cycle
int cycles = 1E1;         // Number of Monte-Carlo cycles

double monte_carlo_1b() {
  int j; double alpha = 1; double beta = 0; double a = 0;
  // Initialize Wave Function
  Psi PDF(alpha, beta, a);
  // Initialize Random Particle Array
  Mat P = random_particles(5, 10);
  // Array of Energies
  // double energy = PDF.energy();
  // Old Probability
  double Psi_old = PDF(P);
  // New Probabilities
  double Psi_new;
  // Monte-Carlo Cycles
  for (int i = 0; i < cycles; i++) {
    // Samples per Cycle
    #pragma omp parallel for
    for (j = 0; j < samples; j++) {
      // Random Walk
      P = random_walk(P, step_size);

    }
  }
  return 0.0;
}

void project_1b_main() {

}
