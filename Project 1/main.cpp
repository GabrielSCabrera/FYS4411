#include <iostream>
#include <chrono>

// #include "./tests/tests_backend.h"
#include "./wavefunctions/Psi.h"
#include "./wavefunctions/Psi_T.h"
#include "./wavefunctions/Psi_OB.h"
//#include "./frontend/adagrad.h"
#include "./backend/monte_carlo_class.h"
#include "./backend/metropolis.h"
#include "./backend/Importance_Sampling.h"
#include "./backend/gradient_descent.h"

// void run_all_tests() {
//   tests_Psi();
// }

/*void run_all_parts() {
  bool one_body = true;  // true: interacting mode, false: one-body mode
  if (one_body) {
    Psi_OB PDF(0, 0, 0, 0);
    adagrad(&PDF);
  } else {
    Psi_T PDF(0, 0, 0, 0);
    adagrad(&PDF);
  }
}*/

void run_Metropolis(bool correlated) {
  double a = 0.0043;          // Atomic Radius
  double gamma = 1.0;         // Potential Elongation Factor
  double alpha = 0.5;
  double beta = 1.0;
  double learning_rate = 5E-3;
  int N = 10;
  int dim = 3;

  if (correlated) {

    Psi_T boson_system(alpha, beta, a, gamma);
    //Importance_Sampling MC(&boson_system, N, dim);
    Metropolis MC(&boson_system, N, dim);

    Mat R = MC.get_initial_R();
    R = MC.equilibriation(R, 100);
    R = MC.sample_energy(R, 1E3);

    MC.print_info();

    boson_system.update_alpha(0.3);
    boson_system.update_beta(0.7);

    gradient_descent(&MC, learning_rate);

  } else {

    Psi_OB boson_system(alpha, beta, a, gamma);
    Importance_Sampling MC(&boson_system, N, dim);
    //Metropolis MC(&boson_system, N, dim);

    Mat R = MC.get_initial_R();
    R = MC.equilibriation(R, 100);
    R = MC.sample_energy(R, 1E3);

    MC.print_info();

    boson_system.update_alpha(0.3);
    boson_system.update_beta(0.7);

    gradient_descent(&MC, learning_rate);
  }
}

int main() {
  srand(1337);
  // run_all_tests();
  // run_all_parts();
  run_Metropolis(false);
}
