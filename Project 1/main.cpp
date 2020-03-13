#include <iostream>
#include <chrono>

// #include "./tests/tests_backend.h"
#include "./wavefunctions/Psi.h"
#include "./wavefunctions/Psi_T.h"
#include "./wavefunctions/Psi_OB.h"
#include "./variational/monte_carlo.h"
#include "./variational/metropolis.h"
#include "./variational/hastings.h"
#include "./variational/gradient_descent.h"


/*write single array (double) to a txt-file
void doubleArrayToFile(double *v , int n, std::string filename, bool zeroPadding = false) {
  std::ofstream myfile(filename + ".txt");
  if (myfile.is_open()) {
    if (zeroPadding) {
      myfile << n+2 << "\n";
      myfile << 0.0 << "\n";
    } else {
      myfile << n << "\n";
    }
    for (int i = 0; i < n; i++) {
      myfile << v[i] << "\n";
    }
    if (zeroPadding) {
      myfile << 0.0 << "\n";
    }
  }
}*/



void run_Metropolis(bool correlated, int N, int dim) {
  //double learning_rate = 5E-4;
  int cycles = 1;
  if (correlated) {
   //Psi_OB boson_system;
    Psi_T boson_system;
    Metropolis MC(&boson_system, N, dim);

    Mat R = MC.get_initial_R_no_overlap();

    //R = MC.equilibriation(R, 100);
    //gradient_descent(&MC, learning_rate);


    R = MC.sample_energy(R, cycles);
    MC.print_info();
    R = MC.sample_energy(R, 10*cycles);
    MC.print_info();
    R = MC.sample_energy(R, 100*cycles);
    MC.print_info();
    R = MC.sample_energy(R, 1000*cycles);
    MC.print_info();


  } else {

    //Psi_OB boson_system;
    Psi_T boson_system;
    Hastings MC(&boson_system, N, dim);

    Mat R = MC.get_initial_R_no_overlap();

    //R = MC.equilibriation(R, 100);
    //gradient_descent(&MC, learning_rate);

    R = MC.sample_energy(R, cycles);
    MC.print_info();
    R = MC.sample_energy(R, 10*cycles);
    MC.print_info();
    R = MC.sample_energy(R, 100*cycles);
    MC.print_info();
    R = MC.sample_energy(R, 1000*cycles);
    MC.print_info();

    //gradient_descent(&MC, learning_rate);
  }
  printf("\n");
}

int main() {
  int N = 25;
  srand(1337);
  printf("\nHastings\n");
  run_Metropolis(false, N, 3);
  run_Metropolis(false, N, 2);
  run_Metropolis(false, N, 1);
  srand(1337);
  printf("\nMetropolis\n");
  run_Metropolis(true, N, 3);
  run_Metropolis(true, N, 2);
  run_Metropolis(true, N, 1);
}
