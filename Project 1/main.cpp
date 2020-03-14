#include <iostream>
#include <chrono>
#include <fstream>
#include <string>
#include <mpi.h>
using namespace std;

// #include "./tests/tests_backend.h"
#include "./wavefunctions/Psi.h"
#include "./wavefunctions/Psi_T.h"
#include "./wavefunctions/Psi_OB.h"
#include "./variational/monte_carlo.h"
#include "./variational/metropolis.h"
#include "./variational/hastings.h"
#include "./variational/gradient_descent.h"

void run(Monte_Carlo* MC) {
  int cycles = 1E5;
  //int equi_cycles = 1E2;
  Mat R = MC->get_initial_R_no_overlap();

  R = MC->sample_energy(R, 1);
  MC->print_info();
  //R = MC->equilibriation(R, equi_cycles);
  R = MC->sample_energy(R, cycles);
  MC->print_info();
  //R = MC->equilibriation(R, equi_cycles);
  R = MC->sample_energy(R, cycles);
  MC->print_info();
  //R = MC->equilibriation(R, equi_cycles);
  R = MC->sample_energy(R, cycles);
  MC->print_info();


  ofstream outfile;
  string filename = MC->filename_E();
  outfile.open(filename);
  MC->write_E_to_file(outfile);
  outfile.close();

  filename = MC->filename_val();
  outfile.open(filename);
  MC->write_val_to_file(outfile);
  outfile.close();
}



// bool correlated is a lying bastard
void run_Metropolis(bool correlated, int N, int dim, double learning_rate=1E-4) {
  if (correlated) {
    Psi_OB boson_system;
    //Psi_T boson_system;
    Metropolis MC(&boson_system, N, dim);
    run(&MC);
  } else {
    Psi_OB boson_system;
    //Psi_T boson_system;
    Hastings MC(&boson_system, N, dim);
    run(&MC);
  }
  printf("\n");
}


int main() {
  int N = 3;
  //srand(1336);
  
  printf("\nHastings\n");
  run_Metropolis(false, N, 3);
  run_Metropolis(false, N, 2);
  run_Metropolis(false, N, 1);

  printf("\nMetropolis\n");
  run_Metropolis(true, N, 3);
  run_Metropolis(true, N, 2);
  run_Metropolis(true, N, 1);
}
