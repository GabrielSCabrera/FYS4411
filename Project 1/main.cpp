#include <iostream>
#include <chrono>
#include <fstream>
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
/*
  gradient_descent(MC, 1E-5, R);
  printf("alpha: %.6lf, beta: %.6lf\n", MC->PDF->get_alpha(), MC->PDF->get_beta());

  cycles = 1E3;
  R = MC->sample_energy(R, cycles);
  MC->print_info();

  MC->PDF->update_alpha(0.5);
  R = MC->equilibriation(R, equi_cycles);
  R = MC->sample_energy(R, cycles);
  MC->print_info();


  int N = R.shape0();
  ofstream myfile;
  string filename = "results/";
  filename.append("test1.dat");
  myfile.open (filename);
  myfile << "Writing this to a file.\n";
  myfile.close();*/
}



void run_Metropolis(bool correlated, int N, int dim) {
  //double learning_rate = 5E-4;
  if (correlated) {
   //Psi_OB boson_system;
    Psi_T boson_system;
    Metropolis MC(&boson_system, N, dim);
    run(&MC);
  } else {
    //Psi_OB boson_system;
    Psi_T boson_system;
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
