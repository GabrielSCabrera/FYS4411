#include <iostream>
#include <chrono>
#include <fstream>
#include <string>
#include <mpi.h>
using namespace std;

#include "./wavefunctions/Psi.h"
#include "./wavefunctions/Psi_T.h"
#include "./wavefunctions/Psi_OB.h"
#include "./variational/monte_carlo.h"
#include "./variational/metropolis.h"
#include "./variational/hastings.h"
#include "./variational/gradient_descent.h"
#include "monte_carlo_simulation.h"

int main(int narg, char** argv) {
  int my_rank, num_procs;
  MPI_Init(&narg, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  //------CHANGE PARAMETERS-------
  int N = 5; 
  int dim = 3; 
  double learning_rate = 1E-4;
  bool importance_samplig = false;
  bool correlated = true;
  //------------------------------

  if (importance_samplig) {
    if (correlated) {
      Psi_T boson_system;
      Hastings MC(&boson_system, N, dim);
      monte_carlo_simulation(&MC, learning_rate, my_rank, num_procs);
    } else {
      Psi_OB boson_system;
      boson_system.update_alpha(0.7);
      Hastings MC(&boson_system, N, dim);
      monte_carlo_simulation(&MC, learning_rate, my_rank, num_procs);
    }
  } else {
    if (correlated) {
      Psi_T boson_system;
      Metropolis MC(&boson_system, N, dim);
      monte_carlo_simulation(&MC, learning_rate, my_rank, num_procs);
    } else {
      Psi_OB boson_system;
      boson_system.update_alpha(0.7);
      Metropolis MC(&boson_system, N, dim);
      monte_carlo_simulation(&MC, learning_rate, my_rank, num_procs);
    }
  }
  MPI_Finalize();
  return 0;
}
