#include <string>
#include <iostream>
#include <mpi.h>
#include <fstream>
#include "./wavefunctions/Psi_T.h"
#include "./wavefunctions/Psi_OB.h"
#include "./variational/metropolis.h"
#include "./variational/gradient_descent.h"
#include "./variational/monte_carlo.h"
#include "./matpak/Mat.h"

int main(int narg, char** argv) {
	double my_alpha;

	int my_rank, num_procs;
  	MPI_Init(&narg, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 

  	//-----change N------------
	int N = 2; 
	//-------------------------

	double alpha = 0.5;
	int cycles = 1e4/num_procs; // cycles per proc
	int equi_cycles = 1e3;

	// set up
	Psi_T bose_system;
	bose_system.update_alpha(alpha);
	Metropolis MC(&bose_system, N, 3);
	Mat R = MC.get_initial_R();
	R = MC.equilibriation(R, equi_cycles);
	
	// Find alpha (choose mean value of all alphas)
	R = gradient_descent(&MC, 0.0, R);
	my_alpha = MC.bose->get_alpha();
	MPI_Reduce(&my_alpha, &alpha, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	alpha /= num_procs;
	MC.bose->update_alpha(alpha);

	// CALCULATE ONE BODY DENSITY....?


	// write to file
	/*
	if (my_rank == 0) {
		std::ofstream outfile;
		std::string filename = "results/part_g/";
		filename.append(".dat");
  		outfile.open(filename);
		outfile.close();
	}
	*/
	MPI_Finalize();
}