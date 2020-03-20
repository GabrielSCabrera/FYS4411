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
	double my_alpha, alpha;
	double* rho = nullptr;

	int my_rank, num_procs;
  	MPI_Init(&narg, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 

  	//-----change N------------
	int N = 10; 
	//-------------------------
	int cycles = 1e6/num_procs; // cycles per proc
	int equi_cycles = 1e3;

	// UNCORRELATED
	// set up
	Psi_OB bose_system_uncorrelated;
	Metropolis MC_uncorrelated(&bose_system_uncorrelated, N, 3);
	Mat R = MC_uncorrelated.get_initial_R();
	R = MC_uncorrelated.equilibriation(R, equi_cycles);
	printf("mkht\n");
	// CALCULATE ONE BODY DENSITY
	R = MC_uncorrelated.equilibriation(R, equi_cycles);
	MC_uncorrelated.one_body_density(R, cycles);

	int N_rho;
	double* my_rho = MC_uncorrelated.get_rho(&N_rho);
	if (my_rank == 0) {
		rho = new double [N_rho];
	}
	MPI_Reduce(my_rho, rho, N_rho, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	printf("m\n");
	if (my_rank == 0) {
		std::ofstream outfile;
		std::string filename = "results/part_g/N_";
		filename.append(std::to_string(N));
		filename.append("_OB.dat");
  		outfile.open(filename);
		outfile << "r in [0, 4]";
		for (int j = 0; j < N_rho; j++) {
			outfile << "\n" << rho[j];
		}
		outfile.close();
	}
	printf("made it\n");
	// CORRELATED

	// set up
	Psi_T bose_system_correlated;
	Metropolis MC_correlated(&bose_system_correlated, N, 3);
	R = MC_correlated.get_initial_R();
	R = MC_correlated.equilibriation(R, equi_cycles);
	
	// Find alpha (choose mean value of all alphas)
	R = gradient_descent(&MC_correlated, 0.0, R);
	my_alpha = MC_correlated.bose->get_alpha();
	MPI_Reduce(&my_alpha, &alpha, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	alpha /= num_procs;
	MC_correlated.bose->update_alpha(alpha);

	// CALCULATE ONE BODY DENSITY
	R = MC_correlated.equilibriation(R, equi_cycles);
	MC_correlated.one_body_density(R, cycles);

	my_rho = MC_correlated.get_rho(&N_rho);
	MPI_Reduce(my_rho, rho, N_rho, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_rank == 0) {
		std::ofstream outfile;
		std::string filename = "results/part_g/N_";
		filename.append(std::to_string(N));
		filename.append("_T.dat");
  		outfile.open(filename);
		outfile << "r in [0, 4]";
		for (int j = 0; j < N_rho; j++) {
			outfile << "\n" << rho[j];
		}
		outfile.close();
		delete[] rho;
	}
	MPI_Finalize();
}