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
	double alpha;
	int** rho = new int*[3];
	int anticipated_max = 5;

	int my_rank, num_procs;
  	MPI_Init(&narg, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  	//-----change N------------
	int N = 10; 
	//-------------------------
	int cycles = 1e7/num_procs; // cycles per proc
	int equi_cycles = 1e3;

	// UNCORRELATED
	// set up
	Psi_OB bose_system_uncorrelated;
	Metropolis MC_uncorrelated(&bose_system_uncorrelated, N, 3);
	Mat R = MC_uncorrelated.get_initial_R();
	R = MC_uncorrelated.equilibriation(R, equi_cycles);
	// CALCULATE ONE BODY DENSITY
	R = MC_uncorrelated.equilibriation(R, equi_cycles);
	MC_uncorrelated.one_body_density(R, cycles, anticipated_max);
	int N_rho;
	int** my_rho = MC_uncorrelated.get_rho(&N_rho);
	if (my_rank == 0) {
		rho[0] = new int [N_rho];
		rho[1] = new int [N_rho];
		rho[2] = new int [N_rho];
	}
	MPI_Reduce(my_rho[0], rho[0], N_rho, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(my_rho[1], rho[1], N_rho, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(my_rho[2], rho[2], N_rho, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (my_rank == 0) {
		std::ofstream outfile;
		std::string filename = "results/part_g/N_";
		filename.append(std::to_string(N));
		filename.append("_OB.dat");
  		outfile.open(filename);
  		outfile << anticipated_max << "\n";
		outfile << "x  y  z";
		for (int j = 0; j < N_rho; j++) {
			outfile << "\n" << rho[0][j] << " " <<rho[1][j] << " " <<rho[2][j];
		}
		outfile.close();
	}
	// CORRELATED

	// set up
	Psi_T bose_system_correlated;
	Metropolis MC_correlated(&bose_system_correlated, N, 3);
	R = MC_correlated.get_initial_R();
	R = MC_correlated.equilibriation(R, equi_cycles);

	alpha = 0.4975;
	MC_correlated.bose->update_alpha(alpha);

	// CALCULATE ONE BODY DENSITY
	R = MC_correlated.equilibriation(R, equi_cycles);
	MC_correlated.one_body_density(R, cycles, anticipated_max);

	my_rho = MC_correlated.get_rho(&N_rho);
	MPI_Reduce(my_rho[0], rho[0], N_rho, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(my_rho[1], rho[1], N_rho, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(my_rho[2], rho[2], N_rho, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_rank == 0) {
		std::ofstream outfile;
		std::string filename = "results/part_g/N_";
		filename.append(std::to_string(N));
		filename.append("_T.dat");
  		outfile.open(filename);
  		outfile << anticipated_max << "\n";
		outfile << "x  y  z";
		for (int j = 0; j < N_rho; j++) {
			outfile << "\n" << rho[0][j] << " " <<rho[1][j] << " " <<rho[2][j];
		}
		outfile.close();
		delete[] rho[0];
		delete[] rho[1];
		delete[] rho[2];
	}
	delete[] rho;
	MPI_Finalize();
}
