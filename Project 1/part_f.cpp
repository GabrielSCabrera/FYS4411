#include <string>
#include <iostream>
#include <mpi.h>
#include <fstream>
#include "./wavefunctions/Psi_T.h"
#include "./variational/metropolis.h"
#include "./variational/gradient_descent.h"
#include "./matpak/Mat.h"

int main(int narg, char** argv) {
	double my_E, my_acceptance_ratio, dE, my_var, my_alpha;
	double E, acceptance_ratio, var;
	double* E_L;

	int my_rank, num_procs;
  	MPI_Init(&narg, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  	//-----change N----------: N = {10, 50, 100}
	int N = 10; 
	//-------------------------

	double alpha = 0.5;

	int cycles = 1e6/num_procs; // cycles per proc
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

	//sample energy
	R = MC.equilibriation(R, 100);
	R = MC.sample_energy(R, cycles);
	MC.print_info();

	// communicate results
	my_E = MC.get_energy();
	my_acceptance_ratio = MC.get_accepted_moves_ratio();
	MPI_Reduce(&my_E, &E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&my_acceptance_ratio, &acceptance_ratio, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&E, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Variance
	E_L = MC.get_E_cycles();
	E /= num_procs;
	my_var = 0.0;
	for (int j = 0; j < cycles; j++) {
		dE = E_L[j] - E;
		my_var += dE*dE;
	}
	my_var /= cycles;
	MPI_Reduce(&my_var, &var, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// write to file
	if (my_rank == 0) {
		std::ofstream outfile;
		std::string filename = "results/part_f/N_";
		filename.append(std::to_string(N));
		filename.append(".dat");
  		outfile.open(filename);
  		outfile << "cycles  " << cycles << "\n";
  		outfile << "workers " << num_procs << "\n";
		outfile << "\n\nalpha   " << alpha;
		outfile << "\nE       " << E;
		outfile << "\nE/Nd    " << E/(3.0*N);
		outfile << "\nvar     " << var;
		outfile << "\naccept  " << acceptance_ratio;
		outfile.close();
	}
	MPI_Finalize();
}
