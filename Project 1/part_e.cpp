#include <string>
#include <iostream>
#include <mpi.h>
#include <fstream>
#include "./wavefunctions/Psi_T.h"
#include "./variational/metropolis.h"
#include "./matpak/Mat.h"

int main(int narg, char** argv) {
	double my_E, my_acceptance_ratio, dE, my_var;
	double E, acceptance_ratio, var;
	double* E_L;
	double* acceptance_ratios = nullptr;
	double* variance = nullptr;
	double* Es = nullptr;

	int my_rank, num_procs;
  	MPI_Init(&narg, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 

  	//-----change N----------: N = {10, 50, 100}
	int N = 50; 
	//-------------------------


	double* alphas = new double [3];
	alphas[0] = 0.45; alphas[1] = 0.5; alphas[2] = 0.55;

	int cycles = 1e6/num_procs; // cycles per proc
	int equi_cycles = 1e3;

	Psi_T bose_system;
	Metropolis MC(&bose_system, N, 3);
	Mat R = MC.get_initial_R();
	R = MC.equilibriation(R, equi_cycles);
	
	if (my_rank == 0) {
	 	Es = new double [3];
	 	variance = new double[3];
	 	acceptance_ratios = new double [3];
	}

	for (int i = 0; i < 3; i++) {
		MC.bose->update_alpha(alphas[i]);
		R = MC.equilibriation(R, 100);
		R = MC.sample_energy(R, cycles);
		MC.print_info();
		my_E = MC.get_energy();
		my_acceptance_ratio = MC.get_accepted_moves_ratio();

		MPI_Reduce(&my_E, &E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&my_acceptance_ratio, &acceptance_ratio, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Bcast(&E, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		E_L = MC.get_E_cycles();
		E /= num_procs;
		my_var = 0.0;
		// Variance
		for (int j = 0; j < cycles; j++) {
			dE = E_L[j] - E;
			my_var += dE*dE;
		}
		my_var /= cycles;
		MPI_Reduce(&my_var, &var, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (my_rank == 0) {
			Es[i] = E;
			acceptance_ratios[i] = acceptance_ratio/num_procs;
			variance[i] = var/num_procs; 
			printf("\n");
		}
	}
	// write to file
	if (my_rank == 0) {
		std::ofstream outfile;
		std::string filename = "results/part_e/N_";
		filename.append(std::to_string(N));
		filename.append(".dat");
  		outfile.open(filename);
  		outfile << "cycles  " << cycles << "\n";
  		outfile << "workers " << num_procs << "\n";
  		for (int i = 0; i < 3; i++) {
			outfile << "\n\nalpha   " << alphas[i];
			outfile << "\nE       " << Es[i];
			outfile << "\nE/Nd    " << Es[i]/(3.0*N);
			outfile << "\nvar     " << variance[i];
			outfile << "\naccept  " << acceptance_ratios[i];
		}
		outfile.close();
 		delete [] Es;
 		delete [] variance;
 		delete [] acceptance_ratios;
	}
	delete[] alphas;
	MPI_Finalize();
}