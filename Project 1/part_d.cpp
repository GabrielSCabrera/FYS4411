#include <iostream>
#include <fstream>
#include <string>
#include <mpi.h>
#include "./wavefunctions/Psi_OB.h"
#include "./variational/hastings.h"


int main(int narg, char** argv) {
	std::ofstream outfile;
	std::string filename;
	std::string path = "results/part_d/alpha";
	std::string end;
	double my_E, my_acceptance_ratio, dE, my_var;
	double E, acceptance_ratio, var;
	double* E_L;
	double* Es = nullptr;

	int my_rank, num_procs;
  	MPI_Init(&narg, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 

  	//-----change N----------: N = {10, 50, 100, 500}
	int N = 10; 
	//-------------------------
	double* alphas = new double [3];
	alphas[0] = 0.45; alphas[1] = 0.5; alphas[2] = 0.55;

	double* delta_t = new double [3];
	delta_t[0] = 0.001; delta_t[1] = 0.05; delta_t[2] = 0.01;

	int dim = 3;
	int cycles = 1e3/num_procs; // cycles per proc
	int equi_cycles = 1e3;
	if (my_rank == 0) {
	 	Es = new double [num_procs*cycles];
	}
	Psi_OB bose_system;
	double dt;
	Hastings MC(&bose_system, N, dim);
	Mat R = MC.get_initial_R();
	R = MC.equilibriation(R, equi_cycles);
	for (int t = 0; t < 3; t++) {
		dt = delta_t[t];
		end = "N_";
		end.append(std::to_string(N));
		end.append("_dt_");
		end.append(std::to_string(t));
		end.append(".dat");
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
			MPI_Gather(E_L, cycles, MPI_DOUBLE, Es, cycles, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (my_rank == 0) {
				filename = path;
				filename.append(std::to_string(i));
				filename.append("/val_");
				filename.append(end);
		  		outfile.open(filename);
		  		outfile << "cycles  " << cycles;
		  		outfile << "\nworkers " << num_procs;
				outfile << "\nalpha   " << alphas[i];
				outfile << "\nE       " << E;
				outfile << "\nE/Nd    " << E/(dim*N);
				outfile << "\nvar     " << var/num_procs;
				outfile << "\naccept  " << acceptance_ratio/num_procs;
				outfile << "\ndt      " << dt;
				outfile.close();

				filename = path;
				filename.append(std::to_string(i));
				filename.append("/E_");
				filename.append(end);
		  		outfile.open(filename);
		  		outfile.open(filename);
		  		outfile << Es[0];
				for (int j = 0; j < cycles*num_procs; j++) {
					outfile << "\n" << Es[j];
				}
				outfile.close();
				printf("\n");
			}
		}
	}
	if (my_rank == 0) {
	 	delete [] Es;
	}
	delete[] alphas;
	delete[] delta_t;
	MPI_Finalize();
	return 0;
}