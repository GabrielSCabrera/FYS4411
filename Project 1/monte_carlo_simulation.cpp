#include <string>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <fstream>
#include <ctime>
#include <chrono>

#include "./wavefunctions/Psi.h"
#include "./variational/monte_carlo.h"
#include "./matpak/Mat.h"
#include "./variational/gradient_descent.h"
#include "monte_carlo_simulation.h"

using namespace std;

void monte_carlo_simulation(Monte_Carlo* MC, double learning_rate, int my_rank, int num_procs) {
	double E, alpha, acceptance_ratio;
	double my_E, my_alpha, my_acceptance_ratio;
	double* my_E_cycles;
	double* E_cycles;

	int cycles = 1e6/num_procs; // cycles per proc
	int equi_cycles = 1e4;

	if (my_rank == 0) {	
		alpha = 0.0;
		E = 0.0;
		acceptance_ratio = 0.0;
		E_cycles = new double [cycles*num_procs];
	}

	Mat R = MC->get_initial_R();

	//--- this should not be here----
	R = MC->equilibriation(R, equi_cycles);
	R = MC->sample_energy(R, cycles);
	MC->print_info();
	//-------------------------

	// Find alpha
	R = gradient_descent(MC, learning_rate, R);
	my_alpha = MC->bose->get_alpha();
	MPI_Reduce(&my_alpha, &alpha, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	alpha /= num_procs;
	MC->bose->update_alpha(alpha);

	// sample energy
	R = MC->equilibriation(R, equi_cycles);
	R = MC->sample_energy(R, cycles);
	MC->print_info();

	// send info to rank zero
	my_E = MC->get_energy();
	my_acceptance_ratio = MC->get_accepted_moves_ratio();
	my_E_cycles = MC->get_E_cycles();
	MPI_Reduce(&my_E, &E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&my_acceptance_ratio, &acceptance_ratio, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Gather(my_E_cycles, cycles, MPI_DOUBLE, E_cycles, cycles, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// write to file
	if (my_rank == 0) {
		// write all energies to file 
		ofstream outfile;
		string filename = MC->filename_E();
		outfile.open(filename);
		outfile << E_cycles[0];
		for (int i = 1; i < cycles*num_procs; i++) {
			outfile << "\n" << E_cycles[i];
		}
		outfile.close();
		delete[] E_cycles;
		// write rest of info to file
		filename = MC->filename_val();
  		outfile.open(filename);
		outfile << "alpha " << alpha << "\n";
		outfile << "beta " << MC->bose->get_beta() << "\n";
		outfile << "E " << E/num_procs << "\n";
		outfile << "accept " << acceptance_ratio/num_procs << "\n";
		outfile << "cycles " << cycles;
		outfile.close();
	}
}