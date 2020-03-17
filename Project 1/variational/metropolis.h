#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"


#ifndef METROPOLIS_CLASS_H
#define METROPOLIS_CLASS_H

class Metropolis : public Monte_Carlo {
protected:
	double step_length = 0.5;
public:
	using Monte_Carlo::Monte_Carlo;
	double acceptance_ratio(Mat* R_new, Mat* R_old, int k);
	void random_walk(Mat *R, int k);
	std::string filename_E();
	std::string filename_val();
};
#endif
