#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include <cmath>


#ifndef HASTINGS_H
#define HASTINGS_H

class Hastings : public Monte_Carlo {
protected:
	double dt = 0.01;
	double dt_sqrt = std::sqrt(dt);
	double D = 0.5;
public:
	using Monte_Carlo::Monte_Carlo;
	double acceptance_ratio(Mat* R_new, Mat* R_old, int k);
	void random_walk(Mat* R, int k);
	std::string filename_E();
	std::string filename_val();
};
#endif
