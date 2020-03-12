#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo_class.h"


#ifndef METROPOLIS_CLASS_H
#define METROPOLIS_CLASS_H

class Metropolis : public Monte_Carlo {
public:
	using Monte_Carlo::Monte_Carlo;
	double acceptance_ratio(double psi_new, double psi_old, Mat R_new, Mat R_old, int index);
	Mat random_walk(Mat R, int index);
};

#endif
