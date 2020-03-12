#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo_class.h"


#ifndef IMPORTANCE_SAMPLING_CLASS_H
#define IMPORTANCE_SAMPLING_CLASS_H

class Importance_Sampling : public Monte_Carlo {
protected:
	double dt = 0.01;
	double dt_sqrt = 0.1;
	double D = 0.5;
public:
	using Monte_Carlo::Monte_Carlo;
	double acceptance_ratio(double psi_new, double psi_old, Mat R_new, Mat R_old, int index);
	Mat random_walk(Mat R, int index);
};

#endif
