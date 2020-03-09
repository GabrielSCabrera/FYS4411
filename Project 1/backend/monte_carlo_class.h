#include "../matpak/Mat.h"
#include "Psi.h"

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

class Monte_Carlo {
protected:
	int N, dim;
	double step_length, x_max;
	Psi* PDF;

public:
	Monte_Carlo::Monte_Carlo;
	virtual Mat random_walk() = 0;

};

#endif