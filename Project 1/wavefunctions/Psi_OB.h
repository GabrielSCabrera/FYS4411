#include "../matpak/Mat.h"
#include "Psi.h"

#ifndef PSI_OB_H
#define PSI_OB_H

class Psi_OB : public Psi  {
  public:
    Psi_OB();
    std::string name();
    bool interaction();
    // CALCULATIONS
    double psi(Mat* R);
    double* drift_force(Mat* R, int k);
    double energy(Mat* R);
    double probability_density_ratio(Mat* R_new, Mat* R_old, int k);
};

#endif
