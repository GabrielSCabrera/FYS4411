#include "../matpak/Mat.h"
#include "Psi.h"

#ifndef PSI_OB_H
#define PSI_OB_H

class Psi_OB : public Psi  {
  public:
    using Psi::Psi;
    // CALLING
    double operator()(Mat R);
    // CALCULATIONS
    double* drift_force(Mat R, int index);
    double energy(Mat R);
    double grad_alpha(Mat R);
    double grad_beta(Mat R);
};

#endif
