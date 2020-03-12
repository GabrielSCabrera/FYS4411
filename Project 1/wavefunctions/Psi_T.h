#include "../matpak/Mat.h"
#include "Psi.h"

#ifndef PSI_T_H
#define PSI_T_H

class Psi_T : public Psi {
  public:
    using Psi::Psi;
    // CALLING
    double operator()(Mat R);
    // CALCULATIONS
    double* drift_force(Mat R, int index);

    double Psi_ob(Mat R, int N);

    double Psi_c(Mat R, int N);

    double energy(Mat R);
    double u_prime(double r_kj);
    double u_double_prime(double r_kj);
};

#endif
