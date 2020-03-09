#include "../matpak/Mat.h"
#include "Psi.h"

#ifndef PSI_T_H
#define PSI_T_H

class Psi_T : public Psi {

  public:

    using Psi::Psi;

    // CALLING

    double operator()(Mat P);

    // CALCULATIONS

    double* drift(double x, double y, double z);

    double Psi_ob(Mat P, int N);

    double Psi_c(Mat P, int N);

    double energy(Mat P);

    double grad_alpha(Mat P);

    double grad_beta(Mat P);

    double grad_alpha_alpha(Mat P);

    double grad_beta_beta(Mat P);

};

#endif
