#include "../matpak/Mat.h"
#include "Psi.h"

#ifndef PSI_OB_H
#define PSI_OB_H

class Psi_OB : public Psi  {

  public:

    // CONSTRUCTORS

    Psi_OB(double alpha, double beta, double a, double omega, double omega_z);

    // DESTRUCTOR

    ~Psi_OB();

    // CALLING

    double operator()(Mat P);

    // CALCULATIONS

    double* drift(double x, double y, double z);

    double Psi_ob(Mat P, int N);

    double energy(Mat P);

    double grad_alpha(Mat P);

    double grad_beta(Mat P);

    double grad_alpha_alpha(Mat P);

    double grad_beta_beta(Mat P);

};

#endif
