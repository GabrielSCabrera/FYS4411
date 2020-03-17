#include "../matpak/Mat.h"
#include "Psi.h"

#ifndef PSI_T_H
#define PSI_T_H

class Psi_T : public Psi {
    public:
    Psi_T();
    // CALLING
    double operator()(Mat R);
    // CALCULATIONS
    double* drift_force(Mat R, int index);

    double Psi_ob(Mat R, int N);
    double Psi_c(Mat R, int N);

    double energy(Mat R);
    double f(double r_ij);
    double u_prime(double r_ij);
    double u_double_prime(double r_ij, double u_prime_ij);
    double probability_density_ratio(Mat R_new, Mat R_old, int k);
    std::string name();
};

#endif
