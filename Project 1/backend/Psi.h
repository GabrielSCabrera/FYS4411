#include "../matpak/Mat.h"

#ifndef PSI_H
#define PSI_H

class Psi {

  protected:

    double alpha;
    double beta;
    double a;
    double omega;
    double omega_z;
    double mass = 1;
    double hbar = 1;
    double c = 1;

  public:

    // CONSTRUCTORS

    Psi(double alpha, double beta, double a, double omega, double omega_z);

    // DESTRUCTOR

    ~Psi();

    // UPDATERS

    void update_alpha(double alpha);

    void update_beta(double beta);

    void update_a(double a);

    void update_omega(double omega);

    void update_omega_z(double omega_z);

    void update_mass(double mass);

    // ATTRIBUTE EXTRACTION

    double get_alpha();

    double get_beta();

    double get_a();

    // CALLING

    double operator()(Mat P);

    // CALCULATIONS

    double V_ext(double x, double y, double z);

    double* drift(double x, double y, double z);

    double greens_ratio(double x0, double y0, double z0,
                        double x1, double y1, double z1,
                        double Ddt);

    double Psi_ob(Mat P, int N);

    double Psi_c(Mat P, int N);

    double energy(Mat P);

    double phi(double x, double y, double z);

    double* grad_phi(double x, double y, double z);

    double laplace_phi(double x, double y, double z);

    double u_prime(double r_jk);

    double u_double_prime(double r_jk);

    double grad_alpha(Mat P);

    double grad_beta(Mat P);

    double grad_alpha_alpha(Mat P);

    double grad_beta_beta(Mat P);

    double grad_beta_alpha(Mat P);

};

#endif
