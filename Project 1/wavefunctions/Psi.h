#include "../matpak/Mat.h"

#ifndef PSI_H
#define PSI_H

class Psi {
  protected:
    double a, a_squared;            // hard radius
    double gamma, gamma_squared;    // elongation factor
    double alpha, alpha_squared;    // variational parameter
    double beta, beta_squared;      // variational parameter
    double D = 0.5;                 // diffusion constant
  public:
    // CONSTRUCTORS
    Psi(double alpha, double beta, double a, double gamma);
    // DESTRUCTOR
    ~Psi();
    // UPDATERS
    void update_alpha(double alpha);
    void update_beta(double beta);
    void update_a(double a);
    void update_gamma(double gamma);
    // ATTRIBUTE EXTRACTION
    double get_alpha();
    double get_beta();
    double get_a();
    double get_gamma();
    // CALLING
    virtual double operator()(Mat R) = 0;
    virtual double energy(Mat R) = 0;
    // CALCULATIONS
    virtual double* drift_force(Mat R, int index) = 0;
    double greens_ratio(Mat R_old, Mat R_new, double dt, int index);

    double phi(double x, double y, double z);
    double* grad_phi(double x, double y, double z);
    double laplace_phi(double x, double y, double z);

    double grad_alpha(Mat R);
    double grad_beta(Mat R);
};

#endif
