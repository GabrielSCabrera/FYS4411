#include "../matpak/Mat.h"

#ifndef PSI_H
#define PSI_H

class Psi {
  protected:
    double D = 0.5;                 // diffusion constant
    double a, a_squared;            // hard radius
    double gamma, gamma_squared;    // elongation factor
    double alpha, alpha_squared;    // variational parameter
    double beta, beta_squared;      // variational parameter
    double minus_two_alpha, minus_two_alpha_beta, 
            minus_four_alpha, minus_four_alpha_beta;
    void update_constants();
  public:
    // CONSTRUCTORS
    Psi();
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
    virtual double probability_density_ratio(Mat R_new, Mat R_old, int k) = 0; 
    double greens_ratio(Mat R_old, Mat R_new, double dt, int index);

    double grad_alpha(Mat R);

    virtual std::string name() = 0;
};

#endif
