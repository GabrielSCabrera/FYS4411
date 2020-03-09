#include "../matpak/Mat.h"

#ifndef PSI_H
#define PSI_H

class Psi {

  protected:

    double alpha;
    double alpha_squared;

    double beta;
    double beta_squared;

    double a;
    double a_squared;

    double gamma;
    double gamma_squared;

    double mass = 1;
    double hbar = 1;
    double c = 1;

  public:

    // CONSTRUCTORS

    Psi(double alpha, double beta, double a, double gamma, double mass);

    // DESTRUCTOR

    ~Psi();

    // UPDATERS

    void update_alpha(double alpha);

    void update_beta(double beta);

    void update_a(double a);

    void update_gamma(double gamma);

    void update_mass(double mass);

    // ATTRIBUTE EXTRACTION

    double get_alpha();

    double get_beta();

    double get_a();

    double get_gamma();

    double get_mass();

    // CALLING

    virtual double operator()(Mat P) = 0;

    // CALCULATIONS

    double V_ext(double x, double y, double z);

    virtual double* drift(double x, double y, double z) = 0;

    double greens_ratio(double x0, double y0, double z0,
                        double x1, double y1, double z1,
                        double Ddt);

    virtual double energy(Mat P) = 0;

    double phi(double x, double y, double z);

    double* grad_phi(double x, double y, double z);

    double laplace_phi(double x, double y, double z);

    double u_prime(double r_jk);

    double u_double_prime(double r_jk);

    virtual double grad_alpha(Mat P) = 0;

    virtual double grad_beta(Mat P) = 0;

    virtual double grad_alpha_alpha(Mat P) = 0;

    virtual double grad_beta_beta(Mat P) = 0;

};

#endif
