#include "../matpak/Mat.h"

#ifndef PSI_H
#define PSI_H

class Psi {

  protected:

    double alpha;
    double beta;
    double a;

  public:

    // CONSTRUCTORS

    Psi(const double& alpha, const double& beta, const double& a);

    // DESTRUCTOR

    ~Psi();

    // UPDATERS

    void update_alpha(const double& alpha);

    void update_beta(const double& beta);

    void update_a(const double& a);

    // ATTRIBUTE EXTRACTION

    double get_alpha();

    double get_beta();

    double get_a();

    // CALLING

    double operator()(Mat particles);

    // CALCULATIONS

    double Psi_ob(Mat P, int N);

    double Psi_c(Mat P, int N);

    double energy();

    double var_phi(double x, double y, double z);

};

#endif
