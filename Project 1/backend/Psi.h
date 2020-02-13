#ifndef PSI_H
#define PSI_H

class Psi {

  protected:

    double alpha;
    double beta;

  public:

    // CONSTRUCTORS

    Psi(const double& alpha, const double& beta);

    // DESTRUCTOR

    ~Psi();

    // UPDATERS

    void update_alpha(const double& alpha);

    void update_beta(const double& beta);

    // ATTRIBUTE EXTRACTION

    double get_alpha();

    double get_beta();

};

#endif
