#include <cmath>

#ifndef PSI_H
#define PSI_H

class Psi {

  protected:

    double alpha;

  public:

    // CONSTRUCTORS

    Psi(const double& alpha);

    // DESTRUCTOR

    ~Psi();

    // UPDATERS

    void update_alpha(const double& alpha);
};

#endif
