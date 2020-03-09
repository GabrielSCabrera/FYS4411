#include <cmath>

#include "Psi.h"
#include "Psi_OB.h"
#include "../matpak/Mat.h"

// CALLING

double Psi_OB::operator()(Mat P) {
  /*
    P – Array of type 'Mat' of shape (N,3)
    N – Number of particles
  */
  double exponent = 0;
  double x; double y; double z;
  for (int i = 0; i < P.shape0(); i++) {
    x = P.get(i,0);
    y = P.get(i,1);
    z = P.get(i,2);
    exponent += x*x + y*y + this->beta*z*z;
  }
  exponent = std::exp(-this->alpha*exponent);
  return exponent;
}

// CALCULATIONS

double* Psi_OB::drift(double x, double y, double z) {
  double* force = new double [3];
  force[0] = -4*this->alpha*x;
  force[1] = -4*this->alpha*y;
  force[2] = -4*this->alpha*this->beta*z;
  return force;
}

double Psi_OB::energy(Mat P) {

  // Kinetic Energy
  double E = 0;

  // Potential
  double V = 0;

  double x_k; double y_k; double z_k;

  for (int k = 0; k < P.shape0(); k++) {
    x_k = P.get(k,0); y_k = P.get(k,1); z_k = P.get(k,2);
    V += V_ext(x_k, y_k, z_k);
    E += x_k*x_k + y_k*y_k + this->beta_squared*z_k*z_k;

  } // END LOOP OVER k

  return P.shape0()*this->alpha*(P.shape1() - 1 + this->beta) + V - 2*this->alpha_squared*E;
}

double Psi_OB::grad_alpha(Mat P) {

  // Kinetic Energy
  double E = 0;

  double x_k; double y_k; double z_k;

  for (int k = 0; k < P.shape0(); k++) {
    x_k = P.get(k,0); y_k = P.get(k,1); z_k = P.get(k,2);
    E += x_k*x_k + y_k*y_k + this->beta_squared*z_k*z_k;
  }

  return P.shape0()*(P.shape1() - 1 + this->beta) - 4*this->alpha*E;
}

double Psi_OB::grad_beta(Mat P) {

    // Kinetic Energy
    double E = 0;

    double z_k;

    for (int k = 0; k < P.shape0(); k++) {
      z_k = P.get(k,2);
      E += z_k*z_k;
    }

    return P.shape0()*this->alpha - 4*this->alpha_squared*this->beta*E;
  }

double Psi_OB::grad_alpha_alpha(Mat P) {

    // Kinetic Energy
    double E = 0;

    for (int k = 0; k < P.shape0(); k++) {
      E += V_ext(P.get(k,0), P.get(k,1), P.get(k,2));
    } // END LOOP OVER k

    return -8*E;
}

double Psi_OB::grad_beta_beta(Mat P) {
  return 0;
}
