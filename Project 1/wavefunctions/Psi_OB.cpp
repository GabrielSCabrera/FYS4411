#include <cmath>
#include "Psi.h"
#include "Psi_OB.h"
#include "../matpak/Mat.h"


// CONSTRUCTOR
Psi_OB::Psi_OB() {
  update_alpha(0.5);
  update_beta(1.0);
  update_a(0.0);
  update_gamma(1.0);
}

// CALLING
double Psi_OB::psi(Mat* R) {
  int M = R->shape1();
  double exponent = 0.0;
  double x;
  for (int i = 0; i < R->shape0(); i++) {
    x = R->get(i, 0);
    exponent += x*x;
    if (M > 1) {
      x = R->get(i, 1);
      exponent += x*x;
      if (M == 3) {
        x = R->get(i, 2);
        exponent += beta*x*x;
      }
    }
  }
  return std::exp(-alpha*exponent);
}

double Psi_OB::probability_density_ratio(Mat* R_new, Mat* R_old, int k) {
  double phi = 0.0;
  double x;
  int M = R_new->shape1();
  // calculate:  phi(r^new_k) - phi(r^old_k) 
  x = R_new->get(k, 0);
  phi -= x*x;
  x = R_old->get(k, 0);
  phi += x*x;
  if (M > 1) {        // y-axis
    x = R_new->get(k, 1);
    phi -= x*x;
    x = R_old->get(k, 1);
    phi += x*x;
    if (M == 3) {     // z-axis
      x = R_new->get(k, 2);
      phi -= beta*x*x;
      x = R_old->get(k, 2);
      phi += beta*x*x;
    }
  }
  phi = std::exp(alpha*phi);
  return phi*phi;
}

double* Psi_OB::drift_force(Mat* R, int k) {
  int M = R->shape1();
  double* force = new double [M];
  force[0] = minus_four_alpha*R->get(k, 0);
  if (M > 1) {
    force[1] = minus_four_alpha*R->get(k, 1);
    if (M == 3) {
      force[2] = minus_four_alpha_beta*R->get(k, 2);
    }
  }
  return force;
}

double Psi_OB::energy(Mat* R) {
  double E = 0.0;   // Kinetic Energy
  double V = 0.0;   // Potential
  int N = R->shape0();
  int M = R->shape1();
  double x, xx;
  for (int k = 0; k < N; k++) {
    x = R->get(k, 0);
    xx = x*x;
    E += xx;
    V += xx;
    if (M > 1) {
      x = R->get(k, 1);
      xx = x*x;
      E += xx;
      V += xx;
      if (M == 3) {
        x = R->get(k, 2);
        xx = x*x;
        V += gamma_squared*xx;
        E += beta_squared*xx;
      }
    }
  }
  double M_beta = (M == 3) ? 2.0 - beta : M;
  return N*alpha*M_beta + 0.5*V - 2*alpha_squared*E;
}

std::string Psi_OB::name() {
  std::string name = "delightful";
  return name;
}

bool Psi_OB::interaction() {
  return false;
}
