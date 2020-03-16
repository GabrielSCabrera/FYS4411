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
double Psi_OB::operator()(Mat R) {
  double exponent = 0.0;
  double x;
  int last_index = R.shape1()-1;
  for (int i = 0; i < R.shape0(); i++) {
    for (int j = 0; j < last_index; j++) {
      x = R.get(i, j);
      exponent += x*x;
    }
    x = R.get(i, last_index);
    exponent += beta*x*x;
  }
  return std::exp(-alpha*exponent);
}

// CALLING
double Psi_OB::probability_density_ratio(Mat R_new, Mat R_old, int k) {
  double phi = 0.0;
  double x;
  int last_index = R_new.shape1() - 1;

  // calculate:  phi(r^new_k) - phi(r^old_k)
  for (int l = 0; l < last_index; l++) {
    x = R_new.get(k, l);
    phi += x*x;

    x = R_old.get(k, l);
    phi -= x*x;
  }
  x = R_new.get(k, last_index);
  phi += beta*x*x;

  x = R_old.get(k, last_index);
  phi -= beta*x*x;

  phi = std::exp(-alpha*phi);
  return phi*phi;
}

double* Psi_OB::drift_force(Mat R, int index) {
  int last_index = R.shape1()-1;
  double* force = new double [R.shape1()];
  for (int i = 0; i < last_index; i++) {
    force[i] = -4*alpha*R.get(index, i);
  }
  force[last_index] = -4*alpha*beta*R.get(index, last_index);
  return force;
}

double Psi_OB::energy(Mat R) {
  double E = 0.0;   // Kinetic Energy
  double V = 0.0;   // Potential
  int N = R.shape0();
  int last_index = R.shape1()-1;
  double x, xx;
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < last_index; j++) {
      x = R.get(k, j);
      xx = x*x;
      E += xx;
      V += xx;
    }
    x = R.get(k, last_index);
    xx = x*x;
    V += gamma_squared*xx;
    E += beta_squared*xx;
  } // END LOOP OVER k
  return N*alpha*(last_index + beta) + 0.5*V - 2*alpha_squared*E;
}

std::string Psi_OB::name() {
  std::string name = "delightful";
  return name;
}
