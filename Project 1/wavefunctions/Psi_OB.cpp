#include <cmath>

#include "Psi.h"
#include "Psi_OB.h"
#include "../matpak/Mat.h"

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

// CALCULATIONS

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
  double x, xx;
  int last_index = R.shape1()-1;
  for (int k = 0; k < R.shape0(); k++) {
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
  return R.shape0()*alpha*(last_index + beta) + 0.5*V - 2*alpha_squared*E;
}


double Psi_OB::grad_alpha(Mat R) {
  double grad_alpha_phi = 0.0;
  double x_i;
  int last_index = R.shape1()-1;

  for (int i = 0; i < R.shape0(); i++) {
    for (int j = 0; j < last_index; j++) {
      x_i = R.get(i, j);
      grad_alpha_phi += x_i*x_i;
    }
    x_i = R.get(i, last_index);
    grad_alpha_phi += beta*x_i*x_i;
  }
  return -grad_alpha_phi;
}


double Psi_OB::grad_beta(Mat R) {
  double grad_beta_phi = 0.0;
  double z_k;
  int last_index = R.shape1()-1;

  for (int k = 0; k < R.shape0(); k++) {
    z_k = R.get(k, last_index);
    grad_beta_phi += z_k*z_k;
  }
  return -alpha*grad_beta_phi;
}

