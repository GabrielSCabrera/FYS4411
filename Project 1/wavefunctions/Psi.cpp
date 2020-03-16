#include <cmath>

#include "Psi.h"
#include "../matpak/Mat.h"

// CONSTRUCTOR
Psi::Psi() {
  update_alpha(0.5);
  update_beta(1.0);
  update_a(0.0043);
  update_gamma(1.0);
}

// DESTRUCTOR
Psi::~Psi() {}

// UPDATERS
void Psi::update_alpha(double alpha) {
  this->alpha = alpha;
  this->alpha_squared = alpha*alpha;
}

void Psi::update_beta(double beta) {
  this->beta = beta;
  this->beta_squared = beta*beta;
}

void Psi::update_a(double a) {
  this->a = a;
  this->a_squared = a*a;
}

void Psi::update_gamma(double gamma) {
  this->gamma = gamma;
  this->gamma_squared = gamma*gamma;
}

// ATTRIBUTE EXTRACTION
double Psi::get_alpha() {
  return this->alpha;
}

double Psi::get_beta() {
  return this->beta;
}

double Psi::get_a() {
  return this->a;
}

double Psi::get_gamma() {
  return this->gamma;
}

// CALCULATIONS
double Psi::greens_ratio(Mat R_old, Mat R_new, double dt, int index) {
  // index : corresponds to particle
  double K = D*dt;
  double* F_old = drift_force(R_old, index);
  double* F_new = drift_force(R_new, index);

  double term;
  double exponent = 0.0;
  for (int i = 0; i < R_new.shape1(); i++) {
      term = R_old.get(index, i) - R_new.get(index, i) - K*F_new[i];
      exponent -= term*term;
      term = R_new.get(index, i) - R_old.get(index, i) - K*F_old[i];
      exponent += term*term;
  }

  delete [] F_old;
  delete [] F_new;
  return std::exp(exponent/(4.0*K));
}

double Psi::phi(double x, double y, double z) {
  return std::exp(-this->alpha*(x*x + y*y + this->beta*z*z));
}

double* Psi::grad_phi(double x, double y, double z) {
  /*
    !! Must be multiplied with phi(x, y, z) !!
  */
  double* grad = new double [3];
  grad[0] = -2*this->alpha*x;
  grad[1] = -2*this->alpha*y;
  grad[2] = -2*this->alpha*z*this->beta;
  return grad;
}

double Psi::laplace_phi(double x, double y, double z) {
  /*
    !! Must be multiplied with phi(x, y, z) !!
  */
  return 2*this->alpha*(2*this->alpha*(x*x+y*y+this->beta*this->beta*z*z) - 2 - this->beta);
}

double Psi::grad_alpha(Mat R) {
  double grad_alpha_phi = 0.0;
  double x_i;
  int last_index = R.shape1()-1;

  for (int i = 0; i < R.shape0(); i++) {
    for (int j = 0; j < last_index; j++) {
      x_i = R.get(i, j);
      grad_alpha_phi += x_i*x_i;
    }
    x_i = R.get(i, last_index);
    grad_alpha_phi += this->beta*x_i*x_i;
  }
  return -grad_alpha_phi;
}

double Psi::grad_beta(Mat R) {
  double grad_beta_phi = 0.0;
  double z_k;
  int last_index = R.shape1()-1;

  for (int k = 0; k < R.shape0(); k++) {
    z_k = R.get(k, last_index);
    grad_beta_phi += z_k*z_k;
  }
  return -this->alpha*grad_beta_phi;
}
