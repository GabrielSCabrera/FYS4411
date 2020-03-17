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
  update_constants();
}

void Psi::update_beta(double beta) {
  this->beta = beta;
  this->beta_squared = beta*beta;
  update_constants();
}

void Psi::update_a(double a) {
  this->a = a;
  this->a_squared = a*a;
}

void Psi::update_gamma(double gamma) {
  this->gamma = gamma;
  this->gamma_squared = gamma*gamma;
}

void Psi::update_constants() {
  // used in grad phi
  minus_two_alpha = -2*alpha;
  minus_two_alpha_beta = minus_two_alpha*beta;
  // used in drag force
  minus_four_alpha = 2*minus_two_alpha;
  minus_four_alpha_beta = minus_four_alpha*beta;
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
double Psi::greens_ratio(Mat* R_new, Mat* R_old, double dt, int k) {
  double K = D*dt;
  double* F_new = drift_force(R_new, k);
  double* F_old = drift_force(R_old, k);

  double term;
  double exponent = 0.0;
  for (int l = 0; l < R_new->shape1(); l++) {
    term = R_old->get(k, l) - R_new->get(k, l) - K*F_new[l];
    exponent -= term*term;
    term = R_new->get(k, l) - R_old->get(k, l) - K*F_old[l];
    exponent += term*term;
  }
  delete [] F_old;
  delete [] F_new;
  return std::exp(exponent/(4.0*K));
}


double Psi::grad_alpha(Mat* R) {
  int M = R->shape1();
  double grad_alpha_phi = 0.0;
  double x;
  for (int i = 0; i < R->shape0(); i++) {
    x = R->get(i, 0);
    grad_alpha_phi += x*x;
    if (M > 1) {
      x = R->get(i, 1);
      grad_alpha_phi += x*x;
      if (M == 3) {
        x = R->get(i, 2);
        grad_alpha_phi += beta*x*x;
      }
    }
  }
  return -grad_alpha_phi;
}
