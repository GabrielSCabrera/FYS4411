#include <cmath>

#include "Psi.h"
#include "../matpak/Mat.h"

// CONSTRUCTOR

Psi::Psi(double alpha, double beta, double a, double gamma, double mass) {
  update_alpha(alpha);
  update_beta(beta);
  update_a(a);
  update_gamma(gamma);
  update_mass(mass);
}

// DESTRUCTOR

Psi::~Psi() {

}

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

void Psi::update_mass(double mass) {
  this->mass = mass;
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

double Psi::get_mass() {
  return this->mass;
}

// CALCULATIONS

double Psi::V_ext(double x, double y, double z) {
  return 0.5*(x*x + y*y + this->gamma_squared*z*z);
}

double Psi::greens_ratio(double x0, double y0, double z0,
                         double x1, double y1, double z1,
                         double K) {
  // K = D*dt

  double* terms = new double[3];

  double* F1 = drift(x1, y1, z1);
  terms[0] = x0-x1-K*F1[0]; terms[1] = y0-y1-K*F1[1]; terms[2] = z0-z1-K*F1[2];
  double exponent = terms[0]*terms[0] + terms[1]*terms[1] + terms[2]*terms[2];
  delete [] F1;

  double* F2 = drift(x0, y0, z0);
  terms[0] = x1-x0-K*F2[0]; terms[1] = y1-y0-K*F2[1]; terms[2] = z1-z0-K*F2[2];
  exponent -= terms[0]*terms[0] + terms[1]*terms[1] + terms[2]*terms[2];

  delete [] F2;
  delete [] terms;
  return std::exp(exponent/(4*K));
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

double Psi::u_prime(double r_kj) {
  if (r_kj <= this-> a) {
    return 0;
  } else {
    return this->a/(r_kj*(r_kj - this->a));
  }
}

double Psi::u_double_prime(double r_kj) {
  /*
    !! Must be multiplied with u_prime(r_kj) !!
  */
  if (r_kj <= this-> a) {
    return 0;
  } else {
    return -(1/(r_kj - this->a) + 1/r_kj);
  }
}
