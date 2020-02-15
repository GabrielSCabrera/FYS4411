#include <cmath>

#include "Psi.h"
#include "../matpak/Mat.h"

// CONSTRUCTOR

Psi::Psi(const double& alpha, const double& beta, const double& a) {
  this->alpha = alpha;
  this->beta = beta;
  this->a = a;
}

// DESTRUCTOR

Psi::~Psi() {

}

// UPDATERS

void Psi::update_alpha(const double& alpha) {
  this->alpha = alpha;
}

void Psi::update_beta(const double& beta) {
  this->beta = beta;
}

void Psi::update_a(const double& a) {
  this->a = a;
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

// CALLING

double Psi::operator()(Mat particles) {
  return 0;
}

// CALCULATIONS

double Psi::Psi_ob(Mat P, int N) {
  /*
    P – Array of type 'Mat' of shape (N,3)
    N – Number of particles
  */
  double exponent = 0;
  double x; double y; double z;
  for (int i = 0; i < N; i++) {
    x = std::pow(P.get(i,0), 2);
    y = std::pow(P.get(i,1), 2);
    z = std::pow(P.get(i,2), 2);
    exponent += x + y + (this->beta)*z;
  }
  exponent = std::exp(-(this->alpha)*exponent);
  return exponent;
}

double Psi::Psi_c(Mat P, int N) {
  /*
    P – Array of type 'Mat' of shape (N,3)
    N – Number of particles
  */
  double product = 1;
  double dx; double dy; double dz;
  double diff; int j;
  for (int i = 0; i < N; i++) {
    for (j = i+1; j < N; j++) {

      dx = std::pow(P.get(i,0) - P.get(j,0), 2);
      dy = std::pow(P.get(i,1) - P.get(j,1), 2);
      dz = std::pow(P.get(i,2) - P.get(j,2), 2);

      diff = std::sqrt(dx + dy + dz);

      if (diff <= 0) {
        product = 0;
        goto stop;
      } else {
        product *= 1 - (this->a)/diff;
      }
    }
  }
  stop:
  return product;
}
