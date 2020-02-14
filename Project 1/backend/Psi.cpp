#include "Psi.h"
#include "../matpak/Mat.h"

// CONSTRUCTOR

Psi::Psi(const double& alpha, const double& beta) {
  this->alpha = alpha;
  this->beta = beta;
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

// ATTRIBUTE EXTRACTION

double Psi::get_alpha() {
  return this->alpha;
}

double Psi::get_beta() {
  return this->beta;
}
