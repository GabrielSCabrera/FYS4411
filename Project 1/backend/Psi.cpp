#include <iostream>
#include <cmath>
#include "Psi.h"
#include "../matpak/Mat.h"

// CONSTRUCTOR

Psi::Psi(const double& alpha) {
  this->alpha = alpha;
}

// DESTRUCTOR

Psi::~Psi() {

}

// UPDATERS

void Psi::update_alpha(const double& alpha) {
  this->alpha = alpha;
}
