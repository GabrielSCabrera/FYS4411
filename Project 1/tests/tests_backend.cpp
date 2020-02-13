#include <cassert>
#include <stdlib.h>

#include "../backend/Psi.h"
#include "../matpak/Mat.h"

void tests_Psi() {
  double alpha_0 = 11;
  double alpha_1 = 32;

  double beta_0 = 11;
  double beta_1 = 32;

  // Initializing Psi instance
  Psi P1(alpha_0, beta_0);

  // Testing get_alpha()
  assert(P1.get_alpha() == alpha_0);

  // Testing get_beta()
  assert(P1.get_beta() == beta_0);

  // Testing update_alpha()
  P1.update_alpha(alpha_1);

  // Testing get_alpha()
  assert(P1.get_alpha() == alpha_1);

  // Testing update_beta()
  P1.update_beta(beta_1);

  // Testing get_beta()
  assert(P1.get_beta() == beta_1);

}
