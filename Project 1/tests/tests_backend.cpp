#include <stdlib.h>
#include <cassert>
#include <cmath>

#include "../backend/Psi.h"
#include "../matpak/Mat.h"

void tests_Psi() {

  double tol = 1E-7;

  double alpha_0 = 1.1;
  double alpha_1 = 1;

  double beta_0 = 13.21;
  double beta_1 = 1;

  double a_0 = 1E-3;
  double a_1 = 1E-5;

  // Initializing Psi instance
  Psi Wf1(alpha_0, beta_0, a_0, 1, 1);

  // Testing get_alpha()
  assert(Wf1.get_alpha() == alpha_0);

  // Testing get_beta()
  assert(Wf1.get_beta() == beta_0);

  // Testing get_a()
  assert(Wf1.get_a() == a_0);

  // Testing update_alpha()
  Wf1.update_alpha(alpha_1);

  // Testing get_alpha()
  assert(Wf1.get_alpha() == alpha_1);

  // Testing update_beta()
  Wf1.update_beta(beta_1);

  // Testing get_beta()
  assert(Wf1.get_beta() == beta_1);

  // Testing update_a()
  Wf1.update_a(a_1);

  // Testing get_a()
  assert(Wf1.get_a() == a_1);

  // Initializing Matrix of Particles
  int N = 5; int d = 3;
  Mat P1(N,d);
  P1 = "[[0,0,0],[0.1,0.1,0.1],[0.2,0.2,0.2],[0.3,0.3,0.3],[0.4,0.4,0.4]]";

  // Testing Psi_ob(Particles, N)
  assert((0.4065696597405991 - Wf1.Psi_ob(P1,N)) <= tol);

  // Testing Psi_c(Particles, N)
  assert((-42.99213492801157 - Wf1.Psi_c(P1,N)) <= tol);

}
