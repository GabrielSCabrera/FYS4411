#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include "gradient_descent.h"
#include <iostream>

// not really gradient decent, but Newton's method
Mat gradient_descent(Monte_Carlo* MC, Mat R) {
  double E_alpha = 42.0;
  int equi_cycles = 100;
  int sample_cycles = 1E3;
  int max_steps = 1E3;
  double tol = 1e-9;

  double alpha = MC->bose->get_alpha();
  int counter = 0;
  while (E_alpha*E_alpha > tol && counter < max_steps) {
    R = MC->equilibriation(R, equi_cycles);
    R = MC->sample_variational_derivatives(R, sample_cycles);

    E_alpha = MC->get_grad_alpha();
    alpha -= E_alpha/MC->E_2alpha;

    MC->bose->update_alpha(alpha);
    counter++;
  }
  if (E_alpha*E_alpha < tol) {printf("stopped because tolerance reached: %d\n", counter);}
  printf("alpha:%.12lf  ->  E': %.3e\n", alpha, E_alpha);
  return R;
}
