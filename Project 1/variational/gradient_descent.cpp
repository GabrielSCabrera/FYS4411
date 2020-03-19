#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include "gradient_descent.h"
#include <iostream>

Mat gradient_descent(Monte_Carlo* MC, double eta, Mat R) {
  // these should be parameters...
  int equi_cycles = 50;
  int sample_cycles = 100;
  int max_steps = 10E3;
  double tol = 1e-16;

  // first iteration
  R = MC->sample_variational_derivatives(R, sample_cycles);

  double E_alpha = MC->get_grad_alpha();
  eta = MC->E_2alpha;
  double alpha = MC->bose->get_alpha() - E_alpha/eta;
  double change = E_alpha*E_alpha;
  
  MC->bose->update_alpha(alpha);
  int counter = 0;
  while (change > tol && counter < max_steps) {
    R = MC->equilibriation(R, equi_cycles);
    R = MC->sample_variational_derivatives(R, sample_cycles);
    //if (counter % 100 == 0) {
    //printf("alpha: %.6lf,  E: %.6lf,  change: %.2e, %.2e\n", alpha, MC->get_energy_mean(), E_alpha, eta);
    //}
    E_alpha = MC->get_grad_alpha();
    eta = MC->E_2alpha;
    alpha -= E_alpha/eta;
    MC->bose->update_alpha(alpha);

    change = E_alpha*E_alpha;
    counter++;
  }
  if (change < tol) {printf("stopped because tolerance reached: %d\n", counter);}
  printf("alpha: %.12lf  ->  change: %.3e\n", alpha, change);
  return R;
}
