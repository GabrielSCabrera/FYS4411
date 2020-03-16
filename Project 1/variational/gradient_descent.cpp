#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include "gradient_descent.h"
#include <iostream>

Mat gradient_descent(Monte_Carlo* MC, double eta, Mat R) {
  // these should be parameters...
  int equi_cycles = 100;
  int sample_cycles = 1E4;
  int max_steps = 50;
  double tol = 1e-10;

  // first iteration
  R = MC->equilibriation(R, equi_cycles);
  R = MC->sample_variational_derivatives(R, sample_cycles);

  double grad_alpha = MC->get_grad_alpha();
  double alpha = MC->PDF->get_alpha() - eta*grad_alpha;
  
  MC->PDF->update_alpha(alpha);
  int counter = 0;
  while (grad_alpha*grad_alpha > tol && counter < max_steps) {
    // Running Monte-Carlo
    R = MC->equilibriation(R, equi_cycles);
    R = MC->sample_variational_derivatives(R, sample_cycles);

    //printf("alpha: %.6lf,  E: %.6lf\n", alpha, MC->get_energy_mean());
    grad_alpha = MC->get_grad_alpha();
    alpha -= eta*grad_alpha;
    MC->PDF->update_alpha(alpha);

    counter++;
  }
  if (grad_alpha*grad_alpha < tol) {printf("stopped because tolerance reached\n");}

  printf("alpha: %.12lf  ->  change: %.3e\n", alpha, grad_alpha*grad_alpha);
  return R;
}
