#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo.h"
#include "gradient_descent.h"
#include <iostream>

Mat gradient_descent(Monte_Carlo* MC, double eta, Mat R) {
  double E_alpha = 42.0;
  double E_2alpha;
  int equi_cycles = 500;
  int sample_cycles = 1E4;
  int max_steps = 1E4;
  double tol = 1e-10;

  double alpha = MC->bose->get_alpha();
  int counter = 0;
  while (E_alpha*E_alpha > tol && counter < max_steps) {
    R = MC->equilibriation(R, equi_cycles);
    R = MC->sample_variational_derivatives(R, sample_cycles);

    E_alpha = MC->get_grad_alpha();
    E_2alpha = MC->E_2alpha;
    alpha -= E_alpha/E_2alpha;
    MC->bose->update_alpha(alpha);

    //if (counter % 100 == 0) {
    //printf("alpha:%.6lf  E:%.6lf  E':%9.2e  E'':%9.2e  E'/E'':%9.2e\n", \
     // alpha, MC->get_energy_mean(), E_alpha, E_2alpha, E_alpha/E_2alpha);
    //}
    counter++;
  }
  if (E_alpha*E_alpha < tol) {printf("stopped because tolerance reached: %d\n", counter);}
  printf("alpha:%.12lf  ->  E': %.3e\n", alpha, E_alpha);
  return R;
}
