#include "../matpak/Mat.h"
#include "../wavefunctions/Psi.h"
#include "monte_carlo_class.h"
#include "gradient_descent.h"
#include <iostream>

void gradient_descent(Monte_Carlo* MC, double eta) {
  double alpha, beta;
  double alpha_prev = MC->PDF->get_alpha();
  double beta_prev = MC->PDF->get_beta();
  double change, change_alpha, change_beta;
  double eta_var;
  int counter = 0;

  // these should be parameters...
  int initial_equi_cycles = 100;
  int equi_cycles = 10;
  int sample_cycles = 1E3;
  int max_steps = 550;
  double tol = 1e-15;

  // first iteration
  Mat R = MC->get_initial_R();
  int N = R.shape0();
  R = MC->equilibriation(R, initial_equi_cycles);
  R = MC->sample_variational_derivatives(R, sample_cycles);

  alpha = alpha_prev - eta*MC->get_grad_alpha();
  beta  = beta_prev - eta*MC->get_grad_beta();

  MC->PDF->update_alpha(alpha);
  MC->PDF->update_beta(beta);

  // absolute change
  change_alpha = alpha - alpha_prev;
  change_beta = beta - beta_prev;
  change = change_alpha*change_alpha + change_beta*change_beta;

  printf("change: %.12lf %.12lf\n", change_alpha, change_beta);
  while (change > tol && counter < max_steps) {

    alpha_prev = alpha;
  	beta_prev = beta;

    // Running Monte-Carlo
    R = MC->equilibriation(R, equi_cycles);
    R = MC->sample_variational_derivatives(R, sample_cycles);

    printf("alpha: %.6lf, beta: %.6lf  E: %.6lf  moves: %.6lf %%\n", alpha, beta, MC->get_energy_mean(), MC->get_accepted_moves_ratio());

    alpha = alpha_prev - eta*MC->get_grad_alpha();
    beta  = beta_prev - eta*MC->get_grad_beta();

    MC->PDF->update_alpha(alpha);
  	MC->PDF->update_beta(beta);

  	change_alpha = alpha - alpha_prev;
	  change_beta = beta - beta_prev;
	  change = change_alpha*change_alpha + change_beta*change_beta;
    counter++;
  }

  if (change < tol) {
    printf("stopped because tolerance reached\n");
  }

  printf("change: %.12lf %.12lf %.12lf\n", change, change_alpha, change_beta);
}
