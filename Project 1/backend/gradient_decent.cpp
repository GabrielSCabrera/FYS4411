
/*
Should we have to separate eta for beta and alpha?
*/
void gradient_decent(Monte_Carlo MC, double eta) {
  double alpha, beta;
  double alpha_prev = 0.8;
  double beta_prev = 0.5;
  double change, change_alpha, change_beta;
  int counter = 0;

  // should be the other way around when we are done debugging
  MC->PDF->update_alpha(alpha_prev);
  MC->PDF->update_beta(beta_prev);

  // these should be parameters...
  int initial_equi_cycles = 1E4;
  int equi_cycles = 1E2;
  int sample_cycles = 104;
  int max_steps = 100;
  double tol = 1e-6;

  // first iteration
  Mat R = MC.get_initial_R();
  R = MC.equilibriation(R, initial_equi_cycles);
  R = MC.sample_variational_derivatives(R, sample_cycles);

  alpha = alpha_prev - eta*MC.get_grad_alpha();
  beta  = beta_prev - eta*MC.get_grad_alpha();

  MC->PDF->update_alpha(alpha);
  MC->PDF->update_beta(beta);

  // absolute change
  change_alpha = alpha - alpha_prev;
  change_beta = beta - beta_prev;
  change = change_alpha*change_alpha + change_beta*change_beta;

  while (change < tol && counter < max_steps) {
  	alpha_prev = alpha;
  	beta_prev = beta;
    // Running Monte-Carlo
    R = MC.equilibriation(R, equi_cycles);
    R = MC.sample_variational_derivatives(R, sample_cycles);

    /* Displaying Stats
        THESE GETTERS DO NOT EXIST YET
    */
    std::cout << "MC-Cycle – alpha = " << alpha 
    		<< ", beta = " << beta
            << ", E = " << MC.get_E()/(N*3) << 
            ", var = " << MC.get_variance() << std::endl;

    /*
    THESE GETTERS DO NOT EXIST YET
    */
    alpha = alpha_prev - eta*MC.get_grad_alpha();
    beta  = beta_prev - eta*MC.get_grad_alpha();

    MC->PDF->update_alpha(alpha);
  	MC->PDF->update_beta(beta);

  	change_alpha = alpha - alpha_prev;
	change_beta = beta - beta_prev;
	change = change_alpha*change_alpha + change_beta*change_beta;
}