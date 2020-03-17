#include <cmath>
#include "Psi.h"
#include "Psi_T.h"
#include "../matpak/Mat.h"

// CONSTRUCTOR
Psi_T::Psi_T() {
  update_alpha(0.5);
  update_beta(2.82843);
  update_a(0.0043);
  update_gamma(2.82843);
}

// CALLING
double Psi_T::psi(Mat* R) {
  double exponent_OB = 0.0;
  double product_C = 1.0;
  double x, dx, r_ij;
  int N = R->shape0();
  int M = R->shape1();
  double* r_i = new double [M];
  for (int i = 0; i < N; i++) {
    // Psi_OB
    x = R->get(i, 0);
    exponent_OB += x*x;
    r_i[0] = x;
    if (M > 1) {
      x = R->get(i, 1);
      exponent_OB += x*x;
      r_i[1] = x;
      if (M == 3) {
        x = R->get(i, 2);
        exponent_OB += beta*x*x;
        r_i[2] = x;
      }
    } 
    // Psi_C
    for (int j = i+1; j < N; j++) {
      r_ij = 0.0;
      for (int l = 0; l < M; l++) {
        dx = r_i[l] - R->get(j, l);
        r_ij += dx*dx;
      }
      r_ij = std::sqrt(r_ij);
      if (r_ij > a) {
        product_C *= f(r_ij);
      } else {
        delete[] r_i;
        return -0.0;
      }
    }
  }
  delete[] r_i;
  return std::exp(-alpha*exponent_OB)*product_C;
}


double Psi_T::probability_density_ratio(Mat* R_new, Mat* R_old, int k) {
  double jastrow = 1.0;
  double phi = 0.0;
  double x, dx, r_ki_old, r_ki_new;
  int M = R_new->shape1();
  double* r_old = new double [M];
  double* r_new = new double [M];
  // calculate:  phi(r^new_k) - phi(r^old_k) 
  x = R_new->get(k, 0);
  phi -= x*x;
  r_new[0] = x;
  x = R_old->get(k, 0);
  phi += x*x;
  r_old[0] = x;
  if (M > 1) {        // y-axis
    x = R_new->get(k, 1);
    phi -= x*x;
    r_new[1] = x;
    x = R_old->get(k, 1);
    phi += x*x;
    r_old[1] = x;
    if (M == 3) {     // z-axis
      x = R_new->get(k, 2);
      phi -= beta*x*x;
      r_new[2] = x;
      x = R_old->get(k, 2);
      phi += beta*x*x;
      r_old[2] = x;
    }
  }
  phi = std::exp(alpha*phi);
  // calculate ratio of jastrow factor
  for (int i = 0; i < R_new->shape0(); i++) {
    if (i == k) {continue;}
    r_ki_old = 0.0; 
    r_ki_new = 0.0;
    for (int l = 0; l < M; l++) {
      x = R_old->get(i, l);
      dx = r_new[l] - x;
      r_ki_new += dx*dx;
      dx = r_old[l] - x;
      r_ki_old += dx*dx;
    }
    r_ki_old = std::sqrt(r_ki_old);
    r_ki_new = std::sqrt(r_ki_new);
    if (r_ki_new <= a) {
      delete[] r_new;
      delete[] r_old;
      return -0.0;
    } 
    if (r_ki_old <= a) {
      delete[] r_new;
      delete[] r_old;
      printf("THIS SHOULD NEVER HAPPEN\n");
      return 1.1;
    }
    jastrow *= f(r_ki_new)/f(r_ki_old);
  }
  delete[] r_new;
  delete[] r_old;
  double ratio = phi*jastrow;
  return ratio*ratio;
}

// CALCULATIONS
double* Psi_T::drift_force(Mat* R, int k) {
  int M = R->shape1();
  double* force = new double [M];
  double* r_k = new double [M];
  double* diff_r_ki = new double [M];
  double x, dx, r_ki, up_ki;
  // grad phi
  x = R->get(k, 0);
  force[0] = minus_four_alpha*x;
  r_k[0] = x;
  if (M > 1) {
    x = R->get(k, 1);
    force[1] = minus_four_alpha*x;
    r_k[1] = x;
    if (M == 3) {
      x = R->get(k, 2);
      force[2] = minus_four_alpha_beta*x;
      r_k[2] = x;
    }
  }
  // grad Psi_T
  for (int i = 0; i < R->shape0(); i++) {
    if (i == k) {continue;}
    r_ki = 0.0;
    for (int l = 0; l < M; l++) {
      dx = r_k[l] - R->get(i, l);
      r_ki += dx*dx;
      diff_r_ki[l] = dx;
    }
    r_ki = std::sqrt(r_ki);
    up_ki = 2*u_prime(r_ki);
    for (int l = 0; l < M; l++) {
      force[l] += diff_r_ki[l]/r_ki*up_ki;
    }
  }
  delete[] r_k;
  delete[] diff_r_ki;
  return force;
}

double Psi_T::energy(Mat* R) {
  int N = R->shape0();
  int M = R->shape1();

  double K;                 // Kinetic Energy
  double V = 0.0;           // Potential Energy (External)
  double laplace_phi = 0.0;
  double laplace_Psi_C = 0.0;
  double grad_phi_grad_Psi_C = 0.0;
  double* grad_phi = new double [M];
  double* grad_Psi_C = new double [M];
  double x, xx, dx; 
  double* r_k = new double [M];
  double r_kj, r_ki, diff_r_kj_r_ki;
  double* diff_r_kj = new double [M];
  double up_kj, up_ki;        // u-prime

  for (int k = 0; k < N; k++) {
    x = R->get(k, 0);
    xx = x*x;
    r_k[0] = x;
    grad_phi[0] = minus_two_alpha*x;
    grad_Psi_C[0] = 0.0;
    laplace_phi += xx;
    V += xx;
    if (M > 1) {// y-axis
      x = R->get(k, 1);
      xx = x*x;
      r_k[1] = x;
      grad_phi[1] = minus_two_alpha*x;
      grad_Psi_C[1] = 0.0;
      laplace_phi += xx;
      V += xx;
      if (M == 3) {// z-axis
        x = R->get(k, 2);
        xx = x*x;
        r_k[2] = x;
        grad_phi[2] = minus_two_alpha_beta*x;
        grad_Psi_C[2] = 0.0;
        laplace_phi += beta_squared*xx;
        V += gamma_squared*xx;
      }
    }
    for (int j = 0; j < N; j++) {
      if (j == k) {continue;}
      r_kj = 0.0;
      for (int l = 0; l < M; l++) {
        dx = r_k[l] - R->get(j, l);
        r_kj += dx*dx;
        diff_r_kj[l] = dx;
      }
      r_kj = std::sqrt(r_kj);
      if (r_kj <= a) {printf("cry\n");}
      up_kj = u_prime(r_kj);
      for (int l = 0; l < M; l++) {
        grad_Psi_C[l] += diff_r_kj[l]/r_kj*up_kj;
      }
      laplace_Psi_C += u_double_prime(r_kj, up_kj) + up_kj*(2.0/r_kj);
      for (int i = 0; i < N; i++) {
        if (i == k) {continue;}
        r_ki = 0.0;
        diff_r_kj_r_ki = 0.0;
        for (int l = 0; l < M; l++) {
          dx = r_k[l] - R->get(i, l);
          r_ki += dx*dx;
          diff_r_kj_r_ki += diff_r_kj[l]*dx;
        }
        r_ki = std::sqrt(r_ki);
        up_ki = u_prime(r_ki);
        if (r_ki <= a) { printf("cry\n");}
        laplace_Psi_C += diff_r_kj_r_ki/(r_kj*r_ki)*up_ki*up_kj;
      }// END LOOP OVER i
    } // END LOOP OVER j
    for (int l = 0; l < M; l++) {
      grad_phi_grad_Psi_C += grad_phi[l]*grad_Psi_C[l];
    }
  }  // END LOOP OVER k
  delete[] grad_phi;
  delete[] grad_Psi_C;
  delete[] r_k;
  delete[] diff_r_kj;
  double M_beta = (M == 3) ? 2.0 + beta : M;
  laplace_phi = 4*alpha_squared*laplace_phi - 2*N*alpha*M_beta;
  K = laplace_phi + laplace_Psi_C + 2*grad_phi_grad_Psi_C;
  return 0.5*(-K + V);
}

double Psi_T::f(double r_ij) {
    return 1.0 - a/r_ij;
}

double Psi_T::u_prime(double r_ij) {
    return a/(r_ij*(r_ij - a));
}

double Psi_T::u_double_prime(double r_ij, double u_prime_ij) {
    return -u_prime_ij*(1.0/(r_ij - a) + 1.0/r_ij);
}

std::string Psi_T::name() {
  std::string name = "repulsive";
  return name;
}

bool Psi_T::interaction() {
  return true;
}

