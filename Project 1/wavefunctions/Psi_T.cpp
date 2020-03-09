#include <cmath>

#include "Psi.h"
#include "Psi_T.h"
#include "../matpak/Mat.h"

// CALLING

double Psi_T::operator()(Mat P) {
  int N = P.shape0();
  return Psi_ob(P, N)*Psi_c(P, N);
}

// CALCULATIONS

double* Psi_T::drift(double x, double y, double z) {
  double* force = new double [3];
  force[0] = -4*this->alpha*x;
  force[1] = -4*this->alpha*y;
  force[2] = -4*this->alpha*z;
  return force;
}

double Psi_T::Psi_ob(Mat P, int N) {
  /*
    P – Array of type 'Mat' of shape (N,3)
    N – Number of particles
  */
  double exponent = 0;
  double x; double y; double z;
  for (int i = 0; i < N; i++) {
    x = P.get(i,0);
    y = P.get(i,1);
    z = P.get(i,2);
    exponent += x*x + y*y + this->beta*z*z;
  }
  exponent = std::exp(-this->alpha*exponent);
  return exponent;
}

double Psi_T::Psi_c(Mat P, int N) {
  /*
    P – Array of type 'Mat' of shape (N,3)
    N – Number of particles
  */
  double product = 1;
  double dx; double dy; double dz;
  double diff;
  for (int i = 0; i < N; i++) {
    for (int j = i+1; j < N; j++) {

      dx = P.get(i,0) - P.get(j,0);
      dy = P.get(i,1) - P.get(j,1);
      dz = P.get(i,2) - P.get(j,2);

      diff = std::sqrt(dx*dx + dy*dy + dz*dz);

      if (diff <= this->a) {
        product = 0;
        goto stop;
      } else {
        product *= 1 - this->a/diff;
      }
    }
  }
  stop:
  return product;
}

double Psi_T::energy(Mat P) {

  int N = P.shape0();

  // Kinetic Energy
  double E = 0;

  // Potential Energy (External)
  double V = 0;

  double x_k; double y_k; double z_k;
  double x_j; double y_j; double z_j;
  double x_i; double y_i; double z_i;

  double r_kj; double r_ki;

  // u-prime
  double up_kj; double up_ki;

  double dx_kj; double dy_kj; double dz_kj;
  double dx_ki; double dy_ki; double dz_ki;

  double* term_2 = new double [3] {0, 0, 0};

  for (int k = 0; k < N; k++) {

    x_k = P.get(k,0); y_k = P.get(k,1); z_k = P.get(k,2);

    // START Term 1
    E += laplace_phi(x_k, y_k, z_k);
    // END Term 1

    term_2[0] = 0; term_2[1] = 0; term_2[2] = 0;

    for (int j = 0; j < N; j++) {

      if (j == k) {
        continue;
      }

      x_j = P.get(j,0); y_j = P.get(j,1); z_j = P.get(j,2);
      dx_kj = x_k - x_j; dy_kj = y_k - y_j; dz_kj = z_k - z_j;
      r_kj = std::sqrt(dx_kj*dx_kj + dy_kj*dy_kj + dz_kj*dz_kj);
      up_kj = u_prime(r_kj);

      // START Term 2
      term_2[0] += dx_kj*up_kj/r_kj;
      term_2[1] += dy_kj*up_kj/r_kj;
      term_2[2] += dz_kj*up_kj/r_kj;
      // END Term 2

      // START Term 5
      E += (u_double_prime(r_kj) + 2/r_kj)*up_kj;
      // END Term 5

      // START Term 3
      for (int i = 0; i < N; i++) {

        if (i == k) {
          continue;
        }

        x_i = P.get(i,0); y_i = P.get(i,1); z_i = P.get(i,2);
        dx_ki = x_k - x_i; dy_ki = y_k - y_i; dz_ki = z_k - z_i;
        r_ki = std::sqrt(dx_ki*dx_ki + dy_ki*dy_ki + dz_ki*dz_ki);
        up_ki = u_prime(r_ki);

        // START Term 4
        E += (dx_kj*dx_ki + dy_kj*dy_ki + dz_kj*dz_ki)*up_ki*up_kj/(r_kj*r_ki);
        // END Term 4

      } // END LOOP OVER i
    } // END LOOP OVER j

    // START Term 2
    double* gp = grad_phi(x_k, y_k, z_k);
    E += 2*(term_2[0]*gp[0] + term_2[1]*gp[1] + term_2[2]*gp[2]);
    // END Term 2

    // START V_ext
    V += V_ext(x_k, y_k, z_k);
    // END V_ext

    delete[] gp;
  } // END LOOP OVER k

  delete[] term_2;

  return -0.5*E + V;

}

double Psi_T::grad_alpha(Mat P) {

  double beta_sq = this->beta*this->beta;

  double x_k; double y_k; double z_k;
  double x_j; double y_j; double z_j;

  double r_kj;

  double dx_kj; double dy_kj; double dz_kj;

  int N = P.shape0();

  double term_1 = 0;
  double term_2 = 0;

  for (int k = 0; k < N; k++) {

    x_k = P.get(k,0); y_k = P.get(k,1); z_k = P.get(k,2);

    term_1 += x_k*x_k + y_k*y_k + beta_sq*z_k*z_k;

    for (int j = 0; j < N; j++) {

      if (j == k) {
        continue;
      }

      x_j = P.get(j,0); y_j = P.get(j,1); z_j = P.get(j,2);
      dx_kj = x_j - x_k; dy_kj = y_j - y_k; dz_kj = z_j - z_k;
      r_kj = std::sqrt(dx_kj*dx_kj + dy_kj*dy_kj + dz_kj*dz_kj);

      term_2 += (x_k*dx_kj + y_k*dy_kj + this->beta*z_k*dz_kj)*u_prime(r_kj)/r_kj;
    }
  }
  return 8*this->alpha*term_1 - 2*term_2 - 4*N - 2*this->beta*N;
}

double Psi_T::grad_beta(Mat P) {

  double x_k; double y_k; double z_k;
  double x_j; double y_j; double z_j;

  double r_kj;

  double dx_kj; double dy_kj; double dz_kj;

  int N = P.shape0();

  double term_1 = 0;
  double term_2 = 0;

  for (int k = 0; k < N; k++) {

    x_k = P.get(k,0); y_k = P.get(k,1); z_k = P.get(k,2);

    term_1 += z_k*z_k;

    for (int j = 0; j < N; j++) {

      if (j == k) {
        continue;
      }

      x_j = P.get(j,0); y_j = P.get(j,1); z_j = P.get(j,2);
      dx_kj = x_j - x_k; dy_kj = y_j - y_k; dz_kj = z_j - z_k;
      r_kj = std::sqrt(dx_kj*dx_kj + dy_kj*dy_kj + dz_kj*dz_kj);

      term_2 += z_k*dz_kj*u_prime(r_kj)/r_kj;
    }
  }
  return 8*this->beta*this->alpha*this->alpha*term_1 - 2*this->alpha*(term_2 + N);
}

double Psi_T::grad_alpha_alpha(Mat P) {

  double beta_sq = this->beta*this->beta;

  double x_k; double y_k; double z_k;

  int N = P.shape0();

  double term_1 = 0;

  for (int k = 0; k < N; k++) {

    x_k = P.get(k,0); y_k = P.get(k,1); z_k = P.get(k,2);

    term_1 += x_k*x_k + y_k*y_k + beta_sq*z_k*z_k;

  }
  return 8*term_1;
}

double Psi_T::grad_beta_beta(Mat P) {

  double z_k;

  int N = P.shape0();

  double term_1 = 0;

  for (int k = 0; k < N; k++) {

    z_k = P.get(k,2);

    term_1 += z_k*z_k;

  }
  return 8*this->alpha*this->alpha*term_1;
}
