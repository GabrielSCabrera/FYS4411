#include <cmath>

#include "Psi.h"
#include "Psi_T.h"
#include "../matpak/Mat.h"

// CALLING
double Psi_T::operator()(Mat R) {
  double exponent_OB = 0.0;
  double product_C = 1.0;
  double x, dx, r_ij;
  int N = R.shape0();
  int M = R.shape1();
  int last_index = M - 1;
  double* r_i = new double [M];

  for (int i = 0; i < N; i++) {
    // Psi_OB
    for (int j = 0; j < last_index; j++) {
      x = R.get(i, j);
      exponent_OB += x*x;
      r_i[j] = x;
    }
    x = R.get(i, last_index);
    exponent_OB += beta*x*x;
    r_i[last_index] = x;
    // Psi_C
    for (int j = i+1; j < N; j++) {
      r_ij = 0.0;
      for (int k = 0; k < M; k++) {
        dx = r_i[k] - R.get(j, k);
        r_ij += dx*dx;
      }
      r_ij = std::sqrt(r_ij);
      if (r_ij > a) {
        product_C *= 1 - a/r_ij;
      } else {
        delete[] r_i;
        return 0;
      }
    }
  }
  delete[] r_i;
  return std::exp(-alpha*exponent_OB)*product_C;
}

// CALCULATIONS
double* Psi_T::drift_force(Mat R, int k) {
  int M = R.shape1();
  int last_index = M - 1;
  double* force = new double [M];
  double* r_k = new double [M];
  double* diff_r_ki = new double [M];
  double x, dx, r_ki, up_ki;
  // grad phi
  for (int l = 0; l < last_index; l++) {
    x = R.get(k, l);
    force[l] = -4*alpha*x;
    r_k[l] = x;
  }
  x = R.get(k, last_index);
  force[last_index] = -4*alpha*beta*x;
  r_k[last_index] = x;
  // grad Psi_T
  for (int i = 0; i < R.shape0(); i++) {
    if (i == k) {continue;}
    r_ki = 0.0;
    for (int l = 0; l < M; l++) {
      dx = r_k[l] - R.get(i, l);
      r_ki += dx*dx;
      diff_r_ki[l] = dx;
    }
    r_ki = std::sqrt(r_ki);
    up_ki = u_prime(r_ki);
    for (int l = 0; l < M; l++) {
      force[l] += diff_r_ki[l]/r_ki*up_ki;
    }
  }
  delete[] r_k;
  delete[] diff_r_ki;
  return force;
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

double Psi_T::energy(Mat R) {
  int N = R.shape0();
  int M = R.shape1();
  int last_index = M - 1;
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
    for (int l = 0; l < last_index; l++) {
      x = R.get(k, l);
      xx = x*x;

      r_k[l] = x;
      grad_phi[l] = -2*alpha*x;
      grad_Psi_C[l] = 0.0;

      laplace_phi += xx;
      V += xx;
    }
    x = R.get(k, last_index);
    xx = x*x;

    r_k[last_index] = x;
    grad_phi[last_index] = -2*alpha*beta*x;
    grad_Psi_C[last_index] = 0.0;

    laplace_phi += beta_squared*xx;
    V += gamma_squared*xx;

    for (int j = 0; j < N; j++) {
      if (j == k) {continue;}
      r_kj = 0.0;
      for (int l = 0; l < M; l++) {
        dx = r_k[l] - R.get(j, l);
        r_kj += dx*dx;
        diff_r_kj[l] = dx;
      }
      r_kj = std::sqrt(r_kj);
      //if (r_kj < a) {printf("cry\n");}
      up_kj = u_prime(r_kj);
      for (int l = 0; l < M; l++) {
        grad_Psi_C[l] += diff_r_kj[l]/r_kj*up_kj;
      }
      laplace_Psi_C += up_kj*(last_index/r_kj - (1.0/(r_kj - a) + 1.0/r_kj));
      for (int i = 0; i < N; i++) {
        if (i == k) {continue;}
        r_ki = 0.0;
        diff_r_kj_r_ki = 0.0;
        for (int l = 0; l < M; l++) {
          dx = r_k[l] - R.get(i, l);
          r_ki += dx*dx;
          diff_r_kj_r_ki += diff_r_kj[l]*dx;
        }
        r_ki = std::sqrt(r_ki);
        up_ki = u_prime(r_ki);
        //if (r_ki < a) { printf("cry\n");}
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
  laplace_phi = 4*alpha_squared*laplace_phi - 2*N*alpha*(last_index + beta);
  K = laplace_phi + laplace_Psi_C + 2*grad_phi_grad_Psi_C;
  return 0.5*(-K + V);

}

double Psi_T::u_prime(double r_kj) {
    return a/(r_kj*(r_kj - a));
}

double Psi_T::u_double_prime(double r_kj) {
    return -(1.0/(r_kj - a) + 1.0/r_kj)*u_prime(r_kj);
}



