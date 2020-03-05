#include <cmath>

#include "Psi.h"
#include "../matpak/Mat.h"

// CONSTRUCTOR

Psi::Psi(double alpha, double beta, double a, double omega, double omega_z) {
  this->alpha = alpha;
  this->beta = beta;
  this->a = a;
  this->omega = omega;
  this->omega_z = omega_z;
  this->mass = mass;
}

// DESTRUCTOR

Psi::~Psi() {

}

// UPDATERS

void Psi::update_alpha(double alpha) {
  this->alpha = alpha;
}

void Psi::update_beta(double beta) {
  this->beta = beta;
}

void Psi::update_a(double a) {
  this->a = a;
}

void Psi::update_omega(double omega) {
  this->omega = omega;
}

void Psi::update_omega_z(double omega_z) {
  this->omega_z = omega_z;
}

void Psi::update_mass(double mass) {
  this->mass = mass;
}

// ATTRIBUTE EXTRACTION

double Psi::get_alpha() {
  return this->alpha;
}

double Psi::get_beta() {
  return this->beta;
}

double Psi::get_a() {
  return this->a;
}

// CALLING

double Psi::operator()(Mat P) {
  int N = P.shape0();
  if (N == 1) {
    return Psi_ob(P, N);
  } else {
    return Psi_ob(P, N)*Psi_c(P, N);
  }
}

// CALCULATIONS

double Psi::V_ext(double x, double y, double z) {
  return 0.5*(this->omega*this->omega*(x*x + y*y)
                      + this->omega_z*this->omega_z*z*z);
}

double* Psi::drift(double x, double y, double z) {
  double* force = new double [3];
  force[0] = -4*this->alpha*x;
  force[1] = -4*this->alpha*y;
  force[2] = -4*this->alpha*z;
  return force;
}

double Psi::greens_ratio(double x0, double y0, double z0,
                         double x1, double y1, double z1,
                         double K) {
  // K = D*dt

  double* terms = new double[3];

  double* F1 = drift(x1, y1, z1);
  terms[0] = x0-x1-K*F1[0]; terms[1] = y0-y1-K*F1[1]; terms[2] = z0-z1-K*F1[2];
  double exponent = terms[0]*terms[0] + terms[1]*terms[1] + terms[2]*terms[2];
  delete [] F1;

  double* F2 = drift(x0, y0, z0);
  terms[0] = x1-x0-K*F2[0]; terms[1] = y1-y0-K*F2[1]; terms[2] = z1-z0-K*F2[2];
  exponent -= terms[0]*terms[0] + terms[1]*terms[1] + terms[2]*terms[2];

  delete [] F2;
  delete [] terms;
  return std::exp(exponent/(4*K));
}

double Psi::Psi_ob(Mat P, int N) {
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

double Psi::Psi_c(Mat P, int N) {
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

double Psi::energy(Mat P) {

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

double Psi::phi(double x, double y, double z) {
  return std::exp(-this->alpha*(x*x + y*y + this->beta*z*z));
}

double* Psi::grad_phi(double x, double y, double z) {
  /*
    !! Must be multiplied with phi(x, y, z) !!
  */
  double* grad = new double [3];
  grad[0] = -2*this->alpha*x;
  grad[1] = -2*this->alpha*y;
  grad[2] = -2*this->alpha*z*this->beta;
  return grad;
}

double Psi::laplace_phi(double x, double y, double z) {
  /*
    !! Must be multiplied with phi(x, y, z) !!
  */
  return 2*this->alpha*(2*this->alpha*(x*x+y*y+this->beta*this->beta*z*z) - 2 - this->beta);
}

double Psi::u_prime(double r_kj) {
  if (r_kj <= this-> a) {
    return 0;
  } else {
    return this->a/(r_kj*(r_kj - this->a));
  }
}

double Psi::u_double_prime(double r_kj) {
  /*
    !! Must be multiplied with u_prime(r_kj) !!
  */
  if (r_kj <= this-> a) {
    return 0;
  } else {
    return -(1/(r_kj - this->a) + 1/r_kj);
  }
}

double Psi::grad_alpha(Mat P) {

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
      dx_kj = x_k - x_j; dy_kj = y_k - y_j; dz_kj = z_k - z_j;
      r_kj = std::sqrt(dx_kj*dx_kj + dy_kj*dy_kj + dz_kj*dz_kj);

      term_2 += (1/r_kj)*(x_k*dx_kj + y_k*dy_kj + this->beta*z_k*dz_kj)*u_prime(r_kj);
    }
  }
  return 8*this->alpha*term_1 - 2*term_2 - 6*N;
}

double Psi::grad_beta(Mat P) {

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
      dx_kj = x_k - x_j; dy_kj = y_k - y_j; dz_kj = z_k - z_j;
      r_kj = std::sqrt(dx_kj*dx_kj + dy_kj*dy_kj + dz_kj*dz_kj);

      term_2 += (1/r_kj)*(z_k*dz_kj)*u_prime(r_kj);
    }
  }
  return 8*this->beta*this->alpha*this->alpha*term_1 - 2*this->alpha*term_2;
}

double Psi::grad_alpha_alpha(Mat P) {

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

double Psi::grad_beta_beta(Mat P) {

  double z_k;

  int N = P.shape0();

  double term_1 = 0;

  for (int k = 0; k < N; k++) {

    z_k = P.get(k,2);

    term_1 += z_k*z_k;

  }
  return 8*this->alpha*this->alpha*term_1;
}

double Psi::grad_beta_alpha(Mat P) {

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
      dx_kj = x_k - x_j; dy_kj = y_k - y_j; dz_kj = z_k - z_j;
      r_kj = std::sqrt(dx_kj*dx_kj + dy_kj*dy_kj + dz_kj*dz_kj);

      term_2 += (1/r_kj)*(z_k*dz_kj)*u_prime(r_kj);
    }
  }
  return 16*this->beta*this->alpha*term_1 - 2*term_2;
}
