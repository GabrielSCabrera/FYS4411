#include <cmath>

#include "Psi.h"
#include "../matpak/Mat.h"

// CONSTRUCTOR

Psi::Psi(double alpha, double beta, double a, double omega, double omega_z,
         double mass) {
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
  int N = P.shape().get(0);
  return Psi_ob(P, N)*Psi_c(P, N);
}

// CALCULATIONS

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

      if (diff <= 0) {
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

  int N = P.shape().get(0);

  double E = 0;
  double V = 0;

  double x1; double y1; double z1;
  double x2; double y2; double z2;
  double x3; double y3; double z3;

  double r12; double r13;

  // u-prime
  double up12; double up13;

  double dx12; double dy12; double dz12;
  double dx13; double dy13; double dz13;

  double* gp = new double [3] {0, 0, 0};
  double* term_2 = new double [3] {0, 0, 0};
  double* term_4 = new double [3] {0, 0, 0};

  for (int k = 0; k < N; k++) {

    x1 = P.get(k,0); y1 = P.get(k,1); z1 = P.get(k,2);

    // START Term 1
    E += laplace_phi(x1, y1, z1);
    // END Term 1

    gp = grad_phi(x1, y1, z1);

    term_2[0] = 0; term_2[1] = 0; term_2[2] = 0;
    term_4[0] = 0; term_4[1] = 0; term_4[2] = 0;

    for (int i = 0; i < N; i++) {
      if (i == k) {
        continue;
      }

      x2 = P.get(i,0); y2 = P.get(i,1); z2 = P.get(i,2);

      dx12 = x2 - x1; dy12 = y2 - y1; dz12 = z2 - z1;

      r12 = std::sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12);

      up12 = u_prime(r12);

      // START Term 2
      term_2[0] += (dx12/r12)*up12;
      term_2[1] += (dy12/r12)*up12;
      term_2[2] += (dz12/r12)*up12;
      // END Term 2

      // START Term 3
      for (int j = 0; j < N; j++) {
        if (j == k) {
          continue;
        }
        x3 = P.get(j,0); y3 = P.get(j,1); z3 = P.get(j,2);

        dx13 = x3 - x1; dy13 = y3 - y1; dz13 = z3 - z1;

        r13 = std::sqrt(dx13*dx13 + dy13*dy13 + dz13*dz13);

        up13 = u_prime(r13);

        term_4[0] += (dx12*dx13)/(r12*r13) * up12 * up13;
        term_4[1] += (dy12*dy13)/(r12*r13) * up12 * up13;
        term_4[2] += (dz12*dz13)/(r12*r13) * up12 * up13;
      }
      // END Term 3

      // START Term 4
      E += u_double_prime(r12)*up12 + (2/r12)*up12;
      // END Term 4

    }
    // MULTIPLY WITH AND SUM OVER TERM 2!!!

    // START Term 2
    term_2[0] *= 2*gp[0];
    term_2[1] *= 2*gp[1];
    term_2[2] *= 2*gp[2];
    E += term_2[0] + term_2[1]+ term_2[2];
    // END Term 2

    // START Term 4
    E += term_4[0] + term_4[1]+ term_4[2];
    // END Term 4

    // START V_ext
    V += V_ext(x1, y1, z1);
    // END V_ext

    // // START V_int
    // for (int i = k+1; i < N; i++) {
    //   if
    // }
    // // END V_int

  }
  delete[] gp;
  delete[] term_2;
  delete[] term_4;
  return -0.5*1.0545718E-34*E/this->mass + V;
}

double Psi::V_ext(double x, double y, double z) {
  return 0.5*this->mass*(this->omega*this->omega*(x*x + y*y)
                      + this->omega_z*this->omega_z*z*z);
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
  return -2*this->alpha*(3-2*this->alpha*(x*x+y*y+this->beta*this->beta*z*z));
}

double Psi::u_prime(double r_jk) {
  if (r_jk <= this-> a) {
    return 0;
  } else {
    return this->a/(r_jk*(r_jk + this->a));
  }
}

double Psi::u_double_prime(double r_jk) {
  /*
    !! Must be multiplied with u_prime(r_jk) !!
  */
  if (r_jk <= this-> a) {
    return 0;
  } else {
    return -1*(1/(r_jk - this->a) + 1/r_jk);
  }
}
