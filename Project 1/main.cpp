#include <iostream>
#include <chrono>

// #include "./tests/tests_backend.h"
#include "./frontend/adagrad.h"
#include "./wavefunctions/Psi.h"
#include "./wavefunctions/Psi_T.h"
#include "./wavefunctions/Psi_OB.h"

// void run_all_tests() {
//   tests_Psi();
// }

void run_all_parts() {
  bool one_body = true;  // true: interacting mode, false: one-body mode
  if (one_body) {
    Psi_OB PDF(0, 0, 0, 0, 0);
    adagrad(&PDF);
  } else {
    Psi_T PDF(0, 0, 0, 0, 0);
    adagrad(&PDF);
  }
}

Mat change(Mat* R) {
  (*R).set(99, 0, 0);
  return (*R);
}

int main() {
  srand(1337);
  // run_all_tests();
  run_all_parts();

  Mat R(2, 3);
  Mat P(2, 3);

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      R.set(i*3 + j, i, j);
    }
  }
  printf("R:\n");
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      std::cout << R.get(i, j) << " ";
    }
    std::cout << std::endl;
  }
  printf("set P to R. P is now:\n");
  P = R;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      std::cout << P.get(i, j) << " ";
    }
    std::cout << std::endl;
  }
  printf("change value in R:\n");
  R = change(&P);
    printf("R:\n");
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      std::cout << R.get(i, j) << " ";
    }
    std::cout << std::endl;
  }
  printf("P is now:\n");
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      std::cout << P.get(i, j) << " ";
    }
    std::cout << std::endl;
  }
  int a = 1;
  int b = 3;
  double c = (double)a/b;
  std::cout << c << std::endl;

}
