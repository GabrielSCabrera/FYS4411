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
    Psi_OB PDF(0, 0, 0, 0);
    adagrad(&PDF);
  } else {
    Psi_T PDF(0, 0, 0, 0);
    adagrad(&PDF);
  }
}

int main() {
  srand(1337);
  // run_all_tests();
  run_all_parts();
}
