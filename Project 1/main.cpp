#include <iostream>
#include <chrono>

#include "./tests/tests_backend.h"
#include "./frontend/adagrad.h"

void run_all_tests() {
  tests_Psi();
}


void run_all_parts() {
  adagrad();
}

int main() {
  srand(1337);
  run_all_tests();
  run_all_parts();
}
