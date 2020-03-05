#include <iostream>
#include <chrono>

#include "./tests/tests_backend.h"
#include "./backend/monte_carlo.h"

void run_all_tests() {
  tests_Psi();
}


void run_all_parts() {
  run_monte_carlo();
}

int main() {
  srand(1337);
  run_all_tests();
  run_all_parts();
}
