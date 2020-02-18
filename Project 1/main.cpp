#include <iostream>
#include <chrono>

#include "./tests/tests_matpak.h"
#include "./tests/tests_backend.h"

#include "./frontend/project_1b.h"

void run_all_tests() {
  tests_Vec();
  tests_Mat();
  tests_Psi();
}

void run_all_parts() {
  project_1b_main();
}

int main() {
  srand(1337);
  run_all_tests();
  run_all_parts();
}
