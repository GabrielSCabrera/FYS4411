#include <iostream>
#include <chrono>

#include "./tests/tests_matpak.h"
#include "./tests/tests_backend.h"

void run_all_tests() {
  tests_Vec();
  tests_Mat();
  tests_Psi();
}

int main() {
  run_all_tests();
}
