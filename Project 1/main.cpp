#include <iostream>
#include <chrono>

#include "./tests/tests_matpak.h"
#include "./tests/tests_backend.h"

#include "./frontend/project_1b.h"

double step_size = 1;       // Step size during random walk
int steps = 1E5;            // Number of Monte-Carle steps per cycle
int cycles = 1E4;           // Number of Monte-Carlo cycles
int N = 1;                  // Number of Particles
int x_max = 1;              // Maximum Initial Distance From Origin
double alpha = 1;           // Hyperparameter
double beta = 1;            // Hyperparameter
double a = 0;               // Atomic Radius [m]
double omega = 1;           // Harmonic Oscillator Frequency
double omega_z = 1;         // Harmonic Oscillator Z-Frequency

void run_all_tests() {
  tests_Vec();
  tests_Mat();
  tests_Psi();
}

void run_all_parts() {
  monte_carlo_1b(step_size, steps, cycles, N, x_max, alpha, beta, a, omega,
                 omega_z);
}

int main() {
  srand(1337);
  run_all_tests();
  run_all_parts();
}
