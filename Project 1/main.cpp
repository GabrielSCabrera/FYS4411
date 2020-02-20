#include <iostream>
#include <chrono>

#include "./tests/tests_matpak.h"
#include "./tests/tests_backend.h"

#include "./frontend/project_1b.h"

double step_size = 1E-1;     // Step size during random walk
int steps = 1E4;            // Number of Monte-Carle steps per cycle
int cycles = 1E3;           // Number of Monte-Carlo cycles
int N = 20;                  // Number of Particles
int x_max = 10;              // Maximum Initial Distance From Origin
double alpha = 1E-4;        // Hyperparameter
double beta = 1E-4;         // Hyperparameter
double a = 2.9e-10;         // Atomic Radius [m]
double omega = 100;         // Harmonic Oscillator Frequency
double omega_z = 100;       // Harmonic Oscillator Z-Frequency
double mass = 1.42e-25;     // Atomic Mass [kg]

void run_all_tests() {
  tests_Vec();
  tests_Mat();
  tests_Psi();
}

void run_all_parts() {
  monte_carlo_1b(step_size, steps, cycles, N, x_max, alpha, beta, a, omega,
                 omega_z, mass);
}

int main() {
  srand(1337);
  run_all_tests();
  run_all_parts();
}
