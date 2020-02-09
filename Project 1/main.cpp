#include <iostream>
#include <chrono>
#include <cmath>
#include "./matpak/Mat.h"
#include "./matpak/tests.h"
#include "./matpak/tools.h"

void tests() {
  auto t0 = std::chrono::high_resolution_clock::now();
  tests_Vec();
  auto t1 = std::chrono::high_resolution_clock::now();
  tests_Mat();
  auto t2 = std::chrono::high_resolution_clock::now();
  tests_MP();
  auto t3 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> vec_elapsed = t1 - t0;
  std::chrono::duration<double> mat_elapsed = t2 - t1;
  std::chrono::duration<double> mp_elapsed = t3 - t2;
  std::chrono::duration<double> tot_elapsed = t3 - t0;
  std::cout << "Elapsed Time (Vectors): " << vec_elapsed.count() << " s\n";
  std::cout << "Elapsed Time (Matrices): " << mat_elapsed.count() << " s\n";
  std::cout << "Elapsed Time (MP): " << mp_elapsed.count() << " s\n";
  std::cout << "Total Elapsed Time: " << tot_elapsed.count() << " s\n";
}

int main() {
  tests();

};
