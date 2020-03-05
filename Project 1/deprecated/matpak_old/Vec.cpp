#include <iostream>
#include <cmath>

#include "tools.h"
#include "Vec.h"

// CONSTRUCTORS
Vec::Vec(const int& a) {
  // Class constructor
  len = a;
  values = new double[len];
}

Vec::Vec(const Vec& u) {
  // Class constructor
  len = u.length();
  values = new double[len];
  for (int i = 0; i < len; i++) {
    values[i] = u.values[i];
  }
}

// DESTRUCTOR

Vec::~Vec() {
  delete [] values;
}

// DATA

int Vec::length() const{
  // Returns the vector's length
  return len;
}

Vec Vec::copy() {
  Vec v(len);
  for (int i = 0; i < len; i++) {
    v.values[i] = values[i];
  }
  return v;
}

// INDEXING

void Vec::set(const double& val, const int& idx) {
  // Sets array value "val" at given index "idx"
  values[idx] = val;
}

double Vec::get(const int& idx) {
  // Returns array value at given index "idx"
  return values[idx];
}

// TESTING

bool Vec::check_dims(const Vec& u, const bool& error){
  if (len != u.length()){
    std::string msg = "Vector length mismatch during operation";
    if (error == true) {
      throw std::invalid_argument(msg);
      return false;
    } else {
      return false;
    }
  } else {
    return true;
  }
}

// VISUALIZATION

std::string Vec::string() {
  std::string str = "[";
  for (int i = 0; i < len; i++){
    str.append(sci_not(values[i], 1, 2));
    if (i < len-1) {
      str.append(", ");
    }
    else {
      str.append("]");
    }
  }
  return str;
}

void Vec::print() {
  // Prints the whole array to the terminal
  std::cout << "[";
  for (int i = 0; i < len; i++){
    std::cout << sci_not(values[i], 1, 2);
    if (i < len-1) {
      std::cout << ", ";
    }
    else {
      std::cout << "]" << std::endl;
    }
  }
}

// ARITHMETIC

Vec Vec::operator+(const Vec& u){
  check_dims(u);
  Vec v(len);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    v.values[i] = values[i] + u.values[i];
  }
  return v;
}

Vec Vec::operator-(const Vec& u){
  check_dims(u);
  Vec v(len);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    v.values[i] = values[i] - u.values[i];
  }
  return v;
}

Vec Vec::operator-(){
  Vec v(len);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    v.values[i] = -values[i];
  }
  return v;
}

Vec Vec::operator*(const double& a) {
  Vec v(len);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    v.values[i] = a*values[i];
  }
  return v;
}

Vec Vec::operator/(const double& a) {
  Vec v(len);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    v.values[i] = values[i]/a;
  }
  return v;
}

// INPLACE OPERATORS

void Vec::operator+=(const Vec& u){
  check_dims(u);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    values[i] = values[i] + u.values[i];
  }
}

void Vec::operator-=(const Vec& u){
  check_dims(u);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    values[i] = values[i] - u.values[i];
  }
}

void Vec::operator*=(const double& a) {
  #pragma omp parallel for
  for (int i = 0; i < len; i++) {
    values[i] = a*values[i];
  }
}

// LINEAR ALGEBRA

double Vec::inner(const Vec& u) {
  check_dims(u);
  double a = 0;
  for (int i = 0; i < len; i++){
    a += values[i]*u.values[i];
  }
  return a;
}

Vec Vec::cross(const Vec& u) {
  check_dims(u);
  if (len != 3) {
    std::string msg = "Cross product must be of 3-element vectors";
    throw std::invalid_argument(msg);
  }
  double a1 = values[0];
  double a2 = values[1];
  double a3 = values[2];
  double u1 = u.values[0];
  double u2 = u.values[1];
  double u3 = u.values[2];
  Vec v(3);
  v.values[0] = a2*u3 - a3*u2;
  v.values[1] = a3*u1 - a1*u3;
  v.values[2] = a1*u2 - a2*u1;
  return v;
}

double Vec::norm() {
  double total = 0;
  for (int i = 0; i < len; i++) {
    total += values[i]*values[i];
  }
  return std::sqrt(total);
}

Vec Vec::unit() {
  Vec v(len);
  double u_norm = norm();
  #pragma omp parallel for
  for (int i = 0; i < len; i++) {
    v.values[i] = values[i]/u_norm;
  }
  return v;
}

// COMPARATORS

bool Vec::operator==(const Vec& u){
  if (check_dims(u, false) == false) {
    return false;
  }
  for (int i = 0; i < len; i++){
    if (values[i] - u.values[i] > tol) {
      return false;
    }
  }
  return true;
}

bool Vec::operator!=(const Vec& u){
  if (check_dims(u, false) == false) {
    return true;
  }
  for (int i = 0; i < len; i++){
    if (values[i] - u.values[i] > tol) {
      return true;
    }
  }
  return false;
}

// ASSIGNMENT OPERATORS

void Vec::operator=(const Vec& u) {
  check_dims(u);
  for (int i = 0; i < len; i++) {
    values[i] = u.values[i];
  }
}
