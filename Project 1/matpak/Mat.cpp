#include <iostream>
#include <cmath>

#include "Mat.h"

// CONSTRUCTOR

Mat::Mat(const int& d1, const int& d2) {

  if (d1 == 0 || d2 == 0) {
    std::string msg = "Matrix cannot have a side-length zero!";
    throw std::invalid_argument(msg);
  }

  size = new int[2] {d1, d2};
  dims = 2;
  len = d1*d2;

  idx_increments = new int[2] {d2, 1};

  values = new double[len];
  for (int i = 0; i < len; i++) {
    values[i] = 0;
  }
}

Mat::Mat(const Mat& u) {
  size = new int[2] {u.size[0], u.size[1]};

  len = u.len;
  dims = u.dims;

  values = new double[len];

  for (int i = 0; i < len; i++) {
    values[i] = u.values[i];
  }

  idx_increments = new int[2] {u.idx_increments[0], u.idx_increments[1]};
}

// DESTRUCTOR

Mat::~Mat() {
  delete[] size;
  delete[] idx_increments;
  delete[] values;
}

// DATA

int Mat::length() {
  // Returns the Matrix's length
  return len;
}

int Mat::dimension() {
  return dims;
}

int* Mat::shape() {
  return size;
}

// INDEXING

int Mat::locate(const int& idx1, const int& idx2) {
  // Converts 2-dimensional indices to 1-D indices
  /*
  [[1,2,3]         [[(0,0) (0,1) (0,2)]
   [4,5,6]]         [(1,0) (1,1) (1,2)]]


   [1,2,3,4,5,6]    [0,1,2,3,4,5]
  */

  return idx1*idx_increments[0] + idx2*idx_increments[1];
}

void Mat::set(const double& val, const int& idx1, const int& idx2) {
  values[locate(idx1, idx2)] = val;
}

double Mat::get(const int& idx1, const int& idx2) {
  // Returns array value at given indices
  return values[locate(idx1, idx2)];
}

void Mat::set_raw(const double& val, const int& idx) {
  // Set an element from the "values" array directly
  values[idx] = val;
}

// ASSIGNMENT OPERATORS

void Mat::operator= (const std::string& s) {
  /* Parses an input string of the form:
    "[[1,2,3], [4,5,6], [7,8,9]]"

    As a matrix, this looks like:
              1 2 3
              4 5 6
              7 8 9

    The shape, dimension, length, and individual
    elements are extracted from the string.

    The user is responsible for making a valid
    shape for the Matrix, irregular rows might
    not be noticed by the parser, for example.
  */
  int N_s = 0;  // Number of Matrix elements
  int size_s[6] = {0, 0, 0, 0, 0, 0}; // Matrix dimensions
  int dims_s = 0; // Number of Matrix dimensions (temporary)
  int dims_max = 0; // Number of Matrix dimensions
  std::string c;

  // Iterating through the string, gathering data
  for (unsigned int i = 0; i < s.length(); i++) {

    c = s.at(i);
    if (c == "[") {
      dims_s++;
    } else if (c == "]") {
      dims_s--;
    } else if (c == ",") {
      if (dims_s > dims_max) {
        dims_max = dims_s;
      }
      size_s[dims_s-1]++;
      N_s++;
    }
    // Checking that the array is of valid form
    if (dims_s < 0) {
      std::string msg = "Invalid Mat assignment string;\n\t\t";
      msg += "Extra closing bracket.";
      throw std::domain_error(msg);
    } else if (dims_s > 5) {
      std::string msg = "Invalid Mat assignment string;\n\t\t";
      msg += "Too many dimensions – 6 is the limit.";
      throw std::domain_error(msg);
    }
  }

  if (dims_s > 0) {
    std::string msg = "Invalid Mat assignment string;\n\t\t";
    msg += "Extra opening bracket.";
    throw std::domain_error(msg);
  }

  // Taking the last comma into consideration
  if (N_s > 0) {
    N_s++;
  }

  // Scaling each dimension's size

  for (int i = dims_max-1; i > 0; i--) {
    size_s[i] /= size_s[i-1];
  }

  // Reiterating through "s" to extract each element
  std::string elements_s[N_s];
  int e = 0;
  bool increment = true;
  for (unsigned int i = 0; i < s.length(); i++) {
    c = s.at(i);
    if (c == "]" || c == ",") {
      if (increment == true) {
        e++;
      }
      increment = false;
    } else if (c != "[" && c != " " && c != "\n") {
      increment = true;
      elements_s[e].append(c);
    }
  }

  /* Converting the string elements to doubles
  and replacing our array with the new values */
  for (int i = 0; i < N_s; i++) {
    values[i] = std::stof(elements_s[i]);
  }
  }

void Mat::operator= (const Mat& u) {
  check_size(u);
  #pragma omp parallel for
  for (int i = 0; i < len; i++) {
    values[i] = u.values[i];
  }
}

// TESTING

bool Mat::check_size(const Mat& u, const bool& error) {
  for (int i = 0; i < 2; i++) {
    if (size[i] != u.size[i]){
      std::string msg = "Matrix length mismatch during operation";
      if (error == true) {
        throw std::domain_error(msg);
        return false;
      }
    }
  }
  return true;
}
