#include <iostream>
#include <cmath>
#include "Vec.h"
#include "Mat.h"
#include "tools.h"
#include <cassert>

// CONSTRUCTOR

Mat::Mat(const int& d1, const int& d2, const int& d3, const int& d4, const int& d5, const int& d6) {
  bool prev_0 = false;
  size = new int[6] {d1, d2, d3, d4, d5, d6};

  for (int i = 0; i < 6; i++) {
    if (size[i] > tol && prev_0 == true) {
      std::string msg = "Cannot have length of zero in a lower dimension!";
      throw std::invalid_argument(msg);
    } else if (size[i] == 0) {
      prev_0 = true;
    } else {
      dims += 1;
      len *= size[i];
    }
  }

  idx_increments = new int[6] {0, 0, 0, 0, 0, 0};
  for (int i = 0; i < dims; i++) {
    idx_increments[i] = 1;
    for (int j = i+1; j < dims; j++) {
      idx_increments[i] *= size[j];
    }
  }

  values = new double[len];
  for (int i = 0; i < len; i++) {
    values[i] = 0;
  }
}

Mat::Mat(const Mat& u) {
  size = new int[6];
  for (int i = 0; i < 6; i++) {
    size[i] = u.size[i];
  }

  len = u.len;
  dims = u.dims;

  values = new double[len];

  for (int i = 0; i < len; i++) {
    values[i] = u.values[i];
  }

  idx_increments = new int[6];
  for (int i = 0; i < 6; i++) {
    idx_increments[i] = u.idx_increments[i];
  }
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

Vec Mat::shape() {
  Vec v(dims);
  for (int i = 0; i < dims; i++) {
    v.set(size[i], i);
  }
  return v;
}

Mat Mat::copy() {
  Mat v(len);
  for (int i = 0; i < len; i++) {
    v.set(get(i), i);
  }
  return v;
}

// INDEXING

int Mat::locate(const int& idx1, const int& idx2, const int& idx3, const int& idx4, const int& idx5, const int& idx6) {
  // Converts N-dimensional indices to 1-D indices
  /*
  [[1,2,3]         [[(0,0) (0,1) (0,2)]
   [4,5,6]]         [(1,0) (1,1) (1,2)]]


   [1,2,3,4,5,6]    [0,1,2,3,4,5]
  */

  return idx1*idx_increments[0] + idx2*idx_increments[1]
       + idx3*idx_increments[2] + idx4*idx_increments[3]
       + idx5*idx_increments[4] + idx6*idx_increments[5];
}

void Mat::set(const double& val, const int& idx1, const int& idx2, const int& idx3, const int& idx4, const int& idx5, const int& idx6) {
  values[locate(idx1, idx2, idx3, idx4, idx5, idx6)] = val;
}

double Mat::get(const int& idx1, const int& idx2, const int& idx3, const int& idx4, const int& idx5, const int& idx6) {
  // Returns array value at given indices
  return values[locate(idx1, idx2, idx3, idx4, idx5, idx6)];
}

void Mat::set_raw(const double& val, const int& idx) {
  // Set an element from the "values" array directly
  values[idx] = val;
}

double Mat::get_raw(const int& idx) {
  // Get an element from the "values" array directly
  return values[idx];
}

Mat Mat::transpose() {
  if (dims > 2) {
    std::string msg = "Must specify transposition axes for greater than two dimensions";
    throw std::domain_error(msg);
  } else if (dims == 1) {
    Mat v(1, size[0]);
    #pragma omp parallel for
    for (int i = 0; i < size[0]; i++) {
      v.values[i] = values[i];
    }
    return v;
  } else {
    Mat v(size[1], size[0]);
    #pragma omp parallel for
    for (int i=0; i < size[1]; i++) {
      for (int j=0; j < size[0]; j++) {
        v.set(get(j,i), i, j);
      }
    }
    return v;
  }
}

Mat Mat::transpose(const int& idx1, const int& idx2, const int& idx3, const int& idx4, const int& idx5, const int& idx6) {
  int indices[6] = {idx1, idx2, idx3, idx4, idx5, idx6};
  int new_idx[6] = {0,0,0,0,0,0};
  bool idx_check[6] = {false, false, false, false, false};
  int check = 0;
  bool slot = false;
  #pragma omp parallel for
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      if (size[i] == indices[j] && idx_check[j] == false) {
        check++;
        idx_check[j] = true;
        new_idx[i] = j;
        break;
      }
      if (slot == false) {
        new_idx[i] = i;
      }
      slot = false;
    }
  }
  #pragma omp parallel for
  for (int i = 0; i < 6; i++) {
    if (indices[i] == 0) {
      indices[i]++;
    }
  }
  int idx[6] = {0,0,0,0,0,0};
  Mat v(idx1, idx2, idx3, idx4, idx5, idx6);
  #pragma omp parallel for
  for (int i6 = 0; i6 < indices[5]; i6++) {
    idx[5] = i6;
    for (int i5 = 0; i5 < indices[4]; i5++) {
      idx[4] = i5;
      for (int i4 = 0; i4 < indices[3]; i4++) {
        idx[3] = i4;
        for (int i3 = 0; i3 < indices[2]; i3++) {
          idx[2] = i3;
          for (int i2 = 0; i2 < indices[1]; i2++) {
            idx[1] = i2;
            for (int i1 = 0; i1 < indices[0]; i1++) {
              idx[0] = i1;
              v.set(get(idx[new_idx[0]], idx[new_idx[1]], idx[new_idx[2]], idx[new_idx[3]], idx[new_idx[4]], idx[new_idx[5]]), i1,i2,i3,i4,i5);
            }
          }
        }
      }
    }
  }
  return v;
}

// TESTING

bool Mat::check_size(const Mat& u, const bool& error) {
  for (int i = 0; i < 6; i++) {
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

bool Mat::is_vector(const Mat& u, const bool& error) {
  std::string msg = "Invalid Mat object;\n";
  msg += "Vector expected";
  for (int i = 3; i < 6; i++) {
    if (u.size[i] > tol) {
      if (error == true) {
        throw std::invalid_argument(msg);
      }
      return false;
    }
  }
  if (size[0] == 1 && size[1] > 1) {
    return true;
  } else if (size[0] > 1 && (size[1] == 1 || size[1] == 0)) {
    return true;
  }
  if (error == true) {
    throw std::invalid_argument(msg);
  }
  return false;
}

// VISUALIZATION

std::string Mat::string() {
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

void Mat::print() {
  // Prints 1-D and 2-D arrays to the terminal
  // For higher dimensions, raises an error
  if (dims > 2) {
    std::string msg = "Attempting to print a higher dimensional array.";
    msg += "\nCan only print 1-D and 2-D arrays.";
    throw std::domain_error(msg);
  }
  else if (dims == 1)
  {
    std::cout << "[";
    for (int i = 0; i < len; i++) {
      std::cout << sci_not(values[i], 1, 2);
      if (i < len - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "]" << std::endl;
  }
  else if (dims == 2)
  {
    std::cout << "[\n";
    for (int i = 0; i < size[0]; i++) {
      std::cout << "[";
      for (int j = 0; j < size[1]; j++) {
        std::cout << sci_not(get(i, j), 1, 2);
        if (j < size[1]-1) {
          std::cout << ", ";
        }
      }
      std::cout << "]" << std::endl;
    }
    std::cout << "]" << std::endl;
  }
}

// ARITHMETIC

Mat Mat::operator+ (const Mat& u){
  check_size(u);
  Mat v(size[0], size[1], size[2], size[3], size[4], size[5]);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    v.values[i] = values[i] + u.values[i];
  }
  return v;
}

Mat Mat::operator- (const Mat& u){
  check_size(u);
  Mat v(size[0], size[1], size[2], size[3], size[4], size[5]);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    v.values[i] = values[i] - u.values[i];
  }
  return v;
}

Mat Mat::operator- (){
  Mat v(size[0], size[1], size[2], size[3], size[4], size[5]);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    v.values[i] = -values[i];
  }
  return v;
}

Mat Mat::operator* (const double& a) {
  Mat v(size[0], size[1], size[2], size[3], size[4], size[5]);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    v.values[i] = a*values[i];
  }
  return v;
}

Mat Mat::operator/ (const double& a) {
  Mat v(size[0], size[1], size[2], size[3], size[4], size[5]);
  double reciprocal = 1/a;
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    v.values[i] = values[i]*reciprocal;
  }
  return v;
}

Mat Mat::prod(const Mat& u) {
  check_size(u);
  Mat v(size[0], size[1], size[2], size[3], size[4], size[5]);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    v.values[i] = values[i]*u.values[i];
  }
  return v;
}

// INPLACE OPERATORS

void Mat::operator+= (const Mat& u) {
  check_size(u);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    values[i] = values[i] + u.values[i];
  }
}

void Mat::operator-= (const Mat& u) {
  check_size(u);
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    values[i] = values[i] - u.values[i];
  }
}

void Mat::operator*= (const double& a) {
  #pragma omp parallel for
  for (int i = 0; i < len; i++) {
    values[i] = a*values[i];
  }
}

// LINEAR ALGEBRA

Mat Mat::mat_mul(Mat u){
  if (dims > 2 || u.dims > 2) {
    std::string msg = "Matrix multiplication is only available for 1-D and 2-D Matrices.";
    throw std::domain_error(msg);
  }
  int d0 = size[0];
  int d1 = size[1];
  int u_d0 = u.size[0];
  int u_d1 = u.size[1];

  if (d1 == 0) {
    d1 = 1;
  }
  if (u_d1 == 0) {
    u_d1 = 1;
  }

  if (d1 != u_d0) {
    std::string msg = "Invalid matrix dimensions for matrix multiplication.";
    throw std::domain_error(msg);
  }

  Mat v(d0, u_d1);
  double total;
  #pragma omp parallel for
  for (int i = 0; i < d0; i++) {
    for (int j = 0; j < u_d1; j++) {
      total = 0;
      for (int k = 0; k < d1; k++) {
        total += get(i,k)*u.get(k,j);
      }
      v.set(total, i, j);
    }
  }
  return v;
}

double Mat::inner(Mat u) {
  check_size(u);
  double a = 0;
  #pragma omp parallel for
  for (int i = 0; i < len; i++){
    a += values[i]*u.values[i];
  }
  return a;
}

double Mat::norm() {
  is_vector(*this);
  double total = 0;
  #pragma omp parallel for
  for (int i = 0; i < len; i++) {
    total += get_raw(i)*get_raw(i);
  }
  return std::sqrt(total);
}

Mat Mat::unit() {
  is_vector(*this);
  Mat v(size[0], size[1], size[2], size[3], size[4], size[5]);
  double u_norm = norm();
  if (u_norm == 0) {
    return v;
  }
  #pragma omp parallel for
  for (int i = 0; i < len; i++) {
    v.set_raw(i, get_raw(i)/u_norm);
  }
  return v;
}

double Mat::determinant() {
  if (dims != 2) {
    std::string msg = "Can only take the determinant of a 2-D matrix";
    throw std::domain_error(msg);
  }

  if (size[0] != size[1]) {
    std::string msg = "Can only take the determinant of a square matrix";
    throw std::domain_error(msg);
  }

  Mat pm(size[0], size[1]);
  int val = 1;
  int val_tmp;
  for (int i = 0; i < size[0]; i++) {
    val_tmp = val;
    for (int j = 0; j < size[1]; j++) {
      pm.set(val_tmp, i, j);
      val_tmp *= -1;
    }
    val *= -1;
  }

  if (size[0] == 2 && size[1] == 2) {
    return get(0,0)*get(1,1) - get(1,0)*get(0,1);
  } else {
    double total = 0;
    Mat v(size[0]-1, size[1]-1);
    double p = 0;
    for (int i = 0; i < size[0]; i++) {
      p = get(0,i);
      v = v*0;
      for (int j = 1; j < size[0]; j++) {
        #pragma omp parallel for
        for (int k = 0; k < size[0]; k++) {
          if (k < i) {
            v.set(get(j,k), j-1, k);
          } else if (k > i) {
            v.set(get(j,k), j-1, k-1);
          }
        }
        total += v.determinant()*p*pm.get(0,i);
      }
    }
    return total;
  }

}

Mat Mat::inverse() {
  if (dims != 2) {
    std::string msg = "Can only take the inverse of a 2-D matrix";
    throw std::domain_error(msg);
  } else if (size[0] != size[1]) {
    std::string msg = "Can only calculate the inverse of a square matrix";
    throw std::domain_error(msg);
  } else if (size[0] == 1) {
    std::string msg = "Can only calculate inverse of a 2x2 matrix or larger";
    throw std::domain_error(msg);
  }

  double det = determinant();
  if (det < tol) {
    std::string msg = "Matrices with a determinant of zero have no inverse";
    throw std::domain_error(msg);
  }

  Mat u = adjugate()/det;
  return u;
}

Mat Mat::adjugate() {
  if (dims != 2) {
    std::string msg = "Can only take the adjoint of a 2-D matrix";
    throw std::domain_error(msg);
  } else if (size[0] != size[1]) {
    std::string msg = "Can only calculate the adjoint of a square matrix";
    throw std::domain_error(msg);
  } else if (size[0] == 1) {
    std::string msg = "Can only calculate adjoint of a 2x2 matrix or larger";
    throw std::domain_error(msg);
  }
  Mat u = cofactor();
  return u.transpose();
}

Mat Mat::cofactor() {
  if (dims != 2) {
    std::string msg = "Can only take the cofactor of a 2-D matrix";
    throw std::domain_error(msg);
  } else if (size[0] != size[1]) {
    std::string msg = "Can only calculate the cofactor of a square matrix";
    throw std::domain_error(msg);
  } else if (size[0] == 1) {
    std::string msg = "Can only calculate cofactor of a 2x2 matrix or larger";
    throw std::domain_error(msg);
  }

  if (size[0] == 2) {
    Mat v(2,2);
    v.values[0] = values[3];
    v.values[1] = -values[1];
    v.values[2] = -values[2];
    v.values[3] = values[0];
    return v;
  } else {
    Mat u(size[0], size[0]);
    Mat v(size[0]-1, size[1]-1);
    int idx1; int idx2; int i; int j; int k; int l;
    int sign = 1; int sign_row = 1;
    for (i = 0; i < size[0]; i++) {
      for (j = 0; j < size[1]; j++) {
        for (k = 0; k < size[0]; k++) {
          #pragma omp parallel for
          for (l = 0; l < size[1]; l++) {
            if(k == i || l == j) {
              continue;
            }
            if (k < i) {
              idx1 = k;
            } else if(k > i) {
              idx1 = k-1;
            }
            if (l < j) {
              idx2 = l;
            } else if(l > j) {
              idx2 = l-1;
            }
            v.set(get(k,l), idx1, idx2);
          }
        }
        u.set(sign*v.determinant(), i,j);
        sign *= -1;
      }
      sign_row *= -1;
      sign = sign_row;
    }
    return u;
  }
}

Mat Mat::outer(const Mat& u) {
  if (dims == 2 && size[0] == 1 && size[1] >= 1 && u.dims == 1) {
    Mat v(size[1], u.size[0]);
    #pragma omp parallel for
    for (int i = 0; i < size[1]; i++) {
      for (int j = 0; j < u.size[0]; j++) {
        v.set(get(0,i)*u.values[j], i, j);
      }
    }
    return v;
  } else if (dims == u.dims && dims == 1) {
    Mat v(size[0], u.size[0]);
    #pragma omp parallel for
    for (int i = 0; i < size[0]; i++) {
      for (int j = 0; j < u.size[0]; j++) {
        v.set(get(i)*u.values[j], i, j);
      }
    }
    return v;
  } else {
    std::string msg = "Can only take the outer product of row & column matrices";
    throw std::domain_error(msg);
  }
}

Mat Mat::rref() {
  if (dims != 2) {
    std::string msg = "Can only rref a 2-D matrix";
    throw std::domain_error(msg);
  } else if (size[0] < 2 || size[1] < 2) {
    std::string msg = "Can only rref a 2x2 matrix or larger";
    throw std::domain_error(msg);
  }
  int rows[size[0]];
  for (int i = 0; i < size[0]; i++) {
    rows[i] = i;
  }
  Mat u(size[0], size[1]);
  int tmp;
  for (int i = 0; i < size[0]-1; i++) {
    if (get(i,i) == 0) {
      for (int j = i+1; j < size[0]; j++) {
        if (get(j,i) != 0) {
          tmp = rows[i];
          rows[i] = rows[j];
          rows[j] = tmp;
          break;
        }
      }
    }
  }
  #pragma omp parallel for
  for (int i = 0; i < size[0]; i++) {
    for (int j = 0; j < size[1]; j++) {
      u.set(get(i,j),rows[i],j);
    }
  }

  double c;
  int i; int j; int k;
  #pragma omp parallel for
  for (i = 0; i < size[0]; i++) {
    for (j = 0; j < size[0]; j++) {
      if (i == j || u.get(i,i) == 0) {
        continue;
      } else {
        c = u.get(j,i)/u.get(i,i);
        for (k = 0; k < size[1]; k++) {
          u.set(u.get(j,k) - u.get(i,k)*c, j, k);
        }
      }
    }
  }
  #pragma omp parallel for
  for (i = 0; i < size[0]; i++) {
    for (j = 0; j < size[1]; j++) {
      if (u.get(i,i) == 0) {
        continue;
      }
      u.set(u.get(i,j)/u.get(i,i), i, j);
    }
  }
  return u;
}

// COMPARATORS

bool Mat::operator== (const Mat& u){
  if (check_size(u, false) == false) {
    return false;
  }
  for (int i = 0; i < len; i++){
    if (values[i] - u.values[i] > tol) {
      std::cout << values[i] << " " << u.values[i] << std::endl;
      std::cout << values[i] - u.values[i] << std::endl;
      return false;
    }
  }
  return true;
}

bool Mat::operator!= (const Mat& u){
  if (check_size(u, false) == false) {
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
  for (int i = 0; i < s.length(); i++) {

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
  for (int i = 0; i < s.length(); i++) {
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
