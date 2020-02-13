#include <iostream>

#include "Vec.h"

#ifndef MAT_H
#define MAT_H

class Mat {

  protected:

    double* values;
    int dims = 0;
    int len = 1;
    int* size;
    int* idx_increments;
    double tol = 1E-6;

  public:

    // CONSTRUCTORS

    Mat(const int& d1, const int& d2 = 0, const int& d3 = 0, const int& d4 = 0, const int& d5 = 0, const int& d6 = 0);

    Mat(const Mat& u);

    // DESTRUCTOR

    ~Mat();

    // DATA

    int length();

    int dimension();

    Vec shape();

    Mat copy();

    // INDEXING
    int locate(const int& idx1, const int& idx2 = 0, const int& idx3 = 0, const int& idx4 = 0, const int& idx5 = 0, const int& idx6 = 0);

    void set(const double& val, const int& idx1, const int& idx2 = 0, const int& idx3 = 0, const int& idx4 = 0, const int& idx5 = 0, const int& idx6 = 0);

    double get(const int& idx1, const int& idx2 = 0, const int& idx3 = 0, const int& idx4 = 0, const int& idx5 = 0, const int& idx6 = 0);

    void set_raw(const double& val, const int& idx);

    double get_raw(const int& idx);

    Mat transpose();

    Mat transpose(const int& idx1, const int& idx2 = 0, const int& idx3 = 0, const int& idx4 = 0, const int& idx5 = 0, const int& idx6 = 0);

    // TESTING

    bool check_size(const Mat& u, const bool& error = true);

    bool is_vector(const Mat& u, const bool& error = true);

    // VISUALIZATION

    std::string string();

    void print();

    Mat operator+(const Mat& u);

    Mat operator-(const Mat& u);

    Mat operator-();

    Mat operator*(const double& a);

    Mat operator/(const double& a);

    Mat prod(const Mat& u);

    // INPLACE OPERATORS

    void operator+=(const Mat& u);

    void operator-=(const Mat& u);

    void operator*=(const double& a);

    // LINEAR ALGEBRA

    Mat mat_mul(Mat u);

    double inner(Mat u);

    double norm();

    Mat unit();

    double determinant();

    Mat inverse();

    Mat adjugate();

    Mat cofactor();

    Mat outer(const Mat& u);

    Mat rref();

    // COMPARATORS;

    bool operator==(const Mat& u);

    bool operator!=(const Mat& u);

    // ASSIGNMENT OPERATORS

    void operator=(const std::string& s);

    void operator=(const Mat& u);

};

#endif
