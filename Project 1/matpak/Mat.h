#include <iostream>

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

    Mat(const int& d1, const int& d2 = 0);

    Mat(const Mat& u);

    // DESTRUCTOR

    ~Mat();

    // DATA

    int length();

    int dimension();

    int shape0();

    int shape1();

    // INDEXING
    int locate(const int& idx1, const int& idx2);

    void set(const double& val, const int& idx1, const int& idx2);

    double get(const int& idx1, const int& idx2);

    void set_raw(const double& val, const int& idx);

    double get_raw(const int& idx);

    // ASSIGNMENT OPERATORS

    void operator=(const std::string& s);

    void operator=(const Mat& u);

    // TESTING

    bool check_size(const Mat& u, const bool& error = true);

};

#endif
