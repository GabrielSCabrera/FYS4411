#ifndef VEC_H
#define VEC_H

class Vec {
  protected:
    double* values;
    int len;
    double tol = 1E-6;
  public:

    // CONSTRUCTORS

    Vec(const int& a);

    Vec(const Vec& v);

    // DESTRUCTOR

    ~Vec();

    // DATA

    int length() const;

    Vec copy();

    // INDEXING

    void set(const double& val, const int& idx);

    double get(const int& idx);

    // TESTING

    bool check_dims(const Vec& u, const bool& error = true);

    // VISUALIZATION

    std::string string();

    void print();

    Vec operator+(const Vec& u);

    Vec operator-(const Vec& u);

    Vec operator-();

    Vec operator*(const double& a);

    Vec operator/(const double& a);

    // INPLACE OPERATORS

    void operator+=(const Vec& u);

    void operator-=(const Vec& u);

    void operator*=(const double& a);

    // LINEAR ALGEBRA

    double inner(const Vec& u);

    Vec cross(const Vec& u);

    double norm();

    Vec unit();

    // COMPARATORS;

    bool operator==(const Vec& u);

    bool operator!=(const Vec& u);

    // ASSIGNMENT OPERATORS

    void operator=(const Vec& u);
};

#endif
