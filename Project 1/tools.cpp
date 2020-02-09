#include <iostream>
#include <cmath>
#include "Vec.h"
#include "Mat.h"

std::string sci_not(double x, int N=1, int n=2) {
  // Creates a string in scientific notation in the Python style from a double "x"

  // Storing the sign of the number
  std::string pre_sign = "";
  if (x < 0) {
    pre_sign = "-";
  }

  // Removing the sign from the number itself

  x = std::abs(x);

  // Finding the exponent for scientific notation
  int l;
  if (x == 0) {
    l = 0;
    pre_sign = "";
  } else {
    l = static_cast<int>(std::log10(x));
  }

  // Finding the sign of the exponent
  std::string post_sign = "+";
  if (0 <= x && x < 1) {
    // l -= 1;  // Why is this here?  Dafuq?
    post_sign = "-";
  }

  // Preparing the double for rounding
  int d = N+n-1;
  x *= std::pow(10, d-l);

  // Determining the modified length and rounding
  int raw_len = 0;
  if (x != 0) {
    raw_len = static_cast<int>(std::log10(x));
  }
  x = round(x);

  // Resets the rounded double to its original power
  int round_len = 0;
  if (x != 0) {
    round_len = static_cast<int>(std::log10(x));
  }

  int diff = round_len - raw_len;
  x /= std::pow(10, d+diff);
  l += diff;
  std::string pre_str = std::to_string(x).substr(0,d+2);

  int padding = N+n+1 - pre_str.length();
  if (padding > 0) {
    pre_str.append(std::string(padding, '0'));
  }

  std::string post_str = "";
  if (l < 10) {
    post_str.append("0");
  }
  post_str.append(std::to_string(abs(l)));

  if (l == 0) {
    post_sign = "+";
  }

  return pre_sign + pre_str + "E" + post_sign + post_str;

}

void print(std::string s) {
  std::cout << s << std::endl;
}

void print(Vec v) {
  v.print();
}

void print(Mat v) {
  v.print();
}

void print(int s) {
  std::cout << std::to_string(s) << std::endl;
}

void print(double s) {
  std::cout << std::to_string(s) << std::endl;
}

std::string bool_to_string(bool b) {
  if (b == true) {
    return "true";
  } else if (b == false) {
    return "false";
  } else {
    std::string msg = "Invalid argument in function \"std::string BoolToString(bool b)\"";
    throw std::invalid_argument(msg);
  }
}
