#include "Vec.h"
#include "Mat.h"

#ifndef TOOLS_H
#define TOOLS_H

std::string sci_not(double x, int N, int n);

void print(std::string s);

void print(Vec v);

void print(Mat v);

void print(int s);

void print(double s);

std::string BoolToString(bool b);

#endif
