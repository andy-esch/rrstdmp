
#ifndef RRMEAN_H
#define RRMEAN_H

#include "rr.h"
#include "stdmp.h"

#include <cmath>
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ios;
using std::setw;

extern const double TWOPI;

double f(const double, const double, const int, const int);
double secantApprox(double, double, const double, const int, \
					int, int, const double);
double rrmean(int&);
void rrtester(int, int, double, int&);

#endif // RRMEAN_H
