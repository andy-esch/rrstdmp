
#ifndef _RRMEAN_H_
#define _RRMEAN_H_

#include "rr.h"
#include "stdmp.h"

double f(const double, const double, const int, const int);
double secantApprox(double, double, const double, const int, \
					int, int, const double);
double rrmean(int&);
void rrtester(int, int, double, int&);

#endif // _RRMEAN_H_