
#ifndef MISC_H
#define MISC_H

#include <cmath>
#include <iostream>
#include <vector>

#include <gsl/gsl_histogram.h>

using std::vector;
using std::cout;
using std::endl;

void logspace(double*, const int, const double, const double);
void linspace(double *, const int, const double, const double);
void cumSumNorm(double*, const int, const double);
double mlefit(double*, const int);
//int fibgen(int, int);
void rrcdf(gsl_histogram*, double*, double*);
void printHistogram(gsl_histogram*);
void printXY(double *, double *, int);

#endif // MISC_H
