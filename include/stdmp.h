/*
 *  stdmp.h
 *  
 *
 *  Created by Andy Eschbacher on 2010-10-12.
 *  
 *
 */
#ifndef STDMP_H
#define STDMP_H

#include <cmath>
#include <iostream>

#include <omp.h>

extern const double TWOPI;

void stdmpInit(double*, double*);

void stdmp(double *__restrict__ x, double *__restrict__ y);

void stdmpBack(double *__restrict__ x, double *__restrict__ y);

double stdmpLifted(double, double, const int);

#endif // STDMP_H
