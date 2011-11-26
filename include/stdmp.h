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

extern const double TWOPI;

void stdmpInit(double*, double*);

inline void stdmp(double *__restrict__ x, double *__restrict__ y)
{
 	extern int globalWindow, globalOverlap;
	extern double k;
	int i;
	static double modCorr = 100.0; // See note in stdmpInit above
	static int diff = globalWindow - globalOverlap;
	
	for ( i = diff; i < globalWindow; i+=2 )
	{ // copies n overlap values -- assumes (w-n) % 2 = 0
		x[i - diff] = x[i];
		y[i - diff] = y[i];
		x[i - diff + 1] = x[i + 1];
		y[i - diff + 1] = y[i + 1];
	}
    
	for ( i = globalOverlap; i < globalWindow; i+=2 )
	{	// calculates new non-overlapping values
		y[i] = fmod(y[i-1] + k*sin( TWOPI * x[i-1] ) + modCorr, 1.0);
		x[i] = fmod(y[i] + x[i-1] + modCorr, 1.0);
		y[i+1] = fmod(y[i] + k*sin( TWOPI * x[i] ) + modCorr, 1.0);
		x[i+1] = fmod(y[i+1] + x[i] + modCorr, 1.0);
	}   
}

double stdmpLifted(double, double, const int);

#endif // _STDMP_H
