/*

 This program generates two vectors (x and y) that form the coordinate pair for
   the phase space of the Chirikov Standard Map.  The first function, without n,
   generates a standard map trajectory for the initial window, while the second
   function stdmp(), with n, captures the overlap values from the previous stdmp 
   call, and calculates from n to w.  The final function, stdmpLifted(), 
   calculates the map lifted off the modulus.  It is used in the rrmean.cpp file 
   for use in finding the golden KAM curve.

  N.B.: All maps assume k2pi is k/(2*pi)

 */

#include "stdmp.h"
#include <cmath>
#include <iostream>

const double TWOPI = 2.0 * M_PI;

void stdmpInit(double *__restrict__ x, double *__restrict__ y)
{
	extern double k;
	extern int globalWindow;
	static double modCorr = 100.0;	// Offset for values that stray below zero
									// ceil(1.75 + k2pi) is sufficient
	for ( int i = 1; i < globalWindow; i++)
	{
		y[i] = fmod(y[i-1] + k*sin( TWOPI * x[i-1] ) + modCorr, 1.0);
		x[i] = fmod(y[i] + x[i-1] + modCorr, 1.0);
 	}
}

void stdmp(double *__restrict__ x, double *__restrict__ y)
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

// Calculates the lifted standard map trajectory, returns the (w-1)st value of x
// Used with secantApprox(...)
//double stdmpLifted(const double k2pi, double x1, double y1, const int w)
double stdmpLifted(double x1, double y1, const int iters)
{
	extern double k;

	for ( int i = 1; i < iters; i++ )
	{
		y1 = y1 + k * sin( TWOPI * x1 );
		x1 = x1 + y1;
	}

	return x1;
}

/*
void noise(double* x, double* y, const int w)
 {
	srand( time(NULL) );
	
	for (int i = 1; i < w; i++)
	{
		x[i] = static_cast<double> (rand())/static_cast<double> (RAND_MAX);
		y[i] = static_cast<double> (rand())/static_cast<double> (RAND_MAX);
	}
}

void noise(double* x, double* y, const int w, const int n)
{
	srand( time(NULL) );
	int i, diff = w-n;
	
	for (i = diff; i < w; i++) // copies n overlap values
	{
		x[i - diff] = x[i];
		y[i - diff] = y[i];
	}
	
	for (int i = n; i < w; i++)
	{
		x[i] = static_cast<double> (rand())/static_cast<double> (RAND_MAX);
		y[i] = static_cast<double> (rand())/static_cast<double> (RAND_MAX);
	}	
	
}


To be logistic map generators, where the x and y components follow
     seperate maps

void logistic(double r, double* x, double* y, const int w)
{

}

void logistic(double r, double* x, double* y, const int w, const int n)
 {

} 
 
 */