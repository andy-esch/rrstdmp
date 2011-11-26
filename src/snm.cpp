
/*
 This program generates two vectors (x and y) that form the coordinate pair for
   the phase space of the Standard Non-Twist Map.  The first function, without n,
   generates a map orbit for the initial window, while the second
   function, with n, captures the overlap values from the previous snm call, 
   and calculates from n to w.
 */

#include "snm.h"

void snm(double a, double b, double* x, double* y, const int w)
{
	for( int i = 1; i < w; i++)
	{
		y[i] = y[i-1] - b*sin(TWOPI*x[i-1]);
		x[i] = fmod(x[i-1] + a*(1 - y[i]*y[i])+100.0,1.0);
	}
}

void snm(double a, double b, double* x, double* y, const int w, const int n)
{
	int i;
	int diff = w - n;
	
	for (i = diff; i < w; i++)
	{ // copies n overlap values
		x[i - diff] = x[i];
		y[i - diff] = y[i];
	} 

	for (i = n; i < w; i++)
	{
		y[i] = y[i-1] - b*sin(TWOPI*x[i-1]);
		x[i] = fmod(x[i-1] + a*(1 - y[i]*y[i])+100.0,1.0);
	}
}
