/*
 logspace creates an array of n elments that are equally spaced on a 
   log_10 scale
 linspace is the value in the exponent of logsum (without ln(10))
 cumSumNorm calculates the cumulative sum from the end of the array backwards
   while simultaneously norming each element along the way
 rrcdf creates the cumulative distribution function of the entered data
 */

#include "misc.h"

void logspace(double *__restrict__ value, const int n, const double min, \
			  const double max)
{
	double delta = (max>min)?(max - min):(min - max);
	double denom = static_cast<double> (n-1);

	for (int i = 0; i < n; i++)
		value[i] = pow(10,min + static_cast<double>(i)*delta/denom);
}

void linspace(double *__restrict__ value, const int n, const double min, \
			  const double max)
{
	double delta = (max>min)?(max-min):(min-max);
	double denom = static_cast<double> (n - 1);

	for (int i = 0; i < n; i++)
		value[i] = min + static_cast<double>(i)*delta/denom;
}

void cumSumNorm(double *__restrict__ xx, const int n, const double norm )
{
	xx[n-1] = xx[n-1]/norm;

	for ( int i = n - 2; i >= 0; i-- )
		xx[i] = xx[i]/norm + xx[i+1];
}

//double mlefit(double *__restrict__ x, int len)
//{
//	double gamm = 0.0;
//	int imin = len/8;	// Do the (p1(x) - s1(x)) analysis to find xmin
//	double xmin = x[imin];
//
//	for (int i = imin; i < len; i++)
//		gamm += x[i];
//
//	gamm = 1.0 + static_cast<double> (len - imin)/(gamm - static_cast<double> (len - imin) * xmin);
//
//	return gamm;
//}

double mlefit(double *__restrict__ x, const int len)
{
	int imin = len/6;
	double gam = 0.0, n = static_cast<double> (len - imin);

	for (int i = imin; i < len; i++)
		gam += x[i];

	gam = gam / n - log10( pow(10.0,x[imin]) - 0.5 );

	return 1.0 + 1.0/gam;
}

// Takes a histogram and returns two arrays (x,y) via pointers
void rrcdf(gsl_histogram* h, double * x, double * y)
{
    double * upper, * lower;
    int nBins = gsl_histogram_bins(h);
    double histSum = gsl_histogram_sum(h);
    vector<int> binTimes;
    
    cout << "Constructing CDF." << endl;
    
    for (int ii = 0; ii < nBins; ii++)
        if ( gsl_histogram_get(h,ii) )
            binTimes.push_back(ii);
    
    x = new double[binTimes.size()];
    y = new double[binTimes.size()];
    
    for ( int jj = binTimes.size() - 1; jj >= 0; jj-- )
    {
        y[jj] = gsl_histogram_get(h,binTimes[jj]);
        gsl_histogram_get_range(h,binTimes.size(),upper,lower);
        x[jj] = ( *upper + *lower ) / 2.0;
    }
}

void printHistogram(gsl_histogram* h)
{
    gsl_histogram_fprintf(stdout, h, "%g", "%g");
}

void printXY(double * x, double * y, int n)
{
    cout << endl;
    cout << "---------------------------" << endl;
    cout << "x\ty" << endl;
    for (int ii = 0; ii < n; ii++)
    {
        cout << x[ii] << '\t' << y[ii] << endl;
    }
    cout << "---------------------------" << endl;
    cout << endl;
}