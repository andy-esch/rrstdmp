
/*
 rrmean:
 This program returns the mean recurrence rate of the golden KAM curve in the
   standard map dependent upon the window size w and threshold value e.  To set 
   a custom rrmean, set v[0] to a specific value, such as 0.5 for a two-period 
   orbit.  The result is weighted by the value multiplying the return quotient.
 */

#include "rrmean.h"

double rrmean(int &__restrict__ rrcntr)
{
	extern double k;
	extern int globalWindow;
	double u[globalWindow], v[globalWindow];
	double rrmean = 0.0;
	const int iters = 20;

	u[0] = 0.0;
	v[0] = secantApprox(0.25,0.75,1E-11,100,233,89,k);

	stdmpInit(u,v);
	rrmean = rrInit(u,v,rrcntr);

	for (int i = 1; i < iters; i++)
	{
		stdmp(u,v);
		rrmean += rr(u,v,rrcntr);
	}

	rrcntr = 0;
	cout.setf(ios::floatfield, ios::fixed); // To return main() to expected output behavior

	return ( 0.98*rrmean/static_cast<double> (iters) );
}

double f(const double k, const double yInit, const int n, const int m)
{
	return ( stdmpLifted( 0.0, yInit, n+1 ) - static_cast<double> (m) );
}

double secantApprox(double xi, double xi_1, const double e, const int iters, \
					int n, int m, const double k)
{
	double d = 0;
	double xlast= xi;
	int pass = 1;
	double fCurr, fLast;
	for (int j = 1; j <= iters; j++)
	{
		fCurr = f(k,xi,n,m);
		fLast = f(k,xi_1,n,m);

		for (int i = 1; i <= iters; i++)
		{ 
			d = ( xi - xi_1 ) / ( fCurr - fLast ) * fCurr;
			if ( fabs(d) < e )
				break;

			xi_1 = xi;
			xi = xi - d;
			fLast = fCurr;
			fCurr = f(k,xi,n,m);
		}

		if ( fabs(xi-xlast) < e )
		{
			return fmod(xi+100.0,1.0);
		}
		else
		{
			xlast = xi;
			m = n - m; // 1/gam^2 math
			n = n + m;
//			n = n + m; // 1/gam math
//			m = n - m;
		}
	}
	cout << "Convergence not found!" << endl;

    return xi;
}

//void rrtester(int w, int n, double e, int& rrcntr)
//{
//	double u[w], v[w];
//	double rrmean = 0.0, testdiff = 0.0, tempdiff = 0.0;
//	const int iters = 200, num = 20;
//
//	cout << "Testing calibration of rr()..." << endl;
//	cout << setw(10) << "Expected" << setw(12) << " Actual  " << setw(10) << "Difference" << endl;
//	cout << setw(10) << "--------" << setw(12) << " ------  " << setw(10) << "----------" << endl;
//
//	for (int i = 1; i < num; i++)
//	{
//		u[0] = 0.0;
//		v[0] = 1.0 / static_cast<double> (i);
//		stdmpInit(0.0,u,v,w);
//		rrmean = rrInit(e,u,v,w,n,rrcntr);
//		for (int j = 0; j < iters; j++)
//		{
//			stdmp(0.0,u,v,w,n);
//			rrmean += rr(e,u,v,w,n,rrcntr);
//		}
//		tempdiff = fabs(1.0/static_cast<double> (i) - rrmean/(static_cast<double> (iters + 1)));
//		testdiff += tempdiff;
//		
//		cout << std::fixed << std::setprecision(5) << setw(10) << 1.0/i << setw(10) << rrmean/(iters+1) << setw(12) << std::scientific << std::setprecision(3) << tempdiff << endl;
//		rrmean = 0.0;
//		rrcntr = 0.0;
//	}
//
//	testdiff /= static_cast<double> (num);
//	cout << endl;
//	cout << "--->  rr() is accurate to within " << std::fixed << testdiff*100 << "% on average." << endl;
//}

