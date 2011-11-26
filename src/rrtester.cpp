/*
 *  rrtester.cpp
 *  
 *
 *  Created by Andy Eschbacher on 4/29/11.
 *
 */

#include "rrtester.h"

void rrtester(int& rrcntr)
{
	extern int globalWindow;
	extern double k;
	double u[globalWindow], v[globalWindow];
	double rrmean = 0.0, testdiff = 0.0, tempdiff = 0.0;
	const int iters = 200, num = 20;

	cout << "Testing calibration of rr()..." << endl;
	cout << setw(10) << "Expected" << setw(12) << " Actual  " << setw(10) << "Difference" << endl;
	cout << setw(10) << "--------" << setw(12) << " ------  " << setw(10) << "----------" << endl;

	for (int i = 1; i < num; i++)
	{
		u[0] = 0.0;
		v[0] = 1.0 / static_cast<double> (i);
//		stdmpInit(0.0,u,v,w);	// These need to be created from scratch, and follow the older definitions of the functions because they are calling a k = 0 option.
		rrmean = rrInit(u,v,rrcntr);
		for (int j = 0; j < iters; j++)
		{
//			stdmp(0.0,u,v,w,n);
			rrmean += rr(u,v,rrcntr);
		}
		tempdiff = fabs(1.0/static_cast<double> (i) - rrmean/(static_cast<double> (iters + 1)));
		testdiff += tempdiff;
		
		cout << std::fixed << std::setprecision(5) << setw(10) << 1.0/i << setw(10) << rrmean/(iters+1) << setw(12) << std::scientific << std::setprecision(3) << tempdiff << endl;
		rrmean = 0.0;
		rrcntr = 0.0;
	}

	testdiff /= static_cast<double> (num);
	cout << endl;
	cout << "--->  rr() is accurate to within " << std::fixed << testdiff*100 << "% on average." << endl;
}
