
#include <fstream>
#include <iostream>
#include <cstdio>
#include "stdmp.h"

using std::ios;
using std::endl;
using std::ofstream;

void stdmpout(const int len, const double xInit, const double yInit, \
			  const double k, const double rrmn)
{
	static int fileNum = 1;
	const char* fileName = "stdout";
	char buffer[40];
	const int chunk = 10000;			// Declare manageable array length
	double xOut[chunk], yOut[chunk];	
	int maxIter = (len/chunk)?((len/chunk>10)?10:len/chunk):0;

	xOut[0] = xInit;
	yOut[0] = yInit;

	int n = sprintf(buffer,"%s%03d-len%d-rrmean%0.3f.txt", \
					fileName, fileNum, len, rrmn);
	ofstream outfile(buffer,ios::out);

	outfile << xOut[0] << '\t' << yOut[0] << endl;

	if ( maxIter > 0 )	// Print either smaller of trajectory
	{					// length 10e5 or (len/chunk)e5 to file
		for ( int j = 0; j < maxIter; j++)	
		{
			stdmp(k,xOut,yOut,chunk);

			for (int i = 1; i < chunk; i++)
				outfile << xOut[i] << '\t' << yOut[i] << endl;

			xOut[0] = xOut[chunk - 1];
			yOut[0] = yOut[chunk - 1];
		}
	}
	else	// If len is smaller than chunk, print 
	{		//   whole trajectory to file
		stdmp(k,xOut,yOut,len);
		for (int i = 0; i < len; i++)
			outfile << xOut[i] << '\t' << yOut[i] << endl;		
	}

	outfile.close();

	fileNum++;
}