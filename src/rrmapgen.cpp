/*
 *  rrmapgen.cpp
 *
 *
 *  Created by Andy Eschbacher on 3/01/11
 *
 */

#include "rr.h"
#include "stdmp.h"
#include "usage.h"
#include "rrmean.h"
#include "summaries.h"
#include "misc.h"
#include "stdmpout.h"

#include <fstream>
#include <cmath>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>

#include <gmpxx.h>
#include <gmp.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_fit.h>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ios;

int main(int argc, char *argv[])
{
	const char* inFilename = "input.txt";
	const char* outFilename = "rr_histogram.txt";
	double* x,* y;								// Will be dynamically declared matrices of length w
	const double TWOPI = 2.0*M_PI;
	// Parameter values
	double e = 0.05, k = 0.97163540631/TWOPI;	// Redefine k to save future computation
	unsigned long int w = 10, n = 1, diff;		// Window size and overlap
	ofstream fid;								// Output file streams

	// Program control variables
	bool silent = false, isSticky = false;
	bool publish = true, doFit = false;

	// Recurrence rate variables
	double thr = -1.0, rrLast, rrCurr, rrMean;		// Sticking Threshold, recurrence rate placeholders
	int rrcntr = 0;								// Counter for redundant rp recurrences
	mpz_t winNumMax, currWin, rrEnter, tau;		// RR time variables
	double tauDoub;	

	bool stdOut = true;
	double stdLenMin = 10.0E4;
	double stdLenMax = 10.0E9;
	double xtemp, ytemp;

	// Initialize multiprecision variables
	mpz_t l;
	mpz_init_set_ui(currWin,1);
	mpz_init(winNumMax);
	mpz_init_set_ui(l,1000);
	mpz_init_set_ui(rrEnter,1);
	mpz_init_set_ui(tau,1);

	// Initialize runtime information
	time_t t1 = time(NULL), t2, rawTime;
	struct tm * timeinfo;
	time( &rawTime );
	timeinfo = localtime( &rawTime );
	
	// Initialize histogram basics
	int nBins = 20;
	 
	/**********************************************************/
	// Read commandline input
	const char* optstring(":i:e:w:n:l:k:r:b:spq:fh?");
	if ( argc == 1 )
	{ 
		usage_rrmapgen(e,k,l,n,w,thr,nBins); 
		exit(0);
	}

    while (optind < argc)
	{
		optarg = NULL;
		int c;
		while ( (c = getopt(argc, argv, optstring)) != -1 )
		{
			switch(c)
			{
				case 'i': //initial values input
					inFilename = argv[optind-1];
					break;
				case 'e':
					e = atof(optarg);
					break;
				case 'w':
					w = atoi(optarg);
					if ( (n >= w) && (n != 1) )
					{
						cerr << "\aError: overlap is ";
						cerr << ( n > w ?"greater than":"equal to" );
						cerr << " window size." << endl;
						exit(0);
					}
					break;
				case 'n':
					n = atoi(optarg);
					if (n<1)
					{
						cerr << "\aError: overlap must be greater than or"
						"equal to 1." << endl;
						exit(0);
					} else if ( (n >= w) && (n != 10) )
					{
						cerr << "\aError: overlap is ";
						cerr << (n>w?"greater than":"equal to");
						cerr << " window size." << endl;
						exit(0);
					}
					break;
				case 'l': 
					mpz_set_str(l,optarg,10);
					break;
				case 'k':
					k = atof(optarg)/TWOPI;
					break;
				case 'r':
					thr = atof(optarg);
					break;
				case 'b':
					nBins = atoi(optarg);
					break;
				case 's':
					silent = true;
					break;
				case 'p':
					publish = false;
					break;
				case 'q':
					stdOut = true;
					stdLenMin = atof(argv[optind-1]);
					//stdLenMax = 10.0*stdLenMin;
					break;
				case 'f':
					doFit = true;
					break;
				case 'h':
				case '?':
					usage_rrmapgen(e,k,l,n,w,thr,nBins);
					exit(0);
				default:
					cerr << "Error: " << optarg << "is not a recognized "
					"option" << endl;
					exit(0);
			}
		}
    }
	
	/**********************************************************/	
	//Initialize variables based on commandline input
	x = new double[w];
	y = new double[w];
	diff = w - n;
	mpz_cdiv_q_ui(winNumMax,l,diff);
	if(thr<0.0) thr = rrmean(k,w,n,e,rrcntr);

	// Initialize histogram
	gsl_histogram * h = gsl_histogram_alloc(nBins);
	double range[nBins+1];
	logspace( range, nBins + 1, 1.999, 9.001 );
	gsl_histogram_set_ranges( h, range, nBins+1 );
	
	// open initial values file and extract initial values
	ifstream input( inFilename, ifstream::in );
	if ( !input.is_open() )
	{
		cerr << "Error: could not open initial value input file " << inFilename << "." << endl;
		cout << "Continue with sample values (x[0] = " << 0.5 << " and y[0] = " << 0.601 << ") or write a sample file to your current working directory? (1 or 2)" << endl;
		int yesno = 1;
		cin >> yesno;
		if ( yesno == 1 )
		{
			x[0] = 0.5;
			y[0] = 0.601;
		} 
		else if ( yesno == 2 )
		{
			ofstream new_input("sample_input.txt",ios::out);
			new_input << 0.5 << '\t' << 0.601 << endl;
			new_input.close();
			cerr << "A sample input file called sample_input.txt was created.\a" << endl;
		}
		else
		{
			cerr << "Error: " << yesno << " is not a valid option."
			"Exiting with no actions taken." << endl;
			exit(0);
		}
	}
	else
	{
		input >> x[0];
		input >> y[0];
		input.close();
	}

	double xInit = x[0];
	double yInit = y[0];

	/**********************************************************/
	// Print summary of input values
	if (!silent)
	{
		cout << endl;
		cout << argv[0] << " summary of inputs: \n" << endl;
		cout << "\tkicking strength (k) = \t " << k*TWOPI << endl;
		cout << "\trp threshold (e) = \t " << e << endl;
		cout << "\twindow size (w) = \t " << w << endl;
		cout << "\twindow overlap (n) = \t " << n << endl;
		cout << "\ttrajectory length (l) =  " << l << endl;
		cout << "\tnumber of rr values = \t " << winNumMax << endl;
		cout << "\tsticking threshold (r) = " << thr << endl;
		cout << "\thistogram bins (b) = \t " << nBins << endl;
		cout << "\tPrinting stdmp trajectories for: " << stdLenMin << " <= t < " << stdLenMax << endl;
		cout << "\tcalculate linear fit: \t " << (doFit?"Yes":"No") << endl;
		cout << "\n\t" << inFilename << " contained the following values: " << endl;
		cout << "\tx[0] = " << x[0] << "\ty[0] = " << y[0] << endl;
		cout << "\nDate and time rrstdmp was initialized: " << asctime(timeinfo) << endl;	
	}

	/**********************************************************/
	// Kernel of program
	
	// Initial data (without overlap)
	stdmp(k,x,y,w);
	rrCurr = rrInit(e,x,y,w,n,rrcntr);

	// Begin analysis before sticky event
	while ( (rrCurr) > thr && mpz_cmp(winNumMax, currWin) )
	{
		stdmp(k,x,y,w,n);
		rrCurr = rr(e,x,y,w,n,rrcntr);
		mpz_add_ui(currWin,currWin,1);
	}

	if (!mpz_cmp(winNumMax,currWin))
	{
		cerr << "\n\aError:\tThreshold (r) larger than largest recurrence rate value.  No data output.\n";
		cerr << "\tConsider a smaller rr threshold (r) or larger data set (l)." << endl;
		publish = false;
		exit(0);
	}

	for (currWin; mpz_cmp(winNumMax, currWin); mpz_add_ui(currWin,currWin,1))
	{
		stdmp(k,x,y,w,n);	// Generate new data
		rrLast = rrCurr;	// Store last rr value
		rrCurr = rr(e,x,y,w,n,rrcntr);	// Assign current rr value
		
		if (rrCurr >= thr && rrLast < thr)	// Entering a sticky region?
		{									// If yes,
			mpz_set(rrEnter,currWin);		// assign time when sticky region entered
			isSticky = true;				// and tell loop you are in a sticky region
			xtemp = x[0];					// Initial values where sticky region begins
			ytemp = y[0];
			
		} else if (rrCurr < thr && rrLast > thr)	// Exiting a sticky region?
		{
			mpz_sub(tau,currWin,rrEnter);						// Assign the length of sticky region to tau
			tauDoub = static_cast<double> (mpz_get_ui(tau));	
			rrMean = (rrMean + rrCurr)/tauDoub;					// Calculate mean recurrence in sticky region
			tauDoub = static_cast<double> (diff*mpz_get_ui(tau) + n);	// Convert tau to type double, 
																		// and convert to time frame of standard map
			gsl_histogram_increment(h,tauDoub); // Insert sticky region length into histogram

			if ( stdOut && (tauDoub > stdLenMin) && ( tauDoub < stdLenMax ) )	// Is this a major sticky event?
				stdmpout(diff*mpz_get_ui(tau)+n,xtemp,ytemp,k,rrMean); // if so, output stdmap trajectory

			rrMean = 0.0;
			isSticky = false;
		}
		
		if (isSticky) rrMean += rrCurr;
	}

	delete [] x;
	delete [] y;
	x = NULL;
	y = NULL;
	/**********************************************************/
	// Construct CDF

	double histSum = gsl_histogram_sum(h);
	std::vector<int> binTimes;
	
	// Find non-zero bins and store their indices in the vector binTimes[]
	for ( int i = 0; i < nBins; i++ )
		if ( gsl_histogram_get(h,i) )
			binTimes.push_back(i);
	
	// Create vectors to store the x and y histogram values that are
	//   modified and copied from the histogram information
	x = new double[binTimes.size()];
	y = new double[binTimes.size()];
	
	for ( int j = binTimes.size()-1; j >= 0; j-- )
	{
		y[j] = gsl_histogram_get(h,binTimes[j]);
		x[j] = ( range[binTimes[j]] + range[binTimes[j]+1] ) / 2.0;
	}	
	
	// Accumulate y component for CDF
	cumSumNorm(y,binTimes.size(),histSum);
	
	// Print to file
	fid.open(outFilename,ios::out);
	if (fid.is_open())
	{
		for (int i = 0; i < binTimes.size(); i++)
			fid << log10(x[i]) << '\t' << log10(y[i]) << endl;
		fid.close();
	} else
	{
		cerr << "Error: Could not open " << outFilename << endl;
		cout << "\tPrinting to screen..." << endl;
		for (int i = 0; i<binTimes.size(); i++)
			fid << log10(x[i]) << '\t' << log10(y[i]) << endl;
	}
	
	/**********************************************************/
	// Perform least square linear fit if chosen
	if (doFit)
	{
		for (int i = 0; i < binTimes.size(); i++)
		{
			x[i] = log10(x[i]);
			y[i] = log10(y[i]);
		}

		double b, m, cov00, cov01, cov11, sumsq;
		gsl_fit_linear(x,1,y,1,binTimes.size(),&b,&m, \
					   &cov00,&cov01,&cov11,&sumsq);
		cout << "Ordinary log-log fit:" << endl;
		cout << '\t' << m << "x + " << b << endl;
		cout << "MLE fit:" << endl;
		cout << '\t' << mlefit(x,binTimes.size()) << endl;
	}
	

	/**********************************************************/
	// End of program

	delete[] x;
	delete[] y;

	t2 = time(NULL) - t1;
	t1 = t2;
	if (!silent) timesummary(t1,t2);

	//Summarize results into log
	if (publish)
	{
		ofstream resultsLog("rrmapgen_results.log",ios::app);
		if (resultsLog.is_open())
		{
			resultsLog << k*TWOPI << "," << e << "," << mpz_get_d(l) << "," << w << "," << n << ",";
			resultsLog << xInit << "," << yInit << "," << t1 << "," << asctime(timeinfo);
			resultsLog.close();
			if (!silent) cout << "Input summary written to rrmapgen_results.log" << endl;
		} else
		{
			cerr << "Error: Could not open rrmapgen_results.log" << endl;
		}
	}

	gsl_histogram_free(h);
	mpz_clear(l);
	
	return 0;
}
