/*
 *  rrstdmp.cpp
 *  
 *
 *  Created by Andy Eschbacher on 3/01/11.
 *
 */

#include "misc.h"
#include "rr.h"
#include "rrmean.h"
#include "rrtester.h"
#include "stdmp.h"
#include "summaries.h"
#include "usage.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#include <gmpxx.h>
#include <gmp.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_fit.h>

using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::vector;

double k = 0.97163540631/(2.0 * M_PI), ge = 0.05;	// redefined k, rr threshold; global/extern
int globalWindow = 100, globalOverlap = 50;		// window, overlap; global/extern

int main(int argc, char **argv)
{
	char inFilename[] = "input1";               // This line and the one below should probably be char var[30] or something like that
	char outFilename[] = "rr_histogram.txt";    //  to avoid memory problems if the file name is longer than 'input1', etc.
	double* x,* y;								// Will be dynamically declared matrices of length w
	unsigned long int diff;						// Window size and overlap
	ofstream fid;								// Output file streams
	
	// Program control variables
	bool silent = false, isSticky = false;
	bool doFit = true, publish = true;
	bool rrtest = false;

	// Recurrence rate variables
	double thr = 0.0, rrLast, rrCurr;		// Sticking Threshold, recurrence rate placeholders
	int rrcntr = 0;							// Counter for redundant rp recurrences
	mpz_t currWin, rrEnter, tau;			// RR time variables
	mpz_t winNumMax, l;						// Trajectory length
	double tauDoub;							// Length of individual sticy event as a double

	// Initialize multiprecision variables
	mpz_init(winNumMax);
	mpz_init_set_ui(currWin,1);
	mpz_init_set_ui(rrEnter,1);
	mpz_init_set_ui(tau,1);
	mpz_init_set_ui(l,25000000);

	// Initialize runtime information
	time_t t1 = time(NULL), t2, rawTime;
	struct tm * timeinfo;
	time( &rawTime );
	timeinfo = localtime( &rawTime );

	// Initialize histogram basics
	int nBins = 50;

	/**********************************************************/
	// Read commandline input
	const char* optstring(":i:o:e:w:n:l:k:r:b:aspfth?");
	if ( argc==1 )
	{
		cout << endl;
		cout << "Error: no options passed." << endl;
		usage_rrstdmp(ge,k,l,globalOverlap,globalWindow,thr,nBins); 
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
				case 'i': // initial values input
					strcpy(inFilename,argv[optind-1]);
					break;
				case 'o': // output file name
					strcpy(outFilename,argv[optind-1]);
					break;
				case 'e':
					ge = atof(optarg);
					break;
				case 'w':
					globalWindow = atoi(optarg);
					if ( (globalOverlap >= globalWindow) && (globalOverlap != 1) )
					{
						cerr << "\aError: window is ";
						cerr << ( globalOverlap > globalWindow ?"smaller than":"equal to" );
						cerr << " overlap size." << endl;
						exit(0);
					}
					break;
				case 'n':
					globalOverlap = atoi(optarg);
					if (globalOverlap<1)
					{
						cerr << "\aError: overlap must be greater than or \
						equal to 1." << endl;
						exit(0);
					}
					else if ( (globalOverlap >= globalWindow) && (globalWindow != 10) )
					{
						cerr << "\aError: overlap is ";
						cerr << (globalOverlap>globalWindow?"greater than":"equal to");
						cerr << " window size." << endl;
						exit(0);
					} 
					break;
				case 'l': 
					mpz_set_str(l,optarg,10);
					break;
				case 'k':
					k = atof(optarg)/(2.0 * M_PI);
					break;
				case 'r':
					thr = atof(optarg);
					if ( thr > 1.0 || thr < 0.0 )
					{
						cerr << "\aError: threshold (r) must be between 0.0 \
						and 1.0" << endl;
						exit(0);
					}
					break;
				case 'b':
					nBins = atoi(optarg);
					if (nBins < 1)
					{
						cerr << "\aError: number of bins must be larger \
						than 0." << endl;
						exit(0);
					}
					break;
				case 'a':
					/* This will take in a file that has the previous run's 
					 * details, including the last data point, so that the
					 * run being called will be a continuation of that. This
					 * function needs to bypass this inputting scheme, though
					 * so hmm..
					 */
					break;
				case 's':
					silent = true;
					break;
				case 'p':
					publish = false;
					break;
				case 'f':
					doFit = false;
					break;
				case 't': // Provisionary feature: to test the output of rr() to ensure it is correct
					rrtest = true;
					break;
				case '?':
				case 'h':
					usage_rrstdmp(ge,k,l,globalOverlap,globalWindow,thr,nBins);
					exit(0);
				default:
					cout << "Invalid entry: '" << static_cast<char> (c) << "'" << endl;
					usage_rrstdmp(ge,k,l,globalOverlap,globalWindow,thr,nBins);
					exit(0);
			}
		}
    }

	/**********************************************************/	
	//Initialize variables based on commandline

	if (rrtest)
	{	// if rrtest = true, run rr() diagnostic; exit program
		cout << "rrtester is currently broken -- exiting program. " << endl;
		//rrtester(w,n,e,rrcntr);
		exit(0);
	}

	if ( (diff = globalWindow - globalOverlap)%2 != 0 ) // Assign diff, check for constraint
	{
		cout << "\aError: w - n must be evenly divisible by 2." << endl;
		exit(0);
	}
	x = new double[globalWindow];
	y = new double[globalWindow];
	mpz_fdiv_q_ui(winNumMax,l,diff);			// Calculate number of windows
	if(thr==0.0) thr = rrmean(rrcntr);			// Calculate RR threshold

	// Initialize histogram
	gsl_histogram * h = gsl_histogram_alloc(nBins);
	double range[nBins+1];
	logspace(range,nBins+1,1.999,9.001);
	gsl_histogram_set_ranges(h,range,nBins+1);

	// open initial values file and extract initial values
	ifstream input(inFilename,ifstream::in);
	if ( !input.is_open() )
	{
		cerr << endl;
		cerr << "ERROR: Could not open initial value input file " \
		<< inFilename << "." << endl;
		cerr << "       A sample input file called 'sample_input' was "
		"created.\a" << endl;
		cerr << endl;
		ofstream new_input("sample_input",ios::out);
		new_input << 0.5 << '\t' << 0.601 << endl;
		new_input.close();
		exit(0);
	}
	input >> x[0];
	input >> y[0];
	input.close();

	/**********************************************************/
	// Print summary of input values
	if (!silent)
	{
		cout << endl;
		cout << argv[0] << " summary of inputs: \n" << endl;
		cout << "\tkicking strength (k) = \t " << k*2.0*M_PI << endl;
		cout << "\trp threshold (e) = \t " << ge << endl;
		cout << "\twindow size (w) = \t " << globalWindow << endl;
		cout << "\twindow overlap (n) = \t " << globalOverlap << endl;
		cout << "\ttrajectory length (l) =  " << l << endl;
		cout << "\tnumber of rr values = \t " << winNumMax << endl;
		cout << "\tsticking threshold (r) = " << thr << endl;
		cout << "\thistogram bins (b) = \t " << nBins << endl;
		cout << "\tcalculate linear fit: \t " << (doFit?"Yes":"No") << endl;
		cout << "\n\t'" << outFilename << "' will contain the final data." << endl;
		cout << "\n\t'" << inFilename << "' contained the following values: " << endl;
		cout << "\tx[0] = " << x[0] << "\ty[0] = " << y[0] << endl;
		cout << "\nDate & time rrstdmp was initialized: " << asctime(timeinfo) << endl;
	}

	/**********************************************************/
	// Kernel of program
	// Initial data (without overlap)
	stdmpInit(x,y);
	rrCurr = rrInit(x,y,rrcntr);

	// If begin in sticky event, exit it
	while ( rrCurr > thr && mpz_cmp(winNumMax, currWin) )
	{
		stdmp(x,y);
		rrCurr = rr(x,y,rrcntr);
		mpz_add_ui(currWin,currWin,1);
	}

	cout << "Exited existing sticky event of length " << currWin << "." << endl;

	if ( !mpz_cmp(winNumMax,currWin) )
	{
		cerr << "\n\aError:\tThreshold (r) larger than largest recurrence" 
		"rate value.  No data output.\n";
		cerr << "\tConsider a smaller rr threshold (r) or larger data set (l)." \
		<< endl;
		publish = false;
		exit(0);
	}

	while ( mpz_cmp(winNumMax, currWin) )	// while currWin < winNumMax
	{
		stdmp(x,y);
		rrLast = rrCurr;
		rrCurr = rr(x,y,rrcntr);

		if (rrCurr > thr && rrLast < thr)
		{
			mpz_set(rrEnter,currWin);
			isSticky = true;
		}
		else if (rrCurr < thr && rrLast > thr)
		{
			mpz_sub(tau,currWin,rrEnter);
			tauDoub = static_cast<double> (diff*mpz_get_ui(tau) + globalOverlap);
			gsl_histogram_increment(h,tauDoub);
			isSticky = false;
		}
		mpz_add_ui(currWin,currWin,1);
	}

	mpz_sub(tau,winNumMax,rrEnter);

	cout << "Uncompleted event of length " << tau << "." << endl;
	cout.precision(15);
	cout.setf(ios::floatfield);
	cout << "Final (x,y) = (" << x[globalWindow-1] << ", " << y[globalWindow-1] << ")" << endl;

	if (1)
	{
		char buff[50];
		// get the thread number from inFilename, increment it by one
		// then put it in %d below
		int globalOverlap = sprintf(buff,"%s_%d",inFilename,1);
		ofstream strand(buff,ios::out);
		strand << x[globalWindow-1] << '\t' << y[globalWindow-1] << endl;
		strand.close();
	}

	/* Output the preceding values into a file called inputX_n.txt, where
	 the X refers to the strain number, and the n refers to the segment of the
	 strain that the initial conditions represent the starting conditions of, 
	 therefore, subsequent calls can load the histogram and keep adding values 
	 to it and your small nightly runs can start to accumulate into hugely big
	 data series with much more reliable statistics.  Appeal to HDF5 files to
	 add this support.  When I do this, though, I need to get rid of the code
	 below that removes the bins with no counts.
	 */

	delete [] x;
	delete [] y;
	x = NULL;
	y = NULL;
	/**********************************************************/
	// Construct CDF
	cout << "Constructing CDF." << endl;
	double histSum = gsl_histogram_sum(h);
	vector<int> binTimes;

	// Find non-zero bins and store their indices in the vector binTimes[]
	for ( int i = 0; i < nBins; i++ )
		if ( gsl_histogram_get(h,i) )	// If non-zero entry, put in binTimes
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

	// Accumulate y component for CDF	--- Look into the accumulate function that's in either numerics or algorithm...
	cumSumNorm(y,binTimes.size(),histSum);

	// Print to file
	fid.open(outFilename,ios::out);
	if (fid.is_open())
	{
		for (int i = 0; i<binTimes.size(); i++)
			fid << log10(x[i]) << '\t' << log10(y[i]) << endl;
		fid.close();
	} else
	{
		cerr << "Error: Could not open " << outFilename << endl;
		cout << "\tPrinting to screen..." << endl;
		for (int i = 0; i<binTimes.size(); i++)
			cout << log10(x[i]) << '\t' << log10(y[i]) << endl;
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

		double b, m, varb, covbm, varm, sumsq;
		gsl_fit_linear(x,1,y,1, binTimes.size(), &b, &m, \
					   &varb, &covbm, &varm, &sumsq);
		double m2 = -1.0 * mlefit(x,binTimes.size());
		cout.precision(5);
		cout.setf(ios::fixed,ios::floatfield);
		cout << "+----------------------------+" << endl;
		cout << "| Ordinary log-log fit:" << endl;
		cout << "| \t(" << m << " +/- " << sqrt(varm) << ")x + (" << b << " +/- " << sqrt(varb) << " ) |" << endl;
		cout << "| MLE fit:" << endl;
		cout << "| \t(" << m2 << " +/- " << (m2-1.0)/sqrt(binTimes.size()) << ")x + " << 0.0 << "  |" << endl; 
		cout << "+----------------------------+" << endl;
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
		double xInit;
		double yInit;		
		input.open(inFilename,ios::in);
		input >> xInit;
		input >> yInit;
		input.close();
		
		ofstream resultsLog("rrstdmp_results.log",ios::app);
		if (resultsLog.is_open())
		{
			resultsLog << k*2.0*M_PI << "," << ge << "," << mpz_get_d(l) << "," \
			<< globalWindow << "," << globalOverlap << "," << xInit << "," << yInit << "," << t1 \
			<< "," << asctime(timeinfo);
/*			if (doFit)
				resultsLog << m << "," << m2 << ",";
 */
			resultsLog.close();
			if (!silent) cout << "Input summary written to rrstdmp_results.log" << endl;
		} else
		{
			cerr << "Error: Could not open rrstdmp_results.log" << endl;
			cerr << k*2.0*M_PI << "," << ge << "," << mpz_get_d(l) << "," \
			<< globalWindow << "," << globalOverlap << "," << xInit << "," << yInit << "," << t1 \
			<< "," << asctime(timeinfo);
		}
	}

	gsl_histogram_free(h);
	mpz_clear(l);
	
	return 0;
}
