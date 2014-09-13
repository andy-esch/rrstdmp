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
#include "stdmp.h"
#include "summaries.h"
#include "usage.h"
#include "cmdLineInput.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_fit.h>
//#include <omp.h>

using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::vector;

typedef unsigned short usInt;


//** Move these down into main() function
double k = 0.97163540631/(2.0 * M_PI), ge = 0.05;	    // redefined k, rr threshold; global/extern //** declare these as constants
const int globalWindow = 100, globalOverlap = 50;       // window, overlap; global/extern //** Declare these as constants
const double TWOPI = 2.0 * M_PI;

int main(int argc, char **argv)
{
	char inFilename[60] = "inputs/input1";
	char outFilename[60] = "rr_histogram.txt";
    ofstream fid;                                       // Output file streams //** Does this need a default initialization?  NULL?

    // Data and variables for run
	double x[globalWindow], y[globalWindow];            //** What about creating pointers to constant locations const double *ptr (or double const* ptr?)
    double * ptrx = x, * ptry = y;
	unsigned int diff = globalWindow - globalOverlap;   // Window size and overlap

	// Program control variables
	bool silent = false, verbose = false;
	bool doFit = true, publish = true;

	// Recurrence rate variables
	double thr = 0.0, rrLast = 0.0, rrCurr = 0.0;		// Sticking Threshold, recurrence rate placeholders
	int rrcntr = 0;                                     // Counter for redundant rp recurrences  //** Declare as pointer?  __restrict__?? //** Change this to usInt
    long int currWin = 1, rrEnter = 1, tau = 1;
    long int l = 25000000, winNumMax = l / diff;
    double wsquared = static_cast<double>(globalWindow * globalWindow);

    double dx = 0.0, dy = 0.0, dmax = 0.0, minusge = 1.0 - ge;
    int RR = 0;

    int ii = 0, jj = 0;
    const double modCorr = 64.0;
    if (verbose)
        cout << "Initialized recurrence rate variables" << endl;

	// Initialize runtime information
	time_t t1 = time(NULL), t2, rawTime;        //** Initialize t2 and rawTime to something?
	struct tm * timeinfo;
	time( &rawTime );
	timeinfo = localtime( &rawTime );
    if (verbose)
        cout << "Initialized runtime variables" << endl;

	// Initialize histogram basics
	int nBins = 50;                             //** Perhaps make this something like static_cast<int>(log10(static_cast<double>(l))) * 8?

	/**********************************************************/
	// Read commandline input
    cmdLineInput(argc, argv, inFilename, outFilename, l, thr, nBins, silent, \
                 publish, doFit);
    if (verbose)
        cout << "Input commandline input" << endl;

	/**********************************************************/	
	//Initialize variables based on commandline
    winNumMax = l / diff;

    //** Change this structure to a try/catch?
	if ( (diff = globalWindow - globalOverlap) % 2 != 0 ) // Assign diff, check for constraint
	{
		cout << "\aError: w - n must be evenly divisible by 2." << endl;
		exit(0);
	}
    else if ( globalOverlap > (globalWindow / 2) )
    {
        cerr << "\aError: n must be equal to or smaller than w / 2" << endl;
        exit(0);
    }

	if(fabs(thr - 0.0) < 1.0e-6) thr = rrmean(rrcntr);			// Calculate RR threshold
    if (verbose)
        cout << "Calculated threshold and number of windows" << endl;

	// Initialize histogram
	gsl_histogram * h = gsl_histogram_alloc(nBins);
	double range[nBins+1];
	logspace(range, nBins+1, 1.999, log10(static_cast<double>(l)) );
	gsl_histogram_set_ranges(h,range,nBins+1);
    if (verbose)
        cout << "Formed histogram struct" << endl;

	// open initial values file and extract initial values
	ifstream input(inFilename,ifstream::in);
	if ( !input.is_open() )
	{
		cerr << endl;
		cerr << "ERROR: Could not open initial value input file " \
		     << inFilename << "." << endl;
		cerr << "       A sample input file called 'sample_input' was " \
             << "created.\a\n" << endl;

		ofstream new_input("sample_input",ios::out);
		new_input << 0.5 << '\t' << 0.601 << endl;
		new_input.close();
		exit(0);
	}
	input >> x[0];
	input >> y[0];
	input.close();
    if (verbose)
        cout << "Input initial values from file" << endl;

	/**********************************************************/
	// Print summary of input values
	if (not silent)
	{
		cout << endl;
		cout << argv[0] << " summary of inputs: \n" << endl;
		cout << "\tkicking strength (k) = \t " << k*2.0*M_PI << endl;
		cout << "\trp threshold (e) = \t " << ge << endl;
		cout << "\twindow size (w) = \t " << globalWindow << endl;
		cout << "\twindow overlap (n) = \t " << globalOverlap << endl;
		cout << "\ttrajectory length (l) =  " << std::scientific << static_cast<double>(l) << std::fixed << endl;
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
	stdmpInit(x,y);                         //** Copy code from stdmp.cpp to here
	rrCurr = rrInit(x,y,rrcntr);            //** Copy code from rr.cpp to here

	// If begin in sticky event, exit event
	while ( rrCurr > thr && winNumMax != currWin )
	{
        /* standard map section */
        // stdmp copy from previous run
        for ( ii = diff; ii < globalWindow; ii += 2 )
        { // copies n overlap values -- assumes (w-n) % 2 == 0
            x[ii - diff] = x[ii];
            y[ii - diff] = y[ii];
            x[ii - diff + 1] = x[ii + 1];
            y[ii - diff + 1] = y[ii + 1];
        }

        // stdmp calculated new values
        for ( ii = globalOverlap; ii < globalWindow; ii+=2 )
        {	// calculates new non-overlapping values
            y[ii] = fmod(y[ii-1] + k * sin( TWOPI * x[ii-1] ) + modCorr, 1.0);
            x[ii] = fmod(y[ii] + x[ii-1] + modCorr, 1.0);
            y[ii+1] = fmod(y[ii] + k*sin( TWOPI * x[ii] ) + modCorr, 1.0);
            x[ii+1] = fmod(y[ii+1] + x[ii] + modCorr, 1.0);
        }

        /* recurrence rate section */
        // Region 1
        RR = rrcntr;
        rrcntr = 0;

        // Region 2
//		rrCurr = rr(x,y,rrcntr);
        for (ii = 0; ii < 50; ii++)
        {
            for (jj = 50; jj < 100; jj++)
            {
                // Can these if statements be amenable to #pragma omp sections?
                //** Are you confident this is sufficient?  (should i and j be swapped?)
                //** Perhaps it is because we're only considering 1/2 the triangle?
                //** Perhaps it's not because we're leaving out x[i] < ge && x[j] > (1.0 - ge)?

                // Calculate recurrences; check if they are on opposite edges
                if ( x[ii] > minusge && ptrx[jj] < ge )
                    dx = fabs(1.0 - x[ii] + ptrx[jj]);
                else
                    dx = fabs(x[ii] - ptrx[jj]);

                if ( y[ii] > minusge && ptry[jj] < ge )
                    dy = fabs(1.0 - y[ii] + ptry[jj]);
                else
                    dy = fabs(y[ii] - ptry[jj]);

                // Maximum Norm
                dx > dy ? (dmax = dx) : (dmax = dy);

                // Apply Threshold
                if (dmax < ge)
                    RR++;
            }
        }

        // Region 3
        for (ii = 50; ii < 99; ii++)
        {
            for (jj = ii+1; jj < 100; jj++)
            {
                if ( x[ii] > minusge && ptrx[jj] < ge )
                    dx = fabs(1.0 - x[ii] + ptrx[jj]);
                else
                    dx = fabs(x[ii] - ptrx[jj]);

                if ( y[ii] > minusge && ptry[jj] < ge )
                    dy = fabs(1.0 - y[ii] + ptry[jj]);
                else
                    dy = fabs(y[ii] - ptry[jj]);

                // Maximum Norm
                dx > dy ? (dmax = dx) : (dmax = dy);

                // Apply Threshold
                if (dmax < ge)
                    rrcntr++;
            }
        }
        
        RR += rrcntr;
        rrCurr = static_cast<double>(2 * RR + globalWindow) / wsquared;

        currWin++;
	}

	cout << "Exited existing sticky event of length " << currWin << "." << endl;

    if ( winNumMax == currWin )
	{
		cerr << "\n\aError:\tThreshold (r) larger than largest recurrence " \
             << "rate value.  No data output.\n" \
		     << "\tConsider a smaller rr threshold (r) or larger data set (l)." \
		     << endl;

		publish = false;
		exit(0);
	}

    while ( currWin < winNumMax )
	{
        /* standard map section */
        // stdmp copy from previous run
        for ( ii = diff; ii < globalWindow; ii += 2 )
        { // copies n overlap values -- assumes (w-n) % 2 = 0
            x[ii - diff] = x[ii];
            y[ii - diff] = y[ii];
            x[ii - diff + 1] = x[ii + 1];
            y[ii - diff + 1] = y[ii + 1];
        }
        // stdmp calculated new values
        for ( ii = globalOverlap; ii < globalWindow; ii+=2 )
        {	// calculates new non-overlapping values
            y[ii] = fmod(y[ii-1] + k * sin( TWOPI * x[ii-1] ) + modCorr, 1.0);
            x[ii] = fmod(y[ii] + x[ii-1] + modCorr, 1.0);
            y[ii+1] = fmod(y[ii] + k*sin( TWOPI * x[ii] ) + modCorr, 1.0);
            x[ii+1] = fmod(y[ii+1] + x[ii] + modCorr, 1.0);
        }

        /* recurrence rate section */
        rrLast = rrCurr;

        // Region 1
        RR = rrcntr;
        rrcntr = 0;

        // Region 2
        for (ii = 0; ii < 50; ii++)
        {
            for (jj = 50; jj < 100; jj++)
            {
                // Can these if statements be amenable to #pragma omp sections?
                //** Are you confident this is sufficient?  (should i and j be swapped?)
                //** Perhaps it is because we're only considering 1/2 the triangle?
                //** Perhaps it's not because we're leaving out x[i] < ge && x[j] > (1.0 - ge)?

                // Calculate recurrences; check if they are on opposite edges
                if ( x[ii] > minusge && ptrx[jj] < ge )
                    dx = fabs(1.0 - x[ii] + ptrx[jj]);
                else
                    dx = fabs(x[ii] - ptrx[jj]);

                if ( y[ii] > minusge && ptry[jj] < ge )
                    dy = fabs(1.0 - y[ii] + ptry[jj]);
                else
                    dy = fabs(y[ii] - ptry[jj]);

                // Maximum Norm
                dx > dy ? (dmax = dx) : (dmax = dy);

                // Apply Threshold
                if (dmax < ge)
                    RR++;
            }
        }

        // Region 3
        for (ii = 50; ii < 99; ii++)
        {
            for (jj = ii+1; jj < 100; jj++)
            {
                if ( x[ii] > minusge && ptrx[jj] < ge )
                    dx = fabs(1.0 - x[ii] + ptrx[jj]);
                else
                    dx = fabs(x[ii] - ptrx[jj]);

                if ( y[ii] > minusge && ptry[jj] < ge )
                    dy = fabs(1.0 - y[ii] + ptry[jj]);
                else
                    dy = fabs(y[ii] - ptry[jj]);

                // Maximum Norm
                dx > dy ? (dmax = dx) : (dmax = dy);

                // Apply Threshold
                if (dmax < ge)
                    rrcntr++;
            }
        }

        RR += rrcntr;
        rrCurr = static_cast<double>(2 * RR + globalWindow) / wsquared;

		if (rrCurr > thr && rrLast <= thr)
		{
            rrEnter = currWin;
		}
		else if (rrCurr <= thr && rrLast > thr)
		{
            tau = currWin - rrEnter;
			gsl_histogram_increment(h,static_cast<double> (diff * tau + globalOverlap));
		}
        currWin++;
	}

    tau = winNumMax - rrEnter;

	cout << "Uncompleted event of length " << tau << "." << endl;
	cout.precision(15);
	cout.setf(ios::floatfield);
	cout << "Final (x,y) = (" << x[globalWindow-1] << ", " << y[globalWindow-1] << ")" << endl;

//	if (false)
//	{
//		char buff[50];
//		// get the thread number from inFilename, increment it by one
//		// then put it in %d below
//		int temp = sprintf(buff,"%s_%d",inFilename,1);
//		ofstream strand(buff,ios::out);
//		strand << x[globalWindow-1] << '\t' << y[globalWindow-1] << endl;
//		strand.close();
//	}

	/* Output the preceding values into a file called inputX_n.txt, where
	 the X refers to the strain number, and the n refers to the segment of the
	 strain that the initial conditions represent the starting conditions of, 
	 therefore, subsequent calls can load the histogram and keep adding values 
	 to it and your small nightly runs can start to accumulate into hugely big
	 data series with much more reliable statistics.  Appeal to HDF5 files to
	 add this support.  When I do this, though, I need to get rid of the code
	 below that removes the bins with no counts.
	 */

	/**********************************************************/
	/*       Construct CDF                                    */

    cout << "Constructing CDF." << endl;
    double histSum = gsl_histogram_sum(h);
    vector<int> binTimes;

    if (verbose)
        printHistogram(h);

    // Find non-zero bins and store their indices in the vector binTimes[]
    for ( int ii = 0; ii < nBins; ii++ )
        if ( gsl_histogram_get(h,ii) )	// If non-zero entry, put in binTimes
            binTimes.push_back(ii);

    if (verbose)
        cout << "binTimes.size() = " << binTimes.size() << endl;

    // Create vectors to store the x and y histogram values that are
    //   modified and copied from the histogram information
    double xTimes[binTimes.size()], yCount[binTimes.size()];

    for ( int j = binTimes.size()-1; j >= 0; j-- )
    {
        yCount[j] = gsl_histogram_get(h,binTimes[j]);
        xTimes[j] = ( range[binTimes[j]] + range[binTimes[j]+1] ) / 2.0;
    }

    if (verbose)
        printXY(xTimes,yCount,binTimes.size());

    // Accumulate y component for CDF
    // --- Look into the accumulate function that's in either numerics or algorithm?
    cumSumNorm(yCount,binTimes.size(),histSum);

    // Print to file
    fid.open(outFilename,ios::app);
    if (fid.is_open())
    {
        fid << asctime(timeinfo);
        for (usInt i = 0; i<binTimes.size(); i++)
            fid << log10(xTimes[i]) << '\t' << log10(yCount[i]) << endl;
        fid << "\n\n";
        fid.close();
    }
    else
    {
        cerr << "Error: Could not open " << outFilename << endl;
        cout << "\tPrinting to screen..." << endl;
        for (usInt i = 0; i<binTimes.size(); i++)
            cout << log10(xTimes[i]) << '\t' << log10(yCount[i]) << endl;
    }
/**********************************************************/
/*    Perform least square linear fit if chosen           */
    if (doFit)
    {
        for (usInt i = 0; i < binTimes.size(); i++)
        {
            xTimes[i] = log10(xTimes[i]);
            yCount[i] = log10(yCount[i]);
        }
        if (verbose)
            printXY(xTimes,yCount,binTimes.size());

        double b, m1, varb, covbm, varm, sumsq;
        gsl_fit_linear(xTimes,1,yCount,1, binTimes.size(), &b, &m1, \
                       &varb, &covbm, &varm, &sumsq);
        double m2 = -1.0 * mlefit(xTimes,binTimes.size());
        cout.precision(5);

        cout.setf(ios::fixed,ios::floatfield);
        cout << "+--------------------------------------------------------+" << endl;
        cout << "| Ordinary log-log fit:                                  |" << endl;
        cout << "| \t(" << m1 << " +/- " << sqrt(varm) << ")x + (" << b << " +/- " << sqrt(varb) << " ) |" << endl;
        cout << "| MLE fit:                                               |" << endl;
        cout << "| \t(" << m2 << " +/- " << (m2-1.0)/sqrt(binTimes.size()) << ")x + " << 0.0 << "               |" << endl; 
        cout << "+--------------------------------------------------------+" << endl;
    }

	/**********************************************************/
	// End of program

	t2 = time(NULL) - t1;
	t1 = t2;
	if (!silent) timesummary(t1,t2);

	//Summarize results into log
	if (publish)
	{
		double xInit = 0.0, yInit = 0.0;
		input.open(inFilename,ios::in);
		input >> xInit;
		input >> yInit;
		input.close();
		
		ofstream resultsLog("rrstdmp_results.log",ios::app);
		if (resultsLog.is_open())
		{
			resultsLog << k*TWOPI << "," << ge << "," << l << "," \
                       << globalWindow << "," << globalOverlap << "," \
                       << xInit << "," << yInit << "," << t1 << "," \
                       << asctime(timeinfo);
//			if (doFit)
//				resultsLog << m1 << "," << m2; //** Find out how to add this to log file

			resultsLog.close();
			if (!silent) cout << "Input summary written to rrstdmp_results.log" << endl;
		} else
		{
			cerr << "Error: Could not open rrstdmp_results.log" << endl;
			cerr << k * TWOPI << "," << ge << "," << l << "," << globalWindow \
                 << "," << globalOverlap << "," << xInit << "," << yInit \
                 << "," << t1 << "," << asctime(timeinfo);
		}
	}

	gsl_histogram_free(h);

	return 0;
}
