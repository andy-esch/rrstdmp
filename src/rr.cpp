/*
 *   This program calculates the recurrence rate of the vectors x[] and y[],
 *     dependent upon the window size w, overlap n, and recurrence threshold e.
 *     Overlap values are stored in rrcntr to reduce redundant computations for
 *     adjacent recurrence plot computations.
 *
 *   Output: recurrence rate (sparsity of the recurrence plot) of the 2-d map (x,y)
 *
 *
 *   Limitation: if n is larger than w/2, the rrmean is incorrect because the 
 *     overlap overlaps with more than one window -- so a different scheme
 *     for capturing the information from other windows needs to be devised.
 *     The program works perfectly well for n <= w/2+1, though.
 */

#include "rr.h"

/*  #include <omp.h>
 *  Currently in serial form (no OpenMP).
 *  To return to OpenMP form, uncomment the following:
 *		#include <omp.h>
 *		#pragma omp parallel...
 *		int rrtemp = 0;
 *
 *	Change rrcntr++ in the if statement in "Applying Threshold" to rrtemp++ instead
 */

double rr(double *__restrict__ x, double *__restrict__ y, \
		  int &__restrict__ rrcntr)
{
	extern int globalWindow, globalOverlap;
	extern double ge;
	double dx, dy, dmax, minusge = 1.0 - ge;
	int RR = rrcntr;	// Transfer previous overlap into current count
	rrcntr = 0;			// Reset overlap for this window
//	int rrtemp = 0;
	static int diff = globalWindow - globalOverlap;
	int i, j;
	static double wsquared = static_cast<double> (globalWindow*globalWindow);

//#pragma omp parallel for private(i,j,dx,dy,dmax) \
                         shared(x,y,ge,minusge,diff,RR,rrcntr) \
                         schedule(guided) \
                         reduction(+:RR)

// Region 1 -- Copied from previous region (copy rrcntr into RR?)

// Regions 2/3
    for (i = 0; i < 50; i++)
	{
		for (j = 50; j < 100; j++)
		{
			// Calculate recurrences; check if they are on opposite edges   //** Are you confident this is sufficient?  (should i and j be swapped?)
                                                                            //** Perhaps it is because we're only considering 1/2 the triangle?
                                                                            //** Perhaps it's not because we're leaving out x[i] < ge && x[j] > (1.0 - ge)?
            // Can these if statements be amenable to #pragma omp sections?
			if ( x[i] > minusge && x[j] < ge )
				dx = fabs(1.0 - x[i] + x[j]);
			else
				dx = fabs(x[i] - x[j]);

			if ( (y[i] > minusge && y[j] < ge) )
				dy = fabs(1.0 - y[i] + y[j]);
			else
				dy = fabs(y[i] - y[j]);

			// Maximum Norm
			dx > dy ? (dmax = dx) : (dmax = dy);

			// Apply Threshold
			if (dmax < ge)
				RR++;
		}
	}


// Region 4
    for (i = 50; i < 100; i++)
    {
        for (j = 51; j < i; j++)
        {
            if ( (x[i] > minusge && x[j] < ge) )
				dx = fabs(1.0 - x[i] + x[j]);
			else
				dx = fabs(x[i] - x[j]);

			if ( (y[i] > minusge && y[j] < ge) )
				dy = fabs(1.0 - y[i] + y[j]);
			else
				dy = fabs(y[i] - y[j]);
        }

        // Maximum Norm
        dx > dy ? (dmax = dx) : (dmax = dy);

        // Apply Threshold
        if (dmax < ge)
            RR++;
    }

	return ( double(2*RR + globalWindow) / wsquared );
    /* many calculations could be avoided if threshold is redefined as renormed 
     value = (thrsOld * wsquared - globalWindow) / 2.0
     If 10^6 runs through rr are given, then 4x as many calcs can be avoided
     (including the double() cast) and rr() could return an int...
     While this may speed up the program some, it would also be harder to under-
     stand.
     */
}


double rrInit(double *__restrict__ x, double *__restrict__ y, \
			  int &__restrict__ rrcntr) 
{
	extern int globalWindow, globalOverlap;
	extern double ge;
	double dx, dy, dmax;
	int RR = 0;
	rrcntr = 0;
	int diff = globalWindow - globalOverlap;

	for (int i = 1; i < globalWindow; i++) 
	{
		for (int j = 0; j < i; j++)
		{
			if ( x[i] > (1.0 - ge) && x[j] < ge )
				dx = fabs(1.0 - x[i] + x[j]);
			else
				dx = fabs(x[i] - x[j]);
			
			if ( y[i] > (1.0 - ge) && y[j] < ge )
				dy = fabs(1.0 - y[i] + y[j]);
			else
				dy = fabs(y[i] - y[j]);

			//maximum norm
			dx > dy ? (dmax = dx) : (dmax = dy);

			// apply threshold
			if (dmax < ge) 
			{
				RR++;
				if ( (i > diff) && (j >= diff) )
					rrcntr++;
			}
		}
	}

	return (double(2*RR+globalWindow)/double(globalWindow*globalWindow));
}
