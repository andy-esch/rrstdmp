
#ifndef RR_H
#define RR_H

#include <cmath> // fabs()
#include <omp.h>

double rr(double *__restrict__ x, double *__restrict__ y, int &__restrict__ rrcntr);
/* { 
 	extern int globalWindow, globalOverlap;
	extern double ge;
	double dx, dy, dmax;
	int RR = rrcntr;	// Transfer previous overlap into current count
	rrcntr = 0;			// Reset overlap for this window
                        //	int rrtemp = 0;
	static int diff = globalWindow - globalOverlap;
	int i, j;
	static double wsquared = static_cast<double> (globalWindow*globalWindow);
    
        //#pragma omp parallel for default(none) private(dx,dy,dmax,i,j,chx,chy) \
         //shared(x,y,diff) schedule(guided) reduction(+:RR,rrtemp)
	for (i = globalOverlap; i < globalWindow; i++)
	{
		for (j = 0; j < i; j++)
		{
            // Calculate recurrences; check if they are on opposite edges
			if ( x[i] > (1.0 - ge) && x[j] < ge )
				dx = fabs(1.0 - x[i] + x[j]);
			else
				dx = fabs(x[i] - x[j]);
            
			if ( y[i] > (1.0 - ge) && y[j] < ge )
				dy = fabs(1.0 - y[i] + y[j]);
			else
				dy = fabs(y[i] - y[j]);
            
                // Maximum Norm
			dx > dy ? (dmax = dx) : (dmax = dy);
            
                // Apply Threshold
			if (dmax < ge)
			{
				RR++;
				if ( (i > diff) && (j >= diff) )
					rrcntr++;
			}
		}
	} // *** End of Parallel Section **
    
        //	rrcntr = rrtemp;
    
	return ( double(2 * RR + globalWindow) / wsquared );
} */

double rrInit(double*, double*, int&);

#endif // RR_H
