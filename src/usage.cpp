
/* usage.cpp -- This program prints the man file-like information for the use of 
                  rrstdmp, as well as printing default values
 */

#include "usage.h"
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>

using namespace std;

void version(void)
{
	string version("$Revision: 2.8.0 $($Date: 2011/05/30 14:04:00 $)");

	int pos = 0;	// remove "$"

	while ( pos >= 0  && pos < version.length())
	{
		pos = version.find("$");
		if ( pos >=0 && pos < version.length())
			version.erase(pos, 1);
	}
	version.erase(0, 10); // remove "Revision:"
	pos = version.find("Date");
	version.erase(pos, 6); // remove "Date:"
	pos = version.find(" )");
	version.erase(pos, 1); // remove space before ")"

	cerr << "Commandline rr calculator, version " << version << endl;
}

void usage_rrstdmp(const double e, const double k, mpz_t l, \
				   const unsigned int n, const unsigned int w, \
				   const double thr, const int b)
{
	cerr << endl;
	version();
	cerr << "  Calculates RR for the Standard Map from given initial \n";
	cerr << "  conditions and kicking parameter.  Output is in the form\n";
	cerr << "  of a log-log histogram of sticking events and linear fits\n";
	cerr << "  for the histogram distribution based on two methods (MLE and\n";
	cerr << "  a standard, though unreliable, least-square fit).\n";
	cerr << endl;
	cerr << "  (c) Peter A. Eschbacher, otto.hasselblad@gmail.com\n";
	cerr << "  Some code borrowed and modified from Norbert Marwan's \n";
	cerr << "  recurrence rate calculator (private correspondence)\n\n";
	cerr << "usage:\n  rrstdmp [options]\n" << endl;
	cerr << "options:" << endl;
	cerr << "    -i <string>   initial values file (input filename)" << endl;
	cerr << "    -o <string>   output file (output filename)" << endl;
	cerr << "    -e <number>   threshold, default = " << e << " (" << e*e*400 << "% of total phase space)" << endl;
	cerr << "    -w <number>   window, default = " << w << endl;
	cerr << "    -n <number>   window overlap, default = " << n << endl;
	cerr << "    -l <number>   length of total stdmap trajectory to be analyzed, default = " << l << endl;
	cerr << "    -k <number>   nonlinearity parameter strength, default = " << k*6.283185307179587 << endl;
	cerr << "    -r <number>   rr sticking event threshold" << endl;
	cerr << "                    If no value is specified, one will be calculated (highly recommended)." << endl;
	cerr << "                    Value is based on orbit closest approximating golden KAM curve for " << endl;
	cerr << "                      the w, n, and k given." << endl;
	cerr << "    -b <number>   number of histogram bins, default = " << b << endl;
	cerr << "    -p            turn off publishing of inputs summary to results.log" << endl;
	cerr << "    -f            perform fits to log-log plot of CDF" << endl;
	cerr << "    -t            perform a calibration test of the recurrence rate calculator" << endl;
	cerr << "    -s            silences output, default off" << endl;
	cerr << "    -h            print this help text" << endl;
	cerr << endl;
}

void usage_rrmapgen(const double e, const double k, mpz_t l, \
					const unsigned int n, const unsigned int w, \
					double thr, const int b)
{
	cerr << endl;
	version();
	cerr << "  Calculates the recurrence rate for the Standard Map from given "
	"initial conditions and parameter k.  Outputs trajectories of StMap for "
	"processing between the values specified with the q flag." << endl;
	cerr << "  (c) Peter A. Eschbacher, otto.hasselblad@gmail.com" << endl;
	cerr << "  Some code borrowed and modified from Norbert Marwan's "
	"recurrence rate calculator (private correspondence)" << endl;
	cerr << endl;
	cerr << "usage:\n  rrmapgen [options]\n\n";
	cerr << "options:" << endl;
	cerr << "    -i <string>      Initial values file (input)" << endl;
	cerr << "    -e <double>      RP threshold, default = " << e << endl;
	cerr << "    -w <integer>     RP window size, default = " << w << endl;
	cerr << "    -n <integer>     RP window overlap, default = " << n << endl;
	cerr << "    -l <integer>     Length of overall stdmap trajectory to be analyzed, default = " << l << endl;
	cerr << "    -k <double>      Standard map parameter strength, default = " << k*6.283185307179587 << endl;
	cerr << "    -r <double>      RR sticking event threshold (not recommended)" << endl;
	cerr << "                       If no value is specified, one will be calculated." << endl;
	cerr << "                       Value is based on golden KAM curve for the w, n, and k given." << endl;
	cerr << "    -b <integer>     Number of histogram bins, default = " << b << endl;
	cerr << "    -p               Turn off publishing of inputs summary to results.log" << endl;
	cerr << "    -f               Perform linear fit to log-log plot of CDF and MLE" << endl;
	cerr << "    -q <integer>     Sticky event values to constrain output of standard map sticking events" << endl;
	cerr << "    -s               Silences stdout (non-verbose), default off" << endl;
	cerr << "    -h               Print this wonderfully informative text" << endl;
	cerr << endl;
}

/*

void changeLog(void)
{	
	cerr << "\n Change Log -- ";
	version();
	cout << "***************************************************************************" << endl;
	cout << " December 15\n";
	cout << "   Full histogram functionality and CDF construction.  Still not 100% confidently implemented";
	cout << "     but pretty much there, even matching gamma values as seen in rrslope.m\n";
	cout << "     Length is now output into results.log in scientific notation, and rrcdf(...)\n";
	cout << "     contains all of the CDF information, where the hist and other variables are passed\n";
	cout << "     in the function.\n";
	cout << " December 13\n";
	cout << "   Histogram functionality, though still in beta testing\n";
	cout << " November 26\n";
	cout << "  Set output precision in rrstdmp to 15 significant figures.  version() became\n";
	cout << "    a function for use in this (changeLog()) and usage(...).  the make file\n";
	cout << "    has been made better, including linkage to GNU Scientific Library,\n";
	cout << "    clean, and stage functionality.\n";
	cout << " November 20\n";
	cout << "  The program is now correctly calibrated (detects n-periods with very good\n";
	cout << "    accuracy).  Reminder: rrmean is output at 98% of calculated value of the \n";
	cout << "    golden KAM curve.\n";
	cout << " November 13/14\n";
	cout << "  Added overlap copying in the rr and rrInit functions with rrcntr.  This \n";
	cout << "    allows far fewer computations and an improved performance.  There is a bug\n";
	cout << "    in rrmean or rr that is dependent on n and outputs the incorrect rrmean if\n";
	cout << "    n is more than half of w.\n";
	cout << " November 12\n";
	cout << "  Added rrmean.cpp to the makefile.  This calculates the recurrence rate of\n";
	cout << "    the golden KAM curve for the w and e parameters enetered if no rr threshold\n";
	cout << "    is entered (the r flag).  The -p flag was also added -- it allows the\n";
	cout << "    option of publishing run details to file.\n";
	cout << " November 10\n";
	cout << "  Added rrmean output to the same file that contains the sticky event length.\n";
	cout << "    The output file size is around 13.6MB per billion trajectory points, but \n";
	cout << "    but this is dependent upon the data analyzed.\n";
	cout << " November 8\n";
	cout << "  Added the sticky threshold kernel.  This now reduces the total compute time.\n";
	cout << " November 4\n";
	cout << "  Changed rr module so that the output is (2*RR+w)/w*w, to give the recurrence \n";
	cout << "    rate  of the full map.  This allows n-periods to be exactly 1/n instead of \n";
	cout << "    an  approximate.  The std module was updated to account for a minor error \n";
	cout << "    in the overlapping of data: a standard error is given if the overlap is \n";
	cout << "    smaller than one.\n";
	cout << " October 28\n";
	cout << "  Added runtime information (time_t and structs) for output and terminal display\n";
	cout << " October 26\n";
	cout << "	Change atoi to atof in k input (wrong type before)\n";
} */
