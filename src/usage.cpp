
/* usage.cpp -- This program prints the man file-like information for the use of 
                  rrstdmp, as well as printing default values
 */

#include "usage.h"

void version(void)
{
	string version("$Revision: 2.8.1 $($Date: 2013/07/23 10:58:00 $)");

	unsigned short pos = 0;	// remove "$"

	while ( pos >= 0  && pos < version.length())
	{
		pos = version.find("$");
		if ( pos >= 0 && pos < version.length())
			version.erase(pos, 1);
	}
	version.erase(0, 10); // remove "Revision:"
	pos = version.find("Date");
	version.erase(pos, 6); // remove "Date:"
	pos = version.find(" )");
	version.erase(pos, 1); // remove space before ")"

	cerr << "Commandline rr calculator, version " << version << endl;
}

void usage_rrstdmp(const double e, const double k, long l, \
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
	cerr << "  recurrence rate calculator (private correspondence)\n";
    cerr << endl;
    cerr << "  View project at http://github.com/ohasselblad/rrstdmp\n\n";
    cerr << endl;
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
	cerr << "    -s            silences output, default off" << endl;
	cerr << "    -h            print this help text" << endl;
	cerr << endl;
}

void usage_rrmapgen(const double e, const double k, long l, \
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
