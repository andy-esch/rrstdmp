/*
 *  cmdLineInput.cpp
 *  
 *  Description:
 *
 *
 *  Created by Peter Eschbacher on 11/10/11.
 *
 */

#include "cmdLineInput.h"

void cmdLineInput(int argc, char **argv, char * inFilename, char * outFilename, \
                  long & l, double & thr, int & nBins, \
                  bool & silent, bool & publish, bool & doFit, bool & rrtest)
{
    const char* optstring(":i:o:e:w:n:l:k:r:b:aspfth?");
    extern int optopt;
    double temp2 = 0.0;

    if ( argc == 1 )
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
                        temp2 = atof(optarg);
                        l = static_cast<long>(temp2);
                        break;
                    case 'k':
                        k = atof(optarg)/TWOPI;
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
                        if (nBins <= 1)
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
                        cout << "\aOption '-" << static_cast<char> (optopt) << "' is not valid." << endl;
                    case 'h':
                        usage_rrstdmp(ge,k,l,globalOverlap,globalWindow,thr,nBins);
                        exit(0);
                    default:
                        usage_rrstdmp(ge,k,l,globalOverlap,globalWindow,thr,nBins);
                        exit(0);
                }
            }
        }
}