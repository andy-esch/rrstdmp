// Compile with g++ -otest -lboost_unit_test_framework TestExample.cpp -L/opt/local/lib -I/opt/local/include
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_misc
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>

// Libraries to be tested
#include "misc.h"
#include "rr.h"
#include "stdmp.h"
#include "rrmean.h"


BOOST_AUTO_TEST_CASE( test1 )
{
    BOOST_CHECK_EQUAL( 2,2 );
}

BOOST_AUTO_TEST_CASE( linspacetest )
{
    double * range;
    int size = 10;
    range = new double[size];
    linspace(range, size, 0, 9); // Should make a set = {0, 1, ..., 9}
    for (int ii = 0; ii < size; ii++)
    {
        BOOST_CHECK_EQUAL( range[ii],ii );
    }
}

BOOST_AUTO_TEST_CASE( logspacetest )
{
    double * range;
    int size = 10;
    range = new double[size];
    logspace(range, size, 0, 9); // Should make a set = {1, 10, 100, 1000, ..., 10^9};
    for (int ii = 0; ii < size; ii++)
        BOOST_CHECK_EQUAL( range[ii],pow(10,ii) );

    delete [] range;
    range = NULL;

    // Test for variables acting as min, max
    double min = 0.0;
    double max = 10.0;
    range = new double[size+1];
    logspace(range,size+1,min,max);
    for (int ii = 0; ii < size; ii++)
        BOOST_CHECK_EQUAL( range[ii],pow(10,ii) );

    delete [] range;
    range = NULL;

    // Test for long int cast as double and then its log10() is taken
    long maxDouble = 1.0e9;
    range = new double[size];
    logspace(range, size, 0, log10(static_cast<double>(maxDouble)));
    for (int ii = 0; ii < size; ii++)
        BOOST_CHECK_EQUAL( range[ii], pow(10,ii) );
}

BOOST_AUTO_TEST_CASE( cumSumNormtest )
{
    double x[10], sum = 0.0;
    double ans[10] = {55.0/55.0, 54.0/55.0, 52.0/55.0, 49.0/55.0, 45.0/55.0, \
                      40.0/55.0, 34.0/55.0, 27.0/55.0, 19.0/55.0, 10.0/55.0};
    for (int ii = 0; ii < 10; ii++)
    {
        x[ii] = static_cast<double>(ii + 1);
        sum += x[ii];
    }
    cumSumNorm(x,10,sum);

    for (int ii = 0; ii < 10; ii++)
        BOOST_CHECK( fabs(ans[ii] - x[ii]) < 1.0e-10 );
}

int globalWindow = 100, globalOverlap = 50;
double ge = 0.01, k = 0.97163540631/(2.0 * M_PI);
const double TWOPI = 2.0 * M_PI;

BOOST_AUTO_TEST_CASE( global_values )
{
    BOOST_CHECK_MESSAGE(fabs(k - 0.97163540631/(2.0 * M_PI)) < 1.0e-10,
                        "k is not equal to set global value");
    BOOST_CHECK_MESSAGE(fabs(ge - 0.01) < 1.0e-10,
                        "global recurrence threshold not set correctly");
    BOOST_CHECK_MESSAGE(globalWindow == 100,
                        "globalWindow not properly set");
    BOOST_CHECK_MESSAGE(globalOverlap == 50,
                        "globalOverlap not properly set");
}

BOOST_AUTO_TEST_CASE( stdmp_period_test )
{
    int size = 100;
    double x[size], y[size], runner = 0.0;
    x[0] = 0.0; y[0] = 0.25;
    k = 0;

    stdmpInit(x,y);
    for (int ii = 0; ii < size; ii++)
    {
        runner = static_cast<double>(ii % 4) * 0.25;
        BOOST_CHECK_MESSAGE(fabs(x[ii] - runner) < 1.0e-5, \
                            "x[" << ii << "]");
        BOOST_CHECK_MESSAGE(fabs(y[ii] - 0.25) < 1.0e-5, \
                            "y[" << ii << "]");
    }

    stdmp(x,y);
    for (int ii = 0; ii < size; ii++)
    {
        runner = static_cast<double>((ii + 2) % 4) * 0.25;
        BOOST_CHECK_MESSAGE(fabs(x[ii] - runner) < 1.0e-5, \
                            "x[" << ii << "]");
        BOOST_CHECK_MESSAGE(fabs(y[ii] - 0.25) < 1.0e-5, \
                            "y[" << ii << "]");
    }
}

BOOST_AUTO_TEST_CASE( rrInit_test )
{
    int size = 100, rrcntr = 0;
    double x[size], y[size];

    // Case 1: k = 0; 4-periodic map
    k = 0;
    x[0] = 0.0; y[0] = 1.0/4.0;
    stdmpInit(x,y);
    BOOST_CHECK_EQUAL(y[0],rrInit(x,y,rrcntr));

    // Case 2: k = k_g; 2-period map
    k = 0.97163540631/(2.0 * M_PI);
    x[0] = 0.0; y[0] = 0.5;
    stdmpInit(x,y);
    BOOST_CHECK_EQUAL(y[0],rrInit(x,y,rrcntr));

    // Case 3: k = 10; randomly chosen initial conditions
    // Expectation: RR = 1% of phase space (bc of chosen threshold ge)
    k = 10.0;
    x[0] = 0.112311; y[0] = M_PI / 10.0;
    stdmpInit(x,y);
    double rrCurr = rrInit(x,y,rrcntr);
    BOOST_CHECK_MESSAGE(fabs(rrCurr - 0.01) < 1.0e-3, \
                        "rrInit 3: (rr - 0.01) = " << (rrCurr - 0.01));
}

BOOST_AUTO_TEST_CASE( rr_test )
{
    int size = 100, rrcntr = 0;
    double x[size], y[size], rrCurr = 0.0;

    // case 1: k = 0; 4-periodic map
    // RR should = 1200 [sum(100 - 4 * n)_n=1->(100/4) = 1200]
    k = 0;
    x[0] = 0.0; y[0] = 1.0 / 4.0;
    stdmpInit(x,y);
    rrCurr = rrInit(x,y,rrcntr);
    BOOST_CHECK_MESSAGE(fabs(rrCurr - 0.25) < 1.0e-5, \
                        "0: (rr - 0.25) = " << (rrCurr - 0.25));

    for (int ii = 0; ii < 10; ii++)
    {
        stdmp(x,y);
        rrCurr = rr(x,y,rrcntr);
        BOOST_CHECK_MESSAGE(fabs(rrCurr - 0.25) < 1.0e-4, \
                            "Case 1." << (ii+1) << ": (rr - 0.25) = " << (rrCurr - 0.25));
    }

    // case 2: k = 10; stochastic(?) -> RR_avg =~ recurrence threshold (ge)
    k = 10;
    x[0] = M_PI / 10.0; y[0] = M_E / 10.0;
    stdmpInit(x,y);
    rrCurr = rrInit(x,y,rrcntr);
    BOOST_CHECK_MESSAGE(fabs(rrCurr - 0.01) < 1.0e-3, \
                        "0: (rr - 0.01) = " << (rrCurr - 0.01));
    for (int ii = 0; ii < 10; ii++)
    {
        stdmp(x,y);
        rrCurr = rr(x,y,rrcntr);
        BOOST_CHECK_MESSAGE(fabs(rrCurr - 0.01) < 1.0e-3, \
                            (ii+1) << ": (rr - 0.01) = " << (rrCurr - 0.01));
    }
}

// Compare the copied values between various runs of stdmp()
BOOST_AUTO_TEST_CASE( stdmpTester )
{
    k = 0.97163540631/(2.0 * M_PI);
    double x1[100], y1[100];
    double x2[100], y2[100];
    x1[0] = 0.05; y1[0] = 0.0625;
    x2[0] = 0.05; y2[0] = 0.0625;

    stdmpInit(x1,y1);
    stdmpInit(x2,y2);
    stdmp(x2,y2);

    for (int ii = 0; ii < 50; ii++)
    {
        BOOST_CHECK_MESSAGE(fabs(x1[ii+50] - x2[ii]) < 1.0e-5, \
                            "Value x[" << ii << "] did not match up\n");
        BOOST_CHECK_MESSAGE(fabs(y1[ii+50] - y2[ii]) < 1.0e-5, \
                            "Value y[" << ii << "] did not match up\n");
    }
}





//EOF
