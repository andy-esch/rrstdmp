// Compile with g++ -otest -lboost_unit_test_framework TestExample.cpp -L/opt/local/lib -I/opt/local/include
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rrstdmp_test
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
    k = 10.0;
    x[0] = 0.112311; y[0] = M_PI / 10.0;
    stdmpInit(x,y);
    double rrCurr = rrInit(x,y,rrcntr);
    BOOST_CHECK_MESSAGE(fabs(rrCurr - 0.01) < 1.0e-10, \
                        "(rr - 0.01) = " << (rrCurr - 0.01));
}

BOOST_AUTO_TEST_CASE( rr_test )
{
    int size = 100, rrcntr = 0;
    double x[size], y[size], rrCurr;

    // case 1: k = 0; 4-periodic map
    k = 0;
    x[0] = 0.0; y[0] = 1.0 / 4.0;
    stdmpInit(x,y);
    rrCurr = rrInit(x,y,rrcntr);
    BOOST_CHECK_MESSAGE(fabs(rrCurr - 0.25) < 1.0e-5, \
                        "(rr - 0.01) = " << (rrCurr - 0.01));

    for (int ii = 0; ii < 10; ii++)
    {
        stdmp(x,y);
        rrCurr = rr(x,y,rrcntr);
        BOOST_CHECK_MESSAGE(fabs(rrCurr - 0.25) < 1.0e-5, \
                            "(rr - 0.01) = " << (rrCurr - 0.01));
    }
}

//EOF
