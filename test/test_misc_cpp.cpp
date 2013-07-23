// Compile with g++ -otest -lboost_unit_test_framework TestExample.cpp -L/opt/local/lib -I/opt/local/include
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rrstdmp test
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>

// Library to be tested
#include "misc.h"


BOOST_AUTO_TEST_CASE( test1 )
{
    BOOST_CHECK_EQUAL( 2,2 );
}

BOOST_AUTO_TEST_CASE( test2 )
{
    BOOST_CHECK_EQUAL( 3,1 );
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
    {
        BOOST_CHECK_EQUAL( range[ii],pow(10,ii) );
    }
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
        BOOST_CHECK( abs(ans[ii] - x[ii]) < 1.0e-20 ); // Is this a sufficient comdition?
}

//EOF
