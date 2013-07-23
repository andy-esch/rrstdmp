// Compile with g++ -otest -lboost_unit_test_framework TestExample.cpp -L/opt/local/lib -I/opt/local/include
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rrstdmp test
#include <boost/test/unit_test.hpp>
#include <cmath>

// Library to be tested
#include "misc.h"


BOOST_AUTO_TEST_CASE( test1 )
{
    BOOST_CHECK( 2==2 );
}

BOOST_AUTO_TEST_CASE( test2 )
{
    BOOST_CHECK( 3==1 );
}

BOOST_AUTO_TEST_CASE( linspacetest )
{
    double * range;
    int size = 10;
    range = new double[size];
    linspace(range, size, 0, 9); // Should make a set = {0, 1, ..., 9}
    for (int ii = 0; ii < size; ii++)
    {
        BOOST_CHECK( range[ii] == ii );
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
        BOOST_CHECK( range[ii] == pow(10,ii) );
    }
}

//EOF
