// Compile with g++ -otest -lboost_unit_test_framework TestExample.cpp -L/opt/local/lib -I/opt/local/include
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Hint test
#include <boost/test/unit_test.hpp>
#include "testlib.h"


BOOST_AUTO_TEST_CASE( test1 )
{
    BOOST_CHECK( 2==2 );
}

BOOST_AUTO_TEST_CASE( test2 )
{
    BOOST_CHECK( 3==1 );
}

BOOST_AUTO_TEST_CASE( test3 )
{
    int num1 = 3, num2 = num1 - 7;
    BOOST_CHECK( sum(num1, num2)==7 );
}

BOOST_AUTO_TEST_CASE( linspacetest )
{
    double * range;
    int size = 10;
    range = new double[size];
    linspace(range, size, 0, 9); // Should make a set = {0, 1, ..., 9}
    for (int ii = 0; ii < 10; ii++)
    {
        BOOST_CHECK( range[ii] == ii );
    }
}

//EOF
