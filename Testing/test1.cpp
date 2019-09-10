#define BOOST_TEST_MODULE My Test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "vector"
#include "utillity.h"
#include "stdlib.h"
#include "stdio.h"


BOOST_AUTO_TEST_CASE(Linspace)
{
    float start = 0.;
    float stop = 1.;
    int steps = 6;
    std::vector<float> simple;
    simple = linspace(start, stop, steps);

    std::vector<float> should{0., 0.2, 0.4, 0.6, 0.8, 1.};

    for(int i = 0; i < steps; i++)
    {
        BOOST_CHECK(AreSame(simple[i], should[i]));
    }

    start = -1;
    simple.clear();
    simple = linspace(start, stop, steps);
    should.clear();
    should = {-1., -0.6, -0.2, 0.2, 0.6, 1.};

    for(int i = 0; i < steps; i++)
    {
        BOOST_CHECK(AreSame(simple[i], should[i]));
    }
}

BOOST_AUTO_TEST_CASE(IBeta)
{
    // All comparisons are from Wolfram Alpha

    float ib = incompleteBeta(0.5, 0.5, 0.1);
    float ib_wolfram = 0.643501;
    BOOST_CHECK(AreSame(ib, ib_wolfram));

    ib = incompleteBeta(0.5, 0.5, 0.1);


}