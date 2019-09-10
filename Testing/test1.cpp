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

    // Compare varying values of u, for a = 0.5, b = 0.5
    int length = 5;
    std::vector<float> z_test_values = linspace((float) 1., (float) 10000., 5);
    for(int i = 0; i < length; i++)
    {
        z_test_values[i] = 1/(z_test_values[i]*z_test_values[i]);
    }
    z_test_values.push_back(0.); // Also test when zero
    length++;
    // 'correct' results from Wolfram Alpha.
    std::vector<float> res_wolfram{3.141592653589793, 0.000799759985303, 0.00039996000066, 0.000266657833977, 0.000200000000333, 0.0};
    float res;
    for(int i = 0; i < length; i++)
    {
        res = incompleteBeta(0.5, 0.5, z_test_values[i]);
        std::cout << "iteration:" << i << " " << res;
        BOOST_CHECK(AreSame(res, res_wolfram[i]));
    }


    //float ib = incompleteBeta(0.5, 0.5, 0.1);
    //float ib_wolfram = 0.643501;
    //BOOST_CHECK(AreSame(ib, ib_wolfram));

    //ib = incompleteBeta(0.5, 0.5, 0.1);


}