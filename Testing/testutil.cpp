#include "testutil.h"

void test_util(void)
{
   test_linspace();
   test_ibeta();
}

void test_linspace(void)
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

void test_ibeta(void)
{
    // All comparisons are from Wolfram Alpha

    // Compare varying values of z, for a = 0.5, b = 0.5
    int length = 5;
    std::vector<float> z_test_values = linspace((float) 1., (float) 10000., 5);
    for(int i = 0; i < length; i++)
    {
        z_test_values[i] = 1/(z_test_values[i]*z_test_values[i]);
    }
    z_test_values.push_back(0.); // Also test when zero
    length++;
    // 'correct' results from Wolfram Alpha.
    std::vector<float> res_wolfram{3.14159265358979,
                                   0.00079975998530,
                                   0.00039996000066,
                                   0.00026665783397,
                                   0.00020000000033,
                                   0.};
    float res;
    for(int i = 0; i < length; i++)
    {
        res = incompleteBeta(0.5, 0.5, z_test_values[i]);
        BOOST_CHECK(AreSame(res, res_wolfram[i]));
    }

    // Compare varying values for a, for z = 0.01, b = 0.5
    // We need to avoid the singularities, which occur at zero and negative integers.
    length = 6;
    std::vector<float> a_test_values = linspace((float) -2.1, (float) 2.5, length);
    std::vector<float> res_wolfram_2{-7625.08886337438196511846,
                                     -200.49433167596152234321,
                                     -12.71336689572775970392,
                                     0.07266455552386126564,
                                     0.00043921575518787742,
                                     4.0143696200406766891e-6};
    for(int i = 0; i < length; i++)
    {
        res = incompleteBeta(a_test_values[i], 0.5, 0.01);
        BOOST_CHECK(AreSame(res, res_wolfram_2[i]));
    }
}

