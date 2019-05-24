#include "testing_utillity.h"

void test_pochhammer(void)
{
   // Test n = 0 results
   test_pochhammer_internal(5., 0, 1.);
   test_pochhammer_internal(0., 0, 1.);
   test_pochhammer_internal(-5, 0, 1.);
   // n = 1 results
   test_pochhammer_internal(5., 1, 5.);
   test_pochhammer_internal(0., 1, 0.);
   test_pochhammer_internal(-5., 1, -5.);
   // n = 2 results
   test_pochhammer_internal(5., 2, 30.);
   test_pochhammer_internal(0., 2, 0.);
   test_pochhammer_internal(-5., 2, 20.);
   // n = 3 results
   test_pochhammer_internal(5., 3, 5.*(5. + 1.)*(5. + 2.));
}

void test_pochhammer_internal(float q, int n, float value)
{
    try
    {
        float res = pochhammer(q, n);

        if(res != value)
        {
            char buffer [50];
            sprintf(buffer, "pochammer (%.1f)_%i returned %.1f, should be %f.", q, n, res, value);
            throw std::invalid_argument(buffer);
        }
    }
    catch (const std::invalid_argument& ia) {
        std::cerr << "Test failed: " << ia.what() << '\n';
    }
}
