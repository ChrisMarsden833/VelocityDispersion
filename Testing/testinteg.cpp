#import "testinteg.h"

void test_integ(void)
{
    test_simpson();
}

void test_simpson(void)
{
    float a = 0., b = 2*PI;

    float In;
    In = SimpsonsRule(cos, 0., 0.5*PI, 3);

    std::cout << "integral" << In << std::endl;
}

