#include "testing_integration.h"

bool test_ARE_sin(void)
{
    std::vector<float> nothing;
    float res = AdaptiveRichardsonExtrapolate(testing_sin, 0, 2*PI, 0.001, nothing);
    std::cout << res << std::endl;
    return false;
}

bool test_RE_exception(void)
{
    std::vector<float> nothing;
    RResult Values;
    Values = RichardsonExtrapolate(testing_sin, 0, 2*PI, 5, nothing);

    return true;
}

float testing_sin(float x, std::vector<float>)
{
    return sin(x);
}

