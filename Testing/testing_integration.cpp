#include "testing_integration.h"

bool test_ARE_sin(void)
{
    std::vector<float> nothing;
    float accuracy = 0.001;

    float res = AdaptiveRichardsonExtrapolate(testing_sin, 0, 2*PI, accuracy, nothing);

    if((int) res != 0)
    {
        std::cout << "Int_cast Failed" << std::endl;
        return false;
    }
    else if(res > accuracy)
    {
        std::cout << "Accuracy test failed" << std::endl;
        return false;
    }

    std::cout << "Integral of sin between zero and 2 PI returns as:" << res << std::endl;

    return false;
}

bool test_RE_exception(void)
{
    // This function should run, but should also throw an exception, complaining about odd numbers
    std::vector<float> nothing;
    RResult Values;
    Values = RichardsonExtrapolate(testing_sin, 0, 2*PI, 5, nothing);

    return true;
}

float testing_sin(float x, std::vector<float>)
{
    return sin(x);
}

