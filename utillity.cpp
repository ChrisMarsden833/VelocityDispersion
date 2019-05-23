#include "utillity.h"

int factorial(int n)
{
    if (n == 0)
        return 1; // base case
    else
        return n * factorial(n-1); // recursive case
}

bool checkEven(int Value)
{
    // Test N to make sure it is even.
    if(Value & 1)
    {
        return false;
    }
    else {
        return true;
    }

}

float incompleteBeta(float a, float b, float z)
{
    float value = pow(-1., a) * boost::math::ibeta(a, 1 - a - b, z/(z-1));
    return value;
}

