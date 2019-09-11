#include "utillity.h"

template<typename T>
T factorial(T n)
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
    float beta;
    if(a <= 0)
    {
        beta = pow(z, a) * pow((1 - z), b) + (a + b) * incompleteBeta(a + 1, b, z);
        beta = beta/a;
        return beta;
    }
    else
    {
        beta = boost::math::beta(a, b, z);
        return beta;
    }
}


void CheckAndMaybeIncrementN(int * N, std::string fname)
{
    try
    {
        if(!checkEven(*N))
        {
            std::string message = std::to_string(*N);
            throw std::invalid_argument("Odd Number supplied for " + fname + " (" + message + ")");
        }
    }
    catch (const std::invalid_argument& ia) {
        std::cerr << ia.what() << "\n" << "Value will be incremented by one." << '\n';
        ++*N;
    }
}

void fastLinspace(float * grid, float * h, float a, float b, int N = 100)
{
    grid = (float *) malloc( (N + 1) * sizeof(float)); // Prep array
    *h = (b - a)/((float) N); // Calculate the spacing of the grid, h
    for(int i = 0; i <= N; i++) grid[i] = a + *h * (float) i; // fill values
}




