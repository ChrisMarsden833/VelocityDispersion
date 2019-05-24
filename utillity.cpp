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
    if(a == -1)
    {
        try
        {
                char buffer [50];
                sprintf(buffer, "incompleteBeta called with Value of a equal to exactly -1 (bad news)");
                throw std::invalid_argument(buffer);
        }
        catch (const std::invalid_argument& ia) {
            std::cerr << "Singularity: " << ia.what() << " Value will be perturbed by 0.0001" << '\n';
        }
        //a += 0.0001;
    }
    float value = (pow(z, a)/a) * hyperGeometricSeries(a, 1.-b, a+1., z);

    return value;
}

double pochhammer(double q, int n)
{
    if(n == 0)
    {
        return 1.;
    }
    else
    {
        double value = 1.;
        for(int i = 0; i <= n - 1; i++)
        {
            value *= q + (double) i;
        }
        return value;
    }
}

float hyperGeometricSeries(float a, float b, float c, float z)
{
    if( ((int) a == a && a < 0))
    {
        // An entirely different process happens in these cases.
        double value = 0;

        for(int n = 0; n <= a; n++)
        {
            value += pow(-1., n) * (factorial(-a)/(factorial(n)*factorial(-a-n))) * (pochhammer(b, n)/pochhammer(c, n)) * pow(z, n);
        }
        return value;

    }
    else if ((int) b == b && b < 0)
    {
        double value = 0;

        for(int n = 0; n <= a; n++)
        {
            value += pow(-1., n) * (factorial(-b)/(factorial(n)*factorial(-b-n))) * (pochhammer(a, n)/pochhammer(c, n)) * pow(z, n);
        }
        return value;
    }
    
    double accuracy = 0.0000001;

    int number_of_iterations = 100000000;

    double value = 0.;
    double oldvalue;

    double term1,
          term2,
          term3,
          problem_term,
          sign,
          factorial_term,
          power_term;


    for(int n = 0; n < number_of_iterations; n++)
    {
        oldvalue = value;
        term1 = pochhammer((double) a, n);
        term2 = pochhammer((double) b, n);
        term3 = pochhammer((double) c, n);
        factorial_term = boost::math::factorial<double>(n);
        power_term = pow(z, n);

        //value += (term1 * term2 * power_term) / (term3 * factorial_term);

        sign = getSign(term1) * getSign(term2) * getSign(term3);

        problem_term = log10(abs(term1)) + log10(abs(term2)) + log10(power_term) - log10(abs(term3)) - log10(factorial_term);

        value += sign * pow(10, problem_term);


        if (abs(oldvalue - value) < accuracy)
        {
            number_of_iterations = 0;
        }
    }
    return (float) value;
}

float getSign(float argument)
{
    if(argument < 0)
    {
        return -1.;
    }
    else
    {
        return 1.;
    }

}




