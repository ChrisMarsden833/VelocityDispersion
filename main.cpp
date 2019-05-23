#include "main.h"


int main()
{


    float a, b;

    a = boost::math::ibeta(0.1, 0.2, 0.2);

    std::cout << "Boost: " << a << std::endl;

    b = incompleteBeta(0.1, 0.2, 0.2);
    std::cout << "Identity: " << b << std::endl;

    return 0;
}