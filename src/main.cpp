#include "main.h"


int main()
{

    Galaxy aGalaxy(12., 0.4, 10., 10., 4.);

    float test = aGalaxy.cumulative_mass(10.);

    std::cout << "Result: " << test << std::endl;


    float old_test = cumSpherMassDistro(10., 10., 4., pow(10., 12.), 10., 10.);
    std::cout << "Result2: " << old_test << std::endl;

    return 0;
}
