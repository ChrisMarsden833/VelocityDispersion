#include "main.h"


int main()
{

    std::cout << "Starting" << std::endl;

    float SM = 11.086;
    float size = 3.2;
    float n = 3.12;
    float beta = 0.115;

    float res = GetUnweightedVelocityDispersion(4.*size, beta, size, n, SM, 0.0);

    std::cout << res << std::endl;

    
    return 0;


}
