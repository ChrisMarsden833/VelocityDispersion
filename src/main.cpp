#include "main.h"


int main()
{

    std::cout << "Starting" << std::endl;

    float SM = 11.1;
    float size = 3.2;
    float n = 3.12;
    float beta = 0.000015;

    float res = GetUnweightedVelocityDispersion(1e-2, beta, size, n, SM, 0.0);

    std::cout << res << std::endl;

    
    return 0;


}
