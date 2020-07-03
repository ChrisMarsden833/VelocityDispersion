#include "main.h"


int main()
{

    std::cout << "Starting" << std::endl;

    float ap_size(10.0);
    float beta(0.15);
    float hlr = ap_size;
    float n = 4.;
    float sm = 10.;
    float z = 0.0;
    float hm = 13.;
    char * name = (char *)"NFW";
    

    float sigma = GetVelocityDispersion(ap_size, beta, hlr, n, sm, z, hm, name);

    std::cout << "Sigma: " << sigma << std::endl;

    
    return 0;


}
