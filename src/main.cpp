#include "main.h"


int main()
{

   float r = 10;
   float HLR = 30;
   float beta = 0.1;
   float sersic_index = 2.;
   float stellar_mass = 10e7;
   float halo_mass = 10e8;
   float halo_size = 2. * HLR;

   float res = rho(r, HLR, sersic_index, stellar_mass, halo_mass, halo_size);

   std::cout << res;

    return 0;
}