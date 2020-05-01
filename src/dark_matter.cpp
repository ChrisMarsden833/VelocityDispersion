#include "dark_matter.h"

float delta_vir(void)
{
    float x = omegam - 1;
    return (18*(PI*PI) + 82*x - 39*x*x)/(x+1);
}

float critical_density(void)
{
    float Hu = H/1000; // Convert H into units of (km s^-1 kPc^-1)
    return (3*Hu*Hu)/(8*PI*GR);
}

float delta_char(float c)
{
    return (delta_vir() * omegam * (c*c*c))/(3 * (log(1+c) + c/(1+c)));
}

float NFW_profile(float r, float rs, float c)
{
    return critical_density() * delta_char(c)/((r/rs)*(1+r/rs)*(1+r/rs));
}


