#include "dark_matter.h"

float delta_vir(float omega_m)
{
    float x = omega_m - 1;
    return (18*(PI*PI) + 82*x - 39*x*x)/(x+1);
}

float critical_density(float H)
{
    H /= 1000; // Convert H into units of (km s^-1 kPc^-1)
    return (3*H*H)/(8*PI*GR);
}

float delta_char(float omega_m, float c)
{
    return (delta_vir(omega_m) * omega_m * (c*c*c))/(3 * (log(1+c) + c/(1+c)));
}

float NFW_profile(float r, float rs, float c, float omega_m, float H)
{
    return critical_density(H) * delta_char(omega_m, c)/((r/rs)*(1+r/rs)*(1+r/rs));
}


