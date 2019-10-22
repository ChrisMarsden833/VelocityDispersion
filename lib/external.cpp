#include <iostream>
#include <math.h>
#include "../src/desmond.h"

extern "C"
{
    float p_nex(float SersicIndex)
    {
        float res = p_n(SersicIndex);

        return res;
    }

    // float innards(float r, float R, float HLR, float beta, float SersicIndex)
    float density(float r, float R, float beta, float HLR, float SersicIndex, float stellar_mass, float dm_rs, float dm_c, float omega_m, float H)
    {
        float res;
        res = rho(r, HLR, SersicIndex, stellar_mass, dm_rs, dm_c, omega_m, H);
        //* cumSpherMassDistro(R, HLR, SersicIndex) / r;
        return res;
    }

    float total_mass(float R, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c, float omega_m,  float H)
    {
        float integral = cumSpherMassDistro(R, Half_Light_radius, SersicIndex, stellar_mass, dm_rs, dm_c, omega_m, H);
        return integral;
    }

    float K_Kernel(float u, float beta)
    {
        return K_Kernel_DW(u, beta);
    }

    float first_integral_internals(float r, float R, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c, float omega_m, float H)
    {
        return full_sigma_integral_internals(r, R, beta, Half_Light_radius, SersicIndex, stellar_mass, dm_rs, dm_c, omega_m, H);
    }

}