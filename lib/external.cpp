#include <iostream>
#include <math.h>
#include "omp.h"
#include "../src/desmond.h"

extern "C"
{
    float p_nex(float SersicIndex)
    {
        float res = p_n(SersicIndex);

        return res;
    }

    // float innards(float r, float R, float HLR, float beta, float SersicIndex)
    float density(float r, float HLR, float SersicIndex, float stellar_mass, float dm_rs, float dm_c, float omega_m, float H)
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

    float sersic_profile(float r, float SersicIndex, float Half_Light_radius, float stellar_mass)
    {
        float value = MassDensityProfile(r, SersicIndex, Half_Light_radius, stellar_mass);
        return value; //log10(value);
    }

    float incbeta(float a, float b, float z)
    {
        return incompleteBeta(a, b, z);
    }

    float rho0(float Half_Light_radius, float SersicIndex, float stellar_mass)
    {
        return rho_0(Half_Light_radius, SersicIndex, stellar_mass);
    }

    float VD_los(float R, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c, float omega_m, float H)
    {
        return sigma_los(R, beta, Half_Light_radius, SersicIndex, stellar_mass, dm_rs, dm_c, omega_m, H);
    }

    float VD(float R_ap, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c, float omega_m, float H)
    {
        return sigma_aperture(R_ap, beta, Half_Light_radius, SersicIndex, stellar_mass, dm_rs, dm_c, omega_m, H);
    }

    void ArrayTest(float * input, int size)
    {
        for(int i = 0; i<size; i++)
        {
            printf("Element %i is %f\n", i, input[i]);
        }

        #pragma omp parallel
        {
            int id = omp_get_thread_num();
            printf("Hello from thread %i\n", id);
        };

    }

    float * ParallelSigma(float * Aperture,
                       float * Beta,
                       float * HalfLightRadius,
                       float * SersicIndex,
                       float * StellarMass,
                       float * HaloRs,
                       float * HaloC,
                       float OmegaM,
                       float H,
                       int size)
    {
        float * res = (float *) std::malloc(sizeof(float) * size);

        #pragma omp parallel for default(none) shared(size, res, Aperture, Beta, HalfLightRadius, SersicIndex, StellarMass, HaloRs, HaloC, OmegaM, H)
        for (int i = 0; i < size; i++)
        {
            res[i] = sigma_aperture(Aperture[i], Beta[i], HalfLightRadius[i], SersicIndex[i], StellarMass[i], HaloRs[i], HaloC[i], OmegaM, H);
        }
        return res;
    }
}