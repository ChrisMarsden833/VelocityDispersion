#include <iostream>
#include <math.h>

#include "../src/desmond.h"

float wrapper(float R, float HLR, float SersicIndex)
{
    std::vector<float> args = {HLR, SersicIndex};

    float a = 1.;
    float b = 2.;
    int N = 100;



    return 0.5;
}


extern "C"
{
    float p_nex(float SersicIndex)
    {
        float res = p_n(SersicIndex);

        return res;
    }

    // float innards(float r, float R, float HLR, float beta, float SersicIndex)
    float innards(float r, float R, float beta, float HLR, float SersicIndex, float stellar_mass, float dm_rho0, float dm_rs)
    {
        float res = K_Kernel_DW(r/R, beta);
        res = rho(r, HLR, SersicIndex, stellar_mass, dm_rho0, dm_rs);
        //* cumSpherMassDistro(R, HLR, SersicIndex) / r;
        return res;
    }
}