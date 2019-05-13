#include "desmond.h"

float MassDensityProfile(float r, float SersicIndex, float Half_Light_radius)
{
    float Sigma_e = 1.;
    float b_value = b_n(SersicIndex);
    float power = 0.;
    float internal_term = 0.;

    if(SersicIndex != 0.)
    {
        float power = 1./SersicIndex;
    }
    if(Half_Light_radius != 0.)
    {
        internal_term = r/Half_Light_radius;
    }
    return Sigma_e * exp(-b_value * (pow(internal_term, power) - 1.));
}

float b_n(float SersicIndex)
{
    if(SersicIndex == 0.)
    {
        return -1./3.; // Zero protection, but might not be reasonable.
    }
    return 2. * SersicIndex - (1./3.) + (.009876/SersicIndex);
}

float rho(float r, float Half_Light_radius, float SersicIndex)
{
    // Rho Zero
    float rho_0_term = rho_0(Half_Light_radius, SersicIndex);

    // First 'power' term
    float p_n_term = -p_n(SersicIndex);
    float radius_ratio = 0.;
    if(Half_Light_radius != 0.)
    {
        radius_ratio = r/Half_Light_radius;
    }

    float power_term = 0;
    if(radius_ratio != 0 || p_n_term > 0)
    {
        power_term = pow(radius_ratio, p_n_term);
    }

    // Exponential Term
    float b_n_term = -b_n(SersicIndex);

    float inverse_n = 0.;
    if(SersicIndex != 0.)
    {
        inverse_n = 1./SersicIndex;
    }

    float exp_power_term = 0;
    if(radius_ratio != 0 || inverse_n > 0)
    {
        exp_power_term = pow(radius_ratio, inverse_n);
    }

    float exp_term = exp(b_n_term * exp_power_term);

    return rho_0_term * power_term * exp_term;
}

float rho_0(float Half_Light_radius, float SersicIndex)
{
    float Sigma_0 = MassDensityProfile(0., SersicIndex, Half_Light_radius);
    float left = Sigma_0 * pow(b_n(SersicIndex), SersicIndex * (1. - p_n(SersicIndex)))/(2.*Half_Light_radius);
    float right = (boost::math::tgamma(2*SersicIndex)/boost::math::tgamma(SersicIndex*(3. - p_n(SersicIndex))));
    return left * right;
}

float p_n(float SersicIndex)
{
    if(SersicIndex == 0.)
    {
        return 1.;
    }
    return 1. - .6097/SersicIndex + .00563/(SersicIndex*SersicIndex);
}

float cumSpherRho(float R, std::vector<float> args)
{
    return 4.*PI*R*R*rho(R, args[0], args[1]);
}

float cumSpherMassDistro(float R, float Half_Light_radius, float SersicIndex)
{
    std::vector<float> args = {Half_Light_radius, SersicIndex};
    return AdaptiveRichardsonExtrapolate(cumSpherRho, 0.0, R, .0001, args);
}