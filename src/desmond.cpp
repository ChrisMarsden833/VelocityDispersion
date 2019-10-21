#include "desmond.h"

// Equation 1 - just a sersic profile
float MassDensityProfile(float r, float SersicIndex, float Half_Light_radius, float stellar_mass)
{
    float Sigma_e = stellar_mass/(PI * Half_Light_radius * Half_Light_radius);
    float b_value = b_n(SersicIndex);

    if(r == 0.)
    {
        // Save a small amount of computations, as this will mostly be called at zero.
        return Sigma_e * exp(b_value);
    }

    float power = 0.;
    float internal_term = 0.;

    if(SersicIndex != 0.)
    {
        power = 1./SersicIndex;
    }
    if(Half_Light_radius != 0.)
    {
        internal_term = r/Half_Light_radius;
    }
    return Sigma_e * exp(-b_value * (pow(internal_term, power) - 1.));
}

float MassDensityProfile_wrapper(float r, std::vector<float> args)
{
    return MassDensityProfile(r, args[0], args[1], args[2]) * r;
}


float b_n(float SersicIndex)
{
    if(SersicIndex == 0.)
    {
        return -1./3.; // Zero protection, but might not be reasonable.
    }
    return 2. * SersicIndex - (1./3.) + (.009876/SersicIndex);
}

float rho(float r, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rho0, float dm_rs)
{
    // Rho Zero
    float rho_0_term = rho_0(Half_Light_radius, SersicIndex, stellar_mass);

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

    return rho_0_term * power_term * exp_term + NFW_profile(r, dm_rho0, dm_rs);
}

float rho_0(float Half_Light_radius, float SersicIndex, float stellar_mass)
{
    float Sigma_0 = MassDensityProfile(0., SersicIndex, Half_Light_radius, stellar_mass);
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
    return 4.*PI*R*R*rho(R, args[0], args[1], args[2], args[3], args[4]);
}

float cumSpherMassDistro(float R, float Half_Light_radius, float SersicIndex, float stellar_mass)
{
    std::vector<float> args = {Half_Light_radius, SersicIndex, stellar_mass};

    float accuracy = SimpsonsRule(cumSpherRho, 0.0, R, Prepass_Subdivisions, args)/10000;

    return AdaptiveRichardsonExtrapolate(cumSpherRho, 0.0, R, accuracy, args);
}

float K_Kernel_DW(float u, float beta)
{
    if(u == 0.)
    {
        return 0.;
    }
    float prefactor = 0.5 * pow(u, (2* beta - 1));
    float term1 = ((3./2.) - beta) * pow(PI, 0.5) * boost::math::tgamma(beta - 0.5) / boost::math::tgamma(beta);
    float term2 = beta * incompleteBeta(beta + 0.5, 0.5, 1./(u*u));
    float term3 = -incompleteBeta(beta - 0.5, 0.5, 1./(u*u));
    return prefactor * (term1 + term2 + term3);
}

float full_sigma_integral_internals(float r, float R, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rho0, float dm_rs)
{
    if(r == 0.)
    {
        return 0.;
    }

    float ratio;

    if(R == 0.)
    {
        ratio = 0.;
    }
    else
    {
        ratio = r/R;
    }

    float Kernal = K_Kernel_DW(ratio, beta);
    float density = rho(r, Half_Light_radius, SersicIndex, stellar_mass, dm_rho0, dm_rs);
    float Mass = cumSpherMassDistro(R, Half_Light_radius, SersicIndex, stellar_mass);

    return Kernal * density * Mass / r;
}

float sigma_internals_wrapper(float r, std::vector<float> args)
{
    return full_sigma_integral_internals(r, args[0], args[1], args[2], args[3], args[4], args[5], args[6]);
}

float sigma_los(float R, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rho0, float dm_rs)
{
    std::vector<float> args = {R, beta, Half_Light_radius, SersicIndex, stellar_mass, dm_rho0, dm_rs};

    float upper_limit = 1000*R;

    float accuracy = SimpsonsRule(sigma_internals_wrapper, R, upper_limit, Prepass_Subdivisions, args)/10000;

    float numerator = 2 * GR * AdaptiveRichardsonExtrapolate(sigma_internals_wrapper, R, upper_limit, accuracy, args);
    float denominator = MassDensityProfile(R, SersicIndex, Half_Light_radius, stellar_mass);
    if(denominator == 0. || numerator == 0.)
    {
        return 0.;
    }
    return pow(numerator/denominator, 0.5);
}

float sigma_los_wrapper(float R, std::vector<float> args)
{
    return sigma_los(R, args[0], args[1], args[2], args[3], args[4], args[5]);
}

float sigma_apature_internals(float r, std::vector<float> args)
{
    float sigma = sigma_los(r, args[0], args[1], args[2], args[3], args[4], args[5]);
    return MassDensityProfile(r, args[1], args[2], args[3]) * sigma * sigma * r;
}

float sigma_aperture(float R_ap, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass)
{
    std::vector<float> args = {beta, Half_Light_radius, SersicIndex, stellar_mass};
    std::vector<float> args2 = {Half_Light_radius, SersicIndex, stellar_mass};

    float accuracy = SimpsonsRule(sigma_apature_internals, 0., R_ap, Prepass_Subdivisions, args)/10000;
    float numerator = AdaptiveRichardsonExtrapolate(sigma_apature_internals, 0., R_ap, accuracy, args);

    accuracy = SimpsonsRule(MassDensityProfile_wrapper, 0, R_ap, Prepass_Subdivisions, args2);
    float denominator = AdaptiveRichardsonExtrapolate(MassDensityProfile_wrapper, 0, R_ap, accuracy, args);

    if(denominator == 0. || numerator == 0.)
    {
        return 0.;
    }

    return pow(numerator/denominator, 0.5);
}

float NFW_profile(float r, float rho0, float Rs)
{
    return rho0/((r/Rs)*(1+(r/Rs))*(1+(r/Rs)));
}




