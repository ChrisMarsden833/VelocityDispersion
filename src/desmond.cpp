#include "desmond.h"

// Equation 1 - just a sersic profile
float MassDensityProfile(float r, float SersicIndex, float Half_Light_radius, float stellar_mass)
{
    float b_value = b_n(SersicIndex);

    //float x = b_value * pow((r/Half_Light_radius), 1./SersicIndex);

    float gamma = boost::math::tgamma_lower(2*SersicIndex, b_value);

    if (gamma == 0)
    {
        gamma = 0.0000001;
    }

    float Sigma_e = stellar_mass / (Half_Light_radius * Half_Light_radius * PI * 2. * 2. * SersicIndex * exp(b_value) * gamma / pow(b_value, 2*SersicIndex) );


    if(r == 0.)
    {
        float result = Sigma_e * exp(b_value);

        if(std::isinf(result))
        {
            printf("-- Triggered inf condition: %f\n", result);
            printf("-- Sigma_e: %f\n", Sigma_e);
            printf("---- HLR: %f\n", Half_Light_radius);
            printf("---- SM: %f\n", stellar_mass);
            printf("-- b_value: %f\n", b_value);
            printf("-- Sersic Index: %f\n", SersicIndex);
            printf("-- gamma: %f\n", gamma);
        }
        // Save a small amount of computations, as this will mostly be called at zero.
        return result;
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
    float result = Sigma_e * exp(-b_value * (pow(internal_term, power) - 1.));

    return result;

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

float rho(float r, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c)
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

    return rho_0_term * power_term * exp_term;// + NFW_profile(r, dm_rs, dm_c, omega_m, H);
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

float cumSpherMassDistro(float R, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c)
{
    std::vector<float> args = {Half_Light_radius, SersicIndex, stellar_mass, dm_rs, dm_c};

    float accuracy = SimpsonsRule(cumSpherRho, 0.0001, R, Prepass_Subdivisions, args)/1000;

    return AdaptiveRichardsonExtrapolate(cumSpherRho, 0.0001, R, accuracy, args);
}


float K_Kernel_DW(float u, float beta)
{
    if(u == 0.)
    {
        return 0.;
    }
    float prefactor = 0.5 * pow(u, (2.* beta - 1.));
    float term1 = ((3./2.) - beta) * pow(PI, 0.5) * boost::math::tgamma(beta - 0.5) / boost::math::tgamma(beta);
    float term2 = beta * incompleteBeta(beta + 0.5, 0.5, 1./(u*u));
    float term3 = -incompleteBeta(beta - 0.5, 0.5, 1./(u*u));

    if(term1 + term2 + term3 < 0)
    {
        return 0; // This can occur for low values of beta, and the values are typically very close to zero, so this is presumably okay.
    }

    float res = prefactor * (term1 + term2 + term3);

    return res;
}

float full_sigma_integral_internals(float r, float R, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c)
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
    float density = rho(r, Half_Light_radius, SersicIndex, stellar_mass, dm_rs, dm_c);
    float Mass = cumSpherMassDistro(r, Half_Light_radius, SersicIndex, stellar_mass, dm_rs, dm_c);
    float res;



    res = Kernal * density * Mass /   r;

    return res;
}


float sigma_internals_wrapper(float r, std::vector<float> args)
{
    return full_sigma_integral_internals(r, args[0], args[1], args[2], args[3], args[4], args[5], args[6]);
}


float sigma_los(float R, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c)
{
    std::vector<float> args = {R, beta, Half_Light_radius, SersicIndex, stellar_mass, dm_rs, dm_c};

    float upper_limit = 100*R;

    float accuracy = SimpsonsRule(sigma_internals_wrapper, R, upper_limit, Prepass_Subdivisions, args)/1000.;

    float numerator = 2 * GR * AdaptiveRichardsonExtrapolate(sigma_internals_wrapper, R, upper_limit, accuracy, args);
    float denominator = MassDensityProfile(R, SersicIndex, Half_Light_radius, stellar_mass);

    float value = pow(numerator/denominator, 0.5);
    
    return value;
}


float sigma_los_wrapper(float R, std::vector<float> args)
{
    return sigma_los(R, args[0], args[1], args[2], args[3], args[4], args[5]);
}

float sigma_apature_internals(float r, std::vector<float> args)
{
    float sigma = sigma_los(r, args[0], args[1], args[2], args[3], args[4], args[5]);

    float MDP = MassDensityProfile(r, args[2], args[1], args[3]);

    return MDP * sigma * sigma * r;
}


float sigma_aperture(float R_ap, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c)
{
    std::vector<float> args = {beta, Half_Light_radius, SersicIndex, stellar_mass, dm_rs, dm_c};
    std::vector<float> args2 = {SersicIndex, Half_Light_radius, stellar_mass};

    float accuracy = SimpsonsRule(sigma_apature_internals, 0., R_ap, Prepass_Subdivisions, args)/1000;
    float numerator = AdaptiveRichardsonExtrapolate(sigma_apature_internals, 0., R_ap, accuracy, args);

    accuracy = SimpsonsRule(MassDensityProfile_wrapper, 0., R_ap, Prepass_Subdivisions, args2)/1000;
    float denominator = AdaptiveRichardsonExtrapolate(MassDensityProfile_wrapper, 0., R_ap, accuracy, args2);

    // Error Handling
    bool nan_fail = (isnan(numerator) + isnan(denominator));
    bool dem_zero_fail = (denominator == 0);
    bool neg_fail = numerator/denominator < 0.;

    if (nan_fail || dem_zero_fail || neg_fail)
    {
        if(nan_fail) std::cout << "ERROR - Aperture Equation (10) sigma_aperture() will return NaN. Numerator:" << numerator << " denominator: " << denominator << std::endl;
        if(dem_zero_fail) std::cout << "ERROR - Aperture Equation (10) sigma_aperture()'s denominator has evaluated to zero, causing a divide by zero error." << std::endl;
        if(neg_fail) std::cout << "ERROR - Aperture Equation (10) sigma_aperture() has evalulated to < 0 before the square root, and complex values are not physical here." << std::endl;
        std::cout << "-- Argument R_ap (aperture size): " << R_ap << std::endl;
        std::cout << "-- Argument beta: " << beta << std::endl;
        std::cout << "-- Argument Half Light Radius: " << Half_Light_radius << std::endl;
        std::cout << "-- Argument SersicIndex: " << SersicIndex << std::endl;
        std::cout << "-- Argument stellar mass: " << stellar_mass << std::endl;
        std::cout << "--- Value will be returned as zero, but you should think about why this occured." << std::endl;
        return 0.;
    }

    return pow(numerator/denominator, 0.5);
}
