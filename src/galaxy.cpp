#include "galaxy.h"

Galaxy::Galaxy(float input_stellar_mass,
	       float input_beta,
	       float input_half_light_radius,
	       float input_aperture_size,
	       float input_sersic_index)
{
	stellar_mass = input_stellar_mass;
	beta = input_beta;
	half_light_radius = input_half_light_radius;
	sersic_index = input_sersic_index;
	
	dark_matter_on = false;

	// Effective constants for divide by zero errors
	sersic_index_eff = sersic_index;
	if(sersic_index == 0.) sersic_index_eff = zero_perturbation;
	half_light_radius_eff = half_light_radius;
	if(half_light_radius == 0.) half_light_radius_eff = zero_perturbation;
	
	// Set the value of b_n 	
	b_n = 2.*sersic_index_eff - (1./3.) + (.009876/sersic_index_eff); 
	
	// find the gamma function used immedately in sigma_e
	float gamma = boost::math::tgamma_lower(2*sersic_index, b_n);
	// Sigma_e, the constant for the sersic profile.
	sigma_e = pow(10, stellar_mass)/( pow(half_light_radius, 2) * PI * 2.* 2.*sersic_index*exp(b_n)*gamma/pow(b_n, 2.*sersic_index));
	// p_n 
	p_n = 1. - .6097/sersic_index_eff  + .00563/(sersic_index_eff*sersic_index_eff);

	// Rho0
	float Sigma0 = this->MassDensity(0.);

	float left = Sigma0 *  pow(b_n, (sersic_index * (1. - p_n))) / (2.*half_light_radius_eff);
	float right = (boost::math::tgamma(2*sersic_index)/boost::math::tgamma(sersic_index*(3. - p_n)));
	rho0 = left * right;

}

void Galaxy::init_dark_matter(string input_profile_name,
			      float input_concentration)
{
	dark_matter_on = true;
	profile_name = input_profile_name;
	concentration = input_concentration;
}

float Galaxy::MassDensity(float r)
{
	float power = 1./sersic_index_eff;
	float internal_term = r/half_light_radius_eff;
	
	float result = sigma_e * exp(-b_n * (pow(internal_term, power) - 1.));	      
	return result;
}

float Galaxy::rho(float r)
{

    float radius_ratio = r/half_light_radius_eff;

    float radius_ratio_eff = radius_ratio;
    //if(radius_ratio == 0. && -p_n <= 0) radius_ratio_eff = zero_perturbation;
    float power_term = pow(radius_ratio_eff, -p_n);

    float inverse_n = 1./sersic_index_eff;
    float exp_power_term = pow(radius_ratio_eff, inverse_n);

    float exp_term = exp(-b_n * exp_power_term);
    
    float dark_matter_term = 0.0;

    if(dark_matter_on)
    {
	    // TODO Call dark matter function
    }
   
    return rho0 * power_term * exp_term + dark_matter_term; // + NFW_profile(r, dm_rs, dm_c, omega_m, H);
}

float Galaxy::mass_shell(float R)
{
	float res = 4.0*PI *R*R* this->rho(R);
	return res;
}


float Galaxy::cumulative_mass(float R)
{
	float accuracy = pow(10., stellar_mass-4);

	auto fp = bind(&Galaxy::mass_shell, this, _1);

	return AdaptiveRichardsonExtrapolate(fp, zero_perturbation, R, accuracy);
}

float Galaxy::K_Kernel_DW(float u)
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

float Galaxy::sigma_integrand(float r)
{ 
    if(R == 0) R = zero_perturbation;

    float ratio = r/R;

    float Kernal = this->K_Kernel_DW(ratio);
    float density = this->rho(r);
    float Mass = this->cumulative_mass(R);
	
    float res;

    if(r == 0) r = zero_perturbation;

    res = Kernal * density * Mass /  r;

    return res;


}
