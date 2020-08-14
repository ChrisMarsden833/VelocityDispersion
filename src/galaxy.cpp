#include "galaxy.h"

Galaxy::Galaxy(float input_stellar_mass,
	       float input_beta,
	       float input_half_light_radius,
	       float input_aperture_size,
	       float input_sersic_index,
	       float z)
{
    /* Constructor. Here we do all the things we only need to do once */

    // First manage the variables into their respective member variables

    stellar_mass = input_stellar_mass;
    beta = input_beta;
	half_light_radius = input_half_light_radius;
	sersic_index = input_sersic_index;
	aperture_size = input_aperture_size;
	redshift = z;

	// Manage beta, in the cases of known sigularities (TODO - generalize this, as  I *think* they occur every 0.5n
    if(beta == 0.0) beta = zero_perturbation;
    if(beta == 0.5) beta += zero_perturbation;

    // Calculate the gamma term within the Kernel, as we use this a lot and it only depends on beta
    gamma_term = (float) boost::math::tgamma(beta - 0.5) / boost::math::tgamma(beta);

    // Check stellar mass is not in units of M_sun. This occurs more often than you might think...
    assert(stellar_mass < 20 &&  "Stellar mass is unexpectedly high (>20) - are the units log10 [M_sun]?");

	// Effective constants for divide by zero errors.
	sersic_index_eff = sersic_index;
	if(sersic_index == 0.) sersic_index_eff = zero_perturbation;
	half_light_radius_eff = half_light_radius;
	if(half_light_radius == 0.) half_light_radius_eff = zero_perturbation;
	
	// Set the value of b_n (equation 2)
	b_n = 2.*sersic_index_eff - (1./3.) + (.009876/sersic_index_eff); 

	// find the gamma function used immedately in sigma_e. This is used by my defn of sigma_e.
	// not sure this is required?
	float gamma = boost::math::tgamma_lower(2*sersic_index, b_n);

	// Sigma_e, the constant for the sersic profile. It's debatable if this should just be M_*/(2*pi*Re^2) or this:
	sigma_e = pow(10, stellar_mass)/( pow(half_light_radius_eff, 2) * PI * 2 * 2. * sersic_index * exp(b_n) * gamma / pow(b_n, 2*sersic_index) );

	// Asserts
    assert(!isnan(sigma_e) || !isinf(sigma_e) || assert_msg("Sigma_e was calculated as: " << sigma_e << std::endl <<
        "====Location: AGalaxy Constructor" << std::endl <<
        "----Numerator (stellar mass, m_sun) = " << pow(10, stellar_mass) << std::endl <<
        "----Half Light Radius (Kpc) = " << half_light_radius_eff << std::endl <<
        "----Sersic Index = " << sersic_index << std::endl <<
        "----b_n = " << b_n << std::endl));

	// p_n - eqn (5)
	p_n = 1. - .6097/sersic_index_eff  + .00563/(sersic_index_eff*sersic_index_eff);

	// Rho0 - eqn (4)
	Sigma0 = this->MassDensity(0.);
	float left = Sigma0 *  pow(b_n, (sersic_index * (1. - p_n))) / (2.*half_light_radius_eff);
	float right = (boost::math::tgamma(2*sersic_index)/boost::math::tgamma(sersic_index*(3. - p_n)));
	rho0 = left * right;

	float as = half_light_radius/pow(p_n, sersic_index);
	float l = rho0 / pow(b_n, (1. - p_n));

	//mass_prefactor = 4.0 * PI * sersic_index * l * pow(as, 3.);
	mass_prefactor = 4.0 * PI * sersic_index * boost::math::tgamma((3.0-p_n)*sersic_index) * l * pow(as, 3.);
    //mass_prefactor =  2. * PI * sersic_index * boost::math::tgamma(2.*sersic_index) * Sigma0 * pow(as, 2.);

    // Asserts
    assert(!isnan(rho0) || !isinf(rho0) || assert_msg("rho0 was calculated as: " << rho0 << std::endl <<
      "====Location: AGalaxy Constructor" << std::endl <<
      "----Sigma0 (MassDensity(0.)) = " << Sigma0 << std::endl <<
      "----LHS = " << left << std::endl <<
      "----RHS = " << right << std::endl));
}

float Galaxy::MassDensity(float r)
{
    // Equation (1)
	float power = 1./sersic_index_eff;
	float internal_term = r/half_light_radius_eff;
	
	float result = sigma_e * exp(-b_n * (pow(internal_term, power) - 1.));

    // Asserts
    assert(!isnan(result) || !isinf(result) || assert_msg("MassDensity (sersic) about to return " << rho0 << std::endl <<
         "----Power (1/n) = " << power << std::endl <<
         "----internal term (r/R_eff) = " << internal_term << std::endl <<
         "----Sigma_e = " << sigma_e << std::endl <<
         "----RHS = " << right << std::endl));

	return result;
}

float Galaxy::MassDensityr(float r)
{
    // Just a wrapper function, sometimes we need it *r for integrals.
    return this->MassDensity(r) * r;
}

float Galaxy::rho4mass(float r)
{
    // Overall Density, sum of all components

    float stars_term = 0.0;
    float dark_matter_term = 0.0;

    if(stars_on)
    {
        // De-projected density, equation (3)
        float radius_ratio = r/half_light_radius_eff;

        float radius_ratio_eff = radius_ratio;
        if(radius_ratio == 0. && -p_n <= 0) radius_ratio_eff = zero_perturbation;
        float power_term = pow(radius_ratio_eff, -p_n);

        float inverse_n = 1./sersic_index_eff;
        float exp_power_term = pow(radius_ratio_eff, inverse_n);

        float exp_term = exp(-b_n * exp_power_term);

        stars_term = rho0 * power_term * exp_term;

        assert(!isnan(stars_term) || !isinf(stars_term) || assert_msg("Stellar Density (eq3) about to return " << stars_term << std::endl <<
           "----rho0 = " << rho0 << std::endl <<
           "----Power term (r/R_eff)^-p_n = " << power_term << std::endl <<
           "----Exp power term (r/R_eff)^1/n = " << exp_power_term << std::endl <<
           "----Exp term exp(-b_n (Exp_power_term)) = " << exp_power_term << std::endl <<
           "-----r = " << r << std::endl));
    }

    if(dark_matter_on)
    {
        dark_matter_term = HaloDensity(r);
        assert(!isnan(dark_matter_term) || !isinf(dark_matter_term) || assert_msg("Dark Matter Density about to contribute " << stars_term << std::endl <<
                                                                                                                             "-----r = " << r << std::endl));
    }

    float res = stars_term + dark_matter_term;
    return res;
}

float Galaxy::rho(float r)
{
    // Overall Density, sum of all components

    float stars_term = 0.0;

    // De-projected density, equation (3)
    float radius_ratio = r/half_light_radius_eff;

    float radius_ratio_eff = radius_ratio;
    if(radius_ratio == 0. && -p_n <= 0) radius_ratio_eff = zero_perturbation;
    float power_term = pow(radius_ratio_eff, -p_n);

    float inverse_n = 1./sersic_index_eff;
    float exp_power_term = pow(radius_ratio_eff, inverse_n);

    float exp_term = exp(-b_n * exp_power_term);

    stars_term = rho0 * power_term * exp_term;

    assert(!isnan(stars_term) || !isinf(stars_term) || assert_msg("Stellar Density (eq3) about to return " << stars_term << std::endl <<
             "----rho0 = " << rho0 << std::endl <<
             "----Power term (r/R_eff)^-p_n = " << power_term << std::endl <<
             "----Exp power term (r/R_eff)^1/n = " << exp_power_term << std::endl <<
             "----Exp term exp(-b_n (Exp_power_term)) = " << exp_power_term << std::endl <<
             "-----r = " << r << std::endl));

    float res = stars_term;
    return res;
}

float Galaxy::mass_shell(float R)
{
    // Mass shell, based on radius or density
	float res = 4.0*PI *R*R * this->rho4mass(R);
    assert(!isnan(res) || !isinf(res) || assert_msg("Mass Shell about to contribute " << res << std::endl <<
      "-----R = " << R << std::endl <<
      "-----Density (function call to rho(R)) = " << this->rho4mass(R) << std::endl));
	return res;
}

float Galaxy::cumulative_mass(float R_arg)
{
    // Calculate the cumulative mass

    // Bind the Mass shell function to a pointer so we can pass it to the integration functions.
    // An awkward consequence of having the class setup
	auto fp = bind(&Galaxy::mass_shell, this, _1);

	//Preliminary pass to get the rough value for the integral. We need it to do Adaptive Richardson Extrapolation properly.
    float mass_accuracy = SimpsonsRule(fp, 0., R_arg, initial_subdiv)/pow(10., precision + cum_mass_precision_modifier);

    //float res = SimpsonsRule(fp, 0., R_arg, 1000);

    // The actual integration
    float res = AdaptiveRichardsonExtrapolate(fp, 0., R_arg, mass_accuracy);

    assert(!isnan(res) || !isinf(res) || assert_msg("Cumulative mass about to return " << res << std::endl <<
      "-----R_arg = " << R_arg << std::endl));
	return res;
}

float Galaxy::analytic_mass(float r)
{
    float mass = pow(10., 8.);

    if(stars_on)
    {
        float as = half_light_radius_eff/pow(b_n, sersic_index);
        float x = r/as;
        float threemp_term = (3.-p_n) * sersic_index;

        mass += pow(10., stellar_mass) * boost::math::tgamma_lower(threemp_term, pow(x, 1./sersic_index) )/ boost::math::tgamma(threemp_term, 0.0);
    }
    if(dark_matter_on)
    {
        if(profile_name == "NFW")
        {
            float rs = HaloRadius/concentration;

            float fc = log(1.+concentration)-concentration/(1.+concentration);
            float logrhos = HaloMass-log10(4.*PI*rs*rs*rs*fc);

            float y = r/rs;

            float numerator = log(y+1.) - y/(y+1.);
            float denominator = log(concentration + 1.) - concentration/(concentration + 1.);

            float res = pow(10., HaloMass) * numerator/denominator;


                    //4. * PI * pow(10., logrhos) * pow(r, 3.) * (log(1. + concentration) - concentration/(1. + concentration));

            mass += res;
        }
        else
        {
            assert(false || assert_msg("Dark Matter profile \"" << profile_name << "\" not recognized"));
        }

    }

    return mass;

}

float Galaxy::K_Kernel_DW(float u)
{
    // Kernel K(U), equation 9.
    if(u == 0.) return 0.;

    float prefactor = 0.5 * pow(u, (2.* beta - 1.));

    float term1 = ((3./2.) - beta) * pow(PI, 0.5) * gamma_term;

    float term2 = beta * incompleteBeta(beta + 0.5, 0.5, 1./(u*u));
    float term3 = -incompleteBeta(beta - 0.5, 0.5, 1./(u*u));

    float res = prefactor * (term1 + term2 + term3);

    if((res < 0.) && (res < zero_perturbation))
    {
        res = zero_perturbation;
    }

    assert((res >= 0.) || assert_msg("K_Kernel about to return negative (" << res << ")." << std::endl <<
       "-----beta = " << beta << std::endl <<
       "-----u = " << u << std::endl <<
       "-----gamma term = " << gamma_term << std::endl <<
       "res = prefactor * (term1 + term2 + term3)" << std::endl <<
       "-----prefactor: ((3./2.) - beta) * pow(PI, 0.5) * gamma_term = " << prefactor << std::endl <<
       "-----term 1: ((3./2.) - beta) * pow(PI, 0.5) * gamma_term = " << term1 << std::endl <<
       "-----term 2: beta * incompleteBeta(beta + 0.5, 0.5, 1./(u*u)) = " << term2 << std::endl <<
       "-----term 3: -incompleteBeta(beta - 0.5, 0.5, 1./(u*u)) = " << term3 << std::endl));


    return res;
}

float Galaxy::sigma_integrand(float r)
{
    // Internals of sigma integrand
    if(R == 0) R = zero_perturbation;

    float ratio = r/R;

    assert(ratio >= 1.0  || assert_msg("sigma_integrand ratio r/R < 0 (value is " << ratio <<
        "). This will violate the incomplete beta function." << std::endl <<
         "-----r = " << r << std::endl <<
         "-----R = " << R << std::endl <<
         "-----HLR = " << half_light_radius << std::endl));

    float Kernel = this->K_Kernel_DW(ratio);
    float density = this->rho(r);
    //float Mass = this->cumulative_mass(r);
    float Mass = this->analytic_mass(r);
    float res;
    if(r == 0) r = zero_perturbation;
    res = Kernel * density * Mass /  r;

    assert((!isnan(res) && !isinf(res)) || assert_msg("Sigma LOS integrand (eq8) about to return " << res << std::endl <<
        "-----r = " << r << std::endl <<
        "-----R = " << R << std::endl <<
        "-----Kernel K(U) = " << Kernel << std::endl <<
        "-----Density rho(r) = " << density << std::endl <<
        "-----Mass M(r) = " << Mass << std::endl));

    if(res < 0)
    {
        #pragma omp critical
        assert((res >= 0.) || assert_msg("Sigma Los integrand (eq8) about to return negative (" << res << ")." << std::endl <<
                                                                                                "-----r = " << r << std::endl <<
                                                                                                "-----R = " << R << std::endl <<
                                                                                                "-----Kernel K(U) = " << Kernel << std::endl <<
                                                                                                "-----Density rho(r) = " << density << std::endl <<
                                                                                                "-----Mass M(r) = " << Mass << std::endl));
    }



    return res;
}

float Galaxy::sigma_los(float R_arg)
{
    assert( R_arg > 0 || assert_msg("R_arg cannot be < 0, got: " << R_arg << std::endl));
    // Full equation to calculate sigma LOS


	float upper_limit = 100.0 * half_light_radius + R_arg; // Hopefully a sufficently high value of r.

	R = R_arg;

	auto fp = bind(&Galaxy::sigma_integrand, this, _1);

	float prepass = 0.0; //SimpsonsRule(fp, R_arg, upper_limit, initial_subdiv);


	float los_accuracy = std::fmin(this->sigma_integrand(R_arg), pow(10., 8.));
    //float los_accuracy = this->sigma_integrand(R_arg); ///10.0; //pow(10., precision); //prepass/pow(10., precision + sigma_los_precision_modifier);

    //float numerator = 2. * GR * SimpsonsRule(fp, R_arg, upper_limit, 10000);

    float integral_term = AdaptiveRichardsonExtrapolate(fp, R_arg, upper_limit, los_accuracy);

	float numerator = 2. * GR *  integral_term;
	float denominator = this->MassDensity(R_arg);
	float value = pow(numerator/denominator, 0.5);


	if(isnan(value) || isinf(value))
	{
    #pragma omp critical
	    {
            assert( (!isnan(value) && !isinf(value)) ||
                    assert_msg(std::endl << "Sigma LOS (eq7) about to return " << value << std::endl <<
                 "-----formula = pow(numerator/denominator, 0.5)" << std::endl <<
                 "-----numerator = " << numerator << std::endl <<
                 "-----denominator = " << denominator << std::endl <<
                 "This occured when:" << std::endl <<
                 "-------R (distance, main argument) = " << R_arg << " kpc" << std::endl <<
                 "-------Integration was performed between " << R_arg << " and " << upper_limit << std::endl <<
                 "-------Integration prepass (simpsons rule with " << initial_subdiv << " elements) return = " << prepass << std::endl <<
                 "-------Accuracy = " << los_accuracy << std::endl <<
                 "-------Full Integration (ARE) = " << integral_term << std::endl));

	    }
	}




	return value;
}

float Galaxy::sigma_ap_integrand(float R_arg)
{
    // Integrand of sigma aperture
	float sigma = this->sigma_los(R_arg);
	float MDP = this->MassDensity(R_arg);
	return MDP * sigma * sigma * R_arg;
}

float Galaxy::sigma_ap(void)
{
    // Full value of halo aperture
	auto numfp = bind(&Galaxy::sigma_ap_integrand, this, _1);

	float accuracy = pow(10, 15-precision);

    float num_accuracy = SimpsonsRule(numfp, 0., aperture_size, initial_subdiv)/pow(10., precision);
  	float numerator = AdaptiveRichardsonExtrapolate(numfp, 0., aperture_size, num_accuracy);

	auto denfp = bind(&Galaxy::MassDensityr, this, _1);

	float den_accuracy = SimpsonsRule(denfp, 0., aperture_size, initial_subdiv)/pow(10., precision);
	float denominator = AdaptiveRichardsonExtrapolate(denfp, 0., aperture_size, den_accuracy);

	return pow(numerator/denominator, 0.5);	
}

void Galaxy::GetHaloC(bool scatter)
{
    // read-off c-Mh from Benedikt's files (taken from http://www.benediktdiemer.com/data/)

    std::vector<int> * IndexesToGrab = new std::vector<int>(0);
    IndexesToGrab->push_back(0); // z
    IndexesToGrab->push_back(2); // M200c
    IndexesToGrab->push_back(4); // c200c

    std::vector<std::vector<float>> * Extracted;
    Extracted = ReadFile(Conc_Path, IndexesToGrab);

    std::vector<float> * Redshift = &Extracted->at(0);
    std::vector<float> * M200c = &Extracted->at(1);
    std::vector<float> * c200c = &Extracted->at(2);

    std::vector<float> * reduced = new std::vector<float>;
    float closest_z = FindClosest(redshift,  Reduce(Redshift, reduced));
    std::vector<bool> * mask = new std::vector<bool>;
    Equals(Redshift, closest_z, mask);

    MaskOut(Redshift, mask);
    MaskOut(M200c, mask);
    MaskOut(c200c, mask);

    // Find the concentration by linearly interpolating the data.

    float logC = LinearInterp(M200c, c200c, pow(10., HaloMass + log10(h)));

    float scatter_magnitude = 0.;

    if(scatter)
    {
	std::default_random_engine generator;
	std::normal_distribution<float> distribution(0., 1.0);
	scatter_magnitude = distribution(generator) * 0.16;
    }
    concentration = logC + scatter_magnitude; //pow(10., logC + scatter_magnitude);

    // Prevent Memory Leaks
    delete IndexesToGrab;
    delete Extracted;
    //delete Redshift;
    //delete M200c;
    //delete c200c;
    //delete mask;

}

void Galaxy::GetHaloR(void)
{
    float OmegaL = 1. - Om;
    float Omegaz = Om*pow((1.+redshift), 3.)/(Om*pow((1.+redshift), 3.) + OmegaL);
    float d = Omegaz-1.;
    float Deltac = 18.*pow(PI, 2.) + 82.*d - 39.*d*d;
    float H0 = 100.*h; //km/s/Mpc
    float HH = H0*pow(Om*pow((1.+redshift), 3.)+(1.-Om), 0.5);
    float GR_Mpc = 4.302e-9;  // Mpc/Msun(km/s)^2
    float rhoc =3.*HH*HH/(8.*PI*GR_Mpc);
    float k = 4.*PI/3.;
    HaloRadius = pow((pow(10., HaloMass)/rhoc/k/200.), (1./3.))*1000.; // kpc
}

float Galaxy::HaloDensity(float r)
{
    float res;

    if(r == 0)
    {
        r = zero_perturbation;
    }

    if (profile_name == "NFW")
    {
        float log10r = log10(r);
        float rs = HaloRadius/concentration;

	    float fc = log(1.+concentration)-concentration/(1.+concentration);
  	    float logrhos = HaloMass-log10(4.*PI*rs*rs*rs*fc);

	    float logrho = logrhos-log10(r/rs) - 2.*log10(1. + r/rs);
	    res = pow(10., logrho);

    }
    else if (profile_name == "DenhenMcLaughlin")
    {
        float log10r = log10(r);
        float rs = HaloRadius/concentration;
        float fc = log(1.+concentration)-concentration/(1.+concentration);
        float logrhos = HaloMass-log10(4.*PI*rs*rs*rs*fc);
        //std::cout << "p0 " << logrhos << std::endl;
        float logrho = log10(2.) + 6.*logrhos - (7./9.)*log10(r/rs) - 6.*log10(1. + pow(r/rs, 4./9.));
        //std::cout << "Lohrho" << logrho << std::endl;
        //exit(1);
        res = pow(10., logrho);
    }
    else if(profile_name == "Burkert")
    {
        // following Cattaneo et al. 2014, I use the same exact rho0 and rs as NFW but different shape:
        float fc = (log(1.+concentration)-concentration/(1.+concentration));
        float rs = HaloRadius/concentration;
        float R0 = log10(rs);
        float xx = r/rs;
        float rho0 = HaloMass - log10(4.*PI*pow(rs, 3.)*fc);
        float logrho = rho0 - log10(1.+xx) - log10(1.+xx*xx);
        res = pow(10., logrho);
    }
    else if(profile_name == "Exponential")
    {
        float xr=r/HaloRadius;
        float logsig0 = HaloMass - log10(2.*PI) - 2.*log10(HaloRadius);
        float expo = exp(-xr);
        float logrho = logsig0 + log10(expo);
        res = pow(10., logrho);
        //std::cout << "Res: " << res << std::endl;
    }
    else if(profile_name == "Hernquist")
    {
        float a= HaloRadius/(1.+sqrt(2.));
        float logrho = HaloMass + log10(a)-log10(2.*PI)-log10(r)-3.*log10(HaloRadius+a);
        res = pow(10., logrho);
    }
    else if(profile_name == "CoredNFW")
    {
        float rs = HaloRadius/concentration;
        float rcore = pow(10., 1.14); //kpc from Newman et al. 2013 paperII
        float b = rs/rcore;
        float fc = log(1. + concentration)-concentration/(1.+concentration);
        float logrhos = HaloMass-log10(4.*PI*pow(rs, 3.)*fc);
        float logrho = log10(b)+logrhos-log10(1.+b*HaloRadius/rs)-2.*log10(1.+HaloRadius/rs);
        res = pow(10.,logrho);
    }
    else
    {
	    std::cout << "HaloDensity() did not recognise profile name: " << profile_name << std::endl;
	    exit(1);
    }

    if(isnan(res) || isinf(res) )
    {
        std::cout << "HaloDensity(r) will return " << res << " when r = " << r << ". Profile = " << profile_name << std::endl;
        exit(1);
    }

    return res;

}

void Galaxy::setDarkMatter(float InputHaloMass, std::string name)
{
	dark_matter_on = true;
	HaloMass = InputHaloMass;
	profile_name = name;
	GetHaloC(false);
	GetHaloR();
}

void Galaxy::setConc_Path(std::string input_path)
{
    Conc_Path = input_path;
}

void Galaxy::set_stars_on(bool new_value) {
    stars_on = new_value;
}


float GetVelocityDispersion(float input_aperture_size,
                            float input_beta,
                            float input_half_light_radius,
                            float input_sersic_index,
                            float input_stellar_mass,
			                float z)
{
	
    Galaxy aGalaxy(input_stellar_mass, input_beta, input_half_light_radius, input_aperture_size, input_sersic_index, z);
    float res =  aGalaxy.sigma_ap();
    return res;
}

float GetVelocityDispersion(float input_aperture_size,
                            float input_beta,
                            float input_half_light_radius,
                            float input_sersic_index,
                            float input_stellar_mass,
                            float z,
                            float halo_mass,
                            char * profile_name,
                            char * c_path)
{
	Galaxy aGalaxy(input_stellar_mass, input_beta, input_half_light_radius, input_aperture_size, input_sersic_index, z);

	std::string path(c_path);
    aGalaxy.setConc_Path(path);

    std::string name(profile_name);
	aGalaxy.setDarkMatter(halo_mass, name);

	return aGalaxy.sigma_ap();
}

float GetUnweightedVelocityDispersion(float R,
                                      float input_beta,
                                      float input_half_light_radius,
                                      float input_sersic_index,
                                      float input_stellar_mass,
                                      float z)
{
    Galaxy aGalaxy(input_stellar_mass, input_beta, input_half_light_radius, 1.0, input_sersic_index, z);
    float res = aGalaxy.sigma_los(R);
    assert(!isnan(res) && "GetUnweightedVelocityDispersion() returned NaN");
    assert(!isinf(res) && "GetUnweightedVelocityDispersion() returned inf");
    return res;
}

float GetUnweightedVelocityDispersion(float R,
                                     float input_beta,
                                     float input_half_light_radius,
                                     float input_sersic_index,
                                     float input_stellar_mass,
                                     float z,
                                     float halo_mass,
                                     char * profile_name,
                                     char * c_path)
{
    Galaxy aGalaxy(input_stellar_mass, input_beta, input_half_light_radius, 1.0, input_sersic_index, z);
    std::string path(c_path);
    aGalaxy.setConc_Path(path);

    std::string name(profile_name);
    aGalaxy.setDarkMatter(halo_mass, name);
    float res = aGalaxy.sigma_los(R);
    assert(!isnan(res) && "GetUnweightedVelocityDispersion() returned NaN");
    assert(!isinf(res) && "GetUnweightedVelocityDispersion() returned inf");
    return res;
}

float GetUnweightedVelocityDispersion(float R,
                                      float z,
                                      float halo_mass,
                                      char * profile_name,
                                      char * c_path)
{
    Galaxy aGalaxy(11., 0.1, 5, 1.0, 4, z);
    std::string path(c_path);
    aGalaxy.setConc_Path(path);

    std::string name(profile_name);
    aGalaxy.setDarkMatter(halo_mass, name);
    aGalaxy.set_stars_on(false);
    float res = aGalaxy.sigma_los(R);
    assert(!isnan(res) && "GetUnweightedVelocityDispersion() returned NaN");
    assert(!isinf(res) && "GetUnweightedVelocityDispersion() returned inf");
    return res;
}


float GetDMrho(float R,
              float input_beta,
              float input_half_light_radius,
              float input_sersic_index,
              float input_stellar_mass,
              float z,
              float halo_mass,
              char * profile_name,
              char * c_path)
{
    Galaxy aGalaxy(input_stellar_mass, input_beta, input_half_light_radius, 1.0, input_sersic_index, z);
    std::string path(c_path);
    aGalaxy.setConc_Path(path);
    std::string name(profile_name);
    aGalaxy.setDarkMatter(halo_mass, name);

    float res = aGalaxy.HaloDensity(R);
    return res;
}

float GetCum(float R,
               float input_beta,
               float input_half_light_radius,
               float input_sersic_index,
               float input_stellar_mass,
               float z,
               float halo_mass,
               char * profile_name,
               char * c_path)
{
    Galaxy aGalaxy(input_stellar_mass, input_beta, input_half_light_radius, 1.0, input_sersic_index, z);
    std::string path(c_path);
    aGalaxy.setConc_Path(path);
    std::string name(profile_name);
    aGalaxy.setDarkMatter(halo_mass, name);

    //float res = aGalaxy.sigma_integrand(R);

    //auto fp = bind(&Galaxy::sigma_integrand, aGalaxy, _1);

    float res = aGalaxy.sigma_integrand(R);

    return res;
}


