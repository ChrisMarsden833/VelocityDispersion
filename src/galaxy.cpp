#include "galaxy.h"

Galaxy::Galaxy(float input_stellar_mass,
	       float input_beta,
	       float input_half_light_radius,
	       float input_aperture_size,
	       float input_sersic_index,
	       float z)
{
	stellar_mass = input_stellar_mass;

    beta = input_beta;

	if(beta == 0.0) beta = zero_perturbation;
	if(beta == 0.5) beta += zero_perturbation;

    gamma_term = (float) boost::math::tgamma(beta - 0.5) / boost::math::tgamma(beta);

	half_light_radius = input_half_light_radius;
	sersic_index = input_sersic_index;
	aperture_size = input_aperture_size;
	redshift = z;
	
	dark_matter_on = false;

	if(stellar_mass > 20)
	{
		#pragma omp critical
		{
			std::cout << "Error, stellar mass is unexpectedly high (" << stellar_mass << ") - are the units log10[M_sun]?" << std::endl;
			exit(1);
		}	
	}

	// Effective constants for divide by zero errors
	sersic_index_eff = sersic_index;
	if(sersic_index == 0.) sersic_index_eff = zero_perturbation;
	half_light_radius_eff = half_light_radius;
	if(half_light_radius == 0.) half_light_radius_eff = zero_perturbation;
	
	// Set the value of b_n.	
	b_n = 2.*sersic_index_eff - (1./3.) + (.009876/sersic_index_eff); 

	// find the gamma function used immedately in sigma_e. This is used by my defn of sigma_e.
	// not sure this is required?
	float gamma = boost::math::tgamma_lower(2*sersic_index, b_n);

	// Sigma_e, the constant for the sersic profile.
	sigma_e = pow(10, stellar_mass)/( pow(half_light_radius_eff, 2) * PI * 2 * 2. * sersic_index * exp(b_n) * gamma / pow(b_n, 2*sersic_index) );

	if(isnan(sigma_e) || isinf(sigma_e))
	{
		#pragma omp critical
		{
			std::cout << "Error: sigma_e = " << sigma_e << std::endl;
			std::cout << "   Numerator (stellar mass, m_sun) = " << pow(10, stellar_mass) << std::endl;
			std::cout << "   Denominator (2*PI*r^2) = " << ( pow(half_light_radius_eff, 2) * PI * 2.) << std::endl; 
			exit(1);
		}
	}


	//* 2.*sersic_index*exp(b_n)*gamma/pow(b_n, 2.*sersic_index)); 
	// I'm not sure if this is correct. Literature makes this out to be just sm/(pi * re^2)

	// p_n - eqn (5)
	p_n = 1. - .6097/sersic_index_eff  + .00563/(sersic_index_eff*sersic_index_eff);

	// Rho0 - eqn (4)
	float Sigma0 = this->MassDensity(0.);
	float left = Sigma0 *  pow(b_n, (sersic_index * (1. - p_n))) / (2.*half_light_radius_eff);
	float right = (boost::math::tgamma(2*sersic_index)/boost::math::tgamma(sersic_index*(3. - p_n)));
	rho0 = left * right;

	if(isnan(rho0) || isinf(rho0))
	{
		#pragma omp critical
		{
			std::cout << "Error: rho0 (equation 4) = " << rho0 << std::endl;
			std::cout << "    left = " << left << std::endl;
			std::cout << "    right = " << right << std::endl;
			exit(1);
		}
	}

	R = 4.;

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

	if(isnan(result) || isinf(result))
	{
		#pragma omp critical
		{
			std::cout << "Error: MassDensity(r) returned " << result << " at r = " << r << std::endl;
			std::cout << "    sigma_e =  " << sigma_e << std::endl;
			std::cout << "    internal_term = " << internal_term << std::endl;
			std::cout << "    power = " <<  power << std::endl;
			exit(1);
		}
	}

	return result;
}

float Galaxy::MassDensityr(float r)
{
	return this->MassDensity(r) * r;
}


float Galaxy::rho(float r)
{

    float radius_ratio = r/half_light_radius_eff;

    float radius_ratio_eff = radius_ratio;
    if(radius_ratio == 0. && -p_n <= 0) radius_ratio_eff = zero_perturbation;
    float power_term = pow(radius_ratio_eff, -p_n);

    float inverse_n = 1./sersic_index_eff;
    float exp_power_term = pow(radius_ratio_eff, inverse_n);

    float exp_term = exp(-b_n * exp_power_term);
    
    float dark_matter_term = 0.0;

    if(dark_matter_on)
    {
	    dark_matter_term = HaloDensity(r);
	    if(isnan(dark_matter_term) || isinf(dark_matter_term))
	    {
		#pragma omp critical
		{
	        	std::cout << "HaloDensity(r) returned " << dark_matter_term << ", at r = " << r << std::endl;
	        	exit(1);
		}
	    }
    }

    float res = rho0 * power_term * exp_term + dark_matter_term;
    if(isnan(res) || isinf(res))
    {
	#pragma omp critical
	{
        	std::cout << "Error: Density profile returned " << res << " at r = " << r << std::endl;
	        std::cout << "    rho0 = " << rho0 << std::endl;
		std::cout << "    power_term = " << power_term << std::endl;
		std::cout << "    exp_term = " << exp_term  << std::endl;
		std::cout << "    Dark_Matter_term = " << dark_matter_term << std::endl;
		exit(1);
	}
    }

    return res;
}

float Galaxy::mass_shell(float R)
{
	float res = 4.0*PI *R*R * this->rho(R);

    if(isnan(res))
    {
        std::cout << "mass_shell returned NaN at R =" << R << std::endl;
        exit(1);
    }

	return res;
}


float Galaxy::cumulative_mass(float R_arg)
{
	float accuracy = pow(10., stellar_mass-precision);

	auto fp = bind(&Galaxy::mass_shell, this, _1);

    float mass_accuracy = SimpsonsRule(fp, 0., R_arg, initial_subdiv)/pow(10., precision + cum_mass_precision_modifier);

	float res = AdaptiveRichardsonExtrapolate(fp, 0., R_arg, mass_accuracy);

    if(isnan(res))
    {
        std::cout << "cumulative_mass returned NaN at R_arg =" << R_arg << std::endl;
        exit(1);
    }

	return res;
}

float Galaxy::K_Kernel_DW(float u)
{
    if(u == 0.)
    {
        return 0.;
    }
    float prefactor = 0.5 * pow(u, (2.* beta - 1.));

    float term1 = ((3./2.) - beta) * pow(PI, 0.5) * gamma_term;

    float term2 = beta * incompleteBeta(beta + 0.5, 0.5, 1./(u*u));
    float term3 = -incompleteBeta(beta - 0.5, 0.5, 1./(u*u));

    float res = prefactor * (term1 + term2 + term3);

    return res;
}

float Galaxy::sigma_integrand(float r)
{ 
    if(R == 0) R = zero_perturbation;

    float ratio = r/R;

    float Kernal = this->K_Kernel_DW(ratio);
    float density = this->rho(r);
    float Mass = this->cumulative_mass(r);
	
    float res;

    if(r == 0) r = zero_perturbation;

    res = Kernal * density * Mass /  r;

    return res;
}


float Galaxy::sigma_los(float R_arg)
{
	float upper_limit = 1000*R_arg;
	R = R_arg;

	auto fp = bind(&Galaxy::sigma_integrand, this, _1);

    float los_accuracy = SimpsonsRule(fp, R_arg, upper_limit, initial_subdiv)/pow(10., precision + sigma_los_precision_modifier);

	float numerator = 2. * GR *  AdaptiveRichardsonExtrapolate(fp, R_arg, upper_limit, los_accuracy);
	float denominator = this->MassDensity(R_arg);

	float value = pow(numerator/denominator, 0.5);

	return value;
}


float Galaxy::sigma_ap_integrand(float R_arg)
{
	float sigma = this->sigma_los(R_arg);
	float MDP = this->MassDensity(R_arg);

	return MDP * sigma * sigma * R_arg;
}


float Galaxy::sigma_ap(void)
{
	auto numfp = bind(&Galaxy::sigma_ap_integrand, this, _1);

	float accuracy = pow(10, 15-precision);

    float num_accuracy = SimpsonsRule(numfp, 0., aperture_size, initial_subdiv)/pow(10., precision);
  	float numerator = AdaptiveRichardsonExtrapolate(numfp, 0., aperture_size, num_accuracy);

	auto denfp = bind(&Galaxy::MassDensityr, this, _1);

	float den_accuracy = SimpsonsRule(denfp, 0., aperture_size, initial_subdiv)/pow(10., precision);
	float denominator = AdaptiveRichardsonExtrapolate(denfp, 0., aperture_size, accuracy);
	

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
        std::cout << "p0 " << logrhos << std::endl;
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


