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
	//float gamma = boost::math::tgamma_lower(2*sersic_index, b_n);

	// Sigma_e, the constant for the sersic profile.
	sigma_e = pow(10, stellar_mass)/( pow(half_light_radius_eff, 2) * PI * 2.);
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

	float res = AdaptiveRichardsonExtrapolate(fp, 0., R_arg, accuracy);

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
    float Mass = this->cumulative_mass(r);
	
    float res;

    if(r == 0) r = zero_perturbation;

    res = Kernal * density * Mass /  r;

    return res;
}


float Galaxy::sigma_los(float R_arg)
{
	float upper_limit = 100*R_arg;
	R = R_arg;

	float accuracy = pow(10., 17 - precision);

	auto fp = bind(&Galaxy::sigma_integrand, this, _1);

	float numerator = 2. * GR *  AdaptiveRichardsonExtrapolate(fp, R_arg, upper_limit, accuracy); 
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

  	float numerator = AdaptiveRichardsonExtrapolate(numfp, 0., aperture_size, accuracy);

	auto denfp = bind(&Galaxy::MassDensityr, this, _1);

	float denominator = AdaptiveRichardsonExtrapolate(denfp, 0., aperture_size, accuracy);
	

	return pow(numerator/denominator, 0.5);	
}

void Galaxy::GetHaloC(bool scatter)
{
    // read-off c-Mh from Benedikt's files (taken from http://www.benediktdiemer.com/data/)
    std::string FilePath = "../data/cM_planck18.txt";

    std::vector<int> * IndexesToGrab = new std::vector<int>(0);
    IndexesToGrab->push_back(0); // z
    IndexesToGrab->push_back(2); // M200c
    IndexesToGrab->push_back(4); // c200c

    std::vector<std::vector<float>> * Extracted;
    Extracted = ReadFile(FilePath, IndexesToGrab);

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
    if(r == 0)
    {
        r = zero_perturbation;
    }

    if (profile_name == "NFW")
    {
        float log10r = log10(r);
        if(isnan(log10r) || isinf(log10r))
        {
            std::cout << "Halo Density(r), Log10(r) = " << log10r << " when r = " << r << std::endl;
            exit(1);
        }

        float rs = HaloRadius/concentration;
	    float fc = log(1.+concentration)-concentration/(1.+concentration);
  	    float logrhos = HaloMass-log10(4.*PI*rs*rs*rs*fc);
	    float logrho = logrhos-log10(pow(10., log10(r))/rs) - 2.*log10(1. + pow(10., log10(r))/rs);
	    float res = pow(10., logrho);

	    if(isnan(res) || isinf(res) )
	    {
	        std::cout << "HaloDensity(r) will return " << res << " when r = " << r << " log10(r) = " << log10r << std::endl;
	        exit(1);
	    }

	    return res;
    }
    else
    {
	    std::cout << "HaloDensity() did not recognise profile name: " << profile_name << std::endl;
	    exit(1);
    }

}

void Galaxy::setDarkMatter(float InputHaloMass, std::string name)
{
	dark_matter_on = true;
	HaloMass = InputHaloMass;
	profile_name = name;
	GetHaloC(false);
	GetHaloR();
	std::cout << "Halo C = " << concentration << std::endl;
	std::cout << "Halo R = " << HaloRadius << std::endl;
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
                            char * profile_name)
{
	Galaxy aGalaxy(input_stellar_mass, input_beta, input_half_light_radius, input_aperture_size, input_sersic_index, z);
	std::string name(profile_name);
	aGalaxy.setDarkMatter(halo_mass, name);

	return aGalaxy.sigma_ap();
}



