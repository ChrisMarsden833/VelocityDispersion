#include "galaxy.h"

Galaxy::Galaxy(float input_aperture_size, int * tracer_flags, int * gravitational_flags)
{
    // Constructor. This all happens when the class is initialized.
	aperture_size = input_aperture_size;
	trace_bulge = bool(tracer_flags[0]);
	grav_bulge = bool(gravitational_flags[0]);
    trace_disk = bool(tracer_flags[1]);
	grav_disk = bool(gravitational_flags[1]);
    grav_halo = bool(gravitational_flags[2]);
    grav_bh = bool(gravitational_flags[3]);
}

// ///////////////////////////
// ///// Bulge Functions /////
// ///////////////////////////

void Galaxy::ConstructBulge(float input_bulge_mass, float input_bulge_beta, float input_bulge_half_light_radius, float input_sersic_index)
{
    // Mass
    assert(input_bulge_mass < 30 || assert_msg("Bulge mass (" << input_bulge_mass << ") is unexpectedly high (>30) - are the units log10 [M_sun]?"));
    bulge_stellar_mass = input_bulge_mass;

    if(bulge_stellar_mass == 0.0)
    {
        // By definition
        trace_bulge = false;
        grav_bulge = false;
    }

    // Beta
    if( fmod(input_bulge_beta, 0.5) == 0.) input_bulge_beta += zero_perturbation;
    bulge_beta = input_bulge_beta;
    bulge_half_light_radius = input_bulge_half_light_radius;

    if((input_sersic_index < 0.13)) input_sersic_index = 0.13;
    if((input_sersic_index > 8.)) input_sersic_index = 8.;
    bulge_sersic_index = input_sersic_index;

    // Calculate the gamma term within the Kernel, often used and only depends on beta
    gamma_term = (float) boost::math::tgamma(bulge_beta - 0.5) / boost::math::tgamma(bulge_beta);

    // Set the value of b_n (D.W. - equation 2)
    b_n = 2.*bulge_sersic_index - (1./3.) + (.009876/bulge_sersic_index);

    float gamma = boost::math::tgamma_lower(2 * bulge_sersic_index, b_n);

    // Sigma_e, the constant for the sersic profile. It's debatable if this should just be M_*/(2*pi*Re^2) or this, but it seems to work as is:
    if(bulge_stellar_mass > 0.)
    {
        sigma_e = pow(10, bulge_stellar_mass) / (pow(bulge_half_light_radius, 2) * PI * 2.  * 2.  * bulge_sersic_index * exp(b_n) * gamma / pow(b_n, 2 * bulge_sersic_index) );
    }


    if (isnan(sigma_e) || isinf(sigma_e)){
        catastrophic_fail = true;
        #pragma omp critical
        {
            std::cout << "Warning: Sigma_e was calculated as: " << sigma_e << std::endl <<
                                  "====Location: Construction of Bulge" << std::endl <<
                                  "----Numerator (stellar mass, m_sun) = " << pow(10, bulge_stellar_mass) << std::endl <<
                                  "----Half Light Radius (Kpc) = " << bulge_half_light_radius << std::endl <<
                                  "----Sersic Index = " << bulge_sersic_index << std::endl <<
                                  "----b_n = " << b_n << std::endl << std::endl;
        }
    };

    // p_n - eqn (5)
    p_n = 1. - .6097/bulge_sersic_index  + .00563/(bulge_sersic_index*bulge_sersic_index);

    // Rho0 - eqn (4)
    Sigma0 = this->BulgeProjectedDensity(0.);
    float left = Sigma0 * pow(b_n, (bulge_sersic_index * (1. - p_n))) / (2. * bulge_half_light_radius);
    float right = (boost::math::tgamma(2 * bulge_sersic_index) / boost::math::tgamma(bulge_sersic_index * (3. - p_n)));
    rho0 = left * right;

    float as = bulge_half_light_radius / pow(p_n, bulge_sersic_index);
    float l = rho0 / pow(b_n, (1. - p_n));
    mass_prefactor = 4.0 * PI * bulge_sersic_index * boost::math::tgamma((3.0 - p_n) * bulge_sersic_index) * l * pow(as, 3.);

    // Asserts
    assert((!isnan(rho0) && !isinf(rho0)) || assert_msg("rho0 was calculated as: " << rho0 << std::endl <<
         "====Location: Construction of bulge" << std::endl <<
	 "----Bulge Mass = " << input_bulge_mass << std::endl <<
	 "----Bulge beta = " << input_bulge_beta << std::endl <<
	 "----Bulge HLR = " << input_bulge_half_light_radius <<
	 "----Bulge n = " << input_sersic_index << std::endl <<
         "----Sigma0 (BulgeMassDensity(0.)) = " << Sigma0 << std::endl <<
         "----LHS = " << left << std::endl <<
         "----RHS = " << right));
}

float Galaxy::BulgeProjectedDensity(float R)
{
    // Equation (1)
	float power = 1./bulge_sersic_index;
	float internal_term = R/bulge_half_light_radius;
	float result = sigma_e * exp(-b_n * (pow(internal_term, power) - 1.));
    	
	assert((!isnan(result) && !isinf(result)) || assert_msg("BulgeProjectedDensity (Sersic Profile) about to return " << result << std::endl <<
         "----Bulge n = " << bulge_sersic_index << std::endl <<
	 "----R = " << R << std::endl <<
       	 "----HLR = " << bulge_half_light_radius  << std::endl <<
	 "result = sigma_e * exp(-b_n * (pow(internal_term, power) - 1.))" << std::endl <<	 
	 "----Power (1/n) = " << power << std::endl <<
         "----internal term (r/R_eff) = " << internal_term << std::endl <<
         "----Sigma_e = " << sigma_e << std::endl));
	return result;
}

float Galaxy::BulgeProjectedDensityxR(float R)
{
    // Just a wrapper function, sometimes we need it *r for integrals.
    return this->BulgeProjectedDensity(R) * R;
}

float Galaxy::BulgeDensity(float r)
{
    // De-projected density, equation (3)
    float radius_ratio = r/bulge_half_light_radius;
    if(radius_ratio == 0. && -p_n <= 0) radius_ratio += zero_perturbation;
    float power_term = pow(radius_ratio, -p_n);
    float inverse_n = 1./bulge_sersic_index;
    float exp_power_term = pow(radius_ratio, inverse_n);
    float exp_term = exp(-b_n * exp_power_term);
    float res = rho0 * power_term * exp_term;

    if(res < 0.) res = 0;

    assert((!isnan(res) && !isinf(res)) || assert_msg("Bulge Density (eq3) to return " << res << std::endl <<
             "----rho0 = " << rho0 << std::endl <<
             "----Power term (r/R_eff)^-p_n = " << power_term << std::endl <<
             "----Exp power term (r/R_eff)^1/n = " << exp_power_term << std::endl <<
             "----Exp term exp(-b_n (Exp_power_term)) = " << exp_power_term << std::endl <<
             "-----r = " << r << std::endl));

    return res;
}

float Galaxy::BulgeMass(float r)
{
    // The total mass of stars in the bulge
    float as = bulge_half_light_radius/pow(b_n, bulge_sersic_index);
    float x = r/as;
    float threemp_term = (3.-p_n) * bulge_sersic_index;
    return pow(10., bulge_stellar_mass) * boost::math::tgamma_lower(threemp_term, pow(x, 1. / bulge_sersic_index) ) / boost::math::tgamma(threemp_term, 0.0);
}

float Galaxy::TotalMassForBulge(float r)
{
    float mass = 0.0;
    if(grav_bulge) mass += rem_prefac * this->BulgeMass(r); //1.72 * this->BulgeMass(r);
    if(grav_halo) mass += this->HaloMass(r);
    if(grav_bh)  mass += pow(10., BHMass);
    // if(grav_disk) mass += this->disk_mass(r);
    
    assert(mass >= 0 || assert_msg("Bulge Total Mass < 0, actual value: " << mass << std::endl <<
                                                          "---BulgeMass(r) = " << this->BulgeMass(r) << std::endl <<
                                                          "---DM(r) = " << this->HaloMass(r) << std::endl <<
                                                          "---Disk Mass(r) = " << this->disk_mass(r) << std::endl <<
                                                          "---r =" << r << std::endl
			    
			    ));
    
    return mass;
}


float Galaxy::K_Kernel_DW(float u)
{
    // Kernel K(U), equation 9.
    if(u == 0.) return 0.;

    if(bulge_beta > 0.5) bulge_beta = 0.5;

    float prefactor = 0.5 * pow(u, (2. * bulge_beta - 1.));
    float term1 = ((3./2.) - bulge_beta) * pow(PI, 0.5) * gamma_term;
    float term2 = bulge_beta * incompleteBeta(bulge_beta + 0.5, 0.5, 1. / (u * u));
    float term3 = -incompleteBeta(bulge_beta - 0.5, 0.5, 1. / (u * u));

    float res = prefactor * (term1 + term2 + term3);

    if((res < 0.) && (res < zero_perturbation)) res = zero_perturbation;

    assert((res >= 0.) || assert_msg("K_Kernel about to return negative (" << res << ")." << std::endl <<
       "-----beta = " << bulge_beta << std::endl <<
       "-----u = " << u << std::endl <<
       "-----gamma term = " << gamma_term << std::endl <<
       "res = prefactor * (term1 + term2 + term3)" << std::endl <<
       "-----prefactor: ((3./2.) - beta) * pow(PI, 0.5) * gamma_term = " << prefactor << std::endl <<
       "-----term 1: ((3./2.) - beta) * pow(PI, 0.5) * gamma_term = " << term1 << std::endl <<
       "-----term 2: beta * incompleteBeta(beta + 0.5, 0.5, 1./(u*u)) = " << term2 << std::endl <<
       "-----term 3: -incompleteBeta(beta - 0.5, 0.5, 1./(u*u)) = " << term3));
    return res;
}

float Galaxy::bulge_sigma_integrand(float r)
{
    // Internals of sigma integrand

    // Manage R, as this can be problematic if R = 0
    if(R == 0) R = zero_perturbation;
    float ratio = r/R;

    assert(ratio >= 1.0  || assert_msg("sigma_integrand ratio r/R < 0 (value is " << ratio <<
      "). This will violate the incomplete beta function." << std::endl <<
      "-----r = " << r << std::endl <<
      "-----R = " << R << std::endl <<
      "-----HLR = " << bulge_half_light_radius));

    float Kernel = this->K_Kernel_DW(ratio);
    float density = this->BulgeDensity(r);
    float Mass = this->TotalMassForBulge(r);
    float res;
    if(r == 0) r = zero_perturbation;

    res = Kernel * density * Mass /  r;

    assert((!isnan(res) && !isinf(res)) || assert_msg("Sigma LOS integrand (eq8) about to return " << res << std::endl <<
        "-----r = " << r << std::endl <<
        "-----R = " << R << std::endl <<
        "-----Kernel K(U) = " << Kernel << std::endl <<
        "-----Density Bulge_rho(r) = " << density << std::endl <<
        "-----Mass M(r) = " << Mass << std::endl));

    assert((res >= 0.) || assert_msg("Sigma Los integrand (eq8) about to return negative (" << res << ")." << std::endl <<
        "-----r = " << r << std::endl <<
        "-----R = " << R << std::endl <<
        "-----Kernel K(U) = " << Kernel << std::endl <<
        "-----Density Bulge_rho(r) = " << density << std::endl <<
        "-----Mass M(r) = " << Mass << std::endl));

    return res;
}

float Galaxy::bulge_sigma_los(float R_arg)
{
    float integral_term = 0., accuracy = 0., upper_limit = 0., los_accuracy = 0.;

    // Full equation to calculate sigma LOS
    if(R_arg == 0.) R_arg = zero_perturbation;

    float bulge_density = this->BulgeProjectedDensity(R_arg);
    if(bulge_density < zero_perturbation) bulge_density = zero_perturbation;

    float K = (2. * GR) / bulge_density;

	R = R_arg;
	auto fp = bind(&Galaxy::bulge_sigma_integrand, this, _1); // An unfortunate consequence of using classes

	if(slow_integrate)
	{
	    upper_limit = 100. * bulge_half_light_radius + R_arg;
	    integral_term = SimpsonsRule(fp, R_arg, upper_limit, 1000);
	}
    else
    {
        accuracy = 0.1; // kms^-1
        los_accuracy = pow(accuracy, 2)/K; // Transform this into the units of the integral

        // Let's find a sensible upper limit for the integral.
        float step = bulge_half_light_radius;
        upper_limit = R_arg + step; // We subdivide into
        float unit_accuracy = 2. * los_accuracy/(upper_limit - R_arg); // Unit accuracy is the accuracy divided by the interval
        while(unit_accuracy > los_accuracy/(upper_limit - R_arg))
        {
            unit_accuracy = this->bulge_sigma_integrand(upper_limit) / (upper_limit - R_arg) ;
            upper_limit += step;
        }
        integral_term = AdaptiveRichardsonExtrapolate(fp, R_arg, upper_limit, los_accuracy);
    }



	float value = pow(K * integral_term, 0.5);

	if(isnan(value) || isinf(value))
	{
    #pragma omp critical
	    {
            assert( (!isnan(value) && !isinf(value)) ||
                    assert_msg(std::endl << "Sigma LOS (eq7) about to return " << value << std::endl <<
                 "-----formula = pow(numerator/denominator, 0.5)" << std::endl <<
                 "-----K (2G/Sigma(r) = " << K << std::endl <<
                 "This occured when:" << std::endl <<
                 "-------R (distance, main argument) = " << R_arg << " kpc" << std::endl <<
                 "-------Bulge Mass " << bulge_stellar_mass << std::endl <<
                 "-------Bulge Density @ r " << bulge_density << std::endl <<
                 "-------Bulge Re " << bulge_half_light_radius << std::endl <<
                 "-------Integration was performed between " << R_arg << " and " << upper_limit << std::endl <<
                 "-------Accuracy = " << los_accuracy << std::endl <<
                 "-------Full Integration (ARE) = " << integral_term << std::endl));
	    }
	}
	return value;
}

float Galaxy::sigma_ap_integrand(float R_arg)
{
    // Integrand of sigma aperture
	float sigma = this->bulge_sigma_los(R_arg);
	float MDP = this->BulgeProjectedDensity(R_arg);
	return MDP * sigma * sigma * R_arg;
}

float Galaxy::bulge_sigma_ap(void)
{
    if(!trace_bulge) return 0.0;

    float numerator, denominator, accuracy;

    // Full value of halo aperture
	auto numfp = bind(&Galaxy::sigma_ap_integrand, this, _1);
    auto denfp = bind(&Galaxy::BulgeProjectedDensityxR, this, _1);

	if(false)
	{
        numerator = SimpsonsRule(numfp, 0, aperture_size, 1000);
	    denominator = SimpsonsRule(denfp, 0, aperture_size, 1000);
	}
	else
    {
        float multiplyer = 1e6;
        int prepass = 12;
        // Do a prepass to get a sense of the value of the integral
        accuracy = SimpsonsRule(numfp, 0, aperture_size, prepass)/multiplyer;
        numerator =  AdaptiveRichardsonExtrapolate(numfp, 0., aperture_size, accuracy);
        accuracy = SimpsonsRule(denfp, 0, aperture_size, prepass)/multiplyer;
        denominator = AdaptiveRichardsonExtrapolate(denfp, 0., aperture_size, accuracy);
    }


    if(denominator == 0.) return 0.;

    float res = pow(numerator/denominator, 0.5);

    if(isnan(res) || isinf(res))
    {
    #pragma omp critical
        {
            assert( (!isnan(res) && !isinf(res)) ||
                    assert_msg(std::endl << "Sigma ap (eq10) about to return " << res << std::endl <<
                                         "-----formula = pow(numerator/denominator, 0.5)" << std::endl <<
                                         "-----numereator: " << numerator <<
                                         "-----denominator: " << denominator <<
                                         "This occured when:" << std::endl <<
                                         "-------Aperture = " << aperture_size << " kpc" << std::endl <<
					 "This occured when:" << std::endl <<
					 "-------Bulge Mass " << bulge_stellar_mass << std::endl <<
			 	         "-------Bulge Re " << bulge_half_light_radius << std::endl ));
        }
    }

    return res;

}


// //////////////////////////
// ///// Halo Functions /////
// //////////////////////////

void Galaxy::ConstructHalo(float haloRs, float haloRhos, const char * halo_profile)
{
    HaloRadius = haloRs;
    HaloRhos = haloRhos;
    halo_profile_name = halo_profile;
}

float Galaxy::HaloDensity(float r)
{
    float density = 0.;
    if(strcmp(halo_profile_name, "Burkert") == 0.)
    {
        density = ( HaloRhos * pow(HaloRadius, 3) ) / ((r + HaloRadius) * (r*r + HaloRadius * HaloRadius));
    }
    else if(strcmp(halo_profile_name, "cNFW") == 0)
    {
        float rc = 5.0; // kpc
        density = ( HaloRhos ) / ((r/HaloRadius + rc/HaloRadius) * (1 + HaloRadius/rc) * (1 + HaloRadius/rc));
    }
    else{
        assert(false >= 0 || assert_msg(" Unknown halo profile type: " << halo_profile_name << std::endl));
    }

    return density;
}

float Galaxy::HaloDensityIntegrand(float r)
{
    return HaloDensity(r) * 4 * PI * r * r;
}

float Galaxy::HaloMass(float r)
{

    if(strcmp(halo_profile_name, "NFW") == 0)
    {
        float x = r/HaloRadius;
        float res = 4. * PI * pow(HaloRadius, 3.) * HaloRhos * (log(1. + x) - x/(1+x));

        if(x/(1+x) > log(1. + x)) res = 0.0; // At very low radii, can return negative.

        if (isinf(abs(res)) || isnan(res))
        {
            assert(false || assert_msg("Halo Analytic Mass, NFW profile, returned " << res << std::endl <<
                                                                                    "---x (r/HaloRadius) = " << x << std::endl <<
                                                                                    "---- r = " << r << std::endl <<
                                                                                    "---- HaloRadius = " << HaloRadius << std::endl <<
                                                                                    "---HaloRhos = " << HaloRhos << std::endl));
        }

        assert(res >= 0 || assert_msg(" Halo Analytic Mass < 0"));
        return res;

    }
    else if(strcmp(halo_profile_name, "Burkert") == 0)
    {
        float x = r/HaloRadius;

        float tanterm = -2.*atan(x);
        float term2 =  2.*log(1+x);
        float term3 = log(1 + x*x);

        float res = PI * HaloRhos * pow(HaloRadius, 3.) * (tanterm + term2 + term3);

        if(res < 0) res = 0.0;

        return res;
        /*
        // Any function that has not been expressed in analytic form
        auto numfp = bind(&Galaxy::HaloDensityIntegrand, this, _1);
        float min = 0.;
        float multiplyer = 1e3;
        int prepass = 10;
        float accuracy = SimpsonsRule(numfp, min, r, prepass)/multiplyer;
        float integral_term = AdaptiveRichardsonExtrapolate(numfp, min, r, accuracy);
        return integral_term;
         */
    }
    else{
        assert(false >= 0 || assert_msg(" Unknown halo profile type: " << halo_profile_name << std::endl));
        return 0.0;
    }


}

float Galaxy::HaloVcirc2(float r)
{
    return GR * HaloMass(r) / r;
}

// //////////////////////////
// ///// Disk Functions /////
// //////////////////////////

void Galaxy::ConstructDisk(float DiskMass, float input_inclination, float disk_scale)
{
    mass_disk = DiskMass;
    disk_inclination = input_inclination;
    disk_scale_length = disk_scale;
}

float Galaxy::disk_projected_density(float R)
{
    float mfactor = (1 - exp(-R/disk_scale_length)*(1+R/disk_scale_length));
    float sigma_0 = (pow(10., mass_disk)/(2. * PI * disk_scale_length * disk_scale_length)) * mfactor;
    return sigma_0 * exp(-R/disk_scale_length);  //(pow(10., mass_disk)/(2. * PI * disk_scale_length * disk_scale_length)) * exp(-R/disk_scale_length);
}

float Galaxy::disk_bessel(float x)
{
    return boost::math::cyl_bessel_i(0., x) * boost::math::cyl_bessel_k(0., x) -
                     boost::math::cyl_bessel_i(1., x) * boost::math::cyl_bessel_k(1., x);
}


float Galaxy::disk_Vcirc2(float R)
{
    if(R == 0.) R += zero_perturbation;
    float x = R/(2.*disk_scale_length);

    float B_term;

    if(x > 80)
    {
        B_term = 1.0; // Insane situations
    }
    else
    {
        B_term = disk_bessel(x); //disk_bessel(1.6);
    }

    float res = 0.0;

    if(grav_disk)
    {
        res += 2. * ((GR * pow(10., mass_disk))/(disk_scale_length)) * x * x * B_term;
    }

    if(grav_halo)
    {
        res += HaloVcirc2(R);
    }

    if(grav_bulge)
    {
        res += GR * BulgeMass(R)/R;
    }
    return res;
}

float Galaxy::disk_mass(float R)
{

    float mfactor = (1 - exp(-R/disk_scale_length)*(1+R/disk_scale_length));
    float sigma_0 = (pow(10., mass_disk)/(2. * PI * disk_scale_length * disk_scale_length)) * mfactor;

    // =Σ0exp (-R/RD)
    float res = 2 * PI * sigma_0 * disk_scale_length * disk_scale_length * mfactor;

            // sigma_0 * exp(-R/disk_scale_length);

    //std::cout << res << std::endl;

    //(pow(10., mass_disk)/(disk_scale_length*disk_scale_length)) * (disk_scale_length * (-exp(-R/disk_scale_length)) * (disk_scale_length + R) + (disk_scale_length*disk_scale_length) );

    if(res < 0. || isnan(res) || isinf(res) ) res = 0.;

    return res;
}

float Galaxy::disk_integrand(float R)
{
    return 2. * PI * R * disk_projected_density(R) * disk_Vcirc2(R) ;
}

float Galaxy::disk_velocity_dispersion2()
{
    if(!trace_disk || mass_disk == 0.0) return 0.0;


    auto numfp = bind(&Galaxy::disk_integrand, this, _1);

    float min = 0.;
    float multiplyer = 1e3;
    int prepass = 10;
    float accuracy = SimpsonsRule(numfp, min, aperture_size, prepass)/multiplyer;
    float integral_term = AdaptiveRichardsonExtrapolate(numfp, min, aperture_size, accuracy);
    //float integral_term = SimpsonsRule(numfp, 0., aperture_size, 10000);

    float mdisk = disk_mass(aperture_size);

    if(mdisk == 0.0)
    {
        return 0.0;
    }

    float constant = (sin(disk_inclination) * sin(disk_inclination))/mdisk;

    float res = constant * integral_term;

    assert( (!isnan(res) && !isinf(res)) ||
        assert_msg("disk_velocity_dispersion2 about to return: " << res << std::endl <<
                         "---- Constant " << constant << std::endl <<
                         "---- Integral Term " << integral_term << std::endl <<
                         "---- aperture_size : " << aperture_size << std::endl <<
                         "---- mass disk : " << mdisk << std::endl <<
                         "---- Disk_inclination: " << disk_inclination << std::endl));

    return res;

}

// //////////////////////////
// // Black Hole Functions //
// //////////////////////////

void Galaxy::Construct_Black_Hole(float bhMass)
{
    BHMass = bhMass;
}

bool Galaxy::getCatFail(void) {
    return catastrophic_fail;
}


// //////////////////////////
// //// Member Functions ////
// //////////////////////////

float GetVelocityDispersion(float Aperture,
                           float bulge_mass,
                           float bulge_radius,
                           float bulge_beta,
                           float bulge_sersicIndex,
                           float disk_mass,
                           float disk_inclination,
                           float disk_scale_length,
                           const char * halo_profile,
                           float haloRs,
                           float haloRhos,
                           float BlackHole_mass,
                           int * tracer_flags,
                           int * gravitational_flags,
                           int mode)
{
    if(Aperture == 0.0) return sqrt(-1);

    if(isnan(Aperture) || isnan(bulge_mass) || isnan(bulge_radius) || isnan(bulge_beta) || isnan(bulge_sersicIndex) ||
       isnan(disk_mass) || isnan(disk_inclination) || isnan(disk_scale_length) || isnan(haloRs) || isnan(haloRhos) ||
       isnan(BlackHole_mass)) return sqrt(-1);

    // Construct galaxy object
    Galaxy aGalaxy(Aperture, tracer_flags, gravitational_flags);

    aGalaxy.ConstructBulge(bulge_mass, bulge_beta, bulge_radius, bulge_sersicIndex);
    aGalaxy.ConstructDisk(disk_mass, disk_inclination, disk_scale_length);
    aGalaxy.ConstructHalo(haloRs, haloRhos, halo_profile);
    aGalaxy.Construct_Black_Hole(BlackHole_mass);

    if(aGalaxy.getCatFail()) return sqrt(-1);

    float res;

    if(mode == 1)
    {
        // Actually compute the disk and bulge velocity dispersion contributions (most of the legwork happens here)
        float bulge_sigma = 0;
        if(bulge_mass != 0)
        {
           bulge_sigma = aGalaxy.bulge_sigma_ap();
        }

        float disk_sigma2 = 0;
        if(disk_mass != 0)
        {
            disk_sigma2 = aGalaxy.disk_velocity_dispersion2();
        }




        // Combine these values in quadrature.
        res = pow(bulge_sigma*bulge_sigma + disk_sigma2, 0.5);

        // Stop in the case that inf or nan is returned, this should never have to happen.
        assert( (!isnan(res) && !isinf(res)) ||
                assert_msg(std::endl << "GetVelocityDispersion() about to return " << res << std::endl <<
                                     "----bulge_sigma: " << bulge_sigma << std::endl <<
                                     "----disk_sigma2: " << disk_sigma2 << std::endl <<
                                     "\t\t Disk_mass : " << disk_mass << std::endl <<
                                     "\t\t Disk_inclination: " << disk_inclination << std::endl));

    }
    else if(mode == 2)
    {
        // Return sigma LOS
        res = aGalaxy.bulge_sigma_los(Aperture);
    }
    else if(mode == 3)
    {
        return sqrt(aGalaxy.disk_Vcirc2(Aperture));
    }
    else
    {
        assert((false) || assert_msg( "Unknown mode : "  << mode));
    }
    return res;
}



