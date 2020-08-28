#include "galaxy.h"

Galaxy::Galaxy(float input_aperture_size, float z)
{
    // Constructor. This all happens when the class is initialized.
	aperture_size = input_aperture_size;
	redshift = z;
}

// ///////////////////////////
// ///// Bulge Functions /////
// ///////////////////////////

void Galaxy::ConstructBulge(float input_bulge_mass, float input_bulge_beta, float input_bulge_half_light_radius, float input_sersic_index)
{
    // Set all the bulge properties.
    bulge_present = true;

    // Mass
    assert(input_bulge_mass < 30 || assert_msg("Bulge mass (" << input_bulge_mass << ") is unexpectedly high (>30) - are the units log10 [M_sun]?"));
    bulge_stellar_mass = input_bulge_mass;

    // Beta
    if( fmod(input_bulge_beta, 0.5) == 0.) input_bulge_beta += zero_perturbation;
    bulge_beta = input_bulge_beta;

    // Radius
    assert(input_bulge_half_light_radius > 0. || assert_msg("The supplied bulge half light radius (" << input_bulge_half_light_radius << ") was not > 0"));
    bulge_half_light_radius = input_bulge_half_light_radius;

    // Sersic Index
    if(input_sersic_index == 0.) input_sersic_index += zero_perturbation;
    bulge_sersic_index = input_sersic_index;

    // Other bulge-related values that we only need to calculate once.

    // Calculate the gamma term within the Kernel, often used and only depends on beta
    gamma_term = (float) boost::math::tgamma(bulge_beta - 0.5) / boost::math::tgamma(bulge_beta);

    // Set the value of b_n (D.W. - equation 2)
    b_n = 2.*bulge_sersic_index - (1./3.) + (.009876/bulge_sersic_index);

    // Find the gamma function used in sigma_e.
    float gamma = boost::math::tgamma_lower(2 * bulge_sersic_index, b_n);
    // Sigma_e, the constant for the sersic profile. It's debatable if this should just be M_*/(2*pi*Re^2) or this, but it seems to work as is:
    sigma_e = pow(10, bulge_stellar_mass) / (pow(bulge_half_light_radius, 2) * PI * 2 * 2. * bulge_sersic_index * exp(b_n) * gamma / pow(b_n, 2 * bulge_sersic_index) );

    assert(!isnan(sigma_e) || !isinf(sigma_e) || assert_msg("Sigma_e was calculated as: " << sigma_e << std::endl <<
          "====Location: Construction of Bulge" << std::endl <<
          "----Numerator (stellar mass, m_sun) = " << pow(10, bulge_stellar_mass) << std::endl <<
          "----Half Light Radius (Kpc) = " << bulge_half_light_radius << std::endl <<
          "----Sersic Index = " << bulge_sersic_index << std::endl <<
          "----b_n = " << b_n));

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
    assert(!isnan(rho0) || !isinf(rho0) || assert_msg("rho0 was calculated as: " << rho0 << std::endl <<
         "====Location: Construction of bulge" << std::endl <<
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
    assert(!isnan(result) || !isinf(result) || assert_msg("BulgeProjectedDensity (Sersic Profile) about to return " << rho0 << std::endl <<
         "----Power (1/n) = " << power << std::endl <<
         "----internal term (r/R_eff) = " << internal_term << std::endl <<
         "----Sigma_e = " << sigma_e << std::endl <<
         "----RHS = " << right));
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

    assert(!isnan(res) || !isinf(res) || assert_msg("Bulge Density (eq3) to return " << res << std::endl <<
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

float Galaxy::BulgeTotalMass(float r)
{
    float mass = 0.0;
    if(bulge_stars_gravitation_on) mass += BulgeMass(r);
    if(dark_matter_gravitation_on) mass += this->HaloAnalyticMass(r);
    if(black_hole_gravitation_on)  mass += pow(10., BHMass);
    return mass;
}

void Galaxy::SetBulgeGravitationalContributions(bool stars, bool dark_matter, bool black_hole)
{
    bulge_stars_gravitation_on = stars;
    dark_matter_gravitation_on = dark_matter;
    black_hole_gravitation_on = black_hole;
}

float Galaxy::K_Kernel_DW(float u)
{
    // Kernel K(U), equation 9.
    if(u == 0.) return 0.;

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
    float Mass = this->BulgeTotalMass(r);
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
    if(R_arg == 0.) R_arg = zero_perturbation;

    // Full equation to calculate sigma LOS
	float upper_limit = 100.0 * bulge_half_light_radius + R_arg; // Hopefully a sufficently high value of r that can never be exceeded by R_arg
	R = R_arg;
	auto fp = bind(&Galaxy::bulge_sigma_integrand, this, _1); // An unfortunate consequence of using classes

	// The accuracy is the smallest of the first value in the integral * 0.01, or 10, 9. This combinationb seems to work well.
	// TODO - Generalize the accuracy in a phyiscal way that presevers performance
	float los_accuracy = std::fmin(this->bulge_sigma_integrand(R_arg) * 0.01, pow(10., 9.));
	float integral_term = AdaptiveRichardsonExtrapolate(fp, R_arg, upper_limit, los_accuracy);
	float numerator = 2. * GR *  integral_term;
	float denominator = this->BulgeProjectedDensity(R_arg);
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

float Galaxy::sigma_ap(void)
{
    float numerator, denominator, accuracy;

    // Full value of halo aperture
	auto numfp = bind(&Galaxy::sigma_ap_integrand, this, _1);

	float multiplyer = 1e3;
	int prepass = 6;
    // Do a prepass to get a sense of the value of the integral
    accuracy = SimpsonsRule(numfp, 0, aperture_size, prepass)/multiplyer;
    numerator = AdaptiveRichardsonExtrapolate(numfp, 0., aperture_size, accuracy);

	auto denfp = bind(&Galaxy::BulgeProjectedDensityxR, this, _1);
    accuracy = SimpsonsRule(denfp, 0, aperture_size, prepass)/multiplyer;
    denominator = AdaptiveRichardsonExtrapolate(denfp, 0., aperture_size, accuracy);

	return pow(numerator/denominator, 0.5);	
}

// //////////////////////////
// ///// Halo Functions /////
// //////////////////////////

void Galaxy::ConstructHalo(float input_halo_mass, std::string input_profile_name, std::string input_conc_path)
{
    HaloMass = input_halo_mass;
    profile_name = input_profile_name;
    Conc_Path = input_conc_path;
    GetHaloC(false);
    GetHaloR();
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

    // Prevent Memory Leaks? Probably not necessary...
    delete IndexesToGrab;
    delete Extracted;
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

    if(r == 0) r = zero_perturbation;

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

float Galaxy::HaloAnalyticMass(float r)
{
    float mass;
    if(profile_name == "NFW")
    {
        float rs = HaloRadius/concentration;
        float fc = log(1.+concentration)-concentration/(1.+concentration);
        float logrhos = HaloMass-log10(4.*PI*rs*rs*rs*fc);
        float prefactor = 4. * PI * pow(10., logrhos) * pow(rs, 3.);
        float ratio = r/rs;
        float rhs = log(1. + ratio) - r/(rs + r);
        float res = prefactor * rhs;
        if(abs(rhs) < zero_perturbation) res = 0.;
        assert(res >= 0. || assert_msg("NFW mass returned negative, res = " << res << std::endl <<
            "---- r = " << r << std::endl <<
            "---- rs = " << rs << std::endl <<
            "---- prefactor 4. * PI * pow(10., logrhos) * pow(rs, 3.) = " << prefactor << std::endl <<
            "---- rhs (log(rs + r) - log(rs) - r/(rs + r) =" << rhs << std::endl <<
            "----- ratio r/rs = " << ratio << std::endl));
        mass = res;
    }
    else assert(false || assert_msg("HaloAnalyticMass() -> Dark Matter profile \"" << profile_name << "\" not recognized/implemented"));
    return mass;
}

float Galaxy::HaloVcirc2(float r)
{
    return GR * HaloAnalyticMass(r) / r;
}

// //////////////////////////
// ///// Disk Functions /////
// //////////////////////////

void Galaxy::ConstructDisk(float DiskMass)
{
    disk_present = true;
    mass_disk = DiskMass;
    float length_log10 = 0.633 + 0.39*(mass_disk - 11.) + 0.069*(mass_disk- 11.)*(mass_disk- 11.);
    disk_scale_length = pow(10., length_log10);
}

float Galaxy::disk_projected_density(float R)
{
    return pow(10., mass_disk)/(2. * PI * disk_scale_length) * exp(-R/disk_scale_length);
}

float Galaxy::disk_Vcirc2(float R)
{
    float x = R/disk_scale_length;
    float B_term = boost::math::cyl_bessel_i(0., x/2.) * boost::math::cyl_bessel_k(0., x/2.) - boost::math::cyl_bessel_i(1., x/2.) * boost::math::cyl_bessel_k(1., x/2.);
    return (GR * pow(10., mass_disk))/(2.*disk_scale_length) * x * B_term;
}

float Galaxy::disk_mass(float R)
{
    return (pow(10., mass_disk)/(disk_scale_length*disk_scale_length)) * (-disk_scale_length * exp(-R/disk_scale_length) * (disk_scale_length + R) + (disk_scale_length*disk_scale_length) );
}

float Galaxy::disk_integrand(float R)
{
    return 2. * PI * R * disk_projected_density(R) * disk_Vcirc2(R);
}

float Galaxy::disk_velocity_dispersion2(float aperture_size, float inclination)
{
    auto numfp = bind(&Galaxy::disk_integrand, this, _1);

    float multiplyer = 1e3;
    int prepass = 6;
    float accuracy = SimpsonsRule(numfp, 0, aperture_size, prepass)/multiplyer;
    float integral_term = AdaptiveRichardsonExtrapolate(numfp, 0., aperture_size, accuracy);
    float constant = (sin(inclination) * sin(inclination))/disk_mass(aperture_size);
    return constant * integral_term;
}

// //////////////////////////
// // Black Hole Functions //
// //////////////////////////

void Galaxy::Construct_Black_Hole(float bhMass)
{
    BlackHolePresent = true;
    BHMass = bhMass;
}

// //////////////////////////
// //// Member Functions ////
// //////////////////////////

float GetVelocityDispersion(float Aperture,
                            float redshift,
                            float bulge_mass,
                            float bulge_radius,
                            float bulge_beta,
                            float bulge_sersicIndex,
                            int * componentFlag,
                            float disk_mass,
                            float Halo_mass,
                            char * profile_name,
                            char * c_path,
                            float BlackHole_mass)
{
    Galaxy aGalaxy(Aperture, redshift);
    aGalaxy.ConstructBulge(bulge_mass, bulge_beta, bulge_radius, bulge_sersicIndex);
    aGalaxy.ConstructDisk(disk_mass);
    aGalaxy.ConstructHalo(Halo_mass, profile_name, c_path);
    aGalaxy.Construct_Black_Hole(BlackHole_mass);

    aGalaxy.SetBulgeGravitationalContributions(componentFlag[0], componentFlag[1], componentFlag[2]);

    float bulge_sigma = aGalaxy.sigma_ap();
    float disk_sigma2 =

    assert(!isnan(res) && "GetVelocityDispersion() returned NaN");
    assert(!isinf(res) && "GetVelocityDispersion() returned inf");

    return res;

}



