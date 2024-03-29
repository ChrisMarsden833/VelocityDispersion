#include "bulge.h"

Bulge::Bulge(float input_stellar_property,
             float input_Re,
             float input_n,
             float input_beta,
             bool input_trace,
             bool input_grav,
             Halo * associatedHalo, 
             BlackHole * associatedBH,
             bool input_using_variable_IMF,
             float input_IMF_R0,
             float input_IMF_R1){
    // Constructor - so assign input values
    stellar_property = input_stellar_property;
    Re = input_Re;
    n = input_n;
    beta = input_beta;
    trace = input_trace;
    grav = input_grav;
    using_variable_IMF = input_using_variable_IMF;
    IMF_R0 = input_IMF_R0;
    IMF_R1 = input_IMF_R1;
    
    // Compute b_n
    b_n = 2.*n - (1./3.) + (.009876/n);

    // Compute I_e
    Ie =  stellar_property / (2. * PI * Re * Re * n * exp(b_n) * boost::math::tgamma(2.*n) / pow(b_n, (2.*n)));

    // Compute p_n
    p_n = 1. - .6097/n  + .00563/(n*n);

    // If beta is an exact multiple of 0.5 there is a singularity. We offset it by a small amount. 
    if( fmod(beta, 0.5) == 0.) beta += 1e-7;

    // Calculate the gamma term within the Kernel, often used and only depends on beta
    gamma_term = (float) boost::math::tgamma(beta - 0.5) / boost::math::tgamma(beta);

    memberHalo = associatedHalo;
    memberBH = associatedBH;

    // Precompute the extra mass from the variable IMF if necessary.
    if(using_variable_IMF) precompute_eXtra_density();
} 

float Bulge::alpha(float r){
    // Alpha, the gradient of the IMF
    float Y = r/Re;
    if(Y > 1) return 0.;
    return (IMF_R0 - IMF_R1)/IMF_R1;
}

float Bulge::ML(float r){
    // The actual projected mass-to-light ratio
    float Y = r/Re;
    float alpha = this->alpha(r);
    if(Y < 1) return IMF_R1;
    return IMF_R1 * (1. + alpha - alpha * Y);
}

float Bulge::sersic_profile(float r){
    // The sersic profile. Returns Mass, not light.
    float prefactor = Ie;
    if(using_variable_IMF) prefactor *= this->ML(r);
    return prefactor * exp(-b_n * ( pow(r/Re, 1./n) - 1 )); 
};

float Bulge::sersic_light_profile(float r){
    // The sersic profile. Returns light, not mass.
    float prefactor = Ie;
    return prefactor * exp(-b_n * ( pow(r/Re, 1./n) - 1 ));
}

float Bulge::prug_sim_profile(float r){
    // Prugniel and Simien profile (deprojected sersic).
    if(r == 0.0) r = 1e-7;
    float prefactor = Ie;
    if(using_variable_IMF) prefactor *= IMF_R1;
    float rho0 = prefactor * exp(b_n) * pow(b_n, n*(1-p_n)) * (boost::math::tgamma(2.*n)/(2.*Re* boost::math::tgamma(n * (3 - p_n))));
    float power_term = pow(r/Re, -p_n);
    float res = rho0 * power_term * exp(-b_n * pow(r/Re, 1/n));

    if(res == 0) return 1e-7;

    return res;
};

float Bulge::prug_sim_profile_light(float r){
    // Prugniel and Simien profile (deprojected sersic).
    if(r == 0.0) r = 1e-7;
    float prefactor = Ie;
    float rho0 = prefactor * exp(b_n) * pow(b_n, n*(1-p_n)) * (boost::math::tgamma(2.*n)/(2.*Re* boost::math::tgamma(n * (3 - p_n))));
    float power_term = pow(r/Re, -p_n);
    float res = rho0 * power_term * exp(-b_n * pow(r/Re, 1/n));
    return res;
};

float Bulge::rhoX_integrand(float Y){
    // The integrand used in the calculation of the 'extra' density when a variable IMF is considered.
    float alpha = this->alpha(Y*Re);

    if(Y == y) Y += 0.00001;

    float res = ( 1/(Y * sqrt(Y*Y - y*y)) ) * exp(-b_n * ( pow(Y, 1./n) - 1 )) \
         * ( (n * alpha * Y) + b_n * (alpha - alpha * Y) * pow(Y, 1./n));

    if(res < 0) return 0.0;
    return res;
}

float Bulge::rhoX_integrand_log(float log_Y){
    // The integrand used in the calculation of the 'extra' density when a variable IMF is considered.
    float Y = pow(10., log_Y);

    float alpha = this->alpha(Y*Re);

    if(abs(Y - y) < 0.001) Y += 0.001;

    float res = ( 1/(Y * sqrt(Y*Y - y*y)) ) * exp(-b_n * ( pow(Y, 1./n) - 1 )) \
         * ( (n * alpha * Y) + b_n * (alpha - alpha * Y) * pow(Y, 1./n)) * Y * log(10);

    if(res < 0) return 0.0;
    return res;
}

float Bulge::rhoX(float r){
    // The extra density when the variable imf is considered.
    y = r/Re;
    if(y < 1){
        //auto fp = bind(&Bulge::rhoX_integrand, this, _1); // An unfortunate consequence of using classes
        auto fp = bind(&Bulge::rhoX_integrand_log, this, _1); // An unfortunate consequence of using classes
        
        float accuracy  = 1e3 * (PI * Re * n) / ((IMF_R1 * Ie));
        
        //float integral = AdaptiveRichardsonExtrapolate(fp, log10(y), 0., accuracy);
        float integral = SimpsonsRule(fp, log10(y), 0, 1000);
        
        float res = ((IMF_R1 * Ie) / (PI * Re * n)) * integral;

        if(res < 0) {
            # pragma omp critical
            {
                assert((res >= 0.) || assert_msg("rhoX was about to return negative (" << res << ")." << std::endl <<
            "-----r = " << r << endl << 
            "-----y (r/Re) = " << y << endl <<
            "-----integral = " << integral << endl <<
            "-----Lum = " << stellar_property << " log10 " << log10(stellar_property) << endl << 
            "-----IMF R0 " << IMF_R0 << endl <<
            "-----IMF R1 " << IMF_R1 << endl));
            }
        }
        return res;
    }
    else return 0.0;
}

void Bulge::precompute_eXtra_density(void){
    // Precompute the extra density
    int len = 500; // Number of subdivisions. More will be more precise, but increase time.
    float eX, massX;
    
    // Bunch of variables used for various housekeeping things
    float lowest = 100;
    float lowest_mass = 1e50;
    float rolling_mass = 0.0;
    float previous_rolling_mass = rolling_mass;
    float previous_r_val = 0.0;
    
    // The grid is spaced in log10 r/re.
    float low = -7; // Log10
    float lowkpc = pow(10., low) * Re;
    r_domain = linspace(low, 1., len);

    //float step = (Re - low)/((float) len);
    
    // Iterate through the grid
    for(int i = 0.; i < len; i++){
        
        // Calculate the extra mass at this step
        eX = this->rhoX( pow(10., (*r_domain)[i]) * Re );

        if(eX < 0) {
            # pragma omp critical
            {
                assert((eX >= 0.) || assert_msg("Precomputed density was determiend to be negative (" << eX << ")." << std::endl <<
            "-----r = " << pow(10., (*r_domain)[i]) * Re << endl << 
            "-----Lum = " << stellar_property << " log10 " << log10(stellar_property) << endl << 
            "-----IMF R0 " << IMF_R0 << endl <<
            "-----IMF R1 " << IMF_R1 << endl));
            }
        }

        assert((!isnan(eX) && !isinf(eX)) || assert_msg("Precomputed density was about to return " << eX << std::endl <<
            "-----r = " << pow(10., (*r_domain)[i]) * Re << std::endl));


        if(isnan(eX)) eX = lowest;
        if(eX < lowest) lowest = eX;
        extra_density->push_back(eX);

        //rolling_mass = this->get_mass_long(pow(10., (*r_domain)[i]) * Re);
        
        // Compute the cumulative mass by 'hand'
        if(i==0){
            rolling_mass = (4./3.) * PI * lowkpc * lowkpc * lowkpc * (*extra_density)[i];
            previous_r_val = pow(10., (*r_domain)[i]) * Re;
        }
        else{
            float r_val = pow(10., (*r_domain)[i]) * Re;
            rolling_mass += 4. * PI * r_val * r_val * (*extra_density)[i] * (r_val - previous_r_val);
            previous_r_val = r_val;
        }
        
        // Checks
        if(rolling_mass <= 0.) {
            #pragma omp critical
            {
                assert((rolling_mass >= 0.) || assert_msg("rolling mass was determiend to be negative (" << rolling_mass << ")." << std::endl <<
                    "-----r = " << pow(10., (*r_domain)[i]) * Re << std::endl <<
                    "-----extra density = " << (*extra_density)[i] << std::endl <<
                    "-----current value = " << rolling_mass << std::endl <<
                    "-----previous value = " << previous_rolling_mass << std::endl ));            }
        }

        if(i > 0) {
            assert(rolling_mass >= (*extra_mass)[i-1] || assert_msg("rolling mass was decreasing " << rolling_mass  << " vs " << (*extra_mass)[i-1] << endl));
        }

        assert((!isnan(rolling_mass) && !isinf(rolling_mass)) || assert_msg("rolling mass was about to return " << eX << std::endl <<
            "-----r = " << pow(10., (*r_domain)[i]) * Re << std::endl));

        extra_mass->push_back(rolling_mass);
    }

}

float Bulge::interpolate_eXtra_density(float r){
    // Interpolate the extra density from the precomputed arrays, when using a variable IMF.
    float l10_r = log10(r/Re);
    if(r < Re & l10_r > r_domain->front()) return LinearInterp(r_domain, extra_density, l10_r);
    else if(l10_r <= r_domain->front()) return extra_density->front();
    return 0.0;
}

float Bulge::interpolate_eXtra_mass(float r){
    // Intepolate the extra mass, in the same way.
    if(!using_variable_IMF) return 0.0;
    // Get it into the same units
    float l10_r = log10(r/Re);

    if(r < Re & l10_r > r_domain->front()){
        float res = LinearInterp(r_domain, extra_mass, l10_r);
        // Checks        
        assert((res >= 0.) || assert_msg("interpolated eXtra mass was determiend to be negative (" << res << ")." << std::endl <<
                    "-----r = " << r << std::endl <<
                    "-----Re = " << Re << std::endl ));
        assert((!isnan(res) && !isinf(res)) || assert_msg("interpolated eXtra mass was about to return " << res << std::endl <<
            "-----r = " << r << std::endl));
        return res;
    } 
    else if(l10_r <= r_domain->front()) return 0.0;
    return extra_mass->back();
}

float Bulge::get_mass_long_integrand(float r){
    // Standard mass integrand, for the extra mass when there is a variable IMF
    return 4. * PI * r * r * this->rhoX(r); 
}

float Bulge::get_mass_long(float r){
    // Full mass calculation. Very slow, so currently unused.
    auto fp = bind(&Bulge::get_mass_long_integrand, this, _1);

    float low = 0.001;

    if(r <= low) low = r/2.;

    //float multiplyer = 1e3;
    //int prepass = 12;
    // Do a prepass to get a sense of the value of the integral
    //float accuracy = SimpsonsRule(fp, 0, r, prepass)/multiplyer;
    //float res = AdaptiveRichardsonExtrapolate(fp, 0.0, r, 1e5);
    float res = SimpsonsRule(fp, low, r, 500.);
    return res;
}

float Bulge::get_density(float r){
    // Function to get the overall density
    float res = this->prug_sim_profile(r);
    
    if(using_variable_IMF) res += this->interpolate_eXtra_density(r);
    
    assert((res != 0.) || assert_msg("Get density about to return 0" << res << endl <<
        "-----r = " << r << endl <<
        "-----Re = " << Re << endl <<
        "-----Prug_sim(r) = " << this->prug_sim_profile(r) << endl));
    return res;
}

float Bulge::analytic_mass_within(float r){
    float threemp_term = (3. - p_n) * n;
    float z = b_n * pow(r/Re, 1/n);
    float prefactor = Ie;
    if(using_variable_IMF) prefactor *= IMF_R1;
    float rho0 = prefactor * exp(b_n) * pow(b_n, n*(1-p_n)) * (boost::math::tgamma(2.*n)/(2.*Re* boost::math::tgamma(n * (3 - p_n))));
    return 4.*PI*n*Re*Re*Re * rho0 * pow(b_n, -(3.-p_n)*n) * boost::math::tgamma_lower(threemp_term, z);
}

float Bulge::projected_mass_within_integrand(float R){
    return 2. * PI * R * this->sersic_profile(R);
}

float Bulge::projected_mass_within(float R){
    auto fp = bind(&Bulge::projected_mass_within_integrand, this, _1); // An unfortunate consequence of using classes
    float res = AdaptiveRichardsonExtrapolate(fp, 0.0, R, 10.);
    return res;
}

float Bulge::mass_within(float r){
    if(!grav) return 0.0;

    float res = this->analytic_mass_within(r); 
    if(using_variable_IMF) res += this->interpolate_eXtra_mass(r); // interpolate extra mass

    assert((!isnan(res) && !isinf(res)) || assert_msg("mass_within (bulge) about to return " << res << std::endl <<
        "-----r = " << r << std::endl <<
        "-----Re = " << Re << std::endl <<
        "-----Analytic mass within = " << this->analytic_mass_within(r)  << std::endl <<
        "-----Mass within variable IMF = " << this->interpolate_eXtra_mass(r) << std::endl));

    return res;
}

float Bulge::get_beta(void)
{
    return beta;
}

float Bulge::K_Kernel(float u)
{
    // Kernel K(U), equation 9.
    if(u == 0.) return 0.;
    if(beta > 0.5) beta = 0.5;

    float prefactor = 0.5 * pow(u, (2. * beta - 1.));
    float term1 = ((3./2.) - beta) * pow(PI, 0.5) * gamma_term;
    float term2 = beta * incompleteBeta(beta + 0.5, 0.5, 1. / (u * u));
    float term3 = -incompleteBeta(beta - 0.5, 0.5, 1. / (u * u));
    float res = prefactor * (term1 + term2 + term3);
    if((res < 0.) && (res < 1e-7)) res = 1e-7;
    assert((res >= 0.) || assert_msg("K_Kernel about to return negative (" << res << ")." << std::endl <<
       "-----beta = " << beta << std::endl <<
       "-----u = " << u << std::endl <<
       "-----gamma term = " << gamma_term << std::endl <<
       "res = prefactor * (term1 + term2 + term3)" << std::endl <<
       "-----prefactor: ((3./2.) - beta) * pow(PI, 0.5) * gamma_term = " << prefactor << std::endl <<
       "-----term 1: ((3./2.) - beta) * pow(PI, 0.5) * gamma_term = " << term1 << std::endl <<
       "-----term 2: beta * incompleteBeta(beta + 0.5, 0.5, 1./(u*u)) = " << term2 << std::endl <<
       "-----term 3: -incompleteBeta(beta - 0.5, 0.5, 1./(u*u)) = " << term3));
    return res;
}

float Bulge::TotalMass(float r) {
    float res = 0.0;
    if(grav) res += this->mass_within(r);
    if(memberHalo) res += memberHalo->mass_within(r);
    if(memberBH) res += memberBH->getMass();
    

    assert((!isnan(res) && !isinf(res)) || assert_msg("TotalMass (bulge) about to return " << res << std::endl <<
        "-----r = " << r << std::endl <<
        "-----Re = " << Re << std::endl <<
        "-----Mass within bulge = " << this->mass_within(r) << std::endl <<
        "-----Mass within Halo = " << memberHalo->mass_within(r) << std::endl <<
        "-----Mass within BH  = " << memberBH->getMass() << std::endl));

    return res;
}

float Bulge::sigma_integrand(float log_r)
{
    float r = pow(10., log_r);
    // Manage R, as this can be problematic if R = 0
    if(R == 0) R = 1e-7;
    float ratio = r/R;

    assert(ratio >= 1.0  || assert_msg("sigma_integrand ratio r/R < 0 (value is " << ratio <<
      "). This will violate the incomplete beta function." << std::endl <<
      "-----r = " << r << std::endl <<
      "-----R = " << R << std::endl <<
      "-----HLR = " << Re));

    float Kernel = this->K_Kernel(ratio);
    float density = this->get_density(r);
    float Mass = this->TotalMass(r);

    float res;
    if(r == 0) r = 1e-7;

    res = Kernel * density * Mass * r * log(10) /  r;

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

float Bulge::sigma_LOS(float R_arg)
{

    float integral_term = 0., accuracy = 0., upper_limit = 0., los_accuracy = 0.;
    if(R_arg == 0.) R_arg = 1e-8;

    float bulge_density = this->sersic_profile(R_arg);
    if(bulge_density < 1e-8) bulge_density = 1e-8;

    float K = (2. * GR) / bulge_density;

	R = R_arg;
	auto fp = bind(&Bulge::sigma_integrand, this, _1); // An unfortunate consequence of using classes

    float approx = sqrt( GR * this->TotalMass(R) / (R_arg * R_arg * 2.5 ) );
    
    /*
    accuracy = 0.1; //min( approx * 0.1, 1.0);
    
    los_accuracy = pow(accuracy, 2)/K; // Transform this into the units of the integral

    // Let's find a sensible upper limit for the integral.
    float step = Re;
    upper_limit = R_arg + step; // We subdivide into
    float unit_accuracy = 2. * los_accuracy/(upper_limit - R_arg); // Unit accuracy is the accuracy divided by the interval
    while(unit_accuracy > los_accuracy/(upper_limit - R_arg))
    {
        unit_accuracy = this->sigma_integrand(upper_limit) / (upper_limit - R_arg) ;
        upper_limit += step;
    }*/
    

    //integral_term = AdaptiveRichardsonExtrapolate(fp, log10(R_arg), log10(upper_limit), los_accuracy);
    integral_term = SimpsonsRule(fp, log10(R_arg), 5, 100);

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
                 "-------Bulge Mass/Luminosity " << stellar_property << std::endl <<
                 "-------Bulge Density @ r " << bulge_density << std::endl <<
                 "-------Bulge Re " << Re << std::endl <<
                 "-------Integration was performed between " << R_arg << " and " << upper_limit << std::endl <<
                 "-------Accuracy = " << los_accuracy << std::endl <<
                 "-------Full Integration (ARE) = " << integral_term << std::endl));
	    }
	}
	return value;
}

float Bulge::Jeans_innerIntegrand_1(float r){
    if(r - R == 0.0) r += 0.0001;
    return (1./this->Jeans_f(r)) * r * this->prug_sim_profile_light(r) / (this->get_density(r) * sqrt(r*r - R*R));
}

float Bulge::Jeans_outerIntegrand_1(float s){
    auto fp = bind(&Bulge::Jeans_innerIntegrand_1, this, _1); // An unfortunate consequence of using classes
    if(R - s == 0.0) s += 0.0001;
    float intergral = SimpsonsRule(fp, R, s, 1000);
    return intergral * this->Jeans_f(s) * this->get_density(s) * this->TotalMass(s)/ (s*s);
}

float Bulge::Jeans_innerIntegrand_2(float r){
    if(r - R == 0.0) r += 0.0001;

    float res =  this->Jeans_fprime(r) * this->prug_sim_profile_light(r)/(this->get_density(r) * this->Jeans_f(r) * this->Jeans_f(r) * sqrt(r*r - R*R));

    assert( (!isnan(res) && !isinf(res)) || assert_msg(std::endl << "Jeans_innerIntegrand_2 about to return " << res << endl <<
                                        "--- this->Jeans_f(r) :" << this->Jeans_f(r) << endl <<
                                        "--- this->get_density(r) :" << this->get_density(r) << endl <<
                                        "--- Re :" << Re << endl <<
                                        "--- Ie : " << Ie << endl <<
                                        "--- r : " << r << endl <<
                                        "--- R : " << R << endl));
    return res;
}

float Bulge::Jeans_outerIntegrand_2(float s){
    auto fp = bind(&Bulge::Jeans_innerIntegrand_2, this, _1); // An unfortunate consequence of using classes
    if(R - s == 0.0) s += 0.0001; 
    float intergral = SimpsonsRule(fp, R, s, 1000);

    assert( (!isnan(intergral) && !isinf(intergral)) || assert_msg(std::endl << "Jeans_outerIntegrand_2 about to return " << intergral << endl <<
                                        "--- s : " << s << endl <<
                                        "--- R : " << R << endl));

    return intergral  * this->Jeans_f(s) * this->get_density(s) * this->TotalMass(s)/(s*s);
}

float Bulge::ChaeIntegrand(float log_u){
    float u = pow(10., log_u);
    float arg = R * sqrt(1 + u * u);
    float res = this->prug_sim_profile_light(arg) * this->sigma_radial2(arg) * (1 + beta * arg / ( 1 + u * u)) * u * log(10);
    return res;
}

float Bulge::sigma_LOS_full(float R_arg){

    //assert((R_arg < 20.*Re) || assert_msg(endl << "R_arg > 10.Re" << endl));

    R = R_arg;
    auto fp = bind(&Bulge::ChaeIntegrand, this, _1); // An unfortunate consequence of using classes

    float l10 = log10(R_arg);
    if(R_arg == 0.0) l10 = -7.;

    float integral = SimpsonsRule(fp, -7, 5, 100);
    float res  = sqrt(integral * 2. * R / this->sersic_light_profile(R_arg));

    assert( (!isnan(res) && !isinf(res)) ||
                    assert_msg(std::endl << "Sigma LOS full full about to return " << res << endl <<
                                        "--- Integral 1: " << integral << endl));
    return res;
}

float Bulge::sigma_radial_integrand(float r){
    float res = GR * this->Jeans_f(r) * this->get_density(r) * this->TotalMass(r) / (r * r);
    return res;
}

float Bulge::sigma_radial_integrand_log(float log_r){
    float r = pow(10., log_r);
    float res = GR * this->Jeans_f(r) * this->get_density(r) * this->TotalMass(r) * r * log(10) / (r * r)  ;
    return res;
}

float Bulge::sigma_radial2(float r){
    auto fp = bind(&Bulge::sigma_radial_integrand_log, this, _1);

    float l10 = log10(r);

    assert( (!isnan(l10)) & (!isinf(l10)) || assert_msg(std::endl << "sigma_radial2 l10 = " << l10 << ", r = " << r << endl));

    float upper = 5.0;
    if(l10 >= upper) upper = l10 + 0.1;
  
    float accuracy = 1.0; //kms^-1
    accuracy *= this->Jeans_f(r) * this->get_density(r);

    float intergral = SimpsonsRule(fp, l10, upper, 50);
    //float intergral = AdaptiveRichardsonExtrapolate(fp, r, pow(10., upper), accuracy); //SimpsonsRule(fp, r, up, 1000);

    assert( (!isnan(intergral) && !isinf(intergral)) ||
                    assert_msg(std::endl << "sigma radial2 intergral " << intergral << endl <<
                                        "--- r: " << r << endl <<
                                        "--- R: " << R << endl << 
                                        "--- l10: " << l10 << endl << 
                                        "--- upper: " << upper << endl)); 
 
    float density = this->get_density(r);
    if(density < 1e-7) density = 1e-7;
    
    float res =  1/(this->Jeans_f(r) * density) * intergral;

    assert( (!isnan(res) && !isinf(res)) ||
                    assert_msg(std::endl << "sigma radial2 about to return " << res << endl <<
                                        "--- r: " << r << endl <<
                                        "--- R: " << R << endl << 
                                        "--- this->Jeans_f(r): " << this->Jeans_f(r) << endl << 
                                        "--- density (potenitally perturbed): " << density << endl <<
                                        "--- this->get_density(r): " << this->get_density(r) << endl)); 

    return res;
}

float Bulge::alt_integrand1(float r){
    if( abs(R  - r) < 1e-7) r += 1e-7;
    float res = this->prug_sim_profile_light(r) * this->sigma_radial2(r) * r / (sqrt(r*r - R*R));
    assert( (!isnan(res) && !isinf(res)) ||
                    assert_msg(std::endl << "alt_integrand1 about to return " << res << endl <<
                                        "--- r: " << r << endl <<
                                        "--- R: " << R << endl <<
                                        "--- this->sigma_radial2(r)" << this->sigma_radial2(r) << endl << 
                                        "--- this->prug_sim_profile_light(r)" << this->prug_sim_profile_light(r) << endl));
    return res;
}

float Bulge::alt_integrand1_log(float log_r){
    float r = pow(10., log_r);

    
    float fac = 1.000001;
    if(r < fac * R ) r = fac * R;
    
    assert(!isnan(r) || assert_msg(endl << "alt_integrand1_log r = nan. Log10r = " << log_r << endl));

    float internal = r*r - R*R;
    
    float res = this->prug_sim_profile_light(r) * this->sigma_radial2(r) * r * r * log(10)/ (sqrt(internal));
    assert( (!isnan(res) && !isinf(res)) ||
                    assert_msg(std::endl << "alt_integrand1 log about to return " << res << endl <<
                                        "--- r: " << r << endl <<
                                        "--- R: " << R << endl <<
                                        "--- r*r - R*R  = " << r*r - R*R << endl <<
                                        "--- internal sqrt(r*r - R*R) = " << internal << endl <<
                                        "--- this->sigma_radial2(r)" << this->sigma_radial2(r) << endl << 
                                        "--- this->prug_sim_profile_light(r)" << this->prug_sim_profile_light(r) << endl));
    return res;
}

float Bulge::alt_integrand2(float r){
    if(R == r) r += 0.001;
    float res = beta * this->prug_sim_profile_light(r) * this->sigma_radial2(r) / (r * sqrt(r*r - R*R));
    assert( (!isnan(res) && !isinf(res)) ||
                    assert_msg(std::endl << "alt_integrand2 about to return " << res << endl <<
                                        "--- r: " << r << endl <<
                                        "--- R: " << R << endl <<
                                        "--- this->sigma_radial2(r)" << this->sigma_radial2(r) << endl << 
                                        "--- this->prug_sim_profile_light(r)" << this->prug_sim_profile_light(r) << endl));
    return res;
}

float Bulge::alt_integrand2_log(float log_r){
    float r = pow(10., log_r);
    
    float fac = 1.000001;
    if(r < fac * R ) r = fac * R;

    assert(!isnan(r) || assert_msg(endl << "alt_integrand2 log r = nan. Log10r = " << log_r << endl));

    float internal = r*r - R*R;
    
    float res = beta * this->prug_sim_profile_light(r) * this->sigma_radial2(r) * r * log(10) / (r * sqrt(internal));
    assert( (!isnan(res) && !isinf(res)) ||
                    assert_msg(std::endl << "alt_integrand2 log about to return " << res << endl <<
                                        "--- r: " << r << endl <<
                                        "--- R: " << R << endl <<
                                        "--- Denominator: " << (r * sqrt(internal)) << endl <<
                                        "--- this->sigma_radial2(r)" << this->sigma_radial2(r) << endl << 
                                        "--- this->prug_sim_profile_light(r)" << this->prug_sim_profile_light(r) << endl));
    return res;
}

float Bulge::Jeans_f(float r){
    return pow(r, 2. * beta);
}

float Bulge::Jeans_fprime(float r){
    return 2. * beta *pow(r, 2. * beta - 1);
}

float Bulge::sigma_ap_integrand(float R_arg){
    // Integrand of sigma aperture
    float sigma;
    if(using_variable_IMF) sigma = this->sigma_LOS_full(R_arg);
    else sigma = this->sigma_LOS(R_arg);

	float MDP = this->sersic_light_profile(R_arg);
	return MDP * sigma * sigma * R_arg;
}

float Bulge::sigma_ap_integrand_log(float log_R_arg){
    float R_arg = pow(10., log_R_arg);

    // Integrand of sigma aperture
    float sigma;
    if(using_variable_IMF) sigma = this->sigma_LOS_full(R_arg);
    else sigma = this->sigma_LOS(R_arg);

	float MDP = this->sersic_light_profile(R_arg);
	return MDP * sigma * sigma * R_arg * R_arg * log(10);
}


float Bulge::sersicxR(float r){
    // Just a wrapper function, sometimes we need it *r for integrals.
    float res = this->sersic_light_profile(r) * r;
    return res;
}

float Bulge::sersicxR_log(float log_r){
    // Just a wrapper function, sometimes we need it *r for integrals.
    float r = pow(10., log_r);
    float res = this->sersic_light_profile(r) * r * r * log(10);
    return res;
}

float Bulge::sersicLightWithin_integrand_log(float log_r){
    // Just a wrapper function, sometimes we need it *r for integrals.
    float r = pow(10., log_r);
    float res = this->sersic_light_profile(r) * r * r * 2 * PI * log(10);
    return res;
}

float Bulge::LightWithin(float R){
    // Get light within for weigting of integral.
    auto numfp = bind(&Bulge::sersicLightWithin_integrand_log, this, _1);
    float res = SimpsonsRule(numfp, -7, log10(R), 50.);
    return res;
}


float Bulge::sigma_ap(float aperture)
{
    if(!trace) return 0.0;

    float numerator, denominator, accuracy;

    // Full value of halo aperture
	auto numfp = bind(&Bulge::sigma_ap_integrand_log, this, _1);
    auto denfp = bind(&Bulge::sersicxR_log, this, _1);

    float multiplyer = 1e5;
    int prepass = 12;
    // Do a prepass to get a sense of the value of the integral
    //accuracy = SimpsonsRule(numfp, -5, log10(aperture), prepass)/multiplyer;
    //numerator =  AdaptiveRichardsonExtrapolate(numfp, -5, log10(aperture), accuracy);
    
    numerator = SimpsonsRule(numfp, -7, log10(aperture), 50.);

    //accuracy = SimpsonsRule(denfp, 0, aperture, prepass)/multiplyer;
    //denominator = AdaptiveRichardsonExtrapolate(denfp, -5, log10aperture, accuracy);
    denominator = SimpsonsRule(denfp, -7, log10(aperture), 50.);

    if(denominator == 0.) return 0.;

    float res = sqrt(numerator/denominator);

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
                                         "-------Aperture = " << aperture << " kpc" << std::endl <<
					 "This occured when:" << std::endl <<
					 "-------Bulge Mass/Luminsoity " << stellar_property << std::endl <<
			 	         "-------Bulge Re " << Re << std::endl ));
        }
    }

    return res;


}