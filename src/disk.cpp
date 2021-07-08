#include "disk.h"

Disk::Disk(float input_mass, float input_h, float input_i, bool input_trace, bool input_grav,
            Bulge * input_bulgeptr, Halo * input_haloptr){
    mass = input_mass;
    h = input_h;
    i = input_i;
    trace = input_trace;
    grav = input_grav;
}

float Disk::projected_density(float R){
    float mfactor = (1 - exp(-R/h)*(1+R/h));
    float sigma_0 = mass/(2. * PI * h * h) * mfactor;
    return sigma_0 * exp(-R/h);
}

float Disk::bessel(float x){
    return boost::math::cyl_bessel_i(0., x) * boost::math::cyl_bessel_k(0., x) -
                     boost::math::cyl_bessel_i(1., x) * boost::math::cyl_bessel_k(1., x);
}

float Disk::Vcirc2(float R){
    if(R == 0.) R += 1e-7;
    float x = R/(2.*h);

    float B_term;
    if(x > 80) B_term = 1.0; // Insane situations
    else B_term = this->bessel(x); //disk_bessel(1.6);
    float res = 0.0;
    if(grav) res += 2. * (GR * mass)/h * x * x * B_term;
    if(AssociatedHalo) res += AssociatedHalo->Vcirc2(R);
    if(AssociatedBulge) res += GR* AssociatedBulge->mass_within(R)/R;
    return res;
}

float Disk::Mass(float R){
    float mfactor = (1 - exp(-R/h)*(1+R/h));
    float sigma_0 = mass/(2. * PI * h * h) * mfactor;
    float res = 2 * PI * sigma_0 * h * h * mfactor;
    if(res < 0. || isnan(res) || isinf(res) ) res = 0.;
    return res;
}

float Disk::disk_integrand(float R){
    return 2. * PI * R * this->projected_density(R) * this->Vcirc2(R);
}

float Disk::VelocityDispersion(float aperture){

    auto numfp = bind(&Disk::disk_integrand, this, _1);

    float min = 0.;
    float multiplyer = 1e6;
    int prepass = 10;
    float accuracy = SimpsonsRule(numfp, min, aperture, prepass)/multiplyer;
    float integral_term = AdaptiveRichardsonExtrapolate(numfp, min, aperture, accuracy);
    float mass_within = this->Mass(aperture);
    if(mass_within == 0.0) return 0.0;
    float constant = (sin(i) * sin(i))/mass_within;
    float res = constant * integral_term;

    assert( (!isnan(res) && !isinf(res)) ||
        assert_msg("disk_velocity_dispersion2 about to return: " << res << std::endl <<
                         "---- Constant " << constant << std::endl <<
                         "---- Integral Term " << integral_term << std::endl <<
                         "---- aperture_size : " << aperture << std::endl <<
                         "---- mass disk : " << mass << std::endl <<
                         "---- Disk_inclination: " << i << std::endl));
    return res;
}

