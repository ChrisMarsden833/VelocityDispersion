#include "halo.h"

Halo::Halo(float input_rho, float input_rs, const char * input_string, bool input_grav){

    assert((!isnan(input_rho) && !isinf(input_rs)) || assert_msg("Halo, input rho = " << input_rho << " input rs = " << input_rs  << std::endl <<
        "-----rho = " << rho << std::endl <<
        "-----rs = " << rs  << std::endl));


    rho = input_rho;
    rs = input_rs;
    name = input_string;
    grav = input_grav;
}

float Halo::mass_within(float r){

    if(!grav) return 0.0;
    
    float res;

    if(strcmp(name, "NFW") == 0){
        float x = r/rs;
        res = 4. * PI * pow(rs, 3.) * rho * (log(1. + x) - x/(1+x));
        if(x/(1+x) > log(1. + x)) res = 0.0; // At very low radii, can return negative.
    }
    else if(strcmp(name, "Burkert") == 0){
        float x = r/rs;
        float tanterm = -2.*atan(x);
        float term2 =  2.*log(1+x);
        float term3 = log(1 + x*x);
        res = PI * rho * pow(rs, 3.) * (tanterm + term2 + term3);
        if(res < 0) res = 0.0;
    }
    else{
        assert(false >= 0 || assert_msg(" Unknown halo profile type: " << name << std::endl));
        res = 0.0;
    }


    assert((!isnan(res) && !isinf(res)) || assert_msg("mass_within (halo) about to return " << res << std::endl <<
        "-----r = " << r << std::endl <<
        "-----rho = " << rho << std::endl <<
        "-----rs = " << rs  << std::endl));

    return res;
}

float Halo::Vcirc2(float r){
    if(!grav) return 0.0;
    return GR * this->mass_within(r)/r;
}
