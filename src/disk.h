#ifndef DISK_H
#define DISK_H

#include "math.h"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include "integration.h"
#include <stdlib.h>
#include <cassert>
#include <functional>
#include "bulge.h"
#include "halo.h"

#define PI 3.14159265
#define GR 4.3009125e-6 // In units of kpc M_sun^-1 (km/s)^2

#define assert_msg(x) !(std::cerr << std::endl << \
    "#######################################################" << std::endl << \
    "################# Assertion failed ####################" << std::endl << \
    x << std::endl << \
    "#######################################################" << std::endl)

using namespace std;
using namespace std::placeholders;

class Disk
{
    public:
        Disk(float input_mass, float input_h, float input_i, bool input_trace, bool input_grav,
            Bulge * input_bulgeptr = nullptr, Halo * input_haloptr = nullptr);
        float projected_density(float R);
        float bessel(float x);
        float Vcirc2(float R);
        float Mass(float R);
        float disk_integrand(float R);
        float VelocityDispersion(float aperture);
        float projected_densityL(float R);
        float LightProfile_integrand_log(float log_r);
        float LightWithin(float r);
    
    private:
        float mass;
        float h;
        float i;
        Bulge * AssociatedBulge = nullptr;
        Halo * AssociatedHalo = nullptr;
        bool trace;
        bool grav;
        float Lum;
    
};

#endif