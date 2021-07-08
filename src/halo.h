#ifndef HALO_H
#define HALO_H

#include "math.h"
#include <boost/math/special_functions/gamma.hpp>
#include "integration.h"
#include <stdlib.h>
#include <cassert>
#include <functional>

#define PI 3.14159265
#define GR 4.3009125e-6 // In units of kpc M_sun^-1 (km/s)^2

#define assert_msg(x) !(std::cerr << std::endl << \
    "#######################################################" << std::endl << \
    "################# Assertion failed ####################" << std::endl << \
    x << std::endl << \
    "#######################################################" << std::endl)

using namespace std;
using namespace std::placeholders;


class Halo
{
    public:
        Halo(float input_rho, float input_rs, const char * input_string, bool input_grav);
        float mass_within(float r);
        float Vcirc2(float r);

    private:
        float rs;
        float rho;
        const char * name;
        bool grav;

};

#endif