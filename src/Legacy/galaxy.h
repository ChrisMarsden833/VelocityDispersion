#ifndef GALAXY_H
#define GALAXY_H

#include <stdlib.h>
#include <string>
#include "math.h"
#include <functional>
#include "integration.h"
#include <cassert>
#include "bulge.h"
#include "halo.h"
#include "disk.h"

#define PI 3.14159265
#define GR 4.3009125e-6 // In units of kpc M_sun^-1 (km/s)^2
#define zero_perturbation 1e-8
#define assert_msg(x) !(std::cerr << std::endl << \
    "#######################################################" << std::endl << \
    "################# Assertion failed ####################" << std::endl << \
    x << std::endl << \
    "#######################################################" << std::endl)

using namespace std;
using namespace std::placeholders;

class Galaxy
{
    public:
        Galaxy(bool input_use_variable_IMF, float magnitude, float logsigma);
        void get_IMF_params(float magnitude, float logsigma);
        void ConstructBulge(float bulge_mag, float bulge_beta, float bulge_re, float bulge_n);
        void ConstructHalo(float input_rho, float input_rs, const char * input_string);
        void ConstructDisk(float disk_mass, float disk_h, float disk_i);

    private:

        bool use_variable_IMF;
        float R0, R8, R1;
        Bulge * memberBulge = nullptr;
        Halo * memberHalo = nullptr;
        Disk * memberDisk = nullptr;

        //
};

#endif