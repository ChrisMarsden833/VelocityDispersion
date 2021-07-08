#ifndef SIGMALIB_H
#define SIGMALIB_H

#include <stdlib.h>
#include <string>
#include "math.h"
#include "bulge.h"
#include "disk.h"
#include "blackhole.h"

/* The main function for calculating sigma */
float Sigma(float aperture,
            float bulge_mass,
            float bulge_radius,
            float bulge_sersic_index,
            float bulge_beta,
            float disk_mass,
            float disk_h,
            float disk_i,
            const char * halo_profile = "None",
            float haloRs = 0.0,
            float haloRhos = 0.0,
            float black_hole_mass = 0.0,
            int * tracer_flags = nullptr,
            int * gravitational_flags = nullptr,
            bool use_variable_IMF = false,
            float varIMF_R1 = 0.0,
            float varIMF_R0 = 0.0,
            int mode = 1);

#endif