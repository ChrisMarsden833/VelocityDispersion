#ifndef BULGE_H
#define BULGE_H

#include "math.h"
#include <boost/math/special_functions/gamma.hpp>
#include "integration.h"
#include <stdlib.h>
#include <functional>
#include "halo.h"
#include "blackhole.h"
#include "omp.h"

#define PI 3.14159265
#define GR 4.3009125e-6 // In units of kpc M_sun^-1 (km/s)^2

using namespace std;
using namespace std::placeholders;

/* The bulge class is undoubtedly the most complex */
class Bulge
{
    public:
        // Constructor 
        Bulge(float input_stellar_property,
              float input_Re,
              float input_n,
              float input_beta = 0.0,
              bool input_trace = true,
              bool input_grav = true,
              Halo * associatedHalo = nullptr, 
              BlackHole * associatedBH = nullptr,
              bool input_using_variable_IMF = false,
              float input_IMF_R0 = 1.,
              float input_IMF_R1 = 1.);

        // Compute alpha, the slope of the ML ratio
        float alpha(float r);
        // Compute the ML ratio
        float ML(float r);

        // The projected sersic profile
        float sersic_profile(float r);
        // The sersic light profile, used for variable IMF
        float sersic_light_profile(float r);
        // The deprojected prugniel and simien profile 
        float prug_sim_profile(float r);
        float prug_sim_profile_light(float r);

        // The integrand for rhoX
        float rhoX_integrand(float Y);
        float rhoX_integrand_log(float log_Y);
        // The eXtra mass from the variable IMF
        float rhoX(float r);

        // Function to precompute the extra density. Massively increases speed!
        void precompute_eXtra_density(void);
        // Interpolation of the eXtra density and mass respectively.
        float interpolate_eXtra_density(float r);
        float interpolate_eXtra_mass(float r);

        // Ingetrand and Function to get the mass 'fully', used in testing and precompution.
        float get_mass_long_integrand(float r);
        float get_mass_long(float r);
        // Generic get_density and get_mass, which have their own controls to detect a variable IMF
        float get_density(float r);
        float mass_within(float r);
        // Mass within according to an analytic integration of Prugniel+Simien profile
        float analytic_mass_within(float r);

        // Integrant and function for calculating the projected mass within. Used in testing.
        float projected_mass_within_integrand(float R);
        float projected_mass_within(float R);

        // Function for computing beta.
        float get_beta(void);
        // The Kernel, K(u)
        float K_Kernel(float u);

        // Total mass at a specific r
        float TotalMass(float r);

        // Integrand and function for LOS sigma.
        float sigma_integrand(float r);
        float sigma_LOS(float R_arg);

        // Integrands for the full Jeans equations solution
        float Jeans_innerIntegrand_1(float r);
        float Jeans_outerIntegrand_1(float s);
        float Jeans_innerIntegrand_2(float r);
        float Jeans_outerIntegrand_2(float s);
        float Jeans_f(float r);
        float Jeans_fprime(float r);
        float sigma_LOS_full(float R_arg);
        float sigma_radial_integrand(float r);
        float sigma_radial_integrand_log(float log_r);
        float sigma_radial2(float r);
        float alt_integrand1(float r);
        float alt_integrand1_log(float log_r);
        float alt_integrand2(float r);
        float alt_integrand2_log(float log_r);

        // (pair of) integrands for aperture sigma
        float sigma_ap_integrand(float R_arg);
        float sigma_ap_integrand_log(float log_R_arg);
        float sersicxR(float r);
        float sersicxR_log(float log_r);
        float sigma_ap(float aperture);

        // Chae's formula for integrating sigma
        float ChaeIntegrand(float log_u);

        float sersicLightWithin_integrand_log(float log_r);
        float LightWithin(float R);

    private:
        float stellar_property; // Stellar mass [M_sun] OR Stellar Luminsoity [L_sun], in the case of a variable IMF
        
        bool using_variable_IMF; // If variable IMF is on or not.
        float IMF_R0; // The IMF parameters
        float IMF_R1;
        float y; // Rolling value used in integration
        float gamma_term; // Gamma term used in sersic profiles.

        bool trace; // If we should 'trace' the bulge
        bool grav; // If we should consider the gravitational influence of the bulge
        
        float Re; // Projected half light radius Re [kpc]
        float n; // Sersic index [dimensionless]
        float beta; // The anisotropy

        float b_n; // 2.*n - (1./3.) + (.009876/n)
        float Ie; // The constant "at the front" of the sersic profile.
        float p_n; // p_n in the PS etc profiles.

        float R = 1.0; // Rolling Constant used in integration.

        int count = 0; // Counter. Not sure this is used. 

        float domain_size;

        // Arrays for the precomputed domains.
        vector<float> * r_domain; 
        vector<float> * extra_density = new vector<float>;
        vector<float> * extra_mass = new vector<float>;

        // Pointers to the halo and black holes that are a part of this galaxy.
        Halo * memberHalo = nullptr;
        BlackHole * memberBH = nullptr;
    
};

#endif