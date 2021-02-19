#ifndef GALAXY_H
#define GALAXY_H

#include <stdlib.h>
#include <string>
#include "math.h"
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <functional>
#include "integration.h"
#include <random>
#include "omp.h"
#include <cassert>
#include <algorithm>

#define PI 3.14159265
#define GR 4.3009125e-6 // In units of kpc M_sun^-1 (km/s)^2
#define h 0.69
#define Om 0.3
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
        // Constructor
		Galaxy(float input_aperture_size);

        // -----------------------
        // --- Bulge functions ---
        // -----------------------

        void ConstructBulge(float input_bulge_mass, float input_bulge_beta, float input_bulge_half_light_radius, float input_sersic_index);

		// Function to return the mass density at radius r - Equation (1). P824
		float BulgeProjectedDensity(float r);

		// Function returning the bulge mass density * r (for integration)
		float BulgeProjectedDensityxR(float r);

		// Function to return the de-projected density at radius r.
		float BulgeDensity(float r);

		// Function to return the analytic expression of the mass at radius r.
        float BulgeMass(float r);

        // Function to return the total (contributing) mass within the bulge ar radius r [kpc]
        float BulgeTotalMass(float r);

        // set the respective contributions to bulge mass
        void SetBulgeGravitationalContributions(bool stars, bool dark_matter, bool black_hole);

		// Function to return the value of the K_Kernal
		float K_Kernel_DW(float u);

		// The Integrand for sigma (internals of the integral)
		float bulge_sigma_integrand(float r);

		// The value of sigma in the LOS
		float bulge_sigma_los(float R_arg);

		// The integrand of sigma aperture (numerator)
		float sigma_ap_integrand(float R_arg);

		// The velocity dispersion within the aperture
		float sigma_ap(void);

		// Set bulge remenant prefactor
		void set_rem_prefac(float input_prefactor);

        // +++++++++++++++++++++++++
        // ++++ Halo Functions +++++
        // +++++++++++++++++++++++++

        // Set up Dark Matter in the galaxy.
        void ConstructHalo(std::string input_profile_name, float haloRs, float haloRhos);

	// Halo analytic mass
	float HaloAnalyticMass(float r);


	// Halo Circular Velocity
	float HaloVcirc2(float r);

        // \\\\\\\\\\\\\\\\\\\\\\\\
        // \\\\ Disk Functions \\\\
        // \\\\\\\\\\\\\\\\\\\\\\\\

        // Set up the galactic disk
        void ConstructDisk(float DiskMass, float input_inclination, float disk_scale);

        // The disk surface density
        float disk_projected_density(float R);

        // The disk mass as a function of r [log10 solar masses]
        float disk_mass(float R);

        // Disk bessel function used in vcirc
        float disk_bessel(float x);

        // The disk circular velocity
        float disk_Vcirc2(float r);

        // disk velocity dispersion integral
        float disk_integrand(float R);

        // Disk velocity dispersion
        float disk_velocity_dispersion2();

        // ]]]]]]]]]]]]]]]]]]]]]]]]]]]]
        // ]]] Black Hole Functions ]]]
        // ]]]]]]]]]]]]]]]]]]]]]]]]]]]]

        // Setter for black hole mass
        void Construct_Black_Hole(float bhMass);

	private:
        // ==========================
        // === General Properties ===
        // ==========================

        // Aperture size [kpc]
        float aperture_size;

        // Take longer integrating, but (potenitally) better accuracy
        bool slow_integrate = false;

        // ++++++++++++++++++++++++++
        // ++++ Halo Properties +++++
        // ++++++++++++++++++++++++++

        bool halo_present = false;

        // Dark Matter Profile Name
        string profile_name;

        // Halo Radius
        float HaloRadius; // [Kpc] R_vir
	    // HaloDensity
	    float HaloRhos; // M_sun/kpc^3

        // --------------------------
        // ---- Bulge Properties ----
        // --------------------------

        // If a bulge exists or not.
        bool bulge_present = false;
		// Stellar Mass of the bulge [log10 M_sun]
		float bulge_stellar_mass = 0.;
        // Beta, The disk anisotropy parameter [dimensionless]
        float bulge_beta = 0.;
        // (bulge) Half Light Radius [kpc]
        float bulge_half_light_radius = 0.;

        // -- Contribution control component influence on bulge --
        // Switch if baryonic Matter is on or not
        bool bulge_stars_gravitation_on = true;
        // Switch if dark matter is on or not.
        bool dark_matter_gravitation_on = true;
        // Switch on black hole or not.
        bool black_hole_gravitation_on = true;

		// Sersic Index [dimensionless]
		float bulge_sersic_index = 0.;
		// The bulge mass density at z.
		float Sigma0 = 0.;
		// Term involving gamma functions, it's best to only calculate once
        float gamma_term = 0.0;
		// R, value to be used as part of the integral
		float R = 0.0;
        // b_n - parameter required for the density profile
        float b_n = 0.;
        // sigma_e - parameter required for mass density profile
        float sigma_e = 0.;
        // p_n - parameter needed for density
        float p_n = 0.;
        // rho0 - constant required for de-projected mass density.
        float rho0 = 0.;
        // Constant for the mass, best saved.
        float mass_prefactor = 0.;
        // Reminant Prefactor  - 1.72
        float rem_prefac  = 1.; // Default is 1, which corresponds to no change.

        // \\\\\\\\\\\\\\\\\\\\\\\\\\\\
        // \\\\ Disk Properties \\\\\\\
        // \\\\\\\\\\\\\\\\\\\\\\\\\\\\

        bool disk_present = false;
        // The mass of the disk [log10 M_sun]
        float mass_disk = 0.;
        // Disk Scale length [kpc]
        float disk_scale_length = 0.;
        // Disk Inclination
        float disk_inclination = 0.;

        // ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
        // ]]] Black Hole Properties ]]]
        // ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

        // Is there a black hole
        bool BlackHolePresent = false;
        // Black Hole Mass
        float BHMass = 0.;

};

float GetVelocityDispersion(float Aperture,
                            float bulge_mass,
                            float bulge_radius,
                            float bulge_beta,
                            float bulge_sersicIndex,
                            float rem_prefactor,
                            int * componentFlag,
                            float disk_mass,
                            float disk_inclination,
                            float disk_scale_length,
                            char * profile_name,
                            float haloRs,
                            float haloRhos,
                            float BlackHole_mass,
                            int mode);

#endif
