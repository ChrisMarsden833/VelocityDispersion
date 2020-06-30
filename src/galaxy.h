#ifndef GALAXY_H
#define GALAXY_H

#include <stdlib.h>
#include <string>
#include "math.h"
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <functional>
#include "integration.h"
#include <random>


#define PI 3.14159265
#define GR 4.3009125e-6 // In units of kpc M_sun^-1 (km/s)^2
#define h 0.69
#define Om 0.3
#define zero_perturbation 0.00000001
#define precision 3


using namespace std;
using namespace std::placeholders;

class Galaxy
{
	public:
		Galaxy(float input_stellar_mass,
		       float input_beta,
		       float input_half_light_radius,
		       float input_aperture_size,
		       float input_sersic_index);

		void init_dark_matter(string input_profile_name,
				      float input_concentration);

		// Function to return the mass density at radius r - Equation (1). P824
		float MassDensity(float r);

		// Function returning the mass density * r (for integration)
		float MassDensityr(float r);

		// Function to return the de-projected density at radius r.
		float rho(float r);

		// Function to return the mass within infintesimal shell at radius R
		float mass_shell(float R);

		// Function to to return the cumulative mass at radius r.
		float cumulative_mass(float R_arg);

		// Function to return the value of the K_Kernal
		float K_Kernel_DW(float u);

		// The Integrand for sigma (internals of the integral)
		float sigma_integrand(float r);

		// The value of sigma in the LOS
		float sigma_los(float R_arg);

		// The integrand of sigma aperture (numerator)
		float sigma_ap_integrand(float R_arg);

		// The velocity dispersion within the aperture
		float sigma_ap(void);

		// The Halo concentration
		void GetHaloC(bool scatter);

		// The Halo Radius
		void GetHaloR(void);

	private:
		// Stellar Mass of the galaxy [log10 M_sun]
		float stellar_mass;
		// Redshift of the galaxy
		float redshift;
		// Beta (anisotropy parameter) [dimensionless]
		float beta;
		// Half Light Radius [kpc]
		float half_light_radius;
		// Effective Half Light Radius perturbed to prevent /0 errors
		float half_light_radius_eff; 
		// Aperture size [kpc]
		float aperture_size;
		// Sersic Index [dimensionless]
		float sersic_index;
		// Effective Sersic Index perturbed to prevent /0 errors 
		float sersic_index_eff;

		// R, value to be used as part of the integral
		float R;
		// Switch if dark matter is on or not.
		bool dark_matter_on;
		// Dark Matter Profile Name
		string profile_name;
		// Halo Mass;
		float HaloMass; // Log10
		// Dark Matter Concentration parameter
		float concentration;
		// Halo Radius
		float HaloRadius; // Kpc


		// b_n - parameter required for the density profile
		float b_n;
		// sigma_e - parameter required for mass density profile
		float sigma_e;
		// p_n - parameter needed for density
		float p_n;

		// rho0 - constant required for de-projected mass density.
		float rho0;

};

float GetVelocityDispersion(float input_aperture_size,
			    float input_beta,
			    float input_half_light_radius,
		   	    float input_sersic_index,
			    float input_stellar_mass);

#endif