#ifndef VELOCITYDISPERSIONS_DESMOND_H
#define VELOCITYDISPERSIONS_DESMOND_H

#include "math.h"
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "integration.h"
#include "dark_matter.h"
#include <iostream>
#include <stdlib.h>

#define PI 3.14159265
#define GR 4.3009125e-6 // In units of kpc M_sun^-1 (km/s)^2
#define Prepass_Subdivisions 6

// Formulae from Desmond and Wechsler 2017

/** ++ MassDensityProfile ++
 * Equation (1). P824
 * @param r : float, the radius
 * @param SersicIndex : float, the Sersic Index
 * @param Half_Light_radius : float, half light radius.
 * @param stellar_mass : float, stellar mass
 * @return
 */
float MassDensityProfile(float r, float SersicIndex, float Half_Light_radius, float stellar_mass);

/** ++ MassDensityProfile_wrapper ++
 * Equation (1). P824. For integrating MassDensityProfile x r, WRT r.
 * @param r : float, the radius
 * @param args : vector of floats, representing the Sersic Index and the Half light Radius.
 * @return The value of the function at r
 */
float MassDensityProfile_wrapper(float r, std::vector<float> args);

/** ++ b_n ++
 * Equation (2), P824, parameter for mass density profile.
 * @param SersicIndex : float, The Sersic Index.
 * @return float, the value of b_n.
 */
float b_n(float SersicIndex);

/** ++ rho ++
 * Equation (3), P824, The de-projected volume density.
 * @param r : float, the radius
 * @param Half_Light_radius : float, the half light radius
 * @param SersicIndex : float, the Sersic Index
 * @param stellar_mass : float, the stellar mass
 * @param dm_rs : float, the dark matter scale radius
 * @param dm_c : float, the dark matter concentration parameter.
 * @return float, the value of the density.
 */
float rho(float r, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c);

/** ++ rho_0 ++
 * Equation (4), P824, Factor for the de-projected volume density.
 * @param Half_Light_radius : float, the half light radius.
 * @param SersicIndex : float, the Sersic Index
 * @param stellar_mass : float, the stellar mass
 * @return float, the value of rho_0
 */
float rho_0(float Half_Light_radius, float SersicIndex, float stellar_mass);


/** ++ p_n ++
 * Equation (5), P824, Factor for de-projected volume density.
 * @param Half_Light_radius : float, the half light radius
 * @param SersicIndex : float, the Sersic Index
 * @return float, the value of p_n
 */
float p_n(float SersicIndex);


/** ++ cumSpherRho ++
 * Guts of the cumulative spherical distribution, used internally for it's integration.
 * @param R : float, Radius
 * @param args : vector<float> of arguments for function, specifically {Half_Light_radius, float SersicIndex}
 * @return float, the value
 */
float cumSpherRho(float R, std::vector<float> args);

/** ++ cumSpherMassDistro ++
 * The Cumulative Spherical Mass distribution
 * @param R : float, the radius. (kpc)
 * @param Half_Light_radius : float, the half light radius. (kpc)
 * @param SersicIndex : float, the Sersic Index (dimensionless)
 * @param stellar_mass : float, the stellar mass (M_sun)
 * @param dm_rs : float, the scale radius r_s of the halo, for a NFW profile (kpc)
 * @param dm_c : float, the NFW halo concentration parameter, c (dimensionless)
 * @param omega_m : float, mass_density cosmological parameter (dimensionless)
 * @return float, the value of the cumulative spherical mass distribution (M_sun)
 */
float cumSpherMassDistro(float R, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c);

/** ++ K_Kernel_DW ++
 * Equation (9), P824, The Kernel function 'u'.
 * @param u : float, the argument for the function (dimensionless).
 * @param beta: float, the anisotropy profile (dimensionless).
 * @return float, the value of the Kernel.
 */
float K_Kernel_DW(float u, float beta);

/** ++ full_sigma_integral_internals ++
 * The internals from the integral part of equation (8), P284.
 * @param r : float, the radius over which we are integrating.
 * @param R : float, the Radius for the larger function.
 * @param beta : float, the anisotropy profile.
 * @param Half_Light_radius : float, the Half_Light_radius.
 * @param SersicIndex : float, the Sersic Index.
 * @param stellar_mass : float, the stellar mass.
 * @param dm_rho0 : float, the dark matter characteristic mass
 * @param dm_rs : float, the dark matter characteristic radius
 * @return float, the returned value.
 */
float full_sigma_integral_internals(float r, float R, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c);

/** ++ sigma_internals_wrapper ++
 * Wrapper function to make sure formatting works for sigma integration
 * @param r : float, the radius
 * @param args : vector of floats containing the other parameters (R, beta, Half_Light_radius and SersicIndex)
 * @return the value of the function an r.
 */
float sigma_internals_wrapper(float r, std::vector<float> args);

/** ++ sigma_internals_wrapper_transformed ++
 * Wrapper function to make sure formatting works for sigma integration, but with transformed coords for integration to +inf
 * @param t : float, the transformed argument
 * @param args : vector of floats containing the other parameters (R, beta, Half_Light_radius and SersicIndex)
 * @return the value of the function an r.
 */
//float sigma_internals_wrapper_transformed(float t, std::vector<float> args);

/** ++ function to calculate the LOS velocity dispersion
 * This is equation (8), P824.
 * @param R : float, the radius at which velocity disperion is calculated
 * @param beta : float, the anisotropy parameter
 * @param Half_Light_radius : float, the half light radius of the galaxy
 * @param SersicIndex : float, the sersic index of the galaxy.
 * @param stellar_mass : float, the stellar mass.
 * @param dm_rho0 : float, the dark matter characteristic mass
 * @param dm_rs : float the dark matter characteristic radius
 * @return float, the velocity dispersion.
 */
float sigma_los(float R, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c);

/** ++ sigma_los_wrapper ++
 * Wrapper function to ensure formatting works for integration.
 * @param R : float, the radius to be integrated over.
 * @param args : vector of floats containing the other parameters (beta, Half_Light_radius and SersicIndex).
 * @return float, the value of the function at R.
 */
float sigma_los_wrapper(float R, std::vector<float> args);

/** ++ sigma_apature_internals ++
 * Wrapper function for the top integral in equation 10.
 * @param r : float, the radius
 * @param args : vector of floats, the arguments
 * @return the value of the function at R.
 */
float sigma_apature_internals(float r, std::vector<float> args);

/** ++ sigma_aperture ++
 * This is equation (10), P824
 * @param R_ap : float, the size of the aperture
 * @param beta : float, the anisotropy profile.
 * @param Half_Light_radius : float, the half light radius.
 * @param SersicIndex : float, the sersic index.
 * @param stellar_mass : float, the stellar mass.
 * @return float, the value of sigma.
 */
float sigma_aperture(float R_ap, float beta, float Half_Light_radius, float SersicIndex, float stellar_mass, float dm_rs, float dm_c);




#endif //VELOCITYDISPERSIONS_DESMOND_H
