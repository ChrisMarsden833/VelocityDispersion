#ifndef VELOCITYDISPERSION_DARK_MATTER_H
#define VELOCITYDISPERSION_DARK_MATTER_H

#include "math.h"

#define PI 3.14159265
#define GR 4.3009125e-6

/** ++ delta_vir ++
 * Function defining delta virial, the average overdensity.
 * @param omega_m : float, mass density cosmological parameter (dimensionless)
 * @return float, the value of the average overdensity (dimensionless)
 */
float delta_vir(float omega_m);

/** ++ critical_density ++
 * function to return the critical density in terms of H.
 * @param H, Hubble's Constant, in (km s^-1 MPc^-1)
 * @return float, the value of the critical density (M_sun kpc^-3)
 */
float critical_density(float H);

/** ++ delta_char ++
 * delta_char, the numerator of the NFW profile
 * @param omega_m : float, mass_density cosmological parameter (dimensionless)
 * @param c : float, concentration parameter (dimensionless)
 * @return float, the value of delta_char (dimensionless)
 */
float delta_char(float omega_m, float c);

/** ++ NFW_profile ++
 * Navarro–Frenk–White (NFW) profile is a spatial mass distribution of dark matter fitted to dark matter halos.
 * @param r : float, radius at which to return the density (kpc)
 * @param rs : float, halo scale radius (kpc)
 * @param c : float, the concentration parameter, r_vir/rs (dimensionless)
 * @param omega_m : float, mass density cosmological parameter (dimensionless)
 * @param H : float, Hubble's Constant ((km/s)/Mpc)
 * @return float, the value of the density (M_sun/kpc
 */
float NFW_profile(float r, float rs, float c, float omega_m, float H);

#endif //VELOCITYDISPERSION_DARK_MATTER_H
