#ifndef VELOCITYDISPERSIONS_DESMOND_H
#define VELOCITYDISPERSIONS_DESMOND_H

#include "math.h"
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "integration.h"

#define PI 3.14159265

// Formulae from Desmond and Wechsler 2017

/** ++ MassDensityProfile ++
 * Equation (1). P824
 * @param r : float, the radius
 * @param SersicIndex : float, the Sersic Index
 * @param Half_Light_radius : float, half light radius.
 * @return
 */
float MassDensityProfile(float r, float SersicIndex, float Half_Light_radius);

/** ++ b_n ++
 * Equation (2), P824, parameter for mass density profile.
 * @param SersicIndex : float, The Sersic Index.
 * @return float, the value of b_n.
 */
float b_n(float SersicIndex);

/** ++ rho ++
 * Equation (3), P284, The de-projected volume density.
 * @param r : float, the radius
 * @param Half_Light_radius : float, the half light radius
 * @param SersicIndex : float, the Sersic Index
 * @return float, the value of the density.
 */
float rho(float r, float Half_Light_radius, float SersicIndex);

/** ++ rho_0 ++
 * Equation (4), P284, Factor for the de-projected volume density.
 * @param Half_Light_radius : float, the half light radius.
 * @param SersicIndex : float, the Sersic Index
 * @return float, the value of rho_0
 */
float rho_0(float Half_Light_radius, float SersicIndex);


/** ++ p_n ++
 * Equation (5), P284, Factor for de-projected volume density.
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
 * @param R : float, the radius.
 * @param Half_Light_radius : float, the half light radius.
 * @param SersicIndex : float, the Sersic Index
 * @return float, the value of the cumulative spherical mass distribution.
 */
float cumSpherMassDistro(float R, float Half_Light_radius, float SersicIndex);

/** ++ K_Kernel_DW ++
 * Equation (9), P284, The Kernel function 'u'.
 * @param u : float, the argument for the function
 */
float K_Kernel_DW(float u);





#endif //VELOCITYDISPERSIONS_DESMOND_H
