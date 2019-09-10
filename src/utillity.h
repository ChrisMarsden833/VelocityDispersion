#ifndef VELOCITYDISPERSIONS_UTILLITY_H
#define VELOCITYDISPERSIONS_UTILLITY_H

#include <iostream>
#include <stdio.h>
#include "stdlib.h"
#include "math.h"
#include "vector"
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/factorials.hpp>

#define PI 3.14159265
#define EU 2.71828182

/** ++ factorial ++
 * Simple recursive function to return the factorial of an integer
 * @param n : template, value
 * @return template, the factorial of value.
 */
template <typename T> T factorial(T n);

/** ++ LargeFactorial ++
 * A more meaty function for safe calculation of large factorials
 * @param n : int, value
 * @return template, the factorial value in log10.
 */
int largeFactoraial(int n);

/** ++ IsEven ++
 * Internal function to check a value is even - error will be thrown otherwise.
 * @param Value : int, Value that will be checked
 */
bool checkEven(int Value);

/** ++ incompleteBeta ++
 * Wrapper function to parameterize the incomplete beta
 * @param a : float, the value of a
 * @param b : float, the value of b
 * @param z : float, the value of z
 * @return float, the value of the function.
 */
float incompleteBeta(float a, float b, float z);

/** ++ pochhammer ++
 * Numerical implementation of the pochhammer symbol
 * @param q : float, the value of q
 * @param n : int, the value of n
 * @return float, the value of the function.
 */
double pochhammer(double q, int n);

/** ++ hyperGeometricSeries ++
 * Numerical implementation of the HyperGeometric Series
 * @param a : float, the value of a
 * @param b : float, the value of b
 * @param c : float, the value of c
 * @param z : float, the value of z
 * @return float, the value of the function
 */
float hyperGeometricSeries(float a, float b, float c, float z);

/** ++ HypergeometricSeries_2 ++
 * Numerical implementation of the HyperGeometric Series 2
 * @param a : float, the value of a
 * @param b : float, the value of b
 * @param c : float, the value of c
 * @param z : float, the value of z
 * @return float, the value of the function
 */
float hyperGeometricSeries_2(float a, float b, float c, float z);

/** ++ getSign ++
 * Simple function to get the sign of a float.
 * @param argument : float, the argument
 * @return : float, 1 or -1
 */
float getSign(float argument);



#endif //VELOCITYDISPERSIONS_UTILLITY_H