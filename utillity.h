#ifndef VELOCITYDISPERSIONS_UTILLITY_H
#define VELOCITYDISPERSIONS_UTILLITY_H

#include <iostream>
#include <stdio.h>
#include "stdlib.h"
#include "math.h"
#include "vector"
#include <boost/math/special_functions/beta.hpp>

/** ++ factorial ++
 * Simple recursive function to return the factorial of an integer
 * @param n : integer value
 * @return integer, the factorial of value.
 */
int factorial(int n);

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


#endif //VELOCITYDISPERSIONS_UTILLITY_H
