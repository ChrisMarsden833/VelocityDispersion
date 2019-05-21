#ifndef VELOCITYDISPERSIONS_UTILLITY_H
#define VELOCITYDISPERSIONS_UTILLITY_H

#include <iostream>
#include <stdio.h>
#include "stdlib.h"
#include "math.h"
#include "vector"

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

#endif //VELOCITYDISPERSIONS_UTILLITY_H
