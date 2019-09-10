#ifndef VELOCITYDISPERSIONS_INTEGRATION_H
#define VELOCITYDISPERSIONS_INTEGRATION_H

// Includes

#include <iostream>
#include <stdio.h>
#include "stdlib.h"
#include "math.h"
#include "vector"
#include "utillity.h"

// General info
/**
 * This file/library uses Richardson extrapolation, with methodology based on Dr Ian Hawke's numerical methods lecture
 * series at the University of Southampton, specifically the section on numerical quadrature.
 */

// Structs
/** ++ RResult ++
 * Wrapper Struct that contains the integration value as well as it's accuracy
 * Contains float values of:
 *  Integral - the value of the integration
 *  Accuracy - the accuracy to which the integral has been calculated.
 * Mainly used for Richardson extrapolation (see below).
 */
struct RResult
{
    float integral;
    float accuracy;
};

// Functions

/** ++ SimpsonsRule ++
 * Performs the Simpsons rule for integration on the specified function.
 * @param f : pointer, to the function in question. The first argument of this function must be the variable over which
           the integration will be performed, the second must accept an array containing any additional arguments. If
           this is a problem, place function in wrapper function and pass that.
 * @param a : float, the lower limit of the integration.
 * @param b : float, the upper limit of the integration.
 * @param N : int, the number of intervals to perform the integration over, or the 'resolution' of the pass.
 * @param extraArguments : vector<float>, containing any extra arguments for f that might be required.
 * @return float, the value of the integration.
 */
float SimpsonsRule(float (*f)(float, std::vector<float>), float a, float b, int N, std::vector<float> extraArguments);

/** ++ RichardsonExtrapolate ++
 * Performs Richardson Extrapolation for the specified function using Simpsons rule.
 * @param f : pointer, to the function in question. The first argument of this function must be the variable over which
           the integration will be performed, the second must accept an array containing any additional arguments. If
           this is a problem, place function in wrapper function and pass that.
 * @param a : float, the lower limit of the integration.
 * @param b : float, the upper limit of the integration.
 * @param steps2 : int, the number of intevals that will be performed on the 'higher' richardson integration, meaning
            that it is a requirement that this number is even.
 * @param extraArguments : vector<float>, containing any extra arguments for f that might be required.
 * @return RResult, the value of the integration coupled with the accuracy of the integration.
 */
RResult RichardsonExtrapolate(float (*f)(float, std::vector<float>), float a, float b, int steps2, std::vector<float> extraArguments);

/** ++ AdaptiveRichardson Extrapolate ++
 * Performs Adaptive Richardson Extrapolation for the specified function using Simpsons rule
 * @param f : pointer, to the function in question. The first argument of this function must be the variable over which
           the integration will be performed, the second must accept an array containing any additional arguments. If
           this is a problem, place function in wrapper function and pass that.
 * @param a : float, the lower limit of the integration.
 * @param b : float, the upper limit of the integration.
 * @param accuracy : float, the required tolerance of the integration. Setting this loop too low can cause this function
            to loop endlesslly. TODO: If this becomes a serious problem, introduce self detection and errors for this.
 * @param extraArguments : vector<float>, containing any extra arguments for f that might be required.
 * @return float, the value of the integration to the specified tolerance.
 */
float AdaptiveRichardsonExtrapolate(float (*f)(float, std::vector<float>), float a, float b, float accuracy, std::vector<float> extraArguments);

/** ++ IterativeRichardsonExtrapolation ++
 * Performs Iterative Richardson Extrapolation for the specified function using Simpsons rule. This is a slightly more
 * 'brute force' approach to Richardson Extrapolation, just increasing the spacing until the error reaches it's desired
 * value.
 * @param f : pointer, to the function in question. The first argument of this function must be the variable over which
           the integration will be performed, the second must accept an array containing any additional arguments. If
           this is a problem, place function in wrapper function and pass that.
 * @param a : float, the lower limit of the integration.
 * @param b : float, the upper limit of the integration.
 * @param accuracy : float, the required tolerance of the integration.
 * @param extraArguments : vector<float>, containing any extra arguments for f that might be required.
 * @return float, the value of the integration to the specified tolerance.
 */
float IterativeRichardsonExtrapolate(float (*f)(float, std::vector<float>), float a, float b, float accuracy, std::vector<float> extraArguments);

#endif //VELOCITYDISPERSIONS_INTEGRATION_H
