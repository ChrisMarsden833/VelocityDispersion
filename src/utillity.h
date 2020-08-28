#ifndef VELOCITYDISPERSIONS_UTILLITY_H
#define VELOCITYDISPERSIONS_UTILLITY_H

#include <iostream>
#include <stdio.h>
#include "stdlib.h"
#include "math.h"
#include "vector"
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include "omp.h"
#include <cassert>

#define PI 3.14159265
#define EU 2.71828182
#define assert_msg(x) !(std::cerr << std::endl << \
    "#######################################################" << std::endl << \
    "################# Assertion failed ####################" << std::endl << \
    x << std::endl << \
    "#######################################################" << std::endl)

/** ++ factorial ++
 * Simple recursive function to return the factorial of an integer
 * @param n : template, value
 * @return template, the factorial of value.
 */
template <typename T> T factorial(T n);

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

/** ++ linspace ++
 * Simple function to mimic numpy's linspace function. This creates
 * a linearly spaced vector of values.
 * @param a : the value at the start of the array
 * @param b : the value at the end of the array
 * @param N : the number of values in the array.
 * @return
 */
template<typename T>
std::vector<T> linspace(T start_in, T end_in, int num = 100)
{
    std::vector<T> linspaced;

    T start = static_cast<T>(start_in);
    T end = static_cast<T>(end_in);
    //T num = static_cast<T>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    T delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspaced;
}

/** ++ AreSame ++
 * A simple function to test that float like variables are the same
 * @param a, first value
 * @param b, second value
 * @return bool, the same or not.
 */
template<typename T>
bool AreSame(T a, T b)
{
    return fabs(a - b) <= fabs(a) * 0.0001;
}

/** ++ CheckAndMaybeIncrementN ++
 * A function that simply checks N is even, and if it's not it will increment it by one.
 * @param N, pointer to integer N.
 * @param fname (optional), the name of the function calling it, for error purposes.
 */
void CheckAndMaybeIncrementN(int * N, std::string fname = "Richardson Extrapolation");


void fastLinspace(float * &grid, float &h, float a, float b, int N = 100);

/** ++ ReadFile ++
 * CSV reader,
 */
std::vector<std::vector<float>> * ReadFile(std::string path, std::vector<int> * IndexesToGrab);

/* ++ InputPrechecks
 * Function to test input indexes
 */
void InputPrechecks(std::vector<int> * IndexesToGrab);

/* ++ checkOutpuConsistency ++ 
 * Function to check output consistency
 */
void checkOutputConsistency(std::vector<int> * IndexesToGrab, std::vector<std::vector<float>> * OutputArrays);

/* ++ Reduce
 * Reduce to single values
 */
std::vector<float> * Reduce(std::vector <float> * input, std::vector<float> * output);

/* ++ FindClosest
 * Find the closest value from a list
 */
float FindClosest(float value, std::vector<float> * data);

/* ++ Equals
 * Generate a mask where the array equals the value 
 */
void Equals(std::vector<float> * array, float value, std::vector<bool> * mask);

/* ++ MaskOut
 * Mask out values in a vector based on a mask vector
 */
void MaskOut(std::vector<float> * array, std::vector<bool> * mask);

/* ++ LinearInterp
 * Linear Interpolation
 */
float LinearInterp(std::vector<float> * X, std::vector<float> * Y, float x);

#endif //VELOCITYDISPERSIONS_UTILLITY_H
