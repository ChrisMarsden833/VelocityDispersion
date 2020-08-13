#include "integration.h"

// Functions

float SimpsonsRule(float (*f)(float, std::vector<float>), float a, float b, int N, std::vector<float> extraArguments)
{
    // Test and adjustment for evenness.
    //std::cout << "Simpsons rule overloaded N: " << N << std::endl;

    CheckAndMaybeIncrementN(&N, "Simpsons Rule Overloaded");

    float * grid = NULL;
    float h;
    fastLinspace(grid, h, a, b, N);

    int N_half = N/2;
    // Term 1 in the formula.
    float sum1 = 0;
    for(int j = 1; j <= (N_half - 1); j++)
    {
        sum1 += (*f)(grid[2*j], extraArguments);
    }
    // Term 2 in the formula
    float sum2 = 0;
    for(int j = 1; j <= N_half; j++)
    {
        sum2 += (*f)(grid[2*j-1], extraArguments);
    }
    // Full formula
    float total = (h/3.0)*((*f)(a, extraArguments) + 2.0 * sum1 + 4 * sum2 + (*f)(b, extraArguments));
    free(grid);
    return total;
}

float SimpsonsRule(float (*f)(float), float a, float b, int N)
{
    CheckAndMaybeIncrementN(&N, "Simpsons Rule, basic");

    float * grid = NULL;
    float h;
    fastLinspace(grid, h, a, b, N);

    // Start the actual Simpsons rule
    int N_half = N/2;

    // Term 1 in the formula.
    float sum1 = 0;
    for(int j = 1; j <= (N_half - 1); j++) sum1 += (*f)(grid[2*j]);

    // Term 2 in the formula
    float sum2 = 0;
    for(int j = 1; j <= N_half; j++) sum2 += (*f)(grid[2*j-1]);

    // Full formula
    float total = (h/3.0)*((*f)(a) + 2.0 * sum1 + 4 * sum2 + (*f)(b));
    free(grid);
    return total;
}

float SimpsonsRule(std::function<float (float)> fun, float a, float b, int N)
{
    CheckAndMaybeIncrementN(&N, "Simpsons Rule, basic");

    assert(b > a || assert_msg("Simpsons rule does not allow b <= a " << std::endl <<
        "----a = " << a << std::endl <<
        "----b = " << b << std::endl));

    float * grid = NULL;
    float h;
    fastLinspace(grid, h, a, b, N);

    // Start the actual Simpsons rule
    int N_half = N/2;

    // Term 1 in the formula.
    float sum1 = 0;
    for(int j = 1; j <= (N_half - 1); j++) sum1 += fun(grid[2*j]);

    // Term 2 in the formula
    float sum2 = 0;
    for(int j = 1; j <= N_half; j++) sum2 += fun(grid[2*j-1]);

    // Full formula
    float total = (h/3.0)*(fun(a) + 2.0 * sum1 + 4 * sum2 + fun(b));

    if(total < 0)
    {
        #pragma omp critical
        {
            std::cout << "######### Simpson's rule full report ##########" << std::endl;
            std::cout << "-Grid: [";
            for(int i = 0; i < N; i++) std::cout << grid[i] << ", ";
            std::cout << "]" << std::endl;

            std::cout << "-Term 1 (sum of below elements):" << std::endl;
            for(int j = 1; j <= (N_half - 1); j++)  std::cout << "-fun(grid[2*" << j << "]) = " << fun(grid[2*j]) << std::endl;
            std::cout << "-Term 2 (sum of below elements):" << std::endl;
            for(int j = 1; j <= N_half; j++) std::cout << "-fun(grid[2*" << j << "-1]) = " << fun(grid[2*j-1]) << std::endl;

            std::cout << "fun(a) [a = " << a << "] = " << fun(a) << std::endl;
            std::cout << "fun(b) [b = " << b << "]= " << fun(b) << std::endl;

            assert( (total >=0 ) || assert_msg(std::endl << "Simpsons Rule about to return negative." << std::endl <<
                "---Value = " << total << std::endl));
        }
    }

    free(grid);
    return total;
}



RResult RichardsonExtrapolate(float (*f)(float, std::vector<float>), float a, float b, int steps2, std::vector<float> extraArguments)
{
    // Struct for the return
    RResult data;

    // Call Simpson's rule twice with different intervals accoring to richardson extrapolation.
    float I2n = SimpsonsRule((*f), a, b, 2*steps2, extraArguments);
    float In = SimpsonsRule((*f), a, b, steps2, extraArguments);

    // Assign the value of the integral according to the formulae for the integral and it's accuracy
    data.integral = (pow(2, 4) * I2n - In)/(pow(2, 4) - 1);
    data.accuracy = (abs(In - I2n)/(pow(2, 4) - 1));
    return data;
}

RResult RichardsonExtrapolate(std::function<float (float)> fun, float a, float b, int steps2)
{
    // Struct for the return
    RResult data;

    // Call Simpsons rule twice at different intervals according to richarson extrapolation.
    float I2n = SimpsonsRule(fun, a, b, 2*steps2);
    float In = SimpsonsRule(fun, a, b, steps2);

    // Assign the value of the integral according to the formulae for the integral and it's accuracy.
    data.integral = (pow(2., 4.) * I2n - In)/(pow(2., 4.) - 1);

    if(data.integral <= 0)
    {
        #pragma omp critical
        {
            assert((data.integral > 0) || assert_msg(std::endl << "Richardson Extrapolate about to return a negative integral" << std::endl <<
                "----Value = " << data.integral << std::endl <<
                "----I2n = " << I2n << std::endl <<
                "----In = " << In << std::endl <<
                "---Domain, a = " << a << ", b = " << b << ". Over (2x) " << steps2 << " steps." << std::endl));
        }
    }

    data.accuracy = abs(In - I2n)/(pow(2., 4.) - 1 );
    return data; 
}

float AdaptiveRichardsonExtrapolate(float (*f)(float, std::vector<float>), float a, float b, float accuracy, std::vector<float> extraArguments)
{
    // Preliminary vector for the boundaries.
    std::vector<float> boundaries = {a, b};

    // Maximum (sane) number of subdivisions).
    //int max_subdivisions = 1000;
    float minimum_acceptable_subdivision_size = abs(b-a)/max_subdivisions;

    // Struct and values for iteration.
    RResult data;
    float integral = 0;
    int calls = 0;

    float new_element;
    float acceptable_local_error;

    // Loop over the boundaries vector elements. Because elements are dynamically added, there remains the possibility
    // of an infinite loop if very low accuracy is used or an unreasonable function is called.
    for(int i = 0; i < boundaries.size()-1; i++)
    {
        bool loop = true; // Stopper
        while(loop)
        {
            // Call Richardson Extrapolation
            data = RichardsonExtrapolate((*f), boundaries[i], boundaries[i+1], 2, extraArguments);
            // Calculate the acceptable local error for this region.
            acceptable_local_error = accuracy*abs(boundaries[i+1] - boundaries[i])/abs(b-a);
            // Based on this, decide if we need to continue.
            loop = (data.accuracy > acceptable_local_error);
            // Work out the size of the new element
            new_element = (boundaries[i] + boundaries[i+1])/2;

            if (new_element - boundaries[i] < minimum_acceptable_subdivision_size and loop)
            {
                // Get the hell out of dodge
                loop = false;

                /*
                try
                {
                    throw std::invalid_argument(" ++ Woah, go easy on the accuracy there; maximum subdivision size reached ++ ");
                }
                catch (const std::invalid_argument& ia) {
                    std::cerr << ia.what() << "\n" << "Integral will move on, but maximum acceptable (local) accuracy was defined as "
                    << acceptable_local_error << " whereas the accuracy of this element is only " << data.accuracy
                    << ". This may or not be a problem for you." << std::endl;
                }
                 */
                // We are loosing some accuracy here in these extreme cases. TODO: Parameterize this, for the user.
            }

            if(loop)
            {
                // Add a new boundaries element exactly between the boundaries used above



                boundaries.insert(boundaries.begin() + i + 1, new_element);
            }
            else
            {
                // Append the integral value to the integral variable.
                integral += data.integral;
            }
        }
    }
    return integral;
}

float AdaptiveRichardsonExtrapolate(std::function<float (float)> fun, float a, float b, float accuracy)
{
    // Preliminary vector for the boundaries.
    std::vector<float> boundaries = {a, b};

    // Maximum (sane) number of subdivisions).
    //int max_subdivisions = 1000;
    float minimum_acceptable_subdivision_size = abs(b-a)/max_subdivisions;

    // Struct and values for iteration.
    RResult data;
    float integral = 0;
    int calls = 0;

    float new_element;
    float acceptable_local_error;

    // Loop over the boundaries vector elements. Because elements are dynamically added, there remains the possibility
    // of an infinite loop if very low accuracy is used or an unreasonable function is called.
    for(int i = 0; i < boundaries.size()-1; i++)
    {
        bool loop = true; // Stopper
        while(loop)
        {
            // Call Richardson Extrapolation
            data = RichardsonExtrapolate(fun, boundaries[i], boundaries[i+1], 2);
            // Calculate the acceptable local error for this region.
            acceptable_local_error = accuracy*abs(boundaries[i+1] - boundaries[i])/abs(b-a);
            // Based on this, decide if we need to continue.
            loop = (data.accuracy > acceptable_local_error);
            // Work out the size of the new element
            new_element = (boundaries[i] + boundaries[i+1])/2;

            if (new_element - boundaries[i] < minimum_acceptable_subdivision_size and loop) loop = false;
            if(loop) boundaries.insert(boundaries.begin() + i + 1, new_element);
            else integral += data.integral;
        }
    }
    return integral;

}
