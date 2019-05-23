#include "integration.h"

// Functions

float SimpsonsRule(float (*f)(float, std::vector<float>), float a, float b, int N, std::vector<float> extraArguments)
{
    // Test and adjustment for evenness.
    try
    {
        if(!checkEven(N))
        {
            std::string message = std::to_string(N);
            throw std::invalid_argument("Odd Number supplied for Simpsons Rule (" + message + ")");
        }
    }
    catch (const std::invalid_argument& ia) {
        std::cerr << ia.what() << "\n" << "Value will be incremented by one." << '\n';
        N++;
    }

    // Create a Grid
    float * grid = (float *) malloc( (N + 1) * sizeof(float)); // Prep array
    float h = (b - a)/((float) N); // Calculate the spacing of the grid, h
    for(int i = 0; i <= N; i++) grid[i] = a + i*h; // fill values

    // Test grid's last element is actually the end
    try
    {
        if(!(grid[N] == b))
        {
            char buffer [50];
            sprintf(buffer, "Last grid element (%f) does not equal the value of b (%f) - it should.", grid[N], b);
            throw std::invalid_argument(buffer);
        }
    }
    catch (const std::invalid_argument& ia) {
        std::cerr << "Critical Error in Simpsons Rule:" << ia.what() << '\n';
        std::exit(1); // This is a critical Fail, so exit.
    }


    // Start the actual Simpsons rule
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


RResult RichardsonExtrapolate(float (*f)(float, std::vector<float>), float a, float b, int steps2, std::vector<float> extraArguments)
{
    // Test and adjustment for evenness.
    try
    {
        if(!checkEven(steps2))
        {
            std::string message = std::to_string(steps2);
            throw std::invalid_argument("Odd Number supplied for Richardson Extrapolation (" + message + ")");
        }
    }
    catch (const std::invalid_argument& ia) {
        std::cerr << ia.what() << "\n" << "Value will be incremented by one." << '\n';
        steps2++;
    }


    // Struct for the return
    RResult data;

    // Call Simpson's rule twice with different intervals accoring to richardson extrapolation.
    float I2n = SimpsonsRule((*f), a, b, steps2, extraArguments);
    float In = SimpsonsRule((*f), a, b, steps2/2, extraArguments);

    // Assign the value of the integral according to the formulae for the integral and it's accuracy
    data.integral = (pow(2, 4) * I2n - In)/(pow(2, 4) - 1);
    data.accuracy = (abs(In - I2n)/(pow(2, 4) - 1));
    return data;
}

float AdaptiveRichardsonExtrapolate(float (*f)(float, std::vector<float>), float a, float b, float accuracy, std::vector<float> extraArguments)
{
    // Preliminary vector for the boundaries.
    std::vector<float> boundaries = {a, b};

    // Maximum (sane) number of subdivisions).
    int max_subdivisions = 1000000;
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
            data = RichardsonExtrapolate((*f), boundaries[i], boundaries[i+1], 4, extraArguments);
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

                try
                {
                    throw std::invalid_argument(" ++ Woah, go easy on the accuracy there; maximum subdivision size reached ++ ");
                }
                catch (const std::invalid_argument& ia) {
                    std::cerr << ia.what() << "\n" << "Integral will move on, but maximum acceptable (local) accuracy was defined as "
                    << acceptable_local_error << " whereas the accuracy of this element is only " << data.accuracy
                    << ". This may or not be a problem for you." << std::endl;
                }
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

float IterativeRichardsonExtrapolate(float (*f)(float, std::vector<float>), float a, float b, float accuracy, std::vector<float> extraArguments)
{
    RResult data;

    // Setter for the number of steps to add each time. Probably should put this as a function argument.
    float steps = 4;
    // Loop until accuracy is okay.
    while(data.accuracy > accuracy)
    {
        data = RichardsonExtrapolate((*f), a, b, steps, extraArguments);
        steps += 4;
    }
    return data.integral;
}
