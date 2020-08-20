#include "utillity.h"

template<typename T>
T factorial(T n)
{
    if (n == 0)
        return 1; // base case
    else
        return n * factorial(n-1); // recursive case
}

bool checkEven(int Value)
{
    // Test N to make sure it is even.
    if(Value & 1)
    {
        return false;
    }
    else {
        return true;
    }

}

float incompleteBeta(float a, float b, float z)
{
    float beta;

    if(a <= 0)
    {
        beta = pow(z, a) * pow((1 - z), b) + (a + b) * incompleteBeta(a + 1, b, z);
        beta = beta/a;
        return beta;
    }
    else
    {
        assert( ((a >= 0) && (a <= 1)) || assert_msg("incompleteBeta recursion finished but a = " << a << std::endl <<
          "====Location: inCompleteBeta" << std::endl <<
          "----a = " << a << std::endl <<
          "----b = " << b << std::endl <<
          "----z = " << z << std::endl));

        assert( ((z >= 0) && (z <= 1)) || assert_msg("incompleteBeta about to be called, but z = " << z << std::endl <<
            "====Location: inCompleteBeta" << std::endl <<
            "----a = " << a << std::endl <<
            "----b = " << b << std::endl <<
            "----z = " << z << std::endl));


        beta = boost::math::beta(a, b, z);


        return beta;
    }
}

void CheckAndMaybeIncrementN(int * N, std::string fname)
{

    if(!checkEven(*N))
    {
        ++*N;
    }
    /*
    try
    {
        if(!checkEven(*N))
        {
            std::string message = std::to_string(*N);
            throw std::invalid_argument("Odd Number supplied for " + fname + " (" + message + ")");
        }
    }
    catch (const std::invalid_argument& ia) {
        std::cerr << ia.what() << "\n" << "Value will be incremented by one." << '\n';
        ++*N;
    }*/
}

void fastLinspace(float * &grid, float &h, float a, float b, int N)
{
    grid = (float *) malloc( (N + 1) * sizeof(float)); // Prep array
    h = (b - a)/((float) N); // Calculate the spacing of the grid, h
    for(int i = 0; i <= N; i++) grid[i] = a + h * (float) i; // fill values
}

std::vector<std::vector<float>> * ReadFile(std::string path, std::vector<int> * IndexesToGrab)
{
    auto start = std::chrono::high_resolution_clock::now(); // Timing

    InputPrechecks(IndexesToGrab);


    // Decleration of output arrays.
    std::vector<std::vector<float>> * OutputArrays = new std::vector<std::vector<float>>(IndexesToGrab->size());

    // Create file pointer and find size
    FILE * fp;
    fp = fopen(path.c_str(), "r");
    if(fp == NULL)
    {
	    #pragma omp critical
	    {
		    std::cout << "Error, file could not be found/opened" << std::endl;
		    std::cout << "    Path:" << path << std::endl;
		    exit(1);
	    }
    }
    // Calculate size of file
    fseek(fp, 0, SEEK_END);
    int fileSize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char * buffer = (char*) malloc(sizeof(char)*fileSize);
    int newFileSize = fread(buffer,1,fileSize, fp); // Read in the entire thing -- insanity.

    // Variables used proper

    // comment syntax
    // TODO implement this as a function argument
    char commentSyntax = '#';
    char delimiter = ' ';

    // Variables used for the processing.
    char dummy;

    std::string dummyStr;

    int commentCount = 0;
    int wordCount = 0;
    int columnIndex = 0;
    int rowIndex = 0;

    bool inWord, inComment;

    float floatBuffer;

    inWord = false;
    inComment = false;

    std::vector<float> ExtractionArray;


    for(int i = 0; i < fileSize; i++)
    {
       dummy = buffer[i];

       if(dummy == commentSyntax) // if we are at the start of a comment
       {
           inComment = true;
       }
       else if(inComment == true) // are we in a comment
       {
           if(dummy == '\n') // and at the end of a line?
           {
               inComment = false; // then the comment ends
               commentCount++;
           }
       }
       else if(inWord == false) // we are not in comment, and not in word
       {
           if(buffer[i] != delimiter)
           {
               inWord = true; // new word start
               dummyStr += buffer[i];

           }
           // else we are in a delimiter, do nothing.
       }
       else if(inWord == true) // if we were in a word...
       {
           if(buffer[i] != delimiter && buffer[i] != '\n') // and still in the word.
           {
               dummyStr += buffer[i];
           }
           else // and are now at the end of said word
           {

               floatBuffer = std::stof(dummyStr);
               for(int j = 0; j < IndexesToGrab->size(); j++)
               {
                   if(columnIndex == IndexesToGrab->at(j))
                   {
                       OutputArrays->at(j).push_back(floatBuffer);
                   }
               }

               dummyStr.clear();
               columnIndex++;

               if(buffer[i] == '\n')
               {
                   columnIndex = 0;
                   rowIndex++;
               }

               wordCount++;
               inWord = false;
           }

       }
    }

    free(buffer);

    // Test to ensure that lengths are consistent
    checkOutputConsistency(IndexesToGrab, OutputArrays);

    fclose(fp);

    auto finish = std::chrono::high_resolution_clock::now(); // Timing
    std::chrono::duration<double> elapsed = finish - start;

    return OutputArrays;

}

void checkOutputConsistency(std::vector<int> *IndexesToGrab, std::vector<std::vector<float>> *OutputArrays)
{
    // function to check output arrays are all the same length

    std::vector<int> sizes;

    int remaining_sizes = IndexesToGrab->size();

    for(int i = 0; i < remaining_sizes; i++)
    {
        sizes.push_back(OutputArrays->at(i).size());
    }
    std::sort(sizes.begin(), sizes.end());

    auto last2 = std::unique(sizes.begin(), sizes.end());

    sizes.erase(last2, sizes.end());

    if(sizes.size() != 1)
    {
        printf("Error - output vectors are different lengths");
    }

}

void InputPrechecks(std::vector<int> *IndexesToGrab)
{
    int previous = 0;

    for(int i = 0; i < IndexesToGrab->size(); i++)
    {
        if(previous > IndexesToGrab->at(i))
        {
            printf("Warning, supplied index array is not in ascending order - this will be reordered\n");
            break;
        }
        previous = IndexesToGrab->at(i);
    }
    std::sort(IndexesToGrab->begin(), IndexesToGrab->end());
    auto last = std::unique(IndexesToGrab->begin(), IndexesToGrab->end());
    int originalLength = IndexesToGrab->size();
    IndexesToGrab->erase(last, IndexesToGrab->end());
    bool duplicatesFlag = (IndexesToGrab->size() < originalLength);
    if(duplicatesFlag) printf("Warning, supplied index array contains duplicates - duplicates will be removed\n");
}

std::vector<float> * Reduce(std::vector<float> * input, std::vector<float> * output)
{
    float previous, current;
    output->clear();

    for(int i=0; i < input->size(); i++)
    {
	if(i==0)
	{
	    previous = input->at(i);
	    output->push_back(previous);
	}
	else
	{
	    current = input->at(i);
	    if(current != previous) output->push_back(current);
	}
    }
    return output;
}

float FindClosest(float value, std::vector<float> * data)
{
    float closest_distance, closest_value;

    for(int i=0; i < data->size(); i++)
    {
	if(i==0)
	{
	    closest_value = data->at(i);
	    closest_distance = abs(closest_value - value);
	}
	else
	{
	    if(abs(data->at(i) - value) < closest_distance)
	    {
	        closest_distance = abs(data->at(i)-value);
		closest_value = data->at(i);
	    }
	}
    }
    delete data;
    return closest_value;
}

void Equals(std::vector<float> * array, float value, std::vector<bool> * mask)
{
	mask->clear();
	for(int i=0; i<array->size(); i++)
	{
		bool input;
		if(array->at(i) == value) input = true;
		else input = false;
		mask->push_back(input);
	}
}

void MaskOut(std::vector<float> * array, std::vector<bool> * mask)
{
	std::vector<float> input(*array);
	array->clear();
	for(int i=0; i<input.size(); i++)
	{
		if(mask->at(i)) array->push_back(input[i]);
	}
}

float LinearInterp(std::vector<float> * X, std::vector<float> * Y, float x)
{
	if(X->size() != Y->size())
	{
		std::cout << "LinearInterp found X and Y of different lengths:" << std::endl;
	      	std::cout << "X: " << X->size() << std::endl;
		std::cout << "Y: " << Y->size() << std::endl;	
		exit(0);
	}


	float x0, x1, y0, y1;
	// Find the above values
	for(int i=1; i<X->size(); i++)
	{
		x0 = X->at(i-1);
		x1 = X->at(i);
		if(x < x1 && x > x0)
		{
			y0 = Y->at(i-1);
			y1 = Y->at(i);
			break;
		}
	}
	
	return y0*(1 - (x-x0)/(x1-x0)) + y1*(x-x0)/(x1-x0);

}

