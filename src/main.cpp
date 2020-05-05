#include "main.h"


int main()
{

    auto start = high_resolution_clock::now();

    Galaxy aGalaxy(12., 0.4, 10., 10., 4.);
    float test = aGalaxy.sigma_ap();

    auto stop = high_resolution_clock::now(); 

    std::cout << "Result: " << test << std::endl;

    auto duration = duration_cast<microseconds>(stop - start);

    cout << "Time taken by function: "  << duration.count() << " microseconds" << endl;


    start = high_resolution_clock::now();

    float old_test = sigma_aperture(10., 0.4, 10., 4, pow(10., 12.), 10., 10.);

    stop = high_resolution_clock::now(); 
    
    std::cout << "Result2: " << old_test << std::endl;
   
    duration = duration_cast<microseconds>(stop - start);

    cout << "Time taken by function: "  << duration.count() << " microseconds" << endl;


    return 0;


}
