#include "main.h"


int main()
{
    float mine;

    float R_ap = 500;
    float beta = 1.;
    float Half_Light_radius = 15000;
    float SersicIndex = 4;

    clock_t start = clock();

    mine = sigma_los(40., 0.51, 30., 3.);

    clock_t end = clock();


    double time = (double) (end-start) / CLOCKS_PER_SEC;

    std::cout << "Time : " << time << std::endl;


    std::cout << "Result :" << mine << std::endl;

    return 0;
}