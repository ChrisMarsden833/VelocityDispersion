#include "main.h"


int main()
{
    float mine;

    float R_ap = 500;
    float beta = 0.01;
    float Half_Light_radius = 15000;
    float SersicIndex = 4;

    clock_t start = clock();

    mine = K_Kernel_DW(1., beta);

    //mine = sigma_aperture(R_ap, beta, Half_Light_radius, SersicIndex);

    clock_t end = clock();


    double time = (double) (end-start) / CLOCKS_PER_SEC;

    std::cout << "Time : " << time << std::endl;


    std::cout << "Result :" << mine << std::endl;

    return 0;
}