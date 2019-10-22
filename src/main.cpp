#include "main.h"


int main()
{
    float u = 20;
    float beta = 0.1;

    float res = K_Kernel_DW(u, beta);

    std::cout << res;

    return 0;
}