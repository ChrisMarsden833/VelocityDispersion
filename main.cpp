#include "main.h"


int main()
{
    float a = -1;
    float b = 0.5;
    float z = 0.001;

    float mine = incompleteBeta(a, b, z);
    std::cout << "Mine : " << mine << std::endl;

    return 0;
}