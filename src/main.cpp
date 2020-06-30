#include "main.h"


int main()
{

    float z = 0.149;
    float mass = 3.15e+08;

    
    std::string FilePath = "/Users/chris/Desktop/cM_planck18.txt";

    std::vector<int> * IndexesToGrab = new std::vector<int>(0);
    IndexesToGrab->push_back(0); // z
    IndexesToGrab->push_back(2); // M200c
    IndexesToGrab->push_back(4); // c200c

    std::vector<std::vector<float>> * Extracted;

    Extracted = ReadFile(FilePath, IndexesToGrab);

    printf("File Read\n");

    std::vector<float> * Redshift = &Extracted->at(0);
    std::vector<float> * M200c = &Extracted->at(1);
    std::vector<float> * c200c = &Extracted->at(2);

    std::vector<float> * reduced = new std::vector<float>;
    float closest_z = FindClosest(z,  Reduce(Redshift, reduced));
    std::vector<bool> * mask = new std::vector<bool>;
    Equals(Redshift, closest_z, mask);
    
    MaskOut(Redshift, mask);
    MaskOut(M200c, mask);
    MaskOut(c200c, mask);

    float logc = LinearInterp(M200c, c200c, mass);

    std::cout << logc << std::endl; 
    
    return 0;


}
