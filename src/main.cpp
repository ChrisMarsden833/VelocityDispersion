#include "main.h"


int main()
{
    Bulge * aBulge = new Bulge(1e11/4., 3.0, 4.0, 0.0, true, true, nullptr, nullptr, true, 8.0, 3.5);

    float mass_num = aBulge->mass_within(1000);
    float mass_ana = aBulge->projected_mass_within(1000);

    cout << "Mass comparison : " << log10(mass_num) << " " << log10(mass_ana) << endl;

    float sigma = aBulge->sigma_LOS(0.12 * 3.0);

    cout << "sigma = " << sigma << " (log10: " << log10(sigma) << ")" << endl;



    return 0;

}
