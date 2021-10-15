#ifndef BH_H
#define BH_H

#include <iostream>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <math.h>
using namespace std;

// Super simple class for the black hole. It's overkill, but in the spirit of consistency, here we are.
class BlackHole
{
    public:
        /* Black hole constructor - we only care about the mass and if it's gravity is switched on.
        The spin and charge options are purely for comedy value - we can now parameterize all possible 
        black hole properties (no hair theorem)!
        */
        BlackHole(float Mass, bool grav, float spin = 0., float charge = 0.);
        float getMass(void);
    private:
        bool internal_grav = false;
        float black_hole_mass;
        
};

#endif