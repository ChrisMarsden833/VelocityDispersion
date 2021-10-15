#include "blackhole.h"

BlackHole::BlackHole(float mass, bool grav, float spin, float charge)
{
    // No Hair
    internal_grav = grav;
    black_hole_mass = mass;
}

float BlackHole::getMass(void){
    if(internal_grav) return pow(10., black_hole_mass);
    return 0.0;
}