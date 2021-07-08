#include "blackhole.h"

BlackHole::BlackHole(float mass, bool grav, float spin, float charge)
{
    // No Hair
    grav = grav;
    mass = mass;
}

float BlackHole::getMass(void){
    if(grav) return mass;
    return 0.0;
}