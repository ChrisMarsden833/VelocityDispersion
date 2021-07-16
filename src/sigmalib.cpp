#include "sigmalib.h"


float Sigma(float aperture,
            float bulge_mass_or_luminosity,
            float bulge_radius,
            float bulge_sersic_index,
            float bulge_beta,
            float disk_mass,
            float disk_h,
            float disk_i,
            const char * halo_profile,
            float haloRs,
            float haloRhos,
            float black_hole_mass,
            int * tracer_flags,
            int * gravitational_flags,
            bool use_variable_IMF,
            float varIMF_R1,
            float varIMF_R0,
            int mode)
{
    // Manage which components contribute gravitationally, and which do not.
    bool trace_bulge = true;
    bool trace_disk = true;

    bool grav_bulge = true;
    bool grav_disk = true;
    bool grav_halo = true;
    bool grav_bh = true;

    if(tracer_flags) {
        trace_bulge = (bool) tracer_flags[0];
        trace_disk = (bool) tracer_flags[1];
    }
    if(gravitational_flags)    {
        grav_bulge = (bool) gravitational_flags[0];
        grav_disk = (bool) gravitational_flags[1];
        grav_halo = (bool) gravitational_flags[2];
        grav_bh = (bool) gravitational_flags[3];
    }

    // Create Halo object
    Halo * aHalo = new Halo(haloRhos, haloRs, halo_profile, grav_halo);
    // Create black hole object
    BlackHole * aBlackHole = new BlackHole(black_hole_mass, grav_bh);

    // Create Bulge object
    Bulge * aBulge = new Bulge(bulge_mass_or_luminosity, bulge_radius, bulge_sersic_index, 
        bulge_beta, trace_bulge, grav_bulge, aHalo, aBlackHole, use_variable_IMF, varIMF_R0, 
        varIMF_R1);
    // Create disk object
    Disk * aDisk = new Disk(disk_mass, disk_h, disk_i, trace_disk, grav_disk, aBulge, aHalo);

    float res;

    if(mode == 1){
        // Get bulge + disk velocity dispersion.
        res = sqrt( pow(aBulge->sigma_ap(aperture), 2.) + aDisk->VelocityDispersion(aperture) ); // add disk in quadrature
    }
    else if(mode == 2)
    {
        res = aBulge->sigma_LOS(aperture); // LOS velocity dispersion
    }
    else if(mode == 3)
    {
        res = aBulge->sigma_LOS_full(aperture);  // Densities, other stuff is for testing.
    }
    else if(mode == 4)
    {
        res = aBulge->prug_sim_profile(aperture) + aBulge->rhoX(aperture); 
    }
    else if(mode == 5)
    {
        res = aBulge->mass_within(aperture);
    }
    else if(mode == 6)
    {
        res = aBulge->get_density(aperture);
    }
    else if(mode == 7)
    {
        res = aBulge->sersic_profile(aperture);
    }
    else{
        cout << "Unknown mode " << mode << endl;
    }
    
    return res;

}