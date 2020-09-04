#include <iostream>
#include <math.h>
#include "omp.h"
#include <thread>
#include "../src/galaxy.h"

extern "C"
{
    float * ParallelSigma(float * Aperture,
                           float * redshift,
                           float * bulge_mass,
                           float * bulge_radius,
                           float * bulge_beta,
                           float * bulge_sersicIndex,
                           int * componentFlag,
                           float * disk_mass,
                           float * disk_inclination,
                           float * Halo_mass,
                           char * profile_name,
                           char * c_path,
                           float * BlackHole_mass,
                           int size)
    {
        printf("\n\n######################################## \n");
	    printf("Chris Marsden's Velocity dispersion code \n");
	    printf("Calculating sigma ap for %i galaxies\n", size);

	    float * res = (float *) std::malloc(sizeof(float) * size);

        int progress = 0;
        #pragma omp parallel for default(none) shared(progress, size, res, Aperture, redshift, bulge_mass, bulge_radius, bulge_beta, bulge_sersicIndex, componentFlag, disk_mass, disk_inclination, Halo_mass, profile_name, c_path, BlackHole_mass, size, stdout) schedule(dynamic, 1)
        for (int i = size-1; i >= 0; i--)
        {

            res[i] = GetVelocityDispersion(Aperture[i], redshift[i], bulge_mass[i], bulge_radius[i], bulge_beta[i], bulge_sersicIndex[i], componentFlag,
                    disk_mass[i], disk_inclination[i], Halo_mass[i], profile_name, c_path, BlackHole_mass[i]);

	    // Print progress
        #pragma omp critical
        {
            progress++;
            printf("\r %i of %i calculated | %2.4f%%", progress, size, 100. * float(progress) / float(size));
            fflush(stdout);
        };
	    }
        return res;
    }


}
