#include <iostream>
#include <math.h>
#include "omp.h"
#include <thread>
#include "../src/galaxy.h"

/* This file serves two purposes:
 * 1. To act as a `wrapper' function, written entirely in C so that it can be called by the ctypes python library
 * 2. To facilitate multiprocessing via OpenMP
 *
 * This code accepts an array for each (unique) parameter, and for each element calls the fully velocity dispersion code.
 * */

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
                           float * disk_scale_length,
                           float * Halo_mass,
                           char * profile_name,
                           float * halo_concentration,
                           float * BlackHole_mass,
                           int mode,
                           int threads,
			   int debug,
                           int size)
    {
        // Little header so we know we've started.
	if(debug == 1)
	{
        	printf("\n\n######################################## \n");
	    	printf(    "Chris Marsden's Velocity dispersion code \n");
	}

        // Reserved memory for results, and progress.
	    float * res = (float *) std::malloc(sizeof(float) * size);
        int progress = 0;

        // Set parallelism up.
        if(threads < 1) threads = omp_get_max_threads();
        if (debug == 1) printf("Preparing to calculate sigma for %i galaxies, using %i threads\n", size, threads);
        omp_set_num_threads(threads);

        // Shedulding is set to schedule(dynamic, 1). This is because tasks can be very asymmetric.
        // Loop back through array as often the last elements can take the longest, so this results in better load management with the above schedule
        #pragma omp parallel for default(none) shared(progress, res, Aperture, redshift, bulge_mass, bulge_radius, bulge_beta, bulge_sersicIndex, componentFlag, disk_mass, disk_inclination, disk_scale_length, Halo_mass, profile_name, halo_concentration, BlackHole_mass, mode, threads, debug, size, stdout) schedule(dynamic, 1)
        for (int i = size-1; i >= 0; i--)
        {
            // Call function itself.
        	res[i] = GetVelocityDispersion(Aperture[i], redshift[i], bulge_mass[i], bulge_radius[i], bulge_beta[i], bulge_sersicIndex[i], componentFlag,
                    disk_mass[i], disk_inclination[i], disk_scale_length[i], Halo_mass[i], profile_name, halo_concentration[i] , BlackHole_mass[i], mode);
	
		if(debug == 1)
		{
		// Print progress
        	#pragma omp critical
        	{
            		progress++;
            		printf("\r %i of %i calculated | %2.4f%%", progress, size, 100. * float(progress) / float(size));
            		fflush(stdout);
        	};
		}
        }
	
	if(debug == 1)
	{
		printf("\nCompleted\n");
        	printf("\n\n######################################## \n");
	}

	// Return the array to the user.
        return res;
    }


}
