#include <iostream>
#include <math.h>
#include "omp.h"
#include <thread>
#include "../src/desmond.h"
#include "../src/galaxy.h"

extern "C"
{
    float * ParallelSigma(float * Aperture,
                       float * Beta,
                       float * HalfLightRadius,
                       float * SersicIndex,
                       float * StellarMass,
                       float * HaloMass,
                       float * z,
                       int size,
                       char * DM,
                       char * c_path)
    {
        printf("\n\n######################################## \n");
	    printf("Chris Marsden's Velocity dispersion code \n");
	    printf("Calculating for %i galaxies\n", size);

	    float * res = (float *) std::malloc(sizeof(float) * size);

        int progress = 0;
        #pragma omp parallel for default(none) shared(size, res, Aperture, Beta, HalfLightRadius, SersicIndex, StellarMass, progress, z, HaloMass, DM, c_path, stdout)
        for (int i = 0; i < size; i++)
        {
            int thread = omp_get_thread_num();
            if(false) {
                printf("%i: Ap: %f\n", thread, Aperture[i]);
                printf("%i, Bta: %f\n", thread, Beta[i]);
                printf("%i, HLR: %f\n", thread, HalfLightRadius[i]);
                printf("%i, n: %f\n", thread, SersicIndex[i]);
                printf("%i, Sm: %f\n", thread, StellarMass[i]);
            }
            if(strcmp(DM, "None") == 0)
            {
                res[i] = GetVelocityDispersion(Aperture[i], Beta[i], HalfLightRadius[i], SersicIndex[i], StellarMass[i], z[i]);
            }
            else
            {
                res[i] = GetVelocityDispersion(Aperture[i], Beta[i], HalfLightRadius[i], SersicIndex[i], StellarMass[i], z[i], HaloMass[i], DM, c_path);
            }


	    // Timer

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
