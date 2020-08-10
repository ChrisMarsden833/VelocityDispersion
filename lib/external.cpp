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
	    printf("Calculating sigma for %i galaxies\n", size);

	    float * res = (float *) std::malloc(sizeof(float) * size);

        int progress = 0;
        #pragma omp parallel for default(none) shared(size, res, Aperture, Beta, HalfLightRadius, SersicIndex, StellarMass, progress, z, HaloMass, DM, c_path, stdout)
        for (int i = 0; i < size; i++)
        {
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

    float * ParallelSigmaLos(float * R,
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
        printf("Calculating sigma_LOS for %i galaxies\n", size);

        float * res = (float *) std::malloc(sizeof(float) * size);

        int progress = 0;
        #pragma omp parallel for default(none) shared(size, res, R, Beta, HalfLightRadius, SersicIndex, StellarMass, progress, z, HaloMass, DM, c_path, stdout)
        for (int i = 0; i < size; i++)
        {
            if(strcmp(DM, "None") == 0)
            {
                res[i] = GetUnweightedVelocityDispersion(R[i], Beta[i], HalfLightRadius[i], SersicIndex[i], StellarMass[i], z[i]);
            }
            else
            {
                res[i] = GetUnweightedVelocityDispersion(R[i], Beta[i], HalfLightRadius[i], SersicIndex[i], StellarMass[i], z[i], HaloMass[i], DM, c_path);
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
