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
                          float * BlackHoleMass,
                          float * z,
                          int size,
                          char * DM,
                          char * c_path,
                          int * componentFlag)
    {
        printf("\n\n######################################## \n");
	    printf("Chris Marsden's Velocity dispersion code \n");
	    printf("Calculating sigma ap for %i galaxies\n", size);

	    float * res = (float *) std::malloc(sizeof(float) * size);

        int progress = 0;
        #pragma omp parallel for default(none) shared(size, res, Aperture, Beta, HalfLightRadius, SersicIndex, StellarMass, progress, z, HaloMass, DM, HaloMass, BlackHoleMass, componentFlag, c_path, stdout) schedule(static, 2)
        for (int i = 0; i < size; i++)
        {

            res[i] = GetVelocityDispersion(Aperture[i], Beta[i], HalfLightRadius[i], SersicIndex[i], StellarMass[i], HaloMass[i], BlackHoleMass[i], z[i], DM, c_path, componentFlag);
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
                             float * BlackHoleMass,
                             float * z,
                             int size,
                             char * DM,
                             char * c_path,
                             int * componentFlag)
    {
        printf("\n\n######################################## \n");
        printf("Chris Marsden's Velocity dispersion code \n");
        printf("Calculating sigma_LOS for %i galaxies\n", size);

        float * res = (float *) std::malloc(sizeof(float) * size);

        int progress = 0;
        #pragma omp parallel for default(none) shared(size, res, R, Beta, HalfLightRadius, SersicIndex, StellarMass, progress, z, HaloMass, BlackHoleMass, DM, c_path, componentFlag, stdout) schedule(static, 2)
        for (int i = 0; i < size; i++)
        {

            res[i] = GetUnweightedVelocityDispersion(R[i], Beta[i], HalfLightRadius[i], SersicIndex[i], StellarMass[i],
                    HaloMass[i], BlackHoleMass[i], z[i], DM, c_path, componentFlag);

            #pragma omp critical
            {
                progress++;
                printf("\r %i of %i calculated | %2.4f%%", progress, size, 100. * float(progress) / float(size));
                fflush(stdout);
            };
        }
        return res;

    }

    float * GetCumMass(float * R,
                             float * Beta,
                             float * HalfLightRadius,
                             float * SersicIndex,
                             float * StellarMass,
                             float * HaloMass,
                             float * z,
                             int size,
                             char * DM,
                             char * c_path,
                             int flag)
    {


        float * res = (float *) std::malloc(sizeof(float) * size);

        int progress = 0;
        #pragma omp parallel for default(none) shared(size, res, R, Beta, HalfLightRadius, SersicIndex, StellarMass, progress, z, HaloMass, DM, c_path, flag, stdout)
        for (int i = 0; i < size; i++)
        {
            res[i] = GetCum(R[i], Beta[i], HalfLightRadius[i], SersicIndex[i], StellarMass[i], z[i], HaloMass[i], DM, c_path);
            #pragma omp critical
            {
                progress++;
                printf("\r %i of %i calculated | %2.4f%%", progress, size, 100. * float(progress) / float(size));
                fflush(stdout);
            };

        }
        return res;

    }

    float * DM_Profile(float * r, float HaloMass, float z, int size, char * DM, char * c_path)
    {
        printf("\n\n######################################## \n");
        printf("Chris Marsden's Velocity dispersion code \n");
        printf("Calculating DM Profile for %i points\n", size);

        float * res = (float *) std::malloc(sizeof(float) * size);

        for (int i = 0; i < size; i++)
        {
            res[i] = GetDMrho(r[i], 0.0014, 5, 3.223, 11.0, z, HaloMass, DM, c_path);
        }
        return res;
    }
}
