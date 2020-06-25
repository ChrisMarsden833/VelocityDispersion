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
                       int size)
    {

	printf("Commencing parallelization of %i elements\n", size);


	float * res = (float *) std::malloc(sizeof(float) * size);

        int progress = 0;

        #pragma omp parallel for default(none) shared(size, res, Aperture, Beta, HalfLightRadius, SersicIndex, StellarMass, progress, stdout)
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

            res[i] = GetVelocityDispersion(Aperture[i], Beta[i], HalfLightRadius[i], SersicIndex[i], StellarMass[i]);

	    // Timer
	    if(false)
	    {
            	#pragma omp atomic
                progress++;

            	if(omp_get_thread_num() == 0)
            	{
                	printf("\r%2.4f %%", 100.*float(progress)/float(size));
                	fflush( stdout );
            	}
            }
	}
        return res;
    }
}
