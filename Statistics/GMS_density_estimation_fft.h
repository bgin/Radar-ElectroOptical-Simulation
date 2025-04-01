
#ifndef __GMS_DENSITY_ESTIMATION_FFT_H__
#define __GMS_DENSITY_ESTIMATION_FFT_H__


#include <stddef.h>

/* Density estimation algorithm AS 176
   Converted to C99 from Fortran 77.
*/


void denest(float * __restrict__,
            const int32_t,
            float,
            float,
            float,
            float * __restrict__,
            float * __restrict__,
            const int32_t,
            int32_t,
            int32_t * __restrict__);







#endif /*__GMS_DENSITY_ESTIMATION_FFT_H__*/
