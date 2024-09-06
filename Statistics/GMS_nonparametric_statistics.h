
#ifndef __GMS_NONPARAMETRIC_STATISTICS_H__
#define __GMS_NONPARAMETRIC_STATISTICS_H__


#include <stddef.h>


void estimate_1st_moment(const float   * __restrict__ rdcyc_data, //ftab
                         const int32_t * __restrict__ xtab, // abscissas
                         float         * __restrict__ ft,  // denset
                         float         * __restrict__ smooth, // denset
                         const float     lo, // denset
                         const float     hi, // denset
                         const float     window,  // denset
                         const float     a, // cubint
                         const float     b, // cubint
                         const int32_t   nsamp, // number of samples must be a power of 2.
                         const int32_t   nft,   // number of FFT sample points
                         float   * __restrict__ mean,
                         int32_t * __restrict__ error,    // cubint error
                         int32_t * __restrict__ denserr); // denset error
                         





















#endif /*__GMS_NONPARAMETRIC_STATISTICS_H__*/
