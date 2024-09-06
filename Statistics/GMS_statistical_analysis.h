
#ifndef __GMS_STATISTICAL_ANALYSIS_H__
#define __GMS_STATISTICAL_ANALYSIS_H__

#include <stddef.h>
            

void run_stat_analysis( const float   * __restrict__  , //ftab
                         const int32_t * __restrict__ , // abscissas
                         float         * __restrict__ ,  // denset
                         float         * __restrict__ , // denset
                         const float     , // denset
                         const float     , // denset
                         const float     ,  // denset
                         const float     , // cubint
                         const float     , // cubint
                         const int32_t   , // number of samples must be a power of 2.
                         const int32_t   , // number of FFT sample points;
                         const char     * __restrict__ );















#endif /*__GMS_STATISTICAL_ANALYSIS_H__*/
