

#ifndef __GMS_DASUM_HPP__
#define __GMS_DASUM_HPP__

//=========================================================================
// Modified and slightly optimized version of OpenBLAS dasum
// Programmer: Bernard Gingold, contact: beniekg@gmail.com
// 06-03-2021 10:25  +00200
// Original copyright below
//========================================================================
/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/


#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"

#ifndef DASUM_ABS_K
#define DASUM_ABS_K(a) ((a) > 0 ? (a) : (-(a)))
#endif

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
__attribute__((no_stack_protector))
static inline
double dasum_kernel_avx2(const int32_t,
                         double * __restrict x1) {

       int32_t i = 0;
       double sum = 0.0;
       if (n >= 256) {
       const int32_t align_256 = ((32 - ((uintptr_t)x1 & (uintptr_t)0x1f)) >> 3) & 0x3;

        for (i = 0; i < align_256; i++) {
            sumf += DASUM_ABS_K(x1[i]);
        }

        n -= align_256;
        x1 += align_256;
    }

    const int32_t tail_index_SSE = n&(~7);
    const int32_t tail_index_AVX2 = n&(~255);

    if (n >= 256) {
        __m256d accum_0, accum_1, accum_2, accum_3;

        accum_0 = _mm256_setzero_pd();
        accum_1 = _mm256_setzero_pd();
        accum_2 = _mm256_setzero_pd();
        accum_3 = _mm256_setzero_pd();

         const __m256i abs_mask = _mm256_set1_epi64x(0x7fffffffffffffff);
        for (i = 0; i < tail_index_AVX2; i += 16) {
            accum_0 += (__m256d)_mm256_and_si256(_mm256_load_si256(&x1[i+ 0]), abs_mask);
            accum_1 += (__m256d)_mm256_and_si256(_mm256_load_si256(&x1[i+ 4]), abs_mask);
            accum_2 += (__m256d)_mm256_and_si256(_mm256_load_si256(&x1[i+ 8]), abs_mask);
            accum_3 += (__m256d)_mm256_and_si256(_mm256_load_si256(&x1[i+12]), abs_mask);
        }

        accum_0 = accum_0 + accum_1 + accum_2 + accum_3;

        __m128d half_accum0;
        half_accum0 = _mm_add_pd(_mm256_extractf128_pd(accum_0, 0), _mm256_extractf128_pd(accum_0, 1));

        half_accum0 = _mm_hadd_pd(half_accum0, half_accum0);

        sum += half_accum0[0];
    }
    
    if (n >= 8) {
        __m128d accum_20, accum_21, accum_22, accum_23;
        accum_20 = _mm_setzero_pd();  
        accum_21 = _mm_setzero_pd(); 
        accum_22 = _mm_setzero_pd(); 
        accum_23 = _mm_setzero_pd(); 

        __m128i abs_mask2 = _mm_set1_epi64x(0x7fffffffffffffff);
        for (i = tail_index_AVX2; i < tail_index_SSE; i += 8) {
            accum_20 += (__m128d)_mm_and_si128(_mm_loadu_si128(&x1[i + 0]), abs_mask2);
            accum_21 += (__m128d)_mm_and_si128(_mm_loadu_si128(&x1[i + 2]), abs_mask2);
            accum_22 += (__m128d)_mm_and_si128(_mm_loadu_si128(&x1[i + 4]), abs_mask2);
            accum_23 += (__m128d)_mm_and_si128(_mm_loadu_si128(&x1[i + 6]), abs_mask2);
        }

        accum_20 = accum_20 + accum_21 + accum_22 + accum_23;
        __m128d half_accum20;
        half_accum20 = _mm_hadd_pd(accum_20, accum_20);

        sum += half_accum20[0];
    }
    
    for (i = tail_index_SSE; i < n; ++i) {
        sum += DASUM_ABS_K(x1[i]);
    }

    return sum;

}


#if defined(__AVX512F__)

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
__attribute__((no_stack_protector))
static inline
double dasum_kernel_avx512(const int32_t n,
                           double * __restrict x1) {

       int32_t i = 0;
       double sum = 0.0;
       if (n >= 256) {
          const int32_t align_512 = ((64 - ((uintptr_t)x1 & (uintptr_t)0x3f)) >> 3) & 0x7;

        for (i = 0; i < align_512; i++) {
            sumf += DASUM_ABS_K(x1[i]);
        }
        
        n -= align_512;
        x1 += align_512;
    }

    const int32_t tail_index_SSE = n&(~7);
    const int32_t tail_index_AVX512 = n&(~255);

    //
    if ( n >= 256 ) {

        __m512d accum_0, accum_1, accum_2, accum_3;
        accum_0 = _mm512_setzero_pd();
        accum_1 = _mm512_setzero_pd();
        accum_2 = _mm512_setzero_pd();
        accum_3 = _mm512_setzero_pd();
        for (i = 0; i < tail_index_AVX512; i += 32) {
            accum_0 += _mm512_abs_pd(_mm512_load_pd(&x1[i + 0]));
            accum_1 += _mm512_abs_pd(_mm512_load_pd(&x1[i + 8]));
            accum_2 += _mm512_abs_pd(_mm512_load_pd(&x1[i +16]));
            accum_3 += _mm512_abs_pd(_mm512_load_pd(&x1[i +24]));
        }

        accum_0 = accum_0 + accum_1 + accum_2 + accum_3;
        sum += _mm512_reduce_add_pd(accum_0);
    }

    if (n >= 8) {
        __m128d accum_20, accum_21, accum_22, accum_23;
        accum_20 = _mm_setzero_pd();  
        accum_21 = _mm_setzero_pd(); 
        accum_22 = _mm_setzero_pd(); 
        accum_23 = _mm_setzero_pd(); 

        __m128i abs_mask2 = _mm_set1_epi64x(0x7fffffffffffffff);
        for (i = tail_index_AVX512; i < tail_index_SSE; i += 8) {
            accum_20 += (__m128d)_mm_and_si128(_mm_loadu_si128(&x1[i + 0]), abs_mask2);
            accum_21 += (__m128d)_mm_and_si128(_mm_loadu_si128(&x1[i + 2]), abs_mask2);
            accum_22 += (__m128d)_mm_and_si128(_mm_loadu_si128(&x1[i + 4]), abs_mask2);
            accum_23 += (__m128d)_mm_and_si128(_mm_loadu_si128(&x1[i + 6]), abs_mask2);
        }

        accum_20 = accum_20 + accum_21 + accum_22 + accum_23;
        __m128d half_accum20;
        half_accum20 = _mm_hadd_pd(accum_20, accum_20);

        sum += half_accum20[0];
    }

    for (i = tail_index_SSE; i < n; ++i) {
        sum += DASUM_ABS_K(x1[i]);
    }

    return sum;
}


#endif


__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
__attribute__((no_stack_protector))
static inline
double dasum(const int32_t n,
             double * __restrict x,
	     const int32_t inc_x) {
 
      if(n<=0 || inc_x<=0) {return(0.0);}
      int32_t i = 0;
      double sum = 0.0;
      if( inc_x == 1 ) {
#if defined(__AVX512F__)
        sum = dasum_kernel_avx512(n, x);
#else
        sum = dasum_kernel_avx2(n,x);
#endif
    } 
    else {
        n *= inc_x;
       
        while(i < n) {
            sum += DASUM_ABS_K(x[i]);
            i += inc_x;
        }
    }
    return(sum);
}




#endif /*__GMS_DASUM_HPP__*/
