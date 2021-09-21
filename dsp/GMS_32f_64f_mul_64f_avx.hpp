

#ifndef __GMS_32F_64F_MUL_64F_AVX_HPP__
#define __GMS_32F_64F_MUL_64F_AVX_HPP__


/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 13-12-2020 09:51AM +00200
    contact: beniekg@gmail.com
    Few modification were added to original
    implementation (ICC pragmas, alignment directives and code layout rescheduled)
    
*/

/*
 * Copyright 2018 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"

#if defined(DSP_32F_64F_MUL_64F_AVX_BLOCK)
    #define DSP_32F_64F_MUL_64F_AVX_BLOCK                                   \
               int32_t idx = 0;                                             \
	       const int32_t len = npoints/8;                               \
	       register __m256 aVal;                                        \
               register  __m128 aVal1, aVal2;                               \
               register __m256d aDbl1, aDbl2, bVal1, bVal2, cVal1, cVal2;
     #endif

 
         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_64f_mul_64f_u_avx_looped(double * __restrict c,
	                                float  * __restrict b,
					float  * __restrict a,
					const int32_t npoints) {
              DSP_32F_64F_MUL_64F_AVX_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
              for(; idx != len; ++idx) {
		    _mm_prefetch((const char*)&a+32,_MM_HINT_T0);
                   aVal = _mm256_loadu_ps(a);
		   aVal1 = _mm256_extractf128_ps(aVal, 0);
		   aDbl1 = _mm256_cvtps_pd(aVal1);
                   aVal2 = _mm256_extractf128_ps(aVal, 1);
		   aDbl2 = _mm256_cvtps_pd(aVal2);
		     _mm_prefetch((const char*)&b+32,_MM_HINT_T0);  
		   bVal1 = _mm256_loadu_pd(b);
		   cVal1 = _mm256_mul_pd(aDbl1, bVal1);
                   bVal2 = _mm256_loadu_pd(b+4);
                   cVal2 = _mm256_mul_pd(aDbl2, bVal2);
		   _mm256_storeu_pd(c,cVal1);
                   _mm256_storeu_pd(c+4,  cVal2);      
                   a += 8;
                   b += 8;
                   c += 8;
	      }
   idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
               for(; idx != npoints; ++idx) {
                     c[i] = ((double)a[i])*b[i];
	       }	      
	 }

	 
         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_64f_mul_64f_a_avx_looped(double * __restrict __ATTR_ALIGN__(32) c,
	                                float  * __restrict __ATTR_ALIGN__(32) b,
					float  * __restrict __ATTR_ALIGN__(32) a,
					const int32_t npoints) {
              DSP_32F_64F_MUL_64F_AVX_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
              __assume_aligned(c,32);
	      __assume_aligned(b,32);
	      __assume_aligned(a,32);
#elif defined __GNUC__ || !defined __INTEL_COMPILER
              c = (double*)__builtin_assume_aligned(c,32);
	      b = (float*) __builtin_assume_aligned(b,32);
	      a = (float*) __builtin_assume_aligned(a,32);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
               for(; idx != len; ++idx) {
		   _mm_prefetch((const char*)&a+32,_MM_HINT_T0);
                   aVal = _mm256_load_ps(a);
		   aVal1 = _mm256_extractf128_ps(aVal,0);
		   aDbl1 = _mm256_cvtps_pd(aVal1);
                   aVal2 = _mm256_extractf128_ps(aVal,1);
		   aDbl2 = _mm256_cvtps_pd(aVal2);
		    _mm_prefetch((const char*)&b+32,_MM_HINT_T0);    
		   bVal1 = _mm256_load_pd(b);
		   cVal1 = _mm256_mul_pd(aDbl1, bVal1);
                   bVal2 = _mm256_load_pd(b+4);
                   cVal2 = _mm256_mul_pd(aDbl2, bVal2);
		   _mm256_store_pd(c,cVal1);
                   _mm256_store_pd(c+4,cVal2);      
                   a += 8;
                   b += 8;
                   c += 8;
               }
               idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
               for(; idx != npoints; ++idx) {
                     c[i] = ((double)a[i])*b[i];
	       }        
	 }


#endif /* __GMS_32F_64F_MUL_64F_AVX_HPP__*/
