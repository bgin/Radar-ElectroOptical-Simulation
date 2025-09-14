
#ifndef __GMS_32F_TAN_32F_H__
#define __GMS_32F_TAN_32F_H__ 091020211027


/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 09-10-2021 10:27AM +00200
    contact: beniekg@gmail.com
    Few modification were added to original
    implementation (ICC pragmas, alignment directives and code layout rescheduled,
    unrolling completely 2-iteration for-loops)
    
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



#if !defined(DSP_32F_TAN_32F_BLOCK)
#define DSP_32F_TAN_32F_BLOCK                                        \
        const __m256  m4pi    = _mm256_set1_ps(1.273239545F);        \
        const __m256  pio4A   = _mm256_set1_ps(0.78515625F);          \
        const __m256  pio4B   = _mm256_set1_ps(0.241876e-3F);         \
	const __m256  feight  = _mm256_set1_ps(0.125F);               \
        const __m256  ffours  = _mm256_set1_ps(4.0F);                 \
        const __m256  ftwos   = _mm256_set1_ps(2.0F);                 \
	const __m256  fones   = _mm256_set1_ps(1.0F);                 \
	const __m256  fzeroes = _mm256_setzero_ps();                  \
	const __m256  cp1     = fones;                                \
        const __m256  cp2     = _mm256_set1_ps(0.83333333e-1F);       \
        const __m256  cp3     = _mm256_set1_ps(0.2777778e-2F);        \
        const __m256  cp4     = _mm256_set1_ps(0.49603e-4F);          \
        const __m256  cp5     = _mm256_set1_ps(0.551e-6F);            \
	const __m256i ones    = _mm256_set1_epi32(1);                 \
	const __m256i twos    = _mm256_set1_epi32(2);                 \
	const __m256  fours   = _mm256_set1_epi32(4);                 \
	__m256 aVal           = fzeroes;                              \
	__m256 s              = fzeroes;                              \
	__m256 sine           = fzeroes;                              \
	__m256 cosine         = fzeroes;                              \
	__m256 tangent        = fzeroes;                              \
	__m256 cond1          = fzeroes;                              \
	__m256 cond2          = fzeroes;                              \
	__m256 cond3          = fzeroes;                              \
	__m256i q;                                                    \
	__m256i r;                                                    \
	int32_t i = 0;                                                \
	const int32_t len = npoints/8;
#endif



       
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	void dsp_32f_tan_32f_u_avx2_looped(float * __restrict b,
	                                    const float * __restrict a,
					    const int32_t npoints); 

	  
     
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	void dsp_32f_tan_32f_a_avx2_looped(float * __restrict __ATTR_ALIGN__(32) b,
	                                    float * __restrict __ATTR_ALIGN__(32) a,
					    const int32_t npoints); 


	  
        
	  __ATTR_HOT__
	  __ATTR_ALIGN__(32)
	  __ATTR_VECTORCALL__
	  __m256 dsp_32f_tan_32f_avx2(__m256 v); 









#endif /*__GMS_32F_TAN_32F_H__*/
