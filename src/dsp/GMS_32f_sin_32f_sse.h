

#ifndef __GMS_32F_SIN_32F_SSE_H__
#define __GMS_32F_SIN_32F_SSE_H__ 300820240744


/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 30-08-2024 07:44PM +00200
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


#if !defined(__FMA__)
#error "Required support of FMA ISA!!"
#endif


#if !defined(DSP_32F_SIN_32F_SSE_BLOCK)
#define DSP_32F_SIN_32F_SSE_BLOCK                                             \
         const __m128 m4pio   =  _mm_set1_ps(1.273239545F);            \
	 const __m128 pio4A   =  _mm_set1_ps(0.78515625F);             \
	 const __m128 pio4B   =  _mm_set1_ps(0.241876e-3F);            \
	 const __m128 finv8   =  _mm_set1_ps(0.125F);                  \
	 const __m128 ffours  = _mm_set1_ps(4.0F);                     \
	 const __m128 fftwos  = _mm_set1_ps(2.0F);                     \
	 const __m128 fones   = _mm_set1_ps(1.0F);                     \
	 const __m128 fzeroes = _mm_setzero_ps();                      \
	 const __m128 finv2   = _mm_set1_ps(0.5f);                     \
	 const __m128 cp1     = _mm_set1_ps(1.0F);                     \
	 const __m128 cp2     = _mm_set1_ps(0.83333333e-1F);           \
	 const __m128 cp3     = _mm_set1_ps(0.2777778e-2F);            \
	 const __m128 cp4     = _mm_set1_ps(0.49603e-4F);              \
         const __m128 cp5     = _mm_set1_ps(0.551e-6F);                \
	 const __m128i ones   = _mm_set1_epi32(1);                     \
	 const __m128i twos   = _mm_set1_epi32(2);                     \
	 const __m128i fours  = _mm_set1_epi32(4);                     \
	 __m128 sine = fzeroes;                                           \
	 __m128 s    = fzeroes;                                           \
	 __m128 cosine = fzeroes;                                         \
	 __m128 t0     = fzeroes;                                         \
	 __m128 condition1;                                               \
	 __m128 condition2;                                               \
	 __m128i q;                                                       \
	 __m128i r;                                                       \
	 int32_t idx = 0;                                                 \
	 int32_t len = npoints/4;
#endif


      
	 void dsp_32f_sin_32f_sse_u_looped(float * __restrict b,
	                                   const float * __restrict a,
					   const int32_t npoints); 




	 void dsp_32f_sin_32f_sse_a_looped(float * __restrict ___ATTR_ALIGN__(16) b,
	                                   const float * __restrict __ATTR_ALIGN__(16) a,
					   const int32_t npoints); 


          __m128 dsp_32f_sin_32_sse(const __m128 v); 

	  







#endif /*__GMS_32F_SIN_32F_SSE_H__*/
