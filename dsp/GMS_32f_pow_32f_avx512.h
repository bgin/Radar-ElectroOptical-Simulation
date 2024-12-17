

#ifndef __GMS_32F_POW_32F_AVX512_H__
#define __GMS_32F_POW_32F_AVX512_H__ 031020211244



/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 03-10-2021 12:44PM +00200
    contact: beniekg@gmail.com
    Few modification were added to original
    implementation (ICC pragmas, alignment directives and code layout rescheduled,
    unrolling completely 2-iteration for-loops)
    Converted to AVX512.
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



#define POW_POLY_DEGREE 3

#define POLY0_AVX512(x, c0) _mm512_set1_ps(c0)
#define POLY1_AVX512(x, c0, c1) \
    _mm256_fmadd_ps(POLY0_AVX512(x, c1), x, _mm512_set1_ps(c0))
#define POLY2_AVX512(x, c0, c1, c2) \
    _mm256_fmadd_ps(POLY1_AVX512(x, c1, c2), x, _mm512_set1_ps(c0))
#define POLY3_AVX512(x, c0, c1, c2, c3) \
    _mm256_fmadd_ps(POLY2_AVX512(x, c1, c2, c3), x, _mm512_set1_ps(c0))
#define POLY4_AVX512(x, c0, c1, c2, c3, c4) \
    _mm256_fmadd_ps(POLY3_AVX512(x, c1, c2, c3, c4), x, _mm512_set1_ps(c0))
#define POLY5_AVX512(x, c0, c1, c2, c3, c4, c5) \
    _mm256_fmadd_ps(POLY4_AVX512(x, c1, c2, c3, c4, c5), x, _mm512_set1_ps(c0))


#if !defined(DSP_32F_POW_32F_AVX512_BLOCK)
#define DSP_32F_POW_32F_AVX512_BLOCK                                                  \
        const __m512  one       = _mm512_set1_ps(1.0F);                        \
        const __m512  exp_hi    = _mm512_set1_ps(88.3762626647949F);           \
        const __m512  exp_lo    = _mm512_set1_ps(-88.3762626647949F);          \
        const __m512  ln2       = _mm512_set1_ps(0.6931471805F);               \
        const __m512  log2EF    = _mm512_set1_ps(1.44269504088896341F);        \
        const __m512  half      = _mm512_set1_ps(0.5F);                        \
        const __m512  exp_C1    = _mm512_set1_ps(0.693359375F);                \
        const __m512  exp_C2    = _mm512_set1_ps(-2.12194440e-4F);             \
        const __m512i  pi32_0x7f = _mm512_set1_epi32(0x7f);                     \
        const __m512  exp_p0    = _mm512_set1_ps(1.9875691500e-4F);            \
        const __m512  exp_p1    = _mm512_set1_ps(1.3981999507e-3F);            \
        const __m512  exp_p2    = _mm512_set1_ps(8.3334519073e-3F);            \
        const __m512  exp_p3    = _mm512_set1_ps(4.1665795894e-2F);            \ 
        const __m512  exp_p4    = _mm512_set1_ps(1.6666665459e-1F);            \
        const __m512  exp_p5    = _mm512_set1_ps(5.0000001201e-1F);            \
        __m512 aVal             = _mm512_setzero_ps();                         \
	__m512 bVal             = aVal;                                        \
	__m512 cVal             = aVal;                                        \
	__m512 logarithm        = aVal;                                        \
	__m512 mantissa         = aVal;                                        \
	__m512 frac             = aVal;                                        \
	__m512 leadingOne       = aVal;                                        \
	__m512 tmp              = aVal;                                        \
	__m512 fx               = aVal;                                        \
	__m512 pow2n            = aVal;                                        \
	__m512 z                = aVal;                                        \
	__m512 y                = aVal;                                        \
	__m512i bias;                                                          \
	__m512i exp;                                                           \
	__m512i emm0;                                                          \
	int32_t idx = 0;                                                       \
	const int32_t len = npoints/16;                                        \
	__mmask16 mask = 0;
#endif



        
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
         void dsp_32f_pow_32f_u_avx512_looped(float * __restrict c,
	                                      const float * __restrict b,
				              const float * __restrict a,
				              const int32_t npoints); 



	
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	void dsp_32f_pow_32f_a_avx512_looped(float * __restrict __ATTR_ALIGN__(64) c,
	                                      float * __restrict __ATTR_ALIGN__(64) b,
				              float * __restrict __ATTR_ALIGN__(64) a,
				              const int32_t npoints); 


	
	  __ATTR_HOT__
	  __ATTR_ALIGN__(32)
	  __ATTR_VECTORCALL__
	  __m512 dsp_32f_pow_32_avx512(const __m512 x,
	                               const __m512 y  ); 

















#endif /*__GMS_32F_POW_32F_AVX512_H__*/
