

#ifndef __GMS_32F_POW_32F_H__
#define __GMS_32F_POW_32F_H__ 031020211027


/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 03-10-2021 10:27AM +00200
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


#define POW_POLY_DEGREE 3

#define POLY0_AVX2(x, c0) _mm256_set1_ps(c0)
#define POLY1_AVX2(x, c0, c1) \
    _mm256_fmadd_ps(POLY0_AVX2(x, c1), x, _mm256_set1_ps(c0))
#define POLY2_AVX2(x, c0, c1, c2) \
    _mm256_fmadd_ps(POLY1_AVX2(x, c1, c2), x, _mm256_set1_ps(c0))
#define POLY3_AVX2(x, c0, c1, c2, c3) \
    _mm256_fmadd_ps(POLY2_AVX2(x, c1, c2, c3), x, _mm256_set1_ps(c0))
#define POLY4_AVX2(x, c0, c1, c2, c3, c4) \
    _mm256_fmadd_ps(POLY3_AVX2(x, c1, c2, c3, c4), x, _mm256_set1_ps(c0))
#define POLY5_AVX2(x, c0, c1, c2, c3, c4, c5) \
    _mm256_fmadd_ps(POLY4_AVX2(x, c1, c2, c3, c4, c5), x, _mm256_set1_ps(c0))


#if !defined(DSP_32F_POW_32F_BLOCK)
#define DSP_32F_POW_32F_BLOCK                                                  \
        const __m256  one       = _mm256_set1_ps(1.0F);                        \
        const __m256  exp_hi    = _mm256_set1_ps(88.3762626647949F);           \
        const __m256  exp_lo    = _mm256_set1_ps(-88.3762626647949F);          \
        const __m256  ln2       = _mm256_set1_ps(0.6931471805F);               \
        const __m256  log2EF    = _mm256_set1_ps(1.44269504088896341F);        \
        const __m256  half      = _mm256_set1_ps(0.5F);                        \
        const __m256  exp_C1    = _mm256_set1_ps(0.693359375F);                \
        const __m256  exp_C2    = _mm256_set1_ps(-2.12194440e-4F);             \
        const __m26i  pi32_0x7f = _mm256_set1_epi32(0x7f);                     \
        const __m256  exp_p0    = _mm256_set1_ps(1.9875691500e-4F);            \
        const __m256  exp_p1    = _mm256_set1_ps(1.3981999507e-3F);            \
        const __m256  exp_p2    = _mm256_set1_ps(8.3334519073e-3F);            \
        const __m256  exp_p3    = _mm256_set1_ps(4.1665795894e-2F);            \ 
        const __m256  exp_p4    = _mm256_set1_ps(1.6666665459e-1F);            \
        const __m256  exp_p5    = _mm256_set1_ps(5.0000001201e-1F);            \
        __m256 aVal             = _mm256_setzero_ps();                         \
	__m256 bVal             = aVal;                                        \
	__m256 cVal             = aVal;                                        \
	__m256 logarithm        = aVal;                                        \
	__m256 mantissa         = aVal;                                        \
	__m256 frac             = aVal;                                        \
	__m256 leadingOne       = aVal;                                        \
	__m256 tmp              = aVal;                                        \
	__m256 fx               = aVal;                                        \
	__m256 mask             = aVal;                                        \
	__m256 pow2n            = aVal;                                        \
	__m256 z                = aVal;                                        \
	__m256 y                = aVal;                                        \
	__m256i bias;                                                          \
	__m256i exp;                                                           \
	__m256i emm0;                                                          \
	int32_t idx = 0;                                                       \
	const int32_t len = npoints/8;
#endif




      
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	void dsp_32f_pow_32f_u_avx2_looped(float * __restrict c,
	                                    const float * __restrict b,
				            const float * __restrict a,
				            const int32_t npoints); 


	
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	void dsp_32f_pow_32f_a_avx2_looped(float * __restrict __ATTR_ALIGN__(32) c,
	                                    float * __restrict __ATTR_ALIGN__(32) b,
				            float * __restrict __ATTR_ALIGN__(32) a,
				            const int32_t npoints); 


	
	  __ATTR_HOT__
	  __ATTR_ALIGN__(32)
	  __ATTR_VECTORCALL__
	  __m256 dsp_32f_pow_32_avx2(const __m256 x,
	                             const __m256 y  ); 


#endif /*__GMS_32F_POW_32F_H__*/
