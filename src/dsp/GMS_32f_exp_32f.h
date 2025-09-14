

#ifndef __GMS_32F_EXP_32F_H__
#define __GMS_32F_EXP_32F_H__


/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 09-10-2021 13:34AM +00200
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

 /* SIMD (SSE4) implementation of exp
   Inspired by Intel Approximate Math library, and based on the
   corresponding algorithms of the cephes math library
*/

/* Copyright (C) 2007  Julien Pommier

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  (this is the zlib license)
*/

#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"



#if !defined(DSP_32F_EXP_32F_BLOCK)
#define DSP_32F_EXP_32F_BLOCK                                                \
         const __m256 one        = _mm256_set1_ps(1.0F);                     \
	 const __m256 exp_hi     = _mm256_set1_ps(88.3762626647949F);        \
	 const __m256 exp_lo     = _mm256_set1_ps(-88.3762626647949F);       \
	 const __m256 log2EF     = _mm256_set1_ps(1.44269504088896341F);     \
	 const __m256 half       = _mm256_set1_ps(0.5F);                     \ 
         const __m256 exp_C1     = _mm256_set1_ps(0.693359375F);             \
         const __m256 exp_C2     = _mm256_set1_ps(-2.12194440e-4F);          \
         const __m256i pi32_0x7f = _mm256_set1_epi32(0x7f);                  \
         const __m256  exp_p0    = _mm256_set1_ps(1.9875691500e-4F);         \
         const __m256  exp_p1    = _mm256_set1_ps(1.3981999507e-3F);         \          
         const __m256  exp_p2    = _mm256_set1_ps(8.3334519073e-3F);         \
         const __m256  exp_p3    = _mm256_set1_ps(4.1665795894e-2F);         \
         const __m256  exp_p4    = _mm256_set1_ps(1.6666665459e-1F);         \
         const __m256  exp_p5    = _mm256_set1_ps(5.0000001201e-1F);         \
	 __m256 aVal  = _mm256_setzero_ps();                                 \
	 __m256 bVal  = aVal;                                                \
	 __m256 tmp   = aVal;                                                \
	 __m256 fx    = aVal;                                                \
	 __m256 mask  = aVal;                                                \
	 __m256 pow2n = aVal;                                                \
	 __m256 z     = aVal;                                                \
	 __m256 y     = aVal;                                                \
	 __m256i emm0;                                                       \
	 int32_ idx = 0;                                                     \
	 const int32_t len = npoints / 8;
#endif


       
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
         void dsp_32f_exp_32f_u_avx_looped(float * __restrict b,
	                                   const float * __restrict a,
					   const int32_t npoints); 
					   
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	void dsp_32f_exp_32f_a_avx_looped(float * __restrict __ATTR_ALIGN__(32) b,
	                                   const float * __restrict __ATTR_ALIGN__(32) a,
					   const int32_t npoints); 
					   
	  __ATTR_HOT__
	  __ATTR_ALIGN__(32)
	  __ATTR_VECTORCALL__
	  __m256 dsp_32f_exp_32f_avx(const __m256 v); 







#endif /*__GMS_32F_EXP_32F_H__*/
