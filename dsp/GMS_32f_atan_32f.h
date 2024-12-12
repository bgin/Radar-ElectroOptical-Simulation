
#ifndef __GMS_32F_ATAN_32F_H__
#define __GMS_32F_ATAN_32F_H__ 241220200320


/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 24-12-2020 3:20PM +00200
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


#include <cstdint>
#include "GMS_config.h"



#if !defined(DSP_32F_ATAN_32F_BLOCK)
    #define DSP_32F_ATAN_32F_BLOCK                               \
       int32_t idx = 0;                                          \
       const int32_t len = npoints/8;                            \
       register __m256 aVal,pio2, x, y, z, arctangent;       \
       register __m256 n_third,p_third,condition,fzeroes,fones;	 \
       register __m256 ftwos,ffours;                             \ 
       pio2 = _mm256_set1_ps(1.5707963267948966192f);            \
       fzeroes = _mm256_setzero_ps();                            \
       fones = _mm256_set1_ps(1.0f);                             \
       ftwos = _mm256_set1_ps(2.0f);                             \
       ffours  = _mm256_set1_ps(4.0f);                           \
       n_third = _mm256_set1_ps(-0.3333333333333333333333333f);  \
       p_third = _mm256_set1_ps(0.3333333333333333333333333f);
#endif



       
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	void dsp_32f_atan_32f_u_avx_looped(float * __restrict b,
	                             float * __restrict a,
				     const int32_t npoints); 
	 }


        
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
        void dsp_32f_atan_32f_a_avx_looped(float * __restrict __ATTR_ALIGN__(64) b,
	                                      float * __restrict __ATTR_ALIGN__(64) a,
				              const int32_t npoints); 


       
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 __ATTR_VECTORCALL__
	__m256 dsp_32f_atan_32f_avx(const __m256 v); 



          





#endif /*__GMS_32F_ATAN_32F_H__*/
