




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
#include "GMS_32f_64f_mul_64f_avx.h"



 
        
	
	 void dsp_32f_64f_mul_64f_u_avx(double * __restrict c,
	                                const float  * __restrict a,
					                const double  * __restrict b,
					                const int32_t npoints) {
              DSP_32F_64F_MUL_64F_AVX_BLOCK

              double * __restrict       ptr_c = c;
              const float * __restrict  ptr_a = a;
              const double * __restrict ptr_b = b;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif 
              for(; idx != len; ++idx) {
                   aVal  = _mm256_loadu_ps(ptr_a);
		           aVal1 = _mm256_extractf128_ps(aVal, 0);
		           aDbl1 = _mm256_cvtps_pd(aVal1);
                   aVal2 = _mm256_extractf128_ps(aVal, 1);
		           aDbl2 = _mm256_cvtps_pd(aVal2);
		           bVal1 = _mm256_loadu_pd(ptr_b);
		           cVal1 = _mm256_mul_pd(aDbl1, bVal1);
                   bVal2 = _mm256_loadu_pd(ptr_b+4);
                   cVal2 = _mm256_mul_pd(aDbl2, bVal2);
		           _mm256_storeu_pd(ptr_c,cVal1);
                   _mm256_storeu_pd(ptr_c+4,cVal2);      
                   ptr_a += 8;
                   ptr_b += 8;
                   ptr_c += 8;
	         }
             idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
               for(; idx != npoints; ++idx) {
                     ptr_c[idx] = ((double)ptr_a[idx])*ptr_b[idx];
	       }	      
	 }

	 
       
	 void dsp_32f_64f_mul_64f_a_avx(double * __restrict __ATTR_ALIGN__(32) c,
	                                const float  * __restrict __ATTR_ALIGN__(32) a,
					                const double  * __restrict __ATTR_ALIGN__(32) b,
					                const int32_t npoints) {
              DSP_32F_64F_MUL_64F_AVX_BLOCK

              double * __restrict       ptr_c = c;
              const float * __restrict  ptr_a = a;
              const double * __restrict ptr_b = b;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
               for(; idx != len; ++idx) {
                   aVal = _mm256_load_ps(ptr_a);
		           aVal1 = _mm256_extractf128_ps(aVal,0);
		           aDbl1 = _mm256_cvtps_pd(aVal1);
                   aVal2 = _mm256_extractf128_ps(aVal,1);
		           aDbl2 = _mm256_cvtps_pd(aVal2);
		           bVal1 = _mm256_load_pd(ptr_b);
		           cVal1 = _mm256_mul_pd(aDbl1, bVal1);
                   bVal2 = _mm256_load_pd(ptr_b+4);
                   cVal2 = _mm256_mul_pd(aDbl2, bVal2);
		           _mm256_store_pd(ptr_c,cVal1);
                   _mm256_store_pd(ptr_c+4,cVal2);      
                   ptr_a += 8;
                   ptr_b += 8;
                   ptr_c += 8;
               }
               idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
               for(; idx != npoints; ++idx) {
                     ptr_c[idx] = ((double)ptr_a[idx])*ptr_b[idx];
	       }        
	 }



