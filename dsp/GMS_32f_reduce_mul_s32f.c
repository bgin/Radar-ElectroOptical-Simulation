



/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 13-12-2020 1:36PM +00200
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



#include "GMS_32f_reduce_mul_s32f.h"


	 float dsp_32f_reduce_mul_s32f_u_avx(float * __restrict data,
	                                     const int32_t npoints) {
               DSP_32F_REDUCE_MUL_S32F_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
               for(; idx != len; ++idx) {
                   ymm0  = _mm256_loadu_ps(data);
		   accum = _mm256_mul_ps(accum,ymm0); // Potential single-precision overflow!!
		   data += 8;
	       }
	       _mm256_storeu_ps(tbuf,accum);
	       // Complete unrolling of horizontal reduction
	       result =  tbuf[0];
	       result *= tbuf[1];
	       result *= tbuf[2];
	       result *= tbuf[3];
	       result *= tbuf[4];
	       result *= tbuf[5];
	       result *= tbuf[6];
	       result *= tbuf[7];
	       idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
               for(; idx != npoints; ++idx) {
                   result *= data[i];
	       }
	       return (result);
	 }

	
	 float dsp_32f_reduce_mul_s32f_a_avx(float * __restrict __ATTR_ALIGN__(32) data,
	                                     const int32_t npoints) {
               DSP_32F_REDUCE_MUL_S32F_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
               __assume_aligned(data,32);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              data = (float*)__builtin_assume_aligned(data,32);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
              for(; idx != len; ++idx) {
                   ymm0  = _mm256_load_ps(data);
		   accum = _mm256_mul_ps(accum,ymm0); // Potential single-precision overflow!!
		   data += 8;
	       }
	       _mm256_store_ps(tbuf,accum);
	       // Complete unrolling of horizontal reduction
	       result =  tbuf[0];
	       result *= tbuf[1];
	       result *= tbuf[2];
	       result *= tbuf[3];
	       result *= tbuf[4];
	       result *= tbuf[5];
	       result *= tbuf[6];
	       result *= tbuf[7];
	       idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
               for(; idx != npoints; ++idx) {
                   result *= data[i];
	       }
	       return (result);
	 }





