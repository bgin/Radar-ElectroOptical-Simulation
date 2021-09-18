
#ifndef __GMS_32F_ASIN_32F_HPP__
#define __GMS_32F_ASIN_32F_HPP__


/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 24-12-2020 11:30AM +00200
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
#include "GMS_cephes.h"



#if !defined(DSP_32F_ASIN_32F_BLOCK)
    #define DSP_32F_ASIN_32F_BLOCK                                   \
        int32_t idx = 0;                                             \
	const int32_t len = npoints/8;                               \
	register __m256 aVal, pio2, x, y, z, arcsine;                \
        register __m256 fzeroes, fones, ftwos, ffours, condition;    \
        register __m256 n_third,p_third;                             \
	pio2 = _mm256_set1_ps(1.5707963267948966192f);               \
	fzeroes = _mm256_setzero_ps();                               \
        fones = _mm256_set1_ps(1.0f);                                \
        ftwos = _mm256_set1_ps(2.0f);                                \
        ffours  = _mm256_set1_ps(4.0f);                              \
        n_third = _mm256_set1_ps(-0.3333333333333333333333333f);     \
        p_third = _mm256_set1_ps(0.3333333333333333333333333f);
#endif
          __ATTR_ALWAYS_INLINE__
          __ATTR_HOT__
          __ATTR_ALIGN__(32)
	  static inline
	  void dsp_32f_asin_32f_u_avx(float * __restrict b,
	                              float * __restrict a,
				      const int32_t npoints) {
             DSP_32F_ASIN_32F_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
              for(; idx != len; ++idx) {
                  aVal = _mm256_loadu_ps(a);
                  aVal = _mm256_div_ps(aVal,
                           _mm256_sqrt_ps(_mm256_mul_ps(_mm256_add_ps(fones, aVal),
                                                          _mm256_sub_ps(fones, aVal))));
                  z = aVal;
                  condition = _mm256_cmp_ps(z, fzeroes, _CMP_LT_OQ);
                  z = _mm256_sub_ps(z, _mm256_and_ps(_mm256_mul_ps(z, ftwos), condition));
                  condition = _mm256_cmp_ps(z, fones, _CMP_LT_OQ);
                  x = _mm256_add_ps(
                           z, _mm256_and_ps(_mm256_sub_ps(_mm256_div_ps(fones, z), z), condition));
		  // Original code contained here a 2-cycle loop
		  /*
                        for (i = 0; i < 2; i++) {
                                x = _mm256_add_ps(x, _mm256_sqrt_ps(_mm256_fmadd_ps(x, x, fones)));
                        }
                  */
		  x = _mm256_add_ps(x, _mm256_sqrt_ps(_mm256_fmadd_ps(x, x, fones)));
		  x = _mm256_add_ps(x, _mm256_sqrt_ps(_mm256_fmadd_ps(x, x, fones)));
		  x = _mm256_div_ps(fones, x);
                  y = fzeroes;
		   // Original code contained here a 2-cycle loop
		   /*
                          for (j = ASIN_TERMS - 1; j >= 0; j--) {
                               y = _mm256_fmadd_ps(
                               y, _mm256_mul_ps(x, x), _mm256_set1_ps(pow(-1, j) / (2 * j + 1)));
                          }
                    */
		  y = _mm256_fmadd_ps(
                         y, _mm256_mul_ps(x, x),n_third); // removing call to pow
		  y = _mm256_fmadd_ps(
                         y, _mm256_mul_ps(x, x), p_third);  // removed call to pow
		  y = _mm256_mul_ps(y, _mm256_mul_ps(x, ffours));
                      condition = _mm256_cmp_ps(z, fones, _CMP_GT_OQ);
                  y = _mm256_add_ps(y, _mm256_and_ps(_mm256_fnmadd_ps(y, ftwos, pio2), condition));
                  arcsine = y;
                  condition = _mm256_cmp_ps(aVal, fzeroes, _CMP_LT_OQ);
                  arcsine = _mm256_sub_ps(arcsine,
                                _mm256_and_ps(_mm256_mul_ps(arcsine, ftwos), condition));
                  _mm256_storeu_ps(b,arcsine);
                  a += 8;
                  b += 8;
	      }
              idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
                   for(; idx != npoints; ++idx) {
                        b[i] = ceph_asinf(a[i]);
	           }	      
	  }

	  
          __ATTR_ALWAYS_INLINE__
          __ATTR_HOT__
          __ATTR_ALIGN__(32)
	  static inline
	  void dsp_32f_asin_32f_a_avx(float * __restrict b,
	                              float * __restrict a,
				      const int32_t npoints) {
             DSP_32F_ASIN_32F_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
              __assume_aligned(b,32);
	      __assume_aligned(a,32);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              b = (float*)__builtin_assume_aligned(b,32);
	      a = (float*)__builtin_assume_aligned(a,32);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
              for(; idx != len; ++idx) {
                  aVal = _mm256_load_ps(a);
                  aVal = _mm256_div_ps(aVal,
                           _mm256_sqrt_ps(_mm256_mul_ps(_mm256_add_ps(fones, aVal),
                                                          _mm256_sub_ps(fones, aVal))));
                  z = aVal;
                  condition = _mm256_cmp_ps(z, fzeroes, _CMP_LT_OQ);
                  z = _mm256_sub_ps(z, _mm256_and_ps(_mm256_mul_ps(z, ftwos), condition));
                  condition = _mm256_cmp_ps(z, fones, _CMP_LT_OQ);
                  x = _mm256_add_ps(
                           z, _mm256_and_ps(_mm256_sub_ps(_mm256_div_ps(fones, z), z), condition));
		  // Original code contained here a 2-cycle loop
		  /*
                        for (i = 0; i < 2; i++) {
                                x = _mm256_add_ps(x, _mm256_sqrt_ps(_mm256_fmadd_ps(x, x, fones)));
                        }
                  */
		  x = _mm256_add_ps(x, _mm256_sqrt_ps(_mm256_fmadd_ps(x, x, fones)));
		  x = _mm256_add_ps(x, _mm256_sqrt_ps(_mm256_fmadd_ps(x, x, fones)));
		  x = _mm256_div_ps(fones, x);
                  y = fzeroes;
		   // Original code contained here a 2-cycle loop
		   /*
                          for (j = ASIN_TERMS - 1; j >= 0; j--) {
                               y = _mm256_fmadd_ps(
                               y, _mm256_mul_ps(x, x), _mm256_set1_ps(pow(-1, j) / (2 * j + 1)));
                          }
                    */
		  y = _mm256_fmadd_ps(
                         y, _mm256_mul_ps(x, x),n_third); // removing call to pow
		  y = _mm256_fmadd_ps(
                         y, _mm256_mul_ps(x, x), p_third);  // removed call to pow
		  y = _mm256_mul_ps(y, _mm256_mul_ps(x, ffours));
                      condition = _mm256_cmp_ps(z, fones, _CMP_GT_OQ);
                  y = _mm256_add_ps(y, _mm256_and_ps(_mm256_fnmadd_ps(y, ftwos, pio2), condition));
                  arcsine = y;
                  condition = _mm256_cmp_ps(aVal, fzeroes, _CMP_LT_OQ);
                  arcsine = _mm256_sub_ps(arcsine,
                                _mm256_and_ps(_mm256_mul_ps(arcsine, ftwos), condition));
                  _mm256_store_ps(b,arcsine);
                  a += 8;
                  b += 8;
	      }
              idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
                   for(; idx != npoints; ++idx) {
                        b[i] = ceph_asinf(a[i]);
	           }	      
	  }





#endif /*__GMS_32F_ASIN_32F_HPP__*/
