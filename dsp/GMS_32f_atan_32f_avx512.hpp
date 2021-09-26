

#ifndef __GMS_32F_ATAN_32F_AVX512_HPP__
#define __GMS_32F_ATAN_32F_AVX512_HPP__ 260920211224

/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 26-09-2021 12:24PM +00200
    contact: beniekg@gmail.com
    Few modification were added to original
    implementation (ICC pragmas, alignment directives and code layout rescheduled,
    unrolling completely 2-iteration for-loops)
    Converted to AVX512 implementation
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

#if !defined(DSP_32F_ATAN_32F_AVX512_BLOCK)
    #define DSP_32F_ATAN_32F_AVX512_BLOCK                        \
       int32_t idx = 0;                                          \
       const int32_t len = npoints/16;                           \
       register __m512 aVal,pio2, x, y, z, arctangent;           \
       register __m512 n_third,p_third,fzeroes,fones;	         \
       register __m512 ftwos,ffours,t0;                          \
       __mmask16 condition = 0;                                  \
       pio2 = _mm512_set1_ps(1.5707963267948966192f);            \
       fzeroes = _mm512_setzero_ps();                            \
       t0 = fzeroes;                                             \
       fones = _mm512_set1_ps(1.0f);                             \
       ftwos = _mm512_set1_ps(2.0f);                             \
       ffours  = _mm512_set1_ps(4.0f);                           \
       n_third = _mm512_set1_ps(-0.3333333333333333333333333f);  \
       p_third = _mm512_set1_ps(0.3333333333333333333333333f);
#endif

         __ATTR_ALWAYS_INLINE__
         __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_atan_32f_u_avx512_looped(float * __restrict b,
	                                       float * __restrict a,
				               const int32_t npoints) {
                DSP_32F_ATAN_32F_AVX512_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                for(; idx != len; ++idx) {
                      _mm_prefetch((const char*)&a+32,_MM_HINT_T0);
		      aVal = _mm512_loadu_ps(a);
		      z = aVal;
		      condition = _mm512_cmp_ps_mask(z,fzeroes,_CMP_LT_OQ);
		      z = _mm512_maskz_sub_ps(condition,z,
		                                    _mm512_mul_ps(z,ftwos));
		      condition = _mm512_cmp_ps_mask(z,fones,_CMP_LT_OQ);
		      x = _mm512_maskz_add_ps(z,
		                            _mm512_sub_ps(_mm512_div_ps(fones,z),z));
		       // Original loop of 2-cycles removed
		      /*
                             for (i = 0; i < 2; i++) {
                                 x = _mm256_add_ps(x, _mm256_sqrt_ps(_mm256_fmadd_ps(x, x, fones)));
                             }
                       */
		      t0 = _mm512_sqrt_ps(_mm512_fmadd_ps(x, x, fones);
		      x  = _mm512_add_ps(x,t0);
		      x  = _mm512_add_ps(x,t0);
		      x  = _mm512_div_ps(fones,x);
		       // Original loop of 2-cycles removed
		      /*
                            for (j = TERMS - 1; j >= 0; j--) {
                                 y = _mm256_fmadd_ps(
                                 y, _mm256_mul_ps(x, x), _mm256_set1_ps(pow(-1, j) / (2 * j + 1)));
                            }
                       */
		        y = _mm512_fmadd_ps(
                             y, _mm512_mul_ps(x, x),n_third); // removing call to pow
		        y = _mm512_fmadd_ps(
                             y, _mm512_mul_ps(x, x), p_third);  // removed call to pow
		        y = _mm512_mul_ps(y, _mm512_mul_ps(x, ffours));
			condition = _mm512_cmp_ps_mask(z,fones,_CMP_GT_OQ);
			y = _mm512_maskz_add_ps(condition,y,_mm512_fnmadd_ps(y,ftwos,pio2));
			arctangent = y;
			condition = _mm512_cmp_ps_mask(aVal,fzeroes,_CMP_LT_OQ);
			arctangent = _mm512_maskz_sub_ps(condition,arctangent,
			                               _mm512_mul_ps(arctangent,ftwos));
			_mm512_storeu_ps(b,arctangent);
			a += 16;
			b += 16;
		}
		 idx = len*16;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
                for(; idx != npoints; ++idx) {
                      b[i] = ceph_atanf(a[i]);
	          }
         }


         __ATTR_ALWAYS_INLINE__
         __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_atan_32f_a_avx512_looped(float * __restrict __ATTR_ALIGN__(64) b,
	                                       float * __restrict __ATTR_ALIGN__(64) a,
				               const int32_t npoints) {
                DSP_32F_ATAN_32F_AVX512_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
              __assume_aligned(b,64);
	      __assume_aligned(a,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              b = (float*)__builtin_assume_aligned(b,64);
	      a = (float*)__builtin_assume_aligned(a,64);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                for(; idx != len; ++idx) {
                      _mm_prefetch((const char*)&a+32,_MM_HINT_T0);
		      aVal = _mm512_load_ps(a);
		      z = aVal;
		      condition = _mm512_cmp_ps_mask(z,fzeroes,_CMP_LT_OQ);
		      z = _mm512_maskz_sub_ps(condition,z,
		                                    _mm512_mul_ps(z,ftwos));
		      condition = _mm512_cmp_ps_mask(z,fones,_CMP_LT_OQ);
		      x = _mm512_maskz_add_ps(z,
		                            _mm512_sub_ps(_mm512_div_ps(fones,z),z));
		       // Original loop of 2-cycles removed
		      /*
                             for (i = 0; i < 2; i++) {
                                 x = _mm256_add_ps(x, _mm256_sqrt_ps(_mm256_fmadd_ps(x, x, fones)));
                             }
                       */
		      t0 = _mm512_sqrt_ps(_mm512_fmadd_ps(x, x, fones);
		      x  = _mm512_add_ps(x,t0);
		      x  = _mm512_add_ps(x,t0);
		      x  = _mm512_div_ps(fones,x);
		       // Original loop of 2-cycles removed
		      /*
                            for (j = TERMS - 1; j >= 0; j--) {
                                 y = _mm256_fmadd_ps(
                                 y, _mm256_mul_ps(x, x), _mm256_set1_ps(pow(-1, j) / (2 * j + 1)));
                            }
                       */
		        y = _mm512_fmadd_ps(
                             y, _mm512_mul_ps(x, x),n_third); // removing call to pow
		        y = _mm512_fmadd_ps(
                             y, _mm512_mul_ps(x, x), p_third);  // removed call to pow
		        y = _mm512_mul_ps(y, _mm512_mul_ps(x, ffours));
			condition = _mm512_cmp_ps_mask(z,fones,_CMP_GT_OQ);
			y = _mm512_maskz_add_ps(condition,y,_mm512_fnmadd_ps(y,ftwos,pio2));
			arctangent = y;
			condition = _mm512_cmp_ps_mask(aVal,fzeroes,_CMP_LT_OQ);
			arctangent = _mm512_maskz_sub_ps(condition,arctangent,
			                               _mm512_mul_ps(arctangent,ftwos));
			_mm512_store_ps(b,arctangent);
			a += 16;
			b += 16;
		}
		 idx = len*16;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
                for(; idx != npoints; ++idx) {
                      b[i] = ceph_atanf(a[i]);
	          }
          }


	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 __ATTR_VECTORCALL__
	 __attribute__((regcall)) // GCC will skip over this attribute!!
	 static inline
	 __m512 dsp_32f_atan_32f_avx512(const __m512 v) {

	             register __m512 aVal,pio2, x, y, z, arctangent;           
                     register __m512 n_third,p_third,fzeroes,fones;	         
                     register __m512 ftwos,ffours,t0;                          
                     __mmask16 condition = 0;                                  
                     pio2 = _mm512_set1_ps(1.5707963267948966192f);            
                     fzeroes = _mm512_setzero_ps();                            
                     t0 = fzeroes;                                             
                     fones = _mm512_set1_ps(1.0f);                             
                     ftwos = _mm512_set1_ps(2.0f);                             
                     ffours  = _mm512_set1_ps(4.0f);                           
                     n_third = _mm512_set1_ps(-0.3333333333333333333333333f);  
                     p_third = _mm512_set1_ps(0.3333333333333333333333333f);
		     aVal = v;
		     z = aVal;
		     condition = _mm512_cmp_ps_mask(z,fzeroes,_CMP_LT_OQ);
		     z = _mm512_maskz_sub_ps(condition,z,
		                                    _mm512_mul_ps(z,ftwos));
		     condition = _mm512_cmp_ps_mask(z,fones,_CMP_LT_OQ);
		     x = _mm512_maskz_add_ps(z,
		                            _mm512_sub_ps(_mm512_div_ps(fones,z),z));
		       // Original loop of 2-cycles removed
		      /*
                             for (i = 0; i < 2; i++) {
                                 x = _mm256_add_ps(x, _mm256_sqrt_ps(_mm256_fmadd_ps(x, x, fones)));
                             }
                       */
		     t0 = _mm512_sqrt_ps(_mm512_fmadd_ps(x, x, fones);
		     x  = _mm512_add_ps(x,t0);
		     x  = _mm512_add_ps(x,t0);
		     x  = _mm512_div_ps(fones,x);
		       // Original loop of 2-cycles removed
		      /*
                            for (j = TERMS - 1; j >= 0; j--) {
                                 y = _mm256_fmadd_ps(
                                 y, _mm256_mul_ps(x, x), _mm256_set1_ps(pow(-1, j) / (2 * j + 1)));
                            }
                       */
		     y = _mm512_fmadd_ps(
                             y, _mm512_mul_ps(x, x),n_third); // removing call to pow
		     y = _mm512_fmadd_ps(
                             y, _mm512_mul_ps(x, x), p_third);  // removed call to pow
		     y = _mm512_mul_ps(y, _mm512_mul_ps(x, ffours));
		     condition = _mm512_cmp_ps_mask(z,fones,_CMP_GT_OQ);
	             y = _mm512_maskz_add_ps(condition,y,_mm512_fnmadd_ps(y,ftwos,pio2));
		     arctangent = y;
		     condition = _mm512_cmp_ps_mask(aVal,fzeroes,_CMP_LT_OQ);
		     arctangent = _mm512_maskz_sub_ps(condition,arctangent,
			                               _mm512_mul_ps(arctangent,ftwos));
		     return (arctangent);
	}










#endif /*__GMS_32F_ATAN_32F_AVX512_HPP__*/
