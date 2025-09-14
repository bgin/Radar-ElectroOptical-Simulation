
#ifndef __GMS_32F_TAN_32F_HPP__
#define __GMS_32F_TAN_32F_HPP__ 091020211027


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
#include "GMS_cephes.h"


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



         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_tan_32f_u_avx2_looped(float * __restrict b,
	                                    const float * __restrict a,
					    const int32_t npoints) {
               DSP_32F_TAN_32F_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                 for(; idx != len; ++idx) {
                       _mm_prefetch((const char*)&a+32,_MM_HINT_T0);
		       aVal = _mm256_loadu_ps(a);
		       s = _mm256_sub_ps(aVal,
                                   _mm256_and_ps(_mm256_mul_ps(aVal, ftwos),
                                                  _mm256_cmp_ps(aVal, fzeroes, _CMP_LT_OS)));
                       q = _mm256_cvtps_epi32(_mm256_floor_ps(_mm256_mul_ps(s, m4pi)));
                       r = _mm256_add_epi32(q, _mm256_and_si256(q, ones));
                       s = _mm256_fnmadd_ps(_mm256_cvtepi32_ps(r), pio4A, s);
                       s = _mm256_fnmadd_ps(_mm256_cvtepi32_ps(r), pio4B, s);
                       s = _mm256_mul_ps(s,feight);
		       s = _mm256_mul_ps(s,s);
		       s = _mm256_mul_ps(
                                     _mm256_fmadd_ps(
                                              _mm256_fmsub_ps(
                                                  _mm256_fmadd_ps(
						           _mm256_fmsub_ps(s, cp5, cp4), s, cp3), s, cp2),
                                                  s,
                                             cp1),
                                     s);
		       const __m256 t0 = _mm256_sub_ps(ffours,s);
		       s = _mm256_mul_ps(s,t0);
		       s = _mm256_mul_ps(s,t0);
		       s = _mm256_mul_ps(s,t0);
		       sine = _mm256_sqrt_ps(_mm256_mul_ps(
		                                _mm256_sub_ps(ftwos, s), s));
                       cosine = _mm256_sub_ps(fones, s);
		       cond1  = _mm256_cmp_ps(
                                  _mm256_cvtepi32_ps(_mm256_and_si256(
				                      _mm256_add_epi32(q, ones), twos)),
                                                                                   fzeroes,
                                                                                      _CMP_NEQ_UQ);
                       cond2 = _mm256_cmp_ps(
                                    _mm256_cmp_ps(
                                            _mm256_cvtepi32_ps(_mm256_and_si256(q, fours)),
					                                  fzeroes, _CMP_NEQ_UQ),
                                                                           _mm256_cmp_ps(aVal, fzeroes, _CMP_LT_OS),
                                                                                                             _CMP_NEQ_UQ);
                       cond3 = _mm256_cmp_ps(
                                     _mm256_cvtepi32_ps(_mm256_and_si256(
				                           _mm256_add_epi32(q, twos), fours)),
                                                                                        fzeroes,
                                                                                            _CMP_NEQ_UQ);
		      __m256 temp = cosine;
                      cosine =
                              _mm256_add_ps(cosine, _mm256_and_ps(
			                                 _mm256_sub_ps(sine, cosine), cond1));
                      sine = _mm256_add_ps(sine, _mm256_and_ps(
		                                         _mm256_sub_ps(temp, sine), cond1));
                      sine = _mm256_sub_ps(
                                     sine, _mm256_and_ps(
				                         _mm256_mul_ps(sine, ftwos), cond2));
                      cosine = _mm256_sub_ps(
                                             cosine,
                                                   _mm256_and_ps(_mm256_mul_ps(cosine, ftwos), cond3));
                      tangent = _mm256_div_ps(sine, cosine);
		      _mm256_storeu_ps(b,tangent);
		      a += 8;
		      b += 8;
		}
		idx = len * 8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
              for(; idx != npoints; ++idx) {
                  b[i] = ceph_tanf(a[i]);
	      }			
	 }


	  
         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_tan_32f_a_avx2_looped(float * __restrict __ATTR_ALIGN__(32) b,
	                                    float * __restrict __ATTR_ALIGN__(32) a,
					    const int32_t npoints) {
               DSP_32F_TAN_32F_BLOCK
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
                       _mm_prefetch((const char*)&a+32,_MM_HINT_T0);
		       aVal = _mm256_load_ps(a);
		       s = _mm256_sub_ps(aVal,
                                   _mm256_and_ps(_mm256_mul_ps(aVal, ftwos),
                                                  _mm256_cmp_ps(aVal, fzeroes, _CMP_LT_OS)));
                       q = _mm256_cvtps_epi32(_mm256_floor_ps(_mm256_mul_ps(s, m4pi)));
                       r = _mm256_add_epi32(q, _mm256_and_si256(q, ones));
                       s = _mm256_fnmadd_ps(_mm256_cvtepi32_ps(r), pio4A, s);
                       s = _mm256_fnmadd_ps(_mm256_cvtepi32_ps(r), pio4B, s);
                       s = _mm256_mul_ps(s,feight);
		       s = _mm256_mul_ps(s,s);
		       s = _mm256_mul_ps(
                                     _mm256_fmadd_ps(
                                              _mm256_fmsub_ps(
                                                  _mm256_fmadd_ps(
						           _mm256_fmsub_ps(s, cp5, cp4), s, cp3), s, cp2),
                                                  s,
                                             cp1),
                                     s);
		       const __m256 t0 = _mm256_sub_ps(ffours,s);
		       s = _mm256_mul_ps(s,t0);
		       s = _mm256_mul_ps(s,t0);
		       s = _mm256_mul_ps(s,t0);
		       sine = _mm256_sqrt_ps(_mm256_mul_ps(
		                                _mm256_sub_ps(ftwos, s), s));
                       cosine = _mm256_sub_ps(fones, s);
		       cond1  = _mm256_cmp_ps(
                                  _mm256_cvtepi32_ps(_mm256_and_si256(
				                      _mm256_add_epi32(q, ones), twos)),
                                                                                   fzeroes,
                                                                                      _CMP_NEQ_UQ);
                       cond2 = _mm256_cmp_ps(
                                    _mm256_cmp_ps(
                                            _mm256_cvtepi32_ps(_mm256_and_si256(q, fours)),
					                                  fzeroes, _CMP_NEQ_UQ),
                                                                           _mm256_cmp_ps(aVal, fzeroes, _CMP_LT_OS),
                                                                                                             _CMP_NEQ_UQ);
                       cond3 = _mm256_cmp_ps(
                                     _mm256_cvtepi32_ps(_mm256_and_si256(
				                           _mm256_add_epi32(q, twos), fours)),
                                                                                        fzeroes,
                                                                                            _CMP_NEQ_UQ);
		      __m256 temp = cosine;
                      cosine =
                              _mm256_add_ps(cosine, _mm256_and_ps(
			                                 _mm256_sub_ps(sine, cosine), cond1));
                      sine = _mm256_add_ps(sine, _mm256_and_ps(
		                                         _mm256_sub_ps(temp, sine), cond1));
                      sine = _mm256_sub_ps(
                                     sine, _mm256_and_ps(
				                         _mm256_mul_ps(sine, ftwos), cond2));
                      cosine = _mm256_sub_ps(
                                             cosine,
                                                   _mm256_and_ps(_mm256_mul_ps(cosine, ftwos), cond3));
                      tangent = _mm256_div_ps(sine, cosine);
		      _mm256_store_ps(b,tangent);
		      a += 8;
		      b += 8;
		}
		idx = len * 8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
              for(; idx != npoints; ++idx) {
                  b[i] = ceph_tanf(a[i]);
	      }			
	 }



	  
          __ATTR_ALWAYS_INLINE__
	  __ATTR_HOT__
	  __ATTR_ALIGN__(32)
	  __ATTR_VECTORCALL__
	  __attribute__((regcall)) // GCC will skip over this attribute!!
	  static inline
          __m256 dsp_32f_tan_32f_avx2(__m256 v) {

                    const __m256  m4pi    = _mm256_set1_ps(1.273239545F);       
                    const __m256  pio4A   = _mm256_set1_ps(0.78515625F);         
                    const __m256  pio4B   = _mm256_set1_ps(0.241876e-3F);       
	            const __m256  feight  = _mm256_set1_ps(0.125F);              
                    const __m256  ffours  = _mm256_set1_ps(4.0F);               
                    const __m256  ftwos   = _mm256_set1_ps(2.0F);               
	            const __m256  fones   = _mm256_set1_ps(1.0F);               
	            const __m256  fzeroes = _mm256_setzero_ps();                
	            const __m256  cp1     = fones;                              
                    const __m256  cp2     = _mm256_set1_ps(0.83333333e-1F);      
                    const __m256  cp3     = _mm256_set1_ps(0.2777778e-2F);       
                    const __m256  cp4     = _mm256_set1_ps(0.49603e-4F);         
                    const __m256  cp5     = _mm256_set1_ps(0.551e-6F);           
	            const __m256i ones    = _mm256_set1_epi32(1);                
	            const __m256i twos    = _mm256_set1_epi32(2);               
	            const __m256  fours   = _mm256_set1_epi32(4);                
	            __m256 aVal           = fzeroes;                            
	            __m256 s              = fzeroes;                             
	            __m256 sine           = fzeroes;                             
	            __m256 cosine         = fzeroes;                             
	            __m256 tangent        = fzeroes;                             
	            __m256 cond1          = fzeroes;                             
	            __m256 cond2          = fzeroes;                             
	            __m256 cond3          = fzeroes;                             
	            __m256i q;                                                  
	            __m256i r;
		    aVal = v;
		    s = _mm256_sub_ps(aVal,
                                   _mm256_and_ps(_mm256_mul_ps(aVal, ftwos),
                                                  _mm256_cmp_ps(aVal, fzeroes, _CMP_LT_OS)));
                    q = _mm256_cvtps_epi32(_mm256_floor_ps(_mm256_mul_ps(s, m4pi)));
                    r = _mm256_add_epi32(q, _mm256_and_si256(q, ones));
                    s = _mm256_fnmadd_ps(_mm256_cvtepi32_ps(r), pio4A, s);
                    s = _mm256_fnmadd_ps(_mm256_cvtepi32_ps(r), pio4B, s);
                    s = _mm256_mul_ps(s,feight);
		    s = _mm256_mul_ps(s,s);
		    s = _mm256_mul_ps(
                                     _mm256_fmadd_ps(
                                              _mm256_fmsub_ps(
                                                  _mm256_fmadd_ps(
						           _mm256_fmsub_ps(s, cp5, cp4), s, cp3), s, cp2),
                                                  s,
                                             cp1),
                                     s);
		     const __m256 t0 = _mm256_sub_ps(ffours,s);
		     s = _mm256_mul_ps(s,t0);
		     s = _mm256_mul_ps(s,t0);
		     s = _mm256_mul_ps(s,t0);
		     sine = _mm256_sqrt_ps(_mm256_mul_ps(
		                                _mm256_sub_ps(ftwos, s), s));
                     cosine = _mm256_sub_ps(fones, s);
		     cond1  = _mm256_cmp_ps(
                                  _mm256_cvtepi32_ps(_mm256_and_si256(
				                      _mm256_add_epi32(q, ones), twos)),
                                                                                   fzeroes,
                                                                                      _CMP_NEQ_UQ);
                     cond2 = _mm256_cmp_ps(
                                    _mm256_cmp_ps(
                                            _mm256_cvtepi32_ps(_mm256_and_si256(q, fours)),
					                                  fzeroes, _CMP_NEQ_UQ),
                                                                           _mm256_cmp_ps(aVal, fzeroes, _CMP_LT_OS),
                                                                                                             _CMP_NEQ_UQ);
                     cond3 = _mm256_cmp_ps(
                                     _mm256_cvtepi32_ps(_mm256_and_si256(
				                           _mm256_add_epi32(q, twos), fours)),
                                                                                        fzeroes,
                                                                                            _CMP_NEQ_UQ);
		     __m256 temp = cosine;
                     cosine =
                              _mm256_add_ps(cosine, _mm256_and_ps(
			                                 _mm256_sub_ps(sine, cosine), cond1));
                     sine = _mm256_add_ps(sine, _mm256_and_ps(
		                                         _mm256_sub_ps(temp, sine), cond1));
                     sine = _mm256_sub_ps(
                                     sine, _mm256_and_ps(
				                         _mm256_mul_ps(sine, ftwos), cond2));
                     cosine = _mm256_sub_ps(
                                             cosine,
                                                   _mm256_and_ps(_mm256_mul_ps(cosine, ftwos), cond3));
                     tangent = _mm256_div_ps(sine, cosine);
		     return (tangent);
	  }











#endif /*__GMS_32F_TAN_32F_HPP__*/
