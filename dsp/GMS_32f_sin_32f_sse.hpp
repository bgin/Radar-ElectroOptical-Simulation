

#ifndef __GMS_32F_SIN_32F_SSE_HPP__
#define __GMS_32F_SIN_32F_SSE_HPP__ 300820240744


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
#include "GMS_cephes.h"

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


         __ATTR_ALWAYS_INLINE__
	 static inline
	 void dsp_32f_sin_32f_sse_u_looped(float * __restrict b,
	                                   const float * __restrict a,
					   const int32_t npoints) {
              DSP_32F_SIN_32F_SSE_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                  for(; idx != len; ++idx) {
                       _mm_prefetch((const char*)&a+32,_MM_HINT_T0);
		       aVal = _mm_loadu_ps(a);
		       s    = _mm_sub_ps(aVal,
                                    _mm_and_ps(_mm_mul_ps(aVal, ftwos),
                                                  _mm_cmp_ps(aVal, fzeroes, _CMP_LT_OQ)));
                       q = _mm_cvtps_epi32(_mm_floor_ps(_mm_mul_ps(s, m4pi)));
                       r = _mm_add_epi32(q, _mm_and_si128(q, ones));
                       s = _mm_fnmadd_ps(_mm_cvtepi32_ps(r), pio4A, s);
                       s = _mm_fnmadd_ps(_mm_cvtepi32_ps(r), pio4B, s);
		       s = _mm_mul_ps(s,finv8);
		       s = _mm_mul_ps(s,s);
		       s = _mm_mul_ps(
                                    _mm_fmadd_ps(
                                             _mm_fmsub_ps(
                                                  _mm_fmadd_ps(
						           _mm_fmsub_ps(s, cp5, cp4), s, cp3), s, cp2),
                                                  s,
                                            cp1),
                                     s);
		        t0 = _mm_sub_ps(ffours,s);
			s  = _mm_mul_ps(s,t0);
			s  = _mm_mul_ps(s,t0);
			s  = _mm_mul_ps(s,t0);
			s  = _mm_mul_ps(s,finv2);
			sine = _mm_sqrt_ps(
			               _mm_mul_ps(_mm_sub_ps(ftwos, s), s));
                        cosine = _mm_sub_ps(fones, s);
                        condition1 = _mm_cmp_ps(
                                        _mm_cvtepi32_ps(
					         _mm_and_si128(_mm_add_epi32(q, ones), twos)),
                                                    fzeroes,
                                                       _CMP_NEQ_UQ);
                        condition2 = _mm_cmp_ps(
                                          _mm_cmp_ps(
                                             _mm_cvtepi32_ps(
					              _mm_and_si128(q, fours)), fzeroes, _CMP_NEQ_UQ),
                                                       _mm_cmp_ps(aVal, fzeroes, _CMP_LT_OS),
                                                                                            _CMP_NEQ_UQ);
			sine =
                              _mm_add_ps(sine,
			               _mm_and_ps(_mm_sub_ps(cosine, sine), condition1));
                        sine = _mm_sub_ps(
                              sine, _mm_and_ps(_mm_mul_ps(sine, ftwos), condition2));
			_mm_storeu_ps(b,sine);
			a += 4;
			b += 4;
		 }
		 idx = len*4;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(2),max(4)
#endif
              for(; idx != npoints; ++idx) {
                  b[i] = ceph_sinf(a[i]);
	      }		 
         }



	 __ATTR_ALWAYS_INLINE__
	 static inline
	 void dsp_32f_sin_32f_sse_a_looped(float * __restrict ___ATTR_ALIGN__(16) b,
	                                   const float * __restrict __ATTR_ALIGN__(16) a,
					   const int32_t npoints) {
              DSP_32F_SIN_32F_SSE_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
              __assume_aligned(b,16);
	      __assume_aligned(a,16);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              b = (float*)__builtin_assume_aligned(b,16);
	      a = (float*)__builtin_assume_aligned(a,16);
#endif		      
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                  for(; idx != len; ++idx) {
                       _mm_prefetch((const char*)&a+32,_MM_HINT_T0);
		       aVal = _mm_load_ps(a);
		       s    = _mm_sub_ps(aVal,
                                    _mm_and_ps(_mm_mul_ps(aVal, ftwos),
                                                  _mm_cmp_ps(aVal, fzeroes, _CMP_LT_OQ)));
                       q = _mm_cvtps_epi32(_mm_floor_ps(_mm_mul_ps(s, m4pi)));
                       r = _mm_add_epi32(q, _mm_and_si128(q, ones));
                       s = _mm_fnmadd_ps(_mm_cvtepi32_ps(r), pio4A, s);
                       s = _mm_fnmadd_ps(_mm_cvtepi32_ps(r), pio4B, s);
		       s = _mm_mul_ps(s,finv8);
		       s = _mm_mul_ps(s,s);
		       s = _mm_mul_ps(
                                    _mm_fmadd_ps(
                                             _mm_fmsub_ps(
                                                  _mm_fmadd_ps(
						           _mm_fmsub_ps(s, cp5, cp4), s, cp3), s, cp2),
                                                  s,
                                            cp1),
                                     s);
		        t0 = _mm_sub_ps(ffours,s);
			s  = _mm_mul_ps(s,t0);
			s  = _mm_mul_ps(s,t0);
			s  = _mm_mul_ps(s,t0);
			s  = _mm_mul_ps(s,finv2);
			sine = _mm_sqrt_ps(
			               _mm_mul_ps(_mm_sub_ps(ftwos, s), s));
                        cosine = _mm_sub_ps(fones, s);
                        condition1 = _mm_cmp_ps(
                                        _mm_cvtepi32_ps(
					         _mm_and_si128(_mm_add_epi32(q, ones), twos)),
                                                    fzeroes,
                                                       _CMP_NEQ_UQ);
                        condition2 = _mm_cmp_ps(
                                          _mm_cmp_ps(
                                             _mm_cvtepi32_ps(
					              _mm_and_si128(q, fours)), fzeroes, _CMP_NEQ_UQ),
                                                       _mm_cmp_ps(aVal, fzeroes, _CMP_LT_OS),
                                                                                            _CMP_NEQ_UQ);
			sine =
                              _mm_add_ps(sine,
			               _mm_and_ps(_mm_sub_ps(cosine, sine), condition1));
                        sine = _mm_sub_ps(
                              sine, _mm_and_ps(_mm_mul_ps(sine, ftwos), condition2));
			_mm_store_ps(b,sine);
			a += 4;
			b += 4;
		 }
		 idx = len*4;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(2),max(4)
#endif
              for(; idx != npoints; ++idx) {
                  b[i] = ceph_sinf(a[i]);
	      }		 
         }



	  __ATTR_ALWAYS_INLINE__
	  static inline
          __m128 dsp_32f_sin_32_sse(const __m128 v) {
                
	          const __m128 m4pio   =  _mm_set1_ps(1.273239545F);           
	          const __m128 pio4A   =  _mm_set1_ps(0.78515625F);            
	          const __m128 pio4B   =  _mm_set1_ps(0.241876e-3F);           
	          const __m128 finv8   =  _mm_set1_ps(0.125F);                  
	          const __m128 ffours  = _mm_set1_ps(4.0F);                    
	          const __m128 fftwos  = _mm_set1_ps(2.0F);                  
	          const __m128 fones   = _mm_set1_ps(1.0F);                    
	          const __m128 fzeroes = _mm_setzero_ps();                     
	          const __m128 finv2   = _mm_set1_ps(0.5f);                   
	          const __m128 cp1     = _mm_set1_ps(1.0F);                    
	          const __m128 cp2     = _mm_set1_ps(0.83333333e-1F);          
	          const __m128 cp3     = _mm_set1_ps(0.2777778e-2F);           
	          const __m128 cp4     = _mm_set1_ps(0.49603e-4F);              
                  const __m128 cp5     = _mm_set1_ps(0.551e-6F);               
	          const __m128i ones   = _mm_set1_epi32(1);                     
	          const __m128i twos   = _mm_set1_epi32(2);                    
	          const __m128i fours  = _mm_set1_epi32(4);                    
	          __m128 sine = fzeroes;                                        
	          __m128 s    = fzeroes;                                          
	          __m128 cosine = fzeroes;                                        
	          __m128 t0     = fzeroes;                                      
	          __m128 condition1;                                              
	          __m128 condition2;                                             
	          __m128i q;                                                     
	          __m128i r;
		  aVal = v;
		  s    = _mm_sub_ps(aVal,
                                    _mm_and_ps(_mm_mul_ps(aVal, ftwos),
                                                  _mm_cmp_ps(aVal, fzeroes, _CMP_LT_OQ)));
                  q = _mm_cvtps_epi32(_mm_floor_ps(_mm_mul_ps(s, m4pi)));
                  r = _mm_add_epi32(q, _mm_and_si128(q, ones));
                  s = _mm_fnmadd_ps(_mm_cvtepi32_ps(r), pio4A, s);
                  s = _mm_fnmadd_ps(_mm_cvtepi32_ps(r), pio4B, s);
		  s = _mm_mul_ps(s,finv8);
		  s = _mm_mul_ps(s,s);
		  s = _mm_mul_ps(
                                    _mm_fmadd_ps(
                                             _mm_fmsub_ps(
                                                  _mm_fmadd_ps(
						           _mm_fmsub_ps(s, cp5, cp4), s, cp3), s, cp2),
                                                  s,
                                            cp1),
                                     s);
		   t0 = _mm_sub_ps(ffours,s);
		   s  = _mm_mul_ps(s,t0);
		   s  = _mm_mul_ps(s,t0);
		   s  = _mm_mul_ps(s,t0);
		   s  = _mm_mul_ps(s,finv2);
		   sine = _mm_sqrt_ps(
			               _mm_mul_ps(_mm_sub_ps(ftwos, s), s));
                   cosine = _mm_sub_ps(fones, s);
                   condition1 = _mm_cmp_ps(
                                        _mm_cvtepi32_ps(
					         _mm_and_si128(_mm_add_epi32(q, ones), twos)),
                                                    fzeroes,
                                                       _CMP_NEQ_UQ);
                    condition2 = _mm_cmp_ps(
                                          _mm_cmp_ps(
                                             _mm_cvtepi32_ps(
					              _mm_and_si128(q, fours)), fzeroes, _CMP_NEQ_UQ),
                                                       _mm_cmp_ps(aVal, fzeroes, _CMP_LT_OS),
                                                                                            _CMP_NEQ_UQ);
		    sine =
                              _mm_add_ps(sine,
			               _mm_and_ps(_mm_sub_ps(cosine, sine), condition1));
                    sine = _mm_sub_ps(
                              sine, _mm_and_ps(_mm_mul_ps(sine, ftwos), condition2));
		    return (sine);
          }
	  

	  







#endif /*__GMS_32F_SIN_32F_SSE_HPP__*/
