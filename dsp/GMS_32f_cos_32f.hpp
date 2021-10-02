
#ifndef __GMS_32F_COS_32F_HPP__
#define __GMS_32F_COS_32F_HPP__ 021020210841


/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 02-10-2021 08:41PM +00200
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



#if !defined(DSP_32F_COS_32F_BLOCK)
#define DSP_32F_COS_32F_BLOCK                                                                                       \
        const __m256 m4pi    = _mm256_set1_ps(1.273239544735162542821171882678754627704620361328125);               \
	const __m256 pio4A   = _mm256_set1_ps(0.7853981554508209228515625);                                         \                            
	const __m256 pio4B   = _mm256_set1_ps(0.794662735614792836713604629039764404296875e-8);                     \
	const __m256 pio4C   = _mm256_set1_ps(0.306161699786838294306516483068750264552437361480769e-16);           \
	const __m256 feight  = _mm256_set1_ps(8.0F);                                                                 \
	const __m256 ffours  = _mm256_set1_ps(4.0F);                                                                \
	const __m256 ftwos   = _mm256_set1_ps(2.0F);                                                                \
	const __m256 fones   = _mm256_set1_ps(1.0F);                                                                \
	const __m256 cp1     = _mm256_set1_ps(1.0F);                                                                \
	const __m256 finv2   = _mm256_set1_ps(0.5F);								    \
	const __m256 cp2     = _mm256_set1_ps(0.08333333333333333F);                                                \
	const __m256 cp2     = _mm256_set1_ps(0.002777777777777778F);                                               \
	const __m256 cp3     = _mm256_set1_ps(0.002777777777777778F);                                               \
	const __m256 cp4     = _mm256_set1_ps(4.96031746031746e-05F);                                               \                                        
	const __m256 cp5     = _mm256_set1_ps(5.511463844797178e-07F);                                              \
	const __m256i zeroes = _mm256_set1_epi32(0);                                                                \
	const __m256i allones = _mm256_set1_epi32(0xFFFFFFFF);                                                      \
	const __m256i twos    = _mm256_set1_epi32(2);                                                               \
	const __m256i fours   = _mm256_set1_epi32(4);                                                               \
	__m256 fzeroes      = _mm256_setzero_ps();                                                                  \
	__m256 aVal         = fzeroes;                                                                              \
	__m256 s;                                                                                                   \
	__m256 r;                                                                                                   \
	__m256 t0 = fzeroes;                                                                                        \
	__m256 sine = fzeroes;                                                                                      \
	__m256 cosine = fzeroes;					                                            \
	__m256i q;					                                                            \
	int32_t idx = 0;                                                                                            \
	const int32_t len = npoints/8;                                                                        
#endif



         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_cos_32f_avx2_u_looped(float * __restrict b,
	                                   const float * __restrict a,
					   const int32_t npoints) {

                  DSP_32F_COS_32F_BLOCK
		  union bit256 {
                       uint8_t i8[32];
                       uint16_t i16[16];
                       uint32_t i[8];
                       float f[8];
                       double d[4];
#ifdef __AVX__
                       __m256 float_vec;
                       __m256i int_vec;
                       __m256d double_vec;
#endif
                  };
		  union bit256 condition1;
		  union bit256 condition2;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                  for(; idx != len; ++idx) {
                      _mm_prefetch((const char*)&a+32,_MM_HINT_T0);
		      aVal = _mm256_loadu_ps(a);
		      s = _mm256_sub_ps(aVal,
                                 _mm256_and_ps(_mm256_mul_ps(aVal, ftwos),
                                                 _mm256_cmp_ps(aVal, fzeroes, _CMP_LT_OS)));
                      q = _mm256_cvtps_epi32(_mm256_floor_ps(
		                                 _mm256_mul_ps(s, m4pi)));
                      r = _mm256_cvtepi32_ps(_mm256_add_epi32(q,
		                                 _mm256_and_si256(q, ones)));
                      s = _mm256_fnmadd_ps(r, pio4A, s);
                      s = _mm256_fnmadd_ps(r, pio4B, s);
                      s = _mm256_fnmadd_ps(r, pio4C, s);
		      s = _mm256_div_ps(s,feight);
		      s = _mm256_mul_ps(s,s);
		      s = _mm256_mul_ps(s,finv2);
                                _mm256_fmadd_ps(
                                       _mm256_fmsub_ps(
                                              _mm256_fmadd_ps(
					               _mm256_fmsub_ps(s, cp5, cp4), s, cp3), s, cp2),
                                              s,
                                       cp1),
                               s);
		      t0 = _mm256_sub_ps(ffours,s);
		      s = _mm256_mul_ps(s,t0);
		      s = _mm256_mul_ps(s,t0);
		      s = _mm256_mul_ps(s,t0);
		      s = _mm256_mul_ps(
		      sine = _mm256_sqrt_ps(_mm256_mul_ps(
		                        _mm256_sub_ps(ftwos, s), s));
                      cosine = _mm256_sub_ps(fones, s);
		      condition1.int_vec =  _mm256_cmpeq_epi32(_mm256_and_si256(
		                            _mm256_add_epi32(q, ones), twos), zeroes);
		      condition1.int_vec = _mm256_xor_si256(allones, condition1.int_vec);
		      condition2.int_vec = _mm256_cmpeq_epi32(
                                         _mm256_and_si256(
					           _mm256_add_epi32(q, twos), fours), zeroes);
		      condition2.int_vec =  _mm256_xor_si256(allones, condition2.int_vec);
                      cosine = _mm256_add_ps(
                                           cosine, _mm256_and_ps(
					             _mm256_sub_ps(sine, cosine),
						              condition1.float_vec));
                      cosine = _mm256_sub_ps(cosine,
                                           _mm256_and_ps(
					             _mm256_mul_ps(cosine,ftwos),
                                                                condition2.float_vec));
                      _mm256_storeu_ps(b,cosine);
		      a += 8;
		      b += 8;
		  }
		  idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
              for(; idx != npoints; ++idx) {
                  b[i] = ceph_cosf(a[i]);
	      }
	 }


         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_cos_32f_avx2_a_looped(float * __restrict __ATTR_ALIGN__(32) b,
	                                   const float * __restrict __ATTR_ALIGN__(32) a,
					   const int32_t npoints) {

                  DSP_32F_COS_32F_BLOCK
		  union bit256 {
                       uint8_t i8[32];
                       uint16_t i16[16];
                       uint32_t i[8];
                       float f[8];
                       double d[4];
#ifdef __AVX__
                       __m256 float_vec;
                       __m256i int_vec;
                       __m256d double_vec;
#endif
                  };
		  union bit256 condition1;
		  union bit256 condition2;
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
                      q = _mm256_cvtps_epi32(_mm256_floor_ps(
		                                 _mm256_mul_ps(s, m4pi)));
                      r = _mm256_cvtepi32_ps(_mm256_add_epi32(q,
		                                 _mm256_and_si256(q, ones)));
                      s = _mm256_fnmadd_ps(r, pio4A, s);
                      s = _mm256_fnmadd_ps(r, pio4B, s);
                      s = _mm256_fnmadd_ps(r, pio4C, s);
		      s = _mm256_div_ps(s,feight);
		      s = _mm256_mul_ps(s,s);
		      s = _mm256_mul_ps(
                                _mm256_fmadd_ps(
                                       _mm256_fmsub_ps(
                                              _mm256_fmadd_ps(
					               _mm256_fmsub_ps(s, cp5, cp4), s, cp3), s, cp2),
                                              s,
                                       cp1),
                               s);
		      t0 = _mm256_sub_ps(ffours,s);
		      s = _mm256_mul_ps(s,t0);
		      s = _mm256_mul_ps(s,t0);
		      s = _mm256_mul_ps(s,t0);
		      s = _mm256_mul_ps(s,finv2);
		      sine = _mm256_sqrt_ps(_mm256_mul_ps(
		                        _mm256_sub_ps(ftwos, s), s));
                      cosine = _mm256_sub_ps(fones, s);
		      condition1.int_vec =  _mm256_cmpeq_epi32(_mm256_and_si256(
		                            _mm256_add_epi32(q, ones), twos), zeroes);
		      condition1.int_vec = _mm256_xor_si256(allones, condition1.int_vec);
		      condition2.int_vec = _mm256_cmpeq_epi32(
                                         _mm256_and_si256(
					           _mm256_add_epi32(q, twos), fours), zeroes);
		      condition2.int_vec =  _mm256_xor_si256(allones, condition2.int_vec);
                      cosine = _mm256_add_ps(
                                           cosine, _mm256_and_ps(
					             _mm256_sub_ps(sine, cosine),
						              condition1.float_vec));
                      cosine = _mm256_sub_ps(cosine,
                                           _mm256_and_ps(
					             _mm256_mul_ps(cosine,ftwos),
                                                                condition2.float_vec));
                      _mm256_store_ps(b,cosine);
		      a += 8;
		      b += 8;
		  }
		  idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
              for(; idx != npoints; ++idx) {
                  b[i] = ceph_cosf(a[i]);
	      }
	 }



          __ATTR_ALWAYS_INLINE__
	  __ATTR_HOT__
	  __ATTR_ALIGN__(32)
	  __ATTR_VECTORCALL__
	  __attribute__((regcall)) // GCC will skip over this attribute!!
	  static inline
          __m256 dsp_32f_cos_32_avx(const __m256 v) {

	           const __m256 m4pi    = _mm256_set1_ps(1.273239544735162542821171882678754627704620361328125);             
	           const __m256 pio4A   = _mm256_set1_ps(0.7853981554508209228515625);                                                                    
	           const __m256 pio4B   = _mm256_set1_ps(0.794662735614792836713604629039764404296875e-8);                    
	           const __m256 pio4C   = _mm256_set1_ps(0.306161699786838294306516483068750264552437361480769e-16);          
	           const __m256 feight  = _mm256_set1_ps(8.0F);                                                                
	           const __m256 ffours  = _mm256_set1_ps(4.0F);                                                               
	           const __m256 ftwos   = _mm256_set1_ps(2.0F);                                                               
	           const __m256 fones   = _mm256_set1_ps(1.0F);                                                               
	           const __m256 cp1     = _mm256_set1_ps(1.0F);                                                              
	           const __m256 finv2   = _mm256_set1_ps(0.5F);								  
	           const __m256 cp2     = _mm256_set1_ps(0.08333333333333333F);                                              
	           const __m256 cp2     = _mm256_set1_ps(0.002777777777777778F);                                              
	           const __m256 cp3     = _mm256_set1_ps(0.002777777777777778F);                                              
	           const __m256 cp4     = _mm256_set1_ps(4.96031746031746e-05F);                                                                                     
	           const __m256 cp5     = _mm256_set1_ps(5.511463844797178e-07F);                                              
	           const __m256i zeroes = _mm256_set1_epi32(0);                                                               
	           const __m256i allones = _mm256_set1_epi32(0xFFFFFFFF);                                                     
 	           const __m256i twos    = _mm256_set1_epi32(2);                                                              
	           const __m256i fours   = _mm256_set1_epi32(4);                                                             
	           __m256 fzeroes      = _mm256_setzero_ps();                                                                 
	           __m256 aVal         = fzeroes;                                                                             
	           __m256 s;                                                                                                  
	           __m256 r;                                                                                                   
	           __m256 t0 = fzeroes;                                                                                       
	           __m256 sine = fzeroes;                                                                                     
	           __m256 cosine = fzeroes;					                                          
	           __m256i q;
		   union bit256 {
                       uint8_t i8[32];
                       uint16_t i16[16];
                       uint32_t i[8];
                       float f[8];
                       double d[4];
#ifdef __AVX__
                       __m256 float_vec;
                       __m256i int_vec;
                       __m256d double_vec;
#endif
                  };
		  union bit256 condition1;
		  union bit256 condition2;
		  aVal = v;
		  s = _mm256_sub_ps(aVal,
                                 _mm256_and_ps(_mm256_mul_ps(aVal, ftwos),
                                                 _mm256_cmp_ps(aVal, fzeroes, _CMP_LT_OS)));
                  q = _mm256_cvtps_epi32(_mm256_floor_ps(
		                                 _mm256_mul_ps(s, m4pi)));
                  r = _mm256_cvtepi32_ps(_mm256_add_epi32(q,
		                                 _mm256_and_si256(q, ones)));
                  s = _mm256_fnmadd_ps(r, pio4A, s);
                  s = _mm256_fnmadd_ps(r, pio4B, s);
                  s = _mm256_fnmadd_ps(r, pio4C, s);
		  s = _mm256_div_ps(s,feight);
		  s = _mm256_mul_ps(s,s);
		  s = _mm256_mul_ps(
                                _mm256_fmadd_ps(
                                       _mm256_fmsub_ps(
                                              _mm256_fmadd_ps(
					               _mm256_fmsub_ps(s, cp5, cp4), s, cp3), s, cp2),
                                              s,
                                       cp1),
                               s);
		   t0 = _mm256_sub_ps(ffours,s);
		   s = _mm256_mul_ps(s,t0);
		   s = _mm256_mul_ps(s,t0);
		   s = _mm256_mul_ps(s,t0);
		   s = _mm256_mul_ps(s,finv2);
		   sine = _mm256_sqrt_ps(_mm256_mul_ps(
		                        _mm256_sub_ps(ftwos, s), s));
                   cosine = _mm256_sub_ps(fones, s);
		   condition1.int_vec =  _mm256_cmpeq_epi32(_mm256_and_si256(
		                            _mm256_add_epi32(q, ones), twos), zeroes);
		   condition1.int_vec = _mm256_xor_si256(allones, condition1.int_vec);
		   condition2.int_vec = _mm256_cmpeq_epi32(
                                         _mm256_and_si256(
					           _mm256_add_epi32(q, twos), fours), zeroes);
		   condition2.int_vec =  _mm256_xor_si256(allones, condition2.int_vec);
                   cosine = _mm256_add_ps(
                                           cosine, _mm256_and_ps(
					             _mm256_sub_ps(sine, cosine),
						              condition1.float_vec));
                   cosine = _mm256_sub_ps(cosine,
                                           _mm256_and_ps(
					             _mm256_mul_ps(cosine,ftwos),
                                                                condition2.float_vec)); 
                   return (cosine);
          }








#endif /*__GMS_32F_COS_32F_HPP__*/
