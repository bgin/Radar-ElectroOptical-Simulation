

#ifndef __GMS_32F_POW_32F_AVX512_HPP__
#define __GMS_32F_POW_32F_AVX512_HPP__ 031020211244



/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 03-10-2021 12:44PM +00200
    contact: beniekg@gmail.com
    Few modification were added to original
    implementation (ICC pragmas, alignment directives and code layout rescheduled,
    unrolling completely 2-iteration for-loops)
    Converted to AVX512.
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


#define POW_POLY_DEGREE 3

#define POLY0_AVX512(x, c0) _mm512_set1_ps(c0)
#define POLY1_AVX512(x, c0, c1) \
    _mm256_fmadd_ps(POLY0_AVX512(x, c1), x, _mm512_set1_ps(c0))
#define POLY2_AVX512(x, c0, c1, c2) \
    _mm256_fmadd_ps(POLY1_AVX512(x, c1, c2), x, _mm512_set1_ps(c0))
#define POLY3_AVX512(x, c0, c1, c2, c3) \
    _mm256_fmadd_ps(POLY2_AVX512(x, c1, c2, c3), x, _mm512_set1_ps(c0))
#define POLY4_AVX512(x, c0, c1, c2, c3, c4) \
    _mm256_fmadd_ps(POLY3_AVX512(x, c1, c2, c3, c4), x, _mm512_set1_ps(c0))
#define POLY5_AVX512(x, c0, c1, c2, c3, c4, c5) \
    _mm256_fmadd_ps(POLY4_AVX512(x, c1, c2, c3, c4, c5), x, _mm512_set1_ps(c0))


#if !defined(DSP_32F_POW_32F_AVX512_BLOCK)
#define DSP_32F_POW_32F_AVX512_BLOCK                                                  \
        const __m512  one       = _mm512_set1_ps(1.0F);                        \
        const __m512  exp_hi    = _mm512_set1_ps(88.3762626647949F);           \
        const __m512  exp_lo    = _mm512_set1_ps(-88.3762626647949F);          \
        const __m512  ln2       = _mm512_set1_ps(0.6931471805F);               \
        const __m512  log2EF    = _mm512_set1_ps(1.44269504088896341F);        \
        const __m512  half      = _mm512_set1_ps(0.5F);                        \
        const __m512  exp_C1    = _mm512_set1_ps(0.693359375F);                \
        const __m512  exp_C2    = _mm512_set1_ps(-2.12194440e-4F);             \
        const __m512i  pi32_0x7f = _mm512_set1_epi32(0x7f);                     \
        const __m512  exp_p0    = _mm512_set1_ps(1.9875691500e-4F);            \
        const __m512  exp_p1    = _mm512_set1_ps(1.3981999507e-3F);            \
        const __m512  exp_p2    = _mm512_set1_ps(8.3334519073e-3F);            \
        const __m512  exp_p3    = _mm512_set1_ps(4.1665795894e-2F);            \ 
        const __m512  exp_p4    = _mm512_set1_ps(1.6666665459e-1F);            \
        const __m512  exp_p5    = _mm512_set1_ps(5.0000001201e-1F);            \
        __m512 aVal             = _mm512_setzero_ps();                         \
	__m512 bVal             = aVal;                                        \
	__m512 cVal             = aVal;                                        \
	__m512 logarithm        = aVal;                                        \
	__m512 mantissa         = aVal;                                        \
	__m512 frac             = aVal;                                        \
	__m512 leadingOne       = aVal;                                        \
	__m512 tmp              = aVal;                                        \
	__m512 fx               = aVal;                                        \
	__m512 pow2n            = aVal;                                        \
	__m512 z                = aVal;                                        \
	__m512 y                = aVal;                                        \
	__m512i bias;                                                          \
	__m512i exp;                                                           \
	__m512i emm0;                                                          \
	int32_t idx = 0;                                                       \
	const int32_t len = npoints/16;                                        \
	__mmask16 mask = 0;
#endif



         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_pow_32f_u_avx512_looped(float * __restrict c,
	                                      const float * __restrict b,
				              const float * __restrict a,
				              const int32_t npoints) {
               DSP_32F_POW_32F_AVX512_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                for(; idx != len; ++idx) {
                      __mm_prefetch((const char*)&a+32,_MM_HINT_T0);
		      aVal = _mm512_loadu_ps(a);
		      bias = _mm512_set1_epi32(127);
		      leadingOne = one;
		      exp = _mm512_sub_epi32(
                                    _mm512_srli_epi32(_mm512_and_si512(_mm512_castps_si512(aVal),
                                                                        _mm512_set1_epi32(0x7f800000)),
                                       23),
                             bias);
                      logarithm = _mm512_cvtepi32_ps(exp);
                      frac = _mm512_or_ps(
                                     leadingOne,
                                         _mm512_and_ps(aVal,
					          _mm512_castsi512_ps(_mm5512_set1_epi32(0x7fffff))));
#if POW_POLY_DEGREE == 6
        mantissa = POLY5_AVX512(frac,
                                  3.1157899f,
                                  -3.3241990f,
                                  2.5988452f,
                                  -1.2315303f,
                                  3.1821337e-1f,
                                  -3.4436006e-2f);
#elif POW_POLY_DEGREE == 5
        mantissa = POLY4_AVX512(frac,
                                  2.8882704548164776201f,
                                  -2.52074962577807006663f,
                                  1.48116647521213171641f,
                                  -0.465725644288844778798f,
                                  0.0596515482674574969533f);
#elif POW_POLY_DEGREE == 4
        mantissa = POLY3_AVX512(frac,
                                  2.61761038894603480148f,
                                  -1.75647175389045657003f,
                                  0.688243882994381274313f,
                                  -0.107254423828329604454f);
#elif POW_POLY_DEGREE == 3
        mantissa = POLY2_AVX512(frac,
                                  2.28330284476918490682f,
                                  -1.04913055217340124191f,
                                  0.204446009836232697516f);
#else
#error
#endif
                       logarithm = _mm512_fmadd_ps(mantissa,
			                   _mm512_sub_ps(frac, leadingOne), logarithm);
                       logarithm = _mm512_mul_ps(logarithm, ln2);
		       _mm_prefetch((const char*)&b+32,_MM_HINT_T0);
		       bVal = _mm512_loadu_ps(b);
		       bVal = _mm512_mul_ps(bVal,logarithm);
		       bVal = _mm512_max_ps(_mm512_min_ps(bVal, exp_hi), exp_lo);
                       fx = _mm512_fmadd_ps(bVal, log2EF, half);
                       emm0 = _mm512_cvttps_epi32(fx);
                       tmp = _mm512_cvtepi32_ps(emm0);
		       mask = _mm512_cmp_ps_mask(tmp,fx,_CMP_GT_OQ);
		       const __m512i m0  = _mm512_set1_epi16(mask);
		       fx = _mm512_sub_ps(tmp,_mm512_castsi512_ps(m0));
		       tmp = _mm512_fnmadd_ps(fx, exp_C1, bVal);
                       bVal = _mm512_fnmadd_ps(fx, exp_C2, tmp);
                       z = _mm512_mul_ps(bVal, bVal);
		       emm0 =
                             _mm512_slli_epi32(_mm512_add_epi32(
			                          _mm512_cvttps_epi32(fx), pi32_0x7f), 23);
			 pow2n = _mm512_castsi512_ps(emm0);
                         y = _mm512_fmadd_ps(exp_p0, bVal, exp_p1);
                         y = _mm512_fmadd_ps(y, bVal, exp_p2);
                         y = _mm512_fmadd_ps(y, bVal, exp_p3);
                         y = _mm512_fmadd_ps(y, bVal, exp_p4);
                         y = _mm512_fmadd_ps(y, bVal, exp_p5);
                         y = _mm512_fmadd_ps(y, z, bVal);
                         y = _mm512_add_ps(y, one);
                         cVal = _mm512_mul_ps(y, pow2n);
			 _mm512_storeu_ps(c,cVal);
			 a += 16;
			 b += 16;
			 c += 16;
		}
	        idx = len*16;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(16)
#endif
              for(; idx != npoints; ++idx) {
                  c[i] = ceph_powf(a[i],b[i]);
	      }			
         }




	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_pow_32f_a_avx512_looped(float * __restrict __ATTR_ALIGN__(64) c,
	                                      float * __restrict __ATTR_ALIGN__(64) b,
				              float * __restrict __ATTR_ALIGN__(64) a,
				              const int32_t npoints) {
               DSP_32F_POW_32F_AVX512_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
              __assume_aligned(c,64);
              __assume_aligned(b,64);
	      __assume_aligned(a,64);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              c = (float*)__builtin_assume_aligned(c,64);
              b = (float*)__builtin_assume_aligned(b,64);
	      a = (float*)__builtin_assume_aligned(a,64);
#endif	
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                for(; idx != len; ++idx) {
                      __mm_prefetch((const char*)&a+32,_MM_HINT_T0);
		      aVal = _mm512_load_ps(a);
		      bias = _mm512_set1_epi32(127);
		      leadingOne = one;
		      exp = _mm512_sub_epi32(
                                    _mm512_srli_epi32(_mm512_and_si512(_mm512_castps_si512(aVal),
                                                                        _mm512_set1_epi32(0x7f800000)),
                                       23),
                             bias);
                      logarithm = _mm512_cvtepi32_ps(exp);
                      frac = _mm512_or_ps(
                                     leadingOne,
                                         _mm512_and_ps(aVal,
					          _mm512_castsi512_ps(_mm5512_set1_epi32(0x7fffff))));
#if POW_POLY_DEGREE == 6
        mantissa = POLY5_AVX512(frac,
                                  3.1157899f,
                                  -3.3241990f,
                                  2.5988452f,
                                  -1.2315303f,
                                  3.1821337e-1f,
                                  -3.4436006e-2f);
#elif POW_POLY_DEGREE == 5
        mantissa = POLY4_AVX512(frac,
                                  2.8882704548164776201f,
                                  -2.52074962577807006663f,
                                  1.48116647521213171641f,
                                  -0.465725644288844778798f,
                                  0.0596515482674574969533f);
#elif POW_POLY_DEGREE == 4
        mantissa = POLY3_AVX512(frac,
                                  2.61761038894603480148f,
                                  -1.75647175389045657003f,
                                  0.688243882994381274313f,
                                  -0.107254423828329604454f);
#elif POW_POLY_DEGREE == 3
        mantissa = POLY2_AVX512(frac,
                                  2.28330284476918490682f,
                                  -1.04913055217340124191f,
                                  0.204446009836232697516f);
#else
#error
#endif
                       logarithm = _mm512_fmadd_ps(mantissa,
			                   _mm512_sub_ps(frac, leadingOne), logarithm);
                       logarithm = _mm512_mul_ps(logarithm, ln2);
		       _mm_prefetch((const char*)&b+32,_MM_HINT_T0);
		       bVal = _mm512_load_ps(b);
		       bVal = _mm512_mul_ps(bVal,logarithm);
		       bVal = _mm512_max_ps(_mm512_min_ps(bVal, exp_hi), exp_lo);
                       fx = _mm512_fmadd_ps(bVal, log2EF, half);
                       emm0 = _mm512_cvttps_epi32(fx);
                       tmp = _mm512_cvtepi32_ps(emm0);
		       mask = _mm512_cmp_ps_mask(tmp,fx,_CMP_GT_OQ);
		       const __m512i m0  = _mm512_set1_epi16(mask);
		       fx = _mm512_sub_ps(tmp,_mm512_castsi512_ps(m0));
		       tmp = _mm512_fnmadd_ps(fx, exp_C1, bVal);
                       bVal = _mm512_fnmadd_ps(fx, exp_C2, tmp);
                       z = _mm512_mul_ps(bVal, bVal);
		       emm0 =
                             _mm512_slli_epi32(_mm512_add_epi32(
			                          _mm512_cvttps_epi32(fx), pi32_0x7f), 23);
			 pow2n = _mm512_castsi512_ps(emm0);
                         y = _mm512_fmadd_ps(exp_p0, bVal, exp_p1);
                         y = _mm512_fmadd_ps(y, bVal, exp_p2);
                         y = _mm512_fmadd_ps(y, bVal, exp_p3);
                         y = _mm512_fmadd_ps(y, bVal, exp_p4);
                         y = _mm512_fmadd_ps(y, bVal, exp_p5);
                         y = _mm512_fmadd_ps(y, z, bVal);
                         y = _mm512_add_ps(y, one);
                         cVal = _mm512_mul_ps(y, pow2n);
			 _mm512_store_ps(c,cVal);
			 a += 16;
			 b += 16;
			 c += 16;
		}
	        idx = len*16;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(16)
#endif
              for(; idx != npoints; ++idx) {
                  c[i] = ceph_powf(a[i],b[i]);
	      }			
         }
	 


	  __ATTR_ALWAYS_INLINE__
	  __ATTR_HOT__
	  __ATTR_ALIGN__(32)
	  __ATTR_VECTORCALL__
	  __attribute__((regcall)) // GCC will skip over this attribute!!
	  static inline
          __m512 dsp_32f_pow_32_avx512(const __m512 x,
	                               const __m512 y  ) {
                    const __m512  one       = _mm512_set1_ps(1.0F);                       
                    const __m512  exp_hi    = _mm512_set1_ps(88.3762626647949F);          
                    const __m512  exp_lo    = _mm512_set1_ps(-88.3762626647949F);         
                    const __m512  ln2       = _mm512_set1_ps(0.6931471805F);              
                    const __m512  log2EF    = _mm512_set1_ps(1.44269504088896341F);       
                    const __m512  half      = _mm512_set1_ps(0.5F);                       
                    const __m512  exp_C1    = _mm512_set1_ps(0.693359375F);              
                    const __m512  exp_C2    = _mm512_set1_ps(-2.12194440e-4F);            
                    const __m512i  pi32_0x7f = _mm512_set1_epi32(0x7f);                   
                    const __m512  exp_p0    = _mm512_set1_ps(1.9875691500e-4F);           
                    const __m512  exp_p1    = _mm512_set1_ps(1.3981999507e-3F);          
                    const __m512  exp_p2    = _mm512_set1_ps(8.3334519073e-3F);          
                    const __m512  exp_p3    = _mm512_set1_ps(4.1665795894e-2F);           
                    const __m512  exp_p4    = _mm512_set1_ps(1.6666665459e-1F);          
                    const __m512  exp_p5    = _mm512_set1_ps(5.0000001201e-1F);           
                    __m512 aVal             = _mm512_setzero_ps();                        
	            __m512 bVal             = aVal;                                      
	            __m512 cVal             = aVal;                                      
	            __m512 logarithm        = aVal;                                      
	            __m512 mantissa         = aVal;                                     
	            __m512 frac             = aVal;                                     
	            __m512 leadingOne       = aVal;                                      
	            __m512 tmp              = aVal;                                      
	            __m512 fx               = aVal;                                      
	            __m512 pow2n            = aVal;                                       
	            __m512 z                = aVal;                                       
	            __m512 y                = aVal;                                      
	            __m512i bias;                                                         
	            __m512i exp;                                                          
	            __m512i emm0;                                                        
                    __mmask16 mask = 0;              
	            aVal = x;
		    bias = _mm512_set1_epi32(127);
		    leadingOne = one;
		    exp = _mm512_sub_epi32(
                                    _mm512_srli_epi32(_mm512_and_si512(_mm512_castps_si512(aVal),
                                                                        _mm512_set1_epi32(0x7f800000)),
                                       23),
                             bias);
                    logarithm = _mm512_cvtepi32_ps(exp);
                    frac = _mm512_or_ps(
                                     leadingOne,
                                         _mm512_and_ps(aVal,
					          _mm512_castsi512_ps(_mm5512_set1_epi32(0x7fffff))));
#if POW_POLY_DEGREE == 6
        mantissa = POLY5_AVX512(frac,
                                  3.1157899f,
                                  -3.3241990f,
                                  2.5988452f,
                                  -1.2315303f,
                                  3.1821337e-1f,
                                  -3.4436006e-2f);
#elif POW_POLY_DEGREE == 5
        mantissa = POLY4_AVX512(frac,
                                  2.8882704548164776201f,
                                  -2.52074962577807006663f,
                                  1.48116647521213171641f,
                                  -0.465725644288844778798f,
                                  0.0596515482674574969533f);
#elif POW_POLY_DEGREE == 4
        mantissa = POLY3_AVX512(frac,
                                  2.61761038894603480148f,
                                  -1.75647175389045657003f,
                                  0.688243882994381274313f,
                                  -0.107254423828329604454f);
#elif POW_POLY_DEGREE == 3
        mantissa = POLY2_AVX512(frac,
                                  2.28330284476918490682f,
                                  -1.04913055217340124191f,
                                  0.204446009836232697516f);
#else
#error
#endif
                       logarithm = _mm512_fmadd_ps(mantissa,
			                   _mm512_sub_ps(frac, leadingOne), logarithm);
                       logarithm = _mm512_mul_ps(logarithm, ln2);
		       bVal = y;
		       bVal = _mm512_mul_ps(bVal,logarithm);
		       bVal = _mm512_max_ps(_mm512_min_ps(bVal, exp_hi), exp_lo);
                       fx = _mm512_fmadd_ps(bVal, log2EF, half);
                       emm0 = _mm512_cvttps_epi32(fx);
                       tmp = _mm512_cvtepi32_ps(emm0);
		       mask = _mm512_cmp_ps_mask(tmp,fx,_CMP_GT_OQ);
		       const __m512i m0  = _mm512_set1_epi16(mask);
		       fx = _mm512_sub_ps(tmp,_mm512_castsi512_ps(m0));
		       tmp = _mm512_fnmadd_ps(fx, exp_C1, bVal);
                       bVal = _mm512_fnmadd_ps(fx, exp_C2, tmp);
                       z = _mm512_mul_ps(bVal, bVal);
		       emm0 =
                             _mm512_slli_epi32(_mm512_add_epi32(
			                          _mm512_cvttps_epi32(fx), pi32_0x7f), 23);
		       pow2n = _mm512_castsi512_ps(emm0);
                       y = _mm512_fmadd_ps(exp_p0, bVal, exp_p1);
                       y = _mm512_fmadd_ps(y, bVal, exp_p2);
                       y = _mm512_fmadd_ps(y, bVal, exp_p3);
                       y = _mm512_fmadd_ps(y, bVal, exp_p4);
                       y = _mm512_fmadd_ps(y, bVal, exp_p5);
                       y = _mm512_fmadd_ps(y, z, bVal);
                       y = _mm512_add_ps(y, one);
                       cVal = _mm512_mul_ps(y, pow2n);
		       return (cVal);
          }
 

















#endif /*__GMS_32F_POW_32F_AVX512_HPP__*/
