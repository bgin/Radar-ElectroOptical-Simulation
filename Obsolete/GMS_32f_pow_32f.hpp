

#ifndef __GMS_32F_POW_32F_HPP__
#define __GMS_32F_POW_32F_HPP__ 031020211027


/*
    Based on VOLK project.
    Modified by Bernard Gingold on:
    Date: 03-10-2021 10:27AM +00200
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

#define POW_POLY_DEGREE 3

#define POLY0_AVX2(x, c0) _mm256_set1_ps(c0)
#define POLY1_AVX2(x, c0, c1) \
    _mm256_fmadd_ps(POLY0_AVX2(x, c1), x, _mm256_set1_ps(c0))
#define POLY2_AVX2(x, c0, c1, c2) \
    _mm256_fmadd_ps(POLY1_AVX2(x, c1, c2), x, _mm256_set1_ps(c0))
#define POLY3_AVX2(x, c0, c1, c2, c3) \
    _mm256_fmadd_ps(POLY2_AVX2(x, c1, c2, c3), x, _mm256_set1_ps(c0))
#define POLY4_AVX2(x, c0, c1, c2, c3, c4) \
    _mm256_fmadd_ps(POLY3_AVX2(x, c1, c2, c3, c4), x, _mm256_set1_ps(c0))
#define POLY5_AVX2(x, c0, c1, c2, c3, c4, c5) \
    _mm256_fmadd_ps(POLY4_AVX2(x, c1, c2, c3, c4, c5), x, _mm256_set1_ps(c0))


#if !defined(DSP_32F_POW_32F_BLOCK)
#define DSP_32F_POW_32F_BLOCK                                                  \
        const __m256  one       = _mm256_set1_ps(1.0F);                        \
        const __m256  exp_hi    = _mm256_set1_ps(88.3762626647949F);           \
        const __m256  exp_lo    = _mm256_set1_ps(-88.3762626647949F);          \
        const __m256  ln2       = _mm256_set1_ps(0.6931471805F);               \
        const __m256  log2EF    = _mm256_set1_ps(1.44269504088896341F);        \
        const __m256  half      = _mm256_set1_ps(0.5F);                        \
        const __m256  exp_C1    = _mm256_set1_ps(0.693359375F);                \
        const __m256  exp_C2    = _mm256_set1_ps(-2.12194440e-4F);             \
        const __m26i  pi32_0x7f = _mm256_set1_epi32(0x7f);                     \
        const __m256  exp_p0    = _mm256_set1_ps(1.9875691500e-4F);            \
        const __m256  exp_p1    = _mm256_set1_ps(1.3981999507e-3F);            \
        const __m256  exp_p2    = _mm256_set1_ps(8.3334519073e-3F);            \
        const __m256  exp_p3    = _mm256_set1_ps(4.1665795894e-2F);            \ 
        const __m256  exp_p4    = _mm256_set1_ps(1.6666665459e-1F);            \
        const __m256  exp_p5    = _mm256_set1_ps(5.0000001201e-1F);            \
        __m256 aVal             = _mm256_setzero_ps();                         \
	__m256 bVal             = aVal;                                        \
	__m256 cVal             = aVal;                                        \
	__m256 logarithm        = aVal;                                        \
	__m256 mantissa         = aVal;                                        \
	__m256 frac             = aVal;                                        \
	__m256 leadingOne       = aVal;                                        \
	__m256 tmp              = aVal;                                        \
	__m256 fx               = aVal;                                        \
	__m256 mask             = aVal;                                        \
	__m256 pow2n            = aVal;                                        \
	__m256 z                = aVal;                                        \
	__m256 y                = aVal;                                        \
	__m256i bias;                                                          \
	__m256i exp;                                                           \
	__m256i emm0;                                                          \
	int32_t idx = 0;                                                       \
	const int32_t len = npoints/8;
#endif




         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_pow_32f_u_avx2_looped(float * __restrict c,
	                                    const float * __restrict b,
				            const float * __restrict a,
				            const int32_t npoints) {
               DSP_32F_POW_32F_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                for(; idx != len; ++idx) {
                      _mm_prefetch((const char*)&a+32,_MM_HINT_T0);
		      aVal = _mm256_loadu_ps(a);
		      bias = _mm256_set1_epi32(127);
		      leadingOne = one;
		      exp = _mm256_sub_epi32(
                                    _mm256_srli_epi32(_mm256_and_si256(_mm256_castps_si256(aVal),
                                                                        _mm256_set1_epi32(0x7f800000)),
                                       23),
                             bias);
                      logarithm = _mm256_cvtepi32_ps(exp);
                      frac = _mm256_or_ps(
                                     leadingOne,
                                         _mm256_and_ps(aVal,
					          _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffff))));
#if POW_POLY_DEGREE == 6
        mantissa = POLY5_AVX2(frac,
                                  3.1157899f,
                                  -3.3241990f,
                                  2.5988452f,
                                  -1.2315303f,
                                  3.1821337e-1f,
                                  -3.4436006e-2f);
#elif POW_POLY_DEGREE == 5
        mantissa = POLY4_AVX2(frac,
                                  2.8882704548164776201f,
                                  -2.52074962577807006663f,
                                  1.48116647521213171641f,
                                  -0.465725644288844778798f,
                                  0.0596515482674574969533f);
#elif POW_POLY_DEGREE == 4
        mantissa = POLY3_AVX2(frac,
                                  2.61761038894603480148f,
                                  -1.75647175389045657003f,
                                  0.688243882994381274313f,
                                  -0.107254423828329604454f);
#elif POW_POLY_DEGREE == 3
        mantissa = POLY2_AVX2(frac,
                                  2.28330284476918490682f,
                                  -1.04913055217340124191f,
                                  0.204446009836232697516f);
#else
#error
#endif
                         logarithm = _mm256_fmadd_ps(mantissa,
			                   _mm256_sub_ps(frac, leadingOne), logarithm);
                         logarithm = _mm256_mul_ps(logarithm, ln2);
			 _mm_prefetch((const char*)&b+32,_MM_HINT_T0);
			 bVal = _mm256_loadu_ps(b);
			 bVal = _mm256_mul_ps(bVal,logarithm);
			 bVal = _mm256_max_ps(_mm256_min_ps(bVal, exp_hi), exp_lo);
                         fx = _mm256_fmadd_ps(bVal, log2EF, half);
                         emm0 = _mm256_cvttps_epi32(fx);
                         tmp = _mm256_cvtepi32_ps(emm0);
                         mask = _mm256_and_ps(_mm256_cmp_ps(tmp, fx, _CMP_GT_OQ), one);
                         fx = _mm256_sub_ps(tmp, mask);
                         tmp = _mm256_fnmadd_ps(fx, exp_C1, bVal);
                         bVal = _mm256_fnmadd_ps(fx, exp_C2, tmp);
                         z = _mm256_mul_ps(bVal, bVal);
			 emm0 =
                             _mm256_slli_epi32(_mm256_add_epi32(
			                          _mm256_cvttps_epi32(fx), pi32_0x7f), 23);
			 pow2n = _mm256_castsi256_ps(emm0);
                         y = _mm256_fmadd_ps(exp_p0, bVal, exp_p1);
                         y = _mm256_fmadd_ps(y, bVal, exp_p2);
                         y = _mm256_fmadd_ps(y, bVal, exp_p3);
                         y = _mm256_fmadd_ps(y, bVal, exp_p4);
                         y = _mm256_fmadd_ps(y, bVal, exp_p5);
                         y = _mm256_fmadd_ps(y, z, bVal);
                         y = _mm256_add_ps(y, one);
                         cVal = _mm256_mul_ps(y, pow2n);
			 _mm256_storeu_ps(c,cVal);
			 a += 8;
			 b += 8;
			 c += 8;
		}
			idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#endif
              for(; idx != npoints; ++idx) {
                  c[i] = ceph_powf(a[i],b[i]);
	      }	
	 }



	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void dsp_32f_pow_32f_a_avx2_looped(float * __restrict __ATTR_ALIGN__(32) c,
	                                    float * __restrict __ATTR_ALIGN__(32) b,
				            float * __restrict __ATTR_ALIGN__(32) a,
				            const int32_t npoints) {
               DSP_32F_POW_32F_BLOCK
#if defined __ICC || defined __INTEL_COMPILER
              __assume_aligned(c,32);
              __assume_aligned(b,32);
	      __assume_aligned(a,32);
#elif defined __GNUC__ && !defined __INTEL_COMPILER
              c = (float*)__builtin_assume_aligned(c,32);
              b = (float*)__builtin_assume_aligned(b,32);
	      a = (float*)__builtin_assume_aligned(a,32);
#endif	
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                for(; idx != len; ++idx) {
                      _mm_prefetch((const char*)&a+32,_MM_HINT_T0);
		      aVal = _mm256_load_ps(a);
		      bias = _mm256_set1_epi32(127);
		      leadingOne = one;
		      exp = _mm256_sub_epi32(
                                    _mm256_srli_epi32(_mm256_and_si256(_mm256_castps_si256(aVal),
                                                                        _mm256_set1_epi32(0x7f800000)),
                                       23),
                             bias);
                      logarithm = _mm256_cvtepi32_ps(exp);
                      frac = _mm256_or_ps(
                                     leadingOne,
                                         _mm256_and_ps(aVal,
					          _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffff))));
#if POW_POLY_DEGREE == 6
        mantissa = POLY5_AVX2(frac,
                                  3.1157899f,
                                  -3.3241990f,
                                  2.5988452f,
                                  -1.2315303f,
                                  3.1821337e-1f,
                                  -3.4436006e-2f);
#elif POW_POLY_DEGREE == 5
        mantissa = POLY4_AVX2(frac,
                                  2.8882704548164776201f,
                                  -2.52074962577807006663f,
                                  1.48116647521213171641f,
                                  -0.465725644288844778798f,
                                  0.0596515482674574969533f);
#elif POW_POLY_DEGREE == 4
        mantissa = POLY3_AVX2(frac,
                                  2.61761038894603480148f,
                                  -1.75647175389045657003f,
                                  0.688243882994381274313f,
                                  -0.107254423828329604454f);
#elif POW_POLY_DEGREE == 3
        mantissa = POLY2_AVX2(frac,
                                  2.28330284476918490682f,
                                  -1.04913055217340124191f,
                                  0.204446009836232697516f);
#else
#error
#endif
                         logarithm = _mm256_fmadd_ps(mantissa,
			                   _mm256_sub_ps(frac, leadingOne), logarithm);
                         logarithm = _mm256_mul_ps(logarithm, ln2);
			 _mm_prefetch((const char*)&b+32,_MM_HINT_T0);
			 bVal = _mm256_load_ps(b);
			 bVal = _mm256_mul_ps(bVal,logarithm);
			 bVal = _mm256_max_ps(_mm256_min_ps(bVal, exp_hi), exp_lo);
                         fx = _mm256_fmadd_ps(bVal, log2EF, half);
                         emm0 = _mm256_cvttps_epi32(fx);
                         tmp = _mm256_cvtepi32_ps(emm0);
                         mask = _mm256_and_ps(_mm256_cmp_ps(tmp, fx, _CMP_GT_OQ), one);
                         fx = _mm256_sub_ps(tmp, mask);
                         tmp = _mm256_fnmadd_ps(fx, exp_C1, bVal);
                         bVal = _mm256_fnmadd_ps(fx, exp_C2, tmp);
                         z = _mm256_mul_ps(bVal, bVal);
			 emm0 =
                             _mm256_slli_epi32(_mm256_add_epi32(
			                          _mm256_cvttps_epi32(fx), pi32_0x7f), 23);
			 pow2n = _mm256_castsi256_ps(emm0);
                         y = _mm256_fmadd_ps(exp_p0, bVal, exp_p1);
                         y = _mm256_fmadd_ps(y, bVal, exp_p2);
                         y = _mm256_fmadd_ps(y, bVal, exp_p3);
                         y = _mm256_fmadd_ps(y, bVal, exp_p4);
                         y = _mm256_fmadd_ps(y, bVal, exp_p5);
                         y = _mm256_fmadd_ps(y, z, bVal);
                         y = _mm256_add_ps(y, one);
                         cVal = _mm256_mul_ps(y, pow2n);
			 _mm256_store_ps(c,cVal);
			 a += 8;
			 b += 8;
			 c += 8;
		}
			idx = len*8;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
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
          __m256 dsp_32f_pow_32_avx2(const __m256 x,
	                             const __m256 y  ) {

                    const __m256  one       = _mm256_set1_ps(1.0F);                       
                    const __m256  exp_hi    = _mm256_set1_ps(88.3762626647949F);          
                    const __m256  exp_lo    = _mm256_set1_ps(-88.3762626647949F);         
                    const __m256  ln2       = _mm256_set1_ps(0.6931471805F);              
                    const __m256  log2EF    = _mm256_set1_ps(1.44269504088896341F);      
                    const __m256  half      = _mm256_set1_ps(0.5F);                       
                    const __m256  exp_C1    = _mm256_set1_ps(0.693359375F);               
                    const __m256  exp_C2    = _mm256_set1_ps(-2.12194440e-4F);            
                    const __m26i  pi32_0x7f = _mm256_set1_epi32(0x7f);                    
                    const __m256  exp_p0    = _mm256_set1_ps(1.9875691500e-4F);           
                    const __m256  exp_p1    = _mm256_set1_ps(1.3981999507e-3F);           
                    const __m256  exp_p2    = _mm256_set1_ps(8.3334519073e-3F);           
                    const __m256  exp_p3    = _mm256_set1_ps(4.1665795894e-2F);           
                    const __m256  exp_p4    = _mm256_set1_ps(1.6666665459e-1F);          
                    const __m256  exp_p5    = _mm256_set1_ps(5.0000001201e-1F);           
                    __m256 aVal             = _mm256_setzero_ps();                        
	            __m256 bVal             = aVal;                                       
	            __m256 cVal             = aVal;                                       
	            __m256 logarithm        = aVal;                                      
	            __m256 mantissa         = aVal;                                       
	            __m256 frac             = aVal;                                      
	            __m256 leadingOne       = aVal;                                       
	            __m256 tmp              = aVal;                                       
	            __m256 fx               = aVal;                                       
	            __m256 mask             = aVal;                                      
	            __m256 pow2n            = aVal;                                      
	            __m256 z                = aVal;                                      
	            __m256 y                = aVal;                                       
	            __m256i bias;                                                         
	            __m256i exp;                                                         
	            __m256i emm0;
		    
		    aVal = x;
		    bias = _mm256_set1_epi32(127);
		    leadingOne = one;
		    exp = _mm256_sub_epi32(
                                    _mm256_srli_epi32(_mm256_and_si256(_mm256_castps_si256(aVal),
                                                                        _mm256_set1_epi32(0x7f800000)),
                                       23),
                             bias);
                    logarithm = _mm256_cvtepi32_ps(exp);
                    frac = _mm256_or_ps(
                                     leadingOne,
                                         _mm256_and_ps(aVal,
					          _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffff))));
#if POW_POLY_DEGREE == 6
        mantissa = POLY5_AVX2(frac,
                                  3.1157899f,
                                  -3.3241990f,
                                  2.5988452f,
                                  -1.2315303f,
                                  3.1821337e-1f,
                                  -3.4436006e-2f);
#elif POW_POLY_DEGREE == 5
        mantissa = POLY4_AVX2(frac,
                                  2.8882704548164776201f,
                                  -2.52074962577807006663f,
                                  1.48116647521213171641f,
                                  -0.465725644288844778798f,
                                  0.0596515482674574969533f);
#elif POW_POLY_DEGREE == 4
        mantissa = POLY3_AVX2(frac,
                                  2.61761038894603480148f,
                                  -1.75647175389045657003f,
                                  0.688243882994381274313f,
                                  -0.107254423828329604454f);
#elif POW_POLY_DEGREE == 3
        mantissa = POLY2_AVX2(frac,
                                  2.28330284476918490682f,
                                  -1.04913055217340124191f,
                                  0.204446009836232697516f);
#else
#error
#endif
                       logarithm = _mm256_fmadd_ps(mantissa,
			                   _mm256_sub_ps(frac, leadingOne), logarithm);
                       logarithm = _mm256_mul_ps(logarithm, ln2);
		       bVal = y;
		       bVal = _mm256_mul_ps(bVal,logarithm);
		       bVal = _mm256_max_ps(_mm256_min_ps(bVal, exp_hi), exp_lo);
                       fx = _mm256_fmadd_ps(bVal, log2EF, half);
                       emm0 = _mm256_cvttps_epi32(fx);
                       tmp = _mm256_cvtepi32_ps(emm0);
                       mask = _mm256_and_ps(_mm256_cmp_ps(tmp, fx, _CMP_GT_OQ), one);
                       fx = _mm256_sub_ps(tmp, mask);
                       tmp = _mm256_fnmadd_ps(fx, exp_C1, bVal);
                       bVal = _mm256_fnmadd_ps(fx, exp_C2, tmp);
                       z = _mm256_mul_ps(bVal, bVal);
		       emm0 =
                             _mm256_slli_epi32(_mm256_add_epi32(
			                          _mm256_cvttps_epi32(fx), pi32_0x7f), 23);
		       pow2n = _mm256_castsi256_ps(emm0);
                       y = _mm256_fmadd_ps(exp_p0, bVal, exp_p1);
                       y = _mm256_fmadd_ps(y, bVal, exp_p2);
                       y = _mm256_fmadd_ps(y, bVal, exp_p3);
                       y = _mm256_fmadd_ps(y, bVal, exp_p4);
                       y = _mm256_fmadd_ps(y, bVal, exp_p5);
                       y = _mm256_fmadd_ps(y, z, bVal);
                       y = _mm256_add_ps(y, one);
                       cVal = _mm256_mul_ps(y, pow2n);
		       return (cVal);
	  }


#endif /*__GMS_32F_POW_32F_HPP__*/
