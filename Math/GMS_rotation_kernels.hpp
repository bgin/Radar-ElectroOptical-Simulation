
#ifndef __GMS_ROTATION_KERNELS_HPP__
#define __GMS_ROTATION_KERNELS_HPP__


namespace file_info {

const unsigned int gGMS_ROTATION_KERNELS_MAJOR = 1U;
const unsigned int gGMS_ROTATION_KERNELS_MINOR = 0U;
const unsigned int gGMS_ROTATION_KERNELS_MICRO = 0U;
const unsigned int gGMS_ROTATION_KERNELS_FULLVER =
       1000U*gGMS_ROTATION_KERNELS_MAJOR+
       100U*gGMS_ROTATION_KERNELS_MINOR +
       10U*gGMS_ROTATION_KERNELS_MICRO;
const char * const pgGMS_ROTATION_KERNELS_CREATION_DATE = "09-08-2021 02:25 PM +00200 (SUN 09 AUG 2021 GMT+2)";
const char * const pgGMS_ROTATION_KERNELS_BUILD_DATE    = __DATE__ ":" __TIME__;
const char * const pgGMS_ROTATION_KERNELS_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
const char * const pgGMS_ROTATION_KERNELS_DESCRIPTION   = "AVX512,AVX/AVX2 vectorized basic rotation operations.";
}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_dcm_avx512.hpp"
#include "GMS_dcm_avx2.hpp"

namespace gms {

         namespace math {

                       /*
                                 This version is *loosely based on the Fortran 90 "rotation.f90" source code
                                 implementation.
                                 
                                 *Many optimizations were applied (precomputation of common subexpression,
                                 constants folding, massive AVX/AVX2 and AVX512 vectorization.)
                                 The original authors copyright statement
                                 Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
                                 Modified      2017-2020, Martin Diehl/Max-Planck-Institut für Eisenforschung GmbH
                                 All rights reserved.
                         */


		    namespace {

		               const __m512 v16f32_2      = _mm512_set1_ps(2.0f);
                               const __m512 v16f32_n1     = _mm512_set1_ps(-1.0f);
			       const __m512 v16f32_spi    = _mm512_set1_ps(1.7724538509055160272982);
			       const __m512 v16f32_s6pi   = _mm512_set1_ps(1.381976597885341917061);
			       const __m512 v16f32_a      = _mm512_set1_ps(1.9257490199582527754939);
			       const __m512 v16f32_ap     = _mm512_set1_ps(2.1450293971110256000775);
			       const __m512 v16f32_sc     = _mm512_set1_ps(0.8977727869612861128953);
			       const __m512 v16f32_beta   = _mm512_set1_ps(0.9628745099791263877469);
			       const __m512 v16f32_r1     = _mm512_set1_ps(1.3306700394914687909256);
			       const __m512 v16f32_r2     = _mm512_set1_ps(1.4142135623730950488017);
			       const __m512 v16f32_pi12   = _mm512_set1_ps(0.2617993877991494365386);
			       const __m512 v16f32_prek   = _mm512_set1_ps(1.6434564029725030125017);
			       const __m512d v8f64_n1     = _mm512_set1_pd(-1.0);
			       const __m512d v8f64_spi    = _mm512_set1_pd(1.7724538509055160272982);
			       const __m512d v8f64_s6pi   = _mm512_set1_pd(1.381976597885341917061);
			       const __m512d v8f64_a      = _mm512_set1_pd(1.9257490199582527754939);
			       const __m512d v8f64_ap     = _mm512_set1_pd(2.1450293971110256000775);
			       const __m512d v8f64_sc     = _mm512_set1_pd(0.8977727869612861128953);
			       const __m512d v8f64_beta   = _mm512_set1_pd(0.9628745099791263877469);
			       const __m512d v8f64_r1     = _mm512_set1_pd(1.3306700394914687909256);
			       const __m512d v8f64_r2     = _mm512_set1_pd(1.4142135623730950488017);
			       const __m512d v8f64_pi12   = _mm512_set1_pd(0.2617993877991494365386);
			       const __m512d v8f64_prek   = _mm512_set1_pd(1.6434564029725030125017);
			       const __m256 v8f32_n1      = _mm256_set1_ps(-1.0f);
			       const __m256 v8f32_spi     = _mm256_set1_ps(1.7724538509055160272982);
			       const __m256 v8f32_s6pi    = _mm256_set1_ps(1.381976597885341917061);
			       const __m256 v8f32_a       = _mm256_set1_ps(1.9257490199582527754939);
			       const __m256 v8f32_ap      = _mm256_set1_ps(2.1450293971110256000775);
			       const __m256 v8f32_sc      = _mm256_set1_ps(0.8977727869612861128953);
			       const __m256 v8f32_beta    = _mm256_set1_ps(0.9628745099791263877469);
			       const __m256 v8f32_r1      = _mm256_set1_ps(1.3306700394914687909256);
			       const __m256 v8f32_r2      = _mm256_set1_ps(1.4142135623730950488017);
			       const __m256 v8f32_pi12    = _mm256_set1_ps(0.2617993877991494365386);
			       const __m256 v8f32_prek    = _mm256_set1_ps(1.6434564029725030125017);
			       const __m256d v4f64_n1     = _mm256_set1_ps(-1.0);
			       const __m256d v4f64_spi    = _mm256_set1_pd(1.7724538509055160272982);
			       const __m256d v4f64_s6pi   = _mm256_set1_pd(1.381976597885341917061);
			       const __m256d v4f64_a      = _mm256_set1_pd(1.9257490199582527754939);
			       const __m256d v4f64_ap     = _mm256_set1_pd(2.1450293971110256000775);
			       const __m256d v4f64_sc     = _mm256_set1_pd(0.8977727869612861128953);
			       const __m256d v4f64_beta   = _mm256_set1_pd(0.9628745099791263877469);
			       const __m256d v4f64_r1     = _mm256_set1_pd(1.3306700394914687909256);
			       const __m256d v4f64_r2     = _mm256_set1_pd(1.4142135623730950488017);
			       const __m256d v4f64_pi12   = _mm256_set1_pd(0.2617993877991494365386);
			       const __m256d v4f64_prek   = _mm256_set1_pd(1.6434564029725030125017);
	       }


	              __ATTR_VECTORCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      DCMatZMM16r4
		      q4x16_to_rmat9x16_zmm16r4(const __m512 q_x,
		                                const __m512 q_y,
						const __m512 q_z,
						const __m512 q_w) {
                            DCMatZMM16r4 rmat{};
			    __m512 t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
			    t0  = _mm512_mul_ps(q_x,q_x);
			    t1  = _mm512_add_ps(_mm512_add_ps(q_x,q_y),q_w);
			    t2  = _mm512_mul_ps(t1,t1);
			    t3  = _mm512_sub_ps(t0,t2); //qq
			    t5  = _mm512_mul_ps(q_y,q_z);
			    t6  = _mm512_mul_ps(q_z,q_w);
			    t7  = _mm512_mul_ps(q_w,q_y);
			    t8  = _mm512_mul_ps(q_x,q_w);
			    t9  = _mm512_mul_ps(q_x,q_y);
			    t10 = _mm512_mul_ps(q_x,q_z);
			    t4 = _mm512_mul_ps(q_y,q_y);
                            rmat.m_vRow1 = _mm512_fmadd_ps(v16f32_2,t4,t3);
			    rmat.m_vRow2 = _mm512_mul_ps(v16f32_2,_mm512_sub_ps(t5,t8));
			    rmat.m_vRow3 = _mm512_mul_ps(v16f32_2,_mm512_add_ps(t7,t10));
			    rmat.m_vRow4 = _mm512_mul_ps(v16f32_2,_mm512_add_ps(t5,t8));
			    rmat.m_vRow5 = _mm512_fmadd_ps(v16f32_2,_mm512_mul_ps(q_z,q_z),t3);
			    rmat.m_vRow6 = _mm512_mul_ps(v16f32_2,_mm512_sub_ps(t6,t9));
			    rmat.m_vRow7 = _mm512_mul_ps(v16f32_2,_mm512_sub_ps(t7,t8));
			    rmat.m_vRow8 = _mm512_mul_ps(v16f32_2,_mm512_add_ps(t6,t9));
			    rmat.m_vRow9 = _mm512_fmadd_ps(v16f32_2,_mm512_mul_ps(q_w,q_w),t3);
			    return (rmat);
		      }
		      
     } // math

} // gms













#endif /* __GMS_ROTATION_KERNELS_HPP__*/
