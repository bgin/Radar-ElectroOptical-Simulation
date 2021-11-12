
#ifndef __GMS_ROTATION_KERNELS_AVX512_HPP__
#define __GMS_ROTATION_KERNELS_AVX512_HPP__


namespace file_info {

const unsigned int gGMS_ROTATION_KERNELS_AVX512_MAJOR = 1U;
const unsigned int gGMS_ROTATION_KERNELS_AVX512_MINOR = 0U;
const unsigned int gGMS_ROTATION_KERNELS_AVX512_MICRO = 0U;
const unsigned int gGMS_ROTATION_KERNELS_AVX512_FULLVER =
       1000U*gGMS_ROTATION_KERNELS_AVX512_MAJOR+
       100U*gGMS_ROTATION_KERNELS_AVX512_MINOR +
       10U*gGMS_ROTATION_KERNELS_AVX512_MICRO;
const char * const pgGMS_ROTATION_KERNELS_AVX512_CREATION_DATE = "09-08-2021 02:25 PM +00200 (SUN 09 AUG 2021 GMT+2)";
const char * const pgGMS_ROTATION_KERNELS_AVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
const char * const pgGMS_ROTATION_KERNELS_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
const char * const pgGMS_ROTATION_KERNELS_AVX512_DESCRIPTION   = "AVX512 vectorized basic rotation operations.";
}


#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_dcm_avx512.hpp"


namespace gms {

         namespace math {

                       /*
                                 This version is *loosely based on the Fortran 90 "rotation.f90" source code
                                 implementation.
                                 
                                 *Many optimizations were applied (precomputation of common subexpression,
                                 constants folding,  AVX512 vectorization.)
                                 The original authors copyright statement
                                 Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
                                 Modified      2017-2020, Martin Diehl/Max-Planck-Institut für Eisenforschung GmbH
                                 All rights reserved.
                         */


		    namespace {

		               const __m512 v16_0      = _mm512_set1_ps(0.0F);
		               const __m512 v16_2      = _mm512_set1_ps(2.0f);
                               const __m512 v16_n1     = _mm512_set1_ps(-1.0f);
			       const __m512 v16_n2     = _mm512_set1_ps(-2.0F);
			       const __m512 v16_spi    = _mm512_set1_ps(1.7724538509055160272982F);
			       const __m512 v16_s6pi   = _mm512_set1_ps(1.381976597885341917061F);
			       const __m512 v16_a      = _mm512_set1_ps(1.9257490199582527754939F);
			       const __m512 v16_ap     = _mm512_set1_ps(2.1450293971110256000775F);
			       const __m512 v16_sc     = _mm512_set1_ps(0.8977727869612861128953F);
			       const __m512 v16_beta   = _mm512_set1_ps(0.9628745099791263877469F);
			       const __m512 v16_r1     = _mm512_set1_ps(1.3306700394914687909256F);
			       const __m512 v16_r2     = _mm512_set1_ps(1.4142135623730950488017F);
			       const __m512 v16_pi12   = _mm512_set1_ps(0.2617993877991494365386F);
			       const __m512 v16_prek   = _mm512_set1_ps(1.6434564029725030125017F);
			       const __m512 v16_pi     = _mm512_set1_pd(3.1415926535897932384626F);
			       const __m512 v16_2pi    = _mm512_set1_pd(6.2831853071795864769253F);
			       const __m512 v16_0      = _mm512_set1_pd(0.0);
		               const __m512 v16_2      = _mm512_set1_pd(2.0);
                               const __m512 v16_n1     = _mm512_set1_pd(-1.0);
			       const __m512 v16_n2     = _mm512_set1_pd(-2.0)
			       const __m512d v8_n1     = _mm512_set1_pd(-1.0);
			       const __m512  v8_pi     = _mm512_set1_pd(3.1415926535897932384626);
			       const __m512  v8_2pi    = _mm512_set1_pd(6.2831853071795864769253);
			       const __m512d v8_spi    = _mm512_set1_pd(1.7724538509055160272982);
			       const __m512d v8_s6pi   = _mm512_set1_pd(1.381976597885341917061);
			       const __m512d v8_a      = _mm512_set1_pd(1.9257490199582527754939);
			       const __m512d v8_ap     = _mm512_set1_pd(2.1450293971110256000775);
			       const __m512d v8_sc     = _mm512_set1_pd(0.8977727869612861128953);
			       const __m512d v8_beta   = _mm512_set1_pd(0.9628745099791263877469);
			       const __m512d v8_r1     = _mm512_set1_pd(1.3306700394914687909256);
			       const __m512d v8_r2     = _mm512_set1_pd(1.4142135623730950488017);
			       const __m512d v8_pi12   = _mm512_set1_pd(0.2617993877991494365386);
			       const __m512d v8f_prek   = _mm512_set1_pd(1.6434564029725030125017);
			     
	       }


	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 fmod_zmm16r4(const __m512 a,
		                          const __m512 b) {

                          __m512 v = _mm512_sub_ps(a,_mm512_mul_ps(
			             _mm512_div_round_ps(a,b,_MM_FROUND_TO_ZERO|_MM_FROUND_NO_EXEC),b));
			  return (v);
			  
		      }

		      
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d fmod_zmm8r8(const __m512d a,
		                          const __m512d b) {

                          __m512d v = _mm512_sub_pd(a,_mm512_mul_pd(
			             _mm512_div_round_pd(a,b,_MM_FROUND_TO_ZERO|_MM_FROUND_NO_EXEC),b));
			  return (v);
			  
		      }

	              __ATTR_REGCALL__
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
                            rmat.m_vRow1 = _mm512_fmadd_ps(v16_2,t4,t3);
			    rmat.m_vRow2 = _mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8));
			    rmat.m_vRow3 = _mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10));
			    rmat.m_vRow4 = _mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8));
			    rmat.m_vRow5 = _mm512_fmadd_ps(v16_2,_mm512_mul_ps(q_z,q_z),t3);
			    rmat.m_vRow6 = _mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9));
			    rmat.m_vRow7 = _mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8));
			    rmat.m_vRow8 = _mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9));
			    rmat.m_vRow9 = _mm512_fmadd_ps(v16_2,_mm512_mul_ps(q_w,q_w),t3);
			    return (rmat);
		      }


	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      DCMatZMM8r4
		      q4x8_to_rmat9x8_zmm8r8(const __m512d q_x,
		                             const __m512d q_y,
					     const __m512d q_z,
					     const __m512d q_w) {
                         
                            DCMatZMM8r8 rmat{};
			    __m512d t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
			    t0  = _mm512_mul_pd(q_x,q_x);
			    t1  = _mm512_add_pd(_mm512_add_pd(q_x,q_y),q_w);
			    t2  = _mm512_mul_pd(t1,t1);
			    t3  = _mm512_sub_pd(t0,t2); //qq
			    t5  = _mm512_mul_pd(q_y,q_z);
			    t6  = _mm512_mul_pd(q_z,q_w);
			    t7  = _mm512_mul_pd(q_w,q_y);
			    t8  = _mm512_mul_pd(q_x,q_w);
			    t9  = _mm512_mul_pd(q_x,q_y);
			    t10 = _mm512_mul_pd(q_x,q_z);
			    t4 = _mm512_mul_pd(q_y,q_y);
                            rmat.m_vRow1 = _mm512_fmadd_pd(v16f32_2,t4,t3);
			    rmat.m_vRow2 = _mm512_mul_pd(v16f32_2,_mm512_sub_pd(t5,t8));
			    rmat.m_vRow3 = _mm512_mul_pd(v16f32_2,_mm512_add_pd(t7,t10));
			    rmat.m_vRow4 = _mm512_mul_pd(v16f32_2,_mm512_add_pd(t5,t8));
			    rmat.m_vRow5 = _mm512_fmadd_pd(v16f32_2,_mm512_mul_pd(q_z,q_z),t3);
			    rmat.m_vRow6 = _mm512_mul_pd(v16f32_2,_mm512_sub_pd(t6,t9));
			    rmat.m_vRow7 = _mm512_mul_pd(v16f32_2,_mm512_sub_pd(t7,t8));
			    rmat.m_vRow8 = _mm512_mul_pd(v16f32_2,_mm512_add_pd(t6,t9));
			    rmat.m_vRow9 = _mm512_fmadd_pd(v16f32_2,_mm512_mul_pd(q_w,q_w),t3);
			    return (rmat);
		    }


		    /*
                            Convert unit quaternion to Euler angles
                      */

		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
	              void
		      q4x16_to_ea3x16_zmm16r4(const __m512 q_x,
		                              const __m512 q_y,
					      const __m512 q_z,
					      const __m512 q_w,
					      __m512 & alpha,
					      __m512 & beta,
					      __m512 & gamma) {

                         const __m512 qxw = _mm512_fmadd_ps(q_x,q_x,_mm512_mul_ps(q_w,q_w));
			 const __m512 qyz = _mm512_fmadd_ps(q_y,q_y,_mm512_mul_ps(q_z,q_z));
			 const __m512 chi = _mm512_sqrt_ps(_mm512_mul_ps(qxw,qyz));
			 __m512 alpha_c;
			 __m512 beta_c;
			 __m512 gamma_c;
			 __mmask16 qyz0 = 0x0;
			 __mmask16 qxw0 = 0x0;
			 __mmask16 k1   = 0x0;
			 __mmask16 k2   = 0x0;
			 __mmask16 k3   = 0x0;
			 qyz0 = _mm512_cmp_ps_mask(qyz0,v16_0,_CMP_EQ_OQ);
			 qxw0 = _mm512_cmp_ps_mask(qxw0,v16_0,_CMP_EQ_OQ);
			 if(1==qyz0) {
			    alpha = _mm512_atan2_ps(_mm512_mul_ps(v16_n2,_mm512_mul_ps(q_x,q_w)),
			                              _mm512_fmsub_ps(q_x,q_x,_mm512_mul_ps(q_w,q_w)));
			    beta  = v16_0;
			    gamma = v16_0;
			 } else if(1==qxw0) {
                            alpha = _mm512_atan2_ps(_mm512_mul_ps(v16_2,_mm512_mul_ps(q_y,q_z)),
			                              _mm512_fmsub_ps(q_y,q_y,_mm512_mul_ps(q_z,q_z)));
			    beta = v16_pi;
			    gamma = v16_0;
			 }
			 else {
			    const __m512 t0 = _mm512_mul_ps(q_y,q_w);
			    const __m512 c0 = _mm512_fmadd_ps(q_x,q_z,t0);
			    const __m512 t1 = _mm512_mul_ps(q_z,q_w);
			    const __m512 c1 = _mm512_fmsub_ps(q_x,q_y,t1);
			 
			    alpha = _mm512_atan2_ps(_mm512_mul_ps(chi,_mm512_mul_ps(c0,v16_n1)),
			                            _mm512_mul_ps(chi,_mm512_mul_ps(c1,v16_n1)));
			    beta  = _mm512_atan2_ps(_mm512_mul_ps(v16_2,chi),_mm512_sub_ps(qxw,qyz));
			    const __m512 c2 = _mm512_fmadd_ps(q_x,_qy,t1);
			    gamma = _mm512_atan2_ps(_mm512_mul_ps(chi,_mm512_mul_ps(c0,v16_1)),
			                            _mm512_mul_ps(chi,_mm512_mul_ps(c2,v16_n1)));
			 }
			 alpha_c = alpha;
			 const __m512 tmp0 = _mm512_fmadd_ps(v16_2,v16_pi,alpha_c);
			 k1 = _mm512_cmp_ps_mask(alpha_c,v16_0,_CMP_LT_OQ);
			 alpha = _mm512_mask_mov_ps(alpha_c,k1,fmod_zmm16r4(tmp0,v16_2pi));
			 beta_c = beta;
			 const _mm512 tmp1 = _mm512_fmadd_ps(v16_2,v16_pi,beta_c);
			 k2 = _mm512_cmp_ps_mask(beta_c,v16_0,k2,_CMP_LT_OQ);
			 beta = _mm512_mask_mov_ps(beta_c,k2,fmod_zmm16r4(tmp1,v16_pi));
			 gamma_c = gamma;
			 const __m512 tmp2 = _mm512_fmadd_ps(v16_2,v16_pi,gamma_c);
			 k3 = _mm512_cmp_ps_mask(gamma_c,v16_0,k3,_CMP_LT_OQ);
			 gamma = _mm512_mask_mov_ps(gamma_c,k3,fmod_zmm16r4(tmp2,v16_2pi));
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
	              void
		      q4x8_to_ea3x8_zmm8r8(   const __m512d q_x,
		                              const __m512d q_y,
					      const __m512d q_z,
					      const __m512d q_w,
					      __m512d & alpha,
					      __m512d & beta,
					      __m512d & gamma) {

                         const __m512d qxw = _mm512_fmadd_pd(q_x,q_x,_mm512_mul_pd(q_w,q_w));
			 const __m512d qyz = _mm512_fmadd_pd(q_y,q_y,_mm512_mul_pd(q_z,q_z));
			 const __m512d chi = _mm512_sqrt_pd(_mm512_mul_pd(qxw,qyz));
			 __m512d alpha_c;
			 __m512d beta_c;
			 __m512d gamma_c;
			 __mmask8 qyz0 = 0x0;
			 __mmask8 qxw0 = 0x0;
			 __mmask8 k1   = 0x0;
			 __mmask8 k2   = 0x0;
			 __mmask8 k3   = 0x0;
			 qyz0 = _mm512_cmp_pd_mask(qyz0,v8_0,_CMP_EQ_OQ);
			 qxw0 = _mm512_cmp_pd_mask(qxw0,v8_0,_CMP_EQ_OQ);
			 if(__builtin_expect(1==qyz0,0)) {
			    alpha = _mm512_atan2_ps(_mm512_mul_pd(v8_n2,_mm512_mul_pd(q_x,q_w)),
			                              _mm512_fmsub_pd(q_x,q_x,_mm512_mul_pd(q_w,q_w)));
			    beta  = v8_0;
			    gamma = v8_0;
			 } else if(__builtin_expect(1==qxw0,0)) {
                            alpha = _mm512_atan2_pd(_mm512_mul_pd(v16_2,_mm512_mul_pd(q_y,q_z)),
			                              _mm512_fmsub_pd(q_y,q_y,_mm512_mul_pd(q_z,q_z)));
			    beta = v8_pi;
			    gamma = v8_0;
			 }
			 else {
			    const __m512d t0 = _mm512_mul_pd(q_y,q_w);
			    const __m512d c0 = _mm512_fmadd_pd(q_x,q_z,t0);
			    const __m512d t1 = _mm512_mul_pd(q_z,q_w);
			    const __m512d c1 = _mm512_fmsub_pd(q_x,q_y,t1);
			 
			    alpha = _mm512_atan2_pd(_mm512_mul_pd(chi,_mm512_mul_pd(c0,v8_n1)),
			                            _mm512_mul_pd(chi,_mm512_mul_ps(c1,v8_n1)));
			    beta  = _mm512_atan2_pd(_mm512_mul_pd(v8_2,chi),_mm512_sub_pd(qxw,qyz));
			    const __m512 c2 = _mm512_fmadd_pd(q_x,_qy,t1);
			    gamma = _mm512_atan2_pd(_mm512_mul_pd(chi,_mm512_mul_pd(c0,v8_1)),
			                            _mm512_mul_pd(chi,_mm512_mul_pd(c2,v8_n1)));
			 }
			 alpha_c = alpha;
			 const __m512d tmp0 = _mm512_fmadd_pd(v8_2,v16_pi,alpha_c);
			 k1 = _mm512_cmp_pd_mask(alpha_c,v8_0,_CMP_LT_OQ);
			 alpha = _mm512_mask_mov_ps(alpha_c,k1,fmod_zmm8r8(tmp0,v8_2pi));
			 beta_c = beta;
			 const _mm512d tmp1 = _mm512_fmadd_ps(v8_2,v8_pi,beta_c);
			 k2 = _mm512_cmp_pd_mask(beta_c,v8_0,k2,_CMP_LT_OQ);
			 beta = _mm512_mask_mov_pd(beta_c,k2,fmod_zmm8r8(tmp1,v8_pi));
			 gamma_c = gamma;
			 const __m512d tmp2 = _mm512_fmadd_pd(v8_2,v8_pi,gamma_c);
			 k3 = _mm512_cmp_pd_mask(gamma_c,v8_0,k3,_CMP_LT_OQ);
			 gamma = _mm512_mask_mov_pd(gamma_c,k3,fmod_zmm8r8(tmp2,v8_2pi));
		    }
		      
     } // math

} // gms













#endif /* __GMS_ROTATION_KERNELS_AVX512_HPP__*/
