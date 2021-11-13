
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
			       const __m512 v16_1o2    = _mm512_set1_ps(0.5F);
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
			       const __m512 v16_pi     = _mm512_set1_ps(3.1415926535897932384626F);
			       const __m512 v16_2pi    = _mm512_set1_ps(6.2831853071795864769253F);
			      
			       const __m512d v8_1o2    = _mm512_set1_pd(0.5);
                               const __m512d v8_0      = _mm512_set1_pd(0.0);
		               const __m512d v8_2      = _mm512_set1_pd(2.0);
			       const __m512d v8_n1     = _mm512_set1_pd(-1.0);
			       const __m512d  v8_pi     = _mm512_set1_pd(3.1415926535897932384626);
			       const __m512d  v8_2pi    = _mm512_set1_pd(6.2831853071795864769253);
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


			       

                                    __ATTR_ALWAYS_INLINE__
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                            static inline
				    __m512d
				    zmm8r8_sign_zmm8r8(const __m512d va,
				                       const __m512d vb) {
				       
				       register __m512d vret = _0;
				       register __m512d t0   = _mm512_abs_pd(va);
                                       __mmask8 gez = 0x0;
				       gez  = _mm512_cmp_pd_mask(vb,v8_0,_CMP_GE_OQ); // Lat=3refc,Thr=1refc
				       vret = _mm512_mask_blend_pd(gez,t0,_mm512_sub_pd(v8_0,t0)); //Lat=1refc,Thr=0.5refc,Lat=4refc,Thr=1refc
				       return (vret);
				                                       
				    }

				    __ATTR_ALWAYS_INLINE__
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                            static inline
				    __m512
				    zmm16r4_sign_zmm16r4(const __m512 va,
				                         const __m512 vb) {
				       
				       register __m512 vret = _0;
				       register __m512 t0   = _mm512_abs_ps(va);
                                       __mmask8 gez = 0x0;
				       gez  = _mm512_cmp_ps_mask(vb,v16_0,_CMP_GE_OQ); // Lat=3refc,Thr=1refc
				       vret = _mm512_mask_blend_ps(gez,t0,_mm512_sub_ps(v16_0,t0)); //Lat=1refc,Thr=0.5refc,Lat=4refc,Thr=1refc
				       return (vret);
				                                       
				    }
			    
			     
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
		      __m512d norm2_zmm8r8(const __m512d y,
		                          const __m512d z,
					  const __m512d w) {

                            const __m512d t0 = _mm512_mul_pd(y,y);
			    const __m512d t1 = _mm512_mul_pd(z,z);
			    const __m512d t2 = _mm512_mul_pd(w,w);
			    const __m512d v  = _mm512_add_pd(t0,_mm512_add_pd(t1,t2));
			    return (_mm512_sqrt_pd(v));
			    
		      }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d norm2_zmm16r4(const __m512 y,
		                            const __m512 z,
					    const __m512 w) {

                            const __m512 t0 = _mm512_mul_ps(y,y);
			    const __m512 t1 = _mm512_mul_ps(z,z);
			    const __m512 t2 = _mm512_mul_ps(w,w);
			    const __m512 v  = _mm512_add_ps(t0,_mm512_add_ps(t1,t2));
			    return (_mm512_sqrt_ps(v));
			    
		      }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 clip_zmm16r4(const __m512 x,
		                          const __m512 lo,
					  const __m512 hi) {

                           return (_mm512_max_ps(lo,_mm512_min_ps(x,hi)));
		      }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d clip_zmm8r8(const __m512d x,
		                          const __m512d lo,
					  const __m512d hi) {

                           return (_mm512_max_pd(lo,_mm512_min_pd(x,hi)));
		      }


		      

	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      DCM9x16
		      q4x16_to_rmat9x16_zmm16r4(const __m512 q_x,
		                                const __m512 q_y,
						const __m512 q_z,
						const __m512 q_w) {
                            DCM9x16 rmat{};
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
                      DCM9x8
		      q4x8_to_rmat9x8_zmm8r8(const __m512d q_x,
		                             const __m512d q_y,
					     const __m512d q_z,
					     const __m512d q_w) {
                         
                            DCM9x8 rmat{};
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



		   /*
                       Convert unit quaternion to axis angle pair
                    */


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void
		      q4x16_to_ax4x16_zmm16r4(const __m512 q_x,
		                              const __m512 q_y,
					      const __m512 q_z,
					      const __m512 q_w,
					      __m512 & ax_1,
					      __m512 & ax_2,
					      __m512 & ax_3,
					      __m512 & ax_4) {

                          const register __m512 t0 = _mm512_mul_ps(q_y,q_y);
			  const register __m512 t1 = _mm512_mul_ps(q_z,q_z);
			  const register __m512 t2 = _mm512_mul_ps(q_w,q_w);
			  const register __m512 v0 = _mm512_add_ps(t0,_mm512_add_ps(t1,t2));
			  __m512 c0 = v16_0;
			  __m512 c1 = v16_0;
			  __mmask16 k1 = 0x0;
			  __mmask16 k2 = 0x0;
			  k1 = _mm512_cmp_ps_mask(v0,v16_0,_CMP_EQ_OQ);
			  if(__builtin_expect(1==k1,0)) {
                             ax_1 = v16_0;
			     ax_2 = v16_0;
			     ax_3 = v16_1;
			     ax_4 = v16_0;
			  }
			  k2 = _mm512_cmp_ps_mask(q_x,_v16_0,_CMP_NEQ_OQ);
			  if(1==k2) {
                             const register __m512 s = _mm512_div_ps(
			                          zmm16r4_sign_zmm16r4(v16_1,q_x),
						  norm2_zmm16r4(q_y,q_z,q_w));
			     ax_1 = _mm512_mul_ps(q_y,s);
			     ax_2 = _mm512_mul_ps(q_z,s);
			     ax_2 = _mm512_mul_ps(q_w,s);
			     const register __m512 omega = _mm512_mul_ps(v16_2,
			                              _mm512_acos_ps(clip_zmm16r4(q_x,v16_n1,v16_1)));
			     ax_4 = omega;
			  }
			  else {
                             ax_1 = q_y;
			     ax_2 = q_z;
			     ax_3 = q_w;
			     ax_4 = v16_pi;
			     
			  }
		     }



		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void
		      q4x8_to_ax4x8_zmm8r8(   const __m512d q_x,
		                              const __m512d q_y,
					      const __m512d q_z,
					      const __m512d q_w,
					      __m512d & ax_1,
					      __m512d & ax_2,
					      __m512d & ax_3,
					      __m512d & ax_4) {

                          const register __m512d t0 = _mm512_mul_pd(q_y,q_y);
			  const register __m512d t1 = _mm512_mul_pd(q_z,q_z);
			  const register __m512d t2 = _mm512_mul_pd(q_w,q_w);
			  const register __m512d v0 = _mm512_add_pd(t0,_mm512_add_pd(t1,t2));
			  __mmask8 k1 = 0x0;
			  __mmask8 k2 = 0x0;
			  k1 = _mm512_cmp_pd_mask(v0,v8_0,_CMP_EQ_OQ);
			  if(__builtin_expect(1==k1,0)) {
                             ax_1 = v8_0;
			     ax_2 = v8_0;
			     ax_3 = v8_1;
			     ax_4 = v8_0;
			  }
			  k2 = _mm512_cmp_pd_mask(q_x,_v8_0,_CMP_NEQ_OQ);
			  if(1==k2) {
                             const register __m512 s = _mm512_div_pd(
			                          zmm8r8_sign_zmm8r8(v8_1,q_x),
						  norm2_zmm8r8(q_y,q_z,q_w));
			     ax_1 = _mm512_mul_pd(q_y,s);
			     ax_2 = _mm512_mul_pd(q_z,s);
			     ax_2 = _mm512_mul_pd(q_w,s);
			     const register __m512 omega = _mm512_mul_pd(v16_2,
			                              _mm512_acos_pd(clip_zmm8r8(q_x,v8_n1,v8_1)));
			     ax_4 = omega;
			  }
			  else {
                             ax_1 = q_y;
			     ax_2 = q_z;
			     ax_3 = q_w;
			     ax_4 = v8_pi;
			     
			  }
		     }


		    /*
                        Convert unit quaternion to Rodrigues vector
                     */
#include <limits>

		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void
		      q4x16_to_rv4x16_zmm16r4(const __m512 q_x,
		                              const __m512 q_y,
					      const __m512 q_z,
					      const __m512 q_w,
					      __m512 & r_x,
					      __m512 & r_y,
					      __m512 & r_z,
					      __m512 & r_w) {

                          const register __m512 thr = _mm512_set1_ps(1.0e-8);
			  const register __m512 inf = _mm512_set1_ps(std::numeric_limits<float>::infinity());
			  register __m512 t0 = v16_0;
			  register __m512 s  = v16_0;
			  __mmask16 k1 = 0x0;
			  __mmask16 k2 = 0x0;
			  k1 = _mm512_cmp_ps_mask(_mm512_abs_ps(q_x),thr,_CMP_LT_OQ);
			  if(__builtin_expect(1==k1,0)) {
                             r_x = q_y;
			     r_y = q_z;
			     r_z = q_w;
			     r_w = inf;
			  }
			  else {
                               s   = norm2_zmm16r4(q_y,q_z,q_w);
			       k2  = _mm512_cmp_ps_mask(s,thr,_CMP_LT_OQ);
			       r_x = _mm512_mask_blend_ps(k2,_mm512_div_ps(q_y,s),v16_0);
			       t0  = _mm512_acos_ps(clip_zmm16r4(q_x,v16_n1,v16_1));
			       r_y = _mm512_mask_blend_ps(k2,_mm512_div_ps(q_z,s),v16_0);
			       r_z = _mm512_mask_blend_ps(k2,_mm512_div_ps(q_w,s),v16_n1);
			       r_w = _mm512_mask_blend_ps(k2,_mm512_tan_ps(t0),v16_0);
			  }
		     }


		     
	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void
		      q4x8_to_rv4x8_zmm8r8(   const __m512d q_x,
		                              const __m512d q_y,
					      const __m512d q_z,
					      const __m512d q_w,
					      __m512d & r_x,
					      __m512d & r_y,
					      __m512d & r_z,
					      __m512d & r_w) {

                          const register __m512d thr = _mm512_set1_pd(1.0e-8);
			  const register __m512d inf = _mm512_set1_pd(std::numeric_limits<double>::infinity());
			  register __m512d t0 = v8_0;
			  register __m512d s  = v8_0;
			  __mmask8 k1 = 0x0;
			  __mmask8 k2 = 0x0;
			  k1 = _mm512_cmp_pd_mask(_mm512_abs_pd(q_x),thr,_CMP_LT_OQ);
			  if(__builtin_expect(1==k1,0)) {
                             r_x = q_y;
			     r_y = q_z;
			     r_z = q_w;
			     r_w = inf;
			  }
			  else {
                               s   = norm2_zmm8r8(q_y,q_z,q_w);
			       k2  = _mm512_cmp_pd_mask(s,thr,_CMP_LT_OQ);
			       r_x = _mm512_mask_blend_pd(k2,_mm512_div_pd(q_y,s),v8_0);
			       t0  = _mm512_acos_pd(clip_zmm8r8(q_x,v8_n1,v8_1));
			       r_y = _mm512_mask_blend_pd(k2,_mm512_div_pd(q_z,s),v8_0);
			       r_z = _mm512_mask_blend_pd(k2,_mm512_div_pd(q_w,s),v8_n1);
			       r_w = _mm512_mask_blend_pd(k2,_mm512_tan_pd(t0),v8_0);
			  }
		     }


		     /*
                           Orientation i.e. (Direct Cosine Matrix)  matrix to Euler angles.
                       */


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void
		      rmat9x16_to_ea3x16_zmm16r4(const DCM9x16 rm,
		                                 __m512 & alpha,
						 __m512 & beta,
						 __m512 & gamma) {

                           const    __m512  thr  = _mm512_set1_ps(1.0e-8);
			   register __m512  t0   = v16_0;
			   register __m512  t1   = v16_0;
			   register __m512  t2   = v16_0;
			   register __m512  zeta = v16_0;
			   register __m512  al_c = v16_0;
			   register __m512  be_c = v16_0;
			   register __m512  ga_c = v16_0;
			   
			   __mmask16 k1 = 0x0;
			   __mmask16 k2 = 0x0;
			   __mmask16 k3 = 0x0;
			   __mmask16 k4 = 0x0;
			   k1 = _mm512_cmp_ps_mask(rm.m_vRow9,v16_1,_CMP_NEQ_OQ);
			   t0 = _mm512_mul_ps(rm.m_vRow9,rm.m_vRow9);
			   zeta = _mm512_div_ps(_mm512_sqrt_ps(_mm512_sub_ps(v16_1,t0)));
			   t0 = _mm512_sub_ps(v16_0,rm.vRow8);
			   alpha = _mm512_mask_blend_ps(k1,_mm512_atan2_ps(rm.vRow2,rm.vRow1),
			                                   _mm512_atan2_ps(_mm512_mul_ps(rm.vRow7,zeta),
			   				                   _mm512_mul_ps(t0,zeta)));
			   t1    = _mm512_mul_ps(v16_1o2,_mm512_mul_ps(v16_pi,_mm512_sub_ps(v16_1,rm.vRow9)));
			   beta  = _mm512_mask_blend_ps(k1,t1,_mm512_acos_ps(rm.vRow9));
			   gamma = _mm512_mask_blend_ps(k1,v16_0,_mm512_atan_ps(_mm512_mul_ps(rm.vRow3,zeta),
			                                                        _mm512_mul_ps(rm.vRow6,zeta)));
			   al_c = alpha;
			   be_c = beta;
			   ga_c = gamma;
			   k2 = _mm512_cmp_ps_mask(_mm512_abs_ps(alpha),thr,_CMP_LT_OQ);
			   alpha = _mm512_mask_mov_ps(al_c,k2,v16_0);
			   k3 = _mm512_cmp_ps_mask(_mm512_abs_ps(beta),thr,_CMP_LT_OQ);
			   beta = _mm512_mask_mov_ps(be_c,k3,v16_0);
			   k4 = _mm512_cmp_ps_mask(_mm512_abs_ps(gamma),thr,_CMP_LT_OQ);
			   gamma = _mm512_mask_mov_ps(ga_c,k4,v16_0);
			   //al_c = alpha;
			   t0 = _mm512_add_ps(alpha,v16_2pi);
			   k2 = _mm512_cmp_ps_mask(alpha,v16_0,_CMP_LT_OQ);
			   alpha = _mm512_mask_mov_ps(al_c,k2,fmod_zmm16r4(t0,v16_2pi));
			   //be_c = beta;
			   t1 = _mm512_add_ps(beta,v16_2pi);
			   k3 = _mm512_cmp_ps_mask(beta,v16_0,_CMP_LT_OQ);
			   beta = _mm512_mask_mov_ps(be_c,k3,fmod_zmm16r4(t1,v16_pi));
			   //ga_c = gamma;
			   t2 = _mm512_add_ps(gamma,v16_2pi);
			   k4 = _mm512_cmp_ps_mask(gamma,v16_0,_CMP_LT_OQ);
			   gamma = _mm512_mask_mov_ps(ga_c,k4,fmod_zmm16r4(t2,v16_2pi));
		     }



		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void
		      rmat9x8_to_ea3x8_zmm8r8(const DCM9x8 rm,
		                                 __m512d & alpha,
						 __m512d & beta,
						 __m512d & gamma) {

                           const    __m512d  thr  = _mm512_set1_pd(1.0e-8);
			   register __m512d  t0   = v8_0;
			   register __m512d  t1   = v8_0;
			   register __m512d  t2   = v8_0;
			   register __m512d  zeta = v8_0;
			   register __m512d  al_c = v8_0;
			   register __m512d  be_c = v8_0;
			   register __m512d  ga_c = v8_0;
			   
			   __mmask8 k1 = 0x0;
			   __mmask8 k2 = 0x0;
			   __mmask8 k3 = 0x0;
			   __mmask8 k4 = 0x0;
			   k1 = _mm512_cmp_pd_mask(rm.m_vRow9,v8_1,_CMP_NEQ_OQ);
			   t0 = _mm512_mul_pd(rm.m_vRow9,rm.m_vRow9);
			   zeta = _mm512_div_pd(_mm512_sqrt_pd(_mm512_sub_pd(v8_1,t0)));
			   t0 = _mm512_sub_pd(v8_0,rm.vRow8);
			   alpha = _mm512_mask_blend_pd(k1,_mm512_atan2_pd(rm.vRow2,rm.vRow1),
			                                   _mm512_atan2_pd(_mm512_mul_pd(rm.vRow7,zeta),
			   				                   _mm512_mul_pd(t0,zeta)));
			   t1    = _mm512_mul_pd(v8_1o2,_mm512_mul_pd(v8_pi,_mm512_sub_pd(v8_1,rm.vRow9)));
			   beta  = _mm512_mask_blend_pd(k1,t1,_mm512_acos_pd(rm.vRow9));
			   gamma = _mm512_mask_blend_pd(k1,v8_0,_mm512_atan_pd(_mm512_mul_pd(rm.vRow3,zeta),
			                                                        _mm512_mul_pd(rm.vRow6,zeta)));
			   al_c = alpha;
			   be_c = beta;
			   ga_c = gamma;
			   k2 = _mm512_cmp_pd_mask(_mm512_abs_pd(alpha),thr,_CMP_LT_OQ);
			   alpha = _mm512_mask_mov_pd(al_c,k2,v8_0);
			   k3 = _mm512_cmp_pd_mask(_mm512_abs_pd(beta),thr,_CMP_LT_OQ);
			   beta = _mm512_mask_mov_pd(be_c,k3,v8_0);
			   k4 = _mm512_cmp_pd_mask(_mm512_abs_pd(gamma),thr,_CMP_LT_OQ);
			   gamma = _mm512_mask_mov_pd(ga_c,k4,v8_0);
			   //al_c = alpha;
			   t0 = _mm512_add_pd(alpha,v8_2pi);
			   k2 = _mm512_cmp_pd_mask(alpha,v8_0,_CMP_LT_OQ);
			   alpha = _mm512_mask_mov_pd(al_c,k2,fmod_zmm8r8(t0,v8_2pi));
			   //be_c = beta;
			   t1 = _mm512_add_pd(beta,v8_2pi);
			   k3 = _mm512_cmp_pd_mask(beta,v8_0,_CMP_LT_OQ);
			   beta = _mm512_mask_mov_pd(be_c,k3,fmod_zmm8r8(t1,v8_pi));
			   //ga_c = gamma;
			   t2 = _mm512_add_pd(gamma,v8_2pi);
			   k4 = _mm512_cmp_pd_mask(gamma,v8_0,_CMP_LT_OQ);
			   gamma = _mm512_mask_mov_pd(ga_c,k4,fmod_zmm8r8(t2,v8_2pi));
		     }
					
		      
     } // math

} // gms













#endif /* __GMS_ROTATION_KERNELS_AVX512_HPP__*/
