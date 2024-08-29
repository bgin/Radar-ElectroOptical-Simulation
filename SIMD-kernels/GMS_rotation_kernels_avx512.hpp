
#ifndef __GMS_ROTATION_KERNELS_AVX512_HPP__
#define __GMS_ROTATION_KERNELS_AVX512_HPP__ 090820210225


namespace file_info {

const unsigned int gGMS_ROTATION_KERNELS_AVX512_MAJOR = 1U;
const unsigned int gGMS_ROTATION_KERNELS_AVX512_MINOR = 0U;
const unsigned int gGMS_ROTATION_KERNELS_AVX512_MICRO = 1U;
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
#include "GMS_rotations_avx512_helpers.hpp"

namespace gms {

         namespace math {

                     


		  

			       

                               

				 


	             


	           


		   


		   


		      

	             
                      __ATTR_ALWAYS_INLINE__
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


	             
                      __ATTR_ALWAYS_INLINE__
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
                       Random rotation matrix 3x3 of sphere.
                       Based on Jim Arvo, 1991 implementation.
                       Original Matrix 3x3 is represented as
                       SIMD Matrix 9x16
                   */

		    
                    __ATTR_ALWAYS_INLINE__
		   static inline
		    DCM9x16
		    random_sphere_rm9x16_zmm16r4(const __m512 vr1,
		                                 const __m512 vr2,
                                                 const __m512 vr3) {

                          DCM9x16 rm;
		          const __m512 theta = _mm512_mul_ps(vr1,v16_2pi);
			  const __m512 phi   = _mm512_mul_ps(vr2,v16_2pi);
			  const __m512 z     = _mm512_add_ps(vr3,vr3);
			  rm.m_vRow9         = _mm512_sub_ps(v16_1,z);
			  const __m512 r     = _mm512_sqrt_ps(z);
			  const __m512 vx    = _mm512_mul_ps(r,_mm512_sin_ps(phi));
			  const __m512 vy    = _mm512_mul_ps(r,_mm512_cos_ps(phi));
			  const __m512 vz    = _mm512_sqrt_ps(_mm512_sub_ps(v16_2,z));
			  rm.m_vRow6         = _mm512_mul_ps(vy,vz);
			  rm.m_vRow3         = _mm512_mul_ps(vx,vz);
			  const __m512 st    = _mm512_sin_ps(theta);
			  const __m512 ct    = _mm512_cos_ps(theta);
			  const __m512 sx    = _mm512_fmsub_ps(vx,ct,_mm512_mul_ps(vy,st));
			  rm.m_vRow7         = _mm512_mul_ps(vz,sx); 
			  rm.m_vRow1         = _mm512_fmsub_ps(vx,sx,ct);
			  const __m512 sy    = _mm512_fmadd_ps(vx,st,_mm512_mul_ps(vy,ct));
			  rm.m_vRow8         = _mm512_mul_ps(vz,sy);
			  rm.m_vRow2         = _mm512_fmsub_ps(vx,sy,st);
			  rm.m_vRow4         = _mm512_fmadd_ps(vy,sx,st);
			  rm.m_vRow5         = _mm512_fmsub_ps(vy,sy,ct);
			  return (rm);
		    }


		      /*
                       Random rotation matrix 3x3 of sphere.
                       Based on Jim Arvo, 1991 implementation.
                       Original Matrix 3x3 is represented as
                       SIMD Matrix 9x8
                   */

		    
                    __ATTR_ALWAYS_INLINE__
		   static inline
		    DCM9x8
		    random_sphere_rm9x8_zmm8r8(  const __m512d vr1,
		                                 const __m512d vr2,
                                                 const __m512d vr3) {

                          DCM9x8 rm;
		          const __m512d theta = _mm512_mul_pd(vr1,v8_2pi);
			  const __m512d phi   = _mm512_mul_pd(vr2,v8_2pi);
			  const __m512d z     = _mm512_add_pd(vr3,vr3);
			  rm.m_vRow9          = _mm512_sub_pd(v8_1,z);
			  const __m512d r     = _mm512_sqrt_pd(z);
			  const __m512d vx    = _mm512_sin_pd(phi,r);
			  const __m512d vy    = _mm512_cos_pd(phi,r);
			  const __m512d vz    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z));
			  rm.m_vRow6          = _mm512_mul_pd(vy,vz);
			  rm.m_vRow3          = _mm512_mul_pd(vx,vz);
			  const __m512d st    = _mm512_sin_pd(theta);
			  const __m512d ct    = _mm512_cos_pd(theta);
			  const __m512d sx    = _mm512_fmsub_pd(vx,ct,_mm512_mul_pd(vy,st));
			  rm.m_vRow7          = _mm512_mul_pd(vz,sx); 
			  rm.m_vRow1          = _mm512_fmsub_pd(vx,sx,ct);
			  const __m512d sy    = _mm512_fmadd_pd(vx,st,_mm512_mul_pd(vy,ct));
			  rm.m_vRow8          = _mm512_mul_pd(vz,sy);
			  rm.m_vRow2          = _mm512_fmsub_pd(vx,sy,st);
			  rm.m_vRow4          = _mm512_fmadd_pd(vy,sx,st);
			  rm.m_vRow5          = _mm512_fmsub_pd(vy,sy,ct);
			  return (rm);
		    }


		    
                            /*  This algorithm generates a gaussian deviate for each coordinate, so
                                *  the total effect is to generate a symmetric 4-D gaussian distribution,
                                *  by separability. Projecting onto the surface of the hypersphere gives
                                *  a uniform distribution.
                                Based on  Ken Shoemake, September 1991 implementation.
                                Manually vectorized.
                            */

		      
                      __ATTR_ALWAYS_INLINE__
		     static inline
	              void
		      urand_q4x16_a_zmm16r4(const __m512 vrx,   // random gaussian vector uniformly distributed [0,1]
		                            const __m512 vry,   // random gaussian vector uniformly distributed [0,1]
					    const __m512 vrz,   // random gaussian vector uniformly distributed [0,1]
					    const __m512 vrw,   // random gaussian vector uniformly distributed [0,1]
		                            float * __restrict __ATTR_ALIGN__(64) q_x,
					    float * __restrict __ATTR_ALIGN__(64) q_y,
					    float * __restrict __ATTR_ALIGN__(64) q_z,
					    float * __restrict __ATTR_ALIGN__(64) q_w) {

		         const register __m512 s1    = _mm512_fmadd_ps(vrx,vrx,_mm512_mul_ps(vry,vry));
			 const register __m512 num1  = _mm512_mul_ps(v16_n2,_mm512_log_ps(s1));
			 const register __m512 s2    = _mm512_fmadd_ps(vrz,vrz,_mm512_mul_ps(vrw,vrw));
			 const register __m512 num2  = _mm512_mul_ps(v16_n2,_mm512_log_ps(s2));
			 const register __m512 r     = _mm512_add_ps(num1,num2);
			 const register __m512 invr  = _mm512_div_ps(v16_1,r);
			 const register __m512 root1 = _mm512_sqrt_ps(_mm512_mul_ps(invr,_mm512_div_ps(num1,s1)));
			 _mm512_store_ps(&q_x[0], _mm512_mul_ps(vrx,root1));
			 _mm512_store_ps(&q_y[0], _mm512_mul_ps(vry,root1));
			 const register __m512 root2 = _mm512_sqrt_ps(_mm512_mul_ps(invr,_mm512_div_ps(num2,s2)));
			 _mm512_store_ps(&q_z[0],  _mm512_mul_ps(vrz,root2));			
			 _mm512_store_ps(&q_w[0],  _mm512_mul_ps(vrw,root2));
		   }


		    
                      __ATTR_ALWAYS_INLINE__
		      static inline
	              void
		      urand_q4x16_u_zmm16r4(const __m512 vrx,   // random gaussian vector uniformly distributed [0,1]
		                            const __m512 vry,   // random gaussian vector uniformly distributed [0,1]
					    const __m512 vrz,   // random gaussian vector uniformly distributed [0,1]
					    const __m512 vrw,   // random gaussian vector uniformly distributed [0,1]
		                            float * __restrict q_x,
					    float * __restrict q_y,
					    float * __restrict q_z,
					    float * __restrict q_w) {

		         const register __m512 s1    = _mm512_fmadd_ps(vrx,vrx,_mm512_mul_ps(vry,vry));
			 const register __m512 num1  = _mm512_mul_ps(v16_n2,_mm512_log_ps(s1));
			 const register __m512 s2    = _mm512_fmadd_ps(vrz,vrz,_mm512_mul_ps(vrw,vrw));
			 const register __m512 num2  = _mm512_mul_ps(v16_n2,_mm512_log_ps(s2));
			 const register __m512 r     = _mm512_add_ps(num1,num2);
			 const register __m512 invr  = _mm512_div_ps(v16_1,r);
			 const register __m512 root1 = _mm512_sqrt_ps(_mm512_mul_ps(invr,_mm512_div_ps(num1,s1)));
			 _mm512_storeu_ps(&q_x[0], _mm512_mul_ps(vrx,root1));
			 _mm512_storeu_ps(&q_y[0], _mm512_mul_ps(vry,root1));
			 const register __m512 root2 = _mm512_sqrt_ps(_mm512_mul_ps(invr,_mm512_div_ps(num2,s2)));
			 _mm512_storeu_ps(&q_z[0],  _mm512_mul_ps(vrz,root2));			
			 _mm512_storeu_ps(&q_w[0],  _mm512_mul_ps(vrw,root2));
		   }



		    
                      __ATTR_ALWAYS_INLINE__
		      static inline
	              void
		      urand_q4x8_a_zmm8r8(  const __m512d vrx,   // random gaussian vector uniformly distributed [0,1]
		                          const __m512d vry,   // random gaussian vector uniformly distributed [0,1]
					  const __m512d vrz,   // random gaussian vector uniformly distributed [0,1]
					  const __m512d vrw,   // random gaussian vector uniformly distributed [0,1]
		                          double * __restrict  __ATTR_ALIGN__(64) q_x,
					  double * __restrict  __ATTR_ALIGN__(64) q_y,
					  double * __restrict  __ATTR_ALIGN__(64) q_z,
					  double * __restrict  __ATTR_ALIGN__(64) q_w) {

		         const register __m512d s1    = _mm512_fmadd_pd(vrx,vrx,_mm512_mul_pd(vry,vry));
			 const register __m512d num1  = _mm512_mul_pd(v8_n2,_mm512_log_pd(s1));
			 const register __m512d s2    = _mm512_fmadd_pd(vrz,vrz,_mm512_mul_pd(vrw,vrw));
			 const register __m512d num2  = _mm512_mul_pd(v8_n2,_mm512_log_pd(s2));
			 const register __m512d r     = _mm512_add_pd(num1,num2);
			 const register __m512d invr  = _mm512_div_pd(v8_1,r);
			 const register __m512d root1 = _mm512_sqrt_pd(_mm512_mul_ps(invr,_mm512_div_pd(num1,s1)));
			 _mm512_store_pd(&q_x[0], _mm512_mul_pd(vrx,root1));
			 _mm512_store_pd(&q_y[0], _mm512_mul_pd(vry,root1));
			 const register __m512d root2 = _mm512_sqrt_pd(_mm512_mul_pd(invr,_mm512_div_pd(num2,s2)));
			 _mm512_store_pd(&q_z[0],  _mm512_mul_pd(vrz,root2));			
			 _mm512_store_pd(&q_w[0],  _mm512_mul_pd(vrw,root2));
		   }


		    
                      __ATTR_ALWAYS_INLINE__
		     static inline
	              void
		      urand_q4x8_u_zmm8r8(  const __m512d vrx,   // random gaussian vector uniformly distributed [0,1]
		                          const __m512d vry,   // random gaussian vector uniformly distributed [0,1]
					  const __m512d vrz,   // random gaussian vector uniformly distributed [0,1]
					  const __m512d vrw,   // random gaussian vector uniformly distributed [0,1]
		                          double * __restrict   q_x,
					  double * __restrict   q_y,
					  double * __restrict   q_z,
					  double * __restrict   q_w) {

		         const register __m512d s1    = _mm512_fmadd_pd(vrx,vrx,_mm512_mul_pd(vry,vry));
			 const register __m512d num1  = _mm512_mul_pd(v8_n2,_mm512_log_pd(s1));
			 const register __m512d s2    = _mm512_fmadd_pd(vrz,vrz,_mm512_mul_pd(vrw,vrw));
			 const register __m512d num2  = _mm512_mul_pd(v8_n2,_mm512_log_pd(s2));
			 const register __m512d r     = _mm512_add_pd(num1,num2);
			 const register __m512d invr  = _mm512_div_pd(v8_1,r);
			 const register __m512d root1 = _mm512_sqrt_pd(_mm512_mul_ps(invr,_mm512_div_pd(num1,s1)));
			 _mm512_storeu_pd(&q_x[0], _mm512_mul_pd(vrx,root1));
			 _mm512_storeu_pd(&q_y[0], _mm512_mul_pd(vry,root1));
			 const register __m512d root2 = _mm512_sqrt_pd(_mm512_mul_pd(invr,_mm512_div_pd(num2,s2)));
			 _mm512_storeu_pd(&q_z[0],  _mm512_mul_pd(vrz,root2));			
			 _mm512_storeu_pd(&q_w[0],  _mm512_mul_pd(vrw,root2));
		   }

	   

		    /*
                            Convert unit quaternion to Euler angles
                      */

		     
                      __ATTR_ALWAYS_INLINE__
		     static inline
	              void
		      q4x16_to_ea3x16_a_zmm16r4(const __m512 q_x,
		                              const __m512 q_y,
					      const __m512 q_z,
					      const __m512 q_w,
					      float * __restrict __ATTR_ALIGN__(64) alpha,
					      float * __restrict __ATTR_ALIGN__(64) beta,
					      float * __restrict __ATTR_ALIGN__(64) gamma) {

                         const __m512 qxw = _mm512_fmadd_ps(q_x,q_x,_mm512_mul_ps(q_w,q_w));
			 const __m512 qyz = _mm512_fmadd_ps(q_y,q_y,_mm512_mul_ps(q_z,q_z));
			 const __m512 chi = _mm512_sqrt_ps(_mm512_mul_ps(qxw,qyz));
			 __m512 alpha_c;
			 __m512 beta_c;
			 __m512 gamma_c;
			 __m512 alp,bet,gam;
			 __mmask16 qyz0 = 0x0;
			 __mmask16 qxw0 = 0x0;
			 __mmask16 k1   = 0x0;
			 __mmask16 k2   = 0x0;
			 __mmask16 k3   = 0x0;
			 qyz0 = _mm512_cmp_ps_mask(qyz,v16_0,_CMP_EQ_OQ);
			 qxw0 = _mm512_cmp_ps_mask(qxw,v16_0,_CMP_EQ_OQ);
			 if(1==qyz0) {
			    _mm512_store_ps(&alpha[0], _mm512_atan2_ps(_mm512_mul_ps(v16_n2,_mm512_mul_ps(q_x,q_w)),
			                              _mm512_fmsub_ps(q_x,q_x,_mm512_mul_ps(q_w,q_w))));
			    _mm512_store_ps(&beta[0],   v16_0);
			    _mm512_store_ps(&gamma[0],  v16_0);
			   // return;
			 } else if(1==qxw0) {
                            _mm512_store_ps(&alpha[0], _mm512_atan2_ps(_mm512_mul_ps(v16_2,_mm512_mul_ps(q_y,q_z)),
			                              _mm512_fmsub_ps(q_y,q_y,_mm512_mul_ps(q_z,q_z))));
			    _mm512_store_ps(&beta[0], v16_pi);
			    _mm512_store_ps(&gama[0], v16_0);
			    //return;
			 }
			 else {
			    const __m512 t0 = _mm512_mul_ps(q_y,q_w);
			    const __m512 c0 = _mm512_fmadd_ps(q_x,q_z,t0);
			    const __m512 t1 = _mm512_mul_ps(q_z,q_w);
			    const __m512 c1 = _mm512_fmsub_ps(q_x,q_y,t1);
			 
			    alp = _mm512_atan2_ps(_mm512_mul_ps(chi,_mm512_mul_ps(c0,v16_n1)),
			                            _mm512_mul_ps(chi,_mm512_mul_ps(c1,v16_n1)));
			    bet  = _mm512_atan2_ps(_mm512_mul_ps(v16_2,chi),_mm512_sub_ps(qxw,qyz));
			    const __m512 c2 = _mm512_fmadd_ps(q_x,_qy,t1);
			    gam = _mm512_atan2_ps(_mm512_mul_ps(chi,_mm512_mul_ps(c0,v16_1)),
			                            _mm512_mul_ps(chi,_mm512_mul_ps(c2,v16_n1)));
			 }
			 alpha_c = alp;
			 const __m512 tmp0 = _mm512_fmadd_ps(v16_2,v16_pi,alpha_c);
			 k1 = _mm512_cmp_ps_mask(alpha_c,v16_0,_CMP_LT_OQ);
			 _mm512_store_ps(&alpha[0],_mm512_mask_mov_ps(alpha_c,k1,fmod_zmm16r4(tmp0,v16_2pi)));
			 beta_c = bet;
			 const _mm512 tmp1 = _mm512_fmadd_ps(v16_2,v16_pi,beta_c);
			 k2 = _mm512_cmp_ps_mask(beta_c,v16_0,k2,_CMP_LT_OQ);
			 _mm512_store_ps(&beta[0], _mm512_mask_mov_ps(beta_c,k2,fmod_zmm16r4(tmp1,v16_pi)));
			 gamma_c = gam;
			 const __m512 tmp2 = _mm512_fmadd_ps(v16_2,v16_pi,gamma_c);
			 k3 = _mm512_cmp_ps_mask(gamma_c,v16_0,k3,_CMP_LT_OQ);
			 _mm512_store_ps(&gamma[0], _mm512_mask_mov_ps(gamma_c,k3,fmod_zmm16r4(tmp2,v16_2pi)));
		    }


		  
                     
                      __ATTR_ALWAYS_INLINE__
		      static inline
	              void
		      q4x16_to_ea3x16_u_zmm16r4(const __m512 q_x,
		                                const __m512 q_y,
					        const __m512 q_z,
					        const __m512 q_w,
					        float * __restrict  alpha,
					        float * __restrict  beta,
					        float * __restrict  gamma) {

                         const __m512 qxw = _mm512_fmadd_ps(q_x,q_x,_mm512_mul_ps(q_w,q_w));
			 const __m512 qyz = _mm512_fmadd_ps(q_y,q_y,_mm512_mul_ps(q_z,q_z));
			 const __m512 chi = _mm512_sqrt_ps(_mm512_mul_ps(qxw,qyz));
			 __m512 alpha_c;
			 __m512 beta_c;
			 __m512 gamma_c;
			 __m512 alp,bet,gam;
			 __mmask16 qyz0 = 0x0;
			 __mmask16 qxw0 = 0x0;
			 __mmask16 k1   = 0x0;
			 __mmask16 k2   = 0x0;
			 __mmask16 k3   = 0x0;
			 qyz0 = _mm512_cmp_ps_mask(qyz,v16_0,_CMP_EQ_OQ);
			 qxw0 = _mm512_cmp_ps_mask(qxw,v16_0,_CMP_EQ_OQ);
			 if(1==qyz0) {
			    _mm512_storeu_ps(&alpha[0], _mm512_atan2_ps(_mm512_mul_ps(v16_n2,_mm512_mul_ps(q_x,q_w)),
			                              _mm512_fmsub_ps(q_x,q_x,_mm512_mul_ps(q_w,q_w))));
			    _mm512_storeu_ps(&beta[0],   v16_0);
			    _mm512_storeu_ps(&gamma[0],  v16_0);
			   // return;
			 } else if(1==qxw0) {
                            _mm512_storeu_ps(&alpha[0], _mm512_atan2_ps(_mm512_mul_ps(v16_2,_mm512_mul_ps(q_y,q_z)),
			                              _mm512_fmsub_ps(q_y,q_y,_mm512_mul_ps(q_z,q_z))));
			    _mm512_storeu_ps(&beta[0], v16_pi);
			    _mm512_storeu_ps(&gama[0], v16_0);
			    //return;
			 }
			 else {
			    const __m512 t0 = _mm512_mul_ps(q_y,q_w);
			    const __m512 c0 = _mm512_fmadd_ps(q_x,q_z,t0);
			    const __m512 t1 = _mm512_mul_ps(q_z,q_w);
			    const __m512 c1 = _mm512_fmsub_ps(q_x,q_y,t1);
			 
			    alp = _mm512_atan2_ps(_mm512_mul_ps(chi,_mm512_mul_ps(c0,v16_n1)),
			                            _mm512_mul_ps(chi,_mm512_mul_ps(c1,v16_n1)));
			    bet  = _mm512_atan2_ps(_mm512_mul_ps(v16_2,chi),_mm512_sub_ps(qxw,qyz));
			    const __m512 c2 = _mm512_fmadd_ps(q_x,_qy,t1);
			    gam = _mm512_atan2_ps(_mm512_mul_ps(chi,_mm512_mul_ps(c0,v16_1)),
			                            _mm512_mul_ps(chi,_mm512_mul_ps(c2,v16_n1)));
			 }
			 alpha_c = alp;
			 const __m512 tmp0 = _mm512_fmadd_ps(v16_2,v16_pi,alpha_c);
			 k1 = _mm512_cmp_ps_mask(alpha_c,v16_0,_CMP_LT_OQ);
			 _mm512_storeu_ps(&alpha[0],_mm512_mask_mov_ps(alpha_c,k1,fmod_zmm16r4(tmp0,v16_2pi)));
			 beta_c = bet;
			 const _mm512 tmp1 = _mm512_fmadd_ps(v16_2,v16_pi,beta_c);
			 k2 = _mm512_cmp_ps_mask(beta_c,v16_0,k2,_CMP_LT_OQ);
			 _mm512_storeu_ps(&beta[0], _mm512_mask_mov_ps(beta_c,k2,fmod_zmm16r4(tmp1,v16_pi)));
			 gamma_c = gam;
			 const __m512 tmp2 = _mm512_fmadd_ps(v16_2,v16_pi,gamma_c);
			 k3 = _mm512_cmp_ps_mask(gamma_c,v16_0,k3,_CMP_LT_OQ);
			 _mm512_storeu_ps(&gamma[0], _mm512_mask_mov_ps(gamma_c,k3,fmod_zmm16r4(tmp2,v16_2pi)));
		    }




		     
                      __ATTR_ALWAYS_INLINE__
		      static inline
	              void
		      q4x8_to_ea3x8_a_zmm8r8(   const __m512d q_x,
		                                const __m512d q_y,
					        const __m512d q_z,
					        const __m512d q_w,
					        double * __restrict __ATTR_ALIGN__(64) alpha,
					        double * __restrict __ATTR_ALIGN__(64) beta,
					        double * __restrict __ATTR_ALIGN__(64) gamma) {

                         const __m512d qxw = _mm512_fmadd_pd(q_x,q_x,_mm512_mul_pd(q_w,q_w));
			 const __m512d qyz = _mm512_fmadd_pd(q_y,q_y,_mm512_mul_pd(q_z,q_z));
			 const __m512d chi = _mm512_sqrt_pd(_mm512_mul_pd(qxw,qyz));
			 __m512d alpha_c;
			 __m512d beta_c;
			 __m512d gamma_c;
			 __m512d alp,bet,gam;
			 __mmask8 qyz0 = 0x0;
			 __mmask8 qxw0 = 0x0;
			 __mmask8 k1   = 0x0;
			 __mmask8 k2   = 0x0;
			 __mmask8 k3   = 0x0;
			 qyz0 = _mm512_cmp_pd_mask(qyz,v8_0,_CMP_EQ_OQ);
			 qxw0 = _mm512_cmp_pd_mask(qxw,v8_0,_CMP_EQ_OQ);
			 if(__builtin_expect(1==qyz0,0)) {
			    _mm512_store_pd(&alpha[0],  _mm512_atan2_ps(_mm512_mul_pd(v8_n2,_mm512_mul_pd(q_x,q_w)),
			                              _mm512_fmsub_pd(q_x,q_x,_mm512_mul_pd(q_w,q_w))));
			    _mm512_store_pd(&beta[0], v8_0);
			    _mm512_store_pd(&gamma[0],v8_0);
			    //return;
			 } else if(__builtin_expect(1==qxw0,0)) {
                            _mm512_store_pd(&alpha[0], _mm512_atan2_pd(_mm512_mul_pd(v16_2,_mm512_mul_pd(q_y,q_z)),
			                              _mm512_fmsub_pd(q_y,q_y,_mm512_mul_pd(q_z,q_z))));
			    _mm512_store_pd(&beta[0], v8_pi);
			    _mm512_store_pd(&gamma[0],v8_0);
			    //return;
			 }
			 else {
			    const __m512d t0 = _mm512_mul_pd(q_y,q_w);
			    const __m512d c0 = _mm512_fmadd_pd(q_x,q_z,t0);
			    const __m512d t1 = _mm512_mul_pd(q_z,q_w);
			    const __m512d c1 = _mm512_fmsub_pd(q_x,q_y,t1);
			 
			    alp = _mm512_atan2_pd(_mm512_mul_pd(chi,_mm512_mul_pd(c0,v8_n1)),
			                            _mm512_mul_pd(chi,_mm512_mul_ps(c1,v8_n1)));
			    bet  = _mm512_atan2_pd(_mm512_mul_pd(v8_2,chi),_mm512_sub_pd(qxw,qyz));
			    const __m512d c2 = _mm512_fmadd_pd(q_x,_qy,t1);
			    gam = _mm512_atan2_pd(_mm512_mul_pd(chi,_mm512_mul_pd(c0,v8_1)),
			                            _mm512_mul_pd(chi,_mm512_mul_pd(c2,v8_n1)));
			 }
			 alpha_c = alp;
			 const __m512d tmp0 = _mm512_fmadd_pd(v8_2,v16_pi,alpha_c);
			 k1 = _mm512_cmp_pd_mask(alpha_c,v8_0,_CMP_LT_OQ);
			 _mm512_store_pd(&alpha[0], _mm512_mask_mov_ps(alpha_c,k1,fmod_zmm8r8(tmp0,v8_2pi)));
			 beta_c = bet;
			 const _mm512d tmp1 = _mm512_fmadd_ps(v8_2,v8_pi,beta_c);
			 k2 = _mm512_cmp_pd_mask(beta_c,v8_0,k2,_CMP_LT_OQ);
			 _mm512_store_pd(&beta[0], _mm512_mask_mov_pd(beta_c,k2,fmod_zmm8r8(tmp1,v8_pi)));
			 gamma_c = gam;
			 const __m512d tmp2 = _mm512_fmadd_pd(v8_2,v8_pi,gamma_c);
			 k3 = _mm512_cmp_pd_mask(gamma_c,v8_0,k3,_CMP_LT_OQ);
			 _mm512_store_pd(&gamma[0], _mm512_mask_mov_pd(gamma_c,k3,fmod_zmm8r8(tmp2,v8_2pi)));
		    }


		     
                      __ATTR_ALWAYS_INLINE__
		      static inline
	              void
		      q4x8_to_ea3x8_u_zmm8r8(   const __m512d q_x,
		                                const __m512d q_y,
					        const __m512d q_z,
					        const __m512d q_w,
					        double * __restrict  alpha,
					        double * __restrict  beta,
					        double * __restrict  gamma) {

                         const __m512d qxw = _mm512_fmadd_pd(q_x,q_x,_mm512_mul_pd(q_w,q_w));
			 const __m512d qyz = _mm512_fmadd_pd(q_y,q_y,_mm512_mul_pd(q_z,q_z));
			 const __m512d chi = _mm512_sqrt_pd(_mm512_mul_pd(qxw,qyz));
			 __m512d alpha_c;
			 __m512d beta_c;
			 __m512d gamma_c;
			 __m512d alp,bet,gam;
			 __mmask8 qyz0 = 0x0;
			 __mmask8 qxw0 = 0x0;
			 __mmask8 k1   = 0x0;
			 __mmask8 k2   = 0x0;
			 __mmask8 k3   = 0x0;
			 qyz0 = _mm512_cmp_pd_mask(qyz,v8_0,_CMP_EQ_OQ);
			 qxw0 = _mm512_cmp_pd_mask(qxw,v8_0,_CMP_EQ_OQ);
			 if(__builtin_expect(1==qyz0,0)) {
			    _mm512_storeu_ps(&alpha[0],  _mm512_atan2_ps(_mm512_mul_pd(v8_n2,_mm512_mul_pd(q_x,q_w)),
			                              _mm512_fmsub_pd(q_x,q_x,_mm512_mul_pd(q_w,q_w))));
			    _mm512_storeu_pd(&beta[0], v8_0);
			    _mm512_storeu_pd(&gamma[0],v8_0);
			    //return;
			 } else if(__builtin_expect(1==qxw0,0)) {
                            _mm512_storeu_pd(&alpha[0], _mm512_atan2_pd(_mm512_mul_pd(v16_2,_mm512_mul_pd(q_y,q_z)),
			                              _mm512_fmsub_pd(q_y,q_y,_mm512_mul_pd(q_z,q_z))));
			    _mm512_storeu_pd(&beta[0], v8_pi);
			    _mm512_storeu_pd(&gamma[0],v8_0);
			    //return;
			 }
			 else {
			    const __m512d t0 = _mm512_mul_pd(q_y,q_w);
			    const __m512d c0 = _mm512_fmadd_pd(q_x,q_z,t0);
			    const __m512d t1 = _mm512_mul_pd(q_z,q_w);
			    const __m512d c1 = _mm512_fmsub_pd(q_x,q_y,t1);
			 
			    alp = _mm512_atan2_pd(_mm512_mul_pd(chi,_mm512_mul_pd(c0,v8_n1)),
			                            _mm512_mul_pd(chi,_mm512_mul_ps(c1,v8_n1)));
			    bet  = _mm512_atan2_pd(_mm512_mul_pd(v8_2,chi),_mm512_sub_pd(qxw,qyz));
			    const __m512d c2 = _mm512_fmadd_pd(q_x,q_y,t1);
			    gam = _mm512_atan2_pd(_mm512_mul_pd(chi,_mm512_mul_pd(c0,v8_1)),
			                            _mm512_mul_pd(chi,_mm512_mul_pd(c2,v8_n1)));
			 }
			 alpha_c = alp;
			 const __m512d tmp0 = _mm512_fmadd_pd(v8_2,v16_pi,alpha_c);
			 k1 = _mm512_cmp_pd_mask(alpha_c,v8_0,_CMP_LT_OQ);
			 _mm512_storeu_pd(&alpha[0], _mm512_mask_mov_ps(alpha_c,k1,fmod_zmm8r8(tmp0,v8_2pi)));
			 beta_c = bet;
			 const _mm512d tmp1 = _mm512_fmadd_ps(v8_2,v8_pi,beta_c);
			 k2 = _mm512_cmp_pd_mask(beta_c,v8_0,k2,_CMP_LT_OQ);
			 _mm512_storeu_pd(&beta[0], _mm512_mask_mov_pd(beta_c,k2,fmod_zmm8r8(tmp1,v8_pi)));
			 gamma_c = gam;
			 const __m512d tmp2 = _mm512_fmadd_pd(v8_2,v8_pi,gamma_c);
			 k3 = _mm512_cmp_pd_mask(gamma_c,v8_0,k3,_CMP_LT_OQ);
			 _mm512_storeu_pd(&gamma[0], _mm512_mask_mov_pd(gamma_c,k3,fmod_zmm8r8(tmp2,v8_2pi)));
		    }




		   /*
                       Convert unit quaternion to axis angle pair
                    */


		    
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void
		      q4x16_to_ax4x16_a_zmm16r4(const __m512 q_x,
		                                const __m512 q_y,
					        const __m512 q_z,
					        const __m512 q_w,
					        float * __restrict __ATTR_ALIGN__(64) ax_1,
					        float * __restrict __ATTR_ALIGN__(64) ax_2,
					        float * __restrict __ATTR_ALIGN__(64) ax_3,
					        float * __restrict __ATTR_ALIGN__(64) ax_4) {

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
                             _mm512_store_ps(&ax_1[0], v16_0);
			     _mm512_store_ps(&ax_2[0], v16_0);
			     _mm512_store_ps(&ax_3[0], v16_1);
			     _mm512_store_ps(&ax_4[0], v16_0);
			     return;
			  }
			  k2 = _mm512_cmp_ps_mask(q_x,_v16_0,_CMP_NEQ_OQ);
			  if(1==k2) {
                             const register __m512 s = _mm512_div_ps(
			                          zmm16r4_sign_zmm16r4(v16_1,q_x),
						  norm2_zmm16r4(q_y,q_z,q_w));
			     _mm512_store_ps(&ax_1[0], _mm512_mul_ps(q_y,s));
			     _mm512_store_ps(&ax_2[0], _mm512_mul_ps(q_z,s));
			     _mm512_store_ps(&ax_3[0],  _mm512_mul_ps(q_w,s));
			     const register __m512 omega = _mm512_mul_ps(v16_2,
			                              _mm512_acos_ps(clip_zmm16r4(q_x,v16_n1,v16_1)));
			     _mm512_store_ps(&ax_4[0] ,omega);
			     return;
			  }
			  else {
                             _mm512_store_ps(&ax_1[0], q_y);
			     _mm512_store_ps(&ax_2[0], q_z);
			     _mm512_store_ps(&ax_3[0], q_w);
			     _mm512_store_ps(&ax_4[0], v16_pi);
			     return;
			  }
		     }


		      
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void
		      q4x16_to_ax4x16_u_zmm16r4(const __m512 q_x,
		                                const __m512 q_y,
					        const __m512 q_z,
					        const __m512 q_w,
					        float * __restrict  ax_1,
					        float * __restrict  ax_2,
					        float * __restrict  ax_3,
					        float * __restrict  ax_4) {

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
                             _mm512_storeu_ps(&ax_1[0], v16_0);
			     _mm512_storeu_ps(&ax_2[0], v16_0);
			     _mm512_storeu_ps(&ax_3[0], v16_1);
			     _mm512_storeu_ps(&ax_4[0], v16_0);
			     return;
			  }
			  k2 = _mm512_cmp_ps_mask(q_x,_v16_0,_CMP_NEQ_OQ);
			  if(1==k2) {
                             const register __m512 s = _mm512_div_ps(
			                          zmm16r4_sign_zmm16r4(v16_1,q_x),
						  norm2_zmm16r4(q_y,q_z,q_w));
			     _mm512_storeu_ps(&ax_1[0], _mm512_mul_ps(q_y,s));
			     _mm512_storeu_ps(&ax_2[0], _mm512_mul_ps(q_z,s));
			     _mm512_storeu_ps(&ax_3[0],  _mm512_mul_ps(q_w,s));
			     const register __m512 omega = _mm512_mul_ps(v16_2,
			                              _mm512_acos_ps(clip_zmm16r4(q_x,v16_n1,v16_1)));
			     _mm512_store_ps(&ax_4[0] ,omega);
			     return;
			  }
			  else {
                             _mm512_storeu_ps(&ax_1[0], q_y);
			     _mm512_storeu_ps(&ax_2[0], q_z);
			     _mm512_storeu_ps(&ax_3[0], q_w);
			     _mm512_storeu_ps(&ax_4[0], v16_pi);
			     return;
			  }
		     }





		     
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void
		      q4x8_to_ax4x8_a_zmm8r8(   const __m512d q_x,
		                              const __m512d q_y,
					      const __m512d q_z,
					      const __m512d q_w,
					      double * __restrict __ATTR_ALIGN__(64) ax_1,
					      double * __restrict __ATTR_ALIGN__(64) ax_2,
					      double * __restrict __ATTR_ALIGN__(64) ax_3,
					      double * __restrict __ATTR_ALIGN__(64) ax_4) {

                          const register __m512d t0 = _mm512_mul_pd(q_y,q_y);
			  const register __m512d t1 = _mm512_mul_pd(q_z,q_z);
			  const register __m512d t2 = _mm512_mul_pd(q_w,q_w);
			  const register __m512d v0 = _mm512_add_pd(t0,_mm512_add_pd(t1,t2));
			  __mmask8 k1 = 0x0;
			  __mmask8 k2 = 0x0;
			  k1 = _mm512_cmp_pd_mask(v0,v8_0,_CMP_EQ_OQ);
			  if(__builtin_expect(1==k1,0)) {
                             _mm512_store_pd(&ax_1[0], v8_0);
			     _mm512_store_pd(&ax_2[0], v8_0);
			     _mm512_store_pd(&ax_3[0], v8_1);
			     _mm512_store_pd(&ax_4[0], v8_0);
			     return;
			  }
			  k2 = _mm512_cmp_pd_mask(q_x,_v8_0,_CMP_NEQ_OQ);
			  if(1==k2) {
                             const register __m512 s = _mm512_div_pd(
			                          zmm8r8_sign_zmm8r8(v8_1,q_x),
						  norm2_zmm8r8(q_y,q_z,q_w));
			     _mm512_store_pd(&ax_1[0], _mm512_mul_pd(q_y,s));
			     _mm512_store_pd(&ax_2[0], _mm512_mul_pd(q_z,s));
			     _mm512_store_pd(&ax_3[0], _mm512_mul_pd(q_w,s));
			     const register __m512 omega = _mm512_mul_pd(v16_2,
			                              _mm512_acos_pd(clip_zmm8r8(q_x,v8_n1,v8_1)));
			     _mm512_store_pd(&ax_4[0], omega);
			     return;
			  }
			  else {
                             _mm512_store_pd(&ax_1[0], q_y);
			     _mm512_store_pd(&ax_2[0], q_z);
			     _mm512_store_pd(&ax_3[0], q_w);
			     _mm512_store_pd(&ax_4[0], v8_pi);
			     return;
			  }
		     }


		      
                      __ATTR_ALWAYS_INLINE__
		     static inline
		      void
		      q4x8_to_ax4x8_u_zmm8r8( const __m512d q_x,
		                              const __m512d q_y,
					      const __m512d q_z,
					      const __m512d q_w,
					      double * __restrict  ax_1,
					      double * __restrict  ax_2,
					      double * __restrict  ax_3,
					      double * __restrict  ax_4) {

                          const register __m512d t0 = _mm512_mul_pd(q_y,q_y);
			  const register __m512d t1 = _mm512_mul_pd(q_z,q_z);
			  const register __m512d t2 = _mm512_mul_pd(q_w,q_w);
			  const register __m512d v0 = _mm512_add_pd(t0,_mm512_add_pd(t1,t2));
			  __mmask8 k1 = 0x0;
			  __mmask8 k2 = 0x0;
			  k1 = _mm512_cmp_pd_mask(v0,v8_0,_CMP_EQ_OQ);
			  if(__builtin_expect(1==k1,0)) {
                             _mm512_storeu_pd(&ax_1[0], v8_0);
			     _mm512_storeu_pd(&ax_2[0], v8_0);
			     _mm512_storeu_pd(&ax_3[0], v8_1);
			     _mm512_storeu_pd(&ax_4[0], v8_0);
			     return;
			  }
			  k2 = _mm512_cmp_pd_mask(q_x,_v8_0,_CMP_NEQ_OQ);
			  if(1==k2) {
                             const register __m512 s = _mm512_div_pd(
			                          zmm8r8_sign_zmm8r8(v8_1,q_x),
						  norm2_zmm8r8(q_y,q_z,q_w));
			     _mm512_storeu_pd(&ax_1[0], _mm512_mul_pd(q_y,s));
			     _mm512_storeu_pd(&ax_2[0], _mm512_mul_pd(q_z,s));
			     _mm512_storeu_pd(&ax_3[0], _mm512_mul_pd(q_w,s));
			     const register __m512 omega = _mm512_mul_pd(v16_2,
			                              _mm512_acos_pd(clip_zmm8r8(q_x,v8_n1,v8_1)));
			     _mm512_storeu_pd(&ax_4[0], omega);
			     return;
			  }
			  else {
                             _mm512_storeu_pd(&ax_1[0], q_y);
			     _mm512_storeu_pd(&ax_2[0], q_z);
			     _mm512_storeu_pd(&ax_3[0], q_w);
			     _mm512_storeu_pd(&ax_4[0], v8_pi);
			     return;
			  }
		     }



		    /*
                        Convert unit quaternion to Rodrigues vector
                     */
#include <limits>

		      
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void
		      q4x16_to_rv4x16_a_zmm16r4(const __m512 q_x,
		                              const __m512 q_y,
					      const __m512 q_z,
					      const __m512 q_w,
					      float * __restrict __ATTR_ALIGN__(64) r_x,
					      float * __restrict __ATTR_ALIGN__(64) r_y,
					      float * __restrict __ATTR_ALIGN__(64) r_z,
					      float * __restrict __ATTR_ALIGN__(64) r_w) {

                          const register __m512 thr = _mm512_set1_ps(1.0e-8);
			  const register __m512 inf = _mm512_set1_ps(std::numeric_limits<float>::infinity());
			  register __m512 t0 = v16_0;
			  register __m512 s  = v16_0;
			  __mmask16 k1 = 0x0;
			  __mmask16 k2 = 0x0;
			  k1 = _mm512_cmp_ps_mask(_mm512_abs_ps(q_x),thr,_CMP_LT_OQ);
			  if(__builtin_expect(1==k1,0)) {
                             _mm512_store_ps(&r_x[0], q_y);
			     _mm512_store_ps(&r_y[0], q_z);
			     _mm512_store_ps(&r_z[0], q_w);
			     _mm512_store_ps(&r_w[0], inf);
			  }
			  else {
                               s   = norm2_zmm16r4(q_y,q_z,q_w);
			       k2  = _mm512_cmp_ps_mask(s,thr,_CMP_LT_OQ);
			       _mm512_store_ps(&r_x[0],_mm512_mask_blend_ps(k2,_mm512_div_ps(q_y,s),v16_0));
			       t0  = _mm512_acos_ps(clip_zmm16r4(q_x,v16_n1,v16_1));
			       _mm512_store_ps(&r_y[0], _mm512_mask_blend_ps(k2,_mm512_div_ps(q_z,s),v16_0));
			       _mm512_store_ps(&r_z[0], _mm512_mask_blend_ps(k2,_mm512_div_ps(q_w,s),v16_n1));
			       _mm512_store_ps(&r_w[0], _mm512_mask_blend_ps(k2,_mm512_tan_ps(t0),v16_0));
			  }
		     }


		      
                      __ATTR_ALWAYS_INLINE__
		     static inline
		      void
		      q4x16_to_rv4x16_u_zmm16r4(const __m512 q_x,
		                              const __m512 q_y,
					      const __m512 q_z,
					      const __m512 q_w,
					      float * __restrict r_x,
					      float * __restrict r_y,
					      float * __restrict r_z,
					      float * __restrict r_w) {

                          const register __m512 thr = _mm512_set1_ps(1.0e-8);
			  const register __m512 inf = _mm512_set1_ps(std::numeric_limits<float>::infinity());
			  register __m512 t0 = v16_0;
			  register __m512 s  = v16_0;
			  __mmask16 k1 = 0x0;
			  __mmask16 k2 = 0x0;
			  k1 = _mm512_cmp_ps_mask(_mm512_abs_ps(q_x),thr,_CMP_LT_OQ);
			  if(__builtin_expect(1==k1,0)) {
                             _mm512_storeu_ps(&r_x[0], q_y);
			     _mm512_storeu_ps(&r_y[0], q_z);
			     _mm512_storeu_ps(&r_z[0], q_w);
			     _mm512_storeu_ps(&r_w[0], inf);
			  }
			  else {
                               s   = norm2_zmm16r4(q_y,q_z,q_w);
			       k2  = _mm512_cmp_ps_mask(s,thr,_CMP_LT_OQ);
			       _mm512_storeu_ps(&r_x[0],_mm512_mask_blend_ps(k2,_mm512_div_ps(q_y,s),v16_0));
			       t0  = _mm512_acos_ps(clip_zmm16r4(q_x,v16_n1,v16_1));
			       _mm512_storeu_ps(&r_y[0], _mm512_mask_blend_ps(k2,_mm512_div_ps(q_z,s),v16_0));
			       _mm512_storeu_ps(&r_z[0], _mm512_mask_blend_ps(k2,_mm512_div_ps(q_w,s),v16_n1));
			       _mm512_storeu_ps(&r_w[0], _mm512_mask_blend_ps(k2,_mm512_tan_ps(t0),v16_0));
			  }
		     }


		     
	              
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void
		      q4x8_to_rv4x8_a_zmm8r8(   const __m512d q_x,
		                              const __m512d q_y,
					      const __m512d q_z,
					      const __m512d q_w,
					      double * __restrict __ATTR_ALIGN__(64) r_x,
					      double * __restrict __ATTR_ALIGN__(64) r_y,
					      double * __restrict __ATTR_ALIGN__(64) r_z,
					      double * __restrict __ATTR_ALIGN__(64) r_w) {

                          const register __m512d thr = _mm512_set1_pd(1.0e-8);
			  const register __m512d inf = _mm512_set1_pd(std::numeric_limits<double>::infinity());
			  register __m512d t0 = v8_0;
			  register __m512d s  = v8_0;
			  __mmask8 k1 = 0x0;
			  __mmask8 k2 = 0x0;
			  k1 = _mm512_cmp_pd_mask(_mm512_abs_pd(q_x),thr,_CMP_LT_OQ);
			  if(__builtin_expect(1==k1,0)) {
                             _mm512_store_pd(&r_x[0], q_y);
			     _mm512_store_pd(&r_y[0], q_z);
			     _mm512_store_pd(&r_z[0], q_w);
			     _mm512_store_pd(&r_w[0], inf);
			  }
			  else {
                               s   = norm2_zmm8r8(q_y,q_z,q_w);
			       k2  = _mm512_cmp_pd_mask(s,thr,_CMP_LT_OQ);
			       _mm512_store_pd(&r_x[0], _mm512_mask_blend_pd(k2,_mm512_div_pd(q_y,s),v8_0));
			       t0  = _mm512_acos_pd(clip_zmm8r8(q_x,v8_n1,v8_1));
			       _mm512_store_pd(&r_y[0], _mm512_mask_blend_pd(k2,_mm512_div_pd(q_z,s),v8_0));
			       _mm512_store_pd(&r_z[0], _mm512_mask_blend_pd(k2,_mm512_div_pd(q_w,s),v8_n1));
			       _mm512_store_pd(&r_w[0], _mm512_mask_blend_pd(k2,_mm512_tan_pd(t0),v8_0));
			  }
		     }

		     
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void
		      q4x8_to_rv4x8_u_zmm8r8(   const __m512d q_x,
		                              const __m512d q_y,
					      const __m512d q_z,
					      const __m512d q_w,
					      double * __restrict  r_x,
					      double * __restrict  r_y,
					      double * __restrict  r_z,
					      double * __restrict  r_w) {

                          const register __m512d thr = _mm512_set1_pd(1.0e-8);
			  const register __m512d inf = _mm512_set1_pd(std::numeric_limits<double>::infinity());
			  register __m512d t0 = v8_0;
			  register __m512d s  = v8_0;
			  __mmask8 k1 = 0x0;
			  __mmask8 k2 = 0x0;
			  k1 = _mm512_cmp_pd_mask(_mm512_abs_pd(q_x),thr,_CMP_LT_OQ);
			  if(__builtin_expect(1==k1,0)) {
                             _mm512_storeu_pd(&r_x[0], q_y);
			     _mm512_storeu_pd(&r_y[0], q_z);
			     _mm512_storeu_pd(&r_z[0], q_w);
			     _mm512_storeu_pd(&r_w[0], inf);
			  }
			  else {
                               s   = norm2_zmm8r8(q_y,q_z,q_w);
			       k2  = _mm512_cmp_pd_mask(s,thr,_CMP_LT_OQ);
			       _mm512_storeu_pd(&r_x[0], _mm512_mask_blend_pd(k2,_mm512_div_pd(q_y,s),v8_0));
			       t0  = _mm512_acos_pd(clip_zmm8r8(q_x,v8_n1,v8_1));
			       _mm512_storeu_pd(&r_y[0], _mm512_mask_blend_pd(k2,_mm512_div_pd(q_z,s),v8_0));
			       _mm512_storeu_pd(&r_z[0], _mm512_mask_blend_pd(k2,_mm512_div_pd(q_w,s),v8_n1));
			       _mm512_storeu_pd(&r_w[0], _mm512_mask_blend_pd(k2,_mm512_tan_pd(t0),v8_0));
			  }
		     }


		     /*
                           Orientation i.e. (Direct Cosine Matrix)  matrix to Euler angles.
                       */


		     
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void
		      rmat9x16_to_ea3x16_a_zmm16r4(const DCM9x16 rm,
		                                 float * __restrict __ATTR_ALIGN__(64) alpha,
						 float * __restrict __ATTR_ALIGN__(64) beta,
						 float * __restrict __ATTR_ALIGN__(64) gamma) {

                           const    __m512  thr  = _mm512_set1_ps(1.0e-8);
			   register __m512  t0   = v16_0;
			   register __m512  t1   = v16_0;
			   register __m512  t2   = v16_0;
			   register __m512  zeta = v16_0;
			   register __m512  al_c = v16_0;
			   register __m512  be_c = v16_0;
			   register __m512  ga_c = v16_0;
			   register __m512  alp  = v16_0;
			   register __m512  bet  = v16_0;
			   register __m512  gam  = v16_0;
			   __mmask16 k1 = 0x0;
			   __mmask16 k2 = 0x0;
			   __mmask16 k3 = 0x0;
			   __mmask16 k4 = 0x0;
			   k1 = _mm512_cmp_ps_mask(rm.m_vRow9,v16_1,_CMP_NEQ_OQ);
			   t0 = _mm512_mul_ps(rm.m_vRow9,rm.m_vRow9);
			   zeta = _mm512_div_ps(_mm512_sqrt_ps(_mm512_sub_ps(v16_1,t0)));
			   t0 = _mm512_sub_ps(v16_0,rm.vRow8);
			   alp = _mm512_mask_blend_ps(k1,_mm512_atan2_ps(rm.vRow2,rm.vRow1),
			                                   _mm512_atan2_ps(_mm512_mul_ps(rm.vRow7,zeta),
			   				                   _mm512_mul_ps(t0,zeta)));
			   t1    = _mm512_mul_ps(v16_1o2,_mm512_mul_ps(v16_pi,_mm512_sub_ps(v16_1,rm.vRow9)));
			   bet  = _mm512_mask_blend_ps(k1,t1,_mm512_acos_ps(rm.vRow9));
			   gam = _mm512_mask_blend_ps(k1,v16_0,_mm512_atan_ps(_mm512_mul_ps(rm.vRow3,zeta),
			                                                        _mm512_mul_ps(rm.vRow6,zeta)));
			   al_c = alp;
			   be_c = bet;
			   ga_c = gam;
			   k2 = _mm512_cmp_ps_mask(_mm512_abs_ps(alp),thr,_CMP_LT_OQ);
			   alp = _mm512_mask_mov_ps(al_c,k2,v16_0);
			   k3 = _mm512_cmp_ps_mask(_mm512_abs_ps(bet),thr,_CMP_LT_OQ);
			   bet = _mm512_mask_mov_ps(be_c,k3,v16_0);
			   k4 = _mm512_cmp_ps_mask(_mm512_abs_ps(gam),thr,_CMP_LT_OQ);
			   gam = _mm512_mask_mov_ps(ga_c,k4,v16_0);
			   //al_c = alpha;
			   t0 = _mm512_add_ps(alp,v16_2pi);
			   k2 = _mm512_cmp_ps_mask(alp,v16_0,_CMP_LT_OQ);
			   alp = _mm512_mask_mov_ps(al_c,k2,fmod_zmm16r4(t0,v16_2pi));
			   //be_c = beta;
			   t1 = _mm512_add_ps(bet,v16_2pi);
			   k3 = _mm512_cmp_ps_mask(bet,v16_0,_CMP_LT_OQ);
			   bet = _mm512_mask_mov_ps(be_c,k3,fmod_zmm16r4(t1,v16_pi));
			   //ga_c = gamma;
			   t2 = _mm512_add_ps(gam,v16_2pi);
			   k4 = _mm512_cmp_ps_mask(gam,v16_0,_CMP_LT_OQ);
			   gam = _mm512_mask_mov_ps(ga_c,k4,fmod_zmm16r4(t2,v16_2pi));
			   _mm512_store_ps(&alpha[0],alp);
			   _mm512_store_ps(&beta[0], bet);
			   _mm512_store_ps(&gamma[0],gam);
		     }


		     
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void
		      rmat9x16_to_ea3x16_u_zmm16r4(const DCM9x16 rm,
		                                 float * __restrict  alpha,
						 float * __restrict  beta,
						 float * __restrict  gamma) {

                           const    __m512  thr  = _mm512_set1_ps(1.0e-8);
			   register __m512  t0   = v16_0;
			   register __m512  t1   = v16_0;
			   register __m512  t2   = v16_0;
			   register __m512  zeta = v16_0;
			   register __m512  al_c = v16_0;
			   register __m512  be_c = v16_0;
			   register __m512  ga_c = v16_0;
			   register __m512  alp  = v16_0;
			   register __m512  bet  = v16_0;
			   register __m512  gam  = v16_0;
			   __mmask16 k1 = 0x0;
			   __mmask16 k2 = 0x0;
			   __mmask16 k3 = 0x0;
			   __mmask16 k4 = 0x0;
			   k1 = _mm512_cmp_ps_mask(rm.m_vRow9,v16_1,_CMP_NEQ_OQ);
			   t0 = _mm512_mul_ps(rm.m_vRow9,rm.m_vRow9);
			   zeta = _mm512_div_ps(_mm512_sqrt_ps(_mm512_sub_ps(v16_1,t0)));
			   t0 = _mm512_sub_ps(v16_0,rm.vRow8);
			   alp = _mm512_mask_blend_ps(k1,_mm512_atan2_ps(rm.vRow2,rm.vRow1),
			                                   _mm512_atan2_ps(_mm512_mul_ps(rm.vRow7,zeta),
			   				                   _mm512_mul_ps(t0,zeta)));
			   t1    = _mm512_mul_ps(v16_1o2,_mm512_mul_ps(v16_pi,_mm512_sub_ps(v16_1,rm.vRow9)));
			   bet  = _mm512_mask_blend_ps(k1,t1,_mm512_acos_ps(rm.vRow9));
			   gam = _mm512_mask_blend_ps(k1,v16_0,_mm512_atan_ps(_mm512_mul_ps(rm.vRow3,zeta),
			                                                        _mm512_mul_ps(rm.vRow6,zeta)));
			   al_c = alp;
			   be_c = bet;
			   ga_c = gam;
			   k2 = _mm512_cmp_ps_mask(_mm512_abs_ps(alp),thr,_CMP_LT_OQ);
			   alp = _mm512_mask_mov_ps(al_c,k2,v16_0);
			   k3 = _mm512_cmp_ps_mask(_mm512_abs_ps(bet),thr,_CMP_LT_OQ);
			   bet = _mm512_mask_mov_ps(be_c,k3,v16_0);
			   k4 = _mm512_cmp_ps_mask(_mm512_abs_ps(gam),thr,_CMP_LT_OQ);
			   gam = _mm512_mask_mov_ps(ga_c,k4,v16_0);
			   //al_c = alpha;
			   t0 = _mm512_add_ps(alp,v16_2pi);
			   k2 = _mm512_cmp_ps_mask(alp,v16_0,_CMP_LT_OQ);
			   alp = _mm512_mask_mov_ps(al_c,k2,fmod_zmm16r4(t0,v16_2pi));
			   //be_c = beta;
			   t1 = _mm512_add_ps(bet,v16_2pi);
			   k3 = _mm512_cmp_ps_mask(bet,v16_0,_CMP_LT_OQ);
			   bet = _mm512_mask_mov_ps(be_c,k3,fmod_zmm16r4(t1,v16_pi));
			   //ga_c = gamma;
			   t2 = _mm512_add_ps(gam,v16_2pi);
			   k4 = _mm512_cmp_ps_mask(gam,v16_0,_CMP_LT_OQ);
			   gam = _mm512_mask_mov_ps(ga_c,k4,fmod_zmm16r4(t2,v16_2pi));
			   _mm512_storeu_ps(&alpha[0],alp);
			   _mm512_storeu_ps(&beta[0], bet);
			   _mm512_storeu_ps(&gamma[0],gam);
		     }



		     
                      __ATTR_ALWAYS_INLINE__
		     static inline
		      void
		      rmat9x8_to_ea3x8_a_zmm8r8(const DCM9x8 rm,
		                              double * __restrict __ATTR_ALIGN__(64) alpha,
					      double * __restrict __ATTR_ALIGN__(64) beta,
					      double * __restrict __ATTR_ALIGN__(64) gamma) {

                           const    __m512d  thr  = _mm512_set1_pd(1.0e-8);
			   register __m512d  t0   = v8_0;
			   register __m512d  t1   = v8_0;
			   register __m512d  t2   = v8_0;
			   register __m512d  zeta = v8_0;
			   register __m512d  al_c = v8_0;
			   register __m512d  be_c = v8_0;
			   register __m512d  ga_c = v8_0;
			   register __m512d  alp  = v8_0;
			   register __m512d  bet  = v8_0;
			   register __m512d  gam  = v8_0;
			   __mmask8 k1 = 0x0;
			   __mmask8 k2 = 0x0;
			   __mmask8 k3 = 0x0;
			   __mmask8 k4 = 0x0;
			   k1 = _mm512_cmp_pd_mask(rm.m_vRow9,v8_1,_CMP_NEQ_OQ);
			   t0 = _mm512_mul_pd(rm.m_vRow9,rm.m_vRow9);
			   zeta = _mm512_div_pd(_mm512_sqrt_pd(_mm512_sub_pd(v8_1,t0)));
			   t0 = _mm512_sub_pd(v8_0,rm.vRow8);
			   alp = _mm512_mask_blend_pd(k1,_mm512_atan2_pd(rm.vRow2,rm.vRow1),
			                                   _mm512_atan2_pd(_mm512_mul_pd(rm.vRow7,zeta),
			   				                   _mm512_mul_pd(t0,zeta)));
			   t1    = _mm512_mul_pd(v8_1o2,_mm512_mul_pd(v8_pi,_mm512_sub_pd(v8_1,rm.vRow9)));
			   bet  = _mm512_mask_blend_pd(k1,t1,_mm512_acos_pd(rm.vRow9));
			   gam = _mm512_mask_blend_pd(k1,v8_0,_mm512_atan_pd(_mm512_mul_pd(rm.vRow3,zeta),
			                                                        _mm512_mul_pd(rm.vRow6,zeta)));
			   al_c = alp;
			   be_c = bet;
			   ga_c = gam;
			   k2 = _mm512_cmp_pd_mask(_mm512_abs_pd(alp),thr,_CMP_LT_OQ);
			   alp = _mm512_mask_mov_pd(al_c,k2,v8_0);
			   k3 = _mm512_cmp_pd_mask(_mm512_abs_pd(bet),thr,_CMP_LT_OQ);
			   bet = _mm512_mask_mov_pd(be_c,k3,v8_0);
			   k4 = _mm512_cmp_pd_mask(_mm512_abs_pd(gam),thr,_CMP_LT_OQ);
			   gam = _mm512_mask_mov_pd(ga_c,k4,v8_0);
			   //al_c = alpha;
			   t0 = _mm512_add_pd(alp,v8_2pi);
			   k2 = _mm512_cmp_pd_mask(alp,v8_0,_CMP_LT_OQ);
			   alp = _mm512_mask_mov_pd(al_c,k2,fmod_zmm8r8(t0,v8_2pi));
			   //be_c = beta;
			   t1 = _mm512_add_pd(bet,v8_2pi);
			   k3 = _mm512_cmp_pd_mask(bet,v8_0,_CMP_LT_OQ);
			   bet = _mm512_mask_mov_pd(be_c,k3,fmod_zmm8r8(t1,v8_pi));
			   //ga_c = gamma;
			   t2 = _mm512_add_pd(gam,v8_2pi);
			   k4 = _mm512_cmp_pd_mask(gam,v8_0,_CMP_LT_OQ);
			   gam = _mm512_mask_mov_pd(ga_c,k4,fmod_zmm8r8(t2,v8_2pi));
			   _mm512_store_pd(&alpha[0],alp);
			   _mm512_store_pd(&beta[0], bet);
			   _mm512_store_pd(&gamma[0],gam);
		     }


		     
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void
		      rmat9x8_to_ea3x8_u_zmm8r8(const DCM9x8 rm,
		                              double * __restrict  alpha,
					      double * __restrict  beta,
					      double * __restrict  gamma) {

                           const    __m512d  thr  = _mm512_set1_pd(1.0e-8);
			   register __m512d  t0   = v8_0;
			   register __m512d  t1   = v8_0;
			   register __m512d  t2   = v8_0;
			   register __m512d  zeta = v8_0;
			   register __m512d  al_c = v8_0;
			   register __m512d  be_c = v8_0;
			   register __m512d  ga_c = v8_0;
			   register __m512d  alp  = v8_0;
			   register __m512d  bet  = v8_0;
			   register __m512d  gam  = v8_0;
			   __mmask8 k1 = 0x0;
			   __mmask8 k2 = 0x0;
			   __mmask8 k3 = 0x0;
			   __mmask8 k4 = 0x0;
			   k1 = _mm512_cmp_pd_mask(rm.m_vRow9,v8_1,_CMP_NEQ_OQ);
			   t0 = _mm512_mul_pd(rm.m_vRow9,rm.m_vRow9);
			   zeta = _mm512_div_pd(_mm512_sqrt_pd(_mm512_sub_pd(v8_1,t0)));
			   t0 = _mm512_sub_pd(v8_0,rm.vRow8);
			   alp = _mm512_mask_blend_pd(k1,_mm512_atan2_pd(rm.vRow2,rm.vRow1),
			                                   _mm512_atan2_pd(_mm512_mul_pd(rm.vRow7,zeta),
			   				                   _mm512_mul_pd(t0,zeta)));
			   t1    = _mm512_mul_pd(v8_1o2,_mm512_mul_pd(v8_pi,_mm512_sub_pd(v8_1,rm.vRow9)));
			   bet  = _mm512_mask_blend_pd(k1,t1,_mm512_acos_pd(rm.vRow9));
			   gam = _mm512_mask_blend_pd(k1,v8_0,_mm512_atan_pd(_mm512_mul_pd(rm.vRow3,zeta),
			                                                        _mm512_mul_pd(rm.vRow6,zeta)));
			   al_c = alp;
			   be_c = bet;
			   ga_c = gam;
			   k2 = _mm512_cmp_pd_mask(_mm512_abs_pd(alp),thr,_CMP_LT_OQ);
			   alp = _mm512_mask_mov_pd(al_c,k2,v8_0);
			   k3 = _mm512_cmp_pd_mask(_mm512_abs_pd(bet),thr,_CMP_LT_OQ);
			   bet = _mm512_mask_mov_pd(be_c,k3,v8_0);
			   k4 = _mm512_cmp_pd_mask(_mm512_abs_pd(gam),thr,_CMP_LT_OQ);
			   gam = _mm512_mask_mov_pd(ga_c,k4,v8_0);
			   //al_c = alpha;
			   t0 = _mm512_add_pd(alp,v8_2pi);
			   k2 = _mm512_cmp_pd_mask(alp,v8_0,_CMP_LT_OQ);
			   alp = _mm512_mask_mov_pd(al_c,k2,fmod_zmm8r8(t0,v8_2pi));
			   //be_c = beta;
			   t1 = _mm512_add_pd(bet,v8_2pi);
			   k3 = _mm512_cmp_pd_mask(bet,v8_0,_CMP_LT_OQ);
			   bet = _mm512_mask_mov_pd(be_c,k3,fmod_zmm8r8(t1,v8_pi));
			   //ga_c = gamma;
			   t2 = _mm512_add_pd(gam,v8_2pi);
			   k4 = _mm512_cmp_pd_mask(gam,v8_0,_CMP_LT_OQ);
			   gam = _mm512_mask_mov_pd(ga_c,k4,fmod_zmm8r8(t2,v8_2pi));
			   _mm512_storeu_pd(&alpha[0],alp);
			   _mm512_storeu_pd(&beta[0], bet);
			   _mm512_storeu_pd(&gamma[0],gam);
		     }
					
		      
     } // math

} // gms













#endif /* __GMS_ROTATION_KERNELS_AVX512_HPP__*/
