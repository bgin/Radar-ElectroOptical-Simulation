

#ifndef __GMS_DCM_AVX512_HPP__
#define __GMS_DCM_AVX512_HPP__

namespace file_info {

const unsigned int gGMS_DCM_AVX512_MAJOR = 1U;
const unsigned int gGMS_DCM_AVX512_MINOR = 0U;
const unsigned int gGMS_DCM_AVX512_MICRO = 0U;
const unsigned int gGMS_DCM_AVX512_FULLVER =
         1000U*gGMS_DCM_AVX512_MAJOR +
	 100U*gGMS_DCM_AVX512_MINOR +
	 10U*gGMS_DCM_AVX512_MICRO;

const char * const pgGMS_DCM_AVX512_CREATION_DATE = "07-08-2021 10:13 AM +00200 (SAT 07 AUG 2021 GMT+2)";
const char * const pgGMS_DCM_AVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
const char * const pgGMS_DCM_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
const char * const pgGMS_DCM_AVX512_DESCRIPTION   = "Earth to Body Axis Direction-Cosine Matrix [AVX512-based implementation]";
}


#include <immintrin.h>
#include "GMS_config.h"

namespace gms {

        namespace math {


	    namespace {

               const static __m512   v16f32_0 = _mm512_setzero_ps();
	       const static __m512d  v8f64_0  = _mm512_setzero_pd();
	    }


	       struct DCM9x16 __ATTR_ALIGN__(64) {

	             /*
                            Each row represent the expanded element of the Direction Cosine Matrix i.e. A[ij]
                            into 1x16 column vector, thus each vector holds the state of consecutive 16 radian
                            arguments i.e. phi, theta, psi.
                        */
                     __m512 m_vRow1;
		     __m512 m_vRow2;
		     __m512 m_vRow3;
		     __m512 m_vRow4;
		     __m512 m_vRow5;
		     __m512 m_vRow6;
		     __m512 m_vRow7;
		     __m512 m_vRow8;
		     __m512 m_vRow9;

		     /*
                          Default Ctor -- zero initialization
                       */
		    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    DCM9x16() {
                       m_vRow1 = v16f32_0;
		       m_vRow2 = v16f32_0;
		       m_vRow3 = v16f32_0;
		       m_vRow4 = v16f32_0;
		       m_vRow5 = v16f32_0;
		       m_vRow6 = v16f32_0;
		       m_vRow7 = v16f32_0;
		       m_vRow8 = v16f32_0;
		       m_vRow9 = v16f32_0;
		    }


		    /*
                         Build from the other rows i.e. of type: __m512
                      */
		    __ATTR_VECTORCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    DCM9x16( const __m512 vRow1,
		             const __m512 vRow2,
			     const __m512 vRow3,
			     const __m512 vRow4,
			     const __m512 vRow5,
			     const __m512 vRow6,
			     const __m512 vRow7,
			     const __m512 vRow8,
			     const __m512 vRow9) {
                        m_vRow1 = vRow1;
			m_vRow2 = vRow2;
			m_vRow3 = vRow3;
			m_vRow4 = vRow4;
			m_vRow5 = vRow5;
			m_vRow6 = vRow6;
			m_vRow7 = vRow7;
			m_vRow8 = vRow8;
			m_vRow9 = vRow9;

		   }


		   /*
                       Build from the 16-elements single precision arrays.
                       Must be aligned on the 64-byte boundaries.
                    */
		    __ATTR_VECTORCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    DCM9x16( const float * __restrict __ATTR_ALIGN__(64) vRow1,
		             const float * __restrict __ATTR_ALIGN__(64) vRow2,
			     const float * __restrict __ATTR_ALIGN__(64) vRow3,
			     const float * __restrict __ATTR_ALIGN__(64) vRow4,
			     const float * __restrict __ATTR_ALIGN__(64) vRow5,
			     const float * __restrict __ATTR_ALIGN__(64) vRow6,
			     const float * __restrict __ATTR_ALIGN__(64) vRow7,
			     const float * __restrict __ATTR_ALIGN__(64) vRow8,
			     const float * __restrict __ATTR_ALIGN__(64) vRow9) {
                         m_vRow1 = _mm512_load_ps(&vRow1[0]);
			 m_vRow2 = _mm512_load_ps(&vRow2[0]);
			 m_vRow3 = _mm512_load_ps(&vRow3[0]);
			 m_vRow4 = _mm512_load_ps(&vRow4[0]);
                         m_vRow5 = _mm512_load_ps(&vRow5[0]);
			 m_vRow6 = _mm512_load_ps(&vRow6[0]);
			 m_vRow7 = _mm512_load_ps(&vRow7[0]);
			 m_vRow8 = _mm512_load_ps(&vRow8[0]);
			 m_vRow9 = _mm512_load_ps(&vRow9[0]);
		    }

		    /*
                        Copy-Ctor
                     */
		    
		    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32) 
		    DCM9x16(const DCM9x16 &x) {
                        m_vRow1 = x.m_vRow1;
			m_vRow2 = x.m_vRow2;
			m_vRow3 = x.m_vRow3;
			m_vRow4 = x.m_vRow4;
			m_vRow5 = x.m_vRow5;
			m_vRow6 = x.m_vRow6;
			m_vRow7 = x.m_vRow7;
			m_vRow8 = x.m_vRow8;
			m_vRow9 = x.m_vRow9;
		    }

                    __ATTR_VECTORCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    DCM9x16 &
		    DCM_compute_zmm16r4(   const __m512 vPhi,
		                           const __m512 vTheta,
				           const __m512 vPsi) {
                      const __m512 sinRoll  = _mm512_sin_ps(vPhi);
		      const __m512 cosRoll  = _mm512_sin_ps(vPhi);
                      const __m512 sinPitch = _mm512_sin_ps(vTheta);
		      const __m512 cosPitch = _mm512_cos_ps(vTheta);
		      const __m512 sinYaw   = _mm512_sin_ps(vPsi);
		      const __m512 cosYaw   = _mm512_cos_ps(vPsi);
		      const __m512 ct1      = _mm512_mul_ps(sinRoll,sinPitch);
		      const __m512 ct2      = _mm512_mul_ps(cosRoll,sinPitch);
		      const __m512 t1       = _mm512_mul_ps(cosRoll,sinYaw);
		      const __m512 t2       = _mm512_mul_ps(cosRoll,cosYaw);
		      const __m512 t3       = _mm512_mul_ps(sinRoll,sinYaw);
		      const __m512 t4       = _mm512_mul_ps(sinRoll,cosYaw);
		      m_vRow1 = _mm512_mul_ps(cosPitch,cosYaw);
		      m_vRow2 = _mm512_mul_ps(cosPitch,sinYaw);
		      m_vRow3 = _mm512_sub_ps(v16f32_0,sinPitch);
		      m_vRow4 = _mm512_fmsub_ps(ct1,cosYaw,t1);
		      m_vRow5 = _mm512_fmadd_ps(ct1,sinYaw,t2);
		      m_vRow6 = _mm512_mul_ps(sinRoll,cosPitch);
		      m_vRow7 = _mm512_fmadd_ps(ct2,cosYaw,t3);
		      m_vRow8 = _mm512_fmsub_ps(ct2,sinYaw,t4);
		      m_vRow9 = _mm512_mul_ps(cosRoll,cosPitch);
		      return (*this);
		  }

       }; // struct


        struct DCM9x8 __ATTR_ALIGN__(64) {

	             /*
                            Each row represent the expanded element of the Direction Cosine Matrix i.e. A[ij]
                            into 1x8 column vector, thus each vector holds the state of consecutive 8 radian
                            arguments i.e. phi, theta, psi.
                        */
                     __m512d m_vRow1;
		     __m512d m_vRow2;
		     __m512d m_vRow3;
		     __m512d m_vRow4;
		     __m512d m_vRow5;
		     __m512d m_vRow6;
		     __m512d m_vRow7;
		     __m512d m_vRow8;
		     __m512d m_vRow9;

		     /*
                          Default Ctor -- zero initialization
                       */
		    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    DCM9x8() {
                       m_vRow1 = v8f64_0;
		       m_vRow2 = v8f64_0;
		       m_vRow3 = v8f64_0;
		       m_vRow4 = v8f64_0;
		       m_vRow5 = v8f64_0;
		       m_vRow6 = v8f64_0;
		       m_vRow7 = v8f64_0;
		       m_vRow8 = v8f64_0;
		       m_vRow9 = v8f64_0;
		    }


		    /*
                         Build from the other rows i.e. of type: __m512
                      */
		    __ATTR_VECTORCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    DCM9x8(      const __m512d vRow1,
		                 const __m512d vRow2,
				 const __m512d vRow3,
				 const __m512d vRow4,
				 const __m512d vRow5,
				 const __m512d vRow6,
				 const __m512d vRow7,
				 const __m512d vRow8,
				 const __m512d vRow9) {
                        m_vRow1 = vRow1;
			m_vRow2 = vRow2;
			m_vRow3 = vRow3;
			m_vRow4 = vRow4;
			m_vRow5 = vRow5;
			m_vRow6 = vRow6;
			m_vRow7 = vRow7;
			m_vRow8 = vRow8;
			m_vRow9 = vRow9;

		   }


		   /*
                       Build from the 8-elements double precision arrays.
                       Must be aligned on the 64-byte boundaries.
                    */
		    __ATTR_VECTORCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    DCM9x8(      const double * __restrict __ATTR_ALIGN__(64) vRow1,
		                 const double * __restrict __ATTR_ALIGN__(64) vRow2,
				 const double * __restrict __ATTR_ALIGN__(64) vRow3,
				 const double * __restrict __ATTR_ALIGN__(64) vRow4,
				 const double * __restrict __ATTR_ALIGN__(64) vRow5,
				 const double * __restrict __ATTR_ALIGN__(64) vRow6,
				 const double * __restrict __ATTR_ALIGN__(64) vRow7,
				 const double * __restrict __ATTR_ALIGN__(64) vRow8,
				 const double * __restrict __ATTR_ALIGN__(64) vRow9) {
                         m_vRow1 = _mm512_load_pd(&vRow1[0]);
			 m_vRow2 = _mm512_load_pd(&vRow2[0]);
			 m_vRow3 = _mm512_load_pd(&vRow3[0]);
			 m_vRow4 = _mm512_load_pd(&vRow4[0]);
                         m_vRow5 = _mm512_load_pd(&vRow5[0]);
			 m_vRow6 = _mm512_load_pd(&vRow6[0]);
			 m_vRow7 = _mm512_load_pd(&vRow7[0]);
			 m_vRow8 = _mm512_load_pd(&vRow8[0]);
			 m_vRow9 = _mm512_load_pd(&vRow9[0]);
		    }

		    /*
                        Copy-Ctor
                     */
		    
		    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32) 
		    DCM9x8(const DCM9x8 &x) {
                        m_vRow1 = x.m_vRow1;
			m_vRow2 = x.m_vRow2;
			m_vRow3 = x.m_vRow3;
			m_vRow4 = x.m_vRow4;
			m_vRow5 = x.m_vRow5;
			m_vRow6 = x.m_vRow6;
			m_vRow7 = x.m_vRow7;
			m_vRow8 = x.m_vRow8;
			m_vRow9 = x.m_vRow9;
		    }

                    __ATTR_VECTORCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    DCM9x8 &
		    DCM_compute_zmm8r8(const __m512d vPhi,
		                           const __m512d vTheta,
				           const __m512d vPsi) {
                      const __m512d sinRoll  = _mm512_sin_pd(vPhi);
		      const __m512d cosRoll  = _mm512_sin_pd(vPhi);
                      const __m512d sinPitch = _mm512_sin_pd(vTheta);
		      const __m512d cosPitch = _mm512_cos_pd(vTheta);
		      const __m512d sinYaw   = _mm512_sin_pd(vPsi);
		      const __m512d cosYaw   = _mm512_cos_pd(vPsi);
		      const __m512d ct1      = _mm512_mul_pd(sinRoll,sinPitch);
		      const __m512d ct2      = _mm512_mul_pd(cosRoll,sinPitch);
		      const __m512d t1       = _mm512_mul_pd(cosRoll,sinYaw);
		      const __m512d t2       = _mm512_mul_pd(cosRoll,cosYaw);
		      const __m512d t3       = _mm512_mul_pd(sinRoll,sinYaw);
		      const __m512d t4       = _mm512_mul_pd(sinRoll,cosYaw);
		      m_vRow1 = _mm512_mul_pd(cosPitch,cosYaw);
		      m_vRow2 = _mm512_mul_pd(cosPitch,sinYaw);
		      m_vRow3 = _mm512_sub_pd(v16f32_0,sinPitch);
		      m_vRow4 = _mm512_fmsub_pd(ct1,cosYaw,t1);
		      m_vRow5 = _mm512_fmadd_pd(ct1,sinYaw,t2);
		      m_vRow6 = _mm512_mul_pd(sinRoll,cosPitch);
		      m_vRow7 = _mm512_fmadd_pd(ct2,cosYaw,t3);
		      m_vRow8 = _mm512_fmsub_pd(ct2,sinYaw,t4);
		      m_vRow9 = _mm512_mul_pd(cosRoll,cosPitch);
		      return (*this);
		  }

       }; // struct

       

   } // math

} // gms



#endif /*__GMS_DCM_AVX512_HPP__*/
