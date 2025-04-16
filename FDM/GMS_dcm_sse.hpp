

#ifndef __GMS_DCM_SSE_HPP__
#define __GMS_DCM_SSE_HPP__ 080820211713

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

namespace file_info {

const unsigned int gGMS_DCM_SSE_MAJOR = 1U;
const unsigned int gGMS_DCM_SSE_MINOR = 0U;
const unsigned int gGMS_DCM_SSE_MICRO = 0U;
const unsigned int gGMS_DCM_SSE_FULLVER =
         1000U*gGMS_DCM_SSE_MAJOR +
	 100U*gGMS_DCM_SSE_MINOR +
	 10U*gGMS_DCM_SSE_MICRO;

const char * const pgGMS_DCM_SSE_CREATION_DATE = "08-08-2021 17:13 AM +00200 (SUN 08 AUG 2021 GMT+2)";
const char * const pgGMS_DCM_SSE_BUILD_DATE    = __DATE__ ":" __TIME__;
const char * const pgGMS_DCM_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
const char * const pgGMS_DCM_SSE_DESCRIPTION   = "Earth to Body Axis Direction-Cosine Matrix [SSE-based implementation]";
}


#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"


namespace gms {

           namespace math {

	           namespace {

                    const __m128   v4f32_0 = _mm_setzero_ps();
		    const __m128d  v2f64_0 = _mm_setzero_pd();
	      }
	      
#define DCMATXMM4R4_BUILD_BLOCK(phi,theta,psi)        \
const __m128 sinRoll  = _mm_sin_ps(phi);              \
const __m128 cosRoll  = _mm_sin_ps(phi);              \
const __m128 sinPitch = _mm_sin_ps(theta);            \
const __m128 cosPitch = _mm_cos_ps(theta);            \
const __m128 sinYaw   = _mm_sin_ps(psi);              \
const __m128 cosYaw   = _mm_cos_ps(psi);              \
const __m128 ct1      = _mm_mul_ps(sinRoll,sinPitch); \
const __m128 ct2      = _mm_mul_ps(cosRoll,sinPitch); \
const __m128 t1       = _mm_mul_ps(cosRoll,sinYaw);   \
const __m128 t2       = _mm_mul_ps(cosRoll,cosYaw);   \
const __m128 t3       = _mm_mul_ps(sinRoll,sinYaw);   \
const __m128 t4       = _mm_mul_ps(sinRoll,cosYaw);   
#endif


#define DCMATXMM2R8_BUILD_BLOCK(phi,theta,psi)         \
const __m128d sinRoll  = _mm_sin_pd(phi);              \
const __m128d cosRoll  = _mm_sin_pd(phi);              \
const __m128d sinPitch = _mm_sin_pd(theta);            \
const __m128d cosPitch = _mm_cos_pd(theta);            \
const __m128d sinYaw   = _mm_sin_pd(psi);              \
const __m128d cosYaw   = _mm_cos_pd(psi);              \
const __m128d ct1      = _mm_mul_pd(sinRoll,sinPitch); \
const __m128d ct2      = _mm_mul_pd(cosRoll,sinPitch); \
const __m128d t1       = _mm_mul_pd(cosRoll,sinYaw);   \
const __m128d t2       = _mm_mul_pd(cosRoll,cosYaw);   \
const __m128d t3       = _mm_mul_pd(sinRoll,sinYaw);   \
const __m128d t4       = _mm_mul_pd(sinRoll,cosYaw);   
#endif

	         struct DCMatXMM4r4 __ATTR_ALIGN__(16) {

                          //  Each row represent the expanded element of the Direction Cosine Matrix i.e. A[ij]
                          //  into 1x4 column vector, thus each vector holds the state of consecutive 4 radian
                          //  arguments i.e. phi, theta, psi.
			__m128 m_vRow1;
			__m128 m_vRow2;
			__m128 m_vRow3;
			__m128 m_vRow4;
			__m128 m_vRow5;
			__m128 m_vRow6;
			__m128 m_vRow7;
			__m128 m_vRow8;
			__m128 m_vRow9;

			// Default Ctor -- zero initialization.
		       __ATTR_ALWAYS_INLINE__
		       __ATTR_HOT__
		       __ATTR_ALIGN__(32)
		       DCMatXMM4r4() {
                           m_vRow1 = v4f32_0;
			   m_vRow2 = v4f32_0;
                           m_vRow3 = v4f32_0;
			   m_vRow4 = v4f32_0;
			   m_vRow5 = v4f32_0;
			   m_vRow6 = v4f32_0;
			   m_vRow7 = v4f32_0;
			   m_vRow8 = v4f32_0;
			   m_vRow8 = v4f32_0;
		       }

		        /*
                         Build from the other rows i.e. of type: __m128
                      */
                      __ATTR_VECTORCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
                      DCMatXMM4r4(    const __m128 vRow1,
		                      const __m128 vRow2,
				      const __m128 vRow3,
				      const __m128 vRow4,
				      const __m128 vRow5,
				      const __m128 vRow6,
				      const __m128 vRow7,
				      const __m128 vRow8,
				      const __m128 vRow9) {
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
                       Build from the 8-elements single precision arrays.
                       Must be aligned on the 32-byte boundaries.
                    */
		    __ATTR_VECTORCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
                    DCMatXMM4r4(    const float * __restrict __ATTR_ALIGN__(16) vRow1,
		                    const float * __restrict __ATTR_ALIGN__(16) vRow2,
				    const float * __restrict __ATTR_ALIGN__(16) vRow3,
				    const float * __restrict __ATTR_ALIGN__(16) vRow4,
				    const float * __restrict __ATTR_ALIGN__(16) vRow5,
				    const float * __restrict __ATTR_ALIGN__(16) vRow6,
				    const float * __restrict __ATTR_ALIGN__(16) vRow7,
				    const float * __restrict __ATTR_ALIGN__(16) vRow8,
				    const float * __restrict __ATTR_ALIGN__(16) vRow9) {
                         m_vRow1 = _mm_load_ps(&vRow1[0]);
			 m_vRow2 = _mm_load_ps(&vRow2[0]);
                         m_vRow3 = _mm_load_ps(&vRow3[0]);
			 m_vRow4 = _mm_load_ps(&vRow4[0]);
			 m_vRow5 = _mm_load_ps(&vRow5[0]);
			 m_vRow6 = _mm_load_ps(&vRow6[0]);
			 m_vRow7 = _mm_load_ps(&vRow7[0]);
			 m_vRow8 = _mm_load_ps(&vRow8[0]);
			 m_vRow9 = _mm_load_ps(&vRow9[0]);
		 }

                  /*
                        Copy-Ctor
                  */
		    
		   __ATTR_ALWAYS_INLINE__
		   __ATTR_HOT__
		   __ATTR_ALIGN__(32)
		   DCMatXMM4r4(const DCMatXMM4r4 &x) {
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
		 DCMatXMM4r4 &
		 DCMat_compute_xmm4r4(const __m128 vPhi,
		                      const __m128 vTheta,
				      const __m128 vPsi) {
	              DCMATXMM4R4_BUILD_BLOCK(vPhi,vTheta,vPsi)
                      m_vRow1 = _mm_mul_ps(cosPitch,cosYaw);
		      m_vRow2 = _mm_mul_ps(cosPitch,sinYaw);
		      m_vRow3 = _mm_sub_ps(v4f32_0,sinPitch);
		      m_vRow4 = _mm_fmsub_ps(ct1,cosYaw,t1);
		      m_vRow5 = _mm_fmadd_ps(ct1,sinYaw,t2);
		      m_vRow6 = _mm_mul_ps(sinRoll,cosPitch);
		      m_vRow7 = _mm_fmadd_ps(ct2,cosYaw,t3);
		      m_vRow8 = _mm_fmsub_ps(ct2,sinYaw,t4);
		      m_vRow9 = _mm_mul_ps(cosRoll,cosPitch);
                      return (*this);
		}
		
		
		
		                           



		 
		      

	     }; // struct
	     
	     
	         __ATTR_VECTORCALL__
                 __ATTR_ALWAYS_INLINE__
		 __ATTR_HOT__
		 __ATTR_ALIGN__(32)
		 void DCMat_compute_xmm4r4(DCMatXMM4r4 &dcmat,
		                           const __m128 vPhi,
		                           const __m128 vTheta,
				           const __m128 vPsi) {
		       
		      DCMATXMM4R4_BUILD_BLOCK(vPhi,vTheta,vPsi)
		      dcmat.m_vRow1 = _mm_mul_ps(cosPitch,cosYaw);
		      dcmat.m_vRow2 = _mm_mul_ps(cosPitch,sinYaw);
		      dcmat.m_vRow3 = _mm_sub_ps(v4f32_0,sinPitch);
		      dcmat.m_vRow4 = _mm_fmsub_ps(ct1,cosYaw,t1);
		      dcmat.m_vRow5 = _mm_fmadd_ps(ct1,sinYaw,t2);
		      dcmat.m_vRow6 = _mm_mul_ps(sinRoll,cosPitch);
		      dcmat.m_vRow7 = _mm_fmadd_ps(ct2,cosYaw,t3);
		      dcmat.m_vRow8 = _mm_fmsub_ps(ct2,sinYaw,t4);
		      dcmat.m_vRow9 = _mm_mul_ps(cosRoll,cosPitch);		           
	        }


	       struct DCMatXMM2r8 __ATTR_ALIGN__(16) {

                          //  Each row represent the expanded element of the Direction Cosine Matrix i.e. A[ij]
                          //  into 1x4 column vector, thus each vector holds the state of consecutive 4 radian
                          //  arguments i.e. phi, theta, psi.
			__m128d m_vRow1;
			__m128d m_vRow2;
			__m128d m_vRow3;
			__m128d m_vRow4;
			__m128d m_vRow5;
			__m128d m_vRow6;
			__m128d m_vRow7;
			__m128d m_vRow8;
			__m128d m_vRow9;

			// Default Ctor -- zero initialization.
		       __ATTR_ALWAYS_INLINE__
		       __ATTR_HOT__
		       __ATTR_ALIGN__(32)
		       DCMatXMM2r8() {
                           m_vRow1 = v2f64_0;
			   m_vRow2 = v2f64_0;
                           m_vRow3 = v2f64_0;
			   m_vRow4 = v2f64_0;
			   m_vRow5 = v2f64_0;
			   m_vRow6 = v2f64_0;
			   m_vRow7 = v2f64_0;
			   m_vRow8 = v2f64_0;
			   m_vRow8 = v2f64_0;
		       }

		        /*
                         Build from the other rows i.e. of type: __m128d
                      */
                      __ATTR_VECTORCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
                      DCMatXMM2r8(    const __m128d vRow1,
		                      const __m128d vRow2,
				      const __m128d vRow3,
				      const __m128d vRow4,
				      const __m128d vRow5,
				      const __m128d vRow6,
				      const __m128d vRow7,
				      const __m128d vRow8,
				      const __m128d vRow9) {
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
                       Build from the 4-elements double precision arrays.
                       Must be aligned on the 32-byte boundaries.
                    */
		    __ATTR_VECTORCALL__
                    __ATTR_ALWAYS_INLINE__
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
                    DCMatXMM2r8(    const double * __restrict __ATTR_ALIGN__(16) vRow1,
		                    const double * __restrict __ATTR_ALIGN__(16) vRow2,
				    const double * __restrict __ATTR_ALIGN__(16) vRow3,
				    const double * __restrict __ATTR_ALIGN__(16) vRow4,
				    const double * __restrict __ATTR_ALIGN__(16) vRow5,
				    const double * __restrict __ATTR_ALIGN__(16) vRow6,
				    const double * __restrict __ATTR_ALIGN__(16) vRow7,
				    const double * __restrict __ATTR_ALIGN__(16) vRow8,
				    const double * __restrict __ATTR_ALIGN__(16) vRow9) {
                         m_vRow1 = _mm_load_pd(&vRow1[0]);
			 m_vRow2 = _mm_load_pd(&vRow2[0]);
                         m_vRow3 = _mm_load_pd(&vRow3[0]);
			 m_vRow4 = _mm_load_pd(&vRow4[0]);
			 m_vRow5 = _mm_load_pd(&vRow5[0]);
			 m_vRow6 = _mm_load_pd(&vRow6[0]);
			 m_vRow7 = _mm_load_pd(&vRow7[0]);
			 m_vRow8 = _mm_load_pd(&vRow8[0]);
			 m_vRow9 = _mm_load_pd(&vRow9[0]);
		 }

                  /*
                        Copy-Ctor
                  */
		    
		   __ATTR_ALWAYS_INLINE__
		   __ATTR_HOT__
		   __ATTR_ALIGN__(32)
		   DCMatXMM2r8(const DCMatXMM2r8 &x) {
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
		 DCMatXMM2r8 &
		 DCM_compute_xmm2r8(  const __m128d vPhi,
		                      const __m128d vTheta,
				      const __m128d vPsi) {
				      
		      DCMATXMM2R8_BUILD_BLOCK(vPhi,vTheta,vPsi)
                      m_vRow1 = _mm_mul_pd(cosPitch,cosYaw);
		      m_vRow2 = _mm_mul_pd(cosPitch,sinYaw);
		      m_vRow3 = _mm_sub_pd(v2f64_0,sinPitch);
		      m_vRow4 = _mm_fmsub_pd(ct1,cosYaw,t1);
		      m_vRow5 = _mm_fmadd_pd(ct1,sinYaw,t2);
		      m_vRow6 = _mm_mul_pd(sinRoll,cosPitch);
		      m_vRow7 = _mm_fmadd_pd(ct2,cosYaw,t3);
		      m_vRow8 = _mm_fmsub_pd(ct2,sinYaw,t4);
		      m_vRow9 = _mm_mul_pd(cosRoll,cosPitch);
                      return (*this);
		}



		 
		      

	     }; // struct
	     
	     
	         __ATTR_VECTORCALL__
                 __ATTR_ALWAYS_INLINE__
		 __ATTR_HOT__
		 __ATTR_ALIGN__(32)
		 void DCMat_compute_xmm2r8(DCMatXMM2r8 &dcmat,
		                           const __m128d vPhi,
		                           const __m128d vTheta,
				           const __m128d vPsi) {
		       
		      DCMATXMM2R8_BUILD_BLOCK(vPhi,vTheta,vPsi)
		      dcmat.m_vRow1 = _mm_mul_ps(cosPitch,cosYaw);
		      dcmat.m_vRow2 = _mm_mul_ps(cosPitch,sinYaw);
		      dcmat.m_vRow3 = _mm_sub_ps(v2f64_0,sinPitch);
		      dcmat.m_vRow4 = _mm_fmsub_ps(ct1,cosYaw,t1);
		      dcmat.m_vRow5 = _mm_fmadd_ps(ct1,sinYaw,t2);
		      dcmat.m_vRow6 = _mm_mul_ps(sinRoll,cosPitch);
		      dcmat.m_vRow7 = _mm_fmadd_ps(ct2,cosYaw,t3);
		      dcmat.m_vRow8 = _mm_fmsub_ps(ct2,sinYaw,t4);
		      dcmat.m_vRow9 = _mm_mul_ps(cosRoll,cosPitch);		           
	        }

   } // math

} // gms






#endif /*__GMS_DCM_SSE_HPP__*/
