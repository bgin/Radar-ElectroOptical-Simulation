

#ifndef __GMS_DCM_AVX2_HPP__
#define __GMS_DCM_AVX2_HPP__

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

const unsigned int gGMS_DCM_AVX2_MAJOR = 1U;
const unsigned int gGMS_DCM_AVX2_MINOR = 0U;
const unsigned int gGMS_DCM_AVX2_MICRO = 0U;
const unsigned int gGMS_DCM_AVX2_FULLVER =
         1000U*gGMS_DCM_AVX2_MAJOR +
	 100U*gGMS_DCM_AVX2_MINOR +
	 10U*gGMS_DCM_AVX2_MICRO;

const char * const pgGMS_DCM_AVX512_CREATION_DATE = "08-08-2021 09:13 AM +00200 (SUN 08 AUG 2021 GMT+2)";
const char * const pgGMS_DCM_AVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
const char * const pgGMS_DCM_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
const char * const pgGMS_DCM_AVX512_DESCRIPTION   = "Earth to Body Axis Direction-Cosine Matrix [AVX2-based implementation]";
}


#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"


namespace gms {

           namespace math {

	           namespace {

                    const __m256   v8f32_0 = _mm256_setzero_ps();
		    const __m256d  v4f32_0 = _mm256_setzero_pd();
	      }


	         struct DCMatYMM8r4 __ATTR_ALIGN__(32) {

                          //  Each row represent the expanded element of the Direction Cosine Matrix i.e. A[ij]
                          //  into 1x8 column vector, thus each vector holds the state of consecutive 8 radian
                          //  arguments i.e. phi, theta, psi.
			__m256 m_vRow1;
			__m256 m_vRow2;
			__m256 m_vRow3;
			__m256 m_vRow4;
			__m256 m_vRow5;
			__m256 m_vRow6;
			__m256 m_vRow7;
			__m256 m_vRow8;
			__m256 m_vRow9;

			// Default Ctor -- zero initialization.
		       __ATTR_ALWAYS_INLINE__
		       __ATTR_HOT__
		       __ATTR_ALIGN__(32)
		       DCMatYMM8r4() {
                           m_vRow1 = v8f32_0;
			   m_vRow2 = v8f32_0;
                           m_vRow3 = v8f32_0;
			   m_vRow4 = v8f32_0;
			   m_vRow5 = v8f32_0;
			   m_vRow6 = v8f32_0;
			   m_vRow7 = v8f32_0;
			   m_vRow8 = v8f32_0;
			   m_vRow8 = v8f32_0;
		       }

		        /*
                         Build from the other rows i.e. of type: __m256
                      */
                      __ATTR_VECTORCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
                      DCMatYMM8r4(    const __m256 vRow1,
		                      const __m256 vRow2,
				      const __m256 vRow3,
				      const __m256 vRow4,
				      const __m256 vRow5,
				      const __m256 vRow6,
				      const __m256 vRow7,
				      const __m256 vRow8,
				      const __m256 vRow9) {
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
                    DCMatYMM8r4(    const float * __restrict __ATTR_ALIGN__(32) vRow1,
		                    const float * __restrict __ATTR_ALIGN__(32) vRow2,
				    const float * __restrict __ATTR_ALIGN__(32) vRow3,
				    const float * __restrict __ATTR_ALIGN__(32) vRow4,
				    const float * __restrict __ATTR_ALIGN__(32) vRow5,
				    const float * __restrict __ATTR_ALIGN__(32) vRow6,
				    const float * __restrict __ATTR_ALIGN__(32) vRow7,
				    const float * __restrict __ATTR_ALIGN__(32) vRow8,
				    const float * __restrict __ATTR_ALIGN__(32) vRow9) {
                         m_vRow1 = _mm256_load_ps(&vRow1[0]);
			 m_vRow2 = _mm256_load_ps(&vRow2[0]);
                         m_vRow3 = _mm256_load_ps(&vRow3[0]);
			 m_vRow4 = _mm256_load_ps(&vRow4[0]);
			 m_vRow5 = _mm256_load_ps(&vRow5[0]);
			 m_vRow6 = _mm256_load_ps(&vRow6[0]);
			 m_vRow7 = _mm256_load_ps(&vRow7[0]);
			 m_vRow8 = _mm256_load_ps(&vRow8[0]);
			 m_vRow9 = _mm256_load_ps(&vRow9[0]);
		 }

                  /*
                        Copy-Ctor
                  */
		    
		   __ATTR_ALWAYS_INLINE__
		   __ATTR_HOT__
		   __ATTR_ALIGN__(32)
		   DCMatYMM8r4(const DCMatYMM8r4 &x) {
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
		 DCMatYMM8r4 &
		 DCMat_compute_ymm8r4(const __m256 vPhi,
		                      const __m256 vTheta,
				      const __m256 vPsi) {
                      const __m256 sinRoll  = _mm256_sin_ps(vPhi);
		      const __m256 cosRoll  = _mm256_sin_ps(vPhi);
                      const __m256 sinPitch = _mm256_sin_ps(vTheta);
		      const __m256 cosPitch = _mm256_cos_ps(vTheta);
		      const __m256 sinYaw   = _mm256_sin_ps(vPsi);
		      const __m256 cosYaw   = _mm256_cos_ps(vPsi);
		      const __m256 ct1      = _mm256_mul_ps(sinRoll,sinPitch);
		      const __m256 ct2      = _mm256_mul_ps(cosRoll,sinPitch);
		      const __m256 t1       = _mm256_mul_ps(cosRoll,sinYaw);
		      const __m256 t2       = _mm256_mul_ps(cosRoll,cosYaw);
		      const __m256 t3       = _mm256_mul_ps(sinRoll,sinYaw);
		      const __m256 t4       = _mm256_mul_ps(sinRoll,cosYaw);
		      m_vRow1 = _mm256_mul_ps(cosPitch,cosYaw);
		      m_vRow2 = _mm256_mul_ps(cosPitch,sinYaw);
		      m_vRow3 = _mm256_sub_ps(v16f32_0,sinPitch);
		      m_vRow4 = _mm256_fmsub_ps(ct1,cosYaw,t1);
		      m_vRow5 = _mm256_fmadd_ps(ct1,sinYaw,t2);
		      m_vRow6 = _mm256_mul_ps(sinRoll,cosPitch);
		      m_vRow7 = _mm256_fmadd_ps(ct2,cosYaw,t3);
		      m_vRow8 = _mm256_fmsub_ps(ct2,sinYaw,t4);
		      m_vRow9 = _mm256_mul_ps(cosRoll,cosPitch);
                      return (*this);
		}



		 
		      

	     }; // struct


	       struct DCMatYMM4r8 __ATTR_ALIGN__(32) {

                          //  Each row represent the expanded element of the Direction Cosine Matrix i.e. A[ij]
                          //  into 1x4 column vector, thus each vector holds the state of consecutive 4 radian
                          //  arguments i.e. phi, theta, psi.
			__m256d m_vRow1;
			__m256d m_vRow2;
			__m256d m_vRow3;
			__m256d m_vRow4;
			__m256d m_vRow5;
			__m256d m_vRow6;
			__m256d m_vRow7;
			__m256d m_vRow8;
			__m256d m_vRow9;

			// Default Ctor -- zero initialization.
		       __ATTR_ALWAYS_INLINE__
		       __ATTR_HOT__
		       __ATTR_ALIGN__(32)
		       DCMatYMM4r8() {
                           m_vRow1 = v4f64_0;
			   m_vRow2 = v4f64_0;
                           m_vRow3 = v4f64_0;
			   m_vRow4 = v4f64_0;
			   m_vRow5 = v4f64_0;
			   m_vRow6 = v4f64_0;
			   m_vRow7 = v4f64_0;
			   m_vRow8 = v4f64_0;
			   m_vRow8 = v4f64_0;
		       }

		        /*
                         Build from the other rows i.e. of type: __m256d
                      */
                      __ATTR_VECTORCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
                      DCMatYMM4r8(    const __m256d vRow1,
		                      const __m256d vRow2,
				      const __m256d vRow3,
				      const __m256d vRow4,
				      const __m256d vRow5,
				      const __m256d vRow6,
				      const __m256d vRow7,
				      const __m256d vRow8,
				      const __m256d vRow9) {
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
                    DCMatYMM4r8(    const double * __restrict __ATTR_ALIGN__(32) vRow1,
		                    const double * __restrict __ATTR_ALIGN__(32) vRow2,
				    const double * __restrict __ATTR_ALIGN__(32) vRow3,
				    const double * __restrict __ATTR_ALIGN__(32) vRow4,
				    const double * __restrict __ATTR_ALIGN__(32) vRow5,
				    const double * __restrict __ATTR_ALIGN__(32) vRow6,
				    const double * __restrict __ATTR_ALIGN__(32) vRow7,
				    const double * __restrict __ATTR_ALIGN__(32) vRow8,
				    const double * __restrict __ATTR_ALIGN__(32) vRow9) {
                         m_vRow1 = _mm256_load_pd(&vRow1[0]);
			 m_vRow2 = _mm256_load_pd(&vRow2[0]);
                         m_vRow3 = _mm256_load_pd(&vRow3[0]);
			 m_vRow4 = _mm256_load_pd(&vRow4[0]);
			 m_vRow5 = _mm256_load_pd(&vRow5[0]);
			 m_vRow6 = _mm256_load_pd(&vRow6[0]);
			 m_vRow7 = _mm256_load_pd(&vRow7[0]);
			 m_vRow8 = _mm256_load_pd(&vRow8[0]);
			 m_vRow9 = _mm256_load_pd(&vRow9[0]);
		 }

                  /*
                        Copy-Ctor
                  */
		    
		   __ATTR_ALWAYS_INLINE__
		   __ATTR_HOT__
		   __ATTR_ALIGN__(32)
		   DCMatYMM4r8(const DCMatYMM4r8 &x) {
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
		 DCMatYMM4r8 &
		 DCM_compute_ymm4r8(  const __m256d vPhi,
		                      const __m256d vTheta,
				      const __m256d vPsi) {
                      const __m256d sinRoll  = _mm256_sin_pd(vPhi);
		      const __m256d cosRoll  = _mm256_sin_pd(vPhi);
                      const __m256d sinPitch = _mm256_sin_pd(vTheta);
		      const __m256d cosPitch = _mm256_cos_pd(vTheta);
		      const __m256d sinYaw   = _mm256_sin_pd(vPsi);
		      const __m256d cosYaw   = _mm256_cos_pd(vPsi);
		      const __m256d ct1      = _mm256_mul_pd(sinRoll,sinPitch);
		      const __m256d ct2      = _mm256_mul_pd(cosRoll,sinPitch);
		      const __m256d t1       = _mm256_mul_pd(cosRoll,sinYaw);
		      const __m256d t2       = _mm256_mul_pd(cosRoll,cosYaw);
		      const __m256d t3       = _mm256_mul_pd(sinRoll,sinYaw);
		      const __m256d t4       = _mm256_mul_pd(sinRoll,cosYaw);

		      m_vRow1 = _mm256_mul_pd(cosPitch,cosYaw);
		      m_vRow2 = _mm256_mul_pd(cosPitch,sinYaw);
		      m_vRow3 = _mm256_sub_pd(v16f32_0,sinPitch);
		      m_vRow4 = _mm256_fmsub_pd(ct1,cosYaw,t1);
		      m_vRow5 = _mm256_fmadd_pd(ct1,sinYaw,t2);
		      m_vRow6 = _mm256_mul_pd(sinRoll,cosPitch);
		      m_vRow7 = _mm256_fmadd_pd(ct2,cosYaw,t3);
		      m_vRow8 = _mm256_fmsub_pd(ct2,sinYaw,t4);
		      m_vRow9 = _mm256_mul_pd(cosRoll,cosPitch);
                      return (*this);
		}



		 
		      

	     }; // struct
	     

   } // math

} // gms






#endif /*__GMS_DCM_AVX2_HPP__*/
