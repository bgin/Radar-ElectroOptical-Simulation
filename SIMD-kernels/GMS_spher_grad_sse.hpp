

#ifndef __GMS_SPHER_GRAD_SSE_HPP__
#define __GMS_SPHER_GRAD_SSE_HPP__ 170520220848


/*LICENSE:
*
*The source code is in the public domain and not licensed or under
*copyright. The information and software may be used freely by the public.
*As required by 17 U.S.C. 403, third parties producing copyrighted works
*consisting predominantly of the material produced by U.S. government
*agencies must provide notice with such work(s) identifying the U.S.
*Government material incorporated and stating that such material is not
*subject to copyright protection.
*
*Derived works shall not identify themselves in a manner that implies an
*endorsement by or an affiliation with the Naval Research Laboratory.
*
*RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
*SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
*RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
*OF RECIPIENT IN THE USE OF THE SOFTWARE.
@@Modified by Bernard Gingold, on 17-05-2022 08:48 +00200 (TUE 17 MAY 2022 08:48 GMT+2)
  contact: beniekg@gmail.com
*/

namespace file_info {

 const unsigned int GMS_SPHER_GRAD_SSE_MAJOR = 1U;
 const unsigned int GMS_SPHER_GRAD_SSE_MINOR = 0U;
 const unsigned int GMS_SPHER_GRAD_SSE_MICRO = 0U;
 const unsigned int GMS_SPHER_GRAD_SSE_FULLVER =
  1000U*GMS_SPHER_GRAD_SSE_MAJOR+100U*GMS_SPHER_GRAD_SSE_MINOR+10U*GMS_SPHER_GRAD_SSE_MICRO;
 const char * const GMS_SPHER_GRAD_SSE_CREATION_DATE = "17-05-2022 08:48 +00200 (TUE 17 MAY 2022 08:48 GMT+2)";
 const char * const GMS_SPHER_GRAD_SSE_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const GMS_SPHER_GRAD_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const GMS_SPHER_GRAD_SSE_SYNOPSIS      = "SSE based spherical coordinates hessians and jacobians functions (vectorized)."


}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_simd_utils.hpp"


namespace gms {

        namespace math {

	
/**CALCSPHERINVJACOBCPP  A C++-only implementations of a function for
 *          computing the Jacobian of of a 3D Cartesian point with respect
 *          to spherical azimuth and elevation. See the Matlab equivalent
 *          for more comments.
 *
 *July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *May  2022 Bernard Gingold, manually vectorized.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
                    
                      __ATTR_ALWAYS_INLINE__
		     static inline
		      void  spher_inv_jac_xmm2r8_a(double * __restrict __ATTR_ALIGN__(16) J0,
		                                   double * __restrict __ATTR_ALIGN__(16) J1,
						   double * __restrict __ATTR_ALIGN__(16) J2,
						   double * __restrict __ATTR_ALIGN__(16) J3,
						   double * __restrict __ATTR_ALIGN__(16) J4,
						   double * __restrict __ATTR_ALIGN__(16) J5,
						   double * __restrict __ATTR_ALIGN__(16) J6,
						   double * __restrict __ATTR_ALIGN__(16) J7,
						   double * __restrict __ATTR_ALIGN__(16) J8,
						   const __m128d x,
						   const __m128d y,
						   const __m128d z,
						   const int32_t sysType) {

                         const __m128d sinaz = _mm_sin_pd(y);
			 const __m128d cosaz = _mm_cos_pd(y);
			 const __m128d sinel = _mm_sin_pd(z);
			 const __m128d cosel = _mm_cos_pd(z);

                         if(sysType==0) {
			    const __m128d nx = xmm2r8_negate(x):
                            //dx/dr
			    _mm_store_pd(&J0[0],_mm_mul_pd(cosaz,cosel));
			    //dy/dr
			    _mm_store_pd(&J1[0], _mm_mul_pd(cosel,sinaz));
			    //dz/dr
			    _mm_store_pd(&J2[0],sinel);
			    //dx/dAz
			    _mm_store_pd(&J3[0],_mm_mul_pd(nx,
			                       _mm_mul_pd(cosel,sinaz)));
			    //dy/dAz
			    _mm_store_pd(&J4[0], _mm_mul_pd(x,_mm_mul_pd(cosaz,cosel)));
			    //dz/dAz
			    _mm_store_pd(&J5[0], _mm_setzero_pd());
			    //dx/dEl
			    _mm_store_pd(&J6[0],_mm_mul_pd(nx,
			                       _mm_mul_pd(cosaz,sinel)));
			    //dy/dEl
			    _mm_store_pd(&J7[0],_mm_mul_pd(nx,
			                       _mm_mul_pd(sinaz,sinel)));
			    //dz/dEl
			    _mm_store_pd(&J8[0],_mm_mul_pd(x,cosel));
			 }
			 else if(sysType==1) {
                             const __m128d nx    = xmm2r8_negate(x):
			     const __m128d cazel = _mm_mul_pd(cosaz,cosel);
			     _mm_store_pd(&J0[0],_mm_mul_pd(cosel,sinaz));
			     _mm_store_pd(&J1[0],sinel);
			     _mm_store_pd(&J2[0],cazel);
			     _mm_store_pd(&J3[0], _mm_mul_pd(x,cazel));
			     _mm_store_pd(&J4[0], _mm_setzero_pd());
			     _mm_store_pd(&J5[0], _mm_mul_pd(nx,_mm_mul_pd(cosel,sinaz)));
			     _mm_store_pd(&J6[0], _mm_mul_pd(nx,_mm_mul_pd(sinaz,sinel)));
			     _mm_store_pd(&J7[0], _mm_mul_pd(x,cosel));
			     _mm_store_pd(&J8[0], _mm_mul_pd(nx,_mm_mul_pd(cosaz,sinel)));
			 }
			 else if(sysType==3) {
                             const __m128d nx    = xmm2r8_negate(x):
			     const __m128d sc    = _mm_mul_pd(sinaz,cosel);
			     const __m128d cc    = _mm_mul_pd(cosaz,cosel);
			     _mm_store_pd(&J0[0],sc);
			     _mm_store_pd(&J1[0],cc);
			     _mm_store_pd(&J2[0],sinel);
			     _mm_store_pd(&J3[0],_mm_mul_pd(x,cc));
			     _mm_store_pd(&J4[0],_mm_mul_pd(nx,sc));
			     _mm_store_pd(&J5[0],_mm_setzero_pd());
			     _mm_store_pd(&J6[0],_mm_mul_pd(nx,_mm_mul_pd(sinaz,sinel)));
			     _mm_store_pd(&J7[0],_mm_mul_pd(nx,_mm_mul_pd(cosaz,sinel)));
			     _mm_store_pd(&J8[0],_mm_mul_pd(x,cosel));
			 }
			 else { // sysType==2
                             const __m128d nx    = xmm2r8_negate(x):
			     const __m128d ss    = _mm_mul_pd(sinaz,sinel);
			     const __m128d cs    = _mm_mul_pd(cosaz,sinel);
			     _mm_store_pd(&J0[0],cs);
			     _mm_store_pd(&J1[0],ss);
			     _mm_store_pd(&J2[0],cosel);
			     _mm_store_pd(&J3[0],_mm_mul_pd(nx,ss));
			     _mm_store_pd(&J4[0],_mm_mul_pd(x,cs));
			     _mm_store_pd(&J5[0],_mm_setzero_pd());
			     _mm_store_pd(&J6[0],_mm_mul_pd(x,_mm_mul_pd(cosaz,cosel)));
			     _mm_store_pd(&J7[0],_mm_mul_pd(x,_mm_mul_pd(cosel,sinaz)));
			     _mm_store_pd(&J8[0],_mm_mul_pd(nx,sinel));
			 }
		   }


		    
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void  spher_inv_jac_xmm2r8_u(double * __restrict J0,
		                                   double * __restrict J1,
						   double * __restrict J2,
						   double * __restrict J3,
						   double * __restrict J4,
						   double * __restrict J5,
						   double * __restrict J6,
						   double * __restrict J7,
						   double * __restrict J8,
						   const __m128d x,
						   const __m128d y,
						   const __m128d z,
						   const int32_t sysType) {

                         const __m128d sinaz = _mm_sin_pd(y);
			 const __m128d cosaz = _mm_cos_pd(y);
			 const __m128d sinel = _mm_sin_pd(z);
			 const __m128d cosel = _mm_cos_pd(z);

                         if(sysType==0) {
			    const __m128d nx = xmm2r8_negate(x):
                            //dx/dr
			    _mm_storeu_pd(&J0[0],_mm_mul_pd(cosaz,cosel));
			    //dy/dr
			    _mm_storeu_pd(&J1[0], _mm_mul_pd(cosel,sinaz));
			    //dz/dr
			    _mm_storeu_pd(&J2[0],sinel);
			    //dx/dAz
			    _mm_storeu_pd(&J3[0],_mm_mul_pd(nx,
			                       _mm_mul_pd(cosel,sinaz)));
			    //dy/dAz
			    _mm_storeu_pd(&J4[0], _mm_mul_pd(x,_mm_mul_pd(cosaz,cosel)));
			    //dz/dAz
			    _mm_storeu_pd(&J5[0], _mm_setzero_pd());
			    //dx/dEl
			    _mm_storeu_pd(&J6[0],_mm_mul_pd(nx,
			                       _mm_mul_pd(cosaz,sinel)));
			    //dy/dEl
			    _mm_storeu_pd(&J7[0],_mm_mul_pd(nx,
			                       _mm_mul_pd(sinaz,sinel)));
			    //dz/dEl
			    _mm_storeu_pd(&J8[0],_mm_mul_pd(x,cosel));
			 }
			 else if(sysType==1) {
                             const __m128d nx    = xmm2r8_negate(x):
			     const __m128d cazel = _mm_mul_pd(cosaz,cosel);
			     _mm_storeu_pd(&J0[0],_mm_mul_pd(cosel,sinaz));
			     _mm_storeu_pd(&J1[0],sinel);
			     _mm_storeu_pd(&J2[0],cazel);
			     _mm_storeu_pd(&J3[0], _mm_mul_pd(x,cazel));
			     _mm_storeu_pd(&J4[0], _mm_setzero_pd());
			     _mm_storeu_pd(&J5[0], _mm_mul_pd(nx,_mm_mul_pd(cosel,sinaz)));
			     _mm_storeu_pd(&J6[0], _mm_mul_pd(nx,_mm_mul_pd(sinaz,sinel)));
			     _mm_storeu_pd(&J7[0], _mm_mul_pd(x,cosel));
			     _mm_storeu_pd(&J8[0], _mm_mul_pd(nx,_mm_mul_pd(cosaz,sinel)));
			 }
			 else if(sysType==3) {
                             const __m128d nx    = xmm2r8_negate(x):
			     const __m128d sc    = _mm_mul_pd(sinaz,cosel);
			     const __m128d cc    = _mm_mul_pd(cosaz,cosel);
			     _mm_storeu_pd(&J0[0],sc);
			     _mm_storeu_pd(&J1[0],cc);
			     _mm_storeu_pd(&J2[0],sinel);
			     _mm_storeu_pd(&J3[0],_mm_mul_pd(x,cc));
			     _mm_storeu_pd(&J4[0],_mm_mul_pd(nx,sc));
			     _mm_storeu_pd(&J5[0],_mm_setzero_pd());
			     _mm_storeu_pd(&J6[0],_mm_mul_pd(nx,_mm_mul_pd(sinaz,sinel)));
			     _mm_storeu_pd(&J7[0],_mm_mul_pd(nx,_mm_mul_pd(cosaz,sinel)));
			     _mm_storeu_pd(&J8[0],_mm_mul_pd(x,cosel));
			 }
			 else { // sysType==2
                             const __m128d nx    = xmm2r8_negate(x):
			     const __m128d ss    = _mm_mul_pd(sinaz,sinel);
			     const __m128d cs    = _mm_mul_pd(cosaz,sinel);
			     _mm_storeu_pd(&J0[0],cs);
			     _mm_storeu_pd(&J1[0],ss);
			     _mm_storeu_pd(&J2[0],cosel);
			     _mm_storeu_pd(&J3[0],_mm_mul_pd(nx,ss));
			     _mm_storeu_pd(&J4[0],_mm_mul_pd(x,cs));
			     _mm_storeu_pd(&J5[0],_mm_setzero_pd());
			     _mm_storeu_pd(&J6[0],_mm_mul_pd(x,_mm_mul_pd(cosaz,cosel)));
			     _mm_storeu_pd(&J7[0],_mm_mul_pd(x,_mm_mul_pd(cosel,sinaz)));
			     _mm_storeu_pd(&J8[0],_mm_mul_pd(nx,sinel));
			 }
		   }




	             


		     
                      __ATTR_ALWAYS_INLINE__
		     static inline
		      void  spher_inv_jac_xmm4r4_a(float * __restrict __ATTR_ALIGN__(16) J0,
		                                    float * __restrict __ATTR_ALIGN__(16) J1,
						    float * __restrict __ATTR_ALIGN__(16) J2,
						    float * __restrict __ATTR_ALIGN__(16) J3,
						    float * __restrict __ATTR_ALIGN__(16) J4,
						    float * __restrict __ATTR_ALIGN__(16) J5,
						    float * __restrict __ATTR_ALIGN__(16) J6,
						    float * __restrict __ATTR_ALIGN__(16) J7,
						    float * __restrict __ATTR_ALIGN__(16) J8,
						    const __m128 x,
						    const __m128 y,
						    const __m128 z,
						    const int32_t sysType) {

                         const __m128 sinaz = _mm_sin_ps(y);
			 const __m128 cosaz = _mm_cos_ps(y);
			 const __m128 sinel = _mm_sin_ps(z);
			 const __m128 cosel = _mm_cos_ps(z);

                         if(sysType==0) {
			    const __m128 nx = xmm4r4_negate(x):
                            //dx/dr
			    _mm_store_ps(&J0[0],_mm_mul_ps(cosaz,cosel));
			    //dy/dr
			    _mm_store_ps(&J1[0],_mm_mul_ps(cosel,sinaz));
			    //dz/dr
			    _mm_store_ps(&J2[0],sinel);
			    //dx/dAz
			    _mm_store_ps(&J3[0],_mm_mul_ps(nx,
			                       _mm_mul_ps(cosel,sinaz)));
			    //dy/dAz
			    _mm_store_ps(&J4[0],_mm_mul_ps(x,_mm_mul_ps(cosaz,cosel)));
			    //dz/dAz
			    _mm_store_ps(&J5[0],_mm_setzero_ps());
			    //dx/dEl
			    _mm_store_ps(&J6[0],_mm_mul_ps(nx,
			                       _mm_mul_ps(cosaz,sinel)));
			    //dy/dEl
			    _mm_store_ps(&J7[0],_mm_mul_ps(nx,
			                       _mm_mul_ps(sinaz,sinel)));
			    //dz/dEl
			    _mm_store_ps(&J8[0],_mm_mul_ps(x,cosel));
			 }
			 else if(sysType==1) {
                             const __m128 nx    = xmm4r4_negate(x):
			     const __m128 cazel = _mm_mul_ps(cosaz,cosel);
			     _mm_store_ps(&J0[0],_mm_mul_ps(cosel,sinaz));
			     _mm_store_ps(&J1[0],sinel);
			     _mm_store_ps(&J2[0],cazel);
			     _mm_store_ps(&J3[0],_mm_mul_ps(x,cazel));
			     _mm_store_ps(&J4[0],_mm_setzero_ps());
			     _mm_store_ps(&J5[0],_mm_mul_ps(nx,_mm_mul_ps(cosel,sinaz)));
			     _mm_store_ps(&J6[0],_mm_mul_ps(nx,_mm_mul_ps(sinaz,sinel)));
			     _mm_store_ps(&J7[0],_mm_mul_ps(x,cosel));
			     _mm_store_ps(&J8[0],_mm_mul_ps(nx,_mm_mul_ps(cosaz,sinel)));
			 }
			 else if(sysType==3) {
                             const __m128 nx    = xmm4r4_negate(x):
			     const __m128 sc    = _mm_mul_ps(sinaz,cosel);
			     const __m128 cc    = _mm_mul_ps(cosaz,cosel);
			     _mm_store_ps(&J0[0],sc);
			     _mm_store_ps(&J1[0],cc);
			     _mm_store_ps(&J2[0],sinel);
			     _mm_store_ps(&J3[0],_mm_mul_ps(x,cc));
			     _mm_store_ps(&J4[0],_mm_mul_ps(nx,sc));
			     _mm_store_ps(&J5[0],_mm_setzero_ps());
			     _mm_store_ps(&J6[0],_mm_mul_ps(nx,_mm_mul_ps(sinaz,sinel)));
			     _mm_store_ps(&J7[0],_mm_mul_ps(nx,_mm_mul_ps(cosaz,sinel)));
			     _mm_store_ps(&J8[0],_mm_mul_ps(x,cosel));
			 }
			 else { // sysType==2
                             const __m128 nx    = xmm4r4_negate(x):
			     const __m128 ss    = _mm_mul_ps(sinaz,sinel);
			     const __m128 cs    = _mm_mul_ps(cosaz,sinel);
			     _mm_store_ps(&J0[0],cs);
			     _mm_store_ps(&J1[0],ss);
			     _mm_store_ps(&J2[0],cosel);
			     _mm_store_ps(&J3[0],_mm_mul_ps(nx,ss));
			     _mm_store_ps(&J4[0],_mm_mul_ps(x,cs));
			     _mm_store_ps(&J5[0],_mm_setzero_ps());
			     _mm_store_ps(&J6[0],_mm_mul_ps(x,_mm_mul_ps(cosaz,cosel)));
			     _mm_store_ps(&J7[0],_mm_mul_ps(x,_mm_mul_ps(cosel,sinaz)));
			     _mm_store_ps(&J8[0],_mm_mul_ps(nx,sinel));
			 }
		   }


		    
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void  spher_inv_jac_xmm4r4_u(float * __restrict J0,
		                                    float * __restrict J1,
						    float * __restrict J2,
						    float * __restrict J3,
						    float * __restrict J4,
						    float * __restrict J5,
						    float * __restrict J6,
						    float * __restrict J7,
						    float * __restrict J8,
						    const __m128 x,
						    const __m128 y,
						    const __m128 z,
						    const int32_t sysType) {

                         const __m128 sinaz = _mm_sin_ps(y);
			 const __m128 cosaz = _mm_cos_ps(y);
			 const __m128 sinel = _mm_sin_ps(z);
			 const __m128 cosel = _mm_cos_ps(z);

                         if(sysType==0) {
			    const __m128 nx = xmm4r4_negate(x):
                            //dx/dr
			    _mm_storeu_ps(&J0[0],_mm_mul_ps(cosaz,cosel));
			    //dy/dr
			    _mm_storeu_ps(&J1[0],_mm_mul_ps(cosel,sinaz));
			    //dz/dr
			    _mm_storeu_ps(&J2[0],sinel);
			    //dx/dAz
			    _mm_storeu_ps(&J3[0],_mm_mul_ps(nx,
			                       _mm_mul_ps(cosel,sinaz)));
			    //dy/dAz
			    _mm_storeu_ps(&J4[0],_mm_mul_ps(x,_mm_mul_ps(cosaz,cosel)));
			    //dz/dAz
			    _mm_storeu_ps(&J5[0],_mm_setzero_ps());
			    //dx/dEl
			    _mm_storeu_ps(&J6[0],_mm_mul_ps(nx,
			                       _mm_mul_ps(cosaz,sinel)));
			    //dy/dEl
			    _mm_storeu_ps(&J7[0],_mm_mul_ps(nx,
			                       _mm_mul_ps(sinaz,sinel)));
			    //dz/dEl
			    _mm_storeu_ps(&J8[0],_mm_mul_ps(x,cosel));
			 }
			 else if(sysType==1) {
                             const __m128 nx    = xmm4r4_negate(x):
			     const __m128 cazel = _mm_mul_ps(cosaz,cosel);
			     _mm_storeu_ps(&J0[0],_mm_mul_ps(cosel,sinaz));
			     _mm_storeu_ps(&J1[0],sinel);
			     _mm_storeu_ps(&J2[0],cazel);
			     _mm_storeu_ps(&J3[0],_mm_mul_ps(x,cazel));
			     _mm_storeu_ps(&J4[0],_mm_setzero_ps());
			     _mm_storeu_ps(&J5[0],_mm_mul_ps(nx,_mm_mul_ps(cosel,sinaz)));
			     _mm_storeu_ps(&J6[0],_mm_mul_ps(nx,_mm_mul_ps(sinaz,sinel)));
			     _mm_storeu_ps(&J7[0],_mm_mul_ps(x,cosel));
			     _mm_storeu_ps(&J8[0],_mm_mul_ps(nx,_mm_mul_ps(cosaz,sinel)));
			 }
			 else if(sysType==3) {
                             const __m128 nx    = xmm4r4_negate(x):
			     const __m128 sc    = _mm_mul_ps(sinaz,cosel);
			     const __m128 cc    = _mm_mul_ps(cosaz,cosel);
			     _mm_storeu_ps(&J0[0],sc);
			     _mm_storeu_ps(&J1[0],cc);
			     _mm_storeu_ps(&J2[0],sinel);
			     _mm_storeu_ps(&J3[0],_mm_mul_ps(x,cc));
			     _mm_storeu_ps(&J4[0],_mm_mul_ps(nx,sc));
			     _mm_storeu_ps(&J5[0],_mm_setzero_ps());
			     _mm_storeu_ps(&J6[0],_mm_mul_ps(nx,_mm_mul_ps(sinaz,sinel)));
			     _mm_storeu_ps(&J7[0],_mm_mul_ps(nx,_mm_mul_ps(cosaz,sinel)));
			     _mm_storeu_ps(&J8[0],_mm_mul_ps(x,cosel));
			 }
			 else { // sysType==2
                             const __m128 nx    = xmm4r4_negate(x):
			     const __m128 ss    = _mm_mul_ps(sinaz,sinel);
			     const __m128 cs    = _mm_mul_ps(cosaz,sinel);
			     _mm_storeu_ps(&J0[0],cs);
			     _mm_storeu_ps(&J1[0],ss);
			     _mm_storeu_ps(&J2[0],cosel);
			     _mm_storeu_ps(&J3[0],_mm_mul_ps(nx,ss));
			     _mm_storeu_ps(&J4[0],_mm_mul_ps(x,cs));
			     _mm_storeu_ps(&J5[0],_mm_setzero_ps());
			     _mm_storeu_ps(&J6[0],_mm_mul_ps(x,_mm_mul_ps(cosaz,cosel)));
			     _mm_storeu_ps(&J7[0],_mm_mul_ps(x,_mm_mul_ps(cosel,sinaz)));
			     _mm_storeu_ps(&J8[0],_mm_mul_ps(nx,sinel));
			 }
		   }


		   /*
                        Different definitions with an array of type __m128d/__m128 i.e. AoSoA
                        instead of SoA.
                    */

		     
                       __ATTR_ALWAYS_INLINE__
		     static inline
		      void  spher_inv_jac_xmm2r8(__m128d * __restrict __ATTR_ALIGN__(16) J, //flatten 9x1 array
		                                 const __m128d x,
						 const __m128d y,
						 const __m128d z,
						 const int32_t sysType) {

                         const __m128d sinaz = _mm_sin_pd(y);
			 const __m128d cosaz = _mm_cos_pd(y);
			 const __m128d sinel = _mm_sin_pd(z);
			 const __m128d cosel = _mm_cos_pd(z);

                         if(sysType==0) {
			    const __m128d nx = xmm2r8_negate(x):
                            //dx/dr
			    J[0] = _mm_mul_pd(cosaz,cosel);
			    //dy/dr
			    J[1] = _mm_mul_pd(cosel,sinaz);
			    //dz/dr
			    J[2] = sinel;
			    //dx/dAz
			    J[3] = _mm_mul_pd(nx,
			                       _mm_mul_pd(cosel,sinaz));
			    //dy/dAz
			    J[4] = _mm_mul_pd(x,_mm_mul_pd(cosaz,cosel));
			    //dz/dAz
			    J[5] = _mm_setzero_pd();
			    //dx/dEl
			    J[6] = _mm_mul_pd(nx,
			                       _mm_mul_pd(cosaz,sinel));
			    //dy/dEl
			    J[7] = _mm_mul_pd(nx,
			                       _mm_mul_pd(sinaz,sinel));
			    //dz/dEl
			    J[8] = _mm_mul_pd(x,cosel);
			 }
			 else if(sysType==1) {
                             const __m128d nx    = xmm2r8_negate(x):
			     const __m128d cazel = _mm_mul_pd(cosaz,cosel);
			     J[0] = _mm_mul_pd(cosel,sinaz);
			     J[1] = sinel;
			     J[2] = cazel;
			     J[3] = _mm_mul_pd(x,cazel);
			     J[4] = _mm_setzero_pd();
			     J[5] = _mm_mul_pd(nx,_mm_mul_pd(cosel,sinaz));
			     J[6] = _mm_mul_pd(nx,_mm_mul_pd(sinaz,sinel));
			     J[7] = _mm_mul_pd(x,cosel);
			     J[8] = _mm_mul_pd(nx,_mm_mul_pd(cosaz,sinel));
			 }
			 else if(sysType==3) {
                             const __m128d nx    = xmm2r8_negate(x):
			     const __m128d sc    = _mm_mul_pd(sinaz,cosel);
			     const __m128d cc    = _mm_mul_pd(cosaz,cosel);
			     J[0] = sc;
			     J[1] = cc;
			     J[2] = sinel;
			     J[3] = _mm_mul_pd(x,cc);
			     J[4] = _mm_mul_pd(nx,sc);
			     J[5] = _mm_setzero_pd();
			     J[6] = _mm_mul_pd(nx,_mm_mul_pd(sinaz,sinel));
			     J[7] = _mm_mul_pd(nx,_mm_mul_pd(cosaz,sinel));
			     J[8] = _mm_mul_pd(x,cosel);
			 }
			 else { // sysType==2
                             const __m128d nx    = xmm2r8_negate(x):
			     const __m128d ss    = _mm_mul_pd(sinaz,sinel);
			     const __m128d cs    = _mm_mul_pd(cosaz,sinel);
			     J[0] = cs;
			     J[1] = ss;
			     J[2] = cosel;
			     J[3] = _mm_mul_pd(nx,ss);
			     J[4] = _mm_mul_pd(x,cs);
			     J[5] = _mm_setzero_pd();
			     J[6] = _mm_mul_pd(x,_mm_mul_pd(cosaz,cosel));
			     J[7] = _mm_mul_pd(x,_mm_mul_pd(cosel,sinaz));
			     J[8] = _mm_mul_pd(nx,sinel);
			 }
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		     static inline
		      void  spher_inv_jac_xmm4r4(__m128 * __restrict __ATTR_ALIGN__(16) J,
		                                  const __m128 x,
						  const __m128 y,
						  const __m128 z,
						  const int32_t sysType) {

                         const __m128 sinaz = _mm_sin_ps(y);
			 const __m128 cosaz = _mm_cos_ps(y);
			 const __m128 sinel = _mm_sin_ps(z);
			 const __m128 cosel = _mm_cos_ps(z);

                         if(sysType==0) {
			    const __m128 nx = xmm4r4_negate(x):
                            //dx/dr
			    J[0] = _mm_mul_ps(cosaz,cosel);
			    //dy/dr
			    J[1] = _mm_mul_ps(cosel,sinaz);
			    //dz/dr
			    J[2] = sinel;
			    //dx/dAz
			    J[3] = _mm_mul_ps(nx,
			                       _mm_mul_ps(cosel,sinaz));
			    //dy/dAz
			    J[4] = _mm_mul_ps(x,_mm_mul_ps(cosaz,cosel));
			    //dz/dAz
			    J[5] = _mm_setzero_ps();
			    //dx/dEl
			    J[6] = _mm_mul_ps(nx,
			                       _mm_mul_ps(cosaz,sinel));
			    //dy/dEl
			    J[7] = _mm_mul_ps(nx,
			                       _mm_mul_ps(sinaz,sinel));
			    //dz/dEl
			    J[8] = _mm_mul_ps(x,cosel);
			 }
			 else if(sysType==1) {
                             const __m128 nx    = xmm4r4_negate(x):
			     const __m128 cazel = _mm_mul_ps(cosaz,cosel);
			     J[0] = _mm_mul_ps(cosel,sinaz);
			     J[1] = sinel;
			     J[2] = cazel;
			     J[3] = _mm_mul_ps(x,cazel);
			     J[4] = _mm_setzero_ps();
			     J[5] = _mm_mul_ps(nx,_mm_mul_ps(cosel,sinaz));
			     J[6] = _mm_mul_ps(nx,_mm_mul_ps(sinaz,sinel));
			     J[7] = _mm_mul_ps(x,cosel);
			     J[8] = _mm_mul_ps(nx,_mm_mul_ps(cosaz,sinel));
			 }
			 else if(sysType==3) {
                             const __m128 nx    = xmm4r4_negate(x):
			     const __m128 sc    = _mm_mul_ps(sinaz,cosel);
			     const __m128 cc    = _mm_mul_ps(cosaz,cosel);
			     J[0] = sc;
			     J[1] = cc;
			     J[2] = sinel;
			     J[3] = _mm_mul_ps(x,cc);
			     J[4] = _mm_mul_ps(nx,sc);
			     J[5] = _mm_setzero_ps();
			     J[6] = _mm_mul_ps(nx,_mm_mul_ps(sinaz,sinel));
			     J[7] = _mm_mul_ps(nx,_mm_mul_ps(cosaz,sinel));
			     J[8] = _mm_mul_ps(x,cosel);
			 }
			 else { // sysType==2
                             const __m128 nx    = xmm4r4_negate(x):
			     const __m128 ss    = _mm_mul_ps(sinaz,sinel);
			     const __m128 cs    = _mm_mul_ps(cosaz,sinel);
			     J[0] = cs;
			     J[1] = ss;
			     J[2] = cosel;
			     J[3] = _mm_mul_ps(nx,ss);
			     J[4] = _mm_mul_ps(x,cs);
			     J[5] = _mm_setzero_ps();
			     J[6] = _mm_mul_ps(x,_mm_mul_ps(cosaz,cosel));
			     J[7] = _mm_mul_ps(x,_mm_mul_ps(cosel,sinaz));
			     J[8] = _mm_mul_ps(nx,sinel);
			 }
		   }

/**SPHERANGHESSIANCPP A C++-only implementation of a function for
 *          computing the Hessian of azimuth and elevation in spherical
 *          coordinates.  See the Matlab equivalent for more comments.
 *
 *June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *May  2022 Bernard Gingold, manually vectorized.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
	              

                   
                      __ATTR_ALWAYS_INLINE__
		     static inline
		      void spher_ang_hess_xmm2r8(__m128d * __restrict __ATTR_ALIGN__(16) H,
		                                 const __m128d G0,
						 const __m128d G1,
						 const __m128d G2,
						 const int32_t sysType) {
                         const __m128d x  = G0;
			 const __m128d x2 = _mm_mul_pd(x,x);
			 const __m128d y  = G1;
			 const __m128d y2 = _mm_mul_pd(y,y);
			 const __m128d z  = G2;
			 const __m128d z2 = _mm_mul_pd(z,z);
			 const __m128d x4 = _mm_mul_pd(x2,x2);
			 
			 if(sysType==0) {
			    const __m128d _2 = _mm_set1_pd(2.0);
			    const __m128d _3 = _mm_set1_pd(3.0);
			    const __m128d _0 = _mm_setzero_pd();
                            __m128d rxy,r4,r4xy,rxy3;
			    __m128d t0,t1,t2,t3,t4,t5,c0;
			    r4xy = _mm_add_pd(x2,y2);
			    rxy  = _mm_sqrt_pd(r4xy);
			    r4   = _mm_add_pd(r4xy,z2);
			    r4   = _mm_mul_pd(r4,r4);
			    r4xy = _mm_mul_pd(r4xy,r4xy);
			    rxy3 = _mm_mul_pd(rxy,_mm_mul_pd(rxy,rxy));
			    t0   = _mm_mul_pd(rxy3,r4);
			    c0   = _mm_mul_pd(rxy,r4);
			    H[0] = _mm_div_pd(_mm_mul_pd(_2,_mm_mul_pd(x,y)),r4xy);
			    H[4] = xmm2r8_negate(H[0]);
			    H[8] = _0;
			    H[3] = _mm_div_pd(_mm_sub_pd(y2,x2),r4xy);
			    H[1] = H[3];
			    H[6] = _0;
			    H[2] = H[6];
			    H[7] = _0;
			    H[5] = H[7];
			    t1   = _mm_mul_pd(y2,_mm_add_pd(y2,z2));
			    t2   = _mm_fmadd_pd(x2,y2,_mm_mul_pd(_2,x4));
			    H[9] = _mm_div_pd(_mm_mul_pd(z,_mm_sub_pd(t1,t2)),t0);
			    t3   = _mm_mul_pd(x2,_mm_sub_pd(z2,y2));
			    t4   = _mm_sub_pd(x4,_mm_mul_pd(_2,_mm_mul_pd(y2,y2)));
			    H[13]= _mm_div_pd(xmm2r8_negate(_mm_mul_pd(z,_mm_add_pd(t3,t4))),t0);
			    H[17]= _mm_div_pd(xmm2r8_negate(_mm_mul_pd(_2,_mm_mul_pd(rxy,z))),r4);
			    t1   = _mm_fmadd_pd(_mm_mul_pd(_3,rxy),rxy,z2);
			    t2   = _mm_mul_pd(x,_mm_mul_pd(y,z));
			    H[12]= _mm_div_pd(_mm_mul_pd(t1,t2),t0);
			    H[10]= H[12];
			    t3   = _mm_mul_pd(x,_mm_sub_pd(_mm_add_pd(x2,y2),z2));
			    H[15]= _mm_div_pd(xmm2r8_negate(t3),c0);
			    H[11]= H[15];
			    t4   = _mm_mul_pd(y,_mm_sub_pd(_mm_add_pd(x2,y2),z2));
			    H[16]= _mm_div_pd(xmm2r8_negate(t4),c0);
			    H[14]= H[16];
			 }
			 else if(sysType==1) {
			    const __m128d _2 = _mm_set1_pd(2.0);
			    const __m128d _3 = _mm_set1_pd(3.0);
			    const __m128d _0 = _mm_setzero_pd();
                            __m128d x4,rxz,r4xz,rxz3,r4;
			    __m128d c0,c1,c2,c3,c4,c5;
			    x4   = _mm_mul_pd(x2,x2);
			    r4xz = _mm_add_pd(x2,z2);
			    rxz  = _mm_sqrt_pd(r4xz);
			    rxz3 = _mm_mul_pd(rxz,_mm_mul_pd(rxz,rxz));
			    r4   = _mm_add_pd(r4xz,y2);
			    r4xz = _mm_mul_pd(r4xz,r4xz);
			    r4   = _mm_mul_pd(r4,r4);
			    c0   = _mm_mul_pd(rxz3,r4);
			    c5   = _mm_mul_pd(rxz,r4);
			    H[0] = _mm_div_pd(xmm2r8_negate(_mm_mul_pd(_2,_mm_mul_pd(x,z))),r4xz);
			    H[4] = _0;
			    H[8] = _mm_div_pd(_mm_mul_pd(_2,_mm_mul_pd(x,z)),r4xz);
			    H[3] = _0;
			    H[1] = H[3];
			    H[6] = _mm_div_pd(_mm_sub_pd(x2,z2),r4xz);
			    H[2] = H[6];
			    H[7] = _0;
			    H[5] = H[7];
			    c1   = _mm_mul_pd(z2,_mm_add_pd(y2,z2));
			    c2   = _mm_fmadd_pd(_2,x4,_mm_mul_pd(x2,z2));
			    H[9] = _mm_div_pd(_mm_mul_pd(y,_mm_sub_pd(c2,c1)),c0);
			    H[13]= xmm2r8_negate(_mm_div_pd(_mm_mul_pd(_2,_mm_mul_pd(y,rxz))),r4);
			    c3   = _mm_mul_pd(x2,_mm_mul_pd(_mm_add_pd(y,z),_mm_sub_pd(y,z)));
			    c4   = _mm_fmsub_pd(x4,_2,_mm_mul_pd(z2,z2));
			    H[17]= xmm2r8_negate(_mm_div_pd(_mm_mul_pd(y,_mm_add_pd(c3,c4))),c0);
			    H[12]= xmm2r8_negate(_mm_div_pd(_mm_mul_pd(_mm_add_pd(_mm_sub_pd(x2,y2),z2)),c5));
			    H[10]= H[12];
			    c1   = _mm_fmadd_pd(_3,x2,_mm_fmadd_pd(_3,z2,y2));
			    c2   = _mm_mul_pd(x,_mm_mul_pd(y,z));
			    H[15]= _mm_div_pd(_mm_mul_pd(c1,c2),c0);
			    H[11]= H[15];
			    c3   = _mm_mul_pd(z,_mm_add_pd(_mm_sub_pd(x2,y2),z2));
			    H[16]= xmm2r8_negate(_mm_div_pd(c3,c5));
			    H[14]= H[16];
			 }
			 else if(sysType==3) {
			    const __m128d _2 = _mm_set1_pd(2.0);
			    const __m128d _3 = _mm_set1_pd(3.0);
			    const __m128d _0 = _mm_setzero_pd();
			    __m128d rxy,r4,r4xy,rxy3;
			    __m128d c0,c1,c2,c3,c4,c5;
			    r4xy  = _mm_add_pd(x2,y2);
			    rxy   = _mm_sqrt_pd(r4xy);
			    r4    = _mm_add_pd(r4xy,z2);
			    r4    = _mm_mul_pd(r4,r4);
			    r4xy  = _mm_mul_pd(r4xy,r4xy);
			    rxy3  = _mm_mul_pd(rxy,_mm_mul_pd(rxy,rxy));
			    c0    = _mm_mul_pd(rxy3,r4);
			    c5    = _mm_mul_pd(rxy,r4);
			    H[0]  = xmm2r8_negate(_mm_div_pd(_mm_mul_pd(_2,_mm_mul_pd(x,y)),r4xy));
			    H[4]  = H[0];
			    H[8]  = _0;
			    H[3]  = _mm_div_pd(_mm_mul_pd(_mm_sub_pd(x,y),
								_mm_add_pd(x,y)),r4xy);
			    H[1]  = H[3];
			    H[6]  = H[0];
			    H[2]  = H[6];
			    H[7]  = _0;
			    H[5]  = H[7];
			    c1    = _mm_mul_pd(y2,_mm_add_pd(y2,z2));
			    c2    = _mm_fmadd_pd(_2,x4,_mm_mul_pd(x2,y2));
			    H[9]  = _mm_div_pd(_mm_mul_pd(z,_mm_sub_pd(c2,c1)),c0);
			    c3    = _mm_mul_pd(x2,_mm_sub_pd(z2,y2));
			    c4    = _mm_fmsub_pd(_2,_mm_mul_pd(y2,y2),x4);
			    H[13] = xmm2r8_negate(_mm_div_pd(
							_mm_mul_pd(z,_mm_add_pd(c3,c4))),c0);
			    H[17] = xmm2r8_negate(_mm_div_pd(_mm_mul_pd(_2,
								_mm_mul_pd(rxy,z)),r4));
			    c1    = _mm_mul_pd(x,_mm_mul_pd(y,z));
			    c2    = _mm_fmadd_pd(_3,_mm_mul_pd(rxy,rxy),z2);
			    H[12] = _mm_div_pd(_mm_mul_pd(c1,c2),c0);
			    H[10] = H[12];
			    c3    = _mm_mul_pd(x,_mm_sub_pd(
			                                  _mm_add_pd(x2,y2),z2));
			    H[15] = xmm2r8_negate(_mm_div_pd(c3,c5));
			    H[11] = H[15];
			    c4    = _mm_mul_pd(y,_mm_sub_pd(
			                                  _mm_add_pd(x2,y2),z2));
			    H[16] = xmm2r8_negate(_mm_div_pd(c4,c5));
			    H[14] = H[16];
				
			 }
			 else { // sysType == 2
                            const __m128d _2 = _mm_set1_pd(2.0);
			    const __m128d _3 = _mm_set1_pd(3.0);
			    const __m128d _0 = _mm_setzero_pd();
			    __m128d r4xy,r4,rxy,x4,y4,rxy3;
			    __m128d c0,c1,c2,c3,c4,c5;
			    r4xy  = _mm_add_pd(x2,y2);
			    r4    = _mm_add_pd(r4xy,z2);
			    r4    = _mm_mul_pd(r4,r4);
			    rxy   = _mm_sqrt_pd(r4xy);
			    r4xy  = _mm_mul_pd(r4xy,r4xy);
			    rxy3  = _mm_mul_pd(rxy,_mm_mul_pd(rxy,rxy));
			    x4    = _mm_mul_pd(x2,x2);
			    y4    = _mm_mul_pd(y2,y2);
			    c0    = _mm_mul_pd(rxy3,r4);
			    c5    = _mm_mul_pd(rxy,r4);
			    H[0]  = _mm_div_pd(_mm_mul_pd(_2,
			                          _mm_mul_pd(x,y)),r4xy);
			    H[4]  = xmm2r8_negate(H[0]);
			    H[8]  = _0;
			    H[3]  = _mm_div_pd(_mm_sub_pd(y2,x2),r4xy);
			    H[1]  = H[0];
			    H[6]  = _0;
			    H[2]  = H[6];
			    H[7]  = _0;
			    H[5]  = H[7];
			    c1    = _mm_fmsub_pd(xmm2r8_negate(_2),x4,_mm_mul_pd(x2,y2));
			    c2    = _mm_fmadd_pd(y2,z2,y4);
			    H[9]  = _mm_div_pd(_mm_mul_pd(z,_mm_add_pd(c1,c2)),c0);
			    c3    = _mm_mul_pd(x2,_mm_sub_pd(z2,y2));
			    c4    = _mm_sub_pd(x4,_mm_mul_pd(_2,y4));
			    H[13] = _mm_div_pd(_mm_mul_pd(z,_mm_add_pd(c3,c4)),c0);
			    H[17] = _mm_div_pd(_mm_mul_pd(_2,_mm_mul_pd(rxy,z)),r4);
			    c1    = _mm_fmadd_pd(_3,_mm_mul_pd(rxy,rxy),z2);
			    c2    = _mm_mul_pd(x,_mm_mul_pd(y,z));
			    H[12] = xmm2r8_negate(_mm_div_pd(_mm_mul_pd(c1,c2),c0));
			    H[10] = H[12];
			    H[15] = _mm_div_pd(_mm_mul_pd(x,_mm_sub_pd(
			                                          _mm_add_pd(x2,y2),z2)),c5);
			    H[11] = H[15];
			    H[16] = _mm_div_pd(_mm_mul_pd(y,_mm_sub_pd(
			                                          _mm_add_pd(x2,y2),z2)),c5);
			    H[14] = H[16];
			}
		   }


		    
                      __ATTR_ALWAYS_INLINE__
		     static inline
		      void spher_ang_hess_xmm4r4(__m128 * __restrict __ATTR_ALIGN__(16) H,
		                                  const __m128 G0,
						  const __m128 G1,
						  const __m128 G2,
						  const int32_t sysType) {
                         const __m128 x  = G0;
			 const __m128 x2 = _mm_mul_ps(x,x);
			 const __m128 y  = G1;
			 const __m128 y2 = _mm_mul_ps(y,y);
			 const __m128 z  = G2;
			 const __m128 z2 = _mm_mul_ps(z,z);
			 const __m128 x4 = _mm_mul_ps(x2,x2);
			 
			 if(sysType==0) {
			    const __m128 _2 = _mm_set1_ps(2.0f);
			    const __m128 _3 = _mm_set1_ps(3.0f);
			    const __m128 _0 = _mm_setzero_ps();
                            __m128 rxy,r4,r4xy,rxy3;
			    __m128 t0,t1,t2,t3,t4,t5,c0;
			    r4xy = _mm_add_ps(x2,y2);
			    rxy  = _mm_sqrt_ps(r4xy);
			    r4   = _mm_add_ps(r4xy,z2);
			    r4   = _mm_mul_ps(r4,r4);
			    r4xy = _mm_mul_ps(r4xy,r4xy);
			    rxy3 = _mm_mul_ps(rxy,_mm_mul_pd(rxy,rxy));
			    t0   = _mm_mul_ps(rxy3,r4);
			    c0   = _mm_mul_ps(rxy,r4);
			    H[0] = _mm_div_ps(_mm_mul_ps(_2,_mm_mul_ps(x,y)),r4xy);
			    H[4] = xmm4r4_negate(H[0]);
			    H[8] = _0;
			    H[3] = _mm_div_ps(_mm_sub_ps(y2,x2),r4xy);
			    H[1] = H[3];
			    H[6] = _0;
			    H[2] = H[6];
			    H[7] = _0;
			    H[5] = H[7];
			    t1   = _mm_mul_ps(y2,_mm_add_ps(y2,z2));
			    t2   = _mm_fmadd_ps(x2,y2,_mm_mul_ps(_2,x4));
			    H[9] = _mm_div_ps(_mm_mul_ps(z,_mm_sub_ps(t1,t2)),t0);
			    t3   = _mm_mul_ps(x2,_mm_sub_ps(z2,y2));
			    t4   = _mm_sub_ps(x4,_mm_mul_ps(_2,_mm_mul_ps(y2,y2)));
			    H[13]= _mm_div_ps(xmm4r4_negate(_mm_mul_ps(z,_mm_add_ps(t3,t4))),t0);
			    H[17]= _mm_div_ps(xmm4r4_negate(_mm_mul_ps(_2,_mm_mul_ps(rxy,z))),r4);
			    t1   = _mm_fmadd_ps(_mm_mul_ps(_3,rxy),rxy,z2);
			    t2   = _mm_mul_ps(x,_mm_mul_ps(y,z));
			    H[12]= _mm_div_ps(_mm_mul_ps(t1,t2),t0);
			    H[10]= H[12];
			    t3   = _mm_mul_ps(x,_mm_sub_ps(_mm_add_ps(x2,y2),z2));
			    H[15]= _mm_div_ps(xmm4r4_negate(t3),c0);
			    H[11]= H[15];
			    t4   = _mm_mul_ps(y,_mm_sub_ps(_mm_add_ps(x2,y2),z2));
			    H[16]= _mm_div_ps(xmm4r4_negate(t4),c0);
			    H[14]= H[16];
			 }
			 else if(sysType==1) {
			    const __m128 _2 = _mm_set1_ps(2.0f);
			    const __m128 _3 = _mm_set1_ps(3.0f);
			    const __m128 _0 = _mm_setzero_ps();
                            __m128 x4,rxz,r4xz,rxz3,r4;
			    __m128 c0,c1,c2,c3,c4,c5;
			    x4   = _mm_mul_ps(x2,x2);
			    r4xz = _mm_add_ps(x2,z2);
			    rxz  = _mm_sqrt_ps(r4xz);
			    rxz3 = _mm_mul_ps(rxz,_mm_mul_ps(rxz,rxz));
			    r4   = _mm_add_ps(r4xz,y2);
			    r4xz = _mm_mul_ps(r4xz,r4xz);
			    r4   = _mm_mul_ps(r4,r4);
			    c0   = _mm_mul_ps(rxz3,r4);
			    c5   = _mm_mul_ps(rxz,r4);
			    H[0] = _mm_div_ps(xmm4r4_negate(_mm_mul_ps(_2,_mm_mul_ps(x,z))),r4xz);
			    H[4] = _0;
			    H[8] = _mm_div_ps(_mm_mul_ps(_2,_mm_mul_ps(x,z)),r4xz);
			    H[3] = _0;
			    H[1] = H[3];
			    H[6] = _mm_div_ps(_mm_sub_ps(x2,z2),r4xz);
			    H[2] = H[6];
			    H[7] = _0;
			    H[5] = H[7];
			    c1   = _mm_mul_ps(z2,_mm_add_ps(y2,z2));
			    c2   = _mm_fmadd_ps(_2,x4,_mm_mul_ps(x2,z2));
			    H[9] = _mm_div_ps(_mm_mul_ps(y,_mm_sub_ps(c2,c1)),c0);
			    H[13]= xmm4r4_negate(_mm_div_ps(_mm_mul_ps(_2,_mm_mul_ps(y,rxz))),r4);
			    c3   = _mm_mul_ps(x2,_mm_mul_ps(_mm_add_ps(y,z),_mm_sub_ps(y,z)));
			    c4   = _mm_fmsub_ps(x4,_2,_mm_mul_ps(z2,z2));
			    H[17]= xmm4r4_negate(_mm_div_ps(_mm_mul_ps(y,_mm_add_ps(c3,c4))),c0);
			    H[12]= xmm4r4_negate(_mm_div_ps(_mm_mul_ps(_mm_add_ps(_mm_sub_ps(x2,y2),z2)),c5));
			    H[10]= H[12];
			    c1   = _mm_fmadd_ps(_3,x2,_mm_fmadd_ps(_3,z2,y2));
			    c2   = _mm_mul_ps(x,_mm_mul_ps(y,z));
			    H[15]= _mm_div_ps(_mm_mul_ps(c1,c2),c0);
			    H[11]= H[15];
			    c3   = _mm_mul_ps(z,_mm_add_ps(_mm_sub_ps(x2,y2),z2));
			    H[16]= xmm4r4_negate(_mm_div_ps(c3,c5));
			    H[14]= H[16];
			 }
			 else if(sysType==3) {
			    const __m128 _2 = _mm_set1_ps(2.0f);
			    const __m128 _3 = _mm_set1_ps(3.0f);
			    const __m128 _0 = _mm_setzero_ps();
			    __m128 rxy,r4,r4xy,rxy3;
			    __m128 c0,c1,c2,c3,c4,c5;
			    r4xy  = _mm_add_ps(x2,y2);
			    rxy   = _mm_sqrt_ps(r4xy);
			    r4    = _mm_add_ps(r4xy,z2);
			    r4    = _mm_mul_ps(r4,r4);
			    r4xy  = _mm_mul_ps(r4xy,r4xy);
			    rxy3  = _mm_mul_ps(rxy,_mm_mul_ps(rxy,rxy));
			    c0    = _mm_mul_ps(rxy3,r4);
			    c5    = _mm_mul_ps(rxy,r4);
			    H[0]  = xmm4r4_negate(_mm_div_ps(_mm_mul_ps(_2,_mm_mul_ps(x,y)),r4xy));
			    H[4]  = H[0];
			    H[8]  = _0;
			    H[3]  = _mm_div_ps(_mm_mul_ps(_mm_sub_ps(x,y),
								_mm_add_ps(x,y)),r4xy);
			    H[1]  = H[3];
			    H[6]  = H[0];
			    H[2]  = H[6];
			    H[7]  = _0;
			    H[5]  = H[7];
			    c1    = _mm_mul_ps(y2,_mm_add_ps(y2,z2));
			    c2    = _mm_fmadd_ps(_2,x4,_mm_mul_ps(x2,y2));
			    H[9]  = _mm_div_ps(_mm_mul_ps(z,_mm_sub_ps(c2,c1)),c0);
			    c3    = _mm_mul_ps(x2,_mm_sub_ps(z2,y2));
			    c4    = _mm_fmsub_ps(_2,_mm_mul_ps(y2,y2),x4);
			    H[13] = xmm4r4_negate(_mm_div_ps(
							_mm_mul_ps(z,_mm_add_ps(c3,c4))),c0);
			    H[17] = xmm4r4_negate(_mm_div_ps(_mm_mul_ps(_2,
								_mm_mul_ps(rxy,z)),r4));
			    c1    = _mm_mul_ps(x,_mm_mul_ps(y,z));
			    c2    = _mm_fmadd_ps(_3,_mm_mul_ps(rxy,rxy),z2);
			    H[12] = _mm_div_ps(_mm_mul_ps(c1,c2),c0);
			    H[10] = H[12];
			    c3    = _mm_mul_ps(x,_mm_sub_ps(
			                                  _mm_add_ps(x2,y2),z2));
			    H[15] = xmm4r4_negate(_mm_div_ps(c3,c5));
			    H[11] = H[15];
			    c4    = _mm_mul_ps(y,_mm_sub_ps(
			                                  _mm_add_ps(x2,y2),z2));
			    H[16] = xmm4r4_negate(_mm_div_ps(c4,c5));
			    H[14] = H[16];
				
			 }
			 else { // sysType == 2
                            const __m128 _2 = _mm_set1_ps(2.0f);
			    const __m128 _3 = _mm_set1_ps(3.0f);
			    const __m128 _0 = _mm_setzero_ps();
			    __m128 r4xy,r4,rxy,x4,y4,rxy3;
			    __m128 c0,c1,c2,c3,c4,c5;
			    r4xy  = _mm_add_ps(x2,y2);
			    r4    = _mm_add_ps(r4xy,z2);
			    r4    = _mm_mul_ps(r4,r4);
			    rxy   = _mm_sqrt_ps(r4xy);
			    r4xy  = _mm_mul_ps(r4xy,r4xy);
			    rxy3  = _mm_mul_ps(rxy,_mm_mul_ps(rxy,rxy));
			    x4    = _mm_mul_ps(x2,x2);
			    y4    = _mm_mul_ps(y2,y2);
			    c0    = _mm_mul_ps(rxy3,r4);
			    c5    = _mm_mul_ps(rxy,r4);
			    H[0]  = _mm_div_ps(_mm_mul_ps(_2,
			                          _mm_mul_ps(x,y)),r4xy);
			    H[4]  = xmm4r4_negate(H[0]);
			    H[8]  = _0;
			    H[3]  = _mm_div_ps(_mm_sub_ps(y2,x2),r4xy);
			    H[1]  = H[0];
			    H[6]  = _0;
			    H[2]  = H[6];
			    H[7]  = _0;
			    H[5]  = H[7];
			    c1    = _mm_fmsub_ps(xmm4r4_negate(_2),x4,_mm_mul_ps(x2,y2));
			    c2    = _mm_fmadd_ps(y2,z2,y4);
			    H[9]  = _mm_div_ps(_mm_mul_ps(z,_mm_add_ps(c1,c2)),c0);
			    c3    = _mm_mul_ps(x2,_mm_sub_ps(z2,y2));
			    c4    = _mm_sub_ps(x4,_mm_mul_ps(_2,y4));
			    H[13] = _mm_div_ps(_mm_mul_ps(z,_mm_add_ps(c3,c4)),c0);
			    H[17] = _mm_div_ps(_mm_mul_ps(_2,_mm_mul_ps(rxy,z)),r4);
			    c1    = _mm_fmadd_ps(_3,_mm_mul_ps(rxy,rxy),z2);
			    c2    = _mm_mul_ps(x,_mm_mul_ps(y,z));
			    H[12] = xmm4r4_negate(_mm_div_ps(_mm_mul_ps(c1,c2),c0));
			    H[10] = H[12];
			    H[15] = _mm_div_ps(_mm_mul_ps(x,_mm_sub_ps(
			                                          _mm_add_ps(x2,y2),z2)),c5);
			    H[11] = H[15];
			    H[16] = _mm_div_ps(_mm_mul_ps(y,_mm_sub_ps(
			                                          _mm_add_ps(x2,y2),z2)),c5);
			    H[14] = H[16];
			}
		   }

/*A C++-only implementations of functions for computing the gradient of
*spherical azimuth and elevation angles. See the Matlab equivalent for
*more comments.
*
*February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
**May  2022 Bernard Gingold, manually vectorized.
**/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/		   

                     
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void spher_ang_grad_xmm2r8_a(double * __restrict __ATTR_ALIGN__(16) Mat0,
		                                   double * __restrict __ATTR_ALIGN__(16) Mat1,
						   double * __restrict __ATTR_ALIGN__(16) Mat2,
						   double * __restrict __ATTR_ALIGN__(16) Mat3,
						   double * __restrict __ATTR_ALIGN__(16) Mat4,
						   double * __restrict __ATTR_ALIGN__(16) Mat5,
						 const __m128d G0,
						 const __m128d G1,
						 const __m128d G2,
						 const __m128d Rx_x,
						 const __m128d Rx_y,
						 const __m128d Rx_z,
						 const __m128d * __restrict __ATTR_ALIGN__(16) M,
						 const int32_t sysType) {
						 
			  const __m128d _0 = _mm_setzero_pd();
                          __m128d J0,J1,J2,J3,J4,J;
			  const __m128d M0    = M[0];
			  const __m128d M1    = M[1];
			  const __m128d temp0 = _mm_sub_pd(G0,Rx_x);
			  const __m128d M3    = M[3];
			  const __m128d M4    = M[4];
			  const __m128d temp1 = _mm_sub_pd(G1,Rx_y);
			  const __m128d M6    = M[6];
			  const __m128d M7    = M[7];
			  const __m128d temp2 = _mm_sub_pd(G2,Rx_z);
			  const __m128d x     = _mm_fmadd_pd(M0,temp0,
			                              _mm_fmadd_pd(M3,temp1,
						            _mm_mul_pd(M6,temp2)));
			  const __m128d M2    = M[2];
			  const __m128d y     = _mm_fmadd_pd(M1,temp0,
			                              _mm_fmadd_pd(M4,temp1,
						            _mm_mul_pd(M7,temp2)));
			  const __m128d M5    = M[5];
			  const __m128d M8    = M[8];
			  const __m128d z     =  _mm_fmadd_pd(M2,temp0,
			                              _mm_fmadd_pd(M5,temp1,
						            _mm_mul_pd(M8,temp2)));
			  if(sysType==0) {
                             const __m128d r2   = _mm_fmadd_pd(x,x,
			                               _mm_fmadd_pd(y,y,
						                 _mm_mul_pd(z,z)));
			     const __m128d sqrv = _mm_fmadd_pd(x,x,
			                                     _mm_mul_pd(y,y));
			     const __m128d sqrtv= _mm_sqrt_pd(sqrv);
			     const __m128d denom= _mm_mul_pd(r2,sqrtv);
			      //Derivatives with respect to x.
			     J0 = xmm2r8_negate(_mm_div_pd(y,sqrv));
			     J1 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(x,z),denom));
			      //Derivatives with respect to y.
			     J2 = _mm_div_pd(x,sqrv);
			     J3 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(y,z),denom));
			     //Derivatives with respect to z.
			     J4 = _0;
			     J5 = _mm_div_pd(sqrtv,r2);
			  }
			  else if(sysType==2) {
                             const __m128d r2   = _mm_fmadd_pd(x,x,
			                               _mm_fmadd_pd(y,y,
						                 _mm_mul_pd(z,z)));
			     const __m128d sqrv = _mm_fmadd_pd(x,x,
			                                     _mm_mul_pd(y,y));
			     const __m128d sqrtv= _mm_sqrt_pd(sqrv);
			     const __m128d denom= _mm_mul_pd(r2,sqrtv);
			     J0 = xmm2r8_negate(_mm_div_pd(y,sqrv));
			     J1 = _mm_div_pd(_mm_mul_pd(x,z),denom);
			     J2 = _mm_div_pd(x,sqrv);
			     J3 = _mm_div_pd(_mm_mul_pd(y,z),denom);
			     J4 = _0;
			     J5 = xmm2r8_negate(_mm_div_pd(sqrtv,r2));
			  }
			  else if(sysType==3) {
                             const __m128d r2   = _mm_fmadd_pd(x,x,
			                               _mm_fmadd_pd(y,y,
						                 _mm_mul_pd(z,z)));
			     const __m128d sqrv = _mm_fmadd_pd(x,x,
			                                     _mm_mul_pd(y,y));
			     const __m128d sqrtv= _mm_sqrt_pd(sqrv);
			     J0 = _mm_div_pd(y,sqrv);
			     J1 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(x,z),denom));
			     J2 = xmm2r8_negate(_mm_div_pd(x,sqrv));
			     J3 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(y,z),denom));
			     J4 = _0;
			     J5 = _mm_div_pd(sqrtv,r2);
			  }
			  else { //sysType==1
                             const __m128d r2   = _mm_fmadd_pd(x,x,
			                               _mm_fmadd_pd(y,y,
						                 _mm_mul_pd(z,z)));
			     const __m128d sqrv = _mm_fmadd_pd(x,x,
			                                     _mm_mul_pd(y,y));
			     const __m128d sqrtv= _mm_sqrt_pd(sqrv);
			     J0 = _mm_div_pd(z,sqrv);
			     J1 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(x,y),denom));
			     J2 = _0;
			     J3 = _mm_div_pd(sqrtv,r2);
			     J4 = xmm2r8_negate(_mm_div_pd(x,sqrv));
			     J5 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(z,y),denom));
			  }

			  //Rotate from local back to global coordinates.
                          //J=J*M;
			  _mm_store_pd(&Mat0[0],_mm_fmadd_pd(J0,M0,
			                   _mm_fmadd_pd(J2,M1,
					              _mm_mul_pd(J4,M2))));
			  _mm_store_pd(&Mat1[0],_mm_fmadd_pd(J1,M0,
			                   _mm_fmadd_pd(J3,M1,
					              _mm_mul_pd(J5,M2))));
			  _mm_store_pd(&Mat2[0],_mm_fmadd_pd(J0,M3,
			                   _mm_fmadd_pd(J2,M4,
					              _mm_mul_pd(J4,M5))));
			  _mm_store_pd(&Mat3[0],_mm_fmadd_pd(J1,M3,
			                   _mm_fmadd_pd(J3,M4,
					              _mm_mul_pd(J5,M5))));
			  _mm_store_pd(&Mat4[0],_mm_fmadd_pd(J0,M6,
			                   _mm_fmadd_pd(J2,M7,
					              _mm_mul_pd(J4,M8))));
			  _mm_store_pd(&Mat5[0],_mm_fmadd_pd(J1,M6,
			                   _mm_fmadd_pd(J3,M7,
					              _mm_mul_pd(J5,M8))));
			  
		    }


		      
                      __ATTR_ALWAYS_INLINE__
		     static inline
		      void spher_ang_grad_xmm2r8_u(double * __restrict  Mat0,
		                                   double * __restrict  Mat1,
						   double * __restrict  Mat2,
						   double * __restrict  Mat3,
						   double * __restrict  Mat4,
						   double * __restrict  Mat5,
						   const __m128d G0,
						   const __m128d G1,
						   const __m128d G2,
						   const __m128d Rx_x,
						   const __m128d Rx_y,
						   const __m128d Rx_z,
						   const __m128d * __restrict __ATTR_ALIGN__(16) M,
						   const int32_t sysType) {
						 
			  const __m128d _0 = _mm_setzero_pd();
                          __m128d J0,J1,J2,J3,J4,J;
			  const __m128d M0    = M[0];
			  const __m128d M1    = M[1];
			  const __m128d temp0 = _mm_sub_pd(G0,Rx_x);
			  const __m128d M3    = M[3];
			  const __m128d M4    = M[4];
			  const __m128d temp1 = _mm_sub_pd(G1,Rx_y);
			  const __m128d M6    = M[6];
			  const __m128d M7    = M[7];
			  const __m128d temp2 = _mm_sub_pd(G2,Rx_z);
			  const __m128d x     = _mm_fmadd_pd(M0,temp0,
			                              _mm_fmadd_pd(M3,temp1,
						            _mm_mul_pd(M6,temp2)));
			  const __m128d M2    = M[2];
			  const __m128d y     = _mm_fmadd_pd(M1,temp0,
			                              _mm_fmadd_pd(M4,temp1,
						            _mm_mul_pd(M7,temp2)));
			  const __m128d M5    = M[5];
			  const __m128d M8    = M[8];
			  const __m128d z     =  _mm_fmadd_pd(M2,temp0,
			                              _mm_fmadd_pd(M5,temp1,
						            _mm_mul_pd(M8,temp2)));
			  if(sysType==0) {
                             const __m128d r2   = _mm_fmadd_pd(x,x,
			                               _mm_fmadd_pd(y,y,
						                 _mm_mul_pd(z,z)));
			     const __m128d sqrv = _mm_fmadd_pd(x,x,
			                                     _mm_mul_pd(y,y));
			     const __m128d sqrtv= _mm_sqrt_pd(sqrv);
			     const __m128d denom= _mm_mul_pd(r2,sqrtv);
			      //Derivatives with respect to x.
			     J0 = xmm2r8_negate(_mm_div_pd(y,sqrv));
			     J1 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(x,z),denom));
			      //Derivatives with respect to y.
			     J2 = _mm_div_pd(x,sqrv);
			     J3 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(y,z),denom));
			     //Derivatives with respect to z.
			     J4 = _0;
			     J5 = _mm_div_pd(sqrtv,r2);
			  }
			  else if(sysType==2) {
                             const __m128d r2   = _mm_fmadd_pd(x,x,
			                               _mm_fmadd_pd(y,y,
						                 _mm_mul_pd(z,z)));
			     const __m128d sqrv = _mm_fmadd_pd(x,x,
			                                     _mm_mul_pd(y,y));
			     const __m128d sqrtv= _mm_sqrt_pd(sqrv);
			     const __m128d denom= _mm_mul_pd(r2,sqrtv);
			     J0 = xmm2r8_negate(_mm_div_pd(y,sqrv));
			     J1 = _mm_div_pd(_mm_mul_pd(x,z),denom);
			     J2 = _mm_div_pd(x,sqrv);
			     J3 = _mm_div_pd(_mm_mul_pd(y,z),denom);
			     J4 = _0;
			     J5 = xmm2r8_negate(_mm_div_pd(sqrtv,r2));
			  }
			  else if(sysType==3) {
                             const __m128d r2   = _mm_fmadd_pd(x,x,
			                               _mm_fmadd_pd(y,y,
						                 _mm_mul_pd(z,z)));
			     const __m128d sqrv = _mm_fmadd_pd(x,x,
			                                     _mm_mul_pd(y,y));
			     const __m128d sqrtv= _mm_sqrt_pd(sqrv);
			     J0 = _mm_div_pd(y,sqrv);
			     J1 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(x,z),denom));
			     J2 = xmm2r8_negate(_mm_div_pd(x,sqrv));
			     J3 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(y,z),denom));
			     J4 = _0;
			     J5 = _mm_div_pd(sqrtv,r2);
			  }
			  else { //sysType==1
                             const __m128d r2   = _mm_fmadd_pd(x,x,
			                               _mm_fmadd_pd(y,y,
						                 _mm_mul_pd(z,z)));
			     const __m128d sqrv = _mm_fmadd_pd(x,x,
			                                     _mm_mul_pd(y,y));
			     const __m128d sqrtv= _mm_sqrt_pd(sqrv);
			     J0 = _mm_div_pd(z,sqrv);
			     J1 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(x,y),denom));
			     J2 = _0;
			     J3 = _mm_div_pd(sqrtv,r2);
			     J4 = xmm2r8_negate(_mm_div_pd(x,sqrv));
			     J5 = xmm2r8_negate(_mm_div_pd(
			                                  _mm_mul_pd(z,y),denom));
			  }

			  //Rotate from local back to global coordinates.
                          //J=J*M;
			  _mm_storeu_pd(&Mat0[0],_mm_fmadd_pd(J0,M0,
			                   _mm_fmadd_pd(J2,M1,
					              _mm_mul_pd(J4,M2))));
			  _mm_storeu_pd(&Mat1[0],_mm_fmadd_pd(J1,M0,
			                   _mm_fmadd_pd(J3,M1,
					              _mm_mul_pd(J5,M2))));
			  _mm_storeu_pd(&Mat2[0],_mm_fmadd_pd(J0,M3,
			                   _mm_fmadd_pd(J2,M4,
					              _mm_mul_pd(J4,M5))));
			  _mm_storeu_pd(&Mat3[0],_mm_fmadd_pd(J1,M3,
			                   _mm_fmadd_pd(J3,M4,
					              _mm_mul_pd(J5,M5))));
			  _mm_storeu_pd(&Mat4[0],_mm_fmadd_pd(J0,M6,
			                   _mm_fmadd_pd(J2,M7,
					              _mm_mul_pd(J4,M8))));
			  _mm_storeu_pd(&Mat5[0],_mm_fmadd_pd(J1,M6,
			                   _mm_fmadd_pd(J3,M7,
					              _mm_mul_pd(J5,M8))));
			  
		    }




		    
                      __ATTR_ALWAYS_INLINE__
		     static inline
		      void spher_ang_grad_xmm4r4_a(float * __restrict __ATTR_ALIGN__(16) Mat0,
		                                    float * __restrict __ATTR_ALIGN__(16) Mat1,
						    float * __restrict __ATTR_ALIGN__(16) Mat2,
						    float * __restrict __ATTR_ALIGN__(16) Mat3,
						    float * __restrict __ATTR_ALIGN__(16) Mat4,
						    float * __restrict __ATTR_ALIGN__(16) Mat5,
						    const __m128 G0,
						    const __m128 G1,
						    const __m128 G2,
						    const __m128 Rx_x,
						    const __m128 Rx_y,
						    const __m128 Rx_z,
						    const __m128 * __restrict __ATTR_ALIGN__(16) M,
						    const int32_t sysType) {
						 
			  const __m128 _0 = _mm_setzero_ps();
                          __m128 J0,J1,J2,J3,J4,J;
			  const __m128 M0    = M[0];
			  const __m128 M1    = M[1];
			  const __m128 temp0 = _mm_sub_ps(G0,Rx_x);
			  const __m128 M3    = M[3];
			  const __m128 M4    = M[4];
			  const __m128 temp1 = _mm_sub_ps(G1,Rx_y);
			  const __m128 M6    = M[6];
			  const __m128 M7    = M[7];
			  const __m128 temp2 = _mm_sub_ps(G2,Rx_z);
			  const __m128 x     = _mm_fmadd_ps(M0,temp0,
			                              _mm_fmadd_ps(M3,temp1,
						            _mm_mul_ps(M6,temp2)));
			  const __m128 M2    = M[2];
			  const __m128 y     = _mm_fmadd_ps(M1,temp0,
			                              _mm_fmadd_ps(M4,temp1,
						            _mm_mul_ps(M7,temp2)));
			  const __m128 M5    = M[5];
			  const __m128 M8    = M[8];
			  const __m128 z     =  _mm_fmadd_ps(M2,temp0,
			                              _mm_fmadd_ps(M5,temp1,
						            _mm_mul_ps(M8,temp2)));
			  if(sysType==0) {
                             const __m128 r2   = _mm_fmadd_ps(x,x,
			                               _mm_fmadd_ps(y,y,
						                 _mm_mul_ps(z,z)));
			     const __m128 sqrv = _mm_fmadd_ps(x,x,
			                                     _mm_mul_ps(y,y));
			     const __m128 sqrtv= _mm_sqrt_ps(sqrv);
			     const __m128 denom= _mm_mul_ps(r2,sqrtv);
			      //Derivatives with respect to x.
			     J0 = xmm4r4_negate(_mm_div_ps(y,sqrv));
			     J1 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(x,z),denom));
			      //Derivatives with respect to y.
			     J2 = _mm_div_ps(x,sqrv);
			     J3 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(y,z),denom));
			     //Derivatives with respect to z.
			     J4 = _0;
			     J5 = _mm_div_ps(sqrtv,r2);
			  }
			  else if(sysType==2) {
                             const __m128 r2   = _mm_fmadd_ps(x,x,
			                               _mm_fmadd_ps(y,y,
						                 _mm_mul_ps(z,z)));
			     const __m128 sqrv = _mm_fmadd_ps(x,x,
			                                     _mm_mul_ps(y,y));
			     const __m128 sqrtv= _mm_sqrt_ps(sqrv);
			     const __m128 denom= _mm_mul_ps(r2,sqrtv);
			     J0 = xmm4r4_negate(_mm_div_ps(y,sqrv));
			     J1 = _mm_div_ps(_mm_mul_ps(x,z),denom);
			     J2 = _mm_div_ps(x,sqrv);
			     J3 = _mm_div_ps(_mm_mul_ps(y,z),denom);
			     J4 = _0;
			     J5 = xmm4r4_negate(_mm_div_ps(sqrtv,r2));
			  }
			  else if(sysType==3) {
                             const __m128 r2   = _mm_fmadd_ps(x,x,
			                               _mm_fmadd_ps(y,y,
						                 _mm_mul_ps(z,z)));
			     const __m128 sqrv = _mm_fmadd_ps(x,x,
			                                     _mm_mul_ps(y,y));
			     const __m128 sqrtv= _mm_sqrt_ps(sqrv);
			     J0 = _mm_div_ps(y,sqrv);
			     J1 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(x,z),denom));
			     J2 = xmm4r4_negate(_mm_div_ps(x,sqrv));
			     J3 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(y,z),denom));
			     J4 = _0;
			     J5 = _mm_div_ps(sqrtv,r2);
			  }
			  else { //sysType==1
                             const __m128 r2   = _mm_fmadd_ps(x,x,
			                               _mm_fmadd_ps(y,y,
						                 _mm_mul_ps(z,z)));
			     const __m128 sqrv = _mm_fmadd_ps(x,x,
			                                     _mm_mul_ps(y,y));
			     const __m128 sqrtv= _mm_sqrt_ps(sqrv);
			     J0 = _mm_div_ps(z,sqrv);
			     J1 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(x,y),denom));
			     J2 = _0;
			     J3 = _mm_div_ps(sqrtv,r2);
			     J4 = xmm4r4_negate(_mm_div_ps(x,sqrv));
			     J5 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(z,y),denom));
			  }

			  //Rotate from local back to global coordinates.
                          //J=J*M;
			  _mm_store_ps(&Mat0[0],_mm_fmadd_ps(J0,M0,
			                   _mm_fmadd_ps(J2,M1,
					              _mm_mul_ps(J4,M2))));
			  _mm_store_ps(&Mat1[0],_mm_fmadd_ps(J1,M0,
			                   _mm_fmadd_ps(J3,M1,
					              _mm_mul_ps(J5,M2))));
			  _mm_store_ps(&Mat2[0],_mm_fmadd_ps(J0,M3,
			                   _mm_fmadd_ps(J2,M4,
					              _mm_mul_ps(J4,M5))));
			  _mm_store_ps(&Mat3[0],_mm_fmadd_ps(J1,M3,
			                   _mm_fmadd_ps(J3,M4,
					              _mm_mul_ps(J5,M5))));
			  _mm_store_ps(&Mat4[0],_mm_fmadd_ps(J0,M6,
			                   _mm_fmadd_ps(J2,M7,
					              _mm_mul_ps(J4,M8))));
			  _mm_store_ps(&Mat5[0],_mm_fmadd_ps(J1,M6,
			                   _mm_fmadd_ps(J3,M7,
					              _mm_mul_ps(J5,M8))));
			  
		    }


		   
                      __ATTR_ALWAYS_INLINE__
		      static inline
		      void spher_ang_grad_xmm4r4_u(float * __restrict  Mat0,
		                                    float * __restrict  Mat1,
						    float * __restrict  Mat2,
						    float * __restrict  Mat3,
						    float * __restrict  Mat4,
						    float * __restrict  Mat5,
						    const __m128 G0,
						    const __m128 G1,
						    const __m128 G2,
						    const __m128 Rx_x,
						    const __m128 Rx_y,
						    const __m128 Rx_z,
						    const __m128 * __restrict __ATTR_ALIGN__(16) M,
						    const int32_t sysType) {
						 
			  const __m128 _0 = _mm_setzero_ps();
                          __m128 J0,J1,J2,J3,J4,J;
			  const __m128 M0    = M[0];
			  const __m128 M1    = M[1];
			  const __m128 temp0 = _mm_sub_ps(G0,Rx_x);
			  const __m128 M3    = M[3];
			  const __m128 M4    = M[4];
			  const __m128 temp1 = _mm_sub_ps(G1,Rx_y);
			  const __m128 M6    = M[6];
			  const __m128 M7    = M[7];
			  const __m128 temp2 = _mm_sub_ps(G2,Rx_z);
			  const __m128 x     = _mm_fmadd_ps(M0,temp0,
			                              _mm_fmadd_ps(M3,temp1,
						            _mm_mul_ps(M6,temp2)));
			  const __m128 M2    = M[2];
			  const __m128 y     = _mm_fmadd_ps(M1,temp0,
			                              _mm_fmadd_ps(M4,temp1,
						            _mm_mul_ps(M7,temp2)));
			  const __m128 M5    = M[5];
			  const __m128 M8    = M[8];
			  const __m128 z     =  _mm_fmadd_ps(M2,temp0,
			                              _mm_fmadd_ps(M5,temp1,
						            _mm_mul_ps(M8,temp2)));
			  if(sysType==0) {
                             const __m128 r2   = _mm_fmadd_ps(x,x,
			                               _mm_fmadd_ps(y,y,
						                 _mm_mul_ps(z,z)));
			     const __m128 sqrv = _mm_fmadd_ps(x,x,
			                                     _mm_mul_ps(y,y));
			     const __m128 sqrtv= _mm_sqrt_ps(sqrv);
			     const __m128 denom= _mm_mul_ps(r2,sqrtv);
			      //Derivatives with respect to x.
			     J0 = xmm4r4_negate(_mm_div_ps(y,sqrv));
			     J1 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(x,z),denom));
			      //Derivatives with respect to y.
			     J2 = _mm_div_ps(x,sqrv);
			     J3 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(y,z),denom));
			     //Derivatives with respect to z.
			     J4 = _0;
			     J5 = _mm_div_ps(sqrtv,r2);
			  }
			  else if(sysType==2) {
                             const __m128 r2   = _mm_fmadd_ps(x,x,
			                               _mm_fmadd_ps(y,y,
						                 _mm_mul_ps(z,z)));
			     const __m128 sqrv = _mm_fmadd_ps(x,x,
			                                     _mm_mul_ps(y,y));
			     const __m128 sqrtv= _mm_sqrt_ps(sqrv);
			     const __m128 denom= _mm_mul_ps(r2,sqrtv);
			     J0 = xmm4r4_negate(_mm_div_ps(y,sqrv));
			     J1 = _mm_div_ps(_mm_mul_ps(x,z),denom);
			     J2 = _mm_div_ps(x,sqrv);
			     J3 = _mm_div_ps(_mm_mul_ps(y,z),denom);
			     J4 = _0;
			     J5 = xmm4r4_negate(_mm_div_ps(sqrtv,r2));
			  }
			  else if(sysType==3) {
                             const __m128 r2   = _mm_fmadd_ps(x,x,
			                               _mm_fmadd_ps(y,y,
						                 _mm_mul_ps(z,z)));
			     const __m128 sqrv = _mm_fmadd_ps(x,x,
			                                     _mm_mul_ps(y,y));
			     const __m128 sqrtv= _mm_sqrt_ps(sqrv);
			     J0 = _mm_div_ps(y,sqrv);
			     J1 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(x,z),denom));
			     J2 = xmm4r4_negate(_mm_div_ps(x,sqrv));
			     J3 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(y,z),denom));
			     J4 = _0;
			     J5 = _mm_div_ps(sqrtv,r2);
			  }
			  else { //sysType==1
                             const __m128 r2   = _mm_fmadd_ps(x,x,
			                               _mm_fmadd_ps(y,y,
						                 _mm_mul_ps(z,z)));
			     const __m128 sqrv = _mm_fmadd_ps(x,x,
			                                     _mm_mul_ps(y,y));
			     const __m128 sqrtv= _mm_sqrt_ps(sqrv);
			     J0 = _mm_div_ps(z,sqrv);
			     J1 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(x,y),denom));
			     J2 = _0;
			     J3 = _mm_div_ps(sqrtv,r2);
			     J4 = xmm4r4_negate(_mm_div_ps(x,sqrv));
			     J5 = xmm4r4_negate(_mm_div_ps(
			                                  _mm_mul_ps(z,y),denom));
			  }

			  //Rotate from local back to global coordinates.
                          //J=J*M;
			  _mm_storeu_ps(&Mat0[0],_mm_fmadd_ps(J0,M0,
			                   _mm_fmadd_ps(J2,M1,
					              _mm_mul_ps(J4,M2))));
			  _mm_storeu_ps(&Mat1[0],_mm_fmadd_ps(J1,M0,
			                   _mm_fmadd_ps(J3,M1,
					              _mm_mul_ps(J5,M2))));
			  _mm_storeu_ps(&Mat2[0],_mm_fmadd_ps(J0,M3,
			                   _mm_fmadd_ps(J2,M4,
					              _mm_mul_ps(J4,M5))));
			  _mm_storeu_ps(&Mat3[0],_mm_fmadd_ps(J1,M3,
			                   _mm_fmadd_ps(J3,M4,
					              _mm_mul_ps(J5,M5))));
			  _mm_storeu_ps(&Mat4[0],_mm_fmadd_ps(J0,M6,
			                   _mm_fmadd_ps(J2,M7,
					              _mm_mul_ps(J4,M8))));
			  _mm_storeu_ps(&Mat5[0],_mm_fmadd_ps(J1,M6,
			                   _mm_fmadd_ps(J3,M7,
					              _mm_mul_ps(J5,M8))));
			  
		    }


	
    } // math

} // gms

















#endif /*__GMS_SPHER_GRAD_SSE_HPP__*/
