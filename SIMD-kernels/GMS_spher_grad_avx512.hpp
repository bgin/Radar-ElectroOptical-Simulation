

#ifndef __GMS_SPHER_GRAD_AVX512_HPP__
#define __GMS_SPHER_GRAD_AVX512_HPP__ 170520220848


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

 const unsigned int gGMS_SPHER_GRAD_AVX512_MAJOR = 1U;
 const unsigned int gGMS_SPHER_GRAD_AVX512_MINOR = 0U;
 const unsigned int gGMS_SPHER_GRAD_AVX512_MICRO = 0U;
 const unsigned int gGMS_SPHER_GRAD_AVX512_FULLVER =
  1000U*gGMS_SPHER_GRAD_AVX512_MAJOR+100U*gGMS_SPHER_GRAD_AVX512_MINOR+10U*gGMS_SPHER_GRAD_AVX512_MICRO;
 const char * const pgGMS_SPHER_GRAD_AVX512_CREATION_DATE = "17-05-2022 08:48 +00200 (TUE 17 MAY 2022 08:48 GMT+2)";
 const char * const pgGMS_SPHER_GRAD_AVX512_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const pgGMS_SPHER_GRAD_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const pgGMS_SPHER_GRAD_AVX512_SYNOPSIS      = "AVX512 based spherical coordinates hessians and jacobians functions (vectorized)."


}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#if (USE_SLEEF_LIB) == 1
#include "GMS_sleefsimddp.hpp"
#include "GMS_sleefsimdsp.hpp"
#endif
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
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void  spher_inv_jac_zmm8r8(__m512d &J0,
		                                 __m512d &J1,
						 __m512d &J2,
						 __m512d &J3,
						 __m512d &J4,
						 __m512d &J5,
						 __m512d &J6,
						 __m512d &J7,
						 __m512d &J8,
						 const __m512d x,
						 const __m512d y,
						 const __m512d z,
						 const int32_t sysType) {
#if (USE_SLEEF_LIB) == 1
                         const __m512d sinaz = xsin(y);
			 const __m512d cosaz = xcos(y);
			 const __m512d sinel = xsin(z);
			 const __m512d cosel = xcos(z);
#else
                         const __m512d sinaz = _mm512_sin_pd(y);
			 const __m512d cosaz = _mm512_cos_pd(y);
			 const __m512d sinel = _mm512_sin_pd(z);
			 const __m512d cosel = _mm512_cos_pd(z);
#endif
                         if(sysType==0) {
			    const __m512d nx = zmm8r8_negate(x):
                            //dx/dr
			    J0 = _mm512_mul_pd(cosaz,cosel);
			    //dy/dr
			    J1 = _mm512_mul_pd(cosel,sinaz);
			    //dz/dr
			    J2 = sinel;
			    //dx/dAz
			    J3 = _mm512_mul_pd(nx,
			                       _mm512_mul_pd(cosel,sinaz));
			    //dy/dAz
			    J4 = _mm512_mul_pd(x,_mm512_mul_pd(cosaz,cosel));
			    //dz/dAz
			    J5 = _mm512_setzero_pd();
			    //dx/dEl
			    J6 = _mm512_mul_pd(nx,
			                       _mm512_mul_pd(cosaz,sinel));
			    //dy/dEl
			    J7 = _mm512_mul_pd(nx,
			                       _mm512_mul_pd(sinaz,sinel));
			    //dz/dEl
			    J8 = _mm512_mul_pd(x,cosel);
			 }
			 else if(sysType==1) {
                             const __m512d nx    = zmm8r8_negate(x):
			     const __m512d cazel = _mm512_mul_pd(cosaz,cosel);
			     J0 = _mm512_mul_pd(cosel,sinaz);
			     J1 = sinel;
			     J2 = cazel;
			     J3 = _mm512_mul_pd(x,cazel);
			     J4 = _mm512_setzero_pd();
			     J5 = _mm512_mul_pd(nx,_mm512_mul_pd(cosel,sinaz));
			     J6 = _mm512_mul_pd(nx,_mm512_mul_pd(sinaz,sinel));
			     J7 = _mm512_mul_pd(x,cosel);
			     J8 = _mm512_mul_pd(nx,_mm512_mul_pd(cosaz,sinel));
			 }
			 else if(sysType==3) {
                             const __m512d nx    = zmm8r8_negate(x):
			     const __m512d sc    = _mm512_mul_pd(sinaz,cosel);
			     const __m512d cc    = _mm512_mul_pd(cosaz,cosel);
			     J0 = sc;
			     J1 = cc;
			     J2 = sinel;
			     J3 = _mm512_mul_pd(x,cc);
			     J4 = _mm512_mul_pd(nx,sc);
			     J5 = _mm512_setzero_pd();
			     J6 = _mm512_mul_pd(nx,_mm512_mul_pd(sinaz,sinel));
			     J7 = _mm512_mul_pd(nx,_mm512_mul_pd(cosaz,sinel));
			     J8 = _mm512_mul_pd(x,cosel);
			 }
			 else { // sysType==2
                             const __m512d nx    = zmm8r8_negate(x):
			     const __m512d ss    = _mm512_mul_pd(sinaz,sinel);
			     const __m512d cs    = _mm512_mul_pd(cosaz,sinel);
			     J0 = cs;
			     J1 = ss;
			     J2 = cosel;
			     J3 = _mm512_mul_pd(nx,ss);
			     J4 = _mm512_mul_pd(x,cs);
			     J5 = _mm512_setzero_pd();
			     J6 = _mm512_mul_pd(x,_mm512_mul_pd(cosaz,cosel));
			     J7 = _mm512_mul_pd(x,_mm512_mul_pd(cosel,sinaz));
			     J8 = _mm512_mul_pd(nx,sinel);
			 }
		   }


	             


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void  spher_inv_jac_zmm16r4(__m512 &J0,
		                                  __m512 &J1,
						  __m512 &J2,
						  __m512 &J3,
						  __m512 &J4,
						  __m512 &J5,
						  __m512 &J6,
						  __m512 &J7,
						  __m512 &J8,
						  const __m512 x,
						  const __m512 y,
						  const __m512 z,
						  const int32_t sysType) {
#if (USE_SLEEF_LIB) == 1
                         const __m512 sinaz = xsinf(y);
			 const __m512 cosaz = xcosf(y);
			 const __m512 sinel = xsinf(z);
			 const __m512 cosel = xcosf(z);
#else
                         const __m512 sinaz = _mm512_sin_ps(y);
			 const __m512 cosaz = _mm512_cos_ps(y);
			 const __m512 sinel = _mm512_sin_ps(z);
			 const __m512 cosel = _mm512_cos_ps(z);
#endif
                         if(sysType==0) {
			    const __m512 nx = zmm16r4_negate(x):
                            //dx/dr
			    J0 = _mm512_mul_ps(cosaz,cosel);
			    //dy/dr
			    J1 = _mm512_mul_ps(cosel,sinaz);
			    //dz/dr
			    J2 = sinel;
			    //dx/dAz
			    J3 = _mm512_mul_ps(nx,
			                       _mm512_mul_ps(cosel,sinaz));
			    //dy/dAz
			    J4 = _mm512_mul_ps(x,_mm512_mul_ps(cosaz,cosel));
			    //dz/dAz
			    J5 = _mm512_setzero_ps();
			    //dx/dEl
			    J6 = _mm512_mul_ps(nx,
			                       _mm512_mul_ps(cosaz,sinel));
			    //dy/dEl
			    J7 = _mm512_mul_ps(nx,
			                       _mm512_mul_ps(sinaz,sinel));
			    //dz/dEl
			    J8 = _mm512_mul_ps(x,cosel);
			 }
			 else if(sysType==1) {
                             const __m512 nx    = zmm16r4_negate(x):
			     const __m512 cazel = _mm512_mul_ps(cosaz,cosel);
			     J0 = _mm512_mul_ps(cosel,sinaz);
			     J1 = sinel;
			     J2 = cazel;
			     J3 = _mm512_mul_ps(x,cazel);
			     J4 = _mm512_setzero_ps();
			     J5 = _mm512_mul_ps(nx,_mm512_mul_ps(cosel,sinaz));
			     J6 = _mm512_mul_ps(nx,_mm512_mul_ps(sinaz,sinel));
			     J7 = _mm512_mul_ps(x,cosel);
			     J8 = _mm512_mul_ps(nx,_mm512_mul_ps(cosaz,sinel));
			 }
			 else if(sysType==3) {
                             const __m512 nx    = zmm16r4_negate(x):
			     const __m512 sc    = _mm512_mul_ps(sinaz,cosel);
			     const __m512 cc    = _mm512_mul_ps(cosaz,cosel);
			     J0 = sc;
			     J1 = cc;
			     J2 = sinel;
			     J3 = _mm512_mul_ps(x,cc);
			     J4 = _mm512_mul_ps(nx,sc);
			     J5 = _mm512_setzero_ps();
			     J6 = _mm512_mul_ps(nx,_mm512_mul_ps(sinaz,sinel));
			     J7 = _mm512_mul_ps(nx,_mm512_mul_ps(cosaz,sinel));
			     J8 = _mm512_mul_ps(x,cosel);
			 }
			 else { // sysType==2
                             const __m512 nx    = zmm16r4_negate(x):
			     const __m512 ss    = _mm512_mul_ps(sinaz,sinel);
			     const __m512 cs    = _mm512_mul_ps(cosaz,sinel);
			     J0 = cs;
			     J1 = ss;
			     J2 = cosel;
			     J3 = _mm512_mul_ps(nx,ss);
			     J4 = _mm512_mul_ps(x,cs);
			     J5 = _mm512_setzero_ps();
			     J6 = _mm512_mul_ps(x,_mm512_mul_ps(cosaz,cosel));
			     J7 = _mm512_mul_ps(x,_mm512_mul_ps(cosel,sinaz));
			     J8 = _mm512_mul_ps(nx,sinel);
			 }
		   }

		   /*
                        Different definitions with an array of type __m512d/__m512 i.e. AoSoA
                        instead of SoA.
                    */

		       __ATTR_REGCALL__
                       __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void  spher_inv_jac_zmm8r8(__m512d * __restrict __ATTR_ALIGN__(64) J, //flatten 9x1 array
		                                 const __m512d x,
						 const __m512d y,
						 const __m512d z,
						 const int32_t sysType) {
#if (USE_SLEEF_LIB) == 1
                         const __m512d sinaz = xsin(y);
			 const __m512d cosaz = xcos(y);
			 const __m512d sinel = xsin(z);
			 const __m512d cosel = xcos(z);
#else
                         const __m512d sinaz = _mm512_sin_pd(y);
			 const __m512d cosaz = _mm512_cos_pd(y);
			 const __m512d sinel = _mm512_sin_pd(z);
			 const __m512d cosel = _mm512_cos_pd(z);
#endif
                         if(sysType==0) {
			    const __m512d nx = zmm8r8_negate(x):
                            //dx/dr
			    J[0] = _mm512_mul_pd(cosaz,cosel);
			    //dy/dr
			    J[1] = _mm512_mul_pd(cosel,sinaz);
			    //dz/dr
			    J[2] = sinel;
			    //dx/dAz
			    J[3] = _mm512_mul_pd(nx,
			                       _mm512_mul_pd(cosel,sinaz));
			    //dy/dAz
			    J[4] = _mm512_mul_pd(x,_mm512_mul_pd(cosaz,cosel));
			    //dz/dAz
			    J[5] = _mm512_setzero_pd();
			    //dx/dEl
			    J[6] = _mm512_mul_pd(nx,
			                       _mm512_mul_pd(cosaz,sinel));
			    //dy/dEl
			    J[7] = _mm512_mul_pd(nx,
			                       _mm512_mul_pd(sinaz,sinel));
			    //dz/dEl
			    J[8] = _mm512_mul_pd(x,cosel);
			 }
			 else if(sysType==1) {
                             const __m512d nx    = zmm8r8_negate(x):
			     const __m512d cazel = _mm512_mul_pd(cosaz,cosel);
			     J[0] = _mm512_mul_pd(cosel,sinaz);
			     J[1] = sinel;
			     J[2] = cazel;
			     J[3] = _mm512_mul_pd(x,cazel);
			     J[4] = _mm512_setzero_pd();
			     J[5] = _mm512_mul_pd(nx,_mm512_mul_pd(cosel,sinaz));
			     J[6] = _mm512_mul_pd(nx,_mm512_mul_pd(sinaz,sinel));
			     J[7] = _mm512_mul_pd(x,cosel);
			     J[8] = _mm512_mul_pd(nx,_mm512_mul_pd(cosaz,sinel));
			 }
			 else if(sysType==3) {
                             const __m512d nx    = zmm8r8_negate(x):
			     const __m512d sc    = _mm512_mul_pd(sinaz,cosel);
			     const __m512d cc    = _mm512_mul_pd(cosaz,cosel);
			     J[0] = sc;
			     J[1] = cc;
			     J[2] = sinel;
			     J[3] = _mm512_mul_pd(x,cc);
			     J[4] = _mm512_mul_pd(nx,sc);
			     J[5] = _mm512_setzero_pd();
			     J[6] = _mm512_mul_pd(nx,_mm512_mul_pd(sinaz,sinel));
			     J[7] = _mm512_mul_pd(nx,_mm512_mul_pd(cosaz,sinel));
			     J[8] = _mm512_mul_pd(x,cosel);
			 }
			 else { // sysType==2
                             const __m512d nx    = zmm8r8_negate(x):
			     const __m512d ss    = _mm512_mul_pd(sinaz,sinel);
			     const __m512d cs    = _mm512_mul_pd(cosaz,sinel);
			     J[0] = cs;
			     J[1] = ss;
			     J[2] = cosel;
			     J[3] = _mm512_mul_pd(nx,ss);
			     J[4] = _mm512_mul_pd(x,cs);
			     J[5] = _mm512_setzero_pd();
			     J[6] = _mm512_mul_pd(x,_mm512_mul_pd(cosaz,cosel));
			     J[7] = _mm512_mul_pd(x,_mm512_mul_pd(cosel,sinaz));
			     J[8] = _mm512_mul_pd(nx,sinel);
			 }
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void  spher_inv_jac_zmm16r4(__m512 * __restrict __ATTR_ALIGN__(64) J,
		                                  const __m512 x,
						  const __m512 y,
						  const __m512 z,
						  const int32_t sysType) {
#if (USE_SLEEF_LIB) == 1
                         const __m512 sinaz = xsinf(y);
			 const __m512 cosaz = xcosf(y);
			 const __m512 sinel = xsinf(z);
			 const __m512 cosel = xcosf(z);
#else
                         const __m512 sinaz = _mm512_sin_ps(y);
			 const __m512 cosaz = _mm512_cos_ps(y);
			 const __m512 sinel = _mm512_sin_ps(z);
			 const __m512 cosel = _mm512_cos_ps(z);
#endif
                         if(sysType==0) {
			    const __m512 nx = zmm16r4_negate(x):
                            //dx/dr
			    J[0] = _mm512_mul_ps(cosaz,cosel);
			    //dy/dr
			    J[1] = _mm512_mul_ps(cosel,sinaz);
			    //dz/dr
			    J[2] = sinel;
			    //dx/dAz
			    J[3] = _mm512_mul_ps(nx,
			                       _mm512_mul_ps(cosel,sinaz));
			    //dy/dAz
			    J[4] = _mm512_mul_ps(x,_mm512_mul_ps(cosaz,cosel));
			    //dz/dAz
			    J[5] = _mm512_setzero_ps();
			    //dx/dEl
			    J[6] = _mm512_mul_ps(nx,
			                       _mm512_mul_ps(cosaz,sinel));
			    //dy/dEl
			    J[7] = _mm512_mul_ps(nx,
			                       _mm512_mul_ps(sinaz,sinel));
			    //dz/dEl
			    J[8] = _mm512_mul_ps(x,cosel);
			 }
			 else if(sysType==1) {
                             const __m512 nx    = zmm16r4_negate(x):
			     const __m512 cazel = _mm512_mul_ps(cosaz,cosel);
			     J[0] = _mm512_mul_ps(cosel,sinaz);
			     J[1] = sinel;
			     J[2] = cazel;
			     J[3] = _mm512_mul_ps(x,cazel);
			     J[4] = _mm512_setzero_ps();
			     J[5] = _mm512_mul_ps(nx,_mm512_mul_ps(cosel,sinaz));
			     J[6] = _mm512_mul_ps(nx,_mm512_mul_ps(sinaz,sinel));
			     J[7] = _mm512_mul_ps(x,cosel);
			     J[8] = _mm512_mul_ps(nx,_mm512_mul_ps(cosaz,sinel));
			 }
			 else if(sysType==3) {
                             const __m512 nx    = zmm16r4_negate(x):
			     const __m512 sc    = _mm512_mul_ps(sinaz,cosel);
			     const __m512 cc    = _mm512_mul_ps(cosaz,cosel);
			     J[0] = sc;
			     J[1] = cc;
			     J[2] = sinel;
			     J[3] = _mm512_mul_ps(x,cc);
			     J[4] = _mm512_mul_ps(nx,sc);
			     J[5] = _mm512_setzero_ps();
			     J[6] = _mm512_mul_ps(nx,_mm512_mul_ps(sinaz,sinel));
			     J[7] = _mm512_mul_ps(nx,_mm512_mul_ps(cosaz,sinel));
			     J[8] = _mm512_mul_ps(x,cosel);
			 }
			 else { // sysType==2
                             const __m512 nx    = zmm16r4_negate(x):
			     const __m512 ss    = _mm512_mul_ps(sinaz,sinel);
			     const __m512 cs    = _mm512_mul_ps(cosaz,sinel);
			     J[0] = cs;
			     J[1] = ss;
			     J[2] = cosel;
			     J[3] = _mm512_mul_ps(nx,ss);
			     J[4] = _mm512_mul_ps(x,cs);
			     J[5] = _mm512_setzero_ps();
			     J[6] = _mm512_mul_ps(x,_mm512_mul_ps(cosaz,cosel));
			     J[7] = _mm512_mul_ps(x,_mm512_mul_ps(cosel,sinaz));
			     J[8] = _mm512_mul_ps(nx,sinel);
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
	              

                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      void spher_ang_hess_zmm8r8(__m512d * __restrict __ATTR_ALIGN__(64) H,
		                                 const __m512d G0,
						 const __m512d G1,
						 const __m512d G2,
						 const int32_t sysType) {
                         const __m512d x  = G0;
			 const __m512d x2 = _mm512_mul_pd(x,x);
			 const __m512d y  = G1;
			 const __m512d y2 = _mm512_mul_pd(y,y);
			 const __m512d z  = G2;
			 const __m512d z2 = _mm512_mul_pd(z,z);
			 const __m512d x4 = _mm512_mul_pd(x2,x2);
			 
			 if(sysType==0) {
			    const __m512d _2 = _mm512_set1_pd(2.0);
			    const __m512d _3 = _mm512_set1_pd(3.0);
			    const __m512d _0 = _mm512_setzero_pd();
                            __m512d rxy,r4,r4xy,rxy3;
			    __m512d t0,t1,t2,t3,t4,t5,c0;
			    r4xy = _mm512_add_pd(x2,y2);
			    rxy  = _mm512_sqrt_pd(r4xy);
			    r4   = _mm512_add_pd(r4xy,z2);
			    r4   = _mm512_mul_pd(r4,r4);
			    r4xy = _mm512_mul_pd(r4xy,r4xy);
			    rxy3 = _mm512_mul_pd(rxy,_mm512_mul_pd(rxy,rxy));
			    t0   = _mm512_mul_pd(rxy3,r4);
			    c0   = _mm512_mul_pd(rxy,r4);
			    H[0] = _mm512_div_pd(_mm512_mul_pd(_2,_mm512_mul_pd(x,y)),r4xy);
			    H[4] = zmm8r8_negate(H[0]);
			    H[8] = _0;
			    H[3] = _mm512_div_pd(_mm512_sub_pd(y2,x2),r4xy);
			    H[1] = H[3];
			    H[6] = _0;
			    H[2] = H[6];
			    H[7] = _0;
			    H[5] = H[7];
			    t1   = _mm512_mul_pd(y2,_mm512_add_pd(y2,z2));
			    t2   = _mm512_fmadd_pd(x2,y2,_mm512_mul_pd(_2,x4));
			    H[9] = _mm512_div_pd(_mm512_mul_pd(z,_mm512_sub_pd(t1,t2)),t0);
			    t3   = _mm512_mul_pd(x2,_mm512_sub_pd(z2,y2));
			    t4   = _mm512_sub_pd(x4,_mm512_mul_pd(_2,_mm512_mul_pd(y2,y2)));
			    H[13]= _mm512_div_pd(zmm8r8_negate(_mm512_mul_pd(z,_mm512_add_pd(t3,t4))),t0);
			    H[17]= _mm512_div_pd(zmm8r8_negate(_mm512_mul_pd(_2,_mm512_mul_pd(rxy,z))),r4);
			    t1   = _mm512_fmadd_pd(_mm512_mul_pd(_3,rxy),rxy,z2);
			    t2   = _mm512_mul_pd(x,_mm512_mul_pd(y,z));
			    H[12]= _mm512_div_pd(_mm512_mul_pd(t1,t2),t0);
			    H[10]= H[12];
			    t3   = _mm512_mul_pd(x,_mm512_sub_pd(_mm512_add_pd(x2,y2),z2));
			    H[15]= _mm512_div_pd(zmm8r8_negate(t3),c0);
			    H[11]= H[15];
			    t4   = _mm512_mul_pd(y,_mm512_sub_pd(_mm512_add_pd(x2,y2),z2));
			    H[16]= _mm512_div_pd(zmm8r8_negate(t4),c0);
			    H[14]= H[16];
			 }
			 else if(sysType==1) {
			    const __m512d _2 = _mm512_set1_pd(2.0);
			    const __m512d _3 = _mm512_set1_pd(3.0);
			    const __m512d _0 = _mm512_setzero_pd();
                            __m512d x4,rxz,r4xz,rxz3,r4;
			    __m512d c0,c1,c2,c3,c4,c5;
			    x4   = _mm512_mul_pd(x2,x2);
			    r4xz = _mm512_add_pd(x2,z2);
			    rxz  = _mm512_sqrt_pd(r4xz);
			    rxz3 = _mm512_mul_pd(rxz,_mm512_mul_pd(rxz,rxz));
			    r4   = _mm512_add_pd(r4xz,y2);
			    r4xz = _mm512_mul_pd(r4xz,r4xz);
			    r4   = _mm512_mul_pd(r4,r4);
			    c0   = _mm512_mul_pd(rxz3,r4);
			    c5   = _mm512_mul_pd(rxz,r4);
			    H[0] = _mm512_div_pd(zmm8r8_negate(_mm512_mul_pd(_2,_mm512_mul_pd(x,z))),r4xz);
			    H[4] = _0;
			    H[8] = _mm512_div_pd(_mm512_mul_pd(_2,_mm512_mul_pd(x,z)),r4xz);
			    H[3] = _0;
			    H[1] = H[3];
			    H[6] = _mm512_div_pd(_mm512_sub_pd(x2,z2),r4xz);
			    H[2] = H[6];
			    H[7] = _0;
			    H[5] = H[7];
			    c1   = _mm512_mul_pd(z2,_mm512_add_pd(y2,z2));
			    c2   = _mm512_fmadd_pd(_2,x4,_mm512_mul_pd(x2,z2));
			    H[9] = _mm512_div_pd(_mm512_mul_pd(y,_mm512_sub_pd(c2,c1)),c0);
			    H[13]= zmm8r8_negate(_mm512_div_pd(_mm512_mul_pd(_2,_mm512_mul_pd(y,rxz))),r4);
			    c3   = _mm512_mul_pd(x2,_mm512_mul_pd(_mm512_add_pd(y,z),_mm512_sub_pd(y,z)));
			    c4   = _mm512_fmsub_pd(x4,_2,_mm512_mul_pd(z2,z2));
			    H[17]= zmm8r8_negate(_mm512_div_pd(_mm512_mul_pd(y,_mm512_add_pd(c3,c4))),c0);
			    H[12]= zmm8r8_negate(_mm512_div_pd(_mm512_mul_pd(_mm512_add_pd(_mm512_sub_pd(x2,y2),z2)),c5));
			    H[10]= H[12];
			    c1   = _mm512_fmadd_pd(_3,x2,_mm512_fmadd_pd(_3,z2,y2));
			    c2   = _mm512_mul_pd(x,_mm512_mul_pd(y,z));
			    H[15]= _mm512_div_pd(_mm512_mul_pd(c1,c2),c0);
			    H[11]= H[15];
			    c3   = _mm512_mul_pd(z,_mm512_add_pd(_mm512_sub_pd(x2,y2),z2));
			    H[16]= zmm8r8_negate(_mm512_div_pd(c3,c5));
			    H[14]= H[16];
			 }
			 else if(sysType==3) {
			    const __m512d _2 = _mm512_set1_pd(2.0);
			    const __m512d _3 = _mm512_set1_pd(3.0);
			    const __m512d _0 = _mm512_setzero_pd();
			    __m512d rxy,r4,r4xy,rxy3;
			    __m512d c0,c1,c2,c3,c4,c5;
			    r4xy  = _mm512_add_pd(x2,y2);
			    rxy   = _mm512_sqrt_pd(r4xy);
			    r4    = _mm512_add_pd(r4xy,z2);
			    r4    = _mm512_mul_pd(r4,r4);
			    r4xy  = _mm512_mul_pd(r4xy,r4xy);
			    rxy3  = _mm512_mul_pd(rxy,_mm512_mul_pd(rxy,rxy));
			    c0    = _mm512_mul_pd(rxy3,r4);
			    c5    = _mm512_mul_pd(rxy,r4);
			    H[0]  = zmm8r8_negate(_mm512_div_pd(_mm512_mul_pd(_2,_mm512_mul_pd(x,y)),r4xy));
			    H[4]  = H[0];
			    H[8]  = _0;
			    H[3]  = _mm512_div_pd(_mm512_mul_pd(_mm512_sub_pd(x,y),
								_mm512_add_pd(x,y)),r4xy);
			    H[1]  = H[3];
			    H[6]  = H[0];
			    H[2]  = H[6];
			    H[7]  = _0;
			    H[5]  = H[7];
			    c1    = _mm512_mul_pd(y2,_mm512_add_pd(y2,z2));
			    c2    = _mm512_fmadd_pd(_2,x4,_mm512_mul_pd(x2,y2));
			    H[9]  = _mm512_div_pd(_mm512_mul_pd(z,_mm512_sub_pd(c2,c1)),c0);
			    c3    = _mm512_mul_pd(x2,_mm512_sub_pd(z2,y2));
			    c4    = _mm512_fmsub_pd(_2,_mm512_mul_pd(y2,y2),x4);
			    H[13] = zmm8r8_negate(_mm512_div_pd(
							_mm512_mul_pd(z,_mm512_add_pd(c3,c4))),c0);
			    H[17] = zmm8r8_negate(_mm512_div_pd(_mm512_mul_pd(_2,
								_mm512_mul_pd(rxy,z)),r4));
			    c1    = _mm512_mul_pd(x,_mm512_mul_pd(y,z));
			    c2    = _mm512_fmadd_pd(_3,_mm512_mul_pd(rxy,rxy),z2);
			    H[12] = _mm512_div_pd(_mm512_mul_pd(c1,c2),c0);
			    H[10] = H[12];
			    c3    = _mm512_mul_pd(x,_mm512_sub_pd(
			                                  _mm512_add_pd(x2,y2),z2));
			    H[15] = zmm8r8_negate(_mm512_div_pd(c3,c5));
			    H[11] = H[15];
			    c4    = _mm512_mul_pd(y,_mm512_sub_pd(
			                                  _mm512_add_pd(x2,y2),z2));
			    H[16] = zmm8r8_negate(_mm512_div_pd(c4,c5));
			    H[14] = H[16];
				
			 }
			 else { // sysType == 2
                            const __m512d _2 = _mm512_set1_pd(2.0);
			    const __m512d _3 = _mm512_set1_pd(3.0);
			    const __m512d _0 = _mm512_setzero_pd();
			    __m512d r4xy,r4,rxy,x4,y4,rxy3;
			    __m512d c0,c1,c2,c3,c4,c5;
			    r4xy  = _mm512_add_pd(x2,y2);
			    r4    = _mm512_add_pd(r4xy,z2);
			    r4    = _mm512_mul_pd(r4,r4);
			    rxy   = _mm512_sqrt_pd(r4xy);
			    r4xy  = _mm512_mul_pd(r4xy,r4xy);
			    rxy3  = _mm512_mul_pd(rxy,_mm512_mul_pd(rxy,rxy));
			    x4    = _mm512_mul_pd(x2,x2);
			    y4    = _mm512_mul_pd(y2,y2);
			    c0    = _mm512_mul_pd(rxy3,r4);
			    c5    = _mm512_mul_pd(rxy,r4);
			    H[0]  = _mm512_div_pd(_mm512_mul_pd(_2,
			                          _mm512_mul_pd(x,y)),r4xy);
			    H[4]  = zmm8r8_negate(H[0]);
			    H[8]  = _0;
			    H[3]  = _mm512_div_pd(_mm512_sub_pd(y2,x2),r4xy);
			    H[1]  = H[0];
			    H[6]  = _0;
			    H[2]  = H[6];
			    H[7]  = _0;
			    H[5]  = H[7];
			    c1    = _mm512_fmsub_pd(zmm8r8_negate(_2),x4,_mm512_mul_pd(x2,y2));
			    c2    = _mm512_fmadd_pd(y2,z2,y4);
			    H[9]  = _mm512_div_pd(_mm512_mul_pd(z,_mm512_add_pd(c1,c2)),c0);
			    c3    = _mm512_mul_pd(x2,_mm512_sub_pd(z2,y2));
			    c4    = _mm512_sub_pd(x4,_mm512_mul_pd(_2,y4));
			    H[13] = _mm512_div_pd(_mm512_mul_pd(z,_mm512_add_pd(c3,c4)),c0);
			    H[17] = _mm512_div_pd(_mm512_mul_pd(_2,_mm512_mul_pd(rxy,z)),r4);
			    c1    = _mm512_fmadd_pd(_3,_mm512_mul_pd(rxy,rxy),z2);
			    c2    = _mm512_mul_pd(x,_mm512_mul_pd(y,z));
			    H[12] = zmm8r8_negate(_mm512_div_pd(_mm512_mul_pd(c1,c2),c0));
			    H[10] = H[12];
			    H[15] = _mm512_div_pd(_mm512_mul_pd(x,_mm512_sub_pd(
			                                          _mm512_add_pd(x2,y2),z2)),c5);
			    H[11] = H[15];
			    H[16] = _mm512_div_pd(_mm512_mul_pd(y,_mm512_sub_pd(
			                                          _mm512_add_pd(x2,y2),z2)),c5);
			    H[14] = H[16];
			}
		   }


		   
	
    } // math

} // gms

















#endif /*__GMS_SPHER_GRAD_AVX512_HPP__*/
