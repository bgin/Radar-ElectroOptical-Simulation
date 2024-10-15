

#ifndef __GMS_SPHER_GRAD_AVX512_H__
#define __GMS_SPHER_GRAD_AVX512_H__ 170520220848


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
                    
                      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void  spher_inv_jac_zmm8r8_a(double * __restrict __ATTR_ALIGN__(64) J0,
		                                   double * __restrict __ATTR_ALIGN__(64) J1,
						   double * __restrict __ATTR_ALIGN__(64) J2,
						   double * __restrict __ATTR_ALIGN__(64) J3,
						   double * __restrict __ATTR_ALIGN__(64) J4,
						   double * __restrict __ATTR_ALIGN__(64) J5,
						   double * __restrict __ATTR_ALIGN__(64) J6,
						   double * __restrict __ATTR_ALIGN__(64) J7,
						   double * __restrict __ATTR_ALIGN__(64) J8,
						   const __m512d x,
						   const __m512d y,
						   const __m512d z,
						   const int32_t sysType); 

		      
                      __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void  spher_inv_jac_zmm8r8_u(double * __restrict J0,
		                                   double * __restrict J1,
						   double * __restrict J2,
						   double * __restrict J3,
						   double * __restrict J4,
						   double * __restrict J5,
						   double * __restrict J6,
						   double * __restrict J7,
						   double * __restrict J8,
						   const __m512d x,
						   const __m512d y,
						   const __m512d z,
						   const int32_t sysType); 



	             


		    
                       __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void  spher_inv_jac_zmm16r4_a(float * __restrict __ATTR_ALIGN__(64) J0,
		                                    float * __restrict __ATTR_ALIGN__(64) J1,
						    float * __restrict __ATTR_ALIGN__(64) J2,
						    float * __restrict __ATTR_ALIGN__(64) J3,
						    float * __restrict __ATTR_ALIGN__(64) J4,
						    float * __restrict __ATTR_ALIGN__(64) J5,
						    float * __restrict __ATTR_ALIGN__(64) J6,
						    float * __restrict __ATTR_ALIGN__(64) J7,
						    float * __restrict __ATTR_ALIGN__(64) J8,
						    const __m512 x,
						    const __m512 y,
						    const __m512 z,
						    const int32_t sysType); 

		    
                       __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void  spher_inv_jac_zmm16r4_u(float * __restrict J0,
		                                    float * __restrict J1,
						    float * __restrict J2,
						    float * __restrict J3,
						    float * __restrict J4,
						    float * __restrict J5,
						    float * __restrict J6,
						    float * __restrict J7,
						    float * __restrict J8,
						    const __m512 x,
						    const __m512 y,
						    const __m512 z,
						    const int32_t sysType); 

		   /*
                        Different definitions with an array of type __m512d/__m512 i.e. AoSoA
                        instead of SoA.
                    */

		      
                       __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void  spher_inv_jac_zmm8r8(__m512d * __restrict __ATTR_ALIGN__(64) J, //flatten 9x1 array
		                                 const __m512d x,
						 const __m512d y,
						 const __m512d z,
						 const int32_t sysType); 

		      
                       __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void  spher_inv_jac_zmm16r4(__m512 * __restrict __ATTR_ALIGN__(64) J,
		                                  const __m512 x,
						  const __m512 y,
						  const __m512 z,
						  const int32_t sysType); 

/**SPHERANGHESSIANCPP A C++-only implementation of a function for
 *          computing the Hessian of azimuth and elevation in spherical
 *          coordinates.  See the Matlab equivalent for more comments.
 *
 *June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *May  2022 Bernard Gingold, manually vectorized.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
	              

                      
                       __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void spher_ang_hess_zmm8r8(__m512d * __restrict __ATTR_ALIGN__(64) H,
		                                 const __m512d G0,
						 const __m512d G1,
						 const __m512d G2,
						 const int32_t sysType); 

		    
                       __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void spher_ang_hess_zmm16r4(__m512 * __restrict __ATTR_ALIGN__(64) H,
		                                  const __m512 G0,
						  const __m512 G1,
						  const __m512 G2,
						  const int32_t sysType); 
/*A C++-only implementations of functions for computing the gradient of
*spherical azimuth and elevation angles. See the Matlab equivalent for
*more comments.
*
*February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
**May  2022 Bernard Gingold, manually vectorized.
**/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/		   

                     
                       __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void spher_ang_grad_zmm8r8_a(double * __restrict __ATTR_ALIGN__(64) Mat0,
		                                   double * __restrict __ATTR_ALIGN__(64) Mat1,
						   double * __restrict __ATTR_ALIGN__(64) Mat2,
						   double * __restrict __ATTR_ALIGN__(64) Mat3,
						   double * __restrict __ATTR_ALIGN__(64) Mat4,
						   double * __restrict __ATTR_ALIGN__(64) Mat5,
						 const __m512d G0,
						 const __m512d G1,
						 const __m512d G2,
						 const __m512d Rx_x,
						 const __m512d Rx_y,
						 const __m512d Rx_z,
						 const __m512d * __restrict __ATTR_ALIGN__(64) M,
						 const int32_t sysType); 

		       
                       __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void spher_ang_grad_zmm8r8_u(double * __restrictMat0,
		                                   double * __restrict  Mat1,
						   double * __restrict  Mat2,
						   double * __restrict  Mat3,
						   double * __restrict  Mat4,
						   double * __restrict  Mat5,
						   const __m512d G0,
						   const __m512d G1,
						   const __m512d G2,
						   const __m512d Rx_x,
						   const __m512d Rx_y,
						   const __m512d Rx_z,
						   const __m512d * __restrict __ATTR_ALIGN__(64) M,
						   const int32_t sysType); 



		     
                       __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void spher_ang_grad_zmm16r4_a(float * __restrict __ATTR_ALIGN__(64) Mat0,
		                                    float * __restrict __ATTR_ALIGN__(64) Mat1,
						    float * __restrict __ATTR_ALIGN__(64) Mat2,
						    float * __restrict __ATTR_ALIGN__(64) Mat3,
						    float * __restrict __ATTR_ALIGN__(64) Mat4,
						    float * __restrict __ATTR_ALIGN__(64) Mat5,
						    const __m512 G0,
						    const __m512 G1,
						    const __m512 G2,
						    const __m512 Rx_x,
						    const __m512 Rx_y,
						    const __m512 Rx_z,
						    const __m512 * __restrict __ATTR_ALIGN__(64) M,
						    const int32_t sysType); 

		      
                       __ATTR_HOT__
                      __ATTR_ALIGN__(32)
		      void spher_ang_grad_zmm16r4_u(float * __restrict  Mat0,
		                                    float * __restrict  Mat1,
						    float * __restrict  Mat2,
						    float * __restrict  Mat3,
						    float * __restrict  Mat4,
						    float * __restrict  Mat5,
						    const __m512 G0,
						    const __m512 G1,
						    const __m512 G2,
						    const __m512 Rx_x,
						    const __m512 Rx_y,
						    const __m512 Rx_z,
						    const __m512 * __restrict __ATTR_ALIGN__(64) M,
						    const int32_t sysType); 

	
    } // math

} // gms

















#endif /*__GMS_SPHER_GRAD_AVX512_H__*/
