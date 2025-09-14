

#ifndef __GMS_SPHER_GRAD_AVX_H__
#define __GMS_SPHER_GRAD_AVX_H__ 170520220848


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

namespace file_info 
{

 const unsigned int GMS_SPHER_GRAD_AVX_MAJOR = 1U;
 const unsigned int GMS_SPHER_GRAD_AVX_MINOR = 0U;
 const unsigned int GMS_SPHER_GRAD_AVX_MICRO = 0U;
 const unsigned int GMS_SPHER_GRAD_AVX_FULLVER =
  1000U*GMS_SPHER_GRAD_AVX_MAJOR+100U*GMS_SPHER_GRAD_AVX_MINOR+10U*GMS_SPHER_GRAD_AVX_MICRO;
 const char * const GMS_SPHER_GRAD_AVX_CREATION_DATE = "17-05-2022 08:48 +00200 (TUE 17 MAY 2022 08:48 GMT+2)";
 const char * const GMS_SPHER_GRAD_AVX_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const GMS_SPHER_GRAD_AVX_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const GMS_SPHER_GRAD_AVX_SYNOPSIS      = "AVX based spherical coordinates hessians and jacobians functions (vectorized).";


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
                    
                     __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void  spher_inv_jac_ymm4r8_a(double * __restrict __ATTR_ALIGN__(32) J0,
		                                   double * __restrict __ATTR_ALIGN__(32) J1,
						   double * __restrict __ATTR_ALIGN__(32) J2,
						   double * __restrict __ATTR_ALIGN__(32) J3,
						   double * __restrict __ATTR_ALIGN__(32) J4,
						   double * __restrict __ATTR_ALIGN__(32) J5,
						   double * __restrict __ATTR_ALIGN__(32) J6,
						   double * __restrict __ATTR_ALIGN__(32) J7,
						   double * __restrict __ATTR_ALIGN__(32) J8,
						   const __m256d x,
						   const __m256d y,
						   const __m256d z,
						   const int32_t sysType); 

		    
                      __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void  spher_inv_jac_ymm4r8_u(double * __restrict J0,
		                                   double * __restrict J1,
						   double * __restrict J2,
						   double * __restrict J3,
						   double * __restrict J4,
						   double * __restrict J5,
						   double * __restrict J6,
						   double * __restrict J7,
						   double * __restrict J8,
						   const __m256d x,
						   const __m256d y,
						   const __m256d z,
						   const int32_t sysType); 

           


		       __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void  spher_inv_jac_ymm8r4_a(float * __restrict __ATTR_ALIGN__(32) J0,
		                                    float * __restrict __ATTR_ALIGN__(32) J1,
						    float * __restrict __ATTR_ALIGN__(32) J2,
						    float * __restrict __ATTR_ALIGN__(32) J3,
						    float * __restrict __ATTR_ALIGN__(32) J4,
						    float * __restrict __ATTR_ALIGN__(32) J5,
						    float * __restrict __ATTR_ALIGN__(32) J6,
						    float * __restrict __ATTR_ALIGN__(32) J7,
						    float * __restrict __ATTR_ALIGN__(32) J8,
						    const __m256 x,
						    const __m256 y,
						    const __m256 z,
						    const int32_t sysType);

		    
                       __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void  spher_inv_jac_ymm8r4_u(float * __restrict J0,
		                                    float * __restrict J1,
						    float * __restrict J2,
						    float * __restrict J3,
						    float * __restrict J4,
						    float * __restrict J5,
						    float * __restrict J6,
						    float * __restrict J7,
						    float * __restrict J8,
						    const __m256 x,
						    const __m256 y,
						    const __m256 z,
						    const int32_t sysType); 

		   /*
                        Different definitions with an array of type __m256d/__m256 i.e. AoSoA
                        instead of SoA.
                    */

#if 0		     
                       __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void  spher_inv_jac_ymm4r8(__m256d * __restrict __ATTR_ALIGN__(32) J, //flatten 9x1 array
		                                 const __m256d x,
						 const __m256d y,
						 const __m256d z,
						 const int32_t sysType); 

		     
                       __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void  spher_inv_jac_ymm8r4(__m256 * __restrict __ATTR_ALIGN__(32) J,
		                                  const __m256 x,
						  const __m256 y,
						  const __m256 z,
						  const int32_t sysType); 
#endif 
/**SPHERANGHESSIANCPP A C++-only implementation of a function for
 *          computing the Hessian of azimuth and elevation in spherical
 *          coordinates.  See the Matlab equivalent for more comments.
 *
 *June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *May  2022 Bernard Gingold, manually vectorized.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
	              

                   
                      __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void spher_ang_hess_ymm4r8(__m256d * __restrict __ATTR_ALIGN__(32) H,
		                                 const __m256d G0,
						                 const __m256d G1,
						                 const __m256d G2,
						                 const int32_t sysType); 

		    
                       __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void spher_ang_hess_ymm8r4(__m256 * __restrict __ATTR_ALIGN__(32) H,
		                                const __m256 G0,
						                const __m256 G1,
						                const __m256 G2,
						                const int32_t sysType);
/*A C++-only implementations of functions for computing the gradient of
*spherical azimuth and elevation angles. See the Matlab equivalent for
*more comments.
*
*February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
**May  2022 Bernard Gingold, manually vectorized.
**/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/		   

#if 0                     
                       __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void spher_ang_grad_ymm4r8_a(double * __restrict __ATTR_ALIGN__(32) Mat0,
		                                   double * __restrict __ATTR_ALIGN__(32) Mat1,
						   double * __restrict __ATTR_ALIGN__(32) Mat2,
						   double * __restrict __ATTR_ALIGN__(32) Mat3,
						   double * __restrict __ATTR_ALIGN__(32) Mat4,
						   double * __restrict __ATTR_ALIGN__(32) Mat5,
						 const __m256d G0,
						 const __m256d G1,
						 const __m256d G2,
						 const __m256d Rx_x,
						 const __m256d Rx_y,
						 const __m256d Rx_z,
						 const __m256d * __restrict __ATTR_ALIGN__(32) M,
						 const int32_t sysType); 

		      
                      __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void spher_ang_grad_ymm4r8_u(double * __restrict  Mat0,
		                                   double * __restrict  Mat1,
						   double * __restrict  Mat2,
						   double * __restrict  Mat3,
						   double * __restrict  Mat4,
						   double * __restrict  Mat5,
						   const __m256d G0,
						   const __m256d G1,
						   const __m256d G2,
						   const __m256d Rx_x,
						   const __m256d Rx_y,
						   const __m256d Rx_z,
						   const __m256d * __restrict __ATTR_ALIGN__(32) M,
						   const int32_t sysType);


		       __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void spher_ang_grad_ymm8r4_a(float * __restrict __ATTR_ALIGN__(32) Mat0,
		                                    float * __restrict __ATTR_ALIGN__(32) Mat1,
						    float * __restrict __ATTR_ALIGN__(32) Mat2,
						    float * __restrict __ATTR_ALIGN__(32) Mat3,
						    float * __restrict __ATTR_ALIGN__(32) Mat4,
						    float * __restrict __ATTR_ALIGN__(32) Mat5,
						    const __m256 G0,
						    const __m256 G1,
						    const __m256 G2,
						    const __m256 Rx_x,
						    const __m256 Rx_y,
						    const __m256 Rx_z,
						    const __m256 * __restrict __ATTR_ALIGN__(32) M,
						    const int32_t sysType);


		   
                       __ATTR_ALIGN__(32)
                      __ATTR_HOT__
		      void spher_ang_grad_ymm8r4_u(float * __restrict  Mat0,
		                                    float * __restrict  Mat1,
						    float * __restrict  Mat2,
						    float * __restrict  Mat3,
						    float * __restrict  Mat4,
						    float * __restrict  Mat5,
						    const __m256 G0,
						    const __m256 G1,
						    const __m256 G2,
						    const __m256 Rx_x,
						    const __m256 Rx_y,
						    const __m256 Rx_z,
						    const __m256 * __restrict __ATTR_ALIGN__(32) M,
						    const int32_t sysType); 
#endif  

	
    } // math

} // gms

















#endif /*__GMS_SPHER_GRAD_AVX_H__*/
