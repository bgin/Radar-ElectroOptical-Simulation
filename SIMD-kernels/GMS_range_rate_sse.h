
#ifndef __GMS_RANGE_RATE_SSE_H__
#define __GMS_RANGE_RATE_SSE_H__ 160820220629

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
@@Modified by Bernard Gingold, on 14-05-2022 11:11 +00200 (SAT 14 MAY 2022 11:11 GMT+2)
  contact: beniekg@gmail.com
*/

namespace file_info {

 const unsigned int gGMS_RANGE_RATE_SSE_MAJOR = 1U;
 const unsigned int gGMS_RANGE_RATE_SSE_MINOR = 0U;
 const unsigned int gGMS_RANGE_RATE_SSE_MICRO = 0U;
 const unsigned int gGMS_RANGE_RATE_SSE_FULLVER =
  1000U*gGMS_RANGE_RATE_SSE_MAJOR+100U*gGMS_RANGE_RATE_SSE_MINOR+10U*gGMS_RANGE_RATE_SSE_MICRO;
 const char * const pgGMS_RANGE_RATE_SSE_CREATION_DATE = "16-08-2024 06:29 +00200 (FRI 16 AUG 2024 06:29 GMT+2)";
 const char * const pgGMS_RANGE_RATE_SSE_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const pgGMS_RANGE_RATE_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const pgGMS_RANGE_RATE_SSE_SYNOPSIS      = "SSE based range-rate functions (vectorized)."


}


#include <immintrin.h>
#include "GMS_config.h"


namespace gms {

         namespace math {

/*GETRANGERATE2DGENCPP A C++ function to convert a Cartesian state in 2D
 *        into a non-relativistic range rate, ignoring atmospheric effects.
 *
 *INPUTS: xTar The 4X1 Cartesian position and velocity vectors
 *             [x;y;xDot;yDot].
 *  useHalfRange A boolean value specifying whether the bistatic (round-
 *             trip) range value (and hence the range rate) has been
 *             divided by two. 
 *         xTx The 4X1 [x;y;xDot;yDot] position and velocity vector of
 *             the transmitter in global Cartesian coordinates.
 *         xRx The 4X1 [x;y;xDot;yDot] position and velocity vector of
 *             the receiver in global Cartesian coordinates.
 *
 *OUTPUTS: rr The range rate as a double.
 *
 *See the comments to the Matlab function getRangeRate for more information
 *on how this function works.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *@@Modified by Bernard Gingold, on May 2022
 **/

                      
                      __ATTR_HOT__
                      __ATTR_VECTORCALL__
                      __m128d range_rate_2d_xmm2r8(const __m128d tar_x,
		                                   const __m128d tar_y,
						   const __m128d tar_xD,
						   const __m128d tar_yD,
						   const __m128d tx_x,
						   const __m128d tx_y,
						   const __m128d tx_xD,
						   const __m128d tx_yD,
						   const __m128d rx_x,
						   const __m128d rx_y,
						   const __m128d rx_xD,
						   const __m128d rx_yD,
						   const bool useHalfRange); 


		   
		        __ATTR_HOT__
                      __ATTR_VECTORCALL__
                      __m128 range_rate_2d_xmm4r4(const __m128 tar_x,
		                                   const __m128 tar_y,
						   const __m128 tar_xD,
						   const __m128 tar_yD,
						   const __m128 tx_x,
						   const __m128 tx_y,
						   const __m128 tx_xD,
						   const __m128 tx_yD,
						   const __m128 rx_x,
						   const __m128 rx_y,
						   const __m128 rx_xD,
						   const __m128 rx_yD,
						   const bool useHalfRange); 
/*GETRANGERATE3DGENCPP A C++ function to convert a Cartesian state in 3D
 *        into a non-relativistic range rate, ignoring atmospheric effects.
 *
 *INPUTS: xTar The 6X1 Cartesian position and velocity vectors
 *             [x;y;z;xDot;yDot;zDot].
 * useHalfRange A boolean value specifying whether the bistatic (round-
 *             trip) range value (and hence the range rate) has been
 *             divided by two. 
 *         xTx The 6X1 [x;y;z;xDot;yDot;zDot] position and velocity
 *             vector of the transmitter in global Cartesian coordinates.
 *         xTx The 6X1 [x;y;z;xDot;yDot;zDot] position and velocity
 *             vector of the receiver in global Cartesian coordinates.
 *
 *OUTPUTS: rr The range rate as a double.
 *
 *See the comments to the Matlab function getRangeRate for more information
 *on how this function works.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *@@Modified by Bernard Gingold, on May 2022
 **/
		   
				   
                    
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
                      __m128d range_rate_3d_xmm2r8(const __m128d tar_x,
		                                   const __m128d tar_y,
						   const __m128d tar_z,
						   const __m128d tar_xD,
						   const __m128d tar_yD,
						   const __m128d tar_zD,
						   const __m128d tx_x,
						   const __m128d tx_y,
						   const __m128d tx_z,
						   const __m128d tx_xD,
						   const __m128d tx_yD,
						   const __m128d tx_zD,
						   const __m128d rx_x,
						   const __m128d rx_y,
						   const __m128d rx_z,
						   const __m128d rx_xD,
						   const __m128d rx_yD,
						   const __m128d rx_zD,
						   const bool useHalfRange); 


		     
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
                      __m128 range_rate_3d_xmm4r4(const __m128 tar_x,
		                                   const __m128 tar_y,
						   const __m128 tar_z,
						   const __m128 tar_xD,
						   const __m128 tar_yD,
						   const __m128 tar_zD,
						   const __m128 tx_x,
						   const __m128 tx_y,
						   const __m128 tx_z,
						   const __m128 tx_xD,
						   const __m128 tx_yD,
						   const __m128 tx_zD,
						   const __m128 rx_x,
						   const __m128 rx_y,
						   const __m128 rx_z,
						   const __m128 rx_xD,
						   const __m128 rx_yD,
						   const __m128 rx_zD,
						   const bool useHalfRange);

/*A C++-only implementations of a functions for computing the gradient of
 *bistatic range.  See the Matlab equivalents for more comments.
 *
*February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
**@@Modified by Bernard Gingold, on May 2022
**/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
		     
                       __ATTR_HOT__
                      __ATTR_VECTORCALL__
                      __m128d range_grad_xmm2r8(const __m128d p,
					        const __m128d Tx,
					        const __m128d Rx,
					        const bool useHalfRange); 
					        
	                  __ATTR_HOT__
                      __ATTR_VECTORCALL__	        
                      __m128 range_grad_xmm4r4(const __m128 p,
					        const __m128 Tx,
					        const __m128 Rx,
					        const bool useHalfRange); 
/**RANGEHESSIANCPP A C++-only implementations of a function for computing
 *the Hessian of bistatic range.  See the Matlab equivalent for more
 *comments.
 *
 *June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *@@Modified by Bernard Gingold, on May 2022
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
		    
                      
                     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline  
                      __m128d range_hessian_1d_xmm2r8() {
                         __m128d H = _mm_setzero_pd();
			 return (H);
		     }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      
		      static inline  
                      __m128 range_hessian_1d_xmm4r4() {
                         __m128 H = _mm_setzero_ps();
			 return (H);
		     }



		     
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hessian_2d_xmm2r8_a(double * __restrict __ATTR_ALIGN__(16) H_0,
		                                     double * __restrict __ATTR_ALIGN__(16) H_1,
						     double * __restrict __ATTR_ALIGN__(16) H_2,
						     double * __restrict __ATTR_ALIGN__(16) H_3,
						     const __m128d x_0,
						     const __m128d x_1,
						     const bool useHalfRange); 


		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hessian_2d_xmm2r8_u(double * __restrict H_0,
		                                     double * __restrict H_1,
						     double * __restrict H_2,
						     double * __restrict H_3,
						     const __m128d x_0,
						     const __m128d x_1,
						     const bool useHalfRange); 


		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hessian_2d_xmm4r4_a(float * __restrict __ATTR_ALIGN__(16) H_0,
		                                      float * __restrict __ATTR_ALIGN__(16) H_1,
						      float * __restrict __ATTR_ALIGN__(16) H_2,
						      float * __restrict __ATTR_ALIGN__(16) H_3,
						      const __m128 x_0,
						      const __m128 x_1,
						      const bool useHalfRange); 


		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hessian_2d_xmm4r4_u(float * __restrict  H_0,
		                                      float * __restrict  H_1,
						      float * __restrict  H_2,
						      float * __restrict  H_3,
						      const __m128 x_0,
						      const __m128 x_1,
						      const bool useHalfRange);
						      
						      
		         __ATTR_HOT__
                      __ATTR_VECTORCALL__			      
		      void range_hessian_3d_xmm2r8_a(double * __restrict __ATTR_ALIGN__(16) H_0,
		                                     double * __restrict __ATTR_ALIGN__(16) H_1,
						     double * __restrict __ATTR_ALIGN__(16) H_2,
						     double * __restrict __ATTR_ALIGN__(16) H_3,
						     double * __restrict __ATTR_ALIGN__(16) H_4,
						     double * __restrict __ATTR_ALIGN__(16) H_5,
						     double * __restrict __ATTR_ALIGN__(16) H_6,
						     double * __restrict __ATTR_ALIGN__(16) H_7,
						     double * __restrict __ATTR_ALIGN__(16) H_8,
						     const __m128d x_0,
						     const __m128d x_1,
						     const __m128d x_2,
						     const bool useHalfRange); 
						     


		     
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hessian_3d_xmm2r8_u(double * __restrict H_0,
		                                     double * __restrict H_1,
						     double * __restrict H_2,
						     double * __restrict H_3,
						     double * __restrict H_4,
						     double * __restrict H_5,
						     double * __restrict H_6,
						     double * __restrict H_7,
						     double * __restrict H_8,
						   const __m128d x_0,
						   const __m128d x_1,
						   const __m128d x_2,
						   const bool useHalfRange); 
						   
				
		            __ATTR_HOT__
                      __ATTR_VECTORCALL__		   
		      void range_hessian_3d_xmm2r8(__m128d * __restrict __ATTR_ALIGN__(16) H,
		                                   const __m128d x_0,
						   const __m128d x_1,
						   const __m128d x_2,
						   const bool useHalfRange); 



		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hessian_3d_xmm4r4_a(float * __restrict __ATTR_ALIGN__(16) H_0,
		                                      float * __restrict __ATTR_ALIGN__(16) H_1,
						      float * __restrict __ATTR_ALIGN__(16) H_2,
						      float * __restrict __ATTR_ALIGN__(16) H_3,
						      float * __restrict __ATTR_ALIGN__(16) H_4,
						      float * __restrict __ATTR_ALIGN__(16) H_5,
						      float * __restrict __ATTR_ALIGN__(16) H_6,
						      float * __restrict __ATTR_ALIGN__(16) H_7,
						      float * __restrict __ATTR_ALIGN__(16) H_8,
						      const __m128 x_0,
						      const __m128 x_1,
						      const __m128 x_2,
						      const bool useHalfRange); 


		     
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hessian_3d_xmm4r4_u(float * __restrict  H_0,
		                                      float * __restrict  H_1,
						      float * __restrict  H_2,
						      float * __restrict  H_3,
						      float * __restrict  H_4,
						      float * __restrict  H_5,
						      float * __restrict  H_6,
						      float * __restrict  H_7,
						      float * __restrict  H_8,
						      const __m128 x_0,
						      const __m128 x_1,
						      const __m128 x_2,
						      const bool useHalfRange); 



		     
                      __ATTR_ALWAYS_INLINE__
		     
		      
		      static inline
		      __m128d range_hess_gen_1d_xmm2r8() {

		             return (_mm_setzero_pd());
		      }


		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hess_gen_2d_xmm2r8(__m128d &H_0,
		                                    __m128d &H_1,
						    __m128d &H_2,
						    __m128d &H_3,
						    const __m128d x_0,
						    const __m128d x_1,
						    const __m128d rx_0,
						    const __m128d rx_1,
						    const __m128d tx_0,
						    const __m128d tx_1,
						    const bool useHalfRange);



		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hess_gen_2d_xmm2r8_a(double * __restrict __ATTR_ALIGN__(16) H_0,
		                                      double * __restrict __ATTR_ALIGN__(16) H_1,
						      double * __restrict __ATTR_ALIGN__(16) H_2,
						      double * __restrict __ATTR_ALIGN__(16) H_3,
						      const __m128d x_0,
						      const __m128d x_1,
						      const __m128d rx_0,
						      const __m128d rx_1,
						      const __m128d tx_0,
						      const __m128d tx_1,
						      const bool useHalfRange); 


		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hess_gen_2d_xmm2r8_u(double * __restrict H_0,
		                                      double * __restrict H_1,
						      double * __restrict H_2,
						      double * __restrict H_3,
						      const __m128d x_0,
						      const __m128d x_1,
						      const __m128d rx_0,
						      const __m128d rx_1,
						      const __m128d tx_0,
						      const __m128d tx_1,
						      const bool useHalfRange); 
						      

	     
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hess_gen_2d_xmm4r4(__m128 &H_0,
		                                     __m128 &H_1,
						     __m128 &H_2,
						     __m128 &H_3,
						     const __m128 x_0,
						     const __m128 x_1,
						     const __m128 rx_0,
						     const __m128 rx_1,
						     const __m128 tx_0,
						     const __m128 tx_1,
						     const bool useHalfRange); 


		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hess_gen_2d_xmm4r4_a(float * __restrict __ATTR_ALIGN__(16) H_0,
		                                       float * __restrict __ATTR_ALIGN__(16) H_1,
						       float * __restrict __ATTR_ALIGN__(16) H_2,
						       float * __restrict __ATTR_ALIGN__(16) H_3,
						       const __m128 x_0,
						       const __m128 x_1,
						       const __m128 rx_0,
						       const __m128 rx_1,
						       const __m128 tx_0,
						       const __m128 tx_1,
						       const bool useHalfRange); 

		     
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hess_gen_2d_xmm4r4_u(float * __restrict  H_0,
		                                       float * __restrict  H_1,
						       float * __restrict  H_2,
						       float * __restrict  H_3,
						       const __m128 x_0,
						       const __m128 x_1,
						       const __m128 rx_0,
						       const __m128 rx_1,
						       const __m128 tx_0,
						       const __m128 tx_1,
						       const bool useHalfRange); 


		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hess_3d_xmm2r8(__m128d &H_0,
		                                __m128d &H_1,
						__m128d &H_2,
						__m128d &H_3,
						__m128d &H_4,
						__m128d &H_5,
						__m128d &H_6,
						__m128d &H_7,
						__m128d &H_8,
						const __m128d x_0,
						const __m128d x_1,
						const __m128d x_2,
						const __m128d rx_0,
						const __m128d rx_1,
						const __m128d rx_2,
						const __m128d tx_0,
						const __m128d tx_1,
						const __m128d tx_2,
						const bool useHalfRange);

		     
                       __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hess_3d_xmm2r8_a( double * __restrict __ATTR_ALIGN__(16) H_0,
		                                   double * __restrict __ATTR_ALIGN__(16) H_1,
						   double * __restrict __ATTR_ALIGN__(16) H_2,
						   double * __restrict __ATTR_ALIGN__(16) H_3,
						   double * __restrict __ATTR_ALIGN__(16) H_4,
						   double * __restrict __ATTR_ALIGN__(16) H_5,
						   double * __restrict __ATTR_ALIGN__(16) H_6,
						   double * __restrict __ATTR_ALIGN__(16) H_7,
						   double * __restrict __ATTR_ALIGN__(16) H_8,
						   const __m128d x_0,
						   const __m128d x_1,
						   const __m128d x_2,
						   const __m128d rx_0,
						   const __m128d rx_1,
						   const __m128d rx_2,
						   const __m128d tx_0,
						   const __m128d tx_1,
						   const __m128d tx_2,
						   const bool useHalfRange); 


		     
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hess_3d_xmm4r4(__m128 &H_0,
		                                __m128 &H_1,
						__m128 &H_2,
						__m128 &H_3,
						__m128 &H_4,
						__m128 &H_5,
						__m128 &H_6,
						__m128 &H_7,
						__m128 &H_8,
						const __m128 x_0,
						const __m128 x_1,
						const __m128 x_2,
						const __m128 rx_0,
						const __m128 rx_1,
						const __m128 rx_2,
						const __m128 tx_0,
						const __m128 tx_1,
						const __m128 tx_2,
						const bool useHalfRange); 
		     
                    
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hess_3d_xmm4r4_a(float * __restrict __ATTR_ALIGN__(16) H_0,
		                                   float * __restrict __ATTR_ALIGN__(16) H_1,
						   float * __restrict __ATTR_ALIGN__(16) H_2,
						   float * __restrict __ATTR_ALIGN__(16) H_3,
						   float * __restrict __ATTR_ALIGN__(16) H_4,
						   float * __restrict __ATTR_ALIGN__(16) H_5,
						   float * __restrict __ATTR_ALIGN__(16) H_6,
						   float * __restrict __ATTR_ALIGN__(16) H_7,
						   float * __restrict __ATTR_ALIGN__(16) H_8,
						   const __m128 x_0,
						   const __m128 x_1,
						   const __m128 x_2,
						   const __m128 rx_0,
						   const __m128 rx_1,
						   const __m128 rx_2,
						   const __m128 tx_0,
						   const __m128 tx_1,
						   const __m128 tx_2,
						   const bool useHalfRange); 


		     
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void range_hess_3d_xmm4r4_u(float * __restrict H_0,
		                                   float * __restrict H_1,
						   float * __restrict H_2,
						   float * __restrict H_3,
						   float * __restrict H_4,
						   float * __restrict H_5,
						   float * __restrict H_6,
						   float * __restrict H_7,
						   float * __restrict H_8,
						   const __m128 x_0,
						   const __m128 x_1,
						   const __m128 x_2,
						   const __m128 rx_0,
						   const __m128 rx_1,
						   const __m128 rx_2,
						   const __m128 tx_0,
						   const __m128 tx_1,
						   const __m128 tx_2,
						   const bool useHalfRange); 



 

 
/*CART2RUVGENCPP A C++ function to convert a Cartesian point into range,
 *           and direction cosines, possibly including the w component.
 *
 *INPUTS: retData A pointer to an array of doubles with 3 elements to
 *                hold the result in [r;u;v] order or with 4 elements if
 *                includeW is true to hold [r;u;v;w].
 *             zC The 3X1 Cartesian points [x;y;z] to be converted.
 *   useHalfRange A boolean value specifying whether the bistatic (round-
 *                trip) range value has been divided by two. 
 *            zTx The 3X1 [x;y;z] location vector of the transmitter in
 *                global Cartesian coordinates.
 *            zRx The 3X1 [x;y;z] location vector of the receiver in global
 *                Cartesian coordinates.
 *             M  A 3X3 rotation matrices to go from the alignment of the
 *                global coordinate system to that at the receiver. It is
 *                stored one columns after the other, consistent with how
 *                Matlab indexes matrices.
 *       includeW A boolean value indicating whether retData has space for
 *                a fourth component and the fourth component should be
 *                included.
 *
 *OUTPUTS: None. The results are placed in retData.
 *
 *See the comments to the Matlab function Cart2ruv for more information
 *on how this function works.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *@@Modified by Bernard Gingold, on May 2022
 **/

	 	    
	             
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_ruv_xmm2r8(__m128d &r,
		                              __m128d &u,
					      __m128d &v,
					      __m128d &w,
					      const __m128d C_x,
					      const __m128d C_y,
					      const __m128d C_z,
					      const __m128d T_x,
					      const __m128d T_y,
					      const __m128d T_z,
					      const __m128d R_x,
					      const __m128d R_y,
					      const __m128d R_z,
					      const __m128d * __restrict __ATTR_ALIGN__(16) M, //flattened 3x3 matrix
					      const bool useHalfRange); 

		     
                        __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_ruv_xmm2r8_a(double * __restrict __ATTR_ALIGN__(16) r,
		                                double * __restrict __ATTR_ALIGN__(16) u,
					        double * __restrict __ATTR_ALIGN__(16) v,
					        double * __restrict __ATTR_ALIGN__(16) w,
					      const __m128d C_x,
					      const __m128d C_y,
					      const __m128d C_z,
					      const __m128d T_x,
					      const __m128d T_y,
					      const __m128d T_z,
					      const __m128d R_x,
					      const __m128d R_y,
					      const __m128d R_z,
					      const __m128d * __restrict __ATTR_ALIGN__(16) M, //flattened 3x3 matrix
					      const bool useHalfRange);

		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_ruv_xmm2r8_u(double * __restrict r,
		                                double * __restrict u,
					        double * __restrict v,
					        double * __restrict w,
					      const __m128d C_x,
					      const __m128d C_y,
					      const __m128d C_z,
					      const __m128d T_x,
					      const __m128d T_y,
					      const __m128d T_z,
					      const __m128d R_x,
					      const __m128d R_y,
					      const __m128d R_z,
					      const __m128d * __restrict __ATTR_ALIGN__(16) M, //flattened 3x3 matrix
					      const bool useHalfRange); 


		     
                      __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_ruv_xmm4r4(__m128 &r,
		                               __m128 &u,
					       __m128 &v,
					       __m128 &w,
					      const __m128 C_x,
					      const __m128 C_y,
					      const __m128 C_z,
					      const __m128 T_x,
					      const __m128 T_y,
					      const __m128 T_z,
					      const __m128 R_x,
					      const __m128 R_y,
					      const __m128 R_z,
					      const __m128 * __restrict __ATTR_ALIGN__(16) M, //flattened 3x3 matrix
					      const bool useHalfRange); 

		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_ruv_xmm4r4_a(float * __restrict __ATTR_ALIGN__(16) r,
		                                 float * __restrict __ATTR_ALIGN__(16) u,
					         float * __restrict __ATTR_ALIGN__(16) v,
					         float * __restrict __ATTR_ALIGN__(16) w,
					         const __m128 C_x,
					         const __m128 C_y,
					         const __m128 C_z,
					         const __m128 T_x,
					         const __m128 T_y,
					         const __m128 T_z,
					         const __m128 R_x,
					         const __m128 R_y,
					         const __m128 R_z,
					         const __m128 * __restrict __ATTR_ALIGN__(16) M, //flattened 3x3 matrix
					         const bool useHalfRange);


		      
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_ruv_xmm4r4_u(float * __restrict r,
		                                 float * __restrict u,
					         float * __restrict v,
					         float * __restrict w,
					         const __m128 C_x,
					         const __m128 C_y,
					         const __m128 C_z,
					         const __m128 T_x,
					         const __m128 T_y,
					         const __m128 T_z,
					         const __m128 R_x,
					         const __m128 R_y,
					         const __m128 R_z,
					         const __m128 * __restrict __ATTR_ALIGN__(16) M, //flattened 3x3 matrix
					         const bool useHalfRange); 


		     
/*CART2SPHEREGENCPP A C++ function to convert Cartesian points to bistatic
 *            range, azimuth and elevation.
 *
 *INPUTS: retData A pointer to an array of doubles with 3 elements to
 *                hold the result in [range;azimuth;elevation]. order.
 *     cartPoints A pointer to the 3X1 Cartesian points [x;y;z] to be
 *                converted.
 *     systemType An integer specifying the axis from which the angles are
 *                measured. Possible values are
 *                0 Azimuth is measured counterclockwise from the x-axis in
 *                  the x-y plane. Elevation is measured up from the x-y
 *                  plane (towards the z-axis). This is consistent with
 *                  common spherical coordinate systems for specifying
 *                  longitude (azimuth) and geocentric latitude
 *                  (elevation).
 *                1 Azimuth is measured counterclockwise from the z-axis in
 *                  the z-x plane. Elevation is measured up from the z-x
 *                  plane (towards the y-axis). This is consistent with
 *                  some spherical coordinate systems that use the z-axis
 *                  as the boresight direction of the radar.
 *                2 This is the same as 0 except instead of being given
 *                  elevation, one desires the angle away from the z-axis,
 *                  which is (pi/2-elevation).
 *                3 This is the same as 0 except azimuth is measured
 *                  clockwise from the y-axis in the x-y plane instead of
 *                  counterclockwise from the x-axis. This coordinate
 *                  system often arises when given "bearings" in a local
 *                  East-North-Up coordinate system, where the bearing
 *                  directions are measured East of North.
 *   useHalfRange A boolean value specifying whether the bistatic (round-
 *                trip) range value has been divided by two. 
 *            zTx The 3X1 [x;y;z] location vector of the transmitter in
 *                global Cartesian coordinates.
 *            zRx The 3X1 [x;y;z] location vector of the receiver in global
 *                Cartesian coordinates.
 *             M  A 3X3  rotation matrices to go from the alignment of the
 *                global coordinate system to that at the receiver. It is
 *                stored one columns after the other, consistent with how
 *                Matlab indexes matrices.
 *
 *OUTPUTS: None. The results are placed in retData.
 *
 *See the comments to the Matlab function Cart2Sphere for more information
 *on how this function works.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *@@Modified by Bernard Gingold, on May 2022
 **/

	             
                          __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_sphere_xmm2r8(__m128d &range,
		                                 __m128d &az,
						 __m128d &elev,
						 const __m128d C_x,
						 const __m128d C_y,
						 const __m128d C_z,
						 const __m128d T_x,
						 const __m128d T_y,
						 const __m128d T_z,
						 const __m128d R_x,
						 const __m128d R_y,
						 const __m128d R_z,
						 const __m128d * __restrict __ATTR_ALIGN__(16) M,
						 const int sysType,
						 const bool useHalfRange); 


		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_sphere_xmm2r8_a(double * __restrict __ATTR_ALIGN__(16) range,
		                                   double * __restrict __ATTR_ALIGN__(16) az,
						   double * __restrict __ATTR_ALIGN__(16) elev,
						   const __m128d C_x,
						   const __m128d C_y,
						   const __m128d C_z,
						   const __m128d T_x,
						   const __m128d T_y,
						   const __m128d T_z,
						   const __m128d R_x,
						   const __m128d R_y,
						   const __m128d R_z,
						   const __m128d * __restrict __ATTR_ALIGN__(16) M,
						   const int sysType,
						   const bool useHalfRange);


		     
                          __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_sphere_xmm2r8_u(double * __restrict  range,
		                                   double * __restrict  az,
						   double * __restrict  elev,
						   const __m128d C_x,
						   const __m128d C_y,
						   const __m128d C_z,
						   const __m128d T_x,
						   const __m128d T_y,
						   const __m128d T_z,
						   const __m128d R_x,
						   const __m128d R_y,
						   const __m128d R_z,
						   const __m128d * __restrict __ATTR_ALIGN__(16) M,
						   const int sysType,
						   const bool useHalfRange); 


		     
                          __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_sphere_xmm4r4(__m128 &range,
		                                  __m128 &az,
						  __m128 &elev,
						  const __m128 C_x,
						  const __m128 C_y,
						  const __m128 C_z,
						  const __m128 T_x,
						  const __m128 T_y,
						  const __m128 T_z,
						  const __m128 R_x,
						  const __m128 R_y,
						  const __m128 R_z,
						  const __m128 * __restrict __ATTR_ALIGN__(16) M,
						  const int sysType,
						  const bool useHalfRange); 

		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_sphere_xmm4r4_a(float * __restrict __ATTR_ALIGN__(16) range,
		                                    float * __restrict __ATTR_ALIGN__(16) az,
						    float * __restrict __ATTR_ALIGN__(16) elev,
						    const __m128 C_x,
						    const __m128 C_y,
						    const __m128 C_z,
						    const __m128 T_x,
						    const __m128 T_y,
						    const __m128 T_z,
						    const __m128 R_x,
						    const __m128 R_y,
						    const __m128 R_z,
						    const __m128 * __restrict __ATTR_ALIGN__(16) M,
						    const int sysType,
						    const bool useHalfRange); 
 

		     
                         __ATTR_HOT__
                      __ATTR_VECTORCALL__
		      void cart_to_sphere_xmm4r4_u(float * __restrict  range,
		                                    float * __restrict  az,
						    float * __restrict  elev,
						    const __m128 C_x,
						    const __m128 C_y,
						    const __m128 C_z,
						    const __m128 T_x,
						    const __m128 T_y,
						    const __m128 T_z,
						    const __m128 R_x,
						    const __m128 R_y,
						    const __m128 R_z,
						    const __m128 * __restrict __ATTR_ALIGN__(16) M,
						    const int sysType,
						    const bool useHalfRange) 
 
 
 
      }

}











#endif /*__GMS_RANGE_RATE_SSE_H__*/
