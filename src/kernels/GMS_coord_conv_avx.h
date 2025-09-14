
#ifndef __GMS_COORD_CONV_AVX_H__
#define __GMS_COORD_CONV_AVX_H__ 120520221511

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
*OF RECIPIENT IN THE USE OF THE SOFTWARE.*/

namespace file_info {

 const unsigned int gGMS_COORD_CONV_AVX_MAJOR = 1U;
 const unsigned int gGMS_COORD_CONV_AVX_MINOR = 0U;
 const unsigned int gGMS_COORD_CONV_AVX_MICRO = 0U;
 const unsigned int gGMS_COORD_CONV_AVX_FULLVER =
  1000U*gGMS_COORD_CONV_AVX_MAJOR+100U*gGMS_COORD_CONV_AVX_MINOR+10U*gGMS_COORD_CONV_AVX_MICRO;
 const char * const pgGMS_COORD_CONV_AVX_CREATION_DATE = "12-05-2022 15:11 +00200 (THR 12 MAY 2022 15:11 GMT+2)";
 const char * const pgGMS_COORD_CONV_AVX_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const pgGMS_COORD_CONV_AVX_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const pgGMS_COORD_CONV_AVX_SYNOPSIS      = "AVX based Coordinate-System conversion functions (vectorized)."


}


#include <immintrin.h>
#include "GMS_config.h"


namespace gms {

        namespace math {

/*
   Adapted from: *GETENUAXISDIRECTIONS A C++ function to compute the unit basis vectors of
 *              an East-North-Up coordinate system.
 *
 *INPUTS: u A pointer to an array in which the unit vectors for the ENU
 *          axes are placed. The first three elements are the first vector,
 *          the next three the second one, etc. If justVertical=true, then
 *          only a single vector, the vertical is returned.
 * plhPoint A length 2 array at which the axes are to be found given in
 *          terms of [latitude;longitude] with the geodetic latitude and
 *          longitude in radians and the height in meters. The latitude
 *          should be between -pi/2 and pi/2.
 * justVertical A boolean parameter. If true then u and c only for the Up
 *          direction will be returned.
 *
 *March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *May   2022 Bernard Gingold, manual vectorization, beniekg@gmail.com
  
*/
                      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
	               void ENU_axis_dir_ymm4r8(__m256d * __restrict __ATTR_ALIGN__(32) u_0, //2-element array of __m256d
                                               __m256d * __restrict __ATTR_ALIGN__(32) u_1, //2-element array of __m256d
                                               __m256d * __restrict __ATTR_ALIGN__(32) u_2, //2-element array of __m256d
                                               const __m256d lat,
                                               const __m256d lon); 


                      __ATTR_REGCALL__
                     __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		       void ENU_axis_dir_ymm8r4(__m256 * __restrict __ATTR_ALIGN__(32) u_0, //2-element array of __m256
                                                __m256 * __restrict __ATTR_ALIGN__(32) u_1, //2-element array of __m256
                                                __m256 * __restrict __ATTR_ALIGN__(32) u_2, //2-element array of __m256
                                                const __m256 lat,
                                                const __m256 lon); 

/**GETENUAXESDIRMAG A C++ function to compute the basis vectors of an East-
 *              North-Up coordinate system as well as the magnitudes of the
 *              unnormalized vectors. The unnormalized vectors are
 *              gradients with respect to ellipsoidal coordinates.
 *
 *INPUTS: u A pointer to an array in which the unit vectors for the ENU
 *          axes are placed. The first three elements are the first vector,
 *          the next three the second one, etc.
 *        c A pointer to a length-3 array in which the magnitudes of the
 *          unnormalized vectors are placed.
 * plhPoint A length 2 array at which the axes are to be found given in
 *          terms of [latitude;longitude] with the geodetic latitude and
 *          longitude in radians and the height in meters. The latitude
 *          should be between -pi/2 and pi/2.
 *        a The semi-major axis of the reference ellipsoid.
 *        f The flattening factor of the reference ellipsoid
 *
 *March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *May   2022 Bernard Gingold, manual vectorization, beniekg@gmail.com
 **/
                      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		       void ENU_axis_dir_mag_ymm4r8(__m256d * __restrict __ATTR_ALIGN__(32) u_0, //3-element array of __m256d
                                                   __m256d * __restrict __ATTR_ALIGN__(32) u_1, //3-element array of __m256d
                                                   __m256d * __restrict __ATTR_ALIGN__(32) u_2, //3-element array of __m256d
                                                   __m256d * __restrict __ATTR_ALIGN__(32) c,   //3-element array of __m256d
                                                   const __m256d lat,
                                                   const __m256d lon,
                                                   const __m256d h,
                                                   const __m256d a,
                                                   const __m256d f); 


                      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void ENU_axis_dir_mag_ymm8r4(__m256 * __restrict __ATTR_ALIGN__(32) u_0, //3-element array of __m256d
                                                   __m256 * __restrict __ATTR_ALIGN__(32) u_1, //3-element array of __m256d
                                                   __m256 * __restrict __ATTR_ALIGN__(32) u_2, //3-element array of __m256d
                                                   __m256 * __restrict __ATTR_ALIGN__(32) c,   //3-element array of __m256d
                                                   const __m256 lat,
                                                   const __m256 lon,
                                                   const __m256 h,
                                                   const __m256 a,
                                                   const __m256 f); 
		     

		      __ATTR_REGCALL__
                     __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void ellips2Cart_ymm4r8(const __m256d lat,
		                              const __m256d lon,
					      const __m256d h,
					      const __m256d a,
					      const __m256d f,
					      __m256d * __restrict __ATTR_ALIGN__(32) cp); 


		      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline     
                      void ellips2Cart_ymm8r4(const __m256 lat,
		                               const __m256 lon,
					       const __m256 h,
					       const __m256 a,
					       const __m256 f,
					       __m256 * __restrict __ATTR_ALIGN__(32) cp); 


		   


         
                                                   
                                             


                 

     } //math

}// gms











#endif /*__GMS_COORDINATES_AVX_H__*/
