
#ifndef __GMS_COORD_CONV_AVX512_H__
#define __GMS_COORD_CONV_AVX512_H__ 120520221311

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

 const unsigned int gGMS_COORD_CONV_AVX512_MAJOR = 1U;
 const unsigned int gGMS_COORD_CONV_AVX512_MINOR = 0U;
 const unsigned int gGMS_COORD_CONV_AVX512_MICRO = 0U;
 const unsigned int gGMS_COORD_CONV_AVX512_FULLVER =
  1000U*gGMS_COORD_CONV_AVX512_MAJOR+100U*gGMS_COORD_CONV_AVX512_MINOR+10U*gGMS_COORD_CONV_AVX512_MICRO;
 const char * const pgGMS_COORD_CONV_AVX512_CREATION_DATE = "12-05-2022 13:11 +00200 (THR 12 MAY 2022 13:11 GMT+2)";
 const char * const pgGMS_COORD_CONV_AVX512_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const pgGMS_COORD_CONV_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const pgGMS_COORD_CONV_AVX512_SYNOPSIS      = "AVX512 based Coordinate-System conversion functions (vectorized)."


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
		      void ENU_axis_dir_zmm8r8(__m512d * __restrict __ATTR_ALIGN__(64) u_0, //2-element array of __m512d
                                               __m512d * __restrict __ATTR_ALIGN__(64) u_1, //2-element array of __m512d
                                               __m512d * __restrict __ATTR_ALIGN__(64) u_2, //2-element array of __m512d
                                               const __m512d lat,
                                               const __m512d lon); 


                      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void ENU_axis_dir_zmm16r4(__m512 * __restrict __ATTR_ALIGN__(64) u_0, //2-element array of __m512
                                                __m512 * __restrict __ATTR_ALIGN__(64) u_1, //2-element array of __m512
                                                __m512 * __restrict __ATTR_ALIGN__(64) u_2, //2-element array of __m512
                                                const __m512 lat,
                                                const __m512 lon); 

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
		      void ENU_axis_dir_mag_zmm8r8(__m512d * __restrict __ATTR_ALIGN__(64) u_0, //3-element array of __m512d
                                                   __m512d * __restrict __ATTR_ALIGN__(64) u_1, //3-element array of __m512d
                                                   __m512d * __restrict __ATTR_ALIGN__(64) u_2, //3-element array of __m512d
                                                   __m512d * __restrict __ATTR_ALIGN__(64) c,   //3-element array of __m512d
                                                   const __m512d lat,
                                                   const __m512d lon,
                                                   const __m512d h,
                                                   const __m512d a,
                                                   const __m512d f); 
                                                   

                      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void ENU_axis_dir_mag_zmm16r4(__m512 * __restrict __ATTR_ALIGN__(64) u_0, //3-element array of __m512d
                                                   __m512 * __restrict __ATTR_ALIGN__(64) u_1, //3-element array of __m512d
                                                   __m512 * __restrict __ATTR_ALIGN__(64) u_2, //3-element array of __m512d
                                                   __m512 * __restrict __ATTR_ALIGN__(64) c,   //3-element array of __m512d
                                                   const __m512 lat,
                                                   const __m512 lon,
                                                   const __m512 h,
                                                   const __m512 a,
                                                   const __m512 f); 

		      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void ellips2Cart_zmm8r8(const __m512d lat,
		                              const __m512d lon,
					      const __m512d h,
					      const __m512d a,
					      const __m512d f,
					      __m512d * __restrict __ATTR_ALIGN__(64) cp); 

		      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void ellips2Cart_zmm16r4(const __m512 lat,
		                               const __m512 lon,
					       const __m512 h,
					       const __m512 a,
					       const __m512 f,
					       __m512 * __restrict __ATTR_ALIGN__(64) cp); 


		   


         
                                                   
                                             


                 

     } //math

}// gms











#endif /*__GMS_COORDINATES_AVX512_H__*/
