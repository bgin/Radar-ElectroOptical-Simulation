
#ifndef __GMS_COORD_CONV_AVX512_HPP__
#define __GMS_COORD_CONV_AVX512_HPP__ 120520221311

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
#if (USE_SLEEF_LIB) == 1
#include "GMS_sleefsimddp.hpp"
#include "GMS_sleefsimdsp.hpp"
#endif
#include "GMS_simd_utils.hpp"

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
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline     
                      void ENU_axis_dir_zmm8r8(__m512d * __restrict __ATTR_ALIGN__(64) u_0, //2-element array of __m512d
                                               __m512d * __restrict __ATTR_ALIGN__(64) u_1, //2-element array of __m512d
                                               __m512d * __restrict __ATTR_ALIGN__(64) u_2, //2-element array of __m512d
                                               const __m512d lat,
                                               const __m512d lon) {
#if (USE_SLEEF_LIB) == 1
                         const __m512d sinp = xsin(lat);
                         const __m512d cosp = xcos(lat);
                         const __m512d sinl = xsin(lon);
                         const __m512d cosl = xcos(lon);
#else
                         const __m512d sinp = _mm512_sin_pd(lat);
                         const __m512d cosp = _mm512_cos_pd(lat);
                         const __m512d sinl = _mm512_sin_pd(lon);
                         const __m512d cosl = _mm512_cos_pd(lon);
#endif
                         const __m512d nsinl = zmm8r8_negate(sinl);
                         const __m512d ncosl = zmm8r8_negate(cosl);
                         // dr/dlambda, normalized (East) -- first element
                         u_0[0] = nsin;
                         u_1[0] = cosl;
                         u_2[0] = _mm512_setzero_pd();
                         // dr/dphi, normalized (North)  -- second element
                         u_0[1] = _mm512_mul_pd(ncosl,sinp);
                         u_1[1] = _mm512_mul_pd(nsinl,sinp);
                         u_2[1] = cosp;
                    }


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline     
                      void ENU_axis_dir_zmm16r4(__m512 * __restrict __ATTR_ALIGN__(64) u_0, //2-element array of __m512
                                                __m512 * __restrict __ATTR_ALIGN__(64) u_1, //2-element array of __m512
                                                __m512 * __restrict __ATTR_ALIGN__(64) u_2, //2-element array of __m512
                                                const __m512 lat,
                                                const __m512 lon) {
#if (USE_SLEEF_LIB) == 1
                         const __m512 sinp = xsinf(lat);
                         const __m512 cosp = xcosf(lat);
                         const __m512 sinl = xsinf(lon);
                         const __m512 cosl = xcosf(lon);
#else
                         const __m512 sinp = _mm512_sin_ps(lat);
                         const __m512 cosp = _mm512_cos_ps(lat);
                         const __m512 sinl = _mm512_sin_ps(lon);
                         const __m512 cosl = _mm512_cos_ps(lon);
#endif
                         const __m512 nsinl = zmm16r4_negate(sinl);
                         const __m512 ncosl = zmm16r4_negate(cosl);
                         // dr/dlambda, normalized (East) -- first element
                         u_0[0] = nsin;
                         u_1[0] = cosl;
                         u_2[0] = _mm512_setzero_ps();
                         // dr/dphi, normalized (North)  -- second element
                         u_0[1] = _mm512_mul_ps(ncosl,sinp);
                         u_1[1] = _mm512_mul_ps(nsinl,sinp);
                         u_2[1] = cosp;
                    }


                 

     } //math

}// gms











#endif /*__GMS_COORDINATES_AVX512_HPP__*/
