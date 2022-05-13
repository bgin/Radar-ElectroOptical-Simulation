
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
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline     
                      void ENU_axis_dir_mag_zmm8r8(__m512d * __restrict __ATTR_ALIGN__(64) u_0, //3-element array of __m512d
                                                   __m512d * __restrict __ATTR_ALIGN__(64) u_1, //3-element array of __m512d
                                                   __m512d * __restrict __ATTR_ALIGN__(64) u_2, //3-element array of __m512d
                                                   __m512d * __restrict __ATTR_ALIGN__(64) c,   //3-element array of __m512d
                                                   const __m512d lat,
                                                   const __m512d lon,
                                                   const __m512d h,
                                                   const __m512d a,
                                                   const __m512d f) {
                         
                         const __m512d _2   = _mm512_set1_pd(2.0);
                         const __m512d _0_02= _mm512_set1_pd(1.0e-2);
                         const __m512d _1   = _mm512_set1_pd(1.0);
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
                         const __m512d e2   = _mm512_sub_pd(_mm512_mul_pd(_2,f),
                                                            _mm512_mul_pd(f,f));
                         const __m512d sarg = _mm512_sqrt_pd(_mm512_mul_pd(_0_02,
                                                             _mm512_mul_pd(sinp,sinp));
                         const __m512d Ne   = _mm512_div_pd(a,sarg);
                         //The derivative of the normal radius of curvature with respect to
                         //phi.
                         const __m512d sarg3= _mm512_mul_pd(sarg,_mm512_mul_pd(sarg,sarg));
                         __m512d t0   = _mm512_mul_pd(_mm512_mul_pd(a,e2),
                                                            _mm512_mul_pd(cosp,sinp));
                         const __m512d dNedPhi = _mm512_div_pd(t0,sarg3);
                          __m512d t1      = _mm512_mul_pd(_mm512_add_pd(Ne,h),sinp);
                          __m512d t2      = _mm512_mul_pd(cosp,dNedPhi);
                         const __m512d temp    = _mm512_sub_pd(t2,t1);
                         // u1  dr/dlambda, normalized (East).
                         u_0[0]  = zmm8r8_negate(sinl);
                         u_1[0]  = cosl;
                         u_2[0]  = _mm512_setzero_pd();
                         // mangitude of the East vector.
                         const __m512d Neh = _mm512_add_pd(Ne,h);
                         const __m512d ca  = _mm512_mul_pd(zmm8r8_negate(Neh),
                                                           _mm512_mul_pd(cosp,sinl));
                         const __m512d cb  = _mm512_mul_pd(Neh,
                                                           _mm512_mul_pd(cosp,cosl)); 
                         t0 = _mm512_add_pd(_mm512_mul_pd(ca,ca),
                                            _mm512_mul_pd(cb,cb));
                         c[0]              = _mm512_sqrt_pd(t0);
                         // u2 dr/dphi, normalized (North)
                         u_0[1] = _mm512_mul_pd(zmm8r8_negate(cosl),sinp);
                         u_1[1] = _mm512_mul_pd(zmm8r8_negate(sinl),sinp); 
                         u_2[1] = cosp;
                         //magnitude of the North vector.
                         const __m512d ca2 = _mm512_mul_pd(temp,cosl);
                         const __m512d cb2 = _mm512_mul_pd(temp,sinl);
                         c[2] = _1;
                         t1   = _mm512_mul_pd(_mm512_mul_pd(Ne,
                                              _mm512_add_pd(_0_02,h)),cosp);
                         t2   = _mm512_mul_pd(_0_02,_mm512_mul_pd(dNedPhi,sinp));
                         const __m512d cc = _mm512_add_pd(t1,t2);
                         const __m512d cc2 = _mm512_mul_pd(cc,cc);
                         const __m512d t3  = _mm512_add_pd(_mm512_mul_pd(ca2,ca2),
                                             _mm512_mul_pd(cb2,cb2));
                         c[1] = _mm512_sqrt_pd(_mm512_add_pd(t3,cc2));
                         // u3 dr/dh (Up)
                         u_0[2] = _mm512_mul_pd(cosp,cosl);
                         u_1[2] = _mm512_mul_pd(cosp,sinl);
                         u_2[2] = sinp;
                     }


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline     
                      void ENU_axis_dir_mag_zmm16r4(__m512 * __restrict __ATTR_ALIGN__(64) u_0, //3-element array of __m512d
                                                   __m512 * __restrict __ATTR_ALIGN__(64) u_1, //3-element array of __m512d
                                                   __m512 * __restrict __ATTR_ALIGN__(64) u_2, //3-element array of __m512d
                                                   __m512 * __restrict __ATTR_ALIGN__(64) c,   //3-element array of __m512d
                                                   const __m512 lat,
                                                   const __m512 lon,
                                                   const __m512 h,
                                                   const __m512 a,
                                                   const __m512 f) {
                         
                         const __m512 _2   = _mm512_set1_ps(2.0f);
                         const __m512 _0_02= _mm512_set1_ps(1.0e-2f);
                         const __m512 _1   = _mm512_set1_ps(1.0f);
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
                         const __m512 e2   = _mm512_sub_ps(_mm512_mul_ps(_2,f),
                                                            _mm512_mul_ps(f,f));
                         const __m512 sarg = _mm512_sqrt_ps(_mm512_mul_ps(_0_02,
                                                             _mm512_mul_ps(sinp,sinp));
                         const __m512 Ne   = _mm512_div_ps(a,sarg);
                         //The derivative of the normal radius of curvature with respect to
                         //phi.
                         const __m512 sarg3= _mm512_mul_ps(sarg,_mm512_mul_ps(sarg,sarg));
                         __m512 t0   = _mm512_mul_ps(_mm512_mul_ps(a,e2),
                                                            _mm512_mul_ps(cosp,sinp));
                         const __m512 dNedPhi = _mm512_div_ps(t0,sarg3);
                          __m512 t1      = _mm512_mul_ps(_mm512_add_ps(Ne,h),sinp);
                          __m512 t2      = _mm512_mul_ps(cosp,dNedPhi);
                         const __m512 temp    = _mm512_sub_ps(t2,t1);
                         // u1  dr/dlambda, normalized (East).
                         u_0[0]  = zmm16r4_negate(sinl);
                         u_1[0]  = cosl;
                         u_2[0]  = _mm512_setzero_ps();
                         // mangitude of the East vector.
                         const __m512 Neh = _mm512_add_ps(Ne,h);
                         const __m512 ca  = _mm512_mul_ps(zmm16r4_negate(Neh),
                                                           _mm512_mul_ps(cosp,sinl));
                         const __m512 cb  = _mm512_mul_ps(Neh,
                                                           _mm512_mul_ps(cosp,cosl)); 
                         t0 = _mm512_add_ps(_mm512_mul_ps(ca,ca),
                                            _mm512_mul_ps(cb,cb));
                         c[0]              = _mm512_sqrt_ps(t0);
                         // u2 dr/dphi, normalized (North)
                         u_0[1] = _mm512_mul_ps(zmm16r4_negate(cosl),sinp);
                         u_1[1] = _mm512_mul_ps(zmm16r4_negate(sinl),sinp); 
                         u_2[1] = cosp;
                         //magnitude of the North vector.
                         const __m512 ca2 = _mm512_mul_ps(temp,cosl);
                         const __m512 cb2 = _mm512_mul_ps(temp,sinl);
                         c[2] = _1;
                         t1   = _mm512_mul_ps(_mm512_mul_pd(Ne,
                                              _mm512_add_pd(_0_02,h)),cosp);
                         t2   = _mm512_mul_ps(_0_02,_mm512_mul_ps(dNedPhi,sinp));
                         const __m512 cc = _mm512_add_ps(t1,t2);
                         const __m512 cc2 = _mm512_mul_ps(cc,cc);
                         const __m512 t3  = _mm512_add_ps(_mm512_mul_pd(ca2,ca2),
                                             _mm512_mul_ps(cb2,cb2));
                         c[1] = _mm512_sqrt_ps(_mm512_add_pd(t3,cc2));
                         // u3 dr/dh (Up)
                         u_0[2] = _mm512_mul_ps(cosp,cosl);
                         u_1[2] = _mm512_mul_ps(cosp,sinl);
                         u_2[2] = sinp;
                     }


         
                                                   
                                             


                 

     } //math

}// gms











#endif /*__GMS_COORDINATES_AVX512_HPP__*/
