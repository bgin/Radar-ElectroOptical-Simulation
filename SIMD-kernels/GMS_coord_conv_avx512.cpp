


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




#include "GMS_coord_conv_avx512.h"
#if (USE_SLEEF_LIB) == 1
#include "GMS_sleefsimddp.hpp"
#include "GMS_sleefsimdsp.hpp"
#endif
#include "GMS_simd_utils.hpp"


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
                     
                      void gms::math::ENU_axis_dir_zmm8r8(__m512d * __restrict __ATTR_ALIGN__(64) u_0, //2-element array of __m512d
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


                  
                      void gms::math::ENU_axis_dir_zmm16r4(__m512 * __restrict __ATTR_ALIGN__(64) u_0, //2-element array of __m512
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
                     
                      void gms::math::ENU_axis_dir_mag_zmm8r8(__m512d * __restrict __ATTR_ALIGN__(64) u_0, //3-element array of __m512d
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


                   
                      void gms::math::ENU_axis_dir_mag_zmm16r4(__m512 * __restrict __ATTR_ALIGN__(64) u_0, //3-element array of __m512d
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
                                                             _mm512_mul_ps(sinp,sinp)));
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
		     

		     
                      void gms::math::ellips2Cart_zmm8r8(const __m512d lat,
		                              const __m512d lon,
					      const __m512d h,
					      const __m512d a,
					      const __m512d f,
					      __m512d * __restrict __ATTR_ALIGN__(64) cp) {

			 const __m512 _2   = _mm512_set1_ps(2.0);
                         const __m512 _0_02= _mm512_set1_ps(1.0e-2);
                         //The square of the first numerical eccentricity
			 const __m512d e2  = _mm512_sub_pd(_mm512_mul_pd(_2,f),
			                                   _mm512_mul_pd(f,f));
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
                         //The normal radius of curvature.
			 const __m512d sarg = _mm512_mul_pd(_0_02,
			                                    _mm512_mul_pd(sinp,sinp));
			 const __m512d Ne   = _mm512_div_pd(a,
			                                    _mm512_sqrt_pd(sarg));
			 const __m512d cparg= _mm512_mul_pd(Ne,
			                                    _mm512_add_pd(_0_02,h));
			 const __m512d Neh  = _mm512_add_pd(Ne,h);
			 cp[0]              = _mm512_mul_pd(Neh,
			                                    _mm512_mul_pd(cosp,cosl));
			 cp[1]              = _mm512_mul_pd(Neh,
			                                    _mm512_mul_pd(cosp,sinl));
			 cp[2]              = _mm512_mul_pd(cparg,sinp);
			 
		   }


		      
                      void gms::math::ellips2Cart_zmm16r4(const __m512 lat,
		                               const __m512 lon,
					       const __m512 h,
					       const __m512 a,
					       const __m512 f,
					       __m512 * __restrict __ATTR_ALIGN__(64) cp) {

			 const __m512 _2   = _mm512_set1_ps(2.0f);
                         const __m512 _0_02= _mm512_set1_ps(1.0e-2f);
                         //The square of the first numerical eccentricity
			 const __m512 e2  = _mm512_sub_pd(_mm512_mul_pd(_2,f),
			                                   _mm512_mul_pd(f,f));
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
                         //The normal radius of curvature.
			 const __m512 sarg = _mm512_mul_ps(_0_02,
			                                    _mm512_mul_ps(sinp,sinp));
			 const __m512 Ne   = _mm512_div_ps(a,
			                                    _mm512_sqrt_ps(sarg));
			 const __m512 cparg= _mm512_mul_ps(Ne,
			                                    _mm512_add_ps(_0_02,h));
			 const __m512 Neh  = _mm512_add_ps(Ne,h);
			 cp[0]              = _mm512_mul_ps(Neh,
			                                    _mm512_mul_ps(cosp,cosl));
			 cp[1]              = _mm512_mul_ps(Neh,
			                                    _mm512_mul_ps(cosp,sinl));
			 cp[2]              = _mm512_mul_ps(cparg,sinp);
			 
		   }


		   


         
                                                   
                                             


                 

 
