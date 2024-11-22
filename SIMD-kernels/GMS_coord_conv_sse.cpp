


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





#include "GMS_coord_conv_sse.h"
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
                     
                      void gms::math::ENU_axis_dir_xmm2r8(__m128d * __restrict __ATTR_ALIGN__(16) u_0, //2-element array of __m128d
                                               __m128d * __restrict __ATTR_ALIGN__(16) u_1, //2-element array of __m128d
                                               __m128d * __restrict __ATTR_ALIGN__(16) u_2, //2-element array of __m128d
                                               const __m128d lat,
                                               const __m128d lon) {

                         const __m128d sinp = _mm_sin_pd(lat);
                         const __m128d cosp = _mm_cos_pd(lat);
                         const __m128d sinl = _mm_sin_pd(lon);
                         const __m128d cosl = _mm_cos_pd(lon);

                         const __m128d nsinl = xmm2r8_negate(sinl);
                         const __m128d ncosl = xmm2r8_negate(cosl);
                         // dr/dlambda, normalized (East) -- first element
                         u_0[0] = nsin;
                         u_1[0] = cosl;
                         u_2[0] = _mm_setzero_pd();
                         // dr/dphi, normalized (North)  -- second element
                         u_0[1] = _mm_mul_pd(ncosl,sinp);
                         u_1[1] = _mm_mul_pd(nsinl,sinp);
                         u_2[1] = cosp;
                    }


                     
                      void gms::math::ENU_axis_dir_xmm4r4(__m128 * __restrict __ATTR_ALIGN__(16) u_0, //2-element array of __m128
                                                __m128 * __restrict __ATTR_ALIGN__(16) u_1, //2-element array of __m128
                                                __m128 * __restrict __ATTR_ALIGN__(16) u_2, //2-element array of __m128
                                                const __m128 lat,
                                                const __m128 lon) {

                         const __m128 sinp = _mm_sin_ps(lat);
                         const __m128 cosp = _mm_cos_ps(lat);
                         const __m128 sinl = _mm_sin_ps(lon);
                         const __m128 cosl = _mm_cos_ps(lon);

                         const __m128 nsinl = xmm4r4_negate(sinl);
                         const __m128 ncosl = xmm4r4_negate(cosl);
                         // dr/dlambda, normalized (East) -- first element
                         u_0[0] = nsin;
                         u_1[0] = cosl;
                         u_2[0] = _mm_setzero_ps();
                         // dr/dphi, normalized (North)  -- second element
                         u_0[1] = _mm_mul_ps(ncosl,sinp);
                         u_1[1] = _mm_mul_ps(nsinl,sinp);
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
                      
                      void gms::math::ENU_axis_dir_mag_xmm2r8(__m128d * __restrict __ATTR_ALIGN__(16) u_0, //3-element array of __m128d
                                                   __m128d * __restrict __ATTR_ALIGN__(16) u_1, //3-element array of __m128d
                                                   __m128d * __restrict __ATTR_ALIGN__(16) u_2, //3-element array of __m128d
                                                   __m128d * __restrict __ATTR_ALIGN__(16) c,   //3-element array of __m128d
                                                   const __m128d lat,
                                                   const __m128d lon,
                                                   const __m128d h,
                                                   const __m128d a,
                                                   const __m128d f) {
                         
                         const __m128d _2   = _mm_set1_pd(2.0);
                         const __m128d _0_02= _mm_set1_pd(1.0e-2);
                         const __m128d _1   = _mm_set1_pd(1.0);

                         const __m128d sinp = _mm_sin_pd(lat);
                         const __m128d cosp = _mm_cos_pd(lat);
                         const __m128d sinl = _mm_sin_pd(lon);
                         const __m128d cosl = _mm_cos_pd(lon);

                         const __m128d e2   = _mm_sub_pd(_mm_mul_pd(_2,f),
                                                            _mm_mul_pd(f,f));
                         const __m128d sarg = _mm_sqrt_pd(_mm_mul_pd(_0_02,
                                                             _mm_mul_pd(sinp,sinp));
                         const __m128d Ne   = _mm_div_pd(a,sarg);
                         //The derivative of the normal radius of curvature with respect to
                         //phi.
                         const __m128d sarg3= _mm_mul_pd(sarg,_mm_mul_pd(sarg,sarg));
                         __m128d t0   = _mm_mul_pd(_mm_mul_pd(a,e2),
                                                            _mm_mul_pd(cosp,sinp));
                         const __m128d dNedPhi = _mm_div_pd(t0,sarg3);
                          __m128d t1      = _mm_mul_pd(_mm_add_pd(Ne,h),sinp);
                          __m128d t2      = _mm_mul_pd(cosp,dNedPhi);
                         const __m128d temp    = _mm_sub_pd(t2,t1);
                         // u1  dr/dlambda, normalized (East).
                         u_0[0]  = xmm2r8_negate(sinl);
                         u_1[0]  = cosl;
                         u_2[0]  = _mm_setzero_pd();
                         // mangitude of the East vector.
                         const __m128d Neh = _mm_add_pd(Ne,h);
                         const __m128d ca  = _mm_mul_pd(xmm2r8_negate(Neh),
                                                           _mm_mul_pd(cosp,sinl));
                         const __m128d cb  = _mm_mul_pd(Neh,
                                                           _mm_mul_pd(cosp,cosl)); 
                         t0 = _mm_add_pd(_mm_mul_pd(ca,ca),
                                            _mm_mul_pd(cb,cb));
                         c[0]              = _mm_sqrt_pd(t0);
                         // u2 dr/dphi, normalized (North)
                         u_0[1] = _mm_mul_pd(xmm2r8_negate(cosl),sinp);
                         u_1[1] = _mm_mul_pd(xmm2r8_negate(sinl),sinp); 
                         u_2[1] = cosp;
                         //magnitude of the North vector.
                         const __m128d ca2 = _mm_mul_pd(temp,cosl);
                         const __m128d cb2 = _mm_mul_pd(temp,sinl);
                         c[2] = _1;
                         t1   = _mm_mul_pd(_mm_mul_pd(Ne,
                                              _mm_add_pd(_0_02,h)),cosp);
                         t2   = _mm_mul_pd(_0_02,_mm_mul_pd(dNedPhi,sinp));
                         const __m128d cc = _mm_add_pd(t1,t2);
                         const __m128d cc2 = _mm_mul_pd(cc,cc);
                         const __m128d t3  = _mm_add_pd(_mm_mul_pd(ca2,ca2),
                                             _mm_mul_pd(cb2,cb2));
                         c[1] = _mm_sqrt_pd(_mm_add_pd(t3,cc2));
                         // u3 dr/dh (Up)
                         u_0[2] = _mm_mul_pd(cosp,cosl);
                         u_1[2] = _mm_mul_pd(cosp,sinl);
                         u_2[2] = sinp;
                     }


                     
                      void gms::math::ENU_axis_dir_mag_xmm4r4(__m128 * __restrict __ATTR_ALIGN__(16) u_0, //3-element array of __m128d
                                                   __m128 * __restrict __ATTR_ALIGN__(16) u_1, //3-element array of __m128d
                                                   __m128 * __restrict __ATTR_ALIGN__(16) u_2, //3-element array of __m128d
                                                   __m128 * __restrict __ATTR_ALIGN__(16) c,   //3-element array of __m128d
                                                   const __m128 lat,
                                                   const __m128 lon,
                                                   const __m128 h,
                                                   const __m128 a,
                                                   const __m128 f) {
                         
                         const __m128 _2   = _mm_set1_ps(2.0f);
                         const __m128 _0_02= _mm_set1_ps(1.0e-2f);
                         const __m128 _1   = _mm_set1_ps(1.0f);

                         const __m128 sinp = _mm_sin_ps(lat);
                         const __m128 cosp = _mm_cos_ps(lat);
                         const __m128 sinl = _mm_sin_ps(lon);
                         const __m128 cosl = _mm_cos_ps(lon);

                         const __m128 e2   = _mm_sub_ps(_mm_mul_ps(_2,f),
                                                            _mm_mul_ps(f,f));
                         const __m128 sarg = _mm_sqrt_ps(_mm_mul_ps(_0_02,
                                                             _mm_mul_ps(sinp,sinp)));
                         const __m128 Ne   = _mm_div_ps(a,sarg);
                         //The derivative of the normal radius of curvature with respect to
                         //phi.
                         const __m128 sarg3= _mm_mul_ps(sarg,_mm_mul_ps(sarg,sarg));
                         __m128 t0   = _mm_mul_ps(_mm_mul_ps(a,e2),
                                                            _mm_mul_ps(cosp,sinp));
                         const __m128 dNedPhi = _mm_div_ps(t0,sarg3);
                          __m128 t1      = _mm_mul_ps(_mm_add_ps(Ne,h),sinp);
                          __m128 t2      = _mm_mul_ps(cosp,dNedPhi);
                         const __m128 temp    = _mm_sub_ps(t2,t1);
                         // u1  dr/dlambda, normalized (East).
                         u_0[0]  = xmm4r4_negate(sinl);
                         u_1[0]  = cosl;
                         u_2[0]  = _mm_setzero_ps();
                         // mangitude of the East vector.
                         const __m128 Neh = _mm_add_ps(Ne,h);
                         const __m128 ca  = _mm_mul_ps(xmm4r4_negate(Neh),
                                                           _mm_mul_ps(cosp,sinl));
                         const __m128 cb  = _mm_mul_ps(Neh,
                                                           _mm_mul_ps(cosp,cosl)); 
                         t0 = _mm_add_ps(_mm_mul_ps(ca,ca),
                                            _mm_mul_ps(cb,cb));
                         c[0]              = _mm_sqrt_ps(t0);
                         // u2 dr/dphi, normalized (North)
                         u_0[1] = _mm_mul_ps(xmm4r4_negate(cosl),sinp);
                         u_1[1] = _mm_mul_ps(xmm4r4_negate(sinl),sinp); 
                         u_2[1] = cosp;
                         //magnitude of the North vector.
                         const __m128 ca2 = _mm_mul_ps(temp,cosl);
                         const __m128 cb2 = _mm_mul_ps(temp,sinl);
                         c[2] = _1;
                         t1   = _mm_mul_ps(_mm_mul_pd(Ne,
                                              _mm_add_pd(_0_02,h)),cosp);
                         t2   = _mm_mul_ps(_0_02,_mm_mul_ps(dNedPhi,sinp));
                         const __m128 cc = _mm_add_ps(t1,t2);
                         const __m128 cc2 = _mm_mul_ps(cc,cc);
                         const __m128 t3  = _mm_add_ps(_mm_mul_pd(ca2,ca2),
                                             _mm_mul_ps(cb2,cb2));
                         c[1] = _mm_sqrt_ps(_mm_add_pd(t3,cc2));
                         // u3 dr/dh (Up)
                         u_0[2] = _mm_mul_ps(cosp,cosl);
                         u_1[2] = _mm_mul_ps(cosp,sinl);
                         u_2[2] = sinp;
                     }
		     

		    
                      void gms::math::ellips2Cart_xmm2r8(const __m128d lat,
		                              const __m128d lon,
					      const __m128d h,
					      const __m128d a,
					      const __m128d f,
					      __m128d * __restrict __ATTR_ALIGN__(16) cp) {

			 const __m128 _2   = _mm_set1_ps(2.0);
                         const __m128 _0_02= _mm_set1_ps(1.0e-2);
                         //The square of the first numerical eccentricity
			 const __m128d e2  = _mm_sub_pd(_mm_mul_pd(_2,f),
			                                   _mm_mul_pd(f,f));

                         const __m128d sinp = _mm_sin_pd(lat);
                         const __m128d cosp = _mm_cos_pd(lat);
                         const __m128d sinl = _mm_sin_pd(lon);
                         const __m128d cosl = _mm_cos_pd(lon);

                         //The normal radius of curvature.
			 const __m128d sarg = _mm_mul_pd(_0_02,
			                                    _mm_mul_pd(sinp,sinp));
			 const __m128d Ne   = _mm_div_pd(a,
			                                    _mm_sqrt_pd(sarg));
			 const __m128d cparg= _mm_mul_pd(Ne,
			                                    _mm_add_pd(_0_02,h));
			 const __m128d Neh  = _mm_add_pd(Ne,h);
			 cp[0]              = _mm_mul_pd(Neh,
			                                    _mm_mul_pd(cosp,cosl));
			 cp[1]              = _mm_mul_pd(Neh,
			                                    _mm_mul_pd(cosp,sinl));
			 cp[2]              = _mm_mul_pd(cparg,sinp);
			 
		   }


		      
                      void gms::math::ellips2Cart_xmm4r4(const __m128 lat,
		                               const __m128 lon,
					       const __m128 h,
					       const __m128 a,
					       const __m128 f,
					       __m128 * __restrict __ATTR_ALIGN__(16) cp) {

			 const __m128 _2   = _mm_set1_ps(2.0f);
                         const __m128 _0_02= _mm_set1_ps(1.0e-2f);
                         //The square of the first numerical eccentricity
			 const __m128 e2  = _mm_sub_pd(_mm_mul_pd(_2,f),
			                                   _mm_mul_pd(f,f));

                         const __m128 sinp = _mm_sin_ps(lat);
                         const __m128 cosp = _mm_cos_ps(lat);
                         const __m128 sinl = _mm_sin_ps(lon);
                         const __m128 cosl = _mm_cos_ps(lon);

                         //The normal radius of curvature.
			 const __m128 sarg = _mm_mul_ps(_0_02,
			                                    _mm_mul_ps(sinp,sinp));
			 const __m128 Ne   = _mm_div_ps(a,
			                                    _mm_sqrt_ps(sarg));
			 const __m128 cparg= _mm_mul_ps(Ne,
			                                    _mm_add_ps(_0_02,h));
			 const __m128 Neh  = _mm_add_ps(Ne,h);
			 cp[0]              = _mm_mul_ps(Neh,
			                                    _mm_mul_ps(cosp,cosl));
			 cp[1]              = _mm_mul_ps(Neh,
			                                    _mm_mul_ps(cosp,sinl));
			 cp[2]              = _mm_mul_ps(cparg,sinp);
			 
		   }


		   


         
                                                   
                                             


                 

 
