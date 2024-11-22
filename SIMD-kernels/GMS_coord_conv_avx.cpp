


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




#include "GMS_coor_conv_avx.h"
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
                      
                      void gms::math::ENU_axis_dir_ymm4r8(__m256d * __restrict __ATTR_ALIGN__(32) u_0, //2-element array of __m256d
                                               __m256d * __restrict __ATTR_ALIGN__(32) u_1, //2-element array of __m256d
                                               __m256d * __restrict __ATTR_ALIGN__(32) u_2, //2-element array of __m256d
                                               const __m256d lat,
                                               const __m256d lon) {

                         const __m256d sinp = _mm256_sin_pd(lat);
                         const __m256d cosp = _mm256_cos_pd(lat);
                         const __m256d sinl = _mm256_sin_pd(lon);
                         const __m256d cosl = _mm256_cos_pd(lon);

                         const __m256d nsinl = ymm4r8_negate(sinl);
                         const __m256d ncosl = ymm4r8_negate(cosl);
                         // dr/dlambda, normalized (East) -- first element
                         u_0[0] = nsin;
                         u_1[0] = cosl;
                         u_2[0] = _mm256_setzero_pd();
                         // dr/dphi, normalized (North)  -- second element
                         u_0[1] = _mm256_mul_pd(ncosl,sinp);
                         u_1[1] = _mm256_mul_pd(nsinl,sinp);
                         u_2[1] = cosp;
                    }


                       
                      void gms::math::ENU_axis_dir_ymm8r4(__m256 * __restrict __ATTR_ALIGN__(32) u_0, //2-element array of __m256
                                                __m256 * __restrict __ATTR_ALIGN__(32) u_1, //2-element array of __m256
                                                __m256 * __restrict __ATTR_ALIGN__(32) u_2, //2-element array of __m256
                                                const __m256 lat,
                                                const __m256 lon) {

                         const __m256 sinp = _mm256_sin_ps(lat);
                         const __m256 cosp = _mm256_cos_ps(lat);
                         const __m256 sinl = _mm256_sin_ps(lon);
                         const __m256 cosl = _mm256_cos_ps(lon);

                         const __m256 nsinl = ymm8r4_negate(sinl);
                         const __m256 ncosl = ymm8r4_negate(cosl);
                         // dr/dlambda, normalized (East) -- first element
                         u_0[0] = nsin;
                         u_1[0] = cosl;
                         u_2[0] = _mm256_setzero_ps();
                         // dr/dphi, normalized (North)  -- second element
                         u_0[1] = _mm256_mul_ps(ncosl,sinp);
                         u_1[1] = _mm256_mul_ps(nsinl,sinp);
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
                   
                      void gms::math::ENU_axis_dir_mag_ymm4r8(__m256d * __restrict __ATTR_ALIGN__(32) u_0, //3-element array of __m256d
                                                   __m256d * __restrict __ATTR_ALIGN__(32) u_1, //3-element array of __m256d
                                                   __m256d * __restrict __ATTR_ALIGN__(32) u_2, //3-element array of __m256d
                                                   __m256d * __restrict __ATTR_ALIGN__(32) c,   //3-element array of __m256d
                                                   const __m256d lat,
                                                   const __m256d lon,
                                                   const __m256d h,
                                                   const __m256d a,
                                                   const __m256d f) {
                         
                         const __m256d _2   = _mm256_set1_pd(2.0);
                         const __m256d _0_02= _mm256_set1_pd(1.0e-2);
                         const __m256d _1   = _mm256_set1_pd(1.0);

                         const __m256d sinp = _mm256_sin_pd(lat);
                         const __m256d cosp = _mm256_cos_pd(lat);
                         const __m256d sinl = _mm256_sin_pd(lon);
                         const __m256d cosl = _mm256_cos_pd(lon);

                         const __m256d e2   = _mm256_sub_pd(_mm256_mul_pd(_2,f),
                                                            _mm256_mul_pd(f,f));
                         const __m256d sarg = _mm256_sqrt_pd(_mm256_mul_pd(_0_02,
                                                             _mm256_mul_pd(sinp,sinp));
                         const __m256d Ne   = _mm256_div_pd(a,sarg);
                         //The derivative of the normal radius of curvature with respect to
                         //phi.
                         const __m256d sarg3= _mm256_mul_pd(sarg,_mm256_mul_pd(sarg,sarg));
                         __m256d t0   = _mm256_mul_pd(_mm256_mul_pd(a,e2),
                                                            _mm256_mul_pd(cosp,sinp));
                         const __m256d dNedPhi = _mm256_div_pd(t0,sarg3);
                          __m256d t1      = _mm256_mul_pd(_mm256_add_pd(Ne,h),sinp);
                          __m256d t2      = _mm256_mul_pd(cosp,dNedPhi);
                         const __m256d temp    = _mm256_sub_pd(t2,t1);
                         // u1  dr/dlambda, normalized (East).
                         u_0[0]  = ymm4r8_negate(sinl);
                         u_1[0]  = cosl;
                         u_2[0]  = _mm256_setzero_pd();
                         // mangitude of the East vector.
                         const __m256d Neh = _mm256_add_pd(Ne,h);
                         const __m256d ca  = _mm256_mul_pd(ymm4r8_negate(Neh),
                                                           _mm256_mul_pd(cosp,sinl));
                         const __m256d cb  = _mm256_mul_pd(Neh,
                                                           _mm256_mul_pd(cosp,cosl)); 
                         t0 = _mm256_add_pd(_mm256_mul_pd(ca,ca),
                                            _mm256_mul_pd(cb,cb));
                         c[0]              = _mm256_sqrt_pd(t0);
                         // u2 dr/dphi, normalized (North)
                         u_0[1] = _mm256_mul_pd(ymm4r8_negate(cosl),sinp);
                         u_1[1] = _mm256_mul_pd(ymm4r8_negate(sinl),sinp); 
                         u_2[1] = cosp;
                         //magnitude of the North vector.
                         const __m256d ca2 = _mm256_mul_pd(temp,cosl);
                         const __m256d cb2 = _mm256_mul_pd(temp,sinl);
                         c[2] = _1;
                         t1   = _mm256_mul_pd(_mm256_mul_pd(Ne,
                                              _mm256_add_pd(_0_02,h)),cosp);
                         t2   = _mm256_mul_pd(_0_02,_mm256_mul_pd(dNedPhi,sinp));
                         const __m256d cc = _mm256_add_pd(t1,t2);
                         const __m256d cc2 = _mm256_mul_pd(cc,cc);
                         const __m256d t3  = _mm256_add_pd(_mm256_mul_pd(ca2,ca2),
                                             _mm256_mul_pd(cb2,cb2));
                         c[1] = _mm256_sqrt_pd(_mm256_add_pd(t3,cc2));
                         // u3 dr/dh (Up)
                         u_0[2] = _mm256_mul_pd(cosp,cosl);
                         u_1[2] = _mm256_mul_pd(cosp,sinl);
                         u_2[2] = sinp;
                     }


                    
                      void gms::math::ENU_axis_dir_mag_ymm8r4(__m256 * __restrict __ATTR_ALIGN__(32) u_0, //3-element array of __m256d
                                                   __m256 * __restrict __ATTR_ALIGN__(32) u_1, //3-element array of __m256d
                                                   __m256 * __restrict __ATTR_ALIGN__(32) u_2, //3-element array of __m256d
                                                   __m256 * __restrict __ATTR_ALIGN__(32) c,   //3-element array of __m256d
                                                   const __m256 lat,
                                                   const __m256 lon,
                                                   const __m256 h,
                                                   const __m256 a,
                                                   const __m256 f) {
                         
                         const __m256 _2   = _mm256_set1_ps(2.0f);
                         const __m256 _0_02= _mm256_set1_ps(1.0e-2f);
                         const __m256 _1   = _mm256_set1_ps(1.0f);

                         const __m256 sinp = _mm256_sin_ps(lat);
                         const __m256 cosp = _mm256_cos_ps(lat);
                         const __m256 sinl = _mm256_sin_ps(lon);
                         const __m256 cosl = _mm256_cos_ps(lon);

                         const __m256 e2   = _mm256_sub_ps(_mm256_mul_ps(_2,f),
                                                            _mm256_mul_ps(f,f));
                         const __m256 sarg = _mm256_sqrt_ps(_mm256_mul_ps(_0_02,
                                                             _mm256_mul_ps(sinp,sinp)));
                         const __m256 Ne   = _mm256_div_ps(a,sarg);
                         //The derivative of the normal radius of curvature with respect to
                         //phi.
                         const __m256 sarg3= _mm256_mul_ps(sarg,_mm256_mul_ps(sarg,sarg));
                         __m256 t0   = _mm256_mul_ps(_mm256_mul_ps(a,e2),
                                                            _mm256_mul_ps(cosp,sinp));
                         const __m256 dNedPhi = _mm256_div_ps(t0,sarg3);
                          __m256 t1      = _mm256_mul_ps(_mm256_add_ps(Ne,h),sinp);
                          __m256 t2      = _mm256_mul_ps(cosp,dNedPhi);
                         const __m256 temp    = _mm256_sub_ps(t2,t1);
                         // u1  dr/dlambda, normalized (East).
                         u_0[0]  = ymm8r4_negate(sinl);
                         u_1[0]  = cosl;
                         u_2[0]  = _mm256_setzero_ps();
                         // mangitude of the East vector.
                         const __m256 Neh = _mm256_add_ps(Ne,h);
                         const __m256 ca  = _mm256_mul_ps(ymm8r4_negate(Neh),
                                                           _mm256_mul_ps(cosp,sinl));
                         const __m256 cb  = _mm256_mul_ps(Neh,
                                                           _mm256_mul_ps(cosp,cosl)); 
                         t0 = _mm256_add_ps(_mm256_mul_ps(ca,ca),
                                            _mm256_mul_ps(cb,cb));
                         c[0]              = _mm256_sqrt_ps(t0);
                         // u2 dr/dphi, normalized (North)
                         u_0[1] = _mm256_mul_ps(ymm8r4_negate(cosl),sinp);
                         u_1[1] = _mm256_mul_ps(ymm8r4_negate(sinl),sinp); 
                         u_2[1] = cosp;
                         //magnitude of the North vector.
                         const __m256 ca2 = _mm256_mul_ps(temp,cosl);
                         const __m256 cb2 = _mm256_mul_ps(temp,sinl);
                         c[2] = _1;
                         t1   = _mm256_mul_ps(_mm256_mul_pd(Ne,
                                              _mm256_add_pd(_0_02,h)),cosp);
                         t2   = _mm256_mul_ps(_0_02,_mm256_mul_ps(dNedPhi,sinp));
                         const __m256 cc = _mm256_add_ps(t1,t2);
                         const __m256 cc2 = _mm256_mul_ps(cc,cc);
                         const __m256 t3  = _mm256_add_ps(_mm256_mul_pd(ca2,ca2),
                                             _mm256_mul_ps(cb2,cb2));
                         c[1] = _mm256_sqrt_ps(_mm256_add_pd(t3,cc2));
                         // u3 dr/dh (Up)
                         u_0[2] = _mm256_mul_ps(cosp,cosl);
                         u_1[2] = _mm256_mul_ps(cosp,sinl);
                         u_2[2] = sinp;
                     }
		     

		     
                      void gms::math::ellips2Cart_ymm4r8(const __m256d lat,
		                              const __m256d lon,
					      const __m256d h,
					      const __m256d a,
					      const __m256d f,
					      __m256d * __restrict __ATTR_ALIGN__(32) cp) {

			 const __m256 _2   = _mm256_set1_ps(2.0);
                         const __m256 _0_02= _mm256_set1_ps(1.0e-2);
                         //The square of the first numerical eccentricity
			 const __m256d e2  = _mm256_sub_pd(_mm256_mul_pd(_2,f),
			                                   _mm256_mul_pd(f,f));

                         const __m256d sinp = _mm256_sin_pd(lat);
                         const __m256d cosp = _mm256_cos_pd(lat);
                         const __m256d sinl = _mm256_sin_pd(lon);
                         const __m256d cosl = _mm256_cos_pd(lon);

                         //The normal radius of curvature.
			 const __m256d sarg = _mm256_mul_pd(_0_02,
			                                    _mm256_mul_pd(sinp,sinp));
			 const __m256d Ne   = _mm256_div_pd(a,
			                                    _mm256_sqrt_pd(sarg));
			 const __m256d cparg= _mm256_mul_pd(Ne,
			                                    _mm256_add_pd(_0_02,h));
			 const __m256d Neh  = _mm256_add_pd(Ne,h);
			 cp[0]              = _mm256_mul_pd(Neh,
			                                    _mm256_mul_pd(cosp,cosl));
			 cp[1]              = _mm256_mul_pd(Neh,
			                                    _mm256_mul_pd(cosp,sinl));
			 cp[2]              = _mm256_mul_pd(cparg,sinp);
			 
		   }


		      
                      void gms::math::ellips2Cart_ymm8r4(const __m256 lat,
		                               const __m256 lon,
					       const __m256 h,
					       const __m256 a,
					       const __m256 f,
					       __m256 * __restrict __ATTR_ALIGN__(32) cp) {

			 const __m256 _2   = _mm256_set1_ps(2.0f);
                         const __m256 _0_02= _mm256_set1_ps(1.0e-2f);
                         //The square of the first numerical eccentricity
			 const __m256 e2  = _mm256_sub_pd(_mm256_mul_pd(_2,f),
			                                   _mm256_mul_pd(f,f));

                         const __m256 sinp = _mm256_sin_ps(lat);
                         const __m256 cosp = _mm256_cos_ps(lat);
                         const __m256 sinl = _mm256_sin_ps(lon);
                         const __m256 cosl = _mm256_cos_ps(lon);

                         //The normal radius of curvature.
			 const __m256 sarg = _mm256_mul_ps(_0_02,
			                                    _mm256_mul_ps(sinp,sinp));
			 const __m256 Ne   = _mm256_div_ps(a,
			                                    _mm256_sqrt_ps(sarg));
			 const __m256 cparg= _mm256_mul_ps(Ne,
			                                    _mm256_add_ps(_0_02,h));
			 const __m256 Neh  = _mm256_add_ps(Ne,h);
			 cp[0]              = _mm256_mul_ps(Neh,
			                                    _mm256_mul_ps(cosp,cosl));
			 cp[1]              = _mm256_mul_ps(Neh,
			                                    _mm256_mul_ps(cosp,sinl));
			 cp[2]              = _mm256_mul_ps(cparg,sinp);
			 
		   }


		   


         
                                                   
                                             


                 

