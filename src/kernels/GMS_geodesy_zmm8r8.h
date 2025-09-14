
#ifndef __GMS_GEODESY_AVX512_H__
#define __GMS_GEODESY_AVX512_H__ 161020210952


namespace file_version {

    const unsigned int gGMS_GEODESY_AVX512_MAJOR = 1U;
    const unsigned int gGMS_GEODESY_AVX512_MINOR = 0U;
    const unsigned int gGMS_GEODESY_AVX512_MICRO = 0U;
    const unsigned int gGMS_GEODESY_AVX512_FULLVER =
      1000U*gGMS_GEODESY_AVX512_MAJOR+
      100U*gGMS_GEODESY_AVX512_MINOR+
      10U*gGMS_GEODESY_AVX512_MICRO;
    const char * const pgGMS_GEODESY_AVX512_CREATION_DATE = "16-10-2021 09:52 AM +00200 (SAT 16 OCT 2021 GMT+2)";
    const char * const pgGMS_GEODESY_AVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_GEODESY_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_GEODESY_AVX512_DESCRIPTION   = "Vectorized (AVX512) geodesic computation implementation."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace  gms {

          namespace math {


	                      namespace {

                                    __ATTR_ALWAYS_INLINE__
                                     static inline
				    __m512d
				    zmm8r8_sign_zmm8r8(const __m512d va,
				                       const __m512d vb) {
				       const register _0 = _mm512_setzero_pd();
				       register __m512d vret = _0;
				       register __m512d t0   = _mm512_abs_pd(va);
                                       __mmask8 gez = 0x0;
				       gez  = _mm512_cmp_pd_mask(vb,_0,_CMP_GE_OQ); // Lat=3refc,Thr=1refc
				       vret = _mm512_mask_blend_pd(gez,t0,_mm512_sub_pd(_0,t0)); //Lat=1refc,Thr=0.5refc,Lat=4refc,Thr=1refc
				       return (vret);
				                                       
				    }
			    }

                          
	                 /*
                              Cartesian to geodetic conversion (kernel).
                              
                          */
	               
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	                void
			cart_to_geodetic_zmm8r8(     const __m512d pos_x, //input position x [km]
			                             const __m512d pos_y, //input position y [km]
						     const __m512d pos_z, //input position z [km]
						     const __m512d a, // semi-minor axis [km]
						     const __m512d b, // semi-major axis [km]
						     __m512d &alt, //output altitude [km]
						     __m512d &lon, //output longtitude [rad]
						     __m512d &lat); 

			
			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			void
			cart_to_geodetic_u_zmm8r8_looped(const double * __restrict pos_x,
			                                 const double * __restrict pos_y,
							 const double * __restrict pos_z,
							 const double a,
							 const double b,
							 double * __restrict alt,
							 double * __restrict lon,
							 double * __restrict lat,
							 const int32_t n); 


			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			void
			cart_to_geodetic_a_zmm8r8_looped(double * __restrict __ATTR_ALIGN__(64) pos_x,
			                                 double * __restrict __ATTR_ALIGN__(64) pos_y,
							 double * __restrict __ATTR_ALIGN__(64) pos_z,
							 const double a,
							 const double b,
							 double * __restrict __ATTR_ALIGN__(64) alt,
							 double * __restrict __ATTR_ALIGN__(64) lon,
							 double * __restrict __ATTR_ALIGN__(64) lat,
							 const int32_t n); 

                    	
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	               	void geodetic_to_cart_zmm8r8(const __m512d a,
			                             const __m512d b,
						     const __m512d lat,
						     const __m512d lon,
						     const __m512d alt,
						     __m512d &pos_x,
						     __m512d &pos_y,
						     __m512d &pos_z);

                    	
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			void
			geodetic_to_cart_u_zmm8r8_looped(const double a,
			                                 const double b,
							 const double * __restrict lat,
							 const double * __restrict lon,
							 const double * __restrict alt,
							 double * __restrict pos_x,
							 double * __restrict pos_y,
							 double * __restrict pos_z,
							 const int32_t n); 

			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
		        void
			geodetic_to_cart_a_zmm8r8_looped(const double a,
			                                 const double b,
							 double * __restrict __ATTR_ALIGN__(64) lat,
							 double * __restrict __ATTR_ALIGN__(64) lon,
							 double * __restrict __ATTR_ALIGN__(64) alt,
							 double * __restrict __ATTR_ALIGN__(64) pos_x,
							 double * __restrict __ATTR_ALIGN__(64) pos_y,
							 double * __restrict __ATTR_ALIGN__(64) pos_z,
							 const int32_t n); 
	              

                        /*
                                Reference: https://www.ngs.noaa.gov/PC_PROD/Inv_Fwd/
                           */
			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	              	void forward_metod_zmm8r8(const __m512d axis,      //ellipsoid semi-maxjor axis
			                         const __m512d flat,      //elipsoid flattening [dimensionless]
						 const __m512d vp1lat,    //vector of 8 starting-points latitude [rad]
						 const __m512d vp1lon,    //vector of 8 starting-points longtitude [rad]
						 const __m512d azvf,      //vector of 8 forward azimutes [rad]
						 const __m512d dstv,      //vector of 8 distances vp1-to-vp2 [m]
						 __m512d &vp2lat,         //vector of 8 endpoints latitude [rad]
						 __m512d &vp2lon,         //vector of 8 endpoints longtitude [rad]
						 __m512d &azvb); 


		     
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
		       	void
			forward_method_u_zmm8r8_looped(const double axis,
			                               const double flat,
						       const double * __restrict plat1,
						       const double * __restrict plon1,
						       const double * __restrict pfaz,
						       const double * __restrict pdst,
						       double * __restrict plat2,
						       double * __restrict plon2,
						       double * __restrict pbaz,
						       const int32_t n); 



		   
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
		       void
			forward_method_a_zmm8r8_looped(const double axis,
			                               const double flat,
						       double * __restrict __ATTR_ALIGN__(64) plat1,
						       double * __restrict __ATTR_ALIGN__(64) plon1,
						       double * __restrict __ATTR_ALIGN__(64) pfaz,
						       double * __restrict __ATTR_ALIGN__(64) pdst,
						       double * __restrict __ATTR_ALIGN__(64) plat2,
						       double * __restrict __ATTR_ALIGN__(64) plon2,
						       double * __restrict __ATTR_ALIGN__(64) pbaz,
						       const int32_t n); 
		      /*
                           @Reference:
                                          http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
                        */
		      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	              	__m512d
			spheroid_distance_zmm8r8(const __m512d vr,
			                         const __m512d vlon1,
						 const __m512d vlat1,
						 const __m512d vlon2,
						 const __m512d vlat2); 


			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			spheroid_distance_u_zmm8r8_looped(const double r,
			                                  double * __restrict plon1,
							  double * __restrict plat1,
							  double * __restrict plon2,
							  double * __restrict plat2,
							  double * __restrict pd,
							  const int32_t n); 

			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			spheroid_distance_a_zmm8r8_looped(const double r,
			                                  double * __restrict __ATTR_ALIGN__(64) plon1,
							  double * __restrict __ATTR_ALIGN__(64) plat1,
							  double * __restrict __ATTR_ALIGN__(64) plon2,
							  double * __restrict __ATTR_ALIGN__(64) plat2,
							  double * __restrict __ATTR_ALIGN__(64) pd,
							  const int32_t n); 

			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	               	__m512d
			geocentric_radius_zmm8r8(const __m512d va,
			                         const __m512d vb,
						 const __m512d vlat); 


		       
		       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			geocentric_radius_u_zmm8r8_looped(const double a,
			                                  const double * __restrict pb,
							  const double * __restrict plat,
							  double * __restrict pr,
							  const int32_t n); 

                      
			__ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			geocentric_radius_a_zmm8r8_looped(const double a,
			                                  const double * __restrict __ATTR_ALIGN__(64) pb,
							  const double * __restrict __ATTR_ALIGN__(64) plat,
							  double * __restrict __ATTR_ALIGN__(64) pr,
							  const int32_t n); 
                        /*
                              Based on: http://www.ngs.noaa.gov/PC_PROD/Inv_Fwd/source/inverse.for
                          */

                       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	              	void inverse_method_zmm8r8(const __m512d va,    // Semi-major axis (equatorial)
			                           const __m512d vrf,   // reciprocal flattening
						   const __m512d vlat1, // Latitude of 8 points [rad, positive north]
						   const __m512d vlon1, // longtitude of 8 points [rad,positive east]
						   const __m512d vlat2, // Latitude of 8 points [rad, positive north]
						   const __m512d vlon2, // Longtitude of 8 points [rad, positive east]
						   __m512d &vfaz,       // Vector of 8 forward azimuths [rad]
						   __m512d &vbaz,       // Vector of 8 backward azimuthes [rad]
						   __m512d &vs,         // Ellipsoidal distance
						   int32_t &icnt,     // iteration count
						   __m512d &sig,       // Spherical distance (auxiliary sphere)
						   __m512d &vld,        // Longtitude difference (auxiliary sphere)
						   int32_t &kind); 




						     
     }

}














#endif /*__GMS_GEODESY_ZMM8R8_H__*/
