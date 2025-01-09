
#ifndef __GMS_GEODESY_YMM4R8_H__
#define __GMS_GEODESY_YMM4R8_H__ 171020211522



namespace file_version {

    const unsigned int GMS_GEODESY_YMM4R8_MAJOR = 1U;
    const unsigned int GMS_GEODESY_YMM4R8_MINOR = 0U;
    const unsigned int GMS_GEODESY_YMM4R8_MICRO = 0U;
    const unsigned int GMS_GEODESY_YMM4R8_FULLVER =
      1000U*GMS_GEODESY_YMM4R8_MAJOR+
      100U*GMS_GEODESY_YMM4R8_MINOR+
      10U*GMS_GEODESY_YMM4R8_MICRO;
    const char * const GMS_GEODESY_YMM4R8_CREATION_DATE = "17-10-2021 15:22  +00200 (SUN 17 OCT 2021 GMT+2)";
    const char * const GMS_GEODESY_YMM4R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_GEODESY_YMM4R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_GEODESY_YMM4R8_DESCRIPTION   = "Vectorized (AVX/AVX2) geodesic computation implementation."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {

          namespace math {


	             /*
                              Cartesian to geodetic conversion (kernel).
                              
                          */
	            
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			void
			cart_to_geodetic_ymm4r8(     const __m256d pos_x, //input position x [km]
			                             const __m256d pos_y, //input position y [km]
						     const __m256d pos_z, //input position z [km]
						     const __m256d a, // semi-minor axis [km]
						     const __m256d b, // semi-major axis [km]
						     __m256d &alt, //output altitude [km]
						     __m256d &lon, //output longtitude [rad]
						     __m256d &lat); 
						     

			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			void
			cart_to_geodetic_u_ymm4r8_looped(const double * __restrict pos_x,
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
			cart_to_geodetic_a_ymm4r8_looped(double * __restrict __ATTR_ALIGN__(32) pos_x,
			                                 double * __restrict __ATTR_ALIGN__(32) pos_y,
							 double * __restrict __ATTR_ALIGN__(32) pos_z,
							 const double a,
							 const double b,
							 double * __restrict __ATTR_ALIGN__(32) alt,
							 double * __restrict __ATTR_ALIGN__(32) lon,
							 double * __restrict __ATTR_ALIGN__(32) lat,
							 const int32_t n); 

	    

		      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			void geodetic_to_cart_ymm4r8(const __m256d a,
			                             const __m256d b,
						     const __m256d lat,
						     const __m256d lon,
						     const __m256d alt,
						     __m256d &pos_x,
						     __m256d &pos_y,
						     __m256d &pos_z); 



			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			void
			geodetic_to_cart_u_ymm4r8_looped(const double a,
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
			geodetic_to_cart_a_ymm4r8_looped(const double a,
			                                 const double b,
							 const double * __restrict  __ATTR_ALIGN__(32) lat,
							 const double * __restrict  __ATTR_ALIGN__(32) lon,
							 const double * __restrict  __ATTR_ALIGN__(32) alt,
							 double * __restrict  __ATTR_ALIGN__(32) pos_x,
							 double * __restrict  __ATTR_ALIGN__(32) pos_y,
							 double * __restrict  __ATTR_ALIGN__(32) pos_z,
							 const int32_t n); 


			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                void forward_method_ymm4r8(const __m256d axis,      //ellipsoid semi-maxjor axis
			                           const __m256d flat,      //elipsoid flattening [dimensionless]
						   const __m256d vp1lat,    //vector of 4 starting-points latitude [rad]
						   const __m256d vp1lon,    //vector of 4 starting-points longtitude [rad]
						   const __m256d azvf,      //vector of 4 forward azimutes [rad]
						   const __m256d dstv,      //vector of 4 distances vp1-to-vp2 [m]
						   __m256d &vp2lat,         //vector of 4 endpoints latitude [rad]
						   __m256d &vp2lon,         //vector of 4 endpoints longtitude [rad]
						   __m256d &azvb); 



			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
		       
			void 
                        forward_method_u_ymm4r8_looped(const double axis,
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
                        forward_method_a_ymm4r8_looped(const double axis,
			                               const double flat,
						       double * __restrict __ATTR_ALIGN__(32) plat1,
						       double * __restrict __ATTR_ALIGN__(32) plon1,
						       double * __restrict __ATTR_ALIGN__(32) pfaz,
						       double * __restrict __ATTR_ALIGN__(32) pdst,
						       double * __restrict __ATTR_ALIGN__(32) plat2,
						       double * __restrict __ATTR_ALIGN__(32) plon2,
						       double * __restrict __ATTR_ALIGN__(32) pbaz,
						       const int32_t n); 


			 /*
                           @Reference:
                                          http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
                        */
		      	
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__m256d
			spheroid_distance_ymm4r8(const __m256d vr,
			                         const __m256d vlon1,
						 const __m256d vlat1,
						 const __m256d vlon2,
						 const __m256d vlat2); 


			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			spheroid_distance_u_ymm4r8_looped(const double r,
			                                  double * __restrict plon1,
							  double * __restrict plat1,
							  double * __restrict plon2,
							  double * __restrict plat2,
							  double * __restrict pd,
							  const int32_t n); 


			
                        
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			spheroid_distance_a_ymm4r8_looped(const double r,
			                                  double * __restrict __ATTR_ALIGN__(32) plon1,
							  double * __restrict __ATTR_ALIGN__(32) plat1,
							  double * __restrict __ATTR_ALIGN__(32) plon2,
							  double * __restrict __ATTR_ALIGN__(32) plat2,
							  double * __restrict __ATTR_ALIGN__(32) pd,
							  const int32_t n); 

			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			
			__m256d
			geocentric_radius_ymm4r8(const __m256d va,
			                         const __m256d vb,
						 const __m256d vlat); 



		     
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			geocentric_radius_u_ymm4r8_looped(const double a,
			                                  const double * __restrict pb,
							  const double * __restrict plat,
							  double * __restrict pr,
							  const int32_t n);


			
			__ATTR_HOT__
                        __ATTR_ALIGN__(32)
                        void
			geocentric_radius_a_ymm4r8_looped(const double a,
			                                  double * __restrict __ATTR_ALIGN__(32) pb,
							  double * __restrict __ATTR_ALIGN__(32) plat,
							  double * __restrict __ATTR_ALIGN__(32) pr,
							  const int32_t n); 





			











     }

}











#endif /*__GMS_GEODESY_YMM4R8_H__*/
