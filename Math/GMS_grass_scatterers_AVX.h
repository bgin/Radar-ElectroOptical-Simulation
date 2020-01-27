

#ifndef __GMS_GRASS_SCATTERERS_AVX_H__
#define __GMS_GRASS_SCATTERERS_AVX_H__


namespace file_info {

  const unsigned int gGMS_GRASS_SCATTERERS_AVX_MAJOR = 1U;
  const unsigned int gGMS_GRASS_SCATTERERS_AVX_MINOR = 0U;
  const unsigned int gGMS_GRASS_SCATTERERS_AVX_MICRO = 0U;
  const unsigned int gGMS_GRASS_SCATTERERS_AVX_FULLVER = 1000U*gGMS_GRASS_SCATTERERS_AVX_MAJOR +
                                                         100U*gGMS_GRASS_SCATTERERS_AVX_MINOR  +
                                                         10U*gGMS_GRASS_SCATTERERS_AVX_MICRO;
  const char * const pgGMS_GRASS_SCATTERER_AVX_CREATE_DATE = "27-01-2020 15:00 +00200 (MON 27 JAN 2020 GMT+2)";
  const char * const pgGMS_GRASS_SCATTERER_AVX_BUILD_DATE  = __DATE__ " " __TIME__;
  const char * const pgGMS_GRASS_SCATTERER_AVX_AUTHOR      = " Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_GRASS_SCATTERER_AVX_SYNOPSYS    = " Model of grass scatterers(AVX implementation) suitable of computation of radar backscatter.";
  
}

#include <cstdint>
#include <complex>
#include "GMS_config.h"
#include "GMS_avxvecf32.h"


namespace gms {

 
          namespace math {

	             namespace {

                          const AVXVec8 VINC0  =   AVXVec8{1.0f,2.0f,3.0f,4.0f,
		                                  5.0f,6.0f,7.0f,8.0f};
		          const AVXVec8 VINC1  =   AVXVec8{9.0f,10.0f,11.0f,12.0f,
		                                  13.0f,14.0f,15.0f,16.0f};
		          const AVXVec8 VINC2  =   AVXVec8{17.0f,18.0f,19.0f,20.0f,
		                                   21.0f,22.0f,23.0f,24.0f};
		          const AVXVec8 VINC3  =   AVXVec8{25.0f,26.0f,27.0f,28.0f,
		                                   29.0f,30.0f,31.0f,32.0f};
                          const AVXVec8 ZERO   =   AVXVec8{};

		          const AVXVec8 TWO    =   AVXVec8{2.0f};

		          const AVXVec8 PI     =   AVXVec8{3.141592653589793f};
	   }

	     // Low temporal access and spatial locality (cold)  data structure (POD) type
	            struct GSColdAVX_t {

                        // Number of grass plants per unit area (1 m^2)
			int32_t m_nplants;
		        //  Number of simulation steps it is equal to Radar PRF (pulse repetetive frequency)
                        int32_t m_nsteps;
			//  Grass plants ordinal number (for the grass field simulation)
			int32_t m_ordinal;
			//   Number of parametric equation evaluation 'points' for the grass plants cylindrical approximation
			int32_t m_grass_param_npts;
			//   Total grass plant area (sum of each plant surface area)
			float   m_tot_area;
                        //   grass plants latitude
			float   m_lat;
			//   grass plants longtitude
			float   m_lon;
			//   Elevation above the sea level
			float   m_elev;
			//   Apparent surface temperature (horizontal polarization)
			float   m_Tah;
			//   Apparent surface temperature (vertical polarization)
			float   m_Tav;
			//  Cross-sectional area of cylinders
			float   * __restrict __ATTR_ALIGN__(8) m_A;
			 // Is water or rather moistness present or not on the branch surface (per n-leaves) (allowed values only 0,1)
			int32_t * __restrict __ATTR_ALIGN__(8) m_moistness;
			//  Complex dielectric constant per each cylinder
			std::complex<float> * __restrict __ATTR_ALIGN__(8) m_eps;
			//  ! Grass parametric equation (approximated as a cylindrical objects)
                         //  ! Parameter x, (r*cos(t))
                         //  ! PAOS type size of arrays is -- npoints/8 1st dim (evaluation of x) ,
                         //  ! nplants/8 2nd dim (number of leaves)
			AVXVec8 * __restrict __ATTR_ALIGN__(8) grass_xparam;
			 //  ! Parameter y, (r*sin(t))
                         //  !2nd dimension is a plant number, 1st dimension evaluation of parameter y
			AVXVec8 * __restrict __ATTR_ALIGN__(8) grass_yparam;
			 //  ! Parameter z, (height)
                         //  !2nd dimension is a plant  number, 1st dimension evaluation of parameter z
			AVXVec8 * __restrict __ATTR_ALIGN__(8) grass_zparam;
		 } __ATTR_ALIGN__(64);

     } // math


} // gms










#endif /*__GMS_GRASS_SCATTERERS_AVX_H__*/
