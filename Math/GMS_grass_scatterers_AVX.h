

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
			int32_t nplants;
		        //  Number of simulation steps it is equal to Radar PRF (pulse repetetive frequency)
                        int32_t nsteps;
			//  Grass plants ordinal number (for the grass field simulation)
			int32_t ordinal;
			//   Number of parametric equation evaluation 'points' for the grass plants cylindrical approximation
			int32_t param_npts;
			//   Total grass plant area (sum of each plant surface area)
			float   tot_area;
                        //   grass plants latitude
			float   lat;
			//   grass plants longtitude
			float   lon;
			//   Elevation above the sea level
			float   elev;
			//   Apparent surface temperature (horizontal polarization)
			float   Tah;
			//   Apparent surface temperature (vertical polarization)
			float   Tav;
		        //  Complex dielectric constant per each cylinder
		        std::complex<float> epsilon;
#if (USE_STRUCT_PADDING) == 1
		      PAD_TO_ALIGNED(4,0,16)
#endif
			//  Cross-sectional area of cylinders
			float   * __restrict __ATTR_ALIGN__(8) A;
			 // Is water or rather moistness present or not on the branch surface (per n-leaves) (allowed values only 0,1)
			int32_t * __restrict __ATTR_ALIGN__(8) moistness;
			
		       
			//  ! Grass parametric equation (approximated as a cylindrical objects)
                         //  ! Parameter x, (r*cos(t))
                         //  ! PAOS type size of arrays is -- npoints/8 1st dim (evaluation of x) ,
                         //  ! nplants/8 2nd dim (number of leaves)
			AVXVec8 * __restrict __ATTR_ALIGN__(8) xparam;
			 //  ! Parameter y, (r*sin(t))
                         //  !2nd dimension is a plant number, 1st dimension evaluation of parameter y
			AVXVec8 * __restrict __ATTR_ALIGN__(8) yparam;
			 //  ! Parameter z, (height)
                         //  !2nd dimension is a plant  number, 1st dimension evaluation of parameter z
			AVXVec8 * __restrict __ATTR_ALIGN__(8) zparam;
#if (USE_STRUCT_PADDING) == 1
		      PAD_TO_ALIGNED(8,1,24)
#endif
		 } __ATTR_ALIGN__(64);


		    	//  ! This is a high termporal and spatial locality data type
                //  ! These data type members characteristics are varying between each sample of Radar PRF.
		    struct GSHotAVX_t  {

                         // Horizontal polarization (full angle sweep [2-90 deg]
			 // The 6 last elements are a 64-byte cache line padding
			 std::complex<float> Polv[96] __ATTR_ALIGN__(64);
			 // Vertical polarization (full angle sweep [2-90 deg]
			 // The 6 last elements are a 64-byte cache line padding
			 std::complex<float> Polh[96] __ATTR_ALIGN__(64);
			 
			 //   ! Grass angle of vibration in x-axis per PRF/s
                         //   ! 1st dimension angle values (rad),  2nd dimension PRF/s,
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) xang;
			 //   ! Grass sine of vibration angle in x-axis per PRF/s
                         //   ! 1st  dimension sine of vibrational angle (rad),  2nd dimension PRF/s,
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) sin_xang;
			 //     ! Grass sine of vibration angle in x-axis per PRF/s
                         //     ! 1st dimension PRF/s, 2nd dimension sine of vibrational angle (rad)
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) cos_xang;
                         //      ! Grass angle of vibration in y-axis per PRF/s
                         //      ! 1st dimension PRF/s, 2nd dimension angle values (rad)
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) yang;
			 //      ! Grass sine of vibration angle in y-axis per PRF/s
                         //      ! 1st dimension PRF/s, 2nd dimension angle of vibrational angle (rad)
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) sin_yang;
			 //        ! Grass sine of vibration angle in y-axis per PRF/s
                         //        ! 1st dimension PRF/s, 2nd dimension sine of vibrational angle (rad)
			 AVXVec8 * __restrict __ATTR_ALIGN__(8) cos_yang;
#if (USE_STRUCT_PADDING) == 1
                         PAD_TO_ALIGNED(8,0,26)
#endif
		 } __ATTR_ALIGN__(64);


	         struct GrassScattererAVX{

                        GSColdAVX_t  m_gsc __ATTR_ALIGN__(64);

			GSHotAVX_t   m_gsh __ATTR_ALIGN__(64);

			GrassScattererAVX() __ATTR_COLD__ __ATTR_ALIGN__(32);

			GrassScattererAVX(const int32_t,
			                  const int32_t,
					  const int32_t,
					  const float,
					  const float,
					  const float,
					  const std::complex<float>) __ATTR_COLD__ __ATTR_ALIGN__(32);

			GrassScattererAVX(const GrassScattererAVX &) = delete;

			GrassScattererAVX(GrassScattererAVX &&) = delete;

			~GrassScattererAVX() __ATTR_COLD__ __ATTR_ALIGN__(32) noexcept(true);

			GrassScattererAVX & operator=(const GrassScattererAVX &) = delete;

			GrassScattererAVX & operator=(GrassScattererAVX &&) = delete;

			void SetGrassMoistnessMask() __ATTR_COLD__ __ATTR_ALIGN__(32);

			void ComputeGrassParamEq_ymm8r4() __ATTR_COLD__ __ATTR_ALIGN__(32);

			void ComputeGrassHVPolarization(
						        const float,
							const float,
							const float) __ATTR_HOT__ __ATTR_ALIGN__(32);			      

			
   

	     } __ATTR_ALIGN__(64);

     } // math


} // gms










#endif /*__GMS_GRASS_SCATTERERS_AVX_H__*/
