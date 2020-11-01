

#ifndef __GMS_GRASS_SCATTERERS_AVX512_H__
#define __GMS_GRASS_SCATTERERS_AVX512_H__


namespace file_info {

  const unsigned int gGMS_GRASS_SCATTERERS_AVX512_MAJOR = 1U;
  const unsigned int gGMS_GRASS_SCATTERERS_AVX512_MINOR = 1U;
  const unsigned int gGMS_GRASS_SCATTERERS_AVX512_MICRO = 0U;
  const unsigned int gGMS_GRASS_SCATTERERS_AVX512_FULLVER = 1000U*gGMS_GRASS_SCATTERERS_AVX512_MAJOR +
                                                         100U*gGMS_GRASS_SCATTERERS_AVX512_MINOR  +
                                                         10U*gGMS_GRASS_SCATTERERS_AVX512_MICRO;
  const char * const pgGMS_GRASS_SCATTERER_AVX512_CREATE_DATE = "30-01-2020 15:28 +00200 (MON 27 JAN 2020 GMT+2)";
  const char * const pgGMS_GRASS_SCATTERER_AVX512_BUILD_DATE  = __DATE__ " " __TIME__;
  const char * const pgGMS_GRASS_SCATTERER_AVX512_AUTHOR      = " Programmer: Bernard Gingold, contact: beniekg@gmail.com";
  const char * const pgGMS_GRASS_SCATTERER_AVX512_SYNOPSYS    = " Model of grass scatterers(AVX512 implementation) suitable of computation of radar backscatter.";
  
}


#include <cstdint>
#include <complex>
#include "GMS_config.h"
#include "GMS_avx512vecf32.h"

namespace gms {

           namespace math {


                    namespace {

                                 const AVX512Vec16 VINC0 = AVX512Vec16{1.0f,2.0f,3.0f,4.0f,
				                                       5.0f,6.0f,7.0f,8.0f,
								       9.0f,10.0f,11.0f,12.0f,
								       13.0f,14.0f,15.0f,16.0f};

				 const AVX512Vec16 VINC1 = AVX512Vec16{17.0f,18.0f,19.0f,20.0f,
				                                       21.0f,22.0f,23.0f,24.0f,
								       25.0f,26.0f,27.0f,28.0f,
								       29.0f,30.0f,31.0f,32.0f};

				 const AVX512Vec16 VINC2 = AVX512Vec16{33.0f,34.0f,35.0f,36.0f,
				                                       37.0f,38.0f,39.0f,40.0f,
								       41.0f,42.0f,43.0f,44.0f,
								       45.0f,46.0f,47.0f,48.0f};

				 const AVX512Vec16 VINC3 = AVX512Vec16{49.0f,50.0f,51.0f,52.0f,
				                                       53.0f,54.0f,55.0f,56.0f,
								       57.0f,58.0f,59.0f,60.0f,
								       61.0f,62.0f,63.0f,64.0f};

				 const AVX512Vec16 ZERO  = AVX512Vec16{0.0f};

				 const AVX512Vec16 n2    = AVX512Vec16{2.0f};

				 const AVX512Vec16 nPI   = AVX512Vec16{3.141592653589793f};
								       
	           } // anon


		  // Low temporal access and spatial locality (cold)  data structure (POD) type
		    struct GSColdAVX512_t {
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
		        // bool: use mlock() function
			bool    mem_lock;
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
			AVX512Vec16 * __restrict __ATTR_ALIGN__(8) xparam;
			 //  ! Parameter y, (r*sin(t))
                         //  !2nd dimension is a plant number, 1st dimension evaluation of parameter y
			AVX512Vec16 * __restrict __ATTR_ALIGN__(8) yparam;
			 //  ! Parameter z, (height)
                         //  !2nd dimension is a plant  number, 1st dimension evaluation of parameter z
			AVX512Vec16 * __restrict __ATTR_ALIGN__(8) zparam;
#if (USE_STRUCT_PADDING) == 1
		      PAD_TO_ALIGNED(8,1,24)
#endif
                 } __ATTR_ALIGN__(64);

                //  ! This is a high termporal and spatial locality data type
                //  ! These data type members characteristics are varying between each sample of Radar PRF.
                      struct GSHotAVX512_t {

                         // Horizontal polarization (full angle sweep [2-90 deg]
			 // The 6 last elements are a 64-byte cache line padding
			 std::complex<float> Polv[96] __ATTR_ALIGN__(64);
			 // Vertical polarization (full angle sweep [2-90 deg]
			 // The 6 last elements are a 64-byte cache line padding
			 std::complex<float> Polh[96] __ATTR_ALIGN__(64);
			 
			 //   ! Grass angle of vibration in x-axis per PRF/s
                         //   ! 1st dimension angle values (rad),  2nd dimension PRF/s,
			 AVX512Vec16 * __restrict __ATTR_ALIGN__(8) xang;
			 //   ! Grass sine of vibration angle in x-axis per PRF/s
                         //   ! 1st  dimension sine of vibrational angle (rad),  2nd dimension PRF/s,
			 AVX512Vec16 * __restrict __ATTR_ALIGN__(8) sin_xang;
			 //     ! Grass sine of vibration angle in x-axis per PRF/s
                         //     ! 1st dimension PRF/s, 2nd dimension sine of vibrational angle (rad)
			 AVX512Vec16 * __restrict __ATTR_ALIGN__(8) cos_xang;
                         //      ! Grass angle of vibration in y-axis per PRF/s
                         //      ! 1st dimension PRF/s, 2nd dimension angle values (rad)
			 AVX512Vec16 * __restrict __ATTR_ALIGN__(8) yang;
			 //      ! Grass sine of vibration angle in y-axis per PRF/s
                         //      ! 1st dimension PRF/s, 2nd dimension angle of vibrational angle (rad)
			 AVX512Vec16 * __restrict __ATTR_ALIGN__(8) sin_yang;
			 //        ! Grass sine of vibration angle in y-axis per PRF/s
                         //        ! 1st dimension PRF/s, 2nd dimension sine of vibrational angle (rad)
			 AVX512Vec16 * __restrict __ATTR_ALIGN__(8) cos_yang;
#if (USE_STRUCT_PADDING) == 1
                         PAD_TO_ALIGNED(8,0,26)
#endif

		 } __ATTR_ALIGN__(64);

	             struct GrassScattererAVX512 {

                        GSColdAVX512_t  m_gsc __ATTR_ALIGN__(64);

			GSHotAVX512_t   m_gsh __ATTR_ALIGN__(64);

			GrassScattererAVX512() __ATTR_COLD__ __ATTR_ALIGN__(32);

			GrassScattererAVX512(const int32_t,
			                     const int32_t,
					     const int32_t,
					     const float,
					     const float,
					     const float,
					     const std::complex<float>,
					     const bool )__ATTR_COLD__ __ATTR_ALIGN__(32);

			GrassScattererAVX512(const GrassScattererAVX512 &) = delete;

			GrassScattererAVX512(GrassScattererAVX512 &&) = delete;

			~GrassScattererAVX512() __ATTR_COLD__ __ATTR_ALIGN__(32) noexcept(true);

			GrassScattererAVX512 & operator=(const GrassScattererAVX512 &) = delete;

			GrassScattererAVX512 & operator=(GrassScattererAVX512 &&) = delete;

			void SetGrassMoistnessMask() __ATTR_COLD__ __ATTR_ALIGN__(32);

		        void ComputeGrassParamEq_zmm16r4() __ATTR_COLD__ __ATTR_ALIGN__(32);

			void ComputeGrassHVPolarization(
						        const float,
							const float,
							const float) __ATTR_HOT__ __ATTR_ALIGN__(32);

                 } __ATTR_ALIGN__(64);


		 
 
		     


    } // math





} // gms









#endif /*__GMS_GRASS_SCATTERERS_AVX512_H__*/
