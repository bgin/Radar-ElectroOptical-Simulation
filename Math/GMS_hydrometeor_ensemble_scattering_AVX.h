
#ifndef __GMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_H__
#define __GMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_H__


namespace file_info {

   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_MAJOR = 1;
   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_MINOR = 0;
   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_MICRO = 0;
   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_FULLVER =
       1000U*gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_MAJOR+
       100U*gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_MINOR+
       10U*gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_MICRO;
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_CREATE_DATE = "08-01-2020 15:34 +00200 (WED 08 JAN 2020 GMT+2)";
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_BUILD_DATE  = __DATE__ " " __TIME__;
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX_SYNOPSIS    = "Scattering (AVX implementation) of Radiation by Moderately Non-Spherical Particles";
}

#include <cstdint>
#include "GMS_config.h"
#include "GMS_avxvecf32.h"


namespace gms {

        namespace  math {


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

                   const AVXVec8 ONE    =   AVXVec8{1.0f};

		   const AVXVec8 TWO    =   AVXVec8{2.0f};

		   const AVXVec8 PI     =   AVXVec8{3.141592653589793f};
                    
                   const AVXVec8 TWO_PI =   AVXVec8{6.283185307179586f};


	  }

	  enum class EnsembleShapesAVX : uint32_t {

	             Cylindrical,
                     Spheroidal,
		     Chebyshev
	  };


	   // Low spatial and temporal frequency data type
	  struct HSColdAVX_t {

                //     ! Number of particles  in aggregated emsemble
		int32_t m_np;
		int32_t m_nshpts;
		//     ! Particles aggregate ID number
		int32_t m_ID;
		//      ! Time evolution steps
                int32_t m_nt;
		//       ! Maximal number of parametric equation points
		int32_t m_nxpts;
		int32_t m_nypts;
		int32_t m_nzpts;
		//    ! Total volume of particles per ensemble
		float   m_tpv;
		//     ! Total particles surface area per ensemble
		float   m_tpsa;
		//     Total particles mass per ensemble
		float   m_tpm;
		//   ! Chebyshev particles shape in aggregated assembly ( (r = r0[1+eTn(cos(theta))])
                //   ! [r,np], where r = parametric shape (cross-section), np =  n-th particle
		AVXVec8 * __restrict __ATTR_ALIGN__(8) m_pcs;
		//    ! Chebyshev particles radii in aggregate ensemble
		AVXVec8 * __restrict __ATTR_ALIGN__(8) m_radii;
		//       ! Chebyshev particles aggregate shape approximated by 3D parametric equations.
                //       ! Components form location of the particle in the ensemble.
                //       ! [3,np], where first dimension represents coordinate components
                //       ! second dimension represent number of particles.
	        float   * __restrict __ATTR_ALIGN__(8) m_pes;
		//        ! Chebyshev particles ensemble( per each particle) parametric equation in x - dimension (non-dimensional)
                //        ! [paramx,np]
		AVXVec8 * __restrict __ATTR_ALIGN__(8) m_ppx;
		//         ! Chebyshev particles ensemble (per each particle)  parametric equation in y  - dimension (non-dimensional)
		AVXVec8 * __restrict __ATTR_ALIGN__(8) m_ppy;
		//        ! Chebyshev particles ensemble (per each particle)  parametric equation in z  - dimension (non-dimensional)
		AVXVec8 * __restrict __ATTR_ALIGN__(8) m_ppz;
		//        ! Hydrometeor type
	        char * __restrict __ATTR_ALIGN__(8) m_type;
		//         ! Ensemble shape  only 2 types supported for now: -- ( spheroid, chebyshev-particle)
	        char * __restrict __ATTR_ALIGN__(8) m_shape;
	   } __ATTR_ALIGN__(64);


	   // ! High temporal and spatial frequency data type.

	   struct  HSHotAVX_t {

               //  ! Yu-lin Xu part of the variables
              
	       double m_dang[1856]    __ATTR_ALIGN__(64);
	       double m_imat[1856]    __ATTR_ALIGN__(64);
	       double m_pol[1856]     __ATTR_ALIGN__(64);
	       double m_i11[1856]     __ATTR_ALIGN__(64);
	       double m_i21[1856]     __ATTR_ALIGN__(64);
	       double m_i12[1856]     __ATTR_ALIGN__(64);
	       double m_i22[1856]     __ATTR_ALIGN__(64);
	       double m_mue[4*4*1856] __ATTR_ALIGN__(64);
	       double m_cext;
               double m_cabs;
               double m_csca;
	       double m_assym;
	       double m_cextv;
	       double m_cabsv;
	       double m_cscav;
	       double m_cbakv;
	       double m_cprv;
	       double m_cexts;
	       double m_cabss;
	       double m_cscas;
	       double m_cbaks;
	       double m_cprs;
	       double * __restrict __ATTR_ALIGN__(8) m_cexti;
	       double * __restrict __ATTR_ALIGN__(8) m_cabsi;
	       double * __restrict __ATTR_ALIGN__(8) m_cscai;
	       double * __restrict __ATTR_ALIGN__(8) m_assymi;
	       double * __restrict __ATTR_ALIGN__(8) m_cpri;
	       //            ! Trajectory of Chebyshev particles ensemble, radial distance component (spherical coordinate system)
               //! [nt]
               float * __restrict __ATTR_ALIGN__(8) m_rdist;
	       //  ! Trajectory of Chebyshev particles ensemble, theta angle component (spherical coordinate system)
               //  ! [nt]
	       float * __restrict __ATTR_ALIGN__(8) m_theta;
	       //  ! Trajectory of Chebyshev particles ensemble, phi angle component (spherical coordinate system)
               //  ! [nt]
	       float * __restrict __ATTR_ALIGN__(8) m_phi;
	       //   ! Chebyshev particles ensemble fall speed 
               //   ![nt]
	       float * __restrict __ATTR_ALIGN__(8) m_vfall;
	  } __ATTR_ALIGN__(64);


	   struct  HMScatterersAVX {

	           HSColdAVX_t m_hsc; //cold

		   HSHotAVX_t  m_hsh; // hot

		   HMScatterersAVX() __ATTR_COLD__ __ATTR_ALIGN__(32);

		   HMScatterersAVX( const int32_t,
		                          const int32_t,
		                          const int32_t,
					  const int32_t,
					  const int32_t,
					  const int32_t,
					  const int32_t,
					  const char *,
					  const char *) __ATTR_COLD__ __ATTR_ALIGN__(32);

		   HMScatterersAVX(const HMScatterersAVX &) = delete;

		   HMScatterersAVX(HMScatterersAVX &&) = delete;

		   ~HMScatterersAVX() noexcept(true);

		   HMScatterersAVX & operator=(const HMScatterersAVX &) = delete;

		   HMScatterersAVX & operator=(HMScatterersAVX &&) = delete;

		   bool ComputeShape_ymm8r4(
#if defined __ICC || defined __INTEL_COMPILER
                                            AVXVec8 * __restrict __ATTR_ALIGN__(64),
		                            AVXVec8 * __restrict __ATTR_ALIGN__(64)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                             float * __restrict __ATTR_ALIGN__(64),
                                             float * __restrict __ATTR_ALIGN__(64)
#endif
)  __ATTR_COLD__ __ATTR_ALIGN__(32);

				         
#if defined __ICC || defined __INTEL_COMPILER
    #if defined __AVX__
	     #pragma intel optimization_parameter target_arch=AVX
    #elif defined __AVX512F__
	     #pragma intel optimization_parameter target_arch=AVX512
    #endif
#endif
		   bool ComputeEnsembleShape(const float,
		                             const float,
					     const EnsembleShapes,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float) __ATTR_COLD__ __ATTR_ALIGN__(64) __ATTR_TCLONES_AVX_AVX512__;

		   bool ComputeXparam_ymm8r4(const AVXVec8 * __restrict __ATTR_ALIGN__(64),
		                             const AVXVec8 * __restrict __ATTR_ALIGN__(64)) __ATTR_COLD__ __ATTR_ALIGN__(32);
					     

		   bool ComputeYparam_ymm8r4(const AVXVec8 * __restrict __ATTR_ALIGN__(64),
		                             const AVXVec8 * __restrict __ATTR_ALIGN__(64)) __ATTR_COLD__ __ATTR_ALIGN__(32);
					    

		   bool ComputeZparam_ymm8r4(const AVXVec8 * __restrict __ATTR_ALIGN__(64),
		                             const AVXVec8 * __restrict __ATTR_ALIGN__(64)) __ATTR_COLD__ __ATTR_ALIGN__(32);
					     

		   bool ComputeEnsembleVolume(const float * __restrict __ATTR_ALIGN__(64),
		                              const float * __restrict __ATTR_ALIGN__(64),
					      const float * __restrict __ATTR_ALIGN__(64)
					      const int32_t ) __ATTR_COLD__ __ATTR_ALIGN__(32) __ATTR_TCLONES_AVX_AVX512__;

		   bool ComputeEnsembleSurface(const float * __restrict __ATTR_ALIGN__(64),
		                               const float * __restrict __ATTR_ALIGN__(64),
					       const float * __restrict __ATTR_ALIGN__(64),
					       const int32_t) __ATTR_COLD__ __ATTR_ALIGN__(32) __ATTR_TCLONES_AVX_AVX512__;

#if defined __ICC || defined __INTEL_COMPILER
    #if defined __AVX__
	     #pragma intel optimization_parameter target_arch=AVX
    #elif defined __AVX512F__
	     #pragma intel optimization_parameter target_arch=AVX512
    #endif
#endif		   
		   void ComputeEnsembleVfall(  const float * __restrict __ATTR_ALIGN__(64),
			               const float * __restrict __ATTR_ALIGN__(64),
				       const float,
				       const float * __restrict __ATTR_ALIGN__(64),
				       const int32_t,
				       const int32_t,
				       const int32_t,
				       const float,
				       const float,
				       const float * __restrict __ATTR_ALIGN__(64),
				       const float,
				       const float,
				       const float * __restrict __ATTR_ALIGN__(64),
				       const float * __restrict __ATTR_ALIGN__(64)) __ATTR_HOT__ __ATTR_ALIGN__(32) __ATTR_TCLONES_AVX_AVX512__;
			                     
	           void ComputeHydroMeteorScattering(
		                                     int32_t,
						     int32_t,
						     int32_t,
						     int32_t,
						     int32_t,
						     int32_t,
						     int32_t,
						     int32_t,
						     int32_t,
						     int32_t,
						     int32_t,
						     double,
						     int32_t,
						     int32_t,
						     int32_t,
						     double,
						     double,
						     int32_t,
						     int32_t,
						     int32_t * __restrict __ATTR_ALIGN__(64),
						     double  * __restrict __ATTR_ALIGN__(64),
						     double  * __restrict __ATTR_ALIGN__(64) __ATTR_HOT__ __ATTR_ALIGN__(32);
		   


	   } __ATTR_ALIGN__(64);

	  

    }  // math


} // gms




#endif /*__GMS_HYDROMETEOR_ENSEMBLE_SCATTERING_H__*/
