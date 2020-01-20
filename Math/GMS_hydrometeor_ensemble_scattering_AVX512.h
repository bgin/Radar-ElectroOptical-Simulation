

#ifndef __GMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_H__
#define __GMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_H__

namespace file_info {

   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_MAJOR = 1U;
   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_MINOR = 0U;
   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_MICRO = 0U;
   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_FULLVER =
       1000U*gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_MAJOR+
       100U*gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_MINOR+
       10U*gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_MICRO;
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_CREATE_DATE = "20-01-2020 10:54 +00200 (MON 20 JAN 2020 GMT+2)";
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_BUILD_DATE  = __DATE__ " " __TIME__;
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_SYNOPSIS    = "Scattering (AVX512 implementation) of Radiation by Moderately Non-Spherical Particles";

}

#include <cstdint>
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

		     const AVX512Vec16 ONE   = AVX512Vec16{1.0f};

		     const AVX512Vec16 TWO   = AVX512Vec16{2.0f};

		     const AVX512Vec16 PI    = AVX512Vec16{3.141592653589793f};

		     const AVX512Vec16 TWOPI = AVX512Vec16{6.283185307179586f};

	    } // anonymous

	    enum class EnsembleShapesAVX512 : uint32_t {

                       Cylindrical,
		       Spheroidal,
		       Chebyshev
	    }

	     // Low spatial and temporal frequency data type
	     struct HSColdAVX512_t {

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
		AVX512Vec16 * __restrict __ATTR_ALIGN__(8) m_pcs;
		//    ! Chebyshev particles radii in aggregate ensemble
		AVX512Vec16 * __restrict __ATTR_ALIGN__(8) m_radii;
		//       ! Chebyshev particles aggregate shape approximated by 3D parametric equations.
                //       ! Components form location of the particle in the ensemble.
                //       ! [3,np], where first dimension represents coordinate components
                //       ! second dimension represent number of particles.
	        float   * __restrict __ATTR_ALIGN__(8) m_pes;
		//        ! Chebyshev particles ensemble( per each particle) parametric equation in x - dimension (non-dimensional)
                //        ! [paramx,np]
		AVX512Vec16 * __restrict __ATTR_ALIGN__(8) m_ppx;
		//         ! Chebyshev particles ensemble (per each particle)  parametric equation in y  - dimension (non-dimensional)
		AVX512Vec16 * __restrict __ATTR_ALIGN__(8) m_ppy;
		//        ! Chebyshev particles ensemble (per each particle)  parametric equation in z  - dimension (non-dimensional)
		AVX512Vec16 * __restrict __ATTR_ALIGN__(8) m_ppz;
		//        ! Hydrometeor type
		char * __restrict __ATTR_ALIGN__(8) m_type;
		//         ! Ensemble shape  only 2 types supported for now: -- ( spheroid, chebyshev-particle)
		char * __restrict __ATTR_ALIGN__(8) m_shape;

	     } __ATTR_ALIGN__(64);

	     // ! High temporal and spatial frequency data type.
	     struct HSHotAVX512_t {
 
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

	     struct HMScatterersAVX512 {

                    HSColdAVX512_t   m_hsc;

		    HSHotAVX512_t    m_hsh;

		    HMScatterersAVX512() __ATTR_COLD__ __ATTR_ALIGN__(32);

		    HMScatterersAVX512(   const int32_t,
		                          const int32_t,
		                          const int32_t,
					  const int32_t,
					  const int32_t,
					  const int32_t,
					  const int32_t,
					  const char *,
					  const char *) __ATTR_COLD__ __ATTR_ALIGN__(32);

		     HMScatterersAVX512(const HMScatterersAVX512 &) = delete;

		     HMScatterersAVX512( HMScatterersAVX512 &&) = delete;

		     ~HMScatterersAVX512() noexcept(true);

		     HMScatterersAVX512 & operator=(const HMScatterersAVX512 &) = delete;

		     HMScatterersAVX512 & operator=( HMScatterersAVX512 &&) = delete;

		     bool ComputeShape_zmm16r4(
#if defined __ICC || defined __INTEL_COMPILER
                                            AVX512Vec16 * __restrict __ATTR_ALIGN__(64),
		                            AVX512Vec16 * __restrict __ATTR_ALIGN__(64)
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
					     const EnsembleShapesAVX512,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float,
					     const float) __ATTR_COLD__ __ATTR_ALIGN__(64) __ATTR_TCLONES_AVX_AVX512__;

	           bool ComputeXparam_zmm16r4(const AVX512Vec16 * __restrict __ATTR_ALIGN__(64),
		                              const AVX512Vec16 * __restrict __ATTR_ALIGN__(64)) __ATTR_COLD__ __ATTR_ALIGN__(32);
					     

		   bool ComputeYparam_zmm16r4(const AVX512Vec16 * __restrict __ATTR_ALIGN__(64),
		                              const AVX512Vec16 * __restrict __ATTR_ALIGN__(64)) __ATTR_COLD__ __ATTR_ALIGN__(32);
					    

		   bool ComputeZparam_zmm16r4(const AVX512Vec16 * __restrict __ATTR_ALIGN__(64),
		                              const AVX512Vec16 * __restrict __ATTR_ALIGN__(64)) __ATTR_COLD__ __ATTR_ALIGN__(32);
					     

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

     } // math


} // gms
















#endif /*__GMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AVX512_H__*/
