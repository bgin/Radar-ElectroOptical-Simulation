
#ifndef __GMS_HYDROMETEOR_ENSEMBLE_SCATTERING_H__
#define __GMS_HYDROMETEOR_ENSEMBLE_SCATTERING_H__


namespace file_info {

   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_MAJOR = 1;
   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_MINOR = 0;
   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_MICRO = 0;
   const unsigned int gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_FULLVER =
       1000U*gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_MAJOR+
       100U*gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_MINOR+
       10U*gGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_MICRO;
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_CREATE_DATE = "08-01-2020 15:34 +00200 (WED 08 JAN 2020 GMT+2)";
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_BUILD_DATE  = __DATE__ " " __TIME__;
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
   const char * const pgGMS_HYDROMETEOR_ENSEMBLE_SCATTERING_SYNOPSIS    = "Scattering of Radiation by Moderately Non-Spherical Particles";
}

#include <cstdint>
#include "GMS_config.h"
#include "GMS_avxvecf32.h"


namespace gms {

        namespace  math {


	   // Low spatial and temporal frequency data type
	   struct HydrometeorScatterersCold_t {

                //     ! Number of particles  in aggregated emsemble
		int32_t m_np;
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
		AVXVec8 * __restrict __ATTR_ALIGN__(64) m_pcs;
		//    ! Chebyshev particles radii in aggregate ensemble
		AVXVec8 * __restrict __ATTR_ALIGN__(64) m_pradii;
		//       ! Chebyshev particles aggregate shape approximated by 3D parametric equations.
                //       ! Components form location of the particle in the ensemble.
                //       ! [3,np], where first dimension represents coordinate components
                //       ! second dimension represent number of particles.
		AVXVec8 * __restrict __ATTR_ALIGN__(64) m_pes;
		//        ! Chebyshev particles ensemble( per each particle) parametric equation in x - dimension (non-dimensional)
                //        ! [paramx,np]
		AVXVec8 * __restrict __ATTR_ALIGN__(64) m_ppx;
		//         ! Chebyshev particles ensemble (per each particle)  parametric equation in y  - dimension (non-dimensional)
		AVXVec8 * __restrict __ATTR_ALIGN__(64) m_ppy;
		//        ! Chebyshev particles ensemble (per each particle)  parametric equation in z  - dimension (non-dimensional)
		AVXVec8 * __restrict __ATTR_ALIGN__(64) m_ppz;
		//        ! Hydrometeor type
		char * m_type;
		//         ! Ensemble shape  only 2 types supported for now: -- ( spheroid, chebyshev-particle)
		char * m_shape;
	   } __ATTR_ALIGN__(64);


	   // ! High temporal and spatial frequency data type.

	   struct  HydrometeorScatterersHot_t {

               //  ! Yu-lin Xu part of the variables
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
	       double m_dang[1856]    __ATTR_ALIGN__(64);
	       double m_imat[1856]    __ATTR_ALIGN__(64);
	       double m_pol[1856]     __ATTR_ALIGN__(64);
	       double m_i11[1856]     __ATTR_ALIGN__(64);
	       double m_i21[1856]     __ATTR_ALIGN__(64);
	       double m_i12[1856]     __ATTR_ALIGN__(64);
	       double m_i22[1856]     __ATTR_ALIGN__(64);
	       double m_mue[4*4*1856] __ATTR_ALIGN__(64);
	       double * __restrict __ATTR_ALIGN__(64) m_cexti;
	       double * __restrict __ATTR_ALIGN__(64) m_cabsi;
	       double * __restrict __ATTR_ALIGN__(64) m_cscai;
	       double * __restrict __ATTR_ALIGN__(64) m_assymi;
	       double * __restrict __ATTR_ALIGN__(64) m_cpri;
	       //            ! Trajectory of Chebyshev particles ensemble, radial distance component (spherical coordinate system)
               //! [nt]
               float * __restrict __ATTR_ALIGN__(64) m_rdist;
	       //  ! Trajectory of Chebyshev particles ensemble, theta angle component (spherical coordinate system)
               //  ! [nt]
	       float * __restrict __ATTR_ALIGN__(64) m_theta;
	       //  ! Trajectory of Chebyshev particles ensemble, phi angle component (spherical coordinate system)
               //  ! [nt]
	       float * __restrict __ATTR_ALIGN__(64) m_phi;
	       //   ! Chebyshev particles ensemble fall speed 
               //   ![nt]
	       float * __restrict __ATTR_ALIGN__(64) m_vfall;
	  } __ATTR_ALIGN__(64);

    }  // math


} // gms




#endif /*__GMS_HYDROMETEOR_ENSEMBLE_SCATTERING_H__*/
