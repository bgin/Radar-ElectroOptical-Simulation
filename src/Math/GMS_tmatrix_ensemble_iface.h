

#ifndef __GMS_TMATRIX_ENSEMBLE_IFACE_H__
#define __GMS_TMATRIX_ENSEMBLE_IFACE_H__


namespace file_info {

       const unsigned int gGMS_TMATRIX_ENSEMBLE_IFACE_MAJOR = 1;
       const unsigned int gGMS_TMATRIX_ENSEMBLE_IFACE_MINOR = 0;
       const unsigned int gGMS_TMATRIX_ENSEMBLE_IFACE_MICRO = 0;
       const unsigned int gGMS_TMATRIX_ENSEMBLE_IFACE_FULLVER =
             1000U*gGMS_TMATRIX_ENSEMBLE_IFACE_MAJOR+
	     100U*gGMS_TMATRIX_ENSEMBLE_IFACE_MINOR+
             10U*gGMS_TMATRIX_ENSEMBLE_IFACE_MICRO;
       const char * const pgGMS_TMATRIX_ENSEMBLE_IFACE_CREATE_DATE = "13-01-2020 17:21 +00200 (MON 13 JAN 2020 GMT+2)";
       const char * const pgGMS_TMATRIX_ENSEMBLE_IFACE_BUILD_DATE  = __DATE__ " " __TIME__;
       const char * const pgGMS_TMATRIX_ENSEMBLE_IFACE_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
       const char * const pgGMS_TMATRIX_ENSEMBLE_IFACE_SYNOPSYS    = "C-interface to Fortran T-Matrix particle ensemble model.";
}

#include "GMS_config.h"

#if (GMS_DEFAULT_CXX_VERSION) == 199711L || (GMS_DEFAULT_CXX_VERSION) == 201103L

extern "C" {


#if defined __ICC || defined __INTEL_COMPILER

     // Fortran version

    /* subroutine tmatrix_mps_driver(&
                                    idMie,small,MXINT,NADD,idscmt,sang,w,irat, &
                                    nL,idshp,shp,r0,cext,cabs,csca,assym,cextv,cabsv, &
                                    cscav,cbakv,cprv,cexts,cabss,cscas,cbaks,cprs, &
                                    dang,inat,pol,i11,i21,i12,i22,cexti,cabsi,cscai, &
                                    assymi,cpri,mue                                    )
    
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: tmatrix_mps_driver
          include 'tmatrix_mps_np.inc'
          
          integer(kind=int4)                   :: idMie
          real(kind=dp)                        :: small
          integer(kind=int4)                   :: MXINT,NADD,idscmt
          real(kind=dp)                        :: sang,w
          integer(kind=int4)                   :: irat,nL
          integer(kind=int4), dimension(nLp)   :: idshp
          real(kind=dp),    dimension(3,nLp)   :: shp
          real(kind=dp),    dimension(9,nLp)   :: r0
          real(kind=dp)                        :: cext,cabs,csca,assym,cextv,cabsv,  &
                                                  cscav,cbakv,cprv,cexts,cabss,cscas,&
                                                  cbaks,cprs
          real(kind=dp),    dimension(NANGMAX) :: dang,inat,pol,i11,i21,i12,i22
          real(kind=dp),    dimension(nLp)     :: cexti,cabsi,cscai,assymi,cpri
          real(kind=dp),    dimension(4,4,NANGMAX) :: mue */

     void  mod_tmatrix_mps_mp_tmatrix_mps_driver_(
						  int *,
						  double *,
						  int *,
						  int *,
						  int *,
						  double *,
						  double *,
						  int *,
						  int *,
						  int * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double *,
						  double *,
						  double *,
						  double *,
						  double *,
                                                  double *,
						  double *,
						  double *,
						  double *,
						  double *,
						  double *,
						  double *,
						  double *,
						  double *,
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64) );

#elif defined __GNUC__  || defined __GFORTRAN__ && !defined __INTEL_COMPILER

      __mod_tmatrix_mps_MOD_tmatrix_mps_driver(   
						  int *,
						  double *,
						  int *,
						  int *,
						  int *,
						  double *,
						  double *,
						  int *,
						  int *,
						  int * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double *,
						  double *,
						  double *,
						  double *,
						  double *,
                                                  double *,
						  double *,
						  double *,
						  double *,
						  double *,
						  double *,
						  double *,
						  double *,
						  double *,
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64),
						  double * __restrict __ATTR_ALIGN__(64) );

#endif








}


#endif












#endif /*__GMS_TMATRIX_ENSEMBLE_IFACE_H__*/
