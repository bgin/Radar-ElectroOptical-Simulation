

#ifndef __GMS_BLAS_IFACE_H__
#define __GMS_BLAS_IFACE_H__



namespace file_info {
#if defined _WIN64
  #include "../GMS_version.h"
#elif defined __linux
  #include "GMS_version.h"
#endif

   const unsigned int gGMS_BLAS_IFACE_MAJOR = gms::common::gVersionInfo.m_VersionMajor;
   const unsigned int gGMS_BLAS_IFACE_MINOR = gms::common::gVersionInfo.m_VersionMinor;
   const unsigned int gGMS_BLAS_IFACE_MICRO = gms::common::gVersionInfo.m_VersionMicro;
   const unsigned int gGMS_BLAS_IFACE_FULLVER =
     1000U*gGMS_BLAS_IFACE_MAJOR+100U*gGMS_BLAS_IFACE_MINOR+10U*gGMS_BLAS_IFACE_MICRO;
   const char * const pgGMS_BLAS_IFACE_CREATE_DATE = "05-12-2019 13:01 +00200 (THR 05 DEC 2019 GMT+2)";
   const char * const pgGMS_BLAS_IFACE_BUILD_DATE  = __DATE__ " " __TIME__;
   const char * const pgGMS_BLAS_IFACE_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
   const char * const pgGMS_BLAS_IFACE_SYNOPSIS    = "C wrappers for explicitly vectorized (type AVX512c8f64_t)  complex BLAS";
}

#if defined (GMS_CXX_98) || defined (GMS_CXX_11) || defined (GMS_CXX_14)

extern "C" {

  // This type is interoperable with Fortran AVX512c8f64_t
  // corresponding type.
  struct VC8F64_t {

    double re[8];
    double im[8];

  } __attribute__((aligned(64)));

  // This type is interoperable with Fortran ZMM8r8_t
  // corresponding type.
  struct V8F64_t {

       double v[8];
  } __attribute__((aligned(64)));

  

  /*  integer(kind=int4),                intent(in)    :: n
      type(AVX512c8f64_t),               intent(in)    :: za
      type(AVX512c8f64_t), dimension(*), intent(in)    :: zx
      !DIR$ ASSUME_ALIGNED zx:64
      integer(kind=int4),                intent(in)    :: incx
      type(AVX512c8f64_t), dimension(*), intent(inout) :: zy
      !DIR$ ASSUME_ALIGNED zy:64
      integer(kind=int4),                intent(in)    :: incy*/
#if defined __INTEL_COMPILER
     void  mod_blas_mp_gms_zaxpy_(const int,
		    const  VC8f64_t * __restrict,
		    const  VC8F64_t * __restrict,
		    const int,
		    VC8F64_t * __restrict,
		    const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zaxpy(const int,
		    const  VC8f64_t * __restrict,
		    const  VC8F64_t * __restrict,
		    const int,
		    VC8F64_t * __restrict,
		    const int);
#endif
//=====================================================================//
    /*  integer(kind=int4),                intent(in)  :: n
      type(AVX512c8f64_t), dimension(*), intent(in)  :: zx
      !DIR$ ASSUME_ALIGNED zx:64
      integer(kind=int4),                intent(in)  :: incx
      type(AVX512c8f64_t), dimension(*), intent(out) :: zy
      integer(kind=dint4),               intent(in)  :: incy*/
#if defined __INTEL_COMPILER     
     void mod_blas_mp_gms_zcopy_(const int,
		                 const VC8F64_t * __restrict,
		                 const int,
		                 VC8F64_t * __restrict,
		                 const int);
#elif defined __GFORTRAN__  || defined __GNUC__
     void __mod_blas_MOD_gms_zcopy(const int,
		                   const  VC8f64_t * __restrict,
		                   const  VC8F64_t * __restrict,
		                   const int,
		                   VC8F64_t * __restrict,
		                   const int);
#endif
//====================================================================//
     /*  integer(kind=int4),                intent(in) :: n
      type(AVX512c8f64_t), dimension(*), intent(in) :: zx
      !DIR$ ASSUME_ALIGNED zx:64
      integer(kind=int4),                intent(in) :: incx
      type(AVX512c8f64_t), dimension(*), intent(in) :: zy
      !DIR$ ASSUME_ALIGNED zy:64
      integer(kind=int4),                intent(in) :: incy*/
#if defined __INTEL_COMPILER
     VC8F64_t mod_blas_mp_gms_zdotc_(const int,
                                     const VC8F64_t * __restrict,
			             const int,
			             const VC8F64_t * __restrict,
			             const int);
#elif defined __GFORTRAN__  || defined __GNUC__
     VC8F64_t __mod_blas_MOD_gms_zdotc(const int,
                                       const VC8F64_t * __restrict,
			               const int,
			               const VC8F64_t * __restrict,
			               const int);
#endif
//===================================================================//
      /*   integer(kind=int4),                intent(in) :: n
      type(AVX512c8f64_t), dimension(*), intent(in) :: zx
      !DIR$ ASSUME_ALIGNED zx:64
      integer(kind=int4),                intent(in) :: incx
      type(AVX512c8f64_t), dimension(*), intent(in) :: zy
      !DIR$ ASSUME_ALIGNED zy:64
      integer(kind=int4),                intent(in) :: incy*/
#if defined __INTEL_COMPILER
     VC8F64_t mod_blas_mp_gms_zdotu_(const int,
                                     const VC8F64_t * __restrict,
			             const int,
			             const VC8F64_t * __restrict,
			             const int);
#elif defined __GFORTRAN__  || defined __GNUC__
     VC8F64_t __mod_blas_MOD_gms_zdot(const int,
                                     const VC8F64_t * __restrict,
			             const int,
			             const VC8F64_t * __restrict,
			             const int);
#endif
 //=================================================================//
     /*  integer(kind=int4),                intent(in)    :: n
      type(AVX512c8f64_t), dimension(*), intent(inout) :: cx
      !DIR$ ASSUME_ALIGNED cx:64
      integer(kind=int4),                intent(in)    :: incx
      type(AVX512c8f64_t), dimension(*), intent(inout) :: cy
      !DIR$ ASSUME_ALIGNED cy:64
      integer(kind=int4),                intent(in)    :: incy
      type(ZMM8r8_t),                    intent(in)    :: c ! scalar extended to vector
      !DIR$ ASSUME_ALIGNED c:64
      type(ZMM8r8_t),                    intent(in)    :: s ! scalar extended to vector	*/
#if defined __INTEL_COMPILER
     void mod_blas_mp_gms_zdrot_(const int,
                                 VC8F64_t * __restrict,
		                 const int,
		                 VC8F64_t * __restrict,
		                 const int,
		                 const V8F64_t * __restrict,
		                 const V8F64_t * __restrict);
#elif defined __GFORTRAN__  || defined __GNUC__
     void __mod_blas_MOD_gms_zdrot(const int,
                                 VC8F64_t * __restrict,
		                 const int,
		                 VC8F64_t * __restrict,
		                 const int,
		                 const V8F64_t * __restrict,
		                 const V8F64_t * __restrict);
#endif
//==================================================================//
     /*  integer(kind=int4),                intent(in)    :: n
       type(AVX512c8f64_t),               intent(in)    :: da
       type(AVX512c8f64_t), dimension(*), intent(inout) :: zx
       !DIR$ ASSUME_ALIGNED zx:64
       integer(kind=int4),                intent(in)    :: incx	*/
#if defined __INTEL_COMPILER
     void mod_blas_mp_gms_zdscal_(const int,
                                  const VC8F64_t * __restrict,
		                  VC8F64_t * __restrict,
		                  const int);
#elif defined __GFORTRAN__ || defined __GNUC__
     void __mod_blas_MOD_gms_zdscal(const int,
                                  const VC8F64_t * __restrict,
		                  VC8F64_t * __restrict,
		                  const int);
#endif
//=================================================================//

      /*  character(len=1),                      intent(in),value :: trans
        integer(kind=int4),                    intent(in),value :: m
        integer(kind=int4),                    intent(in),value :: n
        integer(kind=int4),                    intent(in),value :: kl
        integer(kind=int4),                    intent(in),value :: ku
        type(AVX512c8f64_t),                   intent(in)       :: alpha
        type(AVX512c8f64_t), dimension(lda,*), intent(in)       :: a
        !DIR$ ASSUME_ALIGNED a:64
        integer(kind=int4),                    intent(in),value :: lda
        type(AVX512c8f64_t), dimension(*),     intent(in)       :: x
        !DIR$ ASSUME_ALIGNED x:64
        integer(kind=int4),                    intent(in),value :: incx
        type(AVX512c8f64_t),                   intent(in)       :: beta
        type(AVX512c8f64_t), dimension(*),     intent(inout)    :: y
        integer(kind=int4),                    intent(in),value :: incy */
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zgbmv_(const char,
                                  const int,
		                  const int,
                                  const int,
				  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zgbmv(const char,
                                  const int,
		                  const int,
                                  const int,
				  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#endif
 //==========================================================================//

      
      /* character(len=1),                      intent(in),value    :: transa
       character(len=1),                      intent(in),value    :: transb
       integer(kind=int4),                    intent(in),value    :: m
       integer(kind=int4),                    intent(in),value    :: n
       integer(kind=int4),                    intent(in),value    :: k
       type(AVX512c8f64_t),                   intent(in)          :: alpha
       type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
       !DIR$ ASSUME_ALIGNED a:64
       integer(kind=int4),                    intent(in),value    :: lda
       type(AVX512c8f64_t), dimension(ldb,*), intent(in)          :: b
       !DIR$ ASSUME_ALIGNED b:64
       integer(kind=int4),                    intent(in),value    :: ldb
       type(AVX512c8f64_t),                   intent(in)          :: beta
       type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
       !DIR$ ASSUME_ALIGNED c:64
       integer(kind=int4),                    intent(in),value    :: ldc*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zgemm_(const char,
                                  const char,
				  const int,
				  const int,
				  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zgemm(const char,
                                  const char,
				  const int,
				  const int,
				  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#endif
//==================================================================================//

        /*  character(len=1),                      intent(in),value :: trans
          integer(kind=int4),                    intent(in),value :: m
          integer(kind=int4),                    intent(in),value :: n
          type(AVX512c8f64_t),                   intent(in)       :: alpha
          type(AVX512c8f64_t), dimension(lda,*), intent(in)       :: a
          !DIR$ ASSUME_ALIGNED a:64
          integer(kind=int4),                    intent(in),value :: lda
          type(AVX512c8f64_t), dimension(*),     intent(in)       :: x
          !DIR$ ASSUME_ALIGNED x:64
          integer(kind=int4),                    intent(in),value    :: incx
          type(AVX512c8f64_t),                   intent(in)          :: beta
          type(AVX512c8f64_t), dimension(*),     intent(inout)       :: y
          !DIR$ ASSUME_ALIGNED y:64
          integer(kind=int4),                    intent(in),value    :: incy*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zgemv_(const char,
                                  const int,
				  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
                                  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zgemv(const char,
                                  const int,
				  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
                                  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#endif
//===============================================================================//

       /* integer(kind=int4),                    intent(in),value    :: m
        integer(kind=int4),                    intent(in),value    :: n
        type(AVX512c8f64_t),                   intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
        !DIR$ ASSUME_ALIGNED x:64
        integer(kind=int4),                    intent(in),value    :: incx
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: y
        !DIR$ ASUME_ALIGNED y:64
        integer(kind=int4),                    intent(in),value    :: incy
        type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
        !DIR$ ASSUME_ALIGNED a:64
        integer(kind=int4),                    intent(in),value    :: lda*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zgerc_(const int,
                                  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
                                  const int,
                                  VC8F64_t * __restrict,
				  const int,
				  VC8F64_t * __restrict,
				  const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zgerc(const int,
                                  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
                                  const int,
                                  VC8F64_t * __restrict,
				  const int,
				  VC8F64_t * __restrict,
				  const int);
#endif
//===============================================================================//

       /*  integer(kind=int4),                    intent(in),value    :: m
        integer(kind=int4),                    intent(in),value    :: n
        type(AVX512c8f64_t),                   intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
        !DIR$ ASSUME_ALIGNED x:64
        integer(kind=int4),                    intent(in),value    :: incx
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: y
        !DIR$ ASSUME_ALIGNED y:64
        integer(kind=int4),                    intent(in),value    :: incy
        type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
        integer(kind=int4),                    intent(in),value    :: lda*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zgeru_(const int,
                                  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  VC8F64_t * __restrict,
				  const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zgeru(const int,
                                  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  VC8F64_t * __restrict,
				  const int);
#endif
//================================================================================//

       /*  character(len=1),                      intent(in),value    :: uplo
        integer(kind=int4),                    intent(in),value    :: n
        integer(kind=int4),                    intent(in),value    :: k
        type(AVX512c8f64_t),                   intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
        !DIR$ ASSUME_ALIGNED a:64
        integer(kind=int4),                    intent(in),value    :: lda
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
        !DIR$ ASSUME_ALIGNED x:64
        integer(kind=int4),                    intent(in),value    :: incx
        type(AVX512c8f64_t),                   intent(in)          :: beta
        type(AVX512c8f64_t), dimension(*),     intent(inout)       :: y
        !DIR$ ASSUME_ALIGNED y:64
        integer(kind=int4),                    intent(in),value    :: incy*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zhbmv_(const char,
                                  const int,
				  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zhbmv(const char,
                                  const int,
				  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#endif
//=============================================================================//				  

        /* character(len=1),                      intent(in),value    :: side
          character(len=1),                      intent(in),value    :: uplo
          integer(kind=int4),                    intent(in),value    :: m
          integer(kind=int4),                    intent(in),value    :: n
          type(AVX512c8f64_t),                   intent(in)          :: alpha
          type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
          !DIR$ ASSUME_ALIGNED a:64
          integer(kind=int4),                    intent(in),value    :: lda
          type(AVX512c8f64_t), dimension(ldb,*), intent(in)          :: b
          !DIR$ ASSUME_ALIGNED b:64
          integer(kind=int4),                    intent(in),value    :: ldb
          type(AVX512c8f64_t),                   intent(in)          :: beta
          type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
          !DIR$ ASSUME_ALIGNED c:64
          integer(kind=int4),                    intent(in),value    :: ldc*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zhemm_(const char,
                                  const char,
				  const int,
				  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zhemm(const char,
                                  const char,
				  const int,
				  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#endif
//=============================================================================//

       /*   character(len=1),                       intent(in),value    :: uplo
        integer(kind=int4),                     intent(in),value    :: n
        type(AVX512c8f64_t),                    intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(lda,*),  intent(in)          :: a
        !DIR$ ASSUME_ALIGNED a:64
        integer(kind=int4),                     intent(in),value    :: lda
        type(AVX512c8f64_t), dimension(*),      intent(in)          :: x
        !DIR$ ASSUME_ALIGNED x:64
        integer(kind=int4),                     intent(in),value    :: incx
        type(AVX512c8f64_t),                    intent(in)          :: beta
        type(AVX512c8f64_t), dimension(*),      intent(inout)       :: y
        integer(kind=int4),                     intent(in),value    :: incy*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zhemv_(const char,
                                  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
                                  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zhev(const char,
                                  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
                                  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#endif
//==============================================================================//

      /*   character(len=1),                      intent(in),value    :: uplo
        integer(kind=int4),                    intent(in),value    :: n
        type(ZMM8r4_t),                        intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
        !DIR$ ASSUME_ALIGNED x:64
        integer(kind=int4),                    intent(in),value    :: incx
        type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
        !DIR$ ASSUME_ALIGNED a:64
        integer(kind=int4),                    intent(in),value    :: lda*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zher_(const char,
                                 const int,
				 const V8F64_t * __restrict,
				 const VC8F64_t * __restrict,
				 const int,
				 VC8F64_t * __restrict,
				 const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zher(const char,
                                 const int,
				 const V8F64_t * __restrict,
				 const VC8F64_t * __restrict,
				 const int,
				 VC8F64_t * __restrict,
				 const int);
#endif
 //==============================================================================//

      /*   integer(kind=int4),                     intent(in),value    :: n
       type(AVX512c8f64_t),                    intent(in)          :: za
       type(AVX512c8f64_t), dimension(*),      intent(inout)       :: zx
       !DIR$ ASSUME_ALIGNED zx:64
       integer(kind=int4),                     intent(in),value    :: incx*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zscal_(const int,
                                  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zscal(const int,
                                  const VC8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#endif
//==============================================================================//

       /*  character(len=1),                      intent(in),value    :: uplo
         integer(kind=int4),                    intent(in),value    :: n
         type(AVX512c8f64_t),                   intent(in)          :: alpha
         type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
         !DIR$ ASSUME_ALIGNED x:64
         integer(kind=int4),                    intent(in),value    :: incx
         type(AVX512c8f64_t), dimension(*),     intent(in)          :: y
         !DIR$ ASSUME_ALIGNED y:64
         integer(kind=int4),                    intent(in),value    :: incy
         type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
         !DIR$ ASSUME_ALIGNED a:64
         integer(kind=int4),                    intent(in),value    :: lda*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zher2_(const char,
                                  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  VC8F64_t * __restrict,
				  const int);
#elif defined __GFORTRAN__ || defined __GNUC__
      void __mod_blas_MOD_gms_zher2(const char,
                                  const int,
				  const VC8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
				  const VC8F64_t * __restrict,
				  const int,
				  VC8F64_t * __restrict,
				  const int);
#endif
//===============================================================================//

       /* character(len=1),                      intent(in),value    :: uplo
          character(len=1),                      intent(in),value    :: trans
          integer(kind=int4),                    intent(in),value    :: n
          integer(kind=int4),                    intent(in),value    :: k
          type(AVX512c8f64_t),                   intent(in)          :: alpha
          type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
          !DIR$ ASSUME_ALIGNED a:64
          integer(kind=int4),                    intent(in),value    :: lda
          type(AVX512c8f64_t), dimension(ldb,*), intent(in)          :: b
          !DIR$ ASSUME_ALIGNED b:64
          integer(kind=int4),                    intent(in),value    :: ldb
          type(ZMM8r8_t),                        intent(in)          :: beta
          type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
          !DIR$ ASSUME_ALIGNED c:64
          integer(kind=int4),                    intent(in),value    :: ldc*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zher2k_(const char,
                                   const char,
				   const int,
				   const int,
				   const VC8F64_t * __restrict,
				   const VC8F64_t * __restrict,
				   const int,
				   const VC8F64_t * __restrict,
				   const int,
				   const V8F64_t * __restrict,
				   VC8F64_t * __restrict,
				   const int);
#elif defined __GFORTRAN__  || defined __GNUC__
      void __mod_blas_MOD_gms_zherk2(const char,
                                   const char,
				   const int,
				   const int,
				   const VC8F64_t * __restrict,
				   const VC8F64_t * __restrict,
				   const int,
				   const VC8F64_t * __restrict,
				   const int,
				   const V8F64_t * __restrict,
				   VC8F64_t * __restrict,
				   const int);
#endif
//================================================================================//

       /*  character(len=1),                      intent(in),value    :: uplo
        character(len=1),                      intent(in),value    :: trans
        integer(kind=int4),                    intent(in),value    :: n
        type(ZMM8r8_t),                        intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
        !DIR$ ASSUME_ALIGNED a:64
        integer(kind=int4),                    intent(in),value    :: lda
        type(ZMM8r8_t),                        intent(in)          :: beta
        type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
        !DIR$ ASUME_ALIGNED c:64
        integer(kind=int4),                    intent(in),value    :: ldc*/
#if defined __INTEL_COMPILER
      void mod_blas_mp_gms_zherk_(const char,
                                  const char,
				  const int,
				  const V8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
			          const V8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#elif defined __GFORTRAN__  || defined __GNUC__
      void __mod_blas_MOD_gms_zherk(const char,
                                  const char,
				  const int,
				  const V8F64_t * __restrict,
				  const VC8F64_t * __restrict,
				  const int,
			          const V8F64_t * __restrict,
				  VC8F64_t * __restrict,
				  const int);
#endif
//=================================================================================//			    
				  
} // extern "C"




#endif



















#endif /*__GMS_BLAS_IFACE_H__*/
