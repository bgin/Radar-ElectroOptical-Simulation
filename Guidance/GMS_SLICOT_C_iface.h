

#ifndef __GMS_SLICOT_C_IFACE_H__
#define __GMS_SLICOT_C_IFACE_H__
#include <complex>


#if defined(__cplusplus)

extern "C" {


/*
SUBROUTINE MB3OYZ( M, N, A, LDA, RCOND, SVLMAX, RANK, SVAL, JPVT, &
     TAU, DWORK, ZWORK, INFO )
    INTEGER            INFO, LDA, M, N, RANK
      DOUBLE PRECISION   RCOND, SVLMAX
!C     .. Array Arguments ..

#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            JPVT( * )
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: JPVT
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
#else
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: TAU
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: ZWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   DWORK( * ), SVAL( 3 )
#else
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(3) :: SVAL
#endif
*/
void MB30YZ(int *, int *, std::complex<double> __restrict *,
            int *, double *, double *, int *,
	    double __restrict *, int __restrict *,
	    std::complex<double> __restrict *, double __restrict *,
	    std::complex<double> __restrict *, int *);

/*
SUBROUTINE MB3PYZ( M, N, A, LDA, RCOND, SVLMAX, RANK, SVAL, JPVT,
  TAU, DWORK, ZWORK, INFO )
    INTEGER            INFO, LDA, M, N, RANK
      DOUBLE PRECISION   RCOND, SVLMAX
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            JPVT( * )
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: JPVT
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1      
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
#else
       COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
       COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: TAU
       COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: ZWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
       DOUBLE PRECISION   DWORK( * ), SVAL( 3 )
#else      
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
       DOUBLE PRECISION, DIMENSION(3) :: SVAL
#endif
*/
void MB30YZ(int *,
            int *,
	    std::complex<double> __restrict *,
	    int *,
	    double *,
	    double *,
	    int *,
	    double __restrict *,
	    int __restrict *,
	    std::coomplex<double> __restrict *,
	    double __restrict *,
	    std::complex<double> __restrict *,
	    int *);

/*
SUBROUTINE FB01QD( JOBK, MULTBQ, N, M, P, S, LDS, A, LDA, B,  &
                        LDB, Q, LDQ, C, LDC, R, LDR, K, LDK, TOL, &
                        IWORK, DWORK, LDWORK, INFO )
    CHARACTER         JOBK, MULTBQ
      INTEGER           INFO, LDA, LDB, LDC, LDK, LDQ, LDR, LDS, LDWORK, &
                        M, N, P
      DOUBLE PRECISION  TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           IWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           K(LDK,*), Q(LDQ,*), R(LDR,*), S(LDS,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: K
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S
*/
void FB01QD(const char *,
            const char *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
            int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE FB01SD( JOBX, MULTAB, MULTRC, N, M, P, SINV, LDSINV, &
          AINV, LDAINV, B, LDB, RINV, LDRINV, C, LDC,           &
          QINV, LDQINV, X, RINVY, Z, E, TOL, IWORK,             &
          DWORK, LDWORK, INFO )
   CHARACTER         JOBX, MULTAB, MULTRC
      INTEGER           INFO, LDAINV, LDB, LDC, LDQINV, LDRINV, LDSINV, &
                        LDWORK, M, N, P
      DOUBLE PRECISION  TOL
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           IWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  AINV(LDAINV,*), B(LDB,*), C(LDC,*), DWORK(*), &
                        E(*), QINV(LDQINV,*), RINV(LDRINV,*), RINVY(*), &
                        SINV(LDSINV,*), X(*), Z(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AINV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: E
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: QINV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RINV
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: RINVY
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: SINV
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: X
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: Z
#endif
*/
void FB01SD(const char *,
            const char *,
	    const char *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE FB01VD( N, M, L, P, LDP, A, LDA, B, LDB, C, LDC, Q, &
  LDQ, R, LDR, K, LDK, TOL, IWORK, DWORK, LDWORK,INFO)
    .. Scalar Arguments ..
      INTEGER           INFO, L, LDA, LDB, LDC, LDK, LDP, LDQ, LDR, &
                        LDWORK, M, N
      DOUBLE PRECISION  TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           IWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1      
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           K(LDK,*), P(LDP,*), Q(LDQ,*), R(LDR,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: K
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
#endif
*/
void FB01VD(int *,
            int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);






}











#endif










#endif /*__GMS_SLICOT_C_IFACE_H__*/
