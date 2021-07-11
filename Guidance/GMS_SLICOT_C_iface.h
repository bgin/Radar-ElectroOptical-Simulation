

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

/*
 SUBROUTINE MB02OD( SIDE, UPLO, TRANS, DIAG, NORM, M, N, ALPHA, A,
   LDA, B, LDB, RCOND, TOL, IWORK, DWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER          DIAG, NORM, SIDE, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, M, N
      DOUBLE PRECISION   ALPHA, RCOND, TOL
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            IWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), DWORK(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
#endif
*/
void MB02OD(const char *,
            const char *,
	    const char *,
            const char *,
            const char *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *);

/*
 SUBROUTINE FB01TD( JOBX, MULTRC, N, M, P, SINV, LDSINV, AINV, &
        LDAINV, AINVB, LDAINB, RINV, LDRINV, C, LDC,          &
        QINV, LDQINV, X, RINVY, Z, E, TOL, IWORK,             &
        DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         JOBX, MULTRC
      INTEGER           INFO, LDAINB, LDAINV, LDC, LDQINV, LDRINV, &
                        LDSINV, LDWORK, M, N, P
      DOUBLE PRECISION  TOL
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           IWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  AINV(LDAINV,*), AINVB(LDAINB,*), C(LDC,*), &
                       DWORK(*), E(*), QINV(LDQINV,*), RINV(LDRINV,*), &
                       RINVY(*), SINV(LDSINV,*), X(*), Z(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AINV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AINVB
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: E
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: QINV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RINV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RINVY
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: SINV
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: X
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: Z
#endif
*/
void FB01TD(const char *,
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
SUBROUTINE MB04KD( UPLO, N, M, P, R, LDR, A, LDA, B, LDB, C, LDC, &
     TAU, DWORK )
!C     .. Scalar Arguments ..
      CHARACTER         UPLO
      INTEGER           LDA, LDB, LDC, LDR, M, N, P
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           R(LDR,*), TAU(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
#endif
*/

void MB04KD(const char *,
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
	    double __restrict *);

/*
 SUBROUTINE MB04ID(N, M, P, L, A, LDA, B, LDB, TAU, DWORK, LDWORK,INFO)
!C     .. Scalar Arguments ..
      INTEGER           INFO, L, LDA, LDB, LDWORK, M, N, P
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), TAU(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
#endif
*/
void MB04ID(int *,
            int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB04LD( UPLO, N, M, P, L, LDL, A, LDA, B, LDB, C, LDC &
                  TAU, DWORK )
!C     .. Scalar Arguments ..
      CHARACTER         UPLO
      INTEGER           LDA, LDB, LDC, LDL, M, N, P
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           L(LDL,*), TAU(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: L
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
#endif
*/
void MB04LD(const char *,
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
	    double __restrict *);

/*
SUBROUTINE SB01BD( DICO, N, M, NP, ALPHA, A, LDA, B, LDB, WR, WI, &
   NFP, NAP, NUP, F, LDF, Z, LDZ, TOL, DWORK, &
   LDWORK, IWARN, INFO )
!C     .. Scalar Arguments ..
      CHARACTER        DICO
      INTEGER          INFO, IWARN, LDA, LDB, LDF, LDWORK, LDZ, M, N, &
                       NAP, NFP, NP, NUP
      DOUBLE PRECISION ALPHA, TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
     DOUBLE PRECISION A(LDA,*), B(LDB,*), DWORK(*), F(LDF,*), &
          WI(*), WR(*), Z(LDZ,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WI
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Z
#endif
*/
void SB01BD(const char *,
            int *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    int *,
	    int *);

/*
 SUBROUTINE SB01BX(REIG,N,XR,XI,WR,WI,S,P)
 !C     .. Scalar Arguments ..
      LOGICAL          REIG
      INTEGER          N
      DOUBLE PRECISION P, S, XI ,XR
      !C     .. Array Arguments ..

      DOUBLE PRECISION WI(*), WR(*)
*/
void SB01BD(int *,
            int *,
	    double *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double *,
	    double *);

/*
 SUBROUTINE SB01BY(N,M,S,P,A B,F,TOL,DWORK,INFO)
!C     .. Scalar Arguments ..
      INTEGER           INFO, M, N
      DOUBLE PRECISION  P, S, TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      !DOUBLE PRECISION  A(N,*), B(N,*), DWORK(*), F(M,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F
#endif
*/
void SB01BY(int *,
            int *,
	    double *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double *,
	    double __restrict *,
	    int *);

/*
  SUBROUTINE MB03QD( DICO, STDOM, JOBU, N, NLOW, NSUP, ALPHA, &
       A, LDA, U, LDU, NDIM, DWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER        DICO, JOBU, STDOM
      INTEGER          INFO, LDA, LDU, N, NDIM, NLOW, NSUP
      DOUBLE PRECISION ALPHA
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION A(LDA,*), DWORK(*), U(LDU,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
#endif
*/
void MB03QD(const char *,
            const char *,
	    const char *,
	    int *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int *,
	    double __restrict *,
	    int *);

/*
 SUBROUTINE MB03QY(N,L,A,LDA,U,LDU,E1,E2,INFO)
!C     .. Scalar Arguments ..
      INTEGER          INFO, L, LDA, LDU, N
      DOUBLE PRECISION E1, E2
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION A(LDA,*), U(LDU,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
#endif      
*/
void MB03QY(int *,
            int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double *,
	    int *);

/*
SUBROUTINE TB04CD( JOBD, EQUIL, N, M, P, NPZ, A, LDA, B, LDB, C,  &
LDC, D, LDD, NZ, LDNZ, NP, LDNP, ZEROSR,                          &
ZEROSI, POLESR, POLESI, GAINS, LDGAIN, TOL,                       &
IWORK, DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER          EQUIL, JOBD
      DOUBLE PRECISION   TOL
      INTEGER            INFO, LDA, LDB, LDC, LDD, LDGAIN, LDNP, LDNZ, &
                         LDWORK, M, N, NPZ, P
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*), &
                         DWORK(*), GAINS(LDGAIN,*), POLESI(*),   &
                         POLESR(*), ZEROSI(*), ZEROSR(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GAINS
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: POLESI
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: POLESR
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZEROSI
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZEROSR
#endif    
*/
void TB04CD(const char *,
            const char *,
	    int *,
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
	    int __restrict *,
	    int *,
	    int __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
 SUBROUTINE TB01ID( JOB, N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,  &
      SCALE, INFO)
!C     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            INFO, LDA, LDB, LDC, M, N, P
      DOUBLE PRECISION   MAXRED
!C     ..
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ), &
           SCALE( * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: SCALE
#endif
*/
void TB01ID(const char *,
            int *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *);

/*
SUBROUTINE TB01ZD( JOBZ, N, P, A, LDA, B, C, LDC, NCONT, Z, LDZ,
    TAU, TOL, DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         JOBZ
      INTEGER           INFO, LDA, LDC, LDWORK, LDZ, N, NCONT, P
      DOUBLE PRECISION  TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(*), C(LDC,*), DWORK(*), TAU(*), &
           Z(LDZ,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TAU
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Z
#endif      
*/
void TB01ZD(const char *,
            int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double *,
	    double __restrict *,
	    int *,
	    int *);

/*
 SUBROUTINE MB01PD( SCUN, TYPE, M, N, KL, KU, ANRM, NBL, NROWS, A, &
   LDA, INFO )
!C     .. Scalar Arguments ..
      CHARACTER          SCUN, TYPE
      INTEGER            INFO, KL, KU, LDA, M, MN, N, NBL
      DOUBLE PRECISION   ANRM
      !C     .. Array Arguments ..
   
      INTEGER            NROWS ( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A( LDA, * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
#endif
*/
void MB01PD(const char *,
            const char *,
	    int *,
	    int *,
	    int *,
	    int *,
	    double *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    int * );

/*
 SUBROUTINE MB01QD( TYPE, M, N, KL, KU, CFROM, CTO, NBL, NROWS, A,
   LDA, INFO )
!C     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N, NBL
      DOUBLE PRECISION   CFROM, CTO
!C     ..
      !C     .. Array Arguments ..

      INTEGER            NROWS ( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1      
      DOUBLE PRECISION   A( LDA, * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
#endif
*/
void MB01QD(const char *,
            int *,
	    int *,
	    int *,
	    int *,
	    double *,
	    double *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE TB05AD( BALEIG, INITA, N, M, P, FREQ, A, LDA, B, LDB, &
  C, LDC, RCOND, G, LDG, EVRE, EVIM, HINVB,                      &
  LDHINV, IWORK, DWORK, LDWORK, ZWORK, LZWORK,INFO)
!C     .. Scalar Arguments ..
      CHARACTER         BALEIG, INITA
      INTEGER           INFO, LDA, LDB, LDC, LDG, LDHINV, LDWORK, &
                        LZWORK, M, N, P
      DOUBLE PRECISION  RCOND
      COMPLEX*16        FREQ
      !C     .. Array Arguments ..

      INTEGER           IWORK(*)
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), EVIM(*), &
           EVRE(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: EVIM
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: EVRE
#endif
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      COMPLEX*16        ZWORK(*), G(LDG,*), HINVB(LDHINV,*)
#else
      COMPLEX(16), DIMENSION(:), ALLOCATABLE :: ZWORK
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: G
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: HINVB
#endif
*/
void TB05AD(const char *,
            const char *,
	    int *,
	    int *,
	    int *,
	    std::complex<double> *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    std::complex<double> __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    std::complex<double> __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    std::complex<double> __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB02RZ(TRANS,N,NRHS,H,LDH,IPIV,B,LDB,INFO)
!C     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDB, LDH, N, NRHS
!C     ..
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            IPIV( * )
      COMPLEX*16         B( LDB, * ), H( LDH, * )
#else
      INTEGER,     DIMENSION(:), ALLOCATABLE :: IPIV
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: B,H
#endif
*/
void MB02RZ(const char *,
            int *,
	    int *,
	    std::complex<double> __restrict *,
	    int *,
	    int __restrict *,
	    std::complex<double> __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB02MD( DICO, HINV, UPLO, SCAL, SORT, N, A, LDA, G, &
    LDG, Q, LDQ, RCOND, WR, WI, S, LDS, U, LDU,                &
    IWORK, DWORK, LDWORK, BWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         DICO, HINV, SCAL, SORT, UPLO
      INTEGER           INFO, LDA, LDG, LDQ, LDS, LDU, LDWORK, N
      DOUBLE PRECISION  RCOND
      !C     .. Array Arguments ..

      LOGICAL           BWORK(*)
      INTEGER           IWORK(*)
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), DWORK(*), G(LDG,*), Q(LDQ,*), &
           S(LDS,*), U(LDU,*), WR(*), WI(*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: G
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: WI
#endif

*/
void SB02MD(const char *,
            const char *,
	    const char *,
	    const char *,
	    const char *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *);

/*
SUBROUTINE SB02MU(DICO, HINV, UPLO, N, A, LDA, G, LDG, Q, LDQ, S, &
      LDS, IWORK, DWORK, LDWORK, INFO)
!C     .. Scalar Arguments ..
      CHARACTER         DICO, HINV, UPLO
      INTEGER           INFO, LDA, LDG, LDQ, LDS, LDWORK, N
      !C     .. Array Arguments ..
 
      INTEGER           IWORK(*)
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1      
      DOUBLE PRECISION  A(LDA,*), DWORK(*), G(LDG,*), Q(LDQ,*), &
           S(LDS,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: G
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S
#endif
*/
void SB02MU(const char *,
            const char *,
	    const char *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB02ND( DICO, FACT, UPLO, JOBL, N, M, P, A, LDA, B, &
    LDB, R, LDR, IPIV, L, LDL, X, LDX, RNORM, F,               &
    LDF, OUFACT, IWORK, DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOBL, UPLO
      INTEGER           INFO, LDA, LDB, LDF, LDL, LDR, LDWORK, LDX, M, &
                        N, P
      DOUBLE PRECISION  RNORM
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1 
      INTEGER           IPIV(*), IWORK(*), OUFACT(2)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), F(LDF,*), &
           L(LDL,*), R(LDR,*), X(LDX,*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      INTEGER, DIMENSION(2) :: OUFACT
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: L
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif      
*/
void SB02ND(const char *,
            const char *,
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
	    int __restrict *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);


/*
SUBROUTINE SB02MT( JOBG, JOBL, FACT, UPLO, N, M, A, LDA, B, LDB, &
Q, LDQ, R, LDR, L, LDL, IPIV, OUFACT, G, LDG,    &
IWORK, DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         FACT, JOBG, JOBL, UPLO
      INTEGER           INFO, LDA, LDB, LDG, LDL, LDQ, LDR, LDWORK, M, &
                        N, OUFACT
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      !INTEGER           IPIV(*), IWORK(*)
      !DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), G(LDG,*), &
      !     L(LDL,*), Q(LDQ,*), R(LDR,*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: G
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: L
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
#endif
*/
void SB02MT(const char *,
            const char *,
	    const char *,
	    const char *,
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
	    int __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB01RB( SIDE, UPLO, TRANS, M, N, ALPHA, BETA, R, LDR, &
      A, LDA, B, LDB, INFO)
!C     .. Scalar Arguments ..
      CHARACTER         SIDE, TRANS, UPLO
      INTEGER           INFO, LDA, LDB, LDR, M, N
      DOUBLE PRECISION  ALPHA, BETA
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), R(LDR,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
#endif      
*/
void MB01RB(const char *,
            const char *,
	    const char *,
	    int *,
	    int *,
	    double *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB01RX( SIDE, UPLO, TRANS, M, N, ALPHA, BETA, R, LDR, &
      A, LDA, B, LDB, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         SIDE, TRANS, UPLO
      INTEGER           INFO, LDA, LDB, LDR, M, N
      DOUBLE PRECISION  ALPHA, BETA
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1  
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), R(LDR,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
#endif   
*/
void MB01RX(const char *,
            const char *,
	    const char *,
	    int *,
	    int *,
	    double *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB10ZD( N, M, NP, A, LDA, B, LDB, C, LDC, D, LDD, &
  FACTOR, AK, LDAK, BK, LDBK, CK, LDCK, DK, &
LDDK, RCOND, TOL, IWORK, DWORK, LDWORK, BWORK, &
INFO )
   !C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAK, LDB, LDBK, LDC, LDCK, LDD, &
                         LDDK, LDWORK, M, N, NP
      DOUBLE PRECISION   FACTOR, TOL
!C     ..
!C     .. Array Arguments ..
       
      INTEGER            IWORK( * )
      LOGICAL            BWORK( * )
      DOUBLE PRECISION, DIMENSION(6) :: RCOND
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1  
      DOUBLE PRECISION   A ( LDA,  * ), AK( LDAK, * ), B ( LDB,  * ), &
                         BK( LDBK, * ), C ( LDC,  * ), CK( LDCK, * ), &
                         D ( LDD,  * ), DK( LDDK, * ), DWORK( * )
                         
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
     
#endif
*/
void SB10ZD(int *,
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
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *);

/*
SUBROUTINE SB02OD( DICO, JOBB, FACT, UPLO, JOBL, SORT, N, M, P, A, &
LDA, B, LDB, Q, LDQ, R, LDR, L, LDL, RCOND, X,                     &
LDX, ALFAR, ALFAI, BETA, S, LDS, T, LDT, U,                        &
LDU, TOL, IWORK, DWORK, LDWORK, BWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOBB, JOBL, SORT, UPLO
      INTEGER           INFO, LDA, LDB, LDL, LDQ, LDR, LDS, LDT, LDU, &
                        LDWORK, LDX, M, N, P
      DOUBLE PRECISION  RCOND, TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1 
      LOGICAL           BWORK(*)
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), ALFAI(*), ALFAR(*), B(LDB,*), BETA(*), &
                        DWORK(*), L(LDL,*), Q(LDQ,*), R(LDR,*),          &
                        S(LDS,*), T(LDT,*), U(LDU,*), X(LDX,*)
#else
      LOGICAL, DIMENSION(:), ALLOCATABLE :: BWORK
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ALFAI
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ALFAR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: BETA
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: L
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif     
*/
void SB02OD(const char *,
            const char *,
	    const char *,
	    const char *,
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
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
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
	    int __restrict *,
	    int *);

/*
SUBROUTINE SB02OY( TYPE, DICO, JOBB, FACT, UPLO, JOBL, JOBE, N, M, &
  P, A, LDA, B, LDB, Q, LDQ, R, LDR, L, LDL, E,                    &
  LDE, AF, LDAF, BF, LDBF, TOL, IWORK, DWORK,                      &
  LDWORK, INFO)
!C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOBB, JOBE, JOBL, TYPE, UPLO
      INTEGER           INFO, LDA, LDAF, LDB, LDBF, LDE, LDL, LDQ, LDR, &
                        LDWORK, M, N, P
      DOUBLE PRECISION  TOL
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), AF(LDAF,*), B(LDB,*), BF(LDBF,*),  &
           DWORK(*), E(LDE,*), L(LDL,*), Q(LDQ,*), R(LDR,*)
#else
      INTEGER, DIMENSION(:),ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AF
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BF
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: E
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: L
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
#endif
*/
void SB02OY(const char *,
            const char *,
	    const char *,
	    const char *,
	    const char *,
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
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB03TD( JOB, FACT, TRANA, UPLO, LYAPUN, N, SCALE, A,
LDA, T, LDT, U, LDU, C, LDC, X, LDX, SEP,
RCOND, FERR, WR, WI, IWORK, DWORK, LDWORK,INFO)
!C     .. Scalar Arguments ..
      CHARACTER          FACT, JOB, LYAPUN, TRANA, UPLO
      INTEGER            INFO, LDA, LDC, LDT, LDU, LDWORK, LDX, N
      DOUBLE PRECISION   FERR, RCOND, SCALE, SEP
!C     ..
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), DWORK( * ), &
                         T( LDT, * ), U( LDU, * ), WI( * ), WR( * ), &
                         X( LDX, * )
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WI
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif
*/
void SB03TD(const char *,
            const char *,
	    const char *,
	    const char *,
	    const char *,
	    int *,
	    double *,
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
	    double *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB01RU( UPLO, TRANS, M, N, ALPHA, BETA, R, LDR, A, LDA,
  X, LDX, DWORK, LDWORK, INFO)
!C     .. Scalar Arguments ..
      CHARACTER         TRANS, UPLO
      INTEGER           INFO, LDA, LDR, LDWORK, LDX, M, N
      DOUBLE PRECISION  ALPHA, BETA
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      !DOUBLE PRECISION  A(LDA,*), DWORK(*), R(LDR,*), X(LDX,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif
*/
void MB01RU(const char *,
            const char *,
	    int *,
	    int *,
	    double *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB03MY( TRANA, N, A, LDA, C, LDC, SCALE, INFO)
!C     .. Scalar Arguments ..
      CHARACTER          TRANA
      INTEGER            INFO, LDA, LDC, N
      DOUBLE PRECISION   SCALE
!C     ..
!C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * )
!C     ..
*/
void SB03MY(const char *,
            int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    int *);

/*
SUBROUTINE SB03MW( LTRAN, LUPPER, T, LDT, B, LDB, SCALE, X, LDX,
  XNORM, INFO )
!C     .. Scalar Arguments ..
      LOGICAL            LTRAN, LUPPER
      INTEGER            INFO, LDB, LDT, LDX
      DOUBLE PRECISION   SCALE, XNORM
!C     ..
      !C     .. Array Arguments ..
  
      DOUBLE PRECISION   B( LDB, * ), T( LDT, * ), X( LDX, * )
*/
void SB03MW(int *,
            int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double *,
	    int *);

/*
SUBROUTINE SB03QD( JOB, FACT, TRANA, UPLO, LYAPUN, N, SCALE, A, &
LDA, T, LDT, U, LDU, C, LDC, X, LDX, SEP,  &
RCOND, FERR, IWORK, DWORK, LDWORK, INFO)
!C     .. Scalar Arguments ..
      CHARACTER          FACT, JOB, LYAPUN, TRANA, UPLO
      INTEGER            INFO, LDA, LDC, LDT, LDU, LDWORK, LDX, N
      DOUBLE PRECISION   FERR, RCOND, SCALE, SEP
!C     ..
      !C!     .. Array Arguments ..
  
      INTEGER            IWORK( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), DWORK( * ), &
           T( LDT, * ), U( LDU, * ), X( LDX, * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif
*/
void SB03QD(const char *,
            const char *,
	    const char *,
	    const char *,
	    const char *,
	    int *,
	    double *,
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
	    double *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE MB01UW( SIDE, TRANS, M, N, ALPHA, H, LDH, A, LDA, &
     DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         SIDE, TRANS
      INTEGER           INFO, LDA, LDH, LDWORK, M, N
      DOUBLE PRECISION  ALPHA
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), DWORK(*), H(LDH,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: H
#endif
*/
void MB01UW(const char *,
            const char *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int *);

/*
 SUBROUTINE SB03QX( TRANA, UPLO, LYAPUN, N, XANORM, T, LDT, U, LDU, &
       R, LDR, FERR, IWORK, DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER          LYAPUN, TRANA, UPLO
      INTEGER            INFO, LDR, LDT, LDU, LDWORK, N
      DOUBLE PRECISION   FERR, XANORM
!C     ..
      !C     .. Array Arguments ..

      INTEGER            IWORK( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1      
     DOUBLE PRECISION   DWORK( * ), R( LDR, * ), T( LDT, * ), &
     U( LDU, * )
#else 
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
#endif
*/
void SB03QX(const char *,
            const char *,
	    const char *,
	    int *,
	    double *,
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
SUBROUTINE SB03QY( JOB, TRANA, LYAPUN, N, T, LDT, U, LDU, X, LDX, &
     SEP, THNORM, IWORK, DWORK, LDWORK, INFO)
!C     .. Scalar Arguments ..
      CHARACTER          JOB, LYAPUN, TRANA
      INTEGER            INFO, LDT, LDU, LDWORK, LDX, N
      DOUBLE PRECISION   SEP, THNORM
!C     ..
      !C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER            IWORK( * )
      DOUBLE PRECISION   DWORK( * ), T( LDT, * ), U( LDU, * ), &
           X( LDX, * )
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif
*/
void SB03QY(const char *,
            const char *,
	    const char *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE FD01AD( JP, L, LAMBDA, XIN, YIN, EFOR, XF, EPSBCK, &
CTETA, STETA, YQ, EPOS, EOUT, SALPH, IWARN, &
INFO )
!C     .. Scalar Arguments ..
      CHARACTER         JP
      INTEGER           INFO, IWARN, L
      DOUBLE PRECISION  EFOR, EOUT, EPOS, LAMBDA, XIN, YIN
!C     .. Array Arguments ..
      !DOUBLE PRECISION  CTETA(*), EPSBCK(*), SALPH(*), STETA(*), XF(*), &
      !                  YQ(*)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CTETA
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EPSBCK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SALPH
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: STETA
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XF
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: YQ

*/
void FD01AD(const char *,
            int *,
	    double *,
	    double *,
	    double *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double *,
	    double *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE TD05AD( UNITF, OUTPUT, NP1, MP1, W, A, B, VALR, VALI, &
     INFO )
!C     .. Scalar Arguments ..
      CHARACTER         OUTPUT, UNITF
      INTEGER           INFO, MP1, NP1
      DOUBLE PRECISION  VALI, VALR, W
!C     .. Array Arguments ..
      !DOUBLE PRECISION  A(*), B(*)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: B
*/
void TD05AD(const char *,
            const char *,
	    int *,
	    int *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double *,
	    double *,
	    int *);

/*
SUBROUTINE TF01QD( NC, NB, N, IORD, AR, MA, H, LDH, INFO )
!C     .. Scalar Arguments ..
      INTEGER           INFO, LDH, N, NB, NC
!C     .. Array Arguments ..
      INTEGER           IORD(*)
      DOUBLE PRECISION  AR(*), H(LDH,*), MA(*)
*/
void TF01QD(int *,
            int *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE TF01RD( NA, NB, NC, N, A, LDA, B, LDB, C, LDC, H, LDH, &
     DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LDC, LDH, LDWORK, N, NA, NB, NC
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), H(LDH,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: H
#endif  
*/
void TF01RD(int *,
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
	    int *);

/*
SUBROUTINE TC05AD( LERI, M, P, SVAL, INDEX, PCOEFF, LDPCO1, &
LDPCO2, QCOEFF, LDQCO1, LDQCO2, RCOND, CFREQR,&
LDCFRE, IWORK, DWORK, ZWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         LERI
      INTEGER           INFO, LDCFRE, LDPCO1, LDPCO2, LDQCO1, LDQCO2, M, P
    
      DOUBLE PRECISION  RCOND
      COMPLEX*16        SVAL
!C     .. Array Arguments ..
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      INTEGER           INDEX(*), IWORK(*)
      DOUBLE PRECISION  DWORK(*), PCOEFF(LDPCO1,LDPCO2,*), &
                        QCOEFF(LDQCO1,LDQCO2,*)
      COMPLEX*16        CFREQR(LDCFRE,*), ZWORK(*)
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: INDEX, IWORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: PCOEFF
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QCOEFF
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: CFREQR
      COMPLEX(16), DIMENSION(:), ALLOCATABLE :: ZWORK
#endif
*/
void TC05AD(const char *,
            int *,
	    int *,
	    std::complex<double> *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    int *,
	    double *,
	    std::complex<double> __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    std::complex<double> __restrict *,
	    int *);

/*
SUBROUTINE SB10KD( N, M, NP, A, LDA, B, LDB, C, LDC, FACTOR, &
AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK, RCOND, &
IWORK, DWORK, LDWORK, BWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAK, LDB, LDBK, LDC, LDCK, LDDK, &
                         LDWORK, M, N, NP
      DOUBLE PRECISION   FACTOR
!C     ..
!C     .. Array Arguments ..
       DOUBLE PRECISION, DIMENSION(4) :: RCOND
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
     INTEGER            IWORK( * )
     LOGICAL            BWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AK( LDAK, * ), B( LDB, * ), &
                       BK( LDBK, * ), C( LDC, * ), CK( LDCK, * ), &
                         DK( LDDK, * ), DWORK( * )
#else
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      LOGICAL, DIMENSION(:), ALLOCATABLE :: BWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
#endif     
*/
void SB10KD(int *,
            int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *);

/*
SUBROUTINE SB10AD( JOB, N, M, NP, NCON, NMEAS, GAMMA, A, LDA, &
B, LDB, C, LDC, D, LDD, AK, LDAK, BK, LDBK, CK, &
LDCK, DK, LDDK, AC, LDAC, BC, LDBC, CC, LDCC, &
DC, LDDC, RCOND, GTOL, ACTOL, IWORK, LIWORK, &
DWORK, LDWORK, BWORK, LBWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER            INFO, JOB, LBWORK, LDA, LDAC, LDAK, LDB, LDBC, &
                         LDBK, LDC, LDCC, LDCK, LDD, LDDC, LDDK, LDWORK, &
                         LIWORK, M, N, NCON, NMEAS, NP
      DOUBLE PRECISION   ACTOL, GAMMA, GTOL
!C     ..
!C     .. Array Arguments ..
      LOGICAL            BWORK( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION, DIMENSION(4) :: RCOND
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A( LDA, * ), AC( LDAC, * ), AK( LDAK, * ), &
                        B( LDB, * ), BC( LDBC, * ), BK( LDBK, * ), &
                        C( LDC, * ), CC( LDCC, * ), CK( LDCK, * ), &
                        D( LDD, * ), DC( LDDC, * ), DK( LDDK, * ), &
                        DWORK( * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
#endif
*/
void SB10AD(int *,
            int *,
	    int *,
	    int *,
	    int *,
	    int *,
	    double *,
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
	    double *,
	    double *,
	    int __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
            int *,
	    int *);
	    
/*
SUBROUTINE SB10LD( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC, &
D, LDD, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK, &
AC, LDAC, BC, LDBC, CC, LDCC, DC, LDDC, IWORK, &
DWORK, LDWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAC, LDAK, LDB, LDBC, LDBK, LDC, &
                        LDCC, LDCK, LDD, LDDC, LDDK, LDWORK, M, N, &
                        NCON, NMEAS, NP
!C     ..
!C     .. Array Arguments ..
      INTEGER            IWORK( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION   A( LDA, * ), AC( LDAC, * ), AK( LDAK, * ),
                      B( LDB, * ), BC( LDBC, * ), BK( LDBK, * ),
                      C( LDC, * ), CC( LDCC, * ), CK( LDCK, * ),
                      D( LDD, * ), DC( LDDC, * ), DK( LDDK, * ),
                      DWORK( * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DK
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
#endif
*/
void SB10LD(int *,
            int *,
	    int *,
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
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
            int __restrict *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB10PD( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC, &
D, LDD, TU, LDTU, TY, LDTY, RCOND, TOL, DWORK, &
LDWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LDC, LDD, LDTU, LDTY, LDWORK, &
                         M, N, NCON, NMEAS, NP
      DOUBLE PRECISION   TOL
!C     ..
!C     .. Array Arguments ..
       DOUBLE PRECISION, DIMENSION(2)  :: RCOND
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ), &
                       D( LDD, * ), DWORK( * ), &
                        TU( LDTU, * ), TY( LDTY, * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TU
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TY
#endif      
*/
void SB10PD(int *,
            int *,
	    int *,
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
	    double *,
	    double __restrict *,
	    int *,
	    int *);

/*
SUBROUTINE SB10QD( N, M, NP, NCON, NMEAS, GAMMA, A, LDA, B, LDB,
C, LDC, D, LDD, F, LDF, H, LDH, X, LDX, Y, LDY,
XYCOND, IWORK, DWORK, LDWORK, BWORK, INFO )
!C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LDC, LDD, LDF, LDH, LDWORK, &
                         LDX, LDY, M, N, NCON, NMEAS, NP
      DOUBLE PRECISION   GAMMA
!C     ..
!C     .. Array Arguments ..
      INTEGER            IWORK( * )
      LOGICAL            BWORK( * )
      DOUBLE PRECISION,DIMENSION(2) :: XYCOND
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
     DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ), &
                        D( LDD, * ), DWORK( * ),  F( LDF, * ), &
                        H( LDH, * ), X( LDX, * ), &
                        Y( LDY, * )
#else     
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION,DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: F
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: H
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: X
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: Y
#endif
*/
void SB10QD((int *,
            int *,
	    int *,
	    int *,
	    int *,
	    int *,
	    double *,
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
	    double __restrict *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *);

/*

SUBROUTINE SB02RD( JOB, DICO, HINV, TRANA, UPLO, SCAL, SORT, FACT, &
LYAPUN, N, A, LDA, T, LDT, V, LDV, G, LDG, Q, &
LDQ, X, LDX, SEP, RCOND, FERR, WR, WI, S, LDS, &
IWORK, DWORK, LDWORK, BWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, HINV, JOB, LYAPUN, SCAL, SORT, &
                        TRANA, UPLO
      INTEGER           INFO, LDA, LDG, LDQ, LDS, LDT, LDV, LDWORK, LDX, &
                       N
      DOUBLE PRECISION  FERR, RCOND, SEP
!C     .. Array Arguments ..
      LOGICAL           BWORK(*)
      INTEGER           IWORK(*)
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A(LDA,*), DWORK(*), G(LDG,*), Q(LDQ,*),
                      S(LDS,*), T(LDT,*), V(LDV,*), WI(*), WR(*),
                       X(LDX,*)
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: G
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Q
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: S
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: V
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WI
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif   
*/
void SB02RD(const char *,
            const char *,
	    const char *,
	    const char *,
	    const char *,
	    const char *,
	    const char *,
	    const char *,
	    const char *,
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
	    double *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    int *);

/*
SUBROUTINE MB02PD( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, &
IWORK, DWORK, INFO )
!C     .. Scalar Arguments ..
      CHARACTER         EQUED, FACT, TRANS
      INTEGER           INFO, LDA, LDAF, LDB, LDX, N, NRHS
      DOUBLE PRECISION  RCOND
!C     ..
!C     .. Array Arguments ..
      INTEGER           IPIV( * ), IWORK( * )
#if (GMS_SLICOT_USE_MKL_LAPACK) == 1
      DOUBLE PRECISION  A( LDA, * ), AF( LDAF, * ), B( LDB, * ), &
                      BERR( * ), C( * ), DWORK( * ), FERR( * ), &
                        R( * ), X( LDX, * )
#else
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AF
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BERR
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FERR
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: R
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
#endif
*/
void MB02PD(const char *,
            const char *,
	    int *,
	    int *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    int __restrict *,
	    const char *,
	    double __restrict *,
	    double __restrict *,
	    double __restrict *,
	    int *,
	    double __restrict *,
	    int *,
	    double *,
	    double __restrict *,
	    double __restrict *,
	    int __restrict *,
	    double __restrict *,
	    int *);

/*
*/



}

#endif










#endif /*__GMS_SLICOT_C_IFACE_H__*/
