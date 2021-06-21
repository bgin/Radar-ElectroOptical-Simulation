

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




}

#endif










#endif /*__GMS_SLICOT_C_IFACE_H__*/
