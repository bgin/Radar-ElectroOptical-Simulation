
#ifndef __GMS_LAPACK_C_IFACE_H__
#define __GMS_LAPACK_C_IFACE_H__


#if defined(__cplusplus)


extern "C" {


/*
       SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLARFG
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLARFG
#endif
       implicit none
  
!*
!*  -- LAPACK auxiliary routine (version 3.8.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
!*
!*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*/

void DLARFG(int *, double *, double __restrict *, int *, double *);

/*
    SUBROUTINE DLATZM(SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK)
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLATZM
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLATZM
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION   C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C1
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C2
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: V
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DLATZM(const char *, int *, int *, double __restrict *,
            int *, double *, double __restrict *, double __restrict *,
	    int *, double __restrict *);

/*
     SUBROUTINE DLARFB(SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
   T, LDT, C, LDC, WORK, LDWORK)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLARFB
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLARFB
#endif
       use omp_lib
       implicit none
       
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2013
!*
!*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!*     ..
!*     .. Array Arguments ..
     ! DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
      ! $                   WORK( LDWORK, * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: V
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: WORK
*/
void DLARFB(const char *, const char *, const char *, const char *,
            int *, int *, int *, double __restrict *, int *,
	    double __restrict *, int *, double __restrict *,
	    int *, double __restrict *, int *);

/*
     SUBROUTINE DLARF(SIDE, M, N, V, INCV, TAU, C, LDC, WORK)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLARF
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLARF
#endif
  
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: V
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DLARF(const char *, int *, int *, double __restrict *,
           int *, double *, double __restrict *, int *,
	   double __restrict *);

/*
    SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLARFT
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLARFT
#endif
      use omp_lib
      implicit none
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: V
*/
void DLARFT(const char *, const char *, int *, int *,
            double __restrict *, int *, double *,
	    double __restrict *, int *);

/*
    SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGETRF
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGETRF
#endif
       implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!*     ..
!*     .. Array Arguments ..
      INTEGER            IPIV( * )
      !DOUBLE PRECISION   A( LDA, * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
*/
void DGETRF(int *, int *, double __restrict *, int *,
            int __restrict *, int *);

/*
   SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, INFO)
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGEBRD
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGEBRD
#endif
       implicit none
!*
!*  -- LAPACK computational routine (version 3.8.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ), &
      !                   TAUQ( * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: E
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAUP
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAUQ
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DGEBRD(int *, int *, double __restrict *, int *, double __restrict *,
            double __restrict *, double __restrict *, double __restrict *,
	    double __restrict *, int *, int *);

/*
   SUBROUTINE DGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGEBD2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGEBD2
#endif
       implicit none
!*
!*  -- LAPACK computational routine (version 3.7.1) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2017
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!*     ..
!*     .. Array Arguments ..
    !  DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ),
      !$                   TAUQ( * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: E
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAUP
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAUQ
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DBEBD2(int *, int *, double __restrict *, int *, double __restrict *,
            double __restrict *, double __restrict *, double __restrict *,
	    double __restrict *, int *);

/*
SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, &
     VT, LDVT, WORK, LWORK, INFO )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGESVD
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGESVD
#endif
      implicit none
!*
!*  -- LAPACK driver routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     April 2012
!*
!*     .. Scalar Arguments ..
      CHARACTER          JOBU, JOBVT
      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
!*     ..
!*     .. Array Arguments ..
     ! DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
      !$                   VT( LDVT, * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: S
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: VT
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DGESVD(const char *, const char *, int *, int *, double __restrict *,
            int *, double __restrict *, double __restrict *, int *,
	    double __restrict *, int *, double __restrict *, int *, int *);

/*
SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, &
     LDU, C, LDC, WORK, INFO )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DBDSQR
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DBDQSR
#endif
       implicit none
!*
!*  -- LAPACK computational routine (version 3.7.1) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2017
!*
!*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
!*     ..
!*  !   .. Array Arguments ..
    !  DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ),
      ! $                   VT( LDVT, * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: E
      DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: VT
      DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: WORK
*/
void DBDSQR(const char *, int *, int *, int *, int *, double __restrict *,
            double __restrict *, double __restrict *, int *,
	    double __restrict *, int *, double __restrict *,
	    int *, double __restrict *, int *);

/*
  SUBROUTINE DLASQ1( N, D, E, WORK, INFO )
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASQ1
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASQ1
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, N
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION   D( * ), E( * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
      
*/
void DLASQ1(int *, double __restrict *, double __restrict *,
            double __restrict *, int *);

/*
    SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASCL
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLASCL
#endif
       use omp_lib
       implicit none
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION   A( LDA, * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
*/
void DLASCL(const char *, int *, int , double *, double *,
            int *, int *, double __restrict *, int *, int *);

/*
 SUBROUTINE DLASQ2( N, Z, INFO )
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASQ2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASQ2
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   Z( * )
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z
*/
void DLASQ2(int *, double __restrict *, int *);

/*
SUBROUTINE DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, &
ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, &
DN2, G, TAU )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASQ3
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASQ3
#endif
      implicit none

!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2016
!*
!*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            I0, ITER, N0, NDIV, NFAIL, PP
      DOUBLE PRECISION   DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, &
                         QMAX, SIGMA, TAU
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   Z( * )
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z
*/
void DLASQ3(int *, int *, double __restrict *, int *, double *, double *,
            double *, double *, int *, int *, int *, int *, int *,
	    double *, double *, double *, double *, double *, double *,
	    double *);

/*
SUBROUTINE DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, &
     DN1, DN2, TAU, TTYPE, G )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASQ4
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASQ4
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.1) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            I0, N0, N0IN, PP, TTYPE
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   Z( * )
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z
*/
void DLASQ4(int *, int *, double __restrict *, int *, int *, double *,
            double *, double *, double *, double *, double *, double *,
	    int *, double *);

/*
  SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, &
     DN, DNM1, DNM2, IEEE, EPS )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASQ5
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASQ5
#endif
      implicit none

!*
!*  -- LAPACK computational routine (version 3.7.1) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2017
!*
!*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            I0, N0, PP
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU, &
                        SIGMA, EPS
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION   Z( * )
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z
*/
void DLASQ5(int *, int *, double __restrict *, int *, double *, double *,
            double *, double *, double *, double *, double *, int *, double *);

/*
 SUBROUTINE DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, &
      DNM1, DNM2 )
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASQ6
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASQ6
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            I0, N0, PP
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION   Z( * )
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z
*/
void DLASQ6(int *, int *, double __restrict *, int *, double *, double *,
            double *, double *, double *, double *);

/*
 SUBROUTINE DLASRT( ID, N, D, INFO )
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASRT
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASRT
#endif
       implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   D( * )
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: D
*/
void DLASRT(int *, int *, double __restrict *, int *);

/*
  SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASR
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASR
#endif
    implicit none
    
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: S
*/

void DLASR(const char *, const char *, const char *, int *,
           int *, double __restrict *, double __restrict *,
	   double __restrict *, int *);

/*
 SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGELQF
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGELQF
#endif
    implicit none
    
!*
!*  -- LAPACK computational routine (version 3.9.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2019
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DGELQF(int *, int *, double __restrict *, int *,
            double __restrict *, int *, double __restrict *,
	    int *, int *);

/*
 SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGEQRF
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGEQRF
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.9.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2019
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DGEQRF(int *, int *, double __restrict *, int *,
            double __restrict *, double __restrict *,
	    int *, int *);

/*
 SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGEQR2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGEQR2
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.9.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2019
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DGEQR2(int *, int *, double __restrict *, int *,
            double __restrict *, double __restrict *,
	    int *);

/*
 SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DORGBR
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_av512 :: DORGBR
#endif
      use omp_lib
      implicit none

!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     April 2012
!*
!*     .. Scalar Arguments ..
      CHARACTER          VECT
      INTEGER            INFO, K, LDA, LWORK, M, N
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
         DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DORGBR(const char *, int *, int *, int *, double __restrict *,
            int *, double __restrict *, double __restrict *,
	    int *, int *);

/*
  SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DORGLQ
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_av512 :: DORGLQ
#endif
      use omp_lib
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DORGLQ(int *, int *, int *, double __restrict *, int *,
            double __restrict *, double __restrict *, int *,
	    int *);

/*
  SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )
    !DIR$ ATTRIBUTES FORCEINLINE :: DORGL2
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DORGL2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_av512 :: DORGL2
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DORGL2(int *, int *, int *, double __restrict *, int *,
            double __restrict *, double __restrict *, int *,
	    int *);

/*
  SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DORGQR
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_av512 :: DORGQR
#endif
    implicit none
    
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DORGQR(int *, int *, int *, double __restrict *, int *,
            double __restrict *, double __restrict *, int *,
	    int *);

/*
 SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
    !DIR$ ATTRIBUTES FORCEINLINE :: DORG2R
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DORG2R
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_av512 :: DORG2R
#endif
      use omp_lib
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DORG2R(int *, int *, int *, double __restrict *, int *,
            double __restrict *, double __restrict *, int *);

/*
SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, &
     LDC, WORK, LWORK, INFO )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DORMBR
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DORMBR
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!!*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, VECT
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DORMBR(const char *, const char *, const char *, int *, int *,
            int *, double __restrict *, int *, double __restrict *,
	    double __restrict *, int *, double __restrict *,
	    int *, int *);

/*
SUBROUTINE DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
     WORK, LWORK, INFO )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DORMLQ
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DORMLQ
#endif
       implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DORMLQ(const char *, const char *, int *, int *, int *,
            double __restrict *, int *, double __restrict *,
	    double __restrict *, int *, double __restrict *,
	    int *, int *);

/*
 SUBROUTINE DORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
      WORK, INFO )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DORML2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DORML2
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DORML2(const char *, const char *, int *, int *, int *,
            double __restrict *, int *, double __restrict *,
	    double __restrict *, int *, double __restrict *,
	    int *);

/*
SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
     WORK, LWORK, INFO )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DORMQR
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DORMQR
#endif
  implicit none
  
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!*     ..
!*!     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
          DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DORMQR(const char *, const char *, int *, int *, int *,
            double __restrict *, int *, double __restrict *,
	    double __restrict *, int *, double __restrict *,
	    int *, int *);

/*
 SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
      WORK, INFO )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DORM2R
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DORM2R
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
        DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: TAU
          DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DORM2R(const char *, const char *, int *, int *, int *,
            double __restrict *, int *, double __restrict *,
	    double __restrict *, int *, double __restrict *,
	    int *);

/*
SUBROUTINE DGESVXX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
                         EQUED, R, C, B, LDB, X, LDX, RCOND, RPVGRW, &
                         BERR, N_ERR_BNDS, ERR_BNDS_NORM, &
                         ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, &
                         INFO )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGESVXX
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGESVXX
#endif
       implicit none
!*
!*  -- LAPACK driver routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     April 2012
!*
!*     .. Scalar Arguments ..
      CHARACTER          EQUED, FACT, TRANS
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, &
                        N_ERR_BNDS
      DOUBLE PRECISION   RCOND, RPVGRW
!*     ..
!*     .. Array Arguments ..
      !INTEGER            IPIV( * ), IWORK( * )
      !DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ), &
      !                  X( LDX , * ),WORK( * )
      !DOUBLE PRECISION   R( * ), C( * ), PARAMS( * ), BERR( * ), &
      !                  ERR_BNDS_NORM( NRHS, * ), &
      !                  ERR_BNDS_COMP( NRHS, * )
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AF
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: R
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PARAMS
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BERR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ERR_BNDS_NORM
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ERR_BNDS_COMP
*/
void DGESVXX(const char *, const char *, int *, int *, double __restrict *,
             int *, double __restrict *, int *, int __restrict *, const char *,
	     double __restrict *, double __restrict *, double __restrict *,
	     int *, double __restrict *, int *, double *, double *,
	     double __restrict *, int *, double __restrict *, double __restrict *,
	     int *, double __restrict *, double __restrict *, int __restrict *,
	     int *);

	    
}











#endif


















#endif /*__GMS_LAPACK_C_IFACE_H__*/
