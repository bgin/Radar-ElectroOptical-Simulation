
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

/*
 SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, &
      WORK, LWORK, IWORK, INFO )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGELSD
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGELSD
#endif
   implicit none
   
!*
!*  -- LAPACK driver routine (version 3.7.1) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2017
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
      DOUBLE PRECISION   RCOND
!*     ..
!*     .. Array Arguments ..
      !INTEGER            IWORK( * )
      !DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: S
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
*/
void DGELSD(int *, int *, int *, double __restrict *, int *,
            double __restrict *, int *, double __restrict *,
	    double *, int *, double __restrict *, int *,
	    int __restrict *, int *);

/*
    SUBROUTINE DLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, &
           RANK, WORK, IWORK, INFO )
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLALSD
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLALSD
#endif
        implicit none
        
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS, RANK, SMLSIZ
      DOUBLE PRECISION   RCOND
!*     ..
!*     .. Array Arguments ..
      !INTEGER            IWORK( * )
      !DOUBLE PRECISION   B( LDB, * ), D( * ), E( * ), WORK( * )
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
*/
void DLALSD(const char *, int *, int *, int *, double __restrict *,
            double __restrict *, double __restrict *, int *, double *,
	    int *, double __restrict *, int __restrict *, int *);

/*
SUBROUTINE DLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, &
                        LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, &
                        GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK, &
                        IWORK, INFO )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLALSA
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLALSA
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.1) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2017
!*
!1*     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS, &
                        SMLSIZ
!*     ..
!*     .. Array Arguments ..
      !INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ),
     !$                   K( * ), PERM( LDGCOL, * )
     ! DOUBLE PRECISION   B( LDB, * ), BX( LDBX, * ), C( * ),
    ! $                   DIFL( LDU, * ), DIFR( LDU, * ),
    ! $                   GIVNUM( LDU, * ), POLES( LDU, * ), S( * ),
    ! $                   U( LDU, * ), VT( LDU, * ), WORK( * ),
      ! $                   Z( LDU, * )
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: GIVCOL
      INTEGER, DIMENSION(:), ALLOCATABLE :: GIVPTR
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      INTEGER, DIMENSION(:), ALLOCATABLE :: K
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PERM
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BX
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DIFL
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DIFR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GIVNUM
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: POLES
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: S
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: VT
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Z
*/
void DLALSA(int *, int *, int *, int *, double __restrict *, int *,
            double __restrict *, int *, double __restrict *, int *,
	    double __restrict *, int __restrict *, double __restrict *,
	    double __restrict *, double __restrict *, double __restrict *,
	    int __restrict *, int __restrict *, int *, int __restrict *,
	    double __restrict *, double __restrict *, double __restrict *,
	    double __restrict *, int __restrict *, int *);

/*
  SUBROUTINE DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, &
                        PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, &
                        POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO )                     
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLALS0
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLALS0
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
      INTEGER            GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL, &
                        LDGNUM, NL, NR, NRHS, SQRE
      DOUBLE PRECISION   C, S
!*     ..
!*     .. Array Arguments ..
      !INTEGER            GIVCOL( LDGCOL, * ), PERM( * )
      !DOUBLE PRECISION   B( LDB, * ), BX( LDBX, * ), DIFL( * ),
     !$                   DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ),
      !$                   POLES( LDGNUM, * ), WORK( * ), Z( * )
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: GIVCOL
      INTEGER, DIMENSION(:), ALLOCATABLE :: PERM
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BX
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DIFL
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DIFR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GIVNUM
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: POLES
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z
*/
void DLALS0(int *, int *, int *, int *, int *, double __restrict *,
            int *, double __restrict *, int *, int __restrict *,
	    int *, int __restrict *, int *, double __restrict *,
	    int *, double __restrict *, double __restrict *,
	    double __restrict *, double __restrict *, int *,
	    double *, double *, double __restrict *, int *);

/*
 SUBROUTINE DLASDA( ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K, &
                        DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, &
                        PERM, GIVNUM, C, S, WORK, IWORK, INFO )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASDA
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASDA
#endif
   implicit none
   
!*
!*  -- LAPACK auxiliary routine (version 3.7.1) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2017
!*
!*     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, LDGCOL, LDU, N, SMLSIZ, SQRE
!*     ..
!*     .. Array Arguments ..
     ! INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ),
     !$                   K( * ), PERM( LDGCOL, * )
     ! DOUBLE PRECISION   C( * ), D( * ), DIFL( LDU, * ), DIFR( LDU, * ),
     !$                   E( * ), GIVNUM( LDU, * ), POLES( LDU, * ),
     !$                   S( * ), U( LDU, * ), VT( LDU, * ), WORK( * ),
      !$                   Z( LDU, * )
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: GIVCOL
      INTEGER, DIMENSION(:), ALLOCATABLE :: GIVPTR
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      INTEGER, DIMENSION(:), ALLOCATABLE :: K
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PERM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DIFL
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DIFR
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GIVNUM
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: POLES
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: S
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: VT
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Z
*/
void DLASDA(int *, int *, int *, int *, double __restrict *,
            double __restrict *, double __restrict *, int *,
	    double __restrict *, int __restrict *,
	    double __restrict *, double __restrict *,
	    double __restrict *, double __restrict *,
	    int __restrict *, int __restrict *, int *,
	    int __restrict *, double __restrict *,
	    double __restrict *, double __restrict *,
	    double __restrict *, int __restrict *, int *);

/*
 SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX )
     !DIR$ ATTRIBUTES FORCEINLINE :: DLAMRG
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLAMRG
    !DIR$ OPTIMIZE : 3
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            DTRD1, DTRD2, N1, N2
!*     ..
!*     .. Array Arguments ..
      INTEGER            INDEX( * )
      !DOUBLE PRECISION   A( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
*/
void DLAMRG(int *, int *, double __restrict *, int *,
            int *, int*);

/*
 SUBROUTINE DLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB )
     !DIR$ ATTRIBUTES FORCEINLINE :: DLASDT
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASDT
    !DIR$ OPTIMIZE : 3
#endif
      implicit none
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            LVL, MSUB, N, ND
!*     ..
!*     .. Array Arguments ..
      !INTEGER            INODE( * ), NDIML( * ), NDIMR( * )
      INTEGER, DIMENSION(:), ALLOCATABLE :: INODE,NDIML,NDIMR
*/
void DLASDT(int *, int *, int *, int __restrict *,
            int __restrict *, int __restrict *, int *);

/*
SUBROUTINE DGEEQUB( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, &
     INFO )
   !DIR$ ATTRIBUTES FORCEINLINE :: DGEEQUB
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGEEQUB
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DGEEQUB
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
      INTEGER            INFO, LDA, M, N
      DOUBLE PRECISION   AMAX, COLCND, ROWCND
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION   A( LDA, * ), C( * ), R( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: R
*/
void DGEEQUB(int *, int *, double __restrict *, int *,
             double __restrict *, double __restrict *,
	     double *, double *, double *, int *);

/*
 SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
     !DIR$ ATTRIBUTES FORCEINLINE :: DGETRS
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGETRS
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGETRS
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!*     ..
!*     .. Array Arguments ..
      !INTEGER            IPIV( * )
      !DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
*/
void DGETRS(const char *, int *, int *, double __restrict *,
            int *, int __restrict *, double __restrict *,
	    int *, int *);

/*
SUBROUTINE DLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, &
     EQUED )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLAQGE
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLAQGE
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
      CHARACTER          EQUED
      INTEGER            LDA, M, N
      DOUBLE PRECISION   AMAX, COLCND, ROWCND
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), C( * ), R( * )
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: R
*/
void DLAQGE(int *, int *, double __restrict *, int *,
            double __restrict *, double __restrict *,
	    double *, double *, double *, const char *);

/*
 SUBROUTINE DGERFSX( TRANS, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, &
                         R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, &
                         ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, &
                         WORK, IWORK, INFO )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGERFSX
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGERFSX
#endif
       implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          TRANS, EQUED
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, &
                         N_ERR_BNDS
      DOUBLE PRECISION   RCOND
!*     ..
!*     .. Array Arguments ..
     ! INTEGER            IPIV( * ), IWORK( * )
     ! DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     !$                   X( LDX , * ), WORK( * )
     ! DOUBLE PRECISION   R( * ), C( * ), PARAMS( * ), BERR( * ),
     !$                   ERR_BNDS_NORM( NRHS, * ),
      !$                   ERR_BNDS_COMP( NRHS, * )
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV, IWORK
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
void DGERFSX(const char *, const char *, int *, int *, double __restrict *,
             int *, double __restrict *,  int *, int __restrict *,
	     double __restrict *, double __restrict *, double __restrict *,
	     int *, double __restrict *, int *, double *, double __restrict *,
	     int *, double __restrict *, double __restrict *,  int *,
	     double __restrict *, double __restrict *, double __restrict *,
	     int __restrict *, int *);

/*
SUBROUTINE DLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE, N, NRHS, A, &
                                     LDA, AF, LDAF, IPIV, COLEQU, C, B, &
                                     LDB, Y, LDY, BERR_OUT, N_NORMS, &
                                     ERRS_N, ERRS_C, RES, AYB, DY, &
                                     Y_TAIL, RCOND, ITHRESH, RTHRESH, &
                                     DZ_UB, IGNORE_CWISE, INFO )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLA_DGERFSX_EXTENDED
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLA_DGERFSX_EXTENDED
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.1) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2017
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE, &
                        TRANS_TYPE, N_NORMS, ITHRESH
      LOGICAL            COLEQU, IGNORE_CWISE
      DOUBLE PRECISION   RTHRESH, DZ_UB, RCOND
!*     ..
!*     .. Array Arguments ..
      !INTEGER            IPIV( * )
      !DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     !$                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )
      !DOUBLE PRECISION   C( * ), AYB( * ), RCOND, BERR_OUT( * ),
      !$                   ERRS_N( NRHS, * ), ERRS_C( NRHS, * )
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AF
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Y
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RES
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DY
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Y_TAIL
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AYB
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BERR_OUT
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ERRS_N
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ERRS_C
*/
void DLA_GERFSX_EXTENDED(int *, int *, int *, int *, double __restrict *,
                         int *, double __restrict *, int *, int __restrict *,
			 int *, double __restrict *, double __restrict *,
			 int *, double __restrict , int *, double __restrict *,
			 int *, double __restrict *, double __restrict *,
			 double __restrict *, double __restrict *,
			 double __restrict *, double __restrict *, double *,
			 int *, double *, double *, int *, int *);

/*
SUBROUTINE DLA_GEAMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, &
      Y, INCY )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLA_GEAMV
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLA_GEAMV
#endif
   use omp_lib
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.1) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2017
!*
!*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N, TRANS
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: X
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: Y
*/
void DLA_GEAMV(int *, int *, int *, double *, double __restrict *,
               int *, double __restrict *, int *, double *,
	       double __restrict *, int *);

/*
 SUBROUTINE DPOTRF ( UPLO, N, A, LDA, INFO )
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DPOTRF
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DPOTRF
#emdif
       implicit none
!*
!*  -- LAPACK computational routine (version 3.1) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
*/
void DPOTRF(const char *, int *, double __restrict *, int *, int *);



}











#endif


















#endif /*__GMS_LAPACK_C_IFACE_H__*/
