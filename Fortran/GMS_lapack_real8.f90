



!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!!*> \date December 2016
!*
!*> \ingroup auxOTHERcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DDISNA( JOB, M, N, D, SEP, INFO) !GCC$ ATTRIBUTES inline :: DDISNA !GCC$ ATTRIBUTES aligned(32) :: DDISNA
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DDISNA( JOB, M, N, D, SEP, INFO)
      !DIR$ ATTRIBUTES FORCEINLINE :: DDISNA
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DDISNA
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DDISNA
#endif
   
      implicit none

!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            INFO, M, N
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), SEP( * )
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            DECR, EIGEN, INCR, LEFT, RIGHT, SING
      INTEGER            I, K
      DOUBLE PRECISION   ANORM, EPS, NEWGAP, OLDGAP, SAFMIN, THRESH
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           XERBLA
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      EIGEN = LSAME( JOB, 'E' )
      LEFT = LSAME( JOB, 'L' )
      RIGHT = LSAME( JOB, 'R' )
      SING = LEFT .OR. RIGHT
      IF( EIGEN ) THEN
         K = M
      ELSE IF( SING ) THEN
         K = MIN( M, N )
      END IF
      IF( .NOT.EIGEN .AND. .NOT.SING ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( K.LT.0 ) THEN
         INFO = -3
      ELSE
         INCR = .TRUE.
         DECR = .TRUE.
         DO 10 I = 1, K - 1
            IF( INCR ) &
              INCR = INCR .AND. D( I ).LE.D( I+1 )
            IF( DECR ) &
               DECR = DECR .AND. D( I ).GE.D( I+1 )
   10    CONTINUE
         IF( SING .AND. K.GT.0 ) THEN
            IF( INCR ) &
              INCR = INCR .AND. ZERO.LE.D( 1 )
            IF( DECR ) &
               DECR = DECR .AND. D( K ).GE.ZERO
         END IF
         IF( .NOT.( INCR .OR. DECR ) ) &
            INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
        ! CALL XERBLA( 'DDISNA', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      IF( K.EQ.0 ) &
         RETURN
!*
!*     Compute reciprocal condition numbers
!*
      IF( K.EQ.1 ) THEN
         SEP( 1 ) = DLAMCH( 'O' )
      ELSE
         OLDGAP = ABS( D( 2 )-D( 1 ) )
         SEP( 1 ) = OLDGAP
         DO 20 I = 2, K - 1
            NEWGAP = ABS( D( I+1 )-D( I ) )
            SEP( I ) = MIN( OLDGAP, NEWGAP )
            OLDGAP = NEWGAP
   20    CONTINUE
         SEP( K ) = OLDGAP
      END IF
      IF( SING ) THEN
         IF( ( LEFT .AND. M.GT.N ) .OR. ( RIGHT .AND. M.LT.N ) ) THEN
            IF( INCR ) &
               SEP( 1 ) = MIN( SEP( 1 ), D( 1 ) )
            IF( DECR ) &
               SEP( K ) = MIN( SEP( K ), D( K ) )
         END IF
      END IF
!*
!*     Ensure that reciprocal condition numbers are not less than
!*     threshold, in order to limit the size of the error bound
!*
      EPS = DLAMCH( 'E' )
      SAFMIN = DLAMCH( 'S' )
      ANORM = MAX( ABS( D( 1 ) ), ABS( D( K ) ) )
      IF( ANORM.EQ.ZERO ) THEN
         THRESH = EPS
      ELSE
         THRESH = MAX( EPS*ANORM, SAFMIN )
      END IF
      DO 30 I = 1, K
         SEP( I ) = MAX( SEP( I ), THRESH )
   30 CONTINUE

END SUBROUTINE


!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleGEauxiliary
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK ) !GCC$ ATTRIBUTES HOT :: DLANGE !GCC$ ATTRIBUTES aligned(32) :: DLANGE !GCC$ ATTRIBUTES no_stack_protector :: DLANGE
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLANGE
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLANGE
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
      CHARACTER          NORM
      INTEGER            LDA, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
!*     ..
!*
!* =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SUM, VALUE, TEMP
!*     ..
!*     .. Local Arrays ..
      DOUBLE PRECISION   SSQ( 2 ), COLSSQ( 2 )
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           DLASSQ, DCOMBSSQ
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN, SQRT
!*     ..
!*     .. Executable Statements ..
!*
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!*
!*        Find max(abs(A(i,j))).
!*
         VALUE = ZERO
         
         DO 20 J = 1, N
            DO 10 I = 1, M
               TEMP = ABS( A( I, J ) )
               IF( VALUE.LT.TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!*
!*        Find norm1(A).
!*
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            !$OMP SIMD ALIGNED(A:64) LINEAR(I:1) REDUCTION(+:SUM) UNROLL PARTIAL(6)
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            IF( VALUE.LT.SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!*
!*        Find normI(A).
         !*
          !$OMP SIMD ALIGNED(A:64) LINEAR(I:1)  UNROLL PARTIAL(8)
         DO 50 I = 1, M
            WORK( I ) = ZERO
50       CONTINUE
         !$OMP PARALLEL DO SCHEDULE(STATIC,10) DEFAULT(NONE) SHARED(WORK,A,N,M) PRIVATE(J,I) COLLAPSE(2)    
         DO 70 J = 1, N
                !$OMP SIMD ALIGNED(A:64) LINEAR(I:1)  UNROLL PARTIAL(6)
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
70       CONTINUE
         !$OMP END PARALLEL DO      
         VALUE = ZERO
         !$OMP PARALLEL DO SCHEDULE(STATIC,10) DEFAULT(NONE) SHARED(M,WORK,VALUE) PRIVATE(I,TEMP) IF(M>=1000)     
         DO 80 I = 1, M
            TEMP = WORK( I )
            !$OMP CRITICAL
            IF( VALUE.LT.TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
            !$OMP END CRITICAL
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!*
!*        Find normF(A).
!*        SSQ(1) is scale
!*        SSQ(2) is sum-of-squares
!*        For better accuracy, sum each column separately.
!*
         SSQ( 1 ) = ZERO
         SSQ( 2 ) = ONE
         DO 90 J = 1, N
            COLSSQ( 1 ) = ZERO
            COLSSQ( 2 ) = ONE
            CALL DLASSQ( M, A( 1, J ), 1, COLSSQ( 1 ), COLSSQ( 2 ) )
            CALL DCOMBSSQ( SSQ, COLSSQ )
   90    CONTINUE
         VALUE = SSQ( 1 )*SQRT( SSQ( 2 ) )
      END IF
      DLANGE = VALUE
    
END FUNCTION


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup OTHERauxiliary
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ ) !GCC$ ATTRIBUTES INLINE :: DLASSQ !GCC$ ATTRIBUTES aligned(32) :: DLASSQ
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
     !DIR$ ATTRIBUTES FORCEINLINE :: DLASSQ
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASSQ
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASSQ
#endif
      implicit none
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!*     ..
!*
!* =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
!*     ..
!*     .. External Functions ..
      LOGICAL            DISNAN
      EXTERNAL           DISNAN
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS
!*     ..
!*     .. Executable Statements ..
!*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            ABSXI = ABS( X( IX ) )
            IF( ABSXI.GT.ZERO.OR.DISNAN( ABSXI ) ) THEN
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
END SUBROUTINE


!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*!> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2018
!*
!*> \ingroup OTHERauxiliary
!*
!*  =====================================================================

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DCOMBSSQ( V1, V2 ) !GCC$ ATTRIBUTES INLINE :: DCOMBSSQ !GCC$ ATTRIBUTES aligned(32) :: DCOMBSSQ
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE DCOMBSSQ(V1,V2)
    !DIR$ ATTRIBUTES FORCEINLINE :: DCOMBSSQ
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DCOMBSSQ
    !DIR$ OPTIMIZE : 3
#endif
      implicit none
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2018
!*
!*     .. Array Arguments ..
      DOUBLE PRECISION   V1( 2 ), V2( 2 )
!*     ..
!*
!* =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!*     ..
!*     .. Executable Statements ..
!*
      IF( V1( 1 ).GE.V2( 1 ) ) THEN
         IF( V1( 1 ).NE.ZERO ) THEN
            V1( 2 ) = V1( 2 ) + ( V2( 1 ) / V1( 1 ) )**2 * V2( 2 )
         END IF
      ELSE
         V1( 2 ) = V2( 2 ) + ( V1( 1 ) / V2( 1 ) )**2 * V1( 2 )
         V1( 1 ) = V2( 1 )
      END IF
    
END SUBROUTINE 


!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup OTHERauxiliary
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )  !GCC$ ATTRIBUTES hot :: DLACPY !GCC$ ATTRIBUTES aligned(32) :: DLACPY !GCC$ ATTRIBUTES no_stack_protector :: DLACPY
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
      !DIR$ ATTRIBUTES INLINE :: DLACPY
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLACPY
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLACPY
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
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER            I, J
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MIN
!*     ..
*!     .. Executable Statements ..
!*
      IF( LSAME( UPLO, 'U' ) ) THEN
         !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(B,A,N,M) PRIVATE(J,I)
         DO 20 J = 1, N
            !$OMP SIMD ALIGNED(A:64,B) LINEAR(I:1) UNROLL PARTIAL(8)
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
20       CONTINUE
          !$OMP END PARALLEL DO     
     ELSE IF( LSAME( UPLO, 'L' ) ) THEN
            !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(B,A) COLLAPSE(2) PRIVATE(J,I)    
        DO 40 J = 1, N
              !$OMP SIMD ALIGNED(A:64,B) LINEAR(I:1) UNROLL PARTIAL(8)
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
40       CONTINUE
         !$OMP END PARALLEL DO      
     ELSE
           !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(B,A,N,M) COLLAPSE(2) PRIVATE(J,I)       
        DO 60 J = 1, N
             !$OMP SIMD ALIGNED(A:64,B) LINEAR(I:1) UNROLL PARTIAL(8)
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
60      CONTINUE
          !$OMP END PARALLEL DO     
      END IF
    
END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleOTHERauxiliary
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLAPMT( FORWRD, M, N, X, LDX, K )  !GCC$ ATTRIBUTES hot :: DLAPMT !GCC$ ATTRIBUTES aligned(32) :: DLAPMT !GCC$ ATTRIBUTES no_stack_protector :: DLAPMT
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE DLAPMT(FORWRD, M, N, X, LDX, K)
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLAPMT
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLAPMT
#endif
      implicit none
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      LOGICAL            FORWRD
      INTEGER            LDX, M, N
!*     ..
!*!     .. Array Arguments ..
      INTEGER            K( * )
      DOUBLE PRECISION   X( LDX, * )
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER            I, II, IN, J
      DOUBLE PRECISION   TEMP
!*     ..
!*     .. Executable Statements ..
!*
    !  IF( N.LE.1 )
    ! $   RETURN
!*
      DO 10 I = 1, N
         K( I ) = -K( I )
   10 CONTINUE
!*
      IF( FORWRD ) THEN
!*
!*        Forward permutation
!*
         DO 50 I = 1, N
!*
            IF( K( I ).GT.0 ) &
              GO TO 40
!*
            J = I
            K( J ) = -K( J )
            IN = K( J )
!*
   20       CONTINUE
            IF( K( IN ).GT.0 ) &
              GO TO 40
!*
            DO 30 II = 1, M
               TEMP = X( II, J )
               X( II, J ) = X( II, IN )
               X( II, IN ) = TEMP
   30       CONTINUE
!*
            K( IN ) = -K( IN )
            J = IN
            IN = K( IN )
            GO TO 20
!*!
   40       CONTINUE
!*
   50    CONTINUE
!*
      ELSE
!*
!*        Backward permutation
!*
         DO 90 I = 1, N
!*
            IF( K( I ).GT.0 ) &
               GO TO 80
!*
            K( I ) = -K( I )
            J = K( I )
   60       CONTINUE
            IF( J.EQ.I ) &
               GO TO 80
!*
            DO 70 II = 1, M
               TEMP = X( II, I )
               X( II, I ) = X( II, J )
               X( II, J ) = TEMP
   70       CONTINUE
!*
            K( J ) = -K( J )
            J = K( J )
            GO TO 60
!*
   80       CONTINUE
!*
   90    CONTINUE
!*
      END IF
!*
END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2017
!*
!*> \ingroup doubleOTHERauxiliary
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU ) !GCC$ ATTRIBUTES hot :: DLARFG !GCC$ ATTRIBUTES aligned(32) :: DLARFG !GCC$ ATTRIBUTES no_stack_protector :: DLARFG
#elif defined(__ICC) || defined(__INTEL_COMPILER)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DSCAL
!*     ..
!*     .. Executable Statements ..
!*
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
!*
      XNORM = DNRM2( N-1, X, INCX )
!*
      IF( XNORM.EQ.ZERO ) THEN
!*
!*        H  =  I
!*
         TAU = ZERO
      ELSE
!*
!*        general case
!*
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!*
!*           XNORM, BETA may be inaccurate; scale X and recompute them
!*
            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) ) &
               GO TO 10
!*
!*           New BETA is at most 1, at least SAFMIN
!*
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         END IF
         TAU = ( BETA-ALPHA ) / BETA
         CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
!*
!*        If ALPHA is subnormal, it may lose relative accuracy
!*
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
END SUBROUTINE


!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup OTHERauxiliary
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )  !GCC$ ATTRIBUTES hot :: DLASET !GCC$ ATTRIBUTES aligned(32) :: DLASET !GCC$ ATTRIBUTES no_stack_protector :: DLASET
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE DLASET(UPLO, M, N, ALPHA, BETA, A, LDA)
  !DIR$ ATTRIBUTES FORCEINLINE :: DLASET
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASET
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLASET
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
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
!*     ..
!*     .. Array Arguments ..
     ! DOUBLE PRECISION   A( LDA, * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
!*     ..
!*
!* =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER            I, J
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MIN
!*!     ..
!*     .. Executable Statements ..
!*
      IF( LSAME( UPLO, 'U' ) ) THEN
!*
!*        Set the strictly upper triangular or trapezoidal part of the
!*        array to ALPHA.
         !*
         !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(A,N,M,ALPHA) PRIVATE(J,I)
         DO 20 J = 2, N
            !$OMP SIMD ALIGNED(A:64) LINEAR(I:1) UNROLL PARTIAL(8)
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
20      CONTINUE
        !$OMP END PARALLEL DO       
!*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!*
!*        Set the strictly lower triangular or trapezoidal part of the
!*        array to ALPHA.
         !*
          !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(A,M,N,ALPHA) PRIVATE(J,I) 
         DO 40 J = 1, MIN( M, N )
             !$OMP SIMD ALIGNED(A:64) LINEAR(I:1) 
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
40      CONTINUE
        !$OMP END PARALLEL DO       
!*
      ELSE
!*
!*        Set the leading m-by-n submatrix to ALPHA.
         !*
          !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(A,N,M,ALPHA) COLLAPSE(2) PRIVATE(J,I)
         DO 60 J = 1, N
              !$OMP SIMD ALIGNED(A:64) LINEAR(I:1) 
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
60       CONTINUE
         !$OMP END PARALLEL DO    
      END IF
!*
!*     Set the first min(M,N) diagonal elements to BETA.
      !*
      
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE

END SUBROUTINE


!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleOTHERcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLATZM( SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK )  !GCC$ ATTRIBUTES hot :: DLATZM !GCC$ ATTRIBUTES aligned(32) :: DLATZM !GCC$ ATTRIBUTES no_stack_protector :: DLATZM
#elif defined(__ICC) || defined(__INTEL_COMPILER)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGEMV, DGER
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MIN
!*     ..
!*     .. Executable Statements ..
!*
      IF( ( MIN( M, N ).EQ.0 ) .OR. ( TAU.EQ.ZERO ) ) &
         RETURN
!*
      IF( LSAME( SIDE, 'L' ) ) THEN
!*
!!*        w :=  (C1 + v**T * C2)**T

         CALL DCOPY( N, C1, LDC, WORK, 1 )
         CALL DGEMV( 'Transpose', M-1, N, ONE, C2, LDC, V, INCV, ONE, &
                     WORK, 1 )
!*
!*        [ C1 ] := [ C1 ] - tau* [ 1 ] * w**T
!*        [ C2 ]    [ C2 ]        [ v ]
!*
         CALL DAXPY( N, -TAU, WORK, 1, C1, LDC )
         CALL DGER( M-1, N, -TAU, V, INCV, WORK, 1, C2, LDC )
!*
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!*
!*        w := C1 + C2 * v
!*
         CALL DCOPY( M, C1, 1, WORK, 1 )
         CALL DGEMV( 'No transpose', M, N-1, ONE, C2, LDC, V, INCV, ONE, &
                     WORK, 1 )
!*
!*        [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**T]
!*
         CALL DAXPY( M, -TAU, WORK, 1, C1, 1 )
         CALL DGER( M, N-1, -TAU, WORK, 1, V, INCV, C2, LDC )
      END IF
!*

END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date June 2013
!*
!*> \ingroup doubleOTHERauxiliary
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  The shape of the matrix V and the storage of the vectors which define
!*>  the H(i) is best illustrated by the following example with n = 5 and
!*>  k = 3. The elements equal to 1 are not stored; the corresponding
!*>  array elements are modified but restored on exit. The rest of the
!*>  array is not used.
!*>
!*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!*>
!*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!*>                   ( v1  1    )                     (     1 v2 v2 v2 )
!*>                   ( v1 v2  1 )                     (        1 v3 v3 )
!*>                   ( v1 v2 v3 )
!*>                   ( v1 v2 v3 )
!*>
!*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!*>
!*>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!*>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!*>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!*>                   (     1 v3 )
!*>                   (        1 )
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
T, LDT, C, LDC, WORK, LDWORK )  !GCC$ ATTRIBUTES hot :: DLARFB !GCC$ ATTRIBUTES aligned(32) :: DLARFB !GCC$ ATTRIBUTES no_stack_protector :: DLARFB
#elif defined(__ICC) || defined(__INTEL_COMPILER)
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
!*     ..
!*!
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
!*     ..
!*     .. Executable Statements ..
!*
!*     Quick return if possible
!*
   !   IF( M.LE.0 .OR. N.LE.0 )
  !   $   RETURN
!*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
!*
      IF( LSAME( STOREV, 'C' ) ) THEN
!*
         IF( LSAME( DIRECT, 'F' ) ) THEN
!*
!*           Let  V =  ( V1 )    (first K rows)
!*                     ( V2 )
!*           where  V1  is unit lower triangular.
!*
            IF( LSAME( SIDE, 'L' ) ) THEN
!*
!*              Form  H * C  or  H**T * C  where  C = ( C1 )
!*                                                    ( C2 )
!*
!*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!*
!*              W := C1**T
               !*
               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K) PRIVATE(J)
               DO 10 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
10             CONTINUE
               !$OMP END PARALLEL DO   
!*
!*              W := W * V1
!*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!*
!*                 W := W + C2**T * V2
!*
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                             ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                             ONE, WORK, LDWORK )
               END IF
!*
!*              W := W * T**T  or  W * T
!*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - V * W**T
!*
               IF( M.GT.K ) THEN
!!*
!*                 C2 := C2 - V2 * W**T
!*
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,   &
                             -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, &
                             C( K+1, 1 ), LDC )
               END IF
!*
!!*              W := W * V1**T
!*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                          ONE, V, LDV, WORK, LDWORK )
!*
!*              C1 := C1 - W**T
               !*
               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,N) PRIVATE(J,I) COLLAPSE(2) 
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
30            CONTINUE
              !$OMP END PARALLEL DO       
!*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!*
!*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!*
!*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!*
!*              W := C1
               !!*
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K) PRIVATE(J) 
               DO 40 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!*
!*              W := W * V1
!*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!*
!*                 W := W + C2 * V2
!*
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                             ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                             ONE, WORK, LDWORK )
               END IF
!*
!*              W := W * T  or  W * T**T
!*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - W * V**T
!*
               IF( N.GT.K ) THEN
!*
!*                 C2 := C2 - W * V2**T
!*
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,  &
                             -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                              C( 1, K+1 ), LDC )
               END IF
!*
!*              W := W * V1**T
!*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                           ONE, V, LDV, WORK, LDWORK )
!*
!*              C1 := C1 - W
               !*
                 !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,M) PRIVATE(J,I) COLLAPSE(2) 
               DO 60 J = 1, K
                  !$OMP SIMD ALIGNED(C:64,WORK) LINEAR(I:1) UNROLL PARTIAL(6)
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
60             CONTINUE
                !$OMP END PARALLEL DO     
            END IF

         ELSE
!*
!*           Let  V =  ( V1 )
!*                     ( V2 )    (last K rows)
!*           where  V2  is unit upper triangular.
!*
            IF( LSAME( SIDE, 'L' ) ) THEN
!*
!*              Form  H * C  or  H**T * C  where  C = ( C1 )
!*                                                    ( C2 )
!*
!*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!*
!*              W := C2**T
               !*
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K) PRIVATE(J)   
               DO 70 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
!*
!*              W := W * V2
!*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                          K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!*
!*                 W := W + C1**T * V1
!*
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                             ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!*
!*              W := W * T**T  or  W * T
!*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - V * W**T
!*
               IF( M.GT.K ) THEN
!*
!*                 C1 := C1 - V1 * W**T
!*
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,  &
                             -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!*
!*              W := W * V2**T
!*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                          ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
!*
!*              C2 := C2 - W**T
               !*
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,N) PRIVATE(J,I)  
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
90             CONTINUE
               !$OMP END PARALLEL DO      

            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!*
!*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!*
!*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!*
!*              W := C2
!*           !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K) PRIVATE(J)  
               DO 100 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
!!*
!1*              W := W * V2
!!*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                          K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!*
!*                 W := W + C1 * V1
!*
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,  &
                             ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!*
!*              W := W * T  or  W * T**T
!*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                          ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - W * V**T
!*
               IF( N.GT.K ) THEN
!*
!*                 C1 := C1 - W * V1**T
!*
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,  &
                             -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!*
!*              W := W * V2**T
!*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                           ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
!*
!*              C2 := C2 - W
               !*
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,M) PRIVATE(J,I)  
               DO 120 J = 1, K
                  !$OMP SIMD ALIGNED(C:64,WORK) LINEAR(I:1) UNROLL PARTIAL(6)
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
120            CONTINUE
               !$OMP END PARALLEL DO      
            END IF
         END IF
!*
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!*
         IF( LSAME( DIRECT, 'F' ) ) THEN
!*
!*           Let  V =  ( V1  V2 )    (V1: first K columns)
!*           where  V1  is unit upper triangular.
!*
            IF( LSAME( SIDE, 'L' ) ) THEN
!*
!*              Form  H * C  or  H**T * C  where  C = ( C1 )
!*                                                    ( C2 )
!*
!*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!*
!*              W := C1**T
               !*
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J)   IF(K>=200)
               DO 130 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
!*
!*              W := W * V1**T
!*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                          ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!*
!*                 W := W + C2**T * V2**T
!*
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                             C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, &
                             WORK, LDWORK )
               END IF
!*
!*              W := W * T**T  or  W * T
!*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                          ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - V**T * W**T
!*
               IF( M.GT.K ) THEN
!*
!*                 C2 := C2 - V2**T * W**T
!*
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                             V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
                             C( K+1, 1 ), LDC )
               END IF
!*
!*              W := W * V1
!*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                          K, ONE, V, LDV, WORK, LDWORK ) 
!*
!*              C1 := C1 - W**T
               !*
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,N) PRIVATE(J,I)  
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
!*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!*
!*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!*
!*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!*
!*              W := C1
               !*
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K) PRIVATE(J)  
               DO 160 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!*
!*              W := W * V1**T
!*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                          ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!*
!*                 W := W + C2 * V2**T
!*
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                             ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                             ONE, WORK, LDWORK )
               END IF
!*
!*              W := W * T  or  W * T**T
!*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                          ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - W * V
!*
               IF( N.GT.K ) THEN
!*
!*                 C2 := C2 - W * V2
!*
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                             -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                             C( 1, K+1 ), LDC )
               END IF
!*
!*              W := W * V1
!*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                          K, ONE, V, LDV, WORK, LDWORK )
!*
!*              C1 := C1 - W
               !*
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,M) PRIVATE(J,I)  
               DO 180 J = 1, K
                  !$OMP SIMD ALIGNED(C:64,WORK) LINEAR(I:1) UNROLL PARTIAL(6)
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
!*
            END IF
!*
         ELSE
!*
!*           Let  V =  ( V1  V2 )    (V2: last K columns)
!!*           where  V2  is unit lower triangular.
!*
            IF( LSAME( SIDE, 'L' ) ) THEN
!*
!*              Form  H * C  or  H**T * C  where  C = ( C1 )
!*                                                    ( C2 )
!*
!*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!*
!*              W := C2**T
               !*
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K) PRIVATE(J)  
               DO 190 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
!*
!*              W := W * V2**T
!*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                          ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!*
!*                 W := W + C1**T * V1**T
!*
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                             C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!*
!!*              W := W * T**T  or  W * T
!*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                          ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - V**T * W**T
!*
               IF( M.GT.K ) THEN
!*
!*                 C1 := C1 - V1**T * W**T
!*
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                             V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!*
!*              W := W * V2
!*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
                          K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!*
!*              C2 := C2 - W**T
                          !*
                 !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,N) PRIVATE(J,I)  
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
!*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!*
!*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!*
!*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!*
!*              W := C2
               !*
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K) PRIVATE(J)  
               DO 220 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
!*
!*              W := W * V2**T
!*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                          ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!*
!*                 W := W + C1 * V1**T
!*
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                             ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!*
!*              W := W * T  or  W * T**T
!*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                          ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - W * V
!*
               IF( N.GT.K ) THEN
!*
!*                 C1 := C1 - W * V1
!*
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                             -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!*
!*              W := W * V2
!*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!*
!*              C1 := C1 - W
               !*
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,M) PRIVATE(J,I)  
               DO 240 J = 1, K
                  !$OMP SIMD ALIGNED(C:64,WORK) LINEAR(I:1) UNROLL PARTIAL(6)
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
!*
            END IF
!*
         END IF
      END IF

END SUBROUTINE


!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleOTHERauxiliary
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK) !GCC$ ATTRIBUTES hot :: DLARF !GCC$ ATTRIBUTES aligned(32) :: DLARF !GCC$ ATTRIBUTES no_stack_protector :: DLARF
#elif defined(__ICC) || defined(__INTEL_COMPILER)
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
      
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
!!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILADLR, ILADLC
      EXTERNAL           LSAME, ILADLR, ILADLC
!*     ..
!*     .. Executable Statements ..
!*
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILADLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILADLR(M, LASTV, C, LDC)
         END IF
      END IF
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
      IF( APPLYLEFT ) THEN
!*
!*        Form  H * C
!*
         IF( LASTV.GT.0 ) THEN
!*
!*           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
!*
            CALL DGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, &
                ZERO, WORK, 1 )
!*
!*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
!*
            CALL DGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!*
!*        Form  C * H
!*
         IF( LASTV.GT.0 ) THEN
!*
!*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!*
            CALL DGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC, &
                 V, INCV, ZERO, WORK, 1 )
!*
!*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
!*
            CALL DGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
 
END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleOTHERauxiliary
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  The shape of the matrix V and the storage of the vectors which define
!*>  the H(i) is best illustrated by the following example with n = 5 and
!*>  k = 3. The elements equal to 1 are not stored.
!*>
!*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!*>
!*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!*>                   ( v1  1    )                     (     1 v2 v2 v2 )
!*>                   ( v1 v2  1 )                     (        1 v3 v3 )
!*>                   ( v1 v2 v3 )
!*>                   ( v1 v2 v3 )
!*>
!*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!*>
!*>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!*>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!*>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!*>                   (     1 v3 )
!*>                   (        1 )
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )  !GCC$ ATTRIBUTES hot :: DLARFT !GCC$ ATTRIBUTES aligned(32) :: DLARFT !GCC$ ATTRIBUTES no_stack_protector :: DLARFT
#elif defined(__ICC) || defined(__INTEL_COMPILER)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. Executable Statements ..
!*
!*     Quick return if possible
!*
!      IF( N.EQ.0 )
!     $   RETURN
!*
      IF( LSAME( DIRECT, 'F' ) ) THEN
         PREVLASTV = N
         DO I = 1, K
            PREVLASTV = MAX( I, PREVLASTV )
            IF( TAU( I ).EQ.ZERO ) THEN
!*
!*              H(i)  =  I
               !*
               !$OMP SIMD ALIGNED(T:64) LINEAR(J:1) UNROLL PARTIAL(8)
               DO J = 1, I
                  T( J, I ) = ZERO
               END DO
            ELSE
!*
!*              general case
!*
               IF( LSAME( STOREV, 'C' ) ) THEN
!*                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( LASTV, I ).NE.ZERO ) EXIT
                  END DO
                  
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( I , J )
                  END DO
                  J = MIN( LASTV, PREVLASTV )
!*
!*                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
!*
                  CALL DGEMV( 'Transpose', J-I, I-1, -TAU( I ),  &
                            V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE, &
                             T( 1, I ), 1 )
               ELSE
!*                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( I, LASTV ).NE.ZERO ) EXIT
                  END DO
                   !$OMP SIMD ALIGNED(T:64,TAU,V) LINEAR(J:1) UNROLL PARTIAL(6)
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( J , I )
                  END DO
                  J = MIN( LASTV, PREVLASTV )
!*
!*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T
!*
                  CALL DGEMV( 'No transpose', I-1, J-I, -TAU( I ), &
                             V( 1, I+1 ), LDV, V( I, I+1 ), LDV, ONE, &
                             T( 1, I ), 1 )
               END IF
!*
!*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!*
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
                          LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
               IF( I.GT.1 ) THEN
                  PREVLASTV = MAX( PREVLASTV, LASTV )
               ELSE
                  PREVLASTV = LASTV
               END IF
            END IF
         END DO
      ELSE
         PREVLASTV = 1
         DO I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
!*
!*              H(i)  =  I
               !*
                !$OMP SIMD ALIGNED(T:64) LINEAR(J:1) UNROLL PARTIAL(8)
               DO J = I, K
                  T( J, I ) = ZERO
               END DO
            ELSE
!*
!*              general case
!*
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
!*                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( LASTV, I ).NE.ZERO ) EXIT
                     END DO
                      !$OMP SIMD ALIGNED(T:64,TAU,V) LINEAR(J:1) UNROLL PARTIAL(6)
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( N-K+I , J )
                     END DO
                     J = MAX( LASTV, PREVLASTV )
!*
!*                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)
!*
                     CALL DGEMV( 'Transpose', N-K+I-J, K-I, -TAU( I ),
                                V( J, I+1 ), LDV, V( J, I ), 1, ONE,
                                T( I+1, I ), 1 )
                  ELSE
!*                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( I, LASTV ).NE.ZERO ) EXIT
                     END DO
                      !$OMP SIMD ALIGNED(T:64,TAU,V) LINEAR(J:1) UNROLL PARTIAL(6)
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( J, N-K+I )
                     END DO
                     J = MAX( LASTV, PREVLASTV )
!*
!*                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T
!*
                     CALL DGEMV( 'No transpose', K-I, N-K+I-J, &
                         -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV, &
                         ONE, T( I+1, I ), 1 )
                  END IF
!*
!*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!*
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                             T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                  IF( I.GT.1 ) THEN
                     PREVLASTV = MIN( PREVLASTV, LASTV )
                  ELSE
                     PREVLASTV = LASTV
                  END IF
               END IF
               T( I, I ) = TAU( I )
            END IF
         END DO
      END IF
  
END SUBROUTINE

    
#if 0
*>
*> \verbatim
*>
*> DGETRF computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the right-looking Level 3 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*>                has been completed, but the factor U is exactly
*>                singular, and division by zero will occur if it is used
*>                to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup doubleGEcomputational
*
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO ) !GCC$ ATTRIBUTES hot :: DGETRF !GCC$ ATTRIBUTES aligned(32) :: DGETRF
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*!
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLASWP, DTRSM
!!*     ..
!*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
        ! CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
     ! IF( M.EQ.0 .OR. N.EQ.0 )
     !   RETURN
!*
!*     Determine the block size for this environment.
!*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!*
!*        Use unblocked code.
!*
         CALL DGETRF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!*
!*        Use blocked code.
!*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!*
!*           Factor diagonal and subdiagonal blocks and test for exact
!*           singularity.
!*
            CALL DGETRF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!*
!*           Adjust INFO and the pivot indices.
!*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 ) &
               INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!*
!*           Apply interchanges to columns 1:J-1.
!*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!*
            IF( J+JB.LE.N ) THEN
!*
!*              Apply interchanges to columns J+JB:N.
!*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, &
                           IPIV, 1 )
!*
!*              Compute block row of U.
!*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                          N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
                          LDA )
               IF( J+JB.LE.M ) THEN
!*
!*                 Update trailing submatrix.
!*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                             N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
                             A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), &
                             LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
 
END SUBROUTINE

#if 0  
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*>                has been completed, but the factor U is exactly
*>                singular, and division by zero will occur if it is used
*>                to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2016
*
*> \ingroup doubleGEcomputational
*
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO ) !GCC$ ATTRIBUTES hot :: DGETRF2 !GCC$ ATTRIBUTES aligned(32) :: DGETRF2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGETRF2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGETRF2
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!*     ..
!*     .. Array Arguments ..
      INTEGER            IPIV( * )
      !  DOUBLE PRECISION   A( LDA, * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN, TEMP
      INTEGER            I, IINFO, N1, N2
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IDAMAX
      EXTERNAL           DLAMCH, IDAMAX
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DGEMM, DSCAL, DLASWP, DTRSM
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters
!*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
        ! CALL XERBLA( 'DGETRF2', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
  !    IF( M.EQ.0 .OR. N.EQ.0 )
  !   $   RETURN

      IF ( M.EQ.1 ) THEN
!*
!*        Use unblocked code for one row case
!*        Just need to handle IPIV and INFO
!*
         IPIV( 1 ) = 1
         IF ( A(1,1).EQ.ZERO ) &
           INFO = 1
!*
      ELSE IF( N.EQ.1 ) THEN
!*
!*        Use unblocked code for one column case
!*
!*
!*        Compute machine safe minimum
!*
         SFMIN = DLAMCH('S')
!*
!*        Find pivot and test for singularity
!*
         I = IDAMAX( M, A( 1, 1 ), 1 )
         IPIV( 1 ) = I
         IF( A( I, 1 ).NE.ZERO ) THEN
!*
!*           Apply the interchange
!*
            IF( I.NE.1 ) THEN
               TEMP = A( 1, 1 )
               A( 1, 1 ) = A( I, 1 )
               A( I, 1 ) = TEMP
            END IF
!*
!*           Compute elements 2:M of the column
!*
            IF( ABS(A( 1, 1 )) .GE. SFMIN ) THEN
               CALL DSCAL( M-1, ONE / A( 1, 1 ), A( 2, 1 ), 1 )
            ELSE
               DO 10 I = 1, M-1
                  A( 1+I, 1 ) = A( 1+I, 1 ) / A( 1, 1 )
   10          CONTINUE
            END IF
!*
         ELSE
            INFO = 1
         END IF
!*
      ELSE
!*
!*        Use recursive code
!*
         N1 = MIN( M, N ) / 2
         N2 = N-N1
!*
!*               [ A11 ]
!*        Factor [ --- ]
!*               [ A21 ]
!*
         CALL DGETRF2( M, N1, A, LDA, IPIV, IINFO )

         IF ( INFO.EQ.0 .AND. IINFO.GT.0 ) &
            INFO = IINFO
!*
!*                              [ A12 ]
!*        Apply interchanges to [ --- ]
!*                              [ A22 ]
!*
         CALL DLASWP( N2, A( 1, N1+1 ), LDA, 1, N1, IPIV, 1 )
!*
!*        Solve A12
!*
         CALL DTRSM( 'L', 'L', 'N', 'U', N1, N2, ONE, A, LDA, &
                     A( 1, N1+1 ), LDA )
!*
!*        Update A22
!*
         CALL DGEMM( 'N', 'N', M-N1, N2, N1, -ONE, A( N1+1, 1 ), LDA, &
                     A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA )
!*
!*        Factor A22
!*
         CALL DGETRF2( M-N1, N2, A( N1+1, N1+1 ), LDA, IPIV( N1+1 ), &
                       IINFO )
!*
!*        Adjust INFO and the pivot indices
!*
         IF ( INFO.EQ.0 .AND. IINFO.GT.0 ) &
           INFO = IINFO + N1
         DO 20 I = N1+1, MIN( M, N )
            IPIV( I ) = IPIV( I ) + N1
   20    CONTINUE
!*
!*        Apply interchanges to A21
!*
         CALL DLASWP( N1, A( 1, 1 ), LDA, N1+1, MIN( M, N), IPIV, 1 )
!*
      END IF
 
END SUBROUTINE


!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, INFO) !GCC$ ATTRIBUTES hot :: DGEBRD !GCC$ ATTRIBUTES aligned(32) :: DGEBRD
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
      
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LDWRKX, LDWRKY, LWKOPT, MINMN, NB, &
                         NBMIN, NX, WS
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DGEMM
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
!*     ..
!*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters
!*
      INFO = 0
      NB = MAX( 1, ILAENV( 1, 'DGEBRD', ' ', M, N, -1, -1 ) )
      LWKOPT = ( M+N )*NB
      WORK( 1 ) = DBLE( LWKOPT )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, M, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
      IF( INFO.LT.0 ) THEN
         !CALL XERBLA( 'DGEBRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      MINMN = MIN( M, N )
      IF( MINMN.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!*
      WS = MAX( M, N )
      LDWRKX = M
      LDWRKY = N
!*
      IF( NB.GT.1 .AND. NB.LT.MINMN ) THEN
!*
!*        Set the crossover point NX.
!*
         NX = MAX( NB, ILAENV( 3, 'DGEBRD', ' ', M, N, -1, -1 ) )
!*
!*        Determine when to switch from blocked to unblocked code.
!*
         IF( NX.LT.MINMN ) THEN
            WS = ( M+N )*NB
            IF( LWORK.LT.WS ) THEN
!*
!*              Not enough work space for the optimal NB, consider using
!*              a smaller block size.
!*
               NBMIN = ILAENV( 2, 'DGEBRD', ' ', M, N, -1, -1 )
               IF( LWORK.GE.( M+N )*NBMIN ) THEN
                  NB = LWORK / ( M+N )
               ELSE
                  NB = 1
                  NX = MINMN
               END IF
            END IF
         END IF
      ELSE
         NX = MINMN
      END IF
!*
      DO 30 I = 1, MINMN - NX, NB
!*
!*        Reduce rows and columns i:i+nb-1 to bidiagonal form and return
!*        the matrices X and Y which are needed to update the unreduced
!*        part of the matrix
!*
         CALL DLABRD( M-I+1, N-I+1, NB, A( I, I ), LDA, D( I ), E( I ), &
                     TAUQ( I ), TAUP( I ), WORK, LDWRKX, &
                     WORK( LDWRKX*NB+1 ), LDWRKY )
!*
!*        Update the trailing submatrix A(i+nb:m,i+nb:n), using an update
!*        of the form  A := A - V*Y**T - X*U**T
!*
         CALL DGEMM( 'No transpose', 'Transpose', M-I-NB+1, N-I-NB+1, &
                    NB, -ONE, A( I+NB, I ), LDA, &
                    WORK( LDWRKX*NB+NB+1 ), LDWRKY, ONE, &
                    A( I+NB, I+NB ), LDA )
         CALL DGEMM( 'No transpose', 'No transpose', M-I-NB+1, N-I-NB+1, &
                    NB, -ONE, WORK( NB+1 ), LDWRKX, A( I, I+NB ), LDA, &
                   ONE, A( I+NB, I+NB ), LDA )
!*
!*        Copy diagonal and off-diagonal elements of B back into A
!*
         IF( M.GE.N ) THEN
            DO 10 J = I, I + NB - 1
               A( J, J ) = D( J )
               A( J, J+1 ) = E( J )
   10       CONTINUE
         ELSE
            DO 20 J = I, I + NB - 1
               A( J, J ) = D( J )
               A( J+1, J ) = E( J )
   20       CONTINUE
         END IF
   30 CONTINUE
!*
!*     Use unblocked code to reduce the remainder of the matrix
!*
      CALL DGEBD2( M-I+1, N-I+1, A( I, I ), LDA, D( I ), E( I ), &
                   TAUQ( I ), TAUP( I ), WORK, IINFO )
      WORK( 1 ) = WS
    
END SUBROUTINE
    

!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO ) !GCC$ ATTRIBUTES hot :: DGEBD2 !GCC$ ATTRIBUTES aligned(32) :: DGEBD2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I
!*     ..
!*     .. External Subroutines ..
      !EXTERNAL           DLARF, DLARFG, XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters
!*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.LT.0 ) THEN
         !CALL XERBLA( 'DGEBD2', -INFO )
         RETURN
      END IF
!*
      IF( M.GE.N ) THEN
!*
!*        Reduce to upper bidiagonal form
!*
         DO 10 I = 1, N
!*
!*           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!*
            CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, &
                         TAUQ( I ) )
            D( I ) = A( I, I )
            A( I, I ) = ONE
!*
!*           Apply H(i) to A(i:m,i+1:n) from the left
!*
            IF( I.LT.N ) &
              CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAUQ( I ), &
                          A( I, I+1 ), LDA, WORK )
            A( I, I ) = D( I )

            IF( I.LT.N ) THEN
!*
!!*              Generate elementary reflector G(i) to annihilate
!*              A(i,i+2:n)
!*
               CALL DLARFG( N-I, A( I, I+1 ), A( I, MIN( I+2, N ) ), &
                            LDA, TAUP( I ) )
               E( I ) = A( I, I+1 )
               A( I, I+1 ) = ONE
!*
!*              Apply G(i) to A(i+1:m,i+1:n) from the right
!*
               CALL DLARF( 'Right', M-I, N-I, A( I, I+1 ), LDA, &
                           TAUP( I ), A( I+1, I+1 ), LDA, WORK )
               A( I, I+1 ) = E( I )
            ELSE
               TAUP( I ) = ZERO
            END IF
   10    CONTINUE
      ELSE
!*
!*        Reduce to lower bidiagonal form
!*
         DO 20 I = 1, M
!*
!*           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
!*
            CALL DLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA, &
                         TAUP( I ) )
            D( I ) = A( I, I )
            A( I, I ) = ONE
!*
!*           Apply G(i) to A(i+1:m,i:n) from the right
!*
            IF( I.LT.M ) &
              CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, &
                          TAUP( I ), A( I+1, I ), LDA, WORK )
            A( I, I ) = D( I )

            IF( I.LT.M ) THEN
!*
!*              Generate elementary reflector H(i) to annihilate
!*              A(i+2:m,i)
!*
               CALL DLARFG( M-I, A( I+1, I ), A( MIN( I+2, M ), I ), 1, &
                           TAUQ( I ) )
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE
!*
!*              Apply H(i) to A(i+1:m,i+1:n) from the left
!*
               CALL DLARF( 'Left', M-I, N-I, A( I+1, I ), 1, TAUQ( I ), &
                          A( I+1, I+1 ), LDA, WORK )
               A( I+1, I ) = E( I )
            ELSE
               TAUQ( I ) = ZERO
            END IF
   20    CONTINUE
      END IF
     
END SUBROUTINE


#if 0
*>
*> \verbatim
*>
*> DGESVD computes the singular value decomposition (SVD) of a real
*> M-by-N matrix A, optionally computing the left and/or right singular
*> vectors. The SVD is written
*>
*>      A = U * SIGMA * transpose(V)
*>
*> where SIGMA is an M-by-N matrix which is zero except for its
*> min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
*> V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
*> are the singular values of A; they are real and non-negative, and
*> are returned in descending order.  The first min(m,n) columns of
*> U and V are the left and right singular vectors of A.
*>
*> Note that the routine returns V**T, not V.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBU
*> \verbatim
*>          JOBU is CHARACTER*1
*>          Specifies options for computing all or part of the matrix U:
*>          = 'A':  all M columns of U are returned in array U:
*>          = 'S':  the first min(m,n) columns of U (the left singular
*>                  vectors) are returned in the array U;
*>          = 'O':  the first min(m,n) columns of U (the left singular
*>                  vectors) are overwritten on the array A;
*>          = 'N':  no columns of U (no left singular vectors) are
*>                  computed.
*> \endverbatim
*>
*> \param[in] JOBVT
*> \verbatim
*>          JOBVT is CHARACTER*1
*>          Specifies options for computing all or part of the matrix
*>          V**T:
*>          = 'A':  all N rows of V**T are returned in the array VT;
*>          = 'S':  the first min(m,n) rows of V**T (the right singular
*>                  vectors) are returned in the array VT;
*>          = 'O':  the first min(m,n) rows of V**T (the right singular
*>                  vectors) are overwritten on the array A;
*>          = 'N':  no rows of V**T (no right singular vectors) are
*>                  computed.
*>
*>          JOBVT and JOBU cannot both be 'O'.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the input matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the input matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit,
*>          if JOBU = 'O',  A is overwritten with the first min(m,n)
*>                          columns of U (the left singular vectors,
*>                          stored columnwise);
*>          if JOBVT = 'O', A is overwritten with the first min(m,n)
*>                          rows of V**T (the right singular vectors,
*>                          stored rowwise);
*>          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
*>                          are destroyed.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] S
*> \verbatim
*>          S is DOUBLE PRECISION array, dimension (min(M,N))
*>          The singular values of A, sorted so that S(i) >= S(i+1).
*> \endverbatim
*>
*> \param[out] U
*> \verbatim
*>          U is DOUBLE PRECISION array, dimension (LDU,UCOL)
*>          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
*>          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
*>          if JOBU = 'S', U contains the first min(m,n) columns of U
*>          (the left singular vectors, stored columnwise);
*>          if JOBU = 'N' or 'O', U is not referenced.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of the array U.  LDU >= 1; if
*>          JOBU = 'S' or 'A', LDU >= M.
*> \endverbatim
*>
*> \param[out] VT
*> \verbatim
*>          VT is DOUBLE PRECISION array, dimension (LDVT,N)
*>          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
*>          V**T;
*>          if JOBVT = 'S', VT contains the first min(m,n) rows of
*>          V**T (the right singular vectors, stored rowwise);
*>          if JOBVT = 'N' or 'O', VT is not referenced.
*> \endverbatim
*>
*> \param[in] LDVT
*> \verbatim
*>          LDVT is INTEGER
*>          The leading dimension of the array VT.  LDVT >= 1; if
*>          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
*>          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
*>          superdiagonal elements of an upper bidiagonal matrix B
*>          whose diagonal is in S (not necessarily sorted). B
*>          satisfies A = U * B * VT, so it has the same singular values
*>          as A, and singular vectors related by U and VT.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
*>             - PATH 1  (M much larger than N, JOBU='N')
*>             - PATH 1t (N much larger than M, JOBVT='N')
*>          LWORK >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) for the other paths
*>          For good performance, LWORK should generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  if DBDSQR did not converge, INFO specifies how many
*>                superdiagonals of an intermediate bidiagonal form B
*>                did not converge to zero. See the description of WORK
*>                above for details.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date April 2012
*
*> \ingroup doubleGEsing
*
*  =====================================================================
#endif
SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, &
     VT, LDVT, WORK, LWORK, INFO ) !GCC$ ATTRIBUTES hot :: DGESVD !GCC$ ATTRIBUTES aligned(32) :: DGESVD
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LQUERY, WNTUA, WNTUAS, WNTUN, WNTUO, WNTUS, &
                        WNTVA, WNTVAS, WNTVN, WNTVO, WNTVS
      INTEGER            BDSPAC, BLK, CHUNK, I, IE, IERR, IR, ISCL, &
                        ITAU, ITAUP, ITAUQ, IU, IWORK, LDWRKR, LDWRKU, &
                        MAXWRK, MINMN, MINWRK, MNTHR, NCU, NCVT, NRU, &
                        NRVT, WRKBL
      INTEGER            LWORK_DGEQRF, LWORK_DORGQR_N, LWORK_DORGQR_M, &
                        LWORK_DGEBRD, LWORK_DORGBR_P, LWORK_DORGBR_Q, &
                        LWORK_DGELQF, LWORK_DORGLQ_N, LWORK_DORGLQ_M &
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, SMLNUM
!*     ..
!*     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 )
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DGEMM
                        
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANGE
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      MINMN = MIN( M, N )
      WNTUA = LSAME( JOBU, 'A' )
      WNTUS = LSAME( JOBU, 'S' )
      WNTUAS = WNTUA .OR. WNTUS
      WNTUO = LSAME( JOBU, 'O' )
      WNTUN = LSAME( JOBU, 'N' )
      WNTVA = LSAME( JOBVT, 'A' )
      WNTVS = LSAME( JOBVT, 'S' )
      WNTVAS = WNTVA .OR. WNTVS
      WNTVO = LSAME( JOBVT, 'O' )
      WNTVN = LSAME( JOBVT, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!*
      IF( .NOT.( WNTUA .OR. WNTUS .OR. WNTUO .OR. WNTUN ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( WNTVA .OR. WNTVS .OR. WNTVO .OR. WNTVN ) .OR. &
              ( WNTVO .AND. WNTUO ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDU.LT.1 .OR. ( WNTUAS .AND. LDU.LT.M ) ) THEN
         INFO = -9
      ELSE IF( LDVT.LT.1 .OR. ( WNTVA .AND. LDVT.LT.N ) .OR. &
              ( WNTVS .AND. LDVT.LT.MINMN ) ) THEN
         INFO = -11
      END IF
!*
!*     Compute workspace
!*      (Note: Comments in the code beginning "Workspace:" describe the
!*       minimal amount of workspace needed at that point in the code,
!*       as well as the preferred amount for good performance.
!*       NB refers to the optimal block size for the immediately
!*       following subroutine, as returned by ILAENV.)
!*
      IF( INFO.EQ.0 ) THEN
         MINWRK = 1
         MAXWRK = 1
         IF( M.GE.N .AND. MINMN.GT.0 ) THEN
!*
!*           Compute space needed for DBDSQR
!*
            MNTHR = ILAENV( 6, 'DGESVD', JOBU // JOBVT, M, N, 0, 0 )
            BDSPAC = 5*N
!*           Compute space needed for DGEQRF
            CALL DGEQRF( M, N, A, LDA, DUM(1), DUM(1), -1, IERR )
            LWORK_DGEQRF = INT( DUM(1) )
!*           Compute space needed for DORGQR
            CALL DORGQR( M, N, N, A, LDA, DUM(1), DUM(1), -1, IERR )
            LWORK_DORGQR_N = INT( DUM(1) )
            CALL DORGQR( M, M, N, A, LDA, DUM(1), DUM(1), -1, IERR )
            LWORK_DORGQR_M = INT( DUM(1) )
!*           Compute space needed for DGEBRD
            CALL DGEBRD( N, N, A, LDA, S, DUM(1), DUM(1), &
                        DUM(1), DUM(1), -1, IERR )
            LWORK_DGEBRD = INT( DUM(1) )
!*           Compute space needed for DORGBR P
            CALL DORGBR( 'P', N, N, N, A, LDA, DUM(1), &
                        DUM(1), -1, IERR )
            LWORK_DORGBR_P = INT( DUM(1) )
!*           Compute space needed for DORGBR Q
            CALL DORGBR( 'Q', N, N, N, A, LDA, DUM(1), &
                        DUM(1), -1, IERR )
            LWORK_DORGBR_Q = INT( DUM(1) )
!*
            IF( M.GE.MNTHR ) THEN
               IF( WNTUN ) THEN
!*
!*                 Path 1 (M much larger than N, JOBU='N')
!*
                  MAXWRK = N + LWORK_DGEQRF
                  MAXWRK = MAX( MAXWRK, 3*N + LWORK_DGEBRD )
                  IF( WNTVO .OR. WNTVAS ) &
                    MAXWRK = MAX( MAXWRK, 3*N + LWORK_DORGBR_P )
                  MAXWRK = MAX( MAXWRK, BDSPAC )
                  MINWRK = MAX( 4*N, BDSPAC )
               ELSE IF( WNTUO .AND. WNTVN ) THEN

!*                 Path 2 (M much larger than N, JOBU='O', JOBVT='N')
!*
                  WRKBL = N + LWORK_DGEQRF
                  WRKBL = MAX( WRKBL, N + LWORK_DORGQR_N )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = MAX( N*N + WRKBL, N*N + M*N + N )
                  MINWRK = MAX( 3*N + M, BDSPAC )
               ELSE IF( WNTUO .AND. WNTVAS ) THEN
!*
!*                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
!*                 'A')
!*
                  WRKBL = N + LWORK_DGEQRF
                  WRKBL = MAX( WRKBL, N + LWORK_DORGQR_N )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = MAX( N*N + WRKBL, N*N + M*N + N )
                  MINWRK = MAX( 3*N + M, BDSPAC )
               ELSE IF( WNTUS .AND. WNTVN ) THEN
!*
!*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
!*
                  WRKBL = N + LWORK_DGEQRF
                  WRKBL = MAX( WRKBL, N + LWORK_DORGQR_N )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = N*N + WRKBL
                  MINWRK = MAX( 3*N + M, BDSPAC )
               ELSE IF( WNTUS .AND. WNTVO ) THEN
!*
!*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
!*
                  WRKBL = N + LWORK_DGEQRF
                  WRKBL = MAX( WRKBL, N + LWORK_DORGQR_N )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = 2*N*N + WRKBL
                  MINWRK = MAX( 3*N + M, BDSPAC )
               ELSE IF( WNTUS .AND. WNTVAS ) THEN
!*
!*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
!*                 'A')
!*
                  WRKBL = N + LWORK_DGEQRF
                  WRKBL = MAX( WRKBL, N + LWORK_DORGQR_N )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = N*N + WRKBL
                  MINWRK = MAX( 3*N + M, BDSPAC )
               ELSE IF( WNTUA .AND. WNTVN ) THEN
!*
!*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
!*
                  WRKBL = N + LWORK_DGEQRF
                  WRKBL = MAX( WRKBL, N + LWORK_DORGQR_M )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = N*N + WRKBL
                  MINWRK = MAX( 3*N + M, BDSPAC )
               ELSE IF( WNTUA .AND. WNTVO ) THEN
!*
!*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
!*
                  WRKBL = N + LWORK_DGEQRF
                  WRKBL = MAX( WRKBL, N + LWORK_DORGQR_M )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = 2*N*N + WRKBL
                  MINWRK = MAX( 3*N + M, BDSPAC )
               ELSE IF( WNTUA .AND. WNTVAS ) THEN
!*
!*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
!*                 'A')
!*
                  WRKBL = N + LWORK_DGEQRF
                  WRKBL = MAX( WRKBL, N + LWORK_DORGQR_M )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = N*N + WRKBL
                  MINWRK = MAX( 3*N + M, BDSPAC )
               END IF
            ELSE
!*
!*              Path 10 (M at least N, but not much larger)
!*
               CALL DGEBRD( M, N, A, LDA, S, DUM(1), DUM(1), &
                        DUM(1), DUM(1), -1, IERR )
               LWORK_DGEBRD = INT( DUM(1) )
               MAXWRK = 3*N + LWORK_DGEBRD
               IF( WNTUS .OR. WNTUO ) THEN
                  CALL DORGBR( 'Q', M, N, N, A, LDA, DUM(1), &
                        DUM(1), -1, IERR )
                  LWORK_DORGBR_Q = INT( DUM(1) )
                  MAXWRK = MAX( MAXWRK, 3*N + LWORK_DORGBR_Q )
               END IF
               IF( WNTUA ) THEN
                  CALL DORGBR( 'Q', M, M, N, A, LDA, DUM(1), &
                        DUM(1), -1, IERR )
                  LWORK_DORGBR_Q = INT( DUM(1) )
                  MAXWRK = MAX( MAXWRK, 3*N + LWORK_DORGBR_Q )
               END IF
               IF( .NOT.WNTVN ) THEN
                 MAXWRK = MAX( MAXWRK, 3*N + LWORK_DORGBR_P )
               END IF
               MAXWRK = MAX( MAXWRK, BDSPAC )
               MINWRK = MAX( 3*N + M, BDSPAC )
            END IF
         ELSE IF( MINMN.GT.0 ) THEN
!*
!*           Compute space needed for DBDSQR
!*
            MNTHR = ILAENV( 6, 'DGESVD', JOBU // JOBVT, M, N, 0, 0 )
            BDSPAC = 5*M
!*           Compute space needed for DGELQF
            CALL DGELQF( M, N, A, LDA, DUM(1), DUM(1), -1, IERR )
            LWORK_DGELQF = INT( DUM(1) )
!*           Compute space needed for DORGLQ
            CALL DORGLQ( N, N, M, DUM(1), N, DUM(1), DUM(1), -1, IERR )
            LWORK_DORGLQ_N = INT( DUM(1) )
            CALL DORGLQ( M, N, M, A, LDA, DUM(1), DUM(1), -1, IERR )
            LWORK_DORGLQ_M = INT( DUM(1) )
!*           Compute space needed for DGEBRD
            CALL DGEBRD( M, M, A, LDA, S, DUM(1), DUM(1), &
                        DUM(1), DUM(1), -1, IERR )
            LWORK_DGEBRD = INT( DUM(1) )
!*            Compute space needed for DORGBR P
            CALL DORGBR( 'P', M, M, M, A, N, DUM(1), &
                        DUM(1), -1, IERR )
            LWORK_DORGBR_P = INT( DUM(1) )
!*           Compute space needed for DORGBR Q
            CALL DORGBR( 'Q', M, M, M, A, N, DUM(1), &
                        DUM(1), -1, IERR )
            LWORK_DORGBR_Q = INT( DUM(1) )
            IF( N.GE.MNTHR ) THEN
               IF( WNTVN ) THEN
!*
!*                 Path 1t(N much larger than M, JOBVT='N')
!*
                  MAXWRK = M + LWORK_DGELQF
                  MAXWRK = MAX( MAXWRK, 3*M + LWORK_DGEBRD )
                  IF( WNTUO .OR. WNTUAS ) &
                    MAXWRK = MAX( MAXWRK, 3*M + LWORK_DORGBR_Q )
                  MAXWRK = MAX( MAXWRK, BDSPAC )
                  MINWRK = MAX( 4*M, BDSPAC )
               ELSE IF( WNTVO .AND. WNTUN ) THEN
!*
!*                 Path 2t(N much larger than M, JOBU='N', JOBVT='O')
!*
                  WRKBL = M + LWORK_DGELQF
                  WRKBL = MAX( WRKBL, M + LWORK_DORGLQ_M )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = MAX( M*M + WRKBL, M*M + M*N + M )
                  MINWRK = MAX( 3*M + N, BDSPAC )
               ELSE IF( WNTVO .AND. WNTUAS ) THEN
!*
!*                 Path 3t(N much larger than M, JOBU='S' or 'A',
!*                 JOBVT='O')
!*
                  WRKBL = M + LWORK_DGELQF
                  WRKBL = MAX( WRKBL, M + LWORK_DORGLQ_M )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = MAX( M*M + WRKBL, M*M + M*N + M )
                  MINWRK = MAX( 3*M + N, BDSPAC )
               ELSE IF( WNTVS .AND. WNTUN ) THEN
!*
!*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S')
!*
                  WRKBL = M + LWORK_DGELQF
                  WRKBL = MAX( WRKBL, M + LWORK_DORGLQ_M )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = M*M + WRKBL
                  MINWRK = MAX( 3*M + N, BDSPAC )
               ELSE IF( WNTVS .AND. WNTUO ) THEN
*
*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
*
                  WRKBL = M + LWORK_DGELQF
                  WRKBL = MAX( WRKBL, M + LWORK_DORGLQ_M )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = 2*M*M + WRKBL
                  MINWRK = MAX( 3*M + N, BDSPAC )
               ELSE IF( WNTVS .AND. WNTUAS ) THEN
!*
!*                 Path 6t(N much larger than M, JOBU='S' or 'A',
!*                 JOBVT='S')
!*
                  WRKBL = M + LWORK_DGELQF
                  WRKBL = MAX( WRKBL, M + LWORK_DORGLQ_M )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = M*M + WRKBL
                  MINWRK = MAX( 3*M + N, BDSPAC )
               ELSE IF( WNTVA .AND. WNTUN ) THEN
!*
!*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A')
!*
                  WRKBL = M + LWORK_DGELQF
                  WRKBL = MAX( WRKBL, M + LWORK_DORGLQ_N )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = M*M + WRKBL
                  MINWRK = MAX( 3*M + N, BDSPAC )
               ELSE IF( WNTVA .AND. WNTUO ) THEN
!*
!*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
!*
                  WRKBL = M + LWORK_DGELQF
                  WRKBL = MAX( WRKBL, M + LWORK_DORGLQ_N )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = 2*M*M + WRKBL
                  MINWRK = MAX( 3*M + N, BDSPAC )
               ELSE IF( WNTVA .AND. WNTUAS ) THEN
!*
!!*                 Path 9t(N much larger than M, JOBU='S' or 'A',
!*                 JOBVT='A')
!*
                  WRKBL = M + LWORK_DGELQF
                  WRKBL = MAX( WRKBL, M + LWORK_DORGLQ_N )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_P )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = M*M + WRKBL
                  MINWRK = MAX( 3*M + N, BDSPAC )
               END IF
            ELSE
!*
!!*              Path 10t(N greater than M, but not much larger)
!*
               CALL DGEBRD( M, N, A, LDA, S, DUM(1), DUM(1), &
                        DUM(1), DUM(1), -1, IERR )
               LWORK_DGEBRD = INT( DUM(1) )
               MAXWRK = 3*M + LWORK_DGEBRD
               IF( WNTVS .OR. WNTVO ) THEN
                Compute space needed for DORGBR P
                 CALL DORGBR( 'P', M, N, M, A, N, DUM(1), &
                        DUM(1), -1, IERR )
                 LWORK_DORGBR_P = INT( DUM(1) )
                 MAXWRK = MAX( MAXWRK, 3*M + LWORK_DORGBR_P )
               END IF
               IF( WNTVA ) THEN
                 CALL DORGBR( 'P', N, N, M, A, N, DUM(1), &
                        DUM(1), -1, IERR )
                 LWORK_DORGBR_P = INT( DUM(1) )
                 MAXWRK = MAX( MAXWRK, 3*M + LWORK_DORGBR_P )
               END IF
               IF( .NOT.WNTUN ) THEN
                  MAXWRK = MAX( MAXWRK, 3*M + LWORK_DORGBR_Q )
               END IF
               MAXWRK = MAX( MAXWRK, BDSPAC )
               MINWRK = MAX( 3*M + N, BDSPAC )
            END IF
         END IF
         MAXWRK = MAX( MAXWRK, MINWRK )
         WORK( 1 ) = MAXWRK

         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -13
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGESVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

!*     Quick return if possible
!*
    !  IF( M.EQ.0 .OR. N.EQ.0 ) THEN
    !     RETURN
    !  END IF
!*
!*     Get machine constants
!*
      EPS = DLAMCH( 'P' )
      SMLNUM = SQRT( DLAMCH( 'S' ) ) / EPS
      BIGNUM = ONE / SMLNUM
!*
!*     Scale A if max element outside range [SMLNUM,BIGNUM]
!*
      ANRM = DLANGE( 'M', M, N, A, LDA, DUM )
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ISCL = 1
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR )
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ISCL = 1
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR )
      END IF
!*
      IF( M.GE.N ) THEN
!*
!*        A has at least as many rows as columns. If A has sufficiently
!*        more rows than columns, first reduce using the QR
!*        decomposition (if sufficient workspace available)
!*
         IF( M.GE.MNTHR ) THEN
!*
            IF( WNTUN ) THEN
!*
!*              Path 1 (M much larger than N, JOBU='N')
!*              No left singular vectors to be computed
!*
               ITAU = 1
               IWORK = ITAU + N
!*
!*              Compute A=Q*R
!*              (Workspace: need 2*N, prefer N + N*NB)
!*
               CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                           LWORK-IWORK+1, IERR )
!*
!*              Zero out below R
!*
               IF( N .GT. 1 ) THEN
                  CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), &
                               LDA )
               END IF
               IE = 1
               ITAUQ = IE + N
               ITAUP = ITAUQ + N
               IWORK = ITAUP + N
!*
!*              Bidiagonalize R in A
!*              (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!*
               CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                           WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, &
                           IERR )
               NCVT = 0
               IF( WNTVO .OR. WNTVAS ) THEN
!*
!*                 If right singular vectors desired, generate P'.
!*                 (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!*
                  CALL DORGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  NCVT = N
               END IF
               IWORK = IE + N
!*
!*              Perform bidiagonal QR iteration, computing right
!*              singular vectors of A in A if desired
!*              (Workspace: need BDSPAC)
!*
               CALL DBDSQR( 'U', N, NCVT, 0, 0, S, WORK( IE ), A, LDA, &
                           DUM, 1, DUM, 1, WORK( IWORK ), INFO )
!*
!*              If right singular vectors desired in VT, copy them there
!*
               IF( WNTVAS ) &
                  CALL DLACPY( 'F', N, N, A, LDA, VT, LDVT )

            ELSE IF( WNTUO .AND. WNTVN ) THEN
!*
!*              Path 2 (M much larger than N, JOBU='O', JOBVT='N')
!*              N left singular vectors to be overwritten on A and
!*              no right singular vectors to be computed
!*
               IF( LWORK.GE.N*N+MAX( 4*N, BDSPAC ) ) THEN
!*
!*                 Sufficient workspace for a fast algorithm
!*
                  IR = 1
                  IF( LWORK.GE.MAX( WRKBL, LDA*N + N ) + LDA*N ) THEN
!*
!*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N
!*
                     LDWRKU = LDA
                     LDWRKR = LDA
                  ELSE IF( LWORK.GE.MAX( WRKBL, LDA*N + N ) + N*N ) THEN
!*
!*                    WORK(IU) is LDA by N, WORK(IR) is N by N
!*
                     LDWRKU = LDA
                     LDWRKR = N
                  ELSE
!*
!*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N
!*
                     LDWRKU = ( LWORK-N*N-N ) / N
                     LDWRKR = N
                  END IF
                  ITAU = IR + LDWRKR*N
                  IWORK = ITAU + N
!*
!*                 Compute A=Q*R
!*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!*
                  CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                 Copy R to WORK(IR) and zero out below it
!*
                  CALL DLACPY( 'U', N, N, A, LDA, WORK( IR ), LDWRKR )
                  CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, WORK( IR+1 ), &
                               LDWRKR )
!*
!*                 Generate Q in A
!*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!*
                  CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ), &
                               WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = ITAU
                  ITAUQ = IE + N
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N
!*
!*                 Bidiagonalize R in WORK(IR)
!*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
!*
                  CALL DGEBRD( N, N, WORK( IR ), LDWRKR, S, WORK( IE ), &
                              WORK( ITAUQ ), WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                 Generate left vectors bidiagonalizing R
!*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
!*
                  CALL DORGBR( 'Q', N, N, N, WORK( IR ), LDWRKR, &
                              WORK( ITAUQ ), WORK( IWORK ), &
                              LWORK-IWORK+1, IERR )
                  IWORK = IE + N
!*
!*                 Perform bidiagonal QR iteration, computing left
!*                 singular vectors of R in WORK(IR)
!*                 (Workspace: need N*N + BDSPAC)
!*
                  CALL DBDSQR( 'U', N, 0, N, 0, S, WORK( IE ), DUM, 1, &
                              WORK( IR ), LDWRKR, DUM, 1, &
                              WORK( IWORK ), INFO )
                  IU = IE + N
!*
!*                 Multiply Q in A by left singular vectors of R in
!*                 WORK(IR), storing result in WORK(IU) and copying to A
!!*                 (Workspace: need N*N + 2*N, prefer N*N + M*N + N)
!*
                  DO 10 I = 1, M, LDWRKU
                     CHUNK = MIN( M-I+1, LDWRKU )
                     CALL DGEMM( 'N', 'N', CHUNK, N, N, ONE, A( I, 1 ), &
                                LDA, WORK( IR ), LDWRKR, ZERO, &
                                WORK( IU ), LDWRKU )
                     CALL DLACPY( 'F', CHUNK, N, WORK( IU ), LDWRKU, &
                                 A( I, 1 ), LDA )
   10             CONTINUE

               ELSE
!*
!*                 Insufficient workspace for a fast algorithm
!*
                  IE = 1
                  ITAUQ = IE + N
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N
!*
!*                 Bidiagonalize A
!*                 (Workspace: need 3*N + M, prefer 3*N + (M + N)*NB)
!*
                  CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), &
                              WORK( ITAUQ ), WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                 Generate left vectors bidiagonalizing A
!*                 (Workspace: need 4*N, prefer 3*N + N*NB)
!*
                  CALL DORGBR( 'Q', M, N, N, A, LDA, WORK( ITAUQ ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IWORK = IE + N
!*
!*                 Perform bidiagonal QR iteration, computing left
!*                 singular vectors of A in A
!*                 (Workspace: need BDSPAC)
!*
                  CALL DBDSQR( 'U', N, 0, M, 0, S, WORK( IE ), DUM, 1, &
                               A, LDA, DUM, 1, WORK( IWORK ), INFO )
!*
               END IF
!*
            ELSE IF( WNTUO .AND. WNTVAS ) THEN
!*
!*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
!*              N left singular vectors to be overwritten on A and
!*              N right singular vectors to be computed in VT
!*
               IF( LWORK.GE.N*N+MAX( 4*N, BDSPAC ) ) THEN
!*
!*                 Sufficient workspace for a fast algorithm
!*
                  IR = 1
                  IF( LWORK.GE.MAX( WRKBL, LDA*N + N ) + LDA*N ) THEN
!*
!*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N
!*
                     LDWRKU = LDA
                     LDWRKR = LDA
                  ELSE IF( LWORK.GE.MAX( WRKBL, LDA*N + N ) + N*N ) THEN
!*
!*                    WORK(IU) is LDA by N and WORK(IR) is N by N
!*
                     LDWRKU = LDA
                     LDWRKR = N
                  ELSE
!*
!*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N
!*
                     LDWRKU = ( LWORK-N*N-N ) / N
                     LDWRKR = N
                  END IF
                  ITAU = IR + LDWRKR*N
                  IWORK = ITAU + N
!*
!*                 Compute A=Q*R
!*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!*
                  CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                 Copy R to VT, zeroing out below it
!*
                  CALL DLACPY( 'U', N, N, A, LDA, VT, LDVT )
                  IF( N.GT.1 ) &
                    CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                 VT( 2, 1 ), LDVT )
!*
!*                 Generate Q in A
!*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!*
                  CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ),  &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = ITAU
                  ITAUQ = IE + N
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N
!*
!*                 Bidiagonalize R in VT, copying result to WORK(IR)
!*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
!*
                  CALL DGEBRD( N, N, VT, LDVT, S, WORK( IE ), &
                              WORK( ITAUQ ), WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  CALL DLACPY( 'L', N, N, VT, LDVT, WORK( IR ), LDWRKR )
!!*
!*                 Generate left vectors bidiagonalizing R in WORK(IR)
!*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
!*
                  CALL DORGBR( 'Q', N, N, N, WORK( IR ), LDWRKR, &
                              WORK( ITAUQ ), WORK( IWORK ), &
                              LWORK-IWORK+1, IERR )
!*
!*                 Generate right vectors bidiagonalizing R in VT
!*                 (Workspace: need N*N + 4*N-1, prefer N*N + 3*N + (N-1)*NB)
!*
                  CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IWORK = IE + N
!*
!*                 Perform bidiagonal QR iteration, computing left
!*                 singular vectors of R in WORK(IR) and computing right
!*                 singular vectors of R in VT
!*                 (Workspace: need N*N + BDSPAC)
!*
                  CALL DBDSQR( 'U', N, N, N, 0, S, WORK( IE ), VT, LDVT, &
                              WORK( IR ), LDWRKR, DUM, 1, &
                              WORK( IWORK ), INFO )
                  IU = IE + N
!*
!*                 Multiply Q in A by left singular vectors of R in
!*                 WORK(IR), storing result in WORK(IU) and copying to A
!!*                 (Workspace: need N*N + 2*N, prefer N*N + M*N + N)
!*
                  DO 20 I = 1, M, LDWRKU
                     CHUNK = MIN( M-I+1, LDWRKU )
                     CALL DGEMM( 'N', 'N', CHUNK, N, N, ONE, A( I, 1 ), &
                                LDA, WORK( IR ), LDWRKR, ZERO, &
                                WORK( IU ), LDWRKU )
                     CALL DLACPY( 'F', CHUNK, N, WORK( IU ), LDWRKU, &
                                 A( I, 1 ), LDA )
   20             CONTINUE
!*
               ELSE
!*
!*                 Insufficient workspace for a fast algorithm
!*
                  ITAU = 1
                  IWORK = ITAU + N
!*
!*                 Compute A=Q*R
!*                 (Workspace: need 2*N, prefer N + N*NB)
!*
                  CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                 Copy R to VT, zeroing out below it
!*
                  CALL DLACPY( 'U', N, N, A, LDA, VT, LDVT )
                  IF( N.GT.1 ) &
                    CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                 VT( 2, 1 ), LDVT )
!*
!*                 Generate Q in A
!*                 (Workspace: need 2*N, prefer N + N*NB)
!*
                  CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = ITAU
                  ITAUQ = IE + N
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N

                  
!*                 Bidiagonalize R in VT
!*                 (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!*
                  CALL DGEBRD( N, N, VT, LDVT, S, WORK( IE ), &
                              WORK( ITAUQ ), WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                 Multiply Q in A by left vectors bidiagonalizing R
!*                 (Workspace: need 3*N + M, prefer 3*N + M*NB)
!*
                  CALL DORMBR( 'Q', 'R', 'N', M, N, N, VT, LDVT, &
                              WORK( ITAUQ ), A, LDA, WORK( IWORK ), &
                              LWORK-IWORK+1, IERR )
!*
!*                 Generate right vectors bidiagonalizing R in VT
!*                 (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!*
                  CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IWORK = IE + N
!*
!*                 Perform bidiagonal QR iteration, computing left
!*                 singular vectors of A in A and computing right
!*                 singular vectors of A in VT
!*                 (Workspace: need BDSPAC)
!*
                  CALL DBDSQR( 'U', N, N, M, 0, S, WORK( IE ), VT, LDVT, &
                               A, LDA, DUM, 1, WORK( IWORK ), INFO )

               END IF

            ELSE IF( WNTUS ) THEN

               IF( WNTVN ) THEN
!*
!*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
!*                 N left singular vectors to be computed in U and
!*                 no right singular vectors to be computed
!*
                  IF( LWORK.GE.N*N+MAX( 4*N, BDSPAC ) ) THEN
!*
!*                    Sufficient workspace for a fast algorithm
!*
                     IR = 1
                     IF( LWORK.GE.WRKBL+LDA*N ) THEN
!*
!*                       WORK(IR) is LDA by N
!*
                        LDWRKR = LDA
                     ELSE
!*
!*                       WORK(IR) is N by N
!*
                        LDWRKR = N
                     END IF
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N
!*
!*                    Compute A=Q*R
!*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!*
                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                  WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Copy R to WORK(IR), zeroing out below it
!*
                     CALL DLACPY( 'U', N, N, A, LDA, WORK( IR ), &
                                 LDWRKR )
                     CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                 WORK( IR+1 ), LDWRKR )
!*
!*                    Generate Q in A
!*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!*
                     CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
!*
!*                    Bidiagonalize R in WORK(IR)
!*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
!*
                     CALL DGEBRD( N, N, WORK( IR ), LDWRKR, S, &
                                WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
!*
!*                    Generate left vectors bidiagonalizing R in WORK(IR)
!*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
!*
                     CALL DORGBR( 'Q', N, N, N, WORK( IR ), LDWRKR, &
                                 WORK( ITAUQ ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     IWORK = IE + N
!*
!*                    Perform bidiagonal QR iteration, computing left
!*                    singular vectors of R in WORK(IR)
!*                    (Workspace: need N*N + BDSPAC)
!*
                     CALL DBDSQR( 'U', N, 0, N, 0, S, WORK( IE ), DUM, &
                                 1, WORK( IR ), LDWRKR, DUM, 1, &
                                 WORK( IWORK ), INFO )
!*
!*                    Multiply Q in A by left singular vectors of R in
!*                    WORK(IR), storing result in U
!*                    (Workspace: need N*N)
!*
                     CALL DGEMM( 'N', 'N', M, N, N, ONE, A, LDA, &
                                WORK( IR ), LDWRKR, ZERO, U, LDU )
!*
                  ELSE
!*
!*                    Insufficient workspace for a fast algorithm
!*
                     ITAU = 1
                     IWORK = ITAU + N
!*
!*                    Compute A=Q*R, copying result to U
!*                    (Workspace: need 2*N, prefer N + N*NB)
!*
                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, N, A, LDA, U, LDU )
!*
!*                    Generate Q in U
!*                    (Workspace: need 2*N, prefer N + N*NB)
!*
                     CALL DORGQR( M, N, N, U, LDU, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
!*
!*                    Zero out below R in A
!*!
                     IF( N .GT. 1 ) THEN
                        CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                     A( 2, 1 ), LDA )
                     END IF
!*
!*                    Bidiagonalize R in A
!*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!*
                     CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), &
                                 WORK( ITAUQ ), WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Multiply Q in U by left vectors bidiagonalizing R
!*                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
!*
                     CALL DORMBR( 'Q', 'R', 'N', M, N, N, A, LDA, &
                                 WORK( ITAUQ ), U, LDU, WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     IWORK = IE + N
!*
!*                    Perform bidiagonal QR iteration, computing left
!*                    singular vectors of A in U
!*                    (Workspace: need BDSPAC)
!*
                     CALL DBDSQR( 'U', N, 0, M, 0, S, WORK( IE ), DUM, &
                                 1, U, LDU, DUM, 1, WORK( IWORK ), &
                                 INFO )

                  END IF

               ELSE IF( WNTVO ) THEN
!*
!*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
!*                 N left singular vectors to be computed in U and
!*                 N right singular vectors to be overwritten on A
!*
                  IF( LWORK.GE.2*N*N+MAX( 4*N, BDSPAC ) ) THEN
!*
!*                    Sufficient workspace for a fast algorithm
!*
                     IU = 1
                     IF( LWORK.GE.WRKBL+2*LDA*N ) THEN
!*
!*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N
!*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = LDA
                     ELSE IF( LWORK.GE.WRKBL+( LDA + N )*N ) THEN
!*
!*                       WORK(IU) is LDA by N and WORK(IR) is N by N
!*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     ELSE
!*
!*                       WORK(IU) is N by N and WORK(IR) is N by N
!*
                        LDWRKU = N
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     END IF
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N
!*
!*                    Compute A=Q*R
!*                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB)
!*
                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Copy R to WORK(IU), zeroing out below it
!*
                     CALL DLACPY( 'U', N, N, A, LDA, WORK( IU ), &
                                 LDWRKU )
                     CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                 WORK( IU+1 ), LDWRKU )
!*
!*                    Generate Q in A
!*                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB)
!*
                     CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
!*
!*                    Bidiagonalize R in WORK(IU), copying result to
!*                    WORK(IR)
!*                    (Workspace: need 2*N*N + 4*N,
!*                                prefer 2*N*N+3*N+2*N*NB)
!*
                     CALL DGEBRD( N, N, WORK( IU ), LDWRKU, S, &
                                 WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', N, N, WORK( IU ), LDWRKU, &
                                 WORK( IR ), LDWRKR )
!*
!*                    Generate left bidiagonalizing vectors in WORK(IU)
!*                    (Workspace: need 2*N*N + 4*N, prefer 2*N*N + 3*N + N*NB)
!*
                     CALL DORGBR( 'Q', N, N, N, WORK( IU ), LDWRKU, &
                                 WORK( ITAUQ ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
!*
!*                    Generate right bidiagonalizing vectors in WORK(IR)
!*                    (Workspace: need 2*N*N + 4*N-1,
!*                                prefer 2*N*N+3*N+(N-1)*NB)
!*
                     CALL DORGBR( 'P', N, N, N, WORK( IR ), LDWRKR, &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     IWORK = IE + N
!*
!*                    Perform bidiagonal QR iteration, computing left
!*                    singular vectors of R in WORK(IU) and computing
!*                    right singular vectors of R in WORK(IR)
!*                    (Workspace: need 2*N*N + BDSPAC)
!*
                     CALL DBDSQR( 'U', N, N, N, 0, S, WORK( IE ), &
                                 WORK( IR ), LDWRKR, WORK( IU ), &
                                 LDWRKU, DUM, 1, WORK( IWORK ), INFO )
!*
!*                    Multiply Q in A by left singular vectors of R in
!*                    WORK(IU), storing result in U
!*                    (Workspace: need N*N)
!*
                     CALL DGEMM( 'N', 'N', M, N, N, ONE, A, LDA, &
                                WORK( IU ), LDWRKU, ZERO, U, LDU )
!*
!*                    Copy right singular vectors of R to A
!*                    (Workspace: need N*N)
!*
                     CALL DLACPY( 'F', N, N, WORK( IR ), LDWRKR, A, &
                                 LDA )
!*
                  ELSE
!*
!*                    Insufficient workspace for a fast algorithm
!*
                     ITAU = 1
                     IWORK = ITAU + N
!*
!*                    Compute A=Q*R, copying result to U
!*                    (Workspace: need 2*N, prefer N + N*NB)
!*
                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, N, A, LDA, U, LDU )
!*
!*                    Generate Q in U
!*                    (Workspace: need 2*N, prefer N + N*NB)
!*
                     CALL DORGQR( M, N, N, U, LDU, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
!*
!*                    Zero out below R in A
!*
                     IF( N .GT. 1 ) THEN
                        CALL DLASET( 'L', N-1, N-1, ZERO, ZERO,
                                    A( 2, 1 ), LDA )
                     END IF
!*
!*                    Bidiagonalize R in A
!*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!*
                     CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), &
                                 WORK( ITAUQ ), WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Multiply Q in U by left vectors bidiagonalizing R
!*                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
!*
                     CALL DORMBR( 'Q', 'R', 'N', M, N, N, A, LDA, &
                                 WORK( ITAUQ ), U, LDU, WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
!*
!*                    Generate right vectors bidiagonalizing R in A
!*                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!*
                     CALL DORGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ), &
                                  WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + N
!*
!*                    Perform bidiagonal QR iteration, computing left
!*                    singular vectors of A in U and computing right
!*                    singular vectors of A in A
!*                    (Workspace: need BDSPAC)
!*
                     CALL DBDSQR( 'U', N, N, M, 0, S, WORK( IE ), A, &
                                 LDA, U, LDU, DUM, 1, WORK( IWORK ), &
                                 INFO )
!*
                  END IF
!*
               ELSE IF( WNTVAS ) THEN
!*
!*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S'
!*                         or 'A')
!*                 N left singular vectors to be computed in U and
!*                 N right singular vectors to be computed in VT
!*
                  IF( LWORK.GE.N*N+MAX( 4*N, BDSPAC ) ) THEN
!*
!*                    Sufficient workspace for a fast algorithm
!*
                     IU = 1
                     IF( LWORK.GE.WRKBL+LDA*N ) THEN
!*
!*                       WORK(IU) is LDA by N
!*
                        LDWRKU = LDA
                     ELSE
!*
!*                       WORK(IU) is N by N
!*
                        LDWRKU = N
                     END IF
                     ITAU = IU + LDWRKU*N
                     IWORK = ITAU + N
!*
!*                    Compute A=Q*R
!*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!*
                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Copy R to WORK(IU), zeroing out below it
!*
                     CALL DLACPY( 'U', N, N, A, LDA, WORK( IU ), &
                                 LDWRKU )
                     CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                 WORK( IU+1 ), LDWRKU )
!*
!*                    Generate Q in A
!*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!*
                     CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
!*
!*                    Bidiagonalize R in WORK(IU), copying result to VT
!*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
!*
                     CALL DGEBRD( N, N, WORK( IU ), LDWRKU, S, &
                                 WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', N, N, WORK( IU ), LDWRKU, VT, &
                                LDVT )
!*
!*                    Generate left bidiagonalizing vectors in WORK(IU)
!*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
!*
                     CALL DORGBR( 'Q', N, N, N, WORK( IU ), LDWRKU, &
                                 WORK( ITAUQ ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
!*
!*                    Generate right bidiagonalizing vectors in VT
!*                    (Workspace: need N*N + 4*N-1,
!*                                prefer N*N+3*N+(N-1)*NB)
!*
                     CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + N
!*
!*                    Perform bidiagonal QR iteration, computing left
!*                    singular vectors of R in WORK(IU) and computing
!*                    right singular vectors of R in VT
!*                    (Workspace: need N*N + BDSPAC)
!*
                     CALL DBDSQR( 'U', N, N, N, 0, S, WORK( IE ), VT, &
                                 LDVT, WORK( IU ), LDWRKU, DUM, 1, &
                                 WORK( IWORK ), INFO )
!*
!*                    Multiply Q in A by left singular vectors of R in
!*                    WORK(IU), storing result in U
!*                    (Workspace: need N*N)
!*
                     CALL DGEMM( 'N', 'N', M, N, N, ONE, A, LDA, &
                                WORK( IU ), LDWRKU, ZERO, U, LDU )
!*
                  ELSE
!*
!*                    Insufficient workspace for a fast algorithm
!*
                     ITAU = 1
                     IWORK = ITAU + N
!*
!*                    Compute A=Q*R, copying result to U
!*                    (Workspace: need 2*N, prefer N + N*NB)
!*
                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, N, A, LDA, U, LDU )
!*
!*                    Generate Q in U
!*                    (Workspace: need 2*N, prefer N + N*NB)
!*
                     CALL DORGQR( M, N, N, U, LDU, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Copy R to VT, zeroing out below it
!*
                     CALL DLACPY( 'U', N, N, A, LDA, VT, LDVT )
                     IF( N.GT.1 ) &
                       CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                    VT( 2, 1 ), LDVT )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
!*
!*                    Bidiagonalize R in VT
!*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
!*
                     CALL DGEBRD( N, N, VT, LDVT, S, WORK( IE ), &
                                 WORK( ITAUQ ), WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Multiply Q in U by left bidiagonalizing vectors
!*                    in VT
!*                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
!*
                     CALL DORMBR( 'Q', 'R', 'N', M, N, N, VT, LDVT, &
                                 WORK( ITAUQ ), U, LDU, WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
!*
!*                    Generate right bidiagonalizing vectors in VT
!*                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
!*
                     CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + N
!*
!!*                    Perform bidiagonal QR iteration, computing left
!*!                    singular vectors of A in U and computing right
!*                    singular vectors of A in VT
!*                    (Workspace: need BDSPAC)
!*
                     CALL DBDSQR( 'U', N, N, M, 0, S, WORK( IE ), VT, &
                                 LDVT, U, LDU, DUM, 1, WORK( IWORK ), &
                                 INFO )
!*
                  END IF
!*
               END IF
!*
            ELSE IF( WNTUA ) THEN
!*
               IF( WNTVN ) THEN
!*
!*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
!*                 M left singular vectors to be computed in U and
!*                 no right singular vectors to be computed
!*
                  IF( LWORK.GE.N*N+MAX( N+M, 4*N, BDSPAC ) ) THEN
!*
!*                    Sufficient workspace for a fast algorithm
!*
                     IR = 1
                     IF( LWORK.GE.WRKBL+LDA*N ) THEN
!*
!*                       WORK(IR) is LDA by N
!*
                        LDWRKR = LDA
                     ELSE
!*
!*                       WORK(IR) is N by N
!*
                        LDWRKR = N
                     END IF
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N
!*
!*                    Compute A=Q*R, copying result to U
!*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
!*
                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, N, A, LDA, U, LDU )
!*
!*                    Copy R to WORK(IR), zeroing out below it
!*
                     CALL DLACPY( 'U', N, N, A, LDA, WORK( IR ), &
                                 LDWRKR )
                     CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                 WORK( IR+1 ), LDWRKR )
!*
!*!                    Generate Q in U
!*                    (Workspace: need N*N + N + M, prefer N*N + N + M*NB)
!*
                     CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     CALL DGEBRD( N, N, WORK( IR ), LDWRKR, S, &
                                 WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )

                     CALL DORGBR( 'Q', N, N, N, WORK( IR ), LDWRKR, &
                                 WORK( ITAUQ ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     IWORK = IE + N

                     CALL DBDSQR( 'U', N, 0, N, 0, S, WORK( IE ), DUM, &
                                 1, WORK( IR ), LDWRKR, DUM, 1, &
                                 WORK( IWORK ), INFO )

                     CALL DGEMM( 'N', 'N', M, N, N, ONE, U, LDU, &
                                WORK( IR ), LDWRKR, ZERO, A, LDA )

                     CALL DLACPY( 'F', M, N, A, LDA, U, LDU )

                  ELSE

                     ITAU = 1
                     IWORK = ITAU + N

                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                     CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     IF( N .GT. 1 ) THEN
                        CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                     A( 2, 1 ), LDA )
                     END IF

                     CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), &
                                 WORK( ITAUQ ), WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DORMBR( 'Q', 'R', 'N', M, N, N, A, LDA, &
                                 WORK( ITAUQ ), U, LDU, WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     IWORK = IE + N

                     CALL DBDSQR( 'U', N, 0, M, 0, S, WORK( IE ), DUM, &
                                 1, U, LDU, DUM, 1, WORK( IWORK ), &
                                 INFO )

                  END IF

               ELSE IF( WNTVO ) THEN

                  IF( LWORK.GE.2*N*N+MAX( N+M, 4*N, BDSPAC ) ) THEN

                     IU = 1
                     IF( LWORK.GE.WRKBL+2*LDA*N ) THEN

                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = LDA
                     ELSE IF( LWORK.GE.WRKBL+( LDA + N )*N ) THEN

                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     ELSE

                        LDWRKU = N
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     END IF
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N

                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                     CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DLACPY( 'U', N, N, A, LDA, WORK( IU ), &
                                 LDWRKU )
                     CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                 WORK( IU+1 ), LDWRKU )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     CALL DGEBRD( N, N, WORK( IU ), LDWRKU, S, &
                                 WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', N, N, WORK( IU ), LDWRKU, &
                                 WORK( IR ), LDWRKR )

                     CALL DORGBR( 'Q', N, N, N, WORK( IU ), LDWRKU, &
                                 WORK( ITAUQ ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )

                     CALL DORGBR( 'P', N, N, N, WORK( IR ), LDWRKR, &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     IWORK = IE + N

                     CALL DBDSQR( 'U', N, N, N, 0, S, WORK( IE ), &
                                 WORK( IR ), LDWRKR, WORK( IU ), &
                                 LDWRKU, DUM, 1, WORK( IWORK ), INFO )

                     CALL DGEMM( 'N', 'N', M, N, N, ONE, U, LDU, &
                                WORK( IU ), LDWRKU, ZERO, A, LDA )

                     CALL DLACPY( 'F', M, N, A, LDA, U, LDU )

                     CALL DLACPY( 'F', N, N, WORK( IR ), LDWRKR, A, &
                                 LDA )

                  ELSE

                     ITAU = 1
                     IWORK = ITAU + N

                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                     CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     IF( N .GT. 1 ) THEN
                        CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                     A( 2, 1 ), LDA )
                     END IF

                     CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), &
                                 WORK( ITAUQ ), WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DORMBR( 'Q', 'R', 'N', M, N, N, A, LDA, &
                                 WORK( ITAUQ ), U, LDU, WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )

                     CALL DORGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + N

                     CALL DBDSQR( 'U', N, N, M, 0, S, WORK( IE ), A, &
                                 LDA, U, LDU, DUM, 1, WORK( IWORK ), &
                                 INFO )

                  END IF

               ELSE IF( WNTVAS ) THEN

                  IF( LWORK.GE.N*N+MAX( N+M, 4*N, BDSPAC ) ) THEN

                     IU = 1
                     IF( LWORK.GE.WRKBL+LDA*N ) THEN

                        LDWRKU = LDA
                     ELSE

                        LDWRKU = N
                     END IF
                     ITAU = IU + LDWRKU*N
                     IWORK = ITAU + N

                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                     CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DLACPY( 'U', N, N, A, LDA, WORK( IU ), &
                                 LDWRKU )
                     CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                 WORK( IU+1 ), LDWRKU )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     CALL DGEBRD( N, N, WORK( IU ), LDWRKU, S, &
                                WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', N, N, WORK( IU ), LDWRKU, VT, &
                                LDVT )

                     CALL DORGBR( 'Q', N, N, N, WORK( IU ), LDWRKU, &
                                 WORK( ITAUQ ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )

                     CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + N

                     CALL DBDSQR( 'U', N, N, N, 0, S, WORK( IE ), VT, &
                                 LDVT, WORK( IU ), LDWRKU, DUM, 1, &
                                 WORK( IWORK ), INFO )

                     CALL DGEMM( 'N', 'N', M, N, N, ONE, U, LDU, &
                                WORK( IU ), LDWRKU, ZERO, A, LDA )

                     CALL DLACPY( 'F', M, N, A, LDA, U, LDU )

                  ELSE

                     ITAU = 1
                     IWORK = ITAU + N

                     CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                     CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DLACPY( 'U', N, N, A, LDA, VT, LDVT )
                     IF( N.GT.1 ) &
                       CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, &
                                    VT( 2, 1 ), LDVT )
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     CALL DGEBRD( N, N, VT, LDVT, S, WORK( IE ), &
                                 WORK( ITAUQ ), WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DORMBR( 'Q', 'R', 'N', M, N, N, VT, LDVT, &
                                 WORK( ITAUQ ), U, LDU, WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )

                     CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + N

                     CALL DBDSQR( 'U', N, N, M, 0, S, WORK( IE ), VT, &
                                 LDVT, U, LDU, DUM, 1, WORK( IWORK ), &
                                 INFO )

                  END IF

               END IF

            END IF

         ELSE

            IE = 1
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N

            CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                        WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, &
                        IERR )
            IF( WNTUAS ) THEN

               CALL DLACPY( 'L', M, N, A, LDA, U, LDU )
               IF( WNTUS ) &
                 NCU = N
               IF( WNTUA ) &
                 NCU = M
               CALL DORGBR( 'Q', M, NCU, N, U, LDU, WORK( ITAUQ ), &
                           WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTVAS ) THEN

               CALL DLACPY( 'U', N, N, A, LDA, VT, LDVT )
               CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), &
                           WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTUO ) THEN

               CALL DORGBR( 'Q', M, N, N, A, LDA, WORK( ITAUQ ), &
                           WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTVO ) THEN

               CALL DORGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ), &
                           WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IWORK = IE + N
            IF( WNTUAS .OR. WNTUO ) &
              NRU = M
            IF( WNTUN ) &
              NRU = 0
            IF( WNTVAS .OR. WNTVO ) &
              NCVT = N
            IF( WNTVN ) &
              NCVT = 0
            IF( ( .NOT.WNTUO ) .AND. ( .NOT.WNTVO ) ) THEN

               CALL DBDSQR( 'U', N, NCVT, NRU, 0, S, WORK( IE ), VT, &
                           LDVT, U, LDU, DUM, 1, WORK( IWORK ), INFO )
            ELSE IF( ( .NOT.WNTUO ) .AND. WNTVO ) THEN

               CALL DBDSQR( 'U', N, NCVT, NRU, 0, S, WORK( IE ), A, LDA, &
                           U, LDU, DUM, 1, WORK( IWORK ), INFO )
            ELSE

               CALL DBDSQR( 'U', N, NCVT, NRU, 0, S, WORK( IE ), VT, &
                           LDVT, A, LDA, DUM, 1, WORK( IWORK ), INFO )
            END IF

         END IF

      ELSE

         IF( N.GE.MNTHR ) THEN

            IF( WNTVN ) THEN

               ITAU = 1
               IWORK = ITAU + M
!*
!*              Compute A=L*Q
!*              (Workspace: need 2*M, prefer M + M*NB)
!*
               CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                           LWORK-IWORK+1, IERR )
!*
!*              Zero out above L
!*
               CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA )
               IE = 1
               ITAUQ = IE + M
               ITAUP = ITAUQ + M
               IWORK = ITAUP + M
!*
!*              Bidiagonalize L in A
!*              (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!*
               CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                           WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, &
                           IERR )
               IF( WNTUO .OR. WNTUAS ) THEN
!*
!*                 If left singular vectors desired, generate Q
!                 (Workspace: need 4*M, prefer 3*M + M*NB)
!
                  CALL DORGBR( 'Q', M, M, M, A, LDA, WORK( ITAUQ ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
               END IF
               IWORK = IE + M
               NRU = 0
               IF( WNTUO .OR. WNTUAS ) &
                 NRU = M
!*
!*              Perform bidiagonal QR iteration, computing left singular
!*              vectors of A in A if desired
!*              (Workspace: need BDSPAC)
!*
               CALL DBDSQR( 'U', M, 0, NRU, 0, S, WORK( IE ), DUM, 1, A, &
                           LDA, DUM, 1, WORK( IWORK ), INFO )
!*
!*              If left singular vectors desired in U, copy them there
!*
               IF( WNTUAS ) &
                 CALL DLACPY( 'F', M, M, A, LDA, U, LDU )
!*
            ELSE IF( WNTVO .AND. WNTUN ) THEN
!*
!*              Path 2t(N much larger than M, JOBU='N', JOBVT='O')
!*              M right singular vectors to be overwritten on A and
!*              no left singular vectors to be computed
!*
               IF( LWORK.GE.M*M+MAX( 4*M, BDSPAC ) ) THEN
!*
!*                 Sufficient workspace for a fast algorithm
!*
                  IR = 1
                  IF( LWORK.GE.MAX( WRKBL, LDA*N + M ) + LDA*M ) THEN
!*
!*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
!*
                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = LDA
                  ELSE IF( LWORK.GE.MAX( WRKBL, LDA*N + M ) + M*M ) THEN
!*
!*                    WORK(IU) is LDA by N and WORK(IR) is M by M
!*
                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = M
                  ELSE
!*
!*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M
!*
                     LDWRKU = M
                     CHUNK = ( LWORK-M*M-M ) / M
                     LDWRKR = M
                  END IF
                  ITAU = IR + LDWRKR*M
                  IWORK = ITAU + M
!*
!*                 Compute A=L*Q
!*                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!*
                  CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                 Copy L to WORK(IR) and zero out above it
!*
                  CALL DLACPY( 'L', M, M, A, LDA, WORK( IR ), LDWRKR )
                  CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, &
                              WORK( IR+LDWRKR ), LDWRKR )
!*
!*                 Generate Q in A
!*                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!*
                  CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = ITAU
                  ITAUQ = IE + M
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M
!*
!*                 Bidiagonalize L in WORK(IR)
!*                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
!*
                  CALL DGEBRD( M, M, WORK( IR ), LDWRKR, S, WORK( IE ), &
                              WORK( ITAUQ ), WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                 Generate right vectors bidiagonalizing L
!*                 (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB)
!*
                  CALL DORGBR( 'P', M, M, M, WORK( IR ), LDWRKR, &
                              WORK( ITAUP ), WORK( IWORK ), &
                              LWORK-IWORK+1, IERR )
                  IWORK = IE + M
!*
!*                 Perform bidiagonal QR iteration, computing right
!*                 singular vectors of L in WORK(IR)
!*                 (Workspace: need M*M + BDSPAC)
!*
                  CALL DBDSQR( 'U', M, M, 0, 0, S, WORK( IE ), &
                              WORK( IR ), LDWRKR, DUM, 1, DUM, 1, &
                              WORK( IWORK ), INFO )
                  IU = IE + M
!*
!*                 Multiply right singular vectors of L in WORK(IR) by Q
!*                 in A, storing result in WORK(IU) and copying to A
!*                 (Workspace: need M*M + 2*M, prefer M*M + M*N + M)
!*
                  DO 30 I = 1, N, CHUNK
                     BLK = MIN( N-I+1, CHUNK )
                     CALL DGEMM( 'N', 'N', M, BLK, M, ONE, WORK( IR ), &
                                LDWRKR, A( 1, I ), LDA, ZERO, &
                                WORK( IU ), LDWRKU )
                     CALL DLACPY( 'F', M, BLK, WORK( IU ), LDWRKU, &
                                 A( 1, I ), LDA )
   30             CONTINUE
!*
               ELSE
!*
!*                 Insufficient workspace for a fast algorithm
!*
                  IE = 1
                  ITAUQ = IE + M
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M

                  CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), &
                              WORK( ITAUQ ), WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )

                  CALL DORGBR( 'P', M, N, M, A, LDA, WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IWORK = IE + M

                  CALL DBDSQR( 'L', M, N, 0, 0, S, WORK( IE ), A, LDA, &
                              DUM, 1, DUM, 1, WORK( IWORK ), INFO )

               END IF

            ELSE IF( WNTVO .AND. WNTUAS ) THEN

               IF( LWORK.GE.M*M+MAX( 4*M, BDSPAC ) ) THEN

                  IR = 1
                  IF( LWORK.GE.MAX( WRKBL, LDA*N + M ) + LDA*M ) THEN

                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = LDA
                  ELSE IF( LWORK.GE.MAX( WRKBL, LDA*N + M ) + M*M ) THEN

                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = M
                  ELSE

                     LDWRKU = M
                     CHUNK = ( LWORK-M*M-M ) / M
                     LDWRKR = M
                  END IF
                  ITAU = IR + LDWRKR*M
                  IWORK = ITAU + M

                  CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )

                  CALL DLACPY( 'L', M, M, A, LDA, U, LDU ) &
                  CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, U( 1, 2 ),
                              LDU )
!*
!*                 Generate Q in A
!*                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!*
                  CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = ITAU
                  ITAUQ = IE + M
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M
!*
!*                 Bidiagonalize L in U, copying result to WORK(IR)
!*                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
!*
                  CALL DGEBRD( M, M, U, LDU, S, WORK( IE ), &
                              WORK( ITAUQ ), WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR ) 
                  CALL DLACPY( 'U', M, M, U, LDU, WORK( IR ), LDWRKR )
!*
!*                 Generate right vectors bidiagonalizing L in WORK(IR)
!*                 (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB)
!*
                  CALL DORGBR( 'P', M, M, M, WORK( IR ), LDWRKR, &
                              WORK( ITAUP ), WORK( IWORK ), &
                              LWORK-IWORK+1, IERR )
!*
!*                 Generate left vectors bidiagonalizing L in U
!*                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB)
!*
                  CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IWORK = IE + M
!*
!*                 Perform bidiagonal QR iteration, computing left
!*                 singular vectors of L in U, and computing right
!*                 singular vectors of L in WORK(IR)
!*                 (Workspace: need M*M + BDSPAC)
!*
                  CALL DBDSQR( 'U', M, M, M, 0, S, WORK( IE ), &
                              WORK( IR ), LDWRKR, U, LDU, DUM, 1, &
                             WORK( IWORK ), INFO )
                  IU = IE + M
!*
!*                 Multiply right singular vectors of L in WORK(IR) by Q
!*                 in A, storing result in WORK(IU) and copying to A
!*                 (Workspace: need M*M + 2*M, prefer M*M + M*N + M))
!*
                  DO 40 I = 1, N, CHUNK
                     BLK = MIN( N-I+1, CHUNK )
                     CALL DGEMM( 'N', 'N', M, BLK, M, ONE, WORK( IR ), &
                                LDWRKR, A( 1, I ), LDA, ZERO, &
                                WORK( IU ), LDWRKU )
                     CALL DLACPY( 'F', M, BLK, WORK( IU ), LDWRKU, &
                                 A( 1, I ), LDA )
   40             CONTINUE
!*
               ELSE
!*
!*                 Insufficient workspace for a fast algorithm
!*
                  ITAU = 1
                  IWORK = ITAU + M
!*
!*                 Compute A=L*Q
!*                 (Workspace: need 2*M, prefer M + M*NB)
!*
                  CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                 Copy L to U, zeroing out above it
!*
                  CALL DLACPY( 'L', M, M, A, LDA, U, LDU )
                  CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, U( 1, 2 ), &
                              LDU )
!*
!*                 Generate Q in A
!*                 (Workspace: need 2*M, prefer M + M*NB)
!*
                  CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = ITAU
                  ITAUQ = IE + M
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M
!*
!*                 Bidiagonalize L in U
!*                 (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!*
                  CALL DGEBRD( M, M, U, LDU, S, WORK( IE ), &
                              WORK( ITAUQ ), WORK( ITAUP ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                 Multiply right vectors bidiagonalizing L by Q in A
!*                 (Workspace: need 3*M + N, prefer 3*M + N*NB)
!*
                  CALL DORMBR( 'P', 'L', 'T', M, N, M, U, LDU, &
                              WORK( ITAUP ), A, LDA, WORK( IWORK ), &
                              LWORK-IWORK+1, IERR )
!*
!*                 Generate left vectors bidiagonalizing L in U
!*                 (Workspace: need 4*M, prefer 3*M + M*NB)
!*
                  CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), &
                              WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IWORK = IE + M
!*
!*                 Perform bidiagonal QR iteration, computing left
!*                 singular vectors of A in U and computing right
!*                 singular vectors of A in A
!*                 (Workspace: need BDSPAC)
!*
                  CALL DBDSQR( 'U', M, N, M, 0, S, WORK( IE ), A, LDA, &
                              U, LDU, DUM, 1, WORK( IWORK ), INFO )

               END IF

            ELSE IF( WNTVS ) THEN

               IF( WNTUN ) THEN

!*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S')
!*                 M right singular vectors to be computed in VT and
!*                 no left singular vectors to be computed
!*
                  IF( LWORK.GE.M*M+MAX( 4*M, BDSPAC ) ) THEN
!*
!!*                    Sufficient workspace for a fast algorithm
!*
                     IR = 1
                     IF( LWORK.GE.WRKBL+LDA*M ) THEN
!*
!*                       WORK(IR) is LDA by M
!*
                        LDWRKR = LDA
                     ELSE
!*
!*                       WORK(IR) is M by M
!*
                        LDWRKR = M
                     END IF
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M
!*
!*                    Compute A=L*Q
!*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!*
                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Copy L to WORK(IR), zeroing out above it
!*
                     CALL DLACPY( 'L', M, M, A, LDA, WORK( IR ), &
                                 LDWRKR )
                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, &
                                 WORK( IR+LDWRKR ), LDWRKR )
!*
!*                    Generate Q in A
!*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!*
                     CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
!*
!*                    Bidiagonalize L in WORK(IR)
!*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
!*
                     CALL DGEBRD( M, M, WORK( IR ), LDWRKR, S, &
                                 WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
!*
!*                    Generate right vectors bidiagonalizing L in
!*                    WORK(IR)
!*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB)
!*
                     CALL DORGBR( 'P', M, M, M, WORK( IR ), LDWRKR, &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     IWORK = IE + M
!*
!*                    Perform bidiagonal QR iteration, computing right
!*                    singular vectors of L in WORK(IR)
!*                    (Workspace: need M*M + BDSPAC)
!*
                     CALL DBDSQR( 'U', M, M, 0, 0, S, WORK( IE ), &
                                 WORK( IR ), LDWRKR, DUM, 1, DUM, 1, &
                                 WORK( IWORK ), INFO )
!*
!*                    Multiply right singular vectors of L in WORK(IR) by
!*                    Q in A, storing result in VT
!*                    (Workspace: need M*M)
!*
                     CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IR ), &
                                LDWRKR, A, LDA, ZERO, VT, LDVT )
!*
                  ELSE
!*
!!*                    Insufficient workspace for a fast algorithm
!*
                     ITAU = 1
                     IWORK = ITAU + M
!*
!*                    Compute A=L*Q
!*                    (Workspace: need 2*M, prefer M + M*NB)
!*
                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Copy result to VT
!*
                     CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )
!*
!*                    Generate Q in VT
!*                    (Workspace: need 2*M, prefer M + M*NB)
!*
                     CALL DORGLQ( M, N, M, VT, LDVT, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
!*
!*                    Zero out above L in A
!*
                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), &
                                 LDA )
!*
!*                    Bidiagonalize L in A
!*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!*
                     CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), &
                                 WORK( ITAUQ ), WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Multiply right vectors bidiagonalizing L by Q in VT
!*                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
!*
                     CALL DORMBR( 'P', 'L', 'T', M, N, M, A, LDA, &
                                 WORK( ITAUP ), VT, LDVT, &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + M
!*
!*                    Perform bidiagonal QR iteration, computing right
!*                    singular vectors of A in VT
!*                    (Workspace: need BDSPAC)
!*
                     CALL DBDSQR( 'U', M, N, 0, 0, S, WORK( IE ), VT, &
                                 LDVT, DUM, 1, DUM, 1, WORK( IWORK ), &
                                 INFO )
!*
                  END IF
!*
               ELSE IF( WNTUO ) THEN
!*
!*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
!*                 M right singular vectors to be computed in VT and
!*                 M left singular vectors to be overwritten on A
!*
                  IF( LWORK.GE.2*M*M+MAX( 4*M, BDSPAC ) ) THEN
!*
!*                    Sufficient workspace for a fast algorithm
!*
                     IU = 1
                     IF( LWORK.GE.WRKBL+2*LDA*M ) THEN
!!*
!*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M
!*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = LDA
                     ELSE IF( LWORK.GE.WRKBL+( LDA + M )*M ) THEN
!*
!*                       WORK(IU) is LDA by M and WORK(IR) is M by M
!*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     ELSE
!*
!*                       WORK(IU) is M by M and WORK(IR) is M by M
!*
                        LDWRKU = M
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     END IF
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M
!*
!*                    Compute A=L*Q
!*                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)
!*
                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Copy L to WORK(IU), zeroing out below it
!*
                     CALL DLACPY( 'L', M, M, A, LDA, WORK( IU ), &
                                 LDWRKU )
                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, &
                                 WORK( IU+LDWRKU ), LDWRKU )
!*
!*                    Generate Q in A
!*                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)
!*
                     CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
!*
!*                    Bidiagonalize L in WORK(IU), copying result to
!*                    WORK(IR)
!*                    (Workspace: need 2*M*M + 4*M,
!*                                prefer 2*M*M+3*M+2*M*NB)
!*
                     CALL DGEBRD( M, M, WORK( IU ), LDWRKU, S, &
                                 WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, M, WORK( IU ), LDWRKU, &
                                 WORK( IR ), LDWRKR )
!*
!*                    Generate right bidiagonalizing vectors in WORK(IU)
!*                    (Workspace: need 2*M*M + 4*M-1,
!*                                prefer 2*M*M+3*M+(M-1)*NB)
!*
                     CALL DORGBR( 'P', M, M, M, WORK( IU ), LDWRKU, &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
!*
!*                    Generate left bidiagonalizing vectors in WORK(IR)
!*                    (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB)
!*
                     CALL DORGBR( 'Q', M, M, M, WORK( IR ), LDWRKR, &
                                 WORK( ITAUQ ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     IWORK = IE + M
!*
!*                    Perform bidiagonal QR iteration, computing left
!*                    singular vectors of L in WORK(IR) and computing
!*                    right singular vectors of L in WORK(IU)
!*                    (Workspace: need 2*M*M + BDSPAC)
!*
                     CALL DBDSQR( 'U', M, M, M, 0, S, WORK( IE ), &
                                 WORK( IU ), LDWRKU, WORK( IR ), &
                                 LDWRKR, DUM, 1, WORK( IWORK ), INFO )
!*
!*                    Multiply right singular vectors of L in WORK(IU) by
!*                    Q in A, storing result in VT
!*                    (Workspace: need M*M)
!*
                     CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IU ), &
                                LDWRKU, A, LDA, ZERO, VT, LDVT )
!*
!*                    Copy left singular vectors of L to A
!*                    (Workspace: need M*M)
!*
                     CALL DLACPY( 'F', M, M, WORK( IR ), LDWRKR, A, &
                                 LDA )
!*
                  ELSE
!*
!*                    Insufficient workspace for a fast algorithm
!*
                     ITAU = 1
                     IWORK = ITAU + M
!*
!*                    Compute A=L*Q, copying result to VT
!*                    (Workspace: need 2*M, prefer M + M*NB)
!*
                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                     CALL DORGLQ( M, N, M, VT, LDVT, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), &
                                 LDA )

                     CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), &
                                 WORK( ITAUQ ), WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DORMBR( 'P', 'L', 'T', M, N, M, A, LDA, &
                                 WORK( ITAUP ), VT, LDVT, &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DORGBR( 'Q', M, M, M, A, LDA, WORK( ITAUQ ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + M

                     CALL DBDSQR( 'U', M, N, M, 0, S, WORK( IE ), VT, &
                                 LDVT, A, LDA, DUM, 1, WORK( IWORK ), &
                                INFO )

                  END IF

               ELSE IF( WNTUAS ) THEN

                  IF( LWORK.GE.M*M+MAX( 4*M, BDSPAC ) ) THEN

                     IU = 1
                     IF( LWORK.GE.WRKBL+LDA*M ) THEN

                        LDWRKU = LDA
                     ELSE

                        LDWRKU = M
                     END IF
                     ITAU = IU + LDWRKU*M
                     IWORK = ITAU + M

                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DLACPY( 'L', M, M, A, LDA, WORK( IU ), &
                                 LDWRKU )
                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, &
                                 WORK( IU+LDWRKU ), LDWRKU )

                     CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     CALL DGEBRD( M, M, WORK( IU ), LDWRKU, S, &
                                 WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, M, WORK( IU ), LDWRKU, U, &
                                 LDU )

                     CALL DORGBR( 'P', M, M, M, WORK( IU ), LDWRKU, &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )

                     CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + M

                     CALL DBDSQR( 'U', M, M, M, 0, S, WORK( IE ), &
                                 WORK( IU ), LDWRKU, U, LDU, DUM, 1, &
                                 WORK( IWORK ), INFO )

                     CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IU ), &
                                LDWRKU, A, LDA, ZERO, VT, LDVT )

                  ELSE

                     ITAU = 1
                     IWORK = ITAU + M

                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ),
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                     CALL DORGLQ( M, N, M, VT, LDVT, WORK( ITAU ),
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DLACPY( 'L', M, M, A, LDA, U, LDU )
                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, U( 1, 2 ),
                                 LDU )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     CALL DGEBRD( M, M, U, LDU, S, WORK( IE ),
                                 WORK( ITAUQ ), WORK( ITAUP ),
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DORMBR( 'P', 'L', 'T', M, N, M, U, LDU, &
                                 WORK( ITAUP ), VT, LDVT, &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + M

                     CALL DBDSQR( 'U', M, N, M, 0, S, WORK( IE ), VT, &
                                 LDVT, U, LDU, DUM, 1, WORK( IWORK ), &
                                 INFO )

                  END IF

               END IF

            ELSE IF( WNTVA ) THEN

               IF( WNTUN ) THEN

                  IF( LWORK.GE.M*M+MAX( N + M, 4*M, BDSPAC ) ) THEN

                     IR = 1
                     IF( LWORK.GE.WRKBL+LDA*M ) THEN

                        LDWRKR = LDA
                     ELSE

                        LDWRKR = M
                     END IF
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M

                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                     CALL DLACPY( 'L', M, M, A, LDA, WORK( IR ), &
                                 LDWRKR )
                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, &
                                 WORK( IR+LDWRKR ), LDWRKR )

                     CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     CALL DGEBRD( M, M, WORK( IR ), LDWRKR, S, &
                                 WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )

                     CALL DORGBR( 'P', M, M, M, WORK( IR ), LDWRKR, &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     IWORK = IE + M

                     CALL DBDSQR( 'U', M, M, 0, 0, S, WORK( IE ), &
                                 WORK( IR ), LDWRKR, DUM, 1, DUM, 1, &
                                 WORK( IWORK ), INFO )

                     CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IR ), &
                                LDWRKR, VT, LDVT, ZERO, A, LDA )

                     CALL DLACPY( 'F', M, N, A, LDA, VT, LDVT )

                  ELSE

                     ITAU = 1
                     IWORK = ITAU + M

                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                     CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), &
                                 LDA )

                     CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), &
                                 WORK( ITAUQ ), WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )

                     CALL DORMBR( 'P', 'L', 'T', M, N, M, A, LDA, &
                                 WORK( ITAUP ), VT, LDVT, &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + M

                     CALL DBDSQR( 'U', M, N, 0, 0, S, WORK( IE ), VT, &
                                 LDVT, DUM, 1, DUM, 1, WORK( IWORK ), &
                                 INFO )
!*
                  END IF
!*
               ELSE IF( WNTUO ) THEN
!*
!*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
!*                 N right singular vectors to be computed in VT and
!*                 M left singular vectors to be overwritten on A
!*
                  IF( LWORK.GE.2*M*M+MAX( N + M, 4*M, BDSPAC ) ) THEN
!*
!*                    Sufficient workspace for a fast algorithm
!*
                     IU = 1
                     IF( LWORK.GE.WRKBL+2*LDA*M ) THEN
!*
!*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M
!*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = LDA
                     ELSE IF( LWORK.GE.WRKBL+( LDA + M )*M ) THEN
!*
!*                       WORK(IU) is LDA by M and WORK(IR) is M by M
!*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     ELSE
!*
!*                       WORK(IU) is M by M and WORK(IR) is M by M
!*
                        LDWRKU = M
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     END IF
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M
!*
!*                    Compute A=L*Q, copying result to VT
!*                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)
!*
                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )
!*
!*                    Generate Q in VT
!*                    (Workspace: need 2*M*M + M + N, prefer 2*M*M + M + N*NB)
!*
                     CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Copy L to WORK(IU), zeroing out above it
!*
                     CALL DLACPY( 'L', M, M, A, LDA, WORK( IU ), &
                                 LDWRKU )
                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, &
                                 WORK( IU+LDWRKU ), LDWRKU )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
!*
!*                    Bidiagonalize L in WORK(IU), copying result to
!*                    WORK(IR)
!*                    (Workspace: need 2*M*M + 4*M,
!*                                prefer 2*M*M+3*M+2*M*NB)
!*
                     CALL DGEBRD( M, M, WORK( IU ), LDWRKU, S, &
                                 WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, M, WORK( IU ), LDWRKU, &
                                 WORK( IR ), LDWRKR )
!*
!*                    Generate right bidiagonalizing vectors in WORK(IU)
!*                    (Workspace: need 2*M*M + 4*M-1,
!*                                prefer 2*M*M+3*M+(M-1)*NB)
!*
                     CALL DORGBR( 'P', M, M, M, WORK( IU ), LDWRKU, &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
!*
!*                    Generate left bidiagonalizing vectors in WORK(IR)
!*                    (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB)
!*
                     CALL DORGBR( 'Q', M, M, M, WORK( IR ), LDWRKR, &
                                 WORK( ITAUQ ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     IWORK = IE + M
!*
!*                    Perform bidiagonal QR iteration, computing left
!*                    singular vectors of L in WORK(IR) and computing
!*                    right singular vectors of L in WORK(IU)
!*                    (Workspace: need 2*M*M + BDSPAC)
!*
                     CALL DBDSQR( 'U', M, M, M, 0, S, WORK( IE ), &
                                 WORK( IU ), LDWRKU, WORK( IR ), &
                                 LDWRKR, DUM, 1, WORK( IWORK ), INFO )
!*
!*!                    Multiply right singular vectors of L in WORK(IU) by
!*                    Q in VT, storing result in A
!*                    (Workspace: need M*M)
!*
                     CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IU ), &
                                LDWRKU, VT, LDVT, ZERO, A, LDA )
!*
!*                    Copy right singular vectors of A from A to VT
!*
                     CALL DLACPY( 'F', M, N, A, LDA, VT, LDVT )
!*
!*                    Copy left singular vectors of A from WORK(IR) to A
!*
                     CALL DLACPY( 'F', M, M, WORK( IR ), LDWRKR, A, &
                                 LDA )
!*
                  ELSE
!*
!*                    Insufficient workspace for a fast algorithm
!*
                     ITAU = 1
                     IWORK = ITAU + M
!*
!*                    Compute A=L*Q, copying result to VT
!*                    (Workspace: need 2*M, prefer M + M*NB)
!*
                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )
!*
!*                    Generate Q in VT
!*                    (Workspace: need M + N, prefer M + N*NB)
!*
                     CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
!*
!*                    Zero out above L in A
!*
                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), &
                                 LDA )
!*
!*                    Bidiagonalize L in A
!*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!*
                     CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), &
                                 WORK( ITAUQ ), WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Multiply right bidiagonalizing vectors in A by Q
!*                    in VT
!*                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
!*
                     CALL DORMBR( 'P', 'L', 'T', M, N, M, A, LDA, &
                                 WORK( ITAUP ), VT, LDVT, &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Generate left bidiagonalizing vectors in A
!*                    (Workspace: need 4*M, prefer 3*M + M*NB)
!*
                     CALL DORGBR( 'Q', M, M, M, A, LDA, WORK( ITAUQ ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + M
!*
!!*                    Perform bidiagonal QR iteration, computing left
!*                    singular vectors of A in A and computing right
!*                    singular vectors of A in VT
!*                    (Workspace: need BDSPAC)
!*
                     CALL DBDSQR( 'U', M, N, M, 0, S, WORK( IE ), VT, &
                                 LDVT, A, LDA, DUM, 1, WORK( IWORK ), &
                                 INFO )
!*
                  END IF
!*
               ELSE IF( WNTUAS ) THEN
!*
!*                 Path 9t(N much larger than M, JOBU='S' or 'A',
!*                         JOBVT='A')
!*                 N right singular vectors to be computed in VT and
!*                 M left singular vectors to be computed in U
!*
                  IF( LWORK.GE.M*M+MAX( N + M, 4*M, BDSPAC ) ) THEN
!*
!*                    Sufficient workspace for a fast algorithm
!*
                     IU = 1
                     IF( LWORK.GE.WRKBL+LDA*M ) THEN
!*
!*                       WORK(IU) is LDA by M
!*
                        LDWRKU = LDA
                     ELSE
!*
!*                       WORK(IU) is M by M
!*
                        LDWRKU = M
                     END IF
                     ITAU = IU + LDWRKU*M
                     IWORK = ITAU + M
!*
!*                    Compute A=L*Q, copying result to VT
!*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
!*
                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )
!*
!*                    Generate Q in VT
!*                    (Workspace: need M*M + M + N, prefer M*M + M + N*NB)
!*
                     CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Copy L to WORK(IU), zeroing out above it
!*
                     CALL DLACPY( 'L', M, M, A, LDA, WORK( IU ), &
                                 LDWRKU )
                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO,
                                 WORK( IU+LDWRKU ), LDWRKU )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
!*
!*                    Bidiagonalize L in WORK(IU), copying result to U
!*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
!*
                     CALL DGEBRD( M, M, WORK( IU ), LDWRKU, S, &
                                 WORK( IE ), WORK( ITAUQ ), &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'L', M, M, WORK( IU ), LDWRKU, U, &
                                 LDU )
!*
!*                    Generate right bidiagonalizing vectors in WORK(IU)
!*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB)
!*
                     CALL DORGBR( 'P', M, M, M, WORK( IU ), LDWRKU, &
                                 WORK( ITAUP ), WORK( IWORK ), &
                                 LWORK-IWORK+1, IERR )
!*
!*                    Generate left bidiagonalizing vectors in U
!*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB)
!*
                     CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + M
!*
!*                    Perform bidiagonal QR iteration, computing left
!*                    singular vectors of L in U and computing right
!*                    singular vectors of L in WORK(IU)
!*                    (Workspace: need M*M + BDSPAC)
!*
                     CALL DBDSQR( 'U', M, M, M, 0, S, WORK( IE ), &
                                 WORK( IU ), LDWRKU, U, LDU, DUM, 1, &
                                 WORK( IWORK ), INFO )
!*
!*                    Multiply right singular vectors of L in WORK(IU) by
!*                    Q in VT, storing result in A
!*                    (Workspace: need M*M)
!*
                     CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IU ), &
                                LDWRKU, VT, LDVT, ZERO, A, LDA )
!*
!*                    Copy right singular vectors of A from A to VT
!*
                     CALL DLACPY( 'F', M, N, A, LDA, VT, LDVT )
!*
                  ELSE
!*
!*                    Insufficient workspace for a fast algorithm
!*
                     ITAU = 1
                     IWORK = ITAU + M
!*
!*                    Compute A=L*Q, copying result to VT
!*                    (Workspace: need 2*M, prefer M + M*NB)
!*
                     CALL DGELQF( M, N, A, LDA, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )
!*
!*                    Generate Q in VT
!*                    (Workspace: need M + N, prefer M + N*NB)
!*
                     CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Copy L to U, zeroing out above it
!*
                     CALL DLACPY( 'L', M, M, A, LDA, U, LDU )
                     CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, U( 1, 2 ), &
                                  LDU )
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
!*
!*                    Bidiagonalize L in U
!*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
!*
                     CALL DGEBRD( M, M, U, LDU, S, WORK( IE ), &
                                 WORK( ITAUQ ), WORK( ITAUP ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Multiply right bidiagonalizing vectors in U by Q
!*                    in VT
!*                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
!*
                     CALL DORMBR( 'P', 'L', 'T', M, N, M, U, LDU, &
                                 WORK( ITAUP ), VT, LDVT, &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
!*
!*                    Generate left bidiagonalizing vectors in U
!*                    (Workspace: need 4*M, prefer 3*M + M*NB)
!*
                     CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), &
                                 WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IWORK = IE + M
!*
!*                    Perform bidiagonal QR iteration, computing left
!*                    singular vectors of A in U and computing right
!*                    singular vectors of A in VT
!*                    (Workspace: need BDSPAC)
!*
                     CALL DBDSQR( 'U', M, N, M, 0, S, WORK( IE ), VT, &
                                LDVT, U, LDU, DUM, 1, WORK( IWORK ), &
                                INFO )
!*
                  END IF
!*
               END IF
!*
            END IF
!*
         ELSE
!*
!*           N .LT. MNTHR
!*
!*           Path 10t(N greater than M, but not much larger)
!*           Reduce to bidiagonal form without LQ decomposition
!*
            IE = 1
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
!*
!*           Bidiagonalize A
!*           (Workspace: need 3*M + N, prefer 3*M + (M + N)*NB)
!*
            CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                        WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, &
                        IERR )
            IF( WNTUAS ) THEN
!*
!*              If left singular vectors desired in U, copy result to U
!*              and generate left bidiagonalizing vectors in U
!!*              (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB)
!*
               CALL DLACPY( 'L', M, M, A, LDA, U, LDU )
               CALL DORGBR( 'Q', M, M, N, U, LDU, WORK( ITAUQ ), &
                            WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTVAS ) THEN
!*
!*              If right singular vectors desired in VT, copy result to
!*              VT and generate right bidiagonalizing vectors in VT
!*              (Workspace: need 3*M + NRVT, prefer 3*M + NRVT*NB)
!*
               CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )
               IF( WNTVA ) &
                 NRVT = N
               IF( WNTVS ) &
                 NRVT = M
               CALL DORGBR( 'P', NRVT, N, M, VT, LDVT, WORK( ITAUP ), &
                            WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTUO ) THEN
!*
!*              If left singular vectors desired in A, generate left
!*              bidiagonalizing vectors in A
!*              (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB)
!*
               CALL DORGBR( 'Q', M, M, N, A, LDA, WORK( ITAUQ ), &
                           WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTVO ) THEN
!*
!*              If right singular vectors desired in A, generate right
!*              bidiagonalizing vectors in A
!*              (Workspace: need 4*M, prefer 3*M + M*NB)
!*
               CALL DORGBR( 'P', M, N, M, A, LDA, WORK( ITAUP ), &
                           WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IWORK = IE + M
            IF( WNTUAS .OR. WNTUO ) &
              NRU = M
            IF( WNTUN ) &
              NRU = 0
            IF( WNTVAS .OR. WNTVO ) &
              NCVT = N
            IF( WNTVN ) &
              NCVT = 0
            IF( ( .NOT.WNTUO ) .AND. ( .NOT.WNTVO ) ) THEN
!*
!*              Perform bidiagonal QR iteration, if desired, computing
!*              left singular vectors in U and computing right singular
!*              vectors in VT
!*              (Workspace: need BDSPAC)
!*
               CALL DBDSQR( 'L', M, NCVT, NRU, 0, S, WORK( IE ), VT, &
                          LDVT, U, LDU, DUM, 1, WORK( IWORK ), INFO )
            ELSE IF( ( .NOT.WNTUO ) .AND. WNTVO ) THEN
!*
!*              Perform bidiagonal QR iteration, if desired, computing
!*              left singular vectors in U and computing right singular
!*              vectors in A
!*              (Workspace: need BDSPAC)
!*
               CALL DBDSQR( 'L', M, NCVT, NRU, 0, S, WORK( IE ), A, LDA, &
                           U, LDU, DUM, 1, WORK( IWORK ), INFO )
            ELSE
!*
!*              Perform bidiagonal QR iteration, if desired, computing
!*              left singular vectors in A and computing right singular
!*              vectors in VT
!*              (Workspace: need BDSPAC)
!*
               CALL DBDSQR( 'L', M, NCVT, NRU, 0, S, WORK( IE ), VT, &
                           LDVT, A, LDA, DUM, 1, WORK( IWORK ), INFO )
            END IF
!*
         END IF
!*
      END IF
!*
!*     If DBDSQR failed to converge, copy unconverged superdiagonals
!*     to WORK( 2:MINMN )
!*
      IF( INFO.NE.0 ) THEN
         IF( IE.GT.2 ) THEN
            DO 50 I = 1, MINMN - 1
               WORK( I+1 ) = WORK( I+IE-1 )
   50       CONTINUE
         END IF
         IF( IE.LT.2 ) THEN
            DO 60 I = MINMN - 1, 1, -1
               WORK( I+1 ) = WORK( I+IE-1 )
   60       CONTINUE
         END IF
      END IF
!*
!*     Undo scaling if necessary
!*
      IF( ISCL.EQ.1 ) THEN
         IF( ANRM.GT.BIGNUM )
           CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN,
                        IERR )
         IF( INFO.NE.0 .AND. ANRM.GT.BIGNUM ) &
           CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN-1, 1, WORK( 2 ), &
                        MINMN, IERR )
         IF( ANRM.LT.SMLNUM ) &
           CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, &
                        IERR )
         IF( INFO.NE.0 .AND. ANRM.LT.SMLNUM ) &
           CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN-1, 1, WORK( 2 ), &
                        MINMN, IERR )
      END IF
!*
!*     Return optimal workspace in WORK(1)
!*
      WORK( 1 ) = MAXWRK


END SUBROUTINE

#if 0
*  ===========
*>
*> \verbatim
*>  Bug report from Cezary Dendek.
*>  On March 23rd 2017, the INTEGER variable MAXIT = MAXITR*N**2 is
*>  removed since it can overflow pretty easily (for N larger or equal
*>  than 18,919). We instead use MAXITDIVN = MAXITR*N.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2017
*
*> \ingroup auxOTHERcomputational
*
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, &
     LDU, C, LDC, WORK, INFO ) !GCC$ ATTRIBUTES hot :: DBDSQR !GCC$ ATTRIBUTES aligned(32) :: DBDSQR
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   NEGONE
      PARAMETER          ( NEGONE = -1.0D0 )
      DOUBLE PRECISION   HNDRTH
      PARAMETER          ( HNDRTH = 0.01D0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 10.0D0 )
      DOUBLE PRECISION   HNDRD
      PARAMETER          ( HNDRD = 100.0D0 )
      DOUBLE PRECISION   MEIGTH
      PARAMETER          ( MEIGTH = -0.125D0 )
      INTEGER            MAXITR
      PARAMETER          ( MAXITR = 6 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LOWER, ROTATE
      INTEGER            I, IDIR, ISUB, ITER, ITERDIVN, J, LL, LLL, M, &
                        MAXITDIVN, NM1, NM12, NM13, OLDLL, OLDM
      DOUBLE PRECISION   ABSE, ABSS, COSL, COSR, CS, EPS, F, G, H, MU, &
                        OLDCS, OLDSN, R, SHIFT, SIGMN, SIGMX, SINL, &
                        SINR, SLL, SMAX, SMIN, SMINL, SMINOA, &
                        SN, THRESH, TOL, TOLMUL, UNFL
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DROT,DSCAL,DSWAP
                       
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SIGN, SQRT
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LOWER ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NCVT.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NCC.LT.0 ) THEN
         INFO = -5
      ELSE IF( ( NCVT.EQ.0 .AND. LDVT.LT.1 ) .OR. &
              ( NCVT.GT.0 .AND. LDVT.LT.MAX( 1, N ) ) ) THEN
         INFO = -9
      ELSE IF( LDU.LT.MAX( 1, NRU ) ) THEN
         INFO = -11
      ELSE IF( ( NCC.EQ.0 .AND. LDC.LT.1 ) .OR. &
              ( NCC.GT.0 .AND. LDC.LT.MAX( 1, N ) ) ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DBDSQR', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 ) &
        RETURN
      IF( N.EQ.1 ) &
        GO TO 160
!*
!*     ROTATE is true if any singular vectors desired, false otherwise
!*
      ROTATE = ( NCVT.GT.0 ) .OR. ( NRU.GT.0 ) .OR. ( NCC.GT.0 )
!*
!*     If no singular vectors desired, use qd algorithm
!*
      IF( .NOT.ROTATE ) THEN
         CALL DLASQ1( N, D, E, WORK, INFO )
!*
!*     If INFO equals 2, dqds didn't finish, try to finish
!*
         IF( INFO .NE. 2 ) RETURN
         INFO = 0
      END IF
!*
      NM1 = N - 1
      NM12 = NM1 + NM1
      NM13 = NM12 + NM1
      IDIR = 0
!*
!*     Get machine constants
!*
      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
!*
!*     If matrix lower bidiagonal, rotate to be upper bidiagonal
!*     by applying Givens rotations on the left
!*
      IF( LOWER ) THEN
         DO 10 I = 1, N - 1
            CALL DLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            WORK( I ) = CS
            WORK( NM1+I ) = SN
   10    CONTINUE
!*
!*        Update singular vectors if desired
!*
         IF( NRU.GT.0 ) &
           CALL DLASR( 'R', 'V', 'F', NRU, N, WORK( 1 ), WORK( N ), U, &
                       LDU )
         IF( NCC.GT.0 ) &
           CALL DLASR( 'L', 'V', 'F', N, NCC, WORK( 1 ), WORK( N ), C, &
                       LDC )
      END IF
!*
!*     Compute singular values to relative accuracy TOL
!*     (By setting TOL to be negative, algorithm will compute
!*!     singular values to absolute accuracy ABS(TOL)*norm(input matrix))
!*
      TOLMUL = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) )
      TOL = TOLMUL*EPS
!*
!*     Compute approximate maximum, minimum singular values
!*
      SMAX = ZERO
      DO 20 I = 1, N
         SMAX = MAX( SMAX, ABS( D( I ) ) )
   20 CONTINUE
      DO 30 I = 1, N - 1
         SMAX = MAX( SMAX, ABS( E( I ) ) )
   30 CONTINUE
      SMINL = ZERO
      IF( TOL.GE.ZERO ) THEN
!*
!*        Relative accuracy desired
!*
         SMINOA = ABS( D( 1 ) )
         IF( SMINOA.EQ.ZERO ) &
           GO TO 50
         MU = SMINOA
         DO 40 I = 2, N
            MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) )
            SMINOA = MIN( SMINOA, MU )
            IF( SMINOA.EQ.ZERO ) &
              GO TO 50
   40    CONTINUE
   50    CONTINUE
         SMINOA = SMINOA / SQRT( DBLE( N ) )
         THRESH = MAX( TOL*SMINOA, MAXITR*(N*(N*UNFL)) )
      ELSE
!*
!*        Absolute accuracy desired
!*
         THRESH = MAX( ABS( TOL )*SMAX, MAXITR*(N*(N*UNFL)) )
      END IF
!*
!*     Prepare for main iteration loop for the singular values
!*     (MAXIT is the maximum number of passes through the inner
!*     loop permitted before nonconvergence signalled.)
!*
      MAXITDIVN = MAXITR*N
      ITERDIVN = 0
      ITER = -1
      OLDLL = -1
      OLDM = -1
!*
!*     M points to last element of unconverged part of matrix
!*
      M = N
!*
!*     Begin main iteration loop
!*
   60 CONTINUE
!*
!*     Check for convergence or exceeding iteration count
!*
      IF( M.LE.1 ) &
         GO TO 160
!*
      IF( ITER.GE.N ) THEN
         ITER = ITER - N
         ITERDIVN = ITERDIVN + 1
         IF( ITERDIVN.GE.MAXITDIVN ) &
           GO TO 200
      END IF
!!*
!*     Find diagonal block of matrix to work on
!*
      IF( TOL.LT.ZERO .AND. ABS( D( M ) ).LE.THRESH ) &
        D( M ) = ZERO
      SMAX = ABS( D( M ) )
      SMIN = SMAX
      DO 70 LLL = 1, M - 1
         LL = M - LLL
         ABSS = ABS( D( LL ) )
         ABSE = ABS( E( LL ) )
         IF( TOL.LT.ZERO .AND. ABSS.LE.THRESH ) &
           D( LL ) = ZERO
         IF( ABSE.LE.THRESH ) &
           GO TO 80
         SMIN = MIN( SMIN, ABSS )
         SMAX = MAX( SMAX, ABSS, ABSE )
   70 CONTINUE
      LL = 0
      GO TO 90
   80 CONTINUE
      E( LL ) = ZERO
!*
!*     Matrix splits since E(LL) = 0
!*
      IF( LL.EQ.M-1 ) THEN
!*
!*        Convergence of bottom singular value, return to top of loop
!*
         M = M - 1
         GO TO 60
      END IF
   90 CONTINUE
      LL = LL + 1
!*
!*     E(LL) through E(M-1) are nonzero, E(LL-1) is zero
!*
      IF( LL.EQ.M-1 ) THEN
!*
!*        2 by 2 block, handle separately
!*
         CALL DLASV2( D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR, &
                      COSR, SINL, COSL )
         D( M-1 ) = SIGMX
         E( M-1 ) = ZERO
         D( M ) = SIGMN
!*
!*        Compute singular vectors, if desired
!*
         IF( NCVT.GT.0 ) &
           CALL DROT( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT, COSR, &
                      SINR )
         IF( NRU.GT.0 ) &
           CALL DROT( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL )
         IF( NCC.GT.0 ) &
           CALL DROT( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL, &
                      SINL )
         M = M - 2
         GO TO 60
      END IF
!*
!*     If working on new submatrix, choose shift direction
!*     (from larger end diagonal element towards smaller)
!*
      IF( LL.GT.OLDM .OR. M.LT.OLDLL ) THEN
         IF( ABS( D( LL ) ).GE.ABS( D( M ) ) ) THEN
!*
!*           Chase bulge from top (big end) to bottom (small end)
!*
            IDIR = 1
         ELSE
!*
!*           Chase bulge from bottom (big end) to top (small end)
!*
            IDIR = 2
         END IF
      END IF
!*
!*     Apply convergence tests
!*
      IF( IDIR.EQ.1 ) THEN
!*
!*        Run convergence test in forward direction
!*        First apply standard test to bottom of matrix
!*
         IF( ABS( E( M-1 ) ).LE.ABS( TOL )*ABS( D( M ) ) .OR. &
            ( TOL.LT.ZERO .AND. ABS( E( M-1 ) ).LE.THRESH ) ) THEN
            E( M-1 ) = ZERO
            GO TO 60
         END IF

         IF( TOL.GE.ZERO ) THEN
!*
!*           If relative accuracy desired,
!*           apply convergence criterion forward
!*
            MU = ABS( D( LL ) )
            SMINL = MU
            DO 100 LLL = LL, M - 1
               IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
                  E( LLL ) = ZERO
                  GO TO 60
               END IF
               MU = ABS( D( LLL+1 ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
  100       CONTINUE
         END IF

      ELSE
!*
!*        Run convergence test in backward direction
!*        First apply standard test to top of matrix
!*
         IF( ABS( E( LL ) ).LE.ABS( TOL )*ABS( D( LL ) ) .OR. &
            ( TOL.LT.ZERO .AND. ABS( E( LL ) ).LE.THRESH ) ) THEN
            E( LL ) = ZERO
            GO TO 60
         END IF

         IF( TOL.GE.ZERO ) THEN
!*
!*           If relative accuracy desired,
!*           apply convergence criterion backward
!*
            MU = ABS( D( M ) )
            SMINL = MU
            DO 110 LLL = M - 1, LL, -1
               IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
                  E( LLL ) = ZERO
                  GO TO 60
               END IF
               MU = ABS( D( LLL ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
  110       CONTINUE
         END IF
      END IF
      OLDLL = LL
      OLDM = M
!*
!*     Compute shift.  First, test if shifting would ruin relative
!*     accuracy, and if so set the shift to zero.
!*
      IF( TOL.GE.ZERO .AND. N*TOL*( SMINL / SMAX ).LE. &
         MAX( EPS, HNDRTH*TOL ) ) THEN
!*
!*        Use a zero shift to avoid loss of relative accuracy
!*
         SHIFT = ZERO
      ELSE
!*
!*        Compute the shift from 2-by-2 block at end of matrix
!*
         IF( IDIR.EQ.1 ) THEN
            SLL = ABS( D( LL ) )
            CALL DLAS2( D( M-1 ), E( M-1 ), D( M ), SHIFT, R )
         ELSE
            SLL = ABS( D( M ) )
            CALL DLAS2( D( LL ), E( LL ), D( LL+1 ), SHIFT, R )
         END IF
!*
!*        Test if shift negligible, and if so set to zero
!*
         IF( SLL.GT.ZERO ) THEN
            IF( ( SHIFT / SLL )**2.LT.EPS ) &
              SHIFT = ZERO
         END IF
      END IF
!*
!*     Increment iteration count
!*
      ITER = ITER + M - LL
!*
!*     If SHIFT = 0, do simplified QR iteration
!*
      IF( SHIFT.EQ.ZERO ) THEN
         IF( IDIR.EQ.1 ) THEN
!*
!*           Chase bulge from top to bottom
!*           Save cosines and sines for later singular vector updates
!*
            CS = ONE
            OLDCS = ONE
            DO 120 I = LL, M - 1
               CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
               IF( I.GT.LL ) &
                  E( I-1 ) = OLDSN*R
               CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN, D( I ) )
               WORK( I-LL+1 ) = CS
               WORK( I-LL+1+NM1 ) = SN
               WORK( I-LL+1+NM12 ) = OLDCS
               WORK( I-LL+1+NM13 ) = OLDSN
  120       CONTINUE
            H = D( M )*CS
            D( M ) = H*OLDCS
            E( M-1 ) = H*OLDSN
!*
!*           Update singular vectors
!*
            IF( NCVT.GT.0 ) &
              CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ), &
                          WORK( N ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 ) &
              CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ), &
                          WORK( NM13+1 ), U( 1, LL ), LDU )
            IF( NCC.GT.0 ) &
              CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ), &
                          WORK( NM13+1 ), C( LL, 1 ), LDC )
!*
!*           Test convergence
!*
            IF( ABS( E( M-1 ) ).LE.THRESH ) &
              E( M-1 ) = ZERO
!*
         ELSE
!*
!*           Chase bulge from bottom to top
!*           Save cosines and sines for later singular vector updates
!*
            CS = ONE
            OLDCS = ONE
            DO 130 I = M, LL + 1, -1
               CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
               IF( I.LT.M ) &
                  E( I ) = OLDSN*R
               CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN, D( I ) )
               WORK( I-LL ) = CS
               WORK( I-LL+NM1 ) = -SN
               WORK( I-LL+NM12 ) = OLDCS
               WORK( I-LL+NM13 ) = -OLDSN
  130       CONTINUE
            H = D( LL )*CS
            D( LL ) = H*OLDCS
            E( LL ) = H*OLDSN
!*
!*           Update singular vectors
!*
            IF( NCVT.GT.0 ) &
              CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ), &
                          WORK( NM13+1 ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 ) &
              CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ), &
                          WORK( N ), U( 1, LL ), LDU )
            IF( NCC.GT.0 ) &
              CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ), &
                          WORK( N ), C( LL, 1 ), LDC )
!*
!*           Test convergence
!*
            IF( ABS( E( LL ) ).LE.THRESH ) &
              E( LL ) = ZERO
         END IF
      ELSE
!*
!*        Use nonzero shift
!*
         IF( IDIR.EQ.1 ) THEN
!*
!*           Chase bulge from top to bottom
!*           Save cosines and sines for later singular vector updates
!*
            F = ( ABS( D( LL ) )-SHIFT )* &
               ( SIGN( ONE, D( LL ) )+SHIFT / D( LL ) )
            G = E( LL )
            DO 140 I = LL, M - 1
               CALL DLARTG( F, G, COSR, SINR, R )
               IF( I.GT.LL ) &
                 E( I-1 ) = R
               F = COSR*D( I ) + SINR*E( I )
               E( I ) = COSR*E( I ) - SINR*D( I )
               G = SINR*D( I+1 )
               D( I+1 ) = COSR*D( I+1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( I ) = R
               F = COSL*E( I ) + SINL*D( I+1 )
               D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
               IF( I.LT.M-1 ) THEN
                  G = SINL*E( I+1 )
                  E( I+1 ) = COSL*E( I+1 )
               END IF
               WORK( I-LL+1 ) = COSR
               WORK( I-LL+1+NM1 ) = SINR
               WORK( I-LL+1+NM12 ) = COSL
               WORK( I-LL+1+NM13 ) = SINL
  140       CONTINUE
            E( M-1 ) = F
!*
!*           Update singular vectors
!*
            IF( NCVT.GT.0 ) &
              CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ), &
                          WORK( N ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 ) &
              CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ), &
                          WORK( NM13+1 ), U( 1, LL ), LDU )
            IF( NCC.GT.0 ) &
              CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ), &
                          WORK( NM13+1 ), C( LL, 1 ), LDC )
!*
!*           Test convergence
!*
            IF( ABS( E( M-1 ) ).LE.THRESH ) &
              E( M-1 ) = ZERO
!*
         ELSE
!*
!*           Chase bulge from bottom to top
!*           Save cosines and sines for later singular vector updates
!*
            F = ( ABS( D( M ) )-SHIFT )*( SIGN( ONE, D( M ) )+SHIFT / &
               D( M ) )
            G = E( M-1 )
            DO 150 I = M, LL + 1, -1
               CALL DLARTG( F, G, COSR, SINR, R )
               IF( I.LT.M ) &
                  E( I ) = R
               F = COSR*D( I ) + SINR*E( I-1 )
               E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
               G = SINR*D( I-1 )
               D( I-1 ) = COSR*D( I-1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( I ) = R
               F = COSL*E( I-1 ) + SINL*D( I-1 )
               D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
               IF( I.GT.LL+1 ) THEN
                  G = SINL*E( I-2 )
                  E( I-2 ) = COSL*E( I-2 )
               END IF
               WORK( I-LL ) = COSR
               WORK( I-LL+NM1 ) = -SINR
               WORK( I-LL+NM12 ) = COSL
               WORK( I-LL+NM13 ) = -SINL
  150       CONTINUE
            E( LL ) = F
!*
!*           Test convergence
!*
            IF( ABS( E( LL ) ).LE.THRESH ) &
              E( LL ) = ZERO
!*
!*           Update singular vectors if desired
!*
            IF( NCVT.GT.0 ) &
              CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ), &
                          WORK( NM13+1 ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 ) &
             CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ), &
                          WORK( N ), U( 1, LL ), LDU )
            IF( NCC.GT.0 ) &
              CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ), &
                          WORK( N ), C( LL, 1 ), LDC )
         END IF
      END IF
!*
!*     QR iteration finished, go back and check convergence
!*
      GO TO 60
!*
!*     All singular values converged, so make them positive
!*
  160 CONTINUE
      DO 170 I = 1, N
         IF( D( I ).LT.ZERO ) THEN
            D( I ) = -D( I )
!*
!*           Change sign of singular vectors, if desired
!*
            IF( NCVT.GT.0 ) &
              CALL DSCAL( NCVT, NEGONE, VT( I, 1 ), LDVT )
         END IF
  170 CONTINUE
!*
!*     Sort the singular values into decreasing order (insertion sort on
!*     singular values, but only one transposition per singular vector)
!*
      DO 190 I = 1, N - 1
!*
!*        Scan for smallest D(I)
!*
         ISUB = 1
         SMIN = D( 1 )
         DO 180 J = 2, N + 1 - I
            IF( D( J ).LE.SMIN ) THEN
               ISUB = J
               SMIN = D( J )
            END IF
  180    CONTINUE
         IF( ISUB.NE.N+1-I ) THEN
!*
!*           Swap singular values and vectors
!*
            D( ISUB ) = D( N+1-I )
            D( N+1-I ) = SMIN
            IF( NCVT.GT.0 ) &
              CALL DSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ), &
                          LDVT )
            IF( NRU.GT.0 ) &
              CALL DSWAP( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 )
            IF( NCC.GT.0 ) &
              CALL DSWAP( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC )
         END IF
  190 CONTINUE
      GO TO 220
!*
!*     Maximum number of iterations exceeded, failure to converge
!*
  200 CONTINUE
      INFO = 0
      DO 210 I = 1, N - 1
         IF( E( I ).NE.ZERO ) &
           INFO = INFO + 1
  210 CONTINUE
  220 CONTINUE
  
END SUBROUTINE

#if
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup OTHERauxiliary
*
*  =====================================================================
#endif

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLARTG( F, G, CS, SN, R ) !GCC$ ATTRIBUTES inline :: DLARTG !GCC$ ATTRIBUTES aligned(32) :: DLARTG
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DLARTG( F, G, CS, SN, R )
    !DIR$ ATTRIBUTES FORCEINLINE :: DLARTG
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLARTG
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLARTG
#endif
      implicit none
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
!*     ..
!*     .. Local Scalars ..
!*     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
!*     ..
!*     .. Save statement ..
*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
!*     ..
!*     .. Data statements ..
!*     DATA               FIRST / .TRUE. /
!*     ..
!*     .. Executable Statements ..
!*
!*     IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / &
                  LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
!*        FIRST = .FALSE.
!*     END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 ) &
              GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 ) &
              GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
     
END SUBROUTINE


!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX ) !GCC$ ATTRIBUTES inline :: DLAS2 !GCC$ ATTRIBUTES aligned(32) :: DLAS2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX )
     !DIR$ ATTRIBUTES FORCEINLINE :: DLAS2
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLAS2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLAS2
#endif
       implicit none
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      DOUBLE PRECISION   F, G, H, SSMAX, SSMIN
!*     ..
!*
!*  ====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION   AS, AT, AU, C, FA, FHMN, FHMX, GA, HA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!*     ..
!*     .. Executable Statements ..
!*
      FA = ABS( F )
      GA = ABS( G )
      HA = ABS( H )
      FHMN = MIN( FA, HA )
      FHMX = MAX( FA, HA )
      IF( FHMN.EQ.ZERO ) THEN
         SSMIN = ZERO
         IF( FHMX.EQ.ZERO ) THEN
            SSMAX = GA
         ELSE
            SSMAX = MAX( FHMX, GA )*SQRT( ONE+ &
                   ( MIN( FHMX, GA ) / MAX( FHMX, GA ) )**2 )
         END IF
      ELSE
         IF( GA.LT.FHMX ) THEN
            AS = ONE + FHMN / FHMX
            AT = ( FHMX-FHMN ) / FHMX
            AU = ( GA / FHMX )**2
            C = TWO / ( SQRT( AS*AS+AU )+SQRT( AT*AT+AU ) )
            SSMIN = FHMN*C
            SSMAX = FHMX / C
         ELSE
            AU = FHMX / GA
            IF( AU.EQ.ZERO ) THEN
!*
!*              Avoid possible harmful underflow if exponent range
!*              asymmetric (true SSMIN may not underflow even if
!*              AU underflows)
!*
               SSMIN = ( FHMN*FHMX ) / GA
               SSMAX = GA
            ELSE
               AS = ONE + FHMN / FHMX
               AT = ( FHMX-FHMN ) / FHMX
               C = ONE / ( SQRT( ONE+( AS*AU )**2 )+ &
                  SQRT( ONE+( AT*AU )**2 ) )
               SSMIN = ( FHMN*C )*AU
               SSMIN = SSMIN + SSMIN
               SSMAX = GA / ( C+C )
            END IF
         END IF
      END IF
  
END SUBROUTINE



!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASQ1( N, D, E, WORK, INFO ) !GCC$ ATTRIBUTES hot :: DLASQ1 !GCC$ ATTRIBUTES aligned(32) :: DLASQ1
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
      
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, IINFO
      DOUBLE PRECISION   EPS, SCALE, SAFMIN, SIGMN, SIGMX, TMP0
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DCOPY

!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!*     ..
!*     .. Executable Statements ..
!*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         !CALL XERBLA( 'DLASQ1', -INFO )
         RETURN
      ELSE IF( N.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         D( 1 ) = ABS( D( 1 ) )
         RETURN
      ELSE IF( N.EQ.2 ) THEN
         CALL DLAS2( D( 1 ), E( 1 ), D( 2 ), SIGMN, SIGMX )
         D( 1 ) = SIGMX
         D( 2 ) = SIGMN
         RETURN
      END IF
!*
!*     Estimate the largest singular value.
!*
      SIGMX = ZERO
      DO 10 I = 1, N - 1
         D( I ) = ABS( D( I ) )
         SIGMX = MAX( SIGMX, ABS( E( I ) ) )
   10 CONTINUE
      D( N ) = ABS( D( N ) )
!*
!*     Early return if SIGMX is zero (matrix is already diagonal).
!*
      IF( SIGMX.EQ.ZERO ) THEN
         CALL DLASRT( 'D', N, D, IINFO )
         RETURN
      END IF
!*
      DO 20 I = 1, N
         SIGMX = MAX( SIGMX, D( I ) )
   20 CONTINUE
!*
!*     Copy D and E into WORK (in the Z format) and scale (squaring the
!*     input data makes scaling by a power of the radix pointless).
!*
      EPS = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SCALE = SQRT( EPS / SAFMIN )
      CALL DCOPY( N, D, 1, WORK( 1 ), 2 )
      CALL DCOPY( N-1, E, 1, WORK( 2 ), 2 )
      CALL DLASCL( 'G', 0, 0, SIGMX, SCALE, 2*N-1, 1, WORK, 2*N-1, &
                   IINFO )
!*
!*     Compute the q's and e's.
!*
      DO 30 I = 1, 2*N - 1
         TMP = WORK(I)
         WORK( I ) = TMP*TMP
   30 CONTINUE
      WORK( 2*N ) = ZERO

      CALL DLASQ2( N, WORK, INFO )

      IF( INFO.EQ.0 ) THEN
         DO 40 I = 1, N
            D( I ) = SQRT( WORK( I ) )
   40    CONTINUE
         CALL DLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO )
      ELSE IF( INFO.EQ.2 ) THEN
!*
!*     Maximum number of iterations exceeded.  Move data from WORK
!*     into D and E so the calling subroutine can try to finish
!*
         DO I = 1, N
            D( I ) = SQRT( WORK( 2*I-1 ) )
            E( I ) = SQRT( WORK( 2*I ) )
         END DO
         CALL DLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO )
         CALL DLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, E, N, IINFO )
      END IF

END SUBROUTINE

#if 0
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2016
*
*> \ingroup OTHERauxiliary
*
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO ) !GCC$ ATTRIBUTES hot :: DLASCL !GCC$ ATTRIBUTES aligned(32) :: DLASCL
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH, DISNAN
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!*     ..
!*     .. External Subroutines ..
    !  EXTERNAL           XERBLA
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
!*
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF

      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO .OR. DISNAN(CFROM) ) THEN
         INFO = -4
      ELSE IF( DISNAN(CTO) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR. &
              ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR. &
                 ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) ) &
                  THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR. &
                 ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR. &
                 ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
        ! CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      IF( N.EQ.0 .OR. M.EQ.0 ) &
        RETURN
!*
!*     Get machine parameters
!*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!*
      CFROMC = CFROM
      CTOC = CTO
!*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      IF( CFROM1.EQ.CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         CTO1 = CTOC
      ELSE
         CTO1 = CTOC / BIGNUM
         IF( CTO1.EQ.CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            MUL = CTOC
            DONE = .TRUE.
            CFROMC = ONE
         ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
            MUL = SMLNUM
            DONE = .FALSE.
            CFROMC = CFROM1
         ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
            MUL = BIGNUM
            DONE = .FALSE.
            CTOC = CTO1
         ELSE
            MUL = CTOC / CFROMC
            DONE = .TRUE.
         END IF
      END IF

      IF( ITYPE.EQ.0 ) THEN
!*
!*        Full matrix
         !*
         !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) COLLAPSE(2) PRIVATE(J,I) SHARED(A,N,M,MUL)
         DO 30 J = 1, N
            !$OMP SIMD LINEAR(I:1) UNROLL PARTIAL(6)
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE

      ELSE IF( ITYPE.EQ.1 ) THEN
!*
!*        Lower triangular matrix
         !*
            !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE)  PRIVATE(J,I) SHARED(A,N,M,MUL)
         DO 50 J = 1, N
             !$OMP SIMD LINEAR(I:1)
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE

      ELSE IF( ITYPE.EQ.2 ) THEN
!*
!*        Upper triangular matrix
         !*
           !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE)  PRIVATE(J,I) SHARED(A,N,M,MUL)
         DO 70 J = 1, N
             !$OMP SIMD 
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE

      ELSE IF( ITYPE.EQ.3 ) THEN
!*
!*        Upper Hessenberg matrix
         !*
           !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE)  PRIVATE(J,I) SHARED(A,N,M,MUL)
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE

      ELSE IF( ITYPE.EQ.4 ) THEN
!*
!*        Lower half of a symmetric band matrix
!*
         K3 = KL + 1
         K4 = N + 1
           !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE)  PRIVATE(J,I) SHARED(A,N,K3,K4,MUL)
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE

      ELSE IF( ITYPE.EQ.5 ) THEN
!*
!*        Upper half of a symmetric band matrix
!*
         K1 = KU + 2
         K3 = KU + 1
           !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE)  PRIVATE(J,I) SHARED(A,K1,K3,MUL,N)
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE

      ELSE IF( ITYPE.EQ.6 ) THEN
!*
!*        Band matrix
!*
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
           !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE)  PRIVATE(J,I) SHARED(A,K1,K2,K3,K4,N,MUL)
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE

      END IF

      IF( .NOT.DONE ) &
        GO TO 10

END SUBROUTINE


!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASQ2( N, Z, INFO ) !GCC$ ATTRIBUTES hot :: DLASQ2 !GCC$ ATTRIBUTES aligned(32) :: DLASQ2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   CBIAS
      PARAMETER          ( CBIAS = 1.50D0 )
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO, FOUR, HUNDRD
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, &
                           TWO = 2.0D0, FOUR = 4.0D0, HUNDRD = 100.0D0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            IEEE
      INTEGER            I0, I1, I4, IINFO, IPN4, ITER, IWHILA, IWHILB, &
                        K, KMIN, N0, N1, NBIG, NDIV, NFAIL, PP, SPLT, &
                        TTYPE
      DOUBLE PRECISION   D, DEE, DEEMIN, DESIG, DMIN, DMIN1, DMIN2, DN, &
                        DN1, DN2, E, EMAX, EMIN, EPS, G, OLDEMN, QMAX, &
                        QMIN, S, SAFMIN, SIGMA, T, TAU, TEMP, TOL, &
                        TOL2, TRACE, ZMAX, TEMPE, TEMPQ
!*     ..
!*     .. External Subroutines ..
   !   EXTERNAL           DLASQ3, DLASRT
!*     ..
!*     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, ILAENV
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments.
!*     (in case DLASQ2 is not called by DLASQ1)
!*
      INFO = 0
      EPS = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      TOL = EPS*HUNDRD
      TOL2 = TOL**2

      IF( N.LT.0 ) THEN
         INFO = -1
         !CALL XERBLA( 'DLASQ2', 1 )
         RETURN
      ELSE IF( N.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
!*
!*        1-by-1 case.
!*
         IF( Z( 1 ).LT.ZERO ) THEN
            INFO = -201
            !CALL XERBLA( 'DLASQ2', 2 )
         END IF
         RETURN
      ELSE IF( N.EQ.2 ) THEN
!*
!*        2-by-2 case.
!*
         IF( Z( 2 ).LT.ZERO .OR. Z( 3 ).LT.ZERO ) THEN
            INFO = -2
            !CALL XERBLA( 'DLASQ2', 2 )
            RETURN
         ELSE IF( Z( 3 ).GT.Z( 1 ) ) THEN
            D = Z( 3 )
            Z( 3 ) = Z( 1 )
            Z( 1 ) = D
         END IF
         Z( 5 ) = Z( 1 ) + Z( 2 ) + Z( 3 )
         IF( Z( 2 ).GT.Z( 3 )*TOL2 ) THEN
            T = HALF*( ( Z( 1 )-Z( 3 ) )+Z( 2 ) )
            S = Z( 3 )*( Z( 2 ) / T )
            IF( S.LE.T ) THEN
               S = Z( 3 )*( Z( 2 ) / ( T*( ONE+SQRT( ONE+S / T ) ) ) )
            ELSE
               S = Z( 3 )*( Z( 2 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
            END IF
            T = Z( 1 ) + ( S+Z( 2 ) )
            Z( 3 ) = Z( 3 )*( Z( 1 ) / T )
            Z( 1 ) = T
         END IF
         Z( 2 ) = Z( 3 )
         Z( 6 ) = Z( 2 ) + Z( 1 )
         RETURN
      END IF
!*
!*     Check for negative data and compute sums of q's and e's.
!*
      Z( 2*N ) = ZERO
      EMIN = Z( 2 )
      QMAX = ZERO
      ZMAX = ZERO
      D = ZERO
      E = ZERO
!*
      DO 10 K = 1, 2*( N-1 ), 2
         IF( Z( K ).LT.ZERO ) THEN
            INFO = -( 200+K )
            !CALL XERBLA( 'DLASQ2', 2 )
            RETURN
         ELSE IF( Z( K+1 ).LT.ZERO ) THEN
            INFO = -( 200+K+1 )
            !CALL XERBLA( 'DLASQ2', 2 )
            RETURN
         END IF
         D = D + Z( K )
         E = E + Z( K+1 )
         QMAX = MAX( QMAX, Z( K ) )
         EMIN = MIN( EMIN, Z( K+1 ) )
         ZMAX = MAX( QMAX, ZMAX, Z( K+1 ) )
   10 CONTINUE
      IF( Z( 2*N-1 ).LT.ZERO ) THEN
         INFO = -( 200+2*N-1 )
         !CALL XERBLA( 'DLASQ2', 2 )
         RETURN
      END IF
      D = D + Z( 2*N-1 )
      QMAX = MAX( QMAX, Z( 2*N-1 ) )
      ZMAX = MAX( QMAX, ZMAX )
!*
!*     Check for diagonality.
!*
      IF( E.EQ.ZERO ) THEN
         DO 20 K = 2, N
            Z( K ) = Z( 2*K-1 )
   20    CONTINUE
         CALL DLASRT( 'D', N, Z, IINFO )
         Z( 2*N-1 ) = D
         RETURN
      END IF
!*
      TRACE = D + E
!*
!*     Check for zero data.
!*
      IF( TRACE.EQ.ZERO ) THEN
         Z( 2*N-1 ) = ZERO
         RETURN
      END IF
!*
!*     Check whether the machine is IEEE conformable.
!*
      IEEE = ILAENV( 10, 'DLASQ2', 'N', 1, 2, 3, 4 ).EQ.1 .AND. &
             ILAENV( 11, 'DLASQ2', 'N', 1, 2, 3, 4 ).EQ.1
!*
!*     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
!*
      DO 30 K = 2*N, 2, -2
         Z( 2*K ) = ZERO
         Z( 2*K-1 ) = Z( K )
         Z( 2*K-2 ) = ZERO
         Z( 2*K-3 ) = Z( K-1 )
   30 CONTINUE
!*
      I0 = 1
      N0 = N
!*
!*     Reverse the qd-array, if warranted.
!*
      IF( CBIAS*Z( 4*I0-3 ).LT.Z( 4*N0-3 ) ) THEN
         IPN4 = 4*( I0+N0 )
         DO 40 I4 = 4*I0, 2*( I0+N0-1 ), 4
            TEMP = Z( I4-3 )
            Z( I4-3 ) = Z( IPN4-I4-3 )
            Z( IPN4-I4-3 ) = TEMP
            TEMP = Z( I4-1 )
            Z( I4-1 ) = Z( IPN4-I4-5 )
            Z( IPN4-I4-5 ) = TEMP
   40    CONTINUE
      END IF
!*
!*     Initial split checking via dqd and Li's test.
!*
      PP = 0
!*
      DO 80 K = 1, 2
!*
         D = Z( 4*N0+PP-3 )
         DO 50 I4 = 4*( N0-1 ) + PP, 4*I0 + PP, -4
            IF( Z( I4-1 ).LE.TOL2*D ) THEN
               Z( I4-1 ) = -ZERO
               D = Z( I4-3 )
            ELSE
               D = Z( I4-3 )*( D / ( D+Z( I4-1 ) ) )
            END IF
   50    CONTINUE
!*
!*        dqd maps Z to ZZ plus Li's test.
!*
         EMIN = Z( 4*I0+PP+1 )
         D = Z( 4*I0+PP-3 )
         DO 60 I4 = 4*I0 + PP, 4*( N0-1 ) + PP, 4
            Z( I4-2*PP-2 ) = D + Z( I4-1 )
            IF( Z( I4-1 ).LE.TOL2*D ) THEN
               Z( I4-1 ) = -ZERO
               Z( I4-2*PP-2 ) = D
               Z( I4-2*PP ) = ZERO
               D = Z( I4+1 )
            ELSE IF( SAFMIN*Z( I4+1 ).LT.Z( I4-2*PP-2 ) .AND. &
                    SAFMIN*Z( I4-2*PP-2 ).LT.Z( I4+1 ) ) THEN
               TEMP = Z( I4+1 ) / Z( I4-2*PP-2 )
               Z( I4-2*PP ) = Z( I4-1 )*TEMP
               D = D*TEMP
            ELSE
               Z( I4-2*PP ) = Z( I4+1 )*( Z( I4-1 ) / Z( I4-2*PP-2 ) )
               D = Z( I4+1 )*( D / Z( I4-2*PP-2 ) )
            END IF
            EMIN = MIN( EMIN, Z( I4-2*PP ) )
   60    CONTINUE
         Z( 4*N0-PP-2 ) = D
!*
!*        Now find qmax.
!*
         QMAX = Z( 4*I0-PP-2 )
         DO 70 I4 = 4*I0 - PP + 2, 4*N0 - PP - 2, 4
            QMAX = MAX( QMAX, Z( I4 ) )
   70    CONTINUE
!*
!*        Prepare for the next iteration on K.
!*
         PP = 1 - PP
   80 CONTINUE
!*
!*     Initialise variables to pass to DLASQ3.
!*
      TTYPE = 0
      DMIN1 = ZERO
      DMIN2 = ZERO
      DN    = ZERO
      DN1   = ZERO
      DN2   = ZERO
      G     = ZERO
      TAU   = ZERO
!*
      ITER = 2
      NFAIL = 0
      NDIV = 2*( N0-I0 )
!*
      DO 160 IWHILA = 1, N + 1
         IF( N0.LT.1 ) &
           GO TO 170
!*
!*        While array unfinished do
!*
!*        E(N0) holds the value of SIGMA when submatrix in I0:N0
!*        splits from the rest of the array, but is negated.
!*
         DESIG = ZERO
         IF( N0.EQ.N ) THEN
            SIGMA = ZERO
         ELSE
            SIGMA = -Z( 4*N0-1 )
         END IF
         IF( SIGMA.LT.ZERO ) THEN
            INFO = 1
            RETURN
         END IF
!*
!*        Find last unreduced submatrix's top index I0, find QMAX and
!*        EMIN. Find Gershgorin-type bound if Q's much greater than E's.
!*
         EMAX = ZERO
         IF( N0.GT.I0 ) THEN
            EMIN = ABS( Z( 4*N0-5 ) )
         ELSE
            EMIN = ZERO
         END IF
         QMIN = Z( 4*N0-3 )
         QMAX = QMIN
         DO 90 I4 = 4*N0, 8, -4
            IF( Z( I4-5 ).LE.ZERO ) &
              GO TO 100
            IF( QMIN.GE.FOUR*EMAX ) THEN
               QMIN = MIN( QMIN, Z( I4-3 ) )
               EMAX = MAX( EMAX, Z( I4-5 ) )
            END IF
            QMAX = MAX( QMAX, Z( I4-7 )+Z( I4-5 ) )
            EMIN = MIN( EMIN, Z( I4-5 ) )
   90    CONTINUE
         I4 = 4

  100    CONTINUE
         I0 = I4 / 4
         PP = 0

         IF( N0-I0.GT.1 ) THEN
            DEE = Z( 4*I0-3 )
            DEEMIN = DEE
            KMIN = I0
            DO 110 I4 = 4*I0+1, 4*N0-3, 4
               DEE = Z( I4 )*( DEE /( DEE+Z( I4-2 ) ) )
               IF( DEE.LE.DEEMIN ) THEN
                  DEEMIN = DEE
                  KMIN = ( I4+3 )/4
               END IF
  110       CONTINUE
            IF( (KMIN-I0)*2.LT.N0-KMIN .AND. &
               DEEMIN.LE.HALF*Z(4*N0-3) ) THEN
               IPN4 = 4*( I0+N0 )
               PP = 2
               DO 120 I4 = 4*I0, 2*( I0+N0-1 ), 4
                  TEMP = Z( I4-3 )
                  Z( I4-3 ) = Z( IPN4-I4-3 )
                  Z( IPN4-I4-3 ) = TEMP
                  TEMP = Z( I4-2 )
                  Z( I4-2 ) = Z( IPN4-I4-2 )
                  Z( IPN4-I4-2 ) = TEMP
                  TEMP = Z( I4-1 )
                  Z( I4-1 ) = Z( IPN4-I4-5 )
                  Z( IPN4-I4-5 ) = TEMP
                  TEMP = Z( I4 )
                  Z( I4 ) = Z( IPN4-I4-4 )
                  Z( IPN4-I4-4 ) = TEMP
  120          CONTINUE
            END IF
         END IF
!*
!*        Put -(initial shift) into DMIN.
!*
         DMIN = -MAX( ZERO, QMIN-TWO*SQRT( QMIN )*SQRT( EMAX ) )
!*
!*        Now I0:N0 is unreduced.
!*        PP = 0 for ping, PP = 1 for pong.
!*        PP = 2 indicates that flipping was applied to the Z array and
!*               and that the tests for deflation upon entry in DLASQ3
!*               should not be performed.
!*
         NBIG = 100*( N0-I0+1 )
         DO 140 IWHILB = 1, NBIG
            IF( I0.GT.N0 ) &
              GO TO 150
!*
!*           While submatrix unfinished take a good dqds step.
!*
            CALL DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, &
                        ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, &
                        DN2, G, TAU )
!*
            PP = 1 - PP
!*
!*           When EMIN is very small check for splits.
!*
            IF( PP.EQ.0 .AND. N0-I0.GE.3 ) THEN
               IF( Z( 4*N0 ).LE.TOL2*QMAX .OR. &
                  Z( 4*N0-1 ).LE.TOL2*SIGMA ) THEN
                  SPLT = I0 - 1
                  QMAX = Z( 4*I0-3 )
                  EMIN = Z( 4*I0-1 )
                  OLDEMN = Z( 4*I0 )
                  DO 130 I4 = 4*I0, 4*( N0-3 ), 4
                     IF( Z( I4 ).LE.TOL2*Z( I4-3 ) .OR. &
                        Z( I4-1 ).LE.TOL2*SIGMA ) THEN
                        Z( I4-1 ) = -SIGMA
                        SPLT = I4 / 4
                        QMAX = ZERO
                        EMIN = Z( I4+3 )
                        OLDEMN = Z( I4+4 )
                     ELSE
                        QMAX = MAX( QMAX, Z( I4+1 ) )
                        EMIN = MIN( EMIN, Z( I4-1 ) )
                        OLDEMN = MIN( OLDEMN, Z( I4 ) )
                     END IF
  130             CONTINUE
                  Z( 4*N0-1 ) = EMIN
                  Z( 4*N0 ) = OLDEMN
                  I0 = SPLT + 1
               END IF
            END IF
!*
  140    CONTINUE
!*
         INFO = 2
!!*
!*        Maximum number of iterations exceeded, restore the shift
!*        SIGMA and place the new d's and e's in a qd array.
!*        This might need to be done for several blocks
!*
         I1 = I0
         N1 = N0
 145     CONTINUE
         TEMPQ = Z( 4*I0-3 )
         Z( 4*I0-3 ) = Z( 4*I0-3 ) + SIGMA
         DO K = I0+1, N0
            TEMPE = Z( 4*K-5 )
            Z( 4*K-5 ) = Z( 4*K-5 ) * (TEMPQ / Z( 4*K-7 ))
            TEMPQ = Z( 4*K-3 )
            Z( 4*K-3 ) = Z( 4*K-3 ) + SIGMA + TEMPE - Z( 4*K-5 )
         END DO
!*
!*        Prepare to do this on the previous block if there is one
!*
         IF( I1.GT.1 ) THEN
            N1 = I1-1
            DO WHILE( ( I1.GE.2 ) .AND. ( Z(4*I1-5).GE.ZERO ) )
               I1 = I1 - 1
            END DO
            SIGMA = -Z(4*N1-1)
            GO TO 145
         END IF

         DO K = 1, N
            Z( 2*K-1 ) = Z( 4*K-3 )
!*
!*        Only the block 1..N0 is unfinished.  The rest of the e's
!*        must be essentially zero, although sometimes other data
!*        has been stored in them.
!*
            IF( K.LT.N0 ) THEN
               Z( 2*K ) = Z( 4*K-1 )
            ELSE
               Z( 2*K ) = 0
            END IF
         END DO
         RETURN
!*
!*        end IWHILB
!*
  150    CONTINUE
!*
  160 CONTINUE
!*
      INFO = 3
      RETURN
!*
!*     end IWHILA
!*
  170 CONTINUE
!*
!*     Move q's to the front.
!*
      DO 180 K = 2, N
         Z( K ) = Z( 4*K-3 )
  180 CONTINUE
!*
!*     Sort and compute sum of eigenvalues.
!*
      CALL DLASRT( 'D', N, Z, IINFO )
!*
      E = ZERO
      DO 190 K = N, 1, -1
         E = E + Z( K )
  190 CONTINUE
!*
!*     Store trace, sum(eigenvalues) and information on performance.
!*
      Z( 2*N+1 ) = TRACE
      Z( 2*N+2 ) = E
      Z( 2*N+3 ) = DBLE( ITER )
      Z( 2*N+4 ) = DBLE( NDIV ) / DBLE( N**2 )
      Z( 2*N+5 ) = HUNDRD*NFAIL / DBLE( ITER )
     
END SUBROUTINE

!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, &
ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, &
DN2, G, TAU ) !GCC$ ATTRIBUTES hot :: DLASQ3 !GCC$ ATTRIBUTES aligned(32) :: DLASQ3
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
      INTEGER            I0, ITER, N0, NDIV, NFAIL, PP, TTYPE
      DOUBLE PRECISION   DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, &
                         QMAX, SIGMA, TAU
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   Z( * )
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   CBIAS
      PARAMETER          ( CBIAS = 1.50D0 )
      DOUBLE PRECISION   ZERO, QURTR, HALF, ONE, TWO, HUNDRD
      PARAMETER          ( ZERO = 0.0D0, QURTR = 0.250D0, HALF = 0.5D0, &
                           ONE = 1.0D0, TWO = 2.0D0, HUNDRD = 100.0D0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            IPN4, J4, N0IN, NN
      DOUBLE PRECISION   EPS, S, T, TEMP, TOL, TOL2
!*     ..
!*     .. External Subroutines ..
      !EXTERNAL           DLASQ4, DLASQ5, DLASQ6
!*     ..
!*     .. External Function ..
      DOUBLE PRECISION   DLAMCH
      LOGICAL            DISNAN
      EXTERNAL           DISNAN, DLAMCH
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!*     ..
!*     .. Executable Statements ..
!*
      N0IN = N0
      EPS = DLAMCH( 'Precision' )
      TOL = EPS*HUNDRD
      TOL2 = TOL**2
!*
!*     Check for deflation.
!*
   10 CONTINUE
!*
      IF( N0.LT.I0 ) &
         RETURN
      IF( N0.EQ.I0 ) &
        GO TO 20
      NN = 4*N0 + PP
      IF( N0.EQ.( I0+1 ) ) &
        GO TO 40
!*
!*     Check whether E(N0-1) is negligible, 1 eigenvalue.
!*
      IF( Z( NN-5 ).GT.TOL2*( SIGMA+Z( NN-3 ) ) .AND. &
         Z( NN-2*PP-4 ).GT.TOL2*Z( NN-7 ) ) &
        GO TO 30

   20 CONTINUE

      Z( 4*N0-3 ) = Z( 4*N0+PP-3 ) + SIGMA
      N0 = N0 - 1
      GO TO 10

!*     Check  whether E(N0-2) is negligible, 2 eigenvalues.
!*
   30 CONTINUE
!*
      IF( Z( NN-9 ).GT.TOL2*SIGMA .AND. &
         Z( NN-2*PP-8 ).GT.TOL2*Z( NN-11 ) ) &
        GO TO 50
!*
   40 CONTINUE
!*
      IF( Z( NN-3 ).GT.Z( NN-7 ) ) THEN
         S = Z( NN-3 )
         Z( NN-3 ) = Z( NN-7 )
         Z( NN-7 ) = S
      END IF
      T = HALF*( ( Z( NN-7 )-Z( NN-3 ) )+Z( NN-5 ) )
      IF( Z( NN-5 ).GT.Z( NN-3 )*TOL2.AND.T.NE.ZERO ) THEN
         S = Z( NN-3 )*( Z( NN-5 ) / T )
         IF( S.LE.T ) THEN
            S = Z( NN-3 )*( Z( NN-5 ) /
               ( T*( ONE+SQRT( ONE+S / T ) ) ) ) &
         ELSE
            S = Z( NN-3 )*( Z( NN-5 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
         END IF
         T = Z( NN-7 ) + ( S+Z( NN-5 ) )
         Z( NN-3 ) = Z( NN-3 )*( Z( NN-7 ) / T )
         Z( NN-7 ) = T
      END IF
      Z( 4*N0-7 ) = Z( NN-7 ) + SIGMA
      Z( 4*N0-3 ) = Z( NN-3 ) + SIGMA
      N0 = N0 - 2
      GO TO 10

   50 CONTINUE
      IF( PP.EQ.2 ) &
        PP = 0
!*
!*     Reverse the qd-array, if warranted.
!*
      IF( DMIN.LE.ZERO .OR. N0.LT.N0IN ) THEN
         IF( CBIAS*Z( 4*I0+PP-3 ).LT.Z( 4*N0+PP-3 ) ) THEN
            IPN4 = 4*( I0+N0 )
            DO 60 J4 = 4*I0, 2*( I0+N0-1 ), 4
               TEMP = Z( J4-3 )
               Z( J4-3 ) = Z( IPN4-J4-3 )
               Z( IPN4-J4-3 ) = TEMP
               TEMP = Z( J4-2 )
               Z( J4-2 ) = Z( IPN4-J4-2 )
               Z( IPN4-J4-2 ) = TEMP
               TEMP = Z( J4-1 )
               Z( J4-1 ) = Z( IPN4-J4-5 )
               Z( IPN4-J4-5 ) = TEMP
               TEMP = Z( J4 )
               Z( J4 ) = Z( IPN4-J4-4 )
               Z( IPN4-J4-4 ) = TEMP
   60       CONTINUE
            IF( N0-I0.LE.4 ) THEN
               Z( 4*N0+PP-1 ) = Z( 4*I0+PP-1 )
               Z( 4*N0-PP ) = Z( 4*I0-PP )
            END IF
            DMIN2 = MIN( DMIN2, Z( 4*N0+PP-1 ) )
            Z( 4*N0+PP-1 ) = MIN( Z( 4*N0+PP-1 ), Z( 4*I0+PP-1 ), &
                                 Z( 4*I0+PP+3 ) )
            Z( 4*N0-PP ) = MIN( Z( 4*N0-PP ), Z( 4*I0-PP ), &
                               Z( 4*I0-PP+4 ) )
            QMAX = MAX( QMAX, Z( 4*I0+PP-3 ), Z( 4*I0+PP+1 ) )
            DMIN = -ZERO
         END IF
      END IF
!*
!*     Choose a shift.
!*
      CALL DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1, &
                  DN2, TAU, TTYPE, G )
!*
!*     Call dqds until DMIN > 0.
!*
   70 CONTINUE
!*
      CALL DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN, &
                  DN1, DN2, IEEE, EPS )
!*
      NDIV = NDIV + ( N0-I0+2 )
      ITER = ITER + 1
!*
!*     Check status.
!*
      IF( DMIN.GE.ZERO .AND. DMIN1.GE.ZERO ) THEN
!*
!*        Success.
!*
         GO TO 90
!*
      ELSE IF( DMIN.LT.ZERO .AND. DMIN1.GT.ZERO .AND. &
              Z( 4*( N0-1 )-PP ).LT.TOL*( SIGMA+DN1 ) .AND. &
              ABS( DN ).LT.TOL*SIGMA ) THEN
!*
!*        Convergence hidden by negative DN.
!*
         Z( 4*( N0-1 )-PP+2 ) = ZERO
         DMIN = ZERO
         GO TO 90
      ELSE IF( DMIN.LT.ZERO ) THEN
!*
!*        TAU too big. Select new TAU and try again.
!*
         NFAIL = NFAIL + 1
         IF( TTYPE.LT.-22 ) THEN
!*
!*           Failed twice. Play it safe.
!*
            TAU = ZERO
         ELSE IF( DMIN1.GT.ZERO ) THEN
!*
!*           Late failure. Gives excellent shift.
!*
            TAU = ( TAU+DMIN )*( ONE-TWO*EPS )
            TTYPE = TTYPE - 11
         ELSE
!*
!*           Early failure. Divide by 4.
!*
            TAU = QURTR*TAU
            TTYPE = TTYPE - 12
         END IF
         GO TO 70
      ELSE IF( DISNAN( DMIN ) ) THEN
!*
!*        NaN.
!*
         IF( TAU.EQ.ZERO ) THEN
            GO TO 80
         ELSE
            TAU = ZERO
            GO TO 70
         END IF
      ELSE
!*
!*        Possible underflow. Play it safe.
!*
         GO TO 80
      END IF
!*
!*     Risk of underflow.
!*
   80 CONTINUE
      CALL DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DN1, DN2 )
      NDIV = NDIV + ( N0-I0+2 )
      ITER = ITER + 1
      TAU = ZERO
!*
   90 CONTINUE
      IF( TAU.LT.SIGMA ) THEN
         DESIG = DESIG + TAU
         T = SIGMA + DESIG
         DESIG = DESIG - ( T-SIGMA )
      ELSE
         T = SIGMA + TAU
         DESIG = SIGMA - ( T-TAU ) + DESIG
      END IF
      SIGMA = T

END SUBROUTINE

!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, &
     DN1, DN2, TAU, TTYPE, G ) !!GCC$ ATTRIBUTES hot :: DLASQ4 !GCC$ ATTRIBUTES aligned(32) :: DLASQ4
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   CNST1, CNST2, CNST3
      PARAMETER          ( CNST1 = 0.5630D0, CNST2 = 1.010D0, &
                         CNST3 = 1.050D0 )
      DOUBLE PRECISION   QURTR, THIRD, HALF, ZERO, ONE, TWO, HUNDRD
      PARAMETER          ( QURTR = 0.250D0, THIRD = 0.3330D0, &
                        HALF = 0.50D0, ZERO = 0.0D0, ONE = 1.0D0, &
                        TWO = 2.0D0, HUNDRD = 100.0D0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I4, NN, NP
      DOUBLE PRECISION   A2, B1, B2, GAM, GAP1, GAP2, S
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!*     ..
!*     .. Executable Statements ..
!*
!*     A negative DMIN forces the shift to take that absolute value
!*     TTYPE records the type of shift.
!*
      IF( DMIN.LE.ZERO ) THEN
         TAU = -DMIN
         TTYPE = -1
         RETURN
      END IF
!*
      NN = 4*N0 + PP
      IF( N0IN.EQ.N0 ) THEN
!*!
!*        No eigenvalues deflated.
!*
         IF( DMIN.EQ.DN .OR. DMIN.EQ.DN1 ) THEN
!*
            B1 = SQRT( Z( NN-3 ) )*SQRT( Z( NN-5 ) )
            B2 = SQRT( Z( NN-7 ) )*SQRT( Z( NN-9 ) )
            A2 = Z( NN-7 ) + Z( NN-5 )
!*
!*           Cases 2 and 3.
!*
            IF( DMIN.EQ.DN .AND. DMIN1.EQ.DN1 ) THEN
               GAP2 = DMIN2 - A2 - DMIN2*QURTR
               IF( GAP2.GT.ZERO .AND. GAP2.GT.B2 ) THEN
                  GAP1 = A2 - DN - ( B2 / GAP2 )*B2
               ELSE
                  GAP1 = A2 - DN - ( B1+B2 )
               END IF
               IF( GAP1.GT.ZERO .AND. GAP1.GT.B1 ) THEN
                  S = MAX( DN-( B1 / GAP1 )*B1, HALF*DMIN )
                  TTYPE = -2
               ELSE
                  S = ZERO
                  IF( DN.GT.B1 ) &
                    S = DN - B1
                  IF( A2.GT.( B1+B2 ) ) &
                    S = MIN( S, A2-( B1+B2 ) )
                  S = MAX( S, THIRD*DMIN )
                  TTYPE = -3
               END IF
            ELSE
!*
!*              Case 4.
!*
               TTYPE = -4
               S = QURTR*DMIN
               IF( DMIN.EQ.DN ) THEN
                  GAM = DN
                  A2 = ZERO
                  IF( Z( NN-5 ) .GT. Z( NN-7 ) ) &
                    RETURN
                  B2 = Z( NN-5 ) / Z( NN-7 )
                  NP = NN - 9
               ELSE
                  NP = NN - 2*PP
                  GAM = DN1
                  IF( Z( NP-4 ) .GT. Z( NP-2 ) ) &
                    RETURN
                  A2 = Z( NP-4 ) / Z( NP-2 )
                  IF( Z( NN-9 ) .GT. Z( NN-11 ) ) &
                    RETURN
                  B2 = Z( NN-9 ) / Z( NN-11 )
                  NP = NN - 13
               END IF
!*
!*              Approximate contribution to norm squared from I < NN-1.
!*
               A2 = A2 + B2
               DO 10 I4 = NP, 4*I0 - 1 + PP, -4
                  IF( B2.EQ.ZERO ) &
                    GO TO 20
                  B1 = B2
                  IF( Z( I4 ) .GT. Z( I4-2 ) ) &
                    RETURN
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                  A2 = A2 + B2
                  IF( HUNDRD*MAX( B2, B1 ).LT.A2 .OR. CNST1.LT.A2 ) &
                    GO TO 20
   10          CONTINUE
   20          CONTINUE
               A2 = CNST3*A2
!*
!*              Rayleigh quotient residual bound.
!*
               IF( A2.LT.CNST1 ) &
                 S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
            END IF
         ELSE IF( DMIN.EQ.DN2 ) THEN
!*
!*           Case 5.
!*
            TTYPE = -5
            S = QURTR*DMIN
!*
!*           Compute contribution to norm squared from I > NN-2.
!*
            NP = NN - 2*PP
            B1 = Z( NP-2 )
            B2 = Z( NP-6 )
            GAM = DN2
            IF( Z( NP-8 ).GT.B2 .OR. Z( NP-4 ).GT.B1 ) &
              RETURN
            A2 = ( Z( NP-8 ) / B2 )*( ONE+Z( NP-4 ) / B1 )
!*
!*           Approximate contribution to norm squared from I < NN-2.
!*
            IF( N0-I0.GT.2 ) THEN
               B2 = Z( NN-13 ) / Z( NN-15 )
               A2 = A2 + B2
               DO 30 I4 = NN - 17, 4*I0 - 1 + PP, -4
                  IF( B2.EQ.ZERO ) &
                    GO TO 40
                  B1 = B2
                  IF( Z( I4 ) .GT. Z( I4-2 ) ) &
                    RETURN
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                  A2 = A2 + B2
                  IF( HUNDRD*MAX( B2, B1 ).LT.A2 .OR. CNST1.LT.A2 ) &
                    GO TO 40
   30          CONTINUE
   40          CONTINUE
               A2 = CNST3*A2
            END IF

            IF( A2.LT.CNST1 ) &
              S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
         ELSE
!*
!*           Case 6, no information to guide us.
!*
            IF( TTYPE.EQ.-6 ) THEN
               G = G + THIRD*( ONE-G )
            ELSE IF( TTYPE.EQ.-18 ) THEN
               G = QURTR*THIRD
            ELSE
               G = QURTR
            END IF
            S = G*DMIN
            TTYPE = -6
         END IF
!*
      ELSE IF( N0IN.EQ.( N0+1 ) ) THEN
!*
!*        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
!*
         IF( DMIN1.EQ.DN1 .AND. DMIN2.EQ.DN2 ) THEN
!*
!*           Cases 7 and 8.
!*
            TTYPE = -7
            S = THIRD*DMIN1
            IF( Z( NN-5 ).GT.Z( NN-7 ) ) &
              RETURN
            B1 = Z( NN-5 ) / Z( NN-7 )
            B2 = B1
            IF( B2.EQ.ZERO ) &
              GO TO 60
            DO 50 I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               A2 = B1
               IF( Z( I4 ).GT.Z( I4-2 ) ) &
                 RETURN
               B1 = B1*( Z( I4 ) / Z( I4-2 ) )
               B2 = B2 + B1
               IF( HUNDRD*MAX( B1, A2 ).LT.B2 ) &
                 GO TO 60
   50       CONTINUE
   60       CONTINUE
            B2 = SQRT( CNST3*B2 )
            A2 = DMIN1 / ( ONE+B2**2 )
            GAP2 = HALF*DMIN2 - A2
            IF( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) THEN
               S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
            ELSE
               S = MAX( S, A2*( ONE-CNST2*B2 ) )
               TTYPE = -8
            END IF
         ELSE
!*
!*           Case 9.
!*
            S = QURTR*DMIN1
            IF( DMIN1.EQ.DN1 ) &
             S = HALF*DMIN1
            TTYPE = -9
         END IF
!*
      ELSE IF( N0IN.EQ.( N0+2 ) ) THEN
!*
!*        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
!*
!*        Cases 10 and 11.
!*
         IF( DMIN2.EQ.DN2 .AND. TWO*Z( NN-5 ).LT.Z( NN-7 ) ) THEN
            TTYPE = -10
            S = THIRD*DMIN2
            IF( Z( NN-5 ).GT.Z( NN-7 ) ) &
              RETURN
            B1 = Z( NN-5 ) / Z( NN-7 )
            B2 = B1
            IF( B2.EQ.ZERO ) &
              GO TO 80
            DO 70 I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               IF( Z( I4 ).GT.Z( I4-2 ) ) &
                 RETURN
               B1 = B1*( Z( I4 ) / Z( I4-2 ) )
               B2 = B2 + B1
               IF( HUNDRD*B1.LT.B2 ) &
                 GO TO 80
   70       CONTINUE
   80       CONTINUE
            B2 = SQRT( CNST3*B2 )
            A2 = DMIN2 / ( ONE+B2**2 )
            GAP2 = Z( NN-7 ) + Z( NN-9 ) - &
                  SQRT( Z( NN-11 ) )*SQRT( Z( NN-9 ) ) - A2
            IF( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) THEN
               S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
            ELSE
               S = MAX( S, A2*( ONE-CNST2*B2 ) )
            END IF
         ELSE
            S = QURTR*DMIN2
            TTYPE = -11
         END IF
      ELSE IF( N0IN.GT.( N0+2 ) ) THEN
!*
!*        Case 12, more than two eigenvalues deflated. No information.
!*
         S = ZERO
         TTYPE = -12
      END IF

      TAU = S
   
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, &
     DN, DNM1, DNM2, IEEE, EPS ) !GCC$ ATTRIBUTES hot :: DLASQ5 !GCC$ ATTRIBUTES aligned(32) :: DLASQ5
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameter ..
      DOUBLE PRECISION   ZERO, HALF
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            J4, J4P2
      DOUBLE PRECISION   D, EMIN, TEMP, DTHRESH
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MIN
!*     ..
*!     .. Executable Statements ..
!*
      IF( ( N0-I0-1 ).LE.0 ) &
        RETURN
!*
      DTHRESH = EPS*(SIGMA+TAU)
      IF( TAU.LT.DTHRESH*HALF ) TAU = ZERO
      IF( TAU.NE.ZERO ) THEN
      J4 = 4*I0 + PP - 3
      EMIN = Z( J4+4 )
      D = Z( J4 ) - TAU
      DMIN = D
      DMIN1 = -Z( J4 )

      IF( IEEE ) THEN
!*
!*        Code for IEEE arithmetic.
!*
         IF( PP.EQ.0 ) THEN
            DO 10 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-2 ) = D + Z( J4-1 )
               TEMP = Z( J4+1 ) / Z( J4-2 )
               D = D*TEMP - TAU
               DMIN = MIN( DMIN, D )
               Z( J4 ) = Z( J4-1 )*TEMP
               EMIN = MIN( Z( J4 ), EMIN )
   10       CONTINUE
         ELSE
            DO 20 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-3 ) = D + Z( J4 )
               TEMP = Z( J4+2 ) / Z( J4-3 )
               D = D*TEMP - TAU
               DMIN = MIN( DMIN, D )
               Z( J4-1 ) = Z( J4 )*TEMP
               EMIN = MIN( Z( J4-1 ), EMIN )
   20       CONTINUE
         END IF
!*
!*        Unroll last two steps.
!*
         DNM2 = D
         DMIN2 = DMIN
         J4 = 4*( N0-2 ) - PP
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM2 + Z( J4P2 )
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
         DMIN = MIN( DMIN, DNM1 )
!*
         DMIN1 = DMIN
         J4 = J4 + 4
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM1 + Z( J4P2 )
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
         DMIN = MIN( DMIN, DN )
!*
      ELSE
!*
!*        Code for non IEEE arithmetic.
!*
         IF( PP.EQ.0 ) THEN
            DO 30 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-2 ) = D + Z( J4-1 )
               IF( D.LT.ZERO ) THEN
                  RETURN
               ELSE
                  Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
                  D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU
               END IF
               DMIN = MIN( DMIN, D )
               EMIN = MIN( EMIN, Z( J4 ) )
   30       CONTINUE
         ELSE
            DO 40 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-3 ) = D + Z( J4 )
               IF( D.LT.ZERO ) THEN
                  RETURN
               ELSE
                  Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
                  D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU
               END IF
               DMIN = MIN( DMIN, D )
               EMIN = MIN( EMIN, Z( J4-1 ) )
   40       CONTINUE
         END IF
!*
!*        Unroll last two steps.
!*
         DNM2 = D
         DMIN2 = DMIN
         J4 = 4*( N0-2 ) - PP
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM2 + Z( J4P2 )
         IF( DNM2.LT.ZERO ) THEN
            RETURN
         ELSE
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
         END IF
         DMIN = MIN( DMIN, DNM1 )
!*
         DMIN1 = DMIN
         J4 = J4 + 4
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM1 + Z( J4P2 )
         IF( DNM1.LT.ZERO ) THEN
            RETURN
         ELSE
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
         END IF
         DMIN = MIN( DMIN, DN )
!*
      END IF
      ELSE
!*     This is the version that sets d's to zero if they are small enough
         J4 = 4*I0 + PP - 3
         EMIN = Z( J4+4 )
         D = Z( J4 ) - TAU
         DMIN = D
         DMIN1 = -Z( J4 )
         IF( IEEE ) THEN
!*
!*     Code for IEEE arithmetic.
!*
            IF( PP.EQ.0 ) THEN
               DO 50 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-2 ) = D + Z( J4-1 )
                  TEMP = Z( J4+1 ) / Z( J4-2 )
                  D = D*TEMP - TAU
                  IF( D.LT.DTHRESH ) D = ZERO
                  DMIN = MIN( DMIN, D )
                  Z( J4 ) = Z( J4-1 )*TEMP
                  EMIN = MIN( Z( J4 ), EMIN )
 50            CONTINUE
            ELSE
               DO 60 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-3 ) = D + Z( J4 )
                  TEMP = Z( J4+2 ) / Z( J4-3 )
                  D = D*TEMP - TAU
                  IF( D.LT.DTHRESH ) D = ZERO
                  DMIN = MIN( DMIN, D )
                  Z( J4-1 ) = Z( J4 )*TEMP
                  EMIN = MIN( Z( J4-1 ), EMIN )
 60            CONTINUE
            END IF
!*
!*     Unroll last two steps.
!*
            DNM2 = D
            DMIN2 = DMIN
            J4 = 4*( N0-2 ) - PP
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM2 + Z( J4P2 )
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
            DMIN = MIN( DMIN, DNM1 )
!*
            DMIN1 = DMIN
            J4 = J4 + 4
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM1 + Z( J4P2 )
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
            DMIN = MIN( DMIN, DN )
!*
         ELSE
!*
!*     Code for non IEEE arithmetic.
!*
            IF( PP.EQ.0 ) THEN
               DO 70 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-2 ) = D + Z( J4-1 )
                  IF( D.LT.ZERO ) THEN
                     RETURN
                  ELSE
                     Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
                     D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU
                  END IF
                  IF( D.LT.DTHRESH) D = ZERO
                  DMIN = MIN( DMIN, D )
                  EMIN = MIN( EMIN, Z( J4 ) )
 70            CONTINUE
            ELSE
               DO 80 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-3 ) = D + Z( J4 )
                  IF( D.LT.ZERO ) THEN
                     RETURN
                  ELSE
                     Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
                     D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU
                  END IF
                  IF( D.LT.DTHRESH) D = ZERO
                  DMIN = MIN( DMIN, D )
                  EMIN = MIN( EMIN, Z( J4-1 ) )
 80            CONTINUE
            END IF
!*
!*     Unroll last two steps.
!*
            DNM2 = D
            DMIN2 = DMIN
            J4 = 4*( N0-2 ) - PP
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM2 + Z( J4P2 )
            IF( DNM2.LT.ZERO ) THEN
               RETURN
            ELSE
               Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
               DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
            END IF
            DMIN = MIN( DMIN, DNM1 )
!*
            DMIN1 = DMIN
            J4 = J4 + 4
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM1 + Z( J4P2 )
            IF( DNM1.LT.ZERO ) THEN
               RETURN
            ELSE
               Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
               DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
            END IF
            DMIN = MIN( DMIN, DN )
!*
         END IF
      END IF
!*
      Z( J4+2 ) = DN
      Z( 4*N0-PP ) = EMIN
   
END SUBROUTINE

    
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      SUBROUTINE DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, &
           DNM1, DNM2 ) !GCC$ ATTRIBUTES hot :: DLASQ6 !GCC$ ATTRIBUTES aligned(32) :: DLASQ6
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameter ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            J4, J4P2
      DOUBLE PRECISION   D, EMIN, SAFMIN, TEMP
!*     ..
!*     .. External Function ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MIN
!*     ..
!*     .. Executable Statements ..
!*
      IF( ( N0-I0-1 ).LE.0 ) RETURN
      

      SAFMIN = DLAMCH( 'Safe minimum' )
      J4 = 4*I0 + PP - 3
      EMIN = Z( J4+4 )
      D = Z( J4 )
      DMIN = D

      IF( PP.EQ.0 ) THEN
         DO 10 J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-2 ) = D + Z( J4-1 )
            IF( Z( J4-2 ).EQ.ZERO ) THEN
               Z( J4 ) = ZERO
               D = Z( J4+1 )
               DMIN = D
               EMIN = ZERO
            ELSE IF( SAFMIN*Z( J4+1 ).LT.Z( J4-2 ) .AND. &
                    SAFMIN*Z( J4-2 ).LT.Z( J4+1 ) ) THEN
               TEMP = Z( J4+1 ) / Z( J4-2 )
               Z( J4 ) = Z( J4-1 )*TEMP
               D = D*TEMP
            ELSE
               Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
               D = Z( J4+1 )*( D / Z( J4-2 ) )
            END IF
            DMIN = MIN( DMIN, D )
            EMIN = MIN( EMIN, Z( J4 ) )
   10    CONTINUE
      ELSE
         DO 20 J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-3 ) = D + Z( J4 )
            IF( Z( J4-3 ).EQ.ZERO ) THEN
               Z( J4-1 ) = ZERO
               D = Z( J4+2 )
               DMIN = D
               EMIN = ZERO
            ELSE IF( SAFMIN*Z( J4+2 ).LT.Z( J4-3 ) .AND. &
                    SAFMIN*Z( J4-3 ).LT.Z( J4+2 ) ) THEN
               TEMP = Z( J4+2 ) / Z( J4-3 )
               Z( J4-1 ) = Z( J4 )*TEMP
               D = D*TEMP
            ELSE
               Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
               D = Z( J4+2 )*( D / Z( J4-3 ) )
            END IF
            DMIN = MIN( DMIN, D )
            EMIN = MIN( EMIN, Z( J4-1 ) )
   20    CONTINUE
      END IF
!*
!*     Unroll last two steps.
!*
      DNM2 = D
      DMIN2 = DMIN
      J4 = 4*( N0-2 ) - PP
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM2 + Z( J4P2 )
      IF( Z( J4-2 ).EQ.ZERO ) THEN
         Z( J4 ) = ZERO
         DNM1 = Z( J4P2+2 )
         DMIN = DNM1
         EMIN = ZERO
      ELSE IF( SAFMIN*Z( J4P2+2 ).LT.Z( J4-2 ) .AND. &
              SAFMIN*Z( J4-2 ).LT.Z( J4P2+2 ) ) THEN
         TEMP = Z( J4P2+2 ) / Z( J4-2 )
         Z( J4 ) = Z( J4P2 )*TEMP
         DNM1 = DNM2*TEMP
      ELSE
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) )
      END IF
      DMIN = MIN( DMIN, DNM1 )

      DMIN1 = DMIN
      J4 = J4 + 4
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM1 + Z( J4P2 )
      IF( Z( J4-2 ).EQ.ZERO ) THEN
         Z( J4 ) = ZERO
         DN = Z( J4P2+2 )
         DMIN = DN
         EMIN = ZERO
      ELSE IF( SAFMIN*Z( J4P2+2 ).LT.Z( J4-2 ) .AND. &
              SAFMIN*Z( J4-2 ).LT.Z( J4P2+2 ) ) THEN
         TEMP = Z( J4P2+2 ) / Z( J4-2 )
         Z( J4 ) = Z( J4P2 )*TEMP
         DN = DNM1*TEMP
      ELSE
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) )
      END IF
      DMIN = MIN( DMIN, DN )

      Z( J4+2 ) = DN
      Z( 4*N0-PP ) = EMIN
 
END SUBROUTINE

    
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASRT( ID, N, D, INFO ) !GCC$ ATTRIBUTES hot :: DLASRT !GCC$ ATTRIBUTES aligned(32) :: DLASRT
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
!*     ..
!*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           XERBLA
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DLASRT', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      !IF( N.LE.1 ) 
    
!*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
!*
!*        Do Insertion sort on D( START:ENDD )
!*
         IF( DIR.EQ.0 ) THEN
!*
!*           Sort into decreasing order
!*
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
!*
         ELSE
!*
!*           Sort into increasing order
!*
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
!*
         END IF
!*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
!*
!*        Partition D( START:ENDD ) and stack parts, largest one first
!*
!*        Choose partition entry as median of 3
!*
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
!*
         IF( DIR.EQ.0 ) THEN
!*
!*           Sort into decreasing order
!*
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX ) &
               GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX ) &
              GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
!*
!*           Sort into increasing order
!*
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX ) &
              GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX ) &
              GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 ) &
         GO TO 10
    
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA ) !GCC$ ATTRIBUTES hot :: DLASR !GCC$ ATTRIBUTES aligned(32) :: DLASR
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. External Subroutines ..
     ! EXTERNAL           XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters
!*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT, &
              'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) ) &
               THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
   !   IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
    ! $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
!*
!*        Form  P * A
!*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!*
!*        Form A * P**T
!*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF

END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL ) !GCC$ ATTRIBUTES inline :: DLASV2 !GCC$ ATTRIBUTES aligned(32) :: DLASV2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )
    !DIR$ ATTRIBUTES FORCEINLINE :: DLASV2
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASV2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASV2
#endif
      implicit none
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      DOUBLE PRECISION   CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN
!*     ..
!*
!* =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   FOUR
      PARAMETER          ( FOUR = 4.0D0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            GASMAL, SWAP
      INTEGER            PMAX
      DOUBLE PRECISION   A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M, &
                         MM, R, S, SLT, SRT, T, TEMP, TSIGN, TT
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!*     ..
!*     .. Executable Statements ..
!*
      FT = F
      FA = ABS( FT )
      HT = H
      HA = ABS( H )
!*
!*     PMAX points to the maximum absolute element of matrix
!*       PMAX = 1 if F largest in absolute values
!*       PMAX = 2 if G largest in absolute values
!*       PMAX = 3 if H largest in absolute values
!*
      PMAX = 1
      SWAP = ( HA.GT.FA )
      IF( SWAP ) THEN
         PMAX = 3
         TEMP = FT
         FT = HT
         HT = TEMP
         TEMP = FA
         FA = HA
         HA = TEMP
!*
!*        Now FA .ge. HA
!*
      END IF
      GT = G
      GA = ABS( GT )
      IF( GA.EQ.ZERO ) THEN
!*
!*        Diagonal matrix
!*
         SSMIN = HA
         SSMAX = FA
         CLT = ONE
         CRT = ONE
         SLT = ZERO
         SRT = ZERO
      ELSE
         GASMAL = .TRUE.
         IF( GA.GT.FA ) THEN
            PMAX = 2
            IF( ( FA / GA ).LT.DLAMCH( 'EPS' ) ) THEN
!*
!*              Case of very large GA
!*
               GASMAL = .FALSE.
               SSMAX = GA
               IF( HA.GT.ONE ) THEN
                  SSMIN = FA / ( GA / HA )
               ELSE
                  SSMIN = ( FA / GA )*HA
               END IF
               CLT = ONE
               SLT = HT / GT
               SRT = ONE
               CRT = FT / GT
            END IF
         END IF
         IF( GASMAL ) THEN
!!*
!*           Normal case
!*
            D = FA - HA
            IF( D.EQ.FA ) THEN
!*
!*              Copes with infinite F or H
!*
               L = ONE
            ELSE
               L = D / FA
            END IF
!*
!*           Note that 0 .le. L .le. 1
!*
            M = GT / FT
!*
!*           Note that abs(M) .le. 1/macheps
!*
            T = TWO - L
!*
!*           Note that T .ge. 1
!*
            MM = M*M
            TT = T*T
            S = SQRT( TT+MM )
!*
!*           Note that 1 .le. S .le. 1 + 1/macheps
!*
            IF( L.EQ.ZERO ) THEN
               R = ABS( M )
            ELSE
               R = SQRT( L*L+MM )
            END IF
!*
!*           Note that 0 .le. R .le. 1 + 1/macheps
!*
            A = HALF*( S+R )
!*
!*           Note that 1 .le. A .le. 1 + abs(M)
!*
            SSMIN = HA / A
            SSMAX = FA*A
            IF( MM.EQ.ZERO ) THEN
!*
!*              Note that M is very tiny
!*
               IF( L.EQ.ZERO ) THEN
                  T = SIGN( TWO, FT )*SIGN( ONE, GT )
               ELSE
                  T = GT / SIGN( D, FT ) + M / T
               END IF
            ELSE
               T = ( M / ( S+T )+M / ( R+L ) )*( ONE+A )
            END IF
            L = SQRT( T*T+FOUR )
            CRT = TWO / L
            SRT = T / L
            CLT = ( CRT+SRT*M ) / A
            SLT = ( HT / FT )*SRT / A
         END IF
      END IF
      IF( SWAP ) THEN
         CSL = SRT
         SNL = CRT
         CSR = SLT
         SNR = CLT
      ELSE
         CSL = CLT
         SNL = SLT
         CSR = CRT
         SNR = SRT
      END IF
!*
!*     Correct signs of SSMAX and SSMIN
!*
      IF( PMAX.EQ.1 )&
        TSIGN = SIGN( ONE, CSR )*SIGN( ONE, CSL )*SIGN( ONE, F )
      IF( PMAX.EQ.2 ) &
        TSIGN = SIGN( ONE, SNR )*SIGN( ONE, CSL )*SIGN( ONE, G )
      IF( PMAX.EQ.3 ) &
        TSIGN = SIGN( ONE, SNR )*SIGN( ONE, SNL )*SIGN( ONE, H )
      SSMAX = SIGN( SSMAX, TSIGN )
      SSMIN = SIGN( SSMIN, TSIGN*SIGN( ONE, F )*SIGN( ONE, H ) )
  
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) !GCC$ ATTRIBUTES hot :: DGELQF !GCC$ ATTRIBUTES aligned(32) :: DGELQF
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, &
                        NBMIN, NX
!*     ..
!*     .. External Subroutines ..
     
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      NB = ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
      LWKOPT = M*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DGELQF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!*
      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!*
!*        Determine when to cross over from blocked to unblocked code.
!*
         NX = MAX( 0, ILAENV( 3, 'DGELQF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
!*
!*           Determine if workspace is large enough for blocked code.
!*
            LDWORK = M
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!*
!*              Not enough workspace to use optimal NB:  reduce NB and
!*              determine the minimum value of NB.
!*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DGELQF', ' ', M, N, -1, &
                      -1 ) )
            END IF
         END IF
      END IF

      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!*
!*        Use blocked code initially
!*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
!*
!*           Compute the LQ factorization of the current block
!*           A(i:i+ib-1,i:n)
!*
            CALL DGELQ2( IB, N-I+1, A( I, I ), LDA, TAU( I ), WORK, &
                         IINFO )
            IF( I+IB.LE.M ) THEN
!*
!*              Form the triangular factor of the block reflector
!*              H = H(i) H(i+1) . . . H(i+ib-1)
!*
               CALL DLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ), &
                           LDA, TAU( I ), WORK, LDWORK )
!*
!*              Apply H to A(i+ib:m,i:n) from the right
!*
               CALL DLARFB( 'Right', 'No transpose', 'Forward', &
                           'Rowwise', M-I-IB+1, N-I+1, IB, A( I, I ), &
                           LDA, WORK, LDWORK, A( I+IB, I ), LDA, &
                           WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
!*
!*     Use unblocked code to factor the last or only block.
!*
      IF( I.LE.K ) &
        CALL DGELQ2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, &
                     IINFO )

      WORK( 1 ) = IWS
   
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) !GCC$ ATTRIBUTES hot :: DGEQRF !GCC$ ATTRIBUTES aligned(32) :: DGEQRF
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, &
                        NBMIN, NX
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!*     ..
!*     .. Executable Statements ..
!*
!!*     Test the input arguments
!*
      INFO = 0
      NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DGEQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!*
!*        Determine when to cross over from blocked to unblocked code.
!*
         NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
!*
!*           Determine if workspace is large enough for blocked code.
!*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!*
!*              Not enough workspace to use optimal NB:  reduce NB and
!*              determine the minimum value of NB.
!*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DGEQRF', ' ', M, N, -1, &
                      -1 ) )
            END IF
         END IF
      END IF

      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!*
!*        Use blocked code initially
!*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
!*
!*           Compute the QR factorization of the current block
!*           A(i:m,i:i+ib-1)
!*
            CALL DGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, &
                        IINFO )
            IF( I+IB.LE.N ) THEN
!*
!*              Form the triangular factor of the block reflector
!*              H = H(i) H(i+1) . . . H(i+ib-1)
!*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                           A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!*
!*              Apply H**T to A(i:m,i+ib:n) from the left
!*
               CALL DLARFB( 'Left', 'Transpose', 'Forward', &
                           'Columnwise', M-I+1, N-I-IB+1, IB, &
                           A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                           LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
!*
!*     Use unblocked code to factor the last or only block.
!*
      IF( I.LE.K ) &
        CALL DGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, &
                     IINFO )

      WORK( 1 ) = IWS
   
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO ) !GCC$ ATTRIBUTES hot :: DGEQR2 !GCC$ ATTRIBUTES aligned(32) :: DGEQR2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
!*     ..
!*     .. External Subroutines ..
      !EXTERNAL           DLARF, DLARFG, XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DGEQR2', -INFO )
         RETURN
      END IF
!*
      K = MIN( M, N )
!*
      DO 10 I = 1, K
!*
!*        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!*
         CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, &
                     TAU( I ) )
         IF( I.LT.N ) THEN
!*
!*           Apply H(i) to A(i:m,i+1:n) from the left
!*
            AII = A( I, I )
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                       A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
    
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) !GCC$ ATTRIBUTES hot :: DORGBR !GCC$ ATTRIBUTES aligned(32) :: DORGBR
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, WANTQ
      INTEGER            I, IINFO, J, LWKOPT, MN
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           DORGLQ, DORGQR
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      WANTQ = LSAME( VECT, 'Q' )
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 .OR. ( WANTQ .AND. ( N.GT.M .OR. N.LT.MIN( M, &
              K ) ) ) .OR. ( .NOT.WANTQ .AND. ( M.GT.N .OR. M.LT. &
              MIN( N, K ) ) ) ) THEN
         INFO = -3
      ELSE IF( K.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LWORK.LT.MAX( 1, MN ) .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF

      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = 1
         IF( WANTQ ) THEN
            IF( M.GE.K ) THEN
               CALL DORGQR( M, N, K, A, LDA, TAU, WORK, -1, IINFO )
            ELSE
               IF( M.GT.1 ) THEN
                  CALL DORGQR( M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK, &
                              -1, IINFO )
               END IF
            END IF
         ELSE
            IF( K.LT.N ) THEN
               CALL DORGLQ( M, N, K, A, LDA, TAU, WORK, -1, IINFO )
            ELSE
               IF( N.GT.1 ) THEN
                  CALL DORGLQ( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, &
                              -1, IINFO )
               END IF
            END IF
         END IF
         LWKOPT = WORK( 1 )
         LWKOPT = MAX (LWKOPT, MN)
      END IF

      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DORGBR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         WORK( 1 ) = LWKOPT
         RETURN
      END IF
!*
!*     Quick return if possible
!*
     ! IF( M.EQ.0 .OR. N.EQ.0 ) THEN
     !    WORK( 1 ) = 1
     !    RETURN
     ! END IF
!*
      IF( WANTQ ) THEN
!*
!*        Form Q, determined by a call to DGEBRD to reduce an m-by-k
!*        matrix
!*
         IF( M.GE.K ) THEN
!*
!*           If m >= k, assume m >= n >= k
!*
            CALL DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
!*
         ELSE
!*
!*           If m < k, assume m = n
!*
!*           Shift the vectors which define the elementary reflectors one
!*           column to the right, and set the first row and column of Q
!*           to those of the unit matrix
!*
            DO 20 J = M, 2, -1
               A( 1, J ) = ZERO
               !$OMP SIMD ALIGNED(A:64)
               DO 10 I = J + 1, M
                  A( I, J ) = A( I, J-1 )
   10          CONTINUE
   20       CONTINUE
            A( 1, 1 ) = ONE
             !$OMP SIMD ALIGNED(A:64)  UNROLL PARTIAL(8)   
            DO 30 I = 2, M
               A( I, 1 ) = ZERO
   30       CONTINUE
            IF( M.GT.1 ) THEN
!*
!*              Form Q(2:m,2:m)
!*
               CALL DORGQR( M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK, &
                           LWORK, IINFO )
            END IF
         END IF
      ELSE
!*
!*        Form P**T, determined by a call to DGEBRD to reduce a k-by-n
!*        matrix
!*
         IF( K.LT.N ) THEN
!*
!*           If k < n, assume k <= m <= n
!*
            CALL DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
!*
         ELSE
!*
!*           If k >= n, assume m = n
!*
!*           Shift the vectors which define the elementary reflectors one
!*           row downward, and set the first row and column of P**T to
!*           those of the unit matrix
!*
            A( 1, 1 ) = ONE
             !$OMP SIMD ALIGNED(A:64) UNROLL PARTIAL(8) 
            DO 40 I = 2, N
               A( I, 1 ) = ZERO
40          CONTINUE
              
            DO 60 J = 2, N
                !$OMP SIMD ALIGNED(A:64) 
               DO 50 I = J - 1, 2, -1
                  A( I, J ) = A( I-1, J )
   50          CONTINUE
               A( 1, J ) = ZERO
   60       CONTINUE
            IF( N.GT.1 ) THEN
!*
!*              Form P**T(2:n,2:n)
!*
               CALL DORGLQ( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, &
                           LWORK, IINFO )
            END IF
         END IF
      END IF
      WORK( 1 ) = LWKOPT
   
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) !GCC$ ATTRIBUTES hot :: DORGLQ !GCC$ ATTRIBUTES aligned(32) :: DORGLQ
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, &
                        LWKOPT, NB, NBMIN, NX
!*     ..
!*     .. External Subroutines ..
      !EXTERNAL           DLARFB, DLARFT, DORGL2, XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      NB = ILAENV( 1, 'DORGLQ', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, M )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DORGLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!*
!*     Quick return if possible
!*
     ! IF( M.LE.0 ) THEN
     !    WORK( 1 ) = 1
     !    RETURN
   !   END IF
!*
      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!*
!*        Determine when to cross over from blocked to unblocked code.
!*
         NX = MAX( 0, ILAENV( 3, 'DORGLQ', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
!*
!*           Determine if workspace is large enough for blocked code.
!*
            LDWORK = M
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!*
!*              Not enough workspace to use optimal NB:  reduce NB and
!*              determine the minimum value of NB.
!*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGLQ', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
!*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!*
!*        Use blocked code after the last block.
!*        The first kk rows are handled by the block method.
!!*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
!*
!*        Set A(kk+1:m,1:kk) to zero.
!*
         DO 20 J = 1, KK
            !$OMP SIMD ALIGNED(A:64)
            DO 10 I = KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
!*
!*     Use unblocked code for the last or only block.
!*
      IF( KK.LT.M ) &
        CALL DORGL2( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, &
                     TAU( KK+1 ), WORK, IINFO )
!*
      IF( KK.GT.0 ) THEN
!*
!*        Use blocked code
!*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.M ) THEN
!*
!*              Form the triangular factor of the block reflector
!*              H = H(i) H(i+1) . . . H(i+ib-1)
!*
               CALL DLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ), &
                           LDA, TAU( I ), WORK, LDWORK )
!*
!*              Apply H**T to A(i+ib:m,i:n) from the right
!*
               CALL DLARFB( 'Right', 'Transpose', 'Forward', 'Rowwise', &
                           M-I-IB+1, N-I+1, IB, A( I, I ), LDA, WORK, &
                           LDWORK, A( I+IB, I ), LDA, WORK( IB+1 ), &
                           LDWORK )
            END IF
!*
!*           Apply H**T to columns i:n of current block
!*
            CALL DORGL2( IB, N-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, &
                        IINFO )
!*
!*           Set columns 1:i-1 of current block to zero
!*
            DO 40 J = 1, I - 1
               !$OMP SIMD ALIGNED(A:64)
               DO 30 L = I, I + IB - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF

      WORK( 1 ) = IWS
    
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO ) !GCC$ ATTRIBUTES inline :: DORGL2 !GCC$ ATTRIBUTES aligned(32) :: DORGL2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
      
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, J, L
!*     ..
!*     .. External Subroutines ..
     ! EXTERNAL           DLARF, DSCAL, XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DORGL2', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
   !   IF( M.LE.0 )
    ! $   RETURN
!*
      IF( K.LT.M ) THEN
!*
!*        Initialise rows k+1:m to rows of the unit matrix
!*
         DO 20 J = 1, N
            DO 10 L = K + 1, M
               A( L, J ) = ZERO
   10       CONTINUE
            IF( J.GT.K .AND. J.LE.M ) &
               A( J, J ) = ONE
   20    CONTINUE
      END IF
!*
      DO 40 I = K, 1, -1
!*
!*        Apply H(i) to A(i:m,i:n) from the right
!*
         IF( I.LT.N ) THEN
            IF( I.LT.M ) THEN
               A( I, I ) = ONE
               CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, &
                          TAU( I ), A( I+1, I ), LDA, WORK )
            END IF
            CALL DSCAL( N-I, -TAU( I ), A( I, I+1 ), LDA )
         END IF
         A( I, I ) = ONE - TAU( I )
!*
!*        Set A(i,1:i-1) to zero
!*
         DO 30 L = 1, I - 1
            A( I, L ) = ZERO
   30    CONTINUE
   40 CONTINUE
    
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) !GCC$ ATTRIBUTES hot :: DORGQR !GCC$ ATTRIBUTES aligned(32) :: DORGQR
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, &
                        LWKOPT, NB, NBMIN, NX
!*     ..
!*     .. External Subroutines ..
     ! EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!*
!!*     Quick return if possible
!*
     ! IF( N.LE.0 ) THEN
     !    WORK( 1 ) = 1
     !    RETURN
     ! END IF
!*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!*
!*        Determine when to cross over from blocked to unblocked code.
!*
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
!*
!*           Determine if workspace is large enough for blocked code.
!*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!*
!*              Not enough workspace to use optimal NB:  reduce NB and
!*              determine the minimum value of NB.
!*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
!*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!*
!*        Use blocked code after the last block.
!*        The first kk columns are handled by the block method.
!*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
!*
!*        Set A(1:kk,kk+1:n) to zero.
!*
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
!*
!*     Use unblocked code for the last or only block.
!*
      IF( KK.LT.N ) &
        CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, &
                     TAU( KK+1 ), WORK, IINFO )

      IF( KK.GT.0 ) THEN
!*
!*        Use blocked code
!*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
!*
!*              Form the triangular factor of the block reflector
!*              H = H(i) H(i+1) . . . H(i+ib-1)
!*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                           A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!*
!*              Apply H to A(i:m,i+ib:n) from the left
!*
               CALL DLARFB( 'Left', 'No transpose', 'Forward', &
                           'Columnwise', M-I+1, N-I-IB+1, IB, &
                           A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                           LDA, WORK( IB+1 ), LDWORK )
            END IF
!*
!*           Apply H to rows i:m of current block
!*
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, &
                         IINFO )
!*
!*           Set rows 1:i-1 of current block to zero
!*
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF

      WORK( 1 ) = IWS
    
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO ) !GCC$ ATTRIBUTES inline :: DORG2R !GCC$ ATTRIBUTES aligned(32) :: DORG2R
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, J, L
!*     ..
!*     .. External Subroutines ..
      !EXTERNAL           DLARF, DSCAL, XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DORG2R', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
     ! IF( N.LE.0 )
    ! $   RETURN
!*
!*     Initialise columns k+1:n to columns of the unit matrix
!*
      DO 20 J = K + 1, N
         !$OMP SIMD ALIGNED(A:64) UNROLL PARTIAL(8)
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE

      DO 40 I = K, 1, -1
!*
!*        Apply H(i) to A(i:m,i:n) from the left
!*
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                       A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M ) &
           CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
!*
!*        Set A(1:i-1,i) to zero
         !*
           !$OMP SIMD ALIGNED(A:64) UNROLL PARTIAL(8)
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
  
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, &
     LDC, WORK, LWORK, INFO ) !GCC$ ATTRIBUTES hot :: DORMBR !GCC$ ATTRIBUTES aligned(32) :: DORMBR
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      LOGICAL            APPLYQ, LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           DORMLQ, DORMQR, XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      APPLYQ = LSAME( VECT, 'Q' )
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!*
!*     NQ is the order of Q or P and NW is the minimum dimension of WORK
!*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.APPLYQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( K.LT.0 ) THEN
         INFO = -6
      ELSE IF( ( APPLYQ .AND. LDA.LT.MAX( 1, NQ ) ) .OR. &
              ( .NOT.APPLYQ .AND. LDA.LT.MAX( 1, MIN( NQ, K ) ) ) ) &
             THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF

      IF( INFO.EQ.0 ) THEN
         IF( APPLYQ ) THEN
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M-1, N, M-1,-1)
    
            ELSE
               NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N-1, N-1,-1)
    
            END IF
         ELSE
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'DORMLQ', SIDE // TRANS, M-1, N, M-1,-1)
                
            ELSE
               NB = ILAENV( 1, 'DORMLQ', SIDE // TRANS, M, N-1, N-1,-1)
 
            END IF
         END IF
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMBR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      WORK( 1 ) = 1
   !   IF( M.EQ.0 .OR. N.EQ.0 )
   !  $   RETURN
!*
      IF( APPLYQ ) THEN
!*
!*        Apply Q
!*
         IF( NQ.GE.K ) THEN
!*
!*           Q was determined by a call to DGEBRD with nq >= k
!*
            CALL DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,&
                       WORK, LWORK, IINFO )
         ELSE IF( NQ.GT.1 ) THEN
!*
!*           Q was determined by a call to DGEBRD with nq < k
!*
            IF( LEFT ) THEN
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            END IF
            CALL DORMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, &
                        C( I1, I2 ), LDC, WORK, LWORK, IINFO )
         END IF
      ELSE
!!*
!*        Apply P
!*
         IF( NOTRAN ) THEN
            TRANST = 'T'
         ELSE
            TRANST = 'N'
         END IF
         IF( NQ.GT.K ) THEN

            CALL DORMLQ( SIDE, TRANST, M, N, K, A, LDA, TAU, C, LDC, &
                        WORK, LWORK, IINFO )
         ELSE IF( NQ.GT.1 ) THEN

            IF( LEFT ) THEN
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            END IF
            CALL DORMLQ( SIDE, TRANST, MI, NI, NQ-1, A( 1, 2 ), LDA, &
                       TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
     
END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleOTHERcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
     WORK, LWORK, INFO ) !GCC$ ATTRIBUTES hot :: DORMLQ !GCC$ ATTRIBUTES aligned(32) :: DORMLQ
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1, &
                         TSIZE = LDT*NBMAX )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK, &
                        LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!!*     ..
!*     .. External Subroutines ..
      !EXTERNAL           DLARFB, DLARFT, DORML2, XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!*
!*     NQ is the order of Q and NW is the minimum dimension of WORK
!*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!*
      IF( INFO.EQ.0 ) THEN
!*
!*        Compute the workspace requirements
!*
         NB = MIN( NBMAX, ILAENV( 1, 'DORMLQ', SIDE // TRANS, M, N, K, &
             -1 ) )
         LWKOPT = MAX( 1, NW )*NB + TSIZE
         WORK( 1 ) = LWKOPT
      END IF

      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DORMLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IF( LWORK.LT.NW*NB+TSIZE ) THEN
            NB = (LWORK-TSIZE) / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMLQ', SIDE // TRANS, M, N, K, &
                   -1 ) )
         END IF
      END IF

      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!*
!*        Use unblocked code
!*
         CALL DORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
                     IINFO )
      ELSE
!*
!*        Use blocked code
!*
         IWT = 1 + NW*NB
         IF( ( LEFT .AND. NOTRAN ) .OR. &
            ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF

         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF

         IF( NOTRAN ) THEN
            TRANST = 'T'
         ELSE
            TRANST = 'N'
         END IF

         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!*
!*           Form the triangular factor of the block reflector
!*           H = H(i) H(i+1) . . . H(i+ib-1)
!*
            CALL DLARFT( 'Forward', 'Rowwise', NQ-I+1, IB, A( I, I ), &
                        LDA, TAU( I ), WORK( IWT ), LDT )
            IF( LEFT ) THEN
!*
!*              H or H**T is applied to C(i:m,1:n)
!*
               MI = M - I + 1
               IC = I
            ELSE
!*
!*              H or H**T is applied to C(1:m,i:n)
!*
               NI = N - I + 1
               JC = I
            END IF
!*
!*           Apply H or H**T
!*
            CALL DLARFB( SIDE, TRANST, 'Forward', 'Rowwise', MI, NI, IB, &
                        A( I, I ), LDA, WORK( IWT ), LDT, &
                        C( IC, JC ), LDC, WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
     
END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleOTHERcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      SUBROUTINE DORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
           WORK, INFO ) !GCC$ ATTRIBUTES hot :: DORML2 !GCC$ ATTRIBUTES aligned(32) :: DORML2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           DLARF, XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!*
!*     NQ is the order of Q
!*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DORML2', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
   
!*
      IF( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) &
          THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF

      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF

      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!*
!*           H(i) is applied to C(i:m,1:n)
!*
            MI = M - I + 1
            IC = I
         ELSE
!*
!*           H(i) is applied to C(1:m,i:n)
!*
            NI = N - I + 1
            JC = I
         END IF
!*
!*        Apply H(i)
!*
         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), LDA, TAU( I ), &
                    C( IC, JC ), LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
     
END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleOTHERcomputational
!!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
 SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
      WORK, LWORK, INFO ) !GCC$ ATTRIBUTES hot :: DORMQR !GCC$ ATTRIBUTES aligned(32) :: DORMQR
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1, &
                          TSIZE = LDT*NBMAX )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK, &
                        LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!*     ..
!*     .. External Subroutines ..
     ! EXTERNAL           DLARFB, DLARFT, DORM2R, XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!*
!*     NQ is the order of Q and NW is the minimum dimension of WORK
!*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!*
      IF( INFO.EQ.0 ) THEN
!*
!*        Compute the workspace requirements
!*
         NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N, K, &
             -1 ) )
         LWKOPT = MAX( 1, NW )*NB + TSIZE
         WORK( 1 ) = LWKOPT
      END IF
!*
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DORMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IF( LWORK.LT.NW*NB+TSIZE ) THEN
            NB = (LWORK-TSIZE) / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N, K, &
                   -1 ) )
         END IF
      END IF

      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!*
!*        Use unblocked code
!*
         CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
                    IINFO )
      ELSE
!*
!*        Use blocked code
!*
         IWT = 1 + NW*NB
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. &
            ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF

         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF

         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!*
!*           Form the triangular factor of the block reflector
!*           H = H(i) H(i+1) . . . H(i+ib-1)
!*
            CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), &
                        LDA, TAU( I ), WORK( IWT ), LDT )
            IF( LEFT ) THEN
!*
!*              H or H**T is applied to C(i:m,1:n)
!*
               MI = M - I + 1
               IC = I
            ELSE
!*
!*              H or H**T is applied to C(1:m,i:n)
!*
               NI = N - I + 1
               JC = I
            END IF
!*
!*           Apply H or H**T
!*
            CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, &
                        IB, A( I, I ), LDA, WORK( IWT ), LDT, &
                        C( IC, JC ), LDC, WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
  
END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleOTHERcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
  SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
       WORK, INFO ) !GCC$ ATTRIBUTES hot :: DORM2R !GCC$ ATTRIBUTES aligned(32) :: DORM2R
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. External Subroutines ..
      !EXTERNAL           DLARF, XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!*
!*     NQ is the order of Q
!*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
       RETURN

      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF

      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF

      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!*
!*           H(i) is applied to C(i:m,1:n)
!*
            MI = M - I + 1
            IC = I
         ELSE
!*
!*           H(i) is applied to C(1:m,i:n)
!*
            NI = N - I + 1
            JC = I
         END IF
!*
!*        Apply H(i)
!*
         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ), &
                    LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
    
END SUBROUTINE


!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date April 2012
!*
!*> \ingroup doubleGEsolve
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DGESVXX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
                         EQUED, R, C, B, LDB, X, LDX, RCOND, RPVGRW, &
                         BERR, N_ERR_BNDS, ERR_BNDS_NORM, &
                         ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, &
                         INFO ) !GCC$ ATTRIBUTES hot :: DGESVXX !GCC$ ATTRIBUTES aligned(32) :: DGESVXX
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      INTEGER            FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I
      INTEGER            RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I
      INTEGER            CMP_ERR_I, PIV_GROWTH_I
      PARAMETER          ( FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, &
                        BERR_I = 3 )
      PARAMETER          ( RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 )
      PARAMETER          ( CMP_RCOND_I = 7, CMP_ERR_I = 8, &
                        PIV_GROWTH_I = 9 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU
      INTEGER            INFEQU, J
      DOUBLE PRECISION   AMAX, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, &
                        SMLNUM
!*     ..
!*     .. External Functions ..
      EXTERNAL           LSAME, DLAMCH, DLA_GERPVGRW
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLA_GERPVGRW
!*     ..
!*     .. External Subroutines ..
      !EXTERNAL           DGETRF, DGETRS, DLACPY, DLAQGE,
     !$                   DLASCL2, DGERFSX
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      NOTRAN = LSAME( TRANS, 'N' )
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         ROWEQU = .FALSE.
         COLEQU = .FALSE.
      ELSE
         ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
         COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
      END IF
!*
!*     Default is failure.  If an input parameter is wrong or
!*     factorization fails, make everything look horrible.  Only the
!*     pivot growth is set here, the rest is initialized in DGERFSX.
!*
      RPVGRW = ZERO
!*
!*     Test the input parameters.  PARAMS is not tested until DGERFSX.
!*
      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT. &
          LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
             LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT. &
             ( ROWEQU .OR. COLEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
         INFO = -10
      ELSE
         IF( ROWEQU ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 10 J = 1, N
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
 10         CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -11
            ELSE IF( N.GT.0 ) THEN
               ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
               ROWCND = ONE
            END IF
         END IF
         IF( COLEQU .AND. INFO.EQ.0 ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 20 J = 1, N
               RCMIN = MIN( RCMIN, C( J ) )
               RCMAX = MAX( RCMAX, C( J ) )
 20         CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -12
            ELSE IF( N.GT.0 ) THEN
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
               COLCND = ONE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( LDB.LT.MAX( 1, N ) ) THEN
               INFO = -14
            ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
               INFO = -16
            END IF
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DGESVXX', -INFO )
         RETURN
      END IF

      IF( EQUIL ) THEN
!*
!*     Compute row and column scalings to equilibrate the matrix A.
!*
         CALL DGEEQUB( N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, &
             INFEQU )
         IF( INFEQU.EQ.0 ) THEN
!*
!*     Equilibrate the matrix.
!*
            CALL DLAQGE( N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, &
                EQUED )
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         END IF
!*
!*     If the scaling factors are not applied, set them to 1.0.
!*
         IF ( .NOT.ROWEQU ) THEN
            DO J = 1, N
               R( J ) = 1.0D+0
            END DO
         END IF
         IF ( .NOT.COLEQU ) THEN
            DO J = 1, N
               C( J ) = 1.0D+0
            END DO
         END IF
      END IF
!*
!*     Scale the right-hand side.
!*
      IF( NOTRAN ) THEN
         IF( ROWEQU ) CALL DLASCL2( N, NRHS, R, B, LDB )
      ELSE
         IF( COLEQU ) CALL DLASCL2( N, NRHS, C, B, LDB )
      END IF
!*
      IF( NOFACT .OR. EQUIL ) THEN
!*
!*        Compute the LU factorization of A.
!*
         CALL DLACPY( 'Full', N, N, A, LDA, AF, LDAF )
         CALL DGETRF( N, N, AF, LDAF, IPIV, INFO )
!*
!*        Return if INFO is non-zero.
!*
         IF( INFO.GT.0 ) THEN
!*
!!*           Pivot in column INFO is exactly 0
!*           Compute the reciprocal pivot growth factor of the
!*           leading rank-deficient INFO columns of A.
!*
            RPVGRW = DLA_GERPVGRW( N, INFO, A, LDA, AF, LDAF )
            RETURN
         END IF
      END IF
!*
!*     Compute the reciprocal pivot growth factor RPVGRW.
!*
      RPVGRW = DLA_GERPVGRW( N, N, A, LDA, AF, LDAF )
!*
!*     Compute the solution matrix X.
!*
      CALL DLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL DGETRS( TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO )
!*
!*     Use iterative refinement to improve the computed solution and
!*     compute error bounds and backward error estimates for it.
!*
      CALL DGERFSX( TRANS, EQUED, N, NRHS, A, LDA, AF, LDAF, &
          IPIV, R, C, B, LDB, X, LDX, RCOND, BERR, &
          N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, &
          WORK, IWORK, INFO )
!*
!*     Scale solutions.
!*
      IF ( COLEQU .AND. NOTRAN ) THEN
         CALL DLASCL2 ( N, NRHS, C, X, LDX )
      ELSE IF ( ROWEQU .AND. .NOT.NOTRAN ) THEN
         CALL DLASCL2 ( N, NRHS, R, X, LDX )
      END IF


END SUBROUTINE

!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleGEsolve
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, &
     WORK, LWORK, INFO ) !GCC$ ATTRIBUTES hot :: DGELSS !GCC$ ATTRIBUTES aligned(32) :: DGELSS
#elif defined(__INTEL_COMPILER) || defined(__ICC)
SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, &
     WORK, LWORK, INFO )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGELSS
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DGELSS
#endif
      implicit none
!*
!*  -- LAPACK driver routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
      DOUBLE PRECISION   RCOND
!*     ..
!*!     .. Array Arguments ..
      ! DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: S
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
      
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            BDSPAC, BL, CHUNK, I, IASCL, IBSCL, IE, IL, &
                        ITAU, ITAUP, ITAUQ, IWORK, LDWORK, MAXMN, &
                        MAXWRK, MINMN, MINWRK, MM, MNTHR
      INTEGER            LWORK_DGEQRF, LWORK_DORMQR, LWORK_DGEBRD, &
                        LWORK_DORMBR, LWORK_DORGBR, LWORK_DORMLQ, &
                        LWORK_DGELQF
      DOUBLE PRECISION   ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM, THR
!*     ..
!*     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 )
!*     ..
!*     .. External Subroutines ..
      EXTERNAL            DCOPY, DGEMM, DGEMV
     !$                   DGEQRF, DLABAD, DLACPY, DLASCL, DLASET, DORGBR,
     !$                   DORMBR, DORMLQ, DORMQR, DRSCL
!*     ..
!*     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           ILAENV, DLAMCH, DLANGE
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments
!*
      INFO = 0
      MINMN = MIN( M, N )
      MAXMN = MAX( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, MAXMN ) ) THEN
         INFO = -7
      END IF
!*
!*     Compute workspace
!*      (Note: Comments in the code beginning "Workspace:" describe the
!*       minimal amount of workspace needed at that point in the code,
!*       as well as the preferred amount for good performance.
!*       NB refers to the optimal block size for the immediately
!*       following subroutine, as returned by ILAENV.)
!*
      IF( INFO.EQ.0 ) THEN
         MINWRK = 1
         MAXWRK = 1
         IF( MINMN.GT.0 ) THEN
            MM = M
            MNTHR = ILAENV( 6, 'DGELSS', ' ', M, N, NRHS, -1 )
            IF( M.GE.N .AND. M.GE.MNTHR ) THEN
!*
!*              Path 1a - overdetermined, with many more rows than
!*                        columns
!*
!*              Compute space needed for DGEQRF
               CALL DGEQRF( M, N, A, LDA, DUM(1), DUM(1), -1, INFO )
               LWORK_DGEQRF=DUM(1)
!*              Compute space needed for DORMQR
               CALL DORMQR( 'L', 'T', M, NRHS, N, A, LDA, DUM(1), B, &
                        LDB, DUM(1), -1, INFO )
               LWORK_DORMQR=DUM(1)
               MM = N
               MAXWRK = MAX( MAXWRK, N + LWORK_DGEQRF )
               MAXWRK = MAX( MAXWRK, N + LWORK_DORMQR )
            END IF
            IF( M.GE.N ) THEN
!*
!*              Path 1 - overdetermined or exactly determined
!*
!*              Compute workspace needed for DBDSQR
!*
               BDSPAC = MAX( 1, 5*N )
!*              Compute space needed for DGEBRD
               CALL DGEBRD( MM, N, A, LDA, S, DUM(1), DUM(1), &
                           DUM(1), DUM(1), -1, INFO )
               LWORK_DGEBRD=DUM(1)
!*              Compute space needed for DORMBR
               CALL DORMBR( 'Q', 'L', 'T', MM, NRHS, N, A, LDA, DUM(1), &
                      B, LDB, DUM(1), -1, INFO )
               LWORK_DORMBR=DUM(1)
!*              Compute space needed for DORGBR
               CALL DORGBR( 'P', N, N, N, A, LDA, DUM(1), &
                        DUM(1), -1, INFO )
               LWORK_DORGBR=DUM(1)
!*              Compute total workspace needed
               MAXWRK = MAX( MAXWRK, 3*N + LWORK_DGEBRD )
               MAXWRK = MAX( MAXWRK, 3*N + LWORK_DORMBR )
               MAXWRK = MAX( MAXWRK, 3*N + LWORK_DORGBR )
               MAXWRK = MAX( MAXWRK, BDSPAC )
               MAXWRK = MAX( MAXWRK, N*NRHS )
               MINWRK = MAX( 3*N + MM, 3*N + NRHS, BDSPAC )
               MAXWRK = MAX( MINWRK, MAXWRK )
            END IF
            IF( N.GT.M ) THEN
!*
!*              Compute workspace needed for DBDSQR
!*
               BDSPAC = MAX( 1, 5*M )
               MINWRK = MAX( 3*M+NRHS, 3*M+N, BDSPAC )
               IF( N.GE.MNTHR ) THEN
!*
!*                 Path 2a - underdetermined, with many more columns
!*                 than rows
!*
!*                 Compute space needed for DGELQF
                  CALL DGELQF( M, N, A, LDA, DUM(1), DUM(1), &
                     -1, INFO )
                  LWORK_DGELQF=DUM(1)
                 Compute space needed for DGEBRD
                  CALL DGEBRD( M, M, A, LDA, S, DUM(1), DUM(1), &
                           DUM(1), DUM(1), -1, INFO )
                  LWORK_DGEBRD=DUM(1)
                 Compute space needed for DORMBR
                  CALL DORMBR( 'Q', 'L', 'T', M, NRHS, N, A, LDA, &
                     DUM(1), B, LDB, DUM(1), -1, INFO )
                  LWORK_DORMBR=DUM(1)
                 Compute space needed for DORGBR
                  CALL DORGBR( 'P', M, M, M, A, LDA, DUM(1), &
                        DUM(1), -1, INFO )
                  LWORK_DORGBR=DUM(1)
                 Compute space needed for DORMLQ
                  CALL DORMLQ( 'L', 'T', N, NRHS, M, A, LDA, DUM(1), &
                      B, LDB, DUM(1), -1, INFO )
                  LWORK_DORMLQ=DUM(1)
                 Compute total workspace needed
                  MAXWRK = M + LWORK_DGELQF
                  MAXWRK = MAX( MAXWRK, M*M + 4*M + LWORK_DGEBRD )
                  MAXWRK = MAX( MAXWRK, M*M + 4*M + LWORK_DORMBR )
                  MAXWRK = MAX( MAXWRK, M*M + 4*M + LWORK_DORGBR )
                  MAXWRK = MAX( MAXWRK, M*M + M + BDSPAC )
                  IF( NRHS.GT.1 ) THEN
                     MAXWRK = MAX( MAXWRK, M*M + M + M*NRHS )
                  ELSE
                     MAXWRK = MAX( MAXWRK, M*M + 2*M )
                  END IF
                  MAXWRK = MAX( MAXWRK, M + LWORK_DORMLQ )
               ELSE
!*
!*                 Path 2 - underdetermined
!*
!*                 Compute space needed for DGEBRD
                  CALL DGEBRD( M, N, A, LDA, S, DUM(1), DUM(1), &
                           DUM(1), DUM(1), -1, INFO )
                  LWORK_DGEBRD=DUM(1)
!*                 Compute space needed for DORMBR
                  CALL DORMBR( 'Q', 'L', 'T', M, NRHS, M, A, LDA, &
                     DUM(1), B, LDB, DUM(1), -1, INFO )
                  LWORK_DORMBR=DUM(1)
!*                 Compute space needed for DORGBR
                  CALL DORGBR( 'P', M, N, M, A, LDA, DUM(1), &
                        DUM(1), -1, INFO )
                  LWORK_DORGBR=DUM(1)
                  MAXWRK = 3*M + LWORK_DGEBRD
                  MAXWRK = MAX( MAXWRK, 3*M + LWORK_DORMBR )
                  MAXWRK = MAX( MAXWRK, 3*M + LWORK_DORGBR )
                  MAXWRK = MAX( MAXWRK, BDSPAC )
                  MAXWRK = MAX( MAXWRK, N*NRHS )
               END IF
            END IF
            MAXWRK = MAX( MINWRK, MAXWRK )
         END IF
         WORK( 1 ) = MAXWRK
!*
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) &
           INFO = -12
      END IF
!*
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DGELSS', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         RANK = 0
         RETURN
      END IF
!*
!*     Get machine parameters
!*
      EPS = DLAMCH( 'P' )
      SFMIN = DLAMCH( 'S' )
      SMLNUM = SFMIN / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
!*
!*     Scale A if max element outside range [SMLNUM,BIGNUM]
!*
      ANRM = DLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
!*
!*        Scale matrix norm up to SMLNUM
!*
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
!*
!*        Scale matrix norm down to BIGNUM
!*
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
!*
!*        Matrix all zero. Return zero solution.
!*
         CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         CALL DLASET( 'F', MINMN, 1, ZERO, ZERO, S, MINMN )
         RANK = 0
         GO TO 70
      END IF
!*
!*     Scale B if max element outside range [SMLNUM,BIGNUM]
!*
      BNRM = DLANGE( 'M', M, NRHS, B, LDB, WORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
!*
!*        Scale matrix norm up to SMLNUM
!*
         CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
!*
!*        Scale matrix norm down to BIGNUM
!*
         CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 2
      END IF
!*
!*     Overdetermined case
!*
      IF( M.GE.N ) THEN
!*
!*        Path 1 - overdetermined or exactly determined
!*
         MM = M
         IF( M.GE.MNTHR ) THEN
!*
!*           Path 1a - overdetermined, with many more rows than columns
!*
            MM = N
            ITAU = 1
            IWORK = ITAU + N
!*
!*           Compute A=Q*R
!*           (Workspace: need 2*N, prefer N+N*NB)
!*
            CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                        LWORK-IWORK+1, INFO )
!*
!*           Multiply B by transpose(Q)
!*           (Workspace: need N+NRHS, prefer N+NRHS*NB)
!*
            CALL DORMQR( 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAU ), B, &
                        LDB, WORK( IWORK ), LWORK-IWORK+1, INFO )
!*
!*           Zero out below R
!*
            IF( N.GT.1 ) &
              CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )
         END IF
!*
         IE = 1
         ITAUQ = IE + N
         ITAUP = ITAUQ + N
         IWORK = ITAUP + N
!*
!*        Bidiagonalize R in A
!*        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)
!*
         CALL DGEBRD( MM, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                     WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, &
                     INFO )
!*
!*        Multiply B by transpose of left bidiagonalizing vectors of R
!*        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)
!*
         CALL DORMBR( 'Q', 'L', 'T', MM, NRHS, N, A, LDA, WORK( ITAUQ ), &
                     B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO )
!*
!*        Generate right bidiagonalizing vectors of R in A
!*        (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
!*
         CALL DORGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ), &
                     WORK( IWORK ), LWORK-IWORK+1, INFO )
         IWORK = IE + N
!*
!*        Perform bidiagonal QR iteration
!*          multiply B by transpose of left singular vectors
!*          compute right singular vectors in A
!*        (Workspace: need BDSPAC)
!*
         CALL DBDSQR( 'U', N, N, 0, NRHS, S, WORK( IE ), A, LDA, DUM, &
                     1, B, LDB, WORK( IWORK ), INFO )
         IF( INFO.NE.0 ) &
           GO TO 70
!*
!*        Multiply B by reciprocals of singular values
!*
         THR = MAX( RCOND*S( 1 ), SFMIN )
         IF( RCOND.LT.ZERO ) &
           THR = MAX( EPS*S( 1 ), SFMIN )
         RANK = 0
         DO 10 I = 1, N
            IF( S( I ).GT.THR ) THEN
               CALL DRSCL( NRHS, S( I ), B( I, 1 ), LDB )
               RANK = RANK + 1
            ELSE
               CALL DLASET( 'F', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB )
            END IF
   10    CONTINUE
!*
!*        Multiply B by right singular vectors
!*        (Workspace: need N, prefer N*NRHS)
!*
         IF( LWORK.GE.LDB*NRHS .AND. NRHS.GT.1 ) THEN
            CALL DGEMM( 'T', 'N', N, NRHS, N, ONE, A, LDA, B, LDB, ZERO, &
                       WORK, LDB )
            CALL DLACPY( 'G', N, NRHS, WORK, LDB, B, LDB )
         ELSE IF( NRHS.GT.1 ) THEN
            CHUNK = LWORK / N
            DO 20 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               CALL DGEMM( 'T', 'N', N, BL, N, ONE, A, LDA, B( 1, I ), &
                          LDB, ZERO, WORK, N )
               CALL DLACPY( 'G', N, BL, WORK, N, B( 1, I ), LDB )
   20       CONTINUE
         ELSE
            CALL DGEMV( 'T', N, N, ONE, A, LDA, B, 1, ZERO, WORK, 1 )
            CALL DCOPY( N, WORK, 1, B, 1 )
         END IF

      ELSE IF( N.GE.MNTHR .AND. LWORK.GE.4*M+M*M+ &
              MAX( M, 2*M-4, NRHS, N-3*M ) ) THEN
!*
!*        Path 2a - underdetermined, with many more columns than rows
!*        and sufficient workspace for an efficient algorithm
!*
         LDWORK = M
         IF( LWORK.GE.MAX( 4*M+M*LDA+MAX( M, 2*M-4, NRHS, N-3*M ), &
            M*LDA+M+M*NRHS ) )LDWORK = LDA
         ITAU = 1
         IWORK = M + 1
!*
!*        Compute A=L*Q
!*        (Workspace: need 2*M, prefer M+M*NB)
!*
         CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                     LWORK-IWORK+1, INFO )
         IL = IWORK
!*
!*        Copy L to WORK(IL), zeroing out above it
!*
         CALL DLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWORK )
         CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IL+LDWORK ), &
                     LDWORK )
         IE = IL + LDWORK*M
         ITAUQ = IE + M
         ITAUP = ITAUQ + M
         IWORK = ITAUP + M
!*
!*        Bidiagonalize L in WORK(IL)
!*        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)
!*
         CALL DGEBRD( M, M, WORK( IL ), LDWORK, S, WORK( IE ), &
                     WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), &
                     LWORK-IWORK+1, INFO )
!*
!*        Multiply B by transpose of left bidiagonalizing vectors of L
!*        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
!*
         CALL DORMBR( 'Q', 'L', 'T', M, NRHS, M, WORK( IL ), LDWORK, &
                     WORK( ITAUQ ), B, LDB, WORK( IWORK ), &
                     LWORK-IWORK+1, INFO )
!*
!*        Generate right bidiagonalizing vectors of R in WORK(IL)
!*        (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB)
!*
         CALL DORGBR( 'P', M, M, M, WORK( IL ), LDWORK, WORK( ITAUP ), &
                     WORK( IWORK ), LWORK-IWORK+1, INFO )
         IWORK = IE + M
!*
!*        Perform bidiagonal QR iteration,
!*           computing right singular vectors of L in WORK(IL) and
!*           multiplying B by transpose of left singular vectors
!*        (Workspace: need M*M+M+BDSPAC)
!*
         CALL DBDSQR( 'U', M, M, 0, NRHS, S, WORK( IE ), WORK( IL ), &
                     LDWORK, A, LDA, B, LDB, WORK( IWORK ), INFO )
         IF( INFO.NE.0 ) &
           GO TO 70
!*
!*        Multiply B by reciprocals of singular values
!*
         THR = MAX( RCOND*S( 1 ), SFMIN )
         IF( RCOND.LT.ZERO ) &
           THR = MAX( EPS*S( 1 ), SFMIN )
         RANK = 0
         DO 30 I = 1, M
            IF( S( I ).GT.THR ) THEN
               CALL DRSCL( NRHS, S( I ), B( I, 1 ), LDB )
               RANK = RANK + 1
            ELSE
               CALL DLASET( 'F', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB )
            END IF
   30    CONTINUE
         IWORK = IE
!*
!*        Multiply B by right singular vectors of L in WORK(IL)
!*        (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS)
!*
         IF( LWORK.GE.LDB*NRHS+IWORK-1 .AND. NRHS.GT.1 ) THEN
            CALL DGEMM( 'T', 'N', M, NRHS, M, ONE, WORK( IL ), LDWORK, &
                       B, LDB, ZERO, WORK( IWORK ), LDB )
            CALL DLACPY( 'G', M, NRHS, WORK( IWORK ), LDB, B, LDB )
         ELSE IF( NRHS.GT.1 ) THEN
            CHUNK = ( LWORK-IWORK+1 ) / M
            DO 40 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               CALL DGEMM( 'T', 'N', M, BL, M, ONE, WORK( IL ), LDWORK, &
                          B( 1, I ), LDB, ZERO, WORK( IWORK ), M )
               CALL DLACPY( 'G', M, BL, WORK( IWORK ), M, B( 1, I ), &
                           LDB )
   40       CONTINUE
         ELSE
            CALL DGEMV( 'T', M, M, ONE, WORK( IL ), LDWORK, B( 1, 1 ), &
                       1, ZERO, WORK( IWORK ), 1 )
            CALL DCOPY( M, WORK( IWORK ), 1, B( 1, 1 ), 1 )
         END IF
!*
!*        Zero out below first M rows of B
!*
         CALL DLASET( 'F', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB )
         IWORK = ITAU + M
!*
!*        Multiply transpose(Q) by B
!*        (Workspace: need M+NRHS, prefer M+NRHS*NB)
!*
         CALL DORMLQ( 'L', 'T', N, NRHS, M, A, LDA, WORK( ITAU ), B, &
                     LDB, WORK( IWORK ), LWORK-IWORK+1, INFO )
!*
      ELSE
!*
!*        Path 2 - remaining underdetermined cases
!*
         IE = 1
         ITAUQ = IE + M
         ITAUP = ITAUQ + M
         IWORK = ITAUP + M
!*
!*        Bidiagonalize A
!*        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
!*
         CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ),  &
                     WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, &
                     INFO )
!*
!*        Multiply B by transpose of left bidiagonalizing vectors
!*        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)
!*
         CALL DORMBR( 'Q', 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAUQ ), &
                     B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO )
!*
!*        Generate right bidiagonalizing vectors in A
!*        (Workspace: need 4*M, prefer 3*M+M*NB)
!*
         CALL DORGBR( 'P', M, N, M, A, LDA, WORK( ITAUP ), &
                     WORK( IWORK ), LWORK-IWORK+1, INFO )
         IWORK = IE + M
!*
!*        Perform bidiagonal QR iteration,
!*           computing right singular vectors of A in A and
!*           multiplying B by transpose of left singular vectors
!*        (Workspace: need BDSPAC)
!*
         CALL DBDSQR( 'L', M, N, 0, NRHS, S, WORK( IE ), A, LDA, DUM, &
                     1, B, LDB, WORK( IWORK ), INFO )
         IF( INFO.NE.0 ) &
           GO TO 70
!*
!*        Multiply B by reciprocals of singular values
!*
         THR = MAX( RCOND*S( 1 ), SFMIN )
         IF( RCOND.LT.ZERO ) &
          THR = MAX( EPS*S( 1 ), SFMIN )
         RANK = 0
         DO 50 I = 1, M
            IF( S( I ).GT.THR ) THEN
               CALL DRSCL( NRHS, S( I ), B( I, 1 ), LDB )
               RANK = RANK + 1
            ELSE
               CALL DLASET( 'F', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB )
            END IF
   50    CONTINUE
!*
!*        Multiply B by right singular vectors of A
!*        (Workspace: need N, prefer N*NRHS)
!*
         IF( LWORK.GE.LDB*NRHS .AND. NRHS.GT.1 ) THEN
            CALL DGEMM( 'T', 'N', N, NRHS, M, ONE, A, LDA, B, LDB, ZERO, &
                       WORK, LDB )
            CALL DLACPY( 'F', N, NRHS, WORK, LDB, B, LDB )
         ELSE IF( NRHS.GT.1 ) THEN
            CHUNK = LWORK / N
            DO 60 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               CALL DGEMM( 'T', 'N', N, BL, M, ONE, A, LDA, B( 1, I ), &
                          LDB, ZERO, WORK, N )
               CALL DLACPY( 'F', N, BL, WORK, N, B( 1, I ), LDB )
   60       CONTINUE
         ELSE
            CALL DGEMV( 'T', M, N, ONE, A, LDA, B, 1, ZERO, WORK, 1 )
            CALL DCOPY( N, WORK, 1, B, 1 )
         END IF
      END IF
!*
!*     Undo scaling
!*
      IF( IASCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, &
                     INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, &
                     INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      END IF
      
   70 CONTINUE
      WORK( 1 ) = MAXWRK
  
    END SUBROUTINE

!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*!> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date June 2017
!*
!*> \ingroup doubleGEsolve
!*
!*> \par Contributors:
!*  ==================
!*>
!*>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!!*>       California at Berkeley, USA \n
!*>     Osni Marques, LBNL/NERSC, USA \n
!*
    !*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
 SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, &
      WORK, LWORK, IWORK, INFO ) !GCC$ ATTRIBUTES hot :: DGELSD !GCC$ ATTRIBUTES aligned(32) :: DGELSD
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            IASCL, IBSCL, IE, IL, ITAU, ITAUP, ITAUQ, &
                        LDWORK, LIWORK, MAXMN, MAXWRK, MINMN, MINWRK, &
                       MM, MNTHR, NLVL, NWORK, SMLSIZ, WLALSD
      DOUBLE PRECISION   ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM
!*     ..
!*     .. External Subroutines ..
     ! EXTERNAL           DGEBRD, DGELQF, DGEQRF, DLABAD, DLACPY, DLALSD,
     !$                   DLASCL, DLASET, DORMBR, DORMLQ, DORMQR, XERBLA
!*     ..
!*     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           ILAENV, DLAMCH, DLANGE
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, LOG, MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input arguments.
!*
      INFO = 0
      MINMN = MIN( M, N )
      MAXMN = MAX( M, N )
      MNTHR = ILAENV( 6, 'DGELSD', ' ', M, N, NRHS, -1 )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, MAXMN ) ) THEN
         INFO = -7
      END IF
!*
      SMLSIZ = ILAENV( 9, 'DGELSD', ' ', 0, 0, 0, 0 )
!*
!!*     Compute workspace.
!*     (Note: Comments in the code beginning "Workspace:" describe the
!*     minimal amount of workspace needed at that point in the code,
!*     as well as the preferred amount for good performance.
!*     NB refers to the optimal block size for the immediately
!*     following subroutine, as returned by ILAENV.)
!*
      MINWRK = 1
      LIWORK = 1
      MINMN = MAX( 1, MINMN )
      NLVL = MAX( INT( LOG( DBLE( MINMN ) / DBLE( SMLSIZ+1 ) ) / &
            LOG( TWO ) ) + 1, 0 )

      IF( INFO.EQ.0 ) THEN
         MAXWRK = 0
         LIWORK = 3*MINMN*NLVL + 11*MINMN
         MM = M
         IF( M.GE.N .AND. M.GE.MNTHR ) THEN
!*
!*           Path 1a - overdetermined, with many more rows than columns.
!*
            MM = N
            MAXWRK = MAX( MAXWRK, N+N*ILAENV( 1, 'DGEQRF', ' ', M, N, &
                    -1, -1 ) )
            MAXWRK = MAX( MAXWRK, N+NRHS* &
                    ILAENV( 1, 'DORMQR', 'LT', M, NRHS, N, -1 ) )
         END IF
         IF( M.GE.N ) THEN
!*
!*           Path 1 - overdetermined or exactly determined.
!*
            MAXWRK = MAX( MAXWRK, 3*N+( MM+N )* &
                    ILAENV( 1, 'DGEBRD', ' ', MM, N, -1, -1 ) )
            MAXWRK = MAX( MAXWRK, 3*N+NRHS* &
                    ILAENV( 1, 'DORMBR', 'QLT', MM, NRHS, N, -1 ) )
            MAXWRK = MAX( MAXWRK, 3*N+( N-1 )* &
                    ILAENV( 1, 'DORMBR', 'PLN', N, NRHS, N, -1 ) )
            WLALSD = 9*N+2*N*SMLSIZ+8*N*NLVL+N*NRHS+(SMLSIZ+1)**2
            MAXWRK = MAX( MAXWRK, 3*N+WLALSD )
            MINWRK = MAX( 3*N+MM, 3*N+NRHS, 3*N+WLALSD )
         END IF
         IF( N.GT.M ) THEN
            WLALSD = 9*M+2*M*SMLSIZ+8*M*NLVL+M*NRHS+(SMLSIZ+1)**2
            IF( N.GE.MNTHR ) THEN
!*
!*              Path 2a - underdetermined, with many more columns
!*              than rows.
!*
               MAXWRK = M + M*ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
               MAXWRK = MAX( MAXWRK, M*M+4*M+2*M* &
                       ILAENV( 1, 'DGEBRD', ' ', M, M, -1, -1 ) )
               MAXWRK = MAX( MAXWRK, M*M+4*M+NRHS* &
                        ILAENV( 1, 'DORMBR', 'QLT', M, NRHS, M, -1 ) )
               MAXWRK = MAX( MAXWRK, M*M+4*M+( M-1 )* &
                       ILAENV( 1, 'DORMBR', 'PLN', M, NRHS, M, -1 ) )
               IF( NRHS.GT.1 ) THEN
                  MAXWRK = MAX( MAXWRK, M*M+M+M*NRHS )
               ELSE
                  MAXWRK = MAX( MAXWRK, M*M+2*M )
               END IF
               MAXWRK = MAX( MAXWRK, M+NRHS* &
                       ILAENV( 1, 'DORMLQ', 'LT', N, NRHS, M, -1 ) )
               MAXWRK = MAX( MAXWRK, M*M+4*M+WLALSD )
!     XXX: Ensure the Path 2a case below is triggered.  The workspace
!     calculation should use queries for all routines eventually.
               MAXWRK = MAX( MAXWRK, &
                   4*M+M*M+MAX( M, 2*M-4, NRHS, N-3*M ) )
            ELSE
!*
!*              Path 2 - remaining underdetermined cases.
!*
               MAXWRK = 3*M + ( N+M )*ILAENV( 1, 'DGEBRD', ' ', M, N, &
                       -1, -1 )
               MAXWRK = MAX( MAXWRK, 3*M+NRHS* &
                       ILAENV( 1, 'DORMBR', 'QLT', M, NRHS, N, -1 ) )
               MAXWRK = MAX( MAXWRK, 3*M+M* &
                       ILAENV( 1, 'DORMBR', 'PLN', N, NRHS, M, -1 ) )
               MAXWRK = MAX( MAXWRK, 3*M+WLALSD )
            END IF
            MINWRK = MAX( 3*M+NRHS, 3*M+M, 3*M+WLALSD )
         END IF
         MINWRK = MIN( MINWRK, MAXWRK )
         WORK( 1 ) = MAXWRK
         IWORK( 1 ) = LIWORK

         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGELSD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         GO TO 10
      END IF
!*
!*     Quick return if possible.
!*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         RANK = 0
         RETURN
      END IF
!*
!*     Get machine parameters.
!*
      EPS = DLAMCH( 'P' )
      SFMIN = DLAMCH( 'S' )
      SMLNUM = SFMIN / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
!*
!*     Scale A if max entry outside range [SMLNUM,BIGNUM].
!*
      ANRM = DLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
!*
!*        Scale matrix norm up to SMLNUM.
!*
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
!*
!*        Scale matrix norm down to BIGNUM.
!*
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
!*
!*        Matrix all zero. Return zero solution.
!*
         CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         CALL DLASET( 'F', MINMN, 1, ZERO, ZERO, S, 1 )
         RANK = 0
         GO TO 10
      END IF
!*
!*     Scale B if max entry outside range [SMLNUM,BIGNUM].
!*
      BNRM = DLANGE( 'M', M, NRHS, B, LDB, WORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
!*
!*        Scale matrix norm up to SMLNUM.
!*
         CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
!*
!*        Scale matrix norm down to BIGNUM.
!*
         CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 2
      END IF
!*
!*     If M < N make sure certain entries of B are zero.
!*
      IF( M.LT.N ) &
        CALL DLASET( 'F', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB )
!*
!*     Overdetermined case.
!*
      IF( M.GE.N ) THEN
!*
!*        Path 1 - overdetermined or exactly determined.
!*
         MM = M
         IF( M.GE.MNTHR ) THEN
!*
!*           Path 1a - overdetermined, with many more rows than columns.
!*
            MM = N
            ITAU = 1
            NWORK = ITAU + N
!*
!*           Compute A=Q*R.
!*           (Workspace: need 2*N, prefer N+N*NB)
!*
            CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), &
                        LWORK-NWORK+1, INFO )
!*
!*           Multiply B by transpose(Q).
!*           (Workspace: need N+NRHS, prefer N+NRHS*NB)
!*
            CALL DORMQR( 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAU ), B, &
                        LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!*
!*           Zero out below R.
!*
            IF( N.GT.1 ) THEN
               CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )
            END IF
         END IF
!*
         IE = 1
         ITAUQ = IE + N
         ITAUP = ITAUQ + N
         NWORK = ITAUP + N
!*
!*        Bidiagonalize R in A.
!*        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)
!*
         CALL DGEBRD( MM, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                     WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, &
                     INFO )
!*!
!*        Multiply B by transpose of left bidiagonalizing vectors of R.
!*        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)
!*
         CALL DORMBR( 'Q', 'L', 'T', MM, NRHS, N, A, LDA, WORK( ITAUQ ), &
                     B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!*
!*        Solve the bidiagonal least squares problem.
!*
         CALL DLALSD( 'U', SMLSIZ, N, NRHS, S, WORK( IE ), B, LDB, &
                     RCOND, RANK, WORK( NWORK ), IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF
!*
!*        Multiply B by right bidiagonalizing vectors of R.
!*
         CALL DORMBR( 'P', 'L', 'N', N, NRHS, N, A, LDA, WORK( ITAUP ), &
                     B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!*
      ELSE IF( N.GE.MNTHR .AND. LWORK.GE.4*M+M*M+ &
              MAX( M, 2*M-4, NRHS, N-3*M, WLALSD ) ) THEN
!*
!*        Path 2a - underdetermined, with many more columns than rows
!*        and sufficient workspace for an efficient algorithm.
!*
         LDWORK = M
         IF( LWORK.GE.MAX( 4*M+M*LDA+MAX( M, 2*M-4, NRHS, N-3*M ), &
            M*LDA+M+M*NRHS, 4*M+M*LDA+WLALSD ) )LDWORK = LDA
         ITAU = 1
         NWORK = M + 1
!*
!*        Compute A=L*Q.
!*        (Workspace: need 2*M, prefer M+M*NB)
!*
         CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), &
                     LWORK-NWORK+1, INFO )
         IL = NWORK
!*
!*        Copy L to WORK(IL), zeroing out above its diagonal.
!*
         CALL DLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWORK )
         CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IL+LDWORK ), &
                     LDWORK )
         IE = IL + LDWORK*M
         ITAUQ = IE + M
         ITAUP = ITAUQ + M
         NWORK = ITAUP + M
!*
!*        Bidiagonalize L in WORK(IL).
!*        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)
!*
         CALL DGEBRD( M, M, WORK( IL ), LDWORK, S, WORK( IE ),&
                     WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), &
                     LWORK-NWORK+1, INFO )
!*
!*        Multiply B by transpose of left bidiagonalizing vectors of L.
!*        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
!*
         CALL DORMBR( 'Q', 'L', 'T', M, NRHS, M, WORK( IL ), LDWORK, &
                     WORK( ITAUQ ), B, LDB, WORK( NWORK ), &
                     LWORK-NWORK+1, INFO )
!*
!*        Solve the bidiagonal least squares problem.
!*
         CALL DLALSD( 'U', SMLSIZ, M, NRHS, S, WORK( IE ), B, LDB, &
                     RCOND, RANK, WORK( NWORK ), IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF
!*
!*        Multiply B by right bidiagonalizing vectors of L.
!*
         CALL DORMBR( 'P', 'L', 'N', M, NRHS, M, WORK( IL ), LDWORK, &
                     WORK( ITAUP ), B, LDB, WORK( NWORK ), &
                     LWORK-NWORK+1, INFO )
!*
!*        Zero out below first M rows of B.
!*
         CALL DLASET( 'F', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB )
         NWORK = ITAU + M
!*
!*        Multiply transpose(Q) by B.
!*        (Workspace: need M+NRHS, prefer M+NRHS*NB)
!*
         CALL DORMLQ( 'L', 'T', N, NRHS, M, A, LDA, WORK( ITAU ), B, &
                     LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!*
      ELSE
!*
!*        Path 2 - remaining underdetermined cases.
!*
         IE = 1
         ITAUQ = IE + M
         ITAUP = ITAUQ + M
         NWORK = ITAUP + M
!*
!*        Bidiagonalize A.
!*        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
!*
         CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                     WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, &
                     INFO )
!*
!*        Multiply B by transpose of left bidiagonalizing vectors.
!*        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)
!*
         CALL DORMBR( 'Q', 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAUQ ), &
                     B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!*
!*        Solve the bidiagonal least squares problem.
!*
         CALL DLALSD( 'L', SMLSIZ, M, NRHS, S, WORK( IE ), B, LDB, &
                     RCOND, RANK, WORK( NWORK ), IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF
!*
!*        Multiply B by right bidiagonalizing vectors of A.
!*
         CALL DORMBR( 'P', 'L', 'N', N, NRHS, M, A, LDA, WORK( ITAUP ), &
                     B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
!*
      END IF
!*
!!*     Undo scaling.
!*
      IF( IASCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, &
                     INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, &
                     INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      END IF

   10 CONTINUE
      WORK( 1 ) = MAXWRK
      IWORK( 1 ) = LIWORK
     
END SUBROUTINE 

    
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleOTHERcomputational
!*
!*> \par Contributors:
!*  ==================
!*>
!*>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!*>       California at Berkeley, USA \n
!*>     Osni Marques, LBNL/NERSC, USA \n
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      SUBROUTINE DLALSD( UPLO, SMLSIZ, N, NRHS, D, E, B, LDB, RCOND, &
           RANK, WORK, IWORK, INFO ) !GCC$ ATTRIBUTES hot :: DLALSD !GCC$ ATTRIBUTES hot :: DLALSD
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            BX, BXST, C, DIFL, DIFR, GIVCOL, GIVNUM, &
                        GIVPTR, I, ICMPQ1, ICMPQ2, IWK, J, K, NLVL, &
                        NM1, NSIZE, NSUB, NWORK, PERM, POLES, S, SIZEI, &
                        SMLSZP, SQRE, ST, ST1, U, VT, Z
      DOUBLE PRECISION   CS, EPS, ORGNRM, R, RCND, SN, TOL
!*     ..
!*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           IDAMAX, DLAMCH, DLANST
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLALSA, DLARTG, DLASCL, &
                       DLASDA, DLASDQ, DLASET, DLASRT, DROT
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, SIGN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
!*
      IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.1 ) THEN
         INFO = -4
      ELSE IF( ( LDB.LT.1 ) .OR. ( LDB.LT.N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLALSD', -INFO )
         RETURN
      END IF
!*
      EPS = DLAMCH( 'Epsilon' )
!*
!*     Set up the tolerance.
!*
      IF( ( RCOND.LE.ZERO ) .OR. ( RCOND.GE.ONE ) ) THEN
         RCND = EPS
      ELSE
         RCND = RCOND
      END IF

      RANK = 0
!*
!*     Quick return if possible.
!*
      IF( N.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         IF( D( 1 ).EQ.ZERO ) THEN
            CALL DLASET( 'A', 1, NRHS, ZERO, ZERO, B, LDB )
         ELSE
            RANK = 1
            CALL DLASCL( 'G', 0, 0, D( 1 ), ONE, 1, NRHS, B, LDB, INFO )
            D( 1 ) = ABS( D( 1 ) )
         END IF
         RETURN
      END IF
!*
!*     Rotate the matrix if it is lower bidiagonal.
!*
      IF( UPLO.EQ.'L' ) THEN
         DO 10 I = 1, N - 1
            CALL DLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            IF( NRHS.EQ.1 ) THEN
               CALL DROT( 1, B( I, 1 ), 1, B( I+1, 1 ), 1, CS, SN )
            ELSE
               WORK( I*2-1 ) = CS
               WORK( I*2 ) = SN
            END IF
   10    CONTINUE
         IF( NRHS.GT.1 ) THEN
            DO 30 I = 1, NRHS
               DO 20 J = 1, N - 1
                  CS = WORK( J*2-1 )
                  SN = WORK( J*2 )
                  CALL DROT( 1, B( J, I ), 1, B( J+1, I ), 1, CS, SN )
   20          CONTINUE
   30       CONTINUE
         END IF
      END IF
!*
!*     Scale.
!*
      NM1 = N - 1
      ORGNRM = DLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO ) THEN
         CALL DLASET( 'A', N, NRHS, ZERO, ZERO, B, LDB )
         RETURN
      END IF
!*
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, INFO )
!*
!*     If N is smaller than the minimum divide size SMLSIZ, then solve
!*     the problem with another solver.
!*
      IF( N.LE.SMLSIZ ) THEN
         NWORK = 1 + N*N
         CALL DLASET( 'A', N, N, ZERO, ONE, WORK, N )
         CALL DLASDQ( 'U', 0, N, N, 0, NRHS, D, E, WORK, N, WORK, N, B, &
                     LDB, WORK( NWORK ), INFO )
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
         TOL = RCND*ABS( D( IDAMAX( N, D, 1 ) ) )
         DO 40 I = 1, N
            IF( D( I ).LE.TOL ) THEN
               CALL DLASET( 'A', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB )
            ELSE
               CALL DLASCL( 'G', 0, 0, D( I ), ONE, 1, NRHS, B( I, 1 ), &
                           LDB, INFO )
               RANK = RANK + 1
            END IF
   40    CONTINUE
         CALL DGEMM( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO, &
                    WORK( NWORK ), N )
         CALL DLACPY( 'A', N, NRHS, WORK( NWORK ), N, B, LDB )
!*
!*        Unscale.
!*
         CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
         CALL DLASRT( 'D', N, D, INFO )
         CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO )
!*
         RETURN
      END IF
!*
!*     Book-keeping and setting up some constants.
!*
      NLVL = INT( LOG( DBLE( N ) / DBLE( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1
!*
      SMLSZP = SMLSIZ + 1
!*
      U = 1
      VT = 1 + SMLSIZ*N
      DIFL = VT + SMLSZP*N
      DIFR = DIFL + NLVL*N
      Z = DIFR + NLVL*N*2
      C = Z + NLVL*N
      S = C + N
      POLES = S + N
      GIVNUM = POLES + 2*NLVL*N
      BX = GIVNUM + 2*NLVL*N
      NWORK = BX + N*NRHS
!*
      SIZEI = 1 + N
      K = SIZEI + N
      GIVPTR = K + N
      PERM = GIVPTR + N
      GIVCOL = PERM + NLVL*N
      IWK = GIVCOL + NLVL*N*2
!*
      ST = 1
      SQRE = 0
      ICMPQ1 = 1
      ICMPQ2 = 0
      NSUB = 0
!*
      DO 50 I = 1, N
         IF( ABS( D( I ) ).LT.EPS ) THEN
            D( I ) = SIGN( EPS, D( I ) )
         END IF
   50 CONTINUE
!*
      DO 60 I = 1, NM1
         IF( ( ABS( E( I ) ).LT.EPS ) .OR. ( I.EQ.NM1 ) ) THEN
            NSUB = NSUB + 1
            IWORK( NSUB ) = ST
!*
!*           Subproblem found. First determine its size and then
!*           apply divide and conquer on it.
!*
            IF( I.LT.NM1 ) THEN
!*
!*              A subproblem with E(I) small for I < NM1.
!*
               NSIZE = I - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
            ELSE IF( ABS( E( I ) ).GE.EPS ) THEN
!*
!*              A subproblem with E(NM1) not too small but I = NM1.
!*
               NSIZE = N - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
            ELSE
!*
!*              A subproblem with E(NM1) small. This implies an
!*              1-by-1 subproblem at D(N), which is not solved
!*              explicitly.
!*
               NSIZE = I - ST + 1
               IWORK( SIZEI+NSUB-1 ) = NSIZE
               NSUB = NSUB + 1
               IWORK( NSUB ) = N
               IWORK( SIZEI+NSUB-1 ) = 1
               CALL DCOPY( NRHS, B( N, 1 ), LDB, WORK( BX+NM1 ), N )
            END IF
            ST1 = ST - 1
            IF( NSIZE.EQ.1 ) THEN
!*
!*              This is a 1-by-1 subproblem and is not solved
!*              explicitly.
!*
               CALL DCOPY( NRHS, B( ST, 1 ), LDB, WORK( BX+ST1 ), N )
            ELSE IF( NSIZE.LE.SMLSIZ ) THEN
!*
!*              This is a small subproblem and is solved by DLASDQ.
!*
               CALL DLASET( 'A', NSIZE, NSIZE, ZERO, ONE, &
                           WORK( VT+ST1 ), N )
               CALL DLASDQ( 'U', 0, NSIZE, NSIZE, 0, NRHS, D( ST ),  &
                           E( ST ), WORK( VT+ST1 ), N, WORK( NWORK ), &
                           N, B( ST, 1 ), LDB, WORK( NWORK ), INFO )
               IF( INFO.NE.0 ) THEN
                  RETURN
               END IF
               CALL DLACPY( 'A', NSIZE, NRHS, B( ST, 1 ), LDB, &
                           WORK( BX+ST1 ), N )
            ELSE
!*
!*              A large problem. Solve it using divide and conquer.
!*
               CALL DLASDA( ICMPQ1, SMLSIZ, NSIZE, SQRE, D( ST ),  &
                           E( ST ), WORK( U+ST1 ), N, WORK( VT+ST1 ), &
                           IWORK( K+ST1 ), WORK( DIFL+ST1 ), &
                           WORK( DIFR+ST1 ), WORK( Z+ST1 ), &
                           WORK( POLES+ST1 ), IWORK( GIVPTR+ST1 ), &
                           IWORK( GIVCOL+ST1 ), N, IWORK( PERM+ST1 ), &
                           WORK( GIVNUM+ST1 ), WORK( C+ST1 ), &
                           WORK( S+ST1 ), WORK( NWORK ), IWORK( IWK ), &
                           INFO )
               IF( INFO.NE.0 ) THEN
                  RETURN
               END IF
               BXST = BX + ST1
               CALL DLALSA( ICMPQ2, SMLSIZ, NSIZE, NRHS, B( ST, 1 ), &
                           LDB, WORK( BXST ), N, WORK( U+ST1 ), N, &
                           WORK( VT+ST1 ), IWORK( K+ST1 ), &
                           WORK( DIFL+ST1 ), WORK( DIFR+ST1 ), &
                           WORK( Z+ST1 ), WORK( POLES+ST1 ), &
                           IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, &
                           IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ), &
                           WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ), &
                           IWORK( IWK ), INFO )
               IF( INFO.NE.0 ) THEN
                  RETURN
               END IF
            END IF
            ST = I + 1
         END IF
   60 CONTINUE
!*
!*     Apply the singular values and treat the tiny ones as zero.
!*
      TOL = RCND*ABS( D( IDAMAX( N, D, 1 ) ) )
!*
      DO 70 I = 1, N
!*
!*        Some of the elements in D can be negative because 1-by-1
!*        subproblems were not solved explicitly.
!*
         IF( ABS( D( I ) ).LE.TOL ) THEN
            CALL DLASET( 'A', 1, NRHS, ZERO, ZERO, WORK( BX+I-1 ), N )
         ELSE
            RANK = RANK + 1
            CALL DLASCL( 'G', 0, 0, D( I ), ONE, 1, NRHS, &
                        WORK( BX+I-1 ), N, INFO )
         END IF
         D( I ) = ABS( D( I ) )
   70 CONTINUE
!*
!*     Now apply back the right singular vectors.
!*
      ICMPQ2 = 1
      DO 80 I = 1, NSUB
         ST = IWORK( I )
         ST1 = ST - 1
         NSIZE = IWORK( SIZEI+I-1 )
         BXST = BX + ST1
         IF( NSIZE.EQ.1 ) THEN
            CALL DCOPY( NRHS, WORK( BXST ), N, B( ST, 1 ), LDB )
         ELSE IF( NSIZE.LE.SMLSIZ ) THEN
            CALL DGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, ONE, &
                       WORK( VT+ST1 ), N, WORK( BXST ), N, ZERO, &
                       B( ST, 1 ), LDB )
         ELSE
            CALL DLALSA( ICMPQ2, SMLSIZ, NSIZE, NRHS, WORK( BXST ), N, &
                        B( ST, 1 ), LDB, WORK( U+ST1 ), N, &
                        WORK( VT+ST1 ), IWORK( K+ST1 ), &
                        WORK( DIFL+ST1 ), WORK( DIFR+ST1 ), &
                        WORK( Z+ST1 ), WORK( POLES+ST1 ), &
                        IWORK( GIVPTR+ST1 ), IWORK( GIVCOL+ST1 ), N, &
                        IWORK( PERM+ST1 ), WORK( GIVNUM+ST1 ), &
                        WORK( C+ST1 ), WORK( S+ST1 ), WORK( NWORK ), &
                        IWORK( IWK ), INFO )
            IF( INFO.NE.0 ) THEN
               RETURN
            END IF
         END IF
   80 CONTINUE
!*
!*     Unscale and sort the singular values.
!*
      CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
      CALL DLASRT( 'D', N, D, INFO )
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO )

END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date June 2017
!*
!*> \ingroup doubleOTHERcomputational
!*
!*> \par Contributors:
!*  ==================
!*>
!*>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!*>       California at Berkeley, USA \n
!*>     Osni Marques, LBNL/NERSC, USA \n
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))    
SUBROUTINE DLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, &
                        LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, &
                        GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK, &
                        IWORK, INFO ) !GCC$ ATTRIBUTES hot :: DLALSA !GCC$ ATTRIBUTES aligned(32) :: DLALSA
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
      
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, I1, IC, IM1, INODE, J, LF, LL, LVL, LVL2, &
                        ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL, &
                        NR, NRF, NRP1, SQRE
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM
!*     ..
!!*     .. Executable Statements ..
!*
!*!     Test the input parameters.
!*
      INFO = 0
!*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( SMLSIZ.LT.3 ) THEN
         INFO = -2
      ELSE IF( N.LT.SMLSIZ ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.1 ) THEN
         INFO = -4
      ELSE IF( LDB.LT.N ) THEN
         INFO = -6
      ELSE IF( LDBX.LT.N ) THEN
         INFO = -8
      ELSE IF( LDU.LT.N ) THEN
         INFO = -10
      ELSE IF( LDGCOL.LT.N ) THEN
         INFO = -19
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DLALSA', -INFO )
         RETURN
      END IF
!*
!*!     Book-keeping and  setting up the computation tree.
!*!
      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N

      CALL DLASDT( N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), &
                  IWORK( NDIMR ), SMLSIZ )
!*
!*     The following code applies back the left singular vector factors.
!*     For applying back the right singular vector factors, go to 50.
!*
      IF( ICOMPQ.EQ.1 ) THEN
         GO TO 50
      END IF
!*
!*     The nodes on the bottom level of the tree were solved
!*     by DLASDQ. The corresponding left and right singular vector
!*     matrices are in explicit form. First apply back the left
!*     singular vector matrices.
!*
      NDB1 = ( ND+1 ) / 2
      DO 10 I = NDB1, ND
!*
!*        IC : center row of each node
!*        NL : number of rows of left  subproblem
!*        NR : number of rows of right subproblem
!*        NLF: starting row of the left   subproblem
!*        NRF: starting row of the right  subproblem
!*
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NR = IWORK( NDIMR+I1 )
         NLF = IC - NL
         NRF = IC + 1
         CALL DGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, &
                    B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
         CALL DGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, &
                    B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
   10 CONTINUE
!*
!*     Next copy the rows of B that correspond to unchanged rows
!*     in the bidiagonal matrix to BX.
!*
      DO 20 I = 1, ND
         IC = IWORK( INODE+I-1 )
         CALL DCOPY( NRHS, B( IC, 1 ), LDB, BX( IC, 1 ), LDBX )
   20 CONTINUE
!*
!*     Finally go through the left singular vector matrices of all
!*     the other subproblems bottom-up on the tree.
!*
      J = 2**NLVL
      SQRE = 0
!*
      DO 40 LVL = NLVL, 1, -1
         LVL2 = 2*LVL - 1
!*
!*        find the first node LF and last node LL on
!*        the current level LVL
!*
         IF( LVL.EQ.1 ) THEN
            LF = 1
            LL = 1
         ELSE
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         END IF
         DO 30 I = LF, LL
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            NRF = IC + 1
            J = J - 1
            CALL DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, BX( NLF, 1 ), LDBX, &
                        B( NLF, 1 ), LDB, PERM( NLF, LVL ), &
                        GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, &
                        GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), &
                        DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), &
                        Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK, &
                       INFO )
   30    CONTINUE
   40 CONTINUE
      GO TO 90
!*
!*     ICOMPQ = 1: applying back the right singular vector factors.
!*
   50 CONTINUE
!*
!*     First now go through the right singular vector matrices of all
!!*     the tree nodes top-down.
!*
      J = 0
      DO 70 LVL = 1, NLVL
         LVL2 = 2*LVL - 1
!*
!*        Find the first node LF and last node LL on
!*        the current level LVL.
!*
         IF( LVL.EQ.1 ) THEN
            LF = 1
            LL = 1
         ELSE
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         END IF
         DO 60 I = LL, LF, -1
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            NRF = IC + 1
            IF( I.EQ.LL ) THEN
               SQRE = 0
            ELSE
               SQRE = 1
            END IF
            J = J + 1
            CALL DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B( NLF, 1 ), LDB, &
                        BX( NLF, 1 ), LDBX, PERM( NLF, LVL ), &
                        GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, &
                        GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ), &
                        DIFL( NLF, LVL ), DIFR( NLF, LVL2 ), &
                        Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK, &
                        INFO )
   60    CONTINUE
   70 CONTINUE
!*
!*     The nodes on the bottom level of the tree were solved
!*     by DLASDQ. The corresponding right singular vector
!*     matrices are in explicit form. Apply them back.
!*
      NDB1 = ( ND+1 ) / 2
      DO 80 I = NDB1, ND
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NR = IWORK( NDIMR+I1 )
         NLP1 = NL + 1
         IF( I.EQ.ND ) THEN
            NRP1 = NR
         ELSE
            NRP1 = NR + 1
         END IF
         NLF = IC - NL
         NRF = IC + 1
         CALL DGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, &
                    B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
         CALL DGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, &
                    B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
   80 CONTINUE

   90 CONTINUE

    
END SUBROUTINE

!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleOTHERcomputational
!*
!*> \par Contributors:
!*  ==================
!*>
!*>     Ming Gu and Ren-Cang Li, Computer Science Division, University of
!*>       California at Berkeley, USA \n
!*>     Osni Marques, LBNL/NERSC, USA \n
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))    
      SUBROUTINE DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, &
                        PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, &
                        POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO ) !GCC$ ATTRIBUTES hot :: DLALS0 !GCC$ ATTRIBUTES aligned(32) :: DLALS0
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, NEGONE
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, NEGONE = -1.0D0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, J, M, N, NLP1
      DOUBLE PRECISION   DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, TEMP
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMV,DROT, DSCAL
     
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      N = NL + NR + 1
!*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( NL.LT.1 ) THEN
         INFO = -2
      ELSE IF( NR.LT.1 ) THEN
         INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.1 ) THEN
         INFO = -5
      ELSE IF( LDB.LT.N ) THEN
         INFO = -7
      ELSE IF( LDBX.LT.N ) THEN
         INFO = -9
      ELSE IF( GIVPTR.LT.0 ) THEN
         INFO = -11
      ELSE IF( LDGCOL.LT.N ) THEN
         INFO = -13
      ELSE IF( LDGNUM.LT.N ) THEN
         INFO = -15
      ELSE IF( K.LT.1 ) THEN
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
        ! CALL XERBLA( 'DLALS0', -INFO )
         RETURN
      END IF
!*
      M = N + SQRE
      NLP1 = NL + 1
!*
      IF( ICOMPQ.EQ.0 ) THEN
!*
!*        Apply back orthogonal transformations from the left.
!*
!*        Step (1L): apply back the Givens rotations performed.
!1*
         DO 10 I = 1, GIVPTR
            CALL DROT( NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, &
                      B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), &
                      GIVNUM( I, 1 ) )
   10    CONTINUE
!*
!*        Step (2L): permute rows of B.
!*
         CALL DCOPY( NRHS, B( NLP1, 1 ), LDB, BX( 1, 1 ), LDBX )
         DO 20 I = 2, N
            CALL DCOPY( NRHS, B( PERM( I ), 1 ), LDB, BX( I, 1 ), LDBX )
   20    CONTINUE
!*
!*        Step (3L): apply the inverse of the left singular vector
!*        matrix to BX.
!*
         IF( K.EQ.1 ) THEN
            CALL DCOPY( NRHS, BX, LDBX, B, LDB )
            IF( Z( 1 ).LT.ZERO ) THEN
               CALL DSCAL( NRHS, NEGONE, B, LDB )
            END IF
         ELSE
            !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(DIFL,POLES,DIFR,WORK,Z,K) PRIVATE(J,DIFLJ,DJ,DSIGJ,DIFRJ,DSIGJP,I,TEMP)
            DO 50 J = 1, K
               DIFLJ = DIFL( J )
               DJ = POLES( J, 1 )
               DSIGJ = -POLES( J, 2 )
               IF( J.LT.K ) THEN
                  DIFRJ = -DIFR( J, 1 )
                  DSIGJP = -POLES( J+1, 2 )
               END IF
               IF( ( Z( J ).EQ.ZERO ) .OR. ( POLES( J, 2 ).EQ.ZERO ) ) &
                    THEN
                  WORK( J ) = ZERO
               ELSE
                  WORK( J ) = -POLES( J, 2 )*Z( J ) / DIFLJ / &
                             ( POLES( J, 2 )+DJ )
               END IF
               !$OMP SIMD ALIGNED(Z:64,POLES,WORK) LINEAR(I:1)
               DO 30 I = 1, J - 1
                  IF( ( Z( I ).EQ.ZERO ) .OR. &
                     ( POLES( I, 2 ).EQ.ZERO ) ) THEN
                     WORK( I ) = ZERO
                  ELSE
                     WORK( I ) = POLES( I, 2 )*Z( I ) / &
                                ( DLAMC3( POLES( I, 2 ), DSIGJ )- &
                                DIFLJ ) / ( POLES( I, 2 )+DJ )
                  END IF
    30         CONTINUE
                 !$OMP SIMD ALIGNED(Z:64,POLES,WORK) LINEAR(I:1)  
               DO 40 I = J + 1, K
                  IF( ( Z( I ).EQ.ZERO ) .OR. &
                      ( POLES( I, 2 ).EQ.ZERO ) ) THEN
                     WORK( I ) = ZERO
                  ELSE
                     WORK( I ) = POLES( I, 2 )*Z( I ) / &
                                ( DLAMC3( POLES( I, 2 ), DSIGJP )+ &
                                DIFRJ ) / ( POLES( I, 2 )+DJ )
                  END IF
   40          CONTINUE
               WORK( 1 ) = NEGONE
               TEMP = DNRM2( K, WORK, 1 )
               !$OMP SINGLE
               CALL DGEMV( 'T', K, NRHS, ONE, BX, LDBX, WORK, 1, ZERO, &
                    B( J, 1 ), LDB )
               !$OMP END SINGLE NOWAIT
               !$OMP SINGLE
               CALL DLASCL( 'G', 0, 0, TEMP, ONE, 1, NRHS, B( J, 1 ), &
                    LDB, INFO )
               !$OMP END SINGLE NOWAIT
    50      CONTINUE
            !$OMP END PARALLEL DO   
         END IF
!*
!*!        Move the deflated rows of BX to B also.
!*
         IF( K.LT.MAX( M, N ) ) &
           CALL DLACPY( 'A', N-K, NRHS, BX( K+1, 1 ), LDBX, &
                        B( K+1, 1 ), LDB )
      ELSE
!*
!*        Apply back the right orthogonal transformations.
!*
!*!        Step (1R): apply back the new right singular vector matrix
!*        to B.
!*
         IF( K.EQ.1 ) THEN
            CALL DCOPY( NRHS, B, LDB, BX, LDBX )
         ELSE
            !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(POLES,Z,WORK,DIFR,DIFL,K) PRIVATE(J,DSIGJ,I)
            DO 80 J = 1, K
               DSIGJ = POLES( J, 2 )
               IF( Z( J ).EQ.ZERO ) THEN
                  WORK( J ) = ZERO
               ELSE
                  WORK( J ) = -Z( J ) / DIFL( J ) / &
                             ( DSIGJ+POLES( J, 1 ) ) / DIFR( J, 2 )
               END IF
               !$OMP SIMD ALIGNED(Z:64,WORK,POLES,DIFR) LINEAR(I:1)
               DO 60 I = 1, J - 1
                  IF( Z( J ).EQ.ZERO ) THEN
                     WORK( I ) = ZERO
                  ELSE
                     WORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I+1, &
                                2 ) )-DIFR( I, 1 ) ) / &
                                ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  END IF
60             CONTINUE
                  !$OMP SIMD ALIGNED(Z:64,WORK,POLES,DIFR,DIFL) LINEAR(I:1) 
               DO 70 I = J + 1, K
                  IF( Z( J ).EQ.ZERO ) THEN
                     WORK( I ) = ZERO
                  ELSE
                     WORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I, &
                                2 ) )-DIFL( I ) ) / &
                                ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  END IF
70             CONTINUE
                !$OMP SINGLE  
               CALL DGEMV( 'T', K, NRHS, ONE, B, LDB, WORK, 1, ZERO, &
                    BX( J, 1 ), LDBX )
               !$OMP END SINGLE NOWAIT
   80       CONTINUE
         END IF
!*
!*        Step (2R): if SQRE = 1, apply back the rotation that is
!*        related to the right null space of the subproblem.
!*
         IF( SQRE.EQ.1 ) THEN
            CALL DCOPY( NRHS, B( M, 1 ), LDB, BX( M, 1 ), LDBX )
            CALL DROT( NRHS, BX( 1, 1 ), LDBX, BX( M, 1 ), LDBX, C, S )
         END IF
         IF( K.LT.MAX( M, N ) ) &
           CALL DLACPY( 'A', N-K, NRHS, B( K+1, 1 ), LDB, BX( K+1, 1 ), &
                        LDBX )
!*
!*        Step (3R): permute rows of B.
!*
         CALL DCOPY( NRHS, BX( 1, 1 ), LDBX, B( NLP1, 1 ), LDB )
         IF( SQRE.EQ.1 ) THEN
            CALL DCOPY( NRHS, BX( M, 1 ), LDBX, B( M, 1 ), LDB )
         END IF
         DO 90 I = 2, N
            CALL DCOPY( NRHS, BX( I, 1 ), LDBX, B( PERM( I ), 1 ), LDB )
   90    CONTINUE
!*
!*        Step (4R): apply back the Givens rotations performed.
!*
         DO 100 I = GIVPTR, 1, -1
            CALL DROT( NRHS, B( GIVCOL( I, 2 ), 1 ), LDB, &
                      B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ), &
                      -GIVNUM( I, 1 ) )
  100    CONTINUE
      END IF

END SUBROUTINE

!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date June 2017
!*
!*> \ingroup OTHERauxiliary
!*
!*> \par Contributors:
!*  ==================
!*>
!*>     Ming Gu and Huan Ren, Computer Science Division, University of
!*>     California at Berkeley, USA
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))   
      SUBROUTINE DLASDA( ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K, &
                        DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, &
                        PERM, GIVNUM, C, S, WORK, IWORK, INFO ) !GCC$ ATTRIBUTES hot :: DLASDA !GCC$ ATTRIBUTES aligned(32) :: DLASDA
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, I1, IC, IDXQ, IDXQI, IM1, INODE, ITEMP, IWK, &
                       J, LF, LL, LVL, LVL2, M, NCC, ND, NDB1, NDIML, &
                        NDIMR, NL, NLF, NLP1, NLVL, NR, NRF, NRP1, NRU, &
                        NWORK1, NWORK2, SMLSZP, SQREI, VF, VFI, VL, VLI 
      DOUBLE PRECISION   ALPHA, BETA
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DCOPY
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
!*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( SMLSIZ.LT.3 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
         INFO = -4
      ELSE IF( LDU.LT.( N+SQRE ) ) THEN
         INFO = -8
      ELSE IF( LDGCOL.LT.N ) THEN
         INFO = -17
      END IF
      IF( INFO.NE.0 ) THEN
        ! CALL XERBLA( 'DLASDA', -INFO )
         RETURN
      END IF
!*
      M = N + SQRE
!*
!*     If the input matrix is too small, call DLASDQ to find the SVD.
!*
      IF( N.LE.SMLSIZ ) THEN
         IF( ICOMPQ.EQ.0 ) THEN
            CALL DLASDQ( 'U', SQRE, N, 0, 0, 0, D, E, VT, LDU, U, LDU, &
                        U, LDU, WORK, INFO )
         ELSE
            CALL DLASDQ( 'U', SQRE, N, M, N, 0, D, E, VT, LDU, U, LDU, &
                        U, LDU, WORK, INFO )
         END IF
         RETURN
      END IF
!*
!*     Book-keeping and  set up the computation tree.
!*
      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N
      IDXQ = NDIMR + N
      IWK = IDXQ + N
!*
      NCC = 0
      NRU = 0
!*
      SMLSZP = SMLSIZ + 1
      VF = 1
      VL = VF + M
      NWORK1 = VL + M
      NWORK2 = NWORK1 + SMLSZP*SMLSZP
!*
      CALL DLASDT( N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ), &
                  IWORK( NDIMR ), SMLSIZ )
!*
!*     for the nodes on bottom level of the tree, solve
!*     their subproblems by DLASDQ.
!*
      NDB1 = ( ND+1 ) / 2
      DO 30 I = NDB1, ND
!*
!*        IC : center row of each node
!*        NL : number of rows of left  subproblem
!*        NR : number of rows of right subproblem
!*        NLF: starting row of the left   subproblem
!*        NRF: starting row of the right  subproblem
!*
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NLP1 = NL + 1
         NR = IWORK( NDIMR+I1 )
         NLF = IC - NL
         NRF = IC + 1
         IDXQI = IDXQ + NLF - 2
         VFI = VF + NLF - 1
         VLI = VL + NLF - 1
         SQREI = 1
         IF( ICOMPQ.EQ.0 ) THEN
            CALL DLASET( 'A', NLP1, NLP1, ZERO, ONE, WORK( NWORK1 ), &
                        SMLSZP )
            CALL DLASDQ( 'U', SQREI, NL, NLP1, NRU, NCC, D( NLF ), &
                        E( NLF ), WORK( NWORK1 ), SMLSZP, &
                        WORK( NWORK2 ), NL, WORK( NWORK2 ), NL, &
                        WORK( NWORK2 ), INFO )
            ITEMP = NWORK1 + NL*SMLSZP
            CALL DCOPY( NLP1, WORK( NWORK1 ), 1, WORK( VFI ), 1 )
            CALL DCOPY( NLP1, WORK( ITEMP ), 1, WORK( VLI ), 1 )
         ELSE
            CALL DLASET( 'A', NL, NL, ZERO, ONE, U( NLF, 1 ), LDU )
            CALL DLASET( 'A', NLP1, NLP1, ZERO, ONE, VT( NLF, 1 ), LDU )
            CALL DLASDQ( 'U', SQREI, NL, NLP1, NL, NCC, D( NLF ), &
                        E( NLF ), VT( NLF, 1 ), LDU, U( NLF, 1 ), LDU, &
                        U( NLF, 1 ), LDU, WORK( NWORK1 ), INFO )
            CALL DCOPY( NLP1, VT( NLF, 1 ), 1, WORK( VFI ), 1 )
            CALL DCOPY( NLP1, VT( NLF, NLP1 ), 1, WORK( VLI ), 1 )
         END IF
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
         DO 10 J = 1, NL
            IWORK( IDXQI+J ) = J
   10    CONTINUE
         IF( ( I.EQ.ND ) .AND. ( SQRE.EQ.0 ) ) THEN
            SQREI = 0
         ELSE
            SQREI = 1
         END IF
         IDXQI = IDXQI + NLP1
         VFI = VFI + NLP1
         VLI = VLI + NLP1
         NRP1 = NR + SQREI
         IF( ICOMPQ.EQ.0 ) THEN
            CALL DLASET( 'A', NRP1, NRP1, ZERO, ONE, WORK( NWORK1 ), &
                        SMLSZP )
            CALL DLASDQ( 'U', SQREI, NR, NRP1, NRU, NCC, D( NRF ), &
                        E( NRF ), WORK( NWORK1 ), SMLSZP, &
                        WORK( NWORK2 ), NR, WORK( NWORK2 ), NR, &
                        WORK( NWORK2 ), INFO )
            ITEMP = NWORK1 + ( NRP1-1 )*SMLSZP
            CALL DCOPY( NRP1, WORK( NWORK1 ), 1, WORK( VFI ), 1 )
            CALL DCOPY( NRP1, WORK( ITEMP ), 1, WORK( VLI ), 1 )
         ELSE
            CALL DLASET( 'A', NR, NR, ZERO, ONE, U( NRF, 1 ), LDU )
            CALL DLASET( 'A', NRP1, NRP1, ZERO, ONE, VT( NRF, 1 ), LDU )
            CALL DLASDQ( 'U', SQREI, NR, NRP1, NR, NCC, D( NRF ), &
                        E( NRF ), VT( NRF, 1 ), LDU, U( NRF, 1 ), LDU, &
                       U( NRF, 1 ), LDU, WORK( NWORK1 ), INFO )
            CALL DCOPY( NRP1, VT( NRF, 1 ), 1, WORK( VFI ), 1 )
            CALL DCOPY( NRP1, VT( NRF, NRP1 ), 1, WORK( VLI ), 1 )
         END IF
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
         DO 20 J = 1, NR
            IWORK( IDXQI+J ) = J
   20    CONTINUE
   30 CONTINUE
!*
!*     Now conquer each subproblem bottom-up.
!*
      J = 2**NLVL
      DO 50 LVL = NLVL, 1, -1
         LVL2 = LVL*2 - 1
!*
!*        Find the first node LF and last node LL on
!*        the current level LVL.
!*
         IF( LVL.EQ.1 ) THEN
            LF = 1
            LL = 1
         ELSE
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         END IF
         DO 40 I = LF, LL
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            NRF = IC + 1
            IF( I.EQ.LL ) THEN
               SQREI = SQRE
            ELSE
               SQREI = 1
            END IF
            VFI = VF + NLF - 1
            VLI = VL + NLF - 1
            IDXQI = IDXQ + NLF - 1
            ALPHA = D( IC )
            BETA = E( IC )
            IF( ICOMPQ.EQ.0 ) THEN
               CALL DLASD6( ICOMPQ, NL, NR, SQREI, D( NLF ), &
                           WORK( VFI ), WORK( VLI ), ALPHA, BETA, &
                           IWORK( IDXQI ), PERM, GIVPTR( 1 ), GIVCOL, &
                           LDGCOL, GIVNUM, LDU, POLES, DIFL, DIFR, Z, &
                           K( 1 ), C( 1 ), S( 1 ), WORK( NWORK1 ), &
                           IWORK( IWK ), INFO )
            ELSE
               J = J - 1
               CALL DLASD6( ICOMPQ, NL, NR, SQREI, D( NLF ), &
                           WORK( VFI ), WORK( VLI ), ALPHA, BETA, &
                           IWORK( IDXQI ), PERM( NLF, LVL ), &
                           GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL, &
                           GIVNUM( NLF, LVL2 ), LDU, &
                           POLES( NLF, LVL2 ), DIFL( NLF, LVL ), &
                           DIFR( NLF, LVL2 ), Z( NLF, LVL ), K( J ), &
                           C( J ), S( J ), WORK( NWORK1 ), &
                           IWORK( IWK ), INFO )
            END IF
            IF( INFO.NE.0 ) THEN
               RETURN
            END IF
   40    CONTINUE
   50 CONTINUE

END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))   
      SUBROUTINE DLASD6( ICOMPQ, NL, NR, SQRE, D, VF, VL, ALPHA, BETA, &
                        IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, &
                        LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, WORK, &
                        IWORK, INFO ) !GCC$ ATTRIBUTES hot :: DLASD6 !GCC$ ATTRIBUTES aligned(32) :: DLASD6
#elif defined(__INTEL_COMPILER) || defined(__ICC)
 SUBROUTINE DLASD6( ICOMPQ, NL, NR, SQRE, D, VF, VL, ALPHA, BETA, &
                        IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, &
                        LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, WORK, &
                        IWORK, INFO )
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASD6
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLASD6
#endif
       implicit none
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL, &
                         NR, SQRE
      DOUBLE PRECISION   ALPHA, BETA, C, S
!*     ..
!*     .. Array Arguments ..
     ! INTEGER            GIVCOL( LDGCOL, * ), IDXQ( * ), IWORK( * ),
     !$                   PERM( * )
     ! DOUBLE PRECISION   D( * ), DIFL( * ), DIFR( * ),
     !$                   GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ),
      !$                   VF( * ), VL( * ), WORK( * ), Z( * )
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: GIVCOL
      INTEGER, DIMENSION(:), ALLOCATABLE :: IDXQ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      INTEGER, DIMENSION(:), ALLOCATABLE :: PERM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DIFL
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DIFR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GIVNUM
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: POLES
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VF
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VL
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, IDX, IDXC, IDXP, ISIGMA, IVFW, IVLW, IW, M, &
                        N, N1, N2
      DOUBLE PRECISION   ORGNRM
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DCOPY
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!*     ..
!*     .. Executable Statements ..
!1*!
!*     Test the input parameters.
!*
      INFO = 0
      N = NL + NR + 1
      M = N + SQRE
!*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( NL.LT.1 ) THEN
         INFO = -2
      ELSE IF( NR.LT.1 ) THEN
         INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
         INFO = -4
      ELSE IF( LDGCOL.LT.N ) THEN
         INFO = -14
      ELSE IF( LDGNUM.LT.N ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
        ! CALL XERBLA( 'DLASD6', -INFO )
         RETURN
      END IF
!*
!*     The following values are for bookkeeping purposes only.  They are
!*     integer pointers which indicate the portion of the workspace
!*     used by a particular array in DLASD7 and DLASD8.
!*
      ISIGMA = 1
      IW = ISIGMA + N
      IVFW = IW + M
      IVLW = IVFW + M
!*
      IDX = 1
      IDXC = IDX + N
      IDXP = IDXC + N
!*
!*     Scale.
!*
      ORGNRM = MAX( ABS( ALPHA ), ABS( BETA ) )
      D( NL+1 ) = ZERO
      DO 10 I = 1, N
         IF( ABS( D( I ) ).GT.ORGNRM ) THEN
            ORGNRM = ABS( D( I ) )
         END IF
   10 CONTINUE
      CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
      ALPHA = ALPHA / ORGNRM
      BETA = BETA / ORGNRM
!*
!*     Sort and Deflate singular values.
!*
      CALL DLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, WORK( IW ), VF, &
                  WORK( IVFW ), VL, WORK( IVLW ), ALPHA, BETA, &
                  WORK( ISIGMA ), IWORK( IDX ), IWORK( IDXP ), IDXQ, &
                  PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, C, S, &
                  INFO )
!*
!*     Solve Secular Equation, compute DIFL, DIFR, and update VF, VL.
!*
      CALL DLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDGNUM, &
                  WORK( ISIGMA ), WORK( IW ), INFO )
!*
!*     Report the possible convergence failure.
!*
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
!*
!*     Save the poles if ICOMPQ = 1.
!*
      IF( ICOMPQ.EQ.1 ) THEN
         CALL DCOPY( K, D, 1, POLES( 1, 1 ), 1 )
         CALL DCOPY( K, WORK( ISIGMA ), 1, POLES( 1, 2 ), 1 )
      END IF
!*
!*     Unscale.
!*
      CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
!*
!*     Prepare the IDXQ sorting permutation.
!*
      N1 = K
      N2 = N - K
      CALL DLAMRG( N1, N2, D, 1, -1, IDXQ )

    END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))       
SUBROUTINE DLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR, &
     DSIGMA, WORK, INFO ) !GCC$ ATTRIBUTES hot :: DLASD8 !GCC$ ATTRIBUTES aligned(32) :: DLASD8
#elif defined(__INTEL_COMPILER) || defined(__ICC)
SUBROUTINE DLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR, &
     DSIGMA, WORK, INFO )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASD8
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLASD8
#endif
  use omp_lib
      implicit none
!*
!*  -- LAPACK auxiliary routine (version 3.7.1) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2017
!*
!*     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, K, LDDIFR
!*     ..
!*     .. Array Arguments ..
     ! DOUBLE PRECISION   D( * ), DIFL( * ), DIFR( LDDIFR, * ),
     !$                   DSIGMA( * ), VF( * ), VL( * ), WORK( * ),
      ! $                   Z( * )
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DIFL
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DIFR
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DSIGMA
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VF
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VL
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z
      
!*     ..
!!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, IWK1, IWK2, IWK2I, IWK3, IWK3I, J
      DOUBLE PRECISION   DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, RHO, TEMP
!*     ..
!*     .. External Subroutines ..
     EXTERNAL           DCOPY
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DDOT, DLAMC3, DNRM2
      EXTERNAL           DDOT, DLAMC3, DNRM2
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
!*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( K.LT.1 ) THEN
         INFO = -2
      ELSE IF( LDDIFR.LT.K ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DLASD8', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      IF( K.EQ.1 ) THEN
         D( 1 ) = ABS( Z( 1 ) )
         DIFL( 1 ) = D( 1 )
         IF( ICOMPQ.EQ.1 ) THEN
            DIFL( 2 ) = ONE
            DIFR( 1, 2 ) = ONE
         END IF
         RETURN
      END IF
#if 0
*
*     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
*     be computed with high relative accuracy (barring over/underflow).
*     This is a problem on machines without a guard digit in
*     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
*     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
*     which on any of these machines zeros out the bottommost
*     bit of DSIGMA(I) if it is 1; this makes the subsequent
*     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
*     occurs. On binary machines with a guard digit (almost all
*     machines) it does not change DSIGMA(I) at all. On hexadecimal
*     and decimal machines with a guard digit, it slightly
*     changes the bottommost bits of DSIGMA(I). It does not account
*     for hexadecimal or decimal machines without guard digits
*     (we know of none). We use a subroutine call to compute
*     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
*     this code.
*
#endif
      !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(DSIGMA,K) PRIVATE(I)
      DO 10 I = 1, K
         DSIGMA( I ) = DLAMC3( DSIGMA( I ), DSIGMA( I ) ) - DSIGMA( I )
   10 CONTINUE
!*
!*     Book keeping.
!*
      IWK1 = 1
      IWK2 = IWK1 + K
      IWK3 = IWK2 + K
      IWK2I = IWK2 - 1
      IWK3I = IWK3 - 1
!*
!*     Normalize Z.
!*
      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
!*
!*     Initialize WORK(IWK3).
!*
      CALL DLASET( 'A', K, 1, ONE, ONE, WORK( IWK3 ), K )
!*
!*     Compute the updated singular values, the arrays DIFL, DIFR,
!*     and the updated Z.
!*
      DO 40 J = 1, K
         CALL DLASD4( K, J, DSIGMA, Z, WORK( IWK1 ), RHO, D( J ), &
                     WORK( IWK2 ), INFO )
!*
!*        If the root finder fails, report the convergence failure.
!*
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
         WORK( IWK3I+J ) = WORK( IWK3I+J )*WORK( J )*WORK( IWK2I+J )
         DIFL( J ) = -WORK( J )
         DIFR( J, 1 ) = -WORK( J+1 )
         DO 20 I = 1, J - 1
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* &
                             WORK( IWK2I+I ) / ( DSIGMA( I )- &
                             DSIGMA( J ) ) / ( DSIGMA( I )+ &
                             DSIGMA( J ) )
   20    CONTINUE
         DO 30 I = J + 1, K
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* &
                             WORK( IWK2I+I ) / ( DSIGMA( I )- &
                             DSIGMA( J ) ) / ( DSIGMA( I )+ &
                             DSIGMA( J ) )
   30    CONTINUE
   40 CONTINUE
!*
!*     Compute updated Z.
            !*
     !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(Z,WORK,K) PRIVATE(I)
      DO 50 I = 1, K
         Z( I ) = SIGN( SQRT( ABS( WORK( IWK3I+I ) ) ), Z( I ) )
   50 CONTINUE
!*
!*     Update VF and VL.
         !*
      !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(DILF,D,DSIGMA,DIFR,WORK,Z,K,IWK2I,IWK3I) PRIVATE(J,DIFLJ,DJ,DSIGJ,DIFRJ,DSIGJP,I,TEMP)
      DO 80 J = 1, K
         DIFLJ = DIFL( J )
         DJ = D( J )
         DSIGJ = -DSIGMA( J )
         IF( J.LT.K ) THEN
            DIFRJ = -DIFR( J, 1 )
            DSIGJP = -DSIGMA( J+1 )
         END IF
         WORK( J ) = -Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )
         !$OMP SIMD ALIGNED(WORK:64,Z,DSIGMA)
         DO 60 I = 1, J - 1
            WORK( I ) = Z( I ) / ( DLAMC3( DSIGMA( I ), DSIGJ )-DIFLJ ) &
                         / ( DSIGMA( I )+DJ )
60       CONTINUE
             !$OMP SIMD ALIGNED(WORK:64,Z,DSIGMA)
         DO 70 I = J + 1, K
            WORK( I ) = Z( I ) / ( DLAMC3( DSIGMA( I ), DSIGJP )+DIFRJ ) &
                        / ( DSIGMA( I )+DJ )
   70    CONTINUE
         TEMP = DNRM2( K, WORK, 1 )
         WORK( IWK2I+J ) = DDOT( K, WORK, 1, VF, 1 ) / TEMP
         WORK( IWK3I+J ) = DDOT( K, WORK, 1, VL, 1 ) / TEMP
         IF( ICOMPQ.EQ.1 ) THEN
            DIFR( J, 2 ) = TEMP
         END IF
   80 CONTINUE
!$OMP END PARALLEL DO
      CALL DCOPY( K, WORK( IWK2 ), 1, VF, 1 )
      CALL DCOPY( K, WORK( IWK3 ), 1, VL, 1 )

END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))   
SUBROUTINE DLASD4( N, I, D, Z, DELTA, RHO, SIGMA, WORK, INFO ) !GCC$ ATTRIBUTES Hot :: DLASD4 !GCC$ ATTRIBUTES aligned(32) :: DLASD4
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DLASD4( N, I, D, Z, DELTA, RHO, SIGMA, WORK, INFO )
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASD4
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLASD4
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
      INTEGER            I, INFO, N
      DOUBLE PRECISION   RHO, SIGMA
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   D( * ), DELTA( * ), WORK( * ), Z( * )
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D, DELTA, WORK, Z
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 400 )
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT, TEN
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0, &
                        THREE = 3.0D+0, FOUR = 4.0D+0, EIGHT = 8.0D+0, &
                        TEN = 10.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            ORGATI, SWTCH, SWTCH3, GEOMAVG
      INTEGER            II, IIM1, IIP1, IP1, ITER, J, NITER
      DOUBLE PRECISION   A, B, C, DELSQ, DELSQ2, SQ2, DPHI, DPSI, DTIIM, &
                        DTIIP, DTIPSQ, DTISQ, DTNSQ, DTNSQ1, DW, EPS, &
                        ERRETM, ETA, PHI, PREW, PSI, RHOINV, SGLB, &
                        SGUB, TAU, TAU2, TEMP, TEMP1, TEMP2, W
!*     ..
!*     .. Local Arrays ..
      DOUBLE PRECISION   DD( 3 ), ZZ( 3 )
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DLAED6, DLASD5
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!*     ..
!*     .. Executable Statements ..
!*
!*     Since this routine is called in an inner loop, we do no argument
!*     checking.
!*
!*     Quick return for N=1 and 2.
!*
      INFO = 0
      IF( N.EQ.1 ) THEN
!*
!*        Presumably, I=1 upon entry
!*
         SIGMA = SQRT( D( 1 )*D( 1 )+RHO*Z( 1 )*Z( 1 ) )
         DELTA( 1 ) = ONE
         WORK( 1 ) = ONE
         RETURN
      END IF
      IF( N.EQ.2 ) THEN
         CALL DLASD5( I, D, Z, DELTA, RHO, SIGMA, WORK )
         RETURN
      END IF
!*
!*     Compute machine epsilon
!*
      EPS = DLAMCH( 'Epsilon' )
      RHOINV = ONE / RHO
      TAU2= ZERO
!*
!!*     The case I = N
!*
      IF( I.EQ.N ) THEN
!*
!*        Initialize some basic variables
!*
         II = N - 1
         NITER = 1
!*
!*        Calculate initial guess
!*
         TEMP = RHO / TWO
!*
!*        If ||Z||_2 is not one, then TEMP should be set to
!*        RHO * ||Z||_2^2 / TWO
!*
         TEMP1 = TEMP / ( D( N )+SQRT( D( N )*D( N )+TEMP ) )
         !$OMP PARALLEL DO SCHEDULE(STATIC,2) DEFAULT(NONE) SHARED(WORK,DELTA,D,TEMP1,N) PRIVATE(J)
         DO 10 J = 1, N
            WORK( J ) = D( J ) + D( N ) + TEMP1
            DELTA( J ) = ( D( J )-D( N ) ) - TEMP1
   10    CONTINUE
!*
            PSI = ZERO
         !OMP PARALLEL DO REDUCTION(+:PSI) SHARED(Z,DELTA,WORK,N) FIRSTPRIVATE(PSI) PRIVATE(J)   
         DO 20 J = 1, N - 2
            PSI = PSI + Z( J )*Z( J ) / ( DELTA( J )*WORK( J ) )
   20    CONTINUE

         C = RHOINV + PSI
         W = C + Z( II )*Z( II ) / ( DELTA( II )*WORK( II ) ) + &
            Z( N )*Z( N ) / ( DELTA( N )*WORK( N ) )

         IF( W.LE.ZERO ) THEN
            TEMP1 = SQRT( D( N )*D( N )+RHO )
            TEMP = Z( N-1 )*Z( N-1 ) / ( ( D( N-1 )+TEMP1 )* &
                  ( D( N )-D( N-1 )+RHO / ( D( N )+TEMP1 ) ) ) + &
                  Z( N )*Z( N ) / RHO
!*
!*           The following TAU2 is to approximate
!*           SIGMA_n^2 - D( N )*D( N )
!*
            IF( C.LE.TEMP ) THEN
               TAU = RHO
            ELSE
               DELSQ = ( D( N )-D( N-1 ) )*( D( N )+D( N-1 ) )
               A = -C*DELSQ + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
               B = Z( N )*Z( N )*DELSQ
               IF( A.LT.ZERO ) THEN
                  TAU2 = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
               ELSE
                  TAU2 = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
               END IF
               TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) )
            END IF
!*
!*           It can be proved that
!*               D(N)^2+RHO/2 <= SIGMA_n^2 < D(N)^2+TAU2 <= D(N)^2+RHO
!*
         ELSE
            DELSQ = ( D( N )-D( N-1 ) )*( D( N )+D( N-1 ) )
            A = -C*DELSQ + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DELSQ
!*
!*           The following TAU2 is to approximate
!*           SIGMA_n^2 - D( N )*D( N )
!*
            IF( A.LT.ZERO ) THEN
               TAU2 = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            ELSE
               TAU2 = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            END IF
            TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) )

!*
!*           It can be proved that
!*           D(N)^2 < D(N)^2+TAU2 < SIGMA(N)^2 < D(N)^2+RHO/2
!*
         END IF
!*
!*        The following TAU is to approximate SIGMA_n - D( N )
!*
!*         TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) )
!*
         SIGMA = D( N ) + TAU
          !$OMP PARALLEL DO SCHEDULE(STATIC,10) DEFAULT(NONE) SHARED(WORK,DELTA,D,N,TAU) PRIVATE(J)
         DO 30 J = 1, N
            DELTA( J ) = ( D( J )-D( N ) ) - TAU
            WORK( J ) = D( J ) + D( N ) + TAU
   30    CONTINUE
!*
!*        Evaluate PSI and the derivative DPSI
!*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 40 J = 1, II
            TEMP = Z( J ) / ( DELTA( J )*WORK( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   40    CONTINUE
         ERRETM = ABS( ERRETM )
!*
!*        Evaluate PHI and the derivative DPHI
!*
         TEMP = Z( N ) / ( DELTA( N )*WORK( N ) )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV 
!*    $          + ABS( TAU2 )*( DPSI+DPHI )
!*
         W = RHOINV + PHI + PSI
!*
!*        Test for convergence
!*
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            GO TO 240
         END IF
!*
!*        Calculate the new step
!*
         NITER = NITER + 1
         DTNSQ1 = WORK( N-1 )*DELTA( N-1 )
         DTNSQ = WORK( N )*DELTA( N )
         C = W - DTNSQ1*DPSI - DTNSQ*DPHI
         A = ( DTNSQ+DTNSQ1 )*W - DTNSQ*DTNSQ1*( DPSI+DPHI )
         B = DTNSQ*DTNSQ1*W
         IF( C.LT.ZERO ) &
           C = ABS( C )
         IF( C.EQ.ZERO ) THEN
            ETA = RHO - SIGMA*SIGMA
         ELSE IF( A.GE.ZERO ) THEN
            ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
!*
!*        Note, eta should be positive if w is negative, and
!*        eta should be negative otherwise. However,
!*        if for some reason caused by roundoff, eta*w > 0,
!*        we simply use one Newton step instead. This way
!*        will guarantee eta*w < 0.
!*
         IF( W*ETA.GT.ZERO ) &
           ETA = -W / ( DPSI+DPHI )
         TEMP = ETA - DTNSQ
         IF( TEMP.GT.RHO ) &
           ETA = RHO + DTNSQ

         ETA = ETA / ( SIGMA+SQRT( ETA+SIGMA*SIGMA ) )
         TAU = TAU + ETA
         SIGMA = SIGMA + ETA

         !$OMP PARALLEL DO SCHEDULE(STATIC,10) DEFAULT(NONE) SHARED(WORK,DELTA,ETA,N) PRIVATE(J)
         DO 50 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
            WORK( J ) = WORK( J ) + ETA
   50    CONTINUE
!*
!*        Evaluate PSI and the derivative DPSI
!*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 60 J = 1, II
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   60    CONTINUE
         ERRETM = ABS( ERRETM )
!*
!*        Evaluate PHI and the derivative DPHI
!*
         TAU2 = WORK( N )*DELTA( N )
         TEMP = Z( N ) / TAU2
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV
!*    $          + ABS( TAU2 )*( DPSI+DPHI )
!*
         W = RHOINV + PHI + PSI
!*
!*        Main loop to update the values of the array   DELTA
!*
         ITER = NITER + 1
!*
         DO 90 NITER = ITER, MAXIT
!*
!*           Test for convergence
!*
            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               GO TO 240
            END IF
!*
!*           Calculate the new step
!*
            DTNSQ1 = WORK( N-1 )*DELTA( N-1 )
            DTNSQ = WORK( N )*DELTA( N )
            C = W - DTNSQ1*DPSI - DTNSQ*DPHI
            A = ( DTNSQ+DTNSQ1 )*W - DTNSQ1*DTNSQ*( DPSI+DPHI )
            B = DTNSQ1*DTNSQ*W
            IF( A.GE.ZERO ) THEN
               ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
!*
!*           Note, eta should be positive if w is negative, and
!*           eta should be negative otherwise. However,
!*           if for some reason caused by roundoff, eta*w > 0,
!*           we simply use one Newton step instead. This way
!*           will guarantee eta*w < 0.
!*
            IF( W*ETA.GT.ZERO ) &
              ETA = -W / ( DPSI+DPHI )
            TEMP = ETA - DTNSQ
            IF( TEMP.LE.ZERO ) &
               ETA = ETA / TWO
!*
            ETA = ETA / ( SIGMA+SQRT( ETA+SIGMA*SIGMA ) )
            TAU = TAU + ETA
            SIGMA = SIGMA + ETA
            !*
             !$OMP PARALLEL DO SCHEDULE(STATIC,10) DEFAULT(NONE) SHARED(WORK,DELTA,N,ETA) PRIVATE(J)
            DO 70 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
               WORK( J ) = WORK( J ) + ETA
   70       CONTINUE
!*
!*           Evaluate PSI and the derivative DPSI
!*
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 80 J = 1, II
               TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
   80       CONTINUE
            ERRETM = ABS( ERRETM )
!*
!*           Evaluate PHI and the derivative DPHI
!*
            TAU2 = WORK( N )*DELTA( N )
            TEMP = Z( N ) / TAU2
            PHI = Z( N )*TEMP
            DPHI = TEMP*TEMP
            ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV
!*!    $             + ABS( TAU2 )*( DPSI+DPHI )
!*
            W = RHOINV + PHI + PSI
   90    CONTINUE
!*
!*        Return with INFO = 1, NITER = MAXIT and not converged
!*
         INFO = 1
         GO TO 240
!*
!*        End for the case I = N
!*
      ELSE
!*
!*        The case for I < N
!*
         NITER = 1
         IP1 = I + 1
!*
!*        Calculate initial guess
!*
         DELSQ = ( D( IP1 )-D( I ) )*( D( IP1 )+D( I ) )
         DELSQ2 = DELSQ / TWO
         SQ2=SQRT( ( D( I )*D( I )+D( IP1 )*D( IP1 ) ) / TWO )
         TEMP = DELSQ2 / ( D( I )+SQ2 )
          !$OMP PARALLEL DO SCHEDULE(STATIC,10) DEFAULT(NONE) SHARED(WORK,DELTA,D,N,TEMP) PRIVATE(J)
         DO 100 J = 1, N
            WORK( J ) = D( J ) + D( I ) + TEMP
            DELTA( J ) = ( D( J )-D( I ) ) - TEMP
  100    CONTINUE

            PSI = ZERO
             !OMP PARALLEL DO REDUCTION(+:PSI) SHARED(Z,DELTA,WORK,I) FIRSTPRIVATE(PSI) PRIVATE(J)   
         DO 110 J = 1, I - 1
            PSI = PSI + Z( J )*Z( J ) / ( WORK( J )*DELTA( J ) )
  110    CONTINUE

            PHI = ZERO
             !OMP PARALLEL DO REDUCTION(+:PHI) SHARED(Z,DELTA,WORK) PRIVATE(J)   
         DO 120 J = N, I + 2, -1
            PHI = PHI + Z( J )*Z( J ) / ( WORK( J )*DELTA( J ) )
  120    CONTINUE
         C = RHOINV + PSI + PHI
         W = C + Z( I )*Z( I ) / ( WORK( I )*DELTA( I ) ) + &
            Z( IP1 )*Z( IP1 ) / ( WORK( IP1 )*DELTA( IP1 ) )

         GEOMAVG = .FALSE.
         IF( W.GT.ZERO ) THEN
!*
!*           d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2
!*
!*           We choose d(i) as origin.
!*
            ORGATI = .TRUE.
            II = I
            SGLB = ZERO
            SGUB = DELSQ2  / ( D( I )+SQ2 )
            A = C*DELSQ + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
            B = Z( I )*Z( I )*DELSQ
            IF( A.GT.ZERO ) THEN
               TAU2 = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            ELSE
               TAU2 = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            END IF
!*
!*           TAU2 now is an estimation of SIGMA^2 - D( I )^2. The
!*           following, however, is the corresponding estimation of
!*           SIGMA - D( I ).
!*
            TAU = TAU2 / ( D( I )+SQRT( D( I )*D( I )+TAU2 ) )
            TEMP = SQRT(EPS)
            IF( (D(I).LE.TEMP*D(IP1)).AND.(ABS(Z(I)).LE.TEMP) &
                                    .AND.(D(I).GT.ZERO) ) THEN
               TAU = MIN( TEN*D(I), SGUB )
               GEOMAVG = .TRUE.
            END IF
         ELSE
!*
!*           (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2
!*
!*           We choose d(i+1) as origin.
!*
            ORGATI = .FALSE.
            II = IP1
            SGLB = -DELSQ2  / ( D( II )+SQ2 )
            SGUB = ZERO
            A = C*DELSQ - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
            B = Z( IP1 )*Z( IP1 )*DELSQ
            IF( A.LT.ZERO ) THEN
               TAU2 = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
            ELSE
               TAU2 = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
            END IF
!*
!*           TAU2 now is an estimation of SIGMA^2 - D( IP1 )^2. The
!*           following, however, is the corresponding estimation of
!*           SIGMA - D( IP1 ).
!*
            TAU = TAU2 / ( D( IP1 )+SQRT( ABS( D( IP1 )*D( IP1 )+ &
                 TAU2 ) ) )
         END IF

         SIGMA = D( II ) + TAU
           !$OMP PARALLEL DO SCHEDULE(STATIC,10) DEFAULT(NONE) SHARED(WORK,DELTA,D,N,TAU) PRIVATE(J)
         DO 130 J = 1, N
            WORK( J ) = D( J ) + D( II ) + TAU
            DELTA( J ) = ( D( J )-D( II ) ) - TAU
  130    CONTINUE
         IIM1 = II - 1
         IIP1 = II + 1
!*
!*        Evaluate PSI and the derivative DPSI
!*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 150 J = 1, IIM1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  150    CONTINUE
         ERRETM = ABS( ERRETM )
!*
!*        Evaluate PHI and the derivative DPHI
!*
         DPHI = ZERO
         PHI = ZERO
         DO 160 J = N, IIP1, -1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  160    CONTINUE
!*
         W = RHOINV + PHI + PSI
!*
!*        W is the value of the secular function with
!*        its ii-th element removed.
!*
         SWTCH3 = .FALSE.
         IF( ORGATI ) THEN
            IF( W.LT.ZERO ) &
              SWTCH3 = .TRUE.
         ELSE
            IF( W.GT.ZERO ) &
              SWTCH3 = .TRUE.
         END IF
         IF( II.EQ.1 .OR. II.EQ.N ) &
           SWTCH3 = .FALSE.

         TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = W + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV &
               + THREE*ABS( TEMP ) &
!*              + ABS( TAU2 )*DW
!*
!*        Test for convergence
!*
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            GO TO 240
         END IF
!*
         IF( W.LE.ZERO ) THEN
            SGLB = MAX( SGLB, TAU )
         ELSE
            SGUB = MIN( SGUB, TAU )
         END IF
!*
!*        Calculate the new step
!*
         NITER = NITER + 1
         IF( .NOT.SWTCH3 ) THEN
            DTIPSQ = WORK( IP1 )*DELTA( IP1 )
            DTISQ = WORK( I )*DELTA( I )
            IF( ORGATI ) THEN
               C = W - DTIPSQ*DW + DELSQ*( Z( I ) / DTISQ )**2
            ELSE
               C = W - DTISQ*DW - DELSQ*( Z( IP1 ) / DTIPSQ )**2
            END IF
            A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
            B = DTIPSQ*DTISQ*W
            IF( C.EQ.ZERO ) THEN
               IF( A.EQ.ZERO ) THEN
                  IF( ORGATI ) THEN
                     A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ*( DPSI+DPHI )
                  ELSE
                     A = Z( IP1 )*Z( IP1 ) + DTISQ*DTISQ*( DPSI+DPHI )
                  END IF
               END IF
               ETA = B / A
            ELSE IF( A.LE.ZERO ) THEN
               ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
         ELSE
!*
!*           Interpolation using THREE most relevant poles
!*
            DTIIM = WORK( IIM1 )*DELTA( IIM1 )
            DTIIP = WORK( IIP1 )*DELTA( IIP1 )
            TEMP = RHOINV + PSI + PHI
            IF( ORGATI ) THEN
               TEMP1 = Z( IIM1 ) / DTIIM
               TEMP1 = TEMP1*TEMP1
               C = ( TEMP - DTIIP*( DPSI+DPHI ) ) - &
                  ( D( IIM1 )-D( IIP1 ) )*( D( IIM1 )+D( IIP1 ) )*TEMP1
               ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
               IF( DPSI.LT.TEMP1 ) THEN
                  ZZ( 3 ) = DTIIP*DTIIP*DPHI
               ELSE
                  ZZ( 3 ) = DTIIP*DTIIP*( ( DPSI-TEMP1 )+DPHI )
               END IF
            ELSE
               TEMP1 = Z( IIP1 ) / DTIIP
               TEMP1 = TEMP1*TEMP1
               C = ( TEMP - DTIIM*( DPSI+DPHI ) ) - &
                  ( D( IIP1 )-D( IIM1 ) )*( D( IIM1 )+D( IIP1 ) )*TEMP1
               IF( DPHI.LT.TEMP1 ) THEN
                  ZZ( 1 ) = DTIIM*DTIIM*DPSI
               ELSE
                  ZZ( 1 ) = DTIIM*DTIIM*( DPSI+( DPHI-TEMP1 ) )
               END IF
               ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
            END IF
            ZZ( 2 ) = Z( II )*Z( II )
            DD( 1 ) = DTIIM
            DD( 2 ) = DELTA( II )*WORK( II )
            DD( 3 ) = DTIIP
            CALL DLAED6( NITER, ORGATI, C, DD, ZZ, W, ETA, INFO )

            IF( INFO.NE.0 ) THEN
!*
!*              If INFO is not 0, i.e., DLAED6 failed, switch back
!*              to 2 pole interpolation.
!*
               SWTCH3 = .FALSE.
               INFO = 0
               DTIPSQ = WORK( IP1 )*DELTA( IP1 )
               DTISQ = WORK( I )*DELTA( I )
               IF( ORGATI ) THEN
                  C = W - DTIPSQ*DW + DELSQ*( Z( I ) / DTISQ )**2
               ELSE
                  C = W - DTISQ*DW - DELSQ*( Z( IP1 ) / DTIPSQ )**2
               END IF
               A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
               B = DTIPSQ*DTISQ*W
               IF( C.EQ.ZERO ) THEN
                  IF( A.EQ.ZERO ) THEN
                     IF( ORGATI ) THEN
                        A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ*( DPSI+DPHI )
                     ELSE
                        A = Z( IP1 )*Z( IP1 ) + DTISQ*DTISQ*( DPSI+DPHI)
                     END IF
                  END IF
                  ETA = B / A
               ELSE IF( A.LE.ZERO ) THEN
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               ELSE
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               END IF
            END IF
         END IF
!*
!*        Note, eta should be positive if w is negative, and
!*        eta should be negative otherwise. However,
!*        if for some reason caused by roundoff, eta*w > 0,
!*        we simply use one Newton step instead. This way
!*        will guarantee eta*w < 0.
!*
         IF( W*ETA.GE.ZERO ) &
           ETA = -W / DW

         ETA = ETA / ( SIGMA+SQRT( SIGMA*SIGMA+ETA ) )
         TEMP = TAU + ETA
         IF( TEMP.GT.SGUB .OR. TEMP.LT.SGLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( SGUB-TAU ) / TWO
            ELSE
               ETA = ( SGLB-TAU ) / TWO
            END IF
            IF( GEOMAVG ) THEN
               IF( W .LT. ZERO ) THEN
                  IF( TAU .GT. ZERO ) THEN
                     ETA = SQRT(SGUB*TAU)-TAU
                  END IF
               ELSE
                  IF( SGLB .GT. ZERO ) THEN
                     ETA = SQRT(SGLB*TAU)-TAU
                  END IF
               END IF
            END IF
         END IF

         PREW = W

         TAU = TAU + ETA
         SIGMA = SIGMA + ETA
       !$OMP PARALLEL DO SCHEDULE(STATIC,10) DEFAULT(NONE) SHARED(WORK,DELTA,N,ETA) PRIVATE(J)
         DO 170 J = 1, N
            WORK( J ) = WORK( J ) + ETA
            DELTA( J ) = DELTA( J ) - ETA
  170    CONTINUE
!*
!*        Evaluate PSI and the derivative DPSI
!*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 180 J = 1, IIM1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  180    CONTINUE
         ERRETM = ABS( ERRETM )
!*
!*        Evaluate PHI and the derivative DPHI
!*
         DPHI = ZERO
         PHI = ZERO
         DO 190 J = N, IIP1, -1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  190    CONTINUE
!*
         TAU2 = WORK( II )*DELTA( II )
         TEMP = Z( II ) / TAU2
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = RHOINV + PHI + PSI + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV &
               + THREE*ABS( TEMP )
!*    $          + ABS( TAU2 )*DW
!*
         SWTCH = .FALSE.
         IF( ORGATI ) THEN
            IF( -W.GT.ABS( PREW ) / TEN ) &
              SWTCH = .TRUE.
         ELSE
            IF( W.GT.ABS( PREW ) / TEN ) &
              SWTCH = .TRUE.
         END IF
!*
!*        Main loop to update the values of the array   DELTA and WORK
!*
         ITER = NITER + 1
!*
         DO 230 NITER = ITER, MAXIT
!*
!*           Test for convergence
!*
            IF( ABS( W ).LE.EPS*ERRETM ) THEN 
!*     $          .OR. (SGUB-SGLB).LE.EIGHT*ABS(SGUB+SGLB) ) THEN
               GO TO 240
            END IF
!*
            IF( W.LE.ZERO ) THEN
               SGLB = MAX( SGLB, TAU )
            ELSE
               SGUB = MIN( SGUB, TAU )
            END IF
!*
!*           Calculate the new step
!*
            IF( .NOT.SWTCH3 ) THEN
               DTIPSQ = WORK( IP1 )*DELTA( IP1 )
               DTISQ = WORK( I )*DELTA( I )
               IF( .NOT.SWTCH ) THEN
                  IF( ORGATI ) THEN
                     C = W - DTIPSQ*DW + DELSQ*( Z( I ) / DTISQ )**2
                  ELSE
                     C = W - DTISQ*DW - DELSQ*( Z( IP1 ) / DTIPSQ )**2
                  END IF
               ELSE
                  TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
                  IF( ORGATI ) THEN
                     DPSI = DPSI + TEMP*TEMP
                  ELSE
                     DPHI = DPHI + TEMP*TEMP
                  END IF
                  C = W - DTISQ*DPSI - DTIPSQ*DPHI
               END IF
               A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
               B = DTIPSQ*DTISQ*W
               IF( C.EQ.ZERO ) THEN
                  IF( A.EQ.ZERO ) THEN
                     IF( .NOT.SWTCH ) THEN
                        IF( ORGATI ) THEN &
                           A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ*
                              ( DPSI+DPHI )
                        ELSE
                           A = Z( IP1 )*Z( IP1 ) + &
                              DTISQ*DTISQ*( DPSI+DPHI )
                        END IF
                     ELSE
                        A = DTISQ*DTISQ*DPSI + DTIPSQ*DTIPSQ*DPHI
                     END IF
                  END IF
                  ETA = B / A
               ELSE IF( A.LE.ZERO ) THEN
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               ELSE
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               END IF
            ELSE
!*
!*              Interpolation using THREE most relevant poles
!*
               DTIIM = WORK( IIM1 )*DELTA( IIM1 )
               DTIIP = WORK( IIP1 )*DELTA( IIP1 )
               TEMP = RHOINV + PSI + PHI
               IF( SWTCH ) THEN
                  C = TEMP - DTIIM*DPSI - DTIIP*DPHI
                  ZZ( 1 ) = DTIIM*DTIIM*DPSI
                  ZZ( 3 ) = DTIIP*DTIIP*DPHI
               ELSE
                  IF( ORGATI ) THEN
                     TEMP1 = Z( IIM1 ) / DTIIM
                     TEMP1 = TEMP1*TEMP1
                     TEMP2 = ( D( IIM1 )-D( IIP1 ) )* &
                            ( D( IIM1 )+D( IIP1 ) )*TEMP1
                     C = TEMP - DTIIP*( DPSI+DPHI ) - TEMP2
                     ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                     IF( DPSI.LT.TEMP1 ) THEN
                        ZZ( 3 ) = DTIIP*DTIIP*DPHI
                     ELSE
                        ZZ( 3 ) = DTIIP*DTIIP*( ( DPSI-TEMP1 )+DPHI )
                     END IF
                  ELSE
                     TEMP1 = Z( IIP1 ) / DTIIP
                     TEMP1 = TEMP1*TEMP1
                     TEMP2 = ( D( IIP1 )-D( IIM1 ) )* &
                            ( D( IIM1 )+D( IIP1 ) )*TEMP1
                     C = TEMP - DTIIM*( DPSI+DPHI ) - TEMP2
                     IF( DPHI.LT.TEMP1 ) THEN
                        ZZ( 1 ) = DTIIM*DTIIM*DPSI
                     ELSE
                        ZZ( 1 ) = DTIIM*DTIIM*( DPSI+( DPHI-TEMP1 ) )
                     END IF
                     ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
                  END IF
               END IF
               DD( 1 ) = DTIIM
               DD( 2 ) = DELTA( II )*WORK( II )
               DD( 3 ) = DTIIP
               CALL DLAED6( NITER, ORGATI, C, DD, ZZ, W, ETA, INFO )

               IF( INFO.NE.0 ) THEN
!*
!*                 If INFO is not 0, i.e., DLAED6 failed, switch
!*                 back to two pole interpolation
!*
                  SWTCH3 = .FALSE.
                  INFO = 0
                  DTIPSQ = WORK( IP1 )*DELTA( IP1 )
                  DTISQ = WORK( I )*DELTA( I )
                  IF( .NOT.SWTCH ) THEN
                     IF( ORGATI ) THEN
                        C = W - DTIPSQ*DW + DELSQ*( Z( I )/DTISQ )**2
                     ELSE
                        C = W - DTISQ*DW - DELSQ*( Z( IP1 )/DTIPSQ )**2
                     END IF
                  ELSE
                     TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
                     IF( ORGATI ) THEN
                        DPSI = DPSI + TEMP*TEMP
                     ELSE
                        DPHI = DPHI + TEMP*TEMP
                     END IF
                     C = W - DTISQ*DPSI - DTIPSQ*DPHI
                  END IF
                  A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
                  B = DTIPSQ*DTISQ*W
                  IF( C.EQ.ZERO ) THEN
                     IF( A.EQ.ZERO ) THEN
                        IF( .NOT.SWTCH ) THEN
                           IF( ORGATI ) THEN
                              A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ* &
                                 ( DPSI+DPHI )
                           ELSE
                              A = Z( IP1 )*Z( IP1 ) + &
                                 DTISQ*DTISQ*( DPSI+DPHI )
                           END IF
                        ELSE
                           A = DTISQ*DTISQ*DPSI + DTIPSQ*DTIPSQ*DPHI
                        END IF
                     END IF
                     ETA = B / A
                  ELSE IF( A.LE.ZERO ) THEN
                     ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
                  ELSE
                     ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
                  END IF
               END IF
            END IF
#if 0
*
*           Note, eta should be positive if w is negative, and
*           eta should be negative otherwise. However,
*           if for some reason caused by roundoff, eta*w > 0,
*           we simply use one Newton step instead. This way
*           will guarantee eta*w < 0.
*
#endif
            IF( W*ETA.GE.ZERO ) &
              ETA = -W / DW

            ETA = ETA / ( SIGMA+SQRT( SIGMA*SIGMA+ETA ) )
            TEMP=TAU+ETA
            IF( TEMP.GT.SGUB .OR. TEMP.LT.SGLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( SGUB-TAU ) / TWO
               ELSE
                  ETA = ( SGLB-TAU ) / TWO
               END IF
               IF( GEOMAVG ) THEN
                  IF( W .LT. ZERO ) THEN
                     IF( TAU .GT. ZERO ) THEN
                        ETA = SQRT(SGUB*TAU)-TAU
                     END IF
                  ELSE
                     IF( SGLB .GT. ZERO ) THEN
                        ETA = SQRT(SGLB*TAU)-TAU
                     END IF
                  END IF
               END IF
            END IF

            PREW = W

            TAU = TAU + ETA
            SIGMA = SIGMA + ETA
            !$OMP PARALLEL DO SCHEDULE(STATIC,10) DEFAULT(NONE) SHARED(WORK,DELTA,N,ETA) PRIVATE(J)
            DO 200 J = 1, N
               WORK( J ) = WORK( J ) + ETA
               DELTA( J ) = DELTA( J ) - ETA
  200       CONTINUE
!*
!*           Evaluate PSI and the derivative DPSI
!*
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 210 J = 1, IIM1
               TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
  210       CONTINUE
            ERRETM = ABS( ERRETM )
!*
!*           Evaluate PHI and the derivative DPHI
!*
            DPHI = ZERO
            PHI = ZERO
            DO 220 J = N, IIP1, -1
               TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
               PHI = PHI + Z( J )*TEMP
               DPHI = DPHI + TEMP*TEMP
               ERRETM = ERRETM + PHI
  220       CONTINUE

            TAU2 = WORK( II )*DELTA( II )
            TEMP = Z( II ) / TAU2
            DW = DPSI + DPHI + TEMP*TEMP
            TEMP = Z( II )*TEMP
            W = RHOINV + PHI + PSI + TEMP
            ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV &
                  + THREE*ABS( TEMP )
!*    $             + ABS( TAU2 )*DW
!*
            IF( W*PREW.GT.ZERO .AND. ABS( W ).GT.ABS( PREW ) / TEN ) &
              SWTCH = .NOT.SWTCH

  230    CONTINUE
!*
!*        Return with INFO = 1, NITER = MAXIT and not converged
!*
         INFO = 1
!*
      END IF
!*
  240 CONTINUE
     
END SUBROUTINE


#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))   
      SUBROUTINE DLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, ZW, VF, VFW, VL, &
                        VLW, ALPHA, BETA, DSIGMA, IDX, IDXP, IDXQ, &
                       PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, &
                       C, S, INFO ) !GCC$ ATTRIBUTES hot :: DLASD7 !GCC$ ATTRIBUTES aligned(32) :: DLASD7
#elif defined(__INTEL_COMPILER) || defined(__ICC)
 SUBROUTINE DLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, ZW, VF, VFW, VL, &
                        VLW, ALPHA, BETA, DSIGMA, IDX, IDXP, IDXQ, &
                       PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, &
                       C, S, INFO )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASD7
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLASD7
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
      INTEGER            GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL, &
                        NR, SQRE
      DOUBLE PRECISION   ALPHA, BETA, C, S
!*     ..
!*     .. Array Arguments ..
     ! INTEGER            GIVCOL( LDGCOL, * ), IDX( * ), IDXP( * ),
    ! $                   IDXQ( * ), PERM( * )
    !  DOUBLE PRECISION   D( * ), DSIGMA( * ), GIVNUM( LDGNUM, * ),
    ! $                   VF( * ), VFW( * ), VL( * ), VLW( * ), Z( * ),
      ! $                   ZW( * )
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: GIVCOL
      INTEGER, DIMENSION(:),   ALLOCATABLE :: IDX
      INTEGER, DIMENSION(:),   ALLOCATABLE :: IDXP
      INTEGER, DIMENSION(:),   ALLOCATABLE :: IDXQ
      INTEGER, DIMENSION(:),   ALLOCATABLE :: PERM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DSIGMA
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GIVNUM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VF
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VFW
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VL
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VLW
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZW
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, EIGHT
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0, &
                        EIGHT = 8.0D+0 )
!*     ..
!*     .. Local Scalars ..
!*
      INTEGER            I, IDXI, IDXJ, IDXJP, J, JP, JPREV, K2, M, N, &
                        NLP1, NLP2
      DOUBLE PRECISION   EPS, HLFTOL, TAU, TOL, Z1
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DCOPY, DROT
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      N = NL + NR + 1
      M = N + SQRE
!*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( NL.LT.1 ) THEN
         INFO = -2
      ELSE IF( NR.LT.1 ) THEN
         INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
         INFO = -4
      ELSE IF( LDGCOL.LT.N ) THEN
         INFO = -22
      ELSE IF( LDGNUM.LT.N ) THEN
         INFO = -24
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASD7', -INFO )
         RETURN
      END IF
!*
      NLP1 = NL + 1
      NLP2 = NL + 2
      IF( ICOMPQ.EQ.1 ) THEN
         GIVPTR = 0
      END IF
!*
!*     Generate the first part of the vector Z and move the singular
!*     values in the first part of D one position backward.
!*
      Z1 = ALPHA*VL( NLP1 )
      VL( NLP1 ) = ZERO
      TAU = VF( NLP1 )
      !$OMP SIMD ALIGNED(Z:64,VL,VF,D,IDXQ)
      DO 10 I = NL, 1, -1
         Z( I+1 ) = ALPHA*VL( I )
         VL( I ) = ZERO
         VF( I+1 ) = VF( I )
         D( I+1 ) = D( I )
         IDXQ( I+1 ) = IDXQ( I ) + 1
   10 CONTINUE
      VF( 1 ) = TAU
!*
!*     Generate the second part of the vector Z.
      !*
      !$OMP SIMD ALIGNED(Z:64,VF)
      DO 20 I = NLP2, M
         Z( I ) = BETA*VF( I )
         VF( I ) = ZERO
   20 CONTINUE
!*
!*     Sort the singular values into increasing order
         !*
        !$OMP SIMD ALIGNED(IDXQ:64)   
      DO 30 I = NLP2, N
         IDXQ( I ) = IDXQ( I ) + NLP1
   30 CONTINUE
!*
!*     DSIGMA, IDXC, IDXC, and ZW are used as storage space.
         !*
     
      DO 40 I = 2, N
         DSIGMA( I ) = D( IDXQ( I ) )
         ZW( I ) = Z( IDXQ( I ) )
         VFW( I ) = VF( IDXQ( I ) )
         VLW( I ) = VL( IDXQ( I ) )
   40 CONTINUE

      CALL DLAMRG( NL, NR, DSIGMA( 2 ), 1, 1, IDX( 2 ) )

      DO 50 I = 2, N
         IDXI = 1 + IDX( I )
         D( I ) = DSIGMA( IDXI )
         Z( I ) = ZW( IDXI )
         VF( I ) = VFW( IDXI )
         VL( I ) = VLW( IDXI )
   50 CONTINUE
!*
!*     Calculate the allowable deflation tolerance
!*
      EPS = DLAMCH( 'Epsilon' )
      TOL = MAX( ABS( ALPHA ), ABS( BETA ) )
      TOL = EIGHT*EIGHT*EPS*MAX( ABS( D( N ) ), TOL )
#if 0
*
*     There are 2 kinds of deflation -- first a value in the z-vector
*     is small, second two (or more) singular values are very close
*     together (their difference is small).
*
*     If the value in the z-vector is small, we simply permute the
*     array so that the corresponding singular value is moved to the
*     end.
*
*     If two values in the D-vector are close, we perform a two-sided
*     rotation designed to make one of the corresponding z-vector
*     entries zero, and then permute the array so that the deflated
*     singular value is moved to the end.
*
*     If there are multiple singular values then the problem deflates.
*     Here the number of equal singular values are found.  As each equal
*     singular value is found, an elementary reflector is computed to
*     rotate the corresponding singular subspace so that the
*     corresponding components of Z are zero in this new basis.
*
#endif
      K = 1
      K2 = N + 1
      DO 60 J = 2, N
         IF( ABS( Z( J ) ).LE.TOL ) THEN
!*
!*           Deflate due to small z component.
!*
            K2 = K2 - 1
            IDXP( K2 ) = J
            IF( J.EQ.N ) &
               GO TO 100
         ELSE
            JPREV = J
            GO TO 70
         END IF
   60 CONTINUE
   70 CONTINUE
      J = JPREV
   80 CONTINUE
      J = J + 1
      IF( J.GT.N ) &
        GO TO 90
      IF( ABS( Z( J ) ).LE.TOL ) THEN
!*
!*        Deflate due to small z component.
!*
         K2 = K2 - 1
         IDXP( K2 ) = J
      ELSE
!*
!*        Check if singular values are close enough to allow deflation.
!*
         IF( ABS( D( J )-D( JPREV ) ).LE.TOL ) THEN
!*
!*           Deflation is possible.
!*
            S = Z( JPREV )
            C = Z( J )
!*
!!*           Find sqrt(a**2+b**2) without overflow or
!*           destructive underflow.
!*
            TAU = DLAPY2( C, S )
            Z( J ) = TAU
            Z( JPREV ) = ZERO
            C = C / TAU
            S = -S / TAU
!*
!*           Record the appropriate Givens rotation
!*
            IF( ICOMPQ.EQ.1 ) THEN
               GIVPTR = GIVPTR + 1
               IDXJP = IDXQ( IDX( JPREV )+1 )
               IDXJ = IDXQ( IDX( J )+1 )
               IF( IDXJP.LE.NLP1 ) THEN
                  IDXJP = IDXJP - 1
               END IF
               IF( IDXJ.LE.NLP1 ) THEN
                  IDXJ = IDXJ - 1
               END IF
               GIVCOL( GIVPTR, 2 ) = IDXJP
               GIVCOL( GIVPTR, 1 ) = IDXJ
               GIVNUM( GIVPTR, 2 ) = C
               GIVNUM( GIVPTR, 1 ) = S
            END IF
            CALL DROT( 1, VF( JPREV ), 1, VF( J ), 1, C, S )
            CALL DROT( 1, VL( JPREV ), 1, VL( J ), 1, C, S )
            K2 = K2 - 1
            IDXP( K2 ) = JPREV
            JPREV = J
         ELSE
            K = K + 1
            ZW( K ) = Z( JPREV )
            DSIGMA( K ) = D( JPREV )
            IDXP( K ) = JPREV
            JPREV = J
         END IF
      END IF
      GO TO 80
   90 CONTINUE
!*
!!*     Record the last singular value.
!*
      K = K + 1
      ZW( K ) = Z( JPREV )
      DSIGMA( K ) = D( JPREV )
      IDXP( K ) = JPREV
!*
  100 CONTINUE
!*
!*     Sort the singular values into DSIGMA. The singular values which
!*     were not deflated go into the first K slots of DSIGMA, except
!*     that DSIGMA(1) is treated separately.
      !*
      
      DO 110 J = 2, N
         JP = IDXP( J )
         DSIGMA( J ) = D( JP )
         VFW( J ) = VF( JP )
         VLW( J ) = VL( JP )
  110 CONTINUE
      IF( ICOMPQ.EQ.1 ) THEN
         DO 120 J = 2, N
            JP = IDXP( J )
            PERM( J ) = IDXQ( IDX( JP )+1 )
            IF( PERM( J ).LE.NLP1 ) THEN
               PERM( J ) = PERM( J ) - 1
            END IF
  120    CONTINUE
      END IF
!*
!*     The deflated singular values go back into the last N - K slots of
!*     D.
!*
      CALL DCOPY( N-K, DSIGMA( K+1 ), 1, D( K+1 ), 1 )
!*
!*     Determine DSIGMA(1), DSIGMA(2), Z(1), VF(1), VL(1), VF(M), and
!*     VL(M).
!*
      DSIGMA( 1 ) = ZERO
      HLFTOL = TOL / TWO
      IF( ABS( DSIGMA( 2 ) ).LE.HLFTOL ) &
        DSIGMA( 2 ) = HLFTOL
      IF( M.GT.N ) THEN
         Z( 1 ) = DLAPY2( Z1, Z( M ) )
         IF( Z( 1 ).LE.TOL ) THEN
            C = ONE
            S = ZERO
            Z( 1 ) = TOL
         ELSE
            C = Z1 / Z( 1 )
            S = -Z( M ) / Z( 1 )
         END IF
         CALL DROT( 1, VF( M ), 1, VF( 1 ), 1, C, S )
         CALL DROT( 1, VL( M ), 1, VL( 1 ), 1, C, S )
      ELSE
         IF( ABS( Z1 ).LE.TOL ) THEN
            Z( 1 ) = TOL
         ELSE
            Z( 1 ) = Z1
         END IF
      END IF
!*
!*     Restore Z, VF, and VL.
!*
      CALL DCOPY( K-1, ZW( 2 ), 1, Z( 2 ), 1 )
      CALL DCOPY( N-1, VFW( 2 ), 1, VF( 2 ), 1 )
      CALL DCOPY( N-1, VLW( 2 ), 1, VL( 2 ), 1 )

END SUBROUTINE

    
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER)) 
SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX ) !GCC$ ATTRIBUTES inline :: DLAMRG !GCC$ ATTRIBUTES aligned(32) :: DLAMRG
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER            I, IND1, IND2, N1SV, N2SV
!*     ..
!*     .. Executable Statements ..
!*
      N1SV = N1
      N2SV = N2
      IF( DTRD1.GT.0 ) THEN
         IND1 = 1
      ELSE
         IND1 = N1
      END IF
      IF( DTRD2.GT.0 ) THEN
         IND2 = 1 + N1
      ELSE
         IND2 = N1 + N2
      END IF
      I = 1
!*     while ( (N1SV > 0) & (N2SV > 0) )
   10 CONTINUE
      IF( N1SV.GT.0 .AND. N2SV.GT.0 ) THEN
         IF( A( IND1 ).LE.A( IND2 ) ) THEN
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + DTRD1
            N1SV = N1SV - 1
         ELSE
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + DTRD2
            N2SV = N2SV - 1
         END IF
         GO TO 10
      END IF
!*     end while
      IF( N1SV.EQ.0 ) THEN
         DO 20 N1SV = 1, N2SV
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + DTRD2
   20    CONTINUE
      ELSE
!*     N2SV .EQ. 0
         DO 30 N2SV = 1, N1SV
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + DTRD1
   30    CONTINUE
      END IF

END SUBROUTINE



#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))     
SUBROUTINE DLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB ) !GCC$ ATTRIBUTES inline :: DLASDT
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, IL, IR, LLST, MAXN, NCRNT, NLVL
      DOUBLE PRECISION   TEMP
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, LOG, MAX
!*     ..
!*     .. Executable Statements ..
!*
!*     Find the number of levels on the tree.
!*
      MAXN = MAX( 1, N )
      TEMP = LOG( DBLE( MAXN ) / DBLE( MSUB+1 ) ) / LOG( TWO )
      LVL = INT( TEMP ) + 1
!*
      I = N / 2
      INODE( 1 ) = I + 1
      NDIML( 1 ) = I
      NDIMR( 1 ) = N - I - 1
      IL = 0
      IR = 1
      LLST = 1
      DO 20 NLVL = 1, LVL - 1
!*
!*        Constructing the tree at (NLVL+1)-st level. The number of
!*        nodes created on this level is LLST * 2.
!*
         DO 10 I = 0, LLST - 1
            IL = IL + 2
            IR = IR + 2
            NCRNT = LLST + I
            NDIML( IL ) = NDIML( NCRNT ) / 2
            NDIMR( IL ) = NDIML( NCRNT ) - NDIML( IL ) - 1
            INODE( IL ) = INODE( NCRNT ) - NDIMR( IL ) - 1
            NDIML( IR ) = NDIMR( NCRNT ) / 2
            NDIMR( IR ) = NDIMR( NCRNT ) - NDIML( IR ) - 1
            INODE( IR ) = INODE( NCRNT ) + NDIML( IR ) + 1
   10    CONTINUE
         LLST = LLST*2
   20 CONTINUE
      ND = LLST*2 - 1
    END SUBROUTINE


    

    !*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))    
    SUBROUTINE DLABAD( SMALL, LARGE ) !GCC$ ATTRIBUTES inline :: DLABAD !GCC$ ATTRIBUTES aligned(32) :: DLABAD
#elif defined(__INTEL_COMPILER) || defined(__ICC)
      SUBROUTINE DLABAD( SMALL, LARGE )
     !DIR$ ATTRIBUTES FORCEINLINE :: DLABAD
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLABAD
    !DIR$ OPTIMIZE : 3
   
#endif
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
!*     ..
!*
!!*  =====================================================================
!*
!*     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
!*     ..
!*     .. Executable Statements ..
!*
!*     If it looks like we're on a Cray, take the square root of
!*     SMALL and LARGE to avoid overflow and underflow problems.
!*
      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
!*
    
!*
!*     End of DLABAD
!*
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))    
SUBROUTINE DRSCL( N, SA, SX, INCX ) !GCC$ ATTRIBUTES inline :: DRSCL !GCC$ ATTRIBUTES aligned(32) :: DRSCL
#elif defined(__INTEL_COMPILER) || defined(__ICC)
SUBROUTINE DRSCL( N, SA, SX, INCX )
 !DIR$ ATTRIBUTES FORCEINLINE :: DRSCL
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DRSCL
    !DIR$ OPTIMIZE : 3
#endif
!*
!*  -- LAPACK auxiliary routine (version 3.8.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
!*
!*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SA
!*!     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION   SX( * )
!*     ..
!*
!* =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            DONE
      DOUBLE PRECISION   BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!*     ..
!*     !.. External Subroutines ..
      !EXTERNAL           DSCAL, DLABAD
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS
!*     ..
!*     .. Executable Statements ..
!*
!*     Quick return if possible
!*
     ! IF( N.LE.0 )
   
!*
!*     Get machine parameters
!*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
!*
!*     Initialize the denominator to SA and the numerator to 1.
!*
      CDEN = SA
      CNUM = ONE
!*
   10 CONTINUE
      CDEN1 = CDEN*SMLNUM
      CNUM1 = CNUM / BIGNUM
      IF( ABS( CDEN1 ).GT.ABS( CNUM ) .AND. CNUM.NE.ZERO ) THEN
!*
!*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
!*
         MUL = SMLNUM
         DONE = .FALSE.
         CDEN = CDEN1
      ELSE IF( ABS( CNUM1 ).GT.ABS( CDEN ) ) THEN
!*
!*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
!*
         MUL = BIGNUM
         DONE = .FALSE.
         CNUM = CNUM1
      ELSE
!*!
!*        Multiply X by CNUM / CDEN and return.
!*
         MUL = CNUM / CDEN
         DONE = .TRUE.
      END IF
!*
!*     Scale the vector X by MUL
!*
      CALL DSCAL( N, MUL, SX, INCX )
!*
      IF( .NOT.DONE ) &
         GO TO 10
END SUBROUTINE
     

     

!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleGEcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DGEEQUB( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, &
                      INFO ) !GCC$ ATTRIBUTES INLINE :: DGEEQUB !GCC$ ATTRIBUTES aligned(32) :: DGEEQUB
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   BIGNUM, RCMAX, RCMIN, SMLNUM, RADIX, LOGRDX
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!*     ..
!*     .. External Subroutines ..
      !EXTERNAL           XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, LOG
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DGEEQUB', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible.
!*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         ROWCND = ONE
         COLCND = ONE
         AMAX = ZERO
         RETURN
      END IF
!*
!1!*     Get machine constants.  Assume SMLNUM is a power of the radix.
!*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      RADIX = DLAMCH( 'B' )
      LOGRDX = LOG( RADIX )
!*
!*     Compute row scale factors.
      !*
      !$OMP SIMD ALIGNED(R:64) LINEAR(I:1) UNROLL PARTIAL(8)
      DO 10 I = 1, M
         R( I ) = ZERO
   10 CONTINUE
!*
!*     Find the maximum element in each row.
!*
     DO 30 J = 1, N
          !$OMP SIMD ALIGNED(R:64) LINEAR(I:1) UNROLL PARTIAL(8)  
         DO 20 I = 1, M
            R( I ) = MAX( R( I ), ABS( A( I, J ) ) )
   20    CONTINUE
   30 CONTINUE
      DO I = 1, M
         IF( R( I ).GT.ZERO ) THEN
            R( I ) = RADIX**INT( LOG( R( I ) ) / LOGRDX )
         END IF
      END DO
!*
!*     Find the maximum and minimum scale factors.
!*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 40 I = 1, M
         RCMAX = MAX( RCMAX, R( I ) )
         RCMIN = MIN( RCMIN, R( I ) )
   40 CONTINUE
      AMAX = RCMAX
!*
      IF( RCMIN.EQ.ZERO ) THEN
!*
!*        Find the first zero scale factor and return an error code.
!*
         DO 50 I = 1, M
            IF( R( I ).EQ.ZERO ) THEN
               INFO = I
               RETURN
            END IF
   50    CONTINUE
      ELSE
!*
!*        Invert the scale factors.
         !*
         !$OMP SIMD ALIGNED(R:64) LINEAR(I:1) 
         DO 60 I = 1, M
            R( I ) = ONE / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
   60    CONTINUE
!*
!*        Compute ROWCND = min(R(I)) / max(R(I)).
!*
         ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      END IF
!*
!*     Compute column scale factors
      !*
      !$OMP SIMD ALIGNED(R:64) LINEAR(I:1) UNROLL PARTIAL(8)  
      DO 70 J = 1, N
         C( J ) = ZERO
   70 CONTINUE
!*
!*     Find the maximum element in each column,
!*     assuming the row scaling computed above.
!*
      DO 90 J = 1, N
         !$OMP SIMD ALIGNED(R:64) LINEAR(I:1)  
         DO 80 I = 1, M
            C( J ) = MAX( C( J ), ABS( A( I, J ) )*R( I ) )
   80    CONTINUE
         IF( C( J ).GT.ZERO ) THEN
            C( J ) = RADIX**INT( LOG( C( J ) ) / LOGRDX )
         END IF
   90 CONTINUE
!*
!*     Find the maximum and minimum scale factors.
!*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 100 J = 1, N
         RCMIN = MIN( RCMIN, C( J ) )
         RCMAX = MAX( RCMAX, C( J ) )
  100 CONTINUE
!*
      IF( RCMIN.EQ.ZERO ) THEN
!*
!*        Find the first zero scale factor and return an error code.
!*
         DO 110 J = 1, N
            IF( C( J ).EQ.ZERO ) THEN
               INFO = M + J
               RETURN
            END IF
  110    CONTINUE
      ELSE
!*
!*        Invert the scale factors.
         !*
         !$OMP SIMD ALIGNED(R:64) LINEAR(I:1) UNROLL PARTIAL(8)  
         DO 120 J = 1, N
            C( J ) = ONE / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
  120    CONTINUE
!*
!*        Compute COLCND = min(C(J)) / max(C(J)).
!*
         COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      END IF
END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleGEcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) !GCC$ ATTRIBUTES inline :: DGETRS !GCC$ ATTRIBUTES aligned(32) :: DGETRS
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
      
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            NOTRAN
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
         LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
!      IF( N.EQ.0 .OR. NRHS.EQ.0 )
!     $   RETURN
!*
      IF( NOTRAN ) THEN
!*
!*        Solve A * X = B.
!*
!*        Apply row interchanges to the right hand sides.
!*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!*
!*        Solve L*X = B, overwriting B with X.
!*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
                    ONE, A, LDA, B, LDB )
!*
!*        Solve U*X = B, overwriting B with X.
!*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
                    NRHS, ONE, A, LDA, B, LDB )
      ELSE
!*
!*        Solve A**T * X = B.
!*
!*        Solve U**T *X = B, overwriting B with X.
!*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
                    ONE, A, LDA, B, LDB )
!*
!*        Solve L**T *X = B, overwriting B with X.
!*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, &
                    A, LDA, B, LDB )
!*
!*        Apply row interchanges to the solution vectors.
!*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF

END SUBROUTINE
    
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, &
     EQUED ) !GCC$ ATTRIBUTES hot :: DLAQGE !GCC$ ATTRIBUTES aligned(32) :: DLAQGE
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, THRESH
      PARAMETER          ( ONE = 1.0D+0, THRESH = 0.1D+0 )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   CJ, LARGE, SMALL
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!*     ..
!*     .. Executable Statements ..
!*
!*     Quick return if possible
!*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF
!*
!!*     Initialize LARGE and SMALL.
!*!
      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      LARGE = ONE / SMALL
!*
      IF( ROWCND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE ) &
          THEN
!*
!*        No row scaling
!*
         IF( COLCND.GE.THRESH ) THEN
!*
!*           No column scaling
!*
            EQUED = 'N'
         ELSE
!*
!*           Column scaling
            !*
            !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,A) PRIVATE(J,CJ,I)
            DO 20 J = 1, N
               CJ = C( J )
               !$OMP SIMD ALIGNED(A:64,C) LINEAR(I:1) UNROLL PARTIAL(6)
               DO 10 I = 1, M
                  A( I, J ) = CJ*A( I, J )
   10          CONTINUE
   20       CONTINUE
            EQUED = 'C'
         END IF
      ELSE IF( COLCND.GE.THRESH ) THEN
!*
!*        Row scaling, no column scaling
         !*
          !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,A,R,N,M) PRIVATE(J,I)
         DO 40 J = 1, N
              !$OMP SIMD ALIGNED(A:64,C,R) LINEAR(I:1) UNROLL PARTIAL(6)
            DO 30 I = 1, M
               A( I, J ) = R( I )*A( I, J )
   30       CONTINUE
   40    CONTINUE
         EQUED = 'R'
      ELSE
!*
!*        Row and column scaling
         !*
          !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,A,R,N,M) PRIVATE(J,CJ,I)
         DO 60 J = 1, N
            CJ = C( J )
            !$OMP SIMD ALIGNED(A:64,C,R) LINEAR(I:1) UNROLL PARTIAL(6)
            DO 50 I = 1, M
               A( I, J ) = CJ*R( I )*A( I, J )
   50       CONTINUE
   60    CONTINUE
         EQUED = 'B'
      END IF

END SUBROUTINE


#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLASCL2 ( M, N, D, X, LDX ) !GCC$ ATTRIBUTES inline :: DLASCL2 !GCC$ ATTRIBUTES aligned(32) :: DLASCL2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DLASCL2 ( M, N, D, X, LDX )
    !DIR$ ATTRIBUTES FORCEINLINE :: DLASCL2
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLASCL2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLASCL2
#endif
    use omp_lib
       implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            M, N, LDX
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   D( * ), X( LDX, * )
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X
!*     ..
!*
!*  =====================================================================
!!*
!*     .. Local Scalars ..
      INTEGER            I, J
!*     ..
!*     .. Executable Statements ..
!*
      DO J = 1, N
         !$OMP SIMD ALIGNED(X:64,D) LINEAR(I:1) UNROLL PARTIAL(6)
         DO I = 1, M
            X( I, J ) = X( I, J ) * D( I )
         END DO
      END DO

     
END SUBROUTINE DLASCL2


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
*!> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup doubleGEcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
   SUBROUTINE DGERFSX( TRANS, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, &
                         R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, &
                         ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, &
                         WORK, IWORK, INFO ) !GCC$ ATTRIBUTES hot :: DGERFSX !GCC$ ATTRIBUTES aligned(32) :: DGERFSX
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  ==================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   ITREF_DEFAULT, ITHRESH_DEFAULT
      DOUBLE PRECISION   COMPONENTWISE_DEFAULT, RTHRESH_DEFAULT
      DOUBLE PRECISION   DZTHRESH_DEFAULT
      PARAMETER          ( ITREF_DEFAULT = 1.0D+0 )
      PARAMETER          ( ITHRESH_DEFAULT = 10.0D+0 )
      PARAMETER          ( COMPONENTWISE_DEFAULT = 1.0D+0 )
      PARAMETER          ( RTHRESH_DEFAULT = 0.5D+0 )
      PARAMETER          ( DZTHRESH_DEFAULT = 0.25D+0 )
      INTEGER            LA_LINRX_ITREF_I, LA_LINRX_ITHRESH_I, &
                        LA_LINRX_CWISE_I
      PARAMETER          ( LA_LINRX_ITREF_I = 1, &
                        LA_LINRX_ITHRESH_I = 2 )
      PARAMETER          ( LA_LINRX_CWISE_I = 3 )
      INTEGER            LA_LINRX_TRUST_I, LA_LINRX_ERR_I, &
                        LA_LINRX_RCOND_I
      PARAMETER          ( LA_LINRX_TRUST_I = 1, LA_LINRX_ERR_I = 2 )
      PARAMETER          ( LA_LINRX_RCOND_I = 3 )
!*     ..
!*     .. Local Scalars ..
      CHARACTER(1)       NORM
      LOGICAL            ROWEQU, COLEQU, NOTRAN
      INTEGER            J, TRANS_TYPE, PREC_TYPE, REF_TYPE
      INTEGER            N_NORMS
      DOUBLE PRECISION   ANORM, RCOND_TMP
      DOUBLE PRECISION   ILLRCOND_THRESH, ERR_LBND, CWISE_WRONG
      LOGICAL            IGNORE_CWISE
      INTEGER            ITHRESH
      DOUBLE PRECISION   RTHRESH, UNSTABLE_THRESH
!*     ..
!*     .. External Subroutines ..
     ! EXTERNAL           XERBLA, DGECON, DLA_GERFSX_EXTENDED
!*     ..
*!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!*     ..
!*     .. External Functions ..
      EXTERNAL           LSAME, ILATRANS, ILAPREC
      EXTERNAL           DLAMCH, DLANGE, DLA_GERCOND
      DOUBLE PRECISION   DLAMCH, DLANGE, DLA_GERCOND
      LOGICAL            LSAME
      INTEGER            ILATRANS, ILAPREC
!*     ..
!*     .. Executable Statements ..
!*
!!*     Check the input parameters.
!*
      INFO = 0
      TRANS_TYPE = ILATRANS( TRANS )
      REF_TYPE = INT( ITREF_DEFAULT )
      IF ( NPARAMS .GE. LA_LINRX_ITREF_I ) THEN
         IF ( PARAMS( LA_LINRX_ITREF_I ) .LT. 0.0D+0 ) THEN
            PARAMS( LA_LINRX_ITREF_I ) = ITREF_DEFAULT
         ELSE
            REF_TYPE = PARAMS( LA_LINRX_ITREF_I )
         END IF
      END IF
!*
!*     Set default parameters.
!*
      ILLRCOND_THRESH = DBLE( N ) * DLAMCH( 'Epsilon' )
      ITHRESH = INT( ITHRESH_DEFAULT )
      RTHRESH = RTHRESH_DEFAULT
      UNSTABLE_THRESH = DZTHRESH_DEFAULT
      IGNORE_CWISE = COMPONENTWISE_DEFAULT .EQ. 0.0D+0
!*
      IF ( NPARAMS.GE.LA_LINRX_ITHRESH_I ) THEN
         IF ( PARAMS( LA_LINRX_ITHRESH_I ).LT.0.0D+0 ) THEN
            PARAMS( LA_LINRX_ITHRESH_I ) = ITHRESH
         ELSE
            ITHRESH = INT( PARAMS( LA_LINRX_ITHRESH_I ) )
         END IF
      END IF
      IF ( NPARAMS.GE.LA_LINRX_CWISE_I ) THEN
         IF ( PARAMS( LA_LINRX_CWISE_I ).LT.0.0D+0 ) THEN
            IF ( IGNORE_CWISE ) THEN
               PARAMS( LA_LINRX_CWISE_I ) = 0.0D+0
            ELSE
               PARAMS( LA_LINRX_CWISE_I ) = 1.0D+0
            END IF
         ELSE
            IGNORE_CWISE = PARAMS( LA_LINRX_CWISE_I ) .EQ. 0.0D+0
         END IF
      END IF
      IF ( REF_TYPE .EQ. 0 .OR. N_ERR_BNDS .EQ. 0 ) THEN
         N_NORMS = 0
      ELSE IF ( IGNORE_CWISE ) THEN
         N_NORMS = 1
      ELSE
         N_NORMS = 2
      END IF
!*
      NOTRAN = LSAME( TRANS, 'N' )
      ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
      COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
!*
!*     Test input parameters.
!*
      IF( TRANS_TYPE.EQ.-1 ) THEN
        INFO = -1
      ELSE IF( .NOT.ROWEQU .AND. .NOT.COLEQU .AND. &
              .NOT.LSAME( EQUED, 'N' ) ) THEN
        INFO = -2
      ELSE IF( N.LT.0 ) THEN
        INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
        INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
        INFO = -6
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
        INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
        INFO = -13
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
        INFO = -15
      END IF
      IF( INFO.NE.0 ) THEN
        !CALL XERBLA( 'DGERFSX', -INFO )
        RETURN
      END IF
!*
!*     Quick return if possible.
!*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         RCOND = 1.0D+0
         DO J = 1, NRHS
            BERR( J ) = 0.0D+0
            IF ( N_ERR_BNDS .GE. 1 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
            END IF
            IF ( N_ERR_BNDS .GE. 2 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I) = 0.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 0.0D+0
            END IF
            IF ( N_ERR_BNDS .GE. 3 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_RCOND_I) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = 1.0D+0
            END IF
         END DO
         RETURN
      END IF
!*
!*     Default to failure.
!*
      RCOND = 0.0D+0
      DO J = 1, NRHS
         BERR( J ) = 1.0D+0
         IF ( N_ERR_BNDS .GE. 1 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
         END IF
         IF ( N_ERR_BNDS .GE. 2 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0
         END IF
         IF ( N_ERR_BNDS .GE. 3 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = 0.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = 0.0D+0
         END IF
      END DO
!*
!*     Compute the norm of A and the reciprocal of the condition
!*     number of A.
!*
      IF( NOTRAN ) THEN
         NORM = 'I'
      ELSE
         NORM = '1'
      END IF
      ANORM = DLANGE( NORM, N, N, A, LDA, WORK )
      CALL DGECON( NORM, N, AF, LDAF, ANORM, RCOND, WORK, IWORK, INFO )
!*
!*     Perform refinement on each right-hand side
!*
      IF ( REF_TYPE .NE. 0 ) THEN

         PREC_TYPE = ILAPREC( 'E' )

         IF ( NOTRAN ) THEN
            CALL DLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE,  N, &
                NRHS, A, LDA, AF, LDAF, IPIV, COLEQU, C, B, &
                LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM, &
                ERR_BNDS_COMP, WORK(N+1), WORK(1), WORK(2*N+1), &
                WORK(1), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH, &
                IGNORE_CWISE, INFO )
         ELSE
            CALL DLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE,  N, &
                NRHS, A, LDA, AF, LDAF, IPIV, ROWEQU, R, B, &
                LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM, &
                ERR_BNDS_COMP, WORK(N+1), WORK(1), WORK(2*N+1), &
                WORK(1), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH, &
                IGNORE_CWISE, INFO )
         END IF
      END IF

      ERR_LBND = MAX( 10.0D+0, SQRT( DBLE( N ) ) ) * DLAMCH( 'Epsilon' )
      IF ( N_ERR_BNDS .GE. 1 .AND. N_NORMS .GE. 1 ) THEN
!*
!*     Compute scaled normwise condition number cond(A*C).
!*
         IF ( COLEQU .AND. NOTRAN ) THEN
            RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV, &
                -1, C, INFO, WORK, IWORK )
         ELSE IF ( ROWEQU .AND. .NOT. NOTRAN ) THEN
            RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV, &
                -1, R, INFO, WORK, IWORK )
         ELSE
            RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV, &
                0, R, INFO, WORK, IWORK )
         END IF
         DO J = 1, NRHS
!*
!*     Cap the error at 1.0.
!*
            IF ( N_ERR_BNDS .GE. LA_LINRX_ERR_I &
                .AND. ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) .GT. 1.0D+0 ) &
                ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0
!*
!*     Threshold the error (see LAWN).
!*
            IF ( RCOND_TMP .LT. ILLRCOND_THRESH ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 0.0D+0
               IF ( INFO .LE. N ) INFO = N + J
            ELSE IF ( ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) .LT. ERR_LBND ) &
            THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = ERR_LBND
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0D+0
            END IF
!*
!*     Save the condition number.
!*
            IF ( N_ERR_BNDS .GE. LA_LINRX_RCOND_I ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = RCOND_TMP
            END IF
         END DO
      END IF

      IF ( N_ERR_BNDS .GE. 1 .AND. N_NORMS .GE. 2 ) THEN
!*
!*     Compute componentwise condition number cond(A*diag(Y(:,J))) for
!*     each right-hand side using the current solution as an estimate of
!*     the true solution.  If the componentwise error estimate is too
!*     large, then the solution is a lousy estimate of truth and the
!*     estimated RCOND may be too optimistic.  To avoid misleading users,
!*     the inverse condition number is set to 0.0 when the estimated
!*     cwise error is at least CWISE_WRONG.
!*
         CWISE_WRONG = SQRT( DLAMCH( 'Epsilon' ) )
         DO J = 1, NRHS
            IF ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .LT. CWISE_WRONG ) &
                THEN 
               RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, &
                   IPIV, 1, X(1,J), INFO, WORK, IWORK )
            ELSE
               RCOND_TMP = 0.0D+0
            END IF
!*
!*     Cap the error at 1.0.
!*
            IF ( N_ERR_BNDS .GE. LA_LINRX_ERR_I &
                .AND. ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .GT. 1.0D+0 ) &
                ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0
!*
!*     Threshold the error (see LAWN).
!!*
            IF ( RCOND_TMP .LT. ILLRCOND_THRESH ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 0.0D+0
               IF ( PARAMS( LA_LINRX_CWISE_I ) .EQ. 1.0D+0 &
                   .AND. INFO.LT.N + J ) INFO = N + J
            ELSE IF ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) &
                   .LT. ERR_LBND ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = ERR_LBND
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
            END IF
!*
!*     Save the condition number.
!*
            IF ( N_ERR_BNDS .GE. LA_LINRX_RCOND_I ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = RCOND_TMP
            END IF
         END DO
      END IF

END SUBROUTINE



!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date June 2017
!*
!*> \ingroup doubleGEcomputational
!*
!*  =====================================================================

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE, N, NRHS, A, &
                                     LDA, AF, LDAF, IPIV, COLEQU, C, B, &
                                     LDB, Y, LDY, BERR_OUT, N_NORMS, &
                                     ERRS_N, ERRS_C, RES, AYB, DY, &
                                     Y_TAIL, RCOND, ITHRESH, RTHRESH, &
                                     DZ_UB, IGNORE_CWISE, INFO ) !GCC$ ATTRIBUTES hot :: DLA_GERFSX_EXTENDED !GCC$ ATTRIBUTES aligned(32) :: DLA_GERFSX_EXTENDED
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
      
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      CHARACTER          TRANS
      INTEGER            CNT, I, J, X_STATE, Z_STATE, Y_PREC_STATE
      DOUBLE PRECISION   YK, DYK, YMIN, NORMY, NORMX, NORMDX, DXRAT, &
                        DZRAT, PREVNORMDX, PREV_DZ_Z, DXRATMAX, &
                        DZRATMAX, DX_X, DZ_Z, FINAL_DX_X, FINAL_DZ_Z, &
                        EPS, HUGEVAL, INCR_THRESH
      LOGICAL            INCR_PREC
!*     ..
!*     .. Parameters ..
      INTEGER            UNSTABLE_STATE, WORKING_STATE, CONV_STATE, &
                        NOPROG_STATE, BASE_RESIDUAL, EXTRA_RESIDUAL, &
                        EXTRA_Y
      PARAMETER          ( UNSTABLE_STATE = 0, WORKING_STATE = 1, &
                        CONV_STATE = 2, NOPROG_STATE = 3 )
      PARAMETER          ( BASE_RESIDUAL = 0, EXTRA_RESIDUAL = 1, &
                        EXTRA_Y = 2 )
      INTEGER            FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I
      INTEGER            RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I
      INTEGER            CMP_ERR_I, PIV_GROWTH_I
      PARAMETER          ( FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, &
                          BERR_I = 3 )
      PARAMETER          ( RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 )
      PARAMETER          ( CMP_RCOND_I = 7, CMP_ERR_I = 8, &
                        PIV_GROWTH_I = 9 )
      INTEGER            LA_LINRX_ITREF_I, LA_LINRX_ITHRESH_I, &
                        LA_LINRX_CWISE_I
      PARAMETER          ( LA_LINRX_ITREF_I = 1, &
                        LA_LINRX_ITHRESH_I = 2 )
      PARAMETER          ( LA_LINRX_CWISE_I = 3 )
      INTEGER            LA_LINRX_TRUST_I, LA_LINRX_ERR_I, &
                        LA_LINRX_RCOND_I
      PARAMETER          ( LA_LINRX_TRUST_I = 1, LA_LINRX_ERR_I = 2 )
      PARAMETER          ( LA_LINRX_RCOND_I = 3 )
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGETRS, DGEMV, BLAS_DGEMV_X, &
                        BLAS_DGEMV2_X,  DLAMCH
                       
      DOUBLE PRECISION   DLAMCH
      
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
      IF ( INFO.NE.0 ) RETURN
      TRANS = CHLA_TRANSTYPE(TRANS_TYPE)
      EPS = DLAMCH( 'Epsilon' )
      HUGEVAL = DLAMCH( 'Overflow' )
!*     Force HUGEVAL to Inf
      HUGEVAL = HUGEVAL * HUGEVAL
!*     Using HUGEVAL may lead to spurious underflows.
      INCR_THRESH = DBLE( N ) * EPS
!*
      DO J = 1, NRHS
         Y_PREC_STATE = EXTRA_RESIDUAL
         IF ( Y_PREC_STATE .EQ. EXTRA_Y ) THEN
            DO I = 1, N
               Y_TAIL( I ) = 0.0D+0
            END DO
         END IF

         DXRAT = 0.0D+0
         DXRATMAX = 0.0D+0
         DZRAT = 0.0D+0
         DZRATMAX = 0.0D+0
         FINAL_DX_X = HUGEVAL
         FINAL_DZ_Z = HUGEVAL
         PREVNORMDX = HUGEVAL
         PREV_DZ_Z = HUGEVAL
         DZ_Z = HUGEVAL
         DX_X = HUGEVAL

         X_STATE = WORKING_STATE
         Z_STATE = UNSTABLE_STATE
         INCR_PREC = .FALSE.

         DO CNT = 1, ITHRESH
!*
!*         Compute residual RES = B_s - op(A_s) * Y,
!*             op(A) = A, A**T, or A**H depending on TRANS (and type).
!*
            CALL DCOPY( N, B( 1, J ), 1, RES, 1 )
            IF ( Y_PREC_STATE .EQ. BASE_RESIDUAL ) THEN
               CALL DGEMV( TRANS, N, N, -1.0D+0, A, LDA, Y( 1, J ), 1, &
                   1.0D+0, RES, 1 )
            ELSE IF ( Y_PREC_STATE .EQ. EXTRA_RESIDUAL ) THEN
               CALL BLAS_DGEMV_X( TRANS_TYPE, N, N, -1.0D+0, A, LDA, &
                   Y( 1, J ), 1, 1.0D+0, RES, 1, PREC_TYPE )
            ELSE
               CALL BLAS_DGEMV2_X( TRANS_TYPE, N, N, -1.0D+0, A, LDA, &
                  Y( 1, J ), Y_TAIL, 1, 1.0D+0, RES, 1, PREC_TYPE )
            END IF

!        XXX: RES is no longer needed.
            CALL DCOPY( N, RES, 1, DY, 1 )
            CALL DGETRS( TRANS, N, 1, AF, LDAF, IPIV, DY, N, INFO )
!*
!*         Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT.
!*
            NORMX = 0.0D+0
            NORMY = 0.0D+0
            NORMDX = 0.0D+0
            DZ_Z = 0.0D+0
            YMIN = HUGEVAL
!*
            DO I = 1, N
               YK = ABS( Y( I, J ) )
               DYK = ABS( DY( I ) )

               IF ( YK .NE. 0.0D+0 ) THEN
                  DZ_Z = MAX( DZ_Z, DYK / YK )
               ELSE IF ( DYK .NE. 0.0D+0 ) THEN
                  DZ_Z = HUGEVAL
               END IF

               YMIN = MIN( YMIN, YK )

               NORMY = MAX( NORMY, YK )

               IF ( COLEQU ) THEN
                  NORMX = MAX( NORMX, YK * C( I ) )
                  NORMDX = MAX( NORMDX, DYK * C( I ) )
               ELSE
                  NORMX = NORMY
                  NORMDX = MAX( NORMDX, DYK )
               END IF
            END DO

            IF ( NORMX .NE. 0.0D+0 ) THEN
               DX_X = NORMDX / NORMX
            ELSE IF ( NORMDX .EQ. 0.0D+0 ) THEN
               DX_X = 0.0D+0
            ELSE
               DX_X = HUGEVAL
            END IF

            DXRAT = NORMDX / PREVNORMDX
            DZRAT = DZ_Z / PREV_DZ_Z
!*
!*         Check termination criteria
!*
            IF (.NOT.IGNORE_CWISE &
                .AND. YMIN*RCOND .LT. INCR_THRESH*NORMY &
                .AND. Y_PREC_STATE .LT. EXTRA_Y) &
                INCR_PREC = .TRUE.

            IF ( X_STATE .EQ. NOPROG_STATE .AND. DXRAT .LE. RTHRESH ) &
                 X_STATE = WORKING_STATE
            IF ( X_STATE .EQ. WORKING_STATE ) THEN
               IF ( DX_X .LE. EPS ) THEN
                  X_STATE = CONV_STATE
               ELSE IF ( DXRAT .GT. RTHRESH ) THEN
                  IF ( Y_PREC_STATE .NE. EXTRA_Y ) THEN
                     INCR_PREC = .TRUE.
                  ELSE
                     X_STATE = NOPROG_STATE
                  END IF
               ELSE
                  IF ( DXRAT .GT. DXRATMAX ) DXRATMAX = DXRAT
               END IF
               IF ( X_STATE .GT. WORKING_STATE ) FINAL_DX_X = DX_X
            END IF

            IF ( Z_STATE .EQ. UNSTABLE_STATE .AND. DZ_Z .LE. DZ_UB ) &
                Z_STATE = WORKING_STATE
            IF ( Z_STATE .EQ. NOPROG_STATE .AND. DZRAT .LE. RTHRESH ) &
                Z_STATE = WORKING_STATE
            IF ( Z_STATE .EQ. WORKING_STATE ) THEN
               IF ( DZ_Z .LE. EPS ) THEN
                  Z_STATE = CONV_STATE
               ELSE IF ( DZ_Z .GT. DZ_UB ) THEN
                  Z_STATE = UNSTABLE_STATE
                  DZRATMAX = 0.0D+0
                  FINAL_DZ_Z = HUGEVAL
               ELSE IF ( DZRAT .GT. RTHRESH ) THEN
                  IF ( Y_PREC_STATE .NE. EXTRA_Y ) THEN
                     INCR_PREC = .TRUE.
                  ELSE
                     Z_STATE = NOPROG_STATE
                  END IF
               ELSE
                  IF ( DZRAT .GT. DZRATMAX ) DZRATMAX = DZRAT
               END IF
               IF ( Z_STATE .GT. WORKING_STATE ) FINAL_DZ_Z = DZ_Z
            END IF
!*
!*           Exit if both normwise and componentwise stopped working,
!*           but if componentwise is unstable, let it go at least two
!*           iterations.
!*
            IF ( X_STATE.NE.WORKING_STATE ) THEN
               IF ( IGNORE_CWISE) GOTO 666
               IF ( Z_STATE.EQ.NOPROG_STATE .OR. Z_STATE.EQ.CONV_STATE ) &
                   GOTO 666
               IF ( Z_STATE.EQ.UNSTABLE_STATE .AND. CNT.GT.1 ) GOTO 666
            END IF

            IF ( INCR_PREC ) THEN
               INCR_PREC = .FALSE.
               Y_PREC_STATE = Y_PREC_STATE + 1
               DO I = 1, N
                  Y_TAIL( I ) = 0.0D+0
               END DO
            END IF

            PREVNORMDX = NORMDX
            PREV_DZ_Z = DZ_Z
!*
!*           Update soluton.
!*
            IF ( Y_PREC_STATE .LT. EXTRA_Y ) THEN
               CALL DAXPY( N, 1.0D+0, DY, 1, Y( 1, J ), 1 )
            ELSE
               CALL DLA_WWADDW( N, Y( 1, J ), Y_TAIL, DY )
            END IF

         END DO
!*        Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT.
 666     CONTINUE
!*
!*     Set final_* when cnt hits ithresh.
!*
         IF ( X_STATE .EQ. WORKING_STATE ) FINAL_DX_X = DX_X
         IF ( Z_STATE .EQ. WORKING_STATE ) FINAL_DZ_Z = DZ_Z
!1*
!*     Compute error bounds
!*
         IF (N_NORMS .GE. 1) THEN
            ERRS_N( J, LA_LINRX_ERR_I ) = FINAL_DX_X / (1 - DXRATMAX)
         END IF
         IF ( N_NORMS .GE. 2 ) THEN
            ERRS_C( J, LA_LINRX_ERR_I ) = FINAL_DZ_Z / (1 - DZRATMAX)
         END IF
!*
!*     Compute componentwise relative backward error from formula
!*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
!*     where abs(Z) is the componentwise absolute value of the matrix
!*     or vector Z.
!*
!*         Compute residual RES = B_s - op(A_s) * Y,
!*             op(A) = A, A**T, or A**H depending on TRANS (and type).
!*
         CALL DCOPY( N, B( 1, J ), 1, RES, 1 )
         CALL DGEMV( TRANS, N, N, -1.0D+0, A, LDA, Y(1,J), 1, 1.0D+0, &
          RES, 1 )

         DO I = 1, N
            AYB( I ) = ABS( B( I, J ) )
         END DO
!*
!*     Compute abs(op(A_s))*abs(Y) + abs(B_s).
!*
         CALL DLA_GEAMV ( TRANS_TYPE, N, N, 1.0D+0, &
             A, LDA, Y(1, J), 1, 1.0D+0, AYB, 1 )

         CALL DLA_LIN_BERR ( N, N, 1, RES, AYB, BERR_OUT( J ) )
!*
!*     End of loop for each RHS.
!*
      END DO

END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date June 2017
!*
!*> \ingroup doubleGEcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      SUBROUTINE DLA_GEAMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, &
           Y, INCY ) !GCC$ ATTRIBUTES hot :: DLA_GEAMV !GCC$ ATTRIBUTES aligned(32) :: DLA_GEAMV
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            SYMB_ZERO
      DOUBLE PRECISION   TEMP, SAFE1
      INTEGER            I, INFO, IY, J, JX, KX, KY, LENX, LENY
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DLAMCH
      DOUBLE PRECISION   DLAMCH
!*     ..
!*     .. External Functions ..
      EXTERNAL           ILATRANS
      INTEGER            ILATRANS
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, ABS, SIGN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      IF     ( .NOT.( ( TRANS.EQ.ILATRANS( 'N' ) ) &
                .OR. ( TRANS.EQ.ILATRANS( 'T' ) ) &
                .OR. ( TRANS.EQ.ILATRANS( 'C' )) ) ) THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         !CALL XERBLA( 'DLA_GEAMV ', INFO )
         RETURN
      END IF
!*
!*     Quick return if possible.
!*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
         ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) &
        RETURN
!*
!*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!*     up the start points in  X  and  Y.
!*
      IF( TRANS.EQ.ILATRANS( 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!*
!*     Set SAFE1 essentially to be the underflow threshold times the
!*     number of additions in each row.
!*
      SAFE1 = DLAMCH( 'Safe minimum' )
      SAFE1 = (N+1)*SAFE1
!*
!*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y).
!*
!*     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to
!*     the inexact flag.  Still doesn't help change the iteration order
!*     to per-column.
!*
      IY = KY
      IF ( INCX.EQ.1 ) THEN
         IF( TRANS.EQ.ILATRANS( 'N' ) )THEN
            !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(Y,X,BETA,LENY,LENX,INCX) FIRSTPRIVATE(IY) PRIVATE(I,SYMB_ZERO,J,TEMP)
            DO I = 1, LENY
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. ZERO ) THEN
                  !$OMP SIMD ALIGNED(X:64,Y) LINEAR(J:1)
                  DO J = 1, LENX
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. &
                         ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO ) &
                   Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
            !$OMP END PARALLEL DO
         ELSE
            !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(Y,X,LENY,LENX,INCY,ALPHA) FIRSTPRIVATE(IY) PRIVATE(I,SYMB_ZERO,J,TEMP)
            DO I = 1, LENY
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. ZERO ) THEN
                    !$OMP SIMD ALIGNED(X:64,Y) LINEAR(J:1)
                  DO J = 1, LENX
                     TEMP = ABS( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. &
                          ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( J ) )*TEMP
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO ) &
                   Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
            !$OMP END PARALLEL DO
         END IF
      ELSE
         IF( TRANS.EQ.ILATRANS( 'N' ) )THEN
             !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(Y,X,LENY,LENX,ALPHA) FIRSTPRIVATE(IY) PRIVATE(I,SYMB_ZERO,JX,J,TEMP)
            DO I = 1, LENY
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. ZERO ) THEN
                  JX = KX
                     !$OMP SIMD ALIGNED(X:64,Y) LINEAR(J:1)
                  DO J = 1, LENX
                     TEMP = ABS( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND. &
                        ( X( JX ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               END IF

               IF (.NOT.SYMB_ZERO) &
                   Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
            !$OMP END PARALLEL DO
         ELSE
             !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(Y,X,LENY,LENX,ALPHA) FIRSTPRIVATE(IY) PRIVATE(I,SYMB_ZERO,JX,J,TEMP)
            DO I = 1, LENY
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0D+0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. ZERO ) THEN
                  JX = KX
                   !$OMP SIMD ALIGNED(X:64,Y) LINEAR(J:1)
                  DO J = 1, LENX
                     TEMP = ABS( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND. &
                         ( X( JX ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*ABS( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               END IF

               IF (.NOT.SYMB_ZERO) &
                   Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
            !$OMP END PARALLEL DO
         END IF

      END IF

END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLA_WWADDW( N, X, Y, W ) !GCC$ ATTRIBUTES inline :: DLA_WWADDW !GCC$ ATTRIBUTES aligned(32) :: DLA_WWADDW
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DLA_WWADDW( N, X, Y, W )
    !DIR$ ATTRIBUTES FORCEINLINE :: DLA_WWADDW
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLA_WWADDW
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DLA_WWADDW
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
      INTEGER            N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION   X( * ), Y( * ), W( * )
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X, Y, W
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      DOUBLE PRECISION   S
      INTEGER            I
!*     ..
!*     .. Executable Statements ..
      !*
      !$OMP SIMD ALIGNED(X:64,Y,W) LINEAR(I:1) UNROLL PARTIAL(6)
      DO 10 I = 1, N
        S = X(I) + W(I)
        S = (S + S) - S
        Y(I) = ((X(I) - S) + W(I)) + Y(I)
        X(I) = S
 10   CONTINUE
     
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR ) !GCC$ ATTRIBUTES inline :: DLA_LIN_BERR !GCC$ ATTRIBUTES aligned(32) :: DLA_LIN_BERR
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR )
     !DIR$ ATTRIBUTES FORCEINLINE :: DLA_LIN_BERR
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DLA_LIN_BERR
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DLA_LIN_BERR
#endif
      implicit none
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            N, NZ, NRHS
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION   AYB( N, NRHS ), BERR( NRHS )
      DOUBLE PRECISION   RES( N, NRHS )
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      DOUBLE PRECISION   TMP
      INTEGER            I, J
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!*     ..
!*     .. External Functions ..
      EXTERNAL           DLAMCH
      DOUBLE PRECISION   DLAMCH
      DOUBLE PRECISION   SAFE1
!*     ..
!*     .. Executable Statements ..
!*
!!*     Adding SAFE1 to the numerator guards against spuriously zero
!*     residuals.  A similar safeguard is in the SLA_yyAMV routine used
!*     to compute AYB.
!*
      SAFE1 = DLAMCH( 'Safe minimum' )
      SAFE1 = (NZ+1)*SAFE1

      DO J = 1, NRHS
         BERR(J) = 0.0D+0
         DO I = 1, N
            IF (AYB(I,J) .NE. 0.0D+0) THEN
               TMP = (SAFE1+ABS(RES(I,J)))/AYB(I,J)
               BERR(J) = MAX( BERR(J), TMP )
            END IF
!*
!*     If AYB is exactly 0.0 (and if computed by SLA_yyAMV), then we know
!*     the true residual also must be exactly 0.0.
!*
         END DO
      END DO
END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
CHARACTER*1 FUNCTION CHLA_TRANSTYPE( TRANS ) !GCC$ ATTRIBUTES inline :: CHLA_TRANSTYPE !GCC$ ATTRIBUTES aligned(32) :: CHLA_TRANSTYPE
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  CHARACTER*1 FUNCTION CHLA_TRANSTYPE( TRANS )
     !DIR$ ATTRIBUTES FORCEINLINE :: CHLA_TRANSTYPE
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: CHLA_TRANSTYPE
#endif
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            TRANS
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      INTEGER BLAS_NO_TRANS, BLAS_TRANS, BLAS_CONJ_TRANS
      PARAMETER ( BLAS_NO_TRANS = 111, BLAS_TRANS = 112, &
          BLAS_CONJ_TRANS = 113 )
!*     ..
!*     .. Executable Statements ..
      IF( TRANS.EQ.BLAS_NO_TRANS ) THEN
         CHLA_TRANSTYPE = 'N'
      ELSE IF( TRANS.EQ.BLAS_TRANS ) THEN
         CHLA_TRANSTYPE = 'T'
      ELSE IF( TRANS.EQ.BLAS_CONJ_TRANS ) THEN
         CHLA_TRANSTYPE = 'C'
      ELSE
         CHLA_TRANSTYPE = 'X'
      END IF
   
END FUNCTION 

    
!*  Authors:
!*  ========
!*
!C> \author Univ. of Tennessee
!C> \author Univ. of California Berkeley
!C> \author Univ. of Colorado Denver
!C> \author NAG Ltd.
!*
!C> \date December 2016
!*
!C> \ingroup variantsPOcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DPOTRF ( UPLO, N, A, LDA, INFO ) !GCC$ ATTRIBUTES hot :: DPOTRF !GCC$ ATTRIBUTES aligned(32) :: DPOTRF
#elif defined(__INTEL_COMPILER) || defined(__ICC)
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
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DGEMM, DSYRK, DTRSM
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DPOTRF', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
    !  IF( N.EQ.0 )
     !$   RETURN
!*
!*     Determine the block size for this environment.
!*
      NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
!*
!*        Use unblocked code.
!*
         CALL DPOTF2( UPLO, N, A, LDA, INFO )
      ELSE
!*
!*        Use blocked code.
!*
         IF( UPPER ) THEN
!*
!*           Compute the Cholesky factorization A = U'*U.
!*
            DO 10 J = 1, N, NB
!*
!*              Update and factorize the current diagonal block and test
!*              for non-positive-definiteness.
!*
               JB = MIN( NB, N-J+1 )

               CALL DPOTF2( 'Upper', JB, A( J, J ), LDA, INFO )

               IF( INFO.NE.0 ) &
                 GO TO 30

               IF( J+JB.LE.N ) THEN
!*
!*                 Updating the trailing submatrix.
!*
                  CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', &
                             JB, N-J-JB+1, ONE, A( J, J ), LDA, &
                             A( J, J+JB ), LDA )
                  CALL DSYRK( 'Upper', 'Transpose', N-J-JB+1, JB, -ONE, &
                             A( J, J+JB ), LDA, &
                             ONE, A( J+JB, J+JB ), LDA )
               END IF
   10       CONTINUE

         ELSE
!*
!*           Compute the Cholesky factorization A = L*L'.
!*
            DO 20 J = 1, N, NB
!*
!*              Update and factorize the current diagonal block and test
!*              for non-positive-definiteness.
!*
               JB = MIN( NB, N-J+1 )

               CALL DPOTF2( 'Lower', JB, A( J, J ), LDA, INFO )

               IF( INFO.NE.0 ) &
                 GO TO 30

               IF( J+JB.LE.N ) THEN
!!*
!*                Updating the trailing submatrix.
!*
                 CALL DTRSM( 'Right', 'Lower', 'Transpose', 'Non-unit', &
                            N-J-JB+1, JB, ONE, A( J, J ), LDA, &
                            A( J+JB, J ), LDA )

                 CALL DSYRK( 'Lower', 'No Transpose', N-J-JB+1, JB, &
                            -ONE, A( J+JB, J ), LDA, &
                            ONE, A( J+JB, J+JB ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40

   30 CONTINUE
      INFO = INFO + J - 1

   40 CONTINUE
     
END SUBROUTINE


!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO ) !GCC$ ATTRIBUTES inline :: DPOTF2 !GCC$ ATTRIBUTES aligned(32) :: DPOTF2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )
    !DIR$ ATTRIBUTES FORCEINLINE :: DPOTF2
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DPOTF2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DPOTF2
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
      INTEGER            INFO, LDA, N
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
!!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT, DISNAN
!*     ..
!*     .. External Subroutines ..
      EXTERNAL           DGEMV, DSCAL
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         !CALL XERBLA( 'DPOTF2', -INFO )
         RETURN
      END IF
!*
!*     Quick return if possible
!*
     ! IF( N.EQ.0 )
     !$   RETURN
!*
      IF( UPPER ) THEN
!*
!*        Compute the Cholesky factorization A = U**T *U.
!*
         DO 10 J = 1, N
!*
!*           Compute U(J,J) and test for non-positive-definiteness.
!*
            AJJ = A( J, J ) - DDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            IF( AJJ.LE.ZERO.OR.DISNAN( AJJ ) ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!*
!*           Compute elements J+1:N of row J.
!*
            IF( J.LT.N ) THEN
               CALL DGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ), &
                          LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
!*
!*        Compute the Cholesky factorization A = L*L**T.
!*
         DO 20 J = 1, N
!*
!*           Compute L(J,J) and test for non-positive-definiteness.
!*
            AJJ = A( J, J ) - DDOT( J-1, A( J, 1 ), LDA, A( J, 1 ), &
                 LDA )
            IF( AJJ.LE.ZERO.OR.DISNAN( AJJ ) ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!*
!*           Compute elements J+1:N of column J.
!*
            IF( J.LT.N ) THEN
               CALL DGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ), &
                          LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40

   30 CONTINUE
      INFO = J

   40 CONTINUE
   
END SUBROUTINE



