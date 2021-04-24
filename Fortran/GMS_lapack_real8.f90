


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
    use omp_lib
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
         !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(WORK,A) PRIVATE(J,I) COLLAPSE(2) IF(N>=400)   
         DO 70 J = 1, N
                !$OMP SIMD ALIGNED(A:64) LINEAR(I:1)  UNROLL PARTIAL(6)
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
70       CONTINUE
         !$OMP END PARALLEL DO      
         VALUE = ZERO
         DO 80 I = 1, M
            TEMP = WORK( I )
            IF( VALUE.LT.TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
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
         !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(B,A) PRIVATE(J,I) IF(N>=200)
         DO 20 J = 1, N
            !$OMP SIMD ALIGNED(A:64,B) LINEAR(I:1) UNROLL PARTIAL(8)
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
20       CONTINUE
          !$OMP END PARALLEL DO     
     ELSE IF( LSAME( UPLO, 'L' ) ) THEN
            !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(B,A) COLLAPSE(2) PRIVATE(J,I) IF(N>=200)    
        DO 40 J = 1, N
              !$OMP SIMD ALIGNED(A:64,B) LINEAR(I:1) UNROLL PARTIAL(8)
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
40       CONTINUE
         !$OMP END PARALLEL DO      
     ELSE
           !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(B,A) COLLAPSE(2) PRIVATE(J,I) IF(N>=200)         
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
         !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(A) PRIVATE(J,I) IF(N>=400)
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
          !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(A) PRIVATE(J,I) IF(N>=400)
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
          !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(A) COLLAPSE(2) PRIVATE(J,I) IF(N>=100)
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
       !$OMP SIMD ALIGNED(A:64) LINEAR(I:1) 
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
               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J) IF(K>=200)
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
               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J,I) COLLAPSE(2)  IF(K>=100)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J) IF(K>=200)
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
                 !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J,I) COLLAPSE(2)  IF(K>=100)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J)   IF(K>=200)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J,I)   IF(K>=100)
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
!*           !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J)   IF(K>=100)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J,I)   IF(K>=100)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J,I)   IF(K>=100)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J)   IF(K>=100)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J,I)   IF(K>=100)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J)   IF(K>=100)
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
                 !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J,I)   IF(K>=100)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J)   IF(K>=100)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK) PRIVATE(J,I)   IF(K>=100)
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

    

    
