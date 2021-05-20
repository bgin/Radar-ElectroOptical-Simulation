



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
      EXTERNAL           DBDSQR,DGEMM,DLACPY, &
                         DLASCL, DLASET
                        
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

       
