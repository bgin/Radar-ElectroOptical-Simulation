


  
#if 0

*
*  Arguments:
*  ==========
*
*> \param[in] NORM
*> \verbatim
*>          NORM is CHARACTER*1
*>          Specifies the value to be returned in ZLANGE as described
*>          above.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.  When M = 0,
*>          ZLANGE is set to zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.  When N = 0,
*>          ZLANGE is set to zero.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The m by n matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(M,1).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
*>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
*>          referenced.
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
*> \ingroup complex16GEauxiliary
*
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
DOUBLE PRECISION FUNCTION ZLANGE(NORM, M, N, A, LDA, WORK) !GCC$ ATTRIBUTES HOT :: ZLANGE !GCC$ ATTRIBUTES aligned(32) :: ZLANGE !GCC$ ATTRIBUTES no_stack_protector :: ZLANGE
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  DOUBLE PRECISION FUNCTION ZLANGE(NORM,M,N,A,LDA,WORK)
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZLANGE
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZLANGE
#endif
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
    !*

    use omp_lib

      IMPLICIT NONE
!*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
!*     ..
      !*     .. Array Arguments ..
      
      !DOUBLE PRECISION   WORK( * )
      !COMPLEX*16         A( LDA, * )
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: WORK
      COMPLEX(16),      DIMENSION(:,:), ALLOCATABLE :: A
!#elif defined(__ICC) || defined(__INTEL_COMPILER)
!      DOUBLE PRECISION  WORK(*)
!      COMPLEX*16         A( LDA, * )
!      !DIR$ ASSUME_ALIGNED WORK:64
!      !DIR$ ASSUME_ALIGNED A:64
!#endif
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
!*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           ZLASSQ, DCOMBSSQ
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

         !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(A,N,M) REDUCTION(+:SUM) FIRSTPRIVATE(VALUE) PRIVATE(J,I) COLLAPSE(2)
         DO 40 J = 1, N
            !$OMP SINGLE
            SUM = ZERO
            !$OMP END SINGLE
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
        30  CONTINUE
           !$OMP CRITICAL    
               IF( VALUE.LT.SUM .OR. DISNAN( SUM ) ) VALUE = SUM
           !$OMP END CRITICAL
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!*
!*        Find normI(A).
         !*

         !$OMP SIMD ALIGNED(WORK:64) LINEAR(I:1)
         DO 50 I = 1, M
            WORK( I ) = ZERO
50       CONTINUE

         !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(WORK,A,N,M) PRIVATE(J,I)            
         DO 70 J = 1, N
            !$OMP SIMD  ALIGNED(WORK:64,A) LINEAR(I:1)
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
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
            CALL ZLASSQ( M, A( 1, J ), 1, COLSSQ( 1 ), COLSSQ( 2 ) )
            CALL DCOMBSSQ( SSQ, COLSSQ )
   90    CONTINUE
         VALUE = SSQ( 1 )*SQRT( SSQ( 2 ) )
      END IF

      ZLANGE = VALUE
END FUNCTION

#if 0
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of elements to be used from the vector X.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension (1+(N-1)*INCX)
*>          The vector x as described above.
*>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between successive values of the vector X.
*>          INCX > 0.
*> \endverbatim
*>
*> \param[in,out] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION
*>          On entry, the value  scale  in the equation above.
*>          On exit, SCALE is overwritten with the value  scl .
*> \endverbatim
*>
*> \param[in,out] SUMSQ
*> \verbatim
*>          SUMSQ is DOUBLE PRECISION
*>          On entry, the value  sumsq  in the equation above.
*>          On exit, SUMSQ is overwritten with the value  ssq .
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
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ ) !GCC$ ATTRIBUTES INLINE :: ZLASSQ !GCC$ ATTRIBUTES aligned(32) :: ZLASSQ
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZLASSQ(N,X,INCX,SCALE,SUMSQ)
    !DIR$ ATTRIBUTES FORCEINLINE :: ZLASSQ
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZLASSQ
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
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
!*     ..
!*     .. Array Arguments ..
      COMPLEX*16         X( * )
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
      DOUBLE PRECISION   TEMP1,TMP
!*     ..
!*     .. External Functions ..
      LOGICAL            DISNAN
      EXTERNAL           DISNAN
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG
!*     ..
!*     .. Executable Statements ..
!*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            TEMP1 = ABS( DBLE( X( IX ) ) )
            IF( TEMP1.GT.ZERO.OR.DISNAN( TEMP1 ) ) THEN
               IF( SCALE.LT.TEMP1 ) THEN
                  TMP = SCALE/TEMP1
                  !SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SUMSQ = 1+SUMSQ*(TMP*TMP)
                  SCALE = TEMP1
               ELSE
                  TMP = TEMP1/SCALE
                  !SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
                  SUMSQ = SUMSQ+TMP*TMP
               END IF
            END IF
            TEMP1 = ABS( DIMAG( X( IX ) ) )
            IF( TEMP1.GT.ZERO.OR.DISNAN( TEMP1 ) ) THEN
               IF( SCALE.LT.TEMP1 ) THEN
                  TMP = SCALE/TEMP1
                  !SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SUMSQ = 1+SUMSQ*(TMP*TMP)
                  SCALE = TEMP1
               ELSE
                  TMP = TEMP1/SCALE
                  SUMSQ = SUMSQ + TMP*TMP
               END IF
            END IF
   10    CONTINUE
      END IF

END SUBROUTINE

#if 0
*  Arguments:
*  ==========
*
*> \param[in,out] V1
*> \verbatim
*>          V1 is DOUBLE PRECISION array, dimension (2).
*>          The first scaled sum.
*>          V1(1) = V1_scale, V1(2) = V1_sumsq.
*> \endverbatim
*>
*> \param[in] V2
*> \verbatim
*>          V2 is DOUBLE PRECISION array, dimension (2).
*>          The second scaled sum.
*>          V2(1) = V2_scale, V2(2) = V2_sumsq.
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
*> \date November 2018
*
*> \ingroup OTHERauxiliary
*
*  =====================================================================
#endif
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


#if 0
*> ZLACPY copies all or part of a two-dimensional matrix A to another
*> matrix B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies the part of the matrix A to be copied to B.
*>          = 'U':      Upper triangular part
*>          = 'L':      Lower triangular part
*>          Otherwise:  All of the matrix A
*> \endverbatim
*>
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
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The m by n matrix A.  If UPLO = 'U', only the upper trapezium
*>          is accessed; if UPLO = 'L', only the lower trapezium is
*>          accessed.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,N)
*>          On exit, B = A in the locations specified by UPLO.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,M).
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
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZLACPY(UPLO, M, N, A, LDA, B, LDB) !GCC$ ATTRIBUTES hot :: ZLACPY !GCC$ ATTRIBUTES aligned(32) :: ZLACPY !GCC$ ATTRIBUTES no_stack_protector :: ZLACPY
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZLACPY(UPLO,M,N,A,LDA,B,LDB)
    !DIR$ ATTRIBUTES INLINE :: ZLACPY
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZLACPY
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZLACPY
#endif
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
    !*     .. Scalar Arguments ..

      use omp_lib

      implicit none
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!*     ..
      !*     .. Array Arguments ..
  
      !COMPLEX*16         A( LDA, * ), B( LDB, * )
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: B
!#elif defined(__ICC) || defined(__INTEL_COMPILER)
!      COMPLEX*16         A( LDA, * ), B( LDB, * )
!      !DIR$ ASSUME_ALIGNED A:64
!      !DIR$ ASSUME_ALIGNED B:64
!#endif
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
!*     .. Executable Statements ..
!*
      IF( LSAME( UPLO, 'U' ) ) THEN

         !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(B,A,N,M) PRIVATE(J,I)
         DO 20 J = 1, N
            !$OMP SIMD ALIGNED(B:64,A:64) LINEAR(I:1)
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
!*
        ELSE IF( LSAME( UPLO, 'L' ) ) THEN

         !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(B,A,N,M) PRIVATE(J,I)
         DO 40 J = 1, N
            !$OMP SIMD ALIGNED(B:64,A:64) LINEAR(I:1)
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
!*
        ELSE

         !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(B,A,N,M) PRIVATE(J,I)
         DO 60 J = 1, N
            !$OMP SIMD ALIGNED(B:64,A) LINEAR(I:1)
           DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF

END SUBROUTINE

#if 0
*  Arguments:
*  ==========
*
*> \param[in] FORWRD
*> \verbatim
*>          FORWRD is LOGICAL
*>          = .TRUE., forward permutation
*>          = .FALSE., backward permutation
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix X. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix X. N >= 0.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension (LDX,N)
*>          On entry, the M by N matrix X.
*>          On exit, X contains the permuted matrix X.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X, LDX >= MAX(1,M).
*> \endverbatim
*>
*> \param[in,out] K
*> \verbatim
*>          K is INTEGER array, dimension (N)
*>          On entry, K contains the permutation vector. K is used as
*>          internal workspace, but reset to its original value on
*>          output.
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
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZLAPMT(FORWRD, M, N, X, LDX, K) !GCC$ ATTRIBUTES hot :: ZLAPMT !GCC$ ATTRIBUTES aligned(32) :: ZLAPMT !GCC$ ATTRIBUTES no_stack_protector :: ZLAPMT
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZLAPMT(FORWRD, M, N, X, LDX, K)
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZLAPMT
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: ZLAPMT
#endif
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
    !*
      implicit none
!*     .. Scalar Arguments ..
      LOGICAL            FORWRD
      INTEGER            LDX, M, N
!*     ..
!*     .. Array Arguments ..
      INTEGER            K( * )
      COMPLEX*16         X( LDX, * )
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER            I, II, IN, J
      COMPLEX*16         TEMP
!*     ..
!*     .. Executable Statements ..
!*
      IF( N.LE.1 ) &
         RETURN
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
!*
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

END SUBROUTINE

#if 0
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the elementary reflector.
*> \endverbatim
*>
*> \param[in,out] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>          On entry, the value alpha.
*>          On exit, it is overwritten with the value beta.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension
*>                         (1+(N-2)*abs(INCX))
*>          On entry, the vector x.
*>          On exit, it is overwritten with the vector v.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between elements of X. INCX > 0.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX*16
*>          The value tau.
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
*> \date November 2017
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
#endif

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU ) !GCC$ ATTRIBUTES hot :: ZLARFG !GCC$ ATTRIBUTES aligned(32) :: ZLARFG !GCC$ ATTRIBUTES no_stack_protector :: ZLARFG
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZLARFG
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: ZLARFG
#endif
!*
!*  -- LAPACK auxiliary routine (version 3.8.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
        !*
        implicit none
!*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA, TAU
!*     ..
!*     .. Array Arguments ..
      COMPLEX*16         X( * )
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
      DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
!*     ..
!*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY3, DZNRM2
      COMPLEX*16         ZLADIV
      EXTERNAL           DLAMCH, DLAPY3, DZNRM2, ZLADIV
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, SIGN
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           ZDSCAL, ZSCAL
!*     ..
!*     .. Executable Statements ..
!*
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
!*
      XNORM = DZNRM2( N-1, X, INCX )
      ALPHR = DBLE( ALPHA )
      ALPHI = DIMAG( ALPHA )
!*
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
!*
!*        H  =  I
!*
         TAU = ZERO
      ELSE
!*
!*        general case
!*
         BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN
!*
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!*
!*           XNORM, BETA may be inaccurate; scale X and recompute them
!*
   10       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) ) &
                GO TO 10
!*
!*           New BETA is at most 1, at least SAFMIN
!*
            XNORM = DZNRM2( N-1, X, INCX )
            ALPHA = DCMPLX( ALPHR, ALPHI )
            BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         END IF
         TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
         ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
         CALL ZSCAL( N-1, ALPHA, X, INCX )
!*
!*        If ALPHA is subnormal, it may lose relative accuracy
!*
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF

END SUBROUTINE 


#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZSCAL(N,ZA,ZX,INCX) !GCC$ ATTRIBUTES INLINE :: ZSCAL !GCC$ ATTRIBUTES aligned(32) :: ZSCAL
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZSCAL(N,ZA,ZX,INCX)
    !DIR$ ATTRIBUTES FORCEINLINE :: ZSCAL
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZSCAL
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZSCAL
#endif

    use omp_lib

     implicit none
!*
!*  -- Reference BLAS level1 routine (version 3.8.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
!*
!*     .. Scalar Arguments ..
      COMPLEX*16 ZA
      INTEGER INCX,N
!*     ..
      !*     .. Array Arguments ..

      !COMPLEX*16 ZX(*)
      COMPLEX(16), DIMENSION(:), ALLOCATABLE :: ZX
!#elif defined(__ICC) || defined(__INTEL_COMPILER)
!      COMPLEX(16) ZX(*)
!      !DIR$ ASSUME_ALIGNED ZX:64
!#endif
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER I,NINCX
!*     ..
      !IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
!*
!*        code for increment equal to 1
         !*

         !$OMP SIMD ALIGNED(ZX:64) LINEAR(I:1) UNROLL PARTIAL(6)
         DO I = 1,N
            ZX(I) = ZA*ZX(I)
         END DO
      ELSE
!*
!*        code for increment not equal to 1
!*
         NINCX = N*INCX

         !$OMP SIMD ALIGNED(ZX:64)
         DO I = 1,NINCX,INCX
            ZX(I) = ZA*ZX(I)
         END DO
      END IF
    
END SUBROUTINE


#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZDSCAL(N,DA,ZX,INCX) !GCC$ ATTRIBUTES INLINE :: ZDSCAL !GCC$ ATTRIBUTES aligned(32) :: ZDSCAL
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZDSCAL(N,DA,ZX,INCX)
  !DIR$ ATTRIBUTES FORCEINLINE :: ZDSCAL
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZDSCAL
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZDSCAL
#endif
!*
!*  -- Reference BLAS level1 routine (version 3.8.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
!*
    !*     .. Scalar Arguments ..

      use omp_lib
      implicit none    
      DOUBLE PRECISION DA
      INTEGER INCX,N
!*     ..
      !*     .. Array Arguments ..
    
      !COMPLEX*16 ZX(*)
      COMPLEX(16), DIMENSION(:), ALLOCATABLE :: ZX

!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER I,NINCX
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DCMPLX
!*     ..
!      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
!*
!*        code for increment equal to 1
         !*

         !$OMP SIMD ALIGNED(ZX:64) LINEAR(I:1) UNROLL PARTIAL(6)
         DO I = 1,N
            ZX(I) = DCMPLX(DA,0.0d0)*ZX(I)
         END DO
      ELSE
!*
!*        code for increment not equal to 1
!*
         NINCX = N*INCX

         !$OMP SIMD ALIGNED(ZX:64) UNROLL PARTIAL(6)
         DO I = 1,NINCX,INCX
            ZX(I) = DCMPLX(DA,0.0d0)*ZX(I)
         END DO
      END IF
    
END SUBROUTINE

    
! Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16OTHERauxiliary
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZLASET(UPLO, M, N, ALPHA, BETA, A, LDA) !GCC$ ATTRIBUTES hot :: ZLASET !GCC$ ATTRIBUTES aligned(32) :: ZLASET !GCC$ ATTRIBUTES no_stack_protector :: ZLASET
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZLASET(UPLO, M, N, ALPHA, BETA, A, LDA)
  !DIR$ ATTRIBUTES FORCEINLINE :: ZLASET
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZLASET
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZLASET
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
      COMPLEX*16         ALPHA, BETA
!*     ..
      !*     .. Array Arguments ..

      !COMPLEX*16         A( LDA, * )
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
!#elif defined(__ICC) || defined(__INTEL_COMPILER)
!      COMPLEX(16), DIMENSION(LDA,*) :: A
!      !DIR$ ASSUME_ALIGNED A:64
!#endif
!*     ..
!!*
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
!*     .. Executable Statements ..
!*
      IF( LSAME( UPLO, 'U' ) ) THEN
!*
!*        Set the diagonal to BETA and the strictly upper triangular
!*        part of the array to ALPHA.
         !*

         !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(A,ALPHA,N,M) PRIVATE(J,I)
         DO 20 J = 2, N
            !$OMP SIMD ALIGNED(A:64) LINEAR(I:1)
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
         DO 30 I = 1, MIN( N, M )
            A( I, I ) = BETA
   30    CONTINUE
!*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!*
!*        Set the diagonal to BETA and the strictly lower triangular
!*        part of the array to ALPHA.
         !*

         !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(A,ALPHA,M,N) PRIVATE(J,I)
         DO 50 J = 1, MIN( M, N )
            !$OMP SIMD ALIGNED(A:64) LINEAR(I:1)
            DO 40 I = J + 1, M
               A( I, J ) = ALPHA
   40       CONTINUE
   50    CONTINUE
         DO 60 I = 1, MIN( N, M )
            A( I, I ) = BETA
   60    CONTINUE
!*
      ELSE
!*
!*        Set the array to BETA on the diagonal and ALPHA on the
!*        offdiagonal.
         !*

         !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(A,ALPHA,N,M) PRIVATE(J,I)
         DO 80 J = 1, N
            !$OMP SIMD ALIGNED(A:64) LINEAR(I:1)
            DO 70 I = 1, M
               A( I, J ) = ALPHA
   70       CONTINUE
   80    CONTINUE
         DO 90 I = 1, MIN( M, N )
            A( I, I ) = BETA
   90    CONTINUE
      END IF

END SUBROUTINE

#if 0
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': form P * C
*>          = 'R': form C * P
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension
*>                  (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*>                  (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*>          The vector v in the representation of P. V is not used
*>          if TAU = 0.
*> \endverbatim
*>
*> \param[in] INCV
*> \verbatim
*>          INCV is INTEGER
*>          The increment between elements of v. INCV <> 0
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16
*>          The value tau in the representation of P.
*> \endverbatim
*>
*> \param[in,out] C1
*> \verbatim
*>          C1 is COMPLEX*16 array, dimension
*>                         (LDC,N) if SIDE = 'L'
*>                         (M,1)   if SIDE = 'R'
*>          On entry, the n-vector C1 if SIDE = 'L', or the m-vector C1
*>          if SIDE = 'R'.
*>
*>          On exit, the first row of P*C if SIDE = 'L', or the first
*>          column of C*P if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in,out] C2
*> \verbatim
*>          C2 is COMPLEX*16 array, dimension
*>                         (LDC, N)   if SIDE = 'L'
*>                         (LDC, N-1) if SIDE = 'R'
*>          On entry, the (m - 1) x n matrix C2 if SIDE = 'L', or the
*>          m x (n - 1) matrix C2 if SIDE = 'R'.
*>
*>          On exit, rows 2:m of P*C if SIDE = 'L', or columns 2:m of C*P
*>          if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the arrays C1 and C2.
*>          LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension
*>                      (N) if SIDE = 'L'
*>                      (M) if SIDE = 'R'
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
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZLATZM(SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK) !GCC$ ATTRIBUTES hot :: ZLATZM !GCC$ ATTRIBUTES aligned(32) :: ZLATZM !GCC$ ATTRIBUTES no_stack_protector :: ZLATZM
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE ZLATZM(SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK)
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZLATZM
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: ZLATZM
#endif
!*
!*  -- LAPACK computational routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
  !*     .. Scalar Arguments ..
      implicit none
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
!*     ..
!*     .. Array Arguments ..
      !COMPLEX*16         C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * )
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: C1
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: C2
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: V
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: WORK
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
                         ZERO = ( 0.0D+0, 0.0D+0 ) )
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           ZAXPY, ZCOPY, ZGEMV, ZGERC, ZGERU, ZLACGV
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
!*        w :=  ( C1 + v**H * C2 )**H
!*
         CALL ZCOPY( N, C1, LDC, WORK, 1 )
         CALL ZLACGV( N, WORK, 1 )
         CALL ZGEMV( 'Conjugate transpose', M-1, N, ONE, C2, LDC, V, &
                    INCV, ONE, WORK, 1 )
!*
!*        [ C1 ] := [ C1 ] - tau* [ 1 ] * w**H
!*        [ C2 ]    [ C2 ]        [ v ]
!*
         CALL ZLACGV( N, WORK, 1 )
         CALL ZAXPY( N, -TAU, WORK, 1, C1, LDC )
         CALL ZGERU( M-1, N, -TAU, V, INCV, WORK, 1, C2, LDC )
!*
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!*
!*        w := C1 + C2 * v
!*
         CALL ZCOPY( M, C1, 1, WORK, 1 )
         CALL ZGEMV( 'No transpose', M, N-1, ONE, C2, LDC, V, INCV, ONE, &
                     WORK, 1 )
!*
!*        [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**H]
!*
         CALL ZAXPY( M, -TAU, WORK, 1, C1, 1 )
         CALL ZGERC( M, N-1, -TAU, WORK, 1, V, INCV, C2, LDC )
      END IF
!*
END SUBROUTINE

#if 0
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 4/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY) !GCC$ ATTRIBUTES inline :: ZCOPY !GCC$ ATTRIBUTES aligned(32) :: ZCOPY
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
    !DIR$ ATTRIBUTES FORCEINLINE :: ZCOPY
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZCOPY
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZCOPY

#endif
        use omp_lib
       implicit none
!*
!*  -- Reference BLAS level1 routine (version 3.8.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
!*
!*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!*     ..
!*     .. Array Arguments ..
      ! COMPLEX*16 ZX(*),ZY(*)
      COMPLEX(16), DIMENSION(:), ALLOCATABLE :: ZX
      COMPLEX(16), DIMENSION(:), ALLOCATABLE :: ZY
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER I,IX,IY
!*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!*
!*        code for both increments equal to 1
         !*
         !$OMP SIMD ALIGNED(ZY:64,ZX) LINEAR(I:1) UNROLL PARTIAL(8)
         DO I = 1,N
          ZY(I) = ZX(I)
         END DO
      ELSE
!*
!*        code for unequal increments or equal increments
!*          not equal to 1
!*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         !$OMP SIMD ALIGNED(ZY:64,ZY)  LINEAR(I:1)
         DO I = 1,N
            ZY(IY) = ZX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF     
END SUBROUTINE

#if 0
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY) !GCC$ ATTRIBUTES INLINE :: ZAXPY !GCC$ ATTRIBUTES aligned(32) :: ZAXPY
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
 !DIR$ ATTRIBUTES FORCEINLINE :: ZCOPY
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZCOPY
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZCOPY
#endif
       use omp_lib
       implicit none
!*
!*  -- Reference BLAS level1 routine (version 3.8.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
!*
!*     .. Scalar Arguments ..
      COMPLEX*16 ZA
      INTEGER INCX,INCY,N
!*     ..
!*     .. Array Arguments ..
      !      COMPLEX*16 ZX(*),ZY(*)
      COMPLEX(16), DIMENSION(:), ALLOCATABLE :: ZX
      COMPLEX(16), DIMENSION(:), ALLOCATABLE :: ZY
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER I,IX,IY
!*     ..
!*     .. External Functions ..
      !DOUBLE PRECISION DCABS1
      !EXTERNAL DCABS1
!*     ..
      
      IF (DCABS1(ZA).EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!*
!*        code for both increments equal to 1
         !*
         !$OMP SIMD ALIGNED(ZY:64,ZX) LINEAR(I:1) UNROLL PARTIAL(10)
         DO I = 1,N
            ZY(I) = ZY(I) + ZA*ZX(I)
         END DO
      ELSE
!*
!*        code for unequal increments or equal increments
!*          not equal to 1
!*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         !$OMP SIMD ALIGNED(ZY:64,ZX) LINEAR(I:1) UNROLL PARTIAL(10)
         DO I = 1,N
            ZY(IY) = ZY(IY) + ZA*ZX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
END SUBROUTINE
   
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2017
!*
!*> \ingroup double_blas_level1
!*
!*  =====================================================================
      DOUBLE PRECISION FUNCTION DCABS1(Z)
!*
!*  -- Reference BLAS level1 routine (version 3.8.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
!*
!*     .. Scalar Arguments ..
      COMPLEX*16 Z
!*     ..
!*     ..
!*  =====================================================================
!*
!*     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG
      DCABS1 = ABS(DBLE(Z)) + ABS(DIMAG(Z))
END FUNCTION     

#if 0
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
*> \ingroup complex16_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) !GCC$ ATTRIBUTES hot :: ZGEMV !GCC$ ATTRIBUTES aligned(32) :: ZGEMV !GCC$ ATTRIBUTES no_stack_protector :: ZGEMV
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZGEMV
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZGEMV
#endif
      use omp_lib
      implicit none
!*
!*  -- Reference BLAS level2 routine (version 3.7.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
!*     ..
!*     .. Array Arguments ..
      ! COMPLEX*16 A(LDA,*),X(*),Y(*)
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: X
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: Y
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
!*     ..
!*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL NOCONJ
!*     ..
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
!*     ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
!!          CALL XERBLA('ZGEMV ',INFO)
          RETURN
      END IF
!*
!*     Quick return if possible.
!*
!      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
!     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!*
      NOCONJ = LSAME(TRANS,'T')
!*
!*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!*     up the start points in  X  and  Y.
!*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
!*
!*     Start the operations. In this version the elements of A are
!*     accessed sequentially with one pass through A.
!*
!*     First form  y := beta*y.
!*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
             IF (BETA.EQ.ZERO) THEN

                !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(10)
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
             ELSE
                  !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1)   UNROLL PARTIAL(6)  
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN

                !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1)
                DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1)    
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
!*
!*        Form  y := alpha*A*x + y.
!*
          JX = KX
          IF (INCY.EQ.1) THEN

             !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(Y,A,N,M,ALPHA,N,M) FIRSTPRIVATE(JX) PRIVATE(J,TEMP,I)
              DO 60 J = 1,N
                 TEMP = ALPHA*X(JX)
                  !$OMP SIMD ALIGNED(Y:64,A) LINEAR(I:1) UNROLL PARTIAL(10)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE

             !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(X,Y,A,ALPHA,N,M) FIRSTPRIVATE(JX) PRIVATE(J,TEMP,IY,I)
             DO 80 J = 1,N
                TEMP = ALPHA*X(JX)
                
                   IY = KY
               
                  !$OMP SIMD ALIGNED(Y:64,A) LINEAR(I:1) UNROLL PARTIAL(10)
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
!*
!*        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
!*
          JY = KY
          IF (INCX.EQ.1) THEN

             !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(A,X,Y,N,M,NOCONJ,ALPHA,INCY) REDUCTION(+:TEMP) FIRSTPRIVATE(JY) PRIVATE(J,I)
             DO 110 J = 1,N
                  !$OMP SINGLE
                    TEMP = ZERO
                  !$OMP END SINGLE
                  IF (NOCONJ) THEN
                      !$OMP SIMD ALIGNED(Y:64,A) UNROLL PARTIAL(10)
                      DO 90 I = 1,M
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                  ELSE
                      !$OMP SIMD ALIGNED(Y:64,A) LINEAR(I:1)   UNROLL PARTIAL(10)    
                      DO 100 I = 1,M
                          TEMP = TEMP + DCONJG(A(I,J))*X(I)
  100                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  110         CONTINUE
           ELSE

             !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(A,X,Y,NOCONJ,N,M,ALPHA) REDUCTION(+:TEMP) FIRSTPRIVATE(JY) PRIVATE(J,IX,I,JY)
              DO 140 J = 1,N
                
                  TEMP = ZERO
                  IX = KX
               
                  IF (NOCONJ) THEN
                       !$OMP SIMD ALIGNED(Y:64,A)  UNROLL PARTIAL(10)
                      DO 120 I = 1,M
                          TEMP = TEMP + A(I,J)*X(IX)
                          IX = IX + INCX
  120                 CONTINUE
                  ELSE
                      !$OMP SIMD ALIGNED(Y:64,A)  UNROLL PARTIAL(10)    
                      DO 130 I = 1,M
                          TEMP = TEMP + DCONJG(A(I,J))*X(IX)
                          IX = IX + INCX
  130                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  140         CONTINUE
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
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
!*>     Richard Hanson, Sandia National Labs.
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA) !GCC$ ATTRIBUTES hot :: ZGERC !GCC$ ATTRIBUTES aligned(32) :: ZGERC !GCC$ ATTRIBUTES no_stack_protector :: ZGERC
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    !DIR$ ATTRIBUTES FORCEINLINE :: ZGERC
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZGERC
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZGERC
#endif
!*
!*  -- Reference BLAS level2 routine (version 3.7.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
    !*     .. Scalar Arguments ..
      use omp_lib
      implicit none
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,LDA,M,N
!*     ..
!*     .. Array Arguments ..
      ! COMPLEX*16 A(LDA,*),X(*),Y(*)
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: X
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: Y
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
!*     ..
!*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JY,KX
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
!*     ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
       END IF
       IF (INFO.NE.0) THEN
!          CALL XERBLA('ZGERC ',INFO)
          RETURN
      END IF
!*
!*     Quick return if possible. 
!*
!      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!*
!*     Start the operations. In this version the elements of A are
!*     accessed sequentially with one pass through A.
!*
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN

         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,4) SHARED(Y,A,X,N,M,ALPHA,INCY)  FIRSTPRIVATE(JY) PRIVATE(J,TEMP,I)
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                 TEMP = ALPHA*DCONJG(Y(JY))
                  !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF

         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,4) SHARED(Y,A,X,N,M,ALPHA,KX) FIRSTPRIVATE(JY) PRIVATE(J,TEMP,I,IX)
         DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*DCONJG(Y(JY))
                  IX = KX
                  !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
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
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
!*>     Richard Hanson, Sandia National Labs.
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA) !GCC$ ATTRIBUTES hot :: ZGERU !GCC$ ATTRIBUTES aligned(32) :: ZGERU !GCC$ ATTRIBUTES no_stack_protector
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
 !DIR$ ATTRIBUTES FORCEINLINE :: ZGERU
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZGERU
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZGERU
#endif
!*
!*  -- Reference BLAS level2 routine (version 3.7.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
  !*     .. Scalar Arguments ..
      use omp_lib
      implicit none
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,LDA,M,N
!*     ..
!*     .. Array Arguments ..
      ! COMPLEX*16 A(LDA,*),X(*),Y(*)
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: X
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: Y
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
!*     ..
!*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JY,KX
!*     ..
!*     .. External Subroutines ..
      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MAX
!*     ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
!!          CALL XERBLA('ZGERU ',INFO)
          RETURN
      END IF
!*
!*     Quick return if possible.
!*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!*
!*     Start the operations. In this version the elements of A are
!*     accessed sequentially with one pass through A. 
!*
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN

         !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(Y,A,X,N,M,ALPHA,INCY) FIRSTPRIVATE(JY) PRIVATE(J,TEMP,I)
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                 TEMP = ALPHA*Y(JY)
                  !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF

         !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(Y,A,X,N,M,ALPHA,INCX,KX) FIRSTPRIVATE(JY) PRIVATE(J,TEMP,I,IX)
         DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                 TEMP = ALPHA*Y(JY)
                  IX = KX
                  !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
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
!*> \ingroup complex16OTHERauxiliary
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZLACGV(N, X, INCX) !GCC$ ATTRIBUTES inline :: ZLACGV !GCC$ ATTRIBUTES aligned(32) :: ZLACGV
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZLACGV(N,X,INCX)
 !DIR$ ATTRIBUTES FORCEINLINE :: ZLACGV
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZLACGV
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: ZLACGV
#endif
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER            INCX, N
!*     ..
!*     .. Array Arguments ..
      COMPLEX*16         X( * )
!*     ..
!*
!* =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER            I, IOFF
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!*     ..
!*     .. Executable Statements ..
!*
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = DCONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 ) &
            IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = DCONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
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
!*> \ingroup complex16OTHERcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZUNMQR(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
     WORK, LWORK, INFO) !GCC$ ATTRIBUTES hot :: ZUNMQR !GCC$ ATTRIBUTES aligned(32) :: ZUNMQR !GCC$ ATTRIBUTES no_stack_protector :: ZUNMQR
#elif defined(__ICC) || defined(__INTEL_COMPILER)
 SUBROUTINE ZUNMQR(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
      WORK, LWORK, INFO)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZUNMQR
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: ZUNMQR
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
      !COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: C
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: TAU
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: WORK
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
!      EXTERNAL           ZLARFB, ZLARFT, ZUNM2R
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
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
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
         NB = MIN( NBMAX, ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M, N, K, &
              -1 ) )
         LWKOPT = MAX( 1, NW )*NB + TSIZE
         WORK( 1 ) = LWKOPT
      END IF
!*
      IF( INFO.NE.0 ) THEN
!         CALL XERBLA( 'ZUNMQR', -INFO )
         RETURN
      IF( LQUERY ) THEN
         RETURN
      END IF
!*
!*     Quick return if possible
!*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF

      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IF( LWORK.LT.NW*NB+TSIZE ) THEN
            NB = (LWORK-TSIZE) / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZUNMQR', SIDE // TRANS, M, N, K, &
                    -1 ) )
         END IF
      END IF
!*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!*
!*        Use unblocked code
!*
         CALL ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
                      IINFO )
      ELSE
!*
!*        Use blocked code
!*
         IWT = 1 + NW*NB
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.   &
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
            CALL ZLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), &
                         LDA, TAU( I ), WORK( IWT ), LDT )
            IF( LEFT ) THEN
!*
!*              H or H**H is applied to C(i:m,1:n)
!*
               MI = M - I + 1
               IC = I
            ELSE
!*
!*              H or H**H is applied to C(1:m,i:n)
!*
               NI = N - I + 1
               JC = I
            END IF
!*
!*           Apply H or H**H
!*
            CALL ZLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, &
                         IB, A( I, I ), LDA, WORK( IWT ), LDT,         &
                         C( IC, JC ), LDC, WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
END SUBROUTINE
    

!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date June 2013
!*
!*> \ingroup complex16OTHERauxiliary
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
SUBROUTINE ZLARFB(SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
  T, LDT, C, LDC, WORK, LDWORK) !GCC$ ATTRIBUTES hot :: ZLARFB !GCC$ ATTRIBUTES aligned(32) :: ZLARFB !GCC$ ATTRIBUTES no_stack_protector :: ZLARFB
#elif defined(__ICC) || defined(__INTEL_COMPILER)
 SUBROUTINE ZLARFB(SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
   T, LDT, C, LDC, WORK, LDWORK)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZLARFB
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZLARFB
#endif
       use omp_lib
       implicit none
!      *
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
      !COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ), &
      !     WORK( LDWORK, * )
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: C
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: T
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: V
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: WORK
      
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
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
!      EXTERNAL           ZCOPY, ZGEMM, ZLACGV, ZTRMM
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!*     ..
!*     .. Executable Statements ..
!*
!*     Quick return if possible
!*
!      IF( M.LE.0 .OR. N.LE.0 )
!     $   RETURN
!*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'C'
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
!*              Form  H * C  or  H**H * C  where  C = ( C1 )
!*                                                    ( C2 )
!*
!*              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!*
!*              W := C1**H
!*
               DO 10 J = 1, K
                  CALL ZCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
   10          CONTINUE
!*
!*              W := W * V1
!*
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!*
!*                 W := W + C2**H * V2
!*
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', N, &
                             K, M-K, ONE, C( K+1, 1 ), LDC,            &
                             V( K+1, 1 ), LDV, ONE, WORK, LDWORK )
               END IF
!*
!*              W := W * T**H  or  W * T
!*
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - V * W**H
!*
               IF( M.GT.K ) THEN
!*
!*                 C2 := C2 - V2 * W**H
!*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose',  &
                             M-K, N, K, -ONE, V( K+1, 1 ), LDV, WORK, &
                             LDWORK, ONE, C( K+1, 1 ), LDC )
               END IF
!*
!*              W := W * V1**H
!*
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',  &
                           'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
!*
!*              C1 := C1 - W**H
               !*

               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,N) PRIVATE(J,I) COLLAPSE(2)
               DO 30 J = 1, K
                 
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
   20             CONTINUE
   30          CONTINUE
!*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!*
!*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!*
!*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!*
!*              W := C1
!*
               DO 40 J = 1, K
                  CALL ZCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!*
!*              W := W * V1
!*
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                          K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!*
!*                 W := W + C2 * V2
!*
                  CALL ZGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                             ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,    &
                             ONE, WORK, LDWORK )
               END IF
!*
!*              W := W * T  or  W * T**H
!*
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - W * V**H
!*
               IF( N.GT.K ) THEN
!*
!*                 C2 := C2 - W * V2**H
!*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M,  &
                             N-K, K, -ONE, WORK, LDWORK, V( K+1, 1 ),    &
                              LDV, ONE, C( 1, K+1 ), LDC )
               END IF
!*
!*              W := W * V1**H
!*
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',   &
                           'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
!*
!*              C1 := C1 - W
               !*

               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,M) COLLAPSE(2) PRIVATE(J,I)
            DO 60 J = 1, K
                  !$OMP SIMD ALIGNED(C:64,WORK) LINEAR(I:1) UNROLL PARTIAL(6)
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF

         ELSE
!*
!*           Let  V =  ( V1 )
!*                     ( V2 )    (last K rows)
!*           where  V2  is unit upper triangular.
!*
            IF( LSAME( SIDE, 'L' ) ) THEN
!*
!*              Form  H * C  or  H**H * C  where  C = ( C1 )
!*                                                    ( C2 )
!*
!*              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!*
!*              W := C2**H
!*
               DO 70 J = 1, K
                  CALL ZCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
   70          CONTINUE
!*
!*              W := W * V2
!*
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!*
!*                 W := W + C1**H * V1
!*
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', N,  &
                            K, M-K, ONE, C, LDC, V, LDV, ONE, WORK,      &
                             LDWORK )
               END IF
!*
!*              W := W * T**H  or  W * T
!*
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,  &
                           ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - V * W**H
!*
               IF( M.GT.K ) THEN
!*
!*                 C1 := C1 - V1 * W**H
!*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose',  &
                             M-K, N, K, -ONE, V, LDV, WORK, LDWORK,   &
                             ONE, C, LDC )
               END IF
!*
!*              W := W * V2**H
!*
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',     & 
                          'Unit', N, K, ONE, V( M-K+1, 1 ), LDV, WORK,  &
                          LDWORK )
!*
!*              C2 := C2 - W**H
               !*

               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,M) PRIVATE(J,I) COLLAPSE(2)
              DO 90 J = 1, K
                 
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) -   &
                                     DCONJG( WORK( I, J ) )
   80             CONTINUE
   90          CONTINUE

            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!*
!*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!*
!*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!*
!*              W := C2
!*
               DO 100 J = 1, K
                  CALL ZCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
!*
!*              W := W * V2
!*
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,  &
                           K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!*
!*                 W := W + C1 * V1
!*
                  CALL ZGEMM( 'No transpose', 'No transpose', M, K, N-K,  &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!*
!*              W := W * T  or  W * T**H
!*
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                          ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - W * V**H
!*
               IF( N.GT.K ) THEN
!*
!*                 C1 := C1 - W * V1**H
!*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, &
                             N-K, K, -ONE, WORK, LDWORK, V, LDV, ONE,   &
                             C, LDC )
               END IF
!*
!*              W := W * V2**H
!*
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',   &
                          'Unit', M, K, ONE, V( N-K+1, 1 ), LDV, WORK, &
                           LDWORK )
!*
!*              C2 := C2 - W
               !*

               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,M) PRIVATE(J,I) COLLAPSE(2)
               DO 120 J = 1, K
                  !$OMP SIMD ALIGNED(C:64,WORK) LINEAR(I:1) UNROLL PARTIAL(6)
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF

      ELSE IF( LSAME( STOREV, 'R' ) ) THEN

         IF( LSAME( DIRECT, 'F' ) ) THEN
!*
!*           Let  V =  ( V1  V2 )    (V1: first K columns)
!*           where  V1  is unit upper triangular.
!*
            IF( LSAME( SIDE, 'L' ) ) THEN
!*
!*              Form  H * C  or  H**H * C  where  C = ( C1 )
!*                                                    ( C2 )
!*
!*              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!*
!*              W := C1**H
!*
               DO 130 J = 1, K
                  CALL ZCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
  130          CONTINUE
!*
!*              W := W * V1**H
!*
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',  &
                           'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!*
!*                 W := W + C2**H * V2**H
!*
                  CALL ZGEMM( 'Conjugate transpose',  &
                             'Conjugate transpose', N, K, M-K, ONE, & 
                              C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, & 
                              WORK, LDWORK )
               END IF
!*
!*              W := W * T**H  or  W * T
!*
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,  &
                          ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - V**H * W**H
!*
               IF( M.GT.K ) THEN
!*
!*                 C2 := C2 - V2**H * W**H
!*
                  CALL ZGEMM( 'Conjugate transpose',  &
                             'Conjugate transpose', M-K, N, K, -ONE, &
                             V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
                             C( K+1, 1 ), LDC )
               END IF
!*
!*              W := W * V1
!*
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
!*
!*              C1 := C1 - W**H
               !*

               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,M) PRIVATE(J,I) COLLAPSE(2)
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
  140             CONTINUE
  150          CONTINUE
!*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!*
!*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!*
!*              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!*
!*              W := C1
!*
               DO 160 J = 1, K
                  CALL ZCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!*
!*              W := W * V1**H
!*
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',  & 
                           'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!*
!*                 W := W + C2 * V2**H
!*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M, &
                             K, N-K, ONE, C( 1, K+1 ), LDC,             &
                             V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
               END IF
!*
!*              W := W * T  or  W * T**H
!*
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - W * V
!*
               IF( N.GT.K ) THEN
!*
!*                 C2 := C2 - W * V2
!*
                  CALL ZGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                             -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,  &
                              C( 1, K+1 ), LDC )
               END IF
!*
!*              W := W * V1
!*
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                            K, ONE, V, LDV, WORK, LDWORK )
!*
!*              C1 := C1 - W
               !*

               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,M) COLLAPSE(2) PRIVATE(J,I)
               DO 180 J = 1, K
                  !$OMP SIMD ALIGNED(C:64,WORK) LINEAR(I:1) UNROLL PARTIAL(6)
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE

            END IF

         ELSE
!*
!*           Let  V =  ( V1  V2 )    (V2: last K columns)
!*           where  V2  is unit lower triangular.
!*
            IF( LSAME( SIDE, 'L' ) ) THEN
!*
!*              Form  H * C  or  H**H * C  where  C = ( C1 )
!*                                                    ( C2 )
!*
!*              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!*
!*              W := C2**H
!*
               DO 190 J = 1, K
                  CALL ZCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
  190          CONTINUE
!*
!*              W := W * V2**H
!*
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',     &
                          'Unit', N, K, ONE, V( 1, M-K+1 ), LDV, WORK,  &
                          LDWORK )
               IF( M.GT.K ) THEN
!*
!*                 W := W + C1**H * V1**H
!*
                  CALL ZGEMM( 'Conjugate transpose',  &
                             'Conjugate transpose', N, K, M-K, ONE, C, &
                              LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!*
!*              W := W * T**H  or  W * T
!*
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,  &
                           ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - V**H * W**H
!*
               IF( M.GT.K ) THEN
!*
!*                 C1 := C1 - V1**H * W**H
!*
                  CALL ZGEMM( 'Conjugate transpose',   &
                             'Conjugate transpose', M-K, N, K, -ONE, V, &
                              LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!*
!*              W := W * V2
!*
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,  &
                           K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!*
!*              C2 := C2 - W**H
               !*

               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,M) PRIVATE(J,I) COLLAPSE(2)
               DO 210 J = 1, K
                 DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) -   &
                                    DCONJG( WORK( I, J ) )
  200             CONTINUE
  210          CONTINUE
!*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!*
!*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!*
!*              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!*
!*              W := C2
!*
               DO 220 J = 1, K
                  CALL ZCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
!*
!*              W := W * V2**H
!*
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',  &
                         'Unit', M, K, ONE, V( 1, N-K+1 ), LDV, WORK, &
                           LDWORK )
               IF( N.GT.K ) THEN
!*
!*                 W := W + C1 * V1**H
!*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M,  &
                             K, N-K, ONE, C, LDC, V, LDV, ONE, WORK,  &
                             LDWORK )
               END IF
!*
!*              W := W * T  or  W * T**H
!*
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,  &
                           ONE, T, LDT, WORK, LDWORK )
!*
!*              C := C - W * V
!*
               IF( N.GT.K ) THEN
!*
!*                 C1 := C1 - W * V1
!*
                  CALL ZGEMM( 'No transpose', 'No transpose', M, N-K, K,  &
                             -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!*
!*              W := W * V2
!*
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,  &
                          K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!*
!*              C1 := C1 - W
               !*

               !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(C,WORK,K,M) PRIVATE(J,I) COLLAPSE(2)
               DO 240 J = 1, K
                  !$OMP SIMD ALIGNED(C:64,WORK)
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE

            END IF

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
!*> \ingroup complex16_blas_level3
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 3 Blas routine.
!*>
!*>  -- Written on 8-February-1989.
!*>     Jack Dongarra, Argonne National Laboratory.
!*>     Iain Duff, AERE Harwell.
!*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!*>     Sven Hammarling, Numerical Algorithms Group Ltd.
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) !GCC$ ATTRIBUTES hot :: ZTRMM !GCC$ ATTRIBUTES aligned(32) :: ZTRMM !GCC$ ATTRIBUTES no_stack_protector :: ZTRMM
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZTRMM
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZTRMM
#endif
       use omp_lib
       implicit none
!*
!*  -- Reference BLAS level3 routine (version 3.7.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!*     ..
!*     .. Array Arguments ..
      !COMPLEX*16 A(LDA,*),B(LDB,*)
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: B
!*     ..
!*
!*  =====================================================================
!*
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
!*     ..
!*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
!*     ..
!*     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
!*     ..
!*
!*     Test the input parameters.
!*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOCONJ = LSAME(TRANSA,'T')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
!*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
              (.NOT.LSAME(TRANSA,'T')) .AND. &
              (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
!          CALL XERBLA('ZTRMM ',INFO)
          RETURN
      END IF
!*
!*     Quick return if possible.
!*
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
!*
!!*     And when  alpha.eq.zero.
!*
      IF (ALPHA.EQ.ZERO) THEN
         DO 20 J = 1,N
            !$OMP SIMD ALIGNED(B:64) UNROLL PARTIAL(8) LINEAR(I:1)
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
!*
!*     Start the operations.
!*
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
!*
!*           Form  B := alpha*A*B.
!*
             IF (UPPER) THEN

                !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(B,A,N,M,NOUNIT,ALPHA) PRIVATE(J,K,I,TEMP)
                  DO 50 J = 1,N
                      DO 40 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                             TEMP = ALPHA*B(K,J)
                             !$OMP SIMD ALIGNED(B:64,A) LINEAR(I:1) UNROLL PARTIAL(10)
                             DO 30 I = 1,K - 1
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   30                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP*A(K,K)
                              B(K,J) = TEMP
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE

                !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(B,A,ALPHA,NOUNIT,N,M) PRIVATE(J,K,I,TEMP)
                DO 80 J = 1,N
                      DO 70 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              B(K,J) = TEMP
                              IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                              !$OMP SIMD ALIGNED(B:64,A) LINEAR(I:1) UNROLL PARTIAL(10)
                              DO 60 I = K + 1,M
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
!*
!*           Form  B := alpha*A**T*B   or   B := alpha*A**H*B.
!*
             IF (UPPER) THEN

                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(B,A,N,M,NOCONJ,NOUNIT,ALPHA) PRIVATE(J,K,I,TEMP)
                DO 120 J = 1,N
                      DO 110 I = M,1,-1
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                             IF (NOUNIT) TEMP = TEMP*A(I,I)
                                !$OMP SIMD ALIGNED(B:64,A) LINEAR(I:1) REDUCTION(+:TEMP)
                              DO 90 K = 1,I - 1
                                  TEMP = TEMP + A(K,I)*B(K,J)
   90                         CONTINUE
                          ELSE
                             IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                              !$OMP SIMD ALIGNED(B:64,A) LINEAR(I:1) REDUCTION(+:TEMP)
                              DO 100 K = 1,I - 1
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  100                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  110                 CONTINUE
  120             CONTINUE
             ELSE

                !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(B,A,N,M,NOCONJ,NOUNIT) PRIVATE(J,K,I,TEMP)
                DO 160 J = 1,N
                      DO 150 I = 1,M
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                             IF (NOUNIT) TEMP = TEMP*A(I,I)
                              !$OMP SIMD ALIGNED(B:64,A) LINEAR(I:1) REDUCTION(+:TEMP)
                              DO 130 K = I + 1,M
                                  TEMP = TEMP + A(K,I)*B(K,J)
  130                         CONTINUE
                          ELSE
                             IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                               !$OMP SIMD ALIGNED(B:64,A) LINEAR(I:1) REDUCTION(+:TEMP)
                              DO 140 K = I + 1,M
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  140                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
!*
!*           Form  B := alpha*B*A.
!*
             IF (UPPER) THEN

                !$OMP PARALLEL DO  SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(A,B,N,M,ALPHA)  PRIVATE(J,TEMP,I,K)
                DO 200 J = N,1,-1
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      !$OMP SIMD ALIGNED(B:64) LINEAR(I:1) UNROLL PARTIAL(6)
                      DO 170 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  170                 CONTINUE
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                             TEMP = ALPHA*A(K,J)
                              !$OMP SIMD ALIGNED(B:64) LINEAR(I:1) UNROLL PARTIAL(6)
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
  200             CONTINUE
               ELSE

                !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(A,B,N,ALPHA,NOUNIT) PRIVATE(J,I,TEMP,K)
                DO 240 J = 1,N
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                       !$OMP SIMD ALIGNED(B:64) LINEAR(I:1) UNROLL PARTIAL(6)
                      DO 210 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  210                 CONTINUE
                      DO 230 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                             TEMP = ALPHA*A(K,J)
                               !$OMP SIMD ALIGNED(B:64) LINEAR(I:1) UNROLL PARTIAL(6)
                              DO 220 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  220                         CONTINUE
                          END IF
  230                 CONTINUE
  240             CONTINUE
              END IF
          ELSE
!*
!*           Form  B := alpha*B*A**T   or   B := alpha*B*A**H.
!*
             IF (UPPER) THEN

                !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE)  SHARED(A,B,N,M,NOCONJ,ALPHA,NOUNIT) PRIVATE(K,J,TEMP,I)
                  DO 280 K = 1,N
                      DO 260 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              END IF
                              !$OMP SIMD ALIGNED(B:64) LINEAR(I:1) 
                              DO 250 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*DCONJG(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                           !$OMP SIMD ALIGNED(B:64) LINEAR(I:1) 
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              ELSE
                    !$OMP PARALLEL DO SCHEDULE(STATIC,1) DEFAULT(NONE)  SHARED(A,B) PRIVATE(K,J,TEMP,I)    
                  DO 320 K = N,1,-1
                      DO 300 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              END IF
                                !$OMP SIMD ALIGNED(B:64) LINEAR(I:1) UNROLL PARTIAL(10)
                              DO 290 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  290                         CONTINUE
                          END IF
  300                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*DCONJG(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                          !$OMP SIMD ALIGNED(B:64) LINEAR(I:1) UNROLL PARTIAL(6)
                          DO 310 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  310                     CONTINUE
                      END IF
  320             CONTINUE
              END IF
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
!*> \ingroup complex16OTHERcomputational
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZUNM2R(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
     WORK, INFO) !GCC$ ATTRIBUTES hot :: ZUNM2R !GCC$ ATTRIBUTES aligned(32) :: ZUNM2R !GCC$ ATTRIBUTES no_stack_protector :: ZUNM2R
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE ZUNM2R(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
     WORK, INFO)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZUNM2R
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: ZUNM2R
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
      !COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: C
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: TAU
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: WORK
!!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      COMPLEX*16         AII, TAUI
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           ZLARF
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
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
!      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
!         INFO = -1
!      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
!         INFO = -2
!      ELSE IF( M.LT.0 ) THEN
!         INFO = -3
!      ELSE IF( N.LT.0 ) THEN
!         INFO = -4
!      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
!         INFO = -5
!      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
!         INFO = -7
!      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
!         INFO = -10
!      END IF
!!      IF( INFO.NE.0 ) THEN
!         CALL XERBLA( 'ZUNM2R', -INFO )
!         RETURN
!      END IF
!*
!*     Quick return if possible
!*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
         RETURN

      IF( ( LEFT .AND. .NOT.NOTRAN .OR. .NOT.LEFT .AND. NOTRAN ) ) THEN
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
!*           H(i) or H(i)**H is applied to C(i:m,1:n)
!*
            MI = M - I + 1
            IC = I
         ELSE
!*
!*           H(i) or H(i)**H is applied to C(1:m,i:n)
!*
            NI = N - I + 1
            JC = I
         END IF
!*
!*        Apply H(i) or H(i)**H
!*
         IF( NOTRAN ) THEN
            TAUI = TAU( I )
         ELSE
            TAUI = DCONJG( TAU( I ) )
         END IF
         AII = A( I, I )
         A( I, I ) = ONE
         CALL ZLARF( SIDE, MI, NI, A( I, I ), 1, TAUI, C( IC, JC ), LDC, &
                     WORK )
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
!*> \date December 2016
!*
!!*> \ingroup complex16OTHERauxiliary
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZLARF(SIDE, M, N, V, INCV, TAU, C, LDC, WORK) !GCC$ ATTRIBUTES hot :: ZLARF !GCC$ ATTRIBUTES aligned(32) :: ZLARF !GCC$ ATTRIBUTES no_stack_protector :: ZLARF
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZLARF(SIDE, M, N, V, INCV, TAU, C, LDC, WORK)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZLARF
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: ZLARF
#endif
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
    !*     .. Scalar Arguments ..
      implicit none
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
!*     ..
!*     .. Array Arguments ..
      !COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: C
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: V
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: WORK
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),  &
                         ZERO = ( 0.0D+0, 0.0D+0 ) )
!*     ..
!*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL           ZGEMV, ZGERC
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAZLR, ILAZLC
      EXTERNAL           LSAME, ILAZLR, ILAZLC
!*     ..
!*     .. Executable Statements ..
!*
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!*     Set up variables for scanning V.  LASTV begins pointing to the end
!*     of V.
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
!*     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
!*     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILAZLC(LASTV, N, C, LDC)
         ELSE
!*     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILAZLR(M, LASTV, C, LDC)
         END IF
      END IF
!*     Note that lastc.eq.0 renders the BLAS operations null; no special
!*     case is needed at this level.
      IF( APPLYLEFT ) THEN
!*
!*        Form  H * C
!*
         IF( LASTV.GT.0 ) THEN
!*
!*           w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)
!*
            CALL ZGEMV( 'Conjugate transpose', LASTV, LASTC, ONE, &
                 C, LDC, V, INCV, ZERO, WORK, 1 )
!*
!*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H
!*
            CALL ZGERC( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!*
!*        Form  C * H
!*
         IF( LASTV.GT.0 ) THEN
!*
!*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!*
            CALL ZGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC, &
                 V, INCV, ZERO, WORK, 1 )
!*
!*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H
!*
            CALL ZGERC( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
    
END SUBROUTINE


!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date June 2016
!*
!*> \ingroup complex16OTHERauxiliary
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
SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT) !GCC$ ATTRIBUTES hot :: ZLARFT !GCC$ ATTRIBUTES aligned(32) :: ZLARFT !GCC$ ATTRIBUTES no_stack_protector :: ZLARFT
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZLARFT
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZLARFT
#endif
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     June 2016
!*
    !*     .. Scalar Arguments ..
    use omp_lib
      implicit none
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!*     ..
!*     .. Array Arguments ..
      !COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: T
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: TAU
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: V
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
                         ZERO = ( 0.0D+0, 0.0D+0 ) )
!*     ..
!*     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
!*     ..
!*     .. External Subroutines ..
 !     EXTERNAL           ZGEMV, ZTRMV, ZGEMM
!*     ..
!*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     ..
!*     .. Executable Statements ..
!*
!*     Quick return if possible
!*
      IF( N.EQ.0 ) &
         RETURN

      IF( LSAME( DIRECT, 'F' ) ) THEN
         PREVLASTV = N
         DO I = 1, K
            PREVLASTV = MAX( PREVLASTV, I )
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
                  !$OMP SIMD ALIGNED(T:64,TAU,V) LINEAR(J:1) UNROLL PARTIAL(6)
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * CONJG( V( I , J ) )
                  END DO
                  J = MIN( LASTV, PREVLASTV )
!*
!*                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i)
!*
                  CALL ZGEMV( 'Conjugate transpose', J-I, I-1, &
                             -TAU( I ), V( I+1, 1 ), LDV, &
                             V( I+1, I ), 1, ONE, T( 1, I ), 1 )
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
!*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H
!*
                  CALL ZGEMM( 'N', 'C', I-1, 1, J-I, -TAU( I ), &
                            V( 1, I+1 ), LDV, V( I, I+1 ), LDV, &
                             ONE, T( 1, I ), LDT )
               END IF
!*
!*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!*
               CALL ZTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
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
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * CONJG( V( N-K+I , J ) )
                     END DO
                     J = MAX( LASTV, PREVLASTV )
!*
!*                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i)
!*
                     CALL ZGEMV( 'Conjugate transpose', N-K+I-J, K-I,  &
                                -TAU( I ), V( J, I+1 ), LDV, V( J, I ),
                                1, ONE, T( I+1, I ), 1 )
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
!*                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H
!*
                     CALL ZGEMM( 'N', 'C', K-I, 1, N-K+I-J, -TAU( I ), &
                                V( I+1, J ), LDV, V( I, J ), LDV, &
                                ONE, T( I+1, I ), LDT )
                  END IF
!*
!*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!*
                  CALL ZTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
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
!*> \ingroup complex16_blas_level3
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 3 Blas routine.
!*>
!*>  -- Written on 8-February-1989.
!*>     Jack Dongarra, Argonne National Laboratory.
!*>     Iain Duff, AERE Harwell.
!*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!*>     Sven Hammarling, Numerical Algorithms Group Ltd.
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) !GCC$ ATTRIBUTES hot :: ZGEMM !GCC$ ATTRIBUTES aligned(32) :: ZGEMM !GCC$ ATTRIBUTES no_stack_protector :: ZGEMM
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZGEMM
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZGEMM
#endif
       use omp_lib
       implicit none
!*
!*  -- Reference BLAS level3 routine (version 3.7.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
!*     ..
!*     .. Array Arguments ..
      ! COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: B
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: C
!*     ..
!*
!*  =====================================================================
!*
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
!*     ..
!*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL CONJA,CONJB,NOTA,NOTB
!*     ..
!*     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
!*     ..
!*
!*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!*!     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!*     B  respectively are to be  transposed but  not conjugated  and set
!*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
!*     and the number of rows of  B  respectively.
!*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      CONJA = LSAME(TRANSA,'C')
      CONJB = LSAME(TRANSB,'C')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
!*
!*     Test the input parameters.
!*
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.CONJA) .AND. &
         (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.CONJB) .AND. &
              (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
!          CALL XERBLA('ZGEMM ',INFO)
          RETURN
      END IF
!*
!*     Quick return if possible.
!*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
         (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
!*
!*     And when  alpha.eq.zero.
!*
      IF (ALPHA.EQ.ZERO) THEN
         IF (BETA.EQ.ZERO) THEN
            DO 20 J = 1,N

               !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(8)
                DO 10 I = 1,M
                     C(I,J) = ZERO
   10             CONTINUE
  20         CONTINUE
         ELSE

        DO 40 J = 1,N
           !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)
          DO 30 I = 1,M
               C(I,J) = BETA*C(I,J)
   30      CONTINUE
   40   CONTINUE

        END IF
          RETURN
      END IF
!*
!*     Start the operations.
!*
      IF (NOTB) THEN
          IF (NOTA) THEN
!*
!*           Form  C := alpha*A*B + beta*C.
             !*
             !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,4) SHARED(C,A,N,BETA,M,K,ALPHA) PRIVATE(J,I,L,TEMP) 
              DO 90 J = 1,N
                 IF (BETA.EQ.ZERO) THEN
                      !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(8)
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)    
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                     TEMP = ALPHA*B(L,J)
                      !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)
                      DO 70 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
            !$OMP END PARALLEL DO
                         
          ELSE IF (CONJA) THEN
!*
!*           Form  C := alpha*A**H*B + beta*C.
             !*

             !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,4) SHARED(A,B,C,N,M,K,BETA,ALPHA) PRIVATE(J,I,TEMP,L)
              DO 120 J = 1,N
                  DO 110 I = 1,M
                     TEMP = ZERO
                      !$OMP SIMD ALIGNED(A:64,B) REDUCTION(+:TEMP) LINEAR(L:1)
                      DO 100 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)  
                          C(I,J) = ALPHA*TEMP
                      ELSE
                           !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)  
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120           CONTINUE
            !$OMP END PARALLEL DO
                     
          ELSE
!*
!*           Form  C := alpha*A**T*B + beta*C
             !*

             !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,4) SHARED(A,B,C,N,M,K,BETA,ALPHA) PRIVATE(J,I,TEMP,L) 
              DO 150 J = 1,N
                  DO 140 I = 1,M
                     TEMP = ZERO
                      !$OMP SIMD ALIGNED(A:64,B) LINEAR(L:1) REDUCTION(+:TEMP)
                      DO 130 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  130                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                            !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)     
                          C(I,J) = ALPHA*TEMP
                       ELSE
                            !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)  
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  140             CONTINUE
  150           CONTINUE
            !$OMP END PARALLEL DO
                      
          END IF
      ELSE IF (NOTA) THEN
          IF (CONJB) THEN
!*
!*           Form  C := alpha*A*B**H + beta*C.
             !*

             !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,4) SHARED(C,B,A,BETA,ALPHA,N,M,K) PRIVATE(J,I,L,TEMP) 
             DO 200 J = 1,N
                 IF (BETA.EQ.ZERO) THEN
                       !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(8)    
                      DO 160 I = 1,M
                          C(I,J) = ZERO
  160                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                       !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)        
                      DO 170 I = 1,M
                          C(I,J) = BETA*C(I,J)
  170                 CONTINUE
                  END IF
                  DO 190 L = 1,K
                     TEMP = ALPHA*DCONJG(B(J,L))
                       !$OMP SIMD ALIGNED(C:64,A) LINEAR(I:1) UNROLL PARTIAL(10)    
                      DO 180 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  180                 CONTINUE
  190             CONTINUE
  200           CONTINUE
            !$OMP END PARALLEL DO
                   
          ELSE
!*
!*           Form  C := alpha*A*B**T + beta*C
             !*

             !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,4) SHARED(C,B,A,N,M,ALPHA,BETA) PRIVATE(J,I,L,TEMP) 
              DO 250 J = 1,N
                 IF (BETA.EQ.ZERO) THEN
                      !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(8)    
                      DO 210 I = 1,M
                          C(I,J) = ZERO
  210                 CONTINUE
                 ELSE IF (BETA.NE.ONE) THEN
                        !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)        
                      DO 220 I = 1,M
                          C(I,J) = BETA*C(I,J)
  220                 CONTINUE
                  END IF
                  DO 240 L = 1,K
                     TEMP = ALPHA*B(J,L)
                        !$OMP SIMD ALIGNED(C:64,A) LINEAR(I:1) UNROLL PARTIAL(10)    
                      DO 230 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  230                 CONTINUE
  240             CONTINUE
  250          CONTINUE
            !$OMP END PARALLEL DO
                         
          END IF
      ELSE IF (CONJA) THEN
          IF (CONJB) THEN
!*
!*           Form  C := alpha*A**H*B**H + beta*C.
             !*

             !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,4) SHARED(A,B,C,N,M,K,BETA,ZERO,ALPHA)
             !$OMP PRIVATE(J,I,L) REDUCTION(+:TEMP) COLLAPSE(2)
             DO 280 J = 1,N
                  DO 270 I = 1,M
                     TEMP = ZERO
                       !$OMP SIMD ALIGNED(A:64,B) LINEAR(L:1) 
                      DO 260 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*DCONJG(B(J,L))
  260                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                           !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(10)       
                          C(I,J) = ALPHA*TEMP
                       ELSE
                           !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(10)    
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  270             CONTINUE
  280          CONTINUE
            !$OMP END PARALLEL DO
                    
                      
          ELSE
!*
!*           Form  C := alpha*A**H*B**T + beta*C
             !*

             !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,4) SHARED(A,B,C,N,M,K,BETA,ZERO)
             !$OMP PRIVATE(J,I,L) REDUCTION(+:TEMP) COLLAPSE(2)
             DO 310 J = 1,N
                  DO 300 I = 1,M
                     TEMP = ZERO
                       !$OMP SIMD ALIGNED(A:64,B) LINEAR(L:1) 
                      DO 290 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*B(J,L)
  290                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)       
                          C(I,J) = ALPHA*TEMP
                       ELSE
                          !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(10)       
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  300             CONTINUE
  310           CONTINUE
            !$OMP END PARALLEL DO
                         
          END IF
      ELSE
          IF (CONJB) THEN
!*
!*           Form  C := alpha*A**T*B**H + beta*C
             !*

             !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,4) SHARED(A,B,C,N,M,K,BETA,ZERO) REDUCTION(+:TEMP)
             !$OMP PRIVATE(J,I,L) COLLAPSE(2)
             DO 340 J = 1,N
                  DO 330 I = 1,M
                     TEMP = ZERO
                        !$OMP SIMD ALIGNED(A:64,B) LINEAR(L:1) 
                      DO 320 L = 1,K
                          TEMP = TEMP + A(L,I)*DCONJG(B(J,L))
  320                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                             !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)        
                          C(I,J) = ALPHA*TEMP
                       ELSE
                            !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(10)       
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  330             CONTINUE
  340          CONTINUE
            !$OMP END PARALLEL DO
                       
          ELSE
!*
!*           Form  C := alpha*A**T*B**T + beta*C
             !*

             !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(GUIDED,1) SHARED(A,B,C,N,M,K,BETA,ZERO) PRIVATE(J,I,L)
             !$OMP REDUCTION(+:TEMP) COLLAPSE(2)
             DO 370 J = 1,N
                  DO 360 I = 1,M
                     TEMP = ZERO
                        !$OMP SIMD ALIGNED(A:64,B) LINEAR(L:1)
                      DO 350 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  350                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                               !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)       
                          C(I,J) = ALPHA*TEMP
                       ELSE
                           !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(10)       
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  360             CONTINUE
  370          CONTINUE
            !$OMP END PARALLEL DO
                           
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
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
!*>     Richard Hanson, Sandia National Labs.
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) !GCC$ ATTRIBUTES hot :: ZTRMV !GCC$ ATTRIBUTES aligned(32) :: ZTRMV !GCC$ no_stack_protector :: ZTRMV
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ZTRMV
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ZTRMV
#endif
       use omp_lib
       implicit none
!*
!*  -- Reference BLAS level2 routine (version 3.7.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
!*     ..
!*     .. Array Arguments ..
      ! COMPLEX*16 A(LDA,*),X(*)
      COMPLEX(16), DIMENSION(:,:), ALLOCATABLE :: A
      COMPLEX(16), DIMENSION(:),   ALLOCATABLE :: X
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
!*     ..
!*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
!*     ..
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
!*     ..
!*
!*     Test the input parameters.
!*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
!          CALL XERBLA('ZTRMV ',INFO)
          RETURN
      END IF
!*
!*     Quick return if possible.
!*
      IF (N.EQ.0) RETURN
!*
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
!*
!*     Set up the start point in X if the increment is not unity. This
!*     will be  ( N - 1 )*INCX  too small for descending loops.
!*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
!*
!*     Start the operations. In this version the elements of A are
!*     accessed sequentially with one pass through A.
!*
      IF (LSAME(TRANS,'N')) THEN
!*
!*        Form  x := A*x.
!*
          IF (LSAME(UPLO,'U')) THEN
             IF (INCX.EQ.1) THEN

                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(X,A) PRIVATE(J,TEMP,I)
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                         TEMP = X(J)
                          !$OMP SIMD ALIGNED(X:64,A) LINEAR(I:1) UNROLL PARTIAL(10)
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*A(I,J)
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   20             CONTINUE
                      !$OMP END PARALLEL DO

              ELSE
                 JX = KX

                !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(X,A,INCX,N) FIRSTPRIVATE(JX) PRIVATE(J,TEMP,IX,I)
                DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                           !$OMP SIMD ALIGNED(X:64,A) LINEAR(I:1) UNROLL PARTIAL(10)
                          DO 30 I = 1,J - 1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX + INCX
   40              CONTINUE
                      !$OMP END PARALLEL DO
                      
              END IF
          ELSE
             IF (INCX.EQ.1) THEN

                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(X,A,NOUNIT) PRIVATE(J,TEMP,I)
                DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                         TEMP = X(J)
                          !$OMP SIMD ALIGNED(X:64,A) LINEAR(I:1) UNROLL PARTIAL(10)
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*A(I,J)
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   60             CONTINUE
                      !$OMP END PARALLEL DO
                        
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX

                !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(A,X,N,NOUNIT,KX) FIRSTPRIVATE(JX) PRIVATE(J,TEMP,IX,I) 
                DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                            !$OMP SIMD ALIGNED(X:64,A) LINEAR(I:1) UNROLL PARTIAL(10)
                          DO 70 I = N,J + 1,-1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX - INCX
   80             CONTINUE
                      !$OMP END PARALLEL DO
                         
              END IF
          END IF
      ELSE
!*
!*        Form  x := A**T*x  or  x := A**H*x.
!*
          IF (LSAME(UPLO,'U')) THEN
             IF (INCX.EQ.1) THEN

                !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE)  SHARED(X,A,N,NOCONJ,NOUNIT) PRIVATE(J,I) REDUCTION(+:TEMP)
                DO 110 J = N,1,-1
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                         IF (NOUNIT) TEMP = TEMP*A(J,J)
                          !$OMP SIMD ALIGNED(A:64,X) 
                          DO 90 I = J - 1,1,-1
                              TEMP = TEMP + A(I,J)*X(I)
   90                     CONTINUE
                      ELSE
                         IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                            !$OMP SIMD ALIGNED(A:64,X) 
                          DO 100 I = J - 1,1,-1
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
  100                     CONTINUE
                      END IF
                      X(J) = TEMP
  110              CONTINUE
                      !$OMP END PARALLEL DO
                       
              ELSE
                 JX = KX + (N-1)*INCX

                !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(X,A,N,NOCONJ,NOUNIT) FIRSTPRIVATE(JX) PRIVATE(J,IX,I) REDUCTION(+:TEMP) 
                DO 140 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                         IF (NOUNIT) TEMP = TEMP*A(J,J)
                          !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1)
                          DO 120 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  120                     CONTINUE
                      ELSE
                         IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                           !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1)
                          DO 130 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
  130                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
140               CONTINUE
                      !$OMP END PARALLEL DO
                       
              END IF
          ELSE
             IF (INCX.EQ.1) THEN

                !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) SHARED(X,A,N) PRIVATE(J,I) REDUCTION(+:TEMP)
                DO 170 J = 1,N
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                         IF (NOUNIT) TEMP = TEMP*A(J,J)
                            !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) 
                          DO 150 I = J + 1,N
                              TEMP = TEMP + A(I,J)*X(I)
  150                     CONTINUE
                      ELSE
                         IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                            !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) 
                          DO 160 I = J + 1,N
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
  160                     CONTINUE
                      END IF
                      X(J) = TEMP
  170              CONTINUE
                      !$OMP END PARALLEL DO
                        
              ELSE
                 JX = KX

                !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(NONE) REDUCTION(+:TEMP)  SHARED(X,A,N) FIRSTPRIVATE(JX) PRIVATE(J,IX,I)
                DO 200 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                         IF (NOUNIT) TEMP = TEMP*A(J,J)
                           !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1)
                          DO 180 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  180                     CONTINUE
                      ELSE
                         IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                            !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) 
                          DO 190 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
  190                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
  200             CONTINUE
                      !$OMP END PARALLEL DO
                      
              END IF
          END IF
      END IF

END SUBROUTINE
