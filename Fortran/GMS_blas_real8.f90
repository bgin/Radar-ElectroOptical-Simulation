


!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2017
!*
!*> \ingroup double_blas_level1
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>     jack dongarra, linpack, 3/11/78.
!*>     modified 3/93 to return if incx .le. 0.
!*>     modified 12/3/93, array(1) declarations changed to array(*)
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX) !GCC$ ATTRIBUTES inline :: DASUM !GCC$ ATTRIBUTES aligned(32) :: DASUM
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
    !DIR$ ATTRIBUTES FORCEINLINE :: DASUM
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DASUM
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DASUM
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
      INTEGER INCX,N
!*     ..
!*     .. Array Arguments ..
      !      DOUBLE PRECISION DX(*)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DX
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,M,MP1,NINCX
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DABS,MOD
!*     ..
      DASUM = 0.0d0
      DTEMP = 0.0d0
      !IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
!*        code for increment equal to 1
!*
!!*
!*        clean-up loop
!*
         M = MOD(N,6)
         IF (M.NE.0) THEN
            !$OMP SIMD ALIGNED(DX:64) LINEAR(I:1) REDUCTION(+:DTEMP) UNROLL PARTIAL(6)
            DO I = 1,M
               DTEMP = DTEMP + DABS(DX(I))
            END DO
            IF (N.LT.6) THEN
               DASUM = DTEMP
               RETURN
            END IF
         END IF
         MP1 = M + 1
          !$OMP SIMD ALIGNED(DX:64) LINEAR(I:6) REDUCTION(+:DTEMP) 
         DO I = MP1,N,6
            DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) + &
                    DABS(DX(I+2)) + DABS(DX(I+3)) + &
                    DABS(DX(I+4)) + DABS(DX(I+5))
         END DO
      ELSE
!*
!*        code for increment not equal to 1
!*
         NINCX = N*INCX
         !$OMP SIMD ALIGNED(DX:64) REDUCTION(+:DTEMP) UNROLL PARTIAL(6)
         DO I = 1,NINCX,INCX
            DTEMP = DTEMP + DABS(DX(I))
         END DO
      END IF
      DASUM = DTEMP
     
END FUNCTION 

    
!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2017
!*
!*> \ingroup double_blas_level1
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>     jack dongarra, linpack, 3/11/78.
!*>     modified 12/3/93, array(1) declarations changed to array(*)
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY) !GCC$ ATTRIBUTES inline :: DAXPY !GCC$ ATTRIBUTES aligned(32) :: DAXPY
#elif defined(__INTEL_COMPILER) || defined(__ICC)
SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
   !DIR$ ATTRIBUTES FORCEINLINE :: DAXPY
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DAXPY
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DAXPY
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
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION DX(*),DY(*)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DX
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DY
!!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MOD
!*     ..
!      IF (N.LE.0) RETURN
!      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!*
!*        code for both increments equal to 1
!*
!*
!*        clean-up loop
!*
         M = MOD(N,4)
         IF (M.NE.0) THEN
              DO I = 1,M
               DY(I) = DY(I) + DA*DX(I)
            END DO
         END IF
         IF (N.LT.4) RETURN
         MP1 = M + 1
         !$OMP SIMD ALIGNED(DY,DX) LINEAR(I:4) 
         DO I = MP1,N,4
            DY(I) = DY(I) + DA*DX(I)
            DY(I+1) = DY(I+1) + DA*DX(I+1)
            DY(I+2) = DY(I+2) + DA*DX(I+2)
            DY(I+3) = DY(I+3) + DA*DX(I+3)
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
         DO I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
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
!*> \date November 2017
!*
!*> \ingroup double_blas_level1
!*
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
DOUBLE PRECISION FUNCTION DCABS1(Z) !GCC$ ATTRIBUTES inline :: DCABS1 !GCC$ ATTRIBUTES aligned(32) :: DCABS1
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  DOUBLE PRECISION FUNCTION DCABS1(Z)
 !DIR$ ATTRIBUTES FORCEINLINE :: DCABS1
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DCABS1
    !DIR$ OPTIMIZE : 3
#endif
       implicit none
!*
!*  -- Reference BLAS level1 routine (version 3.8.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
!*
      DCABS1 = ABS(DBLE(Z)) + ABS(DIMAG(Z))
END FUNCTION

    
!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2017
!*
!*> \ingroup double_blas_level1
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>     jack dongarra, linpack, 3/11/78.
!*>     modified 12/3/93, array(1) declarations changed to array(*)
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DCOPY(N,DX,INCX,DY,INCY) !GCC$ ATTRIBUTES inline :: DCOPY !GCC$ ATTRIBUTES aligned(32) :: DCOPY
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
 !DIR$ ATTRIBUTES FORCEINLINE :: DCOPY
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DCOPY
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DCOPY
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
      !DOUBLE PRECISION DX(*),DY(*)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DX
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DY
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MOD
!*     ..
!      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!*
!*        code for both increments equal to 1
!*
!*
!*        clean-up loop
!*
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF
         MP1 = M + 1
         !$OMP SIMD ALIGNED(DY:64,DX) LINEAR(I:7)
         DO I = MP1,N,7
            DY(I) = DX(I)
            DY(I+1) = DX(I+1)
            DY(I+2) = DX(I+2)
            DY(I+3) = DX(I+3)
            DY(I+4) = DX(I+4)
            DY(I+5) = DX(I+5)
            DY(I+6) = DX(I+6)
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
         !$OMP SIMD ALIGNED(DY:64,DX) UNROLL PARTIAL(8)
         DO I = 1,N
            DY(IY) = DX(IX)
            IX = IX + INCX
            IY = IY + INCY
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
!*> \date November 2017
!*
!*> \ingroup double_blas_level1
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>     jack dongarra, linpack, 3/11/78.
!*>     modified 12/3/93, array(1) declarations changed to array(*)
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY) !GCC$ ATTRIBUTES inline :: DDOT !GCC$ ATTRIBUTES aligned(32) :: DDOT
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
  !DIR$ ATTRIBUTES FORCEINLINE :: DDOT
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DDOT
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DDOT
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
      !DOUBLE PRECISION DX(*),DY(*)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DX
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DY
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MOD
!*     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      !IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!*
!*        code for both increments equal to 1
!*
!*
!*        clean-up loop
!*
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DX(I)*DY(I)
            END DO
            IF (N.LT.5) THEN
               DDOT=DTEMP
            RETURN
            END IF
         END IF
         MP1 = M + 1
         !$OMP SIMD ALIGNED(DX:64,DY) LINEAR(I:5) REDUCTION(+:DTEMP)
         DO I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + &
                  DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
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
          !$OMP SIMD ALIGNED(DX:64,DY) REDUCTION(+:DTEMP) UNROLL PARTIAL(6)
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      DDOT = DTEMP
      
END FUNCTION


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
!*> \ingroup double_blas_level2
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
SUBROUTINE DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) !GCC$ ATTRIBUTES hot :: DGBMV !GCC$ ATTRIBUTES aligned(32) :: DGBMV !GCC$ ATTRIBUTES no_stack_protector :: DGBMV
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGBMV
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DGBMV
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
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,KL,KU,LDA,M,N
      CHARACTER TRANS
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION A(LDA,*),X(*),Y(*)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: X
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: Y
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,K,KUP1,KX,KY,LENX,LENY
!*     ..
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
!*     ..
!*
!*     Test the input parameters.
!*
!      INFO = 0
 !     IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
   !  +    .NOT.LSAME(TRANS,'C')) THEN
  !!        INFO = 1
   !   ELSE IF (M.LT.0) THEN
   !       INFO = 2
   !   ELSE IF (N.LT.0) THEN
    !      INFO = 3
     ! ELSE IF (KL.LT.0) THEN
     !     INFO = 4
     ! ELSE IF (KU.LT.0) THEN
    !      INFO = 5
     ! ELSE IF (LDA.LT. (KL+KU+1)) THEN
     !     INFO = 8
     ! ELSE IF (INCX.EQ.0) THEN
     !     INFO = 10
     ! ELSE IF (INCY.EQ.0) THEN
     !     INFO = 13
     ! END IF
     ! IF (INFO.NE.0) THEN
      !    CALL XERBLA('DGBMV ',INFO)
      !    RETURN
     ! END IF
!*
!*     Quick return if possible.
!*
!      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
!     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
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
!*     accessed sequentially with one pass through the band part of A.
!*
!*     First form  y := beta*y.
!*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
             IF (BETA.EQ.ZERO) THEN
                  !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(8)
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
             ELSE
                   !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(6)    
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                   !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(8)
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                   !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(6)    
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      KUP1 = KU + 1
      IF (LSAME(TRANS,'N')) THEN
!*
!*        Form  y := alpha*A*x + y.
!*
          JX = KX
          IF (INCY.EQ.1) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,K,JX) IF(N>=400)
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  K = KUP1 - J
                  DO 50 I = MAX(1,J-KU),MIN(M,J+KL)
                      Y(I) = Y(I) + TEMP*A(K+I,J)
   50             CONTINUE
                  JX = JX + INCX
   60          CONTINUE
               !$OMP END PARALLEL DO   
          ELSE
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,K,IY,JX) IF(N>=400)  
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  K = KUP1 - J
                  DO 70 I = MAX(1,J-KU),MIN(M,J+KL)
                      Y(IY) = Y(IY) + TEMP*A(K+I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
                  IF (J.GT.KU) KY = KY + INCY
80            CONTINUE
              !$OMP END PARALLEL DO    
          END IF
      ELSE
!*
!*        Form  y := alpha*A**T*x + y.
!*
          JY = KY
          IF (INCX.EQ.1) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,K,JY) IF(N>=400)  
              DO 100 J = 1,N
                  TEMP = ZERO
                  K = KUP1 - J
                  DO 90 I = MAX(1,J-KU),MIN(M,J+KL)
                      TEMP = TEMP + A(K+I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
100           CONTINUE
              !$OMP END PARALLEL DO
               ELSE
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,IX,K,JY) IF(N>=400)      
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  K = KUP1 - J
                  DO 110 I = MAX(1,J-KU),MIN(M,J+KL)
                      TEMP = TEMP + A(K+I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
                  IF (J.GT.KU) KX = KX + INCX
120               CONTINUE
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
!*> \ingroup double_blas_level3
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
SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) !GCC$ ATTRIBUTES hot :: DGEMM !GCC$ ATTRIBUTES aligned(32) :: DGEMM !GCC$ ATTRIBUTES no_stack_protector :: DGEMM
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGEMM
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DGEMM
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
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
!*     ..
!*
!*  =====================================================================
!*
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
      !EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MAX
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
!*     ..
!*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!*     ..
!*
!*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!*     and  columns of  A  and the  number of  rows  of  B  respectively.
!*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
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
!      INFO = 0
!      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND.
!     +    (.NOT.LSAME(TRANSA,'T'))) THEN
!          INFO = 1
!      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND.
!     +         (.NOT.LSAME(TRANSB,'T'))) THEN
!          INFO = 2
!      ELSE IF (M.LT.0) THEN
!          INFO = 3
!      ELSE IF (N.LT.0) THEN
!          INFO = 4
!      ELSE IF (K.LT.0) THEN
!          INFO = 5
!      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
!!          INFO = 8
!      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
!          INFO = 10
!      ELSE IF (LDC.LT.MAX(1,M)) THEN
!          INFO = 13
!      END IF
!      IF (INFO.NE.0) THEN
!          CALL XERBLA('DGEMM ',INFO)
!          RETURN
!      END IF
!*
!*     Quick return if possible.
!*
!      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
!     +    (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
!*
!*     And if  alpha.eq.zero.
!*
      IF (ALPHA.EQ.ZERO) THEN
         IF (BETA.EQ.ZERO) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J) IF(N>=400)
            DO 20 J = 1,N
                  !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(8)
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
         ELSE
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J) IF(N>=400)        
            DO 40 J = 1,N
                   !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
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
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J) IF(N>=400)       
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
                        !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(10)   
                      DO 70 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
90             CONTINUE
                !$OMP END PARALLEL DO           
          ELSE
!*
!*           Form  C := alpha*A**T*B + beta*C
             !*
                 !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) COLLAPSE(2) PRIVATE(J,I,TEMP) IF(N>=400)   
              DO 120 J = 1,N
                  DO 110 I = 1,M
                     TEMP = ZERO
                        !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) REDUCTION(+:TEMP) UNROLL PARTIAL(10)   
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
120           CONTINUE
                !$OMP END PARALLEL DO
          END IF
      ELSE
          IF (NOTA) THEN
!*
!*           Form  C := alpha*A*B**T + beta*C
             !*
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J) IF(N>=400)      
              DO 170 J = 1,N
                 IF (BETA.EQ.ZERO) THEN
                         !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(8)   
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                 ELSE IF (BETA.NE.ONE) THEN
                         !$OMP SIMD ALIGNED(C:64) LINEAR(I:1)  UNROLL PARTIAL(6)       
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  END IF
                  DO 160 L = 1,K
                     TEMP = ALPHA*B(J,L)
                         !$OMP SIMD ALIGNED(C:64) LINEAR(I:1)  UNROLL PARTIAL(10)   
                      DO 150 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  150                 CONTINUE
  160             CONTINUE
170            CONTINUE
               !$OMP END PARALLEL DO            
          ELSE
!*
!*           Form  C := alpha*A**T*B**T + beta*C
             !*
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) COLLAPSE(2) PRIVATE(J,I,TEMP) IF(N>=400)   
              DO 200 J = 1,N
                  DO 190 I = 1,M
                     TEMP = ZERO
                       !$OMP SIMD ALIGNED(C:64) LINEAR(I:1)  REDUCTION(+:TEMP) UNROLL PARTIAL(10)   
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  190             CONTINUE
200             CONTINUE
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
!*> \ingroup double_blas_level2
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
SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) !GCC$ ATTRIBUTES hot :: DGEMV !GCC$ ATTRIBUTES aligned(32) :: DGEMV !GCC$ ATTRIBUTES no_stack_protector :: DGEMV
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGEMV
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DGEMV
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
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION A(LDA,*),X(*),Y(*)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: X
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: Y
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
!*     ..
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
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
!      INFO = 0
    !  IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
    ! +    .NOT.LSAME(TRANS,'C')) THEN
    !      INFO = 1
    !  ELSE IF (M.LT.0) THEN
    !      INFO = 2
   !   ELSE IF (N.LT.0) THEN
    !      INFO = 3
    !  ELSE IF (LDA.LT.MAX(1,M)) THEN
   !       INFO = 6
   !   ELSE IF (INCX.EQ.0) THEN
   !       INFO = 8
   !   ELSE IF (INCY.EQ.0) THEN
   !       INFO = 11
   !   END IF
   !   IF (INFO.NE.0) THEN
 !!         CALL XERBLA('DGEMV ',INFO)
 !         RETURN
 !     END IF
!*
!*     Quick return if possible.
!*
!      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
!     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
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
                  !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(8)
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
             ELSE
                    !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(6)   
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                   !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(8)
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                    !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(6)    
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      !IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
!*
!*        Form  y := alpha*A*x + y.
!*
          JX = KX
          IF (INCY.EQ.1) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,JX) IF(N>=400)
              DO 60 J = 1,N
                 TEMP = ALPHA*X(JX)
                    !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(10)    
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
60            CONTINUE
              !$OMP END PARALLEL DO    
           ELSE
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,IY,JX) IF(N>=400)   
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                    !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(10)    
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
80            CONTINUE
              !$OMP END PARALLEL DO    
          END IF
      ELSE
!*
!*        Form  y := alpha*A**T*x + y.
!*
          JY = KY
          IF (INCX.EQ.1) THEN
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,JY) IF(N>=400)   
              DO 100 J = 1,N
                 TEMP = ZERO
                   !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) REDUCTION(+:TEMP) UNROLL PARTIAL(10)    
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
100            CONTINUE
                !$OMP END PARALLEL DO       
           ELSE
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,IX,JY) IF(N>=400)      
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                   !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) REDUCTION(+:TEMP) UNROLL PARTIAL(10) 
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
120          CONTINUE
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
!*> \ingroup double_blas_level2
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
SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA) !GCC$ ATTRIBUTES inline :: DGER !GCC$ ATTRIBUTES aligned(32) :: DGER
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    !DIR$ ATTRIBUTES FORCEINLINE :: DGER
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DGER
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DGER
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
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,M,N
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION A(LDA,*),X(*),Y(*)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: X
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: Y
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MAX
!*     ..
!!*
!*     Test the input parameters.
!*
!      INFO = 0
!      IF (M.LT.0) THEN
!          INFO = 1
 !     ELSE IF (N.LT.0) THEN
!          INFO = 2
!      ELSE IF (INCX.EQ.0) THEN
!          INFO = 5
!      ELSE IF (INCY.EQ.0) THEN
!          INFO = 7
!      ELSE IF (LDA.LT.MAX(1,M)) THEN
!          INFO = 9
!      END IF
!      IF (INFO.NE.0) THEN
!          CALL XERBLA('DGER  ',INFO)
!          RETURN
!      END IF
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
          !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,JY) IF(N>=400)
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                 TEMP = ALPHA*Y(JY)
                 !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
20        CONTINUE
           !$OMP END PARALLEL DO   
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
            !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,IX,JY) IF(N>=400)
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
