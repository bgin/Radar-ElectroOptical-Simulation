


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
40        CONTINUE
          !$OMP END PARALLEL DO    
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
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  -- This version written on 25-October-1982.
!*>     Modified on 14-October-1993 to inline the call to DLASSQ.
!*>     Sven Hammarling, Nag Ltd.
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX) !GCC$ ATTRIBUTES inline :: DNRM2 !GCC$ ATTRIBUTES aligned(32) :: DNRM2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
     !DIR$ ATTRIBUTES FORCEINLINE :: DNRM2
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DNRM2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DNRM2
#endif
      implicit none
!*!
!*  -- Reference BLAS level1 routine (version 3.8.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
!*
!*     .. Scalar Arguments ..
      INTEGER INCX,N
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
!*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
!*        The following loop is equivalent to this call to the LAPACK
!*        auxiliary routine:
!*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
!*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF

      DNRM2 = NORM
    
END FUNCTION


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
SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S) !GCC$ ATTRIBUTES inline :: DROT !GCC$ ATTRIBUTES aligned(32) :: DROT
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
     !DIR$ ATTRIBUTES FORCEINLINE :: DROT
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DROT
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DROT
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
      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION DX(*),DY(*)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DX
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DY
!*     ..
1*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
!*     ..
!      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!*
!*       code for both increments equal to 1
         !*
         !$OMP SIMD ALIGNED(DX:64,DY) LINEAR(I:1) UNROLL(10)
         DO I = 1,N
            DTEMP = C*DX(I) + S*DY(I)
            DY(I) = C*DY(I) - S*DX(I)
            DX(I) = DTEMP
         END DO
      ELSE
!*
!*       code for unequal increments or equal increments not equal
!*         to 1
!*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         !$OMP SIMD ALIGNED(DX:64,DY) LINEAR(I:1) UNROLL(10)
         DO I = 1,N
            DTEMP = C*DX(IX) + S*DY(IY)
            DY(IY) = C*DY(IY) - S*DX(IX)
            DX(IX) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
     
END SUBROUTINE


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
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DROTG(DA,DB,C,S) !GCC$ ATTRIBUTES inline :: DROTG !GCC$ ATTRIBUTES aligned(32) :: DROTG
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DROTG(DA,DB,C,S)
     !DIR$ ATTRIBUTES FORCEINLINE :: DROTG
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DROTG
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DROTG
#endif
       implicit none
!*
!*  -- Reference BLAS level1 routine (version 3.8.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
!*
!*     .. Scalar Arguments ..
      DOUBLE PRECISION C,DA,DB,S
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      DOUBLE PRECISION R,ROE,SCALE,Z
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DABS,DSIGN,DSQRT
!*     ..
      ROE = DB
      IF (DABS(DA).GT.DABS(DB)) ROE = DA
      SCALE = DABS(DA) + DABS(DB)
      IF (SCALE.EQ.0.0d0) THEN
         C = 1.0d0
         S = 0.0d0
         R = 0.0d0
         Z = 0.0d0
      ELSE
         R = SCALE*DSQRT((DA/SCALE)**2+ (DB/SCALE)**2)
         R = DSIGN(1.0d0,ROE)*R
         C = DA/R
         S = DB/R
         Z = 1.0d0
         IF (DABS(DA).GT.DABS(DB)) Z = S
         IF (DABS(DB).GE.DABS(DA) .AND. C.NE.0.0d0) Z = 1.0d0/C
      END IF
      DA = R
      DB = Z
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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM) !GCC$ ATTRIBUTES inline :: DROTM !GCC$ ATTRIBUTES aligned(32) :: DROTM
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM)
      !DIR$ ATTRIBUTES FORCEINLINE :: DROTM
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DROTM
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DROTM
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
      !DOUBLE PRECISION DPARAM(5),DX(*),DY(*)
      DOUBLE PRECISION, DIMENSION(5) :: DPARAM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DX
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DY
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,TWO,W,Z,ZERO
      INTEGER I,KX,KY,NSTEPS
!*     ..
!*     .. Data statements ..
      DATA ZERO,TWO/0.D0,2.D0/
!*     ..
!*
      DFLAG = DPARAM(1)
      IF (DFLAG+TWO.EQ.ZERO) RETURN
      IF (INCX.EQ.INCY.AND.INCX.GT.0) THEN
!*
         NSTEPS = N*INCX
         IF (DFLAG.LT.ZERO) THEN
            DH11 = DPARAM(2)
            DH12 = DPARAM(4)
            DH21 = DPARAM(3)
            DH22 = DPARAM(5)
            !$OMP SIMD ALIGNED(DX:64,DY) UNROLL PARTIAL(6)
            DO I = 1,NSTEPS,INCX
               W = DX(I)
               Z = DY(I)
               DX(I) = W*DH11 + Z*DH12
               DY(I) = W*DH21 + Z*DH22
            END DO
         ELSE IF (DFLAG.EQ.ZERO) THEN
            DH12 = DPARAM(4)
            DH21 = DPARAM(3)
            !$OMP SIMD ALIGNED(DX:64,DY) UNROLL PARTIAL(6)
            DO I = 1,NSTEPS,INCX
               W = DX(I)
               Z = DY(I)
               DX(I) = W + Z*DH12
               DY(I) = W*DH21 + Z
            END DO
         ELSE
            DH11 = DPARAM(2)
            DH22 = DPARAM(5)
            !$OMP SIMD ALIGNED(DX:64,DY) UNROLL PARTIAL(6)
            DO I = 1,NSTEPS,INCX
               W = DX(I)
               Z = DY(I)
               DX(I) = W*DH11 + Z
               DY(I) = -W + DH22*Z
            END DO
         END IF
      ELSE
         KX = 1
         KY = 1
         IF (INCX.LT.0) KX = 1 + (1-N)*INCX
         IF (INCY.LT.0) KY = 1 + (1-N)*INCY
!*
         IF (DFLAG.LT.ZERO) THEN
            DH11 = DPARAM(2)
            DH12 = DPARAM(4)
            DH21 = DPARAM(3)
            DH22 = DPARAM(5)
            !$OMP SIMD ALIGNED(DX:64,DY) LINEAR(I:1) UNROLL PARTIAL(6)
            DO I = 1,N
               W = DX(KX)
               Z = DY(KY)
               DX(KX) = W*DH11 + Z*DH12
               DY(KY) = W*DH21 + Z*DH22
               KX = KX + INCX
               KY = KY + INCY
            END DO
         ELSE IF (DFLAG.EQ.ZERO) THEN
            DH12 = DPARAM(4)
            DH21 = DPARAM(3)
            !$OMP SIMD ALIGNED(DX:64,DY) LINEAR(I:1) UNROLL PARTIAL(6)
            DO I = 1,N
               W = DX(KX)
               Z = DY(KY)
               DX(KX) = W + Z*DH12
               DY(KY) = W*DH21 + Z
               KX = KX + INCX
               KY = KY + INCY
            END DO
         ELSE
             DH11 = DPARAM(2)
             DH22 = DPARAM(5)
             !$OMP SIMD ALIGNED(DX:64,DY) LINEAR(I:1) UNROLL PARTIAL(6)
             DO I = 1,N
                W = DX(KX)
                Z = DY(KY)
                DX(KX) = W*DH11 + Z
                DY(KY) = -W + DH22*Z
                KX = KX + INCX
                KY = KY + INCY
            END DO
         END IF
      END IF
   
END SUBROUTINE


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
SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM) !GCC$ ATTRIBUTES inline :: DROTMG !GCC$ aligned(32) :: DROTMG
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
      !DIR$ ATTRIBUTES FORCEINLINE :: DROTMG
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DROTMG
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: DROTMG
#endif
    
!*
!*  -- Reference BLAS level1 routine (version 3.8.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     November 2017
!*
!*     .. Scalar Arguments ..
      DOUBLE PRECISION DD1,DD2,DX1,DY1
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION DPARAM(5)
!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,DTEMP, &
                       DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DABS
!*     ..
!*     .. Data statements ..
!*
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
      DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
!*     ..

      IF (DD1.LT.ZERO) THEN
!*        GO ZERO-H-D-AND-DX1..
         DFLAG = -ONE
         DH11 = ZERO
         DH12 = ZERO
         DH21 = ZERO
         DH22 = ZERO
!*
         DD1 = ZERO
         DD2 = ZERO
         DX1 = ZERO
      ELSE
!*        CASE-DD1-NONNEGATIVE
         DP2 = DD2*DY1
         IF (DP2.EQ.ZERO) THEN
            DFLAG = -TWO
            DPARAM(1) = DFLAG
            RETURN
         END IF
!*        REGULAR-CASE..
         DP1 = DD1*DX1
         DQ2 = DP2*DY1
         DQ1 = DP1*DX1
!*
         IF (DABS(DQ1).GT.DABS(DQ2)) THEN
            DH21 = -DY1/DX1
            DH12 = DP2/DP1
!*
            DU = ONE - DH12*DH21
!*
           IF (DU.GT.ZERO) THEN
             DFLAG = ZERO
             DD1 = DD1/DU
             DD2 = DD2/DU
             DX1 = DX1*DU
           END IF
         ELSE

            IF (DQ2.LT.ZERO) THEN
!*              GO ZERO-H-D-AND-DX1..
               DFLAG = -ONE
               DH11 = ZERO
               DH12 = ZERO
               DH21 = ZERO
               DH22 = ZERO
!*
               DD1 = ZERO
               DD2 = ZERO
               DX1 = ZERO
            ELSE
               DFLAG = ONE
               DH11 = DP1/DP2
               DH22 = DX1/DY1
               DU = ONE + DH11*DH22
               DTEMP = DD2/DU
               DD2 = DD1/DU
               DD1 = DTEMP
               DX1 = DY1*DU
            END IF
         END IF

!*     PROCEDURE..SCALE-CHECK
         IF (DD1.NE.ZERO) THEN
            DO WHILE ((DD1.LE.RGAMSQ) .OR. (DD1.GE.GAMSQ))
               IF (DFLAG.EQ.ZERO) THEN
                  DH11 = ONE
                  DH22 = ONE
                  DFLAG = -ONE
               ELSE
                  DH21 = -ONE
                  DH12 = ONE
                  DFLAG = -ONE
               END IF
               IF (DD1.LE.RGAMSQ) THEN
                  DD1 = DD1*GAM**2
                  DX1 = DX1/GAM
                  DH11 = DH11/GAM
                  DH12 = DH12/GAM
               ELSE
                  DD1 = DD1/GAM**2
                  DX1 = DX1*GAM
                  DH11 = DH11*GAM
                  DH12 = DH12*GAM
               END IF
            ENDDO
         END IF

         IF (DD2.NE.ZERO) THEN
            DO WHILE ( (DABS(DD2).LE.RGAMSQ) .OR. (DABS(DD2).GE.GAMSQ) )
               IF (DFLAG.EQ.ZERO) THEN
                  DH11 = ONE
                  DH22 = ONE
                  DFLAG = -ONE
               ELSE
                  DH21 = -ONE
                  DH12 = ONE
                  DFLAG = -ONE
               END IF
               IF (DABS(DD2).LE.RGAMSQ) THEN
                  DD2 = DD2*GAM**2
                  DH21 = DH21/GAM
                  DH22 = DH22/GAM
               ELSE
                  DD2 = DD2/GAM**2
                  DH21 = DH21*GAM
                  DH22 = DH22*GAM
               END IF
            END DO
         END IF

      END IF

      IF (DFLAG.LT.ZERO) THEN
         DPARAM(2) = DH11
         DPARAM(3) = DH21
         DPARAM(4) = DH12
         DPARAM(5) = DH22
      ELSE IF (DFLAG.EQ.ZERO) THEN
         DPARAM(3) = DH21
         DPARAM(4) = DH12
      ELSE
         DPARAM(2) = DH11
         DPARAM(5) = DH22
      END IF

      DPARAM(1) = DFLAG
    
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
SUBROUTINE DSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) !GCC$ ATTRIBUTES hot :: DSBMV !GCC$ ATTRIBUTES aligned(32) :: DSBMV !GCC$ ATTRIBUTES no_stack_protector :: DSBMV
#elif defined(__INTEL_COMPILER) || defined(__ICC)
SUBROUTINE DSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSBMV
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSBMV
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
      INTEGER INCX,INCY,K,LDA,N
      CHARACTER UPLO
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
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KPLUS1,KX,KY,L
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
    !  INFO = 0
    !  IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
    !      INFO = 1
    !  ELSE IF (N.LT.0) THEN
    !      INFO = 2
    !  ELSE IF (K.LT.0) THEN
    !      INFO = 3
   !   ELSE IF (LDA.LT. (K+1)) THEN
    !      INFO = 6
   !   ELSE IF (INCX.EQ.0) THEN
    !      INFO = 8
    !  ELSE IF (INCY.EQ.0) THEN
    !      INFO = 11
    !  END IF
   !   IF (INFO.NE.0) THEN
   !       CALL XERBLA('DSBMV ',INFO)
   !       RETURN
   !   END IF
!!*
!*     Quick return if possible.
!*
!      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!*
!*     Set up the start points in  X  and  Y.
!*
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
!*
!*     Start the operations. In this version the elements of the array A
!*     are accessed sequentially with one pass through A.
!*
!*     First form  y := beta*y.
!*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
             IF (BETA.EQ.ZERO) THEN
                 !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(8)
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
             ELSE
                     !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(6)   
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                    !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(6)
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                        !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(6)
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
!      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
!*
!*        Form  y  when upper triangle of A is stored.
!*
          KPLUS1 = K + 1
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
             !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,L) IF(N>=400)
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  L = KPLUS1 - J
                  DO 50 I = MAX(1,J-K),J - 1
                      Y(I) = Y(I) + TEMP1*A(L+I,J)
                      TEMP2 = TEMP2 + A(L+I,J)*X(I)
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*A(KPLUS1,J) + ALPHA*TEMP2
60            CONTINUE
              !$OMP END PARALLEL DO   
          ELSE
              JX = KX
              JY = KY
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,IX,IY,JX,JY,KX,KY,L) IF(N>=400)
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                  L = KPLUS1 - J
                  DO 70 I = MAX(1,J-K),J - 1
                      Y(IY) = Y(IY) + TEMP1*A(L+I,J)
                      TEMP2 = TEMP2 + A(L+I,J)*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*A(KPLUS1,J) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
                  IF (J.GT.K) THEN
                      KX = KX + INCX
                      KY = KY + INCY
                  END IF
80            CONTINUE
              !$OMP END PARALLEL DO    
          END IF
      ELSE
!*
!*        Form  y  when lower triangle of A is stored.
!*
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,L) IF(N>=400)
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*A(1,J)
                  L = 1 - J
                  DO 90 I = J + 1,MIN(N,J+K)
                      Y(I) = Y(I) + TEMP1*A(L+I,J)
                      TEMP2 = TEMP2 + A(L+I,J)*X(I)
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
100           CONTINUE
              !$OMP END PARALLEL DO    
          ELSE
              JX = KX
              JY = KY
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,IX,IY,JX,JY,L) IF(N>=400)
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*A(1,J)
                  L = 1 - J
                  IX = JX
                  IY = JY
                  DO 110 I = J + 1,MIN(N,J+K)
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*A(L+I,J)
                      TEMP2 = TEMP2 + A(L+I,J)*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
120            CONTINUE
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
SUBROUTINE DSCAL(N,DA,DX,INCX) !GCC$ ATTRIBUTES inline :: DSCAL !GCC$ ATTRIBUTES aligned(32) :: DSCAL
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DSCAL(N,DA,DX,INCX)
    !DIR$ ATTRIBUTES FORCEINLINE :: DSCAL
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSCAL
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSCAL
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
      INTEGER I,M,MP1,NINCX
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MOD
!*     ..
!      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
!*
!*        code for increment equal to 1
!*
!*
!*        clean-up loop
!*
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         !$OMP SIMD ALIGNED(DX:64) LINEAR(I:5)
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      ELSE
!*
!*        code for increment not equal to 1
!*
         NINCX = N*INCX
         !$OMP SIMD ALIGNED(DX:64) UNROLL PARTIAL(6)
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
      END IF
    
END SUBROUTINE


!*>  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!*>  Krogh, Basic linear algebra subprograms for Fortran
!*>  usage, Algorithm No. 539, Transactions on Mathematical
!*>  Software 5, 3 (September 1979), pp. 308-323.
!*>
!*>  REVISION HISTORY  (YYMMDD)
!*>
!*>  791001  DATE WRITTEN
!*>  890831  Modified array declarations.  (WRB)
!*>  890831  REVISION DATE from Version 3.2
!*>  891214  Prologue converted to Version 4.0 format.  (BAB)
!*>  920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!*>  920501  Reformatted the REFERENCES section.  (WRB)
!!*>  070118  Reformat to LAPACK style (JL)
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
DOUBLE PRECISION FUNCTION DSDOT(N,SX,INCX,SY,INCY) !GCC$ ATTRIBUTES inline :: DSDOT !GCC$ ATTRIBUTES aligned(32) :: DSDOT
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  DOUBLE PRECISION FUNCTION DSDOT(N,SX,INCX,SY,INCY)
       !DIR$ ATTRIBUTES FORCEINLINE :: DSDOT
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSDOT
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSDOT
#endif
    use omp_lib
      implicit none
!*
!*  -- Reference BLAS level1 routine (version 3.7.0) --
!*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!*     ..
!*     .. Array Arguments ..
      !REAL SX(*),SY(*)
      REAL, DIMENSION(:), ALLOCATABLE :: SX
      REAL, DIMENSION(:), ALLOCATABLE :: SY
      
!*     ..
!*
!*  Authors:
!*  ========
!*  Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
!*  Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      INTEGER I,KX,KY,NS
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DBLE
!*     ..
      DSDOT = 0.0D0
    !  IF (N.LE.0) RETURN
      IF (INCX.EQ.INCY .AND. INCX.GT.0) THEN
!*
!*     Code for equal, positive, non-unit increments.
!*
         NS = N*INCX
         !$OMP SIMD ALIGNED(SX:64,SY) REDUCTION(+:DSDOT) UNROLL PARTIAL(10)
         DO I = 1,NS,INCX
            DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
         END DO
      ELSE
!*
!*     Code for unequal or nonpositive increments.
!*
         KX = 1
         KY = 1
         IF (INCX.LT.0) KX = 1 + (1-N)*INCX
         IF (INCY.LT.0) KY = 1 + (1-N)*INCY
          !$OMP SIMD ALIGNED(SX:64,SY) REDUCTION(+:DSDOT) UNROLL PARTIAL(10)
         DO I = 1,N
            DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
            KX = KX + INCX
            KY = KY + INCY
         END DO
      END IF
     
END SUBROUTINE

    
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
SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY) !GCC$ ATTRIBUTES hot :: DSPMV !GCC$ ATTRIBUTES aligned(32) :: DSPMV
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSPMV
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSPMV
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
      INTEGER INCX,INCY,N
      CHARACTER UPLO
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION AP(*),X(*),Y(*)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AP
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
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY
!*     ..
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*
!*     Test the input parameters.
!*
!      INFO = 0
!      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
!          INFO = 1
!      ELSE IF (N.LT.0) THEN
!          INFO = 2
!      ELSE IF (INCX.EQ.0) THEN
!          INFO = 6
!      ELSE IF (INCY.EQ.0) THEN
!          INFO = 9
 !     END IF
 !     IF (INFO.NE.0) THEN
 !         CALL XERBLA('DSPMV ',INFO)
!          RETURN
!      END IF
!*
!*     Quick return if possible.
!*
!      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!*
!*     Set up the start points in  X  and  Y.
!*
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
!*
!*     Start the operations. In this version the elements of the array AP
!*     are accessed sequentially with one pass through AP.
!*
!*     First form  y := beta*y.
!*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
             IF (BETA.EQ.ZERO) THEN
                 !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(8)
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
             ELSE
                     !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(6)   
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                   !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(8)
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                     !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(6)   
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      !IF (ALPHA.EQ.ZERO) RETURN
      KK = 1
      IF (LSAME(UPLO,'U')) THEN
!*
!*        Form  y  when AP contains the upper triangle.
!*
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,K,KK) IF(N>=400)
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  K = KK
                  !$OMP SIMD ALIGNED(Y:64,X,AP) LINEAR(I:1) REDUCTION(+:TEMP2) UNROLL PARTIAL(10)
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(I)
                      K = K + 1
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2
                  KK = KK + J
60            CONTINUE
              !$OMP END PARALLEL DO    
          ELSE
              JX = KX
              JY = KY
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,IX,IY,JX,JY,KK) IF(N>=400)
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                   !$OMP SIMD ALIGNED(Y:64,X,AP) LINEAR(I:1) REDUCTION(+:TEMP2) UNROLL PARTIAL(10)
                  DO 70 K = KK,KK + J - 2
                      Y(IY) = Y(IY) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + J
80            CONTINUE
                !$OMP END PARALLEL DO  
          END IF
      ELSE
!*
!*        Form  y  when AP contains the lower triangle.
!*
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,K,KK) IF(N>=400)
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*AP(KK)
                  K = KK + 1
                   !$OMP SIMD ALIGNED(Y:64,X,AP) LINEAR(I:1)  REDUCTION(+:TEMP2) UNROLL PARTIAL(10)
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(I)
                      K = K + 1
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
                  KK = KK + (N-J+1)
100           CONTINUE
                !$OMP END PARALLEL DO  
          ELSE
              JX = KX
              JY = KY
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,IX,IY,JX,JY,KK) IF(N>=400)
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*AP(KK)
                  IX = JX
                  IY = JY
                   !$OMP SIMD ALIGNED(Y:64,X,AP) LINEAR(I:1)  REDUCTION(+:TEMP2) UNROLL PARTIAL(10)
                  DO 110 K = KK + 1,KK + N - J
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + (N-J+1)
120           CONTINUE
               !$OMP END PARALLEL DO   
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
SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP) !GCC$ ATTRIBUTES hot :: DSPR !GCC$ ATTRIBUTES aligned(32) :: DSPR !GCC$ ATTRIBUTES no_stack_protector :: DSPR
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP)
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSPR
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSPR
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
      INTEGER INCX,N
      CHARACTER UPLO
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION AP(*),X(*)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AP
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X
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
      INTEGER I,INFO,IX,J,JX,K,KK,KX
!*     ..
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*
!*     Test the input parameters.
!*
!      INFO = 0
!      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
 !         INFO = 1
 !     ELSE IF (N.LT.0) THEN
 !         INFO = 2
 !     ELSE IF (INCX.EQ.0) THEN
   !       INFO = 5
   !   END IF
      IF (INFO.NE.0) THEN
  !        CALL XERBLA('DSPR  ',INFO)
  !        RETURN
   !   END IF
!*
!*     Quick return if possible.
!*
 !     IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!*
!*     Set the start point in X if the increment is not unity.
!*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
!*
!*     Start the operations. In this version the elements of the array AP
!*     are accessed sequentially with one pass through AP.
!*
      KK = 1
      IF (LSAME(UPLO,'U')) THEN
!*
!*        Form  A  when upper triangle is stored in AP.
!*
         IF (INCX.EQ.1) THEN
            !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,K,KK) IF(N>=400)
              DO 20 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      K = KK
                      !$OMP SIMD ALIGNED(AP:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 10 I = 1,J
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   10                 CONTINUE
                  END IF
                  KK = KK + J
20            CONTINUE
              !$OMP END PARALLEL DO    
          ELSE
             JX = KX
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,IX,JX,KK) IF(N>=400)
              DO 40 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = KX
                       !$OMP SIMD ALIGNED(AP:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 30 K = KK,KK + J - 1
                          AP(K) = AP(K) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
                  KK = KK + J
40             CONTINUE
                !$OMP END PARALLEL DO  
          END IF
      ELSE
!*
!*        Form  A  when lower triangle is stored in AP.
!*
         IF (INCX.EQ.1) THEN
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,K,KK) IF(N>=400)
              DO 60 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      K = KK
                       !$OMP SIMD ALIGNED(AP:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 50 I = J,N
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   50                 CONTINUE
                  END IF
                  KK = KK + N - J + 1
60            CONTINUE
                 !$OMP END PARALLEL DO 
          ELSE
             JX = KX
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,IX,JX,KK) IF(N>=400)
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = JX
                      !$OMP SIMD ALIGNED(AP:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 70 K = KK,KK + N - J
                          AP(K) = AP(K) + X(IX)*TEMP
                          IX = IX + INCX
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
                  KK = KK + N - J + 1
   80         CONTINUE
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
SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP) !GCC$ ATTRIBUTES hot :: DSPR2 !GCC$ ATTRIBUTES aligned(32) :: DSPR2 !GCC$ ATTRIBUTES no_stack_protector :; DSPR2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSPR2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSPR2
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
      INTEGER INCX,INCY,N
      CHARACTER UPLO
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION AP(*),X(*),Y(*)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AP
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Y
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY
*!     ..
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*
!*     Test the input parameters.
!*
    !  INFO = 0
    !  IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
    !      INFO = 1
     ! ELSE IF (N.LT.0) THEN
     !     INFO = 2
     ! ELSE IF (INCX.EQ.0) THEN
     !     INFO = 5
     ! ELSE IF (INCY.EQ.0) THEN
     !     INFO = 7
     ! END IF
     ! IF (INFO.NE.0) THEN
     !     CALL XERBLA('DSPR2 ',INFO)
     !     RETURN
     ! END IF
!*
!*     Quick return if possible.
!*
!      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!*
!*     Set up the start points in X and Y if the increments are not both
!*     unity.
!*
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
          IF (INCY.GT.0) THEN
              KY = 1
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
          JX = KX
          JY = KY
      END IF
!*
!*     Start the operations. In this version the elements of the array AP
!*     are accessed sequentially with one pass through AP.
!*
      KK = 1
      IF (LSAME(UPLO,'U')) THEN
!*
!*        Form  A  when upper triangle is stored in AP.
!*
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
             !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,K,KK) IF(N>=400)
              DO 20 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      K = KK
                      !$OMP SIMD ALIGNED(AP:64,X,Y) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 10 I = 1,J
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
   10                 CONTINUE
                  END IF
                  KK = KK + J
20            CONTINUE
              !$OMP END PARALLEL DO    
          ELSE
             DO 40 J = 1,N
                   !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,IX,IY,JX,JY,KK) IF(N>=400)
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = KX
                      IY = KY
                       !$OMP SIMD ALIGNED(AP:64,X,Y) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 30 K = KK,KK + J - 1
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + J
40          CONTINUE
            !$OMP END PARALLEL DO      
          END IF
      ELSE
!*
!*        Form  A  when lower triangle is stored in AP.
!*
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,K,KK) IF(N>=400)
              DO 60 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      K = KK
                        !$OMP SIMD ALIGNED(AP:64,X,Y) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 50 I = J,N
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
   50                 CONTINUE
                  END IF
                  KK = KK + N - J + 1
60            CONTINUE
              !$OMP END PARALLEL DO  
          ELSE
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,IX,IY,JX,JY,KK) IF(N>=400)   
              DO 80 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = JX
                      IY = JY
                        !$OMP SIMD ALIGNED(AP:64,X,Y) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 70 K = KK,KK + N - J
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + N - J + 1
80            CONTINUE
               !$OMP END PARALLEL DO   
          END IF
      END IF

END SUBROUTINE

    
!¬*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!1*> \date November 2017
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
SUBROUTINE DSWAP(N,DX,INCX,DY,INCY) !GCC$ ATTRIBUTES inline :: DSWAP !GCC$ ATTRIBUTES aligned(32) :: DSWAP
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
    !DIR$ ATTRIBUTES FORCEINLINE :: DSWAP
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSWAP
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSWAP
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
      !      DOUBLE PRECISION DX(*),DY(*)
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
!      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!*
!*       code for both increments equal to 1
!*
!*
!*       clean-up loop
!*
         M = MOD(N,3)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DX(I)
               DX(I) = DY(I)
               DY(I) = DTEMP
            END DO
            IF (N.LT.3) RETURN
         END IF
         MP1 = M + 1
         !$OMP SIMD ALIGNED(DX:64,DY) LINEAR(I:3)
         DO I = MP1,N,3
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
            DTEMP = DX(I+1)
            DX(I+1) = DY(I+1)
            DY(I+1) = DTEMP
            DTEMP = DX(I+2)
            DX(I+2) = DY(I+2)
            DY(I+2) = DTEMP
         END DO
      ELSE
!*
!*       code for unequal increments or equal increments not equal
!*         to 1
1*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
           !$OMP SIMD ALIGNED(DX:64,DY) LINEAR(I:1)
         DO I = 1,N
            DTEMP = DX(IX)
            DX(IX) = DY(IY)
            DY(IY) = DTEMP
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
SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC) !GCC$ ATTRIBUTES hot :: DSYMM !GCC$ ATTRIBUTES aligned(32) :: DSYMM !GCC$ ATTRIBUTES no_stack_protector :: DSYMM
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSYMM
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSYMM
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
      INTEGER LDA,LDB,LDC,M,N
      CHARACTER SIDE,UPLO
!*     ..
!*     .. Array Arguments ..
      ! DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
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
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MAX
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,J,K,NROWA
      LOGICAL UPPER
!*     ..
!*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!*     ..
!*
!*     Set NROWA as the number of rows of A.
!*
      IF (LSAME(SIDE,'L')) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      UPPER = LSAME(UPLO,'U')
!*
!*     Test the input parameters.
!*
 !     INFO = 0
      IF ((.NOT.LSAME(SIDE,'L')) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
 !         INFO = 1
  !    ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
 !         INFO = 2
  !    ELSE IF (M.LT.0) THEN
  !!        INFO = 3
 ! 1    ELSE IF (N.LT.0) THEN
  !        INFO = 4
    !  ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
    !      INFO = 7
    !  ELSE IF (LDB.LT.MAX(1,M)) THEN
    !      INFO = 9
    !  ELSE IF (LDC.LT.MAX(1,M)) THEN
    !      INFO = 12
    !  END IF
   !   IF (INFO.NE.0) THEN
    !      CALL XERBLA('DSYMM ',INFO)
   !       RETURN
    !  END IF
!*
!*     Quick return if possible.
!*
   !   IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
   !  +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!*
!*     And when  alpha.eq.zero.
!*
      IF (ALPHA.EQ.ZERO) THEN
         IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
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
      IF (LSAME(SIDE,'L')) THEN
!*
!*        Form  C := alpha*A*B + beta*C.
!*
         IF (UPPER) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) COLLAPSE(2) PRIVATE(J,I,TEMP1,TEMP2) IF(N>=400)
              DO 70 J = 1,N
                  DO 60 I = 1,M
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                      !$OMP SIMD ALIGNED(C:64,B,A) REDUCTION(+:TEMP2) UNROLL PARTIAL(10)
                      DO 50 K = 1,I - 1
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*A(K,I)
   50                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + &
                                   ALPHA*TEMP2
                      END IF
   60             CONTINUE
70             CONTINUE
                !$OMP END PARALLEL DO      
          ELSE
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J) IF(N>=400)        
              DO 100 J = 1,N
                  DO 90 I = M,1,-1
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                       !$OMP SIMD ALIGNED(C:64,B,A) REDUCTION(+:TEMP2) UNROLL PARTIAL(10)
                      DO 80 K = I + 1,M
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*A(K,I)
   80                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + &
                                   ALPHA*TEMP2
                      END IF
   90             CONTINUE
100            CONTINUE
                  !$OMP END PARALLEL DO    
          END IF
      ELSE
!*
!*        Form  C := alpha*B*A + beta*C.
         !*
          !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEAFULT(SHARED) PRIVATE(J,TEMP1) IF(N>=400)
          DO 170 J = 1,N
              TEMP1 = ALPHA*A(J,J)
              IF (BETA.EQ.ZERO) THEN
                   !$OMP SIMD ALIGNED(C:64,B,A) LINEAR(I:1) UNROLL PARTIAL(6)
                  DO 110 I = 1,M
                      C(I,J) = TEMP1*B(I,J)
  110             CONTINUE
              ELSE
                     !$OMP SIMD ALIGNED(C:64,B,A) LINEAR(I:1) UNROLL PARTIAL(6)  
                  DO 120 I = 1,M
                      C(I,J) = BETA*C(I,J) + TEMP1*B(I,J)
  120             CONTINUE
              END IF
              DO 140 K = 1,J - 1
                  IF (UPPER) THEN
                      TEMP1 = ALPHA*A(K,J)
                  ELSE
                      TEMP1 = ALPHA*A(J,K)
                  END IF
                    !$OMP SIMD ALIGNED(C:64,B,A) LINEAR(I:1) UNROLL PARTIAL(10)
                  DO 130 I = 1,M
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
  130             CONTINUE
  140         CONTINUE
              DO 160 K = J + 1,N
                  IF (UPPER) THEN
                      TEMP1 = ALPHA*A(J,K)
                  ELSE
                      TEMP1 = ALPHA*A(K,J)
                  END IF
                    !$OMP SIMD ALIGNED(C:64,B,A) LINEAR(I:1) UNROLL PARTIAL(10)
                  DO 150 I = 1,M
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
  150             CONTINUE
  160         CONTINUE
170        CONTINUE
          !$OMP END PARALLEL DO            
      END IF

END SUBROUTINE


!1*> \author Univ. of Tennessee
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
SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) !GCC$ ATTRIBUTES inline :: DSYMV !GCC$ ATTRIBUTES aligned(32) :: DSYMV !GCC$ ATTRIBUTES no_stack_protector :: DSYMV
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSYMV
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSYMV
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
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
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
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
!*     ..
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MAX
!*     ..
!*
!*     Test the input parameters.
!*
!      INFO = 0
!      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
   !       INFO = 1
   !   ELSE IF (N.LT.0) THEN
   !       INFO = 2
  !    ELSE IF (LDA.LT.MAX(1,N)) THEN
   !       INFO = 5
   !   ELSE IF (INCX.EQ.0) THEN
   !       INFO = 7
   !   ELSE IF (INCY.EQ.0) THEN
  !        INFO = 10
  !    END IF
  !    IF (INFO.NE.0) THEN
  !        CALL XERBLA('DSYMV ',INFO)
  !        RETURN
  !    END IF
!*
!*     Quick return if possible.
!*
!      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!*
!*     Set up the start points in  X  and  Y.
!*
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
!*
!!*     Start the operations. In this version the elements of A are
!*     accessed sequentially with one pass through the triangular part
!*     of A.
!*
!*     First form  y := beta*y.
!*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
             IF (BETA.EQ.ZERO) THEN
                !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(8)
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
             ELSE
                    !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(6)   
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(8)
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                    !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) UNROLL PARTIAL(6)   
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
     ! IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
!*
!*        Form  y  when A is stored in upper triangle.
!*
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
             !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2) IF(N>=400)
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                    !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) REDUCTION(+:TEMP2) UNROLL PARTIAL(10)   
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
60             CONTINUE
                !$OMP END PARALLEL DO  
          ELSE
              JX = KX
              JY = KY
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,IX,IY,JX,JY) IF(N>=400)
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                   !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) REDUCTION(+:TEMP2) UNROLL PARTIAL(10)   
                  DO 70 I = 1,J - 1
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
80             CONTINUE
                !$OMP END PARALLEL DO  
          END IF
      ELSE
!*
!*        Form  y  when A is stored in lower triangle.
!*
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2) IF(N>=400)
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*A(J,J)
                   !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) REDUCTION(+:TEMP2) UNROLL PARTIAL(10)   
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
100            CONTINUE
               !$OMP END PARALLEL DO   
          ELSE
              JX = KX
              JY = KY
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,IX,IY,JX,JY) IF(N>=400)
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*A(J,J)
                  IX = JX
                  IY = JY
                    !$OMP SIMD ALIGNED(Y:64) LINEAR(I:1) REDUCTION(+:TEMP2) UNROLL PARTIAL(10)   
                  DO 110 I = J + 1,N
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
120           CONTINUE
               !$OMP END PARALLEL DO   
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
SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA) !GCC$ ATTRIBUTES inline :: DSYR !GCC$ ATTRIBUTES aligned(32) :: DSYR
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
    !DIR$ ATTRIBUTES FORCEINLINE :: DSYR
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSYR
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSYR
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
      INTEGER INCX,LDA,N
      CHARACTER UPLO
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION A(LDA,*),X(*)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: X
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
      INTEGER I,INFO,IX,J,JX,KX
!*     ..
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MAX
!*     ..
!*
!*     Test the input parameters.
!*
     ! INFO = 0
     ! IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
     !     INFO = 1
     ! ELSE IF (N.LT.0) THEN
     !     INFO = 2
     ! ELSE IF (INCX.EQ.0) THEN
     !     INFO = 5
     ! ELSE IF (LDA.LT.MAX(1,N)) THEN
     !     INFO = 7
     ! END IF
     ! IF (INFO.NE.0) THEN
     !     CALL XERBLA('DSYR  ',INFO)
     !     RETURN
    !  END IF
!*
!*     Quick return if possible.
!*
!      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!*
!*     Set the start point in X if the increment is not unity.
!*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
!*
!*     Start the operations. In this version the elements of A are
!*     accessed sequentially with one pass through the triangular part
!*     of A.
!*
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when A is stored in upper triangle.
*
          IF (INCX.EQ.1) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP) IF(N>=400)
              DO 20 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                     TEMP = ALPHA*X(J)
                      !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 10 I = 1,J
                          A(I,J) = A(I,J) + X(I)*TEMP
   10                 CONTINUE
                  END IF
20            CONTINUE
             !$OMP END PARALLEL DO     
          ELSE
             JX = KX
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,IX,JX) IF(N>=400)
              DO 40 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = KX
                        !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 30 I = 1,J
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
40           CONTINUE
              !$OMP END PARALLEL DO    
          END IF
      ELSE
!*
!*        Form  A  when A is stored in lower triangle.
!*
         IF (INCX.EQ.1) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP) IF(N>=400)
              DO 60 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                     TEMP = ALPHA*X(J)
                       !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 50 I = J,N
                          A(I,J) = A(I,J) + X(I)*TEMP
   50                 CONTINUE
                  END IF
60            CONTINUE
              !$OMP END PARALLEL DO
          ELSE
             JX = KX
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP,IX,JX) IF(N>=400)
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = JX
                        !$OMP SIMD ALIGNED(A:64,X) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 70 I = J,N
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
80           CONTINUE
                !$OMP END PARALLEL DO  
          END IF
      END IF

END SUBROUTINE


!*> \author Univ. of Tennessee
!!*> \author Univ. of California Berkeley
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
SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA) !GCC$ ATTTRIBUTES inline :: DSYR2 !GCC$ ATTRIBUTES aligned(32) :: DSYR2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
     !DIR$ ATTRIBUTES FORCEINLINE :: DSYR2
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSYR2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSYR2
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
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
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
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
!*     ..
!*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!*     ..
!*     .. External Subroutines ..
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MAX
!*     ..
!*
!*     Test the input parameters.
!*
    !  INFO = 0
    !  IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
     !     INFO = 1
    !  ELSE IF (N.LT.0) THEN
    !      INFO = 2
    !  ELSE IF (INCX.EQ.0) THEN
    !      INFO = 5
    !!  ELSE IF (INCY.EQ.0) THEN
     !     INFO = 7
    ! ! ELSE IF (LDA.LT.MAX(1,N)) THEN
    !      INFO = 9
    !  END IF
    !  IF (INFO.NE.0) THEN
   !       CALL XERBLA('DSYR2 ',INFO)
   !       RETURN
   !   END IF
!*
!*     Quick return if possible.
!*
   !   IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!*
!*     Set up the start points in X and Y if the increments are not both
!*     unity.
!*
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
          IF (INCY.GT.0) THEN
              KY = 1
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
          JX = KX
          JY = KY
      END IF
!*
!*     Start the operations. In this version the elements of A are
!*     accessed sequentially with one pass through the triangular part
!*     of A.
!*
      IF (LSAME(UPLO,'U')) THEN
!*
!*        Form  A  when A is stored in the upper triangle.
!*
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2) IF(N>=400)
              DO 20 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      !$OMP SIMD ALIGNED(A:64,X,Y) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 10 I = 1,J
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   10                 CONTINUE
                  END IF
20            CONTINUE
                 !$OMP END PARALLEL DO 
         ELSE
                   !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,IX,IY,JX,JY) IF(N>=400)
              DO 40 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = KX
                      IY = KY
                       !$OMP SIMD ALIGNED(A:64,X,Y) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 30 I = 1,J
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
40            CONTINUE
                !$OMP END PARALLEL DO
          END IF
      ELSE
!*
!*        Form  A  when A is stored in the lower triangle.
!*
         IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
               !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2) IF(N>=400)
              DO 60 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                       !$OMP SIMD ALIGNED(A:64,X,Y) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 50 I = J,N
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   50                 CONTINUE
                  END IF
60            CONTINUE
               !$OMP END PARALLEL DO   
        ELSE
                   !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J,TEMP1,TEMP2,IX,IY,JX,JY) IF(N>=400)
              DO 80 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = JX
                      IY = JY
                       !$OMP SIMD ALIGNED(A:64,X,Y) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 70 I = J,N
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
80            CONTINUE
              !$OMP END PARALLEL DO    
          END IF
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
!*> \ingroup double_blas_level3
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 3 Blas routine.
!*>
!*>
!*>  -- Written on 8-February-1989.
!*>     Jack Dongarra, Argonne National Laboratory.
!*>     Iain Duff, AERE Harwell.
!!*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!*>     Sven Hammarling, Numerical Algorithms Group Ltd.
!*> \endverbatim
!*>
!*  =====================================================================
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) !GCC$ ATTRIBUTES hot :: DSYR2K !GCC$ ATTRIBUTES aligned(32) :: DSYR2K !GCC$ ATTRIBUTES no_stack_protector :: DSYR2K
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSYR2
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSYR2
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
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
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
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MAX
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
!*     ..
!*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!*     ..
!*
!*     Test the input parameters.
!*
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
!*
!      INFO = 0
!      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
!          INFO = 1
!      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
!              (.NOT.LSAME(TRANS,'T')) .AND.
!              (.NOT.LSAME(TRANS,'C'))) THEN
!          INFO = 2
!      ELSE IF (N.LT.0) THEN
!          INFO = 3
 !     ELSE IF (K.LT.0) THEN
!          INFO = 4
!      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
 !         INFO = 7
!!      ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
 !         INFO = 9
 !     ELSE IF (LDC.LT.MAX(1,N)) THEN
 !         INFO = 12
 !     END IF
 !     IF (INFO.NE.0) THEN
 !         CALL XERBLA('DSYR2K',INFO)
 !         RETURN
 !     END IF
!*
!*     Quick return if possible.
!*
!      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
!     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
!*
!*     And when  alpha.eq.zero.
!*
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
!*
!*     Start the operations.
!*
      IF (LSAME(TRANS,'N')) THEN
!*
!*        Form  C := alpha*A*B**T + alpha*B*A**T + C.
!*
         IF (UPPER) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J) IF(N>=400)
              DO 130 J = 1,N
                 IF (BETA.EQ.ZERO) THEN
                     !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                 ELSE IF (BETA.NE.ONE) THEN
                       !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)   
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          !$OMP SIMD ALIGNED(C:64,A,B) LINEAR(I:1) UNROLL PARTIAL(10)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 + &
                                       B(I,L)*TEMP2
  110                     CONTINUE
                      END IF
  120             CONTINUE
130            CONTINUE
               !$OMP END PARALLEL DO        
             ELSE
                    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J) IF(N>=400)    
              DO 180 J = 1,N
                 IF (BETA.EQ.ZERO) THEN
                        !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(10)
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                        !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)  
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                            !$OMP SIMD ALIGNED(C:64,A,B) LINEAR(I:1) UNROLL PARTIAL(10)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 + &
                                       B(I,L)*TEMP2
  160                     CONTINUE
                      END IF
  170             CONTINUE
180            CONTINUE
                !$OMP END PARALLEL DO      
          END IF
      ELSE
!*
!*        Form  C := alpha*A**T*B + alpha*B**T*A + C.
!*
         IF (UPPER) THEN
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) COLLAPSE(2) PRIVATE(J,I,TEMP1,TEMP2) IF(N>=400)    
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                       !$OMP SIMD ALIGNED(C:64,A,B) LINEAR(I:1) REDUCTION(+:TEMP1,TEMP2) UNROLL PARTIAL(10)
                      DO 190 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + &
                                   ALPHA*TEMP2
                      END IF
  200             CONTINUE
210            CONTINUE
                !$OMP END PARALLEL DO      
            ELSE
                  !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED)  PRIVATE(J,TEMP1,TEMP2) IF(N>=400)         
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                        !$OMP SIMD ALIGNED(C:64,A,B) LINEAR(I:1) REDUCTION(+:TEMP1,TEMP2) UNROLL PARTIAL(10)
                      DO 220 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + &
                                   ALPHA*TEMP2
                      END IF
  230             CONTINUE
240            CONTINUE
                !$OMP END PARALLEL DO      
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
SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC) !GCC$ ATTRIBUTES hot :: DSYRK !GCC$ ATTRIBUTES aligned(32) :: DSYRK !GCC$ ATTRIBUTES no_stack_protector :: DSYRK
#elif defined(__INTEL_COMPILER) || defined(__ICC)
  SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DSYRK
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: DSYRK
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
      INTEGER K,LDA,LDC,N
      CHARACTER TRANS,UPLO
!*     ..
!*     .. Array Arguments ..
      !DOUBLE PRECISION A(LDA,*),C(LDC,*)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
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
!      EXTERNAL XERBLA
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MAX
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
!*     ..
!*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!*     ..
!*
!*     Test the input parameters.
!*
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
!*
    !  INFO = 0
    !  IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
    !      INFO = 1
    !  ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
    ! +         (.NOT.LSAME(TRANS,'T')) .AND.
    ! +         (.NOT.LSAME(TRANS,'C'))) THEN
    !      INFO = 2
    !  ELSE IF (N.LT.0) THEN
    !      INFO = 3
    !  ELSE IF (K.LT.0) THEN
    !      INFO = 4
    !  ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
    !      INFO = 7
    !  ELSE IF (LDC.LT.MAX(1,N)) THEN
    !      INFO = 10
    !  END IF
   !   IF (INFO.NE.0) THEN
    !      CALL XERBLA('DSYRK ',INFO)
    !      RETURN
    !  END IF
!*
!*     Quick return if possible.
!*
!      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
!     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
!*
!!*     And when  alpha.eq.zero.
!*
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
!*
!*     Start the operations.
!*
      IF (LSAME(TRANS,'N')) THEN
!*
!*        Form  C := alpha*A*A**T + beta*C.
!*
         IF (UPPER) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J) IF(N>=400)
              DO 130 J = 1,N
                 IF (BETA.EQ.ZERO) THEN
                      !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(8)
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                 ELSE IF (BETA.NE.ONE) THEN
                       !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6) 
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                         TEMP = ALPHA*A(J,L)
                            !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(10)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  110                     CONTINUE
                      END IF
  120             CONTINUE
130           CONTINUE
               !$OMP END PARALLEL DO       
          ELSE
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J) IF(N>=400)        
              DO 180 J = 1,N
                 IF (BETA.EQ.ZERO) THEN
                        !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(8)
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                 ELSE IF (BETA.NE.ONE) THEN
                        !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(6)   
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                         TEMP = ALPHA*A(J,L)
                           !$OMP SIMD ALIGNED(C:64) LINEAR(I:1) UNROLL PARTIAL(10)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  160                     CONTINUE
                      END IF
  170             CONTINUE
180            CONTINUE
                !$OMP END PARALLEL DO      
          END IF
      ELSE
!*
!*        Form  C := alpha*A**T*A + beta*C.
!*
         IF (UPPER) THEN
              !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) COLLAPSE(2) PRIVATE(J,I,TEMP) IF(N>=400)
              DO 210 J = 1,N
                  DO 200 I = 1,J
                     TEMP = ZERO
                      !$OMP SIMD ALIGNED(A:64) LINEAR(L:1) REDUCTION(+:TEMP) UNROLL PARTIAL(10)
                      DO 190 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  200             CONTINUE
210           CONTINUE
               !$OMP END PARALLEL DO       
           ELSE
                !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(J) IF(N>=400)         
              DO 240 J = 1,N
                  DO 230 I = J,N
                     TEMP = ZERO
                       !$OMP SIMD ALIGNED(A:64) LINEAR(L:1) REDUCTION(+:TEMP) UNROLL PARTIAL(10)
                      DO 220 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  230             CONTINUE
240            CONTINUE
                !$OMP END PARALLEL DO       
          END IF
      END IF

END SUBROUTINE
