




      LOGICAL          FUNCTION LSAME( CA, CB)
        !DIR$ ATTRIBUTES FORCEINLINE :: LSAME
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: LSAME
        !DIR$ OPTIMIZE : 3

#if 0
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     March 26, 1990
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  This version of the routine is only correct for ASCII code.
*  Installers must modify the routine for other character-codes.
*
*  For EBCDIC systems the constant IOFF must be changed to -64.
*  For CDC systems using 6-12 bit representations, the system-
*  specific code in comments must be activated.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
*
*     .. Parameters ..
#endif
      INTEGER            IOFF
      PARAMETER        ( IOFF = 32 )
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!*     ..
!*     .. Executable Statements ..
!*
!*     Test if the characters are equal
!*
      LSAME = CA.EQ.CB
!*
!*     Now test for equivalence
!*
      IF( .NOT.LSAME ) THEN
         LSAME = ICHAR( CA ) - IOFF.EQ.ICHAR( CB )
      END IF
      IF( .NOT.LSAME ) THEN
         LSAME = ICHAR( CA ).EQ.ICHAR( CB ) - IOFF
      END IF

END FUNCTION







SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,&
BETA,C,LDC,INFO,MB,NB,NBT,KB)
   
       !dir$ attributes code_align : 32 :: DGEMM
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: DGEMM
!     *     .. Scalar Arguments ..
      use omp_lib
      implicit none
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC, INFO
      INTEGER            MB, NB, NBT, KB 
      DOUBLE PRECISION   ALPHA, BETA
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A,B,C
#if 0
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
!*     op( X ) = X   or   op( X ) = X',
!*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  -- Modified in October-1997.
*     Superscalar GEMM-Based Level 3 BLAS (Version 0.1).
*     Per Ling, Department of Computing Science,
*     Umea University, Sweden.
*
*
*     .. Local Scalars ..
*  -- Modified by Bernard Gingold in 21-02-2021  
#endif
      INTEGER            I, II, ISEC, UISEC, J, JJ, JSEC, UJSEC, &
                         L, LL, LSEC, ULSEC, NROWA, NROWB
      LOGICAL            NOTA, NOTB
      DOUBLE PRECISION   DELTA
      DOUBLE PRECISION   F11, F12, F21, F22, F31, F32, F41, F42
      DOUBLE PRECISION   F13, F14, F23, F24, F33, F34, F43, F44
      DOUBLE PRECISION   TMP0,TMP1,TMP2,TMP3
      DOUBLE PRECISION   C0,C1,C2,C3,C4,C5,C6,C7,C8
                         
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
*     .. External Subroutines ..
      !EXTERNAL           XERBLA
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!*     .. User specified parameters for DGEMM ..
!      Moved to subroutine arguments.
     ! PARAMETER        ( MB = 32, NB = 1024, NBT = 96, KB = 32 )
      DOUBLE PRECISION   T1( KB, MB ), T2( KB, NBT )
!*     ..
!*     .. Executable Statements ..
!*
!*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!*     transposed and set NROWA and NROWB as the number of rows of A  and
!*     the number of rows of B respectively.
!*

      NOTA = LSAME( TRANSA, 'N' )
      NOTB = LSAME( TRANSB, 'N' )
      IF ( NOTA ) THEN
         NROWA = M
      ELSE
         NROWA = K
      END IF
      IF ( NOTB ) THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!*
!*     Test the input parameters.
!*
      INFO = 0
      IF( ( .NOT.NOTA ).AND.( .NOT. LSAME( TRANSA,'C' ) ) &
                               .AND.( .NOT.LSAME( TRANSA, 'T' ) ) )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB ).AND.( .NOT.LSAME( TRANSB, 'C' ) ) &
                               .AND.( .NOT.LSAME( TRANSB, 'T' ) ) )THEN
         INFO = 2
      ELSE IF( M.LT.0 )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( K.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
!CALL XERBLA( 'DGEMM ', INFO )
          RETURN
      END IF
!*
!*     Quick return if possible.
!*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
          ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) RETURN
   
!*
!*     And when alpha.eq.zero.
!*
      IF( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) )THEN
         IF( BETA.EQ.ZERO )THEN
            UISEC = M-MOD( M, 4 )
            DO 30, J = 1, N
   !dir$ assume_aligned C:64
   !dir$ vector aligned 
   !dir$ vector vectorlength(8)
   !dir$ vector always
       
               DO 10, I = 1, UISEC, 4
                  C( I, J ) = ZERO
                  C( I+1, J ) = ZERO
                  C( I+2, J ) = ZERO
                  C( I+3, J ) = ZERO
   10          CONTINUE
   !dir$ assume_aligned C:64
   !dir$ vector aligned 
   !dir$ vector vectorlength(8)
   !dir$ vector always
               DO 20, I = UISEC+1, M
                  C( I, J ) = ZERO
   20          CONTINUE
   30       CONTINUE
         ELSE
            UISEC = M-MOD( M, 4 )
            DO 60, J = 1, N
   !dir$ assume_aligned C:64
   !dir$ assume_aligned BETA:64
   !dir$ vector aligned 
   !dir$ vector vectorlength(8)
   !dir$ vector always
               DO 40, I = 1, UISEC, 4
                  C( I, J ) = BETA*C( I, J )
                  C( I+1, J ) = BETA*C( I+1, J )
                  C( I+2, J ) = BETA*C( I+2, J )
                  C( I+3, J ) = BETA*C( I+3, J )
   40          CONTINUE
   !dir$ assume_aligned C:64
   !dir$ assume_aligned BETA:64
   !dir$ vector aligned 
   !dir$ vector vectorlength(8)
   !dir$ vector always
               DO 50, I = UISEC+1, M
                  C( I, J ) = BETA*C( I, J )
   50          CONTINUE
   60       CONTINUE
         END IF
         RETURN
      END IF
      T1(:,:) = ZERO
      T2(:,:) = ZERO
!*
!*     Start the operations.
!*     
      IF (NOTB) THEN
!*
!*        Form  C := alpha*A*B + beta*C or C := alpha*A'*B + beta*C.
!*
         TMP0 = 0.0D0
         TMP1 = 0.0D0
         TMP2 = 0.0D0
         TMP3 = 0.0D0

         DO 250 JJ = 1, N, NB
            JSEC = MIN( NB, N-JJ+1 )
            UJSEC = JSEC-MOD( JSEC, 4 )
            DO 240 LL = 1, K, KB
               LSEC = MIN( KB, K-LL+1 )
               ULSEC = LSEC-MOD( LSEC, 2 )
!*
!*              Determine if the block of C should be updated with
!*              beta or not.
!*
               DELTA = ONE
               IF( LL.EQ.1 ) DELTA = BETA
!*
               DO 230 II = 1, M, MB
                  ISEC = MIN( MB, M-II+1 )
!*
!*                 T1 := alpha*A' or T1 := alpha*A, copy the transpose
!*                 or the non-transpose of a rectangular block of
!*                 alpha*A to T1.
!*
                  UISEC = ISEC-MOD( ISEC, 2 )
                  IF( NOTA )THEN
                     DO 80, L = LL, LL+ULSEC-1, 2
                        DO 70, I = II, II+UISEC-1, 2
                           T1( L-LL+1, I-II+1 ) = ALPHA*A( I, L )
                           T1( L-LL+2, I-II+1 ) = ALPHA*A( I, L+1 )
                           T1( L-LL+1, I-II+2 ) = ALPHA*A( I+1, L )
                           T1( L-LL+2, I-II+2 ) = ALPHA*A( I+1, L+1 )
   70                   CONTINUE
                        IF( UISEC.LT.ISEC )THEN
                           T1( L-LL+1, ISEC ) = ALPHA*A( II+ISEC-1, L )
                           T1( L-LL+2, ISEC ) = ALPHA*A( II+ISEC-1, L+1)
                                             
                        END IF
   80                CONTINUE
                     IF( ULSEC.LT.LSEC )THEN
                        DO 90, I = II, II+ISEC-1
                           T1( LSEC, I-II+1 ) = ALPHA*A( I, LL+LSEC-1 )
   90                   CONTINUE
                     END IF
                  ELSE
                     DO 110, I = II, II+UISEC-1, 2
   !dir$ assume_aligned T1:64
   !dir$ assume_aligned A:64
   !dir$ vector aligned 
   !dir$ vector vectorlength(8)
   !dir$ vector always
                        DO 100, L = LL, LL+ULSEC-1, 2
                           T1( L-LL+1, I-II+1 )  = ALPHA*A( L, I )
                           T1( L-LL+1, I-II+2 ) = ALPHA*A( L, I+1 )
                           T1( L-LL+2, I-II+1 ) = ALPHA*A( L+1, I )
                           T1( L-LL+2, I-II+2 ) = ALPHA*A( L+1, I+1 )
  100                   CONTINUE
                        IF( ULSEC.LT.LSEC )THEN
                           T1( LSEC, I-II+1 ) = ALPHA*A( LL+LSEC-1, I )
                           T1( LSEC, I-II+2 ) = ALPHA*A( LL+LSEC-1, I+1)
                                             
                        END IF
  110                CONTINUE
                     IF( UISEC.LT.ISEC )THEN
   !dir$ assume_aligned T1:64
   !dir$ assume_aligned A:64
   !dir$ vector aligned 
   !dir$ vector vectorlength(8)
   !dir$ vector always
                        DO 120, L = LL, LL+LSEC-1
                           T1( L-LL+1, ISEC ) = ALPHA*A( L, II+ISEC-1 )
  120                   CONTINUE
                     END IF
                  END IF
!*
!*                 C := T1'*B + beta*C, update a rectangular block
!*                 of C using 4 by 4 unrolling.
!*
                  UISEC = ISEC-MOD( ISEC, 4 )
                  DO 170 J = JJ, JJ+UJSEC-1, 4
                    
                     !$omp simd private(F11,21,F12,F22,F13,F23,F14,F24,F31,F41,F32, &
                     !$omp&             F42,F33,F43,F34,F44) aligned(C:64) linear(J:4)
                     DO 140 I = II, II+UISEC-1, 4
                        F11 = DELTA*C( I,J )
                        F21 = DELTA*C( I+1,J )
                        F12 = DELTA*C( I,J+1 )
                        F22 = DELTA*C( I+1,J+1 )
                        F13 = DELTA*C( I,J+2 )
                        F23 = DELTA*C( I+1,J+2 )
                        F14 = DELTA*C( I,J+3 )
                        F24 = DELTA*C( I+1,J+3 )
                        F31 = DELTA*C( I+2,J )
                        F41 = DELTA*C( I+3,J )
                        F32 = DELTA*C( I+2,J+1 )
                        F42 = DELTA*C( I+3,J+1 )
                        F33 = DELTA*C( I+2,J+2 )
                        F43 = DELTA*C( I+3,J+2 )
                        F34 = DELTA*C( I+2,J+3 )
                        F44 = DELTA*C( I+3,J+3 )
                        !$omp simd private(TMP0,TMP1,TMP2,TMP3) reduction(+:F11,F21,F12,F22,F13, &
                        !$omp&                         F23,F14,F24,F31,F41,F32,F42,F33,F43,F34,F44) &
                        !$omp& aligned(B:64)
                        DO 130 L = LL, LL+LSEC-1
                           TMP0 = B(L,J)
                           TMP1 = B(L,J+1)
                           TMP2 = B(L,J+2)
                           TMP3 = B(L,J+3)
                           F11 = F11 + T1( L-LL+1, I-II+1 )*TMP0
                           F21 = F21 + T1( L-LL+1, I-II+2 )*TMP0
                           F12 = F12 + T1( L-LL+1, I-II+1 )*TMP1
                           F22 = F22 + T1( L-LL+1, I-II+2 )*TMP1
                           F13 = F13 + T1( L-LL+1, I-II+1 )*TMP2
                           F23 = F23 + T1( L-LL+1, I-II+2 )*TMP2
                           F14 = F14 + T1( L-LL+1, I-II+1 )*TMP3
                           F24 = F24 + T1( L-LL+1, I-II+2 )*TMP3
                           F31 = F31 + T1( L-LL+1, I-II+3 )*TMP0
                           F41 = F41 + T1( L-LL+1, I-II+4 )*TMP0
                           F32 = F32 + T1( L-LL+1, I-II+3 )*TMP1
                           F42 = F42 + T1( L-LL+1, I-II+4 )*TMP1
                           F33 = F33 + T1( L-LL+1, I-II+3 )*TMP2
                           F43 = F43 + T1( L-LL+1, I-II+4 )*TMP2
                           F34 = F34 + T1( L-LL+1, I-II+3 )*TMP3
                           F44 = F44 + T1( L-LL+1, I-II+4 )*TMP3
  130                   CONTINUE
                        C( I,J ) = F11
                        C( I+1, J ) = F21
                        C( I, J+1 ) = F12
                        C( I+1, J+1 ) = F22
                        C( I, J+2 ) = F13
                        C( I+1, J+2 ) = F23
                        C( I, J+3 ) = F14
                        C( I+1, J+3 ) = F24
                        C( I+2, J ) = F31
                        C( I+3, J ) = F41
                        C( I+2, J+1 ) = F32
                        C( I+3, J+1 ) = F42
                        C( I+2, J+2 ) = F33
                        C( I+3, J+2 ) = F43
                        C( I+2, J+3 ) = F34
                        C( I+3, J+3 ) = F44
  140                CONTINUE
                     IF( UISEC.LT.ISEC )THEN
                        DO 160 I = II+UISEC, II+ISEC-1
                           F11 = DELTA*C( I, J )
                           F12 = DELTA*C( I, J+1 )
                           F13 = DELTA*C( I, J+2 )
                           F14 = DELTA*C( I, J+3 )
                           !$omp simd private(TMP0,TMP1,TMP2,TMP3) reduction(+:F11,F12,F13,F14) &
                           !$omp& aligned(B:64)
                           DO 150 L = LL, LL+LSEC-1
                              TMP0 = T1( L-LL+1, I-II+1 )*B(L,J)
                              F11 = F11 + TMP0
                              TMP1 =  T1( L-LL+1, I-II+1 )*B( L, J+1 )
                              F12 = F12 + TMP1
                              TMP2 = T1( L-LL+1, I-II+1 )*B( L, J+2 )
                              F13 = F13 + TMP2
                              TMP3 = T1( L-LL+1, I-II+1 )*B( L, J+3 )
                              F14 = F14 + TMP3                             
  150                      CONTINUE
                           C( I, J ) = F11
                           C( I, J+1 ) = F12
                           C( I, J+2 ) = F13
                           C( I, J+3 ) = F14
  160                   CONTINUE
                     END IF
  170             CONTINUE
                  IF( UJSEC.LT.JSEC )THEN
                     TMP0 = 0.0D0
                     TMP1 = 0.0D0
                     TMP2 = 0.0D0
                     TMP3 = 0.0D0
                     DO 220 J = JJ+UJSEC, JJ+JSEC-1
                        DO 190 I = II, II+UISEC-1, 4
                           F11 = DELTA*C( I,J )
                           F21 = DELTA*C( I+1, J )
                           F31 = DELTA*C( I+2, J )
                           F41 = DELTA*C( I+3, J )
                           !$omp simd private(TMP0,TMP1,TMP2,TMP3) reduction(+:F11,F12,F13,F14) &
                           !$omp& aligned(B:64) linear(L:1)
                           DO 180 L = LL, LL+LSEC-1
                              TMP0 = T1( L-LL+1, I-II+1 )*B( L, J )
                              F11  = F11 + TMP0
                              TMP1 = T1( L-LL+1, I-II+2 )*B( L, J )
                              F21  = F21 + TMP1
                              TMP2 = T1( L-LL+1, I-II+3 )*B( L, J )
                              F31  = F31 + TMP2
                              TMP3 = T1( L-LL+1, I-II+4 )*B( L, J )
                              F41  + F41 + TMP3
  180                      CONTINUE
                           C( I,J ) = F11
                           C( I+1, J ) = F21
                           C( I+2, J ) = F31
                           C( I+3, J ) = F41
  190                   CONTINUE
                        DO 210 I = II+UISEC, II+ISEC-1
                           F11 = DELTA*C( I, J )
                           DO 200 L = LL, LL+LSEC-1
                              F11 = F11 + T1( L-LL+1, I-II+1 )*B( L, J )
  200                      CONTINUE
                           C( I, J ) = F11
  210                   CONTINUE
  220                CONTINUE
                  END IF
  230          CONTINUE
  240       CONTINUE
  250    CONTINUE
      ELSE
!*
!*        Form  C := alpha*A*B' + beta*C or C := alpha*A'*B' + beta*C.
!*
         C0  = ZERO
         C1  = ZERO
         C2  = ZERO
         C3  = ZERO
         C4  = ZERO
         C5  = ZERO
         C6  = ZERO
         C7  = ZERO
         C8  = ZERO
         
         DO 470 JJ = 1, N, NBT
            JSEC = MIN( NBT, N-JJ+1 )
            DO 460 LL = 1, K, KB
               LSEC = MIN( KB, K-LL+1 )
!*
!*              Determine if the block of C should be updated with
!*              beta or not.
!*
               DELTA = ONE
               IF( LL.EQ.1 ) DELTA = BETA
!*
!*              T2 := alpha*B', copy the transpose of a rectangular
!*              block of alpha*A to T2.
!*
               ULSEC = LSEC-MOD( LSEC, 2 )
               UJSEC = JSEC-MOD( JSEC, 2 )
               DO 270, L = LL, LL+ULSEC-1, 2
                  DO 260, J = JJ, JJ+UJSEC-1, 2
                     T2( L-LL+1, J-JJ+1 ) = ALPHA*B( J, L )
                     T2( L-LL+2, J-JJ+1 ) = ALPHA*B( J, L+1 )
                     T2( L-LL+1, J-JJ+2 ) = ALPHA*B( J+1, L )
                     T2( L-LL+2, J-JJ+2 ) = ALPHA*B( J+1, L+1 )
  260             CONTINUE
                  IF( UJSEC.LT.JSEC )THEN
                     T2( L-LL+1, JSEC ) = ALPHA*B( JJ+JSEC-1, L )
                     T2( L-LL+2, JSEC ) = ALPHA*B( JJ+JSEC-1, L+1 )
                  END IF
  270          CONTINUE
               IF( ULSEC.LT.LSEC )THEN
                  DO 280, J = JJ, JJ+JSEC-1
                     T2( LSEC, J-JJ+1 ) = ALPHA*B( J, LL+LSEC-1 )
  280             CONTINUE
               END IF

               UJSEC = JSEC-MOD( JSEC, 4 )
               DO 450 II = 1, M, MB
                  ISEC = MIN( MB, M-II+1 )
!*
!*                 T1 := alpha*A' or T1 := alpha*A, copy the transpose
!*                 or the non-transpose of a rectangular block of
!*                 alpha*A to T1.
!*
                  UISEC = ISEC-MOD( ISEC, 2 )
                  IF( NOTA )THEN
                     DO 300, L = LL, LL+ULSEC-1, 2
                        DO 290, I = II, II+UISEC-1, 2
                           T1( L-LL+1, I-II+1 ) = A( I, L )
                           T1( L-LL+2, I-II+1 ) = A( I, L+1 )
                           T1( L-LL+1, I-II+2 ) = A( I+1, L )
                           T1( L-LL+2, I-II+2 ) = A( I+1, L+1 )
  290                   CONTINUE
                        IF( UISEC.LT.ISEC )THEN
                           T1( L-LL+1, ISEC ) = A( II+ISEC-1, L )
                           T1( L-LL+2, ISEC ) = A( II+ISEC-1, L+1 )
                        END IF
  300                CONTINUE
                     IF( ULSEC.LT.LSEC )THEN
                        DO 310, I = II, II+ISEC-1
                           T1( LSEC, I-II+1 ) = A( I, LL+LSEC-1 )
  310                   CONTINUE
                     END IF
                  ELSE
                     DO 330, I = II, II+UISEC-1, 2
                        !$omp simd aligned(A:64) linear(L:2) simdlen(8)
                        DO 320, L = LL, LL+ULSEC-1, 2
                           T1( L-LL+1, I-II+1 ) = A( L, I )
                           T1( L-LL+1, I-II+2 ) = A( L, I+1 )
                           T1( L-LL+2, I-II+1 ) = A( L+1, I )
                           T1( L-LL+2, I-II+2 ) = A( L+1, I+1 )
  320                   CONTINUE
                        IF( ULSEC.LT.LSEC )THEN
                           T1( LSEC, I-II+1 ) = A( LL+LSEC-1, I )
                           T1( LSEC, I-II+2 ) = A( LL+LSEC-1, I+1 )
                        END IF
  330                CONTINUE
                     IF( UISEC.LT.ISEC )THEN
                        !$omp simd aligned(A:64) linear(L:1) simdlen(8) 
                        DO 340, L = LL, LL+LSEC-1
                           T1( L-LL+1, ISEC ) = A( L, II+ISEC-1 )
  340                   CONTINUE
                     END IF
                  END IF
!*
!*                 C := T1'*B + beta*C, update a rectangular block
!*                 of C using 4 by 4 unrolling.
!*
                  UISEC = ISEC-MOD( ISEC, 4 )
                  DO 390 J = JJ, JJ+UJSEC-1, 4
                     DO 360 I = II, II+UISEC-1, 4
                        F11 = DELTA*C( I,J )
                        F21 = DELTA*C( I+1,J )
                        F12 = DELTA*C( I,J+1 )
                        F22 = DELTA*C( I+1,J+1 )
                        F13 = DELTA*C( I,J+2 )
                        F23 = DELTA*C( I+1,J+2 )
                        F14 = DELTA*C( I,J+3 )
                        F24 = DELTA*C( I+1,J+3 )
                        F31 = DELTA*C( I+2,J )
                        F41 = DELTA*C( I+3,J )
                        F32 = DELTA*C( I+2,J+1 )
                        F42 = DELTA*C( I+3,J+1 )
                        F33 = DELTA*C( I+2,J+2 )
                        F43 = DELTA*C( I+3,J+2 )
                        F34 = DELTA*C( I+2,J+3 )
                        F44 = DELTA*C( I+3,J+3 )
                        !$omp simd private(C0,C1,C2,C3,C4,C5,C6,C7) &
                        !$omp& reduction(+:F11,F21,F22,F13,F23,F14,F24,F31,F41,F32,F42,F33,F43,F34,F44) &
                        !$omp& simdlen(8) linear(L:1)
                        DO 350 L = LL, LL+LSEC-1
                           C0  = T1(L-LL+1,I-II+1)
                           C1  = T2(L-LL+1,J-JJ+1)
                           F11 = F11 + C0 * C1
                           C2  = T1(L-LL+1,I-II+2)
                           F21 = F21 + C2 * C1
                           C3 = T2(L-LL+1,J-JJ+2)
                           F12 = F12 + C0 * C3
                           F22 = F22 + C2 * C3
                           C4 = T1(L-LL+1,I-II+3)
                           C5 = T2(L-LL+1,J-JJ+3)
                           F13 = F13 + C0 * C4
                           F23 = F23 + C2 * C4
                           C6 = T1(L-LL+1,I-II+4)
                           C7 = T2(L-LL+1,J-JJ+4)
                           F14 = F14 + C0 * C7
                           F24 = F24 + C2 * C7
                           F31 = F31 + C4 * C1
                           F41 = F41 + C6 * C1
                           F32 = F32 + C4 * C3
                           F42 = F42 + C6 * C3
                           F33 = F33 + C4 * C5
                           F43 = F43 + C6 * C5
                           F34 = F34 + C4 * C7
                           F44 = F44 + C6 * C7
                           
  350                   CONTINUE
                        C( I,J ) = F11
                        C( I+1, J ) = F21
                        C( I, J+1 ) = F12
                        C( I+1, J+1 ) = F22
                        C( I, J+2 ) = F13
                        C( I+1, J+2 ) = F23
                        C( I, J+3 ) = F14
                        C( I+1, J+3 ) = F24
                        C( I+2, J ) = F31
                        C( I+3, J ) = F41
                        C( I+2, J+1 ) = F32
                        C( I+3, J+1 ) = F42
                        C( I+2, J+2 ) = F33
                        C( I+3, J+2 ) = F43
                        C( I+2, J+3 ) = F34
                        C( I+3, J+3 ) = F44
  360                CONTINUE
                     IF( UISEC.LT.ISEC )THEN
                        DO 380 I = II+UISEC, II+ISEC-1
                           F11 = DELTA*C( I, J )
                           F12 = DELTA*C( I, J+1 )
                           F13 = DELTA*C( I, J+2 )
                           F14 = DELTA*C( I, J+3 )
                           !$omp simd reduction(+:F11,F12,F13,F14) &
                           !$omp& simdlen(8) linear(L:1)
                           DO 370 L = LL, LL+LSEC-1
                              F11 = F11 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+1 )
                              F12 = F12 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+2 )
                              F13 = F13 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+3 )
                              F14 = F14 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+4 )
  370                      CONTINUE
                           C( I, J ) = F11
                           C( I, J+1 ) = F12
                           C( I, J+2 ) = F13
                           C( I, J+3 ) = F14
  380                   CONTINUE
                     END IF
  390             CONTINUE
                  IF( UJSEC.LT.JSEC )THEN
                     DO 440 J = JJ+UJSEC, JJ+JSEC-1
                        DO 410 I = II, II+UISEC-1, 4
                           F11 = DELTA*C( I, J )
                           F21 = DELTA*C( I+1, J )
                           F31 = DELTA*C( I+2, J )
                           F41 = DELTA*C( I+3, J )
                           !$omp simd reduction(+:F11,F21,F31,F41) simdlen(8) &
                           !$omp& linear(L:1)
                           DO 400 L = LL, LL+LSEC-1
                              F11 = F11 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+1 )
                              F21 = F21 + T1( L-LL+1, I-II+2 )* &
                                                   T2( L-LL+1, J-JJ+1 )
                              F31 = F31 + T1( L-LL+1, I-II+3 )* &
                                                   T2( L-LL+1, J-JJ+1 )
                              F41 = F41 + T1( L-LL+1, I-II+4 )* &
                                                  T2( L-LL+1, J-JJ+1 )
  400                      CONTINUE
                           C( I,J ) = F11
                           C( I+1, J ) = F21
                           C( I+2, J ) = F31
                           C( I+3, J ) = F41
  410                   CONTINUE
                        DO 430 I = II+UISEC, II+ISEC-1
                           F11 = DELTA*C( I, J )
                           !$omp simd reduction(+:F11) simdlen(8) linear(I:1)
                           DO 420 L = LL, LL+LSEC-1
                              F11 = F11 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+1 )
  420                      CONTINUE
                           C( I, J ) = F11
  430                   CONTINUE
  440                CONTINUE
                  END IF
  450          CONTINUE
  460       CONTINUE
  470    CONTINUE
      END IF

      

      END SUBROUTINE


SUBROUTINE DGEMM_OMP(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,&
                          BETA,C,LDC,INFO,MB,NB,NBT, KB) 
       !dir$ attributes code_align : 32 :: DGEMM_OMP
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: DGEMM_OMP
!     *     .. Scalar Arguments ..
      use omp_lib
      implicit none
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC, INFO
      INTEGER            MB, NB, NBT, KB 
      DOUBLE PRECISION   ALPHA, BETA
!*     .. Array Arguments ..
      !DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A,B,C
#if 0
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
!*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  -- Modified in October-1997.
*     Superscalar GEMM-Based Level 3 BLAS (Version 0.1).
*     Per Ling, Department of Computing Science,
*     Umea University, Sweden.
*
*
*     .. Local Scalars ..
*  -- Modified by Bernard Gingold in 21-02-2021  
#endif
      INTEGER            I, II, ISEC, UISEC, J, JJ, JSEC, UJSEC, &
                         L, LL, LSEC, ULSEC, NROWA, NROWB
      LOGICAL            NOTA, NOTB
      DOUBLE PRECISION   DELTA
      DOUBLE PRECISION   F11, F12, F21, F22, F31, F32, F41, F42
      DOUBLE PRECISION   F13, F14, F23, F24, F33, F34, F43, F44
      DOUBLE PRECISION   TMP0,TMP1,TMP2,TMP3
      DOUBLE PRECISION   C0,C1,C2,C3,C4,C5,C6,C7,C8
                         
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
!*     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!*     .. External Subroutines ..
      !EXTERNAL           XERBLA
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!*     .. User specified parameters for DGEMM ..
!      Moved to subroutine arguments.
     ! PARAMETER        ( MB = 32, NB = 1024, NBT = 96, KB = 32 )
      DOUBLE PRECISION   T1( KB, MB ), T2( KB, NBT )
      !DIR$ ATTRIBUTES ALIGN : 64 :: T1,T2
!*     ..
!*     .. Executable Statements ..
!*
!*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!*     transposed and set NROWA and NROWB as the number of rows of A  and
!*     the number of rows of B respectively.
!*

      NOTA = LSAME( TRANSA, 'N' )
      NOTB = LSAME( TRANSB, 'N' )
      IF ( NOTA ) THEN
         NROWA = M
      ELSE
         NROWA = K
      END IF
      IF ( NOTB ) THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!*
!*     Test the input parameters.
!*
      INFO = 0
      IF( ( .NOT.NOTA ).AND.( .NOT. LSAME( TRANSA,'C' ) ) &
                               .AND.( .NOT.LSAME( TRANSA, 'T' ) ) )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB ).AND.( .NOT.LSAME( TRANSB, 'C' ) ) &
                              .AND.( .NOT.LSAME( TRANSB, 'T' ) ) )THEN
         INFO = 2
      ELSE IF( M.LT.0 )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( K.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
!CALL XERBLA( 'DGEMM ', INFO )
          RETURN
      END IF
!*
!*     Quick return if possible.
!*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
          ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) RETURN
   
!*
!*     And when alpha.eq.zero.
!*
      IF( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) )THEN
         IF( BETA.EQ.ZERO )THEN
            UISEC = M-MOD( M, 4 )
            DO 30, J = 1, N
   !dir$ assume_aligned C:64
   !dir$ vector aligned 
   !dir$ vector vectorlength(8)
   !dir$ vector always
   
               DO 10, I = 1, UISEC, 4
                  C( I, J ) = ZERO
                  C( I+1, J ) = ZERO
                  C( I+2, J ) = ZERO
                  C( I+3, J ) = ZERO
   10          CONTINUE
   !dir$ assume_aligned C:64
   !dir$ vector aligned 
   !dir$ vector vectorlength(8)
   !dir$ vector always
               DO 20, I = UISEC+1, M
                  C( I, J ) = ZERO
   20          CONTINUE
   30       CONTINUE
         ELSE
            UISEC = M-MOD( M, 4 )
            DO 60, J = 1, N
   !dir$ assume_aligned C:64
   !dir$ assume_aligned BETA:64
   !dir$ vector aligned 
   !dir$ vector vectorlength(8)
   !dir$ vector always
               DO 40, I = 1, UISEC, 4
                  C( I, J ) = BETA*C( I, J )
                  C( I+1, J ) = BETA*C( I+1, J )
                  C( I+2, J ) = BETA*C( I+2, J )
                  C( I+3, J ) = BETA*C( I+3, J )
   40          CONTINUE
   !dir$ assume_aligned C:64
   !dir$ assume_aligned BETA:64
   !dir$ vector aligned 
   !dir$ vector vectorlength(8)
   !dir$ vector always
               DO 50, I = UISEC+1, M
                  C( I, J ) = BETA*C( I, J )
   50          CONTINUE
   60       CONTINUE
         END IF
         RETURN
      END IF
      T1(:,:) = ZERO
      T2(:,:) = ZERO
!*
!*     Start the operations.
!*     
      IF (NOTB) THEN
!*
!*        Form  C := alpha*A*B + beta*C or C := alpha*A'*B + beta*C.
!*
         TMP0 = 0.0D0
         TMP1 = 0.0D0
         TMP2 = 0.0D0
         TMP3 = 0.0D0
!$omp parallel do default(none) schedule(dynamic)   &
!$omp private(JJ,JSEC,UJSEC,LL,LSEC,ULSEC,DELTA)    &
!$omp private(II,ISEC,UISEC,L,I)                    &
!$OMP private(F11,F21,F12,F22,F13)                  &
!$omp private(F23,F14,F24,F31,F41,F32,F42)          &
!$omp private(F33,F43,F34,F44)                      &
!$omp private(TMP0,TMP1,TMP2,TMP3)                  &
!$omp shared(N,NB,K,KBONE,BETA,M,MB,ALPHA)          &
!$omp shared(T1,A,C)
         DO 250 JJ = 1, N, NB
            JSEC = MIN( NB, N-JJ+1 )
            UJSEC = JSEC-MOD( JSEC, 4 )
            DO 240 LL = 1, K, KB
               LSEC = MIN( KB, K-LL+1 )
               ULSEC = LSEC-MOD( LSEC, 2 )
!*
!*              Determine if the block of C should be updated with
!*              beta or not.
!*
               DELTA = ONE
               IF( LL.EQ.1 ) DELTA = BETA

               DO 230 II = 1, M, MB
                  ISEC = MIN( MB, M-II+1 )
!*
!*                 T1 := alpha*A' or T1 := alpha*A, copy the transpose
!*                 or the non-transpose of a rectangular block of
!*                 alpha*A to T1.
!*
                  UISEC = ISEC-MOD( ISEC, 2 )
                  IF( NOTA )THEN
                     DO 80, L = LL, LL+ULSEC-1, 2
                        DO 70, I = II, II+UISEC-1, 2
                           T1( L-LL+1, I-II+1 ) = ALPHA*A( I, L )
                           T1( L-LL+2, I-II+1 ) = ALPHA*A( I, L+1 )
                           T1( L-LL+1, I-II+2 ) = ALPHA*A( I+1, L )
                           T1( L-LL+2, I-II+2 ) = ALPHA*A( I+1, L+1 )
   70                   CONTINUE
                        IF( UISEC.LT.ISEC )THEN
                           T1( L-LL+1, ISEC ) = ALPHA*A( II+ISEC-1, L )
                           T1( L-LL+2, ISEC ) = ALPHA*A( II+ISEC-1, L+1)
                                             
                        END IF
   80                CONTINUE
                     IF( ULSEC.LT.LSEC )THEN
                        DO 90, I = II, II+ISEC-1
                           T1( LSEC, I-II+1 ) = ALPHA*A( I, LL+LSEC-1 )
   90                   CONTINUE
                     END IF
                  ELSE
                     DO 110, I = II, II+UISEC-1, 2
                        !dir$ assume_aligned A:64
                        !$omp simd simdlen(4) linear(L:2) aligned(A:64)
                        DO 100, L = LL, LL+ULSEC-1, 2
                           T1( L-LL+1, I-II+1 )  = ALPHA*A( L, I )
                           T1( L-LL+1, I-II+2 ) = ALPHA*A( L, I+1 )
                           T1( L-LL+2, I-II+1 ) = ALPHA*A( L+1, I )
                           T1( L-LL+2, I-II+2 ) = ALPHA*A( L+1, I+1 )
  100                   CONTINUE
                        IF( ULSEC.LT.LSEC )THEN
                           T1( LSEC, I-II+1 ) = ALPHA*A( LL+LSEC-1, I )
                           T1( LSEC, I-II+2 ) = ALPHA*A( LL+LSEC-1, I+1)
                                             
                        END IF
  110                CONTINUE
                     IF( UISEC.LT.ISEC )THEN
                        !dir$ assume_aligned A:64
                        !$omp simd simdlen(8) linear(L:2) aligned(A:64)    
                        DO 120, L = LL, LL+LSEC-1
                           T1( L-LL+1, ISEC ) = ALPHA*A( L, II+ISEC-1 )
  120                   CONTINUE
                     END IF
                  END IF
!*
!*                 C := T1'*B + beta*C, update a rectangular block
!*                 of C using 4 by 4 unrolling.
!*
                  UISEC = ISEC-MOD( ISEC, 4 )
                  DO 170 J = JJ, JJ+UJSEC-1, 4
                     !dir$ assume_aligned C:64
                     !$omp simd simdlen(8) 
                     DO 140 I = II, II+UISEC-1, 4
                        F11 = DELTA*C( I,J )
                        F21 = DELTA*C( I+1,J )
                        F12 = DELTA*C( I,J+1 )
                        F22 = DELTA*C( I+1,J+1 )
                        F13 = DELTA*C( I,J+2 )
                        F23 = DELTA*C( I+1,J+2 )
                        F14 = DELTA*C( I,J+3 )
                        F24 = DELTA*C( I+1,J+3 )
                        F31 = DELTA*C( I+2,J )
                        F41 = DELTA*C( I+3,J )
                        F32 = DELTA*C( I+2,J+1 )
                        F42 = DELTA*C( I+3,J+1 )
                        F33 = DELTA*C( I+2,J+2 )
                        F43 = DELTA*C( I+3,J+2 )
                        F34 = DELTA*C( I+2,J+3 )
                        F44 = DELTA*C( I+3,J+3 )
                        !dir$ assume_aligned T1:64
                        !dir$ assume_aligned B:64
                        !$omp  simd simdlen(8) linear(L:1) private(TMP0,TMP1,TMP2,TMP3) &
                        !$omp& reduction(+:F11,F21,F12,F22,F13,F23,F14,F24,F31,F41,F32, &
                        !$omp&             F42,F33,F43,F34,F44)
                        DO 130 L = LL, LL+LSEC-1
                           call _mm_prefetch(T1(L+32,I),FOR_K_PREFETCH_T1,.false.)
                           TMP0 = B(L,J)
                           TMP1 = B(L,J+1)
                           TMP2 = B(L,J+2)
                           TMP3 = B(L,J+3)
                           F11 = F11 + T1( L-LL+1, I-II+1 )*TMP0
                           F21 = F21 + T1( L-LL+1, I-II+2 )*TMP0
                           F12 = F12 + T1( L-LL+1, I-II+1 )*TMP1
                           F22 = F22 + T1( L-LL+1, I-II+2 )*TMP1
                           F13 = F13 + T1( L-LL+1, I-II+1 )*TMP2
                           F23 = F23 + T1( L-LL+1, I-II+2 )*TMP2
                           F14 = F14 + T1( L-LL+1, I-II+1 )*TMP3
                           F24 = F24 + T1( L-LL+1, I-II+2 )*TMP3
                           F31 = F31 + T1( L-LL+1, I-II+3 )*TMP0
                           F41 = F41 + T1( L-LL+1, I-II+4 )*TMP0
                           F32 = F32 + T1( L-LL+1, I-II+3 )*TMP1
                           F42 = F42 + T1( L-LL+1, I-II+4 )*TMP1
                           F33 = F33 + T1( L-LL+1, I-II+3 )*TMP2
                           F43 = F43 + T1( L-LL+1, I-II+4 )*TMP2
                           F34 = F34 + T1( L-LL+1, I-II+3 )*TMP3
                           F44 = F44 + T1( L-LL+1, I-II+4 )*TMP3
  130                   CONTINUE
                        C( I,J ) = F11
                        C( I+1, J ) = F21
                        C( I, J+1 ) = F12
                        C( I+1, J+1 ) = F22
                        C( I, J+2 ) = F13
                        C( I+1, J+2 ) = F23
                        C( I, J+3 ) = F14
                        C( I+1, J+3 ) = F24
                        C( I+2, J ) = F31
                        C( I+3, J ) = F41
                        C( I+2, J+1 ) = F32
                        C( I+3, J+1 ) = F42
                        C( I+2, J+2 ) = F33
                        C( I+3, J+2 ) = F43
                        C( I+2, J+3 ) = F34
                        C( I+3, J+3 ) = F44
  140                CONTINUE
                     IF( UISEC.LT.ISEC )THEN
                        !dir$ assume_aligned C:64
                        !$omp simd simdlen(8) linear(I:1) 
                        DO 160 I = II+UISEC, II+ISEC-1
                           F11 = DELTA*C( I, J )
                           F12 = DELTA*C( I, J+1 )
                           F13 = DELTA*C( I, J+2 )
                           F14 = DELTA*C( I, J+3 )
                           !dir$ assume_aligned T1:64
                           !dir$ assume_aligned B:64
                           !$omp simd simdlen(8) 
                           !$omp& reduction(+:F11,F12,F13,F14) aligned(B:64)
                           DO 150 L = LL, LL+LSEC-1
                              TMP0 = T1( L-LL+1, I-II+1 )*B(L,J)
                              F11 = F11 + TMP0
                              TMP1 =  T1( L-LL+1, I-II+1 )*B( L, J+1 )
                              F12 = F12 + TMP1
                              TMP2 = T1( L-LL+1, I-II+1 )*B( L, J+2 )
                              F13 = F13 + TMP2
                              TMP3 = T1( L-LL+1, I-II+1 )*B( L, J+3 )
                              F14 = F14 + TMP3                             
  150                      CONTINUE
                           C( I, J ) = F11
                           C( I, J+1 ) = F12
                           C( I, J+2 ) = F13
                           C( I, J+3 ) = F14
  160                   CONTINUE
                     END IF
  170             CONTINUE
                  IF( UJSEC.LT.JSEC )THEN
                     TMP0 = 0.0D0
                     TMP1 = 0.0D0
                     TMP2 = 0.0D0
                     TMP3 = 0.0D0
                     DO 220 J = JJ+UJSEC, JJ+JSEC-1
                        !dir$ assume_aligned C:64
                        !$omp simd simdlen(8) linear(I:4) 
                        DO 190 I = II, II+UISEC-1, 4
                           F11 = DELTA*C( I,J )
                           F21 = DELTA*C( I+1, J )
                           F31 = DELTA*C( I+2, J )
                           F41 = DELTA*C( I+3, J )
                           !dir$ assume_aligned T1:64
                           !dir$ assume_aligned B:64
                           !$omp simd simdlen(8) 
                           !$omp& reduction(+:F11,F12,F13,F14) aligned(B:64)
                           DO 180 L = LL, LL+LSEC-1
                              TMP0 = T1( L-LL+1, I-II+1 )*B( L, J )
                              F11  = F11 + TMP0
                              TMP1 = T1( L-LL+1, I-II+2 )*B( L, J )
                              F21  = F21 + TMP1
                              TMP2 = T1( L-LL+1, I-II+3 )*B( L, J )
                              F31  = F31 + TMP2
                              TMP3 = T1( L-LL+1, I-II+4 )*B( L, J )
                              F41  + F41 + TMP3
  180                      CONTINUE
                           C( I,J ) = F11
                           C( I+1, J ) = F21
                           C( I+2, J ) = F31
                           C( I+3, J ) = F41
  190                   CONTINUE
                        DO 210 I = II+UISEC, II+ISEC-1
                           F11 = DELTA*C( I, J )
                           !dir$ assume_aligned T1:64
                           !dir$ assume_aligned B:64
                           !$omp simd simdlen(8) &
                           !$omp& reduction(+:F11) aligned(B:64)
                           DO 200 L = LL, LL+LSEC-1
                              F11 = F11 + T1( L-LL+1, I-II+1 )*B( L, J )
  200                      CONTINUE
                           C( I, J ) = F11
  210                   CONTINUE
  220                CONTINUE
                  END IF
  230          CONTINUE
  240       CONTINUE
  250    CONTINUE
!$omp end parallel do
      ELSE
!*
!*        Form  C := alpha*A*B' + beta*C or C := alpha*A'*B' + beta*C.
!*
         C0  = ZERO
         C1  = ZERO
         C2  = ZERO
         C3  = ZERO
         C4  = ZERO
         C5  = ZERO
         C6  = ZERO
         C7  = ZERO
         C8  = ZERO
!$omp parallel do default(none) schedule(dynamic) &
!$omp firstprivate(C0,C1,C2,C3,C4,C5,C6,C7,C8)    &
!$omp private(JJ,JSEC,LL,LSEC,DELTA,ULSEC,UJSEC)  &
!$omp private(II,J,ISEC,L,I,UISEC)                &
!$omp private(F11,F21,F12,F22,F13,F23,F14,F24)    &
!$omp private(F31,F41,F32,F42,F33,F43,F34,F44)    &
!$omp shared(N,NBT,K,KB,NOTA,T2,ALPHA,N,C)        &
!$omp shared(M,MB)      
         DO 470 JJ = 1, N, NBT
            JSEC = MIN( NBT, N-JJ+1 )
            DO 460 LL = 1, K, KB
               LSEC = MIN( KB, K-LL+1 )
!*
!*              Determine if the block of C should be updated with
!*              beta or not.
!*
               DELTA = ONE
               IF( LL.EQ.1 ) DELTA = BETA
!*
!*              T2 := alpha*B', copy the transpose of a rectangular
!*              block of alpha*A to T2.
!*
               ULSEC = LSEC-MOD( LSEC, 2 )
               UJSEC = JSEC-MOD( JSEC, 2 )
               DO 270, L = LL, LL+ULSEC-1, 2
                  DO 260, J = JJ, JJ+UJSEC-1, 2
                     T2( L-LL+1, J-JJ+1 ) = ALPHA*B( J, L )
                     T2( L-LL+2, J-JJ+1 ) = ALPHA*B( J, L+1 )
                     T2( L-LL+1, J-JJ+2 ) = ALPHA*B( J+1, L )
                     T2( L-LL+2, J-JJ+2 ) = ALPHA*B( J+1, L+1 )
  260             CONTINUE
                  IF( UJSEC.LT.JSEC )THEN
                     T2( L-LL+1, JSEC ) = ALPHA*B( JJ+JSEC-1, L )
                     T2( L-LL+2, JSEC ) = ALPHA*B( JJ+JSEC-1, L+1 )
                  END IF
  270          CONTINUE
               IF( ULSEC.LT.LSEC )THEN
                  DO 280, J = JJ, JJ+JSEC-1
                     T2( LSEC, J-JJ+1 ) = ALPHA*B( J, LL+LSEC-1 )
  280             CONTINUE
               END IF
!*
               UJSEC = JSEC-MOD( JSEC, 4 )
               DO 450 II = 1, M, MB
                  ISEC = MIN( MB, M-II+1 )
!*
!*                 T1 := alpha*A' or T1 := alpha*A, copy the transpose
!*                 or the non-transpose of a rectangular block of
!*                 alpha*A to T1.
!*
                  UISEC = ISEC-MOD( ISEC, 2 )
                  IF( NOTA )THEN
                     DO 300, L = LL, LL+ULSEC-1, 2
                        DO 290, I = II, II+UISEC-1, 2
                           T1( L-LL+1, I-II+1 ) = A( I, L )
                           T1( L-LL+2, I-II+1 ) = A( I, L+1 )
                           T1( L-LL+1, I-II+2 ) = A( I+1, L )
                           T1( L-LL+2, I-II+2 ) = A( I+1, L+1 )
  290                   CONTINUE
                        IF( UISEC.LT.ISEC )THEN
                           T1( L-LL+1, ISEC ) = A( II+ISEC-1, L )
                           T1( L-LL+2, ISEC ) = A( II+ISEC-1, L+1 )
                        END IF
  300                CONTINUE
                     IF( ULSEC.LT.LSEC )THEN
                        DO 310, I = II, II+ISEC-1
                           T1( LSEC, I-II+1 ) = A( I, LL+LSEC-1 )
  310                   CONTINUE
                     END IF
                  ELSE
                     DO 330, I = II, II+UISEC-1, 2
                        !dir$ assume_aligned T1:64
                        !dir$ assume_aligned A:64
                        !$omp simd simdlen(8) linear(L:2)
                        DO 320, L = LL, LL+ULSEC-1, 2
                           T1( L-LL+1, I-II+1 ) = A( L, I )
                           T1( L-LL+1, I-II+2 ) = A( L, I+1 )
                           T1( L-LL+2, I-II+1 ) = A( L+1, I )
                           T1( L-LL+2, I-II+2 ) = A( L+1, I+1 )
  320                   CONTINUE
                        IF( ULSEC.LT.LSEC )THEN
                           T1( LSEC, I-II+1 ) = A( LL+LSEC-1, I )
                           T1( LSEC, I-II+2 ) = A( LL+LSEC-1, I+1 )
                        END IF
  330                CONTINUE
                     IF( UISEC.LT.ISEC )THEN
                        !dir$ assume_aligned T1:64
                        !dir$ assume_aligned A:64
                        !$omp simd simdlen(8) linear(L:2) 
                        DO 340, L = LL, LL+LSEC-1
                           T1( L-LL+1, ISEC ) = A( L, II+ISEC-1 )
  340                   CONTINUE
                     END IF
                  END IF
!*
!*                 C := T1'*B + beta*C, update a rectangular block
!*                 of C using 4 by 4 unrolling.
!*
                  UISEC = ISEC-MOD( ISEC, 4 )
                  DO 390 J = JJ, JJ+UJSEC-1, 4
                     !dir$ assume_aligned C:64
                     !$omp simd simdlen(8) aligned(C:64) linear(I:4) &
                     
                     DO 360 I = II, II+UISEC-1, 4
                        F11 = DELTA*C( I,J )
                        F21 = DELTA*C( I+1,J )
                        F12 = DELTA*C( I,J+1 )
                        F22 = DELTA*C( I+1,J+1 )
                        F13 = DELTA*C( I,J+2 )
                        F23 = DELTA*C( I+1,J+2 )
                        F14 = DELTA*C( I,J+3 )
                        F24 = DELTA*C( I+1,J+3 )
                        F31 = DELTA*C( I+2,J )
                        F41 = DELTA*C( I+3,J )
                        F32 = DELTA*C( I+2,J+1 )
                        F42 = DELTA*C( I+3,J+1 )
                        F33 = DELTA*C( I+2,J+2 )
                        F43 = DELTA*C( I+3,J+2 )
                        F34 = DELTA*C( I+2,J+3 )
                        F44 = DELTA*C( I+3,J+3 )
                        !dir$ assume_aligned T1:64
                        !dir$ assume_aligned T2:64
                        !dir$ simd simdlen(8)  linear(I:1) 
                        DO 350 L = LL, LL+LSEC-1
                           C0  = T1(L-LL+1,I-II+1)
                           C1  = T2(L-LL+1,J-JJ+1)
                           F11 = F11 + C0 * C1
                           C2  = T1(L-LL+1,I-II+2)
                           F21 = F21 + C2 * C1
                           C3 = T2(L-LL+1,J-JJ+2)
                           F12 = F12 + C0 * C3
                           F22 = F22 + C2 * C3
                           C4 = T1(L-LL+1,I-II+3)
                           C5 = T2(L-LL+1,J-JJ+3)
                           F13 = F13 + C0 * C4
                           F23 = F23 + C2 * C4
                           C6 = T1(L-LL+1,I-II+4)
                           C7 = T2(L-LL+1,J-JJ+4)
                           F14 = F14 + C0 * C7
                           F24 = F24 + C2 * C7
                           F31 = F31 + C4 * C1
                           F41 = F41 + C6 * C1
                           F32 = F32 + C4 * C3
                           F42 = F42 + C6 * C3
                           F33 = F33 + C4 * C5
                           F43 = F43 + C6 * C5
                           F34 = F34 + C4 * C7
                           F44 = F44 + C6 * C7
                           
  350                   CONTINUE
                        C( I,J ) = F11
                        C( I+1, J ) = F21
                        C( I, J+1 ) = F12
                        C( I+1, J+1 ) = F22
                        C( I, J+2 ) = F13
                        C( I+1, J+2 ) = F23
                        C( I, J+3 ) = F14
                        C( I+1, J+3 ) = F24
                        C( I+2, J ) = F31
                        C( I+3, J ) = F41
                        C( I+2, J+1 ) = F32
                        C( I+3, J+1 ) = F42
                        C( I+2, J+2 ) = F33
                        C( I+3, J+2 ) = F43
                        C( I+2, J+3 ) = F34
                        C( I+3, J+3 ) = F44
  360                CONTINUE
                     IF( UISEC.LT.ISEC )THEN
                        !dir$ assume_aligned C:64
                        !$omp simd simdlen(4) linear(I:1) 
                        DO 380 I = II+UISEC, II+ISEC-1
                           F11 = DELTA*C( I, J )
                           F12 = DELTA*C( I, J+1 )
                           F13 = DELTA*C( I, J+2 )
                           F14 = DELTA*C( I, J+3 )
                           !dir$ assume_aligned T1:64
                           !dir$ assume_aligned T2:64
                           !dir$ simd simdlen(8) reduction(+:F11,F12,F13,F14) linear(L:1)
                           DO 370 L = LL, LL+LSEC-1
                              F11 = F11 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+1 )
                              F12 = F12 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+2 )
                              F13 = F13 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+3 )
                              F14 = F14 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+4 )
  370                      CONTINUE
                           C( I, J ) = F11
                           C( I, J+1 ) = F12
                           C( I, J+2 ) = F13
                           C( I, J+3 ) = F14
  380                   CONTINUE
                     END IF
  390             CONTINUE
                  IF( UJSEC.LT.JSEC )THEN
                     DO 440 J = JJ+UJSEC, JJ+JSEC-1
                        !dir$ assume_aligned C:64
                        !$omp simd simdlen(8) linear(I:4) 
                        DO 410 I = II, II+UISEC-1, 4
                           F11 = DELTA*C( I, J )
                           F21 = DELTA*C( I+1, J )
                           F31 = DELTA*C( I+2, J )
                           F41 = DELTA*C( I+3, J )
                           !dir$ assume_aligned T1:64
                           !dir$ assume_aligned T2:64
                           !dir$ simd simdlen(8) reduction(+:F11,F12,F13,F14) linear(L:1)
                           DO 400 L = LL, LL+LSEC-1
                              F11 = F11 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+1 )
                              F21 = F21 + T1( L-LL+1, I-II+2 )* &
                                                   T2( L-LL+1, J-JJ+1 )
                              F31 = F31 + T1( L-LL+1, I-II+3 )* &
                                                   T2( L-LL+1, J-JJ+1 )
                              F41 = F41 + T1( L-LL+1, I-II+4 )* &
                                                  T2( L-LL+1, J-JJ+1 )
  400                      CONTINUE
                           C( I,J ) = F11
                           C( I+1, J ) = F21
                           C( I+2, J ) = F31
                           C( I+3, J ) = F41
  410                   CONTINUE
                        DO 430 I = II+UISEC, II+ISEC-1
                           F11 = DELTA*C( I, J )
                           !dir$ assume_aligned T1:64
                           !dir$ assume_aligned T2:64
                           !dir$ simd simdlen(8) reduction(+:F11) linear(L:1)
                           DO 420 L = LL, LL+LSEC-1
                              F11 = F11 + T1( L-LL+1, I-II+1 )* &
                                                   T2( L-LL+1, J-JJ+1 )
  420                      CONTINUE
                           C( I, J ) = F11
  430                   CONTINUE
  440                CONTINUE
                  END IF
  450          CONTINUE
  460       CONTINUE
  470    CONTINUE
!$omp end parallel do
      END IF
      

      END SUBROUTINE




      

SUBROUTINE DSYMM(SIDE,UPLO M,N,ALPHA,A,LDA,B,LDB, &
BETA,C,LDC,RCB,MB,NB,NBT,KB) 
       !dir$ attributes code_align : 32 :: DSYMM
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: DSYMM
      use omp_lib
      implicit none
      CHARACTER*1        SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      INTEGER            RCB,MB,NB,NBT,KB
      DOUBLE PRECISION   ALPHA, BETA
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
#if 0
*     ..
*
*  Purpose
*  =======
*
*  DSYMM  performs one of the matrix-matrix operations
*
*     C := alpha*A*B + beta*C,
*
*  or
*
*     C := alpha*B*A + beta*C,
*
*  where alpha and beta are scalars,  A is a symmetric matrix and  B and
*  C are  m by n matrices.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE  specifies whether  the  symmetric matrix  A
*           appears on the  left or right  in the  operation as follows:
*
*              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
*
*              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of  the  symmetric  matrix   A  is  to  be
*           referenced as follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of the
*                                  symmetric matrix is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of the
*                                  symmetric matrix is to be referenced.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies the number of rows of the matrix  C.
*           M  must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix C.
*           N  must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
*           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
*           the array  A  must contain the  symmetric matrix,  such that
*           when  UPLO = 'U' or 'u', the leading m by m upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  symmetric matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  m by m  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  symmetric
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
*           the array  A  must contain the  symmetric matrix,  such that
*           when  UPLO = 'U' or 'u', the leading n by n upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  symmetric matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  n by n  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  symmetric
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry, the leading  m by n part of the array  B  must
*           contain the matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n updated
*           matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  -- Rewritten in December-1993.
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*  -- Rewritten in Mars-1995.
*     Superscalar GEMM-Based Level 3 BLAS (Version 0.1).
*     Per Ling, Department of Computing Science,
*     University of Umea, Sweden.
*
*
*     .. Local Scalars ..
#endif
      INTEGER            INFO, NROWA
      LOGICAL            LSIDE, UPPER
      INTEGER            I, J, II, IX, ISEC, UISEC
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
!*     .. External Functions ..
!      LOGICAL            LSAME
      EXTERNAL           LSAME
!*     .. External Subroutines ..
!      EXTERNAL           XERBLA
!      EXTERNAL           DGEMM
!*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION TMP0,TMP1,TMP2,TMP3
!*     .. User specified parameters for DSYMM ..
!      INTEGER            RCB
!      PARAMETER        ( RCB = 96 )
!*     .. Local Arrays ..
      DOUBLE PRECISION   T1( RCB, RCB )
!*     ..
!*     .. Executable Statements ..
!*
!*     Test the input parameters.
!*
      LSIDE = LSAME( SIDE, 'L' )
      UPPER = LSAME( UPLO, 'U' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      INFO = 0
      IF( ( .NOT.LSIDE ).AND.( .NOT.LSAME( SIDE, 'R' ) ) )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER ).AND.( .NOT.LSAME( UPLO, 'L' ) ) )THEN
         INFO = 2
      ELSE IF( M.LT.0 )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, M ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, M ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         !CALL XERBLA( 'DSYMM ', INFO )
         RETURN
      END IF
!*
!*     Quick return if possible.
!*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
          ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) RETURN

      T0 = 0.0D+0
      T1 = 0.0D+0
      T2 = 0.0D+0
      T3 = 0.0D+0
!*
!*     And when alpha.eq.zero.
!*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            UISEC = M-MOD( M, 4 )
            DO 30, J = 1, N
           !dir$ assume_aligned C:64
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector vectorlength(8)
               DO 10, I = 1, UISEC, 4
                  C( I, J ) = ZERO
                  C( I+1, J ) = ZERO
                  C( I+2, J ) = ZERO
                  C( I+3, J ) = ZERO
   10          CONTINUE
           !dir$ assume_aligned C:64
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector vectorlength(8)
               DO 20, I = UISEC+1, M
                  C( I, J ) = ZERO
   20          CONTINUE
   30       CONTINUE
         ELSE
            UISEC = M-MOD( M, 4 )
            DO 60, J = 1, N
           !dir$ assume_aligned C:64
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector vectorlength(8)
               DO 40, I = 1, UISEC, 4
                  C( I, J ) = BETA*C( I, J )
                  C( I+1, J ) = BETA*C( I+1, J )
                  C( I+2, J ) = BETA*C( I+2, J )
                  C( I+3, J ) = BETA*C( I+3, J )
   40          CONTINUE
           !dir$ assume_aligned C:64
           !dir$ vector aligned
           !dir$ vector always
           !dir$ vector vectorlength(8)
               DO 50, I = UISEC+1, M
                  C( I, J ) = BETA*C( I, J )
   50          CONTINUE
   60       CONTINUE
         END IF
         RETURN
      END IF
!*
!*     Start the operations.
!*
      IF( LSIDE )THEN
         IF( UPPER )THEN
!*
!*           Form  C := alpha*A*B + beta*C. Left, Upper.
!1*
!$omp parallel do default(none) schedule(dynamic)   &
!$omp private(II,ISEC,J,UISEC,TMP0,TMP1,TMP2,TMP3)  &
!$omp shared(M,RCB,A,B,ALPHA,T1,N,C,ONE,LDB,BETA)     &
!$omp shared(LDC,MB,NB,NBT,KB)
            DO 90, II = 1, M, RCB
               ISEC = MIN( RCB, M-II+1 )
!*
!*              T1 := A, the upper triangular part of a square diagonal
!*              block of A is copied to upper pangular part of T1 and
!*              the transpose of the strictly upper triangular part of
!*              the block of A is copied to the strictly lower
!1*              triangular part of T1.
!*
               DO 80, J = II+ISEC-2, II-1, -2
                  UISEC = J-II-MOD( J-II, 2 )
                  !dir$ assume_aligned T1:64
                  !dir$ assume_aligned A:64
                  !dir$ vector aligned
                  !dir$ vector vectorlength(8)
                  !dir$ vector always
                  DO 70, I = II, II+UISEC-1, 2
                     TMP0 = A(I,J)
                     TMP1 = A(I+1,J)
                     TMP2 = A(I,J+1)
                     TMP3 = A(I+1,J+1)
                     T1( I-II+1, J-II+1 ) = ALPHA*TMP0
                     T1( I-II+2, J-II+1 ) = ALPHA*TMP1
                     T1( I-II+1, J-II+2 ) = ALPHA*TMP3
                     T1( I-II+2, J-II+2 ) = ALPHA*TMP4
                     T1( J-II+1, I-II+1 ) = ALPHA*TMP0
                     T1( J-II+1, I-II+2 ) = ALPHA*TMP1
                     T1( J-II+2, I-II+1 ) = ALPHA*TMP2
                     T1( J-II+2, I-II+2 ) = ALPHA*TMP3
   70             CONTINUE
                  IF( MOD( J-II, 2 ).EQ.1 )THEN
                     TMP0 = A(J-1,J)
                     TMP1 = A(J-1,J+1)
                     TMP2 = A(J,J+1)
                     T1( J-II, J-II+1 ) = ALPHA*TMP0
                     T1( J-II+1, J-II+1 ) = ALPHA*A( J, J )
                     T1( J-II, J-II+2 ) = ALPHA*TMP1
                     T1( J-II+1, J-II+2 ) = ALPHA*TMP2
                     T1( J-II+1, J-II ) = ALPHA*TMP0
                     T1( J-II+2, J-II ) = ALPHA*TMP1
                     T1( J-II+2, J-II+1 ) = ALPHA*TMP2
                  ELSE IF( J.GE.II )THEN
                     T1( J-II+1, J-II+1 ) = ALPHA*A( J, J )
                     T1( J-II+1, J-II+2 ) = ALPHA*A( J, J+1 )
                     T1( J-II+2, J-II+1 ) = ALPHA*A( J, J+1 )
                  END IF
                  T1( J-II+2, J-II+2 ) = ALPHA*A( J+1, J+1 )
   80          CONTINUE
!*
!*              C := T1'*B + beta*C, general matrix multiplication
!*              involving the symmetric diagonal block of A stored
!*              as a full matrix block in T1.
!*
               CALL DGEMM ( 'T', 'N', ISEC, N, ISEC, ONE, &
                                     T1( 1, 1 ), RCB, B( II, 1 ), LDB, &
                                                BETA, C( II, 1 ), LDC,MB,NB,NBT,KB )
               IF( II.GT.1 )THEN
!*
!*                 C := alpha*A'*B + C, general matrix multiplication
!*                 involving the transpose of a rectangular block of A.
!*
                  CALL DGEMM ( 'T', 'N', ISEC, N, II-1, ALPHA, &
                                       A( 1, II ), LDA, B( 1, 1 ), LDB, &
                                                 ONE, C( II, 1 ), LDC,MB,NB,NBT,KB )
               END IF
               IF( II+ISEC.LE.M )THEN
!*
!*                 C := alpha*A*B + C, general matrix multiplication
!*                 involving a rectangular block of A.
!*
                  CALL DGEMM ( 'N', 'N', ISEC, N, M-II-ISEC+1, ALPHA, &
                           A( II, II+ISEC ), LDA, B( II+ISEC, 1 ), LDB, &
                                                 ONE, C( II, 1 ), LDC MB,NB,NBT,KB )
               END IF
   90       CONTINUE
!$omp end parallel do
         ELSE
!*
!*           Form  C := alpha*A*B + beta*C. Left, Lower.
!*
!$omp parallel do default(none) schedule(dynamic) &
!$omp private(IX,II,ISEC,J,UISEC,TMP0,TMP1,TMP2,TMP3) &
!$omp private(I)                                      &
!$omp shared(M,RCB,ALPHA,T1,A1,N,ONE,LDB,BETA,LDC) &
!$omp shared(MB,NB,NBT,KB,A,B,C)
            DO 120, IX = M, 1, -RCB
               II = MAX( 1, IX-RCB+1 )
               ISEC = IX-II+1
!*
!*              T1 := A, the lower triangular part of a square diagonal
!*              block of A is copied to lower pangular part of T1 and
!*              the transpose of the strictly lower triangular part of
!*              the block of A is copied to the strictly upper
!*              triangular part of T1.
!*
               DO 110, J = II, II+ISEC-1, 2
                  UISEC = II+ISEC-J-2-MOD( II+ISEC-J-2, 2 )
                  T1( J-II+1, J-II+1 ) = ALPHA*A( J, J )
                  IF( MOD( II+ISEC-J-2, 2 ).EQ.0 )THEN
                     TMP0 = A(J+1,J)
                     TMP1 = A(J+1,J+1)
                     T1( J-II+2, J-II+1 ) = ALPHA*TMP0
                     T1( J-II+1, J-II+2 ) = ALPHA*TMP0
                     T1( J-II+2, J-II+2 ) = ALPHA*TMP1
                  ELSE IF( J.LE.II+ISEC-3 )THEN
                     TMP0 = A(J+1,J)
                     TMP1 = A(J+2,J)
                     TMP2 = A(J+1,J+1)
                     TMP3 = A(J+2,J+1)
                     T1( J-II+2, J-II+1 ) = ALPHA*TMP0
                     T1( J-II+1, J-II+2 ) = ALPHA*TMP0
                     T1( J-II+3, J-II+1 ) = ALPHA*TMP1
                     T1( J-II+1, J-II+3 ) = ALPHA*TMP1
                     T1( J-II+2, J-II+2 ) = ALPHA*TMP2
                     T1( J-II+3, J-II+2 ) = ALPHA*TMP3
                     T1( J-II+2, J-II+3 ) = ALPHA*TMP3
                  END IF
                  DO 100 I = II+ISEC-UISEC, II+ISEC-1, 2
                     TMP0 = A(I,J)
                     TMP1 = A(I+1,J)
                     TMP2 = A(I,J+1)
                     TMP3 = A(I+1,J+1)
                     T1( I-II+1, J-II+1 ) = ALPHA*TMP0
                     T1( I-II+2, J-II+1 ) = ALPHA*TMP1
                     T1( I-II+1, J-II+2 ) = ALPHA*TMP2
                     T1( I-II+2, J-II+2 ) = ALPHA*TMP3
                     T1( J-II+1, I-II+1 ) = ALPHA*TMP0
                     T1( J-II+1, I-II+2 ) = ALPHA*TMP1
                     T1( J-II+2, I-II+1 ) = ALPHA*TMP2
                     T1( J-II+2, I-II+2 ) = ALPHA*TMP3
  100             CONTINUE
  110          CONTINUE
!*
!*              C := T1'*B + beta*C, general matrix multiplication
!*              involving the symmetric diagonal block of A stored
!*              as a full matrix block in T1.
!*
               CALL DGEMM ( 'T', 'N', ISEC, N, ISEC, ONE, &
                                      T1( 1, 1 ), RCB, B( II, 1 ), LDB, &
                                                BETA, C( II, 1 ), LDC, MB,NB,NBT,KB )
               IF( II.GT.1 )THEN
!*
!*                 C := alpha*A'*B + C, general matrix multiplication
!*                 involving the transpose of a rectangular block of A.
!*
                  CALL DGEMM ( 'N', 'N', ISEC, N, II-1, ALPHA, &
                                       A( II, 1 ), LDA, B( 1, 1 ), LDB, &
                                                 ONE, C( II, 1 ), LDC, MB,NB,NBT,KB )
               END IF
               IF( II+ISEC.LE.M )THEN
!*
!*                 C := alpha*A*B + C, general matrix multiplication
!*                 involving a rectangular block of A.
!*
                  CALL DGEMM ( 'T', 'N', ISEC, N, M-II-ISEC+1, ALPHA, &
                           A( II+ISEC, II ), LDA, B( II+ISEC, 1 ), LDB, &
                                                 ONE, C( II, 1 ), LDC, MB,NB,NBT,KB )
               END IF
  120       CONTINUE
         END IF
!$omp end parallel do
      ELSE
         IF( UPPER )THEN
!*
!*           Form  C := alpha*B*A + beta*C. Right, Upper.
!*
!$omp parallel do default(none) schedule(dynamic) &
!$omp private(II,ISEC,J,UISEC,I,TMP0,TMP1,TMP2,TMP3)  &
!$omp shared(N,RCB,ALPHA,A,T1,ONE,A,B,C,LDB,BETA,LDC) &
!$omp shared(MB,NB,NBT,KB)
            DO 150, II = 1, N, RCB
               ISEC = MIN( RCB, N-II+1 )
!*
!*              T1 := A, the upper triangular part of a square diagonal
!*              block of A is copied to upper pangular part of T1 and
!*              the transpose of the strictly upper triangular part of
!*              the block of A is copied to the strictly lower
!*              triangular part of T1.
!*
               DO 140, J = II+ISEC-2, II-1, -2
                  UISEC = J-II-MOD( J-II, 2 )
                  !dir$ assume_aligned T1:64
                  !dir$ assume_aligned A:64
                  !dir$ vector aligned
                  !dir$ vector vectorlength(8)
                  !dir$ vector always
                  DO 130, I = II, II+UISEC-1, 2
                     TMP0 = A(I,J)
                     TMP1 = A(I+1,J)
                     TMP2 = A(I,J+1)
                     TMP3 = A(I+1,J+1)
                     T1( I-II+1, J-II+1 ) = ALPHA*TMP0
                     T1( I-II+2, J-II+1 ) = ALPHA*TMP1
                     T1( I-II+1, J-II+2 ) = ALPHA*TMP2
                     T1( I-II+2, J-II+2 ) = ALPHA*TMP3
                     T1( J-II+1, I-II+1 ) = ALPHA*TMP0
                     T1( J-II+1, I-II+2 ) = ALPHA*TMP1
                     T1( J-II+2, I-II+1 ) = ALPHA*TMP2
                     T1( J-II+2, I-II+2 ) = ALPHA*TMP3
  130             CONTINUE
                     IF( MOD( J-II, 2 ).EQ.1 )THEN
                     T1( J-II, J-II+1 ) = ALPHA*A( J-1, J )
                     T1( J-II+1, J-II+1 ) = ALPHA*A( J, J )
                     T1( J-II, J-II+2 ) = ALPHA*A( J-1, J+1 )
                     T1( J-II+1, J-II+2 ) = ALPHA*A( J, J+1 )
                     T1( J-II+1, J-II ) = ALPHA*A( J-1, J )
                     T1( J-II+2, J-II ) = ALPHA*A( J-1, J+1 )
                     T1( J-II+2, J-II+1 ) = ALPHA*A( J, J+1 )
                  ELSE IF( J.GE.II )THEN
                     T1( J-II+1, J-II+1 ) = ALPHA*A( J, J )
                     T1( J-II+1, J-II+2 ) = ALPHA*A( J, J+1 )
                     T1( J-II+2, J-II+1 ) = ALPHA*A( J, J+1 )
                  END IF
                  T1( J-II+2, J-II+2 ) = ALPHA*A( J+1, J+1 )
  140          CONTINUE
!*
!*              C := T1'*B + beta*C, general matrix multiplication
!*              involving the symmetric diagonal block of A stored
!*              as a full matrix block in T1.
!*
               CALL DGEMM ( 'N', 'N', M, ISEC, ISEC, ONE, &
                                      B( 1, II ), LDB, T1( 1, 1 ), RCB, &
                                                BETA, C( 1, II ), LDC,MB,NB,NBT,KB )
              IF( II.GT.1 )THEN
!*
!*                 C := alpha*B*A + C, general matrix multiply
!*                 involving a rectangular block of A.
!*
                  CALL DGEMM ( 'N', 'N', M, ISEC, II-1, ALPHA, &
                                       B( 1, 1 ), LDB, A( 1, II ), LDA, &
                                                 ONE, C( 1, II ), LDC, MB,NB,NBT,KB )
               END IF
               IF( II+ISEC.LE.N )THEN
!*
!*                 C := alpha*B*A' + C, general matrix multiply involving
!*                 the transpose of a rectangular block of A.
!*
                  CALL DGEMM ( 'N', 'T', M, ISEC, N-II-ISEC+1, ALPHA, &
                           B( 1, II+ISEC ), LDB, A( II, II+ISEC ), LDA, &
                                                 ONE, C( 1, II ), LDC,MB,NB,NBT,KB )
               END IF
  150       CONTINUE
!$omp end parallel do
         ELSE
!*
!*           Form  C := alpha*B*A + beta*C. Right, Lower.
!*
!$omp parallel do default(none) schedule(dynamic) &
!$omp private(IX,II,ISEC,J,UISEC,TMP0,TMP1,TMP2,TMP3) &
!$omp private(I)                                      &
!$omp shared(N,RCB,ALPHA,A,T1,B,C,ONE,LDB,LDC)        &
!$omp shared(MB,NB,NBT,KB)
            DO 180, IX = N, 1, -RCB
               II = MAX( 1, IX-RCB+1 )
               ISEC = IX-II+1
!*
!*              T1 := A, the lower triangular part of a square diagonal
!*              block of A is copied to lower pangular part of T1 and
!*              the transpose of the strictly lower triangular part of
!*              the block of A is copied to the strictly upper
!*              triangular part of T1.
!*
               DO 170, J = II, II+ISEC-1, 2
                  UISEC = II+ISEC-J-2-MOD( II+ISEC-J-2, 2 )
                  T1( J-II+1, J-II+1 ) = ALPHA*A( J, J )
                  IF( MOD( II+ISEC-J-2, 2 ).EQ.0 )THEN
                     TMP0 = A(J+1,J)
                     T1( J-II+2, J-II+1 ) = ALPHA*TMP0
                     T1( J-II+1, J-II+2 ) = ALPHA*TMP0
                     T1( J-II+2, J-II+2 ) = ALPHA*A( J+1, J+1 )
                  ELSE IF( J.LE.II+ISEC-3 )THEN
                     TMP0 = A(J+1,J)
                     TMP1 = A(J+2,J)
                     T1( J-II+2, J-II+1 ) = ALPHA*TMP0
                     T1( J-II+1, J-II+2 ) = ALPHA*TMP0
                     T1( J-II+3, J-II+1 ) = ALPHA*TMP1
                     T1( J-II+1, J-II+3 ) = ALPHA*TMP1
                     T1( J-II+2, J-II+2 ) = ALPHA*A( J+1, J+1 )
                     T1( J-II+3, J-II+2 ) = ALPHA*A( J+2, J+1 )
                     T1( J-II+2, J-II+3 ) = ALPHA*A( J+2, J+1 )
                  END IF
                  DO 160 I = II+ISEC-UISEC, II+ISEC-1, 2
                     TMP0 = A(I,J)
                     TMP1 = A(I+1,J)
                     TMP2 = A(I,J+1)
                     TMP3 = A(I+1,J+1)
                     T1( I-II+1, J-II+1 ) = ALPHA*TMP0
                     T1( I-II+2, J-II+1 ) = ALPHA*TMP1
                     T1( I-II+1, J-II+2 ) = ALPHA*TMP2
                     T1( I-II+2, J-II+2 ) = ALPHA*TMP3
                     T1( J-II+1, I-II+1 ) = ALPHA*TMP0
                     T1( J-II+1, I-II+2 ) = ALPHA*TMP1
                     T1( J-II+2, I-II+1 ) = ALPHA*TMP2
                     T1( J-II+2, I-II+2 ) = ALPHA*TMP3
  160             CONTINUE
  170          CONTINUE
!*
!*              C := T1'*B + beta*C, general matrix multiplication
!*              involving the symmetric diagonal block of A stored
!*              as a full matrix block in T1.
!*
               CALL DGEMM ( 'N', 'N', M, ISEC, ISEC, ONE, &
                                      B( 1, II ), LDB, T1( 1, 1 ), RCB, &
                                                BETA, C( 1, II ), LDC,MB,NB,NBT,KB )
               IF( II.GT.1 )THEN
!*
!*                 C := alpha*B*A' + C, general matrix multiply involving
!!*                 the transpose of a rectangular block of A.

                  CALL DGEMM ( 'N', 'T', M, ISEC, II-1, ALPHA, &
                                       B( 1, 1 ), LDB, A( II, 1 ), LDA, &
                                                 ONE, C( 1, II ), LDC,MB,NB,NBT,KB )
               END IF
               IF( II+ISEC.LE.N )THEN
!*
!*                 C := alpha*B*A + C, general matrix multiply
!*                 involving a rectangular block of A.
!*
                  CALL DGEMM ( 'N', 'N', M, ISEC, N-II-ISEC+1, ALPHA, &
                           B( 1, II+ISEC ), LDB, A( II+ISEC, II ), LDA, &
                                                 ONE, C( 1, II ), LDC, MB,NB,NBT,KB )
               END IF
  180       CONTINUE
!$omp end parallel do
         END IF
      END IF

END SUBROUTINE








SUBROUTINE DSYR2K( UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,  &
           BETA,C,LDC,RCB,MB,NB,NBT,KB) 
       !dir$ attributes code_align : 32 :: DSYR2K
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: DSYR2K
      use omp_lib
      implicit none
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC, RCB
      INTEGER            MB,NB,NBT,KB
      DOUBLE PRECISION   ALPHA, BETA

      !DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A,B,C
#if 0
*     ..
*
*  Purpose
*  =======
*
*  DSYR2K  performs one of the symmetric rank 2k operations
*
*     C := alpha*A*B' + alpha*B*A' + beta*C,
*
*  or
*
*     C := alpha*A'*B + alpha*B'*A + beta*C,
*
*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
*  matrices in the second case.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
*                                        beta*C.
*
*              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
*                                        beta*C.
*
*              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
*                                        beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns  of the  matrices  A and B,  and on  entry  with
*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
*           of rows of the matrices  A and B.  K must be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  k by n  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDB must be at least  max( 1, n ), otherwise  LDB must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  symmetric matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  symmetric matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  -- Rewritten in December-1993
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*  -- Rewritten in Mars-1995.
*     Superscalar GEMM-Based Level 3 BLAS (Version 0.1).
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*  -- Modified in October-1997.
*     Superscalar GEMM-Based Level 3 BLAS (Version 0.1).
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*
*     .. Local Scalars ..
#endif
      INTEGER            INFO, NROWA
      INTEGER            I, II, ISEC, J
      INTEGER            UISEC, RISEC
      LOGICAL            UPPER, NOTR
!*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
!*     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!*     .. External Subroutines ..
!      EXTERNAL           XERBLA
!      EXTERNAL           DGEMM
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      DOUBLE PRECISION TMP0,TMP1,TMP2,TMP3
!*     .. User specified parameters for DSYR2K ..
!      INTEGER            RCB
!      PARAMETER        ( RCB = 96 )
!*     .. Local Arrays ..
      DOUBLE PRECISION   T1( RCB, RCB )

      UPPER = LSAME( UPLO, 'U' )
      NOTR = LSAME( TRANS, 'N' )
      IF( NOTR )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      INFO = 0
      IF( ( .NOT.UPPER ).AND.( .NOT.LSAME( UPLO, 'L' ) ) )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTR ).AND.( .NOT.LSAME( TRANS, 'T' ) ).AND. &
                                     ( .NOT.LSAME( TRANS, 'C' ) ) )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( K.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
          RETURN
      END IF

      IF( ( N.EQ.0 ).OR.
         ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) RETURN
     
!*
!*     And when alpha.eq.zero or k.eq.0.
!*
      IF( ALPHA.EQ.ZERO.OR.K.EQ.0 )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 30, J = 1, N
                  UISEC = J-MOD( J, 4 )
!$OMP SIMD
                  DO 10, I = 1, UISEC, 4
                     C( I, J ) = ZERO
                     C( I+1, J ) = ZERO
                     C( I+2, J ) = ZERO
                     C( I+3, J ) = ZERO
   10             CONTINUE
                  DO 20, I = UISEC+1, J
                     C( I, J ) = ZERO
   20             CONTINUE
   30          CONTINUE
            ELSE
               DO 60, J = 1, N
                  UISEC = J-MOD( J, 4 )
!$OMP SIMD
                  DO 40, I = 1, UISEC, 4
                     C( I, J ) = BETA*C( I, J )
                     C( I+1, J ) = BETA*C( I+1, J )
                     C( I+2, J ) = BETA*C( I+2, J )
                     C( I+3, J ) = BETA*C( I+3, J )
40                   CONTINUE
!$OMP SIMD
                  DO 50, I = UISEC+1, J
                     C( I, J ) = BETA*C( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 100, J = 1, N
                  RISEC = MOD( N-J, 4 )+1
!$OMP SIMD
                  DO 80, I = J, J+RISEC-1
                     C( I, J ) = ZERO
80                   CONTINUE
!$OMP SIMD
                  DO 90, I = J+RISEC, N, 4
                     C( I, J ) = ZERO
                     C( I+1, J ) = ZERO
                     C( I+2, J ) = ZERO
                     C( I+3, J ) = ZERO
   90             CONTINUE
  100          CONTINUE
            ELSE
               DO 130, J = 1, N
                  RISEC = MOD( N-J, 4 )+1
!$OMP SIMD
                  DO 110, I = J, J+RISEC-1
                     C( I, J ) = BETA*C( I, J )
110                  CONTINUE
!$OMP SIMD
                  DO 120, I = J+RISEC, N, 4
                     C( I, J ) = BETA*C( I, J )
                     C( I+1, J ) = BETA*C( I+1, J )
                     C( I+2, J ) = BETA*C( I+2, J )
                     C( I+3, J ) = BETA*C( I+3, J )
  120             CONTINUE
  130          CONTINUE
            END IF
         END IF
         RETURN
      END IF
!*
!*     Start the operations.
      !*
      TMP0 = 0.0D+0
      TMP1 = 0.0D+0
      TMP2 = 0.0D+0
      TMP3 = 0.0D+0
      IF( UPPER )THEN
         IF( NOTR )THEN
!*
!*           Form  C := alpha*A*B' + alpha*B*A' + beta*C. Upper, Notr.
!*
!$omp parallel do default(none) schedule(dynamic) &
!$omp firstprivate(TMP0,TMP1,TMP2,TMP3)           &
!$omp private(II,ISEC,J,UISEC,I)                  &
!$omp shared(N,RCB,K,ALPHA,A,B,T1,C,LDA,LDB,ZERO) &
!$omp shared(MB,NB,NBT,KB,BETA)
            DO 160, II = 1, N, RCB
               ISEC = MIN( RCB, N-II+1 )
!*
!*              T1 := alpha*A*B', general matrix multiplication of
!*              rectangular blocks of A and B. An upper triangular
!*              diagonal block of alpha*A*B' + alpha*B*A' can be
!*              constructed from T1. T1 is square.
!*
               CALL DGEMM ( 'N', 'T', ISEC, ISEC, K, ALPHA, A( II, 1 ), &
                        LDA, B( II, 1 ), LDB, ZERO, T1( 1, 1 ),RCB, MB,NB,NBT,KB )
!*
!*              C := T1 + T1' + beta*C, the upper triangular part of C
!*              is updated with beta, the upper triangular part of T1,
!*              and the transpose of the lower triangular part
!*              of T1.
               !
               
               DO 150, J = II+ISEC-2, II-1, -2
                  UISEC = J-II+1-MOD( J-II+1, 2 )
                  !dir$ assume_aligned C:64
                  !dir$ assume_aligned T1:64
                  !dir$ vector aligned
                  !dir$ vector vectorlength(8)
                  !dir$ vector always
                  DO 140, I = II, II+UISEC-1, 2
                     C( I, J ) = BETA*C( I, J ) +
     $                       T1( I-II+1, J-II+1 ) + T1( J-II+1, I-II+1 )
                     C( I+1, J ) = BETA*C( I+1, J ) +
     $                       T1( I-II+2, J-II+1 ) + T1( J-II+1, I-II+2 )
                     C( I, J+1 ) = BETA*C( I, J+1 ) +
     $                       T1( I-II+1, J-II+2 ) + T1( J-II+2, I-II+1 )
                     C( I+1, J+1 ) = BETA*C( I+1, J+1 ) +
     $                       T1( I-II+2, J-II+2 ) + T1( J-II+2, I-II+2 )
  140             CONTINUE
                  IF( MOD( J-II+1, 2 ).EQ.1 )THEN
                     C( J, J ) = BETA*C( J, J ) +
     $                       T1( J-II+1, J-II+1 ) + T1( J-II+1, J-II+1 )
                     C( J, J+1 ) = BETA*C( J, J+1 ) +
     $                       T1( J-II+1, J-II+2 ) + T1( J-II+2, J-II+1 )
                     C( J+1, J+1 ) = BETA*C( J+1, J+1 ) +
     $                       T1( J-II+2, J-II+2 ) + T1( J-II+2, J-II+2 )
                  ELSE
                     C( J+1, J+1 ) = BETA*C( J+1, J+1 ) +
     $                       T1( J-II+2, J-II+2 ) + T1( J-II+2, J-II+2 )
                  END IF
  150          CONTINUE
*
*              C := alpha*A*B' + beta*C  and  C := alpha*B*A' + C,
*              general matrix multiplication of rectangular blocks of C
*              consisting of ISEC columns stretching from 1 to II-1.
*
               IF( II.GT.1 )THEN
                  CALL DGEMM ( 'N', 'T', II-1, ISEC, K, ALPHA,
     $                            A( 1, 1 ), LDA, B( II, 1 ), LDB, BETA,
     $                                                 C( 1, II ), LDC, MB,NB,NBT,KB )
                  CALL DGEMM ( 'N', 'T', II-1, ISEC, K, ALPHA,
     $                             B( 1, 1 ), LDB, A( II, 1 ), LDA, ONE,
     $                                                 C( 1, II ), LDC, MB,NB,NBT,KB )
               END IF
  160       CONTINUE
!$omp end parallel do
         ELSE
*
*           Form  C := alpha*A'*B + alpha*B'*A + beta*C. Upper, Trans.
*
!$omp parallel do default(none) schedule(dynamic) &
!$omp private(II,ISEC,J,UISEC,I)                  &
!$omp shared(N,RCB,K,ALPHA,A,B,T1,C,LDA,LDB,ZERO) &
!$omp shared(MB,NB,NBT,KB,LDC,BETA)
            DO 190, II = 1, N, RCB
               ISEC = MIN( RCB, N-II+1 )
*
*              T1 := alpha*A'*B, general matrix multiplication of
*              rectangular blocks of A and B. An upper triangular
*              diagonal block of alpha*A*B' + alpha*B*A' can be
*              constructed from T1. T1 is square.
*
               CALL DGEMM ( 'T', 'N', ISEC, ISEC, K, ALPHA, A( 1, II ),
     $                     LDA, B( 1, II ), LDB, ZERO, T1( 1, 1 ), RCB , MB,NB,NBT,KB)
*
*              C := T1 + T1' + beta*C, the upper triangular part of C
*              is updated with beta, the upper triangular part of T1,
*              and the transpose of the lower triangular part
*              of T1.
*
               DO 180, J = II+ISEC-2, II-1, -2
                  UISEC = J-II+1-MOD( J-II+1, 2 )
                  !dir$ assume_aligned C:64
                  !dir$ assume_aligned T1:64
                  !dir$ vector aligned
                  !dir$ vector vectorlength(8)
                  !dir$ vector always
                  DO 170, I = II, II+UISEC-1, 2
                     C( I, J ) = BETA*C( I, J ) +
     $                       T1( I-II+1, J-II+1 ) + T1( J-II+1, I-II+1 )
                     C( I+1, J ) = BETA*C( I+1, J ) +
     $                       T1( I-II+2, J-II+1 ) + T1( J-II+1, I-II+2 )
                     C( I, J+1 ) = BETA*C( I, J+1 ) +
     $                       T1( I-II+1, J-II+2 ) + T1( J-II+2, I-II+1 )
                     C( I+1, J+1 ) = BETA*C( I+1, J+1 ) +
     $                       T1( I-II+2, J-II+2 ) + T1( J-II+2, I-II+2 )
  170             CONTINUE
                  IF( MOD( J-II+1, 2 ).EQ.1 )THEN
                     C( J, J ) = BETA*C( J, J ) +
     $                       T1( J-II+1, J-II+1 ) + T1( J-II+1, J-II+1 )
                     C( J, J+1 ) = BETA*C( J, J+1 ) +
     $                       T1( J-II+1, J-II+2 ) + T1( J-II+2, J-II+1 )
                     C( J+1, J+1 ) = BETA*C( J+1, J+1 ) +
     $                       T1( J-II+2, J-II+2 ) + T1( J-II+2, J-II+2 )
                  ELSE
                     C( J+1, J+1 ) = BETA*C( J+1, J+1 ) +
     $                       T1( J-II+2, J-II+2 ) + T1( J-II+2, J-II+2 )
                  END IF
  180          CONTINUE
*
*              C := alpha*A'*B + beta*C  and  C := alpha*B'*A + C,
*              general matrix multiplication of rectangular blocks of C
*              consisting of ISEC columns stretching from 1 to II-1.
*
               IF( II.GT.1 )THEN
                  CALL DGEMM ( 'T', 'N', II-1, ISEC, K, ALPHA,
     $                            A( 1, 1 ), LDA, B( 1, II ), LDB, BETA,
     $                                                 C( 1, II ), LDC ,MB,NB,NBT,KB )
                  CALL DGEMM ( 'T', 'N', II-1, ISEC, K, ALPHA,
     $                             B( 1, 1 ), LDB, A( 1, II ), LDA, ONE,
     $                                                 C( 1, II ), LDC, MB,NB,NBT,KB )
               END IF
  190       CONTINUE
!$omp end parallel do
         END If
      ELSE
         IF( NOTR )THEN
*
*           Form  C := alpha*A*B' + alpha*B*A' + beta*C. Upper, Notr.
*
!$omp parallel do default(none) schedule(dynamic) &
!$omp private(II,ISEC,J,UISEC,I)                  &
!$omp shared(N,RCB,K,ALPHA,A,B,T1,C,LDA,LDB,ZERO) &
!$omp shared(MB,NB,NBT,KB,LDC,BETA)
            DO 220, II = 1, N, RCB
               ISEC = MIN( RCB, N-II+1 )
*
*              T1 := alpha*A*B', general matrix multiplication of
*              rectangular blocks of A and B. An upper triangular
*              diagonal block of alpha*A*B' + alpha*B*A' can be
*              constructed from T1. T1 is square.
*
               CALL DGEMM ( 'N', 'T', ISEC, ISEC, K, ALPHA, A( II, 1 ),
     $                     LDA, B( II, 1 ), LDB, ZERO, T1( 1, 1 ), RCB )
*
*              C := T1 + T1' + beta*C, the lower triangular part of C
*              is updated with beta, the lower triangular part of T1,
*              and the transpose of the lower triangular part
*              of T1.
*
               DO 210, J = II, II+ISEC-1, 2
                  UISEC = II+ISEC-1-J-MOD( II+ISEC-1-J, 2 )
                  !dir$ assume_aligned C:64
                  !dir$ assume_aligned T1:64
                  !dir$ vector aligned
                  !dir$ vector vectorlength(8)
                  !dir$ vector always
                  IF( MOD( ISEC-UISEC, 2 ).EQ.0 )THEN
                     C( J, J ) = BETA*C( J, J ) +
     $                       T1( J-II+1, J-II+1 ) + T1( J-II+1, J-II+1 )
                     C( J+1, J ) = BETA*C( J+1, J ) +
     $                       T1( J-II+2, J-II+1 ) + T1( J-II+1, J-II+2 )
                     C( J+1, J+1 ) = BETA*C( J+1, J+1 ) +
     $                       T1( J-II+2, J-II+2 ) + T1( J-II+2, J-II+2 )
                  ELSE
                     C( J, J ) = BETA*C( J, J ) +
     $                       T1( J-II+1, J-II+1 ) + T1( J-II+1, J-II+1 )
                  END IF
                  DO 200, I = II+ISEC-UISEC, II+ISEC-1, 2
                     C( I, J ) = BETA*C( I, J ) +
     $                       T1( I-II+1, J-II+1 ) + T1( J-II+1, I-II+1 )
                     C( I+1, J ) = BETA*C( I+1, J ) +
     $                       T1( I-II+2, J-II+1 ) + T1( J-II+1, I-II+2 )
                     C( I, J+1 ) = BETA*C( I, J+1 ) +
     $                       T1( I-II+1, J-II+2 ) + T1( J-II+2, I-II+1 )
                     C( I+1, J+1 ) = BETA*C( I+1, J+1 ) +
     $                       T1( I-II+2, J-II+2 ) + T1( J-II+2, I-II+2 )
  200             CONTINUE
  210          CONTINUE
*
*              C := alpha*A*B' + beta*C  and  C := alpha*B*A' + C,
*              general matrix multiplication of rectangular blocks of C
*              consisting of ISEC columns stretching from 1 to II-1.
*
               IF( II+ISEC-1.LT.N )THEN
                  CALL DGEMM ( 'N', 'T', N-II-ISEC+1, ISEC, K, ALPHA,
     $                      A( II+ISEC, 1 ), LDA, B( II, 1 ), LDB, BETA,
     $                                           C( II+ISEC, II ), LDC, MB,NB,NBT,KB )
                  CALL DGEMM ( 'N', 'T', N-II-ISEC+1, ISEC, K, ALPHA,
     $                       B( II+ISEC, 1 ), LDB, A( II, 1 ), LDA, ONE,
     $                                           C( II+ISEC, II ), LDC, MB,NB,NBT,KB )
               END IF
  220       CONTINUE
!$omp end parallel do
         ELSE
*
*           Form  C := alpha*A'*B + alpha*B'*A + beta*C. Upper, Trans.
*
!$omp parallel do default(none) schedule(dynamic) &
!$omp private(II,ISEC,J,UISEC,I)                  &
!$omp shared(N,RCB,K,ALPHA,A,B,T1,C,LDA,LDB,ZERO) &
!$omp shared(MB,NB,NBT,KB,LDC,BETA)
            DO 250, II = 1, N, RCB
               ISEC = MIN( RCB, N-II+1 )
*
*              T1 := alpha*A'*B, general matrix multiplication of
*              rectangular blocks of A and B. An lower triangular
*              diagonal block of alpha*A*B' + alpha*B*A' can be
*              constructed from T1. T1 is square.
*
               CALL DGEMM ( 'T', 'N', ISEC, ISEC, K, ALPHA, A( 1, II ),
     $                     LDA, B( 1, II ), LDB, ZERO, T1( 1, 1 ), RCB )
*
*              C := T1 + T1' + beta*C, the lower triangular part of C
*              is updated with beta, the lower triangular part of T1,
*              and the transpose of the lower triangular part
*              of T1.
*
               DO 240, J = II, II+ISEC-1, 2
                  UISEC = II+ISEC-1-J-MOD( II+ISEC-1-J, 2 )
                  IF( MOD( ISEC-UISEC, 2 ).EQ.0 )THEN
                     C( J, J ) = BETA*C( J, J ) +
     $                       T1( J-II+1, J-II+1 ) + T1( J-II+1, J-II+1 )
                     C( J+1, J ) = BETA*C( J+1, J ) +
     $                       T1( J-II+2, J-II+1 ) + T1( J-II+1, J-II+2 )
                     C( J+1, J+1 ) = BETA*C( J+1, J+1 ) +
     $                       T1( J-II+2, J-II+2 ) + T1( J-II+2, J-II+2 )
                  ELSE
                     C( J, J ) = BETA*C( J, J ) +
     $                       T1( J-II+1, J-II+1 ) + T1( J-II+1, J-II+1 )
                  END IF
                  DO 230, I = II+ISEC-UISEC, II+ISEC-1, 2
                     C( I, J ) = BETA*C( I, J ) +
     $                       T1( I-II+1, J-II+1 ) + T1( J-II+1, I-II+1 )
                     C( I+1, J ) = BETA*C( I+1, J ) +
     $                       T1( I-II+2, J-II+1 ) + T1( J-II+1, I-II+2 )
                     C( I, J+1 ) = BETA*C( I, J+1 ) +
     $                       T1( I-II+1, J-II+2 ) + T1( J-II+2, I-II+1 )
                     C( I+1, J+1 ) = BETA*C( I+1, J+1 ) +
     $                       T1( I-II+2, J-II+2 ) + T1( J-II+2, I-II+2 )
  230             CONTINUE
  240          CONTINUE
*
*              C := alpha*A'*B + beta*C  and  C := alpha*B'*A + C,
*              general matrix multiplication of rectangular blocks of C
*              consisting of ISEC columns stretching from 1 to II-1.
*
               IF( II+ISEC-1.LT.N )THEN
                  CALL DGEMM ( 'T', 'N', N-II-ISEC+1, ISEC, K, ALPHA,
     $                      A( 1, II+ISEC ), LDA, B( 1, II ), LDB, BETA,
     $                                           C( II+ISEC, II ), LDC, MB,NB,NBT,KB )
                  CALL DGEMM ( 'T', 'N', N-II-ISEC+1, ISEC, K, ALPHA,
     $                       B( 1, II+ISEC ), LDB, A( 1, II ), LDA, ONE,
     $                                           C( II+ISEC, II ), LDC, MB,NB,NBT,KB )
               END IF
  250       CONTINUE
!$omp end parallel do
         END IF
      END IF

   
      END SUBROUTINE

