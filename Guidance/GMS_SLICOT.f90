

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE AB04MD( TYPE, N, M, P, ALPHA, BETA, A, LDA, B, LDB, C, &
     LDC, D, LDD, IWORK, DWORK, LDWORK, INFO ) !GCC$ ATTRIBUTES hot :: AB04MD !GCC$ ATTRIBUTES aligned(32) :: AB04MD !GCC$ ATTRIBUTES no_stack_protector :: AB04MD
#elif defined(__INTEL_COMPILER) || defined(__ICC)
SUBROUTINE AB04MD( TYPE, N, M, P, ALPHA, BETA, A, LDA, B, LDB, C, &
     LDC, D, LDD, IWORK, DWORK, LDWORK, INFO )
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: AB04MD
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: AB04MD
  
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To perform a transformation on the parameters (A,B,C,D) of a
C     system, which is equivalent to a bilinear transformation of the
C     corresponding transfer function matrix.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYPE    CHARACTER*1
C             Indicates the type of the original system and the
C             transformation to be performed as follows:
C             = 'D':  discrete-time   -> continuous-time;
C             = 'C':  continuous-time -> discrete-time.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the state matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     ALPHA,  (input) DOUBLE PRECISION
C     BETA    Parameters specifying the bilinear transformation.
C             Recommended values for stable systems: ALPHA = 1,
C             BETA = 1.  ALPHA <> 0, BETA <> 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state matrix A of the original system.
C             On exit, the leading N-by-N part of this array contains
C                              _
C             the state matrix A of the transformed system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B of the original system.
C             On exit, the leading N-by-M part of this array contains
C                              _
C             the input matrix B of the transformed system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the output matrix C of the original system.
C             On exit, the leading P-by-N part of this array contains
C                               _
C             the output matrix C of the transformed system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the input/output matrix D for the original system.
C             On exit, the leading P-by-M part of this array contains
C                                     _
C             the input/output matrix D of the transformed system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,N).
C             For optimum performance LDWORK >= MAX(1,N*NB), where NB
C             is the optimal blocksize.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the matrix (ALPHA*I + A) is exactly singular;
C             = 2:  if the matrix  (BETA*I - A) is exactly singular.
C
C     METHOD
C
C     The parameters of the discrete-time system are transformed into
C     the parameters of the continuous-time system (TYPE = 'D'), or
C     vice-versa (TYPE = 'C') by the transformation:
C
C     1.  Discrete -> continuous
C         _                     -1
C         A = beta*(alpha*I + A)  * (A - alpha*I)
C         _                                     -1
C         B = sqrt(2*alpha*beta) * (alpha*I + A)  * B
C         _                                         -1
C         C = sqrt(2*alpha*beta) * C * (alpha*I + A)
C         _                        -1
C         D = D - C * (alpha*I + A)  * B
C
C     which is equivalent to the bilinear transformation
C
C                       z - alpha
C         z -> s = beta ---------  .
C                       z + alpha
C
C     of one transfer matrix onto the other.
C
C     2.  Continuous -> discrete
C         _                     -1
C         A = alpha*(beta*I - A)  * (beta*I + A)
C         _                                    -1
C         B = sqrt(2*alpha*beta) * (beta*I - A)  * B
C         _                                        -1
C         C = sqrt(2*alpha*beta) * C * (beta*I - A)
C         _                       -1
C         D = D + C * (beta*I - A)  * B
C
C     which is equivalent to the bilinear transformation
C
C                      beta + s
C       s -> z = alpha -------- .
C                      beta - s
C
C     of one transfer matrix onto the other.
C
C     REFERENCES
C
C     [1] Al-Saggaf, U.M. and Franklin, G.F.
C         Model reduction via balanced realizations: a extension and
C         frequency weighting techniques.
C         IEEE Trans. Autom. Contr., AC-33, pp. 687-692, 1988.
C
C     NUMERICAL ASPECTS
C                                                      3
C     The time taken is approximately proportional to N .
C     The accuracy depends mainly on the condition number of the matrix
C     to be inverted.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, and
C                  A. Varga, German Aerospace Research Establishment,
C                  Oberpfaffenhofen, Germany, Nov. 1996.
C     Supersedes Release 2.0 routine AB04AD by W. van der Linden, and
C     A.J. Geurts, Technische Hogeschool Eindhoven, Holland.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Bilinear transformation, continuous-time system, discrete-time
C     system, state-space model.
C
C     ******************************************************************
C
#endif
      use omp_lib
      implicit none
!C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0 )
!C     .. Scalar Arguments ..
      CHARACTER         TYPE
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDWORK, M, N, P
      DOUBLE PRECISION  ALPHA, BETA
!C     .. Array Arguments ..
   !   INTEGER           IWORK(*)
      !  DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*), DWORK(*)
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: C
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
!C     .. Local Scalars ..
      LOGICAL           LTYPE
      INTEGER           I, IP
      DOUBLE PRECISION  AB2, PALPHA, PBETA, SQRAB2
!C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
!C     .. External Subroutines ..
      EXTERNAL          DGEMM, DGETRF, DGETRS, DGETRI, DLASCL, DSCAL, &
                        DSWAP
!C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SIGN, SQRT
!C     .. Executable Statements ..
!C
      INFO = 0
      LTYPE = LSAME( TYPE, 'D' )
!C
!C     Test the input scalar arguments.
!C
      IF( .NOT.LTYPE .AND. .NOT.LSAME( TYPE, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( ALPHA.EQ.ZERO ) THEN
         INFO = -5
      ELSE IF( BETA.EQ.ZERO ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -14
      ELSE IF( LDWORK.LT.MAX( 1, N ) ) THEN
         INFO = -17
      END IF
!C
      IF ( INFO.NE.0 ) THEN
!C
!C        Error return.
!C
       !  CALL XERBLA( 'AB04MD', -INFO )
         RETURN
      END IF
!C
!C     Quick return if possible.
!C
      IF ( MAX( N, M, P ).EQ.0 ) &
          RETURN
!C
!C     (Note: Comments in the code beginning "Workspace:" describe the
!C     minimal amount of real workspace needed at that point in the
!C     code, as well as the preferred amount for good performance.
!C     NB refers to the optimal block size for the immediately
!C     following subroutine, as returned by ILAENV.)
!C
      IF (LTYPE) THEN
!C
!C        Discrete-time to continuous-time with (ALPHA, BETA).
!C
         PALPHA = ALPHA
         PBETA = BETA
      ELSE
!C
!C        Continuous-time to discrete-time with (ALPHA, BETA) is
!C        equivalent with discrete-time to continuous-time with
!C        (-BETA, -ALPHA), if B and C change the sign.
!C
         PALPHA = -BETA
         PBETA = -ALPHA
      END IF
!C
      AB2 = PALPHA*PBETA*TWO
      SQRAB2 = SIGN( SQRT( ABS( AB2 ) ), PALPHA )
!C                          -1
!C     Compute (alpha*I + A)  .
!C
      DO 10 I = 1, N
         A(I,I)  =  A(I,I) + PALPHA
   10 CONTINUE
!C
      CALL DGETRF( N, N, A, LDA, IWORK, INFO )
!C
      IF (INFO.NE.0) THEN
!C
!C        Error return.
!C
         IF (LTYPE) THEN
            INFO = 1
         ELSE
            INFO = 2
         END IF
         RETURN
      END IF
!C                         -1
!C     Compute  (alpha*I+A)  *B.
!C
      CALL DGETRS( 'No transpose', N, M, A, LDA, IWORK, B, LDB, INFO )
!C                               -1
!C     Compute  D - C*(alpha*I+A)  *B.
!C
      CALL DGEMM( 'No transpose', 'No transpose', P, M, N, -ONE, C, &
                  LDC, B, LDB, ONE, D, LDD )
!C
!C     Scale B by  sqrt(2*alpha*beta).
!C
      CALL DLASCL( 'General', 0, 0, ONE, SQRAB2, N, M, B, LDB, INFO )
!C                                                -1
!C     Compute  sqrt(2*alpha*beta)*C*(alpha*I + A)  .
!C
      CALL DTRSM( 'Right', 'Upper', 'No transpose', 'Non-unit', P, N, &
                 SQRAB2, A, LDA, C, LDC )
!C
      CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', P, N, ONE, &
                  A, LDA, C, LDC )
!C
!C     Apply column interchanges to the solution matrix.
      !C
      !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(IWORK,C) PRIVATE(I,IP) IF(N>=100)
      DO 20 I = N-1, 1, -1
         IP = IWORK(I)
         IF ( IP.NE.I ) &
            CALL DSWAP( P, C(1,I), 1, C(1,IP), 1 )
  20  CONTINUE
!C                               -1
!C     Compute beta*(alpha*I + A)  *(A - alpha*I) as
!C                                        -1
!C     beta*I - 2*alpha*beta*(alpha*I + A)  .
!C
!C     Workspace: need N;  prefer N*NB.
!C
      CALL DGETRI( N, A, LDA, IWORK, DWORK, LDWORK, INFO )
      !C
      !$OMP PARALLEL DO SCHEDULE(STATIC,4) DEFAULT(NONE) SHARED(A) PRIVATE(I) IF(N>=100)
      DO 30 I = 1, N
         CALL DSCAL(N, -AB2, A(1,I), 1)
         A(I,I) = A(I,I) + PBETA
   30 CONTINUE

END SUBROUTINE

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE AB01OD( STAGES, JOBU, JOBV, N, M, A, LDA, B, LDB, U, &
 LDU, V, LDV, NCONT, INDCON, KSTAIR, TOL, IWORK, &
 DWORK, LDWORK, INFO ) !GCC$ ATTRIBUTES hot :: AB01OD !GCC$ ATTRIBUTES aligned(32) :: AB01OD !GCC$ ATTRIBUTES no_stack_protector :: AB01OD
#elif defined(__INTEL_COMPILER) || defined(__ICC)
SUBROUTINE AB01OD( STAGES, JOBU, JOBV, N, M, A, LDA, B, LDB, U, &
 LDU, V, LDV, NCONT, INDCON, KSTAIR, TOL, IWORK, &
 DWORK, LDWORK, INFO )

    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: AB01OD
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: AB01OD
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To reduce the matrices A and B using (and optionally accumulating)
C     state-space and input-space transformations U and V respectively,
C     such that the pair of matrices
C
C        Ac = U' * A * U,    Bc = U' * B * V
C
C     are in upper "staircase" form. Specifically,
C
C             [ Acont     *    ]         [ Bcont ]
C        Ac = [                ],   Bc = [       ],
C             [   0    Auncont ]         [   0   ]
C
C        and
C
C                [ A11 A12  . . .  A1,p-1 A1p ]         [ B1 ]
C                [ A21 A22  . . .  A2,p-1 A2p ]         [ 0  ]
C                [  0  A32  . . .  A3,p-1 A3p ]         [ 0  ]
C        Acont = [  .   .   . . .    .     .  ],   Bc = [ .  ],
C                [  .   .     . .    .     .  ]         [ .  ]
C                [  .   .       .    .     .  ]         [ .  ]
C                [  0   0   . . .  Ap,p-1 App ]         [ 0  ]
C
C     where the blocks  B1, A21, ..., Ap,p-1  have full row ranks and
C     p is the controllability index of the pair.  The size of the
C     block Auncont is equal to the dimension of the uncontrollable
C     subspace of the pair (A, B).  The first stage of the reduction,
C     the "forward" stage, accomplishes the reduction to the orthogonal
C     canonical form (see SLICOT library routine AB01ND). The blocks
C     B1, A21, ..., Ap,p-1 are further reduced in a second, "backward"
C     stage to upper triangular form using RQ factorization. Each of
C     these stages is optional.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     STAGES  CHARACTER*1
C             Specifies the reduction stages to be performed as follows:
C             = 'F':  Perform the forward stage only;
C             = 'B':  Perform the backward stage only;
C             = 'A':  Perform both (all) stages.
C
C     JOBU    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix U the state-space transformations as follows:
C             = 'N':  Do not form U;
C             = 'I':  U is internally initialized to the unit matrix (if
C                     STAGES <> 'B'), or updated (if STAGES = 'B'), and
C                     the orthogonal transformation matrix U is
C                     returned.
C
C     JOBV    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix V the input-space transformations as follows:
C             = 'N':  Do not form V;
C             = 'I':  V is initialized to the unit matrix and the
C                     orthogonal transformation matrix V is returned.
C             JOBV is not referenced if STAGES = 'F'.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e. the order of the
C             matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The actual input dimension.  M >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state transition matrix A to be transformed.
C             If STAGES = 'B', A should be in the orthogonal canonical
C             form, as returned by SLICOT library routine AB01ND.
C             On exit, the leading N-by-N part of this array contains
C             the transformed state transition matrix U' * A * U.
C             The leading NCONT-by-NCONT part contains the upper block
C             Hessenberg state matrix Acont in Ac, given by U' * A * U,
C             of a controllable realization for the original system.
C             The elements below the first block-subdiagonal are set to
C             zero.  If STAGES <> 'F', the subdiagonal blocks of A are
C             triangularized by RQ factorization, and the annihilated
C             elements are explicitly zeroed.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B to be transformed.
C             If STAGES = 'B', B should be in the orthogonal canonical
C             form, as returned by SLICOT library routine AB01ND.
C             On exit with STAGES = 'F', the leading N-by-M part of
C             this array contains the transformed input matrix U' * B,
C             with all elements but the first block set to zero.
C             On exit with STAGES <> 'F', the leading N-by-M part of
C             this array contains the transformed input matrix
C             U' * B * V, with all elements but the first block set to
C             zero and the first block in upper triangular form.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,N)
C             If STAGES <> 'B' or JOBU = 'N', then U need not be set
C             on entry.
C             If STAGES = 'B' and JOBU = 'I', then, on entry, the
C             leading N-by-N part of this array must contain the
C             transformation matrix U that reduced the pair to the
C             orthogonal canonical form.
C             On exit, if JOBU = 'I', the leading N-by-N part of this
C             array contains the transformation matrix U that performed
C             the specified reduction.
C             If JOBU = 'N', the array U is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDU = 1 and
C             declare this array to be U(1,1) in the calling program).
C
C     LDU     INTEGER
C             The leading dimension of array U.
C             If JOBU = 'I', LDU >= MAX(1,N);  if JOBU = 'N', LDU >= 1.
C
C     V       (output) DOUBLE PRECISION array, dimension (LDV,M)
C             If JOBV = 'I', then the leading M-by-M part of this array
C             contains the transformation matrix V.
C             If STAGES = 'F', or JOBV = 'N', the array V is not
C             referenced and can be supplied as a dummy array (i.e. set
C             parameter  LDV = 1 and declare this array to be V(1,1) in
C             the calling program).
C
C     LDV     INTEGER
C             The leading dimension of array V.
C             If STAGES <> 'F' and JOBV = 'I', LDV >= MAX(1,M);
C             if STAGES = 'F' or JOBV = 'N', LDV >= 1.
C
C     NCONT   (input/output) INTEGER
C             The order of the controllable state-space representation.
C             NCONT is input only if STAGES = 'B'.
C
C     INDCON  (input/output) INTEGER
C             The number of stairs in the staircase form (also, the
C             controllability index of the controllable part of the
C             system representation).
C             INDCON is input only if STAGES = 'B'.
C
C     KSTAIR  (input/output) INTEGER array, dimension (N)
C             The leading INDCON elements of this array contain the
C             dimensions of the stairs, or, also, the orders of the
C             diagonal blocks of Acont.
C             KSTAIR is input if STAGES = 'B', and output otherwise.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in rank determination when
C             transforming (A, B). If the user sets TOL > 0, then
C             the given value of TOL is used as a lower bound for the
C             reciprocal condition number (see the description of the
C             argument RCOND in the SLICOT routine MB03OD);  a
C             (sub)matrix whose estimated condition number is less than
C             1/TOL is considered to be of full rank.  If the user sets
C             TOL <= 0, then an implicitly computed, default tolerance,
C             defined by  TOLDEF = N*N*EPS,  is used instead, where EPS
C             is the machine precision (see LAPACK Library routine
C             DLAMCH).
C             TOL is not referenced if STAGES = 'B'.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (M)
C             IWORK is not referenced if STAGES = 'B'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             If STAGES <> 'B', LDWORK >= MAX(1, N + MAX(N,3*M));
C             If STAGES =  'B', LDWORK >= MAX(1, M + MAX(N,M)).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     Staircase reduction of the pencil [B|sI - A] is used. Orthogonal
C     transformations U and V are constructed such that
C
C
C                        |B |sI-A      *  . . .  *      *       |
C                        | 1|    11       .      .      .       |
C                        |  |  A    sI-A    .    .      .       |
C                        |  |   21      22    .  .      .       |
C                        |  |        .     .     *      *       |
C     [U'BV|sI - U'AU] = |0 |     0    .     .                  |
C                        |  |            A     sI-A     *       |
C                        |  |             p,p-1    pp           |
C                        |  |                                   |
C                        |0 |         0          0   sI-A       |
C                        |  |                            p+1,p+1|
C
C
!C     where the i-th diagonal block of U'AU has dimension KSTAIR(i),
C     for i = 1,...,p. The value of p is returned in INDCON. The last
C     block contains the uncontrollable modes of the (A,B)-pair which
C     are also the generalized eigenvalues of the above pencil.
C
C     The complete reduction is performed in two stages. The first,
C     forward stage accomplishes the reduction to the orthogonal
C     canonical form. The second, backward stage consists in further
C     reduction to triangular form by applying left and right orthogonal
C     transformations.
C
C     REFERENCES
C
C     [1] Van Dooren, P.
C         The generalized eigenvalue problem in linear system theory.
C         IEEE Trans. Auto. Contr., AC-26, pp. 111-129, 1981.
C
C     [2] Miminis, G. and Paige, C.
C         An algorithm for pole assignment of time-invariant multi-input
C         linear systems.
C         Proc. 21st IEEE CDC, Orlando, Florida, 1, pp. 62-67, 1982.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O((N + M) x N**2) operations and is
C     backward stable (see [1]).
C
C     FURTHER COMMENTS
C
C     If the system matrices A and B are badly scaled, it would be
C     useful to scale them with SLICOT routine TB01ID, before calling
C     the routine.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
C     Supersedes Release 2.0 routine AB01CD by M. Vanbegin, and
C     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium.
C
C     REVISIONS
C
C     January 14, 1997, February 12, 1998, September 22, 2003.
C
C     KEYWORDS
C
C     Controllability, generalized eigenvalue problem, orthogonal
C     transformation, staircase form.
C
C     ******************************************************************
C
#endif
      implicit none
!C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
!C     .. Scalar Arguments ..
      CHARACTER         JOBU, JOBV, STAGES
      INTEGER           INDCON, INFO, LDA, LDB, LDU, LDV, LDWORK, M, N, &
                        NCONT
      DOUBLE PRECISION  TOL
!C     .. Array Arguments ..
      !INTEGER           IWORK(*), KSTAIR(*)
      !DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), U(LDU,*), V(LDV,*)
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
      INTEGER, DIMENSION(:), ALLOCATABLE :: KSTAIR
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: B
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: DWORK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: U
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: V
!C     .. Local Scalars ..
      LOGICAL           LJOBUI, LJOBVI, LSTAGB, LSTGAB
      INTEGER           I, I0, IBSTEP, ITAU, J0, JINI, JWORK, MCRT, MM, &
                        NCRT, WRKOPT
!C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
!C     .. External Subroutines ..
      EXTERNAL          DGERQF, DLACPY, DLASET, DORGRQ, DORMRQ, &
                        DSWAP
!C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
!C     .. Executable Statements ..
!C
      INFO = 0
      LJOBUI = LSAME( JOBU, 'I' )
!C
      LSTAGB = LSAME( STAGES, 'B' )
      LSTGAB = LSAME( STAGES, 'A' ).OR.LSTAGB
!C
      IF ( LSTGAB ) THEN
         LJOBVI = LSAME( JOBV, 'I' )
      END IF
!C
!C     Test the input scalar arguments.
!C
      IF( .NOT.LSTGAB .AND. .NOT.LSAME( STAGES, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LJOBUI .AND. .NOT.LSAME( JOBU, 'N' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDU.LT.1 .OR. ( LJOBUI .AND. LDU.LT.N ) ) THEN
         INFO = -11
      ELSE IF( .NOT.LSTAGB .AND. LDWORK.LT.MAX( 1, N + MAX( N, 3*M ) ) &
              .OR. LSTAGB .AND. LDWORK.LT.MAX( 1, M + MAX( N, M ) ) )  &
              THEN
         INFO = -20
      ELSE IF( LSTAGB .AND. NCONT.GT.N ) THEN
         INFO = -14
      ELSE IF( LSTAGB .AND. INDCON.GT.N ) THEN
         INFO = -15
      ELSE IF( LSTGAB ) THEN
         IF( .NOT.LJOBVI .AND. .NOT.LSAME( JOBV, 'N' ) ) THEN
            INFO = -3
         ELSE IF( LDV.LT.1 .OR. ( LJOBVI .AND. LDV.LT.M ) ) THEN
            INFO = -13
         END IF
      END IF
!C
      IF ( INFO.NE.0 ) THEN
!C
!!C        Error return.
!C
       !  CALL XERBLA( 'AB01OD', -INFO )
         RETURN
      END IF
!C
!C     Quick return if possible.
!C
      IF ( MIN( N, M ).EQ.0 ) THEN
         NCONT  = 0
         INDCON = 0
         IF( N.GT.0 .AND. LJOBUI ) &
            CALL DLASET( 'F', N, N, ZERO, ONE, U, LDU )
         IF( LSTGAB ) THEN
            IF( M.GT.0 .AND. LJOBVI ) &
              CALL DLASET( 'F', M, M, ZERO, ONE, V, LDV )
         END IF
         DWORK(1) = ONE
         RETURN
      END IF
!C
!C     (Note: Comments in the code beginning "Workspace:" describe the
!C     minimal amount of real workspace needed at that point in the
!C     code, as well as the preferred amount for good performance.
!C     NB refers to the optimal block size for the immediately
!C     following subroutine, as returned by ILAENV.)
!C
      ITAU = 1
      WRKOPT = 1
!C
      IF ( .NOT.LSTAGB ) THEN
!C
!C        Perform the forward stage computations of the staircase
!C        algorithm on B and A: reduce the (A, B) pair to orthogonal
!C        canonical form.
!C
!C        Workspace: N + MAX(N,3*M).
!C
         JWORK = N + 1
         CALL AB01ND( JOBU, N, M, A, LDA, B, LDB, NCONT, INDCON, &
                     KSTAIR, U, LDU, DWORK(ITAU), TOL, IWORK,    &
                     DWORK(JWORK), LDWORK-JWORK+1, INFO )
!C
         WRKOPT = INT( DWORK(JWORK) ) + JWORK - 1
      END IF
!C
!C     Exit if no further reduction to triangularize B1 and subdiagonal
!C     blocks of A is required, or if the order of the controllable part
!C     is 0.
!C
      IF ( .NOT.LSTGAB ) THEN
         DWORK(1) = WRKOPT
         RETURN
      ELSE IF ( NCONT.EQ.0 .OR. INDCON.EQ.0 ) THEN
         IF( LJOBVI ) &
            CALL DLASET( 'F', M, M, ZERO, ONE, V, LDV )
         DWORK(1) = WRKOPT
         RETURN
      END IF
!C
!C     Now perform the backward steps except the last one.
!C
      MCRT = KSTAIR(INDCON)
      I0 = NCONT - MCRT + 1
      JWORK = M + 1
!C
      DO 10 IBSTEP = INDCON, 2, -1
         NCRT = KSTAIR(IBSTEP-1)
         J0 = I0 - NCRT
         MM = MIN( NCRT, MCRT )
!C
!C        Compute the RQ factorization of the current subdiagonal block
!C        of A, Ai,i-1 = R*Q (where i is IBSTEP), of dimension
!C        MCRT-by-NCRT, starting in position (I0,J0).
!C        The matrix Q' should postmultiply U, if required.
!C        Workspace: need   M + MCRT;
!C                   prefer M + MCRT*NB.
!C
         CALL DGERQF( MCRT, NCRT, A(I0,J0), LDA, DWORK(ITAU), &
                      DWORK(JWORK), LDWORK-JWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
!C
!C        Set JINI to the first column number in A where the current
!C        transformation Q is to be applied, taking the block Hessenberg
!C!        form into account.
!C
         IF ( IBSTEP.GT.2 ) THEN
            JINI = J0 - KSTAIR(IBSTEP-2)
         ELSE
            JINI = 1
!C
!C           Premultiply the first block row (B1) of B by Q.
!C           Workspace: need   2*M;
!!C                      prefer M + M*NB.
!C
            CALL DORMRQ( 'Left', 'No transpose', NCRT, M, MM, A(I0,J0), &
                        LDA, DWORK(ITAU), B, LDB, DWORK(JWORK), &
                        LDWORK-JWORK+1, INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
         END IF
!C
!C        Premultiply the appropriate block row of A by Q.
!C        Workspace: need   M + N;
!C                   prefer M + N*NB.
!C
         CALL DORMRQ( 'Left', 'No transpose', NCRT, N-JINI+1, MM,  &
                     A(I0,J0), LDA, DWORK(ITAU), A(J0,JINI), LDA, &
                     DWORK(JWORK), LDWORK-JWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
!C
!C        Postmultiply the appropriate block column of A by Q'.
!C        Workspace: need   M +  I0-1;
!C                   prefer M + (I0-1)*NB.
!C
         CALL DORMRQ( 'Right', 'Transpose', I0-1, NCRT, MM, A(I0,J0),  &
                     LDA, DWORK(ITAU), A(1,J0), LDA, DWORK(JWORK),    &
                     LDWORK-JWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
!C
         IF ( LJOBUI ) THEN
!C
!C           Update U, postmultiplying it by Q'.
!C           Workspace: need   M + N;
!C                      prefer M + N*NB.
!C
            CALL DORMRQ( 'Right', 'Transpose', N, NCRT, MM, A(I0,J0),  &
                        LDA, DWORK(ITAU), U(1,J0), LDU, DWORK(JWORK), &
                        LDWORK-JWORK+1, INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
         END IF
!C
!C        Zero the subdiagonal elements of the current subdiagonal block
!C        of A.
!C
         CALL DLASET( 'F', MCRT, NCRT-MCRT, ZERO, ZERO, A(I0,J0), LDA )
         IF ( I0.LT.N ) &
           CALL DLASET( 'L', MCRT-1, MCRT-1, ZERO, ZERO, &
                        A(I0+1,I0-MCRT), LDA )
!C
         MCRT = NCRT
         I0 = J0
!C
   10 CONTINUE
!C
!C     Now perform the last backward step on B, V = Qb'.
!C
!C!     Compute the RQ factorization of the first block of B, B1 = R*Qb.
!C     Workspace: need   M + MCRT;
!C                prefer M + MCRT*NB.
!C
      CALL DGERQF( MCRT, M, B, LDB, DWORK(ITAU), DWORK(JWORK), &
                   LDWORK-JWORK+1, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
!C
      IF ( LJOBVI ) THEN
!C
!C        Accumulate the input-space transformations V.
!C        Workspace: need 2*M;  prefer M + M*NB.
!C
         CALL DLACPY( 'F', MCRT, M-MCRT, B, LDB, V(M-MCRT+1,1), LDV )
         IF ( MCRT.GT.1 ) &
           CALL DLACPY( 'L', MCRT-1, MCRT-1, B(2,M-MCRT+1), LDB, &
                        V(M-MCRT+2,M-MCRT+1), LDV )
         CALL DORGRQ( M, M, MCRT, V, LDV, DWORK(ITAU), DWORK(JWORK), &
                     LDWORK-JWORK+1, INFO )
!C
         DO 20 I = 2, M
            CALL DSWAP( I-1, V(I,1), LDV, V(1,I), 1 )
   20    CONTINUE
!C
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
      END IF
!C
!C     Zero the subdiagonal elements of the submatrix B1.
!C
      CALL DLASET( 'F', MCRT, M-MCRT, ZERO, ZERO, B, LDB )
      IF ( MCRT.GT.1 ) &
        CALL DLASET( 'L', MCRT-1, MCRT-1, ZERO, ZERO, B(2,M-MCRT+1), &
                     LDB )
!C
!C     Set optimal workspace dimension.
!C
      DWORK(1) = WRKOPT
   
END SUBROUTINE
