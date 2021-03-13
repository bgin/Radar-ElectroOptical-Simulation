
! IMPORTANT NOTICE!!
! Set the USE_OMP macro to 0 when linked against multithreaded version of Intel MKL library
!
#if !defined(GMS_SLICOT_OMP_LOOP_PARALLELIZE)
#define GMS_SLICOT_OMP_LOOP_PARALLELIZE 1
#endif

#if defined __GFORTRAN__ && (!defined(__ICC) || !defined(__INTEL_COMPILER))
 SUBROUTINE AB08MZ(EQUIL, N, M, P, A, LDA, B, LDB, C, LDC, D, LDD, &
      RANK, TOL, IWORK, DWORK, ZWORK, LZWORK, INFO ) !GCC$ ATTRIBUTES aligned(32) :: AB08MZ !GCC$ ATTRIBUTES hot :: AB08MZ !GCC$ ATTRIBUTES no_stack_protector :: AB08MZ
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE AB08MZ(EQUIL, N, M, P, A, LDA, B, LDB, C, LDC, D, LDD, &
      RANK, TOL, IWORK, DWORK, ZWORK, LZWORK, INFO )
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: AB08MZ
!DIR$ OPTIMIZE : 3
!DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: AB08MZ
#endif
      implicit none
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To compute the normal rank of the transfer-function matrix of a
C     state-space model (A,B,C,D).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to balance the compound
C             matrix (see METHOD) as follows:
C             = 'S':  Perform balancing (scaling);
C             = 'N':  Do not perform balancing.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of state variables, i.e., the order of the
C             matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     A       (input) COMPLEX*16 array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state dynamics matrix A of the system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) COMPLEX*16 array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input/state matrix B of the system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) COMPLEX*16 array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             state/output matrix C of the system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input) COMPLEX*16 array, dimension (LDD,M)
C             The leading P-by-M part of this array must contain the
C             direct transmission matrix D of the system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     RANK    (output) INTEGER
C             The normal rank of the transfer-function matrix.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             A tolerance used in rank decisions to determine the
C             effective rank, which is defined as the order of the
C             largest leading (or trailing) triangular submatrix in the
C             QR (or RQ) factorization with column (or row) pivoting
C             whose estimated condition number is less than 1/TOL.
C             If the user sets TOL to be less than SQRT((N+P)*(N+M))*EPS
C             then the tolerance is taken as SQRT((N+P)*(N+M))*EPS,
C             where EPS is the machine precision (see LAPACK Library
C             Routine DLAMCH).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (2*N+MAX(M,P)+1)
C
C     DWORK   DOUBLE PRECISION array, dimension (2*MAX(M,P))
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C             On exit, if INFO = 0, ZWORK(1) returns the optimal value
C             of LZWORK.
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.
C             LZWORK >= (N+P)*(N+M) + MAX(MIN(P,M) + MAX(3*M-1,N), 1,
C                                         MIN(P,N) + MAX(3*P-1,N+P,N+M))
C             For optimum performance LZWORK should be larger.
C
C             If LZWORK = -1, then a workspace query is assumed;
C             the routine only calculates the optimal size of the
C             ZWORK array, returns this value as the first entry of
C             the ZWORK array, and no error message related to LZWORK
C             is issued by XERBLA.
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
C     The routine reduces the (N+P)-by-(M+N) compound matrix (B  A)
C                                                            (D  C)
C
C     to one with the same invariant zeros and with D of full row rank.
C     The normal rank of the transfer-function matrix is the rank of D.
C
C     REFERENCES
C
C     [1] Svaricek, F.
C         Computation of the Structural Invariants of Linear
C         Multivariable Systems with an Extended Version of
C         the Program ZEROS.
C         System & Control Letters, 6, pp. 261-266, 1985.
C
C     [2] Emami-Naeini, A. and Van Dooren, P.
C         Computation of Zeros of Linear Multivariable Systems.
C         Automatica, 18, pp. 415-430, 1982.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable (see [2] and [1]).
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001.
C     Complex version: V. Sima, Research Institute for Informatics,
C     Bucharest, Dec. 2008.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2009,
C     Apr. 2009.
C
C     KEYWORDS
C
C     Multivariable system, unitary transformation,
C     structural invariant.
C
C     ******************************************************************
C
#endif
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         EQUIL
      INTEGER           INFO, LDA, LDB, LDC, LDD, LZWORK, M, N, P, RANK
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IWORK(*)
#if defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER           IWORK(*)
      !DIR$ ASSUME_ALIGNED IWORK:64
#endif
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*), ZWORK(*)
      DOUBLE PRECISION  DWORK(*)
#if defined(__ICC) || defined(__INTEL_COMPILER)
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*), ZWORK(*)
!DIR$ ASSUME_ALIGNED A:64
!DIR$ ASSUME_ALIGNED B:64
!DIR$ ASSUME_ALIGNED C:64
!DIR$ ASSUME_ALIGNED D:64
!DIR$ ASSUME_ALIGNED ZWORK:64
      DOUBLE PRECISION  DWORK(*)
!DIR$ ASSUME_ALIGNED DWORK:64
#endif
!C     .. Local Scalars ..
      LOGICAL           LEQUIL, LQUERY
      INTEGER           I, KW, MU, NINFZ, NKROL, NM, NP, NU, RO, &
                        SIGMA, WRKOPT
      DOUBLE PRECISION  MAXRED, SVLMAX, THRESH, TOLER
!C     .. External Functions ..
     
      DOUBLE PRECISION   ZLANGE
      EXTERNAL           ZLANGE
!C     .. External Subroutines ..
      EXTERNAL           ZLACPY
!C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN, SQRT
!C     .. Executable Statements ..
!C
      NP = N + P
      NM = N + M
      INFO = 0
      LEQUIL = LSAME( EQUIL, 'S' )
      LQUERY = ( LZWORK.EQ.-1 )
      WRKOPT = NP*NM
C
C     Test the input scalar arguments.
C
      IF( .NOT.LEQUIL .AND. .NOT.LSAME( EQUIL, 'N' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE
         KW = WRKOPT + MAX( MIN( P, M ) + MAX( 3*M-1, N ), 1, &
                           MIN( P, N ) + MAX( 3*P-1, NP, NM ) )
         IF( LQUERY ) THEN
            SVLMAX = ZERO
            NINFZ  = 0
            CALL AB8NXZ( N, M, P, P, 0, SVLMAX, ZWORK, MAX( 1, NP ), &
                        NINFZ, IWORK, IWORK, MU, NU, NKROL, TOL, IWORK, &
                        DWORK, ZWORK, -1, INFO )
            WRKOPT = MAX( KW, WRKOPT + INT( ZWORK(1) ) )
         ELSE IF( LZWORK.LT.KW ) THEN
            INFO = -17
         END IF
      END IF
!C
      IF ( INFO.NE.0 ) THEN
!C
!C        Error return.
!C
        ! CALL XERBLA( 'AB08MZ', -INFO ) // removing call to xerbla (no need to issue stop)
         RETURN
      ELSE IF( LQUERY ) THEN
         ZWORK(1) = WRKOPT
         RETURN
      END IF
!C
!C     Quick return if possible.
!C
      IF ( MIN( M, P ).EQ.0 ) THEN
         RANK = 0
         ZWORK(1) = ONE
         RETURN
      END IF
!C
      DO 10 I = 1, 2*N+1
         IWORK(I) = 0
 10   CONTINUE
#if 0
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of workspace needed at that point in the code,
C     as well as the preferred amount for good performance.)
C
C     Construct the compound matrix  ( B  A ), dimension (N+P)-by-(M+N).
C                                    ( D  C )
C     Complex workspace: need   (N+P)*(N+M).
C
#endif
      CALL ZLACPY( 'Full', N, M, B, LDB, ZWORK, NP )
      CALL ZLACPY( 'Full', P, M, D, LDD, ZWORK(N+1), NP )
      CALL ZLACPY( 'Full', N, N, A, LDA, ZWORK(NP*M+1), NP )
      CALL ZLACPY( 'Full', P, N, C, LDC, ZWORK(NP*M+N+1), NP )
!C
!C     If required, balance the compound matrix (default MAXRED).
!C     Real Workspace: need   N.
!C
      KW = WRKOPT + 1
      IF ( LEQUIL ) THEN
         MAXRED = ZERO
         CALL TB01IZ( 'A', N, M, P, MAXRED, ZWORK(NP*M+1), NP, ZWORK, &
                      NP, ZWORK(NP*M+N+1), NP, DWORK, INFO )
      END IF
!C
!C     If required, set tolerance.
!C
      THRESH = SQRT( DBLE( NP*NM ) )*DLAMCH( 'Precision' )
      TOLER = TOL
      IF ( TOLER.LT.THRESH ) TOLER = THRESH
      SVLMAX = ZLANGE( 'Frobenius', NP, NM, ZWORK, NP, DWORK )
!C
!C     Reduce this system to one with the same invariant zeros and with
!C     D full row rank MU (the normal rank of the original system).
!C     Real workspace:    need   2*MAX(M,P);
!C     Complex workspace: need   (N+P)*(N+M) +
!C                               MAX( 1, MIN(P,M) + MAX(3*M-1,N),
!C                                       MIN(P,N) + MAX(3*P-1,N+P,N+M) );
!1    C                        prefer larger.
!C     Integer workspace: 2*N+MAX(M,P)+1.
!C
      RO = P
      SIGMA = 0
      NINFZ = 0
      CALL AB8NXZ( N, M, P, RO, SIGMA, SVLMAX, ZWORK, NP, NINFZ, IWORK, &
                  IWORK(N+1), MU, NU, NKROL, TOLER, IWORK(2*N+2),       &
                  DWORK, ZWORK(KW), LZWORK-KW+1, INFO )
      RANK = MU

      ZWORK(1) = MAX( WRKOPT, INT( ZWORK(KW) ) + KW - 1 )


END SUBROUTINE  AB08MZ

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE AB8NXZ( N, M, P, RO, SIGMA, SVLMAX, ABCD, LDABCD, &
                   NINFZ, INFZ, KRONL, MU, NU, NKROL, TOL, IWORK, &
                   DWORK, ZWORK, LZWORK, INFO ) !GCC$ ATTRIBUTES hot :: AB8NXZ !GCC$ ATTRIBUTES aligned(32) :: AB8NXZ !GCC$ ATTRIBUTES no_stack_protector :: AB8NXZ
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE AB8NXZ( N, M, P, RO, SIGMA, SVLMAX, ABCD, LDABCD, &
                   NINFZ, INFZ, KRONL, MU, NU, NKROL, TOL, IWORK, &
                   DWORK, ZWORK, LZWORK, INFO )
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: AB8NXZ
  !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: AB8NZX
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To extract from the (N+P)-by-(M+N) system
C                  ( B  A )
C                  ( D  C )
C     an (NU+MU)-by-(M+NU) "reduced" system
C                  ( B' A')
C                  ( D' C')
C     having the same transmission zeros but with D' of full row rank.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of state variables.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     RO      (input/output) INTEGER
C             On entry,
C             = P     for the original system;
C             = MAX(P-M, 0) for the pertransposed system.
C             On exit, RO contains the last computed rank.
C
C     SIGMA   (input/output) INTEGER
C             On entry,
C             = 0  for the original system;
C             = M  for the pertransposed system.
C             On exit, SIGMA contains the last computed value sigma in
C             the algorithm.
C
C     SVLMAX  (input) DOUBLE PRECISION
C             During each reduction step, the rank-revealing QR
C             factorization of a matrix stops when the estimated minimum
C             singular value is smaller than TOL * MAX(SVLMAX,EMSV),
C             where EMSV is the estimated maximum singular value.
C             SVLMAX >= 0.
C
C     ABCD    (input/output) COMPLEX*16 array, dimension (LDABCD,M+N)
C             On entry, the leading (N+P)-by-(M+N) part of this array
C             must contain the compound input matrix of the system.
C             On exit, the leading (NU+MU)-by-(M+NU) part of this array
C             contains the reduced compound input matrix of the system.
C
C     LDABCD  INTEGER
C             The leading dimension of array ABCD.
C             LDABCD >= MAX(1,N+P).
C
C     NINFZ   (input/output) INTEGER
C             On entry, the currently computed number of infinite zeros.
C             It should be initialized to zero on the first call.
C             NINFZ >= 0.
C             On exit, the number of infinite zeros.
C
C     INFZ    (input/output) INTEGER array, dimension (N)
C             On entry, INFZ(i) must contain the current number of
C             infinite zeros of degree i, where i = 1,2,...,N, found in
C             the previous call(s) of the routine. It should be
C             initialized to zero on the first call.
C             On exit, INFZ(i) contains the number of infinite zeros of
C             degree i, where i = 1,2,...,N.
C
C     KRONL   (input/output) INTEGER array, dimension (N+1)
C             On entry, this array must contain the currently computed
C             left Kronecker (row) indices found in the previous call(s)
C             of the routine. It should be initialized to zero on the
C             first call.
C             On exit, the leading NKROL elements of this array contain
C             the left Kronecker (row) indices.
C
C     MU      (output) INTEGER
C             The normal rank of the transfer function matrix of the
C             original system.
C
C     NU      (output) INTEGER
C             The dimension of the reduced system matrix and the number
C             of (finite) invariant zeros if D' is invertible.
C
C     NKROL   (output) INTEGER
C             The number of left Kronecker indices.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             A tolerance used in rank decisions to determine the
C             effective rank, which is defined as the order of the
C             largest leading (or trailing) triangular submatrix in the
C             QR (or RQ) factorization with column (or row) pivoting
C             whose estimated condition number is less than 1/TOL.
C             NOTE that when SVLMAX > 0, the estimated ranks could be
C             less than those defined above (see SVLMAX).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (MAX(M,P))
C
C     DWORK   DOUBLE PRECISION array, dimension (2*MAX(M,P))
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C             On exit, if INFO = 0, ZWORK(1) returns the optimal value
C             of LZWORK.
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.
C             LZWORK >= MAX( 1, MIN(P,M) + MAX(3*M-1,N),
C                               MIN(P,N) + MAX(3*P-1,N+P,N+M) ).
C             For optimum performance LZWORK should be larger.
C
C             If LZWORK = -1, then a workspace query is assumed;
C             the routine only calculates the optimal size of the
C             ZWORK array, returns this value as the first entry of
C             the ZWORK array, and no error message related to LZWORK
C             is issued by XERBLA.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     REFERENCES
C
C     [1] Svaricek, F.
C         Computation of the Structural Invariants of Linear
C         Multivariable Systems with an Extended Version of
C         the Program ZEROS.
C         System & Control Letters, 6, pp. 261-266, 1985.
C
C     [2] Emami-Naeini, A. and Van Dooren, P.
C         Computation of Zeros of Linear Multivariable Systems.
C         Automatica, 18, pp. 415-430, 1982.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
C     Complex version: V. Sima, Research Institute for Informatics,
C     Bucharest, Nov. 2008 with suggestions from P. Gahinet,
C     The MathWorks.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009,
C     Apr. 2011.
C
C     KEYWORDS
C
C     Generalized eigenvalue problem, Kronecker indices, multivariable
C     system, unitary transformation, structural invariant.
C
C     ******************************************************************
C
#endif
#if(GMS_SLICOT_OMP_LOOP_PARALLELIZE) == 1
       use omp_lib
#endif
       implicit none
!C     .. Parameters ..
      COMPLEX*16        ZERO
      PARAMETER         ( ZERO = ( 0.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION  DZERO
      PARAMETER         ( DZERO = 0.0D0 )
!C     .. Scalar Arguments ..
      INTEGER           INFO, LDABCD, LZWORK, M, MU, N, NINFZ, NKROL, &
                        NU, P, RO, SIGMA
      DOUBLE PRECISION  SVLMAX, TOL
      !C     .. Array Arguments ..
 
      INTEGER           INFZ(*), IWORK(*), KRONL(*)
#if defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER           INFZ(*), IWORK(*), KRONL(*)
      !DIR$ ASSUME_ALIGNED INFZ:64
      !DIR$ ASSUME_ALIGNED IWORK:64
      !DIR$ ASSUME_ALIGNED KRONL:64
#endif
      COMPLEX*16        ABCD(LDABCD,*), ZWORK(*)
#if defined(__ICC) || defined(__INTEL_COMPILER)
      COMPLEX*16        ABCD(LDABCD,*), ZWORK(*)
      !DIR$ ASSUME_ALIGNED ABCD:64
      !DIR$ ASSUME_ALIGNED ZWORK:64
#endif
      DOUBLE PRECISION  DWORK(*)
#if defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION  DWORK(*)
      !DIR$ ASSUME_ALIGNED DWORK:64
#endif
!C     .. Local Scalars ..
      LOGICAL           LQUERY
      INTEGER           I1, IK, IROW, ITAU, IZ, JWORK, MM1, MNTAU, MNU, &
                        MPM, NP, RANK, RO1, TAU, WRKOPT
      COMPLEX*16        TC
!C     .. Local Arrays ..
      DOUBLE PRECISION  SVAL(3)
!C     .. External Subroutines ..
      EXTERNAL          ZLAPMT, ZLARFG, ZLASET,ZLATZM, ZUNMQR, ZUNMRQ
                     
!C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, INT, MAX, MIN
!C     .. Executable Statements ..
!C
      NP   = N + P
      MPM  = MIN( P, M )
      INFO = 0
      LQUERY = ( LZWORK.EQ.-1 )
!C
!C     Test the input scalar arguments.
!C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( RO.NE.P .AND. RO.NE.MAX( P-M, 0 ) ) THEN
         INFO = -4
      ELSE IF( SIGMA.NE.0 .AND. SIGMA.NE.M ) THEN
         INFO = -5
      ELSE IF( SVLMAX.LT.DZERO ) THEN
         INFO = -6
      ELSE IF( LDABCD.LT.MAX( 1, NP ) ) THEN
         INFO = -8
      ELSE IF( NINFZ.LT.0 ) THEN
         INFO = -9
      ELSE
         JWORK = MAX( 1,      MPM + MAX( 3*M - 1, N ), &
                     MIN( P, N ) + MAX( 3*P - 1, NP, N+M ) )
         IF( LQUERY ) THEN
            IF( M.GT.0 ) THEN
               CALL ZUNMQR( 'Left', 'Conjugate', P, N, MPM, ABCD, &
                           LDABCD, ZWORK, ABCD, LDABCD, ZWORK, -1, &
                           INFO )
               WRKOPT = MAX( JWORK, MPM + INT( ZWORK(1) ) )
            ELSE
               WRKOPT = JWORK
            END IF
            CALL ZUNMRQ( 'Right', 'ConjTranspose', NP, N, MIN( P, N ), &
                        ABCD, LDABCD, ZWORK, ABCD, LDABCD, ZWORK, -1,  &
                        INFO )
            WRKOPT = MAX( WRKOPT, MIN( P, N ) + INT( ZWORK(1) ) )
            CALL ZUNMRQ( 'Left', 'NoTranspose', N, M+N, MIN( P, N ),  &
                        ABCD, LDABCD, ZWORK, ABCD, LDABCD, ZWORK, -1, &
                        INFO )
            WRKOPT = MAX( WRKOPT, MIN( P, N ) + INT( ZWORK(1) ) )
         ELSE IF( LZWORK.LT.JWORK ) THEN
            INFO = -19
         END IF
      END IF
!C
      IF ( INFO.NE.0 ) THEN
!C
!C        Error return.
!C
          RETURN
      ELSE IF( LQUERY ) THEN
         ZWORK(1) = WRKOPT
         RETURN
      END IF
!C
      MU = P
      NU = N
!C
      IZ = 0
      IK = 1
      MM1 = M + 1
      ITAU = 1
      NKROL = 0
      WRKOPT = 1
!C
!C     Main reduction loop:
!C
!C            M   NU                  M     NU
!C      NU  [ B   A ]           NU  [ B     A ]
!C      MU  [ D   C ]  -->    SIGMA [ RD   C1 ]   (SIGMA = rank(D) =
!C                             TAU  [ 0    C2 ]    row size of RD)
!C
!C                                    M   NU-RO  RO
!C                            NU-RO [ B1   A11  A12 ]
!C                     -->      RO  [ B2   A21  A22 ]  (RO = rank(C2) =
!C                            SIGMA [ RD   C11  C12 ]   col size of LC)
!C                             TAU  [ 0     0   LC  ]
!C
!C                                     M   NU-RO
!C                            NU-RO [ B1   A11 ]     NU := NU - RO
!C                                  [----------]     MU := RO + SIGMA
!C                     -->      RO  [ B2   A21 ]      D := [B2;RD]
!C                            SIGMA [ RD   C11 ]      C := [A21;C11]
!C
   20 IF ( MU.EQ.0 ) &
          GO TO 80
!C
!C     (Note: Comments in the code beginning "xWorkspace:", where x is
!C     I, D, or C, describe the minimal amount of integer, real and
!C     complex workspace needed at that point in the code, respectively,
!C     as well as the preferred amount for good performance.)
!C
      RO1 = RO
      MNU = M + NU
      IF ( M.GT.0 ) THEN
         IF ( SIGMA.NE.0 ) THEN
            IROW = NU + 1
!C
!C           Compress rows of D.  First exploit triangular shape.
!C           CWorkspace: need   M+N-1.
            !C
#if(GMS_SLICOT_OMP_LOOP_PARALLELIZE) == 1
!$OMP PARALLEL DO PRIVATE(I1,IROW) SCHEDULE(STATIC)            
            DO 40 I1 = 1, SIGMA
               CALL ZLARFG( RO+1, ABCD(IROW,I1), ABCD(IROW+1,I1), 1, &
                           TC )
               CALL ZLATZM( 'L', RO+1, MNU-I1, ABCD(IROW+1,I1), 1, &
                           DCONJG( TC ), ABCD(IROW,I1+1),          &
                           ABCD(IROW+1,I1+1), LDABCD, ZWORK )
               IROW = IROW + 1
40             CONTINUE
!$OMP END PARALLEL DO
#endif
            CALL ZLASET( 'Lower', RO+SIGMA-1, SIGMA, ZERO, ZERO, &
                        ABCD(NU+2,1), LDABCD )
         END IF
!C
!C        Continue with Householder with column pivoting.
!C
!C        The rank of D is the number of (estimated) singular values
!C        that are greater than TOL * MAX(SVLMAX,EMSV). This number
!C        includes the singular values of the first SIGMA columns.
!C        IWorkspace: need   M;
!C        RWorkspace: need   2*M;
!C        CWorkspace: need   min(RO1,M) + 3*M - 1.  RO1 <= P.
!C
         IF ( SIGMA.LT.M ) THEN
            JWORK = ITAU + MIN( RO1, M )
            I1    = SIGMA + 1
            IROW  = NU + I1
            CALL MB3OYZ( RO1, M-SIGMA, ABCD(IROW,I1), LDABCD, TOL,     &
                        SVLMAX, RANK, SVAL, IWORK, ZWORK(ITAU), DWORK, &
                        ZWORK(JWORK), INFO )
            WRKOPT = MAX( WRKOPT, JWORK + 3*M - 2 )
!C
!C           Apply the column permutations to matrices B and part of D.
!C
            CALL ZLAPMT( .TRUE., NU+SIGMA, M-SIGMA, ABCD(1,I1), LDABCD, &
                        IWORK )
!C
            IF ( RANK.GT.0 ) THEN
!C
!C              Apply the Householder transformations to the submatrix C.
!C              CWorkspace:    need   min(RO1,M) + NU;
!C                             prefer min(RO1,M) + NU*NB.
!C
               CALL ZUNMQR( 'Left', 'Conjugate', RO1, NU, RANK,  &
                           ABCD(IROW,I1), LDABCD, ZWORK(ITAU),   &
                           ABCD(IROW,MM1), LDABCD, ZWORK(JWORK), &
                           LZWORK-JWORK+1, INFO )
               WRKOPT = MAX( WRKOPT, INT( ZWORK(JWORK) ) + JWORK - 1 )
               IF ( RO1.GT.1 ) &
                 CALL ZLASET( 'Lower', RO1-1, MIN( RO1-1, RANK ), ZERO, &
                              ZERO, ABCD(IROW+1,I1), LDABCD )
               RO1 = RO1 - RANK
            END IF
         END IF
      END IF

      TAU = RO1
      SIGMA = MU - TAU
!C
!C     Determination of the orders of the infinite zeros.
!C
      IF ( IZ.GT.0 ) THEN
         INFZ(IZ) = INFZ(IZ) + RO - TAU
         NINFZ = NINFZ + IZ*( RO - TAU )
      END IF
      IF ( RO1.EQ.0 ) &
               GO TO 80
      IZ = IZ + 1

      IF ( NU.LE.0 ) THEN
         MU = SIGMA
         NU = 0
         RO = 0
      ELSE
!C
!C        Compress the columns of C2 using RQ factorization with row
!C        pivoting, P * C2 = R * Q.
!C
         I1 = NU + SIGMA + 1
         MNTAU = MIN( TAU, NU )
         JWORK = ITAU + MNTAU
!C
!C        The rank of C2 is the number of (estimated) singular values
!C        greater than TOL * MAX(SVLMAX,EMSV).
!C        IWorkspace: need TAU;
!C        RWorkspace: need 2*TAU;
!C        CWorkspace: need min(TAU,NU) + 3*TAU - 1.
!C
         CALL MB3PYZ( TAU, NU, ABCD(I1,MM1), LDABCD, TOL, SVLMAX, RANK,  &
                     SVAL, IWORK, ZWORK(ITAU), DWORK, ZWORK(JWORK),     &
                     INFO )
         WRKOPT = MAX( WRKOPT, JWORK + 3*TAU - 1 )
         IF ( RANK.GT.0 ) THEN
            IROW = I1 + TAU - RANK
!C
!C           Apply Q' to the first NU columns of [A; C1] from the right.
!C           CWorkspace: need   min(TAU,NU) + NU + SIGMA; SIGMA <= P;
!C                       prefer min(TAU,NU) + (NU  + SIGMA)*NB.
!C
            CALL ZUNMRQ( 'Right', 'ConjTranspose', I1-1, NU, RANK,     &
                        ABCD(IROW,MM1), LDABCD, ZWORK(MNTAU-RANK+1),   &
                        ABCD(1,MM1), LDABCD, ZWORK(JWORK),             &
                        LZWORK-JWORK+1, INFO )
            WRKOPT = MAX( WRKOPT, INT( ZWORK(JWORK) ) + JWORK - 1 )
!C
!C           Apply Q to the first NU rows and M + NU columns of [ B  A ]
!C           from the left.
!C           CWorkspace: need   min(TAU,NU) + M + NU;
!C                       prefer min(TAU,NU) + (M + NU)*NB.
!C
            CALL ZUNMRQ( 'Left', 'NoTranspose', NU, MNU, RANK,       &
                        ABCD(IROW,MM1), LDABCD, ZWORK(MNTAU-RANK+1), &
                        ABCD, LDABCD, ZWORK(JWORK), LZWORK-JWORK+1,  &
                        INFO )
            WRKOPT = MAX( WRKOPT, INT( ZWORK(JWORK) ) + JWORK - 1 )
!C
            CALL ZLASET( 'Full', RANK, NU-RANK, ZERO, ZERO,  &
                        ABCD(IROW,MM1), LDABCD )
            IF ( RANK.GT.1 )  &
              CALL ZLASET( 'Lower', RANK-1, RANK-1, ZERO, ZERO, &
                          ABCD(IROW+1,MM1+NU-RANK), LDABCD )
         END IF

         RO = RANK
      END IF
!C
!C     Determine the left Kronecker indices (row indices).
!C
      KRONL(IK) = KRONL(IK) + TAU - RO
      NKROL = NKROL + KRONL(IK)
      IK = IK + 1
!C
!C     C and D are updated to [A21 ; C11] and [B2 ; RD].
!C
      NU = NU - RO
      MU = SIGMA + RO
      IF ( RO.NE.0 ) &
                GO TO 20
!C
   80 CONTINUE
      ZWORK(1) = WRKOPT
    
END SUBROUTINE  AB8NXZ

    

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE MB3OYZ( M, N, A, LDA, RCOND, SVLMAX, RANK, SVAL, JPVT, &
     TAU, DWORK, ZWORK, INFO ) !GCC$ ATTRIBUTES hot :: MB3OYZ !GCC$ ATTRIBUTES aligned(32) :: MB3OYZ !GCC$ ATTRIBUTES no_stack_protector :: MB3OYZ
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE MB3OYZ( M, N, A, LDA, RCOND, SVLMAX, RANK, SVAL, JPVT, &
     TAU, DWORK, ZWORK, INFO )
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: MB3OYZ
  !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: MB3OYZ
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To compute a rank-revealing QR factorization of a complex general
C     M-by-N matrix  A,  which may be rank-deficient, and estimate its
C     effective rank using incremental condition estimation.
C
C     The routine uses a truncated QR factorization with column pivoting
C                                   [ R11 R12 ]
C        A * P = Q * R,  where  R = [         ],
C                                   [  0  R22 ]
C     with R11 defined as the largest leading upper triangular submatrix
C     whose estimated condition number is less than 1/RCOND.  The order
C     of R11, RANK, is the effective rank of A.  Condition estimation is
C     performed during the QR factorization process.  Matrix R22 is full
C     (but of small norm), or empty.
C
C     MB3OYZ  does not perform any scaling of the matrix A.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix A.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix A.  N >= 0.
C
C     A       (input/output) COMPLEX*16 array, dimension ( LDA, N )
C             On entry, the leading M-by-N part of this array must
C             contain the given matrix A.
C             On exit, the leading RANK-by-RANK upper triangular part
C             of A contains the triangular factor R11, and the elements
C             below the diagonal in the first  RANK  columns, with the
C             array TAU, represent the unitary matrix Q as a product
C             of  RANK  elementary reflectors.
C             The remaining  N-RANK  columns contain the result of the
C             QR factorization process used.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,M).
C
C     RCOND   (input) DOUBLE PRECISION
C             RCOND is used to determine the effective rank of A, which
C             is defined as the order of the largest leading triangular
C             submatrix R11 in the QR factorization with pivoting of A,
C             whose estimated condition number is less than 1/RCOND.
C             0 <= RCOND <= 1.
C             NOTE that when SVLMAX > 0, the estimated rank could be
C             less than that defined above (see SVLMAX).
C
C     SVLMAX  (input) DOUBLE PRECISION
C             If A is a submatrix of another matrix B, and the rank
C             decision should be related to that matrix, then SVLMAX
C             should be an estimate of the largest singular value of B
C             (for instance, the Frobenius norm of B).  If this is not
C             the case, the input value SVLMAX = 0 should work.
C             SVLMAX >= 0.
C
C     RANK    (output) INTEGER
C             The effective (estimated) rank of A, i.e., the order of
C             the submatrix R11.
C
C     SVAL    (output) DOUBLE PRECISION array, dimension ( 3 )
C             The estimates of some of the singular values of the
C             triangular factor R:
C             SVAL(1): largest singular value of R(1:RANK,1:RANK);
C             SVAL(2): smallest singular value of R(1:RANK,1:RANK);
C             SVAL(3): smallest singular value of R(1:RANK+1,1:RANK+1),
C                      if RANK < MIN( M, N ), or of R(1:RANK,1:RANK),
C                      otherwise.
C             If the triangular factorization is a rank-revealing one
C             (which will be the case if the leading columns were well-
C             conditioned), then SVAL(1) will also be an estimate for
C             the largest singular value of A, and SVAL(2) and SVAL(3)
C             will be estimates for the RANK-th and (RANK+1)-st singular
C             values of A, respectively.
C             By examining these values, one can confirm that the rank
C             is well defined with respect to the chosen value of RCOND.
C             The ratio SVAL(1)/SVAL(2) is an estimate of the condition
C             number of R(1:RANK,1:RANK).
C
C     JPVT    (output) INTEGER array, dimension ( N )
C             If JPVT(i) = k, then the i-th column of A*P was the k-th
C             column of A.
C
C     TAU     (output) COMPLEX*16 array, dimension ( MIN( M, N ) )
C             The leading  RANK  elements of TAU contain the scalar
C             factors of the elementary reflectors.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension ( 2*N )
C
C     ZWORK   COMPLEX*16 array, dimension ( 3*N-1 )
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
C     The routine computes a truncated QR factorization with column
C     pivoting of A,  A * P = Q * R,  with  R  defined above, and,
C     during this process, finds the largest leading submatrix whose
C     estimated condition number is less than 1/RCOND, taking the
C     possible positive value of SVLMAX into account.  This is performed
C     using the LAPACK incremental condition estimation scheme and a
C     slightly modified rank decision test.  The factorization process
C     stops when  RANK  has been determined.
C
C     The matrix Q is represented as a product of elementary reflectors
C
C        Q = H(1) H(2) . . . H(k), where k = rank <= min(m,n).
C
C     Each H(i) has the form
C
C        H = I - tau * v * v'
C
C     where tau is a complex scalar, and v is a complex vector with
C     v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
C     A(i+1:m,i), and tau in TAU(i).
C
C     The matrix P is represented in jpvt as follows: If
C        jpvt(j) = i
C     then the jth column of P is the ith canonical unit vector.
C
C     REFERENCES
C
C     [1] Bischof, C.H. and P. Tang.
C         Generalizing Incremental Condition Estimation.
C         LAPACK Working Notes 32, Mathematics and Computer Science
C         Division, Argonne National Laboratory, UT, CS-91-132,
C         May 1991.
C
C     [2] Bischof, C.H. and P. Tang.
C         Robust Incremental Condition Estimation.
C         LAPACK Working Notes 33, Mathematics and Computer Science
C         Division, Argonne National Laboratory, UT, CS-91-133,
C         May 1991.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998.
C     Complex version: V. Sima, Research Institute for Informatics,
C     Bucharest, Nov. 2008.
C
C     REVISIONS
C
C     V. Sima, Jan. 2010, following Bujanovic and Drmac's suggestion.
C
C     KEYWORDS
C
C     Eigenvalue problem, matrix operations, unitary transformation,
C     singular values.
C
C    ******************************************************************
C
#endif

      use omp_lib

      implicit none
!C     .. Parameters ..
      INTEGER            IMAX, IMIN
      PARAMETER          ( IMAX = 1, IMIN = 2 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), &
                         CONE  = ( 1.0D+0, 0.0D+0 ) )
!C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N, RANK
      DOUBLE PRECISION   RCOND, SVLMAX
!C     .. Array Arguments ..

      INTEGER            JPVT( * )
#if defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER            JPVT( * )
      !DIR$ ASSUME_ALIGNED JPVT:64
#endif
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
#if defined(__ICC) || defined(__INTEL_COMPILER)
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
      !DIR$ ASSUME_ALIGNED A:64
      !DIR$ ASSUME_ALIGNED TAU:64
      !DIR$ ASSUME_ALIGNED ZWORK:64
#endif
      DOUBLE PRECISION   DWORK( * ), SVAL( 3 )
#if defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION   DWORK( * ), SVAL( 3 )
      !DIR$ ASSUME_ALIGNED DWORK:64
#endif
!C     ..
!C     .. Local Scalars ..
      INTEGER            I, ISMAX, ISMIN, ITEMP, J, MN, PVT
      COMPLEX*16         AII, C1, C2, S1, S2
      DOUBLE PRECISION   SMAX, SMAXPR, SMIN, SMINPR, TEMP, TEMP2, TOLZ
!C     ..
!C     .. External Functions ..
      INTEGER            IDAMAX
     
      
!C     .. External Subroutines ..
      EXTERNAL           ZLAIC1, ZLARF, ZLARFG, ZSCAL, ZSWAP
!C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCONJG, MAX, MIN, SQRT
!C     ..
!C     .. Executable Statements ..
!C
!C     Test the input scalar arguments.
!C
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( RCOND.LT.ZERO .OR. RCOND.GT.ONE ) THEN
         INFO = -5
      ELSE IF( SVLMAX.LT.ZERO ) THEN
         INFO = -6
      END IF
!C
      IF( INFO.NE.0 ) THEN
          RETURN
      END IF
!C
!C     Quick return if possible.
!C
      MN = MIN( M, N )
      IF( MN.EQ.0 ) THEN
         RANK = 0
         SVAL( 1 ) = ZERO
         SVAL( 2 ) = ZERO
         SVAL( 3 ) = ZERO
         RETURN
      END IF
!C
      TOLZ  = SQRT( DLAMCH( 'Epsilon' ) )
      ISMIN = 1
      ISMAX = ISMIN + N
!C
!C     Initialize partial column norms and pivoting vector. The first n
!C     elements of DWORK store the exact column norms.
      !C

     
!$OMP SIMD ALIGNED(DWORK:64,JPVT:64)
      DO 10 I = 1, N
         DWORK( I ) = DZNRM2( M, A( 1, I ), 1 )
         DWORK( N+I ) = DWORK( I )
         JPVT( I ) = I
   10 CONTINUE
!C
!C     Compute factorization and determine RANK using incremental
!C     condition estimation.
!C
      RANK = 0
!C
   20 CONTINUE
      IF( RANK.LT.MN ) THEN
         I = RANK + 1
!C
!C        Determine ith pivot column and swap if necessary.
!C
         PVT = ( I-1 ) + IDAMAX( N-I+1, DWORK( I ), 1 )
!C
         IF( PVT.NE.I ) THEN
            CALL ZSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( I )
            JPVT( I )   = ITEMP
            DWORK( PVT )   = DWORK( I )
            DWORK( N+PVT ) = DWORK( N+I )
         END IF
!C
!C        Save A(I,I) and generate elementary reflector H(i)
!C        such that H(i)'*[A(i,i);*] = [*;0].
!C
         IF( I.LT.M ) THEN
            AII = A( I, I )
            CALL ZLARFG( M-I+1, A( I, I ), A( I+1, I ), 1, TAU( I ) )
         ELSE
            TAU( M ) = CZERO
         END IF
!C
         IF( RANK.EQ.0 ) THEN
!C
!C           Initialize; exit if matrix is zero (RANK = 0).
!C
            SMAX = ABS( A( 1, 1 ) )
            IF ( SMAX.EQ.ZERO ) THEN
               SVAL( 1 ) = ZERO
               SVAL( 2 ) = ZERO
               SVAL( 3 ) = ZERO
               RETURN
            END IF
            SMIN = SMAX
            SMAXPR = SMAX
            SMINPR = SMIN
            C1 = CONE
            C2 = CONE
         ELSE
!C
!C           One step of incremental condition estimation.
!C
            CALL ZLAIC1( IMIN, RANK, ZWORK( ISMIN ), SMIN, A( 1, I ), &
                         A( I, I ), SMINPR, S1, C1 )
            CALL ZLAIC1( IMAX, RANK, ZWORK( ISMAX ), SMAX, A( 1, I ), &
                         A( I, I ), SMAXPR, S2, C2 )
         END IF

         IF( SVLMAX*RCOND.LE.SMAXPR ) THEN
            IF( SVLMAX*RCOND.LE.SMINPR ) THEN
               IF( SMAXPR*RCOND.LE.SMINPR ) THEN
!!C
!C                 Continue factorization, as rank is at least RANK.
!C
                  IF( I.LT.N ) THEN
!C
!C                    Apply H(i)' to A(i:m,i+1:n) from the left.
!C
                     AII = A( I, I )
                     A( I, I ) = CONE
                     CALL ZLARF( 'Left', M-I+1, N-I, A( I, I ), 1,     &
                                DCONJG( TAU( I ) ), A( I, I+1 ), LDA,  &
                                ZWORK( 2*N+1 ) )
                     A( I, I ) = AII
                  END IF
!C
!C                 Update partial column norms.
                  !C
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(J,TEMP,TEMP2)
                  DO 30 J = I + 1, N
                     IF( DWORK( J ).NE.ZERO ) THEN
                        TEMP = ABS( A( I, J ) ) / DWORK( J )
                        TEMP = MAX( ( ONE + TEMP )*( ONE - TEMP ), ZERO)
                        TEMP2 = TEMP*( DWORK( J ) / DWORK( N+J ) )**2
                        IF( TEMP2.LE.TOLZ ) THEN
                           IF( M-I.GT.0 ) THEN
                              DWORK( J ) = DZNRM2( M-I, A( I+1, J ), 1 )
                              DWORK( N+J ) = DWORK( J )
                           ELSE
                              DWORK( J )   = ZERO
                              DWORK( N+J ) = ZERO
                           END IF
                        ELSE
                           DWORK( J ) = DWORK( J )*SQRT( TEMP )
                        END IF
                     END IF
30                   CONTINUE
!$OMP END PARALLEL DO

                  DO 40 I = 1, RANK
                     ZWORK( ISMIN+I-1 ) = S1*ZWORK( ISMIN+I-1 )
                     ZWORK( ISMAX+I-1 ) = S2*ZWORK( ISMAX+I-1 )
   40             CONTINUE
!C
                  ZWORK( ISMIN+RANK ) = C1
                  ZWORK( ISMAX+RANK ) = C2
                  SMIN = SMINPR
                  SMAX = SMAXPR
                  RANK = RANK + 1
                  GO TO 20
               END IF
            END IF
         END IF
      END IF
!C
!C     Restore the changed part of the (RANK+1)-th column and set SVAL.
!C
      IF ( RANK.LT.N ) THEN
         IF ( I.LT.M ) THEN
            CALL ZSCAL( M-I, -A( I, I )*TAU( I ), A( I+1, I ), 1 )
            A( I, I ) = AII
         END IF
      END IF
      IF ( RANK.EQ.0 ) THEN
         SMIN = ZERO
         SMINPR = ZERO
      END IF
      SVAL( 1 ) = SMAX
      SVAL( 2 ) = SMIN
      SVAL( 3 ) = SMINPR
!C
     
END SUBROUTINE MB3OYZ

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE MB3PYZ( M, N, A, LDA, RCOND, SVLMAX, RANK, SVAL, JPVT,
  TAU, DWORK, ZWORK, INFO ) !GCC$ ATTRIBUTES hot :: MB3PYZ !GCC$ ATTRIBUTES aligned(32) :: MB3PYZ !GCC$ ATTRIBUTES no_stack_protector :: MB3PYZ
#elif defined(__ICC) || defined(__INTEL_COMPILER))
SUBROUTINE MB3PYZ( M, N, A, LDA, RCOND, SVLMAX, RANK, SVAL, JPVT,
  TAU, DWORK, ZWORK, INFO )
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: MB3PYZ
  !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: MB3PYZ
#endif
#if 0
!C
!C     SLICOT RELEASE 5.7.
!C
!C     Copyright (c) 2002-2020 NICONET e.V.
!C
!C     PURPOSE
!C
!C     To compute a rank-revealing RQ factorization of a complex general
!C     M-by-N matrix  A,  which may be rank-deficient, and estimate its
!C     effective rank using incremental condition estimation.
!C
!C     The routine uses a truncated RQ factorization with row pivoting:
C                                   [ R11 R12 ]
C        P * A = R * Q,  where  R = [         ],
C                                   [  0  R22 ]
C     with R22 defined as the largest trailing upper triangular
C     submatrix whose estimated condition number is less than 1/RCOND.
C     The order of R22, RANK, is the effective rank of A.  Condition
C     estimation is performed during the RQ factorization process.
C     Matrix R11 is full (but of small norm), or empty.
C
C     MB3PYZ  does not perform any scaling of the matrix A.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix A.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix A.  N >= 0.
C
C     A       (input/output) COMPLEX*16 array, dimension ( LDA, N )
C             On entry, the leading M-by-N part of this array must
C             contain the given matrix A.
C             On exit, the upper triangle of the subarray
C             A(M-RANK+1:M,N-RANK+1:N) contains the RANK-by-RANK upper
C             triangular matrix R22;  the remaining elements in the last
C             RANK  rows, with the array TAU, represent the unitary
C             matrix Q as a product of  RANK  elementary reflectors
C             (see METHOD).  The first  M-RANK  rows contain the result
C             of the RQ factorization process used.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,M).
C
C     RCOND   (input) DOUBLE PRECISION
C             RCOND is used to determine the effective rank of A, which
C             is defined as the order of the largest trailing triangular
C             submatrix R22 in the RQ factorization with pivoting of A,
C             whose estimated condition number is less than 1/RCOND.
C             0 <= RCOND <= 1.
C             NOTE that when SVLMAX > 0, the estimated rank could be
C             less than that defined above (see SVLMAX).
C
C     SVLMAX  (input) DOUBLE PRECISION
C             If A is a submatrix of another matrix B, and the rank
C             decision should be related to that matrix, then SVLMAX
C             should be an estimate of the largest singular value of B
C             (for instance, the Frobenius norm of B).  If this is not
C             the case, the input value SVLMAX = 0 should work.
C             SVLMAX >= 0.
C
C     RANK    (output) INTEGER
C             The effective (estimated) rank of A, i.e., the order of
C             the submatrix R22.
C
C     SVAL    (output) DOUBLE PRECISION array, dimension ( 3 )
C             The estimates of some of the singular values of the
C             triangular factor R:
C             SVAL(1): largest singular value of
C                      R(M-RANK+1:M,N-RANK+1:N);
C             SVAL(2): smallest singular value of
C                      R(M-RANK+1:M,N-RANK+1:N);
C             SVAL(3): smallest singular value of R(M-RANK:M,N-RANK:N),
C                      if RANK < MIN( M, N ), or of
C                      R(M-RANK+1:M,N-RANK+1:N), otherwise.
C             If the triangular factorization is a rank-revealing one
C             (which will be the case if the trailing rows were well-
C             conditioned), then SVAL(1) will also be an estimate for
C             the largest singular value of A, and SVAL(2) and SVAL(3)
C             will be estimates for the RANK-th and (RANK+1)-st singular
C             values of A, respectively.
C             By examining these values, one can confirm that the rank
C             is well defined with respect to the chosen value of RCOND.
C             The ratio SVAL(1)/SVAL(2) is an estimate of the condition
C             number of R(M-RANK+1:M,N-RANK+1:N).
C
C     JPVT    (output) INTEGER array, dimension ( M )
C             If JPVT(i) = k, then the i-th row of P*A was the k-th row
C             of A.
C
C     TAU     (output) COMPLEX*16 array, dimension ( MIN( M, N ) )
C             The trailing  RANK  elements of TAU contain the scalar
C             factors of the elementary reflectors.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension ( 2*M )
C
C     ZWORK   COMPLEX*16 array, dimension ( 3*M-1 )
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
C     The routine computes a truncated RQ factorization with row
C     pivoting of A,  P * A = R * Q,  with  R  defined above, and,
C     during this process, finds the largest trailing submatrix whose
C     estimated condition number is less than 1/RCOND, taking the
C     possible positive value of SVLMAX into account.  This is performed
C     using an adaptation of the LAPACK incremental condition estimation
C     scheme and a slightly modified rank decision test.  The
C     factorization process stops when  RANK  has been determined.
C
C     The matrix Q is represented as a product of elementary reflectors
C
C        Q = H(k-rank+1)' H(k-rank+2)' . . . H(k)', where k = min(m,n).
C
C     Each H(i) has the form
C
C        H = I - tau * v * v'
C
C     where tau is a complex scalar, and v is a complex vector with
C     v(n-k+i+1:n) = 0 and v(n-k+i) = 1; conjg(v(1:n-k+i-1)) is stored
C     on exit in A(m-k+i,1:n-k+i-1), and tau in TAU(i).
C
C     The matrix P is represented in jpvt as follows: If
C        jpvt(j) = i
C     then the jth row of P is the ith canonical unit vector.
C
C     REFERENCES
C
C     [1] Bischof, C.H. and P. Tang.
C         Generalizing Incremental Condition Estimation.
C         LAPACK Working Notes 32, Mathematics and Computer Science
C         Division, Argonne National Laboratory, UT, CS-91-132,
C         May 1991.
C
C     [2] Bischof, C.H. and P. Tang.
C         Robust Incremental Condition Estimation.
C         LAPACK Working Notes 33, Mathematics and Computer Science
C         Division, Argonne National Laboratory, UT, CS-91-133,
C         May 1991.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998.
C     Complex version: V. Sima, Research Institute for Informatics,
C     Bucharest, Nov. 2008.
C
!C     REVISIONS
!C
!C     V. Sima, Jan. 2010, following Bujanovic and Drmac's suggestion.
!C
!C     KEYWORDS
!C
!C     Eigenvalue problem, matrix operations, unitary transformation,
!C     singular values.
!C
!C    ******************************************************************
!C
#endif

      INTEGER            IMAX, IMIN
      PARAMETER          ( IMAX = 1, IMIN = 2 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
!C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N, RANK
      DOUBLE PRECISION   RCOND, SVLMAX
!C     .. Array Arguments ..
      INTEGER            JPVT( * )
#if defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER            JPVT( * )
!DIR$ ASSUME_ALIGNED JPVT:64
#endif
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
#if defined(__ICC) || defined(__INTEL_COMPILER)
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
!DIR$ ASSUME_ALIGNED A:64
!DIR$ ASSUME_ALIGNED TAU:64
!DIR$ ASSUME_ALIGNED ZWORK:64
#endif
      DOUBLE PRECISION   DWORK( * ), SVAL( 3 )
#if defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION   DWORK( * ), SVAL( 3 )
!DIR$ ASSUME_ALIGNED DWORK:64
#endif
!C     .. Local Scalars ..
      INTEGER            I, ISMAX, ISMIN, ITEMP, J, JWORK, K, MKI, NKI, &
                         PVT
      COMPLEX*16         AII, C1, C2, S1, S2
      DOUBLE PRECISION   SMAX, SMAXPR, SMIN, SMINPR, TEMP, TEMP2, TOLZ
!C     ..
!C     .. External Functions ..
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DLAMCH, DZNRM2
!      EXTERNAL           DLAMCH, DZNRM2, IDAMAX
!C     ..
!C     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZLACGV, ZLAIC1, ZLARF, ZLARFG, &
                         ZSCAL, ZSWAP
!C     ..
!C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!C     ..
!C     .. Executable Statements ..
!C
!C     Test the input scalar arguments.
!C
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( RCOND.LT.ZERO .OR. RCOND.GT.ONE ) THEN
         INFO = -5
      ELSE IF( SVLMAX.LT.ZERO ) THEN
         INFO = -6
      END IF
!C
      IF( INFO.NE.0 ) THEN
          RETURN
      END IF
!C
!C     Quick return if possible.
!C
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         RANK = 0
         SVAL( 1 ) = ZERO
         SVAL( 2 ) = ZERO
         SVAL( 3 ) = ZERO
         RETURN
      END IF
!C
      TOLZ  = SQRT( DLAMCH( 'Epsilon' ) )
      ISMIN = 1
      ISMAX = ISMIN + M
      JWORK = ISMAX + M
!C
!C     Initialize partial row norms and pivoting vector. The first m
!C     elements of DWORK store the exact row norms.
!C
!$OMP SIMD ALIGNED(DWORK:64,JPVT:64)
      DO 10 I = 1, M
         DWORK( I ) = DZNRM2( N, A( I, 1 ), LDA )
         DWORK( M+I ) = DWORK( I )
         JPVT( I ) = I
   10 CONTINUE
!C
!C     Compute factorization and determine RANK using incremental
!C     condition estimation.
!C
      RANK = 0
!C
   20 CONTINUE
      IF( RANK.LT.K ) THEN
         I = K - RANK
!C
!C        Determine ith pivot row and swap if necessary.
!C
         MKI = M - RANK
         NKI = N - RANK
         PVT = IDAMAX( MKI, DWORK, 1 )

         IF( PVT.NE.MKI ) THEN
            CALL ZSWAP( N, A( PVT, 1 ), LDA, A( MKI, 1 ), LDA )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( MKI )
            JPVT( MKI ) = ITEMP
            DWORK( PVT )   = DWORK( MKI )
            DWORK( M+PVT ) = DWORK( M+MKI )
         END IF

         IF( NKI.GT.1 ) THEN
!C
!C           Save A(m-k+i,n-k+i) and generate elementary reflector H(i)
!C           to annihilate A(m-k+i,1:n-k+i-1), k = min(m,n).
!C               A(m-k+i,1:n-k+i) * H(tau,v)        = [0 , *]         <=>
!C               H(conj(tau),v) A(m-k+i,1:n-k+i)^H  = [0 ; *],
!C           using H(tau,v)^H = H(conj(tau),v).
!C
            CALL ZLACGV( NKI, A( MKI, 1 ), LDA )
            AII = A( MKI, NKI )
            CALL ZLARFG( NKI, A( MKI, NKI ), A( MKI, 1 ), LDA, TAU( I ))
                       
         END IF
!C
         IF( RANK.EQ.0 ) THEN
!C
!C           Initialize; exit if matrix is zero (RANK = 0).
!C
            SMAX = ABS( A( M, N ) )
            IF ( SMAX.EQ.ZERO ) THEN
               SVAL( 1 ) = ZERO
               SVAL( 2 ) = ZERO
               SVAL( 3 ) = ZERO
               RETURN
            END IF
            SMIN = SMAX
            SMAXPR = SMAX
            SMINPR = SMIN
            C1 = CONE
            C2 = CONE
         ELSE
!C
!C           One step of incremental condition estimation.
!C
            CALL ZCOPY ( RANK, A( MKI, NKI+1 ), LDA, ZWORK( JWORK ), 1 )
            CALL ZLAIC1( IMIN, RANK, ZWORK( ISMIN ), SMIN, &
                        ZWORK( JWORK ), A( MKI, NKI ), SMINPR, S1, C1 )
            CALL ZLAIC1( IMAX, RANK, ZWORK( ISMAX ), SMAX, &
                        ZWORK( JWORK ), A( MKI, NKI ), SMAXPR, S2, C2 )
         END IF

         IF( SVLMAX*RCOND.LE.SMAXPR ) THEN
            IF( SVLMAX*RCOND.LE.SMINPR ) THEN
               IF( SMAXPR*RCOND.LE.SMINPR ) THEN

                  IF( MKI.GT.1 ) THEN
!C
!C                    Continue factorization, as rank is at least RANK.
!C                    Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right.
!C
                     AII = A( MKI, NKI )
                     A( MKI, NKI ) = CONE
                     CALL ZLARF( 'Right', MKI-1, NKI, A( MKI, 1 ), LDA, &
                                TAU( I ), A, LDA, ZWORK( JWORK ) )
                     A( MKI, NKI ) = AII
!C
!C                    Update partial row norms.
!C
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(J,TEMP,TEMP2)
                     DO 30 J = 1, MKI - 1
                        IF( DWORK( J ).NE.ZERO ) THEN
                           TEMP = ABS( A( J, NKI ) ) / DWORK( J )
                           TEMP = MAX( ( ONE + TEMP )*( ONE - TEMP ), &
                                       ZERO )
                           TEMP2 = TEMP*( DWORK( J ) / DWORK( M+J ) )**2
                           IF( TEMP2.LE.TOLZ ) THEN
                              DWORK( J )   = DZNRM2( NKI-1, A( J, 1 ), &
                                                     LDA )
                              DWORK( M+J ) = DWORK( J )
                           ELSE
                              DWORK( J ) = DWORK( J )*SQRT( TEMP )
                           END IF
                        END IF
   30                CONTINUE
!$OMP END PARALLEL DO

                  END IF

                  DO 40 I = 1, RANK
                     ZWORK( ISMIN+I-1 ) = S1*ZWORK( ISMIN+I-1 )
                     ZWORK( ISMAX+I-1 ) = S2*ZWORK( ISMAX+I-1 )
   40             CONTINUE

                  ZWORK( ISMIN+RANK ) = C1
                  ZWORK( ISMAX+RANK ) = C2
                  SMIN = SMINPR
                  SMAX = SMAXPR
                  RANK = RANK + 1
                  CALL ZLACGV( NKI-1, A( MKI, 1 ), LDA )
                  GO TO 20
               END IF
            END IF
         END IF
      END IF
!C
!C     Restore the changed part of the (M-RANK)-th row and set SVAL.
!C
      IF ( RANK.LT.K .AND. NKI.GT.1 ) THEN
         CALL ZLACGV( NKI-1, A( MKI, 1 ), LDA )
         CALL ZSCAL( NKI-1, -A( MKI, NKI )*TAU( I ), A( MKI, 1 ), LDA )
         A( MKI, NKI ) = AII
      END IF
      SVAL( 1 ) = SMAX
      SVAL( 2 ) = SMIN
      SVAL( 3 ) = SMINPR
END SUBROUTINE MB3PYZ
    
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE FB01QD( JOBK, MULTBQ, N, M, P, S, LDS, A, LDA, B,  &
                        LDB, Q, LDQ, C, LDC, R, LDR, K, LDK, TOL, &
                        IWORK, DWORK, LDWORK, INFO ) !GCC$ ATTRIBUTES hot :: FB01QD !GCC$ ATTRIBUTES aligned(32) :: FB01QD !GCC$ ATTRIBUTES no_stack_protector :: FB01QD
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE FB01QD( JOBK, MULTBQ, N, M, P, S, LDS, A, LDA, B,  &
                        LDB, Q, LDQ, C, LDC, R, LDR, K, LDK, TOL, &
                        IWORK, DWORK, LDWORK, INFO )
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: FB01QD
  !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: FB01QD
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To calculate a combined measurement and time update of one
C     iteration of the time-varying Kalman filter. This update is given
C     for the square root covariance filter, using dense matrices.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBK    CHARACTER*1
C             Indicates whether the user wishes to compute the Kalman
C             filter gain matrix K  as follows:
C                                 i
C             = 'K':  K  is computed and stored in array K;
C                      i
C             = 'N':  K  is not required.
C                      i
C
C     MULTBQ  CHARACTER*1                    1/2
C             Indicates how matrices B  and Q    are to be passed to
C                                     i      i
C             the routine as follows:
C             = 'P':  Array Q is not used and the array B must contain
C                                    1/2
C                     the product B Q   ;
C                                  i i
C             = 'N':  Arrays B and Q must contain the matrices as
C                     described below.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e., the order of the
C             matrices S    and A .  N >= 0.
C                       i-1      i
C
C     M       (input) INTEGER
C             The actual input dimension, i.e., the order of the matrix
C              1/2
C             Q   .  M >= 0.
C              i
C
C     P       (input) INTEGER
C             The actual output dimension, i.e., the order of the matrix
C              1/2
C             R   .  P >= 0.
C              i
C
C     S       (input/output) DOUBLE PRECISION array, dimension (LDS,N)
C             On entry, the leading N-by-N lower triangular part of this
C             array must contain S   , the square root (left Cholesky
C                                 i-1
C             factor) of the state covariance matrix at instant (i-1).
C             On exit, the leading N-by-N lower triangular part of this
C             array contains S , the square root (left Cholesky factor)
C                             i
C             of the state covariance matrix at instant i.
C             The strict upper triangular part of this array is not
C             referenced.
C
C     LDS     INTEGER
C             The leading dimension of array S.  LDS >= MAX(1,N).
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain A ,
C                                                                 i
C             the state transition matrix of the discrete system at
C             instant i.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain B ,
C                                                        1/2      i
C             the input weight matrix (or the product B Q    if
C                                                      i i
C             MULTBQ = 'P') of the discrete system at instant i.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     Q       (input) DOUBLE PRECISION array, dimension (LDQ,*)
C             If MULTBQ = 'N', then the leading M-by-M lower triangular
C                                              1/2
C             part of this array must contain Q   , the square root
C                                              i
C             (left Cholesky factor) of the input (process) noise
C             covariance matrix at instant i.
C             The strict upper triangular part of this array is not
C             referenced.
C             If MULTBQ = 'P', Q is not referenced and can be supplied
C             as a dummy array (i.e., set parameter LDQ = 1 and declare
C             this array to be Q(1,1) in the calling program).
C
C     LDQ     INTEGER
C             The leading dimension of array Q.
C             LDQ >= MAX(1,M) if MULTBQ = 'N';
C             LDQ >= 1        if MULTBQ = 'P'.
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain C , the
C                                                                 i
C             output weight matrix of the discrete system at instant i.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     R       (input/output) DOUBLE PRECISION array, dimension (LDR,P)
C             On entry, the leading P-by-P lower triangular part of this
C                                 1/2
C             array must contain R   , the square root (left Cholesky
C                                 i
C             factor) of the output (measurement) noise covariance
C             matrix at instant i.
C             On exit, the leading P-by-P lower triangular part of this
C                                    1/2
C             array contains (RINOV )   , the square root (left Cholesky
C                                  i
C             factor) of the covariance matrix of the innovations at
C             instant i.
C             The strict upper triangular part of this array is not
C             referenced.
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,P).
C
C     K       (output) DOUBLE PRECISION array, dimension (LDK,P)
C             If JOBK = 'K', and INFO = 0, then the leading N-by-P part
C             of this array contains K , the Kalman filter gain matrix
C                                     i
C             at instant i.
C             If JOBK = 'N', or JOBK = 'K' and INFO = 1, then the
C             leading N-by-P part of this array contains AK , a matrix
C                                                          i
C             related to the Kalman filter gain matrix at instant i (see
C                                                            -1/2
C             METHOD). Specifically, AK  = A P     C'(RINOV')    .
C                                      i    i i|i-1 i      i
C
C     LDK     INTEGER
C             The leading dimension of array K.   LDK >= MAX(1,N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If JOBK = 'K', then TOL is used to test for near
C                                               1/2
C             singularity of the matrix (RINOV )   . If the user sets
C                                             i
C             TOL > 0, then the given value of TOL is used as a
C             lower bound for the reciprocal condition number of that
C             matrix; a matrix whose estimated condition number is less
C             than 1/TOL is considered to be nonsingular. If the user
C             sets TOL <= 0, then an implicitly computed, default
C             tolerance, defined by TOLDEF = P*P*EPS, is used instead,
C             where EPS is the machine precision (see LAPACK Library
C             routine DLAMCH).
C             Otherwise, TOL is not referenced.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK),
C             where LIWORK = P if JOBK = 'K',
C             and   LIWORK = 1 otherwise.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.  If INFO = 0 and JOBK = 'K', DWORK(2) returns
C             an estimate of the reciprocal of the condition number
C                                        1/2
C             (in the 1-norm) of (RINOV )   .
C                                      i
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= MAX(1,N*(P+N)+2*P,N*(N+M+2)),     if JOBK = 'N';
C             LDWORK >= MAX(2,N*(P+N)+2*P,N*(N+M+2),3*P), if JOBK = 'K'.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C                                                        1/2
C             = 1:  if JOBK = 'K' and the matrix (RINOV )   is singular,
C                                                      i           1/2
C                   i.e., the condition number estimate of (RINOV )
C                                                                i
C                   (in the 1-norm) exceeds 1/TOL.  The matrices S, AK ,
C                               1/2                                   i
C                   and (RINOV )    have been computed.
C                             i
C
C     METHOD
C
C     The routine performs one recursion of the square root covariance
C     filter algorithm, summarized as follows:
C
C      |  1/2                      |     |         1/2          |
C      | R      C x S      0       |     | (RINOV )     0     0 |
C      |  i      i   i-1           |     |       i              |
C      |                      1/2  | T = |                      |
C      | 0      A x S    B x Q     |     |     AK       S     0 |
C      |         i   i-1  i   i    |     |       i       i      |
C
C          (Pre-array)                      (Post-array)
C
C     where T is an orthogonal transformation triangularizing the
C     pre-array.
C
C     The state covariance matrix P    is factorized as
C                                  i|i-1
C        P     = S  S'
!C         i|i-1   i  i
!C
!C     and one combined time and measurement update for the state X
C                                                                 i|i-1
C     is given by
C
C        X     = A X      + K (Y - C X     ),
C         i+1|i   i i|i-1    i  i   i i|i-1
C
C                          -1/2
C     where K = AK (RINOV )     is the Kalman filter gain matrix and Y
C            i    i      i                                            i
C     is the observed output of the system.
C
C     The triangularization is done entirely via Householder
C     transformations exploiting the zero pattern of the pre-array.
C
C     REFERENCES
C
C     [1] Anderson, B.D.O. and Moore, J.B.
C         Optimal Filtering.
C         Prentice Hall, Englewood Cliffs, New Jersey, 1979.
C
C     [2] Verhaegen, M.H.G. and Van Dooren, P.
C         Numerical Aspects of Different Kalman Filter Implementations.
C         IEEE Trans. Auto. Contr., AC-31, pp. 907-917, Oct. 1986.
C
!C     [3] Vanbegin, M., Van Dooren, P., and Verhaegen, M.H.G.
!C         Algorithm 675: FORTRAN Subroutines for Computing the Square
!C         Root Covariance Filter and Square Root Information Filter in
C         Dense or Hessenberg Forms.
C         ACM Trans. Math. Software, 15, pp. 243-256, 1989.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires
C
C           3    2                               2   2
C     (7/6)N  + N  x (5/2 x P + M) + N x (1/2 x M + P )
C
C     operations and is backward stable (see [2]).
C
C     CONTRIBUTORS
!C
!C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
!C     Supersedes Release 2.0 routine FB01ED by M. Vanbegin,
!C     P. Van Dooren, and M.H.G. Verhaegen.
!C
!C     REVISIONS
!C
!C     February 20, 1998, November 20, 2003.
!C
!C     KEYWORDS
!C
!C     Kalman filtering, optimal filtering, orthogonal transformation,
!C     recursive estimation, square-root covariance filtering,
!C     square-root filtering.
!C
!C     ******************************************************************
!C
#endif 
      implicit none

!C     .. Parameters ..
      DOUBLE PRECISION  ONE, TWO
      PARAMETER         ( ONE = 1.0D0, TWO = 2.0D0 )
!C     .. Scalar Arguments ..
      CHARACTER         JOBK, MULTBQ
      INTEGER           INFO, LDA, LDB, LDC, LDK, LDQ, LDR, LDS, LDWORK, &
                        M, N, P
      DOUBLE PRECISION  TOL
!C     .. Array Arguments ..
      INTEGER           IWORK(*)
#if defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER           IWORK(*)
!DIR$ ASSUME_ALIGNED IWORK:64
#endif
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
                        K(LDK,*), Q(LDQ,*), R(LDR,*), S(LDS,*)
#if defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
                        K(LDK,*), Q(LDQ,*), R(LDR,*), S(LDS,*)
!DIR$ ASSUME_ALIGNED A:64
!DIR$ ASSUME_ALIGNED B:64
!DIR$ ASSUME_ALIGNED C:64
!DIR$ ASSUME_ALIGNED DWORK:64
!DIR$ ASSUME_ALIGNED K:64
!DIR$ ASSUME_ALIGNED Q:64
!DIR$ ASSUME_ALIGNED R:64
!DIR$ ASSUME_ALIGNED S:64
#endif
!C     .. Local Scalars ..
      LOGICAL           LJOBK, LMULTB
      INTEGER           I12, ITAU, JWORK, N1, PN, WRKOPT
      DOUBLE PRECISION  RCOND

     
!C     .. External Subroutines ..
      EXTERNAL          DGELQF, DLACPY, DTRMM, MB02OD, MB04LD
!C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX
!C     .. Executable Statements ..
!C
      PN = P + N
      N1 = MAX( 1, N )
      INFO = 0
      LJOBK  = LSAME( JOBK, 'K' )
      LMULTB = LSAME( MULTBQ, 'P' )
!C
!C     Test the input scalar arguments.
!C
      IF( .NOT.LJOBK .AND. .NOT.LSAME( JOBK, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LMULTB .AND. .NOT.LSAME( MULTBQ, 'N' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDS.LT.N1 ) THEN
         INFO = -7
      ELSE IF( LDA.LT.N1 ) THEN
         INFO = -9
      ELSE IF( LDB.LT.N1 ) THEN
         INFO = -11
      ELSE IF( LDQ.LT.1 .OR. ( .NOT.LMULTB .AND. LDQ.LT.M ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( LDR.LT.MAX( 1, P ) ) THEN
         INFO = -17
      ELSE IF( LDK.LT.N1 ) THEN
         INFO = -19
      ELSE IF( ( LJOBK .AND. LDWORK.LT.MAX( 2, PN*N + 2*P,  &
                                           N*(N + M + 2), 3*P ) ) .OR. &
           ( .NOT.LJOBK .AND. LDWORK.LT.MAX( 1, PN*N + 2*P, &
                                           N*(N + M + 2) ) ) ) THEN
         INFO = -23
      END IF
!C
      IF ( INFO.NE.0 ) THEN
!C
!C        Error return.
!C
         RETURN
      END IF
!C
!C     Quick return if possible.
!C
      IF ( N.EQ.0 ) THEN
         IF ( LJOBK ) THEN
            DWORK(1) = TWO
            DWORK(2) = ONE
         ELSE
            DWORK(1) = ONE
         END IF
         RETURN
      END IF
!C
!C     Construction of the needed part of the pre-array in DWORK.
!C     To save workspace, only the blocks (1,2), (2,2), and (2,3) will be
!C     constructed as shown below.
!C
!C     Storing A x S and C x S in the (1,1) and (2,1) blocks of DWORK,
!C     respectively.
!C     Workspace: need (N+P)*N.
!C
!C     (Note: Comments in the code beginning "Workspace:" describe the
!C     minimal amount of real workspace needed at that point in the
!C     code, as well as the preferred amount for good performance.
!C     NB refers to the optimal block size for the immediately
!C     following subroutine, as returned by ILAENV.)
!C
      CALL DLACPY( 'Full', N, N, A, LDA, DWORK, PN )
      CALL DLACPY( 'Full', P, N, C, LDC, DWORK(N+1), PN )
      CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Non-unit', PN, N, &
                  ONE, S, LDS, DWORK, PN )
!C
!C     Triangularization (2 steps).
!C
!C     Step 1: annihilate the matrix C x S.
!C     Workspace: need (N+P)*N + 2*P.
!C
      ITAU  = PN*N + 1
      JWORK = ITAU + P
!C
      CALL MB04LD( 'Full', P, N, N, R, LDR, DWORK(N+1), PN, DWORK, PN,  &
                   K, LDK, DWORK(ITAU), DWORK(JWORK) )
      WRKOPT = PN*N + 2*P
!C
!C     Now, the workspace for C x S is no longer needed.
!C     Adjust the leading dimension of DWORK, to save space for the
!C     following computations.
!C
      CALL DLACPY( 'Full', N, N, DWORK, PN, DWORK, N )
      I12 = N*N + 1
!C
!C     Storing B x Q in the (1,2) block of DWORK.
!C     Workspace: need N*(N+M).
!C
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK(I12), N )
      IF ( .NOT.LMULTB )  &
         CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Non-unit', N, M,  &
                     ONE, Q, LDQ, DWORK(I12), N )
      WRKOPT = MAX( WRKOPT, N*( N + M ) )
!C
!C     Step 2: LQ triangularization of the matrix [ A x S  B x Q ], where
!C     A x S was modified at Step 1.
!C     Workspace: need N*(N+M+2);  prefer N*(N+M+1)+N*NB.
!C
      ITAU  = N*( N + M ) + 1
      JWORK = ITAU + N
!C
      CALL DGELQF( N, N+M, DWORK, N, DWORK(ITAU), DWORK(JWORK), &
                  LDWORK-JWORK+1, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
!C
!C     Output S and K (if needed) and set the optimal workspace
!C     dimension (and the reciprocal of the condition number estimate).
!C
      CALL DLACPY( 'Lower', N, N, DWORK, N, S, LDS )
!C
      IF ( LJOBK ) THEN
!C
!C        Compute K.
!C        Workspace: need 3*P.
!C
         CALL MB02OD( 'Right', 'Lower', 'No transpose', 'Non-unit',     &
                     '1-norm', N, P, ONE, R, LDR, K, LDK, RCOND, TOL,   &
                      IWORK, DWORK, INFO )
         IF ( INFO.EQ.0 ) THEN
            WRKOPT = MAX( WRKOPT, 3*P )
            DWORK(2) = RCOND
         END IF
      END IF
!C
      DWORK(1) = WRKOPT
!C
    
C *** Last line of FB01QD ***
END SUBROUTINE FB01QD



!Helpers





DOUBLE PRECISION FUNCTION DZNRM2(N,X,INCX)
  implicit none
!$OMP DECLARE SIMD(DZNRM2) UNIFORM(X), LINEAR(IX:1)
  
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
      COMPLEX*16 X(*)
#if defined(__ICC) || defined(__INTEL_COMPILER)
      COMPLEX*16 X(*)
      !DIR$ ASSUME_ALIGNED X:64
#endif
!*     ..
!*
!*  =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION NORM,SCALE,SSQ,TEMP
      INTEGER IX
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG,SQRT
!*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE
          SCALE = ZERO
          SSQ = ONE
!*        The following loop is equivalent to this call to the LAPACK
!*        auxiliary routine:
!*        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
!*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (DBLE(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(DBLE(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
              IF (DIMAG(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(DIMAG(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF

      DZNRM2 = NORM
END FUNCTION DZNRM2



#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
INTEGER FUNCTION IDAMAX(N,DX,INCX) !GCC$ ATTRIBUTES aligned(32) :: IDAMAX !GCC$ ATTRIBUTES PURE :: IDAMAX !GCC$ ATTRIBUTES INLINE :: IDAMAX
#elif defined(__ICC) || defined(__INTEL_COMPILER)
INTEGER FUNCTION IDAMAX(N,DX,INCX)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: IDAMAX
!DIR$ ATTRIBUTES FORCEINLINE :: IDAMAX
#endif
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
      DOUBLE PRECISION DX(*)
#if defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION DX(*)
!DIR$ ASSUME_ALIGNED DX:64
#endif

!*     ..
!*
!*  =====================================================================
!*
!*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DABS
!*     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
!*
!*        code for increment equal to 1
!*
         DMAX = DABS(DX(1))
         DO I = 2,N
            IF (DABS(DX(I)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(I))
            END IF
         END DO
      ELSE
!*
!*        code for increment not equal to 1
!*
         IX = 1
         DMAX = DABS(DX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DABS(DX(IX)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
     
END FUNCTION IDAMAX


#if defined __GFORTRAN__ && (!defined(__INTEL_COMPILER) || !defined(__ICC))
      LOGICAL          FUNCTION LSAME( CA, CB ) !GCC$ ATTRIBUTES inline :: LSAME
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      LOGICAL          FUNCTION LSAME( CA, CB)
      !DIR$ ATTRIBUTES INLINE :: LSAME
#endif

!*
!*  -- LAPACK auxiliary routine (preliminary version) --
!*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
!*     Courant Institute, NAG Ltd., and Rice University
!*     March 26, 1990
!*
!*     .. Scalar Arguments ..
      CHARACTER          CA, CB
!*     ..
!*
!*  Purpose
!*  =======
!*
!*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!*  case.
!*
!*  This version of the routine is only correct for ASCII code.
!*  Installers must modify the routine for other character-codes.
!*
!*  For EBCDIC systems the constant IOFF must be changed to -64.
!*  For CDC systems using 6-12 bit representations, the system-
!*  specific code in comments must be activated.
!*
!*  Arguments
!*  =========
!*
!*  CA      (input) CHARACTER*1
!*  CB      (input) CHARACTER*1
!*          CA and CB specify the single characters to be compared.
!*
!*
!*     .. Parameters ..
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

END FUNCTION LSAME




*!> DLAMCH determines double precision machine parameters.
!*> \endverbatim
!*
!*  Arguments:
!*  ==========
!*
!*> \param[in] CMACH
!*> \verbatim
!*>          CMACH is CHARACTER*1
!*>          Specifies the value to be returned by DLAMCH:
!*>          = 'E' or 'e',   DLAMCH := eps
!*>          = 'S' or 's ,   DLAMCH := sfmin
!*>          = 'B' or 'b',   DLAMCH := base
!*>          = 'P' or 'p',   DLAMCH := eps*base
!*>          = 'N' or 'n',   DLAMCH := t
!*>          = 'R' or 'r',   DLAMCH := rnd
!*>          = 'M' or 'm',   DLAMCH := emin
!*>          = 'U' or 'u',   DLAMCH := rmin
!*>          = 'L' or 'l',   DLAMCH := emax
!*>          = 'O' or 'o',   DLAMCH := rmax
!*>          where
!*>          eps   = relative machine precision
!*>          sfmin = safe minimum, such that 1/sfmin does not overflow
!*>          base  = base of the machine
!*>          prec  = eps*base
!*>          t     = number of (base) digits in the mantissa
!*>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!*>          emin  = minimum exponent before (gradual) underflow
!*>          rmin  = underflow threshold - base**(emin-1)
!*>          emax  = largest exponent before overflow
!*>          rmax  = overflow threshold  - (base**emax)*(1-eps)
!*> \endverbatim
!*
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
!*> \ingroup auxOTHERauxiliary
!*
!*  =====================================================================

      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!*     December 2016
!*
!*     .. Scalar Arguments ..
      CHARACTER          CMACH
!*     ..
!*
!* =====================================================================
!*
!*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     ..
!*     .. Local Scalars ..
      DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH
!*     ..
!*     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT, &
                         MINEXPONENT, RADIX, TINY
!*     ..
!*     .. Executable Statements ..
!*
!*
!*     Assume rounding, not chopping. Always.
!*
      RND = ONE
!*
      IF( ONE.EQ.RND ) THEN
         EPS = EPSILON(ZERO) * 0.5
      ELSE
         EPS = EPSILON(ZERO)
      END IF
!*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         SFMIN = TINY(ZERO)
         SMALL = ONE / HUGE(ZERO)
         IF( SMALL.GE.SFMIN ) THEN
!*
!*           Use SMALL plus a bit, to avoid the possibility of rounding
!*           causing overflow when computing  1/sfmin.
!*
            SFMIN = SMALL*( ONE+EPS )
         END IF
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = EPS * RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = DIGITS(ZERO)
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = MINEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = tiny(zero)
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = MAXEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = HUGE(ZERO)
      ELSE
         RMACH = ZERO
      END IF

      DLAMCH = RMACH
     

END FUNCTION DLAMCH

DOUBLE PRECISION FUNCTION DLAMC3( A, B )
!*
!*  -- LAPACK auxiliary routine (version 3.7.0) --
!*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!*     November 2010
!*
!*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
!*     ..
!* =====================================================================
!*
!*     .. Executable Statements ..
!*
      DLAMC3 = A + B
!*
   

END FUNCTION DLAMC3

