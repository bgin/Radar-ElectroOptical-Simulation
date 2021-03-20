
! IMPORTANT NOTICE!!
! Set the USE_OMP macro to 0 when linked against multithreaded version of Intel MKL library
!
#if !defined(GMS_SLICOT_OMP_LOOP_PARALLELIZE)
#define GMS_SLICOT_OMP_LOOP_PARALLELIZE 1
#endif

#if !defined(GMS_SLICOT_USE_MKL_LAPACK)
#define GMS_SLICOT_USE_MKL_LAPACK 1
#endif

#if !defined(GMS_SLICOT_USE_REFERENCE_LAPACK)
#define GMS_SLICOT_USE_REFERENCE_LAPACK 0
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
!C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
!C     .. Scalar Arguments ..
      CHARACTER         EQUIL
      INTEGER           INFO, LDA, LDB, LDC, LDD, LZWORK, M, N, P, RANK
      DOUBLE PRECISION  TOL
      !C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      INTEGER           IWORK(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER           IWORK(*)
      !DIR$ ASSUME_ALIGNED IWORK:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*), ZWORK(*)
      DOUBLE PRECISION  DWORK(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
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
!C                  ( D' C')
!C     having the same transmission zeros but with D' of full row rank.
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
!C             The dimension of the reduced system matrix and the number
!C             of (finite) invariant zeros if D' is invertible.
!C
!C     NKROL   (output) INTEGER
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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER)) 
      INTEGER           INFZ(*), IWORK(*), KRONL(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER           INFZ(*), IWORK(*), KRONL(*)
      !DIR$ ASSUME_ALIGNED INFZ:64
      !DIR$ ASSUME_ALIGNED IWORK:64
      !DIR$ ASSUME_ALIGNED KRONL:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      COMPLEX*16        ABCD(LDABCD,*), ZWORK(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      COMPLEX*16        ABCD(LDABCD,*), ZWORK(*)
      !DIR$ ASSUME_ALIGNED ABCD:64
      !DIR$ ASSUME_ALIGNED ZWORK:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION  DWORK(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
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
#endif
            DO 40 I1 = 1, SIGMA
               CALL ZLARFG( RO+1, ABCD(IROW,I1), ABCD(IROW+1,I1), 1, &
                           TC )
               CALL ZLATZM( 'L', RO+1, MNU-I1, ABCD(IROW+1,I1), 1, &
                           DCONJG( TC ), ABCD(IROW,I1+1),          &
                           ABCD(IROW+1,I1+1), LDABCD, ZWORK )
               IROW = IROW + 1
40             CONTINUE

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
!C
!C        H = I - tau * v * v'
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
!C     V. Sima, Jan. 2010, following Bujanovic and Drmac's suggestion.
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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      INTEGER            JPVT( * )
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER            JPVT( * )
      !DIR$ ASSUME_ALIGNED JPVT:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
      !DIR$ ASSUME_ALIGNED A:64
      !DIR$ ASSUME_ALIGNED TAU:64
      !DIR$ ASSUME_ALIGNED ZWORK:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION   DWORK( * ), SVAL( 3 )
#elif defined(__ICC) || defined(__INTEL_COMPILER)
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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      INTEGER            JPVT( * )
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER            JPVT( * )
!DIR$ ASSUME_ALIGNED JPVT:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      COMPLEX*16         A( LDA, * ), TAU( * ), ZWORK( * )
!DIR$ ASSUME_ALIGNED A:64
!DIR$ ASSUME_ALIGNED TAU:64
!DIR$ ASSUME_ALIGNED ZWORK:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION   DWORK( * ), SVAL( 3 )
#elif defined(__ICC) || defined(__INTEL_COMPILER)
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
!C
!C          (Pre-array)                      (Post-array)
!C
!C     where T is an orthogonal transformation triangularizing the
!C     pre-array.
!C
!C     The state covariance matrix P    is factorized as
!C                                  i|i-1
!C        P     = S  S'
!C         i|i-1   i  i
!C
!C     and one combined time and measurement update for the state X
!C                                                                 i|i-1
!C     is given by
!C
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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      INTEGER           IWORK(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER           IWORK(*)
!DIR$ ASSUME_ALIGNED IWORK:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
                        K(LDK,*), Q(LDQ,*), R(LDR,*), S(LDS,*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
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
    

END SUBROUTINE FB01QD

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTRL_COMPILER))
SUBROUTINE FB01SD( JOBX, MULTAB, MULTRC, N, M, P, SINV, LDSINV, &
          AINV, LDAINV, B, LDB, RINV, LDRINV, C, LDC,           &
          QINV, LDQINV, X, RINVY, Z, E, TOL, IWORK,             &
          DWORK, LDWORK, INFO ) !GCC$ ATTRIBUTES hot :: FB01SD !GCC$ ATTRIBUTES aligned(32) :: FB01SD !GCC$ ATTRIBUTES no_stack_protector :: FB01SD
#elif defined(__ICC)  || defined(__INTEL_COMPILER)
SUBROUTINE FB01SD( JOBX, MULTAB, MULTRC, N, M, P, SINV, LDSINV, &
          AINV, LDAINV, B, LDB, RINV, LDRINV, C, LDC,           &
          QINV, LDQINV, X, RINVY, Z, E, TOL, IWORK,             &
          DWORK, LDWORK, INFO )
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: FB01SD
  !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: FB01SD
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
C     for the square root information filter, using dense matrices.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBX    CHARACTER*1
C             Indicates whether X    is to be computed as follows:
C                                i+1
C             = 'X':  X    is computed and stored in array X;
C                      i+1
C             = 'N':  X    is not required.
C                      i+1
C
C     MULTAB  CHARACTER*1             -1
C             Indicates how matrices A   and B  are to be passed to
C                                     i       i
C             the routine as follows:                       -1
C             = 'P':  Array AINV must contain the matrix   A    and the
C                                                       -1  i
C                     array B must contain the product A  B ;
C                                                       i  i
C             = 'N':  Arrays AINV and B must contain the matrices
C                     as described below.
C
C     MULTRC  CHARACTER*1             -1/2
C             Indicates how matrices R     and C    are to be passed to
C                                     i+1       i+1
C             the routine as follows:
C             = 'P':  Array RINV is not used and the array C must
C                                          -1/2
C                     contain the product R    C   ;
C                                          i+1  i+1
C             = 'N':  Arrays RINV and C must contain the matrices
C                     as described below.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e., the order of the
C                       -1      -1
C             matrices S   and A  .  N >= 0.
C                       i       i
C
C     M       (input) INTEGER
C             The actual input dimension, i.e., the order of the matrix
C              -1/2
C             Q    .  M >= 0.
C              i
C
C     P       (input) INTEGER
C             The actual output dimension, i.e., the order of the matrix
C              -1/2
C             R    .  P >= 0.
C              i+1
C
C     SINV    (input/output) DOUBLE PRECISION array, dimension
C             (LDSINV,N)
C             On entry, the leading N-by-N upper triangular part of this
C                                 -1
C             array must contain S  , the inverse of the square root
C                                 i
C             (right Cholesky factor) of the state covariance matrix
C             P    (hence the information square root) at instant i.
C              i|i
C             On exit, the leading N-by-N upper triangular part of this
C                             -1
C             array contains S   , the inverse of the square root (right
C                             i+1
C             Cholesky factor) of the state covariance matrix P
C                                                              i+1|i+1
C             (hence the information square root) at instant i+1.
C             The strict lower triangular part of this array is not
C             referenced.
C
C     LDSINV  INTEGER
C             The leading dimension of array SINV.  LDSINV >= MAX(1,N).
C
C     AINV    (input) DOUBLE PRECISION array, dimension (LDAINV,N)
C                                                                 -1
C             The leading N-by-N part of this array must contain A  ,
C                                                                 i
C             the inverse of the state transition matrix of the discrete
C             system at instant i.
C
C     LDAINV  INTEGER
C             The leading dimension of array AINV.  LDAINV >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain B ,
C                                                      -1         i
C             the input weight matrix (or the product A  B  if
C                                                      i  i
C             MULTAB = 'P') of the discrete system at instant i.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     RINV    (input) DOUBLE PRECISION array, dimension (LDRINV,*)
C             If MULTRC = 'N', then the leading P-by-P upper triangular
C                                              -1/2
C             part of this array must contain R    , the inverse of the
C                                              i+1
C             covariance square root (right Cholesky factor) of the
C             output (measurement) noise (hence the information square
C             root) at instant i+1.
C             The strict lower triangular part of this array is not
C             referenced.
C             Otherwise, RINV is not referenced and can be supplied as a
C             dummy array (i.e., set parameter LDRINV = 1 and declare
C             this array to be RINV(1,1) in the calling program).
C
C     LDRINV  INTEGER
C             The leading dimension of array RINV.
C             LDRINV >= MAX(1,P) if MULTRC = 'N';
C             LDRINV >= 1        if MULTRC = 'P'.
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain C   ,
C                                                       -1/2      i+1
C             the output weight matrix (or the product R    C    if
C                                                       i+1  i+1
C             MULTRC = 'P') of the discrete system at instant i+1.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     QINV    (input/output) DOUBLE PRECISION array, dimension
C             (LDQINV,M)
C             On entry, the leading M-by-M upper triangular part of this
C                                 -1/2
C             array must contain Q    , the inverse of the covariance
C                                 i
C             square root (right Cholesky factor) of the input (process)
C             noise (hence the information square root) at instant i.
C             On exit, the leading M-by-M upper triangular part of this
C                                    -1/2
C             array contains (QINOV )    , the inverse of the covariance
C                                  i
C             square root (right Cholesky factor) of the process noise
C             innovation (hence the information square root) at
C             instant i.
C             The strict lower triangular part of this array is not
C             referenced.
C
C     LDQINV  INTEGER
C             The leading dimension of array QINV.  LDQINV >= MAX(1,M).
C
C     X       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain X , the estimated
C                                                i
C             filtered state at instant i.
C             On exit, if JOBX = 'X', and INFO = 0, then this array
C             contains X   , the estimated filtered state at
C                       i+1
C             instant i+1.
C             On exit, if JOBX = 'N', or JOBX = 'X' and INFO = 1, then
C                                  -1
C             this array contains S   X   .
C                                  i+1 i+1
C
C     RINVY   (input) DOUBLE PRECISION array, dimension (P)
C                                      -1/2
C             This array must contain R    Y   , the product of the
C                                      i+1  i+1
C                                      -1/2
C             upper triangular matrix R     and the measured output
C                                      i+1
C             vector Y    at instant i+1.
C                     i+1
C
C     Z       (input) DOUBLE PRECISION array, dimension (M)
C             This array must contain Z , the mean value of the state
C                                      i
C             process noise at instant i.
C
C     E       (output) DOUBLE PRECISION array, dimension (P)
C             This array contains E   , the estimated error at instant
C                                  i+1
C             i+1.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If JOBX = 'X', then TOL is used to test for near
C                                        -1
C             singularity of the matrix S   . If the user sets
C                                        i+1
C             TOL > 0, then the given value of TOL is used as a
C             lower bound for the reciprocal condition number of that
C             matrix; a matrix whose estimated condition number is less
C             than 1/TOL is considered to be nonsingular. If the user
C             sets TOL <= 0, then an implicitly computed, default
C             tolerance, defined by TOLDEF = N*N*EPS, is used instead,
C             where EPS is the machine precision (see LAPACK Library
C             routine DLAMCH).
C             Otherwise, TOL is not referenced.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             where LIWORK = N if JOBX = 'X',
C             and   LIWORK = 1 otherwise.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.  If INFO = 0 and JOBX = 'X', DWORK(2) returns
C             an estimate of the reciprocal of the condition number
C                                 -1
C             (in the 1-norm) of S   .
C                                 i+1
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= MAX(1,N*(N+2*M)+3*M,(N+P)*(N+1)+2*N),
C                           if JOBX = 'N';
C             LDWORK >= MAX(2,N*(N+2*M)+3*M,(N+P)*(N+1)+2*N,3*N),
C                           if JOBX = 'X'.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;                        -1
C             = 1:  if JOBX = 'X' and the matrix S   is singular,
C                                                 i+1       -1
C                   i.e., the condition number estimate of S    (in the
C                                                           i+1
C                                                         -1    -1/2
C                   1-norm) exceeds 1/TOL.  The matrices S   , Q
C                                                         i+1   i
C                   and E have been computed.
C
C     METHOD
C
C     The routine performs one recursion of the square root information
C     filter algorithm, summarized as follows:
C
C       |    -1/2             -1/2    |     |         -1/2             |
C       |   Q         0      Q    Z   |     | (QINOV )     *     *     |
C       |    i                i    i  |     |       i                  |
C       |                             |     |                          |
C       |  -1 -1     -1 -1    -1      |     |             -1    -1     |
C     T | S  A  B   S  A     S  X     |  =  |    0       S     S   X   |
C       |  i  i  i   i  i     i  i    |     |             i+1   i+1 i+1|
C       |                             |     |                          |
C       |           -1/2      -1/2    |     |                          |
C       |    0     R    C    R    Y   |     |    0         0     E     |
C       |           i+1  i+1  i+1  i+1|     |                     i+1  |
C
C                  (Pre-array)                      (Post-array)
C
C     where T is an orthogonal transformation triangularizing the
C                        -1/2
C     pre-array, (QINOV )     is the inverse of the covariance square
C                      i
C     root (right Cholesky factor) of the process noise innovation
C     (hence the information square root) at instant i, and E    is the
C                                                            i+1
C     estimated error at instant i+1.
C
C     The inverse of the corresponding state covariance matrix P
C                                                               i+1|i+1
C     (hence the information matrix I) is then factorized as
C
!C                   -1         -1     -1
!C        I       = P       = (S   )' S
!C         i+1|i+1   i+1|i+1    i+1    i+1
!C
!C     and one combined time and measurement update for the state is
!C     given by X   .
!C               i+1
!C
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
C     [3] Vanbegin, M., Van Dooren, P., and Verhaegen, M.H.G.
C         Algorithm 675: FORTRAN Subroutines for Computing the Square
C         Root Covariance Filter and Square Root Information Filter in
C         Dense or Hessenberg Forms.
C         ACM Trans. Math. Software, 15, pp. 243-256, 1989.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires approximately
C
C           3    2                              2   2
C     (7/6)N  + N x (7/2 x M + P) + N x (1/2 x P + M )
C
C     operations and is backward stable (see [2]).
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C     Supersedes Release 2.0 routine FB01GD by M. Vanbegin,
C     P. Van Dooren, and M.H.G. Verhaegen.
C
C     REVISIONS
C
C     February 20, 1998, November 20, 2003, February 14, 2004.
C
C     KEYWORDS
C
C     Kalman filtering, optimal filtering, orthogonal transformation,
C     recursive estimation, square-root filtering, square-root
C     information filtering.
C
C     ******************************************************************
C
#endif
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
      use omp_lib
#endif
      implicit none

!C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!C     .. Scalar Arguments ..
      CHARACTER         JOBX, MULTAB, MULTRC
      INTEGER           INFO, LDAINV, LDB, LDC, LDQINV, LDRINV, LDSINV, &
                        LDWORK, M, N, P
      DOUBLE PRECISION  TOL
!C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      INTEGER           IWORK(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER           IWORK(*)
      !DIR$ ASSUME_ALIGNED IWORK:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION  AINV(LDAINV,*), B(LDB,*), C(LDC,*), DWORK(*), &
                        E(*), QINV(LDQINV,*), RINV(LDRINV,*), RINVY(*), &
                        SINV(LDSINV,*), X(*), Z(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
       DOUBLE PRECISION  AINV(LDAINV,*), B(LDB,*), C(LDC,*), DWORK(*), &
                        E(*), QINV(LDQINV,*), RINV(LDRINV,*), RINVY(*), &
                        SINV(LDSINV,*), X(*), Z(*)
      !DIR$ ASSUME_ALIGNED AINV:64
      !DIR$ ASSUME_ALIGNED B:64
      !DIR$ ASSUME_ALIGNED C:64
      !DIR$ ASSUME_ALIGNED DWORK:64
      !DIR$ ASSUME_ALIGNED E:64
      !DIR$ ASSUME_ALIGNED QINV:64
      !DIR$ ASSUME_ALIGNED RINV:64
      !DIR$ ASSUME_ALIGNED RINVY:64
      !DIR$ ASSUME_ALIGNED SINV:64
      !DIR$ ASSUME_ALIGNED X:64
      !DIR$ ASSUME_ALIGNED Z:64
#endif
!C     .. Local Scalars ..
      LOGICAL           LJOBX, LMULTA, LMULTR
      INTEGER           I, I12, I13, I21, I23, IJ, ITAU, JWORK, LDW, M1, &
                        N1, NP, WRKOPT
      DOUBLE PRECISION  RCOND
!C     .. External Functions ..
     
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
!C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMM, DGEQRF, DLACPY, DORMQR,
     $                  DTRMM, DTRMV
!C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX
!C     .. Executable Statements ..
!C
      NP = N + P
      N1 = MAX( 1, N )
      M1 = MAX( 1, M )
      INFO = 0
      LJOBX  = LSAME( JOBX, 'X' )
      LMULTA = LSAME( MULTAB, 'P' )
      LMULTR = LSAME( MULTRC, 'P' )
!C
!C     Test the input scalar arguments.
!C
      IF( .NOT.LJOBX .AND. .NOT.LSAME( JOBX, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LMULTA .AND. .NOT.LSAME( MULTAB, 'N' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LMULTR .AND. .NOT.LSAME( MULTRC, 'N' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDSINV.LT.N1 ) THEN
         INFO = -8
      ELSE IF( LDAINV.LT.N1 ) THEN
         INFO = -10
      ELSE IF( LDB.LT.N1 ) THEN
         INFO = -12
      ELSE IF( LDRINV.LT.1 .OR. ( .NOT.LMULTR .AND. LDRINV.LT.P ) ) THEN
         INFO = -14
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -16
      ELSE IF( LDQINV.LT.M1 ) THEN
         INFO = -18
      ELSE IF( ( LJOBX .AND. LDWORK.LT.MAX( 2, N*(N + 2*M) + 3*M,  &
                                           NP*(N + 1) + 2*N, 3*N ) ) &
                                                                  .OR. &
          ( .NOT.LJOBX .AND. LDWORK.LT.MAX( 1, N*(N + 2*M) + 3*M,  &
                                            NP*(N + 1) + 2*N ) ) ) THEN &
         INFO = -26
      END IF

      IF ( INFO.NE.0 ) THEN
!C
!C        Error return.
!C
          RETURN
      END IF
!C
!C     Quick return if possible.
!C
      IF ( MAX( N, P ).EQ.0 ) THEN
         IF ( LJOBX ) THEN
            DWORK(1) = TWO
            DWORK(2) = ONE
         ELSE
            DWORK(1) = ONE
         END IF
         RETURN
      END IF
!C
!C     Construction of the needed part of the pre-array in DWORK.
!C     To save workspace, only the blocks (1,3), (2,1)-(2,3), (3,2), and
!C     (3,3) will be constructed when needed as shown below.
!C
!C     Storing SINV x AINV and SINV x AINV x B in the (1,1) and (1,2)
!C     blocks of DWORK, respectively.
!C     The variables called Ixy define the starting positions where the
!C     (x,y) blocks of the pre-array are initially stored in DWORK.
!C     Workspace: need N*(N+M).
!C
!C     (Note: Comments in the code beginning "Workspace:" describe the
!C     minimal amount of real workspace needed at that point in the
!C     code, as well as the preferred amount for good performance.
!C     NB refers to the optimal block size for the immediately
!C     following subroutine, as returned by ILAENV.)
!C
      LDW = N1
      I21 = N*N + 1

      CALL DLACPY( 'Full', N, N, AINV, LDAINV, DWORK, LDW )
      IF ( LMULTA ) THEN
         CALL DLACPY( 'Full', N, M, B, LDB, DWORK(I21), LDW )
      ELSE
         CALL DGEMM( 'No transpose', 'No transpose', N, M, N, ONE,  &
                     DWORK, LDW, B, LDB, ZERO, DWORK(I21), LDW )
      END IF
      CALL DTRMM(  'Left', 'Upper', 'No transpose', 'Non-unit', N, N+M,  &
                   ONE, SINV, LDSINV, DWORK, LDW )
!C
!C     Storing the process noise mean value in (1,3) block of DWORK.
!C     Workspace: need N*(N+M) + M.
!C
      I13 = N*( N + M ) + 1
!C
      CALL DCOPY( M, Z, 1, DWORK(I13), 1 )
      CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', M, QINV, LDQINV, &
                 DWORK(I13), 1 )
!C
!C     Computing SINV x X in X.
!C
      CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', N, SINV, LDSINV,  &
                  X, 1 )
!C
!C     Triangularization (2 steps).
!C
!C     Step 1: annihilate the matrix SINV x AINV x B.
!C     Workspace: need N*(N+2*M) + 3*M.
!C
      I12   = I13  + M
      ITAU  = I12  + M*N
      JWORK = ITAU + M
!C
      CALL MB04KD( 'Full', M, N, N, QINV, LDQINV, DWORK(I21), LDW,  &
                  DWORK, LDW, DWORK(I12), M1, DWORK(ITAU), &
                  DWORK(JWORK) )
      WRKOPT = MAX( 1, N*( N + 2*M ) + 3*M )

      IF ( N.EQ.0 ) THEN
         CALL DCOPY( P, RINVY, 1, E, 1 )
         IF ( LJOBX )
     $      DWORK(2) = ONE
         DWORK(1) = WRKOPT
         RETURN
      END IF
!C
!C     Apply the transformations to the last column of the pre-array.
!C     (Only the updated (2,3) block is now needed.)
!C
      IJ = I21
!C
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFULT(SHARED) PRIVATE(I)
#endif
      DO 10 I = 1, M
         CALL DAXPY( N, -DWORK(ITAU+I-1)*( DWORK(I13+I-1) + &
                    DDOT( N, DWORK(IJ), 1, X, 1 ) ), &
                              DWORK(IJ), 1, X, 1 )
         IJ = IJ + N
10       CONTINUE
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1         
!$OMP END PARALLEL DO
#endif
!C
!C     Now, the workspace for SINV x AINV x B, as well as for the updated
!C     (1,2) block of the pre-array, are no longer needed.
!C     Move the computed (2,3) block of the pre-array in the (1,2) block
!C     position of DWORK, to save space for the following computations.
!C     Then, adjust the implicitly defined leading dimension of DWORK,
!C     to make space for storing the (3,2) and (3,3) blocks of the
!C     pre-array.
!C     Workspace: need (N+P)*(N+1).
!C
      CALL DCOPY( N, X, 1, DWORK(I21), 1 )
      LDW = MAX( 1, NP )
!C
      DO 30 I = N + 1, 1, -1
         DO 20 IJ = N, 1, -1
            DWORK(NP*(I-1)+IJ) = DWORK(N*(I-1)+IJ)
   20    CONTINUE
   30 CONTINUE
!C
!C     Copy of RINV x C in the (2,1) block of DWORK.
!C
      CALL DLACPY( 'Full', P, N, C, LDC, DWORK(N+1), LDW )
      IF ( .NOT.LMULTR ) &
         CALL DTRMM(  'Left', 'Upper', 'No transpose', 'Non-unit', P, N, &
                     ONE, RINV, LDRINV, DWORK(N+1), LDW )
!C
!C     Copy the inclusion measurement in the (2,2) block of DWORK.
!C
      I21 = NP*N + 1
      I23 = I21  + N
      CALL DCOPY( P, RINVY, 1, DWORK(I23), 1 )
      WRKOPT = MAX( WRKOPT, NP*( N + 1 ) )
!C
!C     Step 2: QR factorization of the first block column of the matrix
!C
!C        [ SINV x AINV  SINV x X ]
!C        [ RINV x C     RINV x Y ],
!C
!C     where the first block row was modified at Step 1.
!C     Workspace: need   (N+P)*(N+1) + 2*N;
!C                prefer (N+P)*(N+1) + N + N*NB.
!C
      ITAU  = I21  + NP
      JWORK = ITAU + N
!C
      CALL DGEQRF( NP, N, DWORK, LDW, DWORK(ITAU), DWORK(JWORK), &
                   LDWORK-JWORK+1, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
!C
!C     Apply the Householder transformations to the last column.
!C     Workspace: need (N+P)*(N+1) + 1;  prefer (N+P)*(N+1) + NB.
!C
      CALL DORMQR( 'Left', 'Transpose', NP, 1, N, DWORK, LDW,  &
                  DWORK(ITAU), DWORK(I21), LDW, DWORK(JWORK),  &
                   LDWORK-JWORK+1, INFO )  
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
!C
!C     Output SINV, X, and E and set the optimal workspace dimension
!C     (and the reciprocal of the condition number estimate).
!C
      CALL DLACPY( 'Upper', N, N, DWORK, LDW, SINV, LDSINV )
      CALL DCOPY( N, DWORK(I21), 1, X, 1 )
      CALL DCOPY( P, DWORK(I23), 1, E, 1 )
!C
      IF ( LJOBX ) THEN
!C
!C        Compute X.
!C        Workspace: need 3*N.
!C
         CALL MB02OD( 'Left', 'Upper', 'No transpose', 'Non-unit',   &
                     '1-norm', N, 1, ONE, SINV, LDSINV, X, N, RCOND, &
                     TOL, IWORK, DWORK, INFO )
         IF ( INFO.EQ.0 ) THEN
            WRKOPT = MAX( WRKOPT, 3*N )
            DWORK(2) = RCOND
         END IF
      END IF

      DWORK(1) = WRKOPT

END SUBROUTINE FB01SD

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE FB01VD( N, M, L, P, LDP, A, LDA, B, LDB, C, LDC, Q, &
  LDQ, R, LDR, K, LDK, TOL, IWORK, DWORK, LDWORK,INFO) !GCC$ ATTRIBUTES hot :: FB01VD !GCC$ ATTRIBUTES aligned(32) :: FB01VD !GCC$ ATTRIBUTES no_stack_protector :: FB01VD
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE FB01VD( N, M, L, P, LDP, A, LDA, B, LDB, C, LDC, Q, &
  LDQ, R, LDR, K, LDK, TOL, IWORK, DWORK, LDWORK,INFO)
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: FB01VD
  !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: FB01VD
#endif
#if 0                       
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To compute one recursion of the conventional Kalman filter
C     equations. This is one update of the Riccati difference equation
C     and the Kalman filter gain.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e., the order of the
C             matrices P      and A .  N >= 0.
C                       i|i-1      i
C
C     M       (input) INTEGER
C             The actual input dimension, i.e., the order of the matrix
C             Q .  M >= 0.
C              i
C
C     L       (input) INTEGER
C             The actual output dimension, i.e., the order of the matrix
C             R .  L >= 0.
C              i
C
C     P       (input/output) DOUBLE PRECISION array, dimension (LDP,N)
C             On entry, the leading N-by-N part of this array must
C             contain P     , the state covariance matrix at instant
C                      i|i-1
C             (i-1). The upper triangular part only is needed.
C             On exit, if INFO = 0, the leading N-by-N part of this
C             array contains P     , the state covariance matrix at
C                             i+1|i
C             instant i. The strictly lower triangular part is not set.
C             Otherwise, the leading N-by-N part of this array contains
C             P     , its input value.
C              i|i-1
C
C     LDP     INTEGER
C             The leading dimension of array P.  LDP >= MAX(1,N).
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
C                                                                 i
C             the input weight matrix of the discrete system at
C             instant i.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading L-by-N part of this array must contain C ,
C                                                                 i
C             the output weight matrix of the discrete system at
C             instant i.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,L).
C
C     Q       (input) DOUBLE PRECISION array, dimension (LDQ,M)
C             The leading M-by-M part of this array must contain Q ,
C                                                                 i
C             the input (process) noise covariance matrix at instant i.
C             The diagonal elements of this array are modified by the
C             routine, but are restored on exit.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.  LDQ >= MAX(1,M).
C
C     R       (input/output) DOUBLE PRECISION array, dimension (LDR,L)
C             On entry, the leading L-by-L part of this array must
C             contain R , the output (measurement) noise covariance
C                      i
C             matrix at instant i.
C             On exit, if INFO = 0, or INFO = L+1, the leading L-by-L
C                                                                  1/2
C             upper triangular part of this array contains (RINOV )   ,
C                                                                i
C             the square root (left Cholesky factor) of the covariance
C             matrix of the innovations at instant i.
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,L).
C
C     K       (output) DOUBLE PRECISION array, dimension (LDK,L)
C             If INFO = 0, the leading N-by-L part of this array
C             contains K , the Kalman filter gain matrix at instant i.
C                       i
!C             If INFO > 0, the leading N-by-L part of this array
!C             contains the matrix product P     C'.
C                                          i|i-1 i
C
C     LDK     INTEGER
C             The leading dimension of array K.  LDK >= MAX(1,N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used to test for near singularity of
C             the matrix RINOV . If the user sets TOL > 0, then the
C                             i
C             given value of TOL is used as a lower bound for the
C             reciprocal condition number of that matrix; a matrix whose
C             estimated condition number is less than 1/TOL is
C             considered to be nonsingular. If the user sets TOL <= 0,
C             then an implicitly computed, default tolerance, defined by
C             TOLDEF = L*L*EPS, is used instead, where EPS is the
C             machine precision (see LAPACK Library routine DLAMCH).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (L)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, or INFO = L+1, DWORK(1) returns an
C             estimate of the reciprocal of the condition number (in the
C             1-norm) of the matrix RINOV .
C                                        i
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= MAX(1,L*N+3*L,N*N,N*M).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -k, the k-th argument had an illegal
C                   value;
C             = k:  if INFO = k, 1 <= k <= L, the leading minor of order
C                   k of the matrix RINOV  is not positive-definite, and
C                                        i
C                   its Cholesky factorization could not be completed;
C             = L+1: the matrix RINOV  is singular, i.e., the condition
C                                    i
C                   number estimate of RINOV  (in the 1-norm) exceeds
C                                           i
C                   1/TOL.
C
C     METHOD
C
C     The conventional Kalman filter gain used at the i-th recursion
C     step is of the form
C
C                            -1
C        K  = P     C'  RINOV  ,
!C         i    i|i-1 i       i
!C
!C     where RINOV  = C P     C' + R , and the state covariance matrix
C                i    i i|i-1 i    i
C
C     P      is updated by the discrete-time difference Riccati equation
C      i|i-1
C
C        P      = A  (P      - K C P     ) A'  + B Q B'.
C         i+1|i    i   i|i-1    i i i|i-1   i     i i i
C
C     Using these two updates, the combined time and measurement update
C     of the state X      is given by
C                   i|i-1
C
C        X      = A X      + A K (Y  - C X     ),
C         i+1|i    i i|i-1    i i  i    i i|i-1
C
C     where Y  is the new observation at step i.
C            i
C
C     REFERENCES
C
C     [1] Anderson, B.D.O. and Moore, J.B.
C         Optimal Filtering,
C         Prentice Hall, Englewood Cliffs, New Jersey, 1979.
C
C     [2] Verhaegen, M.H.G. and Van Dooren, P.
C         Numerical Aspects of Different Kalman Filter Implementations.
C         IEEE Trans. Auto. Contr., AC-31, pp. 907-917, 1986.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires approximately
C
C             3   2
C      3/2 x N + N  x (3 x L + M/2)
C
C     operations.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C     Supersedes Release 2.0 routine FB01JD by M.H.G. Verhaegen,
C     M. Vanbegin, and P. Van Dooren.
C
C     REVISIONS
C
C     February 20, 1998, November 20, 2003, April 20, 2004.
C
C     KEYWORDS
C
C     Kalman filtering, optimal filtering, recursive estimation.
C
C     ******************************************************************
C
#endif
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
      use omp_lib
#endif
      implicit none
!C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!C     .. Scalar Arguments ..
      INTEGER           INFO, L, LDA, LDB, LDC, LDK, LDP, LDQ, LDR, &
                        LDWORK, M, N
      DOUBLE PRECISION  TOL
      !C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      INTEGER           IWORK(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER           IWORK(*)
      !DIR$ ASSUME_ALIGNED IWORK:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           K(LDK,*), P(LDP,*), Q(LDQ,*), R(LDR,*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           K(LDK,*), P(LDP,*), Q(LDQ,*), R(LDR,*)
      !DIR$ ASSUME_ALIGNED A:64
      !DIR$ ASSUME_ALIGNED B:64
      !DIR$ ASSUME_ALIGNED C:64
      !DIR$ ASSUME_ALIGNED DWORK:64
      !DIR$ ASSUME_ALIGNED K:64
      !DIR$ ASSUME_ALIGNED P:64
      !DIR$ ASSUME_ALIGNED Q:64
      !DIR$ ASSUME_ALIGNED R:64
#endif
!C     .. Local Scalars ..
      INTEGER           J, JWORK, LDW, N1
      DOUBLE PRECISION  RCOND, RNORM, TOLDEF
!C     .. External Functions ..
      DOUBLE PRECISION   DLANSY
      EXTERNAL           DLANSY
!C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DLACPY, DLASET, DPOCON, &
                        DPOTRF, DSCAL, DTRMM, DTRSM, MB01RD
!C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX
!C     .. Executable Statements ..
!C
!C     Test the input scalar arguments.
!C
      INFO = 0
      N1 = MAX( 1, N )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( L.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDP.LT.N1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.N1 ) THEN
         INFO = -7
      ELSE IF( LDB.LT.N1 ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, L ) ) THEN
         INFO = -11
      ELSE IF( LDQ.LT.MAX( 1, M ) ) THEN
         INFO = -13
      ELSE IF( LDR.LT.MAX( 1, L ) ) THEN
         INFO = -15
      ELSE IF( LDK.LT.N1 ) THEN
         INFO = -17
      ELSE IF( LDWORK.LT.MAX( 1, L*N + 3*L, N*N, N*M ) ) THEN
         INFO = -21
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
      IF ( MAX( N, L ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
!C
!C     Efficiently compute RINOV = CPC' + R in R and put CP in DWORK and
!C     PC' in K. (The content of DWORK on exit from MB01RD is used.)
!C     Workspace: need L*N.
!C
!C     (Note: Comments in the code beginning "Workspace:" describe the
!C     minimal amount of real workspace needed at that point in the
!C     code.)
!C
      CALL MB01RD( 'Upper', 'No transpose', L, N, ONE, ONE, R, LDR, C, &
                   LDC, P, LDP, DWORK, LDWORK, INFO )
      LDW = MAX( 1, L )
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(J)
#endif
      DO 10 J = 1, L
         CALL DCOPY( N, DWORK(J), LDW, K(1,J), 1 )
   10 CONTINUE
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
         !$OMP END PARALLEL DO
#endif
      CALL DLACPY( 'Full', L, N, C, LDC, DWORK, LDW )
      CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-unit', L, N, ONE, &
                  P, LDP, DWORK, LDW )
      CALL DSCAL( N, TWO, P, LDP+1 )
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(J)
#endif
      DO 20 J = 1, L
         CALL DAXPY( N, ONE, K(1,J), 1, DWORK(J), LDW )
         CALL DCOPY( N, DWORK(J), LDW, K(1,J), 1 )
20       CONTINUE
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
         !$OMP END PARALLEL DO
#endif
!C
!C     Calculate the Cholesky decomposition U'U of the innovation
!C     covariance matrix RINOV, and its reciprocal condition number.
!C     Workspace: need L*N + 3*L.
!C
      JWORK = L*N + 1
      RNORM = DLANSY( '1-norm', 'Upper', L, R, LDR, DWORK(JWORK) )
!C
      TOLDEF = TOL
      IF ( TOLDEF.LE.ZERO ) &
         TOLDEF = DBLE( L*L )*DLAMCH( 'Epsilon' )
      CALL DPOTRF( 'Upper', L, R, LDR, INFO )
      IF ( INFO.NE.0 ) RETURN
     
!C
      CALL DPOCON( 'Upper', L, R, LDR, RNORM, RCOND, DWORK(JWORK), &
                  IWORK, INFO )
!C
      IF ( RCOND.LT.TOLDEF ) THEN
!C
!C        Error return: RINOV is numerically singular.
!C
         INFO = L+1
         DWORK(1) = RCOND
         RETURN
      END IF

      IF ( L.GT.1 ) &
        CALL DLASET( 'Lower', L-1, L-1, ZERO, ZERO, R(2,1),LDR )
!C                                                          -1
!C     Calculate the Kalman filter gain matrix  K = PC'RINOV .
!C     Workspace: need L*N.
!C
      CALL DTRSM( 'Right', 'Upper', 'No transpose', 'Non-unit', N, L, &
                  ONE, R, LDR, K, LDK )
      CALL DTRSM( 'Right', 'Upper', 'Transpose', 'Non-unit', N, L,  &
                  ONE, R, LDR, K, LDK )
!C
!C     First part of the Riccati equation update: compute A(P-KCP)A'.
!C     The upper triangular part of the symmetric matrix P-KCP is formed.
!C     Workspace: need max(L*N,N*N).
!C
      JWORK = 1
      !C
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(J,JWORK)
#endif
      DO 30 J = 1, N
         CALL DGEMV( 'No transpose', J, L, -ONE, K, LDK, DWORK(JWORK), &
                     1, ONE, P(1,J), 1 )
         JWORK = JWORK + L
30       CONTINUE
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
         !$OMP END PARALLEL DO
#endif
!C
      CALL MB01RD( 'Upper', 'No transpose', N, N, ZERO, ONE, P, LDP, A, &
                   LDA, P, LDP, DWORK, LDWORK, INFO )
!C
!C     Second part of the Riccati equation update: add BQB'.
!C     Workspace: need N*M.
!C
      CALL MB01RD( 'Upper', 'No transpose', N, M, ONE, ONE, P, LDP, B, &
                   LDB, Q, LDQ, DWORK, LDWORK, INFO )
      CALL DSCAL( M, TWO, Q, LDQ+1 )
!C
!C     Set the reciprocal of the condition number estimate.
!C
      DWORK(1) = RCOND
!C
   
END SUBROUTINE FB01VD


    

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE MB02OD( SIDE, UPLO, TRANS, DIAG, NORM, M, N, ALPHA, A,
  LDA, B, LDB, RCOND, TOL, IWORK, DWORK, INFO ) !GCC$ ATTRIBUTES HOT :: MB02OD !GCC$ ATTRIBUTES aligned(32) :: MB020D !GCC$ ATTRIBUTES no_stack_protector :: MB020D
#elif defined(__ICC) || defined(__INTEL_COMPILER)
 SUBROUTINE MB02OD( SIDE, UPLO, TRANS, DIAG, NORM, M, N, ALPHA, A,
   LDA, B, LDB, RCOND, TOL, IWORK, DWORK, INFO )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: NB020D
   !DIR$ OPTIMIZE : 2
   !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: MB02OD
#endif
#if 0 
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To solve (if well-conditioned) one of the matrix equations
C
C        op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
C
C     where alpha is a scalar, X and B are m-by-n matrices, A is a unit,
C     or non-unit, upper or lower triangular matrix and op( A ) is one
C     of
!C
!C        op( A ) = A   or   op( A ) = A'.
C
C     An estimate of the reciprocal of the condition number of the
C     triangular matrix A, in either the 1-norm or the infinity-norm, is
C     also computed as
C
C        RCOND = 1 / ( norm(A) * norm(inv(A)) ).
C
C     and the specified matrix equation is solved only if RCOND is
C     larger than a given tolerance TOL.  In that case, the matrix X is
C     overwritten on B.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     SIDE    CHARACTER*1
C             Specifies whether op( A ) appears on the left or right
C             of X as follows:
C             = 'L':  op( A )*X = alpha*B;
C             = 'R':  X*op( A ) = alpha*B.
C
C     UPLO    CHARACTER*1
C             Specifies whether the matrix A is an upper or lower
C             triangular matrix as follows:
C             = 'U':  A is an upper triangular matrix;
C             = 'L':  A is a lower triangular matrix.
C
C     TRANS   CHARACTER*1
C             Specifies the form of op( A ) to be used in the matrix
C             multiplication as follows:
C             = 'N':  op( A ) = A;
C             = 'T':  op( A ) = A';
C             = 'C':  op( A ) = A'.
C
C     DIAG    CHARACTER*1
C             Specifies whether or not A is unit triangular as follows:
C             = 'U':  A is assumed to be unit triangular;
C             = 'N':  A is not assumed to be unit triangular.
C
C     NORM    CHARACTER*1
C             Specifies whether the 1-norm condition number or the
C             infinity-norm condition number is required:
C             = '1' or 'O':  1-norm;
C             = 'I':         Infinity-norm.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of B.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of B.  N >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar  alpha. When alpha is zero then A is not
C             referenced and B need not be set before entry.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,k),
C             where k is M when SIDE = 'L' and is N when SIDE = 'R'.
C             On entry with UPLO = 'U', the leading k-by-k upper
C             triangular part of this array must contain the upper
C             triangular matrix and the strictly lower triangular part
C             of A is not referenced.
C             On entry with UPLO = 'L', the leading k-by-k lower
!C             triangular part of this array must contain the lower
C             triangular matrix and the strictly upper triangular part
C             of A is not referenced.
C             Note that when DIAG = 'U', the diagonal elements of A are
C             not referenced either, but are assumed to be unity.
!C
!C     LDA     INTEGER
!C             The leading dimension of array A.
!C             LDA >= max(1,M) when SIDE = 'L';
C             LDA >= max(1,N) when SIDE = 'R'.
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading M-by-N part of this array must
C             contain the right-hand side matrix B.
C             On exit, if INFO = 0, the leading M-by-N part of this
C             array contains the solution matrix X.
C             Otherwise, this array is not modified by the routine.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= max(1,M).
C
C     RCOND   (output) DOUBLE PRECISION
C             The reciprocal of the condition number of the matrix A,
C             computed as RCOND = 1/(norm(A) * norm(inv(A))).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used to test for near singularity of
C             the matrix A. If the user sets TOL > 0, then the given
C             value of TOL is used as a lower bound for the reciprocal
C             condition number of that matrix; a matrix whose estimated
!C             condition number is less than 1/TOL is considered to be
!C             nonsingular. If the user sets TOL <= 0, then an implicitly
!C             computed, default tolerance, defined by TOLDEF = k*k*EPS,
!C             is used instead, where EPS is the machine precision (see
!C             LAPACK Library routine DLAMCH).
!C
!C     Workspace
!C
!C     IWORK   INTEGER array, dimension (k)
C
C     DWORK   DOUBLE PRECISION array, dimension (3*k)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
!C             = 1:  the matrix A is numerically singular, i.e. the
!C                   condition number estimate of A (in the specified
!C                   norm) exceeds 1/TOL.
!C
!C     METHOD
!C
!C     An estimate of the reciprocal of the condition number of the
!C     triangular matrix A (in the specified norm) is computed, and if
!C     this estimate is larger then the given (or default) tolerance,
!C     the specified matrix equation is solved using Level 3 BLAS
!C     routine DTRSM.
C
C
C     REFERENCES
C
C     None.
C
C     NUMERICAL ASPECTS
C                             2
C     The algorithm requires k N/2 operations.
C
C     CONTRIBUTORS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C
C     REVISIONS
C
C     February 20, 1998.
C
C     KEYWORDS
C
C     Condition number, matrix algebra, matrix operations.
C
C     ******************************************************************
C
#endif
       implicit none
!C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
!C     .. Scalar Arguments ..
      CHARACTER          DIAG, NORM, SIDE, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, M, N
      DOUBLE PRECISION   ALPHA, RCOND, TOL
!C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      INTEGER            IWORK(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER IWORK(*)
!DIR$ ASSUME_ALIGNED IWORK:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), DWORK(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), DWORK(*)
!DIR$ ASSUME_ALIGNED A:64
!DIR$ ASSUME_ALIGNED B:64
!DIR$ ASSUME_ALIGNED DWORK:64
#endif
!C     .. Local Scalars ..
      LOGICAL            LSIDE, ONENRM
      INTEGER            NROWA
      DOUBLE PRECISION   TOLDEF


!C     .. External Subroutines ..
      EXTERNAL           DTRCON, DTRSM
!C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
!C     .. Executable Statements ..
!C
      LSIDE  = LSAME( SIDE, 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
!C
!C     Test the input scalar arguments.
!C
      INFO = 0
      IF( ( .NOT.LSIDE ).AND.( .NOT.LSAME( SIDE, 'R' ) ) )THEN
         INFO = -1
      ELSE IF( ( .NOT.LSAME( UPLO,  'U' ) ).AND.  &
               ( .NOT.LSAME( UPLO,  'L' ) )      )THEN
         INFO = -2
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND. &
               ( .NOT.LSAME( TRANS, 'T' ) ).AND.
               ( .NOT.LSAME( TRANS, 'C' ) )      )THEN &
         INFO = -3
      ELSE IF( ( .NOT.LSAME( DIAG,  'U' ) ).AND.  &
              ( .NOT.LSAME( DIAG,  'N' ) )      )THEN
         INFO = -4
      ELSE IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 )THEN
         INFO = -6
      ELSE IF( N.LT.0 )THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, M ) )THEN
         INFO = -12
      END IF
!C
      IF( INFO.NE.0 ) THEN
          RETURN
      END IF
!C
!C     Quick return if possible.
!C
      IF( NROWA.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      END IF

      TOLDEF = TOL
      IF ( TOLDEF.LE.ZERO ) &
         TOLDEF = DBLE( NROWA*NROWA )*DLAMCH( 'Epsilon' )

      CALL DTRCON( NORM, UPLO, DIAG, NROWA, A, LDA, RCOND, DWORK, &
                   IWORK, INFO )

      IF ( RCOND.GT.TOLDEF ) THEN
         CALL DTRSM( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, &
                     LDB )
      ELSE
         INFO = 1
      END IF

END SUBROUTINE FB01SD

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE FB01TD( JOBX, MULTRC, N, M, P, SINV, LDSINV, AINV, &
        LDAINV, AINVB, LDAINB, RINV, LDRINV, C, LDC,          &
        QINV, LDQINV, X, RINVY, Z, E, TOL, IWORK,             &
        DWORK, LDWORK, INFO ) !GCC$ ATTRIBUTES HOT :: FB01TD !GCC$ ATTRIBUTES aligned(32) :: FB01TD !GCC$ ATTRIBUTES no_stack_protector :: FB01TD
#elif defined(__ICC) || defined(__INTEL_COMPILER)
 SUBROUTINE FB01TD( JOBX, MULTRC, N, M, P, SINV, LDSINV, AINV, &
        LDAINV, AINVB, LDAINB, RINV, LDRINV, C, LDC,          &
        QINV, LDQINV, X, RINVY, Z, E, TOL, IWORK,             &
        DWORK, LDWORK, INFO )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: FB01TD
   !DIR$ OPTIMIZE : 3
   !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: FB01TD
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
C     iteration of the time-invariant Kalman filter. This update is
C     given for the square root information filter, using the condensed
C     controller Hessenberg form.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBX    CHARACTER*1
C             Indicates whether X    is to be computed as follows:
C                                i+1
C             = 'X':  X    is computed and stored in array X;
C                      i+1
C             = 'N':  X    is not required.
C                      i+1
C
C     MULTRC  CHARACTER*1             -1/2
C             Indicates how matrices R     and C    are to be passed to
C                                     i+1       i+1
C             the routine as follows:
C             = 'P':  Array RINV is not used and the array C must
C                                          -1/2
C                     contain the product R    C   ;
C                                          i+1  i+1
C             = 'N':  Arrays RINV and C must contain the matrices
C                     as described below.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e., the order of the
C                       -1      -1
C             matrices S   and A  .  N >= 0.
C                       i
C
C     M       (input) INTEGER
C             The actual input dimension, i.e., the order of the matrix
C              -1/2
C             Q    .  M >= 0.
C              i
C
C     P       (input) INTEGER
C             The actual output dimension, i.e., the order of the matrix
C              -1/2
C             R    .  P >= 0.
C              i+1
C
C     SINV    (input/output) DOUBLE PRECISION array, dimension
C             (LDSINV,N)
C             On entry, the leading N-by-N upper triangular part of this
C                                 -1
C             array must contain S  , the inverse of the square root
C                                 i
C             (right Cholesky factor) of the state covariance matrix
C             P    (hence the information square root) at instant i.
C              i|i
C             On exit, the leading N-by-N upper triangular part of this
C                             -1
C             array contains S   , the inverse of the square root (right
C                             i+1
C             Cholesky factor) of the state covariance matrix P
C                                                              i+1|i+1
C             (hence the information square root) at instant i+1.
C             The strict lower triangular part of this array is not
C             referenced.
C
C     LDSINV  INTEGER
C             The leading dimension of array SINV.  LDSINV >= MAX(1,N).
C
C     AINV    (input) DOUBLE PRECISION array, dimension (LDAINV,N)
C                                                                 -1
C             The leading N-by-N part of this array must contain A  ,
C             the inverse of the state transition matrix of the discrete
C             system in controller Hessenberg form (e.g., as produced by
C             SLICOT Library Routine TB01MD).
C
C     LDAINV  INTEGER
C             The leading dimension of array AINV.  LDAINV >= MAX(1,N).
C
C     AINVB   (input) DOUBLE PRECISION array, dimension (LDAINB,M)
C                                                                  -1
C             The leading N-by-M part of this array must contain  A  B,
C                             -1
C             the product of A   and the input weight matrix B of the
C             discrete system, in upper controller Hessenberg form
C             (e.g., as produced by SLICOT Library Routine TB01MD).
C
C     LDAINB  INTEGER
C             The leading dimension of array AINVB.  LDAINB >= MAX(1,N).
C
C     RINV    (input) DOUBLE PRECISION array, dimension (LDRINV,*)
C             If MULTRC = 'N', then the leading P-by-P upper triangular
C                                              -1/2
C             part of this array must contain R    , the inverse of the
C                                              i+1
C             covariance square root (right Cholesky factor) of the
C             output (measurement) noise (hence the information square
C             root) at instant i+1.
C             The strict lower triangular part of this array is not
C             referenced.
C             Otherwise, RINV is not referenced and can be supplied as a
C             dummy array (i.e., set parameter LDRINV = 1 and declare
C             this array to be RINV(1,1) in the calling program).
C
C     LDRINV  INTEGER
C             The leading dimension of array RINV.
C             LDRINV >= MAX(1,P) if MULTRC = 'N';
C             LDRINV >= 1        if MULTRC = 'P'.
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain C   ,
C                                                       -1/2      i+1
C             the output weight matrix (or the product R    C    if
C                                                       i+1  i+1
C             MULTRC = 'P') of the discrete system at instant i+1.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     QINV    (input/output) DOUBLE PRECISION array, dimension
C             (LDQINV,M)
C             On entry, the leading M-by-M upper triangular part of this
C                                 -1/2
C             array must contain Q    , the inverse of the covariance
C                                 i
C             square root (right Cholesky factor) of the input (process)
C             noise (hence the information square root) at instant i.
C             On exit, the leading M-by-M upper triangular part of this
C                                    -1/2
C             array contains (QINOV )    , the inverse of the covariance
C                                  i
C             square root (right Cholesky factor) of the process noise
C             innovation (hence the information square root) at
C             instant i.
C             The strict lower triangular part of this array is not
C             referenced.
C
C     LDQINV  INTEGER
C             The leading dimension of array QINV.  LDQINV >= MAX(1,M).
C
C     X       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain X , the estimated
C                                                i
C             filtered state at instant i.
C             On exit, if JOBX = 'X', and INFO = 0, then this array
C             contains X   , the estimated filtered state at
C                       i+1
C             instant i+1.
C             On exit, if JOBX = 'N', or JOBX = 'X' and INFO = 1, then
C                                  -1
C             this array contains S   X   .
C                                  i+1 i+1
C
C     RINVY   (input) DOUBLE PRECISION array, dimension (P)
C                                      -1/2
C             This array must contain R    Y   , the product of the
C                                      i+1  i+1
C                                      -1/2
C             upper triangular matrix R     and the measured output
C                                      i+1
C             vector Y    at instant i+1.
C                     i+1
C
C     Z       (input) DOUBLE PRECISION array, dimension (M)
C             This array must contain Z , the mean value of the state
C                                      i
C             process noise at instant i.
C
C     E       (output) DOUBLE PRECISION array, dimension (P)
C             This array contains E   , the estimated error at instant
C                                  i+1
C             i+1.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If JOBX = 'X', then TOL is used to test for near
C                                        -1
C             singularity of the matrix S   . If the user sets
C                                        i+1
C             TOL > 0, then the given value of TOL is used as a
C             lower bound for the reciprocal condition number of that
C             matrix; a matrix whose estimated condition number is less
C             than 1/TOL is considered to be nonsingular. If the user
C             sets TOL <= 0, then an implicitly computed, default
C             tolerance, defined by TOLDEF = N*N*EPS, is used instead,
C             where EPS is the machine precision (see LAPACK Library
C             routine DLAMCH).
C             Otherwise, TOL is not referenced.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             where LIWORK = N if JOBX = 'X',
C             and   LIWORK = 1 otherwise.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.  If INFO = 0 and JOBX = 'X', DWORK(2) returns
C             an estimate of the reciprocal of the condition number
C                                 -1
C             (in the 1-norm) of S   .
C                                 i+1
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= MAX(1,N*(N+2*M)+3*M,(N+P)*(N+1)+N+MAX(N-1,M+1)),
C                                 if JOBX = 'N';
C             LDWORK >= MAX(2,N*(N+2*M)+3*M,(N+P)*(N+1)+N+MAX(N-1,M+1),
C                           3*N), if JOBX = 'X'.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;                        -1
C             = 1:  if JOBX = 'X' and the matrix S    is singular,
C                                                 i+1       -1
C                   i.e., the condition number estimate of S    (in the
C                                                           i+1
C                                                         -1    -1/2
C                   1-norm) exceeds 1/TOL.  The matrices S   , Q
C                                                         i+1   i
C                   and E have been computed.
C
C     METHOD
C
C     The routine performs one recursion of the square root information
C     filter algorithm, summarized as follows:
C
C       |    -1/2             -1/2    |     |         -1/2             |
C       |   Q         0      Q    Z   |     | (QINOV )     *     *     |
C       |    i                i    i  |     |       i                  |
C       |                             |     |                          |
C       |           -1/2      -1/2    |     |             -1    -1     |
C     T |    0     R    C    R    Y   |  =  |    0       S     S   X   |
C       |           i+1  i+1  i+1  i+1|     |             i+1   i+1 i+1|
C       |                             |     |                          |
C       |  -1 -1     -1 -1    -1      |     |                          |
C       | S  A  B   S  A     S  X     |     |    0         0     E     |
C       |  i         i        i  i    |     |                     i+1  |
C
C                   (Pre-array)                      (Post-array)
C
C     where T is an orthogonal transformation triangularizing the
C                        -1/2
C     pre-array, (QINOV )     is the inverse of the covariance square
C                      i
C     root (right Cholesky factor) of the process noise innovation
C                                                            -1  -1
C     (hence the information square root) at instant i and (A  ,A  B) is
C     in upper controller Hessenberg form.
C
C     An example of the pre-array is given below (where N = 6, M = 2,
C     and P = 3):
C
C         |x x |             | x|
C         |  x |             | x|
C         _______________________
C         |    | x x x x x x | x|
C         |    | x x x x x x | x|
C         |    | x x x x x x | x|
C         _______________________
C         |x x | x x x x x x | x|
C         |  x | x x x x x x | x|
C         |    | x x x x x x | x|
C         |    |   x x x x x | x|
C         |    |     x x x x | x|
C         |    |       x x x | x|
C
C     The inverse of the corresponding state covariance matrix P
C                                                               i+1|i+1
C     (hence the information matrix I) is then factorized as
C
C                    -1         -1     -1
!C         I       = P       = (S   )' S
C          i+1|i+1   i+1|i+1    i+1    i+1
C
C     and one combined time and measurement update for the state is
C     given by X   .
C               i+1
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
C     [2] Van Dooren, P. and Verhaegen, M.H.G.
C         Condensed Forms for Efficient Time-Invariant Kalman Filtering.
C         SIAM J. Sci. Stat. Comp., 9. pp. 516-530, 1988.
C
C     [3] Verhaegen, M.H.G. and Van Dooren, P.
C         Numerical Aspects of Different Kalman Filter Implementations.
C         IEEE Trans. Auto. Contr., AC-31, pp. 907-917, Oct. 1986.
C
C     [4] Vanbegin, M., Van Dooren, P., and Verhaegen, M.H.G.
C         Algorithm 675: FORTRAN Subroutines for Computing the Square
C         Root Covariance Filter and Square Root Information Filter in
C         Dense or Hessenberg Forms.
C         ACM Trans. Math. Software, 15, pp. 243-256, 1989.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires approximately
C
C           3    2                           2          3
C     (1/6)N  + N x (3/2 x M + P) + 2 x N x M  + 2/3 x M
C
C     operations and is backward stable (see [3]).
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C     Supersedes Release 2.0 routine FB01HD by M. Vanbegin,
C     P. Van Dooren, and M.H.G. Verhaegen.
C
C     REVISIONS
C
C     February 20, 1998, November 20, 2003, February 14, 2004.
C
C     KEYWORDS
C
C     Controller Hessenberg form, Kalman filtering, optimal filtering,
C     orthogonal transformation, recursive estimation, square-root
C     filtering, square-root information filtering.
C
C     ******************************************************************
C
#endif
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
     use omp_lib
#endif
       implicit none
!C     .. Parameters ..
      DOUBLE PRECISION  ONE, TWO
      PARAMETER         ( ONE = 1.0D0, TWO = 2.0D0 )
!C     .. Scalar Arguments ..
      CHARACTER         JOBX, MULTRC
      INTEGER           INFO, LDAINB, LDAINV, LDC, LDQINV, LDRINV, &
                        LDSINV, LDWORK, M, N, P
      DOUBLE PRECISION  TOL
!C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      INTEGER           IWORK(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      INTEGER           IWORK(*)
!DIR$ ASSUME_ALIGNED IWORK:64
#endif
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION  AINV(LDAINV,*), AINVB(LDAINB,*), C(LDC,*), &
                       DWORK(*), E(*), QINV(LDQINV,*), RINV(LDRINV,*), &
                       RINVY(*), SINV(LDSINV,*), X(*), Z(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
       DOUBLE PRECISION  AINV(LDAINV,*), AINVB(LDAINB,*), C(LDC,*), &
                       DWORK(*), E(*), QINV(LDQINV,*), RINV(LDRINV,*), &
                       RINVY(*), SINV(LDSINV,*), X(*), Z(*)
!DIR$ ASSUME_ALIGNED ANIV:64
!DIR$ ASSUME_ALIGNED AINVB:64
!DIR$ ASSUME_ALIGNED C:64
!DIR$ ASSUME_ALIGNED DWORK:64
!DIR$ ASSUME_ALIGNED E:64
!DIR$ ASSUME_ALIGNED QINV:64
!DIR$ ASSUME_ALIGNED RINV:64
!DIR$ ASSUME_ALIGNED RINVY:64
!DIR$ ASSUME_ALIGNED SINV:64
!DIR$ ASSUME_ALIGNED X:64
!DIR$ ASSUME_ALIGNED Z:64
#endif

!C     .. Local Scalars ..
      LOGICAL           LJOBX, LMULTR
      INTEGER           I, I12, I13, I23, I32, I33, II, IJ, ITAU, JWORK, &
                        LDW, M1, MP1, N1, NM, NP, WRKOPT
      DOUBLE PRECISION  RCOND
!C     .. External Functions ..
     
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
!C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DLACPY, DTRMM, DTRMV, MB02OD, &
                        MB04ID, MB04KD
!C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
!C     .. Executable Statements ..
!C
      NP = N + P
      NM = N + M
      N1 = MAX( 1, N )
      M1 = MAX( 1, M )
      MP1 = M + 1
      INFO = 0
      LJOBX  = LSAME( JOBX, 'X' )
      LMULTR = LSAME( MULTRC, 'P' )
!C
!C     Test the input scalar arguments.
!C
      IF( .NOT.LJOBX .AND. .NOT.LSAME( JOBX, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LMULTR .AND. .NOT.LSAME( MULTRC, 'N' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDSINV.LT.N1 ) THEN
         INFO = -7
      ELSE IF( LDAINV.LT.N1 ) THEN
         INFO = -9
      ELSE IF( LDAINB.LT.N1 ) THEN
         INFO = -11
      ELSE IF( LDRINV.LT.1 .OR. ( .NOT.LMULTR .AND. LDRINV.LT.P ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( LDQINV.LT.M1 ) THEN
         INFO = -17
      ELSE IF( ( LJOBX .AND. LDWORK.LT.MAX( 2, N*(NM + M) + 3*M, &
                                           NP*(N + 1) + N + &
                                           MAX( N - 1, MP1 ), 3*N ) ) &
                                                                  .OR. &
         ( .NOT.LJOBX .AND. LDWORK.LT.MAX( 1, N*(NM + M) + 3*M, &
                                           NP*(N + 1) + N + &
                                           MAX( N - 1, MP1 ) ) ) ) THEN &
         INFO = -25
      END IF

      IF ( INFO.NE.0 ) THEN
!C
!C        Error return.
!C
         RETURN
      END IF
!C
!C     Quick return if possible.
!C
      IF ( MAX( N, P ).EQ.0 ) THEN
         IF ( LJOBX ) THEN
            DWORK(1) = TWO
            DWORK(2) = ONE
         ELSE
            DWORK(1) = ONE
         END IF
         RETURN
      END IF
!C
!C     Construction of the needed part of the pre-array in DWORK.
!C     To save workspace, only the blocks (1,3), (3,1)-(3,3), (2,2), and
!C     (2,3) will be constructed when needed as shown below.
!C
!C     Storing SINV x AINVB and SINV x AINV in the (1,1) and (1,2)
!C     blocks of DWORK, respectively. The upper trapezoidal structure of
!C     [ AINVB AINV ] is fully exploited. Specifically, if M <= N, the
!C     following partition is used:
!C
!C       [ S1  S2 ] [ B1  A1 A3 ]
!C       [ 0   S3 ] [ 0   A2 A4 ],
!C
!C     where B1, A3, and S1 are M-by-M matrices, A1 and S2 are
!C     M-by-(N-M), A2 and S3 are (N-M)-by-(N-M), A4 is (N-M)-by-M, and
!C     B1, S1, A2, and S3 are upper triangular. The right hand side
!C     matrix above is stored in the workspace. If M > N, the partition
!C     is [ SINV ] [ B1 B2  A ], where B1 is N-by-N, B2 is N-by-(M-N),
!C     and B1 and SINV are upper triangular.
!C     The variables called Ixy define the starting positions where the
!C     (x,y) blocks of the pre-array are initially stored in DWORK.
!C     Workspace: need N*(M+N).
!C
!C     (Note: Comments in the code beginning "Workspace:" describe the
!C     minimal amount of real workspace needed at that point in the
!C     code, as well as the preferred amount for good performance.
!C     NB refers to the optimal block size for the immediately
!C     following subroutine, as returned by ILAENV.)
!C
      LDW = N1
      I32 = N*M + 1

      CALL DLACPY( 'Upper', N, M, AINVB, LDAINB, DWORK, LDW )
      CALL DLACPY( 'Full',  MIN( M, N ), N, AINV, LDAINV, DWORK(I32), &
                  LDW )
      IF ( N.GT.M ) &
        CALL DLACPY( 'Upper', N-M, N, AINV(MP1,1), LDAINV, &
                      DWORK(I32+M), LDW )
!C
!C                    [ B1  A1 ]
!C     Compute SINV x [ 0   A2 ] or SINV x B1 as a product of upper
!C     triangular matrices.
!C     Workspace: need N*(M+N+1).
!C
      II  = 1
      I13 = N*NM + 1
      WRKOPT = MAX( 1, N*NM + N )
!C
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(I,II)
#endif
      DO 10 I = 1, N
         CALL DCOPY( I, DWORK(II), 1, DWORK(I13), 1 )
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I, SINV,
     $               LDSINV, DWORK(I13), 1 )
         CALL DCOPY( I, DWORK(I13), 1, DWORK(II), 1 )
         II = II + N
   10 CONTINUE
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
!$OMP PARALLEL DO END
#endif
!C
!C                    [ A3 ]
!C     Compute SINV x [ A4 ] or SINV x [ B2 A ].
!C
      CALL DTRMM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, M, &
                  ONE, SINV, LDSINV, DWORK(II), LDW ) 
!C
!C     Storing the process noise mean value in (1,3) block of DWORK.
!C     Workspace: need N*(M+N) + M.
!C
      CALL DCOPY( M, Z, 1, DWORK(I13), 1 )
      CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', M, QINV, LDQINV, &
                  DWORK(I13), 1 )
!C
!C     Computing SINV x X in X.
!C
      CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', N, SINV, LDSINV, &
                 X, 1 )
!C
!C     Triangularization (2 steps).
!C
!C     Step 1: annihilate the matrix SINV x AINVB.
!C     Workspace: need N*(N+2*M) + 3*M.
!C
      I12   = I13  + M
      ITAU  = I12  + M*N
      JWORK = ITAU + M
!C
      CALL MB04KD( 'Upper', M, N, N, QINV, LDQINV, DWORK, LDW, &
                  DWORK(I32), LDW, DWORK(I12), M1, DWORK(ITAU), &
                  DWORK(JWORK) )
      WRKOPT = MAX( WRKOPT, N*( NM + M ) + 3*M )
!C
      IF ( N.EQ.0 ) THEN
         CALL DCOPY( P, RINVY, 1, E, 1 )
         IF ( LJOBX ) &
           DWORK(2) = ONE
         DWORK(1) = WRKOPT
         RETURN
      END IF
!C
!C     Apply the transformations to the last column of the pre-array.
!C     (Only the updated (3,3) block is now needed.)
!C
      IJ = 1
      !C
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(I,IJ)
#endif
      DO 20 I = 1, M
         CALL DAXPY( MIN( I, N ), -DWORK(ITAU+I-1)*( DWORK(I13+I-1) + &
                    DDOT( MIN( I, N ), DWORK(IJ), 1, X, 1 ) ), &
                                       DWORK(IJ), 1, X, 1 )
         IJ = IJ + N
20       CONTINUE
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
!$OMP PARALLEL DO END
#endif
!C
!C     Now, the workspace for SINV x AINVB, as well as for the updated
!C!     (1,2) block of the pre-array, are no longer needed.
!C     Move the computed (3,2) and (3,3) blocks of the pre-array in the
!C     (1,1) and (1,2) block positions of DWORK, to save space for the
!C     following computations.
!C     Then, adjust the implicitly defined leading dimension of DWORK,
!C     to make space for storing the (2,2) and (2,3) blocks of the
!C     pre-array.
!C     Workspace: need (P+N)*(N+1).
!C
      CALL DLACPY( 'Full', MIN( M, N ), N, DWORK(I32), LDW, DWORK, LDW )
      IF ( N.GT.M ) &
        CALL DLACPY( 'Upper', N-M, N, DWORK(I32+M), LDW, DWORK(MP1), &
                      LDW )
      LDW = MAX( 1, NP )
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(I)
#endif
      DO 40 I = N, 1, -1
         DO 30 IJ = MIN( N, I+M ), 1, -1
            DWORK(NP*(I-1)+P+IJ) = DWORK(N*(I-1)+IJ)
   30    CONTINUE
   40 CONTINUE
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
!$OMP PARALLEL DO END
#endif
!C     Copy of RINV x C in the (1,1) block of DWORK.
!C
      CALL DLACPY( 'Full', P, N, C, LDC, DWORK, LDW )
      IF ( .NOT.LMULTR ) &
        CALL DTRMM( 'Left', 'Upper', 'No transpose', 'Non-unit', P, N, &
                    ONE, RINV, LDRINV, DWORK, LDW )
!C
!C     Copy the inclusion measurement in the (1,2) block and the updated
!C     X in the (2,2) block of DWORK.
!C
      I23 = NP*N + 1
      I33 = I23  + P
      CALL DCOPY( P, RINVY, 1, DWORK(I23), 1 )
      CALL DCOPY( N, X, 1, DWORK(I33), 1 )
      WRKOPT = MAX( WRKOPT, NP*( N + 1 ) )
!C
!C     Step 2: QR factorization of the first block column of the matrix
!C
!C        [ RINV x C     RINV x Y ],
!C        [ SINV x AINV  SINV x X ]
!C
!C     where the second block row was modified at Step 1.
!C     Workspace: need   (P+N)*(N+1) + N + MAX(N-1,M+1);
!C                prefer (P+N)*(N+1) + N + (M+1)*NB, where NB is the
!C                       optimal block size for DGEQRF called in MB04ID.
!C
      ITAU  = I23  + NP
      JWORK = ITAU + N

      CALL MB04ID( NP, N, MAX( N-MP1, 0 ), 1, DWORK, LDW, DWORK(I23), &
                  LDW, DWORK(ITAU), DWORK(JWORK), LDWORK-JWORK+1, &
                  INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
!C
!C     Output SINV, X, and E and set the optimal workspace dimension
!C     (and the reciprocal of the condition number estimate).
!C
      CALL DLACPY( 'Upper', N, N, DWORK, LDW, SINV, LDSINV )
      CALL DCOPY( N, DWORK(I23), 1, X, 1 )
      IF( P.GT.0 ) &
         CALL DCOPY( P, DWORK(I23+N), 1, E, 1 )
!C
      IF ( LJOBX ) THEN
!C
!C        Compute X.
!C        Workspace: need 3*N.
!C
         CALL MB02OD( 'Left', 'Upper', 'No transpose', 'Non-unit', &
                     '1-norm', N, 1, ONE, SINV, LDSINV, X, N, RCOND, &
                     TOL, IWORK, DWORK, INFO )
         IF ( INFO.EQ.0 ) THEN
            WRKOPT = MAX( WRKOPT, 3*N )
            DWORK(2) = RCOND
         END IF
      END IF

      DWORK(1) = WRKOPT

END SUBROUTINE FB01TD

    
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE MB04KD( UPLO, N, M, P, R, LDR, A, LDA, B, LDB, C, LDC, &
     TAU, DWORK ) !GCC$ ATTRIBUTES hot :: MB04KD !GCC$ ATTRIBUTES aligned(32) :: MB04KD !GCC$ ATTRIBUTES no_stack_protector :: MB04KD
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE MB04KD( UPLO, N, M, P, R, LDR, A, LDA, B, LDB, C, LDC, &
     TAU, DWORK )
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: MB04KD
  !DIR$ OPTIMIZE : 2
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: MB04KD
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To calculate a QR factorization of the first block column and
C     apply the orthogonal transformations (from the left) also to the
C     second block column of a structured matrix, as follows
C                          _
C            [ R   0 ]   [ R   C ]
!C       Q' * [       ] = [       ]
C            [ A   B ]   [ 0   D ]
C                 _
C     where R and R are upper triangular. The matrix A can be full or
C     upper trapezoidal/triangular. The problem structure is exploited.
C     This computation is useful, for instance, in combined measurement
C     and time update of one iteration of the Kalman filter (square
C     root information filter).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Indicates if the matrix A is or not triangular as follows:
C             = 'U':  Matrix A is upper trapezoidal/triangular;
C             = 'F':  Matrix A is full.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER                 _
C             The order of the matrices R and R.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of the matrices B, C and D.  M >= 0.
C
C     P       (input) INTEGER
C             The number of rows of the matrices A, B and D.  P >= 0.
C
C     R       (input/output) DOUBLE PRECISION array, dimension (LDR,N)
C             On entry, the leading N-by-N upper triangular part of this
C             array must contain the upper triangular matrix R.
C             On exit, the leading N-by-N upper triangular part of this
C                                                        _
C             array contains the upper triangular matrix R.
C             The strict lower triangular part of this array is not
C             referenced.
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,N).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, if UPLO = 'F', the leading P-by-N part of this
C             array must contain the matrix A. If UPLO = 'U', the
C             leading MIN(P,N)-by-N part of this array must contain the
C             upper trapezoidal (upper triangular if P >= N) matrix A,
C             and the elements below the diagonal are not referenced.
C             On exit, the leading P-by-N part (upper trapezoidal or
C             triangular, if UPLO = 'U') of this array contains the
C             trailing components (the vectors v, see Method) of the
C             elementary reflectors used in the factorization.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,P).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading P-by-M part of this array must
C             contain the matrix B.
C             On exit, the leading P-by-M part of this array contains
C             the computed matrix D.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,P).
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,M)
C             The leading N-by-M part of this array contains the
C             computed matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,N).
C
C     TAU     (output) DOUBLE PRECISION array, dimension (N)
C             The scalar factors of the elementary reflectors used.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
C
C     METHOD
C
C     The routine uses N Householder transformations exploiting the zero
C     pattern of the block matrix.  A Householder matrix has the form
C
C                                     ( 1 ),
!C        H  = I - tau *u *u',    u  = ( v )
C         i          i  i  i      i   (  i)
C
C     where v  is a P-vector, if UPLO = 'F', or an min(i,P)-vector, if
C            i
C     UPLO = 'U'.  The components of v  are stored in the i-th column
C                                     i
C     of A, and tau  is stored in TAU(i).
C                  i
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTORS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary reflector, QR factorization, orthogonal transformation.
C
C     ******************************************************************
C
#endif
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
use omp_lib
#endif
       implicit none
!C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
!C     .. Scalar Arguments ..
      CHARACTER         UPLO
      INTEGER           LDA, LDB, LDC, LDR, M, N, P
      !C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           R(LDR,*), TAU(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER))
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           R(LDR,*), TAU(*)
      !DIR$ ASSUME_ALIGNED A:64
      !DIR$ ASSUME_ALIGNED B:64
      !DIR$ ASSUME_ALIGNED C:64
      !DIR$ ASSUME_ALIGNED DWORK:64
      !DIR$ ASSUME_ALIGNED R:64
      !DIR$ ASSUME_ALIGNED TAU:64
#endif
!C     .. Local Scalars ..
      LOGICAL           LUPLO
      INTEGER           I, IM

!C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DGER, DLARFG, DSCAL
!C     .. Intrinsic Functions ..
      INTRINSIC         MIN
!C     .. Executable Statements ..
!C
      IF( MIN( N, P ).EQ.0 ) RETURN
       
!C
      LUPLO = LSAME( UPLO, 'U' )
      IM = P
      !C
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(I,IM)
#endif
      DO 10 I = 1, N
!C
!C        Annihilate the I-th column of A and apply the transformations
!C        to the entire block matrix, exploiting its structure.
!C
         IF( LUPLO ) IM = MIN( I, P )
         CALL DLARFG( IM+1, R(I,I), A(1,I), 1, TAU(I) )
         IF( TAU(I).NE.ZERO ) THEN
!C
!C                                      [ R(I,I+1:N)        0     ]
!C           [ w C(I,:) ] := [ 1 v' ] * [                         ]
!C                                      [ A(1:IM,I+1:N) B(1:IM,:) ]
!C
            IF( I.LT.N ) THEN
               CALL DCOPY( N-I, R(I,I+1), LDR, DWORK, 1 )
               CALL DGEMV( 'Transpose', IM, N-I, ONE, A(1,I+1), LDA, &
                           A(1,I), 1, ONE, DWORK, 1 )
            END IF
            CALL DGEMV( 'Transpose', IM, M, ONE, B, LDB, A(1,I), 1, &
                        ZERO, C(I,1), LDC )
!C
!C           [ R(I,I+1:N)      C(I,:)  ]    [ R(I,I+1:N)        0     ]
!C           [                         ] := [                         ]
!C           [ A(1:IM,I+1:N) D(1:IM,:) ]    [ A(1:IM,I+1:N) B(1:IM,:) ]
!C
!C                                                 [ 1 ]
!C                                         - tau * [   ] * [ w C(I,:) ]
!C                                                 [ v ]
!C
            IF( I.LT.N ) THEN
               CALL DAXPY( N-I, -TAU(I), DWORK, 1, R(I,I+1), LDR )
               CALL DGER( IM, N-I, -TAU(I), A(1,I), 1, DWORK, 1, &
                          A(1,I+1), LDA )
            END IF
            CALL DSCAL( M, -TAU(I), C(I,1), LDC )
            CALL DGER( IM, M, ONE, A(1,I), 1, C(I,1), LDC, B, LDB )
         END IF
10       CONTINUE
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
         !$OMP END PARALLEL DO
#endif
     
END SUBROUTINE MB04KD

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE MB04ID(N, M, P, L, A, LDA, B, LDB, TAU, DWORK, LDWORK,INFO) !GCC$ ATTRIBUTES hot :: MB04ID !GCC$ ATTRIBUTES aligned(32) :: MB04ID !GCC$ ATTRIBUTES no_stack_protector :: MB04ID
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE MB04ID(N, M, P, L, A, LDA, B, LDB, TAU, DWORK, LDWORK,INFO)
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: MB04ID
    !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: MB04ID
#endif
#if 0                     
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To compute a QR factorization of an n-by-m matrix A (A = Q * R),
C     having a p-by-min(p,m) zero triangle in the lower left-hand side
C     corner, as shown below, for n = 8, m = 7, and p = 2:
C
C            [ x x x x x x x ]
C            [ x x x x x x x ]
C            [ x x x x x x x ]
C            [ x x x x x x x ]
C        A = [ x x x x x x x ],
C            [ x x x x x x x ]
C            [ 0 x x x x x x ]
C            [ 0 0 x x x x x ]
C
C     and optionally apply the transformations to an n-by-l matrix B
C     (from the left). The problem structure is exploited. This
C     computation is useful, for instance, in combined measurement and
C     time update of one iteration of the time-invariant Kalman filter
C     (square root information filter).
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of rows of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of the matrix A.  M >= 0.
C
C     P       (input) INTEGER
C             The order of the zero triagle.  P >= 0.
C
C     L       (input) INTEGER
C             The number of columns of the matrix B.  L >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
C             On entry, the leading N-by-M part of this array must
C             contain the matrix A. The elements corresponding to the
C             zero P-by-MIN(P,M) lower trapezoidal/triangular part
C             (if P > 0) are not referenced.
C             On exit, the elements on and above the diagonal of this
C             array contain the MIN(N,M)-by-M upper trapezoidal matrix
C             R (R is upper triangular, if N >= M) of the QR
C             factorization, and the relevant elements below the
C             diagonal contain the trailing components (the vectors v,
C             see Method) of the elementary reflectors used in the
C             factorization.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,L)
C             On entry, the leading N-by-L part of this array must
C             contain the matrix B.
C             On exit, the leading N-by-L part of this array contains
C             the updated matrix B.
C             If L = 0, this array is not referenced.
C
C     LDB     INTEGER
C             The leading dimension of array B.
C             LDB >= MAX(1,N) if L > 0;
C             LDB >= 1        if L = 0.
C
C     TAU     (output) DOUBLE PRECISION array, dimension MIN(N,M)
C             The scalar factors of the elementary reflectors used.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= MAX(1,M-1,M-P,L).
C             For optimum performance LDWORK should be larger.
C
C             If LDWORK = -1, then a workspace query is assumed;
C             the routine only calculates the optimal size of the
C             DWORK array, returns this value as the first entry of
C             the DWORK array, and no error message related to LDWORK
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
C     The routine uses min(N,M) Householder transformations exploiting
C     the zero pattern of the matrix.  A Householder matrix has the form
C
!C                                     ( 1 ),
!C        H  = I - tau *u *u',    u  = ( v )
C         i          i  i  i      i   (  i)
C
C     where v  is an (N-P+I-2)-vector.  The components of v  are stored
C            i                                             i
C     in the i-th column of A, beginning from the location i+1, and
C     tau  is stored in TAU(i).
C        i
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTORS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Jan. 2009,
C     Apr. 2009, Apr. 2011.
C
C     KEYWORDS
C
C     Elementary reflector, QR factorization, orthogonal transformation.
C
C     ******************************************************************
C
#endif
#if(GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
      use omp_lib
#endif
      implicit none
!C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
!C     .. Scalar Arguments ..
      INTEGER           INFO, L, LDA, LDB, LDWORK, M, N, P
!C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), TAU(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
       DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), TAU(*)
!DIR$ ASSUME_ALIGNED A:64
!DIR$ ASSUME_ALIGNED B:64
!DIR$ ASSUME_ALIGNED DWORK:64
!DIR$ ASSUME_ALIGNED TAU:64
#endif
!C     .. Local Scalars ..
      LOGICAL           LQUERY
      INTEGER           I, WRKOPT
      DOUBLE PRECISION  FIRST
!C     .. External Subroutines ..
      EXTERNAL          DGEQRF, DLARF, DLARFG, DORMQR
!C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
!C     .. Executable Statements ..
!C
!C     Test the input scalar arguments.
!C
      INFO = 0
      LQUERY = LDWORK.EQ.-1
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( L.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.1 .OR. ( L.GT.0 .AND. LDB.LT.N ) ) THEN
         INFO = -8
      ELSE
         I = MAX( 1, M - 1, M - P, L )
         IF( LQUERY ) THEN
            IF( M.GT.P ) THEN
               CALL DGEQRF( N-P, M-P, A, LDA, TAU, DWORK, -1, INFO )
               WRKOPT = MAX( I, INT( DWORK( 1 ) ) )
               IF ( L.GT.0 ) THEN
                  CALL DORMQR( 'Left', 'Transpose', N-P, L, MIN(N,M)-P, &
                              A, LDA, TAU, B, LDB, DWORK, -1, INFO )
                  WRKOPT = MAX( WRKOPT, INT( DWORK( 1 ) ) )
               END IF
            END IF
         ELSE IF( LDWORK.LT.I ) THEN
            INFO = -11
         END IF
      END IF
!C
      IF( INFO.NE.0 ) THEN
!C
!C        Error return.
!C
          RETURN
      ELSE IF( LQUERY ) THEN
         DWORK(1) = WRKOPT
         RETURN
      END IF
!C
!C     Quick return if possible.
!C
      IF( MIN( M, N ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      ELSE IF( N.LE.P+1 ) THEN
         DO 5 I = 1, MIN( N, M )
            TAU(I) = ZERO
    5    CONTINUE
         DWORK(1) = ONE
         RETURN
      END IF
!C
!C     Annihilate the subdiagonal elements of A and apply the
!C     transformations to B, if L > 0.
!C     Workspace: need MAX(M-1,L).
!C
!C     (Note: Comments in the code beginning "Workspace:" describe the
!C     minimal amount of real workspace needed at that point in the
!C     code, as well as the preferred amount for good performance.
!C     NB refers to the optimal block size for the immediately
!C     following subroutine, as returned by ILAENV.)
!C
#if(GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(I,FIRST)
#endif
      DO 10 I = 1, MIN( P, M )
!C
!C        Exploit the structure of the I-th column of A.
!C
         CALL DLARFG( N-P, A(I,I), A(I+1,I), 1, TAU(I) )
         IF( TAU(I).NE.ZERO ) THEN
!C
            FIRST = A(I,I)
            A(I,I) = ONE
!C
            IF ( I.LT.M ) CALL DLARF( 'Left', N-P, M-I, A(I,I), 1,  &
                                     TAU(I), A(I,I+1), LDA, DWORK )
            IF ( L.GT.0 ) CALL DLARF( 'Left', N-P, L, A(I,I), 1, TAU(I), &
                                     B(I,1), LDB, DWORK )
!C
            A(I,I) = FIRST
         END IF
10       CONTINUE
#if(GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
!$OMP END PARALLEL DO
#endif

      WRKOPT = MAX( 1, M - 1, L )
!C
!C     Fast QR factorization of the remaining right submatrix, if any.
!C     Workspace: need M-P;  prefer (M-P)*NB.
!C
      IF( M.GT.P ) THEN
         CALL DGEQRF( N-P, M-P, A(P+1,P+1), LDA, TAU(P+1), DWORK,  &
                     LDWORK, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
!C
         IF ( L.GT.0 ) THEN
!C
!C           Apply the transformations to B.
!C           Workspace: need L;  prefer L*NB.
!C
            CALL DORMQR( 'Left', 'Transpose', N-P, L, MIN(N,M)-P,  &
                        A(P+1,P+1), LDA, TAU(P+1), B(P+1,1), LDB,  &
                        DWORK, LDWORK, INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
         END IF
      END IF
!C
      DWORK(1) = WRKOPT
    
END SUBROUTINE MB04ID

   

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE MB04LD( UPLO, N, M, P, L, LDL, A, LDA, B, LDB, C, LDC &
                  TAU, DWORK ) !GCC$ ATTRIBUTES hot :: MB04LD !GCC$ ATTRIBUTES aligned(32) :: MB04LD !GCC$ ATTRIBUTES no_stack_protector :: MB04LD
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE MB04LD( UPLO, N, M, P, L, LDL, A, LDA, B, LDB, C, LDC &
                  TAU, DWORK )
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: MB04LD
!DIR$ OPTIMIZE : 2
!DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: MB04LD
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To calculate an LQ factorization of the first block row and apply
C     the orthogonal transformations (from the right) also to the second
C     block row of a structured matrix, as follows
C                        _
C        [ L   A ]     [ L   0 ]
C        [       ]*Q = [       ]
C        [ 0   B ]     [ C   D ]
C                 _
C     where L and L are lower triangular. The matrix A can be full or
C     lower trapezoidal/triangular. The problem structure is exploited.
C     This computation is useful, for instance, in combined measurement
C     and time update of one iteration of the Kalman filter (square
C     root covariance filter).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Indicates if the matrix A is or not triangular as follows:
C             = 'L':  Matrix A is lower trapezoidal/triangular;
C             = 'F':  Matrix A is full.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER                 _
C             The order of the matrices L and L.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of the matrices A, B and D.  M >= 0.
C
C     P       (input) INTEGER
C             The number of rows of the matrices B, C and D.  P >= 0.
C
C     L       (input/output) DOUBLE PRECISION array, dimension (LDL,N)
C             On entry, the leading N-by-N lower triangular part of this
C             array must contain the lower triangular matrix L.
C             On exit, the leading N-by-N lower triangular part of this
C                                                        _
C             array contains the lower triangular matrix L.
C             The strict upper triangular part of this array is not
C             referenced.
C
C     LDL     INTEGER
C             The leading dimension of array L.  LDL >= MAX(1,N).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
C             On entry, if UPLO = 'F', the leading N-by-M part of this
C             array must contain the matrix A. If UPLO = 'L', the
C             leading N-by-MIN(N,M) part of this array must contain the
C             lower trapezoidal (lower triangular if N <= M) matrix A,
C             and the elements above the diagonal are not referenced.
C             On exit, the leading N-by-M part (lower trapezoidal or
C             triangular, if UPLO = 'L') of this array contains the
C             trailing components (the vectors v, see Method) of the
C             elementary reflectors used in the factorization.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading P-by-M part of this array must
C             contain the matrix B.
C             On exit, the leading P-by-M part of this array contains
C             the computed matrix D.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,P).
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array contains the
C             computed matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     TAU     (output) DOUBLE PRECISION array, dimension (N)
C             The scalar factors of the elementary reflectors used.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
C
C     METHOD
C
C     The routine uses N Householder transformations exploiting the zero
C     pattern of the block matrix.  A Householder matrix has the form
C
!C                                     ( 1 ),
!C        H  = I - tau *u *u',    u  = ( v )
C         i          i  i  i      i   (  i)
C
C     where v  is an M-vector, if UPLO = 'F', or an min(i,M)-vector, if
C            i
C     UPLO = 'L'.  The components of v  are stored in the i-th row of A,
C                                     i
C     and tau  is stored in TAU(i).
C            i
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTORS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary reflector, LQ factorization, orthogonal transformation.
C
C     ******************************************************************
C
#endif
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
use omp_lib
#endif
       implicit none
!C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
!C     .. Scalar Arguments ..
      CHARACTER         UPLO
      INTEGER           LDA, LDB, LDC, LDL, M, N, P
      !C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           L(LDL,*), TAU(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), &
           L(LDL,*), TAU(*)
      !DIR$ ASSUME_ALIGNED A:64
      !DIR$ ASSUME_ALIGNED B:64
      !DIR$ ASSUME_ALIGNED C:64
      !DIR$ ASSUME_ALIGNED DWORK:64
      !DIR$ ASSUME_ALIGNED L:64
      !DIR$ ASSUME_ALIGNED TAU:64
#endif
!C     .. Local Scalars ..
      LOGICAL           LUPLO
      INTEGER           I, IM

!C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DGER, DLARFG, DSCAL
!C     .. Intrinsic Functions ..
      INTRINSIC         MIN
!C     .. Executable Statements ..
!C
      IF( MIN( M, N ).EQ.0 ) RETURN
       
!C
      LUPLO = LSAME( UPLO, 'L' )
      IM = M
      !C
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
      !$OMP PARALLEL DO DEAFULT(SHARED) PRIVATE(I,IM) SCHEDULE(STATIC)
#endif
      DO 10 I = 1, N
!C
!C        Annihilate the I-th row of A and apply the transformations to
!C        the entire block matrix, exploiting its structure.
!C
         IF( LUPLO ) IM = MIN( I, M )
         CALL DLARFG( IM+1, L(I,I), A(I,1), LDA, TAU(I) )
         IF( TAU(I).NE.ZERO ) THEN
!C
!C           [    w   ]    [ L(I+1:N,I) A(I+1:N,1:IM) ]   [ 1 ]
!C           [        ] := [                          ] * [   ]
!C           [ C(:,I) ]    [      0        B(:,1:IM)  ]   [ v ]
!C
            IF( I.LT.N ) THEN
               CALL DCOPY( N-I, L(I+1,I), 1, DWORK, 1 )
               CALL DGEMV( 'No transpose', N-I, IM, ONE, A(I+1,1), LDA, &
                           A(I,1), LDA, ONE, DWORK, 1 )
            END IF
            CALL DGEMV( 'No transpose', P, IM, ONE, B, LDB, A(I,1),  &
                        LDA, ZERO, C(1,I), 1 )
!C
!C           [ L(I+1:N,I) A(I+1:N,1:IM) ]    [ L(I+1:N,I) A(I+1:N,1:IM) ]
!C           [                          ] := [                          ]
!C           [   C(:,I)      D(:,1:IM)  ]    [      0        B(:,1:IM)  ]
!C
!C                                                 [    w   ]
!C                                         - tau * [        ] * [ 1 , v']
!C                                                 [ C(:,I) ]
!C
            IF( I.LT.N ) THEN
               CALL DAXPY( N-I, -TAU(I), DWORK, 1, L(I+1,I), 1 )
               CALL DGER( N-I, IM, -TAU(I), DWORK, 1, A(I,1), LDA,  &
                          A(I+1,1), LDA )
            END IF
            CALL DSCAL( P, -TAU(I), C(1,I), 1 )
            CALL DGER( P, IM, ONE, C(1,I), 1, A(I,1), LDA, B, LDB )
         END IF
   10 CONTINUE
#if (GMS_SLICOT_USE_REFERENCE_LAPACK) == 1
         !$OMP END PARALLEL DO
#endif
    
END SUBROUTINE MB04LD

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
SUBROUTINE SB01BD( DICO, N, M, NP, ALPHA, A, LDA, B, LDB, WR, WI, &
   NFP, NAP, NUP, F, LDF, Z, LDZ, TOL, DWORK, &
   LDWORK, IWARN, INFO ) !GCC$ ATTRIBUTES hot :: SB01BD !GCC$ ATTRIBUTES aligned(32) :: SB01BD !GCC$ ATTRIBUTES no_stack_protector :: SB01BD
#elif defined(__ICC) || defined(__INTEL_COMPILER)
SUBROUTINE SB01BD( DICO, N, M, NP, ALPHA, A, LDA, B, LDB, WR, WI, &
   NFP, NAP, NUP, F, LDF, Z, LDZ, TOL, DWORK, &
   LDWORK, IWARN, INFO )
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: SB01BD
!DIR$ OPTIMIZE : 3
!DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: SB01BD
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To determine the state feedback matrix F for a given system (A,B)
C     such that the closed-loop state matrix A+B*F has specified
C     eigenvalues.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the original system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The dimension of the state vector, i.e. the order of the
C             matrix A, and also the number of rows of the matrix B and
C             the number of columns of the matrix F.  N >= 0.
C
C     M       (input) INTEGER
C             The dimension of input vector, i.e. the number of columns
C             of the matrix B and the number of rows of the matrix F.
C             M >= 0.
C
C     NP      (input) INTEGER
C             The number of given eigenvalues. At most N eigenvalues
C             can be assigned.  0 <= NP.
C
C     ALPHA   (input) DOUBLE PRECISION
C             Specifies the maximum admissible value, either for real
C             parts, if DICO = 'C', or for moduli, if DICO = 'D',
C             of the eigenvalues of A which will not be modified by
C             the eigenvalue assignment algorithm.
C             ALPHA >= 0 if DICO = 'D'.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, the leading N-by-N part of this array contains
!C             the matrix Z'*(A+B*F)*Z in a real Schur form.
C             The leading NFP-by-NFP diagonal block of A corresponds
C             to the fixed (unmodified) eigenvalues having real parts
C             less than ALPHA, if DICO = 'C', or moduli less than ALPHA,
C             if DICO = 'D'. The trailing NUP-by-NUP diagonal block of A
C             corresponds to the uncontrollable eigenvalues detected by
C             the eigenvalue assignment algorithm. The elements under
C             the first subdiagonal are set to zero.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input/state matrix.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     WR,WI   (input/output) DOUBLE PRECISION array, dimension (NP)
C             On entry, these arrays must contain the real and imaginary
C             parts, respectively, of the desired eigenvalues of the
C             closed-loop system state-matrix A+B*F. The eigenvalues
C             can be unordered, except that complex conjugate pairs
C             must appear consecutively in these arrays.
C             On exit, if INFO = 0, the leading NAP elements of these
C             arrays contain the real and imaginary parts, respectively,
C             of the assigned eigenvalues. The trailing NP-NAP elements
C             contain the unassigned eigenvalues.
C
C     NFP     (output) INTEGER
C             The number of eigenvalues of A having real parts less than
C             ALPHA, if DICO = 'C', or moduli less than ALPHA, if
C             DICO = 'D'. These eigenvalues are not modified by the
C             eigenvalue assignment algorithm.
C
C     NAP     (output) INTEGER
C             The number of assigned eigenvalues. If INFO = 0 on exit,
C             then NAP = N-NFP-NUP.
C
C     NUP     (output) INTEGER
C             The number of uncontrollable eigenvalues detected by the
C             eigenvalue assignment algorithm (see METHOD).
C
C     F       (output) DOUBLE PRECISION array, dimension (LDF,N)
C             The leading M-by-N part of this array contains the state
C             feedback F, which assigns NAP closed-loop eigenvalues and
C             keeps unaltered N-NAP open-loop eigenvalues.
C
C     LDF     INTEGER
C             The leading dimension of array F.  LDF >= MAX(1,M).
C
C     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N)
C             The leading N-by-N part of this array contains the
C             orthogonal matrix Z which reduces the closed-loop
C             system state matrix A + B*F to upper real Schur form.
C
C     LDZ     INTEGER
C             The leading dimension of array Z.  LDZ >= MAX(1,N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The absolute tolerance level below which the elements of A
C             or B are considered zero (used for controllability tests).
C             If the user sets TOL <= 0, then the default tolerance
C             TOL = N * EPS * max(NORM(A),NORM(B)) is used, where EPS is
C             the machine precision (see LAPACK Library routine DLAMCH)
C             and NORM(A) denotes the 1-norm of A.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of working array DWORK.
C             LDWORK >= MAX( 1,5*M,5*N,2*N+4*M ).
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = K:  K violations of the numerical stability condition
C                   NORM(F) <= 100*NORM(A)/NORM(B) occured during the
C                   assignment of eigenvalues.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the reduction of A to a real Schur form failed;
C             = 2:  a failure was detected during the ordering of the
C                   real Schur form of A, or in the iterative process
!C                   for reordering the eigenvalues of Z'*(A + B*F)*Z
C                   along the diagonal.
C             = 3:  the number of eigenvalues to be assigned is less
C                   than the number of possibly assignable eigenvalues;
C                   NAP eigenvalues have been properly assigned,
C                   but some assignable eigenvalues remain unmodified.
C             = 4:  an attempt is made to place a complex conjugate
C                   pair on the location of a real eigenvalue. This
C                   situation can only appear when N-NFP is odd,
C                   NP > N-NFP-NUP is even, and for the last real
C                   eigenvalue to be modified there exists no available
C                   real eigenvalue to be assigned. However, NAP
C                   eigenvalues have been already properly assigned.
C
C     METHOD
C
C     SB01BD is based on the factorization algorithm of [1].
C     Given the matrices A and B of dimensions N-by-N and N-by-M,
C     respectively, this subroutine constructs an M-by-N matrix F such
C     that A + BF has eigenvalues as follows.
C     Let NFP eigenvalues of A have real parts less than ALPHA, if
C     DICO = 'C', or moduli less then ALPHA, if DICO = 'D'. Then:
C     1) If the pair (A,B) is controllable, then A + B*F has
C        NAP = MIN(NP,N-NFP) eigenvalues assigned from those specified
C        by WR + j*WI and N-NAP unmodified eigenvalues;
C     2) If the pair (A,B) is uncontrollable, then the number of
C        assigned eigenvalues NAP satifies generally the condition
C        NAP <= MIN(NP,N-NFP).
C
C     At the beginning of the algorithm, F = 0 and the matrix A is
C     reduced to an ordered real Schur form by separating its spectrum
C     in two parts. The leading NFP-by-NFP part of the Schur form of
C     A corresponds to the eigenvalues which will not be modified.
C     These eigenvalues have real parts less than ALPHA, if
C     DICO = 'C', or moduli less than ALPHA, if DICO = 'D'.
C     The performed orthogonal transformations are accumulated in Z.
C     After this preliminary reduction, the algorithm proceeds
C     recursively.
C
C     Let F be the feedback matrix at the beginning of a typical step i.
C     At each step of the algorithm one real eigenvalue or two complex
C     conjugate eigenvalues are placed by a feedback Fi of rank 1 or
C     rank 2, respectively. Since the feedback Fi affects only the
C     last 1 or 2 columns of Z'*(A+B*F)*Z, the matrix Z'*(A+B*F+B*Fi)*Z
C     therefore remains in real Schur form. The assigned eigenvalue(s)
C     is (are) then moved to another diagonal position of the real
C     Schur form using reordering techniques and a new block is
C     transfered in the last diagonal position. The feedback matrix F
C     is updated as F <-- F + Fi. The eigenvalue(s) to be assigned at
C     each step is (are) chosen such that the norm of each Fi is
C     minimized.
C
C     If uncontrollable eigenvalues are encountered in the last diagonal
!C     position of the real Schur matrix Z'*(A+B*F)*Z, the algorithm
C     deflates them at the bottom of the real Schur form and redefines
C     accordingly the position of the "last" block.
C
C     Note: Not all uncontrollable eigenvalues of the pair (A,B) are
C     necessarily detected by the eigenvalue assignment algorithm.
C     Undetected uncontrollable eigenvalues may exist if NFP > 0 and/or
C     NP < N-NFP.
C
C     REFERENCES
C
C     [1] Varga A.
C         A Schur method for pole assignment.
C         IEEE Trans. Autom. Control, Vol. AC-26, pp. 517-519, 1981.
C
C     NUMERICAL ASPECTS
C                                            3
C     The algorithm requires no more than 14N  floating point
C     operations. Although no proof of numerical stability is known,
C     the algorithm has always been observed to yield reliable
C     numerical results.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     February 1999. Based on the RASP routine SB01BD.
C
C     REVISIONS
C
C     March 30, 1999, V. Sima, Research Institute for Informatics,
C     Bucharest.
C     April 4, 1999. A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen.
C     May 18, 2003. A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen.
C     Feb. 15, 2004, V. Sima, Research Institute for Informatics,
C     Bucharest.
C     May 12, 2005. A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen.
C     Dec. 29, 2012, V. Sima, Research Institute for Informatics,
C     Bucharest.
C
C     KEYWORDS
C
C     Eigenvalues, eigenvalue assignment, feedback control,
C     pole placement, state-space model.
C
C     ******************************************************************
C
#endif
       implicit none
!C     .. Parameters ..
      DOUBLE PRECISION HUNDR, ONE, TWO, ZERO
      PARAMETER        ( HUNDR = 1.0D2, ONE = 1.0D0, TWO = 2.0D0, &
                        ZERO = 0.0D0 )
!C     .. Scalar Arguments ..
      CHARACTER        DICO
      INTEGER          INFO, IWARN, LDA, LDB, LDF, LDWORK, LDZ, M, N, &
                       NAP, NFP, NP, NUP
      DOUBLE PRECISION ALPHA, TOL
      !C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION A(LDA,*), B(LDB,*), DWORK(*), F(LDF,*), &
           WI(*), WR(*), Z(LDZ,*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION A(LDA,*), B(LDB,*), DWORK(*), F(LDF,*), &
           WI(*), WR(*), Z(LDZ,*)
      !DIR$ ASSUME_ALIGNED A:64
      !DIR$ ASSUME_ALIGNED B:64
      !DIR$ ASSUME_ALIGNED DWORK:64
      !DIR$ ASSUME_ALIGNED F:64
      !DIR$ ASSUME_ALIGNED WI:64
      !DIR$ ASSUME_ALIGNED WR:64
      !DIR$ ASSUME_ALIGNED Z:64
#endif
!C     .. Local Scalars ..
      LOGICAL          CEIG, DISCR, SIMPLB
      INTEGER          I, IB, IB1, IERR, IPC, J, K, KFI, KG, KW, KWI, &
                      KWR, NCUR, NCUR1, NL, NLOW, NMOVES, NPC, NPR,   &
                      NSUP, WRKOPT
      DOUBLE PRECISION ANORM, BNORM, C, P, RMAX, S, X, Y, TOLER, TOLERB
!C     .. Local Arrays ..
      LOGICAL          BWORK(1)
      DOUBLE PRECISION A2(2,2)
!C     .. External Functions ..
      LOGICAL           SELECT
      DOUBLE PRECISION  DLANGE
      EXTERNAL          DLANGE, SELECT
!C     .. External Subroutines ..
      EXTERNAL         DCOPY, DGEES, DGEMM, DLAEXC, DLASET, DROT, &
                       DSWAP
!C     .. Intrinsic Functions ..
      INTRINSIC        DBLE, INT, MAX
!C     ..
!C     .. Executable Statements ..
!C
      DISCR = LSAME( DICO, 'D' )
      IWARN = 0
      INFO  = 0
!C
!C     Check the scalar input parameters.
!C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( NP.LT.0 ) THEN
         INFO = -4
      ELSE IF( DISCR .AND. ( ALPHA.LT.ZERO ) ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
         INFO = -16
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -18
      ELSE IF( LDWORK.LT.MAX( 1, 5*M, 5*N, 2*N + 4*M ) ) THEN
         INFO = -21
      END IF
      IF( INFO.NE.0 )THEN
           RETURN
      END IF
!C
!C     Quick return if possible.
!C
      IF( N.EQ.0 ) THEN
         NFP = 0
         NAP = 0
         NUP = 0
         DWORK(1) = ONE
         RETURN
      END IF
!C
!C     Compute the norms of A and B, and set default tolerances
!C     if necessary.
!C
      ANORM = DLANGE( '1-norm', N, N, A, LDA, DWORK )
      BNORM = DLANGE( '1-norm', N, M, B, LDB, DWORK )
      IF( TOL.LE.ZERO ) THEN
         X = DLAMCH( 'Epsilon' )
         TOLER  = DBLE( N ) * MAX( ANORM, BNORM ) * X
         TOLERB = DBLE( N ) * BNORM * X
      ELSE
         TOLER  = TOL
         TOLERB = TOL
      END IF
!C
!C     Allocate working storage.
!C
      KWR = 1
      KWI = KWR + N
      KW  = KWI + N
!C
!C     Reduce A to real Schur form using an orthogonal similarity
!C     transformation A <- Z'*A*Z and accumulate the transformation in Z.
!C
!C     Workspace:  need   5*N;
!C                 prefer larger.
!C
      CALL DGEES( 'Vectors', 'No ordering', SELECT, N, A, LDA, NCUR, &
                 DWORK(KWR), DWORK(KWI), Z, LDZ, DWORK(KW), &
                 LDWORK-KW+1, BWORK, INFO )
      WRKOPT = KW - 1 + INT( DWORK( KW ) )
      IF( INFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
!C
!C     Reduce A to an ordered real Schur form using an orthogonal
!C     similarity transformation A <- Z'*A*Z and accumulate the
!C     transformations in Z. The separation of the spectrum of A is
!C     performed such that the leading NFP-by-NFP submatrix of A
!C     corresponds to the "good" eigenvalues which will not be
!C     modified. The bottom (N-NFP)-by-(N-NFP) diagonal block of A
!C     corresponds to the "bad" eigenvalues to be modified.
!C
!C     Workspace needed:  N.
!C
      CALL MB03QD( DICO, 'Stable', 'Update', N, 1, N, ALPHA, &
                  A, LDA, Z, LDZ, NFP, DWORK, INFO )
      IF( INFO.NE.0 ) &
          RETURN
!C
!C     Set F = 0.
!C
      CALL DLASET( 'Full', M, N, ZERO, ZERO, F, LDF )
!C
!C     Return if B is negligible (uncontrollable system).
!C
      IF( BNORM.LE.TOLERB ) THEN
         NAP = 0
         NUP = N
         DWORK(1) = WRKOPT
         RETURN
      END IF
!C
!C     Compute the bound for the numerical stability condition.
!C
      RMAX = HUNDR * ANORM / BNORM
!C
!C     Perform eigenvalue assignment if there exist "bad" eigenvalues.
!C
      NAP = 0
      NUP = 0
      IF( NFP.LT.N ) THEN
         KG  = 1
         KFI = KG  + 2*M
         KW  = KFI + 2*M
!C
!C        Set the limits for the bottom diagonal block.
!C
         NLOW = NFP + 1
         NSUP = N
!C
!C        Separate and count real and complex eigenvalues to be assigned.
!C
         NPR = 0
         DO 10 I = 1, NP
            IF( WI(I).EQ.ZERO ) THEN
               NPR = NPR + 1
               K = I - NPR
               IF( K.GT.0 ) THEN
                  S = WR(I)
                  DO 5 J = NPR + K - 1, NPR, -1
                     WR(J+1) = WR(J)
                     WI(J+1) = WI(J)
    5             CONTINUE
                  WR(NPR) = S
                  WI(NPR) = ZERO
               END IF
            END IF
   10    CONTINUE
         NPC = NP - NPR
!C
!C        The first NPR elements of WR and WI contain the real
!C        eigenvalues, the last NPC elements contain the complex
!C        eigenvalues. Set the pointer to complex eigenvalues.
!C
         IPC = NPR + 1
!C
!C        Main loop for assigning one or two eigenvalues.
!C
!C        Terminate if all eigenvalues were assigned, or if there
!C        are no more eigenvalues to be assigned, or if a non-fatal
!C        error condition was set.
!C
!C        WHILE (NLOW <= NSUP and INFO = 0) DO
!C
   20    IF( NLOW.LE.NSUP .AND. INFO.EQ.0 ) THEN
!C
!C           Determine the dimension of the last block.
!C
            IB = 1
            IF( NLOW.LT.NSUP ) THEN
               IF( A(NSUP,NSUP-1).NE.ZERO ) IB = 2
            END IF
!C
!C           Compute G, the current last IB rows of Z'*B.
!C
            NL = NSUP - IB + 1
            CALL DGEMM( 'Transpose', 'NoTranspose', IB, M, N, ONE, &
                       Z(1,NL), LDZ, B, LDB, ZERO, DWORK(KG), IB )
!C
!C           Check the controllability for a simple block.
!C
            IF( DLANGE( '1', IB, M, DWORK(KG), IB, DWORK(KW) )  &
                .LE. TOLERB ) THEN
!C
!C              Deflate the uncontrollable block and resume the
!C              main loop.
!C
               NSUP = NSUP - IB
               NUP = NUP + IB
               GO TO 20
            END IF
!C
!C           Test for termination with INFO = 3.
!C
            IF( NAP.EQ.NP ) THEN
               INFO = 3
!C
!C              Test for compatibility. Terminate if an attempt occurs
!C              to place a complex conjugate pair on a 1x1 block.
!C
            ELSE IF( IB.EQ.1 .AND. NPR.EQ.0 .AND. NLOW.EQ.NSUP ) THEN
               INFO = 4
            ELSE
!C
!C              Set the simple block flag.
!C
               SIMPLB = .TRUE.
!C
!C              Form a 2-by-2 block if necessary from two 1-by-1 blocks.
!C              Consider special case IB = 1, NPR = 1 and
!C              NPR+NPC > NSUP-NLOW+1 to avoid incompatibility.
!C
               IF( ( IB.EQ.1 .AND. NPR.EQ.0 ) .OR.  &
                  ( IB.EQ.1 .AND. NPR.EQ.1 .AND. NSUP.GT.NLOW .AND. &
                    NPR+NPC.GT.NSUP-NLOW+1 ) ) THEN
                  IF( NSUP.GT.2 ) THEN
                     IF( A(NSUP-1,NSUP-2).NE.ZERO ) THEN
!C
!C                       Interchange with the adjacent 2x2 block.
!C
!C                       Workspace needed: N.
!C
                        CALL DLAEXC( .TRUE., N, A, LDA, Z, LDZ, NSUP-2,  &
                                     2, 1, DWORK(KW), INFO )
                        IF( INFO.NE.0 ) THEN
                           INFO = 2
                           RETURN
                        END IF
                     ELSE
!C
!C                       Form a non-simple block by extending the last
!C                       block with a 1x1 block.
!C
                        SIMPLB = .FALSE.
                     END IF
                  ELSE
                     SIMPLB = .FALSE.
                  END IF
                  IB = 2
                  NL = NSUP - IB + 1
!C
!C                 Compute G, the current last IB rows of Z'*B.
!C
                  CALL DGEMM( 'Transpose', 'NoTranspose', IB, M, N, ONE, &
                             Z(1,NL), LDZ, B, LDB, ZERO, DWORK(KG), IB )
                            
!C
!C                 Check the controllability for the current block.
!C
                  IF( DLANGE( '1', IB, M, DWORK(KG), IB, DWORK(KW) ) &
                    .LE.TOLERB ) THEN
!C
!C                    Deflate the uncontrollable block and resume the
!C                    main loop.
!C
                     NSUP = NSUP - IB
                     NUP = NUP + IB
                     GO TO 20
                  END IF
               END IF
!C
               IF( NAP+IB.GT.NP ) THEN
!C
!C                 No sufficient eigenvalues to be assigned.
!C
                  INFO = 3
               ELSE
                  IF( IB.EQ.1 ) THEN
!C
!1C                    A 1-by-1 block.
!C
!C                    Assign the real eigenvalue nearest to A(NSUP,NSUP).
!C
                     X = A(NSUP,NSUP)
                     CALL SB01BX( .TRUE., NPR, X, X, WR, X, S, P )
                     NPR  = NPR - 1
                     CEIG = .FALSE.
                  ELSE
!C
!C                    A 2-by-2 block.
!C
                     IF( SIMPLB ) THEN
!C
!C                       Simple 2-by-2 block with complex eigenvalues.
!C                       Compute the eigenvalues of the last block.
!C
                        CALL MB03QY( N, NL, A, LDA, Z, LDZ, X, Y, INFO )
                        IF( NPC.GT.1 ) THEN
                           CALL SB01BX( .FALSE., NPC, X, Y, &
                                       WR(IPC), WI(IPC), S, P )
                           NPC  = NPC - 2
                           CEIG = .TRUE.
                        ELSE
!C
!C                          Choose the nearest two real eigenvalues.
!C
                           CALL SB01BX( .TRUE., NPR, X, X, WR, X, S, P )
                           CALL SB01BX( .TRUE., NPR-1, X, X, WR, X,  &
                                        Y, P )
                           P = S * Y
                           S = S + Y
                           NPR = NPR - 2
                           CEIG = .FALSE.
                        END IF
                     ELSE
!C
!C                       Non-simple 2x2 block with real eigenvalues.
!C                       Choose the nearest pair of complex eigenvalues.
!C
                        X = ( A(NL,NL) + A(NSUP,NSUP) )/TWO
                        CALL SB01BX( .FALSE., NPC, X, ZERO, WR(IPC), &
                                   WI(IPC), S, P )
                        NPC = NPC - 2
                     END IF
                  END IF
!C
!C                 Form the IBxIB matrix A2 from the current diagonal
!C                 block.
!!C
                  A2(1,1) = A(NL,NL)
                  IF( IB.GT.1 ) THEN
                     A2(1,2) = A(NL,NSUP)
                     A2(2,1) = A(NSUP,NL)
                     A2(2,2) = A(NSUP,NSUP)
                  END IF
!C
!C                 Determine the M-by-IB feedback matrix FI which
!C                 assigns the chosen IB eigenvalues for the pair (A2,G).
!C
!C                 Workspace needed: 5*M.
!C
                  CALL SB01BY( IB, M, S, P, A2, DWORK(KG), DWORK(KFI),  &
                              TOLER, DWORK(KW), IERR )
                  IF( IERR.NE.0 ) THEN
                     IF( IB.EQ.1 .OR. SIMPLB ) THEN
!C
!C                       The simple 1x1 block is uncontrollable.
!C
                        NSUP = NSUP - IB
                        IF( CEIG ) THEN
                           NPC = NPC + IB
                        ELSE
                           NPR = NPR + IB
                        END IF
                        NUP  = NUP + IB
                     ELSE
!C
!C                       The non-simple 2x2 block is uncontrollable.
!C                       Eliminate its uncontrollable part by using
!C                       the information in elements FI(1,1) and F(1,2).
!C
                        C = DWORK(KFI)
                        S = DWORK(KFI+IB)
!C
!C                       Apply the transformation to A and accumulate it
!C                       in Z.
!C
                        CALL DROT( N-NL+1, A(NL,NL), LDA,  &
                                  A(NSUP,NL), LDA, C, S )
                        CALL DROT( N, A(1,NL), 1, A(1,NSUP), 1, C, S )
                        CALL DROT( N, Z(1,NL), 1, Z(1,NSUP), 1, C, S )
!C
!C                       Annihilate the subdiagonal element of the last
!C                       block, redefine the upper limit for the bottom
!C                       block and resume the main loop.
!C
                        A(NSUP,NL) = ZERO
                        NSUP = NL
                        NUP  = NUP + 1
                        NPC  = NPC + 2
                     END IF
                  ELSE
!C
!C                    Successful assignment of IB eigenvalues.
!C
!C                    Update the feedback matrix F <-- F + [0 FI]*Z'.
!C
                     CALL DGEMM( 'NoTranspose', 'Transpose', M, N, &
                                IB, ONE, DWORK(KFI), M, Z(1,NL),   &
                                LDZ, ONE, F, LDF )
!C
!C                    Check for possible numerical instability.
!C
                     IF( DLANGE( '1', M, IB, DWORK(KFI), M, DWORK(KW) ) &
                                .GT. RMAX ) IWARN = IWARN + 1
!C
!C                    Update the state matrix A <-- A + Z'*B*[0 FI].
!C                    Workspace needed: 2*N+4*M.
!C
                     CALL DGEMM( 'NoTranspose', 'NoTranspose', N, IB, &
                                M, ONE, B, LDB, DWORK(KFI), M, ZERO,  &
                                DWORK(KW), N )
                     CALL DGEMM( 'Transpose', 'NoTranspose', NSUP,    & 
                                IB, N, ONE, Z, LDZ, DWORK(KW), N,    &
                                ONE, A(1,NL), LDA )
!C
!C                    Try to split the 2x2 block.
!C
                     IF( IB.EQ.2 ) &
                      CALL MB03QY( N, NL, A, LDA, Z, LDZ, X, Y, &
                                   INFO )
                     NAP = NAP + IB
                     IF( NLOW+IB.LE.NSUP ) THEN
!C
!C                       Move the last block(s) to the leading
!C                       position(s) of the bottom block.
!C
                        NCUR1 = NSUP - IB
                        NMOVES = 1
                        IF( IB.EQ.2 .AND. A(NSUP,NSUP-1).EQ.ZERO ) THEN
                           IB = 1
                           NMOVES = 2
                        END IF
!C
!C                       WHILE (NMOVES > 0) DO
   30                   IF( NMOVES.GT.0 ) THEN
                           NCUR = NCUR1
!C
!C                          WHILE (NCUR >= NLOW) DO
   40                      IF( NCUR.GE.NLOW ) THEN
!C
!C                             Loop for the last block positioning.
!C
                              IB1 = 1
                              IF( NCUR.GT.NLOW ) THEN
                                 IF( A(NCUR,NCUR-1).NE.ZERO ) IB1 = 2
                              END IF
                              CALL DLAEXC( .TRUE., N, A, LDA, Z, LDZ, &
                                          NCUR-IB1+1, IB1, IB,        &
                                          DWORK(KW), INFO )
                              IF( INFO.NE.0 ) THEN
                                 INFO = 2
                                 RETURN
                              END IF
                              NCUR = NCUR - IB1
                              GO TO 40
                           END IF
!C
!C                          END WHILE 40
!C
                           NMOVES = NMOVES - 1
                           NCUR1 = NCUR1 + 1
                           NLOW = NLOW + IB
                           GO TO 30
                        END IF
!C
!C                       END WHILE 30
!C
                     ELSE
                        NLOW = NLOW + IB
                     END IF
                  END IF
               END IF
            END IF
            IF( INFO.EQ.0 ) GO TO 20
!C
!C        END WHILE 20
!C
         END IF
!C
         WRKOPT = MAX( WRKOPT, 5*M, 2*N + 4*M )
      END IF
!C
!C     Annihilate the elements below the first subdiagonal of A.
!C
      IF( N.GT.2) &
        CALL DLASET( 'L', N-2, N-2, ZERO, ZERO, A(3,1), LDA )
      IF( NAP .GT. 0 ) THEN
!C
!C        Move the assigned eigenvalues in the first NAP positions of
!C        WR and WI.
!C
         K = IPC - NPR - 1
         IF( K.GT.0 ) THEN
            IF( K.LE.NPR ) THEN
               CALL DSWAP( K, WR(NPR+1), 1, WR, 1 )
            ELSE
               CALL DCOPY( K, WR(NPR+1), 1, DWORK, 1 )
               CALL DCOPY( NPR, WR, 1, DWORK(K+1), 1 )
               CALL DCOPY( K+NPR, DWORK, 1, WR, 1 )
            END IF
         END IF
         J = NAP - K
         IF( J.GT.0 ) THEN
            CALL DSWAP( J, WR(IPC+NPC), 1, WR(K+1), 1 )
            CALL DSWAP( J, WI(IPC+NPC), 1, WI(K+1), 1 )
         END IF
      END IF
!C
      DWORK(1) = WRKOPT
!C
   
END SUBROUTINE SB01BD

#if defined(__GFORTAN__) && (!defined(__ICC) || defined(__INTEL_COMPILER))
SUBROUTINE SB01BX(REIG,N,XR,XI,WR,WI,S,P) !GCC$ ATTRIBUTES INLINE :: SB01BX !GCC$ ATTRIBUTES ALIGNED(32) :: SB01BX
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE SB01BX(REIG,N,XR,XI,WR,WI,S,P)
 !DIR$ ATTRIBUTES FORCEINLINE :: SB01BX
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: SB01BX
!DIR$ OPTIMIZE : 3
!DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: SB01BX
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To choose a real eigenvalue or a pair of complex conjugate
C     eigenvalues at "minimal" distance to a given real or complex
C     value.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     REIG    LOGICAL
C             Specifies the type of eigenvalues as follows:
C             = .TRUE.,  a real eigenvalue is to be selected;
C             = .FALSE., a pair of complex eigenvalues is to be
C                        selected.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of eigenvalues contained in the arrays WR
C             and WI.  N >= 1.
C
C     XR,XI   (input) DOUBLE PRECISION
C             If REIG = .TRUE., XR must contain the real value and XI
C             is assumed zero and therefore not referenced.
C             If REIG = .FALSE., XR must contain the real part and XI
C             the imaginary part, respectively, of the complex value.
C
C     WR,WI   (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, if REIG = .TRUE., WR must contain the real
C             eigenvalues from which an eigenvalue at minimal distance
C             to XR is to be selected. In this case, WI is considered
C             zero and therefore not referenced.
C             On entry, if REIG = .FALSE., WR and WI must contain the
C             real and imaginary parts, respectively, of the eigenvalues
C             from which a pair of complex conjugate eigenvalues at
C             minimal "distance" to XR + jXI is to be selected.
C             The eigenvalues of each pair of complex conjugate
C             eigenvalues must appear consecutively.
C             On exit, the elements of these arrays are reordered such
C             that the selected eigenvalue(s) is (are) found in the
C             last element(s) of these arrays.
C
C     S,P     (output) DOUBLE PRECISION
C             If REIG = .TRUE., S (and also P) contains the value of
C             the selected real eigenvalue.
C             If REIG = .FALSE., S and P contain the sum and product,
C             respectively, of the selected complex conjugate pair of
C             eigenvalues.
C
C     FURTHER COMMENTS
C
C     For efficiency reasons, |x| + |y| is used for a complex number
C     x + jy, instead of its modulus.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     February 1999. Based on the RASP routine PMDIST.
C
C     REVISIONS
C
C     March 30, 1999, V. Sima, Research Institute for Informatics,
C     Bucharest.
C     Feb. 15, 2004, V. Sima, Research Institute for Informatics,
C     Bucharest.
C
C     ******************************************************************
C
#endif
!C     .. Scalar Arguments ..
      LOGICAL          REIG
      INTEGER          N
      DOUBLE PRECISION P, S, XI ,XR
      !C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION WI(*), WR(*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION WI(*), WR(*)
      !DIR$ ASSUME_ALIGNED WI:64
      !DIR$ ASSUME_ALIGNED WR:64
#endif
!C     .. Local Scalars ..
      INTEGER          I, J, K
      DOUBLE PRECISION X, Y
!C     .. Intrinsic Functions ..
      INTRINSIC        ABS
!C     .. Executable Statements ..
!C
      J = 1
      IF( REIG ) THEN
         Y = ABS( WR(1)-XR )
         DO 10 I = 2, N
            X = ABS( WR(I)-XR )
            IF( X .LT. Y ) THEN
               Y = X
               J = I
            END IF
   10    CONTINUE
         S = WR(J)
         K = N - J
         IF( K .GT. 0 ) THEN
            DO 20 I = J, J + K - 1
               WR(I) = WR(I+1)
   20       CONTINUE
            WR(N) = S
         END IF
         P = S
      ELSE
         Y = ABS( WR(1)-XR ) + ABS( WI(1)-XI )
         DO 30 I = 3, N, 2
            X = ABS( WR(I)-XR ) + ABS( WI(I)-XI )
            IF( X .LT. Y ) THEN
               Y = X
               J = I
            END IF
   30    CONTINUE
         X = WR(J)
         Y = WI(J)
         K = N - J - 1
         IF( K .GT. 0 ) THEN
            DO 40 I = J, J + K - 1
               WR(I) = WR(I+2)
               WI(I) = WI(I+2)
   40       CONTINUE
            WR(N-1) = X
            WI(N-1) = Y
            WR(N) = X
            WI(N) = -Y
         END IF
         S = X + X
         P = X * X + Y * Y
      END IF
     
END SUBROUTINE

#if defined(__GFORTAN__) && (!defined(__ICC) || defined(__INTEL_COMPILER))
SUBROUTINE SB01BY(N,M,S,P,A B,F,TOL,DWORK,INFO) !GCC$ ATTRIBUTES hot :: SB01BY !GCC$ ATTRIBUTES ALIGNED(32) :: SB01BY !GCC$ ATTRIBUTES no_stack_protector :: SB01BY
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE SB01BY(N,M,S,P,A B,F,TOL,DWORK,INFO)
 !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: SB01BY
!DIR$ OPTIMIZE : 3
!DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: SB01BY
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To solve an N-by-N pole placement problem for the simple cases
C     N = 1 or N = 2: given the N-by-N matrix A and N-by-M matrix B,
C     construct an M-by-N matrix F such that A + B*F has prescribed
C     eigenvalues. These eigenvalues are specified by their sum S and
C     product P (if N = 2). The resulting F has minimum Frobenius norm.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A and also the number of rows of
C             the matrix B and the number of columns of the matrix F.
C             N is either 1, if a single real eigenvalue is prescribed
C             or 2, if a complex conjugate pair or a set of two real
C             eigenvalues are prescribed.
C
C     M       (input) INTEGER
C             The number of columns of the matrix B and also the number
C             of rows of the matrix F.  M >= 1.
C
C     S       (input) DOUBLE PRECISION
C             The sum of the prescribed eigenvalues if N = 2 or the
C             value of prescribed eigenvalue if N = 1.
C
C     P       (input) DOUBLE PRECISION
C             The product of the prescribed eigenvalues if N = 2.
C             Not referenced if N = 1.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (N,N)
C             On entry, this array must contain the N-by-N state
C             dynamics matrix whose eigenvalues have to be moved to
C             prescribed locations.
C             On exit, this array contains no useful information.
C
C     B       (input/output) DOUBLE PRECISION array, dimension (N,M)
C             On entry, this array must contain the N-by-M input/state
C             matrix B.
C             On exit, this array contains no useful information.
C
C     F       (output) DOUBLE PRECISION array, dimension (M,N)
C             The state feedback matrix F which assigns one pole or two
C             poles of the closed-loop matrix A + B*F.
C             If N = 2 and the pair (A,B) is not controllable
C             (INFO = 1), then F(1,1) and F(1,2) contain the elements of
C             an orthogonal rotation which can be used to remove the
C             uncontrollable part of the pair (A,B).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The absolute tolerance level below which the elements of A
C             and B are considered zero (used for controllability test).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (M)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 1:  if uncontrollability of the pair (A,B) is detected.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, July 1998.
C     Based on the RASP routine SB01BY.
C
C     REVISIONS
C
C     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest.
C     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven.
C     May  2003, A. Varga, German Aerospace Center.
C
C     KEYWORDS
C
C     Eigenvalue, eigenvalue assignment, feedback control, pole
C     placement, state-space model.
C
C     ******************************************************************
C
#endif
      implicit none
!C     .. Parameters ..
      DOUBLE PRECISION  FOUR, ONE, THREE, ZERO
      PARAMETER         ( FOUR = 4.0D0,  ONE = 1.0D0, THREE = 3.0D0, &
                         ZERO = 0.0D0 )
!C     .. Scalar Arguments ..
      INTEGER           INFO, M, N
      DOUBLE PRECISION  P, S, TOL
!C     .. Array Arguments ..
      DOUBLE PRECISION  A(N,*), B(N,*), DWORK(*), F(M,*)
!C     .. Local Scalars ..
      INTEGER           IR, J
      DOUBLE PRECISION  ABSR, B1, B2, B21, C, C0, C1, C11, C12, C21,    &
                       C22, C3, C4, CS, CU, CV, DC0, DC2, DC3, DIFFR,  &
                       R, RN, S12, S21, SIG, SN, SU, SV, TAU1, TAU2,   &
                       WI, WI1, WR, WR1, X, Y, Z
!C     .. External Functions ..
      DOUBLE PRECISION  DLAMC3
      EXTERNAL          DLAMC3
!C     .. External Subroutines ..
      EXTERNAL          DLANV2, DLARFG, DLASET, DLASV2, DLATZM, DROT
!C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
!C     .. Executable Statements ..
!C
!C     For efficiency reasons, the parameters are not checked.
!C
      INFO = 0
      IF( N.EQ.1 ) THEN
!C
!C        The case N = 1.
!C
         IF( M.GT.1 ) &
            CALL DLARFG( M, B(1,1), B(1,2), N, TAU1 )
         B1 = B(1,1)
         IF( ABS( B1 ).LE.TOL ) THEN
!C
!C           The pair (A,B) is uncontrollable.
!C
            INFO = 1
            RETURN
         END IF
!C
         F(1,1) = ( S - A(1,1) )/B1
         IF( M.GT.1 ) THEN
            CALL DLASET( 'Full', M-1, 1, ZERO, ZERO, F(2,1), M )
            CALL DLATZM( 'Left', M, N, B(1,2), N, TAU1, F(1,1), F(2,1), &
                         M, DWORK )
         END IF
         RETURN
      END IF
!C
!C     In the sequel N = 2.
!C
!C     Compute the singular value decomposition of B in the form
!C
!C                    ( V  0 )                ( B1 0  )
!C     B = U*( G1 0 )*(      )*H2*H1 ,   G1 = (       ),
!C                    ( 0  I )                ( 0  B2 )
!C
!C               ( CU   SU )          ( CV   SV )
!C     where U = (         )  and V = (         )  are orthogonal
!C               (-SU   CU )          (-SV   CV )
!C
!C     rotations and H1 and H2 are elementary Householder reflectors.
!C     ABS(B1) and ABS(B2) are the singular values of matrix B,
!C     with ABS(B1) >= ABS(B2).
!C
!C     Reduce first B to the lower bidiagonal form  ( B1  0  ... 0 ).
!C                                                  ( B21 B2 ... 0 )
      IF( M.EQ.1 ) THEN
!C
!C        Initialization for the case M = 1; no reduction required.
!C
         B1  = B(1,1)
         B21 = B(2,1)
         B2  = ZERO
      ELSE
!C
!C        Postmultiply B with elementary Householder reflectors H1
!C        and H2.
!C
         CALL DLARFG( M, B(1,1), B(1,2), N, TAU1 )
         CALL DLATZM( 'Right', N-1, M, B(1,2), N, TAU1, B(2,1), B(2,2),  &
                      N, DWORK )
         B1  = B(1,1)
         B21 = B(2,1)
         IF( M.GT.2 ) &
            CALL DLARFG( M-1, B(2,2), B(2,3), N, TAU2 )
         B2  = B(2,2)
      END IF
!C
!C     Reduce B to a diagonal form by premultiplying and postmultiplying
!C     it with orthogonal rotations U and V, respectively, and order the
!C     diagonal elements to have decreasing magnitudes.
!C     Note: B2 has been set to zero if M = 1. Thus in the following
!1C     computations the case M = 1 need not to be distinguished.
!1C     Note also that LAPACK routine DLASV2 assumes an upper triangular
!C     matrix, so the results should be adapted.
!C
      CALL DLASV2( B1, B21, B2, X, Y, SU, CU, SV, CV )
      SU = -SU
      B1 =  Y
      B2 =  X
!1C
!C     Compute  A1 = U'*A*U.
!C
      CALL DROT( 2, A(2,1), 2, A(1,1), 2, CU, SU )
      CALL DROT( 2, A(1,2), 1, A(1,1), 1, CU, SU )
!C
!C     Compute the rank of B and check the controllability of the
!C     pair (A,B).
!C
      IR = 0
      IF( ABS( B2 ).GT.TOL ) IR = IR + 1
      IF( ABS( B1 ).GT.TOL ) IR = IR + 1
      IF( IR.EQ.0 .OR. ( IR.EQ.1 .AND. ABS( A(2,1) ).LE.TOL ) ) THEN
         F(1,1) =  CU
         F(1,2) = -SU
!C
!C        The pair (A,B) is uncontrollable.
!C
         INFO = 1
         RETURN
      END IF
!C
!C     Compute F1 which assigns N poles for the reduced pair (A1,G1).
!C
      X = DLAMC3( B1, B2 )
      IF( X.EQ.B1 ) THEN
!C
!C        Rank one G1.
!C
         F(1,1) = ( S - ( A(1,1) + A(2,2) ) )/B1
         F(1,2) = -( A(2,2)*( A(2,2) - S ) + A(2,1)*A(1,2) + P )/ &
                   A(2,1)/B1
         IF( M.GT.1 ) THEN
            F(2,1) = ZERO
            F(2,2) = ZERO
         END IF
      ELSE
!C
!C        Rank two G1.
!C
         Z = ( S - ( A(1,1) + A(2,2) ) )/( B1*B1 + B2*B2 )
         F(1,1) = B1*Z
         F(2,2) = B2*Z
!C
!C        Compute an approximation for the minimum norm parameter
!C        selection.
!C
         X = A(1,1) + B1*F(1,1)
         C = X*( S - X ) - P
         IF( C.GE.ZERO ) THEN
            SIG =  ONE
         ELSE
            SIG = -ONE
         END IF
         S12 = B1/B2
         S21 = B2/B1
         C11 = ZERO
         C12 = ONE
         C21 = SIG*S12*C
         C22 = A(1,2) - SIG*S12*A(2,1)
         CALL DLANV2( C11, C12, C21, C22, WR, WI, WR1, WI1, CS, SN )
         IF( ABS( WR - A(1,2) ).GT.ABS( WR1 - A(1,2) ) ) THEN
            R = WR1
         ELSE
            R = WR
         END IF
!C
!C        Perform Newton iteration to solve the equation for minimum.
!C
         C0 = -C*C
         C1 =  C*A(2,1)
         C4 =  S21*S21
         C3 = -C4*A(1,2)
         DC0 = C1
         DC2 = THREE*C3
         DC3 = FOUR*C4
!C
         DO 10 J = 1, 10
            X  = C0 + R*( C1 + R*R*( C3 + R*C4 ) )
            Y  = DC0 + R*R*( DC2 + R*DC3 )
            IF( Y.EQ.ZERO ) EXIT
            RN = R - X/Y
            ABSR  = ABS( R )
            DIFFR = ABS( R - RN )
            Z = DLAMC3( ABSR, DIFFR )
            IF( Z.EQ.ABSR ) &
              EXIT
            R = RN
   10    CONTINUE
!C
   20    CONTINUE
         IF( R.EQ.ZERO ) R = DLAMCH( 'Epsilon' )
         F(1,2) = (  R  - A(1,2) )/B1
         F(2,1) = ( C/R - A(2,1) )/B2
      END IF
!C
!C     Back-transform F1. Compute first F1*U'.
!C
      CALL DROT( MIN( M, 2 ), F(1,1), 1, F(1,2), 1, CU, SU )
      IF( M.EQ.1 ) &
          RETURN
!C
!C     Compute V'*F1.
!C
      CALL DROT( 2, F(2,1), M, F(1,1), M, CV, SV )
!C
!C               ( F1 )
!C     Form  F = (    ) .
!C               ( 0  )
!C
      IF( M.GT.N ) &
         CALL DLASET( 'Full', M-N, N, ZERO, ZERO, F(N+1,1), M )
!C
!C     Compute H1*H2*F.
!C
      IF( M.GT.2 ) &
        CALL DLATZM( 'Left', M-1, N, B(2,3), N, TAU2, F(2,1), F(3,1), &
                     M, DWORK )
      CALL DLATZM( 'Left', M, N, B(1,2), N, TAU1, F(1,1), F(2,1), M,  &
                   DWORK )

END SUBROUTINE

    
    


#if defined(__GFORTRAN__) && (!defined(__ICC) || defined(__INTEL_COMPILER))
SUBROUTINE MB03QD( DICO, STDOM, JOBU, N, NLOW, NSUP, ALPHA, &
     A, LDA, U, LDU, NDIM, DWORK, INFO ) !GCC$ ATTRIBUTES hot :: MB03QD !GCC$ ATTRIBUTES aligned(32) :: MB03QD !GCC$ ATTRIBUTES no_stack_protector :: MB03QD
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE MB03QD( DICO, STDOM, JOBU, N, NLOW, NSUP, ALPHA, &
       A, LDA, U, LDU, NDIM, DWORK, INFO )
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: MB03QD
!DIR$ OPTIMIZE : 3
!DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: MB03QD
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To reorder the diagonal blocks of a principal submatrix of an
C     upper quasi-triangular matrix A together with their eigenvalues by
C     constructing an orthogonal similarity transformation UT.
C     After reordering, the leading block of the selected submatrix of A
C     has eigenvalues in a suitably defined domain of interest, usually
C     related to stability/instability in a continuous- or discrete-time
C     sense.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the spectrum separation to be
C             performed as follows:
C             = 'C':  continuous-time sense;
C             = 'D':  discrete-time sense.
C
C     STDOM   CHARACTER*1
C             Specifies whether the domain of interest is of stability
C             type (left part of complex plane or inside of a circle)
C             or of instability type (right part of complex plane or
C             outside of a circle) as follows:
C             = 'S':  stability type domain;
C             = 'U':  instability type domain.
C
C     JOBU    CHARACTER*1
C             Indicates how the performed orthogonal transformations UT
C             are accumulated, as follows:
C             = 'I':  U is initialized to the unit matrix and the matrix
C                     UT is returned in U;
C             = 'U':  the given matrix U is updated and the matrix U*UT
C                     is returned in U.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A and U.  N >= 1.
C
C     NLOW,   (input) INTEGER
C     NSUP    NLOW and NSUP specify the boundary indices for the rows
C             and columns of the principal submatrix of A whose diagonal
C             blocks are to be reordered.  1 <= NLOW <= NSUP <= N.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The boundary of the domain of interest for the eigenvalues
C             of A. If DICO = 'C', ALPHA is the boundary value for the
C             real parts of eigenvalues, while for DICO = 'D',
C             ALPHA >= 0 represents the boundary value for the moduli of
C             eigenvalues.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain a matrix in a real Schur form whose 1-by-1 and
C             2-by-2 diagonal blocks between positions NLOW and NSUP
C             are to be reordered.
C             On exit, the leading N-by-N part contains the ordered
!C             real Schur matrix UT' * A * UT with the elements below the
C             first subdiagonal set to zero.
C             The leading NDIM-by-NDIM part of the principal submatrix
C             D = A(NLOW:NSUP,NLOW:NSUP) has eigenvalues in the domain
C             of interest and the trailing part of this submatrix has
C             eigenvalues outside the domain of interest.
C             The domain of interest for lambda(D), the eigenvalues of
C             D, is defined by the parameters ALPHA, DICO and STDOM as
C             follows:
C               For DICO = 'C':
C                  Real(lambda(D)) < ALPHA if STDOM = 'S';
C                  Real(lambda(D)) > ALPHA if STDOM = 'U'.
C               For DICO = 'D':
C                  Abs(lambda(D)) < ALPHA if STDOM = 'S';
C                  Abs(lambda(D)) > ALPHA if STDOM = 'U'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= N.
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,N)
C             On entry with JOBU = 'U', the leading N-by-N part of this
C             array must contain a transformation matrix (e.g. from a
C             previous call to this routine).
C             On exit, if JOBU = 'U', the leading N-by-N part of this
C             array contains the product of the input matrix U and the
C             orthogonal matrix UT used to reorder the diagonal blocks
C             of A.
C             On exit, if JOBU = 'I', the leading N-by-N part of this
C             array contains the matrix UT of the performed orthogonal
C             transformations.
C             Array U need not be set on entry if JOBU = 'I'.
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= N.
C
C     NDIM    (output) INTEGER
C             The number of eigenvalues of the selected principal
C             submatrix lying inside the domain of interest.
C             If NLOW = 1, NDIM is also the dimension of the invariant
C             subspace corresponding to the eigenvalues of the leading
C             NDIM-by-NDIM submatrix. In this case, if U is the
C             orthogonal transformation matrix used to compute and
C             reorder the real Schur form of A, its first NDIM columns
C             form an orthonormal basis for the above invariant
C             subspace.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  A(NLOW,NLOW-1) is nonzero, i.e. A(NLOW,NLOW) is not
C                   the leading element of a 1-by-1 or 2-by-2 diagonal
C                   block of A, or A(NSUP+1,NSUP) is nonzero, i.e.
C                   A(NSUP,NSUP) is not the bottom element of a 1-by-1
C                   or 2-by-2 diagonal block of A;
C             = 2:  two adjacent blocks are too close to swap (the
C                   problem is very ill-conditioned).
C
C     METHOD
C
C     Given an upper quasi-triangular matrix A with 1-by-1 or 2-by-2
C     diagonal blocks, the routine reorders its diagonal blocks along
C     with its eigenvalues by performing an orthogonal similarity
!C     transformation UT' * A * UT. The column transformation UT is also
C     performed on the given (initial) transformation U (resulted from
C     a possible previous step or initialized as the identity matrix).
C     After reordering, the eigenvalues inside the region specified by
C     the parameters ALPHA, DICO and STDOM appear at the top of
C     the selected diagonal block between positions NLOW and NSUP.
C     In other words, lambda(A(NLOW:NSUP,NLOW:NSUP)) are ordered such
C     that lambda(A(NLOW:NLOW+NDIM-1,NLOW:NLOW+NDIM-1)) are inside and
C     lambda(A(NLOW+NDIM:NSUP,NLOW+NDIM:NSUP)) are outside the domain
C     of interest. If NLOW = 1, the first NDIM columns of U*UT span the
C     corresponding invariant subspace of A.
C
C     REFERENCES
C
C     [1] Stewart, G.W.
C         HQR3 and EXCHQZ: FORTRAN subroutines for calculating and
C         ordering the eigenvalues of a real upper Hessenberg matrix.
C         ACM TOMS, 2, pp. 275-280, 1976.
C
C     NUMERICAL ASPECTS
C                                         3
C     The algorithm requires less than 4*N  operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen,
C     April 1998. Based on the RASP routine SEOR1.
C
C     KEYWORDS
C
C     Eigenvalues, invariant subspace, orthogonal transformation, real
C     Schur form, similarity transformation.
C
C    ******************************************************************
C
#endif
!C     .. Parameters ..
      DOUBLE PRECISION ONE, ZERO
      PARAMETER        ( ONE = 1.0D0, ZERO = 0.0D0 )
!C     .. Scalar Arguments ..
      CHARACTER        DICO, JOBU, STDOM
      INTEGER          INFO, LDA, LDU, N, NDIM, NLOW, NSUP
      DOUBLE PRECISION ALPHA
      !C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION A(LDA,*), DWORK(*), U(LDU,*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION A(LDA,*), DWORK(*), U(LDU,*)
      !DIR$ ASSUME_ALIGNED A:64
      !DIR$ ASSUME_ALIGNED DWORK:64
      !DIR$ ASSUME_ALIGNED U:64
#endif
!C     .. Local Scalars ..
      LOGICAL          DISCR, LSTDOM
      INTEGER          IB, L, LM1, NUP
      DOUBLE PRECISION E1, E2, TLAMBD
!C     .. External Functions ..
     
      DOUBLE PRECISION DLAPY2
      EXTERNAL         DLAPY2
!C     .. External Subroutines ..
      EXTERNAL         DLASET, DTREXC
!C     .. Intrinsic Functions ..
      INTRINSIC        ABS
!C     .. Executable Statements ..
!C
      INFO = 0
      DISCR = LSAME( DICO, 'D' )
      LSTDOM = LSAME( STDOM, 'S' )
!C
!C     Check input scalar arguments.
!C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSTDOM .OR. LSAME( STDOM, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( LSAME( JOBU, 'I' ) .OR.  &
                      LSAME( JOBU, 'U' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.1 ) THEN
         INFO = -4
      ELSE IF( NLOW.LT.1 ) THEN
         INFO = -5
      ELSE IF( NLOW.GT.NSUP .OR. NSUP.GT.N ) THEN
         INFO = -6
      ELSE IF( DISCR .AND. ALPHA.LT.ZERO ) THEN
         INFO = -7
      ELSE IF( LDA.LT.N ) THEN
         INFO = -9
      ELSE IF( LDU.LT.N ) THEN
         INFO = -11
      END IF
!C
      IF( INFO.NE.0 ) THEN
!C
!C        Error return.
!C
          RETURN
      END IF

      IF( NLOW.GT.1 ) THEN
         IF( A(NLOW,NLOW-1).NE.ZERO ) INFO = 1
      END IF
      IF( NSUP.LT.N ) THEN
         IF( A(NSUP+1,NSUP).NE.ZERO ) INFO = 1
      END IF
      IF( INFO.NE.0 ) &
          RETURN
!C
!C     Initialize U with an identity matrix if necessary.
!C
      IF( LSAME( JOBU, 'I' ) ) &
        CALL DLASET( 'Full', N, N, ZERO, ONE, U, LDU )

      NDIM = 0
      L = NSUP
      NUP = NSUP
!C
!C     NUP is the minimal value such that the submatrix A(i,j) with
!C     NUP+1 <= i,j <= NSUP contains no eigenvalues inside the domain of
!C     interest. L is such that all the eigenvalues of the submatrix
!C     A(i,j) with L+1 <= i,j <= NUP lie inside the domain of interest.
!C
!C     WHILE( L >= NLOW ) DO
!C
   10 IF( L.GE.NLOW ) THEN
         IB = 1
         IF( L.GT.NLOW ) THEN
            LM1 = L - 1
            IF( A(L,LM1).NE.ZERO ) THEN
               CALL MB03QY( N, LM1, A, LDA, U, LDU, E1, E2, INFO )
               IF( A(L,LM1).NE.ZERO ) IB = 2
            END IF
         END IF
         IF( DISCR ) THEN
            IF( IB.EQ.1 ) THEN
               TLAMBD = ABS( A(L,L) )
            ELSE
               TLAMBD = DLAPY2( E1, E2 )
            END IF
         ELSE
            IF( IB.EQ.1 ) THEN
               TLAMBD = A(L,L)
            ELSE
               TLAMBD = E1
            END IF
         END IF
         IF( (      LSTDOM .AND. TLAMBD.LT.ALPHA ) .OR.    &
             ( .NOT.LSTDOM .AND. TLAMBD.GT.ALPHA ) ) THEN
            NDIM = NDIM + IB
            L = L - IB
         ELSE
            IF( NDIM.NE.0 ) THEN
               CALL DTREXC( 'V', N, A, LDA, U, LDU, L, NUP, DWORK,  &
                           INFO )
               IF( INFO.NE.0 ) THEN
                  INFO = 2
                  RETURN
               END IF
               NUP = NUP - 1
               L = L - 1
            ELSE
               NUP = NUP - IB
               L = L - IB
            END IF
         END IF
         GO TO 10
      END IF

END SUBROUTINE 

#if defined(__GFORTRAN__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
SUBROUTINE MB03QY(N,L,A,LDA,U,LDU,E1,E2,INFO) !GCC$ ATTRIBUTES INLINE :: MB03QY !GCC$ ATTRIBUTES ALIGNED(32) :: MB03QY
#elif defined(__ICC) || defined(__INTEL_COMPILER)
  SUBROUTINE MB03QY(N,L,A,LDA,U,LDU,E1,E2,INFO)
    !DIR$ ATTRIBUTES FORCEINLINE :: MB03QY
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: MB03QY
!DIR$ OPTIMIZE : 3
!DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=Haswell :: MB03QY
#endif
#if 0
C
C     SLICOT RELEASE 5.7.
C
C     Copyright (c) 2002-2020 NICONET e.V.
C
C     PURPOSE
C
C     To compute the eigenvalues of a selected 2-by-2 diagonal block
C     of an upper quasi-triangular matrix, to reduce the selected block
C     to the standard form and to split the block in the case of real
C     eigenvalues by constructing an orthogonal transformation UT.
C     This transformation is applied to A (by similarity) and to
C     another matrix U from the right.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A and UT.  N >= 2.
C
C     L       (input) INTEGER
C             Specifies the position of the block.  1 <= L < N.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper quasi-triangular matrix A whose
C             selected 2-by-2 diagonal block is to be processed.
C             On exit, the leading N-by-N part of this array contains
C             the upper quasi-triangular matrix A after its selected
C             block has been splitt and/or put in the LAPACK standard
C             form.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= N.
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,N)
C             On entry, the leading N-by-N part of this array must
C             contain a transformation matrix U.
C             On exit, the leading N-by-N part of this array contains
C             U*UT, where UT is the transformation matrix used to
C             split and/or standardize the selected block.
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= N.
C
C     E1, E2  (output) DOUBLE PRECISION
C             E1 and E2 contain either the real eigenvalues or the real
C             and positive imaginary parts, respectively, of the complex
C             eigenvalues of the selected 2-by-2 diagonal block of A.
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
C     Let A1 = ( A(L,L)    A(L,L+1)   )
C              ( A(L+1,L)  A(L+1,L+1) )
C     be the specified 2-by-2 diagonal block of matrix A.
C     If the eigenvalues of A1 are complex, then they are computed and
C     stored in E1 and E2, where the real part is stored in E1 and the
C     positive imaginary part in E2. The 2-by-2 block is reduced if
C     necessary to the standard form, such that A(L,L) = A(L+1,L+1), and
C     A(L,L+1) and A(L+1,L) have oposite signs. If the eigenvalues are
C     real, the 2-by-2 block is reduced to an upper triangular form such
C     that ABS(A(L,L)) >= ABS(A(L+1,L+1)).
C     In both cases, an orthogonal rotation U1' is constructed such that
C     U1'*A1*U1 has the appropriate form. Let UT be an extension of U1
!C     to an N-by-N orthogonal matrix, using identity submatrices. Then A
!C     is replaced by UT'*A*UT and the contents of array U is U * UT.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen,
C     March 1998. Based on the RASP routine SPLITB.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Eigenvalues, orthogonal transformation, real Schur form,
C     similarity transformation.
C
C     ******************************************************************
C
#endif
       implicit none
!C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER        ( ZERO = 0.0D0 )
!C     .. Scalar Arguments ..
      INTEGER          INFO, L, LDA, LDU, N
      DOUBLE PRECISION E1, E2
      !C     .. Array Arguments ..
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
      DOUBLE PRECISION A(LDA,*), U(LDU,*)
#elif defined(__ICC) || defined(__INTEL_COMPILER)
      DOUBLE PRECISION A(LDA,*), U(LDU,*)
      !DIR$ ASSUME_ALIGNED A:64
      !DIR$ ASSUME_ALIGNED U:64
#endif
!C     .. Local Scalars ..
      INTEGER          L1
      DOUBLE PRECISION EW1, EW2, CS, SN
!C     .. External Subroutines ..
      EXTERNAL         DLANV2, DROT
!1C     .. Executable Statements ..
!C
      INFO = 0
!C
!C     Test the input scalar arguments.
!C
      IF( N.LT.2 ) THEN
         INFO = -1
      ELSE IF( L.LT.1 .OR. L.GE.N ) THEN
         INFO = -2
      ELSE IF( LDA.LT.N ) THEN
         INFO = -4
      ELSE IF( LDU.LT.N ) THEN
         INFO = -6
      END IF

      IF( INFO.NE.0 ) THEN
!C
!C        Error return.
!C
          RETURN
      END IF
!C
!C     Compute the eigenvalues and the elements of the Givens
!C     transformation.
!C
      L1 = L + 1
      CALL DLANV2( A(L,L), A(L,L1), A(L1,L), A(L1,L1), E1, E2, &
                   EW1, EW2, CS, SN )
      IF( E2.EQ.ZERO ) E2 = EW1
!C
!C     Apply the transformation to A.
!C
      IF( L1.LT.N ) &
         CALL DROT( N-L1, A(L,L1+1), LDA, A(L1,L1+1), LDA, CS, SN )
      CALL DROT( L-1, A(1,L), 1, A(1,L1), 1, CS, SN )
!C
!C     Accumulate the transformation in U.
!C
      CALL DROT( N, U(1,L), 1, U(1,L1), 1, CS, SN )
!C
 
END SUBROUTINE MB03QY





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

