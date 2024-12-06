C
C     file sepx4.f
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 2005 by UCAR                   *
C     *                                                               *
C     *       University Corporation for Atmospheric Research         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                    FISHPACK90  version 1.1                    *
C     *                                                               *
C     *                 A Package of Fortran 77 and 90                *
C     *                                                               *
C     *                Subroutines and Example Programs               *
C     *                                                               *
C     *               for Modeling Geophysical Processes              *
C     *                                                               *
C     *                             by                                *
C     *                                                               *
C     *        John Adams, Paul Swarztrauber and Roland Sweet         *
C     *                                                               *
C     *                             of                                *
C     *                                                               *
C     *         the National Center for Atmospheric Research          *
C     *                                                               *
C     *                Boulder, Colorado  (80307)  U.S.A.             *
C     *                                                               *
C     *                   which is sponsored by                       *
C     *                                                               *
C     *              the National Science Foundation                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE SEPX4 (IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
C    +                  NBDCND,BDC,BDD,COFX,GRHS,USOL,IDMN,PERTRB,
C    +                  IERROR)
C
C
C
C DIMENSION OF           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
C ARGUMENTS              USOL(IDMN,N+1),     GRHS(IDMN,N+1),
C
C
C LATEST REVISION        June 2004
C
C PURPOSE                SEPX4 SOLVES FOR EITHER THE SECOND-ORDER
C                        FINITE DIFFERENCE APPROXIMATION OR A
C                        FOURTH-ORDER APPROXIMATION TO A SEPARABLE
C                        ELLIPTIC EQUATION
C
C                          AF(X)*UXX+BF(X)*UX+CF(X)*U+UYY = G(X,Y)
C
C                        ON A RECTANGLE (X GREATER THAN OR EQUAL TO
C                        A AND LESS THAN OR EQUAL TO B, Y GREATER THAN
C                        OR EQUAL TO C AND LESS THAN OR EQUAL TO D).
C                        ANY COMBINATION OF PERIODIC OR MIXED BOUNDARY
C                        CONDITIONS IS ALLOWED.  IF BOUNDARY
C                        CONDITIONS IN THE X DIRECTION ARE PERIODIC
C                        (SEE MBDCND=0 BELOW) THEN THE COEFFICIENTS
C                        MUST SATISFY
C
C                          AF(X)=C1,BF(X)=0,CF(X)=C2 FOR ALL X.
C
C                        HERE C1,C2 ARE CONSTANTS, C1.GT.0.
C
C                        THE POSSIBLE BOUNDARY CONDITIONS ARE:
C                        IN THE X-DIRECTION:
C                          (0) PERIODIC, U(X+B-A,Y)=U(X,Y) FOR
C                              ALL Y,X
C                          (1) U(A,Y), U(B,Y) ARE SPECIFIED FOR ALL Y
C                          (2) U(A,Y), DU(B,Y)/DX+BETA*U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C                          (3) DU(A,Y)/DX+ALPHA*U(A,Y),DU(B,Y)/DX+
C                              BETA*U(B,Y) ARE SPECIFIED FOR ALL Y
C                          (4) DU(A,Y)/DX+ALPHA*U(A,Y),U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C
C                        IN THE Y-DIRECTION:
C                          (0) PERIODIC, U(X,Y+D-C)=U(X,Y) FOR ALL X,Y
C                          (1) U(X,C),U(X,D) ARE SPECIFIED FOR ALL X
C                          (2) U(X,C),DU(X,D)/DY ARE SPECIFIED FOR
C                              ALL X
C                          (3) DU(X,C)/DY,DU(X,D)/DY ARE SPECIFIED FOR
C                              ALL X
C                          (4) DU(X,C)/DY,U(X,D) ARE SPECIFIED FOR
C                              ALL X
C
C USAGE                  CALL SEPX4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,
C                                   BETA,C,D,N,NBDCND,BDC,BDD,COFX,
C                                   GRHS,USOL,IDMN,W,PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               IORDER
C                          = 2 IF A SECOND-ORDER APPROXIMATION IS
C                              SOUGHT
C                          = 4 IF A FOURTH-ORDER APPROXIMATION IS
C                              SOUGHT
C
c *** caution ***          GRHS SHOULD BE RESET IF SEPX4 WAS FIRST CALLED
C                          WITH IORDER=2 AND WILL BE CALLED AGAIN WITH
C                          IORDER=4.  VALUES IN GRHS ARE DESTROYED BY THE
C                          IORDER=2 CALL.
C
C
C                        A,B
C                          THE RANGE OF THE X-INDEPENDENT VARIABLE,
C                          I.E., X IS GREATER THAN OR EQUAL TO A
C                          AND LESS THAN OR EQUAL TO B.  A MUST BE
C                          LESS THAN B.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (A,B) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE M+1 GRID POINTS IN THE X-
C                          DIRECTION GIVEN BY XI=A+(I-1)*DLX
C                          FOR I=1,2,...,M+1 WHERE DLX=(B-A)/M IS
C                          THE PANEL WIDTH.  M MUST BE LESS THAN
C                          IDMN AND GREATER THAN 5.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT X=A AND X=B
C                          = 0 IF THE SOLUTION IS PERIODIC IN X, I.E.,
C                              U(X+B-A,Y)=U(X,Y)  FOR ALL Y,X
C                          = 1 IF THE SOLUTION IS SPECIFIED AT X=A
C                              AND X=B, I.E., U(A,Y) AND U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C                          = 2 IF THE SOLUTION IS SPECIFIED AT X=A
C                              AND THE BOUNDARY CONDITION IS MIXED AT
C                              X=B, I.E., U(A,Y) AND
C                              DU(B,Y)/DX+BETA*U(B,Y) ARE SPECIFIED
C                              FOR ALL Y
C                          = 3 IF THE BOUNDARY CONDITIONS AT X=A AND
C                              X=B ARE MIXED, I.E.,
C                              DU(A,Y)/DX+ALPHA*U(A,Y) AND
C                              DU(B,Y)/DX+BETA*U(B,Y) ARE SPECIFIED
C                              FOR ALL Y
C                          = 4 IF THE BOUNDARY CONDITION AT X=A IS
C                              MIXED AND THE SOLUTION IS SPECIFIED
C                              AT X=B, I.E., DU(A,Y)/DX+ALPHA*U(A,Y)
C                              AND U(B,Y) ARE SPECIFIED FOR ALL Y
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF
C                          DU(A,Y)/DX+ ALPHA*U(A,Y) AT X=A, WHEN
C                          MBDCND=3 OR 4.
C                          BDA(J) = DU(A,YJ)/DX+ALPHA*U(A,YJ),
C                          J=1,2,...,N+1
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
C                          A DUMMY PARAMETER.
C
C                        ALPHA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN CASE
C                          OF A MIXED BOUNDARY CONDITION AT X=A
C                          (SEE ARGUMENT BDA).  IF MBDCND IS NOT EQUAL
C                          TO EITHER 3 OR 4, THEN ALPHA IS A DUMMY
C                          PARAMETER.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF
C                          DU(B,Y)/DX+ BETA*U(B,Y) AT X=B.
C                          WHEN MBDCND=2 OR 3
C                          BDB(J) = DU(B,YJ)/DX+BETA*U(B,YJ),
C                          J=1,2,...,N+1
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
C                          A DUMMY PARAMETER.
C
C                        BETA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT X=B
C                          (SEE ARGUMENT BDB).  IF MBDCND IS NOT EQUAL
C                          TO 2 OR 3, THEN BETA IS A DUMMY PARAMETER.
C
C                        C,D
C                          THE RANGE OF THE Y-INDEPENDENT VARIABLE,
C                          I.E., Y IS GREATER THAN OR EQUAL TO C AND
C                          LESS THAN OR EQUAL TO D.  C MUST BE LESS
C                          THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (C,D) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE N+1 GRID POINTS IN THE Y-
C                          DIRECTION GIVEN BY YJ=C+(J-1)*DLY FOR
C                          J=1,2,...,N+1 WHERE DLY=(D-C)/N IS THE
C                          PANEL WIDTH.  IN ADDITION, N MUST BE
C                          GREATER THAN 4.
C
C                        NBDCND
C                          INDICATES THE TYPES OF BOUNDARY CONDITIONS
C                          AT Y=C AND Y=D
C                          = 0 IF THE SOLUTION IS PERIODIC IN Y,
C                              I.E., U(X,Y+D-C)=U(X,Y) FOR ALL X,Y
C                          = 1 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND Y = D, I.E., U(X,C)  AND U(X,D)
C                              ARE SPECIFIED FOR ALL X
C                          = 2 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND THE BOUNDARY CONDITION IS MIXED
C                              AT Y=D, I.E., DU(X,C)/DY AND U(X,D)
C                              ARE SPECIFIED FOR ALL X
C                          = 3 IF THE BOUNDARY CONDITIONS ARE MIXED
C                              AT Y=CAND Y=D I.E.,
C                              DU(X,D)/DY AND DU(X,D)/DY ARE
C                              SPECIFIED FOR ALL X
C                          = 4 IF THE BOUNDARY CONDITION IS MIXED
C                              AT Y=C AND THE SOLUTION IS SPECIFIED
C                              AT Y=D, I.E. DU(X,C)/DY+GAMA*U(X,C)
C                              AND U(X,D) ARE SPECIFIED FOR ALL X
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUE DU(X,C)/DY AT Y=C.
C
C                          WHEN NBDCND=3 OR 4
C                            BDC(I) = DU(XI,C)/DY I=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
C                          A DUMMY PARAMETER.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIED THE VALUE OF DU(X,D)/DY AT Y=D.
C
C                          WHEN NBDCND=2 OR 3
C                            BDD(I)=DU(XI,D)/DY I=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
C                          A DUMMY PARAMETER.
C
C                        COFX
C                          A USER-SUPPLIED SUBPROGRAM WITH PARAMETERS
C                          X, AFUN, BFUN, CFUN WHICH RETURNS THE
C                          VALUES OF THE X-DEPENDENT COEFFICIENTS
C                          AF(X), BF(X), CF(X) IN THE ELLIPTIC
C                          EQUATION AT X.  IF BOUNDARY CONDITIONS IN
C                          THE X DIRECTION ARE PERIODIC THEN THE
C                          COEFFICIENTS MUST SATISFY AF(X)=C1,BF(X)=0,
C                          CF(X)=C2 FOR ALL X.  HERE C1.GT.0
C                          AND C2 ARE CONSTANTS.
C
C                          NOTE THAT COFX MUST BE DECLARED EXTERNAL
C                          IN THE CALLING ROUTINE.
C
C                        GRHS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT-HAND SIDE OF THE
C                          ELLIPTIC EQUATION, I.E.,GRHS(I,J)=G(XI,YI),
C                          FOR I=2,...,M, J=2,...,N.  AT THE
C                          BOUNDARIES, GRHS IS DEFINED BY
C
C                          MBDCND   GRHS(1,J)   GRHS(M+1,J)
C                          ------   ---------   -----------
C                            0      G(A,YJ)     G(B,YJ)
C                            1         *           *
C                            2         *        G(B,YJ)  J=1,2,...,N+1
C                            3      G(A,YJ)     G(B,YJ)
C                            4      G(A,YJ)        *
C
C                          NBDCND   GRHS(I,1)   GRHS(I,N+1)
C                          ------   ---------   -----------
C                            0      G(XI,C)     G(XI,D)
C                            1         *           *
C                            2         *        G(XI,D)  I=1,2,...,M+1
C                            3      G(XI,C)     G(XI,D)
C                            4      G(XI,C)        *
C
C                          WHERE * MEANS THESE QUANTITES ARE NOT USED.
C                          GRHS SHOULD BE DIMENSIONED IDMN BY AT LEAST
C                          N+1 IN THE CALLING ROUTINE.
C
c *** caution              GRHS SHOULD BE RESET IF SEPX4 WAS FIRST CALLED
C                          WITH IORDER=2 AND WILL BE CALLED AGAIN WITH
C                          IORDER=4.  VALUES IN GRHS ARE DESTROYED BY THE
C                          IORDER=2 CALL.
C
C                        USOL
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE SOLUTION ALONG THE BOUNDARIES.
C                          AT THE BOUNDARIES, USOL IS DEFINED BY
C
C                          MBDCND   USOL(1,J)   USOL(M+1,J)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(A,YJ)     U(B,YJ)
C                            2      U(A,YJ)        *     J=1,2,...,N+1
C                            3         *           *
C                            4         *        U(B,YJ)
C
C                          NBDCND   USOL(I,1)   USOL(I,N+1)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(XI,C)     U(XI,D)
C                            2      U(XI,C)        *     I=1,2,...,M+1
C                            3         *           *
C                            4         *        U(XI,D)
C
C                          WHERE * MEANS THE QUANTITES ARE NOT USED
C                          IN THE SOLUTION.
C
C                          IF IORDER=2, THE USER MAY EQUIVALENCE GRHS
C                          AND USOL TO SAVE SPACE.  NOTE THAT IN THIS
C                          CASE THE TABLES SPECIFYING THE BOUNDARIES
C                          OF THE GRHS AND USOL ARRAYS DETERMINE THE
C                          BOUNDARIES UNIQUELY EXCEPT AT THE CORNERS.
C                          IF THE TABLES CALL FOR BOTH G(X,Y) AND
C                          U(X,Y) AT A CORNER THEN THE SOLUTION MUST
C                          BE CHOSEN.
C                          FOR EXAMPLE, IF MBDCND=2 AND NBDCND=4,
C                          THEN U(A,C), U(A,D),U(B,D) MUST BE CHOSEN
C                          AT THE CORNERS IN ADDITION TO G(B,C).
C
C                          IF IORDER=4, THEN THE TWO ARRAYS, USOL AND
C                          GRHS, MUST BE DISTINCT.
C
C                          USOL SHOULD BE DIMENSIONED IDMN BY AT LEAST
C                          N+1 IN THE CALLING ROUTINE.
C
C                        IDMN
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAYS
C                          GRHS AND USOL AS IT APPEARS IN THE PROGRAM
C                          CALLING SEPELI.  THIS PARAMETER IS USED
C                          TO SPECIFY THE VARIABLE DIMENSION OF GRHS
C                          AND USOL.  IDMN MUST BE AT LEAST 7 AND
C                          GREATER THAN OR EQUAL TO M+1.
C
C
C ON OUTPUT              USOL
C                          CONTAINS THE APPROXIMATE SOLUTION TO THE
C                          ELLIPTIC EQUATION. USOL(I,J) IS THE
C                          APPROXIMATION TO U(XI,YJ) FOR I=1,2...,M+1
C                          AND J=1,2,...,N+1.  THE APPROXIMATION HAS
C                          ERROR O(DLX**2+DLY**2) IF CALLED WITH
C                          IORDER=2 AND O(DLX**4+DLY**4) IF CALLED
C                          WITH IORDER=4.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS (I.E., ALPHA=BETA=0 IF
C                          MBDCND=3) IS SPECIFIED AND IF CF(X)=0 FOR
C                          ALL X THEN A SOLUTION TO THE DISCRETIZED
C                          MATRIX EQUATION MAY NOT EXIST
C                          (REFLECTING THE NON-UNIQUENESS OF SOLUTIONS
C                          TO THE PDE).
C                          PERTRB IS A CONSTANT CALCULATED AND
C                          SUBTRACTED FROM THE RIGHT HAND SIDE OF THE
C                          MATRIX EQUATION INSURING THE EXISTENCE OF A
C                          SOLUTION.  SEPX4 COMPUTES THIS SOLUTION
C                          WHICH IS A WEIGHTED MINIMAL LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL PROBLEM.  IF
C                          SINGULARITY IS NOT DETECTED PERTRB=0.0 IS
C                          RETURNED BY SEPX4.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS OR FAILURE TO FIND A SOLUTION
C
C                          =  0 NO ERROR
C                          =  1 IF A GREATER THAN B OR C GREATER
C                               THAN D
C                          =  2 IF MBDCND LESS THAN 0 OR MBDCND
C                               GREATER THAN 4
C                          =  3 IF NBDCND LESS THAN 0 OR NBDCND
C                               GREATER THAN 4
C                          =  4 IF ATTEMPT TO FIND A SOLUTION FAILS.
C                               (THE LINEAR SYSTEM GENERATED IS NOT
C                               DIAGONALLY DOMINANT.)
C                          =  5 IF IDMN IS TOO SMALL (SEE DISCUSSION
C                               OF IDMN)
C                          =  6 IF M IS TOO SMALL OR TOO LARGE
C                               (SEE DISCUSSION OF M)
C                          =  7 IF N IS TOO SMALL (SEE DISCUSSION OF N)
C                          =  8 IF IORDER IS NOT 2 OR 4
C                          =  9 IF INTL IS NOT 0 OR 1
C                          = 10 IF AFUN IS LESS THAN OR EQUAL TO ZERO
C                               FOR SOME INTERIOR MESH POINT XI SOME
C                               INTERIOR MESH POINT (XI,YJ)
C                          = 12 IF MBDCND=0 AND AF(X)=CF(X)=CONSTANT
C                               OR BF(X)=0 FOR ALL X IS NOT TRUE.
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C REQUIRED files         fish.f,comf.f,genbun.f,gnbnaux.f,sepaux.f
C
C
C PRECISION              SINGLE
C
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                SEPX4 WAS DEVELOPED AT NCAR BY JOHN C.
C                        ADAMS OF THE SCIENTIFIC COMPUTING DIVISION
C                        IN OCTOBER 1978.  THE BASIS OF THIS CODE IS
C                        NCAR ROUTINE SEPELI.  BOTH PACKAGES WERE
C                        RELEASED ON NCAR'S PUBLIC LIBRARIES IN
C                        JANUARY 1980. SEPX4 was modified in June 2004
c                        incorporating fortran 90 dynamical storage
c                        allocation for work space requirements
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              SEPX4 AUTOMATICALLY DISCRETIZES THE SEPARABLE
C                        ELLIPTIC EQUATION WHICH IS THEN SOLVED BY A
C                        GENERALIZED CYCLIC REDUCTION ALGORITHM IN THE
C                        SUBROUTINE POIS.  THE FOURTH ORDER SOLUTION
C                        IS OBTAINED USING THE TECHNIQUE OF DEFFERRED
C                        CORRECTIONS REFERENCED BELOW.
C
C TIMING                 WHEN POSSIBLE, SEPX4 SHOULD BE USED INSTEAD
C                        OF PACKAGE SEPELI.  THE INCREASE IN SPEED
C                        IS AT LEAST A FACTOR OF THREE.
C
C REFERENCES             KELLER, H.B., NUMERICAL METHODS FOR TWO-POINT
C                        BOUNDARY-VALUE PROBLEMS, BLAISDEL (1968),
C                        WALTHAM, MASS.
C
C                        SWARZTRAUBER, P., AND R. SWEET (1975):
C                        EFFICIENT FORTRAN SUBPROGRAMS FOR THE
C                        SOLUTION OF ELLIPTIC PARTIAL DIFFERENTIAL
C                        EQUATIONS.  NCAR TECHNICAL NOTE
C                          NCAR-TN/IA-109, PP. 135-137.
C***********************************************************************
      SUBROUTINE SEPX4(IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA, C
     1   , D, N, NBDCND, BDC, BDD, COFX, GRHS, USOL, IDMN, PERTRB, 
     2   IERROR)
      USE fish
      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: IORDER
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDMN
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: ALPHA
      REAL  :: BETA
      REAL  :: C
      REAL  :: D
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: GRHS(IDMN,*)
      REAL  :: USOL(IDMN,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: L, K, LOG2N, LENGTH, IRWK, ICWK, I1, I2, I3, I4, I5, I6
     1   , I7, I8, I9, I10, I11, I12, I13
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      external cofx
C-----------------------------------------------
C
C     CHECK INPUT PARAMETERS
C
      CALL C4KPRM(IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,IDMN,IERROR)
      IF (IERROR /= 0) RETURN 
C
C     COMPUTE MINIMUM WORK SPACE AND CHECK WORK SPACE LENGTH INPUT
C
      L = N + 1
      IF (NBDCND == 0) L = N
      K = M + 1
      L = N + 1
C     ESTIMATE LOG BASE 2 OF N
      LOG2N = INT(ALOG(FLOAT(N + 1))/ALOG(2.0) + 0.5)
!     set required work space estimate
      LENGTH = 4*(N + 1) + (10 + LOG2N)*(M + 1)
      IRWK = LENGTH + 6*(K + L) + 1
      ICWK = 0
!     allocate work space
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     return if allocation failed
      IF (IERROR == 20) RETURN 
      IERROR = 0
C
C     SET WORK SPACE INDICES
C
      I1 = LENGTH + 1
      I2 = I1 + L
      I3 = I2 + L
      I4 = I3 + L
      I5 = I4 + L
      I6 = I5 + L
      I7 = I6 + L
      I8 = I7 + K
      I9 = I8 + K
      I10 = I9 + K
      I11 = I10 + K
      I12 = I11 + K
      I13 = 1
      CALL S4ELIP(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     +NBDCND,BDC,BDD,COFX,w%rew(I1),w%rew(I2),w%rew(I3),
     +w%rew(i4),w%rew(i5),w%rew(i6),w%rew(i7),w%rew(i8),
     +w%rew(i9),w%rew(i10),w%rew(i11),w%rew(i12),
     +GRHS,USOL,IDMN,w%rew(i13),PERTRB,IERROR)
!     release dynamically allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE SEPX4


      SUBROUTINE S4ELIP(IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA, 
     1   C, D, N, NBDCND, BDC, BDD, COFX, AN, BN, CN, DN, UN, ZN, AM, BM
     2   , CM, DM, UM, ZM, GRHS, USOL, IDMN, W, PERTRB, IERROR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IORDER
      INTEGER , INTENT(IN) :: M
      INTEGER  :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDMN
      INTEGER , INTENT(INOUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL  :: ALPHA
      REAL  :: BETA
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL  :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL  :: DN(*)
      REAL  :: UN(*)
      REAL  :: ZN(*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL  :: DM(*)
      REAL  :: UM(*)
      REAL  :: ZM(*)
      REAL  :: GRHS(IDMN,*)
      REAL  :: USOL(IDMN,*)
      REAL  :: W(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J, I1, MP, NP, IEROR
      REAL :: XI, AI, BI, CI, AXI, BXI, CXI, DYJ, EYJ, FYJ, AX1, CXM, 
     1   DY1, FYN, GAMA, XNU, PRTRB
      LOGICAL :: SINGLR
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL COFX
C-----------------------------------------------
C
C     S4ELIP SETS UP VECTORS AND ARRAYS FOR INPUT TO BLKTRI
C     AND COMPUTES A SECOND ORDER SOLUTION IN USOL.  A RETURN JUMP TO
C     SEPELI OCCURRS IF IORDER=2.  IF IORDER=4 A FOURTH ORDER
C     SOLUTION IS GENERATED IN USOL.
C
C
C     SET PARAMETERS INTERNALLY
C
      KSWX = MBDCND + 1
      KSWY = NBDCND + 1
      K = M + 1
      L = N + 1
      AIT = A
      BIT = B
      CIT = C
      DIT = D
      DLY = (DIT - CIT)/FLOAT(N)
C
C     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
C     AND NON-SPECIFIED BOUNDARIES.
C
      USOL(2:M,2:N) = DLY**2*GRHS(2:M,2:N)
      IF (KSWX/=2 .AND. KSWX/=3) THEN
         USOL(1,2:N) = DLY**2*GRHS(1,2:N)
      ENDIF
      IF (KSWX/=2 .AND. KSWX/=5) THEN
         USOL(K,2:N) = DLY**2*GRHS(K,2:N)
      ENDIF
      IF (KSWY/=2 .AND. KSWY/=3) THEN
         USOL(2:M,1) = DLY**2*GRHS(2:M,1)
      ENDIF
      IF (KSWY/=2 .AND. KSWY/=5) THEN
         USOL(2:M,L) = DLY**2*GRHS(2:M,L)
      ENDIF
      IF (KSWX/=2 .AND. KSWX/=3 .AND. KSWY/=2 .AND. KSWY/=3) USOL(1,1)
     1    = DLY**2*GRHS(1,1)
      IF (KSWX/=2 .AND. KSWX/=5 .AND. KSWY/=2 .AND. KSWY/=3) USOL(K,1)
     1    = DLY**2*GRHS(K,1)
      IF (KSWX/=2 .AND. KSWX/=3 .AND. KSWY/=2 .AND. KSWY/=5) USOL(1,L)
     1    = DLY**2*GRHS(1,L)
      IF (KSWX/=2 .AND. KSWX/=5 .AND. KSWY/=2 .AND. KSWY/=5) USOL(K,L)
     1    = DLY**2*GRHS(K,L)
      I1 = 1
C
C     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
C
      MP = 1
      IF (KSWX == 1) MP = 0
      NP = NBDCND
C
C     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
C     IN NINT,MINT
C
      DLX = (BIT - AIT)/FLOAT(M)
      MIT = K - 1
      IF (KSWX == 2) MIT = K - 2
      IF (KSWX == 4) MIT = K
      DLY = (DIT - CIT)/FLOAT(N)
      NIT = L - 1
      IF (KSWY == 2) NIT = L - 2
      IF (KSWY == 4) NIT = L
      TDLX3 = 2.0*DLX**3
      DLX4 = DLX**4
      TDLY3 = 2.0*DLY**3
      DLY4 = DLY**4
C
C     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI
C
      IS = 1
      JS = 1
      IF (KSWX==2 .OR. KSWX==3) IS = 2
      IF (KSWY==2 .OR. KSWY==3) JS = 2
      NS = NIT + JS - 1
      MS = MIT + IS - 1
C
C     SET X - DIRECTION
C
      DO I = 1, MIT
         XI = AIT + FLOAT(IS + I - 2)*DLX
         CALL COFX (XI, AI, BI, CI)
         AXI = (AI/DLX - 0.5*BI)/DLX
         BXI = (-2.*AI/DLX**2) + CI
         CXI = (AI/DLX + 0.5*BI)/DLX
         AM(I) = DLY**2*AXI
         BM(I) = DLY**2*BXI
         CM(I) = DLY**2*CXI
      END DO
C
C     SET Y DIRECTION
C
      DYJ = 1.0
      EYJ = -2.0
      FYJ = 1.0
      AN(:NIT) = DYJ
      BN(:NIT) = EYJ
      CN(:NIT) = FYJ
C
C     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
C
      AX1 = AM(1)
      CXM = CM(MIT)
      SELECT CASE (KSWX) 
      CASE (2) 
C
C     DIRICHLET-DIRICHLET IN X DIRECTION
C
         AM(1) = 0.0
         CM(MIT) = 0.0
      CASE (3) 
C
C     DIRICHLET-MIXED IN X DIRECTION
C
         AM(1) = 0.0
         AM(MIT) = AM(MIT) + CXM
         BM(MIT) = BM(MIT) - 2.*BETA*DLX*CXM
         CM(MIT) = 0.0
      CASE (4) 
C
C     MIXED - MIXED IN X DIRECTION
C
         AM(1) = 0.0
         BM(1) = BM(1) + 2.*DLX*ALPHA*AX1
         CM(1) = CM(1) + AX1
         AM(MIT) = AM(MIT) + CXM
         BM(MIT) = BM(MIT) - 2.*DLX*BETA*CXM
         CM(MIT) = 0.0
      CASE (5) 
C
C     MIXED-DIRICHLET IN X DIRECTION
C
         AM(1) = 0.0
         BM(1) = BM(1) + 2.*ALPHA*DLX*AX1
         CM(1) = CM(1) + AX1
         CM(MIT) = 0.0
      END SELECT
c
C     ADJUST IN Y DIRECTION UNLESS PERIODIC
C
      DY1 = AN(1)
      FYN = CN(NIT)
      GAMA = 0.0
      XNU = 0.0
      SELECT CASE (KSWY) 
      CASE (2) 
C
C     DIRICHLET-DIRICHLET IN Y DIRECTION
C
         AN(1) = 0.0
         CN(NIT) = 0.0
      CASE (3) 
C
C     DIRICHLET-MIXED IN Y DIRECTION
C
         AN(1) = 0.0
         AN(NIT) = AN(NIT) + FYN
         BN(NIT) = BN(NIT) - 2.*DLY*XNU*FYN
         CN(NIT) = 0.0
      CASE (4) 
C
C     MIXED - MIXED DIRECTION IN Y DIRECTION
C
         AN(1) = 0.0
         BN(1) = BN(1) + 2.*DLY*GAMA*DY1
         CN(1) = CN(1) + DY1
         AN(NIT) = AN(NIT) + FYN
         BN(NIT) = BN(NIT) - 2.0*DLY*XNU*FYN
         CN(NIT) = 0.0
      CASE (5) 
C
C     MIXED-DIRICHLET IN Y DIRECTION
C
         AN(1) = 0.0
         BN(1) = BN(1) + 2.*DLY*GAMA*DY1
         CN(1) = CN(1) + DY1
         CN(NIT) = 0.0
      END SELECT
      IF (KSWX /= 1) THEN
C
C     ADJUST USOL ALONG X EDGE
C
         IF (KSWX==2 .OR. KSWX==3) THEN
            IF (KSWX==2 .OR. KSWX==5) THEN
               USOL(IS,JS:NS) = USOL(IS,JS:NS) - AX1*USOL(1,JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - CXM*USOL(K,JS:NS)
            ELSE
               USOL(IS,JS:NS) = USOL(IS,JS:NS) - AX1*USOL(1,JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - 2.0*DLX*CXM*BDB(JS:NS)
            ENDIF
         ELSE
            IF (KSWX==2 .OR. KSWX==5) THEN
               USOL(IS,JS:NS) = USOL(IS,JS:NS) + 2.0*DLX*AX1*BDA(JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - CXM*USOL(K,JS:NS)
            ELSE
               USOL(IS,JS:NS) = USOL(IS,JS:NS) + 2.0*DLX*AX1*BDA(JS:NS)
               USOL(MS,JS:NS) = USOL(MS,JS:NS) - 2.0*DLX*CXM*BDB(JS:NS)
            ENDIF
         ENDIF
      ENDIF
      IF (KSWY /= 1) THEN
C
C     ADJUST USOL ALONG Y EDGE
C
         IF (KSWY==2 .OR. KSWY==3) THEN
            IF (KSWY==2 .OR. KSWY==5) THEN
               USOL(IS:MS,JS) = USOL(IS:MS,JS) - DY1*USOL(IS:MS,1)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - FYN*USOL(IS:MS,L)
            ELSE
               USOL(IS:MS,JS) = USOL(IS:MS,JS) - DY1*USOL(IS:MS,1)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - 2.0*DLY*FYN*BDD(IS:MS)
            ENDIF
         ELSE
            IF (KSWY==2 .OR. KSWY==5) THEN
               USOL(IS:MS,JS) = USOL(IS:MS,JS) + 2.0*DLY*DY1*BDC(IS:MS)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - FYN*USOL(IS:MS,L)
            ELSE
               USOL(IS:MS,JS) = USOL(IS:MS,JS) + 2.0*DLY*DY1*BDC(IS:MS)
               USOL(IS:MS,NS) = USOL(IS:MS,NS) - 2.0*DLY*FYN*BDD(IS:MS)
            ENDIF
         ENDIF
      ENDIF
C
C     SAVE ADJUSTED EDGES IN GRHS IF IORDER=4
C
      IF (IORDER == 4) THEN
         GRHS(IS,JS:NS) = USOL(IS,JS:NS)
         GRHS(MS,JS:NS) = USOL(MS,JS:NS)
         GRHS(IS:MS,JS) = USOL(IS:MS,JS)
         GRHS(IS:MS,NS) = USOL(IS:MS,NS)
      ENDIF
      PERTRB = 0.0
C
C     CHECK IF OPERATOR IS SINGULAR
C
      CALL C4KSNG (MBDCND, NBDCND, ALPHA, BETA, COFX, SINGLR)
C
C     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
C     IF SINGULAR
C
      IF (SINGLR) CALL SEPTRI (MIT, AM, BM, CM, DM, UM, ZM)
      IF (SINGLR) CALL SEPTRI (NIT, AN, BN, CN, DN, UN, ZN)
C
C     ADJUST RIGHT HAND SIDE IF NECESSARY
C
      IF (SINGLR) CALL SEPORT (USOL, IDMN, ZN, ZM, PERTRB)
C
C     COMPUTE SOLUTION
C
C     SAVE ADJUSTED RIGHT HAND SIDE IN GRHS
      GRHS(IS:MS,JS:NS) = USOL(IS:MS,JS:NS)
      CALL GENBUN(NP,NIT,MP,MIT,AM,BM,CM,IDMN,USOL(IS,JS),IEROR)
C
C     CHECK IF ERROR DETECTED IN POIS
C     THIS CAN ONLY CORRESPOND TO IERROR=12
      IF (IEROR /= 0) THEN
C       SET ERROR FLAG IF IMPROPER COEFFICIENTS INPUT TO POIS
         IERROR = 12
         RETURN 
      ENDIF
      IF (IERROR /= 0) RETURN 
C
C     SET PERIODIC BOUNDARIES IF NECESSARY
C
      IF (KSWX == 1) THEN
         USOL(K,:L) = USOL(1,:L)
      ENDIF
      IF (KSWY == 1) THEN
         USOL(:K,L) = USOL(:K,1)
      ENDIF
C
C     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
C     NORM IF OPERATOR IS SINGULAR
C
      IF (SINGLR) CALL SEPMIN (USOL, IDMN, ZN, ZM, PRTRB)
C
C     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
C     NOT FLAGGED
C
      IF (IORDER == 2) RETURN 
C
C     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
C
      CALL D4FER (COFX, IDMN, USOL, GRHS)
      IF (SINGLR) CALL SEPORT (USOL, IDMN, ZN, ZM, PERTRB)
C
C     COMPUTE SOLUTION
C
C     SAVE ADJUSTED RIGHT HAND SIDE IN GRHS
      GRHS(IS:MS,JS:NS) = USOL(IS:MS,JS:NS)
      CALL GENBUN(NP,NIT,MP,MIT,AM,BM,CM,IDMN,USOL(IS,JS),IEROR)
C     CHECK IF ERROR DETECTED IN POIS
C     THIS CAN ONLY CORRESPOND TO IERROR=12
      IF (IEROR /= 0) THEN
C       SET ERROR FLAG IF IMPROPER COEFFICIENTS INPUT TO POIS
         IERROR = 12
         RETURN 
      ENDIF
      IF (IERROR /= 0) RETURN 
C
C     SET PERIODIC BOUNDARIES IF NECESSARY
C
      IF (KSWX == 1) THEN
         USOL(K,:L) = USOL(1,:L)
      ENDIF
      IF (KSWY == 1) THEN
         USOL(:K,L) = USOL(:K,1)
      ENDIF
C
C     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
C     NORM IF OPERATOR IS SINGULAR
C
      IF (SINGLR) CALL SEPMIN (USOL, IDMN, ZN, ZM, PRTRB)
      RETURN 
      END SUBROUTINE S4ELIP


      SUBROUTINE C4KPRM(IORDER, A, B, M, MBDCND, C, D, N, NBDCND, COFX, 
     1   IDMN, IERROR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IORDER
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER , INTENT(IN) :: IDMN
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I
      REAL :: DLX, XI, AI, BI, CI
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
C
C     THIS PROGRAM CHECKS THE INPUT PARAMETERS FOR ERRORS
C
C
C
C     CHECK DEFINITION OF SOLUTION REGION
C
      IF (A>=B .OR. C>=D) THEN
         IERROR = 1
         RETURN 
      ENDIF
C
C     CHECK BOUNDARY SWITCHES
C
      IF (MBDCND<0 .OR. MBDCND>4) THEN
         IERROR = 2
         RETURN 
      ENDIF
      IF (NBDCND<0 .OR. NBDCND>4) THEN
         IERROR = 3
         RETURN 
      ENDIF
C
C     CHECK FIRST DIMENSION IN CALLING ROUTINE
C
      IF (IDMN < 7) THEN
         IERROR = 5
         RETURN 
      ENDIF
C
C     CHECK M
C
      IF (M>IDMN - 1 .OR. M<6) THEN
         IERROR = 6
         RETURN 
      ENDIF
C
C     CHECK N
C
      IF (N < 5) THEN
         IERROR = 7
         RETURN 
      ENDIF
C
C     CHECK IORDER
C
      IF (IORDER/=2 .AND. IORDER/=4) THEN
         IERROR = 8
         RETURN 
      ENDIF
C
C     CHECK THAT EQUATION IS ELLIPTIC
C
      DLX = (B - A)/FLOAT(M)
      DO I = 2, M
         XI = A + FLOAT(I - 1)*DLX
         CALL COFX (XI, AI, BI, CI)
         IF (AI > 0.0) CYCLE 
         IERROR = 10
         RETURN 
      END DO
C
C     NO ERROR FOUND
C
      IERROR = 0
      RETURN 
      END SUBROUTINE C4KPRM


      SUBROUTINE C4KSNG(MBDCND, NBDCND, ALPHA, BETA, COFX, SINGLR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: NBDCND
      REAL , INTENT(IN) :: ALPHA
      REAL , INTENT(IN) :: BETA
      LOGICAL , INTENT(OUT) :: SINGLR
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I
      REAL :: XI, AI, BI, CI
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
C
C     THIS SUBROUTINE CHECKS IF THE PDE   SEPELI
C     MUST SOLVE IS A SINGULAR OPERATOR
C
      SINGLR = .FALSE.
C
C     CHECK IF THE BOUNDARY CONDITIONS ARE
C     ENTIRELY PERIODIC AND/OR MIXED
C
      IF(MBDCND/=0.AND.MBDCND/=3.OR.NBDCND/=0.AND.NBDCND/=3)RETURN
C
C     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
C
      IF (MBDCND == 3) THEN
         IF (ALPHA/=0.0 .OR. BETA/=0.0) RETURN 
      ENDIF
C
C     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
C     ARE ZERO
C
      DO I = IS, MS
         XI = AIT + FLOAT(I - 1)*DLX
         CALL COFX (XI, AI, BI, CI)
         IF (CI == 0.0) CYCLE 
         RETURN 
      END DO
C
C     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
C
      SINGLR = .TRUE.
      RETURN 
      END SUBROUTINE C4KSNG


      SUBROUTINE D4FER(COFX, IDMN, USOL, GRHS)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: IDMN
      REAL  :: USOL(IDMN,*)
      REAL , INTENT(INOUT) :: GRHS(IDMN,*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J
      REAL :: XI, AI, BI, CI, UXXX, UXXXX, UYYY, UYYYY, TX, TY
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
C
C     THIS SUBROUTINE FIRST APPROXIMATES THE TRUNCATION ERROR GIVEN BY
C     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY WHERE
C     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 ON THE INTERIOR AND
C     AT THE BOUNDARIES IF PERIODIC(HERE UXXX,UXXXX ARE THE THIRD
C     AND FOURTH PARTIAL DERIVATIVES OF U WITH RESPECT TO X).
C     TX IS OF THE FORM AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
C     AT X=A OR X=B IF THE BOUNDARY CONDITION THERE IS MIXED.
C     TX=0.0 ALONG SPECIFIED BOUNDARIES.  TY HAS SYMMETRIC FORM
C     IN Y WITH X,AFUN(X),BFUN(X) REPLACED BY Y,DFUN(Y),EFUN(Y).
C     THE SECOND ORDER SOLUTION IN USOL IS USED TO APPROXIMATE
C     (VIA SECOND ORDER FINITE DIFFERENCING) THE TRUNCATION ERROR
C     AND THE RESULT IS ADDED TO THE RIGHT HAND SIDE IN GRHS
C     AND THEN TRANSFERRED TO USOL TO BE USED AS A NEW RIGHT
C     HAND SIDE WHEN CALLING BLKTRI FOR A FOURTH ORDER SOLUTION.
C
C
C
C     COMPUTE TRUNCATION ERROR APPROXIMATION OVER THE ENTIRE MESH
C
      DO I = IS, MS
         XI = AIT + FLOAT(I - 1)*DLX
         CALL COFX (XI, AI, BI, CI)
         DO J = JS, NS
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
C
            CALL SEPDX (USOL, IDMN, I, J, UXXX, UXXXX)
            CALL SEPDY (USOL, IDMN, I, J, UYYY, UYYYY)
            TX = AI*UXXXX/12.0 + BI*UXXX/6.0
            TY = UYYYY/12.0
C
C     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
C
            IF (KSWX/=1 .AND. (I==1 .OR. I==K)) TX = AI/3.0*(UXXXX/4.0
     1          + UXXX/DLX)
            IF(KSWY/=1.AND.(J==1.OR.J==L))TY=(UYYYY/4.0+UYYY/DLY)/3.0
            GRHS(I,J) = GRHS(I,J) + DLY**2*(DLX**2*TX + DLY**2*TY)
         END DO
      END DO
C
C     RESET THE RIGHT HAND SIDE IN USOL
C
      USOL(IS:MS,JS:NS) = GRHS(IS:MS,JS:NS)
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 Changes
C-----------------------------------------------------------------------
      END SUBROUTINE D4FER
