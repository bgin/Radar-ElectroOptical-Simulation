C
C     file sepeli.f
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
C     SUBROUTINE SEPELI (INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,
C    +                   D,N,NBDCND,BDC,GAMA,BDD,XNU,COFX,COFY,GRHS,
C    +                   USOL,IDMN,W,PERTRB,IERROR)
C
C DIMENSION OF           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
C ARGUMENTS              USOL(IDMN,N+1),GRHS(IDMN,N+1),
C
C LATEST REVISION        JUNE 2004
C
C PURPOSE                SEPELI SOLVES FOR EITHER THE SECOND-ORDER
C                        FINITE DIFFERENCE APPROXIMATION OR A
C                        FOURTH-ORDER APPROXIMATION TO A SEPARABLE
C                        ELLIPTIC EQUATION
C
C                                 2    2
C                          AF(X)*D U/DX + BF(X)*DU/DX  + CF(X)*U +
C                                 2    2
C                          DF(Y)*D U/DY  + EF(Y)*DU/DY + FF(Y)*U
C
C                          = G(X,Y)
C
C                        ON A RECTANGLE (X GREATER THAN OR EQUAL TO A
C                        AND LESS THAN OR EQUAL TO B; Y GREATER THAN
C                        OR EQUAL TO C AND LESS THAN OR EQUAL TO D).
C                        ANY COMBINATION OF PERIODIC OR MIXED BOUNDARY
C                        CONDITIONS IS ALLOWED.
C
C                        THE POSSIBLE BOUNDARY CONDITIONS ARE:
C                        IN THE X-DIRECTION:
C                        (0) PERIODIC, U(X+B-A,Y)=U(X,Y) FOR ALL
C                            Y,X (1) U(A,Y), U(B,Y) ARE SPECIFIED FOR
C                            ALL Y
C                        (2) U(A,Y), DU(B,Y)/DX+BETA*U(B,Y) ARE
C                            SPECIFIED FOR ALL Y
C                        (3) DU(A,Y)/DX+ALPHA*U(A,Y),DU(B,Y)/DX+
C                            BETA*U(B,Y) ARE SPECIFIED FOR ALL Y
C                        (4) DU(A,Y)/DX+ALPHA*U(A,Y),U(B,Y) ARE
C                            SPECIFIED FOR ALL Y
C
C                        IN THE Y-DIRECTION:
C                        (0) PERIODIC, U(X,Y+D-C)=U(X,Y) FOR ALL X,Y
C                        (1) U(X,C),U(X,D) ARE SPECIFIED FOR ALL X
C                        (2) U(X,C),DU(X,D)/DY+XNU*U(X,D) ARE
C                            SPECIFIED FOR ALL X
C                        (3) DU(X,C)/DY+GAMA*U(X,C),DU(X,D)/DY+
C                            XNU*U(X,D) ARE SPECIFIED FOR ALL X
C                        (4) DU(X,C)/DY+GAMA*U(X,C),U(X,D) ARE
C                            SPECIFIED FOR ALL X
C
C USAGE                  CALL SEPELI (INTL,IORDER,A,B,M,MBDCND,BDA,
C                                     ALPHA,BDB,BETA,C,D,N,NBDCND,BDC,
C                                     GAMA,BDD,XNU,COFX,COFY,GRHS,USOL,
C                                     IDMN,W,PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               INTL
C                          = 0 ON INITIAL ENTRY TO SEPELI OR IF ANY
C                              OF THE ARGUMENTS C,D, N, NBDCND, COFY
C                              ARE CHANGED FROM A PREVIOUS CALL
C                          = 1 IF C, D, N, NBDCND, COFY ARE UNCHANGED
C                              FROM THE PREVIOUS CALL.
C
C                        IORDER
C                          = 2 IF A SECOND-ORDER APPROXIMATION
C                              IS SOUGHT
C                          = 4 IF A FOURTH-ORDER APPROXIMATION
C                              IS SOUGHT
C
C                        A,B
C                          THE RANGE OF THE X-INDEPENDENT VARIABLE,
C                          I.E., X IS GREATER THAN OR EQUAL TO A
C                          AND LESS THAN OR EQUAL TO B.  A MUST BE
C                          LESS THAN B.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL [A,B] IS SUBDIVIDED. HENCE,
C                          THERE WILL BE M+1 GRID POINTS IN THE X-
C                          DIRECTION GIVEN BY XI=A+(I-1)*DLX
C                          FOR I=1,2,...,M+1 WHERE DLX=(B-A)/M IS
C                          THE PANEL WIDTH.  M MUST BE LESS THAN
C                          IDMN AND GREATER THAN 5.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT X=A AND X=B
C
C                          = 0 IF THE SOLUTION IS PERIODIC IN X, I.E.,
C                              U(X+B-A,Y)=U(X,Y)  FOR ALL Y,X
C                          = 1 IF THE SOLUTION IS SPECIFIED AT X=A
C                              AND X=B, I.E., U(A,Y) AND U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C                          = 2 IF THE SOLUTION IS SPECIFIED AT X=A AND
C                              THE BOUNDARY CONDITION IS MIXED AT X=B,
C                              I.E., U(A,Y) AND DU(B,Y)/DX+BETA*U(B,Y)
C                              ARE SPECIFIED FOR ALL Y
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
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
C                          THAT SPECIFIES THE VALUES OF
C                          DU(A,Y)/DX+ ALPHA*U(A,Y) AT X=A, WHEN
C                          MBDCND=3 OR 4.
C                          BDA(J) = DU(A,YJ)/DX+ALPHA*U(A,YJ),
C                          J=1,2,...,N+1. WHEN MBDCND HAS ANY OTHER
C                          OTHER VALUE, BDA IS A DUMMY PARAMETER.
C
C                        ALPHA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT X=A
C                          (SEE ARGUMENT BDA).  IF MBDCND IS NOT
C                          EQUAL TO 3 OR 4 THEN ALPHA IS A DUMMY
C                          PARAMETER.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
C                          THAT SPECIFIES THE VALUES OF
C                          DU(B,Y)/DX+ BETA*U(B,Y) AT X=B.
C                          WHEN MBDCND=2 OR 3
C                          BDB(J) = DU(B,YJ)/DX+BETA*U(B,YJ),
C                          J=1,2,...,N+1. WHEN MBDCND HAS ANY OTHER
C                          OTHER VALUE, BDB IS A DUMMY PARAMETER.
C
C                        BETA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT
C                          X=B (SEE ARGUMENT BDB).  IF MBDCND IS
C                          NOT EQUAL TO 2 OR 3 THEN BETA IS A DUMMY
C                          PARAMETER.
C
C                        C,D
C                          THE RANGE OF THE Y-INDEPENDENT VARIABLE,
C                          I.E., Y IS GREATER THAN OR EQUAL TO C
C                          AND LESS THAN OR EQUAL TO D.  C MUST BE
C                          LESS THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL [C,D] IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS
C                          IN THE Y-DIRECTION GIVEN BY
C                          YJ=C+(J-1)*DLY FOR J=1,2,...,N+1 WHERE
C                          DLY=(D-C)/N IS THE PANEL WIDTH.
C                          IN ADDITION, N MUST BE GREATER THAN 4.
C
C                        NBDCND
C                          INDICATES THE TYPES OF BOUNDARY CONDITIONS
C                          AT Y=C AND Y=D
C
C                          = 0 IF THE SOLUTION IS PERIODIC IN Y,
C                              I.E., U(X,Y+D-C)=U(X,Y)  FOR ALL X,Y
C                          = 1 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND Y = D, I.E., U(X,C) AND U(X,D)
C                              ARE SPECIFIED FOR ALL X
C                          = 2 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND THE BOUNDARY CONDITION IS MIXED
C                              AT Y=D, I.E., U(X,C) AND
C                              DU(X,D)/DY+XNU*U(X,D) ARE SPECIFIED
C                              FOR ALL X
C                          = 3 IF THE BOUNDARY CONDITIONS ARE MIXED
C                              AT Y=C AND Y=D, I.E.,
C                              DU(X,D)/DY+GAMA*U(X,C) AND
C                              DU(X,D)/DY+XNU*U(X,D) ARE SPECIFIED
C                              FOR ALL X
C                          = 4 IF THE BOUNDARY CONDITION IS MIXED
C                              AT Y=C AND THE SOLUTION IS SPECIFIED
C                              AT Y=D, I.E. DU(X,C)/DY+GAMA*U(X,C)
C                              AND U(X,D) ARE SPECIFIED FOR ALL X
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
C                          THAT SPECIFIES THE VALUE OF
C                          DU(X,C)/DY+GAMA*U(X,C) AT Y=C.
C                          WHEN NBDCND=3 OR 4 BDC(I) = DU(XI,C)/DY +
C                          GAMA*U(XI,C), I=1,2,...,M+1.
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC
C                          IS A DUMMY PARAMETER.
C
C                        GAMA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT
C                          Y=C (SEE ARGUMENT BDC).  IF NBDCND IS
C                          NOT EQUAL TO 3 OR 4 THEN GAMA IS A DUMMY
C                          PARAMETER.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
C                          THAT SPECIFIES THE VALUE OF
C                          DU(X,D)/DY + XNU*U(X,D) AT Y=C.
C                          WHEN NBDCND=2 OR 3 BDD(I) = DU(XI,D)/DY +
C                          XNU*U(XI,D), I=1,2,...,M+1.
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD
C                          IS A DUMMY PARAMETER.
C
C                        XNU
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT
C                          Y=D (SEE ARGUMENT BDD).  IF NBDCND IS
C                          NOT EQUAL TO 2 OR 3 THEN XNU IS A
C                          DUMMY PARAMETER.
C
C                        COFX
C                          A USER-SUPPLIED SUBPROGRAM WITH
C                          PARAMETERS X, AFUN, BFUN, CFUN WHICH
C                          RETURNS THE VALUES OF THE X-DEPENDENT
C                          COEFFICIENTS AF(X), BF(X), CF(X) IN THE
C                          ELLIPTIC EQUATION AT X.
C
C                        COFY
C                          A USER-SUPPLIED SUBPROGRAM WITH PARAMETERS
C                          Y, DFUN, EFUN, FFUN WHICH RETURNS THE
C                          VALUES OF THE Y-DEPENDENT COEFFICIENTS
C                          DF(Y), EF(Y), FF(Y) IN THE ELLIPTIC
C                          EQUATION AT Y.
C
C                          NOTE:  COFX AND COFY MUST BE DECLARED
C                          EXTERNAL IN THE CALLING ROUTINE.
C                          THE VALUES RETURNED IN AFUN AND DFUN
C                          MUST SATISFY AFUN*DFUN GREATER THAN 0
C                          FOR A LESS THAN X LESS THAN B, C LESS
C                          THAN Y LESS THAN D (SEE IERROR=10).
C                          THE COEFFICIENTS PROVIDED MAY LEAD TO A
C                          MATRIX EQUATION WHICH IS NOT DIAGONALLY
C                          DOMINANT IN WHICH CASE SOLUTION MAY FAIL
C                          (SEE IERROR=4).
C
C                        GRHS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT-HAND SIDE OF THE
C                          ELLIPTIC EQUATION, I.E.,
C                          GRHS(I,J)=G(XI,YI), FOR I=2,...,M,
C                          J=2,...,N.  AT THE BOUNDARIES, GRHS IS
C                          DEFINED BY
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
C                          WHERE * MEANS THESE QUANTITIES ARE NOT USED.
C                          GRHS SHOULD BE DIMENSIONED IDMN BY AT LEAST
C                          N+1 IN THE CALLING ROUTINE.
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
C                          WHERE * MEANS THE QUANTITIES ARE NOT USED
C                          IN THE SOLUTION.
C
C                          IF IORDER=2, THE USER MAY EQUIVALENCE GRHS
C                          AND USOL TO SAVE SPACE.  NOTE THAT IN THIS
C                          CASE THE TABLES SPECIFYING THE BOUNDARIES
C                          OF THE GRHS AND USOL ARRAYS DETERMINE THE
C                          BOUNDARIES UNIQUELY EXCEPT AT THE CORNERS.
C                          IF THE TABLES CALL FOR BOTH G(X,Y) AND
C                          U(X,Y) AT A CORNER THEN THE SOLUTION MUST
C                          BE CHOSEN.  FOR EXAMPLE, IF MBDCND=2 AND
C                          NBDCND=4, THEN U(A,C), U(A,D), U(B,D) MUST
C                          BE CHOSEN AT THE CORNERS IN ADDITION
C                          TO G(B,C).
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
C                        W
c                          A fortran 90 derived TYPE (fishworkspace) variable
c                          that must be declared by the user.  The first
c                          two declarative statements in the user program
c                          calling SEPELI must be:
c
c                               USE fish
c                               TYPE (fishworkspace) :: W
c
c                          The first statement makes the fishpack module
c                          defined in the file "fish.f" available to the
c                          user program calling SEPELI.  The second statement
c                          declares a derived type variable (defined in
c                          the module "fish.f") which is used internally
c                          in SEPELI to dynamically allocate real and complex
c                          work space used in solution.  An error flag
c                          (IERROR = 20) is set if the required work space
c                          allocation fails (for example if N,M are too large)
c                          Real and complex values are set in the components
c                          of W on a initial (INTL=0) call to SEPELI.  These
c                          must be preserved on non-initial calls (INTL=1)
c                          to SEPELI.  This eliminates redundant calculations
c                          and saves compute time.
c               ****       IMPORTANT!  The user program calling SEPELI should
c                          include the statement:
c
c                               CALL FISHFIN(W)
C
C                          after the final approximation is generated by
C                          SEPELI.  The will deallocate the real and complex
c                          work space of W.  Failure to include this statement
c                          could result in serious memory leakage.
c
C ON OUTPUT              USOL
C                          CONTAINS THE APPROXIMATE SOLUTION TO THE
C                          ELLIPTIC EQUATION.
C                          USOL(I,J) IS THE APPROXIMATION TO U(XI,YJ)
C                          FOR I=1,2...,M+1 AND J=1,2,...,N+1.
C                          THE APPROXIMATION HAS ERROR
C                          O(DLX**2+DLY**2) IF CALLED WITH IORDER=2
C                          AND O(DLX**4+DLY**4) IF CALLED WITH
C                          IORDER=4.
C
C                        W
c                          The derived type (fishworkspace) variable W
c                          contains real and complex values that must not
C                          be destroyed if SEPELI is called again with
C                          INTL=1.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS
C                          (I.E., ALPHA=BETA=0 IF MBDCND=3;
C                          GAMA=XNU=0 IF NBDCND=3) IS SPECIFIED
C                          AND IF THE COEFFICIENTS OF U(X,Y) IN THE
C                          SEPARABLE ELLIPTIC EQUATION ARE ZERO
C                          (I.E., CF(X)=0 FOR X GREATER THAN OR EQUAL
C                          TO A AND LESS THAN OR EQUAL TO B;
C                          FF(Y)=0 FOR Y GREATER THAN OR EQUAL TO C
C                          AND LESS THAN OR EQUAL TO D) THEN A
C                          SOLUTION MAY NOT EXIST.  PERTRB IS A
C                          CONSTANT CALCULATED AND SUBTRACTED FROM
C                          THE RIGHT-HAND SIDE OF THE MATRIX EQUATIONS
C                          GENERATED BY SEPELI WHICH INSURES THAT A
C                          SOLUTION EXISTS. SEPELI THEN COMPUTES THIS
C                          SOLUTION WHICH IS A WEIGHTED MINIMAL LEAST
C                          SQUARES SOLUTION TO THE ORIGINAL PROBLEM.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS OR FAILURE TO FIND A SOLUTION
C                          = 0 NO ERROR
C                          = 1 IF A GREATER THAN B OR C GREATER THAN D
C                          = 2 IF MBDCND LESS THAN 0 OR MBDCND GREATER
C                              THAN 4
C                          = 3 IF NBDCND LESS THAN 0 OR NBDCND GREATER
C                              THAN 4
C                          = 4 IF ATTEMPT TO FIND A SOLUTION FAILS.
C                              (THE LINEAR SYSTEM GENERATED IS NOT
C                              DIAGONALLY DOMINANT.)
C                          = 5 IF IDMN IS TOO SMALL
C                              (SEE DISCUSSION OF IDMN)
C                          = 6 IF M IS TOO SMALL OR TOO LARGE
C                              (SEE DISCUSSION OF M)
C                          = 7 IF N IS TOO SMALL (SEE DISCUSSION OF N)
C                          = 8 IF IORDER IS NOT 2 OR 4
C                          = 9 IF INTL IS NOT 0 OR 1
C                          = 10 IF AFUN*DFUN LESS THAN OR EQUAL TO 0
C                               FOR SOME INTERIOR MESH POINT (XI,YJ)
C                          = 20 If the dynamic allocation of real and
C                               complex work space in the derived type
C                               (fishworkspace) variable W fails (e.g.,
c                               if N,M are too large for the platform used)
c
C                          NOTE (CONCERNING IERROR=4):  FOR THE
C                          COEFFICIENTS INPUT THROUGH COFX, COFY,
C                          THE DISCRETIZATION MAY LEAD TO A BLOCK
C                          TRIDIAGONAL LINEAR SYSTEM WHICH IS NOT
C                          DIAGONALLY DOMINANT (FOR EXAMPLE, THIS
C                          HAPPENS IF CFUN=0 AND BFUN/(2.*DLX) GREATER
C                          THAN AFUN/DLX**2).  IN THIS CASE SOLUTION
C                          MAY FAIL.  THIS CANNOT HAPPEN IN THE LIMIT
C                          AS DLX, DLY APPROACH ZERO.  HENCE, THE
C                          CONDITION MAY BE REMEDIED BY TAKING LARGER
C                          VALUES FOR M OR N.
C
C SPECIAL CONDITIONS     SEE COFX, COFY ARGUMENT DESCRIPTIONS ABOVE.
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED FILES         blktri.f,comf.f,sepaux.f,fish.f
C
C LANGUAGE               Fortran 90
C
C HISTORY                DEVELOPED AT NCAR DURING 1975-76 BY
C                        JOHN C. ADAMS OF THE SCIENTIFIC COMPUTING
C                        DIVISION.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980. Revised in June
C                        2004 using Fortan 90 dynamically allocated work
c                        space and derived data types to eliminate mixed
c                        mode conflicts in the earlier versions. All
c                        statement labels, arithmetic if statements and
c                        computed GO TO statements have been removed from
c                        the current version of SEPELI.
C
C ALGORITHM              SEPELI AUTOMATICALLY DISCRETIZES THE
C                        SEPARABLE ELLIPTIC EQUATION WHICH IS THEN
C                        SOLVED BY A GENERALIZED CYCLIC REDUCTION
C                        ALGORITHM IN THE SUBROUTINE, BLKTRI.  THE
C                        FOURTH-ORDER SOLUTION IS OBTAINED USING
C                        'DEFERRED CORRECTIONS' WHICH IS DESCRIBED
C                        AND REFERENCED IN SECTIONS, REFERENCES AND
C                        METHOD.
C
C TIMING                 THE OPERATIONAL COUNT IS PROPORTIONAL TO
C                        M*N*LOG2(N).
C
C ACCURACY               THE FOLLOWING ACCURACY RESULTS WERE OBTAINED
C                        using 64 bit floating point arithmetic.  Note
C                        THAT THE FOURTH-ORDER accuracy is not realized
C                        UNTIL THE MESH IS sufficiently refined.
C
C                                     SECOND-ORDER  FOURTH-ORDER
C                            M    N     ERROR         ERROR
C
C                             6    6    6.8E-1        1.2E0
C                            14   14    1.4E-1        1.8E-1
C                            30   30    3.2E-2        9.7E-3
C                            62   62    7.5E-3        3.0E-4
C                           126  126    1.8E-3        3.5E-6
C
C
C REFERENCES             KELLER, H.B., NUMERICAL METHODS FOR TWO-POINT
C                        BOUNDARY-VALUE PROBLEMS, BLAISDEL (1968),
C                        WALTHAM, MASS.
C
C                        SWARZTRAUBER, P., AND R. SWEET (1975):
C                        EFFICIENT FORTRAN SUBPROGRAMS FOR THE
C                        SOLUTION OF ELLIPTIC PARTIAL DIFFERENTIAL
C                        EQUATIONS.  NCAR TECHNICAL NOTE
C                        NCAR-TN/IA-109, PP. 135-137.
C***********************************************************************
 
      SUBROUTINE SEPELI(INTL, IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, 
     1   BETA, C, D, N, NBDCND, BDC, GAMA, BDD, XNU, COFX, COFY, GRHS, 
     2   USOL, IDMN, W, PERTRB, IERROR)
      USE fish
      implicit none
      TYPE (fishworkspace) :: w
      external cofx,cofy
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: INTL
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
      REAL  :: GAMA
      REAL  :: XNU
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
      INTEGER::I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,K,L,NP,IRWK,ICWK

      SAVE I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
!     save local variable work space pointers for noninitial call
!     check input arguments
      CALL CHKPRM (INTL, IORDER, A, B, M, MBDCND, C, D, N, NBDCND, COFX
     1   , COFY, IDMN, IERROR)
      IF (IERROR /= 0) RETURN 
      IF (INTL == 0) THEN
!     allocate space and set work space indices on initial call only
         K = M + 1
         L = N + 1
!          compute required blktri work space lengths
         NP = NBDCND
         CALL BLK_SPACE (N, M, IRWK, ICWK)
C
C     SET WORK SPACE INDICES
C
         I1 = IRWK + 1
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
!          set sepeli work space requirements
         IRWK = I12 + K
         ICWK = ICWK + 3*(M + 1)
!          allocate required real and complex work space
         CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!          return if allocation failure
         IF (IERROR == 20) RETURN 
      ENDIF
      IERROR = 0
!     compute second or fourth order solution
      CALL SPELIP (INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     +NBDCND,BDC,GAMA,BDD,XNU,COFX,COFY,w%rew(i1),w%rew(i2),w%rew(i3),
     +w%rew(i4),w%rew(i5),w%rew(i6),w%rew(i7),w%rew(i8),w%rew(i9),
     +w%rew(i10),w%rew(i11),w%rew(i12),grhs,usol,idmn,w%rew,w%cxw,
     +pertrb,ierror)
      RETURN 
      END SUBROUTINE SEPELI


      SUBROUTINE SPELIP(INTL, IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, 
     1   BETA, C, D, N, NBDCND, BDC, GAMA, BDD, XNU, COFX, COFY, AN, BN
     2   , CN, DN, UN, ZN, AM, BM, CM, DM, UM, ZM, GRHS, USOL, IDMN, W, 
     3   WC, PERTRB, IERROR)
      implicit none
      external cofx,cofy
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: INTL
      INTEGER , INTENT(IN) :: IORDER
      INTEGER , INTENT(IN) :: M
      INTEGER  :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDMN
      INTEGER  :: IERROR
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL  :: ALPHA
      REAL  :: BETA
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL  :: GAMA
      REAL  :: XNU
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
      COMPLEX  :: WC(*)
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
      INTEGER :: I, J, I1, MP, NP
      REAL :: XI, AI, BI, CI, AXI, BXI, CXI, YJ, DJ, EJ, FJ, DYJ, EYJ, 
     1   FYJ, AX1, CXM, DY1, FYN, PRTRB
      LOGICAL :: SINGLR
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
C
C     SPELIP SETS UP VECTORS AND ARRAYS FOR INPUT TO BLKTRI
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
C
C     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
C     AND NON-SPECIFIED BOUNDARIES.
C
      USOL(2:M,2:N) = GRHS(2:M,2:N)
      IF (KSWX/=2 .AND. KSWX/=3) THEN
         USOL(1,2:N) = GRHS(1,2:N)
      ENDIF
      IF (KSWX/=2 .AND. KSWX/=5) THEN
         USOL(K,2:N) = GRHS(K,2:N)
      ENDIF
      IF (KSWY/=2 .AND. KSWY/=3) THEN
         USOL(2:M,1) = GRHS(2:M,1)
      ENDIF
      IF (KSWY/=2 .AND. KSWY/=5) THEN
         USOL(2:M,L) = GRHS(2:M,L)
      ENDIF
      IF (KSWX/=2 .AND. KSWX/=3 .AND. KSWY/=2 .AND. KSWY/=3) USOL(1,1)
     1    = GRHS(1,1)
      IF (KSWX/=2 .AND. KSWX/=5 .AND. KSWY/=2 .AND. KSWY/=3) USOL(K,1)
     1    = GRHS(K,1)
      IF (KSWX/=2 .AND. KSWX/=3 .AND. KSWY/=2 .AND. KSWY/=5) USOL(1,L)
     1    = GRHS(1,L)
      IF (KSWX/=2 .AND. KSWX/=5 .AND. KSWY/=2 .AND. KSWY/=5) USOL(K,L)
     1    = GRHS(K,L)
      I1 = 1
C
C     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
C
      MP = 1
      NP = 1
      IF (KSWX == 1) MP = 0
      IF (KSWY == 1) NP = 0
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
         AM(I) = AXI
         BM(I) = BXI
         CM(I) = CXI
      END DO
C
C     SET Y DIRECTION
C
      DO J = 1, NIT
         YJ = CIT + FLOAT(JS + J - 2)*DLY
         CALL COFY (YJ, DJ, EJ, FJ)
         DYJ = (DJ/DLY - 0.5*EJ)/DLY
         EYJ = (-2.*DJ/DLY**2) + FJ
         FYJ = (DJ/DLY + 0.5*EJ)/DLY
         AN(J) = DYJ
         BN(J) = EYJ
         CN(J) = FYJ
      END DO
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
      CASE (5) 
C
C     MIXED-DIRICHLET IN X DIRECTION
C
         AM(1) = 0.0
         BM(1) = BM(1) + 2.*ALPHA*DLX*AX1
         CM(1) = CM(1) + AX1
         CM(MIT) = 0.0
      CASE (3) 
C
C     DIRICHLET-MIXED IN X DIRECTION
C
         AM(1) = 0.0
         AM(MIT) = AM(MIT) + CXM
         BM(MIT) = BM(MIT) - 2.*BETA*DLX*CXM
         CM(MIT) = 0.0
C
C     MIXED - MIXED IN X DIRECTION
C
      CASE (4) 
         AM(1) = 0.0
         BM(1) = BM(1) + 2.*DLX*ALPHA*AX1
         CM(1) = CM(1) + AX1
         AM(MIT) = AM(MIT) + CXM
         BM(MIT) = BM(MIT) - 2.*DLX*BETA*CXM
         CM(MIT) = 0.0
      END SELECT
C
C     ADJUST IN Y DIRECTION UNLESS PERIODIC
C
      DY1 = AN(1)
      FYN = CN(NIT)
      SELECT CASE (KSWY) 
      CASE (2) 
C
C     DIRICHLET-DIRICHLET IN Y DIRECTION
C
         AN(1) = 0.0
         CN(NIT) = 0.0
      CASE (5) 
C
C     MIXED-DIRICHLET IN Y DIRECTION
C
         AN(1) = 0.0
         BN(1) = BN(1) + 2.*DLY*GAMA*DY1
         CN(1) = CN(1) + DY1
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
      CALL CHKSNG(MBDCND,NBDCND,ALPHA,BETA,GAMA,XNU,COFX,COFY,SINGLR)
C
C     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
C     IF SINGULAR
C
      IF (SINGLR) CALL SEPTRI (MIT, AM, BM, CM, DM, UM, ZM)
      IF (SINGLR) CALL SEPTRI (NIT, AN, BN, CN, DN, UN, ZN)
C
C     MAKE INITIALIZATION CALL TO blktrii
C
      IF (INTL == 0) THEN
         CALL BLKTRII (INTL, NP, NIT, AN, BN, CN, MP, MIT, AM, BM, CM, 
     1      IDMN, USOL(IS,JS), IERROR, W, WC)
         IF (IERROR /= 0) RETURN 
      ENDIF
C
C     ADJUST RIGHT HAND SIDE IF NECESSARY
C
      IF (SINGLR) CALL SEPORT (USOL, IDMN, ZN, ZM, PERTRB)
C
C     COMPUTE SOLUTION
C
      CALL BLKTRII (I1, NP, NIT, AN, BN, CN, MP, MIT, AM, BM, CM, IDMN, 
     1   USOL(IS,JS), IERROR, W, WC)
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
      CALL DEFER (COFX, COFY, IDMN, USOL, GRHS)
      IF (SINGLR) CALL SEPORT (USOL, IDMN, ZN, ZM, PERTRB)
C
C     COMPUTE fourth order SOLUTION
C
      CALL BLKTRII (I1, NP, NIT, AN, BN, CN, MP, MIT, AM, BM, CM, IDMN, 
     1   USOL(IS,JS), IERROR, W, WC)
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
      END SUBROUTINE SPELIP


 
      SUBROUTINE CHKPRM(INTL, IORDER, A, B, M, MBDCND, C, D, N, NBDCND, 
     1   COFX, COFY, IDMN, IERROR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: INTL
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
      INTEGER :: I, J
      REAL :: DLX, DLY, XI, AI, BI, CI, YJ, DJ, EJ, FJ
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
C
C     THIS PROGRAM CHECKS THE INPUT arguments FOR ERRORS
C
C
C     CHECK DEFINITION OF SOLUTION REGION
C
      IF (A>=B .OR. C>=D) THEN
         IERROR = 1
         RETURN 
      ENDIF
c
c     check boundary condition arguments
c
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
C     CHECK M,N
C
      IF (M>IDMN - 1 .OR. M<6) THEN
         IERROR = 6
         RETURN 
      ENDIF
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
C     CHECK INTL
C
      IF (INTL/=0 .AND. INTL/=1) THEN
         IERROR = 9
         RETURN 
      ENDIF
C
C     CHECK THAT EQUATION IS ELLIPTIC (only on initial call)
C
      IF (INTL == 0) THEN
         DLX = (B - A)/FLOAT(M)
         DLY = (D - C)/FLOAT(N)
         DO I = 2, M
            XI = A + FLOAT(I - 1)*DLX
            CALL COFX (XI, AI, BI, CI)
            DO J = 2, N
               YJ = C + FLOAT(J - 1)*DLY
               CALL COFY (YJ, DJ, EJ, FJ)
               IF (AI*DJ > 0.0) CYCLE 
               IERROR = 10
               RETURN 
            END DO
         END DO
      ENDIF
C
C     NO ERROR FOUND
C
      IERROR = 0
      RETURN 
      END SUBROUTINE CHKPRM


      SUBROUTINE CHKSNG(MBDCND, NBDCND, ALPHA, BETA, GAMA, XNU, COFX, 
     1   COFY, SINGLR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: NBDCND
      REAL , INTENT(IN) :: ALPHA
      REAL , INTENT(IN) :: BETA
      REAL , INTENT(IN) :: GAMA
      REAL , INTENT(IN) :: XNU
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
      INTEGER :: I, J
      REAL :: XI, AI, BI, CI, YJ, DJ, EJ, FJ
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
 
      IF (NBDCND == 3) THEN
         IF (GAMA/=0.0 .OR. XNU/=0.0) RETURN 
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
      DO J = JS, NS
         YJ = CIT + FLOAT(J - 1)*DLY
         CALL COFY (YJ, DJ, EJ, FJ)
         IF (FJ == 0.0) CYCLE 
         RETURN 
      END DO
C
C     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
C
      SINGLR = .TRUE.
      RETURN 
      END SUBROUTINE CHKSNG


      SUBROUTINE DEFER(COFX, COFY, IDMN, USOL, GRHS)
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
      INTEGER :: J, I
      REAL::YJ,DJ,EJ,FJ,XI,AI,BI,CI,UXXX,UXXXX,UYYY,UYYYY,TX,TY
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
C     COMPUTE TRUNCATION ERROR APPROXIMATION OVER THE ENTIRE MESH
C
      DO J = JS, NS
         YJ = CIT + FLOAT(J - 1)*DLY
         CALL COFY (YJ, DJ, EJ, FJ)
         DO I = IS, MS
            XI = AIT + FLOAT(I - 1)*DLX
            CALL COFX (XI, AI, BI, CI)
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
C
            CALL SEPDX (USOL, IDMN, I, J, UXXX, UXXXX)
            CALL SEPDY (USOL, IDMN, I, J, UYYY, UYYYY)
            TX = AI*UXXXX/12.0 + BI*UXXX/6.0
            TY = DJ*UYYYY/12.0 + EJ*UYYY/6.0
C
C     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
C
            IF (KSWX/=1 .AND. (I==1 .OR. I==K)) TX = AI/3.0*(UXXXX/4.0
     1          + UXXX/DLX)
            IF (KSWY/=1 .AND. (J==1 .OR. J==L)) TY = DJ/3.0*(UYYYY/4.0
     1          + UYYY/DLY)
            GRHS(I,J) = GRHS(I,J) + DLX**2*TX + DLY**2*TY
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
c June      2004    version 5.0, fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE DEFER
