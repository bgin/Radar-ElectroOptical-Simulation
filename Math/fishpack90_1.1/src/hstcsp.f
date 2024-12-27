C
C     file hstcsp.f
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
C     SUBROUTINE HSTCSP (INTL,A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,
C    +                   BDD,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION ON A STAGGERED
C                        GRID TO THE MODIFIED HELMHOLTZ EQUATION IN
C                        SPHERICAL COORDINATES ASSUMING AXISYMMETRY
C                        (NO DEPENDENCE ON LONGITUDE).
C
C                        THE EQUATION IS
C
C                           (1/R**2)(D/DR)(R**2(DU/DR)) +
C                           1/(R**2*SIN(THETA))(D/DTHETA)
C                           (SIN(THETA)(DU/DTHETA)) +
C                           (LAMBDA/(R*SIN(THETA))**2)U  =  F(THETA,R)
C
C                        WHERE THETA IS COLATITUDE AND R IS THE
C                        RADIAL COORDINATE. THIS TWO-DIMENSIONAL
C                        MODIFIED HELMHOLTZ EQUATION RESULTS FROM
C                        THE FOURIER TRANSFORM OF THE THREE-
C                        DIMENSIONAL POISSON EQUATION.
C
C
C USAGE                  CALL HSTCSP (INTL,A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR,W)
C
C ARGUMENTS
C  ON INPUT              INTL
C
C                          = 0  ON INITIAL ENTRY TO HSTCSP OR IF ANY
C                               OF THE ARGUMENTS C, D, N, OR NBDCND
C                               ARE CHANGED FROM A PREVIOUS CALL
C
C                          = 1  IF C, D, N, AND NBDCND ARE ALL
C                               UNCHANGED FROM PREVIOUS CALL TO HSTCSP
C
C                          NOTE:
C                          A CALL WITH INTL = 0 TAKES APPROXIMATELY
C                          1.5 TIMES AS MUCH TIME AS A CALL WITH
C                          INTL = 1.  ONCE A CALL WITH INTL = 0
C                          HAS BEEN MADE THEN SUBSEQUENT SOLUTIONS
C                          CORRESPONDING TO DIFFERENT F, BDA, BDB,
C                          BDC, AND BDD CAN BE OBTAINED FASTER WITH
C                          INTL = 1 SINCE INITIALIZATION IS NOT
C                          REPEATED.
C
C                        A,B
C                          THE RANGE OF THETA (COLATITUDE),
C                          I.E. A .LE. THETA .LE. B.  A
C                          MUST BE LESS THAN B AND A MUST BE
C                          NON-NEGATIVE.  A AND B ARE IN RADIANS.
C                          A = 0 CORRESPONDS TO THE NORTH POLE AND
C                          B = PI CORRESPONDS TO THE SOUTH POLE.
C
C                          * * *  IMPORTANT  * * *
C
C                          IF B IS EQUAL TO PI, THEN B MUST BE
C                          COMPUTED USING THE STATEMENT
C                              B = PIMACH(DUM)
C                          THIS INSURES THAT B IN THE USER'S PROGRAM
C                          IS EQUAL TO PI IN THIS PROGRAM, PERMITTING
C                          SEVERAL TESTS OF THE INPUT PARAMETERS THAT
C                          OTHERWISE WOULD NOT BE POSSIBLE.
C
C                          * * * * * * * * * * * *
C
C                        M
C                          THE NUMBER OF GRID POINTS IN THE INTERVAL
C                          (A,B).  THE GRID POINTS IN THE THETA-
C                          DIRECTION ARE GIVEN BY
C                            THETA(I) = A + (I-0.5)DTHETA
C                          FOR I=1,2,...,M WHERE DTHETA =(B-A)/M.
C                          M MUST BE GREATER THAN 4.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT THETA = A AND THETA = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THETA = B.
C                               (SEE NOTES 1, 2 BELOW)
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = B
C                               (SEE NOTES 1, 2 BELOW).
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = A (SEE NOTES 1, 2 BELOW)
C                               AND THETA = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED AT
C                               THETA = A (SEE NOTES 1, 2 BELOW) AND
C                               THE SOLUTION IS SPECIFIED AT THETA = B.
C
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = A = 0 AND THE SOLUTION IS
C                               SPECIFIED AT THETA = B.
C                               (SEE NOTE 2 BELOW)
C
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = A = 0 AND THE DERIVATIVE OF
C                               THE SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = B
C                               (SEE NOTE 2 BELOW).
C
C                          = 7  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THE SOLUTION IS
C                               UNSPECIFIED AT THETA = B = PI.
C
C                          = 8  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED AT
C                               THETA = A (SEE NOTE 1 BELOW)
C                               AND THE SOLUTION IS UNSPECIFIED AT
C                               THETA = B = PI.
C
C                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
C                                THETA = A = 0 AND THETA = B = PI.
C
C                          NOTE 1:
C                          IF A = 0, DO NOT USE MBDCND = 1,2,3,4,7
C                          OR 8, BUT INSTEAD USE MBDCND = 5, 6, OR 9.
C
C                          NOTE 2:
C                          IF B = PI, DO NOT USE MBDCND = 1,2,3,4,5,
C                          OR 6, BUT INSTEAD USE MBDCND = 7, 8, OR 9.
C
C                          NOTE 3:
C                          WHEN A = 0  AND/OR B = PI THE ONLY
C                          MEANINGFUL BOUNDARY CONDITION IS
C                          DU/DTHETA = 0.   SEE D. GREENSPAN,
C                          'NUMERICAL ANALYSIS OF ELLIPTIC
C                           BOUNDARY VALUE PROBLEMS,'
C                          HARPER AND ROW, 1965, CHAPTER 5.)
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES (IF ANY) OF
C                          THE SOLUTION AT THETA = A.
C
C                          WHEN  MBDCND = 1, 2, OR 7,
C                            BDA(J) = U(A,R(J)),   J=1,2,...,N.
C
C                          WHEN MBDCND = 3, 4, OR 8,
C                            BDA(J) = (D/DTHETA)U(A,R(J)), J=1,2,...,N.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS A
C                          DUMMY VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT THETA = B.
C
C                          WHEN MBDCND = 1, 4, OR 5,
C                            BDB(J) = U(B,R(J)),     J=1,2,...,N.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DTHETA)U(B,R(J)), J=1,2,...,N.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
C                          A DUMMY VARIABLE.
C
C                        C,D
C                          THE RANGE OF R , I.E. C .LE. R .LE. D.
C                          C MUST BE LESS THAN D AND NON-NEGATIVE.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE R-DIRECTION
C                          ARE GIVEN BY R(J) = C + (J-0.5)DR,
C                          J=1,2,...,N, WHERE DR = (D-C)/N.
C                          N MUST BE GREATER THAN 4.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT R = C AND R = D.
C
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               R = C AND R = D.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               R = C AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO R IS
C                               SPECIFIED AT R = D. (SEE NOTE 1 BELOW)
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = C AND R = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS
C                               SPECIFIED AT R = C AND THE SOLUTION
C                               IS SPECIFIED AT R = D.
C
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = C = 0 (SEE NOTE 2 BELOW) AND THE
C                               SOLUTION IS SPECIFIED AT R = D.
C
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = C = 0 (SEE NOTE 2 BELOW)
C                               AND THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = D.
C
C                          NOTE 1:
C                          IF C = 0 AND MBDCND = 3,6,8 OR 9, THE
C                          SYSTEM OF EQUATIONS TO BE SOLVED IS
C                          SINGULAR.  THE UNIQUE SOLUTION IS
C                          DETERMINED BY EXTRAPOLATION TO THE
C                          SPECIFICATION OF U(THETA(1),C).
C                          BUT IN THESE CASES THE RIGHT SIDE OF THE
C                          SYSTEM WILL BE PERTURBED BY THE CONSTANT
C                          PERTRB.
C
C                          NOTE 2:
C                          NBDCND = 5 OR 6 CANNOT BE USED WITH
C                          MBDCND =1, 2, 4, 5, OR 7
C                          (THE FORMER INDICATES THAT THE SOLUTION IS
C                          UNSPECIFIED AT R = 0; THE LATTER INDICATES
C                          SOLUTION IS SPECIFIED).
C                          USE INSTEAD NBDCND = 1 OR 2.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT R = C.  WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(THETA(I),C),    I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DR)U(THETA(I),C), I=1,2,...,M.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
C                          A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT R = D.  WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(THETA(I),D) ,    I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DR)U(THETA(I),D), I=1,2,...,M.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE MODIFIED
C                          HELMHOLTZ EQUATION.  IF LAMBDA IS GREATER
C                          THAN 0, A SOLUTION MAY NOT EXIST.
C                          HOWEVER, HSTCSP WILL ATTEMPT TO FIND A
C                          SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT SIDE OF THE MODIFIED
C                          HELMHOLTZ EQUATION.  FOR I=1,2,...,M AND
C                          J=1,2,...,N
C
C                                F(I,J) = F(THETA(I),R(J)) .
C
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTCSP.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.
C                          IDIMF MUST BE AT LEAST M.
C
C                        W
c                          A fortran 90 derived TYPE (fishworkspace) variable
c                          that must be declared by the user.  The first
c                          two declarative statements in the user program
c                          calling HSTCSP must be:
c
c                               USE fish
c                               TYPE (fishworkspace) :: W
c
c                          The first statement makes the fishpack module
c                          defined in the file "fish.f" available to the
c                          user program calling HSTCSP.  The second statement
c                          declares a derived type variable (defined in
c                          the module "fish.f") which is used internally
c                          in BLKTRI to dynamically allocate real and complex
c                          work space used in solution.  An error flag
c                          (IERROR = 20) is set if the required work space
c                          allocation fails (for example if N,M are too large)
c                          Real and complex values are set in the components
c                          of W on a initial (IFLG=0) call to HSTCSP.  These
c                          must be preserved on non-initial calls (INTL=1)
c                          to HSTCSP.  This eliminates redundant calculations
c                          and saves compute time.
c               ****       IMPORTANT!  The user program calling HSTCSP should
c                          include the statement:
c
c                               CALL FISHFIN(W)
C
C                          after the final approximation is generated by
C                          HSTCSP.  The will deallocate the real and complex
c                          work space of W.  Failure to include this statement
c                          could result in serious memory leakage.
C
C                                                                       
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (THETA(I),R(J)) FOR I=1,2,..,M, J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HSTCSP THEN COMPUTES THIS
C                          SOLUTION, WHICH IS A LEAST SQUARES SOLUTION
C                          TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
C                          A SOLUTION; HENCE, THE SOLUTION IS NOT
C                          UNIQUE.  THE VALUE OF PERTRB SHOULD BE
C                          SMALL COMPARED TO THE RIGHT SIDE F.
C                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
C                          ESSENTIALLY DIFFERENT PROBLEM.
C                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
C                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS. EXCEPT FOR NUMBERS 0 AND 10,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C
C                          =  1  A .LT. 0 OR B .GT. PI
C
C                          =  2  A .GE. B
C
C                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 9
C
C                          =  4  C .LT. 0
C
C                          =  5  C .GE. D
C
C                          =  6  NBDCND .LT. 1 OR NBDCND .GT. 6
C
C                          =  7  N .LT. 5
C
C                          =  8  NBDCND = 5 OR 6 AND
C                                MBDCND = 1, 2, 4, 5, OR 7
C
C                          =  9  C .GT. 0 AND NBDCND .GE. 5
C
C                          = 10  ELMBDA .GT. 0
C
C                          = 11  IDIMF .LT. M
C
C                          = 12  M .LT. 5
C
C                          = 13  A = 0 AND MBDCND =1,2,3,4,7 OR 8
C
C                          = 14  B = PI AND MBDCND .LE. 6
C
C                          = 15  A .GT. 0 AND MBDCND = 5, 6, OR 9
C
C                          = 16  B .LT. PI AND MBDCND .GE. 7
C
C                          = 17  LAMBDA .NE. 0 AND NBDCND .GE. 5
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HSTCSP,
C                          THE USER SHOULD TEST IERROR AFTER THE CALL.
C
C                        = 20 If the dynamic allocation of real and
C                             complex work space in the derived type
C                             (fishworkspace) variable W fails (e.g.,
c                             if N,M are too large for the platform used)
C                                                                       
C                        W
c                             The derived type (fishworkspace) variable W
c                             contains real and complex values that must not
C                             be destroyed if HSTCSP is called again with
C                             IFLG=1.
C                                                                       
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       fish.f,blktri.f,comf.f
C FILES
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
C                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
C                        IN JANUARY 1980. Revised by John Adams in June
C                        2004 using Fortan 90 dynamically allocated work
c                        space and derived data types to eliminate mixed
c                        mode conflicts in the earlier versions.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, ADJUSTS
C                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR
C                        AND CALLS BLKTRI WHICH SOLVES THE LINEAR
C                        SYSTEM OF EQUATIONS.
C
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT IS
C                        ROUGHLY PROPORTIONAL TO M*N*LOG2(N).  THE
C                        TIMING ALSO DEPENDS ON INPUT PARAMETER INTL.
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
C                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
C                        DIGITS FOR N AND M AS LARGE AS 64.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        SUBROUTINE BLKTRI WHICH IS THE ROUTINE
C                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
C
C REFERENCES             P.N. SWARZTRAUBER, "A DIRECT METHOD FOR
C                        THE DISCRETE SOLUTION OF SEPARABLE ELLIPTIC
C                        EQUATIONS",
C                        SIAM J. NUMER. ANAL. 11(1974), PP. 1136-1150.
C
C                        U. SCHUMANN AND R. SWEET, "A DIRECT METHOD FOR
C                        THE SOLUTION OF POISSON'S EQUATION WITH NEUMANN
C                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
C                        ARBITRARY SIZE," J. COMP. PHYS. 20(1976),
C                        PP. 171-182.
C***********************************************************************
      SUBROUTINE HSTCSP(INTL, A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND
     1   , BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      USE fish
      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: INTL
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: A
      REAL  :: B
      REAL  :: C
      REAL  :: D
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDA(*)
      REAL  :: BDB(*)
      REAL  :: BDC(*)
      REAL  :: BDD(*)
      REAL  :: F(IDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IW1, IWBM, IWCM, IWAN, IWBN, IWCN, IWSNTH, IWRSQ, K, L
     1   , NP, IRWK, ICWK, IERR1
      REAL :: PI, DUM

      SAVE IW1, IWBM, IWCM, IWAN, IWBN, IWCN, IWSNTH, IWRSQ
C-----------------------------------------------
c     USE fish
c     TYPE (fishworkspace) :: w
      PI = 4.0*ATAN(1.0)
C
C     CHECK FOR INVALID INPUT PARAMETERS
C
      IERROR = 0
      IF (A<0. .OR. B>PI) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<1 .OR. MBDCND>9) IERROR = 3
      IF (C < 0.) IERROR = 4
      IF (C >= D) IERROR = 5
      IF (NBDCND<1 .OR. NBDCND>6) IERROR = 6
      IF (N < 5) IERROR = 7
      IF ((NBDCND==5 .OR. NBDCND==6) .AND. (MBDCND==1 .OR. MBDCND==2
     1    .OR. MBDCND==4 .OR. MBDCND==5 .OR. MBDCND==7)) IERROR = 8
      IF (C>0. .AND. NBDCND>=5) IERROR = 9
      IF (IDIMF < M) IERROR = 11
      IF (M < 5) IERROR = 12
      IF(A==0..AND.MBDCND/=5.AND.MBDCND/=6.AND.MBDCND/=9)IERROR=13
      IF (B==PI .AND. MBDCND<=6) IERROR = 14
      IF(A>0..AND.(MBDCND==5.OR.MBDCND==6.OR.MBDCND==9))IERROR=15
      IF (B<PI .AND. MBDCND>=7) IERROR = 16
      IF (ELMBDA/=0. .AND. NBDCND>=5) IERROR = 17
      IF (IERROR == 0) THEN
         IF (INTL == 0) THEN
!     allocate required work space
            K = M + 1
            L = N + 1
            NP = NBDCND
!          compute blktri requirements in irwk,icwk
            CALL BLK_SPACE (N, M, IRWK, ICWK)
!     set work space indices
            IW1 = IRWK + 1
            IWBM = IW1 + M
            IWCM = IWBM + M
            IWAN = IWCM + M
            IWBN = IWAN + N
            IWCN = IWBN + N
            IWSNTH = IWCN + N
            IWRSQ = IWSNTH + M
!     allocate hstcsp required work spac
            IRWK = IWRSQ + N
            ICWK = ICWK + 3*K
            CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
            IF (IERROR == 20) RETURN 
         ENDIF
         IERR1 = 0
      CALL HSTCS1 (INTL,A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +ELMBDA,F,IDIMF,PERTRB,IERR1,w%rew(iw1),w%rew(IWBM),w%rew(IWCM),
     +w%rew(IWAN),w%rew(IWBN),w%rew(IWCN),w%rew(IWSNTH),w%rew(IWRSQ),
     +w%rew,w%cxw)
         IERROR = IERR1
      ENDIF
      RETURN 
      END SUBROUTINE HSTCSP


      SUBROUTINE HSTCS1(INTL, A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND
     1   , BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERR1, AM, BM, CM, AN, BN
     2   , CN, SNTH, RSQ, W, WC)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: INTL
      INTEGER  :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER  :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERR1
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDA(*)
      REAL , INTENT(IN) :: BDB(*)
      REAL , INTENT(IN) :: BDC(*)
      REAL , INTENT(IN) :: BDD(*)
      REAL  :: F(IDIMF,*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL , INTENT(INOUT) :: SNTH(*)
      REAL , INTENT(INOUT) :: RSQ(*)
      REAL  :: W(*)
      COMPLEX  :: WC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J, ISW, NB
      REAL :: DTH, DTHSQ, DR, X, Y, A2, A1, A3
C-----------------------------------------------
      DTH = (B - A)/FLOAT(M)
      DTHSQ = DTH*DTH
      DO I = 1, M
         SNTH(I) = SIN(A + (FLOAT(I) - 0.5)*DTH)
      END DO
      DR = (D - C)/FLOAT(N)
      DO J = 1, N
         RSQ(J) = (C + (FLOAT(J) - 0.5)*DR)**2
      END DO
C
C     MULTIPLY RIGHT SIDE BY R(J)**2
C
      DO J = 1, N
         X = RSQ(J)
         F(:M,J) = X*F(:M,J)
      END DO
C
C      DEFINE COEFFICIENTS AM,BM,CM
C
      X = 1./(2.*COS(DTH/2.))
      AM(2:M) = (SNTH(:M-1)+SNTH(2:M))*X
      CM(:M-1) = AM(2:M)
      AM(1) = SIN(A)
      CM(M) = SIN(B)
      DO I = 1, M
         X = 1./SNTH(I)
         Y = X/DTHSQ
         AM(I) = AM(I)*Y
         CM(I) = CM(I)*Y
         BM(I) = ELMBDA*X*X - AM(I) - CM(I)
      END DO
C
C     DEFINE COEFFICIENTS AN,BN,CN
C
      X = C/DR
      DO J = 1, N
         AN(J) = (X + FLOAT(J - 1))**2
         CN(J) = (X + FLOAT(J))**2
         BN(J) = -(AN(J)+CN(J))
      END DO
      ISW = 1
      NB = NBDCND
      IF (C==0. .AND. NB==2) NB = 6
C
C     ENTER DATA ON THETA BOUNDARIES
C
      GO TO (108,108,110,110,112,112,108,110,112) MBDCND
  108 CONTINUE
      BM(1) = BM(1) - AM(1)
      X = 2.*AM(1)
      F(1,:N) = F(1,:N) - X*BDA(:N)
      GO TO 112
  110 CONTINUE
      BM(1) = BM(1) + AM(1)
      X = DTH*AM(1)
      F(1,:N) = F(1,:N) + X*BDA(:N)
  112 CONTINUE
      GO TO (113,115,115,113,113,115,117,117,117) MBDCND
  113 CONTINUE
      BM(M) = BM(M) - CM(M)
      X = 2.*CM(M)
      F(M,:N) = F(M,:N) - X*BDB(:N)
      GO TO 117
  115 CONTINUE
      BM(M) = BM(M) + CM(M)
      X = DTH*CM(M)
      F(M,:N) = F(M,:N) - X*BDB(:N)
  117 CONTINUE
      GO TO (118,118,120,120,122,122) NB
  118 CONTINUE
      BN(1) = BN(1) - AN(1)
      X = 2.*AN(1)
      F(:M,1) = F(:M,1) - X*BDC(:M)
      GO TO 122
  120 CONTINUE
      BN(1) = BN(1) + AN(1)
      X = DR*AN(1)
      F(:M,1) = F(:M,1) + X*BDC(:M)
  122 CONTINUE
      GO TO (123,125,125,123,123,125) NB
  123 CONTINUE
      BN(N) = BN(N) - CN(N)
      X = 2.*CN(N)
      F(:M,N) = F(:M,N) - X*BDD(:M)
      GO TO 127
  125 CONTINUE
      BN(N) = BN(N) + CN(N)
      X = DR*CN(N)
      F(:M,N) = F(:M,N) - X*BDD(:M)
  127 CONTINUE
      PERTRB = 0.
      GO TO (137,137,128,137,137,128,137,128,128) MBDCND
  128 CONTINUE
      GO TO (137,137,129,137,137,129) NB
  129 CONTINUE
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERR1 = 10
         ELSE
            ISW = 2
            DO I = 1, M
               X = 0.
               X = SUM(F(I,:N))
               PERTRB = PERTRB + X*SNTH(I)
            END DO
            X = 0.
            X = SUM(RSQ(:N))
            PERTRB = 2.*(PERTRB*SIN(DTH/2.))/(X*(COS(A) - COS(B)))
            DO J = 1, N
               X = RSQ(J)*PERTRB
               F(:M,J) = F(:M,J) - X
            END DO
         ENDIF
      ENDIF
  137 CONTINUE
      A2 = SUM(F(:M,1))
      A2 = A2/RSQ(1)
C
C     INITIALIZE BLKTRI
C
      IERR1 = 0
      IF (INTL == 0) CALL BLKTRII (0, 1, N, AN, BN, CN, 1, M, AM, BM, CM
     1   , IDIMF, F, IERR1, W, WC)
      CALL BLKTRII(1,1,N,AN,BN,CN,1,M,AM,BM,CM,IDIMF,F,IERR1,W,WC)
      IF (.NOT.(ISW/=2 .OR. C/=0. .OR. NBDCND/=2)) THEN
         A3 = 0.
         A1 = DOT_PRODUCT(SNTH(:M),F(:M,1))
         A3 = SUM(SNTH(:M))
         A1 = A1 + RSQ(1)*A2/2.
         IF(MBDCND==3)A1=A1+(SIN(B)*BDB(1)-SIN(A)*BDA(1))/(2.*(B-A))
         A1 = A1/A3
         A1 = BDC(1) - A1
         F(:M,:N) = F(:M,:N) + A1
      ENDIF
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
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE HSTCS1
