C
C     file hwscsp.f
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
C     SUBROUTINE HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,
C    +                   BDRS,BDRF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C
C DIMENSION OF           BDTS(N+1),     BDTF(N+1), BDRS(M+1), BDRF(M+1),
C ARGUMENTS              F(IDIMF,N+1)
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION
C                        TO THE MODIFIED HELMHOLTZ EQUATION IN
C                        SPHERICAL COORDINATES ASSUMING AXISYMMETRY
C                        (NO DEPENDENCE ON LONGITUDE).  THE EQUATION
C                        IS
C
C                          (1/R**2)(D/DR)((R**2)(D/DR)U) +
C
C                          (1/(R**2)SIN(THETA))(D/DTHETA)
C
C                          (SIN(THETA)(D/DTHETA)U) +
C
C                          (LAMBDA/(RSIN(THETA))**2)U = F(THETA,R).
C
C                        THIS TWO DIMENSIONAL MODIFIED HELMHOLTZ
C                        EQUATION RESULTS FROM THE FOURIER TRANSFORM
C                        OF THE THREE DIMENSIONAL POISSON EQUATION.
C
C USAGE                  CALL HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,
C                                     RS,RF,N,NBDCND,BDRS,BDRF,ELMBDA,
C                                     F,IDIMF,PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT               INTL
C                          = 0  ON INITIAL ENTRY TO HWSCSP OR IF ANY
C                               OF THE ARGUMENTS RS, RF, N, NBDCND
C                               ARE CHANGED FROM A PREVIOUS CALL.
C                          = 1  IF RS, RF, N, NBDCND ARE ALL UNCHANGED
C                               FROM PREVIOUS CALL TO HWSCSP.
C
C                          NOTE:
C                          A CALL WITH INTL=0 TAKES APPROXIMATELY
C                          1.5 TIMES AS MUCH TIME AS A CALL WITH
C                          INTL = 1  .  ONCE A CALL WITH INTL = 0
C                          HAS BEEN MADE THEN SUBSEQUENT SOLUTIONS
C                          CORRESPONDING TO DIFFERENT F, BDTS, BDTF,
C                          BDRS, BDRF CAN BE OBTAINED FASTER WITH
C                          INTL = 1 SINCE INITIALIZATION IS NOT
C                          REPEATED.
C
C                        TS,TF
C                          THE RANGE OF THETA (COLATITUDE), I.E.,
C                          TS .LE. THETA .LE. TF. TS MUST BE LESS
C                          THAN TF.  TS AND TF ARE IN RADIANS. A TS OF
C                          ZERO CORRESPONDS TO THE NORTH POLE AND A
C                          TF OF PI CORRESPONDS TO THE SOUTH POLE.
C
C                          **** IMPORTANT ****
C
C                          IF TF IS EQUAL TO PI THEN IT MUST BE
C                          COMPUTED USING THE STATEMENT
C                          TF = PIMACH(DUM). THIS INSURES THAT TF
C                          IN THE USER'S PROGRAM IS EQUAL TO PI IN
C                          THIS PROGRAM WHICH PERMITS SEVERAL TESTS
C                          OF THE  INPUT PARAMETERS THAT OTHERWISE
C                          WOULD NOT BE POSSIBLE.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (TS,TF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE M+1 GRID POINTS
C                          IN THE THETA-DIRECTION GIVEN BY
C                          THETA(K) = (I-1)DTHETA+TS FOR
C                          I = 1,2,...,M+1, WHERE DTHETA = (TF-TS)/M
C                          IS THE PANEL WIDTH.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT THETA = TS AND  THETA = TF.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THETA = TF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS AND THETA = TF
C                               (SEE NOTES 1,2 BELOW).
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW) AND
C                               SOLUTION IS SPECIFIED AT THETA = TF.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE SOLUTION IS
C                                SPECIFIED AT THETA = TF.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO THETA
C                               IS SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 7  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE SOLUTION IS
C                                UNSPECIFIED AT THETA = TF = PI.
C                          = 8  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW)
C                               AND THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TF = PI.
C                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THETA = TF = PI.
C
C                          NOTE 1:
C                          IF TS = 0, DO NOT USE MBDCND = 3,4, OR 8,
C                          BUT INSTEAD USE MBDCND = 5,6, OR 9  .
C
C                          NOTE 2:
C                          IF TF = PI, DO NOT USE MBDCND = 2,3, OR 6,
C                          BUT INSTEAD USE MBDCND = 7,8, OR 9  .
C
C                        BDTS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TS.  WHEN MBDCND = 3,4, OR 8,
C
C                            BDTS(J) = (D/DTHETA)U(TS,R(J)),
C                            J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTS IS
C                          A DUMMY VARIABLE.
C
C                        BDTF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TF.  WHEN MBDCND = 2,3, OR 6,
C
C                          BDTF(J) = (D/DTHETA)U(TF,R(J)),
C                          J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTF IS
C                          A DUMMY VARIABLE.
C
C                        RS,RF
C                          THE RANGE OF R, I.E., RS .LE. R .LT. RF.
C                          RS MUST BE LESS THAN RF.  RS MUST BE
C                          NON-NEGATIVE.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (RS,RF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS IN THE
C                          R-DIRECTION GIVEN BY R(J) = (J-1)DR+RS
C                          FOR J = 1,2,...,N+1, WHERE DR = (RF-RS)/N
C                          IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 2
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT R = RS AND R = RF.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               R = RS AND R = RF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               R = RS AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO R
C                               IS SPECIFIED AT R = RF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = RS AND R = RF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               RS AND THE SOLUTION IS SPECIFIED AT
C                               R = RF.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = RS = 0 (SEE NOTE BELOW)  AND THE
C                               SOLUTION IS SPECIFIED AT R = RF.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = RS = 0 (SEE NOTE BELOW) AND THE
C                               DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO R IS SPECIFIED AT R = RF.
C
C                          NOTE:
C                          NBDCND = 5 OR 6 CANNOT BE USED WITH
C                          MBDCND = 1,2,4,5, OR 7.  THE FORMER
C                          INDICATES THAT THE SOLUTION IS UNSPECIFIED
C                          AT R = 0, THE LATTER INDICATES THAT THE
C                          SOLUTION IS SPECIFIED).
C                          USE INSTEAD   NBDCND = 1 OR 2  .
C
C                        BDRS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO R AT R = RS.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDRS(I) = (D/DR)U(THETA(I),RS),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDRS IS
C                          A DUMMY VARIABLE.
C
C                        BDRF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
C                          THAT SPECIFIES THE VALUES OF THE
C                          DERIVATIVE OF THE SOLUTION WITH RESPECT
C                          TO R AT R = RF.
C
C                          WHEN NBDCND = 2,3, OR 6,
C                            BDRF(I) = (D/DR)U(THETA(I),RF),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDRF IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSCSP WILL
C                          ATTEMPT TO FIND A SOLUTION.  IF NBDCND = 5
C                          OR 6 OR  MBDCND = 5,6,7,8, OR 9, ELMBDA
C                          MUST BE ZERO.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
C                          LEAST (M+1)*(N+1), SPECIFYING VALUES OF THE
C                          RIGHT SIDE OF THE HELMHOLTZ EQUATION AND
C                          BOUNDARY VALUES (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(THETA(I),R(J)).
C
C                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
C                          FOR J=1,2,...,N+1,  I=1,2,...,M+1,
C
C                          MBDCND   F(1,J)            F(M+1,J)
C                          ------   ----------        ----------
C
C                            1      U(TS,R(J))        U(TF,R(J))
C                            2      U(TS,R(J))        F(TF,R(J))
C                            3      F(TS,R(J))        F(TF,R(J))
C                            4      F(TS,R(J))        U(TF,R(J))
C                            5      F(0,R(J))         U(TF,R(J))
C                            6      F(0,R(J))         F(TF,R(J))
C                            7      U(TS,R(J))        F(PI,R(J))
C                            8      F(TS,R(J))        F(PI,R(J))
C                            9      F(0,R(J))         F(PI,R(J))
C
C                            NBDCND   F(I,1)            F(I,N+1)
C                            ------   --------------    --------------
C
C                              1      U(THETA(I),RS)    U(THETA(I),RF)
C                              2      U(THETA(I),RS)    F(THETA(I),RF)
C                              3      F(THETA(I),RS)    F(THETA(I),RF)
C                              4      F(THETA(I),RS)    U(THETA(I),RF)
C                              5      F(TS,0)           U(THETA(I),RF)
C                              6      F(TS,0)           F(THETA(I),RF)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
C                          U AND THE RIGHT SIDE F AT A CORNER THEN
C                          THE SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSCSP.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M+1  .
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
c                          user program calling HWSCSP.  The second statement
c                          declares a derived type variable (defined in
c                          the module "fish.f") which is used internally
c                          in HWSCSP to dynamically allocate real and complex
c                          work space used in solution.  An error flag
c                          (IERROR = 20) is set if the required work space
c                          allocation fails (for example if N,M are too large)
c                          Real and complex values are set in the components
c                          of W on a initial (INTL=0) call to HWSCSP.  These
c                          must be preserved on non-initial calls (INTL=1)
c                          to HWSCSP.  This eliminates redundant calculations
c                          and saves compute time.
c               ****       IMPORTANT!  The user program calling HWSCSP should
c                          include the statement:
c
c                               CALL FISHFIN(W)
C
C                          after the final approximation is generated by
C                          HWSCSP.  The will deallocate the real and complex
c                          work space of W.  Failure to include this statement
c                          could result in serious memory leakage.
c
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (THETA(I),R(J)),  I = 1,2,...,M+1,
C                                            J = 1,2,...,N+1  .
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
C                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
C                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
C                          CALCULATED AND SUBTRACTED FROM F, WHICH
C                          ENSURES THAT A SOLUTION EXISTS.  HWSCSP
C                          THEN COMPUTES THIS SOLUTION, WHICH IS A
C                          LEAST SQUARES SOLUTION TO THE ORIGINAL
C                          APPROXIMATION. THIS SOLUTION IS NOT UNIQUE
C                          AND IS UNNORMALIZED. THE VALUE OF PERTRB
C                          SHOULD BE SMALL COMPARED TO THE RIGHT SIDE
C                          F. OTHERWISE , A SOLUTION IS OBTAINED TO
C                          AN ESSENTIALLY DIFFERENT PROBLEM. THIS
C                          COMPARISON SHOULD ALWAYS BE MADE TO INSURE
C                          THAT A MEANINGFUL SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 10,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          = 1  TS.LT.0. OR TF.GT.PI
C                          = 2  TS.GE.TF
C                          = 3  M.LT.5
C                          = 4  MBDCND.LT.1 OR MBDCND.GT.9
C                          = 5  RS.LT.0
C                          = 6  RS.GE.RF
C                          = 7  N.LT.5
C                          = 8  NBDCND.LT.1 OR NBDCND.GT.6
C                          = 9  ELMBDA.GT.0
C                          = 10 IDIMF.LT.M+1
C                          = 11 ELMBDA.NE.0 AND MBDCND.GE.5
C                          = 12 ELMBDA.NE.0 AND NBDCND EQUALS 5 OR 6
C                          = 13 MBDCND EQUALS 5,6 OR 9 AND TS.NE.0
C                          = 14 MBDCND.GE.7 AND TF.NE.PI
C                          = 15 TS.EQ.0 AND MBDCND EQUALS 3,4 OR 8
C                          = 16 TF.EQ.PI AND MBDCND EQUALS 2,3 OR 6
C                          = 17 NBDCND.GE.5 AND RS.NE.0
C                          = 18 NBDCND.GE.5 AND MBDCND EQUALS 1,2,4,5 OR
C                          = 20 If the dynamic allocation of real and
C                               complex work space in the derived type
C                               (fishworkspace) variable W fails (e.g.,
c                               if N,M are too large for the platform used)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSLIBY INCORRECT CALL TO HWSCSP, THE
C                          USER SHOULD TEST IERROR AFTER A CALL.
C
C                        W
c                          The derived type (fishworkspace) variable W
c                          contains real and complex values that must not
C                          be destroyed if HWSCSP is called again with
C                          INTL=1.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED files         fish.f,blktri.f,comf.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980. Revised by John
c                        Adams in June 2004 using Fortran 90 dynamically
C                        allocated work space and derived datat types
c                        to eliminate mixed mode conflicts in the earlier
c                        versions.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
C                        AND THEN CALLS BLKTRI TO SOLVE THE SYSTEM.
C
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, "EFFICIENT
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
C                        ELLIPTIC EQUATIONS"
C                          NCAR TN/IA-109, JULY, 1975, 138 PP.
C***********************************************************************
      SUBROUTINE HWSCSP(INTL, TS, TF, M, MBDCND, BDTS, BDTF, RS, RF, N, 
     1   NBDCND, BDRS, BDRF, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
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
      REAL  :: TS
      REAL  :: TF
      REAL  :: RS
      REAL  :: RF
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDTS(*)
      REAL  :: BDTF(*)
      REAL  :: BDRS(*)
      REAL  :: BDRF(*)
      REAL  :: F(IDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, NCK, L, K, NP
     1   , IRWK, ICWK, NP1, MP1
      REAL :: PI, DUM

      SAVE I1, I2, I3, I4, I5, I6, I7, I8, I9, I10
C-----------------------------------------------
      PI = 4.0*ATAN(1.0)
      IERROR = 0
      IF (TS<0. .OR. TF>PI) IERROR = 1
      IF (TS >= TF) IERROR = 2
      IF (M < 5) IERROR = 3
      IF (MBDCND<1 .OR. MBDCND>9) IERROR = 4
      IF (RS < 0.) IERROR = 5
      IF (RS >= RF) IERROR = 6
      IF (N < 5) IERROR = 7
      IF (NBDCND<1 .OR. NBDCND>6) IERROR = 8
      IF (ELMBDA > 0.) IERROR = 9
      IF (IDIMF < M + 1) IERROR = 10
      IF (ELMBDA/=0. .AND. MBDCND>=5) IERROR = 11
      IF (ELMBDA/=0. .AND. (NBDCND==5 .OR. NBDCND==6)) IERROR = 12
      IF((MBDCND==5.OR.MBDCND==6.OR.MBDCND==9).AND.TS/=0.)IERROR=13
      IF (MBDCND>=7 .AND. TF/=PI) IERROR = 14
      IF(TS==0..AND.(MBDCND==4.OR.MBDCND==8.OR.MBDCND==3))IERROR=15
      IF(TF==PI.AND.(MBDCND==2.OR.MBDCND==3.OR.MBDCND==6))IERROR=16
      IF (NBDCND>=5 .AND. RS/=0.) IERROR = 17
      IF (NBDCND>=5 .AND. (MBDCND==1 .OR. MBDCND==2 .OR. MBDCND==5 .OR. 
     1   MBDCND==7)) IERROR = 18
      IF (IERROR/=0 .AND. IERROR/=9) RETURN 
      NCK = N
      GO TO (101,103,102,103,101,103) NBDCND
  101 CONTINUE
      NCK = NCK - 1
      GO TO 103
  102 CONTINUE
      NCK = NCK + 1
  103 CONTINUE
      L = 2
      K = 1
      L = L + L
      K = K + 1
      DO WHILE(NCK - L > 0)
         L = L + L
         K = K + 1
      END DO
      L = L + L
 
      IF (INTL == 0) THEN
!          compute blktri work space lengths
         NP = NBDCND
         CALL BLK_SPACE (N, M, IRWK, ICWK)
         NP1 = N + 1
         MP1 = M + 1
         I1 = (K - 2)*L + K + MAX0(2*N,6*M) + 13
         I2 = I1 + NP1
         I3 = I2 + NP1
         I4 = I3 + NP1
         I5 = I4 + NP1
         I6 = I5 + NP1
         I7 = I6 + MP1
         I8 = I7 + MP1
         I9 = I8 + MP1
         I10 = I9 + MP1
!          set real and complex work space requirements
         IRWK = I10 + MP1
         ICWK = ICWK + 3*M
!          allocate work space
         CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!          return if allocation fails
         IF (IERROR == 20) RETURN 
      ENDIF
      CALL HWSCS1 (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     +BDRF,ELMBDA,F,IDIMF,PERTRB,w%rew,w%cxw,w%rew(i1),w%rew(i2),
     +w%rew(i3),w%rew(i4),w%rew(i5),w%rew(i6),w%rew(i7),w%rew(i8),
     +w%rew(i9),w%rew(i10),ierror)
      RETURN 
      END SUBROUTINE HWSCSP


      SUBROUTINE HWSCS1(INTL, TS, TF, M, MBDCND, BDTS, BDTF, RS, RF, N, 
     1   NBDCND, BDRS, BDRF, ELMBDA, F, IDIMF, PERTRB, W, WC, S, AN, BN
     2   , CN, R, AM, BM, CM, SINT, BMH, IERROR)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: INTL
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL , INTENT(IN) :: TS
      REAL , INTENT(IN) :: TF
      REAL , INTENT(IN) :: RS
      REAL , INTENT(IN) :: RF
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDTS(*)
      REAL , INTENT(IN) :: BDTF(*)
      REAL , INTENT(IN) :: BDRS(*)
      REAL , INTENT(IN) :: BDRF(*)
      REAL  :: F(IDIMF,*)
      REAL  :: W(*)
      REAL , INTENT(INOUT) :: S(*)
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL , INTENT(INOUT) :: R(*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL , INTENT(INOUT) :: SINT(*)
      REAL , INTENT(INOUT) :: BMH(*)
      COMPLEX  :: WC(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MP1, I, NP1, J, MP, NP, ITS, ITF, ITSP, ITFM, ICTR, JRS
     1   , L, JRF, JRSP, JRFM, MUNK, NUNK, ISING, IFLG
      REAL :: PI, DUM, EPS, DTH, TDT, HDTH, SDTS, THETA, T1, DR, HDR, 
     1   TDR, DR2, CZR, AT, CT, WTS, WTF, AR, WTNM, YPS, CR, WRS, WRF, 
     2   WRZ, SUM, R2, HNE, YHLD, RS2, RF2, RSQ, XP, YPH, XPS
      REAL :: EPMACH
C-----------------------------------------------
      PI = 4.0*ATAN(1.0)
      EPS = EPMACH(DUM)
      MP1 = M + 1
      DTH = (TF - TS)/FLOAT(M)
      TDT = DTH + DTH
      HDTH = DTH/2.
      SDTS = 1./(DTH*DTH)
      DO I = 1, MP1
         THETA = TS + FLOAT(I - 1)*DTH
         SINT(I) = SIN(THETA)
         IF (SINT(I) == 0.) CYCLE 
         T1 = SDTS/SINT(I)
         AM(I) = T1*SIN(THETA - HDTH)
         CM(I) = T1*SIN(THETA + HDTH)
         BM(I) = -(AM(I)+CM(I))
      END DO
      NP1 = N + 1
      DR = (RF - RS)/FLOAT(N)
      HDR = DR/2.
      TDR = DR + DR
      DR2 = DR*DR
      CZR = 6.*DTH/(DR2*(COS(TS) - COS(TF)))
      DO J = 1, NP1
         R(J) = RS + FLOAT(J - 1)*DR
         AN(J) = (R(J)-HDR)**2/DR2
         CN(J) = (R(J)+HDR)**2/DR2
         BN(J) = -(AN(J)+CN(J))
      END DO
      MP = 1
      NP = 1
C
C BOUNDARY CONDITION AT PHI=PS
C
      GO TO (104,104,105,105,106,106,104,105,106) MBDCND
  104 CONTINUE
      AT = AM(2)
      ITS = 2
      GO TO 107
  105 CONTINUE
      AT = AM(1)
      ITS = 1
      CM(1) = CM(1) + AM(1)
      GO TO 107
  106 CONTINUE
      ITS = 1
      BM(1) = -4.*SDTS
      CM(1) = -BM(1)
C
C BOUNDARY CONDITION AT PHI=PF
C
  107 CONTINUE
      GO TO (108,109,109,108,108,109,110,110,110) MBDCND
  108 CONTINUE
      CT = CM(M)
      ITF = M
      GO TO 111
  109 CONTINUE
      CT = CM(M+1)
      AM(M+1) = AM(M+1) + CM(M+1)
      ITF = M + 1
      GO TO 111
  110 CONTINUE
      ITF = M + 1
      AM(M+1) = 4.*SDTS
      BM(M+1) = -AM(M+1)
  111 CONTINUE
      WTS = SINT(ITS+1)*AM(ITS+1)/CM(ITS)
      WTF = SINT(ITF-1)*CM(ITF-1)/AM(ITF)
      ITSP = ITS + 1
      ITFM = ITF - 1
C
C BOUNDARY CONDITION AT R=RS
C
      ICTR = 0
      SELECT CASE (NBDCND) 
      CASE DEFAULT
         AR = AN(2)
         JRS = 2
      CASE (3:4) 
         AR = AN(1)
         JRS = 1
         CN(1) = CN(1) + AN(1)
      CASE (5:6) 
         JRS = 2
         ICTR = 1
         S(N) = AN(N)/BN(N)
         DO J = 3, N
            L = N - J + 2
            S(L) = AN(L)/(BN(L)-CN(L)*S(L+1))
         END DO
         S(2) = -S(2)
         DO J = 3, N
            S(J) = -S(J)*S(J-1)
         END DO
         WTNM = WTS + WTF
         DO I = ITSP, ITFM
            WTNM = WTNM + SINT(I)
         END DO
         YPS = CZR*WTNM*(S(2)-1.)
      END SELECT
C
C BOUNDARY CONDITION AT R=RF
C
  118 CONTINUE
      GO TO (119,120,120,119,119,120) NBDCND
  119 CONTINUE
      CR = CN(N)
      JRF = N
      GO TO 121
  120 CONTINUE
      CR = CN(N+1)
      AN(N+1) = AN(N+1) + CN(N+1)
      JRF = N + 1
  121 CONTINUE
      WRS = AN(JRS+1)*R(JRS)**2/CN(JRS)
      WRF = CN(JRF-1)*R(JRF)**2/AN(JRF)
      WRZ = AN(JRS)/CZR
      JRSP = JRS + 1
      JRFM = JRF - 1
      MUNK = ITF - ITS + 1
      NUNK = JRF - JRS + 1
      BMH(ITS:ITF) = BM(ITS:ITF)
      ISING = 0
      GO TO (132,132,123,132,132,123) NBDCND
  123 CONTINUE
      GO TO (132,132,124,132,132,124,132,124,124) MBDCND
  124 CONTINUE
      IF (ELMBDA >= 0.) THEN
         ISING = 1
         SUM = WTS*WRS + WTS*WRF + WTF*WRS + WTF*WRF
         IF (ICTR /= 0) THEN
            SUM = SUM + WRZ
         ENDIF
         DO J = JRSP, JRFM
            R2 = R(J)**2
            DO I = ITSP, ITFM
               SUM = SUM + R2*SINT(I)
            END DO
         END DO
         DO J = JRSP, JRFM
            SUM = SUM + (WTS + WTF)*R(J)**2
         END DO
         DO I = ITSP, ITFM
            SUM = SUM + (WRS + WRF)*SINT(I)
         END DO
         HNE = SUM
      ENDIF
  132 CONTINUE
      GO TO (133,133,133,133,134,134,133,133,134) MBDCND
  133 CONTINUE
      BM(ITS) = BMH(ITS) + ELMBDA/SINT(ITS)**2
  134 CONTINUE
      GO TO (135,135,135,135,135,135,136,136,136) MBDCND
  135 CONTINUE
      BM(ITF) = BMH(ITF) + ELMBDA/SINT(ITF)**2
  136 CONTINUE
      BM(ITSP:ITFM) = BMH(ITSP:ITFM) + ELMBDA/SINT(ITSP:ITFM)**2
      GO TO (138,138,140,140,142,142,138,140,142) MBDCND
  138 CONTINUE
      F(2,JRS:JRF) = F(2,JRS:JRF) - AT*F(1,JRS:JRF)/R(JRS:JRF)**2
      GO TO 142
  140 CONTINUE
      F(1,JRS:JRF) = F(1,JRS:JRF) + TDT*BDTS(JRS:JRF)*AT/R(JRS:JRF)**2
  142 CONTINUE
      GO TO (143,145,145,143,143,145,147,147,147) MBDCND
  143 CONTINUE
      F(M,JRS:JRF) = F(M,JRS:JRF) - CT*F(M+1,JRS:JRF)/R(JRS:JRF)**2
      GO TO 147
  145 CONTINUE
      F(M+1,JRS:JRF)=F(M+1,JRS:JRF)-TDT*BDTF(JRS:JRF)*CT/R(JRS:JRF)**2
  147 CONTINUE
      SELECT CASE (NBDCND) 
      CASE DEFAULT
         IF (MBDCND - 3 /= 0) GO TO 155
         YHLD = F(ITS,1) - CZR/TDT*(SIN(TF)*BDTF(2)-SIN(TS)*BDTS(2))
         F(:MP1,1) = YHLD
      CASE (1:2) 
         RS2 = (RS + DR)**2
         F(ITS:ITF,2) = F(ITS:ITF,2) - AR*F(ITS:ITF,1)/RS2
      CASE (3:4) 
         F(ITS:ITF,1) = F(ITS:ITF,1) + TDR*BDRS(ITS:ITF)*AR/RS**2
      END SELECT
  155 CONTINUE
      GO TO (156,158,158,156,156,158) NBDCND
  156 CONTINUE
      RF2 = (RF - DR)**2
      F(ITS:ITF,N) = F(ITS:ITF,N) - CR*F(ITS:ITF,N+1)/RF2
      GO TO 160
  158 CONTINUE
      F(ITS:ITF,N+1) = F(ITS:ITF,N+1) - TDR*BDRF(ITS:ITF)*CR/RF**2
  160 CONTINUE
      PERTRB = 0.
      IF (ISING /= 0) THEN
         SUM = WTS*WRS*F(ITS,JRS) + WTS*WRF*F(ITS,JRF) + WTF*WRS*F(ITF,
     1      JRS) + WTF*WRF*F(ITF,JRF)
         IF (ICTR /= 0) THEN
            SUM = SUM + WRZ*F(ITS,1)
         ENDIF
         DO J = JRSP, JRFM
            R2 = R(J)**2
            DO I = ITSP, ITFM
               SUM = SUM + R2*SINT(I)*F(I,J)
            END DO
         END DO
         SUM = SUM + DOT_PRODUCT(R(JRSP:JRFM)**2,WTS*F(ITS,JRSP:JRFM)+
     1      WTF*F(ITF,JRSP:JRFM))
         SUM = SUM + DOT_PRODUCT(SINT(ITSP:ITFM),WRS*F(ITSP:ITFM,JRS)+
     1      WRF*F(ITSP:ITFM,JRF))
         PERTRB = SUM/HNE
         F(:MP1,:NP1) = F(:MP1,:NP1) - PERTRB
      ENDIF
      DO J = JRS, JRF
         RSQ = R(J)**2
         F(ITS:ITF,J) = RSQ*F(ITS:ITF,J)
      END DO
      IFLG = INTL
      CALL BLKTRII (IFLG, NP, NUNK, AN(JRS), BN(JRS), CN(JRS), MP, MUNK
     1   , AM(ITS), BM(ITS), CM(ITS), IDIMF, F(ITS,JRS), IERROR, W, WC)
      IFLG = IFLG + 1
      DO WHILE(IFLG - 1 == 0)
         CALL BLKTRII (IFLG, NP, NUNK, AN(JRS), BN(JRS), CN(JRS), MP, 
     1      MUNK, AM(ITS), BM(ITS), CM(ITS), IDIMF, F(ITS,JRS), IERROR, 
     2      W, WC)
         IFLG = IFLG + 1
      END DO
      IF (NBDCND == 0) THEN
         F(:MP1,JRF+1) = F(:MP1,JRS)
      ENDIF
      IF (MBDCND == 0) THEN
         F(ITF+1,:NP1) = F(ITS,:NP1)
      ENDIF
      XP = 0.
      IF (ICTR /= 0) THEN
         IF (ISING == 0) THEN
            SUM = WTS*F(ITS,2) + WTF*F(ITF,2)
            SUM = SUM + DOT_PRODUCT(SINT(ITSP:ITFM),F(ITSP:ITFM,2))
            YPH = CZR*SUM
            XP = (F(ITS,1)-YPH)/YPS
            DO J = JRS, JRF
               XPS = XP*S(J)
               F(ITS:ITF,J) = F(ITS:ITF,J) + XPS
            END DO
         ENDIF
         F(:MP1,1) = XP
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
c June      2004    Version 5.0, Fortran 90 Changes
C-----------------------------------------------------------------------
      END SUBROUTINE HWSCS1
