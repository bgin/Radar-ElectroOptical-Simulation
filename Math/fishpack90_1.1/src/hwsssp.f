C
C     file hwsssp.f
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
C     SUBROUTINE HWSSSP (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,
C    +                   BDPF,ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C DIMENSION OF           BDTS(N+1),    BDTF(N+1), BDPS(M+1), BDPF(M+1),
C ARGUMENTS              F(IDIMF,N+1)
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION TO
C                        THE HELMHOLTZ EQUATION IN SPHERICAL
C                        COORDINATES AND ON THE SURFACE OF THE UNIT
C                        SPHERE (RADIUS OF 1).  THE EQUATION IS
C
C                          (1/SIN(THETA))(D/DTHETA)(SIN(THETA)
C                          (DU/DTHETA)) + (1/SIN(THETA)**2)(D/DPHI)
C                          (DU/DPHI)  + LAMBDA*U = F(THETA,PHI)
C
C                        WHERE THETA IS COLATITUDE AND PHI IS
C                        LONGITUDE.
C
C USAGE                  CALL HWSSSP (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,
C                                     N,NBDCND,BDPS,BDPF,ELMBDA,F,
C                                     IDIMF,PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT               TS,TF
C
C                          THE RANGE OF THETA (COLATITUDE), I.E.,
C                          TS .LE. THETA .LE. TF. TS MUST BE LESS
C                          THAN TF.  TS AND TF ARE IN RADIANS.
C                          A TS OF ZERO CORRESPONDS TO THE NORTH
C                          POLE AND A TF OF PI CORRESPONDS TO
C                          THE SOUTH POLE.
C
C                          * * * IMPORTANT * * *
C
C                          IF TF IS EQUAL TO PI THEN IT MUST BE
C                          COMPUTED USING THE STATEMENT
C                          TF = PIMACH(DUM). THIS INSURES THAT TF
C                          IN THE USER'S PROGRAM IS EQUAL TO PI IN
C                          THIS PROGRAM WHICH PERMITS SEVERAL TESTS
C                          OF THE INPUT PARAMETERS THAT OTHERWISE
C                          WOULD NOT BE POSSIBLE.
C
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (TS,TF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE M+1 GRID POINTS IN THE
C                          THETA-DIRECTION GIVEN BY
C                          THETA(I) = (I-1)DTHETA+TS FOR
C                          I = 1,2,...,M+1, WHERE
C                          DTHETA = (TF-TS)/M IS THE PANEL WIDTH.
C                          M MUST BE GREATER THAN 5
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT THETA = TS AND THETA = TF.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THETA = TF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE DERIVATIVE OF
C                               THE SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               SPECIFIED AT THETA = TS AND
C                               THETA = TF (SEE NOTES 1,2 BELOW).
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW)
C                               AND THE SOLUTION IS SPECIFIED AT
C                               THETA = TF.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE SOLUTION
C                               IS SPECIFIED AT THETA = TF.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO THETA
C                               IS SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 7  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE SOLUTION IS
C                               IS UNSPECIFIED AT THETA = TF = PI.
C                          = 8  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW) AND
C                               THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TF = PI.
C                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THETA = TF = PI.
C
C                          NOTES:
C                          IF TS = 0, DO NOT USE MBDCND = 3,4, OR 8,
C                          BUT INSTEAD USE MBDCND = 5,6, OR 9  .
C
C                          IF TF = PI, DO NOT USE MBDCND = 2,3, OR 6,
C                          BUT INSTEAD USE MBDCND = 7,8, OR 9  .
C
C                        BDTS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TS.  WHEN MBDCND = 3,4, OR 8,
C
C                          BDTS(J) = (D/DTHETA)U(TS,PHI(J)),
C                          J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTS IS
C                          A DUMMY VARIABLE.
C
C                        BDTF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
C                          THAT SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TF.  WHEN MBDCND = 2,3, OR 6,
C
C                          BDTF(J) = (D/DTHETA)U(TF,PHI(J)),
C                          J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTF IS
C                          A DUMMY VARIABLE.
C
C                        PS,PF
C                          THE RANGE OF PHI (LONGITUDE), I.E.,
C                          PS .LE. PHI .LE. PF.  PS MUST BE LESS
C                          THAN PF.  PS AND PF ARE IN RADIANS.
C                          IF PS = 0 AND PF = 2*PI, PERIODIC
C                          BOUNDARY CONDITIONS ARE USUALLY PRESCRIBED.
C
C                          * * * IMPORTANT * * *
C
C                          IF PF IS EQUAL TO 2*PI THEN IT MUST BE
C                          COMPUTED USING THE STATEMENT
C                          PF = 2.*PIMACH(DUM). THIS INSURES THAT
C                          PF IN THE USERS PROGRAM IS EQUAL TO
C                          2*PI IN THIS PROGRAM WHICH PERMITS TESTS
C                          OF THE INPUT PARAMETERS THAT OTHERWISE
C                          WOULD NOT BE POSSIBLE.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (PS,PF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS
C                          IN THE PHI-DIRECTION GIVEN BY
C                          PHI(J) = (J-1)DPHI+PS  FOR
C                          J = 1,2,...,N+1, WHERE
C                          DPHI = (PF-PS)/N IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 4
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT PHI = PS AND PHI = PF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN PHI,
C                               I.U., U(I,J) = U(I,N+J).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               PHI = PS AND PHI = PF
C                               (SEE NOTE BELOW).
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               PHI = PS (SEE NOTE BELOW)
C                               AND THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PHI = PF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PHI = PS AND PHI = PF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PS AND THE SOLUTION IS SPECIFIED
C                               AT PHI = PF
C
C                          NOTE:
C                          NBDCND = 1,2, OR 4 CANNOT BE USED WITH
C                          MBDCND = 5,6,7,8, OR 9.  THE FORMER INDICATES
C                          THAT THE SOLUTION IS SPECIFIED AT A POLE, THE
C                          LATTER INDICATES THAT THE SOLUTION IS NOT
C                          SPECIFIED.  USE INSTEAD  MBDCND = 1 OR 2.
C
C                        BDPS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO PHI AT
C                          PHI = PS.  WHEN NBDCND = 3 OR 4,
C
C                            BDPS(I) = (D/DPHI)U(THETA(I),PS),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDPS IS
C                          A DUMMY VARIABLE.
C
C                        BDPF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO PHI AT
C                          PHI = PF.  WHEN NBDCND = 2 OR 3,
C
C                            BDPF(I) = (D/DPHI)U(THETA(I),PF),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDPF IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSSSP WILL
C                          ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUE OF THE RIGHT SIDE OF THE HELMHOLTZ
C                          EQUATION AND BOUNDARY VALUES (IF ANY).
C                          F MUST BE DIMENSIONED AT LEAST (M+1)*(N+1).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(THETA(I),PHI(J)).
C
C                          ON THE BOUNDARIES F IS DEFINED AS FOLLOWS:
C                          FOR J = 1,2,...,N+1 AND I = 1,2,...,M+1
C
C                          MBDCND   F(1,J)            F(M+1,J)
C                          ------   ------------      ------------
C
C                            1      U(TS,PHI(J))      U(TF,PHI(J))
C                            2      U(TS,PHI(J))      F(TF,PHI(J))
C                            3      F(TS,PHI(J))      F(TF,PHI(J))
C                            4      F(TS,PHI(J))      U(TF,PHI(J))
C                            5      F(0,PS)           U(TF,PHI(J))
C                            6      F(0,PS)           F(TF,PHI(J))
C                            7      U(TS,PHI(J))      F(PI,PS)
C                            8      F(TS,PHI(J))      F(PI,PS)
C                            9      F(0,PS)           F(PI,PS)
C
C                          NBDCND   F(I,1)            F(I,N+1)
C                          ------   --------------    --------------
C
C                            0      F(THETA(I),PS)    F(THETA(I),PS)
C                            1      U(THETA(I),PS)    U(THETA(I),PF)
C                            2      U(THETA(I),PS)    F(THETA(I),PF)
C                            3      F(THETA(I),PS)    F(THETA(I),PF)
C                            4      F(THETA(I),PS)    U(THETA(I),PF)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION U
C                          AND THE RIGHT SIDE F AT A CORNER THEN THE
C                          SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSSSP.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F. IDIMF MUST BE
C                          AT LEAST M+1  .
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (THETA(I),PHI(J)),  I = 1,2,...,M+1  AND
C                          J = 1,2,...,N+1  .
C
C                        PERTRB
C                          IF ONE SPECIFIES A COMBINATION OF PERIODIC,
C                          DERIVATIVE OR UNSPECIFIED BOUNDARY
C                          CONDITIONS FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HWSSSP THEN COMPUTES
C                          THIS SOLUTION, WHICH IS A LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION IS NOT UNIQUE AND IS
C                          UNNORMALIZED. THE VALUE OF PERTRB SHOULD
C                          BE SMALL COMPARED TO THE RIGHT SIDE F.
C                          OTHERWISE , A SOLUTION IS OBTAINED TO AN
C                          ESSENTIALLY DIFFERENT PROBLEM. THIS
C                          COMPARISON SHOULD ALWAYS BE MADE TO INSURE
C                          THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 8,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          = 0  NO ERROR
C                          = 1  TS.LT.0 OR TF.GT.PI
C                          = 2  TS.GE.TF
C                          = 3  MBDCND.LT.1 OR MBDCND.GT.9
C                          = 4  PS.LT.0 OR PS.GT.PI+PI
C                          = 5  PS.GE.PF
C                          = 6  N.LT.5
C                          = 7  M.LT.5
C                          = 8  NBDCND.LT.0 OR NBDCND.GT.4
C                          = 9  ELMBDA.GT.0
C                          = 10 IDIMF.LT.M+1
C                          = 11 NBDCND EQUALS 1,2 OR 4 AND MBDCND.GE.5
C                          = 12 TS.EQ.0 AND MBDCND EQUALS 3,4 OR 8
C                          = 13 TF.EQ.PI AND MBDCND EQUALS 2,3 OR 6
C                          = 14 MBDCND EQUALS 5,6 OR 9 AND TS.NE.0
C                          = 15 MBDCND.GE.7 AND TF.NE.PI
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED files         fish.f,genbun.f,gnbnaux.f,comf.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
C                        AND THEN CALLS GENBUN TO SOLVE THE SYSTEM.
C
C TIMING                 FOR LARGE  M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO
C                          M*N*(LOG2(N)
C                        BUT ALSO DEPENDS ON INPUT PARAMETERS NBDCND
C                        AND MBDCND.
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A LOSS
C                        OF NO MORE THAN THREE SIGNIFICANT DIGITS FOR N
C                        AND M AS LARGE AS 64.  MORE DETAILS ABOUT
C                        ACCURACY CAN BE FOUND IN THE DOCUMENTATION FOR
C                        SUBROUTINE GENBUN WHICH IS THE ROUTINE THAT
C                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
C
C REFERENCES             P. N. SWARZTRAUBER, "THE DIRECT SOLUTION OF
C                        THE DISCRETE POISSON EQUATION ON THE SURFACE OF
C                        A SPHERE", S.I.A.M. J. NUMER. ANAL.,15(1974),
C                        PP 212-215.
C
C                        SWARZTRAUBER,P. AND R. SWEET, "EFFICIENT
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
C                        ELLIPTIC EQUATIONS", NCAR TN/IA-109, JULY,
C                        1975, 138 PP.
C***********************************************************************
      SUBROUTINE HWSSSP(TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, NBDCND
     1   , BDPS, BDPF, ELMBDA, F, IDIMF, PERTRB, IERROR)
      USE fish
      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: TS
      REAL  :: TF
      REAL  :: PS
      REAL  :: PF
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDTS(*)
      REAL  :: BDTF(*)
      REAL  :: BDPS(*)
      REAL  :: BDPF(*)
      REAL  :: F(IDIMF,1)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NBR, IRWK, ICWK
      REAL :: PI, DUM, TPI
C-----------------------------------------------
C
      NBR = NBDCND + 1
      PI = 4.0*ATAN(1.0)
      TPI = 2.*PI
      IERROR = 0
      IF (TS<0. .OR. TF>PI) IERROR = 1
      IF (TS >= TF) IERROR = 2
      IF (MBDCND<1 .OR. MBDCND>9) IERROR = 3
      IF (PS<0. .OR. PF>TPI) IERROR = 4
      IF (PS >= PF) IERROR = 5
      IF (N < 5) IERROR = 6
      IF (M < 5) IERROR = 7
      IF (NBDCND<0 .OR. NBDCND>4) IERROR = 8
      IF (ELMBDA > 0.) IERROR = 9
      IF (IDIMF < M + 1) IERROR = 10
      IF ((NBDCND==1 .OR. NBDCND==2 .OR. NBDCND==4) .AND. MBDCND>=5) 
     1   IERROR = 11
      IF(TS==0..AND.(MBDCND==3.OR.MBDCND==4.OR.MBDCND==8))IERROR=12
      IF(TF==PI.AND.(MBDCND==2.OR.MBDCND==3.OR.MBDCND==6))IERROR=13
      IF((MBDCND==5.OR.MBDCND==6.OR.MBDCND==9).AND.TS/=0.)IERROR=14
      IF (MBDCND>=7 .AND. TF/=PI) IERROR = 15
      IF (IERROR/=0 .AND. IERROR/=9) RETURN 
!     allocate generous work space estimate
      IRWK=4*(N+1)+(16+INT(ALOG(FLOAT(N+1))/ALOG(2.0)))*(M+1)
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hwssspp(TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,
     +             BDPF,ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release dynamically allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HWSSSP


 
      SUBROUTINE HWSSSPP(TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, 
     1   NBDCND, BDPS, BDPF, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      INTEGER  :: IERROR
      REAL  :: TS
      REAL  :: TF
      REAL  :: PS
      REAL  :: PF
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDTS(*)
      REAL  :: BDTF(*)
      REAL  :: BDPS(*)
      REAL  :: BDPF(*)
      REAL  :: F(IDIMF,*)
      REAL  :: W(*)
C-----------------------------------------------
      CALL HWSSS1 (TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, NBDCND, 
     1   BDPS, BDPF, ELMBDA, F, IDIMF, PERTRB, W, W(M+2), W(2*M+3), W(3*
     2   M+4), W(4*M+5), W(5*M+6), W(6*M+7))
      RETURN 
      END SUBROUTINE HWSSSPP


 
      SUBROUTINE HWSSS1(TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, NBDCND
     1   , BDPS, BDPF, ELMBDA, F, IDIMF, PERTRB, AM, BM, CM, SN, SS, 
     2   SINT, D)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER  :: NBDCND
      INTEGER  :: IDIMF
      REAL , INTENT(IN) :: TS
      REAL , INTENT(IN) :: TF
      REAL , INTENT(IN) :: PS
      REAL , INTENT(IN) :: PF
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDTS(*)
      REAL , INTENT(IN) :: BDTF(*)
      REAL , INTENT(IN) :: BDPS(*)
      REAL , INTENT(IN) :: BDPF(*)
      REAL  :: F(IDIMF,*)
      REAL  :: AM(*)
      REAL  :: BM(*)
      REAL  :: CM(*)
      REAL , INTENT(INOUT) :: SN(*)
      REAL , INTENT(INOUT) :: SS(*)
      REAL , INTENT(INOUT) :: SINT(*)
      REAL  :: D(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MP1, NP1, I, INP, ISP, MBR, ITS, ITF, ITSP, ITFM, MUNK
     1   , IID, II, NBR, JPS, JPF, JPSP, JPFM, NUNK, ISING, J, IERROR
      REAL :: PI, DUM, TPI, HPI, FN, FM, DTH, HDTH, TDT, DPHI, TDP, 
     1   DPHI2, EDP2, DTH2, CP, WP, FIM1, THETA, T1, AT, CT, WTS, WTF, 
     2   WPS, WPF, FJJ, CF, SUM, SUM1, HNE, YHLD, SUM2, DFN, DNN, DSN, 
     3   CNP, HLD, DFS, DSS, DNS, CSP, RTN, RTS, DEN
C-----------------------------------------------
C
      PI = 4.0*ATAN(1.0)
      TPI = PI + PI
      HPI = PI/2.
      MP1 = M + 1
      NP1 = N + 1
      FN = N
      FM = M
      DTH = (TF - TS)/FM
      HDTH = DTH/2.
      TDT = DTH + DTH
      DPHI = (PF - PS)/FN
      TDP = DPHI + DPHI
      DPHI2 = DPHI*DPHI
      EDP2 = ELMBDA*DPHI2
      DTH2 = DTH*DTH
      CP = 4./(FN*DTH2)
      WP = FN*SIN(HDTH)/4.
      DO I = 1, MP1
         FIM1 = I - 1
         THETA = FIM1*DTH + TS
         SINT(I) = SIN(THETA)
         IF (SINT(I) == 0.) CYCLE 
         T1 = 1./(DTH2*SINT(I))
         AM(I) = T1*SIN(THETA - HDTH)
         CM(I) = T1*SIN(THETA + HDTH)
         BM(I) = (-AM(I)) - CM(I) + ELMBDA
      END DO
      INP = 0
      ISP = 0
C
C BOUNDARY CONDITION AT THETA=TS
C
      MBR = MBDCND + 1
      GO TO (103,104,104,105,105,106,106,104,105,106) MBR
  103 CONTINUE
      ITS = 1
      GO TO 107
  104 CONTINUE
      AT = AM(2)
      ITS = 2
      GO TO 107
  105 CONTINUE
      AT = AM(1)
      ITS = 1
      CM(1) = AM(1) + CM(1)
      GO TO 107
  106 CONTINUE
      AT = AM(2)
      INP = 1
      ITS = 2
C
C BOUNDARY CONDITION THETA=TF
C
  107 CONTINUE
      GO TO (108,109,110,110,109,109,110,111,111,111) MBR
  108 CONTINUE
      ITF = M
      GO TO 112
  109 CONTINUE
      CT = CM(M)
      ITF = M
      GO TO 112
  110 CONTINUE
      CT = CM(M+1)
      AM(M+1) = AM(M+1) + CM(M+1)
      ITF = M + 1
      GO TO 112
  111 CONTINUE
      ITF = M
      ISP = 1
      CT = CM(M)
C
C COMPUTE HOMOGENEOUS SOLUTION WITH SOLUTION AT POLE EQUAL TO ONE
C
  112 CONTINUE
      ITSP = ITS + 1
      ITFM = ITF - 1
      WTS = SINT(ITS+1)*AM(ITS+1)/CM(ITS)
      WTF = SINT(ITF-1)*CM(ITF-1)/AM(ITF)
      MUNK = ITF - ITS + 1
      IF (ISP > 0) THEN
         D(ITS) = CM(ITS)/BM(ITS)
         DO I = ITSP, M
            D(I) = CM(I)/(BM(I)-AM(I)*D(I-1))
         END DO
         SS(M) = -D(M)
         IID = M - ITS
         DO II = 1, IID
            I = M - II
            SS(I) = -D(I)*SS(I+1)
         END DO
         SS(M+1) = 1.
      ENDIF
      IF (INP > 0) THEN
         SN(1) = 1.
         D(ITF) = AM(ITF)/BM(ITF)
         IID = ITF - 2
         DO II = 1, IID
            I = ITF - II
            D(I) = AM(I)/(BM(I)-CM(I)*D(I+1))
         END DO
         SN(2) = -D(2)
         DO I = 3, ITF
            SN(I) = -D(I)*SN(I-1)
         END DO
      ENDIF
C
C BOUNDARY CONDITIONS AT PHI=PS
C
      NBR = NBDCND + 1
      WPS = 1.
      WPF = 1.
      SELECT CASE (NBR) 
      CASE DEFAULT
         JPS = 1
      CASE (2:3) 
         JPS = 2
      CASE (4:5) 
         JPS = 1
         WPS = 0.5
      END SELECT
C
C BOUNDARY CONDITION AT PHI=PF
C
  124 CONTINUE
      GO TO (125,126,127,127,126) NBR
  125 CONTINUE
      JPF = N
      GO TO 128
  126 CONTINUE
      JPF = N
      GO TO 128
  127 CONTINUE
      WPF = 0.5
      JPF = N + 1
  128 CONTINUE
      JPSP = JPS + 1
      JPFM = JPF - 1
      NUNK = JPF - JPS + 1
      FJJ = JPFM - JPSP + 1
C
C SCALE COEFFICIENTS FOR SUBROUTINE GENBUN
C
      DO I = ITS, ITF
         CF = DPHI2*SINT(I)*SINT(I)
         AM(I) = CF*AM(I)
         BM(I) = CF*BM(I)
         CM(I) = CF*CM(I)
      END DO
      AM(ITS) = 0.
      CM(ITF) = 0.
      ISING = 0
      GO TO (130,138,138,130,138,138,130,138,130,130) MBR
  130 CONTINUE
      GO TO (131,138,138,131,138) NBR
  131 CONTINUE
      IF (ELMBDA >= 0.) THEN
         ISING = 1
         SUM = WTS*WPS + WTS*WPF + WTF*WPS + WTF*WPF
         IF (INP > 0) THEN
            SUM = SUM + WP
         ENDIF
         IF (ISP > 0) THEN
            SUM = SUM + WP
         ENDIF
         SUM1 = 0.
         DO I = ITSP, ITFM
            SUM1 = SUM1 + SINT(I)
         END DO
         SUM = SUM + FJJ*(SUM1 + WTS + WTF)
         SUM = SUM + (WPS + WPF)*SUM1
         HNE = SUM
      ENDIF
  138 CONTINUE
      GO TO (146,142,142,144,144,139,139,142,144,139) MBR
  139 CONTINUE
      IF (NBDCND - 3 /= 0) GO TO 146
      YHLD = F(1,JPS) - 4./(FN*DPHI*DTH2)*(BDPF(2)-BDPS(2))
      F(1,:NP1) = YHLD
      GO TO 146
  142 CONTINUE
      F(2,JPS:JPF) = F(2,JPS:JPF) - AT*F(1,JPS:JPF)
      GO TO 146
  144 CONTINUE
      F(1,JPS:JPF) = F(1,JPS:JPF) + TDT*BDTS(JPS:JPF)*AT
  146 CONTINUE
      GO TO (154,150,152,152,150,150,152,147,147,147) MBR
  147 CONTINUE
      IF (NBDCND - 3 /= 0) GO TO 154
      YHLD = F(M+1,JPS) - 4./(FN*DPHI*DTH2)*(BDPF(M)-BDPS(M))
      F(M+1,:NP1) = YHLD
      GO TO 154
  150 CONTINUE
      F(M,JPS:JPF) = F(M,JPS:JPF) - CT*F(M+1,JPS:JPF)
      GO TO 154
  152 CONTINUE
      F(M+1,JPS:JPF) = F(M+1,JPS:JPF) - TDT*BDTF(JPS:JPF)*CT
  154 CONTINUE
      GO TO (159,155,155,157,157) NBR
  155 CONTINUE
      F(ITS:ITF,2) = F(ITS:ITF,2) - F(ITS:ITF,1)/(DPHI2*SINT(ITS:ITF)*
     1   SINT(ITS:ITF))
      GO TO 159
  157 CONTINUE
      F(ITS:ITF,1) = F(ITS:ITF,1) + TDP*BDPS(ITS:ITF)/(DPHI2*SINT(ITS:
     1   ITF)*SINT(ITS:ITF))
  159 CONTINUE
      GO TO (164,160,162,162,160) NBR
  160 CONTINUE
      F(ITS:ITF,N) = F(ITS:ITF,N) - F(ITS:ITF,N+1)/(DPHI2*SINT(ITS:ITF)*
     1   SINT(ITS:ITF))
      GO TO 164
  162 CONTINUE
      F(ITS:ITF,N+1) = F(ITS:ITF,N+1) - TDP*BDPF(ITS:ITF)/(DPHI2*SINT(
     1   ITS:ITF)*SINT(ITS:ITF))
  164 CONTINUE
      PERTRB = 0.
      IF (ISING /= 0) THEN
         SUM = WTS*WPS*F(ITS,JPS) + WTS*WPF*F(ITS,JPF) + WTF*WPS*F(ITF,
     1      JPS) + WTF*WPF*F(ITF,JPF)
         IF (INP > 0) THEN
            SUM = SUM + WP*F(1,JPS)
         ENDIF
         IF (ISP > 0) THEN
            SUM = SUM + WP*F(M+1,JPS)
         ENDIF
         DO I = ITSP, ITFM
            SUM1 = 0.
            DO J = JPSP, JPFM
               SUM1 = SUM1 + F(I,J)
            END DO
            SUM = SUM + SINT(I)*SUM1
         END DO
         SUM1 = 0.
         SUM2 = 0.
         DO J = JPSP, JPFM
            SUM1 = SUM1 + F(ITS,J)
            SUM2 = SUM2 + F(ITF,J)
         END DO
         SUM = SUM + WTS*SUM1 + WTF*SUM2
         SUM1 = 0.
         SUM2 = 0.
         SUM1 = DOT_PRODUCT(SINT(ITSP:ITFM),F(ITSP:ITFM,JPS))
         SUM2 = DOT_PRODUCT(SINT(ITSP:ITFM),F(ITSP:ITFM,JPF))
         SUM = SUM + WPS*SUM1 + WPF*SUM2
         PERTRB = SUM/HNE
         F(:MP1,:NP1) = F(:MP1,:NP1) - PERTRB
      ENDIF
C
C SCALE RIGHT SIDE FOR SUBROUTINE GENBUN
C
      DO I = ITS, ITF
         CF = DPHI2*SINT(I)*SINT(I)
         F(I,JPS:JPF) = CF*F(I,JPS:JPF)
      END DO
      CALL GENBUNN (NBDCND, NUNK, 1, MUNK, AM(ITS), BM(ITS), CM(ITS), 
     1   IDIMF, F(ITS,JPS), IERROR, D)
      IF (ISING <= 0) GO TO 186
      IF (INP > 0) THEN
         IF (ISP > 0) GO TO 186
         F(1,:NP1) = 0.
         GO TO 209
      ENDIF
      IF (ISP <= 0) GO TO 186
      F(M+1,:NP1) = 0.
      GO TO 209
  186 CONTINUE
      IF (INP > 0) THEN
         SUM = WPS*F(ITS,JPS) + WPF*F(ITS,JPF)
         DO J = JPSP, JPFM
            SUM = SUM + F(ITS,J)
         END DO
         DFN = CP*SUM
         DNN = CP*((WPS + WPF + FJJ)*(SN(2)-1.)) + ELMBDA
         DSN = CP*(WPS + WPF + FJJ)*SN(M)
         IF (ISP > 0) GO TO 194
         CNP = (F(1,1)-DFN)/DNN
         DO I = ITS, ITF
            HLD = CNP*SN(I)
            F(I,JPS:JPF) = F(I,JPS:JPF) + HLD
         END DO
         F(1,:NP1) = CNP
         GO TO 209
      ENDIF
      IF (ISP <= 0) GO TO 209
  194 CONTINUE
      SUM = WPS*F(ITF,JPS) + WPF*F(ITF,JPF)
      DO J = JPSP, JPFM
         SUM = SUM + F(ITF,J)
      END DO
      DFS = CP*SUM
      DSS = CP*((WPS + WPF + FJJ)*(SS(M)-1.)) + ELMBDA
      DNS = CP*(WPS + WPF + FJJ)*SS(2)
      IF (INP <= 0) THEN
         CSP = (F(M+1,1)-DFS)/DSS
         DO I = ITS, ITF
            HLD = CSP*SS(I)
            F(I,JPS:JPF) = F(I,JPS:JPF) + HLD
         END DO
         F(M+1,:NP1) = CSP
      ELSE
         RTN = F(1,1) - DFN
         RTS = F(M+1,1) - DFS
         IF (ISING > 0) THEN
            CSP = 0.
            CNP = RTN/DNN
         ELSE
            IF (ABS(DNN) - ABS(DSN) > 0.) THEN
               DEN = DSS - DNS*DSN/DNN
               RTS = RTS - RTN*DSN/DNN
               CSP = RTS/DEN
               CNP = (RTN - CSP*DNS)/DNN
            ELSE
               DEN = DNS - DSS*DNN/DSN
               RTN = RTN - RTS*DNN/DSN
               CSP = RTN/DEN
               CNP = (RTS - DSS*CSP)/DSN
            ENDIF
         ENDIF
         DO I = ITS, ITF
            HLD = CNP*SN(I) + CSP*SS(I)
            F(I,JPS:JPF) = F(I,JPS:JPF) + HLD
         END DO
         F(1,:NP1) = CNP
         F(M+1,:NP1) = CSP
      ENDIF
  209 CONTINUE
      IF (NBDCND == 0) THEN
         F(:MP1,JPF+1) = F(:MP1,JPS)
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
      END SUBROUTINE HWSSS1
