C
C     file hstssp.f
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
C     SUBROUTINE HSTSSP (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION ON A STAGGERED GRID
C                        TO THE HELMHOLTZ EQUATION IN SPHERICAL
C                        COORDINATES AND ON THE SURFACE OF THE UNIT
C                        SPHERE (RADIUS OF 1).  THE EQUATION IS
C
C                          (1/SIN(THETA))(D/DTHETA)(SIN(THETA)
C                          (DU/DTHETA)) + (1/SIN(THETA)**2)
C                          (D/DPHI)(DU/DPHI) + LAMBDA*U = F(THETA,PHI)
C
C                        WHERE THETA IS COLATITUDE AND PHI IS
C                        LONGITUDE.
C
C USAGE                  CALL HSTSSP (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR)
C
C
C ARGUMENTS
C ON INPUT
C
C                        A,B
C                          THE RANGE OF THETA (COLATITUDE),
C                          I.E. A .LE. THETA .LE. B.
C                          A MUST BE LESS THAN B AND A MUST BE
C                          NON-NEGATIVE.  A AND B ARE IN RADIANS.
C                          A = 0 CORRESPONDS TO THE NORTH POLE AND
C                          B = PI CORRESPONDS TO THE SOUTH POLE.
C
C
C                            * * *  IMPORTANT  * * *
C
C                          IF B IS EQUAL TO PI, THEN B MUST BE
C                          COMPUTED USING THE STATEMENT
C                            B = PIMACH(DUM)
C
C                          THIS INSURES THAT B IN THE USER"S PROGRAM
C                          IS EQUAL TO PI IN THIS PROGRAM WHICH
C                          PERMITS SEVERAL TESTS OF THE INPUT
C                          PARAMETERS THAT OTHERWISE WOULD NOT BE
C                          POSSIBLE.
C
C                            * * * * * * * * * * * *
C                        M
C                          THE NUMBER OF GRID POINTS IN THE INTERVAL
C                          (A,B).  THE GRID POINTS IN THE THETA
C                          DIRECTION ARE GIVEN BY
C                          THETA(I) = A + (I-0.5)DTHETA
C                          FOR I=1,2,...,M WHERE DTHETA =(B-A)/M.
C                          M MUST BE GREATER THAN 2.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT THETA = A AND THETA = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THETA = B.
C                               (SEE NOTE 3 BELOW)
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = B
C                               (SEE NOTES 2 AND 3 BELOW).
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                             WITH RESPECT TO THETA IS SPECIFIED
C                             AT THETA = A
C                             (SEE NOTES 1, 2 BELOW) AND THETA = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = A
C                               (SEE NOTES 1 AND 2 BELOW) AND THE
C                               SOLUTION IS SPECIFIED AT THETA = B.
C
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = A = 0 AND THE SOLUTION IS
C                               SPECIFIED AT THETA = B.
C                               (SEE NOTE 3 BELOW)
C
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = A = 0 AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO THETA
C                               IS SPECIFIED AT THETA = B
C                               (SEE NOTE 2 BELOW).
C
C                          = 7  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = A AND THE SOLUTION IS
C                               UNSPECIFIED AT THETA = B = PI.
C                               (SEE NOTE 3 BELOW)
C
C                          = 8  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED AT
C                               THETA = A (SEE NOTE 1 BELOW)
C                               AND THE SOLUTION IS UNSPECIFIED AT
C                               THETA = B = PI.
C
C                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = A = 0 AND THETA = B = PI.
C
C                          NOTE 1:
C                          IF A = 0, DO NOT USE MBDCND = 3, 4, OR 8,
C                          BUT INSTEAD USE MBDCND = 5, 6, OR 9.
C
C                          NOTE 2:
C                          IF B = PI, DO NOT USE MBDCND = 2, 3, OR 6,
C                          BUT INSTEAD USE MBDCND = 7, 8, OR 9.
C
C                          NOTE 3:
C                          WHEN THE SOLUTION IS SPECIFIED AT
C                          THETA = 0 AND/OR THETA = PI AND THE OTHER
C                          BOUNDARY CONDITIONS ARE COMBINATIONS
C                          OF UNSPECIFIED, NORMAL DERIVATIVE, OR
C                          PERIODICITY A SINGULAR SYSTEM RESULTS.
C                          THE UNIQUE SOLUTION IS DETERMINED BY
C                          EXTRAPOLATION TO THE SPECIFICATION OF THE
C                          SOLUTION AT EITHER THETA = 0 OR THETA = PI.
C                          BUT IN THESE CASES THE RIGHT SIDE OF THE
C                          SYSTEM  WILL BE PERTURBED BY THE CONSTANT
C                          PERTRB.
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES (IF ANY) OF
C                          THE SOLUTION AT THETA = A.
C
C                          WHEN MBDCND = 1, 2, OR 7,
C                            BDA(J) = U(A,PHI(J)) ,      J=1,2,...,N.
C
C                          WHEN MBDCND = 3, 4, OR 8,
C                            BDA(J) = (D/DTHETA)U(A,PHI(J)) ,
C                            J=1,2,...,N.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE,
C                          BDA IS A DUMMY VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT THETA = B.
C
C                          WHEN MBDCND = 1,4, OR 5,
C                            BDB(J) = U(B,PHI(J)) ,       J=1,2,...,N.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DTHETA)U(B,PHI(J)) ,
C                            J=1,2,...,N.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
C                          A DUMMY VARIABLE.
C
C                        C,D
C                          THE RANGE OF PHI (LONGITUDE),
C                          I.E. C .LE. PHI .LE. D.
C                          C MUST BE LESS THAN D.  IF D-C = 2*PI,
C                          PERIODIC BOUNDARY CONDITIONS ARE USUALLY
C                          USUALLY PRESCRIBED.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE PHI-DIRECTION
C                          ARE GIVEN BY PHI(J) = C + (J-0.5)DPHI,
C                          J=1,2,...,N, WHERE DPHI = (D-C)/N.
C                          N MUST BE GREATER THAN 2.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT PHI = C  AND PHI = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN PHI,
C                               I.E.  U(I,J) = U(I,N+J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               PHI = C AND PHI = D
C                               (SEE NOTE BELOW).
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               PHI = C AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO PHI IS
C                               SPECIFIED AT PHI = D
C                               (SEE NOTE BELOW).
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PHI = C AND PHI = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PHI = C AND THE SOLUTION IS
C                               SPECIFIED AT PHI = D
C                               (SEE NOTE BELOW).
C
C                          NOTE:
C                          WHEN NBDCND = 1, 2, OR 4, DO NOT USE
C                          MBDCND = 5, 6, 7, 8, OR 9
C                          (THE FORMER INDICATES THAT THE SOLUTION
C                          IS SPECIFIED AT A POLE; THE LATTER
C                          INDICATES THE SOLUTION IS UNSPECIFIED).
C                          USE INSTEAD MBDCND = 1 OR 2.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT PHI = C.
C
C                          WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(THETA(I),C) ,     I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DPHI)U(THETA(I),C),
C                            I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT PHI = D.
C
C                          WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(THETA(I),D) ,     I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DPHI)U(THETA(I),D) ,
C                            I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
C                          A SOLUTION MAY NOT EXIST.  HOWEVER,
C                          HSTSSP WILL ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE RIGHT SIDE OF THE
C                          HELMHOLTZ EQUATION.
C                          FOR I=1,2,...,M AND J=1,2,...,N
C
C                            F(I,J) = F(THETA(I),PHI(J)) .
C
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTSSP.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.
C                          IDIMF MUST BE AT LEAST M.
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (THETA(I),PHI(J)) FOR
C                          I=1,2,...,M, J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HSTSSP THEN COMPUTES
C                          THIS SOLUTION, WHICH IS A LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL APPROXIMATION.
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
C                          PARAMETERS. EXCEPT TO NUMBERS 0 AND 14,
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
C                          =  4  C .GE. D
C
C                          =  5  N .LE. 2
C
C                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4
C
C                          =  7  A .GT. 0 AND MBDCND = 5, 6, OR 9
C
C                          =  8  A = 0 AND MBDCND = 3, 4, OR 8
C
C                          =  9  B .LT. PI AND MBDCND .GE. 7
C
C                          = 10  B = PI AND MBDCND = 2,3, OR 6
C
C                          = 11  MBDCND .GE. 5 AND NDBCND = 1, 2, OR 4
C
C                          = 12  IDIMF .LT. M
C
C                          = 13  M .LE. 2
C
C                          = 14  LAMBDA .GT. 0
C
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HSTSSP, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED FILES         fish.f,comf.f,genbun.f,gnbnaux.f,poistg.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
C                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
C                        IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90.
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-
C                        DIFFERENCE EQUATIONS, INCORPORATES BOUNDARY
C                        DATA, ADJUSTS THE RIGHT SIDE WHEN THE SYSTEM
C                        IS SINGULAR AND CALLS EITHER POISTG OR GENBUN
C                        WHICH SOLVES THE LINEAR SYSTEM OF EQUATIONS.
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
C                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
C                        DIGITS FOR N AND M AS LARGE AS 64.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        ROUTINE POISTG WHICH IS THE ROUTINE THAT
C                        ACTUALLY SOLVES THE FINITE DIFFERENCE
C                        EQUATIONS.
C
C REFERENCES             U. SCHUMANN AND R. SWEET,"A DIRECT METHOD
C                        FOR THE SOLUTION OF POISSON'S EQUATION WITH
C                        NEUMANN BOUNDARY CONDITIONS ON A STAGGERED
C                        GRID OF ARBITRARY SIZE," J. COMP. PHYS.
C                        20(1976), PP. 171-182.
C***********************************************************************
      SUBROUTINE HSTSSP(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR)
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
      INTEGER :: IRWK, ICWK
      REAL :: PI, DUM
C-----------------------------------------------
C
      IERROR = 0
      PI = 4.0*ATAN(1.0)
      IF (A<0. .OR. B>PI) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<=0 .OR. MBDCND>9) IERROR = 3
      IF (C >= D) IERROR = 4
      IF (N <= 2) IERROR = 5
      IF (NBDCND<0 .OR. NBDCND>=5) IERROR = 6
      IF(A>0..AND.(MBDCND==5.OR.MBDCND==6.OR.MBDCND==9))IERROR=7
      IF(A==0..AND.(MBDCND==3.OR.MBDCND==4.OR.MBDCND==8))IERROR=8
      IF (B<PI .AND. MBDCND>=7) IERROR = 9
      IF(B==PI.AND.(MBDCND==2.OR.MBDCND==3.OR.MBDCND==6))IERROR=10
      IF (MBDCND>=5 .AND. (NBDCND==1 .OR. NBDCND==2 .OR. NBDCND==4)) 
     1   IERROR = 11
      IF (IDIMF < M) IERROR = 12
      IF (M <= 2) IERROR = 13
      IF (IERROR /= 0) RETURN 
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hstsspp(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +                   ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release dynamically allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HSTSSP


 
      SUBROUTINE HSTSSPP(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER  :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: IDIMF
      INTEGER , INTENT(OUT) :: IERROR
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
      REAL  :: F(IDIMF,1)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::NP,ISW,JSW,MB,IWB,IWC,IWR,IWS,I,J,MM1,K,LP,IERR1,I1
      REAL :: DELTAR, DLRSQ, DELTHT, DLTHSQ, PI, A1, A2, A3
C-----------------------------------------------
 
      DELTAR = (B - A)/FLOAT(M)
      DLRSQ = DELTAR**2
      DELTHT = (D - C)/FLOAT(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND + 1
      ISW = 1
      JSW = 1
      MB = MBDCND
      IF (ELMBDA == 0.) THEN
         GO TO (101,102,105,103,101,105,101,105,105) MBDCND
  101    CONTINUE
         IF (A/=0. .OR. B/=PI) GO TO 105
         MB = 9
         GO TO 104
  102    CONTINUE
         IF (A /= 0.) GO TO 105
         MB = 6
         GO TO 104
  103    CONTINUE
         IF (B /= PI) GO TO 105
         MB = 8
  104    CONTINUE
         JSW = 2
      ENDIF
  105 CONTINUE
      IWB = M
      IWC = IWB + M
      IWR = IWC + M
      IWS = IWR + M
      DO I = 1, M
         J = IWR + I
         W(J) = SIN(A + (FLOAT(I) - 0.5)*DELTAR)
         W(I) = SIN(A + FLOAT(I - 1)*DELTAR)/DLRSQ
      END DO
      MM1 = M - 1
      W(IWC+1:MM1+IWC) = W(2:MM1+1)
      W(IWB+1:MM1+IWB) = ELMBDA*W(IWR+1:MM1+IWR) - (W(:MM1)+W(2:MM1+1))
      W(IWR) = SIN(B)/DLRSQ
      W(IWC) = ELMBDA*W(IWS) - (W(M)+W(IWR))
      DO I = 1, M
         J = IWR + I
         A1 = W(J)
         F(I,:N) = A1*F(I,:N)
      END DO
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
      GO TO (110,110,112,112,114,114,110,112,114) MB
  110 CONTINUE
      A1 = 2.*W(1)
      W(IWB+1) = W(IWB+1) - W(1)
      F(1,:N) = F(1,:N) - A1*BDA(:N)
      GO TO 114
  112 CONTINUE
      A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1) + W(1)
      F(1,:N) = F(1,:N) + A1*BDA(:N)
  114 CONTINUE
      GO TO (115,117,117,115,115,117,119,119,119) MB
  115 CONTINUE
      A1 = 2.*W(IWR)
      W(IWC) = W(IWC) - W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
      GO TO 119
  117 CONTINUE
      A1 = DELTAR*W(IWR)
      W(IWC) = W(IWC) + W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
C
C     ENTER BOUNDARY DATA FOR PHI-BOUNDARIES.
C
  119 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (129,120,120,122,122) NP
  120 CONTINUE
      F(:M,1) = F(:M,1) - A1*BDC(:M)/W(IWR+1:M+IWR)
      GO TO 124
  122 CONTINUE
      A1 = 1./DELTHT
      F(:M,1) = F(:M,1) + A1*BDC(:M)/W(IWR+1:M+IWR)
  124 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (129,125,127,127,125) NP
  125 CONTINUE
      F(:M,N) = F(:M,N) - A1*BDD(:M)/W(IWR+1:M+IWR)
      GO TO 129
  127 CONTINUE
      A1 = 1./DELTHT
      F(:M,N) = F(:M,N) - A1*BDD(:M)/W(IWR+1:M+IWR)
  129 CONTINUE
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 14
         ELSE
            GO TO (139,139,132,139,139,132,139,132,132) MB
  132       CONTINUE
            GO TO (133,139,139,133,139) NP
  133       CONTINUE
            ISW = 2
            DO J = 1, N
               PERTRB = PERTRB + SUM(F(:M,J))
            END DO
            A1 = FLOAT(N)*(COS(A) - COS(B))/(2.*SIN(0.5*DELTAR))
            PERTRB = PERTRB/A1
            DO I = 1, M
               J = IWR + I
               A1 = PERTRB*W(J)
               F(I,:N) = F(I,:N) - A1
            END DO
            A2 = SUM(F(1,:N))
            A3 = SUM(F(M,:N))
            A2 = A2/W(IWR+1)
            A3 = A3/W(IWS)
         ENDIF
      ENDIF
  139 CONTINUE
      DO I = 1, M
         J = IWR + I
         A1 = DLTHSQ*W(J)
         W(I) = A1*W(I)
         J = IWC + I
         W(J) = A1*W(J)
         J = IWB + I
         W(J) = A1*W(J)
         F(I,:N) = A1*F(I,:N)
      END DO
      LP = NBDCND
      W(1) = 0.
      W(IWR) = 0.
C
C     CALL POISTG OR GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
C
      IERR1 = 0
      I1 = 1
      IF (NBDCND /= 0) THEN
         CALL POISTGG (LP, N, I1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
      ELSE
         CALL GENBUNN (LP, N, I1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
      ENDIF
      IF (ISW==2 .AND. JSW==2) THEN
         IF (MB == 8) THEN
            A1 = SUM(F(M,:N))
            A1 = (A1 - DLRSQ*A3/16.)/FLOAT(N)
            IF (NBDCND == 3) A1 = A1 + (BDD(M)-BDC(M))/(D - C)
            A1 = BDB(1) - A1
         ELSE
            A1 = SUM(F(1,:N))
            A1 = (A1 - DLRSQ*A2/16.)/FLOAT(N)
            IF (NBDCND == 3) A1 = A1 + (BDD(1)-BDC(1))/(D - C)
            A1 = BDA(1) - A1
         ENDIF
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
      END SUBROUTINE HSTSSPP
