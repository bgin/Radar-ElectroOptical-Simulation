C
C     file hstcyl.f
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
C     SUBROUTINE HSTCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION ON A STAGGERED
C                        GRID TO THE MODIFIED HELMHOLTZ EQUATION
C                        IN CYLINDRICAL COORDINATES. THIS EQUATION
C
C                          (1/R)(D/DR)(R(DU/DR)) + (D/DZ)(DU/DZ)
C
C                          + LAMBDA*(1/R**2)*U = F(R,Z)
C
C                        IS A TWO-DIMENSIONAL MODIFIED HELMHOLTZ
C                        EQUATION RESULTING FROM THE FOURIER TRANSFORM
C                        OF A THREE-DIMENSIONAL POISSON EQUATION.
C
C USAGE                  CALL HSTCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               A,B
C
C                          THE RANGE OF R, I.E. A .LE. R .LE. B.
C                          A MUST BE LESS THAN B AND A MUST BE
C                          BE NON-NEGATIVE.
C
C                        M
C                          THE NUMBER OF GRID POINTS IN THE INTERVAL
C                          (A,B).  THE GRID POINTS IN THE R-DIRECTION
C                          R-DIRECTION ARE GIVEN BY
C                          R(I) = A + (I-0.5)DR FOR I=1,2,...,M
C                          WHERE DR =(B-A)/M.
C                          M MUST BE GREATER THAN 2.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT R = A AND R = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT R = A
C                               (SEE NOTE BELOW) AND R = B.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT R = A
C                               (SEE NOTE BELOW) AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO R IS
C                               SPECIFIED AT R = B.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND R = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND THE
C                               SOLUTION IS SPECIFIED AT R = B.
C
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE SOLUTION IS
C                               SPECIFIED AT R = B.
C
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO R IS SPECIFIED
C                               AT R = B.
C
C                          NOTE:
C                          IF A = 0, DO NOT USE MBDCND = 1,2,3, OR 4,
C                          BUT INSTEAD USE MBDCND = 5 OR 6.
C                          THE RESULTING APPROXIMATION GIVES THE ONLY
C                          MEANINGFUL BOUNDARY CONDITION,
C                          I.E. DU/DR = 0.
C                          (SEE D. GREENSPAN, 'INTRODUCTORY NUMERICAL
C                          ANALYSIS OF ELLIPTIC BOUNDARY VALUE
C                          PROBLEMS,' HARPER AND ROW, 1965, CHAPTER 5.)
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES (IF ANY)
C                          OF THE SOLUTION AT R = A.
C
C                          WHEN MBDCND = 1 OR 2,
C                            BDA(J) = U(A,Z(J)) ,       J=1,2,...,N.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DR)U(A,Z(J)) ,   J=1,2,...,N.
C
C                          WHEN MBDCND = 5 OR 6, BDA IS A DUMMY
C                          VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT R = B.
C
C                          WHEN MBDCND = 1,4,OR 5,
C                            BDB(J) = U(B,Z(J)) ,        J=1,2,...,N.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DR)U(B,Z(J)) ,   J=1,2,...,N.
C
C                        C,D
C                          THE RANGE OF Z, I.E. C .LE. Z .LE. D.
C                          C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE Z-DIRECTION
C                          ARE GIVEN BY Z(J) = C + (J-0.5)DZ,
C                          J=1,2,...,N, WHERE DZ = (D-C)/N.
C                          N MUST BE GREATER THAN 2.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Z = C  AND Z = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.
C                               U(I,J) = U(I,N+J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT Z = C
C                               AND Z = D.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT Z = C
C                               AND THE DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO Z IS SPECIFIED AT Z = D.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO Z IS SPECIFIED AT Z = C
C                               AND Z = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO Z IS SPECIFIED AT Z = C AND
C                               THE SOLUTION IS SPECIFIED AT Z = D.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Z = C.
C
C                          WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(R(I),C) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DZ)U(R(I),C),    I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Z = D.
C
C                          WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(R(I),D) ,       I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DZ)U(R(I),D) ,   I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE MODIFIED
C                          HELMHOLTZ EQUATION.  IF LAMBDA IS GREATER
C                          THAN 0, A SOLUTION MAY NOT EXIST.
C                          HOWEVER, HSTCYL WILL ATTEMPT TO FIND A
C                          SOLUTION.  LAMBDA MUST BE ZERO WHEN
C                          MBDCND = 5 OR 6.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE RIGHT SIDE OF THE
C                          MODIFIED HELMHOLTZ EQUATION.
C                          FOR I=1,2,...,M   AND J=1,2,...,N
C                            F(I,J) = F(R(I),Z(J)) .
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTCYL.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M.
C
C ON OUTPUT
C
C                        F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (R(I),Z(J)) FOR  I=1,2,...,M, J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HSTCYL THEN COMPUTES
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
C                          PARAMETERS. EXCEPT TO NUMBERS 0 AND 11,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C
C                          =  1  A .LT. 0
C
C                          =  2  A .GE. B
C
C                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6
C
C                          =  4  C .GE. D
C
C                          =  5  N .LE. 2
C
C                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4
C
C                          =  7  A = 0 AND MBDCND = 1,2,3, OR 4
C
C                          =  8  A .GT. 0 AND MBDCND .GE. 5
C
C                          =  9  M .LE. 2
C
C                          = 10  IDIMF .LT. M
C
C                          = 11  LAMBDA .GT. 0
C
C                          = 12  A=0, MBDCND .GE. 5, ELMBDA .NE. 0
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HSTCYL, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       fish.f,comf.f,genbun.f,gnbnaux.f,poistg.f
C FILES
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
C                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
C                        IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, ADJUSTS
C                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR AND
C                        CALLS EITHER POISTG OR GENBUN WHICH SOLVES THE
C                        LINEAR SYSTEM OF EQUATIONS.
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               THE SOLUTION PROCESS RESULTS IN A LOSS
C                        OF NO MORE THAN FOUR SIGNIFICANT DIGITS
C                        FOR N AND M AS LARGE AS 64.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        SUBROUTINE POISTG WHICH IS THE ROUTINE THAT
C                        ACTUALLY SOLVES THE FINITE DIFFERENCE
C                        EQUATIONS.
C
C REFERENCES             U. SCHUMANN AND R. SWEET, "A DIRECT METHOD FOR
C                        THE SOLUTION OF POISSON'S EQUATION WITH NEUMANN
C                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
C                        ARBITRARY SIZE," J. COMP. PHYS. 20(1976),
C                        PP. 171-182.
C***********************************************************************
      SUBROUTINE HSTCYL(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
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
C-----------------------------------------------
C
      IERROR = 0
      IF (A < 0.) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<=0 .OR. MBDCND>=7) IERROR = 3
      IF (C >= D) IERROR = 4
      IF (N <= 2) IERROR = 5
      IF (NBDCND<0 .OR. NBDCND>=5) IERROR = 6
      IF (A==0. .AND. MBDCND/=5 .AND. MBDCND/=6) IERROR = 7
      IF (A>0. .AND. MBDCND>=5) IERROR = 8
      IF (IDIMF < M) IERROR = 10
      IF (M <= 2) IERROR = 9
      IF (A==0. .AND. MBDCND>=5 .AND. ELMBDA/=0.) IERROR = 12
      IF (IERROR /= 0) RETURN 
!     allocate real work space
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hstcyll(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +             ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HSTCYL


 
      SUBROUTINE HSTCYLL(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
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
      REAL  :: F(IDIMF,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NP, IWB, IWC, IWR, I, J, K, LP, IERR1
      REAL :: DELTAR, DLRSQ, DELTHT, DLTHSQ, A1
C-----------------------------------------------
      DELTAR = (B - A)/FLOAT(M)
      DLRSQ = DELTAR**2
      DELTHT = (D - C)/FLOAT(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND + 1
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      IWB = M
      IWC = IWB + M
      IWR = IWC + M
      DO I = 1, M
         J = IWR + I
         W(J) = A + (FLOAT(I) - 0.5)*DELTAR
         W(I) = (A + FLOAT(I - 1)*DELTAR)/(DLRSQ*W(J))
         K = IWC + I
         W(K) = (A + FLOAT(I)*DELTAR)/(DLRSQ*W(J))
         K = IWB + I
         W(K) = ELMBDA/W(J)**2 - 2./DLRSQ
      END DO
C
C     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
C
      GO TO (102,102,104,104,106,106) MBDCND
  102 CONTINUE
      A1 = 2.*W(1)
      W(IWB+1) = W(IWB+1) - W(1)
      F(1,:N) = F(1,:N) - A1*BDA(:N)
      GO TO 106
  104 CONTINUE
      A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1) + W(1)
      F(1,:N) = F(1,:N) + A1*BDA(:N)
  106 CONTINUE
      GO TO (107,109,109,107,107,109) MBDCND
  107 CONTINUE
      W(IWC) = W(IWC) - W(IWR)
      A1 = 2.*W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
      GO TO 111
  109 CONTINUE
      W(IWC) = W(IWC) + W(IWR)
      A1 = DELTAR*W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
  111 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (121,112,112,114,114) NP
  112 CONTINUE
      F(:M,1) = F(:M,1) - A1*BDC(:M)
      GO TO 116
  114 CONTINUE
      A1 = 1./DELTHT
      F(:M,1) = F(:M,1) + A1*BDC(:M)
  116 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (121,117,119,119,117) NP
  117 CONTINUE
      F(:M,N) = F(:M,N) - A1*BDD(:M)
      GO TO 121
  119 CONTINUE
      A1 = 1./DELTHT
      F(:M,N) = F(:M,N) - A1*BDD(:M)
  121 CONTINUE
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 11
         ELSE
            GO TO (130,130,124,130,130,124) MBDCND
  124       CONTINUE
            GO TO (125,130,130,125,130) NP
  125       CONTINUE
            DO I = 1, M
               A1 = 0.
               A1 = SUM(F(I,:N))
               J = IWR + I
               PERTRB = PERTRB + A1*W(J)
            END DO
            PERTRB = PERTRB/(FLOAT(M*N)*0.5*(A + B))
            F(:M,:N) = F(:M,:N) - PERTRB
         ENDIF
      ENDIF
  130 CONTINUE
      W(:M) = W(:M)*DLTHSQ
      W(IWC+1:M+IWC) = W(IWC+1:M+IWC)*DLTHSQ
      W(IWB+1:M+IWB) = W(IWB+1:M+IWB)*DLTHSQ
      F(:M,:N) = F(:M,:N)*DLTHSQ
      LP = NBDCND
      W(1) = 0.
      W(IWR) = 0.
C
C     SOLVE THE SYSTEM OF EQUATIONS.
C
      IERR1 = 0
      IF (NBDCND /= 0) THEN
         CALL POISTGG (LP, N, 1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
      ELSE
         CALL GENBUNN (LP, N, 1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
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
      END SUBROUTINE HSTCYLL
