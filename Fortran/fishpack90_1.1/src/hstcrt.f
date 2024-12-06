C
C     file hstcrt.f
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
C     SUBROUTINE HSTCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                 SOLVES THE STANDARD FIVE-POINT FINITE
C                         DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
C                         EQUATION
C                           (D/DX)(DU/DX) + (D/DY)(DU/DY) + LAMBDA*U
C                           = F(X,Y)
C                         ON A STAGGERED GRID IN CARTESIAN COORDINATES.
C
C USAGE                   CALL HSTCRT (A,B,M,MBDCND,BDA,BDB,C,D
C                                      N,NBDCND,BDC,BDD,ELMBDA,
C                                      F,IDIMF,PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT
C
C                        A,B
C                          THE RANGE OF X, I.E. A .LE. X .LE. B.
C                          A MUST BE LESS THAN B.
C
C                        M
C                          THE NUMBER OF GRID POINTS IN THE
C                          INTERVAL (A,B).  THE GRID POINTS
C                          IN THE X-DIRECTION ARE GIVEN BY
C                          X(I) = A + (I-0.5)DX FOR I=1,2,...,M
C                          WHERE DX =(B-A)/M.  M MUST BE GREATER
C                          THAN 2.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT X = A AND X = B.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN X,
C                               U(M+I,J) = U(I,J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               X = A AND X = B.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               X = A AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO X
C                               IS SPECIFIED AT X = B.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED
C                               AT X = A  AND X = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED
C                               AT X = A  AND THE SOLUTION IS
C                               SPECIFIED AT X = B.
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N
C                          THAT SPECIFIES THE BOUNDARY VALUES
C                          (IF ANY) OF THE SOLUTION AT X = A.
C
C                          WHEN MBDCND = 1 OR 2,
C                            BDA(J) = U(A,Y(J)) ,         J=1,2,...,N.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DX)U(A,Y(J)) ,   J=1,2,...,N.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N
C                          THAT SPECIFIES THE BOUNDARY VALUES
C                          OF THE SOLUTION AT X = B.
C
C                          WHEN MBDCND = 1 OR 4
C                            BDB(J) = U(B,Y(J)) ,        J=1,2,...,N.
C
C                          WHEN MBDCND = 2 OR 3
C                            BDB(J) = (D/DX)U(B,Y(J)) ,  J=1,2,...,N.
C
C                        C,D
C                          THE RANGE OF Y, I.E. C .LE. Y .LE. D.
C                          C MUST BE LESS THAN D.
C
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE Y-DIRECTION
C                          ARE GIVEN BY Y(J) = C + (J-0.5)DY,
C                          J=1,2,...,N, WHERE DY = (D-C)/N.
C                          N MUST BE GREATER THAN 2.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Y = C   AND Y = D.
C
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E.
C                               U(I,J) = U(I,N+J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT Y = C
C                               AND Y = D.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT Y = C
C                               AND THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = D.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = C AND Y = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = C AND THE SOLUTION IS SPECIFIED
C                               AT Y = D.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Y = C.
C
C                          WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(X(I),C) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DY)U(X(I),C),   I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Y = D.
C
C                          WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(X(I),D) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DY)U(X(I),D) ,  I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
C                          A SOLUTION MAY NOT EXIST. HOWEVER,
C                          HSTCRT WILL  ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE RIGHT SIDE OF THE
C                          HELMHOLTZ EQUATION.  FOR I=1,2,...,M
C                          AND J=1,2,...,N
C
C                            F(I,J) = F(X(I),Y(J)) .
C
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTCRT.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.
C                          IDIMF MUST BE AT LEAST M.
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (X(I),Y(J)) FOR  I=1,2,...,M, J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
C                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
C                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
C                          CALCULATED AND SUBTRACTED FROM F, WHICH
C                          ENSURES THAT A SOLUTION EXISTS.  HSTCRT
C                          THEN COMPUTES THIS SOLUTION, WHICH IS A
C                          LEAST SQUARES SOLUTION TO THE ORIGINAL
C                          APPROXIMATION.  THIS SOLUTION PLUS ANY
C                          CONSTANT IS ALSO A SOLUTION; HENCE, THE
C                          SOLUTION IS NOT UNIQUE.  THE VALUE OF
C                          PERTRB SHOULD BE SMALL COMPARED TO THE
C                          RIGHT SIDE F.  OTHERWISE, A SOLUTION IS
C                          OBTAINED TO AN ESSENTIALLY DIFFERENT PROBLEM.
C                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
C                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT TO NUMBERS 0 AND  6,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C
C                          =  1  A .GE. B
C
C                          =  2  MBDCND .LT. 0 OR MBDCND .GT. 4
C
C                          =  3  C .GE. D
C
C                          =  4  N .LE. 2
C
C                         =  5  NBDCND .LT. 0 OR NBDCND .GT. 4
C
C                         =  6  LAMBDA .GT. 0
C
C                         =  7  IDIMF .LT. M
C
C                         =  8  M .LE. 2
C
C                         SINCE THIS IS THE ONLY MEANS OF INDICATING
C                         A POSSIBLY INCORRECT CALL TO HSTCRT, THE
C                         USER SHOULD TEST IERROR AFTER THE CALL.
C
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C
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
C                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR
C                        AND CALLS EITHER POISTG OR GENBUN WHICH SOLVES
C                        THE LINEAR SYSTEM OF EQUATIONS.
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A
C                        LOSS OF NO MORE THAN FOUR SIGNIFICANT DIGITS
C                        FOR N AND M AS LARGE AS 64.  MORE DETAILED
C                        INFORMATION ABOUT ACCURACY CAN BE FOUND IN
C                        THE DOCUMENTATION FOR PACKAGE POISTG WHICH
C                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
C
C REFERENCES             U. SCHUMANN AND R. SWEET,"A DIRECT METHOD
C                        FOR THE SOLUTION OF POISSON'S EQUATION WITH
C                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
C                        ARBITRARY SIZE," J. COMP. PHYS. 20(1976),
C                        PP. 171-182.
C***********************************************************************
      SUBROUTINE HSTCRT(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
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
C     CHECK FOR INVALID PARAMETERS.
C
      IERROR = 0
      IF (A >= B) IERROR = 1
      IF (MBDCND<0 .OR. MBDCND>4) IERROR = 2
      IF (C >= D) IERROR = 3
      IF (N <= 2) IERROR = 4
      IF (NBDCND<0 .OR. NBDCND>4) IERROR = 5
      IF (IDIMF < M) IERROR = 7
      IF (M <= 2) IERROR = 8
      IF (IERROR /= 0) RETURN 
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hstcrtt(a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     +             elmbda,f,idimf,pertrb,ierror,w%rew)
      RETURN 
      END SUBROUTINE HSTCRT


 
      SUBROUTINE HSTCRTT(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
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
      INTEGER :: NPEROD, MPEROD, NP, MP, ID2, ID3, ID4, I, J, IERR1
      REAL::DELTAX,TWDELX,DELXSQ,DELTAY,TWDELY,DELYSQ,TWDYSQ,S,ST2
C-----------------------------------------------
 
      NPEROD = NBDCND
      MPEROD = 0
      IF (MBDCND > 0) MPEROD = 1
      DELTAX = (B - A)/FLOAT(M)
      TWDELX = 1./DELTAX
      DELXSQ = 2./DELTAX**2
      DELTAY = (D - C)/FLOAT(N)
      TWDELY = 1./DELTAY
      DELYSQ = DELTAY**2
      TWDYSQ = 2./DELYSQ
      NP = NBDCND + 1
      MP = MBDCND + 1
C
C     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = M
      ID3 = ID2 + M
      ID4 = ID3 + M
      S = (DELTAY/DELTAX)**2
      ST2 = 2.*S
      W(:M) = S
      W(ID2+1:M+ID2) = (-ST2) + ELMBDA*DELYSQ
      W(ID3+1:M+ID3) = S
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      GO TO (111,102,102,104,104) MP
  102 CONTINUE
      F(1,:N) = F(1,:N) - BDA(:N)*DELXSQ
      W(ID2+1) = W(ID2+1) - W(1)
      GO TO 106
  104 CONTINUE
      F(1,:N) = F(1,:N) + BDA(:N)*TWDELX
      W(ID2+1) = W(ID2+1) + W(1)
  106 CONTINUE
      GO TO (111,107,109,109,107) MP
  107 CONTINUE
      F(M,:N) = F(M,:N) - BDB(:N)*DELXSQ
      W(ID3) = W(ID3) - W(1)
      GO TO 111
  109 CONTINUE
      F(M,:N) = F(M,:N) - BDB(:N)*TWDELX
      W(ID3) = W(ID3) + W(1)
  111 CONTINUE
      GO TO (121,112,112,114,114) NP
  112 CONTINUE
      F(:M,1) = F(:M,1) - BDC(:M)*TWDYSQ
      GO TO 116
  114 CONTINUE
      F(:M,1) = F(:M,1) + BDC(:M)*TWDELY
  116 CONTINUE
      GO TO (121,117,119,119,117) NP
  117 CONTINUE
      F(:M,N) = F(:M,N) - BDD(:M)*TWDYSQ
      GO TO 121
  119 CONTINUE
      F(:M,N) = F(:M,N) - BDD(:M)*TWDELY
  121 CONTINUE
      F(:M,:N) = F(:M,:N)*DELYSQ
      IF (MPEROD /= 0) THEN
         W(1) = 0.
         W(ID4) = 0.
      ENDIF
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 6
         ELSE
            GO TO (127,133,133,127,133) MP
  127       CONTINUE
            GO TO (128,133,133,128,133) NP
C
C     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
C     WILL EXIST.
C
  128       CONTINUE
            S = 0.
            DO J = 1, N
               S = S + SUM(F(:M,J))
            END DO
            PERTRB = S/FLOAT(M*N)
            F(:M,:N) = F(:M,:N) - PERTRB
            PERTRB = PERTRB/DELYSQ
C
C     SOLVE THE EQUATION.
C
         ENDIF
      ENDIF
  133 CONTINUE
      IERR1 = 0
      IF (NPEROD /= 0) THEN
         CALL POISTGG (NPEROD, N, MPEROD, M, W(1), W(ID2+1), W(ID3+1), 
     1      IDIMF, F, IERR1, W(ID4+1))
      ELSE
         CALL GENBUNN (NPEROD, N, MPEROD, M, W(1), W(ID2+1), W(ID3+1), 
     1      IDIMF, F, IERR1, W(ID4+1))
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
      END SUBROUTINE HSTCRTT
