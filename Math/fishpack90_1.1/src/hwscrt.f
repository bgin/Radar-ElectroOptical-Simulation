C
C     file hwscrt.f
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
C     SUBROUTINE HWSCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C DIMENSION OF           BDA(N),      BDB(N),   BDC(M),BDD(M),
C ARGUMENTS              F(IDIMF,N)
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
C                        EQUATION IN CARTESIAN COORDINATES.  THIS
C                        EQUATION IS
C
C                          (D/DX)(DU/DX) + (D/DY)(DU/DY)
C                          + LAMBDA*U = F(X,Y).
C
C USAGE                  CALL HWSCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               A,B
C
C                          THE RANGE OF X, I.E., A .LE. X .LE. B.
C                          A MUST BE LESS THAN B.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (A,B) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE M+1 GRID POINTS
C                          IN THE X-DIRECTION GIVEN BY
C                          X(I) = A+(I-1)DX FOR I = 1,2,...,M+1,
C                          WHERE DX = (B-A)/M IS THE PANEL WIDTH.
C                          M MUST BE GREATER THAN 3.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT X = A AND X = B.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN X,
C                               I.E., U(I,J) = U(M+I,J).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               X = A AND X = B.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               X = A AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO X IS
C                               SPECIFIED AT X = B.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED AT
C                               AT X = A AND X = B.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED AT
C                               X = A AND THE SOLUTION IS SPECIFIED
C                               AT X = B.
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO X AT X = A.
C
C                          WHEN MBDCND = 3 OR 4,
C
C                            BDA(J) = (D/DX)U(A,Y(J)), J = 1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
C                          A DUMMY VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
C                          THAT SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO X AT X = B.
C
C                          WHEN MBDCND = 2 OR 3,
C
C                            BDB(J) = (D/DX)U(B,Y(J)), J = 1,2,...,N+1
C
C                          WHEN MBDCND HAS ANY OTHER VALUE BDB IS A
C                          DUMMY VARIABLE.
C
C                        C,D
C                          THE RANGE OF Y, I.E., C .LE. Y .LE. D.
C                          C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (C,D) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE N+1 GRID POINTS IN THE
C                          Y-DIRECTION GIVEN BY Y(J) = C+(J-1)DY
C                          FOR J = 1,2,...,N+1, WHERE
C                          DY = (D-C)/N IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 3.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS AT
C                          Y = C AND Y = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Y,
C                               I.E., U(I,J) = U(I,N+J).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               Y = C AND Y = D.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               Y = C AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO Y IS
C                               SPECIFIED AT Y = D.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = C AND Y = D.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = C AND THE SOLUTION IS SPECIFIED
C                               AT Y = D.
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO Y AT Y = C.
C
C                          WHEN NBDCND = 3 OR 4,
C
C                            BDC(I) = (D/DY)U(X(I),C), I = 1,2,...,M+1
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
C                          A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO Y AT Y = D.
C
C                          WHEN NBDCND = 2 OR 3,
C
C                            BDD(I) = (D/DY)U(X(I),D), I = 1,2,...,M+1
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSCRT WILL
C                          ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
C                          LEAST (M+1)*(N+1), SPECIFYING VALUES OF THE
C                          RIGHT SIDE OF THE HELMHOLTZ  EQUATION AND
C                          BOUNDARY VALUES (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(X(I),Y(J)).
C
C                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
C                          FOR J=1,2,...,N+1,  I=1,2,...,M+1,
C
C                          MBDCND     F(1,J)        F(M+1,J)
C                          ------     ---------     --------
C
C                            0        F(A,Y(J))     F(A,Y(J))
C                            1        U(A,Y(J))     U(B,Y(J))
C                            2        U(A,Y(J))     F(B,Y(J))
C                            3        F(A,Y(J))     F(B,Y(J))
C                            4        F(A,Y(J))     U(B,Y(J))
C
C
C                          NBDCND     F(I,1)        F(I,N+1)
C                          ------     ---------     --------
C
C                            0        F(X(I),C)     F(X(I),C)
C                            1        U(X(I),C)     U(X(I),D)
C                            2        U(X(I),C)     F(X(I),D)
C                            3        F(X(I),C)     F(X(I),D)
C                            4        F(X(I),C)     U(X(I),D)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION U
C                          AND THE RIGHT SIDE F AT A CORNER THEN THE
C                          SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSCRT.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M+1  .
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (X(I),Y(J)), I = 1,2,...,M+1,
C                          J = 1,2,...,N+1  .
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
C                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
C                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
C                          CALCULATED AND SUBTRACTED FROM F, WHICH
C                          ENSURES THAT A SOLUTION EXISTS.  HWSCRT
C                          THEN COMPUTES THIS SOLUTION, WHICH IS A
C                          LEAST SQUARES SOLUTION TO THE ORIGINAL
C                          APPROXIMATION.  THIS SOLUTION PLUS ANY
C                          CONSTANT IS ALSO A SOLUTION.  HENCE, THE
C                          SOLUTION IS NOT UNIQUE.  THE VALUE OF
C                          PERTRB SHOULD BE SMALL COMPARED TO THE
C                          RIGHT SIDE F.  OTHERWISE, A SOLUTION IS
C                          OBTAINED TO AN ESSENTIALLY DIFFERENT
C                          PROBLEM. THIS COMPARISON SHOULD ALWAYS
C                          BE MADE TO INSURE THAT A MEANINGFUL
C                          SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 6,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          = 0  NO ERROR
C                          = 1  A .GE. B
C                          = 2  MBDCND .LT. 0 OR MBDCND .GT. 4
C                          = 3  C .GE. D
C                          = 4  N .LE. 3
C                          = 5  NBDCND .LT. 0 OR NBDCND .GT. 4
C                          = 6  LAMBDA .GT. 0
C                          = 7  IDIMF .LT. M+1
C                          = 8  M .LE. 3
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HWSCRT, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
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
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, "EFFICIENT
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
C                        ELLIPTIC EQUATIONS"
C                          NCAR TN/IA-109, JULY, 1975, 138 PP.
C***********************************************************************
      SUBROUTINE HWSCRT(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
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
      IF (N <= 3) IERROR = 4
      IF (NBDCND<0 .OR. NBDCND>4) IERROR = 5
      IF (IDIMF < M + 1) IERROR = 7
      IF (M <= 3) IERROR = 8
      IF (IERROR /= 0) RETURN 
!     allocate required work space length (generous estimate)
      IRWK=4*(N+1)+(13+INT(ALOG(FLOAT(N+1)/ALOG(2.0))))*(M+1)
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hwscrtt(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +             ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HWSCRT


 
      SUBROUTINE HWSCRTT(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
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
      INTEGER :: NPEROD, MPEROD, NP, NP1, MP, MP1, NSTART, NSTOP, NSKIP
     1   , NUNK, MSTART, MSTOP, MSKIP, J, MUNK, I, ID2, ID3, ID4, MSP1, 
     2   MSTM1, NSP1, NSTM1, IERR1
      REAL::DELTAX,TWDELX,DELXSQ,DELTAY,TWDELY,DELYSQ,S,ST2,A1,A2,S1
C-----------------------------------------------
      NPEROD = NBDCND
      MPEROD = 0
      IF (MBDCND > 0) MPEROD = 1
      DELTAX = (B - A)/FLOAT(M)
      TWDELX = 2./DELTAX
      DELXSQ = 1./DELTAX**2
      DELTAY = (D - C)/FLOAT(N)
      TWDELY = 2./DELTAY
      DELYSQ = 1./DELTAY**2
      NP = NBDCND + 1
      NP1 = N + 1
      MP = MBDCND + 1
      MP1 = M + 1
      NSTART = 1
      NSTOP = N
      NSKIP = 1
      GO TO (104,101,102,103,104) NP
  101 CONTINUE
      NSTART = 2
      GO TO 104
  102 CONTINUE
      NSTART = 2
  103 CONTINUE
      NSTOP = NP1
      NSKIP = 2
  104 CONTINUE
      NUNK = NSTOP - NSTART + 1
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      MSTART = 1
      MSTOP = M
      MSKIP = 1
      GO TO (117,105,106,109,110) MP
  105 CONTINUE
      MSTART = 2
      GO TO 107
  106 CONTINUE
      MSTART = 2
      MSTOP = MP1
      MSKIP = 2
  107 CONTINUE
      F(2,NSTART:NSTOP) = F(2,NSTART:NSTOP) - F(1,NSTART:NSTOP)*DELXSQ
      GO TO 112
  109 CONTINUE
      MSTOP = MP1
      MSKIP = 2
  110 CONTINUE
      F(1,NSTART:NSTOP) = F(1,NSTART:NSTOP) + BDA(NSTART:NSTOP)*TWDELX
  112 CONTINUE
      SELECT CASE (MSKIP) 
      CASE DEFAULT
         F(M,NSTART:NSTOP) = F(M,NSTART:NSTOP) - F(MP1,NSTART:NSTOP)*
     1      DELXSQ
      CASE (2) 
         F(MP1,NSTART:NSTOP) = F(MP1,NSTART:NSTOP) - BDB(NSTART:NSTOP)*
     1      TWDELX
      END SELECT
  117 CONTINUE
      MUNK = MSTOP - MSTART + 1
C
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
C
      GO TO (127,118,118,120,120) NP
  118 CONTINUE
      F(MSTART:MSTOP,2) = F(MSTART:MSTOP,2) - F(MSTART:MSTOP,1)*DELYSQ
      GO TO 122
  120 CONTINUE
      F(MSTART:MSTOP,1) = F(MSTART:MSTOP,1) + BDC(MSTART:MSTOP)*TWDELY
  122 CONTINUE
      SELECT CASE (NSKIP) 
      CASE DEFAULT
         F(MSTART:MSTOP,N) = F(MSTART:MSTOP,N) - F(MSTART:MSTOP,NP1)*
     1      DELYSQ
      CASE (2) 
         F(MSTART:MSTOP,NP1) = F(MSTART:MSTOP,NP1) - BDD(MSTART:MSTOP)*
     1      TWDELY
      END SELECT
C
C    MULTIPLY RIGHT SIDE BY DELTAY**2.
C
  127 CONTINUE
      DELYSQ = DELTAY*DELTAY
      F(MSTART:MSTOP,NSTART:NSTOP) = F(MSTART:MSTOP,NSTART:NSTOP)*DELYSQ
C
C     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = MUNK
      ID3 = ID2 + MUNK
      ID4 = ID3 + MUNK
      S = DELYSQ*DELXSQ
      ST2 = 2.*S
      W(:MUNK) = S
      W(ID2+1:MUNK+ID2) = (-ST2) + ELMBDA*DELYSQ
      W(ID3+1:MUNK+ID3) = S
      IF (MP /= 1) THEN
         W(1) = 0.
         W(ID4) = 0.
      ENDIF
      GO TO (135,135,132,133,134) MP
  132 CONTINUE
      W(ID2) = ST2
      GO TO 135
  133 CONTINUE
      W(ID2) = ST2
  134 CONTINUE
      W(ID3+1) = ST2
  135 CONTINUE
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 6
         ELSE
            IF ((NBDCND==0 .OR. NBDCND==3) .AND. (MBDCND==0 .OR. MBDCND
     1         ==3)) THEN
C
C     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
C     WILL EXIST.
C
               A1 = 1.
               A2 = 1.
               IF (NBDCND == 3) A2 = 2.
               IF (MBDCND == 3) A1 = 2.
               S1 = 0.
               MSP1 = MSTART + 1
               MSTM1 = MSTOP - 1
               NSP1 = NSTART + 1
               NSTM1 = NSTOP - 1
               DO J = NSP1, NSTM1
                  S = 0.
                  S = SUM(F(MSP1:MSTM1,J))
                  S1 = S1 + S*A1 + F(MSTART,J) + F(MSTOP,J)
               END DO
               S1 = A2*S1
               S = 0.
               S = SUM(F(MSP1:MSTM1,NSTART)+F(MSP1:MSTM1,NSTOP))
               S1 = S1 + S*A1 + F(MSTART,NSTART) + F(MSTART,NSTOP) + F(
     1            MSTOP,NSTART) + F(MSTOP,NSTOP)
               S = (2. + FLOAT(NUNK - 2)*A2)*(2. + FLOAT(MUNK - 2)*A1)
               PERTRB = S1/S
               F(MSTART:MSTOP,NSTART:NSTOP) = F(MSTART:MSTOP,NSTART:
     1            NSTOP) - PERTRB
               PERTRB = PERTRB/DELYSQ
C
C     SOLVE THE EQUATION.
C
            ENDIF
         ENDIF
      ENDIF
      IERR1 = 0
      CALL GENBUNN (NPEROD, NUNK, MPEROD, MUNK, W(1), W(ID2+1), W(ID3+1)
     1   , IDIMF, F(MSTART,NSTART), IERR1, W(ID4+1))
C
C     FILL IN IDENTICAL VALUES WHEN HAVE PERIODIC BOUNDARY CONDITIONS.
C
      IF (NBDCND == 0) THEN
         F(MSTART:MSTOP,NP1) = F(MSTART:MSTOP,1)
      ENDIF
      IF (MBDCND == 0) THEN
         F(MP1,NSTART:NSTOP) = F(1,NSTART:NSTOP)
         IF (NBDCND == 0) F(MP1,NP1) = F(1,NP1)
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
      END SUBROUTINE HWSCRTT
