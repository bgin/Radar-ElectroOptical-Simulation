C
C     file hwscyl.f
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
C     SUBROUTINE HWSCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N+1)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION
C                        TO THE HELMHOLTZ EQUATION IN CYLINDRICAL
C                        COORDINATES.  THIS MODIFIED HELMHOLTZ EQUATION
C
C                          (1/R)(D/DR)(R(DU/DR)) + (D/DZ)(DU/DZ)
C
C                          + (LAMBDA/R**2)U = F(R,Z)
C
C                        RESULTS FROM THE FOURIER TRANSFORM OF THE
C                        THREE-DIMENSIONAL POISSON EQUATION.
C
C USAGE                  CALL HWSCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT               A,B
C                          THE RANGE OF R, I.E., A .LE. R .LE. B.
C                          A MUST BE LESS THAN B AND A MUST BE
C                          NON-NEGATIVE.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (A,B) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE M+1 GRID POINTS IN THE
C                          R-DIRECTION GIVEN BY R(I) = A+(I-1)DR,
C                          FOR I = 1,2,...,M+1, WHERE DR = (B-A)/M
C                          IS THE PANEL WIDTH. M MUST BE GREATER
C                          THAN 3.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT R = A AND R = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               R = A AND R = B.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               R = A AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO R IS
C                               SPECIFIED AT R = B.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND R = B.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND THE
C                               SOLUTION IS SPECIFIED AT R = B.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE SOLUTION IS
C                               SPECIFIED AT R = B.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO R IS SPECIFIED
C                               AT R = B.
C
C                          IF A = 0, DO NOT USE MBDCND = 3 OR 4,
C                          BUT INSTEAD USE MBDCND = 1,2,5, OR 6  .
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO R AT R = A.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DR)U(A,Z(J)), J = 1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
C                          A DUMMY VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO R AT R = B.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DR)U(B,Z(J)), J = 1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
C                          A DUMMY VARIABLE.
C
C                        C,D
C                          THE RANGE OF Z, I.E., C .LE. Z .LE. D.
C                          C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (C,D) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE N+1 GRID POINTS IN THE
C                          Z-DIRECTION GIVEN BY Z(J) = C+(J-1)DZ,
C                          FOR J = 1,2,...,N+1,
C                          WHERE DZ = (D-C)/N IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 3.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Z = C AND Z = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Z,
C                               I.E., U(I,1) = U(I,N+1).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               Z = C AND Z = D.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               Z = C AND THE DERIVATIVE OF
C                               THE SOLUTION WITH RESPECT TO Z IS
C                               SPECIFIED AT Z = D.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Z IS
C                               SPECIFIED AT Z = C AND Z = D.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Z IS SPECIFIED AT
C                               Z = C AND THE SOLUTION IS SPECIFIED
C                               AT Z = D.
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO Z AT Z = C.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DZ)U(R(I),C), I = 1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
C                          A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO Z AT Z = D.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DZ)U(R(I),D), I = 1,2,...,M+1
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSCYL WILL
C                          ATTEMPT TO FIND A SOLUTION.  LAMBDA MUST
C                          BE ZERO WHEN MBDCND = 5 OR 6  .
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
C                          LEAST (M+1)*(N+1), SPECIFYING VALUES
C                          OF THE RIGHT SIDE OF THE HELMHOLTZ
C                          EQUATION AND BOUNDARY DATA (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(R(I),Z(J)).
C
C                          ON THE BOUNDARIES F IS DEFINED AS FOLLOWS:
C                          FOR J = 1,2,...,N+1 AND I = 1,2,...,M+1
C
C                          MBDCND   F(1,J)            F(M+1,J)
C                          ------   ---------         ---------
C
C                            1      U(A,Z(J))         U(B,Z(J))
C                            2      U(A,Z(J))         F(B,Z(J))
C                            3      F(A,Z(J))         F(B,Z(J))
C                            4      F(A,Z(J))         U(B,Z(J))
C                            5      F(0,Z(J))         U(B,Z(J))
C                            6      F(0,Z(J))         F(B,Z(J))
C
C                          NBDCND   F(I,1)            F(I,N+1)
C                          ------   ---------         ---------
C
C                            0      F(R(I),C)         F(R(I),C)
C                            1      U(R(I),C)         U(R(I),D)
C                            2      U(R(I),C)         F(R(I),D)
C                            3      F(R(I),C)         F(R(I),D)
C                            4      F(R(I),C)         U(R(I),D)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
C                          U AND THE RIGHT SIDE F AT A CORNER THEN
C                          THE SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSCYL.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M+1  .
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (R(I),Z(J)), I =1,2,...,M+1, J =1,2,...,N+1.
C
C                        PERTRB
C                          IF ONE SPECIFIES A COMBINATION OF PERIODIC,
C                          DERIVATIVE, AND UNSPECIFIED BOUNDARY
C                          CONDITIONS FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HWSCYL THEN COMPUTES
C                          THIS SOLUTION, WHICH IS A LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
C                          A SOLUTION.  HENCE, THE SOLUTION IS NOT
C                          UNIQUE. THE VALUE OF PERTRB SHOULD BE
C                          SMALL COMPARED TO THE RIGHT SIDE F.
C                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
C                          ESSENTIALLY DIFFERENT PROBLEM.  THIS
C                          COMPARISON SHOULD ALWAYS BE MADE TO INSURE
C                          THAT A MEANINGFUL SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG WHICH INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 11,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR.
C                          =  1  A .LT. 0  .
C                          =  2  A .GE. B.
C                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6  .
C                          =  4  C .GE. D.
C                          =  5  N .LE. 3
C                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4  .
C                          =  7  A = 0, MBDCND = 3 OR 4  .
C                          =  8  A .GT. 0, MBDCND .GE. 5  .
C                          =  9  A = 0, LAMBDA .NE. 0, MBDCND .GE. 5  .
C                          = 10  IDIMF .LT. M+1  .
C                          = 11  LAMBDA .GT. 0  .
C                          = 12  M .LE. 3
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HWSCYL, THE
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
      SUBROUTINE HWSCYL(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
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
      IF (A < 0.) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<=0 .OR. MBDCND>=7) IERROR = 3
      IF (C >= D) IERROR = 4
      IF (N <= 3) IERROR = 5
      IF (NBDCND<=(-1) .OR. NBDCND>=5) IERROR = 6
      IF (A==0. .AND. (MBDCND==3 .OR. MBDCND==4)) IERROR = 7
      IF (A>0. .AND. MBDCND>=5) IERROR = 8
      IF (A==0. .AND. ELMBDA/=0. .AND. MBDCND>=5) IERROR = 9
      IF (IDIMF < M + 1) IERROR = 10
      IF (M <= 3) IERROR = 12
      IF (IERROR /= 0) RETURN 
!     allocate real work space
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hwscyll(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +             ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HWSCYL


 
      SUBROUTINE HWSCYLL(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER  :: NBDCND
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
      INTEGER :: MP1, NP1, NP, MSTART, MSTOP, MUNK, NSTART, NSTOP, NUNK
     1   , ID2, ID3, ID4, ID5, ID6, ISTART, IJ, I, J, K, L, NSP1, NSTM1
     2   , IERR1, I1
      REAL::DELTAR,DLRBY2,DLRSQ,DELTHT,DLTHSQ,A1,R,A2,S,S1,S2
C-----------------------------------------------
      MP1 = M + 1
      DELTAR = (B - A)/FLOAT(M)
      DLRBY2 = DELTAR/2.
      DLRSQ = DELTAR**2
      NP1 = N + 1
      DELTHT = (D - C)/FLOAT(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND + 1
C
C     DEFINE RANGE OF INDICES I AND J FOR UNKNOWNS U(I,J).
C
      MSTART = 2
      MSTOP = M
      GO TO (104,103,102,101,101,102) MBDCND
  101 CONTINUE
      MSTART = 1
      GO TO 104
  102 CONTINUE
      MSTART = 1
  103 CONTINUE
      MSTOP = MP1
  104 CONTINUE
      MUNK = MSTOP - MSTART + 1
      NSTART = 1
      NSTOP = N
      GO TO (108,105,106,107,108) NP
  105 CONTINUE
      NSTART = 2
      GO TO 108
  106 CONTINUE
      NSTART = 2
  107 CONTINUE
      NSTOP = NP1
  108 CONTINUE
      NUNK = NSTOP - NSTART + 1
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = MUNK
      ID3 = ID2 + MUNK
      ID4 = ID3 + MUNK
      ID5 = ID4 + MUNK
      ID6 = ID5 + MUNK
      ISTART = 1
      A1 = 2./DLRSQ
      IJ = 0
      IF (MBDCND==3 .OR. MBDCND==4) IJ = 1
      IF (MBDCND > 4) THEN
         W(1) = 0.
         W(ID2+1) = -2.*A1
         W(ID3+1) = 2.*A1
         ISTART = 2
         IJ = 1
      ENDIF
      DO I = ISTART, MUNK
         R = A + FLOAT(I - IJ)*DELTAR
         J = ID5 + I
         W(J) = R
         J = ID6 + I
         W(J) = 1./R**2
         W(I) = (R - DLRBY2)/(R*DLRSQ)
         J = ID3 + I
         W(J) = (R + DLRBY2)/(R*DLRSQ)
         K = ID6 + I
         J = ID2 + I
         W(J) = (-A1) + ELMBDA*W(K)
      END DO
      GO TO (114,111,112,113,114,112) MBDCND
  111 CONTINUE
      W(ID2) = A1
      GO TO 114
  112 CONTINUE
      W(ID2) = A1
  113 CONTINUE
      W(ID3+1) = A1*FLOAT(ISTART)
  114 CONTINUE
      GO TO (115,115,117,117,119,119) MBDCND
  115 CONTINUE
      A1 = W(1)
      F(2,NSTART:NSTOP) = F(2,NSTART:NSTOP) - A1*F(1,NSTART:NSTOP)
      GO TO 119
  117 CONTINUE
      A1 = 2.*DELTAR*W(1)
      F(1,NSTART:NSTOP) = F(1,NSTART:NSTOP) + A1*BDA(NSTART:NSTOP)
  119 CONTINUE
      GO TO (120,122,122,120,120,122) MBDCND
  120 CONTINUE
      A1 = W(ID4)
      F(M,NSTART:NSTOP) = F(M,NSTART:NSTOP) - A1*F(MP1,NSTART:NSTOP)
      GO TO 124
  122 CONTINUE
      A1 = 2.*DELTAR*W(ID4)
      F(MP1,NSTART:NSTOP) = F(MP1,NSTART:NSTOP) - A1*BDB(NSTART:NSTOP)
C
C     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.
C
  124 CONTINUE
      A1 = 1./DLTHSQ
      L = ID5 - MSTART + 1
      GO TO (134,125,125,127,127) NP
  125 CONTINUE
      F(MSTART:MSTOP,2) = F(MSTART:MSTOP,2) - A1*F(MSTART:MSTOP,1)
      GO TO 129
  127 CONTINUE
      A1 = 2./DELTHT
      F(MSTART:MSTOP,1) = F(MSTART:MSTOP,1) + A1*BDC(MSTART:MSTOP)
  129 CONTINUE
      A1 = 1./DLTHSQ
      GO TO (134,130,132,132,130) NP
  130 CONTINUE
      F(MSTART:MSTOP,N) = F(MSTART:MSTOP,N) - A1*F(MSTART:MSTOP,NP1)
      GO TO 134
  132 CONTINUE
      A1 = 2./DELTHT
      F(MSTART:MSTOP,NP1) = F(MSTART:MSTOP,NP1) - A1*BDD(MSTART:MSTOP)
  134 CONTINUE
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 11
         ELSE
            W(ID5+1) = 0.5*(W(ID5+2)-DLRBY2)
            GO TO (146,146,138,146,146,137) MBDCND
  137       CONTINUE
            W(ID5+1) = 0.5*W(ID5+1)
  138       CONTINUE
            GO TO (140,146,146,139,146) NP
  139       CONTINUE
            A2 = 2.
            GO TO 141
  140       CONTINUE
            A2 = 1.
  141       CONTINUE
            K = ID5 + MUNK
            W(K) = 0.5*(W(K-1)+DLRBY2)
            S = 0.
            DO I = MSTART, MSTOP
               S1 = 0.
               NSP1 = NSTART + 1
               NSTM1 = NSTOP - 1
               S1 = SUM(F(I,NSP1:NSTM1))
               K = I + L
               S = S + (A2*S1 + F(I,NSTART)+F(I,NSTOP))*W(K)
            END DO
            S2 = FLOAT(M)*A + (0.75 + FLOAT((M - 1)*(M + 1)))*DLRBY2
            IF (MBDCND == 3) S2 = S2 + 0.25*DLRBY2
            S1 = (2. + A2*FLOAT(NUNK - 2))*S2
            PERTRB = S/S1
            F(MSTART:MSTOP,NSTART:NSTOP) = F(MSTART:MSTOP,NSTART:NSTOP)
     1          - PERTRB
         ENDIF
      ENDIF
  146 CONTINUE
      W(:MSTOP-MSTART+1) = W(:MSTOP-MSTART+1)*DLTHSQ
      W(ID2+1:MSTOP-MSTART+1+ID2) = W(ID2+1:MSTOP-MSTART+1+ID2)*DLTHSQ
      W(ID3+1:MSTOP-MSTART+1+ID3) = W(ID3+1:MSTOP-MSTART+1+ID3)*DLTHSQ
      F(MSTART:MSTOP,NSTART:NSTOP) = F(MSTART:MSTOP,NSTART:NSTOP)*DLTHSQ
      W(1) = 0.
      W(ID4) = 0.
C
C     SOLVE THE SYSTEM OF EQUATIONS.
C
      IERR1 = 0
      I1 = 1
      CALL GENBUNN (NBDCND, NUNK, I1, MUNK, W(1), W(ID2+1), W(ID3+1), 
     1   IDIMF, F(MSTART,NSTART), IERR1, W(ID4+1))
      IF (NBDCND == 0) THEN
         F(MSTART:MSTOP,NP1) = F(MSTART:MSTOP,1)
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
      END SUBROUTINE HWSCYLL
