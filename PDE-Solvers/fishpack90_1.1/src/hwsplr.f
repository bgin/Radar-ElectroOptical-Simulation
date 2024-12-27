C
C     file hwsplr.f
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
C     SUBROUTINE HWSPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N+1)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION TO
C                        THE HELMHOLTZ EQUATION IN POLAR COORDINATES.
C                        THE EQUATION IS
C
C                            (1/R)(D/DR)(R(DU/DR)) +
C                            (1/R**2)(D/DTHETA)(DU/DTHETA) +
C                            LAMBDA*U = F(R,THETA).
C
C USAGE                  CALL HWSPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,
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
C                          FOR I = 1,2,...,M+1,
C                          WHERE DR = (B-A)/M IS THE PANEL WIDTH.
C                          M MUST BE GREATER THAN 3.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT R = A AND R = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               R = A AND R = B.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               R = A AND THE DERIVATIVE OF
C                               THE SOLUTION WITH RESPECT TO R IS
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
C                          NOTE:
C                          IF A = 0, DO NOT USE MBDCND = 3 OR 4, BUT
C                          INSTEAD USE MBDCND = 1,2,5, OR 6  .
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO R AT R = A.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DR)U(A,THETA(J)),
C                            J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
C                          A DUMMY VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO R AT R = B.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DR)U(B,THETA(J)),
C                            J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
C                          A DUMMY VARIABLE.
C
C                        C,D
C                          THE RANGE OF THETA, I.E., C .LE.
C                          THETA .LE. D.  C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (C,D) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE N+1 GRID POINTS IN THE
C                          THETA-DIRECTION GIVEN BY
C                          THETA(J) = C+(J-1)DTHETA FOR
C                          J = 1,2,...,N+1, WHERE
C                          DTHETA = (D-C)/N IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 3.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT THETA = C AND AT THETA = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN THETA,
C                               I.E., U(I,J) = U(I,N+J).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = C AND THETA = D
C                               (SEE NOTE BELOW).
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = C AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = D
C                               (SEE NOTE BELOW).
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = C AND THE SOLUTION IS
C                               SPECIFIED AT THETA = D
C                               (SEE NOTE BELOW).
C
C                          NOTE:
C                          WHEN NBDCND = 1,2, OR 4, DO NOT USE
C                          MBDCND = 5 OR 6
C                          (THE FORMER INDICATES THAT THE SOLUTION
C                          IS SPECIFIED AT R = 0, THE LATTER INDICATES
C                          THE SOLUTION IS UNSPECIFIED AT R = 0).
C                          USE INSTEAD MBDCND = 1 OR 2  .
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = C.  WHEN NBDCND = 3 OR 4,
C
C                            BDC(I) = (D/DTHETA)U(R(I),C),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
C                          A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = D.  WHEN NBDCND = 2 OR 3,
C
C                            BDD(I) = (D/DTHETA)U(R(I),D),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .LT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSPLR WILL
C                          ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
C                          LEAST (M+1)*(N+1), SPECIFYING VALUES
C                          OF THE RIGHT SIDE OF THE HELMHOLTZ
C                          EQUATION AND BOUNDARY DATA (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(R(I),THETA(J)).
C
C                          ON THE BOUNDARIES F IS DEFINED AS FOLLOWS:
C                          FOR J = 1,2,...,N+1 AND I = 1,2,...,M+1
C
C                          MBDCND   F(1,J)            F(M+1,J)
C                          ------   -------------     -------------
C
C                            1      U(A,THETA(J))     U(B,THETA(J))
C                            2      U(A,THETA(J))     F(B,THETA(J))
C                            3      F(A,THETA(J))     F(B,THETA(J))
C                            4      F(A,THETA(J))     U(B,THETA(J))
C                            5      F(0,0)            U(B,THETA(J))
C                            6      F(0,0)            F(B,THETA(J))
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
C                          THEN THE SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSPLR.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M+1.
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (R(I),THETA(J)),
C                          I = 1,2,...,M+1, J = 1,2,...,N+1  .
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HWSPLR THEN COMPUTES
C                          THIS SOLUTION, WHICH IS A LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
C                          A SOLUTION.  HENCE, THE SOLUTION IS NOT
C                          UNIQUE.  PERTRB SHOULD BE SMALL COMPARED
C                          TO THE RIGHT SIDE. OTHERWISE, A SOLUTION
C                          IS OBTAINED TO AN ESSENTIALLY DIFFERENT
C                          PROBLEM.  THIS COMPARISON SHOULD ALWAYS
C                          BE MADE TO INSURE THAT A MEANINGFUL
C                          SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 11,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR.
C                          =  1  A .LT. 0  .
C                          =  2  A .GE. B.
C                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6  .
C                          =  4  C .GE. D.
C                          =  5  N .LE. 3
C                          =  6  NBDCND .LT. 0 OR .GT. 4  .
C                          =  7  A = 0, MBDCND = 3 OR 4  .
C                          =  8  A .GT. 0, MBDCND .GE. 5  .
C                          =  9  MBDCND .GE. 5, NBDCND .NE. 0
C                                AND NBDCND .NE. 3  .
C                          = 10  IDIMF .LT. M+1  .
C                          = 11  LAMBDA .GT. 0  .
C                          = 12  M .LE. 3
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HWSPLR, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
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
      SUBROUTINE HWSPLR(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR)
!     Insert first required declarative statements for dynamically
!     allocated work space
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
      IF (MBDCND>=5 .AND. NBDCND/=0 .AND. NBDCND/=3) IERROR = 9
      IF (IDIMF < M + 1) IERROR = 10
      IF (M <= 3) IERROR = 12
      IF (IERROR /= 0) RETURN 
!     compute and allocate required work space
      IRWK=4*(N+1)+(M+1)*(13+INT(ALOG(FLOAT(N+1))/ALOG(2.0)))
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hwsplrr(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +                   ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HWSPLR


 
      SUBROUTINE HWSPLRR(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
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
     1   , ID2, ID3, ID4, ID5, ID6, IJ, I, J, L, LP, K, I1, IERR1, IP
      REAL :: ALL, DELTAR, DLRBY2, DLRSQ, DELTHT, DLTHSQ, A1, R, S2, A2
     1   , S, S1, YPOLE

      SAVE ALL
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
      MSTOP = MP1
      GO TO (101,105,102,103,104,105) MBDCND
  101 CONTINUE
      MSTOP = M
      GO TO 105
  102 CONTINUE
      MSTART = 1
      GO TO 105
  103 CONTINUE
      MSTART = 1
  104 CONTINUE
      MSTOP = M
  105 CONTINUE
      MUNK = MSTOP - MSTART + 1
      NSTART = 1
      NSTOP = N
      GO TO (109,106,107,108,109) NP
  106 CONTINUE
      NSTART = 2
      GO TO 109
  107 CONTINUE
      NSTART = 2
  108 CONTINUE
      NSTOP = NP1
  109 CONTINUE
      NUNK = NSTOP - NSTART + 1
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = MUNK
      ID3 = ID2 + MUNK
      ID4 = ID3 + MUNK
      ID5 = ID4 + MUNK
      ID6 = ID5 + MUNK
      A1 = 2./DLRSQ
      IJ = 0
      IF (MBDCND==3 .OR. MBDCND==4) IJ = 1
      DO I = 1, MUNK
         R = A + FLOAT(I - IJ)*DELTAR
         J = ID5 + I
         W(J) = R
         J = ID6 + I
         W(J) = 1./R**2
         W(I) = (R - DLRBY2)/(R*DLRSQ)
         J = ID3 + I
         W(J) = (R + DLRBY2)/(R*DLRSQ)
         J = ID2 + I
         W(J) = (-A1) + ELMBDA
      END DO
      GO TO (114,111,112,113,114,111) MBDCND
  111 CONTINUE
      W(ID2) = A1
      GO TO 114
  112 CONTINUE
      W(ID2) = A1
  113 CONTINUE
      W(ID3+1) = A1
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
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
  124 CONTINUE
      A1 = 1./DLTHSQ
      L = ID5 - MSTART + 1
      LP = ID6 - MSTART + 1
      GO TO (134,125,125,127,127) NP
  125 CONTINUE
      F(MSTART:MSTOP,2) = F(MSTART:MSTOP,2) - A1*W(MSTART+LP:MSTOP+LP)*F
     1   (MSTART:MSTOP,1)
      GO TO 129
  127 CONTINUE
      A1 = 2./DELTHT
      F(MSTART:MSTOP,1) = F(MSTART:MSTOP,1) + A1*W(MSTART+LP:MSTOP+LP)*
     1   BDC(MSTART:MSTOP)
  129 CONTINUE
      A1 = 1./DLTHSQ
      GO TO (134,130,132,132,130) NP
  130 CONTINUE
      F(MSTART:MSTOP,N) = F(MSTART:MSTOP,N) - A1*W(MSTART+LP:MSTOP+LP)*F
     1   (MSTART:MSTOP,NP1)
      GO TO 134
  132 CONTINUE
      A1 = 2./DELTHT
      F(MSTART:MSTOP,NP1) = F(MSTART:MSTOP,NP1) - A1*W(MSTART+LP:MSTOP+
     1   LP)*BDD(MSTART:MSTOP)
  134 CONTINUE
      IF (MBDCND>=5 .AND. NBDCND==3) F(1,1) = F(1,1) - (BDD(2)-BDC(2))*
     1   4./(FLOAT(N)*DELTHT*DLRSQ)
C
C     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
C     SOLUTION.
C
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 11
         ELSE
            IF (NBDCND==0 .OR. NBDCND==3) THEN
               S2 = 0.
               GO TO (144,144,137,144,144,138) MBDCND
  137          CONTINUE
               W(ID5+1) = 0.5*(W(ID5+2)-DLRBY2)
               S2 = 0.25*DELTAR
  138          CONTINUE
               A2 = 2.
               IF (NBDCND == 0) A2 = 1.
               J = ID5 + MUNK
               W(J) = 0.5*(W(J-1)+DLRBY2)
               S = 0.
               DO I = MSTART, MSTOP
                  S1 = 0.
                  IJ = NSTART + 1
                  K = NSTOP - 1
                  S1 = SUM(F(I,IJ:K))
                  J = I + L
                  S = S + (A2*S1 + F(I,NSTART)+F(I,NSTOP))*W(J)
               END DO
               S2=FLOAT(M)*A+DELTAR*(FLOAT((M-1)*(M+1))*0.5+0.25)+S2
               S1 = (2. + A2*FLOAT(NUNK - 2))*S2
               IF (MBDCND /= 3) THEN
                  S2 = FLOAT(N)*A2*DELTAR/8.
                  S = S + F(1,1)*S2
                  S1 = S1 + S2
               ENDIF
               PERTRB = S/S1
               F(MSTART:MSTOP,NSTART:NSTOP) = F(MSTART:MSTOP,NSTART:
     1            NSTOP) - PERTRB
            ENDIF
         ENDIF
      ENDIF
  144 CONTINUE
      DO I = MSTART, MSTOP
         K = I - MSTART + 1
         J = I + LP
         A1 = DLTHSQ/W(J)
         W(K) = A1*W(K)
         J = ID2 + K
         W(J) = A1*W(J)
         J = ID3 + K
         W(J) = A1*W(J)
         F(I,NSTART:NSTOP) = A1*F(I,NSTART:NSTOP)
      END DO
      W(1) = 0.
      W(ID4) = 0.
C
C     SOLVE THE SYSTEM OF EQUATIONS.
C
      I1 = 1
      IERR1 = 0
      CALL GENBUNN (NBDCND, NUNK, I1, MUNK, W(1), W(ID2+1), W(ID3+1), 
     1   IDIMF, F(MSTART,NSTART), IERR1, W(ID4+1))
      GO TO (157,157,157,157,148,147) MBDCND
C
C     ADJUST THE SOLUTION AS NECESSARY FOR THE PROBLEMS WHERE A = 0.
C
  147 CONTINUE
      IF (ELMBDA /= 0.) GO TO 148
      YPOLE = 0.
      GO TO 155
  148 CONTINUE
      J = ID5 + MUNK
      W(J) = W(ID2)/W(ID3)
      DO IP = 3, MUNK
         I = MUNK - IP + 2
         J = ID5 + I
         LP = ID2 + I
         K = ID3 + I
         W(J) = W(I)/(W(LP)-W(K)*W(J+1))
      END DO
      W(ID5+1) = -0.5*DLTHSQ/(W(ID2+1)-W(ID3+1)*W(ID5+2))
      DO I = 2, MUNK
         J = ID5 + I
         W(J) = -W(J)*W(J-1)
      END DO
      S = 0.
      S = SUM(F(2,NSTART:NSTOP))
      A2 = NUNK
      IF (NBDCND /= 0) THEN
         S = S - 0.5*(F(2,NSTART)+F(2,NSTOP))
         A2 = A2 - 1.
      ENDIF
      YPOLE = (0.25*DLRSQ*F(1,1)-S/A2)/(W(ID5+1)-1.+ELMBDA*DLRSQ*0.25)
      DO I = MSTART, MSTOP
         K = L + I
         F(I,NSTART:NSTOP) = F(I,NSTART:NSTOP) + YPOLE*W(K)
      END DO
  155 CONTINUE
      F(1,:NP1) = YPOLE
  157 CONTINUE
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
      END SUBROUTINE HWSPLRR
