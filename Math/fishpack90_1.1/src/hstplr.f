C
C     file hstplr.f
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
C     SUBROUTINE HSTPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
C    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
C
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION ON A STAGGERED
C                        GRID TO THE HELMHOLTZ EQUATION IN POLAR
C                        COORDINATES.  THE EQUATION IS
C
C                           (1/R)(D/DR)(R(DU/DR)) +
C                           (1/R**2)(D/DTHETA)(DU/DTHETA) +
C                           LAMBDA*U = F(R,THETA)
C
C USAGE                  CALL HSTPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,
C                                     IDIMF,PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               A,B
C
C                          THE RANGE OF R, I.E. A .LE. R .LE. B.
C                          A MUST BE LESS THAN B AND A MUST BE
C                          NON-NEGATIVE.
C
C                        M
C                          THE NUMBER OF GRID POINTS IN THE INTERVAL
C                          (A,B).  THE GRID POINTS IN THE R-DIRECTION
C                          ARE GIVEN BY R(I) = A + (I-0.5)DR FOR
C                          I=1,2,...,M WHERE DR =(B-A)/M.
C                          M MUST BE GREATER THAN 2.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT R = A AND R = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT R = A
C                               AND R = B.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT R = A
C                               AND THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT R = B.
C                               (SEE NOTE 1 BELOW)
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE 2 BELOW) AND R = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               SPECIFIED AT R = A (SEE NOTE 2 BELOW)
C                               AND THE SOLUTION IS SPECIFIED AT R = B.
C
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
C                          NOTE 1:
C                          IF A = 0, MBDCND = 2, AND NBDCND = 0 OR 3,
C                          THE SYSTEM OF EQUATIONS TO BE SOLVED IS
C                          SINGULAR.  THE UNIQUE SOLUTION IS
C                          IS DETERMINED BY EXTRAPOLATION TO THE
C                          SPECIFICATION OF U(0,THETA(1)).
C                          BUT IN THIS CASE THE RIGHT SIDE OF THE
C                          SYSTEM WILL BE PERTURBED BY THE CONSTANT
C                          PERTRB.
C
C                          NOTE 2:
C                          IF A = 0, DO NOT USE MBDCND = 3 OR 4,
C                          BUT INSTEAD USE MBDCND = 1,2,5, OR 6.
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES (IF ANY) OF
C                          THE SOLUTION AT R = A.
C
C                          WHEN MBDCND = 1 OR 2,
C                            BDA(J) = U(A,THETA(J)) ,     J=1,2,...,N.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DR)U(A,THETA(J)) ,
C                            J=1,2,...,N.
C
C                          WHEN MBDCND = 5 OR 6, BDA IS A DUMMY
C                          VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT R = B.
C
C                          WHEN MBDCND = 1,4, OR 5,
C                            BDB(J) = U(B,THETA(J)) ,     J=1,2,...,N.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DR)U(B,THETA(J)) ,
C                            J=1,2,...,N.
C
C                        C,D
C                          THE RANGE OF THETA, I.E. C .LE. THETA .LE. D.
C                          C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE THETA-
C                          DIRECTION ARE GIVEN BY THETA(J) = C +
C                          (J-0.5)DT,   J=1,2,...,N, WHERE
C                          DT = (D-C)/N.  N MUST BE GREATER THAN 2.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT THETA = C  AND THETA = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN THETA,
C                               I.E. U(I,J) = U(I,N+J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = C AND THETA = D
C                               (SEE NOTE BELOW).
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = C AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = D
C                               (SEE NOTE BELOW).
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = C AND THETA = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = C AND THE SOLUTION IS
C                               SPECIFIED AT THETA = D
C                               (SEE NOTE BELOW).
C
C                          NOTE:
C                          WHEN NBDCND = 1, 2, OR 4, DO NOT USE
C                          MBDCND = 5 OR 6 (THE FORMER INDICATES THAT
C                          THE SOLUTION IS SPECIFIED AT R =  0; THE
C                          LATTER INDICATES THE SOLUTION IS UNSPECIFIED
C                          AT R = 0).  USE INSTEAD MBDCND = 1 OR 2.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT THETA = C.
C
C                          WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(R(I),C) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DTHETA)U(R(I),C),
C                            I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT THETA = D.
C
C                          WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(R(I),D) ,         I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) =(D/DTHETA)U(R(I),D), I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
C                          A SOLUTION MAY NOT EXIST.  HOWEVER, HSTPLR
C                          WILL ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT SIDE OF THE HELMHOLTZ
C                          EQUATION.
C
C                          FOR I=1,2,...,M AND J=1,2,...,N
C                            F(I,J) = F(R(I),THETA(J)) .
C
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTPLR.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.
C                          IDIMF MUST BE AT LEAST M.
C
C
C ON OUTPUT
C
C                        F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (R(I),THETA(J)) FOR I=1,2,...,M,
C                          J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HSTPLR THEN COMPUTES THIS
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
C                          =  7  A = 0 AND MBDCND = 3 OR 4
C
C                          =  8  A .GT. 0 AND MBDCND .GE. 5
C
C                          =  9  MBDCND .GE. 5 AND NBDCND .NE. 0 OR 3
C
C                          = 10  IDIMF .LT. M
C
C                          = 11  LAMBDA .GT. 0
C
C                          = 12  M .LE. 2
C
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HSTPLR, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
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
C PORTABILITY            FORTRAN 90
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
C REFERENCES             U. SCHUMANN AND R. SWEET, "A DIRECT METHOD
C                        FOR THE SOLUTION OF POISSON'S EQUATION WITH
C                        NEUMANN BOUNDARY CONDITIONS ON A STAGGERED
C                        GRID OF ARBITRARY SIZE," J. COMP. PHYS.
C                        20(1976), PP. 171-182.
C***********************************************************************
      SUBROUTINE HSTPLR(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
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
      IERROR = 0
      IF (A < 0.) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<=0 .OR. MBDCND>=7) IERROR = 3
      IF (C >= D) IERROR = 4
      IF (N <= 2) IERROR = 5
      IF (NBDCND<0 .OR. NBDCND>=5) IERROR = 6
      IF (A==0. .AND. (MBDCND==3 .OR. MBDCND==4)) IERROR = 7
      IF (A>0. .AND. MBDCND>=5) IERROR = 8
      IF (MBDCND>=5 .AND. NBDCND/=0 .AND. NBDCND/=3) IERROR = 9
      IF (IDIMF < M) IERROR = 10
      IF (M <= 2) IERROR = 12
      IF (IERROR /= 0) RETURN 
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hstplrr(A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     +                   ELMBDA,F,IDIMF,PERTRB,IERROR,w%rew)
      RETURN 
      END SUBROUTINE HSTPLR


      SUBROUTINE HSTPLRR(A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
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
      INTEGER :: NP, ISW, MB, IWB, IWC, IWR, I, J, K, LP, IERR1
      REAL :: DELTAR, DLRSQ, DELTHT, DLTHSQ, A1, A2
C-----------------------------------------------
      DELTAR = (B - A)/FLOAT(M)
      DLRSQ = DELTAR**2
      DELTHT = (D - C)/FLOAT(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND + 1
      ISW = 1
      MB = MBDCND
      IF (A==0. .AND. MBDCND==2) MB = 6
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      IWB = M
      IWC = IWB + M
      IWR = IWC + M
      DO I = 1, M
         J = IWR + I
         W(J) = A + (FLOAT(I) - 0.5)*DELTAR
         W(I) = (A + FLOAT(I - 1)*DELTAR)/DLRSQ
         K = IWC + I
         W(K) = (A + FLOAT(I)*DELTAR)/DLRSQ
         K = IWB + I
         W(K) = (ELMBDA - 2./DLRSQ)*W(J)
      END DO
      DO I = 1, M
         J = IWR + I
         A1 = W(J)
         F(I,:N) = A1*F(I,:N)
      END DO
C
C     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
C
      GO TO (104,104,106,106,108,108) MB
  104 CONTINUE
      A1 = 2.*W(1)
      W(IWB+1) = W(IWB+1) - W(1)
      F(1,:N) = F(1,:N) - A1*BDA(:N)
      GO TO 108
  106 CONTINUE
      A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1) + W(1)
      F(1,:N) = F(1,:N) + A1*BDA(:N)
  108 CONTINUE
      GO TO (109,111,111,109,109,111) MB
  109 CONTINUE
      A1 = 2.*W(IWR)
      W(IWC) = W(IWC) - W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
      GO TO 113
  111 CONTINUE
      A1 = DELTAR*W(IWR)
      W(IWC) = W(IWC) + W(IWR)
      F(M,:N) = F(M,:N) - A1*BDB(:N)
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
  113 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (123,114,114,116,116) NP
  114 CONTINUE
      F(:M,1) = F(:M,1) - A1*BDC(:M)/W(IWR+1:M+IWR)
      GO TO 118
  116 CONTINUE
      A1 = 1./DELTHT
      F(:M,1) = F(:M,1) + A1*BDC(:M)/W(IWR+1:M+IWR)
  118 CONTINUE
      A1 = 2./DLTHSQ
      GO TO (123,119,121,121,119) NP
  119 CONTINUE
      F(:M,N) = F(:M,N) - A1*BDD(:M)/W(IWR+1:M+IWR)
      GO TO 123
  121 CONTINUE
      A1 = 1./DELTHT
      F(:M,N) = F(:M,N) - A1*BDD(:M)/W(IWR+1:M+IWR)
  123 CONTINUE
      PERTRB = 0.
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 11
         ELSE
            GO TO (133,133,126,133,133,126) MB
  126       CONTINUE
            GO TO (127,133,133,127,133) NP
  127       CONTINUE
            ISW = 2
            DO J = 1, N
               PERTRB = PERTRB + SUM(F(:M,J))
            END DO
            PERTRB = PERTRB/(FLOAT(M*N)*0.5*(A + B))
            DO I = 1, M
               J = IWR + I
               A1 = PERTRB*W(J)
               F(I,:N) = F(I,:N) - A1
            END DO
            A2 = SUM(F(1,:N))
            A2 = A2/W(IWR+1)
         ENDIF
      ENDIF
  133 CONTINUE
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
C     TO SOLVE THE SYSTEM OF EQUATIONS.
C
      IERR1 = 0
      IF (LP /= 0) THEN
         CALL POISTGG (LP, N, 1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
      ELSE
         CALL GENBUNN (LP, N, 1, M, W, W(IWB+1), W(IWC+1), IDIMF, F, 
     1      IERR1, W(IWR+1))
      ENDIF
      IF (.NOT.(A/=0. .OR. MBDCND/=2 .OR. ISW/=2)) THEN
         A1 = SUM(F(1,:N))
         A1 = (A1 - DLRSQ*A2/16.)/FLOAT(N)
         IF (NBDCND == 3) A1 = A1 + (BDD(1)-BDC(1))/(D - C)
         A1 = BDA(1) - A1
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
      END SUBROUTINE HSTPLRR
