C
C     file poistg.f
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
C     SUBROUTINE POISTG (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR)
C
C
C DIMENSION OF           A(M),  B(M),  C(M),  Y(IDIMY,N)
C ARGUMENTS
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE LINEAR SYSTEM OF EQUATIONS
C                        FOR UNKNOWN X VALUES, WHERE I=1,2,...,M
C                        AND J=1,2,...,N
C
C                        A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
C                        + X(I,J-1) - 2.*X(I,J) + X(I,J+1)
C                        = Y(I,J)
C
C                        THE INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
C                        I.E. X(0,J) = X(M,J) AND X(M+1,J) = X(1,J), AND
C                        X(I,0) MAY BE EQUAL TO X(I,1) OR -X(I,1), AND
C                        X(I,N+1) MAY BE EQUAL TO X(I,N) OR -X(I,N),
C                        DEPENDING ON AN INPUT PARAMETER.
C
C USAGE                  CALL POISTG (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,
C                                     IERROR)
C
C ARGUMENTS
C
C ON INPUT
C
C                        NPEROD
C                          INDICATES VALUES WHICH X(I,0) AND X(I,N+1)
C                          ARE ASSUMED TO HAVE.
C                          = 1 IF X(I,0) = -X(I,1) AND X(I,N+1) = -X(I,N
C                          = 2 IF X(I,0) = -X(I,1) AND X(I,N+1) =  X(I,N
C                          = 3 IF X(I,0) =  X(I,1) AND X(I,N+1) =  X(I,N
C                          = 4 IF X(I,0) =  X(I,1) AND X(I,N+1) = -X(I,N
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
C                          N MUST BE GREATER THAN 2.
C
C                        MPEROD
C                          = 0 IF A(1) AND C(M) ARE NOT ZERO
C                          = 1 IF A(1) = C(M) = 0
C
C                        M
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
C                          M MUST BE GREATER THAN 2.
C
C                        A,B,C
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
C                          SPECIFY THE COEFFICIENTS IN THE LINEAR
C                          EQUATIONS GIVEN ABOVE.  IF MPEROD = 0 THE
C                          ARRAY ELEMENTS MUST NOT DEPEND ON INDEX I,
C                          BUT MUST BE CONSTANT.  SPECIFICALLY, THE
C                          SUBROUTINE CHECKS THE FOLLOWING CONDITION
C                            A(I) = C(1)
C                            B(I) = B(1)
C                            C(I) = C(1)
C                          FOR I = 1, 2, ..., M.
C
C                        IDIMY
C                          THE ROW (OR FIRST) DIMENSION OF THE TWO-
C                          DIMENSIONAL ARRAY Y AS IT APPEARS IN THE
C                          PROGRAM CALLING POISTG.  THIS PARAMETER IS
C                          USED TO SPECIFY THE VARIABLE DIMENSION OF Y.
C                          IDIMY MUST BE AT LEAST M.
C
C                        Y
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
C                          OF EQUATIONS GIVEN ABOVE.
C                          Y MUST BE DIMENSIONED AT LEAST M X N.
C
C ON OUTPUT
C
C                        Y
C                          CONTAINS THE SOLUTION X.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBER ZERO, A
C                          SOLUTION IS NOT ATTEMPTED.
C                          = 0  NO ERROR
C                          = 1  IF M .LE. 2
C                          = 2  IF N .LE. 2
C                          = 3  IDIMY .LT. M
C                          = 4  IF NPEROD .LT. 1 OR NPEROD .GT. 4
C                          = 5  IF MPEROD .LT. 0 OR MPEROD .GT. 1
C                          = 6  IF MPEROD = 0 AND A(I) .NE. C(1)
C                               OR B(I) .NE. B(1) OR C(I) .NE. C(1)
C                               FOR SOME I = 1, 2, ..., M.
C                          = 7  IF MPEROD .EQ. 1 .AND.
C                               (A(1).NE.0 .OR. C(M).NE.0)
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING A
C                          POSSIBLY INCORRECT CALL TO POIS3D, THE USER
C                          SHOULD TEST IERROR AFTER THE CALL.
C
C
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED files         fish.f,gnbnaux.f,comf.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY, 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THIS SUBROUTINE IS AN IMPLEMENTATION OF THE
C                        ALGORITHM PRESENTED IN THE REFERENCE BELOW.
C
C TIMING                 FOR LARGE M AND N, THE EXECUTION TIME IS
C                        ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               TO MEASURE THE ACCURACY OF THE ALGORITHM A
C                        UNIFORM RANDOM NUMBER GENERATOR WAS USED TO
C                        CREATE A SOLUTION ARRAY X FOR THE SYSTEM GIVEN
C                        IN THE 'PURPOSE' SECTION ABOVE, WITH
C                          A(I) = C(I) = -0.5*B(I) = 1,    I=1,2,...,M
C                        AND, WHEN MPEROD = 1
C                          A(1) = C(M) = 0
C                          B(1) = B(M) =-1.
C
C                        THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN
C                        SYSTEM AND, USING DOUBLE PRECISION, A RIGHT SID
C                        Y WAS COMPUTED.  USING THIS ARRAY Y SUBROUTINE
C                        POISTG WAS CALLED TO PRODUCE AN APPROXIMATE
C                        SOLUTION Z.  THEN THE RELATIVE ERROR, DEFINED A
C                          E = MAX(ABS(Z(I,J)-X(I,J)))/MAX(ABS(X(I,J)))
C                        WHERE THE TWO MAXIMA ARE TAKEN OVER I=1,2,...,M
C                        AND J=1,2,...,N, WAS COMPUTED.  VALUES OF E ARE
C                        GIVEN IN THE TABLE BELOW FOR SOME TYPICAL VALUE
C                        OF M AND N.
C
C                        M (=N)    MPEROD    NPEROD      E
C                        ------    ------    ------    ------
C
C                          31        0-1       1-4     9.E-13
C                          31        1         1       4.E-13
C                          31        1         3       3.E-13
C                          32        0-1       1-4     3.E-12
C                          32        1         1       3.E-13
C                          32        1         3       1.E-13
C                          33        0-1       1-4     1.E-12
C                          33        1         1       4.E-13
C                          33        1         3       1.E-13
C                          63        0-1       1-4     3.E-12
C                          63        1         1       1.E-12
C                          63        1         3       2.E-13
C                          64        0-1       1-4     4.E-12
C                          64        1         1       1.E-12
C                          64        1         3       6.E-13
C                          65        0-1       1-4     2.E-13
C                          65        1         1       1.E-11
C                          65        1         3       4.E-13
C
C REFERENCES             SCHUMANN, U. AND R. SWEET,"A DIRECT METHOD
C                        FOR THE SOLUTION OF POISSON"S EQUATION WITH
C                        NEUMANN BOUNDARY CONDITIONS ON A STAGGERED
C                        GRID OF ARBITRARY SIZE," J. COMP. PHYS.
C                        20(1976), PP. 171-182.
C *********************************************************************
      SUBROUTINE POISTG(NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, IERROR)
      USE fish
      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: NPEROD
      INTEGER  :: N
      INTEGER  :: MPEROD
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER  :: IERROR
      REAL  :: A(*)
      REAL  :: B(*)
      REAL  :: C(*)
      REAL  :: Y(IDIMY,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, IRWK, ICWK
C-----------------------------------------------
      IERROR = 0
      IF (M <= 2) IERROR = 1
      IF (N <= 2) IERROR = 2
      IF (IDIMY < M) IERROR = 3
      IF (NPEROD<1 .OR. NPEROD>4) IERROR = 4
      IF (MPEROD<0 .OR. MPEROD>1) IERROR = 5
      IF (MPEROD /= 1) THEN
         DO I = 1, M
            IF (A(I) /= C(1)) GO TO 102
            IF (C(I) /= C(1)) GO TO 102
            IF (B(I) /= B(1)) GO TO 102
         END DO
         GO TO 104
  102    CONTINUE
         IERROR = 6
         RETURN 
      ENDIF
      IF (A(1)/=0. .OR. C(M)/=0.) IERROR = 7
  104 CONTINUE
      IF (IERROR /= 0) RETURN 
!     compute and allocate real work space for poistg
      IRWK = M*(9 + INT(ALOG(FLOAT(N))/ALOG(2.0))) + 4*N
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     return if allocation failed (e.g., if n,m are too large)
      IF (IERROR == 20) RETURN 
      call  poistgg(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,w%rew)
!     release work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE POISTG


 
      SUBROUTINE POISTGG(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER  :: N
      INTEGER , INTENT(IN) :: MPEROD
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER  :: IERROR
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL  :: Y(IDIMY,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IWBA, IWBB, IWBC, IWB2, IWB3, IWW1, IWW2, IWW3, IWD, 
     1   IWTCOS, IWP, I, K, J, NP, MP, IPSTOR, IREV, MH, MHM1, MODD, 
     2   MHPI, MHMI, NBY2, MSKIP
      REAL :: A1
C-----------------------------------------------
      IWBA = M + 1
      IWBB = IWBA + M
      IWBC = IWBB + M
      IWB2 = IWBC + M
      IWB3 = IWB2 + M
      IWW1 = IWB3 + M
      IWW2 = IWW1 + M
      IWW3 = IWW2 + M
      IWD = IWW3 + M
      IWTCOS = IWD + M
      IWP = IWTCOS + 4*N
      DO I = 1, M
         K = IWBA + I - 1
         W(K) = -A(I)
         K = IWBC + I - 1
         W(K) = -C(I)
         K = IWBB + I - 1
         W(K) = 2. - B(I)
         Y(I,:N) = -Y(I,:N)
      END DO
      NP = NPEROD
      MP = MPEROD + 1
      GO TO (110,107) MP
  107 CONTINUE
      GO TO (108,108,108,119) NPEROD
  108 CONTINUE
      CALL POSTG2 (NP, N, M, W(IWBA), W(IWBB), W(IWBC), IDIMY, Y, W, W(
     1   IWB2), W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), W
     2   (IWP))
      IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD == 4) GO TO 120
  109 CONTINUE
      GO TO (123,129) MP
  110 CONTINUE
      MH = (M + 1)/2
      MHM1 = MH - 1
      MODD = 1
      IF (MH*2 == M) MODD = 2
      DO J = 1, N
         DO I = 1, MHM1
            W(I) = Y(MH-I,J) - Y(I+MH,J)
            W(I+MH) = Y(MH-I,J) + Y(I+MH,J)
         END DO
         W(MH) = 2.*Y(MH,J)
         GO TO (113,112) MODD
  112    CONTINUE
         W(M) = 2.*Y(M,J)
  113    CONTINUE
         Y(:M,J) = W(:M)
      END DO
      K = IWBC + MHM1 - 1
      I = IWBA + MHM1
      W(K) = 0.
      W(I) = 0.
      W(K+1) = 2.*W(K+1)
      SELECT CASE (MODD) 
      CASE DEFAULT
         K = IWBB + MHM1 - 1
         W(K) = W(K) - W(I-1)
         W(IWBC-1) = W(IWBC-1) + W(IWBB-1)
      CASE (2) 
         W(IWBB-1) = W(K+1)
      END SELECT
  118 CONTINUE
      GO TO 107
  119 CONTINUE
      IREV = 1
      NBY2 = N/2
      NP = 2
  120 CONTINUE
      DO J = 1, NBY2
         MSKIP = N + 1 - J
         DO I = 1, M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
         END DO
      END DO
      GO TO (108,109) IREV
  123 CONTINUE
      DO J = 1, N
         W(MH-1:MH-MHM1:(-1)) = 0.5*(Y(MH+1:MHM1+MH,J)+Y(:MHM1,J))
         W(MH+1:MHM1+MH) = 0.5*(Y(MH+1:MHM1+MH,J)-Y(:MHM1,J))
         W(MH) = 0.5*Y(MH,J)
         GO TO (126,125) MODD
  125    CONTINUE
         W(M) = 0.5*Y(M,J)
  126    CONTINUE
         Y(:M,J) = W(:M)
      END DO
  129 CONTINUE
      W(1) = IPSTOR + IWP - 1
      RETURN 
      END SUBROUTINE POISTGG


      SUBROUTINE POSTG2(NPEROD, N, M, A, BB, C, IDIMQ, Q, B, B2, B3, W, 
     1   W2, W3, D, TCOS, P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: IDIMQ
      REAL  :: A(*)
      REAL  :: BB(*)
      REAL  :: C(*)
      REAL , INTENT(INOUT) :: Q(IDIMQ,*)
      REAL  :: B(*)
      REAL  :: B2(*)
      REAL  :: B3(*)
      REAL  :: W(*)
      REAL  :: W2(*)
      REAL  :: W3(*)
      REAL  :: D(*)
      REAL  :: TCOS(*)
      REAL , INTENT(INOUT) :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER , DIMENSION(4) :: K
      INTEGER :: K1, K2, K3, K4, NP, MR, IP, IPSTOR, I2R, JR, NR, NLAST
     1   , KR, LR, NROD, JSTART, JSTOP, I2RBY2, J, IJUMP, JP1, JP2, JP3
     2   , JM1, JM2, JM3, I, NRODPR, II, NLASTP, JSTEP
      REAL :: FNUM, FNUM2, FI, T
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION ON A STAGGERED GRID.
C
C
      EQUIVALENCE (K(1), K1), (K(2), K2), (K(3), K3), (K(4), K4)
      NP = NPEROD
      FNUM = 0.5*FLOAT(NP/3)
      FNUM2 = 0.5*FLOAT(NP/2)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      IF (NR > 3) THEN
  101    CONTINUE
         JR = 2*I2R
         NROD = 1
         IF ((NR/2)*2 == NR) NROD = 0
         JSTART = 1
         JSTOP = NLAST - JR
         IF (NROD == 0) JSTOP = JSTOP - I2R
         I2RBY2 = I2R/2
         IF (JSTOP < JSTART) THEN
            J = JR
         ELSE
            IJUMP = 1
            DO J = JSTART, JSTOP, JR
               JP1 = J + I2RBY2
               JP2 = J + I2R
               JP3 = JP2 + I2RBY2
               JM1 = J - I2RBY2
               JM2 = J - I2R
               JM3 = JM2 - I2RBY2
               IF (J == 1) THEN
                  CALL COSGEN (I2R, 1, FNUM, 0.5, TCOS)
                  IF (I2R == 1) THEN
                     B(:MR) = Q(:MR,1)
                     Q(:MR,1) = Q(:MR,2)
                     GO TO 112
                  ENDIF
                  B(:MR) = Q(:MR,1) + 0.5*(Q(:MR,JP2)-Q(:MR,JP1)-Q(:MR,
     1               JP3))
                  Q(:MR,1) = Q(:MR,JP2) + Q(:MR,1) - Q(:MR,JP1)
                  GO TO 112
               ENDIF
               GO TO (107,108) IJUMP
  107          CONTINUE
               IJUMP = 2
               CALL COSGEN (I2R, 1, 0.5, 0.0, TCOS)
  108          CONTINUE
               IF (I2R == 1) THEN
                  B(:MR) = 2.*Q(:MR,J)
                  Q(:MR,J) = Q(:MR,JM2) + Q(:MR,JP2)
               ELSE
                  DO I = 1, MR
                     FI = Q(I,J)
                     Q(I,J)=Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
                     B(I) = FI + Q(I,J) - Q(I,JM3) - Q(I,JP3)
                  END DO
               ENDIF
  112          CONTINUE
               CALL TRIX (I2R, 0, MR, A, BB, C, B, TCOS, D, W)
               Q(:MR,J) = Q(:MR,J) + B(:MR)
C
C     END OF REDUCTION FOR REGULAR UNKNOWNS.
C
            END DO
C
C     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
C
            J = JSTOP + JR
         ENDIF
         NLAST = J
         JM1 = J - I2RBY2
         JM2 = J - I2R
         JM3 = JM2 - I2RBY2
         IF (NROD /= 0) THEN
C
C     ODD NUMBER OF UNKNOWNS
C
            IF (I2R == 1) THEN
               B(:MR) = Q(:MR,J)
               Q(:MR,J) = Q(:MR,JM2)
            ELSE
               B(:MR)=Q(:MR,J)+0.5*(Q(:MR,JM2)-Q(:MR,JM1)-Q(:MR,JM3))
               IF (NRODPR == 0) THEN
                  Q(:MR,J) = Q(:MR,JM2) + P(IP+1:MR+IP)
                  IP = IP - MR
               ELSE
                  Q(:MR,J) = Q(:MR,J) - Q(:MR,JM1) + Q(:MR,JM2)
               ENDIF
               IF (LR /= 0) CALL COSGEN (LR, 1, FNUM2, 0.5, TCOS(KR+1))
            ENDIF
            CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS)
            CALL TRIX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,J) = Q(:MR,J) + B(:MR)
            KR = KR + I2R
         ELSE
            JP1 = J + I2RBY2
            JP2 = J + I2R
            IF (I2R == 1) THEN
               B(:MR) = Q(:MR,J)
               TCOS(1) = 0.
               CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
               IP = 0
               IPSTOR = MR
               P(:MR) = B(:MR)
               B(:MR) = B(:MR) + Q(:MR,N)
               TCOS(1) = (-1.) + 2.*FLOAT(NP/2)
               TCOS(2) = 0.
               CALL TRIX (1, 1, MR, A, BB, C, B, TCOS, D, W)
               Q(:MR,J) = Q(:MR,JM2) + P(:MR) + B(:MR)
            ELSE
               B(:MR)=Q(:MR,J)+0.5*(Q(:MR,JM2)-Q(:MR,JM1)-Q(:MR,JM3))
               IF (NRODPR == 0) THEN
                  B(:MR) = B(:MR) + P(IP+1:MR+IP)
               ELSE
                  B(:MR) = B(:MR) + Q(:MR,JP2) - Q(:MR,JP1)
               ENDIF
               CALL COSGEN (I2R, 1, 0.5, 0.0, TCOS)
               CALL TRIX (I2R, 0, MR, A, BB, C, B, TCOS, D, W)
               IP = IP + MR
               IPSTOR = MAX0(IPSTOR,IP + MR)
               P(IP+1:MR+IP) = B(:MR) + 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,
     1            JP1))
               B(:MR) = P(IP+1:MR+IP) + Q(:MR,JP2)
               IF (LR /= 0) THEN
                  CALL COSGEN (LR, 1, FNUM2, 0.5, TCOS(I2R+1))
                  CALL MERGE (TCOS, 0, I2R, I2R, LR, KR)
               ELSE
                  DO I = 1, I2R
                     II = KR + I
                     TCOS(II) = TCOS(I)
                  END DO
               ENDIF
               CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS)
               CALL TRIX (KR, KR, MR, A, BB, C, B, TCOS, D, W)
               Q(:MR,J) = Q(:MR,JM2) + P(IP+1:MR+IP) + B(:MR)
            ENDIF
            LR = KR
            KR = KR + JR
         ENDIF
         NR = (NLAST - 1)/JR + 1
         IF (NR <= 3) GO TO 142
         I2R = JR
         NRODPR = NROD
         GO TO 101
      ENDIF
  142 CONTINUE
      J = 1 + JR
      JM1 = J - I2R
      JP1 = J + I2R
      JM2 = NLAST - I2R
      IF (NR /= 2) THEN
         IF (LR == 0) THEN
            IF (N == 3) THEN
C
C     CASE N = 3.
C
               GO TO (143,148,143) NP
  143          CONTINUE
               B(:MR) = Q(:MR,2)
               B2(:MR) = Q(:MR,1) + Q(:MR,3)
               B3(:MR) = 0.
               SELECT CASE (NP) 
               CASE DEFAULT
                  TCOS(1) = -1.
                  TCOS(2) = 1.
                  K1 = 1
               CASE (1:2) 
                  TCOS(1) = -2.
                  TCOS(2) = 1.
                  TCOS(3) = -1.
                  K1 = 2
               END SELECT
  147          CONTINUE
               K2 = 1
               K3 = 0
               K4 = 0
               GO TO 150
  148          CONTINUE
               B(:MR) = Q(:MR,2)
               B2(:MR) = Q(:MR,3)
               B3(:MR) = Q(:MR,1)
               CALL COSGEN (3, 1, 0.5, 0.0, TCOS)
               TCOS(4) = -1.
               TCOS(5) = 1.
               TCOS(6) = -1.
               TCOS(7) = 1.
               K1 = 3
               K2 = 2
               K3 = 1
               K4 = 1
  150          CONTINUE
               CALL TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
               B(:MR) = B(:MR) + B2(:MR) + B3(:MR)
               GO TO (153,153,152) NP
  152          CONTINUE
               TCOS(1) = 2.
               CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
  153          CONTINUE
               Q(:MR,2) = B(:MR)
               B(:MR) = Q(:MR,1) + B(:MR)
               TCOS(1) = (-1.) + 4.*FNUM
               CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
               Q(:MR,1) = B(:MR)
               JR = 1
               I2R = 0
               GO TO 188
            ENDIF
C
C     CASE N = 2**P+1
C
            B(:MR)=Q(:MR,J)+Q(:MR,1)-Q(:MR,JM1)+Q(:MR,NLAST)-Q(:MR,JM2)
            GO TO (158,160,158) NP
  158       CONTINUE
            B2(:MR) = Q(:MR,1) + Q(:MR,NLAST) + Q(:MR,J) - Q(:MR,JM1) - 
     1         Q(:MR,JP1)
            B3(:MR) = 0.
            K1 = NLAST - 1
            K2 = NLAST + JR - 1
            CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(NLAST))
            TCOS(K2) = 2.*FLOAT(NP - 2)
            CALL COSGEN (JR, 1, 0.5 - FNUM, 0.5, TCOS(K2+1))
            K3 = (3 - NP)/2
            CALL MERGE (TCOS, K1, JR - K3, K2 - K3, JR + K3, 0)
            K1 = K1 - 1 + K3
            CALL COSGEN (JR, 1, FNUM, 0.5, TCOS(K1+1))
            K2 = JR
            K3 = 0
            K4 = 0
            GO TO 162
  160       CONTINUE
            DO I = 1, MR
               FI = (Q(I,J)-Q(I,JM1)-Q(I,JP1))/2.
               B2(I) = Q(I,1) + FI
               B3(I) = Q(I,NLAST) + FI
            END DO
            K1 = NLAST + JR - 1
            K2 = K1 + JR - 1
            CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(K1+1))
            CALL COSGEN (NLAST, 1, 0.5, 0.0, TCOS(K2+1))
            CALL MERGE (TCOS, K1, JR - 1, K2, NLAST, 0)
            K3 = K1 + NLAST - 1
            K4 = K3 + JR
            CALL COSGEN (JR, 1, 0.5, 0.5, TCOS(K3+1))
            CALL COSGEN (JR, 1, 0.0, 0.5, TCOS(K4+1))
            CALL MERGE (TCOS, K3, JR, K4, JR, K1)
            K2 = NLAST - 1
            K3 = JR
            K4 = JR
  162       CONTINUE
            CALL TRI3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
            B(:MR) = B(:MR) + B2(:MR) + B3(:MR)
            IF (NP == 3) THEN
               TCOS(1) = 2.
               CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            ENDIF
            Q(:MR,J) = B(:MR) + 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1))
            B(:MR) = Q(:MR,J) + Q(:MR,1)
            CALL COSGEN (JR, 1, FNUM, 0.5, TCOS)
            CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,1) = Q(:MR,1) - Q(:MR,JM1) + B(:MR)
            GO TO 188
         ENDIF
C
C     CASE OF GENERAL N WITH NR = 3 .
C
         B(:MR) = Q(:MR,1) - Q(:MR,JM1) + Q(:MR,J)
         IF (NROD == 0) THEN
            B(:MR) = B(:MR) + P(IP+1:MR+IP)
         ELSE
            B(:MR) = B(:MR) + Q(:MR,NLAST) - Q(:MR,JM2)
         ENDIF
         DO I = 1, MR
            T = 0.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
            Q(I,J) = T
            B2(I) = Q(I,NLAST) + T
            B3(I) = Q(I,1) + T
         END DO
         K1 = KR + 2*JR
         CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(K1+1))
         K2 = K1 + JR
         TCOS(K2) = 2.*FLOAT(NP - 2)
         K4 = (NP - 1)*(3 - NP)
         K3 = K2 + 1 - K4
         CALL COSGEN(KR+JR+K4,1,FLOAT(K4)/2.,1.-FLOAT(K4),TCOS(K3))
         K4 = 1 - NP/3
         CALL MERGE (TCOS, K1, JR - K4, K2 - K4, KR + JR + K4, 0)
         IF (NP == 3) K1 = K1 - 1
         K2 = KR + JR
         K4 = K1 + K2
         CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS(K4+1))
         K3 = K4 + KR
         CALL COSGEN (JR, 1, FNUM, 0.5, TCOS(K3+1))
         CALL MERGE (TCOS, K4, KR, K3, JR, K1)
         K4 = K3 + JR
         CALL COSGEN (LR, 1, FNUM2, 0.5, TCOS(K4+1))
         CALL MERGE (TCOS, K3, JR, K4, LR, K1 + K2)
         CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS(K3+1))
         K3 = KR
         K4 = KR
         CALL TRI3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
         B(:MR) = B(:MR) + B2(:MR) + B3(:MR)
         IF (NP == 3) THEN
            TCOS(1) = 2.
            CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
         ENDIF
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         B(:MR) = Q(:MR,1) + Q(:MR,J)
         CALL COSGEN (JR, 1, FNUM, 0.5, TCOS)
         CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         IF (JR == 1) THEN
            Q(:MR,1) = B(:MR)
            GO TO 188
         ENDIF
         Q(:MR,1) = Q(:MR,1) - Q(:MR,JM1) + B(:MR)
         GO TO 188
      ENDIF
      B3(:MR) = 0.
      B(:MR) = Q(:MR,1) + P(IP+1:MR+IP)
      Q(:MR,1) = Q(:MR,1) - Q(:MR,JM1)
      B2(:MR) = Q(:MR,1) + Q(:MR,NLAST)
      K1 = KR + JR
      K2 = K1 + JR
      CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(K1+1))
      GO TO (182,183,182) NP
  182 CONTINUE
      TCOS(K2) = 2.*FLOAT(NP - 2)
      CALL COSGEN (KR, 1, 0.0, 1.0, TCOS(K2+1))
      GO TO 184
  183 CONTINUE
      CALL COSGEN (KR + 1, 1, 0.5, 0.0, TCOS(K2))
  184 CONTINUE
      K4 = 1 - NP/3
      CALL MERGE (TCOS, K1, JR - K4, K2 - K4, KR + K4, 0)
      IF (NP == 3) K1 = K1 - 1
      K2 = KR
      CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS(K1+1))
      K4 = K1 + KR
      CALL COSGEN (LR, 1, FNUM2, 0.5, TCOS(K4+1))
      K3 = LR
      K4 = 0
      CALL TRI3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
      B(:MR) = B(:MR) + B2(:MR)
      IF (NP == 3) THEN
         TCOS(1) = 2.
         CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
      ENDIF
      Q(:MR,1) = Q(:MR,1) + B(:MR)
  188 CONTINUE
      J = NLAST - JR
      B(:MR) = Q(:MR,NLAST) + Q(:MR,J)
      JM2 = NLAST - I2R
      IF (JR == 1) THEN
         Q(:MR,NLAST) = 0.
      ELSE
         IF (NROD == 0) THEN
            Q(:MR,NLAST) = P(IP+1:MR+IP)
            IP = IP - MR
         ELSE
            Q(:MR,NLAST) = Q(:MR,NLAST) - Q(:MR,JM2)
         ENDIF
      ENDIF
      CALL COSGEN (KR, 1, FNUM2, 0.5, TCOS)
      CALL COSGEN (LR, 1, FNUM2, 0.5, TCOS(KR+1))
      CALL TRIX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
      Q(:MR,NLAST) = Q(:MR,NLAST) + B(:MR)
      NLASTP = NLAST
  197 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR == 0) GO TO 210
      JSTART = 1 + JR
      KR = KR - JR
      IF (NLAST + JR <= N) THEN
         KR = KR - JR
         NLAST = NLAST + JR
         JSTOP = NLAST - JSTEP
      ELSE
         JSTOP = NLAST - JR
      ENDIF
      LR = KR - JR
      CALL COSGEN (JR, 1, 0.5, 0.0, TCOS)
      DO J = JSTART, JSTOP, JSTEP
         JM2 = J - JR
         JP2 = J + JR
         IF (J == JR) THEN
            B(:MR) = Q(:MR,J) + Q(:MR,JP2)
         ELSE
            B(:MR) = Q(:MR,J) + Q(:MR,JM2) + Q(:MR,JP2)
         ENDIF
         IF (JR == 1) THEN
            Q(:MR,J) = 0.
         ELSE
            JM1 = J - I2R
            JP1 = J + I2R
            Q(:MR,J) = 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1))
         ENDIF
         CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
      END DO
      NROD = 1
      IF (NLAST + I2R <= N) NROD = 0
      IF (NLASTP /= NLAST) GO TO 188
      GO TO 197
  210 CONTINUE
      W(1) = IPSTOR
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
      END SUBROUTINE POSTG2
