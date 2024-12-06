C
C     file genbun.f
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
C     SUBROUTINE GENBUN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR)
C
C
C DIMENSION OF           A(M),B(M),C(M),Y(IDIMY,N)
C ARGUMENTS
C
C LATEST REVISION        JUNE 2004
C
C PURPOSE                THE NAME OF THIS PACKAGE IS A MNEMONIC FOR THE
C                        GENERALIZED BUNEMAN ALGORITHM.
C
C                        IT SOLVES THE REAL LINEAR SYSTEM OF EQUATIONS
C
C                        A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
C                        + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
C
C                        FOR I = 1,2,...,M  AND  J = 1,2,...,N.
C
C                        INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
C                        I.E., X(0,J) = X(M,J) AND X(M+1,J) = X(1,J),
C                        AND X(I,0) MAY EQUAL 0, X(I,2), OR X(I,N),
C                        AND X(I,N+1) MAY EQUAL 0, X(I,N-1), OR X(I,1)
C                        DEPENDING ON AN INPUT PARAMETER.
C
C USAGE                  CALL GENBUN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,
C                                     IERROR)
C
C ARGUMENTS
C
C ON INPUT               NPEROD
C
C                          INDICATES THE VALUES THAT X(I,0) AND
C                          X(I,N+1) ARE ASSUMED TO HAVE.
C
C                          = 0  IF X(I,0) = X(I,N) AND X(I,N+1) =
C                               X(I,1).
C                          = 1  IF X(I,0) = X(I,N+1) = 0  .
C                          = 2  IF X(I,0) = 0 AND X(I,N+1) = X(I,N-1).
C                          = 3  IF X(I,0) = X(I,2) AND X(I,N+1) =
C                               X(I,N-1).
C                          = 4  IF X(I,0) = X(I,2) AND X(I,N+1) = 0.
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
C                          N MUST BE GREATER THAN 2.
C
C                        A,B,C
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
C                          SPECIFY THE COEFFICIENTS IN THE LINEAR
C                          EQUATIONS GIVEN ABOVE.  IF MPEROD = 0
C                          THE ARRAY ELEMENTS MUST NOT DEPEND UPON
C                          THE INDEX I, BUT MUST BE CONSTANT.
C                          SPECIFICALLY, THE SUBROUTINE CHECKS THE
C                          FOLLOWING CONDITION .
C
C                            A(I) = C(1)
C                            C(I) = C(1)
C                            B(I) = B(1)
C
C                          FOR I=1,2,...,M.
C
C                        IDIMY
C                          THE ROW (OR FIRST) DIMENSION OF THE
C                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS
C                          IN THE PROGRAM CALLING GENBUN.
C                          THIS PARAMETER IS USED TO SPECIFY THE
C                          VARIABLE DIMENSION OF Y.
C                          IDIMY MUST BE AT LEAST M.
C
C                        Y
C                          A TWO-DIMENSIONAL COMPLEX ARRAY THAT
C                          SPECIFIES THE VALUES OF THE RIGHT SIDE
C                          OF THE LINEAR SYSTEM OF EQUATIONS GIVEN
C                          ABOVE.
C                          Y MUST BE DIMENSIONED AT LEAST M*N.
C
C
C  ON OUTPUT             Y
C
C                          CONTAINS THE SOLUTION X.
C
C                        IERROR
C                          AN ERROR FLAG WHICH INDICATES INVALID
C                          INPUT PARAMETERS  EXCEPT FOR NUMBER
C                          ZERO, A SOLUTION IS NOT ATTEMPTED.
C
C                          = 0  NO ERROR.
C                          = 1  M .LE. 2  .
C                          = 2  N .LE. 2
C                          = 3  IDIMY .LT. M
C                          = 4  NPEROD .LT. 0 OR NPEROD .GT. 4
C                          = 5  MPEROD .LT. 0 OR MPEROD .GT. 1
C                          = 6  A(I) .NE. C(1) OR C(I) .NE. C(1) OR
C                               B(I) .NE. B(1) FOR
C                               SOME I=1,2,...,M.
C                          = 7  A(1) .NE. 0 OR C(M) .NE. 0 AND
C                                 MPEROD = 1
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C
C SPECIAL CONDITONS      NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED FILES         comf.f,gnbnaux.f,fish.f
C FILES
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN IN 1979 BY ROLAND SWEET OF NCAR'S
C                        SCIENTIFIC COMPUTING DIVISION.  MADE AVAILABLE
C                        ON NCAR'S PUBLIC LIBRARIES IN JANUARY, 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C ALGORITHM              THE LINEAR SYSTEM IS SOLVED BY A CYCLIC
C                        REDUCTION ALGORITHM DESCRIBED IN THE
C                        REFERENCE.
C
C PORTABILITY            FORTRAN 90 --
C                        THE MACHINE DEPENDENT CONSTANT PI IS
C                        DEFINED IN FUNCTION PIMACH.
C
C REFERENCES             SWEET, R., "A CYCLIC REDUCTION ALGORITHM FOR
C                        SOLVING BLOCK TRIDIAGONAL SYSTEMS OF ARBITRARY
C                        DIMENSIONS," SIAM J. ON NUMER. ANAL., 14 (1977)
C                        PP. 706-720.
C
C ACCURACY               THIS TEST WAS PERFORMED ON a platform with
c                        64 bit floating point arithmetic.
C                        A UNIFORM RANDOM NUMBER GENERATOR WAS USED
C                        TO CREATE A SOLUTION ARRAY X FOR THE SYSTEM
C                        GIVEN IN THE 'PURPOSE' DESCRIPTION ABOVE
C                        WITH
C                          A(I) = C(I) = -0.5*B(I) = 1, I=1,2,...,M
C
C                        AND, WHEN MPEROD = 1
C
C                          A(1) = C(M) = 0
C                          A(M) = C(1) = 2.
C
C                        THE SOLUTION X WAS SUBSTITUTED INTO THE
C                        GIVEN SYSTEM  AND, USING DOUBLE PRECISION
C                        A RIGHT SIDE Y WAS COMPUTED.
C                        USING THIS ARRAY Y, SUBROUTINE GENBUN
C                        WAS CALLED TO PRODUCE APPROXIMATE
C                        SOLUTION Z.  THEN RELATIVE ERROR
C                          E = MAX(ABS(Z(I,J)-X(I,J)))/
C                              MAX(ABS(X(I,J)))
C                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
C                        OVER I=1,2,...,M AND J=1,...,N.
C
C                        THE VALUE OF E IS GIVEN IN THE TABLE
C                        BELOW FOR SOME TYPICAL VALUES OF M AND N.
C
C                   M (=N)    MPEROD    NPEROD        E
C                   ------    ------    ------      ------
C
C                     31        0         0         6.E-14
C                     31        1         1         4.E-13
C                     31        1         3         3.E-13
C                     32        0         0         9.E-14
C                     32        1         1         3.E-13
C                     32        1         3         1.E-13
C                     33        0         0         9.E-14
C                     33        1         1         4.E-13
C                     33        1         3         1.E-13
C                     63        0         0         1.E-13
C                     63        1         1         1.E-12
C                     63        1         3         2.E-13
C                     64        0         0         1.E-13
C                     64        1         1         1.E-12
C                     64        1         3         6.E-13
C                     65        0         0         2.E-13
C                     65        1         1         1.E-12
C                     65        1         3         4.E-13
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE GENBUN(NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, IERROR)
      USE fish
      implicit none
      TYPE(fishworkspace) :: w
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
      INTEGER :: IRWK, ICWK
C-----------------------------------------------
      IERROR = 0
!     check input arguments
      IF (M <= 2) then
	 ierror = 1
	 return
      end if
      IF (N <= 2) then
	 ierror = 2
	 return
      end if
      IF (IDIMY < M) then
	 ierror = 3
	 return
      end if
      IF (NPEROD<0 .OR. NPEROD>4) then
	 ierror = 4
	 return
      end if
      IF (MPEROD<0 .OR. MPEROD>1) then
	 ierror = 5
	 return
      end if
!     compute and allocate real work space for genbun
      CALL GEN_SPACE (N, M, IRWK)
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     return if allocation failed (e.g., if n,m are too large)
      IF (IERROR == 20) RETURN 
      call genbunn(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE GENBUN

 
      SUBROUTINE GENBUNN(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER  :: N
      INTEGER , INTENT(IN) :: MPEROD
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER , INTENT(INOUT) :: IERROR
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL  :: Y(IDIMY,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, MP1, IWBA, IWBB, IWBC, IWB2, IWB3, IWW1, IWW2, IWW3
     1   , IWD, IWTCOS, IWP, K, J, MP, NP, IPSTOR, IREV, MH, MHM1, MODD
     2   , MHPI, MHMI, NBY2, MSKIP
      REAL :: ALL, A1

      SAVE ALL
C-----------------------------------------------
      IF (MPEROD /= 1) THEN
         DO I = 2, M
            IF (A(I) /= C(1)) GO TO 103
            IF (C(I) /= C(1)) GO TO 103
            IF (B(I) /= B(1)) GO TO 103
         END DO
         GO TO 104
      ENDIF
      IF (A(1)/=0. .OR. C(M)/=0.) IERROR = 7
      GO TO 104
  103 CONTINUE
      IERROR = 6
  104 CONTINUE
      IF (IERROR /= 0) RETURN 
      MP1 = M + 1
      IWBA = MP1
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
      W(IWBA:M-1+IWBA) = -A(:M)
      W(IWBC:M-1+IWBC) = -C(:M)
      W(IWBB:M-1+IWBB) = 2. - B(:M)
      Y(:M,:N) = -Y(:M,:N)
      MP = MPEROD + 1
      NP = NPEROD + 1
      GO TO (114,107) MP
  107 CONTINUE
      GO TO (108,109,110,111,123) NP
  108 CONTINUE
      CALL POISP2 (M, N, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W(IWB2)
     1   , W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), W(IWP)
     2   )
      GO TO 112
  109 CONTINUE
      CALL POISD2 (M, N, 1, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W(
     1   IWW1), W(IWD), W(IWTCOS), W(IWP))
      GO TO 112
  110 CONTINUE
      CALL POISN2 (M, N, 1, 2, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W
     1   (IWB2), W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), 
     2   W(IWP))
      GO TO 112
  111 CONTINUE
      CALL POISN2 (M, N, 1, 1, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W
     1   (IWB2), W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), 
     2   W(IWP))
  112 CONTINUE
      IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD == 4) GO TO 124
  113 CONTINUE
      GO TO (127,133) MP
  114 CONTINUE
      MH = (M + 1)/2
      MHM1 = MH - 1
      MODD = 1
      IF (MH*2 == M) MODD = 2
      DO J = 1, N
         W(:MHM1) = Y(MH-1:MH-MHM1:(-1),J) - Y(MH+1:MHM1+MH,J)
         W(MH+1:MHM1+MH) = Y(MH-1:MH-MHM1:(-1),J) + Y(MH+1:MHM1+MH,J)
         W(MH) = 2.*Y(MH,J)
         GO TO (117,116) MODD
  116    CONTINUE
         W(M) = 2.*Y(M,J)
  117    CONTINUE
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
  122 CONTINUE
      GO TO 107
C
C     REVERSE COLUMNS WHEN NPEROD = 4.
C
  123 CONTINUE
      IREV = 1
      NBY2 = N/2
  124 CONTINUE
      DO J = 1, NBY2
         MSKIP = N + 1 - J
         DO I = 1, M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
         END DO
      END DO
      GO TO (110,113) IREV
  127 CONTINUE
      DO J = 1, N
         W(MH-1:MH-MHM1:(-1)) = 0.5*(Y(MH+1:MHM1+MH,J)+Y(:MHM1,J))
         W(MH+1:MHM1+MH) = 0.5*(Y(MH+1:MHM1+MH,J)-Y(:MHM1,J))
         W(MH) = 0.5*Y(MH,J)
         GO TO (130,129) MODD
  129    CONTINUE
         W(M) = 0.5*Y(M,J)
  130    CONTINUE
         Y(:M,J) = W(:M)
      END DO
  133 CONTINUE
      W(1) = IPSTOR + IWP - 1
      RETURN 
      END SUBROUTINE GENBUNN


      SUBROUTINE POISD2(MR,NR,ISTAG,BA,BB,BC,Q,IDIMQ,B,W,D,TCOS,P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: MR
      INTEGER , INTENT(IN) :: NR
      INTEGER , INTENT(IN) :: ISTAG
      INTEGER , INTENT(IN) :: IDIMQ
      REAL  :: BA(*)
      REAL  :: BB(*)
      REAL  :: BC(*)
      REAL , INTENT(INOUT) :: Q(IDIMQ,1)
      REAL  :: B(*)
      REAL  :: W(*)
      REAL  :: D(*)
      REAL  :: TCOS(*)
      REAL , INTENT(INOUT) :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: M, N, JSH, IP, IPSTOR, KR, IRREG, JSTSAV, I, LR, NUN, 
     1   JST, JSP, L, NODD, J, JM1, JP1, JM2, JP2, JM3, JP3, NODDPR, IP1
     2   , KRPI, IDEG, JDEG
      REAL :: ALL, FI, T

      SAVE ALL
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION FOR DIRICHLET BOUNDARY
C     CONDITIONS.
C
C     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A.
C     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A+I.
C
      M = MR
      N = NR
      JSH = 0
      FI = 1./FLOAT(ISTAG)
      IP = -M
      IPSTOR = 0
      SELECT CASE (ISTAG) 
      CASE DEFAULT
         KR = 0
         IRREG = 1
         IF (N > 1) GO TO 106
         TCOS(1) = 0.
      CASE (2) 
         KR = 1
         JSTSAV = 1
         IRREG = 2
         IF (N > 1) GO TO 106
         TCOS(1) = -1.
      END SELECT
  103 CONTINUE
      B(:M) = Q(:M,1)
      CALL TRIX (1, 0, M, BA, BB, BC, B, TCOS, D, W)
      Q(:M,1) = B(:M)
      GO TO 183
  106 CONTINUE
      LR = 0
      P(:M) = 0.
      NUN = N
      JST = 1
      JSP = N
C
C     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
C
  108 CONTINUE
      L = 2*JST
      NODD = 2 - 2*((NUN + 1)/2) + NUN
C
C     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
C
      SELECT CASE (NODD) 
      CASE DEFAULT
         JSP = JSP - L
      CASE (1) 
         JSP = JSP - JST
         IF (IRREG /= 1) JSP = JSP - L
      END SELECT
  111 CONTINUE
      CALL COSGEN (JST, 1, 0.5, 0.0, TCOS)
      IF (L <= JSP) THEN
         DO J = L, JSP, L
            JM1 = J - JSH
            JP1 = J + JSH
            JM2 = J - JST
            JP2 = J + JST
            JM3 = JM2 - JSH
            JP3 = JP2 + JSH
            IF (JST == 1) THEN
               B(:M) = 2.*Q(:M,J)
               Q(:M,J) = Q(:M,JM2) + Q(:M,JP2)
            ELSE
               DO I = 1, M
                  T = Q(I,J) - Q(I,JM1) - Q(I,JP1) + Q(I,JM2) + Q(I,JP2)
                  B(I) = T + Q(I,J) - Q(I,JM3) - Q(I,JP3)
                  Q(I,J) = T
               END DO
            ENDIF
            CALL TRIX (JST, 0, M, BA, BB, BC, B, TCOS, D, W)
            Q(:M,J) = Q(:M,J) + B(:M)
         END DO
      ENDIF
C
C     REDUCTION FOR LAST UNKNOWN
C
      SELECT CASE (NODD) 
      CASE DEFAULT
         GO TO (152,120) IRREG
C
C     ODD NUMBER OF UNKNOWNS
C
  120    CONTINUE
         JSP = JSP + L
         J = JSP
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         JM3 = JM2 - JSH
         GO TO (123,121) ISTAG
  121    CONTINUE
         IF (JST /= 1) GO TO 123
         B(:M) = Q(:M,J)
         Q(:M,J) = 0.
         GO TO 130
  123    CONTINUE
         SELECT CASE (NODDPR) 
         CASE DEFAULT
            B(:M) = 0.5*(Q(:M,JM2)-Q(:M,JM1)-Q(:M,JM3)) + P(IP+1:M+IP)
     1          + Q(:M,J)
         CASE (2) 
            B(:M) = 0.5*(Q(:M,JM2)-Q(:M,JM1)-Q(:M,JM3)) + Q(:M,JP2) - Q(
     1         :M,JP1) + Q(:M,J)
         END SELECT
  128    CONTINUE
         Q(:M,J) = 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1))
  130    CONTINUE
         CALL TRIX (JST, 0, M, BA, BB, BC, B, TCOS, D, W)
         IP = IP + M
         IPSTOR = MAX0(IPSTOR,IP + M)
         P(IP+1:M+IP) = Q(:M,J) + B(:M)
         B(:M) = Q(:M,JP2) + P(IP+1:M+IP)
         IF (LR == 0) THEN
            DO I = 1, JST
               KRPI = KR + I
               TCOS(KRPI) = TCOS(I)
            END DO
         ELSE
            CALL COSGEN (LR, JSTSAV, 0., FI, TCOS(JST+1))
            CALL MERGE (TCOS, 0, JST, JST, LR, KR)
         ENDIF
         CALL COSGEN (KR, JSTSAV, 0.0, FI, TCOS)
         CALL TRIX (KR, KR, M, BA, BB, BC, B, TCOS, D, W)
         Q(:M,J) = Q(:M,JM2) + B(:M) + P(IP+1:M+IP)
         LR = KR
         KR = KR + L
C
C     EVEN NUMBER OF UNKNOWNS
C
      CASE (2) 
         JSP = JSP + L
         J = JSP
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         JM3 = JM2 - JSH
         SELECT CASE (IRREG) 
         CASE DEFAULT
            JSTSAV = JST
            IDEG = JST
            KR = L
         CASE (2) 
            CALL COSGEN (KR, JSTSAV, 0.0, FI, TCOS)
            CALL COSGEN (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
            IDEG = KR
            KR = KR + JST
         END SELECT
  139    CONTINUE
         IF (JST == 1) THEN
            IRREG = 2
            B(:M) = Q(:M,J)
            Q(:M,J) = Q(:M,JM2)
         ELSE
            B(:M) = Q(:M,J) + 0.5*(Q(:M,JM2)-Q(:M,JM1)-Q(:M,JM3))
            SELECT CASE (IRREG) 
            CASE DEFAULT
               Q(:M,J) = Q(:M,JM2) + 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1))
               IRREG = 2
            CASE (2) 
               SELECT CASE (NODDPR) 
               CASE DEFAULT
                  Q(:M,J) = Q(:M,JM2) + P(IP+1:M+IP)
                  IP = IP - M
               CASE (2) 
                  Q(:M,J) = Q(:M,JM2) + Q(:M,J) - Q(:M,JM1)
               END SELECT
            END SELECT
         ENDIF
  150    CONTINUE
         CALL TRIX (IDEG, LR, M, BA, BB, BC, B, TCOS, D, W)
         Q(:M,J) = Q(:M,J) + B(:M)
      END SELECT
  152 CONTINUE
      NUN = NUN/2
      NODDPR = NODD
      JSH = JST
      JST = 2*JST
      IF (NUN >= 2) GO TO 108
C
C     START SOLUTION.
C
      J = JSP
      B(:M) = Q(:M,J)
      SELECT CASE (IRREG) 
      CASE DEFAULT
         CALL COSGEN (JST, 1, 0.5, 0.0, TCOS)
         IDEG = JST
      CASE (2) 
         KR = LR + JST
         CALL COSGEN (KR, JSTSAV, 0.0, FI, TCOS)
         CALL COSGEN (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
         IDEG = KR
      END SELECT
  156 CONTINUE
      CALL TRIX (IDEG, LR, M, BA, BB, BC, B, TCOS, D, W)
      JM1 = J - JSH
      JP1 = J + JSH
      SELECT CASE (IRREG) 
      CASE DEFAULT
         Q(:M,J) = 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1)) + B(:M)
      CASE (2) 
         SELECT CASE (NODDPR) 
         CASE DEFAULT
            Q(:M,J) = P(IP+1:M+IP) + B(:M)
            IP = IP - M
         CASE (2) 
            Q(:M,J) = Q(:M,J) - Q(:M,JM1) + B(:M)
         END SELECT
      END SELECT
  164 CONTINUE
      JST = JST/2
      JSH = JST/2
      NUN = 2*NUN
      IF (NUN > N) GO TO 183
      DO J = JST, N, L
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         IF (J <= JST) THEN
            B(:M) = Q(:M,J) + Q(:M,JP2)
         ELSE
            IF (JP2 <= N) GO TO 168
            B(:M) = Q(:M,J) + Q(:M,JM2)
            IF (JST < JSTSAV) IRREG = 1
            GO TO (170,171) IRREG
  168       CONTINUE
            B(:M) = Q(:M,J) + Q(:M,JM2) + Q(:M,JP2)
         ENDIF
  170    CONTINUE
         CALL COSGEN (JST, 1, 0.5, 0.0, TCOS)
         IDEG = JST
         JDEG = 0
         GO TO 172
  171    CONTINUE
         IF (J + L > N) LR = LR - JST
         KR = JST + LR
         CALL COSGEN (KR, JSTSAV, 0.0, FI, TCOS)
         CALL COSGEN (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
         IDEG = KR
         JDEG = LR
  172    CONTINUE
         CALL TRIX (IDEG, JDEG, M, BA, BB, BC, B, TCOS, D, W)
         IF (JST <= 1) THEN
            Q(:M,J) = B(:M)
         ELSE
            IF (JP2 > N) GO TO 177
  175       CONTINUE
            Q(:M,J) = 0.5*(Q(:M,J)-Q(:M,JM1)-Q(:M,JP1)) + B(:M)
            CYCLE 
  177       CONTINUE
            GO TO (175,178) IRREG
  178       CONTINUE
            IF (J + JSH <= N) THEN
               Q(:M,J) = B(:M) + P(IP+1:M+IP)
               IP = IP - M
            ELSE
               Q(:M,J) = B(:M) + Q(:M,J) - Q(:M,JM1)
            ENDIF
         ENDIF
      END DO
      L = L/2
      GO TO 164
  183 CONTINUE
      W(1) = IPSTOR
      RETURN 
      END SUBROUTINE POISD2


      SUBROUTINE POISN2(M, N, ISTAG, MIXBND, A, BB, C, Q, IDIMQ, B, B2, 
     1   B3, W, W2, W3, D, TCOS, P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: ISTAG
      INTEGER , INTENT(IN) :: MIXBND
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
      INTEGER :: K1, K2, K3, K4, MR, IP, IPSTOR, I2R, JR, NR, NLAST, KR
     1   , LR, I, NROD, JSTART, JSTOP, I2RBY2, J, JP1, JP2, JP3, JM1, 
     2   JM2, JM3, NRODPR, II, I1, I2, JR2, NLASTP, JSTEP
      REAL :: ALL, FISTAG, FNUM, FDEN, FI, T

      SAVE ALL
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION WITH NEUMANN BOUNDARY
C     CONDITIONS.
C
C     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS A.
C     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS A-I.
C     MIXBND = 1 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTH BOUNDARIES.
C     MIXBND = 2 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTTOM AND
C     DIRICHLET CONDITION AT TOP.  (FOR THIS CASE, MUST HAVE ISTAG = 1.)
C
      EQUIVALENCE (K(1), K1), (K(2), K2), (K(3), K3), (K(4), K4)
      FISTAG = 3 - ISTAG
      FNUM = 1./FLOAT(ISTAG)
      FDEN = 0.5*FLOAT(ISTAG - 1)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      GO TO (101,103) ISTAG
  101 CONTINUE
      Q(:MR,N) = 0.5*Q(:MR,N)
      GO TO (103,104) MIXBND
  103 CONTINUE
      IF (N <= 3) GO TO 155
  104 CONTINUE
      JR = 2*I2R
      NROD = 1
      IF ((NR/2)*2 == NR) NROD = 0
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         JSTART = 1
      CASE (2) 
         JSTART = JR
         NROD = 1 - NROD
      END SELECT
  107 CONTINUE
      JSTOP = NLAST - JR
      IF (NROD == 0) JSTOP = JSTOP - I2R
      CALL COSGEN (I2R, 1, 0.5, 0.0, TCOS)
      I2RBY2 = I2R/2
      IF (JSTOP < JSTART) THEN
         J = JR
      ELSE
         DO J = JSTART, JSTOP, JR
            JP1 = J + I2RBY2
            JP2 = J + I2R
            JP3 = JP2 + I2RBY2
            JM1 = J - I2RBY2
            JM2 = J - I2R
            JM3 = JM2 - I2RBY2
            IF (J == 1) THEN
               JM1 = JP1
               JM2 = JP2
               JM3 = JP3
            ENDIF
            IF (I2R == 1) THEN
               IF (J == 1) JM2 = JP2
               B(:MR) = 2.*Q(:MR,J)
               Q(:MR,J) = Q(:MR,JM2) + Q(:MR,JP2)
            ELSE
               DO I = 1, MR
                  FI = Q(I,J)
                  Q(I,J)=Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
                  B(I) = FI + Q(I,J) - Q(I,JM3) - Q(I,JP3)
               END DO
            ENDIF
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
            B(:MR) = FISTAG*Q(:MR,J)
            Q(:MR,J) = Q(:MR,JM2)
         ELSE
            B(:MR) = Q(:MR,J) + 0.5*(Q(:MR,JM2)-Q(:MR,JM1)-Q(:MR,JM3))
            IF (NRODPR == 0) THEN
               Q(:MR,J) = Q(:MR,JM2) + P(IP+1:MR+IP)
               IP = IP - MR
            ELSE
               Q(:MR,J) = Q(:MR,J) - Q(:MR,JM1) + Q(:MR,JM2)
            ENDIF
            IF (LR /= 0) THEN
               CALL COSGEN (LR, 1, 0.5, FDEN, TCOS(KR+1))
            ELSE
               B(:MR) = FISTAG*B(:MR)
            ENDIF
         ENDIF
         CALL COSGEN (KR, 1, 0.5, FDEN, TCOS)
         CALL TRIX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         KR = KR + I2R
      ELSE
         JP1 = J + I2RBY2
         JP2 = J + I2R
         IF (I2R == 1) THEN
            B(:MR) = Q(:MR,J)
            CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            IP = 0
            IPSTOR = MR
            SELECT CASE (ISTAG) 
            CASE DEFAULT
               P(:MR) = B(:MR)
               B(:MR) = B(:MR) + Q(:MR,N)
               TCOS(1) = 1.
               TCOS(2) = 0.
               CALL TRIX (1, 1, MR, A, BB, C, B, TCOS, D, W)
               Q(:MR,J) = Q(:MR,JM2) + P(:MR) + B(:MR)
               GO TO 150
            CASE (1) 
               P(:MR) = B(:MR)
               Q(:MR,J) = Q(:MR,JM2) + 2.*Q(:MR,JP2) + 3.*B(:MR)
               GO TO 150
            END SELECT
         ENDIF
         B(:MR) = Q(:MR,J) + 0.5*(Q(:MR,JM2)-Q(:MR,JM1)-Q(:MR,JM3))
         IF (NRODPR == 0) THEN
            B(:MR) = B(:MR) + P(IP+1:MR+IP)
         ELSE
            B(:MR) = B(:MR) + Q(:MR,JP2) - Q(:MR,JP1)
         ENDIF
         CALL TRIX (I2R, 0, MR, A, BB, C, B, TCOS, D, W)
         IP = IP + MR
         IPSTOR = MAX0(IPSTOR,IP + MR)
         P(IP+1:MR+IP) = B(:MR) + 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1))
         B(:MR) = P(IP+1:MR+IP) + Q(:MR,JP2)
         IF (LR /= 0) THEN
            CALL COSGEN (LR, 1, 0.5, FDEN, TCOS(I2R+1))
            CALL MERGE (TCOS, 0, I2R, I2R, LR, KR)
         ELSE
            DO I = 1, I2R
               II = KR + I
               TCOS(II) = TCOS(I)
            END DO
         ENDIF
         CALL COSGEN (KR, 1, 0.5, FDEN, TCOS)
         IF (LR == 0) THEN
            GO TO (146,145) ISTAG
         ENDIF
  145    CONTINUE
         CALL TRIX (KR, KR, MR, A, BB, C, B, TCOS, D, W)
         GO TO 148
  146    CONTINUE
         B(:MR) = FISTAG*B(:MR)
  148    CONTINUE
         Q(:MR,J) = Q(:MR,JM2) + P(IP+1:MR+IP) + B(:MR)
  150    CONTINUE
         LR = KR
         KR = KR + JR
      ENDIF
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         NR = (NLAST - 1)/JR + 1
         IF (NR <= 3) GO TO 155
      CASE (2) 
         NR = NLAST/JR
         IF (NR <= 1) GO TO 192
      END SELECT
  154 CONTINUE
      I2R = JR
      NRODPR = NROD
      GO TO 104
  155 CONTINUE
      J = 1 + JR
      JM1 = J - I2R
      JP1 = J + I2R
      JM2 = NLAST - I2R
      IF (NR /= 2) THEN
         IF (LR /= 0) GO TO 170
         IF (N == 3) THEN
C
C     CASE N = 3.
C
            GO TO (156,168) ISTAG
  156       CONTINUE
            B(:MR) = Q(:MR,2)
            TCOS(1) = 0.
            CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,2) = B(:MR)
            B(:MR) = 4.*B(:MR) + Q(:MR,1) + 2.*Q(:MR,3)
            TCOS(1) = -2.
            TCOS(2) = 2.
            I1 = 2
            I2 = 0
            CALL TRIX (I1, I2, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,2) = Q(:MR,2) + B(:MR)
            B(:MR) = Q(:MR,1) + 2.*Q(:MR,2)
            TCOS(1) = 0.
            CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,1) = B(:MR)
            JR = 1
            I2R = 0
            GO TO 194
         ENDIF
C
C     CASE N = 2**P+1
C
         GO TO (162,170) ISTAG
  162    CONTINUE
         B(:MR) = Q(:MR,J) + 0.5*Q(:MR,1) - Q(:MR,JM1) + Q(:MR,NLAST) - 
     1      Q(:MR,JM2)
         CALL COSGEN (JR, 1, 0.5, 0.0, TCOS)
         CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1)) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,NLAST) + 4.*Q(:MR,J)
         JR2 = 2*JR
         CALL COSGEN (JR, 1, 0.0, 0.0, TCOS)
         TCOS(JR+1:JR*2) = -TCOS(JR:1:(-1))
         CALL TRIX (JR2, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,J)
         CALL COSGEN (JR, 1, 0.5, 0.0, TCOS)
         CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = 0.5*Q(:MR,1) - Q(:MR,JM1) + B(:MR)
         GO TO 194
C
C     CASE OF GENERAL N WITH NR = 3 .
C
  168    CONTINUE
         B(:MR) = Q(:MR,2)
         Q(:MR,2) = 0.
         B2(:MR) = Q(:MR,3)
         B3(:MR) = Q(:MR,1)
         JR = 1
         I2R = 0
         J = 2
         GO TO 177
  170    CONTINUE
         B(:MR) = 0.5*Q(:MR,1) - Q(:MR,JM1) + Q(:MR,J)
         IF (NROD == 0) THEN
            B(:MR) = B(:MR) + P(IP+1:MR+IP)
         ELSE
            B(:MR) = B(:MR) + Q(:MR,NLAST) - Q(:MR,JM2)
         ENDIF
         DO I = 1, MR
            T = 0.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
            Q(I,J) = T
            B2(I) = Q(I,NLAST) + T
            B3(I) = Q(I,1) + 2.*T
         END DO
  177    CONTINUE
         K1 = KR + 2*JR - 1
         K2 = KR + JR
         TCOS(K1+1) = -2.
         K4 = K1 + 3 - ISTAG
         CALL COSGEN (K2 + ISTAG - 2, 1, 0.0, FNUM, TCOS(K4))
         K4 = K1 + K2 + 1
         CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(K4))
         CALL MERGE (TCOS, K1, K2, K1 + K2, JR - 1, 0)
         K3 = K1 + K2 + LR
         CALL COSGEN (JR, 1, 0.5, 0.0, TCOS(K3+1))
         K4 = K3 + JR + 1
         CALL COSGEN (KR, 1, 0.5, FDEN, TCOS(K4))
         CALL MERGE (TCOS, K3, JR, K3 + JR, KR, K1)
         IF (LR /= 0) THEN
            CALL COSGEN (LR, 1, 0.5, FDEN, TCOS(K4))
            CALL MERGE (TCOS, K3, JR, K3 + JR, LR, K3 - LR)
            CALL COSGEN (KR, 1, 0.5, FDEN, TCOS(K4))
         ENDIF
         K3 = KR
         K4 = KR
         CALL TRI3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
         B(:MR) = B(:MR) + B2(:MR) + B3(:MR)
         TCOS(1) = 2.
         CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,J)
         CALL COSGEN (JR, 1, 0.5, 0.0, TCOS)
         CALL TRIX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         IF (JR == 1) THEN
            Q(:MR,1) = B(:MR)
            GO TO 194
         ENDIF
         Q(:MR,1) = 0.5*Q(:MR,1) - Q(:MR,JM1) + B(:MR)
         GO TO 194
      ENDIF
      IF (N == 2) THEN
C
C     CASE  N = 2
C
         B(:MR) = Q(:MR,1)
         TCOS(1) = 0.
         CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = B(:MR)
         B(:MR) = 2.*(Q(:MR,2)+B(:MR))*FISTAG
         TCOS(1) = -FISTAG
         TCOS(2) = 2.
         CALL TRIX (2, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = Q(:MR,1) + B(:MR)
         JR = 1
         I2R = 0
         GO TO 194
      ENDIF
      B3(:MR) = 0.
      B(:MR) = Q(:MR,1) + 2.*P(IP+1:MR+IP)
      Q(:MR,1) = 0.5*Q(:MR,1) - Q(:MR,JM1)
      B2(:MR) = 2.*(Q(:MR,1)+Q(:MR,NLAST))
      K1 = KR + JR - 1
      TCOS(K1+1) = -2.
      K4 = K1 + 3 - ISTAG
      CALL COSGEN (KR + ISTAG - 2, 1, 0.0, FNUM, TCOS(K4))
      K4 = K1 + KR + 1
      CALL COSGEN (JR - 1, 1, 0.0, 1.0, TCOS(K4))
      CALL MERGE (TCOS, K1, KR, K1 + KR, JR - 1, 0)
      CALL COSGEN (KR, 1, 0.5, FDEN, TCOS(K1+1))
      K2 = KR
      K4 = K1 + K2 + 1
      CALL COSGEN (LR, 1, 0.5, FDEN, TCOS(K4))
      K3 = LR
      K4 = 0
      CALL TRI3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
      B(:MR) = B(:MR) + B2(:MR)
      TCOS(1) = 2.
      CALL TRIX (1, 0, MR, A, BB, C, B, TCOS, D, W)
      Q(:MR,1) = Q(:MR,1) + B(:MR)
      GO TO 194
  192 CONTINUE
      B(:MR) = Q(:MR,NLAST)
      GO TO 196
  194 CONTINUE
      J = NLAST - JR
      B(:MR) = Q(:MR,NLAST) + Q(:MR,J)
  196 CONTINUE
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
      CALL COSGEN (KR, 1, 0.5, FDEN, TCOS)
      CALL COSGEN (LR, 1, 0.5, FDEN, TCOS(KR+1))
      IF (LR == 0) THEN
         B(:MR) = FISTAG*B(:MR)
      ENDIF
      CALL TRIX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
      Q(:MR,NLAST) = Q(:MR,NLAST) + B(:MR)
      NLASTP = NLAST
  206 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR == 0) GO TO 222
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         JSTART = 1 + JR
      CASE (2) 
         JSTART = JR
      END SELECT
  209 CONTINUE
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
      IF (NLASTP /= NLAST) GO TO 194
      GO TO 206
  222 CONTINUE
      W(1) = IPSTOR
      RETURN 
      END SUBROUTINE POISN2


      SUBROUTINE POISP2(M,N,A,BB,C,Q,IDIMQ,B,B2,B3,W,W2,W3,D,TCOS,P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER  :: IDIMQ
      REAL  :: A(*)
      REAL  :: BB(*)
      REAL  :: C(*)
      REAL  :: Q(IDIMQ,1)
      REAL  :: B(*)
      REAL  :: B2(*)
      REAL  :: B3(*)
      REAL  :: W(*)
      REAL  :: W2(*)
      REAL  :: W3(*)
      REAL  :: D(*)
      REAL  :: TCOS(*)
      REAL  :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MR, NR, NRM1, J, NRMJ, NRPJ, I, IPSTOR, LH
      REAL :: ALL, S, T

      SAVE ALL
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE POISSON EQUATION WITH PERIODIC BOUNDARY
C     CONDITIONS.
C
      MR = M
      NR = (N + 1)/2
      NRM1 = NR - 1
      IF (2*NR == N) THEN
C
C     EVEN NUMBER OF UNKNOWNS
C
         DO J = 1, NRM1
            NRMJ = NR - J
            NRPJ = NR + J
            DO I = 1, MR
               S = Q(I,NRMJ) - Q(I,NRPJ)
               T = Q(I,NRMJ) + Q(I,NRPJ)
               Q(I,NRMJ) = S
               Q(I,NRPJ) = T
            END DO
         END DO
         Q(:MR,NR) = 2.*Q(:MR,NR)
         Q(:MR,N) = 2.*Q(:MR,N)
         CALL POISD2 (MR, NRM1, 1, A, BB, C, Q, IDIMQ, B, W, D, TCOS, P)
         IPSTOR = W(1)
         CALL POISN2 (MR, NR + 1, 1, 1, A, BB, C, Q(1,NR), IDIMQ, B, B2
     1      , B3, W, W2, W3, D, TCOS, P)
         IPSTOR = MAX0(IPSTOR,INT(W(1)))
         DO J = 1, NRM1
            NRMJ = NR - J
            NRPJ = NR + J
            DO I = 1, MR
               S = 0.5*(Q(I,NRPJ)+Q(I,NRMJ))
               T = 0.5*(Q(I,NRPJ)-Q(I,NRMJ))
               Q(I,NRMJ) = S
               Q(I,NRPJ) = T
            END DO
         END DO
         Q(:MR,NR) = 0.5*Q(:MR,NR)
         Q(:MR,N) = 0.5*Q(:MR,N)
      ELSE
         DO J = 1, NRM1
            NRPJ = N + 1 - J
            DO I = 1, MR
               S = Q(I,J) - Q(I,NRPJ)
               T = Q(I,J) + Q(I,NRPJ)
               Q(I,J) = S
               Q(I,NRPJ) = T
            END DO
         END DO
         Q(:MR,NR) = 2.*Q(:MR,NR)
         LH = NRM1/2
         DO J = 1, LH
            NRMJ = NR - J
            DO I = 1, MR
               S = Q(I,J)
               Q(I,J) = Q(I,NRMJ)
               Q(I,NRMJ) = S
            END DO
         END DO
         CALL POISD2 (MR, NRM1, 2, A, BB, C, Q, IDIMQ, B, W, D, TCOS, P)
         IPSTOR = W(1)
         CALL POISN2 (MR, NR, 2, 1, A, BB, C, Q(1,NR), IDIMQ, B, B2, B3
     1      , W, W2, W3, D, TCOS, P)
         IPSTOR = MAX0(IPSTOR,INT(W(1)))
         DO J = 1, NRM1
            NRPJ = NR + J
            DO I = 1, MR
               S = 0.5*(Q(I,NRPJ)+Q(I,J))
               T = 0.5*(Q(I,NRPJ)-Q(I,J))
               Q(I,NRPJ) = T
               Q(I,J) = S
            END DO
         END DO
         Q(:MR,NR) = 0.5*Q(:MR,NR)
         DO J = 1, LH
            NRMJ = NR - J
            DO I = 1, MR
               S = Q(I,J)
               Q(I,J) = Q(I,NRMJ)
               Q(I,NRMJ) = S
            END DO
         END DO
      ENDIF
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
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE POISP2
