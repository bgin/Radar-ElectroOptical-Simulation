C
C     file cmgnbn.f
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
C     SUBROUTINE CMGNBN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR)
C
C
C DIMENSION OF           A(M),B(M),C(M),Y(IDIMY,N)
C ARGUMENTS
C
C LATEST REVISION        NOVEMBER 2004
C
C PURPOSE                THE NAME OF THIS PACKAGE IS A MNEMONIC FOR THE
C                        COMPLEX GENERALIZED BUNEMAN ALGORITHM.
C                        IT SOLVES THE COMPLEX LINEAR SYSTEM OF EQUATION
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
C USAGE                  CALL CMGNBN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,
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
C                          ONE-DIMENSIONAL COMPLEX ARRAYS OF LENGTH M
C                          THAT SPECIFY THE COEFFICIENTS IN THE LINEAR
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
C                          IN THE PROGRAM CALLING CMGNBN.
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
C SPECIAL CONDITONS      NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       comf.f,fish.f
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
C                        REFERENCE BELOW.
C
C PORTABILITY            FORTRAN 90.  ALL MACHINE DEPENDENT CONSTANTS
C                        ARE DEFINED IN FUNCTION P1MACH.
C
C REFERENCES             SWEET, R., 'A CYCLIC REDUCTION ALGORITHM FOR
C                        SOLVING BLOCK TRIDIAGONAL SYSTEMS OF ARBITRARY
C                        DIMENSIONS,' SIAM J. ON NUMER. ANAL.,
C                          14(SEPT., 1977), PP. 706-720.
C
C ACCURACY               THIS TEST WAS PERFORMED ON A Platform with
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
C                        GIVEN SYSTEM  AND A RIGHT SIDE Y WAS
C                        COMPUTED.  USING THIS ARRAY Y, SUBROUTINE
C                        CMGNBN WAS CALLED TO PRODUCE APPROXIMATE
C                        SOLUTION Z.  THEN RELATIVE ERROR
C                          E = MAX(CABS(Z(I,J)-X(I,J)))/
C                              MAX(CABS(X(I,J)))
C                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
C                        OVER I=1,2,...,M AND J=1,...,N.
C
C                        THE VALUE OF E IS GIVEN IN THE TABLE
C                        BELOW FOR SOME TYPICAL VALUES OF M AND N.
C
C                   M (=N)    MPEROD    NPEROD       E
C                   ------    ------    ------     ------
C
C                     31        0         0        1.E-12
C                     31        1         1        4.E-13
C                     31        1         3        2.E-12
C                     32        0         0        7.E-14
C                     32        1         1        5.E-13
C                     32        1         3        2.E-13
C                     33        0         0        6.E-13
C                     33        1         1        5.E-13
C                     33        1         3        3.E-12
C                     63        0         0        5.E-12
C                     63        1         1        6.E-13
C                     63        1         3        1.E-11
C                     64        0         0        1.E-13
C                     64        1         1        3.E-12
C                     64        1         3        3.E-13
C                     65        0         0        2.E-12
C                     65        1         1        5.E-13
C                     65        1         3        1.E-11
C
C***********************************************************************
      SUBROUTINE CMGNBN(NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, IERROR)
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
      COMPLEX  :: A(*)
      COMPLEX  :: B(*)
      COMPLEX  :: C(*)
      COMPLEX  :: Y(IDIMY,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, ICWK, IRWK
      COMPLEX :: A1
C-----------------------------------------------
      IERROR = 0
      IF (M <= 2) IERROR = 1
      IF (N <= 2) IERROR = 2
      IF (IDIMY < M) IERROR = 3
      IF (NPEROD<0 .OR. NPEROD>4) IERROR = 4
      IF (MPEROD<0 .OR. MPEROD>1) IERROR = 5
      IF (MPEROD /= 1) THEN
         DO I = 2, M
            IF (CABS(A(I)-C(1)) /= 0.) GO TO 103
            IF (CABS(C(I)-C(1)) /= 0.) GO TO 103
            IF (CABS(B(I)-B(1)) /= 0.) GO TO 103
         END DO
         GO TO 104
      ENDIF
      IF (CABS(A(1))/=0. .AND. CABS(C(M))/=0.) IERROR = 7
      GO TO 104
  103 CONTINUE
      IERROR = 6
  104 CONTINUE
      IF (IERROR /= 0) RETURN 
!     allocate required complex work space
      ICWK = (10 + INT(ALOG(FLOAT(N))/ALOG(2.0)))*M + 4*N
      IRWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     return if allocation failed
      IF (IERROR == 20) RETURN 
      call cmgnbnn(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,w%cxw)
!     release dynamically allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE CMGNBN


 
      SUBROUTINE CMGNBNN(NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER  :: N
      INTEGER , INTENT(IN) :: MPEROD
      INTEGER  :: M
      INTEGER  :: IDIMY
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX  :: Y(IDIMY,*)
      COMPLEX  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IWBA, IWBB, IWBC, IWB2, IWB3, IWW1, IWW2, IWW3, IWD, 
     1   IWTCOS, IWP, I, K, J, MP, NP, IPSTOR, IREV, MH, MHM1, MODD, 
     2   MHPI, MHMI, NBY2, MSKIP
      COMPLEX :: A1
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
      MP = MPEROD + 1
      NP = NPEROD + 1
      GO TO (114,107) MP
  107 CONTINUE
      GO TO (108,109,110,111,123) NP
  108 CONTINUE
      CALL CMPOSP (M, N, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W(IWB2)
     1   , W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), W(IWP)
     2   )
      GO TO 112
  109 CONTINUE
      CALL CMPOSD (M, N, 1, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W(
     1   IWW1), W(IWD), W(IWTCOS), W(IWP))
      GO TO 112
  110 CONTINUE
      CALL CMPOSN (M, N, 1, 2, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W
     1   (IWB2), W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), 
     2   W(IWP))
      GO TO 112
  111 CONTINUE
      CALL CMPOSN (M, N, 1, 1, W(IWBA), W(IWBB), W(IWBC), Y, IDIMY, W, W
     1   (IWB2), W(IWB3), W(IWW1), W(IWW2), W(IWW3), W(IWD), W(IWTCOS), 
     2   W(IWP))
  112 CONTINUE
      IPSTOR = REAL(W(IWW1))
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
         DO I = 1, MHM1
            W(I) = Y(MH-I,J) - Y(I+MH,J)
            W(I+MH) = Y(MH-I,J) + Y(I+MH,J)
         END DO
         W(MH) = 2.*Y(MH,J)
         GO TO (117,116) MODD
  116    CONTINUE
         W(M) = 2.*Y(M,J)
  117    CONTINUE
         Y(:M,J) = W(:M)
      END DO
      K = IWBC + MHM1 - 1
      I = IWBA + MHM1
      W(K) = (0.,0.)
      W(I) = (0.,0.)
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
C     REVERSE COLUMNS WHEN NPEROD = 4
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
      W(1) = CMPLX(FLOAT(IPSTOR + IWP - 1),0.)
      RETURN 
      END SUBROUTINE CMGNBNN


      SUBROUTINE CMPOSD(MR,NR,ISTAG,BA,BB,BC,Q,IDIMQ,B,W,D,TCOS,P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: MR
      INTEGER , INTENT(IN) :: NR
      INTEGER , INTENT(IN) :: ISTAG
      INTEGER , INTENT(IN) :: IDIMQ
      COMPLEX  :: BA(*)
      COMPLEX  :: BB(*)
      COMPLEX  :: BC(*)
      COMPLEX , INTENT(INOUT) :: Q(IDIMQ,1)
      COMPLEX  :: B(*)
      COMPLEX  :: W(*)
      COMPLEX  :: D(*)
      COMPLEX  :: TCOS(*)
      COMPLEX , INTENT(INOUT) :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: M, N, IP, IPSTOR, JSH, KR, IRREG, JSTSAV, I, LR, NUN, 
     1   JST, JSP, L, NODD, J, JM1, JP1, JM2, JP2, JM3, JP3, NODDPR, IP1
     2   , KRPI, IDEG, JDEG
      REAL :: FI
      COMPLEX :: T
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
      FI = 1./FLOAT(ISTAG)
      IP = -M
      IPSTOR = 0
      JSH = 0
      SELECT CASE (ISTAG) 
      CASE DEFAULT
         KR = 0
         IRREG = 1
         IF (N > 1) GO TO 106
         TCOS(1) = (0.,0.)
      CASE (2) 
         KR = 1
         JSTSAV = 1
         IRREG = 2
         IF (N > 1) GO TO 106
         TCOS(1) = CMPLX(-1.,0.)
      END SELECT
  103 CONTINUE
      B(:M) = Q(:M,1)
      CALL CMPTRX (1, 0, M, BA, BB, BC, B, TCOS, D, W)
      Q(:M,1) = B(:M)
      GO TO 183
  106 CONTINUE
      LR = 0
      DO I = 1, M
         P(I) = CMPLX(0.,0.)
      END DO
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
      CALL CMPCSG (JST, 1, 0.5, 0.0, TCOS)
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
            CALL CMPTRX (JST, 0, M, BA, BB, BC, B, TCOS, D, W)
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
         DO I = 1, M
            B(I) = Q(I,J)
            Q(I,J) = CMPLX(0.,0.)
         END DO
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
         CALL CMPTRX (JST, 0, M, BA, BB, BC, B, TCOS, D, W)
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
            CALL CMPCSG (LR, JSTSAV, 0., FI, TCOS(JST+1))
            CALL CMPMRG (TCOS, 0, JST, JST, LR, KR)
         ENDIF
         CALL CMPCSG (KR, JSTSAV, 0.0, FI, TCOS)
         CALL CMPTRX (KR, KR, M, BA, BB, BC, B, TCOS, D, W)
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
            CALL CMPCSG (KR, JSTSAV, 0.0, FI, TCOS)
            CALL CMPCSG (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
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
         CALL CMPTRX (IDEG, LR, M, BA, BB, BC, B, TCOS, D, W)
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
         CALL CMPCSG (JST, 1, 0.5, 0.0, TCOS)
         IDEG = JST
      CASE (2) 
         KR = LR + JST
         CALL CMPCSG (KR, JSTSAV, 0.0, FI, TCOS)
         CALL CMPCSG (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
         IDEG = KR
      END SELECT
  156 CONTINUE
      CALL CMPTRX (IDEG, LR, M, BA, BB, BC, B, TCOS, D, W)
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
         CALL CMPCSG (JST, 1, 0.5, 0.0, TCOS)
         IDEG = JST
         JDEG = 0
         GO TO 172
  171    CONTINUE
         IF (J + L > N) LR = LR - JST
         KR = JST + LR
         CALL CMPCSG (KR, JSTSAV, 0.0, FI, TCOS)
         CALL CMPCSG (LR, JSTSAV, 0.0, FI, TCOS(KR+1))
         IDEG = KR
         JDEG = LR
  172    CONTINUE
         CALL CMPTRX (IDEG, JDEG, M, BA, BB, BC, B, TCOS, D, W)
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
      W(1) = CMPLX(FLOAT(IPSTOR),0.)
      RETURN 
      END SUBROUTINE CMPOSD


      SUBROUTINE CMPOSN(M, N, ISTAG, MIXBND, A, BB, C, Q, IDIMQ, B, B2, 
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
      COMPLEX  :: A(*)
      COMPLEX  :: BB(*)
      COMPLEX  :: C(*)
      COMPLEX , INTENT(INOUT) :: Q(IDIMQ,*)
      COMPLEX  :: B(*)
      COMPLEX  :: B2(*)
      COMPLEX  :: B3(*)
      COMPLEX  :: W(*)
      COMPLEX  :: W2(*)
      COMPLEX  :: W3(*)
      COMPLEX  :: D(*)
      COMPLEX  :: TCOS(*)
      COMPLEX , INTENT(INOUT) :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER , DIMENSION(4) :: K
      INTEGER :: K1, K2, K3, K4, MR, IP, IPSTOR, I2R, JR, NR, NLAST, KR
     1   , LR, I, NROD, JSTART, JSTOP, I2RBY2, J, JP1, JP2, JP3, JM1, 
     2   JM2, JM3, NRODPR, II, I1, I2, JR2, NLASTP, JSTEP
      REAL :: FISTAG, FNUM, FDEN
      COMPLEX :: FI, T
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
      CALL CMPCSG (I2R, 1, 0.5, 0.0, TCOS)
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
            CALL CMPTRX (I2R, 0, MR, A, BB, C, B, TCOS, D, W)
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
               CALL CMPCSG (LR, 1, 0.5, FDEN, TCOS(KR+1))
            ELSE
               B(:MR) = FISTAG*B(:MR)
            ENDIF
         ENDIF
         CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS)
         CALL CMPTRX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         KR = KR + I2R
      ELSE
         JP1 = J + I2RBY2
         JP2 = J + I2R
         IF (I2R == 1) THEN
            B(:MR) = Q(:MR,J)
            CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            IP = 0
            IPSTOR = MR
            SELECT CASE (ISTAG) 
            CASE DEFAULT
               P(:MR) = B(:MR)
               B(:MR) = B(:MR) + Q(:MR,N)
               TCOS(1) = CMPLX(1.,0.)
               TCOS(2) = CMPLX(0.,0.)
               CALL CMPTRX (1, 1, MR, A, BB, C, B, TCOS, D, W)
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
         CALL CMPTRX (I2R, 0, MR, A, BB, C, B, TCOS, D, W)
         IP = IP + MR
         IPSTOR = MAX0(IPSTOR,IP + MR)
         P(IP+1:MR+IP) = B(:MR) + 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1))
         B(:MR) = P(IP+1:MR+IP) + Q(:MR,JP2)
         IF (LR /= 0) THEN
            CALL CMPCSG (LR, 1, 0.5, FDEN, TCOS(I2R+1))
            CALL CMPMRG (TCOS, 0, I2R, I2R, LR, KR)
         ELSE
            DO I = 1, I2R
               II = KR + I
               TCOS(II) = TCOS(I)
            END DO
         ENDIF
         CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS)
         IF (LR == 0) THEN
            GO TO (146,145) ISTAG
         ENDIF
  145    CONTINUE
         CALL CMPTRX (KR, KR, MR, A, BB, C, B, TCOS, D, W)
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
            TCOS(1) = CMPLX(0.,0.)
            CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,2) = B(:MR)
            B(:MR) = 4.*B(:MR) + Q(:MR,1) + 2.*Q(:MR,3)
            TCOS(1) = CMPLX(-2.,0.)
            TCOS(2) = CMPLX(2.,0.)
            I1 = 2
            I2 = 0
            CALL CMPTRX (I1, I2, MR, A, BB, C, B, TCOS, D, W)
            Q(:MR,2) = Q(:MR,2) + B(:MR)
            B(:MR) = Q(:MR,1) + 2.*Q(:MR,2)
            TCOS(1) = (0.,0.)
            CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
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
         CALL CMPCSG (JR, 1, 0.5, 0.0, TCOS)
         CALL CMPTRX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1)) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,NLAST) + 4.*Q(:MR,J)
         JR2 = 2*JR
         CALL CMPCSG (JR, 1, 0.0, 0.0, TCOS)
         TCOS(JR+1:JR*2) = -TCOS(JR:1:(-1))
         CALL CMPTRX (JR2, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,J)
         CALL CMPCSG (JR, 1, 0.5, 0.0, TCOS)
         CALL CMPTRX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = 0.5*Q(:MR,1) - Q(:MR,JM1) + B(:MR)
         GO TO 194
C
C     CASE OF GENERAL N WITH NR = 3 .
C
  168    CONTINUE
         B(:MR) = Q(:MR,2)
         Q(:MR,2) = (0.,0.)
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
         TCOS(K1+1) = (-2.,0.)
         K4 = K1 + 3 - ISTAG
         CALL CMPCSG (K2 + ISTAG - 2, 1, 0.0, FNUM, TCOS(K4))
         K4 = K1 + K2 + 1
         CALL CMPCSG (JR - 1, 1, 0.0, 1.0, TCOS(K4))
         CALL CMPMRG (TCOS, K1, K2, K1 + K2, JR - 1, 0)
         K3 = K1 + K2 + LR
         CALL CMPCSG (JR, 1, 0.5, 0.0, TCOS(K3+1))
         K4 = K3 + JR + 1
         CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS(K4))
         CALL CMPMRG (TCOS, K3, JR, K3 + JR, KR, K1)
         IF (LR /= 0) THEN
            CALL CMPCSG (LR, 1, 0.5, FDEN, TCOS(K4))
            CALL CMPMRG (TCOS, K3, JR, K3 + JR, LR, K3 - LR)
            CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS(K4))
         ENDIF
         K3 = KR
         K4 = KR
         CALL CMPTR3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
         B(:MR) = B(:MR) + B2(:MR) + B3(:MR)
         TCOS(1) = (2.,0.)
         CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
         B(:MR) = Q(:MR,1) + 2.*Q(:MR,J)
         CALL CMPCSG (JR, 1, 0.5, 0.0, TCOS)
         CALL CMPTRX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
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
         TCOS(1) = (0.,0.)
         CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = B(:MR)
         B(:MR) = 2.*(Q(:MR,2)+B(:MR))*FISTAG
         TCOS(1) = CMPLX((-FISTAG),0.)
         TCOS(2) = CMPLX(2.,0.)
         CALL CMPTRX (2, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,1) = Q(:MR,1) + B(:MR)
         JR = 1
         I2R = 0
         GO TO 194
      ENDIF
      B3(:MR) = (0.,0.)
      B(:MR) = Q(:MR,1) + 2.*P(IP+1:MR+IP)
      Q(:MR,1) = 0.5*Q(:MR,1) - Q(:MR,JM1)
      B2(:MR) = 2.*(Q(:MR,1)+Q(:MR,NLAST))
      K1 = KR + JR - 1
      TCOS(K1+1) = (-2.,0.)
      K4 = K1 + 3 - ISTAG
      CALL CMPCSG (KR + ISTAG - 2, 1, 0.0, FNUM, TCOS(K4))
      K4 = K1 + KR + 1
      CALL CMPCSG (JR - 1, 1, 0.0, 1.0, TCOS(K4))
      CALL CMPMRG (TCOS, K1, KR, K1 + KR, JR - 1, 0)
      CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS(K1+1))
      K2 = KR
      K4 = K1 + K2 + 1
      CALL CMPCSG (LR, 1, 0.5, FDEN, TCOS(K4))
      K3 = LR
      K4 = 0
      CALL CMPTR3 (MR, A, BB, C, K, B, B2, B3, TCOS, D, W, W2, W3)
      B(:MR) = B(:MR) + B2(:MR)
      TCOS(1) = (2.,0.)
      CALL CMPTRX (1, 0, MR, A, BB, C, B, TCOS, D, W)
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
         Q(:MR,NLAST) = (0.,0.)
      ELSE
         IF (NROD == 0) THEN
            Q(:MR,NLAST) = P(IP+1:MR+IP)
            IP = IP - MR
         ELSE
            Q(:MR,NLAST) = Q(:MR,NLAST) - Q(:MR,JM2)
         ENDIF
      ENDIF
      CALL CMPCSG (KR, 1, 0.5, FDEN, TCOS)
      CALL CMPCSG (LR, 1, 0.5, FDEN, TCOS(KR+1))
      IF (LR == 0) THEN
         B(:MR) = FISTAG*B(:MR)
      ENDIF
      CALL CMPTRX (KR, LR, MR, A, BB, C, B, TCOS, D, W)
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
      CALL CMPCSG (JR, 1, 0.5, 0.0, TCOS)
      DO J = JSTART, JSTOP, JSTEP
         JM2 = J - JR
         JP2 = J + JR
         IF (J == JR) THEN
            B(:MR) = Q(:MR,J) + Q(:MR,JP2)
         ELSE
            B(:MR) = Q(:MR,J) + Q(:MR,JM2) + Q(:MR,JP2)
         ENDIF
         IF (JR == 1) THEN
            Q(:MR,J) = (0.,0.)
         ELSE
            JM1 = J - I2R
            JP1 = J + I2R
            Q(:MR,J) = 0.5*(Q(:MR,J)-Q(:MR,JM1)-Q(:MR,JP1))
         ENDIF
         CALL CMPTRX (JR, 0, MR, A, BB, C, B, TCOS, D, W)
         Q(:MR,J) = Q(:MR,J) + B(:MR)
      END DO
      NROD = 1
      IF (NLAST + I2R <= N) NROD = 0
      IF (NLASTP /= NLAST) GO TO 194
      GO TO 206
  222 CONTINUE
      W(1) = CMPLX(FLOAT(IPSTOR),0.)
      RETURN 
      END SUBROUTINE CMPOSN


      SUBROUTINE CMPOSP(M,N,A,BB,C,Q,IDIMQ,B,B2,B3,W,W2,W3,D,TCOS,P)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER  :: IDIMQ
      COMPLEX  :: A(*)
      COMPLEX  :: BB(*)
      COMPLEX  :: C(*)
      COMPLEX  :: Q(IDIMQ,1)
      COMPLEX  :: B(*)
      COMPLEX  :: B2(*)
      COMPLEX  :: B3(*)
      COMPLEX  :: W(*)
      COMPLEX  :: W2(*)
      COMPLEX  :: W3(*)
      COMPLEX  :: D(*)
      COMPLEX  :: TCOS(*)
      COMPLEX  :: P(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MR, NR, NRM1, J, NRMJ, NRPJ, I, IPSTOR, LH
      COMPLEX :: S, T
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
         CALL CMPOSD (MR, NRM1, 1, A, BB, C, Q, IDIMQ, B, W, D, TCOS, P)
         IPSTOR = REAL(W(1))
         CALL CMPOSN (MR, NR + 1, 1, 1, A, BB, C, Q(1,NR), IDIMQ, B, B2
     1      , B3, W, W2, W3, D, TCOS, P)
         IPSTOR = MAX0(IPSTOR,INT(REAL(W(1))))
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
         CALL CMPOSD (MR, NRM1, 2, A, BB, C, Q, IDIMQ, B, W, D, TCOS, P)
         IPSTOR = REAL(W(1))
         CALL CMPOSN (MR, NR, 2, 1, A, BB, C, Q(1,NR), IDIMQ, B, B2, B3
     1      , W, W2, W3, D, TCOS, P)
         IPSTOR = MAX0(IPSTOR,INT(REAL(W(1))))
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
      W(1) = CMPLX(FLOAT(IPSTOR),0.)
      RETURN 
      END SUBROUTINE CMPOSP


      SUBROUTINE CMPCSG(N, IJUMP, FNUM, FDEN, A)
      implicit none
      REAL PIMACH
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IJUMP
      REAL , INTENT(IN) :: FNUM
      REAL , INTENT(IN) :: FDEN
      COMPLEX , INTENT(OUT) :: A(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: K3, K4, K, K1, K5, I, K2, NP1
      REAL :: PI, DUM, PIBYN, X, Y
C-----------------------------------------------
C
C
C     THIS SUBROUTINE COMPUTES REQUIRED COSINE VALUES IN ASCENDING
C     ORDER.  WHEN IJUMP .GT. 1 THE ROUTINE COMPUTES VALUES
C
C        2*COS(J*PI/L) , J=1,2,...,L AND J .NE. 0(MOD N/IJUMP+1)
C
C     WHERE L = IJUMP*(N/IJUMP+1).
C
C
C     WHEN IJUMP = 1 IT COMPUTES
C
C            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... ,N
C
C     WHERE
C        FNUM = 0.5, FDEN = 0.0,  FOR REGULAR REDUCTION VALUES
C        FNUM = 0.0, FDEN = 1.0, FOR B-R AND C-R WHEN ISTAG = 1
C        FNUM = 0.0, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
C        FNUM = 0.5, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
C                                IN CMPOSN ONLY.
C
C
      PI = 4.0*ATAN(1.0)
      IF (N /= 0) THEN
         IF (IJUMP /= 1) THEN
            K3 = N/IJUMP + 1
            K4 = K3 - 1
            PIBYN = PI/FLOAT(N + IJUMP)
            DO K = 1, IJUMP
               K1 = (K - 1)*K3
               K5 = (K - 1)*K4
               DO I = 1, K4
                  X = K1 + I
                  K2 = K5 + I
                  A(K2) = CMPLX((-2.*COS(X*PIBYN)),0.)
               END DO
            END DO
         ELSE
            NP1 = N + 1
            Y = PI/(FLOAT(N) + FDEN)
            DO I = 1, N
               X = FLOAT(NP1 - I) - FNUM
               A(I) = CMPLX(2.*COS(X*Y),0.)
            END DO
         ENDIF
      ENDIF
      RETURN 
      END SUBROUTINE CMPCSG


      SUBROUTINE CMPMRG(TCOS, I1, M1, I2, M2, I3)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I1
      INTEGER , INTENT(IN) :: M1
      INTEGER , INTENT(IN) :: I2
      INTEGER , INTENT(IN) :: M2
      INTEGER , INTENT(IN) :: I3
      COMPLEX , INTENT(INOUT) :: TCOS(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J11, J3, J1, J2, J, L, K, M
      COMPLEX :: X, Y
C-----------------------------------------------
C
C
C     THIS SUBROUTINE MERGES TWO ASCENDING STRINGS OF NUMBERS IN THE
C     ARRAY TCOS.  THE FIRST STRING IS OF LENGTH M1 AND STARTS AT
C     TCOS(I1+1).  THE SECOND STRING IS OF LENGTH M2 AND STARTS AT
C     TCOS(I2+1).  THE MERGED STRING GOES INTO TCOS(I3+1).
C
C
      J1 = 1
      J2 = 1
      J = I3
      IF (M1 == 0) GO TO 107
      IF (M2 == 0) GO TO 104
  101 CONTINUE
      J11 = J1
      J3 = MAX(M1,J11)
      DO J1 = J11, J3
         J = J + 1
         L = J1 + I1
         X = TCOS(L)
         L = J2 + I2
         Y = TCOS(L)
         IF (REAL(X - Y) > 0.) GO TO 103
         TCOS(J) = X
      END DO
      GO TO 106
  103 CONTINUE
      TCOS(J) = Y
      J2 = J2 + 1
      IF (J2 <= M2) GO TO 101
      IF (J1 > M1) GO TO 109
  104 CONTINUE
      K = J - J1 + 1
      DO J = J1, M1
         M = K + J
         L = J + I1
         TCOS(M) = TCOS(L)
      END DO
      GO TO 109
  106 CONTINUE
      IF (J2 > M2) GO TO 109
  107 CONTINUE
      K = J - J2 + 1
      DO J = J2, M2
         M = K + J
         L = J + I2
         TCOS(M) = TCOS(L)
      END DO
  109 CONTINUE
      RETURN 
      END SUBROUTINE CMPMRG


      SUBROUTINE CMPTRX(IDEGBR, IDEGCR, M, A, B, C, Y, TCOS, D, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDEGBR
      INTEGER , INTENT(IN) :: IDEGCR
      INTEGER , INTENT(IN) :: M
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
      COMPLEX , INTENT(IN) :: TCOS(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MM1, IFB, IFC, L, LINT, K, I, IP
      COMPLEX :: X, XX, Z
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE A SYSTEM OF LINEAR EQUATIONS WHERE THE
C     COEFFICIENT MATRIX IS A RATIONAL FUNCTION IN THE MATRIX GIVEN BY
C     TRIDIAGONAL  ( . . . , A(I), B(I), C(I), . . . ).
C
      MM1 = M - 1
      IFB = IDEGBR + 1
      IFC = IDEGCR + 1
      L = IFB/IFC
      LINT = 1
      DO K = 1, IDEGBR
         X = TCOS(K)
         IF (K == L) THEN
            I = IDEGBR + LINT
            XX = X - TCOS(I)
            W(:M) = Y(:M)
            Y(:M) = XX*Y(:M)
         ENDIF
         Z = 1./(B(1)-X)
         D(1) = C(1)*Z
         Y(1) = Y(1)*Z
         DO I = 2, MM1
            Z = 1./(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y(I) = (Y(I)-A(I)*Y(I-1))*Z
         END DO
         Z = B(M) - X - A(M)*D(MM1)
         IF (CABS(Z) == 0.) THEN
            Y(M) = (0.,0.)
         ELSE
            Y(M) = (Y(M)-A(M)*Y(MM1))/Z
         ENDIF
         DO IP = 1, MM1
            Y(M-IP) = Y(M-IP) - D(M-IP)*Y(M+1-IP)
         END DO
         IF (K /= L) CYCLE 
         Y(:M) = Y(:M) + W(:M)
         LINT = LINT + 1
         L = (LINT*IFB)/IFC
      END DO
      RETURN 
      END SUBROUTINE CMPTRX


      SUBROUTINE CMPTR3(M, A, B, C, K, Y1, Y2, Y3, TCOS, D, W1, W2, W3)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: K(4)
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: Y1(*)
      COMPLEX , INTENT(INOUT) :: Y2(*)
      COMPLEX , INTENT(INOUT) :: Y3(*)
      COMPLEX , INTENT(IN) :: TCOS(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: W1(*)
      COMPLEX , INTENT(INOUT) :: W2(*)
      COMPLEX , INTENT(INOUT) :: W3(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MM1, K1, K2, K3, K4, IF1, IF2, IF3, IF4, K2K3K4, L1, L2
     1   , L3, LINT1, LINT2, LINT3, KINT1, KINT2, KINT3, N, I, IP
      COMPLEX :: X, XX, Z
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE TRIDIAGONAL SYSTEMS
C
      MM1 = M - 1
      K1 = K(1)
      K2 = K(2)
      K3 = K(3)
      K4 = K(4)
      IF1 = K1 + 1
      IF2 = K2 + 1
      IF3 = K3 + 1
      IF4 = K4 + 1
      K2K3K4 = K2 + K3 + K4
      IF (K2K3K4 /= 0) THEN
         L1 = IF1/IF2
         L2 = IF1/IF3
         L3 = IF1/IF4
         LINT1 = 1
         LINT2 = 1
         LINT3 = 1
         KINT1 = K1
         KINT2 = KINT1 + K2
         KINT3 = KINT2 + K3
      ENDIF
      DO N = 1, K1
         X = TCOS(N)
         IF (K2K3K4 /= 0) THEN
            IF (N == L1) THEN
               W1(:M) = Y1(:M)
            ENDIF
            IF (N == L2) THEN
               W2(:M) = Y2(:M)
            ENDIF
            IF (N == L3) THEN
               W3(:M) = Y3(:M)
            ENDIF
         ENDIF
         Z = 1./(B(1)-X)
         D(1) = C(1)*Z
         Y1(1) = Y1(1)*Z
         Y2(1) = Y2(1)*Z
         Y3(1) = Y3(1)*Z
         DO I = 2, M
            Z = 1./(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y1(I) = (Y1(I)-A(I)*Y1(I-1))*Z
            Y2(I) = (Y2(I)-A(I)*Y2(I-1))*Z
            Y3(I) = (Y3(I)-A(I)*Y3(I-1))*Z
         END DO
         DO IP = 1, MM1
            Y1(M-IP) = Y1(M-IP) - D(M-IP)*Y1(M+1-IP)
            Y2(M-IP) = Y2(M-IP) - D(M-IP)*Y2(M+1-IP)
            Y3(M-IP) = Y3(M-IP) - D(M-IP)*Y3(M+1-IP)
         END DO
         IF (K2K3K4 == 0) CYCLE 
         IF (N == L1) THEN
            I = LINT1 + KINT1
            XX = X - TCOS(I)
            Y1(:M) = XX*Y1(:M) + W1(:M)
            LINT1 = LINT1 + 1
            L1 = (LINT1*IF1)/IF2
         ENDIF
         IF (N == L2) THEN
            I = LINT2 + KINT2
            XX = X - TCOS(I)
            Y2(:M) = XX*Y2(:M) + W2(:M)
            LINT2 = LINT2 + 1
            L2 = (LINT2*IF1)/IF3
         ENDIF
         IF (N /= L3) CYCLE 
         I = LINT3 + KINT3
         XX = X - TCOS(I)
         Y3(:M) = XX*Y3(:M) + W3(:M)
         LINT3 = LINT3 + 1
         L3 = (LINT3*IF1)/IF4
      END DO
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
      END SUBROUTINE CMPTR3
