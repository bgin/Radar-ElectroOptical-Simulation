C
C     file gnbnaux.f
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
C
C PACKAGE GNBNAUX
C
C LATEST REVISION        June 2004
C
C PURPOSE                TO PROVIDE AUXILIARY ROUTINES FOR FISHPACK
C                        ENTRIES GENBUN AND POISTG.
C
C USAGE                  THERE ARE NO USER ENTRIES IN THIS PACKAGE.
C                        THE ROUTINES IN THIS PACKAGE ARE NOT INTENDED
C                        TO BE CALLED BY USERS, BUT RATHER BY ROUTINES
C                        IN PACKAGES GENBUN AND POISTG.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN IN 1979 BY ROLAND SWEET OF NCAR'S
C                        SCIENTIFIC COMPUTING DIVISION.  MADE AVAILABLE
C                        ON NCAR'S PUBLIC LIBRARIES IN JANUARY, 1980.
c                        Revised by John Adams in June 2004 incorporating
c                        Fortran 90 features
C
C PORTABILITY            FORTRAN 90
C ********************************************************************
      SUBROUTINE COSGEN(N, IJUMP, FNUM, FDEN, A)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: IJUMP
      REAL , INTENT(IN) :: FNUM
      REAL , INTENT(IN) :: FDEN
      REAL , INTENT(OUT) :: A(*)
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
C                                IN POISN2 ONLY.
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
                  A(K2) = -2.*COS(X*PIBYN)
               END DO
            END DO
         ELSE
            NP1 = N + 1
            Y = PI/(FLOAT(N) + FDEN)
            DO I = 1, N
               X = FLOAT(NP1 - I) - FNUM
               A(I) = 2.*COS(X*Y)
            END DO
         ENDIF
      ENDIF
      RETURN 
      END SUBROUTINE COSGEN


      SUBROUTINE MERGE(TCOS, I1, M1, I2, M2, I3)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I1
      INTEGER , INTENT(IN) :: M1
      INTEGER , INTENT(IN) :: I2
      INTEGER , INTENT(IN) :: M2
      INTEGER , INTENT(IN) :: I3
      REAL , INTENT(INOUT) :: TCOS(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J11, J3, J1, J2, J, L, K, M
      REAL :: X, Y
C-----------------------------------------------
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
         IF (X - Y > 0.) GO TO 103
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
      END SUBROUTINE MERGE


      SUBROUTINE TRIX(IDEGBR, IDEGCR, M, A, B, C, Y, TCOS, D, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDEGBR
      INTEGER , INTENT(IN) :: IDEGCR
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL , INTENT(INOUT) :: Y(*)
      REAL , INTENT(IN) :: TCOS(*)
      REAL , INTENT(INOUT) :: D(*)
      REAL , INTENT(INOUT) :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MM1, IFB, IFC, L, LINT, K, I, IP
      REAL :: X, XX, Z
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
         IF (Z == 0.) THEN
            Y(M) = 0.
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
      END SUBROUTINE TRIX


      SUBROUTINE TRI3(M, A, B, C, K, Y1, Y2, Y3, TCOS, D, W1, W2, W3)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: K(4)
      REAL , INTENT(IN) :: A(*)
      REAL , INTENT(IN) :: B(*)
      REAL , INTENT(IN) :: C(*)
      REAL , INTENT(INOUT) :: Y1(*)
      REAL , INTENT(INOUT) :: Y2(*)
      REAL , INTENT(INOUT) :: Y3(*)
      REAL , INTENT(IN) :: TCOS(*)
      REAL , INTENT(INOUT) :: D(*)
      REAL , INTENT(INOUT) :: W1(*)
      REAL , INTENT(INOUT) :: W2(*)
      REAL , INTENT(INOUT) :: W3(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MM1, K1, K2, K3, K4, IF1, IF2, IF3, IF4, K2K3K4, L1, L2
     1   , L3, LINT1, LINT2, LINT3, KINT1, KINT2, KINT3, N, I, IP
      REAL :: X, Z, XX
C-----------------------------------------------
C
C     SUBROUTINE TO SOLVE THREE LINEAR SYSTEMS WHOSE COMMON COEFFICIENT
C     MATRIX IS A RATIONAL FUNCTION IN THE MATRIX GIVEN BY
C
C                  TRIDIAGONAL (...,A(I),B(I),C(I),...)
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
C OCTOBER   1980    CHANGED SEVERAL DIVIDES OF FLOATING INTEGERS
C                   TO INTEGER DIVIDES TO ACCOMODATE CRAY-1 ARITHMETIC.
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
C-----------------------------------------------------------------------
      END SUBROUTINE TRI3
