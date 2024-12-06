C
C     file sepaux.f
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
C PACKAGE SEPAUX         CONTAINS NO USER ENTRY POINTS.
C
C LATEST REVISION        June 2004
C
C PURPOSE                THIS PACKAGE CONTAINS AUXILIARY ROUTINES FOR
C                        THE FISHPACK SOLVERS SEPELI AND SEPX4.
C
C USAGE                  SINCE THIS PACKAGE CONTAINS NO USER ENTRIES,
C                        NO USAGE INSTRUCTIONS OR ARGUMENT DESCRIPTIONS
C                        ARE GIVEN HERE.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       NONE
C FILES
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                DEVELOPED IN THE LATE 1970'S BY JOHN C. ADAMS
C                        OF NCAR'S SCIENTTIFIC COMPUTING DIVISION.
c                        Revised in June 2004 incorporating fortran 90
c                        features
C
C PORTABILITY            FORTRAN 90
C **********************************************************************
      SUBROUTINE SEPORT(USOL, IDMN, ZN, ZM, PERTRB)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDMN
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(INOUT) :: USOL(IDMN,1)
      REAL , INTENT(IN) :: ZN(*)
      REAL , INTENT(IN) :: ZM(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ISTR, IFNL, JSTR, JFNL, I, II, J, JJ
      REAL :: UTE, ETE
C-----------------------------------------------
C
C     THIS SUBROUTINE ORTHOGANALIZES THE ARRAY USOL WITH RESPECT TO
C     THE CONSTANT ARRAY IN A WEIGHTED LEAST SQUARES NORM
C
      ISTR = IS
      IFNL = MS
      JSTR = JS
      JFNL = NS
C
C     COMPUTE WEIGHTED INNER PRODUCTS
C
      UTE = 0.0
      ETE = 0.0
      DO I = IS, MS
         II = I - IS + 1
         ETE = ETE + SUM(ZM(II)*ZN(:NS-JS+1))
         UTE = UTE + SUM(USOL(I,JS:NS)*ZM(II)*ZN(:NS-JS+1))
      END DO
C
C     SET PERTURBATION PARAMETER
C
      PERTRB = UTE/ETE
C
C     SUBTRACT OFF CONSTANT PERTRB
C
      USOL(ISTR:IFNL,JSTR:JFNL) = USOL(ISTR:IFNL,JSTR:JFNL) - PERTRB
      RETURN 
      END SUBROUTINE SEPORT


      SUBROUTINE SEPMIN(USOL, IDMN, ZN, ZM, PERTB)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDMN
      REAL  :: PERTB
      REAL , INTENT(INOUT) :: USOL(IDMN,1)
      REAL , INTENT(IN) :: ZN(*)
      REAL , INTENT(IN) :: ZM(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ISTR, IFNL, JSTR, JFNL, I, II, J, JJ
      REAL :: UTE, ETE, PERTRB
C-----------------------------------------------
C
C     THIS SUBROUTINE ORHTOGONALIZES THE ARRAY USOL WITH RESPECT TO
C     THE CONSTANT ARRAY IN A WEIGHTED LEAST SQUARES NORM
C
C
C     ENTRY AT SEPMIN OCCURRS WHEN THE FINAL SOLUTION IS
C     TO BE MINIMIZED WITH RESPECT TO THE WEIGHTED
C     LEAST SQUARES NORM
C
      ISTR = 1
      IFNL = K
      JSTR = 1
      JFNL = L
C
C     COMPUTE WEIGHTED INNER PRODUCTS
C
      UTE = 0.0
      ETE = 0.0
      DO I = IS, MS
         II = I - IS + 1
         ETE = ETE + SUM(ZM(II)*ZN(:NS-JS+1))
         UTE = UTE + SUM(USOL(I,JS:NS)*ZM(II)*ZN(:NS-JS+1))
      END DO
C
C     SET PERTURBATION PARAMETER
C
      PERTRB = UTE/ETE
C
C     SUBTRACT OFF CONSTANT PERTRB
C
      USOL(ISTR:IFNL,JSTR:JFNL) = USOL(ISTR:IFNL,JSTR:JFNL) - PERTRB
      RETURN 
      END SUBROUTINE SEPMIN


      SUBROUTINE SEPTRI(N, A, B, C, D, U, Z)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: A(N)
      REAL , INTENT(IN) :: B(N)
      REAL , INTENT(IN) :: C(N)
      REAL , INTENT(INOUT) :: D(N)
      REAL , INTENT(INOUT) :: U(N)
      REAL , INTENT(INOUT) :: Z(N)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: NM2, J, NM1, K
      REAL :: BN, V, DEN, AN
C-----------------------------------------------
C
C     THIS SUBROUTINE SOLVES FOR A NON-ZERO EIGENVECTOR CORRESPONDING
C     TO THE ZERO EIGENVALUE OF THE TRANSPOSE OF THE RANK
C     DEFICIENT ONE MATRIX WITH SUBDIAGONAL A, DIAGONAL B, AND
C     SUPERDIAGONAL C , WITH A(1) IN THE (1,N) POSITION, WITH
C     C(N) IN THE (N,1) POSITION, AND ALL OTHER ELEMENTS ZERO.
C
      BN = B(N)
      D(1) = A(2)/B(1)
      V = A(1)
      U(1) = C(N)/B(1)
      NM2 = N - 2
      DO J = 2, NM2
         DEN = B(J) - C(J-1)*D(J-1)
         D(J) = A(J+1)/DEN
         U(J) = -C(J-1)*U(J-1)/DEN
         BN = BN - V*U(J-1)
         V = -V*D(J-1)
      END DO
      DEN = B(N-1) - C(N-2)*D(N-2)
      D(N-1) = (A(N)-C(N-2)*U(N-2))/DEN
      AN = C(N-1) - V*D(N-2)
      BN = BN - V*U(N-2)
      DEN = BN - AN*D(N-1)
C
C     SET LAST COMPONENT EQUAL TO ONE
C
      Z(N) = 1.0
      Z(N-1) = -D(N-1)
      NM1 = N - 1
      DO J = 2, NM1
         K = N - J
         Z(K) = (-D(K)*Z(K+1)) - U(K)*Z(N)
      END DO
      RETURN 
      END SUBROUTINE SEPTRI


      SUBROUTINE SEPDX(U, IDMN, I, J, UXXX, UXXXX)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDMN
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: J
      REAL , INTENT(OUT) :: UXXX
      REAL , INTENT(OUT) :: UXXXX
      REAL , INTENT(IN) :: U(IDMN,1)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
C
C     THIS PROGRAM COMPUTES SECOND ORDER FINITE DIFFERENCE
C     APPROXIMATIONS TO THE THIRD AND FOURTH X
C     PARTIAL DERIVATIVES OF U AT THE (I,J) MESH POINT
C
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A
C
      IF (I == 1) THEN
         IF (KSWX /= 1) THEN
            UXXX = ((-5.0*U(1,J))+18.0*U(2,J)-24.0*U(3,J)+14.0*U(4,J)-
     1         3.0*U(5,J))/TDLX3
            UXXXX = (3.0*U(1,J)-14.0*U(2,J)+26.0*U(3,J)-24.0*U(4,J)+11.0
     1         *U(5,J)-2.0*U(6,J))/DLX4
            RETURN 
         ELSE
C
C     PERIODIC AT X=A
C
            UXXX = ((-U(K-2,J))+2.0*U(K-1,J)-2.0*U(2,J)+U(3,J))/TDLX3
            UXXXX = (U(K-2,J)-4.0*U(K-1,J)+6.0*U(1,J)-4.0*U(2,J)+U(3,J))
     1         /DLX4
            RETURN 
         ENDIF
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A+DLX
C
      ELSE IF (I == 2) THEN
         IF (KSWX /= 1) THEN
            UXXX = ((-3.0*U(1,J))+10.0*U(2,J)-12.0*U(3,J)+6.0*U(4,J)-U(5
     1         ,J))/TDLX3
            UXXXX = (2.0*U(1,J)-9.0*U(2,J)+16.0*U(3,J)-14.0*U(4,J)+6.0*U
     1         (5,J)-U(6,J))/DLX4
            RETURN 
         ELSE
C
C     PERIODIC AT X=A+DLX
C
            UXXX = ((-U(K-1,J))+2.0*U(1,J)-2.0*U(3,J)+U(4,J))/TDLX3
            UXXXX = (U(K-1,J)-4.0*U(1,J)+6.0*U(2,J)-4.0*U(3,J)+U(4,J))/
     1         DLX4
            RETURN 
         ENDIF
      ELSE IF (I>2 .AND. I<K-1) THEN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
C
         UXXX = ((-U(I-2,J))+2.0*U(I-1,J)-2.0*U(I+1,J)+U(I+2,J))/TDLX3
         UXXXX = (U(I-2,J)-4.0*U(I-1,J)+6.0*U(I,J)-4.0*U(I+1,J)+U(I+2,J)
     1      )/DLX4
         RETURN 
      ELSE IF (I == K - 1) THEN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B-DLX
C
         IF (KSWX /= 1) THEN
            UXXX = (U(K-4,J)-6.0*U(K-3,J)+12.0*U(K-2,J)-10.0*U(K-1,J)+
     1         3.0*U(K,J))/TDLX3
            UXXXX = ((-U(K-5,J))+6.0*U(K-4,J)-14.0*U(K-3,J)+16.0*U(K-2,J
     1         )-9.0*U(K-1,J)+2.0*U(K,J))/DLX4
            RETURN 
         ELSE
C
C     PERIODIC AT X=B-DLX
C
            UXXX = ((-U(K-3,J))+2.0*U(K-2,J)-2.0*U(1,J)+U(2,J))/TDLX3
            UXXXX = (U(K-3,J)-4.0*U(K-2,J)+6.0*U(K-1,J)-4.0*U(1,J)+U(2,J
     1         ))/DLX4
            RETURN 
         ENDIF
      ELSE IF (I == K) THEN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B
C
         UXXX = -(3.0*U(K-4,J)-14.0*U(K-3,J)+24.0*U(K-2,J)-18.0*U(K-1,J)
     1      +5.0*U(K,J))/TDLX3
         UXXXX = ((-2.0*U(K-5,J))+11.0*U(K-4,J)-24.0*U(K-3,J)+26.0*U(K-2
     1      ,J)-14.0*U(K-1,J)+3.0*U(K,J))/DLX4
         RETURN 
      ENDIF
      RETURN 
      END SUBROUTINE SEPDX


      SUBROUTINE SEPDY(U, IDMN, I, J, UYYY, UYYYY)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IDMN
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: J
      REAL , INTENT(OUT) :: UYYY
      REAL , INTENT(OUT) :: UYYYY
      REAL , INTENT(IN) :: U(IDMN,6)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /SPLP/
      COMMON /SPLP/ KSWX, KSWY, K, L, AIT, BIT, CIT, DIT, MIT, NIT, IS, 
     1   MS, JS, NS, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      INTEGER   KSWX, KSWY, K, L, MIT, NIT, IS, MS, JS, NS
      REAL   AIT, BIT, CIT, DIT, DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
C
C     THIS PROGRAM COMPUTES SECOND ORDER FINITE DIFFERENCE
C     APPROXIMATIONS TO THE THIRD AND FOURTH Y
C     PARTIAL DERIVATIVES OF U AT THE (I,J) MESH POINT
C
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C
C
      IF (J == 1) THEN
         IF (KSWY /= 1) THEN
            UYYY = ((-5.0*U(I,1))+18.0*U(I,2)-24.0*U(I,3)+14.0*U(I,4)-
     1         3.0*U(I,5))/TDLY3
            UYYYY = (3.0*U(I,1)-14.0*U(I,2)+26.0*U(I,3)-24.0*U(I,4)+11.0
     1         *U(I,5)-2.0*U(I,6))/DLY4
            RETURN 
         ELSE
C
C     PERIODIC AT X=A
C
            UYYY = ((-U(I,L-2))+2.0*U(I,L-1)-2.0*U(I,2)+U(I,3))/TDLY3
            UYYYY = (U(I,L-2)-4.0*U(I,L-1)+6.0*U(I,1)-4.0*U(I,2)+U(I,3))
     1         /DLY4
            RETURN 
         ENDIF
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C+DLY
C
      ELSE IF (J == 2) THEN
         IF (KSWY /= 1) THEN
            UYYY = ((-3.0*U(I,1))+10.0*U(I,2)-12.0*U(I,3)+6.0*U(I,4)-U(I
     1         ,5))/TDLY3
            UYYYY = (2.0*U(I,1)-9.0*U(I,2)+16.0*U(I,3)-14.0*U(I,4)+6.0*U
     1         (I,5)-U(I,6))/DLY4
            RETURN 
         ELSE
C
C     PERIODIC AT Y=C+DLY
C
            UYYY = ((-U(I,L-1))+2.0*U(I,1)-2.0*U(I,3)+U(I,4))/TDLY3
            UYYYY = (U(I,L-1)-4.0*U(I,1)+6.0*U(I,2)-4.0*U(I,3)+U(I,4))/
     1         DLY4
            RETURN 
         ENDIF
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
C
      ELSE IF (J>2 .AND. J<L-1) THEN
         UYYY = ((-U(I,J-2))+2.0*U(I,J-1)-2.0*U(I,J+1)+U(I,J+2))/TDLY3
         UYYYY = (U(I,J-2)-4.0*U(I,J-1)+6.0*U(I,J)-4.0*U(I,J+1)+U(I,J+2)
     1      )/DLY4
         RETURN 
      ELSE IF (J == L - 1) THEN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D-DLY
C
         IF (KSWY /= 1) THEN
            UYYY = (U(I,L-4)-6.0*U(I,L-3)+12.0*U(I,L-2)-10.0*U(I,L-1)+
     1         3.0*U(I,L))/TDLY3
            UYYYY = ((-U(I,L-5))+6.0*U(I,L-4)-14.0*U(I,L-3)+16.0*U(I,L-2
     1         )-9.0*U(I,L-1)+2.0*U(I,L))/DLY4
            RETURN 
         ELSE
C
C     PERIODIC AT Y=D-DLY
C
            UYYY = ((-U(I,L-3))+2.0*U(I,L-2)-2.0*U(I,1)+U(I,2))/TDLY3
            UYYYY = (U(I,L-3)-4.0*U(I,L-2)+6.0*U(I,L-1)-4.0*U(I,1)+U(I,2
     1         ))/DLY4
            RETURN 
         ENDIF
      ELSE IF (J == L) THEN
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D
C
         UYYY = -(3.0*U(I,L-4)-14.0*U(I,L-3)+24.0*U(I,L-2)-18.0*U(I,L-1)
     1      +5.0*U(I,L))/TDLY3
         UYYYY = ((-2.0*U(I,L-5))+11.0*U(I,L-4)-24.0*U(I,L-3)+26.0*U(I,L
     1      -2)-14.0*U(I,L-1)+3.0*U(I,L))/DLY4
         RETURN 
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
c June      2004    version 5.0, fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE SEPDY
