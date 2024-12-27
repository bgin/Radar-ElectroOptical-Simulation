C
C     file tsepx4.f
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
      PROGRAM TSEPX4
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::M,N,NX,NY,I,J,MBDCND,NBDCND,IDMN,IORDER,IERROR
      REAL , DIMENSION(33,33) :: USOL, GRHS
      REAL , DIMENSION(33) :: BDA, BDB
      REAL :: A, B, C, D, DLX, DLY, X, AF, BF, CF, Y, ALPHA, BETA, DUM, 
     1   PERTRB, ERR, ERR2, ERR4
C-----------------------------------------------
      EXTERNAL COFX4
C
C     DEFINE ARITHMETIC FUNCTIONS GIVING EXACT SOLUTION
C
C
C     SET LIMITS ON REGION
C
      A = 0.0
      B = 1.0
      C = 0.0
      D = 1.0
C
C     SET GRID SIZE
C
      M = 32
      N = 32
      DLX = (B - A)/FLOAT(M)
      DLY = (D - C)/FLOAT(N)
      NX = M + 1
      NY = N + 1
      DO I = 1, NX
         X = A + FLOAT(I - 1)*DLX
C
C     SET SPECIFIED BOUNDARY CONDITIONS AT Y=C,D
C
         USOL(I,1) = UE(X,C)
         USOL(I,NY) = UE(X,D)
         CALL COFX4 (X, AF, BF, CF)
         DO J = 1, NY
            Y = C + FLOAT(J - 1)*DLY
C
C     SET RIGHT HAND SIDE
C
            GRHS(I,J)=AF*UXXE(X,Y)+BF*UXE(X,Y)+CF*UE(X,Y)+UYYE(X,Y)
         END DO
      END DO
C
C     SET MIXED BOUNDARY CONDITIONS AT X=A,B
C
      ALPHA = 1.0
      BETA = 1.0
      DO J = 1, NY
         Y = C + FLOAT(J - 1)*DLY
         BDA(J) = UXE(A,Y) + ALPHA*UE(A,Y)
         BDB(J) = UXE(B,Y) + BETA*UE(B,Y)
      END DO
C
C     SET BOUNDARY SWITHCES
C
      MBDCND = 3
      NBDCND = 1
C
C     SET FIRST DIMENSION OF USOL,GRHS AND WORK SPACE LENGTH
C
      IDMN = 33
C
C     OBTAIN SECOND ORDER APPROXIMATION
C
      IORDER = 2
      CALL SEPX4 (IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA, C, D, 
     1   N, NBDCND, DUM, DUM, COFX4, GRHS, USOL, IDMN, PERTRB, IERROR)
C
C     COMPUTE SECOND ORDER DISCRETIZATION ERROR (RELATIVE)
C     ALSO RESET SPECIFIED BOUNDARIES AND RIGHT HAND SIDE.
C
      ERR = 0.0
      DO I = 1, NX
         X = A + FLOAT(I - 1)*DLX
         USOL(I,1) = UE(X,C)
         USOL(I,NY) = UE(X,D)
         CALL COFX4 (X, AF, BF, CF)
         DO J = 1, NY
            Y = C + FLOAT(J - 1)*DLY
            ERR = AMAX1(ERR,ABS((USOL(I,J)-UE(X,Y))/UE(X,Y)))
C
C     RESET RIGHT HAND SIDE IN GRHS FOR FOURTH ORDER APPROXIMATION CALL
C
            GRHS(I,J)=AF*UXXE(X,Y)+BF*UXE(X,Y)+CF*UE(X,Y)+UYYE(X,Y)
         END DO
      END DO
      ERR2 = ERR
C
C     OBTAIN FOURTH ORDER APPROXIMATION
C
      IORDER = 4
      CALL SEPX4 (IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA, C, D, 
     1   N, NBDCND, DUM, DUM, COFX4, GRHS, USOL, IDMN, PERTRB, IERROR)
C
C     COMPUTE FOURTH ORDER DISCRETIZATION ERROR (RELATIVE)
C
      ERR = 0.0
      DO J = 1, NY
         Y = C + FLOAT(J - 1)*DLY
         DO I = 1, NX
            X = A + FLOAT(I - 1)*DLX
            ERR = AMAX1(ERR,ABS((USOL(I,J)-UE(X,Y))/UE(X,Y)))
         END DO
      END DO
      ERR4 = ERR
      WRITE (*, *) '    SEPEX4 TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0'
      WRITE (*, *) '    Second Order Discretization Error = 1.5985E-4'
      WRITE (*, *) '    Fourth Order Discretization Error = 1.8575E-6'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0'
      WRITE (*, *) '    Second Order Discretization Error = 1.5044E-4'
      WRITE (*, *) '    Fourth Order Discretization Error = 1.5736E-5'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR
      WRITE (*, *) '    Second Order Discretization Error =', ERR2
      WRITE (*, *) '    Fourth Order Discretization Error =', ERR4
      STOP 
      CONTAINS


      REAL FUNCTION UE (S, T)
      REAL, INTENT(IN) :: S
      REAL, INTENT(IN) :: T
      UE = (S*T)**3 + 1.0
      RETURN 
      END FUNCTION UE


      REAL FUNCTION UXE (S, T)
      REAL, INTENT(IN) :: S
      REAL, INTENT(IN) :: T
      UXE = 3.0*S**2*T**3
      RETURN 
      END FUNCTION UXE


      REAL FUNCTION UXXE (S, T)
      REAL, INTENT(IN) :: S
      REAL, INTENT(IN) :: T
      UXXE = 6.0*S*T**3
      RETURN 
      END FUNCTION UXXE


      REAL FUNCTION UYE (S, T)
      REAL, INTENT(IN) :: S
      REAL, INTENT(IN) :: T
      UYE = 3.0*S**3*T**2
      RETURN 
      END FUNCTION UYE


      REAL FUNCTION UYYE (S, T)
      REAL, INTENT(IN) :: S
      REAL, INTENT(IN) :: T
      UYYE = 6.0*S**3*T
      RETURN 
      END FUNCTION UYYE
      END PROGRAM TSEPX4


      SUBROUTINE COFX4(X, AF, BF, CF)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL , INTENT(IN) :: X
      REAL , INTENT(OUT) :: AF
      REAL , INTENT(OUT) :: BF
      REAL , INTENT(OUT) :: CF
C-----------------------------------------------
C
C     SET COEFFICIENTS IN THE X-DIRECTION.
C
      AF = (X + 1.)**2
      BF = 2.0*(X + 1.)
      CF = -X
      RETURN 
      END SUBROUTINE COFX4
