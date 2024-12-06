C
C     file tgenbun.f
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
      PROGRAM TGENBUN
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IDIMF, M, MP1, MPEROD, N, NPEROD, I, J, IERROR
      REAL , DIMENSION(22,40) :: F
      REAL , DIMENSION(20) :: A, B, C
      REAL , DIMENSION(21) :: X
      REAL , DIMENSION(41) :: Y
      REAL :: DX, PI, DY, S, T, TSQ, T4, ERR
C-----------------------------------------------
C
C     FROM THE DIMENSION STATEMENT WE GET THAT IDIMF = 22
C
      IDIMF = 22
      M = 20
      MP1 = M + 1
      MPEROD = 1
      DX = 0.05
      N = 40
      NPEROD = 0
      PI = 4.0*ATAN(1.0)
      DY = PI/20.
C
C     GENERATE GRID POINTS FOR LATER USE.
C
      DO I = 1, MP1
         X(I) = FLOAT(I - 1)*DX
      END DO
      DO J = 1, N
         Y(J) = (-PI) + FLOAT(J - 1)*DY
      END DO
C
C     GENERATE COEFFICIENTS.
C
      S = (DY/DX)**2
      DO I = 2, 19
         T = 1. + X(I)
         TSQ = T**2
         A(I) = (TSQ + T*DX)*S
         B(I) = -2.*TSQ*S
         C(I) = (TSQ - T*DX)*S
      END DO
      A(1) = 0.
      B(1) = -2.*S
      C(1) = -B(1)
      B(20) = -2.*S*(1. + X(20))**2
      A(20) = (-B(20)/2.) + (1. + X(20))*DX*S
      C(20) = 0.
C
C     GENERATE RIGHT SIDE.
C
      DO I = 2, 19
         DO J = 1, N
            F(I,J) = 3.*(1. + X(I))**4*DY**2*SIN(Y(J))
         END DO
      END DO
      T = 1. + X(20)
      TSQ = T**2
      T4 = TSQ**2
      DO J = 1, N
         F(1,J) = (11. + 8./DX)*DY**2*SIN(Y(J))
         F(20,J) = (3.*T4*DY**2 - 16.*TSQ*S + 16.*T*S*DX)*SIN(Y(J))
      END DO
      CALL GENBUN (NPEROD, N, MPEROD, M, A, B, C, IDIMF, F, IERROR)
C
C     COMPUTE DISCRETIAZATION ERROR.  THE EXACT SOLUTION IS
C
C            U(X,Y) = (1+X)**4*SIN(Y) .
C
      ERR = 0.
      DO I = 1, M
         DO J = 1, N
            T = ABS(F(I,J)-(1.+X(I))**4*SIN(Y(J)))
            ERR = AMAX1(T,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    GENBUN TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 9.6406E-3'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 9.6556E-3'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
      STOP 
      END PROGRAM TGENBUN
