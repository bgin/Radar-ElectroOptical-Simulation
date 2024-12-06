C
C     file tpois3d.f
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
      PROGRAM TPOIS3D
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::LDIMF,MDIMF,LPEROD,L,MPEROD,M,NPEROD,N,I,J,K,IERROR
      REAL , DIMENSION(32,33,10) :: F
      REAL , DIMENSION(10) :: A, B, C
      REAL , DIMENSION(30) :: X, Y
      REAL , DIMENSION(10) :: Z
      REAL :: PI, DX, C1, DY, C2, DZ, DZSQ, T, ERR
C-----------------------------------------------
C
C     FROM THE DIMENSION STATEMENT WE GET THAT LDIMF = 32, MDIMF = 33,
C
      LDIMF = 32
      MDIMF = 33
      PI = 4.0*ATAN(1.0)
      LPEROD = 0
      L = 30
      DX = 2.*PI/FLOAT(L)
      C1 = 1./DX**2
      MPEROD = 0
      M = 30
      DY = 2.*PI/FLOAT(M)
      C2 = 1./DY**2
      NPEROD = 1
      N = 10
      DZ = 1./FLOAT(N)
      DZSQ = 1./DZ**2
C
C     GENERATE GRID POINTS FOR LATER USE.
C
      DO I = 1, L
         X(I) = (-PI) + FLOAT(I - 1)*DX
      END DO
      DO J = 1, M
         Y(J) = (-PI) + FLOAT(J - 1)*DY
      END DO
C
C     GENERATE COEFFICIENTS
C
      A(1) = 0.
      B(1) = -2.*DZSQ
      C(1) = -B(1)
      Z(1) = 0.
      DO K = 2, N
         Z(K) = FLOAT(K - 1)*DZ
         T = 1. + Z(K)
         A(K) = T**2*DZSQ + T/DZ
         B(K) = -2.*T**2*DZSQ
         C(K) = T**2*DZSQ - T/DZ
      END DO
C
C     GENERATE RIGHT SIDE OF EQUATION
C
      DO I = 1, L
         DO J = 1, M
            DO K = 2, N
               F(I,J,K) = 2.*SIN(X(I))*SIN(Y(J))*(1. + Z(K))**4
            END DO
         END DO
      END DO
      DO I = 1, L
         DO J = 1, L
            F(I,J,1) = (10. + 8./DZ)*SIN(X(I))*SIN(Y(J))
            F(I,J,N) = F(I,J,N) - C(N)*16.*SIN(X(I))*SIN(Y(J))
         END DO
      END DO
      C(N) = 0.
C
C     CALL POIS3D TO SOLVE EQUATIONS.
C
      CALL POIS3D (LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, C, 
     1   LDIMF, MDIMF, F, IERROR)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C
C              U(X,Y,Z) = SIN(X)*SIN(Y)*(1+Z)**4
C
      ERR = 0.
      DO I = 1, L
         DO J = 1, M
            DO K = 1, N
               T = ABS(F(I,J,K)-SIN(X(I))*SIN(Y(J))*(1.+Z(K))**4)
               ERR = AMAX1(T,ERR)
            END DO
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    POIS3D TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 2.93277E-2'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 2.93390E-2'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
      STOP 
      END PROGRAM TPOIS3D
