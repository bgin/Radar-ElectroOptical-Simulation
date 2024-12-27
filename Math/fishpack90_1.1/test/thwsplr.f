C
C     file thwsplr.f
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
      PROGRAM THWSPLR
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IDIMF, M, MBDCND, N, NBDCND, MP1, NP1, I, J, IERROR
      REAL , DIMENSION(100,50) :: F
      REAL , DIMENSION(51) :: BDC, BDD, R
      REAL , DIMENSION(49) :: THETA
      REAL :: A, B, C, PI, DUM, D, ELMBDA, BDA, BDB, PERTRB, W, ERR, Z
C-----------------------------------------------
C
C          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSPLR TO SOLVE
C     THE EQUATION
C
C     (1/R)(D/DR)(R*(DU/DR)) + (1/R**2)(D/DTHETA)(DU/DTHETA) = 16*R**2
C
C     ON THE QUARTER-DISK 0 .LT. R .LT. 1, 0 .LT. THETA .LT. PI/2 WITH
C     WITH THE BOUNDARY CONDITIONS
C
C     U(1,THETA) = 1 - COS(4*THETA), 0 .LE. THETA .LE. 1
C
C     AND
C
C     (DU/DTHETA)(R,0) = (DU/DTHETA)(R,PI/2) = 0,  0 .LE. R .LE. 1.
C
C     (NOTE THAT THE SOLUTION U IS UNSPECIFIED AT R = 0.)
C          THE R-INTERVAL WILL BE DIVIDED INTO 50 PANELS AND THE
C     THETA-INTERVAL WILL BE DIVIDED INTO 48 PANELS.
C
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.
C
      IDIMF = 100
      A = 0.
      B = 1.
      M = 50
      MBDCND = 5
      C = 0.
      PI = 4.0*ATAN(1.0)
      D = PI/2.
      N = 48
      NBDCND = 3
      ELMBDA = 0.
C
C     AUXILIARY QUANTITIES.
C
      MP1 = M + 1
      NP1 = N + 1
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
C
      DO I = 1, MP1
         R(I) = FLOAT(I - 1)/50.
      END DO
      DO J = 1, NP1
         THETA(J) = FLOAT(J - 1)*PI/96.
      END DO
C
C     GENERATE BOUNDARY DATA.
C
      BDC(:MP1) = 0.
      BDD(:MP1) = 0.
C
C     BDA AND BDB ARE DUMMY VARIABLES.
C
      DO J = 1, NP1
         F(MP1,J) = 1. - COS(4.*THETA(J))
      END DO
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO I = 1, M
         F(I,:NP1) = 16.*R(I)**2
      END DO
      CALL HWSPLR (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD
     1   , ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C                U(R,THETA) = R**4*(1 - COS(4*THETA))
C
      ERR = 0.
      DO I = 1, MP1
         DO J = 1, NP1
            Z = ABS(F(I,J)-R(I)**4*(1.-COS(4.*THETA(J))))
            ERR = AMAX1(Z,ERR)
         END DO
      END DO
      WRITE (*, *) '    HWSPLR TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 6.19134E-4'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 6.20723E-4'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
      STOP 
      END PROGRAM THWSPLR
