C
C     file thstplr.f
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
      PROGRAM THSTPLR
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IDIMF, M, MBDCND, N, NBDCND, I, J, IERROR
      REAL , DIMENSION(51,50) :: F
      REAL , DIMENSION(48) :: BDB
      REAL , DIMENSION(50) :: BDC, BDD, R
      REAL , DIMENSION(48) :: THETA
      REAL :: A, B, C, PI, D, ELMBDA, BDA, PERTRB, ERR, Z
C-----------------------------------------------
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.
C
      IDIMF = 51
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
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
C
      DO I = 1, M
         R(I) = (FLOAT(I) - 0.5)/50.
      END DO
      DO J = 1, N
         THETA(J) = (FLOAT(J) - 0.5)*PI/96.
      END DO
C
C     GENERATE BOUNDARY DATA.
C
      DO J = 1, N
         BDB(J) = 1. - COS(4.*THETA(J))
      END DO
C
C     GENERATE BOUNDARY DATA.
C
      BDC(:M) = 0.
      BDD(:M) = 0.
C
C     BDA IS A DUMMY VARIABLE.
C
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO I = 1, M
         F(I,:N) = 16.*R(I)**2
      END DO
      CALL HSTPLR (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD
     1   , ELMBDA, F, IDIMF, PERTRB, IERROR)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C
C                U(R,THETA) = R**4*(1 - COS(4*THETA))
C
      ERR = 0.
      DO I = 1, M
         DO J = 1, N
            Z = ABS(F(I,J)-R(I)**4*(1.-COS(4.*THETA(J))))
            ERR = AMAX1(Z,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    HSTPLR TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 1.1303E-3'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 1.1300E-3'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
      STOP 
      END PROGRAM THSTPLR
