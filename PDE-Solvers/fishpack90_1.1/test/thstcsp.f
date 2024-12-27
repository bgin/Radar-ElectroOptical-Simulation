C
C     file thstcsp.f
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
      PROGRAM THSTCSP
      USE fish
      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IDIMF, M, MBDCND, I, N, NBDCND, J, INTL, IERROR
      REAL , DIMENSION(47,16) :: F
      REAL , DIMENSION(45) :: BDD, THETA
      REAL , DIMENSION(15) :: R
      REAL , DIMENSION(45) :: COST
      REAL :: A, B, DT, C, D, DR, ELMBDA, BDA, BDB, BDC, PERTRB, ERR, Z
C-----------------------------------------------
C
C     NOTE THAT FROM DIMENSION STATEMENT WE GET THAT IDIMF = 47
C
      IDIMF = 47
      A = 0.
      B = 4.0*ATAN(1.0)
C
C     NOTE THAT B IS SET TO PI USING THE FUNCTION PIMACH AS REQUIRED.
C
      M = 45
      MBDCND = 9
      DT = (B - A)/FLOAT(M)
C
C     DEFINE GRID POINTS THETA(I) AND COS(THETA(I))
C
      DO I = 1, M
         THETA(I) = A + (FLOAT(I) - 0.5)*DT
         COST(I) = COS(THETA(I))
      END DO
      C = 0.
      D = 1.
      N = 15
      NBDCND = 5
      DR = (D - C)/FLOAT(N)
C
C     DEFINE GRID POINTS R(J)
C
      DO J = 1, N
         R(J) = C + (FLOAT(J) - 0.5)*DR
      END DO
C
C     DEFINE BOUNDARY ARRAY BDD.  BDA, BDB, AND BDC ARE DUMMY
C     VARIABLES IN THIS EXAMPLE.
C
      BDD(:M) = COST(:M)**4
      ELMBDA = 0.
C
C     DEFINE RIGHT SIDE F
C
      DO I = 1, M
         F(I,:N) = 12.*(R(:N)*COST(I))**2
      END DO
      INTL = 0
      CALL HSTCSP (INTL, A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC
     1   , BDD, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C
C     U(THETA,R) = (R*COS(THETA))**4
C
      ERR = 0.
      DO I = 1, M
         DO J = 1, N
            Z = ABS(F(I,J)-(R(J)*COST(I))**4)
            ERR = AMAX1(Z,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
!     in this example (contrast with blktri and sepeli) the extra precision
!     does not reduce the discretization error
      WRITE (*, *) '    HSTCSP TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 5.5843E-3'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 5.5845E-3'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
!     release work space allocated by hstcsp (intl=0 call)
      CALL FISHFIN (W)
      STOP 
      END PROGRAM THSTCSP
