C
C     file thstcyl.f
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
      PROGRAM THSTCYL
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IDIMF, M, MBDCND, N, NBDCND, I, J, IERROR
      REAL , DIMENSION(51,52) :: F
      REAL , DIMENSION(52) :: BDB
      REAL , DIMENSION(50) :: BDC, BDD, R
      REAL , DIMENSION(52) :: Z
      REAL :: A, B, C, D, ELMBDA, BDA, PERTRB, X, ERR
C-----------------------------------------------
C
C     PROGRAM TO ILLUSTRATE THE USE OF HSTCYL TO SOLVE THE EQUATION
C
C    (1/R)(D/DR)(R*DU/DR) + (D/DZ)(DU/DZ) = (2*R*Z)**2*(4*Z**2 + 3*R**2)
C
C     ON THE RECTANGLE 0 .LT. R .LT. 1 , 0 .LT. Z .LT. 1 WITH THE
C     BOUNDARY CONDITIONS
C
C     (DU/DR)(1,Z) = 4*Z**2  FOR  0 .LE. Z .LE. 1
C
C     AND
C
C     (DU/DZ)(R,0) = 0 AND (DU/DZ)(R,1) = 4*R**2  FOR  0 .LE. R .LE. 1 .
C
C     THE SOLUTION TO THIS PROBLEM IS NOT UNIQUE.  IT IS A
C     ONE-PARAMETER FAMILY OF SOLUTIONS GIVEN BY
C
C            U(R,Z) = (R*Z)**4 + ARBITRARY CONSTANT .
C
C     THE R-INTERVAL WILL CONTAIN 50 UNKNOWNS AND THE Z-INTERVAL WILL
C     CONTAIN 52 UNKNOWNS.
C
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.
C
      IDIMF = 51
      A = 0.
      B = 1.
      M = 50
      MBDCND = 6
      C = 0.
      D = 1.
      N = 52
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
         Z(J) = (FLOAT(J) - 0.5)/52.
      END DO
C
C     GENERATE BOUNDARY DATA.
C
      BDB(:N) = 4.*Z(:N)**4
C
C     GENERATE BOUNDARY DATA.
C
      BDC(:M) = 0.
      BDD(:M) = 4.*R(:M)**4
C
C     BDA IS A DUMMY VARIABLE.
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO I = 1, M
         F(I,:N) = 4.*R(I)**2*Z(:N)**2*(4.*Z(:N)**2+3.*R(I)**2)
      END DO
      CALL HSTCYL (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD
     1   , ELMBDA, F, IDIMF, PERTRB, IERROR)
C
C     COMPUTE DISCRETIZATION ERROR BY MINIMIZING OVER ALL A THE FUNCTION
C     NORM(F(I,J) - A*1 - U(R(I),Z(J))).  THE EXACT SOLUTION IS
C                U(R,Z) = (R*Z)**4 + ARBITRARY CONSTANT.
C
      X = 0.
      DO I = 1, M
         X = X + SUM(F(I,:N)-(R(I)*Z(:N))**4)
      END DO
      X = X/FLOAT(M*N)
      F(:M,:N) = F(:M,:N) - X
      ERR = 0.
      DO I = 1, M
         DO J = 1, N
            X = ABS(F(I,J)-(R(I)*Z(J))**4)
            ERR = AMAX1(X,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    HSTCYL TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  PERTRB =-4.4311E-4'
      WRITE (*, *) '    Discretization Error = 7.5280E-5 '
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  PERTRB =-4.4321E-4'
      WRITE (*, *) '    Discretization Error = 7.3557E-5'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' PERTRB = ', PERTRB
      WRITE (*, *) '    Discretization Error = ', ERR
      STOP 
      END PROGRAM THSTCYL
