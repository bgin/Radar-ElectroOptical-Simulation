C
C     file thstcrt.f
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
      PROGRAM THSTCRT
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IDIMF, M, MBDCND, N, NBDCND, I, J, IERROR
      REAL , DIMENSION(50,53) :: F
      REAL , DIMENSION(53) :: BDA, BDB
      REAL , DIMENSION(48) :: X
      REAL , DIMENSION(53) :: Y
      REAL::A,B,DX,C,D,DY,ELMBDA,PI,PISQ,T,BDC,BDD,PERTRB,ERR
C-----------------------------------------------
C
C     FROM THE DIMENSION STATEMENT WE GET IDIMF = 50.
C
      IDIMF = 50
      A = 1.
      B = 3.
      M = 48
      DX = (B - A)/FLOAT(M)
      MBDCND = 2
      C = -1.
      D = 1.
      N = 53
      DY = (D - C)/FLOAT(N)
      NBDCND = 0
      ELMBDA = -2.
C
C     AUXILIARY QUANTITIES
C
      PI = 4.0*ATAN(1.0)
      PISQ = PI*PI
C
C     GENERATE AND STORE GRID POINTS FOR COMPUTATION OF BOUNDARY DATA
C     AND THE RIGHT SIDE OF THE HELMHOLTZ EQUATION.
C
      DO I = 1, M
         X(I) = A + (FLOAT(I) - 0.5)*DX
      END DO
      DO J = 1, N
         Y(J) = C + (FLOAT(J) - 0.5)*DY
      END DO
C
C     GENERATE BOUNDARY DATA.
C
      DO J = 1, N
         BDA(J) = 0.
         BDB(J) = -PI*COS(PI*Y(J))
      END DO
C
C     BDC AND BDD ARE DUMMY ARGUMENTS IN THIS EXAMPLE.
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      T = -2.*(PISQ + 1.)
      DO I = 1, M
         DO J = 1, N
            F(I,J) = T*SIN(PI*X(I))*COS(PI*Y(J))
         END DO
      END DO
      CALL HSTCRT (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD
     1   , ELMBDA, F, IDIMF, PERTRB, IERROR)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C
C               U(X,Y) = SIN(PI*X)*COS(PI*Y) .
C
      ERR = 0.
      DO I = 1, M
         DO J = 1, N
            T = ABS(F(I,J)-SIN(PI*X(I))*COS(PI*Y(J)))
            ERR = AMAX1(T,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    HSTCRT TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 1.2600E-3'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 1.2586E-3'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
      STOP 
      END PROGRAM THSTCRT
