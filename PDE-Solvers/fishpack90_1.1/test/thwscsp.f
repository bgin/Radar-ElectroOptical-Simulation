C
C     file thwscsp.f
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
      PROGRAM THWSCSP
      USE fish
      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::INTL,M,MBDCND,N,NBDCND,IDIMF,MP1,I,NP1,J,IERROR
      REAL , DIMENSION(48,33) :: F
      REAL , DIMENSION(33) :: BDTF
      REAL , DIMENSION(48) :: THETA
      REAL , DIMENSION(33) :: R
      REAL :: PI, DUM, TS, TF, RS, RF, ELMBDA, DTHETA, DR, CI4, BDTS, 
     1   BDRS, BDRF, PERTRB, ERR, Z, DPHI, SI
C-----------------------------------------------
C
C     PROGRAM TO ILLUSTRATE THE USE OF HWSCSP
C
C
      PI = 4.0*ATAN(1.0)
      INTL = 0
      TS = 0.
      TF = PI/2.
      M = 36
      MBDCND = 6
      RS = 0.
      RF = 1.
      N = 32
      NBDCND = 5
      ELMBDA = 0.
      IDIMF = 48
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE EQUATION.
C
      MP1 = M + 1
      DTHETA = TF/FLOAT(M)
      DO I = 1, MP1
         THETA(I) = FLOAT(I - 1)*DTHETA
      END DO
      NP1 = N + 1
      DR = 1./FLOAT(N)
      DO J = 1, NP1
         R(J) = FLOAT(J - 1)*DR
      END DO
C
C     GENERATE NORMAL DERIVATIVE DATA AT EQUATOR
C
      BDTF(:NP1) = 0.
C
C     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
C
      DO I = 1, MP1
         F(I,N+1) = COS(THETA(I))**4
      END DO
C
C     COMPUTE RIGHT SIDE OF EQUATION
C
      DO I = 1, MP1
         CI4 = 12.*COS(THETA(I))**2
         F(I,:N) = CI4*R(:N)**2
      END DO
      CALL HWSCSP (INTL, TS, TF, M, MBDCND, BDTS, BDTF, RS, RF, N, 
     1   NBDCND, BDRS, BDRF, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C
C     COMPUTE DISCRETIZATION ERROR
C
      ERR = 0.
      DO I = 1, MP1
         CI4 = COS(THETA(I))**4
         DO J = 1, N
            Z = ABS(F(I,J)-CI4*R(J)**4)
            ERR = AMAX1(Z,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) ' HWSCSP TEST RUN, EXAMPLE 1 *** '
      WRITE (*, *) ' Previous 64 bit floating point arithmetic result '
      WRITE (*, *) ' ierror = 0'
      WRITE (*, *) ' discretization error = 7,9984E-4 '
      WRITE (*, *) ' Previous 32 bit floating point arithmetic result '
      WRITE (*, *) ' ierror = 0'
      WRITE (*, *) ' discretization error = 7.9907E-4 '
      WRITE (*, *) ' The output from your computer is: '
      WRITE (*, *) ' IERROR =', IERROR
      WRITE (*, *) ' Discretization Error =', ERR
C
C     THE FOLLOWING PROGRAM ILLUSTRATES THE USE OF HWSCSP TO SOLVE
C     A THREE DIMENSIONAL PROBLEM WHICH HAS LONGITUDNAL DEPENDENCE
C
      MBDCND = 2
      NBDCND = 1
      DPHI = PI/72.
      ELMBDA = -2.*(1. - COS(DPHI))/DPHI**2
C
C     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
C
      DO I = 1, MP1
         F(I,N+1) = SIN(THETA(I))
      END DO
C
C     COMPUTE RIGHT SIDE OF THE EQUATION
C
      F(:MP1,:N) = 0.
      CALL HWSCSP (INTL, TS, TF, M, MBDCND, BDTS, BDTF, RS, RF, N, 
     1   NBDCND, BDRS, BDRF, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
C
C     COMPUTE DISCRETIZATION ERROR   (FOURIER COEFFICIENTS)
C
      ERR = 0
      DO I = 1, MP1
         SI = SIN(THETA(I))
         DO J = 1, NP1
            Z = ABS(F(I,J)-R(J)*SI)
            ERR = AMAX1(Z,ERR)
         END DO
      END DO
C
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
 
      WRITE (*, *) ' ********** '
      WRITE (*, *) ' ********** '
      WRITE (*, *) ' HWSCSP TEST RUN, EXAMPLE 2 *** '
      WRITE (*, *) ' Previous 64 bit floating point arithmetic result '
      WRITE (*, *) ' ierror = 0'
      WRITE (*, *) ' discretization error = 5.8682E-5 '
      WRITE (*, *) ' Previous 32 bit floating point arithmetic result '
      WRITE (*, *) ' ierror = 0'
      WRITE (*, *) ' discretization error = 5.9962E-5 '
      WRITE (*, *) ' The output from your computer is: '
      WRITE (*, *) ' IERROR =', IERROR
      WRITE (*, *) ' Discretization Error =', ERR
!     release real and complex allocated work space
      CALL FISHFIN (W)
      STOP 
      END PROGRAM THWSCSP
