C
C     file thwsssp.f
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
C     PROGRAM TO ILLUSTRATE THE USE OF HWSSSP
C
      PROGRAM THWSSSP
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: M, MBDCND, N, NBDCND, IDIMF, MP1, I, NP1, J, IERROR
      REAL , DIMENSION(19,73) :: F
      REAL , DIMENSION(73) :: BDTF
      REAL , DIMENSION(19) :: SINT
      REAL , DIMENSION(73) :: SINP
      REAL :: PI, TS, TF, PS, PF, ELMBDA, DTHETA, DPHI, BDTS, BDPS, BDPF
     1   , PERTRB, ERR, Z
C-----------------------------------------------
      PI = 4.0*ATAN(1.0)
      TS = 0
      TF = PI/2.
      M = 18
      MBDCND = 6
      PS = 0
      PF = PI + PI
      N = 72
      NBDCND = 0
      ELMBDA = 0.
      IDIMF = 19
C
C     GENERATE SINES FOR USE IN SUBSEQUENT COMPUTATIONS
C
      DTHETA = TF/FLOAT(M)
      MP1 = M + 1
      DO I = 1, MP1
         SINT(I) = SIN(FLOAT(I - 1)*DTHETA)
      END DO
      DPHI = (PI + PI)/FLOAT(N)
      NP1 = N + 1
      DO J = 1, NP1
         SINP(J) = SIN(FLOAT(J - 1)*DPHI)
      END DO
C
C     COMPUTE RIGHT SIDE OF EQUATION AND STORE IN F
C
      DO J = 1, NP1
         F(:MP1,J) = 2. - 6.*(SINT(:MP1)*SINP(J))**2
      END DO
C
C     STORE DERIVATIVE DATA AT THE EQUATOR
C
      BDTF(:NP1) = 0.
C
      CALL HWSSSP (TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, NBDCND, 
     1   BDPS, BDPF, ELMBDA, F, IDIMF, PERTRB, IERROR)
C
C     COMPUTE DISCRETIZATION ERROR. SINCE PROBLEM IS SINGULAR, THE
C     SOLUTION MUST BE NORMALIZED.
C
      ERR = 0
      DO J = 1, NP1
         DO I = 1, MP1
            Z = ABS(F(I,J)-(SINT(I)*SINP(J))**2-F(1,1))
            ERR = AMAX1(Z,ERR)
         END DO
      END DO
C
      WRITE (*, *) '    HWSSSP TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 3.38107E-3'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 3.3650E-3'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
      STOP 
      END PROGRAM THWSSSP
