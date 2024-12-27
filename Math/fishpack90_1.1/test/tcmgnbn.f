C
C     file tcmgnbn.f
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
      PROGRAM TCMBNGN
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IDIMF, M, MP1, MPEROD, N, NPEROD, I, J, IERROR
      REAL , DIMENSION(21) :: X
      REAL , DIMENSION(41) :: Y
      REAL :: DX, PI, DUM, DY, S, T, TSQ, T4, ERR
      COMPLEX , DIMENSION(22,40) :: F
      COMPLEX, DIMENSION(20) :: A, B, C
C-----------------------------------------------
c
C     PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE CMGNBN TO SOLVE
C     THE EQUATION
C
C     (1+X)**2*(D/DX)(DU/DX) - 2(1+X)(DU/DX) + (D/DY)(DU/DY)
C
C             - SQRT(-1)*U = (3 - SQRT(-1))*(1+X)**4*SIN(Y)         (1)
C
C     ON THE RECTANGLE 0 .LT. X .LT. 1 AND -PI .LT. Y .LT. PI
C     WITH THE BOUNDARY CONDITIONS
C
C     (DU/DX)(0,Y) = 4SIN(Y)                               (2)
C                                -PI .LE. Y .LE. PI
C     U(1,Y) = 16SIN(Y)                                    (3)
C
C     AND WITH U PERIODIC IN Y USING FINITE DIFFERENCES ON A
C     GRID WITH DELTAX (= DX) = 1/20 AND DELTAY (= DY) = PI/20.
C        TO SET UP THE FINITE DIFFERENCE EQUATIONS WE DEFINE
C     THE GRID POINTS
C
C     X(I) = (I-1)DX            I=1,2,...,21
C
C     Y(J) = -PI + (J-1)DY      J=1,2,...,41
C
C     AND LET V(I,J) BE AN APPROXIMATION TO U(X(I),Y(J)).
C     NUMBERING THE GRID POINTS IN THIS FASHION GIVES THE SET
C     OF UNKNOWNS AS V(I,J) FOR I=1,2,...,20 AND J=1,2,...,40.
C     HENCE, IN THE PROGRAM M = 20 AND N = 40.  AT THE INTERIOR
C     GRID POINT (X(I),Y(J)), WE REPLACE ALL DERIVATIVES IN
C     EQUATION (1) BY SECOND ORDER CENTRAL FINITE DIFFERENCES,
C     MULTIPLY BY DY**2, AND COLLECT COEFFICIENTS OF V(I,J) TO
C     GET THE FINITE DIFFERENCE EQUATION
C
C     A(I)V(I-1,J) + B(I)V(I,J) + C(I)V(I+1,J)
C
C     + V(I,J-1) - 2V(I,J) + V(I,J+1) = F(I,J)            (4)
C
C     WHERE S = (DY/DX)**2, AND FOR I=2,3,...,19
C
C     A(I) = (1+X(I))**2*S + (1+X(I))*S*DX
C
C     B(I) = -2(1+X(I))**2*S - SQRT(-1)*DY**2
C
C     C(I) = (1+X(I))**2*S - (1+X(I))*S*DX
C
C     F(I,J) = (3 - SQRT(-1))*(1+X(I))**4*DY**2*SIN(Y(J))
C              FOR J=1,2,...,40.
C
C        TO OBTAIN EQUATIONS FOR I = 1, WE REPLACE THE
C     DERIVATIVE IN EQUATION (2) BY A SECOND ORDER CENTRAL
C     FINITE DIFFERENCE APPROXIMATION, USE THIS EQUATION TO
C     ELIMINATE THE VIRTUAL UNKNOWN V(0,J) IN EQUATION (4)
C     AND ARRIVE AT THE EQUATION
C
C     B(1)V(1,J) + C(1)V(2,J) + V(1,J-1) - 2V(1,J) + V(1,J+1)
C
C                       = F(1,J)
C
C     WHERE
C
C     B(1) = -2S - SQRT(-1)*DY**2 , C(1) = 2S
C
C     F(1,J) = (11-SQRT(-1)+8/DX)*DY**2*SIN(Y(J)),  J=1,2,...,40.
C
C     FOR COMPLETENESS, WE SET A(1) = 0.
C        TO OBTAIN EQUATIONS FOR I = 20, WE INCORPORATE
C     EQUATION (3) INTO EQUATION (4) BY SETTING
C
C     V(21,J) = 16SIN(Y(J))
C
C     AND ARRIVE AT THE EQUATION
C
C     A(20)V(19,J) + B(20)V(20,J)
C
C     + V(20,J-1) - 2V(20,J) + V(20,J+1) = F(20,J)
C
C     WHERE
C
C     A(20) = (1+X(20))**2*S + (1+X(20))*S*DX
C
C     B(20) = -2*(1+X(20))**2*S - SQRT(-1)*DY**2
C
C     F(20,J) = ((3-SQRT(-1))*(1+X(20))**4*DY**2 - 16(1+X(20))**2*S
C                + 16(1+X(20))*S*DX)*SIN(Y(J))
C
C                    FOR J=1,2,...,40.
C
C     FOR COMPLETENESS, WE SET C(20) = 0.  HENCE, IN THE
C     PROGRAM MPEROD = 1.
C        THE PERIODICITY CONDITION ON U GIVES THE CONDITIONS
C
C     V(I,0) = V(I,40) AND V(I,41) = V(I,1) FOR I=1,2,...,20.
C
C     HENCE, IN THE PROGRAM NPEROD = 0.
C
C          THE EXACT SOLUTION TO THIS PROBLEM IS
C
C                  U(X,Y) = (1+X)**4*SIN(Y) .
C
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
      PI = 4.0*atan(1.0)
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
         A(I) = CMPLX((TSQ + T*DX)*S,0.)
         B(I) = (-2.*TSQ*S) - (0.,1.)*DY**2
         C(I) = CMPLX((TSQ - T*DX)*S,0.)
      END DO
      A(1) = (0.,0.)
      B(1) = (-2.*S) - (0.,1.)*DY**2
      C(1) = CMPLX(2.*S,0.)
      B(20) = (-2.*S*(1. + X(20))**2) - (0.,1.)*DY**2
      A(20) = CMPLX(S*(1. + X(20))**2+(1.+X(20))*DX*S,0.)
      C(20) = (0.,0.)
C
C     GENERATE RIGHT SIDE.
C
      DO I = 2, 19
         DO J = 1, N
            F(I,J) = (3.,-1.)*(1. + X(I))**4*DY**2*SIN(Y(J))
         END DO
      END DO
      T = 1. + X(20)
      TSQ = T**2
      T4 = TSQ**2
      DO J = 1, N
         F(1,J) = ((11.,-1.) + 8./DX)*DY**2*SIN(Y(J))
         F(20,J)=((3.,-1.)*T4*DY**2-16.*TSQ*S+16.*T*S*DX)*SIN(Y(J))
      END DO
      CALL CMGNBN (NPEROD, N, MPEROD, M, A, B, C, IDIMF, F, IERROR)
C
C     COMPUTE DISCRETIAZATION ERROR.  THE EXACT SOLUTION IS
C
C            U(X,Y) = (1+X)**4*SIN(Y) .
C
      ERR = 0.
      DO I = 1, M
         DO J = 1, N
            T = CABS(F(I,J)-(1.+X(I))**4*SIN(Y(J)))
            ERR = AMAX1(T,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    CMGNBN TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 9.1620E-3'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 9.1801E-3'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
      STOP 
 
 
 
      END PROGRAM TCMBNGN
