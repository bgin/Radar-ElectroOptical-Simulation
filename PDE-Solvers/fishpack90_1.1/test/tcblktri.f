C
C     file tcblktri.f
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
      PROGRAM TCBLKTRI
      USE fish
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IFLG, NP, N, MP, M, IDIMY, I, J, IERROR
      REAL , DIMENSION(105) :: AN, BN, CN
      REAL , DIMENSION(75) :: S
      REAL , DIMENSION(105) :: T
      REAL::DELTAS,DELTAT,HDS,TDS,TEMP1,TEMP2,TEMP3,HDT,TDT,ERR,Z
      COMPLEX , DIMENSION(75,105) :: Y
      COMPLEX, DIMENSION(75) :: AM, BM, CM
C-----------------------------------------------
C
      IFLG = 0
      NP = 1
      N = 63
      MP = 1
      M = 50
      IDIMY = 75
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
C     COEFFICIENTS AND THE ARRAY Y.
C
      DELTAS = 1./FLOAT(M + 1)
      DO I = 1, M
         S(I) = FLOAT(I)*DELTAS
      END DO
      DELTAT = 1./FLOAT(N + 1)
      DO J = 1, N
         T(J) = FLOAT(J)*DELTAT
      END DO
C
C     COMPUTE THE COEFFICIENTS AM,BM,CM CORRESPONDING TO THE S DIRECTION
C
      HDS = DELTAS/2.
      TDS = DELTAS + DELTAS
      DO I = 1, M
         TEMP1 = 1./(S(I)*TDS)
         TEMP2 = 1./((S(I)-HDS)*TDS)
         TEMP3 = 1./((S(I)+HDS)*TDS)
         AM(I) = CMPLX(TEMP1*TEMP2,0.)
         CM(I) = CMPLX(TEMP1*TEMP3,0.)
         BM(I) = (-(AM(I)+CM(I))) - (0.,1.)
      END DO
C
C     COMPUTE THE COEFFICIENTS AN,BN,CN CORRESPONDING TO THE T DIRECTION
C
      HDT = DELTAT/2.
      TDT = DELTAT + DELTAT
      DO J = 1, N
         TEMP1 = 1./(T(J)*TDT)
         TEMP2 = 1./((T(J)-HDT)*TDT)
         TEMP3 = 1./((T(J)+HDT)*TDT)
         AN(J) = TEMP1*TEMP2
         CN(J) = TEMP1*TEMP3
         BN(J) = -(AN(J)+CN(J))
      END DO
C
C     COMPUTE RIGHT SIDE OF EQUATION
C
      DO J = 1, N
         Y(:M,J) = 3.75*S(:M)*T(J)*(S(:M)**4+T(J)**4) - (0.,1.)*(S(:M)*T
     1      (J))**5
      END DO
C
C     THE NONZERO BOUNDARY CONDITIONS ENTER THE LINEAR SYSTEM VIA
C     THE RIGHT SIDE Y(I,J). IF THE EQUATIONS (3) GIVEN ABOVE
C     ARE EVALUATED AT I=M AND J=1,...,N THEN THE TERM CM(M)*U(M+1,J)
C     IS KNOWN FROM THE BOUNDARY CONDITION TO BE CM(M)*T(J)**5.
C     THEREFORE THIS TERM CAN BE INCLUDED IN THE RIGHT SIDE Y(M,J).
C     THE SAME ANALYSIS APPLIES AT J=N AND I=1,..,M. NOTE THAT THE
C     CORNER AT J=N,I=M INCLUDES CONTRIBUTIONS FROM BOTH BOUNDARIES.
C
      Y(M,:N) = Y(M,:N) - CM(M)*T(:N)**5
      Y(:M,N) = Y(:M,N) - CN(N)*S(:M)**5
      CALL CBLKTRI(IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,IERROR,W)
      IFLG = IFLG + 1
      DO WHILE(IFLG - 1 <= 0)
         CALL CBLKTRI (IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, IDIMY
     1      , Y, IERROR, W)
         IFLG = IFLG + 1
      END DO
      ERR = 0.
      DO J = 1, N
         DO I = 1, M
            Z = CABS(Y(I,J)-(S(I)*T(J))**5)
            ERR = AMAX1(Z,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    CBLKTRI TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 1.6457E-05'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 1.2737E-02'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
!     release dynamically allocated work space
      CALL FISHFIN (W)
      STOP 
      END PROGRAM TCBLKTRI
