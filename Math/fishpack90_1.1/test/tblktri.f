C
C     file tblktri.f
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
C     PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE BLKTRI TO
C     SOLVE THE EQUATION
C
C     .5/S*(D/DS)(.5/S*DU/DS)+.5/T*(D/DT)(.5/T*DU/DT)
C                                                          (1)
C                   = 15/4*S*T*(S**4+T**4)
C
C     ON THE RECTANGLE 0 .LT. S .LT. 1 AND 0 .LT. T .LT. 1
C     WITH THE BOUNDARY CONDITIONS
C
C     U(0,T) = 0
C                            0 .LE. T .LE. 1
C     U(1,T) = T**5
C
C     AND
C
C     U(S,0) = 0
C                            0 .LE. S .LE. 1
C     U(S,1) = S**5
C
C     THE EXACT SOLUTION OF THIS PROBLEM IS U(S,T) = (S*T)**5
C
C     DEFINE THE INTEGERS M = 50 AND N = 63. THEN DEFINE THE
C     GRID INCREMENTS DELTAS = 1/(M+1) AND DELTAT = 1/(N+1).
C
C     THE GRID IS THEN GIVEN BY S(I) = I*DELTAS FOR I = 1,...,M
C     AND T(J) = J*DELTAT FOR J = 1,...,N.
C
C     THE APPROXIMATE SOLUTION IS GIVEN AS THE SOLUTION TO
C     THE FOLLOWING FINITE DIFFERENCE APPROXIMATION OF EQUATION (1).
C
C     .5/(S(I)*DELTAS)*((U(I+1,J)-U(I,J))/(2*S(I+.5)*DELTAS)
C                     -(U(I,J)-U(I-1,J))/(2*S(I-.5)*DELTAS))
C     +.5/(T(I)*DELTAT)*((U(I,J+1)-U(I,J))/(2*T(I+.5)*DELTAT) (2)
C                     -(U(I,J)-U(I,J-1))/(2*T(I-.5)*DELTAT))
C               = 15/4*S(I)*T(J)*(S(I)**4+T(J)**4)
C
C             WHERE S(I+.5) = .5*(S(I+1)+S(I))
C                   S(I-.5) = .5*(S(I)+S(I-1))
C                   T(I+.5) = .5*(T(I+1)+T(I))
C                   T(I-.5) = .5*(T(I)+T(I-1))
C
C     THE APPROACH IS TO WRITE EQUATION (2) IN THE FORM
C
C     AM(I)*U(I-1,J)+BM(I,J)*U(I,J)+CM(I)*U(I+1,J)
C       +AN(J)*U(I,J-1)+BN(J)*U(I,J)+CN(J)*U(I,J+1)      (3)
C           = Y(I,J)
C
C     AND THEN CALL SUBROUTINE BLKTRI TO DETERMINE U(I,J)
C
C
C
      PROGRAM TBLKTRI
      USE fish
      TYPE ( fishworkspace) :: w
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IFLG, NP, N, MP, M, IDIMY, I, J, IERROR
      REAL , DIMENSION(75,105) :: Y
      REAL , DIMENSION(75) :: AM, BM, CM
      REAL , DIMENSION(105) :: AN, BN, CN
      REAL , DIMENSION(75) :: S
      REAL , DIMENSION(105) :: T
      REAL::DELTAS,DELTAT,HDS,TDS,TEMP1,TEMP2,TEMP3,HDT,TDT,ERR,Z
C-----------------------------------------------
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
         AM(I) = TEMP1*TEMP2
         CM(I) = TEMP1*TEMP3
         BM(I) = -(AM(I)+CM(I))
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
         Y(:M,J) = 3.75*S(:M)*T(J)*(S(:M)**4+T(J)**4)
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
C
C     DETERMINE THE APPROXIMATE SOLUTION U(I,J)
C
      CALL BLKTRI(IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,IERROR,W)
      IFLG = IFLG + 1
      DO WHILE(IFLG - 1 <= 0)
         CALL BLKTRI (IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, IDIMY
     1      , Y, IERROR, W)
         IFLG = IFLG + 1
      END DO
      ERR = 0.
      DO J = 1, N
         DO I = 1, M
            Z = ABS(Y(I,J)-(S(I)*T(J))**5)
            ERR = AMAX1(Z,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    BLKTRI TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 1.6478E-05'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 1.2737E-02'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
!     release dynamically allocated work space
      CALL FISHFIN (W)
      STOP 
      END PROGRAM TBLKTRI
