C
C     file hw3crt.f
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
C     SUBROUTINE HW3CRT (XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,MBDCND,BDYS,
C    +                   BDYF,ZS,ZF,N,NBDCND,BDZS,BDZF,ELMBDA,LDIMF,
C    +                   MDIMF,F,PERTRB,IERROR)
C
C
C DIMENSION OF           BDXS(MDIMF,N+1),    BDXF(MDIMF,N+1),
C ARGUMENTS              BDYS(LDIMF,N+1),    BDYF(LDIMF,N+1),
C                        BDZS(LDIMF,M+1),    BDZF(LDIMF,M+1),
C                        F(LDIMF,MDIMF,N+1)
C
C LATEST REVISION        June 2004
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
C                        EQUATION IN CARTESIAN COORDINATES.  THIS
C                        EQUATION IS
C
C                          (D/DX)(DU/DX) + (D/DY)(DU/DY) +
C                          (D/DZ)(DU/DZ) + LAMBDA*U = F(X,Y,Z) .
C
C USAGE                  CALL HW3CRT (XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,
C                                     MBDCND,BDYS,BDYF,ZS,ZF,N,NBDCND,
C                                     BDZS,BDZF,ELMBDA,LDIMF,MDIMF,F,
C                                     PERTRB,IERROR)
C
C ARGUMENTS
C
C ON INPUT               XS,XF
C
C                          THE RANGE OF X, I.E. XS .LE. X .LE. XF .
C                          XS MUST BE LESS THAN XF.
C
C                        L
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (XS,XF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE L+1 GRID POINTS
C                          IN THE X-DIRECTION GIVEN BY
C                          X(I) = XS+(I-1)DX FOR I=1,2,...,L+1,
C                          WHERE DX = (XF-XS)/L IS THE PANEL WIDTH.
C                          L MUST BE AT LEAST 5.
C
C                        LBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT X = XS AND X = XF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN X,
C                               I.E. U(L+I,J,K) = U(I,J,K).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               X = XS AND X = XF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               X = XS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO X IS
C                               SPECIFIED AT X = XF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED AT
C                               X = XS AND X = XF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED AT
C                               X = XS AND THE SOLUTION IS SPECIFIED
C                               AT X=XF.
C
C                        BDXS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE DERIVATIVE OF THE SOLUTION
C                          WITH RESPECT TO X AT X = XS.
C
C                          WHEN LBDCND = 3 OR 4,
C
C                            BDXS(J,K) = (D/DX)U(XS,Y(J),Z(K)),
C                            J=1,2,...,M+1,      K=1,2,...,N+1.
C
C                          WHEN LBDCND HAS ANY OTHER VALUE, BDXS
C                          IS A DUMMY VARIABLE. BDXS MUST BE
C                          DIMENSIONED AT LEAST (M+1)*(N+1).
C
C                        BDXF
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE DERIVATIVE OF THE SOLUTION
C                          WITH RESPECT TO X AT X = XF.
C
C                          WHEN LBDCND = 2 OR 3,
C
C                            BDXF(J,K) = (D/DX)U(XF,Y(J),Z(K)),
C                            J=1,2,...,M+1,      K=1,2,...,N+1.
C
C                          WHEN LBDCND HAS ANY OTHER VALUE, BDXF IS
C                          A DUMMY VARIABLE.  BDXF MUST BE
C                          DIMENSIONED AT LEAST (M+1)*(N+1).
C
C                        YS,YF
C                          THE RANGE OF Y, I.E. YS .LE. Y .LE. YF.
C                          YS MUST BE LESS THAN YF.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (YS,YF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE M+1 GRID POINTS IN
C                          THE Y-DIRECTION GIVEN BY Y(J) = YS+(J-1)DY
C                          FOR J=1,2,...,M+1,
C                          WHERE DY = (YF-YS)/M IS THE PANEL WIDTH.
C                          M MUST BE AT LEAST 5.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Y = YS AND Y = YF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E.
C                               U(I,M+J,K) = U(I,J,K).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               Y = YS AND Y = YF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               Y = YS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO Y IS
C                               SPECIFIED AT Y = YF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = YS AND Y = YF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               AT Y = YS AND THE SOLUTION IS
C                               SPECIFIED AT Y=YF.
C
C                        BDYS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Y AT Y = YS.
C
C                          WHEN MBDCND = 3 OR 4,
C
C                            BDYS(I,K) = (D/DY)U(X(I),YS,Z(K)),
C                            I=1,2,...,L+1,      K=1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDYS
C                          IS A DUMMY VARIABLE. BDYS MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(N+1).
C
C                        BDYF
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Y AT Y = YF.
C
C                          WHEN MBDCND = 2 OR 3,
C
C                            BDYF(I,K) = (D/DY)U(X(I),YF,Z(K)),
C                            I=1,2,...,L+1,      K=1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDYF
C                          IS A DUMMY VARIABLE. BDYF MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(N+1).
C
C                        ZS,ZF
C                          THE RANGE OF Z, I.E. ZS .LE. Z .LE. ZF.
C                          ZS MUST BE LESS THAN ZF.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (ZS,ZF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS
C                          IN THE Z-DIRECTION GIVEN BY
C                          Z(K) = ZS+(K-1)DZ FOR K=1,2,...,N+1,
C                          WHERE DZ = (ZF-ZS)/N IS THE PANEL WIDTH.
C                          N MUST BE AT LEAST 5.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Z = ZS AND Z = ZF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.
C                               U(I,J,N+K) = U(I,J,K).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               Z = ZS AND Z = ZF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               Z = ZS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO Z IS
C                               SPECIFIED AT Z = ZF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Z IS SPECIFIED AT
C                               Z = ZS AND Z = ZF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Z IS SPECIFIED AT
C                               Z = ZS AND THE SOLUTION IS SPECIFIED
C                               AT Z=ZF.
C
C                        BDZS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Z AT Z = ZS.
C
C                          WHEN NBDCND = 3 OR 4,
C
C                            BDZS(I,J) = (D/DZ)U(X(I),Y(J),ZS),
C                            I=1,2,...,L+1,      J=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDZS
C                          IS A DUMMY VARIABLE. BDZS MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(M+1).
C
C                        BDZF
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Z AT Z = ZF.
C
C                          WHEN NBDCND = 2 OR 3,
C
C                            BDZF(I,J) = (D/DZ)U(X(I),Y(J),ZF),
C                            I=1,2,...,L+1,      J=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDZF
C                          IS A DUMMY VARIABLE. BDZF MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(M+1).
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION. IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HW3CRT WILL
C                          ATTEMPT TO FIND A SOLUTION.
C
C                        LDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE
C                          ARRAYS F,BDYS,BDYF,BDZS,AND BDZF AS IT
C                          APPEARS IN THE PROGRAM CALLING HW3CRT.
C                          THIS PARAMETER IS USED TO SPECIFY THE
C                          VARIABLE DIMENSION OF THESE ARRAYS.
C                          LDIMF MUST BE AT LEAST L+1.
C
C                        MDIMF
C                          THE COLUMN (OR SECOND) DIMENSION OF THE
C                          ARRAY F AND THE ROW (OR FIRST) DIMENSION
C                          OF THE ARRAYS BDXS AND BDXF AS IT APPEARS
C                          IN THE PROGRAM CALLING HW3CRT.  THIS
C                          PARAMETER IS USED TO SPECIFY THE VARIABLE
C                          DIMENSION OF THESE ARRAYS.
C                          MDIMF MUST BE AT LEAST M+1.
C
C                        F
C                          A THREE-DIMENSIONAL ARRAY OF DIMENSION AT
C                          AT LEAST (L+1)*(M+1)*(N+1), SPECIFYING THE
C                          VALUES OF THE RIGHT SIDE OF THE HELMHOLZ
C                          EQUATION AND BOUNDARY VALUES (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I=2,3,...,L,  J=2,3,...,M,
C                          AND K=2,3,...,N
C                          F(I,J,K) = F(X(I),Y(J),Z(K)).
C
C                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
C                          FOR J=1,2,...,M+1,  K=1,2,...,N+1,
C                          AND I=1,2,...,L+1
C
C                          LBDCND      F(1,J,K)         F(L+1,J,K)
C                          ------   ---------------   ---------------
C
C                            0      F(XS,Y(J),Z(K))   F(XS,Y(J),Z(K))
C                            1      U(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
C                            2      U(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))
C                            3      F(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))
C                            4      F(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
C
C                          MBDCND      F(I,1,K)         F(I,M+1,K)
C                          ------   ---------------   ---------------
C
C                            0      F(X(I),YS,Z(K))   F(X(I),YS,Z(K))
C                            1      U(X(I),YS,Z(K))   U(X(I),YF,Z(K))
C                            2      U(X(I),YS,Z(K))   F(X(I),YF,Z(K))
C                            3      F(X(I),YS,Z(K))   F(X(I),YF,Z(K))
C                            4      F(X(I),YS,Z(K))   U(X(I),YF,Z(K))
C
C                          NBDCND      F(I,J,1)         F(I,J,N+1)
C                          ------   ---------------   ---------------
C
C                            0      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZS)
C                            1      U(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
C                            2      U(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)
C                            3      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)
C                            4      F(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
C                          U AND THE RIGHT SIDE F ON A BOUNDARY,
C                          THEN THE SOLUTION MUST BE SPECIFIED.
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J,K) OF THE
C                          FINITE DIFFERENCE APPROXIMATION FOR THE
C                          GRID POINT (X(I),Y(J),Z(K)) FOR
C                          I=1,2,...,L+1, J=1,2,...,M+1,
C                          AND K=1,2,...,N+1.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
C                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
C                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
C                          CALCULATED AND SUBTRACTED FROM F, WHICH
C                          ENSURES THAT A SOLUTION EXISTS.  PWSCRT
C                          THEN COMPUTES THIS SOLUTION, WHICH IS A
C                          LEAST SQUARES SOLUTION TO THE ORIGINAL
C                          APPROXIMATION.  THIS SOLUTION IS NOT
C                          UNIQUE AND IS UNNORMALIZED.  THE VALUE OF
C                          PERTRB SHOULD BE SMALL COMPARED TO THE
C                          THE RIGHT SIDE F.  OTHERWISE, A SOLUTION
C                          IS OBTAINED TO AN ESSENTIALLY DIFFERENT
C                          PROBLEM.  THIS COMPARISON SHOULD ALWAYS
C                          BE MADE TO INSURE THAT A MEANINGFUL
C                          SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 12,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C                          =  1  XS .GE. XF
C                          =  2  L .LT. 5
C                          =  3  LBDCND .LT. 0 .OR. LBDCND .GT. 4
C                          =  4  YS .GE. YF
C                          =  5  M .LT. 5
C                          =  6  MBDCND .LT. 0 .OR. MBDCND .GT. 4
C                          =  7  ZS .GE. ZF
C                          =  8  N .LT. 5
C                          =  9  NBDCND .LT. 0 .OR. NBDCND .GT. 4
C                          = 10  LDIMF .LT. L+1
C                          = 11  MDIMF .LT. M+1
C                          = 12  LAMBDA .GT. 0
C                          = 20 If the dynamic allocation of real and
C                               complex work space required for solution
C                               fails (for example if N,M are too large
C                               for your computer)
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HW3CRT, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED Files         fish.f,pois3d.f,fftpack.f,comf.f
C
C LANGUAGE               FORTRAN 90
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980.
c                        Revised in June 2004 by John Adams using
c                        Fortran 90 dynamically allocated work space.
C
C PORTABILITY            FORTRAN 90
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS AND
C                        THEN CALLS POIS3D TO SOLVE THE SYSTEM.
C
C TIMING                 FOR LARGE L, M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO
C                          L*M*N*(LOG2(L)+LOG2(M)+5),
C                        BUT ALSO DEPENDS ON INPUT PARAMETERS LBDCND
C                        AND MBDCND.
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
C                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
C                        DIGITS FOR L, M AND N AS LARGE AS 32.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        ROUTINE POIS3D WHICH IS THE ROUTINE THAT
C                        ACTUALLY SOLVES THE FINITE DIFFERENCE
C                        EQUATIONS.
C
C REFERENCES             NONE
C***********************************************************************
      SUBROUTINE HW3CRT(XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M, MBDCND
     1   , BDYS, BDYF, ZS, ZF, N, NBDCND, BDZS, BDZF, ELMBDA, LDIMF, 
     2   MDIMF, F, PERTRB, IERROR)
      USE fish
      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: L
      INTEGER  :: LBDCND
      INTEGER  :: M
      INTEGER  :: MBDCND
      INTEGER  :: N
      INTEGER  :: NBDCND
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      REAL  :: XS
      REAL  :: XF
      REAL  :: YS
      REAL  :: YF
      REAL  :: ZS
      REAL  :: ZF
      REAL  :: ELMBDA
      REAL  :: PERTRB
      REAL  :: BDXS(MDIMF,*)
      REAL  :: BDXF(MDIMF,*)
      REAL  :: BDYS(LDIMF,*)
      REAL  :: BDYF(LDIMF,*)
      REAL  :: BDZS(LDIMF,*)
      REAL  :: BDZF(LDIMF,*)
      REAL  :: F(LDIMF,MDIMF,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IRWK, ICWK
C-----------------------------------------------
C
C     CHECK FOR INVALID INPUT.
C
      IERROR = 0
      IF (XF <= XS) IERROR = 1
      IF (L < 5) IERROR = 2
      IF (LBDCND<0 .OR. LBDCND>4) IERROR = 3
      IF (YF <= YS) IERROR = 4
      IF (M < 5) IERROR = 5
      IF (MBDCND<0 .OR. MBDCND>4) IERROR = 6
      IF (ZF <= ZS) IERROR = 7
      IF (N < 5) IERROR = 8
      IF (NBDCND<0 .OR. NBDCND>4) IERROR = 9
      IF (LDIMF < L + 1) IERROR = 10
      IF (MDIMF < M + 1) IERROR = 11
c     IF (IERROR .NE. 0) GO TO 188
      IF (IERROR /= 0) RETURN 
!     allocate required work space length (generous estimate)
      IRWK=30+L+M+5*N+MAX0(L,M,N)+7*(INT((L+1)/2)+INT((M+1)/2))
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call hw3crtt(XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,MBDCND,BDYS,
     +             BDYF,ZS,ZF,N,NBDCND,BDZS,BDZF,ELMBDA,LDIMF,
     +             MDIMF,F,PERTRB,IERROR,w%rew)
!     release allocated work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE HW3CRT


 
      SUBROUTINE HW3CRTT(XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M, 
     1   MBDCND, BDYS, BDYF, ZS, ZF, N, NBDCND, BDZS, BDZF, ELMBDA, 
     2   LDIMF, MDIMF, F, PERTRB, IERROR, W)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: L
      INTEGER  :: LBDCND
      INTEGER , INTENT(IN) :: M
      INTEGER  :: MBDCND
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: NBDCND
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER , INTENT(OUT) :: IERROR
      REAL , INTENT(IN) :: XS
      REAL , INTENT(IN) :: XF
      REAL , INTENT(IN) :: YS
      REAL , INTENT(IN) :: YF
      REAL , INTENT(IN) :: ZS
      REAL , INTENT(IN) :: ZF
      REAL , INTENT(IN) :: ELMBDA
      REAL , INTENT(OUT) :: PERTRB
      REAL , INTENT(IN) :: BDXS(MDIMF,*)
      REAL , INTENT(IN) :: BDXF(MDIMF,*)
      REAL , INTENT(IN) :: BDYS(LDIMF,*)
      REAL , INTENT(IN) :: BDYF(LDIMF,*)
      REAL , INTENT(IN) :: BDZS(LDIMF,*)
      REAL , INTENT(IN) :: BDZF(LDIMF,*)
      REAL  :: F(LDIMF,MDIMF,*)
      REAL  :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: MSTART, MSTOP, MP1, MP, MUNK, NP, NP1, NSTART, NSTOP, 
     1   NUNK, LP1, LP, LSTART, LSTOP, J, K, LUNK, I, IWB, IWC, IWW, 
     2   MSTPM1, LSTPM1, NSTPM1, NPEROD, IR
      REAL::DY,TWBYDY,C2,DZ,TWBYDZ,C3,DX,C1,TWBYDX,XLP,YLP,ZLP,S1,S2,S
C-----------------------------------------------
 
      DY = (YF - YS)/FLOAT(M)
      TWBYDY = 2./DY
      C2 = 1./DY**2
      MSTART = 1
      MSTOP = M
      MP1 = M + 1
      MP = MBDCND + 1
      GO TO (104,101,101,102,102) MP
  101 CONTINUE
      MSTART = 2
  102 CONTINUE
      GO TO (104,104,103,103,104) MP
  103 CONTINUE
      MSTOP = MP1
  104 CONTINUE
      MUNK = MSTOP - MSTART + 1
      DZ = (ZF - ZS)/FLOAT(N)
      TWBYDZ = 2./DZ
      NP = NBDCND + 1
      C3 = 1./DZ**2
      NP1 = N + 1
      NSTART = 1
      NSTOP = N
      GO TO (108,105,105,106,106) NP
  105 CONTINUE
      NSTART = 2
  106 CONTINUE
      GO TO (108,108,107,107,108) NP
  107 CONTINUE
      NSTOP = NP1
  108 CONTINUE
      NUNK = NSTOP - NSTART + 1
      LP1 = L + 1
      DX = (XF - XS)/FLOAT(L)
      C1 = 1./DX**2
      TWBYDX = 2./DX
      LP = LBDCND + 1
      LSTART = 1
      LSTOP = L
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      GO TO (122,109,109,112,112) LP
  109 CONTINUE
      LSTART = 2
      F(2,MSTART:MSTOP,NSTART:NSTOP) = F(2,MSTART:MSTOP,NSTART:NSTOP) - 
     1   C1*F(1,MSTART:MSTOP,NSTART:NSTOP)
      GO TO 115
  112 CONTINUE
      F(1,MSTART:MSTOP,NSTART:NSTOP) = F(1,MSTART:MSTOP,NSTART:NSTOP) + 
     1   TWBYDX*BDXS(MSTART:MSTOP,NSTART:NSTOP)
  115 CONTINUE
      GO TO (122,116,119,119,116) LP
  116 CONTINUE
      F(L,MSTART:MSTOP,NSTART:NSTOP) = F(L,MSTART:MSTOP,NSTART:NSTOP) - 
     1   C1*F(LP1,MSTART:MSTOP,NSTART:NSTOP)
      GO TO 122
  119 CONTINUE
      LSTOP = LP1
      F(LP1,MSTART:MSTOP,NSTART:NSTOP) = F(LP1,MSTART:MSTOP,NSTART:NSTOP
     1   ) - TWBYDX*BDXF(MSTART:MSTOP,NSTART:NSTOP)
  122 CONTINUE
      LUNK = LSTOP - LSTART + 1
C
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
C
      GO TO (136,123,123,126,126) MP
  123 CONTINUE
      F(LSTART:LSTOP,2,NSTART:NSTOP) = F(LSTART:LSTOP,2,NSTART:NSTOP) - 
     1   C2*F(LSTART:LSTOP,1,NSTART:NSTOP)
      GO TO 129
  126 CONTINUE
      F(LSTART:LSTOP,1,NSTART:NSTOP) = F(LSTART:LSTOP,1,NSTART:NSTOP) + 
     1   TWBYDY*BDYS(LSTART:LSTOP,NSTART:NSTOP)
  129 CONTINUE
      GO TO (136,130,133,133,130) MP
  130 CONTINUE
      F(LSTART:LSTOP,M,NSTART:NSTOP) = F(LSTART:LSTOP,M,NSTART:NSTOP) - 
     1   C2*F(LSTART:LSTOP,MP1,NSTART:NSTOP)
      GO TO 136
  133 CONTINUE
      F(LSTART:LSTOP,MP1,NSTART:NSTOP) = F(LSTART:LSTOP,MP1,NSTART:NSTOP
     1   ) - TWBYDY*BDYF(LSTART:LSTOP,NSTART:NSTOP)
  136 CONTINUE
      GO TO (150,137,137,140,140) NP
  137 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,2) = F(LSTART:LSTOP,MSTART:MSTOP,2) - 
     1   C3*F(LSTART:LSTOP,MSTART:MSTOP,1)
      GO TO 143
  140 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,1) = F(LSTART:LSTOP,MSTART:MSTOP,1) + 
     1   TWBYDZ*BDZS(LSTART:LSTOP,MSTART:MSTOP)
  143 CONTINUE
      GO TO (150,144,147,147,144) NP
  144 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,N) = F(LSTART:LSTOP,MSTART:MSTOP,N) - 
     1   C3*F(LSTART:LSTOP,MSTART:MSTOP,NP1)
      GO TO 150
  147 CONTINUE
      F(LSTART:LSTOP,MSTART:MSTOP,NP1) = F(LSTART:LSTOP,MSTART:MSTOP,NP1
     1   ) - TWBYDZ*BDZF(LSTART:LSTOP,MSTART:MSTOP)
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
  150 CONTINUE
      IWB = NUNK + 1
      IWC = IWB + NUNK
      IWW = IWC + NUNK
      W(:NUNK) = C3
      W(IWC:NUNK-1+IWC) = C3
      W(IWB:NUNK-1+IWB) = (-2.*C3) + ELMBDA
      GO TO (155,155,153,152,152) NP
  152 CONTINUE
      W(IWC) = 2.*C3
  153 CONTINUE
      GO TO (155,155,154,154,155) NP
  154 CONTINUE
      W(IWB-1) = 2.*C3
  155 CONTINUE
      PERTRB = 0.
C
C     FOR SINGULAR PROBLEMS ADJUST DATA TO INSURE A SOLUTION WILL EXIST.
C
      GO TO (156,172,172,156,172) LP
  156 CONTINUE
      GO TO (157,172,172,157,172) MP
  157 CONTINUE
      GO TO (158,172,172,158,172) NP
  158 CONTINUE
      IF (ELMBDA >= 0.) THEN
         IF (ELMBDA /= 0.) THEN
            IERROR = 12
         ELSE
            MSTPM1 = MSTOP - 1
            LSTPM1 = LSTOP - 1
            NSTPM1 = NSTOP - 1
            XLP = (2 + LP)/3
            YLP = (2 + MP)/3
            ZLP = (2 + NP)/3
            S1 = 0.
            DO K = 2, NSTPM1
               DO J = 2, MSTPM1
                  S1 = S1 + SUM(F(2:LSTPM1,J,K))
                  S1 = S1 + (F(1,J,K)+F(LSTOP,J,K))/XLP
               END DO
               S2 = SUM(F(2:LSTPM1,1,K)+F(2:LSTPM1,MSTOP,K))
               S2 = (S2 + (F(1,1,K)+F(1,MSTOP,K)+F(LSTOP,1,K)+F(LSTOP,
     1            MSTOP,K))/XLP)/YLP
               S1 = S1 + S2
            END DO
            S = (F(1,1,1)+F(LSTOP,1,1)+F(1,1,NSTOP)+F(LSTOP,1,NSTOP)+F(1
     1         ,MSTOP,1)+F(LSTOP,MSTOP,1)+F(1,MSTOP,NSTOP)+F(LSTOP,MSTOP
     2         ,NSTOP))/(XLP*YLP)
            DO J = 2, MSTPM1
               S = S + SUM(F(2:LSTPM1,J,1)+F(2:LSTPM1,J,NSTOP))
            END DO
            S2 = 0.
            S2 = SUM(F(2:LSTPM1,1,1)+F(2:LSTPM1,1,NSTOP)+F(2:LSTPM1,
     1         MSTOP,1)+F(2:LSTPM1,MSTOP,NSTOP))
            S = S2/YLP + S
            S2 = 0.
            S2 = SUM(F(1,2:MSTPM1,1)+F(1,2:MSTPM1,NSTOP)+F(LSTOP,2:
     1         MSTPM1,1)+F(LSTOP,2:MSTPM1,NSTOP))
            S = S2/XLP + S
            PERTRB = (S/ZLP + S1)/((FLOAT(LUNK + 1) - XLP)*(FLOAT(MUNK
     1          + 1) - YLP)*(FLOAT(NUNK + 1) - ZLP))
            F(:LUNK,:MUNK,:NUNK) = F(:LUNK,:MUNK,:NUNK) - PERTRB
         ENDIF
      ENDIF
  172 CONTINUE
      NPEROD = 0
      IF (NBDCND /= 0) THEN
         NPEROD = 1
         W(1) = 0.
         W(IWW-1) = 0.
      ENDIF
      CALL POIS3DD (LBDCND, LUNK, C1, MBDCND, MUNK, C2, NPEROD, NUNK, W
     1   , W(IWB), W(IWC), LDIMF, MDIMF, F(LSTART,MSTART,NSTART), IR, W(
     2   IWW))
C
C     FILL IN SIDES FOR PERIODIC BOUNDARY CONDITIONS.
C
      IF (LP == 1) THEN
         IF (MP == 1) THEN
            F(1,MP1,NSTART:NSTOP) = F(1,1,NSTART:NSTOP)
            MSTOP = MP1
         ENDIF
         IF (NP == 1) THEN
            F(1,MSTART:MSTOP,NP1) = F(1,MSTART:MSTOP,1)
            NSTOP = NP1
         ENDIF
         F(LP1,MSTART:MSTOP,NSTART:NSTOP) = F(1,MSTART:MSTOP,NSTART:
     1      NSTOP)
      ENDIF
      IF (MP == 1) THEN
         IF (NP == 1) THEN
            F(LSTART:LSTOP,1,NP1) = F(LSTART:LSTOP,1,1)
            NSTOP = NP1
         ENDIF
         F(LSTART:LSTOP,MP1,NSTART:NSTOP) = F(LSTART:LSTOP,1,NSTART:
     1      NSTOP)
      ENDIF
      IF (NP == 1) THEN
         F(LSTART:LSTOP,MSTART:MSTOP,NP1) = F(LSTART:LSTOP,MSTART:MSTOP,
     1      1)
      ENDIF
      RETURN 
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
c June      2004    Version 5.0, Fortran 90 changes
C-----------------------------------------------------------------------
      END SUBROUTINE HW3CRTT
