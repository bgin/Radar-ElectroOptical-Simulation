
#if 0
*    This file contains RANRTH and supporting functions and subroutines,
*      and a simple test program RANTST.
*
      PROGRAM RANTST
*
*   Simple Test Program for RANRTH
*
      EXTERNAL FUNSUB
      INTEGER I, NV, RS, IN, NI, NP, N, KEY, NF, LW
      PARAMETER ( N = 3, NI = 20000, NF = 4, LW = 20000 )
      DOUBLE PRECISION VALS(NF), ERRS(NF), WORK(LW)
      PRINT '('' Test for N = '', I2 )', N
      PRINT *, 'Key F Evals.\\ Values'
      DO KEY = 1, 4
         NP = NI
         RS = 0
         DO I = 1, 3
            CALL RANRTH( N, NF, NP, FUNSUB, 1D-5, 1D-5, RS, KEY, 
     &           VALS, ERRS, NV, IN, WORK )
            PRINT '(I3, I9, 1X, 5F8.5)', KEY, NV, VALS
            PRINT '(''     Errors  '', 5F8.5)', ERRS
            NP = NP*4
            RS = 1
          END DO
      END DO
      END
      SUBROUTINE FUNSUB( N, X, NF, F )
      INTEGER I, N, NF
      DOUBLE PRECISION X(N), F(NF), S
      F(1) = 1
      F(2) = 0
      F(3) = 0
      F(4) = 0
      S = SQRT( 2*X(1)**2 + 3*X(2)**2 + 4*X(3)**2 )
      IF ( S .GT. 1E-12 ) THEN
         F(2) = X(1)**2/S
         F(3) = X(2)**2/S
         F(4) = X(3)**2/S
      END IF
      END
*
#endif
      SUBROUTINE RANRTH( M, NF, MXVALS, F, EPSABS, EPSREL, RS, 
     &                   KEY, VALUE, ERROR, INTVLS, INFORM, WK )
****BEGIN PROLOGUE RANRTH
****AUTHOR
*            Alan Genz, Department of Mathematics, Washington State
*            University, Pullman, WA 99164-3113, USA
*              Email: AlanGenz@wsu.edu
*
*        Reference:
*            Genz, A., and Monahan, J. (1998), 
*             Stochastic Integration Rules for Infinite Regions,
*             {\it SIAM J. Sci. Com.}, {\bf 19}, pp. 426--439.
*
****KEYWORDS automatic multidimensional integrator,
*            n-dimensional region ( -infin, infin )^n, Gaussian weight
****PURPOSE  The routine calculates an approximation to a given
*            vector of definite integrals
*
*
*      infin     infin 
*     I    ...  I     w(X)( F ,F ,...,F   ) DX(M)...DX(2)DX(1),
*     -infin    -infin       1  2      NF
*
*       where F = F ( X ,X ,...,X  ), I = 1,2,...,NF,
*              I   I   1  2      M
*
*       w(X) = (2PI)^(-M/2)EXP(-( X(1)**2 + ... + X(M)**2 )/2),
*
*            hopefully satisfying for K = 1, 2, ... NF
*            ABS( I(K)-VALUE(K) ) .LE. MAX( EPSABS, EPSREL*ABS(I(K)) )
*
****DESCRIPTION Computation of integrals over infinite regions with
*               Gaussian weight function.
*
*   ON ENTRY
*
*     M  Integer number of variables, M > 1.
*     NF Integer number of components of the integral NF >= 1.
*     MXVALS Integer maximum number of F calls.
*            When RS > 0, this is the maximum number of new F calls.
*     F Externally declared subroutine for computing all components 
*            of the integrand at the given evaluation point.
*            It must have parameters ( M, X, NF, FUNS )
*            Input parameters:
*              M   Integer number of variables.
*              X   Real array of length M, the evaluation point.
*              NF Integer number of components for I.
*            Output parameter:
*              FUNS Real array of length NF, components of the integrand
*               evaluated at the point X.
*     EPSABS Real requested absolute accuracy.
*     EPSREL Real requested relative accuracy.
*     RS Integer.
*            If RS = 0, this is the first attempt to compute the integral(s).
*            If RS = 1, then a previous calculation is continued. In 
*              this case, the only parameters that may be changed (with 
*              respect to the previous call of the subroutine) are  
*              MXVALS, EPSABS, EPSREL and KEY.
*     KEY  Integer, determines degree of integration rules.
*            If KEY = 1, a degree 1 rule is used, and a minimum of
*                20 calls of F are used.
*            If KEY = 2, a degree 3 rule is used, and a minimum of
*                6( M + 1 ) + 1 calls of F are used.
*            If KEY = 3, a degree 5 rule is used, and a minimum of
*                4( M + 1 )( M + 2 ) + 1 calls of F are used.
*            If KEY = 4, a degree 7 spherical rule is used, 
*                a degree 5 radial rule is used , and a minimum of
*                4( M + 1 )( M^2 +8M + 6 )/3 + 1 calls of F are used.
*     WK   Real work array, with length at least 6NF + 2M for KEY = 1,
*          and with length at least 6NF + M( M + 2 ) for KEY > 1.
*
*   ON RETURN
*
*     VALUE Real array of length NF of approximations to the 
*            components of the integral.
*     ERROR Real array of length NF of estimates of absolute accuracies.
*     INTVLS Integer number of F calls used.
*            When RS > 0, this is the number of new F calls.
*     INFORM  Integer.
*            INFORM = 0 for normal exit, when 
*              ERROR(K) <= MAX( EPSABS, ABS( VALUE(K) )EPSREL ) for
*              0 < K <= NF, with <= MXVALS function values. 
*            INFORM = 1 if MXVALS was too small to obtain the required 
*              accuracy. In this case values of VALUE are returned
*              with estimated absolute accuracies ERROR.
*            INFORM = 2 if MXVALS was less than the minimum required for 
*              the rule type specified by KEY. In this case all 
*              components of VALUE are set to 0 and all components of
*              ERROR are set to -1.
*
****END PROLOGUE RANRTH
*
*   Global variables.
*
      EXTERNAL F
      INTEGER M, NF, MXVALS, RS, KEY
      INTEGER INTVLS, INFORM
      DOUBLE PRECISION EPSABS, EPSREL
      DOUBLE PRECISION VALUE(*), ERROR(*), WK(*)

*   Local variables.
*
      INTEGER I, I1, I2, I3, I4, I5, I6, I7
      DOUBLE PRECISION WEIGHT
*
*   Divide up work space and call RANBAS
*
      I1 =  1 + NF
      I2 = I1 + NF
      I3 = I2 + M
      I4 = I3 + NF
      I5 = I4 + NF
      I6 = I5 + NF
      I7 = I6 + NF
      IF ( RS .EQ. 0 ) THEN
         DO I = 1, NF
            VALUE(I) = 0
            ERROR(I) = -1
         END DO
      END IF
      CALL RANBAS( M, NF, F, MXVALS, KEY, WK, WK(I1), INTVLS,
     &             WK(I2), WK(I3), WK(I4), WK(I5), WK(I6), WK(I7) )
      IF ( INTVLS .LE. MXVALS ) THEN
         INFORM = 0
         DO I = 1, NF
            IF ( ERROR(I) .GT. 0 ) THEN 
               WEIGHT = 1/( 1 + WK(NF+I)/ERROR(I)**2 )
            ELSE 
               WEIGHT = 1
            END IF
            VALUE(I) = VALUE(I) + WEIGHT*( WK(I) - VALUE(I) )
            ERROR(I) = SQRT( WEIGHT*WK(NF+I) )
            IF ( ERROR(I) .GT. MAX( EPSABS, EPSREL*ABS(VALUE(I) ) ) )
     &           INFORM = 1
         END DO
      ELSE
         INTVLS = 0
         INFORM = 2
      END IF
*
****END RANRTH
*
      END
*
      SUBROUTINE RANBAS( M, NF, F, MXVALS, KEY, VALUE, ERROR,
     &                   INTVLS, X, FUNS, FUNC, FUNV, WK, V )
*
*     Stochastic Spherical Radial Rule, for
*
*                         INF                     
*        (2PI)^(-M/2) I  I   exp(-r^2/2) r^(M-1) G(rZ) dr dZ.
*                      U  0       
*
*     U is the surface of unit M-sphere, 
*     Z is an M-vector and G is an NF-vector valued function.
*     In this subroutine, F is a subroutine with calling sequence:
*               CALL F( M, X, NF, G ).
* 
*        Author:
*            Alan Genz, Department of Mathematics, Washington State
*            University, Pullman, WA 99164-3113, USA
*              Email: alangenz@wsu.edu
*
*        References:
*            Genz, A., and Monahan, J. (1997), 
*             A Stochastic Algorithm for High Dimensional Multiple 
*             Integrals over Unbounded Regions with Gaussian Weight,
*             to appear in {\it J. Comp. Appl. Math.}.
*            Genz, A., and Monahan, J. (1998), 
*             Stochastic Integration Rules for Infinite Regions,
*             {\it SIAM J. Sci. Com.}, {\bf 19}, pp. 426--439.
*
*
*   ON ENTRY
*
*     M  Integer number of variables, M > 1.
*     NF Integer number of components of the integral, NF >= 1.
*     MXVALS Integer maximum number of F calls.
*            When RS > 0, this is the maximum number of new F calls.
*     F Externally declared subroutine for computing all components 
*            of the integrand at the given evaluation point.
*            It must have parameters ( M, X, NF, FUNS )
*            Input parameters:
*              M   Integer number of variables.
*              X   Real array of length M, the evaluation point.
*              NF Integer number of components for I.
*            Output parameter:
*              FUNS Real array of length NF, the components of the 
*                      integrand evaluated at the point X.
*     KEY  Integer, determines degree of integration rules.
*            If KEY = 1, a degree 1 rule is used, and a minimum of
*                20 calls of F are used.
*            If KEY = 2, a degree 3 rule is used, and a minimum of
*                6( M + 1 ) + 1 calls of F are used.
*            If KEY = 3, a degree 5 rule is used, and a minimum of
*                4( M + 1 )( M + 2 ) + 1 calls of F are used.
*            If KEY = 4, a degree 7 spherical rule is used, 
*                a degree 5 radial rule is used , and a minimum of
*                4( M + 1 )( M^2 +8M + 6 )/3 + 1 calls of F are used.
*     X, FUNS, FUNC, FUNV, WK and V Real work arrays, with respective 
*                lengths M, NF, NF, NF, NF, and M( M + 1 ).
*
*   ON RETURN
*
*     VALUE Real array of length NF of approximations to the 
*            components of the integral.
*     ERROR Real array of length NF of squares of standard errors.
*     INTVLS Integer number of F calls used.
*
      EXTERNAL F
      INTEGER M, NF, MXVALS, KEY, INTVLS, I, L, NS, NV
      DOUBLE PRECISION VALUE(*), ERROR(*)
      DOUBLE PRECISION X(*), FUNS(*), FUNC(*), FUNV(*), V(*), WK(*) 
      DOUBLE PRECISION R, Q, RM, TH, WC, WV, DIFFER
      DOUBLE PRECISION RNRNOR, BETRAN, GAMRAN
      IF ( KEY .EQ. 2 ) THEN
         NS = MAX( 3, ( MXVALS - 1 )/( 2*( M + 1 ) ) )
         INTVLS = 1 + NS*2*( M + 1 )
      ELSE IF ( KEY .EQ. 3 ) THEN
         NV = 2*( M + 1 )*( M + 2 )
         IF ( M .EQ. 2 ) NV = 12
         IF ( M .EQ. 3 ) NV = 28
         NS = MAX( 2, ( MXVALS - 1 )/NV )
         INTVLS = 1 + NS*NV
      ELSE IF ( KEY .EQ. 4 ) THEN
         NV = 2*( M + 1 )*( M*M +8*M + 6 )/3
         IF ( M .EQ. 2 ) NV = 36
         IF ( M .EQ. 3 ) NV = 76
         IF ( M .EQ. 4 ) NV = 140
         IF ( M .EQ. 5 ) NV = 244
         NS = MAX( 2, ( MXVALS - 1 )/NV )
         INTVLS = 1 + NS*NV
      ELSE
         NS = MAX( 10, MXVALS/2 )
         INTVLS = 2*NS
         WV = 1
      END IF
      IF ( INTVLS .LE. MXVALS ) THEN
         DO I = 1,NF
            VALUE(I) = 0
            ERROR(I) = 0
         END DO
         IF ( 2 .LE. KEY .AND. KEY .LE. 4 ) THEN
            DO I = 1, M
               X(I) = 0
            END DO
            CALL F( M, X, NF, FUNC )
         END IF
         DO L = 1, NS
            IF ( KEY .LE. 1 .OR. 5 .LE. KEY ) THEN
               DO I = 1, M
                  V(I) = RNRNOR()
               END DO
               CALL RNSRUL( KEY, M, NF, F, R, V, X, WK, FUNV )
            ELSE
               RM = M + 2
               IF( KEY .EQ. 2 ) THEN
                  R = SQRT( 2*GAMRAN( RM/2 ) )
                  WV = M/R**2
                  WC = 1 - WV
               ELSE 
                  TH = 3
                  TH = ASIN( BETRAN( RM, TH/2 ) )/2
                  RM = 2*M + 7
                  R = SQRT( 2*GAMRAN( RM/2 ) )
                  Q = R*COS(TH)
                  R = R*SIN(TH)
                  WC = 1 + M*( M + 2 - R**2 - Q**2 )/( R*Q )**2
                  WV = M*( M + 2 - Q**2 )/( R**2*( R**2 - Q**2 ) )
               END IF
               DO I = 1,NF
                  FUNV(I) = WC*FUNC(I)
               END DO
               CALL RNSIMP( M, V, X )
*     
*              Compute integrand average
*     
               CALL RNSRUL( KEY, M, NF, F, R, V, X, WK, FUNS )
               DO I = 1, NF
                  FUNV(I) = FUNV(I) + WV*FUNS(I)
               END DO
               IF ( KEY .GT. 2 ) THEN
                  WV = M*( M + 2 - R**2 )/( Q**2*( Q**2 - R**2 ) )
                  CALL RNSRUL( KEY, M, NF, F, Q, V, X, WK, FUNS )
                  DO I = 1, NF
                     FUNV(I) = FUNV(I) + WV*FUNS(I)
                  END DO
               END IF 
            END IF
            DO I = 1, NF
               DIFFER = ( FUNV(I) - VALUE(I) )/L
               VALUE(I) = VALUE(I) + DIFFER
               ERROR(I) = ( L - 2 )*ERROR(I)/L + DIFFER**2 
            END DO 
         END DO
      END IF
      END
*
      SUBROUTINE RNSRUL( KEY, M, NF, F, R, V, X, RESR, INTV )
*
*     Spherical Radial Rule, for
*            
*         I  G(Z) dZ .
*          U 
*
*     U is the surface of a radius R M-sphere, Z is an M-vector, 
*     G is an NF-vector valued function.
*       In this subroutine, F is a subroutine with calling sequence:
*               CALL F( M, X, NF, G ).
*     KEY is an integer rule degree parameter with 1 <= KEY <= 4. 
*       A degree 2*KEY-1 stochastic spherical rule is used.
*     Output INTV is an NF-vector of integral estimates. 
*     Work vectors V, X and RESR must have respective lengths at 
*      least M*(M+1), M and NF.
* 
      EXTERNAL F
      INTEGER KEY, M, NF, I, IS, J, K, N
      DOUBLE PRECISION R, V( M, * ), X(*), RESR(*), INTV(*)
      DOUBLE PRECISION MP, WV, WM, WC, WT, RM, RC, RT
*
*     Determine Weights
*
      MP = M + 1
      IF ( KEY .EQ. 2 ) THEN
         WV = 1/( 2*MP )
      ELSE IF ( KEY .EQ. 3 ) THEN
         WV =   ( 7 - M )*M /( 2*( M + 2 )*MP**2 )
         WM = 2*( M - 1 )**2/( M*( M + 2 )*MP**2 )
         RM = R/SQRT( 2*( MP - 2 )/M )
         IF ( M .EQ. 2 ) WV = WV + WM
         IF ( M .EQ. 3 ) WM = 2*WM
      ELSE IF ( KEY .EQ. 4 ) THEN
         WV =   1/( 36*M*MP**3*( M + 2 )*( M + 4 ) )
         WM = 144*( M - 1)**3*( 4 - M )*WV
         WC = 486*( M - 2 )**3*WV
         WT =  ( 10*M - 6 )**3*WV
         WV = M**3*( 1800 - 793*M + 9*M*M )*WV
         RM = R/SQRT( 2*( MP - 2 )/M )
         IF ( M .GT. 2 ) RC = R/SQRT( 3*( MP - 3 )/M )
         RT = R/SQRT( ( 10*MP - 16 )/M )
         IF ( M .EQ. 2 ) WV = WV + WM
         IF ( M .EQ. 3 ) WV = WV + WC
         IF ( M .EQ. 3 ) WM = 2*WM
         IF ( M .EQ. 4 ) WM = WC
         IF ( M .EQ. 5 ) WC = 2*WC
      ELSE 
         WV = 1
         WV = WV/2
      ENDIF
      DO I = 1,NF
         INTV(I) = 0
      END DO
*     
*     Compute integrand average
*     
      DO IS = -1, 1, 2
         IF ( KEY .LE. 1 .OR. 5 .LE. KEY ) THEN
            DO I = 1, M
               X(I) = IS*V(I,1)
            END DO
            CALL F( M, X, NF, RESR )
            DO I = 1, NF
               INTV(I) = INTV(I)+ WV*RESR(I)
            END DO
         ELSE 
            DO K = 1, M+1
               DO I = 1, M
                  X(I) = IS*R*V(I,K)
               END DO
               CALL F( M, X, NF, RESR )
               DO I = 1, NF
                  INTV(I) = INTV(I) + WV*RESR(I) 
               END DO
            END DO
         END IF
         IF ( ( KEY .EQ. 3 .OR. KEY .EQ. 4 ) .AND. 
     &        (   M .GT. 3 .OR.   M .EQ. 3 .AND. IS .EQ. 1 ) ) THEN     
            DO K = 1, M
               DO J = K+1, M+1
                  DO I = 1, M
                     X(I) = IS*RM*( V(I,K) + V(I,J) )
                  END DO
                  CALL F( M, X, NF, RESR )
                  DO I = 1, NF
                     INTV(I) = INTV(I) + WM*RESR(I)
                  END DO
               END DO
            END DO
         END IF
         IF ( KEY .EQ. 4 ) THEN
            DO K = 1, M+1
               DO J = 1, M+1
                  IF ( J .NE. K ) THEN
                     DO I = 1, M
                        X(I) = IS*RT*( V(I,K) + 3*V(I,J) )
                     END DO
                     CALL F( M, X, NF, RESR )
                     DO I = 1, NF
                        INTV(I) = INTV(I) + WT*RESR(I)
                     END DO
                  END IF
               END DO
            END DO
            IF ( M .GT. 5 .OR. M .EQ. 5 .AND. IS .EQ. 1 ) THEN
               DO K = 1, M-1
                  DO J = K+1, M
                     DO N = J+1, M+1
                        DO I = 1, M
                           X(I) = IS*RC*( V(I,K) + V(I,J) + V(I,N) )
                        END DO
                        CALL F( M, X, NF, RESR )
                        DO I = 1, NF
                           INTV(I) = INTV(I) + WC*RESR(I)
                        END DO
                     END DO
                  END DO
               END DO
            END IF
         END IF
      END DO
      END
*
      SUBROUTINE RNSIMP( M, V, X )
*
*     Determine random M-simplex with vertices V and work vector X.
*
      INTEGER I, J, K, M
      DOUBLE PRECISION V( M, * ), X(*), AL, BT, RV, MP, RNRNOR
      MP = M + 1
*
*     Determine standard unit simplex centered at origin
*
      DO I = 1, M
         DO J = 1, I-1
            V(I,J) = 0
         END DO
         RV = SQRT( MP/( ( M - I + 1 )*M*( M - I + 2 ) ) )
         V(I,I) = ( M - I + 1 )*RV
         DO J = I+1, M+1
            V(I,J) = -RV
         END DO
      END DO
*     
*     Replace V with (random orthogonal matrix)*V
*     This generates the random orthogonal matrix using a 
*     sequence of random Householder transformations.
*     Reference: Stewart, G.W. (1980),
*     "The Efficient Generation of Random Orthogonal Matrices 
*      with an Application to Condition Estimaors",
*     SIAM J Numer. Anal. 17, pp. 403-408.  
*     
      DO K = M - 1, 1, -1
         AL = 0
         DO I = K, M
            X(I) = RNRNOR()
            AL = AL + X(I)**2
         END DO
         AL = -SQRT(AL)
         BT = 1/( AL*( AL + X(K) ) )
         X(K) = X(K) + AL
         DO J = K, M+1
            AL = 0
            DO I = K, M
               AL = AL + X(I)*V(I,J)
            END DO
            AL = BT*AL
            DO I = K, M
               V(I,J) = V(I,J) - X(I)*AL
            END DO
         END DO
      END DO
      END
*
      DOUBLE PRECISION FUNCTION GAMRAN( ALPHA )
*
*     For Gamma random numbers
*
*     Reference: P Tadikamalla, "Computer Generation of Gamma Random
*                Variables-II", CACM 21(1978), pp. 925-928.
*
      DOUBLE PRECISION U, RNRUNI, X, Y, E
      PARAMETER ( E = 2.71828182845904523536D0 )
      DOUBLE PRECISION ALPHA, A, B, C, D, OLDALP, P, S
      SAVE                    A, B, C, D, OLDALP, P, S
      DATA OLDALP/ 0D0 / 
      IF ( ALPHA .NE. OLDALP ) THEN
         OLDALP = ALPHA
         IF ( ALPHA .LE. 1 ) THEN
            B = ( E + ALPHA )/E
         ELSE
            A = ALPHA - 1
            B = ( 1 + SQRT( 4*ALPHA - 3 ) )/2
            C = A*( 1 + B )/B
            D = ( B - 1 )/( A*B )
            S = A/B
            P = 1/( 2 - EXP(-S) )
         END IF
      END IF
 10   U = RNRUNI( )
      IF ( ALPHA .LE. 1 ) THEN
         Y = B*RNRUNI( )
         IF ( Y .LT. 1 ) THEN
            X = EXP( LOG( Y )/ALPHA )
            IF ( U .GT. EXP( - X ) ) GO TO 10
         ELSE
            X = -LOG( ( B - Y )/ALPHA )
            IF ( LOG( U ) .GT. ( ALPHA - 1 )*LOG( X ) ) GO TO 10
         END IF
      ELSE
         IF ( U .GT. P ) THEN
            Y = -LOG( ( 1 - U )/( 1 - P ) )
            IF ( Y .GT. S ) Y = MOD( Y, S )
            X = A - B*Y
         ELSE
            Y = -LOG( U/P )
            X = A + B*Y
         END IF
         IF ( LOG( RNRUNI( ) ) .GT. A*LOG(D*X) - X + Y + C ) GO TO 10
      END IF
      GAMRAN = X
      END
      DOUBLE PRECISION FUNCTION BETRAN( ALPHA, BETA )
*
*     For random numbers with Beta distribution, uses two Gamma distributed
*     random numbers GAMRAN( ALPHA) and GAMRAN( BETA ) to obtain
*       BETRAN = GAMRAN( ALPHA )/( GAMRAN( ALPHA ) + GAMRAN( BETA ) ).
*
      DOUBLE PRECISION U, RNRUNI, Y, E
      PARAMETER ( E = 2.71828182845904523536D0 )
      DOUBLE PRECISION ALPHA, AA, BA, CA, DA, OLDALP, PA, SA, XA
      DOUBLE PRECISION  BETA, AB, BB, CB, DB, OLDBET, PB, SB, XB
      SAVE                    AA, BA, CA, DA, OLDALP, PA, SA 
      SAVE                    AB, BB, CB, DB, OLDBET, PB, SB 
      DATA OLDALP, OLDBET/ 2*0D0 / 
      IF ( ALPHA .NE. OLDALP ) THEN
         OLDALP = ALPHA
         IF ( ALPHA .LE. 1 ) THEN
            BA = ( E + ALPHA )/E
         ELSE
            AA = ALPHA - 1
            BA = ( 1 + SQRT( 4*ALPHA - 3 ) )/2
            CA = AA*( 1 + BA )/BA
            DA = ( BA - 1 )/( AA*BA )
            SA = AA/BA
            PA = 1/( 2 - EXP(-SA) )
         END IF
      END IF
      IF ( BETA .NE. OLDBET ) THEN
         OLDBET = BETA
         IF ( BETA .LE. 1 ) THEN
            BB = ( E + BETA )/E
         ELSE
            AB = BETA - 1
            BB = ( 1 + SQRT( 4*BETA - 3 ) )/2
            CB = AB*( 1 + BB )/BB
            DB = ( BB - 1 )/( AB*BB )
            SB = AB/BB
            PB = 1/( 2 - EXP(-SB) )
         END IF
      END IF
 10   U = RNRUNI( )
      IF ( ALPHA .LE. 1 ) THEN
         Y = BA*RNRUNI( )
         IF ( Y .LT. 1 ) THEN
            XA = EXP( LOG( Y )/ALPHA )
            IF ( U .GT. EXP( - XA ) ) GO TO 10
         ELSE
            XA = -LOG( ( BA - Y )/ALPHA )
            IF ( LOG( U ) .GT. ( ALPHA - 1 )*LOG( XA ) ) GO TO 10
         END IF
      ELSE
         IF ( U .GT. PA ) THEN
            Y = -LOG( ( 1 - U )/( 1 - PA ) )
            IF ( Y .GT. SA ) Y = MOD( Y, SA )
            XA = AA - BA*Y
         ELSE
            Y = -LOG( U/PA )
            XA = AA + BA*Y
         END IF
         IF ( LOG( RNRUNI( ) ) .GT.  AA*LOG(DA*XA)-XA+Y+CA ) GO TO 10
      END IF
 20   U = RNRUNI( )
      IF ( BETA .LE. 1 ) THEN
         Y = BB*RNRUNI( )
         IF ( Y .LT. 1 ) THEN
            XB = EXP( LOG( Y )/BETA )
            IF ( U .GT. EXP( -XB ) ) GO TO 20
         ELSE
            XB = -LOG( ( BB - Y )/BETA )
            IF ( LOG( U ) .GT. ( BETA - 1 )*LOG( XB ) ) GO TO 20
         END IF
      ELSE
         IF ( U .GT. PB ) THEN
            Y = -LOG( ( 1 - U )/( 1 - PB ) )
            IF ( Y .GT. SB ) Y = MOD( Y, SB )
            XB = AB - BB*Y
         ELSE
            Y = -LOG( U/PB )
            XB = AB + BB*Y
         END IF
         IF ( LOG( RNRUNI( ) ) .GT.  AB*LOG(DB*XB)-XB+Y+CB ) GO TO 20
      END IF
      BETRAN = XA/( XA + XB )
      END
*
      DOUBLE PRECISION FUNCTION RNRNOR()
*
*     RNRNOR generates normal random numbers with zero mean and unit
*     standard deviation, often denoted N(0,1),adapted from G. Marsaglia 
*     and W. W. Tsang: "A Fast, Easily Implemented Method for Sampling
*     from Decreasing or Symmetric Unimodal Density Functions"
*      SIAM J. Sci. Stat. Comput. 5(1984), pp. 349-359.
*
      INTEGER J, N, TN
      DOUBLE PRECISION TWOPIS, AA, B, C, XDN
      PARAMETER ( N = 64, TN = 2*N, TWOPIS = TN/2.506628274631000D0 )
      PARAMETER ( XDN = 0.3601015713011893D0, B = 0.4878991777603940D0 ) 
      PARAMETER (  AA =  12.37586029917064D0, C =  12.67705807886560D0 )
      DOUBLE PRECISION XT, XX, Y, RNRUNI
      DOUBLE PRECISION X(0:N)
      SAVE X
      DATA ( X(J), J = 0, 31 ) /
     &  0.3409450287039653D+00,  0.4573145918669259D+00,
     &  0.5397792816116612D+00,  0.6062426796530441D+00,
     &  0.6631690627645207D+00,  0.7136974590560222D+00,
     &  0.7596124749339174D+00,  0.8020356003555283D+00,
     &  0.8417226679789527D+00,  0.8792102232083114D+00,
     &  0.9148948043867484D+00,  0.9490791137530882D+00,
     &  0.9820004812398864D+00,  0.1013849238029940D+01,
     &  0.1044781036740172D+01,  0.1074925382028552D+01,
     &  0.1104391702268125D+01,  0.1133273776243940D+01,
     &  0.1161653030133931D+01,  0.1189601040838737D+01,
     &  0.1217181470700870D+01,  0.1244451587898246D+01,
     &  0.1271463480572119D+01,  0.1298265041883197D+01,
     &  0.1324900782180860D+01,  0.1351412509933371D+01,
     &  0.1377839912870011D+01,  0.1404221063559975D+01,
     &  0.1430592868502691D+01,  0.1456991476137671D+01,
     &  0.1483452656603219D+01,  0.1510012164318519D+01 /
      DATA ( X(J), J = 32, 64 ) /
     &  0.1536706093359520D+01,  0.1563571235037691D+01,
     &  0.1590645447014253D+01,  0.1617968043674446D+01,
     &  0.1645580218369081D+01,  0.1673525509567038D+01,
     &  0.1701850325062740D+01,  0.1730604541317782D+01,
     &  0.1759842199038300D+01,  0.1789622321566574D+01,
     &  0.1820009890130691D+01,  0.1851077020230275D+01,
     &  0.1882904397592872D+01,  0.1915583051943031D+01,
     &  0.1949216574916360D+01,  0.1983923928905685D+01,
     &  0.2019843052906235D+01,  0.2057135559990095D+01,
     &  0.2095992956249391D+01,  0.2136645022544389D+01,
     &  0.2179371340398135D+01,  0.2224517507216017D+01,
     &  0.2272518554850147D+01,  0.2323933820094302D+01,
     &  0.2379500774082828D+01,  0.2440221797979943D+01,
     &  0.2507511701865317D+01,  0.2583465835225429D+01,
     &  0.2671391590320836D+01,4*0.2776994269662875D+01 /
      Y = RNRUNI()
      J = MOD( INT( TN*RNRUNI() ), N )
      XT = X(J+1)
      RNRNOR = ( Y + Y - 1 )*XT
      IF ( ABS(RNRNOR) .GT. X(J) ) THEN
         XX = B*( XT - ABS(RNRNOR) )/( XT - X(J) )
         Y = RNRUNI()
         IF ( Y .GT. C - AA*EXP( -XX**2/2 ) ) THEN
            RNRNOR = SIGN( XX, RNRNOR )
         ELSE
            IF ( EXP(-XT**2/2)+Y/(TWOPIS*XT).GT.EXP(-RNRNOR**2/2) ) THEN
 10            XX = XDN*LOG( RNRUNI() )
               IF ( -2*LOG( RNRUNI() ) .LE. XX**2 ) GO TO 10
               RNRNOR = SIGN( X(N) - XX, RNRNOR )
            END IF
         END IF
      END IF
      END
*
      DOUBLE PRECISION FUNCTION RNRUNI()
*
*     Uniform (0, 1) random number generator
*
*     Reference:
*     L'Ecuyer, Pierre (1996), 
*     "Combined Multiple Recursive Random Number Generators"
*     Operations Research 44, pp. 816-822.
*
*
      INTEGER A12, A13, A21, A23, P12, P13, P21, P23
      INTEGER Q12, Q13, Q21, Q23, R12, R13, R21, R23
      INTEGER X10, X11, X12, X20, X21, X22, Z, M1, M2, H 
      DOUBLE PRECISION INVMP1
      PARAMETER (  M1 = 2147483647,  M2 = 2145483479 )
      PARAMETER ( A12 =   63308,    Q12 = 33921, R12 = 12979 )
      PARAMETER ( A13 = -183326,    Q13 = 11714, R13 =  2883 )
      PARAMETER ( A21 =   86098,    Q21 = 24919, R21 =  7417 )
      PARAMETER ( A23 = -539608,    Q23 =  3976, R23 =  2071 )
      PARAMETER ( INVMP1 = 4.656612873077392578125D-10 ) 
*                 INVMP1 = 1/(M1+1)
      SAVE X10, X11, X12, X20, X21, X22
      DATA       X10,      X11,      X12,      X20,      X21,      X22  
     &    / 15485857, 17329489, 36312197, 55911127, 75906931, 96210113 /      
*
*     Component 1
*
      H = X10/Q13
      P13 = -A13*( X10 - H*Q13 ) - H*R13
      H = X11/Q12
      P12 =  A12*( X11 - H*Q12 ) - H*R12
      IF ( P13 .LT. 0 ) P13 = P13 + M1
      IF ( P12 .LT. 0 ) P12 = P12 + M1
      X10 = X11 
      X11 = X12
      X12 = P12 - P13
      IF ( X12 .LT. 0 ) X12 = X12 + M1
*
*     Component 2
*
      H = X20/Q23
      P23 = -A23*( X20 - H*Q23 ) - H*R23
      H = X22/Q21
      P21 =  A21*( X22 - H*Q21 ) - H*R21
      IF ( P23 .LT. 0 ) P23 = P23 + M2
      IF ( P21 .LT. 0 ) P21 = P21 + M2
      X20 = X21 
      X21 = X22
      X22 = P21 - P23
      IF ( X22 .LT. 0 ) X22 = X22 + M2
*
*     Combination
*
      Z = X12 - X22
      IF ( Z .LE. 0 ) Z = Z + M1
      RNRUNI = Z*INVMP1
      END
