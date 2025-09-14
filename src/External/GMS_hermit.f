      SUBROUTINE HRMSYM( NDIM, NF, MINPTS, MAXPTS, FUNSUB, EPSABS, 
     &     EPSREL, RESTAR, RESULT, ABSERR, NEVAL, IFAIL, WORK)
****BEGIN PROLOGUE HRMSYM
*
****AUTHOR
*            Alan Genz, 
*            Department of Mathematics 
*            Washington State University
*            Pullman, WA 99164-3113, USA
*            Email: alangenz@wsu.edu
*
*      Reference: Genz, A., and Keister, B. (1996), Fully Symmetric 
*      Interpolatory Rules for Multiple Integrals over Infinite 
*      Regions with Gaussian Weight, J. Comp. Appl. Math., 71, 299-309.
*
****KEYWORDS automatic multidimensional integrator,
*            n-dimensional region (-infin, infin)^n 
*            Gaussian weight
****PURPOSE  The routine calculates an approximation to a given
*            vector of definite integrals
*
*      infin     infin 
*     I    ...  I     w(X)(F ,F ,...,F  ) DX(NDIM)...DX(2)DX(1),
*     -infin    -infin      1  2      NF
*
*       where F = F (X ,X ,...,X    ), I = 1,2,...,NF,
*              I   I  1  2      NDIM
*
*       w(X) = EXP(-( X(1)**2 + ... + X(NDIM)**2 )/2)/SQRT(2*PI)**NDIM.
*
****DESCRIPTION Computation of integrals over infinite regions with
*               Gaussian weight.
*
*   ON ENTRY
*
*     NDIM   Integer, number of variables, with 1 < NDIM <= 1000.
*     NF     Integer, number of components of the integral.
*     MINPTS Integer, minimum number of function evaluations.
*     MAXPTS Integer, maximum number of function evaluations.
*     FUNSUB Externally declared subroutine for computing
*            all components of the integrand at the given
*            evaluation point.
*            It must have parameters (NDIM,X,NF,FUNVLS)
*            Input parameters:
*              NDIM   Integer that defines the dimension of the
*                     integral.
*              X      Real array of dimension NDIM
*                     that defines the evaluation point.
*              NF Integer that defines the number of components of I.
*            Output parameter:
*              FUNVLS Real array of dimension NF
*                     that defines NF components of the integrand.
*
*     EPSABS Real, requested absolute accuracy.
*     EPSREL Real, requested relative accuracy.
*     RESTAR Integer.
*            If RESTAR = 0, this is the first attempt to compute
*            the integral.
*            If RESTAR = 1, then we restart a previous attempt.
*     WORK   Real array of working storage, must have dimension at
*            least 2*NF + (NF+1)*NUMSMS. NUMSMS is the number 
*            of NDIM partitions of the integers 0,...,RULE, after
*            RULE+1 approximations to the integrals have been computed.
*            The required value for NUMSMS will never exceed 10000.
*
*   ON RETURN
*
*     RESULT Real array of dimension NF.
*            Approximations to all components of the integral.
*     ABSERR Real array of dimension NF.
*            Estimates of absolute accuracies.
*     NEVAL  Integer, number of function evaluations used.
*     IFAIL  Integer.
*            IFAIL = 0 for normal exit, when 
*              ABSERR(K) <= MAX( EPSABS, ABS(RESULT(K))*EPSREL ) for
*              all K, 0 < K <= NF, with <= MAXPTS function values. 
*            IFAIL = 1 if MAXPTS was too small to obtain the required 
*              accuracy. In this case values of RESULT are returned
*              with estimated absolute accuracies ABSERR.
*     WORK   Real array of working storage, that contains information
*            that might be needed for subsequent calls of HRMSYM. This
*            array should not be modified. 
*
****ROUTINES CALLED 
****END PROLOGUE HRMSYM
*
*   Global variables.
*
      EXTERNAL FUNSUB
      INTEGER NDIM, NF, MINPTS, MAXPTS, RESTAR
      INTEGER NEVAL, IFAIL
      DOUBLE PRECISION EPSABS, EPSREL
      DOUBLE PRECISION RESULT(*), ABSERR(*), WORK(*)
*
*   Local variables.
*
      INTEGER I, MAXRUL, RULE, MNRULE, INTCLS, NUMSMS
      PARAMETER ( MAXRUL = 25 )
      SAVE RULE, MNRULE
      IFAIL = 1
      IF ( NDIM .LE. 4 ) THEN 
*
*     Call product Gauss-Hermite Rule
*
         CALL HERMIT( NDIM, NF, MINPTS, MAXPTS, FUNSUB, EPSABS, EPSREL, 
     &                RESTAR, RESULT, ABSERR, NEVAL, IFAIL, WORK )
      ELSE
         IF ( RESTAR .EQ. 0 ) THEN
            MNRULE = -1
            RULE = 0
            DO I = 1, NF
               WORK(I) = 0
            END DO
         END IF
         NEVAL = 0
         DO WHILE ( NEVAL .LE. MAXPTS .AND. RULE .LE. MAXRUL 
     &        .AND. ( IFAIL .GT. 0 .OR. NEVAL .LT. MINPTS ) ) 
            CALL HRMTRL( NDIM, NF, FUNSUB, MNRULE, RULE, RESULT, 
     &                   INTCLS, WORK(NF+1), WORK(2*NF+1), NUMSMS ) 
            IFAIL = 0
            DO I = 1, NF
               ABSERR(I) = ABS( RESULT(I) - WORK(I) )
               WORK(I) = RESULT(I)
               IF ( ABSERR(I) .GT. MAX(EPSABS,EPSREL*ABS(RESULT(I))) 
     &           .OR. RULE .EQ. 0 ) IFAIL = 1
            END DO
            NEVAL = NEVAL + INTCLS
            RULE = RULE + 1
         END DO
      END IF
*
****END HRMSYM
*
      END
      SUBROUTINE HRMTRL( S, N, F, MINORD, MAXORD, INTVAL, INTCLS, 
     &                   WORK, FULSMS, NUMSMS )
*  MULTIDIMENSIONAL FULLY SYMMETRIC RULE INTEGRATION SUBROUTINE
*
*   THIS SUBROUTINE COMPUTES A SEQUENCE OF FULLY SYMMETRIC RULE
*   APPROXIMATIONS TO A FULLY SYMMETRIC MULTIPLE INTEGRAL. WRITTEN BY 
*     Alan Genz
*     Department of Mathematics
*     Washington State University
*     Pullman, Washington 99164-3113  USA
*     Telephone: 509-335-2131 
*     Electronic Mail: alangenz@wsu.edu
*     Fax: 509-335-1188 
*
***************  PARAMETERS FOR HRMTRL  ********************************
******INPUT PARAMETERS
*  S       INTEGER number of variables, with 0 < S <= 1000.
*  N       Integer, number of components of the integral.
*  F       EXTERNALly declared user defined integrand subroutine.
*          it must have parameters (S,X,N,F), where X is a REAL S-array
*          and F is a real N-array.
*  MINORD  INTEGER minimum order parameter.  On entry MINORD specifies
*          the current highest order approximation to the integral,
*          available in the array INTVAL.  For the first call
*          MINORD should be set to -1. Otherwise a previous call is
*          assumed that computed INTVAL. On exit MINORD is set to MAXORD.
*  MAXORD  INTEGER maximum order parameter, must be greater than MINORD
*          and not exceed 25. The subroutine computes approximations of 
*          polynomial degree 2*MAXORD+1.
*******OUTPUT PARAMETERS
*  INTVAL  REAL array of length N.  Upon successful exit
*          INTVAL(1),..., INTVAL(N) are approximations to the components
*          of the integral.  These are all approximations of polynomial 
*          degree 2*MAXORD+1.
*  INTCLS  INTEGER number of F values needed for INTVAL
*  WORK    REAL working storage array. 
*  FULSMS  REAL working storage array with dimension (N+1,*). On exit
*          FULSMS(I,J) contains the fully symmetric basic rule sum
*          indexed by the jth s-partition of the integers 0,...,MAXORD, 
*          for the Ith component of the integrand.
*          FULSMS(N+1,J) contains number of points for the fully 
*          symmetric basic rule sum indexed by the Jth S-partition of 
*          the integers 0,...,MAXORD.
*  NUMSMS  INTEGER number of S-partitions of the integers 0,...,MAXORD.
************************************************************************
      EXTERNAL F
      INTEGER S, N, MAXDIM, MINORD, MAXORD, MAXRDM, NUMSMS
      PARAMETER ( MAXDIM = 1000, MAXRDM = 25 ) 
      DOUBLE PRECISION X(MAXDIM), FULWGT, WEIGHT, 
     &                 INTVAL(N), FULSMS(N+1,*), WORK(*) 
      INTEGER D, I, L, MODOFM, M(MAXDIM), K(MAXDIM), INTCLS, PRTCNT
      D = MINORD + 1
      INTCLS = 0
      IF ( D .EQ. 0 ) THEN
         DO I = 1,N
            INTVAL(I) = 0
         END DO
      END IF
*
****  Begin loop for each D
*      for each D find all distinct partitions M with |M| <= D
*
      DO WHILE( D .LE. MIN ( MAXORD, MAXRDM ) )
         PRTCNT = 0
         CALL NXPART( PRTCNT, S, M, MODOFM )
         DO WHILE( MODOFM .LE. D )
*     
****  Calculate updated weight for partition M alnd 
****     fully symmetric sums ( when necessary )
*
            FULWGT = WEIGHT( S, X, M, K, MODOFM, D )
            IF ( D .EQ. MODOFM ) THEN
               DO I = 1,N
                  FULSMS(I,PRTCNT) = 0
               END DO
               FULSMS(N+1,PRTCNT) = 0
            END IF
            IF ( FULSMS(N+1,PRTCNT) .EQ. 0 .AND. FULWGT .NE. 0 ) THEN
               CALL FULSMH( S, M, N, F, FULSMS(1,PRTCNT), X, WORK ) 
     &                      
               INTCLS = INTCLS + FULSMS(N+1,PRTCNT)
            END IF
            DO I = 1,N 
               INTVAL(I) = INTVAL(I) + FULWGT*FULSMS(I,PRTCNT)
            END DO
            CALL NXPART( PRTCNT, S, M, MODOFM )
         END DO
*     
****  End loop for each D
*
         D = D + 1
      END DO
      MINORD = MAXORD
      NUMSMS = PRTCNT - 1
      END
*
      SUBROUTINE FULSMH( S, M, N, F, FULSMS, X, FUNVAL )
*
****  To compute fully symmetric basic rule sums
*
      INTEGER S, M(*), SUMCLS, IX, LX, I,L, MI, ML, IL, N, MX
      DOUBLE PRECISION X(*), INTWGT, FULSMS(*), FUNVAL(*)
*
*        Generators for 1 + 2 + 6 + 10 + 16 = 35 point 
*        degree 51 rule, with degree 1, 5, 15 and 29 imbedded rules.     
*
      PARAMETER ( MX = 17 )
      DOUBLE PRECISION G(0:MX)
      SAVE G
      DATA G( 0),G( 1)/0D0                  , 0.17320508075688773D1/
      DATA G( 2),G( 3)/0.41849560176727319D1, 0.74109534999454084D0/
      DATA G( 4),G( 5)/0.28612795760570581D1, 0.63633944943363700D1/
      DATA G( 6),G( 7)/0.12304236340273060D1, 0.51870160399136561D1/
      DATA G( 8),G( 9)/0.25960831150492022D1, 0.32053337944991945D1/
      DATA G(10),G(11)/0.90169397898903025D1, 0.24899229757996061D0/
      DATA G(12),G(13)/0.79807717985905609D1, 0.22336260616769417D1/
      DATA G(14),G(15)/0.71221067008046167D1, 0.36353185190372782D1/
      DATA G(16),G(17)/0.56981777684881096D1, 0.47364330859522971D1/
      INTWGT = 1
      DO I = 1, S
        IF ( M(I) .NE. 0 ) INTWGT = INTWGT/2
      END DO
      SUMCLS = 0
      DO I = 1,N
         FULSMS(I) = 0
      END DO
*
********  Compute centrally symmetric sum for permutation M
*
 10   DO I = 1, S
         X(I) = -G(M(I))
      END DO
 20   SUMCLS = SUMCLS + 1
      CALL F( S, X, N, FUNVAL )
      DO I = 1,N
         FULSMS(I) = FULSMS(I) + INTWGT*FUNVAL(I)
      END DO
      DO I = 1, S
         X(I) = -X(I)
         IF ( X(I) .GT. 0 ) GO TO 20
      END DO
*
********  END Integration loop for M
*
*
********  Find next distinct permutation of M and loop back
*          to compute next centrally symmetric sum
*
      DO I = 2, S
         IF ( M(I-1) .GT. M(I) ) THEN
            MI = M(I)
            IX = I - 1
            IF ( I .GT. 2 ) THEN
               DO L = 1, IX/2
                  ML = M(L)
                  IL = I - L
                  M(L) = M(IL)
                  M(IL) = ML
                  IF ( ML .LE. MI ) IX = IX - 1
                  IF ( M(L) .GT. MI ) LX = L
               END DO
               IF ( M(IX) .LE. MI ) IX = LX
            END IF
            M(I) = M(IX)
            M(IX) = MI
            GO TO 10
         END IF
      END DO
*
****  END Loop for permutations of M and associated sums
*
*
**** Restore original order to M.
*
      DO I = 1, S/2
         MI = M(I)
         M(I) = M(S-I+1)
         M(S-I+1) = MI
      END DO
      FULSMS(N+1) = SUMCLS
      END
*
      DOUBLE PRECISION FUNCTION WEIGHT( S, INTRPS, M,K, MODOFM, D )
*
****  Subroutine to update weight for partition m
*
      INTEGER S, M(S), K(S), I, L, D, MAXRDM, NZRMAX, MODOFM
      PARAMETER ( MAXRDM = 25 )
      DOUBLE PRECISION INTRPS(S), MOMPRD(0:MAXRDM,0:MAXRDM), MOMNKN,
     &     INTMPA, INTMPB
      SAVE MOMPRD
*
*     Modified moments and generators for 1 + 2 + 6 + 10 + 16 = 35 point 
*        degree 51 rule, with degree 1, 5 and 19 imbedded rules.     
*
      PARAMETER ( NZRMAX = 17 )
      DOUBLE PRECISION A(0:NZRMAX), G(0:NZRMAX)
      DATA A /  2*1D0, 0D0, 6D0, 
     &     -0.48378475125832451D2, 3*0D0, 34020D0, 
     &     -0.98606453173677489D6, 5*0D0, 0.12912054173706603D13, 
     &     -0.11268664521456168D15, 0.29248520348796280D16 /  
      DATA G( 0),G( 1)/0D0                  , 0.17320508075688773D1/
      DATA G( 2),G( 3)/0.41849560176727319D1, 0.74109534999454084D0/
      DATA G( 4),G( 5)/0.28612795760570581D1, 0.63633944943363700D1/
      DATA G( 6),G( 7)/0.12304236340273060D1, 0.51870160399136561D1/
      DATA G( 8),G( 9)/0.25960831150492022D1, 0.32053337944991945D1/
      DATA G(10),G(11)/0.90169397898903025D1, 0.24899229757996061D0/
      DATA G(12),G(13)/0.79807717985905609D1, 0.22336260616769417D1/
      DATA G(14),G(15)/0.71221067008046167D1, 0.36353185190372782D1/
      DATA G(16),G(17)/0.56981777684881096D1, 0.47364330859522971D1/
      DATA MOMPRD(0,0) / 0D0 /
      IF ( MOMPRD(0,0) .EQ. 0 ) THEN
*
****  Calculate moments 
*
         DO L = 0, MAXRDM 
            DO I = 0, MAXRDM 
               MOMPRD(L,I) = 0
            END DO
         END DO
         MOMPRD(0,0) = A(0)
         DO L = 0, NZRMAX
            MOMNKN = 1
            DO I = 1, NZRMAX
               IF ( I .LE. L ) THEN
                  MOMNKN = MOMNKN*( G(L)**2 - G(I-1)**2 )
               ELSE
                  MOMNKN = MOMNKN*( G(L)**2 - G(I)**2 )
               END IF
               IF ( I .GE. L ) MOMPRD(L,I) = A(I)/MOMNKN
            END DO
         END DO
      END IF
*
*     Determine Updated Weight Contribution
*
      DO I = 2,S
        INTRPS(I) = 0
        K(I) = M(I)
      END DO
      K(1) = D - MODOFM + M(1)
 10   INTRPS(1) = MOMPRD( M(1), K(1) )
      DO I = 2, S
        INTRPS(I) = INTRPS(I) + MOMPRD( M(I), K(I) )*INTRPS(I-1)
        INTRPS(I-1) = 0
        K(1) = K(1) - 1
        K(I) = K(I) + 1
        IF ( K(1) .GE. M(1) ) GO TO 10
        K(1) = K(1) + K(I) - M(I)
        K(I) = M(I)
      END DO
      WEIGHT = INTRPS(S)
      END
      SUBROUTINE NXPART( PRTCNT, S, M, MODOFM )
*
****  SUBROUTINE TO COMPUTE THE NEXT S PARTITION
*
      INTEGER S, M(S), PRTCNT, MODOFM, I, MSUM, L
      IF ( PRTCNT .EQ. 0 ) THEN
         DO I = 1, S
            M(I) = 0
         END DO
         PRTCNT = 1
         MODOFM = 0
      ELSE
         PRTCNT = PRTCNT + 1
         MSUM = M(1)
         DO I = 2, S
            MSUM = MSUM + M(I)
            IF ( M(1) .LE. M(I) + 1 ) THEN
               M(I) = 0
            ELSE
               M(1) = MSUM - (I-1)*(M(I)+1)
               DO L = 2, I
                  M(L) = M(I) + 1
               END DO
               RETURN
            END IF
         END DO
         M(1) = MSUM + 1
         MODOFM = M(1)
      END IF
      END
*
      SUBROUTINE MLTRUL( NDIM, NUMFUN, FUNSUB, NP, POINT, WEIGHT, 
     &                   INTVAL, FUNS, X, IC )
*
*     Computes product integration rule
*
      EXTERNAL FUNSUB
      INTEGER NDIM, NP, NUMFUN, ICI, I
      DOUBLE PRECISION POINT(*), WEIGHT(*), INTVAL(*)
      DOUBLE PRECISION WTPROD, X(*), FUNS(*), IC(*)
      DO I = 1, NDIM
         IC(I) = 1
      END DO
      DO I = 1, NUMFUN
         INTVAL(I) = 0
      END DO
 10   WTPROD = 1
      DO I = 1, NDIM
         ICI = IC(I)
         X(I) = POINT(ICI)
         WTPROD = WTPROD*WEIGHT(ICI)
      END DO
      CALL FUNSUB( NDIM, X, NUMFUN, FUNS ) 
      DO I = 1, NUMFUN
         INTVAL(I) = INTVAL(I) + WTPROD*FUNS(I)
      END DO
      DO I = 1, NDIM
         IC(I) = IC(I) + 1
         IF ( IC(I) .LE. NP ) GO TO 10
         IC(I) = 1
      END DO
      END
*
      SUBROUTINE HERMIT( NDIM, NUMFUN, MINPTS, MAXPTS, FUNSUB,
     &     EPSABS, EPSREL, RESTAR, RESULT, ABSERR, NEVAL, IFAIL, WORK ) 
****BEGIN PROLOGUE HERMIT
****CATEGORY NO. H2B1A1
****AUTHOR
*
*     Alan Genz
*     Department of Mathematics
*     Washington State University
*     Pullman, Washington 99164-3113  USA
*     Telephone: 509-335-2131 
*     Electronic Mail: alangenz@wsu.edu
*     Fax: 509-335-1188 
*
****KEYWORDS automatic multidimensional integrator, Gaussian weight
*            n-dimensional region (-infin, infin)^n 
****PURPOSE  The routine calculates an approximation to a vector of
*            definite integrals using Gauss-Hermite product rules.
*
*      infin     infin
*     I    ...  I       w(X) (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
*     -infin    -infin         1  2      NUMFUN
*
*       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN,
*              I   I  1  2      NDIM
*
*       w(X) = EXP(-( X(1)**2 + ... + X(NDIM)**2 )/2)/SQRT(2*PI)**NDIM.
*
****DESCRIPTION Computation of integrals over infinite regions with
*               Gaussian weight using Gauss-Hermite product rules.
*
*   ON ENTRY
*
*     NDIM   Integer, number of variables.
*     NUMFUN Integer, number of components of the integral.
*     MINPTS Integer, minimum number of function evaluations.
*     MAXPTS Integer, maximum number of function evaluations.
*     FUNSUB Externally declared subroutine for computing components of 
*            the integrand at the given evaluation point.
*            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
*            Input parameters:
*              NDIM   Integer, number of dimensions for the integral.
*              X      Real NDIM-array for the evaluation point.
*              NUMFUN Integer, number of components of I.
*            Output parameter:
*              FUNVLS Real NUMFUN-array for integrand components.
*
*     EPSABS Real, requested absolute accuracy.
*     EPSREL Real, requested relative accuracy.
*     RESTAR Integer.
*            If RESTAR = 0, this is the first attempt to compute
*            the integral.
*            If RESTAR = 1, then we restart a previous attempt.
*     WORK   Real array of working storage, must have dimensions at
*            least 2*NUMFUN+2*NDIM
*
*   ON RETURN
*
*     RESULT Real NUMFUN-array of approximations to all components 
*            of the integral.
*     ABSERR Real NUMFUN-array of estimates of absolute accuracies.
*     NEVAL  Integer, number of function evaluations used.
*     IFAIL  Integer.
*            IFAIL = 0 for normal exit, when 
*              ABSERR(K) <= MAX( EPSABS, ABS(RESULT(K))*EPSREL ) for
*              all K, 0 < K <= NUMFUN, with <= MAXPTS function values. 
*            IFAIL = 1 if MAXPTS was too small to obtain the required 
*              accuracy. In this case values of RESULT are returned
*              with estimated absolute accuracies ABSERR.
*
****ROUTINES CALLED MLTRUL
****END PROLOGUE HERMIT
*
*   Global variables.
*
      EXTERNAL FUNSUB
      INTEGER NDIM, NUMFUN, MINPTS, MAXPTS, RESTAR
      INTEGER NEVAL, IFAIL
      DOUBLE PRECISION EPSABS, EPSREL
      DOUBLE PRECISION RESULT(*), ABSERR(*), WORK(*)
*
*   Local variables.
*
      INTEGER I, K, MAXRUL, RULE
      PARAMETER ( MAXRUL = 50 )
      DOUBLE PRECISION POINT(MAXRUL), WEIGHT(MAXRUL)
      SAVE RULE
      DOUBLE PRECISION W(25,50), T(25,50)
      SAVE W, T
      IF ( RESTAR .EQ. 0 ) RULE = 1
      NEVAL = 0
 10   IF ( NEVAL + RULE**NDIM .LE. MAXPTS .AND. RULE .LT. MAXRUL ) THEN
         DO I = 1, RULE/2
            POINT( I) = -T(I,RULE)
            WEIGHT(I) =  W(I,RULE)
            POINT( RULE-I+1) = T(I,RULE)
            WEIGHT(RULE-I+1) = W(I,RULE)
         END DO
         IF ( MOD( RULE, 2 ) .EQ. 1 ) THEN
            POINT(  RULE/2 + 1 ) = 0
            WEIGHT( RULE/2 + 1 ) = W( RULE/2 + 1, RULE ) 
         END IF
         CALL MLTRUL( NDIM, NUMFUN, FUNSUB, RULE, POINT, WEIGHT, 
     &        RESULT, WORK, WORK(NUMFUN+1), WORK(NUMFUN+NDIM+1) )
         NEVAL = NEVAL + RULE**NDIM
         IFAIL = 0
         DO I = 1,NUMFUN
            IF ( RULE .GT. 1 ) THEN
               ABSERR(I) = ABS( RESULT(I) - WORK(2*NDIM+NUMFUN+I)  )
            ELSE
               ABSERR(I) = ABS( RESULT(I) )
            ENDIF
            WORK(2*NDIM+NUMFUN+I) = RESULT(I)
            IF ( ABSERR(I) .GE. MAX( EPSABS, EPSREL*ABS(RESULT(I)) ) 
     &           .OR. RULE .EQ. 1 ) IFAIL = 1
         END DO
         RULE = RULE + 1
         IF ( IFAIL .GT. 0 .OR. NEVAL .LT. MINPTS ) GO TO 10 
      END IF
*
*     Gauss Hermite Weights and Points, N = 1,50
*
      DATA   W(1, 1), T(1, 1) / 1D0, 0D0 /
      DATA ( W(I, 2), T(I, 2), I = 1, 1) /
     &  0.5000000000000001D+00, 0.1000000000000000D+01/
      DATA ( W(I, 3), T(I, 3), I = 1, 2) /
     &  0.1666666666666667D+00, 0.1732050807568877D+01,
     &  0.6666666666666664D+00, 0.1107367643833737D-15/
      DATA ( W(I, 4), T(I, 4), I = 1, 2) /
     &  0.4587585476806855D-01, 0.2334414218338977D+01,
     &  0.4541241452319317D+00, 0.7419637843027258D+00/
      DATA ( W(I, 5), T(I, 5), I = 1, 3) /
     &  0.1125741132772071D-01, 0.2856970013872805D+01,
     &  0.2220759220056126D+00, 0.1355626179974265D+01,
     &  0.5333333333333342D+00, 0.9386691848789097D-16/
      DATA ( W(I, 6), T(I, 6), I = 1, 3) /
     &  0.2555784402056243D-02, 0.3324257433552119D+01,
     &  0.8861574604191447D-01, 0.1889175877753710D+01,
     &  0.4088284695560291D+00, 0.6167065901925933D+00/
      DATA ( W(I, 7), T(I, 7), I = 1, 4) /
     &  0.5482688559722184D-03, 0.3750439717725742D+01,
     &  0.3075712396758645D-01, 0.2366759410734542D+01,
     &  0.2401231786050126D+00, 0.1154405394739968D+01,
     &  0.4571428571428575D+00, 0.2669848554723344D-16/
      DATA ( W(I, 8), T(I, 8), I = 1, 4) /
     &  0.1126145383753679D-03, 0.4144547186125893D+01,
     &  0.9635220120788268D-02, 0.2802485861287542D+01,
     &  0.1172399076617590D+00, 0.1636519042435109D+01,
     &  0.3730122576790775D+00, 0.5390798113513754D+00/
      DATA ( W(I, 9), T(I, 9), I = 1, 5) /
     &  0.2234584400774664D-04, 0.4512745863399781D+01,
     &  0.2789141321231769D-02, 0.3205429002856470D+01,
     &  0.4991640676521780D-01, 0.2076847978677829D+01,
     &  0.2440975028949394D+00, 0.1023255663789133D+01,
     &  0.4063492063492066D+00, 0.0000000000000000D+00/
      DATA ( W(I,10), T(I,10), I = 1, 5) /
     &  0.4310652630718282D-05, 0.4859462828332311D+01,
     &  0.7580709343122187D-03, 0.3581823483551927D+01,
     &  0.1911158050077027D-01, 0.2484325841638954D+01,
     &  0.1354837029802680D+00, 0.1465989094391158D+01,
     &  0.3446423349320194D+00, 0.4849357075154977D+00/
      DATA ( W(I,11), T(I,11), I = 1, 6) /
     &  0.8121849790214922D-06, 0.5188001224374871D+01,
     &  0.1956719302712241D-03, 0.3936166607129977D+01,
     &  0.6720285235537304D-02, 0.2865123160643646D+01,
     &  0.6613874607105794D-01, 0.1876035020154847D+01,
     &  0.2422402998739701D+00, 0.9288689973810635D+00,
     &  0.3694083694083690D+00, 0.0000000000000000D+00/
      DATA ( W(I,12), T(I,12), I = 1, 6) /
     &  0.1499927167637166D-06, 0.5500901704467746D+01,
     &  0.4837184922590630D-04, 0.4271825847932281D+01,
     &  0.2203380687533207D-02, 0.3223709828770096D+01,
     &  0.2911668791236414D-01, 0.2259464451000800D+01,
     &  0.1469670480453302D+00, 0.1340375197151617D+01,
     &  0.3216643615128298D+00, 0.4444030019441390D+00/
      DATA ( W(I,13), T(I,13), I = 1, 7) /
     &  0.2722627642805887D-07, 0.5800167252386502D+01,
     &  0.1152659652733391D-04, 0.4591398448936520D+01,
     &  0.6812363504429268D-03, 0.3563444380281636D+01,
     &  0.1177056050599653D-01, 0.2620689973432215D+01,
     &  0.7916895586044999D-01, 0.1725418379588239D+01,
     &  0.2378715229641365D+00, 0.8566794935194499D+00,
     &  0.3409923409923412D+00, 0.2011511664336819D-15/
      DATA ( W(I,14), T(I,14), I = 1, 7) /
     &  0.4868161257748367D-08, 0.6087409546901291D+01,
     &  0.2660991344067620D-05, 0.4896936397345567D+01,
     &  0.2003395537607445D-03, 0.3886924575059772D+01,
     &  0.4428919106947401D-02, 0.2963036579838668D+01,
     &  0.3865010882425336D-01, 0.2088344745701943D+01,
     &  0.1540833398425136D+00, 0.1242688955485464D+01,
     &  0.3026346268130198D+00, 0.4125904579546022D+00/
      DATA ( W(I,15), T(I,15), I = 1, 8) /
     &  0.8589649899633300D-09, 0.6363947888829836D+01,
     &  0.5975419597920602D-06, 0.5190093591304780D+01,
     &  0.5642146405189029D-04, 0.4196207711269018D+01,
     &  0.1567357503549958D-02, 0.3289082424398766D+01,
     &  0.1736577449213763D-01, 0.2432436827009758D+01,
     &  0.8941779539984458D-01, 0.1606710069028730D+01,
     &  0.2324622936097322D+00, 0.7991290683245483D+00,
     &  0.3182595182595181D+00, 0.0000000000000000D+00/
      DATA ( W(I,16), T(I,16), I = 1, 8) /
     &  0.1497814723161838D-09, 0.6630878198393126D+01,
     &  0.1309473216286842D-06, 0.5472225705949343D+01,
     &  0.1530003216248727D-04, 0.4492955302520013D+01,
     &  0.5259849265739089D-03, 0.3600873624171548D+01,
     &  0.7266937601184742D-02, 0.2760245047630703D+01,
     &  0.4728475235401395D-01, 0.1951980345716333D+01,
     &  0.1583383727509496D+00, 0.1163829100554964D+01,
     &  0.2865685212380120D+00, 0.3867606045005573D+00/
      DATA ( W(I,17), T(I,17), I = 1, 9) /
     &  0.2584314919374932D-10, 0.6889122439895331D+01,
     &  0.2808016117930569D-07, 0.5744460078659410D+01,
     &  0.4012679447979839D-05, 0.4778531589629983D+01,
     &  0.1684914315513387D-03, 0.3900065717198010D+01,
     &  0.2858946062284621D-02, 0.3073797175328194D+01,
     &  0.2308665702571097D-01, 0.2281019440252989D+01,
     &  0.9740637116272111D-01, 0.1509883307796740D+01,
     &  0.2267063084689769D+00, 0.7518426007038956D+00,
     &  0.2995383701266057D+00, 0.0000000000000000D+00/
      DATA ( W(I,18), T(I,18), I = 1, 9) /
     &  0.4416588769358736D-11, 0.7139464849146476D+01,
     &  0.5905488478836554D-08, 0.6007745911359599D+01,
     &  0.1021552397636983D-05, 0.5054072685442739D+01,
     &  0.5179896144116204D-04, 0.4188020231629400D+01,
     &  0.1065484796291652D-02, 0.3374736535778089D+01,
     &  0.1051651775194131D-01, 0.2595833688911239D+01,
     &  0.5489663248022256D-01, 0.1839779921508646D+01,
     &  0.1606853038935128D+00, 0.1098395518091501D+01,
     &  0.2727832346542882D+00, 0.3652457555076979D+00/
      DATA ( W(I,19), T(I,19), I = 1,10) /
     &  0.7482830054057162D-12, 0.7382579024030434D+01,
     &  0.1220370848447449D-08, 0.6262891156513252D+01,
     &  0.2532220032092866D-06, 0.5320536377336039D+01,
     &  0.1535114595466674D-04, 0.4465872626831029D+01,
     &  0.3785021094142701D-03, 0.3664416547450636D+01,
     &  0.4507235420342067D-02, 0.2898051276515753D+01,
     &  0.2866669103011841D-01, 0.2155502761316934D+01,
     &  0.1036036572761442D+00, 0.1428876676078373D+01,
     &  0.2209417121991433D+00, 0.7120850440423796D+00,
     &  0.2837731927515210D+00, 0.4118522463420039D-15/
      DATA ( W(I,20), T(I,20), I = 1,10) /
     &  0.1257800672437914D-12, 0.7619048541679760D+01,
     &  0.2482062362315163D-09, 0.6510590157013660D+01,
     &  0.6127490259983006D-07, 0.5578738805893195D+01,
     &  0.4402121090230841D-05, 0.4734581334046057D+01,
     &  0.1288262799619300D-03, 0.3943967350657311D+01,
     &  0.1830103131080496D-02, 0.3189014816553389D+01,
     &  0.1399783744710099D-01, 0.2458663611172367D+01,
     &  0.6150637206397690D-01, 0.1745247320814126D+01,
     &  0.1617393339840001D+00, 0.1042945348802752D+01,
     &  0.2607930634495551D+00, 0.3469641570813557D+00/
      DATA ( W(I,21), T(I,21), I = 1,11) /
     &  0.2098991219565665D-13, 0.7849382895113822D+01,
     &  0.4975368604121770D-10, 0.6751444718717456D+01,
     &  0.1450661284493093D-07, 0.5829382007304472D+01,
     &  0.1225354836148259D-05, 0.4994963944782024D+01,
     &  0.4219234742551696D-04, 0.4214343981688420D+01,
     &  0.7080477954815349D-03, 0.3469846690475375D+01,
     &  0.6439697051408779D-02, 0.2750592981052372D+01,
     &  0.3395272978654278D-01, 0.2049102468257161D+01,
     &  0.1083922856264195D+00, 0.1359765823211230D+01,
     &  0.2153337156950595D+00, 0.6780456924406435D+00,
     &  0.2702601835728773D+00, 0.0000000000000000D+00/
      DATA ( W(I,22), T(I,22), I = 1,11) /
     &  0.3479460647877136D-14, 0.8074029984021710D+01,
     &  0.9841378982346056D-11, 0.6985980424018808D+01,
     &  0.3366514159458310D-08, 0.6073074951122888D+01,
     &  0.3319853749814059D-06, 0.5247724433714421D+01,
     &  0.1334597712680954D-04, 0.4476361977310866D+01,
     &  0.2622833032559635D-03, 0.3741496350266517D+01,
     &  0.2808761047577212D-02, 0.3032404227831676D+01,
     &  0.1756907288080571D-01, 0.2341759996287707D+01,
     &  0.6719631142889003D-01, 0.1664124839117906D+01,
     &  0.1619062934136754D+00, 0.9951624222712152D+00,
     &  0.2502435965869353D+00, 0.3311793157152742D+00/
      DATA ( W(I,23), T(I,23), I = 1,12) /
     &  0.5732383167802038D-15, 0.8293386027417354D+01,
     &  0.1922935311567786D-11, 0.7214659435051866D+01,
     &  0.7670888862399855D-09, 0.6310349854448401D+01,
     &  0.8775062483861979D-07, 0.5493473986471793D+01,
     &  0.4089977244992140D-05, 0.4730724197451473D+01,
     &  0.9340818609031275D-04, 0.4004775321733304D+01,
     &  0.1167628637497855D-02, 0.3305040021752963D+01,
     &  0.8579678391465647D-02, 0.2624323634059181D+01,
     &  0.3886718370348111D-01, 0.1957327552933424D+01,
     &  0.1120733826026210D+00, 0.1299876468303978D+01,
     &  0.2099596695775429D+00, 0.6484711535344957D+00,
     &  0.2585097408088385D+00, 0.0000000000000000D+00/
      DATA ( W(I,24), T(I,24), I = 1,12) /
     &  0.9390193689041782D-16, 0.8507803519195264D+01,
     &  0.3714974152762395D-12, 0.7437890666021664D+01,
     &  0.1718664927964866D-09, 0.6541675005098631D+01,
     &  0.2267461673480609D-07, 0.5732747175251204D+01,
     &  0.1217659745442582D-05, 0.4978041374639117D+01,
     &  0.3209500565274598D-04, 0.4260383605019904D+01,
     &  0.4647187187793975D-03, 0.3569306764073560D+01,
     &  0.3976608929181313D-02, 0.2897728643223314D+01,
     &  0.2112634440896754D-01, 0.2240467851691752D+01,
     &  0.7206936401717838D-01, 0.1593480429816420D+01,
     &  0.1614595128670001D+00, 0.9534219229321088D+00,
     &  0.2408701155466405D+00, 0.3173700966294525D+00/
      DATA ( W(I,25), T(I,25), I = 1,13) /
     &  0.1530038997998690D-16, 0.8717597678399592D+01,
     &  0.7102103037003980D-13, 0.7656037955393078D+01,
     &  0.3791150000477161D-10, 0.6767464963809719D+01,
     &  0.5738023868899356D-08, 0.5966014690606704D+01,
     &  0.3530152560245470D-06, 0.5218848093644280D+01,
     &  0.1067219490520254D-04, 0.4508929922967284D+01,
     &  0.1777669069265268D-03, 0.3825900569972490D+01,
     &  0.1757850405263803D-02, 0.3162775679388193D+01,
     &  0.1085675599146230D-01, 0.2514473303952205D+01,
     &  0.4337997016764489D-01, 0.1877058369947839D+01,
     &  0.1148809243039517D+00, 0.1247311975616789D+01,
     &  0.2048510256503405D+00, 0.6224622791860757D+00,
     &  0.2481693511764858D+00, 0.0000000000000000D+00/
      DATA ( W(I,26), T(I,26), I = 1,13) /
     &  0.2480694260393664D-17, 0.8923051727828243D+01,
     &  0.1344547649663596D-13, 0.7869426697637738D+01,
     &  0.8242809443163844D-11, 0.6988088770623415D+01,
     &  0.1424293237988014D-08, 0.6193693483796317D+01,
     &  0.9986755573314568D-07, 0.5453615383857833D+01,
     &  0.3443413612308114D-05, 0.4750947483085378D+01,
     &  0.6557558694333818D-04, 0.4075427214412228D+01,
     &  0.7442025763604303D-03, 0.3420156373999979D+01,
     &  0.5302198015682246D-02, 0.2780138499509748D+01,
     &  0.2459766565712125D-01, 0.2151530090121648D+01,
     &  0.7622953220630281D-01, 0.1531215708695402D+01,
     &  0.1605865456137948D+00, 0.9165450413386282D+00,
     &  0.2324707356300776D+00, 0.3051559707592978D+00/
      DATA ( W(I,27), T(I,27), I = 1,14) /
     &  0.4003364766550257D-18, 0.9124421250672931D+01,
     &  0.2522363250873417D-14, 0.8078349274534165D+01,
     &  0.1768236511219616D-11, 0.7203876611910644D+01,
     &  0.3472604702845840D-09, 0.6416154934562174D+01,
     &  0.2761933918447925D-07, 0.5682760761629052D+01,
     &  0.1080581533683832D-05, 0.4986906410679802D+01,
     &  0.2339557671566820D-04, 0.4318417671936682D+01,
     &  0.3028398259361645D-03, 0.3670472986492407D+01,
     &  0.2471872445961970D-02, 0.3038150251871036D+01,
     &  0.1321102584046355D-01, 0.2417683983162542D+01,
     &  0.4748957556274264D-01, 0.1806045213138672D+01,
     &  0.1169962651750102D+00, 0.1200683354549981D+01,
     &  0.2000149701605136D+00, 0.5993548807899852D+00,
     &  0.2389778937255043D+00, 0.0000000000000000D+00/
      DATA ( W(I,28), T(I,28), I = 1,14) /
     &  0.6432547438801930D-19, 0.9321937814408766D+01,
     &  0.4691765569500354D-15, 0.8283069540861421D+01,
     &  0.3745901035176660D-12, 0.7415125286176065D+01,
     &  0.8326609843882241D-10, 0.6633731493950435D+01,
     &  0.7479362584613589D-08, 0.5906656325824994D+01,
     &  0.3304864449926482D-06, 0.5217223673447450D+01,
     &  0.8093584057145153D-05, 0.4555340384596974D+01,
     &  0.1188285381401779D-03, 0.3914253725963635D+01,
     &  0.1104305927857598D-02, 0.3289106970171833D+01,
     &  0.6752459709030160D-02, 0.2676201879526944D+01,
     &  0.2793578476788097D-01, 0.2072582674144621D+01,
     &  0.7977336601159966D-01, 0.1475781736957922D+01,
     &  0.1594181936613094D+00, 0.8836525629929802D+00,
     &  0.2248886297506769D+00, 0.2942517144887133D+00/
      DATA ( W(I,29), T(I,29), I = 1,15) /
     &  0.1029341808721942D-19, 0.9515812006947357D+01,
     &  0.8657491667957282D-16, 0.8483826557846555D+01,
     &  0.7842840425658472D-13, 0.7622102722480985D+01,
     &  0.1965709944734762D-10, 0.6846722135707994D+01,
     &  0.1986123546067053D-08, 0.6125635348243716D+01,
     &  0.9868968543560684D-07, 0.5442271089178505D+01,
     &  0.2721127828058099D-05, 0.4786611062352805D+01,
     &  0.4508394026980976D-04, 0.4151964855100983D+01,
     &  0.4743663738893483D-03, 0.3533533770990993D+01,
     &  0.3297972210833669D-02, 0.2927678154322886D+01,
     &  0.1559400577786720D-01, 0.2331504884065565D+01,
     &  0.5121083528871912D-01, 0.1742616232610662D+01,
     &  0.1185603192669036D+00, 0.1158946149400189D+01,
     &  0.1954459569679010D+00, 0.5786461780331649D+00,
     &  0.2307372767004873D+00, 0.0000000000000000D+00/
      DATA ( W(I,30), T(I,30), I = 1,15) /
     &  0.1640807008117853D-20, 0.9706235997359524D+01,
     &  0.1585560944966296D-16, 0.8680837722732207D+01,
     &  0.1624080129972436D-13, 0.7825051744352813D+01,
     &  0.4573425871326147D-11, 0.7055396866960296D+01,
     &  0.5178459467189710D-09, 0.6339997686869597D+01,
     &  0.2882175154047618D-07, 0.5662381850082873D+01,
     &  0.8909088868621158D-06, 0.5012600596486518D+01,
     &  0.1657998163067346D-04, 0.4384020365898051D+01,
     &  0.1965129439848249D-03, 0.3771894423159236D+01,
     &  0.1544707339866097D-02, 0.3172634639420402D+01,
     &  0.8295747557723240D-02, 0.2583402100229274D+01,
     &  0.3111177018350134D-01, 0.2001858612956431D+01,
     &  0.8278683671562172D-01, 0.1426005658374115D+01,
     &  0.1580469532090208D+00, 0.8540733517109733D+00,
     &  0.2179999718155776D+00, 0.2844387607362094D+00/
      DATA ( W(I,31), T(I,31), I = 1,16) /
     &  0.2605973854893011D-21, 0.9893385708986649D+01,
     &  0.2883352367857899D-17, 0.8874301409488794D+01,
     &  0.3328468324148409D-14, 0.8024193227361653D+01,
     &  0.1049603362311349D-11, 0.7260000488890867D+01,
     &  0.1327251483589731D-09, 0.6550014268765684D+01,
     &  0.8243931619119761D-08, 0.5877855885986261D+01,
     &  0.2845610088162858D-06, 0.5233641511712708D+01,
     &  0.5923202317686233D-05, 0.4610789797323995D+01,
     &  0.7871624069602249D-04, 0.4004600901491224D+01,
     &  0.6960312713792868D-03, 0.3411532415843158D+01,
     &  0.4221717767270697D-02, 0.2828792768157509D+01,
     &  0.1796787584344161D-01, 0.2254095000754410D+01,
     &  0.5456725889447496D-01, 0.1685497905069052D+01,
     &  0.1196831096958545D+00, 0.1121297374047009D+01,
     &  0.1911320047746435D+00, 0.5599475878410030D+00,
     &  0.2232941387424060D+00, 0.9191380905810332D-16/
      DATA ( W(I,32), T(I,32), I = 1,16) /
     &  0.4124607489018384D-22, 0.1007742267422945D+02,
     &  0.5208449591960853D-18, 0.9064399210702408D+01,
     &  0.6755290223670036D-15, 0.8219728765382246D+01,
     &  0.2378064855777808D-12, 0.7460755754121516D+01,
     &  0.3347501239801238D-10, 0.6755930830540704D+01,
     &  0.2312518412074224D-08, 0.6088964309076983D+01,
     &  0.8881290713105934D-07, 0.5450033273623426D+01,
     &  0.2059622103953437D-05, 0.4832604613244488D+01,
     &  0.3055980306089618D-04, 0.4232021109995410D+01,
     &  0.3025570258170642D-03, 0.3644781249880835D+01,
     &  0.2062051051307883D-02, 0.3068135169013122D+01,
     &  0.9903461702320572D-02, 0.2499840415187396D+01,
     &  0.3410984772609194D-01, 0.1938004905925718D+01,
     &  0.8534480827208071D-01, 0.1380980199272144D+01,
     &  0.1565389937575984D+00, 0.8272849037797656D+00,
     &  0.2117055698804795D+00, 0.2755464192302757D+00/
      DATA ( W(I,33), T(I,33), I = 1,17) /
     &  0.6506889970402893D-23, 0.1025849562613868D+02,
     &  0.9349155921250728D-19, 0.9251297851734609D+01,
     &  0.1358447599822066D-15, 0.8411842935668213D+01,
     &  0.5323023871225154D-13, 0.7657866034784416D+01,
     &  0.8316046910555386D-11, 0.6957971061087896D+01,
     &  0.6369261724538402D-09, 0.6295953125159368D+01,
     &  0.2712480030928315D-07, 0.5662046690089215D+01,
     &  0.6982937054010693D-06, 0.5049763451908822D+01,
     &  0.1152282979948381D-04, 0.4454485185983318D+01,
     &  0.1271924628565943D-03, 0.3872747224621657D+01,
     &  0.9695342064240463D-03, 0.3301836743259241D+01,
     &  0.5227605765843533D-02, 0.2739550282026145D+01,
     &  0.2030404475704196D-01, 0.2184038256077331D+01,
     &  0.5758631173578715D-01, 0.1633699795932273D+01,
     &  0.1204510521056565D+00, 0.1087107916669903D+01,
     &  0.1870581852279859D+00, 0.5429533766656418D+00,
     &  0.2165276496896071D+00, 0.0000000000000000D+00/
      DATA ( W(I,34), T(I,34), I = 1,17) /
     &  0.1023327129805438D-23, 0.1043674187100505D+02,
     &  0.1668144375578546D-19, 0.9435150833760799D+01,
     &  0.2708054709262291D-16, 0.8600705233983431D+01,
     &  0.1177937790622538D-13, 0.7851517590699072D+01,
     &  0.2036657679770825D-11, 0.7156339259134842D+01,
     &  0.1724305566745258D-09, 0.6499046353927969D+01,
     &  0.8117409040122166D-08, 0.5869927588596779D+01,
     &  0.2312034264322868D-06, 0.5262536481734633D+01,
     &  0.4227725748387717D-05, 0.4672290691400979D+01,
     &  0.5182712643366873D-04, 0.4095758971954162D+01,
     &  0.4399649667746255D-03, 0.3530261634074169D+01,
     &  0.2650824238310194D-02, 0.2973629650303989D+01,
     &  0.1155073894677711D-01, 0.2424050509756231D+01,
     &  0.3692346295804431D-01, 0.1879964366420186D+01,
     &  0.8751168135862399D-01, 0.1339990625522619D+01,
     &  0.1549420965148600D+00, 0.8028734607837125D+00,
     &  0.2059249366691133D+00, 0.2674391839515571D+00/
      DATA ( W(I,35), T(I,35), I = 1,18) /
     &  0.1604619191790137D-24, 0.1061228847764259D+02,
     &  0.2959542628709020D-20, 0.9616099851106188D+01,
     &  0.5354094198066113D-17, 0.8786471736571588D+01,
     &  0.2578597004420350D-14, 0.8041881508402966D+01,
     &  0.4921158318497817D-12, 0.7351222593955577D+01,
     &  0.4592895709255784D-10, 0.6698448669188526D+01,
     &  0.2383123333705684D-08, 0.6073899909835108D+01,
     &  0.7486460178518943D-07, 0.5471169042499284D+01,
     &  0.1511924078161164D-05, 0.4885706922446520D+01,
     &  0.2050985370117561D-04, 0.4314112806712937D+01,
     &  0.1931447146139865D-03, 0.3753736849678044D+01,
     &  0.1294831077498348D-02, 0.3202440651527510D+01,
     &  0.6300195959720374D-02, 0.2658444295466304D+01,
     &  0.2258156121393648D-01, 0.2120223604836576D+01,
     &  0.6029658086774051D-01, 0.1586437891088980D+01,
     &  0.1209324521970301D+00, 0.1055876792225099D+01,
     &  0.1832085621911518D+00, 0.5274192342262778D+00,
     &  0.2103411454127607D+00, 0.2220535009031490D-16/
      DATA ( W(I,36), T(I,36), I = 1,18) /
     &  0.2509037634634927D-25, 0.1078525331238753D+02,
     &  0.5222366200862934D-21, 0.9794276019583023D+01,
     &  0.1050294193474738D-17, 0.8969286534562617D+01,
     &  0.5587113978649920D-15, 0.8229115367471579D+01,
     &  0.1174029255105150D-12, 0.7542793039211395D+01,
     &  0.1204744586548017D-10, 0.6894347646173911D+01,
     &  0.6871083387212228D-09, 0.6274168326809511D+01,
     &  0.2373822777365743D-07, 0.5675884710106672D+01,
     &  0.5278315192800947D-06, 0.5094978513857614D+01,
     &  0.7896978086723001D-05, 0.4528076990600175D+01,
     &  0.8220025562410442D-04, 0.3972557341929988D+01,
     &  0.6107548335511020D-03, 0.3426308595129129D+01,
     &  0.3304134538435289D-02, 0.2887579695004719D+01,
     &  0.1321657401560215D-01, 0.2354877715992540D+01,
     &  0.3955236976655980D-01, 0.1826896577986744D+01,
     &  0.8934249750438648D-01, 0.1302464954480165D+01,
     &  0.1532910133997116D+00, 0.7805064920524665D+00,
     &  0.2005920064390222D+00, 0.2600079252490002D+00/
      DATA ( W(I,37), T(I,37), I = 1,19) /
     &  0.3912701900272739D-26, 0.1095574594356165D+02,
     &  0.9167997950194993D-22, 0.9969800945691102D+01,
     &  0.2045033189096993D-18, 0.9149282977696753D+01,
     &  0.1198838092763842D-15, 0.8413364679488879D+01,
     &  0.2767202859289950D-13, 0.7731209035776908D+01,
     &  0.3114550633754001D-11, 0.7086915685017345D+01,
     &  0.1947529536776704D-09, 0.6470920475390299D+01,
     &  0.7379429831202629D-08, 0.5876887892594201D+01,
     &  0.1801392729842884D-06, 0.5300328473991422D+01,
     &  0.2963204697508371D-05, 0.4737895299703953D+01,
     &  0.3397941877127311D-04, 0.4186990225627911D+01,
     &  0.2788064134714879D-03, 0.3645526993515297D+01,
     &  0.1670452623119233D-02, 0.3111780274991816D+01,
     &  0.7424836460756516D-02, 0.2584284718210775D+01,
     &  0.2478561700046269D-01, 0.2061764482625703D+01,
     &  0.6272612686101582D-01, 0.1543082026656558D+01,
     &  0.1211816464812620D+00, 0.1027199336691835D+01,
     &  0.1795672590244481D+00, 0.5131472568284307D+00,
     &  0.2046562495907946D+00, 0.0000000000000000D+00/
      DATA ( W(I,38), T(I,38), I = 1,19) /
     &  0.6086019568424894D-27, 0.1112386843494987D+02,
     &  0.1601586834974089D-22, 0.1014278766112967D+02,
     &  0.3953752210235847D-19, 0.9326584757395571D+01,
     &  0.2548653789376282D-16, 0.8594764136386107D+01,
     &  0.6447906524854284D-14, 0.7916616928480361D+01,
     &  0.7941720700519787D-12, 0.7276311665555140D+01,
     &  0.5431469647397207D-10, 0.6664328864414711D+01,
     &  0.2251484789660455D-08, 0.6074366042059564D+01,
     &  0.6017570329946914D-07, 0.5501960756811018D+01,
     &  0.1085199304765951D-05, 0.4943790029802265D+01,
     &  0.1366659848918500D-04, 0.4397278309932714D+01,
     &  0.1234207566124386D-03, 0.3860361738163258D+01,
     &  0.8159914125966156D-03, 0.3331338060544006D+01,
     &  0.4014342501780207D-02, 0.2808766413901917D+01,
     &  0.1488362222069118D-01, 0.2291397563624214D+01,
     &  0.4200051997416258D-01, 0.1778123425136563D+01,
     &  0.9088415555258827D-01, 0.1267939114780749D+01,
     &  0.1516111462601078D+00, 0.7599132481367392D+00,
     &  0.1956519870413640D+00, 0.2531636353359348D+00/
      DATA ( W(I,39), T(I,39), I = 1,20) /
     &  0.9443344575063092D-28, 0.1128971604447632D+02,
     &  0.2784787505225604D-23, 0.1031334144275699D+02,
     &  0.7592450542206611D-20, 0.9501306853782864D+01,
     &  0.5370701458462819D-17, 0.8773438698067656D+01,
     &  0.1486129587733307D-14, 0.8099152213152410D+01,
     &  0.1998726680623669D-12, 0.7462682377866057D+01,
     &  0.1491720010448932D-10, 0.6854552519790964D+01,
     &  0.6748646364787737D-09, 0.6268491549685291D+01,
     &  0.1969872920599359D-07, 0.5700062454271859D+01,
     &  0.3884118283713278D-06, 0.5145964544094186D+01,
     &  0.5356584901373566D-05, 0.4603643074288992D+01,
     &  0.5307742112417270D-04, 0.4071054595947238D+01,
     &  0.3859381697690359D-03, 0.3546517669768644D+01,
     &  0.2093837438880663D-02, 0.3028613314092798D+01,
     &  0.8588083029362219D-02, 0.2516115854339318D+01,
     &  0.2690624148395815D-01, 0.2007943065265498D+01,
     &  0.6490157458316884D-01, 0.1503118900228415D+01,
     &  0.1212421280631244D+00, 0.1000744572356018D+01,
     &  0.1761190277014499D+00, 0.4999751896072681D+00,
     &  0.1994086534474408D+00, 0.1536204102353874D-15/
      DATA ( W(I,40), T(I,40), I = 1,20) /
     &  0.1461839873869467D-28, 0.1145337784154873D+02,
     &  0.4820467940200524D-24, 0.1048156053467427D+02,
     &  0.1448609431551587D-20, 0.9673556366934033D+01,
     &  0.1122275206827074D-17, 0.8949504543855559D+01,
     &  0.3389853443248306D-15, 0.8278940623659475D+01,
     &  0.4968088529197761D-13, 0.7646163764541459D+01,
     &  0.4037638581695192D-11, 0.7041738406453829D+01,
     &  0.1989118526027766D-09, 0.6459423377583766D+01,
     &  0.6325897188548972D-08, 0.5894805675372016D+01,
     &  0.1360342421574886D-06, 0.5344605445720084D+01,
     &  0.2048897436081474D-05, 0.4806287192093873D+01,
     &  0.2221177143247582D-04, 0.4277826156362752D+01,
     &  0.1770729287992397D-03, 0.3757559776168985D+01,
     &  0.1055879016901825D-02, 0.3244088732999869D+01,
     &  0.4773544881823334D-02, 0.2736208340465433D+01,
     &  0.1653784414256937D-01, 0.2232859218634873D+01,
     &  0.4427455520227679D-01, 0.1733090590631720D+01,
     &  0.9217657917006089D-01, 0.1236032004799159D+01,
     &  0.1499211117635710D+00, 0.7408707252859313D+00,
     &  0.1910590096619904D+00, 0.2468328960227240D+00/
      DATA ( W(I,41), T(I,41), I = 1,21) /
     &  0.2257863956583089D-29, 0.1161493725433746D+02,
     &  0.8308558938782992D-25, 0.1064753678631932D+02,
     &  0.2746891228522292D-21, 0.9843433249157988D+01,
     &  0.2326384145587187D-18, 0.9123069907984480D+01,
     &  0.7655982291966812D-16, 0.8456099083269388D+01,
     &  0.1220334874202772D-13, 0.7826882004053867D+01,
     &  0.1077818394935909D-11, 0.7226022663732790D+01,
     &  0.5769853428092003D-10, 0.6647308470747191D+01,
     &  0.1994794756757339D-08, 0.6086349164878472D+01,
     &  0.4667347708107243D-07, 0.5539884440458126D+01,
     &  0.7658186077982435D-06, 0.5005396683404125D+01,
     &  0.9058608622433030D-05, 0.4480878331594004D+01,
     &  0.7894719319504627D-04, 0.3964684028033266D+01,
     &  0.5158014443431912D-03, 0.3455432217780992D+01,
     &  0.2561642428649777D-02, 0.2951937016381193D+01,
     &  0.9777902738208298D-02, 0.2453159345907049D+01,
     &  0.2893721174793441D-01, 0.1958170711977291D+01,
     &  0.6684765935446599D-01, 0.1466125457295967D+01,
     &  0.1211489170115104D+00, 0.9762387671800500D+00,
     &  0.1728495310506020D+00, 0.4877685693194347D+00,
     &  0.1945450277536008D+00, 0.2585532684499631D-15/
      DATA ( W(I,42), T(I,42), I = 1,21) /
     &  0.3479841758734498D-30, 0.1177447255645880D+02,
     &  0.1426197845863333D-25, 0.1081135621818894D+02,
     &  0.5178070329449428D-22, 0.1001103095231325D+02,
     &  0.4785541849652557D-19, 0.9294235815925036D+01,
     &  0.1712817711028008D-16, 0.8630736540442662D+01,
     &  0.2963878417982294D-14, 0.8004954459331870D+01,
     &  0.2839400066530831D-12, 0.7407531683161432D+01,
     &  0.1648408975918176D-10, 0.6832282984221202D+01,
     &  0.6182418956905496D-09, 0.6274839704458262D+01,
     &  0.1570405641380752D-07, 0.5731959941924930D+01,
     &  0.2800444926136691D-06, 0.5201142761234950D+01,
     &  0.3605310211304410D-05, 0.4680396489703305D+01,
     &  0.3425737772910867D-04, 0.4168091525525806D+01,
     &  0.2445266890869028D-03, 0.3662862441045986D+01,
     &  0.1329885905572134D-02, 0.3163540283549010D+01,
     &  0.5573924561218487D-02, 0.2669104116603810D+01,
     &  0.1816809011555160D-01, 0.2178645208762502D+01,
     &  0.4638274179115778D-01, 0.1691339732834912D+01,
     &  0.9325376860459540D-01, 0.1206427277827926D+01,
     &  0.1482345543992958D+00, 0.7231933449391704D+00,
     &  0.1867743488620199D+00, 0.2409545348665914D+00/
      DATA ( W(I,43), T(I,43), I = 1,22) /
     &  0.5352075224079761D-31, 0.1193205730105631D+02,
     &  0.2438509785586889D-26, 0.1097309952495776D+02,
     &  0.9705961174387400D-23, 0.1017643700187913D+02,
     &  0.9772188868027905D-20, 0.9463096735522324D+01,
     &  0.3797453851791020D-17, 0.8802954705722227D+01,
     &  0.7121219095950012D-15, 0.8180490511448660D+01,
     &  0.7386428797509359D-13, 0.7586383052570858D+01,
     &  0.4641647576106650D-11, 0.7014473354182730D+01,
     &  0.1884798249393598D-09, 0.6460413330895460D+01,
     &  0.5186705460818634D-08, 0.5920978461471519D+01,
     &  0.1003006675541783D-06, 0.5393683423125381D+01,
     &  0.1402088827949991D-05, 0.4876551284706595D+01,
     &  0.1448871182961238D-04, 0.4367966934806672D+01,
     &  0.1126824016383807D-03, 0.3866579655692208D+01,
     &  0.6691683916969021D-03, 0.3371235817148052D+01,
     &  0.3070012395787847D-02, 0.2880930777450376D+01,
     &  0.1098373316809114D-01, 0.2394773428842528D+01,
     &  0.3087516724563331D-01, 0.1911959274951712D+01,
     &  0.6858704943109335D-01, 0.1431749364451692D+01,
     &  0.1209303683096624D+00, 0.9534532756297385D+00,
     &  0.1697454597838771D+00, 0.4764148861084949D+00,
     &  0.1900207247825863D+00, 0.0000000000000000D+00/
      DATA ( W(I,44), T(I,44), I = 1,22) /
     &  0.8215250893174055D-32, 0.1208776070905845D+02,
     &  0.4153634860809548D-27, 0.1113284252424393D+02,
     &  0.1809481070701256D-23, 0.1033973350764328D+02,
     &  0.1981515599992759D-20, 0.9629741154666469D+01,
     &  0.8346676981526894D-18, 0.8972848703628870D+01,
     &  0.1693427110357361D-15, 0.8353592294951005D+01,
     &  0.1898518041004563D-13, 0.7762686386163553D+01,
     &  0.1289064815963923D-11, 0.7193997236488432D+01,
     &  0.5656565204102218D-10, 0.6643196399727595D+01,
     &  0.1683045805404500D-08, 0.6107075817029749D+01,
     &  0.3522080139402352D-07, 0.5583164829645622D+01,
     &  0.5334152162808288D-06, 0.5069500234589782D+01,
     &  0.5980506497871435D-05, 0.4564480302199131D+01,
     &  0.5055024343502252D-04, 0.4066767790495557D+01,
     &  0.3269029726329101D-03, 0.3575222999425195D+01,
     &  0.1636879726834301D-02, 0.3088855994142003D+01,
     &  0.6407990218344719D-02, 0.2606791462987426D+01,
     &  0.1976568083382755D-01, 0.2128242119185467D+01,
     &  0.4833422729538417D-01, 0.1652487989479395D+01,
     &  0.9414471679496565D-01, 0.1178859803328855D+01,
     &  0.1465614441714095D+00, 0.7067252381063451D+00,
     &  0.1827650568597313D+00, 0.2354771181719222D+00/
      DATA ( W(I,45), T(I,45), I = 1,23) /
     &  0.1258601266761721D-32, 0.1224164801738294D+02,
     &  0.7049452438028993D-28, 0.1129065655801674D+02,
     &  0.3355898215846036D-24, 0.1050099761933589D+02,
     &  0.3990921764567262D-21, 0.9794252095357150D+01,
     &  0.1819415421781025D-18, 0.9140507651218687D+01,
     &  0.3987388394988911D-16, 0.8524355348620958D+01,
     &  0.4823870286948576D-14, 0.7936544056928108D+01,
     &  0.3532976789571211D-12, 0.7370964332170609D+01,
     &  0.1672373578472200D-10, 0.6823306517591415D+01,
     &  0.5370107061670978D-09, 0.6290378188903940D+01,
     &  0.1213731971819947D-07, 0.5769722503692224D+01,
     &  0.1987368663904533D-06, 0.5259389088450044D+01,
     &  0.2412167676023244D-05, 0.4757788618871960D+01,
     &  0.2210680714129846D-04, 0.4263596248773437D+01,
     &  0.1552889103199588D-03, 0.3775684997351202D+01,
     &  0.8463565255385183D-03, 0.3293078264001083D+01,
     &  0.3614802525284821D-02, 0.2814914962782011D+01,
     &  0.1219644981036976D-01, 0.2340423205074869D+01,
     &  0.3271890357088906D-01, 0.1868899890086888D+01,
     &  0.7014033220973423D-01, 0.1399694432579061D+01,
     &  0.1206095596750258D+00, 0.9321954002025563D+00,
     &  0.1667945553649214D+00, 0.4658191781783683D+00,
     &  0.1857980420096402D+00, 0.0000000000000000D+00/
      DATA ( W(I,46), T(I,46), I = 1,23) /
     &  0.1924664627046275D-33, 0.1239378079201923D+02,
     &  0.1192246272132870D-28, 0.1144660885260523D+02,
     &  0.6192859698073460D-25, 0.1066030193428052D+02,
     &  0.7986127260533585D-22, 0.9956707572498630D+01,
     &  0.3934545090204791D-19, 0.9306015173119082D+01,
     &  0.9300353507700858D-17, 0.8692869193236842D+01,
     &  0.1212244415430154D-14, 0.8108051845054231D+01,
     &  0.9561363845435019D-13, 0.7545477116064648D+01,
     &  0.4874112963243041D-11, 0.7000853362533224D+01,
     &  0.1686106290702903D-09, 0.6471003045526279D+01,
     &  0.4108222480852658D-08, 0.5953482378238384D+01,
     &  0.7258481863429432D-07, 0.5446353016216582D+01,
     &  0.9517578349775918D-06, 0.4948037176754664D+01,
     &  0.9436598821520223D-05, 0.4457221460031489D+01,
     &  0.7183252612530643D-04, 0.3972790546235975D+01,
     &  0.4250585899081316D-03, 0.3493779976317873D+01,
     &  0.1975260370663049D-02, 0.3019341533379349D+01,
     &  0.7268729480993403D-02, 0.2548717171165183D+01,
     &  0.2132401414735712D-01, 0.2081218863808104D+01,
     &  0.5013852830109192D-01, 0.1616212619554741D+01,
     &  0.9487419222577377D-01, 0.1153105444663487D+01,
     &  0.1449090160374568D+00, 0.6913343908693306D+00,
     &  0.1790029030973513D+00, 0.2303570440689780D+00/
      DATA ( W(I,47), T(I,47), I = 1,24) /
     &  0.2937991918931770D-34, 0.1254421721021957D+02,
     &  0.2009633017503842D-29, 0.1160076284239976D+02,
     &  0.1137324853690828D-25, 0.1081771486308368D+02,
     &  0.1588168143815981D-22, 0.1011718100450289D+02,
     &  0.8443830093861547D-20, 0.9469449861345275D+01,
     &  0.2149644881290608D-17, 0.8859217846077781D+01,
     &  0.3014376026444397D-15, 0.8277299513807538D+01,
     &  0.2556494006562142D-13, 0.7717631482085797D+01,
     &  0.1401228890589135D-11, 0.7175939408289193D+01,
     &  0.5213279817493255D-10, 0.6649059958212412D+01,
     &  0.1366926136038253D-08, 0.6134561715294918D+01,
     &  0.2601160240819878D-07, 0.5630517648304252D+01,
     &  0.3677447240564392D-06, 0.5135360748767347D+01,
     &  0.3936411635330397D-05, 0.4647788224870149D+01,
     &  0.3239911617951789D-04, 0.4166695488690111D+01,
     &  0.2076565749155176D-03, 0.3691129181403769D+01,
     &  0.1047278482271871D-02, 0.3220252777152719D+01,
     &  0.4191754574650512D-02, 0.2753320728634848D+01,
     &  0.1340827900837288D-01, 0.2289658541588672D+01,
     &  0.3446881631958411D-01, 0.1828647032406020D+01,
     &  0.7152609366400937D-01, 0.1369709566513547D+01,
     &  0.1202053639683183D+00, 0.9123014234783601D+00,
     &  0.1639855806134269D+00, 0.4559006613049348D+00,
     &  0.1818448921796476D+00, 0.0000000000000000D+00/
      DATA ( W(I,48), T(I,48), I = 1,24) /
     &  0.4477155473876067D-35, 0.1269301231543958D+02,
     &  0.3376456143135106D-30, 0.1175317846161662D+02,
     &  0.2079058971418782D-26, 0.1097330095851353D+02,
     &  0.3139476644799137D-23, 0.1027574158173351D+02,
     &  0.1798854916237414D-20, 0.9630885686949055D+01,
     &  0.4925465109698832D-18, 0.9023480280417658D+01,
     &  0.7419993598065164D-16, 0.8444371322523445D+01,
     &  0.6756677274657091D-14, 0.7887517316500372D+01,
     &  0.3975805958471478D-12, 0.7348660565908911D+01,
     &  0.1588360981246298D-10, 0.6824651320745980D+01,
     &  0.4474287153436318D-09, 0.6313069914905501D+01,
     &  0.9154055746313672D-08, 0.5811999987702533D+01,
     &  0.1392791689559962D-06, 0.5319884620419087D+01,
     &  0.1606393695542611D-05, 0.4835430885878957D+01,
     &  0.1426609232424675D-04, 0.4357544108790293D+01,
     &  0.9881804917611217D-04, 0.3885281117248394D+01,
     &  0.5395865846271863D-03, 0.3417816048684776D+01,
     &  0.2343082111760394D-02, 0.2954414688880120D+01,
     &  0.8149696855351926D-02, 0.2494414743040025D+01,
     &  0.2283821210009132D-01, 0.2037210303259033D+01,
     &  0.5180518354956036D-01, 0.1582239319674964D+01,
     &  0.9546340056142952D-01, 0.1128973231509649D+01,
     &  0.1432824569931762D+00, 0.6769081422063483D+00,
     &  0.1754635418118663D+00, 0.2255570729801638D+00/
      DATA ( W(I,49), T(I,49), I = 1,25) /
     &  0.6811389123116583D-36, 0.1284021824817448D+02,
     &  0.5655218285268972D-31, 0.1190391240788757D+02,
     &  0.3783656387961836D-27, 0.1112712121198839D+02,
     &  0.6170433961110509D-24, 0.1043245459795100D+02,
     &  0.3805275010242669D-21, 0.9790392369513130D+01,
     &  0.1119149962220974D-18, 0.9185730837014976D+01,
     &  0.1808768646225968D-16, 0.8609346484896321D+01,
     &  0.1765991900170034D-14, 0.8055219008757504D+01,
     &  0.1113983833798787D-12, 0.7519106753949869D+01,
     &  0.4771647040124784D-11, 0.6997872987021261D+01,
     &  0.1441769692730330D-09, 0.6489109229832087D+01,
     &  0.3166143293671931D-08, 0.5990909213118121D+01,
     &  0.5175289227690812D-07, 0.5501725495034997D+01,
     &  0.6419581525894787D-06, 0.5020274351287227D+01,
     &  0.6139401305210235D-05, 0.4545470293840825D+01,
     &  0.4586356164456607D-04, 0.4076379533242183D+01,
     &  0.2705406889774859D-03, 0.3612185969516418D+01,
     &  0.1271499268760875D-02, 0.3152165704938866D+01,
     &  0.4796631268346652D-02, 0.2695667488564915D+01,
     &  0.1461268197482865D-01, 0.2242097365824569D+01,
     &  0.3612646737149015D-01, 0.1790906348772901D+01,
     &  0.7276104720165333D-01, 0.1341580271678923D+01,
     &  0.1197332834080540D+00, 0.8936312256129486D+00,
     &  0.1613082628631473D+00, 0.4465901175404579D+00,
     &  0.1781337719310835D+00, 0.0000000000000000D+00/
      DATA ( W(I,50), T(I,50), I = 1,25) /
     &  0.1034607500576990D-36, 0.1298588445541555D+02,
     &  0.9443414659584510D-32, 0.1205301838092448D+02,
     &  0.6856280758924735D-28, 0.1127923332148262D+02,
     &  0.1206044550761014D-24, 0.1058738174919177D+02,
     &  0.7995094477915292D-22, 0.9948035709637500D+01,
     &  0.2522482807168144D-19, 0.9346039593575728D+01,
     &  0.4368171816201588D-17, 0.8772299579514598D+01,
     &  0.4566698246800344D-15, 0.8220815907982127D+01,
     &  0.3083828687005300D-13, 0.7687362406712500D+01,
     &  0.1414228936126661D-11, 0.7168814837853899D+01,
     &  0.4576636712310442D-10, 0.6662775399018720D+01,
     &  0.1077060789389039D-08, 0.6167347388659921D+01,
     &  0.1888225976835208D-07, 0.5680992291033284D+01,
     &  0.2514609880838772D-06, 0.5202434993399912D+01,
     &  0.2584937658949391D-05, 0.4730598550228594D+01,
     &  0.2078485175734569D-04, 0.4264557843038109D+01,
     &  0.1321726328668984D-03, 0.3803505741742012D+01,
     &  0.6708280619787080D-03, 0.3346727774732429D+01,
     &  0.2738160896935348D-02, 0.2893582727707738D+01,
     &  0.9045054154849623D-02, 0.2443487452654017D+01,
     &  0.2430481286424306D-01, 0.1995904709795124D+01,
     &  0.5334352453170102D-01, 0.1550333214338771D+01,
     &  0.9593054035810168D-01, 0.1106299289397183D+01,
     &  0.1416854132499443D+00, 0.6633496795082918D+00,
     &  0.1721258519924433D+00, 0.2210451816445435D+00/
*
*
****END HERMIT
*
      END
