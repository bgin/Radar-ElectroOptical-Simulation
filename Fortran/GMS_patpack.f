#if 0
*   SMPLTST is a simple test driver for PATSYM.
*
*    Sample output.
*
*       PATSYM TEST RESULTS
*
*     FTEST CALLS =  391, IFAIL =  0
*    N   ESTIMATED ERROR    INTEGRAL
*    1     0.00002386     0.13850726
*    2     0.00002056     0.06369358
*    3     0.00002738     0.05861598
*    4     0.00003107     0.05406837
*    5     0.00002793     0.05005375
*    6     0.00000987     0.04654493
*
      EXTERNAL FTEST
      INTEGER N, NF, NW, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL
      PARAMETER ( NDIM = 5, NF = NDIM+1, NW = 5000 )
      DOUBLE PRECISION A(NDIM), B(NDIM), WRKSTR(NW)
      DOUBLE PRECISION ABSEST(NF), FINEST(NF), ABSREQ, RELREQ
      DO N = 1, NDIM
         A(N) = 0
         B(N) = 1
      END DO
      MINCLS = 0
      MAXCLS = 10000
      ABSREQ = 0
      RELREQ = 1E-4
      CALL PATSYM( NDIM, A, B, NF, MINCLS, MAXCLS, FTEST, ABSREQ, 
     *             RELREQ, 0, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR )
      PRINT 9999, NEVAL, IFAIL
 9999 FORMAT (8X, 'PATSYM TEST RESULTS', //'     FTEST CALLS = ', I4,
     * ', IFAIL = ', I2, /'    N   ESTIMATED ERROR    INTEGRAL')
      DO N = 1, NF
         PRINT 9998, N, ABSEST(N), FINEST(N)
 9998    FORMAT (3X, I2, 2F15.8)
      END DO
      END
*
#endif
      SUBROUTINE FTEST( NDIM, Z, NFUN, F ) 
      INTEGER N, NDIM, NFUN
      DOUBLE PRECISION Z(NDIM), F(NFUN), SUM
      SUM = 0
      DO N = 1, NDIM
         SUM = SUM + N*Z(N)**2
      END DO
      F(1) = EXP(-SUM/2)
      DO N = 1, NDIM
         F(N+1) = Z(N)*F(1)
      END DO
      END
*
      SUBROUTINE PATSYM( NDIM, A, B, NUMFUN, MINPTS, MAXPTS, FUNSUB,
     &                   EPSABS, EPSREL, RESTAR, RESULT, ABSERR, NEVAL,
     &                   IFAIL, WORK )
****BEGIN PROLOGUE PATSYM
*
****AUTHOR
*            Alan Genz, 
*            Department of Mathematics 
*            Washington State University
*            Pullman, WA 99164-3113, USA
*            Email: alangenz@wsu.edu
****KEYWORDS automatic multidimensional integrator.
*    
****PURPOSE  The routine calculates an approximation to a given
*            vector of definite integrals
*
*      B(1)      B(NDIM) 
*     I    ...  I        (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
*      A(1)      A(NDIM)   1  2      NUMFUN
*
*       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN.
*              I   I  1  2      NDIM
*
*
****DESCRIPTION Computation of integrals over hyperrectangular region.
*
*   ON ENTRY
*
*     NDIM   Integer.
*            Number of variables. 1 < NDIM <=  20.
*     A      Real array of integration lower limits.
*     B      Real array of integration upper limits.
*     NUMFUN Integer.
*            Number of components of the integral.
*     MINPTS Integer.
*            Minimum number of function evaluations.
*     MAXPTS Integer.
*            Maximum number of function evaluations.
*     FUNSUB Externally declared subroutine for computing
*            all components of the integrand at the given
*            evaluation point.
*            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
*            Input parameters:
*              NDIM   Integer that defines the dimension of the
*                     integral.
*              X      Real array of dimension NDIM
*                     that defines the evaluation point.
*              NUMFUN Integer that defines the number of
*                     components of I.
*            Output parameter:
*              FUNVLS Real array of dimension NUMFUN
*                     that defines NUMFUN components of the integrand.
*
*     EPSABS Real.
*            Requested absolute accuracy.
*     EPSREL Real.
*            Requested relative accuracy.
*     RESTAR Integer.
*            If RESTAR = 0, this is the first attempt to compute
*            the integral.
*            If RESTAR = 1, then we restart a previous attempt.
*     WORK   Real array of working storage, must have dimension at
*            least 2*NUMFUN + (NUMFUN+1)*NUMSMS. NUMSMS is the number 
*            of NDIM partitions of the integers 0,...,RULE, after
*            RULE+1 approximations to the integrals have been computed.
*            The required value for NUMSMS will never exceed 10000.
*
*   ON RETURN
*
*     RESULT Real array of dimension NUMFUN.
*            Approximations to all components of the integral.
*     ABSERR Real array of dimension NUMFUN.
*            Estimates of absolute accuracies, based on differences 
*              from successive calls to PATRUL.
*     NEVAL  Integer.
*            Number of function evaluations used.
*     IFAIL  Integer.
*            IFAIL = 0 for normal exit, when 
*              ABSERR(K) <= MAX( EPSABS, ABS(RESULT(K))*EPSREL ) for
*              all K, 0 < K <= NUMFUN, with <= MAXPTS function values. 
*            IFAIL = 1 if MAXPTS was too small to obtain the required 
*              accuracy. In this case values of RESULT are returned
*              with estimated absolute accuracies ABSERR.
*     WORK   Real array of working storage, that contains information
*            that might be needed for subsequent calls of PATSYM. This
*            array should not be modified. WORK(NUMFUN+1) contains the
*            current value for NUMSMS.
*
****ROUTINES CALLED PATRUL
****END PROLOGUE PATSYM
*
*   Global variables.
*
      EXTERNAL FUNSUB
      INTEGER NDIM, NUMFUN, MINPTS, MAXPTS, RESTAR
      INTEGER NEVAL, IFAIL
      DOUBLE PRECISION A(*), B(*), EPSABS, EPSREL
      DOUBLE PRECISION RESULT(*), ABSERR(*), WORK(*)
*
*   Local variables.
*
      INTEGER I, MAXRUL, RULE, MNRULE, INTCLS, NUMSMS
      PARAMETER ( MAXRUL = 47 )
      SAVE RULE, MNRULE
      IFAIL = 1
      IF ( RESTAR .EQ. 0 ) THEN
         MNRULE = -1
         RULE = 0
      ELSE
         DO I = 1, NUMFUN
            RESULT(I) = WORK(I) 
         END DO
      END IF
      NEVAL = 0
      DO WHILE ( NEVAL .LE. MAXPTS .AND. RULE .LE. MAXRUL .AND. 
     &           ( NEVAL .LE. MINPTS .OR. IFAIL .GT. 0 ) ) 
         CALL PATRUL( NDIM, A, B, NUMFUN, FUNSUB, MNRULE, RULE, RESULT, 
     &        INTCLS, WORK(NUMFUN+1), WORK(2*NUMFUN+1), NUMSMS ) 
         NEVAL = NEVAL + INTCLS
         IFAIL = 0
         DO I = 1, NUMFUN
            IF ( RULE .GT. 0 ) THEN
               ABSERR(I) = ABS( RESULT(I) - WORK(I)  )
            ELSE
               ABSERR(I) = ABS( RESULT(I) )
            ENDIF
            WORK(I) = RESULT(I)
            IF ( ABSERR(I) .GT. MAX(EPSABS,EPSREL*ABS(RESULT(I))) )
     &           IFAIL = 1
         END DO
         RULE = RULE + 1
      END DO
      WORK(NUMFUN+1) = NUMSMS
*
****END PATSYM
*
      END
*
      SUBROUTINE PATRUL( S, A, B, N, F, MINORD, MAXORD, INTVAL, INTCLS, 
     &                   WORK, FULSMS, NUMSMS )
*
*  MULTIDIMENSIONAL FULLY SYMMETRIC RULE INTEGRATION SUBROUTINE
*
*   THIS SUBROUTINE COMPUTES A SEQUENCE OF FULLY SYMMETRIC RULE
*   APPROXIMATIONS TO A FULLY SYMMETRIC MULTIPLE INTEGRAL. WRITTEN BY 
*      Alan Genz
*      Department of Mathematics
*      Washington State University
*      Pullman, Washington 99164-3113  USA
*      Email: alangenz@wsu.edu
*    Reference:
*      Alan Genz, Fully Symmetric Interpolatory Rules for Multiple
*       Integrals, SIAM J. Numer. Anal. 23 (1986), pp. 1273-1283.
*
*
***************  PARAMETERS FOR PATRUL  ********************************
******INPUT PARAMETERS
*  S       Integer number of variables, must exceed 0 but not exceed 20.
*  A       Real array of integration lower limits.
*  B       Real array of integration upper limits.
*  N       Integer number of components for the integral
*  F       Externally declared user defined integrand subroutine.
*          It must have parameters (S,X,N,FUN), where X is a real array
*          with length S and FUN is a real array of length N.
*  MINORD  Integer minimum order parameter; on entry MINORD specifies
*          the current highest order approximation to the integral,
*          available in the array INTVLS.  For the first call of PATRUL
*          MINORD should be set to -1.  Otherwise a previous call is
*          assumed that computed INTVAL. On exit MINORD is set to MAXORD.
*  MAXORD  Integer maximum order parameter, must be greater than MINORD
*          and not exceed 47. 
*******OUTPUT PARAMETERS
*  INTVAL  Real array of length MAXORD; upon successful exit
*          INTVAL(1),..., INTVAL(N) are approximations to the components
*          of the integral.  These are all approximations of polynomial 
*          degree 2*MAXORD+1.
*  INTCLS  Integer number of F values needed for INTVAL
*  WORK    Real work array with length N. 
*  FULSMS  Real work array with dimension (N+1,*). On exit FULSMS(I,J) 
*          contains the fully symmetric basic rule sum indexed by the 
*          jth S-partition of the integers 0, ..., MAXORD, for the ith 
*          component of the integrand. FULSMS(N+1,J) contains number of 
*          points for the fully symmetric basic rule sum indexed by the 
*          jth S-partition of the integers 0, ..., MAXORD.
*  NUMSMS  Integer number of S-partition of the integers 0, ..., MAXORD.
************************************************************************
      EXTERNAL F
      INTEGER S, N, MAXDIM, MINORD, MAXORD, MAXRDM, NUMSMS
      PARAMETER ( MAXDIM = 20, MAXRDM = 47 ) 
      DOUBLE PRECISION A(*), B(*), X(MAXDIM), T(MAXDIM), FULWGT, WGTPAT, 
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
****  BEGIN LOOP FOR EACH D
*      FOR EACH D FIND ALL DISTINCT PARTITIONS M WITH MOD(M)<=D
*
      DO WHILE( D .LE. MIN ( MAXORD, MAXRDM ) )
         PRTCNT = 0
         CALL NXPART( PRTCNT, S, M, MODOFM )
         DO WHILE( MODOFM .LE. D )
*     
****  CALCULATE UPDATED WEIGHT FOR PARTITION M AND 
****     FULLY SYMMETRIC SUMS ( WHEN NECESSARY )
*
            FULWGT = WGTPAT( S, X, M, K, MODOFM, D )
            IF ( D .EQ. MODOFM ) THEN
               DO I = 1,N
                  FULSMS(I,PRTCNT) = 0
               END DO
               FULSMS(N+1,PRTCNT) = 0
            END IF
            IF ( FULSMS(N+1,PRTCNT) .EQ. 0 .AND. FULWGT .NE. 0 ) THEN
               CALL FLSMPT( S, A,B, M, N,F, FULSMS(1,PRTCNT), X,T,WORK ) 
     &                      
               INTCLS = INTCLS + FULSMS(N+1,PRTCNT)
            END IF
            DO I = 1,N 
               INTVAL(I) = INTVAL(I) + FULWGT*FULSMS(I,PRTCNT)
            END DO
            CALL NXPART( PRTCNT, S, M, MODOFM )
         END DO
*     
****  END LOOP FOR EACH D
*
         D = D + 1
      END DO
      MINORD = MAXORD
      NUMSMS = PRTCNT - 1
      END
*
      SUBROUTINE FLSMPT( S, A, B, M, N, F, FULSMS, X, T, FUNVAL )
*
****  To compute fully symmetric basic rule sums
*
      INTEGER S, M(*), SUMCLS, IX, LX, I,L, MI, ML, IL, N, MX
      DOUBLE PRECISION A(*), B(*), X(*), T(*), FULSMS(*), FUNVAL(*)
      PARAMETER ( MX = 31 )
      DOUBLE PRECISION G(0:MX), INTWGT
      SAVE G
*
*        Generators for 1 + 2 + 4 + 8 + 16 + 32 = 63 point degree
*        95 rule, with degree 1, 3, 11, 23, and 47 imbedded rules.     
*          Refererence:
*            T.N.L. Patterson, The Optimum Addition of Points to 
*             Quadrature Formulae, Math. Comp. 22 (1968), pp 847-856.
*
*
      DATA G( 0),G( 1)/0D0                    ,.77459666924148337704D0/
      DATA G( 2),G( 3)/.96049126870802028342D0,.43424374934680255800D0/
      DATA G( 4),G( 5)/.99383196321275502221D0,.22338668642896688163D0/
      DATA G( 6),G( 7)/.62110294673722640294D0,.88845923287225699889D0/
      DATA G( 8),G(11)/.99909812496766759750D0,.98153114955374010698D0/
      DATA G(14),G(13)/.92965485742974005664D0,.83672593816886873551D0/
      DATA G(15),G(10)/.70249620649152707861D0,.53131974364437562397D0/
      DATA G(12),G( 9)/.33113539325797683309D0,.11248894313318662575D0/
      DATA G(16),G(18)/.99987288812035761194D0,.99720625937222195908D0/
      DATA G(20),G(22)/.98868475754742947994D0,.97218287474858179658D0/
      DATA G(24),G(26)/.94634285837340290515D0,.91037115695700429250D0/
      DATA G(28),G(30)/.86390793819369047715D0,.80694053195021761186D0/
      DATA G(31),G(29)/.73975604435269475868D0,.66290966002478059546D0/
      DATA G(27),G(25)/.57719571005204581484D0,.48361802694584102756D0/
      DATA G(23),G(21)/.38335932419873034692D0,.27774982202182431507D0/
      DATA G(19),G(17)/.16823525155220746498D0,.56344313046592789972D-1/
      INTWGT = 1
      DO I = 1, S
        IF ( M(I) .NE. 0 ) INTWGT = INTWGT/2
        INTWGT = INTWGT*( B(I) - A(I) )/2
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
      DO I = 1,S
         T(I) = ( B(I) + A(I) + ( B(I) - A(I) )*X(I) )/2 
      END DO
      CALL F( S, T, N, FUNVAL )
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
      DOUBLE PRECISION FUNCTION WGTPAT( S, INTRPS, M,K, MODOFM, D )
*
****  Subroutine to update weight for partition m
*
      INTEGER S, M(S), K(S), I, L, D, MAXRDM, NZRMAX, MODOFM
      PARAMETER ( MAXRDM = 47 )
      DOUBLE PRECISION INTRPS(S), MOMPRD(0:MAXRDM,0:MAXRDM), MOMNKN
      PARAMETER ( NZRMAX = 31 )
      DOUBLE PRECISION A(0:NZRMAX), G(0:NZRMAX)
      SAVE A, G, MOMPRD
      DATA A / 2D0,  0.66666666666666667D0, 
     &         0D0,  0.45714285714285714D-01,
     &       2*0D0, -0.28065213250398436D-03, 0.13157729695421112D-03,
     &       4*0D0, -0.53949060310550432D-08, 0.16409565802196882D-07,
     &              -0.12211217614373411D-07, 0.52317265561235989D-08,
     &       8*0D0, -0.10141155834616524D-17, 0.14345827598358802D-16,
     &              -0.56865230056143054D-16, 0.11731814910797153D-15,
     &              -0.95580354100927967D-16, 0.64242918014064288D-16,
     &              -0.12072769909636026D-16, 0.19636450073868758D-17/
*
*        Generators for 1 + 2 + 4 + 8 + 16 + 32 = 63 point degree
*        95 rule, with degree 1, 3, 11, 23, and 47 imbedded rules.     
*          Refererence:
*            T.N.L. Patterson, The Optimum Addition of Points to 
*             Quadrature Formulae, Math. Comp. 22 (1968), pp 847-856.
*
      DATA G( 0),G( 1)/0D0                    ,.77459666924148337704D0/
      DATA G( 2),G( 3)/.96049126870802028342D0,.43424374934680255800D0/
      DATA G( 4),G( 5)/.99383196321275502221D0,.22338668642896688163D0/
      DATA G( 6),G( 7)/.62110294673722640294D0,.88845923287225699889D0/
      DATA G( 8),G(11)/.99909812496766759750D0,.98153114955374010698D0/
      DATA G(14),G(13)/.92965485742974005664D0,.83672593816886873551D0/
      DATA G(15),G(10)/.70249620649152707861D0,.53131974364437562397D0/
      DATA G(12),G( 9)/.33113539325797683309D0,.11248894313318662575D0/
      DATA G(16),G(18)/.99987288812035761194D0,.99720625937222195908D0/
      DATA G(20),G(22)/.98868475754742947994D0,.97218287474858179658D0/
      DATA G(24),G(26)/.94634285837340290515D0,.91037115695700429250D0/
      DATA G(28),G(30)/.86390793819369047715D0,.80694053195021761186D0/
      DATA G(31),G(29)/.73975604435269475868D0,.66290966002478059546D0/
      DATA G(27),G(25)/.57719571005204581484D0,.48361802694584102756D0/
      DATA G(23),G(21)/.38335932419873034692D0,.27774982202182431507D0/
      DATA G(19),G(17)/.16823525155220746498D0,.56344313046592789972D-1/
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
      WGTPAT = INTRPS(S)
      END
*
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
