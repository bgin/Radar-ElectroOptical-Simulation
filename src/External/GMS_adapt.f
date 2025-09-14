      SUBROUTINE ADAPT( NDIM, NUMFUN, A, B, MINPTS, MAXPTS, FUNSUB,
     &     EPSABS, EPSREL, KEY, NW, RESTAR, RESULT, ABSERR, NEVAL,
     &     IFAIL, WORK)
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ADAPT
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: ADAPT
#endif
#if 0
****BEGIN PROLOGUE ADAPT
****AUTHOR
*            Alan Genz, Department of Mathematics, Washington State
*            University, Pullman, WA 99164-3113, USA
*              Email: alangenz@wsu.edu
****KEYWORDS automatic multidimensional integrator,
*            n-dimensional hyper-rectangles,
*            general purpose, global adaptive
****PURPOSE  The routine calculates an approximation to a given
*            vector of definite integrals
*
*      B(1) B(2)     B(NDIM)
*     I    I    ... I       (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
*      A(1) A(2)     A(NDIM)  1  2      NUMFUN
*
*       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN.
*              I   I  1  2      NDIM
*
*            hopefully satisfying for each component of I the following
*            claim for accuracy:
*            ABS( I(K)-RESULT(K) ) .LE. MAX( EPSABS, EPSREL*ABS(I(K)) )
****DESCRIPTION Computation of integrals over hyper-rectangular
*            regions.
*            ADAPT is a driver for the integration routine
*            ADBASE, which repeatedly subdivides the region
*            of integration and estimates the integrals and the
*            errors over the subregions with greatest
*            estimated errors until the error request
*            is met or MAXPTS function evaluations have been used.
*
*   ON ENTRY
*
*     NDIM   Integer.
*            Number of variables. 0 < NDIM <=  20.
*     NUMFUN Integer.
*            Number of components of the integral.
*     A      Real array of dimension NDIM.
*            Lower limits of integration.
*     B      Real array of dimension NDIM.
*            Upper limits of integration.
*     MINPTS Integer.
*            Minimum number of function evaluations.
*     MAXPTS Integer.
*            Maximum number of function evaluations.
*            The number of function values for each subregion is NUM.
*            If NDIM = 1 Then NUM = 15
*            ElseIf KEY = 0 Then 
*                  if NDIM < 12 then NUM = 1 + 2*NDIM*(NDIM+3) + 2**NDIM
*                             else NUM = 1 + 2*NDIM*(NDIM+4)
*               Elseif KEY = 1 Then NUM = 1 + 2*NDIM*(NDIM+3) + 2**NDIM
*               Elseif KEY = 2 Then NUM = 1 + 4*NDIM + 6*NDIM*NDIM 
*                                 + 4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
*               Else NUM = 1 + 2*NDIM*(NDIM+4).
*            You must have MAXPTS >= NUM and MAXPTS >= MINPTS.
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
*     KEY    Integer.
*            Key to selected local integration rule.
*            KEY = 0 gives the user a default rule
*            KEY = 1 gives the user a degree 7 integration rule.
*                  This is the recommended general purpose rule.
*            KEY = 2 gives the user a degree 9 integration rule.
*                  This rule is recommended for oscillatory problems.
*            KEY = 3 gives the user a degree 5 integration rule.
*     NW     Integer.
*            Defines the length of the working array WORK.
*            Let MAXSUB denote the maximum allowed number of subregions
*            for the given values of MAXPTS, KEY and NDIM.
*            With MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1, you must have 
*             NW >= MAXSUB*( 2*NDIM + 2*NUMFUN + 2 ) + 7*NUMFUN + NDIM.
*     RESTAR Integer.
*            If RESTAR = 0, this is the first attempt to compute
*            the integral.
*            If RESTAR = 1, then we restart a previous attempt.
*            In this case the only parameters for ADAPT that may
*            be changed (with respect to the previous call of ADAPT)
*            are MINPTS, MAXPTS, EPSABS, EPSREL, KEY and RESTAR.
*
*   ON RETURN
*
*     RESULT Real array of dimension NUMFUN.
*            Approximations to all components of the integral.
*     ABSERR Real array of dimension NUMFUN.
*            Estimates of absolute accuracies.
*     NEVAL  Integer.
*            Number of function evaluations used by ADAPT.
*     IFAIL  Integer.
*            IFAIL = 0 for normal exit, when 
*              ABSERR(K) <= MAX( EPSABS, ABS(RESULT(K))*EPSREL ) for
*              all K, 0 < K <= NUMFUN, with <= MAXPTS function values. 
*            IFAIL = 1 if MAXPTS was too small to obtain the required 
*              accuracy. In this case values of RESULT are returned
*              with estimated absolute accuracies ABSERR.
*            IFAIL = 2 if KEY is less than 0 or KEY > 3.
*            IFAIL = 3 if NDIM is less than 2 or NDIM > 20.
*            IFAIL = 4 if NUMFUN is less than 1.
*            IFAIL = 5 if volume of region of integration is zero.
*            IFAIL = 6 if MAXPTS is less than NUM.
*            IFAIL = 7 if MAXPTS is less than MINPTS.
*            IFAIL = 8 if EPSABS < 0 and EPSREL < 0.
*            IFAIL = 9 if NW is too small.
*            IFAIL = 10 if RESTAR < 0 or RESTAR > 1.
*     WORK   Real array of dimension NW, used as working storage.
*            Let WRKSUB = ( NW - NDIM - 6*NUMFUN )/( 2*NDIM+2*NUMFUN+2 )
*            WORK(1),...,WORK(NUMFUN*MAXSUB) contain
*              the estimated components of the integrals over the
*              subregions.
*            WORK(NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*MAXSUB) contain
*              the estimated errors over the subregions.
*            WORK(2*NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB+NDIM*
*              MAXSUB) contain the centers of the subregions.
*            WORK(2*NUMFUN*WRKSUB+NDIM*WRKSUB+1),...,WORK((2*NUMFUN+
*              NDIM)*WRKSUB+NDIM*MAXSUB) contain subregion half widths.
*            WORK(2*NUMFUN*WRKSUB+2*NDIM*WRKSUB+1),...,WORK(2*NUMFUN*
*              WRKSUB+2*NDIM*WRKSUB+MAXSUB) contain the greatest errors
*              in each subregion.
*            WORK((2*NUMFUN+2*NDIM+1)*WRKSUB+1),...,WORK((2*NUMFUN+
*              2*NDIM+1)*WRKSUB+MAXSUB) contain the heap pointers
*              for the subregions.
*
****ROUTINES CALLED BASCHC, ADBASE
*  REFERENCES
*   P. van Dooren and L. de Ridder Algorithm 6: 
*    An adaptive algorithm for numerical integration over 
*    an n-dimensional cube, J.Comput.Appl. Math. 2(1976)207-217.
*
*   A.C. Genz and A.A. Malik, Algorithm 019. Remarks on algorithm 006:
*    An adaptive algorithm for numerical integration over an
*    N-dimensional rectangular region, J.Comp.Appl.Math. 6(1980)295-302.
*
*   J. Berntsen, T.O. Espelid, and A. Genz, An Adaptive Algorithm 
*    for the Approximate Calculation of Multiple Integrals,
*    ACM Trans.Math.Softw., 17(1991)437-451.
*
*
* Copyright (C) 2013, Alan Genz,  All rights reserved.               
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided the following conditions are met:
*   1. Redistributions of source code must retain the above copyright
*      notice, this list of conditions and the following disclaimer.
*   2. Redistributions in binary form must reproduce the above copyright
*      notice, this list of conditions and the following disclaimer in 
*      the documentation and/or other materials provided with the 
*      distribution.
*   3. The contributor name(s) may not be used to endorse or promote 
*      products derived from this software without specific prior 
*      written permission.
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
* TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
****END PROLOGUE ADAPT
*
*   Global variables.
*
#endif
      EXTERNAL FUNSUB
      integer(kind=4) ::  NDIM, NUMFUN, MINPTS, MAXPTS, KEY, NW, RESTAR
      integer(kind=4) ::  NEVAL, IFAIL
      real(kind=8) ::  A(NDIM), B(NDIM), EPSABS, EPSREL
      real(kind=8) ::  RESULT(NUMFUN), ABSERR(NUMFUN), WORK(NW)
!*
!*   Local variables.
#if 0
*
*   MAXDIM Integer.
*          The maximum allowed value of NDIM.
*   MAXSUB Integer.
*          The maximum allowed number of subdivisions
*          for the given values of KEY, NDIM and MAXPTS.
*   MINSUB Integer.
*          The minimum allowed number of subregions for the given
*          values of MINPTS, KEY and NDIM.
*   WRKSUB Integer.
*          The maximum allowed number of subregions as a function 
*          of NW, NUMFUN, NDIM and NPROC. This determines the length
*          of the main work arrays.
*   NUM    Integer. The number of integrand evaluations needed
*          over each subregion.
*
#endif
      integer(kind=4) ::  MAXDIM, MAXSUB, MINSUB, NUM, NSUB, NEWPTS, NEWCLS, TOTCLS       
      PARAMETER ( MAXDIM = 20 )
      integer(kind=4) ::  WRKSUB, I, J, I1, I2, I3, I4, I5, I6, I7, RS 
      SAVE NSUB
!*
!****FIRST EXECUTABLE STATEMENT ADBAYS
!*
!*   Compute NUM, WTLENG, MAXSUB and MINSUB,
!*   and check the input parameters.
!*
!*
!*   On restart runs the number of subregions from the
!*   previous call is assigned to NSUB.
!*
      IF ( RESTAR .EQ. 0 ) THEN
         NSUB = 1
         DO I = 1, NUMFUN
            WORK( NW - NUMFUN + I ) = 0
         END DO
         RS = 0
      ELSE
         RS = 1
      END IF
      NEWPTS = MAXPTS
      TOTCLS = 0
!*
!*   Split up the work space.
!*
      WRKSUB = ( NW - NDIM - 6*NUMFUN )/( 2*NDIM + 2*NUMFUN + 1 )
      I1 = 1
      I2 = I1 + WRKSUB*NUMFUN
      I3 = I2 + WRKSUB*NUMFUN
      I4 = I3 + WRKSUB*NDIM
      I5 = I4 + WRKSUB*NDIM
      I6 = I5 + WRKSUB
      I7 = I6 + NDIM
      CALL BASCHC( MAXDIM, NDIM, NUMFUN, A,B, MINPTS, NEWPTS, EPSABS,
     &     EPSREL, KEY, NW, RS, NUM, NSUB, MAXSUB,MINSUB, IFAIL )
      IF ( IFAIL .EQ. 0 ) THEN
         IF ( NDIM .EQ. 1 ) THEN
            CALL ADONEV( NUMFUN, A(1), B(1), MINSUB, MAXSUB, FUNSUB,
     &           EPSABS, EPSREL, RS, RESULT, ABSERR, 
     &           NEWCLS, NSUB, IFAIL, WORK(I1), WORK(I2), WORK(I3), 
     &           WORK(I4), WORK(I5), WORK(I6) )
         ELSE
            CALL ADBASE( NDIM, NUMFUN, A, B, MINSUB, MAXSUB, FUNSUB,
     &           EPSABS, EPSREL, KEY, RS, NUM, RESULT, ABSERR, 
     &           NEWCLS, NSUB, IFAIL, WORK(I1), WORK(I2), WORK(I3), 
     &           WORK(I4), WORK(I5), WORK(I6), WORK(I7) )
         END IF
         TOTCLS = TOTCLS + NEWCLS
         DO I = 1,NUMFUN
            WORK( NW - NUMFUN + I ) = RESULT(I)
         END DO
      END IF
      RS = 1
      NEVAL = TOTCLS
!*
!****END ADAPT
1*
      END

      SUBROUTINE BASCHC( MAXDIM, NDIM, NUMFUN, A,B, MINPTS,MAXPTS,
     & EPSABS,EPSREL, KEY, NW, RESTAR, NUM, NSUB, MAXSUB,MINSUB, IFAIL )
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: BASCHC
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: BASCHC
#endif
#if 0
****BEGIN PROLOGUE BASCHC
****PURPOSE  BASCHC checks the validity of the
*            input parameters to ADAPT.
****DESCRIPTION
*            BASCHC computes NUM, MAXSUB, MINSUB and IFAIL as
*            functions of the input parameters to ADAPT,
*            and checks the validity of the input parameters to ADAPT.
*
*   ON ENTRY
*
*     MAXDIM Integer.
*            The maximum allowed number of dimensions.
*     NDIM   Integer.
*            Number of variables. 1 < NDIM <= MAXDIM.
*     NUMFUN Integer.
*            Number of components of the integral.
*     A      Real array of dimension NDIM.
*            Lower limits of integration.
*     B      Real array of dimension NDIM.
*            Upper limits of integration.
*     MINPTS Integer.
*            Minimum number of function evaluations.
*     MAXPTS Integer.
*            Maximum number of function evaluations.
*            The number of function values used in each subregion is NUM.
*            If NDIM = 1 Then NUM = 15
*             ElseIf KEY = 0 Then 
*                  if NDIM < 12 then NUM = 1 + 2*NDIM*(NDIM+3) + 2**NDIM
*                             else NUM = 1 + 2*NDIM*(NDIM+4)
*               Elseif KEY = 1 Then NUM = 1 + 2*NDIM*(NDIM+3) + 2**NDIM
*               Elseif KEY = 2 Then NUM = 1 + 4*NDIM + 6*NDIM*NDIM 
*                                 + 4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
*               Else NUM = 1 + 2*NDIM*(NDIM+4).
*     EPSABS Real.
*            Requested absolute accuracy.
*     EPSREL Real.
*            Requested relative accuracy.
*     KEY    Integer.
*            Key to selected local integration rule.
*            KEY = 0 gives the user a default rule
*            KEY = 1 gives the user a degree 7 integration rule.
*                  This is the recommended general purpose rule.
*            KEY = 2 gives the user a degree 9 integration rule.
*                  This rule is recommended for oscillatory problems.
*            KEY = 3 gives the user a degree 5 integration rule.
*     NW     Integer.
*            Defines the length of the working array WORK.
*            Let MAXSUB denote the maximum allowed number of subregions
*            for the given values of MAXPTS, KEY and NDIM.
*            MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1
*            NW should be greater or equal to
*            MAXSUB*( 2*NDIM + 2*NUMFUN + 1 ) + 6*NUMFUN + NDIM
*     NSUB   Integer.
*            If RESTAR = 1, then NSUB must specify the number
*              of subregions stored in the previous call to ADBASE.
*     RESTAR Integer.
*            If RESTAR = 0, this is the first attempt to compute
*            the integral.
*            If RESTAR = 1, then we restart a previous attempt.
*
*   ON RETURN
*
*     NUM    Integer.
*            The number of function evaluations over each subregion.
*     MAXSUB Integer.
*            The maximum allowed number of subregions for the
*            given values of MAXPTS, KEY and NDIM.
*     MINSUB Integer.
*            The minimum allowed number of subregions for the given
*            values of MINPTS, KEY and NDIM.
*     IFAIL  Integer.
*            IFAIL = 0 for normal exit.
*            IFAIL = 2 if KEY < 0 or KEY > 3.
*            IFAIL = 3 if NDIM < 1 or NDIM > MAXDIM.
*            IFAIL = 4 if NUMFUN less than 1.
*            IFAIL = 5 if volume of region of integration is zero.
*            IFAIL = 6 if MAXPTS is less than NUM.
*            IFAIL = 7 if MAXPTS is less than MINPTS.
*            IFAIL = 8 if EPSABS < 0 and EPSREL < 0.
*            IFAIL = 9 if NW is too small.
*            IFAIL = 10 if illegal RESTAR.
*
****  ROUTINES CALLED NONE
****  END PROLOGUE BASCHC
*
*     Global variables.
*
#endif
      integer(kind=4) ::  NDIM, NUMFUN, MINPTS, MAXPTS, KEY, NW, MAXSUB, RESTAR
      integer(kind=4) ::  NSUB, NUM, IFAIL, MAXDIM, MINSUB
      real(kind=8) ::  A(NDIM), B(NDIM), EPSABS, EPSREL
!*
!*     Local variables.
!*
      integer(kind=4) ::  LIMIT,J
!*
!1****FIRST EXECUTABLE STATEMENT BASCHC
!*
      IFAIL = 0
!*
!*     Check on legal KEY.
!*
      IF ( KEY .LT. 0 .OR. KEY .GT. 3 ) IFAIL = 2
!*
!*     Check on legal NDIM.
!*
      IF ( NDIM .LT. 1 .OR. NDIM .GT. MAXDIM ) IFAIL = 3
!*
!*     Compute NUM as a function of KEY and NDIM.
!*
      IF ( KEY .EQ. 0 ) THEN
         IF ( NDIM .LT. 12 ) THEN
            NUM = 1 + 2*NDIM*(NDIM+3) + 2**NDIM
         ELSE
            NUM = 1 + 2*NDIM*(NDIM+4)
         END IF
      ELSE IF ( KEY .EQ. 1 ) THEN
          NUM = 1 + 2*NDIM*(NDIM+3) + 2**NDIM
      ELSE IF ( KEY .EQ. 2 ) THEN
          NUM = 1 + 4*NDIM + 6*NDIM*NDIM 
     &            + 4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
      ELSE 
          NUM = 1 + 2*NDIM*(NDIM+4)
      END IF
      IF ( NDIM .EQ. 1 ) NUM = 15
!*
!*     Check on positive NUMFUN.
!*
      IF ( NUMFUN .LT. 1 ) IFAIL = 4
!*
!*     Check for legal upper and lower limits of integration.
!!*
      DO J = 1, NDIM
         IF ( A(J) - B(J) .EQ. 0 ) IFAIL = 5
      END DO
!*
!*     Check for MAXPTS < NUM.
!*
      IF ( MAXPTS .LT. NUM ) IFAIL = 6
!*
!*     Check for MAXPTS >= MINPTS.
!*
      IF ( MAXPTS .LT. MINPTS) IFAIL = 7
!*
!*     Compute MAXSUB.
!*
      MAXSUB = MAXPTS/(2*NUM) + NSUB
!*     
!*     Compute MINSUB.
!*
      MINSUB = MAX( 1, MINPTS/(2*NUM) + 1 ) 
!*
!*     Check accuracy requests.
!*
      IF ( EPSABS .LT. 0 .AND. EPSREL .LT. 0 ) IFAIL = 8
!*
!*     Check workspace size.
!*
      LIMIT = MAXSUB* ( 2*NDIM+2*NUMFUN+1 ) + 6*NUMFUN + NDIM
      IF ( NW .LT. LIMIT ) IFAIL = 9
!*     
!*     Check RESTAR.
!*
      IF ( RESTAR .NE. 0 .AND. RESTAR .NE. 1 ) IFAIL = 10
!*
!****END BASCHC
!*
      END
      SUBROUTINE ADBASE( NDIM, NUMFUN, A, B, MINSUB, MAXSUB, FUNSUB,
     & EPSABS, EPSREL, KEY, RESTAR, NUM, RESULT, ABSERR, NEVAL, NSUB,
     & IFAIL, VALUES, ERRORS, CENTRS, HWIDTS, GREATE, X, WORK )
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ADBASE
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: ADBASE
#endif
#if 0
****BEGIN PROLOGUE ADBASE
****KEYWORDS automatic multidimensional integrator,
*            n-dimensional hyper-rectangles,
*            general purpose, global adaptive
****PURPOSE  The routine calculates an approximation to a given
*            vector of definite integrals, I, over a hyper-rectangular
*            region hopefully satisfying for each component of I the
*            following claim for accuracy:
*            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
****DESCRIPTION Computation of integrals over hyper-rectangular
*            regions.
*            ADBASE repeatedly subdivides the region
*            of integration and estimates the integrals and the
*            errors over the subregions with  greatest
*            estimated errors until the error request
*            is met or MAXSUB subregions are stored.
*            The regions are devided in two equally sized parts along
*            the direction with greatest absolute fourth divided
*            difference.
*
*   ON ENTRY
*
*     NDIM   Integer.
*            Number of variables. 1 < NDIM <= MAXDIM.
*     NUMFUN Integer.
*            Number of components of the integral.
*     A      Real array of dimension NDIM.
*            Lower limits of integration.
*     B      Real array of dimension NDIM.
*            Upper limits of integration.
*     MINSUB Integer.
*            The computations proceed until there are at least
*            MINSUB subregions in the data structure.
*     MAXSUB Integer.
*            The computations proceed until there are at most
*            MAXSUB subregions in the data structure.
*
*     FUNSUB Externally declared subroutine for computing
*            all components of the integrand in the given
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
*     KEY    Integer.
*            Key to selected local integration rule.
*     RESTAR Integer.
*            If RESTAR = 0, this is the first attempt to compute
*            the integral.
*            If RESTAR = 1, then we restart a previous attempt.
*              (In this case the output parameters must not be changed 
*               since the last exit.)
*     NUM    Integer.
*            The number of function evaluations over each subregion.
*     NSUB   Integer.
*            If RESTAR = 1, then NSUB must specify the number
*              of subregions stored in the previous call to ADBASE.
*
*   ON RETURN
*
*     RESULT Real array of dimension NUMFUN.
*            Approximations to all components of the integral.
*     ABSERR Real array of dimension NUMFUN.
*            Estimates of absolute accuracies.
*     NEVAL  Integer.
*            Number of function evaluations used.
*     NSUB   Integer.
*            Number of stored subregions.
*     IFAIL  Integer.
*            IFAIL = 0 for normal exit, when 
*              ABSERR(K) <= MAX( EPSABS, ABS(RESULT(K))*EPSREL ) for all
*              K, 1 <= K <= NUMFUN, with <= MAXSUB subregions processed. 
*            IFAIL = 1 if MAXSUB was too small to obtain the required 
*              accuracy. In this case values of RESULT with estimated
*              absolute accuracies ABSERR are returned.
*     VALUES Real array of dimension (NUMFUN,*).
*            Used to store estimated values of the integrals
*            over the subregions.
*     ERRORS Real array of dimension (NUMFUN,*).
*            Used to store the corresponding estimated errors.
*     CENTRS Real array of dimension (NDIM,*).
*            Used to store the centers of the stored subregions.
*     HWIDTS Real array of dimension (NDIM,*).
*            Used to store the half widths of the stored subregions.
*     GREATE Real array of dimension (*).
*            Used to store the greatest estimated errors in
*            all subregions.
*     PONTRS Real array of dimension (*).
*            PONTRS is used to store heap pointers.
*     WORK   Real work array of length at least 5*NUMFUN
*            Used  in BASRUL and SDFFER.
*     X      Real array of length NDIM.
*            Work array used in BASRUL.
*
****REFERENCES
*
*   P. van Dooren and L. de Ridder, Algorithm 6, An adaptive algorithm
*   for numerical integration over an n-dimensional cube, J.Comput.Appl.
*   Math. 2(1976)207-217.
*
*   A.C. Genz and A.A. Malik, Algorithm 019. Remarks on algorithm 006:
*   An adaptive algorithm for numerical integration over an
*   N-dimensional rectangular region,J.Comput.Appl.Math. 6(1980)295-302.
*
****  ROUTINES CALLED TRESTR, BASRUL, SDFFER
****  END PROLOGUE ADBASE
*
*     Global variables.
*
#endif
      EXTERNAL FUNSUB
      integer(kind=4) ::  NDIM, NUMFUN, MINSUB, MAXSUB, KEY, RESTAR
      integer(kind=4) ::  NUM, NEVAL, NSUB, IFAIL
      real(kind=8) ::  A(NDIM),B(NDIM), EPSABS, EPSREL
      real(kind=8) ::  RESULT(NUMFUN), ABSERR(NUMFUN)
      real(kind=8) ::  VALUES(NUMFUN,*), ERRORS(NUMFUN,*)
      real(kind=8) ::  CENTRS(NDIM,*), HWIDTS(NDIM,*)
      real(kind=8) ::  GREATE(*), WORK(*), X(*)
!*
!*     Local variables.
!*
!*   INTSGN is used to get correct sign on the integral.
!*   SBRGNS is the number of stored subregions.
!*   POINTR Pointer to the position in the datastructure where
!*          the new subregions are to be stored.
!*
      integer(kind=4) ::  I, J, INTSGN, SBRGNS, POINTR, DIRECT, WTLENG 
      integer(kind=4) ::  MAXDIM, MAXWTS, NUMNUL
      PARAMETER ( MAXWTS = 9, NUMNUL = 5, MAXDIM = 20)
      real(kind=8) ::  G(MAXWTS*MAXDIM), W(MAXWTS*NUMNUL) 
!*
!****  FIRST PROCESSING STATEMENT for ADBASE
!*
!*     Call BSINIT to compute the weights and abscissas of
!!*     the function evaluation points.
!*
      CALL BSINIT( NDIM, KEY, WTLENG, G, W )
!*
!*     Get the correct sign on the integral.
!*
      INTSGN = 1
      DO J = 1,NDIM
         IF ( B(J) .LT. A(J) ) INTSGN = -INTSGN
      END DO
      NEVAL = 0
      SBRGNS = NSUB
      IF ( RESTAR .EQ. 0 ) THEN
!*     
!*     Initialize the SBRGNS, CENTRS and HWIDTS.
!*
         DO J = 1,NDIM
            CENTRS(J,1) = ( A(J) + B(J) )/2
            HWIDTS(J,1) = ABS( B(J) - A(J) )/2
         END DO
!*     
!*     Apply BASRUL over the whole region.
!*     
         CALL BASRUL( NDIM, CENTRS, HWIDTS, WTLENG, G, W, 
     &        NUMFUN, FUNSUB, X, WORK, VALUES, ERRORS, GREATE(1) )
         NEVAL = NEVAL + NUM
!*     
!****  End initialisation.
!*     
      END IF
!*     
!*     Check for termination.
!*     
 10   IFAIL = 0
      DO J = 1,NUMFUN
         RESULT(J) = 0
         ABSERR(J) = 0
         DO I = 1, SBRGNS
            RESULT(J) = RESULT(J) + VALUES(J,I)*INTSGN
            ABSERR(J) = ABSERR(J) + ERRORS(J,I)
         END DO
         IF ( ABSERR(J) .GT. MAX( EPSABS, EPSREL*ABS(RESULT(J)) ) ) 
     &           IFAIL = 1
      END DO
!*
!****  Begin loop while the error is too great,
!*     and SBRGNS + 1 is less than MAXSUB.
!*     
      IF ( ( IFAIL .NE. 0 .AND. SBRGNS+1 .LE. MAXSUB )
     &                     .OR.  SBRGNS  .LT. MINSUB )   THEN
!!*     
!*     If we are allowed to divide further,
!*     prepare to apply basic rule over each half of the
!*     subregion with greatest error.
!*     
         POINTR = 1
         DO I = 2, SBRGNS
            IF ( GREATE(I) .GE. GREATE(POINTR) ) POINTR = I
         END DO
         CALL SDFFER( NDIM, CENTRS(1,POINTR), HWIDTS(1,POINTR),
     &        NUMFUN, FUNSUB, X, WORK, DIRECT )
!*     
!*     Divide the subregion in two halves. 
!*     
         HWIDTS(DIRECT,POINTR) = HWIDTS(DIRECT,POINTR)/2
         SBRGNS = SBRGNS + 1
         DO J = 1,NDIM
            CENTRS(J,SBRGNS) = CENTRS(J,POINTR)
            HWIDTS(J,SBRGNS) = HWIDTS(J,POINTR)
         END DO
!*
!*     Compute integral and error over first half and store results.
!*
         CENTRS(DIRECT,POINTR) = CENTRS(DIRECT,POINTR) 
     &                         - HWIDTS(DIRECT,POINTR)
         CALL BASRUL( NDIM, CENTRS(1,POINTR), HWIDTS(1,POINTR),
     &        WTLENG, G, W, NUMFUN, FUNSUB, X, WORK, 
     &        VALUES(1,POINTR), ERRORS(1,POINTR), GREATE(POINTR) )
!*
!*     Compute integral and error over second half and store results.
!*
         CENTRS(DIRECT,SBRGNS) = CENTRS(DIRECT,SBRGNS) 
     &                         + HWIDTS(DIRECT,SBRGNS)
         CALL BASRUL( NDIM, CENTRS(1,SBRGNS), HWIDTS(1,SBRGNS),
     &        WTLENG, G, W, NUMFUN, FUNSUB, X, WORK, 
     &        VALUES(1,SBRGNS), ERRORS(1,SBRGNS), GREATE(SBRGNS) )
         NEVAL = NEVAL + 2*NUM
         GO TO 10
      END IF
      NSUB = SBRGNS
!*
!****END ADBASE
!*
      END
      SUBROUTINE BSINIT( NDIM, KEY, WTLENG, G, W )
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: BSINIT
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: BSINIT
#endif
#if 0
****BEGIN PROLOGUE BSINIT
****PURPOSE BSINIT computes abscissas and weights of the integration
*            rule and the null rules to be used in error estimation.
*            These are computed as functions of NDIM and KEY.
****DESCRIPTION BSINIT will for given values of NDIM and KEY compute or
*            select the correct values of the abscissas and 
*            corresponding weights for different integration rules and 
*            null rules and assign them to G and W.
*
*   ON ENTRY
*
*     NDIM   Integer.
*            Number of variables.
*     KEY    Integer.
*            Key to selected local integration rule.
*
*   ON RETURN
*
*     WTLENG Integer.
*            The number of weights in each of the rules.
*     W      Real array of dimension (5,WTLENG).
*            The weights for the basic and null rules.
*            W(1,1), ...,W(1,WTLENG) are weights for the basic rule.
*            W(I,1), ...,W(I,WTLENG), for I > 1 are null rule weights.
*     G      Real array of dimension (NDIM,WTLENG).
*            The fully symmetric sum generators for the rules.
*            G(1,J),...,G(NDIM,J) are the generators for the points
*            associated with the the Jth weights.
*
****ROUTINES CALLED  BSRL05,BSRL07,BSRL09
****END PROLOGUE BSINIT
*
*   Variables.
*
#endif
      integer(kind=4) ::  NDIM, KEY, WTLENG
      real(kind=8) ::  G(*), W(*)
!*
!****FIRST EXECUTABLE STATEMENT BSINIT
!*
!*   Compute W and G.
!*
      IF ( KEY .EQ. 0 ) THEN
         IF ( NDIM .LT. 12 ) THEN
            CALL BSRL07( NDIM, WTLENG, W, G )
         ELSE
            CALL BSRL05( NDIM, WTLENG, W, G )
         ENDIF
      ELSE IF ( KEY .EQ. 1 ) THEN
         CALL BSRL07( NDIM, WTLENG, W, G )
      ELSE IF ( KEY .EQ. 2 ) THEN
         CALL BSRL09( NDIM, WTLENG, W, G )
      ELSE 
         CALL BSRL05( NDIM, WTLENG, W, G )
      END IF
!*
!****END BSINIT
!*
      END

      SUBROUTINE BSRL07( NDIM, WTLENG, W, G )
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: BSRL07
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: BSRL07
#endif
#if 0
****BEGIN PROLOGUE BSRL07
****KEYWORDS basic integration rule, degree 7
****PURPOSE  To initialize a degree 7 basic rule, and null rules.
****AUTHOR   
*            Alan Genz, Department of Mathematics, Washington State
*            University, Pullman, WA 99164-3113, USA
*              Email: genz@gauss.math.wsu.edu
****LAST MODIFICATION 88-05-31
****DESCRIPTION  BSRL07 initializes a degree 7 integration rule,
*            two degree 5 rules, one degree 3 rule and one
*            degree 1 rule for the hypercube [-1,1]**NDIM.
*
*   ON ENTRY
*
*   NDIM   Integer.
*          Number of variables.
*
*   ON EXIT
*
*   WTLENG Integer.
*          The number of weights in each of the rules.
*
*   W      Real array of dimension (5,WTLENG).
*          The weights for the basic and null rules.
*          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
*   G      Real array of dimension (NDIM, WTLENG).
*          The fully symmetric sum generators for the rules.
*          G(1, J), ..., G(NDIM, J) are the are the generators for the
*          points associated with the Jth weights.
*
****REFERENCES A. Genz and A. Malik,
*             "An Imbedded Family of Fully Symmetric Numerical
*              Integration Rules",
*              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
****ROUTINES CALLED NONE
****END PROLOGUE BSRL07
*
*   Global variables
*
#endif
      integer(kind=4) ::  NDIM,WTLENG
      real(kind=8) ::  W(5,*),G(NDIM,*)
!*
!*   Local Variables
!*
      real(kind=8) ::  RATIO,LAM0,LAM1,LAM2,LAMP,TWONDM
      integer(kind=4) ::  RULPTS(6)
      integer(kind=4) ::  I,J
!*
!****FIRST EXECUTABLE STATEMENT BSRL07
!*
!*
!*     Initialize generators, weights and RULPTS
!*
      WTLENG = 6
      DO J = 1,WTLENG
         DO I = 1,NDIM
            G(I,J) = 0
         END DO
         DO I = 1,5
            W(I,J) = 0
         END DO
         RULPTS(J) = 2*NDIM
      END DO
      TWONDM = 2**NDIM
      RULPTS(WTLENG) = TWONDM
      RULPTS(WTLENG-1) = 2*NDIM* (NDIM-1)
      RULPTS(1) = 1
!*
!*     Compute squared generator parameters
!*
      LAM0 = 0.4707_8
      LAMP = 0.5625_8
      LAM1 = 4/ (15-5/LAM0)
      RATIO = (1-LAM1/LAM0)/27
      LAM2 = (5-7*LAM1-35*RATIO)/ (7-35*LAM1/3-35*RATIO/LAM0)
!*
!*     Compute degree 7 rule weights
!*
      W(1,6) = 1/ (3*LAM0)**3/TWONDM
      W(1,5) = (1-5*LAM0/3)/ (60* (LAM1-LAM0)*LAM1**2)
      W(1,3) = (1-5*LAM2/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM2))/
     &         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(1,5)
      W(1,2) = (1-5*LAM1/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM1))/
     &         (10*LAM2* (LAM2-LAM1))
!*
!*     Compute weights for 2 degree 5, 1 degree 3 and 1 degree 1 rules
!*
      W(2,6) = 1/ (36*LAM0**3)/TWONDM
      W(2,5) = (1-9*TWONDM*W(2,6)*LAM0**2)/ (36*LAM1**2)
      W(2,3) = (1-5*LAM2/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM2))/
     &         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(2,5)
      W(2,2) = (1-5*LAM1/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM1))/
     &         (10*LAM2* (LAM2-LAM1))
      W(3,6) = 5/ (108*LAM0**3)/TWONDM
      W(3,5) = (1-9*TWONDM*W(3,6)*LAM0**2)/ (36*LAM1**2)
      W(3,3) = (1-5*LAMP/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAMP))/
     &         (10*LAM1* (LAM1-LAMP)) - 2* (NDIM-1)*W(3,5)
      W(3,4) = (1-5*LAM1/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAM1))/
     &         (10*LAMP* (LAMP-LAM1))
      W(4,6) = 1/ (54*LAM0**3)/TWONDM
      W(4,5) = (1-18*TWONDM*W(4,6)*LAM0**2)/ (72*LAM1**2)
      W(4,3) = (1-10*LAM2/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM2))/
     &         (20*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(4,5)
      W(4,2) = (1-10*LAM1/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM1))/
     &         (20*LAM2* (LAM2-LAM1))
!*
!*     Set generator values
!*
      LAM0 = SQRT(LAM0)
      LAM1 = SQRT(LAM1)
      LAM2 = SQRT(LAM2)
      LAMP = SQRT(LAMP)
      DO I = 1,NDIM
         G(I,WTLENG) = LAM0
      END DO
      G(1,WTLENG-1) = LAM1
      G(2,WTLENG-1) = LAM1
      G(1,WTLENG-4) = LAM2
      G(1,WTLENG-3) = LAM1
      G(1,WTLENG-2) = LAMP
!*
!*     Compute final weight values.
!*
      DO J = 1, 5
         W(J,1) = TWONDM
         DO I = 2, WTLENG
            W(J,I) = TWONDM*W(J,I)
            W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
         END DO
      END DO
      CALL RULNRM ( WTLENG, 5, RULPTS, W )
!*
!****END BSRL07
!*
      END

      SUBROUTINE BSRL05( NDIM, WTLENG, W, G )
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: BSRL05
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: BSRL05
#endif
#if 0
****BEGIN PROLOGUE BSRL05
****KEYWORDS basic integration rule, degree 5.
****PURPOSE  To initialize a degree 5 basic rule.
****AUTHOR   
*            Alan Genz, Department of Mathematics, Washington State
*            University, Pullman, WA 99163-3113, USA
*              Email: genz@gauss.math.wsu.edu
****LAST MODIFICATION 94-01-15
****DESCRIPTION  BSRL05 initializes a degree 5 integration rule,
*            two degree 3 rules and two degree 1 rules
*            for the hypercube [-1,1]**NDIM.
*
*   ON ENTRY
*
*   NDIM   Integer.
*          Number of variables.
*
*   ON EXIT
*
*   WTLENG Integer.
*          The number of weights in each of the rules.
*   W      Real array of dimension (5,WTLENG).
*          The weights for the basic and null rules.
*          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
*          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
*   G      Real array of dimension (NDIM, WTLENG).
*          The fully symmetric sum generators for the rules.
*          G(1, J), ..., G(NDIM, J) are the are the generators for the
*          points associated with the Jth weights.
*
****ROUTINES CALLED NONE
****END PROLOGUE BSRL05
*
*   Global variables
*
#endif
      integer(kind=4) ::  NDIM,WTLENG
      real(kind=8) ::  W(5,*), G(NDIM,*)
!*
!*   Local Variables
!*
      real(kind=8) ::   ONE, THREE, FIVE 
      PARAMETER ( ONE = 1, THREE = 3, FIVE = 5 )
      real(kind=8) ::  TWONDM, LAM1, LAM2, LAM3, LAM4
      integer(kind=4) ::  RULPTS(6), I, J
!*
!****FIRST EXECUTABLE STATEMENT BSRL05
!*
!*
!*     Initialize generators, weights and RULPTS
!*
      WTLENG = 6
      DO J = 1,WTLENG
         DO I = 1,NDIM
            G(I,J) = 0
         END DO
         DO I = 1,5
            W(I,J) = 0
         END DO
         RULPTS(J) = 2*NDIM
      END DO
      TWONDM = 2**NDIM
      RULPTS(WTLENG  ) = 2*NDIM*(NDIM-1)
      RULPTS(1) = 1
!*
!*     Compute squared generator parameters
!*
      LAM1 = 3/FIVE
      LAM2 = 1/THREE
      LAM3 = THREE/4
      LAM4 = 4/FIVE
!*
!*!     Compute degree 5 rule weights
!*
      W(1,WTLENG) = 1/(6*LAM1)**2
      W(1,WTLENG-1) = 1/(6*LAM1) - 2*(NDIM-1)*W(1,WTLENG)
!*
!*     Compute weights for 2 degree 3, and 2 degree 1 rules
!*
      W(2,WTLENG) = 1/(6*LAM1)**2
      W(2,2) = 1/(6*LAM3) - 2*(NDIM-1)*W(2,WTLENG)*LAM1/LAM3
      W(3,2) = ( 1/FIVE - LAM4/3 )/( 2*LAM2*(LAM2-LAM4) )
      W(3,4) = ( 1/FIVE - LAM2/3 )/( 2*LAM4*(LAM4-LAM2) )
      W(4,WTLENG) = ONE/(2*NDIM*(NDIM-1))
      W(5,2) = ONE/(4*NDIM)
      W(5,4) = ONE/(4*NDIM)
!*
!*     Set generator values
!*
      LAM1 = SQRT(LAM1)
      LAM2 = SQRT(LAM2)
      LAM3 = SQRT(LAM3)
      LAM4 = SQRT(LAM4)
      G(1,WTLENG) = LAM1
      G(2,WTLENG) = LAM1
      G(1,WTLENG-1) = LAM1
      G(1,2) = LAM2
      G(1,3) = LAM3
      G(1,4) = LAM4
!*
!*     Compute final weight values.
!*
      DO J = 1,5
         W(J,1) = TWONDM
         DO I = 2,WTLENG
            W(J,I) = TWONDM*W(J,I)
            W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
         END DO
      END DO
      CALL RULNRM ( WTLENG, 5, RULPTS, W )
!*
!****END BSRL05
!*
      END

      SUBROUTINE BSRL09( NDIM, WTLENG, W, G )
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: BSRL09
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: BSRL09
#endif
#if 0
****BEGIN PROLOGUE BSRL09
****KEYWORDS basic integration rule, degree 9
****PURPOSE  To initialize a degree 9 basic rule and null rules.
****AUTHOR   
*            Alan Genz, Department of Mathematics, Washington State
*            University, Pullman, WA 99163-3113, USA
*              Email: alangenz@wsu.edu
****LAST MODIFICATION 88-05-20
****DESCRIPTION  BSRL09 initializes a degree 9 integration rule,
*            two degree 7 rules, one degree 5 rule and one
*            degree 3 rule for the hypercube [-1,1]**NDIM.
*
*   ON ENTRY
*
*   NDIM   Integer.
*          Number of variables.
*
*   ON EXIT
*
*   WTLENG Integer.
*          The number of weights in each of the rules.
*   W      Real array of dimension (5,WTLENG).
*          The weights for the basic rules.
*          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
*   G      Real array of dimension (NDIM, WTLENG).
*          The fully symmetric sum generators for the rules.
*          G(1, J), ..., G(NDIM, J) are the are the generators for the
*          points associated with the Jth weights.
*
****REFERENCES A. Genz and A. Malik,
*             "An Imbedded Family of Fully Symmetric Numerical
*              Integration Rules",
*              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
****ROUTINES CALLED NONE
****END PROLOGUE BSRL09
*
*   Global variables
*
#endif
      integer(kind=4) ::  NDIM,WTLENG
      real(kind=8) ::  W(5,*),G(NDIM,*)
!*
!*   Local Variables
!*
      real(kind=8) ::  RATIO,LAM0,LAM1,LAM2,LAM3,LAMP,TWONDM
      integer(kind=4) ::  RULPTS(9)
      integer(kind=4) ::  I,J
!*
!****FIRST EXECUTABLE STATEMENT BSRL09
!*
!*
!*     Initialize generators, weights and RULPTS
!*
      WTLENG = 9
      IF ( NDIM .EQ. 2 ) WTLENG = 8
      DO J = 1,WTLENG
         DO I = 1,NDIM
            G(I,J) = 0
         END DO
         DO I = 1,5
            W(I,J) = 0
         END DO
         RULPTS(J) = 2*NDIM
      END DO
      TWONDM = 2**NDIM
      RULPTS(WTLENG) = TWONDM
      IF ( NDIM .GT. 2 ) RULPTS(8) = (4*NDIM* (NDIM-1)* (NDIM-2))/3
      RULPTS(7) = 4*NDIM* (NDIM-1)
      RULPTS(6) = 2*NDIM* (NDIM-1)
      RULPTS(1) = 1
!*
!*     Compute squared generator parameters
!*
      LAM0 = 0.4707_8
      LAM1 = 4/ (15-5/LAM0)
      RATIO = (1-LAM1/LAM0)/27
      LAM2 = (5-7*LAM1-35*RATIO)/ (7-35*LAM1/3-35*RATIO/LAM0)
      RATIO = RATIO* (1-LAM2/LAM0)/3
      LAM3 = (7-9* (LAM2+LAM1)+63*LAM2*LAM1/5-63*RATIO)/
     &       (9-63* (LAM2+LAM1)/5+21*LAM2*LAM1-63*RATIO/LAM0)
      LAMP = 0.0625_8
*
*     Compute degree 9 rule weights
*
      W(1,WTLENG) = 1/ (3*LAM0)**4/TWONDM
      IF (NDIM.GT.2) W(1,8) = (1-1/ (3*LAM0))/ (6*LAM1)**3
      W(1,7) = (1-7* (LAM0+LAM1)/5+7*LAM0*LAM1/3)/
     &         (84*LAM1*LAM2* (LAM2-LAM0)* (LAM2-LAM1))
      W(1,6) = (1-7* (LAM0+LAM2)/5+7*LAM0*LAM2/3)/
     &         (84*LAM1*LAM1* (LAM1-LAM0)* (LAM1-LAM2)) -
     &         W(1,7)*LAM2/LAM1 - 2* (NDIM-2)*W(1,8)
      W(1,4) = (1-9* ((LAM0+LAM1+LAM2)/7- (LAM0*LAM1+LAM0*LAM2+
     &         LAM1*LAM2)/5)-3*LAM0*LAM1*LAM2)/
     &         (18*LAM3* (LAM3-LAM0)* (LAM3-LAM1)* (LAM3-LAM2))
      W(1,3) = (1-9* ((LAM0+LAM1+LAM3)/7- (LAM0*LAM1+LAM0*LAM3+
     &         LAM1*LAM3)/5)-3*LAM0*LAM1*LAM3)/
     &         (18*LAM2* (LAM2-LAM0)* (LAM2-LAM1)* (LAM2-LAM3)) -
     &         2* (NDIM-1)*W(1,7)
      W(1,2) = (1-9* ((LAM0+LAM2+LAM3)/7- (LAM0*LAM2+LAM0*LAM3+
     &         LAM2*LAM3)/5)-3*LAM0*LAM2*LAM3)/
     &         (18*LAM1* (LAM1-LAM0)* (LAM1-LAM2)* (LAM1-LAM3)) -
     &         2* (NDIM-1)* (W(1,7)+W(1,6)+ (NDIM-2)*W(1,8))
*
*     Compute weights for 2 degree 7, 1 degree 5 and 1 degree 3 rules
*
      W(2,WTLENG) = 1/ (108*LAM0**4)/TWONDM
      IF (NDIM.GT.2) W(2,8) = (1-27*TWONDM*W(2,9)*LAM0**3)/ (6*LAM1)**3
      W(2,7) = (1-5*LAM1/3-15*TWONDM*W(2,WTLENG)*LAM0**2* (LAM0-LAM1))/
     &          (60*LAM1*LAM2* (LAM2-LAM1))
      W(2,6) = (1-9* (8*LAM1*LAM2*W(2,7)+TWONDM*W(2,WTLENG)*LAM0**2))/
     &         (36*LAM1*LAM1) - 2*W(2,8)* (NDIM-2)
      W(2,4) = (1-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(2,
     &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
     &         (14*LAM3* (LAM3-LAM1)* (LAM3-LAM2))
      W(2,3) = (1-7* ((LAM1+LAM3)/5-LAM1*LAM3/3+TWONDM*W(2,
     &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM3)))/
     &         (14*LAM2* (LAM2-LAM1)* (LAM2-LAM3)) - 2* (NDIM-1)*W(2,7)
      W(2,2) = (1-7* ((LAM2+LAM3)/5-LAM2*LAM3/3+TWONDM*W(2,
     &         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAM3)))/
     &         (14*LAM1* (LAM1-LAM2)* (LAM1-LAM3)) -
     &         2* (NDIM-1)* (W(2,7)+W(2,6)+ (NDIM-2)*W(2,8))
      W(3,WTLENG) = 5/ (324*LAM0**4)/TWONDM
      IF (NDIM.GT.2) W(3,8) = (1-27*TWONDM*W(3,9)*LAM0**3)/ (6*LAM1)**3
      W(3,7) = (1-5*LAM1/3-15*TWONDM*W(3,WTLENG)*LAM0**2* (LAM0-LAM1))/
     &          (60*LAM1*LAM2* (LAM2-LAM1))
      W(3,6) = (1-9* (8*LAM1*LAM2*W(3,7)+TWONDM*W(3,WTLENG)*LAM0**2))/
     &         (36*LAM1*LAM1) - 2*W(3,8)* (NDIM-2)
      W(3,5) = (1-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(3,
     &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
     &         (14*LAMP* (LAMP-LAM1)* (LAMP-LAM2))
      W(3,3) = (1-7* ((LAM1+LAMP)/5-LAM1*LAMP/3+TWONDM*W(3,
     &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAMP)))/
     &         (14*LAM2* (LAM2-LAM1)* (LAM2-LAMP)) - 2* (NDIM-1)*W(3,7)
      W(3,2) = (1-7* ((LAM2+LAMP)/5-LAM2*LAMP/3+TWONDM*W(3,
     &         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAMP)))/
     &         (14*LAM1* (LAM1-LAM2)* (LAM1-LAMP)) -
     &         2* (NDIM-1)* (W(3,7)+W(3,6)+ (NDIM-2)*W(3,8))
      W(4,WTLENG) = 2/ (81*LAM0**4)/TWONDM
      IF (NDIM.GT.2) W(4,8) = (2-27*TWONDM*W(4,9)*LAM0**3)/ (6*LAM1)**3
      W(4,7) = (2-15*LAM1/9-15*TWONDM*W(4,WTLENG)*LAM0* (LAM0-LAM1))/
     &         (60*LAM1*LAM2* (LAM2-LAM1))
      W(4,6) = (1-9* (8*LAM1*LAM2*W(4,7)+TWONDM*W(4,WTLENG)*LAM0**2))/
     &         (36*LAM1*LAM1) - 2*W(4,8)* (NDIM-2)
      W(4,4) = (2-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(4,
     &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
     &         (14*LAM3* (LAM3-LAM1)* (LAM3-LAM2))
      W(4,3) = (2-7* ((LAM1+LAM3)/5-LAM1*LAM3/3+TWONDM*W(4,
     &         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM3)))/
     &         (14*LAM2* (LAM2-LAM1)* (LAM2-LAM3)) - 2* (NDIM-1)*W(4,7)
      W(4,2) = (2-7* ((LAM2+LAM3)/5-LAM2*LAM3/3+TWONDM*W(4,
     &         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAM3)))/
     &         (14*LAM1* (LAM1-LAM2)* (LAM1-LAM3)) -
     &         2* (NDIM-1)* (W(4,7)+W(4,6)+ (NDIM-2)*W(4,8))
      W(5,2) = 1/ (6*LAM1)
!*
!*     Set generator values
!*
      LAM0 = SQRT(LAM0)
      LAM1 = SQRT(LAM1)
      LAM2 = SQRT(LAM2)
      LAM3 = SQRT(LAM3)
      LAMP = SQRT(LAMP)
      DO I = 1,NDIM
         G(I,WTLENG) = LAM0
      END DO
      IF ( NDIM .GT. 2 ) THEN
          G(1,8) = LAM1
          G(2,8) = LAM1
          G(3,8) = LAM1
      END IF
      G(1,7) = LAM1
      G(2,7) = LAM2
      G(1,6) = LAM1
      G(2,6) = LAM1
      G(1,5) = LAMP
      G(1,4) = LAM3
      G(1,3) = LAM2
      G(1,2) = LAM1
!*
!*     Compute final weight values.
!*
      DO J = 1,5
         W(J,1) = TWONDM
         DO I = 2,WTLENG
            W(J,I) = TWONDM*W(J,I)
            W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
         END DO
      END DO
      CALL RULNRM ( WTLENG, 5, RULPTS, W )
!*
!****END BSRL09
!*
      END

      SUBROUTINE RULNRM( LENRUL, NUMNUL, RULPTS, W )
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: RULNRM
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: RULNRM
#endif
      integer(kind=4) ::  LENRUL, NUMNUL, I, J, K, RULPTS(*)
      real(kind=8) ::  ALPHA, NORMCF, NORMNL, W(NUMNUL,*)
!*
!*     Compute orthonormalized null rules.
!*
      NORMCF = 0
      DO I = 1, LENRUL
         NORMCF = NORMCF + RULPTS(I)*W(1,I)*W(1,I)
      END DO
      DO K = 2, NUMNUL
         DO I = 1, LENRUL
            W(K,I) = W(K,I) - W(1,I)
         END DO
         DO J = 2, K-1
            ALPHA = 0
            DO I = 1, LENRUL
               ALPHA = ALPHA + RULPTS(I)*W(J,I)*W(K,I)
            END DO
            ALPHA = -ALPHA/NORMCF
            DO I = 1, LENRUL
               W(K,I) = W(K,I) + ALPHA*W(J,I)
            END DO
         END DO
         NORMNL = 0
         DO I = 1, LENRUL
            NORMNL = NORMNL + RULPTS(I)*W(K,I)*W(K,I)
         END DO
         ALPHA = SQRT( NORMCF/NORMNL )
         DO I = 1, LENRUL
            W(K,I) = ALPHA*W(K,I)
         END DO
      END DO
      END
      SUBROUTINE SDFFER( NDIM, CENTER, HWIDTH, NF, FUNSUB, 
     &                   X, WORK, DIVAXN )
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: SDFERR
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: SDFERR
#endif
!*
!*     Compute second differences and subdivision axis
!*
      EXTERNAL FUNSUB
      integer(kind=4) ::  I, J, NDIM, NF, DIVAXN
      real(kind=8) ::  CENTER(*), HWIDTH(*), X(*), WORK(NF,*)
      real(kind=8) ::  FRTHDF, DIFMAX, DIFSUM
      DIVAXN = 1
      DO I = 1,NDIM 
         IF ( HWIDTH(I) .GT. HWIDTH(DIVAXN) ) DIVAXN = I
         X(I) = CENTER(I)
      END DO
      DIFMAX = 0
      CALL FUNSUB( NDIM, CENTER, NF, WORK(1,3) ) 
      DO I = 1,NDIM
         X(I) = CENTER(I) - 4*HWIDTH(I)/5
         CALL FUNSUB( NDIM, X, NF, WORK(1,1) )
         X(I) = CENTER(I) - 2*HWIDTH(I)/5
         CALL FUNSUB( NDIM, X, NF, WORK(1,2) )
         X(I) = CENTER(I) + 2*HWIDTH(I)/5
         CALL FUNSUB( NDIM, X, NF, WORK(1,4) )
         X(I) = CENTER(I) + 4*HWIDTH(I)/5
         CALL FUNSUB( NDIM, X, NF, WORK(1,5) )
         X(I) = CENTER(I)
         DIFSUM = 0
         DO J = 1,NF
            FRTHDF = ABS( WORK(J,1) - 4*WORK(J,2) + 6*WORK(J,3) 
     &                              - 4*WORK(J,4) +   WORK(J,5) ) 
!*
!*     Ignore differences below roundoff
!*     
            IF ( ABS(WORK(J,3)) + FRTHDF/8 .GT. ABS(WORK(J,3)) ) 
     &           DIFSUM = DIFSUM + FRTHDF
         END DO
         IF ( DIFSUM .GT. DIFMAX ) THEN
            DIFMAX = DIFSUM
            DIVAXN = I
         END IF
      END DO
      END
      SUBROUTINE BASRUL(NDIM, CENTER, HWIDTH, WTLENG, G, W,
     &     NUMFUN, FUNSUB, X, NULL, BASVAL, RGNERR, GREAT)
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: BASRUL
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: BASRUL
#endif
#if 0
****BEGIN PROLOGUE BASRUL
****KEYWORDS basic numerical integration rule
****PURPOSE  To compute basic integration rule values.
****AUTHOR   
*            Alan Genz, Department of Mathematics, Washington State
*            University, Pullman, WA 99164-3113, USA
*              Email: alangenz@wsu.edu
****LAST MODIFICATION 93-08-20
****DESCRIPTION BASRUL computes basic integration rule values for a
*            vector of integrands over a hyper-rectangular region.
*            These are estimates for the integrals. 
*
*   ON ENTRY
*
*   NDIM   Integer.
*          Number of variables.
*   CENTER Real array of dimension NDIM.
*          The coordinates for the center of the region.
*   HWIDTH Real Array of dimension NDIM.
*          HWIDTH(I) is half of the width of dimension I of the region.
*   WTLENG Integer.
*          The number of weights in the basic integration rule.
*   G      Real array of dimension (NDIM,WTLENG).
*          The fully symmetric sum generators for the rules.
*          G(1,J), ..., G(NDIM,J) are the are the generators for the
*          points associated with the Jth weights.
*   W      Real array of dimension (5,WTLENG).
*          The weights for the basic and null rules.
*          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
*          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
*   NUMFUN Integer.
*          Number of components for the vector integrand.
*   FUNSUB Externally declared subroutine.
*          For computing the components of the integrand at a point X.
*          It must have parameters (NDIM,X,NUMFUN,FUNVLS).
*           Input Parameters:
*            X      Real array of dimension NDIM.
*                   Defines the evaluation point.
*            NDIM   Integer.
*                   Number of variables for the integrand.
*            NUMFUN Integer.
*                   Number of components for the vector integrand.
*           Output Parameters:
*            FUNVLS Real array of dimension NUMFUN.
*                   The components of the integrand at the point X.
*   X      Real Array of dimension NDIM.
*          A work array.
*   NULL   Real array of dimension (NUMFUN, 5)
*          A work array.
*
*   ON RETURN
*
*   BASVAL Real array of dimension NUMFUN.
*          The values for the basic rule for each component
*          of the integrand.
*   RGNERR Real array of dimension NUMFUN.
*          The error estimates for each component of the integrand.
*   GREAT  Real maximum error for RGNERR.
*
****ROUTINES CALLED: FULSUM, FUNSUB, TWONRM
*
****END PROLOGUE BASRUL
*
*   Global variables.
*
#endif
      EXTERNAL FUNSUB
      integer(kind=4) ::  WTLENG, NUMFUN, NDIM, NUMNUL
      PARAMETER ( NUMNUL = 4 )
      real(kind=8) ::  GREAT, CENTER(*), X(*), HWIDTH(*), BASVAL(*),
     &     RGNERR(*), NULL(NUMFUN,*), W(5,*), G(NDIM,*)
!*
!*   Local variables.
!*
      real(kind=8) ::  RGNVOL, RGNCMP, RGNCPT, TWONRM
      integer(kind=4) ::  I,J,K
!*
!****FIRST EXECUTABLE STATEMENT BASRUL
!*
!*
      RGNVOL = 1
      DO I = 1,NDIM
         RGNVOL = RGNVOL*HWIDTH(I)
      END DO
      DO J = 1,NUMFUN
         BASVAL(J) = 0
         DO K = 1,NUMNUL
            NULL(J,K) = 0
         END DO
      END DO
!*     
!*    Finish computing the rule values.
!*
      DO I = 1,WTLENG
         CALL FULSUM( NDIM, CENTER, HWIDTH, X, G(1,I), NUMFUN, FUNSUB,
     &                RGNERR, NULL(1,5) )
         DO J = 1,NUMFUN
            BASVAL(J) = BASVAL(J) + W(1,I)*RGNERR(J)
            DO K = 1,NUMNUL
               NULL(J,K) = NULL(J,K) + W(K+1,I)*RGNERR(J)
            END DO
         END DO
      END DO
!*
!*    Compute errors.
!*
      GREAT = 0
      DO J = 1,NUMFUN
         BASVAL(J) = RGNVOL*BASVAL(J)
         RGNERR(J) = TWONRM( NULL(J,1), NULL(J,2) )
         RGNCMP    = TWONRM( NULL(J,2), NULL(J,3) )
         RGNCPT    = TWONRM( NULL(J,3), NULL(J,4) )
         IF ( 4*RGNERR(J) .LT. RGNCMP .AND. 2*RGNCMP .LT. RGNCPT ) 
     &        RGNERR(J) = RGNERR(J)/2 
         IF ( 2*RGNERR(J) .GT. RGNCMP ) 
     &        RGNERR(J) = MAX( RGNERR(J), RGNCMP ) 
         RGNERR(J) = RGNVOL*RGNERR(J)*( 1 + ABS( BASVAL(J) ) )
         GREAT = GREAT + RGNERR(J)
      END DO
!*
!****END BASRUL
!*
      END
*
      real(kind=8)  FUNCTION TWONRM( X, Y )
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: TWONRM
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: TWONORM
      !DIR$ ATTRIBUTES INLINE :: TWONRM
#endif
      real(kind=8) ::  X, Y, ABX, ABY, SQTTWO
      PARAMETER ( SQTTWO = 1.4142135623730950488_8)
      ABX = ABS(X)
      ABY = ABS(Y)
      IF ( ABX .GT. ABY ) THEN
         TWONRM = ABX*SQRT( 1 + ( ABY/ABX )**2 )
      ELSE IF ( ABY .GT. ABX ) THEN
         TWONRM = ABY*SQRT( 1 + ( ABX/ABY )**2 )
      ELSE
         TWONRM = ABY*SQTTWO         
      END IF
      END

      SUBROUTINE FULSUM(NDIM,CENTER,HWIDTH,X,G,NUMFUN,FUNSUB,FULSMS,
     &                  FUNVLS)
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: FULSUM
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: FULSUM
#endif
#if 0
****BEGIN PROLOGUE FULSUM
****KEYWORDS fully symmetric sum
****PURPOSE  To compute fully symmetric basic rule sums
****AUTHOR   
*            Alan Genz, Department of Mathematics, Washington State
*            University, Pullman, WA 99164-3113, USA
*              Email: alangenz@wsu.edu
****LAST MODIFICATION 88-04-08
****DESCRIPTION FULSUM computes a fully symmetric sum for a vector
*            of integrand values over a hyper-rectangular region.
*            The sum is fully symmetric with respect to the center of
*            the region and is taken over all sign changes and
*            permutations of the generators for the sum.
*
*   ON ENTRY
*
*   NDIM   Integer.
*          Number of variables.
*   CENTER Real array of dimension NDIM.
*          The coordinates for the center of the region.
*   HWIDTH Real Array of dimension NDIM.
*          HWIDTH(I) is half of the width of dimension I of the region.
*   X      Real Array of dimension NDIM.
*          A work array.
*   G      Real Array of dimension NDIM.
*          The generators for the fully symmetric sum. These MUST BE
*          non-negative and non-increasing.
*   NUMFUN Integer.
*          Number of components for the vector integrand.
*   FUNSUB Externally declared subroutine.
*          For computing the components of the integrand at a point X.
*          It must have parameters (NDIM, X, NUMFUN, FUNVLS).
*           Input Parameters:
*            X      Real array of dimension NDIM.
*                   Defines the evaluation point.
*            NDIM   Integer.
*                   Number of variables for the integrand.
*            NUMFUN Integer.
*                   Number of components for the vector integrand.
*           Output Parameters:
*            FUNVLS Real array of dimension NUMFUN.
*                   The components of the integrand at the point X.
*   ON RETURN
*
*   FULSMS Real array of dimension NUMFUN.
*          The values for the fully symmetric sums for each component
*          of the integrand.
*   FUNVLS Real array of dimension NUMFUN.
*          A work array.
*
****ROUTINES CALLED: FUNSUB
*
****END PROLOGUE FULSUM
*
*   Global variables.
*
#endif
      EXTERNAL FUNSUB
      integer(kind=4) ::  NDIM,NUMFUN
      real(kind=8) ::  CENTER(NDIM),HWIDTH(NDIM),X(NDIM),G(NDIM),
     &                 FULSMS(NUMFUN),FUNVLS(NUMFUN)
!*
!*   Local variables.
!*
      integer(kind=4) ::  IXCHNG,LXCHNG,I,J,L
      real(kind=8) ::  GL,GI
!*
!****FIRST EXECUTABLE STATEMENT FULSUM
!*
      DO J = 1, NUMFUN
         FULSMS(J) = 0
      END DO
!*     
!*     Compute centrally symmetric sum for permutation of G
!*     
 10   DO I = 1, NDIM
         X(I) = CENTER(I) + G(I)*HWIDTH(I)
      END DO
 20   CALL FUNSUB( NDIM, X, NUMFUN, FUNVLS )
      DO J = 1,NUMFUN
         FULSMS(J) = FULSMS(J) + FUNVLS(J)
      END DO
      DO I = 1, NDIM
         G(I) = - G(I)
         X(I) = CENTER(I) + G(I)*HWIDTH(I)
         IF ( G(I) .LT. 0 ) GO TO 20
      END DO
!*     
!*     Find next distinct permuation of G and loop back for next sum.
!*     Permutations are generated in reverse lexicographic order.
!*     
      DO I = 2, NDIM
         IF ( G(I-1) .GT. G(I) ) THEN
            GI = G(I)
            IXCHNG = I - 1
            DO L = 1, (I-1)/2
               GL = G(L)
               G(L) = G(I-L)
               G(I-L) = GL
               IF ( GL .LE. GI ) IXCHNG = IXCHNG - 1
               IF ( G(L).GT. GI ) LXCHNG = L
            END DO
            IF ( G(IXCHNG) .LE. GI ) IXCHNG = LXCHNG
            G(I) = G(IXCHNG)
            G(IXCHNG) = GI
            GO TO 10
         END IF
      END DO
!*     
!*     Restore original order to generators
!*
      DO I = 1,NDIM/2
          GI = G(I)
          G(I) = G(NDIM-I+1)
          G(NDIM-I+1) = GI
       END DO
!*
!****END FULSUM
!*
      END
      SUBROUTINE ADONEV( NF, A, B, MNS, MXS, F, AB, RE, IR, 
     &     RESULT, ABSERR, IPTS, IM, IFAIL, S, E, C, H, EMX, WORK )
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ADONEV
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: ADONEV
#endif
!*
!*     One Dimensional Adaptive Integration Routine
!*
      EXTERNAL F
      real(kind=8) ::  A, B, AB, RE, RESULT(*), ABSERR(*), WORK(*)
      integer(kind=4) ::  I, IR, NF, IM, IP, IPTS, J, IFAIL, MNS, MXS
      real(kind=8) ::  C(*), H(*), EMX(*), S(NF,*), E(NF,*)
      IPTS = 0
      IF ( IR .EQ. 0 ) THEN
         C(1) = ( B + A )/2
         H(1) = ( B - A )/2
         CALL KRNRDV( C(1), H(1), NF, F, S(1,1), E(1,1), EMX(1), WORK ) 
         DO I = 1,NF
            ABSERR(I) = E(I,1)
            RESULT(I) = S(I,1)
         END DO
         IPTS = IPTS + 15
      END IF
 10   IFAIL = 0
      DO I = 1,NF
         IF ( ABSERR(I) .GT. MAX( AB, RE*ABS( RESULT(I) ) ) ) IFAIL = 1
      END DO
      IF ( IM .LT. MNS .OR. IFAIL .EQ. 1 .AND. IM .LT. MXS ) THEN
         IP = 1
         DO I = 1, IM
            IF ( EMX(I) .GT. EMX(IP) ) IP = I
         END DO
         IM = IM + 1
         H(IP) = H(IP)/2
         H(IM) = H(IP)
         C(IM) = C(IP) + H(IP)
         C(IP) = C(IP) - H(IP)
         CALL KRNRDV( C(IP),H(IP), NF,F, S(1,IP),E(1,IP),EMX(IP), WORK ) 
         CALL KRNRDV( C(IM),H(IM), NF,F, S(1,IM),E(1,IM),EMX(IM), WORK ) 
         IPTS = IPTS + 30
         DO I = 1,NF
            RESULT(I) = 0
            ABSERR(I) = 0
         END DO
         DO J = 1,IM
            DO I = 1,NF
               RESULT(I) = RESULT(I) + S(I,J)
               ABSERR(I) = ABSERR(I) + E(I,J)
            END DO
         END DO
         GO TO 10
      ENDIF
      END


      SUBROUTINE KRNRDV( C, H, NF, F, RESULT, ABSERR, ERROR, FUNS )
#if defined(__INTEL_COMPILER) || defined(__ICC)
      !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: KRNRDV
      !DIR$ OPTIMIZE : 3
      !DIR$ CODE_ALIGN : 32 :: KRNRDV
#endif
!*
!*     Kronrod Rule
!*
      EXTERNAL F
      real(kind=8) ::  C, H, ABSERR(*), FUNS(*), RESULT(*), ERROR, FS
      integer(kind=4) ::  I, J, N, NF
      PARAMETER ( N = 7 )
!*
      real(kind=8) ::  WG(0:N), WGK(0:N), XGK(0:N) 
      SAVE WG, WGK, XGK
!*
!*           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1)
!*           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSE AND THEIR 
!*           CORRESPONDING WEIGHTS ARE GIVEN.
!*
!*           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE 
!*                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
!*                    GAUSS RULE
!*                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!*                    ADDED TO THE 7-POINT GAUSS RULE
!*
!*           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
!*
!*           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
!*
      DATA WG(0) /0.4179591836 7346938775 5102040816 327 D0/
      DATA WG(2) /0.1294849661 6886969327 0611432679 082 D0/
      DATA WG(4) /0.2797053914 8927666790 1467771423 780 D0/
      DATA WG(6) /0.3818300505 0511894495 0369775488 975 D0/
      DATA WG(1), WG(3), WG(5), WG(7) /4*0D0/

      DATA XGK(0) /0.0000000000 0000000000 0000000000 000 D0/
      DATA XGK(1) /0.9914553711 2081263920 6854697526 329 D0/
      DATA XGK(2) /0.9491079123 4275852452 6189684047 851 D0/
      DATA XGK(3) /0.8648644233 5976907278 9712788640 926 D0/
      DATA XGK(4) /0.7415311855 9939443986 3864773280 788 D0/
      DATA XGK(5) /0.5860872354 6769113029 4144838258 730 D0/
      DATA XGK(6) /0.4058451513 7739716690 6606412076 961 D0/
      DATA XGK(7) /0.2077849550 0789846760 0689403773 245 D0/

      DATA WGK(0) /0.2094821410 8472782801 2999174891 714 D0/
      DATA WGK(1) /0.0229353220 1052922496 3732008058 970 D0/
      DATA WGK(2) /0.0630920926 2997855329 0700663189 204 D0/
      DATA WGK(3) /0.1047900103 2225018383 9876322541 518 D0/
      DATA WGK(4) /0.1406532597 1552591874 5189590510 238 D0/
      DATA WGK(5) /0.1690047266 3926790282 6583426598 550 D0/
      DATA WGK(6) /0.1903505780 6478540991 3256402421 014 D0/
      DATA WGK(7) /0.2044329400 7529889241 4161999234 649 D0/
!*
!*           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
!*           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
!*
      CALL F( 1, C, NF, FUNS )
      DO I = 1,NF
         RESULT(I) = WGK(0)*FUNS(I)
         ABSERR(I) =  WG(0)*FUNS(I)
      END DO
      DO J = 1,N
         CALL F( 1, C + H*XGK(J), NF, FUNS )
         CALL F( 1, C - H*XGK(J), NF, FUNS(NF+1) )
         DO I = 1,NF
            FS = FUNS(I) + FUNS(NF+I)
            RESULT(I) = RESULT(I) + WGK(J)*FS 
            ABSERR(I) = ABSERR(I) +  WG(J)*FS
         END DO
      END DO
      ERROR = 0
      DO I = 1,NF
         ABSERR(I) = ABS( H*( RESULT(I) - ABSERR(I) ) )
         ERROR = MAX( ERROR, ABSERR(I) )
         RESULT(I) = H*RESULT(I)
      END DO
      END

!======================================================================
!   Example program
!======================================================================

#if 0


   EXTERNAL FTEST
*   
*        ADAPT EXAMPLE TEST PROGRAM
*
      INTEGER N, MNCLS, MXCLS, IFAIL, NV, ND, NF, LW, KEY
      PARAMETER (ND = 3, NF = 5, LW = 10000 )
      DOUBLE PRECISION A(ND), B(ND), WRKSTR(LW), ABSEST(NF), FINEST(NF)
      DO N = 1, ND
        A(N) = 0
        B(N) = 1
      END DO
      KEY = 1
      MNCLS = 0
      MXCLS = 1000
      CALL ADAPT(ND, NF, A, B, MNCLS, MXCLS, FTEST, 1D-5, 0, 
     *           KEY, LW, 0, FINEST, ABSEST, NV, IFAIL, WRKSTR)
      PRINT 99999, KEY, NV, IFAIL
99999 FORMAT (5X, 'ADAPT TEST RESULTS, WITH KEY =', I2,
     *    /'     FTEST CALLS = ', I5, ', IFAIL = ', I2, 
     *    /'    N  ESTIMATED ERROR   INTEGRAL')
      DO N = 1 , NF
        PRINT '(3X, I2, 2E15.6)', N, ABSEST(N), FINEST(N)
      END DO
      KEY = 2
      CALL ADAPT(ND, NF, A, B, MNCLS, 10*MXCLS, FTEST, 1D-5, 0, 
     *           KEY, LW, 1, FINEST, ABSEST, NV, IFAIL, WRKSTR)
      PRINT *, '      CONTINUATION OF CALCULATION'
      PRINT 99999, KEY, NV, IFAIL
      DO N = 1 , NF
        PRINT '(3X, I2, 2E15.6)', N, ABSEST(N), FINEST(N)
      END DO
      END
      SUBROUTINE FTEST(NDIM, Z, NFUN, F)
      INTEGER N, NDIM, NFUN
      DOUBLE PRECISION Z(NDIM), F(NFUN), SUM
      DO N = 1, NDIM
        F(N) = Z(N)**(N+2)
      END DO
      F(4) = COS(5*Z(3)**2)
      F(5) = EXP(5*Z(2)**2)
      END



#endif
