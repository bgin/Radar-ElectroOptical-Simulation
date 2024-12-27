C
C     file cblktri.f
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
C     SUBROUTINE CBLKTR (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,
C    +                   IERROR)
C
C                                                                       
C DIMENSION OF           AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N)
C ARGUMENTS
C                                                                       
C LATEST REVISION        JUNE 2004
C                                                                       
C PURPOSE                CBLKTR SOLVES A SYSTEM OF LINEAR EQUATIONS     
C                        OF THE FORM                                    
C                                                                       
C                        AN(J)*X(I,J-1) + AM(I)*X(I-1,J) +              
C                        (BN(J)+BM(I))*X(I,J) + CN(J)*X(I,J+1) +        
C                        CM(I)*X(I+1,J) = Y(I,J)                        
C                                                                       
C                        FOR I = 1,2,...,M  AND  J = 1,2,...,N.         
C                                                                       
C                        I+1 AND I-1 ARE EVALUATED MODULO M AND         
C                        J+1 AND J-1 MODULO N, I.E.,                    
C                                                                       
C                        X(I,0) = X(I,N),  X(I,N+1) = X(I,1),           
C                        X(0,J) = X(M,J),  X(M+1,J) = X(1,J).           
C                                                                       
C                        THESE EQUATIONS USUALLY RESULT FROM THE        
C                        DISCRETIZATION OF SEPARABLE ELLIPTIC           
C                        EQUATIONS.  BOUNDARY CONDITIONS MAY BE         
C                        DIRICHLET, NEUMANN, OR PERIODIC.               
C                                                                       
C                        CBLKTRI IS A COMPLEX VERSION OF PACKAGE        
C                        BLKTRI ON ULIB.                                
C                                                                       
C USAGE                  CALL CBLKTR (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,    
C                                     CM,IDIMY,Y,IERROR,W)              
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON INPUT               IFLG                                           
C                                                                       
C                          = 0  INITIALIZATION ONLY.                    
C                               CERTAIN QUANTITIES THAT DEPEND ON NP,   
C                               N, AN, BN, AND CN ARE COMPUTED AND      
C                               STORED IN THE DERIVED DATA TYPE W
C                                                                       
C                          = 1  THE QUANTITIES THAT WERE COMPUTED       
C                               IN THE INITIALIZATION ARE USED          
C                               TO OBTAIN THE SOLUTION X(I,J).          
C                                                                       
C                               NOTE:                                   
C                               A CALL WITH IFLG=0 TAKES                
C                               APPROXIMATELY ONE HALF THE TIME         
C                               AS A CALL WITH IFLG = 1.                
C                               HOWEVER, THE INITIALIZATION DOES        
C                               NOT HAVE TO BE REPEATED UNLESS NP,      
C                               N, AN, BN, OR CN CHANGE.                
C                                                                       
C                        NP                                             
C                          = 0  IF AN(1) AND CN(N) ARE NOT ZERO,        
C                               WHICH CORRESPONDS TO PERIODIC           
C                               BOUNARY CONDITIONS.                     
C                                                                       
C                          = 1  IF AN(1) AND CN(N) ARE ZERO.            
C                                                                       
C                        N                                              
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.   
C                          N MUST BE GREATER THAN 4.                    
C                          THE OPERATION COUNT IS PROPORTIONAL TO       
C                          MNLOG2(N), HENCE N SHOULD BE SELECTED        
C                          LESS THAN OR EQUAL TO M.                     
C                                                                       
C                        AN,BN,CN                                       
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH N           
C                          THAT SPECIFY THE COEFFICIENTS IN THE         
C                          LINEAR EQUATIONS GIVEN ABOVE.                
C                                                                       
C                        MP                                             
C                          = 0  IF AM(1) AND CM(M) ARE NOT ZERO,        
C                               WHICH CORRESPONDS TO PERIODIC           
C                               BOUNDARY CONDITIONS.                    
C                                                                       
C                          = 1  IF AM(1) = CM(M) = 0  .                 
C                                                                       
C                        M                                              
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.   
C                           M MUST BE GREATER THAN 4.                   
C                                                                       
C                        AM,BM,CM                                       
C                          COMPLEX ONE-DIMENSIONAL ARRAYS OF LENGTH M   
C                          THAT SPECIFY THE COEFFICIENTS IN THE LINEAR  
C                          EQUATIONS GIVEN ABOVE.                       
C                                                                       
C                        IDIMY                                          
C                          THE ROW (OR FIRST) DIMENSION OF THE          
C                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS        
C                          IN THE PROGRAM CALLING CBLKTR.               
C                          THIS PARAMETER IS USED TO SPECIFY THE        
C                          VARIABLE DIMENSION OF Y.                     
C                          IDIMY MUST BE AT LEAST M.                    
C                                                                       
C                        Y                                              
C                          A COMPLEX TWO-DIMENSIONAL ARRAY THAT         
C                          SPECIFIES THE VALUES OF THE RIGHT SIDE OF    
C                          THE LINEAR SYSTEM OF EQUATIONS GIVEN ABOVE.  
C                          Y MUST BE DIMENSIONED Y(IDIMY,N) WITH        
C                          IDIMY .GE. M.                                
C                                                                       
C                        W
c                          A fortran 90 derived TYPE (fishworkspace) variable
c                          that must be declared by the user.  The first
c                          two declarative statements in the user program
c                          calling CBLKTRI must be:
c
c                               USE fish
c                               TYPE (fishworkspace) :: W
c
c                          The first statement makes the fishpack module
c                          defined in the file "fish.f" available to the
c                          user program calling CBLKTRI.  The second statement
c                          declares a derived type variable (defined in
c                          the module "fish.f") which is used internally
c                          in CBLKTRI to dynamically allocate real and complex
c                          work space used in solution.  An error flag
c                          (IERROR = 20) is set if the required work space
c                          allocation fails (for example if N,M are too large)
c                          Real and complex values are set in the components
c                          of W on a initial (IFLG=0) call to CBLKTRI.  These
c                          must be preserved on non-initial calls (IFLG=1)
c                          to CBLKTRI.  This eliminates redundant calculations
c                          and saves compute time.
c               ****       IMPORTANT!  The user program calling CBLKTRI should
c                          include the statement:
c
c                               CALL FISHFIN(W)
C
C                          after the final approximation is generated by
C                          CBLKTRI.  The will deallocate the real and complex
c                          work space of W.  Failure to include this statement
c                          could result in serious memory leakage.
c
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON OUTPUT              Y                                              
C                          CONTAINS THE SOLUTION X.                     
C                                                                       
C                        IERROR                                         
C                          AN ERROR FLAG THAT INDICATES INVALID         
C                          INPUT PARAMETERS.  EXCEPT FOR NUMBER ZER0,   
C                          A SOLUTION IS NOT ATTEMPTED.                 
C                                                                       
C                        = 0  NO ERROR.                                 
C                        = 1  M IS LESS THAN 5                          
C                        = 2  N IS LESS THAN 5                          
C                        = 3  IDIMY IS LESS THAN M.                     
C                        = 4  CBLKTR FAILED WHILE COMPUTING RESULTS     
C                             THAT DEPEND ON THE COEFFICIENT ARRAYS     
C                             AN, BN, CN.  CHECK THESE ARRAYS.          
C                        = 5  AN(J)*CN(J-1) IS LESS THAN 0 FOR SOME J.  
C                                                                       
C                             POSSIBLE REASONS FOR THIS CONDITION ARE   
C                             1. THE ARRAYS AN AND CN ARE NOT CORRECT   
C                             2. TOO LARGE A GRID SPACING WAS USED      
C                                IN THE DISCRETIZATION OF THE ELLIPTIC  
C                                EQUATION.                              
C                             3. THE LINEAR EQUATIONS RESULTED FROM A   
C                                PARTIAL DIFFERENTIAL EQUATION WHICH    
C                                WAS NOT ELLIPTIC.                      
C
C                          = 20 If the dynamic allocation of real and
C                               complex work space in the derived type
C                               (fishworkspace) variable W fails (e.g.,
c                               if N,M are too large for the platform used)
c
C                                                                       
C                                                                       
C SPECIAL CONDITIONS     THE ALGORITHM MAY FAIL IF ABS(BM(I)+BN(J))     
C                        IS LESS THAN ABS(AM(I))+ABS(AN(J))+            
C                        ABS(CM(I))+ABS(CN(J))                          
C                        FOR SOME I AND J. THE ALGORITHM WILL ALSO      
C                        FAIL IF AN(J)*CN(J-1) IS LESS THAN ZERO FOR    
C                        SOME J.                                        
C                        SEE THE DESCRIPTION OF THE OUTPUT PARAMETER    
C                        IERROR.                                        
C                                                                       
C                                                                       
C I/O                    NONE                                           
C                                                                       
C PRECISION              SINGLE                                         
C                                                                       
C REQUIRED LIBRARY       comf.f,fish.f
C FILES
C                                                                       
C LANGUAGE               FORTRAN 90
C                                                                       
C HISTORY                WRITTEN BY PAUL SWARZTRAUBER AT NCAR IN        
C                        THE EARLY 1970'S.  REWRITTEN AN RELEASED       
C                        ON NCAR'S PUBLIC SOFTWARE LIBRARIES IN         
C                        JANUARY, 1980. Revised in June 2004 by John
C                        Adams using Fortan 90 dynamically allocated
c                        space and derived data types to eliminate mixed
c                        mode conflicts in the earlier versions.
C                                                                       
C ALGORITHM              GENERALIZED CYCLIC REDUCTION                   
C                        (SEE REFERENCE BELOW)                          
C                                                                       
C PORTABILITY
C                        THE APPROXIMATE MACHINE ACCURACY IS COMPUTED   
C                        IN FUNCTION EPMACH                             
C                                                                       
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, 'EFFICIENT       
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF        
C                        ELLIPTIC EQUATIONS'                            
C                        NCAR TN/IA-109, JULY, 1975, 138 PP.            
C                                                                       
C                        SWARZTRAUBER P. N.,A DIRECT METHOD FOR         
C                        THE DISCRETE SOLUTION OF SEPARABLE             
C                        ELLIPTIC EQUATIONS, S.I.A.M.                   
C                        J. NUMER. ANAL.,11(1974) PP. 1136-1150.        
C                                                                       
C***********************************************************************
      SUBROUTINE CBLKTRI(IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, 
     1   IDIMY, Y, IERROR, W)
      USE fish
      Implicit none
      TYPE (fishworkspace) :: w
      EXTERNAL        PROC       ,PROCP      ,CPROC      ,CPROCP        
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: IFLG
      INTEGER , INTENT(IN) :: NP
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: MP
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER  :: IERROR
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      COMPLEX  :: AM(*)
      COMPLEX  :: BM(*)
      COMPLEX  :: CM(*)
      COMPLEX  :: Y(IDIMY,*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::M2,NH,NL,IWAH,IW1,IWBH,IW2,IW3,IWD,IWW,IWU,IRWK,ICWK
C-----------------------------------------------
C
C TEST M AND N FOR THE PROPER FORM
C
      NM = N
      M2 = M + M
      IERROR = 0
      IF (M - 5 < 0) THEN
         IERROR = 1
      ELSE
         IF (NM - 3 < 0) THEN
            IERROR = 2
         ELSE
            IF (IDIMY - M < 0) THEN
               IERROR = 3
            ELSE
               NH = N
               NPP = NP
               IF (NPP /= 0) THEN
                  NH = NH + 1
               ENDIF
               IK = 2
               K = 1
               IK = IK + IK
               K = K + 1
               DO WHILE(NH - IK > 0)
                  IK = IK + IK
                  K = K + 1
               END DO
               NL = IK
               IK = IK + IK
               NL = NL - 1
               IWAH = (K - 2)*IK + K + 6
               IF (NPP /= 0) THEN
                  IW1 = IWAH
                  IWBH = IW1 + NM
               ELSE
                  IWBH = IWAH + NM + NM
                  IW1 = IWBH
                  NM = NM - 1
               ENDIF
C
C SUBROUTINE COMP B COMPUTES THE ROOTS OF THE B POLYNOMIALS
C
               IF (IERROR == 0) THEN
                  IW2 = IW1 + M
                  IW3 = IW2 + M
                  IWD = IW3 + M
                  IWW = IWD + M
                  IWU = IWW + M
                  IF (IFLG == 0) THEN
                     IRWK = IW1 + 2*N
                     ICWK = IW1 + 6*M
                     CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
                     IF (IERROR /= 0) RETURN 
!     COMPUTE b poly roots (real and complex)
      call ccompb(NL,ierror,an,bn,cn,w%rew,w%cxw,w%rew(iwah),
     +            w%rew(iwbh))
                  ELSE
                     IF (MP /= 0) THEN
      CALL CBLKT1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,w%rew,w%cxw,
     +w%cxw(iw1),w%cxw(iw2),w%cxw(iw3),w%cxw(iwd),w%cxw(iww),
     +w%cxw(iwu),PROC,CPROC)
                     ELSE
      CALL CBLKT1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,w%rew,w%cxw,
     +w%cxw(iw1),w%cxw(iw2),w%cxw(iw3),w%cxw(iwd),w%cxw(iww),
     +w%cxw(iwu),PROCP,CPROCP)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      RETURN 
      END SUBROUTINE CBLKTRI


      SUBROUTINE CBLKT1(N, AN, BN, CN, M, AM, BM, CM, IDIMY, Y, B, BC, 
     1   W1, W2, W3, WD, WW, WU, PRDCT, CPRDCT)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      INTEGER  :: M
      INTEGER , INTENT(IN) :: IDIMY
      REAL  :: AN(*)
      REAL  :: BN(*)
      REAL  :: CN(*)
      REAL  :: B(*)
      COMPLEX  :: AM(*)
      COMPLEX  :: BM(*)
      COMPLEX  :: CM(*)
      COMPLEX  :: Y(IDIMY,*)
      COMPLEX  :: BC(*)
      COMPLEX  :: W1(*)
      COMPLEX  :: W2(*)
      COMPLEX  :: W3(*)
      COMPLEX  :: WD(*)
      COMPLEX  :: WW(*)
      COMPLEX  :: WU(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: KDO, L, IR, I2, I1, I3, I4, IRM1, IM2, NM2, IM3, NM3, 
     1   IM1, NM1, IF, I, IPI1, IPI2, IPI3, IDXC, NC, IDXA, NA, IP2, NP2
     2   , IP1, NP1, IP3, NP3, J, IZ, NZ, IZR, LL, IFD, IP, NP, IMI1, 
     3   IMI2
      REAL :: DUM
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
C
C CBLKT1 SOLVES THE LINEAR SYSTEM
C
C B  CONTAINS THE ROOTS OF ALL THE B POLYNOMIALS
C W1,W2,W3,WD,WW,WU  ARE ALL WORKING ARRAYS
C PRDCT IS EITHER PROCP OR PROC DEPENDING ON WHETHER THE BOUNDARY
C CONDITIONS IN THE M DIRECTION ARE PERIODIC OR NOT
C CPRDCT IS EITHER CPROCP OR CPROC WHICH ARE CALLED IF SOME OF THE ZEROS
C OF THE B POLYNOMIALS ARE COMPLEX.
C
C
C
C BEGIN REDUCTION PHASE
C
      KDO = K - 1
      DO L = 1, KDO
         IR = L - 1
         I2 = 2**IR
         I1 = I2/2
         I3 = I2 + I1
         I4 = I2 + I2
         IRM1 = IR - 1
         CALL CINDXB (I2, IR, IM2, NM2)
         CALL CINDXB (I1, IRM1, IM3, NM3)
         CALL CINDXB (I3, IRM1, IM1, NM1)
         CALL PRDCT (NM2, B(IM2), NM3, B(IM3), NM1, B(IM1), 0, DUM, Y(1,
     1      I2), W3, M, AM, BM, CM, WD, WW, WU)
         IF = 2**K
         DO I = I4, IF, I4
            IF (I - NM > 0) CYCLE 
            IPI1 = I + I1
            IPI2 = I + I2
            IPI3 = I + I3
            CALL CINDXC (I, IR, IDXC, NC)
            IF (I - IF >= 0) CYCLE 
            CALL CINDXA (I, IR, IDXA, NA)
            CALL CINDXB (I - I1, IRM1, IM1, NM1)
            CALL CINDXB (IPI2, IR, IP2, NP2)
            CALL CINDXB (IPI1, IRM1, IP1, NP1)
            CALL CINDXB (IPI3, IRM1, IP3, NP3)
            CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W3, 
     1         W1, M, AM, BM, CM, WD, WW, WU)
            IF (IPI2 - NM > 0) THEN
               W3(:M) = (0.,0.)
               W2(:M) = (0.,0.)
            ELSE
               CALL PRDCT (NP2, B(IP2), NP1, B(IP1), NP3, B(IP3), 0, DUM
     1            , Y(1,IPI2), W3, M, AM, BM, CM, WD, WW, WU)
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W3
     1            , W2, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            Y(:M,I) = W1(:M) + W2(:M) + Y(:M,I)
         END DO
      END DO
      IF (NPP == 0) THEN
         IF = 2**K
         I = IF/2
         I1 = I/2
         CALL CINDXB (I - I1, K - 2, IM1, NM1)
         CALL CINDXB (I + I1, K - 2, IP1, NP1)
         CALL CINDXB (I, K - 1, IZ, NZ)
         CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, Y(1,I)
     1      , W1, M, AM, BM, CM, WD, WW, WU)
         IZR = I
         W2(:M) = W1(:M)
         DO LL = 2, K
            L = K - LL + 1
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I = I2
            CALL CINDXC (I, IR, IDXC, NC)
            CALL CINDXB (I, IR, IZ, NZ)
            CALL CINDXB (I - I1, IR - 1, IM1, NM1)
            CALL CINDXB (I + I1, IR - 1, IP1, NP1)
            CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W1, 
     1         W1, M, AM, BM, CM, WD, WW, WU)
            W1(:M) = Y(:M,I) + W1(:M)
            CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, W1
     1         , W1, M, AM, BM, CM, WD, WW, WU)
         END DO
         L118: DO LL = 2, K
            L = K - LL + 1
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I4 = I2 + I2
            IFD = IF - I2
            DO I = I2, IFD, I4
               IF (I - I2 - IZR /= 0) CYCLE 
               IF (I - NM > 0) CYCLE  L118
               CALL CINDXA (I, IR, IDXA, NA)
               CALL CINDXB (I, IR, IZ, NZ)
               CALL CINDXB (I - I1, IR - 1, IM1, NM1)
               CALL CINDXB (I + I1, IR - 1, IP1, NP1)
               CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W2
     1            , W2, M, AM, BM, CM, WD, WW, WU)
               W2(:M) = Y(:M,I) + W2(:M)
               CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, 
     1            W2, W2, M, AM, BM, CM, WD, WW, WU)
               IZR = I
               IF (I - NM == 0) EXIT  L118
            END DO
         END DO L118
  119    CONTINUE
         Y(:M,NM+1) = Y(:M,NM+1) - CN(NM+1)*W1(:M) - AN(NM+1)*W2(:M)
         CALL CINDXB (IF/2, K - 1, IM1, NM1)
         CALL CINDXB (IF, K - 1, IP, NP)
         IF (NCMPLX /= 0) THEN
            CALL CPRDCT (NM + 1, BC(IP), NM1, B(IM1), 0, DUM, 0, DUM, Y(
     1         1,NM+1), Y(1,NM+1), M, AM, BM, CM, W1, W3, WW)
         ELSE
            CALL PRDCT (NM + 1, B(IP), NM1, B(IM1), 0, DUM, 0, DUM, Y(1,
     1         NM+1), Y(1,NM+1), M, AM, BM, CM, WD, WW, WU)
         ENDIF
         W1(:M) = AN(1)*Y(:M,NM+1)
         W2(:M) = CN(NM)*Y(:M,NM+1)
         Y(:M,1) = Y(:M,1) - W1(:M)
         Y(:M,NM) = Y(:M,NM) - W2(:M)
         DO L = 1, KDO
            IR = L - 1
            I2 = 2**IR
            I4 = I2 + I2
            I1 = I2/2
            I = I4
            CALL CINDXA (I, IR, IDXA, NA)
            CALL CINDXB (I - I2, IR, IM2, NM2)
            CALL CINDXB (I - I2 - I1, IR - 1, IM3, NM3)
            CALL CINDXB (I - I1, IR - 1, IM1, NM1)
            CALL PRDCT (NM2, B(IM2), NM3, B(IM3), NM1, B(IM1), 0, DUM, 
     1         W1, W1, M, AM, BM, CM, WD, WW, WU)
            CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W1, 
     1         W1, M, AM, BM, CM, WD, WW, WU)
            Y(:M,I) = Y(:M,I) - W1(:M)
         END DO
C
         IZR = NM
         L131: DO L = 1, KDO
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I3 = I2 + I1
            I4 = I2 + I2
            IRM1 = IR - 1
            DO I = I4, IF, I4
               IPI1 = I + I1
               IPI2 = I + I2
               IPI3 = I + I3
               IF (IPI2 - IZR /= 0) THEN
                  IF (I - IZR /= 0) CYCLE 
                  CYCLE  L131
               ENDIF
               CALL CINDXC (I, IR, IDXC, NC)
               CALL CINDXB (IPI2, IR, IP2, NP2)
               CALL CINDXB (IPI1, IRM1, IP1, NP1)
               CALL CINDXB (IPI3, IRM1, IP3, NP3)
               CALL PRDCT (NP2, B(IP2), NP1, B(IP1), NP3, B(IP3), 0, DUM
     1            , W2, W2, M, AM, BM, CM, WD, WW, WU)
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W2
     1            , W2, M, AM, BM, CM, WD, WW, WU)
               Y(:M,I) = Y(:M,I) - W2(:M)
               IZR = I
               CYCLE  L131
            END DO
         END DO L131
      ENDIF
C
C BEGIN BACK SUBSTITUTION PHASE
C
      DO LL = 1, K
         L = K - LL + 1
         IR = L - 1
         IRM1 = IR - 1
         I2 = 2**IR
         I1 = I2/2
         I4 = I2 + I2
         IFD = IF - I2
         DO I = I2, IFD, I4
            IF (I - NM > 0) CYCLE 
            IMI1 = I - I1
            IMI2 = I - I2
            IPI1 = I + I1
            IPI2 = I + I2
            CALL CINDXA (I, IR, IDXA, NA)
            CALL CINDXC (I, IR, IDXC, NC)
            CALL CINDXB (I, IR, IZ, NZ)
            CALL CINDXB (IMI1, IRM1, IM1, NM1)
            CALL CINDXB (IPI1, IRM1, IP1, NP1)
            IF (I - I2 <= 0) THEN
               W1(:M) = (0.,0.)
            ELSE
               CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), Y(
     1            1,IMI2), W1, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            IF (IPI2 - NM > 0) THEN
               W2(:M) = (0.,0.)
            ELSE
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), Y(
     1            1,IPI2), W2, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            W1(:M) = Y(:M,I) + W1(:M) + W2(:M)
            CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, W1
     1         , Y(1,I), M, AM, BM, CM, WD, WW, WU)
         END DO
      END DO
      RETURN 
      END SUBROUTINE CBLKT1


      REAL FUNCTION CBSRH (XLL, XRR, IZ, C, A, BH, F, SGN)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: IZ
      REAL , INTENT(IN) :: XLL
      REAL , INTENT(IN) :: XRR
      REAL  :: F
      REAL , INTENT(IN) :: SGN
      REAL  :: C(*)
      REAL  :: A(*)
      REAL  :: BH(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL :: R1, XL, XR, DX, X
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
C-----------------------------------------------
      XL = XLL
      XR = XRR
      DX = 0.5*ABS(XR - XL)
      X = 0.5*(XL + XR)
      R1 = SGN*F(X,IZ,C,A,BH)
      IF (R1 >= 0.) THEN
         IF (R1 == 0.) GO TO 105
         XR = X
      ELSE
         XL = X
      ENDIF
      DX = 0.5*DX
      DO WHILE(DX - CNV > 0.)
         X = 0.5*(XL + XR)
         R1 = SGN*F(X,IZ,C,A,BH)
         IF (R1 >= 0.) THEN
            IF (R1 == 0.) GO TO 105
            XR = X
         ELSE
            XL = X
         ENDIF
         DX = 0.5*DX
      END DO
  105 CONTINUE
      CBSRH = 0.5*(XL + XR)
      RETURN 
      END FUNCTION CBSRH


      SUBROUTINE CCOMPB(N, IERROR, AN, BN, CN, B, BC, AH, BH)
      IMPLICIT NONE
      Real epmach
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER  :: N
      INTEGER  :: IERROR
      REAL  :: AN(*)
      REAL , INTENT(IN) :: BN(*)
      REAL  :: CN(*)
      REAL  :: B(*)
      REAL  :: AH(*)
      REAL  :: BH(*)
      COMPLEX  :: BC(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, IF, KDO, L, IR, I2, I4, IPL, IFD, I, IB, NB, JS, JF
     1   , LS, LH, NMP, L1, L2, J2, J1, N2M2
      REAL :: DUM, BNORM, ARG, D1, D2, D3
C-----------------------------------------------
C
C     CCOMPB COMPUTES THE ROOTS OF THE B POLYNOMIALS USING SUBROUTINE
C     CTEVLS WHICH IS A MODIFICATION THE EISPACK PROGRAM TQLRAT.
C     IERROR IS SET TO 4 IF EITHER CTEVLS FAILS OR IF A(J+1)*C(J) IS
C     LESS THAN ZERO FOR SOME J.  AH,BH ARE TEMPORARY WORK ARRAYS.
C
 
      EPS = EPMACH(DUM)
      BNORM = ABS(BN(1))
      DO J = 2, NM
         BNORM = AMAX1(BNORM,ABS(BN(J)))
         ARG = AN(J)*CN(J-1)
         IF (ARG < 0.) GO TO 119
         B(J) = SIGN(SQRT(ARG),AN(J))
      END DO
      CNV = EPS*BNORM
      IF = 2**K
      KDO = K - 1
      L108: DO L = 1, KDO
         IR = L - 1
         I2 = 2**IR
         I4 = I2 + I2
         IPL = I4 - 1
         IFD = IF - I4
         DO I = I4, IFD, I4
            CALL CINDXB (I, L, IB, NB)
            IF (NB <= 0) CYCLE  L108
            JS = I - IPL
            JF = JS + NB - 1
            LS = 0
            BH(:JF-JS+1) = BN(JS:JF)
            AH(:JF-JS+1) = B(JS:JF)
            CALL CTEVLS (NB, BH, AH, IERROR)
            IF (IERROR /= 0) GO TO 118
            LH = IB - 1
            IF (NB > 0) THEN
               B(LH+1:NB+LH) = -BH(:NB)
               LH = NB + LH
            ENDIF
         END DO
      END DO L108
      B(:NM) = -BN(:NM)
      IF (NPP == 0) THEN
         NMP = NM + 1
         NB = NM + NMP
         DO J = 1, NB
            L1 = MOD(J - 1,NMP) + 1
            L2 = MOD(J + NM - 1,NMP) + 1
            ARG = AN(L1)*CN(L2)
            IF (ARG < 0.) GO TO 119
            BH(J) = SIGN(SQRT(ARG),(-AN(L1)))
            AH(J) = -BN(L1)
         END DO
         CALL CTEVLS (NB, AH, BH, IERROR)
         IF (IERROR /= 0) GO TO 118
         CALL CINDXB (IF, K - 1, J2, LH)
         CALL CINDXB (IF/2, K - 1, J1, LH)
         J2 = J2 + 1
         LH = J2
         N2M2 = J2 + NM + NM - 2
  114    CONTINUE
         D1 = ABS(B(J1)-B(J2-1))
         D2 = ABS(B(J1)-B(J2))
         D3 = ABS(B(J1)-B(J2+1))
         IF (D2>=D1 .OR. D2>=D3) THEN
            B(LH) = B(J2)
            J2 = J2 + 1
            LH = LH + 1
            IF (J2 - N2M2 <= 0) GO TO 114
         ELSE
            J2 = J2 + 1
            J1 = J1 + 1
            IF (J2 - N2M2 <= 0) GO TO 114
         ENDIF
         B(LH) = B(N2M2+1)
         CALL CINDXB (IF, K - 1, J1, J2)
         J2 = J1 + NMP + NMP
         CALL CPPADD (NM + 1, IERROR, AN, CN, B(J1), BC(J1), B(J2))
      ENDIF
      RETURN 
  118 CONTINUE
      IERROR = 4
      RETURN 
  119 CONTINUE
      IERROR = 5
      RETURN 
      END SUBROUTINE CCOMPB


      SUBROUTINE CPROC(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,W,YY)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      COMPLEX , INTENT(IN) :: BD(*)
      COMPLEX , INTENT(IN) :: X(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: W(*)
      COMPLEX  :: YY(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, ID, M1, M2, IA, IFLG, K
      REAL :: RT
      COMPLEX :: CRT, DEN, Y1, Y2
C-----------------------------------------------
C
C PROC APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,W ARE WORK ARRAYS
C ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      Y(:M) = X(:M)
      MM = M - 1
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IFLG = 0
      IF (ID > 0) THEN
         CRT = BD(ID)
         ID = ID - 1
C
C BEGIN SOLUTION TO SYSTEM
C
         D(M) = A(M)/(B(M)-CRT)
         W(M) = Y(M)/(B(M)-CRT)
         DO J = 2, MM
            K = M - J
            DEN = B(K+1) - CRT - C(K+1)*D(K+2)
            D(K+1) = A(K+1)/DEN
            W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
         END DO
         DEN = B(1) - CRT - C(1)*D(2)
         IF (CABS(DEN) /= 0.) THEN
            Y(1) = (Y(1)-C(1)*W(2))/DEN
         ELSE
            Y(1) = (1.,0.)
         ENDIF
         DO J = 2, M
            Y(J) = W(J) - D(J)*Y(J-1)
         END DO
      ENDIF
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 121
         RT = BM2(M2)
         M2 = M2 - 1
      ELSE
         IF (M2 <= 0) THEN
            RT = BM1(M1)
            M1 = M1 - 1
         ELSE
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) > 0.) THEN
               RT = BM1(M1)
               M1 = M1 - 1
            ELSE
               RT = BM2(M2)
               M2 = M2 - 1
            ENDIF
         ENDIF
      ENDIF
      Y1 = (B(1)-RT)*Y(1) + C(1)*Y(2)
      IF (MM - 2 >= 0) THEN
         DO J = 2, MM
            Y2 = A(J)*Y(J-1) + (B(J)-RT)*Y(J) + C(J)*Y(J+1)
            Y(J-1) = Y1
            Y1 = Y2
         END DO
      ENDIF
      Y(M) = A(M)*Y(M-1) + (B(M)-RT)*Y(M)
      Y(M-1) = Y1
      IFLG = 1
      GO TO 102
  121 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IA = IA - 1
         IFLG = 1
C
C SCALAR MULTIPLICATION
C
         Y(:M) = RT*Y(:M)
      ENDIF
      IF (IFLG > 0) GO TO 102
      RETURN 
      END SUBROUTINE CPROC


      SUBROUTINE CPROCP(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,U,YY)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      REAL  :: YY(*)
      COMPLEX , INTENT(IN) :: BD(*)
      COMPLEX , INTENT(IN) :: X(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: U(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, MM2, ID, M1, M2, IA, IFLG, K
      REAL :: RT
      COMPLEX :: V, DEN, BH, YM, AM, Y1, Y2, YH, CRT
C-----------------------------------------------
C
C CPROCP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y
C
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,U ARE WORK ARRAYS
C ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      Y(:M) = X(:M)
      MM = M - 1
      MM2 = M - 2
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IFLG = 0
      IF (ID > 0) THEN
         CRT = BD(ID)
         ID = ID - 1
         IFLG = 1
C
C BEGIN SOLUTION TO SYSTEM
C
         BH = B(M) - CRT
         YM = Y(M)
         DEN = B(1) - CRT
         D(1) = C(1)/DEN
         U(1) = A(1)/DEN
         Y(1) = Y(1)/DEN
         V = C(M)
         IF (MM2 - 2 >= 0) THEN
            DO J = 2, MM2
               DEN = B(J) - CRT - A(J)*D(J-1)
               D(J) = C(J)/DEN
               U(J) = -A(J)*U(J-1)/DEN
               Y(J) = (Y(J)-A(J)*Y(J-1))/DEN
               BH = BH - V*U(J-1)
               YM = YM - V*Y(J-1)
               V = -V*D(J-1)
            END DO
         ENDIF
         DEN = B(M-1) - CRT - A(M-1)*D(M-2)
         D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
         Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2))/DEN
         AM = A(M) - V*D(M-2)
         BH = BH - V*U(M-2)
         YM = YM - V*Y(M-2)
         DEN = BH - AM*D(M-1)
         IF (CABS(DEN) /= 0.) THEN
            Y(M) = (YM - AM*Y(M-1))/DEN
         ELSE
            Y(M) = (1.,0.)
         ENDIF
         Y(M-1) = Y(M-1) - D(M-1)*Y(M)
         DO J = 2, MM
            K = M - J
            Y(K) = Y(K) - D(K)*Y(K+1) - U(K)*Y(M)
         END DO
      ENDIF
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 123
         RT = BM2(M2)
         M2 = M2 - 1
      ELSE
         IF (M2 <= 0) THEN
            RT = BM1(M1)
            M1 = M1 - 1
         ELSE
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) > 0.) THEN
               RT = BM1(M1)
               M1 = M1 - 1
            ELSE
               RT = BM2(M2)
               M2 = M2 - 1
C
C MATRIX MULTIPLICATION
C
            ENDIF
         ENDIF
      ENDIF
      YH = Y(1)
      Y1 = (B(1)-RT)*Y(1) + C(1)*Y(2) + A(1)*Y(M)
      IF (MM - 2 >= 0) THEN
         DO J = 2, MM
            Y2 = A(J)*Y(J-1) + (B(J)-RT)*Y(J) + C(J)*Y(J+1)
            Y(J-1) = Y1
            Y1 = Y2
         END DO
      ENDIF
      Y(M) = A(M)*Y(M-1) + (B(M)-RT)*Y(M) + C(M)*YH
      Y(M-1) = Y1
      IFLG = 1
      GO TO 102
  123 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IA = IA - 1
         IFLG = 1
C
C SCALAR MULTIPLICATION
C
         Y(:M) = RT*Y(:M)
      ENDIF
      IF (IFLG > 0) GO TO 102
      RETURN 
      END SUBROUTINE CPROCP


      SUBROUTINE CINDXA(I, IR, IDXA, NA)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDXA
      INTEGER , INTENT(OUT) :: NA
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
      NA = 2**IR
      IDXA = I - NA + 1
      IF (I - NM > 0) THEN
         NA = 0
      ENDIF
      RETURN 
      END SUBROUTINE CINDXA


      SUBROUTINE CINDXB(I, IR, IDX, IDP)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDX
      INTEGER , INTENT(OUT) :: IDP
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IZH, ID, IPL
C-----------------------------------------------
C
C B(IDX) IS THE LOCATION OF THE FIRST ROOT OF THE B(I,IR) POLYNOMIAL
C
      IDP = 0
      IF (IR >= 0) THEN
         IF (IR <= 0) THEN
            IF (I - NM > 0) GO TO 107
            IDX = I
            IDP = 1
            RETURN 
         ENDIF
         IZH = 2**IR
         ID = I - IZH - IZH
         IDX = ID + ID + (IR - 1)*IK + IR + (IK - I)/IZH + 4
         IPL = IZH - 1
         IDP = IZH + IZH - 1
         IF (I - IPL - NM > 0) THEN
            IDP = 0
            RETURN 
         ENDIF
         IF (I + IPL - NM > 0) THEN
            IDP = NM + IPL - I + 1
         ENDIF
      ENDIF
  107 CONTINUE
      RETURN 
      END SUBROUTINE CINDXB


      SUBROUTINE CINDXC(I, IR, IDXC, NC)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDXC
      INTEGER , INTENT(OUT) :: NC
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
C-----------------------------------------------
      NC = 2**IR
      IDXC = I
      IF (IDXC + NC - 1 - NM > 0) THEN
         NC = 0
      ENDIF
      RETURN 
      END SUBROUTINE CINDXC


      SUBROUTINE CPPADD(N, IERROR, A, C, CBP, BP, BH)
      IMPLICIT NONE
      real psgf,ppspf,ppsgf,cbsrh
      EXTERNAL        PSGF       ,PPSPF      ,PPSGF                     
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(OUT) :: IERROR
      REAL  :: A(*)
      REAL  :: C(*)
      REAL , INTENT(INOUT) :: BP(*)
      REAL  :: BH(*)
      COMPLEX , INTENT(INOUT) :: CBP(*)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /CCBLK/ NPP, K, EPS, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   EPS, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::IZ,IZM,IZM2,J,NT,MODIZ,IS,IF,IG,IT,ICV,I3,I2,NHALF
      REAL :: R4, R5, R6, SCNV, XL, DB, SGN, XR, XM, PSG
      COMPLEX :: CF, CX, FSG, HSG, DD, F, FP, FPP, CDIS, R1, R2, R3
C-----------------------------------------------
C
C     CPPADD COMPUTES THE EIGENVALUES OF THE PERIODIC TRIDIAGONAL
C     MATRIX WITH COEFFICIENTS AN,BN,CN
C
C     N IS THE ORDER OF THE BH AND BP POLYNOMIALS
C     ON OUTPUT BP CONTAINS THE EIGENVALUES
C     CBP IS THE SAME AS BP EXCEPT TYPE COMPLEX
C     BH IS USED TO TEMPORARILY STORE THE ROOTS OF THE B HAT POLYNOMIAL
C       WHICH ENTERS THROUGH BP
C
      SCNV = SQRT(CNV)
      IZ = N
      IZM = IZ - 1
      IZM2 = IZ - 2
      IF (BP(N) - BP(1) <= 0.) THEN
         IF (BP(N) - BP(1) == 0.) GO TO 142
         BH(:N) = BP(N:1:(-1))
      ELSE
         BH(:N) = BP(:N)
      ENDIF
      NCMPLX = 0
      MODIZ = MOD(IZ,2)
      IS = 1
      IF (MODIZ /= 0) THEN
         IF (A(1) < 0.) GO TO 110
         IF (A(1) == 0.) GO TO 142
      ENDIF
      XL = BH(1)
      DB = BH(3) - BH(1)
      XL = XL - DB
      R4 = PSGF(XL,IZ,C,A,BH)
      DO WHILE(R4 <= 0.)
         XL = XL - DB
         R4 = PSGF(XL,IZ,C,A,BH)
      END DO
      SGN = -1.
      CBP(1) = CMPLX(CBSRH(XL,BH(1),IZ,C,A,BH,PSGF,SGN),0.)
      IS = 2
  110 CONTINUE
      IF = IZ - 1
      IF (MODIZ /= 0) THEN
         IF (A(1) > 0.) GO TO 115
         IF (A(1) == 0.) GO TO 142
      ENDIF
      XR = BH(IZ)
      DB = BH(IZ) - BH(IZ-2)
      XR = XR + DB
      R5 = PSGF(XR,IZ,C,A,BH)
      DO WHILE(R5 < 0.)
         XR = XR + DB
         R5 = PSGF(XR,IZ,C,A,BH)
      END DO
      SGN = 1.
      CBP(IZ) = CMPLX(CBSRH(BH(IZ),XR,IZ,C,A,BH,PSGF,SGN),0.)
      IF = IZ - 2
  115 CONTINUE
      DO IG = IS, IF, 2
         XL = BH(IG)
         XR = BH(IG+1)
         SGN = -1.
         XM = CBSRH(XL,XR,IZ,C,A,BH,PPSPF,SGN)
         PSG = PSGF(XM,IZ,C,A,BH)
         IF (ABS(PSG) - EPS <= 0.) GO TO 118
         R6 = PSG*PPSGF(XM,IZ,C,A,BH)
         IF (R6 > 0.) GO TO 119
         IF (R6 == 0.) GO TO 118
         SGN = 1.
         CBP(IG) = CMPLX(CBSRH(BH(IG),XM,IZ,C,A,BH,PSGF,SGN),0.)
         SGN = -1.
         CBP(IG+1) = CMPLX(CBSRH(XM,BH(IG+1),IZ,C,A,BH,PSGF,SGN),0.)
         CYCLE 
C
C     CASE OF A MULTIPLE ZERO
C
  118    CONTINUE
         CBP(IG) = CMPLX(XM,0.)
         CBP(IG+1) = CMPLX(XM,0.)
         CYCLE 
C
C     CASE OF A COMPLEX ZERO
C
  119    CONTINUE
         IT = 0
         ICV = 0
         CX = CMPLX(XM,0.)
  120    CONTINUE
         FSG = (1.,0.)
         HSG = (1.,0.)
         FP = (0.,0.)
         FPP = (0.,0.)
         DO J = 1, IZ
            DD = 1./(CX - BH(J))
            FSG = FSG*A(J)*DD
            HSG = HSG*C(J)*DD
            FP = FP + DD
            FPP = FPP - DD*DD
         END DO
         IF (MODIZ == 0) THEN
            F = (1.,0.) - FSG - HSG
         ELSE
            F = (1.,0.) + FSG + HSG
         ENDIF
         I3 = 0
         IF (CABS(FP) > 0.) THEN
            I3 = 1
            R3 = -F/FP
         ENDIF
         I2 = 0
         IF (CABS(FPP) > 0.) THEN
            I2 = 1
            CDIS = CSQRT(FP**2 - 2.*F*FPP)
            R1 = CDIS - FP
            R2 = (-FP) - CDIS
            IF (CABS(R1) - CABS(R2) > 0.) THEN
               R1 = R1/FPP
            ELSE
               R1 = R2/FPP
            ENDIF
            R2 = 2.*F/FPP/R1
            IF (CABS(R2) < CABS(R1)) R1 = R2
            IF (I3 <= 0) GO TO 133
            IF (CABS(R3) < CABS(R1)) R1 = R3
            GO TO 133
         ENDIF
         R1 = R3
  133    CONTINUE
         CX = CX + R1
         IT = IT + 1
         IF (IT > 50) GO TO 142
         IF (CABS(R1) > SCNV) GO TO 120
         IF (ICV > 0) GO TO 135
         ICV = 1
         GO TO 120
  135    CONTINUE
         CBP(IG) = CX
         CBP(IG+1) = CONJG(CX)
      END DO
      IF (CABS(CBP(N)) - CABS(CBP(1)) <= 0.) THEN
         IF (CABS(CBP(N)) - CABS(CBP(1)) == 0.) GO TO 142
         NHALF = N/2
         DO J = 1, NHALF
            NT = N - J
            CX = CBP(J)
            CBP(J) = CBP(NT+1)
            CBP(NT+1) = CX
         END DO
      ENDIF
      NCMPLX = 1
      DO J = 2, IZ
         IF (AIMAG(CBP(J)) /= 0.) GO TO 143
      END DO
      NCMPLX = 0
      DO J = 2, IZ
         BP(J) = REAL(CBP(J))
      END DO
      GO TO 143
  142 CONTINUE
      IERROR = 4
  143 CONTINUE
      RETURN 
      END SUBROUTINE CPPADD


      SUBROUTINE PROC(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,W,U)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BD(*)
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      COMPLEX , INTENT(IN) :: X(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: W(*)
      COMPLEX  :: U(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, ID, IBR, M1, M2, IA, K
      REAL :: RT
      COMPLEX :: DEN
C-----------------------------------------------
C
C PROC APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,W,U ARE WORKING ARRAYS
C IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      W(:M) = X(:M)
      Y(:M) = W(:M)
      MM = M - 1
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IF (ND == 0) RT = -RT
         IA = IA - 1
C
C SCALAR MULTIPLICATION
C
         Y(:M) = RT*W(:M)
      ENDIF
      IF (ID <= 0) GO TO 125
      RT = BD(ID)
      ID = ID - 1
      IF (ID == 0) IBR = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      D(M) = A(M)/(B(M)-RT)
      W(M) = Y(M)/(B(M)-RT)
      DO J = 2, MM
         K = M - J
         DEN = B(K+1) - RT - C(K+1)*D(K+2)
         D(K+1) = A(K+1)/DEN
         W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
      END DO
      DEN = B(1) - RT - C(1)*D(2)
      W(1) = (1.,0.)
      IF (CABS(DEN) /= 0.) THEN
         W(1) = (Y(1)-C(1)*W(2))/DEN
      ENDIF
      DO J = 2, M
         W(J) = W(J) - D(J)*W(J-1)
      END DO
      IF (NA > 0) GO TO 102
      GO TO 113
  111 CONTINUE
      Y(:M) = W(:M)
      IBR = 1
      GO TO 102
  113 CONTINUE
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 111
      ELSE
         IF (M2 > 0) THEN
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) <= 0.) GO TO 120
         ENDIF
         IF (IBR <= 0) THEN
            IF (ABS(BM1(M1)-BD(ID)) - ABS(BM1(M1)-RT) < 0.) GO TO 111
         ENDIF
         RT = RT - BM1(M1)
         M1 = M1 - 1
         GO TO 123
      ENDIF
  120 CONTINUE
      IF (IBR <= 0) THEN
         IF (ABS(BM2(M2)-BD(ID)) - ABS(BM2(M2)-RT) < 0.) GO TO 111
      ENDIF
      RT = RT - BM2(M2)
      M2 = M2 - 1
  123 CONTINUE
      Y(:M) = Y(:M) + RT*W(:M)
      GO TO 102
  125 CONTINUE
      RETURN 
      END SUBROUTINE PROC


      SUBROUTINE PROCP(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,U,W)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL , INTENT(IN) :: BD(*)
      REAL , INTENT(IN) :: BM1(*)
      REAL , INTENT(IN) :: BM2(*)
      REAL , INTENT(IN) :: AA(*)
      COMPLEX , INTENT(IN) :: X(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
      COMPLEX , INTENT(IN) :: A(*)
      COMPLEX , INTENT(IN) :: B(*)
      COMPLEX , INTENT(IN) :: C(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: U(*)
      COMPLEX , INTENT(INOUT) :: W(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: J, MM, MM2, ID, IBR, M1, M2, IA, K
      REAL :: RT
      COMPLEX :: DEN, YM, V, BH, AM
C-----------------------------------------------
C
C PROCP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
C STORES THE RESULT IN Y        PERIODIC BOUNDARY CONDITIONS
C
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
C NA IS THE LENGTH OF THE ARRAY AA
C X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
C M  IS THE ORDER OF THE MATRIX
C D,U,W ARE WORKING ARRAYS
C IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
C
      Y(:M) = X(:M)
      W(:M) = Y(:M)
      MM = M - 1
      MM2 = M - 2
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IF (ND == 0) RT = -RT
         IA = IA - 1
         Y(:M) = RT*W(:M)
      ENDIF
      IF (ID <= 0) GO TO 128
      RT = BD(ID)
      ID = ID - 1
      IF (ID == 0) IBR = 1
C
C BEGIN SOLUTION TO SYSTEM
C
      BH = B(M) - RT
      YM = Y(M)
      DEN = B(1) - RT
      D(1) = C(1)/DEN
      U(1) = A(1)/DEN
      W(1) = Y(1)/DEN
      V = C(M)
      IF (MM2 - 2 >= 0) THEN
         DO J = 2, MM2
            DEN = B(J) - RT - A(J)*D(J-1)
            D(J) = C(J)/DEN
            U(J) = -A(J)*U(J-1)/DEN
            W(J) = (Y(J)-A(J)*W(J-1))/DEN
            BH = BH - V*U(J-1)
            YM = YM - V*W(J-1)
            V = -V*D(J-1)
         END DO
      ENDIF
      DEN = B(M-1) - RT - A(M-1)*D(M-2)
      D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
      W(M-1) = (Y(M-1)-A(M-1)*W(M-2))/DEN
      AM = A(M) - V*D(M-2)
      BH = BH - V*U(M-2)
      YM = YM - V*W(M-2)
      DEN = BH - AM*D(M-1)
      IF (CABS(DEN) /= 0.) THEN
         W(M) = (YM - AM*W(M-1))/DEN
      ELSE
         W(M) = (1.,0.)
      ENDIF
      W(M-1) = W(M-1) - D(M-1)*W(M)
      DO J = 2, MM
         K = M - J
         W(K) = W(K) - D(K)*W(K+1) - U(K)*W(M)
      END DO
      IF (NA > 0) GO TO 102
      GO TO 116
  114 CONTINUE
      Y(:M) = W(:M)
      IBR = 1
      GO TO 102
  116 CONTINUE
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 114
      ELSE
         IF (M2 > 0) THEN
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) <= 0.) GO TO 123
         ENDIF
         IF (IBR <= 0) THEN
            IF (ABS(BM1(M1)-BD(ID)) - ABS(BM1(M1)-RT) < 0.) GO TO 114
         ENDIF
         RT = RT - BM1(M1)
         M1 = M1 - 1
         GO TO 126
      ENDIF
  123 CONTINUE
      IF (IBR <= 0) THEN
         IF (ABS(BM2(M2)-BD(ID)) - ABS(BM2(M2)-RT) < 0.) GO TO 114
      ENDIF
      RT = RT - BM2(M2)
      M2 = M2 - 1
  126 CONTINUE
      Y(:M) = Y(:M) + RT*W(:M)
      GO TO 102
  128 CONTINUE
      RETURN 
      END SUBROUTINE PROCP


      SUBROUTINE CTEVLS(N, D, E2, IERR)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(OUT) :: IERR
      REAL , INTENT(INOUT) :: D(N)
      REAL , INTENT(INOUT) :: E2(N)
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /CCBLK/
      COMMON /CCBLK/ NPP, K, MACHEP, CNV, NM, NCMPLX, IK
      INTEGER   NPP, K, NM, NCMPLX, IK
      REAL   MACHEP, CNV
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, J, L, M, II, L1, MML, NHALF, NTOP
      REAL :: B, C, F, G, H, P, R, S, DHOLD
C-----------------------------------------------
C
C
C     REAL SQRT,ABS,SIGN
C
C
C     THIS SUBROUTINE IS A MODIFICATION OF THE EISPACK SUBROUTINE TQLRAT
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E2 CONTAINS THE                SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES,
C
C        E2 HAS BEEN DESTROYED,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
C
      IERR = 0
      IF (N /= 1) THEN
C
         E2(:N-1) = E2(2:N)*E2(2:N)
C
         F = 0.0
         B = 0.0
         E2(N) = 0.0
C
         DO L = 1, N
            J = 0
            H = MACHEP*(ABS(D(L))+SQRT(E2(L)))
            IF (B <= H) THEN
               B = H
               C = B*B
            ENDIF
C
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
C
            DO M = L, N
               IF (E2(M) > C) CYCLE 
               EXIT 
C
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
C
            END DO
C
            IF (M /= L) THEN
  105          CONTINUE
               IF (J == 30) GO TO 114
               J = J + 1
C
C     ********** FORM SHIFT **********
C
               L1 = L + 1
               S = SQRT(E2(L))
               G = D(L)
               P = (D(L1)-G)/(2.0*S)
               R = SQRT(P*P + 1.0)
               D(L) = S/(P + SIGN(R,P))
               H = G - D(L)
C
               D(L1:N) = D(L1:N) - H
C
               F = F + H
C
C     ********** RATIONAL QL TRANSFORMATION **********
C
               G = D(M)
               IF (G == 0.0) G = B
               H = G
               S = 0.0
               MML = M - L
C
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
C
               DO II = 1, MML
                  I = M - II
                  P = G*H
                  R = P + E2(I)
                  E2(I+1) = S*R
                  S = E2(I)/R
                  D(I+1) = H + S*(H + D(I))
                  G = D(I) - E2(I)/G
                  IF (G == 0.0) G = B
                  H = G*P/R
               END DO
C
               E2(L) = S*G
               D(L) = H
C
C     ********** GUARD AGAINST UNDERFLOWED H **********
C
               IF (H == 0.0) GO TO 108
               IF (ABS(E2(L)) <= ABS(C/H)) GO TO 108
               E2(L) = H*E2(L)
               IF (E2(L) /= 0.0) GO TO 105
            ENDIF
  108       CONTINUE
            P = D(L) + F
C
C     ********** ORDER EIGENVALUES **********
C
            IF (L /= 1) THEN
C
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
C
               DO II = 2, L
                  I = L + 2 - II
                  IF (P >= D(I-1)) GO TO 111
                  D(I) = D(I-1)
               END DO
            ENDIF
C
            I = 1
  111       CONTINUE
            D(I) = P
         END DO
C
         IF (ABS(D(N)) >= ABS(D(1))) GO TO 115
         NHALF = N/2
         DO I = 1, NHALF
            NTOP = N - I
            DHOLD = D(I)
            D(I) = D(NTOP+1)
            D(NTOP+1) = DHOLD
         END DO
         GO TO 115
C
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
C
  114    CONTINUE
         IERR = L
      ENDIF
  115 CONTINUE
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
      END SUBROUTINE CTEVLS
