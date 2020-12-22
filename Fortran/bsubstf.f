cxx      SUBROUTINE BSUBSTF( ZS,N,IMODEL,LAG,K,IL,LG1,LG2,F,CNST,ZMEAN,SUM,
      SUBROUTINE BSUBSTF( ZS,N,IMODEL,LAG,K,IL,LG1,LG2,ZMEAN,SUM,
     *M,AICM,SDM,A1,SD,AIC,DIC,AICB,SDB,EK,A2,IND,C,C1,C2,B,OEIC,ESUM,
     *OMEAN,OM,E,EMEAN,VARI,SKEW,PEAK,COV,SXX )
C
    
C
cc      PROGRAM BSUBST
C.......................................................................
C.....PLANNED BY H.AKAIKE...............................................
C.....DESIGNED BY H.AKAIKE AND G.KITAGAWA...............................
C.....PROGRAMMED BY G.KITAGAWA AND F.TADA...............................
C.....ADDRESS: THE INSTITUTE OF STATISTICAL MATHEMATICS, 4-6-7 MINAMI-AZ
C..............MINATO-KU, TOKYO 106, JAPAN..............................
C.....DATE OF THE LATEST REVISION:  MAY  14, 1979.......................
C.......................................................................
C.....THIS PROGRAM WAS ORIGINALLY PUBLISHED IN "TIMSAC-78", BY H.AKAIKE,
C.....G.KITAGAWA, E.ARAHATA AND F.TADA, COMPUTER SCIENCE MONOGRAPHS, NO.
C.....THE INSTITUTE OF STATISTICAL MATHEMATICS, TOKYO, 1979.............
C.......................................................................
C     TIMSAC 78.1.3.                                                    
C     _                 ____ _                                          
C     BAYESIAN TYPE ALL SUBSET ANALYSIS OF TIME SERIES BY A MODEL LINEAR
C     IN PARAMETERS                                                     
C                                                                       
C     THIS PROGRAM PRODUCES BAYESIAN ESTIMATES OF TIME SERIES MODELS SUC
C     PURE AR MODELS, AR-MODELS WITH NON-LINEAR TERMS, AR-MODELS WITH PO
C     TYPE MEAN VALUE FUNCTIONS, ETC.  THE GOODNESS OF FIT OF A MODEL IS
C     CHECKED BY THE ANALYSIS OF SEVERAL STEPS AHEAD PREDICTION ERRORS. 
C     BY PREPARING AN EXTERNAL SUBROUTINE SETX PROPERLY, ANY TIME SERIES
C     WHICH IS LINEAR IN PARAMETERS CAN BE TREATED.                     
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
C             REDATA                                                    
C             REDLAG                                                    
C             SETLAG                                                    
C             REDREG                                                    
C             REDUCT                                                    
C             ARMFIT                                                    
C             SBBAYS                                                    
C             CHECK                                                     
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS REQUIRED:                                                
C          MT:    ORIGINAL DATA INPUT DEVICE SPECIFICATION              
C          IMODEL:=1  AUTOREGRESSIVE MODEL                              
C                 =2  POLYNOMIAL TYPE NON-LINEAR MODEL (LAG'S READ IN ) 
C                 =3  POLYNOMIAL TYPE NON-LINEAR MODEL (LAG'S AUTOMATICA
C                 =4  AR-MODEL WITH POLYNOMIAL MEAN VALUE FUNCTION      
C                 =5  ANY NON-LINEAR MODEL                              
C                 =6  POLYNOMIAL TYPE EXPONENTIALLY DAMPED NON-LINEAR MO
C                 =7  THIS MODEL IS RESERVED FOR THE USER'S OPTIONAL USE
C          LAG:   MAXIMUM TIME LAG USED IN THE MODEL                    
C          K:     NUMBER OF REGRESSORS                                  
C          IL:    PREDICTION ERRORS CHECKING (UP TO IL-STEPS AHEAD) IS R
C                 N*IL SHOULD BE LESS THAN OR EQUAL TO 20000            
C                                                                       
C       --   THE FOLLOWING INPUTS ARE REQUIRED AT SUBROUTINE REDATA   --
C                                                                       
C          TITLE:   ORIGINAL DATA SPECIFICATION                         
C          N:       DATA LENGTH                                         
C          DFORM:   INPUT DATA FORMAT SPECIFICATION STATEMENT           
C                   -- EXAMPLE --  (8F10.5 )                            
C          X(I) (I=1,N):   ORIGINAL DATA                                
C       ----------------------------------------------------------------
C                                                                       
cc      !DEC$ ATTRIBUTES DLLEXPORT :: BSUBSTF
C
      PARAMETER  ( MJ2 = 101 )
C
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL*4     Z(10000) , TITLE(20) , TTL(5) , E(20000)               
cc      REAL*4     TITLE(20) , TTL(5)
cc      DIMENSION  Z(10000), E(20000) 
cc      DIMENSION  X(200,101), D(200) , A(100) , B(100)                  
cc      DATA  TTL  / 4H   B,4HAYES,4HIAN ,4HMODE,4HL    /                 
cxx      COMMON     / BBB / L1(50) , L2(50) , L3(50) , SD0 , CONST         
cc      REAL * 4   CONST , SD0                                            
cxx      DIMENSION  ZS(N), Z(N), E(N,IL) 
cxx      DIMENSION  X(N,K+1), A1(K), A2(K), B(K)                  
cxx      DIMENSION  LG1(3,K), LG2(5)
cxx      DIMENSION  SD(K+1), AIC(K+1), DIC(K+1)
cxx      DIMENSION  IND(K), C(K), C1(K+1), C2(K), ESUM(K+1)
cxx      DIMENSION  EMEAN(IL), VARI(IL), SKEW(IL), PEAK(IL), COV(MJ2)
cxx      DIMENSION  SXX(121), FA(N-LAG,IL)
cc      CHARACTER*4 F(20,K+1)
cxx      INTEGER*1  F(80,K+1)
cxx      CHARACTER(4) G
      INTEGER :: N, IMODEL, LAG, K, IL, LG1(3,K), LG2(5), M, IND(K) 
      REAL(8) :: ZS(N), ZMEAN, SUM, AICM, SDM, A1(K), SD(K+1), AIC(K+1),
     1           DIC(K+1), AICB, SDB, EK, A2(K), C(K), C1(K+1), C2(K),
     2           B(K), OEIC, ESUM(K+1), OMEAN, OM, E(N,IL),
     3           EMEAN(IL), VARI(IL), SKEW(IL), PEAK(IL), COV(MJ2),
     4           SXX(121)
      REAL(8) :: Z(N), X(N,K+1), FA(N-LAG,IL), SD0, CONST
      CHARACTER(4) G
C
      COMMON     / BBB / L1(50) , L2(50) , L3(50) , SD0 , CONST         
      COMMON     / EEE /  G(20,31)                                      
      COMMON     / AAA /  NN
C                                                                       
C        EXTERNAL SUBROUTINE DECLARATION:                               
C                                                                       
      EXTERNAL  SETX1                                                   
      EXTERNAL  SETX2                                                   
      EXTERNAL  SETX4                                                   
cx      EXTERNAL  SETX5                                                   
cx      EXTERNAL  SETX6                                                   
cc      EXTERNAL  SETX7                                                   
      EXTERNAL  PRDCT1                                                  
      EXTERNAL  PRDCT2                                                  
cx      EXTERNAL  PRDCT3                                                  
cx      EXTERNAL  PRDCT6                                                  
C
cc	CHARACTER(100) IFLNAM,OFLNAM
cc	CALL FLNAM2( IFLNAM,OFLNAM,NFL )
cc	IF ( NFL.EQ.0 ) GO TO 999
cc	IF ( NFL.EQ.2 ) THEN
cc	   OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR)
cc	ELSE
cc	   CALL SETWND
cc	END IF
C
C
C        PARAMETERS:                                                    
C             MJ1:  ABSOLUTE DIMENSION FOR SUBROUTINE CALL              
C                                                                       
cc      IF ((IMODEL.LE.0) .OR. (IMODEL.GE.8)) GO TO 150
      IF ((IMODEL.LE.0) .OR. (IMODEL.GE.7)) GO TO 150
cc      MJ = 1000                                                         
cc      MJ1 = 200                                                         
      NN = N
      MJ = N
      MJ1 = N
      ISW = 1                                                           
      IPR = 2                                                           
C                                                                       
CC      READ( 5,1 )     MT                                                
cc      MT = 5
cc      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
cc      READ( 5,1 )     IMODEL , LAG , K , IL                             
cc      WRITE( 6,3 )                                                      
cc      IF( IMODEL.EQ.1 )  WRITE( 6,4 )   IMODEL                          
cc      IF( IMODEL.EQ.2 )  WRITE( 6,5 )   IMODEL                          
cc      IF( IMODEL.EQ.3 )  WRITE( 6,5 )   IMODEL                          
cc      IF( IMODEL.EQ.5 )  WRITE( 6,5 )   IMODEL                          
cc      IF( IMODEL.EQ.6 )  WRITE( 6,11 )   IMODEL                         
cc      WRITE( 6,6 )                                                      
cc      IF( IMODEL.EQ.1 )  WRITE( 6,7 )  K                                
cc      IF( IMODEL.EQ.2 )  WRITE( 6,8 )   LAG , K                         
cc      WRITE( 6,2 )     MT                                               
C                                                                       
C                                                                       
C          ---------------------------------------                      
C          ORIGINAL DATA LOADING AND MEAN DELETION                      
C          ---------------------------------------                      
C                                                                       
cc      CALL  REDATA( Z,N,MT,TITLE )                                      
      CALL  REDATA( ZS,Z,N,ZMEAN,SUM )
      NMK = N - LAG                                                     
      LAG1 = LAG + 1                                                    
C                                                                       
C          ---------------------                                        
C          HOUSEHOLDER REDUCTION                                        
C          ---------------------                                        
C                                                                       
cc      GO TO ( 10,20,30,40,50,60,70 ), IMODEL                            
cx      GO TO ( 10,20,30,40,50,60 ), IMODEL                            
cxx      GO TO ( 10,20,30,40 ), IMODEL                            
      IF ( IMODEL .EQ. 1 ) GO TO 10
      IF ( IMODEL .EQ. 2 ) GO TO 20
      IF ( IMODEL .EQ. 3 ) GO TO 30
      IF ( IMODEL .EQ. 4 ) GO TO 40
C                                                                       
   10 K = LAG                                                           
cc      CALL  REDUCT( SETX1,Z,D,NMK,0,K,MJ1,LAG,X )                       
      CALL  REDUCT( SETX1,Z,NMK,0,K,MJ1,LAG,X )                       
      GO TO 100                                                         
C                                                                       
cc   20 CALL  REDLAG( K )                                                 
   20 CONTINUE
      DO 21 I=1,K
         L1(I) = LG1(1,I)
         L2(I) = LG1(2,I)
         L3(I) = LG1(3,I)
   21 CONTINUE
cc      CALL  REDUCT( SETX2,Z,D,NMK,0,K,MJ1,LAG,X )                       
      CALL  REDUCT( SETX2,Z,NMK,0,K,MJ1,LAG,X )                       
      GO TO 100                                                         
C                                                                       
cc   30 CALL  SETLAG( K )                                                 
cc      CALL  REDUCT( SETX2,Z,D,NMK,0,K,MJ1,LAG,X )                       
cxx   30 CALL  SETLAG( K,LG2(1),LG2(2),LG2(3),LG2(4),LG2(5) )
   30 CALL  SETLAG( KK,LG2(1),LG2(2),LG2(3),LG2(4),LG2(5) )
      DO 31 I=1,K
         LG1(1,I) = L1(I)
         LG1(2,I) = L2(I)
         LG1(3,I) = L3(I)
   31 CONTINUE
      CALL  REDUCT( SETX2,Z,NMK,0,K,MJ1,LAG,X )                       
      GO TO 100                                                         
C                                                                       
cc   40 CALL  REDUCT( SETX4,Z,D,NMK,0,K,MJ1,LAG,X )                       
   40 CALL  REDUCT( SETX4,Z,NMK,0,K,MJ1,LAG,X )                       
      GO TO 100                                                         
C                                                                       
cc   50 CALL  REDREG( K )                                                 
cx   50 CONTINUE
cx      DO 51 J=1,K+1
cx         DO 51 I=1,20
cx            II = (I-1)*4+1
cx            G(I,J) = CHAR(F(II,J)) // CHAR(F(II+1,J))
cx     *                             // CHAR(F(II+2,J)) //CHAR(F(II+3,J))
cx   51 CONTINUE
cc      CALL  REDUCT( SETX5,Z,D,NMK,0,K,MJ1,LAG,X )                       
cx      CALL  REDUCT( SETX5,Z,NMK,0,K,MJ1,LAG,X )                       
cx      GO TO 100                                                         
C                                                                       
cx   60 CONTINUE                                                          
cx	SD0 = SUM
cx      CONST = CNST
cc      SD0 = 0.D0                                                        
cc      DO 65  I=1,N                                                      
cc   65 SD0 = SD0 + Z(I)*Z(I)                                             
cc      SD0 = SD0 / N                                                     
cc      READ( 5,9 )     CONST                                             
cc      WRITE( 6,12 )   SD0 , CONST                                       
cc      CALL  REDLAG( K )                                                 
cx      DO 66 I=1,K
cx         L1(I) = LG1(1,I)
cx         L2(I) = LG1(2,I)
cx         L3(I) = LG1(3,I)
cx   66 CONTINUE
cc      CALL  REDUCT( SETX6,Z,D,NMK,0,K,MJ1,LAG,X )                       
cx      CALL  REDUCT( SETX6,Z,NMK,0,K,MJ1,LAG,X )                       
cx      GO TO 100                                                         
C                                                                       
cc   70 CALL  REDUCT( SETX7,Z,D,NMK,0,K,MJ1,LAG,X )                       
cc   70 CALL  REDUCT( SETX7,Z,NMK,0,K,MJ1,LAG,X )                       
C                                                                       
  100 CONTINUE                                                          
cc	CLOSE( MT )
C                                                                       
C          ---------------                                              
C          MAICE PROCEDURE                                              
C          ---------------                                              
C                                                                       
cc      CALL  ARMFIT( X,K,LAG,NMK,ISW,TITLE,MJ1,A,SD,M )                  
cx      IFG=0
cx      CALL ARMFIT( X,K,LAG,NMK,ISW,MJ1,A1,M,SD,AIC,DIC,SDM,AICM,
cx     *             IFG,LU )
      CALL ARMFIT( X,K,LAG,NMK,ISW,MJ1,A1,M,SD,AIC,DIC,SDM,AICM)
C                                                                       
C          ------------------                                           
C          BAYESIAN PROCEDURE                                           
C          ------------------                                           
C                                                                       
cc      CALL  SBBAYS( X,D,K,NMK,IPR,MJ1,A,SD )                            
      CALL  SBBAYS( X,K,NMK,IPR,MJ1,A2,SDB,EK,AICB,IND,C,C1,C2,B,
     *OEIC,ESUM,OMEAN,OM  )
              if( k.eq.32 ) return
C                                                                       
cc      IF( IMODEL .EQ. 1 )  CALL  NRASPE( SD,A,B,K,0,121,TITLE )         
      IF( IMODEL .EQ. 1 )  CALL  NRASPE( SDB,A2,B,K,0,120,SXX )
C                                                                       
      NPS = LAG+1                                                       
C                                                                       
C          -------------------------                                    
C          PREDICTION ERROR CHECKING                                    
C          -------------------------                                    
C                                                                       
cx      GO TO ( 110,120,120,150,130,140 ), IMODEL                         
cxx      GO TO ( 110,120,120,150 ), IMODEL                         
      IF ( IMODEL .EQ. 1 ) GO TO 110 
      IF ( IMODEL .EQ. 2 ) GO TO 120
      IF ( IMODEL .EQ. 3 ) GO TO 120
      GO TO 150                                                         
C                                                                       
cc  110 CALL CHECK( PRDCT1,Z,A,K,0,IL,NPS,N,0,MJ,E )                      
cxx  110 CALL CHECK( PRDCT1,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
  110 CALL CHECK( PRDCT1,Z,A2,K,0,IL,NPS,N,MJ,E,FA,EMEAN,VARI,SKEW,
     *            PEAK,COV,MJ2 )
      GO TO 150                                                         
C                                                                       
cc  120 CALL  CHECK( PRDCT2,Z,A,K,0,IL,NPS,N,0,MJ,E )                     
cxx  120 CALL  CHECK( PRDCT2,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
  120 CALL  CHECK( PRDCT2,Z,A2,K,0,IL,NPS,N,MJ,E,FA,EMEAN,VARI,SKEW,
     *             PEAK,COV,MJ2 )
      GO TO 150                                                         
C                                                                       
cc  130 CALL CHECK( PRDCT3,Z,A,K,0,IL,NPS,N,0,MJ,E )                      
cx  130 CALL CHECK( PRDCT3,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
cx     *            PEAK,COV,MJ2 )
C                                                                       
cc  140 CALL  CHECK( PRDCT6,Z,A,K,0,IL,NPS,N,0,MJ,E )                     
cx  140 CALL  CHECK( PRDCT6,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
cx     *             PEAK,COV,MJ2 )
C                                                                       
  150 CONTINUE                                                          
      RETURN
C                                                                       
cxx    1 FORMAT( 16I5 )                                                    
cxx    2 FORMAT( /1H ,'ORIGINAL DATA INPUT DEVICE  MT =',I4 )              
cxx    3 FORMAT( ' PROGRAM TIMSAC 78.1.3',/,'   SCALAR TIME SERIES MODEL FI
cxx     *TTING;     BAYESIAN PROCEDURE ( ALL SUBSET REGRESSIN TYPE )' )    
cxx    4 FORMAT( //1H ,'MODEL TYPE',I2,/,'   < AUTOREGRESSIVE MODEL >',/,  
cxx     11H ,10X,'Z(I) = A(1)*Z(I-1) + A(2)*Z(I-2) + ... + A(M)*Z(I-M) + E(
cxx     2I)' )                                                             
cxx    5 FORMAT( //1H ,'MODEL TYPE',I2,/,'   < NON-LINEAR MODEL >',/,1H ,10
cxx     1X,'Z(I) = A(1)*Y(I,1) + A(2)*Y(I,2) + ... + A(K)*Y(I,K) + E(I)' ) 
cxx    6 FORMAT( 1H ,2X,'WHERE',/,11X,'M:     ORDER OF THE MODEL',/,11X,'E(
cxx     1I):  GAUSSIAN WHITE NOISE WITH MEAN 0  AND  VARIANCE SD(M).' )    
cxx    7 FORMAT( 1H ,'FITTING UP TO THE ORDER',I3,2X,'IS TRIED' )          
cxx    8 FORMAT( 1H ,'MAXIMUM LAG =',I4,/,' NUMBER OF REGRESSORS =',I4 )   
cxx    9 FORMAT( F10.0 )                                                   
cxx   11 FORMAT( //1H ,'MODEL TYPE',I2,/,'   < EXPONENTIALLY DAMPED NON-LIN
cxx     *EAR MODEL >',/,1H ,10X,'Z(I) = A(1)*Y(I,1) + A(2)*Y(I,2) + ... + A
cxx     *(K)*Y(I,K) + E(I)' )                                              
cxx   12 FORMAT( 1H ,'SIGMA2 =',D15.5,5X,'CONST =',D15.5 )                 
C                                                                       
      E N D                                                             
cc      REAL FUNCTION  BICOEF * 8( K,J )                                  
cxx      REAL*8 FUNCTION  BICOEF( K,J )                                  
      DOUBLE PRECISION FUNCTION BICOEF( K,J )
C                                                                       
C     THIS FUNCTION RETURNS BINOMIAL COEFFICIENTS                       
C                                                                       
C          F(K,J) = K]/(J]*(K-J)])                                      
C                                                                       
C       INPUTS:                                                         
C          K:     NUMBER OF OBJECTS                                     
C          J:     NUMBER OF OBJECTS TAKEN                               
C                                                                       
C       OUTPUT:                                                         
C          F:     NUMBER OF COMBINATIONS OF SELECTING J OBJECTS FROM    
C                 A SET OF K OBJECTS                                    
C                                                                       
cxx      IMPLICIT REAL * 8 ( A-H , O-Z )
      REAL(8) :: SUM, DI
C                                                                       
      KMJ = K-J                                                         
      SUM = 0.D0                                                        
      DO 10   I=1,K                                                     
      DI = I                                                            
cxx   10 SUM = SUM + DLOG( DI )
      SUM = SUM + DLOG( DI )
   10 CONTINUE                                            
C                                                                       
      IF( J .EQ. 0 )   GO TO 30                                         
      DO 20   I=1,J                                                     
      DI = I                                                            
cxx   20 SUM = SUM - DLOG( DI )
      SUM = SUM - DLOG( DI )                                            
   20 CONTINUE
C                                                                       
   30 IF( KMJ .EQ. 0 )   GO TO 50                                       
      DO 40   I=1,KMJ                                                   
      DI = I                                                            
cxx   40 SUM = SUM - DLOG( DI )                                            
      SUM = SUM - DLOG( DI )
   40 CONTINUE
C                                                                       
   50 BICOEF = DEXP( SUM )                                              
      RETURN                                                            
C                                                                       
      END
cc      SUBROUTINE  CHECK( PRDCT,X,A,K,L,IL,NPS,NPE,IPR,MJ,E )            
cxx      SUBROUTINE  CHECK( PRDCT,X,A,K,L,IL,NPS,NPE,IPR,MJ,E,F,EMEAN,VARI,
      SUBROUTINE  CHECK( PRDCT,X,A,K,L,IL,NPS,NPE,MJ,E,F,EMEAN,VARI,
     *                   SKEW,PEAK,COV,MJ2 )            
C                                                                       
C     THIS SUBROUTINE DRAWS HISTGRAMS AND AUTOCOVARIANCE FUNCTION OF ORI
C     DATA OR PREDICTION ERRORS.                                        
C                                                                       
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             GRAPH1                                                    
C             GRAPH2                                                    
C             MOMENT                                                    
C             (PRDCT)                                                   
C       ----------------------------------------------------------------
C       INPUTS:                                                         
C          PRDCT:  EXTERNAL SUBROUTINE DESIGNATION                      
C          X:      ORIGINAL DATA                                        
C          A:      REGRESSION COEFFICIENTS                              
C          K:      NUMBER OF REGRESSORS                                 
C          L:      MA-ORDER ( THIS ARGUMENT IS ONLY USED FOR THE CHECKIN
C                  OF AR-MA MODEL)                                      
C          IL:     MAXIMUM SPAN OF LONG RANGE PREDICTION                
C                  =0     ANALYSIS OF ORIGINAL DATA                     
C                  >0     ANALYSIS OF MULTI-STEP (UP TO IL) PREDICTION E
C          NPS:    PREDICTION STARTING POSITION                         
C          NPE:    PREDICTION ENDING POSITION                           
C          IPR:    =0  MATRIX OF SEVERAL STEP AHEAD PREDICTION ERRORS SU
C                  =1  MATRIX OF SEVERAL STEP AHEAD PREDICTION ERRORS IS
C                      OUT                                              
C          MJ:     ABSOLUTE DIMENSION OF E IN THE MAIN PROGRAM          
C                                                                       
C       OUTPUT:                                                         
C          E:      SEVERAL-STEPS PREDICTION ERRORS                      
C                                                                       
C                                                                       
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL * 4  X , F , E                                               
cx      DIMENSION  X(1) , E(MJ,1)                                         
cc      DIMENSION  A(1) , COV(120)                                        
cxx      DIMENSION  X(NPE) , E(MJ,IL)
cxx      DIMENSION  A(K) , COV(MJ2)                                        
cxx      DIMENSION  EMEAN(IL), VARI(IL), SKEW(IL), PEAK(IL)
cc      COMMON  /COMXX/ F(2000)                                           
cxx      DIMENSION  F(NPE-NPS+1,IL)
      INTEGER :: K, L, IL, NPS, NPE, MJ, MJ2 
      REAL(8) :: X(NPE), A(K), E(MJ,IL), F(NPE-NPS+1,IL), EMEAN(IL),
     1           VARI(IL), SKEW(IL), PEAK(IL), COV(MJ2)
      REAL(8) :: SUM, COV1, SD
C                                                                       
C                                                                       
      ISTEP = 1                                                         
      ISW = IL                                                          
      LAGH = 100                                                        
      N = NPE - NPS - 1                                                 
      IF( LAGH .GE. N )     LAGH = N - 1                                
      LAG1 = LAGH + 1                                                   
      NMK = N - K                                                       
      IF( ISW .GT. 0 )     GO TO 20                                     
C                                                                       
      DO 10  I=NPS,NPE                                                  
cxx   10 E(I,1) = X(I) 
      E(I,1) = X(I)
   10 CONTINUE                                                    
      IL = 1                                                            
      GO TO 36                                                          
C                                                                       
C       ---  SEVERAL STEP AHEAD PREDICITON  ---                         
C                                                                       
   20 CONTINUE                                                          
      CALL  PRDCT( X,A,K,L,IL,NPS,NPE,MJ,E )                            
C                                                                       
C       ---  PREDICTION ERROR  ---                                      
C                                                                       
cxx      DO 30  II=NPS,NPE                                                 
      DO 31  II=NPS,NPE
         I = NPE-II+NPS                                                 
         DO 30  J=1,IL                                                  
         JJ = I-J+1                                                     
cxx   30 E(I,J) = X(I) - E(JJ,J)
      E(I,J) = X(I) - E(JJ,J)
   30 CONTINUE
   31 CONTINUE                                           
      IF( IL .EQ. 1 )     GO TO 34                                      
      DO 35  J=2,IL                                                     
      JJ = J-1                                                          
cxx      DO 35  I=1,JJ                                                     
      DO 33  I=1,JJ 
      II = I+NPS-1                                                      
cxx   35 E(II,J) = 0.D0                                                    
      E(II,J) = 0.D0
   33 CONTINUE
   35 CONTINUE
   34 CONTINUE                                                          
cc      WRITE( 6,690 )                                                    
cc      IF( IPR .EQ. 0 )   GO TO 36                                       
cc      WRITE( 6,670 )     (I,I=1,IL)                                     
cc      DO 37  I=NPS,NPE                                                  
cc   37 WRITE( 6,640 )     I , (E(I,J),J=1,IL)                            
   36 CONTINUE                                                          
C                                                                       
C       ---  MOMENT COMPUTATION  ---                                    
C                                                                       
      DO 50  KK=1,IL                                                    
C                                                                       
      II = NPS+KK-1                                                     
      DO 40  I=II,NPE                                                   
      J = I - II + 1                                                    
cc   40 F(J) = E(I,KK)                                                    
cxx   40 F(J,KK) = E(I,KK)
      F(J,KK) = E(I,KK)
   40 CONTINUE
      NMK = NPE-NPS-(KK-2)                                              
C                                                                       
cc      CALL  MOMENT( F,NMK,EMEAN,VARI,SKEW,PEAK )                        
      CALL  MOMENT( F(1,KK),NMK,EMEAN(KK),VARI(KK),SKEW(KK),PEAK(KK) )
C                                                                       
cc      IF( ISW .GT. 1 )     WRITE( 6,610 )   KK                          
cc      IF( ISW .EQ. 0 )     WRITE( 6,650 )                               
cc      WRITE( 6,620 )     EMEAN , VARI , SKEW , PEAK                     
C                                                                       
C       ---  HISTOGRAM OF F(I)  ---                                     
C                                                                       
cc      SIG = DSQRT(VARI)                                                 
cc      CALL  GRAPH1( F,1,NMK,SIG )                                       
   50 CONTINUE                                                          
C                                                                       
C       ---  AUTOCORRELATION FUNCTION COMPUTATION  ---                  
C                                                                       
      DO 100  KK=1,IL                                                   
C                                                                       
      DO 70   II=1,LAG1                                                 
      JJ = NPS + KK - 1                                                 
      IE = NPE - II + 1                                                 
      SUM = 0.D0                                                        
      DO 60   I=JJ,IE                                                   
      J = I + II- 1                                                     
cxx   60 SUM = SUM + E(I,KK)*E(J,KK)
      SUM = SUM + E(I,KK)*E(J,KK)
   60 CONTINUE                                       
cxx   70 COV(II) = SUM / (NPE-NPS-KK+2)
      COV(II) = SUM / (NPE-NPS-KK+2)
   70 CONTINUE                                    
C                                                                       
      COV1 = COV(1)                                                     
      DO 80   I=1,LAG1                                                  
cxx   80 COV(I) = COV(I) / COV1                                            
      COV(I) = COV(I) / COV1
   80 CONTINUE
C                                                                       
cc      IF( ISW .GT. 0 )   WRITE( 6,630 )     KK                          
cc      IF( ISW .EQ. 0 )   WRITE( 6,660 )                                 
cc      WRITE( 6,680 )     (COV(I),I=1,50)                                
C                                                                       
C       ---  AUTOCORRELATION FUNCTION DISPLAY  ---                      
C                                                                       
      SD = NPE - NPS + 1                                                
      SD = DSQRT( 1.D0/SD )                                             
cc      CALL  GRAPH2( COV,LAG1,SD )                                       
C                                                                       
      IF( ISTEP .EQ. KK )     GO TO 110                                 
C                                                                       
  100 CONTINUE                                                          
C                                                                       
C                                                                       
  110 CONTINUE                                                          
cxx      RETURN                                                            
cxx  610 FORMAT( //1H ,I3,'-STEP AHEAD PREDICTION ERROR' )                 
cxx  620 FORMAT( 1H ,'MEAN       =',D15.8,/,' VARIANCE   =',D15.8,/,' SKEWN
cxx     1ESS   =',D15.8,/,' PEAKEDNESS =',D15.8 )                          
cxx  630 FORMAT( //1H ,'AUTOCORRELATION FUNCTION OF ',I3,'-STEP AHEAD PREDI
cxx     1CTION ERROR' )                                                    
cxx  640 FORMAT( 1H ,I5,5X,10D12.4 )                                       
cxx  650 FORMAT( //,' ORIGINAL DATA' )                                     
cxx  660 FORMAT( //1H ,'AUTOCORRELATION FUNCTION OF ORIGINAL DATA' )       
cxx  670 FORMAT( 1H ,10X,10(4X,'J =',I2,3X) )                              
cxx  680 FORMAT( 1H ,10D13.5 )                                             
cxx  690 FORMAT( ///1H ,45(1H-),2X,'<< J-STEP AHEAD PREDICTION ERROR >>',2X
cxx     1,45(1H-) )                                                        
      END                                                               
      SUBROUTINE  MOMENT( X,N,F1,F2,F3,F4 )                             
C                                                                       
C          +--------------------+                                       
C          ! MOMENT COMPUTATION !                                       
C          +--------------------+                                       
C                                                                       
C     THIS SUBROUTINE COMPUTES MOMENTS.                                 
C                                                                       
C       INPUTS:                                                         
C          X:     ORIGINAL DATA VECTOR                                  
C          N:     DATA LENGTH                                           
C                                                                       
C       OUTPUTS:                                                        
C          F1:    MEAN OF X                                             
C          F2:    VARIANCE OF X                                         
C          F3:    SKEWNESS OF X                                         
C          F4:    PEAKEDNESS OF X                                       
C                                                                       
cxx      IMPLICIT  REAL * 8  ( F )                                         
CC      DIMENSION  X(1)                                                   
cx	REAL*8  X(1)
cxx      REAL*8  X(N)
      INTEGER :: N
      REAL(8) :: X(N), F1, F2, F3, F4
      REAL(8) :: FN, FSUM, FF
C                                                                       
      FN = N                                                            
      FSUM = 0.D0                                                       
      DO  10     I=1,N                                                  
cxx   10 FSUM = FSUM + X(I)                                                
      FSUM = FSUM + X(I)
   10 CONTINUE
C                                                                       
      F1 = FSUM / FN                                                    
C                                                                       
      F2 = 0.D0                                                         
      F3 = 0.D0                                                         
      F4 = 0.D0                                                         
      DO  20     I=1,N                                                  
      FF = X(I) - F1                                                    
      F2 = F2 + FF*FF                                                   
      F3 = F3 + FF**3                                                   
cxx   20 F4 = F4 + FF**4                                                   
      F4 = F4 + FF**4
   20 CONTINUE
C                                                                       
      F2 = F2 / FN                                                      
      F3 = F3 / (FN*F2*DSQRT(F2))                                       
      F4 = F4 / (FN*F2*F2)                                              
C                                                                       
      RETURN                                                            
      END                                                               
      SUBROUTINE  PRDCT1( Z,A,M,L,IL,NPS,NPE,MJ,EZ )                    
C                                                                       
C     THIS SUBROUTINE COMPUTES SEVARAL STEP AHEAD PREDICTION VALUE OF AN
C     AUTOREGRESSIVE MOVING AVERAGE MODEL.                              
C                                                                       
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA VECTOR                                 
C          A:      AR-MA COEFFICIENTS                                   
C          M:      AR-ORDER                                             
C          L:      MA-ORDER                                             
C          IL:     MAXIMUM SPAN OF LONG RANGE PREDICTION                
C          NPS:    PREDICTION STARTING POSITION                         
C          NPE:    PREDICTION ENDING POSITION                           
C          MJ:     ABSOLUTE DIMENSION OF EZ                             
C                                                                       
C       OUTPUT:                                                         
C          EZ:     PREDICTION VALUE MATRIX                              
C                                                                       
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL * 4  Z(1) , EZ(MJ,1)                                         
cx      DIMENSION  Z(1) , EZ(MJ,1)
cx	DIMENSION  A(1)
cxx      DIMENSION  Z(NPE) , EZ(MJ,IL)
cxx      DIMENSION  A(M)                                                   
      INTEGER :: M, L, IL, NPS, NPE, MJ
      REAL(8) :: Z(NPE), A(M), EZ(MJ,IL)
      REAL(8) :: SUM
C                                                                       
C                                                                       
      DO  100     II=NPS,NPE                                            
C                                                                       
cxx      DO  90     KK=1,IL                                                
      DO  91     KK=1,IL 
      KKM1 = KK - 1                                                     
      SUM = 0.D0                                                        
      IF( KK .EQ. 1 )     GO TO 30                                      
      DO  20     I=1,KKM1                                               
      KI = KK - I                                                       
cxx   20 SUM = SUM + A(I)*EZ(II,KI)                                        
      SUM = SUM + A(I)*EZ(II,KI)
   20 CONTINUE
   30 IF( KK .GT. M )     GO TO 50                                      
      DO  40     I=KK,M                                                 
      I1 = II + KKM1 - I                                                
cxx   40 SUM = SUM + A(I)*Z(I1)                                            
      SUM = SUM + A(I)*Z(I1)
   40 CONTINUE
C                                                                       
   50 IF( L .LE. 0 )     GO TO 90                                       
      IF( KK .GT. L )     GO TO 90                                      
      DO  60     I=KK,L                                                 
      I1 = M + I                                                        
      I2 = II + KKM1 - I                                                
      IF( I2 .GE. II )     GO TO 60                                     
      SUM = SUM + A(I1)*(Z(I2)-EZ(I2,1))                                
   60 CONTINUE                                                          
   90 EZ(II,KK) = SUM
   91 CONTINUE                                                   
  100 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END
      SUBROUTINE  PRDCT2( Z,A,K,L,IL,NPS,NPE,MJ1,EZ )                   
C                                                                       
C     THIS SUBROUTINE COMPUTES SEVERAL STEPS AHEAD PREDICTION VALUES OF 
C     NON-LINEAR REGRESSION MODEL.                                      
C                                                                       
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA VECTOR                                 
C          A:      VECTOR OF AR-COEFFICIENTS                            
C          K:      ORDER OF THE AR MODEL                                
C          L:      THIS DUMMY VARIABLE IS NOT REFERENCED IN THIS SUBROUT
C          IL:     MAXIMUM SPAN OF LONG RANGE PREDICTION                
C          NPS:    PREDICTION STARTING POSITION                         
C          NPE:    PREDICTION ENDING POSITION                           
C          MJ1:    ABSOLUTE DIMENSION OF EZ                             
C                                                                       
C       OUTPUT:                                                         
C          EZ:     PREDICTION VALUE MATRIX                              
C                                                                       
CC      IMPLICIT  REAL*8 ( A-D,O-Y )                                      
cxx      IMPLICIT  REAL*8 ( A-H,O-Z ) 
cc      DIMENSION  Z(1) , A(1) , EZ(MJ1,1) , Y(20)                        
cx      DIMENSION  Z(1) , A(1) , EZ(MJ1,1) , Y(IL)                        
cxx      DIMENSION Z(NPE) , A(K) , EZ(MJ1,IL) , Y(IL)
cc      REAL * 4   CSTDMY, SD0DMY
      INTEGER :: K, IL, NPS, NPE, MJ1
      REAL(8) :: Z(NPE), A(K), EZ(MJ1,IL)
      REAL(8) :: Y(IL), CSTDMY, SD0DMY, SUM, XX, X
      COMMON     / BBB /  LAG1(50) , LAG2(50) , LAG3(50), CSTDMY, SD0DMY
CC      INTEGER  RETURN                                                   
C                                                                       
cxx      DO 100  II=NPS,NPE                                                
      DO 110  II=NPS,NPE
         DO 50  J1=1,IL                                                 
            JJ = J1-1                                                   
            SUM = 0.D0                                                  
            DO 40  J=1,K                                                
               XX = 1.D0                                                
               LAG = LAG1(J)                                            
CC               ASSIGN 10 TO RETURN                                      
CC               GO TO 200                                                
CC   10          XX = XX*X
               X = 1.D0
               IF ( LAG .GT. 0 ) THEN
                  I = II+JJ-LAG
                  IF ( I .GE. II ) THEN
                     I = I - II + 1
                     X = Y(I)
                  ELSE
                     X = Z(I)
                  END IF
               END IF
               XX = XX*X                                               
C                                                
               LAG = LAG2(J)                                            
CC               ASSIGN 20 TO RETURN                                      
CC               GO TO 200                                                                                               
CC   20          XX = XX*X
               X = 1.D0
               IF ( LAG .GT. 0 ) THEN
                  I = II+JJ-LAG
                  IF ( I .GE. II ) THEN
                     I = I - II + 1
                     X = Y(I)
                  ELSE
                     X = Z(I)
                  END IF
               END IF
               XX = XX*X
C
               LAG = LAG3(J)                                            
CC               ASSIGN 30 TO RETURN                                      
CC               GO TO 200                                                
CC   30          XX = XX*X
               X = 1.D0
               IF ( LAG .GT. 0 ) THEN
                  I = II+JJ-LAG
                  IF ( I .GE. II ) THEN
                     I = I - II + 1
                     X = Y(I)
                  ELSE
                     X = Z(I)
                  END IF
               END IF
               XX = XX*X
C                                                                       
cxx   40       SUM = SUM + A(J)*XX 
            SUM = SUM + A(J)*XX 
   40       CONTINUE                                        
cxx   50    Y(J1) = SUM
         Y(J1) = SUM
   50    CONTINUE                                                    
C                                                                       
      DO  100  J=1,IL                                                   
      EZ(II,J) = Y(J)                                                   
  100 CONTINUE
  110 CONTINUE                                     
CC      GO TO 300
C        -----  INTERNAL SUBROUTINE  -----
C
CC  200 X = 1.D0
CC      IF ( LAG .LE. 0 )     GO TO 220
CC      I = II+JJ-LAG
CC      IF ( I .GE. II )      GO TO 210
CC      X = Z(I)
CC      GO TO 220
CC  210 I = I - II + 1
CC      X = Y(I)
CC  220 GO TO RETURN, ( 10,20,30 )
C        ----------------------------------
C                                                                                                                                                              
CC  300 RETURN
C    L : DUMMY
      L = L
      RETURN
      END                                                               
cc      SUBROUTINE  SETLAG( K )                                           
      SUBROUTINE  SETLAG( K,LAG1,LAG2,LAG3,LAG4,LAG5 )                                           
C                                                                       
C     THIS SUBROUTINE  PREPARES SPECIFICATION OF REGRESSORS (L1(I),L2(I)
C     (I=1,...,K) FOR THE FITTING OF (POLYNOMIAL TYPE) NON-LINEAR MODEL.
C     THE OUTPUTS ARE USED AS THE INPUTS TO SUBROUTINE SETX2.           
C                                                                       
C       INPUTS:                                                         
C          LAG1:    MAXIMUM TIME LAG OF LINEAR TERM                     
C          LAG2:    MAXIMUM TIME LAG OF SQUARED TERM                    
C          LAG3:    MAXIMUM TIME LAG OF QUADRATIC CROSS TERM            
C          LAG4:    MAXIMUM TIME LAG OF CUBIC TERM                      
C          LAG5:    MAXIMUM TIME LAG OF CUBIC CROSS TERM                
C                                                                       
C       OUTPUTS:                                                        
C          K:       NUMBER OF REGRESSORS                                
C          (L1(I),L2(I),L3(I))  (I=1,K):     SPECIFICATION OF REGRESSORS
C                                                                       
C              ......................................................   
C              I-TH REGRESSOR IS DEFINED BY                             
C                   Z(N-L1(I)) * Z(N-L2(I)) * Z(N-L3(I))                
C              WHERE  0-LAG TERM Z(N-0) IS REPLACED BY THE CONSTANT 1.  
C              ......................................................   
C                                                                       
cc      REAL * 4   CSTDMY,SD0DMY
cxx      REAL * 8   CSTDMY,SD0DMY
      INTEGER :: K, LAG1, LAG2, LAG3, LAG4, LAG5
      REAL(8) :: CSTDMY, SD0DMY
      COMMON     / BBB /  L1(50),L2(50),L3(50),CSTDMY,SD0DMY
C                                                                       
cc      READ( 5,1 )     LAG1,LAG2,LAG3,LAG4,LAG5                          
cc      WRITE( 6,2 )    LAG1,LAG2,LAG3,LAG4,LAG5                          
      IF(LAG1.LE.0)  GO TO 15                                           
      DO 10  I=1,LAG1                                                   
      L1(I) = I                                                         
      L2(I) = 0                                                         
cxx   10 L3(I) = 0
      L3(I) = 0
   10 CONTINUE                                                         
   15 K = LAG1                                                          
C                                                                       
      IF(LAG2.LE.0)   GO TO 30                                          
      DO 20  I=1,LAG2                                                   
      K = K+1                                                           
      L1(K) = I                                                         
      L2(K) = I                                                         
cxx   20 L3(K) = 0                                                         
      L3(K) = 0
   20 CONTINUE
C                                                                       
   30 IF(LAG3.LE.1)   GO TO 50                                          
      LL = LAG3-1                                                       
cxx      DO 40  I=1,LL
      DO 41  I=1,LL                                                     
         I1 = I+1                                                       
         DO 40  J=I1,LAG3                                               
         K = K+1                                                        
         L1(K) = I                                                      
         L2(K) = J                                                      
cxx   40    L3(K) = 0
         L3(K) = 0
   40    CONTINUE
   41 CONTINUE
   50 M  = K                                                            
C                                                                       
      IF(LAG4.LE.0)   GO TO 65                                          
      DO 60  I=1,LAG4                                                   
      K = K+1                                                           
      L1(K) = I                                                         
      L2(K) = I                                                         
cxx   60 L3(K) = I                                                         
      L3(K) = I
   60 CONTINUE
   65 CONTINUE                                                          
C                                                                       
      IF(LAG5.LE.1)   GO TO 80                                          
cxx      DO 70  I=1,LAG5                                                   
cxx         DO 70  J=I,LAG5                                                
      DO 72  I=1,LAG5                                                   
         DO 71  J=I,LAG5
            DO 70  L=J,LAG5                                             
            IF(I.EQ.J .AND. J.EQ.L)  GO TO 70                           
            K = K+1                                                     
            L1(K) = I                                                   
            L2(K) = J                                                   
            L3(K) = L                                                   
   70       CONTINUE
   71    CONTINUE
   72 CONTINUE
C                                                                       
cc   80 WRITE( 6,3 )                                                      
   80 CONTINUE
cxx      IF( LAG1 .EQ. 0 )  GO TO 100                                      
cxx      DO 90  I=1,LAG1                                                   
cc   90 WRITE( 6,4 )     I , L1(I)                                        
cxx   90 CONTINUE
cxx  100 J = LAG1+1                                                        
cxx      IF( LAG2+LAG3 .EQ. 0 )   GO TO 120                                
cxx      DO 110  I=J,M                                                     
cc  110 WRITE( 6,5 )     I , L1(I) , L2(I)                                
cxx  110 CONTINUE
cxx  120 J = M+1                                                           
cxx      IF( LAG4+LAG5 .EQ. 0 )   GO TO 140                                
cxx      DO 130  I=J,K                                                     
cc  130 WRITE( 6,6 )     I ,L1(I) , L2(I) ,L3(I)                          
cxx  130 CONTINUE
cxx  140 CONTINUE                                                          
C                                                                       
      RETURN                                                            
cxx    1 FORMAT( 16I5 )                                                    
cxx    2 FORMAT( /1H ,'LAG1 =',I3,5X,'LAG2 =',I3,5X,'LAG3 =',I3,5X,
cxx     *        'LAG4 =',I3,5X,'LAG5 =',I3 )
cxx    3 FORMAT( 1H ,4X,'M',5X,'REGRESSOR  Y(I,M)' )                       
cxx    4 FORMAT( 1H ,I5,5X,'Z(I-',I2,')' )                                 
cxx    5 FORMAT( 1H ,I5,5X,'Z(I-',I2,') * Z(I-',I2,')' )                   
cxx    6 FORMAT( 1H ,I5,5X,'Z(I-',I2,') * Z(I-',I2,') * Z(I-',I2,')' )     
C                                                                       
      END                                                               
      SUBROUTINE  SETX2( Z,N0,L,K,MJ1,JSW,LAG,X )                       
C                                                                       
C          +----------------------------------------+                   
C          ! MATRIX X SET UP (FOR NON-LINEAR MODEL) !                   
C          +----------------------------------------+                   
C                                                                       
C     THIS SUBROUTINE PREPARES DATA MATRIX X FROM DATA VECTOR Z(I) (I=N0
C     N0+K+LAG) FOR THE FITTING OF NON-LINEAR AUTOREGRESSIVE MODEL.  X I
C     USED AS THE INPUT TO SUBROUTINE HUSHLD.                           
C                                                                       
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA VECTOR                                 
C          N0:     INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATIO
C                  (NEW OBSERVATION STARTS AT N0+LAG+1 AND ENDS AT N0+LA
C          L:      DIMENSION OF THE VECTOR OF NEW OBSERVATIONS          
C          K:      NUMBER OF REGRESSORS                                 
C          MJ1:    ABSOLUTE DIMENSION OF X                              
C          JSW:    =0   TO CONSTRUCT INITIAL L*(K+1) DATA MATRIX        
C                  =1   TO AUGMENT ORIGINAL (K+1)*(K+1) MATRIX X BY AN  
C                       L*(K+1) DATA MATRIX OF ADDITIONAL OBSERVATIONS  
C          LAG:    MAXIMUM TIME LAG                                     
C          KSW:    THIS DUMMY VARIABLE IS NOT REFERENCED IN THIS SUBROUT
C                                                                       
C--  THE FOLLOWING VARIABLE SPECIFICATION IS GIVEN EITHER BY REDLAG OR S
C         (L1(I) , L2(I) , L3(I))  (I=1,K)                              
C                                                                       
C               I-TH REGRESSOR IS DEFINED BY                            
C                    Z(N-L1(I)) * Z(N-L2(I)) * Z(N-L3(I))               
C               WHERE 0-LAG TERM Z(N-0) IS AUTOMATICALLY REPLACED BY CON
C                                                                       
C       OUTPUT:                                                         
C          X:      L*(K+1) MATRIX           IF  JSW = 0                 
C                  (K+1+L)*(K+1) MATRIX     IF  JSW = 1                 
C                                                                       
CC      REAL * 8  X(MJ1,1)                                                
CC      DIMENSION  Z(1)                                                   
cc      REAL * 4   CSTDMY, SD0DMY
cx      REAL * 8  X(MJ1,1) ,  Z(1) , ZTEM
cxx      REAL * 8  X(MJ1,K+1) ,  Z(N0+LAG+L) , ZTEM
cxx      REAL * 8  CSTDMY, SD0DMY
      INTEGER :: N0, L, K, MJ1, JSW, LAG 
      REAL(8) :: Z(N0+LAG+L), X(MJ1,K+1)
      REAL(8) :: CSTDMY, SD0DMY, ZTEM
      COMMON     / BBB /  L1(50) , L2(50) , L3(50), CSTDMY, SD0DMY  
C                                                                       
      K1 = K + 1                                                        
      I0 = K1*JSW                                                       
      DO  10     I=1,L                                                  
      I1 = I + I0                                                       
      J1 = N0 + LAG + I                                                 
cxx   10 X(I1,K1) = Z(J1)                                                  
      X(I1,K1) = Z(J1)
   10 CONTINUE
C                                                                       
      DO  70     II=1,K                                                 
      LL1 = L1(II)                                                      
      LL2 = L2(II)                                                      
      LL3 = L3(II)                                                      
cxx      DO  60     I=1,L
      DO  61     I=1,L
      ZTEM = 1.D0                                                       
      I1 = I + I0                                                       
      J1 = N0 + LAG + I                                                 
      IF( LL1 .EQ. 0 )     GO TO 40                                     
      M1 = J1 - LL1                                                     
      ZTEM = ZTEM * Z(M1)                                               
   40 IF( LL2 .EQ. 0 )     GO TO 50                                     
      M2 = J1 - LL2                                                     
      ZTEM = ZTEM * Z(M2)                                               
   50 IF( LL3 .EQ. 0 )     GO TO 60                                     
      M3 = J1 - LL3                                                     
      ZTEM = ZTEM * Z(M3)                                               
   60 X(I1,II) = ZTEM
   61 CONTINUE                                                   
   70 CONTINUE                                                          
C                                                                       
      RETURN                                                            
C                                                                       
      E N D                                                             
      SUBROUTINE  SETX4( Z,NO,L,K,MJ1,JSW,LAG,X )                       
C                                                                       
C     THIS SUBROUTINE PREPARES DATA MATRIX X FROM DATA VECTOR Z(I) (I=NO
C     NO+K+L) FOR THE FITTING OF AUTOREGRESSIVE MODEL WITH POLYNOMIAL TY
C     VALUE FUNCTION.  X IS THEN USED AS THE INPUT TO SUBROUTINE HUSHLD.
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA VECTOR                                 
C          NO:     INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATIO
C          L:      DIMENSION OF THE VECTOR OF NEW OBSERVATIONS          
C          K:      NUMBER OF REGRESSORS                                 
C          MJ1:    ABSOLUTE DIMENSION OF X                              
C          JSW:    =0   TO CONSTRUCT INITIAL L*(K+1) DATA MATRIX        
C                  =1   TO AUGMENT ORIGINAL (K+1)*(K+1) MATRIX X BY AN  
C                       L*(K+1) DATA MATRIX OF ADDITIONAL OBSERVATIONS  
C          LAG:    MAXIMUM TIME LAG OF THE MODELS                       
C          N:      DATA LENGTH                                          
C                                                                       
C       OUTPUT:                                                         
C          X:      L*(K+1) MATRIX           IF   JSW = 0                
C                  (K+1+L)*(K+1) MATRIX     IF   JSW = 1                
C                                                                       
C                                                                       
cxx      IMPLICIT  REAL  * 8 (A-H,O-Z)                                     
CC      REAL * 4 Z(1)                                                     
cx      DIMENSION  Z(1)
cx      DIMENSION  X(MJ1,1)                                               
cxx      DIMENSION  X(MJ1,K+1) ,  Z(NO+LAG+L)
      INTEGER :: NO, L, K, MJ1, JSW, LAG                       
      REAL(8) :: Z(NO+LAG+L), X(MJ1,K+1)
      REAL(8) :: BN, Y, XX
cc      COMMON     / AAA /  N , M                                         
      COMMON     / AAA /  N
C                                                                       
C          M:      ORDER OF POLYNOMIAL OF MEAN VALUE FUNCTION           
C                                                                       
      M = K - LAG - 1                                                   
      K1 = K + 1                                                        
      I0 = JSW*K1                                                       
      M1 = M + 1                                                        
      LAG=K-M1                                                          
      BN = 2.D0/(N-LAG-1.D0)                                            
cxx      DO 10  I=1,L                                                      
      DO 11  I=1,L
      Y= BN*(NO+I-1)-1.D0                                               
      XX= 1.D0                                                          
      DO 10 J=1,M1                                                      
      II= I+I0                                                          
      X(II,J) = XX                                                      
cxx   10 XX = XX*Y 
      XX = XX*Y 
   10 CONTINUE
   11 CONTINUE                                                        
C                                                                       
cxx      DO 20  I=1,L 
      DO 21  I=1,L                                                     
      II = I+I0                                                         
      NN = NO+LAG+I                                                     
      X(II,K1) = Z(NN)                                                  
      DO 20   J=1,LAG                                                   
      NN = NN-1                                                         
      JJ = J+M1                                                         
cxx   20 X(II,JJ) = Z(NN)                                                  
      X(II,JJ) = Z(NN)
   20 CONTINUE
   21 CONTINUE
C                                                                       
      RETURN                                                            
C                                                                       
      END                                                               
      SUBROUTINE  SRTMIN( X,N,IX )                                      
C                                                                       
C       THIS SUBROUTINE ARRANGES X(I) (I=1,N) IN ORDER OF INCREASING    
C       MAGNITUDE OF X(I)                                               
C                                                                       
C       INPUTS:                                                         
C          X:   VECTOR                                                  
C          N:   DIMENSION OF THE VECTOR                                 
C       OUTPUTS:                                                        
C          X:   ARRANGED VECTOR                                         
C          IND: INDEX OF ARRANGED VECTOR                                
C                                                                       
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cx      DIMENSION  X(1) , IX(1)                                           
cxx      DIMENSION  X(N) , IX(N)
      INTEGER :: N, IX(N)
      REAL(8) :: X(N), XMIN, XT
C                                                                       
      NM1 = N - 1                                                       
      DO  30     I=1,N                                                  
cxx   30 IX(I) = I                                                         
      IX(I) = I 
   30 CONTINUE
      DO  20     II=1,NM1                                               
      XMIN = X(II)                                                      
      MIN = II                                                          
      DO  10     I=II,N                                                 
      IF( XMIN .LT. X(I) )     GO TO 10                                 
      XMIN = X(I)                                                       
      MIN = I                                                           
   10 CONTINUE                                                          
      IF( XMIN .EQ. X(II) )     GO TO 20                                
      XT = X(II)                                                        
      X(II) = X(MIN)                                                    
      X(MIN) = XT                                                       
      IT = IX(II)                                                       
      IX(II) = IX(MIN)                                                  
      IX(MIN) = IT                                                      
   20 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      E N D                                                             
cc      SUBROUTINE  SUBSPC( B,K,N,IPR,EK )                                
cxx      SUBROUTINE  SUBSPC( B,K,N,IPR,EK,IND,C,C1,C2,OEIC,ESUM1,OMEAN,OM )
      SUBROUTINE  SUBSPC( B,K,N,EK,IND,C,C1,C2,OEIC,ESUM1,OMEAN,OM )
C                                                                       
C       THIS SUBROUTINE PRODUCES BAYESIAN ESTIMATES OF PARTIAL CORRELATI
C       BY CHECKING ALL SUBSET REGRESSION MODELS.                       
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             BICOEF                                                    
C             SRTMIN                                                    
C       ----------------------------------------------------------------
C                                                                       
C         INPUTS:                                                       
C           B:   LEAST SQUARES ESTIMATES OF PARTIAL CORRELATIONS        
C           K:   DIMENSION OF VECTOR A                                  
C           N:   NUMBER OF OBSERVATIONS USED FOR THE ESTIMATION OF A(I) 
C           IPR: =0  TO SUPPRESS THE OUTPUTS                            
C                >0  TO PRINT OUT THE OUTPUTS                           
C                                                                       
C         OUTPUTS:                                                      
C           B(I) (I=1,K):   BAYESIAN ESTIMATES OF PARTIAL CORRELATIONS  
C           EK:   EQUIVALENT NUMBER OF FREE PARAMETERS IN THE BAYESIAN M
C                                                                       
cxx      IMPLICIT  REAL * 8 ( A-H , O-Z )                                  
cc      DIMENSION  B(1) , C(50) , D(50,50)                                
cc      DIMENSION  IND(50) , KND(50) , ESUM(50)                           
cc      DIMENSION  C1(50), C2(50), ESUM1(50)
cx      DIMENSION  B(1) , C(K) , D(K+1,K+1)                                
cxx      DIMENSION  B(K) , C(K) , D(K+1,K+1)                                
cxx      DIMENSION  IND(K) , KND(K+1) , ESUM(K+1)                           
cxx      DIMENSION  C1(K+1), C2(K), ESUM1(K+1)
      INTEGER :: K, N, IND(K)
      REAL(8) :: B(K), EK, C(K), C1(K+1), C2(K), OEIC, ESUM1(K+1),
     1           OMEAN, OM
      INTEGER :: KND(K+1)
      REAL(8) :: BICOEF, D(K+1,K+1), ESUM(K+1), CC, DN, SUM, EIC,
     1           SUMC, EXIC, OSUM
C                                                                       
      CC = 1.D0 + DLOG(2.D0)                                            
      K1 = K + 1                                                        
      DN = N                                                            
cxx      DO 10   I=1,K1                                                    
cxx      ESUM(I) = 0.D0                                                    
cxx      DO 10   J=1,K1                                                    
cxx   10 D(I,J) = 0.D0 
      ESUM(1:K1) = 0.D0
      D(1:K1,1:K1) = 0.D0
C                                                                       
C          SQUARE OF PARTIAL CORRELATIONS ( NORMALISED BY MULTIPLYING N 
C                                                                       
      DO 20   I=1,K                                                     
cxx   20 C(I) = B(I)*B(I)*DN                                               
      C(I) = B(I)*B(I)*DN 
   20 CONTINUE
C                                                                       
C          ARRANGEMENT OF C(I) IN ORDER OF INCREASING MAGNITUDE         
C                                                                       
      CALL  SRTMIN( C,K,IND )                                           
cc      IF( IPR .LE. 1 )     GO TO 60                                     
cc      WRITE( 6,7 )                                                      
cc      DO  50     I=1,K                                                  
cc   50 WRITE( 6,6 )     I , IND(I) , C(I)                                
cc   60 CONTINUE                                                          
C                                                                       
C          FIND THE MINIMUM OF EIC                                      
C                                                                       
      OEIC = CC*K                                                       
      SUM = 0.D0                                                        
      DO 30   I=1,K                                                     
      SUM = SUM + C(I)                                                  
      EIC = SUM + CC*(K-I)                                              
cxx   30 IF( OEIC .GT. EIC )   OEIC = EIC                                  
      IF( OEIC .GT. EIC )   OEIC = EIC 
   30 CONTINUE
cc      WRITE( 6,604 )   OEIC                                             
C                                                                       
C--------  COMPUTATION OF EIC'S OF WHOLE SUBSET REGRESSION MODELS  -----
C                                                                       
C          INITIAL SETTING                                              
C                                                                       
      DO 40   I=1,K                                                     
cxx   40 KND(I) = 0                                                        
      KND(I) = 0
   40 CONTINUE
      KND(K1) = 1                                                       
      SUM = 0.D0                                                        
      SUMC = 0.D0
      M = K                                                             
      IP = 0                                                            
      IQ = 0                                                            
C                                                                       
  100 CONTINUE                                                          
C                                                                       
C          -----  SPECIFICATION OF NEXT SUBSET  -----                   
cxx               DO 110   I=1,K
               DO 111   I=1,K
               IF( KND(I) .EQ. 0 )   GO TO 110                          
               KND(I) = 0                                               
               GO TO 120                                                
  110          KND(I) = 1                                               
  111          CONTINUE
  120          CONTINUE                                                 
C          ------------------------------------------                   
C                                                                       
  130 CONTINUE             
      IF( IP .GT. K )   GO TO 200                                       
      IF( KND(IP+1) .EQ. 0 )   GO TO 140                                
C                                                                       
      IF( IQ .EQ. 0 )   GO TO 165                                       
      IF( KND(IQ) .EQ. 1 )   GO TO 150                                  
      IQ = IQ-1                                                         
      SUMC = SUMC + C(IQ+1)                                             
C                                                                       
      IF( SUMC + CC*(K-IP+IQ) .GT. OEIC + 40.D0 )   GO TO 180           
      GO TO 150                                                         
C                                                                       
  140 IP = IP+1                                                         
      IQ = IP-1                                                         
      SUMC = C(IP)                                                      
      IF( SUMC + CC .GT. OEIC + 40.D0 )   GO TO 200                     
C                                                                       
  150 M = K-IP+IQ                                                       
      SUM = SUMC                                                        
      IF( IQ .EQ. 0 )   GO TO 165                                       
      DO 160   I=1,IQ                                                   
      IF( KND(I) .EQ. 1 )   GO TO 160                                   
      M = M-1                                                           
      SUM = SUM + C(I)                                                  
  160 CONTINUE                                                          
  165 CONTINUE                                                          
      EIC = SUM + CC*M - OEIC                                           
      IF( EIC .GT. 40.D0 )   GO TO 100                                  
      EXIC = DEXP( -0.5D0*EIC )                                         
      ESUM(M+1) = ESUM(M+1) + EXIC                                      
      DO 170   I=1,K                                                    
cxx  170 IF( KND(I) .EQ. 1 )   D(I,M+1) = D(I,M+1) + EXIC                  
      IF( KND(I) .EQ. 1 )   D(I,M+1) = D(I,M+1) + EXIC 
  170 CONTINUE
      GO TO 100                                                         
C         --------------------------------------------                  
  180          DO 190   I=1,IP                                          
cxx  190          KND(I) = 1                                               
               KND(I) = 1 
  190          CONTINUE
               KND(IP+1) = 0                                            
               IP = IP+1                                                
               IQ = IP-1                                                
               GO TO 130                                                
C         ---------------------------------------------                 
C                                                                       
C--------------------------  WHOLE SUBSETS CHECKED  --------------------
C                                                                       
  200 CONTINUE                                                          
cc      IF( IPR .GE. 2 )     WRITE( 6,8 )                                 
cc      IF( IPR .GE. 2 )     WRITE( 6,607 )     (ESUM(I),I=1,K1)          
      DO 201 I=1,K1
cxx  201 ESUM1(I) = ESUM(I)
      ESUM1(I) = ESUM(I)
  201 CONTINUE
C                                                                       
C          MEAN OF NUMBER OF PARAMETERS                                 
C                                                                       
      OSUM = 0.D0                                                       
      SUM = ESUM(1)                                                     
      DO 210   I=1,K                                                    
      SUM = SUM + ESUM(I+1)                                             
cxx  210 OSUM = OSUM + I*ESUM(I+1)                                         
      OSUM = OSUM + I*ESUM(I+1)
  210 CONTINUE
      OMEAN = OSUM / SUM                                                
      OM = OMEAN / K                                                    
cc      IF( IPR .GE. 2 )     WRITE( 6,608 )     OMEAN , OM                
C                                                                       
C       --  BINOMIAL TYPE DAMPER  --                                    
C                                                                       
      DO 220   I=1,K1                                                   
      J = I-1                                                           
      KMJ = K-J                                                         
cc  220 C(I) = BICOEF(K,J)*(OM**J)*((1.D0-OM)**KMJ)                       
cc      C(1) = 1.D0 / (1.D0 + K)                                          
cxx  220 C1(I) = BICOEF(K,J)*(OM**J)*((1.D0-OM)**KMJ)                       
      C1(I) = BICOEF(K,J)*(OM**J)*((1.D0-OM)**KMJ)
  220 CONTINUE
      C1(1) = 1.D0 / (1.D0 + K)                                          
      DO 221  I=1,K                                                     
cc  221 C(I+1) = C(I) * I / (1.D0 + K - I)                                
cxx  221 C1(I+1) = C1(I) * I / (1.D0 + K - I)                                
      C1(I+1) = C1(I) * I / (1.D0 + K - I) 
  221 CONTINUE
cc      IF( IPR .GE. 2 )     WRITE( 6,609 )                               
cc      IF( IPR .GE. 2 )     WRITE( 6,607 )     (C(I),I=1,K1)             
C                                                                       
      DO 230   I=1,K1                                                   
cc  230 ESUM(I) = ESUM(I)*C(I)                                            
cxx  230 ESUM(I) = ESUM(I)*C1(I)                                            
      ESUM(I) = ESUM(I)*C1(I)
  230 CONTINUE
C                                                                       
      SUM = 0.D0                                                        
      DO 240   I=1,K1                                                   
cxx  240 SUM = SUM + ESUM(I)                                               
      SUM = SUM + ESUM(I)
  240 CONTINUE
C                                                                       
cxx      DO 250   J=1,K1                                                   
      DO 251   J=1,K1 
      DO 250   I=1,K                                                    
cc  250 D(I,J) = D(I,J)*C(J) / SUM                                        
cxx  250 D(I,J) = D(I,J)*C1(J) / SUM                                        
      D(I,J) = D(I,J)*C1(J) / SUM
  250 CONTINUE
  251 CONTINUE
C                                                                       
C          WEIGHTS OF PARTIAL CORRELATIONS                              
C                                                                       
cxx      DO 260   I=1,K                                                    
cc  260 C(I) = 0.D0                                                       
cxx  260 C2(I) = 0.D0
      C2(1:K) = 0.D0                                                      
cxx      DO 270   I=1,K                                                    
      DO 271   I=1,K
      DO 270   J=1,K1                                                   
cc  270 C(I) = C(I) + D(I,J)                                              
cxx  270 C2(I) = C2(I) + D(I,J)
      C2(I) = C2(I) + D(I,J)
  270 CONTINUE
  271 CONTINUE                                              
cc      IF( IPR .GE. 2 )     WRITE( 6,603 )                               
cc      IF( IPR .GE. 2 )     WRITE( 6,607 )     (C(I),I=1,K)              
C                                                                       
C          AVERAGING AND REARRANGEMENT OF PARTIAL CORRELATIONS          
C                                                                       
      EK = 1.D0                                                         
      DO 280   I=1,K                                                    
      J = IND(I)                                                        
cc      B(J) = B(J)*C(I)                                                  
cc  280 EK = EK + C(I)**2                                                 
      B(J) = B(J)*C2(I)                                                  
cxx  280 EK = EK + C2(I)**2                                                 
      EK = EK + C2(I)**2
  280 CONTINUE
cc      IF( IPR .LE. 1 )     RETURN                                       
cc      WRITE( 6,602 )                                                    
cc      DO  290     I=1,K                                                 
cc  290 WRITE( 6,609 )     I , B(I)                                       
C                                                                       
      RETURN                                                            
C                                                                       
cxx    3 FORMAT( 1H ,20I3 )                                                
cxx    4 FORMAT( 1H ,3D20.10,2I5 )                                         
cxx    6 FORMAT( 1H ,2I7,F13.5 )                                           
cxx    7 FORMAT( 1H ,6X,'I IND(I)',4X,'N*B(I)**2' )                        
cxx    8 FORMAT( 1H ,'ESUM(I) (I=1,M+1)' )                                 
cxx    9 FORMAT( 1H ,'***  BINOMIAL TYPE  ***' )                           
cxx  602 FORMAT( 1H ,'PARTIAL CORRELATIONS OF THE BAYESIAN MODEL',/,1H ,6X,
cxx     1'I',9X,'A(I)' )                                                   
cxx  603 FORMAT( 1H ,'FINAL BAYESIAN WEIGHTS OF PARTIAL CORRELATIONS' )    
cxx  604 FORMAT( 1H ,'  OAIC =',D13.5 )                                    
cxx  606 FORMAT( 1H ,'OSUM =',D20.10,5X,'SUM =',D20.10,5X,'COD =',D20.10 ) 
cxx  607 FORMAT( 1H ,10D13.5 )                                             
cxx  608 FORMAT( 1H ,'OMEAN =',D15.8,5X,'OM =',D15.8 )                     
cxx  609 FORMAT( 1H ,I7,F13.5 )                                            
      END                                                               
cc      SUBROUTINE  SBBAYS( X,D,K,N,IPR,MJ1,A,SD )                        
      SUBROUTINE  SBBAYS( X,K,N,IPR,MJ1,A,SD,EK,AIC,IND,C,C1,C2,B,
     *OEIC,ESUM,OMEAN,OM  )                        
C                                                                       
C     THIS SUBROUTINE PRODUCES BAYESIAN MODEL BASED ON ALL SUBSET       
C     REGRESSION MODELS USING THE OUTPUT OF SUBROUTINE REDUCT.          
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             RECOEF                                                    
C             SDCOMP                                                    
C             SUBSPC                                                    
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS:                                                         
C          X:     N*(K+1) TRIANGULAR MATRIX,OUTPUT OF SUBROUTINE REDUCT 
C          K:     NUMBER OF REGRESSORS OF THE BAYESIAN MODEL            
C          N:     DATA LENGTH                                           
C          IPR:   =0  TO SUPPRESS THE OUTPUTS                           
C                 =1  TO PRINT OUT FINAL RESULT                         
C                 =2  TO PRINT OUT INTERIM AND FINAL RESULTS            
C          MJ1:   ABSOLUTE DIMENSION OF X                               
C                                                                       
C       OUTPUTS:                                                        
C          A(I) (I=1,K):   REGRESSION COEFFICIENTS OF BAYESIAN MODEL    
C          SD:    RESIDUAL VARIANCE                                     
C                                                                       
cxx      IMPLICIT  REAL * 8  (A-H , O-Z )                                  
cc      DIMENSION  X(MJ1,1) , A(1) , D(1)                                 
cc      DIMENSION  B(50) , G(50)                                          
cc      DIMENSION  IND(50), C(50), C1(50), C2(50), ESUM(50)
cx      DIMENSION  X(MJ1,1) , A(1) , D(K)                                 
cxx      DIMENSION  X(MJ1,K+1) , A(K) , D(K)                                 
cxx      DIMENSION  B(K) , G(K)
cxx      DIMENSION  IND(K), C(K), C1(K+1), C2(K), ESUM(K+1)
      INTEGER :: K, N, IPR, MJ1, IND(K)
      REAL(8) :: X(MJ1,K+1), A(K), SD, EK, AIC, C(K), C1(K+1), C2(K),
     1           B(K), OEIC, ESUM(K+1), OMEAN, OM
      REAL(8) :: D(K), G(K), FN, SUM, BB
      K1 = K + 1                                                        
      FN = N                                                            
cc      IF( IPR .GE. 2 )     WRITE( 6,3 )                                 
C                                                                       
C          PARTIAL CORRELATIONS COMPUTATION                             
C                                                                       
      SUM = X(K1,K1)**2                                                 
      DO 10   I=1,K                                                     
      J = K1-I                                                          
      SUM = SUM + X(J,K1)**2                                            
      G(J) = DSQRT( SUM )                                               
cxx   10 B(J) = X(J,K1)*X(J,J) / (G(J)*DABS(X(J,J)))                       
      B(J) = X(J,K1)*X(J,J) / (G(J)*DABS(X(J,J))) 
   10 CONTINUE
C                                                                       
C          PARTIAL CORRELATIONS OF BAYESIAN MODEL COMPUTATION           
C                                                                       
cc      CALL  SUBSPC( B,K,N,IPR,EK )                                      
cxx      CALL  SUBSPC( B,K,N,IPR,EK,IND,C,C1,C2,OEIC,ESUM,OMEAN,OM )
      CALL  SUBSPC( B,K,N,EK,IND,C,C1,C2,OEIC,ESUM,OMEAN,OM )
C                                                                       
C          MODIFICATION OF CROSS-PRODUCTS  X(I,K1) (I=1,K)              
C                                                                       
      DO 30   I=1,K                                                     
cc      B(I) = B(I)*X(I,I)*G(I) / DABS(X(I,I))                            
      BB = B(I)*X(I,I)*G(I) / DABS(X(I,I))                            
      D(I) = X(I,K1)                                                    
cc   30 X(I,K1) = B(I)                                                    
cxx   30 X(I,K1) = BB
      X(I,K1) = BB
   30 CONTINUE
C                                                                       
C          REGRESSION COEFFICIENTS OF BAYSIAN MODEL                     
C                                                                       
      CALL  RECOEF( X,K,K,MJ1,A )                                       
C                                                                       
      DO 40   I=1,K                                                     
cxx   40 X(I,K1) = D(I)                                                    
      X(I,K1) = D(I)
   40 CONTINUE
C                                                                       
C          RESIDUAL VARIANCE AND AIC                                    
C                                                                       
cc      CALL  SDCOMP( X,A,D,N,K,MJ1,SD )                                  
      CALL  SDCOMP( X,A,N,K,MJ1,SD )                                  
C                                                                       
      IF( IPR .EQ. 0 )     RETURN                                       
      AIC = FN*DLOG(SD) + 2.D0*EK                                       
cc      WRITE( 6,6 )     SD , EK , AIC                                    
      RETURN                                                            
C                                                                       
cxx    3 FORMAT( //1H ,18(1H-),/,' BAYESIAN PROCEDURE',/,1H ,18(1H-) )     
cxx    5 FORMAT( 1H ,10D13.5 )                                             
cxx    6 FORMAT( 1H ,'RESIDUAL VARIANCE',16X,'SD =',D19.8,/,1H ,'EQUIVALENT
cxx     1 NUMBER OF PARAMETERS  EK =',F10.3,/,1H ,32X,'AIC =',F15.3 )      
      E N D                                                             

