      SUBROUTINE BLOCARF( ZS,N,LAG,NS0,KMAX,ZMEAN,SUM,AIC,C,B,A,SD,NP,
     *                    NE,SXX )
C
     
C
cc      PROGRAM  BLOCAR                                                   
C.......................................................................
C.....PLANNED BY H.AKAIKE...............................................
C.....DESIGNED BY H.AKAIKE AND G.KITAGAWA...............................
C.....PROGRAMMED BY G.KITAGAWA AND F.TADA...............................
C.....ADDRESS: THE INSTITUTE OF STATISTICAL MATHEMATICS, 4-6-7 MINAMI-AZ
C..............MINATO-KU, TOKYO 106, JAPAN..............................
C.....DATE OF THE LATEST REVISION:  MAR. 6,1979.........................
C.......................................................................
C.....THIS PROGRAM WAS ORIGINALLY PUBLISHED IN "TIMSAC-78", BY H.AKAIKE,
C.....G.KITAGAWA, E.ARAHATA AND F.TADA, COMPUTER SCIENCE MONOGRAPHS, NO.
C.....THE INSTITUTE OF STATISTICAL MATHEMATICS, TOKYO, 1979.............
C.......................................................................
C     TIMSAC 78.3.2.                                                    
C     _                  ___                __                          
C     BAYESIAN METHOD OF LOCALLY STATIONARY AR MODEL FITTING; SCALAR CAS
C                                                                       
C     THIS PROGRAM LOCALLY FITS AUTOREGRESSIVE MODELS TO NON-STATIONARY 
C     SERIES BY A BAYESIAN PROCEDURE.  POWER SPECTRA FOR STATIONARY SPAN
C     ARE GRAPHICALLY PRINTED OUT.  (THIS PROGRAM IS TENTATIVE.)        
C                                                                       
C     INPUTS REQUIRED:                                                  
C             MT:       INPUT DEVICE FOR ORIGINAL DATA (MT=5 : CARD READ
C             LAG:      UPPER LIMIT OF THE ORDER OF AR MODEL, MUST BE LE
C                       OR EQUAL TO 50.                                 
C             NS:       LENGTH OF BASIC LOCAL SPAN                      
C             KSW:      =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESS
C                       =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REG
C                                                                       
C               -- THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE REDA
C             TITLE:    SPECIFICATION OF DATA                           
C             N:        DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 10000
C             DFORM:    INPUT DATA SPECIFICATION STATEMENT.             
C                       -- EXAMPLE  --     (8F10.5)                     
C             (Z(I),I=1,N):  ORIGINAL DATA                              
C               --------------------------------------------------------
C                                                                       
C                                                                       
C     AT EACH STAGE OF MODELLING OF LOCAL AR MODEL, A TWO-STEP BAYESIAN 
C     PROCEDURE IS APPLIED                                              
C        1.   AVERAGING OF THE MODELS WITH DIFFERENT ORDERS FITTED TO TH
C             NEWLY OBTAINED DATA (BY SUBROUTINE ARBAYS).               
C        2.   AVERAGING OF THE MODELS FITTED TO THE PRESENT AND PRECEDIN
C             (AS FAR AS KMAX(<21)) SPANS.                              
C     AS TO THE AVERAGING 1, SEE THE COMMENTS OF PROGRAM UNIBAR.        
C                                                                       
C     AIC OF THE MODEL FITTED TO THE NEW SPAN IS DEFINED BY             
C                                                                       
C          AIC = NS * LOG( SD ) + 2 * EK,                               
C     WHERE                                                             
C          NS:     LENGTH OF NEW DATA                                   
C          SD:     INNOVATION VARIANCE                                  
C          EK:     EQUIVALENT NUMBER OF PARAMETERS, DEFINED AS THE SUM O
C                  SQUARES OF THE BAYESIAN WEIGHTS                      
C                                                                       
C     AIC OF THE MODELS FITTED TO THE PRECEDING SPANS ARE DEFINED BY    
C                                                                       
C          AIC(J+1) = NS * LOG( SD(J) ) + 2     (J=1,KC)                
C     WHERE                                                             
C          SD(J):  PREDICTION ERROR VARIANCE BY THE MODEL FITTED TO J   
C                  PERIODS FORMER SPAN.                                 
C       --------------------------------------------------------------- 
C       REFERENCES:                                                     
C          G.KITAGAWA AND H.AKAIKE(1978), "A PROCEDURE FOR THE MODELING 
C          OF NON-STATIONARY TIME SERIES.",  ANN. INST. STATIST. MATH., 
C          30,B,351-363.                                                
C                                                                       
C          H.AKAIKE(1978), "A BAYESIAN EXTENSION OF THE MINIMUM AIC     
C          PROCEDURE OF AUTOREGRESSIVE MODEL FITTING.",  RESEARCH MEMO. 
C          NO. 126, THE INSTITUTE OF STATISTICAL MATHEMATICS; TOKYO.    
C       --------------------------------------------------------------- 
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
C             REDATA                                                    
C             NONSTB                                                    
C             PRINTA                                                    
C             NRASPE                                                    
C       --------------------------------------------------------------- 
C                                                                       
cc      !DEC$ ATTRIBUTES DLLEXPORT :: BLOCARF
C
CC      IMPLICIT REAL * 8 ( A-H , O-Y )                                   
cxx      IMPLICIT REAL * 8 ( A-H , O-Z )
cc      REAL * 4  TITLE(20) , TTL(13)                                     
cc      DIMENSION Z(10000)                                                
cc      DIMENSION  X(200,51) , D(200) , A(50)                             
cxx      DIMENSION  Z(N), ZS(N)                                            
cxx      DIMENSION  X(NS0,LAG+1), A(LAG,KMAX)
cxx      DIMENSION  AIC(KMAX,KMAX), C(KMAX,KMAX)
cxx      DIMENSION  B(LAG,KMAX), SD(KMAX)
cxx      DIMENSION  NP(KMAX), NE(KMAX)
      INTEGER :: N, LAG, NS0, KMAX, NP(KMAX), NE(KMAX)
      REAL(8) :: ZS(N), ZMEAN, SUM, AIC(KMAX,KMAX), C(KMAX,KMAX),
     1           B(LAG,KMAX), A(LAG,KMAX), SD(KMAX), SXX(121,KMAX)
cc      DATA  TTL / 4H  CU,4HRREN,4HT MO,4HDEL ,4H(AVE,4HRAGE,4H BY ,4HTHE
cc     1 ,4HBAYE,4HSIAN,4H WEI,4HGHTS,4H)    /                            
C
cxx      DIMENSION    SXX(121,KMAX)
cxx      DIMENSION    F(LAG,KMAX)
      REAL(8) :: Z(N), X(NS0,LAG+1), F(LAG,KMAX), BB
C
      EXTERNAL  SETX1
C
cc      CHARACTER(100) IFLNAM,OFLNAM
cc      CALL FLNAM2( IFLNAM,OFLNAM,NFL )
cc      IF (NFL.EQ.0) GO TO 999
cc      IF (NFL.EQ.2) THEN
cc         OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR )
cc      ELSE
cc         CALL SETWND
cc      END IF
C
C
      KSW = 0                                                           
cc      KMAX = 20                                                         
cc      MJ1 = 200                                                         
      MJ1 = NS0
      ISW = 0                                                           
C                                                                       
CC      READ( 5,1 )     MT                                                
cc      MT = 5
cc      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
cc      READ( 5,1 )     LAG , NS                                          
      NS = NS0
C                                                                       
cc      WRITE( 6,3 )                                                      
cc      IF( KSW .EQ. 1 )     WRITE( 6,5 )                                 
cc      IF( KSW .NE. 1 )     WRITE( 6,4 )                                 
cc      WRITE( 6,6 )                                                      
cc      WRITE( 6,2 )   LAG , NS , MT                                      
C                                                                       
C          ---------------------                                        
C          ORIGINAL DATA LOADING                                        
C          ---------------------                                        
C                                                                       
cc      CALL  REDATA( Z,N,MT,TITLE )                                      
cc      SUBROUTINE  REDATA( X,N,MT,TITLE )                                
      CALL  REDATA( ZS,Z,N,ZMEAN,SUM )
cc      CLOSE( MT )
C                                                                       
C     ---  INITIAL SETTING  ---                                         
      LAG1 = LAG+1                                                      
      N0 = 0                                                            
      K = LAG + KSW                                                     
      K3 = K*3                                                          
C                                                                       
      NR = 0
      KC = 0
   10 CONTINUE                                                          
      NR = NR+1
C                                                                       
C          -----------------------------------                          
C          LOCALLY STATIONARY AR-MODEL FITTING                          
C          -----------------------------------                          
C                                                                       
cc      CALL  NONSTB( SETX1,Z,X,D,LAG,N0,NS,KMAX,KSW,ISW,TITLE,MJ1,A,SD ) 
cx      CALL  NONSTB( SETX1,Z,X,LAG,N0,NS,KMAX,KSW,ISW,MJ1,KC,F,
      CALL  NONSTB( SETX1,Z,N,X,LAG,N0,NS,KMAX,KSW,ISW,MJ1,KC,F,
     *              AIC(1,NR),C(1,NR),B(1,NR),A(1,NR),SD(NR) )
C                                                                       
cc      NP = N0 + LAG + 1                                                 
cc      NE = NP+NS-1                                                      
cc      CALL  PRINTA( A,SD,K,TTL,13,TITLE,NP,NE )                         
      NP(NR) = N0 + LAG + 1                                             
      NE(NR) = NP(NR)+NS-1                                              
C                                                                       
C     ---  SPECTRUM DISPLAY  ---                                        
cc      CALL  NRASPE( SD,A,B,K,0,120,TITLE )                              
      CALL  NRASPE( SD(NR),A(1,NR),BB,K,0,120,SXX(1,NR) )
C                                                                       
C                                                                       
C     ---  PREPARATION FOR THE NEXT STEP  ---                           
      N0 = N0 + NS                                                      
      IF(N0+NS+LAG1.GT.N)      NS = N-N0-LAG1                           
      IF(N-N0-NS-LAG1.LT.K3)   NS = N-N0-LAG1                           
      IF( N0+LAG1 .LT. N )     GO TO 10                                 
C                                                                       
C                                                                       
cc      GO TO 999
C
cc  900 CONTINUE
cc      WRITE(6,600) IVAR,OFLNAM
cc      GO TO 999
C
cc  910 CONTINUE
cc      IF (NFL.EQ.2) CLOSE( 6 )
cc#ifdef __linux__
ccC	reopen #6 as stdout
cc      IF (NFL.EQ.2) OPEN(6, FILE='/dev/fd/1')
cc#endif
ccC /* __linux__ */
cc      WRITE(6,610) IVAR,IFLNAM
cc      GO TO 999
C
cxx  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,//5X,100A)
cxx  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//5X,100A)
C
      RETURN
cxx    1 FORMAT( 16I5 )                                                    
cxx    2 FORMAT( ///1H ,'  FITTING UP TO THE ORDER  K =',I3,'  IS TRIED',/,
cxx     1'   BASIC LOCAL SPAN  NS =',I4,/,'   ORIGINAL DATA INPUT DEVICE  M
cxx     2T =',I3 )                                                         
cxx    3 FORMAT( //' PROGRAM TIMSAC 78.3.2',/,'   BAYESIAN METHOD OF LOCALL
cxx     *Y STATIONARY AR MODEL FITTING;   SCALAR CASE',//,'   < BASIC AUTOR
cxx     *EGRESSIVE MODEL >' )                                              
cxx    4 FORMAT( 1H ,10X,'Z(I) = A(1)*Z(I-1) + A(2)*Z(I-2) + ... + A(M)*Z(I
cxx     1-M) + E(I)' )                                                     
cxx    5 FORMAT( 1H ,10X,'Z(I) = A(1) + A(2)*Z(I-1) + ... + A(M+1)*Z(I-M) +
cxx     1 E(I)' )                                                          
cxx    6 FORMAT( 1H ,2X,'WHERE',/,11X,'M:     ORDER OF THE MODEL',/,11X,'E(
cxx     1I):  GAUSSIAN WHITE NOISE WITH MEAN 0  AND  VARIANCE SD(M).' )    
      END                                                               
cc      SUBROUTINE  NONSTB( SETX,Z,X,D,LAG,N0,NS,KMAX,KSW,ISW,TITLE,MJ1,A,
cc     1SD )                                                              
cx      SUBROUTINE  NONSTB( SETX,Z,X,LAG,N0,NS,KMAX1,KSW,ISW,MJ1,KC,F,
      SUBROUTINE  NONSTB( SETX,Z,N,X,LAG,N0,NS,KMAX1,KSW,ISW,MJ1,KC,F,
     1AIC,C,B,A,SD )                                                    
C                                                                       
C     BAYESIAN TYPE NON-STATIONARY AUTOREGRESSIVE MODEL FITTING PROCEDUR
C                                                                       
C     THIS SUBROUTINE FIRST FITS AN AUTOREGRESSIVE MODEL TO THE NEWLY OB
C     DATA SPAN AND THEN PRODUCES A BAYESIAN MODEL BY AVERAGING THE MODE
C     FITTED TO THE PRESENT AND PRECEDING SPANS.                        
C                                                                       
C     BAYESIAN WEIGHT OF THE MODEL FITTED TO I PERIODS FORMER SPAN IS   
C     DEFINED BY                                                        
C          C(I)  =  CONST * P(I) / (I+1)                                
C     WHERE                                                             
C          CONST  =  NORMALIZING CONSTANT                               
C          P(I)  =  EXP( -0.5*AIC(I) )                                  
C          AIC(I)  =  AIC WITH RESPECT TO THE PRESENT DATA.             
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             DMIN                                                      
C             ARBAYS                                                    
C             ARCOEF                                                    
C             BAYSWT                                                    
C             REDUCT                                                    
C             SDCOMP                                                    
C       ----------------------------------------------------------------
C                                                                       
C     INPUTS:                                                           
C        SETX:  EXTERNAL SUBROUTINE DESIGNATION                         
C        Z:     ORIGINAL DATA, OUTPUT OF SUBROUTINE REDATA              
C        X:     WORKING AREA                                            
C        D:     WORKING AREA                                            
C        LAG:   UPPER LIMIT OF ORDER OF AR-MODEL                        
C        N0:    INDEX OF THE END POINT OF THE FORMER SPAN               
C        NS:    LENGTH OF BASIC LOCAL SPAN                              
C        KMAX:  MAXIMUM NUMBER OF PRECEDING MODELS STORED               
C        KSW:   =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR      
C               =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSOR  
C        ISW:   PRINT OUT CONTROL                                       
C        TITLE: TITLE OF DATA                                           
C        MJ1:   ABSOLUTE DIMENSION OF X                                 
C                                                                       
C     OUTPUTS:                                                          
C        A:     AR-COEFFICIENTS OF THE CURRENT MODEL                    
C        SD:    INNOVATION VARIANCE OF THE CURRENT MODEL                
C                                                                       
cxx      IMPLICIT  REAL*8  ( A-H,O-Z )                                     
CC      REAL*4  Z , TITLE                                                 
cc      REAL*4  TITLE
cc      DIMENSION  X(MJ1,1) , D(1) , A(1) , Z(1) , TITLE(1)               
cc      DIMENSION  F(50,20) , AIC(21) , C(50) , B(50)                     
cx      DIMENSION  X(MJ1,1) , D(LAG+KSW+1), A(1) , Z(1)
cxx      DIMENSION  X(MJ1,1) , D(LAG+KSW+1), A(LAG+KSW) , Z(N)
cxx      DIMENSION  F(LAG+KSW,KMAX1) , AIC(KMAX1) , C(KMAX1) , B(LAG+KSW)
cxx      DIMENSION  SDD(LAG+KSW+1), AICC(LAG+KSW+1), DIC(LAG+KSW+1)
cxx      DIMENSION  B1(LAG+KSW), W(LAG+KSW+1)
      INTEGER :: N, LAG, N0, NS, KMAX1, KSW, ISW, MJ1, KC
      REAL(8) :: Z(N), X(MJ1,1), F(LAG+KSW,KMAX1), AIC(KMAX1),
     1           C(KMAX1), B(LAG+KSW), A(LAG+KSW), SD
      REAL(8) :: D(LAG+KSW+1), SDD(LAG+KSW+1), AICC(LAG+KSW+1),
     1           DIC(LAG+KSW+1), B1(LAG+KSW), W(LAG+KSW+1),
     2           AICB, AICM, SDMIN, PN
cc      DATA  KC / 0 /                                                    
      EXTERNAL  SETX
C                                                                       
      K = LAG + KSW                                                     
      KMAX = KMAX1-1
C                                                                       
C     ---  HOUSEHOLDER REDUCTION  ---                                   
cc      CALL  REDUCT( SETX,Z,D,NS,N0,K,MJ1,LAG,X )                        
      CALL  REDUCT( SETX,Z,NS,N0,K,MJ1,LAG,X )                        
C                                                                       
C                                                                       
C     ---  BAYESIAN MODEL FITTED TO THE NEW SPAN  ---                   
cc      CALL  ARBAYS( X,D,K,LAG,NS,ISW,TITLE,MJ1,A,B,SD,AICB )            
cxx      CALL ARBAYS( X,D,K,LAG,NS,ISW,MJ1,SDD,AICC,DIC,AICM,SDMIN,IMIN,
      CALL ARBAYS( X,D,K,NS,ISW,MJ1,SDD,AICC,DIC,AICM,SDMIN,IMIN,
     *             A,B1,B,W,SD,PN,AICB )
C                                                                       
      IF( KC .EQ. 0 )  GO TO 110                                        
C                                                                       
C                                                                       
C     ---  PREDICTION ERROR VARIANCE AND AIC OF THE FORMER MODELS  ---  
      AIC(1) = AICB                                                     
      DO 30  J=1,KC                                                     
         DO 20  I=1,K                                                   
cxx   20    D(I) = F(I,J)                                                  
         D(I) = F(I,J)
   20    CONTINUE
        CALL  ARCOEF( D,K,A )                                           
cc        CALL  SDCOMP( X,A,D,NS,K,MJ1,SD )                               
        CALL  SDCOMP( X,A,NS,K,MJ1,SD )                               
cxx   30 AIC(J+1) = NS*DLOG( SD ) + 2.D0                                   
      AIC(J+1) = NS*DLOG( SD ) + 2.D0
   30 CONTINUE
C                                                                       
C                                                                       
C     ---  BAYESIAN WEIGHTS OF THE MODEL  ---                           
c-------------------------------   06/11/01
ccx      AICM = DMIN( AIC,KC )                                             
      AICM = AIC(1)
      DO 33  I=1,KC
cxx   33 IF( AIC(I) .LT. AICM )  AICM = AIC(I)
      IF( AIC(I) .LT. AICM )  AICM = AIC(I)
   33 CONTINUE
c-------------------------------
      CALL  BAYSWT( AIC,AICM,KC,2,C )                                   
C                                                                       
      KC1 = KC+1                                                        
cc      WRITE( 6,3 )     C(1) , AIC(1)                                    
cc      DO 35  I=2,KC1                                                    
cc         IM1 = I-1                                                      
cc   35 WRITE( 6,4 )     IM1 , C(I) , AIC(I)                              
C                                                                       
C                                                                       
C     ---  AVERAGING OF THE MODELS  ---                                 
      DO 40  I=1,K                                                      
cxx   40 B(I) = B(I)*C(1)                                                  
      B(I) = B(I)*C(1) 
   40 CONTINUE
      DO 70  J=1,KC                                                     
         DO 50  I=1,K                                                   
cxx   50    A(I) = F(I,J)                                                  
         A(I) = F(I,J)
   50    CONTINUE
         DO 60  I=1,K                                                   
cxx   60    B(I) = B(I) + C(J+1)*A(I) 
         B(I) = B(I) + C(J+1)*A(I)
   60    CONTINUE                                     
   70 CONTINUE                                                          
cc      WRITE( 6,5 )     (B(I),I=1,K)                                     
C                                                                       
C                                                                       
C     ---  AR-COEFFICIENTS OF THE CURRENT BAYESIAN MODEL  ---           
      CALL  ARCOEF( B,K,A )                                             
C                                                                       
C                                                                       
C     ---  "PARCOR'S" STORED  ---                                       
cxx      DO 100  J=1,KC                                                    
      DO 101  J=1,KC
         II = KC1-J                                                     
         DO 100  I=1,K                                                  
cxx  100 F(I,II+1) = F(I,II)
         F(I,II+1) = F(I,II)
  100    CONTINUE
  101 CONTINUE                                               
  110 CONTINUE                                                          
      DO 120  I=1,K                                                     
cxx  120 F(I,1) = B(I)                                                     
      F(I,1) = B(I)
  120 CONTINUE
      KC = MIN0( KC+1,KMAX )                                            
C                                                                       
C                                                                       
C     ---  PREDICTION ERROR VARIANCE COMPUTED  ---                      
cc      CALL  SDCOMP( X,A,D,NS,K,MJ1,SD )                                 
      CALL  SDCOMP( X,A,NS,K,MJ1,SD )                                 
C                                                                       
      RETURN                                                            
C                                                                       
cxx    3 FORMAT( ///1H ,13X,'AR-MODEL FITTED TO  !  BAYESIAN WEIGHTS  ! AIC
cxx     1 WITH RESPECT TO THE PRESENT DATA',/,10X,83(1H-),/,1H ,11X,'CURREN
cxx     2T BLOCK',9X,'!',F13.5,7X,'!',F21.3 )                              
cxx    4 FORMAT( 1H ,6X,I5,' PERIOD FORMER BLOCK  !',F13.5,7X,'!',F21.3 )  
cxx    5 FORMAT( //1H ,'PARTIAL AUTOCORRELATION  B(I) (I=1,K)',/,(1X,10D13.
cxx     15) )                                                              
C                                                                       
      END                                                               
