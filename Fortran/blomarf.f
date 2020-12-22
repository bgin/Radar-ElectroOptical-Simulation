      SUBROUTINE  BLOMARF( ZS,N,ID,C,LAG,NS0,KMAX,ZMEAN,ZVARI,BW,AIC,A,
cxx     *                     E,AICB,LKS,LKE )
     *                     E,AICB,LKS,LKE,M )
C
    
C
cc      PROGRAM  BLOMAR                                                   
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
C     TIMSAC 78.3.4.                                                    
C     _                  __                 _            __             
C     BAYESIAN METHOD OF LOCALLY STATIONARY MULTIVARIATE AR MODEL FITTIN
C                                                                       
C     THIS PROGRAM LOCALLY FITS MULTI-VARIATE AUTOREGRESSIVE MODELS TO  
C     NON-STATIONARY TIME SERIES BY A BAYESIAN PROCEDURE.               
C                                                                       
C       --------------------------------------------------------------- 
C       REFERENCES:                                                     
C          G.KITAGAWA AND H.AKAIKE(1978), "A PROCEDURE FOR THE MODELING 
C          OF NON-STATIONARY TIME SERIES.",  ANN. INST. STATIST. MATH., 
C          30,B,351-363.                                                
C          H.AKAIKE(1978), "A BAYESIAN EXTENSION OF THE MINIMUM AIC     
C          PROCEDURE OF AUTOREGRESSIVE MODEL FITTING.",  RESEARCH MEMO. 
C          NO. 126, THE INSTITUTE OF STATISTICAL MATHEMATICS; TOKYO.    
C       --------------------------------------------------------------- 
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
C             MRDATA                                                    
C             MNONSB                                                    
C       --------------------------------------------------------------- 
C       INPUTS REQUIRED;                                                
C          MT:    INPUT DEVICE FOR ORIGINAL DATA (MT=5: CARD READER).   
C          LAG:   UPPER LIMIT OF THE ORDER OF AR-MODEL, MUST BE LESS THA
C                 OR EQUAL TO 50.                                       
C          NS:    LENGTH OF BASIC LOCAL SPAN.                           
C          KSW:   =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR    
C                 =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSOR
C                                                                       
C            -- THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE MRDATA 
C          TITLE: SPECIFICATION OF DATA                                 
C          N:     DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 1000.      
C          ID:    DIMENSION OF DATA,  MUST BE LESS THAN 6               
C                       < ID*(LAG+1)+KSW MUST BE LESS THAN 101 >        
C          IFM:   INPUT FORMAT                                          
C          FORM:  INPUT DATA FORMAT SPECIFICATION STATEMENT.            
C                 -- EXAMPLE --     (8F10.5)                            
C          C(J):  CALIBRATION CONSTANT FOR CHANNEL J (J=1,ID)           
C          Z(I,J): ORIGINAL DATA                                        
C            -----------------------------------------------------------
C                                                                       
cc      !DEC$ ATTRIBUTES DLLEXPORT :: BLOMARF
C
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL * 4  Z                                                       
cc      DIMENSION  Z(1500,10)                                             
cc      DIMENSION  X(200,100) , D(200)                                    
cc      DIMENSION  A(10,10,50) , B(10,10,50)                              
cc      DIMENSION  G(10,10,50) , H(10,10,50) , E(10,10)                   
cxx      DIMENSION  ZS(N,ID), Z(N,ID), C(ID)
cxx      DIMENSION  ZMEAN(ID), ZVARI(ID)
cxx      DIMENSION  A(ID,ID,LAG,KMAX) , B(ID,ID,LAG)
cxx      DIMENSION  G(ID,ID,LAG) , H(ID,ID,LAG) , E(ID,ID,KMAX)
cxx      DIMENSION  BW(KMAX,KMAX), AIC(KMAX,KMAX)
cxx      DIMENSION  AICB(KMAX), LKS(KMAX), LKE(KMAX)
cxx      DIMENSION  F1(LAG*ID,ID,KMAX), F2(LAG*ID,ID,KMAX)
cxxcxx      DIMENSION  X(NS0+LAG+1,(LAG+1)*ID)
cxx      DIMENSION  X((LAG+1)*ID*2,(LAG+1)*ID)
      INTEGER :: N, ID, LAG, NS0, KMAX, LKS(KMAX), LKE(KMAX), M
      REAL(8) :: ZS(N,ID), C(ID), ZMEAN(ID), ZVARI(ID),
     1           BW(KMAX,KMAX), AIC(KMAX,KMAX), A(ID,ID,LAG,KMAX),
     2           E(ID,ID,KMAX), AICB(KMAX)
      REAL(8) :: Z(N,ID), B(ID,ID,LAG), G(ID,ID,LAG), H(ID,ID,LAG),
     1           F1(LAG*ID,ID,KMAX), F2(LAG*ID,ID,KMAX),
     2           X((LAG+1)*ID*2,(LAG+1)*ID)
C
C       PARAMETERS:                                                     
C          MJ:    ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
C          MJ1:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
C          MJ2:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
C          MJ3:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
C                                                                       
cc      CHARACTER(100) IFLNAM,OFLNAM
cc      CALL FLNAM2( IFLNAM,OFLNAM,NFL )
cc      IF (NFL.EQ.0) GO TO 999
cc      IF (NFL.EQ.2) THEN
cc         OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR)
cc      ELSE
cc         CALL SETWND
cc      END IF
C
cc      MJ = 1500                                                         
cc      MJ1 = 200                                                         
cc      MJ3 = 10                                                          
cc      KMAX = 10                                                         
      MJ = N
cxx      MJ1 = NS0+LAG+1
      MJ1 = (LAG+1)*ID*2
      MJ3 = ID
      KSW = 0                                                           
C
      BW(1:KMAX,1:KMAX) = 0.0D0
      AIC(1:KMAX,1:KMAX) = 0.0D0
      A(1:ID,1:ID,1:LAG,1:KMAX) = 0.0D0
      E(1:ID,1:ID,1:KMAX) = 0.0D0
      AICB(1:KMAX) = 0.0D0
      LKS(1:KMAX) = 0
      LKE(1:KMAX) = 0
      F1(1:LAG*ID,1:ID,1:KMAX) = 0.0D0
      F2(1:LAG*ID,1:ID,1:KMAX) = 0.0D0
      X(1:MJ1,1:(LAG+1)*ID) = 0.0D0
C
CC      READ( 5,1 )     MT                                                
cc      MT = 5
cc      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
cc      READ( 5,1 )     LAG , NS                                          
      NS = NS0
C                                                                       
cc      WRITE( 6,2 )                                                      
cc      WRITE( 6,4 )                                                      
cc      WRITE( 6,3 )     LAG , NS , MT                                    
C                                                                       
cc      CALL  MRDATA( MT,MJ,Z,N,ID )                                      
      CALL MRDATA( ZS,Z,N,ID,C,ZMEAN,ZVARI )
cc      CLOSE( MT )
C                                                                       
      L = 0                                                             
      KD = LAG * ID + KSW                                               
      MX = KD * 2                                                       
      KC = 0
C                                                                       
C                                                                       
      M = 0
  111 CONTINUE                                                          
cxx      M = M+1
      LK = L + LAG                                                      
      LK1 = LK + 1                                                      
      IF( LK1 .GE. N )     GO TO 300
      M = M+1
      IF( N-LK1 .LE. NS )     NS = N - LK                               
      IF( N-LK1-NS .LT. MX )     NS = N - LK                            
C                                                                       
cc      CALL MNONSB( Z,X,D,G,H,E,KSW,LAG,L,NS,ID,KMAX,MJ,MJ1,MJ3,A,B,AIC )
cxx      CALL MNONSB( Z,G,H,E(1,1,M),KSW,LAG,L,NS,ID,KMAX,KC,MJ,MJ1,MJ3,
cxx     *             BW(1,M),AIC(1,M),A(1,1,1,M),B,AICB(M),F1,F2 )
      CALL MNONSB( Z,X,G,H,E(1,1,M),KSW,LAG,L,NS,ID,KMAX,KC,MJ,MJ1,
     *             MJ3,BW(1,M),AIC(1,M),A(1,1,1,M),B,AICB(M),F1,F2 )
C                                                                       
      L = L + NS                                                        
C                                                                       
cc      LKE = LK + NS                                                     
      LKE(M) = LK + NS
      MF = LAG                                                          
cc      WRITE( 6,13 )                                                     
cc      WRITE( 6,16 )                                                     
cc      WRITE( 6,14 )     LK1 , LKE                                       
      LKS(M) = LK1
cc      DO 10  I=1,MF                                                     
cc      WRITE( 6,16 )                                                     
cc      DO 10  II=1,ID                                                    
cc      WRITE( 6,15 )     (A(II,JJ,I),JJ=1,ID)                            
cc      IF( II .EQ. 1 )     WRITE( 6,17 )   I                             
cc      IF( II .NE. 1 )     WRITE( 6,21 )                                 
cc   10 CONTINUE                                                          
cc      WRITE( 6,16 )                                                     
cc      WRITE( 6,19 )     MF , AIC                                        
cc      WRITE( 6,16 )                                                     
cc      WRITE( 6,12 )                                                     
cc      DO 20  I=1,ID                                                     
cc      WRITE( 6,11 )     (E(I,J),J=1,ID)                                 
cc   20 WRITE( 6,17 )     I                                               
cc      WRITE( 6,16 )                                                     
cc      WRITE( 6,18 )                                                     
C                                                                       
      GO TO 111                                                         
  300 CONTINUE
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
ccWRITE(6,610) IVAR,IFLNAM
ccC
cxx  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,//,5X,100A)
cxx  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//,5X,100A)
C
cc  999 CONTINUE
      RETURN
cxx    1 FORMAT( 16I5 )                                                    
cxx    2 FORMAT( 1H ,'PROGRAM TIMSAC 78.3.4',//,'   BAYESIAN METHOD OF LOCA
cxx     *LLY STATIONARY MULTIVARIATE AR MODEL FITTING;',//,'  < BASIC AUTOR
cxx     *EGRESSIVE MODEL >' )                                              
cxx    3 FORMAT( //,1H ,'  UPPER LIMIT OF THE ORDER  K =',I3,/,'   BASIC LO
cxx     1CAL SPAN LENGTH  NS =',I4,/,'   ORIGINAL DATA INPUT DEVICE  MT =',
cxx     2I3 )                                                              
cxx    4 FORMAT( //1H ,10X,'Z(N) = A1*Z(N-1) + A2*Z(N-2) + ... + AK*Z(N-K) 
cxx     1+ W(N)',/,1H ,'  WHERE',/,11X,'K:     ORDER OF THE MODEL',/,11X,  
cxx     2'W(N):  INNOVATION' )                                             
cxx   11 FORMAT( 1H ,18X,5D15.5 )                                          
cxx   12 FORMAT( 1H ,10X,1H.,6X,'INNOVATION VARIANCE MATRIX',53X,1H. )     
cxx   13 FORMAT( 1H ,//11X,35(1H.),'  CURRENT MODEL  ',35(1H.) )           
cxx   14 FORMAT( 1H ,10X,1H.,6X,'M',7X,'AM(I,J)',30X,'DATA  Z(K,.); K=',I5,
cxx     11H,,I5,7X,1H. )                                                   
cxx   15 FORMAT( 1H ,18X,5F15.8 )                                          
cxx   16 FORMAT( 1H ,10X,1H.,85X,1H. )                                     
cxx   17 FORMAT( 1H+,10X,1H.,I7,78X,1H. )                                  
cxx   18 FORMAT( 1H ,10X,87(1H.) )                                         
cxx   19 FORMAT( 1H ,10X,1H.,6X,'ORDER =',I5,67X,1H.,/,11X,1H.,6X,'AIC =', 
cxx     1 F15.3,59X,1H. )                                                  
cxx   21 FORMAT( 1H+,10X,1H.,85X,1H. )                                     
      END                                                               
cc      SUBROUTINE  MNONSB( Z,X,D,G,H,E,KSW,LAG,N0,NS,ID,KMAX,MJ,MJ1,MJ3,A
cc     *,B,AICB )                                                         
cxx      SUBROUTINE  MNONSB( Z,G,H,E,KSW,LAG,N0,NS,ID,KMAX1,KC,MJ,MJ1,MJ3,
cxx     *                    C,AIC,A,B,AICB,F1,F2 )
      SUBROUTINE  MNONSB( Z,X,G,H,E,KSW,LAG,N0,NS,ID,KMAX1,KC,MJ,MJ1,
     *                    MJ3,C,AIC,A,B,AICB,F1,F2 )
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             DMIN                                                      
C             BAYSWT                                                    
C             MARCOF                                                    
C             MBYSAR                                                    
C             MREDCT                                                    
C             MSDCOM                                                    
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS:                                                         
C          Z:     ORIGINAL DATA; Z(K,I) (K=1,N) REPRESENTS THE RECORD OF
C                 THE I-TH CHANNEL                                      
C          X:     WORKING AREA                                          
C          D:     WORKING AREA                                          
C          G:     WORKING AREA (PARTIAL AUTOREGRESSION COEFFICIENT MATRI
C                 FORWARD MODEL)                                        
C          H:     WORKING AREA (PARTIAL AUTOREGRESSION COEFFICIENT MATRI
C                 BACKWARD MODEL)                                       
C          E:     WORKING AREA                                          
C          KSW:   =0   CONSTATNT VECTOR IS NOT INCLUDED AS A REGRESSOR  
C                 =1   CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSO
C          LAG:   UPPER LIMIT OF THE ORDER OF AR-MODEL                  
C          N0:    INDEX OF THE END POINT OF THE FORMER SPAN             
C          NS:    LENGTH OF BASIC LOCAL SPAN                            
C          ID:    DIMENSION OF DATA                                     
C          KMAX:  MAXIMUM NUMBER OF PRECEDING MODELS STORED             
C          MJ:    ABSOLUTE DIMENSION OF Z IN THE MAIN PROGRAM           
C          MJ1:   ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM           
C          MJ3:   ABSOLUTE DIMENSION OF A AND B IN THE MAIN PROGRAM     
C                                                                       
C       OUTPUTS:                                                        
C          A:     AUTOREGRESSIVE COEFFICIENT MATRIX OF FORWARD MODEL    
C          B:     AUTOREGRESSIVE COEFFICIENT MATRIX OF BACKWARD MODEL   
C          AICB:  AIC OF THE CURRENT MODEL                              
C                                                                       
cxx      IMPLICIT  REAL *8  ( A-H , O-Z )                                  
CC      REAL*4     Z(MJ,1) , F1(100,10,11) , F2(100,10,11)
cc      DIMENSION  X(MJ1,1) , D(1) , A(MJ3,MJ3,1) , B(MJ3,MJ3,1)          
cx      DIMENSION  Z(MJ,1) , F1(LAG*ID,ID,KMAX1) , F2(LAG*ID,ID,KMAX1)
cx      DIMENSION  X(MJ1,(LAG+1)*ID), A(MJ3,MJ3,1) , B(MJ3,MJ3,1)
cx      DIMENSION  G(MJ3,MJ3,1) , H(MJ3,MJ3,1) , E(MJ3,1)                 
cc      DIMENSION  Y(100,10) , AIC(11) , C(11)                            
cxx      DIMENSION  Z(MJ,ID) , F1(LAG*ID,ID,KMAX1) , F2(LAG*ID,ID,KMAX1)
cxx      DIMENSION  X(MJ1,(LAG+1)*ID), A(MJ3,MJ3,LAG) , B(MJ3,MJ3,LAG)
cxx      DIMENSION  G(MJ3,MJ3,LAG) , H(MJ3,MJ3,LAG) , E(MJ3,ID) 
cxx      DIMENSION  AIC(KMAX1) , C(KMAX1)
      INTEGER :: KSW, LAG, N0, NS, ID, KMAX1, KC, MJ, MJ1, MJ3
      REAL(8) :: Z(MJ,ID), X(MJ1,(LAG+1)*ID), G(MJ3,MJ3,LAG),
     1           H(MJ3,MJ3,LAG), E(MJ3,ID), C(KMAX1), AIC(KMAX1),
     2           A(MJ3,MJ3,LAG), B(MJ3,MJ3,LAG), AICB,
     3           F1(LAG*ID,ID,KMAX1), F2(LAG*ID,ID,KMAX1)
cc      DATA     KC  / 0 /                                                
C
cxx      DIMENSION  SD1(LAG+1), AIC1(LAG+1), DIC1(LAG+1)
cxx      DIMENSION  BW1(LAG+1), BW2(LAG)
      REAL(8) :: SD1(LAG+1), AIC1(LAG+1), DIC1(LAG+1),
     1           BW1(LAG+1), BW2(LAG), AICM, EK, SD, SDMIN
C                                                                       
      KMAX = KMAX1-1
cc      MJ5 = 100                                                         
      IPR = 0                                                           
      KD = LAG * ID                                                     
C          -----------------------------------------------              
C          NEW DATA LOADING AND HOUSEHOLDER TRANSFORMATION              
C          -----------------------------------------------              
cc      CALL  MREDCT( Z,D,NS,N0,LAG,ID,MJ,MJ1,KSW,X )        
      CALL  MREDCT( Z,NS,N0,LAG,ID,MJ,MJ1,KSW,X )
C                                                                       
C          -------------------------------------                        
C          BAYESIAN MODEL FITTED TO THE NEW SPAN                        
C          -------------------------------------                        
cc      CALL  MBYSAR( X,D,NS,LAG,ID,KSW,IPR,MJ1,MJ3,A,B,G,H,E,AICB,EK )   
cxx      CALL  MBYSAR( X,NS,LAG,ID,KSW,IPR,MJ1,MJ3,SD1,AIC1,DIC1,
      CALL  MBYSAR( X,NS,LAG,ID,KSW,MJ1,MJ3,SD1,AIC1,DIC1,
     *              AICM,SDMIN,IMIN,BW1,BW2,A,B,G,H,E,AICB,EK )
C                                                                       
      IF( KC .EQ. 0 )  GO TO 20                                         
C          -----------------------------                                
C          "PARCOR'S" SHIFTED AND STORED                                
C          -----------------------------                                
      KC1 = KC+1                                                        
cxx      DO 10  JJ=1,KC
      DO 12  JJ=1,KC
        II = KC1 - JJ                                                   
cxx        DO 10  I=1,KD                                                   
        DO 11  I=1,KD                                                   
        DO 10  J=1,ID                                                   
      F1(I,J,II+1) = F1(I,J,II)                                         
cxx   10 F2(I,J,II+1) = F2(I,J,II)
      F2(I,J,II+1) = F2(I,J,II)
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE                                         
   20 IM = 0                                                            
cxx      DO 30  II=1,LAG                                                   
cxx      DO 30  I=1,ID
      DO 32  II=1,LAG                                                   
      DO 31  I=1,ID                                                      
      IM = IM+1                                                         
      DO 30  J=1,ID                                                     
      F1(IM,J,1) = G(I,J,II)
cxx   30 F2(IM,J,1) = H(I,J,II)
      F2(IM,J,1) = H(I,J,II)
   30 CONTINUE
   31 CONTINUE
   32 CONTINUE                                            
      IF( KC .EQ. 0 )  GO TO 100                                        
C          ---------------------------------------------------------    
C          PREDICTION ERROR VARIANCES AND AIC'S OF THE FORMER MODELS    
C          ---------------------------------------------------------    
      AIC(1) = AICB                                                     
      DO 50  JJ=1,KC                                                    
        IM = 0                                                          
cxx        DO 40  II=1,LAG                                                 
cxx        DO 40  I=1,ID                                                   
        DO 42  II=1,LAG                                                 
        DO 41  I=1,ID
          IM = IM+1                                                     
          DO 40  J=1,ID                                                 
        G(I,J,II) = F1(IM,J,JJ+1)                                       
cxx   40   H(I,J,II) = F2(IM,J,JJ+1)                                       
        H(I,J,II) = F2(IM,J,JJ+1)
   40     CONTINUE
   41   CONTINUE
   42   CONTINUE 
C                                                                       
        CALL  MARCOF( G,H,ID,LAG,MJ3,A,B )                              
cc      CALL  MSDCOM( X,A,Y,D,NS,LAG,ID,KSW,IPR,MJ1,MJ3,MJ5,E,SD )        
cxx      CALL  MSDCOM( X,A,NS,LAG,ID,KSW,IPR,MJ1,E,SD )        
      CALL  MSDCOM( X,A,NS,LAG,ID,KSW,MJ1,E,SD )        
cxx   50 AIC(JJ+1) = NS*DLOG( SD ) + ID*(ID+1)
      AIC(JJ+1) = NS*DLOG( SD ) + ID*(ID+1)                             
   50 CONTINUE
C          ----------------------------------------                     
C          BAYESIAN WEIGHTS OF THE PRECEDING MODELS                     
C          ----------------------------------------                     
c-----------------------------  06/11/01
ccx      AICM = DMIN( AIC,KC )                                             
      AICM = AIC(1)
      DO 55  I=1,KC
cxx   55 IF( AIC(I) .LT. AICM )  AICM = AIC(I)
      IF( AIC(I) .LT. AICM )  AICM = AIC(I)
   55 CONTINUE
c-----------------------------
      CALL  BAYSWT( AIC,AICM,KC,2,C )                                   
cc      WRITE( 6,3 )     C(1) , AIC(1)                                    
cc      DO 60  I=2,KC1                                                    
cc      IM1 = I-1                                                         
cc   60 WRITE( 6,4 )     IM1 , C(I) , AIC(I)                              
C                                                                       
C          ------------------------                                     
C          AVERAGING OF THE MODELS                                      
C          -----------------------                                      
      EK = EK*C(1)**2                   
cxx      DO 70  II=1,LAG                                                   
cxx      DO 70  I=1,ID                                                     
      DO 72  II=1,LAG                                                   
      DO 71  I=1,ID
      DO 70  J=1,ID                                                     
cxx   70 B(I,J,II) = A(I,J,II)*C(1)                                        
      B(I,J,II) = A(I,J,II)*C(1)
   70 CONTINUE
   71 CONTINUE
   72 CONTINUE
cxx      DO 80  I=1,KD                                                     
      DO 81  I=1,KD
      DO 80  J=1,ID                                                     
      F1(I,J,1) = F1(I,J,1)*C(1)                                        
cxx   80 F2(I,J,1) = F2(I,J,1)*C(1)                                        
      F2(I,J,1) = F2(I,J,1)*C(1)                                        
   80 CONTINUE
   81 CONTINUE
C                                                                       
cxx      DO 90  JJ=1,KC                                                    
cxx        DO 90  I=1,KD                                                   
      DO 92  JJ=1,KC                                                    
        DO 91  I=1,KD                                                   
        DO 90  J=1,ID                                                   
        F1(I,J,1) = F1(I,J,1) + F1(I,J,JJ+1)*C(JJ+1)                    
cxx   90   F2(I,J,1) = F2(I,J,1) + F2(I,J,JJ+1)*C(JJ+1)
        F2(I,J,1) = F2(I,J,1) + F2(I,J,JJ+1)*C(JJ+1)                    
   90   CONTINUE
   91   CONTINUE
   92 CONTINUE
C          -----------------------------------------                    
C          PREDICTION ERROR VARIANCE MATRIX COMPUTED                    
C          -----------------------------------------                    
  100 IM = 0                                                            
      KC = KC + 1                                                       
      IF( KC .GT. KMAX )     KC = KMAX                                  
cxx      DO 110  II=1,LAG                                                  
cxx      DO 110  I=1,ID                                                    
      DO 112  II=1,LAG                                                  
      DO 111  I=1,ID
        IM = IM+1                                                       
        DO 110  J=1,ID                                                  
        G(I,J,II) = F1(IM,J,1)                                          
cxx  110   H(I,J,II) = F2(IM,J,1)                                          
        H(I,J,II) = F2(IM,J,1)
  110   CONTINUE
  111 CONTINUE
  112 CONTINUE
C                                                                       
      CALL  MARCOF( G,H,ID,LAG,MJ3,A,B )                                
cc      CALL  MSDCOM( X,A,Y,D,NS,LAG,ID,KSW,IPR,MJ1,MJ3,MJ5,E,SD )        
cxx      CALL  MSDCOM( X,A,NS,LAG,ID,KSW,IPR,MJ1,E,SD )
      CALL  MSDCOM( X,A,NS,LAG,ID,KSW,MJ1,E,SD )
      AICB = NS*DLOG( SD ) + 2.D0*(EK + KSW*ID) + ID*(ID+1)             
C                                                                       
      RETURN                                                            
cxx    3 FORMAT( ///1H ,13X,'AR-MODEL FITTED TO  !  BAYESIAN WEIGHTS  ! AIC
cxx     1 WITH RESPECT TO THE PRESENT DATA',/,10X,83(1H-),/,1H ,11X,'CURREN
cxx     2T BLOCK',9X,'!',F13.5,7X,'!',F21.3 )
cxx    4 FORMAT( 1H ,6X,I5,' PERIOD FORMER BLOCK  !',F13.5,7X,'!',F21.3 )
      E N D
