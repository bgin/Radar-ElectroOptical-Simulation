      SUBROUTINE UNIBARF( ZS,N,LAG,ZMEAN,SUM,SD,AIC,DIC,IMIN,AICM,SDMIN,
     *                    B1,C,D,B2,AICB,SDB,PN,A,SXX )
C
      INCLUDE 'timsac_f.h'
C
cc      PROGRAM  UNIBAR                                                   
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
C     TIMSAC 78.1.2.                                                    
C     ___        _                  __                                  
C     UNIVARIATE BAYESIAN METHOD OF AR MODEL FITTING                    
C                                                                       
C     THIS PROGRAM FITS AN AUTOREGRESSIVE MODEL BY A BAYESIAN PROCEDURE.
C     THE LEAST SQUARES ESTIMATES OF THE PARAMETERS ARE OBTAINED BY THE 
C     HOUSEHOLDER TRANSFORMATION.                                       
C     ----------------------------------------------------------------- 
C     THE BASIC STATISTIC AIC IS DEFINED BY                             
C                                                                       
C               AIC = N * LOG( SD )  +  2 * M  ,                        
C     WHERE                                                             
C       N:     DATA LENGTH,                                             
C       SD:    ESTIMATE OF INNOVATION VARIANCE, THE AVERAGE OF THE      
C              SQUARED RESIDUALS                                        
C       M:     ORDER OF THE MODEL                                       
C     ----------------------------------------------------------------- 
C     BAYESIAN WEIGHT OF THE M-TH ORDER MODEL IS DEFINED BY             
C            W(M)  =  CONST  *  C(M)  /  ( M + 1 ),                     
C     WHERE                                                             
C            CONST =  NORMALIZING CONSTANT                              
C            C(M)  =  EXP( -0.5*AIC(M) ).                               
C     ----------------------------------------------------------------- 
C     THE EQUIVALENT NUMBER OF FREE PARAMETERS FOR THE BAYESIAN MODEL IS
C     DEFINED BY                                                        
C            EK    =  D(1)**2 + ... + D(K)**2 + 1                       
C     WHERE D(J) IS DEFINED BY                                          
C            D(J) = W(J) + ... + W(K).                                  
C                                                                       
C     M IN THE DEFINITION OF AIC IS REPLACED BY EK TO DEFINE AN EQUIVALE
C     FOR A BAYESIAN MODEL.                                             
C     ----------------------------------------------------------------- 
C       REFERENCES:                                                     
C          H.AKAIKE(1978), "A BAYESIAN EXTENSION OF THE MINIMUM AIC     
C          PROCEDURE OF AUTOREGRESSIVE MODEL FITTING.",  RESEARCH MEMO. 
C          NO. 126, THE INSTITUTE OF STATISTICAL MATHEMATICS; TOKYO.   
C          G.KITAGAWA AND H.AKAIKE(1978), "A PROCEDURE FOR THE MODELING 
C          OF NON-STATIONARY TIME SERIES.",  ANN. INST. STATIST. MATH., 
C          30,B,351-363.                                                
C       --------------------------------------------------------------- 
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
C             REDATA                                                    
C             REDUCT                                                    
C             ARBAYS                                                    
C             NRASPE                                                    
C       --------------------------------------------------------------- 
C       INPUTS REQUIRED:                                                
C             LAG         :  ORDER OF THE AR-MODEL, MUST BE LESS THAN 10
C             MT          :  INPUT DEVICE FOR ORIGINAL DATA (MT=5 : CARD
C                                                                       
C               -- THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE REDA
C                                                                       
C             TITLE:    SPECIFICATION OF DATA                           
C             N:        DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 10000
C             DFORM:    INPUT DATA FORMAT SPECIFICATION STATEMENT.      
C                       -- FOR EXAMPLE --     (8F10.5)                  
C             (Z(I),I=1,N):  ORIGINAL DATA                              
C               --------------------------------------------------------
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: UNIBARF
C                                                                       
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL * 4   Z(10000) , TITLE(20)                                   
cc      REAL * 4   TITLE(20)
cc      DIMENSION  Z(10000)
cc      DIMENSION  X(200,101) , D(200) , A(100) , B(100)                  
cxx      DIMENSION  ZS(N), Z(N)
cxx      DIMENSION  X(N-LAG,LAG+1), D(LAG), A(LAG), B1(LAG), B2(LAG)
cxx      DIMENSION  SD(LAG+1), AIC(LAG+1), DIC(LAG+1), C(LAG+1)
cxx      DIMENSION  SXX(121)
      INTEGER :: N, LAG, IMIN
      REAL(8) :: ZS(N), ZMEAN, SUM, SD(LAG+1), AIC(LAG+1), DIC(LAG+1),
     1           AICM, SDMIN, B1(LAG), C(LAG+1), D(LAG), B2(LAG),
     2           AICB, SDB, PN, A(LAG), SXX(121)
      REAL(8) :: Z(N), X(N-LAG,LAG+1), B
C                                                                       
C        EXTERNAL SUBROUTINE DECLARATION:                               
C                                                                       
      EXTERNAL  SETX1                                                   
C
cc      CHARACTER(100)  IFLNAM,OFLNAM
cc      CALL FLNAM2( IFLNAM,OFLNAM,NFL )
cc      IF ( NFL.EQ.0 ) GO TO 999
cc      IF ( NFL.EQ.2 ) THEN
cc         OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR )
cc      ELSE
cc         CALL SETWND
cc      END IF
C
C        PARAMETERS:                                                    
C             MJ1:  ABSOLUTE DIMENSION FOR SUBROUTINE CALL              
C             ISW:  =0   OUTPUTS ARE SUPPRESSED                         
C                   >0   OUTPUTS ARE PRINTED OUT                        
cc      MJ1 = 200                                                         
      MJ1 = N-LAG                                                       
      ISW = 1                                                           
C                                                                       
cc      WRITE( 6,3 )                                                      
cc      WRITE( 6,4 )                                                      
C
cc      MT = 5
cc      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )                                                                     
CC      READ( 5,1 )     MT                                                
cc      READ( 5,1 )     LAG                                               
cc      WRITE( 6,2 )    LAG , MT                                          
C                                                                       
C                                                                       
C          +---------------------------------------+                  +-
C          ! ORIGINAL DATA INPUT AND MEAN DELETION !                  ! 
C          +---------------------------------------+                  +-
C
cc      CALL  REDATA( Z,N,MT,TITLE )                                      
cc      CLOSE( MT )
      CALL REDATA( ZS,Z,N,ZMEAN,SUM )
      K = LAG                                                           
      NMK = N-LAG                                                       
C                                                                       
C          +-----------------------+                                  +-
C          ! HOUSEHOLDER REDUCTION !                                  ! 
C          +-----------------------+                                  +-
C                                                                       
      N1 = N+1
cc      CALL  REDUCT( SETX1,Z,D,NMK,0,K,MJ1,LAG,X )                       
      CALL  REDUCT( SETX1,Z,NMK,0,K,MJ1,LAG,X )                       
C                                                                       
C         +-----------------------------------------------+           +-
C         ! AUTOREGRESSIVE MODEL FITTING (BAYESIAN MODEL) !           ! 
C         +-----------------------------------------------+           +-
C                                                                       
cc      CALL  ARBAYS( X,D,K,LAG,NMK,ISW,TITLE,MJ1,A,B,SDB,AICB )          
cxx      CALL  ARBAYS( X,D,K,LAG,NMK,ISW,MJ1,SD,AIC,DIC,AICM,SDMIN,IMIN,
      CALL  ARBAYS( X,D,K,NMK,ISW,MJ1,SD,AIC,DIC,AICM,SDMIN,IMIN,
     *              A,B1,B2,C,SDB,PN,AICB )
C                                                                       
C          +------------------------+                                 +-
C          ! POWER SPECTRUM DISPLAY !                                 ! 
C          +------------------------+                                 +-
C                                                                       
cc      CALL  NRASPE( SDB,A,B,K,0,120,TITLE )                             
      CALL  NRASPE( SDB,A,B,K,0,120,SXX )
cc      GO TO 999
C
cc  900 CONTINUE
cc      WRITE(6,600) IVAR,OFLNAM
cc  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,//,5X,100A)
cc       GO TO 999
C
cc  910 CONTINUE
cc      IF ( NFL.EQ.2 ) CLOSE( 6 )
cc#ifdef __linux__
C      reopen #6 as stdout
cc      IF ( NFL.EQ.2 ) OPEN(6, FILE='/dev/fd/1')
cc#endif
C /* __linux__ */
cc      WRITE(6,610) IVAR,IFLNAM
cc  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//,5X,100A)
C                                                                       
cc  999 CONTINUE
      RETURN
cxx    1 FORMAT( 16I5 )                                                    
cxx    2 FORMAT( 1H ,I4,'-TH ORDER BAYESIAN MODEL IS FITTED',/,1H ,2X,     
cxx     1'ORIGINAL DATA INPUT DEVICE  MT =',I4 )                           
cxx    3 FORMAT( ' PROGRAM TIMSAC 78.1.2',/'   EXPONENTIALLY WEIGHTED BAYES
cxx     1IAN AUTOREGRESSIVE MODEL FITTING;  SCALAR CASE')                  
cxx    4 FORMAT(  '   < AUTOREGRESSIVE MODEL >',/,1H ,10X,'Z(I) = A(1)*Z(I-
cxx     11) + A(2)*Z(I-2) +  ...  + A(M)*Z(I-M) + E(I)',/,'   WHERE',/,11X,
cxx     2'M:     ORDER OF THE MODEL',/,11X,'E(I):  GAUSSIAN WHITE NOISE WIT
cxx     3H MEAN 0  AND  VARIANCE SD(M).' )                                 
      E N D                                                             
