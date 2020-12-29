      SUBROUTINE  EXSARF( Z1,N,LAG,ZMEAN,SUM,SD,AIC,DIC,M1,AMIN,SDM1,A1,
cx     *                    SDM2,A2,TMP,IER,JER )
     *                    SDM2,A2,JER )
c----------------   for  isw = 2   --------------------
cxc     *                    SDM22,A22,TMP )
c     *                    SDM22,A22,JER )
C
      INCLUDE 'timsac_f.h'
C
cc      PROGRAM  EXSAR                                                    
C.......................................................................
C.....PLANNED BY H.AKAIKE...............................................
C.....DESIGNED BY H.AKAIKE AND G.KITAGAWA...............................
C.....PROGRAMMED BY G.KITAGAWA AND F.TADA...............................
C.....ADDRESS: THE INSTITUTE OF STATISTICAL MATHEMATICS, 4-6-7 MINAMI-AZ
C..............MINATO-KU, TOKYO 106, JAPAN..............................
C.....DATE OF THE LATEST REVISION:  JUN. 29, 1979.......................
C.......................................................................
C.....THIS PROGRAM WAS ORIGINALLY PUBLISHED IN "TIMSAC-78", BY H.AKAIKE,
C.....G.KITAGAWA, E.ARAHATA AND F.TADA, COMPUTER SCIENCE MONOGRAPHS, NO.
C.....THE INSTITUTE OF STATISTICAL MATHEMATICS, TOKYO, 1979.............
C.......................................................................
C     TIMSAC 78.5.1                                                     
C     __                                 _      __                      
C     EXACT MAXIMUM LIKELIHOOD METHOD OF SCALAR AR-MODEL FITTING        
C                                                                       
C     THIS PROGRAM PRODUCES EXACT MAXIMUM LIKELIHOOD ESTIMATES OF THE   
C     PARAMETERS OF A SCALAR AR-MODEL.                                  
C                                                                       
C     THE AR-MODEL IS GIVEN BY                                          
C                                                                       
C               Z(I) = A(1)*Z(I-1) + ... + A(K)*Z(I-K) + E(I)           
C                                                                       
C     WHERE E(I) IS A ZERO MEAN WHITE NOISE.                            
C                                                                       
C     --------------------------------------------------------------    
C     THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:    
C             REDATA                                                    
C             REDUCT                                                    
C             ARMFIT                                                    
C             RECOEF                                                    
C             ARMLE                                                     
C             PRINTA                                                    
C     --------------------------------------------------------------    
C     INPUTS REQUIRED:                                                  
C          MT:      INPUT DEVICE SPECIFICATION (MT=5: CARD READER)      
C          LAG:     UPPER LIMIT OF AR-ORDER, MUST BE LESS THAN 51       
C                                                                       
C     --  THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE REDATA  --   
C          TITLE:  TITLE OF DATA                                        
C          N:      DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 10000     
C          DFORM:  INPUT DATA FORMAT SPECIFICATION STATEMENT            
C                  -- EXAMPLE --     (8F10.5)                           
C          (Z(I),I=1,N):  ORIGINAL DATA                                 
C     ---------------------------------------------------------------   
C                                                                       
cc      !DEC$ ATTRIBUTES DLLEXPORT :: EXSARF
C
cxx      IMPLICIT  REAL *8  ( A-H , O-Z )                                  
CC      REAL * 4  Z(10000) , TITLE(20)                                    
cc      REAL * 4   TITLE(20)
cc      DIMENSION  Z(10000)
cc      DIMENSION  X(200,51) , D(200) , A(50)                             
cxx      DIMENSION  Z1(N), Z(N)
cxx      DIMENSION  X(N-LAG,LAG+1), A1(LAG), A2(LAG)
c----------------    isw = 2   -----------------
cxx      DIMENSION  A22(LAG,LAG), SDM22(LAG)
c-----------------------------------------------

      INTEGER :: N, LAG, M1, JER
      REAL(8) :: Z1(N), ZMEAN, SUM, SD(LAG+1), AIC(LAG+1), DIC(LAG+1),
     1           AMIN, SDM1, A1(LAG), SDM2, A2(LAG)
      REAL(8) :: Z(N), X(N-LAG,LAG+1)
c----------------    isw = 2   -----------------
      REAL(8) :: A22(LAG,LAG), SDM22(LAG)
c-----------------------------------------------

CC      COMMON     / AAA /  N , Z                                         
cc      COMMON     / AAA /  N
cc      COMMON     / BBB /  Z
cc      COMMON     / CCC /  ISW , IPR                                     
C
cxx      DIMENSION  SD(LAG+1), AIC(LAG+1), DIC(LAG+1)
cx      INTEGER*1  TMP(1)
cx      CHARACTER  CNAME*80
C
cc      CHARACTER(100) IFLNAM,OFLNAM
C                                                                       
C       EXTERNAL SUBROUTINE DECLARATION:                                
C                                                                       
      EXTERNAL   SETX1                                                  
C
cc      CALL FLNAM2( IFLNAM,OFLNAM,NFL )
cc      IF ( NFL.EQ.0 ) GO TO 999
cc      IF ( NFL.EQ.2 ) THEN
cc         OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR )
cc      ELSE
cc         CALL SETWND
cc      END IF
cx      IER=0
cx      LU=3
cx      DO 100 I = 1,80
cx         CNAME(I:I) = ' '
cx  100 CONTINUE
cx      I = 1
cx      IFG = 1
cx      DO WHILE( (IFG.EQ.1) .AND. (I.LE.80) )
cx	   IF ( TMP(I).NE.ICHAR(' ') ) THEN
cx            CNAME(I:I) = CHAR(TMP(I))
cx            I = I+1
cx         ELSE
cx            IFG = 0
cx         END IF
cx      END DO
cx      IF ( I.GT.1 ) THEN
cx         IFG = 1
cx         OPEN (LU,FILE=CNAME,IOSTAT=IVAR)
cx         IF (IVAR .NE. 0) THEN
cxcx            WRITE(*,*) ' ***  exsar temp FILE OPEN ERROR :',CNAME,IVAR
cx            IER=IVAR
cx            IFG=0
cx         END IF
cx      END IF
C                                                                       
C          PARAMETERS:                                                  
C             MJ1:  ABSOLUTE DIMENSION FOR SUBROUTINE CALL              
C             ISW:  =1  PARAMETERS OF MAICE MODEL ONLY ARE REQUESTED    
C                   =2  PARAMETERS OF ALL MODELS ARE REQUESTED          
C             IPR:  PRINT OUT CONTROL IN THE NON-LINEAR OPTIMIZATION PRO
C                                                                       
cc      MJ1 = 200                                                         
      MJ1 = N-LAG
      ISW = 1                                                           
C                                                                       
CC      READ( 5,1 )     MT                                                
cc      MT = 5
cc      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD')
cc      READ( 5,1 )     LAG                                               
cc      WRITE( 6,2 )                                                      
cc      WRITE( 6,3 )   LAG , MT                                           
C                                                                       
C          +-----------------------------------------+                  
C          ! ORIGINAL DATA LOADING AND MEAN DELETION !                  
C          +-----------------------------------------+                  
C                                                                       
cc      CALL  REDATA( Z,N,MT,TITLE )                                      
cc      CLOSE( MT )
      CALL  REDATA( Z1,Z,N,ZMEAN,SUM )
      K = LAG                                                           
      NMK = N - K                                                       
C                                                                       
C          +-----------------------+                                    
C          ! HOUSEHOLDER REDUCTION !                                    
C          +-----------------------+                                    
C                                                                       
cc      CALL  REDUCT( SETX1,Z,D,NMK,0,K,MJ1,LAG,X )                       
      CALL  REDUCT( SETX1,Z,NMK,0,K,MJ1,LAG,X )                       
C                                                                       
C          +------------------+                                         
C          ! AR MODEL FITTING !                                         
C          +------------------+                                         
C                                                                       
cc      CALL  ARMFIT( X,K,LAG,NMK,ISW,TITLE,MJ1,A,SD,M )                  
cx      CALL  ARMFIT( X,K,LAG,NMK,ISW,MJ1,A1,M1,SD,AIC,DIC,SDM1,AMIN,
cx     *              IFG,LU )
      CALL  ARMFIT( X,K,LAG,NMK,ISW,MJ1,A1,M1,SD,AIC,DIC,SDM1,AMIN )
      DO 5  I = 1,K
cxx    5 A2(I) = A1(I)
      A2(I) = A1(I)
    5 CONTINUE
C                                                                       
C          +-------------------------------------------+                
C          ! MAXIMIZATION OF EXACT LIKELIHOOD FUNCTION !                
C          +-------------------------------------------+                
C                                                                       
      JER = 0
      IF( ISW .EQ. 2 )  GO TO 10                                        
C                                                                       
      IPR = 7                                                           
cc      CALL  ARMLE( Z,N,M,K,TITLE,A )                                    
ccx      CALL  ARMLE( Z,N,M1,K,A2,SDM2,ISW,IPR,IFG,LU )
cx      CALL  ARMLE( Z,N,M1,K,A2,SDM2,ISW,IPR,IFG,LU,JER )
cxx      CALL  ARMLE( Z,N,M1,K,A2,SDM2,ISW,IPR,JER )
      CALL  ARMLE( Z,N,M1,K,A2,SDM2,ISW,JER )
cx      IF( JER .NE. 0 ) RETURN
cx      IF (IFG.NE.0) CLOSE(LU)
      RETURN
C                                                                       
   10 DO 20  M=1,K                                                      
cc      CALL  RECOEF( X,M,K,MJ1,A )                                       
      CALL  RECOEF( X,M,K,MJ1,A2 )
      IPR = 5                                                           
cc   20 CALL  ARMLE( Z,N,M,K,TITLE,A )                                    
ccx      CALL  ARMLE( Z,N,M,K,A2,SDM2,ISW,IPR,IFG,LU )
cx      CALL  ARMLE( Z,N,M,K,A2,SDM2,ISW,IPR,IFG,LU,JER )
cxx      CALL  ARMLE( Z,N,M,K,A2,SDM2,ISW,IPR,JER )
      CALL  ARMLE( Z,N,M,K,A2,SDM2,ISW,JER )
      IF( JER .NE. 0 ) RETURN
      DO 15  I = 1,M
cxx   15 A22(I,M) = A2(I)
      A22(I,M) = A2(I)
   15 CONTINUE
      SDM22(M) = SDM2
   20 CONTINUE
cc      GO TO 999
C
cc  900 CONTINUE
cc      WRITE(6,600) IVAR,OFLNAM
cc  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,//,5X,100A)
cc      GO TO 999
C
cc  910 CONTINUE
cc      IF ( NFL.EQ.2 ) CLOSE( 6 )
cc#ifdef __linux__
ccC     reopen #6 as stdout
cc      IF ( NFL.EQ.2 ) OPEN(6, FILE='/dev/fd/1')
cc#endif
ccC /* __linux__ */
cc      WRITE(6,610) IVAR,IFLNAM
cc  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//,5X,100A)
C                                                                       
cxx    1 FORMAT( 16I5 )                                                    
cxx    2 FORMAT( 1H ,'PROGRAM 78.5.1',/,'   EXACT MAXIMUM LIKELIHOOD METHOD
cxx     * OF AUTOREGRESSIVE MODEL FITTING;  SCALAR CASE',/,'   < BASIC AUTO
cxx     *REGRESSIVE MODEL >',/,1H ,10X, 'Z(I) = A(1)*Z(I-1) + ... + A(K)*Z(
cxx     *I-K) + E(I)',/,1H ,'WHERE',/,11X,'K:     ORDER OF THE MODEL',/,11X
cxx     *,'E(I):  ZERO MEAN WHITE NOISE' )                                 
cxx    3 FORMAT( ///1H ,'  FITTING UP TO THE ORDER  K =',I3,'  IS TRIED',/,
cxx     *1H ,'  ORIGINAL DATA INPUT DEVICE  MT =',I3 )                     
C                                                                       
cc  999 CONTINUE
cx      IF (IFG.NE.0) CLOSE(LU)                                                                       
      RETURN
      E N D                                                             
cc      SUBROUTINE  ARMLE( Z,N,K,L,TITLE,A )                              
ccx      SUBROUTINE  ARMLE( Z,N,K,L,A,SDM,ISW,IPR,IFG,LU )
cx      SUBROUTINE  ARMLE( Z,N,K,L,A,SDM,ISW,IPR,IFG,LU,JER )
cxx      SUBROUTINE  ARMLE( Z,N,K,L,A,SDM,ISW,IPR,JER )
      SUBROUTINE  ARMLE( Z,N,K,L,A,SDM,ISW,JER )
C.....DATE OF THE LATEST REVISION:  JUN. 29, 1979.......................
C                                                                       
C     THIS SUBROUTINE PRODUCES EXACT MAXIMUM LIKELIHOOD ESTIMATES OF THE
C     PARAMETERS OF AN AUTOREGRESSIVE MODEL                             
C                                                                       
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             DAVIDN                                                    
C             PRINTA                                                    
C       ----------------------------------------------------------------
C          INPUTS:                                                      
C             Z:       ORIGINAL DATA                                    
C             N:       DATA LENGTH                                      
C             K:       ORDER OF THE AR-MODEL                            
C             L:       HIGHEST POSSIBLE ORDER OF THE MODEL              
C             TITLE:   TITLE OF DATA                                    
C                                                                       
C          OUTPUTS:                                                     
C             A:       MAXIMUM LIKELIHOOD ESTIMATES OF AR-COEFFICIENTS  
C                                                                       
C                                                                       
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL * 4   Z , TITLE                                              
cc      REAL * 4   TITLE
cc      DIMENSION  Z(1) , A(1) , TITLE(1)                                 
cc      DIMENSION  R(31,31) , C(31)                                       
cxx      DIMENSION  Z(N), A(K)
cxx      DIMENSION  R(K+1,K+1), C(L+1)                                       
      INTEGER :: N, K, L, ISW, JER 
      REAL(8) :: Z(N), A(K), SDM
      REAL(8) :: R(K+1,K+1), C(L+1), SUM, F0, AIC, SD
CC      DIMENSION  TTL(4)                                                 
cc      REAL * 4 TTL(8)
cc      COMMON     / DDD /  R , F , AIC , SD                              
CC      DATA         TTL / 8H  MAXIMU,8HM LIKELI,8HHOOD EST,8HIMATE  . /  
cc      DATA         TTL / 4H  MA,4HXIMU,4HM LI,4HKELI,
cc     1					4HHOOD,4H EST,4HIMAT,4HE  . / 
C                                                                       
C       EXTERNAL FUNCTION DECRALATION                                   
C                                                                       
      EXTERNAL  FUNCT                                                   
      EXTERNAL  HESIAN                                                  
C                                                                       
      IHES = 1                                                          
C                                                                       
cc      WRITE( 6,3 )                                                      
cx      IF( IFG.NE.0 ) WRITE( LU,3 )  K
      NMLP1 = N - L + 1                                                 
      K1 = K + 1                                                        
      L1 = L + 1                                                        
      N1 = N + 1                                                        
      NML = N - L                                                       
      DO  20  II=1,L1                                                   
      J = II - 1                                                        
      SUM = 0.D0                                                        
      DO 10  I=L1,NML                                                   
      IMJ = I - J                                                       
cxx   10 SUM = SUM + Z(I)*Z(IMJ)                                           
      SUM = SUM + Z(I)*Z(IMJ)  
   10 CONTINUE
cxx   20 C(II) = SUM                                                       
      C(II) = SUM
   20 CONTINUE
C                                                                       
C       COVARIANCE MATRIX R COMPUTATION                                 
C                                                                       
cxx      DO 70  II=1,K1                                                    
      DO 71  II=1,K1 
      KMI = K1 - II + 1                                                 
      NMI = N1 - II                                                     
      DO 70  JJ=II,K1                                                   
      JMI = JJ - II                                                     
      SUM = C(JMI+1)                                                    
      IF( KMI .GT. L )  GO TO 40                                        
      DO 30  KK=KMI,L                                                   
      KI = KK - JMI                                                     
cxx   30 SUM = SUM + Z(KK)*Z(KI)                                           
      SUM = SUM + Z(KK)*Z(KI) 
   30 CONTINUE
   40 CONTINUE                                                          
      IF( NMLP1 .GT. NMI )  GO TO 60                                    
      DO 50  KK=NMLP1,NMI                                               
      KI = KK - JMI                                                     
cxx   50 SUM = SUM + Z(KK)*Z(KI)                                           
      SUM = SUM + Z(KK)*Z(KI)
   50 CONTINUE
   60 R(II,JJ) = SUM                                                    
      R(JJ,II) = SUM                                                    
   70 CONTINUE
   71 CONTINUE                                                          
C                                                                       
C       DAVIDON'S MINIMIZATION PROCEDURE                                
C                                                                       
      F0 = 1.0D60                                                       
      DO  80   I=1,5                                                    
C                                                                       
cc      CALL  DAVIDN( FUNCT,HESIAN,A,K,IHES )                             
ccx      CALL  DAVIDN( FUNCT,HESIAN,Z,N,A,K,R,IHES,ISW,IPR,AIC,SD,IFG,LU )
cx      CALL  DAVIDN( FUNCT,HESIAN,Z,N,A,K,R,IHES,ISW,IPR,AIC,SD,
cx     *                     IFG,LU,JER )
cxx      CALL  DAVIDN( FUNCT,HESIAN,Z,N,A,K,R,IHES,ISW,IPR,AIC,SD,JER )
      CALL  DAVIDN( FUNCT,HESIAN,Z,N,A,K,R,IHES,ISW,AIC,SD,JER )
      IF( JER .NE. 0 ) RETURN
C                                                                       
      IF( F0-AIC .LT. 0.001 )     GO TO 90                              
cxx   80 F0 = AIC
      F0 = AIC
   80 CONTINUE                                                          
   90 CONTINUE                                                          
C                                                                       
cc      CALL  PRINTA( A,SD,K,TTL,8,TITLE,1,N )                            
      SDM = SD
C                                                                       
      RETURN                                                            
cc    3 FORMAT( //1H ,29(1H-),/,' EXACT LIKELIHOOD MAXIMIZATION',/,1H ,29(
cc     11H-) )                                                            
cxx    3 FORMAT( //1H ,29(1H-),/,' EXACT LIKELIHOOD MAXIMIZATION',10x,
cxx     *'ORDER OF THE AR-MODEL =',i5,/1H ,29(1H-) )                                                            
cxx   65 FORMAT( 1H ,10D13.5 )                                             
      END                                                               
cc      SUBROUTINE  DAVIDN( FUNCT,HESIAN,X,N,IHES )                       
cxx      SUBROUTINE  DAVIDN( FUNCT,HESIAN,Z,NZ,X,N,R,IHES,ISW,IPR,AIC,SD,
ccx     *                    IFG,LU )
cx     *                    IFG,LU,JER )
cxx     *                    JER )
      SUBROUTINE  DAVIDN( FUNCT,HESIAN,Z,NZ,X,N,R,IHES,ISW,AIC,SD,JER )
C.....DATE OF THE LATEST REVISION:  JUN. 29, 1979.......................
C                                                                       
C          MINIMIZATION BY DAVIDON-FLETCHER-POWELL PROCEDURE            
C                                                                       
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             FUNCT                                                     
C             HESIAN                                                    
C             LINEAR                                                    
C       ----------------------------------------------------------------
C          INPUTS:                                                      
C             FUNCT:   EXTERNAL FUNCTION SPECIFICATION                  
C             HESIAN:  EXTERNAL FUNCTION SPECIFICATION                  
C             X:       VECTOR OF INITIAL VALUES                         
C             N:       DIMENSION OF THE VECTOR X                        
C             IHES:    =0   INVERSE OF HESSIAN MATRIX IS NOT AVAILABLE  
C                      =1   INVERSE OF HESSIAN MATRIX IS AVAILABLE      
C                                                                       
C          OUTPUT:                                                      
C             X:       VECTOR OF MINIMIZING SOLUTION                    
C                                                                       
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cc      DIMENSION  X(50) , DX(50) , G(50) , G0(50) , Y(50)                
cc      DIMENSION  H(50,50) , WRK(50) , S(50)                             
cc      DIMENSION  R(31,31)                                               
cxx      DIMENSION  Z(NZ)
cxx      DIMENSION  X(N) , DX(N) , G(N) , G0(N) , Y(N)                
cxx      DIMENSION  H(N,N) , WRK(N) , S(N)                             
cxx      DIMENSION  R(N+1,N+1)
      INTEGER :: NZ, N, IHES, ISW, JER
      REAL(8) :: Z(NZ), X(N), R(N+1,N+1), AIC, SD
      REAL(8) :: DX(N), G(N), G0(N), Y(N), H(N,N), WRK(N), S(N),
     1           TAU1, TAU2, EPS1, EPS2, RAMDA, CONST1, SUM,
     2           S1, S2, STEM, SS, DS2, GTEM, ED, F, XM, XMB                                                           
C
cc      COMMON     / CCC /  ISW, IPR                                      
cc      COMMON     / DDD /  R , F , AIC , SD                              
      DATA  TAU1 , TAU2  /  1.0D-5 , 1.0D-5  /                          
      DATA  EPS1 , EPS2  / 1.0D-5 , 1.0D-5  /                           
      EXTERNAL  FUNCT
C
      RAMDA = 0.5D0                                                     
      CONST1 = 1.0D-70                                                  
C                                                                       
C          INITIAL ESTIMATE OF INVERSE OF HESSIAN                       
C
      H(1:N,1:N) = 0.0D00
      S(1:N) = 0.0D00                                                     
      DX(1:N) = 0.0D00                                                    
      DO  20   I=1,N                                                    
cxx      DO  10   J=1,N                                                    
cxx   10 H(I,J) = 0.0D00
cxx      S(I) = 0.0D00                                                     
cxx      DX(I) = 0.0D00                                                    
cxx   20 H(I,I) = 1.0D00                                                   
      H(I,I) = 1.0D00   
   20 CONTINUE
      ISW = 0                                                           
C                                                                       
cc      CALL  FUNCT( N,X,XM,G,IG )                                        
cx      CALL  FUNCT( Z,NZ,N,X,R,ISW,XM,G,AIC,SD,F,IG )
      CALL  FUNCT( Z,NZ,N,X,R,ISW,XM,G,AIC,SD,F,IG,JER )
      IF( JER .NE. 0 ) RETURN
C                                                                       
cc      IF( IPR .GE. 2 )   WRITE( 6,340 )   XM, SD, AIC                   
cx      IF( (IPR.GE.2) .AND. (IFG.NE.0) )  WRITE( LU,340 )  XM, SD, AIC
C                                                                       
C          INVERSE OF HESSIAN COMPUTATION (IF AVAILABLE)                
C                                                                       
cc      IF( IHES .EQ. 1 )   CALL  HESIAN( X,N,H )                         
      IF( IHES .EQ. 1 )   CALL  HESIAN( X,N,NZ,R,SD,H )
C                                                                       
      ICC = 0                                                           
C      ITERATION                                                        
11110 CONTINUE                                                          
      ICC = ICC + 1                                                     
      DO  11111   IC=1,N                                                
      IF( IC .EQ. 1 .AND. ICC .EQ. 1 )     GO TO 120                    
C                                                                       
      DO  40   I=1,N                                                    
cxx   40 Y(I) = G(I) - G0(I)                                               
      Y(I) = G(I) - G0(I)
   40 CONTINUE
      DO  60   I=1,N                                                    
      SUM = 0.0D00                                                      
      DO  50   J=1,N                                                    
cxx   50 SUM = SUM + Y(J) * H(I,J)                                         
      SUM = SUM + Y(J) * H(I,J)
   50 CONTINUE
cxx   60 WRK(I) = SUM                                                      
      WRK(I) = SUM
   60 CONTINUE
      S1 = 0.0D00                                                       
      S2 = 0.0D00                                                       
      DO  70   I=1,N                                                    
      S1 = S1 + WRK(I) * Y(I)                                           
cxx   70 S2 = S2 + DX(I) * Y(I)                                            
      S2 = S2 + DX(I) * Y(I)
   70 CONTINUE
      IF( S1.LE.CONST1 .OR. S2.LE.CONST1 )  GO TO 900                   
      IF( S1 .LE. S2 )     GO TO 100                                    
C                                                                       
C          UPDATE THE INVERSE OF HESSIAN MATRIX                         
C                                                                       
C               ---  DAVIDON-FLETCHER-POWELL TYPE CORRECTION  ---       
C                                                                       
cxx      DO  90   I=1,N                                                    
      DO  91   I=1,N 
      DO  90   J=I,N                                                    
      H(I,J) = H(I,J) + DX(I)*DX(J)/S2 - WRK(I)*WRK(J)/S1               
cxx   90 H(J,I) = H(I,J)                                                   
      H(J,I) = H(I,J)
   90 CONTINUE
   91 CONTINUE
      GO TO  120                                                        
C                                                                       
C               ---  FLETCHER TYPE CORRECTION  ---                      
C                                                                       
  100 CONTINUE                                                          
      STEM = S1 / S2 + 1.0D00                                           
cxx      DO  110   I=1,N                                                   
      DO  111   I=1,N
      DO  110   J=I,N                                                   
      H(I,J) = H(I,J)- (DX(I)*WRK(J)+WRK(I)*DX(J)-DX(I)*DX(J)*STEM)/S2  
cxx  110 H(J,I) = H(I,J)                                                   
      H(J,I) = H(I,J)
  110 CONTINUE
  111 CONTINUE
C                                                                       
C                                                                       
C                                                                       
  120 CONTINUE                                                          
      SS = 0.0D00                                                       
      DO  150   I=1,N                                                   
      SUM = 0.0D00                                                      
      DO  140   J=1,N                                                   
cxx  140 SUM = SUM + H(I,J)*G(J)                                           
      SUM = SUM + H(I,J)*G(J) 
  140 CONTINUE
      SS = SS + SUM * SUM                                               
cxx  150 S(I) = -SUM                                                       
      S(I) = -SUM
  150 CONTINUE
C                                                                       
C                                                                       
      S1 = 0.0D00                                                       
      S2 = 0.0D00                                                       
      DO  170   I=1,N                                                   
      S1 = S1 + S(I)*G(I)                                               
cxx  170 S2 = S2 + G(I)*G(I)
      S2 = S2 + G(I)*G(I)
  170 CONTINUE                                               
      DS2 = DSQRT(S2)                                                   
      GTEM = DABS(S1) / DS2                                             
      IF( GTEM .LE. TAU1  .AND.  DS2 .LE. TAU2 )     GO TO  900         
      IF( S1 .LT. 0.0D00 )     GO TO  200                               
      H(1:N,1:N) = 0.0D00
      DO  190   I=1,N                                                   
cxx      DO  180   J=1,N                                                   
cxx  180 H(I,J) = 0.0D00                                                   
      H(I,I) = 1.0D00                                                   
cxx  190 S(I) = -S(I)                                                      
      S(I) = -S(I)
  190 CONTINUE
  200 CONTINUE                                                          
C                                                                       
      ED = XM                                                           
C                                                                       
C          LINEAR  SEARCH                                               
C                                                                       
cc      CALL  LINEAR( FUNCT,X,S,RAMDA,ED,N,IG )                           
cxx      CALL  LINEAR( FUNCT,Z,NZ,X,S,RAMDA,ED,N,R,AIC,SD,F,ISW,IPR,IG,
ccx     *              IFG,LU )
cx     *              IFG,LU,JER )
cxx     *              JER )
      CALL  LINEAR( FUNCT,Z,NZ,X,S,RAMDA,ED,N,R,AIC,SD,F,ISW,IG,JER )
      IF( JER .NE. 0 ) RETURN
C                                                                       
cc      IF( IPR .GE. 2 )   WRITE( 6,330 )   RAMDA, F, SD, AIC             
cx      IF( (IPR.GE.2) .AND. (IFG.NE.0) )
cx     *                     WRITE( LU,330 )  RAMDA, F, SD, AIC
C                                                                       
      S1 = 0.0D00                                                       
      DO  210   I=1,N                                                   
      DX(I) = S(I) * RAMDA                                              
      S1 = S1 + DX(I) * DX(I)                                           
      G0(I) = G(I)                                                      
cxx  210 X(I) = X(I) + DX(I)                                               
      X(I) = X(I) + DX(I)
  210 CONTINUE
      XMB = XM                                                          
      ISW = 0                                                           
C                                                                       
cc      CALL  FUNCT( N,X,XM,G,IG )                                        
cx      CALL  FUNCT( Z,NZ,N,X,R,ISW,XM,G,AIC,SD,F,IG )
      CALL  FUNCT( Z,NZ,N,X,R,ISW,XM,G,AIC,SD,F,IG,JER )
      IF( JER .NE. 0 ) RETURN
C                                                                       
      S2 = 0.D0                                                         
      DO  220     I=1,N                                                 
cxx  220 S2 = S2 + G(I)*G(I)
      S2 = S2 + G(I)*G(I)
  220 CONTINUE                                               
      IF( DSQRT(S2) .GT. TAU2 )   GO TO  11111                          
      IF( XMB/XM-1.D0 .LT. EPS1  .AND.  DSQRT(S1) .LT. EPS2 )  GO TO 900
11111 CONTINUE                                                          
      IF( ICC .GE. 5 )     GO TO 900                                    
      GO TO 11110                                                       
  900 CONTINUE                                                          
cx      IF( IPR .LE. 0 )   RETURN                                         
cc      WRITE( 6,600 )                                                    
cc      WRITE( 6,610 )     (X(I),I=1,N)                                   
cc      WRITE( 6,601 )                                                    
cc      WRITE( 6,610 )     (G(I),I=1,N)                                   
cx      IF( IFG .EQ. 0 )   RETURN
cx      WRITE( LU,600 )                                                    
cx      WRITE( LU,610 )     (X(I),I=1,N)                                   
cx      WRITE( LU,601 )                                                    
cx      WRITE( LU,610 )     (G(I),I=1,N)                                   
      RETURN                                                            
cxx  330 FORMAT( 1H ,'LAMBDA =',D15.7,3X,'(-1)LOG LIKELIHOOD =',D23.15,3X, 
cxx     * 'SD =',D22.15,5X,'AIC =',D23.15 )                                
cxx  340 FORMAT( 1H ,26X,'(-1)LOG-LIKELIHOOD =',D23.15,3X,'SD =',D22.15,5X,
cxx     *  'AIC =',D23.15 )                                                
cxx  600 FORMAT( 1H ,'-----  X  -----' )                                   
cxx  601 FORMAT( 1H ,'***  GRADIENT  ***' )                                
cxx  610 FORMAT( 1H ,10D13.5 )                                             
      END
cc      SUBROUTINE  FUNCT( M,A,F,G,IFG )                                  
cx      SUBROUTINE  FUNCT( Z,N,M,A,R,ISW,F,G,AIC,SD,FF,IFG )
      SUBROUTINE  FUNCT( Z,N,M,A,R,ISW,F,G,AIC,SD,FF,IFG,JER )
C                                                                       
C     THIS SUBROUTINE COMPUTES THE EXACT LIKELIHOOD AND ITS GRADIENT OF 
C     THE M-TH ORDER AR-MODEL.                                          
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             ARCOEF                                                    
C             PARCOR                                                    
C             SUBDET                                                    
C       ----------------------------------------------------------------
C                                                                       
C          INPUTS:                                                      
C             M:    ORDER OF THE AR MODEL                               
C             A:    VECTOR OF AR-COEFFICIENTS                           
C                                                                       
C          OUTPUTS:                                                     
C             F:    LIKELIHOOD OF THE AR-MODEL                          
C             G:    GRADIENT OF THE LIKELIHOOD FUNCTION                 
C             IFG:  =0    ;IF MODEL IS STATIONALY                       
C                   =1    ;IF MODEL IS NON-STATIONALY                   
C                                                                       
C                                                                       
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL * 4   Z                                                      
cc      DIMENSION  Z(10000)
CC      COMMON     / AAA /  N , Z(1)                                      
cc      COMMON     / AAA /  N
cc      COMMON     / BBB /  Z
cc      COMMON     / CCC /  ISW, IDUMMY                                          
cc      COMMON     / DDD /  R , FF , AIC , SD                             
cc      DIMENSION   R(31,31) , T(31,31) , U(31,31) , S(31,31)              
cc      DIMENSION   A(20) , B(30) , G(30)                                  
cxx      DIMENSION   Z(N)
cxx      DIMENSION   R(M+1,M+1) , T(M+1,M+1) , U(M+1,M+1) , S(M+1,M+1)              
cxx      DIMENSION   A(M) , B(M) , G(M)                                  
      INTEGER :: N, M, ISW, IFG, JER
      REAL(8) :: Z(N), A(M), R(M+1,M+1), F, G(M), AIC, SD, FF
      REAL(8) :: T(M+1,M+1), U(M+1,M+1), S(M+1,M+1), B(M),
     1           DN, DN1, SUM, DSD, USUM, DETT, UDET
cc      MJ = 31                                                           
      MJ = M+1                                                           
C                                                                       
  280 IFG = 0                                                           
C                                                                       
      DN = N                                                            
      DN1 = 1.0D00 / DN                                                 
      MP1 = M + 1                                                       
      M1 = MP1 / 2                                                      
C                                                                       
C   **  COMPUTATION OF F  **                                            
C                                                                       
C                                                                       
C          INVERSE OF COVARIANCE MATRIX COMPUTATION                     
C                                                                       
C                                                                       
      T(1,1) = 1.0D00 - A(M) * A(M)                                     
      IF( M .EQ. 1 )     GO TO 45                                       
      DO  10     I=2,M                                                  
      IM1 = I - 1                                                       
      II = M - IM1                                                      
cxx   10 T(I,1) = - ( A(IM1)+A(M)*A(II) )                                  
      T(I,1) = - ( A(IM1)+A(M)*A(II) ) 
   10 CONTINUE
      IF( M .EQ. 2 )     GO TO 25                                       
C                                                                       
cxx      DO  20     J=2,M1                                                 
      DO  21     J=2,M1
      JM1 = J - 1                                                       
      JJ = M - JM1                                                      
      DO  20     I=J,JJ                                                 
      IM1 = I - 1                                                       
      II = M - IM1                                                      
cxx   20 T(I,J) = T(IM1,JM1) + A(IM1)*A(JM1) - A(II)*A(JJ)                 
      T(I,J) = T(IM1,JM1) + A(IM1)*A(JM1) - A(II)*A(JJ)
   20 CONTINUE
   21 CONTINUE
C                                                                       
cxx   25 DO  30     J=1,M1                                                 
   25 DO  31     J=1,M1 
      JE = MP1 - J                                                      
      JE1 = JE - 1                                                      
      DO  30     I=J,JE1                                                
      II = MP1 - I                                                      
cxx   30 T(JE,II) = T(I,J)                                                 
      T(JE,II) = T(I,J)
   30 CONTINUE
   31 CONTINUE
cxx      DO  40     I=2,M                                                  
      DO  41     I=2,M
      IM1 = I - 1                                                       
      DO  40     J=1,IM1                                                
cxx   40 T(J,I) = T(I,J)                                                   
      T(J,I) = T(I,J)
   40 CONTINUE
   41 CONTINUE
C                                                                       
   45 CONTINUE                                                          
C                                                                       
C         INNOVATION VARIANCE COMPUTATION                               
C                                                                       
      SD = R(1,1)                                                       
      DO  60     J=1,M                                                  
      J1 = J + 1                                                        
      SUM = 0.0D00                                                      
      DO  50     K=1,M                                                  
cxx   50 SUM = SUM + A(K)*R(J1,K+1)                                        
      SUM = SUM + A(K)*R(J1,K+1)
   50 CONTINUE
cxx   60 SD = SD + SUM*A(J)                                                
      SD = SD + SUM*A(J)
   60 CONTINUE
      SUM = 0.0D00                                                      
      DO  70     K=1,M                                                  
cxx   70 SUM = SUM + A(K)*R(1,K+1)                                         
      SUM = SUM + A(K)*R(1,K+1)
   70 CONTINUE
      SD = SD - 2.0D00*SUM                                              
C                                                                       
      DO  80     J=1,M                                                  
      SUM = 0.0D00                                                      
      DO  85     K=1,M                                                  
cxx   85 SUM = SUM + T(J,K)*Z(K)                                           
      SUM = SUM + T(J,K)*Z(K)
   85 CONTINUE
cxx   80 SD = SD + SUM*Z(J)                                                
      SD = SD + SUM*Z(J)
   80 CONTINUE
      SD = SD * DN1                                                     
C                                                    DS     ) 033,6 (ETI
cxx      DO  90     I=1,M                                                  
      DO  91     I=1,M
      DO  90     J=1,M                                                  
cxx   90 U(I,J) = T(I,J)                                                   
      U(I,J) = T(I,J)
   90 CONTINUE
   91 CONTINUE
      CALL  SUBDET( U,DETT,M,MJ )                                       
C                                                                       
C         CHECK THE STATIONARITY                                        
C                                                                       
C                                                                       
      CALL  PARCOR( A,M,B )                                             
C                                                                       
      DO  210     I=1,M                                                 
      IF( DABS(B(I)) .LT. 1.D0 )     GO TO  210                         
      IF( ISW .EQ. 0 )   B(I) = 0.999D0*DABS(B(I))/B(I)                 
      IFG = 1                                                           
  210 CONTINUE                                                          
      IF( IFG .EQ. 0 )     GO TO  220                                   
      IF( ISW .EQ. 1 )     RETURN                                       
C                                                                       
      CALL  ARCOEF( B,M,A )                                             
      GO TO  280                                                        
  220 CONTINUE                                                          
      IF( SD .GT. 0.0D00 ) GO TO 221
      JER = 11111
      RETURN
  221 CONTINUE
C                                                                       
C         LIKELIHOOD                                                    
C                                                                       
      F = 0.5D00 * ( DN*DLOG(SD) - DLOG(DETT) )                         
      FF = F                                                            
      AIC = 2.D0*F + 2*(M+1)                                            
C                                                                       
      IF( ISW .EQ. 1 )     RETURN                                       
C                                                                       
C       COMPUTATION OF GRADIENT OF F                                    
C                                                                       
      DO  1000     JC=1,M                                               
C                                                                       
      S(1,1) = 0.0D00                                                   
      IF( JC .EQ. M )     S(1,1) = -2.0D00*A(M)                         
      IF( M .EQ. 1 )     GO TO 145                                      
C                                                                       
      DO  100     I=2,M                                                 
      II = MP1 - I                                                      
      IM1 = I - 1                                                       
      S(I,1) = 0.0D00                                                   
      IF( II .EQ. JC )     S(I,1) = -A(M)                               
      IF( M .EQ. JC )     S(I,1) = -A(II)                               
      IF( IM1 .EQ. JC )     S(I,1) = S(I,1) - 1.0D00                    
  100 CONTINUE                                                          
      IF( M .EQ. 2 )     GO TO 120                                      
C                                                                       
cxx      DO  110     K=2,M1                                                
      DO  111     K=2,M1
      KM1 = K - 1                                                       
      KK = MP1 - K                                                      
      DO  110     I=K,KK                                                
      IM1 = I - 1                                                       
      II = MP1 - I                                                      
      S(I,K) = S(IM1,KM1)                                               
      IF( KM1 .EQ. JC )     S(I,K) = S(I,K) + A(IM1)                    
      IF( IM1 .EQ. JC )     S(I,K) = S(I,K) + A(KM1)                    
      IF( II .EQ. JC )     S(I,K) = S(I,K) - A(KK)                      
      IF( KK .EQ. JC )     S(I,K) = S(I,K) - A(II)                      
  110 CONTINUE
  111 CONTINUE                                                          
C                                                                       
cxx  120 DO  130     J=1,M1                                                
  120 DO  131     J=1,M1
      JE=  MP1 - J                                                      
      JE1 = JE - 1                                                      
      DO  130     I=J,JE1                                               
      II = MP1 - I                                                      
cxx  130 S(JE,II) = S(I,J)                                                 
      S(JE,II) = S(I,J)
  130 CONTINUE
  131 CONTINUE
cxx      DO  140     I=2,M                                                 
      DO  141     I=2,M
      IM1 = I - 1                                                       
      DO  140     J=1,IM1                                               
cxx  140 S(J,I) = S(I,J)                                                   
      S(J,I) = S(I,J)
  140 CONTINUE
  141 CONTINUE
  145 CONTINUE                                                          
C                                                                       
C                                                                       
      J1 = JC + 1                                                       
      DSD = R(1,J1)                                                     
      DO  150     K=1,M                                                 
cxx  150 DSD = DSD - A(K) * R(J1,K+1)                                      
      DSD = DSD - A(K) * R(J1,K+1)
  150 CONTINUE
C                                                                       
      DSD = -2.0D00*DSD                                                 
      DO  170     I=1,M                                                 
      SUM = 0.0D00                                                      
      DO  160     J=1,M                                                 
cxx  160 SUM = SUM + S(I,J)*Z(J)
      SUM = SUM + S(I,J)*Z(J)
  160 CONTINUE                                           
cxx  170 DSD = DSD + SUM*Z(I)
      DSD = DSD + SUM*Z(I)
  170 CONTINUE                                              
      DSD = DSD * DN1                                                   
      USUM = 0.0D00                                                     
      DO  200     K=1,M                                                 
cxx      DO  180     I=1,M
      DO  181     I=1,M                                                 
      DO  180     J=1,M                                                 
cxx  180 U(I,J) = T(I,J)
      U(I,J) = T(I,J)
  180 CONTINUE
  181 CONTINUE                                                   
      DO  190     J=1,M                                                 
cxx  190 U(K,J) = S(K,J)                                                   
      U(K,J) = S(K,J)
  190 CONTINUE
C                                                                       
      CALL  SUBDET( U,UDET,M,MJ )                                       
C                                                                       
cxx  200 USUM = USUM + UDET
      USUM = USUM + UDET
  200 CONTINUE                                                
C                                                                       
C         GRADIENT OF LIKELIHOOD FUNCTION                               
C                                                                       
      G(JC) = 0.5D00*( -DN*DSD/SD + USUM/DETT )                         
 1000 CONTINUE                                                          
      DO  300     I=1,M                                                 
cxx  300 G(I) = -G(I)
      G(I) = -G(I)
  300 CONTINUE                                                      
      RETURN                                                            
      END   
cc      SUBROUTINE  HESIAN( X,K,H )                                       
      SUBROUTINE  HESIAN( X,K,N,R,SD,H )
C                                                                       
C     THIS SUBROUTINE RETURNS THE INVERSE OF AN APPROXIMATION TO THE HES
C     OF LOG-LIKELIHOOD FUNCTION OF THE AUTOREGRESSIVE MODEL OF ORDER K.
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
C             INVDET                                                    
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS:                                                         
C          X:     VECTOR OF AR-COEFFICIENTS                             
C          K:     ORDER OF THE AR-MODEL                                 
C          N:     DATA LENGTH                                           
C          Z:     VECTOR OF ORIGINAL DATA                               
C          R:     SQUARE ROOT OF THE COVARIANCE MATRIX, OUTPUT OF SUBROU
C                 ARMLE                                                 
C          SD:    ESTIMATE OF INNOVATION VARIANCE                       
C                                                                       
C       OUTPUT:                                                         
C          H:     INVERSE OF THE HESSIAN                                
cxx      IMPLICIT  REAL*8( A-H,O-Z )                                       
CC      REAL*4  Z                                                         
cc      DIMENSION  Z(10000)
CC      COMMON     / AAA /  N , Z(1)                                      
cc      COMMON     / AAA /  N
cc      COMMON     / BBB /  Z
cc      COMMON     / DDD /  R , F , AIC , SD                              
cc      DIMENSION  R(31,31)                                               
cx      DIMENSION  X(1)                                                   
cc      DIMENSION  H(50,50) , S(50)                                       
cxx      DIMENSION  R(K+1,K+1)                                               
cxx      DIMENSION  X(K)                                                   
cxx      DIMENSION  H(K,K) , S(K)
      INTEGER :: K, N
      REAL(8) :: X(K), R(K+1,K+1), SD, H(K,K)
      REAL(8) :: S(K), SUM, HDET
      DO 10  I=1,K                                                      
      SUM = R(1,I+1)                                                    
      DO 20  II=1,K                                                     
cxx   20 SUM = SUM - X(II)*R(II+1,I+1)                                     
      SUM = SUM - X(II)*R(II+1,I+1)
   20 CONTINUE                                     
cxx   10 S(I) = SUM / SD
      S(I) = SUM / SD                                                   
   10 CONTINUE
cxx      DO 30  I=1,K                                                      
      DO 31  I=1,K 
      DO 30  J=1,K                                                      
cxx   30 H(I,J) = 0.5D0 * (R(I+1,J+1)/SD - S(I)*S(J)/N)                    
      H(I,J) = 0.5D0 * (R(I+1,J+1)/SD - S(I)*S(J)/N)
   30 CONTINUE
   31 CONTINUE
cc      CALL  INVDET( H,HDET,K,50 )                                       
      CALL  INVDET( H,HDET,K,K )                                       
      RETURN                                                            
      E N D                                                             
cc      SUBROUTINE  LINEAR( FUNCT,X,H,RAM,EE,K,IG )                       
cxx      SUBROUTINE  LINEAR( FUNCT,Z,N,X,H,RAM,EE,K,R,AIC,SD,F,ISW,IPR,IG,
ccx     *                    JFG,LU )
cx     *                    JFG,LU,JER )
cxx     *                    
      SUBROUTINE  LINEAR( FUNCT,Z,N,X,H,RAM,EE,K,R,AIC,SD,F,ISW,IG,JER )
C                                                                       
C     THIS SUBROUTINE PERFORMS THE LINEAR SEARCH ALONG THE DIRECTION SPE
C     BY THE VECTOR H                                                   
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
C             FUNCT                                                     
C       ----------------------------------------------------------------
C                                                                       
C     INPUTS:                                                           
C        FUNCT:   EXTERNAL FUNCTION SPECIFICATION                       
C        X:       VECTOR OF POSITION                                    
C        H:       SEARCH DIRECTION                                      
C        K:       DIMENSION OF VECTOR X                                 
C                                                                       
C     OUTPUTS:                                                          
C        RAM:     OPTIMAL STEP WIDTH                                    
C        E2:      MINIMUM FUNCTION VALUE                                
C        IG:      ERROR CODE                                            
C                                                                       
cxx      IMPLICIT  REAL  *8 ( A-H,O-Z )                                    
cxx      INTEGER  RETURN,SUB                                               
cc      DIMENSION  X(1) , H(1) , X1(50)                                   
cc      DIMENSION  G(50)                                                  
cc      COMMON     / CCC /  ISW , IPR                                     
cc      DIMENSION  R(31,31)
cxx      DIMENSION  X(K) , H(K) , X1(K)                                   
cxx      DIMENSION  G(K)                                                  
cxx      DIMENSION  Z(N)
cxx      DIMENSION  R(K+1,K+1)

      INTEGER :: N, K, ISW, IG, JER
      REAL(8) :: Z(N), X(K), H(K), RAM, EE, R(K+1,K+1), AIC, SD, F
      INTEGER :: RETURN, SUB
      REAL(8) :: X1(K), G(K), CONST2, HNORM, E1, E2, E3,
     1           RAM1, RAM2, RAM3, A1, A2, A3, B1, B2
C                                                                       
      ISW = 1                                                           
      IF( RAM .LE. 1.0D-30 )  RAM = 0.01D0                              
      CONST2 = 1.0D-60                                                  
      HNORM = 0.D0                                                      
      DO 10  I=1,K                                                      
cxx   10 HNORM = HNORM + H(I)**2                                           
      HNORM = HNORM + H(I)**2
   10 CONTINUE
      HNORM = DSQRT( HNORM )                                            
C                                                                       
      RAM2 = RAM                                                        
      E1 =EE                                                            
      RAM1 = 0.D0                                                       
C                                                                       
      DO 20  I=1,K                                                      
cxx   20 X1(I) = X(I) + RAM2*H(I)                                          
      X1(I) = X(I) + RAM2*H(I)
   20 CONTINUE
cc      CALL  FUNCT( K,X1,E2,G,IG )                                       
cc      IF(IPR.GE.7)  WRITE(6,2)  RAM2,E2                                 
ccx      CALL  FUNCT( Z,N,K,X1,R,ISW,E2,G,AIC,SD,F,IG )
      CALL  FUNCT( Z,N,K,X1,R,ISW,E2,G,AIC,SD,F,IG,JER )
      IF( JER .NE. 0 ) RETURN
cx      IF( (IPR.GE.7) .AND. (JFG.NE.0) )  WRITE(LU,2)  RAM2,E2
C                                                                       
      IF( IG .EQ. 1 )  GO TO  50                                        
      IF( E2 .GT. E1 )  GO TO 50                                        
   30 RAM3 = RAM2*2.D0                                                  
      DO 40  I=1,K                                                      
cxx   40 X1(I) = X(I) + RAM3*H(I)                                          
      X1(I) = X(I) + RAM3*H(I)
   40 CONTINUE
cc      CALL  FUNCT( K,X1,E3,G,IG )                                       
ccx      CALL  FUNCT( Z,N,K,X1,R,ISW,E3,G,AIC,SD,F,IG )
      CALL  FUNCT( Z,N,K,X1,R,ISW,E3,G,AIC,SD,F,IG,JER )
      IF( JER .NE. 0 ) RETURN
      IF( IG.EQ.1 )  GO TO  500                                         
cc      IF( IPR.GE.7 )  WRITE(6,3)  RAM3,E3                               
cx      IF( (IPR.GE.7) .AND. (JFG.NE.0 ) )  WRITE(LU,3)  RAM3,E3
      IF( E3 .GT. E2 )  GO TO 70                                        
      RAM1 = RAM2                                                       
      RAM2 = RAM3                                                       
      E1 = E2                                                           
      E2 = E3                                                           
      GO TO 30                                                          
C                                                                       
   50 RAM3 = RAM2                                                       
      E3 = E2                                                           
      RAM2 = RAM3*0.1D0                                                 
      IF( RAM2*HNORM .LT. CONST2 )  GO TO  400                          
      DO 60  I=1,K                                                      
cxx   60 X1(I) = X(I) + RAM2*H(I)
      X1(I) = X(I) + RAM2*H(I)
   60 CONTINUE                                          
cc      CALL  FUNCT( K,X1,E2,G,IG )                                       
cc      IF(IPR.GE.7)  WRITE(6,4)  RAM2,E2                                 
ccx      CALL  FUNCT( Z,N,K,X1,R,ISW,E2,G,AIC,SD,F,IG )
      CALL  FUNCT( Z,N,K,X1,R,ISW,E2,G,AIC,SD,F,IG,JER )
      IF( JER .NE. 0 ) RETURN
cx      IF( (IPR.GE.7) .AND. (JFG.NE.0) )  WRITE(LU,4)  RAM2,E2
      IF( E2.GT.E1 )  GO TO 50                                          
C                                                                       
cc   70 ASSIGN 80 TO RETURN                                               
   70 RETURN = 80
      GO TO 200                                                         
C                                                                       
   80 DO 90  I=1,K                                                      
cxx   90 X1(I) = X(I) + RAM*H(I)                                           
      X1(I) = X(I) + RAM*H(I) 
   90 CONTINUE
cc      CALL  FUNCT( K,X1,EE,G,IG )                                       
cc      IF(IPR.GE.7)  WRITE(6,5)  RAM,EE                                  
ccx      CALL  FUNCT( Z,N,K,X1,R,ISW,EE,G,AIC,SD,F,IG )
      CALL  FUNCT( Z,N,K,X1,R,ISW,EE,G,AIC,SD,F,IG,JER )
      IF( JER .NE. 0 ) RETURN
cx      IF( (IPR.GE.7) .AND. (JFG.NE.0) )  WRITE(LU,5)  RAM,EE
C                                                                       
      IFG = 0                                                           
cc      ASSIGN  300 TO  SUB                                               
cc      ASSIGN 200 TO SUB                                                 
cc   95 ASSIGN 130 TO RETURN                                              
      SUB = 300
      SUB = 200
   95 RETURN = 130
      IF( RAM .GT. RAM2 )  GO TO 110                                    
      IF( EE .GE. E2 )  GO TO 100                                       
      RAM3 = RAM2                                                       
      RAM2 = RAM                                                        
      E3 =E2                                                            
      E2 =EE                                                            
cc      GO TO  SUB,( 200,300 )                                            
      IF( SUB .EQ. 200 ) GO TO  200
      IF( SUB .EQ. 300 ) GO TO  300
C                                                                       
  100 RAM1 = RAM                                                        
      E1 = EE                                                           
cc      GO TO  SUB,( 200,300 )                                            
      IF( SUB .EQ. 200 ) GO TO  200
      IF( SUB .EQ. 300 ) GO TO  300
C                                                                       
  110 IF( EE .LE. E2 )  GO TO 120                                       
      RAM3 = RAM                                                        
      E3 = EE                                                           
cc      GO TO  SUB,( 200,300 )                                            
      IF( SUB .EQ. 200 ) GO TO  200
      IF( SUB .EQ. 300 ) GO TO  300
C                                                                       
  120 RAM1 = RAM2                                                       
      RAM2 = RAM                                                        
      E1 = E2                                                           
      E2 = EE                                                           
cc      GO TO  SUB,( 200,300 )                                            
      IF( SUB .EQ. 200 ) GO TO  200
      IF( SUB .EQ. 300 ) GO TO  300
C                                                                       
  130 DO 140  I=1,K                                                     
cxx  140 X1(I) = X(I) + RAM*H(I)                                           
      X1(I) = X(I) + RAM*H(I)
  140 CONTINUE
cc      CALL  FUNCT( K,X1,EE,G,IG )                                       
cc      IF( IPR.GE.7 )  WRITE(6,6)  RAM,EE                                
ccx      CALL  FUNCT( Z,N,K,X1,R,ISW,EE,G,AIC,SD,F,IG )
      CALL  FUNCT( Z,N,K,X1,R,ISW,EE,G,AIC,SD,F,IG,JER )
      IF ( JER .NE. 0 ) RETURN
cx      IF( (IPR.GE.7)  .AND. (JFG.NE.0) )  WRITE(LU,6)  RAM,EE
cc      ASSIGN 200 TO SUB                                                 
      SUB = 200
      IFG = IFG+1                                                       
      IFG = 0                                                           
      IF( IFG .EQ. 1 )  GO TO 95                                        
C                                                                       
      IF( E2 .LT. EE )  RAM = RAM2                                      
      RETURN                                                            
C                                                                       
C      -------  INTERNAL SUBROUTINE SUB1  -------                       
  200 A1 = (RAM3-RAM2)*E1                                               
      A2 = (RAM1-RAM3)*E2                                               
      A3 = (RAM2-RAM1)*E3                                               
      B2 = (A1+A2+A3)*2.D0                                              
      B1 = A1*(RAM3+RAM2) + A2*(RAM1+RAM3) + A3*(RAM2+RAM1)             
      IF( B2 .EQ. 0.D0 )  GO TO 210                                     
      RAM = B1 /B2                                                      
cc      GO TO RETURN ,( 80,130 )                                          
      IF( RETURN .EQ. 80 ) GO TO 80
      IF( RETURN .EQ. 130 ) GO TO 130
C                                                                       
  210 IG = 1                                                            
      RAM = RAM2                                                        
      RETURN                                                            
C                                                                       
C      -------  INTERNAL SUBROUTINE SUB2  -------                       
C                                                                       
  300 IF( RAM3-RAM2 .GT. RAM2-RAM1 )  GO TO 310                         
      RAM = (RAM1+RAM2)*0.5D0                                           
cc      GO TO RETURN ,( 80,130 )                                          
      IF( RETURN .EQ. 80 ) GO TO 80
      IF( RETURN .EQ. 130 ) GO TO 130
C                                                                       
  310 RAM = (RAM2+RAM3)*0.5D0                                           
cc      GO TO RETURN ,( 80,130 )                                          
      IF( RETURN .EQ. 80 ) GO TO 80
      IF( RETURN .EQ. 130 ) GO TO 130
C ------------------------------------------------------------          
C                                                                       
  400 RAM = 0.D0                                                        
      RETURN                                                            
C ------------------------------------------------------------          
C                                                                       
  500 RAM = (RAM2+RAM3)*0.5D0                                           
  510 DO 520  I=1,K                                                     
cxx  520 X1(I) = X(I) + RAM*H(I)
      X1(I) = X(I) + RAM*H(I)
  520 CONTINUE                                           
cc      CALL  FUNCT( K,X1,E3,G,IG )                                       
cc      IF( IPR.GE.7 )  WRITE(6,7)  RAM,E3                                
ccx      CALL  FUNCT( Z,N,K,X1,R,ISW,E3,G,AIC,SD,F,IG )
      CALL  FUNCT( Z,N,K,X1,R,ISW,E3,G,AIC,SD,F,IG,JER )
      IF( JER .NE. 0 ) RETURN
cx      IF( (IPR.GE.7) .AND. (JFG.NE.0) )  WRITE(LU,7)  RAM,E3
      IF( IG.EQ.1 )  GO TO 540                                          
      IF( E3.GT.E2 )  GO TO 530                                         
      RAM1 = RAM2                                                       
      RAM2 = RAM                                                        
      E1 = E2                                                           
      E2 = E3                                                           
      GO TO 500                                                         
C                                                                       
  530 RAM3 = RAM                                                        
      GO TO 70                                                          
C                                                                       
  540 RAM = (RAM2+RAM)*0.5D0                                            
      GO TO 510                                                         
C                                                                       
C ------------------------------------------------------------          
cxx    1 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E1 =',D25.17 )                
cxx    2 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E2 =',D25.17 )                
cxx    3 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E3 =',D25.17 )                
cxx    4 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E4 =',D25.17 )                
cxx    5 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E5 =',D25.17 )                
cxx    6 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E6 =',D25.17 )                
cxx    7 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E7 =',D25.17 )                
      E N D                                                             
