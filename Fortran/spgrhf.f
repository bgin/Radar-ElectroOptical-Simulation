cc      SUBROUTINE SPGRH(Y,N,MODE)
      SUBROUTINE SPGRHF(Y,N,LAGH1,IFPL1,MODE,PERIOD,CXX,CN,XMEAN,SD,AIC,
     *                   PARCOR,PXX,IER)
c
      INCLUDE 'timsac_f.h'
cc      !DEC$ ATTRIBUTES DLLEXPORT::SPGRHF
c
cxx      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      INTEGER H,H1                                                      
cxx      INTEGER PERIOD
CC      REAL*4 X,SXX,BMYA,BMYM                                            
cc      DIMENSION Y(1)                                                    
cc      DIMENSION X(4000),CXX(1001)                                       
cc      DIMENSION A(101),B(101)                                           
cc      DIMENSION SXX(501)                                                
cxx      DIMENSION Y(N)
cxx      DIMENSION CXX(LAGH1),CN(LAGH1)
cxx      DIMENSION A(IFPL1),B(IFPL1)
cxx      DIMENSION PXX(LAGH1)
cxx      DIMENSION SD(IFPL1),AIC(IFPL1),PARCOR(IFPL1-1)
      INTEGER :: N, LAGH1, IFPL1, MODE, PERIOD, IER
      REAL(8) :: Y(N), CXX(LAGH1), CN(LAGH1), XMEAN, SD(IFPL1),
     1           AIC(IFPL1), PARCOR(IFPL1-1), PXX(LAGH1)
      REAL(8) :: A(IFPL1), B(IFPL1), SGME2, OAIC
C
cc      H=60                                                              
cc      H1=H+1                                                            
C                                                                       
cc      DO 10 I=1,N                                                       
cc      X(I)=Y(I)                                                         
cc   10 CONTINUE                                                          
C                                                                       
C     AUTO COVARIANCE COMPUTATION.                                      
cc      LAGH=N-1                                                          
cc      LAGH=MIN0(LAGH,60)                                                
cc      LAGH1=LAGH+1                                                      
cc      CALL SAUTCO(X,CXX,N,LAGH1)
      CALL AUTCORF(Y,N,CXX,CN,LAGH1,XMEAN)
cc      AN=N                                                              
cc      IFPL=3.0D-00*DSQRT(AN)                                            
cc      IFPL=MIN0(IFPL,50,LAGH)                                           
cc      IFPL=MIN0(30,N-1)                                                 
cc      IFPL1=IFPL+1                                                      
cc      CALL SICP2(CXX,IFPL1,N,A,L,SGME2,OAIC)                            
      CALL SICP2(CXX,IFPL1,N,A,L,SGME2,OAIC,SD,AIC,PARCOR,IER)
      IF(MODE .EQ. 0) RETURN                                            
cc      WRITE(6,50)                                                       
      K=0                                                               
cc 1031 CALL SNRASP(A,B,SXX,SGME2,L,K,H1,BMYA,BMYM)                       
cxx 1031 CALL SNRASP(A,B,PXX,SGME2,L,K,LAGH1,PERIOD)
      CALL SNRASP(A,B,PXX,SGME2,L,K,LAGH1,PERIOD)
C
      RETURN                                                            
cxx   50 FORMAT(1H ,'*** HIGH ORDER AR-SPECTRUM AS AN APPROXIMATION TO PERI
cxx     *ODGRAM (ORDER IS FIXED AT 30 ) ***')                              
      END                                                               
cc      SUBROUTINE SICP2(CYY,L1,N,COEF,MO,OSD,OAIC)                       
      SUBROUTINE SICP2(CYY,L1,N,COEF,MO,OSD,OAIC,SD1,AIC1,PARCOR,IER)
C     COMMON SUBROUTINE                                                 
C     THIS SUBROUTINE FITS AUTOREGRESSIVE MODELS OF SUCCESSIVELY        
C     INCREASING ORDER UP TO L(=L1-1).                                  
C     INPUT:                                                            
C     CYY(I),I=0,L1; AUTOCOVARIANCE SEQUENCE                            
C     L1: L1=L+1, L IS THE UPPER LIMIT OF THE MODEL ORDER               
C     N; LENGTH OF ORIGINAL DATA                                        
C     OUT PUT:                                                          
C     COEF; AR-COEFFICIENTS                                             
C     MO: ORDER OF AR                                                   
C     OSD: INNOVATION VARIANCE                                          
C     OAIC: VALUE OF AIC                                                
cxx      IMPLICIT REAL*8(A-H,O-Z)                                          
cxx      DIMENSION CYY(L1),COEF(L1)                                        
cc      DIMENSION A(101),B(101)                                           
cxx      DIMENSION A(L1),B(L1)
cxx      DIMENSION SD1(L1),AIC1(L1),PARCOR(L1-1)
      INTEGER :: L1, N, MO, IER
      REAL(8) :: CYY(L1), COEF(L1), OSD, OAIC, SD1(L1), AIC1(L1),
     1           PARCOR(L1-1)
      REAL(8) :: A(L1), B(L1), CST0, CST1, CST2, CST20, CST05, CST01,
     1           SD, AN, RAN, SCALH, AIC, SE, SDR, D, D2, AM, ANFC
cc      REAL*4 AX,BL,STA,DASH,PLUS                                        
cc      REAL*4 FFFF                                                       
cc      REAL*4  F(41) / 41*1H  /, AMES(41) / 41*1H- /                     
cc      DATA AX,BL,STA,DASH,PLUS/1H!,1H ,1H*,1H-,1H+/                     
c-----
      IER=0
c-----
      CST0=0.0D-00                                                      
      CST1=1.0D-00                                                      
      CST2=2.0D-00                                                      
      CST20=20.0D-00                                                    
      CST05=0.05D-00                                                    
      CST01=0.00001D-00                                                 
      L=L1-1                                                            
      SD=CYY(1)                                                         
      AN=N                                                              
      OAIC=AN*DLOG(SD)                                                  
      OSD=SD                                                            
      MO=0                                                              
C     INITIAL CONDITION PRINT OUT                                       
cc  991 CONTINUE                                                          
C     WRITE(6,1100)                                                     
      RAN=CST1/DSQRT(AN)                                                
      SCALH=CST20                                                       
cxx      JJ0=SCALH+CST1                                                    
cxx      JJL=SCALH*CST2+CST1                                               
      JJ0=INT(SCALH+CST1)
      JJL=INT(SCALH*CST2+CST1)
cc      AMES(1)=PLUS                                                      
cc      AMES(11)=PLUS                                                     
cc      AMES(JJ0)=PLUS                                                    
cc      AMES(JJ0+10)=PLUS                                                 
cc      AMES(JJL)=PLUS                                                    
cxx      IAN=SCALH*(RAN+CST05)                                             
      IAN=INT(SCALH*(RAN+CST05))
      IAN1=IAN+JJ0                                                      
      IAN2=2*IAN+JJ0                                                    
      LAN1=-IAN+JJ0                                                     
      LAN2=-2*IAN+JJ0                                                   
C     WRITE(6,26100)                                                    
cc      WRITE(6,26101)                                                    
cc      WRITE(6,261)                                                      
cc      WRITE(6,26102)                                                    
cc      WRITE(6,262)                                                      
cc      WRITE(6,264) (AMES(J),J=1,JJL)                                    
cc      WRITE(6,859) MO,OSD,OAIC
      SD1(MO+1)=OSD
      AIC1(MO+1)=OAIC
c-----
      AIC=OAIC
c-----
cc      F(JJ0)=AX                                                         
cc      F(IAN1)=AX                                                        
cc      F(IAN2)=AX                                                        
cc      F(LAN1)=AX                                                        
cc      F(LAN2)=AX                                                        
cc      WRITE(6,861) (F(J),J=1,JJL)                                       
      SE=CYY(2)                                                         
C     ITERATION START                                                   
      DO 400 M=1,L                                                      
      SDR=SD/CYY(1)                                                     
      IF(SDR.GE.CST01) GO TO 399                                        
cc      WRITE(6,2600)                                                     
c-----
      IER=2600
c-----
      GO TO 402
  399 MP1=M+1                                                           
      D=SE/SD                                                           
      A(M)=D                                                            
      D2=D*D                                                            
      SD=(CST1-D2)*SD                                                   
      AM=M                                                              
      AIC=AN*DLOG(SD)+CST2*AM                                           
      IF(M.EQ.1) GO TO 410                                              
C     A(I) COMPUTATION                                                  
      LM=M-1                                                            
      DO 420 I=1,LM                                                     
      A(I)=A(I)-D*B(I)                                                  
  420 CONTINUE                                                          
  410 CONTINUE                                                          
      DO 421 I=1,M                                                      
      IM=MP1-I                                                          
cxx  421 B(I)=A(IM)
      B(I)=A(IM)
  421 CONTINUE
C     M,SD,AIC  PRINT OUT                                               
      IF(A(M).LT.CST0) GO TO 300                                        
cxx      NFC=SCALH*(A(M)+CST05)                                            
      NFC=INT(SCALH*(A(M)+CST05))
      GO TO 310                                                         
cxx  300 NFC=SCALH*(A(M)-CST05)                                            
  300 NFC=INT(SCALH*(A(M)-CST05))
  310 ANFC=NFC                                                          
cxx      JJ=ANFC+SCALH+CST1                                                
      JJ=INT(ANFC+SCALH+CST1)                                                
cc      FFFF=F(JJ)                                                        
cc      F(JJ)=STA                                                         
cc      WRITE(6,860) M,SD,AIC,A(M),(F(J),J=1,JJL)                         
      SD1(M+1)=SD
      AIC1(M+1)=AIC
      PARCOR(M)=A(M)
cc      F(JJ)=FFFF                                                        
C                                                                       
cxx  990 IF(OAIC.LT.AIC) GO TO 440                                         
      IF(OAIC.LT.AIC) GO TO 440
      OAIC=AIC                                                          
      OSD=SD                                                            
      MO=M                                                              
      GO TO 440                                                         
C     DO 430 I=1,M                                                      
C 430 COEF(I)=-A(I)                                                     
  440 IF(M.EQ.L) GO TO 400                                              
      SE=CYY(M+2)                                                       
      DO 441 I=1,M                                                      
cxx  441 SE=SE-B(I)*CYY(I+1)                                               
      SE=SE-B(I)*CYY(I+1)
  441 CONTINUE
  400 CONTINUE
  402 CONTINUE
cc      WRITE(6,870) OAIC,MO                                              
      OAIC=AIC                                                         
      OSD=SD                                                            
      MO=L                                                              
      DO 5100 I=1,L                                                     
cxx 5100 COEF(I)=-A(I)
      COEF(I)=-A(I)
 5100 CONTINUE
C                                                                       
C     MO, COEF(I) OUT PUT                                               
C     WRITE(6,1871)                                                     
C     WRITE(6,871)                                                      
C     CALL SUBVCP(COEF,MO)                                              
cc  699 F(JJ0)=BL                                                         
cc      F(IAN1)=BL                                                        
cc      F(IAN2)=BL                                                        
cc      F(LAN1)=BL                                                        
cc      F(LAN2)=BL                                                        
cc      AMES(JJ0)=DASH                                                    
cc      AMES(JJ0+10)=DASH                                                 
cc      AMES(JJL)=DASH                                                    
      RETURN                                                            
cxx26100 FORMAT(1H ,16X,'SD(M)',15X,'AIC(M)',13X,'A(M)')                   
cxx26101 FORMAT(1H ,////,4X,'M',11X,'INNOVATION',10X,'AIC(M)=    ',8X,     
cxx     A      'PARTIAL AUTO-')                                            
cxx  261 FORMAT(1H ,16X,'   VARIANCE',9X,'N*DLOG(SD(M))+2*M',              
cxx     A      '  CORRELATION    ',                                        
cxx     A      'PARTIAL CORRELATION (LINES SHOW +SD AND +2SD)')
cxx26102 FORMAT(1H+,102X,'_',7X,'_')                                       
cxx  262 FORMAT(1H ,69X,'-1',19X,'0',19X,'1')                              
cxx  264 FORMAT(1H ,70X,41A1)                                              
cxx  859 FORMAT(1H ,I5,2X,2D20.5)                                          
cxx  860 FORMAT(1H ,I5,2X,3D20.5,3X,41A1)                                  
cxx  861 FORMAT(1H+,70X,41A1)                                              
cxx  960 FORMAT(1H ,5X,'A(I)')                                             
cxx  870 FORMAT(/,1H ,'MINIMUM AIC(M)=',D12.5,2X,'ATTAINED AT M=',I5/)     
cxx  980 FORMAT(1H ,5X,'COEF(I)')                                          
cxx 1100 FORMAT(1H ,'AIC(M)=N*DLOG(SD)+2*M')                               
cxx  871 FORMAT(1H ,'AR-COEFFICIENTS')                                     
cxx 1871 FORMAT(1H ,'AR MODEL: Y(N)+AR(1)Y(N-1)+...+AR(M)Y(N-M)=X(N)')     
cxx 2600 FORMAT(1H ,'ACCURACY OF COMPUTATION LOST')                        
      END                                                               
cc      SUBROUTINE SNRASP(A,B,SXX,SGME2,L,K,H1,BMYA,BMYM)                 
      SUBROUTINE SNRASP(A,B,PXX,SGME2,L,K,H1,IPPP)
C     THIS PROGRAM COMPUTES POWER SPECTRUM OF AN AR-MA PROCESS DEFINED B
C     X(N)=A(1)X(N-1)+...+A(L)X(N-L)+E(N)+B(1)E(N-1)+...+B(K)E(N-K),    
C     WHERE E(N) IS A WHITE NOISE WITH ZERO MEAN AND VARIANCE EQUAL TO  
C     SGME2.  OUTPUTS PXX(I) ARE GIVEN AT FREQUENCIES I/(2*H)           
C     I=0,1,...,H.                                                      
C     REQUIRED INPUTS ARE:                                              
C     L,K,H,SGME2,(A(I),I=1,L), AND (B(I),I=1,K).                       
C     0 IS ALLOWABLE AS L AND/OR K.                                     
cxx      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      INTEGER H,H1                                                      
cxx      INTEGER H1                                                      
CC      REAL*4 SXX,BMYA,BMYM,FY
cc      REAL*4 FY                                           
cc      INTEGER F,AX,BL,STA                                               
cc      DIMENSION F(61),FY(10)                                            
cc      DIMENSION A(101),B(101)                                           
cc      DIMENSION G(501),GR1(501),GI1(501),GR2(501),GI2(501)              
cc      DIMENSION PXX(501),SXX(501)                                       
cxx      DIMENSION A(L),B(K)                                           
cxx      DIMENSION G(L+K+1),GR1(H1),GI1(H1),GR2(H1),GI2(H1)
cxx      DIMENSION PXX(H1)
      INTEGER :: L, K, H1, IPPP
      REAL(8) :: A(L), B(K), PXX(H1), SGME2
      REAL(8) :: G(L+K+1), GR1(H1), GI1(H1), GR2(H1), GI2(H1),
     1           CST0, CST1, T0
cc      COMMON /ILOGT/IDUMMY(2),IPUNC                                     
cc      COMMON /IDATA/IPPP                                                
cc      DATA AX/1H!/                                                      
cc      DATA BL,STA/1H ,1H*/                                              
      IPPP0=120/IPPP                                                    
      CST0=0.0D-00                                                      
      CST1=1.0D-00
cxx  310 CONTINUE                                                          
      IF(L.LE.0) GO TO 320                                              
      DO 330 I=1,L                                                      
cxx  330 A(I)=-A(I)                                                        
      A(I)=-A(I)
  330 CONTINUE
  320 L1=L+1                                                            
      K1=K+1                                                            
      G(1)=CST1                                                         
      IF(L.LE.0) GO TO 400                                              
      DO 10 I=1,L                                                       
      I1=I+1                                                            
cxx   10 G(I1)=-A(I)
      G(I1)=-A(I)
   10 CONTINUE
C     COMMON SUBROUTINE CALL                                            
  400 CALL FOUGER(G,L1,GR1,GI1,H1)                                      
      G(1)=CST1                                                         
      IF(K.LE.0) GO TO 410                                              
      DO 20 I=1,K                                                       
      I1=I+1                                                            
cxx   20 G(I1)=B(I)                                                        
      G(I1)=B(I)
   20 CONTINUE
C     COMMON SUBROUTINE CALL                                            
  410 CALL FOUGER(G,K1,GR2,GI2,H1)                                      
      DO 30 I=1,H1                                                      
cxx   30 PXX(I)=(GR2(I)**2+GI2(I)**2)/(GR1(I)**2+GI1(I)**2)*SGME2
      PXX(I)=(GR2(I)**2+GI2(I)**2)/(GR1(I)**2+GI1(I)**2)*SGME2
   30 CONTINUE
C     WRITE(6,60)                                                       
C     WRITE(6,160)                                                      
C     WRITE(6,61) L,K,H                                                 
C     WRITE(6,164) SGME2                                                
C     WRITE(6,264)                                                      
C     WRITE(6,265)                                                      
      IF(L.LE.0) GO TO 500                                              
      DO 340 I=1,L                                                      
cxx  340 A(I)=-A(I)
      A(I)=-A(I)
  340 CONTINUE
C     WRITE(6,62)                                                       
C     COMMON SUBROUTINE CALL                                            
C     CALL SUBVCP(A,L)                                                  
  500 IF(K.LE.0) GO TO 510                                              
cc      WRITE(6,63)                                                       
C     COMMON SUBROUTINE CALL                                            
cc      CALL SUBVCP(B,K)                                                  
  510 T0=CST0                                                           
cc      DO 520 I=1,H1                                                     
cc      T=PXX(I)                                                          
cc      IF(T.LT.T0) T=-T                                                  
cc  520 SXX(I)=DLOG10(T)                                                  
C     SEARCH  FOR THE MAXIMUM AND MINIMUM OF (SXX(I),I=0,LAGH)          
cc      CALL DSP3(SXX,H1,BMYA,BMYM)                                       
cc      WRITE(6,266)                                                      
cc  266 FORMAT(1H ,47X,'(DENSITY IN DB)')                                 
cc      FY(1)=BMYM*10.0                                                   
cc      DO 2031 I=2,7                                                     
cc 2031 FY(I)=FY(I-1)+10.0                                                
cc      WRITE(6,2033) (FY(I),I=1,7)                                       
cc 2033 FORMAT(1H ,25X,6(F5.1,5X),F5.1)                                   
cc      WRITE(6,67)                                                       
cc   67 FORMAT(1H ,4X,'I',1X,'RATIONAL SPECTRUM')                         
cc      WRITE(6,2035)                                                     
cc 2035 FORMAT(1H+,28X,6('+',9('-')),'+')                                 
cc      DO 2050 J=1,61                                                    
cc 2050 F(J)=BL                                                           
cc      DO 2040 I=1,H1                                                    
cc      IT=(SXX(I)-BMYM)*10.0+1.0                                         
cc      IM1=I-1                                                           
cc      F(1)=AX                                                           
cc      F(IT)=STA                                                         
cc      WRITE(6,2060) IM1,PXX(I),(F(J),J=1,61)                            
cc      MODP=MOD(IM1,IPPP0)                                               
cc      IF(MODP .EQ. 0) WRITE(6,1600)                                     
cc 1600 FORMAT(1H+,90X,'----X')                                           
cc      IF(IM1 .EQ. 5) WRITE(6,1601)                                      
cc 1601 FORMAT(1H+,98X,'X: IF PEAKS(TROUGHS) APPEAR')                     
cc      IF(IM1 .EQ. 6) WRITE(6,1602)                                      
cc 1602 FORMAT(1H+,101X,'AT THESE FRIQUENCIES,')                          
cc      IF(IM1 .EQ. 7) WRITE(6,1603)                                      
cc 1603 FORMAT(1H+,101X,'TRY LOWER(HIGHER) VALUES')                       
cc      IF(IM1 .EQ. 8) WRITE(6,1604)                                      
cc 1604 FORMAT(1H+,101X,'OF RIGID AND WATCH ABIC')                        
cc      IF(IM1 .EQ. 42 .AND. IPPP .EQ. 12) WRITE(6,1605)                  
cc 1605 FORMAT(1H+,90X,'----- IF A PEAK APPEARS HERE')                    
cc      IF(IM1 .EQ. 43 .AND. IPPP .EQ. 12) WRITE(6,1606)                  
cc 1606 FORMAT(1H+,101X,'TRY TRADING-DAY ADJUSTMENT')                     
cc      F(IT)=BL                                                          
cc 2040 CONTINUE                                                          
cc 2060 FORMAT(1H ,I5,2X,D14.5,7X,61A1)                                   
cc      IF(IPUNC .EQ. 0) RETURN                                           
C                                                                       
cc      WRITE(7,700) H1                                                   
cc      WRITE(7,701)                                                      
cc      WRITE(7,702) (SXX(I),I=1,H1)                                      
cc  700 FORMAT(I5)                                                        
cc  701 FORMAT('(6D12.5)')                                                
cc  702 FORMAT(6D12.5)                                                    
      RETURN                                                            
cxx   60 FORMAT(//1H ,'RATIONAL SPECTRUM')                                 
cxx   61 FORMAT(1H ,2HL=,I5,2X,2HK=,I5,2X,2HH=,I5)                         
cxx   62 FORMAT(1H ,5X,'AR(I)')                                            
cxx   63 FORMAT(1H ,5X,'MA(I)')                                            
cxx  160 FORMAT(1H ,17HINITIAL CONDITION)                                  
cxx  164 FORMAT(1H ,6HSGME2=,D12.5)                                        
cxx  264 FORMAT(1H ,'AR-MA MODEL:')                                        
cxx  265 FORMAT(1H ,'Y(N)+AR(1)Y(N-1)+...+AR(L)Y(N-L)=X(N)+MA(1)X(N-1)+...+
cxx     AMA(K)X(N-K)')                                                     
      END                                                               
