      SUBROUTINE SIMCONF(D,K,H,L,R,B0,BX0,S0,Q0,BC,BD,G,AVY,SI,S2)
C
      INCLUDE 'timsac_f.h'
C
cc      PROGRAM SIMCON                                                    
C     PROGRAM 74.3.2.  OPTIMAL CONTROLLER DESIGN AND SIMULATION         
C-----------------------------------------------------------------------
C     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
C     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
C         TOKYO                                                         
C     ** DATE OF THE LATEST REVISION: MARCH 25, 1977                    
C     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN                       
C        "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(2
C        BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
C        NO.6 MARCH 1976, THE INSTITUTE OF STATISTICAL MATHEMATICS      
C-----------------------------------------------------------------------
C     THIS PROGRAM PRODUCES OPTIMAL CONTROLLER GAIN AND SIMULATES THE   
C     CONTROLLED PROCESS. THE BASIC STATE SPACE MODEL IS OBTAINED       
C     FROM THE AUTOREGRESSIVE MOVING AVERAGE MODEL OF A VECTOR          
C     PROCESS Y(I); Y(I)+B(1)Y(I-1)+..+B(K)Y(I-K)=X(I)+A(1)X(I-1)+..    
C     ..+A(K-1)X(I-K+1).                                                
C                                                                       
C     THE FOLLOWING INPUTS ARE REQUIRED:                                
C      THE OUTPUTS OF PROGRAM MARKOV:                                   
C       (D,K): INTEGERS                                                 
C             D, DIMENSION OF Y(I) (LESS THAN OR EQUAL TO MJD)          
C             K, ORDER OF THE PROCESS (LESS THAN OR EQUAL TO MJK)       
C       (B(I),I=1,K): MATRICES OF AUTOREGRESSIVE COEFFICIENTS           
C       (W(I),I=1,K-1): IMPULSE RESPONCE MATRICES OF THE VECTOR VARIABLE
C                       Y(I) TO THE INNOVATION INPUT                    
C      S: D*D COVARIANCE MATRIX OF INNOVATION X(I)                      
C      (H,L,R): INTEGERS                                                
C              H, SPAN OF CONTROL PERFORMANCE EVALUATION                
C              L, LENGTH OF EXPERIMENTAL OBSERVATION                    
C              R, DIMENSION OF CONTROL INPUT (LESS THAN OR EQUAL TO D)  
C                 THE LAST R COMPONENTS OF Y(I) REPRESENT THE CONTROL IN
C                 AND THE REST REPRESENT THE OUTPUT OF THE SYSTEM.      
C      QX: WEIGHTING MATRIX OF PERFORMANCE, POSITIVE DEFINITE           
C                                                                       
C     BASIC EQUATION OF THE SYSTEM IS AS FOLLOWS;                       
C     V(I+1)=A*V(I)+BD*SI(I)+BC*CI(I+1): STATE VECTOR                   
C     CI(I+1)=G*V(I)                   : CONTROLLER INPUT               
C     Y(I)=C*V(I)                      : OBSERVED OUTPUT                
C     WHERE        +-                              -+                   
C              A = I     0       I       0 ..     0 I ,                 
C                  I     0       0       I ..     0 I                   
C                  I     .       .       . ..     . I                   
C                  I -B(K) -B(K-1) -B(K-2) .. -B(1) I                   
C                  +-                              -+                   
C                                                                       
C                  +-                -+                                 
C              C = I   I 0 0 ... 0    I,                                
C                  +-                -+                                 
C                                                                       
C              +-     -+     +-       -+                                
C              I       I     I  I      I                                
C              I       I     I  W(1)   I                                
C              I BD BC I  =  I  W(2)   I  =  BX,                        
C              I       I     I   .     I                                
C              I       I     I  W(K-1) I                                
C              +-     -+     +-       -+                                
C                                                                       
C              G = INVERSE(-(BC'*Q(H)*BC))*BC'*Q(H)*A,                  
C              Z(I+1) = Q(I)-Q(I)*BC*INVERSE(BC'*Q(I)*BC)*BC'*Q(I) (I=0,
C              Q(I+1) = A'*Z(I+1)*A+Q(0) (I=0,H-1),                     
C                     +-       -+                                       
C              Q(0) = I QX 0..0 I                                       
C                     I  0 0..0 I                                       
C                     I  0 0..0 I                                       
C                     I  . ...0 I                                       
C                     I  0 0..0 I                                       
C                     +-       -+                                       
C              SI(I)= OUTPUT (FIRST D-R COMPONENTS OF Y(I)) INNOVATION. 
C     B(I)(I=1,K) ARE THE AUTOREGRESSIVE COEFFICIENTS OF THE AR-MA      
C     REPRESENTATION OF Y(I).                                           
C     W(J)(J=0,K-1)(W(0)=I) ARE THE IMPULSE RESPONSE MATRICES OF Y(I+J) 
C     THE INNOVATION X(I).                                              
C     BY MODIFYING THE INPUT LOADING PROGRAM, Q(0) MAY BE PUT EQUAL TO A
C     ARBITRARY KD*KD POSITIVE DEFINITE MATRIX.                         
C     THE DIMENSIONS OF THE MATRICES ARE AS FOLLOWES:                   
C     A:  (K*D)*(K*D) MATRIX                                            
C     B:  D*D         MATRIX                                            
C     BD: (K*D)*(D-R) MATRIX                                            
C     BC: (K*D)*R     MATRIX                                            
C     W:  D*D         MATRIX                                            
C     Q:  (K*D)*(K*D) MATRIX                                            
C     QX: D*D         MATRIX                                            
C     Z:  (K*D)*(K*D) MATRIX                                            
C     G:  R*(K*D)     MATRIX                                            
C     S:  D*D         MATRIX                                            
C     V:  (K*D)*1     VECTOR                                            
C     SI: (D-R)*1     VECTOR                                            
C     CI: R*1         VECTOR                                            
C     Y:  D*1         VECTOR                                            
C                                                                       
cc      !DEC$ ATTRIBUTES DLLEXPORT :: SIMCONF
C
cxx      IMPLICIT REAL*8 (A-H,O-Z)                                         
cc      REAL*4 RANDOM                                                     
cxx      REAL*8 RANDM 
cxx      INTEGER D,H,R,H1                                                  
cc      DIMENSION B(10,10,10),BD(100,10),BC(100,10)                       
cc      DIMENSION Q(100,100),Q0(100,100),Z(100,100),G(10,100)             
cc      DIMENSION W1(10)                                                  
cc      DIMENSION S(10,10),WIN(10,10),Y(10,1000),W2(100,100)              
cc      DIMENSION SI(10),CI(10),V(100),AV(100),BCN(100),BDN(100)          
cc      DIMENSION W3(10,10),W4(10,10),S2(10)                              
cc      DIMENSION R1(10,100),R3(10,100)                                   
cc      DIMENSION BX(100,10)                                              
cxx      DIMENSION B0(D,D,K),B(D,D,K),BD(K*D,D-R),BC(K*D,R)
cxx      DIMENSION Q(K*D,K*D),Q0(D,D),Z(K*D,K*D),G(R,K*D)
cxx      DIMENSION W1(D)
cxx      DIMENSION S0(D,D),S(D-R,D-R)
cxx      DIMENSION WIN1(D,D),WIN2(R,D),Y(D,L),W2(K*D,K*D)
cxx      DIMENSION SI(D),CI(R),V(K*D),AV(K*D),AVY(D),BCN(K*D),BDN(K*D)
cxx      DIMENSION W31(R,R),W32(D,D),W41(D,D),W42(R,D),S2(D)
cxx      DIMENSION R1(R,K*D),R3(R,K*D)
cxx      DIMENSION BX0((K-1)*D,D),BX(K*D,D)
      INTEGER :: D, K, H, L, R
      REAL(8) :: B0(D,D,K), BX0((K-1)*D,D), S0(D,D), Q0(D,D), BC(K*D,R),
     1           BD(K*D,D-R), G(R,K*D), AVY(D), SI(D), S2(D)
      INTEGER H1
      REAL(8) :: RANDM, B(D,D,K), Q(K*D,K*D),Z(K*D,K*D), W1(D),
     1           S(D-R,D-R), WIN1(D,D), WIN2(R,D), Y(D,L),
     2           W2(K*D,K*D), CI(R), V(K*D), AV(K*D), BCN(K*D),
     3           BDN(K*D), W31(R,R), W32(D,D), W41(D,D), W42(R,D),
     4           R1(R,K*D), R3(R,K*D), BX(K*D,D), CST0, CST1, CSTM6,
     5           BDMN, AL, SUM, WDET
cc      EQUIVALENCE(Q0(1,1),BCN(1))                                       
cc      EQUIVALENCE(WIN(1,1),W1(1),SI(1),CI(1))                           
cc      EQUIVALENCE(V(1),S2(1))                                           
cc      EQUIVALENCE(Q(1,1),Z(1,1),Y(1,1))                                 
cc      EQUIVALENCE (R1(1,1),G(1,1))                                      
cc      EQUIVALENCE(BX(1,1),R3(1,1),AV(1))                                
cc      DATA B/1000*0.0D-00/,BC/1000*0.0D-00/,BD/1000*0.0D-00/            
cc      DATA Q0/10000*0.0D-00/,Q/10000*0.0D-00/,Z/10000*0.0D-00/          
cc      DATA S/100*0.0D-00/,SI/10*0.0D-00/,CI/10*0.0D-00/                 
cc      DATA W1/10*0.0D-00/,W2/10000*0.0D-00/                             
cc      DATA R1/1000*0.0D-00/                                             
cc      DATA BX/1000*0.0D-00/,WIN/100*0.0D-00/,W3/100*0.0D-00/            
C     ABSOLUTE DIMENSIONS                                               
cc      MJD=10                                                            
cc      MJK=10                                                            
cc      MJKD=MJK*MJD                                                      
C
C     INPUT / OUTPUT DATA FILE OPEN
cc	CALL SETWND
cc	CALL FLOPN2(NFL)
cc	IF (NFL.EQ.0) GO TO 999
C                                                                       
C     INITIAL CONDITION INPUT                                           
cc      READ(5,800) D,K                                                   
cc      WRITE(6,900)                                                      
      KD=K*D                                                            
      KX=KD-D
C     AR-COEFFICIENT MATRICES B(I)(I=1,K) INPUT                         
      ISW=1                                                             
      DO  107 IX=1,K                                                    
C     COMMON SUBROUTINE CALL                                            
cc      CALL REMATX(WIN,D,D,ISW,MJD,MJD)                                  
      DO  106 I=1,D                                                     
      DO  105 J=1,D                                                     
cc  105 B(I,J,IX)=-WIN(I,J)                                               
cxx  105 B(I,J,IX)=-B0(I,J,IX)
      B(I,J,IX)=-B0(I,J,IX)
  105 CONTINUE
  106 CONTINUE                                                          
  107 CONTINUE                                                          
C                                                                       
      CST0=0.0D-00                                                      
      CST1=1.0D-00                                                      
      CSTM6=-6.0D-00                                                    
c
c     initialize data area
      CALL DINIT(Q,KD*KD,CST0)
      CALL DINIT(Z,KD*KD,CST0)
      CALL DINIT(SI,D,CST0)
      CALL DINIT(CI,R,CST0)
      CALL DINIT(W1,D,CST0)
      CALL DINIT(W2,KD*KD,CST0)
      CALL DINIT(R1,R*KD,CST0)
      CALL DINIT(WIN1,D*D,CST0)
      CALL DINIT(WIN2,R*D,CST0)
      CALL DINIT(W31,R*R,CST0)
      CALL DINIT(W32,D*D,CST0)
      CALL DINIT(BX,KD*D,CST0)
c
C     MATRIX BX SET                                                     
      K1=K-1                                                            
      K2=K1*D                                                           
      DO 102 I=1,D                                                      
cxx  102 BX(I,I)=CST1
      BX(I,I)=CST1
  102 CONTINUE
C                                                                       
C     IMPULSE RESPONSE MATRICES W(I),I=1,K-1 INPUT                      
cc      I1=0                                                              
cc      DO  117 IX=1,K1                                                   
C     COMMON SUBROUTINE CALL                                            
C     MATRIX INPUT, ROWWISE                                             
cc      CALL REMATX (WIN,D,D,ISW,MJD,MJD)                                 
cc      I1=I1+D                                                           
cc      DO  116 I=1,D                                                     
cc      I2=I1+I                                                           
cc      DO 115 J=1,D                                                      
cc  115 BX(I2,J)=WIN(I,J)                                                 
cc  116 CONTINUE                                                          
cc  117 CONTINUE                                                          
cxx      DO 116 I=1,(K-1)*D
      DO 126 I=1,(K-1)*D
      DO 116 J=1,D
         BX(I+D,J)=BX0(I,J)
  116 CONTINUE
  126 CONTINUE

C                                                                       
C     INNOVATION COVARIANCE MATRIX INPUT                                
C     COMMON SUBROUTINE CALL                                            
cc      CALL REMATX (S,D,D,ISW,MJD,MJD)                                   
cxx      DO 117 I=1,D-R
      DO 127 I=1,D-R
      DO 117 J=1,D-R
         S(I,J)=S0(I,J)
  117 CONTINUE
  127 CONTINUE
cc      READ (5,800) H,L,R                                                
cc      WRITE(6,901) D,K,H,L,R                                            
C     MATRIX Q(0) SET UP                                                
C     COMMON SUBROUTINE CALL                                            
C     FOR A GENERAL Q(0) REPLACE THE FOLLWING STATEMENT 108-109 BY      
C     CALL REMATX(Q0,KD,KD,ISW,MJKD,MJKD)                               
cc  108 CALL REMATX(W3,D,D,ISW,MJD,MJD)                                   
C                                                                       
C                                                                       
cc      DO  109  I=1,D                                                    
cc      DO  118  J=1,D                                                    
cc  118 Q0(I,J)=W3(I,J)                                                   
cc  109 CONTINUE                                                          
cc      DO 119 I=1,KD
cc      DO 119 J=1,KD
cxx      DO 119 I=1,D
      DO 120 I=1,D
      DO 119 J=1,D
cxx  119 Q(I,J)=Q0(I,J)
      Q(I,J)=Q0(I,J)
  119 CONTINUE
  120 CONTINUE
C                                                                       
C     MATRIX BC AND BD SET UP                                           
      KDR=D-R                                                           
      DO 1120 I=1,KD                                                    
      DO 1119 J=1,KDR
cxx 1119 BD(I,J)=BX(I,J)
      BD(I,J)=BX(I,J)
 1119 CONTINUE
 1120 CONTINUE                                                          
      DO 1122 I=1,KD                                                    
      DO 1121 J=1,R                                                     
      J1=KDR+J                                                          
cxx 1121 BC(I,J)=BX(I,J1) 
      BC(I,J)=BX(I,J1)
 1121 CONTINUE
 1122 CONTINUE                                                          
C                                                                       
C     INITIAL MATRIX PRINT OUT FOR DEBUGGING                            
cc      DO  707 IX=1,K                                                    
cc      DO  708 I=1,D                                                     
cc      DO  709 J=1,D                                                     
cc  709 WIN(I,J)=-B(I,J,IX)                                               
cc  708 CONTINUE                                                          
cc      WRITE (6,700) IX                                                  
cc      CALL SUBMPR (WIN,D,D,MJD,MJD)                                     
cc  707 CONTINUE                                                          
cc      WRITE(6,702)                                                      
cc      CALL SUBMPR (BC,KD,R,MJKD,MJD)                                    
cc      WRITE(6,703)                                                      
cc      CALL SUBMPR (BD,KD,KDR,MJKD,MJD)                                  
cc      WRITE(6,704)                                                      
cc      CALL SUBMPR(W3,D,D,MJD,MJD)                                       
cc      WRITE(6,706)                                                      
cc      CALL SUBMPR(S,D,D,MJD,MJD)                                        
C     INITIAL PRINT FOR DEBUGGING END                                   
C     ITERATIVE COMPUTATION OF G                                        
      H1=H+1                                                            
      DO 200 IX=1,H1                                                    
C     R1=BC'*Q R1 STORE                                                 
cc      CALL MULTRX(BC,KD,R,MJKD,MJD,Q,KD,KD,MJKD,MJKD,R1,R,KD,MJD,MJKD,2)
      CALL MULTRX(BC,KD,R,Q,KD,KD,R1,R,KD,2)
C     W3=R1*BC                                                          
cc      CALL MULTRX(R1,R,KD,MJD,MJKD,BC,KD,R,MJKD,MJD,W3,R,R,MJD,MJD,1)   
      CALL MULTRX(R1,R,KD,BC,KD,R,W31,R,R,1)
C     W3=(INVERSE OF W3)                                                
cc      CALL INVDET(W3,WDET,R,MJD)                                        
      CALL INVDETS(W31,WDET,R)
C     R3=W3*R1 R3 STORE                                                 
cc      CALL MULTRX(W3,R,R,MJD,MJD,R1,R,KD,MJD,MJKD,R3,R,KD,MJD,MJKD,1)   
      CALL MULTRX(W31,R,R,R1,R,KD,R3,R,KD,1)
      IF(IX.EQ.H1) GO TO 210                                            
C     W2=R1'*R3                                                         
cc      CALL MULTRX(R1,R,KD,MJD,MJKD,R3,R,KD,MJD,MJKD,W2,KD,KD,MJKD,MJKD,2
cc     A)                                                                 
      CALL MULTRX(R1,R,KD,R3,R,KD,W2,KD,KD,2)
C     Z=Q-W2                                                            
cc      CALL SUBTAC(Q,W2,Z,KD,KD,MJKD,MJKD)                               
      CALL SUBTAC(Q,W2,Z,KD,KD)
C                                                                       
C     W2=A'*Z COMPUTATION:                                              
C                                                                       
C     +-                               -+     +-                        
C     I     0  0  0  0  ..  0 -B(K)'    I     I     Z(1,1) .. Z(1,K)    
C     I     I  0  0  0  ..  0 -B(K-1)'  I     I      .          .       
C     I     0  I  0  0  ..  0 -B(K-2)'  I     I      .          .       
C     I     0  0  I  0  ..  0 -B(K-3)'  I  *  I      .          .       
C     I     0  0  0  I  ..  0 -B(K-4)'  I     I      .          .       
C     I     .  .  .  .  ..  .    .      I     I      .          .       
C     I     0  0  0  0  ..  I -B(1)'    I     I     Z(K,1) .. Z(K,K)    
C     +-                               -+     +-                        
C                                                                       
C         +-                                                    -+      
C         I   -B(K)'Z(K,1)          .... -B(K)'Z(K,K)            I      
C         I   -B(K-1)'Z(K,1)+Z(1,1) .... -B(K-1)'Z(K,K)+Z(1,K)   I      
C      =  I   -B(K-2)'Z(K,1)+Z(2,1) .... -B(K-2)'Z(K,K)+Z(2,K)   I ,    
C         I   -B(K-3)'Z(K,1)+Z(3,1) .... -B(K-3)'Z(K,K)+Z(3,K)   I      
C         I             .           ....           .             I      
C         I             .           ....           .             I      
C         I             .           ....           .             I      
C         I   -B(1)'Z(K,1)+Z(K-1,1) .... -B(1)'Z(K,K)+Z(K-1,K)   I      
C         +-                                                    -+      
C                                                                       
C     WHERE Z(J,L)(J=1,K;L=1,K) ARE THE D*D SUB MATRICES OF THE         
C     BLOCK MATRIX Z.                                                   
C     A'*Z IS STORED IN W2.                                             
      IZ=K+1                                                            
cc      KX=KD-D                                                           
      J1=D+1                                                            
      KIZ=0                                                             
      DO  408 IY=1,K                                                    
      IZ=IZ-1                                                           
      DO  402 I=1,D                                                     
      DO  401 J=1,D                                                     
cc  401 WIN(I,J)=B(I,J,IZ)                                                
cxx  401 WIN1(I,J)=B(I,J,IZ)
      WIN1(I,J)=B(I,J,IZ)
  401 CONTINUE
  402 CONTINUE                                                          
      KJZ=0                                                             
      DO  407 JY=1,K                                                    
      DO  404 I=1,D                                                     
      KIX=KX+I                                                          
      DO  403 J=1,D                                                     
      KJX=KJZ+J                                                         
cc  403 W3(I,J)=Z(KIX,KJX)                                                
cxx  403 W32(I,J)=Z(KIX,KJX)
      W32(I,J)=Z(KIX,KJX)
  403 CONTINUE
  404 CONTINUE                                                          
cc      CALL MULTRX (WIN,D,D,MJD,MJD,W3,D,D,MJD,MJD,                      
cc     AW4,D,D,MJD,MJD,2)                                                 
      CALL MULTRX (WIN1,D,D,W32,D,D,W41,D,D,2)
      DO  406 I=1,D                                                     
      KIX=KIZ+I                                                         
      DO  405 J=1,D                                                     
      KJX=KJZ+J                                                         
cc  405 W2(KIX,KJX)=W4(I,J)                                               
cxx  405 W2(KIX,KJX)=W41(I,J)
      W2(KIX,KJX)=W41(I,J)
  405 CONTINUE
  406 CONTINUE                                                          
      KJZ=KJZ+D                                                         
  407 CONTINUE                                                          
      KIZ=KIZ+D                                                         
  408 CONTINUE                                                          
      DO  411 IY=J1,KD                                                  
      DO  410 JY=1,KD                                                   
      IYY=IY-D                                                          
cxx  409 W2(IY,JY) = W2(IY,JY) + Z(IYY,JY)                                 
      W2(IY,JY) = W2(IY,JY) + Z(IYY,JY)
  410 CONTINUE                                                          
  411 CONTINUE                                                          
C                                                                       
C     Q=W2*A COMPUTATIONF W2=A'*Z                                       
C                                                                       
C     +-                   -+     +-                                -+  
C     I  W2(1,1)...W2(1,K)  I     I     0     I        0   ..    0   I  
C     I     .          .    I  *  I     0     0        I   ..    0   I  
C     I     .          .    I     I     .     .        .   ..    .   I  
C     I  W2(K,1)...W2(K,K)  I     I  -B(K) -B(K-1) -B(K-2) .. -B(1)  I  
C     +-                   -+     +-                                -+  
C                                                                       
C         +-                                                            
C         I  -W2(1,K)B(K) -W2(1,K)B(K-1)+W2(1,1) ... -W2(1,K)B(1)+W2(1,K
C         I  -W2(2,K)B(K) -W2(2,K)B(K-1)+W2(2,1) ... -W2(2,K)B(1)+W2(2,K
C         I  -W2(3,K)B(K) -W2(3,K)B(K-1)+W2(3,1) ... -W2(3,K)B(1)+W2(3,K
C      =  I  -W2(4,K)B(K) -W2(4,K)B(K-1)+W2(4,1) ... -W2(4,K)B(1)+W2(4,K
C         I         .                  .         ...        .           
C         I         .                  .         ...        .           
C         I         .                  .         ...        .           
C         I  -W2(K,K)B(K) -W2(K,K)B(K-1)+W2(K,1) ... -W2(K,K)B(1)+W2(K,K
C         +-                                                            
C                                                                       
C     WHERE W2(J,L)(J=1,K;L=1,K) ARE THE (J,L) D*D SUB MATRICES         
C     OF THE BLOCK MATRIX W2.                                           
C     W2*A IS STORED IN Q.                                              
      IZ=K+1                                                            
      KJZ=0                                                             
      DO  419 JY=1,K                                                    
      IZ=IZ-1                                                           
      DO  413 I=1,D                                                     
      DO  412 J=1,D                                                     
cxx  412 W32(I,J)=B(I,J,IZ)
      W32(I,J)=B(I,J,IZ)
  412 CONTINUE
  413 CONTINUE                                                          
      KIZ=0                                                             
      DO  418 IY=1,K                                                    
      DO  415 J=1,D                                                     
      KJX=KX+J                                                          
      DO  414 I=1,D                                                     
      KIX=KIZ+I                                                         
cc  414 WIN(I,J)=W2(KIX,KJX)                                              
cxx  414 WIN1(I,J)=W2(KIX,KJX)
      WIN1(I,J)=W2(KIX,KJX)
  414 CONTINUE
  415 CONTINUE                                                          
cc      CALL MULTRX (WIN,D,D,MJD,MJD,W3,D,D,MJD,MJD,                      
cc     AW4,D,D,MJD,MJD,1)                                                 
      CALL MULTRX (WIN1,D,D,W32,D,D,W41,D,D,1)
      DO  417 J=1,D                                                     
      KJX=KJZ+J                                                         
      DO  416 I=1,D                                                     
      KIX=KIZ+I                                                         
cc  416 Q(KIX,KJX)=W4(I,J)                                                
cxx  416 Q(KIX,KJX)=W41(I,J)                                                
      Q(KIX,KJX)=W41(I,J)
  416 CONTINUE
  417 CONTINUE                                                          
      KIZ=KIZ+D                                                         
  418 CONTINUE                                                          
      KJZ=KJZ+D                                                         
  419 CONTINUE                                                          
      DO  422 JY=J1,KD                                                  
      DO  421 IY=1,KD                                                   
      JYY=JY-D
cxx  420 Q(IY,JY)=Q(IY,JY)+W2(IY,JYY)                                      
      Q(IY,JY)=Q(IY,JY)+W2(IY,JYY)
  421 CONTINUE                                                          
  422 CONTINUE                                                          
cc      DO 142 I=1,KD                                                     
cc      DO 142 J=1,KD                                                     
cxx      DO 142 I=1,D
      DO 143 I=1,D
      DO 142 J=1,D
cxx  142 Q(I,J)=Q(I,J)+Q0(I,J)
      Q(I,J)=Q(I,J)+Q0(I,J)
  142 CONTINUE
  143 CONTINUE
  200 CONTINUE                                                          
C                                                                       
C     G=-(INVERSE(BC'*Q*BC))*BC'*Q*A                                    
C     G=R3*A; R3=INVERSE(BC'*Q*BC)*BC'*Q:                               
C                                                                       
C     +-                   -+     +-                               -+   
C     I  W(1) W(2) .. W(K)  I  *  I    0      I       0   ..    0   I   
C     +-                   -+     I    0      0       I   ..    0   I   
C                                 I    0      0       0   ..    0   I   
C                                 I    .      .       .   ..    .   I   
C                                 I -B(K) -B(K-1) -B(K-2) .. -B(1)  I   
C                                 +-                               -+   
C                                                                       
C       +-                                                              
C       I                                                               
C     = I -W(K)B(K) -W(K)B(K-1)+W(1) -W(K)B(K-2)+W(2) ... -W(K)B(1)+W(K-
C       I                                                               
C       +-                                                              
C                                                                       
C     WHERE W(J)(J=1,K) ARE THE J-TH R*D SUB MATRIX OF THE BLOCK MATRIX 
C     R3*A IS STORED IN G.                                              
  210 CONTINUE                                                          
      DO  455 I=1,R                                                     
      DO  454 J=1,D                                                     
      KJX=KX+J                                                          
cc  454 WIN(I,J)=R3(I,KJX)                                                
cxx  454 WIN2(I,J)=R3(I,KJX)
      WIN2(I,J)=R3(I,KJX)
  454 CONTINUE
  455 CONTINUE                                                          
      IZ=K+1                                                            
      DO  461 JY=1,K                                                    
      IZ=IZ-1                                                           
      DO  453 I=1,D                                                     
      DO  452 J=1,D                                                     
cc  452 W3(I,J)=B(I,J,IZ)                                                 
cxx  452 W32(I,J)=B(I,J,IZ)
      W32(I,J)=B(I,J,IZ)
  452 CONTINUE
  453 CONTINUE                                                          
cc      CALL MULTRX (WIN,R,D,MJD,MJD,W3,D,D,MJD,MJD,                      
cc     AW4,R,D,MJD,MJD,1)                                                 
      CALL MULTRX (WIN2,R,D,W32,D,D,W42,R,D,1)
      KJZ=(JY-1)*D                                                      
      DO 457 I=1,R                                                      
      DO 456 J=1,D                                                      
      KJY=KJZ+J                                                         
cc  456 G(I,KJY)=W4(I,J)                                                  
cxx  456 G(I,KJY)=W42(I,J)
      G(I,KJY)=W42(I,J)
  456 CONTINUE
  457 CONTINUE                                                          
cc      IF  (KJZ) 461,461,458                                             
      IF  (KJZ .LE. 0) GO TO 461
      IF  (KJZ .GT. 0) GO TO 458                                             
  458 DO  460 I=1,R                                                     
      DO  459 J=1,D                                                     
      KJX=KJZ+J-D                                                       
      KJY=KJZ+J                                                         
cxx  459 G(I,KJY)=G(I,KJY)+R3(I,KJX)
      G(I,KJY)=G(I,KJY)+R3(I,KJX)
  459 CONTINUE
  460 CONTINUE                                                          
  461 CONTINUE                                                          
      DO 202 I=1,R                                                      
      DO 201 J=1,KD                                                     
cxx  201 G(I,J)=-G(I,J)
      G(I,J)=-G(I,J)
  201 CONTINUE
  202 CONTINUE                                                          
C     G PRINT OUT                                                       
cc      WRITE(6,903) R,KD                                                 
C     COMMON SUBROUTINE CALL                                            
cc      CALL SUBMPR(G,R,KD,MJD,MJKD)                                      
C                                                                       
C                                                                       
C                                                                       
C     ****************                                                  
C     SIMULATION START                                                  
C     ****************                                                  
C     V(I),I=1,KD  CLEAR                                                
cxx      DO  301 I=1,KD                                                    
cxx  301 V(I)=CST0
      V(1:KD)=CST0
C     CONSTANT FOR NOISE GENERATION                                     
C     COMMON SUBROUTINE CALL                                            
cc      SI(1)=RANDOM(1)                                                   
      SI(1)=RANDM(1,KK1,KK2,KK3,KK4)
C     OBSERVED SYSTEM OUTPUT GENERATION I=1,L                           
C     COVARIANCE MATRIX FACTORIZATION                                   
C     COMMON SUBROUTINE CALL                                            
cc      CALL LTINV  (S,KDR,MJD)                                           
      CALL LTINV  (S,KDR)
C     MATRIX S ARRANGEMENT 
      IF(KDR.EQ.1) GO TO 260                                            
cxx      DO 12 I=2,KDR
      DO 13 I=2,KDR
      IM1=I-1                                                           
      DO 12 J=1,IM1                                                     
cxx   12 S(I,J)=S(J,I)
      S(I,J)=S(J,I)
   12 CONTINUE
   13 CONTINUE
  260 CONTINUE                                                          
C                                                                       
      DO  500 J=1,L                                                     
C     AV=A*V COMPUTATION:                                               
C                                                                       
C     +-                               -+     +-    -+                  
C     I     0     I       0   ...    0  I     I V(1) I                  
C     I     0     0       I   ...    0  I     I V(2) I                  
C     I     0     0       0   ...    0  I  *  I V(3) I                  
C     I     .     .       .   ...    .  I     I   .  I                  
C     I     .     .       .   ...    .  I     I   .  I                  
C     I -B(K) -B(K-1) -B(K-2) ... -B(1) I     I V(K) I                  
C     +-                               -+     +-    -+                  
C                                                                       
C       +-                                   -+                         
C       I  V(2)                               I                         
C     = I  V(3)                               I ,                       
C       I  V(4)                               I                         
C       I    .                                I                         
C       I    .                                I                         
C       I  -B(K)V(1)-B(K-1)V(2)-...-B(1)V(K)  I                         
C       +-                                   -+                         
C                                                                       
C     WHERE V(J)(J=1,K) ARE THE D-DIMENSIONAL SUB VECTORS OF V.         
C     A*V IS STORED IN THE VECTOR AV.                                   
cxx      DO  432 I1=1,D                                                    
cxx  432 W1(I1)=CST0
      W1(1:D)=CST0 
      IZ=K+1                                                            
      KJZ=0                                                             
      DO  436 JY=1,K                                                    
      IZ=IZ-1                                                           
      DO  434 I1=1,D                                                    
      DO  433 I2=1,D                                                    
      J2=KJZ+I2                                                         
cxx  433 W1(I1)=W1(I1)+B(I1,I2,IZ)*V(J2)
      W1(I1)=W1(I1)+B(I1,I2,IZ)*V(J2)
  433 CONTINUE
  434 CONTINUE                                                          
      KJZ=KJZ+D                                                         
  436 CONTINUE                                                          
      II1=0                                                             
      DO  439 IY=1,K1                                                   
      I1=II1+1                                                          
      I2=II1+D                                                          
      DO  438 IZ=I1,I2                                                  
      J1=IZ+D                                                           
      AV(IZ)=V(J1)                                                      
  438 CONTINUE                                                          
      II1=II1+D                                                         
  439 CONTINUE                                                          
      DO  441 I1=1,D                                                    
      I2=KX+I1                                                          
      AV(I2)=W1(I1)                                                     
  441 CONTINUE                                                          
C     CI=G*V; CONTROLLER INPUT                                          
C     COMMON SUBROUTINE CALL                                            
cc      CALL MULVER(G,V,CI,R,KD,MJD,MJKD)                                 
      CALL MULVER(G,V,CI,R,KD)
C     BCN=BC*CI                                                         
C     COMMON SUBROUTINE CALL                                            
cc      CALL MULVER(BC,CI,BCN,KD,R,MJKD,MJD)                              
      CALL MULVER(BC,CI,BCN,KD,R)
C     RANDOM VECTOR SI GENERATION                                       
      DO 302 I=1,KDR                                                    
      BDMN=CSTM6                                                        
      DO 303 JRN=1,12                                                   
cc  303 BDMN=BDMN+RANDOM(0)
cxx  303 BDMN=BDMN+RANDM(0,KK1,KK2,KK3,KK4)
      BDMN=BDMN+RANDM(0,KK1,KK2,KK3,KK4)
  303 CONTINUE
cxx  302 BDN(I)=BDMN
      BDN(I)=BDMN
  302 CONTINUE
C     RANDOM NORMAL VECTOR BDN GENERATION                               
C     COMMON SUBROUTINE CALL                                            
cc      CALL LTRVEC(S,BDN,SI,KDR,KDR,MJD,MJD)
      CALL LTRVEC(S,BDN,SI,KDR,KDR)
C     COMMON SUBROUTINE CALL                                            
cc      CALL MULVER(BD,SI,BDN,KD,KDR,MJKD,MJD)
      CALL MULVER(BD,SI,BDN,KD,KDR)
C     V=AV+BCN+BDN                                                      
      DO  309 I=1,KD                                                    
cxx  309 V(I) = AV(I)+BCN(I)+BDN(I)
      V(I) = AV(I)+BCN(I)+BDN(I)
  309 CONTINUE
      DO  310 I=1,D                                                     
cxx  310 Y(I,J)=V(I)
      Y(I,J)=V(I)
  310 CONTINUE
  500 CONTINUE                                                          
C                                                                       
C     ************************************************                  
C     MEAN, VARIANCE AND STANDARD DEVIATION COMPUTATION                 
C     ************************************************                  
      AL=L                                                              
      AL=CST1/AL                                                        
      DO  600 I=1,D                                                     
      SUM=CST0                                                          
      DO  501 J=1,L                                                     
cxx  501 SUM=SUM+Y(I,J)                                                    
      SUM=SUM+Y(I,J)
  501 CONTINUE
cc      AV(I)=AL*SUM                                                      
      AVY(I)=AL*SUM
      SUM=CST0                                                          
      DO  502 J=1,L                                                     
cc  502 SUM=SUM+(Y(I,J)-AV(I))**2                                         
cxx  502 SUM=SUM+(Y(I,J)-AVY(I))**2
      SUM=SUM+(Y(I,J)-AVY(I))**2
  502 CONTINUE
      SI(I)=AL*SUM
      S2(I)= DSQRT(SI(I))                                               
  600 CONTINUE                                                          
cc      WRITE(6,902)                                                      
cc      WRITE(6,904)                                                      
cc      DO  610 I=1,D                                                     
cc  610 WRITE(6,905) I,AV(I),SI(I),S2(I)                                  
C                                                                       
cc      CALL FLCLS2(NFL)
cc  999 CONTINUE
      RETURN
cxx  800 FORMAT(6I5)
cxx  900 FORMAT('1  PROGRAM 74.3.2. OPTIMAL CONTROLLER DESIGN')            
cxx  901 FORMAT(1H ,'D (DIMENSION)=',I2,', K (ORDER)=',I3,', H (HORIZON)=',
cxx     A       I3,', L (LENGTH OF EXPERIMENTAL PERIOD)=',I5,              
cxx     A       ', R (DIMENSION OF CONTROL INPUT)=',I5)                    
cxx  902 FORMAT(//1H ,'OPTIMAL CONTROL SIMULATION')                        
cxx  903 FORMAT(////1H ,'COTROLLER GAIN G(',I3,',',I3,')')                 
cxx  904 FORMAT(1H ,' AVERAGE VALUE OF I-TH COMPONENT OF Y',9X,            
cxx     A'VARIANCE',20X,'STANDARD DEVIATION')                              
cxx  905 FORMAT(1H ,'      I = ',I4,D20.10,8X,D20.10,8X,D20.10)            
cxx  700 FORMAT(' MATRIX  B(',I3,')')                                   
cxx  701 FORMAT(' MATRIX  BX')                                          
cxx  702 FORMAT(' MATRIX  BC')                                          
cxx  703 FORMAT(' MATRIX  BD')                                          
cxx  704 FORMAT(' MATRIX  QX')                                          
cxx  706 FORMAT(' S: COVARIANCE MATRIX OF INNOVATION')                  
      END                                                               
C                                                                       
cc      SUBROUTINE MULTRX(X,MX,NX,MJ1,MJ2,Y,MY,NY,MJ3,MJ4,                
cc     A                  Z,MZ,NZ,MJ5,MJ6,IS)                             
      SUBROUTINE MULTRX(X,MX,NX,Y,MY,NY,Z,MZ,NZ,IS)
C     MATRIX MULTIPLICATION                                             
C     X: MX*NX, ABSOLUTE DIMENSION MJ1*MJ2                              
C     Y: MY*NY, ABSOLUTE DIMENSION MJ3*MJ4                              
C     Z: MZ*NZ, ABSOLUTE DIMENSION MJ5*MJ6                              
C     IS=1: Z=X*Y                                                       
C     IS=2: Z=X'*Y                                                      
C     IS=3: Z=X*Y'                                                      
cxx      IMPLICIT REAL*8 (A-H,O-Z)                                         
cc      DIMENSION X(MJ1,MJ2),Y(MJ3,MJ4),Z(MJ5,MJ6)                        
cxx      DIMENSION X(MX,NX),Y(MY,NY),Z(MZ,NZ)
      INTEGER :: MX, NX, MY, NY, MZ, NZ, IS
      REAL(8) :: X(MX,NX), Y(MY,NY), Z(MZ,NZ)
      REAL(8) :: CST0
      CST0=0.0D-00                                                      
      IF  (IS.EQ.2) GO TO 3050                                          
      IF  (IS.EQ.3) GO TO 3100                                          
cc      MZ=MX                                                             
cc      NZ=NY                                                             
      DO  3009 I=1,MX                                                   
      DO  3008 J=1,NY                                                   
      Z(I,J)=CST0                                                       
      DO  3007 K=1,NX                                                   
cxx 3007 Z(I,J)=Z(I,J)+X(I,K)*Y(K,J)
      Z(I,J)=Z(I,J)+X(I,K)*Y(K,J)
 3007 CONTINUE
 3008 CONTINUE                                                          
 3009 CONTINUE                                                          
      GO TO 3200                                                        
cc 3050 MZ=NX                                                             
cc      NZ=NY                                                             
 3050 CONTINUE
      DO  3059 I=1,NX                                                   
      DO  3058 J=1,NY                                                   
      Z(I,J)=CST0                                                       
      DO  3057 K=1,MX                                                   
cxx 3057 Z(I,J)=Z(I,J)+X(K,I)*Y(K,J)
      Z(I,J)=Z(I,J)+X(K,I)*Y(K,J)
 3057 CONTINUE
 3058 CONTINUE                                                          
 3059 CONTINUE                                                          
      GO TO 3200                                                        
cc 3100 MZ=MX                                                             
cc      NZ=MY                                                             
 3100 CONTINUE
      DO  3159 I=1,MX                                                   
      DO  3158 J=1,MY                                                   
      Z(I,J)=CST0                                                       
      DO  3157 K=1,NX                                                   
cxx 3157 Z(I,J)=Z(I,J)+X(I,K)*Y(J,K)
      Z(I,J)=Z(I,J)+X(I,K)*Y(J,K)
 3157 CONTINUE
 3158 CONTINUE                                                          
 3159 CONTINUE                                                          
 3200 RETURN                                                            
      END                                                               
C                                                                       
cc      SUBROUTINE INVDET(X,XDET,MM,MJ)                                   
      SUBROUTINE INVDETS(X,XDET,MM)
C     COMMON SUBROUTINE                                                 
C     THIS SUBROUTINE COMPUTES THE INVERSE AND DETERMINANT OF           
C     UPPER LEFT MM X MM OF X.                                          
C     X: ORIGINAL MATRIX                                                
C     MM: DIMENSION OF UPPER LEFT OF X (SHOULD BE LESS THAN 11)         
C     XDET: DETERMINANT OF UPPER LEFT MM X MM OF X                      
C     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE                   
C     THE INVERSE MATRIX IS OVERWRITTEN ON THE ORIGINAL.                
C     NEXT STATEMENT SHOULD BE REPLACED BY                              
C     IMPLICIT COMPLEX*16(X)                                            
C     FOR COMPLEX VERSION.  ALSO STATEMENT NO.1 NEEDS MODIFICATION.     
cxx      IMPLICIT REAL*8(X)                                                
cc      DIMENSION X(MJ,MJ)                                                
cc      DIMENSION IDS(10)                                                 
cxx      DIMENSION X(MM,MM)
cxx      DIMENSION IDS(MM)
      INTEGER :: MM
      REAL(8) :: X(MM,MM), XDET
      INTEGER :: IDS(MM)
      REAL(8) :: CST0, CST1, XMAXP, XC
      CST0=0.0D-00                                                      
      CST1=1.0D-00                                                      
      XDET=CST1                                                         
      DO 10 L=1,MM                                                      
C     PIVOTING AT L-TH STAGE                                            
      XMAXP=0.10000D-10                                                 
      MAXI=0                                                            
      DO 110 I=L,MM                                                     
C     FOR COMPLEX VERSION NEXT STATEMENT SHOULD BE REPLACED BY          
C     IF(CDABS(XMAXP).GE.CDABS(X(I,L))) GO TO 110                       
cxx    1 IF(DABS(XMAXP).GE.DABS(X(I,L))) GO TO 110                         
      IF(DABS(XMAXP).GE.DABS(X(I,L))) GO TO 110
      XMAXP=X(I,L)                                                      
      MAXI=I                                                            
  110 CONTINUE                                                          
      IDS(L)=MAXI                                                       
      IF(MAXI.EQ.L) GO TO 120                                           
      IF(MAXI.GT.0) GO TO 121                                           
      XDET=CST0                                                         
      GO TO 140                                                         
C     ROW INTERCHANGE                                                   
  121 DO 14 J=1,MM                                                      
      XC=X(MAXI,J)                                                      
      X(MAXI,J)=X(L,J)                                                  
cxx   14 X(L,J)=XC                                                         
      X(L,J)=XC
   14 CONTINUE
      XDET=-XDET                                                        
  120 XDET=CST1                                                         
      XC=CST1/XMAXP                                                     
      X(L,L)=CST1                                                       
      DO 11 J=1,MM                                                      
cxx   11 X(L,J)=X(L,J)*XC
      X(L,J)=X(L,J)*XC
   11 CONTINUE
      DO 12 I=1,MM                                                      
      IF(I.EQ.L) GO TO 12                                               
      XC=X(I,L)                                                         
      X(I,L)=CST0                                                       
      DO 13 J=1,MM                                                      
cxx   13 X(I,J)=X(I,J)-XC*X(L,J)
      X(I,J)=X(I,J)-XC*X(L,J)
   13 CONTINUE
   12 CONTINUE                                                          
   10 CONTINUE                                                          
      IF(MM.GT.1) GO TO 123                                             
      GO TO 140                                                         
C     COLUMN INTERCHANGE                                                
  123 MM1=MM-1                                                          
      DO 130 J=1,MM1                                                    
      MMJ=MM-J                                                          
      JJ=IDS(MMJ)                                                       
      IF(JJ.EQ.MMJ) GO TO 130                                           
      DO 131 I=1,MM                                                     
      XC=X(I,JJ)                                                        
      X(I,JJ)=X(I,MMJ)                                                  
cxx  131 X(I,MMJ)=XC
      X(I,MMJ)=XC
  131 CONTINUE
  130 CONTINUE                                                          
  140 RETURN                                                            
      END                                                               
