      SUBROUTINE  DECOMPF(DATA,N,IPAR,TREND,SEASNL,AR,
cx     *                    TRAD,NOISE,para,imiss,omaxx )
     *                    TRAD,NOISE,para,imiss,omaxx,ier )
C
      INCLUDE 'timsac_f.h'
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: DECOMPF
C
      PARAMETER (IOPT=1)
      PARAMETER (NIP=9, NPA=26)
cxx      REAL*8    DATA(N), para(NPA), omaxx
cxx      INTEGER   IPAR(NIP)
cxx      REAL*8    TREND(N),SEASNL(N),AR(N),TRAD(N),NOISE(N)
cxx      INTEGER   PERIOD, SORDER
      INTEGER :: N,  IPAR(NIP), imiss, ier
      REAL(8) :: DATA(N), TREND(N), SEASNL(N), AR(N), TRAD(N), NOISE(N),
     1           para(NPA), omaxx
      INTEGER :: PERIOD, SORDER
      COMMON    /COMSM2/  M1, M2, M3, M4, M5, M, L, PERIOD, SORDER,    
     *                    NYEAR, nmonth
cxx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH
      COMMON    /CCC/     ISW, ISMT, IDIF, LOG, MESH
c
      para(1:NPA) = 0.0D0
cxx      CALL  SPARAM0( N,IPAR,NIP,para,NPA )
      CALL  SPARAM0( IPAR,NIP )
      LM1 = L+M+1
c
cxx      call decompff (DATA,N,IPAR,TREND,SEASNL,AR,TRAD,NOISE,
      call decompff (DATA,N,TREND,SEASNL,AR,TRAD,NOISE,
cx     *               para,iopt,imiss,omaxx,LM1)
     *               para,iopt,imiss,omaxx,LM1,ier)
c
      return
      end
c
cxx      SUBROUTINE  DECOMPFF(DATA,N,IPAR,TREND,SEASNL,AR,TRAD,NOISE,
      SUBROUTINE  DECOMPFF(DATA,N,TREND,SEASNL,AR,TRAD,NOISE,
cx     *                     para,iopt,imiss,omaxx,LM1 )
     *                     para,iopt,imiss,omaxx,LM1,ier )
C
C                                                                       
c  Bug fixed (97/10/17)
c
c  99/8/12
c  iopt < 0 -> use para. not optimized. 
c  input 
c  para(1):tau21, para(2):tau22, para(3):tau23,
c  para(4,5,6,...): ARCOEF
c
c  m4 = 0(no trade), 1(2-TDF model), 6(7-TDF model) 
c
c  99/9/2  if iopt < 0 , need ispan (=ipar(9))
c 

C  ...  TIME SERIES DECOMPOSITION (SEASONAL ADJUSTMENT) ...             
C                                                                       
C     THE BASIC MODEL:                                                  
C                                                                       
C          y(n) = T(n) + AR(n) + S(n) + TD(n) + R(n) + W(n)             
C                                                                       
C       where                                                           
C            T(n):       trend component                                
C            AR(n):      AR process                                     
C            S(n):       seasonal component                             
C            TD(n):      trading day factor                             
C            R(n):       any other explanetory variables                
C            W(n):       observational noise                            
C                                                                       
C     COMPONENT MODELS:                                                 
C                                                                       
C       Trend component                                                 
C            T(N) =  T(N-1) + V1(N)                             :M1 = 1 
C            T(N) = 2T(N-1) - T(N-2) + V1(N)                    :M1 = 2 
C            T(N) = 3T(N-1) -3T(N-2) + T(N-2) + V1(N)           :M1 = 3 
C                                                                       
C       AR componet:                                                    
C            AR(n) = a(1)AR(n-1) + ... + a(m2)AR(n-m2) + V2(n)          
C                                                                       
C       Seasonal component:                                             
C            S(N) =  -S(N-1) - ... - S(N-PERIOD+1) + V3(N)    :SORDER=1 
C            S(N) = -2S(N-1) - ... -PERIOB*S(N-PERIOD+1)      :SORDER=2 
C                            - ... - S(n-2PERIOD+2) + V3(n)             
C       Trading day effect:                                             
C            TD(n) = b(1)*TRADE(n,1) + ... + b(7)*TRADE(n,7)            
C                                                                       
C            TRADE(n,i):  number of i-th days of the week in n-th data  
C            b(1) + ... + b(7) = 0                                      
C                                                                       
C     REFERENCES                                                        
C                                                                       
C       G. KITAGAWA (1981), A Nonstationary Time Series Model and Its   
C            Fitting by a Recursive Filter, Journal of Time Series      
C            Analysis, Vol.2, 103-116.                                  
C                                                                       
C       W. GERSCH and G. KITAGAWA (1983), The prediction of time series 
C            with Trends and Seasonalities, Journal of Business and     
C            Economic Statistics, Vol.1, 253-264.                       
C                                                                       
C       G. KITAGAWA (1984), A smoothness priors-state space modeling of 
C            Time Series with Trend and Seasonality, Journal of American
C            Statistical Association, VOL.79, NO.386, 378-389.          
C                                                                       
C     STRUCTURE OF THE PROGRAM                                          
C                                                                       
C       <DECOMP>                                                        
C          |---<SPARAM>                                                 
C          |      |---<ID>                                              
C          |      +---<PARCOR>                                          
C          |---<REDATA>                                                 
C          |---<AREA>                                                   
C          |---<LOGTRF>                                                 
C          |---<TRADE>                                                  
C          |---<EPARAM>                                                 
C          |      |---<SETFGH>                                          
C          |      |---<OPTMIZ>                                          
C          |      |      +---<LINEAR>                                   
C          |      |             +---<FUNCND>                            
C          |      |                    +---<FUNCSA>                     
C          |      |                           |---<ARCOEF>              
C          |      |                           |---<SMOTH3>              
C          |      |                           |      |---<HUSHL7>       
C          |      |                           |      |---<HUSHL4>       
C          |      |                           |      |---<RECOEF>       
C          |      |                           |      +---<ID>           
C          |      |                           +---<STATE>               
C          |      +---<PPARA>                                          
C          |             +---<ARCOEF>                                   
C          |---<FUNCSA>                                                 
C          +---<PLOTDD>                                                 
C                 |---<MAXMIN>                                          
C                 |---<XYAXIS>                                          
C                 +---<PLOTD>                                           
C                                                                       
C     THE FOLLOWING CONTROL PARAMETERS ARE PRESET AS DEFAULT OPTION     
C          M1     = 2    :trend order(0, 1, 2 or 3)                     
C          M2     = 0    :AR order (less than 11, try 2 first)          
C          PERIOD = 12   :number of seasons in one period               
C          SORDER = 1    :seasonal order (0, 1 or 2)                    
C          TRADE  = 0    :trading day adjustment (if TRADE = 1)         
C          MT     = 1    :original data input device, MT = 5:card reader
C          BSPAN  = 300  :maximum data length in filtering              
C          ISPAN  = 100  :number of data for backward filtering         
C          MISING                                                       
C          TAU2(I)       :system noise variances (i=1,2,3)              
C          PAC(I)        :PARCOR (i=1,...,M2)                           
C          IPR    = 7    :print out control                             
C          IDIF   = 1    :numerical differencing (1 sided or 2 sided)   
C          LOG    = 0    :log transformation of data (if LOG = 1)       
C          YEAR          :the first year of the data                    
C          MESH   = 1    :draw mesh on the figure (if MESH > 0)         
C                                                                       
C     These options can be changed by using NAMELIST 'PARAM'            
C                                                                       
C          EXAMPLE:                                                     
C             &PARAM M2=2,LOG=1,&END                                    
C                                                                       
C                                                                       
C     -----  WRITTEN BY GENSHIRO KITAGAWA  ----END S                    
C                                                                       

cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
c
cc      DIMENSION  DATA(N) ,IPAR(11) 
cc      real*4     title(20)
cc      REAL*8     A(40), YMEAN ,para(26)
      PARAMETER (NIP=9, NPA=26)
cxx      DIMENSION  DATA(N) ,IPAR(NIP)
cxx      DIMENSION  DATA(N)
cxx      REAL*8     A(L+M2), YMEAN ,para(NPA)
cxx      REAL*8     TREND(N),SEASNL(N),AR(N),TRAD(N),NOISE(N)
cxx      INTEGER    IMIS(N), PERIOD, SORDER
cc      COMMON     /COMSM1/  WORK(300000)
cxx      DIMENSION  Z(N), E(L,LM1,N), TDAY(N,7)
      INTEGER :: N, iopt, imiss, LM1, ier 
      REAL(8) :: DATA(N), TREND(N), SEASNL(N), AR(N), TRAD(N), NOISE(N),
     1           para(NPA), omaxx
      INTEGER :: IMIS(N), PERIOD, SORDER
      REAL(8) :: A(L+M2), YMEAN, Z(N), E(L,LM1,N), TDAY(N,7), FF
cc      COMMON     /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,
cc     *               NS, NI, MISING, IOUT, LL, NN, NYEAR,nmonth, 
cc     *                     NPREDS, NPREDE, IPRED 
cc      COMMON     /COMSM4/  NP1, NP2, NP3, NP4, NP5, NP6, NP7, NP8       
c      COMMON     /CCC/     ISW, IPR, ISMT, IDIF, LOG                    
cxx      COMMON     /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
cc      common     /cccout/  IMIS(3000)
      COMMON     /COMSM2/  M1, M2, M3, M4, M5, M, L, PERIOD, SORDER,
     *                     NYEAR, nmonth
      COMMON     /CCC/     ISW, ISMT, IDIF, LOG, MESH              
C                                                                       
C           ...  set control parameters  ...                            
C
cc      NN = N
cc      CALL  SPARAM( MT,A,IPAR,para,iopt )
      LLL = L+M2
cxx      CALL  SPARAM( N,A,LLL,IPAR,NIP,para,NPA,iopt )
      CALL  SPARAM( A,LLL,para,NPA,iopt )
c      outmax = omaxx                      
      do 123 i=1,n                      
      imis(i) = 0
cc      if(data(i) .gt. omaxx)  imis(i) = 1
      if (imiss .gt. 0) then
         if(data(i) .gt. omaxx)  imis(i) = 1
      else if (imiss .lt. 0) then
         if(data(i) .lt. omaxx)  imis(i) = 1
      end if
 123  continue
C                                                                       
C           ...  read original data  ...                                
C                                                                       
cc      CALL  REDATA( MT,DATA,ISW,IPR,WORK,N,TITLE,YMEAN )                     
cxx      CALL  REDATAD( DATA,ISW,IPR,Z,N,YMEAN )
      CALL  REDATAD( DATA,ISW,Z,N,YMEAN )
C                                                                       
C           ...  allocation of working area  ...                        
C                                                                       
cc      CALL  AREA( 300000,N )
cc      if(N .le. 0) then
cc      para(3) = -1.0
cc      return
cc      end if                                           
C                                                                       
C           ...  log transformation  ...                                
C                                                                       
c      IF(LOG .EQ. 1)  CALL  LOGTRF( WORK,N )                            
cc      CALL  LOGTRF( WORK,IMIS,N,LOG )
cx      CALL  LOGTRF( Z,IMIS,N,LOG )
      CALL  LOGTRF( Z,IMIS,N,LOG,ier )
      if(ier .ne. 0) return
C                                                                       
C           ...  prepare calender for trading day adjustment  ...       
C                                                                       
      IF( M4 .NE. 0)  then
cc      if( isea .eq. 12 )  CALL  TRADE( NYEAR,nmonth,N,WORK(NP2) )
cc      if( isea .eq. 4 )  CALL  TRADE2( NYEAR,nmonth,N,WORK(NP2) )
      if( period .eq. 12 )  CALL  TRADE( NYEAR,nmonth,N,TDAY )
      if( period .eq. 4 )  CALL  TRADE2( NYEAR,nmonth,N,TDAY )
      endif
C                                                                       
C           ...  estimation of parameters  ...                          
C                                                                           

cc      CALL  EPARAM( A,TITLE ,iopt)
      CALL  EPARAM( Z,E,TDAY,IMIS,N,A,iopt)
C                                                                       
C      write(6,10)
C 10   format( '*****************UJHHBLB*************')
C           ...  smoothing with estimated parameters  ...               
C                                                                       
      ISMT = 1                                                          
cc      LLL = L+M2
cc      CALL  FUNCSA( 1,A,FF,IFG )                                        
       CALL  FUNCSA( Z,E,TDAY,IMIS,N,LM1,LLL,A,FF,IFG)
C                                                                       
C           ...  plot estimated components  ...                         
C                                   

cc      call trpar( A,para)
      call trpar( A,LLL,para,NPA)

C      para(2)= FF
cc      CALL PLOTDD(WORK(NP1),L+M+1,TITLE,A,WORK(NP2),WORK(NP3),WORK(NP7),                                   
cc     *           TREND,SEASNL,AR,TRAD,NOISE)
      CALL PLOTDD(N,Z,E,LM1,TDAY,TREND,SEASNL,AR,TRAD,NOISE)
C     CALL PRINT( WORK(NP1),L+M+1,WORK(NP7) )                           
C                                                                       
c	do 234 i=1,n
c 234	ar(i) = imis(i)
c
c       ar(1) = omaxx

      RETURN                         
      E N D                                                             
cc      SUBROUTINE  ARCOEF( PAC,K,AR )                                    
      SUBROUTINE  ARCOEFD( PAC,K,AR )                                   
C                                                                       
C  ...  TRANSFORMATION FROM PARCOR TO AR COEFFICIENTS  ...              
C                                                                       
C       INPUTS:                                                         
C         PAC:  VECTOR OF PARTIAL AUTOCORRELATIONS                      
C         K:    ORDER OF THE MODEL                                      
C                                                                       
C       OUTPUTS:                                                        
C         AR:   VECTOR OF AR-COEFFICIENTS                               
C                                                                       
C                                                                       
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cc      DIMENSION  AR(K) , PAC(K) , W(100)
cxx      DIMENSION  AR(K) , PAC(K) , W(K)
      INTEGER :: K
      REAL(8) :: PAC(K), AR(K), W(K)
C                                                                       
      DO  30     II=1,K                                                 
      AR(II) = PAC(II)                                                  
      W(II)  = PAC(II)                                                  
      IM1 = II - 1                                                      
      IF( IM1 .LE. 0 )     GO TO 30                                     
      DO  10     J=1,IM1                                                
      JJ = II - J                                                       
cxx   10 AR(J) = W(J) - PAC(II)*W(JJ)                                      
      AR(J) = W(J) - PAC(II)*W(JJ)
   10 CONTINUE
      IF( II .EQ. K )     GO TO 40                                      
      DO  20     J=1,IM1                                                
c   20 W(J) = PAC(J)  # modified  97/10/17                            
cxx   20 W(J) = ar(J)                                                     
      W(J) = ar(J)
   20 CONTINUE
   30 CONTINUE                                                          
   40 CONTINUE                                                          
      RETURN                                                            
      E N D                                                             
cc      SUBROUTINE  EPARAM( A,TITLE ,iopt)                                     
      SUBROUTINE  EPARAM( Z,E,TDAY,IMIS,N,A,iopt)
C                                                                       
C  ...  Estimation of parameters  ...                                   
C                                                                       
cxx      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      DIMENSION  A(40), AI(20)                                          
cxx      DIMENSION  A(L+M2), AI(L+M2)
cxx      DIMENSION  Z(N), E(L,L+M+1,N), TDAY(N,7), IMIS(N)
cxx      INTEGER    PERIOD, SORDER
      INTEGER :: IMIS(N), N, iopt
      REAL(8) :: Z(N), E(L,L+M+1,N), TDAY(N,7), A(L+M2)
      INTEGER :: PERIOD, SORDER
      REAL(8) :: AI(L+M2), 
     1           F1, F2, F3, A1, A2, A3, DI, UI, TDF,
     2           DJACOB, FC, SIG2, AIC, FI, SIG2I, AICI, GI, GC
cc      REAL*4     CPUS, CPUE, TITLE(20), TIME(3)
cc      COMMON     /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,        
cc     *                  NS, NI, MISING, IOUT,LL,N,NYEAR,nmonth,NPREDS,      
cc     *                     NPREDE, NPRED
ccx      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)
cc      COMMON    /CMFUNC/  DJACOB,F,SIG2,AIC,FI,SIG2I,AICI,GI(20),G(20) 
c      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG                   
ccx      COMMON    /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
cxx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH
      COMMON    /COMSM2/  M1, M2, M3, M4, M5, M, L, PERIOD, SORDER,
     *                    NYEAR, nmonth
      COMMON  /COMSM3/  F1(10), F2(10), F3(300), A1(10), A2(10), A3(300)
     *                   ,DI, UI(3), TDF(7)                            
      COMMON /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(200),GC(200)
      COMMON    /CCC/     ISW, ISMT, IDIF, LOG, MESH
      EXTERNAL   FUNCSA                                                 
C                                                                       
      LM = L + M2                                                       
      ISMT = 0                                                          
C      CALL  CLOCK( CPUS,3 )                                             
      DO 10 I=1,LM                                                      
cxx   10 AI(I) = A(I)                                                      
      AI(I) = A(I)
   10 CONTINUE
C                                                                       
      CALL  SETFGH                                                      
      if(iopt .ge. 0) then
cc      CALL  OPTMIZ( FUNCSA,A,L+M2 )
      CALL  OPTMIZ( FUNCSA,Z,E,TDAY,IMIS,N,A,LM,L,L+M+1 )               
      end if
C                                                                       
C      CALL  CLOCK( CPUE,4 )                                             
C      CALL  DATE( DAY )                                                 
C      CALL  CLOCK( TIME,1 )                                             
C                                                                       
c      WRITE(6,650)                                                      
c      WRITE(6,600)  (TITLE(I),I=1,10), N                                
c      IF(LOG .EQ. 1)  WRITE(6,660)  DJACOB                              
c      WRITE(6,610)  M1, M2, M3, M4, M5                                  
C      WRITE(6,620)  DAY, (TIME(I),I=1,2)                                
C                                                                       
c      WRITE(6,630)                                                      
cc      CALL  PPARA( AI,GI,L,M2,FI,AICI,SIG2I )                          
C                                                                       
c      WRITE(6,640)                                                      
cc      CALL  PPARA( A,G,L,M2,F,AIC,SIG2 )                               
cc      CALL  PPARA( A,GC,L,M2,FC,AIC,SIG2 )                               
cc      IF( M4 .EQ. 0 )  GO TO 30                                         
c      WRITE(6,680)                                                      
c      DO 20 I=1,7                                                       
c   20 WRITE(6,690)  I, TDF(I)                                           
cc   30 CONTINUE                                                          
C      CPU = CPUE - CPUS                                                 
C      WRITE(6,670)   CPU                                                
c      WRITE(6,650)                                                      
C                                                                       
      RETURN                                                            
cxx  600 FORMAT( ////10X,'---  DATA  ---',/,10X,10A4,5X,'N =',I4 )
cxx  610 FORMAT( /10X,'---  MODEL  ---',/,10X,'M1 =',I2,5X,'M2 =',I2,5X,   
cxx     *        'M3 =',I3,5X,'M4 =',I2,5X,'M5 =',I2 )                     
cxx  620 FORMAT( /,10X,'---  PROGRAM  DECOMP  ---',/,10X,'DATE: ',A8,5X,   
cxx     *        'TIME: ',3A4 )                                            
cxx  630 FORMAT( ///15X,'<<<  INITIAL ESTIMATES  >>>' )                    
cxx  640 FORMAT( ///15X,'<<<  FINAL ESTIMATES  >>>' )                      
cxx  660 FORMAT( 10X,'***  LOG TRANSFORMED  ***   LOG-JACOBIAN =',         
cxx     *        F12.4 )                                                   
cxx  650 FORMAT( 1H1 )                                                     
cxx  670 FORMAT( //,10X,'CPU TIME =',F10.2 )                               
cxx  680 FORMAT( //,T20,'TRADING DAY',/,T15,'I',T22,'FACTOR' )             
cxx  690 FORMAT( 14X,I1,D16.8 )                                            
      E N D                                                             
cc      SUBROUTINE  FUNCND( FUNCT,M,A,F,G,IFG )                           
      SUBROUTINE  FUNCND( FUNCT,Z,E,TDAY,IMIS,N,M,A,F,G,IFG,L,LM1 )
C                                                                       
C  ...  FUNCTION EVALUATION AND NUMERICAL DIFFERENCING  ...             
C                                                                       
cxx      IMPLICIT   REAL*8( A-H,O-Z )
cxx      DIMENSION  Z(N), E(L,LM1,N), TDAY(N,7), IMIS(N)
cc      DIMENSION  A(M) , G(M) , B(20)
cxx      DIMENSION  A(M) , G(M) , B(M)
      INTEGER :: IMIS(N), N, M, IFG, L, LM1
      REAL(8) :: Z(N), E(L,LM1,N), TDAY(N,7), A(M), F, G(M)
      REAL(8) :: B(M), CONST, FB, FF,
     1           DJACOB, FC, SIG2, AIC, FI, SIG2I, AICI, GI, GC
c      COMMON     / CCC /  ISW , IPR, ISMT, IDIF, log
cxx      COMMON     /CCC/    ISW, IPR, ISMT, IDIF, LOG, MESH              
ccx      COMMON     /CMFUNC/ DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
      COMMON     /CCC/    ISW, ISMT, IDIF, LOG, MESH              
      COMMON   /CMFUNC/ DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(200),GC(200)

      DATA       ICNT /0/                                               
      CONST = 0.0001D0                                                  
C
cc      CALL  FUNCT( M,A,F,IFG )                                          
      CALL  FUNCT( Z,E,TDAY,IMIS,N,LM1,M,A,F,IFG )                      
      FB = F                                                            
      IF( ISW .GE. 1 )   RETURN                                         
C                                                                       
c      WRITE( 6,600 )   (A(I),I=1,M)                                     
      DO 10  I=1,M                                                      
cxx   10 B(I) = A(I)                                                       
      B(I) = A(I)
   10 CONTINUE
C                                                                       
      DO 30  II=1,M                                                     
      B(II) = A(II) + CONST                                             
cc      CALL  FUNCT( M,B,FF,IFG )                                         
      CALL  FUNCT( Z,E,TDAY,IMIS,N,LM1,M,B,FF,IFG )
      IF( IDIF .EQ. 1 )  GO TO 20                                       
      B(II) = A(II) - CONST                                             
cc      CALL  FUNCT( M,B,FB,IFG )
      CALL  FUNCT( Z,E,TDAY,IMIS,N,LM1,M,B,FB,IFG )                    
   20 G(II) = (FF-FB)/(CONST*IDIF)                                      
      IF( G(II) .GT. 1.0D20 )  G(II) = (F-FB)/CONST                     
      IF( G(II) .LT.-1.0D20 )  G(II) = (FF-F)/CONST                     
      IF( FB.GT.F .AND. FF.GT.F )  G(II) = 0.0D0                        
cxx   30 B(II) = A(II)                                                     
      B(II) = A(II)
   30 CONTINUE
C                                                                       
c      WRITE( 6,610 )   (G(I),I=1,M)                                     
      DO 40 I=1,M                                                       
cxx   40 GC(I) = G(I)                                                      
      GC(I) = G(I)
   40 CONTINUE
      ICNT = ICNT + 1                                                   
      IF(ICNT .GT. 1)  RETURN                                           
C                                                                       
      AICI  = AIC                                                       
      SIG2I = SIG2                                                      
      FI    = FC                                                        
      DO 50 I=1,M                                                       
cxx   50 GI(I) = G(I)                                                      
      GI(I) = G(I)
   50 CONTINUE
      RETURN                                                            
cxx  600 FORMAT( 3X,'---  PARAMETER  ---',(/,3X,5D13.5) )                  
cxx  610 FORMAT( 3X,'---  GRADIENT  ---',(/,3X,5D13.5) )                   
      E N D                                                             
cc      SUBROUTINE  FUNCSA( KK,A,FF,IFG )                                 
      SUBROUTINE  FUNCSA( Z,E,TDAY,IMIS,N,LM1,KK,A,FF,IFG )
C                                                                       
C  ...  Initial setting, filtering and smoothing  ...                   
C                                                                       
cxx      IMPLICIT REAL*8( A-H,O-Z )                                        

cc      DIMENSION  Y(40), S(40,40), R(40,40), A(KK)
cc      REAL*8     ZZ(3000)
cc      dimension   imisr(3000)
cxx      DIMENSION  Z(N), E(L,LM1,N), TDAY(N,7), IMIS(N), A(KK)
C
cxx      DIMENSION  IMISR(N), ZZ(N), Y(M), R(LM1+1,LM1), S(M+1,M+1)
cxx      INTEGER    PERIOD, SORDER
      INTEGER :: IMIS(N), N, LM1, KK, IFG 
      REAL(8) :: Z(N), E(L,LM1,N), TDAY(N,7), A(KK), FF
      INTEGER :: PERIOD, SORDER, IMISR(N)
      REAL(8) :: ZZ(N), Y(M), R(LM1+1,LM1), S(M+1,M+1),
     1           F1, F2, F3, A1, A2, A3, DI, UI, TDF, TAU2, SUM
cc      COMMON    /COMSM1/  WORK(300000)
cc      COMMON    /COMSM2/  K1, K2, K3, K4, K5, M, L, ISEA, KSEA,
cc     *              NS, NI, MISING, IOUT, LL, N, NYEAR, nmonth,NPREDS,  
cc     *                    NPREDE, NPRED
ccx      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)
cc      COMMON    /COMSM4/  NP1, NP2, NP3, NP4, NP5, NP6, NP7, NP8
c      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG       
cxx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
cc      common    /cccout/  IMIS(3000)      
      COMMON    /COMSM2/  M1, M2, M3, M4, M5, M, L, PERIOD, SORDER,
     *                    NYEAR, nmonth
      COMMON  /COMSM3/  F1(10), F2(10), F3(300), A1(10), A2(10), A3(300)
     *                   ,DI, UI(3), TDF(7)
      COMMON    /CCC/     ISW, ISMT, IDIF, LOG, MESH              
C                                                                       
      IFG = 0                                                           
cc      K = K1 + K2 + K3 + K4 + K5                                        
cc      K12 = K1 + K2
cc      LK1 = L + K + 1                                    
      M12 = M1 + M2
cc      LM1 = L + M + 1                                                   
cc      MJ = 40                                                           
      MJ=M+1


C                                                                       
C  ...  Backward Filtering  ...                                         
C                                                                       
cc      if(ni .gt. n) ni=n
      ni=n
      DO 10  I=1,NI                                                     
      imisr(i) = imis(NI-I+1)
cc   10 ZZ(I) = WORK(NI-I+1)                                              
cxx   10 ZZ(I) = Z(NI-I+1)                                              
      ZZ(I) = Z(NI-I+1)
   10 CONTINUE
      DO 20  I=1,L                                                      
      TAU2 = 0.5D0*(1.0D0 + DSIN( A(I) ))  + 0.0001
      IF(TAU2 .LT. 1.0D-20)  TAU2 = 1.0D-20                             
cxx   20 UI(I) = 1.0D0/DSQRT( TAU2 )                                     
      UI(I) = 1.0D0/DSQRT( TAU2 )
   20 CONTINUE
cc      IF( K2 .EQ. 0 )  GO TO 40
cc      DO 30  I=1,K2
      IF( M2 .EQ. 0 )  GO TO 40
      DO 30  I=1,M2                                                     
cxx   30 Y(I) = 0.90D0*DSIN( A(L+I) )                                      
      Y(I) = 0.90D0*DSIN( A(L+I) )
   30 CONTINUE
cc      CALL  ARCOEF( Y,K2,A2 )                                           
      CALL  ARCOEFD( Y,M2,A2 )
C                                                                       
cxx   40 DO 50  I=1,MJ
cxx      DO 50  J=1,MJ                                                     
cc      R(I,J) = 0.0D0                                                    
cxx   50 S(I,J) = 0.0D0
cxx      DO 55 I = 1,LM1+1      
cxx         DO 55 J = 1,LM1
cxx   55 R(I,J) = 0.0D0
   40 CONTINUE
      S(1:MJ,1:MJ) = 0.0D0
      R(1:LM1+1,1:LM1) = 0.0D0
      DI = 1.0D0                                                        
cc   70 IF( K2 .EQ. 0 )  GO TO 90                                         
cc      F2(1) = 1.0D0/A2(K2)                                              
cc      IF(K2.EQ.1)  GO TO 90                                             
cc      DO 80  I=2,K2                                                     
cc   80 F2(I) = -A2(I-1)/A2(K2)                                           
      IF( M2 .EQ. 0 )  GO TO 90                                         
      F2(1) = 1.0D0/A2(M2)                                              
      IF(M2.EQ.1)  GO TO 90                                             
      DO 80  I=2,M2                                                     
cxx   80 F2(I) = -A2(I-1)/A2(M2)                                           
      F2(I) = -A2(I-1)/A2(M2)
   80 CONTINUE
   90 CONTINUE                                                          
C                                                                       
cc      CALL  SMOTH3( ZZ,R,S,WORK(NP1),WORK(NP2),WORK(NP3),               
cc     *              WORK(NP4),WORK(NP5),IMISR,NI,LK1,MJ,FF,0,0 )
c     *              WORK(NP4),WORK(NP5),WORK(NP6),NI,LK1,MJ,FF,0,0 )
      CALL  SMOTH3( ZZ,R,S,E,TDAY,IMISR,N,LM1,FF,0,0 )

C                                                                       
C  ...  Transformation of State Vector  ...                             
C                                                                       
cc      CALL  RECOEF( S,K,K,MJ,Y )                                        
cc      CALL  STATE( Y,A1,K1 )                                            
cc      CALL  STATE( Y(K1+1),A2,K2 )                                      
cc      CALL  STATE( Y(K12+1),A3,K3 )                                     
      CALL  RECOEF( S,M,M,M+1,Y )
      CALL  STATE( Y,A1,M1 )
      CALL  STATE( Y(M1+1),A2,M2 )
      CALL  STATE( Y(M12+1),A3,M3 )
C                                                                       
cc      DO 210  I=1,K                                                     
      DO 210  I=1,M
      SUM = 0.0D0                                                       
cc      DO 220  J=I,K                                                     
      DO 220  J=I,M
cxx  220 SUM = SUM + S(I,J)*Y(J)
      SUM = SUM + S(I,J)*Y(J)
  220 CONTINUE
cc  210 S(I,K+1) = SUM
cxx  210 S(I,M+1) = SUM
      S(I,M+1) = SUM
  210 CONTINUE
C                                                                       
C  ...  Forward Filtering and/or Smoothing  ...                         
C                                                                       
cc      CALL  SMOTH3( WORK,R,S,WORK(NP1),WORK(NP2),WORK(NP3),
cc     *              WORK(NP4),WORK(NP5),IMIS,N,LK1,MJ,FF,1,ISMT )  
c     *              WORK(NP4),WORK(NP5),WORK(NP6),N,LK1,MJ,FF,1,ISMT )       
      CALL  SMOTH3( Z,R,S,E,TDAY,IMIS,N,LM1,FF,1,ISMT )
C                                                                       
      FF = -FF                                                          
      RETURN                                                            
C                                                                       
      E N D                                                             
      SUBROUTINE  HUSHL4( X,MJ1,N,K,M,ISW )                             
C                                                                       
C          HOUSEHOLDER TRANSFORMATION;   TYPE 4                         
C                                                                       
cxx      IMPLICIT  REAL*8( A-H,O-Z )                                       
cc      DIMENSION  X(MJ1,K), D(50)                                        
cxx      DIMENSION  X(MJ1,K), D(K)
      INTEGER :: MJ1, N, K, M, ISW
      REAL(8) :: X(MJ1,K)
      REAL(8) :: D(K), TOL, D0, D1, H, G, S
      DATA     TOL /1.0D-30/                                            
C                                                                       
      IF( ISW .EQ. 1 )   GO TO 110                                      
cxx      DO 100  II=M,K                                                    
      DO 101  II=M,K
      D0 = X(II,II)                                                     
      D1 = X(N,II)                                                      
      H = D0**2 + D1**2                                                 
      IF( H .GT. TOL )  GO TO 20                                        
      G = 0.0D0                                                         
      GO TO  100                                                        
   20 G = DSQRT( H )                                                    
      IF( D0 .GE. 0.0D0 )   G = -G                                      
      H = H - D0*G                                                      
      D0 = D0 - G                                                       
      D(II) = D0                                                        
C                                                                       
      IF( II .EQ. K )  GO TO 100                                        
      DO 60  J=II+1,K                                                   
      S = (D0*X(II,J) + D1*X(N,J))/H                                    
      X(II,J) = X(II,J) - D0*S                                          
      X(N,J) = X(N,J) - D1*S                                            
   60 CONTINUE                                                          
  100 X(II,II) = G                                                      
  101 CONTINUE
      RETURN                                                            
C                                                                       
C                                                                       
  110 CONTINUE                                                          
      DO 120  J=M,K-1                                                   
      S = (D(J)*X(J,K)+X(N,J)*X(N,K))                                   
      H = -X(J,J)*D(J)                                                  
      S = S/H                                                           
      X(J,K) = X(J,K)-D(J)*S                                            
cxx  120 X(N,K) = X(N,K)-X(N,J)*S
      X(N,K) = X(N,K)-X(N,J)*S
  120 CONTINUE
      RETURN                                                            
C                                                                       
      E N D                                                             
      SUBROUTINE  HUSHL7( X,D,MJ1,K,M,KE )                              
C                                                                       
C     Householder Transformation,  TYPE 7                               
C                                                                       
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cxx      DIMENSION  X(MJ1,K) , D(MJ1)                                      
      INTEGER :: MJ1, K, M, KE
      REAL(8) :: X(MJ1,K), D(MJ1)
      REAL(8) :: TOL, H, G, F, S
C                                                                       
           TOL = 1.0D-30                                                
C                                                                       
cxx      DO 100  II=1,KE                                                   
      DO 101  II=1,KE
      NE = MAX(M,II) + 1                                                
         H = 0.0D00                                                     
         DO 10  I=II,NE                                                 
            D(I) = X(I,II)                                              
cxx   10       H = H + D(I)*D(I)
            H = H + D(I)*D(I)
   10    CONTINUE
         IF( H .GT. TOL )  GO TO 20                                     
         G = 0.0D00                                                     
         GO TO 100                                                      
   20    G = DSQRT( H )                                                 
         F = X(II,II)                                                   
         IF( F .GE. 0.0D00 )   G = -G                                   
         D(II) = F - G                                                  
         H = H - F*G                                                    
C                                                                       
C          FORM  (I - D*D'/H) * X, WHERE H = D'D/2                      
C                                                                       
         II1 = II+1                                                     
      DO 30 I=II1,NE                                                    
cxx   30 X(I,II) = 0.0D0
      X(I,II) = 0.0D0
   30 CONTINUE
         IF( II .EQ. K )  GO TO 100                                     
C                                                                       
         DO 60  J=II1,K                                                 
            S = 0.0D00                                                  
            DO 40  I=II,NE                                              
cxx   40       S = S + D(I)*X(I,J)                                         
            S = S + D(I)*X(I,J)
   40       CONTINUE
            S = S/H                                                     
            DO 50  I=II,NE                                              
cxx   50      X(I,J) = X(I,J) - D(I)*S                                     
            X(I,J) = X(I,J) - D(I)*S
   50       CONTINUE
   60    CONTINUE                                                       
  100 X(II,II) = G                                                      
  101 CONTINUE
C                                                                       
      RETURN                                                            
C                                                                       
      E N D                                                             
      FUNCTION  ID( K )                                                 
C                                                                       
C  ...  ID = 1:    IF K > 0                                             
C       ID = 0:    OTHERWISE                                            
C                                                                       
      ID = 0                                                            
      IF( K .GT. 0 )  ID = 1                                            
      RETURN                                                            
      E N D                                                             
cc      SUBROUTINE  LINEA1( FUNCT,X,H,RAM,EE,G,K,IG )                     
      SUBROUTINE  LINEA1( FUNCT,Z,E,TDAY,IMIS,N,L,LM1,
     *                    X,H,RAM,EE,G,K,IG )
C                                                                       
C  ...  LINE SEARCH (WOLFE'S ALGORITHM)  ...                            
C                                                                       
cxx      IMPLICIT  REAL  *8 ( A-H,O-Z )                                    
cxx      INTEGER  RETURN,SUB                                               
cc      DIMENSION  X(K) , H(K) , X1(20)                                   
cc      DIMENSION  G(20)                                                  
cxx      DIMENSION  Z(N), E(L,LM1,N), TDAY(N,7), IMIS(N)
cxx      DIMENSION  X(K) , H(K) , G(K), X1(K)
      INTEGER :: IMIS(N), N, L, LM1, K, IG
      REAL(8) :: Z(N), E(L,LM1,N), TDAY(N,7), X(K), H(K), RAM,
     1           EE, G(K)
      INTEGER :: RETURN, SUB
      REAL(8) :: X1(K), C1, C2, F0, SUM0, CONST2, HNORM,
     1           RAM1, RAM2, RAM3, A1, A2, A3, B1, B2, E1, E2,
     2           E3, H1, H2, SUM
c      COMMON     / CCC /  ISW , IPR                          
c      COMMON     /CCC/     ISW, IPR, ISMT, IDIF, LOG
cxx      COMMON     /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH
      COMMON     /CCC/     ISW, ISMT, IDIF, LOG, MESH
      EXTERNAL   FUNCT              
      C1 = 0.01D0                                                       
      C2 = 0.5D0                                                        
C
cc
      IFG = 0
cc
      F0 = EE                                                           
      SUM0 = 0.0D0                                                      
      DO 15 I=1,K                                                       
cxx   15 SUM0 = SUM0 + G(I)*H(I)                                           
      SUM0 = SUM0 + G(I)*H(I)
   15 CONTINUE
      ISW = 1                                                           
      RAM = 0.5D0                                                       
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
      IF( RAM2*HNORM .GT. 5.0D0 )  GO TO 48                             
      DO 20  I=1,K                                                      
cxx   20 X1(I) = X(I) + RAM2*H(I)                                          
      X1(I) = X(I) + RAM2*H(I)
   20 CONTINUE
cc      CALL  FUNCND( FUNCT,K,X1,E2,G,IG )                                
      CALL  FUNCND( FUNCT,Z,E,TDAY,IMIS,N,K,X1,E2,G,IG,L,LM1 )
c      IF(IPR.GE.7)  WRITE(6,2)  RAM2,E2                                 
C                                                                       
      IF( IG .EQ. 1 )  GO TO  50                                        
      IF( E2 .GT. E1 )  GO TO 50                                        
CCCCCCCCCCCCCCCCCCCCCCC                                                 
      H1 = E2 - F0 - C1*RAM2*SUM0                                       
      SUM = 0.0D0                                                       
      DO 25 I=1,K                                                       
cxx   25 SUM = SUM + G(I)*H(I)                                             
      SUM = SUM + G(I)*H(I)
   25 CONTINUE
      H2 = SUM - C2*SUM0                                                
      IF(H1.LE.0.0D0 .AND. H2.GE.0.0D0)  RETURN                         
CCCCCCCCCCCCCCCCCCCCCCC                                                 
   30 RAM3 = RAM2*4.D0                                                  
      DO 40  I=1,K                                                      
cxx   40 X1(I) = X(I) + RAM3*H(I)
      X1(I) = X(I) + RAM3*H(I)
   40 CONTINUE
cc      CALL  FUNCND( FUNCT,K,X1,E3,G,IG )                                
      CALL  FUNCND( FUNCT,Z,E,TDAY,IMIS,N,K,X1,E3,G,IG,L,LM1 )
      IF( IG.EQ.1 )  GO TO  500                                         
c      IF( IPR.GE.7 )  WRITE(6,3)  RAM3,E3                               
      IF( E3 .GT. E2 )  GO TO 70                                        
      IF(RAM3.GT.1.0D10 .AND. E3.LT.E1)  GO TO 45                       
      IF(RAM3.GT.1.0D10 .AND. E3.GE.E1)  GO TO 46                       
      RAM1 = RAM2                                                       
      RAM2 = RAM3                                                       
      E1 = E2                                                           
      E2 = E3                                                           
      GO TO 30                                                          
C                                                                       
   45 RAM = RAM3                                                        
      EE = E3                                                           
      RETURN                                                            
C                                                                       
   46 RAM = 0.0D0                                                       
      RETURN                                                            
C                                                                       
   48 E2 = 1.0D30                                                       
C                                                                       
   50 RAM3 = RAM2                                                       
      E3 = E2                                                           
      RAM2 = RAM3*0.1D0                                                 
      IF( RAM2*HNORM .LT. CONST2 )  GO TO  400                          
      DO 60  I=1,K                                                      
cxx   60 X1(I) = X(I) + RAM2*H(I)                                          
      X1(I) = X(I) + RAM2*H(I)
   60 CONTINUE
cc      CALL  FUNCND( FUNCT,K,X1,E2,G,IG )                                
      CALL  FUNCND( FUNCT,Z,E,TDAY,IMIS,N,K,X1,E2,G,IG,L,LM1 )
c      IF(IPR.GE.7)  WRITE(6,4)  RAM2,E2                                 
      IF( E2.GT.E1 )  GO TO 50                                          
CCCCCCCCCCCCCCCCCCCCCCC                                                 
      H1 = E2 - F0 - C1*RAM2*SUM0                                       
      SUM = 0.0D0                                                       
      DO 65 I=1,K                                                       
cxx   65 SUM = SUM + G(I)*H(I)                                             
      SUM = SUM + G(I)*H(I)
   65 CONTINUE
      H2 = SUM - C2*SUM0                                                
      IF(H1.GT.0.0D0 .OR. H2.LT.0.0D0)  GO TO 70                        
      RAM = RAM2                                                        
      RETURN                                                            
CCCCCCCCCCCCCCCCCCCCCCC                                                 
C                                                                       
cc   70 ASSIGN 80 TO RETURN                                               
   70 RETURN = 80
      GO TO 200                                                         
C                                                                       
   80 DO 90  I=1,K                                                      
cxx   90 X1(I) = X(I) + RAM*H(I)                                           
      X1(I) = X(I) + RAM*H(I)
   90 CONTINUE
cc      CALL  FUNCND( FUNCT,K,X1,EE,G,IG )                                
      CALL  FUNCND( FUNCT,Z,E,TDAY,IMIS,N,K,X1,EE,G,IG,L,LM1 )
c      IF(IPR.GE.7)  WRITE(6,5)  RAM,EE                                  
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
      IF( SUB .EQ. 200 ) GO TO 200
      IF( SUB .EQ. 300 ) GO TO 300
C                                                                       
  100 RAM1 = RAM                                                        
      E1 = EE                                                           
cc      GO TO  SUB,( 200,300 )                                            
      IF( SUB .EQ. 200 ) GO TO 200
      IF( SUB .EQ. 300 ) GO TO 300
C                                                                       
  110 IF( EE .LE. E2 )  GO TO 120                                       
      RAM3 = RAM                                                        
      E3 = EE                                                           
cc      GO TO  SUB,( 200,300 )                                            
      IF( SUB .EQ. 200 ) GO TO 200
      IF( SUB .EQ. 300 ) GO TO 300
C                                                                       
  120 RAM1 = RAM2                                                       
      RAM2 = RAM                                                        
      E1 = E2                                                           
      E2 = EE                                                           
cc      GO TO  SUB,( 200,300 )                                            
      IF( SUB .EQ. 200 ) GO TO 200
      IF( SUB .EQ. 300 ) GO TO 300
C                                                                       
  130 DO 140  I=1,K                                                     
cxx  140 X1(I) = X(I) + RAM*H(I)                                           
      X1(I) = X(I) + RAM*H(I)
  140 CONTINUE
cc      CALL  FUNCND( FUNCT,K,X1,EE,G,IG )                                
      CALL  FUNCND( FUNCT,Z,E,TDAY,IMIS,N,K,X1,EE,G,IG,L,LM1 )
c      IF( IPR.GE.7 )  WRITE(6,6)  RAM,EE                                
cc      ASSIGN 200 TO SUB                                                 
      SUB = 200
      IFG = IFG+1                                                       
C     IFG = 0                                                           
C                                                                       
      H1 = EE - F0 - C1*RAM*SUM0                                        
      SUM = 0.0D0                                                       
      DO 145 I=1,K                                                      
cxx  145 SUM = SUM + G(I)*H(I)                                             
      SUM = SUM + G(I)*H(I)
  145 CONTINUE
      H2 = SUM - C2*SUM0                                                
      IF( H1.LE.0.0D0 .AND. H2.LE.0.0D0 )  GO TO 150                    
      IF( IFG .LE. 20 )  GO TO 95                                       
C                                                                       
  150 IF( E2 .LT. EE )  RAM = RAM2                                      
      RETURN                                                            
C                                                                       
C      -------  INTERNAL SUBROUTINE SUB1  -------                       
  200 IF( RAM3-RAM2 .GT. 5.0D0*(RAM2-RAM1) )  GO TO 202                 
      IF( RAM2-RAM1 .GT. 5.0D0*(RAM3-RAM2) )  GO TO 204                 
      A1 = (RAM3-RAM2)*E1                                               
      A2 = (RAM1-RAM3)*E2                                               
      A3 = (RAM2-RAM1)*E3                                               
      B2 = (A1+A2+A3)*2.D0                                              
      B1 = A1*(RAM3+RAM2) + A2*(RAM1+RAM3) + A3*(RAM2+RAM1)             
      IF( B2 .EQ. 0.D0 )  GO TO 210                                     
      RAM = B1 /B2                                                      
cc      GO TO RETURN ,( 80,130 )                                          
      IF( RETURN .EQ. 80 ) GO TO 80
      IF( RETURN .EQ. 130 ) GO TO 130
  202 RAM = (4.0D0*RAM2 + RAM3)/5.0D0                                   
cc      GO TO RETURN, (80,130)                                            
      IF( RETURN .EQ. 80 ) GO TO 80
      IF( RETURN .EQ. 130 ) GO TO 130
  204 RAM = (RAM1 + 4.0D0*RAM2)/5.0D0                                   
cc      GO TO RETURN, (80,130)                                            
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
      IF( RETURN .EQ. 80 ) GO TO 130
C ---------------------------------
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
cc      CALL  FUNCND( FUNCT,K,X1,E3,G,IG )                                
      CALL  FUNCND( FUNCT,Z,E,TDAY,IMIS,N,K,X1,E3,G,IG,L,LM1 )
c      IF( IPR.GE.7 )  WRITE(6,7)  RAM,E3                                
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
cxx    1 FORMAT( 5X ,'RAM =',D10.3, 5X,'E1 =',F10.3 )                      
cxx    2 FORMAT( 5X ,'RAM =',D10.3, 5X,'E2 =',F10.3 )                      
cxx    3 FORMAT( 5X ,'RAM =',D10.3, 5X,'E3 =',F10.3 )                      
cxx    4 FORMAT( 5X ,'RAM =',D10.3, 5X,'E4 =',F10.3 )                      
cxx    5 FORMAT( 5X ,'RAM =',D10.3, 5X,'E5 =',F10.3 )                      
cxx    6 FORMAT( 5X ,'RAM =',D10.3, 5X,'E6 =',F10.3 )                      
cxx    7 FORMAT( 5X ,'RAM =',D10.3, 5X,'E7 =',F10.3 )                      
      E N D                                                             
cc      SUBROUTINE  LOGTRF( Z,N,ilog )                                         
cx      SUBROUTINE  LOGTRF( Z,IMIS,N,ilog )
      SUBROUTINE  LOGTRF( Z,IMIS,N,ilog,ier )
C                                                                       
C ... LOG TRANSFORMATION ...                                            
C                                                                       
cxx      IMPLICIT REAL*8(A-H,O-Z)                                          
cxx      REAL*8   Z(N)
cxx      DIMENSION   IMIS(N)                                               
      INTEGER :: IMIS(N), N, ilog, ier
      REAL(8) :: Z(N), DJACOB, FC, SIG2, AIC, FI, SIG2I, AICI, GI, GC
C      COMMON     /CMFUNC/  DJACOB                                       
cc      COMMON     /CMFUNC/  DJACOB,F,SIG2,AIC,FI,SIG2I,AICI,GI(20),G(20)
cc      common     /cccout/  IMIS(3000)
ccx      COMMON    /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
      COMMON  /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(200),GC(200)
C                                                                       
      ier = 0
      DJACOB = 0.0D0                                                    
      if(ilog .eq. 0) return
      ier = -1
      DO 10 I=1,N                                                       
c      DJACOB = DJACOB - DLOG10(Z(I))                                    
cx	if(IMIS(i) .ne. 1)  DJACOB = DJACOB - DLOG(Z(I)) 
c   10 Z(I) = DLOG10( Z(I) )   
cx	Z(I) = DLOG( Z(I) )                                             
         if(IMIS(i) .ne. 1) then
            if(Z(I) .le. 0) return
            DJACOB = DJACOB - DLOG(Z(I)) 
            Z(I) = DLOG( Z(I) )
         end if
 10   continue
      ier = 0
C                                                                       
c      WRITE(6,600)  DJACOB                                              
C     WRITE(6,610) (Z(I),I=1,N)                                         
      RETURN                                                            
C                                                                       
cxx  600 FORMAT(1H0,'OPTION = 1, (LOG TRANSFORMATION),     LOG-JACOBIAN =',
cxx     1       D15.7)                                                     
cxx  610 FORMAT(1H ,10D13.6)                                               
      END
cc      SUBROUTINE  OPTMIZ( FUNCT,X,N )                                   
      SUBROUTINE  OPTMIZ( FUNCT,Z,E,TDAY,IMIS,NN,X,N,L,LM1 )
C                                                                       
C  ...  NUMERICAL OPTIMIZATION  ...                                     
C            LATEST REVISION:  JUNE 20, 1983                            
C                                                                       
cxx      IMPLICIT  REAL*8 (A-H,O-Z)
cc      DIMENSION  X(20) , DX(20) , G(20) , G0(20) , Y(20)                
cc      DIMENSION  H(20,20) , WRK(20) , S(20)
cxx      DIMENSION  Z(NN), E(L,LM1,NN), TDAY(NN,7), IMIS(NN), X(N)
cxx      DIMENSION  DX(N) ,G(N) ,G0(N) ,Y(N) ,H(N,N), WRK(N), S(N)
      INTEGER :: IMIS(NN), NN, N, L, LM1
      REAL(8) :: Z(NN), E(L,LM1,NN), TDAY(NN,7), X(N)
      REAL(8) :: DX(N), G(N), G0(N), Y(N), H(N,N), WRK(N), S(N),
     1           DJACOB, FC, SIG2, AIC, FI, SIG2I, AICI, GI, GC,
     2           TAU2, EPS1, EPS2, CONST1, XM, SUM, S1, S2, STEM,
     3           SS, DS2, GTEM, ED, RAMDA, XMB
c      COMMON     / CCC /  ISW, IPR                                      
c      COMMON     /CMFUNC/  DJACOB, F, SD, AIC                           
cc      COMMON     /CMFUNC/  DJACOB,F,SIG2,AIC,FI,SIG2I,AICI,
cc     *                   GI(20),GDUM(20) 
cxx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
ccx      COMMON    /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
      COMMON    /CCC/     ISW, ISMT, IDIF, LOG, MESH              
      COMMON  /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(200),GC(200)
      EXTERNAL   FUNCT
cc      DATA  TAU1 , TAU2  /  1.0D-1 , 1.0D-1  /             
      DATA  TAU2  /  1.0D-1  /
      DATA  EPS1 , EPS2  / 1.0D-2 , 1.0D-2  /                           
      CONST1 = 1.0D-50                                                  
      ISW = 0                                                           
C                                                                       
cc      CALL  FUNCND( FUNCT,N,X,XM,G,IG )
      XM = 0
      CALL  FUNCND( FUNCT,Z,E,TDAY,IMIS,NN,N,X,XM,G,IG,L,LM1 )
C                                                                       
C          INITIAL ESTIMATE OF INVERSE OF HESSIAN                       
C                                                                       
      ICOUNT=0

cxx 1000 DO 20  I=1,N                                                      
cxx      DO 10  J=1,N                                                      
cxx   10 H(I,J) = 0.0D0                                                    
cxx      S(I)   = 0.0D0                                                    
cxx      DX(I)  = 0.0D0                                                    
cxx   20 H(I,I) = 1.0D0
 1000 CONTINUE
      H(1:N,1:N) = 0.0D0                                                    
      S(1:N)     = 0.0D0                                                    
      DX(1:N)    = 0.0D0
      DO 20  I=1,N
      H(I,I) = 1.0D0
   20 CONTINUE
      ICC = 0                                                           

      ICOUNT=ICOUNT+1
      IF(ICOUNT .GT. 4) THEN
c         WRITE(6,700)
c         WRITE(6,701)
         RETURN
      ENDIF
C                                                                       
c      IF( IPR .GE. 2 )   WRITE( 6,340 )   XM                            
C                                                                       
C  ...  QUASI NEWTON ALGORITHM  ...                                     
C                                                                       
 2000 ICC = ICC + 1                                                     
      IF( ICC .EQ. 1 )   GO TO 120                                      
C                                                                       
      DO 40  I=1,N                                                      
cxx   40 Y(I) = G(I) - G0(I)                                               
      Y(I) = G(I) - G0(I)
   40 CONTINUE
      DO 60  I=1,N                                                      
      SUM = 0.0D0                                                       
      DO 50  J=1,N                                                      
cxx   50 SUM = SUM + Y(J)*H(I,J)                                           
      SUM = SUM + Y(J)*H(I,J)
   50 CONTINUE
cxx   60 WRK(I) = SUM                                                      
      WRK(I) = SUM
   60 CONTINUE
      S1 = 0.0D0                                                        
      S2 = 0.0D0                                                        
      DO 70  I=1,N                                                      
      S1 = S1 + WRK(I)*Y(I)                                             
cxx   70 S2 = S2 + DX(I) *Y(I)                                             
      S2 = S2 + DX(I) *Y(I)
   70 CONTINUE
      IF( S1.LE.CONST1 .OR. S2.LE.CONST1 )  GO TO 900                   
C                                                                       
C  ...  BFGS FORMULA FOR UPDATING INVERSE OF HESSIAN MATRIX  ...        
C                                                                       
      STEM = S1 / S2 + 1.0D00                                           
cxx      DO 110  I=1,N                                                     
      DO 111  I=1,N
      DO 110  J=I,N                                                     
      H(I,J) = H(I,J)- (DX(I)*WRK(J)+WRK(I)*DX(J)-DX(I)*DX(J)*STEM)/S2  
cxx  110 H(J,I) = H(I,J)                                                   
      H(J,I) = H(I,J)
  110 CONTINUE
  111 CONTINUE
C                                                                       
  120 SS = 0.0D0                                                        
      DO 150  I=1,N                                                     
      SUM = 0.0D0                                                       
      DO 140  J=1,N                                                     
cxx  140 SUM = SUM + H(I,J)*G(J)                                           
      SUM = SUM + H(I,J)*G(J)
  140 CONTINUE
      SS = SS + SUM * SUM                                               
cxx  150 S(I) = -SUM                                                       
      S(I) = -SUM
  150 CONTINUE
C                                                                       
      S1 = 0.0D0                                                        
      S2 = 0.0D0                                                        
      DO 170  I=1,N                                                     
      S1 = S1 + S(I)*G(I)                                               
cxx  170 S2 = S2 + G(I)**2                                                 
      S2 = S2 + G(I)**2
  170 CONTINUE
      DS2 = DSQRT(S2)                                                   
      GTEM = DABS(S1) / DS2                                             
C     IF( GTEM .LE. TAU1  .AND.  DS2 .LE. TAU2 )     GO TO  900         
      IF( S1 .GE. 0.0D0 )  GO TO  1000                                  
      ED = XM                                                           
C                                                                       
C  ...  LINE SEARCH  ...                                                
C                                                                       
cc      CALL  LINEA1( FUNCT,X,S,RAMDA,ED,G,N,IG )                         
      CALL LINEA1(FUNCT,Z,E,TDAY,IMIS,NN,L,LM1,X,S,RAMDA,ED,G,N,IG)
C                                                                       
c      IF( IPR .GE. 2 )   WRITE( 6,330 )   RAMDA, ED                     
C                                                                       
      S1 = 0.0D0                                                        
      DO 210  I=1,N                                                     
      DX(I) = S(I)*RAMDA                                                
      S1 = S1 + DX(I)**2                                                
      G0(I) = G(I)                                                      
cxx  210 X(I) = X(I) + DX(I)                                               
      X(I) = X(I) + DX(I)
  210 CONTINUE
      XMB  = XM                                                         
      ISW  = 0                                                          
C                                                                       
cc      CALL  FUNCND( FUNCT,N,X,XM,G,IG )
      CALL  FUNCND( FUNCT,Z,E,TDAY,IMIS,NN,N,X,XM,G,IG,L,LM1 )
C                                                                       
      S2 = 0.D0                                                         
      DO 220  I=1,N                                                     
cxx  220 S2 = S2 + G(I)**2
      S2 = S2 + G(I)**2
  220 CONTINUE
      IF( DSQRT(S2) .LT. TAU2 )   GO TO 900                             
      IF( XMB-XM .LT. EPS1  .AND.  DSQRT(S1) .LT. EPS2 )  GO TO 900     
      IF( XMB-XM .LT. 0.0001 .AND. ICC .GT. N ) GO TO 900
C      IF( ICC .GE. N*2 )  GO TO 1000                                    
      GO TO 2000                                                        
C                                                                       
  900 CONTINUE                                                          
      S2 = 0.0D0                                                        
      DO 230 I=1,N                                                      
cxx  230 S2 = S2 + G(I)**2                                                 
      S2 = S2 + G(I)**2
  230 CONTINUE
      IF( DSQRT(S2) .GT. 1.0D0 )  GO TO 1000                            
cxx      IF( IPR .LE. 0 )   RETURN                                         
c      WRITE(6,620)                                                      
c      WRITE( 6,600 )                                                    
c      WRITE( 6,610 )     (X(I),I=1,N)                                   
c      WRITE( 6,601 )                                                    
c      WRITE( 6,610 )     (G(I),I=1,N)                                   
      RETURN                                                            
cxx  330 FORMAT( 5X ,'RAM =',D10.3,5X,'E0 =',F10.3 )                       
cxx  340 FORMAT( 25X,'EE =',F10.3 )                                        
cxx  600 FORMAT( 1H0,'--  PARAMETER  ---' )                                
cxx  601 FORMAT( 1H0,'--  GRADIENT  --' )                                  
cxx  610 FORMAT( 1H ,10D13.5 )                                             
cxx  620 FORMAT( 10X,'<<<  FINAL ESTIMATES  >>>' )                         
cxx 700  FORMAT( ' OPTIMIZATION DOES NOT SUCCEED ' )
cxx 701  FORMAT( ' CHECK THE MODEL ! ' )
      E N D                                                             

cc      SUBROUTINE  PLOTDD( E,LM1,TITLE,A,TRADE,REG,Z ,COMP1,COMP2,COMP3,
      SUBROUTINE  PLOTDD( N,Z,E,LM1,TRADE,COMP1,COMP2,COMP3,
     *                    COMP4,COMP5)
cxx      IMPLICIT REAL*8(A-H,O-Z)                                          

C  ...  PLOT ESTIMATED PARAMETRS, ORIGINAL DATA AND ESTIMATED PARAMETERS
C
cxx      INTEGER    PERIOD, SORDER
      INTEGER :: N, LM1
      REAL(8) :: Z(N), E(L,LM1,N), TRADE(N,7),  COMP1(N), COMP2(N),
     1           COMP3(N), COMP4(N), COMP5(N)
      INTEGER :: PERIOD, SORDER
      REAL(8) :: REG(N,M5), SUM, tmp
C      COMMON   /COMSM1/  nwww, WORK                                        
cc      COMMON     /COMSM1/  WORK(300000)
cc      COMMON     /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,                  
cc     *               NS,NI,MISING,IOUT,LL,N,NYEAR,nmonth,NPS,NPE,NPRED
cxx      COMMON     /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
      COMMON     /COMSM2/  M1, M2, M3, M4, M5, M, L, PERIOD, SORDER,
     *                     NYEAR, nmonth
      COMMON     /CCC/     ISW, ISMT, IDIF, LOG, MESH              
cc      DIMENSION  E(L,LM1,N), Z(N), TRADE(N,7), REG(N,M5)
cxx      DIMENSION  E(L,LM1,N), Z(N), TRADE(N,7), REG(N,M5)
cxx      REAL*8     COMP1(N),COMP2(N),COMP3(N),COMP4(N),COMP5(N)
cc      real*4     TITLE(20)
cc      REAL*8     A(40),COMP1(N),COMP2(N),COMP3(N) 
cc      REAL*8     COMP4(N),COMP5(N)


C      CALL  PLTCK                                                       
C      CALL  PLOTS( 8HDATA,0.0 )                                         
C      CALL  NEWPEN( 2 )                                                 
C                                                                       
C      CALL  MAXMIN( WORK,N,Y0,Y1,DY )                                   
C                                                                       
C      XX = 1.5                                                          
C      YY = 15.0                                                         
C      X0 = 0                                                            
C      X1 = N                                                            
C      DX = 12.0                                                         
C      IF( N .GT. 200 )  DX = 24.0                                       
C      IF( N .GT. 400 )  DX = 48.0                                       
C      NS = 1                                                            
C      WX = N*0.1                                                        
C  150 IF(WX .LE. 13.0)  GO TO 160                                       
C      WX = WX *0.8                                                      
C      GO TO 150                                                         
C  160 CONTINUE                                                          
C      WY = 6.0                                                          
C      WC = 0.22                                                         
      M12 = M1 + M2                                                     
      M123= M1 + M2 + M3                                                
      M1234=M1 + M2 + M3 + M4                                           
C                                                                       
c      NNY = (Y1-Y0)/DY                                                  
c      IF( NNY .GE. 4 )  DY1 = DY                                        
c      IF( NNY .EQ. 2 .OR. NNY .EQ. 3 )  DY1 = DY*0.5                    
c      IF( NNY .EQ. 1 )  DY1 = DY*0.25                                   
c      WY1 = DY1*2.0/(Y1-Y0)*WY                                          
c      Y01 =-DY1                                                         
c      Y11 = DY1                                                         
C                                                                       
C      CALL  TITLEP( TITLE,A )                                           
C      CALL  PLOT( XX,YY,-3 )                                            
C                                                                       
      DO 5 I=1,N
         COMP1(I) = 0.0D0
         COMP2(I) = 0.0D0
         COMP3(I) = 0.0D0
         COMP4(I) = 0.0D0
    5 CONTINUE
C
C  ...  PLOT ORIGINAL DATA AND TREND  ...                               
C                                                                       
C      CALL  XYAXIS( 0.0,0.0,WX,WY,X0,X1,Y0,Y1,DX,DY,2,MESH )            
C      CALL  PLOTD( WORK,N,Y0,Y1,WX,WY,1,1,2 )                           
      DO 10 I=1,N                                                       
cxx   10 COMP1(I) = E(1,1,I)                                               
      COMP1(I) = E(1,1,I)
   10 CONTINUE
C      CALL  PLOTD( Z,N,Y0,Y1,WX,WY,1,1,2 )                              
C      CALL  SYMBOL( WX-WC*18,WY+0.1,WC,'ORIGINAL AND TREND',0.0,18 )    
C                                                                       
C  ...  PLOT SEASONAL COMPONENT  ...                                    
C                                                                       
      IF( SORDER .NE. 0) THEN
      DO 20 I=1,N                                                       
cxx   20 COMP2(I) = E(1,M1+M2+1,I)
      COMP2(I) = E(1,M1+M2+1,I)
   20 CONTINUE
      END IF
C      CALL  MAXMIN( Z,N,ZMIN,ZMAX,DZ )                                  
C      DY2 = DY                                                          
C      IF( NNY .EQ. 1 )  DY2 = DY*0.5                                    
C      Y02 =-DY2                                                         
C      Y12 = DY2                                                         
C      NNY = 2                                                           
C  170 IF( Y02 .LE. ZMIN )  GO TO 180                                    
C      Y02 = Y02 - DY2                                                   
C      NNY = NNY + 1                                                     
C      GO TO 170                                                         
C  180 IF( Y12 .GE. ZMAX )  GO TO 190                                    
C      Y12 = Y12 + DY2                                                   
C      NNY = NNY + 1                                                     
C      GO TO 180                                                         
C  190 CONTINUE                                                          
C      WY2 = DY2*NNY*WY/(Y1-Y0)                                          
C      CALL  PLOT( 0.0,-WY2-1.5,-3 )                                     
C      CALL  XYAXIS( 0.0,0.0,WX,WY2,X0,X1,Y02,Y12,DX,DY2,2,MESH )        
C      CALL  PLOTD( Z,N,Y02,Y12,WX,WY2,1,1,2 )                           
C      CALL  SYMBOL( WX-WC*8,WY2+0.1,WC,'SEASONAL',0.0,8 )               
C                                                                       
C  ...  PLOT OBSERVATIONAL NOISE  ...                                   
C                                                                       

      IF( M4 .EQ. 6 )  then
      DO 35 I=1,N                                                       
      SUM = 0.0D0                                                       
      DO 30 J=1,6                                                       
cxx   30 SUM = SUM + E(1,M123+J,N)*(TRADE(I,J)-TRADE(I,7))                 
      SUM = SUM + E(1,M123+J,N)*(TRADE(I,J)-TRADE(I,7))
   30 CONTINUE
cx   35 E(2,1,I) = SUM                                                    
      E(2,1,I) = SUM
   35 CONTINUE
      end if
      if(m4 .eq. 1) then
      DO 37 I=1,N
      tmp=trade(I,2)+trade(I,3)+trade(I,4)+trade(I,5)+trade(I,6)
cxx 37   e(2,1,i) = (trade(I,1)+trade(I,7)-0.4*tmp)*e(1,M123+1,n)
      e(2,1,i) = (trade(I,1)+trade(I,7)-0.4*tmp)*e(1,M123+1,n)
 37   CONTINUE
      end if


cxx   40 IF( M5 .EQ. 0 )  GO TO 60                                         
      IF( M5 .EQ. 0 )  GO TO 60                                         
      DO 55 I=1,N                                                       
      SUM = 0.0D0                                                       
      DO 50 J=1,M5                                                      
cxx   50 SUM = SUM + E(1,M1234+J,N)*REG(I,J)                               
      SUM = SUM + E(1,M1234+J,N)*REG(I,J)
   50 CONTINUE
cxx   55 E(2,2,I) = SUM                                                    
      E(2,2,I) = SUM
   55 CONTINUE
   60 CONTINUE                                                          
C      CALL  PLOT( 0.0,-WY1-1.5,-3 )                                     
C      CALL  XYAXIS( 0.0,0.0,WX,WY1,X0,X1,Y01,Y11,DX,DY1,2,MESH )        
      DO 70 I=1,N                                                       
cc      Z(I) = WORK(I) - E(1,1,I)*ID(M1) - E(1,M1+1,I)*ID(M2)             
cc     *     - E(1,M12+1,I)*ID(M3) - E(2,1,I)*ID(M4) - E(2,2,I)*ID(M5)    
cc 70   COMP5(I)=Z(I)
cxx 70   COMP5(I) = Z(I) - E(1,1,I)*ID(M1) - E(1,M1+1,I)*ID(M2)
      COMP5(I) = Z(I) - E(1,1,I)*ID(M1) - E(1,M1+1,I)*ID(M2)            
     *     - E(1,M12+1,I)*ID(M3) - E(2,1,I)*ID(M4) - E(2,2,I)*ID(M5)
   70 CONTINUE
C                                                                       
C      CALL  PLOTD( Z,N,Y01,Y11,WX,WY1,1,1,2 )                           
C      CALL  SYMBOL( WX-WC*5,WY1+0.1,WC,'NOISE',0.0,5 )                  
C                                                                       
      IF( M2 .EQ. 0 )  GO TO 100                                        
C     CALL  TITLEP( TITLE,A )                                           
C      CALL  PLOT( 16.0,WY1+WY2+3.0,-3 )                                 
C                                                                       
C  ...  PLOT AR PROCESS  ...                                            
C                                                                       
C      CALL  XYAXIS( 0.0,0.0,WX,WY2,X0,X1,Y02,Y12,DX,DY2,2,MESH )        
      DO 80 I=1,N                                                       
cxx   80 COMP3(I) = E(1,M1+1,I)                                            
      COMP3(I) = E(1,M1+1,I)
   80 CONTINUE
C      CALL  PLOTD( Z,N,Y02,Y12,WX,WY2,1,1,2 )                           
C      CALL  SYMBOL( WX-WC*10,WY2+0.1,WC,'AR PROCESS',0.0,10 )           
C                                                                       
C  ...  PLOT DATA AND TREND PLUS AR  ...                                
C                                                                       
C      CALL  PLOT( 0.0,-8.0,-3 )                                         
C      CALL  XYAXIS( 0.0,0.0,WX,WY,X0,X1,Y0,Y1,DX,DY,2,MESH )            
C      CALL  PLOTD( WORK,N,Y0,Y1,WX,WY,1,1,2 )                           
C      DO 90 I=1,N                                                       
C   90 Z(I) = E(1,1,I) + E(1,M1+1,I)                                     
C      CALL  PLOTD( Z,N,Y0,Y1,WX,WY,1,1,2 )                              
C      CALL  SYMBOL(WX-WC*26,WY+0.1,WC,'ORIGINAL AND TREND PLUS AR',0.0, 
C     *             28 )                                                 
C                                                                       
  100 IF( M4 .EQ. 0 )  GO TO 130                                        
C      IF( M2 .NE. 0 )  CALL  PLOTI                                      
C      IF( M2 .NE. 0 )  CALL  TITLEP( TITLE,A )                          
C      IF( M2 .NE. 0 )  CALL  PLOT( XX,YY,-3 )                           
C      IF( M2 .EQ. 0 )  CALL  PLOT( 16.0,WY1+WY2+3.0,-3 )                
C                                                                       
C  ...  PLOT TRADING DAY EFFECT  ...                                    
C                                                                       
C      CALL  XYAXIS( 0.0,0.0,WX,WY2,X0,X1,Y02,Y12,DX,DY2,2,MESH )        
      DO 110 I=1,N                                                      
cxx  110 COMP4(I) = E(2,1,I)                                                   
      COMP4(I) = E(2,1,I)
  110 CONTINUE
C      CALL  PLOTD( Z,N,Y02,Y12,WX,WY2,1,1,2 )                           
C      CALL  SYMBOL( WX-WC*18,WY2+0.1,WC,'TRADING DAY EFFECT',0.0,18 )   
C                                                                       
C  ...  PLOT SEASONAL PLUS TRADING DAY EFFECT  ...                      
C                                                                       
C      CALL  PLOT( 0.0,-WY2-2.0,-3 )                                     
C      CALL  XYAXIS( 0.0,0.0,WX,WY2,X0,X1,Y02,Y12,DX,DY2,2,MESH )        
C      DO 120 I=1,N                                                      
C  120 Z(I) = E(1,M12+1,I) + E(2,1,I)                                    
C      CALL  PLOTD( Z,N,Y02,Y12,WX,WY2,1,1,2 )                           
C      CALL  SYMBOL( WX-WC*32,WY2+0.1,WC,'SEASONAL PLUS TRADING DAY EFFEC
C     *T',0.0,32 )                                                       
  130 CONTINUE                                                          
C      CALL  PLOTE                                                       
C      CALL  PLTCE                                                       
      RETURN                                                            
cxx  77  format(' ^*&$^$%#&^&*nvkdnfvafv')
      E N D                                                             
cc      SUBROUTINE  REDATA( MT,DATA,ISW,IPR,X,N,TITLE,XM )                     
cxx      SUBROUTINE  REDATAD( DATA,ISW,IPR,X,N,XM )                     
      SUBROUTINE  REDATAD( DATA,ISW,X,N,XM )                     
C                                                                       
C     THIS SUBROUTINE IS USED FOR THE LOADING OF ORIGINAL DATA,         
C     THE DATA IS LOADED THROUGH THE DEVICE SPECIFIED BY @MT@.          
C     EACH DATA SET IS COMPOSED OF TITLE, DATA LENGTH, DATA FORMAT      
C     AND ORIGINAL DATA.  IF ISW IS SET EQUAL TO 0, THE MEAN VALUE      
C     IS SUBTRACTED FROM THE ORIGINAL DATA.                             
C                                                                       
C       INPUTS:                                                         
C         MT:     INPUT DEVICE SPECIFICATION                            
C         ISW:    =0 MEAN VALUE IS SUBTRACTED FROM THE ORIGINAL DATA    
C                 =1 MEAN VALUE IS NOT SUBTRACTED                       
C         IPR:    =1 TO PRINT OUT DATA BY D-FORMAT                      
C                 =2 TO PRINT OUT DATA BY F-FORMAT                      
C         TITLE:  TITLE OF DATA                                         
C         N:      DATA LENGTH                                           
C         FORM:   INPUT DATA FORMAT SPECIFICATION                       
C         X(I) (I=1,N):  ORIGINAL DATA                                  
C                                                                       
C       OUTPUTS:                                                        
C         X:      ORIGINAL DATA                                         
C         N:      DATA LENGTH                                           
C         TITLE:  TITLE OF DATA                                         
C         XM:     MEAN                                                  
C                                                                       
cxx      IMPLICIT REAL*8( A-H,O-Z )                                        
cc      REAL * 4     FORM(20), TITLE(20)
cx      REAL*8       X(1) ,DATA(1)         
cxx      REAL*8       X(N) ,DATA(N)
      INTEGER :: ISW, N         
      REAL(8) :: DATA(N), X(N), XM
      REAL(8) :: FN, S1, S2, S3, S4, XX
C                                                                       
C       LOADING OF TITLE, DATA LENGTH, FORMAT SPECIFICATION AND DATA    
C                                                                       
C      READ( MT,5 )      TITLE                                           
C      READ( MT,1 )      N                                               
C      READ( MT,5 )      FORM                                            
C      READ( MT,FORM )   (X(I),I=1,N)                                    
C                                                                       
C       ORIGINAL DATA PRINT OUT                                         
C                                                                       
      DO 53 I=1,N
cxx 53      X(I)=DATA(I)
         X(I)=DATA(I)
 53   CONTINUE

c      WRITE( 6,9 )     N , (FORM(I),I=1,17)                             
c      WRITE( 6,8 )     TITLE                                            
c      IF(IPR .EQ. 1)  WRITE( 6,3 )  (X(I),I=1,N)                        
c      IF(IPR .EQ. 2)  WRITE( 6,4 )  (X(I),I=1,N)                        
C                                                                       
C       SAMPLE MOMENTS COMPUTATION                                      
C   by Sato
cc	return



      FN = N                                                            
      S1 = 0.0D0                                                        
      DO 10     I=1,N                                                   
cxx   10 S1 = S1 + X(I)                                                    
      S1 = S1 + X(I)
   10 CONTINUE
      XM = S1 / FN                                                      
      S2 = 0.0D0                                                        
      S3 = 0.0D0                                                        
      S4 = 0.0D0                                                        
      DO 20  I=1,N                                                      
      XX = X(I) - XM                                                    
      S2 = S2 + XX**2                                                   
      S3 = S3 + XX**3                                                   
cxx   20 S4 = S4 + XX**4                                                   
      S4 = S4 + XX**4
   20 CONTINUE
C                                                                       
      S2 = S2 / FN                                                      
      S3 = S3 / (FN*S2*DSQRT(S2))                                       
      S4 = S4 / (FN*S2*S2)                                              
C                                                                       
C       PRINT OUT SAMPLE MOMENTS                                        
C                                                                       
c      WRITE( 6,7 )   XM, S2, S3, S4                                     
      IF( ISW .EQ. 1 )  RETURN                                          
C                                                                       
C       MEAN DELETION                                                   
C                                                                       
      DO 30  I=1,N                                                      
cxx   30 X(I) = X(I) - XM                                                  
      X(I) = X(I) - XM
   30 CONTINUE
C                                                                       
      RETURN                                                            
C                                                                       
cxx    1 FORMAT( 16I5 )                                                    
cxx    3 FORMAT( 1H ,10D13.5 )                                             
cxx    4 FORMAT( 1H ,10F10.4 )                                             
                                                                        
cxx    5 FORMAT ( 20A4 )                                                   
cxx    7 FORMAT( 1H0,'MEAN      =',D17.8,/,' VARIANCE  =',D17.8,/,         
cxx     *        ' SKEWNESS  =',D17.8,/,' KURTOSIS  =',D17.8 )             
cxx    8 FORMAT( 1H ,20A4 )                                                
cxx    9 FORMAT( 1H0,'<<  ORIGINAL DATA  X(I) (I=1,N)  >>',5X,'N =',I5,4X, 
cxx     *       'FORMAT =',17A4 )                                          
C                                                                       
      E N D                                                             
      SUBROUTINE  SETFGH                                                
C                                                                       
C  ...  SET F, G AND H MATRICES OF STATE SPACE MODEL  ...               
C                                                                       
cxx      IMPLICIT  REAL*8(A-H,O-Z)                                         
cxx      INTEGER   PERIOD, SORDER
      INTEGER :: PERIOD, SORDER
      REAL(8) :: F1, F2, F3, A1, A2, A3, DI, UI, TDF
cc      COMMON    /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,        
cc     *              NS, NI, MISING, IOUT, LL, N, NYEAR, nmonth,NPREDS,  
cc     *                    NPREDE, NPRED
ccx      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)
      COMMON    /COMSM2/  M1, M2, M3, M4, M5, M, L, PERIOD, SORDER,
     *                    NYEAR, nmonth
      COMMON  /COMSM3/  F1(10), F2(10), F3(300), A1(10), A2(10), A3(300)
     *                   ,DI, UI(3), TDF(7)
C     *                   ,DI, TAU(3)                                   
C                                                                       
      DI = 1.0D0                                                        
C                                                                       
      IF( M1 .EQ. 0 )  GO TO 40                                         
cxx      GO TO (10,20,30), M1                                              
      IF( M1 .EQ. 1 ) GO TO 10
      IF( M1 .EQ. 2 ) GO TO 20
      IF( M1 .EQ. 3 ) GO TO 30                                               
C                                                                       
   10 F1(1) =  1.0D0                                                    
      A1(1) =  1.0D0                                                    
      GO TO 40                                                          
   20 F1(1) = -1.0D0                                                    
      F1(2) =  2.0D0                                                    
      A1(1) =  2.0D0                                                    
      A1(2) = -1.0D0                                                    
      GO TO 40                                                          
   30 F1(1) =  1.0D0                                                    
      F1(2) = -3.0D0                                                    
      F1(3) =  3.0D0                                                    
      A1(1) =  3.0D0                                                    
      A1(2) = -3.0D0                                                    
      A1(3) =  1.0D0                                                    
C                                                                       
cc   40 IF( KSEA .EQ. 0 )  GO TO 160                                      
cc      IF( KSEA .EQ. 2 )  GO TO 130                                      
cc      IF( KSEA .EQ. -1 ) GO TO 80                                       
   40 IF( SORDER .EQ. 0 )  GO TO 160                                    
      IF( SORDER .EQ. 2 )  GO TO 130                                    
      IF( SORDER .EQ. -1 ) GO TO 80                                     
      DO 70  I=1,M3                                                     
      A3(I) = -1.0D0                                                    
cxx   70 F3(I) = -1.0D0
      F3(I) = -1.0D0
   70 CONTINUE
      GO TO 150                                                         
cxx   80 DO 90 I=1,M3                                                      
cxx      A3(I) = 0.0D0                                                     
cxx   90 F3(I) = 0.0D0
   80 CONTINUE                                                      
      A3(1:M3) = 0.0D0                                                     
      F3(1:M3) = 0.0D0                                                     
      A3(M3)= 1.0D0                                                     
      F3(1) = 1.0D0                                                     
      GO TO 150                                                         
cc  130 DO 140  I=1,ISEA-1                                                
  130 DO 140  I=1,PERIOD-1
      F3(I) = -I                                                        
      A3(I) = -I - 1                                                    
cc      F3(ISEA+I-1) = I - ISEA - 1                                       
cc  140 A3(ISEA+I-1) = I - ISEA                                           
      F3(PERIOD+I-1) = I - PERIOD - 1                                   
cxx  140 A3(PERIOD+I-1) = I - PERIOD
      A3(PERIOD+I-1) = I - PERIOD
  140 CONTINUE
  150 CONTINUE                                                          
  160 CONTINUE                                                          
C                                                                       
      RETURN                                                            
C                                                                       
C                                                                       
      E N D                                                             
      SUBROUTINE  SMOTH3( Z,R,S,T,TRADE,IMIS,N,LM1,F,ILKF,ISMT )
cc      SUBROUTINE  SMOTH3( Z,R,S,T,TRADE,REG,EPRED,VPRED,IMIS,N,LM1,
cc     *                    F,ILKF,ISMT )
CC     *                    MJ,F,ILKF,ISMT )                              
C                                                                       
C     ...  INFORMATION SQUARE ROOT FILTER & SMOOTHER  ...               
C          FOR SEASONAL ADJUSTMENT                                      
C                                                                       
cxx      IMPLICIT  REAL*8  ( A-H,O-Z )                                     
cc      REAL*8     Z(N), T(L,LM1,N), REG(N,M5), EPRED(N), VPRED(N), TRADE
cc      DIMENSION  S(MJ,MJ), R(MJ,MJ), IMIS(N), TRADE(N,7)                
cc      DIMENSION  D(100), WI(10,10), X(40), TT(40), E(40)
cxx      REAL*8     Z(N), T(L,LM1,N), REG(N,M5)
cxx      DIMENSION  S(M+1,M+1), R(LM1+1,LM1), IMIS(N), TRADE(N,7)            
cxx      DIMENSION  D(LM1+1), WI(3,3), X(LM1), TT(L), E(M)                
cxx      INTEGER    PERIOD, SORDER
      INTEGER :: IMIS(N), N, LM1, ILKF, ISMT
      REAL(8) :: Z(N), R(LM1+1,LM1), S(M+1,M+1), T(L,LM1,N),
     1           TRADE(N,7), F
      INTEGER :: PERIOD, SORDER
      REAL(8) :: REG(N,M5), D(LM1+1), WI(3,3), X(LM1), TT(L), E(M),
     1           F1, F2, F3, A1, A2, A3, DI, UI, TDF,
     2           DJACOB, FC, SIG2, AIC, FI, SIG2I, AICI, GI, GC,
     3           PAI, SDET, SUM, tmpp, tmp
cc      COMMON    /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,        
cc     *             NS, NI, MISING, IOUT, LL, NN, NYEAR,nmonth, NPREDS, 
cc     *                    NPREDE, IPRED
ccx      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)
c      COMMON    /CMFUNC/  DJACOB, FF, SIG2, AIC                        
cc      COMMON    /CMFUNC/  DJACOB,FF,
cc     *            SIG2,AIC,FI,SIG2I,AICI,GI(20),G(20) 
ccx      COMMON    /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
      COMMON    /COMSM2/  M1, M2, M3, M4, M5, M, L, PERIOD, SORDER,
     *                    NYEAR, nmonth
      COMMON  /COMSM3/  F1(10), F2(10), F3(300), A1(10), A2(10), A3(300)
     *                   ,DI, UI(3), TDF(7)                            
      COMMON  /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(200),GC(200)

      DATA  PAI/3.1415926535D0/                                         
C                                                                       
cc      M    = M1 + M2 + M3 + M4 + M5                                     
      M12  = M1 + M2                                                    
      M123 = M1 + M2 + M3                                               
      M45 = M4 + M5                                                     
      LM   = L + M                                                      
      L1   = L + M1                                                     
      L12  = L + M1 + M2                                                
      L123 = L + M1 + M2 + M3                                           
      L1234= L + M1 + M2 + M3 + M4                                      
      LL2  = ID(M1) + ID(M2)                                            
      SIG2 = 0.0D0                                                      
      SDET = 0.0D0
      NNN = 0                                                      
C
      DO 500  II=1,N
C                                                                       
C     ...  TIME  UPDATE (PREDICTION)  ...                               
C                                                                       
      IF( M1 .EQ. 0 )  GO TO 140                                        
      IF(M1.EQ.1)  GO TO 115                                            
cxx      DO 110  J=2,M1                                                    
      DO 111  J=2,M1
      DO 110  I=1,J                                                     
cxx  110 R(L+I,L+J-1) = S(I,J)                                             
      R(L+I,L+J-1) = S(I,J)
  110 CONTINUE
  111 CONTINUE
  115 DO 130  I=1,M1                                                    
      SUM = 0.0D0                                                       
      DO 120  J=I,M1                                                    
cxx  120 SUM = SUM + S(I,J)*F1(J)                                          
      SUM = SUM + S(I,J)*F1(J)
  120 CONTINUE
cxx  130 R(L+I,L1) = SUM                                                   
      R(L+I,L1) = SUM
  130 CONTINUE
C                                                                       
  140 IF( M2 .EQ. 0 )  GO TO 5                                          
      IF(M2 .EQ. 1)  GO TO 155                                          
cxx      DO 150  J=2,M2                                                    
      DO 151  J=2,M2
      JJ = M1 + J                                                       
      DO 150  I=1,JJ                                                    
cxx  150 R(L+I,L1+J-1) = S(I,M1+J)                                         
      R(L+I,L1+J-1) = S(I,M1+J)
  150 CONTINUE
  151 CONTINUE
  155 DO 170  I=1,M2                                                    
      SUM = 0.0D0                                                       
      DO 160  J=I,M2                                                    
cxx  160 SUM = SUM + S(M1+I,M1+J)*F2(J)                                    
      SUM = SUM + S(M1+I,M1+J)*F2(J)
  160 CONTINUE
cxx  170 R(L1+I,L12) = SUM                                                 
      R(L1+I,L12) = SUM
  170 CONTINUE
      DO 190  I=1,M1                                                    
      SUM = 0.0D0                                                       
      DO 180  J=1,M2                                                    
cxx  180 SUM = SUM + S(I,M1+J)*F2(J)                                       
      SUM = SUM + S(I,M1+J)*F2(J)
  180 CONTINUE
cxx  190 R(L+I,L12) = SUM                                                  
      R(L+I,L12) = SUM
  190 CONTINUE
cxx  195 CONTINUE                                                          
C                                                                       
    5 IF( M3 .EQ. 0 )  GO TO 55                                         
cxx      DO 10  J=2,M3                                                     
      DO 11  J=2,M3
      JJ = M12 + J                                                      
      DO 10  I=1,JJ                                                     
cxx   10 R(L+I,L12+J-1) = S(I,M12+J)                                       
      R(L+I,L12+J-1) = S(I,M12+J)
   10 CONTINUE
   11 CONTINUE
      DO 30  I=1,M3                                                     
      SUM = 0.0D0                                                       
      DO 20  J=I,M3                                                     
cxx   20 SUM = SUM + S(M12+I,M12+J)*F3(J)                                  
      SUM = SUM + S(M12+I,M12+J)*F3(J)
   20 CONTINUE
cxx   30 R(L12+I,L123) = SUM                                               
      R(L12+I,L123) = SUM
   30 CONTINUE
C
      DO 50  I=1,M12                                                    
      SUM = 0.0D0                                                       
      DO 40  J=1,M3                                                     
cxx   40 SUM = SUM + S(I,M12+J)*F3(J)                                      
      SUM = SUM + S(I,M12+J)*F3(J)
   40 CONTINUE
cxx   50 R(L+I,L123) = SUM                                                 
      R(L+I,L123) = SUM
   50 CONTINUE
C                                                                       
   55 IF( M45 .EQ. 0 )  GO TO 70                                        
cxx      DO 60  J=1,M45
      DO 61  J=1,M45
      JJ = M123 + J                                                     
      DO 60  I=1,JJ                                                     
cxx   60 R(L+I,L123+J) = S(I,M123+J)
      R(L+I,L123+J) = S(I,M123+J)
   60 CONTINUE
   61 CONTINUE
   70 CONTINUE                                                          
C                                                                       
      IF( M1 .EQ. 0 )  GO TO 215                                        
      DO 210  I=1,M1                                                    
cxx  210 R(L+I,1) = -R(L+I,L+1)                                            
      R(L+I,1) = -R(L+I,L+1)
  210 CONTINUE
  215 IF( M2 .EQ. 0 ) GO TO 225                                         
      DO 220  I=1,M12                                                   
cxx  220 R(L+I,LL2) = -R(L+I,L+M1+1)
      R(L+I,LL2) = -R(L+I,L+M1+1)
  220 CONTINUE
C                                                                       
  225 IF( M3 .EQ. 0 )  GO TO 235                                        
      DO 230  I=1,M123                                                  
cxx  230 R(L+I,L) = -R(L+I,L+M12+1)                                        
      R(L+I,L) = -R(L+I,L+M12+1)
  230 CONTINUE
C                                                                       
cxx  235 DO 250  I=1,L                                                     
cxx      DO 250  J=1,LM+1                                                  
cxx  250 R(I,J) = 0.0D0
  235 CONTINUE
      R(1:L,1:LM+1) = 0.0D0
      DO 240  I=1,M                                                     
cxx  240 R(L+I,L+M+1) = S(I,M+1)                                           
      R(L+I,L+M+1) = S(I,M+1)
  240 CONTINUE
C                                                                       
      DO 260 I=1,L                                                      
cxx  260 R(I,I) = UI(I)                                                    
      R(I,I) = UI(I)
  260 CONTINUE
C                                                                       
cc      CALL  HUSHL7( R,D,MJ,L+M+1,L12+1,L+M123 )                         
      CALL  HUSHL7( R,D,LM1+1,LM1,L12+1,L+M123 )
C                                                                       
cxx      DO 300  I=1,L                                                     
      DO 301  I=1,L
      DO 300  J=1,LM+1                                                  
cxx  300 T(I,J,II) = R(I,J)
      T(I,J,II) = R(I,J)
  300 CONTINUE
  301 CONTINUE
C                                                                       
C     ...  MEASUREMENT  UPDATE (FILTERING)  ...                         
C                                                                       
c    modified at 96.10 by S.S.
c 400  continue
cxx 400    IF( IMIS(II) .EQ. 1 )  GO TO 420
        IF( IMIS(II) .EQ. 1 )  GO TO 420
c
      nnn=nnn+1
C     WRITE(6,997) M,UI(1)                                              
cxx      DO 410  J=1,M+1                                                   
      DO 411  J=1,M+1
      S(M+1,J) = 0.0D0                                                  
      DO 410  I=1,J                                                     
cxx  410 S(I,J) = R(L+I,L+J)                                               
      S(I,J) = R(L+I,L+J)
  410 CONTINUE
  411 CONTINUE
C                                                                       
      S(M+1,1)    = DI                                                  
      S(M+1,M1+1) = DI                                                  
      S(M+1,M12+1)= DI                                                  
      S(M+1,M+1)  = DI*Z(II)                                            
      IF( M4 .EQ. 6 )  then
      DO 414  I=1,6                                                     
      JJ = II                                                           
      IF( ILKF .EQ. 0 )  JJ = N - II + 1                                
cxx  414 S(M+1,M123+I) = DI*(TRADE(JJ,I) - TRADE(JJ,7))                    
      S(M+1,M123+I) = DI*(TRADE(JJ,I) - TRADE(JJ,7))
  414 CONTINUE
      end if

      if(m4 .eq. 1) then

      JJ = II                                                           
      IF( ILKF .EQ. 0 )  JJ = N - II + 1                                
      tmpp=trade(jj,2)+trade(jj,3)
      tmpp=tmpp+trade(jj,4)+trade(jj,5)+trade(jj,6)
      s(m+1,m123+1)=di*(trade(jj,1)+trade(jj,7)-0.4*tmpp)
      end if

cxx  415 IF( M5 .EQ. 0 )  GO TO 418                                        
      IF( M5 .EQ. 0 )  GO TO 418
      DO 417  I=1,M5
cxx  417 S(M+1,M123+M4+I) = REG(I,II)                                      
      S(M+1,M123+M4+I) = REG(I,II)
  417 CONTINUE
  418 CONTINUE                                                          
C                                                                       
C     CALL  PRINT2( S,M+1,M+1,MJ,1 )                                    
cc      CALL  HUSHL4( S,MJ,M+1,M+1,1,0 )                                  
      CALL  HUSHL4( S,M+1,M+1,M+1,1,0 )
C     CALL  PRINT2( S,M+1,M+1,MJ,1 )                                    
C                                                                       
C     ...  PREDICTION ERROR & DETERMINANT  ...                          
C                                                                       
cc      IF( IPRED.EQ.1 .AND. II.GE.NPREDS )  GO TO 430                    
      IF(ILKF .EQ. 0)  GO TO 500                                        
      SIG2 = SIG2 + S(M+1,M+1)**2                                       
      DO 440  I=1,M                                                     
cxx  440 SDET = SDET + DLOG( S(I,I)**2 ) - DLOG( R(L+I,L+I)**2 )           
      SDET = SDET + DLOG( S(I,I)**2 ) - DLOG( R(L+I,L+I)**2 )
  440 CONTINUE
      GO TO 500                                                         
C                                                                       
C  ...  LONG RANGE PREDICTION  ...                                      
C                                                                       
cc  430 CALL  RECOEF( R,LM,LM,MJ,X )                                      
cxx  430 CALL  RECOEF( R,LM,LM,LM1+1,X )
      CALL  RECOEF( R,LM,LM,LM1+1,X )
      SUM = ID(M1)*X(L+1) + ID(M2)*X(L1+1) + ID(M3)*X(L12+1)            
      IF( M4 .EQ. 6 )  then
      DO 431 J=1,6                                                      
cxx  431 SUM = SUM + (TRADE(II,J)-TRADE(II,7))*X(L123+J)
      SUM = SUM + (TRADE(II,J)-TRADE(II,7))*X(L123+J)
  431 CONTINUE
      end if
      if(m4 .eq. 1) then
      tmp=trade(II,2)+trade(II,3)+trade(II,4)+trade(II,5)+trade(II,6)
      sum=sum+(trade(II,1)+trade(II,7)-0.4*tmp)*x(L123+1)
      end if
cc  432 IF( M5 .EQ. 0 )  GO TO 434                                        
cc      DO 433 J=1,M5                                                     
cc  433 SUM = SUM + REG(II,J)*X(L1234+J)                                  
cc  434 EPRED(II) = SUM                                                   
cc      SUM = 1.0D0                                                       
cc      DO 445 I=1,M                                                      
cc  445 SUM = SUM * S(I,I)/R(L+I,L+I)                                     
cc      VPRED(II) = SUM**2                                                
ccc      WRITE(6,977)  II, EPRED(II), VPRED(II)                            
cc  977 FORMAT( 1H ,I5,2F15.7 )                                           
C                                                                       
C  ...  MISSING OBSERVATION  ...                                        
C                                                                       
cxx  420 DO 425  I=1,M
  420 DO 426  I=1,M
      DO 425  J=1,M+1                                                   
cxx  425 S(I,J) = R(L+I,L+J)
      S(I,J) = R(L+I,L+J)
  425 CONTINUE
  426 CONTINUE
  500 CONTINUE                                                          
C                                                                       
C     ... LIKELIHOOD COMPUTATION  ...                                   
C                                                                       

C   99/8/12
      IF(ILKF .EQ. 0)  GO TO 600                                        

      SIG2 = SIG2/(2.0D0*nnn)                                           
      F = -0.5D0*(nnn*( DLOG(PAI*2.D0) + DLOG(SIG2) + 1.0D0 ) + SDET)   
cc      FF = F                                                            
      FC = F                                                            
c    modified 97/10/17
      AIC = -2.0D0*(F+DJACOB) + 2.0D0*(M2 + L + 1 + m)                  
cxx  901 FORMAT(1H ,4D16.9)                                                
C 900 FORMAT(1H ,'F =',D15.8,5X,'SIG2 =',D15.8,5X,'DET =',D15.8,        
C    1     5X,'FF =',D15.8,5X,'SD =',D15.8,5X,'LOG V =',D15.8 )         
cxx  900 FORMAT(1H ,3D20.10)                                               
C                                                                       
C                                                                       
C                                                                       
C     ...  SMOOTHING  ...                                               
C                                                                       
  600 IF(ISMT .EQ.  0) RETURN                                           
cc      CALL  RECOEF( S,M,M,MJ,X )
      CALL  RECOEF( S,M,M,M+1,X )                                      
C                                                                       
      DO 610 I=1,M                                                      
cxx  610 E(I) = X(I)
      E(I) = X(I)
  610 CONTINUE
C                                                                       
      DO 1200  III=1,N-1                                                
C                                                                       
      II = N-III+1                                                      
C                                                                       
C  ...  WI = INVERSE OF T(I,J)  ...                                     
C                                                                       
      WI(1,1) = 1.0D0/T(1,1,II)                                         
      IF( L .EQ. 1 )  GO TO 615                                         
      WI(2,2) = 1.0D0/T(2,2,II)                                         
      WI(1,2) = -T(1,2,II)*WI(2,2)*WI(1,1)                              
      IF( L .EQ. 2 )  GO TO 615                                         
      WI(3,3) = 1.0D0/T(3,3,II)                                         
      WI(2,3) = -T(2,3,II)*WI(2,2)*WI(3,3)                              
      WI(1,3) = -(T(1,2,II)*WI(2,3) + T(1,3,II)*WI(3,3))*WI(1,1)        
  615 DO 630  I=1,L                                                     
      SUM = T(I,LM1,II)                                                 
      DO 620  J=1,M                                                     
cxx  620 SUM = SUM - T(I,L+J,II)*X(J)
      SUM = SUM - T(I,L+J,II)*X(J)
  620 CONTINUE
cxx  630 TT(I) = SUM
      TT(I) = SUM
  630 CONTINUE
C                                                                       
      DO 650  I=1,L                                                     
      SUM = 0.0D0                                                       
      DO 640  J=I,L                                                     
cxx  640 SUM = SUM + WI(I,J)*TT(J)                                         
      SUM = SUM + WI(I,J)*TT(J)
  640 CONTINUE
cxx  650 D(I) = SUM
      D(I) = SUM
  650 CONTINUE
      X(1) = X(1) - D(1)                                                
      IF(M1.GT.0 .AND. M2.GT.0)  X(M1+1) = X(M1+1) - D(2)               
      X(M12+1) = X(M12+1) - D(L)                                        
      DO 660  I=1,M                                                     
cxx  660 D(I) = X(I)
      D(I) = X(I)
  660 CONTINUE
C                                                                       
      IF( M1 .EQ. 0 )  GO TO 715                                        
      DO 710  I=1,M1                                                    
cxx  710 X(I) = F1(I)*D(M1)                                                
      X(I) = F1(I)*D(M1)
  710 CONTINUE
  715 IF(M1.LE.1)  GO TO  730                                           
      DO 720  I=2,M1                                                    
cxx  720 X(I) = X(I) + D(I-1)                                              
      X(I) = X(I) + D(I-1)
  720 CONTINUE
C                                                                       
  730 IF( M2 .EQ. 0 )  GO TO 760                                        
      DO 740  I=1,M2                                                    
cxx  740 X(M1+I) = F2(I)*D(M12)
      X(M1+I) = F2(I)*D(M12)
  740 CONTINUE
      IF(M2 .LE. 1)  GO TO 760                                          
      DO 750  I=2,M2                                                    
cxx  750 X(M1+I) = X(M1+I) + D(M1+I-1)                                     
      X(M1+I) = X(M1+I) + D(M1+I-1)
  750 CONTINUE
C                                                                       
  760 IF( M3 .LE. 0 )  GO TO 775                                        
      DO 770  I=1,M3                                                    
cxx  770 X(M12+I) = F3(I)*D(M123)                                          
      X(M12+I) = F3(I)*D(M123)
  770 CONTINUE
  775 IF( M3 .LE. 1 )  GO TO 785                                        
      DO 780  I=2,M3                                                    
cxx  780 X(M12+I) = X(M12+I) + D(M12+I-1)
      X(M12+I) = X(M12+I) + D(M12+I-1)
  780 CONTINUE
C                                                                       
  785 CONTINUE                                                          
cxx 2000 CONTINUE                                                          
      DO 790 I=1,M                                                      
      T(1,I,II) = E(I)                                                  
cxx  790 E(I) = X(I)                                                       
      E(I) = X(I)
  790 CONTINUE
C                                                                       
 1200 CONTINUE                                                          
C                                                                       
      DO 1210 I=1,M                                                     
cxx 1210 T(1,I,1) = E(I)                                                   
      T(1,I,1) = E(I)
 1210 CONTINUE
      IF( M4 .EQ. 6 ) then
      SUM = 0.0D0                                                       
      DO 1220 I=1,6                                                     
      TDF(I) = X(M123+I)                                                
cxx 1220 SUM = SUM + TDF(I)                                                
      SUM = SUM + TDF(I)
 1220 CONTINUE
      TDF(7) = -SUM                                                     
      end if
      IF( M4 .EQ. 1 ) then
      TDF(1) = X(M123+1)
      TDF(7) = X(M123+1)
      do 1222 i=2,6
cxx 1222    TDF(i) = -0.4*TDF(1)
         TDF(i) = -0.4*TDF(1)
 1222 continue
      end if
C                                                                       
      RETURN                                                            
      E N D                                                             
cxx      SUBROUTINE  SPARAM0( N,IPAR,NIP,para,NPA )
      SUBROUTINE  SPARAM0( IPAR,NIP )
C                                                                       
C  ...  Set or read control parameters ...                              
C                                                                       
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  IPAR(NIP), para(NPA)
cxx      INTEGER    PERIOD, SORDER, TRADE
      INTEGER :: IPAR(NIP), NIP
      INTEGER :: PERIOD, SORDER, TRADE
      COMMON     /COMSM2/  M1, M2, M3, M4, M5, M, L, PERIOD, SORDER,    
     *                     NYEAR, nmonth
cxx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
      COMMON    /CCC/     ISW, ISMT, IDIF, LOG, MESH              
C                                                                       
C  ...  SET DEFAULT VALUES  ...                                         
C                                                                       
      M1 = IPAR(1)       
      M2 = IPAR(2)
      m4 = 0
      M5 = 0                        
      LOG  = IPAR(5)                            
cxx      IPR  = 7                                            
      ISW  = 1                                                          
      IDIF = IPAR(7)
      if(idif .gt. 2) idif=1
      if(idif .lt. 1) idif=1
      MESH = 1                                                          
      PERIOD = IPAR(3)                               
      SORDER = IPAR(4)                                                 
      TRADE  = IPAR(6)                
      if(TRADE .EQ. 1) TRADE = 7
      IF(TRADE .GE. 1)  NYEAR=IPAR(8)
      IF(TRADE .GE. 1)  nmonth=IPAR(9)
      M3 = (PERIOD - 1)*SORDER                                          
      IF(SORDER.EQ.-1)  M3 = PERIOD                                     
      IF(TRADE .GE. 1)  M4 = TRADE-1                                   
      M = M1 + M2 + M3 + M4 + M5                                        
      L = ID(M1) + ID(M2) + ID(M3)                                      
      L = MAX(2,L)
C                                                                       
      RETURN                                                            
      END                                                             
cc      SUBROUTINE  SPARAM( MT,A,IPAR,para,iopt )                                        
cxx      SUBROUTINE  SPARAM( N,A,NA,IPAR,NIP,para,NPA,iopt )
      SUBROUTINE  SPARAM( A,NA,para,NPA,iopt )
C                                                                       
C  ...  Set or read control parameters ...                              
C                                                                       
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      REAL*8     TAU2(3), A(40), PAC(10), F1,F2,F3,A1,AR,A3,para(26)
cc      REAL*8     ARCC(10)
cc      INTEGER   PERIOD, SORDER, BSPAN, OUTLIR, TRADE, YEAR, PRED        
cc      INTEGER   PREDS, PREDE ,IPAR(11),month
cc      DIMENSION  TAU2(3), PAC(10), ARCC(10)
cxx      DIMENSION  A(NA), IPAR(NIP), para(NPA)
cxx      DIMENSION  TAU2(3), PAC(M2), ARCC(M2)
cxx      INTEGER    PERIOD, SORDER
      INTEGER :: NA, NPA, iopt
      REAL(8) :: A(NA), para(NPA)
      INTEGER :: PERIOD, SORDER
      REAL(8) :: TAU2(3), PAC(M2), ARCC(M2),
     1           F1, F2, F3, A1, A2, A3, DI, UI, TDF
      COMMON     /COMSM2/  M1, M2, M3, M4, M5, M, L, PERIOD, SORDER,    
     *                     NYEAR, nmonth
cc     *              BSPAN, ISPAN, MISING, OUTLIR, LL, N, YEAR, month,
cc     *                    PREDS, PREDE, PRED
cc      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), AR(10), A3(30)    
ccx      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)    
cxx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
      COMMON  /COMSM3/  F1(10), F2(10), F3(300), A1(10), A2(10), A3(300)    
     *                   ,DI, UI(3), TDF(7)
      COMMON    /CCC/     ISW, ISMT, IDIF, LOG, MESH              
C      NAMELIST  /PARAM/  M1, M2, M5, PERIOD, SORDER, TRADE, MT, BSPAN,  
C     *                   ISPAN, MISING, OUTLIR, TAU2, PAC, IPR, IDIF,   
C     *                   LOG, YEAR, PRED, MESH                          
C                                                                       
C  ...  SET DEFAULT VALUES  ...                                         
C                                                                       
cc      M1 = IPAR(1)       
cc      M2 = IPAR(2)
cc      m4 = 0
cc      M5 = 0                        
cc      MT = 1                                                            
cc      LOG  = IPAR(5)                            
cc      IPR  = 7                                            
cc      ISW  = 1                                                          
cc      IDIF = IPAR(7)
cc      if(idif .gt. 2) idif=1
cc      if(idif .lt. 1) idif=1
cc      MESH = 1                                                          
cc      PERIOD = IPAR(3)                               
cc      SORDER = IPAR(4)                                                 
cc	TRADE  = IPAR(6)                
cc      if(TRADE .EQ. 1) TRADE = 7
cc      IF(TRADE .GE. 1)  YEAR=IPAR(10)
cc      IF(TRADE .GE. 1)  month=IPAR(11)
c      BSPAN  = IPAR(8)                                                 
c      ISPAN  = IPAR(9)          
c      if(n .lt. ispan) ispan = n
cc      bspan = n
cc      ispan = n
cc      if(iopt .lt. 0) ispan = ipar(9)
cc      PRED = 0                                                          
cc      MISING = 0                                                        
cc      OUTLIR = 0                                                        
      TAU2(1) = 0.005D0                                                 
      TAU2(2) = 0.800D0                                                 

C 99/8/12    
      TAU2(3) = 1.0D-3                                                  
C      TAU2(3) = 1.0D-5    
      if(m2 .eq. 0) tau2(2) = 0.001D0

cxx      do 18 i=1,7
cxx 18      TDF(i) = 0.0D00
      TDF(1:7) = 0.0D00
cc      DO 20  I=1,10                                                     
      DO 20  I=1,M2
cxx   20 PAC(I) = 0.88D0*(-0.6D0)**(I-1)                                   
      PAC(I) = 0.88D0*(-0.6D0)**(I-1)
   20 CONTINUE

      if(iopt .lt. 0) then
         do 21 i=1,3
          tau2(i)=para(i) - 0.0001
          if(tau2(i) .ge. 1.0D0) tau2(i)=1.0D0 - 0.1D-20
          if(tau2(i) .le. 0.0D0) tau2(i)=0.0D0 + 0.1D-20
 21       continue

c            write(*,*) tau2(1),tau2(2),tau2(3)
         if(m2 .gt. 0) then
            do 22 i=1,m2
cxx 22            arcc(i) = para(3+i)
               arcc(i) = para(3+i)
 22         continue
            call PARCOR( ARCC,m2,PAC )
c            write(*,*) PAC(1),PAC(2)
         end if
       end if
C                                                                       
C  ...  READ IN CONTROL PARAMETERS  ...                                 
C                                                                       
C      READ(5,PARAM)                                                     
C                                                                       
cc      M3 = (PERIOD - 1)*SORDER                                          
cc      IF(SORDER.EQ.-1)  M3 = PERIOD                                     
cc      IF(TRADE .GE. 1)  M4 = TRADE-1                                          
cc      M = M1 + M2 + M3 + M4 + M5                                        
cc      L = ID(M1) + ID(M2) + ID(M3)                                      
cc      LL= ID(M1) + ID(M2) + ID(M3) + ID(M4) + ID(M5)                    
      DO 30 I=1,L                                                       
cxx   30 A(I) = DASIN( TAU2(I)*2.0D0 - 1.0D0 )                             
      A(I) = DASIN( TAU2(I)*2.0D0 - 1.0D0 )
   30 CONTINUE
c      WRITE(6,610)                                                      
c      WRITE(6,600) M1,M2,M3,M4,M5,M,L,PERIOD,SORDER,BSPAN,ISPAN,MISING, 
c     *           OUTLIR, LL                                             
      IF(M2.EQ.0)  RETURN                                               
C                                                                       
      DO 40 I=1,M2                                                      
cxx   40 A(L+I) = DASIN( PAC(I)/0.90D0 )
      A(L+I) = DASIN( PAC(I)/0.90D0 )
   40 CONTINUE
C                                                                       
      RETURN                                                            
cxx  600 FORMAT( 10X,'M1     =',I3,/,10X,'M2     =',I3,/,10X,'M3     =',   
cxx     *   I3,/,10X,'M4     =',I3,/,10X,'M5     =',I3,/,10X,'M      =',   
cxx     *   I3,/,10X,'L      =',I3,/,10X,'PERIOD =',I3,/,10X,'SORDER =',   
cxx     *   I3,/,10X,'BSPAN  =',I3,/,10X,'ISPAN  =',I3,/,10X,'MISING =',   
cxx     *   I3,/,10X,'OUTLIR =',I3,/,10X,'LL     =',I3 )                   
cxx  610 FORMAT( 5X,'---  PROGRAM  DECOMP  ---' )                          
      E N D                                                             
      SUBROUTINE  STATE( X,A,K )                                        
C                                                                       
C  ...  TRANSFORMATION OF STATE VECTOR FOR TIME REVERSED MODEL  ...     
C                                                                       
cxx      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      DIMENSION  X(K), A(K), Y(30)                                      
cxx      DIMENSION  X(K), A(K), Y(K)
      INTEGER :: K
      REAL(8) :: X(K), A(K), Y(K), SUM
C                                                                       
      IF( K .EQ. 0 )  RETURN                                            
cxx      DO 20  I=1,K
      DO 21  I=1,K
      SUM = A(I)*X(1)                                                   
      IF( I .LT. K )  SUM = SUM + X(I+1)                                
      IF( I .EQ. 1 )  GO TO 20                                          
      DO 10  J=1,I-1                                                    
cxx   10 SUM = SUM + A(J)*Y(I-J)
      SUM = SUM + A(J)*Y(I-J)
   10 CONTINUE
   20 Y(I) = SUM                                                        
   21 CONTINUE
C                                                                       
      X(1) = Y(1)                                                       
      IF( K .EQ. 1 ) RETURN                                             
      DO 40  I=2,K                                                      
      SUM = 0.0D0                                                       
      DO 30  J=I,K                                                      
cxx   30 SUM = SUM + A(J)*Y(J-I+2)                                         
      SUM = SUM + A(J)*Y(J-I+2)
   30 CONTINUE
cxx   40 X(I) = SUM
      X(I) = SUM
   40 CONTINUE
C                                                                       
      RETURN                                                            
      E N D                                                             

      SUBROUTINE  TRADE( JSYEAR,nmonth,N,TDAY )

cxx      IMPLICIT REAL*8(A-H,O-Z)                                          

C
C  ...  This subroutine computes the number of days of the week
C       in each month, Nov.1981
C    modified at '96 by S.S.
C    This subroutine should not be used after 2099.
C
cxx      DIMENSION  TDAY(N,7), IX(12)
      INTEGER :: JSYEAR, nmonth, N
      REAL(8) :: TDAY(N,7)
      INTEGER :: IX(12)
      DATA   IX  /3,0,3,2,3,2,3,3,2,3,2,3/
C
c      open(1,file='tmp.dat')
      JS = JSYEAR - 1900
cc      I0 = MOD( JS+(JS-1)/4,7 ) + 1
      I2 = MOD( JS+(JS-1)/4,7 ) + 1
      JJ = 2-nmonth
      II = 0
 5     II = II + 1
      I1 = II + JS - 1
      IX(2) = 0
      IF( MOD(I1,4) .EQ. 0 )  IX(2) = 1
      IF( MOD(I1+1900,100).EQ.0 )  IX(2) = 0
      IF( MOD(I1+1900,400).EQ.0 )  IX(2) = 1
      DO 30  J=1,12
      DO 10  I=1,7
cxx 10    if( jj.gt.0 ) TDAY(JJ,I) = 4.0
       if( jj.gt.0 ) TDAY(JJ,I) = 4.0
 10   CONTINUE
C
      IE = IX(J)
      IF( IE .EQ. 0 )  GO TO 28
      I0 = I2
      DO 20  I=1,IE
      I2 = I0 + I
      IF( I2 .GT. 7 ) I2 = I2 - 7
cxx 20    if( jj.gt.0 ) TDAY(JJ,I2) = 5.0
       if( jj.gt.0 ) TDAY(JJ,I2) = 5.0
 20   CONTINUE
cc      I0 = I2

c       if(jj .gt.0) WRITE(1,*)  (TDAY(JJ,I),I=1,7)
 28   JJ = JJ + 1
      IF( JJ .GT.N ) then 

c      do 29 i9=1,n
c 29      write(1,*) (tDAY(i9,j9),j9=1,7)

         RETURN
      end if
 30   continue
      GO TO 5
C
      E N D

      SUBROUTINE  TRADE2( JSYEAR,nquart,N,TDAY )

cxx      IMPLICIT REAL*8(A-H,O-Z)                                          

C
C  ...  This subroutine computes the number of days of the week
C       in each quarter.
C    modified at '96 by S.S.
C    This subroutine should not be used after 2099.
C
cxx      DIMENSION  TDAY(N,7), IX(4)
      INTEGER :: JSYEAR, nquart, N
      REAL(8) :: TDAY(N,7)
      INTEGER :: IX(4)
      DATA   IX  /6,7,8,8/
C
c      write(6,*) nquart, jsyear
      JS = JSYEAR - 1900
cc      I0 = MOD( JS+(JS-1)/4,7 ) + 1
      I2 = MOD( JS+(JS-1)/4,7 ) + 1
      JJ = 2-nquart
      II = 0
 5     II = II + 1
      I1 = II + JS - 1
      IX(1) = 6
      IF( MOD(I1,4) .EQ. 0 )  IX(1) = 7
      IF( MOD(I1+1900,100).EQ.0 )  IX(1) = 6
      IF( MOD(I1+1900,400).EQ.0 )  IX(1) = 7
      DO 30  J=1,4
      DO 10  I=1,7
cxx 10    if( jj.gt.0 ) TDAY(JJ,I) = 12.0
       if( jj.gt.0 ) TDAY(JJ,I) = 12.0
 10   CONTINUE
C
      IE = IX(J)
c      write(6,*) IE, jj
      IF( IE .EQ. 0 )  GO TO 28
      I0 = I2
      DO 20  I=1,IE
      I2 = I0 + I
      IF( I2 .GT. 7 ) I2 = I2- 7
      IF( I2 .GT. 7 ) I2 = I2- 7
c      write(6,*) i2, tday(jj,i2)
cxx 20    if( jj.gt.0 ) TDAY(JJ,I2) = TDAY(JJ,I2) + 1.0
       if( jj.gt.0 ) TDAY(JJ,I2) = TDAY(JJ,I2) + 1.0
 20   CONTINUE
cc      I0 = I2
c      WRITE(6,*)  (TDAY(JJ,I),I=1,7)
 28   JJ = JJ + 1
      IF( JJ .GT.N ) RETURN
 30   continue
      GO TO 5
C
      E N D

cc      subroutine trpar( a,para )
      subroutine trpar( a,na,para,npa )
c      SUBROUTINE  TITLEP( TITLE,A )                                     
C                                                                       
C  ...  PLOT DATA ID AND ESTIMATED PARAMETERS  ...                      
C                                                                       
cxx      IMPLICIT REAL*8(A-H,O-Z)                                          
c      REAL*4   FM1, FM2, FM3, FM4, FM5, FAIC,FSIG2,TAU1,TAU2,TAU3       
c      REAL*4     TITLE(20), DAY(2), TIME(3)                             
cc      DIMENSION  A(40),para(26),a2(10),atmp(10)
cxx      DIMENSION  A(na),para(npa),ar(M2),atmp(M2)
cxx      INTEGER    PERIOD, SORDER
      INTEGER :: na, npa
      REAL(8) :: A(na), para(npa)
      INTEGER :: PERIOD, SORDER
      REAL(8) :: ar(M2), atmp(M2), tau1, tau2, tau3,
     1             DJACOB, FC, SIG2, AIC, FI, SIG2I, AICI, GI, GC,
     2             F1, F2, F3, A1, A2, A3, DI, UI, TDF
cc      COMMON  /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,                  
cc     *         NS,NI,MISING,IOUT,LL,N,NYEAR,nmonth,NPS,NPE,NPRED    
c      COMMON  /CMFUNC/  DJACOB, F, SIG2, AIC                           
cc      COMMON    /CMFUNC/  DJACOB,F,SIG2,AIC,FI,SIG2I,AICI,GI(20),G(20) 
cc      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), AR(10), A3(30)    
ccx      COMMON    /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
ccx      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)
      COMMON    /COMSM2/  M1, M2, M3, M4, M5, M, L, PERIOD, SORDER,
     *                    NYEAR, nmonth
      COMMON  /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(200),GC(200)
      COMMON  /COMSM3/  F1(10), F2(10), F3(300), A1(10), A2(10), A3(300)
     *                   ,DI, UI(3), TDF(7)
C                                                                       
       tau1 = 0.0D00
       tau2 = 0.0D00
       tau3 = 0.0D00


c      FM1 = M1                                                          
c      FM2 = M2                                                          
c      FM3 = M3                                                          
c      FM4 = M4                                                          
c      FM5 = M5
      para(1) = AIC                                                     
cc      para(2) = F                                                        
      para(2) = FC
      para(3) = SIG2                                                    
c      CALL  DATE( DAY )                                                 
c      CALL  CLOCK( TIME,1 )                                             
C     TAU1 = DEXP( A(1) )/(1.0D0 + DEXP(A(1)))                          
C     IF(L.GE.2)  TAU2 = DEXP( A(2) )/(1.0D0 + DEXP( A(2) ))            
C     IF(L.GE.3)  TAU3 = DEXP( A(3) )/(1.0D0 + DEXP( A(3) ))            
      TAU1 = 0.5D0*(1.0D0 + DSIN( A(1) )) + 0.0001

      IF(L.GE.2)  TAU2 = 0.5D0*(1.0D0 + DSIN( A(2) )) + 0.0001
C     IF(L.EQ.3)  TAU2 = DEXP( A(2) )                                   
      IF(L.GE.3)  TAU3 = 0.5D0*(1.0D0 + DSIN( A(3) )) + 0.0001          
C     TAU1 = A(1)**2                                                    
C     TAU2 = A(2)**2                                                    
C     TAU3 = A(3)**2                                                    
c  610 FORMAT(1H ,10F10.5 )                                              
C                                                                       
      para(4) = tau1
      para(5) = tau2
      para(6) = tau3
      IF( M2 .EQ. 0 )  GO TO 40                                         
      DO 30  I=1,M2                                                     
cxx   30 atmp(I) = 0.90D0*DSIN( A(L+I) )                                   
      atmp(I) = 0.90D0*DSIN( A(L+I) )
   30 CONTINUE
cc      CALL  ARCOEF( atmp,M2,A2 )
      CALL  ARCOEFD( atmp,M2,AR )
      do 35 I=1,M2
cc 35      para(i+6) = A2(i)
cxx 35      para(i+6) = AR(i)
         para(i+6) = AR(i)
 35   continue
 40   continue
      do 41 i=1,7
cxx 41      para(i+6+m2) = TDF(i)
         para(i+6+m2) = TDF(i)
 41   continue
C                                                                       
      RETURN                                                            
      E N D                                                             
