cx      SUBROUTINE  UNIMARF( ZS,N,LAG,ZMEAN,SUM,SD,AIC,DIC,M,AICM,SDM,A,
cx     *                     TMP,IER )
      SUBROUTINE  UNIMARF( ZS,N,LAG,ZMEAN,SUM,SD,AIC,DIC,M,AICM,SDM,A )
C
      INCLUDE 'timsac_f.h'
C
cc      PROGRAM  UNIMAR                                                   
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
C     TIMSAC 78.1.1.                                                    
C     ___                _                     __                       
C     UNIVARIATE CASE OF MINIMUM AIC METHOD OF AR MODEL FITTING.        
C                                                                       
C       THIS IS THE BASIC PROGRAM FOR THE FITTING OF AUTOREGRESSIVE MODE
C       OF SUCCESSIVELY HIGHER ORDERS BY THE METHOD OF LEAST SQUARES    
C       REALIZED THROUGH HOUSEHOLDER TRANSFORMATION.  THE OUTPUTS ARE   
C       THE ESTIMATES OF THE COEFFICIENTS, THE INNOVATION VARIANCES AND 
C       CORRESPONDING AIC STATISTICS.  AIC IS DEFINED BY                
C                                                                       
C               AIC  =  N * LOG( SD )  +  2 * ( NUMBER OF PARAMETERS )  
C       WHERE                                                           
C             N:    DATA LENGTH,                                        
C             SD:   ESTIMATE OF THE INNOVATION VARIANCE.                
C                                                                       
C       --------------------------------------------------------------- 
C       REFERENCE:                                                      
C          G.KITAGAWA AND H.AKAIKE(1978), "A PROCEDURE FOR THE MODELING 
C          OF NON-STATIONARY TIME SERIES.",  ANN. INST. STATIST. MATH., 
C          30,B,351-363.                                                
C       --------------------------------------------------------------- 
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
C             REDATA                                                    
C             REDUCT                                                    
C             ARMFIT                                                    
C       --------------------------------------------------------------- 
C       INPUTS REQUIRED:                                                
C             MT:    INPUT DEVICE SPECIFICATION (MT=5 : CARD READER)    
C             LAG:   UPPER LIMIT OF AR-ORDER, MUST BE LESS THAN 101     
C                                                                       
C       --  THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE REDATA  -- 
C                                                                       
C             TITLE: TITLE OF DATA                                      
C             N:     DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 10000   
C             DFORM: INPUT DATA FORMAT SPECIFICATION STATEMENT          
C                    -- FOR EXAMPLE --     (8F10.5)                     
C             (Z(I),I=1,N):  ORIGINAL DATA                              
C       --------------------------------------------------------------- 
C                                                                       
cc      !DEC$ ATTRIBUTES DLLEXPORT :: UNIMARF
C
cxx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL * 4   Z(10000) , TITLE(20)                                   
cc      REAL * 4   TITLE(20)
cc      DIMENSION  Z(10000)
cc      DIMENSION  X(200,101) , D(200) , A(100)                           
cxx      DIMENSION  ZS(N), Z(N)
cxx      DIMENSION  X(N+1,LAG+1), A(LAG)
cxx      DIMENSION  SD(LAG+1), AIC(LAG+1), DIC(LAG+1)
      INTEGER :: N, LAG, M
      REAL(8) :: ZS(N), ZMEAN, SUM, SD(LAG+1), AIC(LAG+1), DIC(LAG+1),
     1           AICM, SDM, A(LAG)
      REAL(8) :: Z(N), X(N+1,LAG+1)
cx      INTEGER*1  TMP(1)
cx      CHARACTER  CNAME*80
C                                                                       
C        EXTERNAL SUBROUTINE DECLARATION:                               
C                                                                       
      EXTERNAL  SETX1                                                   
C
cc      CHARACTER(100) IFLNAM,OFLNAM
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
cxcx            WRITE(*,*) ' ***  unimar temp FILE OPEN ERROR :',CNAME,IVAR
cx            IER=IVAR
cx            IFG=0
cx         END IF
cx      END IF
C                                                                       
C        PARAMETERS:                                                    
C             MJ1:  ABSOLUTE DIMENSION FOR SUBROUTINE CALL              
C             ISW:  =1   PARAMETERS OF MAICE MODEL ONLY ARE REQUESTED   
C                   =2   PARAMETERS OF ALL MODELS ARE REQUESTED         
C                                                                       
cc      MJ1 = 200                                                         
      MJ1 = N+1
      ISW = 2                                                           
C                                                                       
CC      READ( 5,1 )     MT                                                
cc      MT = 5
cc      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
cc      READ( 5,1 )     LAG                                               
cc      WRITE( 6,3 )                                                      
cc      WRITE( 6,2 )     LAG                                              
cc      WRITE( 6,4 )     MT                                               
C                                                                       
C          +-----------------------------------------+                +-
C          ! ORIGINAL DATA LOADING AND MEAN DELETION !                ! 
C          +-----------------------------------------+                +-
C                                                                       
cc      CALL  REDATA( Z,N,MT,TITLE )                                      
cc	CLOSE( MT )
      CALL  REDATA( ZS,Z,N,ZMEAN,SUM )                                      
      NMK = N - LAG                                                     
      K = LAG                                                           
C                                                                       
C          +-----------------------+                                  +-
C          ! HOUSEHOLDER REDUCTION !                                  ! 
C          +-----------------------+                                  +-
C                                                                       
cc      CALL  REDUCT( SETX1,Z,D,NMK,0,K,MJ1,LAG,X )                       
      CALL  REDUCT( SETX1,Z,NMK,0,K,MJ1,LAG,X )                       
C                                                                       
C          +------------------+                                       +-
C          ! AR MODEL FITTING !                                       ! 
C          +------------------+                                       +-
C                                                                       
cc      CALL  ARMFIT( X,K,LAG,NMK,ISW,TITLE,MJ1,A,SD,M )                  
cx      CALL  ARMFIT( X,K,LAG,NMK,ISW,MJ1,A,M,SD,AIC,DIC,SDM,AICM,
cx     *              IFG,LU )
      CALL  ARMFIT( X,K,LAG,NMK,ISW,MJ1,A,M,SD,AIC,DIC,SDM,AICM )
C
cc      GO TO 999
C                                                                      +
cc  900 CONTINUE
cc      WRITE(6,600) IVAR,OFLNAM
cc  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,//,5X,100A)
cc      GO TO 999
C                                                                      !
cc  910 CONTINUE
cc      IF ( NFL.EQ.2 ) CLOSE( 6 )
cc#ifdef __linux__
ccC      reopen #6 as stdout
cc      IF ( NFL.EQ.2 ) OPEN(6, FILE='/dev/fd/1')
cc#endif
ccC /* __linux__ */
cc      WRITE(6,610) IVAR,IFLNAM
cc  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//,5X,100A)
C                                                                      +
cc  999 CONTINUE
cx      IF (IFG.NE.0) CLOSE(LU)                                           
      RETURN
C                                                                       
cxx    1 FORMAT( 16I5 )                                                    
cxx    2 FORMAT( 1H ,'FITTING UP TO THE ORDER  K =',I3,'  IS TRIED' )      
cxx    3 FORMAT( ' PROGRAM TIMSAC 78.1.1',/'   AUTOREGRESSIVE MODEL FITTING
cxx     1  (SCALAR CASE)  ;   LEAST SQUARES METHOD BY HOUSEHOLDER TRANSFORM
cxx     2ATION',/,'   < AUTOREGRESSIVE MODEL >',/,1H ,10X,'Z(I) = A(1)*Z(I-
cxx     31) + A(2)*Z(I-2) +  ...  + A(M)*Z(I-M) + E(I)',/,'   WHERE',/,11X,
cxx     4'M:     ORDER OF THE MODEL',/,11X,'E(I):  GAUSSIAN WHITE NOISE WIT
cxx     5H MEAN 0  AND  VARIANCE SD(M).' )                                 
cxx    4 FORMAT( 1H ,'ORIGINAL DATA INPUT DEVICE  MT =',I3 )               
C                                                                       
      E N D                                                             
