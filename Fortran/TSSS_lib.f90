

!C     PROGRAM 7.1  ARFIT
!ccx      SUBROUTINE ARFITF( Y,N,LAG,NF,ISW,SIG2,AIC,MAR,A,PAR,SP )

      SUBROUTINE ARFITF( Y,N,LAG,NF,MJ2,ISW,SIG2,AIC,MAR,A,PAR,SP ) !GCC$ ATTRIBUTES HOT :: ARFITT !GCC$ ATTRIBUTES ALIGNED(32) :: ARFIT
        use omp_lib
        implicit none
!C
!      INCLUDE 'TSSS_f.h'
#if 0
!C
C  ...  AR model fitting  ...
C
C     Inputs:
C        IDEV:    Input device for time series
C        LAG:     Highest order of AR model
C        ISW:     Estimation procedure
C                 = 1:  Yule-Walker method
C                 = 2:  Least squares (Householder) method
C                 = 3:  Partial autoregression method
C                 = 4:  PARCOR method
C                 = 5:  Burg's algorithm (MEM)
C     The following inputs are required in the subroutine READTS.
C        TITLE:   Caption of the data set
C        N:       Data length
C        FORMAT:  Reading format
C        Y(I):    Time series, (i=1,...,N)
C     Parameters:
C        NMAX:    Adjustable dimension of Y (NMAX.GE.N)
C        MJ,MJ2:  Adjustable dimensions
C        NF:      Number of frequencies for computing spectrum
C     @TEST.PN71:  12/12/90,1/7/91 Y.I. and G.K. 9/4/91
C
cc      PARAMETER( NMAX=200,MJ=20,MJ2=100,NF=200,IDEV=1 )
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION Y(NMAX), COV(0:MJ,4), SIG2(0:MJ), SP(0:NF)
cc      DIMENSION D(MJ2), X(MJ2,MJ+1)
cc      DIMENSION AIC(0:MJ), A(MJ,MJ), B(MJ), PAR(MJ)
cc      DIMENSION FE(NMAX), BE(NMAX)
cxx      DIMENSION Y(N), COV(0:LAG,4), SIG2(0:LAG), SP(0:NF)
cxx      DIMENSION X(MJ2,LAG+1)
cxx      DIMENSION AIC(0:LAG), A(LAG,LAG), B(LAG+1), PAR(LAG)
cxx      DIMENSION FE(N), BE(N)
!c
#endif
      INTEGER :: N, LAG, NF, MJ2, ISW
      REAL(8) :: Y(N), SIG2(0:LAG), AIC(0:LAG), A(LAG,LAG), PAR(LAG), &
                SP(0:NF)
#if defined __ICC
      !DIR$ ASSUME_ALIGNED Y:64,SIG2:64,AIC:64,A:64,PAR:64,SP:64
#endif
      REAL(8) :: COV(0:LAG,4), X(MJ2,LAG+1), B(LAG+1), FE(N), BE(N),  &
                OUTMIN, OUTMAX, YMEAN
#if defined __ICC
     !DIR$ ATTRIBUTES ALIGN : 64 :: COV
     !DIR$ ATTRIBUTES ALIGN : 64 :: X
     !DIR$ ATTRIBUTES ALIGN : 64 :: B
     !DIR$ ATTRIBUTES ALIGN : 64 :: FE
     !DIR$ ATTRIBUTES ALIGN : 64 :: BE
#endif
!c
      EXTERNAL  SETXAR
      DATA  OUTMIN/-1.0e+30_8/, OUTMAX/1.0e+30_8/
!c
      X(1:MJ2,1:LAG+1) = 0.0_8
      PAR(1:LAG) = 0.0_8
!cc      READ(5,*)  LAG, ISW
!cc      CALL  READTS( IDEV,Y,N )
      CALL  MEAN( Y,N,-1.0D30,1.0D30,NSUM,YMEAN )
!C
!C  ...  Yule-Walker method  ...
!C
      IF( ISW.EQ.1 )  THEN
!cc         CALL  UNICOR( Y,N,LAG,OUTMIN,OUTMAX,COV )
         CALL  UNICOR( Y,N,LAG,OUTMIN,OUTMAX,COV,YMEAN )
         CALL  ARYULE( COV,N,LAG,SIG2,AIC,PAR,A,MAR)
      END IF
!C
!C  ...  Householder method  ...
!C
      IF( ISW.EQ.2 )  THEN
!cc         CALL  REDUCT( SETXAR,Y,D,N-LAG,0,LAG,MJ2,X )
!cc         CALL  REGRES( X,LAG,N-LAG,MJ2,MJ,A,SIG2,AIC,MAR )
         CALL  REDUCT( SETXAR,Y,N-LAG,0,LAG,MJ2,X )
         CALL  REGRES( X,LAG,N-LAG,MJ2,A,SIG2,AIC,MAR )
!cxx         CALL  PARCOR( A(1,MAR),MAR,PAR )
         CALL  PARCOR( A(1,LAG),LAG,PAR )
      END IF
!C
!C  ...  PARCOR method  ...
!C
      IF( ISW.GE.3 )  THEN
         CALL  ARPCOR( Y,FE,BE,SIG2,AIC,LAG,N,PAR,ISW-2,MAR )
         DO 10 I=1,LAG
!cxx   10    CALL  ARCOEF( PAR,I,A(1,I) )
            CALL  ARCOEF( PAR,I,A(1,I) )
   10    CONTINUE
      END IF
!C
!C  ...  Power spectrum  ...
!C
      CALL  ARMASP( A(1,MAR),MAR,B,0,SIG2(MAR),NF,SP )
!C
!cc      CALL  PRARSP( N,LAG,A,MAR,PAR,SIG2,AIC,SP,MJ,NF,ISW )
!C      CALL  PTAR( PAR,AIC,SP,LAG,NF,MAR,SIG2,ISW )
!cc      STOP
      RETURN
      E N D

      SUBROUTINE  ARPCOR( Y,FE,BE,SIG2,AIC,K,N,PARCOR,ISW,MAR ) !GCC$ ATTRIBUTES ALIGNED(32) :: ARPCOR !GCC$ ATTRIBUTES HOT :: ARPCOR !GCC$ ATTRIBUTES INLINE :: ARPCOR
        use omp_lib
        implicit none
#if 0
C
C  ...  PARCOR method  ...
C
C     Inputs:
C        Y(I):    Time series
C        N:       Data length
C        K:       Highest AR order
C        ISW:     Estimation procedure
C                 = 1:   Partial autoregression
C                 = 2:   PARCOR
C                 = 3:   Burg's algorithm (MEM)
C        FE,BE:   Working area
C     Outputs:
C        SIG2(I): Innovation variance
C        AIC(I):  AIC
C        PARCOR(I):  PARCOR
C        MAR:     MAICE order
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION  Y(N), FE(N), BE(N)
cxx      DIMENSION  SIG2(0:K), AIC(0:K), PARCOR(K)
cc      DIMENSION  A(50), B(50), FA(50), BA(50)
cxx      DIMENSION  A(K), B(K), FA(K), BA(K)
C
#endif
      INTEGER :: K, N, MAR
      REAL(8) :: Y(N), FE(N), BE(N), SIG2(0:K), AIC(0:K), PARCOR(K)
#if defined __ICC
      !DIR$ ASSUME_ALIGNED Y:64,FE:64,BE:64,SIG2:64,AIC:64,PARCOR:64
#endif
      REAL(8) :: A(K), B(K), FA(K), BA(K), PI, SUM, AICM, FB, FF, BB, &
               FE0, X
#if defined __ICC
      !DIR$ ATTRIBUTES ALIGN : 64 :: A,B,FA,BA
#endif
!C
      PI = 3.1415926535_8
!C
      SUM = 0.0_8
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
!$OMP SIMD REDUCTION(+:SUM) LINEAR(I:1)
#endif
      DO 5 I=K+1,N
!cxx    5 SUM = SUM + Y(I)**2
      SUM = SUM + Y(I)**2
    5 CONTINUE
      SIG2(0) = SUM/(N-K)
      AIC(0) = (N-K)*( DLOG(2*PI) + 1+ DLOG(SIG2(0))) + 2
      MAR = 0
      AICM = AIC(0)
!C
#if defined __ICC
      !DIR$ VECTOR ALIGNED
      !DIR$ SIMD 
#elif defined __GFORTRAN__
      !$OMP SIMD  LINEAR(I:1)
#endif
      DO 10 I=1,N
      FE(I) = Y(I)
!cxx   10 BE(I) = Y(I)
      BE(I) = Y(I)
   10 CONTINUE
!C
      DO 100 M=1,K
      FB = 0.0_8
      FF = 0.0_8
      BB = 0.0_8
#if defined __ICC
      !DIR$ VECTOR ALIGNED
      !DIR$ SIMD REDUCTION(+:FB,FF)
#elif defined __GFORTRAN__
      !$OMP SIMD REDUCTION(+:FB,FF) LINEAR(I:1)
#endif
      DO 20 I=M+1,N
      FB = FB + FE(I)*BE(I-M)
      FF = FF + BE(I-M)**2
!cxx   20 BB = BB + FE(I)**2
      BB = BB + FE(I)**2
   20 CONTINUE
!C
      IF( ISW.EQ.1) THEN
        A(M) = FB/FF
        B(M) = FB/BB
      END IF
      IF( ISW.EQ.2 ) THEN
        A(M) = FB/(DSQRT(FF*BB))
        B(M) = FB/(DSQRT(FF*BB))
      END IF
      IF( ISW.EQ.3 ) THEN
        A(M) = FB/((FF+BB)/2)
        B(M) = FB/((FF+BB)/2)
      END IF
!C
      
      DO 50 L=1,M-1
      A(L) = FA(L)-A(M)*BA(M-L)
!cxx   50 B(L) = BA(L)-B(M)*FA(M-L)
      B(L) = BA(L)-B(M)*FA(M-L)
   50 CONTINUE
#if defined __ICC
   !DIR$ VECTOR ALIGNED
   !DIR$ SIMD 
#elif defined __GFORTRAN__
   !$OMP SIMD  LINEAR(I:1)
#endif
      DO 60 L=1,M
      FA(L) = A(L)
!cxx   60 BA(L) = B(L)
      BA(L) = B(L)
   60 CONTINUE
      DO 70 I=M+1,N
      FE0 = FE(I)
      FE(I) = FE(I)-A(M)*BE(I-M)
!cxx   70 BE(I-M) = BE(I-M)-B(M)*FE0
      BE(I-M) = BE(I-M)-B(M)*FE0
   70 CONTINUE
!C
      SUM = 0.0_8
#if defined __ICC
      !DIR$ VECTOR ALIGNED
      !DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
      !$OMP SIMD REDUCTION(+:SUM) LINEAR(I:1)
#endif
      DO 80 I=K+1,N
!cxx   80 SUM = SUM + FE(I)**2
      SUM = SUM + FE(I)**2
   80 CONTINUE
      X = SUM
      PARCOR(M) = A(M)
      SIG2(M) = X/(N-K)
      AIC(M) = (N-K)*( DLOG(2*PI) + 1+ DLOG(SIG2(M)))+ 2*(M+1)
      IF( AIC(M).LT.AICM )  THEN
         MAR  = M
         AICM = AIC(M)
      END IF
!C
  100 CONTINUE
!C
      RETURN
      E N D


      SUBROUTINE ARMAF(M,L,A,B,SIG2,K,KMAX,NF,G,COV,PAR,SP, &
         ROOTA,ROOTB,IER,JER) !GCC$ ATTRIBUTES ALIGNED(32) :: ARMAF  !GCC$ ATTRIBUTES HOT :: ARMAF
                     implicit none
   !C
    !     INCLUDE 'TSSS_f.h'
#if 0
   C
   C     PROGRAM 6.1  ARMA
   C
   C  ...  This program analyses time series via ARMA modeling  ...
   C
   C     Inputs:
   C        IDEV:  Input devie (=5: CONSOLE)
   C        M:     AR order
   C        L:     MA order
   C        SIG2:  Innovation variance
   C        A(I):  AR coefficients
   C        B(I):  MA coefficients
   C        K:     Maximum lag of autocovariance function (K<MJ+1)
   C     Parameter:
   C        MJ:    Adjustable dimension of A, B, etc..
   C        NF:    Number of frequencies in evaluating spectrum
   C     Programmed by Y.I and G.K.
   C
   C     Outputs:
   C         IER:  =1 : MATRIX WITH ZERO ROW IN DECOMPOSE
   C               =2 : SINGULAR MATRIX IN DECOMPOSE.ZERO DIDIVIDE IN SOLVE
   C               =3 : CONVERGENCE IN IMPRUV.MATRIX IS NEARLY SINGULAR
   C         JER:  =1 : NON-CONVERGENCE AT POLYRT
   C
   cc      !DEC$ ATTRIBUTES DLLEXPORT::ARMAF
   C
   cc      PARAMETER( NF=200,MJ=100,IDEV=1 )
   cxx      IMPLICIT REAL*8( A-H,O-Z )
   cc      DIMENSION  A(MJ), B(MJ), PAR(MJ), ROOTA(MJ,2), ROOTB(MJ,2)
   cc      DIMENSION  G(0:MJ), COV(0:MJ), SP(0:NF)
   cc      DIMENSION  WRK1(0:MJ), WRK2(0:MJ), WRK3(MJ)
   cc      DATA  PAR/MJ*0.0D0/
   cxx      DIMENSION  A(M), B(L), PAR(K), ROOTA(M,2), ROOTB(L,2)
   cxx      DIMENSION  G(0:KMAX), COV(0:K), SP(0:NF)
   cxx      DIMENSION  WRK1(0:K), WRK2(0:K), WRK3(K,K)
#endif
         INTEGER :: M, L, K, KMAX, NF, IER, JER
         REAL(8) :: A(M), B(L), SIG2, G(0:KMAX), COV(0:K), PAR(K), &
                    SP(0:NF), ROOTA(M,2), ROOTB(L,2)
#if defined __ICC
         !DIR$ ASSUME_ALIGNED A:64,B:64,G:64,COV:64,PAR:64,SP:64,ROOTA:64,ROOTB:64
#endif
         REAL(8) :: WRK1(0:K), WRK2(0:K), WRK3(K,K)
#if defined __ICC 
         !DIR$ ATTRIBUTES ALIGN : 64 :: WRK1, WRK2, WRK3
#endif
#if 0
   C
   cc      WRITE( 6,* )  'K = ?'
   cc      READ( 5,* )   K
   cc      OPEN( IDEV,FILE='arma.dat' )
   cc      READ(IDEV,*)  M, L
   cc      READ(IDEV,*)  SIG2
   cc      READ(IDEV,*)  (A(I),I=1,M)
   cc      READ(IDEV,*)  (B(I),I=1,L)
   cc      CLOSE( IDEV )
   C
   cc      KMAX = MAX(M,L,K)
#endif
         CALL  IMPULS( M,L,A,B,K,G )
   !cc      CALL  ARMCOV( M,L,A,B,SIG2,K,COV )
         CALL  ARMCOV( M,L,A,B,SIG2,K,COV,KMAX,IER )
         if( ier.ne.0 ) return
   !c------------
         PAR(1:K) = 0.0D0
   !c------------
         CALL  PARCOR( A,M,PAR )
         CALL  ARCOEF( PAR,M,A )
   !cc      IF( L.GT.0 )  CALL  ARYULE( COV,1000,K,WRK1,WRK2,PAR,WRK3,MAR )
         IF( L.GT.0 )  CALL  ARYULE( COV,0,K,WRK1,WRK2,PAR,WRK3,MAR )
         CALL  ARMASP( A,M,B,L,SIG2,NF,SP )
  ! cc      CALL  CHROOT( A,M,ROOTA,MJ )
  ! cc      CALL  CHROOT( B,L,ROOTB,MJ )
         CALL  CHROOT( A,M,ROOTA,M,JER1 )
         CALL  CHROOT( B,L,ROOTB,L,JER2 )
         JER = JER1
         IF( JER2 .NE. 0 ) JER = JER + JER2 + 1
   !cc      CALL  PRARMA( M,L,A,B,G,K,COV,K,PAR,SP,NF,ROOTA,ROOTB,MJ )
   !C      CALL  PTARMA( G,K,COV,K,PAR,SP,NF,ROOTA,M,ROOTB,L,MJ )
         RETURN
         E N D

  ! cc      SUBROUTINE  CHROOT( A,M,ROOT,MJ )

         SUBROUTINE  CHROOT( A,M,ROOT,MJ,IER ) !GCC$ ATTRIBUTES INLINE :: CHROOT !GCC$ ATTRIBUTES ALIGNED(32) :: CHROOT
#if 0
   C
   C  ...  Characteristic roots of the AR or MA operator  ...
   C
   C     Inputs:
   C        A:     Vector of AR or MA coefficients
   C        M:     Order of the model
   C        MJ:    Adjustable dimension of ROOT
   C     Output:
   C        ROOT:  Characteristic roots (real part,imaginary part)
   C
   cxx      IMPLICIT  REAL*8(A-H,O-Z)
   cxx      DIMENSION  A(M), ROOT(MJ,2)
   cc      DIMENSION  C(50), CW(50)
   cxx      DIMENSION  C(M+1)
   C
#endif
         INTEGER :: M, MJ, IER
#if defined __ICC
         !DIR$ ASSUME_ALIGNED : 64 :: A, M
#endif
         REAL(8) :: A(M), ROOT(MJ,2)
         REAL(8) :: C(M+1)
#if defined __ICC
         !DIR$ ATTRIBUTES ALIGN : 64 :: C
#endif
  ! C
         IER = 0
  ! C
         IF( M.EQ.0 )  RETURN
         DO 10  I=1,M
   !cxx   10 C(I) = -A(M-I+1)
         C(I) = -A(M-I+1)
      10 CONTINUE
         C(M+1) = 1.0D0
         MMAX = M
  ! C
   !C  ... Characteristic roots of operator 1-A(1)*S- ... -A(M)*S**M=0  ...
  ! C
  ! cc      CALL  POLYRT( C,CW,MMAX,ROOT(1,1),ROOT(1,2),IER )
         CALL  POLYRT( C,MMAX,ROOT(1,1),ROOT(1,2),IER )
   !cc      IF( IER.NE.0 )   WRITE(6,600)
   !C
         RETURN
   !cxx  600 FORMAT( 1H0,'*****  NON-CONVERGENCE AT POLYRT  *****' )
         E N D

   
         SUBROUTINE  POLYRT( A,M,ROOTR,ROOTI,IER ) !GCC$ ATTRIBUTES ALIGNED(32) :: POLYRT !GCC$ ATTRIBUTES HOT :: POLYRT
            use omp_lib
            implicit none
#if 0
   C
   C  ...  This subroutine finds the roots of the equation
   C            A(1) + A(2)*S + ... + A(M)*(S**M) = 0
   C       by Newton-Raphson method
   C
   C     Inputs:
   C        A:     Coefficients of the equation
   C        B:     Working area
   C        M:     Degree of the polynomial
   C     Outputs:
   C        ROOTR:   Real parts of the roots
   C        ROOTI:   Imaginary parts of the roots
   C        IER:     Error code to indicate non-convergence
   C
   cxx      IMPLICIT  REAL*8(A-H,O-Z)
   cc      DIMENSION A(1), B(1), ROOTR(1), ROOTI(1)
   cxx      DIMENSION A(M+1), B(M+3), ROOTR(M), ROOTI(M)
   C
#endif
         INTEGER :: M, IER 
         REAL(8) :: A(M+1), ROOTR(M), ROOTI(M)
#if defined __ICC
         !DIR$ ASSUME_ALIGNED : 64 :: A,ROOTR,ROOTI
#endif
         REAL(8) :: B(M+3), XPR, YPR, X, XO, Y, YO, UX, UY, U, V, XT, YT, &
                   XT2, YT2, FI, SUM, DX, DY, TEM, ALPH
#if defined __ICC
          !DIR$ ATTRIBUTES ALIGN : 64 :: B
#endif
   
         IFIT = 0
         ISW = 0
         JSW = 0
         K = M
         IER = 0
         KX = K
         KXX = K+1
         KJI = KXX
         K2 = 1
   !c------------
         XPR = 0.0_8
         YPR = 0.0_8
   !c------------
         DO 10  I=1,KXX
         A(I) = A(I)/A(K+1)
         J = KXX - I+1
   !cxx   10 B(J) = A(I)
         B(J) = A(I)
      10 CONTINUE
      20 XO = 0.5e-02_8
         YO = 0.1e-01_8
         IN = 0
      30 X = XO
   !C
         XO = -10.0_8*YO
         YO = -10.0_8*X
         X = XO
         Y = YO
         ICT = 0
         IN = IN+1
         GO TO 50
   !C
      40 IFIT = 1
         XPR = X
         YPR = Y
   !C
      50 UX = 0.0_8
         UY = 0.0_8
         V  = 0._8
         YT = 0.0_8
         XT = 1.0_8
         U = B(K+1)
         IF( U .EQ. 0.0_8 )  GO TO 140
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ SIMD REDUCTION(+:U,V,UX) REDUCTION(-:UY)
#elif defined __GFORTRAN__
!$OMP SIMD REDUCTION(+:U,V,UX) REDUCTION(-:UY) LINEAR(I:1,L)
#endif
         DO 60  I=1,K
         L = K-I+1
         XT2 = X*XT - Y*YT
         YT2 = X*YT+Y*XT
         U = U + B(L)*XT2
         V = V + B(L)*YT2
         FI = I
         UX = UX + FI*XT*B(L)
         UY = UY - FI*YT*B(L)
         XT = XT2
         YT = YT2
      60 CONTINUE
         SUM = UX**2 + UY**2
         IF( SUM .EQ. 0.D0 )  GO TO 100
         DX = (V*UY - U*UX)/SUM
         X = X + DX
         DY = -(U*UY + V*UX)/SUM
         Y = Y + DY
         IF( DABS(DY)+DABS(DX) .LT. 1.0D-10 )  GO TO 80
   !C
         ICT = ICT+1
         IF( ICT .LT. 500 )  GO TO 50
         ISW = 1
         IF( IN .GE. 5 )  GO TO 70
         IF( IFIT .NE. 0 )  GO TO 80
      65 ISW = 0
         IFIT = 0
         GO TO 30
   !C
      70 IF( IFIT .EQ. 0 )  GO TO 300
         JSW = 1
      80 DO 90  L=1,KXX
         J = KJI-L+1
         TEM = A(J)
         A(J) = B(L)
   !cxx   90 B(L) = TEM
         B(L) = TEM
      90 CONTINUE
         ITEM = K
         K = KX
         KX = ITEM
     100 IF( IFIT .EQ. 0 )  GO TO 40
         IF( JSW .EQ. 1 )   GO TO 110
   !cc      IF( ISW-1 )  120,65,120
         IF( ISW-1 .LT. 0 )  GO TO 120
         IF( ISW-1 .EQ. 0 )  GO TO 65
         IF( ISW-1 .GT. 0 )  GO TO 120
     110 X = XPR
         Y = YPR
         ISW = 0
         JSW = 0
     120 IFIT = 0
         IF( X .EQ. 0.0_8 )  GO TO 130
         IF( DABS( Y/X ) .LT. 1.0D-08 )  GO TO 150
     130 ALPH = 2*X
         SUM = X*X + Y*Y
         K = K-2
         GO TO 160
   !C
     140 X = 0.0_8
         KX = KX-1
         KXX = KXX-1
     150 Y = 0.0_8
         SUM = 0.0_8
         ALPH = X
         K = K-1
     160 B(2) = B(2) + ALPH*B(1)
   !cxx  170 DO 180  L=2,K
         DO 180  L=2,K
   !cxx  180 B(L+1) = B(L+1) + ALPH*B(L) - SUM*B(L-1)
         B(L+1) = B(L+1) + ALPH*B(L) - SUM*B(L-1)
     180 CONTINUE
   !cc  190 ROOTI(K2) = Y
   !cc      ROOTR(K2) = X
     190 IF( K2.LE.M ) ROOTI(K2) = Y
         IF( K2.LE.M ) ROOTR(K2) = X
         K2 = K2+1
         IF( SUM .EQ. 0.D0 )  GO TO 200
         Y = -Y
         SUM = 0.0_8
         GO TO 190
     200 IF (K .LE. 0 )  RETURN
         GO TO 20
   !C
     300 IER = 1
         RETURN
         E N D

  
         SUBROUTINE ARMAFT2F( Y0,N,MMAX,LMAX,MLMAX,TSIG2,TFF,TAIC,AR,CMA,IER ) !GCC$ ATTRIBUTES HOT :: ARMAFT2F !GCC$ ATTRIBUTES ALIGNED(32) :: ARMAFT2F
#if defined __GFORTRAN__
            use omp_lib
#endif
            implicit none
 #if 0                            
   C
   C  ...  ARMA MODEL FITTING  ...
   C
   C     Inputs:
   C        M:       AR Order (M <= 20)
   C        L:       MA Order (M <= 20)
   C        IPARAM:  =0    Use defalt initail values
   C                 =1    Read intial values
   C        AR(I):   AR coefficients (I=1,M3)
   C     Parameters:
   C        NMAX:    Adjustable dimension of Y
   C        MJ:      Adjustable dimension of XF, VF, etc.
   C     @TEST.FILTER2O    NOV.29,1990, SEP.02,1992
   C
   cc      PARAMETER( NMAX=1000,MJ=40 )
   cc      IMPLICIT REAL*8(A-H,O-Z)
   cc      DIMENSION  AA(40), AR(40), PAR(40), CMA(40)
   cc      DIMENSION  TAIC(0:20,0:20), TFF(0:20,0:20)
   cc      DIMENSION  SAA(40,0:20,0:20)
   cc      COMMON  /C92825/  OUTMIN, OUTMAX
   cc      COMMON  /C92826/  Y(NMAX)
   cc      COMMON  /C92907/  ALIMIT
   cc      COMMON  /C92908/  M, L, N
   cc      COMMON  /C92909/  FLK, SIG2
   cc      COMMON  / CCC /  ISW, IPR, ISMT, IDIF
#endif
         INTEGER :: N, MMAX, LMAX, MLMAX, IER
         REAL(8) :: Y0(N), TSIG2(0:MMAX, 0:LMAX), TFF(0:MMAX,0:LMAX), &
                   TAIC(0:MMAX, 0:LMAX), AR(MMAX, 0:MMAX, 0:LMAX), &
                   CMA(LMAX, 0:MMAX, 0:LMAX)
#if defined __ICC
        !DIR$ ASSUME_ALIGNED Y0:64,TSIG2:64,TFF:64,TAIC:64,AR:64,CMA:64
#endif
         INTEGER :: IPRAM, IOPT, IDIF
         REAL(8) :: Y(N), PI, SUM, YMEAN, YVAR, ALIMIT, OUTMIN, OUTMAX, &
                   PAR(MLMAX), AA(MMAX+LMAX), SIG2, FLK, AIC, &
                   SAA(MMAX+LMAX, 0:MMAX, 0:LMAX)
#if defined __ICC
        !DIR$ ATTRIBUTES ALIGN : 64 :: Y,PAR,AA,SAA
#endif
   
         !EXTERNAL  FFARMA
         PI = 3.1415926535897932384626433_8
   !cc      IPR = 0
   !c
         TSIG2(0:MMAX, 0:LMAX) = 0.0_8
         TFF(0:MMAX,0:LMAX) = 0.0_8
         TAIC(0:MMAX, 0:LMAX) = 0.0_8
         AR(1:MMAX, 0:MMAX, 0:LMAX) = 0.0_8
         CMA(1:LMAX, 0:MMAX, 0:LMAX) = 0.0_8
         PAR(1:MLMAX) = 0.0_8
         AA(1:(MMAX+LMAX)) = 0.0_8
         SAA(1:(MMAX+LMAX), 0:MMAX, 0:LMAX) = 0.0_8
   !C
   !C  ...  Read Time Series  ...
   !C
   !cc      CALL  READTS( 1,Y,N )
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
!$OMP SIMD LINEAR(I:1)
#endif
         DO 5 I = 1,N
         Y(I) = Y0(I)
       5 CONTINUE
   !C
   !C  ...  Subtrac Mean Value  ...
   !C
         SUM = 0.0_8
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
!$OMP SIMD REDUCTION(+:SUM) LINEAR(I:1)
#endif
         DO 10 I=1,N
  ! cc   10 SUM = SUM + Y(I)
         SUM = SUM + Y(I)
      10 CONTINUE
         YMEAN = SUM/N
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#elif defined __GFORTRAN__
         !$OMP SIMD  LINEAR(I:1)
#endif
         DO 20 I=1,N
  ! cc   20 Y(I) = Y(I) - YMEAN
         Y(I) = Y(I) - YMEAN
      20 CONTINUE
   
   !cc      DO 25 I=0,20
   !cc      DO 25 J=0,20
   !cc   25 TFF(I,J) = -1.0D10
         DO 26 I=0,MMAX
         DO 25 J=0,LMAX
         TFF(I,J) = -1.0E+10_8
      25 CONTINUE
      26 CONTINUE
   !C
   !C  ...  Read Model Orders  ...
   !C
   !C      READ( 5,* )  M, L, IPARAM
   !cc      open( 3,file='temp.dat' )
   !cc      MMAX = 10
   !ccc     LMAX = 10
  ! cc      DO 100  M=0,MMAX
         DO 101  M=0,MMAX
         DO 100  L=0,LMAX
         
   !cc      write(3,*) M, L
         IF( M.EQ.0 .and. L.EQ.0 )  THEN
            SUM = 0.0_8
#if defined __ICC
            !DIR$ VECTOR ALIGNED
            !DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
            !$OMP SIMD REDUCTION(+:SUM) LINEAR(I:1)
#endif
            DO 110 I=1,N
   !cc  110    SUM = SUM + Y(I)**2
            SUM = SUM + Y(I)**2
     110    CONTINUE
            YVAR = SUM/N
            TFF(0,0) = -0.5_8*N*(DLOG(2*PI*YVAR) + 1)
            TAIC(0,0) = -2*TFF(0,0) + 2
            GO TO 100
         END IF
   !C
   !C  ...  Set Defalt Parameters  ...
   !C
   !cc      IPR   = 0
         IPRAM = 0
         IOPT  = 1
         IDIF  = 2
         ALIMIT = 0.95_8
         
         IF( L.EQ.0 .OR. M.EQ.0 ) THEN 
         
  ! cc      CALL  SPARA1( M,L,AR,CMA,OUTMIN,OUTMAX,IOPT )
         CALL  SPARA1( M,L,MLMAX,AR(1,M,L),CMA(1,M,L),OUTMIN,OUTMAX,IOPT )
#if 0
   C      IF( IPARAM.EQ.1 )  THEN
   C         READ( 5,* )  (AR(I),I=1,M)
   C         READ( 5,* )  (CMA(I),I=1,L)
   C      END IF
   C
   C      WRITE(6,*) M, L
   C      WRITE(6,*) (AR(I),I=1,M)
   C      WRITE(6,*) (CMA(I),I=1,L)
   C
   cc      CALL  PARCOR( AR,M,PAR )
#endif
         CALL  PARCOR( AR(1,M,L),M,PAR )
  ! C      write(6,*) (PAR(I),I=1,M)
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
!$OMP SIMD LINEAR(I:1)
#endif
         DO 30 I=1,M
  ! cc   30 AA(I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
         AA(I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
      30 CONTINUE
  ! cc      CALL  PARCOR( CMA,L,PAR )
         CALL  PARCOR( CMA(1,M,L),L,PAR )
  ! C      write(6,*) (PAR(I),I=1,L)
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 40 I=1,L
  ! cc   40 AA(M+I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
         AA(M+I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
      40 CONTINUE
      
         ELSE
            AA(M+L) = 0.0000_8
  ! C         IF( M.EQ.0.AND.L.EQ.1 )  AA(1) = 0.18D0
  ! C         IF( M.GT.0 )  THEN
            IF( TFF(M-1,L).GT.TFF(M,L-1) )  THEN
               DO 70 I=1,M-1
  ! cc   70       AA(I) = SAA(I,M-1,L)
               AA(I) = SAA(I,M-1,L) 
      70       CONTINUE
              AA(M) = 0.0000_8
               DO 80 I=1,L
  ! cc   80       AA(M+I) = SAA(M-1+I,M-1,L)         
               AA(M+I) = SAA(M-1+I,M-1,L)         
      80       CONTINUE         
            END IF
  ! C         END IF
         END IF
         
  ! cc      WRITE(6,*) M, L
   !cc      WRITE(6,666) (AA(I),I=1,M+L)
   !cc  666 FORMAT( 10F10.5 )
     
  ! C
  ! C  ...  Maximum Likelihood Method  ...
  ! C
  ! cc      write(3,*)  (AA(I),I=1,M+L)
  ! cc      IF( IOPT.EQ.1 )  CALL  DAVIDN( FFARMA,AA,M+L,2 )
         ier = 0
         IF( IOPT.EQ.1 )  CALL  DAVIDN( FFARMA,AA,M+L,2, &
            Y,N,M,L,MLMAX,OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
         if( ier.ne.0 ) return
  ! C
  ! C
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 50 I=1,M
  ! cc   50 PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0D0)/(DEXP(AA(I))+1.0D0)
         PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0_8)/(DEXP(AA(I))+1.0_8)
      50 CONTINUE
  ! cc      CALL  ARCOEF( PAR,M,AR )
         CALL  ARCOEF( PAR,M,AR(1,M,L) )
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 60 I=1,L
  ! cc   60 PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0D0)/(DEXP(AA(M+I))+1.0D0)
         PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0_8)/(DEXP(AA(M+I))+1.0_8)
      60 CONTINUE
  ! cc      CALL  ARCOEF( PAR,L,CMA )
         CALL  ARCOEF( PAR,L,CMA(1,M,L) )
  ! C
         AIC = -2*FLK + 2*(M+L+1)
  ! C
  ! C  ...  Print out the Maximum Likelihood Estimates  ...
  ! C
  ! cc      WRITE(6,666) (AA(I),I=1,M+L)
   
  ! cc      CALL  PRARMA( M,L,AR,CMA,SIG2,FLK,AIC )
         TFF(M,L) = FLK
         TAIC(M,L) = AIC
         TSIG2(M,L) = SIG2
         DO 90 I=1,M+L
  ! cc   90 SAA(I,M,L) = AA(I)
         SAA(I,M,L) = AA(I)
      90 CONTINUE
  ! cc      WRITE(3,33)  FLK, AIC
  ! cc   33 FORMAT( 5f12.4 )
  ! cc      WRITE(3,*)  (AA(I),I=1,M+L)
     100 CONTINUE
     101 CONTINUE
#if 0
   C
   cc      open( 2,file='temp_table.dat' )
   cc      WRITE(6,*)  'log-Likelihood(M,L)'
   cc      WRITE(2,*)  'log-Likelihood(M,L)'
   cc      DO 200 I=0,MMAX
   cc      WRITE(2,600)  (TFF(I,J),J=0,LMAX)
   cc  200 WRITE(6,600)  (TFF(I,J),J=0,LMAX)
   cc      WRITE(2,*)  'AIC(M,L)'
   cc      WRITE(6,*)  'AIC(M,L)'
   cc      DO 210 I=0,MMAX
   cc      WRITE(2,600)  (TAIC(I,J),J=0,LMAX)
   cc  210 WRITE(6,600)  (TAIC(I,J),J=0,LMAX)
   cc      close( 2 )
   cc      close( 3 )
   cc  600 FORMAT( 11F10.4 )
   cc      STOP
#endif
         RETURN
         E N D
    
   
         
    SUBROUTINE ARMAFTF( Y0,N,M,L,MLMAX,IPARAM,AR0,CMA0, &
                            SIG2,FLK,AIC,AR,CMA,IER ) !GCC$ ATTRIBUTES HOT :: ARMAFTF !GCC$ ATTRIBUTES ALIGNED(32) :: ARMAFTF
            use omp_lib
            implicit none
   !C
    !     INCLUDE 'TSSS_f.h'
#if 0
   C
   C  ...  ARMA MODEL FITTING  ...
   C
   C     Inputs:
   C        L:       AR Order (M <= 20)
   C        M:       MA Order (M <= 20)
   C        IPARAM:  =0    Use defalt initail values
   C                 =1    Read intial values
   C        AR(I):   AR coefficients (I=1,M3)
   C     Parameters:
   C        NMAX:    Adjustable dimension of Y
   C        MJ:      Adjustable dimension of XF, VF, etc.
   C     @TEST.FILTER2O    NOV.29,1990, SEP.02,1992
   C
   cc      PARAMETER( NMAX=1000,MJ=20 )
   cxx      IMPLICIT REAL*8(A-H,O-Z)
   cc      DIMENSION  AA(20), AR(20), PAR(20), CMA(20)
   cxx      DIMENSION  AR0(M), CMA0(L)
   cxx      DIMENSION  AA(M+L), AR(M), PAR(MLMAX), CMA(L)
   cxx      DIMENSION  Y0(N), Y(N)
   cc      COMMON  /C92825/  OUTMIN, OUTMAX
   cc      COMMON  /C92826/  Y(NMAX)
   cc      COMMON  /C92907/  ALIMIT
   cc      COMMON  /C92908/  M, L, N
   cc      COMMON  /C92909/  FLK, SIG2
   cc      COMMON  / CCC /  ISW, IPR, ISMT, IDIF
   C
#endif
         INTEGER :: N, M, L, MLMAX, IPARAM, IER 
         REAL(8) :: Y0(N), AR0(M), CMA0(L), SIG2, FLK, AIC, AR(M), CMA(L)
#if defined __ICC
         !DIR$ ASSUME_ALIGNED Y0:64,AR0:64,CMA0:64,AR:64,CMA:64
#endif
         REAL(8) :: Y(N), AA(M+L), PAR(MLMAX), ALIMIT, OUTMIN, OUTMAX, SUM, &
                   YMEAN
#if defined __ICC
         !DIR$ ATTRIBUTES ALIGN : 64 :: Y,AA,PAR
#endif
        ! EXTERNAL  FFARMA
#if 0
   C
   C  ...  Read Model Orders  ...
   C
   cc      READ( 5,* )  M, L, IPARAM
   C
   C  ...  Set Defalt Parameters  ...
   C
   cc      IPR  = 7
#endif
         ALIMIT = 0.95_8
   !cc      CALL  SPARA1( M,L,AR,CMA,OUTMIN,OUTMAX,IOPT )
         CALL  SPARA1( M,L,MLMAX,AR,CMA,OUTMIN,OUTMAX,IOPT )
         IF( IPARAM.EQ.1 )  THEN
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
!$OMP SIMD LINEAR(I:1)
#endif
            DO 5 I=1,M
   !cxx    5    AR(I) = AR0(I)
            AR(I) = AR0(I)
       5    CONTINUE
#if defined __ICC
       !DIR$ VECTOR ALIGNED
       !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
       !$OMP SIMD LINEAR(I:1)
#endif
            DO 6 I=1,L
   !cxx    6    CMA(I) = CMA0(I)
            CMA(I) = CMA0(I)
       6    CONTINUE
  ! cc         READ( 5,* )  (AR(I),I=1,M)
  ! cc         READ( 5,* )  (CMA(I),I=1,L)
         END IF
  ! C
  ! C  ...  Read Time Series  ...
  ! C
  ! cc      CALL  READTS( 1,Y,N )
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 7 I = 1,N
   !cxx    7 Y(I) = Y0(I)
         Y(I) = Y0(I)
       7 CONTINUE
   !C
   !C  ...  Subtrac Mean Value  ...
   !C
         SUM = 0.0_8
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
!$OMP SIMD REDUCTION(+:SUM) LINEAR(I:1)
#endif
         DO 10 I=1,N
   !cxx   10 SUM = SUM + Y(I)
         SUM = SUM + Y(I)
      10 CONTINUE
         YMEAN = SUM/N
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 20 I=1,N
   !cxx   20 Y(I) = Y(I) - YMEAN
         Y(I) = Y(I) - YMEAN
      20 CONTINUE
   !C
  ! cc      WRITE(6,*) M, L
  ! cc      WRITE(6,*) (AR(I),I=1,M)
  ! C
         CALL  PARCOR( AR,M,PAR )
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 30 I=1,M
  ! cxx   30 AA(I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
         AA(I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
      30 CONTINUE
         CALL  PARCOR( CMA,L,PAR )
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 40 I=1,L
   !cxx   40 AA(M+I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
         AA(M+I) = DLOG( (ALIMIT+PAR(I))/(ALIMIT-PAR(I)) )
      40 CONTINUE
  ! C
  ! C  ...  Maximum Likelihood Method  ...
  ! C
  ! cc      IF( IOPT.EQ.1 )  CALL  DAVIDN( FFARMA,AA,M+L,2 )
         ier = 0
         IF( IOPT.EQ.1 )  CALL  DAVIDN( FFARMA,AA,M+L,2, &
            Y,N,M,L,MLMAX,OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
         if( ier.ne.0 ) return
  ! C
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 50 I=1,M
  ! cxx   50 PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0D0)/(DEXP(AA(I))+1.0D0)
         PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0D0)/(DEXP(AA(I))+1.0D0)
      50 CONTINUE
         CALL  ARCOEF( PAR,M,AR )
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 60 I=1,L
  ! cxx   60 PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0D0)/(DEXP(AA(M+I))+1.0D0)
         PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0D0)/(DEXP(AA(M+I))+1.0D0)
      60 CONTINUE
         CALL  ARCOEF( PAR,L,CMA )
  ! C
         AIC = -2*FLK + 2*(M+L+1)
  ! C
  ! C  ...  Print out the Maximum Likelihood Estimates  ...
  ! C
  ! cc      CALL  PRARMA( M,L,AR,CMA,SIG2,FLK,AIC )
  ! C
  ! cc      STOP
         RETURN
         E N D
   
         SUBROUTINE  SPARA1( M,L,MLMAX,AR,CMA,OUTMIN,OUTMAX,IOPT ) !GCC$ ATTRIBUTES INLINE :: SPARA1 !GCC$ ATTRIBUTES ALIGNED(32) :: SPARA1
#if defined __GFORTRAN__
            use omp_lib
#endif
            implicit none
#if 0
   C
   C  ...  Set Default Parameters  ...
   C
   C     Inputs:
   C       M:       AR order
   C       L:       MA order
   C     Outputs:
   C       TAU2:    System noise variance
   C       AR:      AR coefficients
   C       OUTMIN:  Lower bound for outliers
   C       OUTMAX:  Upper bound for outliers
   C       IOPT:    (=1  MLE by numerical optimization)
   C
   cxx      IMPLICIT REAL*8(A-H,O-Z)
   cc      DIMENSION  AR(20), PAR(20), CMA(20)
   cxx      DIMENSION  AR(M), PAR(MLMAX), CMA(L)
   C
#endif
         INTEGER :: M, L, MLMAX, IOPT
         REAL(8) :: AR(M), CMA(L), OUTMIN, OUTMAX
#if defined __ICC
         !DIR$ ASSUME_ALIGNED AR:64,CMA:64
#endif
         REAL(8) :: PAR(MLMAX)
#if defined __ICC
         !DIR$ ATTRIBUTES ALIGN : 64 :: PAR
#endif
   !C
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
!$OMP SIMD LINEAR(I:1)
#endif
         DO 10 I=1,M
   !cxx   10 PAR(I) = -(-0.6D0)**I
         PAR(I) = -(-0.6_8)**I
      10 CONTINUE
         CALL  ARCOEF( PAR,M,AR )
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 20 I=1,L
  ! cxx   20 PAR(I) = -(-0.6D0)**I
  ! ccc      PAR(I) = -(-0.6D0)**I
         PAR(I) = -(-0.5_8)**I
      20 CONTINUE
         CALL  ARCOEF( PAR,L,CMA )
   !C
         OUTMIN = -1.0e+30_8
         OUTMAX =  1.0e+30_8
         IOPT = 1
   !C
         RETURN
         E N D

   
         SUBROUTINE  FILTR3( Y,XF,VF,A,B,M,NS,N,OUTMIN,OUTMAX, &
                            FF,OVAR ) !GCC$ ATTRIBUTES ALIGNED(32) :: FILTR3 !GCC$ ATTRIBUTES HOT :: FILTR3
 #if defined __GFORTRAN__
            use omp_lib
#endif
            implicit none
#if 0               
   C
   C  ...  Kalman filter  ...
   C
   C     Inputs:
   C        Y:      time series
   C        N:      data length
   C        NS:     Start position of filtering
   C        XF:     Initial state vector
   C        VF:     Initial covariance matrix
   C        A:      Parameters of the matrix F
   C        B:      Parameters of the matrix G
   C        C:      Parameters of the matrix H
   C        K:      Dimension of the state vector
   C        MJ:     Adjustable dimension of XF, VF
   C        OUTMIN: Lower limit for detecting outliers
   C        OUTMAX: Upper limit for detecting outliers
   C     Outputs:
   C        FF:     Log likelihood
   C        SIG2:   Estimated variance
   C
   cxx      IMPLICIT REAL*8(A-H,O-Z)
   cc      DIMENSION  A(MJ), B(MJ), C(MJ), Y(N)
   cc      DIMENSION  XF(MJ), VF(MJ,MJ), XP(40), VP(40,40)
   cc      DIMENSION  WRK(40,40), VH(40), GAIN(40)
   cxx      DIMENSION  A(M), B(M), C(M), Y(N)
   cxx      DIMENSION  XF(M), VF(M,M), XP(M), VP(M,M)
   cxx      DIMENSION  WRK(M,M), VH(M), GAIN(M)
   C
#endif
         INTEGER :: M, NS, N
         REAL(8) :: Y(N), VF(M,M), A(M), B(M), OUTMIN, OUTMAX, FF, OVAR
#if defined __ICC
         !DIR$ ASSUME_ALIGNED Y:64,VF:64,A:64,B:64
#endif
         REAL(8) :: XF(M), XP(M), VP(M,M), WRK(M,M), VH(M), GAIN(M), PI, &
                    SDET, PVAR, PERR
#if defined __ICC
         !DIR$ ATTRIBUTES ALIGN : 64 XF,XP,VP,WRK,VH,GAIN
#endif
    
   !C
         DATA   PI  /3.1415926535897932384626433_8/
   !C
         OVAR = 0.0_8
         SDET = 0.0_8
         NSUM = 0
   !C
         DO 300  II=NS,N
   !C
   !C  ...  ONE STEP AHEAD PREDICTION  ...
  ! C
         XP(M) = A(M)*XF(1)
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
!$OMP SIMD LINEAR(I:1)
#endif
         DO 100  I=1,M-1
  ! cxx  100 XP(I) = A(I)*XF(1) + XF(I+1)
         XP(I) = A(I)*XF(1) + XF(I+1)
     100 CONTINUE
  ! C
  ! cxx      DO 110  J=1,M
         DO 111  J=1,M
         WRK(M,J) = A(M)*VF(1,J)
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 110  I=1,M-1
   !cxx  110 WRK(I,J) = A(I)*VF(1,J) + VF(I+1,J)
         WRK(I,J) = A(I)*VF(1,J) + VF(I+1,J)
     110 CONTINUE
     111 CONTINUE
   !C
   !cxx      DO 120  I=1,M
         DO 121  I=1,M
         VP(I,M) = WRK(I,1)*A(M)
         DO 120  J=1,M-1
  ! cxx  120 VP(I,J) = WRK(I,1)*A(J) + WRK(I,J+1)
         VP(I,J) = WRK(I,1)*A(J) + WRK(I,J+1)
     120 CONTINUE
     121 CONTINUE
   !C
  ! cxx      DO 140  I=1,M
         DO 141  I=1,M
#if defined __ICC
            !DIR$ VECTOR ALIGNED
            !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
            !$OMP SIMD LINEAR(I:1)
#endif
         DO 140  J=1,M
  ! cxx  140 VP(I,J) = VP(I,J) + B(I)*B(J)
         VP(I,J) = VP(I,J) + B(I)*B(J)
     140 CONTINUE
     141 CONTINUE
  ! C
  ! C  ...  FILTERING  ...
  ! C
         IF( Y(II).GT.OUTMIN .AND. Y(II).LT.OUTMAX .AND. II.LE.N )  THEN
  ! C
#if defined __ICC
            !DIR$ VECTOR ALIGNED
            !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
            !$OMP SIMD LINEAR(I:1)
#endif
         DO 210  I=1,M
   !cxx  210 VH(I) = VP(I,1)
         VH(I) = VP(I,1)
     210 CONTINUE
  ! C
         PVAR = VH(1)
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         IF( PVAR.LE.1.0e-30_8 )  THEN
            FF = -1.0e+20_8
            RETURN
         END IF
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         PERR = Y(II) - XP(1)
  ! C
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 230  I=1,M
  ! cxx  230 GAIN(I) = VH(I)/PVAR
         GAIN(I) = VH(I)/PVAR
     230 CONTINUE
  ! C
#if defined __ICC
     !DIR$ VECTOR ALIGNED
     !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
     !$OMP SIMD LINEAR(I:1)
#endif
         DO 250  I=1,M
  ! cxx  250 XF(I) = XP(I) + GAIN(I)*PERR
         XF(I) = XP(I) + GAIN(I)*PERR
     250 CONTINUE
  ! C
  ! cxx      DO 260  J=1,M
         DO 261  J=1,M
#if defined __ICC
            !DIR$ VECTOR ALIGNED
            !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
            !$OMP SIMD LINEAR(I:1)
#endif
         DO 260  I=1,M
  ! cxx  260 VF(I,J) = VP(I,J) - GAIN(I)*VH(J)
         VF(I,J) = VP(I,J) - GAIN(I)*VH(J)
     260 CONTINUE
     261 CONTINUE
  ! C
         OVAR = OVAR + PERR**2/PVAR
         SDET = SDET + DLOG(PVAR)
         NSUM = NSUM + 1
   !C
   !C  ...  MISSING OBSERVATION  ...
   !C
         ELSE
   !cxx      DO 270  I=1,M
         DO 271  I=1,M
         XF(I) = XP(I)
         DO 270  J=1,M
   !cxx  270 VF(I,J) = VP(I,J)
         VF(I,J) = VP(I,J)
     270 CONTINUE
     271 CONTINUE
         END IF
  ! C
     300 CONTINUE
         OVAR = OVAR/NSUM
         FF = -0.5_8*(NSUM*DLOG(PI*2*OVAR) + SDET + NSUM)
   !C
         RETURN
         E N D
   
         SUBROUTINE  ISTAT3( M,L,MM,AR,CMA,XF,VF,IER ) !GCC$ ATTRIBUTES ALIGNED(32) :: ISTAT3 !GCC$ ATTRIBUTES HOT :: ISTAT3
#if defined __GFORTRAN__
          use omp_lib
#endif
          implicit none
#if 0
   C
   C  ...  Initial state ...
   C
   C     Inputs:
   C        M:     AR order
   C        L:     MA order
   C        AR:    AR coefficient
   C        CMA:   MA coefficient
   C        MJ:    Adjustable dimension of F
   C     Outputs:
   C         XF:   State vector, X(0|0)
   C         VF:   State covarance matrix, V(0|0)
   C
   cxx      IMPLICIT REAL*8(A-H,O-Z)
   cc      DIMENSION  AR(*), CMA(*)
   cc      DIMENSION  XF(MJ), VF(MJ,MJ), COV(0:20), G(0:20)
   cxx      DIMENSION  AR(M), CMA(L)
   cxx      DIMENSION  XF(MM), VF(MM,MM), COV(0:MM), G(0:MM)
   C
#endif

         INTEGER :: M, L, MM, IER 
         REAL(8) :: AR(M), CMA(L), XF(MM), VF(MM,MM)
#if defined __ICC
         !DIR$ ASSUME_ALIGNED AR:64,CMA:64,XF:64,VF:64
#endif
         REAL(8) :: COV(0:MM), G(0:MM), SUM
#if defined __ICC
         !DIR$ ATTRIBUTES ALIGN : 64 :: COV,G
#endif
  ! C
  ! cc      MM = MAX0( M,L+1 )
  ! cc      DO 10  I=1,MJ
  ! cxx      DO 10  I=1,MM
         DO 11  I=1,MM
         XF(I) = 0.0_8
  ! cc      DO 10  J=1,MJ
         DO 10  J=1,MM
  ! cxx   10 VF(I,J) = 0.0D0
         VF(I,J) = 0.0_8
      10 CONTINUE
      11 CONTINUE
  ! C
  ! cc      CALL  ARMCOV( M,L,AR,CMA,1.0D0,MM,COV )
         CALL  ARMCOV( M,L,AR,CMA,1.0D0,MM,COV,MM,IER )
         if( ier.ne.0 ) return
         CALL  IMPULS( M,L,AR,CMA,MM,G )
         VF(1,1) = COV(0)
         DO 50 I=2,MM
         SUM = 0.0_8
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
!$OMP SIMD REDUCTION(+:SUM) LINEAR(J:1)
#endif
         DO 30 J=I,M
   !cxx   30 SUM = SUM + AR(J)*COV(J-I+1)
         SUM = SUM + AR(J)*COV(J-I+1)
      30 CONTINUE
#if defined __ICC
      !DIR$ VECTOR ALIGNED
      !DIR$ SIMD REDUCTION(-:SUM)
#elif defined __GFORTRAN__
      !$OMP SIMD REDUCTION(-:SUM) LINEAR(J:1)
#endif
         DO 40 J=I-1,L
   !cxx   40 SUM = SUM - CMA(J)*G(J-I+1)
         SUM = SUM - CMA(J)*G(J-I+1)
      40 CONTINUE
         VF(1,I) = SUM
   !cxx   50 VF(I,1) = SUM
         VF(I,1) = SUM
      50 CONTINUE
   !cxx      DO 100 I=2,MM
         DO 101 I=2,MM
         DO 100 J=I,MM
         SUM = 0.0_8
  ! cxx      DO 60 I1=I,M
         DO 61 I1=I,M
#if defined __ICC
            !DIR$ VECTOR ALIGNED
            !DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
            !$OMP SIMD REDUCTION(+:SUM) LINEAR(J1:1)
#endif
         DO 60 J1=J,M
  ! cxx   60 SUM = SUM + AR(I1)*AR(J1)*COV(IABS(J1-J-I1+I))
         SUM = SUM + AR(I1)*AR(J1)*COV(IABS(J1-J-I1+I))
      60 CONTINUE
      61 CONTINUE
   !cxx      DO 70 I1=I,M
         DO 71 I1=I,M
  ! c  modified 2019/12/09 =============
  ! ccc   DO 70 J1=J-I+I1,L
         JMIN = MAX(J-1,J-I+I1)
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD REDUCTION(-:SUM)
   #elif defined __GFORTRAN__
         !$OMP SIMD REDUCTION(-:SUM) LINEAR(J1:1)
   #endif
         DO 70 J1=JMIN,L
   !c ==================================
   !cxx   70 SUM = SUM - AR(I1)*CMA(J1)*G(IABS(J1-J-I1+I))
         SUM = SUM - AR(I1)*CMA(J1)*G(IABS(J1-J-I1+I))
      70 CONTINUE
      71 CONTINUE
   !cxx      DO 80 I1=J,M
         DO 81 I1=J,M
   !c  modified 2019/12/09 =============
   !ccc   DO 80 J1=I-J+I1,L
         JMIN = MAX( I-1,I-J+I1)
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD REDUCTION(-:SUM)
#elif defined __GFORTRAN__
         !$OMP SIMD REDUCTION(-:SUM) LINEAR(J1:1)
#endif
         DO 80 J1=JMIN,L
   !c ==================================
   !cxx   80 SUM = SUM - AR(I1)*CMA(J1)*G(IABS(J1-I-I1+J))
         SUM = SUM - AR(I1)*CMA(J1)*G(IABS(J1-I-I1+J))
      80 CONTINUE
      81 CONTINUE
  ! c  modified 2019/12/09 =============
   !ccc      DO 90 I1=I-1,L
         DO 90 I1=I-1,L+I-J
  ! ccc      SUM = SUM + CMA(I1)*CMA(J1)
         SUM = SUM + CMA(I1)*CMA(J-I+I1)
      90 CONTINUE
   !c ==================================
         VF(I,J) = SUM
   !cxx  100 VF(J,I) = SUM
         VF(J,I) = SUM
     100 CONTINUE
     101 CONTINUE
   !C
         RETURN
         E N D

   
         SUBROUTINE  SETABC( M,L,AR,CMA,A,B,C,MM ) !GCC$ ATTRIBUTES INLINE :: SETABC !GCC$ ATTRIBUTES ALIGNED(32) :: SETABC
            implicit none
#if 0
   C
   C  ...  State space model for Seasonal Adjustment  ...
   C
   C     Input:
   C       M:     AR Order
   C       L:     MA Order
   C       AR:    AR coefficient
   C       CMA:   MA coefficient
   C     Outputs:
   C       A,B,C: Parameters of F, G, H
   C
   cxx      IMPLICIT REAL*8(A-H,O-Z)
   cc      DIMENSION  A(*), B(*), C(*)
   cc      DIMENSION  AR(*), CMA(*)
   cxx      DIMENSION  A(MM), B(MM), C(MM)
   cxx      DIMENSION  AR(M), CMA(L)
   C
#endif
         INTEGER :: M, L, MM
         REAL(8) :: AR(M), CMA(L), A(MM), B(MM), C(MM)
#if defined __ICC
         !DIR$ ASSUME_ALIGNED AR:64,CMA:64,A:64,B:64,C:64
#endif
   !C
   !cc      MM = MAX0( M,L+1 )
   !cxx      DO 10  I=1,MM
   !cxx      C(I) = 0.0D0
   !cxx      A(I) = 0.0D0
   !cxx   10 B(I) = 0.0D0
         C(1:MM) = 0.0_8
         A(1:MM) = 0.0_8
         B(1:MM) = 0.0_8
   
   !C
         DO 20  I=1,M
   !cxx   20 A(I) = AR(I)
         A(I) = AR(I)
      20 CONTINUE
         B(1) = 1.0_8
         DO 30 I=1,L
   !cxx   30 B(I+1) =-CMA(I)
         B(I+1) =-CMA(I)
      30 CONTINUE
         C(1) = 1.0_8
  ! C
         RETURN
         E N D

   
         SUBROUTINE  FFARMA( K,AA,FF,IFG, &
                     Y,N,M,L,MM,OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER ) !GCC$ ATTRIBUTES HOT :: FFARMA !GCC$ ATTRIBUTES ALIGNED(32) :: FFARMA
#if defined __GFORTRAN__
                use omp_lib
#endif
            implicit none
#if 0
   C
   C  ... Function for Maximum likelihood estimation (seasonal model) ...
   C
   C     Inputs:
   C        K:      Number of parameters
   C        AA(I):  Parameter vector
   C     Outputs:
   C        FF:     -(Log likelihood)
   C        IFG:    =1 if some conditions are violated
   C
   cxx      IMPLICIT REAL*8(A-H,O-Z)
   cc      DIMENSION  PAR(20), GDUMMY(20)
   cc      DIMENSION  AA(20), AR(20), CMA(20)
   cc      DIMENSION  A(20), B(20), C(20)
   cc      DIMENSION  XF(40), VF(40,40)
   cc      COMMON  /C92826/   Y(1000)
   cc      COMMON  /C92825/   OUTMIN, OUTMAX
   cc      COMMON  /C92907/  ALIMIT
   cc      COMMON  /C92908/   M, L, N
   cc      COMMON  /C92909/   FLK, SIG2
   cxx      DIMENSION  PAR(MM)
   cxx      DIMENSION  AA(K), AR(M), CMA(L)
   cxx      DIMENSION  A(MM), B(MM), C(MM)
   cxx      DIMENSION  XF(MM), VF(MM,MM)
   cxx      DIMENSION  Y(N)
   C
#endif
         INTEGER :: K, IFG, N, M, L, MM, IER
         REAL(8) :: AA(K), FF, Y(N), OUTMIN, OUTMAX, ALIMIT, FLK, SIG2
#if defined __ICC
         !DIR$ ASSUME_ALIGNED AA:64,Y:64
#endif
         REAL(8) :: PAR(MM), AR(M), CMA(L), A(MM), B(MM), C(MM), XF(MM), &
                    VF(MM,MM)
#if defined __ICC
         !DIR$ ATTRIBUTES ALIGN : 64 :: PAR,AR,CMA,A,B,C,XF,VF
#endif
  ! C
  ! cc      MJ  = 20
  ! C
 !  cc      MM = MAX0( M,L+1 )
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
!$OMP SIMD LINEAR(I:1)
#endif
         DO 10 I=1,M
  ! cxx   10 PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0D0)/(DEXP(AA(I))+1.0D0)
         PAR(I)  = ALIMIT*(DEXP(AA(I))-1.0D0)/(DEXP(AA(I))+1.0D0)
      10 CONTINUE
         CALL  ARCOEF( PAR,M,AR )
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 20 I=1,L
  ! cxx   20 PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0D0)/(DEXP(AA(M+I))+1.0D0)
         PAR(I)  = ALIMIT*(DEXP(AA(M+I))-1.0D0)/(DEXP(AA(M+I))+1.0D0)
      20 CONTINUE
         CALL  ARCOEF( PAR,L,CMA )
         IFG = 0
  ! C
  ! cc      CALL  SETABC( M,L,AR,CMA,A,B,C )
  ! cc      CALL  ISTAT3( M,L,AR,CMA,MJ,XF,VF )
  ! cc      CALL  FILTR3( Y,XF,VF,A,B,C,MM,MJ,1,N,OUTMIN,OUTMAX,FLK,SIG2 )
         CALL  SETABC( M,L,AR,CMA,A,B,C,MM )
         CALL  ISTAT3( M,L,MM,AR,CMA,XF,VF,IER )
         if( ier.ne.0 ) return
  ! cxx      CALL  FILTR3( Y,XF,VF,A,B,C,MM,1,N,OUTMIN,OUTMAX,FLK,SIG2 )
         CALL  FILTR3( Y,XF,VF,A,B,MM,1,N,OUTMIN,OUTMAX,FLK,SIG2 )
         FF = -FLK
         RETURN
  ! C
  ! cx  100 IFG = 1
  ! cx      FF = 1.0D20
  ! cx      RETURN
         E N D

   
         SUBROUTINE  DAVIDN( FUNCT,X,N,NDIF, YY,NN,M,L,MLMAX, &
                           OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER ) !GCC$ ATTRIBUTES HOT :: DAVIDN !GCC$ ATTRIBUTES ALIGNED(32) :: DAVIDN
#if defined __GFORTRAN__
                    use omp_lib
#endif
                implicit none
#if 0
   C
   C  ...  6/20/83, 12/19/92
   C
   cxx      IMPLICIT  REAL*8( A-H,O-Z )
   cc      DIMENSION  X(40), DX(40), G(40), G0(40), Y(40)
   cc      DIMENSION  H(40,40), WRK(40), S(40)
   cc      COMMON     / CCC /  ISW, IPR, ISMT, IDIF
   cc      COMMON     / DDD /  XM , AIC , SD
   cxx      DIMENSION  X(N), DX(N), G(N), G0(N), Y(N)
   cxx      DIMENSION  H(N,N), WRK(N), S(N)
   cxx      DIMENSION  YY(NN)
   C
#endif
         INTEGER :: N, NDIF, NN, M, L, MLMAX, IER
         REAL(8) :: X(N), YY(NN), OUTMIN, OUTMAX, ALIMIT, FLK, SIG2
#if defined __ICC
         !DIR$ ASSUME_ALIGNED X:64,YY:64
#endif
         REAL(8) :: DX(N), G(N), G0(N), Y(N), H(N,N), WRK(N), S(N), TAU2, &
                   EPS1, EPS2, RAMDA, CONST1, SUM, S1, S2, STEM, SS, ED,  &
                   XM, XMB
#if defined __ICC
         !DIR$ ATTRIBUTES ALIGN : 64 :: DX,G,G0,Y,H,WRK,S
#endif
   
         EXTERNAL  FUNCT
         DATA        TAU2  /          1.0e-6_8  /
         DATA  EPS1, EPS2  / 1.0e-6_8 , 1.0e-6_8  /
         RAMDA  = -0.5_8
         CONST1 = 1.0e-30_8
         IPR = 0
 !  c-------   
         IDIF = NDIF
 !  c-------
  ! C
  ! C          INITIAL ESTIMATE OF INVERSE OF HESSIAN
 !  C
         ICOUNT = 0
    1000 CONTINUE
#if 0
   cxx 1000 DO 20  I=1,N
   cxx      DO 10  J=1,N
   cxx   10 H(I,J) = 0.0D0
   cxx      S(I)   = 0.0D0
   cxx      DX(I)  = 0.0D0
   cxx   20 H(I,I) = 1.0D0
#endif
         H(1:N,1:N) = 0.0_8
         S(1:N)   = 0.0_8
         DX(1:N)  = 0.0_8
         DO 20  I=1,N
         H(I,I) = 1.0_8
      20 CONTINUE
         ISW = 0
 !  C
 !  cc      IF( NDIF.EQ.0 )  CALL  FUNCT( N,X,XM,G,IG )
  ! cc      IF( NDIF.GE.1 )  CALL  FUNCND( FUNCT,N,X,XM,G,IG )
  ! cxx      IF( NDIF.EQ.0 ) CALL  FUNCT ( N,X,XM,G,IG, YY,NN,M,L,MLMAX,
         IF( NDIF.EQ.0 ) CALL  FUNCT ( N,X,XM,IG, YY,NN,M,L,MLMAX,
                            OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER ) &
         IF( NDIF.GE.1 ) CALL  FUNCND( FUNCT,N,X,XM,G,IG,YY,NN,M,L,MLMAX, &
                           OUTMIN,OUTMAX,ALIMIT,ISW,IDIF,FLK,SIG2,IER )
         if( ier.ne.0 ) return
   !C
   !cc      IF( IPR .GE. 2 )   WRITE( 6,640 )   XM, SD, AIC
   !cc      WRITE(6,650) (X(I),I=1,N), XM, (G(I),I=1,N)
   !C
         ICC = 0
  ! C      ITERATION
    2000 CONTINUE
         ICC = ICC + 1
         IF( ICC .EQ. 1 )  GO TO 120
   !C
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
!$OMP SIMD LINEAR(I:1)
#endif
         DO 40  I=1,N
  ! cxx   40 Y(I) = G(I) - G0(I)
         Y(I) = G(I) - G0(I)
      40 CONTINUE
         DO 60  I=1,N
         SUM = 0.0_8
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
         !$OMP SIMD REDUCTION(+:SUM) LINEAR(I:1)
#endif
         DO 50  J=1,N
  ! cxx   50 SUM = SUM + Y(J)*H(I,J)
         SUM = SUM + Y(J)*H(I,J)
      50 CONTINUE
   !cxx   60 WRK(I) = SUM
         WRK(I) = SUM
      60 CONTINUE
         S1 = 0.0_8
         S2 = 0.0_8
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD REDUCTION(+:S1,S2)
#elif defined __GFORTRAN__
         !$OMP SIMD REDUCTION(+:S1,S2) LINEAR(I:1)
#endif
         DO 70  I=1,N
         S1 = S1 + WRK(I)*Y(I)
   !cxx   70 S2 = S2 + DX(I) *Y(I)
         S2 = S2 + DX(I) *Y(I)
      70 CONTINUE
         IF( S1.LE.CONST1 .OR. S2.LE.CONST1 )  GO TO 900
#if 0
   C     IF( S1 .LE. S2 )   GO TO 100
   C
   C          UPDATE THE INVERSE OF HESSIAN MATRIX
   C
   C               ---  BROYDEN-FLETCHER-GOLDFARB-SHANNO TYPE CORRECTION  -
   C
   cxx  100 CONTINUE
#endif
         STEM = S1 / S2 + 1.0D0
   !cxx      DO 110  I=1,N
         DO 111  I=1,N
         DO 110  J=I,N
         H(I,J) = H(I,J)-(DX(I)*WRK(J)+WRK(I)*DX(J)-DX(I)*DX(J)*STEM)/S2
   !cxx  110 H(J,I) = H(I,J)
         H(J,I) = H(I,J)
     110 CONTINUE
     111 CONTINUE
   !C
     120 CONTINUE
         SS = 0.0_8
         DO 140  I=1,N
         SUM = 0.0_8
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
         !$OMP SIMD REDUCTION(+:SUM) LINEAR(J:1)
#endif
         DO 130  J=1,N
   !cxx  130 SUM = SUM + H(I,J)*G(J)
         SUM = SUM + H(I,J)*G(J)
     130 CONTINUE
         SS  = SS + SUM**2
  ! cxx  140 S(I) = -SUM
         S(I) = -SUM
     140 CONTINUE
   !C
         S1 = 0.0_8
         S2 = 0.0_8
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD REDUCTION(+:S1,S2)
#elif defined __GFORTRAN__
         !$OMP SIMD REDUCTION(+:S1,S2) LINEAR(I:1)
#endif
         DO 150  I=1,N
         S1 = S1 + S(I)*G(I)
   !cxx  150 S2 = S2 + G(I)*G(I)
         S2 = S2 + G(I)*G(I)
     150 CONTINUE
  ! C     DS2 = DSQRT(S2)
  ! C     GTEM = DABS(S1)/DS2
  ! C     IF( GTEM.LE.TAU1 .AND. DS2.LE.TAU2 )  GO TO 900
         IF( S1.LT.0.0_8 )  GO TO 200
         DO 170  I=1,N
         DO 160  J=1,N
   !cxx  160 H(I,J) = 0.0D0
         H(I,J) = 0.0_8
     160 CONTINUE
         H(I,I) = 1.0_8
  ! cxx  170 S(I) = -S(I)
         S(I) = -S(I)
     170 CONTINUE
     200 CONTINUE
   !C
         ED = XM
  ! C
   !C          LINEAR  SEARCH
   !C
   !cc      CALL  LINEAR( FUNCT,X,S,RAMDA,ED,N,IG )
         CALL  LINEAR( FUNCT,X,S,RAMDA,ED,N,IG, YY,NN,M,L, MLMAX,OUTMIN,OUTMAX,ALIMIT,ISW,FLK,SIG2,IER )
   !cxx     *         MLMAX,OUTMIN,OUTMAX,ALIMIT,ISW,IDIF,FLK,SIG2,IER )
               
         if( ier.ne.0 ) return
   !C
   !cc      IF( IPR .GE. 2 )  WRITE( 6,630 )  RAMDA, ED, SD, AIC
  ! C
         S1 = 0.0_8
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD REDUCTION(+:S1)
#elif defined __GFORTRAN__
         !$OMP SIMD REDUCTION(+:S1) LINEAR(I:1)
#endif
         DO 210  I=1,N
         DX(I) = S(I)*RAMDA
         S1 = S1 + DX(I)*DX(I)
         G0(I) = G(I)
   !cxx  210 X(I) = X(I) + DX(I)
         X(I) = X(I) + DX(I)
     210 CONTINUE
         XMB = XM
         ISW = 0
  ! C
  ! cc      IF( NDIF.EQ.0 )  CALL  FUNCT( N,X,XM,G,IG )
  ! cc      IF( NDIF.GE.1 )  CALL  FUNCND( FUNCT,N,X,XM,G,IG )
  ! cxx      IF( NDIF.EQ.0 )  CALL  FUNCT ( N,X,XM,G,IG, YY,NN,M,L,MLMAX,
         IF( NDIF.EQ.0 )  CALL  FUNCT ( N,X,XM,IG, YY,NN,M,L,MLMAX, &
                             OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
         IF( NDIF.GE.1 )  CALL  FUNCND( FUNCT,N,X,XM,G,IG,YY,NN,M,L,MLMAX, &
                             OUTMIN,OUTMAX,ALIMIT,ISW,IDIF,FLK,SIG2,IER )
         if( ier.ne.0 ) return
  ! cc      WRITE(6,650) (X(I),I=1,N), XM, (G(I),I=1,N)
   !C
         S2 = 0.0_8
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD REDUCTION(+:S2)
#elif defined __GFORTRAN__
         !$OMP SIMD REDUCTION(+:S2) LINEAR(I:1)
#endif
         DO 220  I=1,N
   !cxx  220 S2 = S2 + G(I)*G(I)
         S2 = S2 + G(I)*G(I)
     220 CONTINUE
         IF( DSQRT(S2) .LT. TAU2 )  GO TO 900
         IF( XMB/XM-1._8.LT.EPS1 .AND. DSQRT(S1).LT.EPS2 )  GO TO 900
   !C     IF( ICC .GE. 5 )  GO TO 900
         GO TO 2000
     900 CONTINUE
         IF( IPR .LE. 0 )  RETURN
   !cc      WRITE( 6,600 )
   !cc      WRITE( 6,610 )  (X(I),I=1,N)
   !cc      WRITE( 6,620 )
   !cc      WRITE( 6,610 )  (G(I),I=1,N)
   !C
         ICOUNT  = ICOUNT + 1
         S2 = 0.0_8
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD REDUCTION(+:S2)
#elif defined __GFORTRAN__
         !$OMP SIMD REDUCTION(+:S2) LINEAR(I:1)
#endif
         DO 910  I=1,N
   !cxx  910 S2 = S2 + G(I)*G(I)
         S2 = S2 + G(I)*G(I)
     910 CONTINUE
         IF( S2.GT.1.0D0.AND.ICOUNT.LE.1 )   GO TO 1000
         RETURN
#if 0
   cxx  600 FORMAT( 1H0,'-----  X  -----' )
   cxx  610 FORMAT( 1H ,10D13.5 )
   cxx  620 FORMAT( 1H0,'***  GRADIENT  ***' )
   cxx  630 FORMAT( 1H ,'LAMBDA =',D15.7,3X,'(-1)LOG LIKELIHOOD =',D23.15,
   cxx     *        3X,'SD =',D22.15,5X,'AIC =',D23.15 )
   cxx  640 FORMAT( 1H ,26X,'(-1)LOG-LIKELIHOOD =',D23.15,3X,'SD =',D22.15,
   cxx     *        5X,'AIC =',D23.15 )
   cxx  650 FORMAT( 10X,5F12.7 )
#endif
         E N D
   
         SUBROUTINE  FUNCND( FUNCT,M,A,F,G,IFG, Y,N,MM,L,MLMAX,  &
                            OUTMIN,OUTMAX,ALIMIT,ISW,IDIF,FLK,SIG2,IER ) !GCC$ ATTRIBUTES HOT :: FUNCND !GCC$ ATTRIBUTES ALIGNED(32) :: FUNCND
#if 0
   C
   C  ...  FUNCTION EVALUATION AND NUMERICAL DIFFERENCING  ...
   C
   cxx      IMPLICIT   REAL*8( A-H,O-Z )
   cc      DIMENSION  A(M) , G(M) , B(20),GD(5)
   cc      COMMON  / CCC /  ISW , IPR, ISMT, IDIF
   cc      COMMON  /CMFUNC/ DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
   cxx      DIMENSION  A(M) , G(M) , B(M)
   cxx      DIMENSION  Y(N)
   C
#endif
         INTEGER :: M, IFG, N, MM, L, MLMAX, ISW, IDIF, IER
         REAL(8) :: A(M), F, G(M), Y(N), OUTMIN, OUTMAX, ALIMIT, FLK, SIG2
         REAL(8) :: B(M), CONST, FB, FF
         EXTERNAL FUNCT
  ! C
  ! C     DATA       ICNT /0/
         CONST = 0.00001_8
  ! C
  ! cc      CALL  FUNCT( M,A,F,GD,IFG )
  ! cxx      CALL  FUNCT( M,A,F,GD,IFG, Y,N,MM,L,MLMAX,
         CALL  FUNCT( M,A,F,IFG, Y,N,MM,L,MLMAX, &
                            OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
         if( ier.ne.0 ) return
         FB = F
         IF( ISW .GE. 1 )   RETURN
   !C
   !C     WRITE( 6,600 )   (A(I),I=1,M)
         DO 10  I=1,M
   !cxx   10 B(I) = A(I)
         B(I) = A(I)
      10 CONTINUE
  ! C
         DO 30  II=1,M
         B(II) = A(II) + CONST
  ! cc      CALL  FUNCT( M,B,FF,GD,IFG )
  ! cxx      CALL  FUNCT( M,B,FF,GD,IFG, Y,N,MM,L,MLMAX,
         CALL  FUNCT( M,B,FF,IFG, Y,N,MM,L,MLMAX, &
                            OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
         if( ier.ne.0 ) return
         IF( IDIF .EQ. 1 )  GO TO 20
         B(II) = A(II) - CONST
   !cc      CALL  FUNCT( M,B,FB,GD,IFG )
  ! cxx      CALL  FUNCT( M,B,FB,GD,IFG, Y,N,MM,L,MLMAX,
         CALL  FUNCT( M,B,FB,IFG, Y,N,MM,L,MLMAX, &
                            OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
      20 G(II) = (FF-FB)/(CONST*IDIF)
         IF( G(II) .GT. 1.0e+20_8 )  G(II) = (F-FB)/CONST
         IF( G(II) .LT.-1.0e+20_8 )  G(II) = (FF-F)/CONST
         IF( FB.GT.F .AND. FF.GT.F )  G(II) = 0.0_8
   cxx   30 B(II) = A(II)
         B(II) = A(II)
      30 CONTINUE
   
         RETURN
   
  ! cxx  600 FORMAT( 3X,'---  PARAMETER  ---',(/,3X,5D13.5) )
  ! cxx  610 FORMAT( 3X,'---  GRADIENT  ---',(/,3X,5D13.5) )
         E N D
  
         SUBROUTINE  LINEAR( FUNCT,X,H,RAM,EE,K,IG, Y,N,M,L, &
            MLMAX,OUTMIN,OUTMAX,ALIMIT,ISW,FLK,SIG2,IER ) !GCC$ ATTRIBUTES HOT :: LINEAR !GCC$ ATTRIBUTES ALIGNED(32) :: LINEAR
#if defined __GFORTRAN__
              use omp_lib
#endif
              implicit none
 #if 0  
   C
   C  ...  LINEAR SEARCH  ...
   C
   cxx      IMPLICIT  REAL*8( A-H,O-Z )
   cc      INTEGER  RETURN,SUB
   cc      DIMENSION  X(1), H(1), X1(K)
   cc      DIMENSION  G(40)
   cc      COMMON     / CCC /  ISW , IPR, ISMT, IDIF
   cxx      DIMENSION  X(K), H(K), X1(K)
   cxx      DIMENSION Y(N)
   C
#endif
         INTEGER :: K, IG, N, M, L, MLMAX, ISW, IER
         REAL(8) :: X(K), H(K), RAM, EE, Y(N), OUTMIN, OUTMAX, ALIMIT, FLK, &
                    SIG2
#if defined __ICC
        !DIR$ ASSUME_ALIGNED X:64,H:64,Y:64
         INTEGER :: ire510
         REAL(8) :: X1(K), CONST2, HNORM, RAM1, RAM2, RAM3, E1, E2, E3, &
                   A1, A2, A3, B1, B2
#if defined __ICC
        !DIR$ ATTRIBUTES ALIGN : 64 :: X1
#endif
         EXTERNAL FUNCT
  
         ISW = 1
         IG  = 0
         RAM = 0.1_8
         CONST2 = 1.0e-60_8
      14 DO 15 I=1,K
   !cxx   15 IF( DABS(H(I)) .GT. 1.0D10 )  GO TO 16
         IF( DABS(H(I)) .GT. 1.0e+10_8 )  GO TO 16
      15 CONTINUE
         GO TO 18
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR ALIGNED
#elif defined __GFORTRAN__
!$OMP SIMD LINEAR(I:1)
#endif
      16 DO 17 I=1,K
   !cxx   17 H(I) = H(I)*1.0D-10
         H(I) = H(I)*1.0e-10_8
      17 CONTINUE
         GO TO 14
      18 CONTINUE
         HNORM = 0.0_8
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ SIMD REDUCTION(+:HNORM)
#elif defined __GFORTRAN__
!$OMP SIMD REDUCTION(+:HNORM)
#endif
         DO 10  I=1,K
   !cxx   10 HNORM = HNORM + H(I)**2
         HNORM = HNORM + H(I)**2
      10 CONTINUE
         HNORM = DSQRT( HNORM )
   !C
         RAM2 = RAM
         E1 =EE
         RAM1 = 0.0_8
   !C
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALIGNED
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 20  I=1,K
  ! cxx   20 X1(I) = X(I) + RAM2*H(I)
         X1(I) = X(I) + RAM2*H(I)
      20 CONTINUE
  ! cc      CALL  FUNCT( K,X1,E2,G,IG )
   !cxx      CALL  FUNCT( K,X1,E2,G,IG, Y,N,M,L,MLMAX,
         CALL  FUNCT( K,X1,E2,IG, Y,N,M,L,MLMAX, &
                           OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
         if( ier.ne.0 ) return
   !cc      IF(IPR.GE.7)  WRITE(6,2)  RAM2,E2
   !C
         IF( IG .EQ. 1  )  GO TO 50
         IF( E2 .GT. E1 )  GO TO 50
      30 RAM3 = RAM2*4.0_8
#if defined __ICC
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALIGNED
#elif defined __GFORTRAN__
      !$OMP SIMD LINEAR(I:1)
#endif
         DO 40  I=1,K
  ! cxx   40 X1(I) = X(I) + RAM3*H(I)
         X1(I) = X(I) + RAM3*H(I)
      40 CONTINUE
   !cc      CALL  FUNCT( K,X1,E3,G,IG )
   !cxx      CALL  FUNCT( K,X1,E3,G,IG, Y,N,M,L,MLMAX,
         CALL  FUNCT( K,X1,E3,IG, Y,N,M,L,MLMAX, &
                            OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
         if( ier.ne.0 ) return
          IF( IG.EQ.1 )  GO TO  500
  ! cc      IF( IPR.GE.7 )  WRITE(6,3)  RAM3,E3
         IF( E3 .GT. E2 )  GO TO 70
         IF(RAM3.GT.1.0D10 .AND. E3.LT.E1)  GO TO 45
         IF(RAM3.GT.1.0D10 .AND. E3.GE.E1)  GO TO 46
         RAM1 = RAM2
         RAM2 = RAM3
         E1 = E2
         E2 = E3
         GO TO 30
   
      45 RAM = RAM3
         EE = E3
         RETURN
   
      46 RAM = 0.0_8
         RETURN
   
      50 RAM3 = RAM2
         E3 = E2
         RAM2 = RAM3*0.1_8
         IF( RAM2*HNORM .LT. CONST2 )  GO TO  400
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALIGNED
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
         DO 60  I=1,K
   !cxx   60 X1(I) = X(I) + RAM2*H(I)
         X1(I) = X(I) + RAM2*H(I)
      60 CONTINUE
   !cc      CALL  FUNCT( K,X1,E2,G,IG )
  ! cxx      CALL  FUNCT( K,X1,E2,G,IG, Y,N,M,L,MLMAX,
         CALL  FUNCT( K,X1,E2,IG, Y,N,M,L,MLMAX, &
                           OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
         if( ier.ne.0 ) return
   !cc      IF(IPR.GE.7)  WRITE(6,4)  RAM2,E2
         IF( E2.GT.E1 )  GO TO 50
   !C
   !cc   70 ASSIGN 80 TO RETURN
      70 CONTINUE
         IRET = 80
         GO TO 200
  ! C
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALIGNED
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
      80 DO 90  I=1,K
   !cxx   90 X1(I) = X(I) + RAM*H(I)
         X1(I) = X(I) + RAM*H(I)
      90 CONTINUE
  ! cc      CALL  FUNCT( K,X1,EE,G,IG )
  ! cxx      CALL  FUNCT( K,X1,EE,G,IG, Y,N,M,L,MLMAX,
         CALL  FUNCT( K,X1,EE,IG, Y,N,M,L,MLMAX, &
                            OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
         if( ier.ne.0 ) return
   !cc      IF(IPR.GE.7)  WRITE(6,5)  RAM,EE
   !C
         IFG = 0
  ! cc      ASSIGN  300 TO  SUB
   !cc      ASSIGN 200 TO SUB
   !cc   95 ASSIGN 130 TO RETURN
         ISUB = 300
         ISUB = 200
      95 CONTINUE
         IRET = 130
         IF( RAM .GT. RAM2 )  GO TO 110
         IF( EE .GE. E2 )  GO TO 100
         RAM3 = RAM2
         RAM2 = RAM
         E3 =E2
         E2 =EE
  ! cc      GO TO  SUB,( 200,300 )
         IF( ISUB.EQ.200 ) GO TO 200
         IF( ISUB.EQ.300 ) GO TO 300
   !C
     100 RAM1 = RAM
         E1 = EE
  ! cc      GO TO  SUB,( 200,300 )
         IF( ISUB.EQ.200 ) GO TO 200
         IF( ISUB.EQ.300 ) GO TO 300
  ! C
     110 IF( EE .LE. E2 )  GO TO 120
         RAM3 = RAM
         E3 = EE
   !cc      GO TO  SUB,( 200,300 )
         IF( ISUB.EQ.200 ) GO TO 200
         IF( ISUB.EQ.300 ) GO TO 300
  ! C
     120 RAM1 = RAM2
         RAM2 = RAM
         E1 = E2
         E2 = EE
  ! cc      GO TO  SUB,( 200,300 )
         IF( ISUB.EQ.200 ) GO TO 200
         IF( ISUB.EQ.300 ) GO TO 300
  ! C
#if defined __ICC
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALIGNED
#elif defined __GFORTRAN__
         !$OMP SIMD LINEAR(I:1)
#endif
     130 DO 140  I=1,K
   !cxx  140 X1(I) = X(I) + RAM*H(I)
         X1(I) = X(I) + RAM*H(I)
     140 CONTINUE
  ! cc      CALL  FUNCT( K,X1,EE,G,IG )
  ! cxx      CALL  FUNCT( K,X1,EE,G,IG, Y,N,M,L,MLMAX,
         CALL  FUNCT( K,X1,EE,IG, Y,N,M,L,MLMAX, &
                            OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
         if( ier.ne.0 ) return
   !cc      IF( IPR.GE.7 )  WRITE(6,6)  RAM,EE
  ! cc      ASSIGN 200 TO SUB
         ISUB = 200
         IFG = IFG+1
  ! cxx  600 FORMAT( 1H ,6D20.13 )
   !C     IF( HNORM*(RAM3-RAM1) .GT. 1.0D-12 )  GO TO 95
         IF( RAM3-RAM1 .GT. RAM2*1.0D-1 )  GO TO 95
   !C     IF( E1-EE .GT. 1.0D-15 )  GO TO 95
   !C     IF( E3-EE .GT. 1.0D-15 )  GO TO 95
   !C     IF(DX.LE.C1*XN+C2 .AND. DF.LE.C1*DABS(EE)+C2)  RETURN
   !C     IF( IFG .LE. 2 )  GO TO 95
   !C
         IF( E2 .LT. EE )  RAM = RAM2
         RETURN
   !C
  ! C      -------  INTERNAL SUBROUTINE SUB1  -------
     200 IF( RAM3-RAM2 .GT. 5.0D0*(RAM2-RAM1) )  GO TO 202
         IF( RAM2-RAM1 .GT. 5.0D0*(RAM3-RAM2) )  GO TO 204
         A1 = (RAM3-RAM2)*E1
         A2 = (RAM1-RAM3)*E2
         A3 = (RAM2-RAM1)*E3
         B2 = (A1+A2+A3)*2.D0
         B1 = A1*(RAM3+RAM2) + A2*(RAM1+RAM3) + A3*(RAM2+RAM1)
         IF( B2 .EQ. 0.D0 )  GO TO 210
         RAM = B1 /B2
   !C
         IF( RAM .LE. RAM1 )  RAM = (RAM1 + RAM2)/2.0D0
         IF( RAM .GE. RAM3 )  RAM = (RAM2 + RAM3)/2.0D0
         IF( DABS(RAM-RAM2) .LE. 1.0D-15 )  RAM = (RAM2*4.D0 + RAM3)/5.0D0
  ! cc      GO TO RETURN ,( 80,130 )
         IF( IRET.EQ.80 ) GO TO 80
         IF( IRET.EQ.130 ) GO TO 130
     202 RAM = (4.0D0*RAM2 + RAM3)/5.0D0
  ! cc      GO TO RETURN, (80,130)
         IF( IRET.EQ.80 ) GO TO 80
         IF( IRET.EQ.130 ) GO TO 130
     204 RAM = (RAM1 + 4.0D0*RAM2)/5.0D0
  ! cc      GO TO RETURN, (80,130)
         IF( IRET.EQ.80 ) GO TO 80
         IF( IRET.EQ.130 ) GO TO 130
  ! C
     210 IG = 1
         RAM = RAM2
         RETURN
  ! C
  ! C      -------  INTERNAL SUBROUTINE SUB2  -------
  ! C
     300 IF( RAM3-RAM2 .GT. RAM2-RAM1 )  GO TO 310
         RAM = (RAM1+RAM2)*0.5D0
   !cc      GO TO RETURN ,( 80,130 )
         IF( IRET.EQ.80 ) GO TO 80
         IF( IRET.EQ.130 ) GO TO 130
   
     310 RAM = (RAM2+RAM3)*0.5D0
   !cc      GO TO RETURN ,( 80,130 )
         IF( IRET.EQ.80 ) GO TO 80
         IF( IRET.EQ.130 ) GO TO 130
   !C
     400 RAM = 0.0_8
         RETURN
   !C
     500 RAM = (RAM2+RAM3)*0.5_8
   !cxx 19/12/07
            ire510 = 0
   !cxx
#if defined __ICC
            !DIR$ VECTOR ALIGNED
            !DIR$ VECTOR ALIGNED
#elif defined __GFORTRAN__
            !$OMP SIMD LINEAR(I:1)
#endif
     510 DO 520  I=1,K
   !cxx  520 X1(I) = X(I) + RAM*H(I)
         X1(I) = X(I) + RAM*H(I)
     520 CONTINUE
  ! cc      CALL  FUNCT( K,X1,E3,G,IG )
  ! cxx      CALL  FUNCT( K,X1,E3,G,IG, Y,N,M,L,MLMAX,
         CALL  FUNCT( K,X1,E3,IG, Y,N,M,L,MLMAX, &
                            OUTMIN,OUTMAX,ALIMIT,FLK,SIG2,IER )
         if( ier.ne.0 ) return
   !cc      IF( IPR.GE.7 )  WRITE(6,7)  RAM,E3
         IF( IG.EQ.1 )  GO TO 540
         IF( E3.GT.E2 )  GO TO 530
         RAM1 = RAM2
         RAM2 = RAM
         E1 = E2
         E2 = E3
         GO TO 500
   !C
     530 RAM3 = RAM
         GO TO 70
   !C
     540 RAM = (RAM2+RAM)*0.5_8
   !cxx      GO TO 510
   !cxx 19/12/07
            if (RAM .EQ. RAM2) return
            ire510 = ire510 + 1
            if (ire510 .le. 100) GO TO 510
#if 0
   C
   C ------------------------------------------------------------
   cxx    1 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E1 =',D25.17 )
   cxx    2 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E2 =',D25.17 )
   cxx    3 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E3 =',D25.17 )
   cxx    4 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E4 =',D25.17 )
   cxx    5 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E5 =',D25.17 )
   cxx    6 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E6 =',D25.17 )
   cxx    7 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E7 =',D25.17 )
#endif
         E N D
      
        
         SUBROUTINE BOXCOXF(Y, N, AICZT, FFZT, AICZ, FFZ, ZMEAN, ZVAR, ZZ) !GCC$ ATTRIBUTES HOT :: BOXCOXF !GCC$ ATTRIBUTES ALIGNED(32) :: BOXCOXF
           implicit none
#if 0
   C
   C  ...  Box-Cox transformation
   C
   C     The following inputs are required in READTS
   C        TITLE:   title of the data set
   C        N:       data length
   C        Y(I):    data
   C     12/27/90 Y.I. 8/7/91 G.K.
   C
   cc      PARAMETER(MJ=1000)
   cxx      IMPLICIT REAL*8(A-H,O-Z )
   cc      CHARACTER  TITLE*72
   cc      DIMENSION  Y(MJ), Z(MJ)
   cc      COMMON  /CMDATA/  TITLE
   cxx      DIMENSION Y(N), Z(N), ZZ(N)
   cxx      DIMENSION AICZT(21), FFZT(21), AICZ(21), FFZ(21)
   cxx      DIMENSION ZMEAN(21), ZVAR(21)
   C
#endif
         INTEGER :: N
         REAL(8) :: Y(N), AICZT(21), FFZT(21), AICZ(21), FFZ(21), &
                    ZMEAN(21), ZVAR(21), ZZ(N)

         REAL(8) :: Z(N), YMEAN, YVAR, FFY, AICY, A, ZJACOB, AICM
#if defined __ICC
         !DIR$ ATTRIBUTES ALIGN : 64 :: Z
#endif
  ! C
  ! cc      CALL  READTS( 1,Y,N )
         CALL  GAUSSM( Y,N,YMEAN,YVAR,FFY,AICY )
  ! cc      WRITE(6,600)
  ! cc      WRITE(6,610)  TITLE
  ! cc      WRITE(6,620)
         I = 0
         AICM=AICZT(1)
         DO 200 II=10,-10,-1
            I = I+1
         A = II/10.0_8
  ! cc      CALL  BOXCOX( Y,N,A,Z,ZJACOB )
         CALL  BOXCOX( Y,N,A,Z,ZJACOB )
   !C
  ! cc      CALL  GAUSSM( Z,N,ZMEAN,ZVAR,FFZ,AICZ )
         CALL  GAUSSM( Z,N,ZMEAN(I),ZVAR(I),FFZ(I),AICZ(I) )
  ! C
  ! cc      FFZT = FFZ + ZJACOB
  ! cc      AICZT = AICZ-2*ZJACOB
         FFZT(I) = FFZ(I) + ZJACOB
         AICZT(I) = AICZ(I)-2*ZJACOB
  ! cc      WRITE(6,630)  A, AICZT, FFZT, AICZ, FFZ, ZMEAN, ZVAR
  ! c-----
         IF( I.EQ.1 ) AICM=AICZT(1)
         IF( AICZT(I).LE.AICM ) THEN
            DO 100 J=1,N
  ! cxx  100    ZZ(J) = Z(J)
            ZZ(J) = Z(J)
     100    CONTINUE
            AICM = AICZT(I)
         END IF
  ! c-----
     200 CONTINUE
  ! cc      STOP
         RETURN
  ! cxx  600 FORMAT( 1H ,'PROGRAM 4.4:   BOX-COX TRANSFORMATION' )
  ! cxx  610 FORMAT( 1H ,A72 )
  ! cxx  620 FORMAT( 1H ,'LAMBDA',5X,'AIC''',8X,'LL''',7X,'AIC',9X,'LL',
  ! cxx     *            9X,'MEAN',9X,'VARIANCE' )
  ! cxx  630 FORMAT( 1H ,F5.2,4F11.2,2D15.6 )
         E N D

         SUBROUTINE BOXCOX( Y,N,A,Z,ZJACOB ) !GCC$ ATTRIBUTES HOT :: BOXCOX !GCC$ ATTRIBUTES ALIGNED(32) :: BOXCOX
#if defined __GFORTRAN__ 
                use omp_lib
#endif
           implicit none
#if 0           
   C
   C  ...  Box-Cox transformation:  Z = (Y**A - 1)/A  ...
   C
   C     Inputs:
   C        Y(I):   data
   C        N:      data length
   C        A:      lambda
   C     Outputs:
   C        Z(I):   transformed data
   C        ZJACOB: log(Jacobian), log |dz/dy|
   C
   cxx      IMPLICIT REAL*8(A-H,O-Z )
   cxx      DIMENSION  Y(N), Z(N)
   C
#endif
         INTEGER :: N
         REAL(8) :: Y(N), A, Z(N), ZJACOB 
#if defined __ICC
         !DIR$ ASSUME_ALIGNED Y:64,Z:64
#endif
         REAL(8) :: SUM
   
         SUM = 0.0_8
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
!$OMP SIMD REDUCTION(+:SUM)
#endif
         DO 10 I=1,N
         IF( A .NE. 0.0_8) THEN
            Z(I) = (Y(I)**A - 1)/A
            SUM = SUM + (A-1)*DLOG( DABS( Y(I) ) )
         ELSE
            Z(I) = DLOG( Y(I) )
            SUM = SUM - DLOG( DABS( Y(I) ) )
         END IF
      10 CONTINUE
         ZJACOB = SUM
   
         RETURN
         E N D


         SUBROUTINE GAUSSM( Y,N,YMEAN,YVAR,FF,AIC ) !GCC$ ATTRIBUTES HOT :: GAUSSM !GCC$ ATTRIBUTES ALIGNED(32) :: GAUSSM
#if defined __GFORTRAN__
               use omp_lib
#endif
            implicit none
#if 0
   C
   C  ...  This subroutine fits Gaussian distribution model  ...
   C
   C     Inputs:
   C        Y(I):   data
   C        N:      data length
   C     Outputs:
   C        YMEAN:  mean
   C        YVAR:   variance
   C        FF:     log-likelihood
   C        AIC:    AIC = -2FF + 4
   C
   cxx      IMPLICIT REAL*8( A-H,O-Z )
   cxx      DIMENSION Y(N)
   C
#endif
         INTEGER :: N
         REAL(8) :: Y(N), YMEAN, YVAR, FF, AIC 
#if defined __ICC
         !DIR$ ASSUME_ALIGNED Y:64
#endif
         REAL(8) :: PI, SUM
   
         DATA  PI/3.1415926535897932384626433_8/
   
         SUM = 0.0_8
#if defined __ICC
!DIR$ VECTOR ALWAYS
!DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
!$OMP SIMD REDUCTION(+:SUM) LINEAR(I:1)
#endif
         DO 10 I=1,N
   !cxx   10 SUM = SUM + Y(I)
         SUM = SUM + Y(I)
      10 CONTINUE
         YMEAN = SUM/N
   
         SUM = 0.0_8
#if defined __ICC
         !DIR$ VECTOR ALWAYS
         !DIR$ SIMD REDUCTION(+:SUM)
#elif defined __GFORTRAN__
         !$OMP SIMD REDUCTION(+:SUM) LINEAR(I:1)
#endif
         DO 20 I=1,N
   !cxx   20 SUM = SUM + (Y(I)-YMEAN)**2
         SUM = SUM + (Y(I)-YMEAN)**2
      20 CONTINUE
         YVAR = SUM/N
         FF  = -0.5_8*N*(DLOG( 2*PI*YVAR ) + 1)
         AIC = -2*FF + 4
   
         RETURN
         E N D

#if 0         
         C------------------------------------------------- 
         C     ARCOEF  ---  arfit, armaft, season, tvar
         C     FOUGER  ---  spgrh, tvar
         C     PARCOR  ---  arfit, armaft, season, tvar
         C     REGRES  ---  lsqr, arfit, lsar1, lsar2, polreg
         C     COMAIC  ---  lsqr, arfit, lsar1, lsar2, polreg
         C     HUSHLD  ---  lsqr, arfit, marlsq, lsar1, lsar2 
         C     RECOEF  ---  lsqr, arfit, lsar1, lsar2, polreg
         C     REDUCT  ---  lsqr, arfit, lsar1, lsar2, polreg
         C     SETXAR  ---  arfit, lsar1, lsar2
         C     ARYULE  ---  arma, arfit
         C     ARMASP  ---  arma, arfit, lsar1, tvspc
         C     FOURIE  ---  period, arma, marspc, arfit, tvspc
         C     MOMENT ---  trend, season, tvvar, ngsmth
         C     SMOOTH  ---  tvvar, smooth
         C     GINVRS  ---  tvvar, smooth, season
         C     SETSEA  ---  simssm. ngsim
         C     CHOLES  ---  simssm. ngsim
         C     CRSCOR  ---  crscor, marfit
         C     MEAN  ---  unicor, crscor, period, arfit, marfit
         C     AUTCOR  --- unicor, arfit
         C     AUTCOV  ---  unicor, period, arfit
         C     UNICOR  ---  unicor, arfit
         C     ERRACF  ---  unicor, arfit
         C     WINDOW  ---  period, fftper
         C     ARMCOV  ---  arma, armaft, season
         C     DECOM  ---  arma, armaft, season
         C     IMPULS  ---  arma, armaft, season
         C     SOLVE  ---  arma, armaft, season
         C     SMOTH1  ---  seasom, tvar
         C
         C     ID  ---  simssm, ngsim, season
         C     INIT  ---  simssm, ngsim
         C     GAUSS  ---  densty, klinfo, ngsmth, ngsim
         C     PEARSN  ---  densty, ngsmth, ngsim
         C     DBLEXP  ---  densty, ngsim, ngsmth
         C     CAUCHY  ---  densty, klinfo
         C-------------------------------------------------
         C
#endif
               SUBROUTINE  ARCOEF( PAR,K,A ) !GCC$ ATTRIBUTES HOT :: ARCOEF !GCC$ ATTRIBUTES ALIGNED(32) :: ARCOEF

                  implicit none
#if 0
         C
         C  ...  This subroutine computes AR coefficients from PARCOR  ...
         C
         C     Inputs:
         C        PAR(I):   PARCOR
         C        K:        Order of AR model
         C     Output:
         C        A(I):     AR coefficient
         C
         cxx      IMPLICIT  REAL*8(A-H,O-Z)
         cc      DIMENSION  PAR(K), A(K), AA(50)
         cxx      DIMENSION  PAR(K), A(K), AA(K)
#endif
               INTEGER :: K
               REAL(8) :: PAR(K), A(K), AA(K)
#if defined __ICC
               !DIR$ ASSUME_ALIGNED PAR:64,A:64,AA:64
#endif         
               DO 30  II=1,K
               A(II)  = PAR(II)
               AA(II) = PAR(II)
               IF( II-1.LE.0 )  GO TO 30

               DO 10  J=1,II-1
         !cxx   10 A(J) = AA(J) - PAR(II)*AA(II-J)
               A(J) = AA(J) - PAR(II)*AA(II-J)
            10 CONTINUE
               IF( II.LT.K )  THEN
                 DO 20  J=1,II-1
         !cxx   20   AA(J) = A(J)
                 AA(J) = A(J)
            20   CONTINUE
               END IF
            30 CONTINUE
               RETURN
               E N D
         
         
               SUBROUTINE FOUGER(G,LGP1,FC,FS,LF1) !GCC$ ATTRIBUTES ALIGNED(32) :: FOUGER !GCC$ ATTRIBUTES HOT :: FOUGER

                    implicit none
#if 0
         C
         C     FOURIER TRANSFORM (GOERTZEL METHOD)
         C     THIS SUBROUTINE COMPUTES FOURIER TRANSFORM OF G(I),I=0,1,...,LG AT
         C     FREQUENCIES K/(2*LF), K=0,1,...,LF AND RETURNS COSINE TRANSFORM IN
         C     FC(K) AND SINE TRANSFORM IN FS(K).
         C
         C       INPUTS:
         C          G(I):   ORIGINAL DATA (I=0,1,...,LG)
         C          LG1:    = LG+1
         C          LF1:    = LF+1
         C
         C       OUTPUTS:
         C          FC(I):  COSINE TRANSFORM OF G  (I=0,1,...,LF)
         C          FB(I):  SINE TRANSFORM OF G  (I=0,1,...,LF)
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cxx      DIMENSION G(LGP1),FC(LF1),FS(LF1)
#endif
               INTEGER :: LGP1, LF1
               REAL(8) :: G(LGP1),FC(LF1),FS(LF1)
#if defined __ICC
                !DIR$ ASSUME_ALIGNED G:64,FC:64,FS:64
#endif
               REAL(8) :: PI, ALF, T, AK, TK, CK, SK, CK2, UM0, UM1, UM2
               LG=LGP1-1
               LF=LF1-1
         !C     REVERSAL OF G(I),I=1,...,LGP1 INTO G(LG3-I)   LG3=LGP1+1
               IF(LGP1.LE.1) GO TO 110
               LG3=LGP1+1
               LG4=LGP1/2
               DO 100 I=1,LG4
               I2=LG3-I
               T=G(I)
               G(I)=G(I2)
         !cxx  100 G(I2)=T
               G(I2)=T
           100 CONTINUE
           110 PI=3.1415926536D0
               ALF=LF
               T=PI/ALF
               DO 10 K=1,LF1
               AK=K-1
               TK=T*AK
               CK=DCOS(TK)
               SK=DSIN(TK)
               CK2=CK+CK
               UM1=0.0_8
               UM2=0.0_8
               IF(LG.EQ.0) GO TO 12
               DO 11 I=1,LG
               UM0=CK2*UM1-UM2+G(I)
               UM2=UM1
         !cxx   11 UM1=UM0
               UM1=UM0
            11 CONTINUE
            12 FC(K)=CK*UM1-UM2+G(LGP1)
               FS(K)=SK*UM1
            10 CONTINUE
               RETURN
               END
        
               SUBROUTINE  PARCOR( A,K,PAR ) !GCC$ ATTRIBUTES HOT :: PARCOR !GCC$ ATTRIBUTES ALIGNED(32) :: PARCOR
#if defined __GFORTRAN__
                       use omp_lib
#endif
                    implicit none
#if 0
         C  ...  This subriutine computes PARCOR for AR coefficients  ...
         C
         C       Inputs:
         C          A:   Vector of AR coefficients
         C          K:   Order of the model
         C
         C       Output:
         C          PAR: Vector of partial autocorrelations (PARCOR)
         C
         cxx      IMPLICIT  REAL*8(A-H,O-Z)
         cc      DIMENSION  A(K), PAR(K), G(50)
         cxx      DIMENSION  A(K), PAR(K), G(K)
#endif
               INTEGER :: K
               REAL(8) :: A(K), PAR(K)
#if defined __ICC
              !DIR$ ASSUME_ALIGNED A:64,PAR:64
#endif
               REAL(8) :: G(K), S
#if defined __ICC
              !DIR$ ATTRIBUTES ALIGN : 64 :: G
#endif   
#if defined __ICC
!DIR$ VECTOR ALIGNED
!DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
!$OMP SIMD LINEAR(I:1)
#endif      
               DO 10  I=1,K
         !cxx   10 PAR(I) = A(I)
               PAR(I) = A(I)
            10 CONTINUE
        
               IF( K .EQ. 1 )   RETURN                                           
         
               DO 40  II=K-1,1,-1
               S = 1.0_8 - PAR(II+1)**2
#if defined __ICC
               !DIR$ VECTOR ALIGNED
               !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
               !$OMP SIMD 
#endif  
               DO 20  I=1,II
        ! cxx   20 G(I) = (PAR(I) + PAR(II+1)*PAR(II-I+1))/S
               G(I) = (PAR(I) + PAR(II+1)*PAR(II-I+1))/S
            20 CONTINUE
               I2 = (II+1)/2
               IF( MOD( II,2 ).EQ.1 )  G(I2) = PAR(I2)/(1.D0 - PAR(II+1))
#if defined __ICC
               !DIR$ VECTOR ALIGNED
               !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
               !$OMP SIMD LINEAR(I:1)
#endif  
               DO 30  I=1,II
        ! cxx   30 PAR(I) = G(I)
               PAR(I) = G(I)
            30 CONTINUE
            40 CONTINUE
         
               RETURN
               E N D
        
               SUBROUTINE  REGRES( X,K,N,MJ1,A,SIG2,AIC,IMIN )
         C
         C  ...  Regression model fitting  ...
         C  ...  Order of the model is selected by AIC  ...
         C
         C     Inputs:
         C        X:      Householder reduced form (upper triangular matrix)
         C        K:      Maximum number of regressors
         C        N:      Number of data
         C        MJ1:    Adjustable dimension of X
         C        MJ2:    Adjustable dimension of A, ISG2 and AIC
         C     Outputs:
         C        A(I,M): Regression coefficients of the model with order M
         C        SIG2:   Residual variances
         C        AIC:    AIC's
         C        IMIN:   MAICE order
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cc      DIMENSION  X(MJ1,1), A(MJ2,MJ2), SIG2(0:MJ2), AIC(0:MJ2)
         cxx      DIMENSION  X(MJ1,K+1), A(K,K), SIG2(0:K), AIC(0:K)
               INTEGER :: K, N, MJ1, IMIN
               REAL(8) :: X(MJ1,K+1), A(K,K), SIG2(0:K), AIC(0:K)
               REAL(8) :: AICM
         C
               A(1:K, 1:K) = 0.0D0
         C
               CALL  COMAIC( X,N,K,MJ1,SIG2,AIC )
         C
               IMIN = 0
               AICM = AIC(0)
               DO 10  M=1,K
               IF( AIC(M).LT.AICM )  THEN
                  IMIN = M
                  AICM = AIC(M)
               END IF
               CALL  RECOEF( X,M,K,MJ1,A(1,M) )
            10 CONTINUE
         C
               RETURN
               E N D
         C
         C
               SUBROUTINE  COMAIC( X,N,K,MJ1,SIG2,AIC )
         C
         C  ...  This subroutine computes residual variances and AIC's  ...
         C
         C     Inputs:
         C        X(I,J):  Householder reduced form
         C        N:       Data length
         C        K:       Highest order
         C        MJ1:     Adjustable dimension of X
         C     Outputs:
         C        SIG2(I): residual variance of the model with order I
         C        AIC(I):  AIC of the model with order I
         C
         cxx      IMPLICIT  REAL*8(A-H,O-Z)
         cc      DIMENSION  X(MJ1,1), AIC(0:K), SIG2(0:K)
         cxx      DIMENSION  X(MJ1,K+1), AIC(0:K), SIG2(0:K)
               INTEGER :: N, K, MJ1
               REAL(8) :: X(MJ1,K+1), SIG2(0:K), AIC(0:K)
               REAL(8) :: PI2, PVAR
               DATA  PI2/6.28318531D0/
         C
               PVAR = 0.0D0
               DO 10 I=K,0,-1
               PVAR = PVAR + X(I+1,K+1)**2
               SIG2(I) = PVAR / N
         cxx   10 AIC(I)  = N*DLOG( PI2*SIG2(I) ) + N + 2*(I+1)
               AIC(I)  = N*DLOG( PI2*SIG2(I) ) + N + 2*(I+1)
            10 CONTINUE
         C
               RETURN
               E N D
         C
         C
         cc      SUBROUTINE  HUSHLD( X,D,MJ1,N,K )
               SUBROUTINE  HUSHLD( X,MJ1,N,K )
         C
         C  ...  Householder transformation  ...
         C
         C     Inputs:
         C        X(I,J):  Original matrix
         C        D(I):    Working area
         C        MJ1:     Adjustable dimension of X
         C        N:       Number of rows of X
         C        K:       Number of columns of X
         C     Output:
         C        X(I,J):  Householder reduced form (upper triangular form)
         C
         cxx      IMPLICIT  REAL*8(A-H,O-Z)
         cc      DIMENSION  X(MJ1,1), D(MJ1)
         cxx      DIMENSION  X(MJ1,K), D(MJ1)
               INTEGER :: MJ1, N, K
               REAL(8) :: X(MJ1,K)
               REAL(8) :: D(MJ1), TOL, H, G, F, S
         C
               TOL = 1.0D-60
         C
         cxx      DO 100  II=1,K
               DO 101  II=1,K
                 H = 0.0D0
                 DO 10  I=II,N
                   D(I) = X(I,II)
         cxx   10     H = H + D(I)**2
                   H = H + D(I)**2
            10 CONTINUE
                 IF( H .GT. TOL )  GO TO 20
                 G = 0.0D0
                 GO TO 100
            20   G = DSQRT( H )
                 F = X(II,II)
                 IF( F .GE. 0.0D0 )   G = -G
                 D(II) = F - G
                 H = H - F*G
         C
                 DO 30  I=II+1,N
         cxx   30   X(I,II) = 0.0D0
                 X(I,II) = 0.0D0
            30   CONTINUE
                 DO 60  J=II+1,K
                   S = 0.0D0
                   DO 40  I=II,N
         cxx   40     S = S + D(I)*X(I,J)
                     S = S + D(I)*X(I,J)
            40     CONTINUE
                   S = S/H
                   DO 50  I=II,N
         cxx   50     X(I,J) = X(I,J) - D(I)*S
                     X(I,J) = X(I,J) - D(I)*S
            50     CONTINUE
            60   CONTINUE
           100 X(II,II) = G
           101 CONTINUE
         C
               RETURN
               E N D
         C
         C
               SUBROUTINE  RECOEF( X,M,K,MJ,A )
         C
         C  ...  Regression coefficients  ...
         C
         C     Inputs:
         C        X(I,J):  Householder reduced form
         C        M:       Number of actually used regressors
         C        K:       Heighest order
         C        MJ:      Adjustable dimension of X
         C     Output:
         C        A(I):    Vector of regression coefficients
         C
         cxx      IMPLICIT REAL*8 (A-H,O-Z)
         cc      DIMENSION  X(MJ,1), A(1)
         cxx      DIMENSION  X(MJ,K+1), A(M)
               INTEGER :: M, K, MJ
               REAL(8) :: X(MJ,K+1), A(M)
               REAL(8) :: SUM
         C
               A(M) = X(M,K+1)/X(M,M)
         c-----
               IF( M .EQ. 1 ) RETURN
         c-----
               DO 20 I=M-1,1,-1
               SUM = X(I,K+1)
               DO 10 J=I+1,M
         cxx   10 SUM  = SUM - A(J)*X(I,J)
               SUM  = SUM - A(J)*X(I,J)
            10 CONTINUE
         cxx   20 A(I) = SUM/X(I,I)
               A(I) = SUM/X(I,I)
            20 CONTINUE
         C
               RETURN
               E N D
         C
         C
         cc      SUBROUTINE  REDUCT( SETX,Z,D,NMK,N0,K,MJ1,X )
               SUBROUTINE  REDUCT1( SETX,Z,NMK,N0,K,MJ1,X )
         C
         C  ...  Successive Householder reduction  ...
         C
         C     Inputs:
         C        SETX:    Name of the subroutine for making X(I,J)
         C        Z(I):    Data vector
         C        D(I):    Working area
         C        NMK:     Number of actually used observations
         C        N0:      Time point of the previous set ofobservations
         C        K:       Heighest order of the model
         C        MJ1:     Adjustable dimension of X
         C     Output:
         C        X(I,J):  data matrix
         C
         cxx      IMPLICIT  REAL*8( A-H,O-Z )
         cc      DIMENSION  X(MJ1,1) , D(1), Z(1)
         cx      DIMENSION  X(MJ1,1) , Z(1)
         cxx      DIMENSION  X(MJ1,K+1) , Z(N0+NMK)
               INTEGER :: NMK, N0, K, MJ1
               REAL(8) :: Z(N0+NMK), X(MJ1,K+1)
         C
               L = MIN0( NMK,MJ1 )
               K1 = K + 1
               N1 = L
         C
               CALL  SETX( Z,N0,L,K,MJ1,0,X )
         cc      CALL  HUSHLD( X,D,MJ1,L,K1 )
               CALL  HUSHLD( X,MJ1,L,K1 )
               IF( N1 .GE. NMK )  RETURN
         C
            10 L = MIN0( NMK-N1,MJ1-K1 )
         C
               LK = L + K1
               N2 = N0 + N1
               CALL  SETX( Z,N2,L,K,MJ1,1,X )
         cc      CALL  HUSHLD( X,D,MJ1,LK,K1 )
               CALL  HUSHLD( X,MJ1,LK,K1 )
               N1 = N1 + L
               IF( N1.LT.NMK )  GO TO 10
         C
               RETURN
         C
               E N D
         C
         C
               SUBROUTINE  REDUCT( SETX,Z,NMK,N0,K,MJ1,X )
         C
         C  ...  Successive Householder reduction  ...
         C
         C     Inputs:
         C        SETX:    Name of the subroutine for making X(I,J)
         C        Z(I):    Data vector
         C        D(I):    Working area
         C        NMK:     Number of actually used observations
         C        N0:      Time point of the previous set ofobservations
         C        K:       Heighest order of the model
         C        MJ1:     Adjustable dimension of X
         C     Output:
         C        X(I,J):  data matrix
         C
         cxx      IMPLICIT  REAL*8( A-H,O-Z )
         cc      DIMENSION  X(MJ1,1) , D(1), Z(1)
         cx      DIMENSION  X(MJ1,1) , Z(1)
         cxx      DIMENSION  X(MJ1,K+1) , Z(N0+NMK+K)
               INTEGER :: NMK, N0, K, MJ1
               REAL(8) :: Z(N0+NMK+K), X(MJ1,K+1)
         C
               L = MIN0( NMK,MJ1 )
               K1 = K + 1
               N1 = L
         C
         cdd      CALL  SETX( Z,N0,L,K,MJ1,0,X )
               CALL SETXAR( Z,N0,L,K,MJ1,0,X )
         cc      CALL  HUSHLD( X,D,MJ1,L,K1 )
               CALL  HUSHLD( X,MJ1,L,K1 )
               IF( N1 .GE. NMK )  RETURN
         C
            10 L = MIN0( NMK-N1,MJ1-K1 )
         C
               LK = L + K1
               N2 = N0 + N1
               CALL  SETX( Z,N2,L,K,MJ1,1,X )
         cc      CALL  HUSHLD( X,D,MJ1,LK,K1 )
               CALL  HUSHLD( X,MJ1,LK,K1 )
               N1 = N1 + L
               IF( N1.LT.NMK )  GO TO 10
         C
               RETURN
         C
               E N D
         C
         C
               SUBROUTINE  SETXAR( Z,N0,L,K,MJ1,JSW,X )
         C
         C  ...  Data matrix for AR model  ...
         C
         C     Inputs:
         C        Z(I):    Time series
         C        N0:      Origin of the current observations
         C        L:       Number of current observations
         C        K:       Number of regressors
         C        MJ1:     Adjustable dimension of X
         C        JSW=0:   Make initial data matrix
         C           =1:   Apend L*(K+1) data matrix below the triangular one
         C     Output:
         C        X(I,J):  Data matrix
         C
         cc      REAL*8  X(MJ1,1), Z(1)
         cxx      REAL*8  X(MJ1,K+1), Z(N0+L+K)
               INTEGER :: N0, L, K, MJ1, JSW
               REAL(8) :: Z(N0+L+K), X(MJ1,K+1)
         C
               I0 = 0
               IF( JSW .EQ. 1 )     I0 = K+1
         cxx      DO 10  I=1,L
               DO 11  I=1,L
                 II = I + I0
                 JJ = N0 + K + I
                 X(II,K+1) = Z(JJ)
               DO 10  J=1,K
                 JJ = JJ - 1
         cxx   10 X(II,J) = Z(JJ)
               X(II,J) = Z(JJ)
            10 CONTINUE
            11 CONTINUE
         C
               RETURN
               E N D
         C
         C
               SUBROUTINE ARYULE( C,N,MAXM,SIG2,AIC,PARCOR,A,MAR )
         C
         C  ...  Yule-Walker method  ...
         C
         C     Inputs:
         C        C(I):    Autocovariance function
         C        N:       Data length
         C        MAXM:    Highest AR order
         C     Outputs:
         C        SIG2(I): Innovation variance
         C        AIC(I):  AIC
         C        PARCOR(I):  PARCOR
         C        AMIN:     AR coefficients of the best model
         C        MAR:      Selected order of the model
         C
         cxx      IMPLICIT REAL*8 (A-H,O-Z)
         cxx      DIMENSION  C(0:MAXM), SIG2(0:MAXM), AIC(0:MAXM)
         cxx      DIMENSION  PARCOR(MAXM), A(MAXM,MAXM)
               INTEGER :: N, MAXM, MAR
               REAL(8) :: C(0:MAXM), SIG2(0:MAXM), AIC(0:MAXM), PARCOR(MAXM),
              1           A(MAXM,MAXM)
               REAL(8) :: CONST, AICMIN, SUM
         C
               CONST = N*(DLOG(2*3.1415926535D0) + 1)
         C
               SIG2(0) = C(0)
               AIC(0) = CONST + N*DLOG(SIG2(0)) + 2
               AICMIN = AIC(0)
               MAR = 0
         C
               DO 50 M=1,MAXM
               SUM  = C(M)
               IF (M.GE.2) THEN
               DO 10 I=1,M-1
         cxx   10 SUM  = SUM - A(I,M-1)*C(M-I)
               SUM  = SUM - A(I,M-1)*C(M-I)
            10 CONTINUE
               END IF
               A(M,M) = SUM /SIG2(M-1)
               IF (M.GE.2) THEN
               DO 20 J=1,M-1
         cxx   20 A(J,M) = A(J,M-1)-A(M,M)*A(M-J,M-1)
               A(J,M) = A(J,M-1)-A(M,M)*A(M-J,M-1)
            20 CONTINUE
               END IF
               SIG2(M) = SIG2(M-1)*(1.0D0-A(M,M)**2)
               AIC(M) = CONST + N*DLOG(SIG2(M)) + 2*(M+1)
               PARCOR(M) = A(M,M)
               IF( AIC(M).LT.AICMIN )  THEN
                  AICMIN = AIC(M)
                  MAR = M
               END IF
            50 CONTINUE
               RETURN
               E N D
         C
         C
               SUBROUTINE ARMASP( A,M,B,L,SIG2,NF,SP )
         C
         C  ...  Logarithm of the power spectrum of the ARMA model  ...
         C
         C     Inputs:
         C        M:     AR order
         C        L:     MA order
         C        A(I):  AR coefficient
         C        B(I):  MA coefficient
         C        SIG2:  Innovation variance
         C        NF:    Number of frequencies
         C     Output:
         C        SP(I): Power spectrum (in log scale)
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cxx      DIMENSION A(M), B(L)
         cc      DIMENSION SP(0:NF), H(0:500), FR(0:500), FI(0:500)
         cxx      DIMENSION SP(0:NF), H(0:M+L), FR(0:NF), FI(0:NF)
               INTEGER :: M, L, NF
               REAL(8) :: A(M), B(L), SIG2, SP(0:NF)
               REAL(8) :: H(0:M+L), FR(0:NF), FI(0:NF)
         C
               H(0) = 1.0D0
               DO 10 I=1,M
         cxx   10 H(I) = -A(I)
               H(I) = -A(I)
            10 CONTINUE
         C
               CALL  FOURIE( H,M+1,NF+1,FR,FI )
         C
               DO 20 I=0,NF
         cxx   20 SP(I) = SIG2/( FR(I)**2 + FI(I)**2 )
               SP(I) = SIG2/( FR(I)**2 + FI(I)**2 )
            20 CONTINUE
         C
               IF (L .EQ. 0) GO TO 41
               H(0) = 1.0D0
               DO 30 I=1,L
         cxx   30 H(I) = -B(I)
               H(I) = -B(I)
            30 CONTINUE
               CALL  FOURIE( H,L+1,NF+1,FR,FI )
               DO 40 I=0,NF
         cxx   40 SP(I) = SP(I)*( FR(I)**2 + FI(I)**2 )
               SP(I) = SP(I)*( FR(I)**2 + FI(I)**2 )
            40 CONTINUE
            41 CONTINUE
         C
               DO 50 I=0,NF
         cxx   50 SP(I) = DLOG10( SP(I) )
               SP(I) = DLOG10( SP(I) )
            50 CONTINUE
         C
               RETURN
               E N D
         C
         C
               SUBROUTINE FOURIE( X,N,M,FC,FS )
         C
         C  ...  Discrete Fourier transformation by Goertzel method  ...
         C
         C     Inputs:
         C        X(I):   data (I=1,N)
         C        N:      data length
         C        M:      number of Fourier components
         C        FC(J):  Fourier cosine transform (J=1,M)
         C        FS(J):  Fourier sine transform   (J=1,M)
         C
         cxx      IMPLICIT REAL*8 (A-H,O-Z)
         cxx      DIMENSION  X(N), FC(M), FS(M)
               INTEGER :: N, M
               REAL(8) :: X(N), FC(M), FS(M)
               REAL(8) :: PI, W, CI, SI, T0, T1, T2
               DATA  PI/3.14159265358979D0/
         C
               W = PI/(M-1)
               DO 20 I=1,M
               CI = DCOS(W*(I-1))
               SI = DSIN(W*(I-1))
               T1 = 0.0
               T2 = 0.0
               DO 10 J=N,2,-1
                 T0 = 2*CI*T1 - T2 + X(J)
                 T2 = T1
         cxx   10   T1 = T0
               T1 = T0
            10 CONTINUE
               FC(I) = CI*T1 - T2 + X(1)
         cxx   20 FS(I) = SI*T1
               FS(I) = SI*T1
            20 CONTINUE
         C
               RETURN
               E N D
               SUBROUTINE  MOMENT( Y,N,YMEAN,VAR )
         C
         C  ...  Mean and variance of the data  ...
         C
         C     Inputs:
         C       Y(I):   data
         C       N:      data length
         C     Outputs:
         C       YMEAN:  mean
         C       VAR:    variance
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cxx      DIMENSION  Y(N)
               INTEGER :: N
               REAL(8) :: Y(N), YMEAN, VAR, SUM
         C
               SUM = 0.0D0
               DO 10 I=1,N
         cxx   10 SUM = SUM + Y(I)
               SUM = SUM + Y(I)
            10 CONTINUE
               YMEAN = SUM/N
               SUM = 0.0D0
               DO 20 I=1,N
         cxx   20 SUM = SUM + (Y(I)-YMEAN)**2
               SUM = SUM + (Y(I)-YMEAN)**2
            20 CONTINUE
               VAR = SUM/N
         C
               RETURN
               E N D
         cc      SUBROUTINE  SMOOTH( F,M,MJ,NDIM,NS,NFE,NPE,VFS,VPS,XFS,XPS,
         cc     *                    VSS,XSS )
               SUBROUTINE  SMOOTH( F,M,NDIM,NS,NFE,NPE,VFS,VPS,XFS,XPS,VSS,XSS )
         C
         C  ...  Fixed Interval Smoother (General Form)  ...
         C
         C     Inputs:
         C        NS:     Start position of filtering
         C        NFE:    End position of filtering
         C        NPE:    End position of prediction
         C        M:      Dimension of the state vector
         C        F:      M*M matrix
         C        MJ:     Adjustable dimension of XF, VF
         C        NDIM:   Adjustable dimension of XFS, XPS, VFS, VPS
         C        NMAX    Adjustable dimension of Y
         C        VFS:    Covariance matrices of the filter
         C        VPS:    Covariance matrices of the predictor
         C        XFS:    Mean vectors of the filter
         C        XPS:    Mean vectors of the predictor
         C     Outputs:
         C        VSS:    Covariance matrices of the smoother
         C        XSS:    Mean vectors of the smoother
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cc      DIMENSION  F(MJ,MJ)
         cc      DIMENSION  XS(40), VS(40,40), VP(40,40)
         cc      DIMENSION  XFS(MJ,NDIM), XPS(MJ,NDIM), XSS(MJ,NDIM)
         cc      DIMENSION  VFS(MJ,MJ,NDIM), VPS(MJ,MJ,NDIM), VSS(MJ,MJ,NDIM)
         cc      DIMENSION  WRK(40,40), SGAIN(40,40)
         cxx      DIMENSION  F(M,M)
         cxx      DIMENSION  XS(M), VS(M,M), VP(M,M)
         cxx      DIMENSION  XFS(M,NDIM), XPS(M,NDIM), XSS(M,NDIM)
         cxx      DIMENSION  VFS(M,M,NDIM), VPS(M,M,NDIM), VSS(M,M,NDIM)
         cxx      DIMENSION  WRK(M,M), SGAIN(M,M)
               INTEGER :: M, NDIM, NS, NFE, NPE
               REAL(8) :: F(M,M), VFS(M,M,NDIM), VPS(M,M,NDIM), XFS(M,NDIM),
              1           XPS(M,NDIM), VSS(M,M,NDIM), XSS(M,NDIM)
               REAL(8) :: XS(M), VS(M,M), VP(M,M), WRK(M,M), SGAIN(M,M), VDET,
              1           SUM
         C
         C
         C  ...  SMOOTHING  ...
         C
         cxx      DO 10 II=NFE,NPE
         cxx      DO 10  I=1,M
               DO 12 II=NFE,NPE
               DO 11  I=1,M
               XSS(I,II)   = XFS(I,II)
               DO 10  J=1,M
         cxx   10 VSS(I,J,II) = VFS(I,J,II)
               VSS(I,J,II) = VFS(I,J,II)
            10 CONTINUE
            11 CONTINUE
            12 CONTINUE
         cxx      DO 20  I=1,M
               DO 21  I=1,M
               XS(I)   = XFS(I,NFE)
               DO 20  J=1,M
         cxx   20 VS(I,J) = VFS(I,J,NFE)
               VS(I,J) = VFS(I,J,NFE)
            20 CONTINUE
            21 CONTINUE
         C
               DO 500  II=NFE-1,NS,-1
         C
               NZERO = 0
               DO 100 I=1,M
         cxx  100 IF( VFS(I,I,II).GT.1.0D-12 )  NZERO = NZERO + 
               IF( VFS(I,I,II).GT.1.0D-12 )  NZERO = NZERO + 1
           100 CONTINUE
         C
               IF( NZERO.EQ.0 )  THEN
         cxx         DO 110  I=1,M
                  DO 111  I=1,M
                  XS(I)     = XFS(I,II)
                  XSS(I,II) = XFS(I,II)
                  DO 110  J=1,M
                  VS(I,J)     = VFS(I,J,II)
         cxx  110    VSS(I,J,II) = VFS(I,J,II)
                  VSS(I,J,II) = VFS(I,J,II)
           110    CONTINUE
           111    CONTINUE
         C
               ELSE
         cxx      DO 410  I=1,M
               DO 411  I=1,M
               DO 410  J=1,M
         cxx  410 VP(I,J) = VPS(I,J,II+1)
               VP(I,J) = VPS(I,J,II+1)
           410 CONTINUE
           411 CONTINUE
         C
         cc      CALL  GINVRS( VP,VDET,M,40 )
               CALL  GINVRS( VP,VDET,M )
         C
         cxx      DO 425  I=1,M
               DO 426  I=1,M
               DO 425  J=1,M
               SUM = 0.0D0
               DO 420  IJ=1,M
         cxx  420 SUM = SUM + VFS(I,IJ,II)*F(J,IJ)
               SUM = SUM + VFS(I,IJ,II)*F(J,IJ)
           420 CONTINUE
         cxx  425 WRK(I,J) = SUM
               WRK(I,J) = SUM
           425 CONTINUE
           426 CONTINUE
         C
         cxx      DO 440  I=1,M
               DO 441  I=1,M
               DO 440  J=1,M
               SUM = 0.0D0
               DO 430 IJ=1,M
         cxx  430 SUM = SUM + WRK(I,IJ)*VP(IJ,J)
               SUM = SUM + WRK(I,IJ)*VP(IJ,J)
           430 CONTINUE
         cxx  440 SGAIN(I,J) = SUM
               SGAIN(I,J) = SUM
           440 CONTINUE
           441 CONTINUE
         C
         cxx      DO 450  I=1,M
               DO 451  I=1,M
               XS(I) = XFS(I,II)
               DO 450  J=1,M
               WRK(I,J) = 0.0D0
         cxx  450 VS (I,J) = VFS(I,J,II)
               VS (I,J) = VFS(I,J,II)
           450 CONTINUE
           451 CONTINUE
         C
         cxx      DO 460  J=1,M
               DO 461  J=1,M
               DO 460  I=1,M
         cxx  460 XS(I) = XS(I) + SGAIN(I,J)*(XSS(J,II+1) - XPS(J,II+1))
               XS(I) = XS(I) + SGAIN(I,J)*(XSS(J,II+1) - XPS(J,II+1))
           460 CONTINUE
           461 CONTINUE
         C
         cxx      DO 470  J=1,M
         cxx      DO 470 IJ=1,M
               DO 472  J=1,M
               DO 471 IJ=1,M
               DO 470  I=1,M
         cxx  470 WRK(I,J)=WRK(I,J) + SGAIN(I,IJ)*(VSS(IJ,J,II+1)-VPS(IJ,J,II+1))
               WRK(I,J)=WRK(I,J) + SGAIN(I,IJ)*(VSS(IJ,J,II+1)-VPS(IJ,J,II+1))
           470 CONTINUE
           471 CONTINUE
           472 CONTINUE
         C
         cxx      DO 480  J=1,M
         cxx      DO 480 IJ=1,M
               DO 482  J=1,M
               DO 481 IJ=1,M
               DO 480  I=1,M
         cxx  480 VS(I,J) = VS(I,J) + WRK(I,IJ)*SGAIN(J,IJ)
               VS(I,J) = VS(I,J) + WRK(I,IJ)*SGAIN(J,IJ)
           480 CONTINUE
           481 CONTINUE
           482 CONTINUE
               DO 485 I=1,M
         cxx  485 IF( VS(I,I).LT.0.0D0 )  VS(I,I) = 0.0D0
               IF( VS(I,I).LT.0.0D0 )  VS(I,I) = 0.0D0
           485 CONTINUE
         C
         cxx      DO 490  I=1,M
               DO 491  I=1,M
               XSS(I,II) = XS(I)
               DO 490  J=1,M
         cxx  490 VSS(I,J,II) = VS(I,J)
               VSS(I,J,II) = VS(I,J)
           490 CONTINUE
           491 CONTINUE
               END IF
         C
           500 CONTINUE
         C
               RETURN
               E N D
         cc      SUBROUTINE  GINVRS( A,DET,M,MJ )
               SUBROUTINE  GINVRS( A,DET,M )
         C
         C  ...  Generalized inverse of a square matrix A  ...
         C
         C     Inputs:
         C        A:     M*M matrix
         C        M:     Dimension of A
         C        MJ:    Adjustable dimension of A
         C     Outputs:
         C        A:     Generalize inverse of A
         C        DET:   Determinant of A
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cc      DIMENSION  A(MJ,MJ), IND(50)
         cxx      DIMENSION  A(M,M), IND(M+1)
               INTEGER :: M, IND(M+1)
               REAL(8) :: A(M,M), DET, EPS, AMAX, SUM
         C
               EPS = 1.0D-10
         C
               DO 10  I=1,M
         cxx   10 IND(I) = I
                IND(I) = I
            10 CONTINUE
         cc------------
               I0 = 0
               LMAX = 0
         cc------------
               DO 60  L=1,M
               AMAX = 0.0D0
               DO 20  I=L,M
               IF( A(IND(I),IND(I)).GT.AMAX )  THEN
                 AMAX = A(IND(I),IND(I))
                 I0 = I
               END IF
            20 CONTINUE
               IF( AMAX.GT.EPS*A(IND(1),IND(1)) )  THEN
                  IMAX = IND(I0)
                  DO 30  I=I0,L+1,-1
         cxx   30    IND(I) = IND(I-1)
                  IND(I) = IND(I-1)
            30    CONTINUE
                  IND(L) = IMAX
                  LMAX   = L
         cxx         DO 40  I=L+1,M
                  DO 41  I=L+1,M
                  A(IND(I),IMAX) = -A(IND(I),IMAX)/A(IMAX,IMAX)
                  DO 40  J=L+1,M
         cxx   40    A(IND(I),IND(J)) = A(IND(I),IND(J))
                  A(IND(I),IND(J)) = A(IND(I),IND(J))
              *                    + A(IND(I),IMAX)*A(IMAX,IND(J))
            40    CONTINUE
            41    CONTINUE
               ELSE
         cxx         DO 50  I=L,M
                  DO 51  I=L,M
                  DO 50  J=L,M
         cxx   50    A(IND(I),IND(J)) = 0.0D0
                  A(IND(I),IND(J)) = 0.0D0
            50    CONTINUE
            51    CONTINUE
                  GO TO 70
               END IF
            60 CONTINUE
         C
            70 DET = 1.0D0
               DO 80  I=1,M
         C     DET = DET*A(IND(I),IND(I))
         C     IF( A(IND(I),IND(I)).GT.EPS*AMAX )  THEN
               IF( A(IND(I),IND(I)).GT.0.0D0 )  THEN
                 A(IND(I),IND(I)) = 1.0D0/A(IND(I),IND(I))
               ELSE
                 A(IND(I),IND(I)) = 0.0D0
               END IF
            80 CONTINUE
         C
               MS = MIN0( M-1,LMAX )
               DO 200  L=MS,1,-1
               DO 100 J=L+1,M
               SUM = 0.0D0
               DO 90  I=L+1,M
         cxx   90 SUM = SUM + A(IND(I),IND(L))*A(IND(I),IND(J))
               SUM = SUM + A(IND(I),IND(L))*A(IND(I),IND(J))
            90 CONTINUE
         cxx  100 A(IND(L),IND(J)) = SUM
               A(IND(L),IND(J)) = SUM
           100 CONTINUE
               SUM = A(IND(L),IND(L))
               DO 110  I=L+1,M
         cxx  110 SUM = SUM + A(IND(L),IND(I))*A(IND(I),IND(L))
               SUM = SUM + A(IND(L),IND(I))*A(IND(I),IND(L))
           110 CONTINUE
               A(IND(L),IND(L)) = SUM
               DO 120  I=L+1,M
         cxx  120 A(IND(I),IND(L)) = A(IND(L),IND(I))
               A(IND(I),IND(L)) = A(IND(L),IND(I))
           120 CONTINUE
               IMAX = IND(L)
               DO 130 I=L+1,M
               IF( IMAX.GT.IND(I) )  THEN
                  IND(I-1) = IND(I)
                  IND(I)   = IMAX
               END IF
           130 CONTINUE
           200 CONTINUE
         C
               RETURN
               E N D
         
         cc      SUBROUTINE  CRSCOR( Y,N,ID,LAG,MJ,MJ1,OUTMIN,OUTMAX,C,R,YMEAN )
                SUBROUTINE  CRSCOR( Y,N,ID,LAG,OUTMIN,OUTMAX,C,R,YMEAN )
         C
         C ... cross correlation function computation ...
         C
         C     Inputs:
         C        Y(I,J):   multi-variate time series
         C        N:        data length
         C        ID:       dimension of the observation
         C        LAG:      maximum lag of cross-covariance
         C        MJ:       adjustable dimension (MJ.GE.N)
         C        MJ1:      adjustable dimension (MJ1.GE.ID)
         C        OUTMIN:   bound for outliers in low side
         C        OUTMAX:   bound for outliers in high side
         C     Outputs:
         C        C:        sample cross-covariance function
         C        R:        sample cross-correlation function
         C        YMEAN:    sample mean vector
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cc      DIMENSION  Y(MJ,MJ1), OUTMIN(MJ1), OUTMAX(MJ1)
         cc      DIMENSION  C(0:LAG,MJ1,MJ1), R(0:LAG,MJ1,MJ1)
         cc      DIMENSION  YMEAN(MJ1), NSUM(10)
         cxx      DIMENSION  Y(N,ID), OUTMIN(ID), OUTMAX(ID)
         cxx      DIMENSION  C(0:LAG,ID,ID), R(0:LAG,ID,ID)
         cxx      DIMENSION  YMEAN(ID), NSUM(ID)
               INTEGER :: N, ID, LAG, NSUM(ID)
               REAL(8) :: Y(N,ID), OUTMIN(ID), OUTMAX(ID), C(0:LAG,ID,ID),
              1           R(0:LAG,ID,ID), YMEAN(ID), SUM
         C
               DO 10 J=1,ID
         cxx   10 CALL  MEAN( Y(1,J),N,OUTMIN(J),OUTMAX(J),NSUM(J),YMEAN(J) )
               CALL  MEAN( Y(1,J),N,OUTMIN(J),OUTMAX(J),NSUM(J),YMEAN(J) )
            10 CONTINUE
         C
         cxx      DO 30 I=1,ID
         cxx      DO 30 J=1,ID
               DO 32 I=1,ID
               DO 31 J=1,ID
         cc      WRITE(6,*)  I,J
               DO 30 L=0,LAG
               SUM = 0.0D0
               NNN = 0
               DO 20 II=L+1,N
               IF( Y(II,I).GT.OUTMIN(I).AND.Y(II,I).LT.OUTMAX(I) )  THEN
                 IF( Y(II-L,J).GT.OUTMIN(J).AND.Y(II-L,J).LT.OUTMAX(J) )  THEN
                   SUM = SUM + (Y(II,I)-YMEAN(I))*(Y(II-L,J)-YMEAN(J))
                   NNN = NNN + 1
                 END IF
               END IF
            20 CONTINUE
         cxx   30 C(L,I,J)=SUM/DSQRT( DBLE( NSUM(I)*NSUM(J) ) )
               C(L,I,J)=SUM/DSQRT( DBLE( NSUM(I)*NSUM(J) ) )
            30 CONTINUE
            31 CONTINUE
            32 CONTINUE
         C
         cxx      DO 40 I=1,ID
         cxx      DO 40 J=1,ID
               DO 42 I=1,ID
               DO 41 J=1,ID
               DO 40 L=0,LAG
         cxx   40 R(L,I,J) = C(L,I,J)/DSQRT(C(0,I,I)*C(0,J,J))
               R(L,I,J) = C(L,I,J)/DSQRT(C(0,I,I)*C(0,J,J))
            40 CONTINUE
            41 CONTINUE
            42 CONTINUE
               RETURN
         C
               E N D
         
               SUBROUTINE  MEAN( Y,N,OUTMIN,OUTMAX,NSUM,YMEAN )
         C
         C  ...  This subroutine computes sample mean  ...
         C
         C     Inputs:
         C        Y(I):    time series
         C        N:       data length
         C        OUTMIN:  bound for outliers in low side
         C        OUTMAX:  bound for outliers in high side
         C     Outputs:
         C        NSUM:    number of non-outlier observations
         C        YMEAN:   sample mean
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cxx      DIMENSION Y(N)
               INTEGER :: N, NSUM
               REAL(8) :: Y(N), OUTMIN, OUTMAX, YMEAN, YSUM
         C
               NSUM = 0
               YSUM = 0.0D0
               DO 10 I=1,N
               IF( Y(I) .GT.OUTMIN .AND. Y(I).LT.OUTMAX ) THEN
                  NSUM = NSUM + 1
                  YSUM = YSUM + Y(I)
               END IF
            10 CONTINUE
               YMEAN = YSUM/NSUM
         C
               RETURN
               E N D
               SUBROUTINE AUTCOR( C,MAXLAG,R )
         C
         C  ...  This subroutine compute sample autocorrelation  ...
         C
         C     Inputs:
         C        C(I):    sample autocovariance
         C        MAXLAG:  maximum lag of autocovariance
         C     Output:
         C        R(I):    sample autocorrelation
         C
         cxx      IMPLICIT REAL*8 (A-H,O-Z)
         cxx      DIMENSION  C(0:MAXLAG), R(0:MAXLAG)
               INTEGER :: MAXLAG
               REAL(8) :: C(0:MAXLAG), R(0:MAXLAG)
         C
               DO 10 LAG=0,MAXLAG
         cxx   10 R(LAG) = C(LAG)/C(0)
               R(LAG) = C(LAG)/C(0)
            10 CONTINUE
               RETURN
               E N D
         cc      SUBROUTINE AUTCOV( Y,N,MAXLAG,OUTMIN,OUTMAX,C )
               SUBROUTINE AUTCOV( Y,N,MAXLAG,OUTMIN,OUTMAX,C,YMEAN )
         C
         C  ...  This subroutine computes sample autocovariance  ...
         C
         C     Inputs:
         C        Y(I):    time series
         C        N:       data length
         C        MAXLAG:  maximum lag of autocovariance
         C        OUTMIN:  bound for outliers in low side
         C        OUTMAX:  bound for outliers in high side
         C     OUTPUT:
         C        C(I):    autocovariance function
         C
         cxx      IMPLICIT REAL*8( A-H,O-Z )
         cxx      DIMENSION Y(N), C(0:MAXLAG )
               INTEGER :: N, MAXLAG
               REAL(8) :: Y(N), OUTMIN, OUTMAX, C(0:MAXLAG ), YMEAN, SUM
         C
         C  ...  sample mean  ...
         C
               CALL  MEAN( Y,N,OUTMIN,OUTMAX,NSUM,YMEAN )
         C
         C  ...  sample autocovariance  ...
         C
               DO 20 LAG = 0,MAXLAG
               SUM = 0.0D0
               DO 10 I=LAG+1,N
               IF( Y(I).GT.OUTMIN .AND. Y(I).LT.OUTMAX )   THEN
                  IF( Y(I-LAG).GT.OUTMIN .AND. Y(I-LAG).LT.OUTMAX ) THEN
                  SUM = SUM + ( Y(I)-YMEAN)*( Y(I-LAG) - YMEAN )
                  END IF
               END IF
            10 CONTINUE
         cxx   20 C(LAG) = SUM/NSUM
               C(LAG) = SUM/NSUM
            20 CONTINUE
         C
               RETURN
               E N D
         cc      SUBROUTINE  UNICOR( Y,N,MAXLAG,OUTMIN,OUTMAX,COV )
               SUBROUTINE  UNICOR( Y,N,MAXLAG,OUTMIN,OUTMAX,COV,YMEAN )
         C
         C  ...  This subroutine computes sample autocovariance function,
         C       sample autocorrelationfunction and their error bounds  ...
         C
         C     Inputs:
         C        Y(I):    time series
         C        N:       data length
         C        MAXLAG:  maximum lag of autocovariance
         C        OUTMIN:  bound for outliers in low side
         C        OUTMAX:  bound for outliers in high side
         C     Output:
         C        COV(I,J):  I=0,...,MAXLAG
         C                   J = 1   sample autocovariance
         C                     = 2   sample autocorrelation
         C                     = 3   error bound for sample autocovariance
         C                     = 4   error bound for sample autocorrelation
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cxx      DIMENSION  Y(N), COV(0:MAXLAG,4)
               INTEGER :: N, MAXLAG
               REAL(8) :: Y(N), OUTMIN, OUTMAX, COV(0:MAXLAG,4), YMEAN
         C
         C  ...  sample autocovariance function  ...
         C
         cc      CALL  AUTCOV( Y,N,MAXLAG,OUTMIN,OUTMAX,COV )
               CALL  AUTCOV( Y,N,MAXLAG,OUTMIN,OUTMAX,COV,YMEAN )
         C
         C  ...  sample autocorrelation function  ...
         C
               CALL  AUTCOR( COV,MAXLAG,COV(0,2) )
         C
         C  ...  error bounds  ...
         C
               CALL  ERRACF( COV,N,MAXLAG,COV(0,3),COV(0,4) )
         C
               RETURN
               E N D
               SUBROUTINE ERRACF( C,N,MAXLAG,CERR,RERR )
         C
         C  ...  This subroutine computes error bounds for sample autocovariance
         C       and autocorrelation function  ...
         C
         C     Inputs:
         C        C(I):     sample autocovariance
         C        N:        data length
         C        MAXLAG:   maximum lag of autocovariance
         C     Outputs:
         C        CERR(I):  error bound for I-th autocovariance
         C        RERR(I):  error bound for I-th autocorrelation
         C
         cxx      IMPLICIT REAL*8 (A-H,O-Z)
         cxx      DIMENSION  C(0:MAXLAG), CERR(0:MAXLAG), RERR(0:MAXLAG)
               INTEGER :: N, MAXLAG
               REAL(8) :: C(0:MAXLAG), CERR(0:MAXLAG), RERR(0:MAXLAG), SUM
         C
               SUM = C(0)**2
               CERR(0)   = DSQRT( 2*SUM/N )
               DO 10 LAG=1,MAXLAG
               IF(LAG.GE.2)  SUM = SUM + 2*C(LAG-1)**2
         cxx   10 CERR(LAG) = DSQRT( SUM/N )
               CERR(LAG) = DSQRT( SUM/N )
            10 CONTINUE
         C
               RERR(0) = 0.0D0
               DO 20 LAG=1,MAXLAG
         cxx   20 RERR(LAG) = CERR(LAG)/C(0)
               RERR(LAG) = CERR(LAG)/C(0)
            20 CONTINUE
               RETURN
               E N D
         C
         C
         cc      SUBROUTINE  WINDOW( PE,NP,IWINDW,SPE )
               SUBROUTINE  WINDOW( PE,NP,IWINDW,SPE,IFG )
         C
         C  ...  Smoothing by spectral window and log-transformation  ...
         C
         C     Inputs:
         C        PE(I):    raw specrum
         C        NP:       number of frequencies
         C        IWINDW:   window type (0: box-car, 1: Hanning, 2: Hamming)
         C     Outputs:
         C        SPE(I):   logarithm of smoothed periodogram
         C
         cxx      IMPLICIT  REAL*8(A-H,O-Z)
         cxx      DIMENSION  PE(0:NP), SPE(0:NP)
         cxx      DIMENSION  W(0:1,2)
               INTEGER :: NP, IWINDW, IFG
               REAL(8) ::PE(0:NP), SPE(0:NP), W(0:1,2), PMIN
               DATA  W/0.5D0, 0.25D0, 0.54D0, 0.23D0/
         C
               IF( IWINDW.EQ.0 )  THEN
                 PMIN = 1.0D30
                 DO 10 I=0,NP
         cxx   10   IF( PE(I).GT.0.0D0 .AND. PE(I).LT.PMIN )  PMIN = PE(I)
                 IF( PE(I).GT.0.0D0 .AND. PE(I).LT.PMIN )  PMIN = PE(I)
            10   CONTINUE
                 DO 20 I=0,NP
         cxx   20   SPE(I) = DMAX1( PE(I),PMIN )
                 SPE(I) = DMAX1( PE(I),PMIN )
            20   CONTINUE
               ELSE
                 SPE(0) = W(0,IWINDW)*PE(0) + W(1,IWINDW)*PE(1)*2
                 SPE(NP)= W(0,IWINDW)*PE(NP)+ W(1,IWINDW)*PE(NP-1)*2
                 DO 30  I=1,NP-1
         cxx   30   SPE(I) = W(0,IWINDW)*PE(I) + W(1,IWINDW)*(PE(I-1) + PE(I+1))
                 SPE(I) = W(0,IWINDW)*PE(I) + W(1,IWINDW)*(PE(I-1) + PE(I+1))
            30   CONTINUE
               END IF
         c---------- 2013/07/03
               IF ( MINVAL(SPE) .GT. 0 ) THEN
               IFG = 0
               DO 50 I=0,NP
         c   50 SPE(I) = DLOG10( SPE(I) )
                SPE(I) = DLOG10( SPE(I) )
            50 CONTINUE
               ELSE
                  IFG = -1
               END IF
         C
               RETURN
               E N D
         C
         C
               SUBROUTINE  SETSEA( M1,M2,M3,IPER,AR,TAU1,TAU2,TAU3,SIG2,
         cc     *                    F,G,H,Q,R,M,K,MJ )
              *                    F,G,H,Q,R,M,K )
         C
         C  ...  SET STATE SPACE MODEL FOR SEASONAL ADJUSTMENT  ...
         C
         C     INPUTS:
         C        M1,M2,M3:   MODEL ORDERS
         C        IPER:       NUMBER OF SEASONS IN ONE PERIOD
         C        AR(I):      AR COEFFICIENTS
         C        TAU1-TAU3:  SYSTEM NOISE VARIANCES
         C        MJ:         ADJUSTABLE DIMENSION
         C     OUTPUTS:
         C        F,G,H:      MATRICES OF STATE SPACE MODEL
         C        Q,R:        SYSTEM NOISE AND OBS NOISE VARIANCES
         C        M:          STATE DIMENSION
         C        K:          DIMENSION OF SYSTEM NOISE
         C     MODIFIED  2/16/93
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cc      DIMENSION  F(MJ,MJ), G(MJ,MJ), H(MJ), Q(MJ,MJ), AR(*), R(1,1)
         cxx      DIMENSION  F(M,M), G(M,K), H(M), Q(K,K), AR(M3), R(1,1)
               INTEGER :: M1, M2, M3, IPER, M, K 
               REAL(8) :: AR(M3), TAU1, TAU2, TAU3, SIG2, F(M,M), G(M,K), H(M),
              1           Q(K,K), R(1,1)
         C
         cc      M = M1 + M2*(IPER-1) + M3
         cc      K = ID(M1) + ID(M2) + ID(M3)
         cxx      DO 10 I=1,M
         cxx      H(I) = 0.0D0
         cxx      DO 10 J=1,M
         cxx   10 F(I,J) = 0.0D0
         cxx      DO 20 I=1,M
         cxx      DO 20 J=1,K
         cxx   20 G(I,J) = 0.0D0
         cxx      DO 30 I=1,K
         cc   30 Q(I,I) = 0.0D0
         cxx      DO 30 J=1,K
         cxx   30 Q(I,J) = 0.0D0
               H(1:M) = 0.0D0
               F(1:M,1:M) = 0.0D0
               G(1:M,1:K) = 0.0D0
               Q(1:K,1:K) = 0.0D0
         C
               IF( M1.GT.0 )  THEN
                 IF( M1.EQ.1 )  F(1,1) = 1.0D0
                 IF( M1.EQ.2 )  THEN
                   F(1,1) = 2.0D0
                   F(1,2) =-1.0D0
                   F(2,1) = 1.0D0
                 END IF
                 G(1,1) = 1.0D0
                 H(1)   = 1.0D0
                 Q(1,1) = TAU1
               END IF
         C
         C  ...  SEASONAL COMPONENT  ...
         C
               IF( M2.GT.0 )  THEN
                 L1 = ID(M1) + 1
                 DO 40 I=1,IPER-1
         cxx   40   F(M1+1,M1+I) = -1.0D0
                 F(M1+1,M1+I) = -1.0D0
            40   CONTINUE
                 DO 50 I=2,IPER-1
         cxx   50   F(M1+I,M1+I-1) = 1.0D0
                 F(M1+I,M1+I-1) = 1.0D0
            50   CONTINUE
                 G(M1+1,L1) = 1.0D0
                 H(M1+1)    = 1.0D0
                 Q(L1,L1)   = TAU2
               END IF
         C
         C  ...  AR COMPONENT  ...
         C
               IF( M3.GT.0 )  THEN
                 M12 = M1 + M2*(IPER-1)
                 L12 = ID(M1) + ID(M2) + 1
                 DO 60 I=1,M3
         cxx   60   F(M12+1,M12+I) = AR(I)
                 F(M12+1,M12+I) = AR(I)
            60   CONTINUE
                 DO 70 I=2,M3
         cxx   70   F(M12+I,M12+I-1) = 1.0D0
                 F(M12+I,M12+I-1) = 1.0D0
            70   CONTINUE
                 G(M12+1,L12) = 1.0D0
                 H(M12+1)     = 1.0D0
                 Q(L12,L12)   = TAU3
               END IF
         C
               R(1,1) = SIG2
         C
               RETURN
               E N D
               SUBROUTINE  CHOLES ( X,MJ,N,Y,NJ )
         C
         C  ...  CHOLESKY DECOMPOSITION  ...
         C
         C     INPUTS:
         C        X(I,J):   SYMMETRIC MATRIX
         C        N:        DIMENSION OF Z
         C        MJ:       ADJUSTABLE DIMENSION OF X
         C        NJ:       ADJUSTABLE DIMENSION OF Y
         C     OUTPUT:
         C        Y(I,J):   LOWER TRIANGULAR MATRIX; X = Y*Y'
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cxx      DIMENSION  X(MJ,MJ), Y(NJ,NJ)
               INTEGER :: MJ, N, NJ
               REAL(8) :: X(MJ,MJ), Y(NJ,NJ), SUM1, SUM2
         C
         cxx      DO 10 I=1,N
         cxx      DO 10 J=1,N
         cxx   10 Y(I,J) = 0.0D0
               Y(1:N,1:N) = 0.0D0
         C
               DO 100 J = 1,N
               SUM1 = X(J,J)
               DO 20 K=1,J-1
         cxx   20 SUM1 = SUM1 - Y(J,K)**2
               SUM1 = SUM1 - Y(J,K)**2
            20 CONTINUE
               IF( SUM1.GT.0.0D0 ) Y(J,J) = DSQRT( SUM1 )
               IF( SUM1.EQ.0.0D0 ) Y(J,J) = 0.0D0
               DO 40 I=J+1,N
               SUM2 = 0.0D0
               DO 30 K=1,J-1
         cxx   30 SUM2 = SUM2 + Y(I,K)*Y(J,K)
               SUM2 = SUM2 + Y(I,K)*Y(J,K)
            30 CONTINUE
         cxx   40 Y(I,J) = ( X(I,J) - SUM2 ) / Y(J,J)
               Y(I,J) = ( X(I,J) - SUM2 ) / Y(J,J)
            40 CONTINUE
           100 CONTINUE
         C
               RETURN
               E N D
         cc      SUBROUTINE ARMCOV( M,L,A,B,SIG2,K,COV )
               SUBROUTINE ARMCOV( M,L,A,B,SIG2,K,COV,KMAX,IER )
         C
         C ...  Autocovariance Function of ARMA model  ...
         C
         C     Inputs:
         C        M:     AR order
         C        L:     MA order
         C        A(I):  AR coefficient
         C        B(I):  MA coefficient
         C        SIG2:  innovation variance
         C        K:     Required maximum lag of autocovariance
         C     Output:
         C        COV(I):  Autocovariance function
         C     Y.I.
         cxx      IMPLICIT REAL*8( A-H,O-Z )
         cc      DIMENSION  A(*), B(*), COV(0:K), G(0:100), X(30,30)
         cc      DIMENSION  Z(100), UL(30,30), IPS(100)
         cxx      DIMENSION  A(M), B(L), COV(0:K), G(0:KMAX), X(M+1,M+1)
         cxx      DIMENSION  Z(M+1), UL(M+1,M+1), IPS(M+1)
               INTEGER :: M, L, K, KMAX, IER, IPS(M+1)
               REAL(8) :: A(M), B(L), SIG2, COV(0:K)
               REAL(8) :: G(0:KMAX), X(M+1,M+1), Z(M+1), UL(M+1,M+1), SUM
         C
         cc      KMAX = MAX(M,L,K)
               CALL  IMPULS( M,L,A,B,KMAX,G )
         C
         cxx      DO 10 I=1,M+1
         cxx      DO 10 J=1,M+1
         cxx   10 X(I,J) = 0.0D0
               X(1:M+1,1:M+1) = 0.0D0
               DO 20 I=1,M+1
         cxx   20 X(I,I) = 1.0D0
               X(I,I) = 1.0D0
            20 CONTINUE
         c
               IF( M.GT.0 ) THEN
         
         cxx      DO 30 I=1,M
               DO 31 I=1,M
               DO 30 J=2,M-I+2
         cxx   30 X(I,J) = X(I,J) - A(I+J-2)
               X(I,J) = X(I,J) - A(I+J-2)
            30 CONTINUE
            31 CONTINUE
         cxx      DO 40 I=2,M+1
               DO 41 I=2,M+1
               DO 40 J=1,I-1
         cxx   40 X(I,J) = X(I,J) - A(I-J)
               X(I,J) = X(I,J) - A(I-J)
            40 CONTINUE
            41 CONTINUE
         
               END IF
         c
         C
         cc      CALL  DECOM( M+1,X,30,UL,IPS )
               CALL  DECOM( M+1,X,UL,IPS,IER )
               if( ier.ne.0 ) return
         C
               SUM = 1.0D0
         c
               IF( L.GT.0 ) THEN
         
               DO 50 J=1,L
         cxx   50 SUM = SUM - B(J)*G(J)
               SUM = SUM - B(J)*G(J)
            50 CONTINUE
         
               END IF
         c
               Z(1)= SIG2*SUM
         c
               IF( M.GT.0 ) THEN
         
               DO 70 I=2,M+1
               SUM = 0.0D0
               DO 60 J=I-1,L
         cc   60 SUM = SUM - B(J)*G(J-I+1)
         cxx   60 IF( L.GT. 0) SUM = SUM - B(J)*G(J-I+1)
               IF( L.GT. 0) SUM = SUM - B(J)*G(J-I+1)
            60 CONTINUE
         cxx   70 Z(I) = SIG2*SUM
               Z(I) = SIG2*SUM
            70 CONTINUE
         
               END IF
         c
         C
         cc      CALL  SOLVE( M+1,UL,30,Z,COV,IPS)
               CALL  SOLVE( M+1,UL,Z,COV,IPS )
         C
               DO 100 J=M+1,K
               SUM = 0.0D0
               DO 80 I=1,M
         cc   80 SUM = SUM + A(I)*COV(J-I)
         cxx   80 IF( M.GT.0 ) SUM = SUM + A(I)*COV(J-I)
               IF( M.GT.0 ) SUM = SUM + A(I)*COV(J-I)
            80 CONTINUE
               DO 90 I=J,L
         cc   90 SUM = SUM - B(I)*G(I-J)*SIG2
         cxx   90 IF( L.GT.0 ) SUM = SUM - B(I)*G(I-J)*SIG2
               IF( L.GT.0 ) SUM = SUM - B(I)*G(I-J)*SIG2
            90 CONTINUE
         cxx  100 COV(J) = SUM
               COV(J) = SUM
           100 CONTINUE
         C
               RETURN
               E N D
         C
         cc      SUBROUTINE DECOM( N,A,MJ,UL,IPS )
               SUBROUTINE DECOM( N,A,UL,IPS,IER )
         C
         C  ...  UL decomposition:  A = L*U  ...
         C
         C     Inputs:
         C        N:      Dimension of the matrix A
         C        A(I,J): N*N positive definite matrix
         C        MJ:     Adjustable dimension of A and UL
         C     Outputs:
         C        UL(I,J):  L-I and U
         C        IPS:    Index vector
         C     Y.I.
         C        IER:    Error code
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z )
         cc      DIMENSION A(MJ,*),UL(MJ,*),SCALES(100),IPS(100)
         cxx      DIMENSION A(N,N),UL(N,N),SCALES(N),IPS(N)
               INTEGER :: N, IPS(N), IER
               REAL(8) :: A(N,N), UL(N,N), SCALES(N), RNORM, BIG, SIZE, PIVOT, TM
         C
               IER = 0
         C
               DO 20 I=1,N
               IPS(I) = I
               RNORM = 0.0D0
               DO 10 J=1,N
               UL(I,J) = A(I,J)
               IF(RNORM.LT.ABS(UL(I,J)) ) RNORM = ABS( UL(I,J) )
            10 CONTINUE
               IF( RNORM .NE. 0.0D0 )  THEN
                   SCALES(I) = 1/RNORM
               ELSE
                   SCALES(I) = 0.0D0
         cc          CALL  SING(0)         
                   IER = 1
               END IF
            20 CONTINUE
               if( ier.ne.0 ) return
         C
         cc-------------
               INDEX = 0
         cc-------------
               DO 60 K=1,N-1
               BIG = 0.0D0
               DO 30 I=K,N
               SIZE = ABS( UL(IPS(I),K) )*SCALES( IPS(I) )
               IF( BIG.LT.SIZE ) THEN
                   BIG = SIZE
                   INDEX = I
               END IF
            30 CONTINUE
               IF( BIG.EQ. 0.0D0 )  THEN
         cc          CALL  SING(1)
                   IER = 2
                   GO TO 60
               END IF
               IF( INDEX.NE.K ) THEN
               J = IPS(K)
               IPS(K) = IPS(INDEX)
               IPS(INDEX) = J
               END IF
         C
               PIVOT = UL(IPS(K),K)
               DO 50 I=K+1,N
               TM = UL( IPS(I),K)/PIVOT
               UL( IPS(I),K) = TM
               IF( TM.NE. 0.0D0 )  THEN
               DO 40 J = K+1,N
         cxx   40 UL( IPS(I),J ) = UL( IPS(I),J)-TM*UL( IPS(K),J)
               UL( IPS(I),J ) = UL( IPS(I),J)-TM*UL( IPS(K),J)
            40 CONTINUE
         C     WRITE(6,*) (UL(IPS(I),J),J=1,N)
               END IF
            50 CONTINUE
            60 CONTINUE
               if( ier.ne.0 ) return
         C
         cc      IF( UL(IPS(N),N) .EQ. 0.0D0 )   CALL  SING(2)
               IF( UL(IPS(N),N) .EQ. 0.0D0 )   IER = 3
               RETURN
               E N D
         C
               SUBROUTINE IMPULS( M,L,A,B,K,G )
         C
         C ...  Impulse Response Function  ...
         C
         C     Inputs:
         C        M:     AR order
         C        L:     MA order
         C        A(I):  AR coefficient
         C        B(I):  MA coefficient
         C        K:     Required maximum lag of impulse respose
         C     Output:
         C        G(I):  Impulse response function
         C     Y.I.
         cxx      IMPLICIT REAL*8( A-H,O-Z )
         cc      DIMENSION A(*), B(*), G(0:K)
         cxx      DIMENSION A(M), B(L), G(0:K)
               INTEGER :: M, L, K
               REAL(8) :: A(M), B(L), G(0:K), SUM
         C
               G(0) = 1.0
               DO  20 I=1,K
               SUM = 0.0D0
               IF(I.LE.L) SUM = -B(I)
               DO  10 J=1,I
         cxx   10 IF(J.LE.M) SUM = SUM + A(J)*G(I-J)
               IF(J.LE.M) SUM = SUM + A(J)*G(I-J)
            10 CONTINUE
         cxx   20 G(I) = SUM
               G(I) = SUM
            20 CONTINUE
         C
               RETURN
               E N D
         C
         cc      SUBROUTINE  SOLVE( N,UL,MJ,B,X,IPS )
               SUBROUTINE  SOLVE( N,UL,B,X,IPS )
         C
         C  ...  Solve Ax=b using UL obtained by DECOM  ...
         C
         C     Inputs:
         C        N:     Dimension of UL and B
         C        UL:    LU decomposition of A
         C        MJ:    Adjustable dimension of A
         C        B:
         C        IPS:   index vector
         C     Output:
         C        X:     Solution
         C     Y.I.
         cxx      IMPLICIT REAL*8( A-H,O-Z )
         cc      DIMENSION UL(MJ,*),B(*),X(*),IPS(100)
         cxx      DIMENSION UL(N,N),B(N),X(N),IPS(N)
               INTEGER :: N, IPS(N)
               REAL(8) :: UL(N,N), B(N),X(N), SUM
         C
               DO 20 I=1,N
               SUM = 0.0D0
               DO 10 J=1,I-1
         cxx   10 SUM = SUM + UL(IPS(I),J)*X(J)
               SUM = SUM + UL(IPS(I),J)*X(J)
            10 CONTINUE
         cxx   20 X(I) = B(IPS(I)) - SUM
               X(I) = B(IPS(I)) - SUM
            20 CONTINUE
         C
               DO 40 I=N,1,-1
               SUM = 0.0D0
               DO 30 J=I+1,N
         cxx   30 SUM = SUM + UL(IPS(I),J)*X(J)
               SUM = SUM + UL(IPS(I),J)*X(J)
            30 CONTINUE
         cxx   40 X(I) = ( X(I)-SUM )/UL(IPS(I),I)
               X(I) = ( X(I)-SUM )/UL(IPS(I),I)
            40 CONTINUE
               RETURN
         cxx  600 FORMAT(1H ,'N=',I10,/,(5X,'IPS=',I10 ) )
               E N D
         C
         C
         cc      SUBROUTINE  SMOTH1( A,M,MMAX,NC,NS,N,NE,MJ,
               SUBROUTINE  SMOTH1( A,M,MMAX,NC,NS,N,NE,NMAX,MJ,
              *                    VFS,VPS,VSS,XFS,XPS,XSS )
         C
         C  ...  Fixed Interval Smoother (General Form)  ...
         C
         C     Inputs:
         C        A:      Parameter of matrix F
         C        M:      Dimension of the state vector
         C        MMAX:   Adjustable dimension of A
         C        NC:     Number of components
         C        N:      Data length
         C        NS:     Start position of filtering
         C        NE:     End position of filtering
         C        MJ:     Adjustable dimension of XF, VF
         C        VFS:    Covariance matrices of the filter
         C        VPS:    Covariance matrices of the predictor
         C        XFS:    Mean vectors of the filter
         C        XPS:    Mean vectors of the predictor
         C     Outputs:
         C        VSS:    Covariance matrices of the smoother
         C        XSS:    Mean vectors of the smoother
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cxx      DIMENSION  A(MMAX,NC), M(NC)
         cc      DIMENSION  XS(40), VS(40,40), VP(40,40)
         cc      DIMENSION  XFS(MJ,N), XPS(MJ,N), XSS(MJ,N)
         cc      DIMENSION  VFS(MJ,MJ,N), VPS(MJ,MJ,N), VSS(MJ,MJ,N)
         cc      DIMENSION  WRK(40,40), SGAIN(40,40)
         cc      DIMENSION  I0(10)
         cxx      DIMENSION  XS(MJ), VS(MJ,MJ), VP(MJ,MJ)
         cxx      DIMENSION  XFS(MJ,NMAX), XPS(MJ,NMAX), XSS(MJ,NMAX)
         cxx      DIMENSION  VFS(MJ,MJ,NMAX), VPS(MJ,MJ,NMAX), VSS(MJ,MJ,NMAX)
         cxx      DIMENSION  WRK(MJ,MJ), SGAIN(MJ,MJ)
         cxx      DIMENSION  I0(NC)
               INTEGER :: MMAX, NC, N S, N, NE, NMAX, MJ, M(NC), I0(NC)
               REAL(8) :: A(MMAX,NC), VFS(MJ,MJ,NMAX), VPS(MJ,MJ,NMAX),
              1           VSS(MJ,MJ,NMAX), XFS(MJ,NMAX), XPS(MJ,NMAX),
              2           XSS(MJ,NMAX)
               REAL(8) :: XS(MJ), VS(MJ,MJ), VP(MJ,MJ), WRK(MJ,MJ), SGAIN(MJ,MJ),
              1           SUM, VDET
         C
               I0(1) = 0
               DO 10 I=2,NC
         cxx   10 I0(I) = I0(I-1) + M(I-1)
               I0(I) = I0(I-1) + M(I-1)
            10 CONTINUE
               MM = I0(NC) + M(NC)
         C
         C  ...  SMOOTHING  ...
         C
               NSS = MIN0( N,NE )
         cxx      DO 20  I=1,MM
               DO 21  I=1,MM
               XS(I)      = XFS(I,NSS)
               XSS(I,NSS) = XFS(I,NSS)
               DO 20  J=1,MM
               VS(I,J)      = VFS(I,J,NSS)
         cxx   20 VSS(I,J,NSS) = VFS(I,J,NSS)
               VSS(I,J,NSS) = VFS(I,J,NSS)
            20 CONTINUE
            21 CONTINUE
         cxx      DO 30 II=NSS+1,NE
         cxx      DO 30 I=1,MM
               DO 32 II=NSS+1,NE
               DO 31 I=1,MM
               XSS(I,II) = XFS(I,II)
               DO 30 J=1,MM
         cxx   30 VSS(I,J,II) = VFS(I,J,II)
               VSS(I,J,II) = VFS(I,J,II)
            30 CONTINUE
            31 CONTINUE
            32 CONTINUE
         C
               DO 500  II=NSS-1,NS,-1
         C
               NZERO = 0
               DO 100 I=1,MM
         cxx  100 IF( VFS(I,I,II).GT.1.0D-12 )  NZERO = NZERO + 1
               IF( VFS(I,I,II).GT.1.0D-12 )  NZERO = NZERO + 1
           100 CONTINUE
         C
               IF( NZERO.EQ.0 )  THEN
         cxx         DO 110  I=1,MM
                  DO 111  I=1,MM
                  XS(I)     = XFS(I,II)
                  XSS(I,II) = XFS(I,II)
                  DO 110  J=1,MM
                  VS(I,J)     = VFS(I,J,II)
         cxx  110    VSS(I,J,II) = VFS(I,J,II)
                  VSS(I,J,II) = VFS(I,J,II)
           110 CONTINUE
           111 CONTINUE
         C
               ELSE
         cxx      DO 410  I=1,MM
               DO 411  I=1,MM
               DO 410  J=1,MM
         cxx  410 VP(I,J) = VPS(I,J,II+1)
               VP(I,J) = VPS(I,J,II+1)
           410 CONTINUE
           411 CONTINUE
         C
         cc      CALL  GINVRS( VP,VDET,MM,40 )
               CALL  GINVRS( VP,VDET,MM )
         C
         cxx      DO 420  I=1,MM
         cxx      DO 420  L=1,NC
               DO 422  I=1,MM
               DO 421  L=1,NC
               WRK(I,I0(L)+M(L)) = VFS(I,I0(L)+1,II)*A(M(L),L)
               DO 420  J=1,M(L)-1
         cxx  420 WRK(I,I0(L)+J) = VFS(I,I0(L)+1,II)*A(J,L) + VFS(I,I0(L)+J+1,II)
               WRK(I,I0(L)+J) = VFS(I,I0(L)+1,II)*A(J,L) + VFS(I,I0(L)+J+1,II)
           420 CONTINUE
           421 CONTINUE
           422 CONTINUE
         C
         cxx      DO 440  I=1,MM
               DO 441  I=1,MM
               DO 440  J=1,MM
               SUM = 0.0D0
               DO 430 IJ=1,MM
         cxx  430 SUM = SUM + WRK(I,IJ)*VP(IJ,J)
               SUM = SUM + WRK(I,IJ)*VP(IJ,J)
           430 CONTINUE
         cxx  440 SGAIN(I,J) = SUM
               SGAIN(I,J) = SUM
           440 CONTINUE
           441 CONTINUE
         C
         cxx      DO 450  I=1,MM
               DO 451  I=1,MM
               XS(I) = XFS(I,II)
               DO 450  J=1,MM
               WRK(I,J) = 0.0D0
         cxx  450 VS (I,J) = VFS(I,J,II)
               VS (I,J) = VFS(I,J,II)
           450 CONTINUE
           451 CONTINUE
         C
         cxx      DO 460  J=1,MM
               DO 461  J=1,MM
               DO 460  I=1,MM
         cxx  460 XS(I) = XS(I) + SGAIN(I,J)*(XSS(J,II+1) - XPS(J,II+1))
               XS(I) = XS(I) + SGAIN(I,J)*(XSS(J,II+1) - XPS(J,II+1))
           460 CONTINUE
           461 CONTINUE
         C
         cxx      DO 470  J=1,MM
         cxx      DO 470 IJ=1,MM
               DO 472  J=1,MM
               DO 471 IJ=1,MM
               DO 470  I=1,MM
         cxx  470 WRK(I,J) = WRK(I,J) + SGAIN(I,IJ)*(VSS(IJ,J,II+1)-VPS(IJ,J,II+1))
               WRK(I,J) = WRK(I,J) + SGAIN(I,IJ)*(VSS(IJ,J,II+1)-VPS(IJ,J,II+1))
           470 CONTINUE
           471 CONTINUE
           472 CONTINUE
         C
         cxx      DO 480  J=1,MM
         cxx      DO 480 IJ=1,MM
               DO 482  J=1,MM
               DO 481 IJ=1,MM
               DO 480  I=1,MM
         cxx  480 VS(I,J) = VS(I,J) + WRK(I,IJ)*SGAIN(J,IJ)
               VS(I,J) = VS(I,J) + WRK(I,IJ)*SGAIN(J,IJ)
           480 CONTINUE
           481 CONTINUE
           482 CONTINUE
               DO 485 I=1,MM
         cxx  485 IF( VS(I,I).LT.0.0D0 )  VS(I,I) = 0.0D0
               IF( VS(I,I).LT.0.0D0 )  VS(I,I) = 0.0D0
           485 CONTINUE
         C
         cxx      DO 490  I=1,MM
               DO 491  I=1,MM
               XSS(I,II) = XS(I)
               DO 490  J=1,MM
         cxx  490 VSS(I,J,II) = VS(I,J)
               VSS(I,J,II) = VS(I,J)
           490 CONTINUE
           491 CONTINUE
               END IF
         C
           500 CONTINUE
         C
               RETURN
               E N D
         C------------------------------------------------------------------
         C
               FUNCTION  ID( K )                                                 
         C                                                                       
         C  ...  ID = 1:    IF K > 0                                             
         C       ID = 0:    OTHERWISE                                            
         C                                                                       
               ID = 0                                                            
               IF( K .GT. 0 )  ID = 1                                            
               RETURN                                                            
               E N D                                                             
         
               SUBROUTINE INIT(IX)
               INTEGER :: IX, v(8)
         C
               if ( IX .ge. 0 ) then
                   call init_genrand64(IX)
               else
                  call date_and_time( values=v )
                  call init_genrand64( sum(v) )
               endif
               RETURN
               END
         
               DOUBLE PRECISION FUNCTION  GAUSS( X,PARAM )
         C
         C  ...  Gaussian (normal) distribution  ...
         C
         C     Inputs:
         C        X:
         C        PARAM(1):  mean
         C        PARAM(2):  variance
         C     Output:
         C        GAUSS:     density at X
         C
         cxx      IMPLICIT  REAL*8(A-H,O-Z)
         cxx      DIMENSION  PARAM(2)
               REAL(8) :: X, PARAM(2), C1
               DATA  C1  /2.506628275D0/
         C
               GAUSS = DEXP( -(X-PARAM(1))**2/(2*PARAM(2)) )/(C1*DSQRT(PARAM(2)))
               RETURN
               E N D
               DOUBLE PRECISION FUNCTION  PEARSN( X,PARAM )
         C
         C  ...  Pearson family of  distributions  ...
         C
         C     Inputs:
         C        X:
         C        PARAM(1):  location parameter, mu
         C        PARAM(2):  dispersion parameter, tau square
         C        PARAM(3):  shape parameter
         C     Output:
         C        PEARSN:    density at X
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cxx      DIMENSION  PARAM(3)
               REAL(8) :: X, PARAM(3), PI, dgammafn
               DATA  PI/3.1415926535D0/
         C
         CXX      PEARSN = DGAMMA(PARAM(3))/DGAMMA(PARAM(3)-0.5D0)
               PEARSN = dgammafn(PARAM(3))/dgammafn(PARAM(3)-0.5D0)
              *                  /DSQRT(PI)*PARAM(2)**(PARAM(3)-0.5D0)
              *                  /((X-PARAM(1))**2 + PARAM(2))**PARAM(3)
               RETURN
         C
               END
         
               DOUBLE PRECISION FUNCTION  DBLEXP( X,PARAM )
         C
         C  ...  double exponential distribution  f(x) = exp(x - exp(x))  ...
         C
         C     Inputs:
         C        X:
         C     Output:
         C        DBLEXP:    density at X
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cxx      dimension PARAM(3)
               REAL(8) :: X, PARAM(3)
         C
         cxx 2018/07/02      DBLEXP = DEXP( X-DEXP(X) )
               DBLEXP = DEXP( (X-PARAM(1))-DEXP(X-PARAM(1)) )
               RETURN
         C
               E N D
         
               DOUBLE PRECISION FUNCTION  CAUCHY( X,PARAM )
         C
         C  ...  Cauchy distribution  ...
         C
         C     Inputs:
         C        X:
         C        PARAM(1):  location parameter, mu
         C        PARAM(2):  dispersion parameter, tau square
         C     Output:
         C        CAUCHY:    density at X
         C
         cxx      IMPLICIT REAL*8(A-H,O-Z)
         cxx      DIMENSION  PARAM(2)
               REAL(8) :: X, PARAM(2), PI
               DATA  PI /3.1415926535D0/
         C
               CAUCHY = DSQRT( PARAM(2) )/(PARAM(2) + (X-PARAM(1))**2)/PI
               RETURN
         C
               E N D
         
   
