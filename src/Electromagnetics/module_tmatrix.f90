




module mod_tmatrix
      implicit none
      
 !include "mkl_lapack.fi"
 
   
     
   
     
    integer, parameter :: I32P = 4   ! returns  4
     
    integer, parameter :: R32P = 4 !selected_real_kind(6,37)
             
    integer, parameter :: R64P = 8 !selected_real_kind(15,307) 
    
    
    !=====================================================59
    !  File and module information:
    !  version,creation and build date, author,description
    !=====================================================59
    
    ! Major version
    !integer(I32P), parameter, public :: MOD_TMATRIX_MAJOR = 1_I32P
    
    ! Minor version
    !integer(I32P), parameter, public :: MOD_TMATRIX_MINOR = 0_I32P
    
    ! Micro version
   ! integer(I32P), parameter, public :: MOD_TMATRIX_MICRO = 0_I32P
    
    ! Module full version
   ! integer(I32P), parameter, public :: MOD_TMATRIX_FULLVER = 1000_I32P*MOD_TMATRIX_MAJOR + &
   !                                                           100_I32P*MOD_TMATRIX_MINOR  + &
   !                                                           10_I32P*MOD_TMATRIX_MICRO
    ! Module creation date
    !character(*),  parameter, public :: MOD_TMATRIX_CREATION_DATE = "27-05-2018 11:26 +00200 (SUN 27 MAY 2018 GMT+2)"
    
    ! Module build date   (should be set after successful compilation)
   ! character(*),  parameter, public :: MOD_TMATRIX_BUILD_DATE = " "
    
    ! Module author info
   ! character(*),  parameter, public :: MOD_TMATRIX_AUTHOR = "Programmer: Michael Mishchenko NASA, modified by: Bernard Gingold, e-mail: beniekg@gmail.com"      
                      
    ! Module short description
    !character(*),  parameter, public :: MOD_TMATRIX_DESCRIPT = " CALCULATION OF LIGHT SCATTERING BY POLYDISPERSE, RANDOMLY ORIENTED PARTICLES OF IDENTICAL AXIALLY SYMMETRIC SHAPE "
    
   
    
    
    integer(I32P), parameter, public :: NPN1  = 100_I32P
    
    integer(I32P), parameter, public :: NPNG1 = 300_I32P
    
    integer(I32P), parameter, public :: NPNG2 = 2_I32P*NPNG1
    
    integer(I32P), parameter, public :: NPN2  = 2_I32P*NPN1
    
    integer(I32P), parameter, public :: NPL   = NPN2+1_I32P
    
    integer(I32P), parameter, public :: NPN3  = NPN1+1_I32P
    
    integer(I32P), parameter, public :: NPN4  = 80_I32P
    
    integer(I32P), parameter, public :: NPN5  = 2_I32P*NPN4
    
    integer(I32P), parameter, public :: NPN6  = NPN4+1_I32P
    
    integer(I32P), parameter, public :: NPL1  = NPN5+1_I32P
    
    real(R64P),    parameter, private :: PI    = 3.1415926535897932384626433832795_R64P
    
   
    
   
    
    contains
    
    
    
    
    subroutine TMATRIX_Driver(RAT,NDISTR,AXMAX,NPNAX,B,GAM, &
                               NKMAX,EPS,NP,LAM,MRR,MRI,DDELT,NPNA,NDGS, &
                               CEXT,CSCA,W,COSPH,REFF,VEFF,ALPHA1,ALPHA2,ALPHA3,ALPHA4, &
                               BETA1,BETA2,SCATMAT)
          implicit none
          real(R64P),    intent(in)      :: RAT
          integer(I32P), intent(in)      :: NDISTR
          real(R64P),    intent(in)      :: AXMAX
          integer(I32P), intent(in)      :: NPNAX
          real(R64P),    intent(in)      :: B
          real(R64P),    intent(in)      :: GAM
          integer(I32P), intent(in)      :: NKMAX
          real(R64P),    intent(in)      :: EPS
          integer(I32P), intent(in)      :: NP
          real(R64P),    intent(in)      :: LAM
          real(R64P), intent(in)         :: MRR
          real(R64P), intent(in)         :: MRI
          real(R64P),    intent(inout)      :: DDELT
          integer(I32P), intent(in)      :: NPNA  ! Number of scattering angles
          integer(I32P), intent(in)      :: NDGS
          real(R64P),    intent(inout)   :: CEXT  ! Result
          real(R64P),    intent(inout)   :: CSCA  ! Result
          real(R64P),    intent(inout)   :: W     ! Result
          real(R64P),    intent(inout)   :: COSPH ! Result
          real(R64P),    intent(inout)   :: REFF  ! Result
          real(R64P),    intent(inout)   :: VEFF  ! Result
          real(R64P), dimension(NPL), intent(inout) :: ALPHA1
!DIR$     ASSUME_ALIGNED ALPHA1:64
          real(R64P), dimension(NPL), intent(inout) :: ALPHA2
!DIR$     ASSUME_ALIGNED ALPHA2:64
          real(R64P), dimension(NPL), intent(inout) :: ALPHA3
!DIR$     ASSUME_ALIGNED ALPHA3:64
          real(R64P), dimension(NPL), intent(inout) :: ALPHA4
!DIR$     ASSUME_ALIGNED ALPHA4:64
          real(R64P), dimension(NPL), intent(inout) :: BETA1
!DIR$     ASSUME_ALIGNED BETA1:64
          real(R64P), dimension(NPL), intent(inout) :: BETA2
!DIR$     ASSUME_ALIGNED BETA2:64
          real(R64P), dimension(6,NPNA), intent(inout) :: SCATMAT   ! contains elements of the normalized scattering matrix
!DIR$     ASSUME_ALIGNED SCATMAT:64
          ! Locals
          real(R64P), dimension(NPNG2) :: X,WA,S,SS,R,DR,DDR,DRR,DRI
!DIR$     ATTRIBUTES ALIGN : 64 :: X
!DIR$     ATTRIBUTES ALIGN : 64 :: WA
!DIR$     ATTRIBUTES ALIGN : 64 :: S
!DIR$     ATTRIBUTES ALIGN : 64 :: SS
!DIR$     ATTRIBUTES ALIGN : 64 :: R
!DIR$     ATTRIBUTES ALIGN : 64 :: DR
!DIR$     ATTRIBUTES ALIGN : 64 :: DDR
!DIR$     ATTRIBUTES ALIGN : 64 :: DRR
!DIR$     ATTRIBUTES ALIGN : 64 :: DRI
          real(R64P), dimension(NPN1) :: AN
!DIR$     ATTRIBUTES ALIGN : 64 :: AN
          real(R64P), dimension(NPN1,NPN1) :: ANN
!DIR$     ATTRIBUTES ALIGN : 64 ::  ANN
          real(R64P), dimension(1000)      :: XG,WG
!DIR$     ATTRIBUTES ALIGN : 64 ::  XG
!DIR$     ATTRIBUTES ALIGN : 64 ::  WG
          real(R64P), dimension(2000)      :: XG1,WG1
!DIR$     ATTRIBUTES ALIGN : 64 ::  XG1
!DIR$     ATTRIBUTES ALIGN : 64 ::  WG1
          real(R64P), dimension(NPN2,NPN2) :: TR1,TI1
          real(R64P), dimension(NPL)       :: AL1,AL2,AL3,AL4,BE1,BE2
!DIR$     ATTRIBUTES ALIGN : 64 ::  AL1
!DIR$     ATTRIBUTES ALIGN : 64 ::  AL2
!DIR$     ATTRIBUTES ALIGN : 64 ::  AL3
!DIR$     ATTRIBUTES ALIGN : 64 ::  AL4
!DIR$     ATTRIBUTES ALIGN : 64 ::  BE1
!DIR$     ATTRIBUTES ALIGN : 64 ::  BE2
          real(R64P), dimension(NPN6,NPN4,NPN4) :: RT11,RT12,RT21,RT22, &
                                                   IT11,IT12,IT21,IT22
          integer(I32P) :: NCHECK,IAX,I,NK,II,L1MAX,INK,NNM,N2,M,NMAX, &
                           NM,NM1,N11,NN2,IXXX,INM1,NMA,NGAUSS,MMAX, &
                           N,N1,NMIN,NMAX1,NNNGGG,NGAUS,NGGG,NN1,M1,N22,   &
                           L1M,LMAX,L1
          real(R64P)    :: DAX,R1,R2,AXI,Z1,Z2,Z3,CSCAT,CEXTIN,   &
                           A,XEV,QEXT1,QSCA1,QEXT,QSCA,TR1NN,DN1, &
                           TI1NN,TR1NN1,TI1NN1,DSCA,DEXT,DQSCA,   &
                           DQEXT,ZZ1,ZZ2,ZZ3,ZZ4,ZZ5,ZZ6,ZZ7,ZZ8, &
                           QSC,QXT,COEFF1,WGII,WGI,WALB,ASYMM,P,  &
                           PPI,PIR,PII
          
          COMMON /CT/ TR1,TI1
          COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
!!DIR$     ATTRIBUTES ALIGN : 64 :: /CT/
!!DIR$     ATTRIBUTES ALIGN : 64 :: /TMAT/
          !Exec code ...
          
          X   = 0._R64P
          WA  = 0._R64P 
          S   = 0._R64P 
          SS  = 0._R64P 
          R   = 0._R64P 
          DR  = 0._R64P 
          DDR = 0._R64P 
          DRR = 0._R64P 
          DRI = 0._R64P
          AN  = 0._R64P
          ANN = 0._R64P
          XG  = 0._R64P
          WG  = 0._R64P
          XG1 = 0._R64P
          WG1 = 0._R64P
          TR1 = 0._R64P
          TI1 = 0._R64P
          AL1 = 0._R64P
          AL2 = 0._R64P
          AL3 = 0._R64P
          AL4 = 0._R64P
          BE1 = 0._R64P
          BE2 = 0._R64P
          RT11 = 0._R64P
          RT12 = 0._R64P
          RT21 = 0._R64P
          RT22 = 0._R64P
          IT11 = 0._R64P
          IT12 = 0._R64P
          IT21 = 0._R64P
          IT22 = 0._R64P
          P = DACOS(-1._R64P)
          NCHECK = 0_I32P
          IF (NP.EQ.-1.OR.NP.EQ.-2) NCHECK=1
          IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK=1
          WRITE (6,5454) NCHECK
 5454     FORMAT ('NCHECK=',I1)
          DAX=AXMAX/NPNAX
          IF (DABS(RAT-1._R64P).GT.1.0E-8_R64P.AND.NP.EQ.-1) CALL SAREA (EPS,RAT)
          IF (DABS(RAT-1._R64P).GT.1.0E-8_R64P.AND.NP.GE.0) CALL SURFCH(NP,EPS,RAT)
          IF (DABS(RAT-1._R64P).GT.1.0E-8_R64P.AND.NP.EQ.-2) CALL SAREAC (EPS,RAT)
          WRITE(6,8000) RAT
8000      FORMAT ('RAT=',F8.6)
          IF(NP.EQ.-1.AND.EPS.GE.1._R64P) THEN
             WRITE(6,7000) EPS
7000         FORMAT('RANDOMLY ORIENTED OBLATE SPHEROIDS, A/B=',F11.7)  
          END IF
          IF(NP.EQ.-1.AND.EPS.LT.1._R64P) THEN 
              WRITE(6,7001) EPS
7001          FORMAT('RANDOMLY ORIENTED PROLATE SPHEROIDS, A/B=',F11.7)
          END IF
          IF(NP.GE.0) THEN 
             WRITE(6,7100) NP,EPS
7100         FORMAT('RANDOMLY ORIENTED CHEBYSHEV PARTICLES, T',  I1,'(',F5.2,')')
          END IF
          
         IF(NP.EQ.-2.AND.EPS.GE.1.0_R64P) THEN 
            WRITE(6, 7150) EPS
7150        FORMAT('RANDOMLY ORIENTED OBLATE CYLINDERS, D/L=',F11.7)
         END IF
         IF(NP.EQ.-2.AND.EPS.LT.1._R64P) THEN 
             WRITE(6, 7151) EPS
7151         FORMAT('RANDOMLY ORIENTED PROLATE CYLINDERS, D/L=',F11.7)
         END IF
         WRITE(6,7200) DDELT
7200     FORMAT ('ACCURACY OF COMPUTATIONS DDELT = ',D8.2)  
         WRITE(6,7400) LAM,MRR,MRI
7400     FORMAT('LAM=',F10.6,3X,'MRR=',D10.4,3X,'MRI=',D10.4)
         
        
       
! 
! 

 !
! 
! 
! 
          DDELT=0.1_R64P*DDELT
          DO 600 IAX=1,NPNAX
              AXI=AXMAX-DAX*DFLOAT(IAX-1)
              R1=0.89031_R64P*AXI
              R2=1.56538_R64P*AXI
              NK=IDINT(AXI*NKMAX/AXMAX+2)                       
              IF (NK.GT.1000) THEN 
                 WRITE(6,8001) NK
8001             FORMAT ('NK=',I4,' I.E., IS GREATER THAN 1000. ',  &
                       'EXECUTION TERMINATED.')
                 STOP
              END IF
             
              IF (NDISTR.EQ.3) CALL POWER (AXI,B,R1,R2)
!
              CALL GAUSS (NK,0,0,XG,WG)
              Z1=(R2-R1)*0.5_R64P
              Z2=(R1+R2)*0.5_R64P
              Z3=R1*0.5_R64P
              IF (NDISTR.EQ.5) GO TO 3
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
              DO I=1,NK
                 XG1(I)=Z1*XG(I)+Z2
                 WG1(I)=WG(I)*Z1
              ENDDO
              GO TO 4
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
    3         DO I=1,NK
                 XG1(I)=Z3*XG(I)+Z3
                 WG1(I)=WG(I)*Z3
              ENDDO
              DO I=NK+1,2*NK
                 II=I-NK
                 XG1(I)=Z1*XG(II)+Z2
                 WG1(I)=WG(II)*Z1
             ENDDO
             NK=NK*2
4            CALL DISTRB (NK,XG1,WG1,NDISTR,AXI,B,GAM,R1,R2,   &
                          REFF,VEFF,P)
             !PRINT 8002,R1,R2
 !8002    FORMAT('R1=',F10.6,'   R2=',F10.6)
 !        IF (DABS(RAT-1D0).LE.1D-6) PRINT 8003, REFF,VEFF
 !        IF (DABS(RAT-1D0).GT.1D-6) PRINT 8004, REFF,VEFF
 !8003    FORMAT('EQUAL-VOLUME-SPHERE REFF=',F8.4,'   VEFF=',F7.4)
 !8004    FORMAT('EQUAL-SURFACE-AREA-SPHERE REFF=',F8.4,
 !    &          '   VEFF=',F7.4)
 !        PRINT 7250,NK
 !7250    FORMAT('NUMBER OF GAUSSIAN QUADRATURE POINTS ',
 !    &          'IN SIZE AVERAGING =',I4)

             ! DO I=1,NPL
             !    ALPHA1(I)=0._R64P
             !    ALPHA2(I)=0._R64P
             !    ALPHA3(I)=0._R64P
              !   ALPHA4(I)=0._R64P
             !    BETA1(I)=0._R64P
             !    BETA2(I)=0._R64P
             !ENDDO      
             CSCAT=0._R64P
             CEXTIN=0._R64P
             L1MAX=0
             DO 500 INK=1,NK
                I=NK-INK+1
                A=RAT*XG1(I)
                XEV=2._R64P*P*A/LAM
                IXXX=XEV+4.05_R64P*XEV**0.333333_R64P
                INM1=MAX0(4,IXXX)
            IF (INM1.GE.NPN1) THEN 
                
                WRITE(6, 7333) NPN1
7333            FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,  &
                         '.  EXECUTION TERMINATED')
                STOP
            END IF
           
! 
            QEXT1=0._R64P
            QSCA1=0._R64P
            DO 50 NMA=INM1,NPN1
               NMAX=NMA
               MMAX=1
               NGAUSS=NMAX*NDGS
               IF (NGAUSS.GT.NPNG1) THEN 
                   WRITE(6,7340) NGAUSS
7340               FORMAT('NGAUSS =',I3,' I.E. IS GREATER THAN NPNG1.',  &
                    '  EXECUTION TERMINATED')
                   STOP
               END IF
               
             
 !7340         
! 7334          FORMAT(' NMAX =', I3,'  DC2=',D8.2,'   DC1=',D8.2)
! 7335          FORMAT('               NMAX1 =', I3,'  DC2=',D8.2,    &
 !                     '  DC1=',D8.2)
               CALL CONST(NGAUSS,NMAX,MMAX,P,X,WA,AN,ANN,S,SS,NP,EPS)
               CALL VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R, &
                       DR,DDR,DRR,DRI,NMAX)
               CALL TMATR0 (NGAUSS,X,WA,AN,ANN,S,SS,PPI,PIR,PII,R,DR,    &
                         DDR,DRR,DRI,NMAX,NCHECK)
               QEXT=0._R64P
               QSCA=0._R64P
               DO N=1_I32P,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2_I32P*N+1_I32P)
                  QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN+ &
                                 TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  QEXT=QEXT+(TR1NN+TR1NN1)*DN1
               ENDDO
               DSCA=DABS((QSCA1-QSCA)/QSCA)
               DEXT=DABS((QEXT1-QEXT)/QEXT)
!              PRINT 7334, NMAX,DSCA,DEXT
               QEXT1=QEXT
               QSCA1=QSCA
               NMIN=DFLOAT(NMAX)/2._R64P+1._R64P
               DO 10 N=NMIN,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  DQSCA=DN1*(TR1NN*TR1NN+TI1NN*TI1NN+ &
                             TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  DQEXT=(TR1NN+TR1NN1)*DN1
                  DQSCA=DABS(DQSCA/QSCA)
                  DQEXT=DABS(DQEXT/QEXT)
                  NMAX1=N
                  IF (DQSCA.LE.DDELT.AND.DQEXT.LE.DDELT) GO TO 12
   10          CONTINUE              
   12          CONTINUE
!c              PRINT 7335, NMAX1,DQSCA,DQEXT
               IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 55
               IF (NMA.EQ.NPN1) THEN 
                   WRITE(6,7336) NPN1
7336               FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,   &
           '.  EXECUTION TERMINATED') 
                   STOP
               END IF
                  
   50       CONTINUE
   55       NNNGGG=NGAUSS+1
!            IF (NGAUSS.EQ.NPNG1) PRINT 7336
            MMAX=NMAX1
!            DO 150 NGAUS=NNNGGG,NPNG1
               NGAUSS=NGAUS
               NGGG=2*NGAUSS
! 7336          FORMAT('WARNING: NGAUSS=NPNG1')
! 7337          FORMAT(' NG=',I3,'  DC2=',D8.2,'   DC1=',D8.2)
               CALL CONST(NGAUSS,NMAX,MMAX,P,X,WA,AN,ANN,S,SS,NP,EPS)
               CALL VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R,  &
                      DR,DDR,DRR,DRI,NMAX)
               CALL TMATR0 (NGAUSS,X,WA,AN,ANN,S,SS,PPI,PIR,PII,R,DR,      &
                          DDR,DRR,DRI,NMAX,NCHECK)
               QEXT=0._R64P
               QSCA=0._R64P
               DO 104 N=1,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN+ &
                               TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  QEXT=QEXT+(TR1NN+TR1NN1)*DN1
  104          CONTINUE
!               DSCA=DABS((QSCA1-QSCA)/QSCA)
!               DEXT=DABS((QEXT1-QEXT)/QEXT)
!c              PRINT 7337, NGGG,DSCA,DEXT
               QEXT1=QEXT
               QSCA1=QSCA
               IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 155
!               IF (NGAUS.EQ.NPNG1) PRINT 7336
 ! 150       CONTINUE
  155       CONTINUE
            QSCA=0._R64P
            QEXT=0._R64P
            NNM=NMAX*2

            DO 204 N=1,NNM
               QEXT=QEXT+TR1(N,N)
  204       CONTINUE
            IF (NMAX1.GT.NPN4) THEN 
                WRITE(6,7550) NMAX1
7550            FORMAT ('NMAX1 = ',I3, ', i.e. greater than NPN4.',  &
                   ' Execution terminated')
                STOP
            END IF
                      
            DO 213 N2=1,NMAX1
               NN2=N2+NMAX
               DO 213 N1=1,NMAX1
                  NN1=N1+NMAX
                  ZZ1=TR1(N1,N2)
                  RT11(1,N1,N2)=ZZ1
                  ZZ2=TI1(N1,N2)
                  IT11(1,N1,N2)=ZZ2
                  ZZ3=TR1(N1,NN2)
                  RT12(1,N1,N2)=ZZ3
                  ZZ4=TI1(N1,NN2)
                  IT12(1,N1,N2)=ZZ4
                  ZZ5=TR1(NN1,N2)
                  RT21(1,N1,N2)=ZZ5
                  ZZ6=TI1(NN1,N2)
                  IT21(1,N1,N2)=ZZ6
                  ZZ7=TR1(NN1,NN2)
                  RT22(1,N1,N2)=ZZ7
                  ZZ8=TI1(NN1,NN2)
                  IT22(1,N1,N2)=ZZ8
                  QSCA=QSCA+ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4+ &
                       ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8
  213       CONTINUE
!            PRINT 7800,0,DABS(QEXT),QSCA,NMAX
            DO 220 M=1,NMAX1
               CALL TMATR(M,NGAUSS,X,WA,AN,ANN,S,SS,PPI,PIR,PII,R,DR, &
                         DDR,DRR,DRI,NMAX,NCHECK)
               NM=NMAX-M+1
               NM1=NMAX1-M+1
               M1=M+1
               QSC=0._R64P
               DO 214 N2=1,NM1
                  NN2=N2+M-1
                  N22=N2+NM
                  DO 214 N1=1,NM1
                     NN1=N1+M-1
                     N11=N1+NM
                     ZZ1=TR1(N1,N2)
                     RT11(M1,NN1,NN2)=ZZ1
                     ZZ2=TI1(N1,N2)
                     IT11(M1,NN1,NN2)=ZZ2
                     ZZ3=TR1(N1,N22)
                     RT12(M1,NN1,NN2)=ZZ3
                     ZZ4=TI1(N1,N22)
                     IT12(M1,NN1,NN2)=ZZ4
                     ZZ5=TR1(N11,N2)
                     RT21(M1,NN1,NN2)=ZZ5
                     ZZ6=TI1(N11,N2)
                     IT21(M1,NN1,NN2)=ZZ6
                     ZZ7=TR1(N11,N22)
                     RT22(M1,NN1,NN2)=ZZ7
                     ZZ8=TI1(N11,N22)
                     IT22(M1,NN1,NN2)=ZZ8
                     QSC=QSC+(ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4 + &
                           ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8)*2._R64P
  214          CONTINUE
               NNM=2*NM
               QXT=0._R64P
               DO 215 N=1,NNM
                  QXT=QXT+TR1(N,N)*2._R64P
  215          CONTINUE
               QSCA=QSCA+QSC
               QEXT=QEXT+QXT
!               PRINT 7800,M,DABS(QXT),QSC,NMAX
! 7800          FORMAT(' m=',I3,'  qxt=',d12.6,'  qsc=',d12.6,
!     &                '  nmax=',I3)
  220       CONTINUE
            COEFF1=LAM*LAM*0.5_R64P/PI
            CSCA=QSCA*COEFF1
            CEXT=-QEXT*COEFF1
!c           PRINT 7880, NMAX,NMAX1
! 7880       FORMAT ('nmax=',I3,'   nmax1=',I3)
            CALL GSP (NMAX1,CSCA,LAM,AL1,AL2,AL3,AL4,BE1,BE2,LMAX)
            L1M=LMAX+1
            L1MAX=MAX(L1MAX,L1M)
            WGII=WG1(I)
            WGI=WGII*CSCA
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
            DO 250 L1=1,L1M
               ALPHA1(L1)=ALPHA1(L1)+AL1(L1)*WGI
               ALPHA2(L1)=ALPHA2(L1)+AL2(L1)*WGI
               ALPHA3(L1)=ALPHA3(L1)+AL3(L1)*WGI
               ALPHA4(L1)=ALPHA4(L1)+AL4(L1)*WGI
               BETA1(L1)=BETA1(L1)+BE1(L1)*WGI
               BETA2(L1)=BETA2(L1)+BE2(L1)*WGI
  250       CONTINUE
            CSCAT=CSCAT+WGI
            CEXTIN=CEXTIN+CEXT*WGII
!C           PRINT 6070, I,NMAX,NMAX1,NGAUSS
! 6070       FORMAT(4I6)
500          CONTINUE
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
             DO 510 L1=1,L1MAX
                ALPHA1(L1)=ALPHA1(L1)/CSCAT
                ALPHA2(L1)=ALPHA2(L1)/CSCAT
                ALPHA3(L1)=ALPHA3(L1)/CSCAT
                ALPHA4(L1)=ALPHA4(L1)/CSCAT
                BETA1(L1)=BETA1(L1)/CSCAT
                BETA2(L1)=BETA2(L1)/CSCAT
  510    CONTINUE
         WALB=CSCAT/CEXTIN
         CALL HOVENR(L1MAX,ALPHA1,ALPHA2,ALPHA3,ALPHA4,BETA1,BETA2)
         ASYMM=ALPHA1(2)*0.33333333333333333333_R64P
!         PRINT 9100,CEXTIN,CSCAT,WALB,ASYMM
! 9100    FORMAT('CEXT=',D12.6,2X,'CSCA=',D12.6,2X,
!     &          2X,'W=',D12.6,2X,'<COS>=',D12.6)
!         IF (WALB.GT.1._R64P) PRINT 9111
! 9111    FORMAT ('WARNING: W IS GREATER THAN 1')
!         WRITE (10,580) WALB,L1MAX
!         DO L=1,L1MAX
!            WRITE (10,575) ALPH1(L),ALPH2(L),ALPH3(L),ALPH4(L),
!     &                     BET1(L),BET2(L)           
!         ENDDO   
!  575    FORMAT(6D14.7)
!  580    FORMAT(D14.8,I8)
         LMAX=L1MAX-1
         CALL MATR (ALPHA1,ALPHA2,ALPHA3,ALPHA4,BETA1,BETA2,SCATMAT,LMAX,NPNA)
600      CONTINUE
      !ITIME=MCLOCK()
      !TIME=DFLOAT(ITIME)/6000D0
     ! PRINT 1001,TIME
 !1001 FORMAT (' time =',F8.2,' min')
       
    end subroutine
    
    subroutine CONST(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
          implicit none
          integer(I32P) :: NGAUSS,NMAX,MMAX
          real(R64P)    :: P
          real(R64P), dimension(NPNG2) :: X,W
!DIR$     ASSUME_ALIGNED X:64,W:64
          real(R64P), dimension(NPN1)  :: AN
!DIR$     ASSUME_ALIGNED AN:64
          real(R64P), dimension(NPN1,NPN1) :: ANN
!DIR$     ASSUME_ALIGNED ANN:64
          real(R64P), dimension(NPNG2)     :: S,SS
!DIR$     ASSUME_ALIGNED S:64,SS:64
          integer(I32P)                    :: NP
          real(R64P)                       :: EPS
          ! Locals
          real(R64P), dimension(NPNG1) :: X1,W1,X2,W2
!DIR$     ATTRIBUTES ALIGN : 64 :: X1
!DIR$     ATTRIBUTES ALIGN : 64 :: W1
!DIR$     ATTRIBUTES ALIGN : 64 :: X2
!DIR$     ATTRIBUTES ALIGN : 64 :: W2
          real(R64P), dimension(NPN1) :: DD
!DIR$     ATTRIBUTES ALIGN : 64 :: DD
          integer(I32P) :: N,NN,N1,NG,NG1,NG2,I
          real(R64P)    :: D,DDD,XX,Y
          !Exec code ...
          X1 = 0._R64P
          W1 = 0._R64P
          X2 = 0._R64P
          W2 = 0._R64P
          DD = 0._R64P
           DO 10 N=1,NMAX
                 NN=N*(N+1)
                 AN(N)=DFLOAT(NN)
                 D=DSQRT(DFLOAT(2*N+1)/DFLOAT(NN))
                 DD(N)=D
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
                 DO 10 N1=1,N
                       DDD=D*DD(N1)*0.5D0
                       ANN(N,N1)=DDD
                       ANN(N1,N)=DDD
        10 CONTINUE
           NG=2*NGAUSS
           IF (NP.EQ.-2) GO  TO 11
           CALL GAUSS(NG,0,0,X,W)
           GO TO 19
   11      NG1=DFLOAT(NGAUSS)/2D0
           NG2=NGAUSS-NG1
           XX=-DCOS(DATAN(EPS))
           CALL GAUSS(NG1,0,0,X1,W1)
           CALL GAUSS(NG2,0,0,X2,W2)
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
           DO 12 I=1,NG1
              W(I)=0.5_R64P*(XX+1D0)*W1(I)
              X(I)=0.5_R64P*(XX+1D0)*X1(I)+0.5_R64P*(XX-1D0)
        12 CONTINUE
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
           DO 14 I=1,NG2
              W(I+NG1)=-0.5_R64P*XX*W2(I)
              X(I+NG1)=-0.5_R64P*XX*X2(I)+0.5_R64P*XX
        14 CONTINUE
        DO 16 I=1,NGAUSS
              W(NG-I+1)=W(I)
              X(NG-I+1)=-X(I)
     16 CONTINUE
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
     19 DO 20 I=1,NGAUSS
              Y=X(I)
              Y=1._R64P/(1._R64P-Y*Y)
              SS(I)=Y
              SS(NG-I+1)=Y
              Y=DSQRT(Y)
              S(I)=Y
              S(NG-I+1)=Y
20      CONTINUE
        
    end subroutine
    
    subroutine VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII, &
                    R,DR,DDR,DRR,DRI,NMAX                   )
          implicit none
          real(R64P) :: LAM,MRR,MRI,A,EPS
          integer(I32P) :: NP,NGAUSS
          real(R64P), dimension(NPNG2) :: X
!DIR$     ASSUME_ALIGNED X:64
          real(R64P)    :: P,PPI,PIR,PII
          real(R64P), dimension(NPNG2) :: R,DR,DDR,DRR,DRI
!DIR$     ASSUME_ALIGNED R:64,DR:64,DDR:64,DRR:64,DRI:64
          integer(I32P) :: NMAX
          ! Locals
          real(R64P), dimension(NPNG2) :: Z,ZR,ZI
!DIR$     ATTRIBUTES ALIGN : 64 :: Z
!DIR$     ATTRIBUTES ALIGN : 64 :: ZR
!DIR$     ATTRIBUTES ALIGN : 64 :: ZI
          real(R64P), dimension(NPNG2,NPN1) :: J,Y,JR,JI,DJ,DJR, &
                                               DJI,DY
          integer(I32P) :: NG,I,NNMAX1,NNMAX2
          real(R64P)    :: V,PRR,PRI,TA,VV,V1,V2,TB,PI
          common  /CBESS/  J,Y,JR,JI,DJ,DY,DJR,DJI
!!DIR$     ATTRIBUTES ALIGN : 64 /CBESS/
          ! Exec code ....
          Z   = 0._R64P
          ZR  = 0._R64P
          ZI  = 0._R64P
          J   = 0._R64P
          Y   = 0._R64P
          JR  = 0._R64P
          JI  = 0._R64P
          DJ  = 0._R64P
          DJR = 0._R64P
          DJI = 0._R64P
          DY  = 0._R64P
          NG=NGAUSS*2
          IF (NP.EQ.-1) CALL RSP1(X,NG,NGAUSS,A,EPS,NP,R,DR)
          IF (NP.GE.0) CALL RSP2(X,NG,A,EPS,NP,R,DR)
          IF (NP.EQ.-2) CALL RSP3(X,NG,NGAUSS,A,EPS,R,DR)
          PI=P*2._R64P/LAM
          PPI=PI*PI
          PIR=PPI*MRR
          PII=PPI*MRI
          V=1._R64P/(MRR*MRR+MRI*MRI)
          PRR=MRR*V
          PRI=-MRI*V
          TA=0._R64P
!DIR$   VECTOR ALIGNED
!DIR$   SIMD VECTORLENGTHFOR(REAL(KIND=8))
          DO 10 I=1,NG
             VV=DSQRT(R(I))
             V=VV*PI
             TA=MAX(TA,V)
             VV=1D0/V
             DDR(I)=VV
             DRR(I)=PRR*VV
             DRI(I)=PRI*VV
             V1=V*MRR
             V2=V*MRI
             Z(I)=V
             ZR(I)=V1
             ZI(I)=V2
     10 CONTINUE
        IF (NMAX.GT.NPN1) PRINT 9000,NMAX,NPN1
        IF (NMAX.GT.NPN1) STOP
 9000   FORMAT(' NMAX = ',I2,', i.e., greater than ',I3)
        TB=TA*DSQRT(MRR*MRR+MRI*MRI)
        TB=DMAX1(TB,DFLOAT(NMAX))
        NNMAX1=1.2_R64P*DSQRT(DMAX1(TA,DFLOAT(NMAX)))+3._R64P
        NNMAX2=(TB+4._R64P*(TB**0.33333_R64P)+1.2_R64P*DSQRT(TB))
        NNMAX2=NNMAX2-NMAX+5
        CALL BESS(Z,ZR,ZI,NG,NMAX,NNMAX1,NNMAX2)
    end subroutine
    
    subroutine RSP1(X,NG,NGAUSS,REV,EPS,NP,R,DR)
          implicit none
          real(R64P), dimension(NG) :: X,R,DR
!DIR$     ASSUME_ALIGNED X:64,R:64,DR:64
          integer(I32P) :: NG,NGAUSS
          real(R64P)    :: REV,EPS
          integer(I32P) :: NP
          ! Locals
          real(R64P) :: A,AA,EE,EE1,C,CC,SS,S,RR
          integer(I32P) :: I
          ! Exec code ...
          A = REV*EPS**(0.333333333333333333_R64P)
          AA = A*A
          EE = EPS*EPS
          EE1 = EE-1._R64P
!DIR$  VECTOR ALIGNED
!DIR$  SIMD VECTORLENGTHFOR(REAL(KIND=8))
          DO 50 I=1,NGAUSS
                C=X(I)
                CC=C*C
                SS=1._R64P-CC
                S=DSQRT(SS)
                RR=1._R64P/(SS+EE*CC)
                R(I)=AA*RR
                R(NG-I+1)=R(I)
                DR(I)=RR*C*S*EE1
                DR(NG-I+1)=-DR(I)
50        CONTINUE
          
    end subroutine
    
    subroutine RSP2(X,NG,REV,EPS,N,R,DR)
          implicit none
          real(R64P), dimension(NG) :: X,R,DR
!DIR$     ASSUME_ALIGNED X:64,R:64,DR:64
          integer(I32P) :: NG
          real(R64P)    :: REV,EPS
          integer(I32P) :: N
          ! Locals
          real(R64P) :: DNP,DN,DN4,EP,A,R0,XI,RI
          integer(I32P) :: I
          ! Exec code ...
           DNP=DFLOAT(N)
           DN=DNP*DNP
           DN4=DN*4._R64P
           EP=EPS*EPS
           A=1._R64P+1.5_R64P*EP*(DN4-2._R64P)/(DN4-1._R64P)
           I=(DNP+0.1_R64P)*0.5_R64P
           I=2*I
           IF (I.EQ.N) A=A-3._R64P*EPS*(1._R64P+0.25_R64P*EP)/   &
                      (DN-1._R64P)-0.25_R64P*EP*EPS/(9._R64P*DN-1._R64P)
           R0=REV*A**(-0.33333333333333333333_R64P)
!DIR$      VECTOR ALIGNED
!DIR$      SIMD VECTORLENGTHFOR(REAL(KIND=8))
           DO 50 I=1,NG
                 XI=DACOS(X(I))*DNP
                 RI=R0*(1D0+EPS*DCOS(XI))
                 R(I)=RI*RI
                 DR(I)=-R0*EPS*DNP*DSIN(XI)/RI
50         CONTINUE
           
    end subroutine
    
    subroutine RSP3 (X,NG,NGAUSS,REV,EPS,R,DR)
          implicit none
          integer(I32P) :: NG
          real(R64P), dimension(NG) :: X,R,DR
!DIR$     ASSUME_ALIGNED X:64,R:64,DR:64
          integer(I32P) :: NGAUSS
          real(R64P)    :: REV,EPS
          ! Locals
          real(R64P) :: H,A,CO,SI,RAD,RTHET
          integer(I32P) :: I
          ! Exec code
          H=REV*( (2._R64P/(3._R64P*EPS*EPS))**(0.3333333333333333333333333_R64P) )
          A=H*EPS
          DO 50 I=1,NGAUSS
                CO=-X(I)
                SI=DSQRT(1._R64P-CO*CO)
                IF (SI/CO.GT.A/H) GO TO 20
                    RAD=H/CO
                    RTHET=H*SI/(CO*CO)
                    GO TO 30
   20          RAD=A/SI
               RTHET=-A*CO/(SI*SI)
   30          R(I)=RAD*RAD   
               R(NG-I+1)=R(I)
               DR(I)=-RTHET/RAD
               DR(NG-I+1)=-DR(I)
   50 CONTINUE
     
    end subroutine
    
    subroutine BESS(X,XR,XI,NG,NMAX,NNMAX1,NNMAX2)
          implicit none
          real(R64P), dimension(NG) :: X,XR,XI
!DIR$     ASSUME_ALIGNED X:64,XR:64,XI:64
          integer(I32P) :: NG,NMAX,NNMAX1,NNMAX2
          real(R64P), dimension(NPNG2,NPN1) :: J,Y,JR,JI,DJ,DY,  &
                                               DJR,DJI
          ! Locals
          real(R64P), dimension(NPN1) :: AJ,AY,AJR,AJI,ADJ,ADY,ADJR,ADJI
!DIR$     ATTRIBUTES ALIGN : 64 :: AJ
!DIR$     ATTRIBUTES ALIGN : 64 :: AY
!DIR$     ATTRIBUTES ALIGN : 64 :: AJR
!DIR$     ATTRIBUTES ALIGN : 64 :: AJI
!DIR$     ATTRIBUTES ALIGN : 64 :: ADJ
!DIR$     ATTRIBUTES ALIGN : 64 :: ADY
!DIR$     ATTRIBUTES ALIGN : 64 :: ADJR
!DIR$     ATTRIBUTES ALIGN : 64 :: ADJI
          integer(I32P) :: I,N
          real(R64P)    :: XX,YR,YI
          common /CBESS/  J,Y,JR,JI,DJ,DY,DJR,DJI
!!DIR$     ATTRIBUTES ALIGN : 64 :: /CBESS/
          ! Exec code ...
          J  = 0._R64P
          Y  = 0._R64P
          JR = 0._R64P
          JI = 0._R64P
          DJ = 0._R64P
          DY = 0._R64P
          DJR = 0._R64P
          DJI = 0._R64P
          AJ  = 0._R64P
          AY  = 0._R64P
          AJR = 0._R64P
          AJI = 0._R64P
          ADJ = 0._R64P
          ADY = 0._R64P
          ADJR = 0._R64P
          ADJI = 0._R64P
          DO 10 I=1,NG
             XX=X(I)
             CALL RJB(XX,AJ,ADJ,NMAX,NNMAX1)
             CALL RYB(XX,AY,ADY,NMAX)
             YR=XR(I)
             YI=XI(I)
             CALL CJB(YR,YI,AJR,AJI,ADJR,ADJI,NMAX,NNMAX2)

             DO 10 N=1,NMAX
                   J(I,N)=AJ(N)
                   Y(I,N)=AY(N)
                   JR(I,N)=AJR(N)
                   JI(I,N)=AJI(N)
                   DJ(I,N)=ADJ(N)
                   DY(I,N)=ADY(N)
                   DJR(I,N)=ADJR(N)
                   DJI(I,N)=ADJI(N)
10        CONTINUE
          
    end subroutine
    
      subroutine RJB(X,Y,U,NMAX,NNMAX)
          implicit none
          real(R64P) :: X
          real(R64P), dimension(NMAX) :: Y
!DIR$     ASSUME_ALIGNED Y:64
          real(R64P), dimension(NMAX)  :: U
!DIR$     ASSUME_ALIGNED U:64
          real(R64P), dimension(800) :: Z
!DIR$     ATTRIBUTES ALIGN : 64 :: Z
          integer(I32P) :: NMAX,NNMAX
          ! Locals
          integer(I32P) :: L,L1,I,I1
          real(R64P)    :: XX,Z0,Y0,Y1,YI,YI1
          ! Exec code ....
          Z = 0._R64P
          L=NMAX+NNMAX
          XX=1._R64P/X
          Z(L)=1D0/(DFLOAT(2*L+1)*XX)
          L1=L-1
          DO 5 I=1,L1
               I1=L-I
               Z(I1)=1._R64P/(DFLOAT(2*I1+1)*XX-Z(I1+1))
       5 CONTINUE
         Z0=1._R64P/(XX-Z(1))
         Y0=Z0*DCOS(X)*XX
         Y1=Y0*Z(1)
         U(1)=Y0-Y1*XX
         Y(1)=Y1
         DO 10 I=2,NMAX
               YI1=Y(I-1)
               YI=YI1*Z(I)
               U(I)=YI1-DFLOAT(I)*YI*XX
               Y(I)=YI
     10 CONTINUE
     
    end subroutine
    
      subroutine RYB(X,Y,V,NMAX)
          implicit none
          real(R64P) :: X
          real(R64P), dimension(NMAX) :: Y,V
!DIR$     ASSUME_ALIGNED Y:64,V:64
          integer(I32P) :: NMAX
          ! Locals
          real(R64P) :: C,S,X1,X2,X3,Y1
          integer(I32P) :: I,NMAX1
          ! Exec code ...
          C=DCOS(X)
          S=DSIN(X)
          X1=1._R64P/X
          X2=X1*X1
          X3=X2*X1
          Y1=-C*X2-S*X1
          Y(1)=Y1
          Y(2)=(-3._R64P*X3+X1)*C-3._R64P*X2*S
          NMAX1=NMAX-1
          DO 5 I=2,NMAX1
    5            Y(I+1)=DFLOAT(2*I+1)*X1*Y(I)-Y(I-1)
                 V(1)=-X1*(C+Y1)
                 DO 10 I=2,NMAX
  10                   V(I)=Y(I-1)-DFLOAT(I)*X1*Y(I)
     
    end subroutine
    
!C**********************************************************************
!C                                                                     *
!C   CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND       *
!C   J=JR+I*JI OF COMPLEX ARGUMENT X=XR+I*XI OF ORDERS FROM 1 TO NMAX  *
!C   BY USING BACKWARD RECURSION. PARAMETR NNMAX DETERMINES NUMERICAL  *
!C   ACCURACY. U=UR+I*UI - FUNCTION (1/X)(D/DX)(X*J(X))                *
!C                                                                     *
!C**********************************************************************
 
      subroutine CJB (XR,XI,YR,YI,UR,UI,NMAX,NNMAX)
          implicit none
          real(R64P) :: XR,XI
          real(R64P), dimension(NMAX) ::  YR,YI,UR,UI
!DIR$     ASSUME_ALIGNED YR:64,YI:64,UR:64,UI:64
          integer(I32P) :: NMAX,NNMAX
          ! Locals
          real(R64P), dimension(NPN1) :: CYR,CYI,CUR,CUI
!DIR$     ATTRIBUTES ALIGN : 64 :: CYR
!DIR$     ATTRIBUTES ALIGN : 64 :: CYI
!DIR$     ATTRIBUTES ALIGN : 64 :: CUR
!DIR$     ATTRIBUTES ALIGN : 64 :: CUI
          real(R64P), dimension(1200) :: CZR,CZI
!DIR$     ATTRIBUTES ALIGN : 64 :: CZR
!DIR$     ATTRIBUTES ALIGN : 64 :: CZI
          integer(I32P) :: L,L1,I,I1
          real(R64P)    :: XRXI,CXXR,CXXI,QF,AR,ARI,AI,  &
                           CZ0R,CZ0I,CR,CI,CY0R,CY1R,CY0I,CY1I , &
                           CU1R,CU1I,QI,CYI1R,CYI1I,CYII,CYIR, &
                           CUIR,CUII
          CYR = 0._R64P
          CYI = 0._R64P
          CUR = 0._R64P
          CUI = 0._R64P
          CZR = 0._R64P
          CZI = 0._R64P
          L=NMAX+NNMAX
          XRXI=1._R64P/(XR*XR+XI*XI)
          CXXR=XR*XRXI
          CXXI=-XI*XRXI 
          QF=1._R64P/DFLOAT(2*L+1)
          CZR(L)=XR*QF
          CZI(L)=XI*QF
          L1=L-1
          DO I=1,L1
             I1=L-I
             QF=DFLOAT(2*I1+1)
             AR=QF*CXXR-CZR(I1+1)
             AI=QF*CXXI-CZI(I1+1)
             ARI=1D0/(AR*AR+AI*AI)
             CZR(I1)=AR*ARI
             CZI(I1)=-AI*ARI
          ENDDO   
          AR=CXXR-CZR(1)
          AI=CXXI-CZI(1)
          ARI=1._R64P/(AR*AR+AI*AI)
          CZ0R=AR*ARI
          CZ0I=-AI*ARI
          CR=DCOS(XR)*DCOSH(XI)
          CI=-DSIN(XR)*DSINH(XI)
          AR=CZ0R*CR-CZ0I*CI
          AI=CZ0I*CR+CZ0R*CI
          CY0R=AR*CXXR-AI*CXXI
          CY0I=AI*CXXR+AR*CXXI
          CY1R=CY0R*CZR(1)-CY0I*CZI(1)
          CY1I=CY0I*CZR(1)+CY0R*CZI(1)
          AR=CY1R*CXXR-CY1I*CXXI
          AI=CY1I*CXXR+CY1R*CXXI
          CU1R=CY0R-AR
          CU1I=CY0I-AI
          CYR(1)=CY1R
          CYI(1)=CY1I
          CUR(1)=CU1R
          CUI(1)=CU1I
          YR(1)=CY1R
          YI(1)=CY1I
          UR(1)=CU1R
          UI(1)=CU1I
          DO I=2,NMAX
             QI=DFLOAT(I)
             CYI1R=CYR(I-1)
             CYI1I=CYI(I-1)
             CYIR=CYI1R*CZR(I)-CYI1I*CZI(I)
             CYII=CYI1I*CZR(I)+CYI1R*CZI(I)
             AR=CYIR*CXXR-CYII*CXXI
             AI=CYII*CXXR+CYIR*CXXI
             CUIR=CYI1R-QI*AR
             CUII=CYI1I-QI*AI
             CYR(I)=CYIR
             CYI(I)=CYII
             CUR(I)=CUIR
             CUI(I)=CUII
             YR(I)=CYIR
             YI(I)=CYII
             UR(I)=CUIR
             UI(I)=CUII
        ENDDO   
       
    end subroutine
    
    subroutine TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,    &
                       DRR,DRI,NMAX,NCHECK)
          implicit none
          integer(I32P) :: NGAUSS
          real(R64P), dimension(NPNG2) :: X,W,S,SS
!DIR$     ASSUME_ALIGNED X:64,W:64,AN:64,S:64,SS:64
          real(R64P), dimension(NPN1) :: AN
!DIR$     ASSUME_ALIGNED AN:64
          real(R64P), dimension(NPN1,NPN1) :: ANN
!DIR$     ASSUME_ALIGNED ANN:64
          real(R64P) :: PPI,PIR,PII
          real(R64P), dimension(NPNG2) :: R,DR,DDR,DRR,DRI
!DIR$     ASSUME_ALIGNED R:64,DR:64,DDR:64,DRR:64,DRI:64
          integer(I32P) :: NMAX,NCHECK
          ! Locals
          real(R64P), dimension(NPNG2) :: SIG 
!DIR$     ATTRIBUTES ALIGN : 64 :: SIG
          real(R64P), dimension(NPNG2,NPN1) :: J,Y,JR,JI,DJ,DY,DJR,  &
                                               DJI,D1,D2
!DIR$     ATTRIBUTES ALIGN : 64 :: D1
!DIR$     ATTRIBUTES ALIGN : 64 :: D2
          real(R64P), dimension(NPNG2) :: DS,DSS,RR,DV1,DV2
!DIR$     ATTRIBUTES ALIGN : 64 :: DS
!DIR$     ATTRIBUTES ALIGN : 64 :: DSS
!DIR$     ATTRIBUTES ALIGN : 64 :: RR
!DIR$     ATTRIBUTES ALIGN : 64 :: DV1
!DIR$     ATTRIBUTES ALIGN : 64 :: DV2          
          real(R64P), dimension(NPN1,NPN1) :: R11,R12,R21,R22,I11,I12,  &
                                              I21,I22,RG11,RG12,RG21,  &
                                              RG22,IG11,IG12,IG21,IG22
          real(R64P), dimension(NPN2,NPN2) :: QR,QI,RGQR,RGQI,TQR,TQI, &
                                              TRGQR,TRGQI,TR1,TI1
!DIR$     ATTRIBUTES ALIGN : 64 :: TQR          
!DIR$     ATTRIBUTES ALIGN : 64 :: TQI
!DIR$     ATTRIBUTES ALIGN : 64 :: TRGQR
!DIR$     ATTRIBUTES ALIGN : 64 :: TRGQI
          real(R32P), dimension(NPN6*NPN4*NPN4*8) :: PLUS
!DIR$     ATTRIBUTES ALIGN : 64 :: PLUS
          integer(I32P) :: MM1,NNMAX,NG,NGSS,N,I,I1,I2,   &
                           N1,N2,K1,KK1,K2,KK2,NM
          real(R64P)    :: FACTOR,SI,DD1,DD2,AN1,AN2,     &
                           AR12,AR21,AI12,AI21,GR12,      &
                           GR21,GI12,GI21,D1N1,D2N1,      &
                           A12,A21,A22,AA1,QJ1,QY1,       &
                           QJR2,QJI2,QDJR2,QDJI2,QDJ1,    &
                           QDY1,C1R,C1I,B1R,B1I,      &
                           C2R,C2I,B2R,B2I,DDRI,C3R,      &
                           C3I,B3R,B3I,C4R,C4I,B4R,       &
                           B4I,DRRI,DRII,C5R,C5I,B5R,     &
                           URI,RRI,F1,F2,AN12,TPIR,       &
                           TPII,TPPI,TAR12,TAI12,TGR12,   &
                           TGI12,TAR21,TAI21,TGR21,TGI21, &
                           D1N2,D2N2,B5I
                           
      COMMON /TMAT/ PLUS,   &
                 R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,  &
                 IG11,IG12,IG21,IG22
!!DIR$   ATTRIBUTES ALIGN : 64 :: /TMAT/
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
!!DIR$   ATTRIBUTES ALIGN : 64 :: /CBESS/
      COMMON /CT/ TR1,TI1
!!DIR$   ATTRIBUTES ALIGN : 64 :: /CT/
      COMMON /CTT/ QR,QI,RGQR,RGQI
!!DIR$   ATTRIBUTES ALIGN : 64 :: /CTT/  
      SIG = 0._R64P
      J   = 0._R64P
      Y   = 0._R64P
      JR  = 0._R64P
      JI  = 0._R64P
      DJ  = 0._R64P
      DY  = 0._R64P
      DJR = 0._R64P
      DJI = 0._R64P
      D1  = 0._R64P
      D2  = 0._R64P
      DS  = 0._R64P 
      DSS = 0._R64P 
      RR  = 0._R64P 
      DV1 = 0._R64P
      DV2 = 0._R64P
      R11 = 0._R64P 
      R12 = 0._R64P
      R21 = 0._R64P 
      R22 = 0._R64P 
      I11 = 0._R64P 
      I12 = 0._R64P
      I21 = 0._R64P
      I22 = 0._R64P 
      RG11 = 0._R64P 
      RG12 = 0._R64P 
      RG21 = 0._R64P
      RG22 = 0._R64P 
      IG11 = 0._R64P 
      IG12 = 0._R64P 
      IG21 = 0._R64P 
      IG22 = 0._R64P
      QR   = 0._R64P 
      QI   = 0._R64P 
      RGQR = 0._R64P 
      RGQI = 0._R64P 
      TQR  = 0._R64P 
      TQI  = 0._R64P
      TRGQR = 0._R64P 
      TRGQI = 0._R64P 
      TR1   = 0._R64P
      TI1   = 0._R64P
      PLUS  = 0._R64P
      MM1=1
      NNMAX=NMAX+NMAX
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1D0
      IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2D0
         ELSE
            CONTINUE
      ENDIF
      SI=1D0
      DO 5 N=1,NNMAX
           SI=-SI
           SIG(N)=SI
    5 CONTINUE
   20 DO 25 I=1,NGAUSS
         I1=NGAUSS+I
         I2=NGAUSS-I+1
         CALL VIG ( X(I1), NMAX, 0, DV1, DV2)
!DIR$    VECTOR ALIGNED
!DIR$    SIMD VECTORLENGTHFOR(REAL(KIND=8))
         DO 25 N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI
25    CONTINUE
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
   30 DO 40 I=1,NGSS
           RR(I)=W(I)*R(I)
   40 CONTINUE
 
      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR12=0._R64P
                AR21=0._R64P
                AI12=0._R64P
                AI21=0._R64P
                GR12=0._R64P
                GR21=0._R64P
                GI12=0._R64P
                GI21=0._R64P
                IF (NCHECK.EQ.1.AND.SIG(N1+N2).LT.0._R64P) GO TO 205
!DIR$   FMA
                DO 200 I=1,NGSS
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21
 
                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)
 
                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1
 
                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1
 
                    DDRI=DDR(I)
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
 
                    DRRI=DRR(I)
                    DRII=DRI(I)
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII
 
                    URI=DR(I)
                    RRI=RR(I)
 
                    F1=RRI*A22
                    F2=RRI*URI*AN1*A12
                    AR12=AR12+F1*B2R+F2*B3R
                    AI12=AI12+F1*B2I+F2*B3I
                    GR12=GR12+F1*C2R+F2*C3R
                    GI12=GI12+F1*C2I+F2*C3I
 
                    F2=RRI*URI*AN2*A21
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
  200           CONTINUE
 
  205           AN12=ANN(N1,N2)*FACTOR
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
  300 CONTINUE
 
      TPIR=PIR
      TPII=PII
      TPPI=PPI
 
      NM=NMAX
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
!DIR$   FMA
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=0._R64P
                TQI(K1,KK2)=0._R64P
                TRGQR(K1,KK2)=0._R64P
                TRGQI(K1,KK2)=0._R64P
 
                TQR(KK1,K2)=0._R64P
                TQI(KK1,K2)=0._R64P
                TRGQR(KK1,K2)=0._R64P
                TRGQI(KK1,K2)=0._R64P
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX

           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
      CALL TT(NMAX,NCHECK)
      
    end subroutine
    
    subroutine TMATR (M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,    &
                       DRR,DRI,NMAX,NCHECK)
          implicit none
          integer(I32P) :: M,NGAUSS
          real(R64P), dimension(NPNG2) :: X,W,S,SS
!DIR$     ASSUME_ALIGNED X:64,W:64,S:64,SS:64
          real(R64P), dimension(NPN1) :: AN
!DIR$     ASSUME_ALIGNED AN:64
          real(R64P), dimension(NPN1,NPN1) :: ANN
!DIR$     ASSUME_ALIGNED ANN:64
          real(R64P) :: PPI,PIR,PII
          real(R64P), dimension(NPNG2) :: R,DR,DDR,DRR,DRI
!DIR$     ASSUME_ALIGNED R:64,DR:64,DDR:64,DRR:64,DRI:64
          integer(I32P) :: NMAX,NCHECK
          ! Locals
          real(R64P), dimension(NPNG2) :: SIG 
!DIR$     ATTRIBUTES ALIGN : 64 :: SIG
          real(R64P), dimension(NPNG2,NPN1) :: J,Y,JR,JI,DJ,DY,DJR,  &
                                               DJI,D1,D2
!DIR$     ATTRIBUTES ALIGN : 64 :: D1
!DIR$     ATTRIBUTES ALIGN : 64 :: D2
          real(R64P), dimension(NPNG2) :: DS,DSS,RR,DV1,DV2
!DIR$     ATTRIBUTES ALIGN : 64 :: DS
!DIR$     ATTRIBUTES ALIGN : 64 :: DSS
!DIR$     ATTRIBUTES ALIGN : 64 :: RR
!DIR$     ATTRIBUTES ALIGN : 64 :: DV1
!DIR$     ATTRIBUTES ALIGN : 64 :: DV2          
          real(R64P), dimension(NPN1,NPN1) :: R11,R12,R21,R22,I11,I12,  &
                                              I21,I22,RG11,RG12,RG21,  &
                                              RG22,IG11,IG12,IG21,IG22
          real(R64P), dimension(NPN2,NPN2) :: QR,QI,RGQR,RGQI,TQR,TQI, &
                                              TRGQR,TRGQI,TR1,TI1
!DIR$     ATTRIBUTES ALIGN : 64 :: TQR          
!DIR$     ATTRIBUTES ALIGN : 64 :: TQI
!DIR$     ATTRIBUTES ALIGN : 64 :: TRGQR
!DIR$     ATTRIBUTES ALIGN : 64 :: TRGQI
          real(R32P), dimension(NPN6*NPN4*NPN4*8) :: PLUS
!DIR$     ATTRIBUTES ALIGN : 64 :: PLUS
          integer(I32P) :: MM1,NNMAX,NG,NGSS,N,I,I1,I2,   &
                           N1,N2,K1,KK1,K2,KK2,NM
          real(R64P)    :: FACTOR,QM,QMM,WR,SI,DD1,DD2,AN1,AN2,     &
                           AR11,AR12,AR21,AR22,AI11,A11,AA2,        &
                           AI12,AI21,AI22,GR11,GR12,GR22,           &
                           GR21,GI11,GI12,GI21,GI22,D1N1,D2N1,      &
                           A12,A21,A22,AA1,QJ1,QY1,       &
                           QJR2,QJI2,QDJR2,QDJI2,QDJ1,    &
                           QDY1,C1R,C1I,B1R,B1I,      &
                           C2R,C2I,B2R,B2I,DDRI,C3R,      &
                           C3I,B3R,B3I,C4R,C4I,B4R,       &
                           B4I,DRRI,DRII,C5R,C5I,B5R,     &
                           URI,RRI,F1,F2,AN12,TPIR,       &
                           TPII,TPPI,TAR12,TAI12,TGR12,   &
                           TGI12,TAR21,TAI21,TGR21,TGI21, &
                           D1N2,D2N2,B5I,C6R,C6I,B6R,B6I, &
                           C7R,C7I,B7R,B7I,C8R,C8I,B8R,   &
                           B8I,DSI,DSSI,E1,E2,E3, TAR11,  &
                           TAI11,TGR11,TGI11,TAR22,TAI22, &
                           TGR22,TGI22
                           
      COMMON /TMAT/ PLUS,   &
                 R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,  &
                 IG11,IG12,IG21,IG22
!!DIR$   ATTRIBUTES ALIGN : 64 :: /TMAT/
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
!!DIR$   ATTRIBUTES ALIGN : 64 :: /CBESS/
      COMMON /CT/ TR1,TI1
!!DIR$   ATTRIBUTES ALIGN : 64 :: /CT/
      COMMON /CTT/ QR,QI,RGQR,RGQI
!!DIR$   ATTRIBUTES ALIGN : 64 :: /CTT/ 
      SIG = 0._R64P
      J   = 0._R64P
      Y   = 0._R64P
      JR  = 0._R64P
      JI  = 0._R64P
      DJ  = 0._R64P
      DY  = 0._R64P
      DJR = 0._R64P
      DJI = 0._R64P
      D1  = 0._R64P
      D2  = 0._R64P
      DS  = 0._R64P 
      DSS = 0._R64P 
      RR  = 0._R64P 
      DV1 = 0._R64P
      DV2 = 0._R64P
      R11 = 0._R64P 
      R12 = 0._R64P
      R21 = 0._R64P 
      R22 = 0._R64P 
      I11 = 0._R64P 
      I12 = 0._R64P
      I21 = 0._R64P
      I22 = 0._R64P 
      RG11 = 0._R64P 
      RG12 = 0._R64P 
      RG21 = 0._R64P
      RG22 = 0._R64P 
      IG11 = 0._R64P 
      IG12 = 0._R64P 
      IG21 = 0._R64P 
      IG22 = 0._R64P
      QR   = 0._R64P 
      QI   = 0._R64P 
      RGQR = 0._R64P 
      RGQI = 0._R64P 
      TQR  = 0._R64P 
      TQI  = 0._R64P
      TRGQR = 0._R64P 
      TRGQI = 0._R64P 
      TR1   = 0._R64P
      TI1   = 0._R64P
      PLUS  = 0._R64P
      MM1=M
      QM=DFLOAT(M)
      QMM=QM*QM
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1._R64P
       IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2D0
         ELSE
            CONTINUE
      ENDIF
      SI=1._R64P
      NM=NMAX+NMAX
      DO 5 N=1,NM
           SI=-SI
           SIG(N)=SI
    5 CONTINUE
   20 DO 25 I=1,NGAUSS
         I1=NGAUSS+I
         I2=NGAUSS-I+1
         CALL VIG (X(I1),NMAX,M,DV1,DV2)
         DO 25 N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI
25    CONTINUE
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
   30 DO 40 I=1,NGSS
           WR=W(I)*R(I)
           DS(I)=S(I)*QM*WR
           DSS(I)=SS(I)*QMM
           RR(I)=WR
   40 CONTINUE
 
      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR11=0._R64P
                AR12=0._R64P
                AR21=0._R64P
                AR22=0._R64P
                AI11=0._R64P
                AI12=0._R64P
                AI21=0._R64P
                AI22=0._R64P
                GR11=0._R64P
                GR12=0._R64P
                GR21=0._R64P
                GR22=0._R64P
                GI11=0._R64P
                GI12=0._R64P
                GI21=0._R64P
                GI22=0._R64P
                SI=SIG(N1+N2)
!DIR$   FMA 
                DO 200 I=1,NGSS
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A11=D1N1*D1N2
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21
                    AA2=A11*DSS(I)+A22
                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)
 
                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1
 
                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1
 
                    DDRI=DDR(I)
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
 
                    DRRI=DRR(I)
                    DRII=DRI(I)
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII
 
                    C6R=QDJR2*QDJ1
                    C6I=QDJI2*QDJ1
                    B6R=C6R-QDJI2*QDY1
                    B6I=C6I+QDJR2*QDY1
 
                    C7R=C4R*DDRI
                    C7I=C4I*DDRI
                    B7R=B4R*DDRI
                    B7I=B4I*DDRI
 
                    C8R=C2R*DRRI-C2I*DRII
                    C8I=C2I*DRRI+C2R*DRII
                    B8R=B2R*DRRI-B2I*DRII
                    B8I=B2I*DRRI+B2R*DRII
 
                    URI=DR(I)
                    DSI=DS(I)
                    DSSI=DSS(I)
                    RRI=RR(I)
 
                    IF (NCHECK.EQ.1.AND.SI.GT.0._R64P) GO TO 150
 
                    E1=DSI*AA1
                    AR11=AR11+E1*B1R
                    AI11=AI11+E1*B1I
                    GR11=GR11+E1*C1R
                    GI11=GI11+E1*C1I
                    IF (NCHECK.EQ.1) GO TO 160
 
  150               F1=RRI*AA2
                    F2=RRI*URI*AN1*A12
                    AR12=AR12+F1*B2R+F2*B3R
                    AI12=AI12+F1*B2I+F2*B3I
                    GR12=GR12+F1*C2R+F2*C3R
                    GI12=GI12+F1*C2I+F2*C3I
 
                    F2=RRI*URI*AN2*A21
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
                    IF (NCHECK.EQ.1) GO TO 200
 
  160               E2=DSI*URI*A11
                    E3=E2*AN2
                    E2=E2*AN1
                    AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                    AI22=AI22+E1*B6I+E2*B7I+E3*B8I
                    GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                    GI22=GI22+E1*C6I+E2*C7I+E3*C8I
  200           CONTINUE
                AN12=ANN(N1,N2)*FACTOR
                R11(N1,N2)=AR11*AN12
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                R22(N1,N2)=AR22*AN12
                I11(N1,N2)=AI11*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                I22(N1,N2)=AI22*AN12
                RG11(N1,N2)=GR11*AN12
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                RG22(N1,N2)=GR22*AN12
                IG11(N1,N2)=GI11*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
                IG22(N1,N2)=GI22*AN12
 
  300 CONTINUE
      TPIR=PIR
      TPII=PII
      TPPI=PPI
      NM=NMAX-MM1+1
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
!DIR$  FMA
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR11=-R11(N1,N2)
                TAI11=-I11(N1,N2)
                TGR11=-RG11(N1,N2)
                TGI11=-IG11(N1,N2)
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TAR22=-R22(N1,N2)
                TAI22=-I22(N1,N2)
                TGR22=-RG22(N1,N2)
                TGI22=-IG22(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=TPIR*TAR11-TPII*TAI11+TPPI*TAR22
                TQI(K1,KK2)=TPIR*TAI11+TPII*TAR11+TPPI*TAI22
                TRGQR(K1,KK2)=TPIR*TGR11-TPII*TGI11+TPPI*TGR22
                TRGQI(K1,KK2)=TPIR*TGI11+TPII*TGR11+TPPI*TGI22
 
                TQR(KK1,K2)=TPIR*TAR22-TPII*TAI22+TPPI*TAR11
                TQI(KK1,K2)=TPIR*TAI22+TPII*TAR22+TPPI*TAI11
                TRGQR(KK1,K2)=TPIR*TGR22-TPII*TGI22+TPPI*TGR11
                TRGQI(KK1,K2)=TPIR*TGI22+TPII*TGR22+TPPI*TGI11
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
 
      CALL TT(NM,NCHECK)
 
     
    end subroutine
    
    subroutine VIG (X, NMAX, M, DV1, DV2)
      implicit none
      real(R64P) :: X
      integer(I32P) :: NMAX,M
      real(R64P), dimension(NPN1) :: DV1,DV2
!DIR$ ASSUME_ALIGNED VD1:64, DV2:64 
      ! Locals
      real(R64P) :: A,QS,QS1,D1,D2,QN,QN1,QN2,D3,DER,  &
                    QMM,QNM,QNM1
      integer(I32P) :: I,N,I2
      ! Exec code ...
      A=1._R64P
      QS=DSQRT(1._R64P-X*X)
      QS1=1._R64P/QS
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
      DO N=1,NMAX
         DV1(N)=0._R64P
         DV2(N)=0._R64P
      ENDDO   
      IF (M.NE.0) GO TO 20
      D1=1._R64P
      D2=X  
      DO N=1,NMAX  
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1 
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO   
      RETURN
   20 QMM=DFLOAT(M*M)
      DO I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS
      ENDDO   
      D1=0._R64P
      D2=A 
      DO N=M,NMAX
         QN=DFLOAT(N)
         QN2=DFLOAT(2*N+1)
         QN1=DFLOAT(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO   
      
    end subroutine
    
!    C**********************************************************************
!C                                                                     *
!C   CALCULATION OF THE MATRIX    T = - RG(Q) * (Q**(-1))              *
!C                                                                     *
!C   INPUT INFORTMATION IS IN COMMON /CTT/                             *
!C   OUTPUT INFORMATION IS IN COMMON /CT/                              *
!C                                                                     *
!C**********************************************************************
 
      subroutine TT(NMAX,NCHECK)
      implicit none
      integer(I32P) :: NMAX,NCHECK
      real(R64P), dimension(NPN2,NPN2) :: F,B,QR,QI,RGQR,RGQI,  &
                                          A,C,D,E
!DIR$ ATTRIBUTES ALIGN : 64 :: F
!DIR$ ATTRIBUTES ALIGN : 64 :: B
!DIR$ ATTRIBUTES ALIGN : 64 :: A
!DIR$ ATTRIBUTES ALIGN : 64 :: C
!DIR$ ATTRIBUTES ALIGN : 64 :: D
!DIR$ ATTRIBUTES ALIGN : 64 :: E  
      real(R64P), dimension(NPN2) :: WORK
!DIR$ ATTRIBUTES ALIGN : 64 :: WORK
      real(R64P), dimension(NPN2,NPN2) :: TR1,TI1
      COMPLEX*16 ZQ(NPN2,NPN2),ZW(NPN2)
      integer(I32P), dimension(NPN2) :: IPIV,IPVT
      integer(I32P) :: NDIM,NNMAX,I,J,INFO,K
      real(R64P)    :: TR,TI,ARR,ARI,AR,AI
      COMMON /CT/ TR1,TI1
      COMMON /CTT/ QR,QI,RGQR,RGQI
!!DIR$ ATTRIBUTES ALIGN : 64 :: /CT/
!!DIR$ ATTRIBUTES ALIGN : 64 :: /CTT/
      ! Exec code ....
      F = 0._R64P
      B = 0._R64P
      QR = 0._R64P
      QI = 0._R64P
      RGQR = 0._R64P
      RGQI = 0._R64P
      A = 0._R64P
      C = 0._R64P
      D = 0._R64P
      E = 0._R64P
      WORK = 0._R64P
      TR1 = 0._R64P
      TI1 = 0._R64P
      ZW = DCMPLX(0._R64P,0._R64P)
      ZQ = DCMPLX(0._R64P,0._R64P)
      IPIV = 0
      IPVT = 0
      NDIM=NPN2
      NNMAX=2*NMAX
 
!C     Matrix inversion from LAPACK 
 
      DO I=1,NNMAX
	   DO J=1,NNMAX
	      ZQ(I,J)=DCMPLX(QR(I,J),QI(I,J))
	   ENDDO
      ENDDO
      INFO=0
      CALL ZGETRF(NNMAX,NNMAX,ZQ,NPN2,IPIV,INFO)
      IF (INFO.NE.0) WRITE (6,1100) INFO
      CALL ZGETRI(NNMAX,ZQ,NPN2,IPIV,ZW,NPN2,INFO)
      IF (INFO.NE.0) WRITE (6,1100) INFO

 1100 FORMAT ('WARNING:  info=', i2)
      DO I=1,NNMAX
         DO J=1,NNMAX
            TR=0._R64P
            TI=0._R64P
!DIR$   FMA
	    DO K=1,NNMAX
                 ARR=RGQR(I,K)
                 ARI=RGQI(I,K)
                 AR=ZQ(K,J)
                 AI=DIMAG(ZQ(K,J))
                 TR=TR-ARR*AR+ARI*AI
                 TI=TI-ARR*AI-ARI*AR
            ENDDO
	        TR1(I,J)=TR
            TI1(I,J)=TI
         ENDDO
      ENDDO
     
    end subroutine
    
!    C********************************************************************
!C                                                                   *
!C   CALCULATION OF THE EXPANSION COEFFICIENTS FOR (I,Q,U,V) -       *
!C   REPRESENTATION.                                                 *
!C                                                                   *
!C   INPUT PARAMETERS:                                               *
!C                                                                   *
!C      LAM - WAVELENGTH OF LIGHT                                    *
!C      CSCA - SCATTERING CROSS SECTION                              *
!C      TR AND TI - ELEMENTS OF THE T-MATRIX. TRANSFERRED THROUGH    *
!C                  COMMON /CTM/                                     *
!C      NMAX - DIMENSION OF T(M)-MATRICES                            *
!C                                                                   *
!C   OUTPUT INFORTMATION:                                            *
!C                                                                   *
!C      ALF1,...,ALF4,BET1,BET2 - EXPANSION COEFFICIENTS             *
!C      LMAX - NUMBER OF COEFFICIENTS MINUS 1                        *
!C                                                                   *
!C********************************************************************
 
    subroutine GSP(NMAX,CSCA,LAM,ALF1,ALF2,ALF3,ALF4,BET1,BET2,LMAX)
          implicit none
          integer(I32P) :: NMAX
          real(R64P)    :: CSCA,LAM
          real(R64P), dimension(NPL) :: ALF1,ALF2,ALF3,ALF4,BET1,BET2
!DIR$     ASSUME_ALIGNED ALF1:64,ALF2:64,ALF3:64,ALF4:64,BET1:64,BET2:64
          integer(I32P) :: LMAX
          ! Locals
          real(R64P), dimension(900) :: SSIGN
          real(R64P), dimension(NPL) :: SSI
!DIR$     ATTRIBUTES ALIGN : 64 :: SSI
          real(R64P), dimension(NPN1) :: SSJ
!DIR$     ATTRIBUTES ALIGN : 64 :: SSJ
          real(R64P), dimension(NPL1,NPN4) :: TR1,TR2,TI1,TI2
          real(R64P), dimension(NPL1,NPN6) :: G1,G2
          real(R64P), dimension(NPN4) :: AR1,AR2,AI1,AI2
!DIR$     ATTRIBUTES ALIGN : 64 :: AR1     
!DIR$     ATTRIBUTES ALIGN : 64 :: AR2
!DIR$     ATTRIBUTES ALIGN : 64 :: AI1
!DIR$     ATTRIBUTES ALIGN : 64 :: AI2          
          real(R64P), dimension(NPN4,NPN4) :: FR,FI,FF
!DIR$     ATTRIBUTES ALIGN : 64 :: FR
!DIR$     ATTRIBUTES ALIGN : 64 :: FI
!DIR$     ATTRIBUTES ALIGN : 64 :: FF
          real(R32P), dimension(NPL1,NPL1,NPN4) :: B1R,B1I,B2R,B2I
          real(R32P), dimension(NPL1,NPN4,NPN4) :: D1,D2,D3,D4,D5R,D5I
!!DIR$     ATTRIBUTES ALIGN : 64 :: D1
!!DIR$     ATTRIBUTES ALIGN : 64 :: D2
!!DIR$     ATTRIBUTES ALIGN : 64 :: D2
!!DIR$     ATTRIBUTES ALIGN : 64 :: D4
!!DIR$     ATTRIBUTES ALIGN : 64 :: D5R
!!DIR$     ATTRIBUTES ALIGN : 64 :: D5I
          real(R32P), dimension(NPN6*NPN4*NPN4*8) :: PLUS1
!DIR$     ATTRIBUTES ALIGN : 64 :: PLUS1
          real(R32P), dimension(NPN6,NPN4,NPN4) :: TR11,TR12,TR21,TR22, & 
                                                   TI11,TI12,TI21,TI22
          complex(16), dimension(NPN1) :: CIM
          complex(16) :: CI,CCI,CCJ
         
          integer(I32P) :: L1MAX,I,I1,J,NMAX1,         &
                           K1,K2,K3,K4,K5,K6,NN,       &
                           M1,M,M1MAX,L1,NN1,NN1MAX,   &
                           N1,NNMAX,NNMIN,KN,N,NNN,M2, &
                           M1MIN,NN1MIN,MMAX,MMIN,  &
                           NL,L
          real(R64P) :: SI,SJ,TT1,TT2,TT3,TT4,TT5,TT6,     &
                        TT7,TT8,T1,T2,T3,T4,SIG,AAR1,AAR2, &
                        AAI1,AAI2,SSS,RR1,RI1,RR2,RI2,     &
                        XR,XI,BBR1,BBR2,BBI1,BBI2,XX,X1,   &
                        X2,X3,X4,X5,X6,X7,X8,DD1,DD2,DD3,  &
                        DD4,DD5R,DD5I,DK,G1L,G2L,G3L,G4L,  &
                        SL,G5LR,G5LI,DM1,DM2,SSS1,FFN,DM3, &
                        DM4,DM5R,DM5I
 
      COMMON /TMAT/ TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22
      COMMON /CBESS/ B1R,B1I,B2R,B2I    
      COMMON /SS/ SSIGN
      EQUIVALENCE ( PLUS1(1),TR11(1,1,1) )
      EQUIVALENCE (D1(1,1,1),PLUS1(1)),                   &     
                 (D2(1,1,1),PLUS1(NPL1*NPN4*NPN4+1)),     &
                 (D3(1,1,1),PLUS1(NPL1*NPN4*NPN4*2+1)),   &
                 (D4(1,1,1),PLUS1(NPL1*NPN4*NPN4*3+1)),   &
                 (D5R(1,1,1),PLUS1(NPL1*NPN4*NPN4*4+1)) 
      
      SSIGN = 0._R64P
      SSI = 0._R64P
      SSJ = 0._R64P
      TR1 = 0._R64P
      TR2 = 0._R64P
      TI1 = 0._R64P
      TI2 = 0._R64P
      G1 = 0._R64P
      G2 = 0._R64P
      AR1 = 0._R64P
      AR2 = 0._R64P
      AI1 = 0._R64P
      AI2 = 0._R64P
      FR = 0._R64P
      FI = 0._R64P
      FF = 0._R64P
      B1R = 0._R64P 
      B1I = 0._R64P 
      B2R = 0._R64P 
      B2I = 0._R64P
      D1 = 0._R64P 
      D2 = 0._R64P 
      D3 = 0._R64P 
      D4 = 0._R64P 
      D5R = 0._R64P 
      D5I = 0._R64P
      PLUS1 = 0._R64P
      TR11 = 0._R64P
      TR12 = 0._R64P
      TR21 = 0._R64P
      TR22 = 0._R64P
      TI11 = 0._R64P
      TI12 = 0._R64P
      TI21 = 0._R64P
      TI22 = 0._R64P
      CIM = DCMPLX(0._R64P,0._R64P)
      CALL FACT
      CALL SIGNUM
      LMAX=2*NMAX
      L1MAX=LMAX+1
      CI=(0D0,1D0)
      CIM(1)=CI
      DO 2 I=2,NMAX
         CIM(I)=CIM(I-1)*CI
    2 CONTINUE
      SSI(1)=1._R64P
      DO 3 I=1,LMAX
         I1=I+1
         SI=DFLOAT(2*I+1)
         SSI(I1)=SI
         IF(I.LE.NMAX) SSJ(I)=DSQRT(SI)
    3 CONTINUE
      CI=-CI
      DO 5 I=1,NMAX
         SI=SSJ(I)
         CCI=CIM(I)
         DO 4 J=1,NMAX
            SJ=1._R64P/SSJ(J)
            CCJ=CIM(J)*SJ/CCI
            FR(J,I)=CCJ
            FI(J,I)=CCJ*CI
            FF(J,I)=SI*SJ
    4    CONTINUE
    5 CONTINUE
      NMAX1=NMAX+1
 
!C *****  CALCULATION OF THE ARRAYS B1 AND B2  *****
 
      K1=1
      K2=0
      K3=0
      K4=1
      K5=1
      K6=2
 
!C     PRINT 3300, B1,B2
 3300 FORMAT (' B1 AND B2')
      DO 100 N=1,NMAX
 
!C *****  CALCULATION OF THE ARRAYS T1 AND T2  *****
 
 
         DO 10 NN=1,NMAX
            M1MAX=MIN0(N,NN)+1
            DO 6 M1=1,M1MAX
               M=M1-1
               L1=NPN6+M
               TT1=TR11(M1,N,NN)
               TT2=TR12(M1,N,NN)
               TT3=TR21(M1,N,NN)
               TT4=TR22(M1,N,NN)
               TT5=TI11(M1,N,NN)
               TT6=TI12(M1,N,NN)
               TT7=TI21(M1,N,NN)
               TT8=TI22(M1,N,NN)
               T1=TT1+TT2
               T2=TT3+TT4
               T3=TT5+TT6
               T4=TT7+TT8
               TR1(L1,NN)=T1+T2
               TR2(L1,NN)=T1-T2
               TI1(L1,NN)=T3+T4
               TI2(L1,NN)=T3-T4
               IF(M.EQ.0) GO TO 6
               L1=NPN6-M
               T1=TT1-TT2
               T2=TT3-TT4
               T3=TT5-TT6
               T4=TT7-TT8
               TR1(L1,NN)=T1-T2
               TR2(L1,NN)=T1+T2
               TI1(L1,NN)=T3-T4
               TI2(L1,NN)=T3+T4
    6       CONTINUE
   10    CONTINUE
 
!C  *****  END OF THE CALCULATION OF THE ARRAYS T1 AND T2  *****
 
         NN1MAX=NMAX1+N
         DO 40 NN1=1,NN1MAX
            N1=NN1-1
 
!C  *****  CALCULATION OF THE ARRAYS A1 AND A2  *****
 
            CALL CCG(N,N1,NMAX,K1,K2,G1)
            NNMAX=MIN0(NMAX,N1+N)
            NNMIN=MAX0(1,IABS(N-N1))
            KN=N+NN1
            DO 15 NN=NNMIN,NNMAX
               NNN=NN+1
               SIG=SSIGN(KN+NN)
               M1MAX=MIN0(N,NN)+NPN6
               AAR1=0._R64P
               AAR2=0._R64P
               AAI1=0._R64P
               AAI2=0._R64P
!DIR$     FMA
               DO 13 M1=NPN6,M1MAX
                  M=M1-NPN6
                  SSS=G1(M1,NNN)
                  RR1=TR1(M1,NN)
                  RI1=TI1(M1,NN)
                  RR2=TR2(M1,NN)
                  RI2=TI2(M1,NN)
                  IF(M.EQ.0) GO TO 12
                  M2=NPN6-M
                  RR1=RR1+TR1(M2,NN)*SIG
                  RI1=RI1+TI1(M2,NN)*SIG
                  RR2=RR2+TR2(M2,NN)*SIG
                  RI2=RI2+TI2(M2,NN)*SIG
   12             AAR1=AAR1+SSS*RR1
                  AAI1=AAI1+SSS*RI1
                  AAR2=AAR2+SSS*RR2
                  AAI2=AAI2+SSS*RI2
   13          CONTINUE
               XR=FR(NN,N)
               XI=FI(NN,N)
               AR1(NN)=AAR1*XR-AAI1*XI
               AI1(NN)=AAR1*XI+AAI1*XR
               AR2(NN)=AAR2*XR-AAI2*XI
               AI2(NN)=AAR2*XI+AAI2*XR
   15       CONTINUE
 
!C  *****  END OF THE CALCULATION OF THE ARRAYS A1 AND A2 ****
 
            CALL CCG(N,N1,NMAX,K3,K4,G2)
            M1=MAX0(-N1+1,-N)
            M2=MIN0(N1+1,N)
            M1MAX=M2+NPN6
            M1MIN=M1+NPN6
            DO 30 M1=M1MIN,M1MAX
               BBR1=0._R64P
               BBI1=0._R64P
               BBR2=0._R64P
               BBI2=0._R64P
!DIR$     FMA
               DO 25 NN=NNMIN,NNMAX
                  NNN=NN+1
                  SSS=G2(M1,NNN)
                  BBR1=BBR1+SSS*AR1(NN)
                  BBI1=BBI1+SSS*AI1(NN)
                  BBR2=BBR2+SSS*AR2(NN)
                  BBI2=BBI2+SSS*AI2(NN)
   25          CONTINUE
               B1R(NN1,M1,N)=BBR1
               B1I(NN1,M1,N)=BBI1
               B2R(NN1,M1,N)=BBR2
               B2I(NN1,M1,N)=BBI2
   30       CONTINUE
   40    CONTINUE
  100 CONTINUE
 
!C  *****  END OF THE CALCULATION OF THE ARRAYS B1 AND B2 ****
 
!C  *****  CALCULATION OF THE ARRAYS D1,D2,D3,D4, AND D5  *****
 
!C     PRINT 3301
 3301 FORMAT(' D1, D2, ...')
      DO 200 N=1,NMAX
         DO 190 NN=1,NMAX
            M1=MIN0(N,NN)
            M1MAX=NPN6+M1
            M1MIN=NPN6-M1
            NN1MAX=NMAX1+MIN0(N,NN)
            DO 180 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD1=0._R64P
               DD2=0._R64P
!DIR$   FMA
               DO 150 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B1R(NN1,M1,NN)
                  X4=B1I(NN1,M1,NN)
                  X5=B2R(NN1,M1,N)
                  X6=B2I(NN1,M1,N)
                  X7=B2R(NN1,M1,NN)
                  X8=B2I(NN1,M1,NN)
                  DD1=DD1+XX*(X1*X3+X2*X4)
                  DD2=DD2+XX*(X5*X7+X6*X8)
  150          CONTINUE
               D1(M1,NN,N)=DD1
               D2(M1,NN,N)=DD2
  180       CONTINUE
            MMAX=MIN0(N,NN+2)
            MMIN=MAX0(-N,-NN+2)
            M1MAX=NPN6+MMAX
            M1MIN=NPN6+MMIN
            DO 186 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD3=0._R64P
               DD4=0._R64P
               DD5R=0._R64P
               DD5I=0._R64P
               M2=-M+2+NPN6
               DO 183 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B2R(NN1,M1,N)
                  X4=B2I(NN1,M1,N)
                  X5=B1R(NN1,M2,NN)
                  X6=B1I(NN1,M2,NN)
                  X7=B2R(NN1,M2,NN)
                  X8=B2I(NN1,M2,NN)
                  DD3=DD3+XX*(X1*X5+X2*X6)
                  DD4=DD4+XX*(X3*X7+X4*X8)
                  DD5R=DD5R+XX*(X3*X5+X4*X6)
                  DD5I=DD5I+XX*(X4*X5-X3*X6)
  183          CONTINUE
               D3(M1,NN,N)=DD3
               D4(M1,NN,N)=DD4
               D5R(M1,NN,N)=DD5R
               D5I(M1,NN,N)=DD5I
  186       CONTINUE
  190    CONTINUE
  200 CONTINUE
 
!C  *****  END OF THE CALCULATION OF THE D-ARRAYS *****
 
!C  *****  CALCULATION OF THE EXPANSION COEFFICIENTS *****
 
!C     PRINT 3303
 3303 FORMAT (' G1, G2, ...')
 
      DK=LAM*LAM/(4._R64P*CSCA*DACOS(-1._R64P))
      DO 300 L1=1,L1MAX
         G1L=0._R64P
         G2L=0._R64P
         G3L=0._R64P
         G4L=0._R64P
         G5LR=0._R64P
         G5LI=0._R64P
         L=L1-1
         SL=SSI(L1)*DK
         DO 290 N=1,NMAX
            NNMIN=MAX0(1,IABS(N-L))
            NNMAX=MIN0(NMAX,N+L)
            IF(NNMAX.LT.NNMIN) GO TO 290
            CALL CCG(N,L,NMAX,K1,K2,G1)
            IF(L.GE.2) CALL CCG(N,L,NMAX,K5,K6,G2)
            NL=N+L
            DO 280  NN=NNMIN,NNMAX
               NNN=NN+1
               MMAX=MIN0(N,NN)
               M1MIN=NPN6-MMAX
               M1MAX=NPN6+MMAX
               SI=SSIGN(NL+NNN)
               DM1=0._R64P
               DM2=0._R64P
!DIR$   FMA
               DO 270 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  IF(M.GE.0) SSS1=G1(M1,NNN)
                  IF(M.LT.0) SSS1=G1(NPN6-M,NNN)*SI
                  DM1=DM1+SSS1*D1(M1,NN,N)
                  DM2=DM2+SSS1*D2(M1,NN,N)
  270          CONTINUE
               FFN=FF(NN,N)
               SSS=G1(NPN6+1,NNN)*FFN
               G1L=G1L+SSS*DM1
               G2L=G2L+SSS*DM2*SI
               IF(L.LT.2) GO TO 280
               DM3=0._R64P
               DM4=0._R64P
               DM5R=0._R64P
               DM5I=0._R64P
               MMAX=MIN0(N,NN+2)
               MMIN=MAX0(-N,-NN+2)
               M1MAX=NPN6+MMAX
               M1MIN=NPN6+MMIN
!DIR$ FMA
               DO 275 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  SSS1=G2(NPN6-M,NNN)
                  DM3=DM3+SSS1*D3(M1,NN,N)
                  DM4=DM4+SSS1*D4(M1,NN,N)
                  DM5R=DM5R+SSS1*D5R(M1,NN,N)
                  DM5I=DM5I+SSS1*D5I(M1,NN,N)
  275          CONTINUE
               G5LR=G5LR-SSS*DM5R
               G5LI=G5LI-SSS*DM5I
               SSS=G2(NPN4,NNN)*FFN
               G3L=G3L+SSS*DM3
               G4L=G4L+SSS*DM4*SI
  280       CONTINUE
  290    CONTINUE
         G1L=G1L*SL
         G2L=G2L*SL
         G3L=G3L*SL
         G4L=G4L*SL
         G5LR=G5LR*SL
         G5LI=G5LI*SL
         ALF1(L1)=G1L+G2L
         ALF2(L1)=G3L+G4L
         ALF3(L1)=G3L-G4L
         ALF4(L1)=G1L-G2L
         BET1(L1)=G5LR*2._R64P
         BET2(L1)=G5LI*2._R64P
         LMAX=L
         IF(DABS(G1L).LT.1.0E-6_R64P) GO TO 500
  300 CONTINUE
  500 CONTINUE
     
    end subroutine
    
!    C   CALCULATION OF THE QUANTITIES F(N+1)=0.5*LN(N!)
!C   0.LE.N.LE.899
 
      subroutine FACT
          implicit none
          real(R64P), dimension(900) :: F
          ! Locals
          integer(I32P) :: I,I1
          COMMON /FAC/ F
!!DIR$     ATTRIBUTES ALIGN : 64 :: /FAC/
          ! Exec code ... 
          F(1)=0._R64P
          F(2)=0._R64P
          DO 2 I=3,900
               I1=I-1
               F(I)=F(I1)+0.5_R64P*DLOG(DFLOAT(I1))
     2 CONTINUE
     
    end subroutine
    
!    C   CALCULATION OF THE ARRAY SSIGN(N+1)=SIGN(N)
!C   0.LE.N.LE.899
 
      subroutine SIGNUM
          implicit none
          real(R64P), dimension(900) :: SSIGN
          ! Locals
          integer(I32P) :: N
          COMMON /SS/ SSIGN
!!DIR$     ATTRIBUTES ALIGN : 64 :: /SS/
          ! Exec code ...
          SSIGN(1)=1._R64P
          DO 2 N=2,899 
               SSIGN(N)=-SSIGN(N-1)
       2 CONTINUE
     
    end subroutine
    
!    C******************************************************************
!C
!C   CALCULATION OF CLEBSCH-GORDAN COEFFICIENTS
!C   (N,M:N1,M1/NN,MM)
!C   FOR GIVEN N AND N1. M1=MM-M, INDEX MM IS FOUND FROM M AS
!C   MM=M*K1+K2
!C
!C   INPUT PARAMETERS :  N,N1,NMAX,K1,K2
!C                               N.LE.NMAX
!C                               N.GE.1
!C                               N1.GE.0
!C                               N1.LE.N+NMAX
!C   OUTPUT PARAMETERS : GG(M+NPN6,NN+1) - ARRAY OF THE CORRESPONDING
!C                                       COEFFICIENTS
!C                               /M/.LE.N
!C                               /M1/=/M*(K1-1)+K2/.LE.N1
!C                               NN.LE.MIN(N+N1,NMAX)
!C                               NN.GE.MAX(/MM/,/N-N1/)
!C   IF K1=1 AND K2=0, THEN 0.LE.M.LE.N
 
 
    subroutine CCG(N,N1,NMAX,K1,K2,GG)
          implicit none
          integer(I32P) :: N,N1,NMAX,K1,K2
          real(R64P), dimension(NPL1,NPN6) :: GG
!DIR$     ASSUME_ALIGNED GG:64
          ! Locals
          real(R64P), dimension(0:NPN5) :: CD,CU
!DIR$     ATTRIBUTES ALIGN : 64 :: CD
!DIR$     ATTRIBUTES ALIGN : 64 :: CU
          integer(I32P) :: NNF,MIN,MF,MIND,MM,M1,NNL,  &
                           NNU,NNM,NN,M
          real(R64P)    :: C,C2,C1,A,B,D
          ! Exec code ... 
          CD = 0._R64P
          CU = 0._R64P
          IF(NMAX.LE.NPN4.AND.0.LE.N1.AND.  &
              N1.LE.NMAX+N.AND.             &
              N.GE.1.AND.N.LE.NMAX) GO TO 1
               PRINT 5001
               STOP
 5001          FORMAT(' ERROR IN SUBROUTINE CCG')
          
     
    1     NNF=MIN0(N+N1,NMAX)
          MIN=NPN6-N
          MF=NPN6+N
          IF(K1.EQ.1.AND.K2.EQ.0) MIN=NPN6
          DO 100 MIND=MIN,MF
              M=MIND-NPN6
              MM=M*K1+K2
              M1=MM-M
              IF(IABS(M1).GT.N1) GO TO 90
                 NNL=MAX0(IABS(MM),IABS(N-N1))
                 IF(NNL.GT.NNF) GO TO 90
                    NNU=N+N1
                    NNM=(NNU+NNL)*0.5_R64P
                    IF (NNU.EQ.NNL) NNM=NNL
                    CALL CCGIN(N,N1,M,MM,C)
                    CU(NNL)=C  
                    IF (NNL.EQ.NNF) GO TO 50
                    C2=0D0
                    C1=C
           DO 7 NN=NNL+1,MIN0(NNM,NNF)
                A=DFLOAT((NN+MM)*(NN-MM)*(N1-N+NN))
                A=A*DFLOAT((N-N1+NN)*(N+N1-NN+1)*(N+N1+NN+1))
                A=DFLOAT(4*NN*NN)/A
                A=A*DFLOAT((2*NN+1)*(2*NN-1))
                A=DSQRT(A)
                B=0.5_R64P*DFLOAT(M-M1)
                D=0._R64P
                IF(NN.EQ.1) GO TO 5
                B=DFLOAT(2*NN*(NN-1))
                B=DFLOAT((2*M-MM)*NN*(NN-1)-MM*N*(N+1)+    &
                    MM*N1*(N1+1))/B
                D=DFLOAT(4*(NN-1)*(NN-1))
                D=D*DFLOAT((2*NN-3)*(2*NN-1))
                D=DFLOAT((NN-MM-1)*(NN+MM-1)*(N1-N+NN-1))/D
                D=D*DFLOAT((N-N1+NN-1)*(N+N1-NN+2)*(N+N1+NN))
                D=DSQRT(D)
    5           C=A*(B*C1-D*C2)
                C2=C1
                C1=C
                CU(NN)=C
    7    CONTINUE
         IF (NNF.LE.NNM) GO TO 50
         CALL DIRECT(N,M,N1,M1,NNU,MM,C)
         CD(NNU)=C
         IF (NNU.EQ.NNM+1) GO TO 50
         C2=0._R64P
         C1=C
         DO 12 NN=NNU-1,NNM+1,-1
            A=DFLOAT((NN-MM+1)*(NN+MM+1)*(N1-N+NN+1))
            A=A*DFLOAT((N-N1+NN+1)*(N+N1-NN)*(N+N1+NN+2))
            A=DFLOAT(4*(NN+1)*(NN+1))/A
            A=A*DFLOAT((2*NN+1)*(2*NN+3))
            A=DSQRT(A)
            B=DFLOAT(2*(NN+2)*(NN+1))
            B=DFLOAT((2*M-MM)*(NN+2)*(NN+1)-MM*N*(N+1)+ &
                                    MM*N1*(N1+1))/B
            D=DFLOAT(4*(NN+2)*(NN+2))
            D=D*DFLOAT((2*NN+5)*(2*NN+3))
            D=DFLOAT((NN+MM+2)*(NN-MM+2)*(N1-N+NN+2))/D
            D=D*DFLOAT((N-N1+NN+2)*(N+N1-NN-1)*(N+N1+NN+3))
            D=DSQRT(D)
            C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CD(NN)=C
   12    CONTINUE
   50    DO 9 NN=NNL,NNF
            IF (NN.LE.NNM) GG(MIND,NN+1)=CU(NN)
            IF (NN.GT.NNM) GG(MIND,NN+1)=CD(NN)
!c           WRITE (6,*) N,M,N1,M1,NN,MM,GG(MIND,NN+1)
    9    CONTINUE
   90    CONTINUE
  100 CONTINUE
     
    end subroutine
    
    subroutine DIRECT (N,M,N1,M1,NN,MM,C)
          implicit none
          integer(I32P) :: N,M,N1,M1,NN,MM
          real(R64P), dimension(900) ::  F
          real(R64P) :: C
          COMMON /FAC/ F
!!DIR$     ATTRIBUTES ALIGN : 64 :: /FAC/
          ! Exec code ...
          F = 0._R64P
          C=F(2*N+1)+F(2*N1+1)+F(N+N1+M+M1+1)+F(N+N1-M-M1+1)    
          C=C-F(2*(N+N1)+1)-F(N+M+1)-F(N-M+1)-F(N1+M1+1)-F(N1-M1+1)
          C=DEXP(C)
     
    end subroutine
    
!    C*********************************************************************
!C
!C   CALCULATION OF THE CLEBCSH-GORDAN COEFFICIENTS
!C   G=(N,M:N1,MM-M/NN,MM)
!C   FOR GIVEN N,N1,M,MM, WHERE NN=MAX(/MM/,/N-N1/)
!C                               /M/.LE.N
!C                               /MM-M/.LE.N1
!C                               /MM/.LE.N+N1
 
    subroutine CCGIN(N,N1,M,MM,G)
          implicit none
          integer(I32P) :: N,N1,M,MM
          real(R64P) :: G
           ! Locals
          real(R64P), dimension(900) :: F,SSIGN
          integer(I32P) :: L1,L2,L3,K,M1,N2,M2,N12,M12
          real(R64P)    :: A
          COMMON /SS/ SSIGN
          COMMON /FAC/ F
!!DIR$     ATTRIBUTES ALIGN : 64 :: /SS/
!!DIR$     ATTRIBUTES ALIGN : 64 :: /FAC/
          F = 0._R64P
          SSIGN = 0._R64P
          M1=MM-M
          IF(N.GE.IABS(M).AND.N1.GE.IABS(M1).AND. &
             IABS(MM).LE.(N+N1)) GO TO 1
          PRINT 5001
          STOP
 5001     FORMAT(' ERROR IN SUBROUTINE CCGIN')
    1     IF (IABS(MM).GT.IABS(N-N1)) GO TO 100
              L1=N
              L2=N1
              L3=M
              IF(N1.LE.N) GO TO 50
              K=N
              N=N1
              N1=K
              K=M
              M=M1
              M1=K
   50         N2=N*2
              M2=M*2
              N12=N1*2
              M12=M1*2
              G=SSIGN(N1+M1+1)* &
                DEXP(F(N+M+1)+F(N-M+1)+F(N12+1)+F(N2-N12+2)-F(N2+2)- &
                     F(N1+M1+1)-F(N1-M1+1)-F(N-N1+MM+1)-F(N-N1-MM+1))
              N=L1
              N1=L2
              M=L3
              RETURN
  100         A=1._R64P
              L1=M
              L2=MM
              IF(MM.GE.0) GO TO 150
              MM=-MM
              M=-M
              M1=-M1
              A=SSIGN(MM+N+N1+1)
  150         G=A*SSIGN(N+M+1)*  &
                DEXP(F(2*MM+2)+F(N+N1-MM+1)+F(N+M+1)+F(N1+M1+1)-        &
                     F(N+N1+MM+2)-F(N-N1+MM+1)-F(-N+N1+MM+1)-F(N-M+1)-  &
                     F(N1-M1+1))
             M=L1
             MM=L2
      
    end subroutine
    
    subroutine SAREA (D,RAT)
          implicit none
          real(R64P) :: D,RAT
          ! Locals
          real(R64P) :: E,R
          IF (D.GE.1._R64P) GO TO 10
          E=DSQRT(1D0-D*D)
          R=0.5_R64P*(D**(2._R64P/3._R64P) + D**(-1._R64P/3._R64P)*DASIN(E)/E)
          R=DSQRT(R)
          RAT=1._R64P/R
          RETURN
   10     E=DSQRT(1._R64P-1._R64P/(D*D))
          R=0.25D0*(2D0*D**(2D0/3D0) + D**(-4D0/3D0)*DLOG((1D0+E)/(1D0-E)) /E)
          R=DSQRT(R)
          RAT=1._R64P/R
    end subroutine  
     
    subroutine SURFCH (N,E,RAT)
          implicit none
          integer(I32P) :: N
          real(R64P)    :: E,RAT
          real(R64P), dimension(60) :: X,W
          real(R64P) :: DN,E2,EN,S,V,XI,DX,DXN,DS,DSN,DCN,  &
                        A,A2,ENS,RS,RV
          integer(I32P) :: NG,I
          X = 0._R64P
          W = 0._R64P
          DN=DFLOAT(N)
          E2=E*E
          EN=E*DN
          NG=60
          CALL GAUSS (NG,0,0,X,W)
          S=0._R64P
          V=0._R64P
          DO 10 I=1,NG
             XI=X(I)
             DX=DACOS(XI)
             DXN=DN*DX
             DS=DSIN(DX)
             DSN=DSIN(DXN)
             DCN=DCOS(DXN)
             A=1._R64P+E*DCN
             A2=A*A
             ENS=EN*DSN
             S=S+W(I)*A*DSQRT(A2+ENS*ENS)
             V=V+W(I)*(DS*A+XI*ENS)*DS*A2
   10 CONTINUE
        RS=DSQRT(S*0.5_R64P)
        RV=(V*3._R64P/4._R64P)**(1._R64P/3._R64P)
        RAT=RV/RS
     
    end subroutine
    
    subroutine SAREAC (EPS,RAT)
          implicit none
          real(R64P) :: EPS,RAT
          RAT=(1.5_R64P/EPS)**(0.333333333333333333333333333_R64P)
          RAT=RAT/DSQRT( (EPS+2._R64P)/(2._R64P*EPS) )
     
    end subroutine
    
!    C  COMPUTATION OF R1 AND R2 FOR A POWER LAW SIZE DISTRIBUTION WITH
!C  EFFECTIVE RADIUS A AND EFFECTIVE VARIANCE B
 
    subroutine POWER (A,B,R1,R2)
          implicit none
          real(R64P) :: A,B,R1,R2
          !external F
          !interface 
          !    double precision function F(R1)
          !      double precision, intent(in) :: R1
          !    end function F
          !end interface
          real(R64P) :: AX,BX,AA,BB
          COMMON /Loc/ AA,BB
          AA=A
          BB=B
          AX=1.0E-5_R64P
          BX=A-1.0E-5_R64P
          R1=ZEROIN (AX,BX,0._R64P)
          R2=(1._R64P+B)*2._R64P*A-R1
     
    end subroutine
    
    real(R64P) function ZEROIN (AX,BX,TOL)
          implicit none
          real(R64P) :: AX,BX,TOL
         ! interface 
         !     double precision function F(x)
         !       double precision, intent(in) :: x
         !     end function F
         ! end interface
          real(R64P) :: TOL1,A,B,FA,FB,C,FC,D,E,XM,  &
                        S,P,Q,R,EPS
          EPS=1D0
       10 EPS=0.5_R64P*EPS
          TOL1=1._R64P+EPS
          IF (TOL1.GT.1D0) GO TO 10
       15 A=AX
          B=BX
          FA=F(A)
          FB=F(B)
       20 C=A
          FC=FA
          D=B-A
          E=D
       30 IF (DABS(FC).GE.DABS(FB)) GO TO 40
       35 A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
       40 TOL1=2.0_R64P*EPS*DABS(B)+0.5_R64P*TOL
          XM=0.5_R64P*(C-B)
          IF (DABS(XM).LE.TOL1) GO TO 90
       44 IF (FB.EQ.0._R64P) GO TO 90
       45 IF (DABS(E).LT.TOL1) GO TO 70
       46 IF (DABS(FA).LE.DABS(FB)) GO TO 70
       47 IF (A.NE.C) GO TO 50
       48 S=FB/FA
          P=2._R64P*XM*S
          Q=1._R64P-S
          GO TO 60
       50 Q=FA/FC
          R=FB/FC
          S=FB/FA
          P=S*(2.0_R64P*XM*Q*(Q-R)-(B-A)*(R-1.0_R64P))
          Q=(Q-1.0_R64P)*(R-1.0_R64P)*(S-1.0_R64P)
       60 IF (P.GT.0.0_R64P) Q=-Q
          P=DABS(P)
          IF ((2.0_R64P*P).GE.(3.0_R64P*XM*Q-DABS(TOL1*Q))) GO TO 70
       64 IF (P.GE.DABS(0.5_R64P*E*Q)) GO TO 70
       65 E=D
          D=P/Q
          GO TO 80
       70 D=XM
          E=D
       80 A=B
          FA=FB
          IF (DABS(D).GT.TOL1) B=B+D
          IF (DABS(D).LE.TOL1) B=B+DSIGN(TOL1,XM)
          FB=F(B)
          IF ((FB*(FC/DABS(FC))).GT.0D0) GO TO 20
       85 GO TO 30
       90 ZEROIN=B
       
    end function
    
    real(R64P) function F(R1)
          implicit none
          real(R64P) :: R1
          real(R64P) :: R2 ,A,B
          COMMON /Loc/ A,B
          R2=(1.0_R64P+B)*2.0_R64P*A-R1
          F=(R2-R1)/DLOG(R2/R1)-A
      
    end function
    
!    C**********************************************************************
!C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         *
!C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      *
!C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.               *
!C    N - NUMBER OF POINTS                                             *
!C    Z - DIVISION POINTS                                              *
!C    W - WEIGHTS                                                      *
!C**********************************************************************
 
    subroutine  GAUSS ( N,IND1,IND2,Z,W )
          implicit none
          integer(I32P) :: N,IND1,IND2
          real(R64P), dimension(N) :: Z,W
!DIR$     ASSUME_ALIGNED Z:64,W:64
          integer(I32P) :: IND,K,M,I,NITER,J
          real(R64P)    :: A,B,C,F,X,PB,CHECK,PC,DJ,  &
                           PA,ZZ
          A=1._R64P
          B=2._R64P
          C=3._R64P
          IND=MOD(N,2)
          K=N/2+IND
          F=DFLOAT(N)
          DO 100 I=1,K
                 M=N+1-I
                 IF(I.EQ.1) X=A-B/((F+A)*F)
                 IF(I.EQ.2) X=(Z(N)-A)*4.0_R64P+Z(N)
                 IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6_R64P+Z(N-1)
                 IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
                 IF(I.EQ.K.AND.IND.EQ.1) X=0._R64P
                 NITER=0
                 CHECK=1D-16
   10            PB=1._R64P
                 NITER=NITER+1
                 IF (NITER.LE.100) GO TO 15
                 CHECK=CHECK*10D0
   15            PC=X
                 DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
              PA=A/((PB-X*PC)*F)
              PB=PA*PC*(A-X*X)
              X=X-PB
              IF(DABS(PB).GT.CHECK*DABS(X)) GO TO 10
              Z(M)=X
              W(M)=PA*PA*(A-X*X)
              IF(IND1.EQ.0) W(M)=B*W(M)
              IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
              Z(I)=-Z(M)
              W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA', ' OF ',I4,'-TH ORDER')
      
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
!C     PRINT 1300,N
 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
     
    end subroutine
    
    subroutine DISTRB (NNK,YY,WY,NDISTR,AA,BB,GAM,R1,R2,REFF,    &             
                        VEFF,P1)                                               
          implicit none
          integer(I32P) :: NNK,NDISTR
          real(R64P), dimension(NNK) ::YY,WY
!DIR$     ASSUME_ALIGNED YY:64,WW:64
          real(R64P) :: AA,BB,GAM,R1,R2,REFF,VEFF,P1
          ! Locals
          integer(I32P) :: I
          real(R64P)    :: A2,B2,DB,X,Y,DA,DAB,SUM,G,XI
          IF (NDISTR.EQ.2) GO TO 100                                                
          IF (NDISTR.EQ.3) GO TO 200                                                
          IF (NDISTR.EQ.4) GO TO 300                                                
          IF (NDISTR.EQ.5) GO TO 360
          PRINT 1001,AA,BB,GAM                                                      
 1001 FORMAT('MODIFIED GAMMA DISTRIBUTION, alpha=',F6.4,'  r_c=',F6.4,'  gamma=',F6.4)                  
                                                     
          A2=AA/GAM                                                                 
          DB=1D0/BB
       DO 50 I=1,NNK                                                             
             X=YY(I)                                                             
             Y=X**AA                                                                
             X=X*DB
             Y=Y*DEXP(-A2*(X**GAM))                                                 
             WY(I)=WY(I)*Y                                                       
    50 CONTINUE                                                                  
       GO TO 400                                                                 
  100 PRINT 1002,AA,BB                                                          
 1002 FORMAT('LOG-NORMAL DISTRIBUTION, r_g=',F8.4,'  [ln(sigma_g)]**2=', F6.4)                   
            
      DA=1.0_R64P/AA                                                                 
      DO 150 I=1,NNK                                                            
             X=YY(I)                                                                
             Y=DLOG(X*DA)                                                          
             Y=DEXP(-Y*Y*0.5_R64P/BB)/X                                             
             WY(I)=WY(I)*Y                                                          
  150 CONTINUE                                                                  
      GO TO 400                                                                 
  200 PRINT 1003                                                                
 1003 FORMAT('POWER LAW DISTRIBUTION OF HANSEN & TRAVIS 1974')                 
      DO 250 I=1,NNK                                                            
             X=YY(I)                                                                
             WY(I)=WY(I)/(X*X*X)                                                 
  250 CONTINUE                                                                  
      GO TO 400                                                                 
  300 PRINT 1004,AA,BB                                                          
 1004 FORMAT ('GAMMA DISTRIBUTION,  a=',F6.3,'  b=',F6.4)
      B2=(1.0_R64P-3.0_R64P*BB)/BB                                                        
      DAB=1.0_R64P/(AA*BB)                                                          
      DO 350 I=1,NNK                                                            
             X=YY(I)                                                                
             X=(X**B2)*DEXP(-X*DAB)                                                 
             WY(I)=WY(I)*X                                                       
  350 CONTINUE                                                                  
      GO TO 400                                                                 
  360 PRINT 1005,BB
 1005 FORMAT ('MODIFIED POWER LAW DISTRIBUTION,  alpha=',D10.4)
      DO 370 I=1,NNK
             X=YY(I)
             IF (X.LE.R1) WY(I)=WY(I)
             IF (X.GT.R1) WY(I)=WY(I)*(X/R1)**BB
  370 CONTINUE
  400 CONTINUE                                                                  
      SUM=0._R64P
!DIR$ SIMD REDUCTION(+,SUM)
      DO 450 I=1,NNK
             SUM=SUM+WY(I)
  450 CONTINUE
      SUM=1.0_R64P/SUM
!DIR$ VECTOR ALIGNED
!DIR$ SIMD VECTORLENGTHFOR(REAL(KIND=8))
      DO 500 I=1,NNK
             WY(I)=WY(I)*SUM
  500 CONTINUE
      G=0._R64P
      DO 550 I=1,NNK
              X=YY(I)
              G=G+X*X*WY(I)
  550 CONTINUE
      REFF=0.0_R64P
      DO 600 I=1,NNK
         X=YY(I)
         REFF=REFF+X*X*X*WY(I)
  600 CONTINUE
      REFF=REFF/G
      VEFF=0._R64P
      DO 650 I=1,NNK
         X=YY(I)
         XI=X-REFF
         VEFF=VEFF+XI*XI*X*X*WY(I)
  650 CONTINUE
      VEFF=VEFF/(G*REFF*REFF)
                                                                     
    end subroutine
    
    subroutine HOVENR(L1,A1,A2,A3,A4,B1,B2)
          implicit none
          integer(I32P) :: L1
          real(R64P), dimension(L1) :: A1,A2,A3,A4,B1,B2
!DIR$     ASSUME_ALIGNED A1:64,A2:64,A3:64,A4:64,B1:64,B2:64
          integer(I32P) :: L,KONTR,LL,I
          real(R64P)    :: DL,DDL,AA1,AA2,AA3,AA4,BB1,BB2,C,  &
                           CC,C1,C2,C3
          DO 100 L=1,L1
                 KONTR=1
                 LL=L-1
                 DL=DFLOAT(LL)*2.0_R64P+1.0_R64P
                 DDL=DL*0.48_R64P
                 AA1=A1(L)
                 AA2=A2(L)
                 AA3=A3(L)
                 AA4=A4(L)
                 BB1=B1(L)
                 BB2=B2(L)
                 IF(LL.GE.1.AND.DABS(AA1).GE.DL) KONTR=2
                 IF(DABS(AA2).GE.DL) KONTR=2
                 IF(DABS(AA3).GE.DL) KONTR=2
                 IF(DABS(AA4).GE.DL) KONTR=2
                 IF(DABS(BB1).GE.DDL) KONTR=2
                 IF(DABS(BB2).GE.DDL) KONTR=2
                 IF(KONTR.EQ.2) PRINT 3000,LL
                 C=-0.1_R64P
                 DO 50 I=1,11
                       C=C+0.1D0
                       CC=C*C
                       C1=CC*BB2*BB2
                       C2=C*AA4
                       C3=C*AA3
                       IF((DL-C*AA1)*(DL-C*AA2)-CC*BB1*BB1.LE.-1D-4) KONTR=2
                       IF((DL-C2)*(DL-C3)+C1.LE.-1D-4) KONTR=2
                       IF((DL+C2)*(DL-C3)-C1.LE.-1D-4) KONTR=2
                       IF((DL-C2)*(DL+C3)-C1.LE.-1D-4) KONTR=2
                       IF(KONTR.EQ.2) PRINT 4000,LL,C
          50    CONTINUE
     100 CONTINUE
      IF(KONTR.EQ.1) PRINT 2000
 2000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS SATISFIED')
 3000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=',I3)
 4000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=',I3,'   A=',D9.2)
    end subroutine
    
!    C****************************************************************
 
!C    CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION
!C    COEFFICIENTS
! 
!C    A1,...,B2 - EXPANSION COEFFICIENTS
!C    LMAX - NUMBER OF COEFFICIENTS MINUS 1
!C    N - NUMBER OF SCATTERING ANGLES
!C        THE CORRESPONDING SCATTERING ANGLES ARE GIVEN BY
!C        180*(I-1)/(N-1) (DEGREES), WHERE I NUMBERS THE ANGLES
 
    SUBROUTINE MATR(A1,A2,A3,A4,B1,B2,scatmat,LMAX,NPNA)
          implicit none
          real(R64P), dimension(NPL) :: A1,A2,A3,A4,B1,B2
!DIR$     ASSUME_ALIGNED A1:64,A2:64,A3:64,A4:64,B1:64,B2:64
          real(R64P), dimension(6,NPNA) :: scatmat
          integer(I32P) :: LMAX,NPNA
          ! Locals
          integer(I32P) :: N,L1MAX,I1,L1,L
          real(R64P)    :: DN,DA,DB,TB,TAA,U,F11,F12,F33,F34,F22,F44, &
                           P1,P2,P3,P4,PP1,PP2,PP3,PP4,D6,DL,DL1,PL1, &
                           P,F2,F3,PL2,PL3,PL4
          ! Exec code ...                 
          N=NPNA
          DN=1.0_R64P/DFLOAT(N-1)
          DA=PI*DN
          DB=180.0_R64P*DN
          L1MAX=LMAX+1
!      PRINT 1000
! 1000 FORMAT(' ')
!      PRINT 1001
! 1001 FORMAT(' ',2X,'S',6X,'ALPHA1',6X,'ALPHA2',6X,'ALPHA3',
!     &       6X,'ALPHA4',7X,'BETA1',7X,'BETA2')
!      DO 10 L1=1,L1MAX
 !        L=L1-1
!         PRINT 1002,L,A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
!   10 CONTINUE
 !1002 FORMAT(' ',I3,6F12.5)
          TB=-DB
          TAA=-DA
!          PRINT 1000
!      PRINT 1003
! 1003 FORMAT(' ',5X,'<',8X,'F11',8X,'F22',8X,'F33',
!     & 8X,'F44',8X,'F12',8X,'F34')
          D6=DSQRT(6.0_R64P)*0.25_R64P
      DO 500 I1=1,N
             TAA=TAA+DA
             TB=TB+DB
             U=DCOS(TAA)
             F11=0._R64P
             F2=0._R64P
             F3=0._R64P
             F44=0._R64P
             F12=0._R64P
             F34=0._R64P
             P1=0._R64P
             P2=0._R64P
             P3=0._R64P
             P4=0._R64P
             PP1=1._R64P
             PP2=0.25_R64P*(1.0_R64P+U)*(1.0_R64P+U)
             PP3=0.25_R64P*(1.0_R64P-U)*(1.0_R64P-U)
             PP4=D6*(U*U-1.0_R64P)
         DO 400 L1=1,L1MAX
            L=L1-1
            DL=DFLOAT(L)
            DL1=DFLOAT(L1)
            F11=F11+A1(L1)*PP1
            F44=F44+A4(L1)*PP1
            IF(L.EQ.LMAX) GO TO 350
            PL1=DFLOAT(2*L+1)
            P=(PL1*U*PP1-DL*P1)/DL1
            P1=PP1
            PP1=P
  350       IF(L.LT.2) GO TO 400
            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            F12=F12+B1(L1)*PP4
            F34=F34+B2(L1)*PP4
            IF(L.EQ.LMAX) GO TO 400
            PL2=DFLOAT(L*L1)*U
            PL3=DFLOAT(L1*(L*L-4))
            PL4=1D0/DFLOAT(L*(L1*L1-4))
            P=(PL1*(PL2-4.0_R64P)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4._R64P)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DFLOAT(L*L-4))*P4)/DSQRT(DFLOAT(L1*L1-4))
            P4=PP4
            PP4=P
  400    CONTINUE
         F22=(F2+F3)*0.5_R64P
         F33=(F2-F3)*0.5_R64P
!C        F22=F22/F11
!C        F33=F33/F11
!C        F44=F44/F11
!C        F12=-F12/F11
!C        F34=F34/F11
         scatmat(1,I1) = F11
         scatmat(2,I1) = F12
         scatmat(3,I1) = F33
         scatmat(4,I1) = F34
         scatmat(5,I1) = F22
         scatmat(6,I1) = F44
!         PRINT 1004,TB,F11,F22,F33,F44,F12,F34
  500 CONTINUE
  !    PRINT 1000 
 !1004 FORMAT(' ',F6.2,6F11.4)
      
    end subroutine
    
    
    
     
    
    
    
    
 
 
    
    
    
    
    
                     

end module mod_tmatrix