
include 'GMS_config.fpp'

module  mod_tmatrix_mps

 !===================================================================================85
 !---------------------------- DESCRIPTION ------------------------------------------85
 !
 !
 !
 !          Module  name:
 !                         'mod_tmatrix_mps'
 !          
 !          Purpose:
 !                       This module is a modified version of  Yu-lin Xu 'gmm03Trd.f'
 !                       program.
 !          History:
 !                        Date: 05-06-2018
 !                        Time: 18:46 GMT+2
 !                        Date: 27-08-2019 (Modified)
  !                        Time: 19:12 GMT+2
  !                       Date: 12-09-2019
  !                       Time: 17:32 GMT+2
 !          Version:
 !
 !                      Major: 2
 !                      Minor: 0
 !                      Micro: 0
 !
 !          Author:  
 !                        The code has been developed by Yu-lin Xu.
 !                        The subroutine tm0d.f together with its auxiliaries in 
 !                        this gmm03Trd.f code are a part of Mishchenko's public-domain code
 !                        ampld.new.f available online at http://www.giss.nasa.gov/~crmim.
 !
 !                        For questions/comments/suggestions/bugs/problems please contact 
!C                        Yu-lin Xu at yu-lin.xu1@nasa.gov
 !          Modified:
 !                   Bernard Gingold on 05-06-2018 (@see history)
 !                 
 !          References:
 !         
 !                         (1) The multiparticle-scattering formulation used in this code can be 
!C                             found in the papers
!C                               Bruning and Lo, IEEE Trans. Anten. Prop. AP-19, 378 (1971)
!C                               Xu, Appl. Opt. 34, 4573 (1995)  
!C                               Appl. Opt. 36, 9496 (1997) 
!C                               Phys. Lett. A 249, 30 (1998)
!C                               Phys. Rev. E 67, 046620 (2003)
!C                               J. Opt. Soc. Am. A 20, 2093 (2003)
!C                               Xu and Wang, Phys. Rev. E 58, 3931 (1998)
!C                               Xu, Gustafson, Giovane, Blum, and Tehranian,
!C                               Phys. Rev. E 60, 2347 (1999)
!C                               Xu and Gustafson, JQSRT 70, 395 (2001) 
!C                               Xu and Khlebtsov, JQSRT 79-80, 1121 (2003)
!C                               Xu, J. Opt. Soc. Am. A 20, 2093 (2003)
!C                               Xu, Phys. Rev. E 67, 046620 (2003)
!C                               Xu and Gustafson, in "Recent Research Development in Optics" 
!C                               (Kerala, India: Research Signpost, 2003), pp.599-648
!C                               Xu, JQSRT 89, 385 (2004)
!C                          (2)  Numerical techniques used in this code can be found in the papers
!C                               Stein, Q. Appl. Math. 19, 15 (1961)
!C                               Cruzan, Q. Appl. Math. 20, 33 (1962)
!C                               Mackowski, Proc. R. Soc. Lond. A 433, 599 (1991)
!C                               Wang and van de Hulst, Appl. Opt. 30, 106 (1991)
!C                               H.A. van der Vorst, SIAM J. Sci. Stat. Comput. 13, 631 (1992)
!C                               Gutknecht, SIAM J. Sci. Comput. 14, 1020 (1993)
!C                               Xu, J. Comput. Appl. Math. 85, 53 (1997)
!C                               J. Comput. Phys. 139, 137 (1998)
!C                            (3) In the input of this code individual particles can be an 
!C                                arbitrary mixture of homogeneous spheres, core-mantle spheres, 
!C                                and certain types of axially symmetric shapes including spheroids, 
!C                                finite circular cylinders, Chebyshev particles, and generalized 
!C                                Chebyshev particles (simulating the shapes of distorted water 
!C                                drops). The subroutine tm0d.f together with its auxiliaries in 
!C                                this gmm03Trd.f code are a part of Mishchenko's public-domain code
!C                                ampld.new.f available online at http://www.giss.nasa.gov/~crmim. 
!C                                It computes the T-matrix of an individual nonspherical particle 
!C                                with an axially symmetric shape in single-body scattering, based 
!C                                on Waterman's extended boundary condition method (or called the 
!C                                null field method).                                 
 !    
 !         
 !          E-mail:
 !                  
 !                      beniekg@gmail.com
 !==================================================================================85
    ! Tab:5 col - Type and etc.. definitions
    ! Tab:10,11 col - Type , function and subroutine code blocks.
        
   
    use mod_kinds,    only : int4, sp, dp
#if defined __INTEL_COMPILER
    use IFCORE,       only : TRACEBACKQQ
    use IFPORT,       only : SYSTEMQQ,GETLASTERRORQQ,DCLOCK
#endif
    use mod_lnf
    implicit none
    !=====================================================59
    !  File and module information:
    !  version,creation and build date, author,description
    !=====================================================59
    
    ! Major version
    integer(kind=int4), parameter, public :: MOD_TMATRIX_MPS_MAJOR = 2
    
    ! Minor version
    integer(kind=int4), parameter, public :: MOD_TMATRIX_MPS_MINOR = 0
    
    ! Micro version
    integer(kind=int4), parameter, public :: MOD_TMATRIX_MPS_MICRO = 0
    
    ! Module full version
    integer(kind=int4), parameter, public :: MOD_TMATRIX_MPS_FULLVER = 1000*MOD_TMATRIX_MPS_MAJOR + &
                                                                  100*MOD_TMATRIX_MPS_MINOR  + &
                                                                  10*MOD_TMATRIX_MPS_MICRO
    
    ! Module creation date
    character(*),  parameter, public :: MOD_TMATRIX_MPS_CREATE_DATE = "05-06-2018 18:46 +00200 (TUE 05 JUN 2018 GMT+2)"
    
    ! Module build date ( should be set after successful compilation)
    character(*),  parameter, public :: MOD_TMATRIX_MPS_BUILD_DATE = __DATE__ " " __TIME__
    
    ! Module author info
    character(*),  parameter, public :: MOD_TMATRIX_MPS_AUTHOR = "Programmer: Yu-lin Xu, yu-lin.xu1@nasa.gov, modified by Bernard Gingold, beniekg@gmail.com"
    
    ! Module short description
    character(*),  parameter, public :: MOD_TMATRIX_MPS_DESCRIPT = " Calculation of random- orientation averaged radiative-scattering by a single or an ensemble of variously shaped small particles" 
                         
    ! Constants
    integer(kind=int4), parameter, private :: NPN1  = 100
    
    integer(kind=int4), parameter, private :: NPNG1 = 500
    
    integer(kind=int4), parameter, private :: NPNG2 = 2*NPNG1
    
    integer(kind=int4), parameter, private :: NPN2  =  2*NPN1
    
    integer(kind=int4), parameter, private :: NPL   = NPN2+1
    
    integer(kind=int4), parameter, private :: NPN3  = NPN1+1
    
    integer(kind=int4), parameter, private :: NPN4  = NPN1
    
    integer(kind=int4), parameter, private :: NPN5  = 2*NPN4
    
    integer(kind=int4), parameter, private :: NPN6  = NPN4+1
    
    integer(kind=int4), parameter, public  :: NANGMAX = 1808
    
    complex(16),   parameter, private :: CZERO = DCMPLX(0._dp,0._dp)
    
    real(kind=dp),    parameter, private :: ZERO = 0._dp
    
    contains
    
      subroutine tmatrix_mps_driver(analysis_type,pmc_loop744,pmc_loop965,pmc_loop1097,pmc_loop1179, &
                                    pmc_loop1778,pmc_loop1828,pmc_loop2039,pmc_loop2092,pmc_loop2323, &
                                    idMie,small,MXINT,NADD,idscmt,sang,w,irat, &
                                    nL,idshp,shp,r0,cext,cabs,csca,assym,cextv,cabsv, &
                                    cscav,cbakv,cprv,cexts,cabss,cscas,cbaks,cprs, &
                                    dang,inat,pol,i11,i21,i12,i22,cexti,cabsi,cscai, &
                                    assymi,cpri,mue                                    )
    
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: tmatrix_mps_driver
          include 'tmatrix_mps_np.inc'
          integer(kind=int4)                   :: analysis_type
          integer(kind=int4)                   :: pmc_loop744,pmc_loop965,pmc_loop1097,pmc_loop1179, &
                                                  pmc_loop1778,pmc_loop1828,pmc_loop2039,pmc_loop2092, &
                                                  pmc_loop2323
          integer(kind=int4)                   :: idMie
          real(kind=dp)                        :: small
          integer(kind=int4)                   :: MXINT,NADD,idscmt
          real(kind=dp)                        :: sang,w
          integer(kind=int4)                   :: irat,nL
          integer(kind=int4), dimension(nLp)   :: idshp
          real(kind=dp),    dimension(3,nLp)   :: shp
          real(kind=dp),    dimension(9,nLp)   :: r0
          real(kind=dp)                        :: cext,cabs,csca,assym,cextv,cabsv,  &
                                                  cscav,cbakv,cprv,cexts,cabss,cscas,&
                                                  cbaks,cprs
          real(kind=dp),    dimension(NANGMAX) :: dang,inat,pol,i11,i21,i12,i22
          real(kind=dp),    dimension(nLp)     :: cexti,cabsi,cscai,assymi,cpri
          real(kind=dp),    dimension(4,4,NANGMAX) :: mue
          ! Locals
          integer(kind=int4), parameter :: nmp   =   np*(np+2)       !
          integer(kind=int4), parameter :: nmp0  =   (np+1)*(np+4)/2 
          integer(kind=int4), parameter :: np2   =   2*np
          integer(kind=int4), parameter :: ni0   =   np*(np+1)*(np2+1)/3+np*np
          integer(kind=int4), parameter :: ng0   =   np*(2*np**3+10*np**2+19*np+5)/6
          integer(kind=int4), parameter :: nrc   =   4*np*(np+1)*(np+2)/3+np
          integer(kind=int4), parameter :: nij   =   nLp*(nLp-1)/2
          real(kind=dp), automatic :: start,end,duration
          integer(kind=int4) :: u,v,u0
          logical(kind=int4), automatic :: result
          integer(kind=int4), automatic :: ret
          integer(kind=int4), automatic :: line
          integer(kind=int4), automatic :: aerr
          character(len=256), automatic :: emsg
          integer(kind=int4), dimension(nLp) :: nmax,uvmax,ind
          real(kind=dp), dimension(nLp), automatic :: x,xc
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: x,xc
#endif
          real(kind=dp), dimension(3,nLp) :: r00
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: R00
#endif
          real(kind=dp), dimension(0:np2+1) :: besj,besy
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: besj,besy
#endif
          real(kind=dp), allocatable, dimension(:,:) :: drot ! dynamic alloc  ! real(kind=dp), dimension(nrc,nij)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: drot
#endif
          real(kind=dp), dimension(nLp), automatic     :: c0i,c1i
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: c0i,c1i
#endif
          real(kind=dp), dimension(nij)     :: bes0
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: bes0
#endif
          real(kind=dp), dimension(5,nij)   :: confg
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: confg
#endif
          real(kind=dp), dimension(2)      :: taup,taupj,taupg,taupjg,tau0p, &
                                            tau1p,tau2p,tau0pg,tau1pg,     &
                                            tau2pg,tau0pj,tau1pj,tau2pj,   &
                                            tau0pjg,tau1pjg,tau2pjg
          real(kind=dp), dimension(2,2)    :: tau20,tau11,tau02,tau20g,tau11g, &
                                            tau02g,taum,taumg
          real(kind=dp), dimension(np2+1)   :: w01s,wcf
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: w01s,wcf
#endif
          real(kind=dp), dimension(0:4*(np+1)) :: fnr
          real(kind=dp), allocatable, dimension(:,:,:,:) :: wmf1,wm1,wsdt !dynamic alloc   real(kind=dp), dimension(np,np,0:np,0:np2)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: wmf1,wm1,wsdt
#endif
          real(kind=dp), dimension(0:np+2)       :: bcof
          real(kind=dp), dimension(-np:np,0:nmp) ::  dc
          real(kind=dp), dimension(nmp0) :: pi,tau
          integer(kind=int4), dimension(ni0) :: iga0
          real(kind=dp),    dimension(ng0) :: ga0
          real(kind=dp),    dimension(ni0) :: cof0
          real(kind=dp),    dimension(nmp) :: cofsr
          real(kind=dp), dimension(NPN6,NPN4,NPN4)   :: RT11,RT12,RT21,RT22, &
                                                     IT11,IT12,IT21,IT22
          complex(16), dimension(2,2,2,2)         :: A0p,A1p,B0p,B1p,A0pg,&
                                                     A1pg,B0pg,B1pg
          complex(16), allocatable, dimension(:,:)         :: atr0,btr0  ! dynamic alloc  complex(16), dimension(ni0,nij)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: atr0,btr0
#endif
          complex(16), dimension(2,np,nmp)        :: atr

          complex(16), dimension(nmp)            :: at,bt
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: at,bt
#endif
          complex(16), allocatable, dimension(:,:)         :: atr1,btr1  ! dynamic alloc  complex(16), dimension(ni0,nij)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: atr1,btr1
#endif
          complex(16), allocatable, dimension(:,:)          :: ek ! dynamic alloc   complex(16), dimension(np,nij)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 ::  ek
#endif
          complex(16), dimension(nLp)           :: ref,refc
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 ::  ref,refc
#endif
          complex(16), allocatable, dimension(:,:)         :: p0,q0   !  complex(16), dimension(nLp,nmp)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 ::  p0,q0
#endif
          complex(16), dimension(np)            :: an,bn
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 ::  an,bn
#endif
          complex(16), dimension(nLp)           :: B2i
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 ::  B2i
#endif
          complex(16), allocatable, dimension(:,:)         :: at0,bt0 !  complex(16), dimension(nLp,nmp)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 ::  at0,bt0
#endif
          complex(16), dimension(nmp)             :: at1,bt1
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 ::  at1,bt1
#endif
          complex(16), allocatable, dimension(:,:)         :: as,bs,as2,bs2  !  complex(16), dimension(nLp,nmp)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: as,bs,as2,bs2
#endif
          complex(16), allocatable, dimension(:,:,:)        :: tta,ttb,tta0,ttb0,asr,bsr, &  ! dynamic alloc   complex(16), dimension(nLp,nLp,nmp)
                                                        as0,bs0,asc,bsc,as1,bs1,   &
                                                        ast,bst,asp,bsp,asv,bsv
#if defined __INTEL_COMPILER
!DIR$     ATTRIBUTES ALIGN : 64 :: tta,ttb,tta0,ttb0,asr,bsr,as0,bs0,asc,bsc
          !DIR$     ATTRIBUTES ALIGN : 64 :: as1,bs1,ast,bst,asp,bsp,asv,bsv
#endif
          complex(16), dimension(nmp)                 :: atj,btj
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: atj,btj
#endif
          complex(16), allocatable, dimension(:,:,:,:,:,:) :: pct ! dynamic alloc   complex(16), dimension(nLp,nLp,2,2,nmp,nmp)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: pct
#endif
          complex(16), allocatable, dimension(:,:,:,:,:)      :: A1m,A2m  ! dynamic alloc   complex(16), dimension(0:np,np,np,2,2)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: A1m,A2m
#endif
          complex(16), allocatable, dimension(:,:,:,:)           :: B11n,B12n,B13n,B21n,B22n,B23n  !  complex(16), dimension(np,np,2,2)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: B11n,B12n,B13n,B21n,B22n,B23n
#endif
          complex(16), allocatable, dimension(:,:,:,:)    ::   fhmf1,fmf1,fhm1v,fm1v,fhm1q, &   ! dynamic alloc   complex(16), dimension(0:np2,np,0:np2,2) 
               fm1q,fhmf1vq,fmf1vq,fhmf4vq
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 ::  fhmf1,fmf1,fhm1v,fm1v,fhm1q,fm1q,fhmf4vq,fmf1vq
#endif
          complex(16), dimension(0:np2,2,2)           ::   fhas,fnhs
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: fhas,fnhs
#endif
          complex(16), allocatable, dimension(:,:,:,:,:)     ::   tbar  ! dynamic alloc    complex(16), dimension(nLp,2,2,nmp,nmp)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: tbar
#endif
          complex(16), allocatable, dimension(:,:,:,:)         ::   tbar0 ! dynamic alloc   complex(16), dimension(2,2,nmp,nmp)
#if defined __INTEL_COMPILER
          !DIR$     ATTRIBUTES ALIGN : 64 :: tbar0
#endif
          complex(16), dimension(2,2)      ::   bar
          complex(16), dimension(-2*np:2*np)          :: ekt
          
          
        
#if defined __INTEL_COMPILER              
               !DIR$ ATTRIBUTES ALIGN : 8 :: k,pione,gcs,gcv,eps
               !DIR$ ATTRIBUTES ALIGN : 8 :: fint,temp,temp0,x0,y0,z0,gcsr
               !DIR$ ATTRIBUTES ALIGN : 8 :: gcvr,xv,xs,ratio,RAT,DDELT,alph
               !DIR$ ATTRIBUTES ALIGN : 8 :: beta, s,t,ca,sa
               !DIR$ ATTRIBUTES ALIGN : 8 :: sb,cb,xd,d,cz,cext0,cext1,fuv1,fuv2
               !DIR$ ATTRIBUTES ALIGN : 8 :: fuv3,guv,guv3,guv5,guv4,guv6,rn,rm
               !DIR$ ATTRIBUTES ALIGN : 8 :: p,fnp,fn,gmn,gmn1,gmn2,gmn3,gmn5,gmn6
               !DIR$ ATTRIBUTES ALIGN : 8 ::  gmn4,cwmf1,theta,gn,cbak,lnfacd,gmnj
               !DIR$ ATTRIBUTES ALIGN : 8 :: guv1,guv2,sphi,cphi,cpr,fv,fv0,cv,cwm1,cq,gt,xt,cscax,cextx,cprx,cabsx,assymx,cbakx
#endif              
               real(kind=dp) :: k,pih,twopi,pione,gcs,gcv,eps,fint,  &
                         temp,temp0,x0,y0,z0,gcsr,gcvr,xv,xs, &
                        ratio,RAT,DDELT,alph,beta,s,t,ca,sa, &
                        sb,cb,xd,d,cz,cext0,cext1,fuv1,fuv2, &
                        fuv3,guv,guv3,guv5,guv4,guv6,rn,rm,  &
                        p,fnp,fn,gmn,gmn1,gmn2,gmn3,gmn5,gmn6,    &
                        gmn4,cwmf1,theta,gn,cbak,lnfacd,gmnj,&
                        guv1,guv2,xt,sphi,cphi,cpr,fv,fv0,cv, &
                        cwm1,cq,gt,cscax,cextx,cprx,cabsx,assymx,cbakx
               
         
          
          
#if defined __INTEL_COMPILER             
             !DIR$ ATTRIBUTES ALIGN : 4 ::  i,j1,j2,n,in0,iv0,m,imn,iuv,imn1
             !DIR$ ATTRIBUTES ALIGN : 4 ::  iuv1,nmax0,imax,n0,ii,ij,iv,inn,n1          !! in original code was: in, here it was changed to: inn
             !DIR$ ATTRIBUTES ALIGN : 4 ::  ip,iq,is,isn,isv,nlarge,nbes,irc      
             !DIR$ ATTRIBUTES ALIGN : 4 ::  itrc,iuvc,niter,ijmax,ijmin,iuvp        
             !DIR$ ATTRIBUTES ALIGN : 4 ::  jv1,iuv2,iuv3,juv1,juv2,juv3,juv4       
             !DIR$ ATTRIBUTES ALIGN : 4 ::  juv5,juv6,ntemp,n2,imn2,imn3,jmn1,jmn2
             !DIR$ ATTRIBUTES ALIGN : 4 ::  jmn3,jmn4,jmn5,jmn6,nmax2,it,nmf1
             !DIR$ ATTRIBUTES ALIGN : 4 ::  jn,jp,ids,ms,mw,isf,iwv,iwf,nvs,itmin
             !DIR$ ATTRIBUTES ALIGN : 4 ::  itmax,ia,iang,ik,jk,itau,itau0,jtau0
             !DIR$ ATTRIBUTES ALIGN : 4 ::  jtau,nang,nang2,nsmall,j,NDGS,in,jj,nm1,nj
#endif              
             integer(kind=int4) :: i,j1,j2,n,in0,iv0,m,imn,iuv,imn1,        &
                           iuv1,nmax0,imax,n0,ii,ij,iv,inn,n1,      &     !! in original code was: in, here it was changed to: inn
                           ip,iq,is,isn,isv,nlarge,nbes,irc,        &
                           itrc,iuvc,niter,ijmax,ijmin,iuvp,        &
                           jv1,iuv2,iuv3,juv1,juv2,juv3,juv4,       &
                           juv5,juv6,ntemp,n2,imn2,imn3,jmn1,jmn2,  &
                           jmn3,jmn4,jmn5,jmn6,nmax2,it,nmf1,    &
                           jn,jp,ids,ms,mw,isf,iwv,iwf,nvs,itmin,   &
                           itmax,ia,iang,ik,jk,itau,itau0,jtau0,    &
                           jtau,nang,nang2,nsmall,j,NDGS,in,jj,nm1, &
                           nj
             
        
          
         
#if defined __INTEL_COMPILER              
              !DIR$ ATTRIBUTES ALIGN : 16 ::   A,B,cmz,Aj,Bj,A2,B2,Aj2,Bj2,A0,B0,ephi,ci,cin,ci0
             !DIR$ ATTRIBUTES ALIGN : 16 ::   A1,B1,Aj1,Bj1,Aj0,Bj0,cmzj,cmzg,Ag,Bg
#endif
               complex(16)  :: A,B,cmz,Aj,Bj,A2,B2,Aj2,Bj2,A0,B0,ephi,ci,cin,ci0,    &
                           A1,B1,Aj1,Bj1,Aj0,Bj0,cmzj,cmzg,Ag,Bg
               
         
          
          
           common/MIESUB/ twopi,pih
           common/rot/bcof,dc
           common/fnr/fnr
           common/pitau/pi,tau
           common/tran/atr
           common/ig0/iga0
           common/g0/ga0
           common/cofmnv0/cof0
           common/crot/cofsr
           common /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
#if defined __INTEL_COMPILER
           !DIR$ ATTRIBUTES ALIGN : 64 :: /MIESUB/
           !DIR$ ATTRIBUTES ALIGN : 64 :: /rot/
           
           !DIR$ ATTRIBUTES ALIGN : 64 :: /pitau/
           !DIR$ ATTRIBUTES ALIGN : 64 :: / tran/
           !DIR$ ATTRIBUTES ALIGN : 64 :: /ig0/
           !DIR$ ATTRIBUTES ALIGN : 64 :: /g0/
           !DIR$ ATTRIBUTES ALIGN : 64 :: /cofmnv0/
           !DIR$ ATTRIBUTES ALIGN : 64 :: / crot/
           !DIR$ ATTRIBUTES ALIGN : 64 :: /TMAT/
#endif          
           !Exec code ....
           aerr = -9999
           emsg = " "
           include 'tmatrix_mps_vtune_analysis.inc'
           ! Begin allocation
           allocate(drot(nrc,nij),STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(wmf1(np,np,0:np,0:np2),    &
                    wm1(np,np,0:np,0:np2),     &
                    wsdt(np,np,0:np,0:np2),     &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(atr0(ni0,nij),             &
                    btr0(ni0,nij),             &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(atr1(ni0,nij),             &
                    btr1(ni0,nij),             &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(ek(np,nij),                &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(p0(nLp,nmp),               &
                    q0(nLp,nmp),               &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(at0(nLp,nmp),              &
                    bt0(nLp,nmp),              &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(as(nLp,nmp),               &
                    bs(nLp,nmp),               &
                    as2(nLp,nmp),              &
                    bs2(nLp,nmp),              &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(tta(nLp,nLp,nmp),          &
                    ttb(nLp,nLp,nmp),          &
                    tta0(nLp,nLp,nmp),         &
                    ttb0(nLp,nLp,nmp),         &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(asr(nLp,nLp,nmp),          &
                    bsr(nLp,nLp,nmp),          &
                    as0(nLp,nLp,nmp),          &
                    bs0(nLp,nLp,nmp),          &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(asc(nLp,nLp,nmp),          &
                    bsc(nLp,nLp,nmp),          &
                    as1(nLp,nLp,nmp),          &
                    bs1(nLp,nLp,nmp),          &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(ast(nLp,nLp,nmp),          &
                    bst(nLp,nLp,nmp),          &
                    asp(nLp,nLp,nmp),          &
                    bsp(nLp,nLp,nmp),          &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(asv(nLp,nLp,nmp),          &
                    bsv(nLp,nLp,nmp),          &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(pct(nLp,nLp,2,2,nmp,nmp),  &
                    A1m(0:np,np,np,2,2),       &
                    A2m(0:np,np,np,2,2),       &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(B11n(np,np,2,2),           &
                    B12n(np,np,2,2),           &
                    B13n(np,np,2,2),           &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(B21n(np,np,2,2),           &
                    B22n(np,np,2,2),           &
                    B23n(np,np,2,2),           &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(fhmf1(0:np2,np,0:np2,2),   &
                    fmf1(0:np2,np,0:np2,2),    &
                    fhm1v(0:np2,np,0:np2,2),   &
                    fm1v(0:np2,np,0:np2,2),    &
                    fhm1q(0:np2,np,0:np2,2),   &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(fm1q(0:np2,np,0:np2,2),    &
                    fhmf1vq(0:np2,np,0:np2,2), &
                    fmf1vq(0:np2,np,0:np2,2),  &
                    fhmf4vq(0:np2,np,0:np2,2), &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999
           allocate(tbar(nLp,2,2,nmp,nmp),     &
                    tbar0(2,2,nmp,nmp),        &
                    STAT=aerr,ERRMSG=emsg)
           if(aerr /= 0) goto 9999     
           ! Initialization   of arrays
#if defined __INTEL_COMPILER
           !DIR$ VECTOR ALIGNED
           !DIR$ SIMD
#endif
           do i = 1, nLp
               nmax(i)    = 0
               uvmax(i)   = 0
              ! ind(i)     = 0
               x(i)       = ZERO
               xc(i)      = ZERO
               c0i(i)     = ZERO
               c1i(i)     = ZERO
           end do
           
          
           besj = ZERO
           besy = ZERO
          
          
          
           bes0 = ZERO
           confg = ZERO
           taup = ZERO
           taupj = ZERO
           taupg = ZERO
           taupjg = ZERO
           tau0p = ZERO
           tau1p = ZERO
           tau2p = ZERO
           tau0pg = ZERO
           tau1pg = ZERO
           tau2pg = ZERO
           tau0pj = ZERO
           tau1pj = ZERO
           tau2pj = ZERO
           tau0pjg = ZERO
           tau1pjg = ZERO
           tau2pjg = ZERO
           tau20 = ZERO
           tau11 = ZERO
           tau02 = ZERO
           tau20g = ZERO
           tau11g = ZERO
           tau02g = ZERO
           taum = ZERO
           taumg = ZERO
           w01s = ZERO
           wcf = ZERO
         
           wsdt = ZERO
           bcof = ZERO
           dc = ZERO
           pi = ZERO
           tau = ZERO
           cof0 = ZERO
           cofsr = ZERO
         
        
           atr0 = CZERO
           btr0 = CZERO
           atr  = CZERO
           at   = CZERO
           bt   = CZERO
           atr1 = CZERO
           btr1 = CZERO
           
           ref  = CZERO
           refc = CZERO
           p0   = CZERO
           q0   = CZERO
           an   = CZERO
           an   = CZERO
           bn   = CZERO
           B2i  = CZERO
           at0  = CZERO
           bt0  = CZERO
           at1  = CZERO
           bt1  = CZERO
           atj  = CZERO
           btj  = CZERO
          
          
          
          
          
           
          
         
           
          
           
           
          
          
           
          
           
           pih   = dacos(0.0_dp)
           twopi  = 4.0_dp*pih
           pione  = 2.0_dp*pih
           ci=dcmplx(0._dp,1._dp)
           cin=dcmplx(0._dp,-1._dp)
           ci0=dcmplx(0._dp,0._dp)
           gcs=0._dp
           gcv=0._dp
           eps=1.0E-20_dp
           fint=0._dp
           result = .false.
           ret    = -1
     ! OPEN(UNIT=1,FILE='gmm03Tr.in',STATUS='OLD')
     ! READ(1,'(a20)') FLNAME
!C-----------------------------------------------------------------------
!C  FLNAME: input file name for the configuration of an ensemble
!C          the first line in this file is the incident wavelength, 
!C          the second line is the total number of particles, 
!C          rest lines specify shape, orientation, position of
!C             particle-center, size, and complex refractive index for
!C             each component particles; one line for a particle and 
!C             each line must contain thirteen numbers, including 
!C                idshp, shp(1), shp(2), shp(3), 
!C                x, y, z, r, Re(m), Im(m), rc, Rec(m), Imc(m)
!C             idshp -- an integer number to specify particle shape 
!C                (the same as NP in Mishchenko's code)   
!C                 0 for a sphere 
!C                -1 for a spheroid 
!C                -2 for a finite cylinder
!C                -3 for generalized Chebyshev particles (describing 
!C                   the shape of distorted water drops)
!C                 idshp>0 for Chebyshev particles, which is the degree
!C                    of the Chebyshev polynomial
!C             shp(1) -- to specify the aspect ratio of a nonspherical
!C                particle (the parameter EPS in Mishchenko's 
!C               "ampld.new.f" code), 
!C                for a spheroid, it is the ratio of the horizontal to 
!C                    rotational axises 
!C                for a finite circular cylinder, it is the 
!C                    diameter-to-length ratio
!C                for Chebyshev particles, it is the deformation 
!C                    parameter
!C                for spheres, it is ineffective and can be set to 0
!C             shp(2) and shp(3) -- the azimuth and zenith angles of the 
!C                axis of rotational symmetry of a nonspherical particle
!C                to specify its orientation, i.e., the Euler angles of 
!C                alpha and beta to rotate the coordinate system so that 
!C                the z axis after rotation is along the axis of 
!C                rotational symmetry of the particle
!C                for spheres, simply put shp(2)=shp(3)=0 
!C             x,y,z -- the Cartesian coordinates of a particle-center
!C             r -- equivalent-sphere radius to specify particle size, 
!C                for spherical particles, it is the real sphere radius 
!C             Re(m), Im(m) -- the real and imaginary parts of the 
!C               refractive index of a particle 
!C                (for a core-mantle sphere, these refer to the mantle)
!C             rc, Rec(m), Imc(m) -- effective for a core-mantle sphere 
!C                only, which are the radius, real and imaginary parts
!C                of the refractive index of the core of a core-mantle  
!C                sphere; for all other types of particles these three 
!C                parameters are ineffective and can simply be set to 0
!C    (see the sample file given above)
!C
       !   write(6,'(a12,a20)') 'input file: ',FLNAME    
       !   READ(1,*) idMie
!C-----------------------------------------------------------------------
!C  default: idMie=0
!C  idMie=1: calculating only coherent scattering, no interaction 
!C-----------------------------------------------------------------------
        write(6,'(a7,i3)') 'idMie: ',idMie
     ! READ(1,*) small,MXINT
!C----------------------------------------------------------------------- 
!C  small: error tolerance for the iterative solution process for solving 
!C         the interacting scattering coefficients (1.d-6)
!C  MXINT: the maximum number of iterations allowed in the iterative 
!C  solution process
!C  This code uses BI-CGSTAB, the stabilized Bi-Conjugate Gradient method
!C  [see H.A. van der Vorst, SIAM J. Sci. Stat. Comput. 13, 631, (1992); 
!C  M.H. Gutknecht, SIAM J. Sci. comput. 14, pp.1020-1033 (1993)].
!C-----------------------------------------------------------------------
      write(6,'(a22,e10.2)')    'Convergence criterion:',small
      write(6,'(a37,i5)')       'Maximum number of iteration allowed:',MXINT
      
     ! READ(1,*) NADD
!C-----------------------------------------------------------------------
!C  NADD is the number of terms to add to the scattering orders required  
!C  by Wiscombe's criterion, which is for spherical component particles 
!C  only and can be negative or positive in the range of [-9,99]   
!C  (Normally, set NADD=0)
!C-----------------------------------------------------------------------
      write(6,'(a41,i3)')   'Scat. orders added to Wiscombe criterion:', NADD 
    ! +         
!      READ(1,*) idscmt,sang
!C----------------------------------------------------------------------- 
!C  idscmt: switch for calculation of scattering matrix
!C          (idscmt<0 for not calculating scattering matrix)
!C  sang: the scattering angle interval for output  
!C  example: when sang=1, the results in output will be written for every 
!C  degree of scattering angle from 0 to 180 
!C-----------------------------------------------------------------------  
      if(sang.lt.0.001_dp) sang=90._dp
         nang=90._dp/sang+1._dp
         nang2=2*nang-1
      if(nang2.gt.nangmax) then
         write(6,*) '  sang too small'
         write(6,*) '  please increase sang in tmatrix_mps_driver'
         write(6,*) '  and try again, or'
         write(6,*) '  increase nangmax in the parameter line of the'
         write(6,*) '  main code, recompile, then try again'
#if defined __INTEL_COMPILER
         call TRACEBACKQQ(STRING="Angle: sang too small!!", USER_EXIT_CODE= -1)
#elif defined __GFORTRAN__
         call backtrace()
#endif
         !
         stop
      endif
     ! close(1)
    !  write(6,'(/)')

     ! OPEN(UNIT=2,FILE=FLNAME,STATUS='OLD')
     ! READ(2,*) w,irat
!      C-----------------------------------------------------------------------
!C  w -- incident wavelength
!C  irat -- the same as the parameter RAT in Mishchenko's code
!C     when irat=1, particle size is specified in terms of the 
!C     equal-volume-sphere radius, otherwise, particle size is specified 
!C     in terms of the equal-surface-area-sphere radius
!C-----------------------------------------------------------------------
!      READ(2,*) nL
!C-----------------------------------------------------------------------
!C  nL -- number of spheres in the aggregate
!C-----------------------------------------------------------------------
      if(nL.gt.nLp) then
         write(6,*) 'Parameter nLp too small, must be >', nL 
         write(6,*) 'Change nLp in gmm01f.par, recompile, then try again'
         stop  
      end if
      if(nL.eq.1) idMie=1
!      C-----------------------------------------------------------------------
!C  input the configuration and particle parameters for an snsemble
!C  as mentioned above, one line for a particle and each line contains 
!C  13 numbers:
!C  idshp(i) - to specify particle shape 
!C     0 for a sphere, -1 for a spheroid, -2 for a finite cylinder, ... 
!C  shp(1,i) - aspect ratio for a nonspherical particle
!C  shp(2,i) and shp(3,i) - the azimuth and zenith angles of the 
!C     axis of rotational symmetry of a nonspherical particle
!C  the rest of the line is the same as in the input data for the 
!C  codes gmm02f.f, gmm02s.f, or gmm02TrA.f, which includes 9 numbers:
!C     x-, y-, z-coordinates of a particle-center, the radius of  
!C     an (equivalent) sphere in the same unit of the incident 
!C     wavelength, the real and imaginary parts of the refractive 
!C     index, the last three numbers are for a core-mantle sphere, 
!C     which are the radius, the real and imaginary parts of 
!C     refractive index of the core, for all other shapes of particles, 
!C     simply set these three to 0. 
      !C-----------------------------------------------------------------------
         r00 = ZERO
         do 1 i=1,nL
        ! read(2,*,err=10) idshp(i),(shp(j,i),j=1,3),
    ! +                             (r0(j,i),j=1,9)
           if(idshp(i).eq.0) then
              do j=1,3
                  shp(j,i)=0._dp
              enddo
         endif
         temp=dabs(shp(1,i)-1._dp)
         temp0=0.0000001_dp
         if(idshp(i).eq.-1.and.temp.lt.temp0) then
            idshp(i)=0
            do j=1,3
               shp(j,i)=0._dp
            enddo
         endif
         x0=r0(1,i)
         y0=r0(2,i)
         z0=r0(3,i)
         r00(1,i)=x0
         r00(2,i)=y0
         r00(3,i)=z0
         r0(6,i)=dabs(r0(6,i))
         r0(9,i)=dabs(r0(9,i))
         ref(i)=dcmplx(r0(5,i),r0(6,i))
         refc(i)=dcmplx(r0(8,i),r0(9,i)) 
         gcs=gcs+r0(4,i)*r0(4,i)
         gcv=gcv+r0(4,i)*r0(4,i)*r0(4,i)
1        continue
         !close(2)
         gcsr=dsqrt(gcs)
         gcvr=gcv**0.333333333333333333333333333333_dp
         goto 11
 !10   write(6,*) 'fatal error in the input file'
 !     stop
 11      k=twopi/w
        xv=k*gcvr
         xs=k*gcsr
         if(irat.eq.1) then 
             write(6,'(a,f7.3)') ' volume-equiv. xv:  ', xv
         else 
             write(6,'(a,f7.3)') ' surface-equiv. xs: ', xs
         endif
       !  write(6,'(/)')
    !  fileout='gmm03Trd.out'
     ! fileout1=fileout
      if(idMie.eq.1) then
         write(6,*)  &
           '*** Calculating only coherent scattering ***'
         write(6,*) &
           '*** No interaction included ****************'
      endif
      do i=1,nL
         x(i)=k*r0(4,i)
         if(r0(7,i).lt.0.000000000001_dp) then
            xc(i)=0._dp
            refc(i)=dcmplx(0._dp,0._dp)
         else
            xc(i)=k*r0(7,i)
         endif
      enddo
!      C ----------------------------------------------------------------------
!C  Calculating individual T-matrix for each component particles 
!C  For spherical particles, it only needs to calculate the "Mie" 
!C  scattering coefficients
!C  (1) the subroutine "scoatabd.f" used here for calculating scattering
!C      coefficients of homogeneous and core-mantle spheres is originally 
!C      the code "SCSERD.FOR" written by R.T. Wang and W.X. Wang  
!C      [see R.T. Wang and W.X. Wang, "Scattering by Arbitrarily Large 
!C      Homogeneous/Concentric Speres - Exact Theory with Use of New 
!C      Efficient Algorithms," in Proc. of the 1985 CRDC Scientific 
!C      Conference on Obscuration and Aerosol Research, R. Kohl, ed. 
!C      (Army CRDEC-SP-86019, Aberdeen, MD 1986), pp. 381-409], which 
!C      uses ratio method of Wang and van de Hulst in the calculation of 
!C      Riccati-Bessel functions [see Wang and van de Hulst, Appl. Opt. 
!C      30, 106 (1991), Xu, Gustafson, Giovane, Blum, and Tehranian, 
!C      Phys. Rev. E 60, 2347 (1999)]
!C  (2) for nonspherical particles, the individual-particle T-matrices  
!C      in the particle reference system are calculated by the subroutine 
!C      "tmod.f", which is a part of Mishchenko's public domain code 
!C      "ampld.new.f" 
!C ----------------------------------------------------------------------
    
      nmax0=1
      tbar = CZERO

      ! Guard this with conditional compilation
#if defined __INTEL_COMPILER
#if (USE_PERF_PROFILER) == 1
      include 'tmatrix_mps_perf_loop744.inc'    
#endif
#elif defined __GFORTRAN__
#if (USE_PERF_PROFILER) == 1
      ! ADD gfortran version
#endif
#endif
      line     = __LINE__       
      start    = 0.0_dp
      end      = 0.0_dp
      duration = 0.0_dp
      start = DCLOCK()
#if defined __INTEL_COMPILER
      !DIR$ LOOP COUNT (10,100,300)
#endif
      do i=1,nL
         do j1=1,nmp
#if defined __INTEL_COMPILER
            !DIR$       UNROLL (4)
#elif defined __GFORTRAN__
            !GCC$ UNROLL 4
#endif
            do j2=1,nmp
               tbar(i,1,1,j1,j2)=0._dp
               tbar(i,1,2,j1,j2)=0._dp
               tbar(i,2,1,j1,j2)=0._dp
               tbar(i,2,2,j1,j2)=0._dp
            enddo
         enddo
         if(i.eq.1) goto  12

#if defined __INTEL_COMPILER           
!DIR$    DISTRIBUTE POINT
#endif            
         do 121 j=i-1,1,-1
            if(idshp(i).eq.idshp(j).and.shp(1,i).eq.shp(1,j)) then
               if(xc(i).eq.xc(j).and.refc(i).eq.refc(j)) then
                  if(x(i).eq.x(j).and.ref(i).eq.ref(j)) then
                     nmax(i)=nmax(j)
                     uvmax(i)=uvmax(j)
                     do n=1,nmax(j)
                        do v=1,nmax(j)
			   in0=n*(n+1)
			   iv0=v*(v+1)
			   tbar(i,1,1,in0,iv0)=tbar(j,1,1,in0,iv0)
                           tbar(i,1,2,in0,iv0)=tbar(j,1,2,in0,iv0)
                           tbar(i,2,1,in0,iv0)=tbar(j,2,1,in0,iv0)
                           tbar(i,2,2,in0,iv0)=tbar(j,2,2,in0,iv0)
			               do  m=1,min(n,v)
			                    imn=m+in0
			                    iuv=m+iv0
			                    A=tbar(j,1,1,imn,iuv)
			                    tbar(i,1,1,imn,iuv)=A
			                    A=tbar(j,1,2,imn,iuv)
			                    tbar(i,1,2,imn,iuv)=A
			                    A=tbar(j,2,1,imn,iuv)
			                    tbar(i,2,1,imn,iuv)=A
			                    A=tbar(j,2,2,imn,iuv)
			                    tbar(i,2,2,imn,iuv)=A
                                            imn=-m+in0
			                    iuv=-m+iv0
			                    A=tbar(j,1,1,imn,iuv)
			                    tbar(i,1,1,imn,iuv)=A
                                            A=tbar(j,1,2,imn,iuv)
			                    tbar(i,1,2,imn,iuv)=A
			                    A=tbar(j,2,1,imn,iuv)
			                    tbar(i,2,1,imn,iuv)=A
			                    A=tbar(j,2,2,imn,iuv)
			                    tbar(i,2,2,imn,iuv)=A
			               enddo
                        enddo
                     enddo
!c		     write(6,*) 
!c     +                  'same type of particles: ',i,' and ',j
                     goto 15
                  endif
               endif
            endif
 121     continue
 12      if(idshp(i).eq.0) then
           ratio=xc(i)/x(i)
            if(ratio.gt.1._dp) then
               write(6,*) 'size of core > mantle for particle ',i
#if defined __INTEL_COMPILER              
               call TRACEBACKQQ(STRING="size of core > mantle for particle", USER_EXIT_CODE= -1)
#elif defined __GFORTRAN__
               call backtrace()
#endif
	          stop
            endif
            if(ratio.lt.0.00000000000001_dp) ratio=0._dp	   
            call scoatabd(x(i),ratio,r0(8,i),-r0(9,i),r0(5,i),     &
                         -r0(6,i),np,an,bn,NADD,nmax(i))
            if(nmax(i).gt.np) then
               write(6,*) '  Parameter np too small, must be >',  &
                             nmax(i)
               write(6,*) '  Please change np in tmatrix_mps_driver,' 
               write(6,*) '  recompile, then try again'
#if defined __INTEL_COMPILER
               call TRACEBACKQQ(STRING="Parameter np too small",USER_EXIT_CODE= -1)
#elif defined __GFORTRAN__
               call backtrace()
#endif
               stop
            endif
            uvmax(i)=nmax(i)*(nmax(i)+2)
            write(6,'(a,1x,i4)')     &
              'Actual single-particle expansion truncation:',   &
               nmax(i)
            do j=1,uvmax(i)
               v=dsqrt(dble(j))
               tbar(i,1,1,j,j)=an(v)
               tbar(i,2,2,j,j)=bn(v)   
            enddo
	        write(6,'(a,i3,a,f10.4)') 'sphere #',i,    &
           '  individual size parameter: ',x(i)
	        goto 15
        endif
         RAT=irat
         DDELT=0.001_dp
         NDGS=2
#if defined __INTEL_COMPILER
         !DIR$   LOOP COUNT ( 101)
#endif
	 do m=1,NPN6
#if defined __INTEL_COMPILER
     !DIR$   LOOP COUNT (100)
#endif
     do n=1,NPN4
#if defined __INTEL_COMPILER
        !DIR$ IVDEP:LOOP
#endif
#if defined __INTEL_COMPILER
        !DIR$ UNROLL (4)
#elif defined __GFORTRAN__
        !GCC$ UNROLL 4
#endif
	       do v=1,NPN4
	              RT11(m,n,v)=0._dp
		      RT12(m,n,v)=0._dp
		      RT21(m,n,v)=0._dp
		      RT22(m,n,v)=0._dp
		      IT11(m,n,v)=0._dp
		      IT12(m,n,v)=0._dp
		      IT21(m,n,v)=0._dp
		      IT22(m,n,v)=0._dp
	       enddo
	    enddo
	 enddo
         call tm0d(w,idshp(i),shp(1,i),r0(4,i),RAT,  &
                  r0(5,i),r0(6,i),DDELT,NDGS,nmax(i))
         write(6,'(a,i3,a,f10.4)') 'particle #',i,  &
           '   individual size parameter: ',x(i)
	 write(6,*) 'individual T-matrix: ', i, '   ', nmax(i) 
         if(nmax(i).gt.np) then
            write(6,*) ' Parameter np too small, must be >',     &
                      nmax(i)
            write(6,*) ' Please change np in gmm01f.par,' 
            write(6,*) '   recompile,then try again'
#if defined __INTEL_COMPILER
            call TRACEBACKQQ(STRING="Parameter np too small",USER_EXIT_CODE = -1)
#elif defined __GFORTRAN__
            call backtrace()
#endif
            stop
         endif
         if(nmax(i).gt.NPN1) then
            write(6,*) ' Parameter NPN1 too small, must be >',  &
                        nmax(i)
            write(6,*) '   Please change NPN1 in ampld.par.f,' 
            write(6,*) '   recompile,then try again'
#if defined __INTEL_COMPILER
            call TRACEBACKQQ(STRING="Parameter NPN1 too small",USER_EXIT_CODE = -1)
#elif defined __GFORTRAN__
            call backtrace()
#endif
            stop
         endif
         uvmax(i)=nmax(i)*(nmax(i)+2)

         do n=1,nmax(i)
            in0=n*(n+1)
            do v=1,nmax(i)
               iv0=v*(v+1)
               nsmall=min(n,v)
	       A=-ci**(v-n)
	       imn=in0
               iuv=iv0
               B=dcmplx(RT11(1,n,v),IT11(1,n,v))
               tbar(i,1,1,imn,iuv)=A*B
               B=dcmplx(RT12(1,n,v),IT12(1,n,v))
               tbar(i,1,2,imn,iuv)=A*B
               B=dcmplx(RT21(1,n,v),IT21(1,n,v))
               tbar(i,2,1,imn,iuv)=A*B
               B=dcmplx(RT22(1,n,v),IT22(1,n,v))
               tbar(i,2,2,imn,iuv)=A*B
               do m=1,nsmall
                  imn=m+in0
                  iuv=m+iv0
                  imn1=-m+in0
                  iuv1=-m+iv0
                  B=dcmplx(RT11(m+1,n,v),IT11(m+1,n,v))
                  tbar(i,1,1,imn,iuv)=A*B
                  tbar(i,1,1,imn1,iuv1)=tbar(i,1,1,imn,iuv)
                  B=dcmplx(RT12(m+1,n,v),IT12(m+1,n,v))
                  tbar(i,1,2,imn,iuv)=A*B
                  tbar(i,1,2,imn1,iuv1)=-tbar(i,1,2,imn,iuv)
                  B=dcmplx(RT21(m+1,n,v),IT21(m+1,n,v))
                  tbar(i,2,1,imn,iuv)=A*B
                  tbar(i,2,1,imn1,iuv1)=-tbar(i,2,1,imn,iuv)
                  B=dcmplx(RT22(m+1,n,v),IT22(m+1,n,v))
                  tbar(i,2,2,imn,iuv)=A*B
                  tbar(i,2,2,imn1,iuv1)=tbar(i,2,2,imn,iuv)
               enddo
            enddo
         enddo  
15       if(nmax(i).gt.nmax0) then
            nmax0=nmax(i)
	    imax=i
	 endif
      enddo
      end = DCLOCK()
      duration = end - start
      print*, "Loop at line: ", line+7, " executed in: ", duration, " microseconds."
      write(6,*) 'maximum scattering order: ',imax,'   ',nmax0
      write(6,'(/)')
      write(6,*) 'input particle-positions: '
      i=1
      write(6,'(i5,3f14.5)') i,r0(1,1),r0(2,1),r0(3,i)
      i=nL
      write(6,'(i5,3f14.5)') i,r0(1,i),r0(2,i),r0(3,i)
!C-----------------------------------------------------------------------
!C  calculating constants and Gaunt coefficients
!C-----------------------------------------------------------------------
      n0=nmax0+2
      fnr(0)=0._dp
      do n=1,4*(nmax0+1)
         fnr(n)=dsqrt(dble(n))
      enddo
      bcof(0)=1._dp
      do n=0,n0-1
         bcof(n+1)=fnr(n+n+2)*fnr(n+n+1)*bcof(n)/fnr(n+1)/fnr(n+1)
      enddo
!C-----------------------------------------------------------------------
!C  calculating T-matrices of individual particles in their respective
!C  specified orientations
!C-----------------------------------------------------------------------
          tbar0 = CZERO
          ekt = CZERO
          

          ! wrap it conditional compilation block
#if defined __INTEL_COMPILER
#if (USE_PERF_PROFILER) == 1
        include 'tmatrix_mps_perf_loop965.inc' 
#endif
#elif defined __GFORTRAN__
#if (USE_PERF_PROFILER) == 1

#endif
#endif
        line  = __LINE__     
        start = 0.0_dp
        end   = 0.0_dp
        duration = 0.0_dp
        start = DCLOCK()
#if defined __INTEL_COMPILER
        !DIR$   LOOP COUNT (10,100,300)
#endif
        do 17 i=1,nL
             if(idshp(i).eq.0) goto 17
                alph=shp(2,i)*pih*0.01111111111111111111111111111111_dp
                beta=shp(3,i)*pih*0.01111111111111111111111111111111_dp
                s=dabs(alph)
                t=dabs(beta)
                if(s.lt.0.0000000001_dp.and.t.lt.0.00000000001_dp) goto 17
                if(i.eq.1) goto  172
          do 171 j=i-1,1,-1
                ii=idshp(i)
                ij=idshp(j)
	        ca=shp(2,j)
	        s=shp(1,i)
	        t=shp(1,j)
            if(ii.eq.ij.and.s.eq.t) then
               if(x(i).eq.x(j).and.ref(i).eq.ref(j)) then
                  sa=shp(2,i)
                  ca=shp(2,j)
                  sb=shp(3,i)
                  cb=shp(3,j)
                  if(sa.eq.ca.and.sb.eq.cb) then
                     do inn=1,uvmax(j)
                        do iv=1,uvmax(j)
			   tbar(i,1,1,inn,iv)=tbar(j,1,1,inn,iv)
                           tbar(i,1,2,inn,iv)=tbar(j,1,2,inn,iv)
                           tbar(i,2,1,inn,iv)=tbar(j,2,1,inn,iv)
                           tbar(i,2,2,inn,iv)=tbar(j,2,2,inn,iv)
                        enddo
                     enddo
                     goto 17
                  endif
               endif
            endif
 171     continue
 172     do in=1,uvmax(i)
            do iv=1,uvmax(i)
               tbar0(1,1,inn,iv)=tbar(i,1,1,inn,iv)
               tbar0(1,2,inn,iv)=tbar(i,1,2,inn,iv)
               tbar0(2,1,inn,iv)=tbar(i,2,1,inn,iv)
               tbar0(2,2,inn,iv)=tbar(i,2,2,inn,iv)
            enddo
         enddo
         sa=dsin(alph)
         ca=dcos(alph)
         cb=dcos(beta)
         A=dcmplx(ca,sa)
         n1=nmax(i)
	 ekt(0)=1._dp
         do m=1,2*n1
            ekt(m)=A**m
	    ekt(-m)=qconjg(ekt(m))
         enddo
         call rotcoef(cb,n1)
         do m=-n1,n1
            do u=-n1,n1
               A=ekt(u-m)
               do 173 n=1,n1
                  if(iabs(m).gt.n) go to 173
                    in0=n*(n+1)
	            imn=in0+m
                  do 174 v=1,n1
                     if(iabs(u).gt.v) go to 173
                     iv0=v*(v+1)
                     iuv=iv0+u
                     nsmall=min(n,v)
                     do ip=1,2
                        do iq=1,2
                           bar(ip,iq)=0._dp
                        enddo
                     enddo
                     do is=-nsmall,nsmall
                        isn=in0+is
                        isv=iv0+is
                        t=dc(m,isn)*dc(u,isv)
                        B=A*t
                        do ip=1,2
                           do iq=1,2
                              Bj=tbar0(ip,iq,isn,isv)
                              bar(ip,iq)=bar(ip,iq)+B*Bj
                           enddo
                        enddo	 
                     enddo      
                     do ip=1,2
                        do iq=1,2
                           tbar(i,ip,iq,imn,iuv)=bar(ip,iq)
                        enddo
                     enddo
 174              continue
 173           continue
            enddo
         enddo
17     continue
       end = DCLOCK()
       duration = end - start
       print*, "Loop at line: ", line+4, "executed in: ", duration, " microseconds."
!        C
!C  calculating Gaunt coefficients
!C  the formulation used here for the calculation of Gaunt coefficients 
!C  can be found in Bruning and Lo, IEEE Trans. Anten. Prop. Ap-19, 378 
!C  (1971) and Xu, J. Comput. Appl. Math. 85, 53 (1997), J. Comput. Phys. 
!C  139, 137 (1998)
!C
      call cofsrd(nmax0)	
      call cofd0(nmax0)
      call cofnv0(nmax0)
      call gau0(nmax0)
!      C-----------------------------------------------------------------------
!C  calculating rotational and translation coefficients
!C-----------------------------------------------------------------------
      ek = CZERO
      drot = ZERO

     
      include 'tmatrix_mps_perf_loop1097.inc' 
    

      line  = __LINE__       
      start = 0.0_dp
      end   = 0.0_dp
      duration = 0.0_dp
      start = DCLOCK()
!DIR$ LOOP COUNT (9,99,299)
      do i=1,nL-1
!DIR$ LOOP COUNT (10,100,300)
         do j=i+1,nL
            ij=(j-1)*(j-2)/2+j-i
            x0=r0(1,i)-r0(1,j)
            y0=r0(2,i)-r0(2,j)
            z0=r0(3,i)-r0(3,j)
            call carsphd(x0,y0,z0,d,xt,sphi,cphi)
            temp=(r0(4,i)+r0(4,j))/d
            confg(1,ij)=x0
            confg(2,ij)=y0
            confg(3,ij)=z0
            confg(4,ij)=d
            confg(5,ij)=temp
            ephi=dcmplx(cphi,sphi)
            nlarge=max(nmax(i),nmax(j))
            do m=1,nlarge
               ek(m,ij)=ephi**m
            enddo
            xd=k*d
            nbes=2*nlarge+1
            call besseljd(nbes,xd,besj)
            call besselyd(nbes,xd,besy)
!C
!C  calculating reduced rotation matrix elements (the subroutine 
!C  rotcoef used here is originally written by Mackowski)
!C
            call rotcoef(xt,nlarge)
            irc=0
            do n=1,nlarge
               n1=n*(n+1)
               do u=-n,n
                  do m=-n,n
                     imn=n1+m
                     irc=irc+1
                     drot(irc,ij)=dc(u,imn)
                  enddo
               enddo
            enddo
            itrc=0
            nsmall=min(nmax(i),nmax(j))
!C
!C  the formulation used here for the calculation of vector translation 
!C  coefficients are from Cruzan, Q. Appl. Math. 20, 33 (1962) and  
!C  Xu, J. Comput. Phys. 139, 137 (1998)
!C
            do m=-nsmall,nsmall
               n1=max(1,iabs(m))
               do n=n1,nlarge
                  do v=n1,nlarge
                     itrc=itrc+1
                     call cofxuds0(nmax0,m,n,v,besj,besy,         &
                                 atr0(itrc,ij),btr0(itrc,ij),     &
                                 atr1(itrc,ij),btr1(itrc,ij))     
                  enddo
               enddo
            enddo
         enddo
      enddo
      end = DCLOCK()
      duration = end  - start
      print*, "Loop at line: ", line+4, " executed in: ", duration, " microseconds."
      if(idMie.eq.1) then
         do j=1,nL
            do iuv=1,uvmax(j)
               do i=1,nL
                  do imn=1,uvmax(i)
                     pct(j,i,1,1,imn,iuv)=0.0_dp
                     pct(j,i,1,2,imn,iuv)=0.0_dp
                     pct(j,i,2,1,imn,iuv)=0.0_dp
                     pct(j,i,2,2,imn,iuv)=0.0_dp
                     if(j.eq.i) then
                        pct(j,i,1,1,imn,iuv)=tbar(i,1,1,imn,iuv)
                        pct(j,i,1,2,imn,iuv)=tbar(i,1,2,imn,iuv)
                        pct(j,i,2,1,imn,iuv)=tbar(i,2,1,imn,iuv)
                        pct(j,i,2,2,imn,iuv)=tbar(i,2,2,imn,iuv)
                     endif     
                  enddo
               enddo
            enddo
         enddo
         goto 1800
      endif
!C-----------------------------------------------------------------------
!c  begins iteration process to solve T-matrix columns of (u,v,q)
!c  BI-CGSTAB [see van der Vorst, SIAM J. Sci. Stat. Comput. 13, 631 
!c  (1992); Gutknecht, SIAM J. Sci. comput. 14, pp. 1020-1033 (1993)]   
!C-----------------------------------------------------------------------
      do i=1,nL-1
         do j=i+1,nL
            ij=(j-1)*(j-2)/2+j-i
            d=confg(4,ij)
            xd=k*d
            bes0(ij)=dsin(xd)/xd	      
         enddo
      enddo      
      write(6,*) 'Starting Bi-CGSTAB to solve T-matrix'         
      n0=nmax0*(nmax0+2)
      asr = CZERO
      bsr = CZERO
      as  = CZERO
      bs  = CZERO
      tta0 = CZERO
      ttb0 = CZERO
      tta  = CZERO
      ttb  = CZERO
      as1  = CZERO
      bs1  = CZERO
      ast  = CZERO
      bst  = CZERO
      as0  = CZERO
      bs0  = CZERO
      asp  = CZERO
      bsp  = CZERO
      asv  = CZERO
      bsv  = CZERO

  
      include 'tmatrix_mps_perf_loop1179.inc'    
      

      line  = __LINE__      
      start = 0.0_dp
      end = 0.0_dp
      duration = 0.0_dp
      start = DCLOCK()
      do 1001 iuv=1,n0
         v=dsqrt(dble(iuv))
         iuvc=v*v
!c         iuvc=1
         u=iuv-v*(v+1)
         do 1002 iq=1,2
            iuv1=-u+v+v*v
            do imn=1,iuvc-1
               n=dsqrt(dble(imn))
               m=imn-n*n-n
               imn1=-m+n*n+n
               cz=(-1)**(m+u+n+v)
               do j=1,nL
                  do i=1,nL
                     A=pct(i,j,iq,1,iuv1,imn1)
                     asr(j,i,imn)=cz*A
                     A=pct(i,j,iq,2,iuv1,imn1)
                     bsr(j,i,imn)=cz*A
                  enddo
               enddo
            enddo
            do j=1,nL
	        ind(j)=0
	       if(iuvc.eq.1) goto 10021
               do i=1,nL
                  do imn=1,iuvc-1
                     as(i,imn)=asr(j,i,imn)
                     bs(i,imn)=bsr(j,i,imn)
                  enddo
                  do imn=iuvc,uvmax(i)
                     as(i,imn)=ci0
                     bs(i,imn)=ci0
                  enddo
               enddo
               call transT(nL,r0,nmax,uvmax,fint,atr0,btr0,     &
                   ek,drot,as,bs,as2,bs2,ind,confg,iuvc,2)
               do i=1,nL
                  do imn=iuvc,uvmax(i)
                     tta0(j,i,imn)=as2(i,imn)
                     ttb0(j,i,imn)=bs2(i,imn)
                  enddo
               enddo
	       goto 10022
10021	       do i=1,nL
                  do imn=iuvc,uvmax(i)
                     tta0(j,i,imn)=ci0
                     ttb0(j,i,imn)=ci0
                  enddo
               enddo
10022          continue
            enddo	      	      	      	      
            do j=1,nL
               do i=1,nL
                  do imn=iuvc,n0
                     asr(j,i,imn)=ci0
                     bsr(j,i,imn)=ci0
                  enddo
               enddo
            enddo
            do j=1,nL
               if(iuv.le.uvmax(j)) then
	             do imn=iuvc,uvmax(j)
                        asr(j,j,imn)=tbar(j,1,iq,imn,iuv)
		        bsr(j,j,imn)=tbar(j,2,iq,imn,iuv)
		     enddo
               endif
            enddo
            do j=1,nL
               temp=0._dp
               do imn=iuvc,uvmax(j)
                  A=asr(j,j,imn)
                  temp0=A*qconjg(A)
                  temp=temp+temp0
                  A=bsr(j,j,imn)
                  temp0=A*qconjg(A)
                  temp=temp+temp0
               enddo
               c0i(j)=temp
            enddo
            niter=1
            do j=1,nL
               do i=1,nL
                  do imn=1,iuvc-1
                     as(i,imn)=ci0
                     bs(i,imn)=ci0
                  enddo
                  do imn=iuvc,uvmax(i)
                     as(i,imn)=asr(j,i,imn)
                     bs(i,imn)=bsr(j,i,imn)
                  enddo
               enddo
               call transT(nL,r0,nmax,uvmax,fint,atr0,btr0,  &
                   ek,drot,as,bs,as2,bs2,ind,confg,iuvc,1)
               do i=1,nL
                  do imn=iuvc,uvmax(i)
                     tta(j,i,imn)=as2(i,imn)+tta0(j,i,imn)
                     ttb(j,i,imn)=bs2(i,imn)+ttb0(j,i,imn)
                  enddo
               enddo
            enddo
            do 611 i=1,nL
	           c1i(i)=0._dp
	        if(iuv.gt.uvmax(i)) goto 611
               do 6111 imn=iuvc,uvmax(i)
                  n=dsqrt(dble(imn))
                  do 6112 j=1,nL
                     if(iuv.gt.uvmax(j)) goto 6112
                     as1(j,i,imn)=0._dp
                     bs1(j,i,imn)=0._dp
                     ijmax=max(i,j)
                     ijmin=min(i,j)
                     if(ijmax.eq.ijmin) then
                        cz=1._dp
                     else
                        ij=(ijmax-1)*(ijmax-2)/2
			ij=ij+ijmax-ijmin
                        cz=bes0(ij)
                     endif
                     A=tbar(i,1,iq,imn,iuv)
                     as1(j,i,imn)=as1(j,i,imn)+cz*A
                     A=tbar(i,2,iq,imn,iuv)
                     bs1(j,i,imn)=bs1(j,i,imn)+cz*A
                     do 6113 jj=1,nL
		        if(iuv.gt.uvmax(jj)) goto 6113
                        ijmax=max(jj,j)
                        ijmin=min(jj,j)
                        if(ijmax.eq.ijmin) then
                           cz=1._dp
                        else
                           ij=(ijmax-1)*(ijmax-2)/2
			   ij=ij+ijmax-ijmin
                           cz=bes0(ij)
                        endif
                        A0=asr(jj,i,imn)
                        B0=bsr(jj,i,imn)
			do iuvp=iuvc,uvmax(i)
			   A=tbar(i,1,1,imn,iuvp)
                           A0=A0+A*tta(jj,i,iuvp)
                           B=tbar(i,1,2,imn,iuvp)
                           A0=A0+B*ttb(jj,i,iuvp)
                           A=tbar(i,2,1,imn,iuvp)
                           B0=B0+A*tta(jj,i,iuvp)
                           B=tbar(i,2,2,imn,iuvp)
                           B0=B0+B*ttb(jj,i,iuvp)
                        enddo
                        as1(j,i,imn)=as1(j,i,imn)-cz*A0
                        bs1(j,i,imn)=bs1(j,i,imn)-cz*B0
 6113                continue
                     A=as1(j,i,imn)
                     B=bs1(j,i,imn)
                     c1i(i)=c1i(i)+A*qconjg(A)
                     c1i(i)=c1i(i)+B*qconjg(B)
 6112	          continue
 6111 	       continue
 611        continue
            temp=0._dp
            B0=0._dp
            do 612 i=1,nL
               cext0=c1i(i)/c0i(i)
               if(cext0.lt.small) ind(i)=1
               if(ind(i).gt.0) goto 612
               if(cext0.gt.temp) temp=cext0
               B0=B0+c1i(i)
 612	    continue
            if(temp.lt.small) then
               do i=1,nL
                  do imn=1,uvmax(i)
                     do j=1,nL
                        pct(j,i,1,iq,imn,iuv)=asr(j,i,imn)
                        pct(j,i,2,iq,imn,iuv)=bsr(j,i,imn)
                     enddo
                  enddo
               enddo
               goto 1002
            endif
            do 613 i=1,nL
               if(ind(i).gt.0) goto 613
               do imn=1,iuvc-1
                  do j=1,nL
                     asp(j,i,imn)=ci0
                     bsp(j,i,imn)=ci0
                     as0(j,i,imn)=ci0
                     bs0(j,i,imn)=ci0
                  enddo
               enddo
               do imn=iuvc,uvmax(i)
                  do j=1,nL
                     asp(j,i,imn)=as1(j,i,imn)
                     bsp(j,i,imn)=bs1(j,i,imn)
                     as0(j,i,imn)=as1(j,i,imn)
                     bs0(j,i,imn)=bs1(j,i,imn)
                  enddo
               enddo
 613        continue
            do j=1,nL
               do i=1,nL
                  do imn=1,iuvc-1
                     as(i,imn)=ci0
                     bs(i,imn)=ci0
                  enddo
                  do imn=iuvc,uvmax(i)
                     as(i,imn)=asp(j,i,imn)
                     bs(i,imn)=bsp(j,i,imn)
                  enddo
               enddo
               call transT(nL,r0,nmax,uvmax,fint,atr0,btr0,   &
                   ek,drot,as,bs,as2,bs2,ind,confg,iuvc,1)
               do i=1,nL
                  do imn=iuvc,uvmax(i)
                     tta(j,i,imn)=as2(i,imn)
                     ttb(j,i,imn)=bs2(i,imn)
                  enddo
               enddo
            enddo
            A0=0.d0
            do 614 i=1,nL
               if(ind(i).gt.0) goto 614
               do 6141 imn=iuvc,uvmax(i)
                  n=dsqrt(dble(imn))
                  do 6142 j=1,nL
                     if(iuv.gt.uvmax(j)) goto 6142
                     ast(j,i,imn)=0._dp
                     bst(j,i,imn)=0._dp
                     do 6143 jj=1,nL
		            if(iuv.gt.uvmax(jj)) goto 6143
                        ijmax=max(jj,j)
                        ijmin=min(jj,j)
                        if(ijmax.eq.ijmin) then
                           cz=1._dp
                        else
                           ij=(ijmax-1)*(ijmax-2)/2+ijmax-ijmin
                           cz=bes0(ij)
                        endif
                        Aj2=asp(jj,i,imn)
                        Bj2=bsp(jj,i,imn)
			do iuvp=iuvc,uvmax(i)
			   A=tbar(i,1,1,imn,iuvp)
                           Aj2=Aj2+A*tta(jj,i,iuvp)
                           B=tbar(i,1,2,imn,iuvp)
                           Aj2=Aj2+B*ttb(jj,i,iuvp)
                           A=tbar(i,2,1,imn,iuvp)
                           Bj2=Bj2+A*tta(jj,i,iuvp)
                           B=tbar(i,2,2,imn,iuvp)
                           Bj2=Bj2+B*ttb(jj,i,iuvp)
                        enddo
                        ast(j,i,imn)=ast(j,i,imn)+cz*Aj2
                        bst(j,i,imn)=bst(j,i,imn)+cz*Bj2
 6143                continue	          
                     A0=A0+conjg(as0(j,i,imn))*ast(j,i,imn)
                     A0=A0+conjg(bs0(j,i,imn))*bst(j,i,imn)
 6142             continue
 6141          continue
 614        continue
            if(abs(A0).lt.1.0e-200_dp) then
               do i=1,nL
                  do imn=1,uvmax(i)
                     do j=1,nL
                        pct(j,i,1,iq,imn,iuv)=asr(j,i,imn)
                        pct(j,i,2,iq,imn,iuv)=bsr(j,i,imn)
                     enddo
                  enddo
               enddo  
               goto 1002
            endif
            Aj=B0/A0
  62        do 621 i=1,nL
               if(ind(i).gt.0) goto 621
               do imn=iuvc,uvmax(i)
                  do 6211 j=1,nL
                     if(iuv.gt.uvmax(j)) then
                        asv(j,i,imn)=0._dp
                        bsv(j,i,imn)=0._dp
                        go to 6211
                     endif
                     asv(j,i,imn)=asp(j,i,imn)-Aj*ast(j,i,imn)
                     bsv(j,i,imn)=bsp(j,i,imn)-Aj*bst(j,i,imn)
 6211             continue 
               enddo
 621        continue
            do j=1,nL
               do i=1,nL
                  do imn=1,iuvc-1
                     as(i,imn)=ci0
                     bs(i,imn)=ci0
                  enddo
                  do imn=iuvc,uvmax(i)
                     as(i,imn)=asv(j,i,imn)
                     bs(i,imn)=bsv(j,i,imn)
                  enddo
               enddo
               call transT(nL,r0,nmax,uvmax,fint,atr0,btr0,  &
                   ek,drot,as,bs,as2,bs2,ind,confg,iuvc,1)
               do i=1,nL
                  do imn=iuvc,uvmax(i)
                     tta(j,i,imn)=as2(i,imn)
                     ttb(j,i,imn)=bs2(i,imn)
                  enddo
               enddo
            enddo
            A2=0.0_dp
            B2=0.0_dp
            do 622 i=1,nL
               if(ind(i).gt.0) goto 622
               do 6221 imn=iuvc,uvmax(i)
                  n=dsqrt(dble(imn))
                  do 6222 j=1,nL
                     if(iuv.gt.uvmax(j)) goto 6222
                         asc(j,i,imn)=0.0_dp
                         bsc(j,i,imn)=0.0_dp
                        do 6223 jj=1,nL
		             if(iuv.gt.uvmax(jj)) goto 6223
                              ijmax=max(jj,j)
                              ijmin=min(jj,j)
                             if(ijmax.eq.ijmin) then
                               cz=1.0_dp
                             else
                               ij=(ijmax-1)*(ijmax-2)/2+ijmax-ijmin
                               cz=bes0(ij)
                             endif
                        Aj2=asv(jj,i,imn)
                        Bj2=bsv(jj,i,imn)
		       	do iuvp=iuvc,uvmax(i)
			   A=tbar(i,1,1,imn,iuvp)
                           Aj2=Aj2+A*tta(jj,i,iuvp)
                           B=tbar(i,1,2,imn,iuvp)
                           Aj2=Aj2+B*ttb(jj,i,iuvp)
                           A=tbar(i,2,1,imn,iuvp)
                           Bj2=Bj2+A*tta(jj,i,iuvp)
                           B=tbar(i,2,2,imn,iuvp)
                           Bj2=Bj2+B*ttb(jj,i,iuvp)
                        enddo
                        asc(j,i,imn)=asc(j,i,imn)+cz*Aj2
                        bsc(j,i,imn)=bsc(j,i,imn)+cz*Bj2
 6223                continue	            
                     A2=A2+conjg(asc(j,i,imn))*asv(j,i,imn)
                     A2=A2+conjg(bsc(j,i,imn))*bsv(j,i,imn)
                     B2=B2+conjg(asc(j,i,imn))*asc(j,i,imn)
                     B2=B2+conjg(bsc(j,i,imn))*bsc(j,i,imn)
 6222             continue
 6221          continue
 622        continue
            if(abs(B2).lt.1.0e-200_dp) then
               do i=1,nL
                  do imn=1,uvmax(i)
                     do j=1,nL
                        pct(j,i,1,iq,imn,iuv)=asr(j,i,imn)
                        pct(j,i,2,iq,imn,iuv)=bsr(j,i,imn)
                     enddo
                  enddo
               enddo  
               goto 1002
            endif
            Bj=A2/B2
            do 623 i=1,nL
               if(ind(i).gt.0) goto 623
               do imn=iuvc,uvmax(i)
                  do 6231 j=1,nL
                     if(iuv.gt.uvmax(j)) then
                        asp(j,i,imn)=ci0
                        bsp(j,i,imn)=ci0
                        go to 6231
                     endif   
                     asp(j,i,imn)=asv(j,i,imn)-Bj*asc(j,i,imn)
                     bsp(j,i,imn)=bsv(j,i,imn)-Bj*bsc(j,i,imn)
 6231             continue
               enddo
 623        continue
            do 624 i=1,nL
               c1i(i)=0.0_dp
               if(ind(i).gt.0) goto 624
               do 6241 imn=iuvc,uvmax(i)
                  do 6242 j=1,nL
                     if(iuv.gt.uvmax(j)) goto 6242
                     Aj2=Aj*as1(j,i,imn)+Bj*asv(j,i,imn)
                     Bj2=Aj*bs1(j,i,imn)+Bj*bsv(j,i,imn)
                     asr(j,i,imn)=asr(j,i,imn)+Aj2
                     bsr(j,i,imn)=bsr(j,i,imn)+Bj2
                     c1i(i)=c1i(i)+Aj2*conjg(Aj2)
                     c1i(i)=c1i(i)+Bj2*conjg(Bj2)
 6242             continue
 6241          continue
 624        continue
            cext0=0.0_dp
            cext1=0.0_dp
            do 625 i=1,nL
               if(ind(i).gt.0) goto 625
               cext0=cext0+c0i(i)
               cext1=cext1+c1i(i)
 625        continue
            temp=cext1/cext0
            if(iuv.eq.1.and.iq.eq.1) then
               if(niter/20*20.eq.niter.or.niter.eq.1) then
                  write(6,'(a11,i4,2x,e15.7)')   &
                    'iteration #',niter,temp
               endif
            endif
            if(temp.lt.small) then
               if(u.eq.v) then
                  write(6,'(3i4,a12,i4,2x,e15.5)') u,v,iq,   &
                     '  iteration ',niter,temp	         
               endif
               do i=1,nL
                  do imn=1,uvmax(i)
                     do j=1,nL
                        pct(j,i,1,iq,imn,iuv)=asr(j,i,imn)
                        pct(j,i,2,iq,imn,iuv)=bsr(j,i,imn)
                     enddo
                  enddo
               enddo  
               goto 1002
            endif
            if(niter.gt.MXINT) then
               write(6,'(a9,i5,i4,i3)') 'Caution: ',u,v,iq
               write(6,*) '*** Maximum iterations exceeded ***'
               write(6,*) '*** Solution may be inaccurate  ***'	
               goto 1002
            endif
            B2=0.0_dp
            do 626 i=1,nL 	   
               if(ind(i).gt.0) goto 626
               B2i(i)=0.0_dp
               do imn=iuvc,uvmax(i)
                  do 6261 j=1,nL
                     if(iuv.gt.uvmax(j)) goto 6261
                     Aj2=conjg(as0(j,i,imn))*asp(j,i,imn)
                     Bj2=conjg(bs0(j,i,imn))*bsp(j,i,imn)
                     B2i(i)=B2i(i)+Aj2
                     B2i(i)=B2i(i)+Bj2
 6261             continue
               enddo
 	       B2=B2+B2i(i)
 626        continue	
            A0=B0*Bj
            if(abs(A0).lt.1.0E-200_dp) then
               do i=1,nL
                  do imn=1,uvmax(i)
                     do j=1,nL
                        pct(j,i,1,iq,imn,iuv)=asr(j,i,imn)
                        pct(j,i,2,iq,imn,iuv)=bsr(j,i,imn)
                     enddo
                  enddo
               enddo  
               goto 1002
            endif
            A0=-Aj*B2/A0
            do 627 i=1,nL
               if(ind(i).gt.0) goto 627	   
               do imn=iuvc,uvmax(i)
                  do 6271 j=1,nL
                     if(iuv.gt.uvmax(j)) goto 6271
                     as1(j,i,imn)=as1(j,i,imn)-Bj*ast(j,i,imn)
                     bs1(j,i,imn)=bs1(j,i,imn)-Bj*bst(j,i,imn)
                     as1(j,i,imn)=asp(j,i,imn)-A0*as1(j,i,imn)
                     bs1(j,i,imn)=bsp(j,i,imn)-A0*bs1(j,i,imn)
 6271             continue
               enddo
 627        continue
            B0=B2
            do j=1,nL
               do i=1,nL
                  do imn=1,iuvc-1
                     as(i,imn)=ci0
                     bs(i,imn)=ci0
                  enddo
                  do imn=iuvc,uvmax(i)
                     as(i,imn)=as1(j,i,imn)
                     bs(i,imn)=bs1(j,i,imn)
                  enddo
               enddo
               call transT(nL,r0,nmax,uvmax,fint,atr0,btr0,    &
                   ek,drot,as,bs,as2,bs2,ind,confg,iuvc,1)
               do i=1,nL
                  do imn=iuvc,uvmax(i)
                     tta(j,i,imn)=as2(i,imn)
                     ttb(j,i,imn)=bs2(i,imn)
                  enddo
               enddo
            enddo
            A0=0.0_dp
            do 629 i=1,nL
               if(ind(i).gt.0) goto 629
               do 6291 imn=iuvc,uvmax(i)
                  n=dsqrt(dble(imn))
                  do 6292 j=1,nL
                     if(iuv.gt.uvmax(j)) goto 6292
                     ast(j,i,imn)=0._dp
                     bst(j,i,imn)=0._dp
                     do 6293 jj=1,nL
                        if(iuv.gt.uvmax(jj)) goto 6293
		   	   ijmax=max(jj,j)
                           ijmin=min(jj,j)
                        if(ijmax.eq.ijmin) then
                           cz=1._dp
                        else
                           ij=(ijmax-1)*(ijmax-2)/2+ijmax-ijmin
                           cz=bes0(ij)
                        endif
                        Aj2=as1(jj,i,imn)
                        Bj2=bs1(jj,i,imn)
                        do iuvp=iuvc,uvmax(i)
			   A=tbar(i,1,1,imn,iuvp)
                           Aj2=Aj2+A*tta(jj,i,iuvp)
                           B=tbar(i,1,2,imn,iuvp)
                           Aj2=Aj2+B*ttb(jj,i,iuvp)
                           A=tbar(i,2,1,imn,iuvp)
                           Bj2=Bj2+A*tta(jj,i,iuvp)
                           B=tbar(i,2,2,imn,iuvp)
                           Bj2=Bj2+B*ttb(jj,i,iuvp)
                        enddo
                        ast(j,i,imn)=ast(j,i,imn)+cz*Aj2
                        bst(j,i,imn)=bst(j,i,imn)+cz*Bj2
 6293                continue	 
                     A0=A0+conjg(as0(j,i,imn))*ast(j,i,imn)
                     A0=A0+conjg(bs0(j,i,imn))*bst(j,i,imn)
 6292             continue
 6291          continue
 629        continue
            if(abs(A0).lt.1.0E-200_dp) then
               do i=1,nL
                  do imn=1,uvmax(i)
                     do j=1,nL
                        pct(j,i,1,iq,imn,iuv)=asr(j,i,imn)
                        pct(j,i,2,iq,imn,iuv)=bsr(j,i,imn)
                     enddo
                  enddo
               enddo  
               goto 1002
            endif
            Aj=B0/A0
            niter=niter+1
            goto 62
 1002    continue
1001  continue

      end = DCLOCK()
      duration = end - start
      print*, "Loop at line: ", line+4, " executed in: ", duration, " microseconds."
            
!      C-----------------------------------------------------------------------
!C  calculating random-orientation averaged total and 
!C  individual-particle extinction cross-sections
!C-----------------------------------------------------------------------
 1800 do i=1,nL
         ind(i)=0
         cexti(i)=0.0_dp
         cscai(i)=0.0_dp
         cpri(i)=0.0_dp
      enddo
      cext=0.0_dp
      csca=0.0_dp
      cpr=0.0_dp

      n0=nmax0*(nmax0+2)
      as2 = CZERO
      bs2 = CZERO
      
      include 'tmatrix_mps_perf_loop1778.inc'

      
      line = __LINE__
      start = 0.0_dp
      end   = 0.0_dp
      duration = 0.0_dp
      start = DCLOCK()
      do iuv=1,n0
         do iq=1,2
            do 1801 j=1,nL
               if(iuv.gt.uvmax(j)) goto 1801
               do i=1,nL
                  do imn=1,uvmax(i)
                     as(i,imn)=pct(j,i,1,iq,imn,iuv)
                     bs(i,imn)=pct(j,i,2,iq,imn,iuv)
                  enddo
               enddo
               call transT(nL,r0,nmax,uvmax,fint,atr1,btr1,   &
                          ek,drot,as,bs,as2,bs2,ind,confg,1,1)
               do i=1,nL
                  do imn=1,uvmax(i)
                     A=as2(i,imn)+as(i,imn)
                     B=bs2(i,imn)+bs(i,imn)
                     pct(j,i,1,iq,imn,iuv)=A
                     pct(j,i,2,iq,imn,iuv)=B
                  enddo
               enddo
               cz=pct(j,j,iq,iq,iuv,iuv)
               cext=cext+cz
               cexti(j)=cexti(j)+cz
 1801	    continue
         enddo
      enddo
      end = DCLOCK()
      duration = end - start
      print*, "Loop at line: ", line+4, " executedin: ", duration, " microseconds."
!      C-----------------------------------------------------------------------
!C  calculating random-orientation averaged asymmetry parameter and  
!C  total and individual-particle scattering cross-sections 
      !C-----------------------------------------------------------------------
      

      
        include 'tmatrix_mps_perf_loop1828.inc' 
      

        line  = __LINE__
        start = 0.0_dp
        end   = 0.0_dp
        duration = 0.0_dp
        start = DCLOCK()
        do j=1,nL
          do iuv=1,uvmax(j)
            v=dsqrt(dble(iuv))
            u=iuv-v*v-v
            iuv1=iuv-2*u
            cb=(-1)**(v+u)
            jv1=v+1
            iuv2=-u+jv1*(jv1+1)
            iuv3=-u+(v-1)*v
            fuv1=v*jv1
            fuv1=dble(u)/fuv1
            fuv2=fnr(v)*fnr(v+2)/fnr(2*v+1)
            fuv2=fuv2/fnr(2*v+3)/dble(jv1)
            fuv2=fuv2*fnr(jv1-u)*fnr(jv1+u)
            fuv3=fnr(v-1)*fnr(jv1)/fnr(2*v-1)
            fuv3=fuv3/fnr(2*v+1)/dble(v)
            fuv3=fuv3*fnr(v-u)*fnr(v+u)
            juv1=iuv1-1
            juv2=iuv1+1
            juv3=iuv2-1
            juv4=iuv3-1
            juv5=iuv2+1
            juv6=iuv3+1
            guv=v*(v+1)
            guv1=-fnr(v-u)*fnr(v+u+1)/guv
            guv2=-fnr(v+u)*fnr(v-u+1)/guv
            guv=dble(v+1)*fnr(v)*fnr(v+2)*fnr(2*v+1)*fnr(2*v+3)
            guv3=u*(u+1)+((v-u)*(v+u+3)+(v+u)*(v-u+1))/2
            guv3=guv3*fnr(v+u+1)*fnr(v+u+2)/guv
            guv5=u*(u-1)+((v-u)*(v+u+1)+(v+u)*(v-u+3))/2
            guv5=-guv5*fnr(v-u+1)*fnr(v-u+2)/guv
            if(v.gt.1) then 
               guv=dble(v)*fnr(v-1)*fnr(v+1)*fnr(2*v-1)*fnr(2*v+1)
               guv4=u*(u+1)+((v-u-2)*(v+u+1)+(v+u)*(v-u+1))/2
               ntemp=v-u-1
               if(ntemp.lt.0) then
                  temp=0._dp
               else
                  temp=fnr(ntemp)
               endif
               guv4=-guv4*fnr(v-u)*temp/guv
               guv6=u*(u-1)+((v-u)*(v+u+1)+(v+u-2)*(v-u+1))/2
               ntemp=v+u-1
               if(ntemp.lt.0) then
                  temp=0._dp
               else
                  temp=fnr(ntemp)
               endif
               guv6=guv6*fnr(v+u)*temp/guv
            endif
            do i=1,nL
               do imn=1,uvmax(i)
                  n=dsqrt(dble(imn))
                  m=imn-n*n-n
                  imn1=-m+n*n+n
                  sb=(-1)**(m+n)
                  sb=cb*sb
                  n1=n+1
                  n2=2*n
                  rn=1._dp/dble(n*n1)
                  p=fnr(n)*fnr(n+2)/fnr(n2+1)
                  p=p/fnr(n2+3)/dble(n1)
                  t=fnr(n-1)*fnr(n+1)/fnr(n2-1)
                  t=t/fnr(n2+1)/dble(n)
                  rm=dble(m)*rn
                  imn2=(n+1)*(n+2)+m
                  fnp=fnr(n+m+1)*fnr(n-m+1)*p	
                  imn3=(n-1)*n+m
                  fn=fnr(n+m)*fnr(n-m)*t
                  jmn1=imn+1
                  jmn2=imn-1
                  jmn3=imn2+1
                  jmn4=imn3+1
                  jmn5=imn2-1
                  jmn6=imn3-1
                  gmn=n*(n+1)
                  gmn1=-fnr(n-m)*fnr(n+m+1)/gmn
                  gmn2=-fnr(n+m)*fnr(n-m+1)/gmn
                  gmn=dble(n+1)*fnr(n)*fnr(n+2)
                  gmn=gmn*fnr(2*n+1)*fnr(2*n+3)
                  gmn3=m*(m+1)+((n-m)*(n+m+3)+(n+m)*(n-m+1))/2
                  gmn3=gmn3*fnr(n+m+1)*fnr(n+m+2)/gmn
                  gmn5=m*(m-1)+((n-m)*(n+m+1)+(n+m)*(n-m+3))/2
                  gmn5=-gmn5*fnr(n-m+1)*fnr(n-m+2)/gmn
                  if(n.gt.1) then
                     gmn=dble(n)*fnr(n-1)*fnr(n+1)
                     gmn=gmn*fnr(2*n-1)*fnr(2*n+1)
                     gmn4=m*(m+1)+((n-m-2)*(n+m+1)+(n+m)*(n-m+1))/2
                     ntemp=n-m-1
                     if(ntemp.lt.0) then
                        temp=0.0_dp
                     else
                        temp=fnr(ntemp)
                     endif
                     gmn4=-gmn4*fnr(n-m)*temp/gmn
                     gmn6=m*(m-1)+((n-m)*(n+m+1)+(n+m-2)*(n-m+1))/2
                     ntemp=n+m-1
                     if(ntemp.lt.0) then
                        temp=0.0_dp
                     else
                        temp=fnr(ntemp)
                     endif
                     gmn6=gmn6*fnr(n+m)*temp/gmn
                  endif
                  do ip=1,2
                     do iq=1,2
                        A=pct(j,i,ip,iq,imn,iuv)
                        B=sb*pct(i,j,iq,ip,iuv1,imn1)
                        B=conjg(B)
                        cz=B*A
                        csca=csca+cz
                        cscai(j)=cscai(j)+cz	                  
                        A=rm*pct(j,i,3-ip,iq,imn,iuv)
                        if(n.eq.nmax(i)) goto 1951
                        A=A+fnp*pct(j,i,ip,iq,imn2,iuv)
 1951                   if(n.eq.1.or.iabs(m).gt.n-1) goto 1952
                        A=A+fn*pct(j,i,ip,iq,imn3,iuv)
 1952                   B=fuv1*pct(i,j,3-iq,ip,iuv1,imn1)
                        if(v.eq.nmax(j)) goto 1953
                        B=B-fuv2*pct(i,j,iq,ip,iuv2,imn1)
 1953                   if(v.eq.1.or.iabs(u).gt.v-1) goto 1954
                        B=B-fuv3*pct(i,j,iq,ip,iuv3,imn1)
 1954                   B=sb*conjg(B)
                        cz=B*A
                        cpr=cpr+cz
                        cpri(j)=cpri(j)+cz
                        A=ci0
                        if(iabs(m+1).gt.n) goto 1961
                        A=gmn1*pct(j,i,3-ip,iq,jmn1,iuv)
 1961                   if(n.eq.nmax(i)) goto 1962
                        A=A+gmn3*pct(j,i,ip,iq,jmn3,iuv)
 1962                   if(n.eq.1) goto 1963
                        if(iabs(m+1).gt.n-1) goto 1963
                        A=A+gmn4*pct(j,i,ip,iq,jmn4,iuv)
 1963                   B=ci0
                        if(iabs(u-1).gt.v) goto 1964
                        B=-guv2*pct(i,j,3-iq,ip,juv2,imn1)
 1964                   if(v.eq.nmax(j)) goto 1965
                        B=B+guv5*pct(i,j,iq,ip,juv5,imn1)
 1965                   if(v.eq.1) goto 1966
                        if(iabs(u-1).gt.v-1) goto 1966
                        B=B+guv6*pct(i,j,iq,ip,juv6,imn1)
 1966                   cz=conjg(B)*A
                        temp=0.5_dp*sb*cz
                        cpr=cpr+temp
                        cpri(j)=cpri(j)+temp
                        A=ci0
                        if(iabs(m-1).gt.n) goto 1971
                        A=gmn2*pct(j,i,3-ip,iq,jmn2,iuv)
 1971                   if(n.eq.nmax(i)) goto 1972
                        A=A+gmn5*pct(j,i,ip,iq,jmn5,iuv)
 1972                   if(n.eq.1) goto 1973
                        if(iabs(m-1).gt.n-1) goto 1973
                        A=A+gmn6*pct(j,i,ip,iq,jmn6,iuv)
 1973                   B=ci0
                        if(iabs(u+1).gt.v) goto 1974
                        B=-guv1*pct(i,j,3-iq,ip,juv1,imn1)
 1974                   if(v.eq.nmax(j)) goto 1975
                        B=B+guv3*pct(i,j,iq,ip,juv3,imn1)
 1975                   if(v.eq.1) goto 1976
                        if(iabs(u+1).gt.v-1) goto 1976
                        B=B+guv4*pct(i,j,iq,ip,juv4,imn1)
 1976                   cz=conjg(B)*A
                        temp=0.5_dp*sb*cz
                        cpr=cpr+temp
                        cpri(j)=cpri(j)+temp
                     enddo
                  enddo     
               enddo
            enddo
         enddo
      enddo
      end = DCLOCK()
      duration = end - start
      print*, "Loop at line: ", line+5, " executed in: ", duration, " microseconds."
      if(idscmt.lt.0) goto 2000
!C-----------------------------------------------------------------------
!C  calculating random-orientation averaged Mueller matrix elements
!C-----------------------------------------------------------------------
         write(6,'(/)')
         write(6,*) 'Calculating (Mueller) scattering matrix'
         write(6,'(/)') 
         n0=nmax0*(nmax0+2)
         nmax2=2*nmax0
         wmf1 = ZERO
         wm1  = ZERO

     
         include 'tmatrix_mps_perf_loop2039.inc'
    

         line = __LINE__
         start = 0.0_dp
         end   = 0.0_dp
         duration = 0.0_dp
         start = dclock()
         do n=1,nmax0
            do v=1,nmax0
               do m=0,n
                  it=n+v+1
                  call xuwigd(n,v,m,-1,w01s,wcf,it,nmf1)
                do ii=1,nmf1
                   it=it-1
                   wmf1(n,v,m,it)=w01s(ii)
                enddo
               it=n+v+1
               call xuwigd(n,v,m,1,w01s,wcf,it,nm1)
               do ii=1,nm1
                  it=it-1
                  wm1(n,v,m,it)=w01s(ii)
               enddo
            enddo
         enddo
      enddo
      end = dclock()
      duration = end - start
      print*, "Loop at line: ", line+4, " executed in: ", duration, " microseconds."
      do n=1,nmax0
         do jn=1,nmax0
            do ip=1,2
               do jp=1,2
                  B11n(n,jn,ip,jp)=ci0
                  B12n(n,jn,ip,jp)=ci0
                  B13n(n,jn,ip,jp)=ci0
                  B21n(n,jn,ip,jp)=ci0
                  B22n(n,jn,ip,jp)=ci0
                  B23n(n,jn,ip,jp)=ci0
                  do m=0,nmax0
                     A1m(m,n,jn,ip,jp)=ci0
                     A2m(m,n,jn,ip,jp)=ci0
                  enddo
               enddo
            enddo
         enddo
      enddo
      wsdt = ZERO

    
      include 'tmatrix_mps_perf_loop2092.inc'
      

      line = __LINE__
      start = 0.0_dp
      end   = 0.0_dp
      duration = 0.0_dp
      start = dclock()
      do 1900 ids=-nmax2,nmax2         
         do n=1,nmax0          
            do 19001 v=1,nmax0
               if(iabs(ids).gt.n+v) goto 19001
               fv0=0.5_dp*fnr(2*v+1)
               do 19002 ms=-n,n
                  mw=ms+ids
                  if(iabs(mw).gt.v) goto 19002
                  it=n+v+1
                  call xuwigd(n,v,ms,-mw,w01s,wcf,it,nvs)
                  fv=((-1)**ms)*fv0
                  do ii=1,nvs
                     it=it-1
                     wsdt(n,v,ms,it)=fv*w01s(ii)
                  enddo
19002          continue
19001       continue
         enddo
         
         do i=1,nL
            do j=1,nL
               do it=0,nmax2
                  do n=1,nmax0
                     do m=0,n
                        do ip=1,2
                           fhmf1(it,n,m,ip)=ci0
                           fmf1(it,n,m,ip)=ci0
                           fhm1v(it,n,m,ip)=ci0
                           fm1v(it,n,m,ip)=ci0
                           fhm1q(it,n,m,ip)=ci0
                           fm1q(it,n,m,ip)=ci0
                           fhmf1vq(it,n,m,ip)=ci0
                           fmf1vq(it,n,m,ip)=ci0
                        enddo
                     enddo
                  enddo
               enddo
               do 19003 n=1,nmax0
                  if(n.gt.nmax(i)) goto 19003
                  in0=n*(n+1)
                  cv=1._dp
                  do 19004 v=1,nmax0
                     if(v.gt.nmax(j)) goto 19004
                     if(iabs(ids).gt.n+v) goto 19004
                     cv=-cv
                     iv0=v*(v+1)
                     do it=0,n+v
                        do ip=1,2
                           do iq=1,2
                              fhas(it,ip,iq)=ci0
                              fnhs(it,ip,iq)=ci0
                           enddo
                        enddo
                     enddo
                     do 1901 ms=-n,n
                        mw=ms+ids
                        if(iabs(mw).gt.v) goto 1901
                        it=n+v+1
                        fn=(-1)**(n+v+ids)
                        isn=in0+ms
                        isf=in0-ms
                        iwv=iv0+mw
                        iwf=iv0-mw
                        nvs=n+v-max(iabs(n-v),iabs(ids))+1
                        do ii=1,nvs
                           it=it-1
                           do ip=1,2
                              do iq=1,2
                                 A=pct(i,j,iq,ip,iwf,isf)
                                 A=fn*wsdt(n,v,ms,it)*A
                                 fhas(it,ip,iq)=      &
                                    fhas(it,ip,iq)+A
                                 B=pct(j,i,ip,iq,isn,iwv)
                                 B=wsdt(n,v,ms,it)*B
                                 fnhs(it,ip,iq)=      &
                                    fnhs(it,ip,iq)+B
                              enddo
                           enddo
                        enddo
 1901                continue
                     itmin=max(iabs(ids),iabs(n-v))
                     do it=itmin,n+v
                        do m=0,n
                           if(iabs(m-1).gt.it) then
                              cwmf1=0._dp
                           else
                              cwmf1=wmf1(n,v,m,it)
                           endif
                           if(m+1.gt.it) then
                              cwm1=0._dp
                           else
                              cwm1=wm1(n,v,m,it)
                           endif
                           do ip=1,2
                              do iq=1,2
                                 cq=(-1)**iq
                                 A=cwmf1*fhas(it,ip,iq)
                                 A=conjg(A)
                                 B=cv*A
                                 fhmf1(it,n,m,ip)= &
                                 fhmf1(it,n,m,ip)+A
                                 fhmf1vq(it,n,m,ip)= &
                                 fhmf1vq(it,n,m,ip)+cq*B
                                 A=cwm1*fhas(it,ip,iq)
                                 A=conjg(A)
                                 B=cv*A
                                 fhm1v(it,n,m,ip)=   &
                                   fhm1v(it,n,m,ip)+B
                                 fhm1q(it,n,m,ip)=   &
                                    fhm1q(it,n,m,ip)+cq*A
                                 A=cwmf1*fnhs(it,ip,iq)
                                 B=cv*A
                                 fmf1(it,n,m,ip)=    &
                                   fmf1(it,n,m,ip)+A
                                 fmf1vq(it,n,m,ip)=   &
                                    fmf1vq(it,n,m,ip)+cq*B
                                 A=cwm1*fnhs(it,ip,iq)
                                 B=cv*A
                                 fm1v(it,n,m,ip)=    &
                                   fm1v(it,n,m,ip)+B
                                 fm1q(it,n,m,ip)=    &
                                   fm1q(it,n,m,ip)+cq*A
                              enddo
                           enddo
                        enddo
                     enddo
19004             continue
19003          continue
               do 1904 n=1,nmax0
                  if(n.gt.nmax(i)) goto 1904
                  itmin=iabs(ids)
                  if(itmin.gt.n+nmax0) goto 1904
                  do 1905 jn=1,nmax0
                     if(jn.gt.nmax(j)) goto 1905
                     if(itmin.gt.jn+nmax0) goto 1905
                     itmax=min(n,jn)+nmax0
                     do 1906 it=itmin,itmax
                        gt=2*it+1
                        do m=0,min(n,jn)
                           do ip=1,2
                              do jp=1,2
                                 A=fhmf1(it,n,m,ip)
                                 B=fmf1(it,jn,m,jp)
                                 A0=A*B
                                 A=fhm1q(it,n,m,ip)
                                 B=fm1q(it,jn,m,jp)
                                 A=gt*(A0+A*B)
                                 A1m(m,n,jn,ip,jp)=     &
                                 A1m(m,n,jn,ip,jp)+A
                                 fn=(-1)**(ip+jp+n+jn)
                                 A=fhmf1vq(it,n,m,ip)
                                 B=fmf1vq(it,jn,m,jp)
                                 A0=A*B
                                 A=fhm1v(it,n,m,ip)
                                 B=fm1v(it,jn,m,jp)
                                 A=fn*gt*(A0+A*B)
                                 A2m(m,n,jn,ip,jp)=     &
                                 A2m(m,n,jn,ip,jp)+A
                                 if(m.gt.0) goto 19061
                                 temp=((-1)**ip)*gt
                                 A=fhm1q(it,n,0,ip)
                                 B=fmf1(it,jn,2,jp)
                                 A=temp*A*B
                                 B11n(n,jn,ip,jp)=    &
                                 B11n(n,jn,ip,jp)+A
                                 nj=ip+n+it
                                 temp=((-1)**nj)*gt
                                 A=fhmf1vq(it,n,1,ip)
                                 B=fmf1(it,jn,1,jp)
                                 A=temp*A*B
                                 B12n(n,jn,ip,jp)=    &
                                 B12n(n,jn,ip,jp)+A
                                 nj=ip+n+jn
                                 temp=((-1)**nj)*gt
                                 A=fhmf1vq(it,n,2,ip)
                                 B=fm1v(it,jn,0,jp)
                                 A=temp*A*B
                                 B13n(n,jn,ip,jp)=    &
                                 B13n(n,jn,ip,jp)+A
                                 nj=jp+n+jn
                                 temp=((-1)**nj)*gt
                                 A=fhm1v(it,n,0,ip)
                                 B=fmf1vq(it,jn,2,jp)
                                 A=temp*A*B
                                 B21n(n,jn,ip,jp)=   &
                                 B21n(n,jn,ip,jp)+A
                                 nj=jp+jn+it
                                 temp=((-1)**nj)*gt
                                 A=fhmf1(it,n,1,ip)
                                 B=fmf1vq(it,jn,1,jp)
                                 A=temp*A*B
                                 B22n(n,jn,ip,jp)=   &
                                 B22n(n,jn,ip,jp)+A
                                 temp=((-1)**jp)*gt
                                 A=fhmf1(it,n,2,ip)
                                 B=fm1q(it,jn,0,jp)
                                 A=temp*A*B
                                 B23n(n,jn,ip,jp)=   &
                                 B23n(n,jn,ip,jp)+A
19061                            continue
                              enddo
                           enddo
                        enddo
 1906                continue
 1905             continue
 1904          continue
            enddo
         enddo
1900  continue
      end =  dclock()
      duration = end - start
      print*, "Loop at line: ", line+4, " executed in: ", duration, " microseconds."

      include 'tmatrix_mps_perf_loop2323.inc'
      
      line = __LINE__
      start = 0.0_dp
      end   = 0.0_dp
      duration = 0.0_dp
      start = dclock()
      do 1910 ia=1,nang
         iang=2*nang-ia
         dang(ia)=sang*dble(ia-1)
         dang(iang)=180.0_dp-dang(ia)                  
         theta=dang(ia)*pione/180.0_dp
         xt=dcos(theta)
         call tipitaud(nmax0,xt)
         do i=1,4
            do j=1,4
               mue(i,j,ia)=0.0_dp
               mue(i,j,iang)=0.0_dp
            enddo
         enddo
         do i=1,2
            do j=1,2
               do ik=1,2
                  do jk=1,2
                     A0p(i,j,ik,jk)=ci0
                     A1p(i,j,ik,jk)=ci0
                     A0pg(i,j,ik,jk)=ci0
                     A1pg(i,j,ik,jk)=ci0
                     B0p(i,j,ik,jk)=ci0
                     B1p(i,j,ik,jk)=ci0
                     B0pg(i,j,ik,jk)=ci0
                     B1pg(i,j,ik,jk)=ci0
                  enddo
               enddo
            enddo
         enddo
         do n=1,nmax0
            gn=(-1)**n
            itau0=(n-1)*(n+2)/2+1
            tau0p(1)=0.5_dp*tau(itau0)
            tau0p(2)=0._dp
            tau0pg(1)=-gn*tau0p(1)
            tau0pg(2)=0._dp
            itau=itau0+1
            gn=-gn
            tau1p(1)=tau(itau)
            tau1p(2)=pi(itau)
            tau1pg(1)=-gn*tau1p(1)
            tau1pg(2)=gn*tau1p(2)
            itau=itau0+2
            gn=-gn
            tau2p(1)=tau(itau)
            tau2p(2)=pi(itau)
            tau2pg(1)=-gn*tau2p(1)
            tau2pg(2)=gn*tau2p(2)
            do jn=1,nmax0
               jtau0=(jn-1)*(jn+2)/2+1
               gn=(-1)**jn
               tau0pj(1)=0.5_dp*tau(jtau0)
               tau0pj(2)=0._dp
               tau0pjg(1)=-gn*tau0pj(1)
               tau0pjg(2)=0._dp
               jtau=jtau0+1
               gn=-gn
               tau1pj(1)=tau(jtau)
               tau1pj(2)=pi(jtau)
               tau1pjg(1)=-gn*tau1pj(1)
               tau1pjg(2)=gn*tau1pj(2)
               jtau=jtau0+2
               gn=-gn
               tau2pj(1)=tau(jtau)
               tau2pj(2)=pi(jtau)
               tau2pjg(1)=-gn*tau2pj(1)
               tau2pjg(2)=gn*tau2pj(2)
               do ik=1,2
                  do jk=1,2
                     tau20(ik,jk)=tau2p(ik)*tau0pj(jk)
                     tau11(ik,jk)=tau1p(ik)*tau1pj(jk)
                     tau02(ik,jk)=tau0p(ik)*tau2pj(jk)
                     tau20g(ik,jk)=tau2pg(ik)*tau0pjg(jk)
                     tau11g(ik,jk)=tau1pg(ik)*tau1pjg(jk)
                     tau02g(ik,jk)=tau0pg(ik)*tau2pjg(jk)
                  enddo
               enddo
               do ip=1,2
                  do jp=1,2
                     A=B11n(n,jn,ip,jp)
                     B=B12n(n,jn,ip,jp)
                     cmz=B13n(n,jn,ip,jp)
                     Aj=B21n(n,jn,ip,jp)
                     Bj=B22n(n,jn,ip,jp)
                     cmzj=B23n(n,jn,ip,jp)
                     do ik=1,2
                        do jk=1,2
                           A1=A*tau02(ik,jk)
                           A2=B*tau11(ik,jk)
                           Aj2=cmz*tau20(ik,jk)
                           B1=A1-A2+Aj2                           
                           A1=Aj*tau02(ik,jk)
                           A2=Bj*tau11(ik,jk)
                           Aj2=cmzj*tau20(ik,jk)
                           B2=A1-A2+Aj2
                           A0=B1+B2
                           B0p(ip,jp,ik,jk)=B0p(ip,jp,ik,jk)+A0
                           A0=B1-B2
                           B1p(ip,jp,ik,jk)=B1p(ip,jp,ik,jk)+A0
                           A1=A*tau02g(ik,jk)
                           A2=B*tau11g(ik,jk)
                           Aj2=cmz*tau20g(ik,jk)
                           B1=A1-A2+Aj2
                           A1=Aj*tau02g(ik,jk)
                           A2=Bj*tau11g(ik,jk)
                           Aj2=cmzj*tau20g(ik,jk)
                           B2=A1-A2+Aj2
                           A0=B1+B2
                           B0pg(ip,jp,ik,jk)=B0pg(ip,jp,ik,jk)+A0
                           A0=B1-B2
                           B1pg(ip,jp,ik,jk)=B1pg(ip,jp,ik,jk)+A0
                        enddo
                     enddo
                  enddo
               enddo
               do m=0,min(n,jn)
                  itau=itau0+m
                  gmn=(-1)**(m+n)
                  taup(1)=tau(itau)
                  taup(2)=pi(itau)
                  if(m.eq.0) then
                     taup(1)=0.5_dp*taup(1)
                     taup(2)=0._dp
                  endif
                  taupg(1)=-gmn*taup(1)
                  taupg(2)=gmn*taup(2)
                  gmnj=(-1)**(jn+m)
                  jtau=jtau0+m
                  taupj(1)=tau(jtau)
                  taupj(2)=pi(jtau)
                  if(m.eq.0) then
                     taupj(1)=0.5_dp*taupj(1)
                     taupj(2)=0._dp
                  endif
                  taupjg(1)=-gmnj*taupj(1)
                  taupjg(2)=gmnj*taupj(2)
                  do ik=1,2
                     do jk=1,2
                        taum(ik,jk)=taup(ik)*taupj(jk)
                        taumg(ik,jk)=taupg(ik)*taupjg(jk)
                     enddo
                  enddo            
                  do ip=1,2
                     do jp=1,2
                        A=A1m(m,n,jn,ip,jp)
                        Aj=A
                        B=A2m(m,n,jn,ip,jp)
                        A=A+B
                        Aj=Aj-B
                        do ik=1,2
                           do jk=1,2
                              B0=A*taum(ik,jk)
                              B1=Aj*taum(ik,jk)
                              A0p(ip,jp,ik,jk)=   &
                              A0p(ip,jp,ik,jk)+B0
                              A1p(ip,jp,ik,jk)=    &
                              A1p(ip,jp,ik,jk)+B1
                              B0=A*taumg(ik,jk)
                              B1=Aj*taumg(ik,jk)
                              A0pg(ip,jp,ik,jk)=   &
                              A0pg(ip,jp,ik,jk)+B0
                              A1pg(ip,jp,ik,jk)=   &
                              A1pg(ip,jp,ik,jk)+B1
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo         
         do ip=1,2
            do jp=1,2
               temp=A0p(ip,jp,ip,jp)+A0p(ip,jp,3-ip,3-jp)
               mue(1,1,ia)=mue(1,1,ia)+temp
               temp=B0p(ip,jp,3-ip,3-jp)-B0p(ip,jp,ip,jp)
               mue(1,2,ia)=mue(1,2,ia)+temp
               A0=A1p(ip,jp,ip,jp)+A1p(ip,jp,3-ip,3-jp)
               A0=A0-B1p(ip,jp,ip,jp)+B1p(ip,jp,3-ip,3-jp)
               temp=-imag(A0)
               mue(1,3,ia)=mue(1,3,ia)+temp
               temp=-A0
               mue(1,4,ia)=mue(1,4,ia)+temp
               temp=A0p(ip,jp,ip,jp)-A0p(ip,jp,3-ip,3-jp)
               mue(2,1,ia)=mue(2,1,ia)+temp
               temp=B0p(ip,jp,ip,jp)+B0p(ip,jp,3-ip,3-jp)
               mue(2,2,ia)=mue(2,2,ia)-temp
               A0=A1p(ip,jp,ip,jp)+A1p(ip,jp,3-ip,3-jp)
               A0=A0+B1p(ip,jp,ip,jp)-B1p(ip,jp,3-ip,3-jp)
               temp=-imag(A0)
               mue(2,3,ia)=mue(2,3,ia)+temp
               temp=-A0
               mue(2,4,ia)=mue(2,4,ia)+temp
               temp=0.5_dp*imag(A1p(ip,jp,ip,3-jp))
               mue(3,1,ia)=mue(3,1,ia)+temp
               temp=-0.5_dp*imag(B1p(ip,jp,ip,3-jp))
               mue(3,2,ia)=mue(3,2,ia)+temp
               A0=A0p(ip,jp,ip,3-jp)-A0p(ip,jp,3-ip,jp)
               A0=A0-B0p(ip,jp,ip,3-jp)-B0p(ip,jp,3-ip,jp)
               temp=A0
               mue(3,3,ia)=mue(3,3,ia)+temp
               temp=-imag(A0)
               mue(3,4,ia)=mue(3,4,ia)+temp
               temp=-0.5_dp*A1p(ip,jp,ip,3-jp)
               mue(4,1,ia)=mue(4,1,ia)+temp
               temp=0.5_dp*B1p(ip,jp,ip,3-jp)
               mue(4,2,ia)=mue(4,2,ia)+temp
               A0=A0p(ip,jp,ip,3-jp)+A0p(ip,jp,3-ip,jp)
               A0=A0-B0p(ip,jp,ip,3-jp)+B0p(ip,jp,3-ip,jp)
               temp=imag(A0)
               mue(4,3,ia)=mue(4,3,ia)+temp
               temp=A0
               mue(4,4,ia)=mue(4,4,ia)+temp
               if(ia.eq.iang) goto 19101
               temp=A0pg(ip,jp,ip,jp)+A0pg(ip,jp,3-ip,3-jp)
               mue(1,1,iang)=mue(1,1,iang)+temp
               temp=B0pg(ip,jp,3-ip,3-jp)-B0pg(ip,jp,ip,jp)
               mue(1,2,iang)=mue(1,2,iang)+temp
               A0=A1pg(ip,jp,ip,jp)+A1pg(ip,jp,3-ip,3-jp)
               A0=A0-B1pg(ip,jp,ip,jp)+B1pg(ip,jp,3-ip,3-jp)
               temp=-imag(A0)
               mue(1,3,iang)=mue(1,3,iang)+temp
               temp=-A0
               mue(1,4,iang)=mue(1,4,iang)+temp
               temp=A0pg(ip,jp,ip,jp)-A0pg(ip,jp,3-ip,3-jp)
               mue(2,1,iang)=mue(2,1,iang)+temp
               temp=B0pg(ip,jp,ip,jp)+B0pg(ip,jp,3-ip,3-jp)
               mue(2,2,iang)=mue(2,2,iang)-temp
               A0=A1pg(ip,jp,ip,jp)+A1pg(ip,jp,3-ip,3-jp)
               A0=A0+B1pg(ip,jp,ip,jp)-B1pg(ip,jp,3-ip,3-jp)
               temp=-imag(A0)
               mue(2,3,iang)=mue(2,3,iang)+temp
               temp=-A0
               mue(2,4,iang)=mue(2,4,iang)+temp
               temp=0.5_dp*imag(A1pg(ip,jp,ip,3-jp))
               mue(3,1,iang)=mue(3,1,iang)+temp
               temp=-0.5_dp*imag(B1pg(ip,jp,ip,3-jp))
               mue(3,2,iang)=mue(3,2,iang)+temp
               A0=A0pg(ip,jp,ip,3-jp)-A0pg(ip,jp,3-ip,jp)
               A0=A0-B0pg(ip,jp,ip,3-jp)-B0pg(ip,jp,3-ip,jp)
               temp=A0
               mue(3,3,iang)=mue(3,3,iang)+temp
               temp=-imag(A0)
               mue(3,4,iang)=mue(3,4,iang)+temp
               temp=-0.5_dp*A1pg(ip,jp,ip,3-jp)
               mue(4,1,iang)=mue(4,1,iang)+temp
               temp=0.5_dp*B1pg(ip,jp,ip,3-jp)
               mue(4,2,iang)=mue(4,2,iang)+temp
               A0=A0pg(ip,jp,ip,3-jp)+A0pg(ip,jp,3-ip,jp)
               A0=A0-B0pg(ip,jp,ip,3-jp)+B0pg(ip,jp,3-ip,jp)
               temp=imag(A0)
               mue(4,3,iang)=mue(4,3,iang)+temp
               temp=A0
               mue(4,4,iang)=mue(4,4,iang)+temp
19101          continue
            enddo
         enddo
         temp=mue(1,1,ia)+mue(1,2,ia)+mue(2,1,ia)+mue(2,2,ia)
         i22(ia)=0.5_dp*temp
         temp=mue(1,1,ia)-mue(1,2,ia)-mue(2,1,ia)+mue(2,2,ia)
         i11(ia)=0.5_dp*temp
         temp=mue(1,1,ia)+mue(1,2,ia)-mue(2,1,ia)-mue(2,2,ia)
         i21(ia)=0.5_dp*temp
         temp=mue(1,1,ia)-mue(1,2,ia)+mue(2,1,ia)-mue(2,2,ia)
         i12(ia)=0.5_dp*temp
         if(ia.eq.iang) goto 1910
         temp=mue(1,1,iang)+mue(1,2,iang)+mue(2,1,iang)+mue(2,2,iang)
         i22(iang)=0.5_dp*temp
         temp=mue(1,1,iang)-mue(1,2,iang)-mue(2,1,iang)+mue(2,2,iang)
         i11(iang)=0.5_dp*temp
         temp=mue(1,1,iang)+mue(1,2,iang)-mue(2,1,iang)-mue(2,2,iang)
         i21(iang)=0.5_dp*temp
         temp=mue(1,1,iang)-mue(1,2,iang)+mue(2,1,iang)-mue(2,2,iang)
         i12(iang)=0.5_dp*temp
1910 continue
      end = dclock()
      duration = end - start
      print*, "Loop at line: ", line+4, " executed in: ", duration, " microseconds."
      cbak=i11(2*nang-1)
      do i=1,nang2
         inat(i)=i11(i)+i22(i)+i12(i)+i21(i)
         pol(i)=(i11(i)+i12(i)-i22(i)-i21(i))/inat(i)
      enddo
 2000 cz=twopi/k**2
      csca=csca*cz
      cext=cext*cz
      cabs=cext-csca
      cpr=cpr*cz
      assym=cpr/csca
      cbak=2.d0*cbak*cz
      write(6,'(5x,a6,7x,a6,7x,a6,7x,a6,7x,a5,6x,a12)')    &
      '<Cext>','<Cabs>','<Csca>','<Cbak>','<Cpr>','<cos(theta)>'
      write(6,'(2x,6e13.5)') cext,cabs,csca,cbak,cext-cpr,assym
      cscax=0.0_dp
      cextx=0.0_dp
      cprx=0.0_dp
      do i=1,nL
         cscax=cscax+cscai(i)
         cextx=cextx+cexti(i)
         cprx=cprx+cpri(i)
      enddo
      cscax=cscax*cz
      cextx=cextx*cz
      cprx=cprx*cz
      cabsx=cextx-cscax
      assymx=cprx/cscax
      cbakx=cbak	
      write(6,'(2x,6e13.5)')   &
          cextx,cabsx,cscax,cbakx,cextx-cprx,assymx
      do i=1,nL
         cexti(i)=cexti(i)*cz
         cscai(i)=cscai(i)*cz
         cabsi(i)=cexti(i)-cscai(i)
         cpri(i)=cpri(i)*cz
         assymi(i)=cpri(i)/csca
         write(6,'(i5,5e15.6)') i,cexti(i),cabsi(i),   &
                              cscai(i),cpri(i),assymi(i)
      enddo
!c	flout='cr'//fileout
!c	open(33,file=flout,status='unknown')
!c	write(33,'(a20,a47)') 
!c     +     flout,'(Total and individual-particle cross sections)'
!c	write(33,'(a32,2x,a22)') 
!c     +     'input sphere-aggregate filename:',FLNAME
!c	write(33,'(12x,a4,11x,a4,11x,a4,11x,a3,8x,a12)')
!c     +    'Cext','Cabs','Csca','Cpr','<cos(theta)>'
!c	write(33,'(a5,5e15.6)') 'total',cext,cabs,csca,cext-cpr,
!c     +     assym
!c	do i=1,nL
!c	   write(33,'(i5,5e15.6)') i,cexti(i),cabsi(i),cscai(i),
!c     +           cpri(i),assymi(i)
!c	enddo
!c	close(33)       
      cz=pione*gcvr*gcvr
      cscav=csca/cz
      cextv=cext/cz
      cprv=cpr/cz
      cbakv=cbak/cz
      cabsv=cextv-cscav
      temp=gcvr*gcvr/gcs
      cscas=cscav*temp
      cexts=cextv*temp
      cabss=cabsv*temp
      cbaks=cbakv*temp
      cprs=cprv*temp
 222  format(6e13.5)
 221  format(3x,a7,6x,a7,6x,a7,6x,a7,6x,a6,5x,a12) 
      if(irat.eq.1) then        
         write(6,221)    &
           '<Qextv>','<Qabsv>','<Qscav>','<Qbakv>','<Qprv>',   &
           '<cos(theta)>'	
         write(6,222) cextv,cabsv,cscav,cbakv,cprv,assym
      else
         write(6,221)                                         &
           '<Qexts>','<Qabss>','<Qscas>','<Qbaks>','<Qprs>', &
           '<cos(theta)>'
         write(6,222) cexts,cabss,cscas,cbaks,cprs,assym
      endif
!      open(12,file=fileout,status='unknown')
!      write(12,'(/)')
      if(irat.eq.1) then
         write(6,'(1x,a5,f8.3)')  'xv:',xv     
                
      else
         write(6,'(1x,a5,f8.3)')   'xs:',xs         
                 
      endif
      if(idscmt.lt.0) then
         write(6,'(/)')
         write(12,*) '*** backscattering and scatttering matrix',  &
                    ' are not calculated ***'
      endif
      write(6,'(/)')
      write(6,221)                                                  &
         '<Cext>','<Cabs>','<Csca>','<Cbak>','<Cpr>','<cos(theta)>' 
      write(6,222)                                                  &
         cext,cabs,csca,cbak,cext-cpr,assym
      if(irat.eq.1) then 
         write(6,221)                                              &
           '<Qextv>','<Qabsv>','<Qscav>','<Qbakv>','<Qprv>',       &
           '<cos(theta)>'
         write(6,222) cextv,cabsv,cscav,cbakv,cprv,assym
      else
         write(6,221)                                             &
           '<Qexts>','<Qabss>','<Qscas>','<Qbaks>','<Qprs>',      &
            '<cos(theta)>'
         write(6,222) cexts,cabss,cscas,cbaks,cprs,assym
      endif
      !if(idscmt.lt.0) goto 2001
      write(6,'(/)')
      write(6,'(2x,a4,4x,a7,4x,a6,4x,a7,6x,a7,6x,a7,6x,a7)')       &
        's.a.','<total>','<pol.>','<S1*S1>','<S4*S4>','<S3*S3>',   &
        '<S2*S2>'               
      do i=1,nang2
         write(12,'(f6.1,e13.5,f8.4,4e13.5)')                      &
           dang(i),inat(i),pol(i),i11(i),i21(i),i12(i),i22(i) 
      enddo
      write(6,'(/)')	
      write(6,'(1x,a50)')                                          &
         'Scattering matrix (4X4 for each scattering angle):'
      do i=1,nang2
         write(6,'(f7.1,4e16.7)')                                  &
           dang(i),mue(1,1,i),mue(1,2,i),mue(1,3,i),mue(1,4,i)
         write(6,'(7x,4e16.7)')                                    &
                   mue(2,1,i),mue(2,2,i),mue(2,3,i),mue(2,4,i)
         write(6,'(7x,4e16.7)')                                    &
                   mue(3,1,i),mue(3,2,i),mue(3,3,i),mue(3,4,i)
         write(6,'(7x,4e16.7)')                                    &
                   mue(4,1,i),mue(4,2,i),mue(4,3,i),mue(4,4,i)
      enddo
9999  print*, "*****FATAL-ERROR*****"
      print*, emsg
      ERROR STOP 
    end subroutine
                                  
                                  
!    C
!C  subroutine scoatabd.f 
!C  slightly edited from SCSERD.FOR BY R.T. Wang 
!C  SUB SCSERD: DP CORE-MANTLE-SPHERE MIE COEFFICIENTS 
!C  6-7-91   Ru T. Wang
!C  ED. OF SCSE1W(1974(UNIVAC) & 10-84(VAX)) & W.X.WANG(5-7-85).
!C  DOUBLE PRECISION,DOWNWARD RECURSION RATIO FUNCTION ALGORITHMS.
!C  Q=(CORE RADIUS)/(MANTLE RADIUS); XB=SIZE PARAMETER OF MANTLE
!C  COMPLEX REFRACTIVE INDEX OF CORE = (XM1,YM1)
!C  COMPLEX REFRACTIVE INDEX OF MANTLE = (XM2,YM2)
!C
    subroutine scoatabd(XB,Q,XM1,YM1,XM2,YM2,np,an,bn,NADD,NSTOP)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: scoatabd
          real(kind=dp) :: XB,Q,XM1,YM1,XM2,YM2
          integer(kind=int4) :: np
          complex(16), dimension(np) :: an,bn
          integer(kind=int4) :: NADD,NSTOP
          ! Locals
          integer(kind=int4),parameter :: nab = 500
          integer(kind=int4), parameter :: ndx = 5000
          real(kind=dp), dimension(nab) :: AR,AI,BR,BI
!DIR$     ATTRIBUTES ALIGN : 64 :: AR,AI,BR,BI
          real(kind=dp), dimension(ndx) ::  AM1AR,AM1AI,AM2AR,AM2AI,   &
                                         AM2BR,AM2BI,AB,SM2AR,      &
                                         SM2AI,SM2BR,SM2BI,SB,      &
                                         BM2AR,BM2AI,BM2BR,BM2BI,   &
                                         BDBR,BDBI,BB,CM2AR,        &
                                         CM2AI,CM2BR,CM2BI,CB
!DIR$     ATTRIBUTES ALIGN : 64 :: AM2BR,AM2BI,AB,SM2AR,SM2AI
!DIR$     ATTRIBUTES ALIGN : 64 :: SM2BR,SM2BI,SB,BM2AR,BM2AI
!DIR$     ATTRIBUTES ALIGN : 64 :: BM2BR,BM2BI,BDBR,BDBI,BB
!DIR$     ATTRIBUTES ALIGN : 64 :: CM2AR,CM2AI,CM2BR,CM2BI,CB
          real(kind=dp), dimension(4) ::     U,V,D1,EV,SHV,CHV,SU,     &
                                          CU,SN1R,SN1I,CN1R,CN1I
          real(kind=dp) :: YM1P,YM2P,oneth,CTST,CNX,FCT0,SM1A0R,SM1A0I,     &
                        SM2A0R,SM2A0I,SM2B0R,SM2B0I,SB0,CNN,FCT, SM1A1R, &
                        SM1A1I,SM2A1R,SM2A1I,QSM1A,QSM2A,QSM2B,BM2A0R,   &
                        BM2A0I, QSM2B0, BM2B0R,BM2B0I,BB0,BDB0R, BDB0I,  &
                        UM2R0,UM2I0,VM2R0,VM2I0,CN,CM2A0R,QCM2A,DBDB0R,  &
                        DBDB0I,QCM2B,QDB,CM2B1R,CM2B1I,CB1,DBDB1R,DBDB1I,&
                        SNB,CNB,SSB1R,SSB1I,SSB2R,SSB2I,SZB2R,SZB2I,ANNR,&
                        ANNI,ANDR,ANDI,BNNR,BNNI,BNDR,BNDI,UM2R1,UM2I1,  &
                        UM2R,UM2I,VM2RI,VM2I1,VM2R,VM2I,UVR,UVI,SSA1R,   &
                        SSA1I,CSB1R,CSB1I,SCS1R,SCS1I,UVS1R,UVS1I,CSA1R, &
                        CSA1I,CSS1R,CSS1I,CZB1R,CZB1I,SCZ1R,SCZ1I,UVZ1R, &
                        UVZ1I,CSZ1R,CSZ1I,SSA2R,SSA2I,CSB2R,CSB2I,SCS2R, &
                        SCS2I,UVS2R,UVS2I,CSA2R,CSA2I,CSS2R,CSS2I,CZB2R, &
                        CZB2I,SCZ2R,SCZ2I,UVZ2R,UVZ2I,CSZ2R,CSZ2I,AND,   &
                        ABANDR,ABANDI,ABANNR,ABANNI,AAAA,ABBNDR, ABBNDI, &
                        ABBNNR,ABBNNI,BBBB,TI,YM,XM,XN,XA,QCM2A0,QCM2B0, &
                        CM2A0I,SZB1R,SZB1I,VM2R1,BND
          integer(kind=int4) :: K,NX,I,J,N,ii
          ! Exec code....
!DIR$     VECTOR ALIGNED
!DIR$     SIMD VECTORLENGTHFOR(REAL(KIND=8))
          do ii = 1, nab
              AR(ii) = ZERO
              AI(ii) = ZERO
              BR(ii) = ZERO
              BI(ii) = ZERO
          end do
!DIR$     VECTOR ALIGNED
!DIR$     SIMD VECTORLENGTHFOR(REAL(KIND=8))          
          do ii = 1, ndx
             AM1AR(ii) = ZERO
             AM1AI(ii) = ZERO
             AM2AR(ii) = ZERO
             AM2AI(ii) = ZERO
             AM2BR(ii) = ZERO
             AM2BI(ii) = ZERO
             AB(ii)    = ZERO
             SM2AR(ii) = ZERO
             SM2AI(ii) = ZERO
             SM2BR(ii) = ZERO
             SM2BI(ii) = ZERO
             SB(ii)    = ZERO
             BM2AR(ii) = ZERO
             BM2AI(ii) = ZERO
             BM2BR(ii) = ZERO
             BM2BI(ii) = ZERO
             BDBR(ii)  = ZERO
             BDBI(ii)  = ZERO
             BB(ii)    = ZERO
             CM2AR(ii) = ZERO
             CM2AI(ii) = ZERO
             CM2BR(ii) = ZERO
             CM2BI(ii) = ZERO
             CB(ii)    = ZERO
          end do
           do i=1,np
              an(i)=0.0_dp
              bn(i)=0.0_dp
          enddo
          XM=DMAX1(XM1,XM2)
          YM1P=DABS(YM1)
          YM2P=DABS(YM2)
          YM=DMAX1(YM1P,YM2P)
          XN=XB*DSQRT(XM**2+YM**2)
          NX=1.1_dp*XN+10._dp
      if(NX.gt.ndx) then
         write(6,*) 'parameter (ndx) in sub. scoatabd too small'
         write(6,*) 'please change ndx to ',NX
         write(6,*) 'recompile, and then try again'
         call TRACEBACKQQ(STRING="parameter (ndx) in sub. scoatabd too small",USER_EXIT_CODE = -1)
         stop
      endif
      oneth=0.3333333333333333333333333333_dp
      NSTOP=XB+4._dp*XB**oneth
      NSTOP=NSTOP+2+NADD   
      if(NSTOP.gt.ndx) then 
         write(6,*) 'particle size too large'
         call TRACEBACKQQ(STRING="particle size too large",USER_EXIT_CODE = -1)
         stop   
      endif   	       
      XA=XB*Q
      CTST=0.00000000000001_dp           
      U(1)=XM1*XA
      V(1)=YM1*XA
      U(2)=XM2*XA
      V(2)=YM2*XA
      U(3)=XM2*XB
      V(3)=YM2*XB
      U(4)=XB
      V(4)=0.0_dp
      K=1
      IF(Q.EQ.0.0_dp) K=3
      DO 10 J=K,4
         D1(J)=U(J)*U(J)+V(J)*V(J)
         EV(J)=DEXP(V(J))
         SHV(J)=0.5_dp*(EV(J)*EV(J)-1.0_dp)
         CHV(J)=0.5_dp*(EV(J)*EV(J)+1.0_dp)
         SU(J)=DSIN(U(J))
         CU(J)=DCOS(U(J))
         SN1R(J)=SU(J)*CHV(J)
         SN1I(J)=CU(J)*SHV(J)
         CN1R(J)=CU(J)*CHV(J)
10    CN1I(J)=-SU(J)*SHV(J)
      CNX=NX
      FCT0=2.0_dp*CNX+3.0_dp
      IF(Q.EQ.0.0_dp) GO TO 12
      SM1A0R=U(1)/FCT0
      SM1A0I=V(1)/FCT0
      SM2A0R=U(2)/FCT0
      SM2A0I=V(2)/FCT0
12    SM2B0R=U(3)/FCT0
      SM2B0I=V(3)/FCT0
      SB0=U(4)/FCT0
      DO 15 I=1,NX
         N=NX-I+1
         CNN=N
         FCT=2.0_dp*CNN+1.0_dp
         IF(Q.EQ.0.0D0) GO TO 18
         SM1A1R=+U(1)*FCT/D1(1)-SM1A0R
         SM1A1I=-V(1)*FCT/D1(1)-SM1A0I
         SM2A1R=+U(2)*FCT/D1(2)-SM2A0R
         SM2A1I=-V(2)*FCT/D1(2)-SM2A0I
         AM1AR(N)=-CNN*U(1)/D1(1)+SM1A1R
         AM1AI(N)=+CNN*V(1)/D1(1)+SM1A1I
         AM2AR(N)=-CNN*U(2)/D1(2)+SM2A1R
         AM2AI(N)=+CNN*V(2)/D1(2)+SM2A1I
         QSM1A=SM1A1R**2+SM1A1I**2
         QSM2A=SM2A1R**2+SM2A1I**2
         SM1A0R=+SM1A1R/QSM1A
         SM1A0I=-SM1A1I/QSM1A
         SM2A0R=+SM2A1R/QSM2A
         SM2A0I=-SM2A1I/QSM2A
         SM2AR(N)=SM2A0R
         SM2AI(N)=SM2A0I
18       SM2BR(N)=+U(3)*FCT/D1(3)-SM2B0R
         SM2BI(N)=-V(3)*FCT/D1(3)-SM2B0I
         AM2BR(N)=-CNN*U(3)/D1(3)+SM2BR(N)
         AM2BI(N)=+CNN*V(3)/D1(3)+SM2BI(N)
         QSM2B=SM2BR(N)**2+SM2BI(N)**2
         SM2B0R=+SM2BR(N)/QSM2B
         SM2B0I=-SM2BI(N)/QSM2B
         AB(N)=(CNN+1.0_dp)/U(4)-SB0
         SB0=U(4)/(FCT-SB0*U(4))
15    SB(N)=SB0
      IF(Q.EQ.0.0_dp) GO TO 20
      QCM2A0=CN1R(2)**2+CN1I(2)**2
      BM2A0R=-SN1R(2)*CN1R(2)-SN1I(2)*CN1I(2)
      BM2A0R=BM2A0R/QCM2A0
      BM2A0I=+SN1R(2)*CN1I(2)-SN1I(2)*CN1R(2)
      BM2A0I=BM2A0I/QCM2A0
20    QSM2B0=SN1R(3)**2+SN1I(3)**2
      QCM2B0=CN1R(3)**2+CN1I(3)**2
      BM2B0R=-SN1R(3)*CN1R(3)-SN1I(3)*CN1I(3)
      BM2B0R=BM2B0R/QCM2B0
      BM2B0I=+SN1R(3)*CN1I(3)-SN1I(3)*CN1R(3)
      BM2B0I=BM2B0I/QCM2B0
      BB0=-SN1R(4)/CN1R(4)
      BDB0R=0.0_dp
      BDB0I=-1.0_dp
      IF(Q.EQ.0.0_dp) GO TO 22
      UM2R0=+SN1R(2)*SN1R(3)+SN1I(2)*SN1I(3)
      UM2R0=UM2R0/QSM2B0
      UM2I0=-SN1R(2)*SN1I(3)+SN1I(2)*SN1R(3)
      UM2I0=UM2I0/QSM2B0
      VM2R0=+CN1R(2)*CN1R(3)+CN1I(2)*CN1I(3)
      VM2R0=VM2R0/QCM2A0
      VM2I0=+CN1R(2)*CN1I(3)-CN1I(2)*CN1R(3)
      VM2I0=VM2I0/QCM2A0
22    CONTINUE
      DO 25 N=1,NX
         CN=N
         IF(Q.EQ.0.0_dp) GO TO 24
         CM2A0R=+CN*U(2)/D1(2)-BM2A0R
         CM2A0I=-CN*V(2)/D1(2)-BM2A0I
         QCM2A=CM2A0R**2+CM2A0I**2
         CM2AR(N)=+CM2A0R/QCM2A
         CM2AI(N)=-CM2A0I/QCM2A
         BM2AR(N)=-CN*U(2)/D1(2)+CM2AR(N)
         BM2AI(N)=+CN*V(2)/D1(2)+CM2AI(N)
         BM2A0R=BM2AR(N)
         BM2A0I=BM2AI(N)
24       CM2BR(N)=+CN*U(3)/D1(3)-BM2B0R
         CM2BI(N)=-CN*V(3)/D1(3)-BM2B0I
         CB(N)=CN/U(4)-BB0
         DBDB0R=+CN/U(4)-BDB0R
         DBDB0I=-BDB0I
         QCM2B=CM2BR(N)**2+CM2BI(N)**2
         QDB=DBDB0R**2+DBDB0I**2
         CM2B1R=+CM2BR(N)/QCM2B
         CM2B1I=-CM2BI(N)/QCM2B
         CB1=1.0_dp/CB(N)
         DBDB1R=DBDB0R/QDB
         DBDB1I=-DBDB0I/QDB
         BM2BR(N)=-CN*U(3)/D1(3)+CM2B1R
         BM2BI(N)=+CN*V(3)/D1(3)+CM2B1I
         BM2B0R=BM2BR(N)
         BM2B0I=BM2BI(N)
         BB(N)=-CN/U(4)+CB1
         BB0=BB(N)
         BDBR(N)=-CN/U(4)+DBDB1R
         BDBI(N)=DBDB1I
         BDB0R=BDBR(N)
         BDB0I=BDBI(N)
         SNB=SB(N)*SN1R(4)
         CNB=CB(N)*CN1R(4)
         SSB1R=AM2BR(N)-XM2*AB(N)
         SSB1I=AM2BI(N)-YM2*AB(N)
         SZB1R=AM2BR(N)-XM2*BDBR(N)+YM2*BDBI(N)
         SZB1I=AM2BI(N)-YM2*BDBR(N)-XM2*BDBI(N)
         SSB2R=XM2*AM2BR(N)-YM2*AM2BI(N)-AB(N)
         SSB2I=XM2*AM2BI(N)+YM2*AM2BR(N)
         SZB2R=XM2*AM2BR(N)-YM2*AM2BI(N)-BDBR(N)
         SZB2I=XM2*AM2BI(N)+YM2*AM2BR(N)-BDBI(N)
         IF(Q.NE.0.0_dp) GO TO 45
         ANNR=SNB*SSB1R
         ANNI=SNB*SSB1I
         ANDR=SNB*SZB1R-CNB*SZB1I
         ANDI=CNB*SZB1R+SNB*SZB1I
         BNNR=SNB*SSB2R
         BNNI=SNB*SSB2I
         BNDR=SNB*SZB2R-CNB*SZB2I
         BNDI=CNB*SZB2R+SNB*SZB2I
         GO TO 65
45       UM2R1=SM2AR(N)*SM2BR(N)-SM2AI(N)*SM2BI(N)
         UM2I1=SM2AR(N)*SM2BI(N)+SM2AI(N)*SM2BR(N)
         UM2R=UM2R1*UM2R0-UM2I1*UM2I0
         UM2I=UM2R1*UM2I0+UM2I1*UM2R0
         VM2R1=CM2AR(N)*CM2BR(N)-CM2AI(N)*CM2BI(N)
         VM2I1=CM2AR(N)*CM2BI(N)+CM2AI(N)*CM2BR(N)
         VM2R=VM2R1*VM2R0-VM2I1*VM2I0
         VM2I=VM2R1*VM2I0+VM2I1*VM2R0
         UVR=UM2R*VM2R-UM2I*VM2I
         UVI=UM2R*VM2I+UM2I*VM2R
         SSA1R=XM1*AM2AR(N)-YM1*AM2AI(N)- &
              XM2*AM1AR(N)+YM2*AM1AI(N)
         SSA1I=YM1*AM2AR(N)+XM1*AM2AI(N)- &
               YM2*AM1AR(N)-XM2*AM1AI(N)
         CSB1R=BM2BR(N)-XM2*AB(N)
         CSB1I=BM2BI(N)-YM2*AB(N)
         SCS1R=SSA1R*CSB1R-SSA1I*CSB1I
         SCS1I=SSA1R*CSB1I+SSA1I*CSB1R
         UVS1R=UVR*SCS1R-UVI*SCS1I
         UVS1I=UVR*SCS1I+UVI*SCS1R
         CSA1R=XM1*BM2AR(N)-YM1*BM2AI(N)- &
              XM2*AM1AR(N)+YM2*AM1AI(N)
         CSA1I=YM1*BM2AR(N)+XM1*BM2AI(N)- &
              YM2*AM1AR(N)-XM2*AM1AI(N)
         CSS1R=CSA1R*SSB1R-CSA1I*SSB1I
         CSS1I=CSA1R*SSB1I+CSA1I*SSB1R
         CZB1R=BM2BR(N)-XM2*BDBR(N)+YM2*BDBI(N)
         CZB1I=BM2BI(N)-YM2*BDBR(N)-XM2*BDBI(N)
         SCZ1R=SSA1R*CZB1R-SSA1I*CZB1I
         SCZ1I=SSA1R*CZB1I+SSA1I*CZB1R
         UVZ1R=UVR*SCZ1R-UVI*SCZ1I
         UVZ1I=UVI*SCZ1R+UVR*SCZ1I
         CSZ1R=CSA1R*SZB1R-CSA1I*SZB1I
         CSZ1I=CSA1R*SZB1I+CSA1I*SZB1R
         SSA2R=XM2*AM2AR(N)-YM2*AM2AI(N)- &
             XM1*AM1AR(N)+YM1*AM1AI(N)
         SSA2I=YM2*AM2AR(N)+XM2*AM2AI(N)- &
             YM1*AM1AR(N)-XM1*AM1AI(N)
         CSB2R=XM2*BM2BR(N)-YM2*BM2BI(N)-AB(N)
         CSB2I=XM2*BM2BI(N)+YM2*BM2BR(N)
         SCS2R=SSA2R*CSB2R-SSA2I*CSB2I
         SCS2I=SSA2R*CSB2I+SSA2I*CSB2R
         UVS2R=UVR*SCS2R-UVI*SCS2I
         UVS2I=UVR*SCS2I+UVI*SCS2R
         CSA2R=XM2*BM2AR(N)-YM2*BM2AI(N)- &
             XM1*AM1AR(N)+YM1*AM1AI(N)
         CSA2I=YM2*BM2AR(N)+XM2*BM2AI(N)- &
             YM1*AM1AR(N)-XM1*AM1AI(N)
         CSS2R=CSA2R*SSB2R-CSA2I*SSB2I
         CSS2I=CSA2R*SSB2I+CSA2I*SSB2R
         CZB2R=XM2*BM2BR(N)-YM2*BM2BI(N)-BDBR(N)
         CZB2I=XM2*BM2BI(N)+YM2*BM2BR(N)-BDBI(N)
         SCZ2R=SSA2R*CZB2R-SSA2I*CZB2I
         SCZ2I=SSA2R*CZB2I+SSA2I*CZB2R
         UVZ2R=UVR*SCZ2R-UVI*SCZ2I
         UVZ2I=UVI*SCZ2R+UVR*SCZ2I
         CSZ2R=CSA2R*SZB2R-CSA2I*SZB2I
         CSZ2I=CSA2R*SZB2I+CSA2I*SZB2R
55       ANNR=SNB*(UVS1R-CSS1R)
         ANNI=SNB*(UVS1I-CSS1I)
         ANDR=SNB*(UVZ1R-CSZ1R)-CNB*(UVZ1I-CSZ1I)
         ANDI=CNB*(UVZ1R-CSZ1R)+SNB*(UVZ1I-CSZ1I)
         BNNR=SNB*(UVS2R-CSS2R)
         BNNI=SNB*(UVS2I-CSS2I)
         BNDR=SNB*(UVZ2R-CSZ2R)-CNB*(UVZ2I-CSZ2I)
         BNDI=CNB*(UVZ2R-CSZ2R)+SNB*(UVZ2I-CSZ2I)
65       AND=ANDR*ANDR+ANDI*ANDI
         IF(AND.NE.0.0_dp) GO TO 70 
         ABANDR=DABS(ANDR)
         ABANDI=DABS(ANDI)
         ABANNR=DABS(ANNR)
         ABANNI=DABS(ANNI)
         AAAA=DMAX1(ABANDR,ABANDI,ABANNR,ABANNI)
         ANDR=ANDR/AAAA
         ANDI=ANDI/AAAA
         ANNR=ANNR/AAAA
         ANNI=ANNI/AAAA
         GO TO 65
70       CONTINUE
         AR(N)=(ANNR*ANDR+ANNI*ANDI)/AND
         AI(N)=(ANNI*ANDR-ANNR*ANDI)/AND
75       BND=BNDR*BNDR+BNDI*BNDI
         IF(BND.NE.0.0_dp) GO TO 80
         ABBNDR=DABS(BNDR)
         ABBNDI=DABS(BNDI)
         ABBNNR=DABS(BNNR)
         ABBNNI=DABS(BNNI)
         BBBB=DMAX1(ABBNDR,ABBNDI,ABBNNR,ABBNNI)
         BNDR=BNDR/BBBB
         BNDI=BNDI/BBBB
         BNNR=BNNR/BBBB
         BNNI=BNNI/BBBB
         GO TO 75
80       CONTINUE
         BR(N)=(BNNR*BNDR+BNNI*BNDI)/BND
         BI(N)=(BNNI*BNDR-BNNR*BNDI)/BND
         TI=AR(N)*AR(N)+AI(N)*AI(N)+BR(N)*BR(N)+ &
     	 BI(N)*BI(N)
         TI=TI/(AR(1)*AR(1)+AI(1)*AI(1)+BR(1)*BR(1)+ &
         BI(1)*BI(1))
         SN1R(4)=SNB
         CN1R(4)=CNB
         UM2R0=UM2R
         UM2I0=UM2I
         VM2R0=VM2R
         VM2I0=VM2I
         IF(NSTOP-N) 50,50,25
!C         IF(TI-CTST) 50,50,25
25    CONTINUE
50    CONTINUE
      do i=1,NSTOP
         an(i)=dcmplx(AR(i),-AI(i))
         bn(i)=dcmplx(BR(i),-BI(i))
      enddo
             
                        
    end subroutine
    
    subroutine cofsrd(nmax)
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: cofsrd
          include 'tmatrix_mps_np.inc'
          integer(kind=int4) :: nmax
          ! Locals
          integer(kind=int4), parameter :: nmp = np*(np+2)
          real(kind=dp), dimension(nmp) :: cofsr
          real(kind=dp) :: c
          integer(kind=int4) :: i,n,m
          common/crot/cofsr
          ! Exec code ....
          
          i=0
          do n=1,nmax
             do m=-n,n
                i=i+1       
                c=lnf(dble(n-m))-lnf(dble(n+m))
                cofsr(i)=0.5_dp*c
!c               c=0.5d0*c
!c            cofsr(i)=dexp(c)
         enddo
      enddo
    end subroutine
    
    subroutine cofd0(nmax)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: cofd0
          include 'tmatrix_mps_np.inc'
          integer(kind=int4) :: nmax
          ! Locals
          integer(kind=int4), parameter :: nmp = np*(np+2)
          integer(kind=int4), parameter :: ni0 = np*(np+1)*(2*np+1)/3+np*np
        
          integer(kind=int4) :: v,i,m,ns,inm,ivm,n	
          real(kind=dp) :: sm,c,c0,c1
          real(kind=dp), dimension(ni0)        :: cof0
          real(kind=dp), dimension(nmp)        :: cofsr
          real(kind=dp), dimension(0:4*(np+1)) :: fnr
          common/cofmnv0/cof0
          common/crot/cofsr	
          common/fnr/fnr
          ! Exec code ...
          
          i=0
          sm=-0.5_dp*dble((-1)**nmax)
          do m=-nmax,nmax
             ns=max(1,iabs(m))
             sm=-sm
             do n=ns,nmax
                inm=n*(n+1)-m
                do v=ns,nmax
                   i=i+1
                   ivm=v*(v+1)+m
                   c=cofsr(inm)+cofsr(ivm)	            
                   c=sm*dexp(c)
                   c0=fnr(2*n+1)*fnr(2*v+1)
                   c1=fnr(n)*fnr(v)*fnr(n+1)*fnr(v+1)
                   c0=c0/c1
                   cof0(i)=c*c0
               enddo
            enddo
          enddo
    end subroutine
    
    subroutine cofnv0(nmax)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: cofnv0
          include 'tmatrix_mps_np.inc'
          integer(kind=int4) :: nmax
          ! Locals
          integer(kind=int4) :: n,v
          real(kind=dp) :: c1
          real(kind=dp), dimension(np,np) :: cnv
          common/cfnv/cnv
          ! Exec code ...
          do n=1,nmax
             do v=n,nmax
                c1=lnf(dble(2*n))+lnf(dble(2*v))
                c1=c1-lnf(dble(2*n+2*v))
                c1=c1+2.0_dp*lnf(dble(n+v))
                c1=c1-lnf(dble(n))-lnf(dble(v))
                cnv(n,v)=c1
            enddo
          enddo
    end subroutine
    
!    C  subroutine gau0.f generates tabulated values for  
!C  Gaunt coefficients up to n=v=n_max
      subroutine gau0(nmax)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gau0
          include 'tmatrix_mps_np.inc'
          integer(kind=int4) :: nmax
          integer(kind=int4), parameter :: ni0 =   np*(np+1)*(2*np+1)/3+np*np
          integer(kind=int4), parameter :: ng0 =   np*(2*np**3+10*np**2+19*np+5)/6
          integer(kind=int4) :: v,qmax,uvmax,i,na,m,ns,n
          integer(kind=int4), dimension(ni0) :: iga0
          real(kind=dp), dimension(ng0) :: ga0
          common/g0/ga0
          common/ig0/iga0
          ! Exec code ....
          
          na=0
          uvmax=nmax*(nmax+2)
          i=0
          do m=-nmax,nmax
             ns=max(1,iabs(m))
             do n=ns,nmax
                do v=ns,nmax
                   call gxurcd0(-m,n,v,qmax,na)
                   i=i+1
                   iga0(i)=na
                   na=na+qmax+1
               enddo
           enddo
          enddo
      end subroutine
      
!      c  transforms the rectangular coordinates (x,y,z)
!c  to spherical coordinates (r,theta,phi)
      subroutine carsphd(x,y,z,r,xt,sphi,cphi)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: carsphd
          real(kind=dp) :: x,y,z,r,xt,sphi,cphi
          r=dsqrt(x*x+y*y+z*z)
          if(r.eq.0._dp) then
             xt=1._dp
             sphi=0._dp
             cphi=1._dp
             return
          endif
          xt=z/r
          if(y.eq.0._dp.and.x.eq.0._dp) then
             sphi=0._dp
             cphi=1._dp
             return
          endif
      sphi=dsqrt(x*x+y*y)
      cphi=x/sphi
      sphi=y/sphi
     
      end  subroutine
      
!      c  subroutine besseljd.f  (in double precision arithmetic)
!c  returns an array of the spherical Bessel function of the  
!c  first kind with a real argument: j_0,j_1,j_2,...,j_{NC}
!c  uses Ru Wang's ratio method for the downward recursive 
!c  calculation of the Riccati-Bessel function Psi_n(z)=z j_n(z) 
!c  [see Xu et al., Physical Review E, v.60, 2347-2365 (1999)]
    SUBROUTINE besseljd(NC,X,BESJ)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: besseljd
          INTEGER(kind=int4) ::  NC,NX,K,N
          real(KIND=DP) :: X,PN,CN,X2
          real(KIND=DP), dimension(0:NC) ::  BESJ
          DO K=1,NC
             BESJ(K)=0._dp
          ENDDO
          IF(DABS(X).LT.1.0E-12_dp) THEN
             BESJ(0)=1.0_dp
             BESJ(1)=0.3333333333333333333_dp*X
             RETURN
          ENDIF
!c  down-recursively calculating an array of the ratio functions  
!c  P_{NC},P_{NC-1},...,P(2) stored in the same array for j_n, 
!c  starting with an asymptotic value P_{NX+1}=X/(2 NX+3) at the 
!c  highest order NX+1, where NX=NC+1.1X+10 
          NX=1.1_dp*X+10._dp
          NX=NC+NX   
          PN=X/DBLE(2*NX+3)
          DO 5 K=1,NX-1
               N=NX-K+1
               CN=DBLE(N)
               PN=X/(DBLE(2*N+1)-PN*X)
               IF(N.GT.NC) GOTO 5
                  BESJ(N)=PN
       5 CONTINUE
!C  calculating j_0(x) and j_1(x)
         IF(DABS(X)-0.1_dp) 10,11,11
    10   X2=X*X
         BESJ(0)=1._dp-X2/72._dp
         BESJ(0)=1._dp-X2/42._dp*BESJ(0)
         BESJ(0)=1._dp-X2/20._dp*BESJ(0)
         BESJ(0)=1._dp-X2/6._dp*BESJ(0)
         BESJ(1)=1._dp/45360._dp-X2/3991680._dp
         BESJ(1)=1._dp/840._dp-X2*BESJ(1)
         BESJ(1)=1._dp/30._dp-X2*BESJ(1)
         BESJ(1)=X*(0.3333333333333333_dp-X2*BESJ(1))
      GOTO 12
   11    BESJ(0)=DSIN(X)/X
         BESJ(1)=(DSIN(X)/X-DCOS(X))/X      
!c  calculating j_2,...,j_{NC} 
   12 DO 20 N=2,NC
         BESJ(N)=BESJ(N)*BESJ(N-1)
   20 CONTINUE
    end subroutine
    
!    c  sub. besselyd.f  (in double precision arithmetic)
!c  returns an array of the spherical Bessel function of
!c  the second kind with a real argument: y_0,y_1,...,y_n

    subroutine besselyd(n,x,besy)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: besselyd
          integer(kind=int4) ::  i,n
          real(kind=dp) :: x,besyn,x2
          real(kind=dp), dimension(0:n) :: besy
          ! Exec code ...
          if(x.eq.0.0_dp) then
             write(6,'(A32 F22.15)') 'bad argument in sub. besselyd', x
             call TRACEBACKQQ(STRING="bad argument in sub. besselyd",USER_EXIT_CODE = -1)
             stop
          endif
         if(dabs(x)-0.1_dp)10,11,11
  10        x2=x*x
            besyn=1._dp-x2/72._dp
            besyn=1._dp-x2/42._dp*besyn
            besyn=1._dp-x2/20._dp*besyn
            besyn=1._dp-x2/6._dp*besyn
            besy(0)=1._dp-x2/56._dp
            besy(0)=1._dp-x2/30._dp*besy(0)
            besy(0)=1._dp-x2/12._dp*besy(0)
            besy(0)=1._dp-0.5_dp*x2*besy(0)
            besy(0)=-besy(0)/x
            goto 12
  11        besyn=dsin(x)/x
            besy(0)=-dcos(x)/x
  12        besy(1)=besy(0)/x-besyn
            do i=2,n
               besy(i)=dble(2*i-1)/x*besy(i-1)-besy(i-2)
            enddo
     
    end  subroutine
    
 ! c  subroutine rotcoef is originally written by Mackowski, Fuller, and 
!c  Mishchenko (taken from the code scsmfo1b.for released to public by 
!c  the authors at 8/1999)
!c  the rotational coefficients are required for the implementation of 
!c  Mackowski's three-step numerical technique in subroutine rtr.f for 
!c  decomposition of translation matrix into rotational and axial 
!c  translational parts [see Mackowski, Proc. R. Soc. Lond. A 433, 599 
!c  (1991)]
!c  Yu-lin Xu   12/2000   
    subroutine rotcoef(cbe,nmax)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: rotcoef
          include 'tmatrix_mps_np.inc'
          real(kind=dp) :: cbe
          integer(kind=int4) :: nmax
          ! Locals
          integer(kind=int4), parameter ::  nmp =  np*(np+2)
          real(kind=dp) :: sbe,cbe2,sben,dkt,fmn,dkm0,dkm1, &
                        dkn1,sbe2
          integer(kind=int4) :: inn,n,nn1,k,im1,m,m1,kn,IM
          real(kind=dp), dimension(-2*np:2*np) :: dk0,dk01
!DIR$     ATTRIBUTES ALIGN : 64 :: dk0,dk01    
          real(kind=dp), dimension(0:np+2) :: bcof
          real(kind=dp), dimension(-np:np,0:nmp) :: dc
          real(kind=dp), dimension(0:4*(np+1)) :: fnr
          common/rot/bcof,dc
          common/fnr/fnr
          ! Exec code ....
         ! dk0  = 0._R64P
          ! dk01 = 0._R64P
          
          sbe=dsqrt((1._dp+cbe)*(1._dp-cbe))
          cbe2=.5_dp*(1._dp+cbe)
          sbe2=.5_dp*(1._dp-cbe)
          inn=1           ! was 'in'
          dk0(0)=1._dp
          sben=1._dp
          dc(0,0)=1._dp
          dk01(0)=0._dp
          do n=1,nmax
             nn1=n*(n+1)
             inn=-inn
             sben=sben*sbe/2._dp
             dk0(n)=dble(inn)*sben*bcof(n)
             dk0(-n)=dble(inn)*dk0(n)
             dk01(n)=0._dp
             dk01(-n)=0._dp
             dc(0,nn1+n)=dk0(n)
             dc(0,nn1-n)=dk0(-n)
             do k=-n+1,n-1
                kn=nn1+k
                dkt=dk01(k)
                dk01(k)=dk0(k)
                dk0(k)=(cbe*dble(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt)/ &
                       (fnr(n+k)*fnr(n-k))
                dc(0,kn)=dk0(k)
            enddo
            im=1
            do m=1,n
               im=-im
               fmn=1._dp/fnr(n-m+1)/fnr(n+m)
               m1=m-1
               dkm0=0._dp
               do k=-n,n
                  kn=nn1+k
                  dkm1=dkm0
                  dkm0=dc(m1,kn)
                  if(k.eq.n) then
                     dkn1=0._dp
                  else
                     dkn1=dc(m1,kn+1)
                  endif
                  dc(m,kn)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1- &
                           fnr(n-k)*fnr(n+k+1)*sbe2*dkn1-  &
                           dble(k)*sbe*dc(m1,kn))*fmn
                  dc(-m,nn1-k)=dble((-1)**(k)*im)*dc(m,kn)
               enddo
            enddo
          enddo
    end subroutine
    
!    c  subroutine cofxuds0.f returns the two classes of vector 
!c  (axial) translation coefficients for a given combination of 
!c  (m,n,m,v) and a given dimensionless translation distance kd 
!cu uses subroutine gid0.f 
    subroutine cofxuds0(nmax,m,n,v,sja,sya,A,B,Aj,Bj)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: cofxuds0
      implicit double precision (a-h,o-z)
      include 'tmatrix_mps_np.inc'
      integer(kind=int4) :: nmax,m,n
      integer(kind=int4), parameter :: ni0 = np*(np+1)*(2*np+1)/3+np*np
      integer(kind=int4), parameter :: ng0 = np*(2*np**3+10*np**2+19*np+5)/6
      integer(kind=int4) :: v,p,qmax,ig,i,nv2,id
      real(kind=dp), dimension(0:n+v+1) ::  sja,sya
      integer(kind=int4), dimension(ni0) :: iga0
      complex(16) ::  A,B,Aj,Bj,signz
      common/ig0/iga0
      common/g0/ga0(ng0)
      common/cofmnv0/cof0(ni0)
     
      
      fa(m,p)=dble(-2*m*p*(p-1))
      fb(n,v,p)=dble(p*p-(n+v+1)*(n+v+1))* &
                dble(p*p-(n-v)*(n-v))/dble(4*p*p-1)
      if(iabs(m).gt.n.or.iabs(m).gt.v) then 
         write(6,*) '|m|>n or v in subroutine cofxuds0.f'
         stop
      endif
      A=0.0_dp
      B=0.0_dp
      Aj=0.0_dp
      Bj=0.0_dp
      call gid0(nmax,m,n,v,id)
      c=cof0(id)
      ig=iga0(id)
      nv2=v*(v+1)+n*(n+1)
      signz=dcmplx(0.0_dp,1.0_dp)**(n+v)
      p=n+v+2
      qmax=min(n,v)
      do i=1,qmax+1
         p=p-2
         cp=dble(nv2-p*(p+1))*ga0(ig+i)
         sj=sja(p)
         sy=sya(p)
         A=A+dcmplx(sj,sy)*signz*cp
         Aj=Aj+sj*signz*cp
         signz=-signz
      enddo
      A=A*c
      Aj=Aj*c
      if(m.eq.0) return
      signz=dcmplx(0.0_dp,1.0_dp)**(n+v+1)
      p=n+v
      do 20 i=1,qmax
         p=p-2
         signz=-signz
         if(i.eq.1) then
            cp=dble(2*p+3)*fa(m,p+3)      
            cp=cp*ga0(ig+1)/dble((p+3)*(p+2))
            goto 21
         endif
         if(i.eq.qmax) then
            if(p.eq.0) goto 22
            nv2=p*(p+1)
            cp=dble(2*p+3)*fa(m,p+1)
            cp=-cp*ga0(ig+i+1)/dble(nv2)
            goto 21
         endif
 22      c4=fa(m,p+2)
         cp=-dble((p+1)*(p+2))*fb(n,v,p+2)*ga0(ig+i)
         cp=cp+dble((p+2)*(p+1))*fb(n,v,p+1)*ga0(ig+i+1)
         cp=cp*dble(2*p+3)/c4
 21      sj=sja(p+1)
         sy=sya(p+1)
         B=B+dcmplx(sj,sy)*signz*cp
         Bj=Bj+sj*signz*cp
 20   continue
      B=B*c
      Bj=Bj*c
      
    end  subroutine
    
!    c  subroutine tipitaud.f 
!c  calculates pi(m,n) & tau(m,n) up to a specified nmax for all 
!c  m=0,1,...n at a given x
!c  pi(m,n) & tau(m,n) calculated are normalized by  
!c         C_mn=[(2n+1)(n-m)!/n/(n+1)/(n+m)!]^(1/2)
!c  Yu-lin Xu    12/2000

    subroutine tipitaud(nmax,x)
      ! include 'gmm01f.par'
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: tipitaud
     
      implicit double precision (a-h,o-z)
      include 'tmatrix_mps_np.inc'
      integer(kind=int4), parameter :: nmp0 = (np+1)*(np+4)/2
      integer(kind=int4) :: nmax
      integer(kind=int4) :: nt,imn,m,n,i,j,i1,i2
     
      common/fnr/fnr(0:4*(np+1))
      common/pitau/pi(nmp0),tau(nmp0)
      
      nt=(nmax+1)*(nmax+4)/2         ! calculates pi up to nmax+1
      if(nt.gt.nmp0.or.dabs(x).gt.1.d0) then
         write(6,*) 'dimension or argument wrong in sub. tipitaud'
         write(6,*) 'argument: ',x
         stop
      endif
      sx=dsqrt(1.d0-x*x)
      pi(1)=0.d0                     ! pi(0,1)  pi(0,n)=0 when m=0
      pi(2)=dsqrt(.75d0)             ! pi(1,1)
      pi(3)=0.d0                     ! pi(0,2)
      t125=dsqrt(1.25d0) 
      pi(4)=t125*x                      ! pi(1,2)  
      pi(5)=t125*sx                     ! pi(2,2)
      imn=5
      do i=3,nmax+1
         n=i
         imn=imn+1
         pi(imn)=0.d0                ! pi(0,n)=0
         do 11 j=2,n
            m=j-1
            imn=imn+1
            i1=(n-2)*(n+1)/2+m+1
            if(m.eq.n-1) then
               pi(imn)=fnr(n-1)*fnr(2*n+1)/fnr(n+1)*x*pi(i1)
               goto 11
            endif
            i2=(n-3)*n/2+m+1
            t=fnr(n)*fnr(2*n-3)
            t=fnr(n-2)*fnr(n-m-1)*fnr(n+m-1)/t
            pi(imn)=fnr(2*n-1)*x*pi(i1)-t*pi(i2)
            t=fnr(n+1)*fnr(n-m)*fnr(n+m)
            t=fnr(n-1)*fnr(2*n+1)/t
            pi(imn)=t*pi(imn)
 11      continue
         imn=imn+1
         i1=(n-2)*(n+1)/2+n
         t=fnr(n-1)*fnr(n+1)
         t=dsqrt(.5d0)*fnr(n)*fnr(2*n+1)/t
         pi(imn)=t*sx*pi(i1)
      enddo
      tx=x*sx
      tau(1)=-dsqrt(1.5d0)*sx          ! tau(0,1)
      tau(2)=pi(2)*x                   ! tau(1,1)
      tau(3)=-dsqrt(7.5d0)*tx          ! tau(0,2)   
      tau(4)=t125*(2.d0*x*x-1.d0)      ! tau(1,2)
      tau(5)=t125*tx                   ! tau(2,2)
      imn=5
      do i=3,nmax
         n=i
         do 21 j=1,n+1
            m=j-1
            imn=imn+1
            if(m.eq.0) then
               i1=(n-2)*(n+1)/2+1
               i2=(n-3)*n/2+1
               t=fnr(2*n-3)
               t=fnr(n-2)*fnr(n)/t
               tau(imn)=fnr(2*n-1)*x*tau(i1)-t*tau(i2)
               t=fnr(n-1)*fnr(n+1)
               t=fnr(2*n+1)/t
               tau(imn)=t*tau(imn)
               goto 21
            endif
            i1=n*(n+3)/2+m+1
            t=fnr(n)*fnr(2*n+3)
            t=fnr(n+2)*fnr(2*n+1)*fnr(n-m+1)*fnr(n+m+1)/t
            tau(imn)=t*pi(i1)-dble(n+1)*x*pi(imn)
            tau(imn)=tau(imn)/dble(m)
 21      continue
      enddo
   
    end subroutine
    
!    c  lnfacd.f  (double precision function)
!c  returns ln(z!)  z>-1.0
!c  based on Lanczos' method [see Xu, Journal of Computational 
!c  Physics, v.139, 137-165 (1998)]
    
     ! integer(kind=int4) :: i
     ! double precision z,a,b,cp,c0(11)
     ! data c0/0.16427423239836267d5, -0.48589401600331902d5, &
     !         0.55557391003815523d5, -0.30964901015912058d5, &
     !         0.87287202992571788d4, -0.11714474574532352d4, &
     !         0.63103078123601037d2, -0.93060589791758878d0, &
     !         0.13919002438227877d-2,-0.45006835613027859d-8,&
     !         0.13069587914063262d-9/ 
     ! a=1.d0
     ! cp=2.5066282746310005d0
    !  b=z+10.5d0
    !  b=(z+0.5d0)*dlog(b)-b
    !  do i=1,11
    !    z=z+1.d0
    !    a=a+c0(i)/z
    !  enddo
    !  lnfacd=b+dlog(cp*a)
    !  
   ! end  function lnfacd

  
      
!      c  gxurcd0.f to compute Gaunt coefficients a(-m,n,m,v,p)
!cu uses lnfacd.f to compute ln(z!)
      subroutine gxurcd0(m,n,v,qmax,na)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gxurcd0
        implicit double precision (a-h,o-z)
        include 'tmatrix_mps_np.inc'
        integer(kind=int4), parameter ::  ng0 = np*(2*np**3+10*np**2+19*np+5)/6
       
        integer(kind=int4) :: m, v,qmax,p,n,na,nq,  &
                              i
      double precision cnv(np,np),ga0(ng0)
      common/cfnv/cnv
      common/g0/ga0
      ! Exec code .....
     
      fb(n,v,p)=dble(p-(n+v+1))*dble(p+(n+v+1))* &
                dble(p-(n-v))*dble(p+(n-v))/ &
                (dble(2*p+1)*dble(2*p-1))	
      if(iabs(m).gt.n.or.iabs(m).gt.v) then
         write(6,*) 'warning: |m|>n or v in gxurcd0'
         qmax=-1	   
         return
      endif
      qmax=min(n,v)
      nq=qmax+1
      if(n.le.v) then 
         c1=cnv(n,v)
      else
         c1=cnv(v,n)
      endif
      c1=c1-lnf(dble(n-m))-lnf(dble(v+m))
      ga0(na+1)=dexp(c1)
      if(qmax.lt.1) return	
      p=n+v
      do 8 i=2,nq
         p=p-2
         if(m.eq.0) then
            c1=fb(n,v,p+1)
            c2=fb(n,v,p+2)
            goto 2
         endif
         c1=fb(n,v,p+1)
         c2=dble(4*m*m)+fb(n,v,p+2)+fb(n,v,p+3)
         if(i.eq.2) goto 2
         c3=-fb(n,v,p+4)
         goto 4
  2      ga0(na+i)=c2*ga0(na+i-1)/c1
         goto 8
  4      ga0(na+i)=(c2*ga0(na+i-1)+c3*ga0(na+i-2))/c1
  8   continue
     
      end subroutine
      
     subroutine gid0(nmax,m,n,iv,id)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gid0
           integer(kind=int4) :: nmax,m,n,iv,id
           ! Locals
           integer(kind=int4) :: nt,ns,nc0
           nt=nmax*(nmax+1)*(2*nmax+1)/3+nmax*nmax
           ns=max(1,iabs(m))
           nc0=nmax-iabs(m)
           id=nc0*(nc0+1)*(2*nc0+1)/6
           if(m) 10,11,12
 10        id=id+(n-ns)*(nc0+1)+iv-ns+1  
           return
 11        id=id+(n-ns)*nmax+iv
           return
 12        id=id+(nmax-n)*(nc0+1)+nmax-iv
           id=nt-id
           
     end  subroutine
     
     function plgndrd(l,mr,x)  
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: plgndrd
      integer(kind=int4) ::  l,mr,m,index,i,ll
      real(kind=dp) :: x
     
      real(kind=dp) :: fact,pll,pmm,pmmp1,somx2
      real(kind=dp) ::  plgndrd
      m=mr
      index=0
      if(m.lt.0) then
         index=1
         m=-m
      endif
      if(dabs(x).gt.1.0_dp) then
	     write(6,*) 'n,m,x: ', l,(-2*index+1)*m,x
         stop 'bad arguments in plgndrd'
      end if
      if(m.gt.l) then
	     plgndrd=0._dp
         return
      end if
      pmm=1.0_dp
      if(m.gt.0) then
        somx2=dsqrt((1._dp-x)*(1._dp+x))
        fact=1._dp
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2._dp
11      continue
      endif
      if(l.eq.m) then
        plgndrd=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
           plgndrd=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndrd=pll
        endif
      endif
      plgndrd=-plgndrd
      if(m/2*2.eq.m) plgndrd=-plgndrd
      if(index.gt.0) then
         fact=lnf(dble(l-m))-lnf(dble(l+m))
         fact=dexp(fact)
          plgndrd=-plgndrd*fact
         if(m/2*2.eq.m) plgndrd=-plgndrd
      endif 
     
     end  function
     
    subroutine rtrT(anpt,nodrj,nodri,ekt,drot,ij1,ij2,ii1,ii2)      
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: rtrT

          include 'tmatrix_mps_np.inc'
          integer(kind=int4), parameter ::  nmp = np*(np+2), &
                                            nrc = 4*np*(np+1)*(np+2)/3+np, &
                                            nij =  nLp*(nLp-1)/2
          integer(kind=int4) :: nodrj,nodri
          
          
          complex(16), dimension(2,nmp) :: anpt
          complex(16), dimension(np)    :: ekt
          real(kind=dp),  dimension(nrc)   :: drot
          integer(kind=int4) :: ij1,ij2,ii1,ii2
          ! Locals
          complex(16), dimension(2,2*np) :: ant
!DIR$     ATTRIBUTES ALIGN : 64 :: ant
          complex(16), dimension(2,-np:np) :: amt
!DIR$     ATTRIBUTES ALIGN : 64 :: amt
          complex(16), dimension(-np:np) :: ek
!DIR$     ATTRIBUTES ALIGN : 64 :: ek
          complex(16), dimension(2,np,nmp) :: atr
          complex(16) :: a,b
          integer(kind=int4) :: nj1,nj2,ni1,ni2,nmax,irc,n,k,kn,mmax, &
                           l,ml,inn,imn,m,n1,n1j,ip,n1i
          real(kind=dp) :: sik
          common/tran/atr
          ! Exec code ...
         
        
         
          nj1=dsqrt(dble(ij1))
          nj2=dsqrt(dble(ij2))
          ni1=dsqrt(dble(ii1))
          ni2=dsqrt(dble(ii2))
          ek(0)=1._dp
          nmax=max(nj2,ni2)
         do m=1,nmax
            ek(m)=ekt(m)
            ek(-m)=conjg(ek(m))
        enddo
        irc=0
        do n=1,nj1-1
           do k=-n,n
              do m=-n,n
                 irc=irc+1
              enddo
         enddo
       enddo
       do n=nj1,nj2
           n1=n*(n+1)
            do m=-n,n
               amt(1,m)=0._dp
               amt(2,m)=0._dp
           enddo
           do 10 k=-n,n
               kn=n1+k
               if(kn.lt.ij1) then
                  irc=irc+2*n+1
                 goto 10
              endif
            if(kn.gt.ij2) goto 11
            a=ek(k)*anpt(1,kn)
            b=ek(k)*anpt(2,kn)
            do m=-n,n
               irc=irc+1
               amt(1,m)=amt(1,m)+a*drot(irc)
               amt(2,m)=amt(2,m)+b*drot(irc)
            enddo
 10      continue
 11      do m=-n,n
            imn=n1+m
            anpt(1,imn)=amt(1,m)
            anpt(2,imn)=amt(2,m)
         enddo
      enddo
      mmax=min(nj2,ni2)
      do m=-mmax,mmax
         n1=max(1,iabs(m))
         n1j=max(n1,nj1)
         do n=n1j,nj2
            imn=n*(n+1)+m
            do ip=1,2
               ant(ip,n)=anpt(ip,imn)
            enddo
         enddo
         n1i=max(n1,ni1)
         do n=n1i,ni2
            imn=n*(n+1)+m
            a=0.0_dp
            b=0.0-dp
            do 20 l=n1j,nj2
               ml=l*(l+1)+m
               a=a+atr(1,n,ml)*ant(1,l)+ &
                 atr(2,n,ml)*ant(2,l)
               b=b+atr(1,n,ml)*ant(2,l)+ &
                 atr(2,n,ml)*ant(1,l)
 20         continue
            anpt(1,imn) = a
            anpt(2,imn) = b
         enddo
      enddo
      inn=1
      irc=0
      do n=1,ni1-1
         do k=-n,n
            do m=-n,n
               irc=irc+1
            enddo
         enddo
      enddo
      do n=ni1,ni2
         inn=-inn
         n1=n*(n+1)
         do m=-n,n
            amt(1,m)=0.0_dp
            amt(2,m)=0.0_dp
         enddo
         sik=-inn
         do 30 k=-n,n
	    sik=-sik
	    kn=n1+k
	    a=sik*anpt(1,kn)
            b=sik*anpt(2,kn)
            do m=-n,n
               irc=irc+1
               amt(1,m)=amt(1,m)+a*drot(irc)
               amt(2,m)=amt(2,m)+b*drot(irc)
            enddo
 30      continue
         sik=-inn
         do 31 m=-n,n
            sik=-sik
            imn=n1+m
            if(imn.lt.ii1.or.imn.gt.ii2) goto 31
            anpt(1,imn)=amt(1,m)*ek(-m)*sik
            anpt(2,imn)=amt(2,m)*ek(-m)*sik
 31      continue
      enddo
     
    end subroutine
    
    subroutine transT(nL,r0,nmax,uvmax,fint,atr0,btr0,ek,      &
                      drot,as,bs,as1,bs1,ind,confg,iuvc,isw)
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: transT
          include 'tmatrix_mps_np.inc'
          integer(kind=int4), parameter ::  nmp = np*(np+2), &
                                            ni0 = np*(np+1)*(2*np+1)/3+np*np, &
                                            nrc = 4*np*(np+1)*(np+2)/3+np, &
                                            nij = nLp*(nLp-1)/2
          integer(kind=int4) :: nL
          real(kind=dp), dimension(6,nLp) :: r0
          integer(kind=int4), dimension(nLp) :: nmax,uvmax,ind 
          real(kind=dp) :: fint
          
                                     
          complex(16), dimension(2,np,nmp) :: atr
          complex(16), dimension(ni0,nij)  :: atr0,btr0
          complex(16), dimension(np,nij)   :: ek
          complex(16), dimension(nLp,nmp)  :: as,bs,as1,bs1
          
          real(kind=dp),  dimension(nrc,nij)  :: drot
          real(kind=dp),  dimension(5,nij)    :: confg
          integer(kind=int4) :: iuvc,isw
          ! Locals
          complex(16), dimension(2,nmp) :: at1
          integer(kind=int4) :: i,imn,j,ij,nlarge,itrc,nsmall,m,n1,n, &
                           v,ij1,ij2,ii1,ii2,iuv
          real(kind=dp) :: x0,y0,z0,temp,sic
          common/tran/atr
          ! Exec code ....
         
          do i=1,nL
             do imn=1,uvmax(i)
                 as1(i,imn)=dcmplx(0._dp,0._dp)
                 bs1(i,imn)=dcmplx(0._dp,0._dp)
             enddo
          enddo
          do 11 i=1,nL-1
             if(ind(i).gt.0) goto 11
             do 10 j=i+1,nL
             if(ind(j).gt.0) goto 10
                ij=(j-1)*(j-2)/2+j-i
                x0=confg(1,ij)
                y0=confg(2,ij)
                z0=confg(3,ij)
                temp=confg(5,ij)
                if(temp.le.fint) goto 10
                nlarge=max(nmax(i),nmax(j))
                itrc=0
                nsmall=min(nmax(i),nmax(j))
            do m=-nsmall,nsmall
               n1=max(1,iabs(m))
               do n=n1,nlarge
                  do v=n1,nlarge
                     itrc=itrc+1
                     iuv=v*(v+1)+m
                     atr(1,n,iuv)=atr0(itrc,ij)
                     atr(2,n,iuv)=btr0(itrc,ij)
                     if(x0.eq.0._dp.and.y0.eq.0._dp) then
                        if(z0.lt.0._dp) goto 20
                     endif
                     goto 21
 20                  sic=dble((-1)**(n+v))
                     atr(1,n,iuv)=sic*atr(1,n,iuv)
                     atr(2,n,iuv)=-sic*atr(2,n,iuv)
 21                  continue
                  enddo
               enddo
            enddo
            do iuv=1,uvmax(j)
               at1(1,iuv)=as(j,iuv)
               at1(2,iuv)=bs(j,iuv)
            enddo
            If(nmax(j).lt.nlarge) then
               do n=nmax(j)+1,nlarge
                  do m=-n,n
                     iuv=n*(n+1)+m
                     at1(1,iuv)=dcmplx(0._dp,0._dp)
                     at1(2,iuv)=dcmplx(0._dp,0._dp)
                  enddo
               enddo
            endif
            goto(22,23),isw
 22         ij1=iuvc
            ij2=uvmax(j)
            ii1=iuvc
            ii2=uvmax(i)
            goto 24
 23         ij1=1
            ij2=iuvc-1
            ii1=iuvc
            ii2=uvmax(i)
 24         if(x0.eq.0._dp.and.y0.eq.0._dp) then
               call trvT(at1,nmax(j),nmax(i),ij1,ij2,ii1,ii2)
            else
               call rtrT(at1,nmax(j),nmax(i),ek(1,ij),  &
                        drot(1,ij),ij1,ij2,ii1,ii2)
            endif
            do imn=ii1,ii2
               as1(i,imn)=as1(i,imn)+at1(1,imn)
               bs1(i,imn)=bs1(i,imn)+at1(2,imn)
            enddo
            do m=-nsmall,nsmall
               n1=max(1,iabs(m))
               do n=n1,nlarge
                  do v=n1,nlarge
                     iuv=v*(v+1)+m
                     sic=dble((-1)**(n+v))
                     atr(1,n,iuv)=sic*atr(1,n,iuv)
                     atr(2,n,iuv)=-sic*atr(2,n,iuv)
                  enddo
               enddo
            enddo
            do iuv=1,uvmax(i)
               at1(1,iuv)=as(i,iuv)
               at1(2,iuv)=bs(i,iuv)
            enddo
            If(nmax(i).lt.nlarge) then
               do n=nmax(i)+1,nlarge
                  do m=-n,n
                     iuv=n*(n+1)+m
                     at1(1,iuv)=dcmplx(0._dp,0._dp)
                     at1(2,iuv)=dcmplx(0._dp,0._dp)
                  enddo
               enddo
            endif
            goto(25,26),isw
 25         ii1=iuvc
            ii2=uvmax(i)
            ij1=iuvc
            ij2=uvmax(j)
            goto 27
 26         ii1=1
            ii2=iuvc-1
            ij1=iuvc
            ij2=uvmax(j)
 27         if(x0.eq.0._dp.and.y0.eq.0._dp) then
               call trvT(at1,nmax(i),nmax(j),ii1,ii2,ij1,ij2)
            else 
               call rtrT(at1,nmax(i),nmax(j),ek(1,ij),   &
                        drot(1,ij),ii1,ii2,ij1,ij2)
            endif
            do imn=ij1,ij2
               as1(j,imn)=as1(j,imn)+at1(1,imn)
               bs1(j,imn)=bs1(j,imn)+at1(2,imn)
            enddo
 10      continue
 11   continue
     
    end  subroutine
                     
    subroutine trvT(anpt,nodrj,nodri,ij1,ij2,ii1,ii2)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: trvT
          include 'tmatrix_mps_np.inc'
          integer(kind=int4), parameter :: nmp = np*(np+2)
          complex(16), dimension(2,nmp) :: anpt
          integer(kind=int4) :: nodrj,nodri,ij1,ij2,ii1,ii2
          ! Locals
          complex(16), dimension(2,2*np) :: ant
          complex(16) :: a,b
          complex(16), dimension(2,np,nmp) :: atr
          integer(kind=int4) :: nji,nj2,ni1,ni2,mmax,m,n1,nj1,imn,  &
                           ip,n1i,l,ml,n,n1j
          common/tran/atr
          ! Exec code ...
       
          nj1=dsqrt(dble(ij1))
          nj2=dsqrt(dble(ij2))
          ni1=dsqrt(dble(ii1))
          ni2=dsqrt(dble(ii2))
          mmax=min(nj2,ni2)
          do m=-mmax,mmax
             n1=max(1,iabs(m))
             n1j=max(n1,nj1)
             do 10 n=n1j,nj2
                   imn=n*(n+1)+m
                   if(imn.lt.ij1.or.imn.gt.ij2) goto 10
               do ip=1,2
                  ant(ip,n)=anpt(ip,imn)
              enddo
 10      continue
         n1i=max(n1,ni1)
         do 20 n=n1i,ni2
            imn=n*(n+1)+m
            if(imn.lt.ii1.or.imn.gt.ii2) goto 20
            a=0._dp
            b=0._dp
            do 21 l=n1j,nodrj
               ml=l*(l+1)+m
               if(ml.lt.ij1.or.ml.gt.ij2) goto 21               
               a=a+atr(1,n,ml)*ant(1,l)+ &
                atr(2,n,ml)*ant(2,l)
               b=b+atr(1,n,ml)*ant(2,l)+ &
                atr(2,n,ml)*ant(1,l)
 21         continue
            anpt(1,imn) = a
            anpt(2,imn) = b
 20      continue
      enddo
     
    end subroutine
    
!    c
!c  subroutine xuwigd.f
!c  returns the entire set of Wigner 3jm symbols in double 
!c  precision for a given integer group (j1,j2,m1,m2) 
!c  needs lnfacd.f to compute ln(z!)
!c  using the algorithm described in [Xu, Journal of Computational
!c  Physics 139, 137-165 (1998)]
!c 
      subroutine xuwigd(j1,j2,m1,m2,c,cf,n,kmax)
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: xuwigd
          include 'tmatrix_mps_np.inc'
          integer(kind=int4) :: j1,j2,m1,m2,n,kmax
          real(kind=dp), dimension(n) :: c,cf
          real(kind=dp), dimension(0:4*(np+1)) :: fnr
          ! Locals
          integer(kind=int4) :: i,m3,j,j3f,j3,k
          real(kind=dp) :: s,t,cr,a,b
          real(kind=dp), parameter :: small = 0.000000000001_dp
          common/fnr/fnr
          !data small/1.0E-12_dp/
         ! code ....
         a(j1,j2,j3,m1,m2)=fnr(j3+j1-j2)*fnr(j1+j2+1+j3)* &
                           fnr(j3+m1+m2)*fnr(j3-j1+j2)*   &
                           fnr(j1+j2+1-j3)*fnr(j3-m1-m2)
         b(j1,j2,j3,m1,m2)=dble(2*j3+1)*(dble((m1+m2)*  &
                           (j1*j1+j1-j2*j2-j2))-        &
                           dble((m1-m2)*(j3*j3+j3)))
        do i=1,n
           c(i)=0.0_dp
        enddo
        if(iabs(m1).gt.j1.or.iabs(m2).gt.j2) return
        kmax=j1+j2-max(iabs(j1-j2),iabs(m1+m2))+1
        if(kmax.gt.n) stop
        j3=j1+j2
        m3=-m1-m2
        s=lnf(dble(2*j1))+lnf(dble(2*j2))
        s=s-lnf(dble(2*j1+2*j2+1))+lnf(dble(j1+j2-m1-m2))
        s=s-lnf(dble(j1-m1))-lnf(dble(j1+m1))
        s=s+lnf(dble(j1+j2+m1+m2))-lnf(dble(j2-m2))
        s=s-lnf(dble(j2+m2))
        j=j1-j2+m1+m2
        c(1)=((-1)**j)*dexp(0.5_dp*s)
        if(kmax.eq.1) return
        c(2)=-b(j1,j2,j3,m1,m2)*c(1)/   &
              (dble(j3+1)*a(j1,j2,j3,m1,m2))
        if(kmax.eq.2) return
        j3f=max(iabs(j1-j2),iabs(m1+m2))
        j3=j3f
        if(j3.eq.0) then
            c(kmax)=dble((-1)**(j1-m1))/fnr(2*j1+1)
        else
	        s=lnf(dble(j3+j1-j2))+lnf(dble(j3-j1+j2))
	        s=s+lnf(dble(j1+j2-j3))+lnf(dble(j1-m1))
	        s=s+lnf(dble(j1+m1))+lnf(dble(j2-m2))
	        s=s+lnf(dble(j2+m2))-lnf(dble(j1+j2+j3+1))
	        s=s+lnf(dble(j3-m3))+lnf(dble(j3+m3))
	        t=0.5_dp*s
	        if(j3.eq.j1-j2.or.j3.eq.m3) k=j2+m2
	        if(j3.eq.j2-j1.or.j3.eq.-m3) k=j1-m1
	        s=lnf(dble(j1+j2-j3-k))+lnf(dble(j1-m1-k))
	        s=s+lnf(dble(k))+lnf(dble(j2+m2-k))
	        s=s+lnf(dble(j3-j2+m1+k))+lnf(dble(j3-j1-m2+k))
	        j=j1-j2-m3
	        c(kmax)=((-1)**(k+j))*dexp(t-s)           
       endif
       if(kmax.eq.3) return
       if(j3.eq.0) then
         c(kmax-1)=c(kmax)*dble(m1)/(fnr(j1)*fnr(j1+1))
       else
         c(kmax-1)=-b(j1,j2,j3,m1,m2)*c(kmax)/(dble(j3)*   &
                   a(j1,j2,j3+1,m1,m2))
       endif
       if(kmax.eq.4) return
         cf(kmax)=c(kmax)
         cf(kmax-1)=c(kmax-1)
         j3=j3f
       do i=kmax-2,3,-1
         j3=j3+1
         cf(i)=dble(j3+1)*a(j1,j2,j3,m1,m2)*cf(i+2)+   &
                    b(j1,j2,j3,m1,m2)*cf(i+1)
         cf(i)=-cf(i)/(dble(j3)*a(j1,j2,j3+1,m1,m2))
      enddo
      j3=j1+j2
      do i=3,kmax-2
         j3=j3-1
         c(i)=dble(j3)*a(j1,j2,j3+1,m1,m2)*c(i-2)+    &
              b(j1,j2,j3,m1,m2)*c(i-1)
         c(i)=-c(i)/(dble(j3+1)*a(j1,j2,j3,m1,m2))
         cr=dabs(c(i)-cf(i))
         if(dabs(c(i)).gt.1.0_dp) cr=cr/dabs(c(i))
         if(cr.lt.small) then
            do j=i,kmax-2
               c(j)=cf(j)
            enddo
            return
         endif
      enddo
     
      end  subroutine
      
    
    
!    C  The following is a part of Mishchenko's public domain code of 
!C  ampld.new.f, which is available at URL: www.giss.nasa.gov/~crmin. 
!C  For detailed information on the formulation and algorithms,  
!C  please see the original code by Mishchenko.

    subroutine tm0d(LAM,NP,EPS,AXI,RAT,MRR,MRI,DDELT,NDGS,NMAX)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: tm0d
          real(kind=dp) :: LAM
          integer(kind=int4) :: NP
          real(kind=dp)    :: EPS,AXI,RAT,MRR,MRI,DDELT
          integer(kind=int4) :: NDGS,NMAX
          ! Locals
          real(kind=dp), dimension(NPNG2) :: X,W,S,SS,R,DR,DDR,DRR,DRI
          real(kind=dp), dimension(NPN1)  :: AN
          real(kind=dp), dimension(NPN1,NPN1) :: ANN
!DIR$     ATTRIBUTES ALIGN : 64 :: X,W,S,SS,R,DR,DDR,DRR,DRI
!DIR$     ATTRIBUTES ALIGN : 64 :: AN,ANN
          real(kind=dp), dimension(NPN2,NPN2) :: TR1,TI1
          real(kind=dp), dimension(NPN6,NPN4,NPN4) :: RT11,RT12,RT21,RT22, &
                                                   IT11,IT12,IT21,IT22
          integer(kind=int4) :: ICHOICE,NCHECK, IXXX,INM1,NMA,MMAX,NGAUSS,  &
                           N,N1,NNNGGG,NGGG,NGAUS,NNM,N2,NN2,NN1,M,NM, &
                           N11,N22,M1
          real(kind=dp) :: P,A,XEV,QEXT1,QSCA1,QEXT,QSCA,TR1NN,TI1NN, TR1NN1, &
                        TI1NN1,DN1,DSCA,DEXT,PPI,PIR,PII,ZZ1,ZZ2,ZZ3,ZZ4,ZZ5, &
                        ZZ6,ZZ7,ZZ8,QSC,QXT,WALB
 
         COMMON /CT/ TR1,TI1
         COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
         COMMON /CHOICE/ ICHOICE
 
!C  OPEN FILES *******************************************************
 
      P=DACOS(-1._dp)
      ICHOICE=2
      NCHECK=0
      IF (NP.EQ.-1.OR.NP.EQ.-2) NCHECK=1
      IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK=1
      WRITE (6,5454) ICHOICE,NCHECK
 5454 FORMAT ('ICHOICE=',I1,'  NCHECK=',I1)
      IF (DABS(RAT-1.0_dp).GT.1.0E-8_dp.AND.NP.EQ.-1) CALL SAREA (EPS,RAT)
      IF (DABS(RAT-1.0_dp).GT.1.0E-8_dp.AND.NP.GE.0) CALL SURFCH(NP,EPS,RAT)
      IF (DABS(RAT-1.0_dp).GT.1.0E-8_dp.AND.NP.EQ.-2) CALL SAREAC (EPS,RAT)
      IF (NP.EQ.-3) CALL DROP (RAT)
!C     PRINT 8000, RAT
 8000 FORMAT ('RAT=',F8.6)
      IF(NP.EQ.-1.AND.EPS.GE.1._dp) PRINT 7000,EPS
      IF(NP.EQ.-1.AND.EPS.LT.1._dp) PRINT 7001,EPS
      IF(NP.GE.0) PRINT 7100,NP,EPS
      IF(NP.EQ.-2.AND.EPS.GE.1._dp) PRINT 7150,EPS
      IF(NP.EQ.-2.AND.EPS.LT.1._dp) PRINT 7151,EPS
      IF(NP.EQ.-3) PRINT 7160
      PRINT 7400, LAM,MRR,MRI
      PRINT 7200,DDELT
 7000 FORMAT('OBLATE SPHEROIDS, A/B=',F11.7)
 7001 FORMAT('PROLATE SPHEROIDS, A/B=',F11.7)
 7100 FORMAT('CHEBYSHEV PARTICLES, T',      &
           I1,'(',F5.2,')')
 7150 FORMAT('OBLATE CYLINDERS, D/L=',F11.7)
 7151 FORMAT('PROLATE CYLINDERS, D/L=',F11.7)
 7160 FORMAT('GENERALIZED CHEBYSHEV PARTICLES')
 7200 FORMAT ('ACCURACY OF COMPUTATIONS DDELT = ',D8.2)
 7400 FORMAT('LAM=',F10.6,3X,'MRR=',D10.4,3X,'MRI=',D10.4)
      DDELT=0.1_dp*DDELT
      IF (DABS(RAT-1._dp).LE.0.000001_dp) PRINT 8003, AXI
      IF (DABS(RAT-1._dp).GT.0.000001_dp) PRINT 8004, AXI
 8003 FORMAT('EQUAL-VOLUME-SPHERE RADIUS=',F8.4)
 8004 FORMAT('EQUAL-SURFACE-AREA-SPHERE RADIUS=',F8.4)
      A=RAT*AXI
      XEV=2.0_dp*P*A/LAM
      IXXX=XEV+4.05_dp*XEV**0.33333333333333333_dp
      INM1=MAX0(4,IXXX)
      IF (INM1.GE.NPN1) PRINT 7333, NPN1
      IF (INM1.GE.NPN1) STOP
 7333 FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,      &
           '.  EXECUTION TERMINATED')
      QEXT1=0.0_dp
      QSCA1=0.0_dp
      DO 50 NMA=INM1,NPN1
         NMAX=NMA
         MMAX=1
         NGAUSS=NMAX*NDGS
         IF (NGAUSS.GT.NPNG1) PRINT 7340, NGAUSS
         IF (NGAUSS.GT.NPNG1) STOP
 7340    FORMAT('NGAUSS =',I3,' I.E. IS GREATER THAN NPNG1.',    &
             '  EXECUTION TERMINATED')
 7334    FORMAT(' NMAX =', I3,'  DC2=',D8.2,'   DC1=',D8.2)
         CALL CONST(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
         CALL VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R,    &
               DR,DDR,DRR,DRI,NMAX)
         CALL TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,      &
                  DDR,DRR,DRI,NMAX,NCHECK)
         QEXT=0._dp
         QSCA=0._dp
         DO 4 N=1,NMAX
            N1=N+NMAX
            TR1NN=TR1(N,N)
            TI1NN=TI1(N,N)
            TR1NN1=TR1(N1,N1)
            TI1NN1=TI1(N1,N1)
            DN1=2*N+1
            QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN+ &
                          TR1NN1*TR1NN1+TI1NN1*TI1NN1)
            QEXT=QEXT+(TR1NN+TR1NN1)*DN1
    4    CONTINUE
         DSCA=DABS((QSCA1-QSCA)/QSCA)
         DEXT=DABS((QEXT1-QEXT)/QEXT)
         QEXT1=QEXT
         QSCA1=QSCA
!c        PRINT 7334, NMAX,DSCA,DEXT
         IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 55
         IF (NMA.EQ.NPN1) PRINT 7333, NPN1
         IF (NMA.EQ.NPN1) STOP      
   50 CONTINUE
   55 NNNGGG=NGAUSS+1
      MMAX=NMAX
      IF (NGAUSS.EQ.NPNG1) PRINT 7336
      IF (NGAUSS.EQ.NPNG1) GO TO 155 
      DO 150 NGAUS=NNNGGG,NPNG1
         NGAUSS=NGAUS
         NGGG=2*NGAUSS
 7336    FORMAT('WARNING: NGAUSS=NPNG1')
 7337    FORMAT(' NG=',I3,'  DC2=',D8.2,'   DC1=',D8.2)
         CALL CONST(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
         CALL VARY(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R,   &
                 DR,DDR,DRR,DRI,NMAX)
         CALL TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,      &
                    DDR,DRR,DRI,NMAX,NCHECK)
         QEXT=0.0_dp
         QSCA=0.0_dp
         DO 104 N=1,NMAX
            N1=N+NMAX
            TR1NN=TR1(N,N)
            TI1NN=TI1(N,N)
            TR1NN1=TR1(N1,N1)
            TI1NN1=TI1(N1,N1)
            DN1=2*N+1
            QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN+ &
                          TR1NN1*TR1NN1+TI1NN1*TI1NN1)
            QEXT=QEXT+(TR1NN+TR1NN1)*DN1
  104    CONTINUE
         DSCA=DABS((QSCA1-QSCA)/QSCA)
         DEXT=DABS((QEXT1-QEXT)/QEXT)
!c        PRINT 7337, NGGG,DSCA,DEXT
         QEXT1=QEXT
         QSCA1=QSCA
         IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 155
         IF (NGAUS.EQ.NPNG1) PRINT 7336
  150 CONTINUE
  155 CONTINUE
      QSCA=0._dp
      QEXT=0._dp
      NNM=NMAX*2
      DO 204 N=1,NNM
         QEXT=QEXT+TR1(N,N)
  204 CONTINUE
      DO 213 N2=1,NMAX
         NN2=N2+NMAX
         DO 213 N1=1,NMAX
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
  213 CONTINUE
!c     PRINT 7800,0,DABS(QEXT),QSCA,NMAX
      DO 220 M=1,NMAX
         CALL TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,     &
                    DDR,DRR,DRI,NMAX,NCHECK)
         NM=NMAX-M+1
         M1=M+1
         QSC=0._dp
         DO 214 N2=1,NM
            NN2=N2+M-1
            N22=N2+NM
            DO 214 N1=1,NM
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
               QSC=QSC+(ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4+   &
                      ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8)*2._dp
  214    CONTINUE
         NNM=2*NM
         QXT=0._dp
         DO 215 N=1,NNM
            QXT=QXT+TR1(N,N)*2._dp
  215    CONTINUE
         QSCA=QSCA+QSC
         QEXT=QEXT+QXT
!c        PRINT 7800,M,DABS(QXT),QSC,NMAX
 7800    FORMAT(' m=',I3,'  qxt=',D12.6,'  qsc=',D12.6,    &
             '  nmax=',I3)
  220 CONTINUE
      WALB=-QSCA/QEXT
      IF (WALB.GT.1_dp+DDELT) PRINT 9111
 9111 FORMAT ('WARNING: W IS GREATER THAN 1')

!c      ITIME=MCLOCK()
!c      TIME=DFLOAT(ITIME)/6000D0
!c      PRINT 1001,TIME
! 1001 FORMAT (' time =',F8.2,' min')
   
    end subroutine
    
!    C*****************************************************************
!C
!C     Calculation of the functions
!C     DV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)
!C     and
!C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x)
!C     1.LE.N.LE.NMAX
!C     0.LE.X.LE.1

    SUBROUTINE VIGAMPL (X, NMAX, M, DV1, DV2)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: VIGAMPL
          real(kind=dp) :: X
          integer(kind=int4) :: NMAX,M
          real(kind=dp), dimension(NPN6) ::  DV1, DV2
          ! Locals
          integer(kind=int4) :: N,I,I2
          real(kind=dp)    :: DX,A,QS,QS1,DSI,D1,D2,QN,QN1,QN2,  &
                           D3,DER,QNM,QNM1,QMM,DN
         DO 1 N=1,NMAX
             DV1(N)=0._dp
             DV2(N)=0._dp
       1 CONTINUE
         DX=DABS(X)
         IF (DABS(1._dp-DX).LE.0.0000000001_dp) GO TO 100
         A=1._dp
         QS=DSQRT(1._dp-X*X)
         QS1=1._dp/QS
         DSI=QS1
         IF (M.NE.0) GO TO 20
         D1=1.0-dp
         D2=X  
         DO 5 N=1,NMAX  
             QN=DFLOAT(N)
             QN1=DFLOAT(N+1)
             QN2=DFLOAT(2*N+1)
             D3=(QN2*X*D2-QN*D1)/QN1 
             DER=QS1*(QN1*QN/QN2)*(-D1+D3)
             DV1(N)=D2*DSI
             DV2(N)=DER
             D1=D2
             D2=D3
    5 CONTINUE
      RETURN
   20 QMM=DFLOAT(M*M)
      DO 25 I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS
   25 CONTINUE
      D1=0._dp
      D2=A 
      DO 30 N=M,NMAX
         QN=DFLOAT(N)
         QN2=DFLOAT(2*N+1)
         QN1=DFLOAT(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2
         DV1(N)=D2*DSI
         DV2(N)=DER
         D1=D2
         D2=D3
   30 CONTINUE
      RETURN
  100 IF (M.NE.1) RETURN
      DO 110 N=1,NMAX
         DN=DFLOAT(N*(N+1))
         DN=0.5_dp*DSQRT(DN)
         IF (X.LT.0._dp) DN=DN*(-1)**(N+1)
         DV1(N)=DN
         IF (X.LT.0._dp) DN=-DN
         DV2(N)=DN
  110 CONTINUE
    
    end subroutine
    
    SUBROUTINE CONST (NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: CONST
          integer(kind=int4) :: NGAUSS,NMAX,MMAX,NP
          real(kind=dp) :: P,EPS
          real(kind=dp), dimension(NPNG2) ::  X,W,S,SS
          real(kind=dp), dimension(NPNG1) ::  X1,W1,X2,W2
!DIR$     ATTRIBUTES ALIGN : 64 :: X1,W1,X2,W2
          real(kind=dp), dimension(NPN1) :: AN
          real(kind=dp), dimension(NPN1,NPN1) :: ANN
          real(kind=dp), dimension(NPN1) :: DD
!DIR$     ATTRIBUTES ALIGN : 64 :: DD          
          integer(kind=int4) :: N,NN,N1,NG,NG1,NG2,I
          real(kind=dp)    :: D,DDD,XX,Y
          ! Exec code ...
          DO 10 N=1,NMAX
                NN=N*(N+1)
                AN(N)=DFLOAT(NN)
                D=DSQRT(DFLOAT(2*N+1)/DFLOAT(NN))
                DD(N)=D
                DO 10 N1=1,N
                     DDD=D*DD(N1)*0.5_dp
                     ANN(N,N1)=DDD
                     ANN(N1,N)=DDD
        10 CONTINUE
           NG=2*NGAUSS
           IF (NP.EQ.-2) GO  TO 11
           CALL GAUSS(NG,0,0,X,W)
           GO TO 19
        11 NG1=DFLOAT(NGAUSS)/2._dp
           NG2=NGAUSS-NG1
           XX=-DCOS(DATAN(EPS))
           CALL GAUSS(NG1,0,0,X1,W1)
           CALL GAUSS(NG2,0,0,X2,W2)
!DIR$      SIMD
           DO 12 I=1,NG1
                 W(I)=0.5_dp*(XX+1.0_dp)*W1(I)
                 X(I)=0.5_dp*(XX+1.0_dp)*X1(I)+0.5_dp*(XX-1._dp)
        12 CONTINUE
         DO 14 I=1,NG2
              W(I+NG1)=-0.5_dp*XX*W2(I)
              X(I+NG1)=-0.5_dp*XX*X2(I)+0.5_dp*XX
       14 CONTINUE
       DO 16 I=1,NGAUSS
             W(NG-I+1)=W(I)
             X(NG-I+1)=-X(I)
      16 CONTINUE
      19 DO 20 I=1,NGAUSS
               Y=X(I)
               Y=1D0/(1.0_dp-Y*Y)
               SS(I)=Y
               SS(NG-I+1)=Y
               Y=DSQRT(Y)
               S(I)=Y
               S(NG-I+1)=Y
      20 CONTINUE
     
    end subroutine
    
    SUBROUTINE VARY (LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,    &
                      R,DR,DDR,DRR,DRI,NMAX)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: VARY
          real(kind=dp) :: LAM,MRR,MRI,A,EPS
          integer(kind=int4) :: NP,NGAUSS
          real(kind=dp), dimension(NPNG2) :: X
          real(kind=dp) :: P,PPI,PIR,PII
          real(kind=dp), dimension(NPNG2) :: R,DR,DDR,DRR,DRI
          integer(kind=int4) :: NMAX
          ! Locals
          real(kind=dp), dimension(NPNG2) :: Z,ZR,ZI
!DIR$     ATTRIBUTES ALIGN : 64 :: Z,ZR,ZI
          real(kind=dp), dimension(NPNG2,NPN1) :: J,Y,JR,JI,DJ,DJR,DJI,DY
!!DIR$     ATTRIBUTES ALIGN : 64 :: J,Y,JR,JI,DJ,DJR,DJI,DY      
          integer(kind=int4) :: NG,I,NNMAX1, NNMAX2
          real(kind=dp) :: PI,V,PRR,PRI,TA,VV,V1,V2,TB
          COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
          NG=NGAUSS*2
          IF (NP.GT.0) CALL RSP2(X,NG,A,EPS,NP,R,DR)
          IF (NP.EQ.-1) CALL RSP1(X,NG,NGAUSS,A,EPS,NP,R,DR)
          IF (NP.EQ.-2) CALL RSP3(X,NG,NGAUSS,A,EPS,R,DR)
          IF (NP.EQ.-3) CALL RSP4(X,NG,A,R,DR)
          PI=P*2._dp/LAM
          PPI=PI*PI
          PIR=PPI*MRR
          PII=PPI*MRI
          V=1._dp/(MRR*MRR+MRI*MRI)
          PRR=MRR*V
          PRI=-MRI*V
          TA=0._dp
          DO 10 I=1,NG
                VV=DSQRT(R(I))
                V=VV*PI
                TA=MAX(TA,V)
                VV=1._dp/V
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
   9000 FORMAT(' NMAX = ',I2,', i.e., greater than ',I3)
      TB=TA*DSQRT(MRR*MRR+MRI*MRI)
      TB=DMAX1(TB,DFLOAT(NMAX))
      NNMAX1=1.2_dp*DSQRT(DMAX1(TA,DFLOAT(NMAX)))+3._dp
      NNMAX2=(TB+4._dp*(TB**0.33333_dp)+1.2_dp*DSQRT(TB))
      NNMAX2=NNMAX2-NMAX+5
      CALL BESS(Z,ZR,ZI,NG,NMAX,NNMAX1,NNMAX2)
    
    end subroutine
                      
    SUBROUTINE RSP1 (X,NG,NGAUSS,REV,EPS,NP,R,DR)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: RSP1
          real(kind=dp), dimension(NG)  :: X
          integer(kind=int4) :: NG,NGAUSS
          real(kind=dp) :: REV,EPS
          integer(kind=int4) :: NP
          real(kind=dp), dimension(NG) :: R,DR
          ! Locals
          real(kind=dp) :: A,AA,EE,EE1,C,CC,SS,S,RR
          integer(kind=int4) :: I
          ! Exec code ...
          A=REV*EPS**(0.3333333333333333333333_dp)
          AA=A*A
          EE=EPS*EPS
          EE1=EE-1._dp
          DO 50 I=1,NGAUSS
                C=X(I)
                CC=C*C
                SS=1.0_dp-CC
                S=DSQRT(SS)
                RR=1.0_dp/(SS+EE*CC)
                R(I)=AA*RR
                R(NG-I+1)=R(I)
                DR(I)=RR*C*S*EE1
                DR(NG-I+1)=-DR(I)
       50 CONTINUE
     
    end subroutine
    
    SUBROUTINE RSP2 (X,NG,REV,EPS,N,R,DR)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: RSP2
          real(kind=dp), dimension(NG) :: X
          integer(kind=int4) :: NG
          real(kind=dp) :: REV,EPS
          integer(kind=int4) :: N
          real(kind=dp), dimension(NG) :: R,DR
          ! Locals
          real(kind=dp) :: DNP,DN,DN4,EP,A,R0,XI,RI
          integer(kind=int4) :: I
          ! Exec code ....          
          DNP=DFLOAT(N)
          DN=DNP*DNP
          DN4=DN*4._dp
          EP=EPS*EPS
          A=1._dp+1.5_dp*EP*(DN4-2._dp)/(DN4-1._dp)
          I=(DNP+0.1_dp)*0.5_dp
          I=2*I
          IF (I.EQ.N) A=A-3._dp*EPS*(1._dp+0.25_dp*EP)/     &
                  (DN-1._dp)-0.25_dp*EP*EPS/(9._dp*DN-1._dp)
          R0=REV*A**(-0.3333333333333333333333_dp)
          DO 50 I=1,NG
                XI=DACOS(X(I))*DNP
                RI=R0*(1._dp+EPS*DCOS(XI))
                R(I)=RI*RI
                DR(I)=-R0*EPS*DNP*DSIN(XI)/RI
!c        WRITE (6,*) I,R(I),DR(I)
         50 CONTINUE
    
    end subroutine
    
      SUBROUTINE RSP3 (X,NG,NGAUSS,REV,EPS,R,DR)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: RSP3
          real(kind=dp), dimension(NG) :: X
          integer(kind=int4) :: NG,NGAUSS
          real(kind=dp) :: REV,EPS
          real(kind=dp), dimension(NG) :: R,DR
          ! Locals
          real(kind=dp) :: H,A,CO,SI,RAD,RTHET
          integer(kind=int4) :: I
          ! Exec code ....
          H=REV*( (2._dp/(3._dp*EPS*EPS))**(0.3333333333333333333333_dp) )
          A=H*EPS
          DO 50 I=1,NGAUSS
                CO=-X(I)
                SI=DSQRT(1.0_dp-CO*CO)
                IF (SI/CO.GT.A/H) GO TO 20
                RAD=H/CO
                RTHET=H*SI/(CO*CO)
                GO TO 30
   20           RAD=A/SI
                RTHET=-A*CO/(SI*SI)
   30           R(I)=RAD*RAD
                R(NG-I+1)=R(I)
                DR(I)=-RTHET/RAD
                DR(NG-I+1)=-DR(I)
     50 CONTINUE
     
      end subroutine
      
!      C**********************************************************************
!C                                                                     *
!C   Calculation of the functions R(I)=r(y)**2 and                     *
!C   DR(I)=((d/dy)r(y))/r(y) for a distorted                           *
!C   droplet specified by the parameters r_ev (equal-volume-sphere     *
!C   radius) and c_n (Chebyshev expansion coefficients)                *
!C   Y(I)=arccos(X(I))                                                 *
!C   1.LE.I.LE.NG                                                      *
!C   X - arguments                                                     *
!C                                                                     *
!C**********************************************************************

    SUBROUTINE RSP4 (X,NG,REV,R,DR)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: RSP4
          real(kind=dp), dimension(NG) :: X
          integer(kind=int4) :: NG
          real(kind=dp) :: REV
          real(kind=dp), dimension(NG) :: DR,R
          ! Locals
          integer(kind=int4), parameter :: NC = 10
          real(kind=dp), dimension(0:NC) :: C
          real(kind=dp) :: R0,R0V,XI,RI,DRI,XIN
          integer(kind=int4) :: I,N
      
          COMMON /CDROP/ C,R0V
          R0=REV*R0V
          DO I=1,NG
             XI=DACOS(X(I))
             RI=1.0_dp+C(0)
             DRI=0._dp
             DO N=1,NC
                XIN=XI*N
                RI=RI+C(N)*DCOS(XIN)
                DRI=DRI-C(N)*N*DSIN(XIN)
             ENDDO
         RI=RI*R0
         DRI=DRI*R0
         R(I)=RI*RI
         DR(I)=DRI/RI
!c        WRITE (6,*) I,R(I),DR(I)
      ENDDO
    
    end subroutine
    
    SUBROUTINE BESS (X,XR,XI,NG,NMAX,NNMAX1,NNMAX2)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: BESS
          real(kind=dp), dimension(NG) :: X,XR,XI
          integer(kind=int4) :: NG,NMAX,NNMAX1,NNMAX2
          ! lOCALS
          real(kind=dp), dimension(NPNG2,NPN1) :: J,Y,JR,JI,DJ,DY,DJR,DJI
          real(kind=dp), dimension(NPN1) :: AJ,AY,AJR,AJI,ADJ,ADY,ADJR,ADJI
!DIR$     ATTRIBUTES ALIGN : 64 ::  AJ,AY,AJR,AJI,ADJ,ADY,ADJR,ADJI  
          integer(kind=int4) :: I,N
          real(kind=dp) :: XX,YR,YI
          COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
 
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
       10 CONTINUE
    
    end subroutine
    
    SUBROUTINE RJB(X,Y,U,NMAX,NNMAX)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: RJB
          real(kind=dp)  ::  X
          real(kind=dp), dimension(NMAX) :: Y,U
          integer(kind=int4) :: NMAX,NNMAX
          ! Locals
          real(kind=dp), dimension(800) :: Z
!DIR$     ATTRIBUTES ALIGN : 64 :: Z          
          integer(kind=int4) :: L,LI,L1,I1,I
          real(kind=dp) :: XX,Z0,Y0,Y1,YI1,YI
          ! Exec code ....
          L=NMAX+NNMAX
          XX=1._dp/X
          Z(L)=1._dp/(DFLOAT(2*L+1)*XX)
          L1=L-1
          DO 5 I=1,L1
               I1=L-I
               Z(I1)=1._dp/(DFLOAT(2*I1+1)*XX-Z(I1+1))
        5 CONTINUE
        Z0=1._dp/(XX-Z(1))
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
    
    SUBROUTINE RYB(X,Y,V,NMAX)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: RYB
          real(kind=dp) :: X
          real(kind=dp), dimension(NMAX) :: Y,V
          integer(kind=int4) :: NMAX
          ! Locals
          integer(kind=int4) :: I,NMAX1
          real(kind=dp) :: C,S,X1,X2,X3,Y1
          C=DCOS(X)
          S=DSIN(X)
          X1=1._dp/X
          X2=X1*X1
          X3=X2*X1
          Y1=-C*X2-S*X1
          Y(1)=Y1
          Y(2)=(-3._dp*X3+X1)*C-3._dp*X2*S
          NMAX1=NMAX-1
          DO 5 I=2,NMAX1
        5      Y(I+1)=DFLOAT(2*I+1)*X1*Y(I)-Y(I-1)
          V(1)=-X1*(C+Y1)
          DO 10 I=2,NMAX
  10         V(I)=Y(I-1)-DFLOAT(I)*X1*Y(I)
     
    END SUBROUTINE
    
!    C**********************************************************************
!C                                                                     *
!C   CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND       *
!C   J=JR+I*JI OF COMPLEX ARGUMENT X=XR+I*XI OF ORDERS FROM 1 TO NMAX  *
!C   BY USING BACKWARD RECURSION. PARAMETR NNMAX DETERMINES NUMERICAL  *
!C   ACCURACY. U=UR+I*UI - FUNCTION (1/X)(D/DX)(X*J(X))                *
!C                                                                     *
!C**********************************************************************
 
    SUBROUTINE CJB (XR,XI,YR,YI,UR,UI,NMAX,NNMAX)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: CJB
          real(kind=dp) :: XR,XI
          real(kind=dp), dimension(NMAX) :: YR,YI,UR,UI
          integer(kind=int4) :: NMAX,NNMAX
          real(kind=dp), dimension(NPN1) :: CYR,CYI,CUR,CUI
!DIR$     ATTRIBUTES ALIGN : 64 :: CYR,CYI,CUR,CUI
          real(kind=dp), dimension(1200) :: CZR,CZI
!DIR$     ATTRIBUTES ALIGN : 64 :: CZR,CZI
          integer(kind=int4) :: L,L1,I,I1
          real(kind=dp) :: XRXI,CXXR,CXXI,QF,AR,AI,ARI,CZ0R,CZ0I, &
                        CR,CI,CY0R,CY0I,CY1R,CY1I,CU1R,CU1I,   &
                        QI,CYI1R, CYI1I,CUII,CYIR,CUIR,CYII
          L=NMAX+NNMAX
          XRXI=1._dp/(XR*XR+XI*XI)
          CXXR=XR*XRXI
          CXXI=-XI*XRXI 
          QF=1._dp/DFLOAT(2*L+1)
          CZR(L)=XR*QF
          CZI(L)=XI*QF
          L1=L-1
          DO I=1,L1
             I1=L-I
             QF=DFLOAT(2*I1+1)
             AR=QF*CXXR-CZR(I1+1)
             AI=QF*CXXI-CZI(I1+1)
             ARI=1._dp/(AR*AR+AI*AI)
             CZR(I1)=AR*ARI
             CZI(I1)=-AI*ARI
         ENDDO   
         AR=CXXR-CZR(1)
         AI=CXXI-CZI(1)
         ARI=1._dp/(AR*AR+AI*AI)
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
     
    END  SUBROUTINE
    
    SUBROUTINE TMATR0 (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,   &
                       DRR,DRI,NMAX,NCHECK)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: TMATR0        
          integer(kind=int4) :: NGAUSS
          real(kind=dp), dimension(NPNG2) :: X,W
          real(kind=dp), dimension(NPN1)  :: AN
          real(kind=dp), dimension(NPN1,NPN1) :: ANN
          real(kind=dp), dimension(NPNG2) :: S,SS
          real(kind=dp) :: PPI,PIR,PII
          real(kind=dp), dimension(NPNG2) :: R,DR,DDR,DRR,DRI
          integer(kind=int4) :: NMAX,NCHECK
          ! Locals
          real(kind=dp), dimension(NPN2) :: SIG
!DIR$     ATTRIBUTES ALIGN : 64 :: SIG
          real(kind=dp), dimension(NPNG2,NPN1) :: J,Y,JR,JI,DJ,DY,DJR,DJI,D1,D2
          real(kind=dp), dimension(NPNG2) :: DS,DSS,RR
          real(kind=dp), dimension(NPN1) :: DV1,DV2
!DIR$     ATTRIBUTES ALIGN : 64 :: DV1,DV2   
          real(kind=dp), dimension(NPN1,NPN1) :: R11,R12,R21,R22,I11,I12,I21,I22, &
                                              RG11,RG12,RG21,RG22,IG11,IG12,   &
                                              IG21,IG22
          real(kind=dp), dimension(NPN2,NPN2) :: QR,QI,RGQR,RGQI,TQR,TQI,TRGQR,TRGQI
          real(kind=dp), dimension(NPN2,NPN2) :: TR1,TI1
          integer(kind=int4) :: MM1,NNMAX,NG,NGSS,N,I,I1,I2,N1,NM,K1,KK1,N2,K2,KK2, &
                                AN12
          real(kind=dp) :: FACTOR,SI,DD1,DD2,AN2,AR12,AR21,AI12,AI21,GR12,GR21, &
                        GI12,GI21,D1N1,D2N1,D1N2,D2N2,AA1,QJ1,QY1,QJR2,QJI2,QDJR2, &
                        QDJI2,QDJ1,QDY1,C1R,C1I,B1R,B1I,C2R,C2I,B2R,B2I,C3R,C3I,   &
                        B3R,B3I,C4R,C4I,B4R,B4I,DRRI,DRII,C5R,C5I,B5R,B5I,URI,RRI, &
                        F1,F2,AN1,TPIR,TPII,TPPI,TAR12,TAI12,TGR12,TGI12,TAR21,TAI21, &
                        TGR21,TGI21,A12,A21,A22,DDRI
          COMMON /TMAT99/           &
                 R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,      &
                 IG11,IG12,IG21,IG22
          COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
          COMMON /CT/ TR1,TI1
          COMMON /CTT/ QR,QI,RGQR,RGQI
          MM1=1
          NNMAX=NMAX+NMAX
          NG=2*NGAUSS
          NGSS=NG
          FACTOR=1._dp
          IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2._dp
          ELSE
            CONTINUE
         ENDIF
         SI=1._dp
         DO 5 N=1,NNMAX
           SI=-SI
           SIG(N)=SI
       5 CONTINUE
      20 DO 25 I=1,NGAUSS
               I1=NGAUSS+I
               I2=NGAUSS-I+1
               CALL VIG ( X(I1), NMAX, 0, DV1, DV2)
          DO 25 N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI
       25 CONTINUE
     30 DO 40 I=1,NGSS
           RR(I)=W(I)*R(I)
   40 CONTINUE
 
      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR12=0._dp
                AR21=0._dp
                AI12=0._dp
                AI21=0._dp
                GR12=0._dp
                GR21=0._dp
                GI12=0._dp
                GI21=0._dp
                IF (NCHECK.EQ.1.AND.SIG(N1+N2).LT.0._dp) GO TO 205

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
 
                TQR(K1,KK2)=0._dp
                TQI(K1,KK2)=0._dp
                TRGQR(K1,KK2)=0._dp
                TRGQI(K1,KK2)=0._dp
 
                TQR(KK1,K2)=0._dp
                TQI(KK1,K2)=0._dp
                TRGQR(KK1,K2)=0._dp
                TRGQI(KK1,K2)=0._dp
 
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
     
    END SUBROUTINE
                       
      subroutine TMATR (M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,    &
                       DRR,DRI,NMAX,NCHECK)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: TMATR
          integer(kind=int4) :: M,NGAUSS
          real(kind=dp), dimension(NPNG2) :: X,W,S,SS
          real(kind=dp), dimension(NPN1)  :: AN
          real(kind=dp), dimension(NPN1,NPN1) :: ANN

          real(kind=dp) :: PPI,PIR,PII
          real(kind=dp), dimension(NPNG2) :: R,DR,DDR,DRR,DRI

          integer(kind=int4) :: NMAX,NCHECK
          ! Locals
          real(kind=dp), dimension(NPNG2) :: SIG 
!DIR$     ATTRIBUTES ALIGN : 64 :: SIG
          real(kind=dp), dimension(NPNG2,NPN1) :: J,Y,JR,JI,DJ,DY,DJR,  &
                                               DJI,D1,D2
!DIR$     ATTRIBUTES ALIGN : 64 :: D1
!DIR$     ATTRIBUTES ALIGN : 64 :: D2
          real(kind=dp), dimension(NPNG2) :: DS,DSS,RR,DV1,DV2
!DIR$     ATTRIBUTES ALIGN : 64 :: DS
!DIR$     ATTRIBUTES ALIGN : 64 :: DSS
!DIR$     ATTRIBUTES ALIGN : 64 :: RR
!DIR$     ATTRIBUTES ALIGN : 64 :: DV1
!DIR$     ATTRIBUTES ALIGN : 64 :: DV2          
          real(kind=dp), dimension(NPN1,NPN1) :: R11,R12,R21,R22,I11,I12,  &
                                              I21,I22,RG11,RG12,RG21,  &
                                              RG22,IG11,IG12,IG21,IG22
          real(kind=dp), dimension(NPN2,NPN2) :: QR,QI,RGQR,RGQI,TQR,TQI, &
                                              TRGQR,TRGQI,TR1,TI1
!DIR$     ATTRIBUTES ALIGN : 64 :: TQR          
!DIR$     ATTRIBUTES ALIGN : 64 :: TQI
!DIR$     ATTRIBUTES ALIGN : 64 :: TRGQR
!DIR$     ATTRIBUTES ALIGN : 64 :: TRGQI
          real(kind=sp), dimension(NPN6*NPN4*NPN4*8) :: PLUS

          integer(kind=int4) :: MM1,NNMAX,NG,NGSS,N,I,I1,I2,   &
                           N1,N2,K1,KK1,K2,KK2,NM
          real(kind=dp)    :: FACTOR,QM,QMM,WR,SI,DD1,DD2,AN1,AN2,     &
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
                           
      COMMON /TMAT99/ PLUS,   &
                 R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,  &
                 IG11,IG12,IG21,IG22
!DIR$   ATTRIBUTES ALIGN : 64 :: /TMAT99/
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
!DIR$   ATTRIBUTES ALIGN : 64 :: /CBESS/
      COMMON /CT/ TR1,TI1
!DIR$   ATTRIBUTES ALIGN : 64 :: /CT/
      COMMON /CTT/ QR,QI,RGQR,RGQI
!DIR$   ATTRIBUTES ALIGN : 64 :: /CTT/      
      MM1=M
      QM=DFLOAT(M)
      QMM=QM*QM
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1._dp
       IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2.0_dp
         ELSE
            CONTINUE
      ENDIF
      SI=1._dp
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
                AR11=0._dp
                AR12=0._dp
                AR21=0._dp
                AR22=0._dp
                AI11=0._dp
                AI12=0._dp
                AI21=0._dp
                AI22=0._dp
                GR11=0._dp
                GR12=0._dp
                GR21=0._dp
                GR22=0._dp
                GI11=0._dp
                GI12=0._dp
                GI21=0._dp
                GI22=0._dp
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
 
                    IF (NCHECK.EQ.1.AND.SI.GT.0._dp) GO TO 150
 
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
                       
    SUBROUTINE VIG (X, NMAX, M, DV1, DV2)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: VIG
          real(kind=dp) :: X
          integer(kind=int4) :: NMAX,M
          real(kind=dp), dimension(NPN1) :: DV1,DV2
          ! Locals
          real(kind=dp) :: A,QS,QS1,D1,D2,QN,QN1,QN2,D3,DER,  &
                        QMM,QNM,QNM1
          integer(kind=int4) :: N,I,I2
           ! Exec code ....
          A=1._dp
          QS=DSQRT(1._dp-X*X)
          QS1=1._dp/QS
          DO N=1,NMAX
                DV1(N)=0._dp
                DV2(N)=0._dp
         ENDDO   
      IF (M.NE.0) GO TO 20
      D1=1._dp
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
      D1=0._dp
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
      
    END SUBROUTINE
    
!    C**********************************************************************
!C                                                                     *
!C   CALCULATION OF THE MATRIX    T = - RG(Q) * (Q**(-1))              *
!C                                                                     *
!C   INPUT INFORTMATION IS IN COMMON /CTT/                             *
!C   OUTPUT INFORMATION IS IN COMMON /CT/                              *
!C                                                                     *
!C**********************************************************************
 
    SUBROUTINE TT(NMAX,NCHECK)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: TT
          integer(kind=int4) :: NMAX,NCHECK
          real(kind=dp), dimension(NPN2,NPN2) :: F
!DIR$     ATTRIBUTES ALIGN : 64 :: F
          real(kind=dp), dimension(NPN2) :: B,WORK
!DIR$     ATTRIBUTES ALIGN : 64 :: B,WORK
          real(kind=dp), dimension(NPN2,NPN2) :: QR,QI,RGQR,RGQI,  &
                                              A,C,D,E
!DIR$     ATTRIBUTES ALIGN : 64 ::  A,C,D,E
          real(kind=dp), dimension(NPN2,NPN2) :: TR1,TI1
          COMPLEX(16), dimension(NPN2,NPN2) ::  ZQ
          COMPLEX(16), dimension(NPN2) :: ZW
!DIR$     ATTRIBUTES ALIGN : 64 :: ZQ,ZW
          COMPLEX(16), dimension(NPN2,NPN2) :: ZQR,ZAFAC,ZT,       &
                                               ZTHETA
!DIR$     ATTRIBUTES ALIGN : 64 :: ZQR,ZAFAC,ZT,ZTHETA     
     
         INTEGER(kind=int4), dimension(NPN2) ::  IPIV,IPVT
         integer(kind=int4) :: I,NDIM,NNMAX,J,ICHOICE,INFO,K,IFAIL,N1,N2
         real(kind=dp) :: TR,TI,ARR,ARI,AR,AI,COND
         COMMON /CHOICE/ ICHOICE
         COMMON /CT/ TR1,TI1
         COMMON /CTT/ QR,QI,RGQR,RGQI
         NDIM=NPN2
         NNMAX=2*NMAX
         IF (ICHOICE.EQ.2) GO TO 5
 
!C	Inversion from NAG-LIB or Waterman's method
 
         DO I=1,NNMAX
!DIR$ SIMD
	   DO J=1,NNMAX
	      ZQ(I,J)=DCMPLX(QR(I,J),QI(I,J))
	      ZAFAC(I,J)=ZQ(I,J)
	   ENDDO
	ENDDO
	IF (ICHOICE.EQ.1) THEN
	   INFO=0
!c           CALL F07ARF(NNMAX,NNMAX,ZQ,NPN2,IPIV,INFO)
           IF (INFO.NE.0) WRITE (6,1100) INFO
!c           CALL F07AWF(NNMAX,ZQ,NPN2,IPIV,ZW,NPN2,INFO)
           IF (INFO.NE.0) WRITE (6,1100) INFO
 1100      FORMAT ('WARNING:  info=', i2)
	   DO I=1,NNMAX
	      DO J=1,NNMAX
	         TR=0._dp
	         TI=0._dp
	         DO K=1,NNMAX
                    ARR=RGQR(I,K)
                    ARI=RGQI(I,K)
                    AR=ZQ(K,J)
                    AI=IMAG(ZQ(K,J))
                    TR=TR-ARR*AR+ARI*AI
                    TI=TI-ARR*AI-ARI*AR
                 ENDDO
	         TR1(I,J)=TR
	         TI1(I,J)=TI
	      ENDDO
	   ENDDO
 
	   ELSE
	   IFAIL=0
!C          CALL F01RCF(NNMAX,NNMAX,ZAFAC,NPN2,ZTHETA,IFAIL)
!C          CALL F01REF('S',NNMAX,NNMAX,NNMAX,ZAFAC,
!C    &                 NPN2,ZTHETA,ZW,IFAIL)
    DO I=1,NNMAX
!DIR$ SIMD
	      DO J=1,NNMAX
	         ZQ(I,J)=CMPLX(REAL(ZAFAC(I,J)),-    &
                      IMAG(ZAFAC(I,J)))
	      ENDDO
	   ENDDO
	   DO I=1,NNMAX
	      DO J=1,NNMAX
                 IF (I.LE.NNMAX/2.AND.I.EQ.J) THEN
	            D(I,J)=-1.0_dp
                    ELSE IF (I.GT.NNMAX/2.AND.I.EQ.J) THEN
	               D(I,J)=1._dp
	               ELSE
	               D(I,J)=0._dp
	         ENDIF
              ENDDO
           ENDDO
	   DO I=1,NNMAX
	      DO J=1,NNMAX
	         ZT(I,J)=CMPLX(0.0_dp,0.0_dp)
	         DO K=1,NNMAX
	            ZT(I,J)=ZT(I,J)+D(I,I)* &
                       ZQ(I,K)*D(K,K)*ZQ(J,K)
	         ENDDO
	         ZT(I,J)=0.5_dp*(ZT(I,J)-D(I,J)**2)
	         TR1(I,J)=REAL(ZT(I,j))
	         TI1(I,J)=IMAG(ZT(i,j))
              ENDDO
	   ENDDO
	ENDIF
		
	GOTO 70
 
!C  Gaussian elimination

    5 DO 10 N1=1,NNMAX
         DO 10 N2=1,NNMAX
            F(N1,N2)=QI(N1,N2)
   10 CONTINUE
      IF (NCHECK.EQ.1) THEN
          CALL INV1(NMAX,F,A)
        ELSE
          CALL INVERT(NDIM,NNMAX,F,A,COND,IPVT,WORK,B) 
      ENDIF
      CALL PROD(QR,A,C,NDIM,NNMAX)
      CALL PROD(C,QR,D,NDIM,NNMAX)
      DO 20 N1=1,NNMAX
!DIR$ SIMD
           DO 20 N2=1,NNMAX
                C(N1,N2)=D(N1,N2)+QI(N1,N2)
   20 CONTINUE
      IF (NCHECK.EQ.1) THEN
          CALL INV1(NMAX,C,QI)
        ELSE
          CALL INVERT(NDIM,NNMAX,C,QI,COND,IPVT,WORK,B) 
      ENDIF
      CALL PROD(A,QR,D,NDIM,NNMAX)
      CALL PROD(D,QI,QR,NDIM,NNMAX)
 
      CALL PROD(RGQR,QR,A,NDIM,NNMAX)
      CALL PROD(RGQI,QI,C,NDIM,NNMAX)
      CALL PROD(RGQR,QI,D,NDIM,NNMAX)
      CALL PROD(RGQI,QR,E,NDIM,NNMAX)
      DO 30 N1=1,NNMAX
           DO 30 N2=1,NNMAX
                TR1(N1,N2)=-A(N1,N2)-C(N1,N2)
                TI1(N1,N2)= D(N1,N2)-E(N1,N2)
   30 CONTINUE
   70 RETURN
    END SUBROUTINE
    
    SUBROUTINE PROD(A,B,C,NDIM,N)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: PROD
          real(kind=dp), dimension(NDIM,N) :: A,B,C
          integer(kind=int4) :: NDIM,N
          ! Locals
          real(kind=dp) :: CIJ
          integer(kind=int4) :: I,J,K
      DO 10 I=1,N
           DO 10 J=1,N
                CIJ=0._dp
                DO 5 K=1,N
                     CIJ=CIJ+A(I,K)*B(K,J)
    5           CONTINUE
                C(I,J)=CIJ
   10 CONTINUE
     
    END   SUBROUTINE
    
    SUBROUTINE INV1 (NMAX,F,A)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: INV1
          integer(kind=int4) :: NMAX
          real(kind=dp), dimension(NPN2,NPN2) :: F,A
          ! Locals
          real(kind=dp) :: COND
          real(kind=dp), dimension(NPN1) :: WORK
          real(kind=dp), dimension(NPN1) :: B
!DIR$     ATTRIBUTES ALIGN : 64 :: WORK,B
          real(kind=dp), dimension(NPN1,NPN1) :: Q1,Q2,P1,P2
!DIR$     ATTRIBUTES ALIGN : 64 :: Q1,Q2,P1,P2
          integer(kind=int4), dimension(NPN1) ::  IPVT,IND1,IND2
          integer(kind=int4) :: NDIM,NN1,NN2,I,NNMAX,J,I1,I2,J1,J2
          ! Exec code ....
          NDIM=NPN1
          NN1=(DFLOAT(NMAX)-0.1_dp)*0.5_dp+1._dp 
          NN2=NMAX-NN1
          DO 5 I=1,NMAX
                IND1(I)=2*I-1
                IF(I.GT.NN1) IND1(I)=NMAX+2*(I-NN1)
                   IND2(I)=2*I
                IF(I.GT.NN2) IND2(I)=NMAX+2*(I-NN2)-1
         5 CONTINUE
           NNMAX=2*NMAX
           DO 15 I=1,NMAX
                 I1=IND1(I)
                 I2=IND2(I)
              DO 15 J=1,NMAX
                    J1=IND1(J)
                    J2=IND2(J)
                    Q1(J,I)=F(J1,I1)
                    Q2(J,I)=F(J2,I2)
         15 CONTINUE
            CALL INVERT(NDIM,NMAX,Q1,P1,COND,IPVT,WORK,B)
            CALL INVERT(NDIM,NMAX,Q2,P2,COND,IPVT,WORK,B)
         DO 30 I=1,NNMAX
             DO 30 J=1,NNMAX
                A(J,I)=0._dp
       30 CONTINUE
        DO 40 I=1,NMAX
              I1=IND1(I)
              I2=IND2(I)
              DO 40 J=1,NMAX
                   J1=IND1(J)
                   J2=IND2(J)
                   A(J1,I1)=P1(J,I)
                   A(J2,I2)=P2(J,I)
      40 CONTINUE
     
    END  SUBROUTINE
    
    SUBROUTINE INVERT (NDIM,N,A,X,COND,IPVT,WORK,B)
!DIR4 ATTRIBUTES CODE_ALIGN : 32 :: INVERT
          integer(kind=int4) :: NDIM,N
          real(kind=dp), dimension(NDIM,N) :: A,X
          real(kind=dp) :: COND
          integer(kind=int4), dimension(N) ::  IPVT
          real(kind=dp), dimension(N) :: WORK,B
          ! Locals
          integer(kind=int4) :: I,J
          ! Exec code ...
          CALL DECOMP (NDIM,N,A,COND,IPVT,WORK)
          IF (COND+1.0_dp.EQ.COND) PRINT 5,COND
!C     IF (COND+1D0.EQ.COND) STOP
    5 FORMAT(' THE MATRIX IS SINGULAR FOR THE GIVEN NUMERICAL ACCURACY '  , &
                'COND = ',D12.6)
      DO 30 I=1,N
           DO 10 J=1,N
                B(J)=0._dp
                IF (J.EQ.I) B(J)=1._dp
  10       CONTINUE
           CALL SOLVE (NDIM,N,A,B,IPVT)
           DO 30 J=1,N
                X(J,I)=B(J)
   30 CONTINUE
     
    END SUBROUTINE
    
    SUBROUTINE DECOMP (NDIM,N,A,COND,IPVT,WORK)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DECOMP
          integer(kind=int4) :: NDIM,N
          real(kind=dp), dimension(NDIM,N) :: A
          real(kind=dp) :: COND
          integer(kind=int4), dimension(N) ::  IPVT
          real(kind=dp), dimension(N) :: WORK
          ! Locals
          integer(kind=int4) :: NM1,J,I,K,KP1,M,KB,KM1,KP
          real(kind=dp) :: ANORM,T,EK,YNORM,ZNORM
          ! Exec code ...
          IPVT(N)=1
          IF(N.EQ.1) GO TO 80
          NM1=N-1
          ANORM=0._dp
          DO 10 J=1,N
             T=0._dp
!DIR$ SIMD
              DO 5 I=1,N
                 T=T+DABS(A(I,J))
    5     CONTINUE
          IF (T.GT.ANORM) ANORM=T
   10 CONTINUE
      DO 35 K=1,NM1
          KP1=K+1
          M=K
          DO 15 I=KP1,N
              IF (DABS(A(I,K)).GT.DABS(A(M,K))) M=I
   15     CONTINUE
          IPVT(K)=M
          IF (M.NE.K) IPVT(N)=-IPVT(N)
          T=A(M,K)
          A(M,K)=A(K,K)
          A(K,K)=T
          IF (T.EQ.0._dp) GO TO 35
          DO 20 I=KP1,N
              A(I,K)=-A(I,K)/T
   20     CONTINUE
          DO 30 J=KP1,N
              T=A(M,J)
              A(M,J)=A(K,J)
              A(K,J)=T
              IF (T.EQ.0._dp) GO TO 30
              DO 25 I=KP1,N
                  A(I,J)=A(I,J)+A(I,K)*T
   25         CONTINUE
   30     CONTINUE
   35 CONTINUE
      DO 50 K=1,N
          T=0._dp
          IF (K.EQ.1) GO TO 45
          KM1=K-1
          DO 40 I=1,KM1
              T=T+A(I,K)*WORK(I)
   40     CONTINUE
   45     EK=1._dp
          IF (T.LT.0.0_dp) EK=-1._dp
          IF (A(K,K).EQ.0._dp) GO TO 90
          WORK(K)=-(EK+T)/A(K,K)
   50 CONTINUE
      DO 60 KB=1,NM1
          K=N-KB
          T=0._dp
          KP1=K+1
          DO 55 I=KP1,N
              T=T+A(I,K)*WORK(K)
   55     CONTINUE
          WORK(K)=T
          M=IPVT(K)
          IF (M.EQ.K) GO TO 60
          T=WORK(M)
          WORK(M)=WORK(K)
          WORK(K)=T
   60 CONTINUE
      YNORM=0._dp
      DO 65 I=1,N
          YNORM=YNORM+DABS(WORK(I))
   65 CONTINUE
      CALL SOLVE (NDIM,N,A,WORK,IPVT)
      ZNORM=0._dp
      DO 70 I=1,N
          ZNORM=ZNORM+DABS(WORK(I))
   70 CONTINUE
      COND=ANORM*ZNORM/YNORM
      IF (COND.LT.1.0_dp) COND=1._dp
      RETURN
   80 COND=1._dp
      IF (A(1,1).NE.0._dp) RETURN
   90 COND=HUGE(1._dp)
    
    END  SUBROUTINE
    
    SUBROUTINE SOLVE (NDIM,N,A,B,IPVT)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: SOLVE
          integer(kind=int4) :: NDIM,N
          real(kind=dp), dimension(NDIM,N) :: A
          real(kind=dp), dimension(N)      :: B
          integer(kind=int4), dimension(N) :: IPVT
          ! Locals
          integer(kind=int4) :: NM1,K,KP1,I,KB,KM1,M
          real(kind=dp) :: T
          IF (N.EQ.1) GO TO 50
             NM1=N-1
          DO 20 K=1,NM1
                KP1=K+1
                M=IPVT(K)
                T=B(M)
                B(M)=B(K)
                B(K)=T
            !DIR4 SIMD
            DO 10 I=KP1,N
              B(I)=B(I)+A(I,K)*T
   10     CONTINUE
   20 CONTINUE
      DO 40 KB=1,NM1
          KM1=N-KB
          K=KM1+1
          B(K)=B(K)/A(K,K)
          T=-B(K)
          !DIR$ SIMD
          DO 30 I=1,KM1
              B(I)=B(I)+A(I,K)*T
   30     CONTINUE
   40 CONTINUE
   50 B(1)=B(1)/A(1,1)
     
    END  SUBROUTINE
!DIR$ ATTRIBUTES INLINE :: SAREA    
    SUBROUTINE SAREA (D,RAT)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: SAREA
          real(kind=dp) :: D,RAT
          ! Locals
          real(kind=dp) :: E,R
          IF (D.GE.1) GO TO 10
          E=DSQRT(1._dp-D*D)
          R=0.5_dp*(D**(2._dp/3._dp) + D**(-0.33333333333333333_dp)*DASIN(E)/E)
          R=DSQRT(R)
          RAT=1._dp/R
          RETURN
       10 E=DSQRT(1._dp-1._dp/(D*D))
          R=0.25_dp*(2._dp*D**(2._dp/3._dp) + D**(-4._dp/3._dp)*DLOG((1._dp)/(0.1_dp))/E)
    
          R=DSQRT(R)
          RAT=1._dp/R
      
    END  SUBROUTINE
    
    SUBROUTINE SURFCH (N,E,RAT)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: SURFCH
          integer(kind=int4) :: N
          real(kind=dp) :: E,RAT
          ! Locals
          real(kind=dp), dimension(60) :: X,W
          integer(kind=int4) :: I,NG
          real(kind=dp) :: DN,E2,EN,S,V,XI,DX,DXI,DS,DSN,DCN,A2, &
                        A,ENS,RS,RV,DXN
          DN=DFLOAT(N)
          E2=E*E
          EN=E*DN
          NG=60
          CALL GAUSS (NG,0,0,X,W)
          S=0._dp
          V=0._dp
          DO 10 I=1,NG
                XI=X(I)
                DX=DACOS(XI)
                DXN=DN*DX
                DS=DSIN(DX)
                DSN=DSIN(DXN)
                DCN=DCOS(DXN)
                A=1._dp+E*DCN
                A2=A*A
                ENS=EN*DSN
                S=S+W(I)*A*DSQRT(A2+ENS*ENS)
                V=V+W(I)*(DS*A+XI*ENS)*DS*A2
      10 CONTINUE
        RS=DSQRT(S*0.5_dp)
        RV=(V*3._dp/4._dp)**(0.3333333333333333333333_dp)
        RAT=RV/RS
    
    END  SUBROUTINE
!DIR$ ATTRIBUTES INLINE :: SAREAC    
    SUBROUTINE SAREAC (EPS,RAT)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: SAREAC
        real(kind=dp) :: EPS,RAT
        RAT=(1.5_dp/EPS)**(0.33333333333333333333333_dp)
        RAT=RAT/DSQRT( (EPS+2._dp)/(2._dp*EPS) )
      
    END SUBROUTINE
    
    SUBROUTINE DROP (RAT)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: DROP
          real(kind=dp) :: RAT
          ! Locls
          integer(kind=int4), parameter :: NC = 10, NG = 60
          real(kind=dp), dimension(NG) :: X,W
          real(kind=dp), dimension(0:NC) :: C
          real(kind=dp) :: V,S,XI,WI,RI,DRI,XIN,CI,RISI,RS,RV,R0V,SI
          integer(kind=int4) :: I,N
          COMMON /CDROP/ C,R0V
          C(0)=-0.0481_dp
          C(1)= 0.0359_dp
          C(2)=-0.1263_dp
          C(3)= 0.0244_dp
          C(4)= 0.0091_dp
          C(5)=-0.0099_dp
          C(6)= 0.0015_dp
          C(7)= 0.0025_dp
          C(8)=-0.0016_dp
          C(9)=-0.0002_dp
          C(10)= 0.0010_dp
          CALL GAUSS (NG,0,0,X,W)
          S=0._dp
          V=0._dp
          DO I=1,NG
               XI=DACOS(X(I))
               WI=W(I)
               RI=1._dp+C(0)
               DRI=0._dp
               DO N=1,NC
                    XIN=XI*N
                    RI=RI+C(N)*DCOS(XIN)
                    DRI=DRI-C(N)*N*DSIN(XIN)
               ENDDO
               SI=DSIN(XI)
               CI=X(I)
               RISI=RI*SI
               S=S+WI*RI*DSQRT(RI*RI+DRI*DRI)
               V=V+WI*RI*RISI*(RISI-DRI*CI)
         ENDDO
         RS=DSQRT(S*0.5_dp)
         RV=(V*3._dp*0.25_dp)**(0.33333333333333333333_dp)
         IF (DABS(RAT-1._dp).GT.0.00000001_dp) RAT=RV/RS
         R0V=1._dp/RV
         WRITE (6,1000) R0V
         DO N=0,NC
             WRITE (6,1001) N,C(N)
         ENDDO
 1000 FORMAT ('r_0/r_ev=',F7.4)
 1001 FORMAT ('c_',I2,'=',F7.4)
    
    END  SUBROUTINE
    

!C**********************************************************************
!C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         *
!C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      *
!C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.               *
!C    N - NUMBER OF POINTS                                             *
!C    Z - DIVISION POINTS                                              *
!C    W - WEIGHTS                                                      *
!C**********************************************************************
 
    SUBROUTINE GAUSS (N,IND1,IND2,Z,W)
!DIR$ ATTRIBUTES CODE_ALIGN : 32 :: GAUSS
          integer(kind=int4) :: N,IND1,IND2
          real(kind=dp), dimension(N) :: Z,W
          real(kind=dp), parameter :: A = 1._dp, B = 2._dp, C = 3._dp
        !  DATA A,B,C /1._R64P,2D0,3D0/
          integer(kind=int4) :: IND,K,I,M,NITER,J
          real(kind=dp) :: F,CHECK,PB,PC,X,DJ,PA,ZZ
          IND=MOD(N,2)
          K=N/2+IND
          F=DFLOAT(N)
          DO 100 I=1,K
                 M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4.0_dp+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6_dp+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0.0_dp
          NITER=0
          CHECK=0.0000000000000001_dp
   10     PB=1._dp
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10._dp
   15     PC=X
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
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',   &
       ' OF ',I4,'-TH ORDER')
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
    
    END  SUBROUTINE
    
    

    
    
    
    
      
      
    
    
    
    
    
    
    
    
    
    
    
    
    
                     
                     
     
    
    
      
     

    
    
      
      
     
          
      
    
     
      
     
          
     
    
    


end module mod_tmatrix_mps
