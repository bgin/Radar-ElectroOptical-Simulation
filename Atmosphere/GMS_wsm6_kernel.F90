
#ifdef IINSIDE
#define _NISLFV_RAIN_PLM_ nislfv_rain_plm_ii
#define _NISLFV_RAIN_PLM6_ nislfv_rain_plm6_ii
#else
#define _NISLFV_RAIN_PLM_ nislfv_rain_plm
#define _NISLFV_RAIN_PLM6_ nislfv_rain_plm6
#endif

MODULE module_mp_wsm6
  !
   use mod_kinds, int4,sp
   IMPLICIT NONE
!
   REAL(kind=sp), PARAMETER, PUBLIC :: dtcldcr     = 120._sp ! maximum time step for minor loops
   REAL(kind=sp), PARAMETER, PUBLIC :: n0r = 8.e6_sp         ! intercept parameter rain
   REAL(kind=sp), PARAMETER, PUBLIC :: n0g = 4.e6_sp         ! intercept parameter graupel
   REAL(kind=sp), PARAMETER, PUBLIC :: avtr = 841.9_sp       ! a constant for terminal velocity of rain
   REAL(kind=sp), PARAMETER, PUBLIC :: bvtr = 0.8_sp         ! a constant for terminal velocity of rain
   REAL(kind=sp), PARAMETER, PUBLIC :: r0 = .8e-5_sp         ! 8 microm  in contrast to 10 micro m
   REAL(kind=sp), PARAMETER, PUBLIC :: peaut = .55_sp        ! collection efficiency
   REAL(kind=sp), PARAMETER, PUBLIC :: xncr = 3.e8_sp        ! maritime cloud in contrast to 3.e8 in tc80
   REAL(kind=sp), PARAMETER, PUBLIC :: xmyu = 1.718e-5_sp    ! the dynamic viscosity kgm-1s-1
   REAL(kind=sp), PARAMETER, PUBLIC :: avts = 11.72_sp       ! a constant for terminal velocity of snow
   REAL(kind=sp), PARAMETER, PUBLIC :: bvts = .41_sp         ! a constant for terminal velocity of snow
   REAL(kind=sp), PARAMETER, PUBLIC :: avtg = 330._sp        ! a constant for terminal velocity of graupel
   REAL(kind=sp), PARAMETER, PUBLIC :: bvtg = 0.8_sp         ! a constant for terminal velocity of graupel
   REAL(kind=sp), PARAMETER, PUBLIC :: deng = 500._sp        ! density of graupel
   REAL(kind=sp), PARAMETER, PUBLIC :: n0smax =  1.e11_sp    ! maximum n0s (t=-90C unlimited)
   REAL(kind=sp), PARAMETER, PUBLIC :: lamdarmax = 8.e4_sp   ! limited maximum value for slope parameter of rain
   REAL(kind=sp), PARAMETER, PUBLIC :: lamdasmax = 1.e5_sp   ! limited maximum value for slope parameter of snow
   REAL(kind=sp), PARAMETER, PUBLIC :: lamdagmax = 6.e4_sp   ! limited maximum value for slope parameter of graupel
   REAL(kind=sp), PARAMETER, PUBLIC :: dicon = 11.9_sp       ! constant for the cloud-ice diamter
   REAL(kind=sp), PARAMETER, PUBLIC :: dimax = 500.e-6_sp    ! limited maximum value for the cloud-ice diamter
   REAL(kind=sp), PARAMETER, PUBLIC :: n0s = 2.e6_sp         ! temperature dependent intercept parameter snow
   REAL(kind=sp), PARAMETER, PUBLIC :: alpha = .12_sp        ! .122 exponen factor for n0s
   REAL(kind=sp), PARAMETER, PUBLIC :: pfrz1 = 100._sp       ! constant in Biggs freezing
   REAL(kind=sp), PARAMETER, PUBLIC :: pfrz2 = 0.66_sp       ! constant in Biggs freezing
   REAL(kind=sp), PARAMETER, PUBLIC :: qcrmin = 1.e-9_sp     ! minimun values for qr, qs, and qg
   REAL(kind=sp), PARAMETER, PUBLIC :: eacrc = 1.0_sp        ! Snow/cloud-water collection efficiency
   REAL(kind=sp), PARAMETER, PUBLIC :: dens  =  100.0_sp     ! Density of snow
   REAL(kind=sp), PARAMETER, PUBLIC :: qs0   =  6.e-4_sp     ! threshold amount for aggretion to occur
   REAL(kind=sp), SAVE ::                                      &
             qc0, qck1,bvtr1,bvtr2,bvtr3,bvtr4,g1pbr, &
             g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr,    &
             bvtr6,g6pbr,                             &
             precr1,precr2,roqimax,bvts1,             &
             bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs,     &
             g5pbso2,pvts,pacrs,precs1,precs2,pidn0r, &
             pidn0s,xlv1,pacrc,pi,                    &
             bvtg1,bvtg2,bvtg3,bvtg4,g1pbg,           &
             g3pbg,g4pbg,g5pbgo2,pvtg,pacrg,          &
             precg1,precg2,pidn0g,                    &
             rslopermax,rslopesmax,rslopegmax,        &
             rsloperbmax,rslopesbmax,rslopegbmax,     &
             rsloper2max,rslopes2max,rslopeg2max,     &
             rsloper3max,rslopes3max,rslopeg3max
CONTAINS
!===================================================================
!
!===================================================================
  !
#if defined __GFORTRAN__
  SUBROUTINE wsm62D(t, q                                          &   
                   ,qci, qrs, den, p, delz                        &
                   ,delt,g, cpd, cpv, rd, rv, t0c                 &
                   ,ep1, ep2, qmin                                &
                   ,XLS, XLV0, XLF0, den0, denr                   &
                   ,cliq,cice,psat                                &
                   ,lat                                           &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte                     &
                   ,snow,snowncv                                  &
                   ,graupel,graupelncv                            &
                   )    !GCC$ ATTRIBUTES hot :: wsm62D   !GCC$ ATTRIBUTES aligned(32) :: wsm62D
#elif defined __INTEL_COMPILER
  SUBROUTINE wsm62D(t, q                                          &   
                   ,qci, qrs, den, p, delz                        &
                   ,delt,g, cpd, cpv, rd, rv, t0c                 &
                   ,ep1, ep2, qmin                                &
                   ,XLS, XLV0, XLF0, den0, denr                   &
                   ,cliq,cice,psat                                &
                   ,lat                                           &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte                     &
                   ,snow,snowncv                                  &
                   ,graupel,graupelncv                            &
                   )
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: wsm62D
#endif
!-------------------------------------------------------------------
 ! IMPLICIT NONE
!-------------------------------------------------------------------
!
!  This code is a 6-class GRAUPEL phase microphyiscs scheme (WSM6) of the 
!  Single-Moment MicroPhyiscs (WSMMP). The WSMMP assumes that ice nuclei
!  number concentration is a function of temperature, and seperate assumption
!  is developed, in which ice crystal number concentration is a function
!  of ice amount. A theoretical background of the ice-microphysics and related
!  processes in the WSMMPs are described in Hong et al. (2004).
!  All production terms in the WSM6 scheme are described in Hong and Lim (2006).
!  All units are in m.k.s. and source/sink terms in kgkg-1s-1.
!
!  WSM6 cloud scheme
!
!  Coded by Song-You Hong and Jeong-Ock Jade Lim (Yonsei Univ.)
!           Summer 2003
!
!  Implemented by Song-You Hong (Yonsei Univ.) and Jimy Dudhia (NCAR)
!           Summer 2004
!
!  History :  semi-lagrangian scheme sedimentation(JH), and clean up
!             Hong, August 2009
!
!  Reference) Hong, Dudhia, Chen (HDC, 2004) Mon. Wea. Rev.
!             Hong and Lim (HL, 2006) J. Korean Meteor. Soc.
!             Dudhia, Hong and Lim (DHL, 2008) J. Meteor. Soc. Japan
!             Lin, Farley, Orville (LFO, 1983) J. Appl. Meteor.
!             Rutledge, Hobbs (RH83, 1983) J. Atmos. Sci.
!             Rutledge, Hobbs (RH84, 1984) J. Atmos. Sci.
!             Juang and Hong (JH, 2010) Mon. Wea. Rev.
!
  INTEGER(kind=int4),      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte,  &
                                      lat
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ),                           &
        INTENT(INOUT) ::                                          &
                                                               t
  REAL(kind=sp), DIMENSION( its:ite , kts:kte, 2 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qci
  REAL(kind=sp), DIMENSION( its:ite , kts:kte, 3 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qrs
  REAL(kind=sp), DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                               q
  REAL(kind=sp), DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(IN   ) ::                                          &
                                                             den, &
                                                               p, &
                                                            delz
  REAL(kind=sp), INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                             cpd, &
                                                             cpv, &
                                                             t0c, &
                                                            den0, &
                                                              rd, &
                                                              rv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL(kind=sp), DIMENSION( ims:ime ),                                     &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr
  REAL(kind=sp), DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv
  REAL(kind=sp), DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                 graupel, &
                                                      graupelncv
  ! LOCAL VAR
#if defined __INTEL_COMPILER
  REAL(kind=sp), DIMENSION( its:ite , kts:kte , 3) ::             &
                                                              rh, &
                                                              qs, &
                                                          rslope, &
                                                         rslope2, &
                                                         rslope3, &
                                                         rslopeb, &
                                                         qrs_tmp, & 
                                                            falk, &
                                                            fall, &
                                                            work1
  !DIR$ ATTRIBUTES ALIGN : 64 :: rh,qs,rslope,rslope2
  !DIR$ ATTRIBUTES ALIGN : 64 :: rslope3,rslopeb,qrs_tmp
  !DIR$ ATTRIBUTES ALIGN : 64 :: falk,fall,work1
#elif defined __GFORTRAN__
  REAL(kind=sp), DIMENSION( its:ite, kts:kte, 3) :: rh        !GCC$ ATTRIBUTES aligned(64) :: rh
  REAL(kind=sp), DIMENSION( its:ite, kts:kte, 3) :: qs        !GCC$ ATTRIBUTES aligned(64) :: qs
  REAL(kind=sp), DIMENSION( its:ite, kts:kte, 3) :: rslope    !GCC$ ATTRIBUTES aligned(64) :: rslope
  REAL(kind=sp), DIMENSION( its:ite, kts:kte, 3) :: rslope2   !GCC$ ATTRIBUTES aligned(64) :: rslope2
  REAL(kind=sp), DIMENSION( its:ite, kts:kte, 3) :: rslope3   !GCC$ ATTRIBUTES aligned(64) :: rslope3
  REAL(kind=sp), DIMENSION( its:ite, kts:kte, 3) :: rslopeb   !GCC$ ATTRIBUTES aligned(64) :: rslopeb
  REAL(kind=sp), DIMENSION( its:ite, kts:kte, 3) :: qrs_tmp   !GCC$ ATTRIBUTES aligned(64) :: qrs_tmp
  REAL(kind=sp), DIMENSION( its:ite, kts:kte, 3) :: falk      !GCC$ ATTRIBUTES aligned(64) :: falk
  REAL(kind=sp), DIMENSION( its:ite, kts:kte, 3) :: fall      !GCC$ ATTRIBUTES aligned(64) :: fall
  REAL(kind=sp), DIMENSION( its:ite, kts:kte, 3) :: work1     !GCC$ ATTRIBUTES aligned(64) :: work1
#endif
#if defined __INTEL_COMPILER
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) ::                &
                                                           fallc, &
                                                           falkc, &
                                                          work1c, &
                                                          work2c, &
                                                           workr, &
                                                           worka
  !DIR$ ATTRIBUTES ALIGN : 64 :: fallc,falkc,work1c
  !DIR$ ATTRIBUTES ALIGN : 64 :: work2c,workr,worka
#elif defined __GFORTRAN__
  REAL(kind=sp), DIMENSION( its:ite, kts:kte) :: fallc       !GCC$ ATTRIBUTES aligned(64) :: fallc
  REAL(kind=sp), DIMENSION( its:ite, kts:kte) :: falkc       !GCC$ ATTRIBUTES aligned(64) :: falkc
  REAL(kind=sp), DIMENSION( its:ite, kts:kte) :: work1c      !GCC$ ATTRIBUTES aligned(64) :: work1c
  REAL(kind=sp), DIMENSION( its:ite, kts:kte) :: work2c      !GCC$ ATTRIBUTES aligned(64) :: work2c
  REAL(kind=sp), DIMENSION( its:ite, kts:kte) :: workr       !GCC$ ATTRIBUTES aligned(64) :: workr
  REAL(kind=sp), DIMENSION( its:ite, kts:kte) :: worka       !GCC$ ATTRIBUTES aligned(64) :: worka
#endif
#if defined __INTEL_COMPILER
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) ::                         &
                                                         den_tmp, &
                                                         delz_tmp
  !DIR$ ATTRIBUTES ALIGN : 64 :: den_tmp,delz_tmp
#elif defined __GFORTRAN__
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: den_tmp   !GCC$ ATTRIBUTES aligned(64) :: den_tmp
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: delz_tmp  !GCC$ ATTRIBUTES aligned(64) :: delz_tmp
#endif
#if defined __INTEL_COMPILER
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) ::                &
                                                           pigen, &
                                                           pidep, &
                                                           pcond, &
                                                           prevp, &
                                                           psevp, &
                                                           pgevp, &
                                                           psdep, &
                                                           pgdep, &
                                                           praut, &
                                                           psaut, &
                                                           pgaut, &
                                                           piacr, &
                                                           pracw, &
                                                           praci, &
                                                           pracs, &
                                                           psacw, &
                                                           psaci, &
                                                           psacr, &
                                                           pgacw, &
                                                           pgaci, &
                                                           pgacr, &
                                                           pgacs, &
                                                           paacw, &
                                                           psmlt, &
                                                           pgmlt, &
                                                           pseml, &
                                                           pgeml
  !DIR$ ATTRIBUTES ALIGN : 64 :: pigen,pidep,pcond
  !DIR$ ATTRIBUTES ALIGN : 64 :: prevp,psevp,pgevp
  !DIR$ ATTRIBUTES ALIGN : 64 :: praut,psaut,pgaut
  !DIR$ ATTRIBUTES ALIGN : 64 :: piacr,pracw,praci
  !DIR$ ATTRIBUTES ALIGN : 64 :: pracs,psacw,psaci
  !DIR$ ATTRIBUTES ALIGN : 64 :: psacr,pgacw,pgaci
  !DIR$ ATTRIBUTES ALIGN : 64 :: pgacr,pgacs,paacw
  !DIR$ ATTRIBUTES ALIGN : 64 :: psmlt,pgmlt,pseml
  !DIR$ ATTRIBUTES ALIGN : 64 :: pgeml
#elif defined __GFORTRAN__
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pigen    !GCC$ ATTRIBUTES aligned(64) :: pigen
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pidep    !GCC$ ATTRIBUTES aligned(64) :: pidep
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pcond    !GCC$ ATTRIBUTES aligned(64) :: pcond
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: prevp    !GCC$ ATTRIBUTES aligned(64) :: prevp
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: psevp    !GCC$ ATTRIBUTES aligned(64) :: psevp
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pgevp    !GCC$ ATTRIBUTES aligned(64) :: pgevp
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: psdep    !GCC$ ATTRIBUTES aligned(64) :: psdep
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pgdep    !GCC$ ATTRIBUTES aligned(64) :: pgdep
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: praut    !GCC$ ATTRIBUTES aligned(64) :: praut
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: psaut    !GCC$ ATTRIBUTES aligned(64) :: psaut
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pgaut    !GCC$ ATTRIBUTES aligned(64) :: pgaut
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: piacr    !GCC$ ATTRIBUTES aligned(64) :: piacr
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pracw    !GCC$ ATTRIBUTES aligned(64) :: pracw
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: psacw    !GCC$ ATTRIBUTES aligned(64) :: psacw
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: psaci    !GCC$ ATTRIBUTES aligned(64) :: psaci
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: psacr    !GCC$ ATTRIBUTES aligned(64) :: psacr
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pgacw    !GCC$ ATTRIBUTES aligned(64) :: pgacw
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pgaci    !GCC$ ATTRIBUTES aligned(64) :: pgaci
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pgacr    !GCC$ ATTRIBUTES aligned(64) :: pgacr
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pgacs    !GCC$ ATTRIBUTES aligned(64) :: pgacs
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: paacw    !GCC$ ATTRIBUTES aligned(64) :: paacw
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: psmlt    !GCC$ ATTRIBUTES aligned(64) :: psmlt
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pgmlt    !GCC$ ATTRIBUTES aligned(64) :: pgmlt
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pseml    !GCC$ ATTRIBUTES aligned(64) :: pseml
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: pgeml    !GCC$ ATTRIBUTES aligned(64) :: pgeml
#endif
#if defined __INTEL_COMPILER
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) ::                         &
                                                            qsum, &
                                                              xl, &
                                                             cpm, &
                                                           work2, &
                                                          denfac, &
                                                             xni, &
                                                         denqrs1, &
                                                         denqrs2, &
                                                         denqrs3, &
                                                          denqci, & 
                                                          delta2, &
                                                          delta3, &
                                                          n0sfac
  !DIR$ ATTRIBUTES ALIGN : 64 :: qsum,xl,cpm,work2,denfac
  !DIR$ ATTRIBUTES ALIGN : 64 :: xni,denqrs1,denqrs2,denqrs3
  !DIR$ ATTRIBUTES ALIGN : 64 :: denqci,delta2,delta3,n0sfac
#elif defined __GFORTRAN__
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: qsum    !GCC$ ATTRIBUTES aligned(64) :: qsum
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: xl      !GCC$ ATTRIBUTES aligned(64) :: xl
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: cpm     !GCC$ ATTRIBUTES aligned(64) :: cpm
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: work2   !GCC$ ATTRIBUTES aligned(64) :: work2
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: denfac  !GCC$ ATTRIBUTES aligned(64) :: xni
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: denqrs1 !GCC$ ATTRIBUTES aligned(64) :: denqrs1
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: denqrs2 !GCC$ ATTRIBUTES aligned(64) :: denqrs2
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: denqrs3 !GCC$ ATTRIBUTES aligned(64) :: denqrs3
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: denqci  !GCC$ ATTRIBUTES aligned(64) :: denqci
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: delta2  !GCC$ ATTRIBUTES aligned(64) :: delta2
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: delta3  !GCC$ ATTRIBUTES aligned(64) :: delta3
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) :: n0sfac  !GCC$ ATTRIBUTES aligned(64) :: n0sfac
#endif
#if defined __INTEL_COMPILER
  REAL(kind=sp), DIMENSION( its:ite ) ::                 delqrs1, &
                                                         delqrs2, &
                                                         delqrs3, &
                                                         delqi
  !DIR$ ATTRIBUTES ALIGN : 64 :: delqrs1,delqrs2,delqrs3,delqi
#elif defined __GFORTRAN__
  REAL(kind=sp), DIMENSION( its:ite ) ::  delqrs1   !GCC$ ATTRIBUTES aligned(64) :: delqrs1
  REAL(kind=sp), DIMENSION( its:ite ) ::  delqrs2   !GCC$ ATTRIBUTES aligned(64) :: delqrs2
  REAL(kind=sp), DIMENSION( its:ite ) ::  delqrs3   !GCC$ ATTRIBUTES aligned(64) :: delqrs3
  REAL(kind=sp), DIMENSION( its:ite ) ::  delqi     !GCC$ ATTRIBUTES aligned(64) :: delqi
#endif
#if defined __INTEL_COMPILER 
  REAL(kind=sp), DIMENSION( its:ite ) ::              tstepsnow, &
                                                      tstepgraup
  !DIR$ ATTRIBUTES ALIGN : 64 :: tstepsnow,tstepgraup
#elif defined __GFORTRAN__
  REAL(kind=sp), DIMENSION( its:ite ) :: tstepsnow  !GCC$ ATTRIBUTES aligned(64) :: tstepsnow
  REAL(kind=sp), DIMENSION( its:ite ) :: tstepgraup !GCC$ ATTRIBUTES aligned(64) :: tstepgraup
#endif
#if defined __INTEL_COMPILER
  INTEGER(kind=int4), DIMENSION( its:ite ) ::              mstep, &
                                                           numdt
  !DIR$ ATTRIBUTES ALIGN : 64 :: mstep,numdt
#elif defined __GFORTRAN__
  INTEGER(kind=int4), DIMENSION( its:ite ) :: mstep  !GCC$ ATTRIBUTES aligned(64) :: mstep
  INTEGER(kind=int4), DIMENSION( its:ite ) :: numdt  !GCC$ ATTRIBUTES aligned(64) :: numdt
#endif
#if defined __INTEL_COMPILER
  LOGICAL(kind=int4), DIMENSION( its:ite ) :: flgcld 
#endif
  !DIR$ ATTRIBUTES ALIGN : 64 :: flgcld
  REAL(kind=sp)  ::                                                        &
            cpmcal, xlcal, diffus,                                &
            viscos, xka, venfac, conden, diffac,                  &
            x, y, z, a, b, c, d, e,                               &
            qdt, holdrr, holdrs, holdrg, supcol, supcolt, pvt,    &
            coeres, supsat, dtcld, xmi, eacrs, satdt,             &
            qimax, diameter, xni0, roqi0,                         &
            fallsum, fallsum_qsi, fallsum_qg,                     &
            vt2i,vt2r,vt2s,vt2g,acrfac,egs,egi,                   &
            xlwork2, factor, source, value,                       &
            xlf, pfrzdtc, pfrzdtr, supice, alpha2
  REAL(kind=sp)  :: vt2ave
  REAL(kind=sp)  :: holdc, holdci
  INTEGER(kind=int4) :: i, j, k, mstepmax,                                   &
            iprt, latd, lond, loop, loops, ifsat, n, idim, kdim
  INTEGER(kind=int4) :: itest,ktest
! Temporaries used for inlining fpvs function
  REAL(kind=sp)  :: dldti, xb, xai, tr, xbi, xa, hvap, cvap, hsub, dldt, ttp
  ! variables for optimization
#if defined __INTEL_COMPILER
  REAL(kind=sp), DIMENSION( its:ite ) ::                             tvec1
  !DIR$ ATTRIBUTES ALIGNED : 64 :: tvec1
#elif defined __GFORTRAN__
  REAL(kind=sp), DIMENSION( its:ite ) ::                             tvec1  !GCC$ ATTRIBUTES aligned(64) :: tvec1
#endif
  REAL(kind=sp)                       ::                              temp
!
!=================================================================
!   compute internal functions
!
      cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
      xlcal(x) = xlv0-xlv1*(x-t0c)
!----------------------------------------------------------------
!     diffus: diffusion coefficient of the water vapor
!     viscos: kinematic viscosity(m2s-1)
!     Optimizatin : A**B => exp(log(A)*(B))
!
      diffus(x,y) = 8.794e-5 * exp(log(x)*(1.81)) / y        ! 8.794e-5*x**1.81/y
      viscos(x,y) = 1.496e-6 * (x*sqrt(x)) /(x+120.)/y  ! 1.496e-6*x**1.5/(x+120.)/y
      xka(x,y) = 1.414e3*viscos(x,y)*y
      diffac(a,b,c,d,e) = d*a*a/(xka(c,d)*rv*c*c)+1./(e*diffus(c,b))
      venfac(a,b,c) = exp(log((viscos(b,c)/diffus(b,a)))*((.3333333)))         &
                     /sqrt(viscos(b,c))*sqrt(sqrt(den0/c))
      conden(a,b,c,d,e) = (max(b,qmin)-c)/(1.+d*d/(rv*e)*c/(a*a))
!

!!DIR$ ASSUME_ALIGNED t:64,qci:64,qrs:64,q:64,den:64,p:64,delz:64,rain:64,rainncv:64,sr:64
!!DIR$ ASSUME_ALIGNED snow:64,snowncv:64,graupel:64,graupelncv:64,mstep:64,flgcld:64
!!DIR$ ASSUME_ALIGNED worka:64,pgacs:64,n0sfac:64,work2:64,psmlt:64,rslope:64
!!DIR$ ASSUME_ALIGNED rslope2:64,rslope3:64,rslopeb:64,cpm:64,pgmlt:64
!!DIR$ ASSUME_ALIGNED tstepsnow:64,tstepgraup:64

!
      idim = ite-its+1
      kdim = kte-kts+1
itest=979
ktest=1
!
!----------------------------------------------------------------
!     padding 0 for negative values generated by dynamics
!
     do k = kts, kte
#if defined __INTEL_COMPILER
        !DIR$ VECTOR ALIGNED
        !DIR$ SIMD
        !DIR$ UNROLL(4)
#elif defined __GFORTRAN__
        !GCC$ VECTOR
        !GCC$ UNROLL 4
#endif
        do i = its, ite
#if defined __INTEL_COMPILER
           !DIR$ ASSUME_ALIGNED qci:64
           !DIR$ ASSUME_ALIGNED qrs:64
#endif         
          qci(i,k,1) = max(qci(i,k,1),0.0)
          qrs(i,k,1) = max(qrs(i,k,1),0.0)
          qci(i,k,2) = max(qci(i,k,2),0.0)
          qrs(i,k,2) = max(qrs(i,k,2),0.0)
          qrs(i,k,3) = max(qrs(i,k,3),0.0)
        enddo
      enddo
!
!----------------------------------------------------------------
!     latent heat for phase changes and heat capacity. neglect the
!     changes during microphysical process calculation
!     emanuel(1994)
!
      do k = kts, kte
        do i = its, ite
          cpm(i,k) = cpmcal(q(i,k))
          xl(i,k) = xlcal(t(i,k))
        enddo
      enddo
      do k = kts, kte
#if defined __INTEL_COMPILER
        !DIR$ VECTOR ALIGNED
        !DIR$ SIMD
        !DIR$ UNROLL(4)
#elif defined __GFORTRAN__
        !GCC$ VECTOR
        !GCC$ UNROLL 4
#endif         
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED delz:64,delz_tmp:64
            !DIR$ ASSUME_ALIGNED den:64,den_tmp:64
#endif
          delz_tmp(i,k) = delz(i,k)
          den_tmp(i,k) = den(i,k)
        enddo
      enddo
!
!----------------------------------------------------------------
!    initialize the surface rain, snow, graupel
!
      do i = its, ite
        rainncv(i) = 0.
        if(PRESENT (snowncv) .AND. PRESENT (snow)) snowncv(i,lat) = 0.
        if(PRESENT (graupelncv) .AND. PRESENT (graupel)) graupelncv(i,lat) = 0.
        sr(i) = 0.
! new local array to catch step snow and graupel
        tstepsnow(i) = 0.
        tstepgraup(i) = 0.
      enddo
!
!----------------------------------------------------------------
!     compute the minor time steps.
!
      loops = max(nint(delt/dtcldcr),1)
      dtcld = delt/loops
      if(delt.le.dtcldcr) dtcld = delt
!
      do loop = 1,loops
!
!----------------------------------------------------------------
!     initialize the large scale variables
!
      
      do i = its, ite
        mstep(i) = 1
        flgcld(i) = .true.
      enddo

      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
         !DIR$ UNROLL(7)
#elif defined __GFORTRAN__
         !GCC$ VECTOR
         !GCC$ UNROLL 7
#endif
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED denfac:64
#endif
          denfac(i,k) = sqrt(den0/den(i,k))
        enddo
      enddo
!
! Inline expansion for fpvs
!         qs(i,k,1) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!         qs(i,k,2) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ IVDEP
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#elif defined __GFORTRAN__
         !GCC$ IVDEP
         !GCC$ VECTOR
#endif
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED t:64
            !DIR$ ASSUME_ALIGNED p:64
            !DIR$ ASSUME_ALIGNED qs:64
            !DIR$ ASSUME_ALIGNED q:64
            !DIR$ ASSUME_ALIGNED rh:64
#endif
#if defined __GFORTRAN__
            !GCC$ builtin (exp) attributes simd
            !GCC$ builtin (log) attributes simd
#endif
          tr=ttp/t(i,k)
          qs(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          qs(i,k,1) = min(qs(i,k,1),0.99*p(i,k))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          rh(i,k,1) = max(q(i,k) / qs(i,k,1),qmin)
          tr=ttp/t(i,k)
          if(t(i,k).lt.ttp) then
            qs(i,k,2)=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
          else
            qs(i,k,2)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          endif
          qs(i,k,2) = min(qs(i,k,2),0.99*p(i,k))
          qs(i,k,2) = ep2 * qs(i,k,2) / (p(i,k) - qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
          rh(i,k,2) = max(q(i,k) / qs(i,k,2),qmin)
        enddo
      enddo
!
!----------------------------------------------------------------
!     initialize the variables for microphysical physics
!
!
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
         !DIR$ UNROLL(4)
#elif defined __GFORTRAN__
         !GCC$ VECTOR
         !GCC$ UNROLL 4
#endif
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED prevp:64
            !DIR$ ASSUME_ALIGNED psdep:64
            !DIR$ ASSUME_ALIGNED pgdep:64
            !DIR$ ASSUME_ALIGNED praut:64
            !DIR$ ASSUME_ALIGNED psaut:64
            !DIR$ ASSUME_ALIGNED pgaut:64
            !DIR$ ASSUME_ALIGNED pracw:64
            !DIR$ ASSUME_ALIGNED praci:64
            !DIR$ ASSUME_ALIGNED piacr:64
            !DIR$ ASSUME_ALIGNED psaci:64
            !DIR$ ASSUME_ALIGNED psacw:64
            !DIR$ ASSUME_ALIGNED pracs:64
            !DIR$ ASSUME_ALIGNED psacr:64
            !DIR$ ASSUME_ALIGNED pgacw:64
            !DIR$ ASSUME_ALIGNED paacw:64
            !DIR$ ASSUME_ALIGNED pgaci:64
            !DIR$ ASSUME_ALIGNED pgacr:64
            !DIR$ ASSUME_ALIGNED pgacs:64
            !DIR$ ASSUME_ALIGNED pigen:64
            !DIR$ ASSUME_ALIGNED pidep:64
            !DIR$ ASSUME_ALIGNED pcond:64
            !DIR$ ASSUME_ALIGNED psmlt:64
            !DIR$ ASSUME_ALIGNED pgmlt:64
            !DIR$ ASSUME_ALIGNED pseml:64
            !DIR$ ASSUME_ALIGNED pgeml:64
            !DIR$ ASSUME_ALIGNED psevp:64
            !DIR$ ASSUME_ALIGNED pgevp:64
            !DIR$ ASSUME_ALIGNED falk:64
            !DIR$ ASSUME_ALIGNED fall:64
            !DIR$ ASSUME_ALIGNED fallc:64
            !DIR$ ASSUME_ALIGNED falkc:64
#endif
          prevp(i,k) = 0.
          psdep(i,k) = 0.
          pgdep(i,k) = 0.
          praut(i,k) = 0.
          psaut(i,k) = 0.
          pgaut(i,k) = 0.
          pracw(i,k) = 0.
          praci(i,k) = 0.
          piacr(i,k) = 0.
          psaci(i,k) = 0.
          psacw(i,k) = 0.
          pracs(i,k) = 0.
          psacr(i,k) = 0.
          pgacw(i,k) = 0.
          paacw(i,k) = 0.
          pgaci(i,k) = 0.
          pgacr(i,k) = 0.
          pgacs(i,k) = 0.
          pigen(i,k) = 0.
          pidep(i,k) = 0.
          pcond(i,k) = 0.
          psmlt(i,k) = 0.
          pgmlt(i,k) = 0.
          pseml(i,k) = 0.
          pgeml(i,k) = 0.
          psevp(i,k) = 0.
          pgevp(i,k) = 0.
          falk(i,k,1) = 0.
          falk(i,k,2) = 0.
          falk(i,k,3) = 0.
          fall(i,k,1) = 0.
          fall(i,k,2) = 0.
          fall(i,k,3) = 0.
          fallc(i,k) = 0.
          falkc(i,k) = 0.
          xni(i,k) = 1.e3_sp
        enddo
      enddo
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(1)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(2)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(1)
#endif
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL 1
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL 2
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL 1
#endif
#endif
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED qci:64
            !DIR$ ASSUME_ALIGNED xni:64
#endif
          temp = (den(i,k)*max(qci(i,k,2),qmin))
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7_sp*temp,1.e3_sp),1.e6_sp)
        enddo
      enddo
!
!----------------------------------------------------------------
!     compute the fallout term:
!     first, vertical terminal velosity for minor loops
!----------------------------------------------------------------
      do k = kts, kte
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
        enddo
      enddo
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, & 
                     work1,its,ite,kts,kte)
!
      do k = kte, kts, -1
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
         !DIR$ UNROLL(4)
#elif defined __GFORTRAN__
         !GCC$ VECTOR
         !GCC$ UNROLL 4
#endif
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED workr:64
            !DIR$ ASSUME_ALIGNED work1:64
            !DIR$ ASSUME_ALIGNED worka:64
            !DIR$ ASSUME_ALIGNED qsum:64
#endif
          workr(i,k) = work1(i,k,1)
          qsum(i,k) = max( (qrs(i,k,2)+qrs(i,k,3)), 1.E-15_sp)
          IF ( qsum(i,k) .gt. 1.e-15_sp ) THEN
            worka(i,k) = (work1(i,k,2)*qrs(i,k,2) + work1(i,k,3)*qrs(i,k,3)) &
                      /qsum(i,k)
          ELSE
            worka(i,k) = 0.
          ENDIF
          denqrs1(i,k) = den(i,k)*qrs(i,k,1)
          denqrs2(i,k) = den(i,k)*qrs(i,k,2)
          denqrs3(i,k) = den(i,k)*qrs(i,k,3)
          if(qrs(i,k,1).le.0.0) workr(i,k) = 0.0_sp
        enddo
      enddo
      call _NISLFV_RAIN_PLM_(its,ite,kts,kte,den_tmp,denfac,t,delz_tmp,workr,denqrs1,  &
                           delqrs1,dtcld,1,1)
      call _NISLFV_RAIN_PLM6_(its,ite,kts,kte,den_tmp,denfac,t,delz_tmp,worka,         & 
                           denqrs2,denqrs3,delqrs2,delqrs3,dtcld,1,1)
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL 5
#endif
#endif
        do i = its, ite
          qrs(i,k,1) = max(denqrs1(i,k)/den(i,k),0.)
          qrs(i,k,2) = max(denqrs2(i,k)/den(i,k),0.)
          qrs(i,k,3) = max(denqrs3(i,k)/den(i,k),0.)
          fall(i,k,1) = denqrs1(i,k)*workr(i,k)/delz(i,k)
          fall(i,k,2) = denqrs2(i,k)*worka(i,k)/delz(i,k)
          fall(i,k,3) = denqrs3(i,k)*worka(i,k)/delz(i,k)
        enddo
     enddo
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL 5
#endif
#endif     
      do i = its, ite
        fall(i,1,1) = delqrs1(i)/delz(i,1)/dtcld
        fall(i,1,2) = delqrs2(i)/delz(i,1)/dtcld
        fall(i,1,3) = delqrs3(i)/delz(i,1)/dtcld
      enddo
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
         !DIR$ UNROLL(4)
#elif defined __GFORTRAN__
         !GCC$ VECTOR
         !GCC$ UNROLL 4
#endif
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
        enddo
      enddo
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                     work1,its,ite,kts,kte)
!
      do k = kte, kts, -1 
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ IVDEP
         !DIR$ VECTOR ALWAYS
         !DIR$ FMA
#elif defined __GFORTRAN__
         !GCC$ IVDEP
         !GCC$ VECTOR
#endif
         do i = its, ite
#if defined __GFORTRAN__
            !GCC$ builtin (exp) attributes simd
#endif
          supcol = t0c-t(i,k)
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(t(i,k).gt.t0c) then
!---------------------------------------------------------------
! psmlt: melting of snow [HL A33] [RH83 A25]
!       (T>T0: S->R)
!---------------------------------------------------------------
            xlf = xlf0
            work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
            if(qrs(i,k,2).gt.0.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psmlt(i,k) = xka(t(i,k),den(i,k))/xlf*(t0c-t(i,k))*pi/2.       &
                         *n0sfac(i,k)*(precs1*rslope2(i,k,2)                 &
                         +precs2*work2(i,k)*coeres)
              psmlt(i,k) = min(max(psmlt(i,k)*dtcld/mstep(i),                &
                          -qrs(i,k,2)/mstep(i)),0.)
              qrs(i,k,2) = qrs(i,k,2) + psmlt(i,k)
              qrs(i,k,1) = qrs(i,k,1) - psmlt(i,k)
              t(i,k) = t(i,k) + xlf/cpm(i,k)*psmlt(i,k)
            endif
!---------------------------------------------------------------
! pgmlt: melting of graupel [HL A23]  [LFO 47]
!       (T>T0: G->R)
!---------------------------------------------------------------
            if(qrs(i,k,3).gt.0.) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgmlt(i,k) = xka(t(i,k),den(i,k))/xlf                          &
                         *(t0c-t(i,k))*(precg1*rslope2(i,k,3)                &
                         +precg2*work2(i,k)*coeres)
              pgmlt(i,k) = min(max(pgmlt(i,k)*dtcld/mstep(i),                &
                          -qrs(i,k,3)/mstep(i)),0.)                          
              qrs(i,k,3) = qrs(i,k,3) + pgmlt(i,k)
              qrs(i,k,1) = qrs(i,k,1) - pgmlt(i,k)
              t(i,k) = t(i,k) + xlf/cpm(i,k)*pgmlt(i,k)
            endif
          endif
        enddo
      enddo
!---------------------------------------------------------------
! Vice [ms-1] : fallout of ice crystal [HDC 5a]
!---------------------------------------------------------------
      do k = kte, kts, -1
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !GCC$ VECTOR
         !GCC$ builtin  (exp) attributes simd
         !GCC$ builtin  (log) attributes simd
#endif
         do i = its, ite

          if(qci(i,k,2).le.0.) then
            work1c(i,k) = 0.
          else
            xmi = den(i,k)*qci(i,k,2)/xni(i,k)
            diameter  = max(min(dicon * sqrt(xmi),dimax), 1.e-25_sp)
            work1c(i,k) = 1.49e4_sp*exp(log(diameter)*(1.31_sp))
          endif
        enddo
      enddo
!
!  forward semi-laglangian scheme (JH), PCM (piecewise constant),  (linear)
!
      do k = kte, kts, -1
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(3)
#endif
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL 5
#endif
#endif
        do i = its, ite
          denqci(i,k) = den(i,k)*qci(i,k,2)
        enddo
      enddo
      call _NISLFV_RAIN_PLM_(its,ite,kts,kte,den_tmp,denfac,t,delz_tmp,work1c,denqci,  &
                           delqi,dtcld,1,0)
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL 5
#endif
#endif              
        do i = its, ite
          qci(i,k,2) = max(denqci(i,k)/den(i,k),0.)
        enddo
      enddo
      do i = its, ite
        fallc(i,1) = delqi(i)/delz(i,1)/dtcld
      enddo
!
!----------------------------------------------------------------

      !      rain (unit is mm/sec;kgm-2s-1: /1000*delt ===> m)==> mm for wrf
!
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ SIMD
#elif defined __GFORTRAN__
      !GCC$ VECTOR
#endif
      do i = its, ite
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED rainncv:64
         !DIR$ ASSUME_ALIGNED rain:64
         !DIR$ ASSUME_ALIGNED tstepsnow:64
         !DIR$ ASSUME_ALIGNED snowncv:64
         !DIR$ ASSUME_ALIGNED snow:64
         !DIR$ ASSUME_ALIGNED tstepgraup:64
         !DIR$ ASSUME_ALIGNED graupelncv:64
         !DIR$ ASSUME_ALIGNED graupel:64
#endif
        fallsum = fall(i,kts,1)+fall(i,kts,2)+fall(i,kts,3)+fallc(i,kts)
        fallsum_qsi = fall(i,kts,2)+fallc(i,kts)
        fallsum_qg = fall(i,kts,3)
        if(fallsum.gt.0.) then
          rainncv(i) = fallsum*delz(i,kts)/denr*dtcld*1000. + rainncv(i)
          rain(i) = fallsum*delz(i,kts)/denr*dtcld*1000. + rain(i)
        endif
        if(fallsum_qsi.gt.0.) then
          tstepsnow(i)   = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.            &
                           +tstepsnow(i)
        IF ( PRESENT (snowncv) .AND. PRESENT (snow)) THEN
          snowncv(i,lat) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.            & 
                           +snowncv(i,lat)
          snow(i,lat) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000. + snow(i,lat)
        ENDIF
        endif
        if(fallsum_qg.gt.0.) then
          tstepgraup(i)  = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.            &
                           +tstepgraup(i)
        IF ( PRESENT (graupelncv) .AND. PRESENT (graupel)) THEN
          graupelncv(i,lat) = fallsum_qg*delz(i,kts)/denr*dtcld*1000.          &   
                              + graupelncv(i,lat)
          graupel(i,lat) = fallsum_qg*delz(i,kts)/denr*dtcld*1000. + graupel(i,lat)
        ENDIF
        endif
!       if(fallsum.gt.0.)sr(i)=(snowncv(i,lat) + graupelncv(i,lat))/(rainncv(i)+1.e-12)
        if(fallsum.gt.0.)sr(i)=(tstepsnow(i) + tstepgraup(i))/(rainncv(i)+1.e-12)
      enddo
!
!---------------------------------------------------------------
! pimlt: instantaneous melting of cloud ice [HL A47] [RH83 A28]
!       (T>T0: I->C)
!---------------------------------------------------------------
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ IVDEP
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !GCC$ IVDEP
         !GCC$ VECTOR
#endif
        do i = its, ite
          supcol = t0c-t(i,k)
          xlf = xls-xl(i,k)
          if(supcol.lt.0.) xlf = xlf0
          if(supcol.lt.0.and.qci(i,k,2).gt.0.) then
            qci(i,k,1) = qci(i,k,1) + qci(i,k,2)
            t(i,k) = t(i,k) - xlf/cpm(i,k)*qci(i,k,2)
            qci(i,k,2) = 0.
          endif
!---------------------------------------------------------------
! pihmf: homogeneous freezing of cloud water below -40c [HL A45]
!        (T<-40C: C->I)
!---------------------------------------------------------------
          if(supcol.gt.40..and.qci(i,k,1).gt.0.) then
            qci(i,k,2) = qci(i,k,2) + qci(i,k,1)
            t(i,k) = t(i,k) + xlf/cpm(i,k)*qci(i,k,1)
            qci(i,k,1) = 0.
          endif
!---------------------------------------------------------------
! pihtf: heterogeneous freezing of cloud water [HL A44]
!        (T0>T>-40C: C->I)
!---------------------------------------------------------------
          if(supcol.gt.0..and.qci(i,k,1).gt.qmin) then
!           pfrzdtc = min(pfrz1*(exp(pfrz2*supcol)-1.)                         &
!              *den(i,k)/denr/xncr*qci(i,k,1)**2*dtcld,qci(i,k,1))
            supcolt=min(supcol,50.)
            pfrzdtc = min(pfrz1*(exp(pfrz2*supcolt)-1.)                        &
            *den(i,k)/denr/xncr*qci(i,k,1)*qci(i,k,1)*dtcld,qci(i,k,1))
            qci(i,k,2) = qci(i,k,2) + pfrzdtc
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtc
            qci(i,k,1) = qci(i,k,1)-pfrzdtc
          endif
!---------------------------------------------------------------
! pgfrz: freezing of rain water [HL A20] [LFO 45]
!        (T<T0, R->G)
!---------------------------------------------------------------
          if(supcol.gt.0..and.qrs(i,k,1).gt.0.) then
!           pfrzdtr = min(20.*pi**2*pfrz1*n0r*denr/den(i,k)                    &
!                 *(exp(pfrz2*supcol)-1.)*rslope3(i,k,1)**2                    &
!                 *rslope(i,k,1)*dtcld,qrs(i,k,1))
            temp = rslope3(i,k,1)
            temp = temp*temp*rslope(i,k,1)
            supcolt=min(supcol,50.)
            pfrzdtr = min(20.*(pi*pi)*pfrz1*n0r*denr/den(i,k)                  &
                  *(exp(pfrz2*supcolt)-1.)*temp*dtcld,                         &
                  qrs(i,k,1))
            qrs(i,k,3) = qrs(i,k,3) + pfrzdtr
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtr
            qrs(i,k,1) = qrs(i,k,1)-pfrzdtr
          endif
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     update the slope parameters for microphysics computation
!
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
         !DIR$ UNROLL(4)
#elif defined __GFORTRAN__
         !GCC$ VECTOR
         !GCC$ UNROLL 4
#endif
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
        enddo
      enddo
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                     work1,its,ite,kts,kte)
!------------------------------------------------------------------
!     work1:  the thermodynamic term in the denominator associated with
!             heat conduction and vapor diffusion
!             (ry88, y93, h85)
!     work2: parameter associated with the ventilation effects(y93)
!
      do k = kts, kte
        do i = its, ite
          work1(i,k,1) = diffac(xl(i,k),p(i,k),t(i,k),den(i,k),qs(i,k,1))
          work1(i,k,2) = diffac(xls,p(i,k),t(i,k),den(i,k),qs(i,k,2))
          work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
        enddo
      enddo
!
!===============================================================
!
! warm rain processes
!
! - follows the processes in RH83 and LFO except for autoconcersion
!
!===============================================================
!
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ IVDEP
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !GCC$ IVDEP
         !GCC$ VECTOR
#endif
        do i = its, ite
          supsat = max(q(i,k),qmin)-qs(i,k,1)
          satdt = supsat/dtcld
!---------------------------------------------------------------
! praut: auto conversion rate from cloud to rain [HDC 16]
!        (C->R)
!---------------------------------------------------------------
          if(qci(i,k,1).gt.qc0) then
            praut(i,k) = qck1*qci(i,k,1)**(7./3.)
            praut(i,k) = min(praut(i,k),qci(i,k,1)/dtcld)
          endif
!---------------------------------------------------------------
! pracw: accretion of cloud water by rain [HL A40] [LFO 51]
!        (C->R)
!---------------------------------------------------------------
          if(qrs(i,k,1).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            pracw(i,k) = min(pacrr*rslope3(i,k,1)*rslopeb(i,k,1)               &
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!---------------------------------------------------------------
! prevp: evaporation/condensation rate of rain [HDC 14]
!        (V->R or R->V)
!---------------------------------------------------------------
          if(qrs(i,k,1).gt.0.) then
            coeres = rslope2(i,k,1)*sqrt(rslope(i,k,1)*rslopeb(i,k,1))
            prevp(i,k) = (rh(i,k,1)-1.)*(precr1*rslope2(i,k,1)                 &
                         +precr2*work2(i,k)*coeres)/work1(i,k,1)
            if(prevp(i,k).lt.0.) then
              prevp(i,k) = max(prevp(i,k),-qrs(i,k,1)/dtcld)
              prevp(i,k) = max(prevp(i,k),satdt/2)
            else
              prevp(i,k) = min(prevp(i,k),satdt/2)
            endif
          endif
        enddo
      enddo
!
!===============================================================
!
! cold rain processes
!
! - follows the revised ice microphysics processes in HDC
! - the processes same as in RH83 and RH84  and LFO behave
!   following ice crystal hapits defined in HDC, inclduing
!   intercept parameter for snow (n0s), ice crystal number
!   concentration (ni), ice nuclei number concentration
!   (n0i), ice diameter (d)
!
!===============================================================
!
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ IVDEP
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !DIR$ IVDEP
         !DIR$ VECTOR
#endif
         do i = its, ite
#if defined __GFORTRAN__
            !GCC$ builtin (exp) attributes simd
#endif           
          supcol = t0c-t(i,k)
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          supsat = max(q(i,k),qmin)-qs(i,k,2)
          satdt = supsat/dtcld
          ifsat = 0
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
!         xni(i,k) = min(max(5.38e7*(den(i,k)                                  &
!                      *max(qci(i,k,2),qmin))**0.75,1.e3),1.e6)
          temp = (den(i,k)*max(qci(i,k,2),qmin))
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
          eacrs = exp(0.07*(-supcol))
!
          xmi = den(i,k)*qci(i,k,2)/xni(i,k)
          diameter  = min(dicon * sqrt(xmi),dimax)
          vt2i = 1.49e4*diameter**1.31
          vt2r=pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt2s=pvts*rslopeb(i,k,2)*denfac(i,k)
          vt2g=pvtg*rslopeb(i,k,3)*denfac(i,k)
          qsum(i,k) = max( (qrs(i,k,2)+qrs(i,k,3)), 1.E-15)
          if(qsum(i,k) .gt. 1.e-15) then
          vt2ave=(vt2s*qrs(i,k,2)+vt2g*qrs(i,k,3))/(qsum(i,k))
          else
          vt2ave=0.
          endif
          if(supcol.gt.0.and.qci(i,k,2).gt.qmin) then
            if(qrs(i,k,1).gt.qcrmin) then
!-------------------------------------------------------------
! praci: Accretion of cloud ice by rain [HL A15] [LFO 25]
!        (T<T0: I->R)
!-------------------------------------------------------------
              acrfac = 2.*rslope3(i,k,1)+2.*diameter*rslope2(i,k,1)            &
                      +diameter**2*rslope(i,k,1)
              praci(i,k) = pi*qci(i,k,2)*n0r*abs(vt2r-vt2i)*acrfac/4.
              praci(i,k) = min(praci(i,k),qci(i,k,2)/dtcld)
!-------------------------------------------------------------
! piacr: Accretion of rain by cloud ice [HL A19] [LFO 26]
!        (T<T0: R->S or R->G)
!-------------------------------------------------------------
              piacr(i,k) = pi**2*avtr*n0r*denr*xni(i,k)*denfac(i,k)            &
                          *g6pbr*rslope3(i,k,1)*rslope3(i,k,1)                 &
                          *rslopeb(i,k,1)/24./den(i,k)
              piacr(i,k) = min(piacr(i,k),qrs(i,k,1)/dtcld)
            endif
!-------------------------------------------------------------
! psaci: Accretion of cloud ice by snow [HDC 10]
!        (T<T0: I->S)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.qcrmin) then
              acrfac = 2.*rslope3(i,k,2)+2.*diameter*rslope2(i,k,2)            &
                      +diameter**2*rslope(i,k,2)
              psaci(i,k) = pi*qci(i,k,2)*eacrs*n0s*n0sfac(i,k)                 &
                          *abs(vt2ave-vt2i)*acrfac/4.
              psaci(i,k) = min(psaci(i,k),qci(i,k,2)/dtcld)
            endif
!-------------------------------------------------------------
! pgaci: Accretion of cloud ice by graupel [HL A17] [LFO 41]
!        (T<T0: I->G)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.qcrmin) then
              egi = exp(0.07*(-supcol))
              acrfac = 2.*rslope3(i,k,3)+2.*diameter*rslope2(i,k,3)            &
                      +diameter**2*rslope(i,k,3)
              pgaci(i,k) = pi*egi*qci(i,k,2)*n0g*abs(vt2ave-vt2i)*acrfac/4.
              pgaci(i,k) = min(pgaci(i,k),qci(i,k,2)/dtcld)
            endif
          endif
!-------------------------------------------------------------
! psacw: Accretion of cloud water by snow  [HL A7] [LFO 24]
!        (T<T0: C->S, and T>=T0: C->R)
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            psacw(i,k) = min(pacrc*n0sfac(i,k)*rslope3(i,k,2)*rslopeb(i,k,2)   &    
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pgacw: Accretion of cloud water by graupel [HL A6] [LFO 40]
!        (T<T0: C->G, and T>=T0: C->R)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            pgacw(i,k) = min(pacrg*rslope3(i,k,3)*rslopeb(i,k,3)               &
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! paacw: Accretion of cloud water by averaged snow/graupel 
!        (T<T0: C->G or S, and T>=T0: C->R) 
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qrs(i,k,3).gt.qcrmin) then
            paacw(i,k) = (qrs(i,k,2)*psacw(i,k)+qrs(i,k,3)*pgacw(i,k))         & 
                        /(qsum(i,k))
           endif      
!-------------------------------------------------------------
! pracs: Accretion of snow by rain [HL A11] [LFO 27]
!         (T<T0: S->G)
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qrs(i,k,1).gt.qcrmin) then
            if(supcol.gt.0) then
              acrfac = 5.*rslope3(i,k,2)*rslope3(i,k,2)*rslope(i,k,1)          &
                      +2.*rslope3(i,k,2)*rslope2(i,k,2)*rslope2(i,k,1)         &
                      +.5*rslope2(i,k,2)*rslope2(i,k,2)*rslope3(i,k,1)
              pracs(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2r-vt2ave)          &
                          *(dens/den(i,k))*acrfac
              pracs(i,k) = min(pracs(i,k),qrs(i,k,2)/dtcld)
            endif
!-------------------------------------------------------------
! psacr: Accretion of rain by snow [HL A10] [LFO 28]
!         (T<T0:R->S or R->G) (T>=T0: enhance melting of snow)
!-------------------------------------------------------------
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,2)            &
                    +2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,2)           &
                    +.5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,2)
            psacr(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2ave-vt2r)            &
                        *(denr/den(i,k))*acrfac
            psacr(i,k) = min(psacr(i,k),qrs(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pgacr: Accretion of rain by graupel [HL A12] [LFO 42]
!         (T<T0: R->G) (T>=T0: enhance melting of graupel)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qrs(i,k,1).gt.qcrmin) then
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,3)            &
                    +2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,3)           &
                    +.5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,3)
            pgacr(i,k) = pi**2*n0r*n0g*abs(vt2ave-vt2r)*(denr/den(i,k))        &
                        *acrfac
            pgacr(i,k) = min(pgacr(i,k),qrs(i,k,1)/dtcld)
          endif
!
!-------------------------------------------------------------
! pgacs: Accretion of snow by graupel [HL A13] [LFO 29]
!        (S->G): This process is eliminated in V3.0 with the 
!        new combined snow/graupel fall speeds
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qrs(i,k,2).gt.qcrmin) then
            pgacs(i,k) = 0.
          endif
          if(supcol.le.0) then
            xlf = xlf0
!-------------------------------------------------------------
! pseml: Enhanced melting of snow by accretion of water [HL A34]
!        (T>=T0: S->R)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0.)                                               &
              pseml(i,k) = min(max(cliq*supcol*(paacw(i,k)+psacr(i,k))         &
                          /xlf,-qrs(i,k,2)/dtcld),0.)
!-------------------------------------------------------------
! pgeml: Enhanced melting of graupel by accretion of water [HL A24] [RH84 A21-A22]
!        (T>=T0: G->R)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0.)                                               &
              pgeml(i,k) = min(max(cliq*supcol*(paacw(i,k)+pgacr(i,k))         &
                          /xlf,-qrs(i,k,3)/dtcld),0.)
          endif
          if(supcol.gt.0) then
!-------------------------------------------------------------
! pidep: Deposition/Sublimation rate of ice [HDC 9]
!       (T<T0: V->I or I->V)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.and.ifsat.ne.1) then
              pidep(i,k) = 4.*diameter*xni(i,k)*(rh(i,k,2)-1.)/work1(i,k,2)
              supice = satdt-prevp(i,k)
              if(pidep(i,k).lt.0.) then
                pidep(i,k) = max(max(pidep(i,k),satdt/2),supice)
                pidep(i,k) = max(pidep(i,k),-qci(i,k,2)/dtcld)
              else
                pidep(i,k) = min(min(pidep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)).ge.abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! psdep: deposition/sublimation rate of snow [HDC 14]
!        (T<T0: V->S or S->V)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psdep(i,k) = (rh(i,k,2)-1.)*n0sfac(i,k)*(precs1*rslope2(i,k,2)   &    
                           + precs2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)
              if(psdep(i,k).lt.0.) then
                psdep(i,k) = max(psdep(i,k),-qrs(i,k,2)/dtcld)
                psdep(i,k) = max(max(psdep(i,k),satdt/2),supice)
              else
                psdep(i,k) = min(min(psdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)).ge.abs(satdt))          &
                ifsat = 1
            endif
!-------------------------------------------------------------
! pgdep: deposition/sublimation rate of graupel [HL A21] [LFO 46]
!        (T<T0: V->G or G->V)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgdep(i,k) = (rh(i,k,2)-1.)*(precg1*rslope2(i,k,3)               &
                              +precg2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)
              if(pgdep(i,k).lt.0.) then
                pgdep(i,k) = max(pgdep(i,k),-qrs(i,k,3)/dtcld)
                pgdep(i,k) = max(max(pgdep(i,k),satdt/2),supice)
              else
                pgdep(i,k) = min(min(pgdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)+pgdep(i,k)).ge.          &
                abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! pigen: generation(nucleation) of ice from vapor [HL 50] [HDC 7-8]
!       (T<T0: V->I)
!-------------------------------------------------------------
            if(supsat.gt.0.and.ifsat.ne.1) then
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)-pgdep(i,k)
              xni0 = 1.e3*exp(0.1*supcol)
              roqi0 = 4.92e-11_sp*xni0**1.33#_sp
              pigen(i,k) = max(0.,(roqi0/den(i,k)-max(qci(i,k,2),0.))/dtcld)
              pigen(i,k) = min(min(pigen(i,k),satdt),supice)
            endif
!
!-------------------------------------------------------------
! psaut: conversion(aggregation) of ice to snow [HDC 12]
!        (T<T0: I->S)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.) then
              qimax = roqimax/den(i,k)
              psaut(i,k) = max(0.,(qci(i,k,2)-qimax)/dtcld)
            endif
!
!-------------------------------------------------------------
! pgaut: conversion(aggregation) of snow to graupel [HL A4] [LFO 37]
!        (T<T0: S->G)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0.) then
              alpha2 = 1.e-3*exp(0.09_sp*(-supcol))
              pgaut(i,k) = min(max(0.,alpha2*(qrs(i,k,2)-qs0)),qrs(i,k,2)/dtcld)
            endif
          endif
!
!-------------------------------------------------------------
! psevp: Evaporation of melting snow [HL A35] [RH83 A27]
!       (T>=T0: S->V)
!-------------------------------------------------------------
          if(supcol.lt.0.) then
            if(qrs(i,k,2).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psevp(i,k) = (rh(i,k,1)-1.)*n0sfac(i,k)*(precs1                  &
                           *rslope2(i,k,2)+precs2*work2(i,k)                   &
                           *coeres)/work1(i,k,1)
              psevp(i,k) = min(max(psevp(i,k),-qrs(i,k,2)/dtcld),0.)
            endif
!-------------------------------------------------------------
! pgevp: Evaporation of melting graupel [HL A25] [RH84 A19]
!       (T>=T0: G->V)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgevp(i,k) = (rh(i,k,1)-1.)*(precg1*rslope2(i,k,3)               &
                         +precg2*work2(i,k)*coeres)/work1(i,k,1)
              pgevp(i,k) = min(max(pgevp(i,k),-qrs(i,k,3)/dtcld),0.)
            endif
          endif
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     check mass conservation of generation terms and feedback to the
!     large scale
!
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#endif
        do i = its, ite
          delta2(i,k)=0._sp
          delta3(i,k)=0._sp
          if(qrs(i,k,1).lt.1.e-4.and.qrs(i,k,2).lt.1.e-4) delta2(i,k)=1.
          if(qrs(i,k,1).lt.1.e-4) delta3(i,k)=1.
        enddo
      enddo
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(3)
#endif
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL 5
#endif
#endif

        do i = its, ite
!
!     cloud water
!
          value = max(qmin,qci(i,k,1))
          source = (praut(i,k)+pracw(i,k)+paacw(i,k)+paacw(i,k))*dtcld
          if(t(i,k).le.t0c) then
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
            endif
          else
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
            endif
          endif
        enddo
      enddo
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(3)
#endif
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL 5
#endif
#endif
        do i = its, ite
!
!     cloud ice
!
          if(t(i,k).le.t0c) then
            value = max(qmin,qci(i,k,2))
            source = (psaut(i,k)-pigen(i,k)-pidep(i,k)+praci(i,k)+psaci(i,k)   &
                     +pgaci(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              psaut(i,k) = psaut(i,k)*factor
              pigen(i,k) = pigen(i,k)*factor
              pidep(i,k) = pidep(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
            endif
          endif
        enddo
      enddo
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(3)
#endif
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL 5
#endif
#endif
        do i = its, ite
!
!     rain
!
          value = max(qmin,qrs(i,k,1))
          if(t(i,k).le.t0c) then
            source = (-praut(i,k)-prevp(i,k)-pracw(i,k)+piacr(i,k)+psacr(i,k)  &
                      +pgacr(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
            endif
          else
            source = (-paacw(i,k)-praut(i,k)+pseml(i,k)+pgeml(i,k)-pracw(i,k)  &
                      -paacw(i,k)-prevp(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif
          endif
        enddo
      enddo
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(3)
#endif
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL 5
#endif
#endif
        do i = its, ite
!
!     snow
!
          if(t(i,k).le.t0c) then
            value = max(qmin,qrs(i,k,2))
            source = -(psdep(i,k)+psaut(i,k)-pgaut(i,k)+paacw(i,k)+piacr(i,k)  &
                     *delta3(i,k)+praci(i,k)*delta3(i,k)                       &
                     -pracs(i,k)*(1.-delta2(i,k))                              &
                     +psacr(i,k)*delta2(i,k)+psaci(i,k)-pgacs(i,k) )*dtcld
            if (source.gt.value) then
              factor = value/source
              psdep(i,k) = psdep(i,k)*factor
              psaut(i,k) = psaut(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif
          else
            value = max(qcrmin,qrs(i,k,2))
            source=(pgacs(i,k)-pseml(i,k)-psevp(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              psevp(i,k) = psevp(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
            endif
          endif
        enddo
      enddo
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(6)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(3)
#endif
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL 6
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL 5
#endif
#endif
        do i = its, ite
!
!     graupel
!
          if(t(i,k).le.t0c) then
            value = max(qmin,qrs(i,k,3))
            source = -(pgdep(i,k)+pgaut(i,k)                                   &
                     +piacr(i,k)*(1.-delta3(i,k))+praci(i,k)*(1.-delta3(i,k))  &
                     +psacr(i,k)*(1.-delta2(i,k))+pracs(i,k)*(1.-delta2(i,k))  &
                     +pgaci(i,k)+paacw(i,k)+pgacr(i,k)+pgacs(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgdep(i,k) = pgdep(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif
          else
            value = max(qcrmin,qrs(i,k,3))
            source=-(pgacs(i,k)+pgevp(i,k)+pgeml(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              pgevp(i,k) = pgevp(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif
          endif
        enddo
      enddo
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#endif
        do i = its, ite
!
!     update
!
          if(t(i,k).le.t0c) then
            work2(i,k)=-(prevp(i,k)+psdep(i,k)+pgdep(i,k)+pigen(i,k)+pidep(i,k))
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)                 &
                           +paacw(i,k)+paacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)                 &
                           +prevp(i,k)-piacr(i,k)-pgacr(i,k)                   &
                           -psacr(i,k))*dtcld,0.)
            qci(i,k,2) = max(qci(i,k,2)-(psaut(i,k)+praci(i,k)                 &
                           +psaci(i,k)+pgaci(i,k)-pigen(i,k)-pidep(i,k))       &
                           *dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psdep(i,k)+psaut(i,k)+paacw(i,k)      &
                           -pgaut(i,k)+piacr(i,k)*delta3(i,k)                  &
                           +praci(i,k)*delta3(i,k)+psaci(i,k)-pgacs(i,k)       &
                           -pracs(i,k)*(1.-delta2(i,k))+psacr(i,k)*delta2(i,k))&
                           *dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgdep(i,k)+pgaut(i,k)                 &
                           +piacr(i,k)*(1.-delta3(i,k))                        &
                           +praci(i,k)*(1.-delta3(i,k))                        &
                           +psacr(i,k)*(1.-delta2(i,k))                        &
                           +pracs(i,k)*(1.-delta2(i,k))+pgaci(i,k)+paacw(i,k)  &
                           +pgacr(i,k)+pgacs(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xls*(psdep(i,k)+pgdep(i,k)+pidep(i,k)+pigen(i,k))       &
                      -xl(i,k)*prevp(i,k)-xlf*(piacr(i,k)+paacw(i,k)           &
                      +paacw(i,k)+pgacr(i,k)+psacr(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          else
            work2(i,k)=-(prevp(i,k)+psevp(i,k)+pgevp(i,k))
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)                 &
                    +paacw(i,k)+paacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)                 &
                    +prevp(i,k)+paacw(i,k)+paacw(i,k)-pseml(i,k)               &
                    -pgeml(i,k))*dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psevp(i,k)-pgacs(i,k)                 &
                    +pseml(i,k))*dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgacs(i,k)+pgevp(i,k)                 &
                    +pgeml(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xl(i,k)*(prevp(i,k)+psevp(i,k)+pgevp(i,k))              &
                      -xlf*(pseml(i,k)+pgeml(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          endif
        enddo
      enddo
!
! Inline expansion for fpvs
!         qs(i,k,1) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!         qs(i,k,2) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !GCC$ VECTOR
         !GCC$ builtin (exp) attributes simd
         !GCC$ builtin (log) attributes simd
#endif
        do i = its, ite
          tr=ttp/t(i,k)
          qs(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          qs(i,k,1) = min(qs(i,k,1),0.99*p(i,k))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          tr=ttp/t(i,k)
          if(t(i,k).lt.ttp) then
            qs(i,k,2)=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
          else
            qs(i,k,2)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          endif
          qs(i,k,2) = min(qs(i,k,2),0.99*p(i,k))
          qs(i,k,2) = ep2 * qs(i,k,2) / (p(i,k) - qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
        enddo
      enddo
!
!----------------------------------------------------------------
!  pcond: condensational/evaporational rate of cloud water [HL A46] [RH83 A6]
!     if there exists additional water vapor condensated/if
!     evaporation of cloud water is not enough to remove subsaturation
!
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#endif
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED work1:64
            !DIR$ ASSUME_ALIGNED work2:64
#endif
          work1(i,k,1) = conden(t(i,k),q(i,k),qs(i,k,1),xl(i,k),cpm(i,k))
          work2(i,k) = qci(i,k,1)+work1(i,k,1)
          pcond(i,k) = min(max(work1(i,k,1)/dtcld,0.),max(q(i,k),0.)/dtcld)
          if(qci(i,k,1).gt.0..and.work1(i,k,1).lt.0.)                          &
            pcond(i,k) = max(work1(i,k,1),-qci(i,k,1))/dtcld
          q(i,k) = q(i,k)-pcond(i,k)*dtcld
          qci(i,k,1) = max(qci(i,k,1)+pcond(i,k)*dtcld,0.)
          t(i,k) = t(i,k)+pcond(i,k)*xl(i,k)/cpm(i,k)*dtcld
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     padding for small values
!
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
        do i = its, ite
          if(qci(i,k,1).le.qmin) qci(i,k,1) = 0.0
          if(qci(i,k,2).le.qmin) qci(i,k,2) = 0.0
        enddo
      enddo
      enddo                  ! big loops
  END SUBROUTINE wsm62d
  !--------------------------------------------------------------------------
#if defined __GFORTRAN__
      subroutine slope_wsm6(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
           vt,its,ite,kts,kte) !GCC$ ATTRIBUTES hot :: slope_wsm6 !GCC$ ATTRIBUTES always_inline :: slope_wsm6 !GCC$ ATTRIBUTES aligned(32) :: slope_wsm6 !GCC$ ATTRIBUTES target_clones ("avx,avx512")
#elif defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES inline :: slope_wsm6
      subroutine slope_wsm6(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
           vt,its,ite,kts,kte)
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: slope_wsm6
        !DIR$ ATTRIBUTES VECTOR :: slope_wsm6
#endif
  IMPLICIT NONE
  INTEGER       ::               its,ite, jts,jte, kts,kte
  REAL(kind=sp), DIMENSION( its:ite , kts:kte,3) ::                                     &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &                                                 
                                                                      rslope2, &                                                 
                                                                      rslope3, &                                                 
                                                                           vt
  REAL(kind=sp), DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL(kind=sp), PARAMETER  :: t0c = 273.15_sp
#if defined __INTEL_COMPILER
  !DIR$ ATTRIBUTES ALIGN : 64 :: n0sfac
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) ::                                      &
                                                                        n0sfac
#elif defined __GFORTRAN__
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) ::                                      &
                                                        n0sfac  !GCC$ ATTRIBUTES aligned(64) :: n0sfac
#endif
  REAL(kind=sp)       ::  lamdar, lamdas, lamdag, x, y, z, supcol
  integer(kind=int4)  :: i, j, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))      ! (pidn0g/(x*y))**.25
!



      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ VECTOR ALIGNED
         !DIR$ SIMD
#elif defined __GFORTRAN__
         !GCC$ VECTOR
#endif
        do i = its, ite
           supcol = t0c-t(i,k)
#if defined __INTEL_COMPILER
           !DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
#endif
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(i,k,1).le.qcrmin)then
            rslope(i,k,1) = rslopermax
            rslopeb(i,k,1) = rsloperbmax
            rslope2(i,k,1) = rsloper2max
            rslope3(i,k,1) = rsloper3max
          else
            rslope(i,k,1) = 1./lamdar(qrs(i,k,1),den(i,k))
            rslopeb(i,k,1) = rslope(i,k,1)**bvtr
            rslope2(i,k,1) = rslope(i,k,1)*rslope(i,k,1)
            rslope3(i,k,1) = rslope2(i,k,1)*rslope(i,k,1)
          endif
          if(qrs(i,k,2).le.qcrmin)then
            rslope(i,k,2) = rslopesmax
            rslopeb(i,k,2) = rslopesbmax
            rslope2(i,k,2) = rslopes2max
            rslope3(i,k,2) = rslopes3max
          else
            rslope(i,k,2) = 1./lamdas(qrs(i,k,2),den(i,k),n0sfac(i,k))
            rslopeb(i,k,2) = rslope(i,k,2)**bvts
            rslope2(i,k,2) = rslope(i,k,2)*rslope(i,k,2)
            rslope3(i,k,2) = rslope2(i,k,2)*rslope(i,k,2)
          endif
          if(qrs(i,k,3).le.qcrmin)then
            rslope(i,k,3) = rslopegmax
            rslopeb(i,k,3) = rslopegbmax
            rslope2(i,k,3) = rslopeg2max
            rslope3(i,k,3) = rslopeg3max
          else
            rslope(i,k,3) = 1./lamdag(qrs(i,k,3),den(i,k))
            rslopeb(i,k,3) = rslope(i,k,3)**bvtg
            rslope2(i,k,3) = rslope(i,k,3)*rslope(i,k,3)
            rslope3(i,k,3) = rslope2(i,k,3)*rslope(i,k,3)
          endif
          vt(i,k,1) = pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt(i,k,2) = pvts*rslopeb(i,k,2)*denfac(i,k)
          vt(i,k,3) = pvtg*rslopeb(i,k,3)*denfac(i,k)
          if(qrs(i,k,1).le.0.0) vt(i,k,1) = 0.0
          if(qrs(i,k,2).le.0.0) vt(i,k,2) = 0.0
          if(qrs(i,k,3).le.0.0) vt(i,k,3) = 0.0
        enddo
      enddo
  END subroutine slope_wsm6
!-----------------------------------------------------------------------------
#ifndef IINSIDE
#if defined __GFORTRAN__ 
      subroutine slope_rain(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   & 
       vt,kts,kte) !GCC$ ATTRIBUTES hot :: slope_rain !GCC$ ATTRIBUTES always_inline :: slope_rain !GCC$ ATTRIBUTES aligned(32) :: slope_rain !GCC$ ATTRIBUTES target_clones("avx,avx512") :: slope_rain
#elif defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES INLINE :: slope_rain
      subroutine slope_rain(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   & 
           vt,kts,kte)
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: slope_rain
        !DIR$ ATTIRBUTES VECTOR :: slope_rain
#endif
  IMPLICIT NONE
  INTEGER(kind=int4)       ::               kts,kte
  REAL(kind=sp), DIMENSION( kts:kte) ::                                                 &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &      
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL(kind=sp), PARAMETER  :: t0c = 273.15_sp
  REAL(kind=sp), DIMENSION( kts:kte ) ::                                                &
                                                                       n0sfac
  REAL(kind=sp)       ::  lamdar, x, y, z, supcol
  integer(kind=int4)  :: k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
      !
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
      !GCC$ VECTOR
#endif
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
#endif
          if(qrs(k).le.qcrmin)then
            rslope(k) = rslopermax
            rslopeb(k) = rsloperbmax
            rslope2(k) = rsloper2max
            rslope3(k) = rsloper3max
          else
            rslope(k) = 1./lamdar(qrs(k),den(k))
            rslopeb(k) = rslope(k)**bvtr
            rslope2(k) = rslope(k)*rslope(k)
            rslope3(k) = rslope2(k)*rslope(k)
          endif
          vt(k) = pvtr*rslopeb(k)*denfac(k)
          if(qrs(k).le.0.0) vt(k) = 0.0
      enddo
  END subroutine slope_rain
  !------------------------------------------------------------------------------
#if defined __GFORTRAN__
      subroutine slope_snow(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
      vt,kts,kte) !GCC$ ATTRIBUTES hot :: slope_snow !GCC$ ATTRIBUTES always_inline :: slope_snow !GCC$ ATTRIBUTES aligned(32) :: slope_snow !GCC$ ATTRIBUTES target_clones("avx,avx512") :: slope_snow
#elif defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES INLINE :: slope_snow 
      subroutine slope_snow(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3, &
           vt,kts,kte)
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: slope_snow
       !DIR$ ATTIRBUTES VECTOR :: slope_snow
#endif
  IMPLICIT NONE
  INTEGER(kind=int4)       ::               kts,kte
  REAL(kind=sp), DIMENSION( kts:kte) ::                                                 &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &  
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL(kind=sp), PARAMETER  :: t0c = 273.15_sp
  REAL(kind=sp), DIMENSION( kts:kte ) ::                                                &
                                                                       n0sfac
  REAL(kind=sp)       ::  lamdas, x, y, z, supcol
  integer(kind=int4) :: k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
!
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
      !GCC$ VECTOR
#endif
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
#endif
          supcol = t0c-t(k)
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(k).le.qcrmin)then
            rslope(k) = rslopesmax
            rslopeb(k) = rslopesbmax
            rslope2(k) = rslopes2max
            rslope3(k) = rslopes3max
          else
            rslope(k) = 1./lamdas(qrs(k),den(k),n0sfac(k))
            rslopeb(k) = rslope(k)**bvts
            rslope2(k) = rslope(k)*rslope(k)
            rslope3(k) = rslope2(k)*rslope(k)
          endif
          vt(k) = pvts*rslopeb(k)*denfac(k)
          if(qrs(k).le.0.0) vt(k) = 0.0
      enddo
  END subroutine slope_snow
  !----------------------------------------------------------------------------------
#if defined __GFORTRAN__
      
      subroutine slope_graup(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,  &
   vt,kts,kte) !GCC$ ATTRIBUTES hot :: slope_graup !GCC$ ATTRIBUTES always_inline :: slope_graup !GCC$ ATTRIBUTES aligned(32) :: slope_graup !GCC$ ATTRIBUTES target_clones("avx,avx512") :: slope_graup
#elif defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES INLINE :: slope_graup
      subroutine slope_graup(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3, &
           vt,kts,kte)
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: slope_graup
        !DIR$ ATTRIBUTES VECTOR :: slope_graup
#endif
  IMPLICIT NONE
  INTEGER(kind=int4)       :: kts,kte
  REAL(kind=sp), DIMENSION( kts:kte) ::                                                 &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &  
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL(kind=sp), PARAMETER  :: t0c = 273.15_sp
  REAL(kind=sp), DIMENSION( kts:kte ) ::                                                &
                                                                       n0sfac
  REAL(kind=sp)       ::  lamdag, x, y, z, supcol
  integer(kind=int4) :: j, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))      ! (pidn0g/(x*y))**.25
!
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
      !GCC$ VECTOR
#endif
      do k = kts, kte
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
#endif
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          if(qrs(k).le.qcrmin)then
            rslope(k) = rslopegmax
            rslopeb(k) = rslopegbmax
            rslope2(k) = rslopeg2max
            rslope3(k) = rslopeg3max
          else
            rslope(k) = 1./lamdag(qrs(k),den(k))
            rslopeb(k) = rslope(k)**bvtg
            rslope2(k) = rslope(k)*rslope(k)
            rslope3(k) = rslope2(k)*rslope(k)
          endif
          vt(k) = pvtg*rslopeb(k)*denfac(k)
          if(qrs(k).le.0.0) vt(k) = 0.0
      enddo
  END subroutine slope_graup
!---------------------------------------------------------------------------------
  !-------------------------------------------------------------------
#if defined __GFORTRAN__
  
  SUBROUTINE nislfv_rain_plm(its,ite,kts,kte,denl,denfacl,tkl,dzl,wwl,rql,precip,dt,id,iter)
    !GCC$ ATTRIBUTES hot :: nislfv_rain_plm !GCC$ ATTRIBUTES aligned(32) :: nislfv_rain_plm !GCC$ ATTRIBUTES target_clones("avx,avx512") :: nislfv_rain_plm
#elif defined __INTEL_COMPILER
  SUBROUTINE nislfv_rain_plm(its,ite,kts,kte,denl,denfacl,tkl,dzl,wwl,rql,precip,dt,id,iter)
    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: nislfv_rain_plm
    !DIR$ ATTRIBUTES VECTOR :: nislfv_rain_plm
#endif
!-------------------------------------------------------------------
!
! for non-iteration semi-Lagrangain forward advection for cloud
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise linear method
! this routine is under assumption of decfl < 1 for semi_Lagrangian
!
! dzl    depth of model layer in meter
! wwl    terminal velocity at model layer m/s
! rql    cloud density*mixing ration
! precip precipitation
! dt     time step
! id     kind of precip: 0 test case; 1 raindrop
! iter   how many time to guess mean terminal velocity: 0 pure forward.
!        0 : use departure wind for advection
!        1 : use mean wind for advection
!        > 1 : use mean wind after iter-1 iterations
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         implemented by song-you hong
!
   
      integer(kind=int4) ::   its,ite,kts,kte,id
      real(kind=sp)      ::  dt
      real(kind=sp),  dimension(its:ite,kts:kte) ::  dzl,wwl,rql
      real(kind=sp),  dimension(its:ite) :: precip
      real(kind=sp),  dimension(its:ite,kts:kte) ::  denl,denfacl,tkl
      !

      integer(kind=int4) ::  i,k,n,m,kk,kb,kt,iter
      real(kind=sp)      ::  tl,tl2,qql,dql,qqd
      real(kind=sp)      ::  th,th2,qqh,dqh
      real(kind=sp)      ::  zsum,qsum,dim,dip,c1,con1,fa1,fa2
      real(kind=sp)      ::  allold, allnew, zz, dzamin, cflmax, decfl
#if defined __GFORTRAN__
      real(kind=sp), dimension(kts:kte)   ::   dz  !GCC$ ATTRIBUTES aligned(64)    :: dz
      real(kind=sp), dimension(kts:kte)   ::   ww  !GCC$ ATTRIBUTES aligned(64)    :: ww
      real(kind=sp), dimension(kts:kte)   ::   qq  !GCC$ ATTRIBUTES aligned(64)    :: qq
      real(kind=sp), dimension(kts:kte)   ::   wd  !GCC$ ATTRIBUTES aligned(64)    :: wd
      real(kind=sp), dimension(kts:kte)   ::   wa  !GCC$ ATTRIBUTES aligned(64)    :: wa
      real(kind=sp), dimension(kts:kte)   ::   was !GCC$ ATTRIBUTES aligned(64)    :: was
      real(kind=sp), dimension(kts:kte)   ::   den !GCC$ ATTRIBUTES aligned(64)    :: den
      real(kind=sp), dimension(kts:kte)   ::   denfac !GCC$ ATTRIBUTES aligned(64) :: denfac
      real(kind=sp), dimension(kts:kte)   ::   tk  !GCC$ ATTRIBUTES aligned(64)    :: tk
      real(kind=sp), dimension(kts:kte+1) ::   wi  !GCC$ ATTRIBUTES aligned(64)    :: wi
      real(kind=sp), dimension(kts:kte+1) ::   zi  !GCC$ ATTRIBUTES aligned(64)    :: zi
      real(kind=sp), dimension(kts:kte+1) ::   za  !GCC$ ATTRIBUTES aligned(64)    :: za
      real(kind=sp), dimension(kts:kte)   ::   qn  !GCC$ ATTRIBUTES aligned(64)    :: qn
      real(kind=sp), dimension(kts:kte)   ::   qr  !GCC$ ATTRIBUTES aligned(64)    :: qr
      real(kind=sp), dimension(kts:kte)   ::   tmp1 !GCC$ ATTRIBUTES aligned(64)   :: tmp1
      real(kind=sp), dimension(kts:kte)   ::   tmp2 !GCC$ ATTRIBUTES aligned(64)   :: tmp2
      real(kind=sp), dimension(kts:kte)   ::   tmp3 !GCC$ ATTRIBUTES aligned(64)   :: tmp3
      real(kind=sp), dimension(kts:kte+1) ::   dza  !GCC$ ATTRIBUTES aligned(64)   :: dza
      real(kind=sp), dimension(kts:kte+1) ::   qa   !GCC$ ATTRIBUTES aligned(64)   :: qa
      real(kind=sp), dimension(kts:kte+1) ::   qmi  !GCC$ ATTRIBUTES aligned(64)   :: qmi
      real(kind=sp), dimension(kts:kte+1) ::   qpi  !GCC$ ATTRIBUTES aligned(64)   :: qpi
#elif defined __INTEL_COMPILER
      real(kind=sp), dimension(kts:kte)   ::   dz, ww, qq, wd, wa, was
      !DIR$ ATTRIBUTES ALIGN : 64 :: dz
      !DIR$ ATTRIBUTES ALIGN : 64 :: ww
      !DIR$ ATTRIBUTES ALIGN : 64 :: qq
      !DIR$ ATTRIBUTES ALIGN : 64 :: wd
      !DIR$ ATTRIBUTES ALIGN : 64 :: wa
      !DIR$ ATTRIBUTES ALIGN : 64 :: was
      real(kind=sp), dimension(kts:kte)   ::   den, denfac, tk
      !DIR$ ATTRIBUTES ALIGN : 64 :: den
      !DIR$ ATTRIBUTES ALIGN : 64 :: denfac
      !DIR$ ATTRIBUTES ALIGN : 64 :: tk
      real(kind=sp), dimension(kts:kte+1) ::   wi, zi, za
      !DIR$ ATTRIBUTES ALIGN : 64 :: wi
      !DIR$ ATTRIBUTES ALIGN : 64 :: zi
      !DIR$ ATTRIBUTES ALIGN : 64 :: za
      real(kind=sp), dimension(kts:kte)   ::   qn, qr,tmp,tmp1,tmp2,tmp3
      !DIR$ ATTRIBUTES ALIGN : 64 :: qn
      !DIR$ ATTRIBUTES ALIGN : 64 :: qr
      !DIR$ ATTRIBUTES ALIGN : 64 :: tmp1
      !DIR$ ATTRIBUTES ALIGN : 64 :: tmp2
      !DIR$ ATTRIBUTES ALIGN : 64 :: tmp3
      real(kind=sp), dimension(kts:kte+1) ::   dza, qa, qmi, qpi
      !DIR$ ATTRIBUTES ALIGN : 64 :: dza
      !DIR$ ATTRIBUTES ALIGN : 64 :: qa
      !DIR$ ATTRIBUTES ALIGN : 64 :: qmi
      !DIR$ ATTRIBUTES ALIGN : 64 :: qpi
#endif
!

      precip(:) = 0.0
!
      i_loop : do i=its,ite
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED denl:64,denfacl:64,tkl:64,dzl:64,wwl:64,rql:64,precip:64
#endif
! -----------------------------------
      dz(:) = dzl(i,:)
      qq(:) = rql(i,:)
      ww(:) = wwl(i,:)
      den(:) = denl(i,:)
      denfac(:) = denfacl(i,:)
      tk(:) = tkl(i,:)
! skip for no precipitation for all layers
      allold = 0.0
#if defined __INTEL_COMPILER
      !DIR$ SIMD REDUCTION(+:allold)
      do k=kts,kte
        allold = allold + qq(k)
      enddo
      if(allold.le.0.0) then
        cycle i_loop
      endif
!
! compute interface values
      zi(kts)=0.0
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
#elif defined __GFORTRAN__
      !GCC$ VECTOR
#endif
      do k=kts,kte
        zi(k+1) = zi(k)+dz(k)
      enddo
!
! save departure wind
      wd(:) = ww(:)
      n=1
 100  continue
! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! 2nd order interpolation to get wi
      wi(kts) = ww(kts)
      wi(kte+1) = ww(kte)
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif      
      do k=kts+1,kte
        wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
      enddo
! 3rd order interpolation to get wi
      fa1 = 9./16.
      fa2 = 1./16.
      wi(kts) = ww(kts)
      wi(kts+1) = 0.5*(ww(kts+1)+ww(kts))
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif
      do k=kts+2,kte-1
        wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
      enddo
      wi(kte) = 0.5*(ww(kte)+ww(kte-1))
      wi(kte+1) = ww(kte)
!
! terminate of top of raingroup
      do k=kts+1,kte
        if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
      enddo
!
! diffusivity of wi
      con1 = 0.05
      do k=kte,kts,-1
        decfl = (wi(k+1)-wi(k))*dt/dz(k)
        if( decfl .gt. con1 ) then
          wi(k) = wi(k+1) - con1*dz(k)/dt
        endif
      enddo
      ! compute arrival point
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif      
      do k=kts,kte+1
        za(k) = zi(k) - wi(k)*dt
      enddo
!
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif         
        dza(k) = za(k+1)-za(k)
      enddo
      dza(kte+1) = zi(kte+1) - za(kte+1)
!
      ! computer deformation at arrival point
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif      
      do k=kts,kte
        qa(k) = qq(k)*dz(k)/dza(k)
        qr(k) = qa(k)/den(k)
      enddo
      qa(kte+1) = 0.0
!     call maxmin(kte-kts+1,1,qa,' arrival points ')
!
! compute arrival terminal velocity, and estimate mean terminal velocity
! then back to use mean terminal velocity
      if( n.le.iter ) then
        call slope_rain(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,kts,kte)
        if( n.ge.2 ) wa(kts:kte)=0.5*(wa(kts:kte)+was(kts:kte))
        do k=kts,kte
!#ifdef DEBUG
!        print*,' slope_wsm3 ',qr(k)*1000.,den(k),denfac(k),tk(k),tmp(k),tmp1(k),tmp2(k),ww(k),wa(k)
!#endif
! mean wind is average of departure and new arrival winds
          ww(k) = 0.5* ( wd(k)+wa(k) )
        enddo
        was(:) = wa(:)
        n=n+1
        go to 100
      endif
!
      ! estimate values at arrival cell interface with monotone
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif      
      do k=kts+1,kte
        dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
        dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
        if( dip*dim.le.0.0 ) then
          qmi(k)=qa(k)
          qpi(k)=qa(k)
        else
          qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
          qmi(k)=2.0*qa(k)-qpi(k)
          if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          endif
        endif
      enddo
      qpi(kts)=qa(kts)
      qmi(kts)=qa(kts)
      qmi(kte+1)=qa(kte+1)
      qpi(kte+1)=qa(kte+1)
!
! interpolation to regular point
      qn = 0.0
      kb=kts
      kt=kts
      intp : do k=kts,kte
             kb=max(kb-1,kts)
             kt=max(kt-1,kts)
! find kb and kt
             if( zi(k).ge.za(kte+1) ) then
               exit intp
             else
               find_kb : do kk=kb,kte
                         if( zi(k).le.za(kk+1) ) then
                           kb = kk
                           exit find_kb
                         else
                           cycle find_kb
                         endif
               enddo find_kb
               find_kt : do kk=kt,kte
                         if( zi(k+1).le.za(kk) ) then
                           kt = kk
                           exit find_kt
                         else
                           cycle find_kt
                         endif
               enddo find_kt
               kt = kt - 1
! compute q with piecewise constant method
               if( kt.eq.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 th=(zi(k+1)-za(kb))/dza(kb)
                 tl2=tl*tl
                 th2=th*th
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qqh=qqd*th2+qmi(kb)*th
                 qql=qqd*tl2+qmi(kb)*tl
                 qn(k) = (qqh-qql)/(th-tl)
               else if( kt.gt.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 tl2=tl*tl
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qql=qqd*tl2+qmi(kb)*tl
                 dql = qa(kb)-qql
                 zsum  = (1.-tl)*dza(kb)
                 qsum  = dql*dza(kb)
                 if( kt-kb.gt.1 ) then
                 do m=kb+1,kt-1
                   zsum = zsum + dza(m)
                   qsum = qsum + qa(m) * dza(m)
                 enddo
                 endif
                 th=(zi(k+1)-za(kt))/dza(kt)
                 th2=th*th
                 qqd=0.5*(qpi(kt)-qmi(kt))
                 dqh=qqd*th2+qmi(kt)*th
                 zsum  = zsum + th*dza(kt)
                 qsum  = qsum + dqh*dza(kt)
                 qn(k) = qsum/zsum
               endif
               cycle intp
             endif
!
       enddo intp
!
! rain out
      sum_precip: do k=kts,kte
                    if( za(k).lt.0.0 .and. za(k+1).lt.0.0 ) then
                      precip(i) = precip(i) + qa(k)*dza(k)
                      cycle sum_precip
                    else if ( za(k).lt.0.0 .and. za(k+1).ge.0.0 ) then
                      precip(i) = precip(i) + qa(k)*(0.0-za(k))
                      exit sum_precip
                    endif
                    exit sum_precip
      enddo sum_precip
!
! replace the new values

      rql(i,:) = qn(:)
!
! ----------------------------------
      enddo i_loop
!
    END SUBROUTINE nislfv_rain_plm
    
    !-------------------------------------------------------------------
#if defined __GFORTRAN__
    SUBROUTINE nislfv_rain_plm6(its,ite,kts,kte,denl,denfacl,tkl,dzl,wwl,rql,rql2, &
         precip1, precip2,dt,id,iter) !GCC$ ATTRIBUTES hot :: nislfv_rain_plm6 !GCC$ ATTRIBUTES aligned(32) :: nislfv_rain_plm6 !GCC$ ATTRIBUTES target_clones("avx,avx512") :: nislfv_rain_plm6
#elif defined __INTEL_COMPILER
   SUBROUTINE nislfv_rain_plm6(its,ite,kts,kte,denl,denfacl,tkl,dzl,wwl,rql,rql2, &
        precip1, precip2,dt,id,iter)
     !DIR$ ATTRIBUTE CODE_ALIGN : 32 :: nislfv_rain_plm6
     !DIR$ ATTRIBUTES VECTOR :: nislfv_rain_plm6
#endif
!-------------------------------------------------------------------
!
! for non-iteration semi-Lagrangain forward advection for cloud
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise linear method
! this routine is under assumption of decfl < 1 for semi_Lagrangian
!
! dzl    depth of model layer in meter
! wwl    terminal velocity at model layer m/s
! rql    cloud density*mixing ration
! precip precipitation
! dt     time step
! id     kind of precip: 0 test case; 1 raindrop
! iter   how many time to guess mean terminal velocity: 0 pure forward.
!        0 : use departure wind for advection
!        1 : use mean wind for advection
!        > 1 : use mean wind after iter-1 iterations
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         implemented by song-you hong
!
      
      integer(kind=int4) ::   its,ite,kts,kte,id
      real(kind=sp)      ::  dt
      real(kind=sp), dimension(its:ite,kts:kte) ::  dzl,wwl,rql,rql2
      real(kind=sp), dimension(its:ite)         ::  precip,precip1,precip2
      real(kind=sp), dimension(its:ite,kts:kte) ::  denl,denfacl,tkl
!
      integer(kind=int4) ::  i,k,n,m,kk,kb,kt,iter,ist
      real(kind=sp)      ::  tl,tl2,qql,dql,qqd
      real(kind=sp)      ::  th,th2,qqh,dqh
      real(kind=sp)      ::  zsum,qsum,dim,dip,c1,con1,fa1,fa2
      real(kind=sp)      ::  allold, allnew, zz, dzamin, cflmax, decfl
#if defined __GFORTRAN__
      real(kind=sp), dimension(kts:kte)   ::  dz  !GCC$ ATTRIBUTES aligned(64) :: dz
      real(kind=sp), dimension(kts:kte)   ::  ww  !GCC$ ATTRIBUTES aligned(64) :: ww
      real(kind=sp), dimension(kts:kte)   ::  qq  !GCC$ ATTRIBUTES aligned(64) :: qq
      real(kind=sp), dimension(kts:kte)   ::  qq2 !GCC$ ATTRIBUTES aligned(64) :: qq2
      real(kind=sp), dimension(kts:kte)   ::  wd  !GCC$ ATTRIBUTES aligned(64) :: wd
      real(kind=sp), dimension(kts:kte)   ::  wa  !GCC$ ATTRIBUTES aligned(64) :: wa
      real(kind=sp), dimension(kts:kte)   ::  wa2 !GCC$ ATTRIBUTES aligned(64) :: wa2
      real(kind=sp), dimension(kts:kte)   ::  was !GCC$ ATTRIBUTES aligned(64) :: was
      real(kind=sp), dimension(kts:kte)   ::  den !GCC$ ATTRIBUTES aligned(64) :: den
      real(kind=sp), dimension(kts:kte)   ::  denfac  !GCC$ ATTRIBUTES aligned(64) :: denfac
      real(kind=sp), dimension(kts:kte)   ::  tk  !GCC$ ATTRIBUTES aligned(64) :: tk
      real(kind=sp), dimension(kts:kte+1) ::  wi  !GCC$ ATTRIBUTES aligned(64) :: wi
      real(kind=sp), dimension(kts:kte+1) ::  zi  !GCC$ ATTRIBUTES aligned(64) :: zi
      real(kind=sp), dimension(kts:kte+1) ::  za  !GCC$ ATTRIBUTES aligned(64) :: za
      real(kind=sp), dimension(kts:kte)   ::  qn  !GCC$ ATTRIBUTES aligned(64) :: qn
      real(kind=sp), dimension(kts:kte)   ::  qr  !GCC$ ATTRIBUTES aligned(64) :: qr
      real(kind=sp), dimension(kts:kte)   ::  qr2 !GCC$ ATTRIBUTES aligned(64) :: qr2
      real(kind=sp), dimension(kts:kte)   ::  tmp !GCC$ ATTRIBUTES aligned(64) :: tmp
      real(kind=sp), dimension(kts:kte)   ::  tmp1 !GCC$ ATTRIBUTES aligned(64) :: tmp1
      real(kind=sp), dimension(kts:kte)   ::  tmp2 !GCC$ ATTRIBUTES aligned(64) :: tmp2
      real(kind=sp), dimension(kts:kte)   ::  tmp3 !GCC$ ATTRIBUTES aligned(64) :: tmp3
      real(kind=sp), dimension(kts:kte+1) ::  dza  !GCC$ ATTRIBUTES aligned(64) :: dza
      real(kind=sp), dimension(kts:kte)   ::  qa   !GCC$ ATTRIBUTES aligned(64) :: qa
      real(kind=sp), dimension(kts:kte)   ::  qa2  !GCC$ ATTRIBUTES aligned(64) :: qa2
      real(kind=sp), dimension(kts:kte)   ::  qmi  !GCC$ ATTRIBUTES aligned(64) :: qmi
      real(kind=sp), dimension(kts:kte)   ::  qpi  !GCC$ ATTRIBUTES aligned(64) :: qpi
#elif defined __INTEL_COMPILER
      real(kind=sp), dimension(kts:kte)   ::  dz, ww, qq, qq2, wd, wa, wa2, was
      !DIR$ ATTRIBUTES ALIGN : 64 :: dz
      !DIR$ ATTRIBUTES ALIGN : 64 :: ww
      !DIR$ ATTRIBUTES ALIGN : 64 :: qq
      !DIR$ ATTRIBUTES ALIGN : 64 :: qq2
      !DIR$ ATTRIBUTES ALIGN : 64 :: wd
      !DIR$ ATTRIBUTES ALIGN : 64 :: wa
      !DIR$ ATTRIBUTES ALIGN : 64 :: wa2
      !DIR$ ATTRIBUTES ALIGN : 64 :: was
      real(kind=sp), dimension(kts:kte)   ::  den, denfac, tk
      !DIR$ ATTRIBUTES ALIGN : 64 :: den
      !DIR$ ATTRIBUTES ALIGN : 64 :: denfac
      !DIR$ ATTRIBUTES ALIGN : 64 :: tk
      real(kind=sp), dimension(kts:kte+1) ::  wi, zi, za
      !DIR$ ATTRIBUTES ALIGN : 64 :: wi
      !DIR$ ATTRIBUTES ALIGN : 64 :: zi
      !DIR$ ATTRIBUTES ALIGN : 64 :: za
      real(kind=sp), dimension(kts:kte)   ::  qn, qr,qr2,tmp,tmp1,tmp2,tmp3
      !DIR$ ATTRIBUTES ALIGN : 64 :: qn
      !DIR$ ATTRIBUTES ALIGN : 64 :: qr
      !DIR$ ATTRIBUTES ALIGN : 64 :: qr2
      !DIR$ ATTRIBUTES ALIGN : 64 :: tmp
      !DIR$ ATTRIBUTES ALIGN : 64 :: tmp1
      !DIR$ ATTRIBUTES ALIGN : 64 :: tmp2
      !DIR$ ATTRIBUTES ALIGN : 64 :: tmp3
      real(kind=sp), dimension(kts:kte+1) ::  dza, qa, qa2,qmi, qpi
      !DIR$ ATTRIBUTES ALIGN : 64 :: dza
      !DIR$ ATTRIBUTES ALIGN : 64 :: qa
      !DIR$ ATTRIBUTES ALIGN : 64 :: qa2
      !DIR$ ATTRIBUTES ALIGN : 64 :: qmi
      !DIR$ ATTRIBUTES ALIGN : 64 :: qpi
#endif
!

      precip(:) = 0.0_sp
      precip1(:) = 0.0_sp
      precip2(:) = 0.0_sp
!
      i_loop : do i=its,ite
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED denl:64,denfacl:64,tkl:64,dzl:64,wwl:64,rql:64,rql2:64,precip1:64,precip2:64
#endif
! -----------------------------------
      dz(:) = dzl(i,:)
      qq(:) = rql(i,:)
      qq2(:) = rql2(i,:)
      ww(:) = wwl(i,:)
      den(:) = denl(i,:)
      denfac(:) = denfacl(i,:)
      tk(:) = tkl(i,:)
! skip for no precipitation for all layers
      allold = 0.0
#if defined __INTEL_COMPILER
      !DIR$ SIMD REDUCTION(+:allold)
#endif
      do k=kts,kte
        allold = allold + qq(k)
      enddo
      if(allold.le.0.0_sp) then
        cycle i_loop
      endif
!
! compute interface values
      zi(kts)=0.0_sp
      
      do k=kts,kte
        zi(k+1) = zi(k)+dz(k)
      enddo
!
! save departure wind
      wd(:) = ww(:)
      n=1
 100  continue
! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! 2nd order interpolation to get wi
      wi(kts) = ww(kts)
      wi(kte+1) = ww(kte)
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif 
      do k=kts+1,kte
        wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
      enddo
! 3rd order interpolation to get wi
      fa1 = 9./16.
      fa2 = 1./16.
      wi(kts) = ww(kts)
      wi(kts+1) = 0.5*(ww(kts+1)+ww(kts))
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif       
      do k=kts+2,kte-1
        wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
      enddo
      wi(kte) = 0.5*(ww(kte)+ww(kte-1))
      wi(kte+1) = ww(kte)
!
! terminate of top of raingroup
      do k=kts+1,kte
        if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
      enddo
!
! diffusivity of wi
      con1 = 0.05_sp
      do k=kte,kts,-1
        decfl = (wi(k+1)-wi(k))*dt/dz(k)
        if( decfl .gt. con1 ) then
          wi(k) = wi(k+1) - con1*dz(k)/dt
        endif
      enddo
      ! compute arrival point
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif       
      do k=kts,kte+1
        za(k) = zi(k) - wi(k)*dt
      enddo
      !
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif       
      do k=kts,kte
        dza(k) = za(k+1)-za(k)
      enddo
      dza(kte+1) = zi(kte+1) - za(kte+1)
!
      ! computer deformation at arrival point
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif       
      do k=kts,kte
        qa(k) = qq(k)*dz(k)/dza(k)
        qa2(k) = qq2(k)*dz(k)/dza(k)
        qr(k) = qa(k)/den(k)
        qr2(k) = qa2(k)/den(k)
      enddo
      qa(kte+1) = 0.0_sp
      qa2(kte+1) = 0.0_sp
!     call maxmin(kte-kts+1,1,qa,' arrival points ')
!
! compute arrival terminal velocity, and estimate mean terminal velocity
! then back to use mean terminal velocity
      if( n.le.iter ) then
        call slope_snow(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,kts,kte)
        call slope_graup(qr2,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa2,kts,kte)
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ SIMD IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif         
        do k = kts, kte
          tmp(k) = max((qr(k)+qr2(k)), 1.E-15_sp)
          IF ( tmp(k) .gt. 1.e-15_sp ) THEN
            wa(k) = (wa(k)*qr(k) + wa2(k)*qr2(k))/tmp(k)
          ELSE
            wa(k) = 0._sp
          ENDIF
        enddo
        if( n.ge.2 ) wa(kts:kte)=0.5*(wa(kts:kte)+was(kts:kte))
        do k=kts,kte
!#ifdef DEBUG
!        print*,' slope_wsm3 ',qr(k)*1000.,den(k),denfac(k),tk(k),tmp(k),tmp1(k),tmp2(k), &
!           ww(k),wa(k)
!#endif
! mean wind is average of departure and new arrival winds
          ww(k) = 0.5* ( wd(k)+wa(k) )
        enddo
        was(:) = wa(:)
        n=n+1
        go to 100
      endif
      ist_loop : do ist = 1, 2
      if (ist.eq.2) then
       qa(:) = qa2(:)
      endif
!
      precip(i) = 0.
!
! estimate values at arrival cell interface with monotone
      do k=kts+1,kte
        dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
        dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
        if( dip*dim.le.0.0 ) then
          qmi(k)=qa(k)
          qpi(k)=qa(k)
        else
          qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
          qmi(k)=2.0*qa(k)-qpi(k)
          if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          endif
        endif
      enddo
      qpi(kts)=qa(kts)
      qmi(kts)=qa(kts)
      qmi(kte+1)=qa(kte+1)
      qpi(kte+1)=qa(kte+1)
!
! interpolation to regular point
      qn = 0.0
      kb=kts
      kt=kts
      intp : do k=kts,kte
             kb=max(kb-1,kts)
             kt=max(kt-1,kts)
! find kb and kt
             if( zi(k).ge.za(kte+1) ) then
               exit intp
             else
               find_kb : do kk=kb,kte
                         if( zi(k).le.za(kk+1) ) then
                           kb = kk
                           exit find_kb
                         else
                           cycle find_kb
                         endif
               enddo find_kb
               find_kt : do kk=kt,kte
                         if( zi(k+1).le.za(kk) ) then
                           kt = kk
                           exit find_kt
                         else
                           cycle find_kt
                         endif
               enddo find_kt
               kt = kt - 1
! compute q with piecewise constant method
               if( kt.eq.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 th=(zi(k+1)-za(kb))/dza(kb)
                 tl2=tl*tl
                 th2=th*th
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qqh=qqd*th2+qmi(kb)*th
                 qql=qqd*tl2+qmi(kb)*tl
                 qn(k) = (qqh-qql)/(th-tl)
               else if( kt.gt.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 tl2=tl*tl
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qql=qqd*tl2+qmi(kb)*tl
                 dql = qa(kb)-qql
                 zsum  = (1.-tl)*dza(kb)
                 qsum  = dql*dza(kb)
                 if( kt-kb.gt.1 ) then
                 do m=kb+1,kt-1
                   zsum = zsum + dza(m)
                   qsum = qsum + qa(m) * dza(m)
                 enddo
                 endif
                 th=(zi(k+1)-za(kt))/dza(kt)
                 th2=th*th
                 qqd=0.5*(qpi(kt)-qmi(kt))
                 dqh=qqd*th2+qmi(kt)*th
                 zsum  = zsum + th*dza(kt)
                 qsum  = qsum + dqh*dza(kt)
                 qn(k) = qsum/zsum
               endif
               cycle intp
             endif
!
       enddo intp
!
! rain out
      sum_precip: do k=kts,kte
                    if( za(k).lt.0.0 .and. za(k+1).lt.0.0 ) then
                      precip(i) = precip(i) + qa(k)*dza(k)
                      cycle sum_precip
                    else if ( za(k).lt.0.0 .and. za(k+1).ge.0.0 ) then
                      precip(i) = precip(i) + qa(k)*(0.0-za(k))
                      exit sum_precip
                    endif
                    exit sum_precip
      enddo sum_precip
!
! replace the new values
      if(ist.eq.1) then
        rql(i,:) = qn(:)
        precip1(i) = precip(i)
      else
        rql2(i,:) = qn(:)
        precip2(i) = precip(i)
      endif
      enddo ist_loop
!
! ----------------------------------
      enddo i_loop
!
  END SUBROUTINE nislfv_rain_plm6
!---------------------------------------------------------------------------------
#else
  !-------------------------------------------------------------------
#if defined __GFORTRAN__
                 subroutine slope_rain_ii(qrs,den,denfac,t,rslope,rslopeb,&
 rslope2,rslope3, vt,its,ite,kts,kte,lmask) !GCC$ ATTRIBUTES hot :: slope_rain_ii !GCC$ ATTRIBUTES always_inline :: slope_rain_ii !GCC$ ATTRIBUTES aligned(32) :: slope_rain_ii !GCC$ ATTRIBUTES target_clones("avx,avx512") :: slope_rain_ii
#elif defined __INTEL_COMPILER
                  !DIR$ ATTRIBUTES INLINE :: slope_rain_ii
                  subroutine slope_rain_ii(qrs,den,denfac,t,rslope,rslopeb,&
rslope2,rslope3, vt,its,ite,kts,kte,lmask)
                    !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: slope_rain_ii
                    !DIR$ ATTRIBUTES VECTOR :: slope_rain_ii
#endif
                    
 
  INTEGER(kind=int4)       :: its,ite,kts,kte
  REAL(kind=sp), DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &      
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL(kind=sp), PARAMETER  :: t0c = 273.15_sp
  LOGICAL(kind=int4) :: lmask(its:ite)
  REAL(kind=sp)       ::  lamdar, x, y, z, supcol
  integer :: i, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
!

      do k = kts, kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif

         do i = its, ite
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
#endif            
         if (lmask(i)) then
          if(qrs(i,k).le.qcrmin)then
            rslope(i,k) = rslopermax
            rslopeb(i,k) = rsloperbmax
            rslope2(i,k) = rsloper2max
            rslope3(i,k) = rsloper3max
          else
            rslope(i,k) = 1./lamdar(qrs(i,k),den(i,k))
            rslopeb(i,k) = rslope(i,k)**bvtr
            rslope2(i,k) = rslope(i,k)*rslope(i,k)
            rslope3(i,k) = rslope2(i,k)*rslope(i,k)
          endif
          vt(i,k) = pvtr*rslopeb(i,k)*denfac(i,k)
          if(qrs(i,k).le.0.0) vt(i,k) = 0.0
         endif
        enddo
      enddo
  END subroutine slope_rain_ii
  !------------------------------------------------------------------------------
#if defined __GFORTRAN__
      subroutine slope_snow_ii(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,&
vt,its,ite,kts,kte,lmask) !GCC$ ATTRIBUTES hot :: slope_snow_ii !GCC$ ATTRIBUTES always_inline :: slope_snow_ii !GCC$ ATTRIBUTES aligned(32) :: slope_snow_ii !GCC$ ATTRIBUTES target_clones("avx,avx512") :: slope_snow_ii
#elif defined __INTEL_COMPILER
        subroutine slope_snow_ii(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,&
             vt,its,ite,kts,kte,lmask)
#endif
  
  INTEGER(kind=int4)       :: its,ite,kts,kte
  REAL(kind=sp), DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &  
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL(kind=sp), PARAMETER  :: t0c = 273.15_sp
#if defined __GFORTRAN__
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) ::       n0sfac      !GCC$ ATTRIBUTES aligned(64) :: n0sfac
#elif defined __INTEL_COMPILER
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ) ::       n0sfac 
  !DIR$ ATTRIBUTES ALIGN : 64 :: n0sfac
#endif
  LOGICAL(kind=int4) :: lmask(its:ite)
  REAL(kind=sp)       ::  lamdas, x, y, z, supcol
  integer(kind=int4) :: i, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
!

      do k = kts, kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif         
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64,n0sfac:64
#endif
         if (lmask(i)) then
          supcol = t0c-t(i,k)
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(i,k).le.qcrmin)then
            rslope(i,k) = rslopesmax
            rslopeb(i,k) = rslopesbmax
            rslope2(i,k) = rslopes2max
            rslope3(i,k) = rslopes3max
          else
            rslope(i,k) = 1./lamdas(qrs(i,k),den(i,k),n0sfac(i,k))
            rslopeb(i,k) = rslope(i,k)**bvts
            rslope2(i,k) = rslope(i,k)*rslope(i,k)
            rslope3(i,k) = rslope2(i,k)*rslope(i,k)
          endif
          vt(i,k) = pvts*rslopeb(i,k)*denfac(i,k)
          if(qrs(i,k).le.0.0) vt(i,k) = 0.0
         endif
        enddo
      enddo
  END subroutine slope_snow_ii
  !----------------------------------------------------------------------------------
#if defined __GFORTRAN__  
      subroutine slope_graup_ii(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,&
vt,its,ite,kts,kte,lmask) !GCC$ ATTRIBUTES hot :: slope_graup_ii !GCC$ ATTRIBUTES always_inline :: slope_graup_ii !GCC$ ATTRIBUTES aligned(32) :: slope_graup_ii !GCC$ ATTRIBUTES target_clones("avx,avx512") :: slope_graup_ii

  INTEGER(kind=int4)      :: its,ite,kts,kte
  REAL(kind=sp), DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &  
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL(kind=sp), PARAMETER  :: t0c = 273.15_sp
  LOGICAL(kind=int4) :: lmask(its:ite)
  REAL(kind=sp)       ::  lamdag, x, y, z, supcol
  integer(kind=int4) :: i, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))      ! (pidn0g/(x*y))**.25
!

      do k = kts, kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif           
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
#endif
         if (lmask(i)) then
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          if(qrs(i,k).le.qcrmin)then
            rslope(i,k) = rslopegmax
            rslopeb(i,k) = rslopegbmax
            rslope2(i,k) = rslopeg2max
            rslope3(i,k) = rslopeg3max
          else
            rslope(i,k) = 1./lamdag(qrs(i,k),den(i,k))
            rslopeb(i,k) = rslope(i,k)**bvtg
            rslope2(i,k) = rslope(i,k)*rslope(i,k)
            rslope3(i,k) = rslope2(i,k)*rslope(i,k)
          endif
          vt(i,k) = pvtg*rslopeb(i,k)*denfac(i,k)
          if(qrs(i,k).le.0.0) vt(i,k) = 0.0
         endif
        enddo
      enddo
  END subroutine slope_graup_ii
!---------------------------------------------------------------------------------
  !-------------------------------------------------------------------
#if defined __GFORTRAN__
  SUBROUTINE nislfv_rain_plm_ii(its,ite,kts,kte,denl,denfacl, &
       tkl,dzl,wwl,rql,precip,dt,id,iter) !GCC$ ATTRIBUTES hot :: nislfv_rain_plm_ii !GCC$ ATTRIBUTES aligned(32) :: nislfv_rain_plm_ii !GCC$ ATTRIBUTES
#elif defined __INTEL_COMPILER
   SUBROUTINE nislfv_rain_plm_ii(its,ite,kts,kte,denl,denfacl, &
        tkl,dzl,wwl,rql,precip,dt,id,iter)
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: nislfv_rain_plm_ii
#endif
!-------------------------------------------------------------------
!
! for non-iteration semi-Lagrangain forward advection for cloud
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise linear method
! this routine is under assumption of decfl < 1 for semi_Lagrangian
!
! dzl    depth of model layer in meter
! wwl    terminal velocity at model layer m/s
! rql    cloud density*mixing ration
! precip precipitation
! dt     time step
! id     kind of precip: 0 test case; 1 raindrop
! iter   how many time to guess mean terminal velocity: 0 pure forward.
!        0 : use departure wind for advection
!        1 : use mean wind for advection
!        > 1 : use mean wind after iter-1 iterations
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         implemented by song-you hong
!
      
      integer(kind=int4) ::   its,ite,kts,kte,id
      real(kind=sp)      ::   dt
      real(kind=sp), dimension(its:ite,kts:kte) ::  dzl,wwl,rql
      real(kind=sp), dimension(its:ite)         ::  precip
      real(kind=sp), dimension(its:ite,kts:kte) ::  denl,denfacl,tkl
!
      integer(kind=int4) ::   i,k,n,m,kk,iter
#ifdef MASK_HISTOGRAM
#if defined __GFORTRAN__
      integer(kind=int4), dimension(kts:kte)     :: intp_count   !GCC$ ATTRIBUTES aligned(64) :: intp_count
      integer(kind=int4), dimension(0:ite-its+1) :: intp_hist    !GCC$ ATTRIBUTES aligned(64) :: intp_hist
#elif defined __INTEL_COMPILER
      integer(kind=int4), dimension(kts:kte)     :: intp_count
      !DIR$ ATTRIBUTES ALIGN : 64 :: intp_count
      integer(kind=int4), dimension(0:ite-its+1) :: intp_hist
#endif
#endif
      real(kind=sp) ::  dim,dip,con1,fa1,fa2,decfl
#if defined __GFORTRAN__
      real(kind=sp),  dimension(its:ite) ::  allold            !GCC$ ATTRIBUTES aligned(64) :: allold
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  dz      !GCC$ ATTRIBUTES aligned(64) :: dz
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  ww      !GCC$ ATTRIBUTES aligned(64) :: ww
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  qq      !GCC$ ATTRIBUTES aligned(64) :: qq
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  wd      !GCC$ ATTRIBUTES aligned(64) :: wd
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  wa      !GCC$ ATTRIBUTES aligned(64) :: wa
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  was     !GCC$ ATTRIBUTES aligned(64) :: was
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  den     !GCC$ ATTRIBUTES aligned(64) :: den
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  denfac  !GCC$ ATTRIBUTES aligned(64) :: denfac
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  tk      !GCC$ ATTRIBUTES aligned(64) :: tk
      real(kind=sp),  dimension(its:ite,kts:kte+1) ::  wi      !GCC$ ATTRIBUTES aligned(64) :: wi
      real(kind=sp),  dimension(its:ite,kts:kte+1) ::  zi      !GCC$ ATTRIBUTES aligned(64) :: zi
      real(kind=sp),  dimension(its:ite,kts:kte+1) ::  za      !GCC$ ATTRIBUTES aligned(64) :: za
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  qn      !GCC$ ATTRIBUTES aligned(64) :: qn
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  qr      !GCC$ ATTRIBUTES aligned(64) :: qr
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  tmp     !GCC$ ATTRIBUTES aligned(64) :: tmp
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  tmp1    !GCC$ ATTRIBUTES aligned(64) :: tmp1
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  tmp2    !GCC$ ATTRIBUTES aligned(64) :: tmp2
      real(kind=sp),  dimension(its:ite,kts:kte)   ::  tmp3    !GCC$ ATTRIBUTES aligned(64) :: tmp3
      real(kind=sp),  dimension(its:ite,kts:kte+1) ::  dza     !GCC$ ATTRIBUTES aligned(64) :: dza
      real(kind=sp),  dimension(its:ite,kts:kte+1) ::  qa      !GCC$ ATTRIBUTES aligned(64) :: qa
      real(kind=sp),  dimension(its:ite,kts:kte+1) ::  qmi     !GCC$ ATTRIBUTES aligned(64) :: qmi
      real(kind=sp),  dimension(its:ite,kts:kte+1) ::  qpi     !GCC$ ATTRIBUTES aligned(64) :: qpi
      logical(kind=int4), dimension(its:ite)       ::  lmask   !GCC$ ATTRIBUTES aligned(64) :: lmask
#elif defined __INTEL_COMPILER
      real(kind=sp), dimension(its:ite) ::  allold
      !DIR$ ATTRIBUTES ALIGN : 64 :: allold
      real(kind=sp), dimension(its:ite,kts:kte)   ::   dz,ww, qq, wd, wa, was
      !DIR$ ATTRIBUTES ALIGN : 64 :: dz,ww, qq, wd, wa, was
      real(kind=sp), dimension(its:ite,kts:kte)   ::   den, denfac, tk
      !DIR$ ATTRIBUTES ALIGN : 64 :: den, denfac, tk
      real(kind=sp), dimension(its:ite,kts:kte+1) ::   wi, zi, za
      !DIR$ ATTRIBUTES ALIGN : 64 :: wi, zi, za
      real(kind=sp), dimension(its:ite,kts:kte)   ::   qn, qr,tmp,tmp1,tmp2,tmp3
      !DIR$ ATTRIBUTES ALIGN : 64 :: qn, qr,tmp,tmp1,tmp2,tmp3
      real(kind=sp), dimension(its:ite,kts:kte)   ::   dza, qa, qmi, qpi
      !DIR$ ATTRIBUTES ALIGN : 64 :: dza, qa, qmi, qpi
      logical(kind=int4), dimension(its:ite)      ::   lmask
      !DIR$ ATTRIBUTES ALIGN : 64 :: lmask
#endif
      !
#if defined __GFORTRAN__
      INTEGER(kind=int4) ::  minkb, minkt
      LOGICAL(kind=int4), DIMENSION(its:ite) :: intp_mask  !GCC$ ATTRIBUTES aligned(64) :: intp_mask
      LOGICAL(kind=int4), DIMENSION(its:ite) :: tmask      !GCC$ ATTRIBUTES aligned(64) :: tmask
      INTEGER(kind=int4), DIMENSION(its:ite) :: kb         !GCC$ ATTRIBUTES aligned(64) :: kb
      INTEGER(kind=int4), DIMENSION(its:ite) :: kt         !GCC$ ATTRIBUTES aligned(64) :: kt
      REAL(kind=sp),      DIMENSION(its:ite) :: tl         !GCC$ ATTRIBUTES aligned(64) :: tl
      REAL(kind=sp),      DIMENSION(its:ite) :: tl2        !GCC$ ATTRIBUTES aligned(64) :: tl2
      REAL(kind=sp),      DIMENSION(its:ite) :: th         !GCC$ ATTRIBUTES aligned(64) :: th
      REAL(kind=sp),      DIMENSION(its:ite) :: th2        !GCC$ ATTRIBUTES aligned(64) :: th2
      REAL(kind=sp),      DIMENSION(its:ite) :: qqd        !GCC$ ATTRIBUTES aligned(64) :: qqd
      REAL(kind=sp),      DIMENSION(its:ite) :: qqh        !GCC$ ATTRIBUTES aligned(64) :: qqh
      REAL(kind=sp),      DIMENSION(its:ite) :: qql        !GCC$ ATTRIBUTES aligned(64) :: qql
      REAL(kind=sp),      DIMENSION(its:ite) :: zsum       !GCC$ ATTRIBUTES aligned(64) :: zsum
      REAL(kind=sp),      DIMENSION(its:ite) :: qsum       !GCC$ ATTRIBUTES aligned(64) :: qsum
      REAL(kind=sp),      DIMENSION(its:ite) :: dql        !GCC$ ATTRIBUTES aligned(64) :: dql
      REAL(kind=sp),      DIMENSION(its:ite) :: dqh        !GCC$ ATTRIBUTES aligned(64) :: dqh
      REAL(kind=sp),      DIMENSION(its:ite) :: za_gath_t  !GCC$ ATTRIBUTES aligned(64) :: za_gath_t
      REAL(kind=sp),      DIMENSION(its:ite) :: za_gath_b  !GCC$ ATTRIBUTES aligned(64) :: za_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: qa_gath_b  !GCC$ ATTRIBUTES aligned(64) :: qa_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: dza_gath_t !GCC$ ATTRIBUTES aligned(64) :: dza_gath_t
      REAL(kind=sp),      DIMENSION(its:ite) :: dza_gath_b !GCC$ ATTRIBUTES aligned(64) :: dza_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: qpi_gath_t !GCC$ ATTRIBUTES aligned(64) :: qpi_gath_t
      REAL(kind=sp),      DIMENSION(its:ite) :: qpi_gath_b !GCC$ ATTRIBUTES aligned(64) :: qpi_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: qmi_gath_t !GCC$ ATTRIBUTES aligned(64) :: qmi_gath_t
      REAL(kind=sp),      DIMENSION(its:ite) :: qmi_gath_b !GCC$ ATTRIBUTES aligned(64) :: qmi_gath_b
#elif defined __INTEL_COMPILER
      INTEGER(kind=int4) ::  minkb, minkt
      LOGICAL(kind=int4), DIMENSION(its:ite) :: intp_mask, tmask
      !DIR$ ATTRIBUTES ALIGN : 64 :: intp_mask, tmask
      INTEGER(kind=int4), DIMENSION(its:ite) :: kb, kt
      !DIR$ ATTRIBUTES ALIGN : 64 :: kb, kt
      REAL(kind=sp),      DIMENSION(its:ite) :: tl,tl2,th,th2,qqd,qqh,qql,zsum,qsum,dql,dqh
      !DIR$ ATTRIBUTES ALIGN : 64 :: tl,tl2,th,th2,qqd,qqh,qql,zsum,qsum,dql,dqh
      REAL(kind=sp),      DIMENSION(its:ite) :: za_gath_t,za_gath_b
      !DIR$ ATTRIBUTES ALIGN : 64 :: za_gath_t,za_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: qa_gath_b
      !DIR$ ATTRIBUTES ALIGN : 64 :: qa_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: dza_gath_t,dza_gath_b
      !DIR$ ATTRIBUTES ALIGN : 64 :: dza_gath_t,dza_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: qpi_gath_t,qpi_gath_b
      !DIR$ ATTRIBUTES ALIGN : 64 :: qpi_gath_t,qpi_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: qmi_gath_t,qmi_gath_b
      !DIR$ ATTRIBUTES ALIGN : 64 ::  qmi_gath_t,qmi_gath_b
#endif
!


#if defined __INTEL_COMPILER
!DIR$ ASSUME_ALIGNED denl:64,denfacl:64,tkl:64,dzl:64,wwl:64,rql:64,precip:64,lmask:64
!DIR$ ASSUME_ALIGNED intp_mask:64,tmask:64,kb:64,kt:64,tl:64,tl2:64,qqd:64,qqh:64
!DIR$ ASSUME_ALIGNED qql:64,zsum:64,dql:64,dqh:64,za_gath_t:64,za_gath_b:64,qa_gath_b:64
!DIR$ ASSUME_ALIGNED dza_gath_t:64,dza_gath_b:64,qpi_gath_t:64,qpi_gath_b:64
!DIR$ ASSUME_ALIGNED qmi_gath_t:64,qmi_gath_b:64,wi:64,ww:64,za:64,zi:64,dza:64
#endif
      precip(:) = 0.0_sp
!
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif          
         do i=its,ite

! -----------------------------------
      dz(i,k) = dzl(i,k)
      qq(i,k) = rql(i,k)
      ww(i,k) = wwl(i,k)
      den(i,k) = denl(i,k)
      denfac(i,k) = denfacl(i,k)
      tk(i,k) = tkl(i,k)
      enddo
      enddo
! skip for no precipitation for all layers
      do i=its,ite
      allold(i) = 0.0
      enddo
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif          
      do i=its,ite
        allold(i) = allold(i) + qq(i,k)
      enddo
      enddo
      if (maxval(allold).le.0.0) return
      lmask = allold .gt. 0.0
!
! compute interface values
      do i=its,ite
      if (lmask(i)) then
      zi(i,kts)=0.0
      endif
      enddo
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif          
      do i=its,ite
      if (lmask(i)) then
        zi(i,k+1) = zi(i,k)+dz(i,k)
      endif
      enddo
      enddo
!
! save departure wind
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif          
      do i=its,ite
      if (lmask(i)) then
      wd(i,k) = ww(i,k)
      endif
      enddo
      enddo
      n=1
      do while (n.le.(iter+1))
      do i=its,ite
      if (lmask(i)) then
! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! 2nd order interpolation to get wi
      wi(i,kts) = ww(i,kts)
      wi(i,kte+1) = ww(i,kte)
      endif
      enddo
      do k=kts+1,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif          
      do i=its,ite
      if (lmask(i)) then
        wi(i,k) = (ww(i,k)*dz(i,k-1)+ww(i,k-1)*dz(i,k))/(dz(i,k-1)+dz(i,k))
      endif
      enddo
      enddo
! 3rd order interpolation to get wi
      fa1 = 9./16.
      fa2 = 1./16.
      do i=its,ite
      if (lmask(i)) then
      wi(i,kts) = ww(i,kts)
      wi(i,kts+1) = 0.5*(ww(i,kts+1)+ww(i,kts))
      endif
      enddo
      do k=kts+2,kte-1
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif           
      do i=its,ite
      if (lmask(i)) then
        wi(i,k) = fa1*(ww(i,k)+ww(i,k-1))-fa2*(ww(i,k+1)+ww(i,k-2))
      endif
      enddo
      enddo
      do i=its,ite
      if (lmask(i)) then
      wi(i,kte) = 0.5*(ww(i,kte)+ww(i,kte-1))
      wi(i,kte+1) = ww(i,kte)
      endif
      enddo
!
! terminate of top of raingroup
      do k=kts+1,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif  
      do i=its,ite
      if (lmask(i)) then
        if( ww(i,k).eq.0.0 ) wi(i,k)=ww(i,k-1)
      endif
      enddo
      enddo
!
! diffusivity of wi
      con1 = 0.05_sp
      do k=kte,kts,-1
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif           
      do i=its,ite
      if (lmask(i)) then
        decfl = (wi(i,k+1)-wi(i,k))*dt/dz(i,k)
        if( decfl .gt. con1 ) then
          wi(i,k) = wi(i,k+1) - con1*dz(i,k)/dt
        endif
      endif
      enddo
      enddo
! compute arrival point
      do k=kts,kte+1
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif           
      do i=its,ite
      if (lmask(i)) then
        za(i,k) = zi(i,k) - wi(i,k)*dt
      endif
      enddo
      enddo
!
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif           
      do i=its,ite
      if (lmask(i)) then
        dza(i,k) = za(i,k+1)-za(i,k)
      endif
      enddo
      enddo
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif     
      do i=its,ite
      if (lmask(i)) then
      dza(i,kte+1) = zi(i,kte+1) - za(i,kte+1)
      endif
      enddo
!
! computer deformation at arrival point
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif           
      do i=its,ite
      if (lmask(i)) then
        qa(i,k) = qq(i,k)*dz(i,k)/dza(i,k)
        qr(i,k) = qa(i,k)/den(i,k)
      endif
      enddo
      enddo
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif  
      do i=its,ite
      if (lmask(i)) then
      qa(i,kte+1) = 0.0
!     call maxmin(kte-kts+1,1,qa(i,:),' arrival points ')
      endif
      enddo
!
! compute arrival terminal velocity, and estimate mean terminal velocity
! then back to use mean terminal velocity
      if( n.le.iter ) then
        call slope_rain_ii(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,its,ite,kts,kte,lmask)
        if( n.ge.2 ) then
        do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif  
        do i=its,ite
        if (lmask(i)) then
          wa(i,k)=0.5*(wa(i,k)+was(i,k))
        endif
        enddo
        enddo
        endif
        do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif  
        do i=its,ite
        if (lmask(i)) then
!#ifdef DEBUG
!        print*,' slope_wsm3 ',qr(i,k)*1000.,den(i,k),denfac(i,k),tk(i,k),tmp(i,k),tmp1(i,k),tmp2(i,k),ww(i,k),wa(i,k)
!#endif
! mean wind is average of departure and new arrival winds
          ww(i,k) = 0.5* ( wd(i,k)+wa(i,k) )
          was(i,k) = wa(i,k)
        endif
        enddo
        enddo
      endif
      n=n+1
      enddo
!
! estimate values at arrival cell interface with monotone
      do k=kts+1,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif           
      do i=its,ite
      if (lmask(i)) then
        dip=(qa(i,k+1)-qa(i,k))/(dza(i,k+1)+dza(i,k))
        dim=(qa(i,k)-qa(i,k-1))/(dza(i,k-1)+dza(i,k))
        if( dip*dim.le.0.0 ) then
          qmi(i,k)=qa(i,k)
          qpi(i,k)=qa(i,k)
        else
          qpi(i,k)=qa(i,k)+0.5*(dip+dim)*dza(i,k)
          qmi(i,k)=2.0*qa(i,k)-qpi(i,k)
          if( qpi(i,k).lt.0.0 .or. qmi(i,k).lt.0.0 ) then
            qpi(i,k) = qa(i,k)
            qmi(i,k) = qa(i,k)
          endif
        endif
      endif
      enddo
      enddo
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif     
      do i=its,ite
      if (lmask(i)) then
      qpi(i,kts)=qa(i,kts)
      qmi(i,kts)=qa(i,kts)
      qmi(i,kte+1)=qa(i,kte+1)
      qpi(i,kte+1)=qa(i,kte+1)
      endif
      enddo
!
! interpolation to regular point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGIN WSM5 CODE FROM John Michalakes !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      qn = 0.0
      kb=kts  ! kb is a vector
      kt=kts  ! kt is a vector
#ifdef MASK_HISTOGRAM
      intp_hist = 0
#endif
      INTP : do k=kts,kte
             kb=max(kb-1,kts)
             kt=max(kt-1,kts)
! find kb and kt
             intp_mask = ( zi(:,k).lt.za(:,kte+1) .AND. lmask )
             tmask = intp_mask
             minkb = 999
             minkt = 999
             DO i=its,ite
               IF ( tmask(i) .AND. kb(i) .lt. minkb ) minkb = kb(i)
               IF ( tmask(i) .AND. kt(i) .lt. minkt ) minkt = kt(i)
             ENDDO
             find_kb : do kk=minkb,kte
               DO i=its,ite
               IF ( tmask(i) .AND. zi(i,k).le.za(i,kk+1) ) THEN
                 kb(i) = kk
                 tmask(i) = .FALSE.
               ENDIF
               ENDDO
             enddo find_kb

             tmask = intp_mask
             find_kt : do kk=minkt,kte
               DO i=its,ite
               IF ( tmask(i) .AND. zi(i,k+1).le.za(i,kk) ) THEN
                 kt(i) = kk
                 tmask(i) = .FALSE.
               ENDIF
               ENDDO
             enddo find_kt
             kt = max(kt - 1,kts)

!#define RANGE_CHECKING
#ifndef RANGE_CHECKING
# define DX1 (i+(kb(i)-1)*(ite-its+1)),1
# define DX2 (i+(kt(i)-1)*(ite-its+1)),1
#else
# define DX1 i,kb(i)
# define DX2 i,kt(i)
#endif
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif  
             DO i = its,ite
               qa_gath_b(i) = qa(DX1)
               za_gath_b(i) = za(DX1)
               dza_gath_b(i) = dza(DX1)
               qpi_gath_b(i) = qpi(DX1)
               qmi_gath_b(i) = qmi(DX1)
             ENDDO
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif  
             DO i = its,ite
               za_gath_t(i) = za(DX2)
               dza_gath_t(i) = dza(DX2)
               qpi_gath_t(i) = qpi(DX2)
               qmi_gath_t(i) = qmi(DX2)
             ENDDO
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif  
             DO i = its,ite
             IF ( kt(i) .eq. kb(i) .AND. intp_mask(i) ) THEN
               tl(i)=(zi(i,k)-za_gath_b(i))/dza_gath_b(i)
               th(i)=(zi(i,k+1)-za_gath_b(i))/dza_gath_b(i)
               tl2(i) = tl(i)*tl(i)
               th2(i) = th(i)*th(i)
               qqd(i)=0.5*(qpi_gath_b(i)-qmi_gath_b(i))
               qqh(i)=qqd(i)*th2(i)+qmi_gath_b(i)*th(i)
               qql(i)=qqd(i)*tl2(i)+qmi_gath_b(i)*tl(i)
               qn(i,k) = (qqh(i)-qql(i))/(th(i)-tl(i))
             ELSE IF ( kt(i) .gt. kb(i) .AND. intp_mask(i) ) THEN
               tl(i)=(zi(i,k)-za_gath_b(i))/dza_gath_b(i)
               tl2(i)=tl(i)*tl(i)
               qqd(i)=0.5*(qpi_gath_b(i)-qmi_gath_b(i))
               qql(i)=qqd(i)*tl2(i)+qmi_gath_b(i)*tl(i)
               dql(i) = qa_gath_b(i)-qql(i)
               zsum(i)  = (1.-tl(i))*dza_gath_b(i)
               qsum(i)  = dql(i)*dza_gath_b(i)
             ENDIF
             ENDDO
#ifdef MASK_HISTOGRAM
             intp_count(k) = 0
             DO i = its,ite
               IF ( kt(i) .ge. kb(i) .AND. intp_mask(i) ) THEN
                 intp_count(k) = intp_count(k) + 1
               ENDIF
             ENDDO
#endif
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif               
             DO i = its,ite
               if( kt(i)-kb(i).gt.1 .AND. intp_mask(i) ) then
                 do m=kb(i)+1,kt(i)-1
                     zsum(i) = zsum(i) + dza(i,m)
                     qsum(i) = qsum(i) + qa(i,m) * dza(i,m)
                 enddo
               endif
            ENDDO
 #if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif             
             DO i = its,ite
             IF ( kt(i) .gt. kb(i) .AND. intp_mask(i) ) THEN
               th(i)=(zi(i,k+1)-za_gath_t(i))/dza_gath_t(i)
               th2(i) = th(i)*th(i)
               qqd(i)=0.5*(qpi_gath_t(i)-qmi_gath_t(i))
               dqh(i)=qqd(i)*th2(i)+qmi_gath_t(i)*th(i)
               zsum(i)  = zsum(i) + th(i)*dza_gath_t(i)
               qsum(i)  = qsum(i) + dqh(i)*dza_gath_t(i)
               qn(i,k) = qsum(i)/zsum(i)
             ENDIF
             ENDDO
       ENDDO intp
#ifdef MASK_HISTOGRAM
       do k=kts,kte
!print *,'DEBUG:  intp_count(',k,') = ',intp_count(k)
         IF ((intp_count(k) < 0) .OR. (intp_count(k) > (ite-its+1))) THEN
           print *,'ERROR:  intp_count(',k,') = ',intp_count(k)
           stop
         ENDIF
         intp_hist(intp_count(k)) = intp_hist(intp_count(k)) + 1
       enddo
       if (ite-its+1 == 8) then
         write (6,110) intp_hist
110   format ('intp_hist =  ',9i3)
       else
         do i=its,ite
           print *,'intp_hist(',i,') = ',intp_hist(i)
         enddo
       endif
#endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END WSM5 CODE FROM John Michalakes !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! rain out
      intp_mask = lmask
      sum_precip: do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif           
             DO i = its,ite
             IF (za(i,k).lt.0.0.and.za(i,k+1).lt.0.0.AND.intp_mask(i)) THEN
               precip(i) = precip(i) + qa(i,k)*dza(i,k)
             ELSE IF (za(i,k).lt.0.0.and.za(i,k+1).ge.0.0.AND.intp_mask(i)) THEN
               precip(i) = precip(i) + qa(i,k)*(0.0-za(i,k))
               intp_mask(i) = .FALSE.
             ENDIF
             ENDDO
      enddo sum_precip
!
! replace the new values
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif  
      do i=its,ite
      if (lmask(i)) then
      rql(i,k) = qn(i,k)
      endif
      enddo
      enddo
!
! ----------------------------------
!
  END SUBROUTINE nislfv_rain_plm_ii
  !-------------------------------------------------------------------
#if defined __GFORTRAN__
  SUBROUTINE nislfv_rain_plm6_ii(its,ite,kts,kte,denl,denfacl,tkl,&
       dzl,wwl,rql,rql2, precip1, precip2,dt,id,iter) !GCC$ ATTRIBUTES hot :: nislfv_rain_plm6_ii !GCC$ ATTRIBUTES aligned(32) :: nislfv_rain_plm6_ii
#elif defined __INTEL_COMPILER
   SUBROUTINE nislfv_rain_plm6_ii(its,ite,kts,kte,denl,denfacl,tkl,&
        dzl,wwl,rql,rql2, precip1, precip2,dt,id,iter)
     !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: nislfv_rain_plm6_ii
#endif
!-------------------------------------------------------------------
!
! for non-iteration semi-Lagrangain forward advection for cloud
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise linear method
! this routine is under assumption of decfl < 1 for semi_Lagrangian
!
! dzl    depth of model layer in meter
! wwl    terminal velocity at model layer m/s
! rql    cloud density*mixing ration
! precip precipitation
! dt     time step
! id     kind of precip: 0 test case; 1 raindrop
! iter   how many time to guess mean terminal velocity: 0 pure forward.
!        0 : use departure wind for advection
!        1 : use mean wind for advection
!        > 1 : use mean wind after iter-1 iterations
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         implemented by song-you hong
!
     
      integer(kind=int4) ::   its,ite,kts,kte,id
      real(kind=sp)      ::   dt
      real(kind=sp), dimension(its:ite,kts:kte) ::  dzl,wwl,rql,rql2
      real(kind=sp), dimension(its:ite)         ::  precip,precip1,precip2
      real(kind=sp), dimension(its:ite,kts:kte) ::  denl,denfacl,tkl
!
      integer(kind=int4) ::   i,k,n,m,kk,iter,ist
      real(kind=sp)      ::   dim,dip,con1,fa1,fa2,defcl
#if defined __GFORTRAN__
      real(kind=sp), dimension(its:ite) ::   allold
      real(kind=sp), dimension(its:ite,kts:kte)   ::  dz      !GCC$ ATTRIBUTES aligned(64) :: dz
      real(kind=sp), dimension(its:ite,kts:kte)   ::  ww      !GCC$ ATTRIBUTES aligned(64) :: ww
      real(kind=sp), dimension(its:ite,kts:kte)   ::  qq      !GCC$ ATTRIBUTES aligned(64) :: qq
      real(kind=sp), dimension(its:ite,kts:kte)   ::  qq2     !GCC$ ATTRIBUTES aligned(64) :: qq2
      real(kind=sp), dimension(its:ite,kts:kte)   ::  wd      !GCC$ ATTRIBUTES aligned(64) :: wd
      real(kind=sp), dimension(its:ite,kts:kte)   ::  wa      !GCC$ ATTRIBUTES aligned(64) :: wa
      real(kind=sp), dimension(its:ite,kts:kte)   ::  wa2     !GCC$ ATTRIBUTES aligned(64) :: wa2
      real(kind=sp), dimension(its:ite,kts:kte)   ::  was     !GCC$ ATTRIBUTES aligned(64) :: was
      real(kind=sp), dimension(its:ite,kts:kte)   ::  den     !GCC$ ATTRIBUTES aligned(64) :: den
      real(kind=sp), dimension(its:ite,kts:kte)   ::  denfac  !GCC$ ATTRIBUTES aligned(64) :: denfac
      real(kind=sp), dimension(its:ite,kts:kte)   ::  tk      !GCC$ ATTRIBUTES aligned(64) :: tk
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  wi      !GCC$ ATTRIBUTES aligned(64) :: wi
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  zi      !GCC$ ATTRIBUTES aligned(64) :: zi
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  za      !GCC$ ATTRIBUTES aligned(64) :: za
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  qn      !GCC$ ATTRIBUTES aligned(64) :: qn
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  qr      !GCC$ ATTRIBUTES aligned(64) :: qr
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  qr2     !GCC$ ATTRIBUTES aligned(64) :: qr2
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  tmp     !GCC$ ATTRIBUTES aligned(64) :: tmp
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  tmp1    !GCC$ ATTRIBUTES aligned(64) :: tmp1
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  tmp2    !GCC$ ATTRIBUTES aligned(64) :: tmp2
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  tmp3    !GCC$ ATTRIBUTES aligned(64) :: tmp3
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  dza     !GCC$ ATTRIBUTES aligned(64) :: dza
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  qa      !GCC$ ATTRIBUTES aligned(64) :: qa
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  qa2     !GCC$ ATTRIBUTES aligned(64) :: qa2
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  qmi     !GCC$ ATTRIBUTES aligned(64) :: qmi
      real(kind=sp), dimension(its:ite,kts:kte+1) ::  qpi     !GCC$ ATTRIBUTES aligned(64) :: qpi
      logical(kind=int4), dimension(its:ite)      ::  lmask   !GCC$ ATTRIBUTES aligned(64) :: lmask
#elif defined __INTEL_COMPILER
      real(kind=sp), dimension(its:ite)           ::   allold
      !DIR$ ATTRIBUTES ALIGN : 64 :: allold
      real(kind=sp), dimension(its:ite,kts:kte)   ::   dz, ww, qq, qq2, wd, wa, wa2, was
      !DIR$ ATTRIBUTES ALIGN : 64 :: dz, ww, qq, qq2, wd, wa, wa2, was
      real(kind=sp), dimension(its:ite,kts:kte)   ::   den, denfac, tk
      !DIR$ ATTRIBUTES ALIGN : 64 ::  den, denfac, tk
      real(kind=sp), dimension(its:ite,kts:kte+1) ::   wi zi, za
      !DIR$ ATTRIBUTES ALIGN : 64 ::   wi zi, za
      real(kind=sp), dimension(its:ite,kts:kte)   ::   qn,qr,qr2,tmp,tmp1,tmp2,tmp3
      !DIR$ ATTRIBUTES ALIGN : 64 ::    qn,qr,qr2,tmp,tmp1,tmp2,tmp3
      real(kind=sp), dimension(its:ite,kts:kte+1) ::   dza, qa,qa2,qmi,qpi
      !DIR$ ATTRIBUTES ALIGN : 64 ::     dza, qa,qa2,qmi,qpi
      logical(kind=int4), dimension(its:ite)      ::   lmask
      !DIR$ ATTRIBUTES ALIGN : 64 ::     lmask
#endif
!
      INTEGER(kind=int4) ::  minkb, minkt
#if defined __GFORTRAN__
      LOGICAL(kind=int4), DIMENSION(its:ite) :: intp_mask     !GCC$ ATTRIBUTES aligned(64) :: intp_mask
      LOGICAL(kind=int4), DIMENSION(its:ite) :: tmask         !GCC$ ATTRIBUTES aligned(64) :: tmask
      INTEGER(kind=int4), DIMENSION(its:ite) :: kb            !GCC$ ATTRIBUTES aligned(64) :: kb
      INTEGER(kind=int4), DIMENSION(its:ite) :: kt            !GCC$ ATTRIBUTES aligned(64) :: kt
      REAL(kind=sp),      DIMENSION(its:ite) :: tl            !GCC$ ATTRIBUTES aligned(64) :: tl
      REAL(kind=sp),      DIMENSION(its:ite) :: tl2           !GCC$ ATTRIBUTES aligned(64) :: tl2
      REAL(kind=sp),      DIMENSION(its:ite) :: th            !GCC$ ATTRIBUTES aligned(64) :: th
      REAL(kind=sp),      DIMENSION(its:ite) :: th2           !GCC$ ATTRIBUTES aligned(64) :: th2
      REAL(kind=sp),      DIMENSION(its:ite) :: qqd           !GCC$ ATTRIBUTES aligned(64) :: qqd
      REAL(kind=sp),      DIMENSION(its:ite) :: qqh           !GCC$ ATTRIBUTES aligned(64) :: qqh
      REAL(kind=sp),      DIMENSION(its:ite) :: qql           !GCC$ ATTRIBUTES aligned(64) :: qql
      REAL(kind=sp),      DIMENSION(its:ite) :: zsum          !GCC$ ATTRIBUTES aligned(64) :: zsum
      REAL(kind=sp),      DIMENSION(its:ite) :: qsum          !GCC$ ATTRIBUTES aligned(64) :: qsum
      REAL(kind=sp),      DIMENSION(its:ite) :: dql           !GCC$ ATTRIBUTES aligned(64) :: dql
      REAL(kind=sp),      DIMENSION(its:ite) :: dqh           !GCC$ ATTRIBUTES aligned(64) :: dqh
      REAL(kind=sp),      DIMENSION(its:ite) :: za_gath_t     !GCC$ ATTRIBUTES aligned(64) :: za_gath_t 
      REAL(kind=sp),      DIMENSION(its:ite) :: za_gath_b     !GCC$ ATTRIBUTES aligned(64) :: za_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: qa_gath_b     !GCC$ ATTRIBUTES aligned(64) :: qa_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: dza_gath_t    !GCC$ ATTRIBUTES aligned(64) :: dza_gath_t
      REAL(kind=sp),      DIMENSION(its:ite) :: dza_gath_b    !GCC$ ATTRIBUTES aligned(64) :: dza_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: qpi_gath_t    !GCC$ ATTRIBUTES aligned(64) :: qpi_gath_t
      REAL(kind=sp),      DIMENSION(its:ite) :: qpi_gath_b    !GCC$ ATTRIBUTES aligned(64) :: qpi_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: qmi_gath_t    !GCC$ ATTRIBUTES aligned(64) :: qmi_gath_t
      REAL(kind=sp),      DIMENSION(its:ite) :: qmi_gath_b    !GCC$ ATTRIBUTES aligned(64) :: qmi_gath_b
#elif defined __INTEL_COMPILER
      LOGICAL(kind=int4), DIMENSION(its:ite) :: intp_mask, tmask
      !DIR$ ATTRIBUTES ALIGN : 64 :: intp_mask, tmask
      INTEGER(kind=int4), DIMENSION(its:ite) :: kb, kt
      !DIR$ ATTRIBUTES ALIGN : 64 :: kb, kt
      REAL(kind=sp),      DIMENSION(its:ite) :: tl,tl2,th,th2,qqd,qqh,qql,zsum,qsum,dql,dqh
      !DIR$ ATTRIBUTES ALIGN : 64 ::  tl,tl2,th,th2,qqd,qqh,qql,zsum,qsum,dql,dqh
      REAL(kind=sp),      DIMENSION(its:ite) :: za_gath_t,za_gath_b
      !DIR$ ATTRIBUTES ALIGN : 64 ::   za_gath_t,za_gath_b 
      REAL(kind=sp),      DIMENSION(its:ite) :: qa_gath_b
      !DIR$ ATTRIBUTES ALIGN : 64 ::   qa_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: dza_gath_t,dza_gath_b
      !DIR$ ATTRIBUTES ALIGN : 64 ::   dza_gath_t,dza_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: qpi_gath_t,qpi_gath_b
      !DIR$ ATTRIBUTES ALIGN : 64 ::   qpi_gath_t,qpi_gath_b
      REAL(kind=sp),      DIMENSION(its:ite) :: qmi_gath_t,qmi_gath_b
      !DIR$ ATTRIBUTES ALIGN : 64 :: qmi_gath_t,qmi_gath_b
#endif
!
#if defined __INTEL_COMPILER
!DIR$ ASSUME_ALIGNED denl:64,denfacl:64,tkl:64,dzl:64,wwl:64,rql:64,precip:64,lmask:64
!DIR$ ASSUME_ALIGNED intp_mask:64,tmask:64,kb:64,kt:64,tl:64,tl2:64,qqd:64,qqh:64
!DIR$ ASSUME_ALIGNED qql:64,zsum:64,dql:64,dqh:64,za_gath_t:64,za_gath_b:64,qa_gath_b:64
!DIR$ ASSUME_ALIGNED dza_gath_t:64,dza_gath_b:64,qpi_gath_t:64,qpi_gath_b:64
!DIR$ ASSUME_ALIGNED qmi_gath_t:64,qmi_gath_b:64,wi:64,ww:64,za:64,zi:64,dza:64
!DIR$ ASSUME_ALIGNED precip1:64,precip2:64,qmi:64,qa:64,rql2:64
#endif
      precip(:) = 0.0
      precip1(:) = 0.0
      precip2(:) = 0.0
!
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif 
      do i=its,ite
! -----------------------------------
      dz(i,k) = dzl(i,k)
      qq(i,k) = rql(i,k)
      qq2(i,k) = rql2(i,k)
      ww(i,k) = wwl(i,k)
      den(i,k) = denl(i,k)
      denfac(i,k) = denfacl(i,k)
      tk(i,k) = tkl(i,k)
      enddo
      enddo
! skip for no precipitation for all layers
      do i=its,ite
      allold(i) = 0.0
      enddo
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif          
      do i=its,ite
        allold(i) = allold(i) + qq(i,k)
      enddo
      enddo
      if (maxval(allold).le.0.0) return
      lmask = allold .gt. 0.0
!
! compute interface values
      do i=its,ite
      if(lmask(i)) then
      zi(i,kts)=0.0
      endif
      enddo
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif          
      do i=its,ite
      if(lmask(i)) then
        zi(i,k+1) = zi(i,k)+dz(i,k)
      endif
      enddo
      enddo
!
! save departure wind
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
    
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      
#endif          
      do i=its,ite
      if(lmask(i)) then
      wd(i,k) = ww(i,k)
      endif
      enddo
      enddo
      n=1
      do while (n.le.(iter+1))
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif          
      do i=its,ite
      if(lmask(i)) then
! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! 2nd order interpolation to get wi
      wi(i,kts) = ww(i,kts)
      wi(i,kte+1) = ww(i,kte)
      endif
      enddo
      do k=kts+1,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif          
      do i=its,ite
      if(lmask(i)) then
        wi(i,k) = (ww(i,k)*dz(i,k-1)+ww(i,k-1)*dz(i,k))/(dz(i,k-1)+dz(i,k))
      endif
      enddo
      enddo
! 3rd order interpolation to get wi
      fa1 = 9./16.
      fa2 = 1./16.
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif       
      do i=its,ite
      if(lmask(i)) then
      wi(i,kts) = ww(i,kts)
      wi(i,kts+1) = 0.5*(ww(i,kts+1)+ww(i,kts))
      endif
      enddo
      do k=kts+2,kte-1
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif          
      do i=its,ite
      if(lmask(i)) then
        wi(i,k) = fa1*(ww(i,k)+ww(i,k-1))-fa2*(ww(i,k+1)+ww(i,k-2))
      endif
      enddo
      enddo
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif    
      do i=its,ite
      if(lmask(i)) then
      wi(i,kte) = 0.5*(ww(i,kte)+ww(i,kte-1))
      wi(i,kte+1) = ww(i,kte)
      endif
      enddo
!
! terminate of top of raingroup
      do k=kts+1,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif          
      do i=its,ite
      if(lmask(i)) then
        if( ww(i,k).eq.0.0 ) wi(i,k)=ww(i,k-1)
      endif
      enddo
      enddo
!
! diffusivity of wi
      con1 = 0.05
      do k=kte,kts,-1
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif          
      do i=its,ite
      if(lmask(i)) then
        decfl = (wi(i,k+1)-wi(i,k))*dt/dz(i,k)
        if( decfl .gt. con1 ) then
          wi(i,k) = wi(i,k+1) - con1*dz(i,k)/dt
        endif
      endif
      enddo
      enddo
! compute arrival point
      do k=kts,kte+1
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif          
      do i=its,ite
      if(lmask(i)) then
        za(i,k) = zi(i,k) - wi(i,k)*dt
      endif
      enddo
      enddo
!
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif
      do i=its,ite
      if(lmask(i)) then
        dza(i,k) = za(i,k+1)-za(i,k)
      endif
      enddo
      enddo
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif
      do i=its,ite
      if(lmask(i)) then
      dza(i,kte+1) = zi(i,kte+1) - za(i,kte+1)
      endif
      enddo
!
! computer deformation at arrival point
      do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif
      do i=its,ite
      if(lmask(i)) then
        qa(i,k) = qq(i,k)*dz(i,k)/dza(i,k)
        qa2(i,k) = qq2(i,k)*dz(i,k)/dza(i,k)
        qr(i,k) = qa(i,k)/den(i,k)
        qr2(i,k) = qa2(i,k)/den(i,k)
      endif
      enddo
      enddo
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif   
      do i=its,ite
      if(lmask(i)) then
      qa(i,kte+1) = 0.0
      qa2(i,kte+1) = 0.0
!     call maxmin(kte-kts+1,1,qa(i,:),' arrival points ')
      endif
      enddo
!
! compute arrival terminal velocity, and estimate mean terminal velocity
! then back to use mean terminal velocity
      if( n.le.iter ) then
        call slope_snow_ii(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,its,ite,kts,kte,lmask)
        call slope_graup_ii(qr2,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa2,its,ite,kts,kte,lmask)
        do k = kts,kte
        do i=its,ite
        if(lmask(i)) then
          tmp(i,k) = max((qr(i,k)+qr2(i,k)), 1.E-15_sp)
          IF ( tmp(i,k) .gt. 1.e-15_sp ) THEN
            wa(i,k) = (wa(i,k)*qr(i,k) + wa2(i,k)*qr2(i,k))/tmp(i,k)
          ELSE
            wa(i,k) = 0.
          ENDIF
        endif
        enddo
        enddo
        if( n.ge.2 ) then
        do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif              
        do i=its,ite
        if(lmask(i)) then
          wa(i,k)=0.5*(wa(i,k)+was(i,k))
        endif
        enddo
        enddo
        endif
        do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif
        do i=its,ite
        if(lmask(i)) then
!#ifdef DEBUG
!        print*,' slope_wsm3 ',qr(i,k)*1000.,den(i,k),denfac(i,k),tk(i,k),tmp(i,k),tmp1(i,k),tmp2(i,k), &
!           ww(i,k),wa(i,k)
!#endif
! mean wind is average of departure and new arrival winds
          ww(i,k) = 0.5* ( wd(i,k)+wa(i,k) )
          was(i,k) = wa(i,k)
        endif
        enddo
        enddo
      endif
      n=n+1
      enddo
      ist_loop : do ist = 1, 2
      if (ist.eq.2) then
      do k=kts,kte+1
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif
      do i=its,ite
      if(lmask(i)) then
       qa(i,k) = qa2(i,k)
      endif
      enddo
      enddo
      endif
!
      do i=its,ite
      if(lmask(i)) then
      precip(i) = 0.
      endif
      enddo
!
! estimate values at arrival cell interface with monotone
      do k=kts+1,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif
      do i=its,ite
      if(lmask(i)) then
        dip=(qa(i,k+1)-qa(i,k))/(dza(i,k+1)+dza(i,k))
        dim=(qa(i,k)-qa(i,k-1))/(dza(i,k-1)+dza(i,k))
        if( dip*dim.le.0.0 ) then
          qmi(i,k)=qa(i,k)
          qpi(i,k)=qa(i,k)
        else
          qpi(i,k)=qa(i,k)+0.5*(dip+dim)*dza(i,k)
          qmi(i,k)=2.0*qa(i,k)-qpi(i,k)
          if( qpi(i,k).lt.0.0 .or. qmi(i,k).lt.0.0 ) then
            qpi(i,k) = qa(i,k)
            qmi(i,k) = qa(i,k)
          endif
        endif
      endif
      enddo
      enddo
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif
      do i=its,ite
      if(lmask(i)) then
      qpi(i,kts)=qa(i,kts)
      qmi(i,kts)=qa(i,kts)
      qmi(i,kte+1)=qa(i,kte+1)
      qpi(i,kte+1)=qa(i,kte+1)
      endif
      enddo
!
! interpolation to regular point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGIN WSM5 CODE FROM John Michalakes !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      qn = 0.0
      kb=kts  ! kb is a vector
      kt=kts  ! kt is a vector
      INTP : do k=kts,kte
             kb=max(kb-1,kts)
             kt=max(kt-1,kts)
! find kb and kt
             intp_mask = ( zi(:,k).lt.za(:,kte+1) .AND. lmask )
             tmask = intp_mask
             minkb = 999
             minkt = 999
             DO i=its,ite
               IF ( tmask(i) .AND. kb(i) .lt. minkb ) minkb = kb(i)
               IF ( tmask(i) .AND. kt(i) .lt. minkt ) minkt = kt(i)
             ENDDO
             find_kb : do kk=minkb,kte
               DO i=its,ite
               IF ( tmask(i) .AND. zi(i,k).le.za(i,kk+1) ) THEN
                 kb(i) = kk
                 tmask(i) = .FALSE.
               ENDIF
               ENDDO
             enddo find_kb

             tmask = intp_mask
             find_kt : do kk=minkt,kte
               DO i=its,ite
               IF ( tmask(i) .AND. zi(i,k+1).le.za(i,kk) ) THEN
                 kt(i) = kk
                 tmask(i) = .FALSE.
               ENDIF
               ENDDO
             enddo find_kt
             kt = max(kt - 1,kts)

!#define RANGE_CHECKING
#ifndef RANGE_CHECKING
# define DX1 (i+(kb(i)-1)*(ite-its+1)),1
# define DX2 (i+(kt(i)-1)*(ite-its+1)),1
#else
# define DX1 i,kb(i)
# define DX2 i,kt(i)
#endif
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif
             DO i = its,ite
               qa_gath_b(i) = qa(DX1)
               za_gath_b(i) = za(DX1)
               dza_gath_b(i) = dza(DX1)
               qpi_gath_b(i) = qpi(DX1)
               qmi_gath_b(i) = qmi(DX1)
             ENDDO
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
    
#endif
             DO i = its,ite
               za_gath_t(i) = za(DX2)
               dza_gath_t(i) = dza(DX2)
               qpi_gath_t(i) = qpi(DX2)
               qmi_gath_t(i) = qmi(DX2)
             ENDDO
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
    
#endif
             DO i = its,ite
             IF ( kt(i) .eq. kb(i) .AND. intp_mask(i) ) THEN
               tl(i)=(zi(i,k)-za_gath_b(i))/dza_gath_b(i)
               th(i)=(zi(i,k+1)-za_gath_b(i))/dza_gath_b(i)
               tl2(i) = tl(i)*tl(i)
               th2(i) = th(i)*th(i)
               qqd(i)=0.5*(qpi_gath_b(i)-qmi_gath_b(i))
               qqh(i)=qqd(i)*th2(i)+qmi_gath_b(i)*th(i)
               qql(i)=qqd(i)*tl2(i)+qmi_gath_b(i)*tl(i)
               qn(i,k) = (qqh(i)-qql(i))/(th(i)-tl(i))
             ELSE IF ( kt(i) .gt. kb(i) .AND. intp_mask(i) ) THEN
               tl(i)=(zi(i,k)-za_gath_b(i))/dza_gath_b(i)
               tl2(i)=tl(i)*tl(i)
               qqd(i)=0.5*(qpi_gath_b(i)-qmi_gath_b(i))
               qql(i)=qqd(i)*tl2(i)+qmi_gath_b(i)*tl(i)
               dql(i) = qa_gath_b(i)-qql(i)
               zsum(i)  = (1.-tl(i))*dza_gath_b(i)
               qsum(i)  = dql(i)*dza_gath_b(i)
             ENDIF
             ENDDO
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif          
             DO i = its,ite
               if( kt(i)-kb(i).gt.1 .AND. intp_mask(i) ) then
                 do m=kb(i)+1,kt(i)-1
                     zsum(i) = zsum(i) + dza(i,m)
                     qsum(i) = qsum(i) + qa(i,m) * dza(i,m)
                 enddo
               endif
            ENDDO
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif            
             DO i = its,ite
             IF ( kt(i) .gt. kb(i) .AND. intp_mask(i) ) THEN
               th(i)=(zi(i,k+1)-za_gath_t(i))/dza_gath_t(i)
               th2(i) = th(i)*th(i)
               qqd(i)=0.5*(qpi_gath_t(i)-qmi_gath_t(i))
               dqh(i)=qqd(i)*th2(i)+qmi_gath_t(i)*th(i)
               zsum(i)  = zsum(i) + th(i)*dza_gath_t(i)
               qsum(i)  = qsum(i) + dqh(i)*dza_gath_t(i)
               qn(i,k) = qsum(i)/zsum(i)
             ENDIF
             ENDDO
       ENDDO intp
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END WSM5 CODE FROM John Michalakes !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! rain out
      intp_mask = lmask
      sum_precip: do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
      !DIR$ IVDEP
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      !GCC$ IVDEP
#endif         
             DO i = its,ite
             IF (za(i,k).lt.0.0.and.za(i,k+1).lt.0.0.AND.intp_mask(i)) THEN
               precip(i) = precip(i) + qa(i,k)*dza(i,k)
             ELSE IF (za(i,k).lt.0.0.and.za(i,k+1).ge.0.0.AND.intp_mask(i)) THEN
               precip(i) = precip(i) + qa(i,k)*(0.0-za(i,k))
               intp_mask(i) = .FALSE.
             ENDIF
             ENDDO
      enddo sum_precip
!
! replace the new values
      if(ist.eq.1) then
        do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
    
#elif defined __GFORTRAN__
      !GCC$ VECTOR
    
#endif
        do i=its,ite
        if(lmask(i)) then
        rql(i,k) = qn(i,k)
        endif
        enddo
        enddo
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
    
#endif
        do i=its,ite
        if(lmask(i)) then
        precip1(i) = precip(i)
        endif
        enddo
      else
        do k=kts,kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
      
#endif
        do i=its,ite
        if(lmask(i)) then
        rql2(i,k) = qn(i,k)
        endif
        enddo
        enddo
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif
        do i=its,ite
        if(lmask(i)) then
        precip2(i) = precip(i)
        endif
        enddo
      endif
      enddo ist_loop
!
! ----------------------------------
!
  END SUBROUTINE nislfv_rain_plm6_ii
#endif
!+---+-----------------------------------------------------------------+
  ! Read array from unitno and convert from "no-chunk" to "chunk"
#if defined __GFORTRAN__
  SUBROUTINE readarray2(arr,arrname,unitno,ips,ipe) !GCC$ ATTRIBUTES cold :: readarray2 !GCC$ ATTRIBUTES aligned(32) :: readarray2
#elif defined __INTEL_COMPILER
    SUBROUTINE readarray2(arr,arrname,unitno,ips,ipe)
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: readarray2
#endif
    REAL(kind=sp), contiguous,   INTENT(OUT) :: arr(:,:)
    CHARACTER(LEN=*),            INTENT(IN)  :: arrname
    INTEGER(kind=int4),          INTENT(IN)  :: unitno
    INTEGER(kind=int4),          INTENT(IN)  :: ips,ipe
    REAL(kind=sp), ALLOCATABLE :: tmparr(:)
    INTEGER(kind=int4) :: i,j,ij,jsize,CHUNK,ipn
    CHUNK = size(arr,1)
    jsize = size(arr,2)
    ALLOCATE(tmparr(ips:ipe))
    read(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do i = 1,CHUNK
        ! replicate last column if needed
        ipn = min(ij+i+ips-1,ipe)
        arr(i,j) = tmparr(ipn)
      enddo
    enddo
    DEALLOCATE(tmparr)
  END SUBROUTINE readarray2

!+---+-----------------------------------------------------------------+
  ! Read array from unitno and convert from "no-chunk" to "chunk"
#if defined __GFORTRAN__
  SUBROUTINE readarray3(arr,arrname,unitno,ips,ipe) !GCC$ ATTRIBUTES cold :: readarray3 !GCC$ ATTRIBUTES aligned(32) :: readarray3
#elif defined __INTEL_COMPILER
    SUBROUTINE readarray3(arr,arrname,unitno,ips,ipe)
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: readarray3
#endif
    REAL(kind=sp),  contiguous   INTENT(OUT) :: arr(:,:,:)
    CHARACTER(LEN=*),            INTENT(IN)  :: arrname
    INTEGER(kind=int4),          INTENT(IN)  :: unitno
    INTEGER(kind=int4),          INTENT(IN)  :: ips,ipe
    REAL(kind=sp), ALLOCATABLE :: tmparr(:,:)
    INTEGER(kind=int4) :: i,j,k,ij,ksize,jsize,CHUNK,ipn
    CHUNK = size(arr,1)
    ksize = size(arr,2)
    jsize = size(arr,3)
    ALLOCATE(tmparr(ips:ipe,ksize))
    read(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do i = 1,CHUNK
        ! replicate last column if needed
        ipn = min(ij+i+ips-1,ipe)
        do k = 1,ksize
          arr(i,k,j) = tmparr(ipn,k)
        enddo
      enddo
    enddo
    DEALLOCATE(tmparr)
  END SUBROUTINE readarray3

!+---+-----------------------------------------------------------------+
  ! Read array from unitno and convert from "no-chunk" to "chunk"
  SUBROUTINE readarray4(arr,arrname,unitno,ips,ipe)
    REAL(kind=sp),  contiguous           INTENT(OUT) :: arr(:,:,:,:)
    CHARACTER(LEN=*),                    INTENT(IN)  :: arrname
    INTEGER(kind=int4),                  INTENT(IN)  :: unitno
    INTEGER(kind=int4),                  INTENT(IN)  :: ips,ipe
    REAL(kind=sp), ALLOCATABLE :: tmparr(:,:,:)
    INTEGER(kind=int4) :: i,j,k,ij,ksize,jsize,CHUNK,ipn,m,msize
    CHUNK = size(arr,1)
    ksize = size(arr,2)
    msize = size(arr,3)
    jsize = size(arr,4)
    ALLOCATE(tmparr(ips:ipe,ksize,msize))
    read(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do m = 1,msize
        do i = 1,CHUNK
          ! replicate last column if needed
          ipn = min(ij+i+ips-1,ipe)
          do k = 1,ksize
            arr(i,k,m,j) = tmparr(ipn,k,m)
          enddo
        enddo
      enddo
    enddo
    DEALLOCATE(tmparr)
  END SUBROUTINE readarray4

!+---+-----------------------------------------------------------------+
  ! Convert array from "chunk" to "no-chunk" and write to unitno.  
  SUBROUTINE writearray2(arr,arrname,unitno,ips,ipe)
    REAL(kind=sp),    contiguous,           INTENT(IN) :: arr(:,:)
    CHARACTER(LEN=*),                       INTENT(IN) :: arrname
    INTEGER(kind=int4),                     INTENT(IN) :: unitno
    INTEGER(kind=int4) ,                    INTENT(IN) :: ips,ipe
    REAL(kind=sp), ALLOCATABLE :: tmparr(:)
    INTEGER(kind=sp) :: i,j,ij,jsize,CHUNK,ipn
    CHUNK = size(arr,1)
    jsize = size(arr,2)
    ALLOCATE(tmparr(ips:ipe))
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do i = 1,CHUNK
        ipn = ij+i+ips-1
        ! skip any replicated columns
        if ((ips<=ipn).and.(ipn<=ipe)) then
          tmparr(ipn) = arr(i,j)
        endif
      enddo
    enddo
    write(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    DEALLOCATE(tmparr)
  END SUBROUTINE writearray2

!+---+-----------------------------------------------------------------+
  ! Convert array from "chunk" to "no-chunk" and write to unitno.  
  SUBROUTINE writearray3(arr,arrname,unitno,ips,ipe)
    REAL(kind=sp),    contiguous,         INTENT(IN) :: arr(:,:,:)
    CHARACTER(LEN=*),                     INTENT(IN) :: arrname
    INTEGER(kind=int4),                   INTENT(IN) :: unitno
    INTEGER(kind=int4),                   INTENT(IN) :: ips,ipe
    REAL(kind=sp), ALLOCATABLE :: tmparr(:,:)
    INTEGER(kind=int4) :: i,j,k,ij,ksize,jsize,CHUNK,ipn
    CHUNK = size(arr,1)
    ksize = size(arr,2)
    jsize = size(arr,3)
    ALLOCATE(tmparr(ips:ipe,ksize))
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do i = 1,CHUNK
        ipn = ij+i+ips-1
        ! skip any replicated columns
        if ((ips<=ipn).and.(ipn<=ipe)) then
          do k = 1,ksize
            tmparr(ipn,k) = arr(i,k,j)
          enddo
        endif
      enddo
    enddo
    write(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    DEALLOCATE(tmparr)
  END SUBROUTINE writearray3

!+---+-----------------------------------------------------------------+
  ! Convert array from "chunk" to "no-chunk" and write to unitno.  
  SUBROUTINE writearray4(arr,arrname,unitno,ips,ipe)
    REAL(kind=sp),          contiguous,   INTENT(IN) :: arr(:,:,:,:)
    CHARACTER(LEN=*),                     INTENT(IN) :: arrname
    INTEGER(kind=int4),                   INTENT(IN) :: unitno
    INTEGER(kind=int4),                   INTENT(IN) :: ips,ipe
    REAL(kind=sp), ALLOCATABLE :: tmparr(:,:,:)
    INTEGER :: i,j,k,ij,ksize,jsize,CHUNK,ipn,m,msize
    CHUNK = size(arr,1)
    ksize = size(arr,2)
    msize = size(arr,3)
    jsize = size(arr,4)
    ALLOCATE(tmparr(ips:ipe,ksize,msize))
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do m = 1,msize
        do i = 1,CHUNK
          ipn = ij+i+ips-1
          ! skip any replicated columns
          if ((ips<=ipn).and.(ipn<=ipe)) then
            do k = 1,ksize
              tmparr(ipn,k,m) = arr(i,k,m,j)
            enddo
          endif
        enddo
      enddo
    enddo
    write(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    DEALLOCATE(tmparr)
  END SUBROUTINE writearray4
#if defined __GFORTRAN__
  SUBROUTINE firstTouch(t, q                                      &   
                   ,qci, qrs, den, p, delz                        &
                   ,lat                                           &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte                     &
                   ,snow,snowncv                                  &
                   ,graupel,graupelncv                            &
                   )  !GCC$ ATTRIBUTES cold :: firstTouch !GCC$ ATTRIBUTES aligned(32) :: firstTouch
#elif defined __INTEL_COMPILER
    SUBROUTINE firstTouch(t, q                                      &   
                   ,qci, qrs, den, p, delz                        &
                   ,lat                                           &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte                     &
                   ,snow,snowncv                                  &
                   ,graupel,graupelncv                            &
                   )
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: firstTouch
#endif
!-------------------------------------------------------------------
 
!-------------------------------------------------------------------
!
!  Mimic memory access patterns of wsm62D() while setting physics arrays 
!  to zero.  
!
  INTEGER(kind=int4),      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte,  &
                                      lat
  REAL(kind=sp), DIMENSION( its:ite , kts:kte ),                           &
        INTENT(INOUT) ::                                          &
                                                               t
  REAL(kind=sp), DIMENSION( its:ite , kts:kte, 2 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qci
  REAL(kind=sp), DIMENSION( its:ite , kts:kte, 3 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qrs
  REAL(kind=sp), DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                               q
  REAL(kind=sp), DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                             den, &
                                                               p, &
                                                            delz
  REAL(kind=sp), DIMENSION( ims:ime ),                                     &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr
  REAL(kind=sp), DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv
  REAL(kind=sp), DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                 graupel, &
                                                      graupelncv
! LOCAL VAR
  INTEGER :: i, k
!
  do k = kts, kte
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
    
#endif     
     do i = its, ite
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED t:64,q:64,qci:64,qrs:64,den:64,p:64,delz:64
#endif
          t(i,k) = 0.
          q(i,k) = 0.
          qci(i,k,1) = 0.
          qci(i,k,2) = 0.
          qrs(i,k,1) = 0.
          qrs(i,k,2) = 0.
          qrs(i,k,3) = 0.
          den(i,k) = 0.
          p(i,k) = 0.
          delz(i,k) = 0.
        enddo
     enddo
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
    
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif     
     do i = its, ite
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED rain:64,rainncv:64,sr:64
#endif
        rain(i) = 0.
        rainncv(i) = 0.
        sr(i) = 0.
      enddo
      if (PRESENT(snow)) then
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED snow:64
#endif            
          snow(i,lat) = 0.
        enddo
      endif
      if (PRESENT(snowncv)) then
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
   
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif         
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED snowncv:64
#endif
          snowncv(i,lat) = 0.
        enddo
      endif
      if (PRESENT(graupel)) then
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
     
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif         
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED graupel:64
#endif
          graupel(i,lat) = 0.
        enddo
      endif
      if (PRESENT(graupelncv)) then
#if defined __INTEL_COMPILER
      !DIR$ VECTOR ALIGNED
      !DIR$ VECTOR ALWAYS
    
#elif defined __GFORTRAN__
      !GCC$ VECTOR
     
#endif         
         do i = its, ite
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED graupelncv:64
#endif
          graupelncv(i,lat) = 0.
        enddo
      endif
  END SUBROUTINE firstTouch

END MODULE module_mp_wsm6




 ! PROGRAM wsm6
!-------------------------------------------------------------------
 ! USE module_mp_wsm6
!-------------------------------------------------------------------
!  IMPLICIT NONE
!-------------------------------------------------------------------
 ! INTEGER ::   ids, ide,  jds, jde, kds,kde , &
 !              ims, ime,  jms, jme, kms,kme , &
 !              its, ite,  jts, jte, kts,kte
  ! "chunk" indices
!  INTEGER ::  iids,iide, jjds,jjde,           &
!              iims,iime, jjms,jjme,           &
 !             iits,iite, jjts,jjte
! timers
!#include <gptl.inc>
!  INTEGER :: ret
!  REAL*8  :: totaltime
!  INTEGER :: handle = 0
! Contains definition of _CHUNK_
!#include <chunk.h>
! Contains definition of _NZ_
!#include <literalk.h>
 ! INTEGER :: CHUNK, NZ_BUILD
 ! INTEGER :: num_tiles_C
 ! INTEGER ::  itimestep
 ! REAL, ALLOCATABLE :: t(:,:,:)
 ! REAL, ALLOCATABLE :: qci(:,:,:,:)
 ! REAL, ALLOCATABLE :: qrs(:,:,:,:)
 ! REAL, ALLOCATABLE :: q(:,:,:)
 ! REAL, ALLOCATABLE :: den(:,:,:)
 ! REAL, ALLOCATABLE :: p(:,:,:)
 ! REAL, ALLOCATABLE :: delz(:,:,:)
!  REAL :: delt, g, rd, rv, t0c, den0, cpd, cpv, ep1, ep2, qmin, &
 !         XLS, XLV0, XLF0, cliq, cice, psat, denr
 ! REAL, ALLOCATABLE :: rain(:,:)
 ! REAL, ALLOCATABLE :: rainncv(:,:)
 ! REAL, ALLOCATABLE :: sr(:,:)
 ! REAL, ALLOCATABLE :: snow(:,:)
 ! REAL, ALLOCATABLE :: snowncv(:,:)
 ! REAL, ALLOCATABLE :: graupel(:,:)
 ! REAL, ALLOCATABLE :: graupelncv(:,:)

!+---+-----------------------------------------------------------------+
! LOCAL VAR
 ! INTEGER ::               i,j,k
 ! INTEGER :: ios, unitno
 ! CHARACTER(LEN=64) :: fn  ! file name

!+---+-----------------------------------------------------------------+
! EXTERNAL FUNCTIONS
!#ifdef _OPENMP
!  integer, external :: omp_get_max_threads
!#endif
!+---+-----------------------------------------------------------------+

 ! call gptlprocess_namelist ('GPTLnamelist', 77, ret)
 ! ret = gptlinitialize ()
 ! ret = gptlstart('Total')

  ! read constants
 ! PRINT *,'wsm6init():  read constants'
 ! fn = 'wsm6_constants.dat'
 ! unitno=31
 ! open (unitno,file=trim(fn),form="unformatted",action='read', &
 !       iostat=ios)
 ! if (ios /= 0) then
 !   write(6,*) 'ERROR: failed to open constant file ',trim(fn), &
 !              ' . stopping'
 !   stop
!  endif
!  read(unitno) pi, xlv1
!  read(unitno) qc0, qck1
 ! read(unitno) bvtr1, bvtr2, bvtr3, bvtr4, bvtr6, g1pbr, g3pbr, &
 !              g4pbr, g6pbr, g5pbro2, pvtr, eacrr, pacrr, &
 !              precr1, precr2, roqimax
 ! read(unitno) bvts1, bvts2, bvts3, bvts4, g1pbs, g3pbs, g4pbs,  &
 !              g5pbso2, pvts, pacrs, precs1, precs2, pidn0r, pidn0s
 ! read(unitno) pacrc
 ! read(unitno) bvtg1, bvtg2, bvtg3, bvtg4, g1pbg, g3pbg, g4pbg,  &
 !              pacrg, g5pbgo2, pvtg, precg1, precg2, pidn0g
 !! read(unitno) rslopermax, rslopesmax, rslopegmax, rsloperbmax,  &
 !              rslopesbmax, rslopegbmax, rsloper2max, rslopes2max,  &
 !              rslopeg2max, rsloper3max, rslopes3max, rslopeg3max
 ! close(unitno)

  ! read input data
  !PRINT *,'wsm62D():  read input state'
  !fn = 'wsm6_input.dat'
 ! unitno=31
 ! open (unitno,file=trim(fn),form="unformatted",action='read', &
 !       iostat=ios)
 ! if (ios /= 0) then
 !   write(6,*) 'ERROR: failed to open input file ',trim(fn), &
 !              ' . stopping'
 !   stop
 ! endif
 ! read(unitno) itimestep
 ! PRINT *,'wsm62D():  itimestep == ',itimestep
  ! read serial versions of indices
 ! read(unitno) ids,ide, jds,jde, kds,kde, &
 !              ims,ime, jms,jme, kms,kme, &
 !              its,ite, jts,jte, kts,kte

  ! assert constraints on indices
  ! NOTE that [ikj]d[se] are ignored
 ! if ((ims/=its).or.(ime/=ite).or.(jms/=jts).or.(jme/=jte)) then
 !   print *,'ERROR:  index mismatch reading file ',trim(fn), &
 !           ' . stopping'
 !   stop
!  endif
 ! if ((ims/=1).or.(jms/=1)) then
 !   print *,'ERROR:  incorrect start index reading file ',trim(fn), &
 !           ' . stopping'
 !   stop
 ! endif

  ! set default values of "chunk" indices
 ! iids = ids
 ! iide = ide
 ! iims = ims
 ! iime = ime
 ! iits = its
 ! iite = ite
 ! jjds = jds
 ! jjde = jde
 ! jjms = jms
 ! jjme = jme
 ! jjts = jts
 ! jjte = jte

  ! set up optional "i" chunking and optional fixed vertical extent
!#ifdef _CHUNK_
 !  CHUNK = _CHUNK_
!#else
!   CHUNK = iite-iits+1
!#endif
!#ifdef _NZ_
!   NZ_BUILD = _NZ_
   ! if specified, NZ_BUILD must match namelist setting for nvl
 !  if (NZ_BUILD/=kte-kts+1) then
 !    print *, 'ERROR:  Build-time-specified NZ must equal namelist nz, values are: ',NZ_BUILD,kte-kts+1
 !    call flush(6)
!     stop
! !  endif
!#else
!   NZ_BUILD = kte-kts+1
!#endif
 !  num_tiles_C = (iite-iits+1) / CHUNK
 !  if (mod((iite-iits+1),CHUNK) > 0) then
 !    num_tiles_C = num_tiles_C + 1
  ! endif
  ! iime = CHUNK
 !  iite = CHUNK
 !  jjme = num_tiles_C
 !  jjte = num_tiles_C
 !  PRINT *,'ims,ime,iims,iime',ims,ime,iims,iime
!#ifdef _OPENMP
!   PRINT "('omp_get_max_threads() returned:  ',I9)",omp_get_max_threads()
!#endif
!   PRINT *,'CHUNK =', CHUNK
!   PRINT *,'NUMBER OF CHUNKS =', num_tiles_C
!#ifdef _NZ_
!   PRINT *,'NZ_BUILD =', NZ_BUILD
!#endif

  ! allocate arrays
!  ALLOCATE(t(iits:iite,kts:kte,jjts:jjte))
!  ALLOCATE(qci(iits:iite,kts:kte,2,jjts:jjte))
!  ALLOCATE(qrs(iits:iite,kts:kte,3,jjts:jjte))
!  ALLOCATE(q(iims:iime,kms:kme,jjms:jjme))
!  ALLOCATE(den(iims:iime,kms:kme,jjms:jjme))
!  ALLOCATE(p(iims:iime,kms:kme,jjms:jjme))
!  ALLOCATE(delz(iims:iime,kms:kme,jjms:jjme))
!  ALLOCATE(rain(iims:iime,jjms:jjme))
!  ALLOCATE(rainncv(iims:iime,jjms:jjme))
!  ALLOCATE(sr(iims:iime,jjms:jjme))
!  ALLOCATE(snow(iims:iime,jjms:jjme))
 ! ALLOCATE(snowncv(iims:iime,jjms:jjme))
!  ALLOCATE(graupel(iims:iime,jjms:jjme))
 ! ALLOCATE(graupelncv(iims:iime,jjms:jjme))

!!$OMP PARALLEL DO &
!!$OMP PRIVATE ( j ) &
!!$OMP SCHEDULE(runtime)
!  do j = jjts,jjte
!    CALL firstTouch(t(iits,kts,j), q(iims,kms,j)              &
 !              ,qci(iits,kts,1,j), qrs(iits,kts,1,j)          &
 !              ,den(iims,kms,j)                               &
 !              ,p(iims,kms,j), delz(iims,kms,j)               &
 !              ,j                                             &
 !              ,rain(iims,j),rainncv(iims,j)                  &
 !              ,sr(iims,j)                                    &
 !              ,iids,iide, jjds,jjde, kds,kde                 &
 !              ,iims,iime, jjms,jjme, kms,kme                 &
 !              ,iits,iite, jjts,jjte, kts,kte                 &
 !              ,snow,snowncv                                  &
 !              ,graupel,graupelncv                            &
                                                            !  )
 ! enddo  ! j loop
!!$OMP END PARALLEL DO

  ! read remaining input data
!  call readarray3(t,'t',unitno,its,ite)
 ! call readarray4(qci,'qci',unitno,its,ite)
 ! call readarray4(qrs,'qrs',unitno,its,ite)
 ! call readarray3(q,'q',unitno,its,ite)
 ! call readarray3(den,'den',unitno,its,ite)
 ! call readarray3(p,'p',unitno,its,ite)
 ! call readarray3(delz,'delz',unitno,its,ite)
 ! read(unitno) delt
 ! print *,' delt = ',delt
 ! read(unitno) g
 ! print *,' g = ',g
 ! read(unitno) cpd
 ! print *,' cpd = ',cpd
 ! read(unitno) cpv
 ! print *,' cpv = ',cpv
 ! read(unitno) t0c
!  print *,' t0c = ',t0c
 ! read(unitno) den0
 ! print *,' den0 = ',den0
 ! read(unitno) rd
 ! print *,' rd = ',rd
 ! read(unitno) rv
!  print *,' rv = ',rv
!  read(unitno) ep1
!  print *,' ep1 = ',ep1
!  read(unitno) ep2
!  print *,' ep2 = ',ep2
!  read(unitno) qmin
!  print *,' qmin = ',qmin
!  read(unitno) XLS
 ! print *,' XLS = ',XLS
 ! read(unitno) XLV0
 ! print *,' XLV0 = ',XLV0
 ! read(unitno) XLF0
 ! print *,' XLF0 = ',XLF0
 ! read(unitno) cliq
!  print *,' cliq = ',cliq
 ! read(unitno) cice
!  print *,' cice = ',cice
!  read(unitno) psat
!  print *,' psat = ',psat
!  read(unitno) denr
!  print *,' denr = ',denr
!  call readarray2(rain,'rain',unitno,its,ite)
!  call readarray2(rainncv,'rainncv',unitno,its,ite)
!  call readarray2(sr,'sr',unitno,its,ite)
!  call readarray2(snow,'snow',unitno,its,ite)
!  call readarray2(snowncv,'snowncv',unitno,its,ite)
!  call readarray2(graupel,'graupel',unitno,its,ite)
!  call readarray2(graupelncv,'graupelncv',unitno,its,ite)
!  close(unitno)

  ! minimize timer overhead inside OpenMP loop
 ! ret = gptlinit_handle ('WSM62D', handle)
  ! call WSM6
 ! ret = gptlstart('WSM62D+OpenMP')
!!$OMP PARALLEL DO &
!!$OMP PRIVATE ( j,ret ) &
!!$OMP SCHEDULE(runtime)
!  do j = jjts,jjte
!    ret = gptlstart_handle ('WSM62D', handle)
!    CALL wsm62D(t(iits,kts,j), q(iims,kms,j)                  &
          !     ,qci(iits,kts,1,j), qrs(iits,kts,1,j)          &
          !     ,den(iims,kms,j)                               &
          !     ,p(iims,kms,j), delz(iims,kms,j)               &
          !     ,delt,g, cpd, cpv, rd, rv, t0c                 &
          !     ,ep1, ep2, qmin                                &
          !     ,XLS, XLV0, XLF0, den0, denr                   &
          !     ,cliq,cice,psat                                &
          !     ,j                                             &
         !      ,rain(iims,j),rainncv(iims,j)                  &
         !      ,sr(iims,j)                                    &
         !      ,iids,iide, jjds,jjde, kds,kde                 &
         !      ,iims,iime, jjms,jjme, kms,kme                 &
  !             ,iits,iite, jjts,jjte, kts,kte                 &
  !             ,snow,snowncv                                  &
  !             ,graupel,graupelncv                            &
  !!                                                            )
 !   ret = gptlstop_handle ('WSM62D', handle)
!  enddo  ! j loop
!!$OMP END PARALLEL DO
  ! ret = gptlstop('WSM62D+OpenMP')

  ! write output data
 ! PRINT *,'wsm62D():  itimestep == ',itimestep,', write output state'
  !fn = 'wsm6_output.dat'
 ! unitno=31
 ! open (unitno,file=trim(fn),form="unformatted",action='write', &
 !       iostat=ios)
 ! if (ios /= 0) then
 !   write(6,*) 'ERROR: failed to open output file ',trim(fn), &
 !              ' . stopping'
 !   stop
 ! endif
 ! call writearray3(t,'t',unitno,its,ite)
 ! call writearray4(qci,'qci',unitno,its,ite)
 ! call writearray4(qrs,'qrs',unitno,its,ite)
 ! call writearray3(q,'q',unitno,its,ite)
 ! call writearray2(rain,'rain',unitno,its,ite)
 ! call writearray2(rainncv,'rainncv',unitno,its,ite)
 ! call writearray2(sr,'sr',unitno,its,ite)
!  call writearray2(snow,'snow',unitno,its,ite)
 ! call writearray2(snowncv,'snowncv',unitno,its,ite)
 ! call writearray2(graupel,'graupel',unitno,its,ite)
 ! call writearray2(graupelncv,'graupelncv',unitno,its,ite)
 ! close(unitno)

  ! deallocate arrays
!  DEALLOCATE(t)
 ! DEALLOCATE(qci)
 ! DEALLOCATE(qrs)
!  DEALLOCATE(q)
!  DEALLOCATE(den)
 ! DEALLOCATE(p)
 ! DEALLOCATE(delz)
 ! DEALLOCATE(rain)
!  DEALLOCATE(rainncv)
!  DEALLOCATE(sr)
!  DEALLOCATE(snow)
!  DEALLOCATE(snowncv)
 ! DEALLOCATE(graupel)
 ! DEALLOCATE(graupelncv)

!  ret = gptlstop('Total')
 !! ret = gptlget_wallclock ('Total', 0, totaltime)  ! The "0" is thread number
!  print*,''
 ! print*,'Total time =' , totaltime

!+---+-----------------------------------------------------------------+

  ! print timing info
 ! ret = gptlpr (0)
 ! ret = gptlpr_summary (0)

 ! END PROGRAM wsm6

