
#include "GMS_config.fpp"

!MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.

module mod_wsm6_driver

    !===================================================================================85
 !---------------------------- DESCRIPTION ------------------------------------------85
 !
 !
 !
 !          Module  name:
 !                         'mod_wsm6_driver'
 !          
 !          Purpose:
 !                       This is a driver for --> 
 !                       6-class GRAUPEL phase microphyiscs scheme (WSM6) of the 
 !                       Single-Moment MicroPhyiscs (WSMMP)
 !          History:
 !                        Date: 26-05-2018
 !                        Time: 10:26 GMT+2
 !
 !          Version:
 !
 !                      Major: 1
 !                      Minor: 1
 !                      Micro: 0
 !
 !          Author:  
 !                   Song-You Hong and Jeong-Ock Jade Lim (Yonsei Univ.)
 !          Modified:
 !                   Bernard Gingold on 26-05-2018
 !                 
 !          References:
 !         
 !                
 !    
 !         
 !          E-mail:
 !                  
 !                      beniekg@gmail.com
 !==================================================================================85
    ! Tab:5 col - Type and etc.. definitions
    ! Tab:10,11 col - Type , function and subroutine code blocks.
    
    
    use mod_kinds,   only : i4,sp
#if defined(USE_PROFILING) && (SAMPLE_PMC) == 1
    use pmc_samples, only : core_counters1D, &
                            init_core_counters1D
    use fast_pmc_access
    use IFCORE, only : FOR_LFENCE
#endif  
    use module_mp_wsm6, only : firstTouch,wsm62D,
                               readarray2,readarray3,readarray4,     &
                               pi, xlv1,                             &
                               qc0, qck1,                            &
                               bvtr1, bvtr2, bvtr3,                  &
                               bvtr4, bvtr6, g1pbr, g3pbr,           &
                               g4pbr, g6pbr, g5pbro2, pvtr,          &
                               eacrr, pacrr,                         &
                               precr1, precr2, roqimax,              &
                               bvts1, bvts2, bvts3, bvts4,           &
                               g1pbs, g3pbs, g4pbs,                  &
                               g5pbso2, pvts, pacrs, precs1,         & 
                               precs2, pidn0r, pidn0s, pacrc,        &
                               bvtg1, bvtg2, bvtg3, bvtg4,           &
                               g1pbg, g3pbg, g4pbg,                  &
                               pacrg, g5pbgo2, pvtg, precg1,         &
                               precg2, pidn0g ,rslopermax,           &
                               rslopesmax, rslopegmax, rsloperbmax,  &
                               rslopesbmax, rslopegbmax,             &
                               rsloper2max, rslopes2max,             &
                               rslopeg2max, rsloper3max,             &
                               rslopes3max, rslopeg3max

    
    implicit none
    
    public :: wsm6D_driver
    !=====================================================59
    !  File and module information:
    !  version,creation and build date, author,description
    !=====================================================59
    
    ! Major version
    integer(kind=i4), parameter, public :: MOD_WSM6_DRIVER_MAJOR = 1
    
    ! Minor version
    integer(kind=i4), parameter, public :: MOD_WSM6_DRIVER_MINOR = 1
    
    ! Micro version
    integer(kind=i4), parameter, public :: MOD_WSM6_DRIVER_MICRO = 0
    
    ! Module full version
    integer(kind=i4), parameter, public :: MOD_WSM6_DRIVER_FULLVER = 1000*MOD_WSM6_DRIVER_MAJOR + &
                                                                  100*MOD_WSM6_DRIVER_MINOR  + &
                                                                  10*MOD_WSM6_DRIVER_MICRO
    ! Module creation date
    character(*),  parameter, public :: MOD_WSM6_DRIVER_CREATE_DATE = "26-05-2018 10:26 +00200 (SAT 26 MAY 2018 GMT+2)"
    
    ! Module build date (should be set after successful compilation)
    character(*),  parameter, public :: MOD_WSM6_DRIVER_BUILD_DATE = __DATE__ " " __TIME__
    
    ! Module author info
    character(*),  parameter, public :: MOD_WSM6_DRIVER_AUTHOR = "Song-You Hong and Jeong-Ock Jade Lim, modified by: Bernard Gingold, e-mail: beniekg@gmail.com"
    
    ! Module short description
    character(*),  parameter, public :: MOD_WSM6_DRIVER_DESCRIPT = "Driver for -class GRAUPEL phase microphyiscs scheme (WSM6) "
    




    
    contains


    subroutine wsm6D_driver(  ids,ide,jds,jde,kds,kde      &
                             ,ims,ime,jms,jme,kms,kme      &
                             ,its,ite,jts,jte,kts,kte      &
                             ,t,qci,qrs,q,den,p,delz       &
                             ,rain,rainncv,sr,snow         &
                             ,snowncv,graupel,graupelncv   &
#if defined(USE_PROFILING) && (SAMPLE_PMC) == 1
                             core_counters                                
#endif                                     )

                             

                             
                             

     
#if defined(USE_PROFILING) && (SAMPLE_PMC) == 1
           use ISO_C_BINDING, only : c_size_t
#endif
           use omp_lib
        


          integer(kind=int4), intent(inout) :: ids,ide,jds,jde,kds,kde,  &
                                               ims,ime,jms,jme,kms,kme,  &
                                               its,ite,jts,jte,kts,kte  
          real(kind=sp), allocatable, dimension(:,:,:),   intent(inout)   :: t 
          real(kind=sp), allocatable, dimension(:,:,:,:), intent(inout)   :: qci
          real(kind=sp), allocatable, dimension(:,:,:,:), intent(inout)   :: qrs
          real(kind=sp), allocatable, dimension(:,:,:),   intent(inout)   :: q
          real(kind=sp), allocatable, dimension(:,:,:),   intent(inout)   :: den
          real(kind=sp), allocatable, dimension(:,:,:),   intent(inout)   :: p
          real(kind=sp), allocatable, dimension(:,:,:),   intent(inout)   :: delz
          real(kind=sp), allocatable, dimension(:,:),     intent(inout)   :: rain
          real(kind=sp), allocatable, dimension(:,:),     intent(inout)   :: rainncv
          real(kind=sp), allocatable, dimension(:,:),     intent(inout)   :: sr
          real(kind=sp), allocatable, dimension(:,:),     intent(inout)   :: snow
          real(kind=sp), allocatable, dimension(:,:),     intent(inout)   :: snowncv
          real(kind=sp), allocatable, dimension(:,:),     intent(inout)   :: graupel
          real(kind=sp), allocatable, dimension(:,:),     intent(inout)   :: graupelncv
#if defined(USE_PROFILING) && (SAMPLE_PMC) == 1
          type(core_counters1D),                          intent(inout)   :: core_counters
#endif
         
!       
          real(kind=sp) :: delt,g,rd,rv,t0c,den0,cpd,cpv,ep1,   &
               ep2,qmin,XLS,XLV0,cliq,cice,psat,denr
#if defined(USE_PROFILING) && (SAMPLE_PMC) == 1
          integer(kind=c_size_t), volatile  :: FCRef_start,FCRef_end,  &
                                               FCAct_start,FCAct_end,  &
                                               FCIns_start,FCIns_end
          integer(kind=c_size_t), volatile  :: PMC0_start,PMC0_end,   &
                                               PMC1_start,PMC1_end,   &
                                               PMC2_start,PMC2_end,   &
                                               PMC3_start,PMC3_end,   &
                                               PMC4_start,PMC4_end,   &
                                               PMC5_start,PMC5_end,   &
                                               PMC6_start,PMC6_end,   &
                                               PMC7_start,PMC7_end,   &
                                               TSC_start, TSC_end
          integer(kind=c_size_t), volatile  :: dummy0,dummy1,dummy2,  &
                                               dummy3,dummy4,dummy5,  &
                                               dummy6,dummy7,dummy8,  &
                                               dummy9,dummy10,dummy11,&
                                               dummy12
          integer(kind=c_int),    volatile  :: socket,core
          integer(kind=c_int),    volatile  :: core_ctr_width,fix_ctr_width
#endif
          integer(kind=int4) :: i,j,k,CHUNK,num_tiles_C
          integer(kind=int4) :: ios,unitno
         
          ! So called 'chunk indices'
          integer(kind=int4) :: iids,iide,jjds,jjde, &
                           iims,iime,jjms,jjme, &
                           iits,iite,jjts,jjte
          
          ! Exec code ....
         
          ! Check allocation status, if allocated return immediately
          if( ALLOCATED(t)       .OR.  &
              ALLOCATED(qci)     .OR.  &
              ALLOCATED(qrs)     .OR.  &
              ALLOCATED(q)       .OR.  &
              ALLOCATED(den)     .OR.  &
              ALLOCATED(p)       .OR.  &
              ALLOCATED(delz)    .OR.  &
              ALLOCATED(rain)    .OR.  &
              ALLOCATED(rainncv) .OR.  &
              ALLOCATED(sr)      .OR.  &
              ALLOCATED(snow)    .OR.  &
              ALLOCATED(snowncv) .OR.  &
              ALLOCATED(graupel) .OR.  &
              ALLOCATED(graupelncv)     ) then
                print*, "In File: ", __FILE__," at line: ", __LINE__, " Non-Fatal error: -- allocated array in wsm6D_driver!!"
                err = -1
                return
          end if
          fn = "wsm6_constants.dat"
          unitno = 31
          open(unitno,file = trim(fn),form = "unformatted",action = "read", &
               iostat=ios )
          if(ios /= 0) then
              print*, "In File: ", __FILE__, " at line: ",__LINE__, & 
                      " FATAL-ERROR: Failed to open file: ",trim(fn)
              !err = ios
              ERROR STOP "FATAL-ERROR: Failed to open file: ",trim(fn)
          end if
          read(unitno) pi,xlv1
          read(unitno) qc0,qck1
          read(unitno) bvtr1, bvtr2, bvtr3, bvtr4, bvtr6, g1pbr, g3pbr, &
                       g4pbr, g6pbr, g5pbro2, pvtr, eacrr, pacrr, &
                       precr1, precr2, roqimax
          read(unitno) bvts1, bvts2, bvts3, bvts4, g1pbs, g3pbs, g4pbs,  &
                       g5pbso2, pvts, pacrs, precs1, precs2, pidn0r, pidn0s
          read(unitno) pacrc
          read(unitno) bvtg1, bvtg2, bvtg3, bvtg4, g1pbg, g3pbg, g4pbg,  &
                       pacrg, g5pbgo2, pvtg, precg1, precg2, pidn0g
          read(unitno) rslopermax, rslopesmax, rslopegmax, rsloperbmax,  &
                       rslopesbmax, rslopegbmax, rsloper2max, rslopes2max,  &
                       rslopeg2max, rsloper3max, rslopes3max, rslopeg3max
          close(unitno) 
          ! Read input data
          fn = "wsm6_input.dat"
          unitno = 31
          open(unitno,file=trim(fn),form="unformatted",action="read", &
               iostat=ios )
          if(ios /= 0) then
              print*, "In File: ",__FILE__, "at line: ",__LINE__, &
                      "FATAL-ERROR: Failed to open file: ",trim(fn)
              ERROR STOP "FATAL-ERROR: Failed to open file: ",trim(fn)
          end if
          read(unitno) ids,ide,jds,jde,kds,kde, &
                       ims,ime,jms,jme,kms,kme, &
                       its,ite,jts,jte,kts,kte
          ! Check indices
          if((ims/=its) .OR.   &
             (ime/=ite) .OR.   &
             (jms/=jts) .OR.   &
             (jme/=jte)         ) then
                print*, " In File: ",__FILE__, "at line: ",__LINE__, &
                        " FATAL-ERROR: Index mismatch found in file: ", trim(fn)
                ERROR STOP " FATAL-ERROR: Index mismatch found in file: ", trim(fn)
          end if
          if((ims/=1) .OR. (jms/=1)) then
               print*, " In File: ",__FILE__, "at line: ",__LINE__, &
                        " FATAL-ERROR: Incorrect start index found in file: ", trim(fn)
               ERROR STOP " FATAL-ERROR: Incorrect start index found in file: ", trim(fn)
          end if
           ! set default values of "chunk" indices
          iids = ids
          iide = ide
          iims = ims
          iime = ime
          iits = its
          iite = ite
          jjds = jds
          jjde = jde
          jjms = jms
          jjme = jme
          jjts = jts
          jjte = jte
          CHUNK = iite-iits+1
          num_tiles_C = (iite-iits+1) / CHUNK
          if(mod((iite-iits+1),CHUNK) > 0) then
              num_tiles_C = num_tiles_C + 1
          end if
          iime = CHUNK
          iite = CHUNK
          jjme = num_tiles_C
          jjte = num_tiles_C
          !omp_set_num_threads(nthreads)
          ! Array allocation
          allocate(t(iits:iite,kts:kte,jjts:jjte),     &
                   qci(iits:iite,kts:kte,2,jjts:jjte), &
                   qrs(iits:iite,kts:kte,3,jjts:jjte), &
                   q(iims:iime,kms:kme,jjms:jjme),     &
                   den(iims:iime,kms:kme,jjms:jjme),   &
                   p(iims:iime,kms:kme,jjms:jjme),     &
                   delz(iims:iime,kms:kme,jjms:jjme),  &
                   rain(iims:iime,jjms:jjme),          &
                   rainncv(iims:iime,jjms:jjme),       &
                   sr(iims:iime,jjms:jjme),            &
                   snow(iims:iime,jjms:jjme),          &
                   snowncv(iims:iime,jjms:jjme),       &
                   graupel(iims:iime,jjms:jjme),       &
                   graupelncv(iims:iime,jjms:jjme))   
#if defined(USE_PROFILING) && (SAMPLE_PMC) == 1
          call init_core_counters1D(core_counters,jjts,jjte)
#endif
          
            

          
!$OMP PARALLEL DO
!$OMP PRIVATE(j)
!$OMP SCHEDULE(runtime)
          do j = jjts, jjte
              call firstTouch(t(iits,kts,j),     q(iims,kms,j),     &
                              qci(iits,kts,1,j), qrs(iits,kts,1,j), &
                              den(iims,kms,j),                      &
                              p(iims,kms,j),     delz(iims,kms,j),  &
                              j,                                    &
                              rain(iims,j),      rainncv(iims,j),   &
                              sr(iims,j),                           &
                              iids,iide, jjds,jjde, kds,kde,        &
                              iims,iime, jjms,jjme, kms,kme,        &
                              iits,iite, jjts,jjte, kts,kte,        &
                              snow,snowncv,                         &
                              graupel,graupelncv                    )

          end do
!$OMP END PARALLEL DO
           ! read remaining input data
           call readarray3(t,'t',unitno,its,ite)
           call readarray4(qci,'qci',unitno,its,ite)
           call readarray4(qrs,'qrs',unitno,its,ite)
           call readarray3(q,'q',unitno,its,ite)
           call readarray3(den,'den',unitno,its,ite)
           call readarray3(p,'p',unitno,its,ite)
           call readarray3(delz,'delz',its,ite)
           read(unitno) delt,g,cpd,cpv,t0c,den0,  &
                        rd,rv,ep1,ep2,qmin,XLS,   &
                        XLV0,XLF0,cliq,cice,psat, &
                        denr
           call readarray2(rain,'rain',unitno,its,ite)
           call readarray2(rainncv,'rainncv',unitno,its,ite)
           call readarray2(sr,'sr',unitno,its,ite)
           call readarray2(snow,'snow',unitno,its,ite)
           call readarray2(snowncv,'snowncv',unitno,its,ite)
           call readarray2(graupel,'graupel',unitno,its,ite)
           call readarray2(graupelncv,'graupelncv',unitno,its,ite)
           close(unitno)
      
#if defined(USE_PROFILING) && (SAMPLE_PMC) == 1
           ! Kernel warmp-up run
!$OMP PARALLEL DO
!$OMP PRIVATE(j)
!$OMP SCHEDULE(runtime)
           do j = jjts, jjte
             
               call wsm62D( t(iits,kts,j),     q(iims,kms,j),      &
                            qci(iits,kts,1,j), qrs(iits,kts,1,j),  &
                            den(iims,kms,j),                       &
                            p(iims,kms,j),     delz(iims,kms,j),   &
                            delt,g,cpd,cpv,rd,rv,t0c,              &
                            ep1,ep2,qmin,                          &
                            XLS,XLV0,XLF0,den0,denr,               &
                            cliq,cice,psat,                        &
                            j,                                     &
                            rain(iims,j),      rainncv(iims,j),    &
                            sr(iims,j),                            &
                            iids,iide, jjds,jjde, kds,kde,         &
                            iims,iims, jjms,jjme, kms,kme,         &
                            iits,iite, jjts,jjte, kts,kte,         &
                            snow,snowncv,                          &
                            graupel,graupelncv                    )
           
           end do
!$OMP END PARALLEL DO           
#endif

           !Now the profiled run
#if defined(USE_PROFILING) && (SAMPLE_PMC) == 1
           !Minimize the L1I and ITLB noise pollution
           !Counter measurement functions probably
           !will not be inlined.
           dummy0 = rdpmc_instructions()
           dummy1 = rdpmc_actual_cycles()
           dummy2 = rdpmc_reference_cycles()
           dummy4 = full_rdtscp(socket,core)
           dummy5 = rdpmc(0)
           dummy6 = rdpmc(1)
           dummy7 = rdpmc(2)
           dummy8 = rdpmc(3)
           dummy9 = rdpmc(4)
           dummy10= rdpmc(5)
           dummy11= rdpmc(6)
           dummy12= rdpmc(7)
           core_ctr_width = get_core_counter_width()
           fix_ctr_width  = get_fixed_counter_width()
           dummy13= corrected_pmc_delta(dummy1,dummy2,core_ctr_width)
           if(core_ctr_width /= fix_ctr_width) then
              print*, "Warning -- programmable counters width does not match the fixed counters width!!"
           endif
!$OMP PARALLEL DO
!$OMP PRIVATE(j,TSC_start,TSC_end,FCRef_start,FCRef_end,FCAct_start,FCAct_end, &
!$OMP&        FCIns_start,FCIns_end,PMC0_start,PMC0_end,PMC1_start,PMC1_end,   &
!$OMP&        PMC2_start,PMC2_end,PMC3_start,PMC3_end,PMC4_start,PMC4_end,     &
!$OMP&        PMC5_start,PMC5_end,PMC6_start,PMC6_end,PMC7_start,PMC7_end)
!$OMP SHARED(core_counters)           
!$OMP SCHEDULE(runtime)
           do j = jjts, jjte
               call FOR_LFENCE()
               TSC_start   = rdtscp()
               FCRef_start = rdpmc_reference_cycles()
               FCAct_start = rdpmc_actual_cycles()
               FCIns_start = rdpmc_instructions()
               PMC0_start  = rdpmc(0)
               PMC1_start  = rdpmc(1)
               PMC2_start  = rdpmc(2)
               PMC3_start  = rdpmc(3)
               PMC4_start  = rdpmc(4)
               PMC5_start  = rdpmc(5)
               PMC6_start  = rdpmc(6)
               PMC7_start  = rdpmc(7)
               call FOR_LFENCE()
               call wsm62D( t(iits,kts,j),     q(iims,kms,j),      &
                            qci(iits,kts,1,j), qrs(iits,kts,1,j),  &
                            den(iims,kms,j),                       &
                            p(iims,kms,j),     delz(iims,kms,j),   &
                            delt,g,cpd,cpv,rd,rv,t0c,              &
                            ep1,ep2,qmin,                          &
                            XLS,XLV0,XLF0,den0,denr,               &
                            cliq,cice,psat,                        &
                            j,                                     &
                            rain(iims,j),      rainncv(iims,j),    &
                            sr(iims,j),                            &
                            iids,iide, jjds,jjde, kds,kde,         &
                            iims,iims, jjms,jjme, kms,kme,         &
                            iits,iite, jjts,jjte, kts,kte,         &
                            snow,snowncv,                          &
                            graupel,graupelncv                    )
               call FOR_LFENCE()
               TSC_end   = rdtscp()
               FCRef_end = rdpmc_reference_cycles()
               FCAct_end = rdpmc_actual_cycles()
               FCIns_end = rdpmc_instructions()
               PMC0_end  = rdpmc(0)
               PMC1_end  = rdpmc(1)
               PMC2_end  = rdpmc(2)
               PMC3_end  = rdpmc(3)
               PMC4_end  = rdpmc(4)
               PMC5_end  = rdpmc(5)
               PMC6_end  = rdpmc(6)
               PMC7_end  = rdpmc(7)
               call FOR_LFENCE()
               core_counters.tsc(j) = corrected_pmc_delta(TSC_end,TSC_start,core_ctr_width)
               core_counters.r30b(j)= corrected_pmc_delta(FCRef_end,FCRef_start,core_ctr_width)
               core_counters.r30a(j)= corrected_pmc_delta(FCAct_end,FCAct_start,core_counter_width)
               core_counters.c0(j)  = corrected_pmc_delta(PMC0_end,PMC0_start,core_counter_width)
               core_counters.c1(j)  = corrected_pmc_delta(PMC1_end,PMC1_start,core_counter_width)
               core_counters.c2(j)  = corrected_pmc_delta(PMC2_end,PMC2_start,core_counter_width)
               core_counters.c3(j)  = corrected_pmc_delta(PMC3_end,PMC3_start,core_counter_width)
               core_counters.c4(j)  = corrected_pmc_delta(PMC4_end,PMC4_start,core_counter_width)
               core_counters.c5(j)  = corrected_pmc_delta(PMC5_end,PMC5_start,core_counter_width)
               core_counters.c6(j)  = corrected_pmc_delta(PMC6_end,PMC6_start,core_counter_width)
               core_counters.c7(j)  = corrected_pmc_delta(PMC7_end,PMC7_start,core_counter_width)
           end do
!$OMP END PARALLEL DO           
#else
!$OMP PARALLEL DO
!$OMP PRIVATE(j)
!$OMP SCHEDULE(runtime)
           do j = jjts, jjte
             
               call wsm62D( t(iits,kts,j),     q(iims,kms,j),      &
                            qci(iits,kts,1,j), qrs(iits,kts,1,j),  &
                            den(iims,kms,j),                       &
                            p(iims,kms,j),     delz(iims,kms,j),   &
                            delt,g,cpd,cpv,rd,rv,t0c,              &
                            ep1,ep2,qmin,                          &
                            XLS,XLV0,XLF0,den0,denr,               &
                            cliq,cice,psat,                        &
                            j,                                     &
                            rain(iims,j),      rainncv(iims,j),    &
                            sr(iims,j),                            &
                            iids,iide, jjds,jjde, kds,kde,         &
                            iims,iims, jjms,jjme, kms,kme,         &
                            iits,iite, jjts,jjte, kts,kte,         &
                            snow,snowncv,                          &
                            graupel,graupelncv                    )
           
           end do
!$OMP END PARALLEL DO
#endif        
    end subroutine
           
  

end module wsm6_driver
