
#include "Config.fpp"

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
    
    !======================================
    ! type: WSM6Type_t
    !======================================
    type, public :: WSM6Type_t
        
          public
           ! WRF indexing
          integer(kind=i4) :: ids,ide,jds,jde,kds,kde,  &
                                ims,ime,jms,jme,kms,kme,  &
                                its,ite,jts,jte,kts,kte
          ! Meteo fields
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:,:)   :: t
          !DIR$ ATTRIBUTES ALIGN : 64 :: t
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:,:)   :: t    !GCC$ ATTRIBUTES aligned(64) :: t
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:,:,:) :: qci
          !DIR$ ATTRIBUTES ALIGN : 64 :: qci
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:,:,:) :: qci  !GCC$ ATTRIBUTES aligned(64) :: qci
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:,:,:) :: qrs
          !DIR$ ATTRIBUTES ALIGN : 64 :: qrs
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:,:,:) :: qrs  !GCC$ ATTRIBUTES aligned(64) :: qrs
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:,:)   :: q
          !DIR$ ATTRIBUTES ALIGN : 64 :: q
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:,:)   :: q    !GCC$ ATTRIBUTES aligned(64) :: q
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:,:)   :: den
          !DIR$ ATTRIBUTES ALIGN : 64 :: den
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:,:)   :: den  !GCC$ ATTRIBUTES aligned(64) :: den
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:,:)   :: p
          !DIR$ ATTRIBUTES ALIGN : 64 :: p
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:,:)   :: p    !GCC$ ATTRIBUTES aligned(64) :: p
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:,:)   :: delz
          !DIR$ ATTRIBUTES ALIGN : 64 :: delz
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:,:)   :: delz !GCC$ ATTRIBUTES aligned(64) :: delz
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:)     :: rain
          !DIR$ ATTRIBUTES ALIGN : 64 :: rain
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:)     :: rain !GCC$ ATTRIBUTES aligned(64) :: rain
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:)     :: rainncv
          !DIR$ ATTRIBUTES ALIGN : 64 :: rainncv
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:)     :: rainncv !GCC$ ATTRIBUTES aligned(64) :: rainncv
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:)     :: sr
          !DIR$ ATTRIBUTES ALIGN : 64 :: sr
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:)     :: sr     !GCC$ ATTRIBUTES aligned(64) :: sr
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:)     :: snow
          !DIR$ ATTRIBUTES ALIGN : 64 :: snow
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:)     :: snow   !GCC$ ATTRIBUTES aligned(64) :: snow
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:)     :: snowncv
          !DIR$ ATTRIBUTES ALIGN : 64 :: snowncv
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:)     :: snowncv !GCC$ ATTRIBUTES aligned(64) :: snowcv
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:)     :: graupel
          !DIR$ ATTRIBUTES ALIGN : 64 :: graupel
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:)     :: graupel !GCC$ ATTRIBUTES aligned(64) :: graupel
#endif
#if defined __INTEL_COMPILER
          real(kind=sp), allocatable, dimension(:,:)     :: graupelnvc
          !DIR$ ATTRIBUTES ALIGN : 64 :: graupelnvc
#elif defined __GFORTRAN__
          real(kind=sp), allocatable, dimension(:,:)     :: graupelnvc !GCC$ ATTRIBUTES aligned(64) :: graupelnvc
#endif
    end type WSM6Type_t
    
    contains

#if defined __GFORTRAN__
    subroutine wsm6D_driver(  ids,ide,jds,jde,kds,kde,     &
                             ims,ime,jms,jme,kms,kme,      &
                             its,ite,jts,jte,kts,kte,      &
                             t,qci,qrs,q,den,p,delz,       &
                             rain,rainncv,sr,snow,         &
                             snowncv,graupel,graupelncv,   &
                             logging,verbose,fname,append, &
                             nthreads,err                   )    !GCC$ ATTRIBUTES cold :: wsm6D_driver  !GCC$ ATTRIBUTES aligned(32) :: wsm6D_driver
#elif defined __INTEL_COMPILER
    subroutine wsm6D_driver( ids,ide,jds,jde,kds,kde,      &
                             ims,ime,jms,jme,kms,kme,      &
                             its,ite,jts,jte,kts,kte,      &
                             t,qci,qrs,q,den,p,delz,       &
                             rain,rainncv,sr,snow,         &
                             snowncv,graupel,graupelncv,   &
                             logging,verbose,fname,append, &
                             nthreads,err                )

      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: wsm6D_driver
#endif
          use omp_lib
          use mod_print_error, only : handle_fatal_memory_error
!          use mod_code_timing

          integer(kind=i4), intent(inout) :: ids,ide,jds,jde,kds,kde,  &
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
         
          logical(kind=i4),                               intent(in)      :: logging,verbose
          character(len=*),                                 intent(in)      :: fname
          logical(kind=i4),                               intent(in)      :: append
          integer(kind=i4),                               intent(in)      :: nthreads
          integer(kind=i4),                               intent(inout)   :: err
          
!          type(QPCTimer_t),allocatable, dimension(:),  intent(inout)   :: timers -- dummy argument
          !Locals
          real(kind=sp) :: delt,g,rd,rv,t0c,den0,cpd,cpv,ep1,   &
                        ep2,qmin,XLS,XLV0,cliq,cice,psat,denr
          integer(kind=i4) :: i,j,k,CHUNK,num_tiles_C
          integer(kind=i4) :: ios,unitno
          character(len=256) :: emsg
          character(len=64)  :: fn ! filename
          integer(kind=i4)      :: aerr
          ! So called 'chunk indices'
          integer(kind=i4) :: iids,iide,jjds,jjde, &
                           iims,iime,jjms,jjme, &
                           iits,iite,jjts,jjte
          integer(kind=1)      :: ifail
          logical(kind=i4)      :: bfail
          ! Exec code ....
          if(err < 0) err = 0
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
          omp_set_num_threads(nthreads)
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
                   graupelncv(iims:iime,jjms:jjme),    &
                   timers(jjts:jjte),                  &
                   STAT=aerr,                          &
                   ERRMSG=emsg                  )
          if(aerr /= 0) then
              call handle_fatal_memory_error(logging,verbose,append,fname,   &
                                             "logger:320 --> mod_wsm6_driver/wsm6D_driver: -- Memory Allocation Failure !!!", &
                                             " mod_wsm6_driver/wsm6D_driver:320 -- Memory Allocation Failure !!!" ,  &
                                             emsg,__LINE__  )
          end if
        
            
         
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
          ! do j = jjts, jjte
          !     call qpctimer_init(timers(j),"WSM6D thread timer: #"//trim(to_str(j)), &
          !                        "module_wsm6_driver/wsm6D_driver"              )
          ! end do
           
!$OMP PARALLEL DO
!$OMP PRIVATE(j)
!$OMP SCHEDULE(runtime)
           do j = jjts, jjte
              ! call qpctimer_start(timers(j),ifail)
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
              ! call qpctimer_stop(timers(j),ifail)
           end do
!$OMP END PARALLEL DO
          ! print*, " Crude results of QueryPerformanceCounter measurements."
          ! do j = jjts, jjte
          !     call qpctimer_delta(timers(j),ifail)
          !     call qpctimer_print(timers(j))
          ! end do
    end subroutine
           
    character(len=20) function to_str(ivalue)
          implicit none
          integer(kind=i4), intent(in) :: ivalue
          ! Exec code ...
          write(to_str,*) ivalue
          to_str = adjustl(to_str)
    end function

end module wsm6_driver
