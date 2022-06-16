
#include "Config.fpp"

module mod_wind_wrf


 !===================================================================================85
 !---------------------------- DESCRIPTION ------------------------------------------85
 !
 !
 !
 !          Module  name:
 !                         'mod_wind_wrf'
 !          
 !          Purpose:
 !                       This module is responsible for initialization
 !                       3D Wind fields, i.e. U,V,W vectors.
 !                       The data is computed by the WRF numerical simulation
 !                       This module is based on CR-SIM v3.0 implementation
 !                       on WRF NetCDF output reading subroutines.
 !                       
 !          History:
 !                        Date: 17-06-2018
 !                        Time: 11:19 GMT+2
  !                        Modified
  !                       Date: 16-12-2019
  !                       Time: 21:04 GMT+2
 !          Version:
 !
 !                      Major: 1
 !                      Minor: 1
 !                      Micro: 0
 !
 !          Author:  
 !                      Copyright (C) 2014 McGill Clouds Research Group (//radarscience.weebly.com//)
!!                      Contact email address: aleksandra.tatarevic@mcgill.ca pavlos.kollias@mcgill.ca
 !          Modified:
 !                   Bernard Gingold on 17-06-2018
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

   
    use mod_kinds, only : i4,i8,sp,dp
    use netcdf
    implicit none
    !=====================================================59
    !  File and module information:
    !  version,creation and build date, author,description
    !=====================================================59
    
    ! Major version
    integer(kind=i4), parameter, public :: MOD_WIND_WRF_MAJOR = 1
    
    ! Minor version
    integer(kind=i4), parameter, public :: MOD_WIND_WRF_MINOR = 1
    
    ! Micro version
    integer(kind=i4), parameter, public :: MOD_WIND_WRF_MICRO = 0
    
    ! Module full version
    integer(kind=i4), parameter, public :: MOD_WIND_WRF_FULLVER = 1000*MOD_WIND_WRF_MAJOR+  &
                                                               100*MOD_WIND_WRF_MINOR+   &
                                                               10*MOD_WIND_WRF_MICRO
    
    ! Module creation date
    character(*),  parameter, public :: MOD_WIND_WRF_CREATE_DATE = "17-06-2018 11:26 +00200 (SUN 17 JUN 2018 GMT+2)"
    
    
    character(*),  parameter, public :: MOD_WIND_WRF_BUILD_DATE =  __DATE__ " " __TIME__

    ! Module author
    character(*),  parameter, public :: MOD_WIND_WRF_AUTHOR = "Based on CR-SIM v3.0, modified by Bernard Gingold, beniekg@gmail.com"
    
    ! Module short description
    character(*),  parameter, public :: MOD_WIND_WRF_DESCRIPT = "Wind 3D fields."
    
    !================================
    !  type: Wind3D_t
    !================================
    type, public :: Wind3D_t
        
          public
          
          integer(kind=i4) :: nt   ! time evolution
          
          integer(kind=i4) :: nx,nxp1  ! dimension x, staggered nxp1
          
          integer(kind=i4) :: ny,nyp1  ! dimension y, staggered nyp1
          
          integer(kind=i4) :: nz,nzp1  ! dimension z, staggered nzp1
      
          real(kind=dp), allocatable, dimension(:,:,:,:) :: U   ! U-component
!DIR$     ATTRIBUTES ALIGN : 64 :: U


          real(kind=dp), allocatable, dimension(:,:,:,:) :: V   ! V-component
!DIR$     ATTRIBUTES ALIGN : 64 :: V

#if defined __INTEL_COMPILE
          real(kind=dp), allocatable, dimension(:,:,:,:) :: W   ! W-component
!DIR$     ATTRIBUTES ALIGN : 64 :: W        

    end type Wind3D_t
    
    contains


     
      subroutine init_wind3D(wind3d,fwrfout,statd,stata,logging,verbose,append,fname)

        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: init_wind3D
        use IFPORT,          only : TRACEBACKQQ

          use mod_print_error, only : handle_fatal_memory_error,  &
                                      print_fatal_error
          type(Wind3D_t),   intent(inout) :: wind3d
          character(len=*), intent(in)    :: fwrfout
          integer(kind=i4),    intent(out)   :: statd,stata
          logical(kind=i4),    intent(in)    :: logging,verbose,append
          character(len=*), intent(in)    :: fname
          ! Locals
          character(len=nf90_max_name) :: name
          character(len=256) :: emsg,nf90err
          integer(kind=i4) ::  nDimsd,nVarsd,length,iDim,ncidd,  &  ! dimension variables
                            nDimsv,nVarsv,ncidv,iVarv,aerr
          integer(kind=i4), allocatable, dimension(:)                :: dim_lengths
          character(len=nf90_max_name), allocatable, dimension(:) :: dim_names
          ! Exec code .....
          wind3D%nt   = 0
          wind3D%nx   = 0
          wind3D%nxp1 = 0
          wind3D%ny   = 0
          wind3D%nyp1 = 0
          wind3D%nz   = 0
          wind3D%nzp1 = 0
          statd = nf90_open(trim(fwrfout), NF90_NOWRITE,ncidd)
          if(statd /= 0) then 
              call print_fatal_error("!!!Fatal-Error!!! -- in module: mod_wind_wrf in subroutine:  int_wind3D",    &
                                     " nf90_open failed to open file: fwrfout",    &
                                     "No system message available!!",              &
                                     __FILE__, __LINE__  )
#if defined __INTEL_COMPILER
              call TRACEBACKQQ(STRING="nf90_open failed",USER_EXIT_CODE = -1)
#elif defined __GFORTRAN__
              call backtrace()
#endif
              ERROR STOP
          end if
          statd = nf90_inquire(ncidd,nDimsd,nVarsd)
          if(statd /= 0) then
              call print_fatal_error("!!!Fatal-Error!!! -- in module: mod_wind_wrf in subroutine: init_wind3D", &
                                     " nf90_inquire failed!!",  &
                                     "No system message available!!",  &
                                     __FILE__,__LINE__ )
#if defined __INTEL_COMPILER
              call TRACEBACKQQ(STRING="nf90_inquire failed",USER_EXIT_CODE = -1)
#elif defined __GFORTRAN__
              call backtrace()
#endif
              ERROR STOP
          end if
          allocate(dim_names(1:nDimsd),  &
                   dim_lengths(1:nDimsd))
          
          do iDim = 1, nDimsd
              
              statd = nf90_inquire_dimension(ncidd,iDim,name=name,len=length)
              dim_names(iDim) = name
              dim_lengths(iDim) = length
              select case (name)
                 case ('Time')             ;   wind3D%nt   = length
                 case ('west_east')        ;   wind3D%nx   = length
                 case ('west_east_stag')   ;   wind3D%nxp1 = length
                 case ('south_north')      ;   wind3D%ny   = length
                 case ('south_north_stag') ;   wind3D%nyp1 = length
                 case ('bottom_top')       ;   wind3D%nz   = length
                 case ('bottom_top_stag')  ;   wind3D%nzp1 = length
             end select
          end do
          statd = nf90_close(ncidd)
          if(wind3D%nx+1 /= wind3D%nxp1) then
              print*, " Invalid indices: ", wind3D%nx+1 /= wind3Dnxp1
              ERROR STOP "!!! Fatal-Error !!! in init_wind3D ... exiting!!"
          end if
          if(wind3D%ny+1 /= wind3D%nyp1) then
              print*, " Invalid indices: ", wind3D%ny+1 /= wind3D%nyp1
              ERROR STOP "!!! Fatal-Error !!! in init_wind3D ... exiting!!"
          end if
         if(allocated(wind3D%U)) then
             deallocate(wind3D%U)
         end if
         if(allocated(wind3D%V)) then
             deallocate(wind3D%V)
         end if
         if(allocated(wind3D%W)) then
             deallocate(wind3D%W)
         end if
         associate(nt=>wind3D%nt,nx=>wind3D%nx,         &
                   ny=>wind3D%ny,nz=>wind3D%nz,         &
                   nxp1=>wind3D%nxp1,nyp1=>wind3D%nyp1, &
                   nzp1=>wind3D%nzp1      )
             allocate(wind3D%U(nxp1,ny,nz,nt),      &
                      wind3D%V(nx,nyp1,nz,nt),      &
                      wind3D%W(nx,ny,nzp1,nt),      &
                      STAT=aerr,                    &
                      ERRMSG=emsg          )
         end associate
         if(aerr /= 0) then
             call handle_fatal_memory_error(logging,verbose,append,fname,   &
                                             "logger:201 --> mod_wind_wrf/init_wind3D: Memory Allocation Failure!!" , &
                                             "mod_wind_wrf/init_wind3D:201 -- Memory Allocation Failure!!" , &
                                             emsg,__LINE__)
         end if
         stata = nf90_open(trim(fwrfout), NF90_NOWRITE,ncidv)
         if(stata /= 0) then
              call print_fatal_error("!!!Fatal-Error!!! -- in module: mod_wind_wrf in subroutine:  int_wind3D",    &
                                     "nf90_open failed to open file: fwrfout",    &
                                     "No system message available!!",              &
                                     __FILE__, __LINE__  )
#if defined __INTEL_COMPILER
              call TRACEBACKQQ(STRING="nf90_open failed",USER_EXIT_CODE = -1)
#elif defined __GFORTRAN__
              call backtrace()
#endif
              ERROR STOP
         end if
         stata = nf90_inquire(ncidv,nDimsv,nVarsv)
         if(stata /= 0) then
              call print_fatal_error("!!!Fatal-Error!!! -- in module: mod_wind_wrf in subroutine:  int_wind3D",    &
                                     "nf90_inquire failed !!",    &
                                     "No system message available!!",              &
                                     __FILE__, __LINE__  )
#if defined __INTEL_COMPILER
              call TRACEBACKQQ(STRING="nf90_inquire failed",USER_EXIT_CODE = -1)
#elif defined __GFORTRAN_
              call backtrace()
#endif
              ERROR STOP
         end if
         
         do iVarv = 1, nVarsv
             
             stata = nf90_inquire_variable(ncidv,iVarv,name=name)
             if(stata /= 0) then
                 nf90err = "nf90_inquire_variable -- failed!!"
                 goto 9999
             end if
             select case (trim(adjust(name)))
             case ('U')
                 stata = nf90_get_var(ncidv,iVarv,wind3D%U)
                 if(stata /= 0) then
                     nf90err = "nf90_get_var -- failed to copy U-component!!"
                     goto 9999
                 end if
             case ('V')
                 stata = nf90_get_var(ncidv,iVarv,wind3D%V)
                 if(stata /= 0) then
                     nf90err = "nf90_get_var -- failed to copy V-component!!"
                     goto 9999
                 end if
             case ('W')
                 stata = nf90_get_var(ncidv,iVarv,wind3D%W)
                  if(stata /= 0) then
                      nf90err = "nf90_get_var -- failed to copy W-component!!"
                      goto 9999
                  end if
             end select
             
         end do
         
9999     if(stata /= 0) then
            write(6,*) nf90err
#if defined __INTEL_COMPILER
            call TRACEBACKQQ(STRING=nf90err,USER_EXIT_CODE = -1)
#elif defined  __GFORTRAN__
            call backtrace()
#endif
             ERROR STOP nf90err
         else
             stata = nf90_close(ncidv)
         end if
         
    end subroutine


end module mod_wind_wrf
