module GMS_mod_spectrum

    !
    ! wavy - A spectral ocean wave modeling and development framework
    ! Copyright (c) 2017, Wavebit Scientific LLC
    ! All rights reserved.
    !
    ! Licensed under the BSD-3 clause license. See LICENSE for details.
    ! Modified by Bernard Gingold (beniekg@gmail.com) on 16/03/2019
    !

    use ,intrinsic :: ISO_C_BINDING, Only :  c_int, c_double
    use GMS_mod_kinds,     only : int32_t, dp
    use GMS_mod_utility,   only : diff, diff_periodic
    use GMS_mod_datetime,  only : datetime
    use GMS_mod_timedelta, only : timedelta
    implicit none

    private

    type, public :: spectrum_type

          public

          type(datetime)  :: start_time
          type(datetime)  :: end_time
          type(timedelta) :: time_step
!DIR$   ATTRIBUTES ALIGN : 64 :: spec
          real(kind=dp), allocatable, dimension(:,:) :: spec ! 2-d spectrum
!DIR$   ATTRIBUTES ALIGN : 64 :: f
          real(kind=dp), allocatable, dimension(:)   :: f    ! Frequency [Hz]
!DIR$   ATTRIBUTES ALIGN : 64 :: df
          real(kind=dp), allocatable, dimension(:)   :: df   ! Frequency spacing [Hz]
!DIR$   ATTRIBUTES ALIGN : 64 :: k
          real(kind=dp), allocatable, dimension(:)   :: k    ! Wavenumber [rad/m]
!DIR$   ATTRIBUTES ALIGN : 64 :: dk
          real(kind=dp), allocatable, dimension(:)   :: dk   ! Wavenumber spacing [rad/m]
!DIR$   ATTRIBUTES ALIGN : 64 ::  th
          real(kind=dp), allocatable, dimension(:)   :: th   ! Direction [rad]
!DIR$   ATTRIBUTES ALIGN : 64 :: dth
          real(kind=dp), allocatable, dimension(:)   :: dth  ! Directional spacing [rad]
!DIR$   ATTRIBUTES ALIGN : 64 :: cp
          real(kind=dp), allocatable, dimension(:)   :: cp   ! Phase speed [m/s]
!DIR$   ATTRIBUTES ALIGN : 64 :: cg
          real(kind=dp), allocatable, dimension(:)   :: cg   ! Group speed [m/s]
!DIR$   ATTRIBUTES ALIGN : 64 :: u
          real(kind=dp), allocatable, dimension(:)   :: u    ! Mean current velocity in x-direction [m/s]
!DIR$   ATTRIBUTES ALIGN : 64 :: v
          real(kind=dp), allocatable, dimension(:)   :: v    !  Mean current velocity in y-direction [m/s]
!DIR$   ATTRIBUTES ALIGN : 64 :: z
          real(kind=dp), allocatable, dimension(:)   :: z    ! Depth levels for current array [m]
          real(kind=dp) :: air_density     !! Air density [kg/m^3]
          real(kind=dp) :: depth           !! Mean water depth [m]
          real(kind=dp) :: elevation       !! Mean surface elevation [m]
          real(kind=dp) :: grav            !! Gravitational acceleration [m/s^2]
          real(kind=dp) :: surface_tension !! Surface tension [N/m]
          real(kind=dp) :: water_density   !! Water density [kg/m^3]

    end type spectrum_type

    contains

    subroutine constructor(spectrum,fmin,fmax,df,ndirs,depth,grav,    &
                 air_density,water_density,surface_tension,           &
                 errstate,iounit,logging,verbose,append,fname)
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: constructor
          use GMS_mod_print_error,  only : handle_fatal_memory_error,  &
                                       print_non_fatal_error
          type(spectrum_type),      intent(inout) :: spectrum
          real(kind=dp),            intent(in)    :: fmin
          real(kind=dp),            intent(in)    :: fmax
          real(kind=dp),            intent(in)    :: df
          integer(kind=int32_t),    intent(in)    :: ndirs
          real(kind=dp),            intent(in)    :: depth
          real(kind=dp),            intent(in), optional :: grav
          real(kind=dp),            intent(in), optional :: air_density
          real(kind=dp),            intent(in), optional :: water_density
          real(kind=dp),            intent(in), optional :: surface_tension
          logical(kind=int32_t),    intent(inout) :: errstate
          integer(kind=int32_t),    intent(in)    :: iounit
          logical(kind=int32_t),    intent(in)    :: logging
          logical(kind=int32_t),    intent(in)    :: verbose
          logical(kind=int32_t),    intent(in)    :: append
          character(len=*),         intent(in)    :: fname
          ! Locals
!DIR$   ATTRIBUTES ALIGN : 64 :: tmp_f,tmp_k,th_tmp
          real(kind=dp), allocatable, dimension(:) :: tmp_f,tmp_k,th_tmp
          character(len=256),    automatic :: emsg

          integer(kind=int32_t), automatic :: aerr
          integer(kind=int32_t), automatic :: n
          integer(kind=int32_t), automatic :: nfreqs
          ! EXec code.....
          if(present(grav)) then
                spectrum.grav = grav
          else
                spectrum.grav = 9.8_dp
          end if
          if(present(air_density)) then
                spectrum.air_density = air_density
          else
                spectrum.air_density = 1.2_dp
          end if
          if(present(water_density)) then
                spectrum.water_density = water_density
          else
                spectrum.water_density = 1.0e3_dp
          end if
          if(present(surface_tension)) then
                spectrum.surface_tension = surface_tension
          else
                spectrum.surface_tension = 0.07_dp
          end if
          spectrum.depth = depth
          if(fmin == fmax) then
                nfreqs = 1
          else
                nfreqs = int((log(fmax)-log(fmin))/log(df))
          end if
          if(allocated(spectrum.spec)) then
                deallocate(spectrum.spec)
                allocate(spectrum.spec(nfreqs,ndirs),  &
                            STAT=aerr,ERRMSG=emsg)
                if(aerr /= 0) goto 9999
          else
                allocate(spectrum.spec(nfreqs,ndirs),  &
                            STAT=aerr,ERRMSG=emsg)
                if(aerr /= 0) goto 9999
          end if
          spectrum.spec = 0.0_dp
          if(allocated(spectrum.f)) then
                 deallocate(spectrum.f)
                 allocate(spectrum.f(nfreqs),  &
                             STAT=aerr,ERRMSG=emsg)
                  if(aerr /= 0) goto 9999
          else
                  allocate(spectrum.f(nfreqs),  &
                             STAT=aerr,ERRMSG=emsg)
                  if(aerr /= 0) goto 9999
           end if
           if(allocated(spectrum.df)) then
                   deallocate(spectrum.df)
                   allocate(spectrum.df(nfreqs), &
                               STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           else
                   allocate(spectrum.df(nfreqs), &
                               STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           end if
           if(allocated(spectrum.k)) then
                   deallocate(spectrum.k))
                   allocate(spectrum.k(nfreqs), &
                               STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           else
                   allocate(spectrum.k(nfreqs), &
                               STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           end if
           if(allocated(spectrum.dk)) then
                   deallocate(spectrum.dk)
                   allocate(spectrum.dk(nfreqs), &
                                STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           else
                   allocate(spectrum.dk(nfreqs), &
                                STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           end if
           if(allocated(spectrum.cp)) then
                   deallocate(spectrum.cp)
                   allocate(spectrum.cp(nfreqs), &
                                STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           else
                   allocate(spectrum.cp(nfreqs), &
                                STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           end if
           if(allocated(spectrum.cg)) then
                   deallocate(spectrum.cg)
                   allocate(spectrum.cg(nfreqs), &
                                STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           else
                   allocate(spectrum.cg(nfreqs), &
                                STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           end if
           if(allocated(spectrum.th)) then
                   deallocate(spectrum.th)
                   allocate(spectrum.th(ndirs), &
                                STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           else
                   allocate(spectrum.th(nfreqs), &
                                STAT=aerr,ERRMSG=emsg)
                   if(aerr /= 0) goto 9999
           end if
           allocate(tmp_f(nfreqs))
           allocate(tmp_k(nfreqs))
           tmp_f = 0.0_dp; tmp_k = 0.0_dp
           if(nfreqs == 1) then
                spectrum.f(1) = fmin
           else
                do concurrent(n=1:nfreqs)
                    spectrum.f(n) = exp(log(fmin)+(n-1)*log(df))
                end do
           end if
           do concurrent(n=1:nfreqs)
                spectrum.k(n) = wavenumber(spectrum.f(n),          &
                                           spectrum.depth,         &
                                           spectrum.water_density, &
                                           spectrum.grav,          &
                                           spectrum.surface_tension)
           end do
           spectrum.cp = pi2_const*spectrum.f/spectrum.k
           call diff(spectrum.f,tmp_f,nfreqs)
           call diff(spectrum.k,tmp_k,nfreqs)
           spectrum.cg = pi2_const*tmp_f/tmp_k

           do concurrent(n=1:ndirs)
                spectrum.th(n) = (n-0.5_dp*(ndirs+1))*2pi_const/ndirs
           end do
           spectrum.df = tmp_f
           spectrum.dk = tmp_k
           if(ndirs > 1) then
                allocate(th_tmp(ndirs))
                call diff_periodic(spectrum.th.th_tmp,ndirs)
                spectrum.dth = th_tmp
                spectrum.dth(1) = spectrum.dth(2)
                spectrum.dth(ndirs) = spectrum.dth(ndirs-1)
           else
                spectrum.dth = [1]
           end if
           call setCurrent2d(spectrum,[0.0_dp],[0.0_dp],[0.0_dp])
           return
9999
           call handle_fatal_memory_error(iounit, logging,verbose,append,fname, &
                     "logger: "// __FILE__ // "module: mod_spectrum, subroutine: constructor -- Memory Allocation Failure !!", &                                                        &
                              "module: mod_spectrum, subroutine: constructor -- Memory Allocation Failure !!", &
                              emsg,243)
    end subroutine constructor

!DIR$   ATTRIBUTES INLINE :: isAllocated
    pure elemental logical function isAllocated(spectrum)
!! Returns the allocation status of the spectrum array.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: isAllocated
          type(spectrum_type),      intent(in) :: spectrum
          ! Exec code ....
          isAllocated = allocated(spectrum.spec)
    end function isAllocated

!DIR$   ATTRIBUTES INLINE :: isMonochromatic
    pure elemental function isMonochromatic(spectrum) result(monochromatic)
  !! Returns `.true.` if only one frequency bin is allocated,
  !! and `.false.` otherwise.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: isMonochromatic
          type(spectrum_type),      intent(in) :: spectrum
          ! Locals
          logical(kind=int32_t) :: monochromatic
          ! Exec code ....
          if(size(spectrum.f) == 1) then
                monochromatic = .true.
          else
                monochromatic = .false.
          end if
    end function isMonochromatic

!DIR$   ATTRIBUTES INLINE :: isOmnidirectional
    pure elemental function isOmnidirectional(spectrum) result(omnidirect)
 !! Returns `.true.` if only one direction bin is allocated,
  !! and `.false.` otherwise.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: isOmnidirectional
          type(spectrum_type),      intent(in) :: spectrum
          ! LOcals
          logical(kind=int32_t) :: omnidirect
          ! Exec code ...
          if(size(spectrum.th) == 1) then
                omnidirect = .true.
          else
                omnidirect = .false.
          end if
    end function isOmnidirectional

    subroutine getFrequency2d(spectrum,f_out,nfreqs,ndirs)
  !!Extract frequency [Hz] array of the spectrum instance, reshaped to
  !! match the spectrum array shape. This method is most useful for conforming
  !! shape array in 2-d spectrum computations.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: getFRequency2d
          type(spectrum_type),      intent(in) :: spectrum
 !DIR$  ASSUME_ALIGNED f_out:64
          real(kind=dp),    dimension(nfreqs,ndirs), intent(out) :: f_out
          integer(kind=int32_t),    intent(in) :: nfreqs
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Locals
          integer(kind=int32_t), automatic :: ndir
          do concurrent(ndir=1:ndirs)
                f_out(:,ndir) = spectrum.f
          end do
    end subroutine getFrequency2d

    subroutine getWavenumber2d(spectrum,k_out,nfreqs,ndirs)
!! Returns the wavenumber [rad/m] array of the spectrum instance, reshaped to
  !! match the spectrum array shape. This method is most useful for conforming
  !! shape array in 2-d spectrum computations.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: getWavenumber2d
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED k_out:64
          real(kind=dp),    dimension(nfreqs,ndirs), intent(out) :: k_out
          integer(kind=int32_t),    intent(in) :: nfreqs
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Locals
          integer(kind=int32_t), automatic :: ndir
          ! EXec code ....
          do concurrent(ndir=1:ndirs)
                k_out(:,ndir) = spectrum.k
          end do
    end subroutine getWavenumber2d

!DIR$   ATTRIBUTES INLINE :: getWavelength
    subroutine getWavelength(spectrum,lambda,nfreqs)
          use mod_constants, only : pi2_const
!! Returns the wavelength [m] array of the spectrum instance.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: getWavelength
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED lambda:64
          real(kind=dp),    dimension(nfreqs), intent(out) :: lambda
          integer(kind=int32_t),    intent(in) ::nfreqs
          ! Exec code ....
          lambda = pi2_const/spectrum.k
    end subroutine getWavelength

    subroutine getDirection2d(spectrum,th_out,nfreqs,ndirs)
 !! Returns the directions [rad] array of the spectrum instance, reshaped to
  !! match the spectrum array shape. This method is most useful for conforming
  !! shape array in 2-d spectrum computations.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: getDirection2d
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED th_out:64
          real(kind=dp),    dimension(nfreqs,ndirs), intent(out) :: th_out
          integer(kind=int32_t),    intent(in) :: nfreqs
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Locals
          integer(kind=int32_t), automatic :: ndir
          ! EXec code ....
          do concurrent(ndir=1:ndirs)
                th_out(nfreq,:) = spectrum.th
          end do
    end subroutine getDirection2d

    subroutine getPhaseSpeed2d(spectrum,cp_out,nfreqs,ndirs)
  !! Returns the phase speed [m/s] array of the spectrum instance.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: getPhaseSpeed2d
          type(spectrum_type),      intent(in) :: spectrum
!DIR$ ASSUME_ALIGNED cp_out:64
          real(kind=dp),    dimension(nfreqs,ndirs),  intent(out) :: cp_out
          integer(kind=int32_t),    intent(in) :: nfreqs
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Locals
          integer(kind=int32_t), automatic :: ndir
          ! EXec code ....
          do concurrent(ndir=1:ndirs)
                cp_out(:,ndir) = spectrum.cp
          end do
    end subroutine getPhaseSpeed2d

    subroutine getGroupSpeed2d(spectrum,cg_out,nfreqs,ndirs)
 !! Returns the group speed [m/s] array of the spectrum instance.
 !DIR$  ATTRIBUTES CODE_ALIGN:32 :: getGroupSpeed2d
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED cg:64
          real(kind=dp), dimension(nfreqs,ndirs),   intent(out) :: cg_out
          integer(kind=int32_t),    intent(in) :: nfreqs
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Locals
          integer(kind=int32_t), automatic :: ndir
          ! EXec code ....
          do concurrent(ndir=1:ndirs)
                cg_out(:,ndir) = spectrum.cg
          end do
    end subroutine getGroupSpeed2d

!DIR$ ATTRIBUTES INLINE :: getWaveAction
    subroutine getWaveAction(spectrum,wave_action,freq2d,nfreqs,ndirs)
!! Returns the wave action spectrum, which corresponds to the the wave
  !! variance spectrum normalized by the intrinsic frequency.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: getWaveAction
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED wave_action:64
          real(kind=dp), dimension(nfreqs,ndirs), intent(in) :: wave_action
!DIR$   ASSUME_ALIGNED freq2d:64
          real(kind=dp), dimension(nfreqs,ndirs), intent(in) :: freq2d
          integer(kind=int32_t),    intent(in) :: nfreqs
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Exec code
          wave_action = spectrum.spec/freq2d
    end subroutine getWaveAction

!DIR$   ATTRIBUTES INLINE :: getAmplitude
    subroutine getAmplitude(spectrum,a,nfreqs,ndirs)
  !! Returns the amplitude array.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: getAmplitude
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED a:64
          real(kind=dp),    dimension(nfreqs,ndirs), intent(out) :: a
          integer(kind=int32_t),    intent(in) :: nfreqs
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Locals
          integer(kind=int32_t) :: ndir,ndirs2
          ! Exec code ...
          if(isMonochromatic(spectrum)) then
                a = sqrt(spectrum.spec+spectrum.spec)
          else
                ndirs2 = size(spectrum.th)
                do concurrent(ndir=1:ndirs2)
                         a(:,ndir) = sqrt((spectrum.spec(:,ndir)+spectrum.spec(:,ndir)) * &
                                           spectrum.df)
                end do
          end if

    end subroutine getAmplitude

!DIR$   ATTRIBUTES INLINE :: omnidirectionalSpectrum
    subroutine omnidirectionalSpectrum(spectrum,spec_out,nfreqs,ndirs)
 !! Returns the omnidirectional spectrum that corresponds to the input
  !! directional spectrum, integrated over all directions.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: omnidirectionalSpectrum
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED spec_out:64
          real(kind=dp),    dimension(nfreqs),  intent(out) :: spec_out
          integer(kind=int32_t),    intent(in) :: nfreqs
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Locals
          integer(kind=int32_t) :: ndir
!DIR$   VECTOR ALIGNED
!DIR$   SIMD
          do ndir = 1, ndirs
                spec_out(:) = spec_out(:)+spectrum.spec(:,ndir)*spectrum.dth(ndir)
          end do
    end subroutine omnidirectionalSpectrum



    subroutine setCurrent1d(spectrum,u,z)
!! Sets the 1-d current velocity field.
!DIR$   ATTRIBUTES CODE_ALIGN:32 ::  setCurrent1d
          type(spectrum_type),      intent(inout) :: spectrum
!DIR$   ASSUME_ALIGNED u:64
          real(kind=dp),    dimension(:), intent(in) :: u
!DIR$   ASSUME_ALIGNED z:64
          real(kind=dp),    dimension(:), intent(in) :: z
          ! Exec code ....
          spectrum.u = u
          spectrum.z = z
          allocate(spectrum.v(size(u))
          spectrum.v = 0.0_dp
    end subroutine setCurrent1d


    subroutine setCurrent2d(spectrum,u,v,z)
!! Sets the 2-d current velocity field.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: setCurrent2d
          type(spectrum_type),      intent(inout) :: spectrum
!DIR$   ASSUME_ALIGNED u:64
          real(kind=dp),    dimension(:), intent(in) :: u
!DIR$   ASSUME_ALIGNED v:64
          real(kind=dp),    dimension(:), intent(in) :: v
!DIR$   ASSUME_ALIGNED z:64
          real(kind=dp),    dimension(:), intent(in) :: z
          ! EXec code ....
          spectrum.u = u
          spectrum.v = v
          spectrum.z = z
    end subroutine setCurrent2d

!DIR$   ATTRIBUTES INLINE :: setDepth
    subroutine setDepth(spectrum,depth)
 !! Sets the mean surface elevation value.
 !DIR$  ATTRIBUTES CODE_ALIGN:32 :: setDepth
          type(spectrum_type),  intent(inout) :: spectrum
          real(kind=dp),        intent(in)    :: depth
          ! Exec code ...
          spectrum.depth = depth
    end subroutine setDepth

!DIR$   ATTRIBUTES INLINE :: setElevation
    subroutine setElevation(spectrum,elevation)
 !! Sets the mean surface elevation value.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: setElevation
          type(spectrum_type),   intent(inout) :: spectrum
          real(kind=dp),         intent(in)    :: elevation
          ! Exec code ...
          spectrum.elevation = elevation
    end subroutine setElevation

!DIR$   ATTRIBUTES INLINE :: setGravity
    subroutine setGravity(spectrum,gravity)
!! Sets the gravitational acceleration [m/s^2].
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: setGravity
          type(spectrum_type),    intent(inout) :: spectrum
          real(kind=dp),          intent(in)    :: gravity
          ! Exec code ....
          spectrum.gravity = gravity
    end subroutine setGravity

!DIR$   ATTRIBUTES INLINE :: setSurfaceTension
    subroutine setSurfaceTension(spectrum,surface_tension)
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: setSurfaceTension
          type(spectrum_type),  intent(inout) :: spectrum
          real(kind=dp),        intent(in)    :: surface_tension
          ! Exec code ...
          spectrum.surface_tension = surface_tension
    end subroutine setSurfaceTension

!DIR$   ATTRIBUTES INLINE :: setAirDensity
    subroutine setAirDensity(spectrum,air_density)
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: setAirDensity
          type(spectrum_type),  intent(inout) :: spectrum
          real(kind=dp),        intent(in)    :: air_density
          ! Exec code ....
          spectrum.air_density = air_density
    end subroutine setAirDensity

!DIR$   ATTRIBUTES INLINE :: setWaterDensity
    subroutine setWaterDensity(spectrum,water_density)
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: setWaterDensity
          type(spectrum_type),      intent(inout) :: spectrum
          real(kind=dp),            intent(in)    :: water_density
          ! Exec code ...
          spectrum.water_density = water_density
    end subroutine setWaterDensity

    subroutine meanSquareSlopeDirectional(spectrum,mss,wavenum_spec,dir_proj,  &
                                          nfreqs,ndirs)
!! For each directional frequency bin, computes the mean square slope of all
  !! all waves longer than that bin, projected to the direction of that bin.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: meanSquareSlopeDirectional
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED mss:64
          real(kind=dp),    dimension(nfreqs,ndirs), intent(out) :: mss
!DIR$   ASSUME_ALIGNED wavenum_spec:64
          real(kind=dp),    dimension(nfreqs,ndirs), intent(in)  :: wavenum_spec
!DIR$   ASSUME_ALIGNED dir_proj:64
          real(kind=dp),    dimension(ndirs,ndirs),  intent(in)  :: dir_proj
          integer(kind=int32_t),                     intent(in)  :: nfreqs
          integer(kind=int32_t),                     intent(in)  :: ndirs
          ! Locals
          integer(kind=int32_t), automatic :: ndir
          integer(kind=int32_t), automatic :: nfreq
          ! EXec code ....
          do concurrent(ndir=1:ndirs)
                dir_proj(:,ndir) = abs(cos(spectrum.th(ndir)-spectrum.th(ndir)))
          end do
          do ndir=1, ndirs
!DIR$   SIMD VECTORLENGTHFOR(REAL(KIND=dp))
             do nfreq=2,nfreqs
                 mss(nfreq,ndir) = mss(nfreq-1,ndir)                &
                 + sum(wavenum_spec(nfreq-1,:)*dir_proj(:,ndir)     &
                 * spectrum.k(nfreq-1)**2*spectrum.dk(nfreq-1)
             end do
          end do
    end subroutine meanSquareSlopeDirectional

!DIR$   ATTRIBUTES INLINE :: momentum_x
    pure elemental real(kind=dp) function momentum_x(spectrum)
!! Returns total wave momentum [kg/m/s] in x-direction.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: momentum_x
          type(spectrum_type),  intent(in) :: spectrum
          ! LOcals
          integer(kind=int32_t), automatic :: n,nfreqs
          ! Exec code ....
          nfreqs = size(spectrum.f)
          momentum_x = 0.0_dp
!DIR$ SIMD REDUCTION(+:momentum_x)
          do n=1, nfreqs
              momentum_x = momentum_x+sum(spectrum.spec(n,:)*spectrum.dth*cos(spectrum.th)) &
                * spectrum.df(n)/spectrum.cp(n)
          end do
          momentum_x = mometum_x*spectrum.water_density*spectrum.grav
    end function momentum_x

!DIR$   ATTRIBUTES INLINE :: momentum_y
    pure elemental real(kind=dp) function momentum_y(spectrum)
 !! Returns total wave momentum [kg/m/s] in y-direction.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: momentum_y
          type(spectrum_type),  intent(in) :: spectrum
          ! Locals
          integer(kind=int32_t), automatic :: n,nfreqs
          ! Exec code ....
          nfreqs = size(spectrum.f)
          momentum_y = 0.0_dp
!DIR$   SIMD REDUCTION(+:momentum_y)
          do n=1, nfreqs
               momentum_y = momentum_y+sum(spectrum.spec(n,:)*spectrum.dth*sin(spectrum.th)) &
                            * spectrum.df(n)/spectrum.cp(n)
          end do
          momentum_y = momentum_y*spectrum.water_density*spectrum.grav
    end function  momentum_y

!DIR$   ATTRIBUTES INLINE :: momentumFlux_xx
    pure elemental real(kind=dp) function momentumFlux_xx(spectrum)
!! Returns total advective flux [kg/m^2/s^2] in y-direction of momentum in
  !! y-direction.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: momentumFlux_xx
          type(spectrum_type),    intent(in) :: spectrum
          ! Locals
          integer(kind=int32_t), automatic :: n.nfreqs
          ! Exec code ....
          nfreqs = size(spectrum.f)
          momentumFlux_xx = 0.0_dp
!DIR$   SIMD REDUCTION(+:momentumFlux_xx)
          do n=1, nfreqs
               momentumFlux_xx = momentumFlux_xx            &
                + sum(spectrum.spec(n,:)*spectrum.dth*cos(spectrum.th)**2) &
                * spectrum.df(n)*spectrum.cg(n)/spectrum.cp(n)
          end do
          momentumFlux_xx = momentumFlux_xx*spectrum.water_density*spectrum.grav
    end function momentumFlux_xx

!DIR$   ATTRIBUTES INLINE :: momentumFlux_xy
    pure elemental real(kind=dp) function momentumFlux_xy(spectrum)
 !! Returns total advective flux [kg/m^2/s^2] in x-direction of momentum in
  !! y-direction and vice versa (flux in y-direction of momentum in
  !! y-direction), because \int{Cgx*My} == \int{Cgy*Mx}.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: momentumFlux_xy
          type(spectrum_type),      intent(in) :: spectrum
          ! Locals
          integer(kind=int32_t), automatic :: n,nfreqs
          ! Exec code ....
          nfreqs = size(spectrum.f)
          momentumFlux_xy = 0.0_dp
!DIR$   SIMD  REDUCTION(+:momentumFlux_xy)
          do n=1, nfreqs
              momentumFlux_xy = momentumFlux_xy             &
              + sum(spectrum.spec(n,:)*spectrum.dth*cos(spectrum.th)*sin(spectrum.th) &
              * spectrum.df(n)*spectrum.cg(n)/spectrum.cp(n)
          end do
          momentumFlux_xy = momentumFlux_xy*spectrum.water_density*spectrum.grav
    end function momentumFlux_xy

!DIR$   ATTRIBUTES INLINE :: frequencyMoment
    pure elemental real(kind=dp) function frequencyMoment(spectrum,n)
 !! Returns the spectral frequency moment of order n.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: frequencyMoment
          type(spectrum_type),   intent(in) :: spectrum
          integer(kind=int32_t), intent(in) :: n
          ! Exec code ....
          frequencyMoment = sum(spectrum.f**n*sum(spectrum.spec,dim=2)*spectrum.df)
    end function frequencyMoment

!DIR$   ATTRIBUTES INLINE :: peakedness
    pure elemental real(kind=dp) function peakedness(spectrum,spec_in,ndirs)
  !! Returns the peakedness parameter that quantifies the sharpness of the
  !! spectral peak, following Goda (1970).
  !!
  !! References:
  !!
  !! Goda, Y., 1970. Numerical experiments on waves statistics with spectral
  !! simulation. *Report. Port and Harbour Research Institute*, Japan, **9**,
  !! 3-57.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: peakedness
          use GMS_mod_constants, only : eps
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED spec_in:64
          real(kind=dp), dimension(ndirs), intent(in) :: spec_in
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Exec code ....
          peakedness = 2*sum(spectrum.f*spec_in**2*spectrum.df) &
                       / (frequencyMoment(spectrum,0)**2+eps)
    end function peakedness

!DIR$   ATTRIBUTES INLINE :: peakFrequency
    pure elemental real(kind=dp) function peakFrequency(spectrum,spec_in,ndirs)
 !! Returns the peak frequency based on Young (1995).
  !!
  !! References:
  !!
  !! Young, I, 1995. The determination of confidence limits associated with
  !! estimates of the spectral peak frequency. *Ocean Engng.*, **22**, 669-686.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: peakFrequency
          use GMS_mod_constants, only : eps
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED spec_in:64
          real(kind=dp),    dimension(ndirs), intent(in) :: spec_in
          integer(kind=int32_t),    intent(in) :: ndirs
          ! EXec code
          peakFrequency = sum(spectrum.f*spec_in**4*spectrum.df) &
                          /(sum(spec_in**4*spectrum.df)+eps)
    end function peakFrequency

!DIR$   ATTRIBUTES INLINE :: peakFrequencyDiscrete
    pure elemental real(kind=dp) function peakFrequencyDiscrete(spectrum,spec_in,ndirs)
!! Returns the peak frequency based on simple discrete maximum location of
  !! the spectrum array.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: peakFrequencyDiscrete
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED spec_in:64
          real(kind=dp),    dimension(ndirs), intent(in) :: spec_in
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Exec code
          peakFrequencyDiscrete = &
                spectrum.f(maxloc(spec_in,dim=1))
    end function peakFrequencyDiscrete



!DIR$   ATTRIBUTES INLINE :: momentumFlux_yy
    pure elemental real(kind=dp) function momentumFlux_yy(spectrum)
 !! Returns total advective flux [kg/m^2/s^2] in y-direction of momentum in
  !! y-direction.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: momentumFlux_yy
          type(spectrum_type),  intent(in) :: spectrum

          ! Locals
          integer(kind=int32_t), automatic :: n,nfreqs
          ! Exec code ....
          nfreqs = size(spectrum.f)
          momentumFlux_yy = 0.0_dp
!DIR$   SIMD  REDUCTION(+:momentumFlux_yy)
          do n=1, nfreqs
               momentumFlux_yy = momentumFlux_yy            &
               + sum(spectrum.spec(n,:)*spectrum.dth*sin(spectrum.th)**2) &
               * spectrum.df(n)*spectrum.cg(n)/spectrum.cp(n)
          end do
          momentumFlux_yy = momentumFlux_yy*spectrum.water_density*spectrum.grav
    end function momentumFlux_yy

    pure elemental function wavenumber(f,depth,water_density,grav, &
                                        surface_tension)
          use mod_constants, only : pi2_const,eps
!! Solves the linear water wave dispersion relationship using a
 !! Newton-Raphson iteration loop.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: wavenumber
          real(kind=dp),    intent(in) :: f
          real(kind=dp),    intent(in), optional :: depth
          real(kind=dp),    intent(in), optional :: water_density
          real(kind=dp),    intent(in), optional :: grav
          real(kind=dp),    intent(in), optional :: surface_tension
          ! LOcals
          real(kind=dp), automatic :: dk,b,fnd,t
          real(kind=dp) :: wavenumber
          integer(kind=int32_t), automatic :: counter
          ! EXec code ....
          associate(k=>wavenumber)
                fnd = pi2_const*f*sqrt(depth/grav)
                k   = fnd*fnd
                b = surface_tension/(water_density*grav*depth**2)
                counter = 1
                dk = 2.0e-3_dp
                newton_raphson: do
                        t = tanh(k)
                        dk = -(fnd**2-k*t*(1.0_dp+b*k**2)) &
                           /(3.0_dp*b*k**2*t+t+k*(1.0_dp+b*k**2)*(1.0_dp-t**2))
                        k = k-dk
                        if(abs(dk) < eps .or. counter > 100) then
                                exit newton_raphson
                        end if
                        counter = counter + 1
               end do newton_raphson
               k = k/depth
          end associate
    end function wavenumber

!DIR$ ATTRIBUTES INLINE :: wavenumberMoment
    pure elemental real(kind=dp) function wavenumberMoment(spectrum,spec_in,nfreqs,ndirs,n)
!! Returns the spectral wavenumber moment of order n.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: wavenumberMoment
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED spec_in:64
          real(kind=dp),    dimension(nfreqs,ndirs), intent(in) :: spec_in
          integer(kind=int32_t),    intent(in) :: nfreqs
          integer(kind=int32_t),    intent(in) :: ndirs
          integer(kind=int32_t),    intent(in) :: n
          ! Exec code ...
          wavenumberMoment = sum(spectrum.k**n*sum(spec_in,dim=2) &
                             * spectrum.dk
    end function wavenumberMoment

!DIR$   ATTRIBUTES INLINE :: significantWaveHeight
    pure elemental real(kind=dp) function significantWaveHeight(spectrum)
 !! Returns the significant wave height [m].
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: significantWaveHeight
          type(spectrum_type),      intent(in) :: spectrum
          ! Exec code ....
          significantWaveHeight =  4.0_dp*sqrt(frequencyMoment(spectrum,0))
    end function significantWaveHeight

!DIR$   ATTRIBUTES INLINE :: significantSurfaceOrbitalVelocity
    pure elemental real(kind=dp) function significantSurfaceOrbitalVelocity(spectrum) &
            result(uorb)
 !! Returns the significant surface orbital velocity [m/s].
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: significantSurfaceOrbitalVelocity
          GMS_use mod_constants, only : pi2_const
          type(spectrum_type),      intent(in) :: spectrum
          ! Exec code
          uorb = 2*sqrt(sum((pi2_const*spectrum.f)**2*sum(spectrum.spec,dim=2)*spectrum.df))
    end function significantSurfaceOrbitalVelocity

!DIR$   ATTRIBUTES INLINE :: meanPeriod
    pure elemental real(kind=p) function meanPeriod(spectrum)
 !! Returns the mean wave period [s].
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: meanPeriod
          use mod_constants, only : eps
          type(spectrum_type),      intent(in) :: spectrum
          ! Exec code ...
          meanPeriod = frequencyMoment(spectrum,0)/(frequencyMoment(spectrum,1)+eps)
    end function meanPeriod

!DIR$   ATTRIBUTES INLINE :: meanPeriodZeroCrossing
    pure elemental real(kind=dp) function meanPeriodZeroCrossing(spectrum)
 !! Returns the zero-crossing mean wave period [s]:
  !!
  !! Tm02 = \sqrt(m_0 / m_2)
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: meanPeriodZeroCrossing
          use GMS_mod_constants, only : eps
          type(spectrum_type),  intent(in) :: spectrum
          ! Exec code ...
          meanPeriodZeroCrossing = sqrt(frequencyMoment(spectrum,0) &
                                    /(frequencyMoment(2)+eps)
    end function meanPeriodZeroCrossing


    subroutine wavenumberSpectrum(spectrum,spec_out,nfreqs,ndirs)
          use GMS_mod_constants, only : v1_over_pi
!! Returns the wavenumber spectrum array of the spectrum instance.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: wavenumberSpectrum
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED spec_out:64
          real(kind=dp),    dimension(nfreqs,ndirs), intent(out) :: spec_out
          integer(kind=int32_t),    intent(in) :: nfreqs
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Locals
          integer(kind=int32_t) :: ndir
          ! Exec code
          do concurrent(ndir=1:ndirs)
                spec_out(:,ndir) =  spectrum.spec(nfreq,:)*spectrum.cg(nfreq)*v1_over_pi
          end do
    end subroutine wavenumberSpectrum

    subroutine saturationSpectrum(spectrum,spec_out,wavenum_spec,nfreqs,ndirs)
!! Returns the saturation spectrum B(k) = F(k)k^4.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: saturationSpectrum
          type(spectrum_type),      intent(in) :: spectrum
!DIR$   ASSUME_ALIGNED spec_out:64
          real(kind=dp),    dimension(nfreqs,ndirs), intent(out) :: spec_out
!DIR$   ASSUME_ALIGNED wavenum_spec:64
          real(kind=dp),    dimension(nfreqs,ndirs), intent(in)  :: wavenum_spec
          integer(kind=int32_t),    intent(in) :: nfreqs
          integer(kind=int32_t),    intent(in) :: ndirs
          ! Locals
          integer(kind=int32_t) :: n
          ! Exec code...
          do concurrent(n=1:ndirs)
                spec_out(:,n) = wavenum_spec(:,n)*spectrum.k**4
          end do
    end subroutine saturationSpectrum



end module GMS_mod_spectrum
