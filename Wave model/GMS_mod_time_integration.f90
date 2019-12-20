module GMS_mod_time_integration

!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD 3-clause license. See LICENSE for details.
! Modified by Bernard Gingold (beniekg@gmail.com) on 25/03/2019
    use GMS_mod_kinds,    only : int32_t, dp
    use GMS_mod_spectrum, only : spectrum_type
    use GMS_mod_domain,   only : domain_type
    implicit none

    private

    public :: backward_euler
    public :: exact_exponential
    public :: forward_euler
    public :: integrate

    interface integrate
         module procedure :: integrate_spectrum
         module procedure :: integrate_domain
    end interface integrate

   ! interface exact_exponential
   !     module procedure :: exact_exponential_spectrum
   !     module procedure :: exact_exponential_domain
   ! end interface exact_exponential
!
   ! interface backward_euler
   !     module procedure :: backward_euler_spectrum
   !     module procedure :: backward_euler_domain
   ! end interface backward_euler

   ! interface forward_euler
   !     module procedure :: forward_euler_spectrum
   !     module procedure :: forward_euler_domain
   ! end interface forward_euler

    contains

    subroutine integrate_spectrum(func,initial,tendency,spectrum,dt)
!! Integrates spectrum forward in time using a time integration method
  !! provided as the argument `func`.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: integrate_spectrum
          interface
                subroutine func(initial,tendency,spectrum,dt)
                    import :: spectrum_type, dp
                    type(spectrum_type),    intent(in)  :: initial
                    type(spectrum_type),    intent(in)  :: tendency
                    type(spectrum_type),    intent(out) :: spectrum
                    real(kind=dp),          intent(in)  :: dt
                end subroutine func
          end interface
          type(spectrum_type),      intent(in)  :: initial
          type(spectrum_type),      intent(in)  :: tendency
          type(spectrum_type),      intent(out) :: spectrum
          real(kind=dp),            intent(in)  :: dt
          ! EXec code ....
          call func(initial,tendency,spectrum,dt)
    end subroutine integrate_spectrum

    subroutine integrate_domain(func,initial,tendency,domain,dt)
!! Integrates domain forward in time using a time integration method
  !! provided as the argument `func`.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: integrate_domain
          interface
               subroutine func(initial,tendency,domain,dt)
                    import :: domain_type, dp
                    type(domain_type),      intent(in)  :: initial
                    type(domain_type),      intent(in)  :: tendency
                    type(domain_type),      intent(out) :: domain
                    real(kind=dp),          intent(in)  :: dt
               end subroutine func
          end interface
          type(domain_type),        intent(in)  :: initial
          type(domain_type),        intent(in)  :: tendency
          type(domain_type),        intent(out) :: domain
          real(kind=dp),            intent(in)  :: dt
          ! Exec code ....
          call func(initial,tendency,domain,dt)
    end subroutine integrate_domain


end module GMS_mod_time_integration
