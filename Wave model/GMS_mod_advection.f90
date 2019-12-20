module GMS_mod_advection

!
! wavy - A spectral ocean wave modeling and development framework
! Copyright (c) 2017, Wavebit Scientific LLC
! All rights reserved.
!
! Licensed under the BSD-3 clause license. See LICENSE for details.
! Modified by Bernard Gingold (beniekg@gmail.com) on 24/03/2019
! Removing Object-Oriented crap.
!===============================================================================
    use GMS_mod_kinds, only : int32_t, dp
    implicit none

    private
    ! Generic procedures
    public :: advectUpwind1stOrder
    public :: advectCentered2ndOrder

! Specific procedures
    public :: advectUpwind1stOrder1dRank1
    public :: advectUpwind1stOrder1dRank2
    public :: advectUpwind1stOrder2dRank1
    public :: advectUpwind1stOrder2dRank2
    public :: advectCentered2ndOrder1dRank1
    public :: advectCentered2ndOrder1dRank2

    interface advectUpwind1stOrder
        module procedure :: advectUpwind1stOrder1dRank0
        module procedure :: advectUpwind1stOrder1dRank1
        module procedure :: advectUpwind1stOrder1dRank2
        module procedure :: advectUpwind1stOrder2dRank0
        module procedure :: advectUpwind1stOrder2dRank1
        module procedure :: advectUpwind1stOrder2dRank2
    end interface advectUpwind1stOrder

    interface advectCentered2ndOrder
        module procedure :: advectCentered2ndOrder1dRank0
        module procedure :: advectCentered2ndOrder1dRank1
        module procedure :: advectCentered2ndOrder1dRank2
        module procedure :: advectCentered2ndOrder2dRank0
        module procedure :: advectCentered2ndOrder2dRank1
        module procedure :: advectCentered2ndOrder2dRank2
    end interface advectCentered2ndOrder

    contains


    subroutine advectUpwind1stOrder1dRank0(f,u,dx,dfdt)
  !! Computes the advective tendency of an input field f given the advective
  !! velocity field u [m/s] and grid spacing dx [m], using a first order,
  !! positive-definite upwind differencing. Fields f and u are defined on a
  !! semi-staggered Arakawa C-grid:
  !!
  !! --u---f---u---f---u---f--
  !!   |  i-1  |   i   |  i+1
  !!
  !! The differentiation range is [2,size(f)-1], so 1 halo cell on each end
  !! needs to be set priod to calling this function.
  !!
  !! This function is for 1-dimensional input arrays. It is overloaded by the
  !! `advectUpwind1stOrder1d` procedure.
!DIR$   ATTRIBUTES CODE_ALIGN:32 ::  advectUpwind1stOrder1dRank0
!DIR$   ASSUME_ALIGNED f:64
          real(kind=dp), dimension(:), intent(in) :: f
!DIR$   ASSUME_ALIGNED u:64
          real(kind=dp), dimension(:), intent(in) :: u
!DIR$   ASSUME_ALIGNED dx:64
          real(kind=dp), dimension(:), intent(in) :: dx
!DIR$   ASSUME_ALIGNED dfdt:64
          real(kind=dp), dimension(:), intent(out) :: dfdt

          ! Locals
          integer(kind=dint32_t), automatic :: i,idm
          ! EXec code ....
           idm = size(f)
           allocate(dfdt(2:idm-1))
           dfdt = 0.0_dp

          do concurrent(i=2:idm-1)
                  dfdt(i) =  - 0.5_dp*((u(i+1)+abs(u(i+1)))*f(i)  &
                             + (u(i+1)-abs(u(i+1)))*f(i+1)&
                             - (u(i)+abs(u(i)))*f(i-1)    &
                             - (u(i)-abs(u(i)))*f(i))     &
                             / dx(i)
           end do
    end subroutine advectUpwind1stOrder1dRank0

    subroutine  advectUpwind1stOrder1dRank1(f,u,dx,dfdt)
 !! Computes the advective tendency of an input field f given the advective
  !! velocity field u [m/s] and grid spacing dx [m], using a first order,
  !! positive-definite upwind differencing. Fields f and u are defined on a
  !! semi-staggered Arakawa C-grid:
  !!
  !! --u---f---u---f---u---f--
  !!   |  i-1  |   i   |  i+1
  !!
  !! The differentiation range is [2,size(f)-1], so 1 halo cell on each end
  !! needs to be set priod to calling this function.
  !!
  !! This function is for 2-dimensional input arrays. It is overloaded by the
  !! `advectUpwind1stOrder1d` procedure.
!DIR$   ATTRIBUTES CODE_ALIGN:32 ::  advectUpwind1stOrder1dRank1
!DIR$   ASSUME_ALIGNED f:64
          real(kind=dp), dimension(:,:), intent(in) :: f
!DIR$   ASSUME_ALIGNED u:64
          real(kind=dp), dimension(:,:), intent(in) :: u
!DIR$   ASSUME_ALIGNED dx:64
          real(kind=dp), dimension(:), intent(in) :: dx
!DIR$   ASSUME_ALIGNED dfdt:64
          real(kind=dp), dimension(:,:), intent(out) :: dfdt

          ! Locals
          integer(kind=int32_t), automatic ::  i,m,idm,mdm
          ! EXec code ....
          mdm = size(f,dim=1)
          idm = size(f,dim=2)

          allocate(dfdt(mdm,2:idm-1))
          dfdt = 0.0_dp ! Memset called here!
          do concurrent(i=2:idm-1,m=1:mdm)
                         dfdt(m,i) = - 0.5_dp*((u(m,i+1)+abs(u(m,i+1)))*f(m,i)   &
                                     + (u(m,i+1)-abs(u(m,i+1)))*f(m,i+1) &
                                     - (u(m,i)+abs(u(m,i)))*f(m,i-1)     &
                                     - (u(m,i)-abs(u(m,i)))*f(m,i))      &
                                     / dx(i)
         end do
    end subroutine advectUpwind1stOrder1dRank1

    subroutine advectUpwind1stOrder1dRank2(f,u,dx,dfdt)
  !! Computes the advective tendency of an input field f given the advective
  !! velocity field u [m/s] and grid spacing dx [m], using a first order,
  !! positive-definite upwind differencing. Fields f and u are defined on a
  !! semi-staggered Arakawa C-grid:
  !!
  !! The differentiation range is [2,size(f)-1], so 1 halo cell on each end
  !! needs to be set priod to calling this function.
  !!
  !! --u---f---u---f---u---f--
  !!   |  i-1  |   i   |  i+1
  !!
  !! This function is for 3-dimensional input arrays. It is overloaded by the
  !! `advectUpwind1stOrder1d` procedure.
!DIR$   ATTRIBUTES CODE_ALIGN:32 ::  advectUpwind1stOrder1dRank2
!DIR$   ASSUME_ALIGNED f:64
          real(kind=dp), dimension(:,:,:), intent(in) :: f
!DIR$   ASSUME_ALIGNED u:64
          real(kind=dp), dimension(:,:,:), intent(in) :: u
!DIR$   ASSUME_ALIGNED dx:64
          real(kind=dp), dimension(:), intent(in) :: dx
!DIR$   ASSUME_ALIGNED dfdt:64
          real(kind=dp), dimension(:,:,:), intent(out) :: dfdt
          ! LOcals
          integer(kind=int4) :: i,m,n

          !Exec code ...
          mdm = size(f,dim=1)
          ndm = size(f,dim=2)
          idm = size(f,dim=3)


          allocate(dfdt(mdm,ndm,2:idm-1))
          dfdt = 0.0_dp

          do concurrent(i = 2:idm-1,m = 1:mdm,n = 1:ndm)

                 dfdt(m,n,i) = - 0.5_dp*((u(m,n,i+1)+abs(u(m,n,i+1)))*f(m,n,i)  &
                               + (u(m,n,i+1)-abs(u(m,n,i+1)))*f(m,n,i+1)&
                               - (u(m,n,i)+abs(u(m,n,i)))*f(m,n,i-1)    &
                               - (u(m,n,i)-abs(u(m,n,i)))*f(m,n,i))     &
                               / dx(i)

           end do

    end subroutine advectUpwind1stOrder1dRank2

    subroutine advectUpwind1stOrder2dRank0(f,u,v,dx,dy,dfdt)
!! Computes the 2-d advective tendency of an input field f given the advective
  !! velocity field u and v [m/s] and grid spacing dx and dy [m], using a first
  !! order, positive-definite upwind differencing. Fields f, u, and v are
  !! defined on a semi-staggered Arakawa C-grid:
  !!
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j+1   u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !!  j    u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j-1   u   f   u   f   u   f
  !!       |  i-1  |   i   |  i+1
  !!
  !! The differentiation range is [2:size(f,dim=1)-1,2:size(f,dim=2)-1], so 1
  !! halo cell on each end needs to be set priod to calling this function.
  !!
  !! This function is for 1-dimensional input arrays. It is overloaded by the
  !! `advectUpwind1stOrder2d` procedure.
!DIR$   ATTRIBUTES CODE_ALIGN:32 ::  advectUpwind1stOrder2dRank0
!DIR$   ASSUME_ALIGNED f:64
    real(kind=dp),dimension(:,:),intent(in) :: f
    !! Input field to be advected
!DIR$   ASSUME_ALIGNED u:64
    real(kind=dp),dimension(:,:),intent(in) :: u
    !! Advective velocity in x-direction [m/s]
!DIR$   ASSUME_ALIGNED v:64
    real(kind=dp),dimension(:,:),intent(in) :: v
    !! Advective velocity in y-direction [m/s]
!DIR$   ASSUME_ALIGNED dx:64
    real(kind=dp),dimension(:,:),intent(in) :: dx
    !! Grid spacing in x-direction [m]
!DIR$   ASSUME_ALIGNED dy:64
    real(kind=dp),dimension(:,:),intent(in) :: dy
    !! Grid spacing in x-direction [m]

    real(kind=dp),dimension(:,:),allocatable :: dfdt
    !! Advective tendency

          integer(kind=int32_t) :: i,j
          integer(kind=int32_t) :: idm,jdm
!       EXEC CODE ....
          idm = size(f,dim=1)
          jdm = size(f,dim=2)

          allocate(dfdt(2:idm-1,2:jdm-1))
          dfdt = 0.0_dp

          do concurrent(i = 2:idm-1,j = 2:jdm-1)

                    dfdt(i,j) = - 0.25_dp                                     &
                         *((u(i+1,j)+abs(u(i+1,j)))*(dy(i,j)+dy(i+1,j))*f(i,j)   &
                         + (u(i+1,j)-abs(u(i+1,j)))*(dy(i,j)+dy(i+1,j))*f(i+1,j) &
                         - (u(i,j)+abs(u(i,j)))*(dy(i-1,j)+dy(i,j))*f(i-1,j)     &
                         - (u(i,j)-abs(u(i,j)))*(dy(i-1,j)+dy(i,j))*f(i,j)       &
                         + (v(i,j)+abs(v(i,j)))*(dx(i,j)+dx(i,j+1))*f(i,j)       &
                         + (v(i,j)-abs(v(i,j)))*(dx(i,j)+dx(i,j+1))*f(i+1,j)     &
                         - (v(i,j-1)+abs(v(i,j-1)))*(dx(i,j-1)+dx(i,j))*f(i-1,j) &
                         - (v(i,j-1)-abs(v(i,j-1)))*(dx(i,j-1)+dx(i,j))*f(i,j))  &
                         / (dx(i,j)*dy(i,j))
           end do

    end subroutine advectUpwind1stOrder2dRank0

   subroutine advectUpwind1stOrder2dRank1(f,u,v,dx,dy,dfdt)

  !! Computes the 2-d advective tendency of an input field f given the advective
  !! velocity field u and v [m/s] and grid spacing dx and dy [m], using a first
  !! order, positive-definite upwind differencing. Fields f, u, and v are
  !! defined on a semi-staggered Arakawa C-grid:
  !!
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j+1   u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !!  j    u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j-1   u   f   u   f   u   f
  !!       |  i-1  |   i   |  i+1
  !!
  !! The differentiation range is [2:size(f,dim=1)-1,2:size(f,dim=2)-1], so 1
  !! halo cell on each end needs to be set priod to calling this function.
  !!
  !! This function is for 1-dimensional input arrays. It is overloaded by the
  !! `advectUpwind1stOrder2d` procedure.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: advectUpwind1stOrder2dRank1
!DIR$   ASSUME_ALIGNED f:64
    real(kind=dp),dimension(:,:,:),intent(in) :: f
    !! Input field to be advected
!DIR$   ASSUME_ALIGNED u:64
    real(kind=dp),dimension(:,:,:),intent(in) :: u
    !! Advective velocity in x-direction [m/s]
!DIR$   ASSUME_ALIGNED v:64
    real(kind=dp),dimension(:,:,:),intent(in) :: v
    !! Advective velocity in y-direction [m/s]
!DIR$   ASSUME_ALIGNED dx:64
    real(kind=dp),dimension(:,:),intent(in) :: dx
    !! Grid spacing in x-direction [m]
!DIR$   ASSUME_ALIGNED dy:64
    real(kind=dp),dimension(:,:),intent(in) :: dy
    !! Grid spacing in x-direction [m]
!DIR$ ATTRIBUTES ALIGN : 64 :: dfdt
    real(kind=dp),dimension(:,:,:),allocatable :: dfdt
    !! Advective tendency

          integer(kind=int32_t) :: i,j,m
          integer(kind=int32_t) :: idm,jdm,mdm

          mdm = size(f,dim=1)
          idm = size(f,dim=2)
          jdm = size(f,dim=3)

          allocate(dfdt(mdm,2:idm-1,2:jdm-1))
          dfdt = 0.0_dp

          do concurrent(i = 2:idm-1,j = 2:jdm-1,m = 1:mdm)

                                     dfdt(m,i,j) = - 0.25_dp                                        &
                     *((u(m,i+1,j)+abs(u(m,i+1,j)))*(dy(i,j)+dy(i+1,j))*f(m,i,j)  &
                     + (u(m,i+1,j)-abs(u(m,i+1,j)))*(dy(i,j)+dy(i+1,j))*f(m,i+1,j)&
                     - (u(m,i,j)+abs(u(m,i,j)))*(dy(i-1,j)+dy(i,j))*f(m,i-1,j)    &
                     - (u(m,i,j)-abs(u(m,i,j)))*(dy(i-1,j)+dy(i,j))*f(m,i,j)      &
                     + (v(m,i,j)+abs(v(m,i,j)))*(dx(i,j)+dx(i,j+1))*f(m,i,j)      &
                     + (v(m,i,j)-abs(v(m,i,j)))*(dx(i,j)+dx(i,j+1))*f(m,i+1,j)    &
                     - (v(m,i,j-1)+abs(v(m,i,j-1)))*(dx(i,j-1)+dx(i,j))*f(m,i-1,j)&
                     - (v(m,i,j-1)-abs(v(m,i,j-1)))*(dx(i,j-1)+dx(i,j))*f(m,i,j)) &
                    / (dx(i,j)*dy(i,j))

        end do

    end subroutine  advectUpwind1stOrder2dRank1

   subroutine  advectUpwind1stOrder2dRank2(f,u,v,dx,dy,dfdt)

  !! Computes the 2-d advective tendency of an input field f given the advective
  !! velocity field u and v [m/s] and grid spacing dx and dy [m], using a first
  !! order, positive-definite upwind differencing. Fields f, u, and v are
  !! defined on a semi-staggered Arakawa C-grid:
  !!
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j+1   u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !!  j    u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j-1   u   f   u   f   u   f
  !!       |  i-1  |   i   |  i+1
  !!
  !! The differentiation range is [2:size(f,dim=1)-1,2:size(f,dim=2)-1], so 1
  !! halo cell on each end needs to be set priod to calling this function.
  !!
  !! This function is for 1-dimensional input arrays. It is overloaded by the
  !! `advectUpwind1stOrder2d` procedure.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: advectUpwind1stOrder2dRank2
!DIR$   ASSUME_ALIGNED f:64
    real(kind=dp),dimension(:,:,:,:),intent(in) :: f
    !! Input field to be advected
!DIR$   ASSUME_ALIGNED u:64
    real(kind=dp),dimension(:,:,:,:),intent(in) :: u
    !! Advective velocity in x-direction [m/s]
!DIR$   ASSUME_ALIGNED v:64
    real(kind=dp),dimension(:,:,:,:),intent(in) :: v
    !! Advective velocity in y-direction [m/s]
!DIR$   ASSUME_ALIGNED dx:64
    real(kind=dp),dimension(:,:),intent(in) :: dx
    !! Grid spacing in x-direction [m]
!DIR$   ASSUME_ALIGNED dy:64
    real(kind=dp),dimension(:,:),intent(in) :: dy
    !! Grid spacing in x-direction [m]

    real(kind=dp),dimension(:,:,:,:),allocatable :: dfdt
    !! Advective tendency

          integer(kind=int32_t) :: i,j,m,n
          integer(kind=int32_t) :: idm,jdm,mdm,ndm

          mdm = size(f,dim=1)
          ndm = size(f,dim=2)
          idm = size(f,dim=3)
          jdm = size(f,dim=4)

          allocate(dfdt(mdm,ndm,2:idm-1,2:jdm-1))
          dfdt = 0.0_dp

          do concurrent(i = 2:idm-1,j = 2:jdm-1,m = 1:mdm,n = 1:ndm)

                 dfdt(m,n,i,j) = - 0.25_dp                                            &
                      *((u(m,n,i+1,j)+abs(u(m,n,i+1,j)))*(dy(i,j)+dy(i+1,j))*f(m,n,i,j)  &
                      + (u(m,n,i+1,j)-abs(u(m,n,i+1,j)))*(dy(i,j)+dy(i+1,j))*f(m,n,i+1,j)&
                      - (u(m,n,i,j)+abs(u(m,n,i,j)))*(dy(i-1,j)+dy(i,j))*f(m,n,i-1,j)    &
                      - (u(m,n,i,j)-abs(u(m,n,i,j)))*(dy(i-1,j)+dy(i,j))*f(m,n,i,j)      &
                      + (v(m,n,i,j)+abs(v(m,n,i,j)))*(dx(i,j)+dx(i,j+1))*f(m,n,i,j)      &
                      + (v(m,n,i,j)-abs(v(m,n,i,j)))*(dx(i,j)+dx(i,j+1))*f(m,n,i+1,j)    &
                      - (v(m,n,i,j-1)+abs(v(m,n,i,j-1)))*(dx(i,j-1)+dx(i,j))*f(m,n,i-1,j)&
                      - (v(m,n,i,j-1)-abs(v(m,n,i,j-1)))*(dx(i,j-1)+dx(i,j))*f(m,n,i,j)) &
                         / (dx(i,j)*dy(i,j))

          end do

    end subroutine  advectUpwind1stOrder2dRank2

    subroutine advectCentered2ndOrder1dRank0(f,u,dx,dfdt)

  !! Computes the advective tendency of an input field f given the advective
  !! velocity field u [m/s] and grid spacing dx [m], using a second order
  !! centered differencing. Fields f and u are defined on a
  !! semi-staggered Arakawa C-grid:
  !!
  !! --u---f---u---f---u---f--
  !!   |  i-1  |   i   |  i+1
  !!
  !! The differentiation range is [2,size(f)-1], so 1 halo cell on each end
  !! needs to be set priod to calling this function.
  !!
  !! This function is for 1-dimensional input arrays. It is overloaded by the
  !! `advectCentered2ndOrder1d` procedure.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: advectCentered2ndOrder1dRank
!DIR$   ASSUME)ALIGNED f:64
          real(kind=dp),dimension(:),intent(in) :: f
    !! Input field to be advected
!DIR$   ASSUME_ALIGNED u:64
          real(kind=dp),dimension(:),intent(in) :: u
    !! Advective velocity [m/s]
!DIR$   ASSUME_ALIGNED dx:64
          real(kind=dp),dimension(:),intent(in) :: dx
    !! Grid spacing [m]
!DIR$ ATTRIBUTES ALIGN : 64 :: dfdt
          real(kind=dp),dimension(:),allocatable :: dfdt
    !! Advective tendency

            integer(kind=int32_t) :: i
            integer(kind=int32_t) :: idm

            idm = size(f)

            allocate(dfdt(2:idm-1))
            dfdt = 0.0_dp

            do concurrent(i = 2:idm-1)

                  dfdt(i) = - 0.5_dp*(u(i+1)*(f(i)+f(i+1))&
                            - u(i)*(f(i-1)+f(i))) &
                            / dx(i)

            end do

    end subroutine advectCentered2ndOrder1dRank0

   subroutine  advectCentered2ndOrder1dRank1(f,u,dx,dfdt)

  !! Computes the advective tendency of an input field f given the advective
  !! velocity field u [m/s] and grid spacing dx [m], using a second order
  !! centered differencing. Fields f and u are defined on a
  !! semi-staggered Arakawa C-grid:
  !!
  !! --u---f---u---f---u---f--
  !!   |  i-1  |   i   |  i+1
  !!
  !! The differentiation range is [2,size(f)-1], so 1 halo cell on each end
  !! needs to be set priod to calling this function.
  !!
  !! This function is for 2-dimensional input arrays. It is overloaded by the
  !! `advectCentered2ndOrder1d` procedure.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: advectCentered2ndOrder1dRank1
!DIR$   ASSUME_ALIGNED f:64
    real(kind=dp),dimension(:,:),intent(in) :: f
    !! Input field to be advected
!DIR$   ASSUME_ALIGNED u:64
    real(kind=dp),dimension(:,:),intent(in) :: u
    !! Advective velocity [m/s]
!DIR$   ASSUME_ALIGNED dx:64
    real(kind=dp),dimension(:),intent(in) :: dx
    !! Grid spacing [m]
!DIR$ ATTRIBUTES ALIGN : 64 :: dfdt
    real(kind=dp),dimension(:,:),allocatable :: dfdt
    !! Advective tendency

          integer(kind=int32_t) :: i,m
          integer(kind=int32_y) :: idm,mdm

          mdm = size(f,dim=1)
          idm = size(f,dim=2)

          allocate(dfdt(mdm,2:idm-1))
          dfdt = 0.0_dp

          do concurrent(i = 2:idm-1,m = 1:mdm)

                dfdt(m,i) = - 0.5_dp*(u(m,i+1)*(f(m,i)+f(m,i+1))&
                            - u(m,i)*(f(m,i-1)+f(m,i))) &
                           / dx(i)

          end do

    end subroutine  advectCentered2ndOrder1dRank1

   subroutine  advectCentered2ndOrder1dRank2(f,u,dx,dfdt)

  !! Computes the advective tendency of an input field f given the advective
  !! velocity field u [m/s] and grid spacing dx [m], using a second order
  !! centered differencing. Fields f and u are defined on a
  !! semi-staggered Arakawa C-grid:
  !!
  !! --u---f---u---f---u---f--
  !!   |  i-1  |   i   |  i+1
  !!
  !! The differentiation range is [2,size(f)-1], so 1 halo cell on each end
  !! needs to be set priod to calling this function.
  !!
  !! This function is for 3-dimensional input arrays. It is overloaded by the
  !! `advectCentered2ndOrder1d` procedure.
!DIR$ ATTRIBUTES CODE_ALIGN:32 :: advectCentered2ndOrder1dRank2

    real(kind=dp),dimension(:,:,:),intent(in) :: f
    !! Input field to be advected
    real(kind=dp),dimension(:,:,:),intent(in) :: u
    !! Advective velocity [m/s]
    real(kind=dp),dimension(:),intent(in) :: dx
    !! Grid spacing [m]
!DIR$ ATTRIBUTES ALIGN : 64 :: dfdt
    real(kind=dp),dimension(:,:,:),allocatable :: dfdt
    !! Advective tendency

          integer(kind=int32_t) :: i,m,n
          integer(kind=int32_t) :: idm,mdm,ndm

          mdm = size(f,dim=1)
          ndm = size(f,dim=2)
          idm = size(f,dim=3)

          allocate(dfdt(mdm,ndm,2:idm-1))
          dfdt = 0.0_dp

          do concurrent(i = 2:idm-1,m = 1:mdm,n = 1:ndm)

               dfdt(m,n,i) = - 0.5_dp*(u(m,n,i+1)*(f(m,n,i)+f(m,n,i+1))&
                             - u(m,n,i)*(f(m,n,i-1)+f(m,n,i))) &
                             / dx(i)

          end  do

    end subroutine advectCentered2ndOrder1dRank2

   subroutine  advectCentered2ndOrder2dRank0(f,u,v,dx,dy,dfdt)

  !! Computes the 2-d advective tendency of an input field f given the advective
  !! velocity field u and v [m/s] and grid spacing dx and dy [m], using a second
  !! order centered differencing. Fields f, u, and v are defined on a
  !! semi-staggered Arakawa C-grid:
  !!
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j+1   u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !!  j    u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j-1   u   f   u   f   u   f
  !!       |  i-1  |   i   |  i+1
  !!
  !! The differentiation range is [2:size(f,dim=1)-1,2:size(f,dim=2)-1], so 1
  !! halo cell on each end needs to be set priod to calling this function.
  !!
  !! This function is for 1-dimensional input arrays. It is overloaded by the
  !! `advectCentered2ndOrder2d` procedure.
!DIR$   ATTRIBUTES CODE_ALIGN:32 :: advectCentered2ndOrder2dRank0
!DIR$   ASSUME_ALIGNED f:64
          real(kind=dp),dimension(:,:),intent(in) :: f
    !! Input field to be advected
!DIR$   ASSUME_ALIGNED u:64
          real(kind=dp),dimension(:,:),intent(in) :: u
    !! Advective velocity in x-direction [m/s]
!DIR$   ASSUME_ALIGNED v:64
          real(kind=dp),dimension(:,:),intent(in) :: v
    !! Advective velocity in y-direction [m/s]
!DIR$   ASSUME_ALIGNED dx:64
          real(kind=dp),dimension(:,:),intent(in) :: dx
    !! Grid spacing in x-direction [m]
!DIR$   ASSUME_ALIGNED dy:64
          real(kind=dp),dimension(:,:),intent(in) :: dy
    !! Grid spacing in y-direction [m]
!DIR$ ATTRIBUTES ALIGN : 64 ::  dfdt
          real(kind=dp),dimension(:,:),allocatable :: dfdt
    !! Advective tendency

          integer(kind=int32_t) :: i,j
          integer(kind=int32_t) :: idm,jdm

          idm = size(f,dim=1)
          jdm = size(f,dim=2)

          allocate(dfdt(2:idm-1,2:jdm-1))
          dfdt = 0.0_dp

          do concurrent(i = 2:idm-1,j = 2:jdm-1)

                 dfdt(i,j) = - 0.25_dp*(u(i+1,j)*(dy(i,j)+dy(i+1,j))*(f(i,j)+f(i+1,j)) &
                             - u(i,j)*(dy(i-1,j)+dy(i,j))*(f(i-1,j)+f(i,j))   &
                             + v(i,j)*(dy(i,j)+dy(i,j+1))*(f(i,j)+f(i,j+1))   &
                             - v(i,j-1)*(dy(i,j-1)+dy(i,j))*(f(i,j-1)+f(i,j)))&
                             / (dx(i,j)*dy(i,j))

           end do
    end subroutine  advectCentered2ndOrder2dRank0


   subroutine  advectCentered2ndOrder2dRank1(f,u,v,dx,dy,dfdt)

  !! Computes the 2-d advective tendency of an input field f given the advective
  !! velocity field u and v [m/s] and grid spacing dx and dy [m], using a second
  !! order centered differencing. Fields f, u, and v are defined on a
  !! semi-staggered Arakawa C-grid:
  !!
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j+1   u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !!  j    u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j-1   u   f   u   f   u   f
  !!       |  i-1  |   i   |  i+1
  !!
  !! The differentiation range is [2:size(f,dim=1)-1,2:size(f,dim=2)-1], so 1
  !! halo cell on each end needs to be set priod to calling this function.
  !!
  !! This function is for 1-dimensional input arrays. It is overloaded by the
  !! `advectCentered2ndOrder2d` procedure.
!DIR$ ATTRIBUTES CODE_ALIGN:32 :: advectCentered2ndOrder2dRank1
!DIR$ ASSUME_ALIGNED f:64
          real(kind=dp),dimension(:,:,:),intent(in) :: f
    !! Input field to be advected
!DIR$ ASSUME_ALIGNED u:64
          real(kind=dp),dimension(:,:,:),intent(in) :: u
    !! Advective velocity in x-direction [m/s]
!DIR$ ASSUME_ALIGNED v:64
          real(kind=dp),dimension(:,:,:),intent(in) :: v
    !! Advective velocity in y-direction [m/s]
!DIR$ ASSUME_ALIGNED dx:64
          real(kind=dp),dimension(:,:),intent(in) :: dx
    !! Grid spacing in x-direction [m]
!DIR$ ASSUME_ALIGNED dy:64
          real(kind=dp),dimension(:,:),intent(in) :: dy
    !! Grid spacing in y-direction [m]
!DIR$ ATTRIBUTES ALIGN : 64 :: dfdt
          real(kind=dp),dimension(:,:,:),allocatable :: dfdt
    !! Advective tendency

          integer(kind=int32_t) :: i,j,m
          integer(kind=int32_t) :: idm,jdm,mdm

          mdm = size(f,dim=1)
          idm = size(f,dim=2)
          jdm = size(f,dim=3)

          allocate(dfdt(mdm,2:idm-1,2:jdm-1))
          dfdt = 0.0_dp

          do concurrent(i = 2:idm-1,j = 2:jdm-1,m = 1:mdm)

                dfdt(m,i,j) = - 0.25_dp                                  &
                              *(u(m,i+1,j)*(dy(i,j)+dy(i+1,j))*(f(m,i,j)+f(m,i+1,j)) &
                              - u(m,i,j)*(dy(i-1,j)+dy(i,j))*(f(m,i-1,j)+f(m,i,j))   &
                              + v(m,i,j)*(dy(i,j)+dy(i,j+1))*(f(m,i,j)+f(m,i,j+1))   &
                              - v(m,i,j-1)*(dy(i,j-1)+dy(i,j))*(f(m,i,j-1)+f(m,i,j)))&
                              / (dx(i,j)*dy(i,j))

          end do

    end subroutine advectCentered2ndOrder2dRank1

    subroutine  advectCentered2ndOrder2dRank2(f,u,v,dx,dy,dfdt)

  !! Computes the 2-d advective tendency of an input field f given the advective
  !! velocity field u and v [m/s] and grid spacing dx and dy [m], using a second
  !! order centered differencing. Fields f, u, and v are defined on a
  !! semi-staggered Arakawa C-grid:
  !!
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j+1   u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !!  j    u   f   u   f   u   f
  !!       |       |       |
  !!     --+---v---+---v---+---v--
  !!       |       |       |
  !! j-1   u   f   u   f   u   f
  !!       |  i-1  |   i   |  i+1
  !!
  !! The differentiation range is [2:size(f,dim=1)-1,2:size(f,dim=2)-1], so 1
  !! halo cell on each end needs to be set priod to calling this function.
  !!
  !! This function is for 1-dimensional input arrays. It is overloaded by the
  !! `advectCentered2ndOrder2d` procedure.
!DIR$ ATTRIBUTES CODE_ALIGN:32 :: advectCentered2ndOrder2dRank2
!DIR$   ASSUME_ALIGNED f:64
          real(kind=dp),dimension(:,:,:,:),intent(in) :: f
    !! Input field to be advected
!DIR$   ASSUME_ALIGNED u:64
          real(kind=dp),dimension(:,:,:,:),intent(in) :: u
    !! Advective velocity in x-direction [m/s]
!DIR$   ASSUME_ALIGNED v:64
          real(kind=dp),dimension(:,:,:,:),intent(in) :: v
    !! Advective velocity in y-direction [m/s]
!DIR$   ASSUME_ALIGNED dx:64
          real(kind=dp),dimension(:,:),intent(in) :: dx
    !! Grid spacing in x-direction [m]
!DIR$   ASSUME_ALIGNED dy:64
          real(kind=dp),dimension(:,:),intent(in) :: dy
    !! Grid spacing in y-direction [m]
!DIR$ ATTRIBUTES ALIGN : 64 :: dfdt
          real(kind=dp),dimension(:,:,:,:),allocatable :: dfdt
    !! Advective tendency

          integer(kind=int32_t) :: i,j,m,n
          integer(kind=int32_t) :: idm,jdm,mdm,ndm

          mdm = size(f,dim=1)
          ndm = size(f,dim=2)
          idm = size(f,dim=3)
          jdm = size(f,dim=4)

          allocate(dfdt(mdm,ndm,2:idm-1,2:jdm-1))
          dfdt = 0.0_dp

          do concurrent(i = 2:idm-1,j = 2:jdm-1,m = 1:mdm,n = 1:ndm)

                     dfdt(m,n,i,j) = - 0.25_dp                                      &
                          *(u(m,n,i+1,j)*(dy(i,j)+dy(i+1,j))*(f(m,n,i,j)+f(m,n,i+1,j)) &
                          - u(m,n,i,j)*(dy(i-1,j)+dy(i,j))*(f(m,n,i-1,j)+f(m,n,i,j))   &
                          + v(m,n,i,j)*(dy(i,j)+dy(i,j+1))*(f(m,n,i,j)+f(m,n,i,j+1))   &
                          - v(m,n,i,j-1)*(dy(i,j-1)+dy(i,j))*(f(m,n,i,j-1)+f(m,n,i,j)))&
                         / (dx(i,j)*dy(i,j))

          enddo

    end subroutine  advectCentered2ndOrder2dRank2

end module GMS_mod_advection
