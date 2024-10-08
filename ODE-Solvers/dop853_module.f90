!*****************************************************************************************
!> author: Jacob Williams
!
!  Modern Fortran Edition of the DOP853 ODE Solver.
!
!### See also
!   * [DOP853.f](http://www.unige.ch/~hairer/prog/nonstiff/dop853.f) -- The original code
!
!### History
!   * Jacob Williams : December 2015 : Created module from the DOP853 Fortran 77 code.
!   * Development continues at [GitHub](https://github.com/jacobwilliams/dop853).
!
!### License
!
!### Original DOP853 license:
!
!        Copyright (c) 2004, Ernst Hairer
!
!        Redistribution and use in source and binary forms, with or without
!        modification, are permitted provided that the following conditions are
!        met:
!
!        - Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.
!
!        - Redistributions in binary form must reproduce the above copyright
!        notice, this list of conditions and the following disclaimer in the
!        documentation and/or other materials provided with the distribution.
!
!        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!        IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!        TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!        PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR
!        CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!        EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!        PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!        PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!        LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!        NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!### License for updated version:
!
!        Modern Fortran Edition of the DOP853 ODE Solver
!        https://github.com/jacobwilliams/dop853
!
!        Copyright (c) 2015-2017, Jacob Williams
!        All rights reserved.
!
!        Redistribution and use in source and binary forms, with or without modification,
!        are permitted provided that the following conditions are met:
!
!        * Redistributions of source code must retain the above copyright notice, this
!          list of conditions and the following disclaimer.
!
!        * Redistributions in binary form must reproduce the above copyright notice, this
!          list of conditions and the following disclaimer in the documentation and/or
!          other materials provided with the distribution.
!
!        * The names of its contributors may not be used to endorse or promote products
!          derived from this software without specific prior written permission.
!
!        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!        ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!        WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
!        ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!        (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!        LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
!        ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!*****************************************************************************************

    module dop853_module
    use mod_kinds, only : i4,dp
    use dop853_constants
    use iso_fortran_env, only: output_unit,error_unit

    implicit none

    type,public :: dop853_class

        private

        !internal variables:
        integer(i4) :: n      = 0   !! the dimension of the system
        integer(i4) :: nfcn   = 0   !! number of function evaluations
        integer(i4) :: nstep  = 0   !! number of computed steps
        integer(i4) :: naccpt = 0   !! number of accepted steps
        integer(i4) :: nrejct = 0   !! number of rejected steps (due to error test),
                                !! (step rejections in the first step are not counted)
        integer(i4) :: nrdens = 0   !! number of components, for which dense output
                                !! is required. for `0 < nrdens < n` the components
                                !! (for which dense output is required) have to be
                                !! specified in `icomp(1),...,icomp(nrdens)`.
        real(dp) :: h = 0.0_dp  !! predicted step size of the last accepted step

        !input paramters:
        !  these parameters allow
        !  to adapt the code to the problem and to the needs of
        !  the user. set them on class initialization.

        integer(i4) :: iprint = output_unit   !! switch for printing error messages
                                          !! if `iprint==0` no messages are being printed
                                          !! if `iprint/=0` messages are printed with
                                          !! `write (iprint,*)` ...

        integer(i4) :: nmax = 100000       !! the maximal number of allowed steps.
        integer(i4) :: nstiff = 1000       !! test for stiffness is activated after step number
                                       !! `j*nstiff` (`j` integer(i4)), provided `nstiff>0`.
                                       !! for negative `nstiff` the stiffness test is
                                       !! never activated.
        real(dp) :: hinitial = 0.0_dp  !! initial step size, for `hinitial=0` an initial guess
                                       !! is computed with help of the function [[hinit]].
        real(dp) :: hmax = 0.0_dp      !! maximal step size, defaults to `xend-x` if `hmax=0`.
        real(dp) :: safe = 0.9_dp      !! safety factor in step size prediction
        real(dp) :: fac1 = 0.333_dp    !! parameter for step size selection.
                                       !! the new step size is chosen subject to the restriction
                                       !! `fac1 <= hnew/hold <= fac2`
        real(dp) :: fac2 = 6.0_dp      !! parameter for step size selection.
                                       !! the new step size is chosen subject to the restriction
                                       !! `fac1 <= hnew/hold <= fac2`
        real(dp) :: beta = 0.0_dp      !! is the `beta` for stabilized step size control
                                       !! (see section iv.2). positive values of beta ( <= 0.04 )
                                       !! make the step size control more stable.

        integer(i4),dimension(:),allocatable  :: icomp      !! `dimension(nrdens)`
                                                        !! the components for which dense output is required
        real(dp),dimension(:),allocatable :: cont       !! `dimension(8*nrdens)`
        !dir$ attributes align : 64 :: icomp
        !dir$ attributes align : 64 :: cont

        !formerly in the condo8 common block:
        real(dp) :: xold = 0.0_dp
        real(dp) :: hout = 0.0_dp

        !user-defined procedures:
        procedure(deriv_func),pointer  :: fcn    => null() !! subroutine computing the value of `f(x,y)`
        procedure(solout_func),pointer :: solout => null() !! subroutine providing the
                                                           !! numerical solution during integration.
                                                           !! if `iout>=1`, it is called during integration.

        contains

        private

        procedure,public :: initialize => set_parameters    !! initialization routine.
        procedure,public :: integrate  => dop853            !! main integration routine.
        procedure,public :: destroy    => destroy_dop853    !! destructor.
        procedure,public :: info       => get_dop853_info   !! to get info after a run.

        procedure :: dp86co
        procedure :: hinit
        procedure,public :: contd8  !! can be called in user's [[solout_func]] for dense output.

    end type dop853_class

    abstract interface

        subroutine deriv_func(me,x,y,f)
            !! subroutine computing the value of \( dy/dx = f(x,y) \)
            import :: dp,dop853_class
            implicit none
            class(dop853_class),intent(inout)   :: me
            real(dp),intent(in)                 :: x    !! independent variable \(x\)
            real(dp),dimension(:),intent(in)    :: y    !! state vector \( y(x) \) [size n]
            real(dp),dimension(:),intent(out)   :: f    !! derivative vector \( f(x,y) = dy/dx \) [size n]
        end subroutine deriv_func

        subroutine solout_func(me,nr,xold,x,y,irtrn,xout)
            !! `solout` furnishes the solution `y` at the `nr`-th
            !! grid-point `x` (thereby the initial value is
            !! the first grid-point).
            import :: dp,dop853_class
            implicit none
            class(dop853_class),intent(inout) :: me
            integer(i4),intent(in)                :: nr    !! grid point (0,1,...)
            real(dp),intent(in)               :: xold  !! the preceding grid point
            real(dp),intent(in)               :: x     !! current grid point
            real(dp),dimension(:),intent(in)  :: y     !! state vector \( y(x) \) [size n]
            integer(i4),intent(inout)             :: irtrn !! serves to interrupt the integration. if
                                                       !! `irtrn` is set `<0`, [[dop853]] will return to
                                                       !! the calling program. if the numerical solution
                                                       !! is altered in `solout`, set `irtrn = 2`.
            real(dp),intent(out)              :: xout  !! `xout` can be used for efficient intermediate output
                                                       !! if one puts `iout=3`. when `nr=1` define the first
                                                       !! output point `xout` in `solout`. the subroutine
                                                       !! `solout` will be called only when `xout` is in the
                                                       !! interval `[xold,x]`; during this call
                                                       !! a new value for `xout` can be defined, etc.
        end subroutine solout_func

    end interface

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Get info from a [[dop853_class]].

    subroutine get_dop853_info(me,n,nfcn,nstep,naccpt,nrejct,h)

    implicit none

    class(dop853_class),intent(in) :: me
    integer(i4),intent(out),optional   :: n      !! dimension of the system
    integer(i4),intent(out),optional   :: nfcn   !! number of function evaluations
    integer(i4),intent(out),optional   :: nstep  !! number of computed steps
    integer(i4),intent(out),optional   :: naccpt !! number of accepted steps
    integer(i4),intent(out),optional   :: nrejct !! number of rejected steps (due to error test),
                                             !! (step rejections in the first step are not counted)
    real(dp),intent(out),optional  :: h      !! predicted step size of the last accepted step

    if (present(n     )) n      = me%n
    if (present(nfcn  )) nfcn   = me%nfcn
    if (present(nstep )) nstep  = me%nstep
    if (present(naccpt)) naccpt = me%naccpt
    if (present(nrejct)) nrejct = me%nrejct
    if (present(h))      h      = me%h

    end subroutine get_dop853_info
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[dop853_class]].

    subroutine destroy_dop853(me)

    implicit none

    class(dop853_class),intent(out) :: me

    end subroutine destroy_dop853
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set the optional inputs for [[dop853]].
!
!@note In the original code, these were part of the `work` and `iwork` arrays.

    subroutine set_parameters(me,n,fcn,solout,iprint,nstiff,nmax,hinitial,&
                              hmax,safe,fac1,fac2,beta,icomp,status_ok)

    implicit none

    class(dop853_class),intent(inout)        :: me
    integer(i4),intent(in)                       :: n         !! the dimension of the system (size of \(y\) and \(y'\) vectors)
    procedure(deriv_func)                    :: fcn       !! subroutine computing the value of \( y' = f(x,y) \)
    procedure(solout_func),optional          :: solout    !! subroutine providing the
                                                          !! numerical solution during integration.
                                                          !! if `iout>=1`, it is called during integration.
                                                          !! supply a dummy subroutine if `iout=0`.
    integer(i4),intent(in),optional              :: iprint    !! switch for printing error messages
                                                          !! if `iprint==0` no messages are being printed
                                                          !! if `iprint/=0` messages are printed with
                                                          !! `write (iprint,*)` ...
    integer(i4),intent(in),optional              :: nstiff    !! test for stiffness is activated after step number
                                                          !! `j*nstiff` (`j` integer(i4)), provided `nstiff>0`.
                                                          !! for negative `nstiff` the stiffness test is
                                                          !! never activated.
    integer(i4),intent(in),optional              :: nmax      !! the maximal number of allowed steps.
    real(dp),intent(in),optional             :: hinitial  !! initial step size, for `hinitial=0` an initial guess
                                                          !! is computed with help of the function [[hinit]].
    real(dp),intent(in),optional             :: hmax      !! maximal step size, defaults to `xend-x` if `hmax=0`.
    real(dp),intent(in),optional             :: safe      !! safety factor in step size prediction
    real(dp),intent(in),optional             :: fac1      !! parameter for step size selection.
                                                          !! the new step size is chosen subject to the restriction
                                                          !! `fac1 <= hnew/hold <= fac2`
    real(dp),intent(in),optional             :: fac2      !! parameter for step size selection.
                                                          !! the new step size is chosen subject to the restriction
                                                          !! `fac1 <= hnew/hold <= fac2`
    real(dp),intent(in),optional             :: beta      !! is the `beta` for stabilized step size control
                                                          !! (see section iv.2). positive values of `beta` ( <= 0.04 )
                                                          !! make the step size control more stable.
    integer(i4),dimension(:),intent(in),optional :: icomp     !! the components for which dense output is required (size from 0 to `n`).
    !dir$ assume_aligned icomp:64
    logical,intent(out)                      :: status_ok !! will be false for invalid inputs.

    call me%destroy()

    status_ok = .true.

    !required inputs:
    me%n = n
    me%fcn => fcn

    !optional inputs:

    if (present(solout)) me%solout => solout

    if (present(iprint))   me%iprint   = iprint
    if (present(nstiff))   me%nstiff   = nstiff
    if (present(hinitial)) me%hinitial = hinitial
    if (present(hmax))     me%hmax     = hmax
    if (present(fac1))     me%fac1     = fac1
    if (present(fac2))     me%fac2     = fac2

    if (present(nmax)) then
        if ( nmax<=0 ) then
            if ( me%iprint/=0 ) &
                write (me%iprint,*) ' wrong input nmax=', nmax
            status_ok = .false.
        else
            me%nmax = nmax
        end if
    end if

    if (present(safe)) then
        if ( safe>=1.0_dp .or. safe<=1.0e-4_dp ) then
            if ( me%iprint/=0 ) &
                write (me%iprint,*) ' curious input for safety factor safe:', &
                                    safe
            status_ok = .false.
        else
            me%safe = safe
        end if
    end if

    if (present(beta)) then
        if ( beta<=0.0_dp ) then
           me%beta = 0.0_dp
        else
           if ( beta>0.2_dp ) then
              if ( me%iprint/=0 ) write (me%iprint,*) &
                                    ' curious input for beta: ', beta
              status_ok = .false.
           else
              me%beta = beta
           end if
        end if
    end if

    if (present(icomp)) then
        me%nrdens = size(icomp)
        !check validity of icomp array:
        if (size(icomp)<=me%n .and. all(icomp>0 .and. icomp<=me%n)) then
            allocate(me%icomp(me%nrdens));  me%icomp = icomp
            allocate(me%cont(8*me%nrdens)); me%cont = 0.0_dp
        else
            if ( me%iprint/=0 ) write (me%iprint,*) &
                                  ' invalid icomp array: ',icomp
            status_ok = .false.
        end if
    end if

    end subroutine set_parameters
!*****************************************************************************************

!*****************************************************************************************
!>
!  Numerical solution of a system of first order
!  ordinary differential equations \( y'=f(x,y) \).
!  This is an explicit Runge-Kutta method of order 8(5,3)
!  due to Dormand & Prince (with stepsize control and
!  dense output).
!
!### Authors
!  * E. Hairer and G. Wanner
!    Universite de Geneve, Dept. De Mathematiques
!    ch-1211 geneve 24, switzerland
!    e-mail:  ernst.hairer@unige.ch
!             gerhard.wanner@unige.ch
!  * Version of October 11, 2009
!    (new option `iout=3` for sparse dense output)
!  * Jacob Williams, Dec 2015: significant refactoring into modern Fortran.
!
!### Reference
!  * E. Hairer, S.P. Norsett and G. Wanner, [Solving Ordinary
!    Differential Equations I. Nonstiff Problems. 2nd Edition](http://www.unige.ch/~hairer/books.html).
!    Springer Series in Computational Mathematics, Springer-Verlag (1993)

      subroutine dop853(me,x,y,xend,rtol,atol,iout,idid)

      implicit none

      class(dop853_class),intent(inout)       :: me
      real(dp),intent(inout)                  :: x      !! *input:* initial value of independent variable.
                                                        !! *output:* `x` for which the solution has been computed
                                                        !! (after successful return `x=xend`).
      real(dp),dimension(:),intent(inout)     :: y      !! *input:* initial values for `y`. [size n]
      !dir$ assume_aligned y:64                                                  
                                                        !! *output:* numerical solution at `x`.
      real(dp),intent(in)                     :: xend   !! final x-value (`xend`-`x` may be positive or negative)
      real(dp),dimension(:),intent(in)        :: rtol   !! relative error tolerance. `rtol` and `atol`
                                                        !! can be both scalars or else both vectors of length `n`.
      !dir$ assume_aligned rtol:64
      real(dp),dimension(:),intent(in)        :: atol   !! absolute error tolerance. `rtol` and `atol`
                                                        !! can be both scalars or else both vectors of length `n`.
                                                        !! `atol` should be strictly positive (possibly very small)
      !dir$ assume_aligned atol:64
      integer(i4),intent(in)                      :: iout   !! switch for calling the subroutine `solout`:
                                                        !!
                                                        !! * `iout=0`: subroutine is never called
                                                        !! * `iout=1`: subroutine is called after every successful step
                                                        !! * `iout=2`: dense output is performed after every successful step
                                                        !! * `iout=3`: dense output is performed in steps defined by the user
                                                        !!          (see `xout` above)
      integer(i4),intent(out)                     :: idid   !! reports on successfulness upon return:
                                                        !!
                                                        !! * `idid=1`  computation successful,
                                                        !! * `idid=2`  comput. successful (interrupted by [[solout]]),
                                                        !! * `idid=-1` input is not consistent,
                                                        !! * `idid=-2` larger `nmax` is needed,
                                                        !! * `idid=-3` step size becomes too small.
                                                        !! * `idid=-4` problem is probably stiff (interrupted).

      real(dp) :: beta,fac1,fac2,h,hmax,safe
      integer(i4)  :: i,ieco,iprint,istore,nrdens,nstiff,nmax
      logical  :: arret
      integer(i4)  :: itol    !! switch for `rtol` and `atol`:
                          !!
                          !! * `itol=0`: both `rtol` and `atol` are scalars.
                          !!    the code keeps, roughly, the local error of
                          !!    `y(i)` below `rtol*abs(y(i))+atol`.
                          !!
                          !! * `itol=1`: both `rtol` and `atol` are vectors.
                          !!    the code keeps the local error of `y(i)` below
                          !!    `rtol(i)*abs(y(i))+atol(i)`.

      iprint = me%iprint
      arret = .false.

      !check procedures:
      if (.not. associated(me%fcn)) then
          if ( iprint/=0 ) &
                write (iprint,*) &
                'Error in dop853: procedure FCN is not associated.'
          idid = -1
          return
      end if

      if (iout/=0 .and. .not. associated(me%solout)) then
          if ( iprint/=0 ) &
                write (iprint,*) &
                'Error in dop853: procedure SOLOUT must be associated if IOUT/=0.'
          idid = -1
          return
      end if

      !scalar or vector tolerances:
      if (size(rtol)==1 .and. size(atol)==1) then
          itol = 0
      elseif (size(rtol)==me%n .and. size(atol)==me%n) then
          itol = 1
      else
          if ( iprint/=0 ) &
                write (iprint,*) &
                'Error in dop853: improper dimensions for rtol and/or atol.'
          idid = -1
          return
      end if

      ! setting the parameters
      me%nfcn   = 0
      me%nstep  = 0
      me%naccpt = 0
      me%nrejct = 0

      nmax = me%nmax
      nrdens = me%nrdens  !number of dense output components

      ! nstiff parameter for stiffness detection
      if ( me%nstiff<=0 ) then
          nstiff = nmax + 10  !no stiffness check
      else
          nstiff = me%nstiff
      end if

      if ( nrdens<0 .or. me%nrdens>me%n ) then
         if ( iprint/=0 ) write (iprint,*) ' curious input nrdens=' , nrdens
         arret = .true.
      else
         if ( nrdens>0 .and. iout<2 .and. iprint/=0 ) &
                write (iprint,*) ' warning: set iout=2 or iout=3 for dense output '
      end if

      if (size(y)/=me%n) then
          if ( iprint/=0 ) &
            write (iprint,*) ' error: y must have n elements: size(y)= ',size(y)
          arret = .true.
      end if

      safe = me%safe
      fac1 = me%fac1
      fac2 = me%fac2
      beta = me%beta

      if ( me%hmax==0.0_dp ) then
          hmax = xend - x
      else
          hmax = me%hmax
      end if

      h = me%hinitial     ! initial step size
      me%h = h

      ! when a fail has occurred, we return with idid=-1
      if ( arret ) then

         idid = -1

      else

        ! call to core integrator
        call me%dp86co(x,y,xend,hmax,h,rtol,atol,itol,iprint, &
                       iout,idid,nmax,nstiff,safe,beta,fac1,fac2, &
                       me%nfcn,me%nstep,me%naccpt,me%nrejct)

        me%h = h  ! may have been updated

      end if

      end subroutine dop853
!*****************************************************************************************

!*****************************************************************************************
!>
!  Core integrator for [[dop853]].
!  parameters same as in [[dop853]] with workspace added.

    subroutine dp86co(me,x,y,xend,hmax,h,rtol,atol,itol,iprint, &
                        iout,idid,nmax,nstiff,safe, &
                        beta,fac1,fac2, &
                        nfcn,nstep,naccpt,nrejct)

    implicit none

    class(dop853_class),intent(inout)   :: me
    real(dp),intent(inout)              :: x
    real(dp),dimension(:),intent(inout) :: y
    !dir$ assume_aligned y:64
    real(dp),intent(in)                 :: xend
    real(dp),intent(inout)              :: hmax
    real(dp),intent(inout)              :: h
    real(dp),dimension(:),intent(in)    :: rtol
    real(dp),dimension(:),intent(in)    :: atol
    !dir$ assume_aligned rtol:64
    !dir$ assume_aligned atol:64
    integer(i4),intent(in)                  :: itol
    integer(i4),intent(in)                  :: iprint
    integer(i4),intent(in)                  :: iout
    integer(i4),intent(out)                 :: idid
    integer(i4),intent(in)                  :: nmax
    integer(i4),intent(in)                  :: nstiff
    real(dp),intent(in)                 :: safe
    real(dp),intent(in)                 :: beta
    real(dp),intent(in)                 :: fac1
    real(dp),intent(in)                 :: fac2
    integer(i4),intent(inout)               :: nfcn
    integer(i4),intent(inout)               :: nstep
    integer(i4),intent(inout)               :: naccpt
    integer(i4),intent(inout)               :: nrejct

    real(dp),dimension(me%n) :: y1,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10
    real(dp) :: atoli,bspl,deno,err,err2,erri,expo1,fac,fac11,&
                facc1,facc2,facold,hlamb,hnew,posneg,rtoli,&
                sk,stden,stnum,xout,xph,ydiff
    integer(i4) :: i,iasti,iord,irtrn,j,nonsti,nrd
    logical :: reject,last,event,abort

    ! initialisations
    nrd = me%nrdens
    facold = 1.0e-4_dp
    expo1 = 1.0_dp/8.0_dp - beta*0.2_dp
    facc1 = 1.0_dp/fac1
    facc2 = 1.0_dp/fac2
    posneg = sign(1.0_dp,xend-x)

    ! initial preparations
    atoli = atol(1)
    rtoli = rtol(1)
    last = .false.
    hlamb = 0.0_dp
    iasti = 0
    call me%fcn(x,y,k1)
    hmax = abs(hmax)
    iord = 8
    if ( h==0.0_dp ) then
        h = me%hinit(x,y,posneg,k1,iord,hmax,atol,rtol,itol)
    end if
    nfcn = nfcn + 2
    reject = .false.
    me%xold = x
    if ( iout/=0 ) then
        irtrn = 1
        me%hout = 1.0_dp
        call me%solout(naccpt+1,me%xold,x,y,irtrn,xout)
        abort = ( irtrn<0 )
    else
        abort = .false.
    end if

    if (.not. abort) then
        do
            ! basic integration step
            if ( nstep>nmax ) then
                if ( iprint/=0 ) &
                        write (iprint,'(A,E18.4)') ' exit of dop853 at x=', x
                if ( iprint/=0 ) &
                        write (iprint,*) ' more than nmax =' , nmax , 'steps are needed'
                idid = -2
                return
            elseif ( 0.1_dp*abs(h)<=abs(x)*uround ) then
                if ( iprint/=0 ) &
                        write (iprint,'(A,E18.4)') ' exit of dop853 at x=', x
                if ( iprint/=0 ) &
                        write (iprint,*) ' step size too small, h=' , h
                idid = -3
                return
            else
                if ( (x+1.01_dp*h-xend)*posneg>0.0_dp ) then
                    h = xend - x
                    last = .true.
                end if
                nstep = nstep + 1
                ! the twelve stages
                if ( irtrn>=2 ) call me%fcn(x,y,k1)
                y1 = y + h*a21*k1
                call me%fcn(x+c2*h,y1,k2)
                y1 = y + h*(a31*k1+a32*k2)
                call me%fcn(x+c3*h,y1,k3)
                y1 = y + h*(a41*k1+a43*k3)
                call me%fcn(x+c4*h,y1,k4)
                y1 = y + h*(a51*k1+a53*k3+a54*k4)
                call me%fcn(x+c5*h,y1,k5)
                y1 = y + h*(a61*k1+a64*k4+a65*k5)
                call me%fcn(x+c6*h,y1,k6)
                y1 = y + h*(a71*k1+a74*k4+a75*k5+a76*k6)
                call me%fcn(x+c7*h,y1,k7)
                y1 = y + h*(a81*k1+a84*k4+a85*k5+a86*k6+a87*k7)
                call me%fcn(x+c8*h,y1,k8)
                y1 = y + h*(a91*k1+a94*k4+a95*k5+a96*k6+a97*k7+a98*k8)
                call me%fcn(x+c9*h,y1,k9)
                y1 = y + h*(a101*k1+a104*k4+a105*k5+a106*k6+a107*k7+&
                            a108*k8+a109*k9)
                call me%fcn(x+c10*h,y1,k10)
                y1 = y + h*(a111*k1+a114*k4+a115*k5+a116*k6+a117*k7+&
                            a118*k8+a119*k9+a1110*k10)
                call me%fcn(x+c11*h,y1,k2)
                xph = x + h
                y1 = y + h*(a121*k1+a124*k4+a125*k5+a126*k6+a127*k7+&
                            a128*k8+a129*k9+a1210*k10+a1211*k2)
                call me%fcn(xph,y1,k3)
                nfcn = nfcn + 11
                k4 = b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k2+b12*k3
                k5 = y + h*k4
                ! error estimation
                err = 0.0_dp
                err2 = 0.0_dp
                if ( itol==0 ) then
                    do i = 1 , me%n
                        sk = atoli + rtoli*max(abs(y(i)),abs(k5(i)))
                        erri = k4(i) - bhh1*k1(i) - bhh2*k9(i) - bhh3*k3(i)
                        err2 = err2 + (erri/sk)**2
                        erri = er1*k1(i) + er6*k6(i) + er7*k7(i) + er8*k8(i) &
                                + er9*k9(i) + er10*k10(i) + er11*k2(i) &
                                + er12*k3(i)
                        err = err + (erri/sk)**2
                    end do
                else
                    do i = 1 , me%n
                        sk = atol(i) + rtol(i)*max(abs(y(i)),abs(k5(i)))
                        erri = k4(i) - bhh1*k1(i) - bhh2*k9(i) - bhh3*k3(i)
                        err2 = err2 + (erri/sk)**2
                        erri = er1*k1(i) + er6*k6(i) + er7*k7(i) + er8*k8(i) &
                                + er9*k9(i) + er10*k10(i) + er11*k2(i) &
                                + er12*k3(i)
                        err = err + (erri/sk)**2
                    end do
                end if
                deno = err + 0.01_dp*err2
                if ( deno<=0.0_dp ) deno = 1.0_dp
                err = abs(h)*err*sqrt(1.0_dp/(me%n*deno))
                ! computation of hnew
                fac11 = err**expo1
                ! lund-stabilization
                fac = fac11/facold**beta
                ! we require  fac1 <= hnew/h <= fac2
                fac = max(facc2,min(facc1,fac/safe))
                hnew = h/fac
                if ( err<=1.0_dp ) then
                    ! step is accepted
                    facold = max(err,1.0e-4_dp)
                    naccpt = naccpt + 1
                    call me%fcn(xph,k5,k4)
                    nfcn = nfcn + 1
                    ! stiffness detection
                    if ( mod(naccpt,nstiff)==0 .or. iasti>0 ) then
                        stnum = 0.0_dp
                        stden = 0.0_dp
                        do i = 1 , me%n
                            stnum = stnum + (k4(i)-k3(i))**2
                            stden = stden + (k5(i)-y1(i))**2
                        end do
                        if ( stden>0.0_dp ) hlamb = abs(h)*sqrt(stnum/stden)
                        if ( hlamb>6.1_dp ) then
                            nonsti = 0
                            iasti = iasti + 1
                            if ( iasti==15 ) then
                                if ( iprint/=0 ) &
                                        write (iprint,*) &
                                        ' the problem seems to become stiff at x = ', x
                                if ( iprint==0 ) then
                                    idid = -4      ! fail exit
                                    return
                                end if
                            end if
                        else
                            nonsti = nonsti + 1
                            if ( nonsti==6 ) iasti = 0
                        end if
                    end if
                    ! final preparation for dense output
                    event = (iout==3) .and. (xout<=xph)
                    if ( iout==2 .or. event ) then
                        ! save the first function evaluations
                        do j = 1 , nrd
                            i = me%icomp(j)
                            me%cont(j) = y(i)
                            ydiff = k5(i) - y(i)
                            me%cont(j+nrd) = ydiff
                            bspl = h*k1(i) - ydiff
                            me%cont(j+nrd*2) = bspl
                            me%cont(j+nrd*3) = ydiff - h*k4(i) - bspl
                            me%cont(j+nrd*4) = d41*k1(i)+d46*k6(i)+d47*k7(i)+&
                                               d48*k8(i)+d49*k9(i)+d410*k10(i)+&
                                               d411*k2(i)+d412*k3(i)
                            me%cont(j+nrd*5) = d51*k1(i)+d56*k6(i)+d57*k7(i)+&
                                               d58*k8(i)+d59*k9(i)+d510*k10(i)+&
                                               d511*k2(i)+d512*k3(i)
                            me%cont(j+nrd*6) = d61*k1(i)+d66*k6(i)+d67*k7(i)+&
                                               d68*k8(i)+d69*k9(i)+d610*k10(i)+&
                                               d611*k2(i)+d612*k3(i)
                            me%cont(j+nrd*7) = d71*k1(i)+d76*k6(i)+d77*k7(i)+&
                                               d78*k8(i)+d79*k9(i)+d710*k10(i)+&
                                               d711*k2(i)+d712*k3(i)
                        end do
                        ! the next three function evaluations
                        y1 = y + h*(a141*k1+a147*k7+a148*k8+a149*k9+&
                                    a1410*k10+a1411*k2+a1412*k3+a1413*k4)
                        call me%fcn(x+c14*h,y1,k10)
                        y1 = y + h*(a151*k1+a156*k6+a157*k7+a158*k8+&
                                    a1511*k2+a1512*k3+a1513*k4+a1514*k10)
                        call me%fcn(x+c15*h,y1,k2)
                        y1 = y + h*(a161*k1+a166*k6+a167*k7+a168*k8+a169*k9+&
                                    a1613*k4+a1614*k10+a1615*k2)
                        call me%fcn(x+c16*h,y1,k3)
                        nfcn = nfcn + 3
                        ! final preparation
                        do j = 1 , nrd
                            i = me%icomp(j)
                            me%cont(j+nrd*4) = h*(me%cont(j+nrd*4)+d413*k4(i)+&
                                               d414*k10(i)+d415*k2(i)+d416*k3(i))
                            me%cont(j+nrd*5) = h*(me%cont(j+nrd*5)+d513*k4(i)+&
                                               d514*k10(i)+d515*k2(i)+d516*k3(i))
                            me%cont(j+nrd*6) = h*(me%cont(j+nrd*6)+d613*k4(i)+&
                                               d614*k10(i)+d615*k2(i)+d616*k3(i))
                            me%cont(j+nrd*7) = h*(me%cont(j+nrd*7)+d713*k4(i)+&
                                               d714*k10(i)+d715*k2(i)+d716*k3(i))
                        end do
                        me%hout = h
                    end if
                    k1 = k4
                    y = k5
                    me%xold = x
                    x = xph
                    if ( iout==1 .or. iout==2 .or. event ) then
                        call me%solout(naccpt+1,me%xold,x,y,irtrn,xout)
                        if ( irtrn<0 ) exit !abort
                    end if
                    ! normal exit
                    if ( last ) then
                        h = hnew
                        idid = 1
                        return
                    end if
                    if ( abs(hnew)>hmax ) hnew = posneg*hmax
                    if ( reject ) hnew = posneg*min(abs(hnew),abs(h))
                    reject = .false.
                else
                    ! step is rejected
                    hnew = h/min(facc1,fac11/safe)
                    reject = .true.
                    if ( naccpt>=1 ) nrejct = nrejct + 1
                    last = .false.
                end if
                h = hnew
            end if
        end do
    end if

    if ( iprint/=0 ) write (iprint,'(A,E18.4)') ' exit of dop853 at x=', x
    idid = 2

    end subroutine dp86co
!*****************************************************************************************

!*****************************************************************************************
!>
!  computation of an initial step size guess

    function hinit(me,x,y,posneg,f0,iord,hmax,atol,rtol,itol)
    use omp_lib
    implicit none

    class(dop853_class),intent(inout) :: me
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y       !! dimension(n)
    !dir$ assume_aligned y:64
    real(dp),intent(in)               :: posneg
    real(dp),dimension(:),intent(in)  :: f0      !! dimension(n)
    !dir$ assume_aligned f0:64
    integer(i4),intent(in)                :: iord
    real(dp),intent(in)               :: hmax
    real(dp),dimension(:),intent(in)  :: atol
    real(dp),dimension(:),intent(in)  :: rtol
    !dir$ assume_aligned atol:64
    !dir$ assume_aligned rtol:64
    integer(i4),intent(in)                :: itol

    real(dp) :: atoli,der12,der2,dnf,dny,h,h1,hinit,rtoli,sk
    integer(i4) :: i,n
    real(dp),dimension(me%n)  :: f1,y1

    ! compute a first guess for explicit euler as
    !   h = 0.01 * norm (y0) / norm (f0)
    ! the increment for explicit euler is small
    ! compared to the solution
    n=m%n
    dnf = 0.0_dp
    dny = 0.0_dp
    atoli = atol(1)
    rtoli = rtol(1)
    if ( itol==0 ) then
        !$omp simd reduction(+:dnf) reduction(+:dny) if(n>=16)
        do i = 1 , n
            sk = atoli + rtoli*abs(y(i))
            dnf = dnf + (f0(i)/sk)**2
            dny = dny + (y(i)/sk)**2
        end do
    else
         !$omp simd reduction(+:dnf) reduction(+:dny) if(n>=16)
        do i = 1 , n
            sk = atol(i) + rtol(i)*abs(y(i))
            dnf = dnf + (f0(i)/sk)**2
            dny = dny + (y(i)/sk)**2
        end do
    end if
    if ( dnf<=1.0e-10_dp .or. dny<=1.0e-10_dp ) then
        h = 1.0e-6_dp
    else
        h = sqrt(dny/dnf)*0.01_dp
    end if
    h = min(h,hmax)
    h = sign(h,posneg)
    ! perform an explicit euler step
    do i = 1 , me%n
        y1(i) = y(i) + h*f0(i)
    end do
    call me%fcn(x+h,y1,f1)
    ! estimate the second derivative of the solution
    der2 = 0.0_dp
    if ( itol==0 ) then
         !$omp simd reduction(+:der2) if(n>=16)
        do i = 1 , n
            sk = atoli + rtoli*abs(y(i))
            der2 = der2 + ((f1(i)-f0(i))/sk)**2
        end do
    else
        !$omp simd reduction(+:der2) if(n>=16)
        do i = 1 , n
            sk = atol(i) + rtol(i)*abs(y(i))
            der2 = der2 + ((f1(i)-f0(i))/sk)**2
        end do
    end if
    der2 = sqrt(der2)/h
    ! step size is computed such that
    !  h**iord * max ( norm (f0), norm (der2)) = 0.01
    der12 = max(abs(der2),sqrt(dnf))
    if ( der12<=1.0e-15_dp ) then
        h1 = max(1.0e-6_dp,abs(h)*1.0e-3_dp)
    else
        h1 = (0.01_dp/der12)**(1.0_dp/iord)
    end if

    h = min(100.0_dp*abs(h),h1,hmax)
    hinit = sign(h,posneg)

    end function hinit
!*****************************************************************************************

!*****************************************************************************************
!>
!  this function can be used for continuous output in connection
!  with the output-subroutine for [[dop853]]. it provides an
!  approximation to the `ii`-th component of the solution at `x`.

    function contd8(me,ii,x) result(y)

    implicit none

    class(dop853_class),intent(in) :: me
    integer(i4),intent(in)             :: ii
    real(dp),intent(in)            :: x
    real(dp)                       :: y

    real(dp) :: conpar, s, s1
    integer(i4) :: i,j,nd,ierr

    ! compute place of ii-th component
    i = 0
    do j = 1, me%nrdens
        if ( me%icomp(j)==ii ) i = j
    end do
    if ( i==0 ) then
        !always report this message, since it is an invalid use of the code.
        if (me%iprint==0) then
            ierr = error_unit
        else
            ierr = me%iprint
        end if
        write (ierr,*) &
            ' Error in contd8: no dense output available for component:', ii
        y = 0.0_dp
    else
        nd = me%nrdens
        s = (x-me%xold)/me%hout
        s1 = 1.0_dp - s
        conpar = me%cont(i+nd*4) + &
                 s*(me%cont(i+nd*5)+ &
                 s1*(me%cont(i+nd*6)+s*me%cont(i+nd*7)))
        y = me%cont(i) + &
            s*(me%cont(i+nd)+ &
            s1*(me%cont(i+nd*2)+&
            s*(me%cont(i+nd*3)+s1*conpar)))
    end if

    end function contd8
!*****************************************************************************************

!*****************************************************************************************
    end module dop853_module
!*****************************************************************************************
