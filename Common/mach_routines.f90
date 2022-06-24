!*****************************************************************************************
!>
!  Fortran machine constants.
!
!### Reference
!  * P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!    a portable library, ACM Transactions on Mathematical
!    Software 4, 2 (June 1978), pp. 177-188.
!
!### Author
!  * Jacob Williams, June 2022

   module mach_routines

   !use iso_fortran_env
   use mod_kinds, only : i4,sp,dp
   implicit none

   private

   !integer, parameter :: sp = kind(1.0)    !! single precision kind
   !i!nteger, parameter :: dp = kind(1.0d0)  !! double precision kind

   public :: r1mach, d1mach, i1mach

   contains
!*****************************************************************************************

!*******************************************************************
!>
!  Return floating point machine dependent constants.
!
!  D1MACH can be used to obtain machine-dependent parameters for the
!  local machine environment.  It is a function subprogram with one
!  (input) argument, and can be referenced as follows:
!
!```fortran
!       A = D1MACH(I)
!```
!
!  where I=1,...,5.  The (output) value of A above is determined by
!  the (input) value of I.  The results for various values of I are
!  discussed below.
!
!   * `D1MACH(1) = B**(EMIN-1)`, the smallest positive magnitude.
!   * `D1MACH(2) = B**EMAX*(1 - B**(-T))`, the largest magnitude.
!   * `D1MACH(3) = B**(-T)`, the smallest relative spacing.
!   * `D1MACH(4) = B**(1-T)`, the largest relative spacing.
!   * `D1MACH(5) = LOG10(B)`
!
!  Assume single precision numbers are represented in the T-digit,
!  base-B form
!
!             sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!  where 0 <= X(I) < B for I=1,...,T, 0 < X(1), and
!  EMIN <= E <= EMAX.
!
!  The values of B, T, EMIN and EMAX are provided in I1MACH as
!  follows:
!
!  * `I1MACH(10) = B`, the base.
!  * `I1MACH(11) = T`, the number of base-B digits.
!  * `I1MACH(12) = EMIN`, the smallest exponent E.
!  * `I1MACH(13) = EMAX`, the largest exponent E.
!
!### Author
!  * Fox, P. A., (Bell Labs)
!  * Hall, A. D., (Bell Labs)
!  * Schryer, N. L., (Bell Labs)
!
!### Revision history
!  * 790101  DATE WRITTEN
!  * 960329  Modified for Fortran 90 (BE after suggestions by EHG)
!  * 220619  Modified for Fortran 2008 (JW)

   pure real(dp) function d1mach(i)

      integer(i4), intent(in) :: i

      real(dp), parameter :: x = 1.0_dp
      real(dp), parameter :: b = real(radix(x), dp)

      select case (i)
      case (1); d1mach = b**(minexponent(x) - 1) ! the smallest positive magnitude.
      case (2); d1mach = huge(x)                 ! the largest magnitude.
      case (3); d1mach = b**(-digits(x))         ! the smallest relative spacing.
      case (4); d1mach = b**(1 - digits(x))      ! the largest relative spacing.
      case (5); d1mach = log10(b)
      case default
         error stop 'Error in d1mach - i out of bounds'
      end select

   end function d1mach
!*******************************************************************

!*******************************************************************
!>
!  Return integer machine dependent constants.
!
!   I1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument and can be referenced as follows:
!
!```fortran
!        K = I1MACH(I)
!```
!
!   where I=1,...,16.  The (output) value of K above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   I/O unit numbers:
!     * `I1MACH( 1)`   = the standard input unit.
!     * `I1MACH( 2)`   = the standard output unit.
!     * `I1MACH( 3)`   = the standard punch unit.
!     * `I1MACH( 4)`   = the standard error message unit.
!
!   Words:
!     * `I1MACH( 5)` = the number of bits per integer storage unit.
!     * `I1MACH( 6)` = the number of characters per integer storage unit.
!
!   Integers:
!     assume integers are represented in the S-digit, base-A form
!```
!                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!                where 0 <= X(I) < A for I=0,...,S-1.
!```
!
!     * `I1MACH( 7) = A`, the base.
!     * `I1MACH( 8) = S`, the number of base-A digits.
!     * `I1MACH( 9) = A**S - 1`, the largest magnitude.
!
!   Floating-Point Numbers:
!     Assume floating-point numbers are represented in the T-digit, base-B form
!```
!                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!                where 0 <= X(I) < B for I=1,...,T,
!                0 < X(1), and EMIN <= E <= EMAX.
!```
!     `I1MACH(10) = B`, the base.
!
!   Single-Precision:
!     * `I1MACH(11) = T`, the number of base-B digits.
!     * `I1MACH(12) = EMIN`, the smallest exponent E.
!     * `I1MACH(13) = EMAX`, the largest exponent E.
!
!   Double-Precision:
!     * `I1MACH(14) = T`, the number of base-B digits.
!     * `I1MACH(15) = EMIN`, the smallest exponent E.
!     * `I1MACH(16) = EMAX`, the largest exponent E.
!
!### Author
!  * Fox, P. A., (Bell Labs)
!  * Hall, A. D., (Bell Labs)
!  * Schryer, N. L., (Bell Labs)
!
!### Revision history
!  * 750101  DATE WRITTEN
!  * 960411  Modified for Fortran 90 (BE after suggestions by EHG).
!  * 980727  Modified value of I1MACH(6) (BE after suggestion by EHG).
!  * 220619  Modified for Fortran 2008 (JW)

   pure integer(i4) function i1mach(i)

      integer(i4), intent(in) :: i

      real(sp), parameter :: x = 1.0_sp
      real(dp), parameter :: xx = 1.0_dp

      select case (i)
      case (1);  i1mach = input_unit
      case (2);  i1mach = output_unit
      case (3);  i1mach = 0              ! Punch unit is no longer used
      case (4);  i1mach = error_unit
      case (5);  i1mach = numeric_storage_size
      case (6);  i1mach = numeric_storage_size / character_storage_size ! characters per integer
      case (7);  i1mach = radix(1)
      case (8);  i1mach = numeric_storage_size - 1
      case (9);  i1mach = huge(1)
      case (10); i1mach = radix(x)
      case (11); i1mach = digits(x)
      case (12); i1mach = minexponent(x)
      case (13); i1mach = maxexponent(x)
      case (14); i1mach = digits(xx)
      case (15); i1mach = minexponent(xx)
      case (16); i1mach = maxexponent(xx)
      case default
         error stop 'Error in i1mach - i out of bounds'
      end select

   end function i1mach
!*******************************************************************

!*******************************************************************
!>
!  Return floating point machine dependent constants.
!!
!   R1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!```fortran
!        A = R1MACH(I)
!```
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   * `R1MACH(1) = B**(EMIN-1)`, the smallest positive magnitude.
!   * `R1MACH(2) = B**EMAX*(1 - B**(-T))`, the largest magnitude.
!   * `R1MACH(3) = B**(-T)`, the smallest relative spacing.
!   * `R1MACH(4) = B**(1-T)`, the largest relative spacing.
!   * `R1MACH(5) = LOG10(B)`
!
!   Assume single precision numbers are represented in the T-digit, base-B form
!```
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 <= X(I) < B for I=1,...,T, 0 < X(1), and
!   EMIN <= E <= EMAX.
!```
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   * `I1MACH(10) = B`, the base.
!   * `I1MACH(11) = T`, the number of base-B digits.
!   * `I1MACH(12) = EMIN`, the smallest exponent E.
!   * `I1MACH(13) = EMAX`, the largest exponent E.
!
!### Author
!  * Fox, P. A., (Bell Labs)
!  * Hall, A. D., (Bell Labs)
!  * Schryer, N. L., (Bell Labs)
!
!### Revision history
!  * 790101  DATE WRITTEN
!  * 960329  Modified for Fortran 90 (BE after suggestions by EG)
!  * 220619  Modified for Fortran 2008 (JW)

   pure real(sp) function r1mach(i)

      integer(i4), intent(in) :: i

      real(sp), parameter :: x = 1.0_sp
      real(sp), parameter :: b = real(radix(x), sp)

      select case (i)
      case (1); r1mach = b**(minexponent(x) - 1) ! the smallest positive magnitude.
      case (2); r1mach = huge(x)                 ! the largest magnitude.
      case (3); r1mach = b**(-digits(x))         ! the smallest relative spacing.
      case (4); r1mach = b**(1 - digits(x))      ! the largest relative spacing.
      case (5); r1mach = log10(b)
      case default
         error stop 'Error in r1mach - i out of bounds'
      end select

   end function r1mach
!*******************************************************************

!*****************************************************************************************
   end module mach_routines
!*****************************************************************************************
