
! KGEN-generated Fortran source file
!
! Filename    : shr_spfn_mod.F90
! Generated at: 2015-03-30 21:09:53
! KGEN version: 0.4.5



    MODULE shr_spfn_mod
        ! Module for common mathematical functions
        ! This #ifdef is to allow the module to be compiled with no dependencies,
        ! even on shr_kind_mod.
        USE shr_kind_mod, ONLY: r8 => shr_kind_r8
        IMPLICIT NONE
        PRIVATE
        ! Error functions



        ! Gamma functions
        ! Note that we lack an implementation of log_gamma, but we do have an
        ! implementation of the upper incomplete gamma function, which is not in
        ! Fortran 2008.
        ! Note also that this gamma function is only for double precision. We
        ! haven't needed an r4 version yet.
        PUBLIC shr_spfn_gamma

        INTERFACE shr_spfn_gamma
            MODULE PROCEDURE shr_spfn_gamma_r8
        END INTERFACE shr_spfn_gamma
        ! Mathematical constants
        ! sqrt(pi)
        ! Define machine-specific constants needed in this module.
        ! These were used by the original gamma and calerf functions to guarantee
        ! safety against overflow, and precision, on many different machines.
        ! By defining the constants in this way, we assume that 1/xmin is
        ! representable (i.e. does not overflow the real type). This assumption was
        ! not in the original code, but is valid for IEEE single and double
        ! precision.
        ! Double precision
        !---------------------------------------------------------------------
        ! Machine epsilon
        ! "Huge" value is returned when actual value would be infinite.
        ! Smallest normal value.
        ! Largest number that, when added to 1., yields 1.
        ! Largest argument for which erfcx > 0.
        ! Single precision
        !---------------------------------------------------------------------
        ! Machine epsilon
        ! "Huge" value is returned when actual value would be infinite.
        ! Smallest normal value.
        ! Largest number that, when added to 1., yields 1.
        ! Largest argument for which erfcx > 0.
        ! For gamma/igamma
        ! Approximate value of largest acceptable argument to gamma,
        ! for IEEE double-precision.
        CONTAINS

        ! write subroutines
        ! No subroutines
        ! No module extern variables
        ! Wrapper functions for erf


        ! Wrapper functions for erfc


        ! Wrapper functions for erfc_scaled



        elemental FUNCTION shr_spfn_gamma_r8(x) RESULT ( res )
            REAL(KIND=r8), intent(in) :: x
            REAL(KIND=r8) :: res
            ! Call intrinsic gamma.
            INTRINSIC gamma
            res = gamma(x)
        END FUNCTION shr_spfn_gamma_r8
        !------------------------------------------------------------------
        !
        ! 6 December 2006 -- B. Eaton
        ! The following comments are from the original version of CALERF.
        ! The only changes in implementing this module are that the function
        ! names previously used for the single precision versions have been
        ! adopted for the new generic interfaces.  To support these interfaces
        ! there is now both a single precision version (calerf_r4) and a
        ! double precision version (calerf_r8) of CALERF below.  These versions
        ! are hardcoded to use IEEE arithmetic.
        !
        !------------------------------------------------------------------
        !
        ! This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
        !   for a real argument  x.  It contains three FUNCTION type
        !   subprograms: ERF, ERFC, and ERFCX (or ERF_R8, ERFC_R8, and ERFCX_R8),
        !   and one SUBROUTINE type subprogram, CALERF.  The calling
        !   statements for the primary entries are:
        !
        !                   Y=ERF(X)     (or   Y=ERF_R8(X)),
        !
        !                   Y=ERFC(X)    (or   Y=ERFC_R8(X)),
        !   and
        !                   Y=ERFCX(X)   (or   Y=ERFCX_R8(X)).
        !
        !   The routine  CALERF  is intended for internal packet use only,
        !   all computations within the packet being concentrated in this
        !   routine.  The function subprograms invoke  CALERF  with the
        !   statement
        !
        !          CALL CALERF(ARG,RESULT,JINT)
        !
        !   where the parameter usage is as follows
        !
        !      Function                     Parameters for CALERF
        !       call              ARG                  Result          JINT
        !
        !     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
        !     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
        !     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
        !
        !   The main computation evaluates near-minimax approximations
        !   from "Rational Chebyshev approximations for the error function"
        !   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
        !   transportable program uses rational functions that theoretically
        !   approximate  erf(x)  and  erfc(x)  to at least 18 significant
        !   decimal digits.  The accuracy achieved depends on the arithmetic
        !   system, the compiler, the intrinsic functions, and proper
        !   selection of the machine-dependent constants.
        !
        !*******************************************************************
        !*******************************************************************
        !
        ! Explanation of machine-dependent constants
        !
        !   XMIN   = the smallest positive floating-point number.
        !   XINF   = the largest positive finite floating-point number.
        !   XNEG   = the largest negative argument acceptable to ERFCX;
        !            the negative of the solution to the equation
        !            2*exp(x*x) = XINF.
        !   XSMALL = argument below which erf(x) may be represented by
        !            2*x/sqrt(pi)  and above which  x*x  will not underflow.
        !            A conservative value is the largest machine number X
        !            such that   1.0 + X = 1.0   to machine precision.
        !   XBIG   = largest argument acceptable to ERFC;  solution to
        !            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
        !            W(x) = exp(-x*x)/[x*sqrt(pi)].
        !   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
        !            machine precision.  A conservative value is
        !            1/[2*sqrt(XSMALL)]
        !   XMAX   = largest acceptable argument to ERFCX; the minimum
        !            of XINF and 1/[sqrt(pi)*XMIN].
        !
        !   Approximate values for some important machines are:
        !
        !                          XMIN       XINF        XNEG     XSMALL
        !
        !  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
        !  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
        !  IEEE (IBM/XT,
        !    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
        !  IEEE (IBM/XT,
        !    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
        !  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
        !  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
        !  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
        !  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
        !
        !
        !                          XBIG       XHUGE       XMAX
        !
        !  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
        !  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
        !  IEEE (IBM/XT,
        !    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
        !  IEEE (IBM/XT,
        !    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
        !  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
        !  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
        !  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
        !  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
        !
        !*******************************************************************
        !*******************************************************************
        !
        ! Error returns
        !
        !  The program returns  ERFC = 0      for  ARG .GE. XBIG;
        !
        !                       ERFCX = XINF  for  ARG .LT. XNEG;
        !      and
        !                       ERFCX = 0     for  ARG .GE. XMAX.
        !
        !
        ! Intrinsic functions required are:
        !
        !     ABS, AINT, EXP
        !
        !
        !  Author: W. J. Cody
        !          Mathematics and Computer Science Division
        !          Argonne National Laboratory
        !          Argonne, IL 60439
        !
        !  Latest modification: March 19, 1990
        !
        !------------------------------------------------------------------

        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !! Incomplete Gamma function
        !!
        !! @author  Tianyi Fan
        !! @version August-2010

    END MODULE shr_spfn_mod
