module Logical_Constants
! ##############################################################################
!
! Copyright 2018 IRD
!
! This file is part of statpack.
!
! statpack is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of 
! the License, or (at your option) any later version.
!
! statpack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You can find a copy of the GNU Lesser General Public License
! in the statpack/doc directory.
!
! ##############################################################################
!                                                                              *
! ******************************************************************************
! MODULE EXPORTING LOGICAL CONSTANTS OF KIND 'lgl'.                            *
!                                                                              *
! BY ONLY USING LOGICAL VALUES AS DEFINED WITHIN THIS MODULE (e.g. THE         *
! LOGICAL CONSTANTS true AND false OF KIND lgl), ALL PROBLEMS ASSOCIATED WITH  *
! THE CONVERSION OF LOGICAL LITERAL VALUES CAN BE TOTALLY AVOIDED.             *
!                                                                              *
! LATEST REVISION : 21/03/2018                                                 *
!                                                                              *
! ##############################################################################
!                                                                              *
! ******************************************************************************
!
!
! USED MODULES
! ============
!
    use Select_Parameters, only : lgl
!
! STRONG TYPING IMPOSED 
! =====================
!    
    implicit none
!
! PUBLIC ENTITIES 
! ===============
!
    private
    public :: true, false
!
! MODULE PARAMETERS 
! =================
!
! SOME CONVENIENT TRUE FALSE ABBREVIATIONS.
!
    logical(lgl), parameter  :: true  = .true._lgl
    logical(lgl), parameter  :: false = .false._lgl
!    
!    
! *********************************
! END OF MODULE Logical_Constants *
! *********************************
!    
end module Logical_Constants
!    
! ***********************************************************************************************
!
module Reals_Constants   
! ##############################################################################
!
! Copyright 2018 IRD
!
! This file is part of statpack.
!
! statpack is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of 
! the License, or (at your option) any later version.
!
! statpack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You can find a copy of the GNU Lesser General Public License
! in the statpack/doc directory.
!
! ##############################################################################
!                                                                              *
! ******************************************************************************
! THIS MODULE PROVIDES NAMES FOR ALL REQUIRED LITERAL REAL VALUES OF KIND      *
! 'stnd' AND 'extd' USED IN STATPACK.                                          *
!                                                                              *
! BY ONLY USING REAL VALUES AS DEFINED WITHIN THIS MODULE,                     *
! ALL PROBLEMS ASSOCIATED WITH THE PRECISION OF REAL                           *
! LITERAL VALUES CAN BE TOTALLY AVOIDED.                                       *
!                                                                              *
! LATEST REVISION : 31/10/2018                                                 *
!                                                                              *
! ##############################################################################
!                                                                              *
! ******************************************************************************
!
!
! USED MODULES 
! ============
!    
    use Select_Parameters, only : stnd, extd
!
! STRONG TYPING IMPOSED 
! =====================
!    
    implicit none
!
! PUBLIC ENTITIES 
! ===============
!    
    public
    private :: stnd, extd
!
! MODULE PARAMETERS 
! =================
!
! REAL CONSTANTS AT PRECISION stnd.
!
!   DIGITS
!
    real(stnd), parameter ::      &
        zero    = 0,              & 
        one     = 1,              & 
        two     = 2,              & 
        three   = 3,              &
        four    = 4,              &
        five    = 5,              &
        six     = 6,              &
        seven   = 7,              &
        eight   = 8,              &
        nine    = 9,              &
        ten     =10               
!
!   TENTHS
!
    real(stnd), parameter ::      &
        c0_1    = 0.1_stnd,       &
        c0_2    = 0.2_stnd,       &
        c0_3    = 0.3_stnd,       &
        c0_4    = 0.4_stnd,       &
        c0_5    = 0.5_stnd,       &
        c0_6    = 0.6_stnd,       &
        c0_7    = 0.7_stnd,       &
        c0_8    = 0.8_stnd,       &
        c0_9    = 0.9_stnd
!
!   RECIPROCALS
!
    real(stnd), parameter ::      &
        tenth   = 0.1_stnd,       &
        ninth   = one/nine,       &
        eighth  = 0.125_stnd,     &
        seventh = one/seven,      &
        sixth   = one/six,        &
        fifth   = 0.2_stnd,       &
        quarter = 0.25_stnd,      &
        third   = one/three,      &
        half    = 0.5_stnd
!
!   FRACTIONS a/b NAMED AS fa_b
!
    real(stnd), parameter ::      &
        f1_29   = one/29,         &
        f2_3    = two/three,      &
        f4_3    = four/three,     &
        f7_3    = seven/three
!
!   INTEGRAL VALUES TO 99
!
    real(stnd), parameter ::                                            &
        c10     = 10,           c40     = 40,           c70     = 70,   &
        c11     = 11,           c41     = 41,           c71     = 71,   &
        c12     = 12,           c42     = 42,           c72     = 72,   &
        c13     = 13,           c43     = 43,           c73     = 73,   &
        c14     = 14,           c44     = 44,           c74     = 74,   &
        c15     = 15,           c45     = 45,           c75     = 75,   &
        c16     = 16,           c46     = 46,           c76     = 76,   &
        c17     = 17,           c47     = 47,           c77     = 77,   &
        c18     = 18,           c48     = 48,           c78     = 78,   &
        c19     = 19,           c49     = 49,           c79     = 79,   &
        c20     = 20,           c50     = 50,           c80     = 80,   &
        c21     = 21,           c51     = 51,           c81     = 81,   &
        c22     = 22,           c52     = 52,           c82     = 82,   &
        c23     = 23,           c53     = 53,           c83     = 83,   &
        c24     = 24,           c54     = 54,           c84     = 84,   &
        c25     = 25,           c55     = 55,           c85     = 85,   &
        c26     = 26,           c56     = 56,           c86     = 86,   &
        c27     = 27,           c57     = 57,           c87     = 87,   &
        c28     = 28,           c58     = 58,           c88     = 88,   &
        c29     = 29,           c59     = 59,           c89     = 89,   &
        c30     = 30,           c60     = 60,           c90     = 90,   &
        c31     = 31,           c61     = 61,           c91     = 91,   &
        c32     = 32,           c62     = 62,           c92     = 92,   &
        c33     = 33,           c63     = 63,           c93     = 93,   &
        c34     = 34,           c64     = 64,           c94     = 94,   &
        c35     = 35,           c65     = 65,           c95     = 95,   &
        c36     = 36,           c66     = 66,           c96     = 96,   &
        c37     = 37,           c67     = 67,           c97     = 97,   &
        c38     = 38,           c68     = 68,           c98     = 98,   &
        c39     = 39,           c69     = 69,           c99     = 99
!
!   MISCELLANEOUS INTEGRAL VALUES
!
    real(stnd), parameter ::      &
        c100    = 100,            &
        c120    = 120,            &
        c180    = 180,            &
        c200    = 200,            &
        c256    = 256,            &
        c300    = 300,            &
        c360    = 360,            &
        c400    = 400,            &
        c500    = 500,            &
        c600    = 600,            &
        c681    = 681,            &
        c700    = 700,            &
        c800    = 800,            &
        c900    = 900,            &
        c991    = 991,            &
        c1000   = 1000,           &
        c1162   = 1162,           &
        c2324   = 2324,           &
        c2000   = 2000,           & 
        c3000   = 3000,           & 
        c4000   = 4000,           & 
        c5000   = 5000,           & 
        c10000  = 10000,          &
        c20700  = 20700,          &
        c40000  = 40000
!
!    MISCELLANEOUS REAL VALUES
!
!    FORM: d.dd          NAMED AS  cd_dd
!                   nn
!    FORM: d.dd x 10     NAMED AS  cd_ddenn
!                   -nn
!    FORM: d.dd x 10     NAMED AS  cd_ddmnn
!
    real(stnd), parameter ::      &
        c1_e6   = 1.0e6_stnd,     &
        c1_m1   = 1.0e-1_stnd,    &
        c1_m2   = 1.0e-2_stnd,    &
        c1_m3   = 1.0e-3_stnd,    & 
        c1_m4   = 1.0e-4_stnd,    &
        c1_m5   = 1.0e-5_stnd,    &
        c1_m6   = 1.0e-6_stnd,    &
        c1_m7   = 1.0e-7_stnd,    &
        c1_m8   = 1.0e-8_stnd,    &
        c1_m9   = 1.0e-9_stnd,    &
        c1_m10  = 1.0e-10_stnd,   &
        c1_m11  = 1.0e-11_stnd,   &
        c1_m12  = 1.0e-12_stnd,   &
        c1_m13  = 1.0e-13_stnd,   &
        c1_m14  = 1.0e-14_stnd,   &
        c1_m15  = 1.0e-15_stnd,   &
        c1_m16  = 1.0e-16_stnd,   &
        c2_m3   = 2.0e-3_stnd,    &
        c2_m6   = 2.0e-6_stnd,    &
        c4_m2   = 4.0e-2_stnd,    &
        c5_m2   = 5.0e-2_stnd,    &
        c0_08   = 0.08_stnd,      &
        c0_089  = 0.089_stnd,     &
        c0_1136 = 0.1136_stnd,    &
        c0_717  = 0.717_stnd,     &
        c0_822  = 0.822_stnd,     &
        c0_999  = 0.999_stnd,     &
        c1_0001 = 1.0001_stnd,    &
        c1_01   = 1.01_stnd,      &
        c1_2    = 1.2_stnd,       &
        c1_3    = 1.3_stnd,       &
        c1_5    = 1.5_stnd,       &
        c1_6    = 1.6_stnd,       &
        c2_25   = 2.25_stnd,      &
        c2_5    = 2.5_stnd,       &
        c2_625  = 2.625_stnd,     &
        c3_3    = 3.3_stnd,       &
        c4_5    = 4.5_stnd,       &
        c6_3    = 6.3_stnd,       &
        c7_5    = 7.5_stnd,       &
        c10_1   = 10.1_stnd,      &
        c18_8055 = 18.8055_stnd,  &
        c19_8   = 19.8_stnd,      &
        c20_2   = 20.2_stnd,      &
        c85_5   = 85.5_stnd,      &
        c94_5   = 94.5_stnd,      &
        c96_36  = 96.36_stnd
!
!   FREQUENTLY USED MATHEMATICAL CONSTANTS
!
    real(stnd), parameter ::                                                &
        pi       = 3.1415926535897932384626433832795028841971693993751_stnd,&
        pio2     = 1.57079632679489661923132169163975144209858_stnd,        &
        twopi    = 6.283185307179586476925286766559005768394_stnd,          &  
        sqrt2    = 1.41421356237309504880168872420969807856967_stnd,        &
        euler    = 0.5772156649015328606065120900824024310422_stnd
!
!
! REAL CONSTANTS AT PRECISION extd.
!
!   DIGITS IN EXTENDED PRECISION
!
    real(extd), parameter ::           &
        zero_extd    = 0,              & 
        one_extd     = 1,              & 
        two_extd     = 2,              & 
        three_extd   = 3,              &
        four_extd    = 4,              &
        five_extd    = 5,              &
        six_extd     = 6,              &
        seven_extd   = 7,              &
        eight_extd   = 8,              &
        nine_extd    = 9,              &
        ten_extd     =10               
!
!   RECIPROCALS IN EXTENDED PRECISION
!
    real(extd), parameter ::                     &
        tenth_extd   = 0.1_extd,                 &
        ninth_extd   = one_extd/nine_extd,       &
        eighth_extd  = 0.125_extd,               &
        seventh_extd = one_extd/seven_extd,      &
        sixth_extd   = one_extd/six_extd,        &
        fifth_extd   = 0.2_extd,                 &
        quarter_extd = 0.25_extd,                &
        third_extd   = one_extd/three_extd,      &
        half_extd    = 0.5_extd
!
!   MISCELLANEOUS INTEGRAL VALUES IN EXTENDED PRECISION
!
    real(extd), parameter ::           &
        c37_extd     = 37,             &
        c1000_extd   = 1000
!
!    MISCELLANEOUS REAL VALUES IN EXTENDED PRECISION
!
!    FORM: d.dd          NAMED AS  cd_dd_extd
!                   nn
!    FORM: d.dd x 10     NAMED AS  cd_ddenn_extd
!                   -nn
!    FORM: d.dd x 10     NAMED AS  cd_ddmnn_extd
!
    real(extd), parameter ::           &
        c1_m15_extd    = 1.0e-15_extd, &
        c1_m16_extd    = 1.0e-16_extd, &
        c0_65_extd     = 0.65_extd,    &
        c6_5_extd      = 6.5_extd,     &
        c365_25_extd   = 365.25_extd,  &
        c0_0004_extd   = .0004_extd
!
!   FREQUENTLY USED MATHEMATICAL CONSTANTS IN EXTENDED PRECISION
!
    real(extd), parameter ::                                                      &
        pi_extd        = 3.1415926535897932384626433832795028841971693993751_extd,&
        pio2_extd      = 1.57079632679489661923132169163975144209858_extd,        &
        twopi_extd     = 6.283185307179586476925286766559005768394_extd,          &  
        sqrt2_extd     = 1.41421356237309504880168872420969807856967_extd,        &
        euler_extd     = 0.5772156649015328606065120900824024310422_extd,         &
        lnsqrt2pi_extd = 0.9189385332046727_extd
!
!    
!    
! ******************************* 
! END OF MODULE Reals_Constants *
! *******************************
!    
end module Reals_Constants
!    
! ***********************************************************************************************
!
module Num_Constants   
! ##############################################################################
!
! Copyright 2018 IRD
!
! This file is part of statpack.
!
! statpack is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of 
! the License, or (at your option) any later version.
!
! statpack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You can find a copy of the GNU Lesser General Public License
! in the statpack/doc directory.
!
! ##############################################################################
!                                                                              *
! ******************************************************************************
! THIS MODULE PROVIDES SIMPLE NAMES AND ROUTINES FOR THE VARIOUS               *
! MACHINE DEPENDENT CONSTANTS. ALL ARE FOR PRECISION 'stnd'.                   *
!                                                                              *
! LATEST REVISION : 28/09/2018                                                 *
!                                                                              *
! ##############################################################################
!                                                                              *
! ******************************************************************************
!                                                                              *
! THE C PROCESSOR MACROS USED IN THIS MODULE ARE:                              *
!                                                                              *
!  _F2003       FOR ACTIVATING FORTRAN 2003 CONSTRUCTS                         *
!                                                                              *
! ******************************************************************************
!
!
! USED MODULES 
! ============
!    
    use Select_Parameters,     only : stnd, i4b, lgl
    use Reals_Constants,       only : zero, one, two, c2000
    use Logical_Constants
!
! USED MODULES IN SPECIFIC PROCEDURES
! ===================================
!   
!   use ieee_arithmetic  (if the CPP macro _F2003 is activated)
!
! STRONG TYPING IMPOSED 
! =====================
!    
    implicit none
!
! PUBLIC ENTITIES 
! ===============
!    
    private
    public ::                   &
       maxexp,                  &
       minexp,                  &
       base,                    &
       nbasedigits,             &
       decprec,                 &
       decexpr,                 &
       machmaxexp,              &
       machminexp,              &
       machbase,                &
       machnbasedigits,         &
       machdecprec,             &
       machdecexpr,             &
       macheps,                 & 
       machulp,                 & 
       machtiny,                &
       machhuge,                &
       machsmlnum,              &
       machbignum,              &
       lamch,                   &
       mach,                    &
       test_ieee,               &
       test_nan,                &
       is_nan,                  &
       replace_nan,             &
       true_nan,                &
       nan
!
! GENERIC INTERFACES FOR ROUTINES WITH OVERLOADED VERSIONS
! ========================================================
!
    interface is_nan
        module procedure    is_nan_r, is_nan_rv, is_nan_rm
    end interface
!
    interface replace_nan
        module procedure    replace_nan_r, replace_nan_rv, replace_nan_rm
    end interface
!
! MODULE VARIABLES
! ================
!    
! MACHINE NUMERIC CHARACTERISTICS WITHOUT INTRINSIC VALUES 
! IN FORTRAN90. PRIVATE ENTITIES.
!
    real(stnd), save ::   &
        unitrnd             ! UNIT ROUNDOFF OF THE MACHINE
!
    integer(i4b), save :: &
        nrnd,             & ! ROUNDING INDICATOR FOR FLOATING-POINT ADDITION
        ngrd                ! GUARD DIGITS INDICATOR FOR CHOPPING ARITHMETIC
!
!
! logical INDICATOR FOR ROUTINES lamch, mach AND test_ieee . PRIVATE ENTITY.
!
    logical(lgl) :: first_mach = true
!
! MODULE PARAMETERS 
! =================
!
! VALUES OF MACHINE NUMERIC CHARACTERISTICS
! WITH INTRINSIC VALUES IN FORTRAN90.
!
    integer(i4b), parameter ::             &
       maxexp      = maxexponent(unitrnd), & ! LARGEST EXPONENT BEFORE OVERFLOW
       minexp      = minexponent(unitrnd), & ! MINIMUM EXPONENT BEFORE (GRADUAL) UNDERFLOW
       base        = radix(unitrnd),       & ! BASE OF THE MACHINE
       nbasedigits = digits(unitrnd),      & ! NUMBER OF (base) DIGITS IN THE MANTISSA
       decprec     = precision(unitrnd),   & ! NUMBER OF EQUIVALENT DECIMAL DIGITS IN THE MANTISSA
       decexpr     = range(unitrnd)          ! EQUIVALENT DECIMAL EXPONENT RANGE
!
    real(stnd), parameter ::                &
       machmaxexp      = maxexp,            & ! LARGEST EXPONENT BEFORE OVERFLOW
       machminexp      = minexp,            & ! MINIMUM EXPONENT BEFORE (GRADUAL) UNDERFLOW
       machbase        = base,              & ! BASE OF THE MACHINE
       machnbasedigits = nbasedigits,       & ! NUMBER OF (base) DIGITS IN THE MANTISSA
       machdecprec     = decprec,           & ! NUMBER OF EQUIVALENT DECIMAL DIGITS IN THE MANTISSA
       machdecexpr     = decexpr,           & ! EQUIVALENT DECIMAL EXPONENT RANGE
       macheps         = epsilon(unitrnd),  & ! MACHINE EPSILON     :  base**(1-nbasedigits)
       machulp         = macheps*machbase,  & ! MACHINE PRECISION   :  base*macheps
       machtiny        = tiny(unitrnd),     & ! UNDERFLOW THRESHOLD :  base**(minexp-1)
       machhuge        = huge(unitrnd),     & ! OVERFLOW THRESHOLD  : (base**maxexp)*(1-base**(-nbasedigits))
       machsmlnum      = machtiny/machulp,  & ! SCALED MINIMUM
       machbignum      = one/machsmlnum       ! SCALED MAXIMUM
!
!
! =============================================================================================
!
                                  contains
!                                 ========
!
!
! =============================================================================================
!
!
    function lamch( cmach )
!
! Purpose
! _______
!
!   LAMCH determines machine parameters for precision STND.
!
!
! Arguments
! _________
!
!   CMACH   (INPUT) character*1
!           Specifies the value to be returned by LAMCH. If:
!
!           - CMACH = 'S' or 's',   LAMCH := sfmin
!           - CMACH = 'T' or 't',   LAMCH := t
!           - CMACH = 'R' or 'r',   LAMCH := rnd
!           - CMACH = 'G' or 'g',   LAMCH := grd
!           - CMACH = 'U' or 'u',   LAMCH := unitrnd
!           - CMACH = 'P' or 'p',   LAMCH := prec
!
!           where:
!
!           - sfmin = safe minimum, such that 1/sfmin does not overflow.
!
!           - t = number of (base) digits in the floating-point significand.
!
!           - rnd = 0.0 when floating-point addition rounds upward, downward 
!                   or toward zero;
!
!                 = 1.0 when floating-point addition rounds to nearest,
!                   but not in the IEEE style;
!
!                 = 2.0 when floating-point addition rounds in the IEEE style.
!
!           - grd = 1. if floating-point arithmetic chops (rnd = 0.) and more 
!                   than t digits participate in the post-normalization shift 
!                   of the floating-point significand in multiplication,
!
!                 = 0.0 otherwise.
!
!           - unitrnd = unit roundoff of the machine.
!
!           - prec = unitrnd*machbase.
!
!
! Further Details
! _______________
!
!   The routine is based on the routine  DLAMCH  in LAPACK77 (version 3).
!
!   For any other characters, LAMCH returns the bit pattern corresponding to a quiet NaN.
!
!
! _____________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    character, intent(in) :: cmach
!
    real(stnd) :: lamch
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd), save :: rnd, grd, prec, sfmin=machtiny
    real(stnd)       :: a
!
    logical(lgl) :: first=true
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!$OMP CRITICAL
!
    if ( first ) then
!
        first = false
!
!       comp_mach RETURNS MACHINE-SPECIFIC PARAMETERS IN ITS ARGUMENT LIST.
!
        if ( first_mach ) call comp_mach( nrnd, ngrd, unitrnd )
!
        rnd = real( nrnd, stnd )
        grd = real( ngrd, stnd )
!
!       DETERMINE PRECISION.
!
        prec = machbase*unitrnd
!
!       FIND sfmin.
!
        a = one/machhuge
!
!       IF a>=machtiny, USE a PLUS A BIT, TO AVOID THE POSSIBILITY
!       OF ROUNDING CAUSING OVERFLOW WHEN COMPUTING 1/sfmin.
!
        if ( a>=machtiny ) sfmin = nearest( a, one )
!
    end if
!
!$OMP END CRITICAL
!
!   RETURN REQUIRED FUNCTION VALUE.
!
    select case( cmach )
!
        case( achar(83), achar(115) ) ; lamch = sfmin 
        case( achar(84), achar(116) ) ; lamch = machnbasedigits
        case( achar(82), achar(114) ) ; lamch = rnd 
        case( achar(71), achar(103) ) ; lamch = grd 
        case( achar(85), achar(117) ) ; lamch = unitrnd 
        case( achar(80), achar(112) ) ; lamch = prec
        case default                  ; lamch = true_nan()
!
    end select
!
!
! END OF FUNCTION lamch
! _____________________
!
    end function lamch
!
! =============================================================================================
!
    subroutine mach( basedigits, irnd, iuflow, igrd, iexp, ifloat, expepspos, expepsneg,    &
                     minexpbase, maxexpbase, epspos, epsneg, epsilpos, epsilneg, rndunit )
!
! Purpose
! _______
!
!   This subroutine is intended to determine the parameters of the floating-point
!   arithmetic system specified below.
!                                                                              
!
! Arguments
! _________
!
!   BASEDIGITS  (OUTPUT, OPTIONAL) integer(i4b)         
!               The number of base digits in the floating-point significand.
!
!   IRND        (OUTPUT, OPTIONAL) integer(i4b)
!               A parameter indicating whether proper rounding or chopping 
!               (rounding upward, downward, toward zero) occurs in addition. If:
!
!               - IRND = 0 if floating-point addition rounds upward, downward
!                 or toward zero;
!               - IRND = 1 if floating-point addition rounds to nearest,
!                 but not in the IEEE style;
!               - IRND = 2 if floating-point addition rounds in the IEEE style.
!
!   IUFLOW      (OUTPUT, OPTIONAL) integer(i4b)
!               A parameter indicating whether underflow is full or partial:
!
!               - IUFLOW = 0   if there is full underflow (flush to zero, etc);
!               - IUFLOW = 1   if there is partial underflow.
!
!   IGRD        (OUTPUT, OPTIONAL) integer(i4b)
!               The number of guard digits for multiplication with chopping
!               arithmetic (IRND = 0). If:
!
!               - IGRD = 0   if floating-point arithmetic rounds, or if it
!                 chops and only  BASEDIGITS  digits participate
!                 in the post-normalization shift of the
!                 floating-point significand in multiplication;
!               - IGRD = 1   if floating-point arithmetic chops and more
!                 than  BASEDIGITS  digits participate in the
!                 post-normalization shift of the floating-point
!                 significand in multiplication.
!
!   IEXP        (OUTPUT, OPTIONAL) integer(i4b)
!               A guess for the number of bits dedicated to the representation            
!               of the exponent of a floating point number if BASE is
!               a power of two and -1 otherwise.
!
!   IFLOAT      (OUTPUT, OPTIONAL) integer(i4b)
!               A guess for the number of bits dedicated to the representation            
!               of a floating point number if BASE is a power of two 
!               and -1 otherwise.                    
!                                                                              
!   EXPEPSPOS   (OUTPUT, OPTIONAL) integer(i4b)
!               The largest in magnitude negative integer such that
!
!                   1.0 + float(base)**(expepspos) /= 1.         
!
!   EXPEPSNEG   (OUTPUT, OPTIONAL) integer(i4b)         
!               The largest in magnitude negative integer such that            
!
!                   1.0 - float(base)**(expepsneg) /= 1.         
!
!   MINEXPBASE  (OUTPUT, OPTIONAL) integer(i4b)
!               The largest in magnitude negative integer such that         
!               float(base)**minexpbase is positive and normalized.  
!
!   MAXEXPBASE  (OUTPUT, OPTIONAL) integer(i4b)
!               The largest in magnitude positive integer such that         
!               float(base)**(maxexpbase) is positive and normalized. 
!
!   EPSPOS      (OUTPUT, OPTIONAL) real(stnd)         
!               The smallest power of BASE whose sum with 1. is         
!               greater than 1. That is, float(base)**(expepspos).                    
!
!   EPSNEG      (OUTPUT, OPTIONAL) real(stnd)          
!               The smallest power of BASE whose difference with 1.           
!               is less than 1. That is, float(base)**(expepsneg).                 
!
!   EPSILPOS    (OUTPUT, OPTIONAL) real(stnd)
!               The smallest positive floating point number whose           
!               sum with 1. is greater than 1.             
!
!   EPSILNEG    (OUTPUT, OPTIONAL) real(stnd)
!               The smallest positive floating point number whose           
!               difference with 1. is less than 1.             
!
!   RNDUNIT     (OUTPUT, OPTIONAL) real(stnd)
!               Unit roundoff of the machine.
!
!                                                                              
! Further Details
! _______________
!
!   This subroutine is based on the routines  MACHAR  by Cody and  DLAMCH  in
!   LAPACK77 (version 3). For further details, See:
!
!   (1) Malcolm M.A., 1972: Algorithms to reveal properties of
!          floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!   (2) Gentleman, W.M., and Marovich, S.B., 1974: More on algorithms
!          that reveal properties of floating point arithmetic units.
!          Comms. of the ACM, 17, 276-277.
!
!   (3) Cody, W.J., 1988: MACHAR: A subroutine to dynamically 
!           determine machine parameters,TOMS 14, No. 4, 303-311.                                 
!                                                                              
!                                                                              
! _____________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), optional :: basedigits, irnd, iuflow, igrd, iexp, ifloat,    &
                                           expepspos, expepsneg, minexpbase, maxexpbase
!
    real(stnd), intent(out),   optional :: epspos, epsneg, epsilpos, epsilneg, rndunit
!   
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b), save :: nuflow, nexp, nfloat, expposeps, expnegeps
    integer(i4b)       :: ia, ip, iexpsum
!
    real(stnd), save   :: poseps, negeps, posepsil, negepsil
    real(stnd)         :: baseinv, a, b, c, eps1, eps2
!
    logical(lgl)       :: first=true
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   THROUGHOUT THIS ROUTINE WE USE THE FUNCTION add TO ENSURE
!   THAT RELEVANT VALUES ARE STORED AND NOT HELD IN REGISTERS,
!   OR ARE NOT AFFECTED BY OPTIMIZERS.
!
!$OMP CRITICAL
!
    if ( first ) then
!
        first = false
!
!       comp_mach RETURNS MACHINE-SPECIFIC PARAMETERS IN ITS ARGUMENT LIST.
!
        if ( first_mach ) call comp_mach( nrnd, ngrd, unitrnd )
!
!       DETERMINE negeps AND expnegeps.
!
        baseinv = one/machbase
!
        negeps = spacing( one )*machbase
!
        do
!
            b = negeps*baseinv
            c = add( one, -b )
!
            if ( c==one ) exit
!
            negeps = b
!
        end do
!
        expnegeps = exponent( negeps ) - 1_i4b
!
!       DETERMINE negepsil, THE SMALLEST POSITIVE FLOATING POINT NUMBER
!       WHOSE DIFFERENCE WITH ONE IS LESS THAN ONE.
!
        eps1 = negeps
        eps2 = b
!
        do
!
            b = add( eps1, eps2 )/two
            c = add( one, -b )
!
            if ( c<one ) then
!
                if ( b==eps1 ) exit
!
                eps1 = b
!
            else
!
                if ( b==eps2 ) exit
!
                eps2 = b
!
            end if
!
        end do
!
        do
!
            b = nearest( eps1, -one )
            c = add( one, -b )
!
            if ( c>=one ) exit
!
            eps1 = b
!
            eps2 = nearest( eps2, one )
            c = add( one, -eps2 )
!
            if ( c<one ) exit
!
        end do
!
        negepsil = merge( eps1, eps2, c>=one )
!
!       DETERMINE poseps AND expposeps.
!
        poseps = spacing( two )*machbase
!
        do
!
            b = poseps*baseinv
            c = add( one, b )
!
            if ( c==one ) exit
!
            poseps = b
!
        end do
!
        expposeps = exponent( poseps ) - 1_i4b
!
!       DETERMINE posepsil, THE SMALLEST POSITIVE FLOATING POINT NUMBER
!       WHOSE SUM WITH ONE IS GREATER THAN ONE.
!
        eps1 = poseps
        eps2 = b
!
        do
!
            b = add( eps1, eps2 )/two
            c = add( one, b )
!
            if ( c>one ) then
!
                if ( b==eps1 ) exit
!
                eps1 = b
!
            else
!
                if ( b==eps2 ) exit
!
                eps2 = b
!
            end if
!
        end do
!
        do
!
            b = nearest( eps1, -one )
            c = add( one, b )
!
            if ( c<=one ) exit
!
            eps1 = b
!
            eps2 = nearest( eps2, one )
            c = add( one, eps2 )
!
            if ( c>one ) exit
!
        end do
!
        posepsil = merge( eps1, eps2, c<=one )
!
!       CHECK FOR PARTIAL UNDERFLOW. WE ASSUME PARTIAL UNDERFLOW
!       IF WE FOUND NON VANISHING DENORMALIZED FLOATING POINT NUMBERS.
!
        a = add( machtiny/two , zero )
        b = a*one
        c = add( b, b )
!
        nuflow = merge( 1_i4b, 0_i4b , (c/=zero) .and. (zero<a) .and. (a<machtiny) )
!
!       CHECK IF base IS A POWER OF 2.
!
        ia = base
        ip = 0_i4b
!
        do
!
            if ( ia==1_i4b ) exit
!
            ia = ia/2_i4b
            ip = ip + 1_i4b
!
        end do
!
        if ( mod( base, 2_i4b**ip )==0_i4b ) then
!
!           DETERMINE nexp THE NUMBER OF BITS NEEDED TO STORE THE EXPONENT.
!           WE ASSUME THAT THERE ARE nexp BINARY DIGITS SET ASIDE FOR 
!           THE REPRESENTATION OF THE EXPONENT. THAT IS  maxexp+abs(minexp)+1  
!           SUM APPROXIMATELY TO A POWER OF 2.
!
            iexpsum = maxexp - minexp + 1_i4b
!
            if ( nuflow==1_i4b .and. nrnd==2_i4b )  then
!
!               CORRECT iexpsum FOR IEEE MACHINES.
!
!               WE ASSUME IEEE ARITHMETIC IF WE FOUND DENORMALIZED NUMBERS,
!               IF ARITHMETIC SEEMS TO ROUND IN THE IEEE STYLE AND IF QUIET NaNs EXIST.
!               IF THIS IS TRUE, WE NEED TO AUGMENT iexpsum BY 2 (ONE EXPONENT FOR 
!               infinity AND NaNs AND ONE EXPONENT FOR zero AND DENORMALIZED NUMBERS).
!
!               CHECK FOR THE EXISTENCE OF NaNs.
!
                if ( test_nan() ) iexpsum = iexpsum + 2_i4b
!
            end if
!
            nexp = nint( log( real( iexpsum, stnd ) )/log( two ) , i4b )
!
            if ( base==2_i4b .and. log10(machhuge)>c2000 ) then
!
!               CORRECT nexp FOR CRAY MACHINES WITH ONE RESERVED EXPONENT 
!               (BINARY) DIGIT. WE ASSUME THAT WE ARE ON A CRAY IF THE OVERFLOW
!               TRESHOLD IS SUFFICIENTLY LARGE AND IF base=2 .
!
                nexp = nexp + 1_i4b
!
            end if
!
!           DETERMINE nfloat THE TOTAL NUMBER OF BITS NEEDED TO STORE A
!           FLOATING-POINT NUMBER.
!
            nfloat = 1_i4b + nexp + ( ip*nbasedigits )
!
!           EITHER THERE ARE AN ODD NUMBER OF BITS USED TO STORE A
!           FLOATING-POINT NUMBER, WHICH IS UNLIKELY, OR SOME BITS ARE
!           NOT USED IN THE REPRESENTATION OF NUMBERS, WHICH IS POSSIBLE
!           (E.G. CRAY MACHINES) OR THE MANTISSA HAS AN IMPLICIT BIT,
!           (E.G. IEEE MACHINES, DEC VAX MACHINES), WHICH IS PERHAPS THE
!           MOST LIKELY. WE ASSUME THAT WE ARE ON A CRAY IF THE OVERFLOW
!           TRESHOLD IS SUFFICIENTLY LARGE. IF THIS IS FALSE, WE ASSUME
!           THE LAST ALTERNATIVE AND WE NEED TO REDUCE nfloat BY ONE.
!
            if ( mod(nfloat,2_i4b)==1_i4b .and. base==2_i4b .and. log10(machhuge)<=c2000 ) then
                nfloat = nfloat - 1_i4b
            end if
!
!           CORRECT nfloat IF SOME BITS ARE NOT USED IN THE REPRESENTATION 
!           OF NUMBERS AS IN CRAY OR RS6000 MACHINES.
!
            select case( nfloat )
!
                case( 33_i4b:64_i4b )
                    nfloat = 64_i4b
                case( 65_i4b:128_i4b )
                    nfloat = 128_i4b
                case(129_i4b:)
                    nfloat = -1_i4b
!
            end select
!
        else
!
            nexp = -1_i4b
            nfloat = -1_i4b
!
        end if
!
    end if
!
!$OMP END CRITICAL
!
!   RETURN VALUES.
!
    if ( present( basedigits ) ) basedigits = nbasedigits
    if ( present( irnd ) )       irnd       = nrnd
    if ( present( iuflow ) )     iuflow     = nuflow
    if ( present( igrd ) )       igrd       = ngrd
    if ( present( iexp ) )       iexp       = nexp
    if ( present( ifloat ) )     ifloat     = nfloat
    if ( present( expepspos ) )  expepspos  = expposeps
    if ( present( expepsneg ) )  expepsneg  = expnegeps
    if ( present( minexpbase ) ) minexpbase = minexp - 1_i4b
    if ( present( maxexpbase ) ) maxexpbase = maxexp - 1_i4b
    if ( present( epspos ) )     epspos     = poseps
    if ( present( epsneg ) )     epsneg     = negeps
    if ( present( epsilpos ) )   epsilpos   = posepsil
    if ( present( epsilneg ) )   epsilneg   = negepsil
    if ( present( rndunit ) )    rndunit    = unitrnd
!
!
! END OF SUBROUTINE mach
! ______________________
!
    end subroutine mach
!
! ===========================================================================================
!
    function test_ieee(  )
!
! Purpose
! _______
!
!   TEST_IEEE try to determine if the computer follows the IEEE standard 754 
!   for binary floating-point arithmetic.
!
! Arguments
! _________
!
!   none
!
!                                                                              
! Further Details
! _______________
!
!   If the compiler follows the Fortran 2003 standard, the facilities provided by the
!   IEEE_ARITHMETIC module are used to determine if the computer follows the IEEE standard
!   754 for binary floating-point arithmetic.
!
!
! _____________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _F2003
    use ieee_arithmetic, only : ieee_support_datatype
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    logical(lgl) :: test_ieee
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
#ifdef _F2003
    real(stnd)   :: a
!   
    logical(lgl) :: ieee
#else
    real(stnd)   :: a, b, c
!
    integer(i4b)   :: ia, ip
!   
    logical(lgl) :: first=true, ieee=false
#endif
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
#ifdef _F2003
!
    ieee = ieee_support_datatype( a )
!
#else
!
!$OMP CRITICAL
!
    if ( first ) then
!
        first = false
!
!       FIRST CHECK IF base IS A POWER OF 2.
!
        ia = base
        ip = 0_i4b
!
        do
!
            if ( ia==1_i4b ) exit
!
            ia = ia/2_i4b
            ip = ip + 1_i4b
!
        end do
!
        if ( mod( base, 2_i4b**ip )==0_i4b ) then
!
!           SECOND CHECK FOR PARTIAL UNDERFLOW. WE ASSUME PARTIAL UNDERFLOW
!           IF WE FOUND NON VANISHING DENORMALIZED FLOATING POINT NUMBERS.
!
            a = add( machtiny/two , zero )
            b = a*one
            c = add( b, b )
!
            if ( (c/=zero) .and. (zero<a) .and. (a<machtiny) ) then
!
!               THIRD CHECK FOR THE EXISTENCE OF NaNs.
!
                if ( test_nan() ) then
!
!                   comp_mach RETURNS MACHINE-SPECIFIC PARAMETERS IN ITS ARGUMENT LIST.
!
                    if ( first_mach ) call comp_mach( nrnd, ngrd, unitrnd )
!
!                   WE ASSUME IEEE ARITHMETIC IF base IS A POWER OF 2,
!                   IF WE FOUND DENORMALIZED NUMBERS, IF QUIET NaNs EXIST
!                   AND IF ARITHMETIC SEEMS TO ROUND IN THE IEEE STYLE.
!
                    ieee = nrnd==2_i4b
!
                end if
!
            end if
!
        end if
!
    end if
!
!$OMP END CRITICAL
!
#endif
!
!   RETURN FUNCTION VALUE.
!
    test_ieee = ieee
!
!
! END OF FUNCTION test_ieee
! _________________________
!
    end function test_ieee
!
! =============================================================================================
! 
    function test_nan(  ) result( testnan )
!
! Purpose
! _______
!
!   TEST_NAN returns TRUE if NaNs exist, and FALSE otherwise.
!
!
! Arguments
! _________
!
!   None
!
!
! Further Details
! _______________
!
!   If the compiler follows the Fortran 2003 standard, the facilities provided by the
!   IEEE_ARITHMETIC module are used to determine if NaNs exist as defined in the IEEE standard
!   754 for binary floating-point arithmetic.
!
!   Otherwise, the routine exploits the IEEE requirement that NaNs
!   compare as unequal to all values, including themselves.
!
!   For further details, see:
!
!   (1) Cody, W.J., and Coonen, J.T., 1993: Algorithm 722,
!           TOMS 19, No. 4, 443-451.
!
!
! _____________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _F2003
    use ieee_arithmetic, only : ieee_support_datatype, ieee_support_nan
#endif
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd) :: a, b
!   
    logical(lgl) :: testnan
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
#ifdef _F2003
    if ( ieee_support_datatype( a ) ) then
        testnan = ieee_support_nan( a )
    else
        a = true_nan()
        b = true_nan()
!
        testnan = a/=b
    end if
#else
    a = true_nan()
    b = true_nan()
!
    testnan = a/=b
#endif
!
!
! END OF FUNCTION test_nan
! ________________________
!
    end function test_nan
!
! =============================================================================================
!
    function is_nan_r( x ) result( isnan )
!
! Purpose
! _______
!
!   This function returns TRUE if the scalar X is a NaN, and FALSE otherwise.
!
!
! Arguments
! _________
!
!   X   (INPUT) real(stnd)
!       The floating point number to be tested.
!
!
! Further Details
! _______________
!
!   If the compiler follows the Fortran 2003 standard, the facilities provided by the
!   IEEE_ARITHMETIC module are used to detect NaNs as defined in the IEEE standard
!   754 for binary floating-point arithmetic.
!
!   Otherwise, the routine exploits the IEEE requirement that NaNs
!   compare as unequal to all values, including themselves.
!
!   Finally, if the computer does not follow the IEEE standard 754 for binary floating-point
!   arithmetic, this function returns TRUE if the scalar X is equal to huge(X).
!
!   For further details, see:
!
!   (1) Cody, W.J., and Coonen, J.T., 1993: Algorithm 722,
!           TOMS 19, No. 4, 443-451.
!
!
! _____________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _F2003
    use ieee_arithmetic, only : ieee_support_datatype, ieee_is_nan
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in) :: x
!   
    logical(lgl) :: isnan
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( test_nan() ) then
!
#ifdef _F2003
        if ( ieee_support_datatype( x ) ) then
            isnan = ieee_is_nan( x )
        else
            isnan = diff_r( x, x )
        end if
#else
        isnan = diff_r( x, x )
!        isnan = x/=x
#endif
!
    else
!
        isnan = x==machhuge
!
    end if
!
!
! END OF FUNCTION is_nan_r
! ________________________
!
    end function is_nan_r
!
! =============================================================================================
!
    function is_nan_rv( x )  result( isnan )
!
! Purpose
! _______
!
!   This function returns the value TRUE if any of the elements of the vector X is a NaN, 
!   and FALSE otherwise.
!
!
! Arguments
! _________
!
!   X   (INPUT) real(stnd), dimension(:)
!       The floating point vector to be tested.
!
!
! Further Details
! _______________
!
!   If the compiler follows the Fortran 2003 standard, the facilities provided by the
!   IEEE_ARITHMETIC module are used to detect NaNs as defined in the IEEE standard
!   754 for binary floating-point arithmetic.
!
!   Otherwise, the routine exploits the IEEE requirement that NaNs
!   compare as unequal to all values, including themselves.
!
!   If the computer does not follow the IEEE standard 754 for binary floating-point
!   arithmetic, this function returns TRUE if any of the elements of the vector X is equal
!   to huge(X).
!
!   For further details, see:
!
!   (1) Cody, W.J., and Coonen, J.T., 1993: Algorithm 722,
!           TOMS 19, No. 4, 443-451.
!
!
! _____________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _F2003
    use ieee_arithmetic, only : ieee_support_datatype, ieee_is_nan
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:),intent(in) :: x
!   
    logical(lgl) :: isnan
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( test_nan() ) then
!
#ifdef _F2003
        if ( ieee_support_datatype( x ) ) then
            isnan = any( ieee_is_nan( x(:) ) )
        else
            isnan = diff_rv( x(:), x(:) )
        end if
#else
        isnan = diff_rv( x(:), x(:) )
!        isnan = any( x(:)/=x(:) )
#endif
!
    else
!
        isnan = any( x(:)==machhuge )
!
    end if
!
!
! END OF FUNCTION is_nan_rv
! _________________________
!
    end function is_nan_rv
!
! =============================================================================================
!
    function is_nan_rm( x )  result( isnan )
!
! Purpose
! _______
!
!   This function returns the value TRUE if any of the elements of the matrix X is a NaN, 
!   and FALSE otherwise.
!
!
! Arguments
! _________
!
!   X   (INPUT) real(stnd), dimension(:,:)
!       The floating point matrix to be tested.
!
!
! Further Details
! _______________
!
!   If the compiler follows the Fortran 2003 standard, the facilities provided by the
!   IEEE_ARITHMETIC module are used to detect NaNs as defined in the IEEE standard
!   754 for binary floating-point arithmetic.
!
!   Otherwise, the routine exploits the IEEE requirement that NaNs
!   compare as unequal to all values, including themselves.
!
!   If the computer does not follow the IEEE standard 754 for binary floating-point
!   arithmetic, this function returns TRUE if any of the elements of the matrix X is
!   equal to huge(X).
!
!   For further details, see:
!
!   (1) Cody, W.J., and Coonen, J.T., 1993: Algorithm 722,
!           TOMS 19, No. 4, 443-451.
!
!
! _____________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _F2003
    use ieee_arithmetic, only : ieee_support_datatype, ieee_is_nan
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:,:), intent(in) :: x
!   
    logical(lgl) :: isnan
!   
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( test_nan() ) then
!
#ifdef _F2003
        if ( ieee_support_datatype( x ) ) then
            isnan = any( ieee_is_nan( x(:,:) ) )
        else
            isnan = diff_rm( x(:,:), x(:,:) )
        end if
#else
        isnan = diff_rm( x(:,:), x(:,:) )
!        isnan = any( x(:,:)/=x(:,:) )
#endif
!
    else
!
        isnan = any( x(:,:)==machhuge )
!
    end if
!
!
! END OF FUNCTION is_nan_rm
! _________________________
!
    end function is_nan_rm
!
! =============================================================================================
!
    subroutine replace_nan_r( x, missing )
!
! Purpose
! _______
!
!   This subroutine replaces the scalar X with the scalar MISSING, if X is a NaN on input.
!
! Arguments
! _________
!
!   X           (INPUT/OUTPUT) real(stnd)
!               The floating point number to be tested.
!
!   MISSING     (INPUT) real(stnd)
!               The floating point number used to replace NaNs.
!
!
! Further Details
! _______________
!
!   If the compiler follows the Fortran 2003 standard, the facilities provided by the
!   IEEE_ARITHMETIC module are used to detect NaNs as defined in the IEEE standard
!   754 for binary floating-point arithmetic.
!
!   Otherwise, the routine exploits the IEEE requirement that NaNs
!   compare as unequal to all values, including themselves.
!
!   If the computer does not follow the IEEE standard 754 for binary floating-point
!   arithmetic, this subroutine replaces X with the scalar MISSING, if X is equal to
!   huge(X).
!
!   For further details, see:
!
!   (1) Cody, W.J., and Coonen, J.T., 1993: Algorithm 722,
!           TOMS 19, No. 4, 443-451.
!
!
! _____________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _F2003
    use ieee_arithmetic, only : ieee_support_datatype, ieee_is_nan
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout) :: x
    real(stnd), intent(in)    :: missing
!   
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( test_nan() ) then
!
#ifdef _F2003
        if ( ieee_support_datatype( x ) ) then
            if ( ieee_is_nan( x ) )  x = missing
        else
            if ( diff_r( x, x) )  x = missing
        end if
#else
        if ( diff_r( x, x) )  x = missing
!        if ( x/=x )  x = missing
#endif
!
    else
!
        if ( x==machhuge )  x = missing
!
    end if
!
!
! END OF SUBROUTINE replace_nan_r
! _______________________________
!
    end subroutine replace_nan_r
!
! =============================================================================================
!
    subroutine replace_nan_rv( x, missing )
!
! Purpose
! _______
!
!   This subroutine replaces the elements of the vector X which are NaNs by the scalar MISSING.
!
! Arguments
! _________
!
!   X        (INPUT/OUTPUT) real(stnd), dimension(:)
!            The floating point vector to be tested.
!
!   MISSING  (INPUT) real(stnd)
!            The floating point number used to replace the NaNs.
!
!
! Further Details
! _______________
!
!   If the compiler follows the Fortran 2003 standard, the facilities provided by the
!   IEEE_ARITHMETIC module are used to detect NaNs as defined in the IEEE standard
!   754 for binary floating-point arithmetic.
!
!   Otherwise, the routine exploits the IEEE requirement that NaNs
!   compare as unequal to all values, including themselves.
!
!   Finally, if the computer does not follow the IEEE standard 754 for binary floating-point
!   arithmetic, this subroutine replaces the elements of the vector X which are equal
!   to huge(X) with the scalar MISSING.
!
!   For further details, see:
!
!   (1) Cody, W.J., and Coonen, J.T., 1993: Algorithm 722,
!           TOMS 19, No. 4, 443-451.
!
!
! _____________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _F2003
    use ieee_arithmetic, only : ieee_support_datatype, ieee_is_nan
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout) :: x(:)
    real(stnd), intent(in)    :: missing
!   
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( test_nan() ) then
!
#ifdef _F2003
        if ( ieee_support_datatype( x ) ) then
            where( ieee_is_nan( x(:) ) )  x(:) = missing
        else
            where( maskdiff_rv(x(:),x(:)) )  x(:) = missing
        end if
#else
        where( maskdiff_rv(x(:),x(:)) )  x(:) = missing
!        where( x(:)/=x(:) )  x(:) = missing
#endif
!
    else
!
        where( x(:)==machhuge )  x(:) = missing
!
    end if
!
!
! END OF SUBROUTINE replace_nan_rv
! ________________________________
!
    end subroutine replace_nan_rv
!
! =============================================================================================
!
    subroutine replace_nan_rm( x, missing )
!
! Purpose
! _______
!
!   This subroutine replaces the elements of the matrix X which are NaNs by the scalar MISSING.
!
! Arguments
! _________
!
!   X        (INPUT/OUTPUT) real(stnd), dimension(:,:)
!            The floating point matrix to be tested.
!
!   MISSING  (INPUT) real(stnd)
!            The floating point number used to replace the NaNs.
!
!
! Further Details
! _______________
!
!   If the compiler follows the Fortran 2003 standard, the facilities provided by the
!   IEEE_ARITHMETIC module are used to detect NaNs as defined in the IEEE standard
!   754 for binary floating-point arithmetic.
!
!   Otherwise, the routine exploits the IEEE requirement that NaNs
!   compare as unequal to all values, including themselves.
!
!   Finally, if the computer does not follow the IEEE standard 754 for binary floating-point
!   arithmetic, this subroutine replaces the elements of the matrix X which are equal
!   to huge(X) with the scalar MISSING.
!
!   For further details, see:
!
!   (1) Cody, W.J., and Coonen, J.T., 1993: Algorithm 722,
!           TOMS 19, No. 4, 443-451.
!
!
! _____________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _F2003
    use ieee_arithmetic, only : ieee_support_datatype, ieee_is_nan
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout) :: x(:,:)
    real(stnd), intent(in)    :: missing
!   
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( test_nan() ) then
!
#ifdef _F2003
        if ( ieee_support_datatype( x ) ) then
            where( ieee_is_nan( x(:,:) ) )  x(:,:) = missing
        else
            where( maskdiff_rm(x(:,:),x(:,:)) )  x(:,:) = missing
        end if
#else
        where( maskdiff_rm(x(:,:),x(:,:)) )  x(:,:) = missing
!        where( x(:,:)/=x(:,:) )  x(:,:) = missing
#endif
!
    else
!
        where( x(:,:)==machhuge )  x(:,:) = missing
!
    end if
!
!
! END OF SUBROUTINE replace_nan_rm
! ________________________________
!
    end subroutine replace_nan_rm
!
! =============================================================================================
!
    function nan( ) result( x )
!
! Purpose
! _______
!
!   NAN returns as a scalar function, the bit pattern corresponding to a quiet NaN
!   in the IEEE standard 754 for binary floating-point arithmetic if the machine
!   recognizes NaNs or the maximum floating point number of kind STND otherwise.
!
!
! Arguments
! _________
!
!   None
!
!
! Further Details
! _______________
!
!   If the compiler follows the Fortran 2003 standard, the facilities provided by the
!   IEEE_ARITHMETIC module are used to create a quiet NaN as defined in the IEEE standard
!   754 for binary floating-point arithmetic.
!
!   Otherwise, the routine exploits the IEEE requirement that NaNs
!   compare as unequal to all values, including themselves.
!
!   Finally, NAN returns the maximum floating point number of kind STND, if the computer
!   does not follow the IEEE standard 754 for binary floating-point arithmetic.
!
!
! _____________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _F2003
    use ieee_arithmetic, only : ieee_support_datatype, ieee_quiet_nan, ieee_value
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!   
  real(stnd) :: x
!   
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: a
!   
!
! EXECUTABLE STATEMENTS
! _____________________
!
!
!   TEST THE EXISTENCE OF NaNs.
!
    if ( test_nan() ) then
!
!       RETURN A QUIET NaN.
!
#ifdef _F2003
        if ( ieee_support_datatype( x ) ) then
            a = ieee_value( a, ieee_quiet_nan )
        else
            a = true_nan()
        end if
#else
        a = true_nan()
#endif
!
        x = a
!
    else
!
!BUG:   SOME COMPILERS (e.g. PGF90) DO NOT IMPLEMENT NaNs
!       CORRECTLY. THUS, IT IS SAFER TO USE THE VALUE huge()
!       INSTEAD OF A TRUE QUIET NaN.
!
!       RETURN huge(a).
!
        x = machhuge
!
    end if
!
!
! END OF FUNCTION nan
! ___________________
!
    end function nan
!
! =============================================================================================
!
    function true_nan( ) result( x )
!
! Purpose
! _______
!
!   TRUE_NAN returns as a scalar function, the bit pattern corresponding to a quiet NaN
!   in the IEEE standard 754 for binary floating-point arithmetic.
!
!
! Arguments
! _________
!
!   None
!
!
! _____________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!   
  real(stnd) :: x
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: idum(8)
!
    real(stnd)   :: nan_tmp
!   
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   SET THE BIT PATTERN FOR A QUIET NaN.
!
    idum(:) = not( 0_i4b )
    nan_tmp = transfer( idum, nan_tmp )
!
!   RETURN A QUIET NaN.
!
    x = nan_tmp
!
!
! END OF FUNCTION true_nan
! ________________________
!
    end function true_nan
!
!
! =========================================================================================
!                                 PRIVATE SUBROUTINES
! =========================================================================================
!
!
    subroutine comp_mach( nrnd, ngrd, unitrnd )
!
! Purpose
! _______
!
!   COMP_MACH is a service routine for LAMCH, MACH and TEST_IEEE. COMP_MACH
!   determines the machine parameters specified in its argument list.
!  
!
! Arguments
! _________
!
!   NRND        (OUTPUT) integer(i4b)
!               A parameter indicating whether proper rounding or chopping 
!               (rounding upward, downward, toward zero) occurs in addition: 
!
!               - NRND = 0 if floating-point addition rounds upward, downward
!                 or toward zero;
!               - NRND = 1 if floating-point addition rounds to nearest,
!                 but not in the IEEE style;
!               - NRND = 2 if floating-point addition rounds in the IEEE style.
!
!   NGRD        (OUTPUT) integer(i4b)
!               The number of guard digits for multiplication with chopping
!               arithmetic (IRND = 0). If:
!
!               - NGRD = 0    if floating-point arithmetic rounds, or if it
!                 chops and only  NBASEDIGITS  digits participate
!                 in the post-normalization shift of the
!                 floating-point significand in multiplication;
!               - NGRD = 1    if floating-point arithmetic chops and more
!                 than  NBASEDIGITS  digits participate in the
!                 post-normalization shift of the floating-point
!                 significand in multiplication.
!
!   UNITRND     (OUTPUT) real(stnd)
!               unit roundoff of the machine.
!
!
! Further Details
! _______________
!
!   This subroutine is based on the routines  MACHAR  by Cody and  DLAMCH  in
!   LAPACK77 (version 3). For further details, see:
!
!   (1) Malcolm, M.A., 1972: Algorithms to reveal properties of
!           floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!   (2) Gentleman, W.M., and Marovich, S.B., 1974: More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
!   (3) Cody, W.J., 1988: MACHAR: A subroutine to dynamically 
!         determine machine parameters,TOMS 14, No. 4, 303-311.                                 
!                                                                              
!                                                                              
! _____________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out) :: nrnd, ngrd
!
    real(stnd),   intent(out) :: unitrnd
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd) :: a, b, c, baseh, tempa
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!
!   THROUGHOUT THIS ROUTINE WE USE THE FUNCTION add TO ENSURE
!   THAT RELEVANT VALUES ARE STORED AND NOT HELD IN REGISTERS,
!   OR ARE NOT AFFECTED BY OPTIMIZERS.
!
    first_mach = false
!   
!   NOW DETERMINE THE MANTISSA, nbasedigits. IT MAY BE DETERMINE BY
!   POWERING. nbasedigits IS THE SMALLEST POSITIVE INTEGER
!   FOR WHICH
!
!          real( base**nbasedigits + 1 ) - real( base**nbasedigits ) /= 1
!   
!   HOWEVER, IT IS SAFER TO USE THE INTRINSIC FUNCTION digits() .
!
!    nbasedigits = 0_i4b
!    a = one
!    c = one
!    do while( c==one )
!        nbasedigits = nbasedigits + 1_i4b
!        a = a * machbase
!        b = add( a, one )
!        c = add( b, -a )
!    end do
!
    a = machbase**int(nbasedigits)
!
!   DETERMINE WETHER ROUNDING OR CHOPPING OCCURS, BY ADDING
!   A BIT LESS THAN base/2 AND A BIT MORE THAN base/2 TO a.
!
    baseh = machbase/two
    b = nearest( baseh, -one )
    c = add( a, b )
    nrnd = merge( 1_i4b, 0_i4b, c==a )
    b = nearest( baseh, one )
    c = add( a, b )
!   
    if ( nrnd==1_i4b .and. c==a ) nrnd = 0_i4b
!   
    if ( nrnd==1_i4b ) then
!
        ngrd = 0_i4b
!
!       COMPUTE UNIT ROUNDOFF IF ROUNDING OCCURS.
!
        unitrnd = macheps/two
!
!       TRY AND DECIDE WHETHER ROUNDIND ID DONE IN THE IEEE 'ROUND
!       TO NEAREST' STYLE. base/2 IS HALF A UNIT IN THE LAST PLACE
!       OF THE TWO NUMBERS a AND a+base. FURTHERMORE, a IS EVEN,
!       I.E. HAS LAST BIT zero, AND a+base IS ODD, I.E. HAS LAST
!       BIT one. THUS ADDING base/2 TO a SHOULD NOT CHANGE a, BUT
!       ADDING base/2 TO a+base SHOULD CHANGE a+base.
!
        tempa = add( a, machbase )
        b = add( a, baseh )
        c = add( tempa, baseh)
!   
        if ( a==b .and. c>tempa ) nrnd = 2_i4b
!
    else
!
!       COMPUTE UNIT ROUNDOFF IF CHOPPING OCCURS.
!
        unitrnd = macheps
!
!       DETERMINE THE NUMBER OF GUARD DIGITS FOR
!       MULTIPLICATION WITH TRUNCATING ARITHMETIC.
!
        b = nearest( one, one )
        c = add( b*one, zero )
        ngrd = merge( 1_i4b, 0_i4b, c/=one )
!
    end if
!
!
! END OF SUBROUTINE comp_mach
! ___________________________
!
    end subroutine comp_mach
!
! ===========================================================================================
!
    real(stnd) function add( y, z )
!
! Purpose
! _______
!
!   ADD  is intended to force  Y  and  Z  to be stored prior to doing
!   the addition of  Y  and  Z ,  for use in situations where optimizers
!   might hold one of these in a register.
!
!
! Arguments
! _________
!
!   Y, Z    (INPUT) real(stnd)
!           The values Y and Z.
!                                                 
!
! _____________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in) :: y, z
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    add = y + z
!
!
! END OF FUNCTION add
! ___________________
!
    end function add
!
! =============================================================================================
! 
    function diff_r( y, z )
!
! Purpose
! _______
!
!   DIFF_R is intended to force the scalars  Y  and  Z  to be stored prior to doing
!   the comparison of  Y  and  Z ,  for use in situations where optimizers
!   might hold one of these in a register.
!
!
! Arguments
! _________
!
!   Y, Z    (INPUT) real(stnd)
!           The values y and z to be compared.
!                                                                              
!                                                                              
! _____________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in) :: y, z
!
    logical(lgl) :: diff_r
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    diff_r = y/=z
!
!
! END OF FUNCTION diff_r
! ______________________
!
    end function diff_r
!
! =============================================================================================
! 
    function diff_rv( y, z )
!
! Purpose
! _______
!
!   DIFF_RV is intended to force the vectors  Y(:)  and  Z(:)  to be stored prior to doing
!   the comparison of  Y  and  Z ,  for use in situations where optimizers
!   might hold one of these in a register.
!
!
! Arguments
! _________
!
!   Y, Z    (INPUT) real(stnd), dimension(:)
!           The vectors X and Z to be compared.
!                                                                              
!                                                                              
! _____________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:), intent(in) :: y, z
!
    logical(lgl) :: diff_rv
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: ny, nz, n
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    ny = size( y )
    nz = size( z )
    n = min( ny, nz )
!
    diff_rv = any( y(:n)/=z(:n) )
!
!
! END OF FUNCTION diff_rv
! _______________________
!
    end function diff_rv
!
! =============================================================================================
! 
    function diff_rm( y, z )
!
! Purpose
! _______
!
!   DIFF_RM is intended to force the matrices  Y(:,:)  and  Z(:,:)  to be stored prior to doing
!   the comparison of  Y  and  Z ,  for use in situations where optimizers
!   might hold one of these in a register.
!
!
! Arguments
! _________
!
!   Y, Z    (INPUT) real(stnd), dimension(:,:)
!           The matrices Y and Z to be compared.
!                                                                              
!
! _____________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:,:), intent(in) :: y, z
!
    logical(lgl) :: diff_rm
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: ny, nz, n1, n2
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    ny = size( y, 1 )
    nz = size( z, 1 )
    n1 = min( ny, nz )
!
    ny = size( y, 2 )
    nz = size( z, 2 )
    n2 = min( ny, nz )
!
    diff_rm = any( y(:n1,:n2)/=z(:n1,:n2) )
!
!
! END OF FUNCTION diff_rm
! _______________________
!
    end function diff_rm
!
! =============================================================================================
! 
    function maskdiff_rv( y, z )
!
! Purpose
! _______
!
!   MASKDIFF_RV  is intended to force  Y(:)  and  Z(:)  to be stored prior to doing
!   the comparison of  Y  and  Z ,  for use in situations where optimizers
!   might hold one of these in a register.
!
!
! Arguments
! _________
!
!   Y, Z    (INPUT) real(stnd), dimension(:)
!           The vectors Y and Z to be compared. It is assumed that Y and Z have
!           the same size.
!                                                                              
!                                                                              
! _____________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:), intent(in) :: y, z
!
    logical(lgl), dimension(size(y)) :: maskdiff_rv
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    n = size( y )
!
    maskdiff_rv(:n) = y(:n)/=z(:n)
!
!
! END OF FUNCTION maskdiff_rv
! ___________________________
!
    end function maskdiff_rv
!
! =============================================================================================
! 
    function maskdiff_rm( y, z )
!
! Purpose
! _______
!
!   MASKDIFF_RM  is intended to force  Y(:,:)  and  Z(:,:)  to be stored prior to doing
!   the comparison of  Y  and  Z ,  for use in situations where optimizers
!   might hold one of these in a register.
!
!
! Arguments
! _________
!
!   Y, Z    (INPUT) real(stnd), dimension(:,:)
!           The matrices Y and Z to be compared. It is assumed that Y and Z have
!           the same shape.
!
!                                                                              
! _____________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:,:), intent(in) :: y, z
!
    logical(lgl), dimension(size(y,1),size(y,2)) :: maskdiff_rm
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n1, n2
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    n1 = size( y, 1 )
    n2 = size( y, 2 )
!
    maskdiff_rm(:n1,:n2) = y(:n1,:n2)/=z(:n1,:n2)
!
!
! END OF FUNCTION maskdiff_rm
! ___________________________
!
    end function maskdiff_rm
!
! =============================================================================================
!
!
! *****************************
! END OF MODULE Num_Constants *
! *****************************
!    
end module Num_Constants
!    
! ***********************************************************************************************
!
module Char_Constants
! ##############################################################################
!
! Copyright 2020 IRD
!
! This file is part of statpack.
!
! statpack is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of 
! the License, or (at your option) any later version.
!
! statpack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You can find a copy of the GNU Lesser General Public License
! in the statpack/doc directory.
!
! ##############################################################################
!                                                                              *
! ******************************************************************************
! MODULE EXPORTING CHARACTER CONSTANTS, STRINGS AND ERROR MESSAGES FOR         *
! ROUTINES AVAILABLE IN STATPACK.                                              *
!                                                                              *
! LATEST REVISION : 04/12/2020                                                 *
!                                                                              *
! ##############################################################################
!                                                                              *
! ******************************************************************************
!
!
! STRONG TYPING IMPOSED 
! =====================
!    
    implicit none
!
! PUBLIC ENTITIES 
! ===============
!
    public
!    
! MODULE PARAMETERS 
! =================
!
! NAMES FOR COMMON CHARACTERS.
!
    character(len=1), parameter ::    &
        ampersand   = achar(38)  ,    &
        apostrophe  = achar(39)  ,    &
        atSign      = achar(64)  ,    &
        backslash   = achar(92)  ,    &
        backquote   = achar(96)  ,    &
        bang        = achar(33)  ,    &
        blank       = achar(32)  ,    &
        caret       = achar(94)  ,    &
        cbrace      = achar(125) ,    &
        cbracket    = achar(93)  ,    &
        cparen      = achar(41)  ,    &
        colon       = achar(58)  ,    &
        comma       = achar(44)  ,    &
        dash        = achar(45)  ,    &
        dollar      = achar(36)  ,    &
        equals      = achar(61)  ,    &
        exclamation = achar(33)  ,    &
        greaterthan = achar(62)  ,    &
        hash        = achar(35) 
!
    character(len=1), parameter ::    &
        lessthan    = achar(60)  ,    &
        minus       = achar(45)  ,    &
        obrace      = achar(123) ,    &
        obracket    = achar(91)  ,    &
        oparen      = achar(40)  ,    &
        percent     = achar(37)  ,    &
        period      = achar(46)  ,    &
        plus        = achar(43)  ,    &
        quesmark    = achar(63)  ,    &
        quote       = achar(34)  ,    &
        semicolon   = achar(59)  ,    &
        slash       = achar(47)  ,    &
        star        = achar(42)  ,    &
        tilde       = achar(126) ,    &
        vertBar     = achar(124) ,    &
        underscore  = achar(95)
!
! NAME FOR NULL STRING.
!
    character(len=0), parameter :: null = ''
!
! NAMES FOR ASCII COMMAND CHARACTERS.
!
    character(len=1), parameter ::    &
        bell = achar(7) ,             & ! BELL
        bs   = achar(8) ,             & ! BACK SPACE
        ht   = achar(9) ,             & ! HORIZONTAL TABULATION
        lf   = achar(10) ,            & ! LINE FEED
        vt   = achar(11) ,            & ! VERTICAL TABULATION
        ff   = achar(12) ,            & ! FORM FEED
        cr   = achar(13) ,            & ! CARRIAGE RETURN
        so   = achar(14) ,            & ! SHIFT OUT
        si   = achar(15) ,            & ! SHIFT IN
        esc  = achar(27) ,            & ! ESCAPE
        del  = achar(127)               ! DELETE
!
! ERROR MESSAGE FOR allocate STATEMENT .
!
    character(len=*), parameter ::                                   &
    allocate_error  = ' : problem in attempt to allocate memory !'
!
! ERROR MESSAGES FOR MODULE Utilities .
!
    character(len=*), parameter ::                                              &
    utilities_error1  = 'an assertion failed with this tag : ',                 &
    utilities_error2  = 'an assert_eq failed with this tag : ',                 &
    utilities_error3  = 'program terminated by MERROR !',                       &
    utilities_error4  = ' : CFROM==zero !'
!
! ERROR MESSAGES FOR MODULE Random .
!
    character(len=*), parameter ::                                                       &
    random_error1  = ' : 32-bits integers not available !',                              &
    random_error2  = ' : Default RNG must be 1, 2 or 3!',                                &
    random_error3  = ' : Too many arguments !',                                          &
    random_error4  = ' : No system_clock !',                                             &
    random_error5  = ' : Size(PUT) < seed_size !',                                       &
    random_error6  = ' : Size(GET) < seed_size !',                                       &
    random_error7  = ' : Nonstandard integer representation !',                          &
    random_error8  = ' : Integer overflows are unsafe !',                                &
    random_error9  = ' : Default RNG must be 3 if 32-bit integers are not available !',  &
    random_error10 = ' : KISS RNGs cannot be used if 32-bit integers are not available !'
!
! ERROR MESSAGES FOR MODULE Time_Procedures .
!
    character(len=*), parameter ::                                                  &
    time_error1  = ' : only valid dates in the  Gregorian  calendar are allowed !', &
    time_error2  = ' : min(size(T1),size(T0)) < 7  !'
!
! ERROR MESSAGES FOR MODULE QR_Procedures .
!
    character(len=*), parameter ::                                                  &
    qr_error1  = ' : estimated rank of MAT < min(KRANK,size(MAT,1),size(MAT,2)) !', &
    qr_error2  = ' : CHUNK <= 0  !',                                                &
    qr_error3  = ' : Optional IP array argument is missing !'
!
! ERROR MESSAGES FOR MODULE Lin_Procedures .
!
    character(len=*), parameter ::                                                  &
    lin_error1  = ' : MAT is a singular matrix !'
!
! ERROR MESSAGES FOR MODULE Eig_Procedures .
!
    character(len=*), parameter ::                                                     &
    eig_error1  = ' : MAT do not contain elementary reflectors !',                     &
    eig_error2  = ' : The eigenvalues are not in decreasing order !',                  &
    eig_error3  = ' : Only one of the optional arguments LE and THETA must be used !'
!
! ERROR MESSAGES FOR MODULE SVD_Procedures .
!
    character(len=*), parameter ::                                                         &
    svd_error1  = ' : A or B have more columns than rows !',                               &
    svd_error2  = ' : Only one of the optional arguments LS and THETA must be used !',     &
    svd_error3  = ' : The singular values are not in decreasing order !',                  &
    svd_error4  = ' : Error in the parallel path !',                                       &
    svd_error5  = ' : size(mat,1) < size(mat,2) is not allowed !',                         &
    svd_error6  = ' : Initial tridiagonal matrix is not positive-definite',                &
    svd_error7  = ' : The singular values must be positive or zero !',                     &
    svd_error8  = ' : abs(p(1,1)) is not equal to one !',                                  &
    svd_error9  = ' : RELERR must be greater than 4.*epsilon( RELERR ) and less than one !'
!
! ERROR MESSAGES FOR MODULE LLSQ_Procedures .
!
    character(len=*), parameter ::                                                    &
    llsq_error1  = ' : size(MAT)=0 !',                                                &
    llsq_error2  = ' : size(MAT)=0 or size(B)=0 !',                                   &
    llsq_error3  = ' : size(VEC)=0 !',                                                &
    llsq_error4  = ' : size(VEC)=0 or size(B)=0 !'
!
! ERROR MESSAGES FOR MODULE Prob_Procedures .
!
    character(len=*), parameter ::                                                    &
    prob_error1  = ' : P outside range (0, 1) !',                                     &
    prob_error2  = ' : Number of degrees of freedom < 1 !',                           &
    prob_error3  = ' : Number of degrees of freedom < 1 or P outside range (0, 1) !', &
    prob_error4  = ' : Number of degrees of freedom < 1 or X2 < 0  !',                &
    prob_error5  = ' : Numbers of degrees of freedom < 1 or F < 0  !',                &
    prob_error6  = ' : X <= 0  !',                                                    &
    prob_error7  = ' : Either A or B <= 0 !',                                         &
    prob_error8  = ' : X outside range (0, 1) !',                                     &
    prob_error9  = ' : ACU <= 0 !',                                                   &
    prob_error10 = ' : Either DF1 or DF2 <= 0 !',                                     &
    prob_error11 = ' : F < 0  !',                                                     &
    prob_error12 = ' : X < 0  !',                                                     &
    prob_error13 = ' : GAMP<=0 !',                                                    &
    prob_error14 = ' : X2<0 !',                                                       &
    prob_error15 = ' : DF<0.5 !',                                                     &
    prob_error16 = ' : DF<0.5 or DF>2.E5 !'  ,                                        &
    prob_error17 = ' : Maximum number of iterations exceeded !',                      &
    prob_error18 = ' : GAMP<0.25 !',                                                  &
    prob_error19 = ' : Sample size <= 1 or X <= 0  !',                                &
    prob_error20 = ' : N <= 0  !',                                                    &
    prob_error21 = ' : K < 0 or K>N  !',                                              &
    prob_error22 = ' : Either A or B <= 0.1 !',                                       &
    prob_error23 = ' : Either DF1 or DF2 <= 0.2 !'
!
! ERROR MESSAGES FOR MODULE Stat_Procedures AND Mul_Stat_Procedures.
!
    character(len=*), parameter ::                                      &
    stat_error1  = ' : Number of variables < 1 !',                      &
    stat_error2  = ' : CORTEST<=0 or CORTEST>=1 !',                     &
    stat_error3  = ' : UTEST<=0 or UTEST>=1 !',                         &
    stat_error4  = ' : Some groups are empty !',                        &
    stat_error5  = ' : DIMVAR>2 or DIMVAR<0 !',                         &
    stat_error6  = ' : Some singular values are <= 0 !',                &
    stat_error7  = ' : The number of observations is zero !',           &
    stat_error8  = ' : The number of shuffles is <=0 !',                &
    stat_error9  = ' : METHOD=1 or METHOD=2 !',                         &
    stat_error10 = ' : PERIODICITY<=0 or PERIODICITY>=size(Y) !',       &
    stat_error11 = ' : BLOCK_SIZE<=0 or BLOCK_SIZE>=size(Y) !',         &
    stat_error12 = ' : BLOCK_SIZE<PERIODICITY !',                       &
    stat_error13 = ' : SEASON<=0 !',                                    &
    stat_error14 = ' : size(Y) is not a multiple of SEASON !',          &
    stat_error15 = ' : BLOCK_SIZE>SEASON !',                            &
    stat_error16 = ' : PERIODICITY>SEASON !'
!
! ERROR MESSAGES FOR MODULE FFT_Procedures.
!
    character(len=*), parameter ::                                                      &
    fft_error1  = ' : values in SHAP must be greater than 0 !',                         &
    fft_error2  = ' : DIM must be greater than 0 and less or equal to size(SHAP) !',    &
    fft_error3  = ' : size(SHAP) must be less or equal to 3 !',                         &
    fft_error4  = ' : LENGTH1<=0 !',                                                    &
    fft_error5  = ' : LENGTH1<=0 or LENGTH2<=0 ! ',                                     &
    fft_error6  = ' : LENGTH1<=0 or LENGTH2<=0 or LENGTH3<=0 ! ',                       &
    fft_error7  = ' : size(VEC) is not an even integer !',                              &
    fft_error8  = ' : Wrong initialization, use INIT_FFT first !',                      &
    fft_error9  = ' : size(MAT,2) is not an even integer !',                            &
    fft_error10 = ' : DIM>2 or DIM<0 !',                                                &
    fft_error11 = ' : DIM must be equal to 1 or 2 !',                                   &
    fft_error12 = ' : DIM must be equal to 1, 2 or 3 !'
!
! ERROR MESSAGES FOR MODULE Time_Series_Procedures.
!
    character(len=*), parameter ::                                                     &
    tseries_error1  = ' : The smoothing factor must be >0 !',                          &
    tseries_error2  = ' : DIMVAR>2 or DIMVAR<=0 !',                                    &
    tseries_error3  = ' : The ITDEG parameter must be 0, 1 or 2 !',                    &
    tseries_error4  = ' : The NP parameter must be >1 !',                              &
    tseries_error5  = ' : The ISDEG parameter must be 0 or 1 !',                       &
    tseries_error6  = ' : The ILDEG parameter must be 0, 1 or 2 !',                    &
    tseries_error7  = ' : The NO parameter should be a nonnegative integer !',         &
    tseries_error8  = ' : The NI parameter must be >0 !',                              &
    tseries_error9  = ' : The LEN parameter must be >=1 and < size(X) !',              &
    tseries_error10 = ' : Number of points must be at least 4 !',                      &
    tseries_error11 = ' : Filter periods cannot be negative !',                        &
    tseries_error12 = ' : Periods shorter than 2 sample points cannot be observed !',  &
    tseries_error13 = ' : Periods longer than signal length cannot be observed !',     &
    tseries_error14 = ' : size(VEC) is not an even integer !',                         &
    tseries_error15 = ' : Unimplemented window function !',                            &
    tseries_error16 = ' : size(MAT,2) is not an even integer !',                       &
    tseries_error17 = ' : L is greater than size(VEC) !',                              &
    tseries_error18 = ' : L is not a positive even integer !',                         &
    tseries_error19 = ' : L is greater than size(MAT,2) !',                            &
    tseries_error20 = ' : Optional parameter SYM must be equal to -one, one or zero !'
!    
    character(len=*), parameter ::                                                                             &
    tseries_error21 = ' : size(MAT,1) is less or equal to 0 !',                                                &
    tseries_error22 = ' : MAX_ALLOC is less than 1 or greater than size(MAT,1) !',                             &
    tseries_error23 = ' : Length of data window must be a positive even integer !',                            &
    tseries_error24 = ' : L0 is not a positive integer !',                                                     &
    tseries_error25 = ' : L+L0 is not a positive even integer !',                                              &
    tseries_error26 = ' : The number K of filter coefficients must be greater or equal to 3 !',                &
    tseries_error27 = ' : The number K of filter coefficients must be odd !',                                  &
    tseries_error28 = ' : The cutoff frequency FC must be greater than zero and less than half !',             &
    tseries_error29 = ' : The cutoff frequency FCL must be greater than zero and less than half !',            &
    tseries_error30 = ' : The cutoff frequency FCH must be greater than zero and less than half !',            &
    tseries_error31 = ' : The cutoff period PL must be greater than 2 !',                                      &
    tseries_error32 = ' : The cutoff period PH must be greater than 2 !',                                      &
    tseries_error33 = ' : FC, minus one over the number of filter terms, K, must be greater or equal to 0 !',  &
    tseries_error34 = ' : FCH minus one over the number of filter terms, K, must be greater or equal to 0 !',  &
    tseries_error35 = ' : FC plus one over the number of filter terms, K, must be less than 0.5 !',            &
    tseries_error36 = ' : FCL plus one over the number of filter terms, K, must be less than 0.5 !',           &
    tseries_error37 = ' : FCL minus 1.3/(K+1) must be greater or equal to FCH plus 1.3/(K+1) !',               &
    tseries_error38 = ' : The frequency FREQ must be greater than zero and less than half !',                  &
    tseries_error39 = ' : The frequency FREQ must be greater or equal to ( 1.3/(K+1) + 1/K ) !',               &
    tseries_error40 = ' : The frequency FREQ must be less than 0.5 - ( 1.3/(K+1) + 1/K ) !'
!    
    character(len=*), parameter ::                                                                                   &
    tseries_error41 = ' : The frequency FREQ must be greater or equal to 1.3/(K+1) !',                               &
    tseries_error42 = ' : The frequency FREQ must be less than 0.5 - 1.3/(K+1) !',                                   &
    tseries_error43 = ' : size(COEF) must be greater or equal to 3 and less or equal to size(VEC) !',                &
    tseries_error44 = ' : size(COEF) must be greater or equal to 3 and less or equal to size(MAT,2) !',              &
    tseries_error45 = ' : size(COEF) must be odd !',                                                                 &
    tseries_error46 = ' : COEF(:)==zero !',                                                                          &
    tseries_error47 = ' : The input filter COEF(:) must be symmetric !',                                             &
    tseries_error48 = ' : size(COEF) must be greater or equal to 3 !',                                               &
    tseries_error49 = ' : FCL must be greater than FCH !',                                                           &
    tseries_error50 = ' : size(SMOOTH_PARAM) must be greater or equal to 1 !',                                       &
    tseries_error51 = ' : Elements in SMOOTH_PARAM(:) must be greater than zero !',                                  &
    tseries_error52 = ' : Elements in SMOOTH_PARAM(:) must be greater than zero and less than size(VEC) !',          &
    tseries_error53 = ' : Elements in SMOOTH_PARAM(:) must be greater than zero and less than size(MAT,2) !',        &
    tseries_error54 = ' : Elements in SMOOTH_PARAM(:) must be greater than zero and less than size(VEC)/2 + 1 !',    &
    tseries_error55 = ' : Elements in SMOOTH_PARAM(:) must be greater than zero and less than size(MAT,2)/2 + 1 !',  &
    tseries_error56 = ' : Elements in SMOOTH_PARAM(:) must be greater than zero and less than ((L+L0)/2) + 1  !',    &
    tseries_error57 = ' : NSEG is a strictly positive integer  !',                                                   &
    tseries_error58 = ' : Computation of dof with overlapped segments and this data window is not implemented  !',   &
    tseries_error59 = ' : PROBTEST<=0 or PROBTEST>=1 !',                                                             &
    tseries_error60 = ' : NSMOOTH must be greater than 2 and less or equal to size(VEC) !'
!    
    character(len=*), parameter ::                                                                                   &
    tseries_error61 = ' : NSMOOTH must be odd !',                                                                    &
    tseries_error62 = ' : NSMOOTH must be greater than 2 and less or equal to size(MAT,2) !',                        &
    tseries_error63 = ' : size(WK)+L0 must be even !',                                                               &
    tseries_error64 = ' : NSMOOTH must be greater than 2 and less or equal to (size(WK)+L0)/2 + 1 !',                &
    tseries_error65 = ' : NSMOOTH must be greater than 2 and less or equal to (size(VEC)/2)+1 !',                    &
    tseries_error66 = ' : NSMOOTH must be greater than 2 and less or equal to (size(MAT,2)/2)+1 !',                  &
    tseries_error67 = ' : NSMOOTH must be greater than 2 and less or equal to ((L+L0)/2)+1 !',                       &
    tseries_error68 = ' : PINTERVAL<=0 or PINTERVAL>=1 !',                                                           &
    tseries_error69 = ' : Either EDOFN or EDOFD <= 0 !',                                                             &
    tseries_error70 = ' : Either EDOFN or EDOFD <= 1 !',                                                             &
    tseries_error71 = ' : size(PSVEC1) must be greater or equal to 10 !',                                            &
    tseries_error72 = ' : size(PSMAT1,2) must be greater or equal to 10 !',                                          &
    tseries_error73 = ' : Invalid spectral estimates !',                                                             &
    tseries_error74 = ' : size(PSVECN) must be greater or equal to 2 !',                                             &
    tseries_error75 = ' : size(PSMATN,2) must be greater or equal to 2 !',                                           &
    tseries_error76 = ' : Elements in SMOOTH_PARAM(:) must be greater than zero and less than (size(WK)+L0)/2 + 1 !'
!
! ERROR MESSAGES FOR MODULE SVD_Weight .
!
    character(len=*), parameter ::                                                     &
    svdw_error1  = ' : irank>min(nvar,nobs) or irank<=0 !',                            &
    svdw_error2  = ' : The statement ALLOCATE causes STATUS = ',                       &
    svdw_error3  = ' : Some of the weights are negative !',                            &
    svdw_error4  = ' : Some variables have only missing values !',                     &
    svdw_error5  = ' : Some observations have only missing values !',                  &
    svdw_error6  = ' : Some variables have less than irank values !',                  &
    svdw_error7  = ' : Some observations have less than irank values !',               &
    svdw_error8  = ' : Some missing values have a non-zero weight !',                  &
    svdw_error9  = ' : Weighted data matrix all zero !',                               &
    svdw_error10 = ' : One column of initialization matrix for right_mat all zero !',  &
    svdw_error11 = ' : a deficient matrix was encountered, use pivoting !',            &
    svdw_error12 = ' : a deficient matrix was encountered, use Marquardt !',           &
    svdw_error13 = ' : B = 0 during iterations !',                                     &
    svdw_error14 = ' : Marquardt parameter is too large !   ',                         &
    svdw_error15 = ' : Gauss-Newton parameter is too small !'
!
! COMMON STRINGS FOR MODULE SVD_Weight .
!
    character(len=*), parameter ::                      &
    svdw_str1='Norm of gradient of F at (A,B) = ',      &
    svdw_str2='Test of gradient of F at (A,B) = '
!    
!    
! ******************************
! END OF MODULE Char_Constants *
! ******************************
!    
end module Char_Constants
