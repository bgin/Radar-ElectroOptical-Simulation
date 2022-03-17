
! ####################################################################################
! Modified by Bernard Gingold, beniekg@gmail.com
! Copyright 2020 IRD
!
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
! ####################################################################################
!                                                                                    *
! ************************************************************************************
! MODULE EXPORTING SUBROUTINES AND FUNCTIONS FOR TIME SERIES ANALYSIS                *
!                                                                                    *
! LATEST REVISION : 16/09/2020                                                       *
!                                                                                    *
! ####################################################################################
!                                                                                    *
! ************************************************************************************
!                                                                                    *
! THE C PROCESSOR MACROS USED IN THIS MODULE ARE:                                    *
!                                                                                    *
!  _OPENMP         FOR ACTIVATING OPENMP PARALLELIZATION                             *
!  _DOT_PRODUCT    FOR REPLACING THE dot_product INTRINSIC FUNCTION WITH STATPACK    *
!                  dot_product2 FUNCTION                                             *
!  _MATMUL         FOR REPLACING THE matmul INTRINSIC FUNCTION WITH STATPACK         *
!                  matmul2 FUNCTION                                                  *
!  _INTERNAL_PROC  FOR DESACTIVATING OPENMP PARALLELIZATION FOR SECTIONS OF CODES    *
!                  WHICH CONTAIN CALLS TO INTERNAL PROCEDURES (USEFUL FOR SOME       *
!                  COMPILERS LIKE MIPSpro f90 OR NAG FORTRAN COMPILER WHEN OPENMP    *
!                  IS USED)                                                          *
!  _USE_NAGWARE    FOR DESACTIVATING SOME OPENMP CONSTRUCTS FOR THE NAG FORTRAN      *
!                  COMPILER. HERE, THIS INCLUDES ONLY THE ACTIVATION OF THE          *
!                  _INTERNAL_PROC MACRO                                              *
!                                                                                    *
! ************************************************************************************
!
#ifdef _USE_NAGWARE
#define _INTERNAL_PROC
#endif
!
!
! USED MODULES
! ============
!
    use Select_Parameters,  only : lgl, i4b, stnd
#ifdef _MATMUL
    use Utilities,          only : matmul=>matmul2
#endif
#ifdef _DOT_PRODUCT
    use Utilities,             only : dot_product=>dot_product2
#endif
#ifdef _OPENMP
    use Select_Parameters,  only : omp_limit
    use omp_lib,            only : omp_get_num_procs, omp_get_max_threads, omp_in_parallel
#endif
!
! OTHER MODULES USED IN SPECIFIC SUBROUTINES
! ==========================================
!
!   use Logical_Constants
!   use Char_Constants
!   use Reals_Constants
!   use Num_Constants
!   use Utilities
!   use Random
!   use Stat_Procedures
!   use FFT_Procedures
!   use Prob_Procedures
!
! STRONG TYPING IMPOSED 
! =====================
!    
!   implicit none
!
! PUBLIC ENTITIES 
! ===============
!    
!   ALL SUBROUTINES, FUNCTIONS, VARIABLES AND PARAMETERS ARE PRIVATE BY DEFAULT.
!
  ! private
  !  public ::   comp_smooth, comp_trend, comp_stlez, comp_stl,             &
  !              ma, detrend, hwfilter,  hwfilter2,                         &
  !              lp_coef, lp_coef2, hp_coef, hp_coef2, bd_coef, bd_coef2,   &
  !              pk_coef, moddan_coef, freq_func, extend, taper,            &
  !              symlin_filter, symlin_filter2, moddan_filter, dan_filter,  &
  !              data_window, estim_dof, estim_dof2, comp_conflim,          &
  !              spctrm_ratio, spctrm_ratio2, spctrm_ratio3, spctrm_ratio4, &
  !              spctrm_diff, spctrm_diff2,                                 &
   !             power_spctrm,   cross_spctrm,                              &
  !              power_spctrm2,   cross_spctrm2,                            &
  !              power_spectrum, cross_spectrum,                            &
  !              power_spectrum2, cross_spectrum2
!
! GENERIC INTERFACES FOR ROUTINES WITH OVERLOADED VERSIONS
! ========================================================
!
 !   interface comp_smooth
 !       module procedure    comp_smooth_rv, comp_smooth_rm, comp_smooth_rt
 !   end interface
!
 !   interface comp_trend
 !       module procedure    comp_trend_rv, comp_trend_rm
 !   end interface
!
  !  interface comp_stlez
  !      module procedure    comp_stlez_rv, comp_stlez_rm
  !  end interface
!
  !  interface comp_stl
  !      module procedure    comp_stl_rv, comp_stl_rm
  !  end interface
!
  !  interface hwfilter
  !      module procedure    hwfilter_rv, hwfilter_rm
  !  end interface
!
  !  interface hwfilter2
  !      module procedure    hwfilter2_rv, hwfilter2_rm
  !  end interface
!
  !  interface symlin_filter
  !      module procedure    symlin_filter_rv, symlin_filter_rm
  !  end interface
!
  !  interface symlin_filter2
  !      module procedure    symlin_filter2_rv, symlin_filter2_rm
  !  end interface
!
  !  interface dan_filter
  !      module procedure    dan_filter_rv, dan_filter_rm
  !  end interface
!
 !   interface moddan_filter
 !       module procedure    moddan_filter_rv, moddan_filter_rm
 !   end interface
!
  !  interface extend
  !      module procedure    extend_rv, extend_rm
  !  end interface
!
  !  interface taper
  !      module procedure    taper_rv, taper_rm
  !  end interface
!
  !  interface detrend
  !      module procedure    detrend_rv, detrend_rm
  !  end interface
!
  ! interface comp_conflim
  !      module procedure    comp_conflim_r, comp_conflim_rv
  !  end interface
!
  !  interface spctrm_ratio
  !      module procedure    spctrm_ratio_r, spctrm_ratio_rv
  !  end interface
!
   ! interface spctrm_ratio2
  !      module procedure    spctrm_ratio2_rv, spctrm_ratio2_rm
  !  end interface
!
   ! interface spctrm_ratio3
   !     module procedure    spctrm_ratio3_rv, spctrm_ratio3_rm
   ! end interface
!
    !interface spctrm_ratio4
   !     module procedure    spctrm_ratio4_rv, spctrm_ratio4_rm
   ! end interface
!
   ! interface spctrm_diff
   !     module procedure    spctrm_diff_rv, spctrm_diff_rm
   ! end interface
!
   ! interface spctrm_diff2
   !     module procedure    spctrm_diff2_rv, spctrm_diff2_rm
   ! end interface
!
   ! interface power_spctrm
   !     module procedure    power_spctrm_rv, power_spctrm_rm
   ! end interface
!
   ! interface cross_spctrm
   !     module procedure    cross_spctrm_rv, cross_spctrm_rm
   ! end interface
!
   ! interface power_spctrm2
   !     module procedure    power_spctrm2_rv, power_spctrm2_rm
  !  end interface
!
   ! interface cross_spctrm2
   !     module procedure    cross_spctrm2_rv, cross_spctrm2_rm
   ! end interface
!
  !  interface power_spectrum
  !      module procedure    power_spectrum_rv, power_spectrum_rm
  !  end interface
!
  !  interface cross_spectrum
  !      module procedure    cross_spectrum_rv, cross_spectrum_rm
  !  end interface
!
  !  interface power_spectrum2
  !      module procedure    power_spectrum2_rv, power_spectrum2_rm
  !  end interface
!
 !   interface cross_spectrum2
 !       module procedure    cross_spectrum2_rv, cross_spectrum2_rm
 !   end interface
!
!
! =========================================================================================
!
!                                  contains
!                                 ========
!
! =========================================================================================
!                               SMOOTHING AND TIME SERIES DECOMPOSITIONS
! =========================================================================================
!
!
    subroutine comp_smooth_rv( x, smooth_factor )
       
#if defined (__INTEL_COMPILER) || defined(__ICC)
       !dir$ attributes code_align : 32 :: comp_smooth_rv
       !dir$ attributes forceinline :: comp_smooth_rv
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: comp_smooth_rv
#endif
!
!
! Purpose
! _______
!                                                                              
!   Smooth a time series.
!
!
! Arguments
! _________
!                                                                              
!   X               (INPUT/OUTPUT) real(stnd), dimension(:)
!                   On entry, input vector containing size(X) observations
!                   which must be smoothed with a smoothing factor of SMOOTH_FACTOR.
!                   On exit, the smoothed vector.
!
!
!   SMOOTH_FACTOR   (INPUT) integer(i4b)
!                   On entry, the smoothing factor. The smoothing factor must be greater
!                   than 0 and less than size(X).
!
!
! Further Details
! _______________
!
!   The input vector is smoothed with a moving average of, approximately, (2 * SMOOTH_FACTOR) + 1
!   terms.
!
!   For further details, see:
!
!   (1) Olagnon, M., 1996:
!            Traitement de donnees numeriques avec Fortran 90.
!            Masson, 264 pages, Chapter 11.1.2,  ISBN 2-225-85259-6.
!
!                                                                                                 
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use Reals_Constants,    only : zero, one, half
    use Char_Constants,     only : tseries_error1
    use Utilities,          only : merror
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:), intent(inout) :: x
!
    integer(i4b),             intent(in)    :: smooth_factor
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd), dimension(size(x)) :: xtmp
    !dir$ attributes align : 64 :: xtmp
    real(stnd)                     :: wgt1, wgt2, xi, xj, tmp
!
    integer(i4b) :: nx, i, j
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='comp_smooth'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nx = size( x )
!
!   TEST THE ARGUMENTS.
!
    if ( smooth_factor<=0_i4b ) then
       call merror( name_proc//tseries_error1 )
    end if
!
    if ( nx<=smooth_factor ) return    
!
!   INITIALIZATIONS.
!
    xtmp(1_i4b:nx) = zero
!
    i  = 1_i4b
    j  = nx
    xi = x(i)
    xj = x(j)
!
!   FIRST LOOP TO AVOID SIDE EFFECTS.
!
    !dir$ assume_aligned xtmp:64
    !dir$ assume_aligned x:64
    do
!
        xtmp(i) = xtmp(i)+ xi
        xtmp(j) = xtmp(j)+ xj
!
        i = i + 1_i4b
        j = j - 1_i4b
!
        tmp  = real( i, stnd )
        wgt2 = one/tmp
        wgt1 = ( tmp - one )*wgt2
!
        xi = wgt1*xi + wgt2*x(i)
        xj = wgt1*xj + wgt2*x(j)
!
        if ( i>smooth_factor ) exit
!
    end do        
!
!   MAIN LOOP.
!
    !dir$ assume_aligned xtmp:64
    !dir$ assume_aligned x:64
    do
!
        xtmp(i) = xtmp(i)+ xi
        xtmp(j) = xtmp(j)+ xj
!
        if ( i>=nx ) exit
!
        i = i + 1_i4b
        j = j - 1_i4b
        xi = wgt1*xi + wgt2*x(i)
        xj = wgt1*xj + wgt2*x(j)
!
    end do    
!
    x(1_i4b:nx) = half*xtmp(1_i4b:nx)
!
!
! END OF SUBROUTINE comp_smooth_rv
! ________________________________
!
    end subroutine comp_smooth_rv
!
! ===========================================================================================
!
    subroutine comp_smooth_rm( x, smooth_factor, dimvar )

       !dir$ attributes code_align : 32 :: comp_smooth_rm
       !dir$ attributes forceinline :: comp_smooth_rm
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: comp_smooth_rm

!
! Purpose
! _______
!                                              
!   Smooth the rows or the columns of a matrix.
!
!
! Arguments
! _________
!                                                                  
!   X               (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                   On entry, input matrix containing size(X,3-DIMVAR)
!                   observations on size(X,DIMVAR) variables which must be smoothed with
!                   a smoothing factor of SMOOTH_FACTOR.
!                   By default, DIMVAR is equal to 1. See description of optional DIMVAR
!                   argument for details.
!
!   SMOOTH_FACTOR   (INPUT) integer(i4b)
!                   On entry, the smoothing factor. The smoothing factor must be greater
!                   than 0 and less than size(X,3-DIMVAR).
!
!   DIMVAR          (INPUT, OPTIONAL) integer(i4b)
!                   On entry, if DIMVAR is present, DIMVAR is used as follows:
!
!                   - DIMVAR = 1, the input matrix X contains size(X,2) observations
!                     on size(X,1) variables and the rows of X will be smoothed.
!                   - DIMVAR = 2, the input submatrix X contains size(X,1) observations
!                     on size(X,2) variables and the columns of X will be smoothed.
!
!                   The default is DIMVAR = 1.
!
!
! Further Details
! _______________
!
!   The input matris is smoothed along the specified dimension with a moving average of,
!   approximately, (2 * SMOOTH_FACTOR) + 1 terms.
!
!   For further details, see:
!
!   (1) Olagnon, M., 1996:
!            Traitement de donnees numeriques avec Fortran 90.
!            Masson, 264 pages, Chapter 11.1.2,  ISBN 2-225-85259-6.
!
!                                                                                                 
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use Reals_Constants,   only : zero, one, half
    use Char_Constants,    only : tseries_error1, tseries_error2
    use Utilities,         only : merror
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:,:), intent(inout) :: x
!
    integer(i4b), intent(in)             :: smooth_factor
    integer(i4b), intent(in),  optional  :: dimvar
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd), dimension(size(x,1),size(x,2))       :: xtmp
    !dir$ attributes align : 64 :: xtmp
    real(stnd), dimension(max(size(x,1),size(x,2)))  :: xi, xj
    !dir$ attributes align : 64 :: xi
    !dir$ attributes align : 64 :: xj
    real(stnd)                                       :: wgt1, wgt2, tmp
!
    integer(i4b) :: nx, mx, i, j
!
    logical(lgl) :: dim_is_one
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='comp_smooth'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    if ( present(dimvar) ) then
        if ( dimvar==1_i4b .or. dimvar==2_i4b ) then
            i = dimvar
        else
           call merror( name_proc//tseries_error2 )
        end if
    else
        i = 1_i4b
    end if
!
    dim_is_one = i==1_i4b
!
    if ( smooth_factor<=0_i4b ) then
       call merror( name_proc//tseries_error1 )
    end if
!
    nx = size( x, 1 )
    mx = size( x, 2 )
!
    if ( size( x, 3-int(i) )<=smooth_factor .or. size( x, int(i) )< 1 ) return        
!
!   INITIALIZATIONS.
!
    xtmp(1_i4b:nx,1_i4b:mx) = zero
    i = 1_i4b
!
    if ( dim_is_one ) then
!
!       INITIALIZATIONS.
!
        j = mx
        xi(1_i4b:nx) = x(1_i4b:nx,i)
        xj(1_i4b:nx) = x(1_i4b:nx,j)
!
!       FIRST LOOP TO AVOID SIDE EFFECTS.
!
        !dir$ assume_aligned xtmp:64
        !dir$ assume_aligned xi:64
        !dir$ assume_aligned xj:64
        !dir$ assume_aligned x:64
        do
            xtmp(1_i4b:nx,i) = xtmp(1_i4b:nx,i)+ xi(1_i4b:nx)
            xtmp(1_i4b:nx,j) = xtmp(1_i4b:nx,j)+ xj(1_i4b:nx)
!
            i = i + 1_i4b
            j = j - 1_i4b
!
            tmp = real( i, stnd )
            wgt2 = one/tmp
            wgt1 = (tmp - one)*wgt2
!
            xi(1_i4b:nx) = wgt1*xi(1_i4b:nx) + wgt2*x(1_i4b:nx,i)
            xj(1_i4b:nx) = wgt1*xj(1_i4b:nx) + wgt2*x(1_i4b:nx,j)
!
            if ( i>smooth_factor ) exit
!
        end do        
!
!       MAIN LOOP.
!
        !dir$ assume_aligned xtmp:64
        !dir$ assume_aligned xi:64
        !dir$ assume_aligned xj:64
        !dir$ assume_aligned x:64
        do
!
            xtmp(1_i4b:nx,i) = xtmp(1_i4b:nx,i)+ xi(1_i4b:nx)
            xtmp(1_i4b:nx,j) = xtmp(1_i4b:nx,j)+ xj(1_i4b:nx)
!
            if ( i>=mx ) exit
!
            i = i + 1_i4b
            j = j - 1_i4b
!
            xi(1_i4b:nx) = wgt1*xi(1_i4b:nx) + wgt2*x(1_i4b:nx,i)
            xj(1_i4b:nx) = wgt1*xj(1_i4b:nx) + wgt2*x(1_i4b:nx,j)
!
        end do        
!
    else
!
!       INITIALIZATIONS.
!
        j = nx
!
        xi(1_i4b:mx) = x(i,1_i4b:mx)
        xj(1_i4b:mx) = x(j,1_i4b:mx)
!
!       FIRST LOOP TO AVOID SIDE EFFECTS.
!
        !dir$ assume_aligned xtmp:64
        !dir$ assume_aligned xi:64
        !dir$ assume_aligned xj:64
        !dir$ assume_aligned x:64
        do
!
            xtmp(i,1_i4b:mx) = xtmp(i,1_i4b:mx)+ xi(1_i4b:mx)
            xtmp(j,1_i4b:mx) = xtmp(j,1_i4b:mx)+ xj(1_i4b:mx)
!
            i = i + 1_i4b
            j = j - 1_i4b
!
            tmp = real( i, stnd )
            wgt2 = one/tmp
            wgt1 = (tmp - one)*wgt2
!
            xi(1_i4b:mx) = wgt1*xi(1_i4b:mx) + wgt2*x(i,1_i4b:mx)
            xj(1_i4b:mx) = wgt1*xj(1_i4b:mx) + wgt2*x(j,1_i4b:mx)
!
            if ( i>smooth_factor ) exit
!
        end do        
!
!       MAIN LOOP.
!
        !dir$ assume_aligned xtmp:64
        !dir$ assume_aligned xi:64
        !dir$ assume_aligned xj:64
        !dir$ assume_aligned x:64
        do
            xtmp(i,1_i4b:mx) = xtmp(i,1_i4b:mx)+ xi(1_i4b:mx)
            xtmp(j,1_i4b:mx) = xtmp(j,1_i4b:mx)+ xj(1_i4b:mx)
!
            if ( i>=nx ) exit
!
            i = i + 1_i4b
            j = j - 1_i4b
            xi(1_i4b:mx) = wgt1*xi(1_i4b:mx) + wgt2*x(i,1_i4b:mx)
            xj(1_i4b:mx) = wgt1*xj(1_i4b:mx) + wgt2*x(j,1_i4b:mx)
!
        end do        
!
    end if
!
    x(1_i4b:nx,1_i4b:mx) = half*xtmp(1_i4b:nx,1_i4b:mx)
!
!
! END OF SUBROUTINE comp_smooth_rm
! ________________________________
!
    end subroutine comp_smooth_rm
!
! ===========================================================================================
!
    subroutine comp_smooth_rt( x, smooth_factor )
       !dir$ attributes code_align : 32 :: comp_smooth_rt
       !dir$ attributes forceinline :: comp_smooth_rt
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: comp_smooth_rt
!
! Purpose
! _______
!                                                                              
!   Smooth a tridimensional array along the third dimension.
!
!
! Arguments
! _________
!                                                                              
!   X               (INPUT/OUTPUT) real(stnd), dimension(:,:,:)
!                   On entry, input tridimensional array containing size(X,3)
!                   observations on size(X,1) by size(X,2) variables which must be smoothed with
!                   a smoothing factor of SMOOTH_FACTOR.
!
!   SMOOTH_FACTOR   (INPUT) integer(i4b)
!                   On entry, the smoothing factor. The smoothing factor must be greater
!                   than 0 and less than size(X,3).
!
!
! Further Details
! _______________
!
!   The input tridimensional array is smoothed along the third dimension with a moving average of,
!   approximately, (2 * SMOOTH_FACTOR) + 1 terms.
!
!   For further details, see:
!
!   (1) Olagnon, M., 1996:
!            Traitement de donnees numeriques avec Fortran 90.
!            Masson, 264 pages, Chapter 11.1.2,  ISBN 2-225-85259-6.
!
!                                                                                                 
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use Reals_Constants,   only : zero, one, half
    use Char_Constants,    only : tseries_error1
    use Utilities,         only : merror
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:,:,:), intent(inout) :: x
!
    integer(i4b),                 intent(in)    :: smooth_factor
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd), dimension(size(x,1),size(x,2),size(x,3)) :: xtmp
    !dir$ attributes align : 64 :: xtmp
    real(stnd), dimension(size(x,1),size(x,2))           :: xi, xj
    !dir$ attributes align : 64 :: xi
    !dir$ attributes align : 64 :: xj
    real(stnd)                                           :: wgt1, wgt2, tmp
!
    integer(i4b) :: nx, mx, px, i, j
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='comp_smooth'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    if ( smooth_factor<=0_i4b ) then
        call merror( name_proc//tseries_error1 )
    end if
!
    nx = size( x, 1 )
    mx = size( x, 2 )
    px = size( x, 3 )
!
    if ( px<=smooth_factor .or. min(nx,mx)< 1 ) return        
!
!   INITIALIZATIONS.
!
    xtmp(1_i4b:nx,1_i4b:mx,1_i4b:px) = zero
!
    i = 1_i4b
    j = px
!
    xi(1_i4b:nx,1_i4b:mx) = x(1_i4b:nx,1_i4b:mx,i)
    xj(1_i4b:nx,1_i4b:mx) = x(1_i4b:nx,1_i4b:mx,j)
!
!   FIRST LOOP TO AVOID SIDE EFFECTS.
!
    !dir$ assume_aligned xtmp:64
    !dir$ assume_aligned xi:64
    !dir$ assume_aligned xj:64
    do
        xtmp(1_i4b:nx,1_i4b:mx,i) = xtmp(1_i4b:nx,1_i4b:mx,i)+ xi(1_i4b:nx,1_i4b:mx)
        xtmp(1_i4b:nx,1_i4b:mx,j) = xtmp(1_i4b:nx,1_i4b:mx,j)+ xj(1_i4b:nx,1_i4b:mx)
!
        i = i + 1_i4b
        j = j - 1_i4b
!
        tmp = real( i, stnd )
        wgt2 = one/tmp
        wgt1 = (tmp - one)*wgt2
!
        xi(1_i4b:nx,1_i4b:mx) = wgt1*xi(1_i4b:nx,1_i4b:mx) + wgt2*x(1_i4b:nx,1_i4b:mx,i)
        xj(1_i4b:nx,1_i4b:mx) = wgt1*xj(1_i4b:nx,1_i4b:mx) + wgt2*x(1_i4b:nx,1_i4b:mx,j)
!
        if ( i>smooth_factor ) exit
!
    end do        
!
!   MAIN LOOP.
!
    !dir$ assume_aligned xtmp:64
    !dir$ assume_aligned xi:64
    !dir$ assume_aligned xj:64
    do
!
        xtmp(1_i4b:nx,1_i4b:mx,i) = xtmp(1_i4b:nx,1_i4b:mx,i)+ xi(1_i4b:nx,1_i4b:mx)
        xtmp(1_i4b:nx,1_i4b:mx,j) = xtmp(1_i4b:nx,1_i4b:mx,j)+ xj(1_i4b:nx,1_i4b:mx)
!
        if ( i>=px ) exit
!
        i = i + 1_i4b
        j = j - 1_i4b
!
        xi(1_i4b:nx,1_i4b:mx) = wgt1*xi(1_i4b:nx,1_i4b:mx) + wgt2*x(1_i4b:nx,1_i4b:mx,i)
        xj(1_i4b:nx,1_i4b:mx) = wgt1*xj(1_i4b:nx,1_i4b:mx) + wgt2*x(1_i4b:nx,1_i4b:mx,j)
!
    end do        
!
    x(1_i4b:nx,1_i4b:mx,1_i4b:px) = half*xtmp(1_i4b:nx,1_i4b:mx,1_i4b:px)
!
!
! END OF SUBROUTINE comp_smooth_rt
! ________________________________
!
    end subroutine comp_smooth_rt
!
! =========================================================================================
!
    subroutine comp_trend_rv( y, nt, itdeg, robust, trend, ntjump, maxiter, rw, no, ok )
       !dir$ attributes code_align : 32 :: comp_trend_rv
       !dir$ attributes forceinline :: comp_trend_rv
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: comp_trend_rv
!
! Purpose
! _______
!                                                                              
!   COMP_TREND extracts a smoothed component from a time series using a LOESS method.
!   It returns the smoothed component (e.g. the trend) and, optionally, the robustness
!   weights.
!
!
! Arguments
! _________
!                                                                              
!   Y        (INPUT) real(stnd), dimension(:)
!            On entry, the time series to be decomposed.
!
!   NT       (INPUT) integer(i4b)
!            On entry, the length of the trend smoother. The value of
!            NT should be an odd integer greater than or equal to 3.
!            As NT increases the values of the trend component become smoother.
!
!   ITDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted  polynomial  in  trend
!            smoothing.  The value must be 0, 1 or 2.
!
!   ROBUST   (INPUT) logical(lgl)
!            On entry, TRUE if robustness iterations are to be used,
!            FALSE otherwise.
!
!            Robustness iterations are  carried out until convergence of
!            the trend component, with MAXITER iterations maximum.
!            Convergence occurs if the maximum changes in the
!            trend fit is less than 1% of the component's range
!            after the previous iteration.
!
!   TREND    (OUTPUT) real(stnd), dimension(:)
!            On output, the smoothed (e.g. trend) component.
!
!            TREND must verify:  size(TREND) = size(Y).
!
!   NTJUMP   (INPUT, OPTIONAL) integer(i4b)
!            On entry, the  skipping value for trend smoothing.
!            By default, NTJUMP is set to NT/10.
!
!   MAXITER  (INPUT, OPTIONAL) integer(i4b)
!            On entry, the maximum number of robustness iterations.
!
!            The default is 15. This argument is not used if ROBUST=FALSE.
!
!   RW       (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!            On output, final  robustness  weights. All RW elements are 1 if
!            ROBUST=FALSE .
!
!            RW must verify:  size(RW) = size(Y).
!
!   NO       (OUTPUT, OPTIONAL) integer(i4b)
!            On output, if:
!
!            - ROBUST=TRUE : the number of robustness iterations.  The  iterations
!              end if a convergence criterion is met or if the number is MAXITER.
!            - ROBUST=FALSE : NO is set to 0.
!
!   OK       (OUTPUT, OPTIONAL) logical(lgl)
!            On output, if:
!
!            - ROBUST=TRUE : OK is set to TRUE if the convergence criterion is met
!              and to FALSE otherwise.
!            - ROBUST=FALSE : OK is set to TRUE.
!
!
! Further Details
! _______________
!
!   This subroutine is adapted from subroutine STL developped by
!   Cleveland and coworkers at AT&T Bell Laboratories.
!
!   This subroutine decomposes a time series into trend and residual components,
!   assuming that the time series has no seasonal cycle or other harmonic components.
!   The algorithm uses LOESS interpolation to smooth the time series and find the trend.
!
!   The LOESS smoother for estimating the trend is specified with three parameters: a width
!   (e.g. NT), a degree (e.g. ITDEG) and a jump (e.g. NTJUMP). The width specifies the number
!   of data points that the local interpolation uses to smooth each point, the degree specifies
!   the degree of the local polynomial that is fit to the data, and the jump specifies how many
!   points are skipped between Loess interpolations, with linear interpolation being done between
!   these points.
!
!   If the optional ROBUST argument is set to true, the process is iterative and includes robustness
!   iterations that take advandages of the weighted-least-squares underpinnings of LOESS to remove the
!   effects of outliers.
!
!   Note that, finally, that this subroutine expects equally spaced data with no missing values.
!
!   For further details, see:
!   
!   (1)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I.: 
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           Statistics Research Report, AT&T Bell Laboratories.
!   
!   (2)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I., 1990:
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           J. Official Stat., 6, 3-73.
!
!   (3)  Crotinger, J., 2017:
!            Java implementation of Seasonal-Trend-Loess time-series decomposition algorithm.
!            https://github.com/ServiceNow/stl-decomp-4j
!                                                                                                 
!                                                                                                 
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use Logical_Constants, only : true, false
    use Reals_Constants,   only : one, c0_9, c1_m2, ten
    use Char_Constants,    only : tseries_error3
    use Utilities,         only : assert, assert_eq, merror
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:), intent(in)             :: y
    real(stnd), dimension(:), intent(out)            :: trend
    real(stnd), dimension(:), intent(out), optional  :: rw
!
    integer(i4b), intent(in)             :: nt, itdeg
    integer(i4b), intent(in),   optional :: ntjump, maxiter
    integer(i4b), intent(out),  optional :: no
!
    logical(lgl), intent(in)             :: robust
    logical(lgl), intent(out),  optional :: ok
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd), dimension(size(y))     :: trend2, rw2, work
    !dir$ attributes align : 64 :: trend2
    !dir$ attributes align : 64 :: rw2
    !dir$ attributes align : 64 :: work
    real(stnd)                         :: maxt, mint, maxdt
!
    integer(i4b) :: n, k, nt2, ntjump2, maxiter2
!
    logical(lgl) :: ok2
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='comp_trend'
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = assert_eq(  int(size(y),i4b),        &
                    int(size(trend),i4b),    &
                    name_proc )
!
    if ( present(rw) ) then
        call assert( logical( n==int(size(rw),i4b), lgl ), name_proc )
    end if
!
    nt2 = max( 3_i4b, nt )
    if ( mod(nt2,2_i4b)==0_i4b) nt2 = nt2 + 1_i4b
!
    if ( itdeg<0_i4b .or. itdeg>2_i4b  ) then
        call merror( name_proc//tseries_error3 )
    end if
!
    if ( present(ntjump) ) then
        ntjump2 = min( max( 1_i4b, ntjump ), n-1_i4b )
    else
        ntjump2 = max( 1_i4b, int(real(nt2,stnd)/ten + c0_9,i4b) )
    end if
!
    if ( present(maxiter) ) then
        maxiter2 = max( maxiter, 1_i4b )
    else
        maxiter2 = 15_i4b
    end if
!
!   INITIALIZATIONS.
!
    rw2(1_i4b:n)  = one
    k             = 0_i4b
    ok2           = true
!
!   FIRST ITERATION.
!
    work(1_i4b:n) = y(1_i4b:n)
!
    call ess( work(1_i4b:n), nt2, itdeg, ntjump2, false, rw2(1_i4b:n), trend2(1_i4b:n) )
!
!   ROBUST ITERATIONS.
!
    if ( robust ) then
!
        !dir$ assume_aligned trend:64
        !dir$ assume_aligned work:64
        !dir$ assume_aligned rw2:64
        !dir$ assume_aligned trend2:64
        do
!
            trend(1_i4b:n)  = trend2(1_i4b:n)
!
            k = k + 1_i4b
!
!           COMPUTE ROBUST WEIGHTS.
!
            call rwts( work(1_i4b:n), trend2(1_i4b:n), rw2(1_i4b:n) )
!
            call ess( work(1_i4b:n), nt2, itdeg, ntjump2, true, rw2(1_i4b:n), trend2(1_i4b:n) )
!
!           TEST CONVERGENCE OF THE ROBUST ITERATIONS.
!
            maxt  = maxval( trend(1_i4b:n) )
            mint  = minval( trend(1_i4b:n) )
            maxdt = maxval( abs( trend2(1_i4b:n) - trend(1_i4b:n) ) )
!
            ok2 =(maxdt/(maxt-mint))<c1_m2
!
            if ( ok2 .or. k>=maxiter2 ) exit
!
        end do
!
    end if
!
!   SAVE FINAL RESULTS.
!
    trend(1_i4b:n)  = trend2(1_i4b:n)
!
    if ( present(rw) ) then
        rw(1_i4b:n) = rw2(1_i4b:n)
    end if
!
    if ( present(no) ) then
        no = k
    end if
!
    if ( present(ok) ) then
        ok = ok2
    end if
!
!
! END OF SUBROUTINE comp_trend_rv
! _______________________________
!
    end subroutine comp_trend_rv
!
! =========================================================================================
!
    subroutine comp_trend_rm( y, nt, itdeg, robust, trend, ntjump, maxiter, rw,  no, ok )
       !dir$ attributes code_align : 32 :: comp_trend_rm
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: comp_trend_rm
!
!
! Purpose
! _______
!                                                                              
!   COMP_TREND extracts smoothed components from the (time series) columns of a matrix
!   using a LOESS method. It returns the smoothed components (e.g. the trends) and,
!   optionally, the robustness weights.
!
!
! Arguments
! _________
!                                                                              
!   Y        (INPUT) real(stnd), dimension(:,:)
!            On entry, the matrix to be decomposed.
!
!   NT       (INPUT) integer(i4b)
!            On entry, the length of the trend smoother. The value of
!            NT should be an odd integer greater than or equal to 3.
!            As NT increases the values of the trend component become smoother.
!
!   ITDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted  polynomial  in  trend
!            smoothing.  The value must be 0, 1 or 2.
!
!   ROBUST   (INPUT) logical(lgl)
!            On entry, TRUE if robustness iterations are to be used,
!            FALSE otherwise.
!
!            Robustness iterations are  carried out until convergence of
!            the trend component, with MAXITER iterations maximum.
!            Convergence occurs if the maximum changes in the
!            trend fit is less than 1% of the component's range
!            after the previous iteration.
!
!   TREND    (OUTPUT) real(stnd), dimension(:,:)
!            On output, the smoothed (e.g. trend) components.
!
!            TREND must verify: size(TREND,1) = size(Y,1) and size(TREND,2) = size(Y,2).
!
!   NTJUMP   (INPUT, OPTIONAL) integer(i4b)
!            On entry, the  skipping value for trend smoothing.
!
!            By default, NTJUMP is set to NT/10.
!
!   MAXITER  (INPUT, OPTIONAL) integer(i4b)
!            On entry, the maximum number of robustness iterations.
!
!            The default is 15. This argument is not used if ROBUST=FALSE.
!
!   RW       (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!            On output, final  robustness  weights. All RW elements are 1 if
!            ROBUST=FALSE .
!
!            RW must verify: size(RW,1) = size(Y,1) and size(RW,2) = size(Y,2).
!
!   NO       (OUTPUT, OPTIONAL) integer(i4b), dimension(:)
!            On output, if
!
!            - ROBUST=TRUE : NO(i) is the number of robustness iterations for each
!              time series Y(:,i). The  iterations end if a convergence criterion is met 
!              or if the number is MAXITER for each time series Y(:,i).
!            - ROBUST=FALSE : NO(:) is set to 0.
!
!            NO must verify:  size(NO) = size(Y,2).
!
!   OK       (OUTPUT, OPTIONAL) logical(lgl), dimension(:)
!            On output, if
!
!            - ROBUST=TRUE : OK(i) is set to TRUE if the convergence criterion is met
!              for time series Y(:,i) and to FALSE otherwise.
!            - ROBUST=FALSE : OK(:) is set to TRUE.
!
!            OK must verify:  size(OK) = size(Y,2).
!
!
! Further Details
! _______________
!
!   This subroutine is adapted from subroutine STL developped by
!   Cleveland and coworkers at AT&T Bell Laboratories.
!
!   This subroutine decomposes a multi-channel time series into trend and residual components,
!   assuming that the multi-channel time series has no seasonal cycle or other harmonic components.
!   The algorithm uses LOESS interpolation to smooth the multi-channel time series and find the trends.
!
!   The LOESS smoother for estimating the trends is specified with three parameters: a width
!   (e.g. NT), a degree (e.g. ITDEG) and a jump (e.g. NTJUMP). The width specifies the number
!   of data points that the local interpolation uses to smooth each point, the degree specifies
!   the degree of the local polynomial that is fit to the data, and the jump specifies how many
!   points are skipped between Loess interpolations, with linear interpolation being done between
!   these points.
!
!   If the optional ROBUST argument is set to true, the process is iterative and includes robustness
!   iterations that take advandages of the weighted-least-squares underpinnings of LOESS to remove the
!   effects of outliers.
!
!   Note that, finally, that this subroutine expects equally spaced data with no missing values.
!
!   For further details, see description of COMP_STL and :
!   
!   (1)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I.: 
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           Statistics Research Report, AT&T Bell Laboratories.
!   
!   (2)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I., 1990:
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           J. Official Stat., 6, 3-73.
!
!   (3)  Crotinger, J., 2017:
!            Java implementation of Seasonal-Trend-Loess time-series decomposition algorithm.
!            https://github.com/ServiceNow/stl-decomp-4j
!                                                                                                 
!                                                                                                 
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    
    use Select_Parameters,  only : lgl, i4b, stnd
    use Select_Parameters,  only : omp_limit
    use omp_lib,            only : omp_get_num_procs, omp_get_max_threads, omp_in_parallel
    use Logical_Constants, only : true, false
    use Reals_Constants,   only : one, c0_9, c1_m2, ten
    use Char_Constants,    only : tseries_error3
    use Utilities,         only : assert, assert_eq, merror
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
 
    real(stnd), dimension(:,:), intent(in)             :: y
    real(stnd), dimension(:,:), intent(out)            :: trend
    real(stnd), dimension(:,:), intent(out),  optional :: rw
!
    integer(i4b),               intent(in)             :: nt, itdeg
    integer(i4b),               intent(in),   optional :: ntjump, maxiter
    integer(i4b), dimension(:), intent(out),  optional :: no
!
    logical(lgl),               intent(in)             :: robust
    logical(lgl), dimension(:), intent(out),  optional :: ok
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd), dimension(size(y,1)) :: trend2, rw2, work
    !dir$ attributes align : 64 :: trend2
    !dir$ attributes align : 64 :: rw2
    !dir$ attributes align : 64 :: work
    real(stnd)                       :: maxt, mint, maxdt
!
    integer(i4b) :: n, m, i, k, nt2, ntjump2, maxiter2
!
    logical(lgl) :: ok2, use_no, use_ok, use_rw

    logical      :: test_par

!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='comp_trend'
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = assert_eq(  int(size(y,1),i4b),        &
                    int(size(trend,1),i4b),    &
                    name_proc )
!
    m = assert_eq(  int(size(y,2),i4b),        &
                    int(size(trend,2),i4b),    &
                    name_proc )
!
    use_rw = false
!
    if ( present(rw) ) then
        call assert( logical( n==int(size(rw,1),i4b), lgl ),     &
                     logical( m==int(size(rw,2),i4b), lgl ),     &
                     name_proc )
        use_rw = true
    end if
!
    use_no = false
!
    if ( present(no) ) then
        call assert( logical( m==int(size(no),i4b), lgl ),  name_proc )
        use_no = true
    end if
!
    use_ok = false
!
    if ( present(ok) ) then
        call assert( logical( m==int(size(ok),i4b), lgl ),  name_proc )
        use_ok = true
    end if
!
    nt2 = max( 3_i4b, nt )
    if ( mod(nt2,2_i4b)==0_i4b) nt2 = nt2 + 1_i4b
!
    if ( itdeg<0_i4b .or. itdeg>2_i4b  ) then
       call merror( name_proc//tseries_error3 )
    end if
!
    if ( present(ntjump) ) then
        ntjump2 = min( max( 1_i4b, ntjump ), n-1_i4b )
    else
        ntjump2 = max( 1_i4b, int(real(nt2,stnd)/ten + c0_9,i4b) )
    end if
!
    if ( present(maxiter) ) then
        maxiter2 = max( maxiter, 1_i4b )
    else
        maxiter2 = 15_i4b
    end if
!

    i = omp_get_num_procs()
    k = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )      .and.      &
               (m*n)>=omp_limit                .and.      &
               i>1_i4b                         .and.      &
               k>1_i4b                         .and.      &
               m>1_i4b

     
!
!$OMP PARALLEL IF(test_par)                                      &
!$OMP         ,PRIVATE(i,k,ok2,rw2,trend2,work,maxt,mint,maxdt)  &
!$OMP         ,SHARED(nt2,itdeg,ntjump2)
!
!   INITIALIZATIONS.
!
    rw2(1_i4b:n)  = one
    ok2           = true
!
!$OMP DO SCHEDULE(STATIC,8) 
!
    !dir$ assume_aligned work:64
    !dir$ assume_aligned y:64
    !dir$ assume_aligned rw:64
    !dir$ assume_aligned trend2:64
    !dir$ assume_aligned trend:64
    do i = 1_i4b, m
!
        k = 0_i4b
!
!       FIRST ITERATION.
!
        work(1_i4b:n) = y(1_i4b:n,i)
!
        call ess( work(1_i4b:n), nt2, itdeg, ntjump2, false, rw2(1_i4b:n), trend2(1_i4b:n) )
!
!       ROBUST ITERATIONS.
!
        if ( robust ) then
!
            do
!
                trend(1_i4b:n,i)  = trend2(1_i4b:n)
!
                k = k + 1_i4b
!
!               COMPUTE ROBUST WEIGHTS.
!
                call rwts( work(1_i4b:n), trend2(1_i4b:n), rw2(1_i4b:n) )
!
                call ess( work(1_i4b:n), nt2, itdeg, ntjump2, true, rw2(1_i4b:n), trend2(1_i4b:n) )
!
!               TEST CONVERGENCE OF THE ROBUST ITERATIONS.
!
                maxt  = maxval( trend(1_i4b:n,i) )
                mint  = minval( trend(1_i4b:n,i) )
                maxdt = maxval( abs( trend2(1_i4b:n) - trend(1_i4b:n,i) ) )
!
                ok2 = (maxdt/(maxt-mint))<c1_m2
!
                if ( ok2 .or. k>=maxiter2 ) exit
!
            end do
!
        end if
!
!       SAVE FINAL RESULTS.
!
        trend(1_i4b:n,i)  = trend2(1_i4b:n)
!
        if ( use_rw ) then
            rw(1_i4b:n,i) = rw2(1_i4b:n)
        end if
!
        if ( use_no ) then
            no(i) = k
        end if
!
        if ( use_ok) then
            ok(i) = ok2
        end if
!
    end do
!
!$OMP END DO
!
!$OMP END PARALLEL
!
!
! END OF SUBROUTINE comp_trend_rm
! _______________________________
!
    end subroutine comp_trend_rm
!
! =========================================================================================
!
    subroutine comp_stlez_rv( y, np, ns, isdeg, itdeg, robust, season, trend, ni, nt, nl, &
                              ildeg, nsjump, ntjump, nljump, maxiter, rw, no, ok )
       !dir$ attributes code_align : 32 :: comp_stlez_rv
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: comp_stlez_rv
!
!
! Purpose
! _______
!                                              
!   COMP_STLEZ decomposes a time series into seasonal and trend  components.
!   It returns the components and, optionally, the robustness weights.
!
!   COMP_STLEZ offers an easy to use version of COMP_STL subroutine, also
!   included in STATPACK, by defaulting most parameters values associated
!   with the three LOESS smoothers used in COMP_STL.
!
!
! Arguments
! _________
!                                          
!   Y        (INPUT) real(stnd), dimension(:)
!            On entry, the time series to be decomposed.
!
!   NP       (INPUT) integer(i4b)
!            On entry, the period of the seasonal component. For example,
!            if  the  time series is monthly with a yearly cycle, then
!            NP=12 should be used. NP must be greater than 1.
!
!   NS       (INPUT) integer(i4b)
!            On entry, the length of the seasonal smoother.  The value of
!            NS should be an odd integer greater than or equal to 3; NS>6
!            is recommended.   As NS  increases  the  values  of  the
!            seasonal component at a given point in the seasonal cycle
!            (e.g., January values of a monthly series with  a  yearly
!            cycle) become smoother.
!
!   ISDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted  polynomial in seasonal
!            smoothing.  The value must be 0 or 1.
!
!   ITDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted  polynomial in trend
!            smoothing.  The value must be 0, 1 or 2.
!
!   ROBUST   (INPUT) logical(lgl)
!            On entry, TRUE if robustness iterations are to be used,
!            FALSE otherwise.
!
!            Robustness iterations are  carried out until convergence of
!            both seasonal and trend components, with MAXITER iterations
!            maximum.  Convergence occurs if the maximum changes in individual
!            seasonal  and  trend fits are less than 1% of the component's range
!            after the previous iteration.
!
!   SEASON   (OUTPUT) real(stnd), dimension(:)
!            On output, the seasonal component.
!
!            SEASON must verify:  size(SEASON) = size(Y).
!
!   TREND    (OUTPUT) real(stnd), dimension(:)
!            On output, the trend component.
!
!            TREND must verify:  size(TREND) = size(Y).
!
!   NI       (INPUT, OPTIONAL) integer(i4b)
!            On entry, the number of loops for updating the seasonal and trend components.
!            The value of NI  should be a positive integer.
!
!            By default, NI=2 if ROBUST=FALSE and NI=1 if ROBUST=TRUE.
!   
!   NT       (INPUT, OPTIONAL) integer(i4b)
!            On entry, the length of the trend smoother. The value of
!            NT should be an odd integer greater than or equal to 3.
!            A value of NT between 1.5 * NP and 2 * NP is recommended.
!            As NT increases the values of the trend component become smoother.
!
!            By default, NT is set to the smallest odd integer greater
!            than or equal to  (1.5 * NP) / (1-(1.5/NS)).
!
!   NL       (INPUT, OPTIONAL) integer(i4b)
!            On entry, the length of the low-pass filter. The value 
!            of NL should be an odd integer greater than or equal to 3.
!            The smallest odd integer greater than or equal to NP is recommended.
!
!            By default, NL is set to the smallest odd integer greater
!            than or equal to NP.
!
!   ILDEG    (INPUT, OPTIONAL) integer(i4b)
!            On entry, the degree of locally-fitted polynomial in low-pass smoothing.
!
!            By default, ILDEG is set to ITDEG.
!
!   NSJUMP   (INPUT, OPTIONAL) integer(i4b)
!            On entry, the skipping value for seasonal smoothing.
!            The seasonal  smoother  skips  ahead NSJUMP points and then
!            linearly interpolates in between.  The  value of NSJUMP
!            should  be  a  positive  integer; if NSJUMP=1, a seasonal
!            smooth is calculated  at all size(Y) points. To  make  the
!            procedure  run  faster, a reasonable choice for NSJUMP is
!            10% or 20% of NS.
!
!            By default, NSJUMP is set to NS/10.
!
!   NTJUMP   (INPUT, OPTIONAL) integer(i4b)
!            On entry, the  skipping value for trend smoothing.
!
!            By default, NTJUMP is set to NT/10.
!
!   NLJUMP   (INPUT, OPTIONAL) integer(i4b)
!            On entry, the skipping value for the low-pass filter.
!
!            By default, NLJUMP is set to NL/10.
!
!   MAXITER  (INPUT, OPTIONAL) integer(i4b)
!            On entry, the maximum number of robustness iterations.
!
!            The default is 15. This argument is not used if ROBUST=FALSE.
!
!   RW       (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!            On output, final  robustness  weights. All RW elements are 1 if
!            ROBUST=FALSE .
!
!            RW must verify:  size(RW) = size(Y).
!
!   NO       (OUTPUT, OPTIONAL) integer(i4b)
!            On output, if
!
!            - ROBUST=TRUE : the number of robustness iterations.  The  iterations
!              end if a convergence criterion is met or if the number is MAXITER.
!            - ROBUST=FALSE : NO is set to 0.
!
!   OK       (OUTPUT, OPTIONAL) logical(lgl)
!            On output, if
!
!            - ROBUST=TRUE : OK is set to TRUE if the convergence criterion is met
!              and to FALSE otherwise.
!            - ROBUST=FALSE : OK is set to TRUE.
!
!
! Further Details
! _______________
!
!   This subroutine is a FORTRAN90 implementation of subroutine STLEZ developped by
!   Cleveland and coworkers at AT&T Bell Laboratories.
!
!   At a minimum, COMP_STLEZ requires specifying the periodicity of the data
!   (e.g. NP, 12 for monthly), the width of the LOESS smoother used to smooth the
!   cyclic seasonal sub-series (e.g. NS) and the degree of the locally-fitted polynomial
!   in seasonal (e.g. ISDEG) and trend (e.g. ITDEG) smoothing.
!
!   COMP_STLEZ sets, by default, others parameters of the STL procedure to the values recommended
!   in Cleveland et al. (1990). It also includes tests of convergence if robust iterations
!   are carried out. Otherwise, COMP_STLEZ is similar to COMP_STL.
!
!   For further details, see description of COMP_STL and:
!   
!   (1)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I.: 
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           Statistics Research Report, AT&T Bell Laboratories.
!   
!   (2)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I., 1990:
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           J. Official Stat., 6, 3-73.
!
!   (3)  Crotinger, J., 2017:
!            Java implementation of Seasonal-Trend-Loess time-series decomposition algorithm.
!            https://github.com/ServiceNow/stl-decomp-4j
!
!
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use Logical_Constants, only : true, false
    use Reals_Constants,   only : zero, one, half, c1_5, c0_9, c1_m2, ten
    use Char_Constants,    only : tseries_error3, tseries_error4, tseries_error5, tseries_error6
    use Utilities,         only : assert, assert_eq, merror
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:), intent(in)             :: y
    real(stnd), dimension(:), intent(out)            :: season, trend
    real(stnd), dimension(:), intent(out), optional  :: rw
!
    integer(i4b), intent(in)             :: np, ns, isdeg, itdeg
    integer(i4b), intent(in),   optional :: ni, nt, nl, ildeg, nsjump, ntjump, nljump, maxiter
    integer(i4b), intent(out),  optional :: no
!
    logical(lgl), intent(in)             :: robust
    logical(lgl), intent(out),  optional :: ok
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd), dimension(size(y)+2*max(np,1_i4b),3) :: work
    !dir$ attributes align : 64 ::  work
    real(stnd), dimension(size(y))                   :: trend2, season2, rw2
    !dir$ attributes align : 64 :: trend2
    !dir$ attributes align : 64 :: season2
    !dir$ attributes align : 64 :: rw2
    real(stnd)                                       :: maxs, mins, maxt, mint, maxds, maxdt
!
    integer(i4b) :: n, j, k, nt2, nl2, ni2, ns2, ildeg2, nsjump2, ntjump2, nljump2, maxiter2
!
    logical(lgl) :: ok2
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='comp_stlez'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = assert_eq(  int(size(y),i4b),        &
                    int(size(season),i4b),   &
                    int(size(trend),i4b),    &
                    name_proc )
!
    if ( present(rw) ) then
        call assert( logical( n==int(size(rw),i4b), lgl ), name_proc )
    end if
!
    if ( np<=1_i4b ) then
        call merror( name_proc//tseries_error4 )
    end if
!
    if ( isdeg/=1_i4b .and. isdeg/=0_i4b  ) then
        call merror( name_proc//tseries_error5 )
    end if
!
    if ( itdeg<0_i4b .or. itdeg>2_i4b  ) then
        call merror( name_proc//tseries_error3 )
    end if
!
    if ( present(ildeg) ) then
        if ( ildeg<0_i4b .or. ildeg>2_i4b  ) then
            call merror( name_proc//tseries_error6 )
        end if
        ildeg2 = ildeg
    else
        ildeg2 = itdeg
    end if
!
    ns2 = max( 3_i4b, ns )
    if ( mod(ns2,2_i4b)==0_i4b) ns2 = ns2 + 1_i4b
!
    if ( present(nt) ) then
        nt2 = nt
    else
        nt2 = ( c1_5*real(np,stnd) )/( one - c1_5/real(ns2,stnd) ) + half
    end if
!
    nt2 = max( 3_i4b, nt2 )
!
    if ( mod(nt2,2_i4b)==0_i4b) nt2 = nt2 + 1_i4b
!
    if ( present(nl) ) then
        nl2 = max( 3_i4b, nl )
    else
        nl2 = np
    end if
    if ( mod(nl2,2_i4b)==0_i4b) nl2 = nl2 + 1_i4b
!
    if ( present(ni) ) then
        ni2 = max( ni, 1_i4b )
    else
        if ( robust ) then
            ni2 = 1_i4b
        else
            ni2 = 2_i4b
        end if
    end if
!
    if ( present(nsjump) ) then
        nsjump2 = min( max( 1_i4b, nsjump ), n-1_i4b )
    else
        nsjump2 = max( 1_i4b, int(real(ns2,stnd)/ten + c0_9,i4b) )
    end if
!
    if ( present(ntjump) ) then
        ntjump2 = min( max( 1_i4b, ntjump ), n-1_i4b )
    else
        ntjump2 = max( 1_i4b, int(real(nt2,stnd)/ten + c0_9,i4b) )
    end if
!
    if ( present(nljump) ) then
        nljump2 = min( max( 1_i4b, nljump ), n-1_i4b )
    else
        nljump2 = max( 1_i4b, int(real(nl2,stnd)/ten + c0_9,i4b) )
    end if
!
    if ( present(maxiter) ) then
        maxiter2 = max( maxiter, 1_i4b )
    else
        maxiter2 = 15_i4b
    end if
!
!   INITIALIZATIONS.
!
    trend2(1_i4b:n) = zero
    rw2(1_i4b:n)    = one
    k               = 0_i4b
    ok2             = true
!
!   FIRST ITERATIONS.
!
     !dir$ assume_aligned work:64
     !dir$ assume_aligned y:64
     !dir$ assume_aligned trend2:64
     !dir$ assume_aligned rw2:64
     !dir$ assume_aligned season2:64
    do j = 1_i4b, ni2
!
!       STEP1 : DETRENDING.
!
        work(1_i4b:n,1_i4b) = y(1_i4b:n) - trend2(1_i4b:n)
!
!       STEP2 : CYCLE-SUBSERIES SMOOTHING.
!
        call ss( work(1_i4b:n,1_i4b), np, ns2, isdeg, nsjump2, false,     &
                 rw2(1_i4b:n), work(1_i4b:n+2_i4b*np,2_i4b) )
!
!       call fts( work(:n+2*np,2), np, work(:n+2*np,3), work(:n+2*np,1) )
!
!       STEP3 : LOW-PASS FILTERING OF SMOOTHED CYCLE-SUBSERIES.
!
        call ma( work(1_i4b:n+2_i4b*np,2_i4b), np, work(1_i4b:n+2_i4b*np,3_i4b) )
        call ma( work(1_i4b:n+np+1_i4b,3_i4b), np, work(1_i4b:n+np+1_i4b,1_i4b) )
        call ma( work(1_i4b:n+2_i4b,1_i4b), 3_i4b, work(1_i4b:n+2_i4b,3_i4b) )

        call ess( work(1_i4b:n,3_i4b), nl2, ildeg2, nljump2, false,        &
                  rw2(1_i4b:n), work(1_i4b:n,1_i4b) )
!
!       STEP4 : DETRENDING OF SMOOTHED CYCLE-SUBSERIES.
!
        season2(1_i4b:n) = work(np+1_i4b:n+np,2_i4b) - work(1_i4b:n,1_i4b)
!
!       STEP5 : DESEASONALIZING.
!
        work(1_i4b:n,1_i4b) = y(1_i4b:n) - season2(1_i4b:n)
!
!       STEP6 : TREND SMOOTHING.
!
        call ess( work(1_i4b:n,1_i4b), nt2, itdeg, ntjump2, false, rw2(1_i4b:n), trend2(1_i4b:n) )
!
    end do
!
!   ROBUST ITERATIONS.
!
    if ( robust ) then
!
     !dir$ assume_aligned work:64
     !dir$ assume_aligned y:64
     !dir$ assume_aligned trend2:64
     !dir$ assume_aligned rw2:64
     !dir$ assume_aligned season2:64
        do
!
!           UPDATE SEASONAL AND TREND COMPONENTS.
!
            season(1_i4b:n) = season2(1_i4b:n)
            trend(1_i4b:n)  = trend2(1_i4b:n)
!
            k = k + 1_i4b
!
!           COMPUTE ROBUST WEIGHTS.
!
            work(1_i4b:n,1_i4b) = trend2(1_i4b:n) + season2(1_i4b:n)
            call rwts( y(1_i4b:n), work(1_i4b:n,1_i4b), rw2(1_i4b:n) )
!
!           INNER LOOP.
!
            do j = 1_i4b, ni2
!
!               STEP1 : DETRENDING STEP.
!
                work(1_i4b:n,1_i4b) = y(1_i4b:n) - trend2(1_i4b:n)
!
!               STEP2 : CYCLE-SUBSERIES SMOOTHING.
!
                call ss( work(1_i4b:n,1_i4b), np, ns2, isdeg, nsjump2, true,    &
                         rw2(1_i4b:n), work(1_i4b:n+2_i4b*np,2_i4b) )
!
!               call fts( work(:n+2*np,2), np, work(:n+2*np,3), work(:n+2*np,1) )
!
!               STEP3 : LOW-PASS FILTERING OF SMOOTHED CYCLE-SUBSERIES.
!
                call ma( work(1_i4b:n+2_i4b*np,2_i4b), np, work(1_i4b:n+2_i4b*np,3_i4b) )
                call ma( work(1_i4b:n+np+1,3_i4b), np, work(1_i4b:n+np+1,1_i4b) )
                call ma( work(1_i4b:n+2_i4b,1_i4b), 3_i4b, work(1_i4b:n+2_i4b,3_i4b) )

                call ess( work(1_i4b:n,3_i4b), nl2, ildeg2, nljump2, false,       &
                          rw2(1_i4b:n), work(1_i4b:n,1_i4b) )
!
!               STEP4 : DETRENDING OF SMOOTHED CYCLE-SUBSERIES.
!
                season2(1_i4b:n) = work(np+1_i4b:n+np,2_i4b) - work(1_i4b:n,1_i4b)
!
!               STEP5 : DESEASONALIZING.
!
                work(1_i4b:n,1_i4b) = y(1_i4b:n) - season2(1_i4b:n)
!
!               STEP6 : TREND SMOOTHING.
!
                call ess( work(1_i4b:n,1_i4b), nt2, itdeg, ntjump2, true,            &
                          rw2(1_i4b:n), trend2(1_i4b:n) )
!
            end do
!
!           TEST CONVERGENCE OF THE ROBUST ITERATIONS.
!
            maxs  = maxval( season(1_i4b:n) )
            mins  = minval( season(1_i4b:n) )
            maxt  = maxval( trend(1_i4b:n) )
            mint  = minval( trend(1_i4b:n) )
            maxds = maxval( abs( season2(1_i4b:n) - season(1_i4b:n) ) )
            maxdt = maxval( abs( trend2(1_i4b:n) - trend(1_i4b:n) ) )
!
            ok2 = (maxds/(maxs-mins)<c1_m2)  .and.  (maxdt/(maxt-mint)<c1_m2)
!
            if ( ok2 .or. k>=maxiter2 ) exit
!
        end do
!
    end if
!
!   SAVE FINAL RESULTS.
!
    season(1_i4b:n) = season2(1_i4b:n)
    trend(1_i4b:n)  = trend2(1_i4b:n)
!
    if ( present(rw) ) then
        rw(1_i4b:n) = rw2(1_i4b:n)
    end if
!
    if ( present(no) ) then
        no = k
    end if
!
    if ( present(ok) ) then
        ok = ok2
    end if
!
!
! END OF SUBROUTINE comp_stlez_rv
! _______________________________
!
    end subroutine comp_stlez_rv
!
! =========================================================================================
!
    subroutine comp_stlez_rm( y, np, ns, isdeg, itdeg, robust, season, trend, ni, nt, nl, &
                              ildeg, nsjump, ntjump, nljump, maxiter, rw, no, ok )
       !dir$ attributes code_align : 32 :: comp_stlez_rm
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: comp_stlez_rm
!
! Purpose
! _______
!                                          
!   COMP_STLEZ decomposes the (time series) columns of a matrix into seasonal and 
!   trend components. It returns the components and, optionally, the robustness weights.
!
!   COMP_STLEZ offers an easy to use version of COMP_STL subroutine, also
!   included in STATPACK, by defaulting most parameters values associated
!   with the three LOESS smoothers used in COMP_STL.
!
!
! Arguments
! _________
!                                                  
!   Y        (INPUT) real(stnd), dimension(:,:)
!            On entry, the matrix to be decomposed.
!
!   NP       (INPUT) integer(i4b)
!            On entry, the period of the seasonal component. For example,
!            if  the  time series is monthly with a yearly cycle, then
!            NP=12 should be used. NP must be greater than 1.
!
!   NS       (INPUT) integer(i4b)
!            On entry, the length of the seasonal smoother.  The value of
!            NS should be an odd integer greater than or equal to 3; NS>6
!            is recommended. As NS increases  the  values  of  the
!            seasonal component at a given point in the seasonal cycle
!            (e.g., January values of a monthly series with  a  yearly
!            cycle) become smoother.
!
!   ISDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted  polynomial  in  seasonal
!            smoothing.  The value must be 0 or 1.
!
!   ITDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted  polynomial  in  trend
!            smoothing.  The value must be 0, 1 or 2.
!
!   ROBUST   (INPUT) logical(lgl)
!            On entry, TRUE if robustness iterations are to be used,
!            FALSE otherwise.
!            Robustness iterations are  carried out until convergence of
!            both seasonal and trend components, with MAXITER iterations
!            maximum.  Convergence occurs if the maximum changes in individual
!            seasonal  and  trend fits are less than 1% of the component's range
!            after the previous iteration.
!
!   SEASON   (OUTPUT) real(stnd), dimension(:,:)
!            On output, the seasonal components.
!
!            SEASON must verify: size(SEASON,1) = size(Y,1) and size(SEASON,2) = size(Y,2).
!
!   TREND    (OUTPUT) real(stnd), dimension(:,:)
!            On output, the trend components.
!
!            TREND must verify: size(TREND,1) = size(Y,1) and size(TREND,2) = size(Y,2).
!
!   NI       (INPUT, OPTIONAL) integer(i4b)
!            On entry, the number of loops for updating the seasonal and trend components.
!            The value of NI  should be a positive integer.
!
!            By default, NI=2 if ROBUST=FALSE and NI=1 if ROBUST=TRUE.
!   
!   NT       (INPUT, OPTIONAL) integer(i4b)
!            On entry, the length of the trend smoother. The value of
!            NT should be an odd integer greater than or equal to 3.
!
!            A value of NT between 1.5 * NP and 2 * NP is recommended.
!            As NT increases the values of the trend component become smoother.
!
!            By default, NT is set to the smallest odd integer greater
!            than or equal to  (1.5 * NP) / (1-(1.5/NS)).
!
!   NL       (INPUT, OPTIONAL) integer(i4b)
!            On entry, the length of the low-pass filter. The value 
!            of NL should be an odd integer greater than or equal to 3.
!
!            The smallest odd integer greater than or equal to NP is recommended.
!
!            By default, NL is set to the smallest odd integer greater
!            than or equal to NP.
!
!   ILDEG    (INPUT, OPTIONAL) integer(i4b)
!            On entry, the degree of locally-fitted polynomial in low-pass smoothing.
!
!            By default, ILDEG is set to ITDEG.
!
!
!   NSJUMP   (INPUT, OPTIONAL) integer(i4b)
!            On entry, the skipping value for seasonal smoothing.
!            The seasonal  smoother  skips  ahead NSJUMP points and then
!            linearly interpolates in between.  The  value of NSJUMP
!            should  be  a  positive  integer; if NSJUMP=1, a seasonal
!            smooth is calculated  at all size(Y) points. To  make  the
!            procedure  run  faster, a reasonable choice for NSJUMP is
!            10% or 20% of NS.
!
!            By default, NSJUMP is set to NS/10.
!
!   NTJUMP   (INPUT, OPTIONAL) integer(i4b)
!            On entry, the  skipping value for trend smoothing.
!
!            By default, NTJUMP is set to NT/10.
!
!   NLJUMP   (INPUT, OPTIONAL) integer(i4b)
!            On entry, the skipping value for the low-pass filter.
!
!            By default, NLJUMP is set to NL/10.
!
!   MAXITER  (INPUT, OPTIONAL) integer(i4b)
!            On entry, the maximum number of robustness iterations.
!
!            The default is 15. This argument is not used if ROBUST=FALSE.
!
!   RW       (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!            On output, final  robustness  weights. All RW elements are 1 if
!            ROBUST=FALSE .
!
!            RW must verify: size(RW,1) = size(Y,1) and size(RW,2) = size(Y,2).
!
!   NO       (OUTPUT, OPTIONAL) integer(i4b), dimension(:)
!            On output, if
!
!            - ROBUST=TRUE : NO(i) is the number of robustness iterations for each
!              time series Y(:,i). The  iterations end if a convergence criterion is met 
!              or if the number is MAXITER for each time series Y(:,i).
!            - ROBUST=FALSE : NO(:) is set to 0.
!
!            NO must verify:  size(NO) = size(Y,2).
!
!   OK       (OUTPUT, OPTIONAL) logical(lgl), dimension(:)
!            On output, if
!
!            - ROBUST=TRUE : OK(i) is set to TRUE if the convergence criterion is met
!              for time series Y(:,i) and to FALSE otherwise.
!            - ROBUST=FALSE : OK(:) is set to TRUE.
!
!            OK must verify:  size(OK) = size(Y,2).
!
!
! Further Details
! _______________
!
!   This subroutine is a FORTRAN90 implementation of subroutine STLEZ developped by
!   Cleveland and coworkers at AT&T Bell Laboratories.
!
!   At a minimum, COMP_STLEZ requires specifying the periodicity of the data
!   (e.g. NP, 12 for monthly), the width of the LOESS smoother used to smooth the
!   cyclic seasonal sub-series (e.g. NS) and the degree of the locally-fitted polynomial
!   in seasonal (e.g. ISDEG) and trend (e.g. ITDEG) smoothing.
!
!   COMP_STLEZ sets, by default, others parameters of the STL procedure to the values recommended
!   in Cleveland et al. (1990). It also includes tests of convergence if robust iterations
!   are carried out. Otherwise, COMP_STLEZ is similar to COMP_STL.
!
!   For further details, see:
!   
!   (1)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I.: 
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           Statistics Research Report, AT&T Bell Laboratories.
!   
!   (2)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I., 1990:
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           J. Official Stat., 6, 3-73.
!
!   (3)  Crotinger, J., 2017:
!            Java implementation of Seasonal-Trend-Loess time-series decomposition algorithm.
!            https://github.com/ServiceNow/stl-decomp-4j
!
!                                                                                                 
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use Select_Parameters,  only : omp_limit
    use omp_lib,            only : omp_get_num_procs, omp_get_max_threads, omp_in_parallel
    use Logical_Constants, only : true, false
    use Reals_Constants,   only : zero, one, half, c1_5, c0_9, c1_m2, ten
    use Char_Constants,    only : tseries_error3, tseries_error4, tseries_error5, tseries_error6
    use Utilities,         only : assert, assert_eq, merror
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:,:), intent(in)             :: y
    real(stnd), dimension(:,:), intent(out)            :: season, trend
    real(stnd), dimension(:,:), intent(out),  optional :: rw
!
    integer(i4b),               intent(in)             :: np, ns, isdeg, itdeg
    integer(i4b),               intent(in),   optional :: ni, nt, nl, ildeg, nsjump, ntjump, nljump, maxiter
    integer(i4b), dimension(:), intent(out),  optional :: no
!
    logical(lgl),               intent(in)             :: robust
    logical(lgl), dimension(:), intent(out),  optional :: ok
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd), dimension(size(y,1)+2*max(np,1_i4b),3) :: work
    !dir$ attributes align : 64 :: work
    real(stnd), dimension(size(y,1))                   :: trend2, season2, rw2
    !dir$ attributes align : 64 :: trend2
    !dir$ attributes align : 64 :: season2
    !dir$ attributes align : 64 :: rw2
    real(stnd)                                         :: maxs, mins, maxt, mint, maxds, maxdt
!
    integer(i4b) :: n, m, i, j, k, nt2, nl2, ni2, ns2, ildeg2, nsjump2, ntjump2, nljump2, maxiter2
!
    logical(lgl) :: ok2, use_no, use_ok, use_rw

    logical      :: test_par

!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='comp_stlez'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = assert_eq(  int(size(y,1),i4b),        &
                    int(size(season,1),i4b),   &
                    int(size(trend,1),i4b),    &
                    name_proc )
!
    m = assert_eq(  int(size(y,2),i4b),        &
                    int(size(season,2),i4b),   &
                    int(size(trend,2),i4b),    &
                    name_proc )
!
    use_rw = false
!
    if ( present(rw) ) then
        call assert( logical( n==int(size(rw,1),i4b), lgl ),     &
                     logical( m==int(size(rw,2),i4b), lgl ),     &
                     name_proc )
        use_rw = true
    end if
!
    use_no = false
!
    if ( present(no) ) then
        call assert( logical( m==int(size(no),i4b), lgl ),  name_proc )
        use_no = true
    end if
!
    use_ok = false
!
    if ( present(ok) ) then
        call assert( logical( m==int(size(ok),i4b), lgl ),  name_proc )
        use_ok = true
    end if
!
    if ( np<=1_i4b ) then
        call merror( name_proc//tseries_error4 )
    end if
!
    if ( isdeg/=1_i4b .and. isdeg/=0_i4b  ) then
        call merror( name_proc//tseries_error5 )
    end if
!
    if ( itdeg<0_i4b .or. itdeg>2_i4b  ) then
        call merror( name_proc//tseries_error3 )
    end if
!
    if ( present(ildeg) ) then
        if ( ildeg<0_i4b .or. ildeg>2_i4b  ) then
            call merror( name_proc//tseries_error6 )
        end if
        ildeg2 = ildeg
    else
        ildeg2 = itdeg
    end if
!
    ns2 = max( 3_i4b, ns )
!
    if ( mod(ns2,2_i4b)==0_i4b) ns2 = ns2 + 1_i4b
!
    if ( present(nt) ) then
        nt2 = nt
    else
        nt2 = ( c1_5*real(np,stnd) )/( one - c1_5/real(ns2,stnd) ) + half
    end if
!
    nt2 = max( 3_i4b, nt2 )
!
    if ( mod(nt2,2_i4b)==0_i4b) nt2 = nt2 + 1_i4b
!
    if ( present(nl) ) then
        nl2 = max( 3_i4b, nl )
    else
        nl2 = np
    end if
!
    if ( mod(nl2,2_i4b)==0_i4b) nl2 = nl2 + 1_i4b
!
    if ( present(ni) ) then
        ni2 = max( ni, 1_i4b )
    else
        if ( robust ) then
            ni2 = 1_i4b
        else
            ni2 = 2_i4b
        end if
    end if
!
    if ( present(nsjump) ) then
        nsjump2 = min( max( 1_i4b, nsjump ), n-1_i4b )
    else
        nsjump2 = max( 1_i4b, int(real(ns2,stnd)/ten + c0_9,i4b) )
    end if
!
    if ( present(ntjump) ) then
        ntjump2 = min( max( 1_i4b, ntjump ), n-1_i4b )
    else
        ntjump2 = max( 1_i4b, int(real(nt,stnd)/ten + c0_9,i4b) )
    end if
!
    if ( present(nljump) ) then
        nljump2 = min( max( 1_i4b, nljump ), n-1_i4b )
    else
        nljump2 = max( 1_i4b, int(real(nl,stnd)/ten + c0_9,i4b) )
    end if
!
    if ( present(maxiter) ) then
        maxiter2 = max( maxiter, 1_i4b )
    else
        maxiter2 = 15_i4b
    end if
!

    i = omp_get_num_procs()
    k = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )      .and.      &
               (m*n)>=omp_limit                .and.      &
               i>1_i4b                         .and.      &
               k>1_i4b                         .and.      &
               m>1_i4b

!
!$OMP PARALLEL IF(test_par)                                                               &
!$OMP         ,PRIVATE(i,j,k,ok2,rw2,trend2,season2,work,maxs,mins,maxt,mint,maxds,maxdt) &
!$OMP         ,SHARED(np,ns2,isdeg,nsjump2,nl2,ildeg2,nljump2,nt2,itdeg,ntjump2)
!
!   INITIALIZATIONS.
!
    rw2(1_i4b:n)  = one
    ok2      = true
  
!
!$OMP DO SCHEDULE(STATIC,8) 
!
    !dir$ assume_aligned rw2:64
    !dir$ assume_aligned trend2:64
    !dir$ assume_aligned y:64
    !dir$ assume_aligned work:64
    !dir$ assume_aligned season2:64
    !dir$ assume_aligned season:64
    do i = 1_i4b, m
!
        trend2(1_i4b:n) = zero
        k          = 0_i4b
!
!       FIRST ITERATIONS.
!
        do j = 1_i4b, ni2
!
!           STEP1 : DETRENDING.
!
            work(1_i4b:n,1_i4b) = y(1_i4b:n,i) - trend2(1_i4b:n)
!
!           STEP2 : CYCLE-SUBSERIES SMOOTHING.
!
            call ss( work(1_i4b:n,1_i4b), np, ns2, isdeg, nsjump2, false,      &
                     rw2(1_i4b:n), work(1_i4b:n+2_i4b*np,2_i4b) )
!
!           call fts( work(:n+2*np,2), np, work(:n+2*np,3), work(:n+2*np,1) )
!
!           STEP3 : LOW-PASS FILTERING OF SMOOTHED CYCLE-SUBSERIES.
!
            call ma( work(1_i4b:n+2_i4b*np,2_i4b), np, work(1_i4b:n+2_i4b*np,3_i4b) )
            call ma( work(1_i4b:n+np+1,3_i4b), np, work(1_i4b:n+np+1,1_i4b) )
            call ma( work(1_i4b:n+2_i4b,1_i4b), 3_i4b, work(1_i4b:n+2_i4b,3_i4b) )

            call ess( work(1_i4b:n,3_i4b), nl2, ildeg2, nljump2, false,     &
                      rw2(1_i4b:n), work(1_i4b:n,1_i4b) )
!
!           STEP4 : DETRENDING OF SMOOTHED CYCLE-SUBSERIES.
!
            season2(1_i4b:n) = work(np+1_i4b:n+np,2_i4b) - work(1_i4b:n,1_i4b)
!
!           STEP5 : DESEASONALIZING.
!
            work(1_i4b:n,1_i4b) = y(1_i4b:n,i) - season2(1_i4b:n)
!
!           STEP6 : TREND SMOOTHING.
!
            call ess( work(1_i4b:n,1_i4b), nt2, itdeg, ntjump2, false,           &
                      rw2(1_i4b:n), trend2(1_i4b:n) )
!
        end do
!
!       ROBUST ITERATIONS.
!
        if ( robust ) then
!
              !dir$ assume_aligned rw2:64
              !dir$ assume_aligned trend2:64
              !dir$ assume_aligned y:64
              !dir$ assume_aligned work:64
              !dir$ assume_aligned season2:64
              !dir$ assume_aligned season:64
            do
!
!               UPDATE SEASONAL AND TREND COMPONENTS.
!
                season(1_i4b:n,i) = season2(1_i4b:n)
                trend(1_i4b:n,i)  = trend2(1_i4b:n)
!
                k = k + 1_i4b
!
!               COMPUTE ROBUST WEIGHTS.
!
                work(1_i4b:n,1_i4b) = trend2(1_i4b:n) + season2(1_i4b:n)
                call rwts( y(1_i4b:n,i), work(1_i4b:n,1_i4b), rw2(1_i4b:n) )
!
!               INNER LOOP.
!
                do j = 1_i4b, ni2
!
!                   STEP1 : DETRENDING STEP.
!
                    work(1_i4b:n,1_i4b) = y(1_i4b:n,i) - trend2(1_i4b:n)
!
!                   STEP2 : CYCLE-SUBSERIES SMOOTHING.
!
                    call ss( work(1_i4b:n,1_i4b), np, ns2, isdeg, nsjump2, true,      &
                             rw2(1_i4b:n), work(1_i4b:n+2_i4b*np,2_i4b) )
!
!                   call fts( work(:n+2*np,2), np, work(:n+2*np,3), work(:n+2*np,1) )
!
!                   STEP3 : LOW-PASS FILTERING OF SMOOTHED CYCLE-SUBSERIES.
!
                    call ma( work(1_i4b:n+2_i4b*np,2_i4b), np, work(1_i4b:n+2_i4b*np,3_i4b) )
                    call ma( work(1_i4b:n+np+1,3_i4b), np, work(1_i4b:n+np+1,1_i4b) )
                    call ma( work(1_i4b:n+2_i4b,1_i4b), 3_i4b, work(1_i4b:n+2_i4b,3_i4b) )

                    call ess( work(1_i4b:n,3_i4b), nl2, ildeg2, nljump2, false,    &
                              rw2(1_i4b:n), work(1_i4b:n,1_i4b) )
!
!                   STEP4 : DETRENDING OF SMOOTHED CYCLE-SUBSERIES.
!
                    season2(1_i4b:n) = work(np+1_i4b:n+np,2_i4b) - work(1_i4b:n,1_i4b)
!
!                   STEP5 : DESEASONALIZING.
!
                    work(1_i4b:n,1_i4b) = y(1_i4b:n,i) - season2(1_i4b:n)
!
!                   STEP6 : TREND SMOOTHING.
!
                    call ess( work(1_i4b:n,1_i4b), nt2, itdeg, ntjump2, true,       &
                              rw2(1_i4b:n), trend2(1_i4b:n) )
!
                end do
!
!               TEST CONVERGENCE OF THE ROBUST ITERATIONS.
!
                maxs  = maxval( season(1_i4b:n,i) )
                mins  = minval( season(1_i4b:n,i) )
                maxt  = maxval( trend(1_i4b:n,i) )
                mint  = minval( trend(1_i4b:n,i) )
                maxds = maxval( abs( season2(1_i4b:n) - season(1_i4b:n,i) ) )
                maxdt = maxval( abs( trend2(1_i4b:n) - trend(1_i4b:n,i) ) )
!
                ok2 = (maxds/(maxs-mins)<c1_m2)  .and.  (maxdt/(maxt-mint)<c1_m2)
!
                if ( ok2 .or. k>=maxiter2 ) exit
!
            end do
!
        end if
!
!       SAVE FINAL RESULTS.
!
        season(1_i4b:n,i) = season2(1_i4b:n)
        trend(1_i4b:n,i)  = trend2(1_i4b:n)
!
        if ( use_rw ) then
            rw(1_i4b:n,i) = rw2(1_i4b:n)
        end if
!
        if ( use_no ) then
            no(i) = k
        end if
!
        if ( use_ok) then
            ok(i) = ok2
        end if
!
    end do
!
!$OMP END DO
!
!$OMP END PARALLEL
!
!
! END OF SUBROUTINE comp_stlez_rm
! _______________________________
!
    end subroutine comp_stlez_rm
!
! =========================================================================================
!
    subroutine comp_stl_rv( y, np, ni, no, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, &
                            ns, nt, nl, rw, season, trend )
       !dir$ attributes code_align : 32 :: comp_stl_rv
       !dir$ attributes forceinline :: comp_stl_rv
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: comp_stl_rv
!
! Purpose
! _______
!                                                                              
!   COMP_STL decomposes a time series into seasonal and trend components using LOESS smoothers.
!   It returns the components and robustness weights.
!
!
! Arguments
! _________
!                                                                              
!   Y        (INPUT) real(stnd), dimension(:)
!            On entry, the time series to be decomposed.
!
!   NP       (INPUT) integer(i4b)
!            On entry, the period of the seasonal component. For example,
!            if  the  time series is monthly with a yearly cycle, then
!            NP=12 should be used. NP must be greater than 1.
!
!   NI       (INPUT) integer(i4b)
!            On entry, the number of loops for updating the seasonal and trend components.
!            The value of NI  should be a strictly positive integer.
!
!   NO       (INPUT) integer(i4b)
!            On entry, the number of robustness iterations.
!            The value of NO  should be a positive integer.
!
!   ISDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted  polynomial  in  seasonal
!            smoothing.  The value must be 0 or 1.
!
!   ITDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted  polynomial  in  trend
!            smoothing.  The value must be 0, 1 or 2.
!
!   ILDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted polynomial in low-pass smoothing.
!            The value must be 0, 1 or 2.
!
!   NSJUMP   (INPUT/OUTPUT) integer(i4b)
!            On entry, the skipping value for seasonal smoothing.
!            The seasonal  smoother  skips  ahead NSJUMP points and then
!            linearly interpolates in between.  The  value of NSJUMP
!            should  be  a  positive  integer; if NSJUMP=1, a seasonal
!            smooth is calculated  at all size(Y) points. To  make  the
!            procedure  run  faster, a reasonable choice for NSJUMP is
!            10% or 20% of NS.
!
!   NTJUMP   (INPUT/OUTPUT) integer(i4b)
!            On entry, the  skipping value for trend smoothing.
!
!   NLJUMP   (INPUT/OUTPUT) integer(i4b)
!            On entry, the skipping value for the low-pass filter.
!
!   NS       (INPUT/OUTPUT) integer(i4b)
!            On entry, the length of the seasonal smoother.  The value of
!            NS should be an odd integer greater than or equal to 3; NS>6
!            is recommended.   As NS  increases  the  values  of  the
!            seasonal component at a given point in the seasonal cycle
!            (e.g., January values of a monthly series with  a  yearly
!            cycle) become smoother.
!
!   NT       (INPUT/OUTPUT) integer(i4b)
!            On entry, the length of the trend smoother. The value of
!            NT should be an odd integer greater than or equal to 3.
!            A value of NT between 1.5 * NP and 2 * NP is recommended.
!            As NT increases the values of the trend component become smoother.
!
!   NL       (INPUT/OUTPUT) integer(i4b)
!            On entry, the length of the low-pass filter. The value 
!            of NL should be an odd integer greater than or equal to 3.
!            The smallest odd integer greater than or equal to NP is recommended.
!
!   RW       (OUTPUT) real(stnd), dimension(:)
!            On output, final  robustness  weights. All RW elements are 1 if
!            NO=0 .
!
!            RW must verify:  size(RW) = size(Y).
!
!   SEASON   (OUTPUT) real(stnd), dimension(:)
!            On output, the seasonal component.
!
!            SEASON must verify:  size(SEASON) = size(Y).
!
!   TREND    (OUTPUT) real(stnd), dimension(:)
!            On output, the trend component.
!
!            TREND must verify:  size(TREND) = size(Y).
!
!
! Further Details
! _______________
!
!   This subroutine is a FORTRAN90 implementation of subroutine STL developped by
!   Cleveland and coworkers at AT&T Bell Laboratories.
!
!   This subroutine decomposes a time series into seasonal, trend and residual components.
!   The algorithm uses LOESS interpolation and smoothers to smooth the time series and estimate
!   the seasonal (or harmonic) component and the trend. This process is iterative with many steps
!   and may include robustness iterations that take advantage of the weighted-least-squares underpinnings
!   of LOESS to remove the effects of outliers.
!
!   There are three LOESS smoothers in COMP_STL and each require three parameters: a width, a degree, and a jump.
!   The width specifies the number of data points that the local interpolation uses to smooth each point, the degree
!   specifies the degree of the local polynomial that is fit to the data, and the jump specifies how many points are
!   skipped between LOESS interpolations, with linear interpolation being done between these points.
!
!   The LOESS smoother for estimating the trend is specified with the following parameters: a width
!   (e.g. NT), a degree (e.g. ITDEG) and a jump (e.g. NTJUMP).
!
!   The LOESS smoother for estimating the seasonal component is specified with the following parameters:
!   a width (e.g. NS), a degree (e.g. ISDEG) and a jump (e.g. NSJUMP).
!
!   The LOESS smoother for low-pass filtering is specified with the following parameters: a width
!   (e.g. NL), a degree (e.g. ILDEG) and a jump (e.g. NLJUMP).
!
!   If the NO argument is set to an integer value greater than 0, the process includes also robustness
!   iterations that take advandages of the weighted-least-squares underpinnings of LOESS to
!   remove the effects of outliers.
!
!   Note that, finally, that this subroutine expects equally spaced data with no missing values.
!
!   For further details, see:
!   
!   (1)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I.: 
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           Statistics Research Report, AT&T Bell Laboratories.
!   
!   (2)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I., 1990:
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           J. Official Stat., 6, 3-73.
!
!   (3)  Crotinger, J., 2017:
!            Java implementation of Seasonal-Trend-Loess time-series decomposition algorithm.
!            https://github.com/ServiceNow/stl-decomp-4j
!
!                                                                                                 
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use Logical_Constants, only : true, false
    use Reals_Constants,   only : zero, one
    use Char_Constants,    only : tseries_error3, tseries_error4, tseries_error5, tseries_error6,   &
                                  tseries_error7, tseries_error8
    use Utilities,         only : assert_eq, merror
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:), intent(in)  :: y
    real(stnd), dimension(:), intent(out) :: rw, season, trend
!
    integer(i4b), intent(in)    :: np, ni, no, isdeg, itdeg, ildeg
!
    integer(i4b), intent(inout) :: ns, nt, nl, nsjump, ntjump, nljump
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd), dimension(size(y)+2*max(np,1_i4b),3) :: work
    !dir$ attributes align : 64 :: work
!
    integer(i4b) :: n, k, j
!
    logical(lgl) :: userw
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='comp_stl'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = assert_eq(  int(size(y),i4b),        &
                    int(size(rw),i4b),       &
                    int(size(season),i4b),   &
                    int(size(trend),i4b),    &
                    name_proc )
!
    if ( np<=1_i4b ) then
        call merror( name_proc//tseries_error4 )
    end if
!
    if ( no<0_i4b ) then
        call merror( name_proc//tseries_error7 )
    end if
!
    if ( ni<=0_i4b ) then
        call merror( name_proc//tseries_error8 )
    end if
!
    if ( isdeg/=1_i4b .and. isdeg/=0_i4b  ) then
        call merror( name_proc//tseries_error5 )
    end if
!
    if ( itdeg<0_i4b .or. itdeg>2_i4b  ) then
        call merror( name_proc//tseries_error3 )
    end if
!
    if ( ildeg<0_i4b .or. ildeg>2_i4b  ) then
        call merror( name_proc//tseries_error6 )
    end if
!
    ns = max( 3_i4b, ns )
    if ( mod(ns,2_i4b)==0_i4b) ns = ns + 1_i4b
!
    nt = max( 3_i4b, nt )
    if ( mod(nt,2_i4b)==0_i4b) nt = nt + 1_i4b
!
    nl = max( 3_i4b, nl )
    if ( mod(nl,2_i4b)==0_i4b) nl = nl + 1_i4b
!
    nsjump = min( max( 1_i4b, nsjump ), n-1_i4b )
    ntjump = min( max( 1_i4b, ntjump ), n-1_i4b )
    nljump = min( max( 1_i4b, nljump ), n-1_i4b )
!
!   INITIALIZATIONS.
!
    userw          = false
    trend(1_i4b:n) = zero
!
    if ( no==0_i4b ) rw(1_i4b:n) = one
!
!   OUTER LOOP.
!
    !dir$ assume_aligned work:64
    !dir$ assume_aligned y:64
    !dir$ assume_aligned trend:64
    !dir$ assume_aligned rw:64
    !dir$ assume_aligned season:64
    do k = 0_i4b, no
!
!       call onestp( y, np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump,     &
!                    ntjump, nljump, ni, userw, rw, season, trend, work )
!
!       INNER LOOP.
!
        do j = 1_i4b, ni
!
!           STEP1 : DETRENDING STEP.
!
            work(1_i4b:n,1_i4b) = y(1_i4b:n) - trend(1_i4b:n)
!
!           STEP2 : CYCLE-SUBSERIES SMOOTHING.
!
            call ss( work(1_i4b:n,1_i4b), np, ns, isdeg, nsjump, userw,    &
                     rw(1_i4b:n), work(1_i4b:n+2_i4b*np,2_i4b) )
!
!           call fts( work(:n+2*np,2), np, work(:n+2*np,3), work(:n+2*np,1) )
!
!           STEP3 : LOW-PASS FILTERING OF SMOOTHED CYCLE-SUBSERIES.
!
            call ma( work(1_i4b:n+2_i4b*np,2_i4b), np, work(1_i4b:n+2_i4b*np,3_i4b) )
            call ma( work(1_i4b:n+np+1,3_i4b), np, work(1_i4b:n+np+1,1_i4b) )
            call ma( work(1_i4b:n+2_i4b,1_i4b), 3_i4b, work(1_i4b:n+2_i4b,3_i4b) )
!
            call ess( work(1_i4b:n,3_i4b), nl, ildeg, nljump, false,     &
                      rw(1_i4b:n), work(1_i4b:n,1_i4b) )
!
!           STEP4 : DETRENDING OF SMOOTHED CYCLE-SUBSERIES.
!
            season(1_i4b:n) = work(np+1_i4b:n+np,2_i4b) - work(1_i4b:n,1_i4b)
!
!           STEP5 : DESEASONALIZING.
!
            work(1_i4b:n,1_i4b) = y(1_i4b:n) - season(1_i4b:n)
!
!           STEP6 : TREND SMOOTHING.
!
            call ess( work(1_i4b:n,1_i4b), nt, itdeg, ntjump, userw,         &
                      rw(1_i4b:n), trend(1_i4b:n) )
!
        end do
!
        if ( k<no ) then
!
!           COMPUTE ROBUST WEIGHTS.
!
            work(1_i4b:n,1_i4b) = trend(1_i4b:n) + season(1_i4b:n)
            call rwts( y(1_i4b:n), work(1_i4b:n,1_i4b), rw(1_i4b:n) )
            userw = true
!
        end if
!
    end do        
!
!
! END OF SUBROUTINE comp_stl_rv
! _____________________________
!
    end subroutine comp_stl_rv
!
! ===========================================================================================
!
    subroutine comp_stl_rm( y, np, ni, no, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, &
                            ns, nt, nl, rw, season, trend )
       !dir$ attributes code_align : 32 :: comp_stl_rm
       !dir$ attributes forceinline :: comp_stl_rm
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: comp_stl_rm
!
! Purpose
! _______
!                                          
!   COMP_STL decomposes the (time series) columns of a matrix into seasonal and trend
!   components using LOESS smoothers. It returns the components and robustness weights.
!
!
! Arguments
! _________
!                                          
!   Y        (INPUT) real(stnd), dimension(:,:)
!            On entry, the time series to be decomposed.
!
!   NP       (INPUT) integer(i4b)
!            On entry, the period of the seasonal component. For example,
!            if  the  time series is monthly with a yearly cycle, then
!            NP=12 should be used. NP must be greater than 1.
!
!   NI       (INPUT) integer(i4b)
!            On entry, the number of loops for updating the seasonal and trend components.
!            The value of NI  should be a strictly positive integer.
!
!   NO       (INPUT) integer(i4b)
!            On entry, the number of robustness iterations.
!            The value of NO  should be a positive integer.
!
!   ISDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted  polynomial  in  seasonal
!            smoothing.  The value must be 0 or 1.
!
!   ITDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted  polynomial  in  trend
!            smoothing.  The value must be 0, 1 or 2.
!
!   ILDEG    (INPUT) integer(i4b)
!            On entry, the degree of locally-fitted polynomial in low-pass smoothing.
!            The value must be 0, 1 or 2.
!
!   NSJUMP   (INPUT/OUTPUT) integer(i4b)
!            On entry, the skipping value for seasonal smoothing.
!            The seasonal  smoother  skips  ahead NSJUMP points and then
!            linearly interpolates in between.  The  value of NSJUMP
!            should  be  a  positive  integer; if NSJUMP=1, a seasonal
!            smooth is calculated  at all size(Y) points. To  make  the
!            procedure  run  faster, a reasonable choice for NSJUMP is
!            10% or 20% of NS.
!
!   NTJUMP   (INPUT/OUTPUT) integer(i4b)
!            On entry, the  skipping value for trend smoothing.
!
!   NLJUMP   (INPUT/OUTPUT) integer(i4b)
!            On entry, the skipping value for the low-pass filter.
!
!   NS       (INPUT/OUTPUT) integer(i4b)
!            On entry, the length of the seasonal smoother.  The value of
!            NS should be an odd integer greater than or equal to 3; NS>6
!            is recommended.   As NS  increases  the  values  of  the
!            seasonal component at a given point in the seasonal cycle
!            (e.g., January values of a monthly series with  a  yearly
!            cycle) become smoother.
!
!   NT       (INPUT/OUTPUT) integer(i4b)
!            On entry, the length of the trend smoother. The value of
!            NT should be an odd integer greater than or equal to 3.
!            A value of NT between 1.5 * NP and 2 * NP is recommended.
!            As NT increases the values of the trend component become smoother.
!
!   NL       (INPUT/OUTPUT) integer(i4b)
!            On entry, the length of the low-pass filter. The value 
!            of NL should be an odd integer greater than or equal to 3.
!            The smallest odd integer greater than or equal to NP is recommended.
!
!   RW       (OUTPUT) real(stnd), dimension(:,:)
!            On output, final  robustness  weights. All RW elemets are 1 if
!            NO=0 .
!
!            RW must verify: size(RW,1) = size(Y,1) and size(RW,2) = size(Y,2).
!
!   SEASON   (OUTPUT) real(stnd), dimension(:,:)
!            On output, the seasonal components.
!
!            SEASON must verify: size(SEASON,1) = size(Y,1) and size(SEASON,2) = size(Y,2).
!
!   TREND    (OUTPUT) real(stnd), dimension(:,:)
!            On output, the trend components.
!
!            TREND must verify: size(TREND,1) = size(Y,1) and size(TREND,2) = size(Y,2).
!
!
! Further Details
! _______________
!
!   This subroutine is a FORTRAN90 implementation of subroutine STL developped by
!   Cleveland and coworkers at AT&T Bell Laboratories.
!
!   This subroutine decomposes a multi-channel time series into seasonal, trend and residual components.
!   The algorithm uses LOESS interpolation and smoothers to smooth the multi-channel time series and estimate
!   the seasonal (or harmonic) components and the trends. This process is iterative with many steps
!   and may include robustness iterations that take advantage of the weighted-least-squares underpinnings
!   of LOESS to remove the effects of outliers.
!
!   There are three LOESS smoothers in COMP_STL and each require three parameters: a width, a degree, and a jump.
!   The width specifies the number of data points that the local interpolation uses to smooth each point, the degree
!   specifies the degree of the local polynomial that is fit to the data, and the jump specifies how many points are
!   skipped between LOESS interpolations, with linear interpolation being done between these points.
!
!   The LOESS smoother for estimating the trend is specified with the following parameters: a width
!   (e.g. NT), a degree (e.g. ITDEG) and a jump (e.g. NTJUMP).
!
!   The LOESS smoother for estimating the seasonal component is specified with the following parameters:
!   a width (e.g. NS), a degree (e.g. ISDEG) and a jump (e.g. NSJUMP).
!
!   The LOESS smoother for low-pass filtering is specified with the following parameters: a width
!   (e.g. NL), a degree (e.g. ILDEG) and a jump (e.g. NLJUMP).
!
!   If the NO argument is set to an integer value greater than 0, the process includes also robustness
!   iterations that take advandages of the weighted-least-squares underpinnings of LOESS to
!   remove the effects of outliers.
!
!   Note that, finally, that this subroutine expects equally spaced data with no missing values.
!
!   For further details, see:
!   
!   (1)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I.: 
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           Statistics Research Report, AT&T Bell Laboratories.
!   
!   (2)  Cleveland, R.B., Cleveland, W.S., McRae, J.E., and Terpenning, I., 1990:
!           STL: A Seasonal-Trend Decomposition  Procedure Based on Loess.
!           J. Official Stat., 6, 3-73.
!
!   (3)  Crotinger, J., 2017:
!            Java implementation of Seasonal-Trend-Loess time-series decomposition algorithm.
!            https://github.com/ServiceNow/stl-decomp-4j
!
!                                                                                                 
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use Select_Parameters,  only : omp_limit
    use omp_lib,            only : omp_get_num_procs, omp_get_max_threads, omp_in_parallel
    use Logical_Constants, only : true, false
    use Reals_Constants,   only : zero, one
    use Char_Constants,    only : tseries_error3, tseries_error4, tseries_error5, tseries_error6,   &
                                  tseries_error7, tseries_error8
    use Utilities,         only : assert_eq, merror
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:,:), intent(in)  :: y
    real(stnd), dimension(:,:), intent(out) :: rw, season, trend
!
    integer(i4b), intent(in)    :: np, ni, no, isdeg, itdeg, ildeg
!
    integer(i4b), intent(inout) :: ns, nt, nl, nsjump, ntjump, nljump
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd), dimension(size(y)+2*max(np,1_i4b),3) :: work
    !dir$ attributes align : 64 :: work
!
    integer(i4b) :: m, n, k, j, i
!
    logical(lgl) :: userw

    logical      :: test_par

!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='comp_stl'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = assert_eq(  int(size(y,1),i4b),        &
                    int(size(rw,1),i4b),       &
                    int(size(season,1),i4b),   &
                    int(size(trend,1),i4b),    &
                    name_proc )
!
    m = assert_eq(  int(size(y,2),i4b),        &
                    int(size(rw,2),i4b),       &
                    int(size(season,2),i4b),   &
                    int(size(trend,2),i4b),    &
                    name_proc )
!
    if ( np<=1_i4b ) then
        call merror( name_proc//tseries_error4 )
    end if
!
    if ( no<0_i4b ) then
        call merror( name_proc//tseries_error7 )
    end if
!
    if ( ni<=0_i4b ) then
        call merror( name_proc//tseries_error8 )
    end if
!
    if ( isdeg/=1_i4b .and. isdeg/=0_i4b  ) then
        call merror( name_proc//tseries_error5 )
    end if
!
    if ( itdeg<0_i4b .or. itdeg>2_i4b  ) then
        call merror( name_proc//tseries_error3 )
    end if
!
    if ( ildeg<0_i4b .or. ildeg>2_i4b  ) then
        call merror( name_proc//tseries_error6 )
    end if
!
    ns = max( 3_i4b, ns )
    if ( mod(ns,2_i4b)==0_i4b) ns = ns + 1_i4b
!
    nt = max( 3_i4b, nt )
    if ( mod(nt,2_i4b)==0_i4b) nt = nt + 1_i4b
!
    nl = max( 3_i4b, nl )
    if ( mod(nl,2_i4b)==0_i4b) nl = nl + 1_i4b
!
    nsjump = min( max( 1_i4b, nsjump ), n-1_i4b )
    ntjump = min( max( 1_i4b, ntjump ), n-1_i4b )
    nljump = min( max( 1_i4b, nljump ), n-1_i4b )
!
!   INITIALIZATIONS.
!
    trend(1_i4b:n,:m) = zero
!
    if ( no==0_i4b ) rw(1_i4b:n,:m) = one
!

    i = omp_get_num_procs()
    k = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )      .and.      &
               (m*n)>=omp_limit                .and.      &
               i>1_i4b                         .and.      &
               k>1_i4b                         .and.      &
               m>1_i4b

    !dir$ assume_aligned work:64
    !dir$ assume_aligned y:64
    !dir$ assume_aligned trend:64
    !dir$ assume_aligned rw:64
    !dir$ assume_aligned season:64
!
!$OMP PARALLEL DO IF(test_par)                &
!$OMP            ,PRIVATE(i,j,k,userw,work)   &
!$OMP            ,SHARED(np,ns,nl,nt,isdeg,ildeg,itdeg,nsjump,nljump,ntjump)
!
    do i = 1_i4b, m
!
        userw = false
!
!       OUTER LOOP.
!
        do k = 0_i4b, no
!
!           INNER LOOP.
!
            do j = 1_i4b, ni
!
!               STEP1 : DETRENDING STEP.
!
                work(1_i4b:n,1_i4b) = y(1_i4b:n,i) - trend(1_i4b:n,i)
!
!               STEP2 : CYCLE-SUBSERIES SMOOTHING.
!
                call ss( work(1_i4b:n,1_i4b), np, ns, isdeg, nsjump, userw,      &
                         rw(1_i4b:n,i), work(1_i4b:n+2_i4b*np,2_i4b) )
!
!               STEP3 : LOW-PASS FILTERING OF SMOOTHED CYCLE-SUBSERIES.
!
                call ma( work(1_i4b:n+2_i4b*np,2_i4b), np, work(1_i4b:n+2_i4b*np,3_i4b) )
                call ma( work(1_i4b:n+np+1,3_i4b), np, work(1_i4b:n+np+1,1_i4b) )
                call ma( work(1_i4b:n+2_i4b,1_i4b), 3_i4b, work(1_i4b:n+2_i4b,3_i4b) )

                call ess( work(1_i4b:n,3_i4b), nl, ildeg, nljump, false ,    &
                          rw(1_i4b:n,i), work(1_i4b:n,1_i4b) )
!
!               STEP4 : DETRENDING OF SMOOTHED CYCLE-SUBSERIES.
!
                season(1_i4b:n,i) = work(np+1_i4b:n+np,2_i4b) - work(1_i4b:n,1_i4b)
!
!               STEP5 : DESEASONALIZING.
!
                work(1_i4b:n,1_i4b) = y(1_i4b:n,i) - season(1_i4b:n,i)
!
!               STEP6 : TREND SMOOTHING.
!
                call ess( work(1_i4b:n,1_i4b), nt, itdeg, ntjump, userw,         &
                          rw(1_i4b:n,i), trend(1_i4b:n,i) )
!
            end do
!
            if ( k<no ) then
!
!               COMPUTE ROBUST WEIGHTS.
!
                work(1_i4b:n,1_i4b) = trend(1_i4b:n,i) + season(1_i4b:n,i)
                call rwts( y(1_i4b:n,i), work(1_i4b:n,1_i4b), rw(1_i4b:n,i) )
                userw = true
!
            end if
!
         end do
!
    end do
!
!$OMP END PARALLEL DO
!
!
! END OF SUBROUTINE comp_stl_rm
! _____________________________
!
    end subroutine comp_stl_rm
!
! ===========================================================================================
!
    subroutine ma( x, len, ave)
       !dir$ attributes code_align : 32 :: ma
       !dir$ attributes forceinline :: ma
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: ma
!
! Purpose
! _______
!                                                  
!   Smooth the vector X with a moving average of length LEN and output the result
!   in the vector AVE.
!
!
! Arguments
! _________
!                                                                              
!   X    (INPUT) real(stnd), dimension(:)
!        On entry, the vector to smooth.
!
!   LEN  (INPUT) integer(i4b)
!        On entry, the length of the moving average.
!        The argument LEN must be >=1 and < size(x).
!
!
!   AVE  (OUTPUT) real(stnd), dimension(size(x))
!        On output, AVE(1:size(X)-LEN+1) contains the smoothed values and
!        AVE(size(X)-LEN+2:size(X)) is unchanged.
!
!
! Further Details
! _______________
!
!   This subroutine is a low-level subroutine used by subroutines COMP_STLEZ and COMP_STL.
!
!                                                                                                 
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use omp_lib
    use Reals_Constants,   only : one
    use Char_Constants,    only : tseries_error9
    use Utilities,         only : merror
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:),       intent(in)  :: x
    real(stnd), dimension(size(x)), intent(out) :: ave
!
    integer(i4b), intent(in)  :: len
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd)   :: v
!
    integer(i4b) :: n, j, k, m, newn
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='ma'
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( x )
!
    if ( len<1_i4b .or. len>=n ) then
       call merror( name_proc//tseries_error9 )
    end if
!
    newn = n - len + 1_i4b
!
    v = sum( x(1_i4b:len) )
    ave(1_i4b) = v
    k = len
    m = 0_i4b
!
    !dir$ assume_aligned x:64
    !dir$ assume_aligned ave:64
    !$omp simd simdlen(8) reduction(-:v)
    do j=2_i4b, newn
!
        k = k + 1_i4b
        m = m + 1_i4b
!
        v = v - x(m) + x(k)
        ave(j) = v
!
    end do
!
    ave(1_i4b:newn) = ave(1_i4b:newn)*( one/real(len, stnd) )
!
!
! END OF SUBROUTINE ma
! ____________________
!
    end subroutine ma
!
!
! =========================================================================================
!                                    DETRENDING PROCEDURES
! =========================================================================================
!
!
    subroutine detrend_rv( vec, trend, orig, slope )
       !dir$ attributes code_align : 32 :: detrend_rv
       !dir$ attributes forceinline :: detrend_rv
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: detrend_rv
!
! Purpose
! _______
!
!   Subroutine DETREND detrends a time series (e.g. the argument VEC).
!
!
! Arguments
! _________
!
!   VEC          (INPUT/OUTPUT) real(stnd), dimension(:)
!                The time series vector to be detrended.
!
!   TREND        (INPUT) integer(i4b)
!                If:
!
!                - TREND=1 The mean of the time series is removed
!                - TREND=2 The drift from the time series is removed by using the formula:
!
!                                  drift = (VEC(size(VEC)) - VEC(1))/(size(VEC) - 1)
!
!                - TREND=3 The least-squares line from the time series is removed.
!
!                For other values of TREND nothing is done.
!
!   ORIG         (OUTPUT, OPTIONAL) real(stnd)
!                On exit, the constant term if TREND=1 or 3.
!
!   SLOPE        (OUTPUT, OPTIONAL) real(stnd)
!                On exit, the linear term if TREND=2 or 3.
!
!
! Further Details
! _______________
!
!   On exit, the original time series may be recovered with the formula
!
!             VEC(i) = VEC(i) + ORIG + SLOPE * real(i-1,stnd)
!
!   for i=1, size(vec), in all the cases.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use omp_lib
    use Utilities,         only : arth
    use Reals_Constants,   only : zero, one, three, six
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:)  :: vec
    real(stnd),   intent(out), optional        :: orig, slope
!
    integer(i4b), intent(in) :: trend
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, k
!
    real(stnd)   :: orig2, slope2
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='detrend'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m = size( vec )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( m<=0_i4b ) then
!
        if ( present( orig) )  orig  = zero
        if ( present( slope) ) slope = zero
!
        return
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    orig2  = zero
    slope2 = zero
!
    select case( abs(trend) )
!
       case( 1_i4b )
!
!           REMOVE MEAN.
!
            !dir$ assume_aligned vec:64
            orig2    = sum( vec(:m) )/real( m, stnd )
            vec(:m) = vec(:m) - orig2
!
        case( 2_i4b )
!
!           REMOVE DRIFT.
!
            !dir$ assume_aligned vec:64
            slope2 = (vec(m) - vec(1_i4b))/real( m-1_i4b, stnd )
            vec(:m) = vec(:m) - slope2*arth( zero, one, m )           
!
        case( 3_i4b )
!
!           REMOVE LINEAR LEAST SQUARES LINE.
!
             !dir$ assume_aligned vec:64
             !$omp simd simdlen(8) reduction(+:orig2)
            do k = 1_i4b, m
                slope2 = real( k, stnd )*vec(k) - orig2 + slope2
                orig2 = orig2 + vec(k)
            end do
!
            orig2  = (orig2 - three*slope2/real(m+1_i4b,stnd))/real(m,stnd)
            slope2 = slope2/(real( m*(m-1_i4b)*(m+1_i4b), stnd)/six)
            vec(:m) = vec(:m) - ( orig2 + slope2*arth( zero, one, m ) )       
!
    end select
!   
    if ( present(orig) ) then
        orig = orig2
    end if
!   
    if ( present(slope) ) then
       slope = slope2
    end if
!
!
! END OF SUBROUTINE detrend_rv
! ____________________________
!
    end subroutine detrend_rv
!
! =========================================================================================
!
    subroutine detrend_rm( mat, trend, orig, slope )
       !dir$ attributes code_align : 32 :: detrend_rm
       !dir$ attributes forceinline :: detrend_rm
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: detrend_rm
!
! Purpose
! _______
!
!   Subroutine DETREND detrends a multi-channel time series (e.g. the argument MAT).
!   Each row of matrix MAT is a real time series
!
!
! Arguments
! _________
!
!   MAT          (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                The multi-channel time series matrix to be detrended.
!
!   TREND        (INPUT) integer(i4b)
!                If:
!
!                - TREND=1 The means of the time series are removed
!                - TREND=2 The drifts from the time series are removed by using the formula:
!
!                            drift(:) = (MAT(:,size(MAT,2)) - MAT(:,1))/(size(MAT,2) - 1)
!
!                - TREND=3 The least-squares lines from the time series are removed.
!
!                For other values of TREND nothing is done.
!
!   ORIG         (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                On exit, the constant terms if TREND=1 or 3.
!
!                The size of ORIG must verify:  size(ORIG) = size(MAT,1) .
!                
!
!   SLOPE        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                On exit, the linear terms if TREND=2 or 3.
!
!                The size of SLOPE must verify:  size(SLOPE) = size(MAT,1) .
!
!
! Further Details
! _______________
!
!   On exit, the original time series may be recovered with the formula
!
!             MAT(j,i) = MAT(j,i) + ORIG(j) + SLOPE(j) * real(i-1,stnd)
!
!   for i=1, size(MAT,2) and j=1, size(MAT,1), in all the cases.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use omp_lib
    use Utilities,         only : assert
    use Reals_Constants,   only : one, zero, three, six
    use Logical_Constants, only : true, false
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout),  dimension(:,:)         :: mat
    real(stnd),   intent(out),    dimension(:), optional :: orig, slope
!
    integer(i4b), intent(in) :: trend
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, n, k
!
    real(stnd)                         :: c
    real(stnd), dimension(size(mat,1)) :: orig2, slope2
    !dir$ attributes align : 64 :: orig2,slope2
!
    logical(lgl) :: out_orig, out_slope
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='detrend'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
    out_orig = false
!
    if ( present(orig) )  then
!
        call assert( logical(int(size(orig),i4b)==n,lgl),     &
                     name_proc )
!
        orig(:n) = zero
        out_orig = true
!
    end if
!
    out_slope = false
!
    if ( present(slope) )  then
!
        call assert( logical(int(size(slope),i4b)==n,lgl),     &
                     name_proc )
!
        slope(:n) = zero
        out_slope = true
!
    end if
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( m<=0_i4b .or. n<=0_i4b  ) then
!
        return
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    orig2(:n)  = zero
    slope2(:n) = zero
!
    select case( abs(trend) )
!
        case( 1_i4b )
!
!           REMOVE MEANS.
!
             !dir$ assume_aligned orig2:64
             !dir$ assume_aligned mat:64
            orig2(:n)    = sum( mat(:n,:m),dim=2 )*( one/real( m, stnd ) )
!
            !$omp simd simdlen(8) reduction(-:mat)
            do k = 1_i4b, m
                mat(:n,k) = mat(:n,k) - orig2(:n)
            end do
!
        case( 2_i4b )
!
!           REMOVE DRIFTS.
!
            !dir$ assume_aligned slope2:64
            !dir$ assume_aligned mat:64
            slope2(:n) = (mat(:n,m) - mat(:n,1_i4b))*( one/real( m-1_i4b, stnd ) )
!
            !$omp simd simdlen(8) reduction(-:mat)
            do k = 1_i4b, m
                c         = real( k, stnd ) - one
                mat(:n,k) = mat(:n,k) - c*slope2(:n)           
            end do
!
        case( 3_i4b )
!
!           REMOVE LINEAR LEAST SQUARES LINES.
!
            !dir$ assume_aligned slope2:64
            !dir$ assume_aligned mat:64
            !dir$ assume_aligned orig2:64
            do k = 1_i4b, m
                slope2(:n) = real( k, stnd )*mat(:n,k) - orig2(:n) + slope2(:n)
                orig2(:n)  = orig2(:n) + mat(:n,k)
            end do
!
            orig2(:n)  = ( orig2(:n) - (three/real(m+1_i4b,stnd))*slope2(:n) )*( one/real(m,stnd) )
            slope2(:n) = slope2(:n)*( one/(real( m*(m-1_i4b)*(m+1_i4b), stnd)/six) )
!
             !$omp simd simdlen(8) reduction(-:mat)
            do k = 1_i4b, m
                c         = real( k, stnd ) - one
                mat(:n,k) = mat(:n,k) - ( c*slope2(:n) + orig2(:n) )         
            end do
!
    end select
!   
    if ( out_orig ) then
        orig(:n) = orig2(:n)
    end if
!   
    if ( out_slope ) then
       slope(:n) = slope2(:n)
    end if
!
!
! END OF SUBROUTINE detrend_rm
! ____________________________
!
    end subroutine detrend_rm
!
!
! =========================================================================================
!                               FREQUENCY FILTERING PROCEDURES
! =========================================================================================
!
!
    subroutine hwfilter_rv( vec, pl, ph, initfft, trend, win )
       !dir$ attributes code_align : 32 :: hwfilter_rv
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: hwfilter_rv
!
! Purpose
! _______
!
!   Subroutine HWFILTER filters a time series (e.g. the argument VEC) in the frequency band
!   limited by periods PL and PH by windowed filtering (PL and PH are expressed in number
!   of points, i.e. PL=6(18) and PH=32(96) selects periods between 1.5 yrs and 8 yrs for
!   quarterly (monthly) data).
!
!
! Arguments
! _________
!
!   VEC          (INPUT/OUTPUT) real(stnd), dimension(:)
!                The time series vector to be filtered.
!
!                Size(VEC) must be greater or equal to 4.
!
!   PL           (INPUT) integer(i4b)
!                Minimum period of oscillation of desired component.
!                Use PL=0 for high-pass filtering frequencies corresponding
!                to periods shorter than PH, 
!                PL must be equal to 0 or greater or equal to 2.
!                Moreover, PL must be less or equal to size(VEC).
!
!   PH           (INPUT) integer(i4b)
!                Maximum period of oscillation of desired component.
!                USE PH=0 for low-pass filtering frequencies corresponding
!                to periods longer than PL.
!                PH must be equal to 0 or greater or equal to 2.
!                Moreover, PH must be less or equal to size(VEC).
!
!   INITFFT      (INPUT, OPTIONAL) logical(lgl)
!                On entry, if INITFFT is set to false, it is assumed that a call to subroutine
!                INIT_FFT has been done before calling subroutine HWFILTER in order to 
!                sets up constants and functions for use by subroutine FFT_ROW which is called inside
!                subroutine HWFILTER (the call to INIT_FFT must have the following form: 
!
!                     call init_fft( size(VEC) )
!
!                If INITFFT is set to true, the call to INIT_FFT is done inside subroutine
!                HWFILTER and a call to END_FFT is also done before leaving
!                subroutine HWFILTER.
!
!                The default is INITFFT=true.
!
!   TREND        (INPUT, OPTIONAL) integer(i4b)
!                If:
!
!                - TREND=+/-1 The mean of the time series is removed before time filtering
!                - TREND=+/-2 The drift from the time series is removed before time filtering
!                  by using the formula: drift = (VEC(size(VEC)) - VEC(1))/(size(VEC) - 1)
!                - TREND=+/-3 The least-squares line from the time series is removed before
!                  time filtering.
!
!                IF TREND=-1,-2 or -3, the mean, drift or least-squares line is reintroduced
!                post-filtering, respectively.
!                For other values of TREND nothing is done before or after filtering.
!
!   WIN          (INPUT, OPTIONAL) real(stnd)
!                By default, Hamming window filtering is used (i.e. WIN=0.54).
!                SET WIN=0.5 for Hanning window or WIN=1 for rectangular window.
!
!                WIN must be greater or equal to O.5 and less or equal to 1.
!
!
! Further Details
! _______________
!
!   Use PL=0 for high-pass filtering frequencies corresponding to periods shorter than PH, 
!   or PH=0 for low-pass filtering frequencies corresponding to periods longer than PL.
!
!   Setting PH<PL is also allowed and performs band rejection of periods between PH and PL
!   (i.e. in that case the meaning of the PL and PH arguments are reversed).
!
!   Examples: 
!
!   For quarterly data:
!      call hwfilter( vec, pl=6, ph=32) returns component with periods between 1.5 and 8 yrs.
!
!   For monthly data:
!      call hwfilter( vec, pl=0, ph=24) returns component with all periods less than 2 yrs.
!
!   For more details and algorithm, see:
!
!   (1) Iacobucci, A., and Noullez, A., 2005:
!           A Frequency Selective Filter for Short-Length Time Series.
!           Computational Economics, 25,75-102.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    !use omp_lib
    use Utilities,         only : merror, arth
    use Char_Constants,    only : tseries_error10, tseries_error11, tseries_error12,  &
                                  tseries_error13
    use Reals_Constants,   only : one, zero, half
    use Logical_Constants, only : true, false
    use FFT_Procedures,    only : fft_row, init_fft, end_fft
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:)  :: vec
    real(stnd),   intent(in), optional         :: win
!
    integer(i4b), intent(in)                :: pl, ph
    integer(i4b), intent(in), optional      :: trend
!
    logical(lgl),  intent(in), optional :: initfft
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, pl2, ph2, k, km, kodd, trend2
!
    real(stnd)   :: comp, win2, hwin, orig, slope, a, b, d
    real(stnd), dimension((size(vec)/2)+1) :: h2
    real(stnd), dimension(size(vec))       :: h
    !dir$ attributes align : 64 :: h2
    !dir$ attributes align : 64 :: h
!
    complex(stnd), dimension(size(vec)) :: vect
    !dir$ attributes align : 64 ::  vect
!
    logical(lgl) :: odd, initfft2
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='hwfilter'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m = size( vec )
!
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    km = m/2_i4b
!
    if ( m==2_i4b*km ) then
        odd  = false
        kodd = km
    else
        odd  = true
        kodd = km + 1_i4b
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
!   CHECK PERIODS ORDERING TO SELECT FREQUENCIES.
!
    if ( ph<pl ) then
        pl2   = ph
        ph2   = pl
        comp  = one
    else
        pl2   = pl
        ph2   = ph
        comp  = zero
    end if
!
    if ( pl2<0_i4b)        &
    call merror( name_proc//tseries_error11 )
!
    if ( ph==1_i4b .or. pl==1_i4b )       &
    call merror( name_proc//tseries_error12 )
!
    if ( m<ph2)            &
    call merror( name_proc//tseries_error13 )
!
!   USE HAMMING WINDOW BY DEFAULT.
!
    win2 = 0.54_stnd
!
    if ( present(win) ) then
!
!       IF WINDOW IS NOT POSITIVE DEFINITE, USE HAMMING WINDOW.
!
        if ( win>=half .and. win<=one ) win2 = win
    end if    
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
!
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
!   BUILD UP IDEAL FILTER GAIN VECTOR.
!
    do k=0_i4b, km
!
!       SELECT FREQUENCY BAND.
!
        if ( k*ph2>=m .and. k*pl2<=m ) then
            h(k+1_i4b) = one - comp
        else
            h(k+1_i4b) = comp
        end if
!
    end do
!
!   CONVOLVE WITH WINDOW TRANSFORM TO GET FILTER FREQUENCY RESPONSE IF NEEDED.
!
    
    if ( win2/=one ) then
!
        hwin = (one-win2)*half
        a    = h(2_i4b)
!
        if ( odd ) then
            b = h(km+1_i4b)
        else
            b = h(km)
        end if
!
         !dir$ assume_aligned h:64
         !$omp simd simdlen(8) linear(k:1)
        do k = 1_i4b, km
            d    = h(k)
            h(k) = win2*d + hwin*( a+h(k+1_i4b) )
            a    = d
        end do
!
        h(km+1_i4b) = win2*h(km+1_i4b) + hwin*(a+b)
!
    end if
!
    h2(2_i4b:kodd)       = h(2_i4b:kodd)
    h(m:km+2_i4b:-1_i4b) = h2(2_i4b:kodd)
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    k = abs( trend2 )
!
    if ( k>=1_i4b .and. k<=3_i4b ) then
        call detrend_rv( vec(:m), k, orig, slope )
    end if
!
!   INITIALIZE THE FFT SUBROUTINE IF REQUIRED.
!
    if ( initfft2 ) then
        call init_fft( m )
    end if
!
!   TRANSFORM THE REAL SEQUENCE.
!
    vect(:m) = cmplx( vec(:m), zero, kind=stnd )
!
    call fft_row( vect(:m), forward=true  )
!
!   NOW, FILTER THE TIME SERIES.
!   COMPUTE FILTERED FOURIER TRANSFORM OF THE SERIES.
!
    vect(:m) = h(:m)*vect(:m)
!
!   INVERT THE SEQUENCE BACK.
!
    call fft_row( vect(:m), forward=false )
!   
    vec(:m) = real( vect(:m) )
!
!   DEALLOCATE WORK ARRAYS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!   ADD MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEAN.
!
            !dir$ assume_aligned vec:64
            vec(:m) = vec(:m) + orig
!
        case( -2_i4b )
!
!           ADD DRIFT.
!
            !dir$ assume_aligned vec:64
            vec(:m) = vec(:m) + slope*arth( zero, one, m )           
!
!            do k = 1_i4b, m
!                a      = real( k, stnd ) - one
!                vec(k) = vec(k) + a*slope           
!            end do
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINE.
!         
            !dir$ assume_aligned vec:64
            vec(:m) = vec(:m) + ( orig + slope*arth( zero, one, m ) )       
!         
!            do k = 1_i4b, m
!                a      = real( k, stnd ) - one
!                vec(k) = vec(k) + a*slope + orig      
!            end do
!
    end select
!
!
! END OF SUBROUTINE hwfilter_rv
! _____________________________
!
    end subroutine hwfilter_rv
!
! =========================================================================================
!
    subroutine hwfilter_rm( mat, pl, ph, initfft, trend, win, max_alloc )
       !dir$ attributes code_align : 32 :: hwfilter_rm
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: hwfilter_rm
!
! Purpose
! _______
!
!   Subroutine HWFILTER filters a multi-channel time series (e.g. the argument MAT) in the frequency
!   band limited by periods PL and PH by windowed filtering (PL and PH are expressed in number
!   of points, i.e. PL=6(18) and PH=32(96) selects periods between 1.5 yrs and 8 yrs for
!   quarterly (monthly) data).
!
!
! Arguments
! _________
!
!   MAT          (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                The multi-channel time series matrix to be filtered.
!                Each column of MAT corresponds to one observation.
!
!                Size(MAT,2) must be greater or equal to 4.
!
!   PL           (INPUT) integer(i4b)
!                Minimum period of oscillation of desired component.
!                Use PL=0 for high-pass filtering frequencies corresponding
!                to periods shorter than PH, 
!                PL must be equal to 0 or greater or equal to 2.
!                Moreover, PL must be less or equal to size(MAT,2).
!
!   PH           (INPUT) integer(i4b)
!                Maximum period of oscillation of desired component.
!                USE PH=0 for low-pass filtering frequencies corresponding
!                to periods longer than PL.
!                PH must be equal to 0 or greater or equal to 2.
!                Moreover, PH must be less or equal to size(MAT,2).
!
!   INITFFT      (INPUT, OPTIONAL) logical(lgl)
!                On entry, if INITFFT is set to false, it is assumed that a call to subroutine
!                INIT_FFT has been done before calling subroutine HWFILTER in order to 
!                sets up constants and functions for use by subroutine FFT_ROW which is called inside
!                subroutine HWFILTER (the call to INIT_FFT must have the following form: 
!
!                     call init_fft( shape(MAT), dim=2_i4b )
!
!                If INITFFT is set to true, the call to INIT_FFT is done inside subroutine
!                HWFILTER and a call to END_FFT is also done before leaving
!                subroutine HWFILTER.
!
!                The default is INITFFT=true.
!
!   TREND        (INPUT, OPTIONAL) integer(i4b)
!                If:
!
!                - TREND=+/-1 The means of the time series are removed before time filtering
!                - TREND=+/-2 The drifts from the time series are removed before time filtering
!                  by using the formula: drift(:) = (MAT(:,size(MAT,2)) - MAT(:,1))/(size(MAT,2) - 1)
!                - TREND=+/-3 The least-squares lines from the time series are removed before
!                  time filtering.
!
!                IF TREND=-1,-2 or -3, the means, drifts or least-squares lines are reintroduced
!                post-filtering, respectively.
!                For other values of TREND nothing is done before or after filtering.
!
!   WIN          (INPUT, OPTIONAL) real(stnd)
!                By default, Hamming window filtering is used (i.e. WIN=0.54).
!                SET WIN=0.5 for Hanning window or WIN=1 for rectangular window.
!
!                WIN must be greater or equal to O.5 and less or equal to 1.
!
!   MAX_ALLOC    (INPUT, OPTIONAL) integer(i4b)
!                MAX_ALLOC is a factor which allows to reduce the workspace used to compute the 
!                Fourier transform of the data if necessary at the expense of increasing the computing time.
!                MAX_ALLOC must be greater or equal to 1 and less or equal to size(MAT,1).
!
!                The default is MAX_ALLOC= size(MAT,1).
!
!
! Further Details
! _______________
!
!   Use PL=0 for high-pass filtering frequencies corresponding to periods shorter than PH, 
!   or PH=0 for low-pass filtering frequencies corresponding to periods longer than PL.
!
!   Setting PH<PL is also allowed and performs band rejection of periods between PH and PL
!   (i.e. in that case the meaning of the PL and PH arguments are reversed).
!
!   Examples: 
!
!   For quarterly data:
!      call hwfilter( mat, pl=6, ph=32) returns components with periods between 1.5 and 8 yrs.
!
!   For monthly data:
!      call hwfilter( mat, pl=0, ph=24) returns components with all periods less than 2 yrs.
!
!   For more details and algorithm, see :
!
!   (1) Iacobucci, A., and Noullez, A., 2005:
!           A Frequency Selective Filter for Short-Length Time Series.
!           Computational Economics, 25,75-102.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use omp_lib
    use Utilities,         only : merror
    use Char_Constants,    only : tseries_error10, tseries_error11, tseries_error12,  &
                                  tseries_error13, tseries_error22, allocate_error
    use Reals_Constants,   only : one, zero, half
    use Logical_Constants, only : true, false
    use FFT_Procedures,    only : fft_row, init_fft, end_fft
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:,:) :: mat
    real(stnd),   intent(in), optional          :: win
!
    integer(i4b), intent(in)                :: pl, ph
    integer(i4b), intent(in), optional      :: trend, max_alloc
!
    logical(lgl),  intent(in), optional :: initfft
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, m, p, pl2, ph2, k, km, kodd, i, imax, i1, i2, trend2, max_alloc2
    integer      :: iok
!
    real(stnd)                                      :: comp, win2, hwin, c, d, e
    real(stnd), dimension(size(mat,1))              :: orig, slope
    !dir$ attributes align : 64 :: orig
    !dir$ attributes align : 64 :: slope
    real(stnd), dimension((size(mat,2)/2)+1)        :: h2
    !dir$ attributes align : 64 :: h2
    real(stnd), dimension(size(mat,2))              :: h
    !dir$ attributes align : 64 :: h
!
    complex(stnd), dimension(:,:), allocatable  :: matt
    !dir$ attributes align : 64 :: matt
!
    logical(lgl) :: odd, initfft2
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='hwfilter'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n =  size( mat, 1 )
    if ( n<=0_i4b )  return
!
    m = size( mat, 2 )
!
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    km = m/2_i4b
!
    if ( m==2_i4b*km ) then
        odd  = false
        kodd = km
    else
        odd  = true
        kodd = km + 1_i4b
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    max_alloc2 = n
!
    if ( present(max_alloc) ) then
        if ( max_alloc<1_i4b .or. max_alloc>n ) then
            call merror( name_proc//tseries_error22 )
        else
            max_alloc2 = max_alloc
        end if
    end if
!
!   CHECK PERIODS ORDERING TO SELECT FREQUENCIES.
!
    if ( ph<pl ) then
        pl2   = ph
        ph2   = pl
        comp  = one
    else
        pl2   = pl
        ph2   = ph
        comp  = zero
    end if
!
    if ( pl2<0_i4b)        &
    call merror( name_proc//tseries_error11 )
!
    if ( ph==1_i4b .or. pl==1_i4b )       &
    call merror( name_proc//tseries_error12 )
!
    if ( m<ph2)            &
    call merror( name_proc//tseries_error13 )
!
!   INITIALIZE THE FFT SUBROUTINE IF REQUIRED.
!
    if ( initfft2 ) then
        call init_fft( (/ n, m /), dim=2_i4b )
    end if
!
    imax = n/max_alloc2
!
    i = imax*max_alloc2
!
    if ( i==n - 1_i4b ) then
        max_alloc2 = max_alloc2 + 1_i4b
    else if ( i/=n ) then
        imax = imax + 1_i4b
    end if
!
!    if ( (imax*max_alloc2)/=n ) imax = imax + 1_i4b
!
!   ALLOCATE WORK SPACE FOR THE FOURIER TRANSFORM OF THE DATA.
!
    allocate( matt(max_alloc2,m), stat=iok )
!
    if ( iok/=0 ) then
        call merror( name_proc//allocate_error )
    end if
!
!   USE HAMMING WINDOW BY DEFAULT.
!
    win2 = 0.54_stnd
!
    if ( present(win) ) then
!
!       IF WINDOW IS NOT POSITIVE DEFINITE, USE HAMMING WINDOW.
!
        if ( win>=half .and. win<=one ) win2 = win
    end if    
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
!
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
!   BUILD UP IDEAL FILTER GAIN VECTOR.
!
    do k=0_i4b, km
!
!       SELECT FREQUENCY BAND.
!
        if ( k*ph2>=m .and. k*pl2<=m ) then
            h(k+1_i4b) = one - comp
        else
            h(k+1_i4b) = comp
        end if
!
    end do
!
!   CONVOLVE WITH WINDOW TRANSFORM TO GET FILTER FREQUENCY RESPONSE IF NEEDED.
!
    if ( win2/=one) then
!
        hwin = (one-win2)*half
        c    = h(2_i4b)
!
        if ( odd ) then
            e = h(km+1_i4b)
        else
            e = h(km)
        end if
!
        !dir$ assume_aligned h:64
        !$omp simd simdlen(8) linear(k:1)
        do k = 1_i4b, km
!
            d    = h(k)
            h(k) = win2*d + hwin*( c+h(k+1_i4b) )
            c    = d
!
        end do
!
        h(km+1_i4b) = win2*h(km+1_i4b) + hwin*(c+e)
!
    end if
!
    !dir$ assume_aligned h2:64
    !dir$ assume_aligned h:64
    h2(2_i4b:kodd)       = h(2_i4b:kodd)
    h(m:km+2_i4b:-1_i4b) = h2(2_i4b:kodd)
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    k = abs( trend2 )
!
    if ( k>=1_i4b .and. k<=3_i4b ) then
        call detrend_rm( mat(:n,:m), k, orig(:n), slope(:n) )
    end if
!
    i2 = 0_i4b
!
     !dir$ assume_aligned matt:64
    do i=1_i4b, imax
!
        i1 = i2 + 1_i4b
        i2 = min( n, i2+max_alloc2 )
        p  = i2 - i1 + 1_i4b
!
!       TRANSFORM THE REAL SEQUENCE.
!
        matt(:p,:m) = cmplx( mat(i1:i2,:m), zero, kind=stnd )
        call fft_row( matt(:p,:m), forward=true  )
!
!       NOW, FILTER THE TIME SERIES.
!       COMPUTE FILTERED FOURIER TRANSFORM OF THE SERIES.
!
        !dir$ assume_aligned h:64
        !$omp simd simdlen(8) reduction(*:matt)
        do k=1_i4b,m
            matt(:p,k) = h(k)*matt(:p,k)
        end do
!
!       INVERT THE SEQUENCE BACK.
!
        call fft_row( matt(:p,:m), forward=false )
!   
        mat(i1:i2,:m) = real( matt(:p,:m) )
!
    end do
!
!   DEALLOCATE WORK ARRAYS.
!
    deallocate( matt )
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!   ADD MEANS, DRIFTS OR LINEAR LEAST SQUARES LINES FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEANS.
!
            !dir$ assume_aligned mat:64
            !dir$ assume_aligned orig:64
            !$omp simd simdlen(8) reduction(+:mat) linear(k:1)
            do k = 1_i4b, m
                mat(:n,k) = mat(:n,k) + orig(:n)
            end do
!
        case( -2_i4b )
!
!           ADD DRIFTS.
!
            !dir$ assume_aligned mat:64
            !dir$ assume_aligned slope:64
            !$omp simd simdlen(8) reduction(+:mat) linear(k:1) 
            do k = 1_i4b, m
                c         = real( k, stnd ) - one
                mat(:n,k) = mat(:n,k) + c*slope(:n)           
            end do
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINES.
!         
            !dir$ assume_aligned mat:64
            !dir$ assume_aligned slope:64 
            !dir$ assume_aligned orig:64
            !$omp simd simdlen(8) reduction(+:mat) linear(k:1)
            do k = 1_i4b, m
                c         = real( k, stnd ) - one
                mat(:n,k) = mat(:n,k) + c*slope(:n) + orig(:n)      
            end do
!
    end select
!
!
! END OF SUBROUTINE hwfilter_rm
! _____________________________
!
    end subroutine hwfilter_rm
!
! =========================================================================================
!
    subroutine hwfilter2_rv( vec, pl, ph, trend, win )
       !dir$ attributes code_align : 32 :: hwfilter2_rv
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: hwfilter2_rv
!
! Purpose
! _______
!
!   Subroutine HWFILTER2 filters a time series (e.g. the argument VEC) in the frequency band
!   limited by periods PL and PH by windowed filtering (PL and PH are expressed in number
!   of points, i.e. PL=6(18) and PH=32(96) selects periods between 1.5 yrs and 8 yrs for
!   quarterly (monthly) data).
!
!
! Arguments
! _________
!
!   VEC          (INPUT/OUTPUT) real(stnd), dimension(:)
!                The time series vector to be filtered.
!                Size(VEC) must be greater or equal to 4.
!
!   PL           (INPUT) integer(i4b)
!                Minimum period of oscillation of desired component.
!                Use PL=0 for high-pass filtering frequencies corresponding
!                to periods shorter than PH, 
!                PL must be equal to 0 or greater or equal to 2.
!                Moreover, PL must be less or equal to size(VEC).
!
!   PH           (INPUT) integer(i4b)
!                Maximum period of oscillation of desired component.
!                USE PH=0 for low-pass filtering frequencies corresponding
!                to periods longer than PL.
!                PH must be equal to 0 or greater or equal to 2.
!                Moreover, PH must be less or equal to size(VEC).
!
!   TREND        (INPUT, OPTIONAL) integer(i4b)
!                If:
!
!                - TREND=+/-1 The mean of the time series is removed before time filtering
!                - TREND=+/-2 The drift from the time series is removed before time filtering
!                  by using the formula: drift = (VEC(size(VEC)) - VEC(1))/(size(VEC) - 1)
!                - TREND=+/-3 The least-squares line from the time series is removed before
!                  time filtering.
!
!                IF TREND=-1,-2 or -3, the mean, drift or least-squares line is reintroduced
!                post-filtering, respectively.
!                For other values of TREND nothing is done before or after filtering.
!
!   WIN          (INPUT, OPTIONAL) real(stnd)
!                By default, Hamming window filtering is used (i.e. WIN=0.54).
!                SET WIN=0.5 for Hanning window or WIN=1 for rectangular window.
!
!                WIN must be greater or equal to O.5 and less or equal to 1.
!
!
! Further Details
! _______________
!
!   Use PL=0 for high-pass filtering frequencies corresponding to periods shorter than PH, 
!   or PH=0 for low-pass filtering frequencies corresponding to periods longer than PL.
!
!   Setting PH<PL is also allowed and performs band rejection of periods between PH and PL
!   (i.e. in that case the meaning of the PL and PH arguments are reversed).
!
!   The unique difference between HWFILTER2 and HWFILTER is the use of the Goertzel method
!   for computing the Fourier transform of the data instead of a Fast Fourier Transform algorithm.
!
!   Examples: 
!
!   For quarterly data:
!      call hwfilter2( vec, pl=6, ph=32) returns component with periods between 1.5 and 8 yrs.
!
!   For monthly data:
!      call hwfilter2( vec, pl=0, ph=24) returns component with all periods less than 2 yrs.
!
!   For more details and algorithm, see :
!
!   (1) Iacobucci, A., and Noullez, A., 2005:
!           A Frequency Selective Filter for Short-Length Time Series.
!           Computational Economics, 25,75-102.
!
!   (2) Goertzel, G., 1958:
!           An Algorithm for the Evaluation of Finite Trigonometric Series.
!           The American Mathematical Monthly, Vol. 65, No. 1, pp. 34-35
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : omp_limit
    use omp_lib,            only : omp_get_num_procs, omp_get_max_threads, omp_in_parallel
    use Utilities,         only : merror, arth
    use Char_Constants,    only : tseries_error10, tseries_error11, tseries_error12,  &
                                  tseries_error13
    use Reals_Constants,   only : one, zero, half, two, four, pi
    use Logical_Constants, only : true, false
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:)  :: vec
    real(stnd),   intent(in), optional         :: win
!
    integer(i4b), intent(in)                :: pl, ph
    integer(i4b), intent(in), optional      :: trend
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, pl2, ph2, i, k, km, trend2
!
    real(stnd)   :: comp, win2, hwin, orig, slope, a, b, c, d, e,    &
                    ar, ai, br, bi
    real(stnd), dimension((size(vec)/2)+1)  :: h, vecr, veci
    !dir$ attributes align : 64 :: h
    !dir$ attributes align : 64 :: vecr
    !dir$ attributes align : 64 :: veci
!
    logical(lgl) :: odd
!

    logical      :: test_par

!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='hwfilter2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m = size( vec )
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
!   CHECK PERIODS ORDERING TO SELECT FREQUENCIES.
!
    if ( ph<pl ) then
        pl2   = ph
        ph2   = pl
        comp  = one
    else
        pl2   = pl
        ph2   = ph
        comp  = zero
    end if
!
    if ( pl2<0_i4b)        &
    call merror( name_proc//tseries_error11 )
!
    if ( pl==1_i4b .or. ph==1_i4b )        &
    call merror( name_proc//tseries_error12 )
!
    if ( m<ph2)            &
    call merror( name_proc//tseries_error13 )
!
!   USE HAMMING WINDOW BY DEFAULT.
!
    win2 = 0.54_stnd
!
    if ( present(win) ) then
!
!       IF WINDOW IS NOT POSITIVE DEFINITE, USE HAMMING WINDOW.
!
        if ( win>=half .and. win<=one ) win2 = win
    end if    
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
!
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
    km = m/2_i4b
!
    if ( m==2_i4b*km ) then
        odd = false
    else
        odd = true
    end if
!
!   BUILD UP IDEAL FILTER GAIN VECTOR.
!
    !dir$ assume_aligned h:64
    !$omp simd simdlen(8) linear(i:1)
    do i=0_i4b, km
!
!       SELECT FREQUENCY BAND.
!
        if ( i*ph2>=m .and. i*pl2<=m ) then
            h(i+1_i4b) = one - comp
        else
            h(i+1_i4b) = comp
        end if
!
    end do
!
!   CONVOLVE WITH WINDOW TRANSFORM TO GET FILTER FREQUENCY RESPONSE IF NEEDED.
!
    if ( win2/=one ) then
!
        hwin = (one-win2)*half
        a    = h(2_i4b)
!
        if ( odd ) then
            b = h(km+1_i4b)
        else
            b = h(km)
        end if
!
        !dir$ assume_aligned h:64
        !$omp simd simdlen(8) linear(k:1)
        do k = 1_i4b, km
!
            d    = h(k)
            h(k) = win2*d + hwin*( a+h(k+1_i4b) )
            a    = d
!
        end do
!
        h(km+1_i4b) = win2*h(km+1_i4b) + hwin*(a+b)
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    k = abs( trend2 )
!
    if ( k>=1_i4b .and. k<=3_i4b ) then
        call detrend_rv( vec(:m), k, orig, slope )
    end if
!   

    i = omp_get_num_procs()
    k = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               m>=omp_limit                 .and.      &
               i>1_i4b                      .and.      &
               k>1_i4b

    !dir$ assume_aligned h:64
    !dir$ assume_aligned vecr:64
    !dir$ assume_aligned veci:64
!
!$OMP PARALLEL IF(test_par)                             &
!$OMP         ,PRIVATE(i,k,a,b,c,d,e,ar,br,ai,bi)       &
!$OMP         ,SHARED(m,km,vec,vecr,veci,h,odd)
!
!   NOW, FILTER THE TIME SERIES.
!   COMPUTE FILTERED FOURIER TRANSFORM OF THE SERIES USING GOERTZEL METHOD.
!
    c = pi/real( m, stnd )
!
!$OMP DO SCHEDULE(STATIC,8) 
!
    do k = 0_i4b, km
!
        if ( h(k+1_i4b)==zero ) then
!
            vecr(k+1_i4b) = zero
            veci(k+1_i4b) = zero
!
        else
!
            e = c*real( k, stnd )
            d = -four*sin( e )**2
            a = zero
            b = zero
!
            !$omp simd simdlen(8) private(a) reduction(+:b)
            do i = m, 1_i4b, -1_i4b
!
                a = a + b
                b = b + d*a+ vec(i)
!
            end do
!
            vecr(k+1_i4b) = h(k+1_i4b)*(b-half*d*a)
            veci(k+1_i4b) = -h(k+1_i4b)*a*sin( two*e )
!
        end if
!
    end do
!
!$OMP END DO
!
!   SET ODD/EVEN NUMBER OF POINTS LAST SPECTRUM VALUE.
!
    if ( odd ) then
        b = vecr(km+1_i4b)
        a = veci(km+1_i4b)
    else
        b = half*vecr(km+1_i4b)
        a = zero
    end if
!
!   COMPUTE INVERSE FOURIER TRANSFORM TO GET BACK THE DATA.
!
!$OMP DO SCHEDULE(STATIC,8) 
!
    do i = 0_i4b, m-1_i4b
!
        e  = c*real( i, stnd )
        d  = -four*sin( e )**2
        ar = b
        br = b
        ai = a
        bi = a
!
        !$omp simd simdlen(8) reduction(+:br) reduction(+:ar) &
        !$omp&                reduction(+:bi) reduction(+:ai)
        do k = km, 2_i4b, -1_i4b
!
            br = br + d*ar + vecr(k)
            ar = ar + br
            bi = bi + d*ai + veci(k)
            ai = ai + bi
!
        end do
!
        vec(i+1_i4b) = two*(br-sin(two*e)*ai) + d*ar + vecr(1_i4b)
!
    end do
!
!$OMP END DO
!
!$OMP END PARALLEL
!
    vec(:m) = (one/real(m,stnd))*vec(:m)
!
!   ADD MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEAN.
!
            vec(:m) = vec(:m) + orig
!
        case( -2_i4b )
!
!           ADD DRIFT.
!
            vec(:m) = vec(:m) + slope*arth( zero, one, m )           
!
!            do k = 1_i4b, m
!                a      = real( k, stnd ) - one
!                vec(k) = vec(k) + a*slope           
!            end do
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINE.
!         
            vec(:m) = vec(:m) + ( orig + slope*arth( zero, one, m ) )       
!         
!            do k = 1_i4b, m
!                a      = real( k, stnd ) - one
!                vec(k) = vec(k) + a*slope + orig      
!            end do
!
    end select
!
!
! END OF SUBROUTINE hwfilter2_rv
! ______________________________
!
    end subroutine hwfilter2_rv
!
! =========================================================================================
!
    subroutine hwfilter2_rm( mat, pl, ph, trend, win )
       !dir$ attributes code_align : 32 :: hwfilter2_rm
       !dir$ optimize: 3
       !dir$ attributes optimization_parameter:TARGET_ARCH=skylake_avx512 :: hwfilter2_rm
!
! Purpose
! _______
!
!   Subroutine HWFILTER2 filters a multi-channel time series (e.g. the argument MAT) in the frequency
!   band limited by periods PL and PH by windowed filtering (PL and PH are expressed in number
!   of points, i.e. PL=6(18) and PH=32(96) selects periods between 1.5 yrs and 8 yrs for
!   quarterly (monthly) data).
!
!
! Arguments
! _________
!
!   MAT          (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                The multi-channel time series matrix to be filtered.
!                Each column of MAT corresponds to one observation.
!
!                Size(MAT,2) must be greater or equal to 4.
!
!   PL           (INPUT) integer(i4b)
!                Minimum period of oscillation of desired component.
!                Use PL=0 for high-pass filtering frequencies corresponding
!                to periods shorter than PH, 
!                PL must be equal to 0 or greater or equal to 2.
!                Moreover, PL must be less or equal to size(MAT,2).
!
!   PH           (INPUT) integer(i4b)
!                Maximum period of oscillation of desired component.
!                USE PH=0 for low-pass filtering frequencies corresponding
!                to periods longer than PL.
!                PH must be equal to 0 or greater or equal to 2.
!                Moreover, PH must be less or equal to size(MAT,2).
!
!   TREND        (INPUT, OPTIONAL) integer(i4b)
!                If:
!
!                - TREND=+/-1 The means of the time series are removed before time filtering
!                - TREND=+/-2 The drifts from the time series are removed before time filtering
!                  by using the formula: drift(:) = (MAT(:,size(MAT,2)) - MAT(:,1))/(size(MAT,2) - 1)
!                - TREND=+/-3 The least-squares lines from the time series are removed before
!                  time filtering.
!
!                IF TREND=-1,-2 or -3, the means, drifts or least-squares lines are reintroduced
!                post-filtering, respectively.
!                For other values of TREND nothing is done before or after filtering.
!
!   WIN          (INPUT, OPTIONAL) real(stnd)
!                By default, Hamming window filtering is used (i.e. WIN=0.54).
!                SET WIN=0.5 for Hanning window or WIN=1 for rectangular window.
!
!                WIN must be greater or equal to O.5 and less or equal to 1.
!
!
! Further Details
! _______________
!
!   Use PL=0 for high-pass filtering frequencies corresponding to periods shorter than PH, 
!   or PH=0 for low-pass filtering frequencies corresponding to periods longer than PL.
!
!   Setting PH<PL is also allowed and performs band rejection of periods between PH and PL
!   (i.e. in that case the meaning of the PL and PH arguments are reversed).
!
!   The unique difference between HWFILTER2 and HWFILTER is the use of the Goertzel method
!   for computing the Fourier transform of the data instead of a Fast Fourier Transform algorithm.
!
!   Examples: 
!
!   For quarterly data:
!      call hwfilter2( mat, pl=6, ph=32) returns components with periods between 1.5 and 8 yrs.
!
!   For monthly data:
!      call hwfilter2( mat, pl=0, ph=24) returns components with all periods less than 2 yrs.
!
!   For more details and algorithm, see :
!
!   (1) Iacobucci, A., and Noullez, A., 2005:
!           A Frequency Selective Filter for Short-Length Time Series.
!           Computational Economics, 25,75-102.
!
!   (2) Goertzel, G., 1958:
!           An Algorithm for the Evaluation of Finite Trigonometric Series.
!           The American Mathematical Monthly, Vol. 65, No. 1, pp. 34-35
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Select_Parameters,  only : lgl, i4b, stnd
    use Select_Parameters,  only : omp_limit
    use omp_lib,            only : omp_get_num_procs, omp_get_max_threads, omp_in_parallel
    use Utilities,         only : merror
    use Char_Constants,    only : tseries_error10, tseries_error11, tseries_error12,  &
                                  tseries_error13
    use Reals_Constants,   only : one, zero, half, two, four, pi
    use Logical_Constants, only : true, false
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:,:) :: mat
    real(stnd),   intent(in), optional          :: win
!
    integer(i4b), intent(in)                :: pl, ph
    integer(i4b), intent(in), optional      :: trend
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, m, pl2, ph2, i, k, km, trend2
!
    real(stnd)                                            :: comp, win2, hwin, c, d, e
    real(stnd), dimension(size(mat,1))                    :: orig, slope, a, b, ar, ai, br, bi
    !dir$ attributes align : 64 :: orig
    !dir$ attributes align : 64 :: slope
    !dir$ attributes align : 64 :: a
    !dir$ attributes align : 64 :: b
    !dir$ attributes align : 64 :: ar
    !dir$ attributes align : 64 :: ai
    !dir$ attributes align : 64 :: br
    !dir$ attributes align : 64 :: bi
    real(stnd), dimension((size(mat,2)/2)+1)              :: h
    !dir$ attributes align : 64 :: h
    real(stnd), dimension(size(mat,1),(size(mat,2)/2)+1)  :: matr, mati
    !dir$ attributes align : 64 :: matr
    !dir$ attributes align : 64 :: mati
!
    logical(lgl) :: odd
!

    logical      :: test_par

!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='hwfilter2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n =  size( mat, 1 )
    if ( n<=0_i4b )  return
!
    m = size( mat, 2 )
!
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
!   CHECK PERIODS ORDERING TO SELECT FREQUENCIES.
!
    if ( ph<pl ) then
        pl2   = ph
        ph2   = pl
        comp  = one
    else
        pl2   = pl
        ph2   = ph
        comp  = zero
    end if
!
    if ( pl2<0_i4b)        &
    call merror( name_proc//tseries_error11 )
!
    if ( pl==1_i4b .or. ph==1_i4b )        &
    call merror( name_proc//tseries_error12 )
!
    if ( m<ph2)            &
    call merror( name_proc//tseries_error13 )
!
!   USE HAMMING WINDOW BY DEFAULT.
!
    win2 = 0.54_stnd
!
    if ( present(win) ) then
!
!       IF WINDOW IS NOT POSITIVE DEFINITE, USE HAMMING WINDOW.
!
        if ( win>=half .and. win<=one ) win2 = win
    end if    
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
!
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
    km = m/2_i4b
!
    if ( m==2_i4b*km ) then
        odd = false
    else
        odd = true
    end if
!
!   BUILD UP IDEAL FILTER GAIN VECTOR.
!
    !dir$ assume_aligned h:64
     !$omp simd simdlen(8) linear(i:1)
    do i=0_i4b, km
!
!       SELECT FREQUENCY BAND.
!
        if ( i*ph2>=m .and. i*pl2<=m ) then
            h(i+1_i4b) = one - comp
        else
            h(i+1_i4b) = comp
        end if
!
    end do
!
!   CONVOLVE WITH WINDOW TRANSFORM TO GET FILTER FREQUENCY RESPONSE IF NEEDED.
!
    if ( win2/=one) then
!
        hwin = (one-win2)*half
        c    = h(2_i4b)
!
        if ( odd ) then
            e = h(km+1_i4b)
        else
            e = h(km)
        end if
!
        !dir$ assume_aligned h:64
        !$omp simd simdlen(8) private(d,c) linear(k:1)
        do k = 1_i4b, km
!
            d    = h(k)
            h(k) = win2*d + hwin*( c+h(k+1_i4b) )
            c    = d
!
        end do
!
        h(km+1_i4b) = win2*h(km+1_i4b) + hwin*(c+e)
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    k = abs( trend2 )
!
    if ( k>=1_i4b .and. k<=3_i4b ) then
        call detrend_rm( mat(:n,:m), k, orig(:n), slope(:n) )
    end if
!   

    i = omp_get_num_procs()
    k = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               m*n>=omp_limit               .and.      &
               i>1_i4b                      .and.      &
               k>1_i4b

!
!$OMP PARALLEL IF(test_par)                             &
!$OMP         ,PRIVATE(i,k,a,b,c,d,e,ar,br,ai,bi)       &
!$OMP         ,SHARED(n,m,km,mat,matr,mati,h,odd)
!
!   NOW, FILTER THE TIME SERIES.
!   COMPUTE FILTERED FOURIER TRANSFORM OF THE SERIES USING GOERTZEL METHOD.
!
    c = pi/real( m, stnd )
!
!$OMP DO SCHEDULE(STATIC,8) 
!
    !dir$ assume_aligned matr:64
    !dir$ assume_aligned mati:64
    !dir$ assume_aligned a:64
    !dir$ assume_aligned b:64
    do k = 0_i4b, km
!
        if ( h(k+1_i4b)==zero ) then
!
            matr(:n,k+1_i4b) = zero
            mati(:n,k+1_i4b) = zero
!
        else
!
            e     = c*real( k, stnd )
            d     = -four*sin( e )**2
            a(:n) = zero
            b(:n) = zero
!
            !$omp simd simdlen(8) reduction(+:a) reduction(+:b)
            do i = m, 1_i4b, -1_i4b
!
                a(:n) = a(:n) + b(:n)
                b(:n) = b(:n) + d*a(:n)+ mat(:n,i)
!
            end do
!
            matr(:n,k+1_i4b) = h(k+1_i4b)*( b(:n)-half*d*a(:n) )
            mati(:n,k+1_i4b) = -h(k+1_i4b)*sin( two*e )*a(:n)
!
        end if
!
    end do
!
!$OMP END DO
!
!   SET ODD/EVEN NUMBER OF POINTS LAST SPECTRUM VALUE.
!
    if ( odd ) then
        b(:n) = matr(:n,km+1_i4b)
        a(:n) = mati(:n,km+1_i4b)
    else
        b(:n) = half*matr(:n,km+1_i4b)
        a(:n) = zero
    end if
!
!   COMPUTE INVERSE FOURIER TRANSFORM TO GET BACK THE FILTERED DATA.
!
!$OMP DO SCHEDULE(STATIC,8) 
!
    !dir$ assume_aligned ar:64
    !dir$ assume_aligned br:64
    !dir$ assume_aligned ai:64
    !dir$ assume_aligned bi:64
    do i = 0_i4b, m-1_i4b
!
        e  = c*real( i, stnd )
        d  = -four*sin( e )**2
        ar(:n) = b(:n)
        br(:n) = b(:n)
        ai(:n) = a(:n)
        bi(:n) = a(:n)
!
        !$omp simd simdlen(8) reduction(+:br) reduction(+:ar) &
        !$omp&                reduction(+:bi) reduction(+:ai)
        do k = km, 2_i4b, -1_i4b
!
            br(:n) = br(:n) + d*ar(:n) + matr(:n,k)
            ar(:n) = ar(:n) + br(:n)
            bi(:n) = bi(:n) + d*ai(:n) + mati(:n,k)
            ai(:n) = ai(:n) + bi(:n)
!
        end do
!
        mat(:n,i+1_i4b) = two*(br(:n)-sin(two*e)*ai(:n)) + d*ar(:n) + matr(:n,1_i4b)
!
    end do
!
!$OMP END DO
!
!$OMP END PARALLEL
!
    mat(:n,:m) = ( one/real(m,stnd) )*mat(:n,:m)
!
!   ADD MEANS, DRIFTS OR LINEAR LEAST SQUARES LINES FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEANS.
!
            !dir$ assume_aligned mat:64
            !dir$ assume_aligned orig:64
            !$omp simd simdlen(8) reduction(+:mat)
            do k = 1_i4b, m
                mat(:n,k) = mat(:n,k) + orig(:n)
            end do
!
        case( -2_i4b )
!
!           ADD DRIFTS.
!
            !dir$ assume_aligned mat:64
            !dir$ assume_aligned slope:64
            !$omp simd simdlen(8) private(c) reduction(+:mat)
            do k = 1_i4b, m
                c         = real( k, stnd ) - one
                mat(:n,k) = mat(:n,k) + c*slope(:n)           
            end do
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINES.
!         
            !dir$ assume_aligned mat:64
            !dir$ assume_aligned slope:64
            !dir$ assume_aligned orig:64
            !$omp simd simdlen(8) private(c) reduction(+:mat)
            do k = 1_i4b, m
                c         = real( k, stnd ) - one
                mat(:n,k) = mat(:n,k) + c*slope(:n) + orig(:n)      
            end do
!
    end select
!
!
! END OF SUBROUTINE hwfilter2_rm
! ______________________________
!
    end subroutine hwfilter2_rm
!
! =========================================================================================
!
    function lp_coef( pl, k, fc, notest_fc ) result( lpcoef )
!
! Purpose
! _______
!
!   Function LP_COEF computes the K-term least squares approximation to an 
!   -ideal- low pass filter with cutoff period PL (e.g. cutoff frequency FC =  1/PL).
!   
!   This filter has a transfer function with a transition band of width delta surrounding FC,
!   where delta = 4 * pi/K when FC is expressed in radians.
! 
! 
! Arguments
! _________
!
!   PL          (INPUT) integer(i4b)
!               Minimum period of oscillation of desired component. The corresponding
!               cutoff frequency is FC=1/PL (i.e. filter has zero response in the interval
!               [FC+1/K,Nyquist] and one response in the  interval [0,FC-1/K].
!
!               PL must be greater than 2 and FC must verify the following inequalities:
!
!               - FC - 1/K >= 0
!               - FC + 1/K <  0.5
!
!   K           (INPUT) integer(i4b)
!               The number of filter terms to be computed.
!               K must be greater or equal to 3 and odd.
!
!   FC          (INPUT, OPTIONAL) real(stnd)
!               The user chosen cutoff frequency in cycles per sample interval.
!               If the optional argument FC is used, the PL argument is not used
!               to determine the cutoff frequency.
!
!               FC must verify the following inequalities:
!
!               - FC - 1/K >= 0
!               - FC + 1/K <  0.5
!
!   NOTEST_FC   (INPUT, OPTIONAL) logical(lgl)
!               On input, if this optional logical argument is set to true, the two tests on the
!               cutoff frequency (e.g. FC - 1/K >= 0 and FC + 1/K <  0.5) are bypassed. However, in
!               that case, the cutoff frequency FC must still verify the inequalities 0 < FC < 0.5.
!
!
! Further Details
! _______________
!
!   Function LP_COEF computes symmetric linear low-pass  filter coefficients
!   using a  least squares  approximation  to an  ideal low-pass  filter with
!   convergence factors  (i.e. Lanczos window) which reduce  overshoot and  
!   ripple  (Bloomfield, 1976).
!
!   This low-pass filter has a transfer function which changes from approximately
!   one to zero in a transition band about the ideal cutoff frequency FC (FC=1/PL), 
!   that is from (FC - 1/K) to (FC + 1/K),  as discussed  in  section  6.4  of  Bloomfield
!   (1976).  The  user  must  specify  the  cutoff period (or the cutoff frequency) and the number
!   of filter coefficients, which must be odd.
!
!   The  user must also choose the number of filter terms, K, so that 
!   (FC - 1/K) >= 0 and (FC + 1/K) < 0.5 if the optional logical argument NOTEST_FC is not used
!   or is not set to true.
!
!   In addition, K must be chosen as a compromise between:
!
!      1)  A sharp cutoff, that is, 1/K small; and
!
!      2)  Minimizing the number of data points lost by the filtering operations
!          (e.g. (K-1)/2 data points will be lost from each end of the series).
!
!   The subroutine returns the normalized low-pass filter coefficients.   
!
!   This function is adapted from the STARPAC software developed by the National Institute
!   of Standards and Technology (NIST). For more details and algorithm, see
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Logical_Constants, only : true
    use Reals_Constants,   only : zero, half, one, two, pi
    use Char_Constants,    only : tseries_error26, tseries_error27, tseries_error28, tseries_error31,     &
                                  tseries_error33, tseries_error35
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: pl, k
!
    real(stnd), intent(in), optional :: fc
!
    real(stnd), dimension(k) :: lpcoef
!
    logical(lgl), intent(in), optional :: notest_fc
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)               :: fc2, sum, arg, con, con2, tmp
    real(stnd), dimension(k) :: coef
!
    integer(i4b) :: i, kmid, khalf
!
    logical(lgl) :: test_fc
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='lp_coef'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    if ( k<3_i4b )              &
    call merror( name_proc//tseries_error26 )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error27 )
!
    if ( present(fc) ) then
!
        if ( fc<=zero .or. fc>=half )      &
        call merror( name_proc//tseries_error28 )
!
        fc2  = fc
!
    else
!
        if ( pl<=2_i4b )            &
        call merror( name_proc//tseries_error31 )
!
        fc2  = one/real( pl, stnd )
!
    end if
!
    test_fc = true
!
    if ( present( notest_fc ) ) then
        test_fc = .not.notest_fc
    end if
!
!   BYPASS THE TESTS ON THE CUTOFF FREQUENCY IF REQUIRED.
!
    if ( test_fc ) then
!
        tmp = one/real( k, stnd )
!
        if ( (fc2 - tmp) < zero  )            &
        call merror( name_proc//tseries_error33 )
!
        if ( (fc2 + tmp) >= half )            &
        call merror( name_proc//tseries_error35 )
!
    end if
!
!   CALCULATE THE LOW-PASS FILTER COEFFICIENTS.
!
    kmid  = ( k + 1_i4b )/2_i4b
    khalf = ( k - 1_i4b )/2_i4b
!
    con  = pi/real( kmid, stnd )
    con2 = fc2*two*pi
!
    sum  = two*fc2
    arg  = zero
!
    coef(kmid) = sum
!
    do i=1_i4b, khalf
!
        arg = arg + con
        tmp = sin( real(i,stnd)*con2 )*sin( arg )/( real(i,stnd)*pi*arg )
!
        coef(kmid-i) = tmp
        coef(kmid+i) = tmp
!
        sum = sum + tmp + tmp
!
    end do
!
!   NORMALIZE THE LOW-PASS COEFFICIENTS.
!
    lpcoef(:k) = coef(:k)/sum
!
!
! END OF FUNCTION lp_coef
! _______________________
!
    end function lp_coef
!
! =========================================================================================
!
    function lp_coef2( pl, k, fc, win, notest_fc  ) result( lpcoef2 )
!
! Purpose
! _______
!
!   Function LP_COEF2 computes the K-term least squares approximation to an 
!   -ideal- low pass filter with cutoff period PL (e.g. cutoff frequency FC =  1/PL)
!   by windowed filtering (e.g. Hamming window is used).
!   
!   This filter has a transfer function with a transition band of width delta surrounding FC,
!   where delta = 4 * pi/K when FC is expressed in radians.
!   
! 
! Arguments
! _________
!
!   PL          (INPUT) integer(i4b)
!               Minimum period of oscillation of desired component. The corresponding
!               cutoff frequency is FC=1/PL (i.e. filter has zero response in the interval
!               [FC+1/K,Nyquist] and one response in the interval [0,FC-1/K].

!               PL must be greater than two and FC must verify the following inequalities:
!
!               - FC - 1/K >= 0
!               - FC + 1/K <  0.5
!
!   K           (INPUT) integer(i4b)
!               The number of filter terms to be computed.
!               K must be greater or equal to 3 and odd.
!
!   FC          (INPUT, OPTIONAL) real(stnd)
!               The user chosen cutoff frequency in cycles per sample interval.
!               If the optional argument FC is used, the PL argument is not used
!               to determine the cutoff frequency.
!
!               FC must verify the following inequalities:
!
!               - FC - 1/K >= 0
!               - FC + 1/K <  0.5
!
!   WIN         (INPUT, OPTIONAL) real(stnd)
!               By default, Hamming window filtering is used (i.e. WIN=0.54).
!               Set WIN=0.5 for Hanning window or WIN=1 for rectangular window.
!
!               WIN must be greater or equal to O.5 and less or equal to 1, otherwise
!               WIN is reset to 0.54.
!
!   NOTEST_FC   (INPUT, OPTIONAL) logical(lgl)
!               On input, if this optional logical argument is set to true, the two tests on the
!               cutoff frequency (e.g. FC - 1/K >= 0 and FC + 1/K <  0.5) are bypassed. However, in
!               that case, the cutoff frequency FC must still verify the inequalities 0 < FC < 0.5.
!
!
! Further Details
! _______________
!
!   Function LP_COEF2  computes symmetric linear low-pass filter coefficients
!   using a  least squares  approximation  to an  ideal low-pass  filter. The Hamming
!   window is used to reduce overshoot and ripple in the transfert function of the 
!   ideal low-pass filter.
!
!   This low-pass filter has a transfer function which changes from approximately
!   one to zero in a transition band about the ideal cutoff frequency FC (FC=1/PL), 
!   that is from (FC - 1/K) to (FC + 1/K),  as discussed  in  section  6.4  of  Bloomfield
!   (1976).  The  user  must  specify  the  cutoff period (or the cutoff frequency) and the number
!   of filter coefficients, which must be odd.
!
!   The  user must also choose the number of filter terms, K, so that 
!   (FC - 1/K) >= 0 and (FC + 1/K) < 0.5 if the optional logical argument NOTEST_FC is not used
!   or is not set to true.
!
!   The overshoot and the associated ripples in the ideal transfert function are reduced by the use
!   of the Hamming window.
!
!   In addition, K must be chosen as a compromise between:
!
!      1)  A sharp cutoff, that is, 1/K small; and
!
!      2)  Minimizing the number of data points lost by the filtering operations
!          (e.g. (K-1)/2 data points will be lost from each end of the series).
!
!   The subroutine returns the normalized low-pass filter coefficients.   
!
!   For more details and algorithm, see
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Logical_Constants, only : true
    use Reals_Constants,   only : zero, half, one, two, pi
    use Char_Constants,    only : tseries_error26, tseries_error27, tseries_error28, tseries_error31,     &
                                  tseries_error33, tseries_error35
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: pl, k
!
    real(stnd), intent(in), optional :: fc, win
!
    real(stnd), dimension(k) :: lpcoef2
!
    logical(lgl), intent(in), optional :: notest_fc
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)               :: fc2, sum, arg, con, con2, tmp, win2
    real(stnd), dimension(k) :: coef
!
    integer(i4b) :: i, kmid, khalf
!
    logical(lgl) :: test_fc
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='lp_coef2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    if ( k<3_i4b )              &
    call merror( name_proc//tseries_error26 )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error27 )
!
    if ( present(fc) ) then
!
        if ( fc<=zero .or. fc>=half )      &
        call merror( name_proc//tseries_error28 )
!
        fc2  = fc
!
    else
!
        if ( pl<=2_i4b )            &
        call merror( name_proc//tseries_error31 )
!
        fc2  = one/real( pl, stnd )
!
    end if
!
    test_fc = true
!
    if ( present( notest_fc ) ) then
        test_fc = .not.notest_fc
    end if
!
!   BYPASS THE TESTS ON THE CUTOFF FREQUENCY IF REQUIRED.
!
    if ( test_fc ) then
!
        tmp = one/real( k, stnd )
!
        if ( (fc2 - tmp) < zero  )            &
        call merror( name_proc//tseries_error33 )
!
        if ( (fc2 + tmp) >= half )            &
        call merror( name_proc//tseries_error35 )
!
    end if
!
!   USE HAMMING WINDOW BY DEFAULT.
!
    win2 = 0.54_stnd
!
    if ( present(win) ) then
!
!       IF WINDOW IS NOT POSITIVE DEFINITE, USE HAMMING WINDOW.
!
        if ( win>=half .and. win<=one ) win2 = win
!
    end if    
!
!   CALCULATE THE LOW-PASS FILTER COEFFICIENTS.
!
    kmid  = ( k + 1_i4b )/2_i4b
    khalf = ( k - 1_i4b )/2_i4b
!
    con  = pi/real( khalf, stnd )
    con2 = fc2*two*pi
!
    sum        = two*fc2
    coef(kmid) = sum
!
    do i=1_i4b, khalf
!
        arg  = real( i, stnd )
        tmp = win2 + (one-win2)*cos( arg*con )
        tmp = tmp*( sin( arg*con2 )/( arg*pi ) )
!
        coef(kmid-i) = tmp
        coef(kmid+i) = tmp
!
        sum = sum + tmp + tmp
!
    end do
!
!   NORMALIZE THE LOW-PASS COEFFICIENTS.
!
    lpcoef2(:k) = coef(:k)/sum
!
!
! END OF FUNCTION lp_coef2
! ________________________
!
    end function lp_coef2
!
! =========================================================================================
!
    function hp_coef( ph, k, fc, notest_fc ) result( hpcoef )
!
! Purpose
! _______
!
!   Function HP_COEF computes the K-term least squares approximation to an 
!   -ideal- high pass filter with cutoff period PH (e.g. cutoff frequency FC =  1/PH).
!
!   This filter has a transfer function with a transition band of width delta surrounding FC,
!   where delta = 4 * pi/K when FC is expressed in radians.
!
! Arguments
! _________
!
!   PH           (INPUT) integer(i4b)
!                Maximum period of oscillation of desired component. The corresponding
!                cutoff frequency is FC=1/PH (i.e. filter has one response in the interval
!                [FC+1/K,Nyquist] and zero response in the interval [0,FC-1/K].
!
!                PH must be greater than two and FC must verify the following inequalities:
!
!                - FC - 1/K >= 0
!                - FC + 1/K <  0.5
!
!   K            (INPUT) integer(i4b)
!                The number of filter terms to be computed.
!                K must be greater or equal to 3 and odd.
!
!   FC           (INPUT, OPTIONAL) real(stnd)
!                The user chosen cutoff frequency in cycles per sample interval.
!                If the optional argument FC is used, the PH argument is not used
!                to determine the cutoff frequency.
!
!                FC must verify the following inequalities:
!
!                - FC - 1/K >= 0
!                - FC + 1/K <  0.5
!
!   NOTEST_FC    (INPUT, OPTIONAL) logical(lgl)
!                On input, if this optional logical argument is set to true, the two tests on the
!                cutoff frequency (e.g. FC - 1/K >= 0 and FC + 1/K <  0.5) are bypassed. However, in
!                that case, the cutoff frequency FC must still verify the inequalities 0 < FC < 0.5.
!
!
! Further Details
! _______________
!
!   Function HP_COEF computes symmetric linear high-pass filter coefficients
!   from the corresponding low-pass filter as given by function LP_COEF. This is equivalent
!   to subtracting  the low-pass filtered series from the original time series.
!
!   This high-pass filter has a transfer function which changes from approximately
!   zero to one in a transition band about the ideal cutoff frequency FC (FC=1/PH), 
!   that is from (FC - 1/K) to (FC + 1/K),  as discussed  in  section  6.4  of  Bloomfield
!   (1976).  The  user  must  specify  the  cutoff period (or the cutoff frequency) and the number
!   of filter coefficients, which must be odd.
!
!   The  user must also choose the number of filter terms, K, so that (FC - 1/K) >= 0 
!   and (FC + 1/K) < 0.5 if the optional logical argument NOTEST_FC is not used
!   or is not set to true.
!
!   In addition, K must be chosen as a compromise between:
!
!      1)  A sharp cutoff, that is, 1/K small; and
!
!      2)  Minimizing the number of data points lost by the filtering operations
!          (e.g. (K-1)/2 data points will be lost from each end of the series).
!
!   The subroutine returns the high-pass filter coefficients.   
!
!   This function is adapted from the STARPAC software developed by the National Institute
!   of Standards and Technology (NIST).
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Logical_Constants, only : true
    use Reals_Constants,   only : zero, half, one, two, pi
    use Char_Constants,    only : tseries_error26, tseries_error27, tseries_error28, tseries_error32,   &
                                  tseries_error33, tseries_error35
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: ph, k
!
    real(stnd), intent(in), optional :: fc
!
    real(stnd), dimension(k) :: hpcoef
!
    logical(lgl), intent(in), optional :: notest_fc
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)               :: fc2, sum, arg, con, con2, tmp
    real(stnd), dimension(k) :: coef
!
    integer(i4b) :: i, kmid, khalf
!
    logical(lgl) :: test_fc
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='hp_coef'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    if ( k<3_i4b )              &
    call merror( name_proc//tseries_error26  )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error27 )
!
    if ( present(fc) ) then
!
        if ( fc<=zero .or. fc>=half )            &
        call merror( name_proc//tseries_error28 )
!
        fc2  = fc
!
    else
!
        if ( ph<=2_i4b )            &
        call merror( name_proc//tseries_error32 )
!
        fc2  = one/real( ph, stnd )
!
    end if
!
    test_fc = true
!
    if ( present( notest_fc ) ) then
        test_fc = .not.notest_fc
    end if
!
!   BYPASS THE TESTS ON THE CUTOFF FREQUENCY IF REQUIRED.
!
    if ( test_fc ) then
!
        tmp = one/real( k, stnd )
!
        if ( (fc2 - tmp) < zero  )            &
        call merror( name_proc//tseries_error33 )
!
        if ( (fc2 + tmp) >= half )            &
        call merror( name_proc//tseries_error35 )
!
    end if
!
!   CALCULATE THE LOW-PASS FILTER COEFFICIENTS.
!
    kmid  = ( k + 1_i4b )/2_i4b
    khalf = ( k - 1_i4b )/2_i4b
!
    con  = pi/real( kmid, stnd )
    con2 = fc2*two*pi
!
    sum  = two*fc2
    arg  = zero
!
    coef(kmid) = sum
!
    do i=1_i4b, khalf
!
        arg = arg + con
        tmp = sin( real(i,stnd)*con2 )*sin( arg )/( real(i,stnd)*pi*arg )
!
        coef(kmid-i) = tmp
        coef(kmid+i) = tmp
!
        sum = sum + tmp + tmp
!
    end do
!
!   NORMALIZE THE COEFFICIENTS.
!
    hpcoef(:k) = -coef(:k)/sum
!
    hpcoef(kmid) = hpcoef(kmid) + one
!
!
! END OF FUNCTION hp_coef
! _______________________
!
    end function hp_coef
!
! =========================================================================================
!
    function hp_coef2( ph, k, fc, win, notest_fc ) result( hpcoef2 )
!
! Purpose
! _______
!
!   Function HP_COEF2 computes the K-term least squares approximation to an 
!   -ideal- high pass filter with cutoff period PH (e.g. cutoff frequency FC =  1/PH)
!   by windowed filtering (e.g. Hamming window is used).
!
!   This filter has a transfer function with a transition band of width delta surrounding FC,
!   where delta = 4 * pi/K when FC is expressed in radians.
!
! Arguments
! _________
!
!   PH           (INPUT) integer(i4b)
!                Maximum period of oscillation of desired component. The corresponding
!                cutoff frequency is FC=1/PH (i.e. filter has one response in the interval
!                [FC+1/K,Nyquist] and zero response in the interval [0,FC-1/K].
!
!                PH must be greater than two and FC must verify the following inequalities:
!
!                - FC - 1/K >= 0
!                - FC + 1/K <  0.5
!
!   K            (INPUT) integer(i4b)
!                The number of filter terms to be computed.
!                K must be greater or equal to 3 and odd.
!
!   FC           (INPUT, OPTIONAL) real(stnd)
!                The user chosen cutoff frequency in cycles per sample interval.
!                If the optional argument FC is used, the PH argument is not used
!                to determine the cutoff frequency.
!
!                FC must verify the following inequalities:
!
!                - FC - 1/K >= 0
!                - FC + 1/K <  0.5
!
!   WIN          (INPUT, OPTIONAL) real(stnd)
!                By default, Hamming window filtering is used (i.e. WIN=0.54).
!                Set WIN=0.5 for Hanning window or WIN=1 for rectangular window.
!
!                WIN must be greater or equal to O.5 and less or equal to 1, otherwise
!                WIN is reset to 0.54.
!
!   NOTEST_FC    (INPUT, OPTIONAL) logical(lgl)
!                On input, if this optional logical argument is set to true, the two tests on the
!                cutoff frequency (e.g. FC - 1/K >= 0 and FC + 1/K <  0.5) are bypassed. However, in
!                that case, the cutoff frequency FC must still verify the inequalities 0 < FC < 0.5.
!
!
! Further Details
! _______________
!
!   Function HP_COEF2 computes symmetric linear high-pass filter coefficients
!   from the corresponding low-pass filter as given by function LP_COEF2. This is equivalent
!   to subtracting  the low-pass filtered series from the original time series.
!
!   This high-pass filter has a transfer function which changes from approximately
!   zero to one in a transition band about the ideal cutoff frequency FC (FC=1/PH), 
!   that is from (FC - 1/K) to (FC + 1/K),  as discussed  in  section  6.4  of  Bloomfield
!   (1976).  The  user  must  specify  the  cutoff period (or the cutoff frequency) and the number
!   of filter coefficients, which must be odd.
!
!   The  user must also choose the number of filter terms, K, so that (FC - 1/K) >= 0 
!   and (FC + 1/K) < 0.5 if the optional logical argument NOTEST_FC is not used
!   or is not set to true.
!
!   The overshoot and the associated ripples in the ideal transfert function are reduced by the use
!   of the Hamming window.
!
!   In addition, K must be chosen as a compromise between:
!
!      1)  A sharp cutoff, that is, 1/K small; and
!
!      2)  Minimizing the number of data points lost by the filtering operations
!          (e.g. (K-1)/2 data points will be lost from each end of the series).
!
!   The subroutine returns the high-pass filter coefficients.   
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Logical_Constants, only : true
    use Reals_Constants,   only : zero, half, one, two, pi
    use Char_Constants,    only : tseries_error26, tseries_error27, tseries_error28, tseries_error32,   &
                                  tseries_error33, tseries_error35
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: ph, k
!
    real(stnd), intent(in), optional :: fc, win
!
    real(stnd), dimension(k) :: hpcoef2
!
    logical(lgl), intent(in), optional :: notest_fc
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)               :: fc2, sum, arg, con, con2, tmp, win2
    real(stnd), dimension(k) :: coef
!
    integer(i4b) :: i, kmid, khalf
!
    logical(lgl) :: test_fc
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='hp_coef2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    if ( k<3_i4b )              &
    call merror( name_proc//tseries_error26  )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error27 )
!
    if ( present(fc) ) then
!
        if ( fc<=zero .or. fc>=half )            &
        call merror( name_proc//tseries_error28 )
!
        fc2  = fc
!
    else
!
        if ( ph<=2_i4b )            &
        call merror( name_proc//tseries_error32 )
!
        fc2  = one/real( ph, stnd )
!
    end if
!
    test_fc = true
!
    if ( present( notest_fc ) ) then
        test_fc = .not.notest_fc
    end if
!
!   BYPASS THE TESTS ON THE CUTOFF FREQUENCY IF REQUIRED.
!
    if ( test_fc ) then
!
        tmp = one/real( k, stnd )
!
        if ( (fc2 - tmp) < zero  )            &
        call merror( name_proc//tseries_error33 )
!
        if ( (fc2 + tmp) >= half )            &
        call merror( name_proc//tseries_error35 )
!
    end if
!
!   USE HAMMING WINDOW BY DEFAULT.
!
    win2 = 0.54_stnd
!
    if ( present(win) ) then
!
!       IF WINDOW IS NOT POSITIVE DEFINITE, USE HAMMING WINDOW.
!
        if ( win>=half .and. win<=one ) win2 = win
!
    end if    
!
!   CALCULATE THE LOW-PASS FILTER COEFFICIENTS.
!
    kmid  = ( k + 1_i4b )/2_i4b
    khalf = ( k - 1_i4b )/2_i4b
!
    con  = pi/real( khalf, stnd )
    con2 = fc2*two*pi
!
    sum        = two*fc2
    coef(kmid) = sum
!
    do i=1_i4b, khalf
!
        arg  = real( i, stnd )
        tmp = win2 + (one-win2)*cos( arg*con )
        tmp = tmp*( sin( arg*con2 )/( arg*pi ) )
!
        coef(kmid-i) = tmp
        coef(kmid+i) = tmp
!
        sum = sum + tmp + tmp
!
    end do
!
!   NORMALIZE THE COEFFICIENTS.
!
    hpcoef2(:k) = -coef(:k)/sum
!
    hpcoef2(kmid) = hpcoef2(kmid) + one
!
!
! END OF FUNCTION hp_coef2
! ________________________
!
    end function hp_coef2
!
! =========================================================================================
!
    function bd_coef( pl, ph, k, fch, fcl, notest_fc ) result( bdcoef )
!
! Purpose
! _______
!
!   Function BD_COEF computes the K-term least squares approximation to an 
!   -ideal- band pass filter with cutoff periods PL and PH (e.g. cutoff frequencies 1/PL and 1/PH).
!   
!   PL and PH are expressed in number of points, i.e. PL=6(18) and PH=32(96) selects periods
!   between 1.5 yrs and 8 yrs for quarterly(monthly) data).
!
!   Alternatively, the user can directly specify the two cutoff frequencies, FCL and FCH, 
!   corresponding to PL and PH.
!   
!
! Arguments
! _________
!
!   PL          (INPUT) integer(i4b)
!               Minimum period of oscillation of desired component. The corresponding
!               cutoff frequency is 1/PL.
!               PL must be greater than two and must verify the following inequalities:
!
!               - 1/PH +  1.3/(K+1) <= 1/PL - 1.3/(K+1)
!               - 1/PL + 1/K        <  0.5
!
!   PH          (INPUT) integer(i4b)
!               Maximum period of oscillation of desired component. The corresponding
!               cutoff frequency is 1/PH.
!               PH must be greater than two and 1/PH must verify the following inequalities:
!
!               - 0                <= 1/PH - 1/K 
!               - 1/PH + 1.3/(K+1) <= 1/PL -  1.3/(K+1)
!
!   K           (INPUT) integer(i4b)
!               The number of filter terms to be computed.
!               K must be greater or equal to 3 and odd.
!
!   FCH         (INPUT, OPTIONAL) real(stnd)
!               The user chosen (low) cutoff frequency in cycles per sample interval.
!               If the optional argument FCH is used, the PH argument is not used to
!               determine the (low) cutoff frequency.
!
!               FCH must verify the following inequalities:
!
!               - 0                 <= FCH - 1/K 
!               - FCH +  1.3/(K+1)  <= FCL -  1.3/(K+1)
!
!   FCL         (INPUT, OPTIONAL) real(stnd)
!               The user chosen (high) cutoff frequency in cycles per sample interval.
!               If the optional argument FCL is used, the PL argument is not used to
!               determine the cutoff (high) frequency.
!
!               FCL must verify the following inequalities:
!
!               - FCH +  1.3/(K+1) <= FCL -  1.3/(K+1)
!               - FCL + 1/K        <  0.5
!
!   NOTEST_FC   (INPUT, OPTIONAL) logical(lgl)
!               On input, if this optional logical argument is set to true, the tests on the
!               cutoff frequencies FCH and FCL  (e.g. FCH - 1/K >= 0 and FCL + 1/K <  0.5) are bypassed.
!               However, in that case FCH and FCL must still verify the inequalities FCH > 0, FCL < 0.5
!               and  FCH +  1.3/(K+1) <= FCL -  1.3/(K+1) .
!
!
! Further Details
! _______________
!
!   Function BD_COEF computes symmetric linear band-pass  filter coefficients
!   using a  least squares  approximation  to an  ideal band-pass  filter  that has
!   convergence  factors  which reduce  overshoot and  ripple  (Bloomfield, 1976).
!
!   This band-pass filter is computed as the difference between two low-pass filters
!   with cutoff frequencies  1/PH and 1/PL, respectively (or FCH and FCL).
!
!   This band-pass filter has a transfer function which changes  from approximately
!   zero to one and one to zero in the transition bands about the ideal cutoff frequencies 1/PH and 1/PL), 
!   that is from (1/PH - 1/K) to (1/PH + 1/K) and (1/PL - 1/K) to (1/PL + 1/K), respectively. 
!   The  user  must  specify  the  two cutoff periods and the number of filter coefficients,
!   which must be odd.
!
!   The  user must also choose the number of filter terms, K, so that:
!
!   - 0<=(1/PH - 1/K)
!   - (1/PH +  1.3/(K+1)) <= (1/PL -  1.3/(K+1))
!   - (1/PL + 1/K) < 0.5
!
!   However, if the optional logical argument NOTEST_FC is used and is set to true, the two tests
!
!   - 0<=(1/PH - 1/K)
!   - (1/PL + 1/K) < 0.5
!
!   are bypassed.
!
!   In addition, K must be chosen as a compromise between:
!
!      1)  A sharp cutoff, that is, 1/K small; and
!
!      2)  Minimizing the number of data points lost by the filtering operations
!          (e.g. (K-1)/2 data points will be lost from each end of the series).
!
!   The subroutine returns the difference between the two corresponding normalized low-pass filter coefficients
!   as computed by function LP_COEF.   
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York, Chapter 6.
!
!   (2) Duchon, C., 1979:
!             Lanczos filtering in one and two dimensions.
!             Journal of applied meteorology, vol. 18, 1016-1022.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Logical_Constants, only : true
    use Reals_Constants,   only : zero, half, one, c1_3
    use Char_Constants,    only : tseries_error26, tseries_error27, tseries_error29, tseries_error30,  &
                                  tseries_error31, tseries_error32, tseries_error34, tseries_error36,  &
                                  tseries_error37
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: pl, ph, k
!
    real(stnd), intent(in), optional :: fch, fcl
!
    real(stnd), dimension(k) :: bdcoef
!
    logical(lgl), intent(in), optional :: notest_fc
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)               :: fcl2, fch2, tmp
    real(stnd), dimension(k) :: coefl, coefh
!
    logical(lgl) :: test_fc
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='bd_coef'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    if ( k<3_i4b )              &
    call merror( name_proc//tseries_error26 )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error27 )
!
    if ( present(fcl) ) then
!
        if ( fcl<=zero .or. fcl>=half )            &
        call merror( name_proc//tseries_error29 )
!
        fcl2  = fcl
!
    else
!
        if ( pl<=2_i4b )            &
        call merror( name_proc//tseries_error31 )
!
        fcl2  = one/real( pl, stnd )
!
    end if
!
    if ( present(fch) ) then
!
        if ( fch<=zero .or. fch>=half )            &
        call merror( name_proc//tseries_error30 )
!
        fch2  = fch
!
    else
!
        if ( ph<=2_i4b )            &
        call merror( name_proc//tseries_error32 )
!
        fch2  = one/real( ph, stnd )
!
    end if
!
    test_fc = true
!
    if ( present( notest_fc ) ) then
        test_fc = .not.notest_fc
    end if
!
!   BYPASS THE TESTS ON THE CUTOFF FREQUENCY IF REQUIRED.
!
    if ( test_fc ) then
!
        tmp = one/real( k, stnd )
!
        if ( (fch2 - tmp) < zero  )            &
        call merror( name_proc//tseries_error34 )
!
        if ( (fcl2 + tmp) >= half )            &
        call merror( name_proc//tseries_error36 )
!
    end if
!
!   CHECK THAT k IS LARGE ENOUGH TO ENSURE UNIT RESPONSE AT BAND CENTER.
!
    tmp  = c1_3/real( k + 1_i4b, stnd )
!
    if ( (fch2 + tmp)  >  (fcl2 - tmp) )    &
    call merror( name_proc//tseries_error37 )
!
!   CALCULATE THE LOW-PASS FILTER COEFFICIENTS CORRESPONDING TO THE TWO CUTOFF FREQUENCIES fcl2 and fch2.
!
    coefl(:k) = lp_coef( pl, k, fc=fcl2, notest_fc=notest_fc )
    coefh(:k) = lp_coef( ph, k, fc=fch2, notest_fc=notest_fc )
!
!   COMPUTE THE BAND-PASS COEFFICIENTS AS THE DIFFERENCE IN WEIGHT FUNCTIONS
!   FOR THE TWO NORMALIZED LOW-PASS FILTERS WITH CUTOFFS fcl2 AND fch2.
!
    bdcoef(:k) = coefl(:k) - coefh(:k)
!
!
! END OF FUNCTION bd_coef
! _______________________
!
    end function bd_coef
!
! =========================================================================================
!
    function bd_coef2( pl, ph, k, fch, fcl, win, notest_fc ) result( bdcoef2 )
!
! Purpose
! _______
!
!   Function BD_COEF2 computes the K-term least squares approximation to an 
!   -ideal- band pass filter with cutoff periods PL and PH (e.g. cutoff frequencies 1/PL and 1/PH)
!   by windowed filtering (e.g. Hamming window is used).
!   
!   PL and PH are expressed in number of points, i.e. PL=6(18) and PH=32(96) selects periods
!   between 1.5 yrs and 8 yrs for quarterly(monthly) data).
!
!   Alternatively, the user can directly specify the two cutoff frequencies, FCL and FCH, 
!   corresponding to PL and PH.
!   
!
! Arguments
! _________
!
!   PL          (INPUT) integer(i4b)
!               Minimum period of oscillation of desired component. The corresponding
!               cutoff frequency is 1/PL.
!               PL must be greater than two and must verify the following inequalities:
!
!               - 1/PH        < 1/PL
!               - 1/PL + 1/K  <  0.5
!
!   PH          (INPUT) integer(i4b)
!               Maximum period of oscillation of desired component. The corresponding
!               cutoff frequency is 1/PH.
!               PH must be greater than two and 1/PH must verify the following inequalities:
!
!               - 0    <= 1/PH - 1/K 
!               - 1/PH <  1/PL
!
!   K           (INPUT) integer(i4b)
!               The number of filter terms to be computed.
!               K must be greater or equal to 3 and odd.
!
!   FCH         (INPUT, OPTIONAL) real(stnd)
!               The user chosen (low) cutoff frequency in cycles per sample interval.
!               If the optional argument FCH is used, the PH argument is not used to
!               determine the (low) cutoff frequency.
!
!               FCH must verify the following inequalities:
!
!               - 0   <= FCH - 1/K 
!               - FCH <  FCL
!
!   FCL         (INPUT, OPTIONAL) real(stnd)
!               The user chosen (high) cutoff frequency in cycles per sample interval.
!               If the optional argument FCL is used, the PL argument is not used to
!               determine the cutoff (high) frequency.
!
!               FCL must verify the following inequalities:
!
!               - FCH        < FCL
!               - FCL + 1/K  < 0.5
!
!   WIN         (INPUT, OPTIONAL) real(stnd)
!               By default, Hamming window filtering is used (i.e. WIN=0.54).
!               Set WIN=0.5 for Hanning window or WIN=1 for rectangular window.
!
!               WIN must be greater or equal to O.5 and less or equal to 1, otherwise
!               WIN is reset to 0.54.
!
!   NOTEST_FC   (INPUT, OPTIONAL) logical(lgl)
!               On input, if this optional logical argument is set to true, the tests on the
!               cutoff frequencies FCH and FCL  (e.g. FCH - 1/K >= 0 and FCL + 1/K <  0.5) are bypassed.
!               However, in that case FCH and FCL must still verify the inequalities FCH > 0, FCL < 0.5
!               and  FCH < FCL .
!
!
! Further Details
! _______________
!
!   Function BD_COEF2 computes symmetric linear band-pass  filter coefficients
!   using a  least squares  approximation  to an  ideal band-pass  filter. The Hamming
!   window is used to reduce  overshoot and  ripple in the transfert function of the 
!   ideal low-pass  filter.
!
!   This band-pass filter is computed as the difference between two low-pass filters
!   with cutoff frequencies  1/PH and 1/PL, respectively (or FCH and FCL).
!
!   This band-pass filter has a transfer function which changes  from approximately
!   zero to one and one to zero in the transition bands about the ideal cutoff frequencies 1/PH and 1/PL), 
!   that is from (1/PH - 1/K) to (1/PH + 1/K) and (1/PL - 1/K) to (1/PL + 1/K), respectively. 
!   The  user  must  specify  the  two cutoff periods and the number of filter coefficients,
!   which must be odd. The  user must also choose the number of filter terms, K, so that:
!
!   - 0            <= (1/PH - 1/K)
!   - 1/PH         <  1/PL
!   - (1/PL + 1/K) <  0.5
!
!   However, if the optional logical argument NOTEST_FC is used and is set to true, the two tests
!
!   - 0            <= (1/PH - 1/K)
!   - (1/PL + 1/K) <  0.5
!
!   are bypassed.
!
!   The overshoot and the associated ripples in the ideal transfert function are reduced by the use
!   of the Hamming window.
!
!   In addition, K must be chosen as a compromise between:
!
!      1)  A sharp cutoff, that is, 1/K small; and
!
!      2)  Minimizing the number of data points lost by the filtering operations
!          (e.g. (K-1)/2 data points will be lost from each end of the series).
!
!   The subroutine returns the difference between the two corresponding normalized low-pass filter coefficients
!   as computed by function LP_COEF2.   
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Logical_Constants, only : true
    use Reals_Constants,   only : zero, half, one
    use Char_Constants,    only : tseries_error26, tseries_error27, tseries_error29, tseries_error30,  &
                                  tseries_error31, tseries_error32, tseries_error34, tseries_error36,  &
                                  tseries_error49
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: pl, ph, k
!
    real(stnd), intent(in), optional :: fch, fcl, win
!
    real(stnd), dimension(k) :: bdcoef2
!
    logical(lgl), intent(in), optional :: notest_fc
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)               :: fcl2, fch2, tmp
    real(stnd), dimension(k) :: coefl, coefh
!
    logical(lgl) :: test_fc
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='bd_coef2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    if ( k<3_i4b )              &
    call merror( name_proc//tseries_error26 )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error27 )
!
    if ( present(fcl) ) then
!
        if ( fcl<=zero .or. fcl>=half )            &
        call merror( name_proc//tseries_error29 )
!
        fcl2  = fcl
!
    else
!
        if ( pl<=2_i4b )            &
        call merror( name_proc//tseries_error31 )
!
        fcl2  = one/real( pl, stnd )
!
    end if
!
    if ( present(fch) ) then
!
        if ( fch<=zero .or. fch>=half )            &
        call merror( name_proc//tseries_error30 )
!
        fch2  = fch
!
    else
!
        if ( ph<=2_i4b )            &
        call merror( name_proc//tseries_error32 )
!
        fch2  = one/real( ph, stnd )
!
    end if
!
    test_fc = true
!
    if ( present( notest_fc ) ) then
        test_fc = .not.notest_fc
    end if
!
!   BYPASS THE TESTS ON THE CUTOFF FREQUENCY IF REQUIRED.
!
    if ( test_fc ) then
!
        tmp = one/real( k, stnd )
!
        if ( (fch2 - tmp) < zero  )            &
        call merror( name_proc//tseries_error34 )
!
        if ( (fcl2 + tmp) >= half )            &
        call merror( name_proc//tseries_error36 )
!
    end if
!
    if ( fch2  >=  fcl2 )   &
    call merror( name_proc//tseries_error49 )
!
!   CALCULATE THE LOW-PASS FILTER COEFFICIENTS CORRESPONDING TO THE TWO CUTOFF FREQUENCIES fcl2 and fch2.
!
    coefl(:k) = lp_coef2( pl, k, fc=fcl2, win=win, notest_fc=notest_fc )
    coefh(:k) = lp_coef2( ph, k, fc=fch2, win=win, notest_fc=notest_fc )
!
!   COMPUTE THE BAND-PASS COEFFICIENTS AS THE DIFFERENCE IN WEIGHT FUNCTIONS
!   FOR THE TWO NORMALIZED LOW-PASS FILTERS WITH CUTOFFS fcl2 AND fch2.
!
    bdcoef2(:k) = coefl(:k) - coefh(:k)
!
!
! END OF FUNCTION bd_coef2
! ________________________
!
    end function bd_coef2
!
! =========================================================================================
!
    function pk_coef( freq, k, notest_freq ) result( pkcoef )
!
! Purpose
! _______
!
!   Function PK_COEF computes the K-term least squares approximation to an 
!   -ideal- band pass filter with peak response near one at the single frequency
!   FREQ (e.g. the peak response is at period=1/FREQ).
!   
!
! Arguments
! _________
!
!   FREQ         (INPUT) real(stnd)
!                The band pass filter will have unit response at the single frequency FREQ.
!                FREQ is expressed in cycles per sample interval.
!
!                The frequency FREQ must also be greater or equal to ( 1.3/(K+1) + 1/K ) and less than
!                0.5 - ( 1.3/(K+1) + 1/K ).
!
!   K            (INPUT) integer(i4b)
!                The number of filter terms to be computed.
!                K must be greater or equal to 3 and odd.
!
!   NOTEST_FREQ  (INPUT, OPTIONAL) logical(lgl)
!                On input, if this optional logical argument is set to true, the frequency FREQ must
!                only be greater or equal to 1.3/(K+1) and less than 0.5 - 1.3/(K+1).
!
!
! Further Details
! _______________
!
!   Function PK_COEF computes symmetric linear band-pass  filter coefficients
!   using a  least squares  approximation  to an  ideal band-pass  filter  that has
!   convergence  factors  which reduce  overshoot and  ripple  (Bloomfield, 1976).
!
!   This band-pass filter is computed as the difference between two low-pass filters
!   with cutoff frequencies FCL and FCH, respectively (Duchon, 1979).
!
!   This band-pass filter has a transfer function which changes  from approximately
!   zero to one and one to zero in the transition bands about the cutoff frequencies FCH and FCL, 
!   that is from (FCH - 1/K) to FREQ and FREQ to (FCL + 1/K), respectively. 
!   The  user  must  specify  the frequency FREQ with unit response and the number of filter
!   coefficients, which must be odd. The  user must also choose the number of filter terms,
!   K, as a compromise between:
!
!      1)  A sharp cutoff, that is, 1/K small; and
!
!      2)  Minimizing the number of data points lost by the filtering operations
!          (e.g. (K-1)/2 data points will be lost from each end of the series).
!
!   The subroutine computes the two cutoff frequencies FCL and FCH as described by Duchon (1979) 
!   and returns the difference between the two corresponding normalized low-pass filter coefficients
!   as computed by function LP_COEF.
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York, Chapter 6.
!
!   (2) Duchon, C., 1979:
!             Lanczos filtering in one and two dimensions.
!             Journal of applied meteorology, vol. 18, 1016-1022.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Logical_Constants, only : true
    use Reals_Constants,   only : zero, half, one, c1_3
    use Char_Constants,    only : tseries_error26, tseries_error27, tseries_error38, tseries_error39,  &
                                  tseries_error40, tseries_error41, tseries_error42
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: k
!
    real(stnd), intent(in)   :: freq
!
    real(stnd), dimension(k) :: pkcoef
!
    logical(lgl), intent(in), optional :: notest_freq
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)               :: fcl, fch, tmp
    real(stnd), dimension(k) :: coefl, coefh
!
    integer(i4b) :: i
!
    logical(lgl) :: test_freq
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='pk_coef'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    if ( k<3_i4b )              &
    call merror( name_proc//tseries_error26 )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error27 )
!
    if ( freq<=zero .or. freq>=half )            &
    call merror( name_proc//tseries_error38 )
!
    test_freq = true
!
    if ( present( notest_freq ) ) then
        test_freq = .not.notest_freq
    end if
!
!   COMPUTE THE TWO CUTOFF FREQUENCIES fcl AND fch .
!
    tmp  = c1_3/real( k + 1_i4b, stnd )
    fch  = freq - tmp
    fcl  = freq + tmp
!
    if ( test_freq ) then
!
        tmp = one/real( k, stnd )
!
        if ( (fch - tmp) < zero  )            &
        call merror( name_proc//tseries_error39 )
!
        if ( (fcl + tmp) >= half )            &
        call merror( name_proc//tseries_error40 )
!
    else
!
        if ( fch < zero  )            &
        call merror( name_proc//tseries_error41 )
!
        if ( fcl >= half )            &
        call merror( name_proc//tseries_error42 )
!
    end if
!
!   CALCULATE THE LOW-PASS FILTER COEFFICIENTS CORRESPONDING TO THE TWO CUTOFF FREQUENCIES fcl and fch.
!
    i = 0_i4b
!
    coefl(:k) = lp_coef( i, k , fc=fcl, notest_fc=notest_freq )
    coefh(:k) = lp_coef( i, k , fc=fch, notest_fc=notest_freq )
!
!   COMPUTE THE BAND-PASS COEFFICIENTS AS THE DIFFERENCE IN WEIGHT FUNCTIONS
!   FOR THE TWO NORMALIZED LOW-PASS FILTERS WITH CUTOFFS fcl AND fch.
!
    pkcoef(:k) = coefl(:k) - coefh(:k)
!
!
! END OF FUNCTION pk_coef
! _______________________
!
    end function pk_coef
!
! =========================================================================================
!
    function moddan_coef( k, smooth_param ) result( moddancoef )
!
! Purpose
! _______
!
!   This function computes the impulse response function (e.g. weights) corresponding to a number
!   of applications of modified Daniell filters as done in subroutine MODDAN_FILTER.
! 
! 
! Arguments
! _________
!
!   K             (INPUT) integer(i4b)
!                 The number of filter weights to be computed.
!                 K must be equal to 2 * (2+sum(SMOOTH_PARAM(:)))- 1
!
!   SMOOTH_PARAM  (INPUT) integer(i4b), dimension(:)
!                 The array of the half-lengths of the modified Daniell filters to be applied.
!                 All the values in SMOOTH_PARAM(:) must be greater than 0.
!
!                 Size(SMOOTH_PARAM) must be greater or equal to 1.
!
!
! Further Details
! _______________
!
!   For definition, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror, assert
    use Reals_Constants,   only : zero, one
    use Char_Constants,    only : tseries_error50, tseries_error51
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in)                :: k
    integer(i4b), intent(in), dimension(:)  :: smooth_param
!
    real(stnd), dimension(k) :: moddancoef
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd), dimension(k) :: coef
!
    integer(i4b) :: k2, index, nparam
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='moddan_coef'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    nparam = size( smooth_param )
!
    if ( nparam<=0_i4b  )    &
    call merror( name_proc//tseries_error50 )
!
!   CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIEL FILTERS.
!
    if ( any( smooth_param(:nparam)<=0_i4b ) )   &
    call merror( name_proc//tseries_error51 )
!
    index = 2_i4b + sum( smooth_param(:nparam) )
    k2    = 2_i4b*index - 1_i4b
!
    call assert( logical( k==k2, lgl ), name_proc )
!
!   COMPUTE THE FILTER WEIGHTS.
!
    coef(:k)    = zero
    coef(index) = one
!
    call moddan_filter_rv( coef(:k), smooth_param(:), sym=one )
!
    moddancoef(:k) = coef(:k)
!
!
! END OF FUNCTION moddan_coef
! ___________________________
!
    end function moddan_coef
!
! =========================================================================================
!
    subroutine freq_func( nfreq, coef, freqr, four_freq, freq )
!
! Purpose
! _______
!
!   Subroutine FREQ_FUNC computes the frequency response function (e.g. the transfer function)
!   of the symmetric linear filter given by the argument COEF(:).
!   
!   The frequency response function is computed at NFREQ frequencies regularly sampled
!   between 0 and the Nyquist frequency if the optional logical argument FOUR_FREQ is not used
!   or at the NFREQ Fourier frequencies 2 * pi * j/nfreq for j=0 to NFREQ-1 if this argument is used
!   and set to true.
! 
! Arguments
! _________
!
!   NFREQ       (INPUT) integer(i4b)
!               The number of frequencies at which the frequency response function
!               must be evaluated.
!
!   COEF        (INPUT) real(stnd), dimension(:)
!               The array of symmetric linear filter coefficients.
!
!               Size(COEF) must be greater or equal to 3 and odd.
!
!   FREQR       (OUTPUT) real(stnd), dimension(NFREQ)
!               On output, the frequency response function.
!
!   FOUR_FREQ   (INPUT, OPTIONAL) logical(lgl)
!               On input, if this argument is set to true the frequency response function
!               is evaluated at the Fourier frequencies 2 * pi * j/nfreq for j=0 to NFREQ-1.
!
!   FREQ        (OUTPUT, OPTIONAL) real(stnd), dimension(NFREQ)
!               The NFREQ frequencies, in cycles per sample interval, at which
!               the frequency response function are evaluated.
!
!
! Further Details
! _______________
!
!   For more details, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York, Chapter 6.
!
!   (2) Oppenheim, A.V., and Schafer, R.W., 1999:
!            Discrete-Time Signal Processing.
!            Second Edition. Prentice-Hall, New Jersey.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror, assert, arth
    use Logical_Constants, only : true, false
    use Reals_Constants,   only : zero, one, two, pi
    use Char_Constants,    only : tseries_error45, tseries_error46, tseries_error47, tseries_error48
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: nfreq
!
    real(stnd), intent(in),  dimension(:)           :: coef
    real(stnd), intent(out), dimension(:)           :: freqr
    real(stnd), intent(out), dimension(:), optional :: freq
!
    logical(lgl), intent(in),              optional :: four_freq
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)                   :: tmp
    real(stnd), dimension(nfreq) :: tmpvec, freq2
!
    integer(i4b) :: k, khalf, kmid, i, i2, n
!
    logical(lgl) :: four_freq2, store_freq
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='freq_func'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    call assert( logical( nfreq==int(size(freqr),i4b), lgl ),     &
                 name_proc )
!
    store_freq = false
!
    if ( present(freq) ) then
        call assert( logical( nfreq==int(size(freq),i4b), lgl ),  &
                     name_proc )
        store_freq = true
    end if
!
    four_freq2 = false
!
    if ( present(four_freq) ) then
        four_freq2 = four_freq
    end if
!
    k = size( coef )
!
    if ( k<3_i4b )              &
    call merror( name_proc//tseries_error48 )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error45 )
!
!   CHECK THE INPUT VALUES OF THE SYMMETRIC LINEAR FILTER.
!
    if ( all( coef(:k)==zero ) )   &
    call merror( name_proc//tseries_error46 )
!
    khalf = ( k - 1_i4b )/2_i4b
!
    do i = 1_i4b, khalf
!
        i2 = k + 1_i4b - i
        tmp = coef(i) - coef(i2)
        if ( abs(tmp)>abs(coef(i2))*epsilon(zero) ) exit
!
    end do
!
!    do i = 1_i4b, khalf
!        i2 = k + 1_i4b - i
!        if ( coef(i)/=coef(i2) ) exit
!    end do
!
    if ( i<=khalf )       &
    call merror( name_proc//tseries_error47 )
!
!   COMPUTE THE TRANSFER FUNCTION.
!
    kmid  = ( k + 1_i4b )/2_i4b
!
    freqr(:nfreq)  = coef(kmid)
    tmpvec(:nfreq) = zero
!
    if ( four_freq2 ) then
!
        tmp = (two*pi)/real( nfreq, stnd )
!
        i   = nfreq/2
        i2  = i + 1_i4b
        n   = merge( i, i2, nfreq==i*2_i4b )
!
    else
!
        tmp = pi/real( nfreq-1_i4b, stnd )
!
        i2 = nfreq
!
    end if
!
    if ( store_freq ) then
        i = nfreq
    else
        i = i2
    end if
!
    freq2(:i)  = arth( zero, tmp, i )
!
    do i = 1_i4b, khalf
!
        tmpvec(:i2) = tmpvec(:i2) + freq2(:i2)
        freqr(:i2)  = freqr(:i2)  + two*coef(kmid+i)*cos( tmpvec(:i2) )
!
    end do
!
    if ( four_freq2 ) then
        tmpvec(2_i4b:n)              = freqr(2_i4b:n)
        freqr(nfreq:i2+1_i4b:-1_i4b) = tmpvec(2_i4b:n)
    end if
!
!   COMPUTE THE ASSOCIATED FREQUENCIES IF REQUIRED.
!
    if ( store_freq ) then
        freq(:nfreq)   = freq2(:nfreq)*( one/(two*pi) )
    end if
!
!
! END OF SUBROUTINE freq_func
! ___________________________
!
    end subroutine freq_func
!
! =========================================================================================
!
    subroutine symlin_filter_rv( vec, coef, trend, nfilt )
!
! Purpose
! _______
!
!   Subroutine SYMLIN_FILTER performs a symmetric filtering operation on an input
!   time series (e.g. the argument VEC).
!
!
! Arguments
! _________
!
!   VEC          (INPUT/OUTPUT) real(stnd), dimension(:)
!                On input, the vector containing the time series to be filtered.
!                On output, the filtered time series is returned in VEC(:NFILT).
!                Note that (size(COEF)-1)/2 data points will be lost from each end of the series,
!                so that NFILT (NFILT= size(VEC) - size(COEF) + 1) time observations are returned
!                and the remainig and ending part of VEC(:) is set to zero.
!
!                Size(VEC) must be greater or equal to 4.
!
!   COEF         (INPUT) real(stnd), dimension(:)
!                The array of symmetric linear filter coefficients.
!
!                Size(COEF) must be odd, greater or equal to 3 and less or equal to size(VEC).
!
!
!   TREND        (INPUT, OPTIONAL) integer(i4b)
!                If:
!
!                - TREND=+/-1 The mean of the time series is removed before time filtering
!                - TREND=+/-2 The drift from the time series is removed before time filtering
!                  by using the formula: drift = (VEC(size(VEC)) - VEC(1))/(size(VEC) - 1)
!                - TREND=+/-3 The least-squares line from the time series is removed before
!                  time filtering.
!
!                IF TREND=-1,-2 or -3, the mean, drift or least-squares line is reintroduced
!                post-filtering, respectively.
!                For other values of TREND nothing is done before or after filtering.
!
!   NFILT        (OUTPUT, OPTIONAL) integer(i4b)
!                The number of time observations in the filtered time series.
!                On output, NFILT= size(VEC) - size(COEF) + 1.
!
!
! Further Details
! _______________
!
!   The filtering is done in place and (size(COEF)-1)/2 observations will be lost from
!   each end of the time series.
!
!   Note, also, that the filtered time series is shifted in time and is stored in VEC(1:NFILT)
!   on output, with NFILT= size(VEC) - size(COEF) + 1.
!
!   The symmetric linear filter coefficients (e.g. the array COEF) can be computed with
!   the help of functions LP_COEF, LP_COEF2, HP_COEF, HP_COEF2, BD_COEF and BD_COEF2.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror, arth
    use Reals_Constants,   only : zero, one
    use Char_Constants,    only : tseries_error10, tseries_error43, tseries_error45, tseries_error46,    &
                                  tseries_error47
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:)  :: vec
    real(stnd),   intent(in),    dimension(:)  :: coef
!
    integer(i4b), intent(in), optional         :: trend
    integer(i4b), intent(out), optional        :: nfilt
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, k, khalf, i, i2, trend2, nvecf
!
    real(stnd)   :: orig, slope, tmp, offset
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='symlin_filter'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m = size( vec )
!
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    k = size( coef )
!
    if ( k<3_i4b .or. k>m )   &
    call merror( name_proc//tseries_error43 )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error45 )
!
!   CHECK THE INPUT VALUES OF THE SYMMETRIC LINEAR FILTER.
!
    if ( all( coef(:k)==zero ) )   &
    call merror( name_proc//tseries_error46 )
!
    khalf = ( k - 1_i4b )/2_i4b
!
    do i = 1_i4b, khalf
!
        i2 = k + 1_i4b - i
        tmp = coef(i) - coef(i2)
!
        if ( abs(tmp)>abs(coef(i2))*epsilon(zero) ) exit
!
    end do
!
!    do i = 1_i4b, khalf
!        i2 = k + 1_i4b - i
!        if ( coef(i)/=coef(i2) ) exit
!    end do
!
    if ( i<=khalf )       &
    call merror( name_proc//tseries_error47 )
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
!
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    i = abs( trend2 )
!
    if ( i>=1_i4b .and. i<=3_i4b ) then
        call detrend_rv( vec(:m), i, orig, slope )
    end if
!
!   NOW, FILTER THE TIME SERIES IN PLACE.
!
    nvecf = m - k + 1_i4b
    i2    = k - 1_i4b
!
    do i = 1_i4b, nvecf
!
        i2    = i2 + 1_i4b
!
        tmp = dot_product( coef(:k), vec(i:i2) )
        vec(i) = tmp
!
    end do
!
!   SET THE REMAINING OBSERVATIONS TO ZERO.
!
    vec(nvecf+1_i4b:m) =zero
!
!   ADD MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEAN.
!
            vec(:nvecf) = vec(:nvecf) + orig
!
        case( -2_i4b )
!
!           ADD DRIFT.
!         
            offset = real( (k - 1_i4b)/2_i4b, stnd )
!
            vec(:nvecf) = vec(:nvecf) + slope*arth( offset, one, nvecf )           
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINE.
!         
            offset = real( (k - 1_i4b)/2_i4b, stnd )
!         
            vec(:nvecf) = vec(:nvecf) + ( orig + slope*arth( offset, one, nvecf ) )       
!         
    end select
!
!   SET THE OPTIONAL ARGUMENT nfilt, IF REQUIRED.
!
    if ( present(nfilt) ) then
        nfilt = nvecf
    end if
!
!
! END OF SUBROUTINE symlin_filter_rv
! __________________________________
!
    end subroutine symlin_filter_rv
!
! =========================================================================================
!
    subroutine symlin_filter_rm( mat, coef, trend, nfilt )
!
! Purpose
! _______
!
!   Subroutine SYMLIN_FILTER performs a symmetric filtering operation on an input multi-channel
!   time series (e.g. the argument MAT).
!
!
! Arguments
! _________
!
!   MAT          (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                The multi-channel time series matrix to be filtered.
!                Each column of MAT corresponds to one observation.
!                On output, the multi-channel filtered time series are returned in MAT(:,:NFILT).
!                Note that (size(COEF)-1)/2 observations will be lost from each end of the multi-channel
!                series, so that NFILT (NFILT= size(MAT,2) - size(COEF) + 1) time observations are returned
!                and the remainig part of MAT(:,:) is set to zero.
!
!                Size(MAT,2) must be greater or equal to 4.
!
!   COEF         (INPUT) real(stnd), dimension(:)
!                The array of symmetric linear filter coefficients.
!
!                Size(COEF) must be odd, greater or equal to 3 and less or equal to size(MAT,2).
!
!   TREND        (INPUT, OPTIONAL) integer(i4b)
!                If:
!
!                - TREND=+/-1 The means of the time series are removed before time filtering
!                - TREND=+/-2 The drifts from the time series are removed before time filtering
!                  by using the formula: drift(:) = (MAT(:,size(MAT,2)) - MAT(:,1))/(size(MAT,2) - 1)
!                - TREND=+/-3 The least-squares lines from the time series are removed before
!                  time filtering.
!
!                IF TREND=-1,-2 or -3, the means, drifts or least-squares lines are reintroduced
!                post-filtering, respectively.
!                For other values of TREND nothing is done before or after filtering.
!
!   NFILT        (OUTPUT, OPTIONAL) integer(i4b)
!                The number of time observations in the filtered multi-channel time series.
!                On output, NFILT= size(MAT,2) - size(COEF) + 1.
!
!
! Further Details
! _______________
!
!   The filtering is done in place and (size(COEF)-1)/2 observations will be lost from
!   each end of the multi-channel series.
!
!   Note, also, that the filtered multi-channel time series is shifted in time and is stored
!   in MAT(:,1:NFILT) on output, with NFILT= size(MAT,2) - size(COEF) + 1.
!
!   The symmetric linear filter coefficients (e.g. the array COEF) can be computed with
!   the help of functions LP_COEF, LP_COEF2, HP_COEF, HP_COEF2, BD_COEF and BD_COEF2.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Reals_Constants,   only : zero
    use Char_Constants,    only : tseries_error10, tseries_error44, tseries_error45, tseries_error46,   &
                                  tseries_error47
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:,:)  :: mat
    real(stnd),   intent(in),    dimension(:)    :: coef
!
    integer(i4b), intent(in), optional           :: trend
    integer(i4b), intent(out), optional        :: nfilt
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, m, k, khalf, i, i2, trend2, nmatf
!
    real(stnd), dimension(size(mat,1)) :: orig, slope, tmpvec
    real(stnd)                         :: tmp, offset
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='symlin_filter'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( mat, 1 )
!
    if ( n<=0_i4b )  return
!
    m = size( mat, 2 )
!
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    k = size( coef )
!
    if ( k<3_i4b .or. k>m )   &
    call merror( name_proc//tseries_error44 )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error45 )
!
!   CHECK THE INPUT VALUES OF THE SYMMETRIC LINEAR FILTER.
!
    if ( all( coef(:k)==zero ) )   &
    call merror( name_proc//tseries_error46 )
!
    khalf = ( k - 1_i4b )/2_i4b
!
    do i = 1_i4b, khalf
!
        i2 = k + 1_i4b - i
        tmp = coef(i) - coef(i2)
        if ( abs(tmp)>abs(coef(i2))*epsilon(zero) ) exit
!
    end do
!
!    do i = 1_i4b, khalf
!        i2 = k + 1_i4b - i
!        if ( coef(i)/=coef(i2) ) exit
!    end do
!
    if ( i<=khalf )       &
    call merror( name_proc//tseries_error47 )
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
!
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    i = abs( trend2 )
!
    if ( i>=1_i4b .and. i<=3_i4b ) then
        call detrend_rm( mat(:n,:m), i, orig(:n), slope(:n) )
    end if
!
!   NOW, FILTER THE TIME SERIES IN PLACE.
!
    nmatf = m - k + 1_i4b
    i2    = k - 1_i4b
!
    do i = 1_i4b, nmatf
!
        i2    = i2 + 1_i4b
!
        tmpvec(:n) = matmul( mat(:n,i:i2), coef(:k) )
        mat(:n,i)  = tmpvec(:n)
!
    end do
!
!   SET THE REMAINING OBSERVATIONS TO ZERO.
!
    mat(:n,nmatf+1_i4b:m) = zero
!
!   ADD MEANS, DRIFTS OR LINEAR LEAST SQUARES LINES FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEANS.
!
            do i = 1_i4b, nmatf
                mat(:n,i) = mat(:n,i) + orig(:n)
            end do
!
        case( -2_i4b )
!
!           ADD DRIFTS.
!
            offset = real( (k - 3_i4b)/2_i4b, stnd )
!
            do i = 1_i4b, nmatf
                tmp       = real( i, stnd ) + offset
                mat(:n,i) = mat(:n,i) + tmp*slope(:n)           
            end do
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINES.
!         
            offset = real( (k - 3_i4b)/2_i4b, stnd )
!
            do i = 1_i4b, nmatf
                tmp       = real( i, stnd ) + offset
                mat(:n,i) = mat(:n,i) + tmp*slope(:n) + orig(:n)      
            end do
!
    end select
!
!   SET THE OPTIONAL ARGUMENT nfilt, IF REQUIRED.
!
    if ( present(nfilt) ) then
        nfilt = nmatf
    end if
!
!
! END OF SUBROUTINE symlin_filter_rm
! __________________________________
!
    end subroutine symlin_filter_rm
!
! =========================================================================================
!
    subroutine symlin_filter2_rv( vec, coef, trend, usefft, initfft )
!
! Purpose
! _______
!
!   Subroutine SYMLIN_FILTER2 performs a symmetric filtering operation on
!   an input time series (e.g. the argument VEC).
!
!
! Arguments
! _________
!
!   VEC          (INPUT/OUTPUT) real(stnd), dimension(:)
!                On input, the vector containing the time series to be filtered.
!                On output, the filtered time series is returned.
!
!                Size(VEC) must be greater or equal to 4.
!
!   COEF         (INPUT) real(stnd), dimension(:)
!                The array of symmetric linear filter coefficients.
!
!                Size(COEF) must be odd, greater or equal to 3 and less or equal to size(VEC).
!
!   TREND        (INPUT, OPTIONAL) integer(i4b)
!                If:
!
!                - TREND=+/-1 The mean of the time series is removed before time filtering
!                - TREND=+/-2 The drift from the time series is removed before time filtering
!                  by using the formula: drift = (VEC(size(VEC)) - VEC(1))/(size(VEC) - 1)
!                - TREND=+/-3 The least-squares line from the time series is removed before
!                  time filtering.
!
!                IF TREND=-1,-2 or -3, the mean, drift or least-squares line is reintroduced
!                post-filtering, respectively.
!                For other values of TREND nothing is done before or after filtering.
!
!   USEFFT       (INPUT, OPTIONAL) logical(lgl)
!                On input, if USEFFT is used and is set to true, the symmetric linear filter is applied
!                to the argument VEC by using a Fast Fourier Transform and the convolution theorem.
!
!   INITFFT      (INPUT, OPTIONAL) logical(lgl)
!                On entry, if INITFFT is set to false, it is assumed that a call to subroutine
!                INIT_FFT has been done before calling subroutine SYMLIN_FILTER2 in order to 
!                sets up constants and functions for use by subroutine FFT_ROW which is called inside
!                subroutine SYMLIN_FILTER2 (the call to INIT_FFT must have the following form: 
!
!                       call init_fft( shape(VEC) )
!
!                If INITFFT is set to true, the call to INIT_FFT is done inside subroutine
!                SYMLIN_FILTER2 and a call to END_FFT is also done before leaving
!                subroutine SYMLIN_FILTER2.
!                This optional argument has an effect only if argument USEFFT is used with the value true.
!
!                The default is INITFFT=true .
!
!
! Further Details
! _______________
!
!   No time observations will be lost, however the first and last (size(COEF)-1)/2 time
!   observations are affected by end effects.
!
!   If USEFFT is used with the value true, the values at both ends of the output series
!   are computed by assuming that the input series is part of a periodic sequence of
!   period size(VEC). Otherwise, each end of the filtered time series is estimated by truncated
!   the symmetric linear filter coefficients array.
!
!   The symmetric linear filter coefficients (e.g. the array COEF) can be computed with
!   the help of functions LP_COEF, LP_COEF2, HP_COEF, HP_COEF2, BD_COEF and BD_COEF2.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror, arth
    use Logical_Constants, only : true, false
    use Reals_Constants,   only : zero, one
    use Char_Constants,    only : tseries_error10, tseries_error43, tseries_error45, tseries_error46,    &
                                  tseries_error47, allocate_error
    use FFT_Procedures,    only : init_fft, fft_row, end_fft
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:)  :: vec
    real(stnd),   intent(in),    dimension(:)  :: coef
!
    integer(i4b), intent(in), optional         :: trend
!
    logical(lgl), intent(in), optional         :: usefft, initfft
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd), dimension(size(vec)) :: vec2
    real(stnd)                       :: orig, slope, tmp
!
    complex(stnd), dimension(:), allocatable :: vecc
!
    integer(i4b) :: m, k, i, i1, i2, trend2, kmid, khalf
!   
    integer :: iok
!
    logical(lgl) :: usefft2, initfft2
!
#ifdef _OPENMP
    logical      :: test_par
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='symlin_filter2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m = size( vec )
!
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    k = size( coef )
!
    if ( k<3_i4b .or. k>m )   &
    call merror( name_proc//tseries_error43 )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error45 )
!
!   CHECK THE INPUT VALUES OF THE SYMMETRIC LINEAR FILTER.
!
    if ( all( coef(:k)==zero ) )   &
    call merror( name_proc//tseries_error46 )
!
    khalf = ( k - 1_i4b )/2_i4b
!
    do i = 1_i4b, khalf
!
        i2 = k + 1_i4b - i
        tmp = coef(i) - coef(i2)
        if ( abs(tmp)>abs(coef(i2))*epsilon(zero) ) exit
!
    end do
!
!    do i = 1_i4b, khalf
!        i2 = k + 1_i4b - i
!        if ( coef(i)/=coef(i2) ) exit
!    end do
!
    if ( i<=khalf )  &
    call merror( name_proc//tseries_error47 )
!
    usefft2 = false
!
    if ( present(usefft) ) then
        usefft2 = usefft
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
!
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    i = abs( trend2 )
!
    if ( i>=1_i4b .and. i<=3_i4b ) then
        call detrend_rv( vec(:m), i, orig, slope )
    end if
!
    if ( usefft2 ) then
!
!       ALLOCATE WORK ARRAY.
!
        allocate( vecc(m), stat=iok )
!
        if ( iok/=0 ) then
            call merror( name_proc//allocate_error )
        end if
!
!       COMPUTE THE TRANSFERT FUNCTION OF THE SYMMETRIC LINEAR FILTER AT THE FOURIER FREQUENCIES.
!
        call freq_func( m, coef(:k), vec2(:m), four_freq=true )
!
!       NOW, APPLY THE FILTER USING THE FFT AND THE CONVOLUTION THEOREM.
!
!       INITIALIZE THE FFT SUBROUTINE.
!
        if ( initfft2 ) then
            call init_fft( m )
        end if
!
!       TRANSFORM THE TIME SERIES.
!
        vecc(:m) = cmplx( vec(:m), zero, kind=stnd )
!
        call fft_row( vecc(:m), forward=true  )
!
!       MULTIPLY THE FOURIER TRANSFORM OF THE TIME SERIES
!       BY THE TRANSFERT FUNCTION OF THE FILTER.
!
        vecc(:m) = vecc(:m)*vec2(:m)
!
!       INVERT THE SEQUENCE BACK TO GET THE FILTERED TIME SERIES.
!
        call fft_row( vecc(:m), forward=false )
!
        vec(:m) = real( vecc(:m),  kind=stnd )
!
!       DEALLOCATE WORK ARRAYS.
!
        if ( initfft2 ) then
            call end_fft()
        end if
!
        deallocate( vecc )
!
    else
!
!       MAKE A COPY OF THE TIME SERIES.
!
        vec2(:m) = vec(:m)
!
!       NOW, FILTER THE TIME SERIES.
!
        kmid  = ( k + 1_i4b )/2_i4b
!   
#ifdef _OPENMP
        i1 = omp_get_num_procs()
        i2 = omp_get_max_threads()
        test_par = .not.( omp_in_parallel() )   .and.      &
                   m>=omp_limit                 .and.      &
                   i1>1_i4b                     .and.      &
                   i2>1_i4b
!
        if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(i,i1,i2)                      &
!$OMP         ,SHARED(m,k,kmid,khalf,vec,vec2)
!
!
!$OMP DO SCHEDULE(STATIC) 
!
            do i = 1_i4b, khalf
!
                i1 = kmid + 1_i4b - i
                i2 = khalf + i
!
                vec(i) = dot_product( coef(i1:k), vec2(:i2) )
!
            end do
!
!$OMP END DO NOWAIT
!
!
!$OMP DO SCHEDULE(STATIC) 
!
            do i = kmid, m-khalf
!
                i1 = i - kmid + 1_i4b
                i2 = k  + i1  - 1_i4b
!
                vec(i) = dot_product( coef(:k), vec2(i1:i2) )
!
            end do
!
!$OMP END DO NOWAIT
!
!
!$OMP DO SCHEDULE(STATIC) 
!
            do i =  m-khalf+1_i4b, m
!
                i1 = k - khalf - i + m
                i2 = m - i1 + 1_i4b
!
                vec(i) = dot_product( coef(:i1), vec2(i2:m) )
!
            end do
!
!$OMP END DO
!
!$OMP END PARALLEL
!
        else
!
#endif
            i1 = kmid + 1_i4b
            i2 = khalf
!
            do i = 1_i4b, khalf
!
                i1 = i1 - 1_i4b
                i2 = i2 + 1_i4b
!
                vec(i) = dot_product( coef(i1:k), vec2(:i2) )
!
            end do
!
            i1 = 0_i4b
            i2 = k - 1_i4b
!
            do i = kmid, m-khalf
!
                i1 = i1 + 1_i4b
                i2 = i2 + 1_i4b
!
                vec(i) = dot_product( coef(:k), vec2(i1:i2) )
!
            end do
!
            i1 = k
            i2 = m - k + 1_i4b
!
            do i =  m-khalf+1_i4b, m
!
                i1 = i1 - 1_i4b
                i2 = i2 + 1_i4b
!
                vec(i) = dot_product( coef(:i1), vec2(i2:m) )
!
            end do
!
#ifdef _OPENMP
        end if
#endif
!
    end if
!
!   ADD MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEAN.
!
            vec(:m) = vec(:m) + orig
!
        case( -2_i4b )
!
!           ADD DRIFT.
!         
            vec(:m) = vec(:m) + slope*arth( zero, one, m )           
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINE.
!         
            vec(:m) = vec(:m) + ( orig + slope*arth( zero, one, m ) )       
!         
    end select
!
!
! END OF SUBROUTINE symlin_filter_rv
! __________________________________
!
    end subroutine symlin_filter2_rv
!
! =========================================================================================
!
    subroutine symlin_filter2_rm( mat, coef, trend, usefft, initfft )
!
! Purpose
! _______
!
!   Subroutine SYMLIN_FILTER2 performs a symmetric filtering operation on an input multi-channel
!   time series (e.g. the argument MAT).
!
!
! Arguments
! _________
!
!   MAT          (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                The multi-channel time series matrix to be filtered.
!                Each column of MAT corresponds to one observation.
!                On output, the multi-channel filtered time series are returned.
!
!                Size(MAT,2) must be greater or equal to 4.
!
!   COEF         (INPUT) real(stnd), dimension(:)
!                The array of symmetric linear filter coefficients.
!
!                Size(COEF) must be odd, greater or equal to 3 and less or equal to size(MAT,2).
!
!   TREND        (INPUT, OPTIONAL) integer(i4b)
!                If:
!
!                - TREND=+/-1 The means of the time series are removed before time filtering
!                - TREND=+/-2 The drifts from the time series are removed before time filtering
!                  by using the formula: drift(:) = (MAT(:,size(MAT,2)) - MAT(:,1))/(size(MAT,2) - 1)
!                - TREND=+/-3 The least-squares lines from the time series are removed before
!                  time filtering.
!
!                IF TREND=-1,-2 or -3, the means, drifts or least-squares lines are reintroduced
!                post-filtering, respectively.
!                For other values of TREND nothing is done before or after filtering.
!
!   USEFFT       (INPUT, OPTIONAL) logical(lgl)
!                On input, if USEFFT is used and is set to true, the symmetric linear filter is applied
!                to the argument VEC by using a Fast Fourier Transform and the convolution theorem.
!
!   INITFFT      (INPUT, OPTIONAL) logical(lgl)
!                On entry, if INITFFT is set to false, it is assumed that a call to subroutine
!                INIT_FFT has been done before calling subroutine SYMLIN_FILTER2 in order to 
!                sets up constants and functions for use by subroutine FFT_ROW which is called inside
!                subroutine SYMLIN_FILTER2 (the call to INIT_FFT must have the following form: 
!
!                      call init_fft( shape(MAT), dim=2_i4b )
!
!                If INITFFT is set to true, the call to INIT_FFT is done inside subroutine
!                SYMLIN_FILTER2 and a call to END_FFT is also done before leaving
!                subroutine SYMLIN_FILTER2.
!                This optional argument has an effect only if argument USEFFT is used with the value true.
!
!                The default is INITFFT=true .
!
!
! Further Details
! _______________
!
!   No time observations will be lost, however the first and last (size(COEF)-1)/2 time
!   observations are affected by end effects.
!
!   If USEFFT is used with the value true, the values at both ends of the output multi-channel
!   time series are computed by assuming that the input multi-channel series is part of a periodic
!   sequence of period size(VEC). Otherwise, each end of the filtered multi-channel time series is
!   estimated by truncated the symmetric linear filter coefficients array.
!
!   The symmetric linear filter coefficients (e.g. the array COEF) can be computed with
!   the help of functions LP_COEF, LP_COEF2, HP_COEF, HP_COEF2, BD_COEF and BD_COEF2.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Logical_Constants, only : true, false
    use Reals_Constants,   only : zero, one
    use Char_Constants,    only : tseries_error10, tseries_error44, tseries_error45, tseries_error46,   &
                                  tseries_error47, allocate_error
    use FFT_Procedures,    only : init_fft, fft_row, end_fft
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:,:)  :: mat
    real(stnd),   intent(in),    dimension(:)    :: coef
!
    integer(i4b), intent(in), optional           :: trend
!
    logical(lgl), intent(in), optional           :: usefft, initfft
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd), dimension(:,:),        allocatable :: mat2
    real(stnd), dimension(:),          allocatable :: freqr
    real(stnd), dimension(size(mat,1))             :: orig, slope
    real(stnd)                                     :: tmp
!
    complex(stnd), dimension(:,:),     allocatable :: matc
!
    integer(i4b) :: n, m, k, i, i1, i2, trend2, kmid, khalf
!   
    integer :: iok
!
    logical(lgl) :: usefft2, initfft2
!
#ifdef _OPENMP
!
    integer(i4b)                       :: j
!
    real(stnd), dimension(size(coef))  :: coef2
    real(stnd), dimension(size(mat,1)) :: vec
!
    logical                            :: test_par
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='symlin_filter2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( mat, 1 )
!
    if ( n<=0_i4b )  return
!
    m = size( mat, 2 )
!
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    k = size( coef )
!
    if ( k<3_i4b .or. k>m )   &
    call merror( name_proc//tseries_error44 )
!
    if ( (k/2_i4b)*2_i4b==k )   &
    call merror( name_proc//tseries_error45 )
!
!   CHECK THE INPUT VALUES OF THE SYMMETRIC LINEAR FILTER.
!
    if ( all( coef(:k)==zero ) )   &
    call merror( name_proc//tseries_error46 )
!
    khalf = ( k - 1_i4b )/2_i4b
!
    do i = 1_i4b, khalf
!
        i2 = k + 1_i4b - i
        tmp = coef(i) - coef(i2)
!
        if ( abs(tmp)>abs(coef(i2))*epsilon(zero) ) exit
!
    end do
!
!    do i = 1_i4b, khalf
!        i2 = k + 1_i4b - i
!        if ( coef(i)/=coef(i2) ) exit
!    end do
!
    if ( i<=khalf )       &
    call merror( name_proc//tseries_error47 )
!
    usefft2 = false
!
    if ( present(usefft) ) then
        usefft2 = usefft
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
!
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    i = abs( trend2 )
!
    if ( i>=1_i4b .and. i<=3_i4b ) then
!
        call detrend_rm( mat(:n,:m), i, orig(:n), slope(:n) )
!
    end if
!
    if ( usefft2 ) then
!
!       ALLOCATE WORK ARRAY.
!
        allocate( matc(n,m), freqr(m), stat=iok )
!
        if ( iok/=0 ) then
            call merror( name_proc//allocate_error )
        end if
!
!       COMPUTE THE TRANSFERT FUNCTION OF THE SYMMETRIC LINEAR FILTER AT THE FOURIER FREQUENCIES.
!
        call freq_func( m, coef(:k), freqr(:m), four_freq=true )
!
!       NOW, APPLY THE FILTER USING THE FFT AND THE CONVOLUTION THEOREM.
!
!       INITIALIZE THE FFT SUBROUTINE.
!
        if ( initfft2 ) then
            call init_fft( (/ n, m /), dim=2_i4b )
        end if
!
!       TRANSFORM THE TIME SERIES.
!
        matc(:n,:m) = cmplx( mat(:n,:m), zero, kind=stnd )
!
        call fft_row( matc(:n,:m), forward=true  )
!
!       MULTIPLY THE FOURIER TRANSFORM OF THE TIME SERIES
!       BY THE TRANSFERT FUNCTION OF THE FILTER.
!
        do i = 1_i4b, m
!
            matc(:n,i) = matc(:n,i)*freqr(i)
!
        end do
!
!       INVERT THE SEQUENCE BACK TO GET THE FILTERED TIME SERIES.
!
        call fft_row( matc(:n,:m), forward=false )
!
        mat(:n,:m) = real( matc(:n,:m),  kind=stnd )
!
!       DEALLOCATE WORK ARRAYS.
!
        if ( initfft2 ) then
            call end_fft( )
        end if
!
        deallocate( matc, freqr )
!
    else
!
!       ALLOCATE WORK ARRAY.
!
        allocate( mat2(n,m), stat=iok )
!
        if ( iok/=0 ) then
            call merror( name_proc//allocate_error )
        end if
!
!       MAKE A COPY OF THE MULTICHANNEL TIME SERIES.
!
        mat2(:n,:m) = mat(:n,:m)
!
!       NOW, FILTER THE TIME SERIES.
!
        kmid  = ( k + 1_i4b )/2_i4b
!   
#ifdef _OPENMP
        i1 = omp_get_num_procs()
        i2 = omp_get_max_threads()
        test_par = .not.( omp_in_parallel() )   .and.      &
                   m*n>=omp_limit               .and.      &
                   i1>1_i4b                     .and.      &
                   i2>1_i4b
!
        if ( test_par ) then       
!
!$OMP PARALLEL PRIVATE(i,i1,i2,vec,coef2)                     &
!$OMP         ,SHARED(m,n,k,kmid,khalf,mat,mat2,coef)
!
            coef2(:k) = coef(:k)
!
!$OMP DO SCHEDULE(STATIC) 
!
            do i = 1_i4b, khalf
!
!                i1 = kmid + 1_i4b - i
!                i2 = khalf + i
!
!                mat(:n,i) = matmul( mat2(:n,1_i4b:i2), coef2(i1:k) )
!
!BUG: For compiling STATPACK with OPENMP and PGI compiler
!     replace the preceding section with
!
                i1 = kmid  - i
                i2 = khalf + i
!
                vec(:n) = zero
!
                do j=1_i4b, i2
!
                    i1 = i1 + 1_i4b
                    vec(:n) = vec(:n) + coef2(i1)*mat2(:n,j)
!
                end do
!
                mat(:n,i) = vec(:n)
!
            end do
!
!$OMP END DO NOWAIT
!
!
!$OMP DO SCHEDULE(STATIC) 
!
            do i = kmid, m-khalf
!
!                i1 = i - kmid + 1_i4b
!                i2 = k  + i1  - 1_i4b
!
!                mat(:n,i) = matmul( mat2(:n,i1:i2), coef2(1_i4b:k) )
!
!BUG: For compiling STATPACK with OPENMP and PGI compiler
!     replace the preceding section with
!
                i1 = i - kmid
!
                vec(:n) = zero
!
                do j=1_i4b, k
!
                    i1 = i1 + 1_i4b
                    vec(:n) = vec(:n) + coef2(j)*mat2(:n,i1)
!
                end do
!
                mat(:n,i) = vec(:n)
!
            end do
!
!$OMP END DO NOWAIT
!
!
!$OMP DO SCHEDULE(STATIC) 
!
            do i =  m-khalf+1_i4b, m
!
!                i1 = k - khalf - i + m
!                i2 = m - i1 + 1_i4b
!
!                mat(:n,i) = matmul( mat2(:n,i2:m), coef2(1_i4b:i1) )
!
!BUG: For compiling STATPACK with OPENMP and PGI compiler
!     replace the preceding section with
!
                i1 = k - khalf - i + m
                i2 = m - i1
!
                vec(:n) = zero
!
                do j=1_i4b, i1
!
                    i2 = i2 + 1_i4b
                    vec(:n) = vec(:n) + coef2(j)*mat2(:n,i2)
!
                end do
!
                mat(:n,i) = vec(:n)
!
            end do
!
!$OMP END DO
!
!$OMP END PARALLEL
!
        else
!
#endif
!
            i1 = kmid + 1_i4b
            i2 = khalf
!
            do i = 1_i4b, khalf
!
                i1 = i1 - 1_i4b
                i2 = i2 + 1_i4b
!
                mat(:n,i) = matmul( mat2(:n,:i2), coef(i1:k) )
!
            end do
!
            i1 = 0_i4b
            i2 = k - 1_i4b
!
            do i = kmid, m-khalf
!
                i1 = i1 + 1_i4b
                i2 = i2 + 1_i4b
!
                mat(:n,i) = matmul( mat2(:n,i1:i2), coef(:k) )
!
            end do
!
            i1 = k
            i2 = m - k + 1_i4b
!
            do i =  m-khalf+1_i4b, m
!
                i1 = i1 - 1_i4b
                i2 = i2 + 1_i4b
!
                mat(:n,i) = matmul( mat2(:n,i2:m), coef(:i1) )
!
            end do
!
#ifdef _OPENMP
        end if
#endif
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( mat2 )
!
    end if
!
!   ADD MEANS, DRIFTS OR LINEAR LEAST SQUARES LINES FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEANS.
!
            do i = 1_i4b, m
                mat(:n,i) = mat(:n,i) + orig(:n)
            end do
!
        case( -2_i4b )
!
!           ADD DRIFTS.
!
            do i = 1_i4b, m
                tmp       = real( i, stnd ) - one
                mat(:n,i) = mat(:n,i) + tmp*slope(:n)           
            end do
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINES.
!         
            do i = 1_i4b, m
                tmp       = real( i, stnd ) - one
                mat(:n,i) = mat(:n,i) + tmp*slope(:n) + orig(:n)      
            end do
!
        end select
!
!
! END OF SUBROUTINE symlin_filter2_rm
! ___________________________________
!
    end subroutine symlin_filter2_rm
!
! =========================================================================================
!
    subroutine dan_filter_rv( vec, nsmooth, sym, trend )
!
! Purpose
! _______
!
!   Subroutine DAN_FILTER smooths an input time series (e.g. the argument VEC) by applying a
!   Daniell filter (e.g. a simple moving average) of length NSMOOTH.
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On input, the vector containing the time series to be filtered.
!                 On output, the filtered time series is returned.
!
!                 Size(VEC) must be greater or equal to 4.
!
!   NSMOOTH       (INPUT) integer(i4b)
!                 The length of the  Daniell filter to be applied to the time series. NSMOOTH must be odd.
!                 Moreover, NSMOOTH must be greater or equal to 3 and less or equal to size(VEC).
!
!   SYM           (INPUT, OPTIONAL) real(stnd)
!                 An optional indictor variable used to designate wether the series has an even symmetry
!                 (SYM = one), an odd symmetry (SYM = -one) or no symmetry (SYM = zero). Other values
!                 than -one, one or zero are not allowed for the optional argument SYM.
!
!                 The default value for SYM is one.
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+/-1 The mean of the time series is removed before time filtering
!                 - TREND=+/-2 The drift from the time series is removed before time filtering
!                   by using the formula: drift = (VEC(size(VEC)) - VEC(1))/(size(VEC) - 1)
!                 - TREND=+/-3 The least-squares line from the time series is removed before
!                   time filtering.
!
!                 IF TREND=-1,-2 or -3, the mean, drift or least-squares line is reintroduced
!                 post-filtering, respectively.
!                 For other values of TREND nothing is done before or after filtering.
!
!
! Further Details
! _______________
!
!   Subroutine DAN_FILTER smooths an input time series by applying a Daniell filter
!   as discussed in chapter 7 of Bloomfield (1976).
!
!   This subroutine use the hypothesis of the (even or odd) symmetry of the input time series
!   to avoid losing values from the ends of the series.  
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P.,1976:
!             Fourier analysis of time series- An introduction, 
!             John Wiley and Sons, New York, Chapter 7.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror, arth
    use Reals_Constants,   only : zero, one
    use Char_Constants,    only : tseries_error10, tseries_error20, tseries_error60, tseries_error61
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:) :: vec
    real(stnd),   intent(in), optional        :: sym
!
    integer(i4b), intent(in)           :: nsmooth
    integer(i4b), intent(in), optional :: trend
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, i, j, j1, j2, lim, trend2
!
    real(stnd)                       :: sym2, con, orig, slope, veci
    real(stnd), dimension(size(vec)) :: vec2
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    logical      :: test_par
#endif
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='dan_filter'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m = size( vec )
!
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
!   CHECK THE INPUT VALUE FOR THE LENGTH OF THE DANIELL FILTER.
!
    if ( nsmooth<3_i4b .or. nsmooth>m )      &
    call merror( name_proc//tseries_error60 )
!
    if ( (nsmooth/2_i4b)*2_i4b==nsmooth )   &
    call merror( name_proc//tseries_error61 )
!
    sym2 = one
!
    if ( present(sym) ) then
!
        if ( abs(sym)/=one .and. sym/=zero )   &
        call merror( name_proc//tseries_error20 )
!
        sym2 = sym
!
    end if
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
!
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    i = abs( trend2 )
!
    if ( i>=1_i4b .and. i<=3_i4b ) then
        call detrend_rv( vec(:m), i, orig, slope )
    end if
!
!   NOW, FILTER THE TIME SERIES .
!
    vec2(:m) = vec(:m)
!
    lim = ( nsmooth - 1_i4b )/2_i4b
    con = one/real( nsmooth , stnd )
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    j1 = omp_get_num_procs()
    j2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               m*nsmooth>=omp_limit         .and.      &
               j1>1_i4b                     .and.      &
               J2>1_i4b
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(i,j,j1,j2,veci)                     &
!$OMP        , SHARED(m,lim,con,sym2,vec,vec2)
!
!$OMP DO SCHEDULE(STATIC) 
!
        do i = 1_i4b, lim
!
            veci = vec2(i)
!
            do j = 1_i4b, lim
!
                j1 = i - j
                j2 = i + j
!
                veci = veci + extd_rv( j1, sym2 ) + vec2(j2)
!
            end do
!
            vec(i) = veci*con
!
        end do
!
!$OMP END DO
!
!$OMP DO SCHEDULE(STATIC) 
!
        do i = lim+1_i4b, m-lim
!
            j1 = i - lim
            j2 = i + lim
!
            vec(i) = sum( vec2(j1:j2) )*con
!
        end do
!
!$OMP END DO
!
!$OMP DO SCHEDULE(STATIC) 
!
        do i = m-lim+1_i4b, m
!
            veci = vec2(i)
!
            do j = 1_i4b, lim
!
                j1 = i - j
                j2 = i + j
!
                veci = veci + vec2(j1) + extd_rv( j2, sym2 )
!
            end do
!
            vec(i) = veci*con
!
        end do
!
!$OMP END DO
!
!
! $OMP DO SCHEDULE(STATIC) 
!
!        do i = 1_i4b, m
!
!            veci = vec2(i)
!
!            do j = 1_i4b, lim
!
!                j1 = i - j
!                j2 = i + j
!
!                if ( j1>=1_i4b .and. j2<=m ) then
!                    veci = veci + vec2(j1) + vec2(j2)
!                else
!                    veci = veci + extd_rv( j1, sym2 ) +  extd_rv( j2, sym2 )
!                end if
!
!            end do
!
!            vec(i) = veci*con
!
!        end do
!
! $OMP END DO
!
!$OMP END PARALLEL
!
    else
!
#endif
#endif
!
        do i = 1_i4b, lim
!
            veci = vec2(i)
!
            do j = 1_i4b, lim
!
                j1 = i - j
                j2 = i + j
!
                veci = veci + extd_rv( j1, sym2 ) + vec2(j2)
!
            end do
!
            vec(i) = veci*con
!
        end do
!
        do i = lim+1_i4b, m-lim
!
            j1 = i - lim
            j2 = i + lim
!
            vec(i) = sum( vec2(j1:j2) )*con
!
        end do
!
        do i = m-lim+1_i4b, m
!
            veci = vec2(i)
!
            do j = 1_i4b, lim
!
                j1 = i - j
                j2 = i + j
!
                veci = veci + vec2(j1) + extd_rv( j2, sym2 )
!
            end do
!
            vec(i) = veci*con
!
        end do
!
!        do i = 1_i4b, m
!
!            veci = vec2(i)
!
!            do j = 1_i4b, lim
!
!                j1 = i - j
!                j2 = i + j
!
!                if ( j1>=1_i4b .and. j2<=m ) then
!                    veci = veci + vec2(j1) + vec2(j2)
!                else
!                    veci = veci + extd_rv( j1, sym2 ) + extd_rv( j2, sym2 )
!                end if
!
!            end do
!
!            vec(i) = veci*con
!
!        end do
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    end if
#endif
#endif
!
!   ADD MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEAN.
!
            vec(:m) = vec(:m) + orig
!
        case( -2_i4b )
!
!           ADD DRIFT.
!         
            vec(:m) = vec(:m) + slope*arth( zero, one, m )           
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINE.
!         
            vec(:m) = vec(:m) + ( orig + slope*arth( zero, one, m ) )       
!         
    end select
!
!----------------------------------------------------------------------
                         contains
!----------------------------------------------------------------------
!
    function extd_rv( index, sym ) result( vecj )
!
! Purpose
! _______
!
!   This internal function returns the INDEX-th term in the series VEC,
!   extending it if necessary with even or odd symmetry according
!   to the sign of SYM, which should be either plus or minus one.
!   (Note: the value zero will result in the extended value being zero).
!
! 
! Arguments
! _________
!
!   INDEX         (INPUT) integer(i4b)
!                 On input, the index of the desired term in the time series.
!                 INDEX may be any integer.
!
!   SYM           (INPUT) real(stnd)
!                 An indictor variable used to designate wether the series has an even symmetry
!                 (SYM = one), an odd symmetry (SYM = -one) or no symmetry (SYM = zero). Other values
!                 than -one, one or zero are not allowed, however no checking is done on the SYM
!                 argument.
!
!
! Further Details
! _______________
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: index
!
    real(stnd), intent(in) :: sym
!
    real(stnd)  :: vecj
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: tmp
!
    integer(i4b) :: index2
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    index2 = index
    tmp    = one
!
    do
!
        if ( index2<1_i4b ) then
            index2 = 2_i4b - index2
            tmp    = tmp*sym
        end if
!
        if ( index2<=m ) exit
!
        index2 = 2_i4b*m - index2
        tmp    = tmp*sym
!
    end do
!
    vecj = tmp*vec2(index2)
!
!
! END OF FUNCTION extd_rv
! _______________________
!
    end function extd_rv
!
!----------------------------------------------------------------------
!
! END OF SUBROUTINE dan_filter_rv
! _______________________________
!
    end subroutine dan_filter_rv
!
! =========================================================================================
!
    subroutine dan_filter_rm( mat, nsmooth, sym, trend )
!
! Purpose
! _______
!
!   Subroutine DAN_FILTER smooths an input multi-channel time series (the argument MAT) by
!   applying a Daniell filter (e.g. a simple moving average) of length NSMOOTH.
!
!
! Arguments
! _________
!
!   MAT           (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                 The multi-channel time series matrix to be filtered.
!                 Each column of MAT corresponds to one observation.
!                 On output, the multi-channel filtered time series are returned.
!
!                 Size(MAT,2) must be greater or equal to 4.
!
!   NSMOOTH       (INPUT) integer(i4b)
!                 The length of the  Daniell filter to be applied to the time series. NSMOOTH must be odd.
!                 Moreover, NSMOOTH must be greater or equal to 3 and less or equal to size(MAT,2).
!
!   SYM           (INPUT, OPTIONAL) real(stnd)
!                 An optional indictor variable used to designate wether the series has an even symmetry
!                 (SYM = one), an odd symmetry (SYM = -one) or no symmetry (SYM = zero). Other values
!                 than -one, one or zero are not allowed for the optional argument SYM.
!
!                 The default value for SYM is one.
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+/-1 The means of the time series are removed before time filtering
!                 - TREND=+/-2 The drifts from the time series are removed before time filtering
!                   by using the formula: drift(:) = (MAT(:,size(MAT,2)) - MAT(:,1))/(size(MAT,2) - 1)
!                 - TREND=+/-3 The least-squares lines from the time series are removed before
!                   time filtering.
!
!                 IF TREND=-1,-2 or -3, the means, drifts or least-squares lines are reintroduced
!                 post-filtering, respectively.
!                 For other values of TREND nothing is done before or after filtering.
!
!
! Further Details
! _______________
!
!   Subroutine DAN_FILTER smooths an input multi-channel time series by applying a
!   Daniell filter as discussed in chapter 7 of Bloomfield (1976).
!
!   This subroutine may use the hypothesis of the (even or odd) symmetry of the input time series
!   to avoid losing values from the ends of the series.  
!
!   For more details and algorithm, see
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 7.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Reals_Constants,   only : zero, one
    use Char_Constants,    only : tseries_error10, tseries_error20, tseries_error61, tseries_error62
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:,:)  :: mat
    real(stnd),   intent(in), optional           :: sym
!
    integer(i4b), intent(in)           :: nsmooth
    integer(i4b), intent(in), optional :: trend
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, m, i, j, j1, j2, k, lim, trend2
!
    real(stnd)                                     :: sym2, con
    real(stnd), dimension(size(mat,1))             :: orig, slope, mati
    real(stnd), dimension(size(mat,1),size(mat,2)) :: mat2
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    logical      :: test_par
#endif
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='dan_filter'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( mat, 1 )
    if ( n<=0_i4b )  return
!
    m = size( mat, 2 )
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
!   CHECK THE INPUT VALUE FOR THE LENGTH OF THE DANIEL FILTER.
!
    if ( nsmooth<3_i4b .or. nsmooth>m )     &
    call merror( name_proc//tseries_error62 )
!
    if ( (nsmooth/2_i4b)*2_i4b==nsmooth )   &
    call merror( name_proc//tseries_error61 )
!
    sym2 = one
    if ( present(sym) ) then
!
        if ( abs(sym)/=one .and. sym/=zero )   &
        call merror( name_proc//tseries_error20 )
!
        sym2 = sym
!
    end if
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    i = abs( trend2 )
    if ( i>=1_i4b .and. i<=3_i4b ) then
        call detrend_rm( mat(:n,:m), i, orig(:n), slope(:n) )
    end if
!
!   NOW, FILTER THE TIME SERIES .
!
    mat2(:n,:m) = mat(:n,:m)
!
    lim = ( nsmooth - 1_i4b )/2_i4b
    con = one/real( nsmooth , stnd )
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    j1 = omp_get_num_procs()
    j2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               m*n*nsmooth>=omp_limit       .and.      &
               j1>1_i4b                     .and.      &
               J2>1_i4b
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(i,j,j1,j2,mati)                 &
!$OMP          , SHARED(m,n,lim,con,sym2,mat,mat2)
!
!
!$OMP DO SCHEDULE(STATIC) 
!
        do i = 1_i4b, lim
!
            mati(:n) = mat2(:n,i)
!
            do j = 1_i4b, lim
!
                j1 = i - j
                j2 = i + j
!
                mati(:n) = mati(:n) + extd_rm( j1, sym2, n ) + mat2(:n,j2)
!
            end do
!
            mat(:n,i) = mati(:n)*con
!
        end do
!
!$OMP END DO
!
!$OMP DO SCHEDULE(STATIC) 
!
        do i = lim+1_i4b, m-lim
!
            j1 = i - lim
            j2 = i + lim
!
            mat(:n,i) = sum( mat2(:n,j1:j2), dim=2 )*con
!
        end do
!
!$OMP END DO
!
!$OMP DO SCHEDULE(STATIC) 
!
        do i = m-lim+1_i4b, m
!
            mati(:n) = mat2(:n,i)
!
            do j = 1_i4b, lim
!
                j1 = i - j
                j2 = i + j
!
                mati(:n) = mati(:n) + mat2(:n,j1) + extd_rm( j2, sym2, n )
!
            end do
!
            mat(:n,i) = mati(:n)*con
!
        end do
!
!$OMP END DO
!
!$OMP END PARALLEL
!
    else
!
#endif
#endif
!
        do i = 1_i4b, lim
!
            mati(:n) = mat2(:n,i)
!
            do j = 1_i4b, lim
!
                j1 = i - j
                j2 = i + j
!
                mati(:n) = mati(:n) + extd_rm( j1, sym2, n ) + mat2(:n,j2)
!
            end do
!
            mat(:n,i) = mati(:n)*con
!
        end do
!
        do i = lim+1_i4b, m-lim
!
            j1 = i - lim
            j2 = i + lim
!
            mat(:n,i) = sum( mat2(:n,j1:j2), dim=2 )*con
!
        end do
!
        do i = m-lim+1_i4b, m
!
            mati(:n) = mat2(:n,i)
!
            do j = 1_i4b, lim
!
                j1 = i - j
                j2 = i + j
!
                mati(:n) = mati(:n) + mat2(:n,j1) + extd_rm( j2, sym2, n )
!
            end do
!
            mat(:n,i) = mati(:n)*con
!
        end do
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    end if
#endif
#endif
!
!   ADD MEANS, DRIFTS OR LINEAR LEAST SQUARES LINES FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEANS.
!
            do k = 1_i4b, m
                mat(:n,k) = mat(:n,k) + orig(:n)
            end do
!
        case( -2_i4b )
!
!           ADD DRIFTS.
!
            do k = 1_i4b, m
                con       = real( k, stnd ) - one
                mat(:n,k) = mat(:n,k) + con*slope(:n)           
            end do
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINES.
!         
            do k = 1_i4b, m
                con       = real( k, stnd ) - one
                mat(:n,k) = mat(:n,k) + con*slope(:n) + orig(:n)      
            end do
!
    end select
!
!----------------------------------------------------------------------
                         contains
!----------------------------------------------------------------------
!
    function extd_rm( index, sym, nt ) result( matj )
!
! Purpose
! _______
!
!   This internal function returns the INDEX-th term in the multi-channel time series MAT,
!   extending it if necessary with even or odd symmetry according
!   to the sign of SYM, which should be either plus or minus one.
!   (Note: the value zero will result in the extended value being zero).
!
! 
! Arguments
! _________
!
!   INDEX         (INPUT) integer(i4b)
!                 On input, the index of the desired term in the multi-channel time series.
!                 INDEX may be any integer.
!
!   SYM           (INPUT) real(stnd)
!                 An indictor variable used to designate wether the series has an even symmetry
!                 (SYM = one), an odd symmetry (SYM = -one) or no symmetry (SYM = zero). Other values
!                 than -one, one or zero are not allowed, however no checking is done on the SYM
!                 argument.
!
!   NT            (INPUT) integer(i4b)
!                 On input, the number of time series to extend.
!
!
! Further Details
! _______________
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: index, nt
!
    real(stnd), intent(in) :: sym
!
    real(stnd), dimension(nt)  :: matj
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: tmp
!
    integer(i4b) :: index2
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    index2 = index
    tmp    = one
!
    do
!
        if ( index2<1_i4b ) then
            index2 = 2_i4b - index2
            tmp    = tmp*sym
        end if
!
        if ( index2<=m ) exit
!
        index2 = 2_i4b*m - index2
        tmp    = tmp*sym
!
    end do
!
    matj(:nt) = tmp*mat2(:nt,index2)
!
!
! END OF FUNCTION extd_rm
! _______________________
!
    end function extd_rm
!
!----------------------------------------------------------------------
!
!
! END OF SUBROUTINE dan_filter_rm
! _______________________________
!
    end subroutine dan_filter_rm
!
! =========================================================================================
!
    subroutine moddan_filter_rv( vec, smooth_param, sym, trend )
!
! Purpose
! _______
!
!   Subroutine MODDAN_FILTER smooths an input time series (e.g. the argument VEC) by
!   applying a sequence of modified Daniell filters.
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On input, the vector containing the time series to be filtered.
!                 On output, the filtered time series is returned.
!                 Size(VEC) must be greater or equal to 4.
!
!   SMOOTH_PARAM  (INPUT) integer(i4b), dimension(:)
!                 The array of the half-lengths of the modified Daniell filters to be applied
!                 to the time series. All the values in SMOOTH_PARAM(:) must be greater than 0
!                 and less than size(VEC) .
!
!                 Size(SMOOTH_PARAM) must be greater or equal to 1.
!
!   SYM           (INPUT, OPTIONAL) real(stnd)
!                 An optional indictor variable used to designate wether the series has an even symmetry
!                 (SYM = one), an odd symmetry (SYM = -one) or no symmetry (SYM = zero). Other values
!                 than -one, one or zero are not allowed for the optional argument SYM.
!
!                 The default value for SYM is one.
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+/-1 The mean of the time series is removed before time filtering
!                 - TREND=+/-2 The drift from the time series is removed before time filtering
!                   by using the formula: drift = (VEC(size(VEC)) - VEC(1))/(size(VEC) - 1)
!                 - TREND=+/-3 The least-squares line from the time series is removed before
!                   time filtering.
!
!                 IF TREND=-1,-2 or -3, the mean, drift or least-squares line is reintroduced
!                 post-filtering, respectively.
!                 For other values of TREND nothing is done before or after filtering.
!
!
! Further Details
! _______________
!
!   Subroutine MODDAN_FILTER smooths an input time series by applying a sequence of
!   modified Daniell filters as discussed in chapter 7 of Bloomfield (1976). This
!   subroutine use the hypothesis of the (even or odd) symmetry of the input time series
!   to avoid losing values from the ends of the series.  
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 7.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror, arth
    use Reals_Constants,   only : zero, half, one
    use Char_Constants,    only : tseries_error10, tseries_error20, tseries_error52
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:) :: vec
    real(stnd),   intent(in), optional        :: sym
!
    integer(i4b), intent(in), dimension(:)  :: smooth_param
    integer(i4b), intent(in), optional      :: trend
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, nparam, l, k, i, j, j1, j2, lim, trend2
!
    real(stnd)                       :: sym2, con, orig, slope, veci
    real(stnd), dimension(size(vec)) :: vec2
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    logical      :: test_par
#endif
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='moddan_filter'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m = size( vec )
!
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    nparam = size( smooth_param )
!
    if ( nparam<=0_i4b  ) then
        return
    end if
!
!   CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIEL FILTERS.
!
    if ( any( smooth_param(:nparam)<=0_i4b .or. smooth_param(:nparam)>=m ) )   &
    call merror( name_proc//tseries_error52 )
!
    sym2 = one
!
    if ( present(sym) ) then
!
        if ( abs(sym)/=one .and. sym/=zero )   &
        call merror( name_proc//tseries_error20 )
!
        sym2 = sym
!
    end if
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
!
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    i = abs( trend2 )
!
    if ( i>=1_i4b .and. i<=3_i4b ) then
        call detrend_rv( vec(:m), i, orig, slope )
    end if
!
!   NOW, FILTER THE TIME SERIES.
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    j1 = omp_get_num_procs()
    j2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               m*nparam>=omp_limit          .and.      &
               j1>1_i4b                     .and.      &
               J2>1_i4b
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(l,k,lim,i,j,j1,j2,con,veci)              &
!$OMP          , SHARED(m,nparam,smooth_param,sym2,vec,vec2)
!
        do l = 1_i4b, nparam
!
            k   = smooth_param(l)
            lim = k - 1_i4b
            con = one/real( 2_i4b*k , stnd )
!
!$OMP SINGLE
!
            vec2(:m) = vec(:m)
!
!$OMP END SINGLE
!
!$OMP DO SCHEDULE(STATIC) 
!
            do i = 1_i4b, m
!
                veci = vec2(i)
!
                do j = 1_i4b, lim
!
                    j1 = i - j
                    j2 = i + j
!
                    if ( j1>=1_i4b .and. j2<=m ) then
                        veci = veci + vec2(j1) + vec2(j2)
                    else
                        veci = veci + extd_rv( j1, sym2 ) +     &
                                      extd_rv( j2, sym2 )
                    end if
!
                end do
!
                j1 = i - k
                j2 = i + k
!
                if ( j1>=1_i4b .and. j2<=m ) then
                    vec(i) = ( veci + ( vec2(j1) + vec2(j2) )*half )*con
                else
                    vec(i) = ( veci + ( extd_rv( j1, sym2 ) +            &
                                        extd_rv( j2, sym2 ) )*half )*con
                end if
!
            end do
!
!$OMP END DO
!
        end do
!
!$OMP END PARALLEL
!
    else
!
#endif
#endif
!
        do l = 1_i4b, nparam
!
            k   = smooth_param(l)
            lim = k - 1_i4b
            con = one/real( 2_i4b*k , stnd )
!
            vec2(:m) = vec(:m)
!
            do i = 1_i4b, m
!
                veci = vec2(i)
!
                do j = 1_i4b, lim
!
                    j1 = i - j
                    j2 = i + j
!
                    if ( j1>=1_i4b .and. j2<=m ) then
                        veci = veci + vec2(j1) + vec2(j2)
                    else
!                        veci = veci + extend_rv( vec2(:m), j1, sym2 ) +     &
!                                      extend_rv( vec2(:m), j2, sym2 )
                        veci = veci + extd_rv( j1, sym2 ) +     &
                                      extd_rv( j2, sym2 )
                    end if
!
                end do
!
                j1 = i - k
                j2 = i + k
!
                if ( j1>=1_i4b .and. j2<=m ) then
                    vec(i) = ( veci + ( vec2(j1) + vec2(j2) )*half )*con
                else
!                    vec(i) = ( veci + ( extend_rv( vec2(:m), j1, sym2 ) +            &
!                                        extend_rv( vec2(:m), j2, sym2 ) )*half )*con
                    vec(i) = ( veci + ( extd_rv( j1, sym2 ) +            &
                                        extd_rv( j2, sym2 ) )*half )*con
                end if
!
            end do
!
        end do
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    end if
#endif
#endif
!
!   ADD MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEAN.
!
            vec(:m) = vec(:m) + orig
!
        case( -2_i4b )
!
!           ADD DRIFT.
!         
            vec(:m) = vec(:m) + slope*arth( zero, one, m )           
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINE.
!         
            vec(:m) = vec(:m) + ( orig + slope*arth( zero, one, m ) )       
!         
    end select
!
!----------------------------------------------------------------------
                         contains
!----------------------------------------------------------------------
!
    function extd_rv( index, sym ) result( vecj )
!
! Purpose
! _______
!
!   This internal function returns the INDEX-th term in the series VEC,
!   extending it if necessary with even or odd symmetry according
!   to the sign of SYM, which should be either plus or minus one.
!   (Note: the value zero will result in the extended value being zero).
!
! 
! Arguments
! _________
!
!   INDEX         (INPUT) integer(i4b)
!                 On input, the index of the desired term in the time series.
!                 INDEX may be any integer.
!
!   SYM           (INPUT) real(stnd)
!                 An indictor variable used to designate wether the series has an even symmetry
!                 (SYM = one), an odd symmetry (SYM = -one) or no symmetry (SYM = zero). Other values
!                 than -one, one or zero are not allowed, however no checking is done on the SYM
!                 argument.
!
!
! Further Details
! _______________
!
!   For more details and algorithm, see
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: index
!
    real(stnd), intent(in) :: sym
!
    real(stnd)  :: vecj
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: tmp
!
    integer(i4b) :: index2
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    index2 = index
    tmp    = one
!
    do
!
        if ( index2<1_i4b ) then
            index2 = 2_i4b - index2
            tmp    = tmp*sym
        end if
!
        if ( index2<=m ) exit
!
        index2 = 2_i4b*m - index2
        tmp    = tmp*sym
!
    end do
!
    vecj = tmp*vec2(index2)
!
!
! END OF FUNCTION extd_rv
! _______________________
!
    end function extd_rv
!
!----------------------------------------------------------------------
!
! END OF SUBROUTINE moddan_filter_rv
! __________________________________
!
    end subroutine moddan_filter_rv
!
! =========================================================================================
!
    subroutine moddan_filter_rm( mat, smooth_param, sym, trend )
!
! Purpose
! _______
!
!   Subroutine MODDAN_FILTER smooths an input multi-channel time series (the argument MAT) by
!   applying a sequence of modified Daniell filters.
!
!
! Arguments
! _________
!
!   MAT           (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                 The multi-channel time series matrix to be filtered.
!                 Each column of MAT corresponds to one observation.
!                 On output, the multi-channel filtered time series are returned.
!
!                 Size(MAT,2) must be greater or equal to 4.
!
!   SMOOTH_PARAM  (INPUT) integer(i4b), dimension(:)
!                 The array of the half-lengths of the modified Daniell filters to be applied
!                 to the time series. All the values in SMOOTH_PARAM(:) must be greater than 0
!                 and less than size(MAT,2) .
!
!                 Size(SMOOTH_PARAM) must be greater or equal to 1.
!
!   SYM           (INPUT, OPTIONAL) real(stnd)
!                 An optional indictor variable used to designate wether the series has an even symmetry
!                 (SYM = one), an odd symmetry (SYM = -one) or no symmetry (SYM = zero). Other values
!                 than -one, one or zero are not allowed for the optional argument SYM.
!
!                 The default value for SYM is one.
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+/-1 The means of the time series are removed before time filtering
!                 - TREND=+/-2 The drifts from the time series are removed before time filtering
!                   by using the formula: drift(:) = (MAT(:,size(MAT,2)) - MAT(:,1))/(size(MAT,2) - 1)
!                 - TREND=+/-3 The least-squares lines from the time series are removed before
!                   time filtering.
!
!                 IF TREND=-1,-2 or -3, the means, drifts or least-squares lines are reintroduced
!                 post-filtering, respectively.
!                 For other values of TREND nothing is done before or after filtering.
!
!
! Further Details
! _______________
!
!   Subroutine MODDAN_FILTER smooths an input multi-channel time series by applying a
!   sequence of modified Daniell filters as discussed in chapter 7 of Bloomfield (1976). This
!   subroutine may use the hypothesis of the (even or odd) symmetry of the input time series
!   to avoid losing values from the ends of the series.  
!
!   For more details and algorithm see
!
!   (1) Bloomfield, P., 1976: Fourier analysis of time series- An introduction, 
!              John Wiley and Sons, New York, Chapter 7.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Reals_Constants,   only : zero, half, one
    use Char_Constants,    only : tseries_error10, tseries_error20, tseries_error53
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:,:)  :: mat
    real(stnd),   intent(in), optional           :: sym
!
    integer(i4b), intent(in), dimension(:)  :: smooth_param
    integer(i4b), intent(in), optional      :: trend
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, m, nparam, l, k, i, j, j1, j2, lim, trend2
!
    real(stnd)                                     :: sym2, con
    real(stnd), dimension(size(mat,1))             :: orig, slope, mati
    real(stnd), dimension(size(mat,1),size(mat,2)) :: mat2
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    logical      :: test_par
#endif
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='moddan_filter'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( mat, 1 )
    if ( n<=0_i4b )  return
!
    m = size( mat, 2 )
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    nparam = size( smooth_param )
!
    if ( nparam<=0_i4b  ) then
        return
    end if
!
!   CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIEL FILTERS.
!
    if ( any( smooth_param(:nparam)<=0_i4b .or. smooth_param(:nparam)>=m ) )   &
    call merror( name_proc//tseries_error53 )
!
    sym2 = one
    if ( present(sym) ) then
!
        if ( abs(sym)/=one .and. sym/=zero )   &
        call merror( name_proc//tseries_error20 )
!
        sym2 = sym
!
    end if
!
!   NO TRANSFORMATION BEFORE OR AFTER FILTERING BY DEFAULT.
!
    trend2 = 0_i4b
    if ( present(trend) ) then
        if ( abs(trend)>=1_i4b .and. abs(trend)<=3_i4b ) trend2 = trend
    end if    
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF NEEDED.
!
    i = abs( trend2 )
    if ( i>=1_i4b .and. i<=3_i4b ) then
        call detrend_rm( mat(:n,:m), i, orig(:n), slope(:n) )
    end if
!
!   NOW, FILTER THE TIME SERIES.
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    j1 = omp_get_num_procs()
    j2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               m*n*nparam>=omp_limit        .and.      &
               j1>1_i4b                     .and.      &
               J2>1_i4b
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(l,k,lim,i,j,j1,j2,con,mati)              &
!$OMP          , SHARED(m,n,nparam,smooth_param,sym2,mat,mat2)
!
        do l = 1_i4b, nparam
!
            k   = smooth_param(l)
            lim = k - 1_i4b
            con = one/real( 2_i4b*k , stnd )
!
!$OMP SINGLE
!
            mat2(:n,:m) = mat(:n,:m)
!
!$OMP END SINGLE
!
!$OMP DO SCHEDULE(STATIC) 
!
            do i = 1_i4b, m
!
                mati(:n) = mat2(:n,i)
!
                do j = 1_i4b, lim
!
                    j1 = i - j
                    j2 = i + j
!
                    if ( j1>=1_i4b .and. j2<=m ) then
                        mati(:n) = mati(:n) + mat2(:n,j1) + mat2(:n,j2)
                    else
                        mati(:n) = mati(:n) + extd_rm( j1, sym2, n ) +     &
                                              extd_rm( j2, sym2, n )
                    end if
!
                end do
!
                j1 = i - k
                j2 = i + k
!
                if ( j1>=1_i4b .and. j2<=m ) then
                    mat(:n,i) = ( mati(:n) + ( mat2(:n,j1) + mat2(:n,j2) )*half )*con
                else
                    mat(:n,i) = ( mati(:n) + ( extd_rm( j1, sym2, n ) +              &
                                               extd_rm( j2, sym2, n ) )*half )*con
               end if
!
            end do
!
!$OMP END DO
!
        end do
!
!$OMP END PARALLEL
!
    else
!
#endif
#endif
!
        do l = 1_i4b, nparam
!
            k   = smooth_param(l)
            lim = k - 1_i4b
            con = one/real( 2_i4b*k , stnd )
!
            mat2(:n,:m) = mat(:n,:m)
!
            do i = 1_i4b, m
!
                mati(:n) = mat2(:n,i)
!
                do j = 1_i4b, lim
!
                    j1 = i - j
                    j2 = i + j
!
                    if ( j1>=1_i4b .and. j2<=m ) then
                        mati(:n) = mati(:n) + mat2(:n,j1) + mat2(:n,j2)
                    else
                        mati(:n) = mati(:n) + extd_rm( j1, sym2, n ) +     &
                                              extd_rm( j2, sym2, n )
                    end if
!
                end do
!
                j1 = i - k
                j2 = i + k
!
                if ( j1>=1_i4b .and. j2<=m ) then
                    mat(:n,i) = ( mati(:n) + ( mat2(:n,j1) + mat2(:n,j2) )*half )*con
                else
                    mat(:n,i) = ( mati(:n) + ( extd_rm( j1, sym2, n ) +              &
                                               extd_rm( j2, sym2, n ) )*half )*con
               end if
!
            end do
!
        end do
!
#ifndef _INTERNAL_PROC
#ifdef _OPENMP
    end if
#endif
#endif
!
!   ADD MEANS, DRIFTS OR LINEAR LEAST SQUARES LINES FROM THE SERIES IF NEEDED.
!
    select case( trend2 )
!
        case( -1_i4b )
!
!           ADD MEANS.
!
            do k = 1_i4b, m
                mat(:n,k) = mat(:n,k) + orig(:n)
            end do
!
        case( -2_i4b )
!
!           ADD DRIFTS.
!
            do k = 1_i4b, m
                con       = real( k, stnd ) - one
                mat(:n,k) = mat(:n,k) + con*slope(:n)           
            end do
!
        case( -3_i4b )
!
!           ADD LINEAR LEAST SQUARES LINES.
!         
            do k = 1_i4b, m
                con       = real( k, stnd ) - one
                mat(:n,k) = mat(:n,k) + con*slope(:n) + orig(:n)      
            end do
!
    end select
!
!----------------------------------------------------------------------
                         contains
!----------------------------------------------------------------------
!
    function extd_rm( index, sym, nt ) result( matj )
!
! Purpose
! _______
!
!   This internal function returns the INDEX-th term in the multi-channel time series MAT,
!   extending it if necessary with even or odd symmetry according
!   to the sign of SYM, which should be either plus or minus one.
!   (Note: the value zero will result in the extended value being zero).
!
! 
! Arguments
! _________
!
!   INDEX         (INPUT) integer(i4b)
!                 On input, the index of the desired term in the multi-channel time series.
!                 INDEX may be any integer.
!
!   SYM           (INPUT) real(stnd)
!                 An indictor variable used to designate wether the series has an even symmetry
!                 (SYM = one), an odd symmetry (SYM = -one) or no symmetry (SYM = zero). Other values
!                 than -one, one or zero are not allowed, however no checking is done on the SYM
!                 argument.
!
!   NT            (INPUT) integer(i4b)
!                 On input, the number of time series to extend.
!
!
! Further Details
! _______________
!
!   For more details and algorithm, see
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: index, nt
!
    real(stnd), intent(in) :: sym
!
    real(stnd), dimension(nt)  :: matj
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: tmp
!
    integer(i4b) :: index2
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    index2 = index
    tmp    = one
!
    do
!
        if ( index2<1_i4b ) then
            index2 = 2_i4b - index2
            tmp    = tmp*sym
        end if
!
        if ( index2<=m ) exit
!
        index2 = 2_i4b*m - index2
        tmp    = tmp*sym
!
    end do
!
    matj(:nt) = tmp*mat2(:nt,index2)
!
!
! END OF FUNCTION extd_rm
! _______________________
!
    end function extd_rm
!
!----------------------------------------------------------------------
!
!
! END OF SUBROUTINE moddan_filter_rm
! __________________________________
!
    end subroutine moddan_filter_rm
!
! =========================================================================================
!
    function extend_rv( vec, index, sym ) result( vecj )
!
! Purpose
! _______
!
!   This function returns the INDEX-th term in the series VEC,
!   extending it if necessary with even or odd symmetry according
!   to the sign of SYM, which should be either plus or minus one.
!   (Note: the value zero will result in the extended value being zero).
!
! 
! Arguments
! _________
!
!   VEC           (INPUT) real(stnd), dimension(:)
!                 On input, the vector containing the time series.
!                 If size(VEC) is zero, the extended value returned is zero.
!
!   INDEX         (INPUT) integer(i4b)
!                 On input, the index of the desired term in the time series.
!                 INDEX may be any integer.
!
!   SYM           (INPUT) real(stnd)
!                 An indictor variable used to designate wether the series has an even symmetry
!                 (SYM = one), an odd symmetry (SYM = -one) or no symmetry (SYM = zero). Other values
!                 than -one, one or zero are not allowed, however no checking is done on the SYM
!                 argument.
!
!
! Further Details
! _______________
!
!   For more details and algorithm, see
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Reals_Constants,   only : zero, one
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: index
!
    real(stnd), intent(in)               :: sym
    real(stnd), intent(in), dimension(:) :: vec
!
    real(stnd)  :: vecj
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: con
!
    integer(i4b) :: m, j
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m = size( vec )
!
    if ( m<=0_i4b ) then
        vecj = zero
    else
!
        j   = index
        con = one
!
        do
!
            if ( j<1_i4b ) then
                j   = 2_i4b - j
                con = con*sym
            end if
!
            if ( j<=m ) exit
!
            j   = 2_i4b*m - j
            con = con*sym
!
        end do
!
        vecj = con*vec(j)
!
    end if
!
!
! END OF FUNCTION extend_rv
! _________________________
!
    end function extend_rv
!
! =========================================================================================
!
    function extend_rm( mat, index, sym ) result( matj )
!
! Purpose
! _______
!
!   This function returns the INDEX-th term in the multi-channel series MAT,
!   extending it if necessary with even or odd symmetry according
!   to the sign of SYM, which should be either plus or minus one.
!   Note: the value zero will result in the extended value being zero.
!
! 
! Arguments
! _________
!
!   MAT           (INPUT) real(stnd), dimension(:,:)
!                 On input, the matrix containing the multi-channel time series.
!                 Each column of MAT corresponds to one observation.
!                 If size(MAT,2) is zero, the extended vector (which is dimensionned as
!                 size(MAT,1)) returned is zero.
!
!   INDEX         (INPUT) integer(i4b)
!                 On input, the index of the desired term in the multi-channel time series.
!
!   SYM           (INPUT) real(stnd)
!                 An indictor variable used to designate wether the series has an even symmetry
!                 (SYM = one), an odd symmetry (SYM = -one) or no symmetry (SYM = zero). Other values
!                 than -one, one or zero are not allowed, however no checking is done on the SYM
!                 argument.
!
!
! Further Details
! _______________
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 6.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Reals_Constants,   only : zero, one
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: index
!
    real(stnd), intent(in)                 :: sym
    real(stnd), intent(in), dimension(:,:) :: mat
!
    real(stnd), dimension(size(mat,1))  :: matj
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: con
!
    integer(i4b) :: n, m, j
!

! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
    if ( m<=0_i4b ) then
!
        matj(:n) = zero
!
    else
!
        j   = index
        con = one
!
        do
!
            if ( j<1_i4b ) then
                j   = 2_i4b - j
                con = con*sym
            end if
!
            if ( j<=m ) exit
!
            j   = 2_i4b*m - j
            con = con*sym
!
        end do
!
        matj(:n) = con*mat(:n,j)
!
    end if
!
!
! END OF FUNCTION extend_rm
! _________________________
!
    end function extend_rm
!
! =========================================================================================
!
    subroutine taper_rv( vec, taperp )
!
! Purpose
! _______
!
!   Subroutine TAPER applies a split-cosine-bell taper on an input
!   time series (e.g. the argument VEC).
!
! Arguments
! _________
!
!   VEC          (INPUT/OUTPUT) real(stnd), dimension(:)
!                On input, the vector containing the time series to be tapered.
!                On output, the tapered time series is returned.
!
!   TAPERP       (INPUT) real(stnd)
!                The total percentage of the data to be tapered.
!                TAPERP must be greater than zero and less or equal to one,
!                otherwise the series is not tapered.
!
!
! Further Details
! _______________
!
!   This subroutine is adapted from Bloomfield (1976).
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 5.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Reals_Constants,   only : zero, one, half, pi
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout), dimension(:)  :: vec
!
    real(stnd), intent(in) :: taperp
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, n, i, j
!
    real(stnd)   :: con, weight
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m = size( vec )
!
    if ( (m<=0_i4b) .or. (taperp<=zero) .or. (taperp>one) ) return
!
    n    = int( taperp*real( m, stnd ) + half, i4b )/2_i4b
!
    if ( n>0_i4b ) then
!
        con  = one/real( n, stnd )
!
        do i = 1_i4b, n
            weight = half - half*cos( pi*(real(i,stnd)-half)*con )
            vec(i) = weight*vec(i)
            j = m + 1_i4b - i
            vec(j) = weight*vec(j)
        end do
!
    end if
!
! END OF SUBROUTINE taper_rv
! __________________________
!
    end subroutine taper_rv
!
! =========================================================================================
!
    subroutine taper_rm( mat, taperp )
!
! Purpose
! _______
!
!   Subroutine TAPER applies a split-cosine-bell taper on an input
!   multi-channel time series (the argument MAT).
!
!
! Arguments
! _________
!
!   MAT           (INPUT) real(stnd), dimension(:,:)
!                 On input, the matrix containing the multi-channel time series.
!                 Each column of MAT corresponds to one observation.
!
!   TAPERP       (INPUT) real(stnd)
!                The total percentage of the data to be tapered.
!                TAPERP must be greater than zero and less or equal to one,
!                otherwise the series is not tapered.
!
!
! Further Details
! _______________
!
!   This subroutine is adapted from Bloomfield (1976).
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 5.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Reals_Constants,   only : zero, one, half, pi
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout), dimension(:,:)  :: mat
!
    real(stnd), intent(in) :: taperp
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: p, m, n, i, j
!
    real(stnd)   :: con, weight
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    p = size( mat, 1 )
    m = size( mat, 2 )
!
    if ( (m*p<=0_i4b) .or. (taperp<=zero) .or. (taperp>one) ) return
!
    n    = int( taperp*real( m, stnd ) + half, i4b )/2_i4b
!
    if ( n>0_i4b ) then
!
        con  = one/real( n, stnd )
!
        do i = 1_i4b, n
!
            weight = half - half*cos( pi*(real(i,stnd)-half)*con )
            mat(:p,i) = weight*mat(:p,i)
!
            j = m + 1_i4b - i
!
            mat(:p,j) = weight*mat(:p,j)
!
        end do
!
    end if
!
! END OF SUBROUTINE taper_rm
! __________________________
!
    end subroutine taper_rm
!
!
! =========================================================================================
!                               POWER SPECTRA AND CROSS SPECTRA PROCEDURES
! =========================================================================================
!
!
    function data_window( n, win, taperp ) result( wk )
!
! Purpose
! _______
!
!   Function DATA_WINDOW computes data windows used in spectral computations.
!
!
! Arguments
! _________
!
!   N        (INPUT) integer(i4b)
!            The size of the data window. N must be an even positive integer.
!
!   WIN      (INPUT) integer(i4b)
!            On entry, this argument specify the form of the data window.
!            If:
!
!            -  WIN=+1 The Bartlett window is used
!            -  WIN=+2 The square window is used
!            -  WIN=+3 The Welch window is used
!            -  WIN=+4 The Hann window is used
!            -  WIN=+5 The Hamming window is used
!            -  WIN=+6 A split-cosine-bell window is used
!   
!            For other values of WIN, a square window is returned.   
!   
!   TAPERP   (INPUT, OPTIONAL) real(stnd)
!            The total percentage of the data to be tapered if WIN=6.
!            TAPERP must be greater than zero and less or equal to one,
!            otherwise the default value is used.
!
!            The default is 0.2 .
!
!
! Further Details
! _______________
!
!   For more details, see:
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 5.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : arth, merror
    use Reals_Constants,   only : zero, half, one, pi, c0_2
    use Char_Constants,    only : tseries_error23
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in) :: n, win
!
    real(stnd), intent(in), optional :: taperp
!
    real(stnd), dimension(n) :: wk
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)                 :: c1, c2, taperp2
    real(stnd), dimension(n/2) :: wk2
!
    integer(i4b) :: nd2, i1
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='data_window'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<=0_i4b ) return
!
    nd2 = n/2_i4b
!
    if ( n/=2*nd2 )   &
    call merror( name_proc//tseries_error23 )
!
!   CALCULATE DATA WINDOW.
!
    c1 = real( nd2, stnd )
    c2 = one/c1
!
    select case( win )
!
        case( 1_i4b )
!
!           BARTLETT WINDOW.
!
            wk2(:nd2)  = one - abs( (arth(zero, one, nd2)-c1)*c2 )
!
        case( 2_i4b )
!
!           SQUARE WINDOW.
!
            wk2(:nd2) = one
!
        case( 3_i4b )
!
!           WELCH WINDOW.
!
            wk2(:nd2)   = one - ( (arth(zero, one, nd2)-c1)*c2 )**2
!
        case( 4_i4b )
!
!           HANN WINDOW.
!
            c1 = pi*c2
            wk2(:nd2) = half*(one - cos( arth(zero, c1, nd2) ))
!
        case( 5_i4b )
!
!           HAMMING WINDOW.
!
            c1 = pi*c2
            wk2(:nd2) = 0.54_stnd - 0.46_stnd*cos( arth(zero, c1, nd2) )
!
        case( 6_i4b )
!
!           SPLIT-COSINE-BELL WINDOW.
!
            taperp2 = c0_2
            if ( present( taperp ) ) then
                if ( (taperp2>zero) .and. (taperp2<=one) ) then
                    taperp2 = taperp
                end if
            end if
!
            i1  = int( taperp2*real( n, stnd ) + half, i4b )/2_i4b
!
            if ( i1>0_i4b ) then
!
                c2 = pi/real( i1, stnd )
                c1 = (pi-half)/real( i1, stnd )
!
                wk2(1_i4b:i1)     = half - half*cos( arth(c1, c2, i1) )
                wk2(i1+1_i4b:nd2) = one 
!
            else
                wk2(:nd2) = one
            end if
!
        case default
!
            wk2(:nd2) = one
!
    end select
!
    wk(:nd2)        = wk2(:nd2)
    wk(nd2+1_i4b:n) = wk2(nd2:1_i4b:-1_i4b)
!
!
! END OF FUNCTION data_window
! ___________________________
!
    end function data_window
!
! =========================================================================================
!
    function estim_dof( wk, win, smooth_param, l0, nseg, overlap  ) result( dof )
!
! Purpose
! _______
!
!   Function ESTIM_DOF computes "the equivalent number of degrees of freedom" of power
!   and cross spectrum estimates as calculated by subroutines POWER_SPECTRUM, CROSS_SPECTRUM,
!   POWER_SPECTRUM2 and CROSS_SPECTRUM2.
!
!
! Arguments
! _________
!
!   WK            (INPUT)  real(stnd), dimension(:)
!                 On entry, this argument specify the data window used in the computations of the
!                 power and/or cross spectra.
!
!                 Spectral computations are at (Size(WK)/2)+1 frequencies if the optional
!                 argument L0 is absent and are at ((Size(WK)+L0)/2)+1 frequencies if L0 is
!                 present (L0 is the number of zeros added to each segment).

!                 Size(WK) must be greater or equal to 4 and Size(WK)+L0 must be even.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the form of the data window given in argument WK.
!                 If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!   
!                 For other values of WIN, a message error is issued and the program is stopped.
!
!                 The default is WIN=+3, e.g. the Welch window.
!
!   SMOOTH_PARAM  (INPUT, OPTIONAL) integer(i4b), dimension(:)
!                 if SMOOTH_PARAM is used, the power and/or cross spectrum have been estimated
!                 by repeated smoothing of the periodogram with modified Daniell weights.
!
!                 On entry, SMOOTH_PARAM(:) gives the array of the half-lengths of the
!                 modified Daniell filters that have been applied.
!
!                 All the values in SMOOTH_PARAM(:) must be greater than 0 and less
!                 than ((size(WK)+L0)/2) + 1 .
!
!                 Size(SMOOTH_PARAM) must be greater or equal to 1.
!
!   L0            (INPUT, OPTIONAL) integer(i4b)
!                 The number of zeros added to the time series (or segment) in order to
!                 obtain more finely spaced spectral estimates. L0 must be a positive integer.
!                 Moreover, Size(VEC)+L0 must be even.
!
!                 The default is L0=0, e.g. no zeros are added to the time series.
!
!   NSEG          (INPUT, OPTIONAL) integer(i4b)
!                 The number of segments if the spectra have been computed by 
!                 POWER_SPECTRUM2 and CROSS_SPECTRUM2 . NSEG must be a positive integer.
!
!                 The segments are assumed to be independent or to overlap by one half of
!                 their length if the optional argument OVERLAP is used and is set to true.
!                 Let L = size(WK). Then, the number of segments may be computed as follows:
!
!                 - N/L      if OVERLAP=false
!                 - (2N/L)-1 if OVERLAP=true
!
!                 where N is equal to:
!
!                 - the length of the original time series (call it M) if
!                   this length is evenly divisible by L,
!                 - M+L-mod(M,L) if M is not evenly divisible L.
!
!                 The default is NSEG=1, e.g. the time series is not segmented.
!
!   OVERLAP       (INPUT, OPTIONAL) logical(lgl)
!                 If OVERLAP is set to false, the spectrum estimates have been computed
!                 from nonoverlapping segments.
!
!                 If OVERLAP is set to true, the spectrum estimates have been computed
!                 from overlapped segments (subroutines POWER_SPECTRUM2 and CROSS_SPECTRUM2
!                 may overlap the segments by one half of their length.
!
!                 The default is OVERLAP=false .
!
!
! Further Details
! _______________
!
!   The computed equivalent number of degrees of freedom must be divided by two for the zero
!   and Nyquist frequencies.
!
!   Furthermore, the computed equivalent number of degrees of freedom is not right near the 
!   zero and Nyquist frequencies if the PSD estimates have been smoothed by modified
!   Daniell filters.
!   The reason is that ESTIM_DOF assumes that smoothing involves averaging independent frequency
!   ordinates. This is true except near the zero and Nyquist frequencies where an average
!   may contain contributions from negative frequencies, which are identical to and hence not
!   independent of positive frequency spectral values. Thus, the number of degrees of freedom
!   in PSD estimates near the 0 and Nyquist frequencies are as little as half the number
!   of degrees of freedom of the spectral estimates away from these frequency extremes if
!   the optional argument SMOOTH_PARAM is used.
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 8.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Reals_Constants,   only : two, nine, c11, c16, c18
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error15,   &
                                  tseries_error24, tseries_error50, tseries_error57,  &
                                  tseries_error58, tseries_error63, tseries_error76
    use Logical_Constants, only : true, false
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),    intent(in), dimension(:)           :: wk
!
    integer(i4b),  intent(in),               optional :: win, l0, nseg
    integer(i4b),  intent(in), dimension(:), optional :: smooth_param
!
    logical(lgl),  intent(in), optional ::  overlap
!
    real(stnd)          :: dof
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd) :: g2
    real(stnd), dimension(:), allocatable :: coef
!
    integer(i4b) :: m, sl, nparam, ncoef, nseg2, win2
    integer      :: iok
!
    logical(lgl) :: smooth, overlap2
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='estim_dof'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m = size( wk )
!
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
!   DETERMINE THE LENGTH OF THE SEGMENTS.
!
    if ( present(l0) ) then
!
        if ( l0<0_i4b )     &
        call merror( name_proc//tseries_error24 )
!
        sl = m + l0
!
    else
!
        sl = m
!
    end if
!
!   CHECK IF sl IS AN EVEN INTEGER.
!
    if ( sl/=2_i4b*(sl/2_i4b) )     &
    call merror( name_proc//tseries_error63 )
!
!   DETERMINE THE FORM OF THE DATA WINDOW.
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
!
    else
!
        win2 = 3_i4b
!
    end if
!
    smooth = false
!
    if ( present(smooth_param) ) then
!
!       THE SPECTRAL ESTIMATES HAVE BEEN SMOOTHED WITH DANIELL FILTERS.
!
        nparam = size( smooth_param )
!
        if ( nparam<=0_i4b  )    &
        call merror( name_proc//tseries_error50 )
!
!       CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIELL FILTERS.
!
        if ( any( smooth_param(:nparam)<=0_i4b .or. smooth_param(:nparam)>=(sl/2_i4b)+1_i4b ) )   &
        call merror( name_proc//tseries_error76 )
!
        smooth = true
!
    end if
!
!   DETERMINE THE NUMBER OF THE SEGMENTS.
!
    overlap2 = false
!
    if ( present(nseg) ) then
!
        if ( nseg<=0_i4b )       &
        call merror( name_proc//tseries_error57 )
!
        nseg2 = nseg
!
        if ( present(overlap) ) then
            overlap2 = overlap
        end if
!
        if ( nseg2/=1_i4b .and. overlap2 .and. win2/=1_i4b .and. win2/=3_i4b  )    &
        call merror( name_proc//tseries_error58 )
!
    else
!
        nseg2 = 1_i4b
!
    end if
!
!   COMPUTE THE EQUIVALENT NUMBER OF DEGREES OF FREEDOM.
!
    dof = two
!
    if ( nseg2==1_i4b ) then
!
!        SCALE THE NUMBER OF DEGREES OF FREEDOM TO ADJUST FOR
!        FRACTIONAL PIECES OF INFORMATION THAT OCCUR AT EACH 
!        FOURIER SERIES FREQUENCY BECAUSE OF LEAKAGE WHEN THE 
!        TIME SERIES IS PADDED WITH ZEROS OR/AND TAPERED BEFORE
!        FOURIER TRANSFORMING.
!
        g2 = real( sl, stnd)*(sum(wk(:m)**4)/(sum(wk(:m)**2))**2)
!
        dof = dof/g2
!
    else
!
!       WELCH METHOD OF AVERAGING PERIODOGRAMS.
!
        if ( overlap2 ) then
!
            select case( win2 )
!
                case( 1_i4b )
!
                    dof = (dof*real(nseg2, stnd )*c16)/c18
!
                case( 3_i4b )
!
                    dof = (dof*real(nseg2, stnd )*nine)/c11
!
            end select
!
        else
!
            dof = dof*real(nseg2, stnd )
!
        end if
!
!        dof = dof*( real( m, stnd)/real( sl, stnd) )
!
    end if
!
    if ( smooth ) then
!
!       THE SPECTRUM IS ESTIMATED BY SMOOTHING THE PERIODOGRAM.
!
!       COMPUTE THE NUMBER OF FILTER COEFFICIENTS.
!
        ncoef = 2_i4b*( 2_i4b + sum( smooth_param(:nparam) )) - 1_i4b
!
!       ALLOCATE WORK ARRAY FOR THE FILTER COEFFICIENTS.
!
        allocate( coef(ncoef),  stat=iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       FUNCTION moddan_coef COMPUTES THE IMPULSE RESPONSE FUNCTION CORRESPONDING TO A NUMBER
!       OF APPLICATIONS MODIFIED DANIELL FILTERS .
!
        coef(:ncoef) = moddan_coef( ncoef, smooth_param(:nparam) )
!
        dof = dof/dot_product( coef(:ncoef), coef(:ncoef) )
!
!       DEALLOCATE WORK ARRAY FOR THE FILTER COEFFICIENTS.
!
        deallocate( coef )
!
    end if
!
!
! END OF FUNCTION estim_dof
! _________________________
!
    end function estim_dof
!
! =========================================================================================
!
    function estim_dof2( wk, l0, win, nsmooth, nseg, overlap  ) result( dof )
!
! Purpose
! _______
!
!   Function ESTIM_DOF2 computes "the equivalent number of degrees of freedom" of power
!   and cross spectrum estimates as calculated by subroutines POWER_SPCTRM, CROSS_SPCTRM,
!   POWER_SPCTRM2 and CROSS_SPCTRM2.
!
!
! Arguments
! _________
!
!   WK            (INPUT)  real(stnd), dimension(:)
!                 On entry, this argument specifies the data window used in the computations of the
!                 power and/or cross spectra.
!
!                 Spectral computations are at ((Size(WK)+L0)/2)+1 frequencies (L0 is the number of
!                 zeros added to each segment).
!
!                 Size(WK) must be greater or equal to 4 and Size(WK)+L0 must be even.
!
!   L0            (INPUT) integer(i4b)
!                 The number of zeros added to the time series (or segment) in order to
!                 obtain more finely spaced spectral estimates. L0 must be a positive integer.
!                 Moreover, Size(VEC)+L0 must be even.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the form of the data window given in argument WK.
!                 If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!   
!                 For other values of WIN, a message error is issued and the program is stopped.
!
!                 The default is WIN=+3, e.g. the Welch window.
!
!   NSMOOTH       (INPUT, OPTIONAL) integer(i4b)
!                 if NSMOOTH is used, the power and/or cross spectra have been estimated
!                 by smoothing the periodogram with Daniell weights.
!
!                 On entry, NSMOOTH gives the length of the Daniell filter
!                 that has been applied.
!
!                 Setting NSMOOTH=0 on entry is equivalent to omit the optional argument
!                 NSMOOTH. Otherwise, NSMOOTH must be odd, greater than 2 and less or
!                 equal to ((size(WK)+L0)/2) + 1 .
!
!   NSEG          (INPUT, OPTIONAL) integer(i4b)
!                 The number of segments if the spectra have been computed by 
!                 POWER_SPCTRM2 and CROSS_SPCTRM2 . NSEG must be a positive integer.
!
!                 The segments are assumed to be independent or to overlap by one half of
!                 their length if the optional argument OVERLAP is used and is set to true.
!                 Let L = size(WK). Then, the number of segments may be computed as follows:
!
!                 - N/L      if OVERLAP=false
!                 - (2N/L)-1 if OVERLAP=true
!
!                 where N is equal to:
!
!                 - the length of the original time series (call it M) if
!                   this length is evenly divisible by L,
!                 - M+L-mod(M,L) if M is not evenly divisible L.
!
!                 The default is NSEG=1, e.g. the time series is not segmented.
!
!   OVERLAP       (INPUT, OPTIONAL) logical(lgl)
!                 If OVERLAP is set to false, the spectrum estimates have been computed
!                 from nonoverlapping segments.
!
!                 If OVERLAP is set to true, the spectrum estimates have been computed
!                 from overlapped segments (subroutines POWER_SPCTRUM2 and CROSS_SPCTRUM2
!                 may overlap the segments by one half of their length.
!
!                 The default is OVERLAP=false .
!
!
! Further Details
! _______________
!
!   For more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!           Fourier analysis of time series- An introduction. 
!           John Wiley and Sons, New York, Chapter 8.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Reals_Constants,   only : one, two, nine, c11, c16, c18
    use Char_Constants,    only : tseries_error10, tseries_error15,  &
                                  tseries_error24, tseries_error57,  &
                                  tseries_error58, tseries_error63,  &
                                  tseries_error61, tseries_error64
    use Logical_Constants, only : true, false
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),    intent(in), dimension(:) :: wk
!
    integer(i4b),  intent(in)            :: l0
    integer(i4b),  intent(in),  optional :: win, nseg, nsmooth
!
    logical(lgl),  intent(in), optional ::  overlap
!
    real(stnd), dimension( (size(wk)+l0)/2 + 1) :: dof
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd) :: g2, dof2
!
    integer(i4b) :: m, sl, nf, nseg2, win2, i, j, k, df, lim
!
    logical(lgl) :: smooth, overlap2
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='estim_dof2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m = size( wk )
!
    if ( m<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
!   DETERMINE THE LENGTH OF THE SEGMENTS.
!
    if ( l0<0_i4b )     &
    call merror( name_proc//tseries_error24 )
!
    sl = m + l0
!
!   CHECK IF sl IS AN EVEN INTEGER.
!
    if ( sl/=2_i4b*(sl/2_i4b) )     &
    call merror( name_proc//tseries_error63 )
!
    nf = ( sl/2_i4b )+1_i4b
!
!   DETERMINE THE FORM OF THE DATA WINDOW.
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
!
    else
!
        win2 = 3_i4b
!
    end if
!
    smooth = false
!
    if ( present(nsmooth) ) then
!
        if ( nsmooth/=0_i4b ) then
!
!           THE SPECTRAL ESTIMATES HAVE BEEN SMOOTHED WITH A DANIELL FILTER.
!           CHECK THE INPUT VALUE FOR THE LENGTH OF THE DANIELL FILTER.
!
            if ( (nsmooth/2_i4b)*2_i4b==nsmooth )   &
            call merror( name_proc//tseries_error61 )
!
            if ( nsmooth<2_i4b .or. nsmooth>nf  )   &
            call merror( name_proc//tseries_error64 )
!
            smooth = true
!
        end if
!
    end if
!
!   DETERMINE THE NUMBER OF THE SEGMENTS.
!
    overlap2 = false
!
    if ( present(nseg) ) then
!
        if ( nseg<=0_i4b )       &
        call merror( name_proc//tseries_error57 )
!
        nseg2 = nseg
!
        if ( present(overlap) ) then
            overlap2 = overlap
        end if
!
        if ( nseg2/=1_i4b .and. overlap2 .and. win2/=1_i4b .and. win2/=3_i4b  )    &
        call merror( name_proc//tseries_error58 )
!
    else
!
        nseg2 = 1_i4b
!
    end if
!
!   COMPUTE THE EQUIVALENT NUMBER OF DEGREES OF FREEDOM.
!
    if ( nseg2==1_i4b ) then
!
!        SCALE THE NUMBER OF DEGREES OF FREEDOM TO ADJUST FOR
!        FRACTIONAL PIECES OF INFORMATION THAT OCCUR AT EACH 
!        FOURIER SERIES FREQUENCY BECAUSE OF LEAKAGE WHEN THE 
!        TIME SERIES IS PADDED WITH ZEROS OR/AND TAPERED BEFORE
!        FOURIER TRANSFORMING.
!
        g2 = real( sl, stnd)*(sum(wk(:m)**4)/(sum(wk(:m)**2))**2)
!
        dof2 = one/g2
!
    else
!
!       WELCH METHOD OF AVERAGING PERIODOGRAMS.
!
        if ( overlap2 ) then
!
            select case( win2 )
!
                case( 1_i4b )
!
                    dof2 = (real(nseg2, stnd )*c16)/c18
!
                case( 3_i4b )
!
                    dof2 = (real(nseg2, stnd )*nine)/c11
!
            end select
!
        else
!
            dof2 = real(nseg2, stnd )
!
        end if
!
!        dof2 = dof2*( real( m, stnd)/real( sl, stnd) )
!
    end if
!
    if ( smooth ) then
!
!       THE SPECTRUM IS ESTIMATED BY SMOOTHING THE PERIODOGRAM.

        lim = ( nsmooth - 1_i4b )/2_i4b
!
        do i = 1_i4b, nf
!
!           i IS THE INDEX OF CENTER OF INTERVAL TO BE SMOOTHED
!           (E.G. THE I-TH SPECTRAL ESTIMATES).
!
            df = 0_i4b
!
            do j = -lim, lim
!
                k = i + j
!
!               COUNT THE NUMBER OF DEGREES OF FREEDOM ADDED BY THE
!               PERIODOGRAM VALUE OF INDEX k .
!
                if ( k==1_i4b .or. k==nf ) then
!
!                   k IS THE ZERO OR NYQUIST FREQUENCIES.
!
                    df = df + 1_i4b
!
                else if ( k>1_i4b .and. k<nf ) then
!
!                   k LIES BETWEEN THE ZERO OR NYQUIST FREQUENCIES.
!
                    df = df + 2_i4b
!
                end if
!
            end do
!
            dof(i) = dof2*real( df, stnd )
!
        end do
!
    else
!
!       THE SPECTRUM HAS NOT BEEN SMOOTHED WITH A DANIELL FILTER.
!       HOWEVER, THE COMPUTED EQUIVALENT NUMBER OF DEGREES OF DEGREES
!       OF FREEDOM MUST BE DIVIDED BY TWO FOR THE ZERO AND NYQUIST FREQUENCIES.
!
        dof(1_i4b)          = dof2        
        dof(2_i4b:nf-1_i4b) = dof2*two   
        dof(nf)             = dof2       
!
    end if
!
!
! END OF FUNCTION estim_dof2
! __________________________
!
    end function estim_dof2
!
! =========================================================================================
!
    subroutine comp_conflim_r( edof, probtest, conlwr, conupr, testcoher )
!
! Purpose
! _______
!
!   Subroutine COMP_CONFLIM estimates confidence limit factors for spectral estimates
!   and, optionally, critical value for testing the null hypothesis that squared coherency
!   is zero.
!
!
! Arguments
! _________
!
!   EDOF          (INPUT) real(stnd)
!                 On entry, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ) and the critical value for
!                 testing the null hypothesis that the squared coherency is zero (e.g.
!                 the TESTCOHER optional argument).
!                 PROBTEST must verify   0. < P < 1.
!
!                 The default is 0.05 .
!
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd)
!
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates by these constants to get the
!                 lower and upper limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!   
!   TESTCOHER     (OUTPUT, OPTIONAL) real(stnd)
!                 On output, this argument specifies the critical value for testing the null
!                 hypothesis that the squared coherency is zero at the PROBTEST * 100% significance
!                 level (e.g. squared coherencies less than TESTCOHER should be regarded as not
!                 significantly different from zero at the PROBTEST * 100% significance level).
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Reals_Constants,   only : zero, half, one, two, c5_m2
    use Num_Constants,     only : nan
    use Char_Constants,    only : tseries_error59
    use Prob_Procedures,   only : pinvf2, pinvq2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in) :: edof
!
    real(stnd), intent(in),  optional :: probtest
    real(stnd), intent(out), optional :: conlwr, conupr, testcoher
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: c1, c2, miss, probtest2, con
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='comp_conflim'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    probtest2 = c5_m2
!
    if ( present(probtest) ) then
!
        if ( zero>=probtest .or. probtest>=one ) then
            call merror( name_proc//tseries_error59 )
        else
            probtest2 = probtest
        end if
!
    end if
!
!   OUTPUT CONFIDENCE LIMIT FACTORS AND
!   PROBABILITIES FOR SQUARED COHERENCIES IF REQUIRED.
!
    c1 = probtest2*half
    c2 = one - c1
!
    miss = nan()
!
    if ( present(conlwr) ) then
!
        if ( edof>=half ) then
            conlwr = edof/pinvq2( c2, edof )
        else
            conlwr = miss
        end if
!
    end if
!
    if ( present(conupr) ) then
!
        if ( edof>=half ) then
            conupr = edof/pinvq2( c1, edof )
        else
            conupr = miss
        end if
!
    end if
!
    if ( present(testcoher) ) then
!
        if ( edof>two ) then
!
!            c1 = two/(edof-two)
!            c2 = one - probtest2**c1
!            testcoher = max( min( one, c2 ), zero )
!
            con = one - probtest2
            c2  = edof - two
            c1  = two*pinvf2( con, two, c2 )
            testcoher = c1/( c2 + c1 )
        else
            testcoher = miss
        end if
!
    end if
!
!
! END OF SUBROUTINE comp_conflim_r
! ________________________________
!
    end subroutine comp_conflim_r
!
! =========================================================================================
!
    subroutine comp_conflim_rv( edof, probtest, conlwr, conupr, testcoher )
!
! Purpose
! _______
!
!   Subroutine COMP_CONFLIM estimates confidence limit factors for spectral estimates
!   and, optionally, critical values for testing the null hypothesis that squared coherencies
!   are zero.
!
!
! Arguments
! _________
!
!   EDOF          (INPUT) real(stnd), dimension(:)
!                 On entry, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ) and the critical value for
!                 testing the null hypothesis that the squared coherency is zero (e.g.
!                 the TESTCOHER optional argument).
!                 PROBTEST must verify   0. < P < 1.
!
!                 The default is 0.05 .
!
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates by these constants to get the
!                 lower and upper limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!                 CONLWR must verify:  size(CONLWR) = size(EDOF) .
!
!                 CONUPR must verify:  size(CONUPR) = size(EDOF) .
!   
!   TESTCOHER     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, this argument specifies the critical values for testing the null
!                 hypothesis that the squared coherencies are zero at the PROBTEST * 100% significance
!                 level (e.g. squared coherencies less than TESTCOHER(:) should be regarded as not
!                 significantly different from zero at the PROBTEST * 100% significance level).
!
!                 TESTCOHER must verify:  size(TESTCOHER) = size(EDOF) .
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror
    use Reals_Constants,   only : zero, half, one, two, c5_m2
    use Num_Constants,     only : nan
    use Char_Constants,    only : tseries_error59
    use Prob_Procedures,   only : pinvf2, pinvq2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),  dimension(:) :: edof
!
    real(stnd), intent(in),                optional :: probtest
    real(stnd), intent(out), dimension(:), optional :: conlwr, conupr, testcoher
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: nf, i
!
    real(stnd)   :: c1, c2, miss, chi2, probtest2, df, df2, con
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='comp_conflim'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    nf = size( edof )
!
    if ( nf<1_i4b ) return
!
    probtest2 = c5_m2
!
    if ( present(probtest) ) then
!
        if ( zero>=probtest .or. probtest>=one ) then
            call merror( name_proc//tseries_error59 )
        else
            probtest2 = probtest
        end if
!
    end if
!
    if ( present( conlwr ) ) then
        call assert( logical(int(size(conlwr),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    if ( present( conupr ) ) then
        call assert( logical(int(size(conupr),i4b)==nf,lgl),    &
        name_proc )
    end if
!
    if ( present( testcoher ) ) then
        call assert( logical(int(size(testcoher),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
!   OUTPUT CONFIDENCE LIMIT FACTORS AND
!   PROBABILITIES FOR SQUARED COHERENCIES IF REQUIRED.
!
    c1 = probtest2*half
    c2 = one - c1
!
    miss = nan()
!
    if ( present(conlwr) ) then
!
        df2 = -one
!
        do i = 1_i4b, nf 
!
               df = edof(i)
!
               if ( df>=half ) then
!
                   if ( df/=df2 ) then
                       df2  = df
                       chi2 = pinvq2( c2, df )
                       con  = df/chi2
                   end if
!
                   conlwr(i) = con
!
               else
!
                   conlwr(i) = miss
!
               end if
!
       end do
!
    end if
!
    if ( present(conupr) ) then
!
        df2 = -one
!
        do i = 1_i4b, nf 
!
                df = edof(i)
!
                if ( df>=half ) then
!
                    if ( df/=df2 ) then
                        df2  = df
                        chi2 = pinvq2( c1, df )
                        con  = df/chi2
                    end if
!
                    conupr(i) = con
!
                else
!
                    conupr(i) = miss
!
                end if
!
        end do
!
    end if
!
    if ( present(testcoher) ) then
!
        df2 = -one
        con = one - probtest2
!
        do i = 1_i4b, nf 
!
                df = edof(i)
!
                if ( df>two ) then
!
                    if ( df/=df2 ) then
!
                        df2 = df
!
!                        c1 = two/(df-two)
!                        c2 = max( min( one, one - probtest2**c1 ), zero )
!
                        c2 = df - two
                        c1  = two*pinvf2( con, two, c2 )
                        c2 = c1/( c2 + c1 )
!
                    end if

                    testcoher(i) = c2
!
                else
!
                    testcoher(i) = miss
!
                end if
!
        end do
!
    end if
!
!
! END OF SUBROUTINE comp_conflim_rv
! _________________________________
!
    end subroutine comp_conflim_rv
!
! =========================================================================================
!
    subroutine spctrm_ratio_r( edofn, edofd, lwr_ratio, upr_ratio, pinterval )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_RATIO calculates a pointwise tolerance interval for the ratio of two 
!   estimated spectra under the assumption that the two "true" underlying spectra
!   are the same.
!
!
! Arguments
! _________
!
!   EDOFN         (INPUT) real(stnd)
!                 On exit, the equivalent number of degrees of freedom of the first
!                 estimated spectrum (e.g. the numerator of the ratio of the
!                 two estimated spectra).
!
!                 EDOFN must be greater than zero.
!
!   EDOFD        (INPUT) real(stnd)
!                 On exit, the equivalent number of degrees of freedom of the second
!                 estimated spectrum (e.g. the denominator of the ratio of the
!                 two estimated spectra).
!
!                 EDOFD must be greater than zero.
!
!
!   LWR_RATIO     (OUTPUT) real(stnd)
!
!   UPR_RATIO     (OUTPUT) real(stnd)
!                 On output, these arguments specify the lower and upper critical ratios of the
!                 computed PINTERVAL * 100% tolerance interval for the ratio of the power spectral
!                 density estimates.
!
!          
!   PINTERVAL     (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. This probability is used to determine the upper and lower
!                 critical ratios of the computed tolerance interval. A PINTERVAL * 100% tolerance interval
!                 is computed and output in the two arguments LWR_RATIO and  UPR_RATIO.
!                 PINTERVAL must verify:   0. < PINTERVAL < 1.
!
!                 The default value is 0.90, e.g. a 90% tolerance interval is computed.
!
!
! Further Details
! _______________
!
!   For more details, see:
!
!   (1) Diggle, P.J., 1990:
!         Time series: a biostatistical introduction.
!         Clarendon Press, Oxford, Chapter 4.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : merror
    use Reals_Constants,   only : zero, half, one, c5_m2
    use Char_Constants,    only : tseries_error68, tseries_error69
    use Prob_Procedures,   only : pinvf2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in) :: edofn, edofd
!
    real(stnd), intent(in),  optional :: pinterval
    real(stnd), intent(out)           :: lwr_ratio, upr_ratio
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: c1, c2
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_ratio'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    c1 = c5_m2
!
    if ( present(pinterval) ) then
!
        if ( zero>=pinterval .or. pinterval>=one ) then
            call merror( name_proc//tseries_error68 )
        else
            c1 = half*( one - pinterval )
        end if
!
    end if
!
    if ( zero>=edofn .or. zero>=edofd )    &
    call merror( name_proc//tseries_error69 )
!
!   COMPUTE THE UPPER AND LOWER RATIOS OF THE pinterval*100% 
!   TOLERANCE INTERVAL.
!
    c2   = one - c1
!
    lwr_ratio = pinvf2( c1, edofn, edofd )
    upr_ratio = pinvf2( c2, edofn, edofd )
!
!
! END OF SUBROUTINE spctrm_ratio_r
! ________________________________
!
    end subroutine spctrm_ratio_r
!
! =========================================================================================
!
    subroutine spctrm_ratio_rv( edofn, edofd, lwr_ratio, upr_ratio, pinterval )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_RATIO calculates pointwise tolerance intervals for the ratio of two 
!   estimated spectra under the assumption that the two "true" underlying spectra are the same.
!
!
! Arguments
! _________
!
!   EDOFN         (INPUT) real(stnd), dimension(:)
!                 On exit, the equivalent number of degrees of freedom of the first
!                 estimated spectrum (e.g. the numerator of the ratio of the
!                 two estimated spectra).
!
!                 Elements of EDOFN(:) must be greater than zero.
!
!   EDOFD        (INPUT) real(stnd), dimension(:)
!                 On exit, the equivalent number of degrees of freedom of the second
!                 estimated spectrum (e.g. the denominator of the ratio of the
!                 two estimated spectra).
!                 Elements of EDOFD(:) must be greater than zero.
!
!                 EDOFD must verify:  size(EDOFD) = size(EDOFN) .
!
!   LWR_RATIO     (OUTPUT) real(stnd), dimension(:)
!
!   UPR_RATIO     (OUTPUT) real(stnd), dimension(:)
!                 On output, these arguments specify the lower and upper critical ratios of the
!                 computed PINTERVAL * 100% tolerance interval for the ratio of the  power spectral
!                 density estimates.
!
!                 LWR_RATIO must verify:  size(LWR_RATIO) = size(EDOFN) .
!
!                 UPR_RATIO must verify:  size(UPR_RATIO) = size(EDOFN) .
!          
!   PINTERVAL     (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. This probability is used to determine the upper and lower
!                 critical ratios of the computed tolerance interval. A PINTERVAL * 100% tolerance interval
!                 is computed and output in the two arguments LWR_RATIO and  UPR_RATIO.
!                 PINTERVAL must verify:   0. < PINTERVAL < 1.
!
!                 The default value is 0.90, e.g. a 90% tolerance interval is computed .
!
!
! Further Details
! _______________
!
!   For more details, see:
!
!   (1) Diggle, P.J., 1990:
!         Time series: a biostatistical introduction.
!         Clarendon Press, Oxford, Chapter 4.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert_eq, merror
    use Reals_Constants,   only : zero, half, one, c5_m2
    use Char_Constants,    only : tseries_error68, tseries_error69
    use Prob_Procedures,   only : pinvf2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),  dimension(:) :: edofn, edofd
!
    real(stnd), intent(in),                optional :: pinterval
    real(stnd), intent(out), dimension(:)           :: lwr_ratio, upr_ratio
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: nf, i
!
    real(stnd)   :: c1, c2, con1, con2, dfn, df2n, dfd, df2d
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_ratio'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    nf =  assert_eq( int(size(edofn),i4b) ,          &
                     int(size(edofd),i4b) ,          &
                     int(size(lwr_ratio),i4b) ,      &
                     int(size(upr_ratio),i4b) ,      &
                     name_proc )
!
    if ( nf<1_i4b ) return
!
    c1 = c5_m2
!
    if ( present(pinterval) ) then
!
        if ( zero>=pinterval .or. pinterval>=one ) then
            call merror( name_proc//tseries_error68 )
        else
            c1 = half*( one - pinterval )
        end if
!
    end if
!
    if ( any( zero>=edofn(:nf) .or. zero>=edofd(:nf) ) )    &
    call merror( name_proc//tseries_error69 )
!
!   COMPUTE THE UPPER AND LOWER RATIOS OF THE pinterval*100% 
!   TOLERANCE INTERVAL.
!
    c2 = one - c1
!
    df2n = -one
    df2d = -one
!
    do i = 1_i4b, nf 
!
        dfn = edofn(i)
        dfd = edofd(i)
!
        if ( dfn/=df2n .or. dfd/=df2d ) then
!
            df2n  = dfn
            df2d  = dfd
!
            con1 = pinvf2( c1, dfn, dfd )
            con2 = pinvf2( c2, dfn, dfd )
!
        end if
!
        lwr_ratio(i) = con1
        upr_ratio(i) = con2
!
    end do
!
!
! END OF SUBROUTINE spctrm_ratio_rv
! _________________________________
!
    end subroutine spctrm_ratio_rv
!
! =========================================================================================
!
    subroutine spctrm_ratio2_rv( psvecn, psvecd, edofn, edofd, prob, min_ratio, max_ratio, &
                                 prob_min_ratio, prob_max_ratio )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_RATIO2 calculates a conservative critical probability value (e.g. p-value)
!   for testing the hypothesis of a common spectrum for two estimated spectra (e.g. the arguments
!   PSVECN, PSVECD). This conservative critical probability value is computed from the minimim and
!   maximum values of the ratio of the two estimated spectra and the associated probabilities of 
!   obtaining, respectively, a value less (for the minimum ratio) and higher (for the maximum ratio) 
!   than attained under the null hypothesis of a common spectrum for the two time series.
!
!
! Arguments
! _________
!
!
!   PSVECN             (INPUT) real(stnd), dimension(:)
!                      On entry, a real vector containing the Power Spectral Density (PSD)
!                      estimates of the first time series (e.g. the numerator of the ratio of
!                      the two estimated spectra).  
!
!                      All elements in PSVECN(:) must be greater or equal to zero and size(PSVECN)
!                      must be greater or equal to 2.
!
!   PSVECD             (INPUT) real(stnd), dimension(:)
!                      On entry, a real vector containing the Power Spectral Density (PSD)
!                      estimates of the second time series (e.g. the denominator of the ratio of
!                      the two estimated spectra).  
!                      All elements in PSVECD(:) must be greater than zero and size(PSVECD)
!                      must be greater or equal to 2.
!
!                      PSVECD must also verify:  size(PSVECD) = size(PSVECN) .
!
!   EDOFN              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the first
!                      estimated spectrum (e.g. the numerator of the ratio of the
!                      two estimated spectra).
!
!                      EDOFN must be greater than zero.
!
!   EDOFD              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the second
!                      estimated spectrum (e.g. the denominator of the ratio of the
!                      two estimated spectra).
!
!                      EDOFD must be greater than zero.
!
!   PROB               (OUTPUT) real(stnd)
!                      On exit, the conservative critical probability value (e.g. p-value) computed
!                      under the hypothesis that the two "true" underlying spectra are the same.
!                      See the description of the PROB_MIN_RATIO and PROB_MAX_RATIO optional
!                      arguments for more details.
!
!   MIN_RATIO          (OUTPUT, OPTIONAL) real(stnd)
!
!   MAX_RATIO          (OUTPUT, OPTIONAL) real(stnd)
!                      On output, these arguments give, respectively, the minimum and maximum values
!                      of the ratio of the two PSD estimates.
!
!   PROB_MIN_RATIO     (OUTPUT, OPTIONAL) real(stnd)
!
!   PROB_MAX_RATIO     (OUTPUT, OPTIONAL) real(stnd)
!                      On output, these arguments give, respectively, the probabilities of obtaining
!                      a smaller value of the minimum ratio (e.g. the argument MIN_RATIO) and a greater
!                      value of the maximum ratio (e.g. the argument MAX_RATIO) under the null hypothesis
!                      that the two "true" underlying spectra are the same.
!
!                      The PROB argument is computed as 2 * min(PROB_MIN_RATIO,PROB_MAX_RATIO).
!
!
! Further Details
! _______________
!
!   This statistical test relies on the assumptions that the different spectral ordinates
!   have the same sampling distribution and are independent of each other for each series.
!   This means, in particular, that the spectral ordinates corresponding to the zero and
!   Nyquist frequencies must be excluded from the PSVECN and PSVECD vectors before calling
!   SPCTRM_RATIO2 and that the two estimated spectra have not been obtained by smoothing
!   the periodogram in the frequency domain.
!
!   It is also assumed that the PSVECN and PSVECD realizations are independent.
!
!   For more details, see:
!
!   (1) Diggle, P.J., 1990:
!         Time series: a biostatistical introduction.
!         Clarendon Press, Oxford, Chapter 4.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert_eq, merror
    use Reals_Constants,   only : zero, one, two
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : tseries_error69, tseries_error74
    use Prob_Procedures,   only : probf2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),  dimension(:) :: psvecn, psvecd
    real(stnd), intent(in)                :: edofn, edofd
!
    real(stnd), intent(out)             :: prob
    real(stnd), intent(out), optional   :: min_ratio, max_ratio, prob_min_ratio, prob_max_ratio
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: nf, nmask
!
    real(stnd)                            :: minr, maxr, pminr, pmaxr, miss
    real(stnd),  dimension(size(psvecn))  :: ratio
!
    logical(lgl),  dimension(size(psvecn)) :: mask
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_ratio2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    nf =  assert_eq( int(size(psvecn),i4b) ,          &
                     int(size(psvecd),i4b) ,          &
                     name_proc )
!
    if ( nf<=1_i4b )     &
    call merror( name_proc//tseries_error74  )
!
    if ( zero>=edofn .or. zero>=edofd )           &
    call merror( name_proc//tseries_error69 )
!
    mask(:nf) = psvecd(:nf)>zero .and. psvecn(:nf)>=zero
!
    nmask = count( mask(:nf) )
!
    if ( nmask>=2_i4b ) then
!
!       COMPUTE THE RATIO OF THE TWO SPECTRAL ESTIMATES.
!
        where( mask(:nf) )
            ratio(:nf) = psvecn(:nf)/psvecd(:nf)
        elsewhere
            ratio(:nf) = one
        end where
!
!       COMPUTE THE MIN AND MAX OF THE RATIO OF THE TWO SPECTRAL ESTIMATES.
!
        minr = minval( ratio(:nf), mask=mask(:nf) )
        maxr = maxval( ratio(:nf), mask=mask(:nf) )
!
!       COMPUTE THE ATTAINED SIGNIFICANCE LEVELS FOR THE STATISTICS minr AND maxr .
!
        pminr = one - ( probf2( minr, edofn, edofd, true  ) )**nmask
        pmaxr = one - ( probf2( maxr, edofn, edofd, false ) )**nmask
!
!       COMPUTE THE CONSERVATIVE PROBABILITY.
!
        prob = two*min( pminr, pmaxr )
!
    else
!
        miss = nan()
!
        minr  = miss
        maxr  = miss
        pminr = miss
        pmaxr = miss
        prob  = miss
!       
    end if
!
!   OUTPUT THE REQUIRED OPTIONAL ARGUMENTS.
!
    if ( present( min_ratio ) ) then
        min_ratio = minr
    end if
!
    if ( present( max_ratio ) ) then
        max_ratio = maxr
    end if
!
    if ( present( prob_min_ratio ) ) then
        prob_min_ratio = pminr
    end if
!
    if ( present( prob_max_ratio ) ) then
        prob_max_ratio = pmaxr
    end if
!
!
! END OF SUBROUTINE spctrm_ratio2_rv
! __________________________________
!
    end subroutine spctrm_ratio2_rv
!
! =========================================================================================
!
    subroutine spctrm_ratio2_rm( psmatn, psmatd, edofn, edofd, prob, min_ratio, max_ratio, &
                                 prob_min_ratio, prob_max_ratio )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_RATIO2 calculates conservative critical probability values (e.g. p-values)
!   for testing the hypothesis of a common spectrum for the elements of two estimated multi-channel
!   spectra (e.g. the arguments PSMATN, PSMATD).
!
!   These conservative critical probability values are computed from the minimim and maximum 
!   values of the ratio of the two estimated multi-channel spectra and the associated probabilities
!   of obtaining, respectively, a value less (for the minimum ratio) and higher (for the maximum ratio)
!   than attained under the null hypothesis of a common spectrum for the two multi-channel time series.
!
!
! Arguments
! _________
!
!
!   PSMATN             (INPUT) real(stnd), dimension(:,:)
!                      On entry, a real matrix containing the Power Spectral Density (PSD)
!                      estimates of the first multi-channel time series (e.g. the numerator of
!                      the ratio of the two estimated multi-channel spectra).  
!                      Each row of the real matrix PSMATN contains the estimated spectrum
!                      of the corresponding "row" of the first multi-channel times series.  
!
!                      All elements in PSMATN(:,:) must be greater or equal to zero and size(PSMATN,2)
!                      must be greater or equal to 2.
!
!   PSMATD             (INPUT) real(stnd), dimension(:,:)
!                      On entry, a real matrix containing the Power Spectral Density (PSD)
!                      estimates of the second multi-channel time series (e.g. the denominator
!                      of the ratio of the two estimated multi-channel spectra).  
!                      Each row of the real matrix PSMATD contains the estimated spectrum
!                      of the corresponding "row" of the second multi-channel times series.  
!                      All elements in PSMATD(:,:) must be greater than zero and size(PSMATD,2)
!                      must be greater or equal to 2.
!
!                      PSMATD must also verify:
!
!                      - size(PSMATD,1) = size(PSMATN,1) ,
!                      - size(PSMATD,2) = size(PSMATN,2) .
!
!   EDOFN              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the first
!                      estimated spectrum (e.g. the numerator of the ratio of the
!                      two estimated spectra).
!
!                      EDOFN must be greater than zero.
!
!   EDOFD              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the second
!                      estimated spectrum (e.g. the denominator of the ratio of the
!                      two estimated spectra).
!
!                      EDOFD must be greater than zero.
!
!   PROB               (OUTPUT) real(stnd), dimension(:)
!                      On exit, the conservative critical  probability values (e.g. p-values) computed
!                      under the hypothesis that the two "true" underlying multi-channel spectra are 
!                      the same.
!                      See the description of the PROB_MIN_RATIO and PROB_MAX_RATIO optional
!                      arguments for more details.
!
!                      PROB must verify:  size(PROB) = size(PSMATN,1) .
!
!   MIN_RATIO          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!
!   MAX_RATIO          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                      On output, these arguments give, respectively, the minimum and maximum values
!                      of the ratio of the two multi-channel PSD estimates.
!
!                      MIN_RATIO must verify:  size(MIN_RATIO) = size(PSMATN,1) .
!
!                      MAX_RATIO must verify:  size(MAX_RATIO) = size(PSMATN,1) .
!
!   PROB_MIN_RATIO     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!
!   PROB_MAX_RATIO     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                      On output, these arguments give, respectively, the probabilities of obtaining
!                      a smaller value of the minimum ratio (e.g. the argument MIN_RATIO) and a greater
!                      value of the maximum ratio (e.g. the argument MAX_RATIO) under the null hypothesis
!                      that the two "true" underlying multi-channel spectra are the same.
!                      The PROB(:) argument is calculated as 2 * min(PROB_MIN_RATIO(:),PROB_MAX_RATIO(:)).
!
!                      PROB_MIN_RATIO must verify:  size(PROB_MIN_RATIO) = size(PSMATN,1) .
!
!                      PROB_MAX_RATIO must verify:  size(PROB_MAX_RATIO) = size(PSMATN,1) .
!
!
! Further Details
! _______________
!
!   This statistical test relies on the assumptions that the different spectral ordinates
!   have the same sampling distribution and are independent of each other for each series.
!   This means, in particular, that the spectral ordinates corresponding to the zero and
!   Nyquist frequencies must be excluded from the PSMATN and PSMATD matrices before calling
!   SPCTRM_RATIO2 and that the two estimated multi-channel spectra have not been obtained
!   by smoothing the periodogram in the frequency domain.
!
!   It is also assumed that the PSMATN and PSMATD realizations are independent.
!
!   For more details, see:
!
!   (1) Diggle, P.J., 1990:
!         Time series: a biostatistical introduction.
!         Clarendon Press, Oxford, Chapter 4.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert_eq, assert, merror
    use Reals_Constants,   only : zero, one, two
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : tseries_error69, tseries_error75
    use Prob_Procedures,   only : probf2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),  dimension(:,:) :: psmatn, psmatd
    real(stnd), intent(in)                  :: edofn, edofd
!
    real(stnd), intent(out),  dimension(:)             :: prob
    real(stnd), intent(out),  dimension(:), optional   :: min_ratio, max_ratio, prob_min_ratio, prob_max_ratio
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b)                             :: nf, n
    integer(i4b),  dimension(size(psmatn,1)) :: nmask
!
    real(stnd)                                             :: miss
    real(stnd),  dimension(size(psmatn,1))                 :: minr, maxr, pminr, pmaxr
    real(stnd),  dimension(size(psmatn,1),size(psmatn,2))  :: ratio
!
    logical(lgl),  dimension(size(psmatn,1),size(psmatn,2)) :: mask
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_ratio2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n =  assert_eq( int(size(psmatn,1),i4b) ,          &
                    int(size(psmatd,1),i4b) ,          &
                    int(size(prob),i4b) ,              &
                    name_proc )
!
    if ( n<1_i4b ) return
!
    nf =  assert_eq( int(size(psmatn,2),i4b) ,          &
                     int(size(psmatd,2),i4b) ,          &
                     name_proc )
!
    if ( nf<=1_i4b )     &
    call merror( name_proc//tseries_error75 )
!
    if ( present( min_ratio ) ) then
        call assert( logical(int(size(min_ratio),i4b)==n,lgl),    &
                     name_proc )
    end if
!
    if ( present( max_ratio ) ) then
        call assert( logical(int(size(max_ratio),i4b)==n,lgl),    &
                     name_proc )
    end if
!
    if ( present( prob_min_ratio ) ) then
        call assert( logical(int(size(prob_min_ratio),i4b)==n,lgl),    &
                     name_proc )
    end if
!
    if ( present( prob_max_ratio ) ) then
        call assert( logical(int(size(prob_max_ratio),i4b)==n,lgl),    &
                     name_proc )
    end if
!
    if ( zero>=edofn .or. zero>=edofd )   &
    call merror( name_proc//tseries_error69  )
!
    mask(:n,:nf) = psmatd(:n,:nf)>zero .and. psmatn(:n,:nf)>=zero
    nmask(:n)    = count( mask(:n,:nf), dim=2 )
!
!   COMPUTE THE RATIO OF THE TWO SPECTRAL ESTIMATES.
!
    where( mask(:n,:nf) )
         ratio(:n,:nf) = psmatn(:n,:nf)/psmatd(:n,:nf)
    elsewhere
         ratio(:n,:nf) = one
    end where
!
!   COMPUTE THE MIN AND MAX OF THE RATIO OF THE TWO SPECTRAL ESTIMATES.
!
    minr(:n) = minval( ratio(:n,:nf), mask=mask(:n,:nf), dim=2 )
    maxr(:n) = maxval( ratio(:n,:nf), mask=mask(:n,:nf), dim=2 )
!
    miss = nan()
!
    pminr(:n) = probf2( minr(:n), edofn, edofd, true  )
    pmaxr(:n) = probf2( maxr(:n), edofn, edofd, false )
!
    where( nmask(:n)>=2_i4b )
!
!       COMPUTE THE ATTAINED SIGNIFICANCE LEVELS FOR THE STATISTICS minr AND maxr .
!
        pminr(:n) = one - pminr(:n)**nmask(:n)
        pmaxr(:n) = one - pmaxr(:n)**nmask(:n)
!
!       COMPUTE THE CONSERVATIVE PROBABILITY.
!
        prob(:n) = two*min( pminr(:n), pmaxr(:n) )
!
    elsewhere
!
        minr(:n)  = miss
        maxr(:n)  = miss
        pminr(:n) = miss
        pmaxr(:n) = miss
        prob(:n)  = miss
!       
    end where
!
!   OUTPUT THE REQUIRED OPTIONAL ARGUMENTS.
!
    if ( present( min_ratio ) ) then
        min_ratio(:n) = minr(:n)
    end if
!
    if ( present( max_ratio ) ) then
        max_ratio(:n) = maxr(:n)
    end if
!
    if ( present( prob_min_ratio ) ) then
        prob_min_ratio(:n) = pminr(:n)
    end if
!
    if ( present( prob_max_ratio ) ) then
        prob_max_ratio(:n) = pmaxr(:n)
    end if
!
!
! END OF SUBROUTINE spctrm_ratio2_rm
! __________________________________
!
    end subroutine spctrm_ratio2_rm
!
! =========================================================================================
!
    subroutine spctrm_ratio3_rv( psvecn, psvecd, edofn, edofd, chi2_stat, prob )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_RATIO3 calculates an approximate critical probability value (e.g. p-value)
!   for testing the hypothesis of a common spectrum for two estimated spectra (e.g. the arguments
!   PSVECN, PSVECD). This approximate critical probability value is derived from the following
!   CHI2 statistic :
!
!       CHI2_STAT  =  ( 2/EDOFN + 2/EDOFD )**(-1)  [ sum k=1 to nf ]  log( PSVECN(k) / PSVECD(k) )**(2)
!
!   where nf = size(PSVECN) = size(PSVECD). In order to derive an approximate critical probability
!   value, it is assumed that CHI2_STAT has an approximate CHI2 distribution with nf degrees of freedom.
!
!
! Arguments
! _________
!
!
!   PSVECN             (INPUT) real(stnd), dimension(:)
!                      On entry, a real vector containing the Power Spectral Density (PSD)
!                      estimates of the first time series (e.g. the numerator of the ratio of
!                      the two estimated spectra).  
!
!                      All elements in PSVECN(:) must be greater than zero and size(PSVECN)
!                      must be greater or equal to 2.
!
!   PSVECD             (INPUT) real(stnd), dimension(:)
!                      On entry, a real vector containing the Power Spectral Density (PSD)
!                      estimates of the second time series (e.g. the denominator of the ratio of
!                      the two estimated spectra).  
!
!                      All elements in PSVECD(:) must be greater than zero and size(PSVECD)
!                      must be greater or equal to 2.
!
!                      PSVECD must verify:  size(PSVECD) = size(PSVECN) .
!
!   EDOFN              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the first
!                      estimated spectrum (e.g. the numerator of the ratio of the
!                      two estimated spectra).
!
!                      EDOFN must be greater than one.
!
!   EDOFD              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the second
!                      estimated spectrum (e.g. the denominator of the ratio of the
!                      two estimated spectra).
!
!                      EDOFD must be greater than one.
!
!   CHI2_STAT          (OUTPUT) real(stnd)
!                      On output, the CHI2 statistic which is assumed to follow a CHI2 distribution
!                      with size(PSVECN) degrees of freedom under the null hypothesis of a common
!                      spectrum.
!
!   PROB               (OUTPUT) real(stnd)
!                      On exit, the aproximate critical probability value (e.g. p-value) computed
!                      under the hypothesis that the two "true" underlying spectra are the same.
!                      PROB is calculated as the probability of obtaining a value greater or equal
!                      to CHI2_STAT under the hypothesis of a common spectrum for the two series.
!
!
! Further Details
! _______________
!
!   This statistical test relies on the assumptions that the different spectral ordinates
!   have the same sampling distribution and are independent of each other for each time series.
!   This means, in particular, that the spectral ordinates corresponding to the zero and Nyquist
!   frequencies must be excluded from the PSVECN and PSVECD vectors before calling SPCTRM_RATIO3
!   and that the two estimated spectra have not been obtained by smoothing the periodogram in the
!   frequency domain.
!
!   It is also assumed that the PSVECN and PSVECD realizations are independent.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert_eq, merror
    use Reals_Constants,   only : zero, one, two
    use Num_Constants,     only : nan
    use Logical_Constants, only : true
    use Char_Constants,    only : tseries_error70, tseries_error74
    use Prob_Procedures,   only : probq
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),  dimension(:) :: psvecn, psvecd
    real(stnd), intent(in)                :: edofn, edofd
!
    real(stnd), intent(out)             :: chi2_stat, prob
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: nf, nmask
!
    real(stnd)                            :: miss, con
    real(stnd),  dimension(size(psvecn))  :: lnratio
!
    logical(lgl),  dimension(size(psvecn)) :: mask
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_ratio3'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    nf =  assert_eq( int(size(psvecn),i4b) ,          &
                     int(size(psvecd),i4b) ,          &
                     name_proc )
!
    if ( nf<=1_i4b )     &
    call merror( name_proc//tseries_error74 )
!
    if ( one>=edofn .or. one>=edofd )    &
    call merror( name_proc//tseries_error70 )
!
    mask(:nf) = psvecn(:nf)>zero .and. psvecd(:nf)>zero
    nmask = count( mask(:nf) )
!
    if ( nmask>=2_i4b ) then
!
!       COMPUTE THE SQUARE OF THE LOG RATIO OF THE TWO SPECTRAL ESTIMATES.
!
        con = (one/(edofn-one)) - (one/(edofd-one))
!
        where( mask(:nf) )
            lnratio(:nf) = ( log( psvecn(:nf)/psvecd(:nf) ) - con )**2
        elsewhere
            lnratio(:nf) = zero
        end where
!
!       COMPUTE THE CHI2 STATISTIC AND THE ASSOCIATED PROBABILITY.
!
        con = (two/(edofn-one)) + (two/(edofd-one))
!
        chi2_stat = (one/con )*sum( lnratio(:nf) )
        prob      = probq( chi2_stat, nmask, true )
!
    else
!
        miss = nan()
!
        chi2_stat = miss
        prob      = miss
!       
    end if
!
!
! END OF SUBROUTINE spctrm_ratio3_rv
! __________________________________
!
    end subroutine spctrm_ratio3_rv
!
! =========================================================================================
!
    subroutine spctrm_ratio3_rm( psmatn, psmatd, edofn, edofd, chi2_stat, prob )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_RATIO3 calculates approximate critical probability values (e.g. p-values)
!   for testing the hypothesis of a common spectrum for two estimated multi-channel spectra (eg
!   the arguments PSMATN, PSMATD). Thess approximate critical probability values are derived
!   from the following CHI2 statistics :
!
!        CHI2_STAT(:n) =  ( 2/EDOFN + 2/EDOFD )**(-1)  [ sum k=1 to nf ]  log( PSMATN(:n,k) / PSMATD(:n,k) )**(2)
!
!   where n = size(PSMATN,1) = size(PSMATD,1) = size(CHI2_STAT) and  nf = size(PSMATN,2) = size(PSMATD,2).
!   In order to derive approximate critical probability values, it is assumed that each element of
!   CHI2_STAT(:n) has an approximate CHI2 distribution with nf degrees of freedom.
!
!
! Arguments
! _________
!
!
!   PSMATN             (INPUT) real(stnd), dimension(:,:)
!                      On entry, a real matrix containing the Power Spectral Density (PSD)
!                      estimates of the first multi-channel time series (e.g. the numerator of
!                      the ratio of the two estimated multi-channel spectra).  
!                      Each row of the real matrix PSMATN contains the estimated spectrum
!                      of the corresponding "row" of the first multi-channel times series.  
!
!                      All elements in PSMATN(:,:) must be greater than zero and size(PSMATN,2)
!                      must be greater or equal to 2.
!
!   PSMATD             (INPUT) real(stnd), dimension(:,:)
!                      On entry, a real matrix containing the Power Spectral Density (PSD)
!                      estimates of the second multi-channel time series (e.g. the denominator
!                      of the ratio of the two estimated multi-channel spectra).  
!                      Each row of the real matrix PSMATD contains the estimated spectrum
!                      of the corresponding "row" of the second multi-channel times series.  
!                      All elements in PSMATD(:,:) must be greater than zero and size(PSMATD,2)
!                      must be greater or equal to 2.
!
!                      PSMATD must also verify:
!
!                      - size(PSMATD,1) = size(PSMATN,1) ,
!                      - size(PSMATD,2) = size(PSMATN,2) .
!
!   EDOFN              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the first
!                      estimated spectrum (e.g. the numerator of the ratio of the
!                      two estimated spectra).
!
!                      EDOFN must be greater than one.
!
!   EDOFD              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the second
!                      estimated spectrum (e.g. the denominator of the ratio of the
!                      two estimated spectra).
!
!                      EDOFD must be greater than one.
!
!   CHI2_STAT          (OUTPUT) real(stnd), dimension(:)
!                      On output, the CHI2 statistics which are assumed to follow a CHI2 distribution
!                      with size(PSMATN,2) degrees of freedom under the null hypothesis of a common
!                      spectrum.
!
!                      CHI2_STAT must verify:  size(CHI2_STAT) = size(PSMATN,1) .
!
!   PROB               (OUTPUT) real(stnd), dimension(:)
!                      On exit, the aproximate critical probability values (e.g. p-values) computed
!                      under the hypothesis that the two "true" underlying multi-channel spectra are
!                      the same. Each element of PROB(:) is calculated as the probability of obtaining
!                      a value greater or equal to the corresponding element of CHI2_STAT(:) under the
!                      hypothesis of a common spectrum for the two (single) series.
!
!                      PROB must verify:  size(PROB) = size(PSMATN,1) .
!
!
! Further Details
! _______________
!
!   This statistical test relies on the assumptions that the different spectral ordinates
!   have the same sampling distribution and are independent of each other for each time series.
!   This means, in particular, that the spectral ordinates corresponding to the zero and Nyquist
!   frequencies must be excluded from the PSMATN and PSMATD matrices before calling SPCTRM_RATIO3
!   and that the two estimated multi-channel spectra have not been obtained by smoothing the
!   periodogram in the frequency domain.
!
!   It is also assumed that the PSMATN and PSMATD realizations are independent.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert_eq, merror
    use Reals_Constants,   only : zero, one, two
    use Num_Constants,     only : nan
    use Logical_Constants, only : true
    use Char_Constants,    only : tseries_error70, tseries_error75
    use Prob_Procedures,   only : probq
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),  dimension(:,:) :: psmatn, psmatd
    real(stnd), intent(in)                  :: edofn, edofd
!
    real(stnd), intent(out),  dimension(:)  :: chi2_stat, prob
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b)                             :: nf, n, i
    integer(i4b),  dimension(size(psmatn,1)) :: nmask
!
    real(stnd)                                             :: miss, con
    real(stnd),  dimension(size(psmatn,1),size(psmatn,2))  :: lnratio
!
    logical(lgl),  dimension(size(psmatn,1),size(psmatn,2)) :: mask
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_ratio3'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n =  assert_eq( int(size(psmatn,1),i4b) ,          &
                    int(size(psmatd,1),i4b) ,          &
                    int(size(chi2_stat),i4b) ,         &
                    int(size(prob),i4b) ,              &
                    name_proc )
!
    if ( n<1_i4b ) return
!
    nf =  assert_eq( int(size(psmatn,2),i4b) ,          &
                     int(size(psmatd,2),i4b) ,          &
                     name_proc )
!
    if ( nf<=1_i4b )     &
    call merror( name_proc//tseries_error75 )
!
    if ( one>=edofn .or. one>=edofd )   &
    call merror( name_proc//tseries_error70  )
!
    mask(:n,:nf) = psmatd(:n,:nf)>zero .and. psmatn(:n,:nf)>zero
    nmask(:n)    = count( mask(:n,:nf), dim=2 )
!
!   COMPUTE THE SQUARE OF THE LOG RATIO OF THE TWO SPECTRAL ESTIMATES.
!
    con = (one/(edofn-one)) - (one/(edofd-one))
!
    where( mask(:n,:nf) )
         lnratio(:n,:nf) = ( log( psmatn(:n,:nf)/psmatd(:n,:nf) ) - con )**2
    elsewhere
         lnratio(:n,:nf) = zero
    end where
!
!   COMPUTE THE CHI2 STATISTICS AND THE ASSOCIATED PROBABILITIES.
!
    con = (two/(edofn-one)) + (two/(edofd-one))
!
    chi2_stat(:n) = (one/con )*sum( lnratio(:n,:nf), dim=2 )
!
    miss = nan()
!
    i = nmask(1_i4b)
!
    if ( all( i==nmask(:n) ) ) then
!
        if ( i>1_i4b ) then
            prob(:n) = probq( chi2_stat(:n), i, true )
        else
            chi2_stat(:n) = miss
            prob(:n)      = miss
        end if
!
    else
!
        do i=1_i4b, n
!
            if ( nmask(i)>1_i4b ) then
                prob(i) = probq( chi2_stat(i), nmask(i), true )
            else
                chi2_stat(i) = miss
                prob(i)      = miss
            end if
!
        end do
!
    end if
!
!
! END OF SUBROUTINE spctrm_ratio3_rm
! __________________________________
!
    end subroutine spctrm_ratio3_rm
!
! =========================================================================================
!
    subroutine spctrm_ratio4_rv( psvecn, psvecd, edofn, edofd, range_stat, prob )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_RATIO4 calculates an approximate critical probability value (e.g. p-value)
!   for testing the hypothesis of a common shape for two estimated spectra (e.g. the arguments
!   PSVECN, PSVECD). This approximate critical probability value is derived from the following
!   RANGE statistic :
!                                                  
!     RANGE_STAT  =  ( 2/EDOFN + 2/EDOFD )**(-1/2) * ( maxval( logratio(:nf) ) - minval( logratio(:nf) ) )
!                                                  
!   where nf = size(PSVECN) = size(PSVECD) and logratio(:nf) is given as
!
!              logratio(:nf) = log( PSVECN(:nf) / PSVECD(:nf) )
!
!   In order to derive an approximate critical probability value, it is assumed that the elements of
!   the vector logratio(:nf) are independent and follow approximately a normal distribution with mean
!   (1/EDOFN) - (1/EDOFD) and variance  (2/EDOFN) + (2/EDOFD). Than, the distribution of the statistic
!   RANGE_STAT may be approximated by the distribution function of the range of nf independent normal
!   random variables with mean and variance as specified above.
!
!
! Arguments
! _________
!
!
!   PSVECN             (INPUT) real(stnd), dimension(:)
!                      On entry, a real vector containing the Power Spectral Density (PSD)
!                      estimates of the first time series (e.g. the numerator of the ratio of
!                      the two estimated spectra).  
!
!                      All elements in PSVECN(:) must be greater than zero and size(PSVECN)
!                      must be greater or equal to 2.
!
!   PSVECD             (INPUT) real(stnd), dimension(:)
!                      On entry, a real vector containing the Power Spectral Density (PSD)
!                      estimates of the second time series (e.g. the denominator of the ratio of
!                      the two estimated spectra).  
!
!                      All elements in PSVECD(:) must be greater than zero and size(PSVECD)
!                      must be greater or equal to 2.
!
!                      PSVECD must also verify:  size(PSVECD) = size(PSVECN) .
!
!   EDOFN              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the first
!                      estimated spectrum (e.g. the numerator of the ratio of the
!                      two estimated spectra).
!
!                      EDOFN must be greater than one.
!
!   EDOFD              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the second
!                      estimated spectrum (e.g. the denominator of the ratio of the
!                      two estimated spectra).
!
!                      EDOFD must be greater than one.
!
!   RANGE_STAT         (OUTPUT) real(stnd)
!                      On output, the range statistic which is assumed to follow the distribution
!                      of the range of nf=size(PSVECN) independent standard normal variates under
!                      the null hypothesis of a common shape spectrum.
!
!   PROB               (OUTPUT) real(stnd)
!                      On exit, the aproximate critical probability value (e.g. p-value) computed
!                      under the hypothesis that the two "true" underlying spectra have the same shape.
!                      PROB is calculated as the probability of obtaining a value greater or equal
!                      to RANGE_STAT under the hypothesis of a common shape spectrum for the two series.
!
!
! Further Details
! _______________
!
!   This statistical test relies on the assumptions that the different spectral ordinates
!   have the same sampling distribution and are independent of each other for each time series.
!   This means, in particular, that the spectral ordinates corresponding to the zero and Nyquist
!   frequencies must be excluded from the PSVECN and PSVECD vectors before calling SPCTRM_RATIO4
!   and that the two estimated spectra have not been obtained by smoothing the periodogram in the
!   frequency domain.
!
!   It is also assumed that the PSVECN and PSVECD realizations are independent.
!
!   For more details, see:
!
!   (1) Coates, D.S., and Diggle, P.J., 1986:
!             Tests for comparing two estimated spectral densities.
!             J. Time series Analysis, Vol. 7, pp. 7-20 .
!
!   (2) Potscher, B.,M., and Reschenhofer, E., 1988:
!             Discriminating between two spectral densities in case of replicated observations.
!             J. Time series Analysis, Vol. 9, pp. 221-224 .
!
!   (3) Potscher, B.,M., and Reschenhofer, E., 1989:
!             Distribution of the Coates-Diggle test statistic in case of replicated observations.
!             Statistics, Vol. 20, pp. 417-421 .
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert_eq, merror
    use Reals_Constants,   only : zero, one, two
    use Num_Constants,     only : nan
    use Char_Constants,    only : tseries_error70, tseries_error74
    use Prob_Procedures,   only : rangen
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),  dimension(:) :: psvecn, psvecd
    real(stnd), intent(in)                :: edofn, edofd
!
    real(stnd), intent(out)               :: range_stat, prob
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: nf, nmask
!
    real(stnd)                            :: con, maxr, minr
    real(stnd),  dimension(size(psvecn))  :: ratio
!
    logical(lgl),  dimension(size(psvecn)) :: mask
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_ratio4'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    nf =  assert_eq( int(size(psvecn),i4b) ,          &
                     int(size(psvecd),i4b) ,          &
                     name_proc )
!
    if ( nf<=1_i4b )     &
    call merror( name_proc//tseries_error74 )
!
    if ( one>=edofn .or. one>=edofd )    &
    call merror( name_proc//tseries_error70 )
!
    mask(:nf) = psvecn(:nf)>zero .and. psvecd(:nf)>zero
    nmask = count( mask(:nf) )
!
    if ( nmask>=2_i4b ) then
!
!       COMPUTE THE RATIO OF THE TWO SPECTRAL ESTIMATES.
!
        where( mask(:nf) )
            ratio(:nf) = psvecn(:nf)/psvecd(:nf)
        elsewhere
            ratio(:nf) = one
        end where
!
!       COMPUTE THE MIN AND MAX OF THE RATIO OF THE TWO SPECTRAL ESTIMATES.
!
        minr = minval( ratio(:nf), mask=mask(:nf) )
        maxr = maxval( ratio(:nf), mask=mask(:nf) )
!
!       COMPUTE THE RANGE STATISTIC AND THE ASSOCIATED PROBABILITY.
!
        con        = sqrt( (two/(edofn-one)) + (two/(edofd-one)) )
        range_stat = ( log( maxr ) - log( minr ) )/con
!
        if ( range_stat/=zero ) then
            prob  = one - rangen( range_stat, nmask )
        else
            prob  = one
        end if
!
    else
!
        con = nan()
!
        range_stat = con
        prob       = con
!       
    end if
!
!
! END OF SUBROUTINE spctrm_ratio4_rv
! __________________________________
!
    end subroutine spctrm_ratio4_rv
!
! =========================================================================================
!
    subroutine spctrm_ratio4_rm( psmatn, psmatd, edofn, edofd, range_stat, prob )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_RATIO4 calculates approximate critical probability values (e.g. p-values)
!   for testing the hypothesis of a common shape for two estimated multi-channel spectra (e.g. the
!   arguments PSMATN, PSMATD). These approximate critical probability values are derived from
!   the following range statistics :
!                                                  
!     RANGE_STAT(:n)  =  ( 2/EDOFN + 2/EDOFD )**(-1/2) * ( maxval( logratio(:n,:nf), dim=2 ) - minval( logratio(:n,:nf), dim=2  ) )
!                                                  
!   where n = size(PSMATN,1) = size(PSMATD,1), nf = size(PSMATN,2) = size(PSMATD,2) and
!   logratio(:n,:nf) is given as
!
!              logratio(:n,:nf) = log( PSMATN(:n,:nf) / PSMATD(:n,:nf) )
!
!   In order to derive approximate critical probability values, it is assumed that the elements of
!   the vectors logratio(i,:nf), for i=1 to n, are independent and follow approximately a normal 
!   distribution with mean (1/EDOFN) - (1/EDOFD) and variance  (2/EDOFN) + (2/EDOFD).
!   Than, the distribution of the statistics RANGE_STAT(:n) may be approximated by the distribution
!   function of the range of nf independent normal random variables with mean and variance
!   as specified above.
!
!
! Arguments
! _________
!
!
!   PSMATN             (INPUT) real(stnd), dimension(:,:)
!                      On entry, a real matrix containing the Power Spectral Density (PSD)
!                      estimates of the first multi-channel time series (e.g. the numerator of
!                      the ratio of the two estimated multi-channel spectra).  
!                      Each row of the real matrix PSMATN contains the estimated spectrum
!                      of the corresponding "row" of the first multi-channel times series.  
!
!                      All elements in PSMATN(:,:) must be greater than zero and size(PSMATN,2)
!                      must be greater or equal to 2.
!
!   PSMATD             (INPUT) real(stnd), dimension(:,:)
!                      On entry, a real matrix containing the Power Spectral Density (PSD)
!                      estimates of the second multi-channel time series (e.g. the denominator
!                      of the ratio of the two estimated multi-channel spectra).  
!                      Each row of the real matrix PSMATD contains the estimated spectrum
!                      of the corresponding "row" of the second multi-channel times series.  
!
!                      All elements in PSMATD(:,:) must be greater than zero and size(PSMATD,2)
!                      must be greater or equal to 2.
!
!                      PSMATD must also verify:
!
!                      - size(PSMATD,1) = size(PSMATN,1) ,
!                      - size(PSMATD,2) = size(PSMATN,2) .
!
!   EDOFN              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the first
!                      estimated spectrum (e.g. the numerator of the ratio of the
!                      two multi-channel estimated spectra).
!
!                      EDOFN must be greater than one.
!
!   EDOFD              (INPUT) real(stnd)
!                      On exit, the equivalent number of degrees of freedom of the second
!                      estimated spectrum (e.g. the denominator of the ratio of the
!                      two multi-channel estimated spectra).
!
!                      EDOFD must be greater than one.
!
!   RANGE_STAT         (OUTPUT) real(stnd), dimension(:)
!                      On output, the range statistics which are assumed to follow the distribution
!                      of the range of nf=size(PSMATN,2) independent standard normal variates under
!                      the null hypothesis of a common shape spectrum.
!
!                      RANGE_STAT must verify:  size(RANGE_STAT) = size(PSMATN,1) .
!
!   PROB               (OUTPUT) real(stnd), dimension(:)
!                      On exit, the aproximate critical probability values (e.g. p-values) computed
!                      under the hypothesis that the two "true" underlying multi-channel spectra have
!                      the same shape.
!
!                      Each element of PROB(:) is calculated as the probability of obtaining a value greater or equal
!                      to the corresponding element of RANGE_STAT(:) under the hypothesis of a common shape spectrum
!                      for the two multi-channel series.
!
!
! Further Details
! _______________
!
!   This statistical test relies on the assumptions that the different spectral ordinates
!   have the same sampling distribution and are independent of each other for each multi-channel
!   time series.
!   This means, in particular, that the spectral ordinates corresponding to the zero and Nyquist
!   frequencies must be excluded from the PSMATN and PSMATD matrices before calling SPCTRM_RATIO4
!   and that the two estimated multi-channel spectra have not been obtained by smoothing the periodogram
!   in the frequency domain.
!
!   It is also assumed that the PSMATN and PSMATD multi-channel realizations are independent.
!
!   For more details, see:
!
!   (1) Coates, D.S., and Diggle, P.J., 1986:
!             Tests for comparing two estimated spectral densities.
!             J. Time series Analysis, Vol. 7, pp. 7-20 .
!
!   (2) Potscher, B.,M., and Reschenhofer, E., 1988:
!             Discriminating between two spectral densities in case of replicated observations.
!             J. Time series Analysis, Vol. 9, pp. 221-224 .
!
!   (3) Potscher, B.,M., and Reschenhofer, E., 1989:
!             Distribution of the Coates-Diggle test statistic in case of replicated observations.
!             Statistics, Vol. 20, pp. 417-421 .
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert_eq, merror
    use Reals_Constants,   only : zero, one, two
    use Num_Constants,     only : nan
    use Char_Constants,    only : tseries_error70, tseries_error75
    use Prob_Procedures,   only : rangen
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),  dimension(:,:) :: psmatn, psmatd
    real(stnd), intent(in)                  :: edofn, edofd
!
    real(stnd), intent(out),  dimension(:)  :: range_stat, prob
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b)                             :: nf, n, i
    integer(i4b),  dimension(size(psmatn,1)) :: nmask
!
    real(stnd)                                             :: con, miss
    real(stnd),  dimension(size(psmatn,1))                 :: maxr, minr
    real(stnd),  dimension(size(psmatn,1),size(psmatn,2))  :: ratio
!
    logical(lgl),  dimension(size(psmatn,1),size(psmatn,2)) :: mask
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_ratio4'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n =  assert_eq( int(size(psmatn,1),i4b) ,          &
                    int(size(psmatd,1),i4b) ,          &
                    int(size(range_stat),i4b) ,        &
                    int(size(prob),i4b) ,              &
                    name_proc )
!
    nf =  assert_eq( int(size(psmatn,2),i4b) ,          &
                     int(size(psmatd,2),i4b) ,          &
                     name_proc )
!
    if ( nf<=1_i4b )     &
    call merror( name_proc//tseries_error75 )
!
    if ( one>=edofn .or. one>=edofd )    &
    call merror( name_proc//tseries_error70 )
!
    mask(:n,:nf) = psmatd(:n,:nf)>zero .and. psmatn(:n,:nf)>zero
    nmask(:n)    = count( mask(:n,:nf), dim=2 )
!
    con  = sqrt( (two/(edofn-one)) + (two/(edofd-one)) )
!
    if ( all( nmask(:n)==nf ) ) then
!
!       COMPUTE THE RATIO OF THE TWO SPECTRAL ESTIMATES.
!
        ratio(:n,:nf) = psmatn(:n,:nf)/psmatd(:n,:nf)
!
!       COMPUTE THE MIN AND MAX OF THE RATIO OF THE TWO SPECTRAL ESTIMATES.
!
        minr(:n) = minval( ratio(:n,:nf) , dim=2 )
        maxr(:n) = maxval( ratio(:n,:nf) , dim=2 )
!
!       COMPUTE THE RANGE STATISTICS AND THE ASSOCIATED PROBABILITIES.
!
        range_stat(:n) =  (one/con)*( log( maxr(:n) ) - log( minr(:n) ) )
!
        where( range_stat(:n)/=zero )
            maxr(:n) = range_stat(:n)
        elsewhere
            maxr(:n) = one
        end where
!
        prob(:n)  = one - rangen( maxr(:n), nf )
!
        where( range_stat(:n)==zero ) prob(:n)  = one
!
    else
!
!       COMPUTE THE RATIO OF THE TWO SPECTRAL ESTIMATES.
!
         where( mask(:n,:nf) )
             ratio(:n,:nf) = psmatn(:n,:nf)/psmatd(:n,:nf)
        elsewhere
             ratio(:n,:nf) = one
        end where    
!
!       COMPUTE THE MIN AND MAX OF THE RATIO OF THE TWO SPECTRAL ESTIMATES.
!
        minr(:n) = minval( ratio(:n,:nf) , dim=2, mask=mask(:n,:nf)  )
        maxr(:n) = maxval( ratio(:n,:nf) , dim=2, mask=mask(:n,:nf)  )
!
!       COMPUTE THE RANGE STATISTICS AND THE ASSOCIATED PROBABILITIES.
!
        miss = nan()
!
        do i=1_i4b, n
!
            if ( nmask(i)>1_i4b ) then
!
                range_stat(i) =  ( log( maxr(i) ) - log( minr(i) ) )/con
                if ( range_stat(i)/=zero ) then
                    prob(i) = one - rangen( range_stat(i), nmask(i) )
                else
                    prob(i) = one
                end if
!
            else
!
                range_stat(i) = miss
                prob(i)       = miss
!
            end if
!
        end do
!       
    end if
!
!
! END OF SUBROUTINE spctrm_ratio4_rm
! __________________________________
!
    end subroutine spctrm_ratio4_rm
!
! =========================================================================================
!
    subroutine spctrm_diff_rv( psvec1, psvec2, ks_stat, prob, nrep, norm, initseed  )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_DIFF calculates an approximate critical probability value (e.g. p-value)
!   for testing the hypothesis of a common shape for two estimated spectra (e.g. the arguments PSVEC1
!   and PSVEC2). This approximate critical probability value is derived from the following
!   Kolmogorov-Smirnov statistic (e.g. the KS_STAT argument) :
!
!                           D =   [ sup m=1 to nf ]  | F1(m) - F2(m) |
!
!   where nf = size(PSVEC1) = size(PSVEC2), F1(:) and F2(:) are the normalized cumulative periodograms
!   computed from the estimated spectra PSVEC1(:) and PSVEC2(:). The distribution of D under the null 
!   hypothesis of a common shape for the spectra of the two series is approximated by calculating
!   D for some large number (e.g. the NREP argument) of random interchanges of periodogram ordinates
!   at each frequency for the two estimated spectra (e.g. the arguments PSVEC1(:) and PSVEC2(:)).
!
! Arguments
! _________
!
!
!   PSVEC1   (INPUT) real(stnd), dimension(:)
!            On entry, a real vector containing the Power Spectral Density (PSD)
!            estimates of the first time series.  
!
!            size(PSVECN) must be greater or equal to 10.
!
!   PSVEC2   (INPUT) real(stnd), dimension(:)
!            On entry, a real vector containing the Power Spectral Density (PSD)
!            estimates of the second time series.  
!
!            PSVEC2 must verify:  size(PSVEC2) = size(PSVEC1) .
!
!   KS_STAT  (OUTPUT) real(stnd)
!            On output, the Kolmogorov-Smirnov statistic.
!
!   PROB     (OUTPUT) real(stnd)
!            On exit, the aproximate critical probability value (e.g. p-value) computed
!            under the hypothesis that the two "true" underlying spectra have the same shape.
!            PROB is calculated as the probability of obtaining a value greater or equal
!            to KS_STAT under the hypothesis of a common shape for the spectra of the 
!            two series.
!
!   NREP     (INPUT, OPTIONAL) integer(i4b)
!            On entry, when argument NREP is present, NREP specifies the number of random interchanges
!            of periodogram ordinates at each frequency in order to approximate the randomization 
!            distribution of KS_STAT under the null hypothesis.
!
!            The default is 99.
!
!   NORM     (INPUT, OPTIONAL) logical(lgl)
!            On entry, if:
!
!            - NORM=true, KS_STAT is calculated from normalized cumulative periodograms.
!            - NORM=false KS_STAT is calculated from non-normalized cumulative periodograms.
!              In that case the null hypothesis is that the spectra of the two time series is
!              the same.
!                                    
!            The default is NORM=true.
!
!   INITSEED (INPUT, OPTIONAL) logical(lgl)
!            On entry, if INITSEED=true, a call to RANDOM_SEED_() without arguments
!            is done in the subroutine, in order to initiate a non-repeatable reset
!            of the seed used by the STATPACK random generator.
!
!            The default is INITSEED=false.
!
!
! Further Details
! _______________
!
!   This statistical randomization test relies on the assumptions that the different spectral ordinates
!   have the same sampling distribution and are independent of each other. This means, in
!   particular, that the spectral ordinates corresponding to the zero and Nyquist frequencies
!   must be excluded from the PSVEC1 and PSVEC2 vectors before calling SPCTRM_DIFF and
!   that the two estimated spectra have not been obtained by smoothing the periodogram
!   in the frequency domain.
!
!   This randomization test could only be used to compare two periodograms or two spectral
!   estimates computed as the the average of, say, r periodograms for each time series.
!
!   For more details, see:
!
!   (1) Diggle, P.J., and Fisher, N.I., 1991:
!               Nonparametric comparison of cumulative periodograms,
!               Applied Statistics, Vol. 40, No 3, pp. 423-434.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert_eq, merror, cumsum
    use Reals_Constants,   only : zero, half, one
    use Logical_Constants, only : true
    use Char_Constants,    only : stat_error8, tseries_error71, tseries_error73
    use Random,            only : random_seed_, random_number_
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),  dimension(:) :: psvec1, psvec2
    real(stnd), intent(out)               :: ks_stat, prob
!
    integer(i4b), intent(in), optional :: nrep
!
    logical(lgl), intent(in),  optional :: norm, initseed
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: nf, nfm1, nrep2, nge, i
!
    real(stnd)                            :: ks
    real(stnd),  dimension(size(psvec1))  :: psveca, psvecb, f1, f2
!
    logical(lgl)                           :: norm2
!
#ifdef _OPENMP
    integer(i4b) :: nge2, j
!
    logical      :: test_par
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_diff'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    nf =  assert_eq( int(size(psvec1),i4b) ,          &
                     int(size(psvec2),i4b) ,          &
                     name_proc )
!
    if ( nf<10_i4b  )    &
    call merror( name_proc//tseries_error71  )
!
    if ( all( psvec1(:nf)<=zero ) .or. all( psvec2(:nf)<=zero )  )    &
    call merror( name_proc//tseries_error73 )
!
    if ( present(initseed) ) then
!
        if ( initseed ) then
!
!           INITIALIZE THE RANDOM GENERATOR IF REQUIRED.
!
!$OMP SINGLE
            call random_seed_()
!$OMP END SINGLE
!
        end if
!
    end if
!
    nrep2 = 99_i4b
!
    if ( present(nrep) ) then
!
        if ( nrep<=0_i4b )        &
        call merror( name_proc//stat_error8 )
!
        nrep2 = nrep
!
    end if
!
    norm2 = true
!
    if ( present(norm) ) then
        norm2 = norm
    end if
!
!   COMPUTE THE TWO CUMULATIVE PERIODOGRAMS.
!
    f1(:nf) = cumsum( psvec1(:nf) )
    f2(:nf) = cumsum( psvec2(:nf) )
!
!   NORMALIZED THE CUMULATIVE PERIODOGRAMS IF REQUIRED.
!
    if ( norm2 ) then
!
        nfm1 = nf - 1_i4b
!
        ks = one/f1(nf)
        f1(:nfm1) = f1(:nfm1)*ks
!
        ks = one/f2(nf)
        f2(:nfm1) = f2(:nfm1)*ks
!
    else
!
        nfm1 = nf
!
    end if
!
!   COMPUTE THE KOLMOGOROV-SMIRNOV STATISTIC.
!
    ks_stat = maxval( abs( f1(:nfm1) - f2(:nfm1) ) )
!
    nge = 0_i4b
!
#ifdef _OPENMP
!
    i = omp_get_num_procs()
    j = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() ) .and. i>1 .and. j>1
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(psveca,psvecb,f1,f2,ks,nge2,i)
!
        nge2 = 0_i4b
!
!$OMP DO SCHEDULE(STATIC) 
!
        do i=1_i4b, nrep2
!
!           GENERATE A RANDOM NUMBER UNIFORM SEQUENCE.
!
            call random_number_(psveca(:nf))
!
!           RANDOM INTERCHANGES OF PERIODOGRAM ORDINATES.
!
            where( psveca(:nf)<=half )
                psveca(:nf) = psvec1(:nf)
                psvecb(:nf) = psvec2(:nf)
            elsewhere
                psveca(:nf) = psvec2(:nf)
                psvecb(:nf) = psvec1(:nf)
            end where
!
!           COMPUTE THE CUMULATIVE PERIODOGRAMS
!
            f1(:nf) = cumsum( psveca(:nf) )
            f2(:nf) = cumsum( psvecb(:nf) )
!
!           NORMALIZED THE CUMULATIVE PERIODOGRAMS IF REQUIRED.
!
            if ( norm2 ) then
!
                ks = one/f1(nf)
                f1(:nfm1) = f1(:nfm1)*ks
!
                ks = one/f2(nf)
                f2(:nfm1) = f2(:nfm1)*ks
!
            end if
!
!           COMPUTE THE KOLMOGOROV-SMIRNOV STATISTIC.
!
            ks = maxval( abs( f1(:nfm1) - f2(:nfm1) ) )
!
            if ( ks>=ks_stat ) nge2 = nge2 + 1_i4b
!
        end do
!
!$OMP END DO NOWAIT
!
!$OMP ATOMIC
        nge = nge + nge2
!
!$OMP END PARALLEL
!
    else
!
#endif
!
        do i=1_i4b, nrep2
!
!           GENERATE A RANDOM NUMBER SEQUENCE.
!
            call random_number_(psveca(:nf))
!
!           RANDOM INTERCHANGES OF PERIODOGRAM ORDINATES.
!
            where( psveca(:nf)<=half )
                psveca(:nf) = psvec1(:nf)
                psvecb(:nf) = psvec2(:nf)
            elsewhere
                psveca(:nf) = psvec2(:nf)
                psvecb(:nf) = psvec1(:nf)
            end where
!
!           COMPUTE THE CUMULATIVE PERIODOGRAMS
!
            f1(:nf) = cumsum( psveca(:nf) )
            f2(:nf) = cumsum( psvecb(:nf) )
!
!           NORMALIZED THE CUMULATIVE PERIODOGRAMS IF REQUIRED.
!
            if ( norm2 ) then
!
                ks = one/f1(nf)
                f1(:nfm1) = f1(:nfm1)*ks
!
                ks = one/f2(nf)
                f2(:nfm1) = f2(:nfm1)*ks
!
            end if
!
!           COMPUTE THE KOLMOGOROV-SMIRNOV STATISTIC.
!
            ks = maxval( abs( f1(:nfm1) - f2(:nfm1) ) )
!
!           CALCULATE THE SIGNIFICANCE PROBABILITY OF THE OBSERVED
!           KOLMOGOROV-SMIRNOV STATISTIC.
!
            if ( ks>=ks_stat ) nge = nge + 1_i4b
!
        end do
!
#ifdef _OPENMP
    end if
#endif
!
!   COMPUTE THE SIGNIFICANCE LEVEL.
!
    prob = real( nge+1_i4b, stnd )/real( nrep2+1_i4b, stnd )
!
!
! END OF SUBROUTINE spctrm_diff_rv
! ________________________________
!
    end subroutine spctrm_diff_rv
!
! =========================================================================================
!
    subroutine spctrm_diff_rm( psmat1, psmat2, ks_stat, prob, nrep, norm, initseed  )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_DIFF calculates approximate critical probability values (e.g. p-value)
!   for testing the hypothesis of a common shape for two estimated multi-channel spectra (e.g. the
!   arguments PSMAT1 and PSMAT2). These approximate critical probability values are derived from
!   the following Kolmogorov-Smirnov statistics (e.g. the KS_STAT(:) argument) :
!
!           D(j) =   [ sup m=1 to nf ] | F1(j,m) - F2(j,m) |      for j=1, ... , size(PSMAT1,1)
!
!   where nf = size(PSMAT1,2) = size(PSMAT2,2), F1(:,:) and F2(:,:) are the normalized cumulative periodograms
!   computed from the estimated spectra PSMAT1(:,:) and PSMAT2(:,:). The distribution of D under the null 
!   hypothesis of a common shape for the spectra of the two multi-channel series is approximated by calculating
!   D for some large number (e.g. the NREP argument) of random interchanges of periodogram ordinates
!   at each frequency for the two estimated multi-channel spectra (e.g. the arguments PSMAT1(:,:) and PSMAT2(:,:)).
!
!
! Arguments
! _________
!
!
!   PSMAT1   (INPUT) real(stnd), dimension(:,:)
!            On entry, a real matrix containing the Power Spectral Density (PSD)
!            estimates of the first multi-channel time series.  
!            Each row of the real matrix PSMAT1 contains the estimated spectrum
!            of the corresponding "row" of the first multi-channel times series.  
!
!            size(PSMATN,2) must be greater or equal to 10.
!
!   PSMAT2   (INPUT) real(stnd), dimension(:,:)
!            On entry, a real matrix containing the Power Spectral Density (PSD)
!            estimates of the second multi-channel time series.  
!            Each row of the real matrix PSMAT2 contains the estimated spectrum
!            of the corresponding "row" of the second multi-channel times series.  
!
!            PSMAT2 must verify:
!
!            - size(PSMAT2,1) = size(PSMAT1,1) ,
!            - size(PSMAT2,2) = size(PSMAT1,2) .
!
!   KS_STAT  (OUTPUT) real(stnd), dimension(:)
!            On output, the Kolmogorov-Smirnov statistics for the multi-channel 
!            times series .
!
!            KS_STAT must verify:  size(KS_STAT) = size(PSMAT1,1) .
!
!   PROB     (OUTPUT) real(stnd), dimension(:)
!            On exit, the aproximate critical probability values (e.g. p-values) computed
!            under the hypothesis that the two "true" underlying multi-channel spectra have
!            the same shape.
!
!            PROB is calculated as the probability of obtaining a value greater or equal
!            to KS_STAT under the hypothesis of a common shape for the spectra of the 
!            two series.
!
!            PROB must verify:  size(PROB) = size(PSMAT1,1) .
!
!   NREP     (INPUT, OPTIONAL) integer(i4b)
!            On entry, when argument NREP is present, NREP specifies the number of random interchanges
!            of periodogram ordinates at each frequency in order to approximate the randomization 
!            distribution of KS_STAT under the null hypothesis.
!
!            The default is 99.
!
!   NORM     (INPUT, OPTIONAL) logical(lgl)
!            On entry, if:
!
!            - NORM=true, KS_STAT is calculated from normalized cumulative periodograms.
!            - NORM=false KS_STAT is calculated from non-normalized cumulative periodograms.
!              In that case the null hypothesis is that the spectra of the two multi-channel
!              time series is the same.
!
!            The default is NORM=true.
!
!   INITSEED (INPUT, OPTIONAL) logical(lgl)
!            On entry, if INITSEED=true, a call to RANDOM_SEED_() without arguments
!            is done in the subroutine, in order to initiate a non-repeatable reset
!            of the seed used by the STATPACK random generator.
!
!            The default is INITSEED=false.
!
!
! Further Details
! _______________
!
!   This statistical randomization test relies on the assumptions that the different spectral ordinates
!   have the same sampling distribution and are independent of each other. This means, in
!   particular, that the spectral ordinates corresponding to the zero and Nyquist frequencies
!   must be excluded from the PSMAT1 and PSMAT2 matrices before calling SPCTRM_DIFF and
!   that the two estimated multi-channel spectra have not been obtained by smoothing the
!   periodograms in the frequency domain.
!
!   This randomization test could only be used to compare two periodograms or two spectral
!   estimates computed as the the average of, say, r periodograms for each time series.
!
!   For more details see:
!
!   (1) Diggle, P.J., and Fisher, N.I., 1991:
!               Nonparametric comparison of cumulative periodograms,
!               Applied Statistics, Vol. 40, No 3, pp. 423-434.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert_eq, merror
    use Reals_Constants,   only : zero, half, one
    use Logical_Constants, only : true
    use Char_Constants,    only : stat_error8, tseries_error72, tseries_error73
    use Random,            only : random_seed_, random_number_
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),   dimension(:,:) :: psmat1, psmat2
    real(stnd), intent(out),  dimension(:)   :: ks_stat, prob
!
    integer(i4b), intent(in), optional :: nrep
!
    logical(lgl), intent(in),  optional :: norm, initseed
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b)                             :: n, nf, nfm1, nrep2, i, j
    integer(i4b),  dimension(size(psmat1,1)) :: nge
!
    real(stnd),  dimension(size(psmat1,1))                 :: ks
    real(stnd),  dimension(size(psmat1,1),size(psmat1,2))  :: f1, f2
!
    logical(lgl) :: norm2
!
#ifdef _OPENMP
    integer(i4b),  dimension(size(psmat1,1)) :: nge2
!
    logical      :: test_par
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_diff'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n =  assert_eq( int(size(psmat1,1),i4b) ,          &
                    int(size(psmat2,1),i4b) ,          &
                    int(size(ks_stat),i4b) ,           &
                    int(size(prob),i4b) ,              &
                    name_proc )
!
    if ( n<1_i4b ) return
!
    nf =  assert_eq( int(size(psmat1,2),i4b) ,          &
                     int(size(psmat2,2),i4b) ,          &
                     name_proc )
!
    if ( nf<10_i4b  )    &
    call merror( name_proc//tseries_error72 )
!
    if ( any( all( psmat1(:n,:nf)<=zero, dim=2 ) .or. all( psmat2(:n,:nf)<=zero, dim=2 )  ) )    &
    call merror( name_proc//tseries_error73 )
!
    if ( present(initseed) ) then
!
        if ( initseed ) then
!
!           INITIALIZE THE RANDOM GENERATOR IF REQUIRED.
!
!$OMP SINGLE
            call random_seed_()
!$OMP END SINGLE
!
        end if
!
    end if
!
    nrep2 = 99_i4b
!
    if ( present(nrep) ) then
        if ( nrep<=0_i4b )        &
        call merror( name_proc//stat_error8 )
        nrep2 = nrep
    end if
!
    norm2 = true
!
    if ( present(norm) ) then
        norm2 = norm
    end if
!
!   COMPUTE THE TWO CUMULATIVE PERIODOGRAMS.
!
    f1(:n,1_i4b) = psmat1(:n,1_i4b)
!
    do j = 2_i4b, nf
        f1(:n,j) = f1(:n,j-1_i4b) + psmat1(:n,j)
    end do
!
    f2(:n,1_i4b) = psmat2(:n,1_i4b)
!
    do j = 2_i4b, nf
        f2(:n,j) = f2(:n,j-1_i4b) + psmat2(:n,j)
    end do
!
!   NORMALIZED THE CUMULATIVE PERIODOGRAMS IF REQUIRED.
!
    if ( norm2 ) then
!
        nfm1 = nf - 1_i4b
!
        ks(:n) = one/f1(:n,nf)
!
        do j = 1_i4b, nfm1
            f1(:n,j) = f1(:n,j)*ks(:n)
        end do
!
        ks(:n) = one/f2(:n,nf)
!
        do j = 1_i4b, nfm1
            f2(:n,j) = f2(:n,j)*ks(:n)
        end do
!
    else
!
        nfm1 = nf
!
    end if
!
!   COMPUTE THE KOLMOGOROV-SMIRNOV STATISTICS.
!
    ks_stat(:n) = maxval( abs( f1(:n,:nfm1) - f2(:n,:nfm1) ), dim=2 )
!
    nge(:n) = 0_i4b
!
#ifdef _OPENMP
!
    i = omp_get_num_procs()
    j = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() ) .and. i>1 .and. j>1
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(f1,f2,ks,nge2,i,j)
!
        nge2(:n) = 0_i4b
!
!$OMP DO SCHEDULE(STATIC) 
!
        do i=1_i4b, nrep2
!
!           GENERATE A RANDOM NUMBER UNIFORM SEQUENCE.
!
            call random_number_( f1(:n,:nf) )
!
!           RANDOM INTERCHANGES OF PERIODOGRAM ORDINATES.
!
            where( f1(:n,:nf)<=half )
                f1(:n,:nf) = psmat1(:n,:nf)
                f2(:n,:nf) = psmat2(:n,:nf)
            elsewhere
                f1(:n,:nf) = psmat2(:n,:nf)
                f2(:n,:nf) = psmat1(:n,:nf)
            end where
!
!           COMPUTE THE TWO CUMULATIVE PERIODOGRAMS.
!
            do j = 2_i4b, nf
                f1(:n,j) = f1(:n,j-1_i4b) + f1(:n,j)
            end do
!
            do j = 2_i4b, nf
                f2(:n,j) = f2(:n,j-1_i4b) + f2(:n,j)
            end do
!
!           NORMALIZED THE CUMULATIVE PERIODOGRAMS IF REQUIRED.
!
            if ( norm2 ) then
!
                ks(:n) = one/f1(:n,nf)
!
                do j = 1_i4b, nfm1
                    f1(:n,j) = f1(:n,j)*ks(:n)
                end do
!
                ks(:n) = one/f2(:n,nf)
!
                do j = 1_i4b, nfm1
                    f2(:n,j) = f2(:n,j)*ks(:n)
                end do
!
            end if
!
!           COMPUTE THE KOLMOGOROV-SMIRNOV STATISTICS.
!
            ks(:n) = maxval( abs( f1(:n,:nfm1) - f2(:n,:nfm1) ), dim=2 )
!
!           CALCULATE THE SIGNIFICANCE PROBABILITIES OF THE OBSERVED
!           KOLMOGOROV-SMIRNOV STATISTICS.
!
            where( ks(:n)>=ks_stat(:n) ) nge2(:n) = nge2(:n) + 1_i4b
!
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL  (diff)
!
        nge(:n) = nge(:n) + nge2(:n)
!
!$OMP END CRITICAL  (diff)
!
!$OMP END PARALLEL
!
    else
!
#endif
!
        do i=1_i4b, nrep2
!
!           GENERATE A RANDOM NUMBER SEQUENCE.
!
            call random_number_( f1(:n,:nf) )
!
!           RANDOM INTERCHANGES OF PERIODOGRAM ORDINATES.
!
            where( f1(:n,:nf)<=half )
                f1(:n,:nf) = psmat1(:n,:nf)
                f2(:n,:nf) = psmat2(:n,:nf)
            elsewhere
                f1(:n,:nf) = psmat2(:n,:nf)
                f2(:n,:nf) = psmat1(:n,:nf)
            end where
!
!           COMPUTE THE TWO CUMULATIVE PERIODOGRAMS.
!
            do j = 2_i4b, nf
                f1(:n,j) = f1(:n,j-1_i4b) + f1(:n,j)
            end do
!
            do j = 2_i4b, nf
                f2(:n,j) = f2(:n,j-1_i4b) + f2(:n,j)
            end do
!
!           NORMALIZED THE CUMULATIVE PERIODOGRAMS IF REQUIRED.
!
            if ( norm2 ) then
!
                ks(:n) = one/f1(:n,nf)
!
                do j = 1_i4b, nfm1
                    f1(:n,j) = f1(:n,j)*ks(:n)
                end do
!
                ks(:n) = one/f2(:n,nf)
!
                do j = 1_i4b, nfm1
                    f2(:n,j) = f2(:n,j)*ks(:n)
                end do
!
            end if
!
!           COMPUTE THE KOLMOGOROV-SMIRNOV STATISTICS.
!
            ks(:n) = maxval( abs( f1(:n,:nfm1) - f2(:n,:nfm1) ), dim=2 )
!
!           CALCULATE THE SIGNIFICANCE PROBABILITY OF THE OBSERVED
!           KOLMOGOROV-SMIRNOV STATISTICS.
!
            where( ks(:n)>=ks_stat(:n) ) nge(:n) = nge(:n) + 1_i4b
!
        end do
!
#ifdef _OPENMP
    end if
#endif
!
!   COMPUTE THE SIGNIFICANCE LEVELS.
!
    prob(:n) = real( nge(:n)+1_i4b, stnd )*(one/real( nrep2+1_i4b, stnd ))
!
!
! END OF SUBROUTINE spctrm_diff_rm
! ________________________________
!
    end subroutine spctrm_diff_rm
!
! =========================================================================================
!
    subroutine spctrm_diff2_rv( psvec1, psvec2, chi2_stat, prob, nrep, initseed  )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_DIFF2 calculates an approximate critical probability value (e.g. p-value)
!   for testing the hypothesis of a common underlying spectrum for the two estimated spectra (e.g. the
!   arguments PSVEC1 and PSVEC2). This approximate critical probability value is derived from the following
!   CHI2 statistic (e.g. the CHI2_STAT argument) :
!
!       CHI2_STAT  =  (1/nf) [ sum k=1 to nf ]  log( PSVEC1(k) / PSVEC2(k) )**(2)
!
!   where nf = size(PSVEC1) = size(PSVEC2). The distribution of CHI2_STAT under the null 
!   hypothesis of a common spectrum for the two series is approximated by calculating
!   CHI2_STAT for some large number (e.g. the NREP argument) of random interchanges of periodogram ordinates
!   at each frequency for the two estimated spectra (e.g. the arguments PSVEC1(:) and PSVEC2(:)).
!
! Arguments
! _________
!
!
!   PSVEC1     (INPUT) real(stnd), dimension(:)
!              On entry, a real vector containing the Power Spectral Density (PSD)
!              estimates of the first time series.
!
!              All elements in PSVEC1(:) must be greater than zero.
!              size(PSVECN) must be greater or equal to 10.
!
!   PSVEC2     (INPUT) real(stnd), dimension(:)
!              On entry, a real vector containing the Power Spectral Density (PSD)
!              estimates of the second time series.  
!
!              All elements in PSVEC2(:) must be greater than zero.
!
!              PSVEC2 must verify:  size(PSVEC2) = size(PSVEC1) .
!
!   CHI2_STAT  (OUTPUT) real(stnd)
!              On output, the CHI2 statistic.
!
!   PROB       (OUTPUT) real(stnd)
!              On exit, the aproximate critical probability value (e.g. p-value) computed
!              under the hypothesis that the two "true" underlying spectra are the same.
!
!              PROB is calculated as the probability of obtaining a value greater or equal
!              to CHI2_STAT under the hypothesis of a common spectrum for the 
!              two series.
!
!   NREP       (INPUT, OPTIONAL) integer(i4b)
!              On entry, when argument NREP is present, NREP specifies the number of random interchanges
!              of periodogram ordinates at each frequency in order to approximate the randomization 
!              distribution of CHI2_STAT under the null hypothesis.
!
!              The default is 99.
!
!   INITSEED   (INPUT, OPTIONAL) logical(lgl)
!              On entry, if INITSEED=true, a call to RANDOM_SEED_() without arguments
!              is done in the subroutine, in order to initiates a non-repeatable reset
!              of the seed used by the STATPACK random generator.
!
!              The default is INITSEED=false.
!
!
! Further Details
! _______________
!
!   This statistical randomization test relies on the assumptions that the different spectral ordinates
!   have the same sampling distribution and are independent of each other. This means, in
!   particular, that the spectral ordinates corresponding to the zero and Nyquist frequencies
!   must be excluded from the PSVEC1 and PSVEC2 vectors before calling SPCTRM_DIFF2 and
!   that the two estimated spectra have not been obtained by smoothing the periodograms
!   in the frequency domain.
!
!   This randomization test could only be used to compare two periodograms or two spectral
!   estimates computed as the the average of, say, r periodograms for each time series.
!
!   Finally, none of the spectral estimates must be zero.
!
!   For more details see:
!
!   (1) Diggle, P.J., and Fisher, N.I., 1991:
!               Nonparametric comparison of cumulative periodograms,
!               Applied Statistics, Vol. 40, No 3, pp. 423-434.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert_eq, merror
    use Reals_Constants,   only : zero, half, one
    use Char_Constants,    only : stat_error8, tseries_error71, tseries_error73
    use Random,            only : random_seed_, random_number_
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),  dimension(:) :: psvec1, psvec2
    real(stnd), intent(out)               :: chi2_stat, prob
!
    integer(i4b), intent(in), optional :: nrep
!
    logical(lgl), intent(in),  optional :: initseed
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: nf, nrep2, nge, i
!
    real(stnd)                            :: chi2
    real(stnd),  dimension(size(psvec1))  :: lnratio
!
#ifdef _OPENMP
    integer(i4b) :: nge2, j
!
    logical      :: test_par
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_diff2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    nf =  assert_eq( int(size(psvec1),i4b) ,          &
                     int(size(psvec2),i4b) ,          &
                     name_proc )
!
    if ( nf<10_i4b  )    &
    call merror( name_proc//tseries_error71 )
!
    if ( any( psvec1(:nf)<=zero .or. psvec2(:nf)<=zero )  )       &
    call merror( name_proc//tseries_error73 )
!
    if ( present(initseed) ) then
!
        if ( initseed ) then
!
!           INITIALIZE THE RANDOM GENERATOR IF REQUIRED.
!
!$OMP SINGLE
            call random_seed_()
!$OMP END SINGLE
!
        end if
!
    end if
!
    nrep2 = 99_i4b
!
    if ( present(nrep) ) then
        if ( nrep<=0_i4b )        &
        call merror( name_proc//stat_error8 )
        nrep2 = nrep
    end if
!
!   COMPUTE THE SQUARE OF THE LOG RATIO OF THE TWO SPECTRAL ESTIMATES.
!
    lnratio(:nf) = ( log( psvec1(:nf)/psvec2(:nf) ) )**2
!
!   COMPUTE THE CHI2 STATISTIC AND THE ASSOCIATED PROBABILITY.
!
    chi2_stat = sum( lnratio(:nf) )
!
    nge = 0_i4b
!
#ifdef _OPENMP
!
    i = omp_get_num_procs()
    j = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() ) .and. i>1 .and. j>1
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(lnratio,chi2,nge2,i)
!
        nge2 = 0_i4b
!
!$OMP DO SCHEDULE(STATIC) 
!
        do i=1_i4b, nrep2
!
!           GENERATE A RANDOM NUMBER UNIFORM SEQUENCE.
!
            call random_number_( lnratio(:nf) )
!
!           RANDOM INTERCHANGES OF PERIODOGRAM ORDINATES.
!
            where( lnratio(:nf)<=half )
                lnratio(:nf) = ( log( psvec1(:nf)/psvec2(:nf) ) )**2
            elsewhere
                lnratio(:nf) = ( log( psvec2(:nf)/psvec1(:nf) ) )**2
            end where
!
!           COMPUTE THE CHI2 STATISTIC.
!
            chi2 = sum( lnratio(:nf) )
!
            if ( chi2>=chi2_stat ) nge2 = nge2 + 1_i4b
!
        end do
!
!$OMP END DO NOWAIT
!
!$OMP ATOMIC
        nge = nge + nge2
!
!$OMP END PARALLEL
!
    else
!
#endif
!
        do i=1_i4b, nrep2
!
!           GENERATE A RANDOM NUMBER SEQUENCE.
!
            call random_number_( lnratio(:nf) )
!
!           RANDOM INTERCHANGES OF PERIODOGRAM ORDINATES.
!
            where( lnratio(:nf)<=half )
                lnratio(:nf) = ( log( psvec1(:nf)/psvec2(:nf) ) )**2
            elsewhere
                lnratio(:nf) = ( log( psvec2(:nf)/psvec1(:nf) ) )**2
            end where
!
!           COMPUTE THE CHI2 STATISTIC.
!
            chi2 = sum( lnratio(:nf) )
!
            if ( chi2>=chi2_stat ) nge = nge + 1_i4b
!
        end do
!
#ifdef _OPENMP
    end if
#endif
!
    chi2_stat = ( one/real( nf, stnd ) )*chi2_stat
!
!   COMPUTE THE SIGNIFICANCE LEVEL.
!
    prob = real( nge+1_i4b, stnd )/real( nrep2+1_i4b, stnd )
!
!
! END OF SUBROUTINE spctrm_diff2_rv
! _________________________________
!
    end subroutine spctrm_diff2_rv
!
! =========================================================================================
!
    subroutine spctrm_diff2_rm( psmat1, psmat2, chi2_stat, prob, nrep, initseed  )
!
! Purpose
! _______
!
!   Subroutine SPCTRM_DIFF2 calculates approximate critical probability values (e.g. p-value)
!   for testing the hypothesis of a common spectrum for two estimated multi-channel spectra (e.g. the
!   arguments PSMAT1 and PSMAT2). These approximate critical probability values are derived from
!   the following CHI2 statistics (e.g. the CHI2_STAT(:) argument) :
!
!       CHI2_STAT(:n) =  ( 1/nf ) [ sum k=1 to nf ]  log( PSMAT1(:n,k) / PSMAT2(:n,k) )**(2)
!
!   where n = size(PSMAT1,1) = size(PSMAT2,1) = size(CHI2_STAT) and  nf = size(PSMAT2,2) = size(PSMAT2,2).
!
!   The distribution of CHI2_STAT under the null hypothesis of a common spectrum for the spectra of the
!   two multi-channel series is approximated by calculating CHI2_STAT for some large number (e.g. the NREP
!   argument) of random interchanges of periodogram ordinates at each frequency for the two 
!   estimated multi-channel spectra (e.g. the arguments PSMAT1(:,:) and PSMAT2(:,:)).
!
! Arguments
! _________
!
!
!   PSMAT1     (INPUT) real(stnd), dimension(:,:)
!              On entry, a real matrix containing the Power Spectral Density (PSD)
!              estimates of the first multi-channel time series.  
!              Each row of the real matrix PSMAT1 contains the estimated spectrum
!              of the corresponding "row" of the first multi-channel times series.  
!
!              All elements in PSMAT1(:,:) must be greater than zero.
!              size(PSMATN,2) must be greater or equal to 10.
!
!   PSMAT2     (INPUT) real(stnd), dimension(:,:)
!              On entry, a real matrix containing the Power Spectral Density (PSD)
!              estimates of the second multi-channel time series.  
!              Each row of the real matrix PSMAT2 contains the estimated spectrum
!              of the corresponding "row" of the second multi-channel times series.  
!
!              All elements in PSMAT2(:,:) must be greater than zero.
!
!              PSMAT2 must verify:
!
!              - size(PSMAT2,1) = size(PSMAT1,1) ,
!              - size(PSMAT2,2) = size(PSMAT2,2) .
!
!   CHI2_STAT  (OUTPUT) real(stnd), dimension(:)
!              On output, the CHI2 statistics computed from the multi-channel 
!              times series .
!
!              CHI2_STAT must verify:  size(CHI2_STAT) = size(PSMAT1,1) .
!
!   PROB       (OUTPUT) real(stnd), dimension(:)
!              On exit, the aproximate critical probability values (e.g. p-values) computed
!              under the hypothesis that the two "true" underlying multi-channel spectra are
!              the same.
!
!              PROB is calculated as the probability of obtaining a value greater or equal
!              to CHI2_STAT under the hypothesis of a common spectrum for the two
!              multi-channel series.
!
!              PROB must verify:  size(PROB) = size(PSMAT1,1) .
!
!   NREP       (INPUT, OPTIONAL) integer(i4b)
!              On entry, when argument NREP is present, NREP specifies the number of random interchanges
!              of periodogram ordinates at each frequency in order to approximate the randomization 
!              distribution of CHI2_STAT under the null hypothesis.
!
!              The default is 99.
!
!   INITSEED   (INPUT, OPTIONAL) logical(lgl)
!              On entry, if INITSEED=true, a call to RANDOM_SEED_() without arguments
!              is done in the subroutine, in order to initiate a non-repeatable reset
!              of the seed used by the STATPACK random generator.
!
!              The default is INITSEED=false.
!
!
! Further Details
! _______________
!
!   This statistical randomization test relies on the assumptions that the different spectral ordinates
!   have the same sampling distribution and are independent of each other. This means, in
!   particular, that the spectral ordinates corresponding to the zero and Nyquist frequencies
!   must be excluded from the PSMAT1 and PSMAT2 matrices before calling SPCTRM_DIFF2 and
!   that the two estimated multi-channel spectra have not been obtained by smoothing the
!   periodograms in the frequency domain.
!
!   This randomization test could only be used to compare two periodograms or two spectral
!   estimates computed as the the average of, say, r periodograms for each time series.
!
!   Finally, none of the spectral estimates must be zero.
!
!   For more details see:
!
!   (1) Diggle, P.J., and Fisher, N.I., 1991:
!               Nonparametric comparison of cumulative periodograms,
!               Applied Statistics, Vol. 40, No 3, pp. 423-434.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert_eq, merror
    use Reals_Constants,   only : zero, half, one
    use Char_Constants,    only : stat_error8, tseries_error72, tseries_error73
    use Random,            only : random_seed_, random_number_
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(in),   dimension(:,:) :: psmat1, psmat2
    real(stnd), intent(out),  dimension(:)   :: chi2_stat, prob
!
    integer(i4b), intent(in), optional :: nrep
!
    logical(lgl), intent(in),  optional :: initseed
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b)                             :: n, nf, nrep2, i
    integer(i4b),  dimension(size(psmat1,1)) :: nge
!
    real(stnd),  dimension(size(psmat1,1))                 :: chi2
    real(stnd),  dimension(size(psmat1,1),size(psmat1,2))  :: lnratio
!
#ifdef _OPENMP
    integer(i4b)                             :: j
    integer(i4b),  dimension(size(psmat1,1)) :: nge2
!
    logical      :: test_par
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='spctrm_diff2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n =  assert_eq( int(size(psmat1,1),i4b) ,          &
                    int(size(psmat2,1),i4b) ,          &
                    int(size(chi2_stat),i4b) ,         &
                    int(size(prob),i4b) ,              &
                    name_proc )
!
    if ( n<1_i4b ) return
!
    nf =  assert_eq( int(size(psmat1,2),i4b) ,          &
                     int(size(psmat2,2),i4b) ,          &
                     name_proc )
!
    if ( nf<10_i4b  )    &
    call merror( name_proc//tseries_error72 )
!
    if ( any( psmat1(:n,:nf)<=zero .or. psmat2(:n,:nf)<=zero )  )       &
    call merror( name_proc//tseries_error73 )
!
    if ( present(initseed) ) then
!
        if ( initseed ) then
!
!           INITIALIZE THE RANDOM GENERATOR IF REQUIRED.
!
!$OMP SINGLE
            call random_seed_()
!$OMP END SINGLE
!
        end if
!
    end if
!
    nrep2 = 99_i4b
!
    if ( present(nrep) ) then
!
        if ( nrep<=0_i4b )        &
        call merror( name_proc//stat_error8 )
!
        nrep2 = nrep
!
    end if
!
!   COMPUTE THE SQUARE OF THE LOG RATIO OF THE TWO SPECTRAL ESTIMATES.
!
    lnratio(:n,:nf) = ( log( psmat1(:n,:nf)/psmat2(:n,:nf) ) )**2
!
!   COMPUTE THE CHI2 STATISTIC AND THE ASSOCIATED PROBABILITY.
!
    chi2_stat(:n) = sum( lnratio(:n,:nf), dim=2 )
!
    nge(:n) = 0_i4b
!
#ifdef _OPENMP
!
    i = omp_get_num_procs()
    j = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() ) .and. i>1 .and. j>1
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(lnratio,chi2,nge2,i)
!
        nge2(:n) = 0_i4b
!
!$OMP DO SCHEDULE(STATIC) 
!
        do i=1_i4b, nrep2
!
!           GENERATE A RANDOM NUMBER UNIFORM SEQUENCE.
!
            call random_number_( lnratio(:n,:nf) )
!
!           RANDOM INTERCHANGES OF PERIODOGRAM ORDINATES.
!
            where( lnratio(:n,:nf)<=half )
                lnratio(:n,:nf) = ( log( psmat1(:n,:nf)/psmat2(:n,:nf) ) )**2
            elsewhere
                lnratio(:n,:nf) = ( log( psmat2(:n,:nf)/psmat1(:n,:nf) ) )**2
            end where
!
!           COMPUTE THE CHI2 STATISTICS.
!
            chi2(:n) = sum( lnratio(:n,:nf), dim=2 )
!
!           CALCULATE THE SIGNIFICANCE PROBABILITY OF THE OBSERVED
!           CHI2 STATISTICS.
!
            where( chi2(:n)>=chi2_stat(:n) ) nge2(:n) = nge2(:n) + 1_i4b
!
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL  (diff2)
!
        nge(:n) = nge(:n) + nge2(:n)
!
!$OMP END CRITICAL  (diff2)
!
!$OMP END PARALLEL
!
    else
!
#endif
!
        do i=1_i4b, nrep2
!
!           GENERATE A RANDOM NUMBER SEQUENCE.
!
            call random_number_( lnratio(:n,:nf) )
!
!           RANDOM INTERCHANGES OF PERIODOGRAM ORDINATES.
!
            where( lnratio(:n,:nf)<=half )
                lnratio(:n,:nf) = ( log( psmat1(:n,:nf)/psmat2(:n,:nf) ) )**2
            elsewhere
                lnratio(:n,:nf) = ( log( psmat2(:n,:nf)/psmat1(:n,:nf) ) )**2
            end where
!
!           COMPUTE THE CHI2 STATISTICS.
!
            chi2(:n) = sum( lnratio(:n,:nf), dim=2 )
!
!           CALCULATE THE SIGNIFICANCE PROBABILITY OF THE OBSERVED
!           CHI2 STATISTICS.
!
            where( chi2(:n)>=chi2_stat(:n) ) nge(:n) = nge(:n) + 1_i4b
!
        end do
!
#ifdef _OPENMP
    end if
#endif
!
    chi2_stat(:n) = ( one/real( nf, stnd ) )*chi2_stat(:n)
!
!   COMPUTE THE SIGNIFICANCE LEVELS.
!
    prob(:n) = real( nge(:n)+1_i4b, stnd )*(one/real( nrep2+1_i4b, stnd ))
!
!
! END OF SUBROUTINE spctrm_diff2_rm
! _________________________________
!
    end subroutine spctrm_diff2_rm
!
! =========================================================================================
!
    subroutine power_spctrm_rv( vec, psvec, freq, fftvec, edof, bandwidth, conlwr, conupr, &
                                initfft, normpsd, nsmooth, trend, win, taperp, probtest    )
!
! Purpose
! _______
!
!   Subroutine POWER_SPCTRM computes a Fast Fourier Transform (FFT) estimate
!   of the power spectrum of a real time series, VEC. The real valued sequence
!   VEC must be of even length.
!
!   The Power Spectral Density (PSD) estimates are returned in units which are
!   the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the real time series for which the power spectrum
!                 must be estimated.
!                 If WIN/=2 or TREND=1, 2 or 3,  VEC is used
!                 as workspace and is transformed.
!
!                 Size(VEC) must be an even (positive) integer greater or equal to 4.
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = size(VEC)/2 + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/size(VEC)
!                 for i=1,2, ... , (size(VEC)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = size(VEC)/2 + 1 .
!
!   FFTVEC        (OUTPUT, OPTIONAL) complex(stnd), dimension(:)
!                 On exit, a complex vector of length (size(VEC)/2)+1 containing
!                 the Fast Fourier Transform of the product of the (detrended, e.g. 
!                 the TREND argument) real time series VEC with the choosen window 
!                 function (e.g. The WIN argument).
!
!                 FFTVEC must verify:  size(FFTVEC) = size(VEC)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!                 EDOF must verify:  size(EDOF) = size(VEC)/2 + 1 .
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the bandwidth of the power spectrum estimates.
!
!                 BANDWIDTH must verify:  size(BANDWIDTH) = size(VEC)/2 + 1 .
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100%
!                 confidence limit factors, respectively. Multiply the PSD estimates (e.g. the
!                 PSVEC(:) argument) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!                 CONLWR must verify:  size(CONLWR) = size(VEC)/2 + 1 .
!
!                 CONUPR must verify:  size(CONUPR) = size(VEC)/2 + 1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine POWER_SPCTRM in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine POWER_SPCTRM. This call to INITFFT must have the following form: 
!
!                       call init_fft( size(VEC)/2 )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   POWER_SPCTRM and a call to END_FFT is also done before leaving
!                   subroutine POWER_SPCTRM.
!
!                 The default is INITFFT=true .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the PSD estimates are normalized in such
!                   a way that the total area under the power spectrum is equal to the variance of
!                   the time series VEC.
!                 - NORMPSD = false, the sum of the PSD estimates
!                   (e.g. sum( PSVEC(2:) ) is equal to the variance of the time series.
!
!                 The default is NORMPSD=true .
!
!   NSMOOTH       (INPUT, OPTIONAL) integer(i4b)
!                 if NSMOOTH is used, the PSD estimates are computed by smoothing 
!                 the periodogram with Daniell weights (e.g. a simple moving average).
!
!                 On entry, NSMOOTH gives the length of the Daniell filter to be applied.
!
!                 Setting NSMOOTH=0 on entry is equivalent to omit the optional argument
!                 NSMOOTH. Otherwise, NSMOOTH must be odd, greater than 2 and less or
!                 equal to size(VEC)/2+1 .
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The mean of the time series is removed before computing the spectrum
!                 - TREND=+2 The drift from the time series is removed before computing the spectrum
!                   by using the formula: drift = (VEC(size(VEC)) - VEC(1))/(size(VEC) - 1)
!                 - TREND=+3 The least-squares line from the time series is removed before
!                   computing the spectrum.
!
!                 For other values of TREND nothing is done before estimating the power spectrum.
!
!                 The default is TREND=1, e.g. the mean of the time series is removed before the
!                 computations.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power spectrum. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the selected
!   data window (e.g. WIN=1,2,3,4,5,6) is applied to the time series and the PSD estimates are
!   computed by the FFT of this transformed time series. Optionally, theses PSD estimates
!   may then be smoothed in the frequency domain by a Daniell filter (e.g. if
!   NSMOOTH is used).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!      
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, arth
    use Reals_Constants,   only : zero, one, two, pi
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error14, tseries_error15,    &
                                  tseries_error59, tseries_error61, tseries_error65
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:) :: vec
    real(stnd), intent(out),    dimension(:) :: psvec
!
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out), dimension(:), optional :: freq, edof, conlwr, conupr, bandwidth
!
    complex(stnd), intent(out), dimension(:), optional :: fftvec
!
    integer(i4b),  intent(in),               optional :: trend, win, nsmooth
!
    logical(lgl),  intent(in), optional :: normpsd, initfft
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, nd2, nf, win2, trend2
    integer      :: iok
!
    real(stnd)                            :: c1, c2, sumw
    real(stnd), dimension(size(vec))      :: wk
    real(stnd), dimension(:), allocatable :: edof2
!
    complex(stnd), dimension(size(vec)/2+1)  :: cwk
!
    logical(lgl)  :: normpsd2, initfft2, smooth, out_freq, out_fftvec, out_dof
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='power_spctrm'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n   = size( vec )
!
    if ( n<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    nd2 = n/2_i4b
!
    if ( n/=2*nd2 )   &
    call merror( name_proc//tseries_error14 )
!
    nf  = nd2 + 1_i4b
!
    call assert( logical(int(size(psvec),i4b)==nf,lgl),     &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_fftvec = present( fftvec )
!
    if ( out_fftvec ) then
        call assert( logical(int(size(fftvec),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth )        &
              .or. present( conlwr ) .or. present( conupr ) 
!
    if ( out_dof ) then
!
        if ( present(probtest) ) then
            if ( zero>=probtest .or. probtest>=one ) then
                call merror( name_proc//tseries_error59 )
            end if
        end if
!
        if ( present( edof ) ) then
            call assert( logical(int(size(edof),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( bandwidth ) ) then
            call assert( logical(int(size(bandwidth),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conlwr ) ) then
            call assert( logical(int(size(conlwr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conupr ) ) then
            call assert( logical(int(size(conupr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
    end if
!
    if ( present(trend) ) then
        trend2 = trend
    else
        trend2 = 1_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
    else
        win2 = 3_i4b
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( nd2 )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    smooth = false
!
    if ( present(nsmooth) ) then
!
        if ( nsmooth/=0_i4b ) then
!
!           CHECK THE INPUT VALUE FOR THE LENGTH OF THE DANIELL FILTER.
!
            if ( nsmooth<3_i4b .or. nsmooth>nf )      &
            call merror( name_proc//tseries_error65 )
!
            if ( (nsmooth/2_i4b)*2_i4b==nsmooth )   &
            call merror( name_proc//tseries_error61 )
!
            smooth = true
!
        end if
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trend2>=1_i4b .and. trend2<=3_i4b ) then
        call detrend_rv( vec(:n), trend2 )
    end if
!
!   TAPER THE TIME SERIES IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:n)  = data_window( n, win2, taperp=taperp )
!
!       APPLY DATA WINDOW.
!
        vec(:n) = vec(:n)*wk(:n)
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:n), wk(:n) )
!
    else    
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:n)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( n, stnd )
!
    end if
!
!   COMPUTE PSD ESTIMATES
!
    call real_fft( vec(:n), cwk(:nf), true)
!
    psvec(:nd2) = real(cwk(:nd2),stnd)**2 + aimag(cwk(:nd2))**2
    psvec(nf)   = real(cwk(nf),stnd)**2
!
!   NORMALIZE THE PSD ESTIMATES.
!
    c1 = one/(sumw*real( n, stnd))
!
    if ( normpsd2 ) then
!
!       CHATFIELD DEFINITION OF THE SPECTRAL DENSITY FUNCTION.
!
        c2 = real( n, stnd)/pi
        psvec(1_i4b:nf) = (c1*c2)*psvec(1_i4b:nf)
!
    else
!
!       POLLOCK DEFINITION OF THE PERIODOGRAM. IF THE ESTIMATES
!       ARE DIVIDED BY 4*PI, THIS IS THE POLLOCK DEFINITION OF
!       THE SPECTRAL DENSITY FUNCTION.
!
        psvec(1_i4b)     = c1*psvec(1_i4b)
        psvec(2_i4b:nd2) = (two*c1)*psvec(2_i4b:nd2)
        psvec(nf)        = c1*psvec(nf)
!
    end if
!
!   SMOOTH THE PSD ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
        call dan_filter_rv( psvec(:nf), nsmooth, sym=one )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( n, stnd )
        freq(:nf) = arth(zero, c1, nf)
    end if
!
!   OUTPUT FFT OF THE DATA IF REQUIRED.
!
    if ( out_fftvec ) then
        fftvec(:nf) = cwk(:nf)
    end if
!
!   OUTPUT DEGREES OF FREEDOM, BANDWIDTH OR CONFIDENCE LIMIT FACTORS IF REQUIRED.
!
    if ( out_dof ) then
!
        allocate( edof2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
        edof2(:nf) = estim_dof2( wk(:n), 0_i4b, nsmooth=nsmooth )
!
        if ( present(conlwr) .or. present(conupr) ) then
            call comp_conflim_rv( edof2(:nf), probtest=probtest, conlwr=conlwr, conupr=conupr )
        end if
!
        if ( present(edof) ) then
            edof(:nf) = edof2(:nf)
        end if
!
        if ( present(bandwidth) ) then
            bandwidth(:nf) = edof2(:nf)/( two*real( n, stnd ) )
        end if
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( edof2 )
!
    end if
!
!   DEALLOCATE WORK ARRAYS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE power_spctrm_rv
! _________________________________
!
    end subroutine power_spctrm_rv
!
! =========================================================================================
!
    subroutine power_spctrm_rm( mat, psmat, freq, fftmat, edof, bandwidth, conlwr, conupr, &
                                initfft, normpsd, nsmooth, trend, win, taperp, probtest    )
!
! Purpose
! _______
!
!   Subroutine POWER_SPCTRM computes a Fast Fourier Transform (FFT) estimate
!   of the power spectra of the rows of the real matrix, MAT. size(MAT,2) must
!   be of even length.
!
!   The Power Spectral Density (PSD) estimates are returned in units which are
!   the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   MAT           (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                 On entry, the real time series for which power spectra
!                 must be estimated. Each row of MAT is a real time series.
!
!                 If WIN/=2 or TREND=1, 2 or 3,  MAT is used
!                 as workspace and is transformed.
!
!                 Size(MAT,2) must be an even (positive) integer greater or equal to 4.
!
!   PSMAT         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix containing the Power Spectral Density
!                 (PSD) estimates for each row of the real matrix MAT.  
!
!                 The shape of PSMAT must verify:
!
!                 - size(PSMAT,1) = size(MAT,1) ;
!                 - size(PSMAT,2) = size(MAT,2)/2 + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(MAT,2)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/size(VEC)
!                 for i=1,2, ... , (size(MAT,2)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = size(MAT,2)/2 + 1 .
!
!   FFTMAT        (OUTPUT, OPTIONAL) complex(stnd), dimension(:,:)
!                 On exit, a complex matrix containing the Fast Fourier Transform 
!                 of the product of the (detrended, e.g. the TREND argument) 
!                 real time series in each row of MAT with the choosen window 
!                 function (e.g. The WIN argument).
!
!                 The shape of FFTMAT must verify:
!
!                 - size(FFTMAT,1) = size(MAT,1) ;
!                 - size(FFTMAT,2) = size(MAT,2)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!                 EDOF must verify:  size(EDOF) = size(MAT,2)/2 + 1 .
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the bandwidth of the power spectrum estimates.
!
!                 BANDWIDTH must verify:  size(BANDWIDTH) = size(MAT,2)/2 + 1 .
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the
!                 PSMAT(:,:) argument) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!                 CONLWR must verify:  size(CONLWR) = size(MAT,2)/2 + 1 .
!
!                 CONUPR must verify:  size(CONUPR) = size(MAT,2)/2 + 1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine POWER_SPCTRM in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine POWER_SPCTRM. This call to INITFFT must have the following form: 
!
!                         call init_fft( (/ size(MAT,1), size(MAT,2)/2 /), dim=2_i4b )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   POWER_SPCTRM and a call to END_FFT is also done before leaving
!                   subroutine POWER_SPCTRM
!
!                 The default is INITFFT=true .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD is set to true, the PSD estimates are normalized in such
!                   a way that the total area under the power spectrum is equal to the variance of
!                   the time series MAT.
!                 - NORMPSD = false, the sum of the PSD estimates
!                   for each row of MAT (e.g. sum( PSMAT(:,2:), dim=2 ) is equal to the variance of
!                   the corresponding time series.
!
!                 The default is NORMPSD=true .
!
!   NSMOOTH       (INPUT, OPTIONAL) integer(i4b)
!                 if NSMOOTH is used, the PSD estimates are computed by smoothing 
!                 the periodogram with Daniell weights (e.g. a simple moving average).
!
!                 On entry, NSMOOTH gives the length of the Daniell filter to be applied.
!
!                 Setting NSMOOTH=0 on entry is equivalent to omit the optional argument
!                 NSMOOTH. Otherwise, NSMOOTH must be odd, greater than 2 and less or
!                 equal to size(MAT,2)/2+1 .
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The means of the time series are removed before computing the spectra
!                 - TREND=+2 The drifts from the time series are removed before computing the spectra
!                   by using the formula: drift(i) = (MAT(i,size(MAT,2)) - MAT(i,1))/(size(MAT,2) - 1)
!                 - TREND=+3 The least-squares lines from the time series are removed before
!                   computing the spectra.
!
!                 For other values of TREND nothing is done before estimating the power spectra.
!
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power spectrum. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!                   WIN=+1 The Bartlett window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the selected
!   data window (e.g. WIN=1,2,3,4,5,6) is applied to the time series and the PSD estimates are
!   computed by the FFT of these transformed time series. Optionally, theses PSD estimates
!   may then be smoothed in the frequency domain by modified Daniell filters (e.g. if 
!   SMOOTH_PARAM is used).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, arth
    use Reals_Constants,   only : zero, one, two, pi
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error15, tseries_error16,    &
                                  tseries_error59, tseries_error61, tseries_error66
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:,:) :: mat
!
    real(stnd), intent(out), dimension(:,:)         :: psmat
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out), dimension(:), optional :: freq, edof, bandwidth, conlwr, conupr
!
    complex(stnd), intent(out), dimension(:,:), optional :: fftmat
!
    integer(i4b),  intent(in),               optional :: trend, win, nsmooth
!
    logical(lgl),  intent(in), optional :: normpsd, initfft
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, n, nd2, nf, k, win2, trend2
    integer      :: iok
!
    real(stnd)                               :: c1, c2, sumw
    real(stnd), dimension(size(mat,2))       :: wk
    real(stnd), dimension(:),    allocatable :: edof2
!
    complex(stnd), dimension(size(mat,1),size(mat,2)/2+1)  :: cwk
!
    logical(lgl)  :: normpsd2, initfft2, smooth, out_freq, out_fftmat, out_dof
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='power_spctrm'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m   = size( mat, 1 )
    n   = size( mat, 2 )
!
    if ( m<=0_i4b ) return
!
    if ( n<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    nd2 = n/2_i4b
!
    if ( n/=2*nd2 )   &
    call merror( name_proc//tseries_error16 )
!
    nf  = nd2 + 1_i4b
!
    call assert( logical(int(size(psmat,1),i4b)==m,lgl),     &
                 logical(int(size(psmat,2),i4b)==nf,lgl),    &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_fftmat = present( fftmat )
!
    if ( out_fftmat ) then
        call assert( logical(int(size(fftmat,1),i4b)==m,lgl),     &
                     logical(int(size(fftmat,2),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth )        &
              .or. present( conlwr ) .or. present( conupr ) 
!
    if ( out_dof ) then
!
        if ( present(probtest) ) then
            if ( zero>=probtest .or. probtest>=one ) then
                call merror( name_proc//tseries_error59 )
            end if
        end if
!
        if ( present( edof ) ) then
            call assert( logical(int(size(edof),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( bandwidth ) ) then
            call assert( logical(int(size(bandwidth),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conlwr ) ) then
            call assert( logical(int(size(conlwr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conupr ) ) then
            call assert( logical(int(size(conupr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
    end if
!
    if ( present(trend) ) then
        trend2 = trend
    else
        trend2 = 1_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
    else
        win2 = 3_i4b
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( (/ m, nd2/), dim=2_i4b )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    smooth = false
!
    if ( present(nsmooth) ) then
!
        if ( nsmooth/=0_i4b ) then
!
!           CHECK THE INPUT VALUE FOR THE LENGTH OF THE DANIELL FILTER.
!
            if ( nsmooth<3_i4b .or. nsmooth>nf )      &
            call merror( name_proc//tseries_error66 )
!
            if ( (nsmooth/2_i4b)*2_i4b==nsmooth )   &
            call merror( name_proc//tseries_error61 )
!
            smooth = true
!
        end if
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trend2>=1_i4b .and. trend2<=3_i4b ) then
        call detrend_rm( mat(:m,:n), trend2 )
    end if
!
!   TAPER THE TIME SERIES IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:n)  = data_window( n, win2, taperp=taperp )
!
!       APPLY DATA WINDOW.
!
        do k = 1_i4b, n
            mat(:m,k) = wk(k)*mat(:m,k)
        end do
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:n), wk(:n) )
!
    else    
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:n)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( n, stnd )
!
    end if
!
!   COMPUTE THE PSD ESTIMATES.
!
    call real_fft( mat(:m,:n), cwk(:m,:nf), true)
!
    psmat(:m,:nd2) = real(cwk(:m,:nd2),stnd)**2 + aimag(cwk(:m,:nd2))**2
    psmat(:m,nf)   = real(cwk(:m,nf),stnd)**2
!
!   NORMALIZE THE PSD ESTIMATES.
!
    c1 = one/(sumw*real( n, stnd))
!
    if ( normpsd2 ) then
!
!       CHATFIELD DEFINITION OF THE SPECTRAL DENSITY FUNCTION.
!
        c2 = real( n, stnd)/pi
        psmat(:m,:nf) = (c1*c2)*psmat(:m,:nf)
!
    else
!
!       POLLOCK DEFINITION OF THE PERIODOGRAM. IF THE ESTIMATES
!       ARE DIVIDED BY 4*PI, THIS IS THE POLLOCK DEFINITION OF
!       THE SPECTRAL DENSITY FUNCTION.
!
        psmat(:m,1_i4b)     = c1*psmat(:m,1_i4b)
        psmat(:m,2_i4b:nd2) = (two*c1)*psmat(:m,2_i4b:nd2)
        psmat(:m,nf)        = c1*psmat(:m,nf)
!
    end if
!
!   SMOOTH THE PSD ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
        call dan_filter_rm( psmat(:m,:nf), nsmooth, sym=one )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( n, stnd )
        freq(:nf) = arth(zero, c1, nf)
    end if
!
!   OUTPUT FFT OF THE DATA IF REQUIRED.
!
    if ( out_fftmat ) then
        fftmat(:m,:nf) = cwk(:m,:nf)
    end if
!
!   OUTPUT DEGREES OF FREEDOM, BANDWIDTH OR CONFIDENCE LIMIT FACTORS IF REQUIRED.
!
    if ( out_dof ) then
!
        allocate( edof2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
        edof2(:nf) = estim_dof2( wk(:n), 0_i4b, nsmooth=nsmooth )
!
        if ( present(conlwr) .or. present(conupr) ) then
            call comp_conflim_rv( edof2(:nf), probtest=probtest, conlwr=conlwr, conupr=conupr )
        end if
!
        if ( present(edof) ) then
            edof(:nf) = edof2(:nf)
        end if
!
        if ( present(bandwidth) ) then
            bandwidth(:nf) = edof2(:nf)/( two*real( n, stnd ) )
        end if
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( edof2 )
!
    end if
!
!   DEALLOCATE WORK ARRAYS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE power_spctrm_rm
! _________________________________
!
    end subroutine power_spctrm_rm
!
! =========================================================================================
!
    subroutine cross_spctrm_rv( vec, vec2, psvec, psvec2, phase, coher, freq, edof, bandwidth, &
                                conlwr, conupr, testcoher, ampli, co_spect, quad_spect,        &
                                prob_coher, initfft, normpsd, nsmooth, trend, win, taperp,     &
                                probtest )
!
! Purpose
! _______
!
!   Subroutine CROSS_SPCTRM computes Fast Fourier Transform (FFT) estimates
!   of the power and cross spectra of two real time series, VEC and VEC2. The real valued
!   sequences VEC and VEC2  must be of even length.
!
!   The Power Spectral Density (PSD) and Cross Spectral Density (CSD) estimates are returned
!   in units which are the square of the data (if NORMPSD=false) or in spectral density units
!   (if NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the first real time series for which the power and cross spectra
!                 must be estimated.
!                 If WIN/=2 or TREND=1, 2 or 3,  VEC is used as workspace and is transformed.
!
!                 Size(VEC) must be an even (positive) integer greater or equal to 4.
!
!   VEC2          (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the second real time series for which the power and cross spectra
!                 must be estimated.
!                 If WIN/=2 or TREND=1, 2 or 3,  VEC2 is used as workspace and is transformed.
!
!                 VEC2 must verify:  size(VEC2) = size(VEC).
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = size(VEC)/2 + 1 .
!
!   PSVEC2        (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC2)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC2.  
!
!                 PSVEC2 must verify:  size(PSVEC2) = size(VEC)/2 + 1 .
!
!   PHASE         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the phase of the cross spectrum, given in fractions
!                 of a circle (e.g. on the closed interval (0,1) ).
!
!                 PHASE must verify:  size(PHASE) = size(VEC)/2 + 1 .
!
!   COHER         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the squared coherency  estimates for all frequencies.
!
!                 COHER must verify:  size(COHER) = size(VEC)/2 + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/size(VEC)
!                 for i=1,2, ... , (size(VEC)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = size(VEC)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!                 EDOF must verify:  size(EDOF) = size(VEC)/2 + 1 .
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the bandwidth of the power spectrum estimates.
!
!                 BANDWIDTH must verify:  size(BANDWIDTH) = size(VEC)/2 + 1 .
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the PSVEC(:) and
!                 PSVEC2(:) arguments) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!                 CONLWR must verify:  size(CONLWR) = size(VEC)/2 + 1 .
!
!                 CONUPR must verify:  size(CONUPR) = size(VEC)/2 + 1 .
!   
!   TESTCOHER     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, this argument specifies the critical value for testing the null
!                 hypothesis that the squared coherency is zero at the PROBTEST * 100% significance
!                 level (e.g. elements of COHER(:) less than TESTCOHER(:) should be regarded as not
!                 significantly different from zero at the PROBTEST * 100% significance level).
!
!                 TESTCOHER must verify:  size(TESTCOHER) = size(VEC)/2 + 1 .
!
!   AMPLI         (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the cross-amplitude spectrum.  
!
!                 AMPLI must verify:  size(AMPLI) = (size(VEC)/2) + 1 .
!
!   CO_SPECT      (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1  containing
!                 the co-spectrum (e.g. the real part of cross-spectrum).  
!
!                 CO_SPECT must verify:  size(CO_SPECT) = (size(VEC)/2) + 1 .
!
!   QUAD_SPECT    (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the quadrature spectrum (e.g. the imaginary part of
!                 cross-spectrum with a minus sign).  
!
!                 QUAD_SPECT must verify:  size(QUAD_SPECT) = (size(VEC)/2) + 1 .
!
!   PROB_COHER    (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the probabilities that the computed sample squared coherencies
!                 came from an ergodic stationary bivariate process with (corresponding)
!                 squared coherencies equal to zero.
!
!                 PROB_COHER  must verify:  size(PROB_COHER) = (size(VEC)/2) + 1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine CROSS_SPCTRM in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine CROSS_SPCTRM. This call to INITFFT must have the following form: 
!
!                        call init_fft( size(VEC)/2 )
!
!                 - INITFFT is set to true, the call to INIT_FFT is done inside subroutine
!                   CROSS_SPCTRM and a call to END_FFT is also done before leaving
!                   subroutine CROSS_SPCTRM.
!
!                 The default is INITFFT=true .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the power and cross spectra estimates are normalized
!                   in such a way that the total area under the power spectrum is equal to the variance
!                   of the time series VEC and VEC2.
!                 - NORMPSD = false, the sum of the PSD estimates
!                   (e.g. sum(PSVEC(2:)) and sum(PSVEC2(2:)) ) is equal to the variance of the corresponding
!                   time series.
!
!                 The default is NORMPSD=true .
!
!   NSMOOTH       (INPUT, OPTIONAL) integer(i4b)
!                 If NSMOOTH is used, the PSD and CSD estimates are computed by smoothing 
!                 the periodogram with Daniell weights (e.g. a simple moving average).
!
!                 On entry, NSMOOTH gives the length of the Daniell filter to be applied.
!
!                 Setting NSMOOTH=0 on entry is equivalent to omit the optional argument
!                 NSMOOTH. Otherwise, NSMOOTH must be odd, greater than 2 and less or
!                 equal to size(VEC)/2+1 .
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 -  TREND=+1 The mean of the two time series is removed before computing 
!                    the power and cross spectra.
!                 -  TREND=+2 The drift from the two time series is removed before computing
!                    the power and cross spectra.
!                 -  TREND=+3 The least-squares line from the two time series is removed before
!                    computing the power and cross spectra.
!
!                 For other values of TREND nothing is done before estimating the power and cross spectra.
!
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power and cross spectra. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ) and the critical value for
!                 testing the null hypothesis that the squared coherency is zero (e.g.
!                 the TESTCOHER optional argument).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the selected
!   data window (e.g. WIN=1,2,3,4,5,6) is applied to the time series and the PSD and CSD estimates are
!   computed by the FFT of these transformed time series. Optionally, theses PSD and CSD estimates
!   may then be smoothed in the frequency domain by modified Daniell filters (e.g. if argument 
!   NSMOOTH is used).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, assert_eq, arth
    use Reals_Constants,   only : zero, one, two, pi, twopi
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error14, tseries_error15,    &
                                  tseries_error59, tseries_error61, tseries_error65
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : probf2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:) :: vec, vec2
    real(stnd), intent(out),    dimension(:) :: psvec, psvec2, phase, coher
!
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out), dimension(:), optional :: freq, ampli, co_spect, quad_spect, prob_coher, &
                                                       edof, bandwidth, conlwr, conupr, testcoher
!
    integer(i4b),  intent(in),               optional :: trend, win, nsmooth
!
    logical(lgl),  intent(in), optional :: normpsd, initfft
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, nd2, nf, k1, k2, win2, trend2, i
    integer      :: iok
!
    real(stnd)                            :: c1, c2, sumw, df, df2
    real(stnd), dimension(size(vec))      :: wk
    real(stnd), dimension(size(vec)/2+1)  :: rcospect, icospect
    real(stnd), dimension(:), allocatable :: edof2
!
    complex(stnd), dimension(size(vec)/2+1)  :: cwk1, cwk2
!
    logical(lgl)  :: normpsd2, initfft2, smooth, out_freq, out_ampli,      &
                     out_co_spect, out_quad_spect, out_prob_coher, out_dof
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='cross_spctrm'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n =  assert_eq( int(size(vec),i4b) ,           &
                    int(size(vec2),i4b) ,          &
                    name_proc )
!
    if ( n<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    nd2 = n/2_i4b
!
    if ( n/=2*nd2 )   &
    call merror( name_proc//tseries_error14 )
!
    nf  = nd2 + 1_i4b
!
    call assert( logical(int(size(psvec),i4b)==nf,lgl),   &
                 logical(int(size(psvec2),i4b)==nf,lgl),  &
                 logical(int(size(phase),i4b)==nf,lgl),   &
                 logical(int(size(coher),i4b)==nf,lgl),   &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_ampli = present( ampli )
!
    if ( out_ampli ) then
        call assert( logical(int(size(ampli),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_co_spect = present( co_spect )
!
    if ( out_co_spect ) then
        call assert( logical(int(size(co_spect),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_quad_spect = present( quad_spect )
!
    if ( out_quad_spect ) then
        call assert( logical(int(size(quad_spect),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_prob_coher = present( prob_coher )
!
    if ( out_prob_coher ) then
        call assert( logical(int(size(prob_coher),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth ) .or. out_prob_coher   &
              .or. present( conlwr ) .or. present( conupr ) .or. present( testcoher ) 
!
    if ( out_dof ) then
!
        if ( present(probtest) ) then
            if ( zero>=probtest .or. probtest>=one ) then
                call merror( name_proc//tseries_error59 )
            end if
        end if
!
        if ( present( edof ) ) then
            call assert( logical(int(size(edof),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( bandwidth ) ) then
            call assert( logical(int(size(bandwidth),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conlwr ) ) then
            call assert( logical(int(size(conlwr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conupr ) ) then
            call assert( logical(int(size(conupr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( testcoher ) ) then
            call assert( logical(int(size(testcoher),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
    end if
!
    if ( present(trend) ) then
        trend2 = trend
    else
        trend2 = 1_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
!
    else
!
        win2 = 3_i4b
!
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( nd2 )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    smooth = false
!
    if ( present(nsmooth) ) then
!
        if ( nsmooth/=0_i4b ) then
!
!           CHECK THE INPUT VALUE FOR THE LENGTH OF THE DANIELL FILTER.
!
            if ( nsmooth<3_i4b .or. nsmooth>nf )      &
            call merror( name_proc//tseries_error65 )
!
            if ( (nsmooth/2_i4b)*2_i4b==nsmooth )   &
            call merror( name_proc//tseries_error61 )
!
            smooth = true
!
        end if
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trend2>=1_i4b .and. trend2<=3_i4b ) then
!
        call detrend_rv( vec(:n),  trend2 )
        call detrend_rv( vec2(:n), trend2 )
!
    end if
!
!   TAPER THE TIME SERIES IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:n)  = data_window( n, win2, taperp=taperp )
!
!       APPLY DATA WINDOW.
!
        vec(:n) = vec(:n)*wk(:n)
!
        vec2(:n) = vec2(:n)*wk(:n)
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:n), wk(:n) )
!
    else    
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:n)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( n, stnd )
!
    end if
!
!   COMPUTE DEGREES OF FREEDOM IF REQUIRED.
!
    if ( out_dof ) then
!
        allocate( edof2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
        edof2(:nf) = estim_dof2( wk(:n), 0_i4b, nsmooth=nsmooth )
!
    end if
!
!   COMPUTE PSD AND CSD ESTIMATES .
!
!   COMPUTE FFT FOR THE FIRST SERIES.
!
    call real_fft( vec(:n), cwk1(:nf), true)
!
!   COMPUTE FFT FOR THE SECOND SERIES.
!
    call real_fft( vec2(:n), cwk2(:nf), true)
!
!   COMPUTE POWER SPECTRUM ESIMATES.
!
    psvec(:nd2) = real(cwk1(:nd2),stnd)**2 + aimag(cwk1(:nd2))**2
    psvec(nf)   = real(cwk1(nf),stnd)**2
!
    psvec2(:nd2) = real(cwk2(:nd2),stnd)**2 + aimag(cwk2(:nd2))**2
    psvec2(nf)   = real(cwk2(nf),stnd)**2
!
!   COMPUTE CO-SPECTRUM ESIMATES.
!
    rcospect(:nd2)  = real(cwk1(:nd2),stnd)*real(cwk2(:nd2),stnd) +   &
                      aimag(cwk1(:nd2))*aimag(cwk2(:nd2))
    rcospect(nf)    = real(cwk1(nf),stnd)*real(cwk2(nf),stnd)
!
!   COMPUTE QUADRATURE-SPECTRUM ESIMATES.
!
    icospect(:nd2)  = aimag(cwk1(:nd2))*real(cwk2(:nd2),stnd)    -   &
                      real(cwk1(:nd2),stnd)*aimag(cwk2(:nd2))
    icospect(nf)    = zero
!
!   NORMALIZE THE POWER AND CROSS SPECTRA ESTIMATES.
!
    c1 = one/(sumw*real( n, stnd))
!
    if ( normpsd2 ) then
!
!       CHATFIELD DEFINITION OF THE SPECTRAL DENSITY FUNCTION.
!
        c2 = (real( n, stnd)/pi)*c1
        k1 = 1_i4b
        k2 = nf
!
    else
!
!       POLLOCK DEFINITION OF THE PERIODOGRAM. IF THE ESTIMATES
!       ARE DIVIDED BY 4*PI, THIS IS THE POLLOCK DEFINITION OF
!       THE SPECTRAL DENSITY FUNCTION.
!
        c2 = two*c1
        k1 = 2_i4b
        k2 = nd2
!
        psvec(1_i4b)  = c1*psvec(1_i4b)
        psvec2(1_i4b) = c1*psvec2(1_i4b)
!
        rcospect(1_i4b) = c1*rcospect(1_i4b)
        icospect(1_i4b) = c1*icospect(1_i4b)
!
        psvec(nf)  = c1*psvec(nf)
        psvec2(nf) = c1*psvec2(nf)
!
        rcospect(nf)  = c1*rcospect(nf)
        icospect(nf)  = c1*icospect(nf)
!
    end if
!
    psvec(k1:k2)  = c2*psvec(k1:k2)
    psvec2(k1:k2) = c2*psvec2(k1:k2)
!
    rcospect(k1:k2) = c2*rcospect(k1:k2)
    icospect(k1:k2) = c2*icospect(k1:k2)
!
!   SMOOTH THE POWER AND CROSS SPECTRA ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
!       SMOOTH POWER SPECTRUM OF THE FIRST SERIES.
!
        call dan_filter_rv( psvec(:nf), nsmooth, sym=one )
!
!       SMOOTH POWER SPECTRUM OF THE SECOND SERIES.
!
        call dan_filter_rv( psvec2(:nf), nsmooth, sym=one )
!
!       SMOOTH CO-SPECTRUM.
!
        call dan_filter_rv( rcospect(:nf), nsmooth, sym=one )
!
!       SMOOTH QUADRATURE-SPECTRUM.
!
        call dan_filter_rv( icospect(:nf), nsmooth, sym=-one )
!
    end if
!
    if ( out_co_spect ) then
!
!       OUTPUT THE CO-SPECTRUM IF REQUIRED.
!
        co_spect(:nf) = rcospect(:nf)
!
    end if
!
    if ( out_quad_spect ) then
!
!       OUTPUT THE QUADRATURE-SPECTRUM IF REQUIRED.
!
        quad_spect(:nf) = -icospect(:nf)
!
    end if
!
!   COMPUTE PHASE OF THE CROSS-SPECTRUM.
!
    where ( rcospect(:nf)==zero .and. icospect(:nf)==zero )
        phase(:nf) = zero
    elsewhere
        wk(:nf)   = atan2( icospect(:nf), rcospect(:nf) )
        phase(:nf) = (one/twopi)*wk(:nf) + merge( zero, one, wk(:nf)>=zero ) 
    end where
!
!   COMPUTE MAGNITUDE OF THE CROSS-SPECTRUM.
!
    wk(:nf) = rcospect(:nf)**2 + icospect(:nf)**2
!
!   COMPUTE THE SQUARED COHERENCIES FOR ALL FREQUENCIES.
!
    where ( psvec(:nf)/=zero .and. psvec2(:nf)/=zero )
        rcospect(:nf)  = psvec(:nf)*psvec2(:nf)
        coher(:nf)     = min( wk(:nf)/rcospect(:nf), one )
    elsewhere
        coher(:nf) = zero
    end where 
!
    if ( out_ampli ) then
!
!       COMPUTE AND OUTPUT AMPLITUDE OF THE CROSS-SPECTRUM IF REQUIRED.
!
        ampli(:nf) = sqrt( wk(:nf) )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( n, stnd )
        freq(:nf) = arth(zero, c1, nf)
    end if
!
!   OUTPUT DEGREES OF FREEDOM, BANDWIDTH, CONFIDENCE LIMIT FACTORS 
!   OR PROBABILITIES FOR SQUARED COHERENCIES IF REQUIRED.
!
    if ( out_dof ) then
!
        if ( present(edof) ) then
            edof(:nf) = edof2(:nf)
        end if
!
        if ( present(bandwidth) ) then
            bandwidth(:nf) = edof2(:nf)/( two*real( n, stnd ) )
        end if
!
        if ( present(conlwr) .or. present(conupr) .or. present(testcoher) ) then
            call comp_conflim_rv( edof2(:nf), probtest=probtest, conlwr=conlwr,       &
                                  conupr=conupr, testcoher=testcoher )
        end if
!
!       COMPUTE AND OUTPUT SIGNIFICANCE PROBABILITY FOR SQUARED COHERENCIES.
!
        if ( out_prob_coher ) then
!
            sumw = nan()
!
            do i = 1_i4b, nf
!            
                df = edof2(i)  
!            
                if ( df>two ) then
!            
                    c1 = coher(i)
!            
                    if ( c1/=one ) then
!            
                        c2 = ( (df/two - one)*c1 )/( one - c1 )
                        df2 = df - two
!            
                        prob_coher(i) = probf2( c2, two, df2, true )
!            
                    else
!
                        prob_coher(i) = zero
!
                    end if 
!
                else
!
                    prob_coher(i)    = sumw
!
                end if
!
            end do
!
        end if
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( edof2 )
!
    end if
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE cross_spctrm_rv
! _________________________________
!
    end subroutine cross_spctrm_rv
!
! =========================================================================================
!
    subroutine cross_spctrm_rm( vec, mat, psvec, psmat, phase, coher, freq, edof, bandwidth,  &
                                conlwr, conupr, testcoher, ampli, co_spect, quad_spect,       &
                                prob_coher, initfft, normpsd, nsmooth, trend, win, taperp,    &
                                probtest )
!
! Purpose
! _______
!
!   Subroutine CROSS_SPCTRM computes Fast Fourier Transform (FFT) estimates
!   of the power and cross spectra of the real time series, VEC, and the multi-channel
!   real time series MAT.
!
!   The Power Spectral Density (PSD) and Cross Spectral Density (CSD) estimates are returned
!   in units which are the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the real time series for which the power and cross spectra
!                 must be estimated.
!                 If WIN/=2 or TREND=1, 2 or 3,  VEC is used as workspace and is transformed.
!
!                 Size(VEC) must be an even (positive) integer greater or equal to 4.
!
!   MAT           (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                 On entry, the multi-channel real time series for which the power 
!                 and cross spectra must be estimated. Each row of MAT is a real time series.
!                 If WIN/=2 or TREND=1, 2 or 3,  MAR is used as workspace and is transformed.
!
!                 The shape of MAT must verify:  size(MAT,2) = size(VEC).
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = size(VEC)/2 + 1 .
!
!   PSMAT         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the Power Spectral Density (PSD) estimates of each row of MAT.
!
!                 The shape of PSMAT must verify:
!
!                 - size(PSMAT,1) =  size(MAT,1) ;
!                 - size(PSMAT,2) =  size(VEC)/2 + 1 .
!
!   PHASE         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the phase of the cross spectrum, given in fractions
!                 of a circle (e.g. on the closed interval (0,1) ).
!
!                 The shape of PHASE must verify:
!
!                 - size(PHASE,1) =  size(MAT,1) ;
!                 - size(PHASE,2) =  size(VEC)/2 + 1 .
!
!   COHER         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the squared coherency  estimates for all frequencies.
!
!                 The shape of COHER must verify:
!
!                 - size(COHER,1) =  size(MAT,1) ;
!                 - size(COHER,2) =  size(VEC)/2 + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/size(VEC)
!                 for i=1,2, ... , (size(VEC)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = size(VEC)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!                 EDOF must verify:  size(EDOF) = size(VEC)/2 + 1 .
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the bandwidth of the power spectrum estimates.
!
!                 BANDWIDTH must verify:  size(BANDWIDTH) = size(VEC)/2 + 1 .
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the PSVEC(:) and
!                 PSMAT(:,:) arguments) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!                 CONLWR must verify:  size(CONLWR) = size(VEC)/2 + 1 .
!
!                 CONUPR must verify:  size(CONUPR) = size(VEC)/2 + 1 .
!   
!   TESTCOHER     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, this argument specifies the critical value for testing the null
!                 hypothesis that the squared coherency is zero at the PROBTEST * 100% significance
!                 level (e.g. elements of COHER(:,:) less than TESTCOHER(:) should be regarded as not
!                 significantly different from zero at the PROBTEST * 100% significance level).
!
!                 TESTCOHER must verify:  size(TESTCOHER) = size(VEC)/2 + 1 .
!
!   AMPLI         (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the cross-amplitude spectra.
!
!                 The shape of AMPLI must verify:
!
!                 - size(AMPLI,1) =  size(MAT,1) ;
!                 - size(AMPLI,2) =  (size(VEC)/2) + 1 .
!
!   CO_SPECT      (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the co-spectra (e.g. the real part of cross-spectra).  
!
!                 The shape of CO_SPECT must verify:
!
!                 - size(CO_SPECT,1) =  size(MAT,1) ;
!                 - size(CO_SPECT,2) =  (size(VEC)/2) + 1 .
!
!   QUAD_SPECT    (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the quadrature spectrum (e.g. the imaginary part of
!                 cross-spectrum with a minus sign).  
!
!                 The shape of QUAD_SPECT must verify:
!
!                 - size(QUAD_SPECT,1) =  size(MAT,1) ;
!                 - size(QUAD_SPECT,2) =  (size(VEC)/2) + 1 .
!
!   PROB_COHER    (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the probabilities that the computed sample squared coherencies
!                 came from an ergodic stationary bivariate process with (corresponding)
!                 squared coherencies equal to zero.
!
!                 The shape of PROB_COHER must verify:
!
!                 - size(PROB_COHER,1) =  size(MAT,1) ;
!                 - size(PROB_COHER,2) =  (size(VEC)/2) + 1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine CROSS_SPCTRM in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine CROSS_SPCTRM. This call to INITFFT must have the following form: 
!
!                      call init_fft( (/ size(MAT,1), size(MAT,2)/2 /), dim=2_i4b )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   CROSS_SPCTRM and a call to END_FFT is also done before leaving
!                   subroutine CROSS_SPCTRM.
!
!                 The default is INITFFT=true .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the power and cross spectra estimates are
!                   normalized in such a way that the total area under the power spectra is equal to
!                   the variance of the time series contained in VEC and in each row of MAT.
!                 - NORMPSD = false, the sum of the PSD estimates
!                   (e.g. sum(PSVEC(2:)) and sum(PSMAT(:,2:),dim=2) ) is equal to the variance of
!                   the corresponding time series.
!
!                 The default is NORMPSD=true .
!
!   NSMOOTH       (INPUT, OPTIONAL) integer(i4b)
!                 If NSMOOTH is used, the PSD and CSD estimates are computed by smoothing 
!                 the periodogram with Daniell weights (e.g. a simple moving average).
!
!                 On entry, NSMOOTH gives the length of the Daniell filter to be applied.
!
!                 Setting NSMOOTH=0 on entry is equivalent to omit the optional argument
!                 NSMOOTH. Otherwise, NSMOOTH must be odd, greater than 2 and less or
!                 equal to size(VEC)/2+1 .
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The means of the time series are removed before computing the power and cross spectra
!                 - TREND=+2 The drifts from time series are removed before computing the power and cross spectra
!                 - TREND=+3 The least-squares lines from time series are removed before
!                   computing the power and cross spectra.
!
!                 For other values of TREND nothing is done before estimating the power and cross spectra.
!
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power and cross spectra. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ) and the critical value for
!                 testing the null hypothesis that the squared coherency is zero (e.g.
!                 the TESTCOHER optional argument).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the selected
!   data window (e.g. WIN=1,2,3,4,5,6) is applied to the time series and the PSD and CSD estimates are
!   computed by the FFT of these transformed time series. Optionally, theses PSD and CSD estimates
!   may then be smoothed in the frequency domain by modified Daniell filters (e.g. if argument
!   NSMOOTH is used).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, assert_eq, arth
    use Reals_Constants,   only : zero, one, two, pi, twopi
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error14, tseries_error15,    &
                                  tseries_error21, tseries_error59, tseries_error61, tseries_error65
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : probf2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:)   :: vec
    real(stnd), intent(inout),  dimension(:,:) :: mat
!
    real(stnd), intent(out), dimension(:)             :: psvec
    real(stnd), intent(out), dimension(:,:)           :: psmat, phase, coher
    real(stnd), intent(in),                  optional :: taperp, probtest
    real(stnd), intent(out), dimension(:),   optional :: freq, edof, bandwidth, conlwr, conupr, testcoher
    real(stnd), intent(out), dimension(:,:), optional :: ampli, co_spect, quad_spect, prob_coher
!
    integer(i4b),  intent(in),               optional :: trend, win, nsmooth
!
    logical(lgl),  intent(in), optional :: normpsd, initfft
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, n, nd2, nf, k, k1, k2, win2, trend2, i
    integer      :: iok
!
    real(stnd)                                :: c1, c2, sumw, df, df2
    real(stnd), dimension(size(vec))          :: wk
    real(stnd), dimension(:),   allocatable   :: temp, edof2
    real(stnd), dimension(:,:), allocatable   :: magni
!
    complex(stnd), dimension(:),   allocatable :: cwk1
    complex(stnd), dimension(:,:), allocatable :: cwk2
!
    logical(lgl)  :: normpsd2, initfft2, smooth, out_freq, out_ampli, out_co_spect,    &
                     out_quad_spect, out_prob_coher, out_dof
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='cross_spctrm'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m =  assert_eq( int(size(mat,1),i4b) ,          &
                    int(size(psmat,1),i4b) ,        &
                    int(size(phase,1),i4b) ,        &
                    int(size(coher,1),i4b) ,        &
                    name_proc )
!
    if ( m<=0_i4b )     &
    call merror( name_proc//tseries_error21 )
!
    n =  assert_eq( int(size(vec),i4b) ,          &
                    int(size(mat,2),i4b) ,        &
                    name_proc )
!
    if ( n<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    nd2 = n/2_i4b
!
    if ( n/=2*nd2 )   &
    call merror( name_proc//tseries_error14 )
!
    nf  = nd2 + 1_i4b
!
    call assert( logical(int(size(psvec),i4b)==nf,lgl),     &
                 logical(int(size(psmat,2),i4b)==nf,lgl),   &
                 logical(int(size(phase,2),i4b)==nf,lgl),   &
                 logical(int(size(coher,2),i4b)==nf,lgl),   &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_ampli = present( ampli )
!
    if ( out_ampli ) then
        call assert( logical(int(size(ampli,1),i4b)==m,lgl),   &
                     logical(int(size(ampli,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_co_spect = present( co_spect )
!
    if ( out_co_spect ) then
        call assert( logical(int(size(co_spect,1),i4b)==m,lgl),   &
                     logical(int(size(co_spect,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_quad_spect = present( quad_spect )
!
    if ( out_quad_spect ) then
        call assert( logical(int(size(quad_spect,1),i4b)==m,lgl),   &
                     logical(int(size(quad_spect,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_prob_coher = present( prob_coher )
    if ( out_prob_coher ) then
        call assert( logical(int(size(prob_coher,1),i4b)==m,lgl),   &
                     logical(int(size(prob_coher,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth ) .or. out_prob_coher   &
              .or. present( conlwr ) .or. present( conupr ) .or. present( testcoher ) 
!
    if ( out_dof ) then
!
        if ( present(probtest) ) then
            if ( zero>=probtest .or. probtest>=one ) then
                call merror( name_proc//tseries_error59 )
            end if
        end if
!
        if ( present( edof ) ) then
            call assert( logical(int(size(edof),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( bandwidth ) ) then
            call assert( logical(int(size(bandwidth),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conlwr ) ) then
            call assert( logical(int(size(conlwr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conupr ) ) then
            call assert( logical(int(size(conupr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( testcoher ) ) then
            call assert( logical(int(size(testcoher),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
    end if
!
    if ( present(trend) ) then
        trend2 = trend
    else
        trend2 = 1_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
!
    else
!
        win2 = 3_i4b
!
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( (/ m, nd2/), dim=2_i4b )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    smooth = false
!
    if ( present(nsmooth) ) then
!
        if ( nsmooth/=0_i4b ) then
!
!           CHECK THE INPUT VALUE FOR THE LENGTH OF THE DANIELL FILTER.
!
            if ( nsmooth<3_i4b .or. nsmooth>nf )      &
            call merror( name_proc//tseries_error65 )
!
            if ( (nsmooth/2_i4b)*2_i4b==nsmooth )   &
            call merror( name_proc//tseries_error61 )
!
            smooth = true
!
        end if
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trend2>=1_i4b .and. trend2<=3_i4b ) then
        call detrend_rv( vec(:n),    trend2 )
        call detrend_rm( mat(:m,:n), trend2 )
    end if
!
!   TAPER THE TIME SERIES IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:n)  = data_window( n, win2, taperp=taperp )
!
!       APPLY DATA WINDOW.
!
        vec(:n) = vec(:n)*wk(:n)
!
        mat(:m,:n) = mat(:m,:n)*spread( wk(:n), dim=1, ncopies=m )
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:n), wk(:n) )
!
    else    
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:n)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( n, stnd )
!
    end if
!
!   COMPUTE PSD AND CSD ESTIMATES.
!
!   ALLOCATE WORK MATRIX.
!
    allocate( cwk1(nf), cwk2(m,nf), stat = iok )
!
    if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!   COMPUTE FFT FOR THE SERIES.
!
    call real_fft( vec(:n), cwk1(:nf), true)
!
    call real_fft( mat(:m,:n), cwk2(:m,:nf), true )
!
!   COMPUTE POWER SPECTRUM ESIMATES.
!
    psvec(:nd2) = real(cwk1(:nd2),stnd)**2 + aimag(cwk1(:nd2))**2
    psvec(nf)   = real(cwk1(nf),stnd)**2
!
    psmat(:m,:nd2) = real(cwk2(:m,:nd2),stnd)**2  + aimag(cwk2(:m,:nd2))**2
    psmat(:m,nf)   = real(cwk2(:m,nf),stnd)**2
!
!   COMPUTE CO-SPECTRUM ESIMATES.
!
    coher(:m,:nd2)  = spread(real(cwk1(:nd2),stnd),dim=1,ncopies=m)*real(cwk2(:m,:nd2),stnd) +   &
                      spread(aimag(cwk1(:nd2)),    dim=1,ncopies=m)*aimag(cwk2(:m,:nd2))
    coher(:m,nf)    = real(cwk1(nf),stnd)*real(cwk2(:m,nf),stnd)
!
!   COMPUTE QUADRATURE-SPECTRUM ESIMATES.
!
    phase(:m,:nd2)  = spread(aimag(cwk1(:nd2)),    dim=1,ncopies=m)*real(cwk2(:m,:nd2),stnd) -   &
                      spread(real(cwk1(:nd2),stnd),dim=1,ncopies=m)*aimag(cwk2(:m,:nd2))
    phase(:m,nf)    = zero
!
!   DEALLOCATE WORK ARRAYS.
!
    deallocate( cwk1, cwk2 )
!
!   NORMALIZE THE POWER AND CROSS SPECTRA ESTIMATES.
!
    c1 = one/(sumw*real( n, stnd))
!
    if ( normpsd2 ) then
!
!       CHATFIELD DEFINITION OF THE SPECTRAL DENSITY FUNCTION.
!
        c2 = (real( n, stnd)/pi)*c1
        k1 = 1_i4b
        k2 = nf
!
    else
!
!       POLLOCK DEFINITION OF THE PERIODOGRAM. IF THE ESTIMATES
!       ARE DIVIDED BY 4*PI, THIS IS THE POLLOCK DEFINITION OF
!       THE SPECTRAL DENSITY FUNCTION.
!
        c2 = two*c1
        k1 = 2_i4b
        k2 = nd2
!
        psvec(1_i4b)    = c1*psvec(1_i4b)
        psmat(:m,1_i4b) = c1*psmat(:m,1_i4b)
!
        coher(:m,1_i4b)  = c1*coher(:m,1_i4b)
        phase(:m,1_i4b)  = c1*phase(:m,1_i4b)
!
        psvec(nf)    = c1*psvec(nf)
        psmat(:m,nf) = c1*psmat(:m,nf)
!
        coher(:m,nf)  = c1*coher(:m,nf)
        phase(:m,nf)  = c1*phase(:m,nf)
!
    end if
!
    psvec(k1:k2)    = c2*psvec(k1:k2)
    psmat(:m,k1:k2) = c2*psmat(:m,k1:k2)
!
    coher(:m,k1:k2)   = c2*coher(:m,k1:k2)
    phase(:m,k1:k2)   = c2*phase(:m,k1:k2)
!
!   SMOOTH THE POWER AND CROSS SPECTRA ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
!       SMOOTH POWER SPECTRUM OF THE FIRST SERIES.
!
        call dan_filter_rv( psvec(:nf), nsmooth, sym=one )
!
!       SMOOTH POWER SPECTRUM OF THE OTHER SERIES.
!
        call dan_filter_rm( psmat(:m,:nf), nsmooth, sym=one )
!
!       SMOOTH CO-SPECTRUM.
!
        call dan_filter_rm( coher(:m,:nf), nsmooth, sym=one )
!
!       SMOOTH QUADRATURE-SPECTRUM.
!
        call dan_filter_rm( phase(:m,:nf), nsmooth, sym=-one )
!
    end if
!
    if ( out_co_spect ) then
!
!       OUTPUT THE CO-SPECTRUM IF REQUIRED.
!
        co_spect(:m,:nf) = coher(:m,:nf)
!
    end if
!
    if ( out_quad_spect ) then
!
!       OUTPUT THE QUADRATURE-SPECTRUM IF REQUIRED.
!
        quad_spect(:m,:nf) = -phase(:m,:nf)
!
    end if
!
!   ALLOCATE WORK MATRICES.
!
    allocate( temp(m), magni(m,nf), stat = iok )
!
    if ( iok/=0 ) call merror( name_proc//allocate_error )
!
    do k = 1_i4b, nf
!
!       COMPUTE MAGNITUDE OF THE CROSS-SPECTRUM.
!
        magni(:m,k) = coher(:m,k)**2 + phase(:m,k)**2
!
!       COMPUTE PHASE OF THE CROSS-SPECTRUM.
!
        where ( coher(:m,k)/=zero .or. phase(:m,k)/=zero )
            temp(:m)   = atan2( phase(:m,k), coher(:m,k) )
            phase(:m,k) = (one/twopi)*temp(:m) + merge( zero, one, temp(:m)>=zero ) 
        end where
!
!       COMPUTE THE SQUARED COHERENCY.
!
        if ( psvec(k)/=zero ) then
!
            where( psmat(:m,k)/=zero )
                coher(:m,k) = min( magni(:m,k)/(psvec(k)*psmat(:m,k)), one )
            elsewhere
                coher(:m,k) = zero
            end where
!
        else
!
            coher(:m,k) = zero
!
        end if
!
    end do   
!
!   OUTPUT AMPLITUDE OF THE CROSS-SPECTRUM IF REQUIRED.
!
    if ( out_ampli ) then
!
        ampli(:m,:nf) = sqrt( magni(:m,:nf) )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
!
        c1 = one/real( n, stnd )
        freq(:nf) = arth(zero, c1, nf)
!
    end if
!
!   OUTPUT DEGREES OF FREEDOM, BANDWIDTH, CONFIDENCE LIMIT FACTORS 
!   OR PROBABILITIES FOR SQUARED COHERENCIES IF REQUIRED.
!
    if ( out_dof ) then
!
        allocate( edof2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
        edof2(:nf) = estim_dof2( wk(:n), 0_i4b, nsmooth=nsmooth )
!
        if ( present(edof) ) then
            edof(:nf) = edof2(:nf)
        end if
!
        if ( present(bandwidth) ) then
            bandwidth(:nf) = edof2(:nf)/( two*real( n, stnd ) )
        end if
!
        if ( present(conlwr) .or. present(conupr) .or. present(testcoher) ) then
            call comp_conflim_rv( edof2(:nf), probtest=probtest, conlwr=conlwr,       &
                                  conupr=conupr, testcoher=testcoher )
        end if
!
!       COMPUTE AND OUTPUT SIGNIFICANCE PROBABILITY FOR SQUARED COHERENCIES.
!
        if ( out_prob_coher ) then
!
            sumw = nan()
!
            do i = 1_i4b, nf
!            
                df = edof2(i)   
!            
                if ( df>two ) then
!            
                    where ( coher(:m,i)/=one )
                        temp(:m) = ( (df/two - one)*coher(:m,i) )/( one - coher(:m,i) )
                    elsewhere
                        temp(:m) = zero
                    end where
!            
                    df2 = df - two
!            
                    prob_coher(:m,i) = probf2( temp(:m), two, df2, true )
!            
                    where( coher(:m,i)==one ) prob_coher(:m,i) = zero
!
                else
!
                    prob_coher(:m,i) = sumw
!
                end if
!
            end do
!
        end if
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( edof2 )
!
    end if
!
!   DEALLOCATE WORK ARRAYS.
!
    deallocate( temp, magni )
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE cross_spctrm_rm
! _________________________________
!
    end subroutine cross_spctrm_rm
!
! =========================================================================================
!
    subroutine power_spctrm2_rv( vec, l, psvec, freq, edof, bandwidth, conlwr, conupr, initfft,      &
                                 overlap, normpsd, nsmooth, trend, trend2, win, taperp, l0, probtest )
!
! Purpose
! _______
!
!   Subroutine POWER_SPCTRM2 computes a Fast Fourier Transform (FFT) estimate
!   of the power spectrum of a real time series.
!
!   The Power Spectral Density (PSD) estimates are returned in units which are
!   the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the real time series for which the power spectrum
!                 must be estimated.
!                 If TREND=1, 2 or 3,  VEC is used as workspace and is transformed.
!
!                 Size(VEC) must be greater or equal to 4.
!
!   L             (INPUT) integer(i4b)
!                 On entry, an integer used to segment the time series. L is the
!                 length of the segments. L must be a positive even integer, less
!                 or equal to size(VEC), but greater or equal to 4.
!
!                 Spectral computations are at (L/2)+1 frequencies if the optional
!                 argument L0 is absent and are at ((L+L0)/2)+1 frequencies if L0 is
!                 present (L0 is the number of zeros added to each segment).
!
!                 Suggested values for L+L0 are 16, 32, 64 or 128 (e.g. an integer power of two,
!                 in order to speed the computations).
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = ((L+L0)/2) + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/(L+L0)
!                 for i=1,2, ... , ((L+L0)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = (L+L0)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!                 EDOF must verify:  size(EDOF) = ((L+L0)/2) + 1 .
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the bandwidth of the power spectrum estimates.
!
!                 BANDWIDTH must verify:  size(BANDWIDTH) = ((L+L0)/2) + 1 .
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the
!                 PSVEC(:) argument) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!                 CONLWR must verify:  size(CONLWR) = ((L+L0)/2) + 1 .
!
!                 CONUPR must verify:  size(CONUPR) = ((L+L0)/2) + 1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine POWER_SPCTRM2 in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine POWER_SPCTRM2. This call to INITFFT must have the following form: 
!
!                      call init_fft( (L+L0)/2 )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   POWER_SPCTRM2 and a call to END_FFT is also done before leaving
!                   subroutine POWER_SPCTRM2.
!
!                 The default is INITFFT=true .
!
!   OVERLAP       (INPUT, OPTIONAL) logical(lgl)
!                 If:
!
!                 - OVERLAP = false, the subroutine segments the data 
!                   without any overlapping.
!                 - OVERLAP=true, the subroutine overlaps the segments
!                   by one half of their length (which is equal to L).
!
!                 In both cases, zeros are eventually added to each segment (if argument L0 is present)
!                 and each segment will be FFT'd, and the resulting periodograms
!                 will averaged together to obtain a Power Spectrum Density estimate at the
!                 ((L+L0)/2)+1 frequencies.
!
!                 The default is OVERLAP=false .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the PSD estimates are normalized in such
!                   a way that the total area under the power spectrum is equal to the variance of
!                   the time series VEC.
!                 - NORMPSD is set to false, the sum of the PSD estimates
!                   (e.g. sum( PSVEC(2:) ) is equal to the variance of the time series.
!
!                 The default is NORMPSD=true .
!
!   NSMOOTH       (INPUT, OPTIONAL) integer(i4b)
!                 If NSMOOTH is used, the PSD estimates are computed by smoothing 
!                 the periodogram with Daniell weights (e.g. a simple moving average).
!
!                 On entry, NSMOOTH gives the length of the Daniell filter to be applied.
!
!                 Setting NSMOOTH=0 on entry is equivalent to omit the optional argument
!                 NSMOOTH. Otherwise, NSMOOTH must be odd, greater than 2 and less or
!                 equal to (L+L0)/2+1 .
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 -  TREND=+1 The mean of the time series is removed before computing the spectrum
!                 -  TREND=+2 The drift from the time series is removed before computing the spectrum
!                    by using the formula: drift = (VEC(size(VEC)) - VEC(1))/(size(VEC) - 1)
!                 -  TREND=+3 The least-squares line from the time series is removed before
!                    computing the spectrum.
!
!                 For other values of TREND nothing is done before estimating the power spectrum.
!
!                 The default is TREND=1, e.g. the mean of the time series is removed before the
!                 computations.
!
!   TREND2        (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 -  TREND2=+1 The mean of the time segment is removed before computing the spectrum
!                    on this segment.
!                 -  TREND2=+2 The drift from the time segment is removed before computing the spectrum
!                    on this segment.
!                 -  TREND2=+3 The least-squares line from the time segment is removed before
!                    computing the spectrum on this segment.
!
!                 For other values of TREND2 nothing is done before estimating the power spectrum
!                 on each segment.
!
!                 The default is TREND2=0, e.g. nothing is done before estimating the power spectrum
!                 on each segment.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power spectrum. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!
!   L0            (INPUT, OPTIONAL) integer(i4b)
!                 The number of zeros added to each time segment in order to obtain more finely
!                 spaced spectral estimates. L+L0 must be a positive even integer.
!
!                 The default is L0=0, e.g. no zeros are added to each time segment.
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the series
!   is padded with zero on the right such that the length of the resulting time series is 
!   evenly divisible by L (a positive even integer). The length, N, of this resulting time
!   series is the first integer greater than or equal to size(VEC) which is evenly divisible
!   by L. If size(VEC) is not evenly divisible by L, N is equal to size(VEC)+L-mod(size(VEC),L).
!   
!   Optionally, the mean or the trend may also be removed from each time segment (e.g. TREND2=1,2,3).
!   Optionally, zeros may be added to each time segment (e.g. the optional arguemnt L0) if more
!   finely spaced spectral esimates are desired.
!   
!   The stability of the PSD estimates depends on the averaging process. That is, the greater
!   the number of segments ( N/L if OVERLAP=false and (2N/L)-1 if OVERLAP=true), the more
!   stable the resulting PSD estimates.
!   
!   Optionally, theses PSD estimates may then be smoothed again in the frequency
!   domain by a Daniell filter (e.g. if argument NSMOOTH is used).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, arth
    use Reals_Constants,   only : zero, one, two, pi
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error15, tseries_error17,    &
                                  tseries_error18, tseries_error24, tseries_error25, tseries_error58,   &
                                  tseries_error59, tseries_error61, tseries_error67
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout), dimension(:) :: vec
    real(stnd), intent(out),   dimension(:) :: psvec
!
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out), dimension(:), optional :: freq, edof, bandwidth, conlwr, conupr
!
    integer(i4b),  intent(in)                         :: l
    integer(i4b),  intent(in),               optional :: trend, trend2, win, l0, nsmooth
!
    logical(lgl),  intent(in), optional :: initfft, overlap, normpsd
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, n2, nf, k, l2, ld2, m, win2, trendb, trend2b,    &
                    i, i1, i2, step, sl, sld2
    integer      :: iok
!
    real(stnd)                            :: c1, c2, sumw
    real(stnd), dimension(l)              :: wk
    real(stnd), dimension(:), allocatable :: seg, edof2
!
    complex(stnd), dimension(:), allocatable  :: cwk
!
    logical(lgl)  :: normpsd2, overlap2, initfft2, smooth, out_freq, out_dof
!
#ifdef _OPENMP
!
    real(stnd), dimension(:), allocatable :: psvecb
!
    logical  :: test_par
!
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='power_spctrm2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n   = size( vec )
!
    if ( l<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    if ( l>n )     &
    call merror( name_proc//tseries_error17 )
!
!   CHECK IF l IS A POSITIVE EVEN INTEGER.
!
    ld2 = l/2_i4b
!
    if ( l/=2_i4b*ld2 )     &
    call merror( name_proc//tseries_error18 )
!
!   DETERMINE THE LENGTH OF THE SEGMENTS.
!
    if ( present(l0) ) then
!
        if ( l0<0_i4b )     &
        call merror( name_proc//tseries_error24 )
!
        sl = l + l0
!
!       CHECK IF sl IS A POSITIVE EVEN INTEGER.
!
        sld2 = sl/2_i4b
!
        if ( sl/=2_i4b*sld2 )     &
        call merror( name_proc//tseries_error25  )
!
    else
!
        sl   = l
        sld2 = ld2
!
    end if
!
    nf  = sld2 + 1_i4b
!
    call assert( logical(int(size(psvec),i4b)==nf,lgl),    &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth )        &
              .or. present( conlwr ) .or. present( conupr ) 
!
    if ( out_dof ) then
!
        if ( present(probtest) ) then
            if ( zero>=probtest .or. probtest>=one ) then
                call merror( name_proc//tseries_error59 )
            end if
        end if
!
        if ( present( edof ) ) then
            call assert( logical(int(size(edof),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( bandwidth ) ) then
            call assert( logical(int(size(bandwidth),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conlwr ) ) then
            call assert( logical(int(size(conlwr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conupr ) ) then
            call assert( logical(int(size(conupr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
    end if
!
    if ( present(trend) ) then
        trendb = trend
    else
        trendb = 1_i4b
    end if
!
    if ( present(trend2) ) then
        trend2b = trend2
    else
        trend2b = 0_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
!
    else
!
        win2 = 3_i4b
!
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( sld2 )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    overlap2 = false
!
    if ( present(overlap) ) then
        overlap2 = overlap
    end if
!
    smooth = false
!
    if ( present(nsmooth) ) then
!
        if ( nsmooth/=0_i4b ) then
!
!           CHECK THE INPUT VALUE FOR THE LENGTH OF THE DANIELL FILTER.
!
            if ( nsmooth<3_i4b .or. nsmooth>nf )      &
            call merror( name_proc//tseries_error67 )
!
            if ( (nsmooth/2_i4b)*2_i4b==nsmooth )   &
            call merror( name_proc//tseries_error61 )
!
            smooth = true
!
        end if
!
    end if
!
!   FIND THE FIRST INTEGER GREATER THAN OR EQUAL TO n THAT IS EVENLY
!   DIVISIBLE BY l.
!
    i = mod( n, l )
!
    if ( i==0_i4b ) then
        n2 = n
    else
        n2 = n + l - i
    end if
!
!   DETERMINE THE NUMBER OF SEGMENTS.
!
    if ( overlap2 ) then
        m    = (2*n2/l) - 1_i4b
        step = ld2
    else
        m    = n2/l
        step = l
    end if
!
    if ( out_dof ) then
!
        if ( m/=1_i4b .and. overlap2 .and. win2/=1_i4b .and. win2/=3_i4b  )    &
        call merror( name_proc//tseries_error58 )
!
    end if
!
!   CALCULATE DATA WINDOW IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:l)  = data_window( l, win2, taperp=taperp )
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:l), wk(:l) )
!
    else
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:l)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( l, stnd )
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trendb>=1_i4b .and. trendb<=3_i4b ) then
        call detrend_rv( vec(:n), trendb )
    end if
!
!   ZERO OUTPUT PSD ESTIMATES FOR SUMMING OVER THE m SEGMENTS.
!
    psvec(:nf) = zero
!
!   COMPUTE PSD ESTIMATES BY SUMMING OVER THE m SEGMENTS.
!
#ifdef _OPENMP
    i1 = omp_get_num_procs()
    i2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               i1>1_i4b                     .and.      &
               i2>1_i4b                     .and.      &
               m>=i2
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(k,i,i1,i2,l2,iok,psvecb,seg,cwk)
!
!       ALLOCATE WORK VECTORS.
!
        allocate( psvecb(nf), seg(sl), cwk(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ZERO PARTIAL OUTPUT PSD ESTIMATES FOR SUMMING OVER EACH THREAD.
!
        psvecb(:nf) = zero
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(l+1_i4b:sl) = zero
        end if    
!
!$OMP DO SCHEDULE(STATIC)
!
        do k=1_i4b, m
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg(:l2) = vec(i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENT IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rv( seg(:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:l2) = seg(:l2)*wk(:l2)
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT OF THE TIME SEGMENT.
!
            call real_fft( seg(:sl), cwk(:nf), true)
!
!           UPDATE SPECTRAL ESTIMATES.
!
            psvecb(:sld2) = psvecb(:sld2) + real(cwk(:sld2),stnd)**2 + aimag(cwk(:sld2))**2
            psvecb(nf)    = psvecb(nf)    + real(cwk(nf),stnd)**2
!        
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (updatepsd)
        psvec(:nf) = psvec(:nf) + psvecb(:nf)
!$OMP END CRITICAL (updatepsd)
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( psvecb, seg, cwk )
!
!$OMP END PARALLEL
!
    else
!
#endif
!
!       ALLOCATE WORK VECTORS.
!
        allocate( seg(sl), cwk(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(l+1_i4b:sl) = zero
        end if    
!
        do k=1_i4b, m
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
            l2 = i2 - i1 + 1_i4b
!
            seg(:l2) = vec(i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENT IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rv( seg(:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:l2) = seg(:l2)*wk(:l2)
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT OF THE TIME SEGMENT.
!
            call real_fft( seg(:sl), cwk(:nf), true)
!
            psvec(:sld2) = psvec(:sld2) + real(cwk(:sld2),stnd)**2 + aimag(cwk(:sld2))**2
            psvec(nf)    = psvec(nf)    + real(cwk(nf),stnd)**2
!        
        end do
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( seg, cwk )
!        
#ifdef _OPENMP
!
    end if
!
#endif
!
!   NORMALIZE THE PSD ESTIMATES.
!
    c1 = one/( sumw*real(l*m, stnd) )
!
    if ( normpsd2 ) then
!
!       CHATFIELD DEFINITION OF THE SPECTRAL DENSITY FUNCTION.
!
        c2 = real( l, stnd)/pi
        psvec(1_i4b:nf) = (c1*c2)*psvec(1_i4b:nf)
!
    else
!
!       POLLOCK DEFINITION OF THE PERIODOGRAM. IF THE ESTIMATES
!       ARE DIVIDED BY 4*PI, THIS IS THE POLLOCK DEFINITION OF
!       THE SPECTRAL DENSITY FUNCTION.
!
        psvec(1_i4b)      = c1*psvec(1_i4b)
        psvec(2_i4b:sld2) = (two*c1)*psvec(2_i4b:sld2)
        psvec(nf)         = c1*psvec(nf)
!
    end if
!
!   SMOOTH THE PSD ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
        call dan_filter_rv( psvec(:nf), nsmooth, sym=one )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( sl, stnd )
        freq(:nf) = arth(zero, c1, nf)
    end if
!
!   OUTPUT DEGREES OF FREEDOM OR BANDWIDTH IF REQUIRED.
!
    if ( out_dof ) then
!
        allocate( edof2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
        i = sl - l
!
        edof2(:nf) = estim_dof2( wk(:l), i, win=win, nsmooth=nsmooth, nseg=m, overlap=overlap )
!
        if ( present(conlwr) .or. present(conupr) ) then
            call comp_conflim_rv( edof2(:nf), probtest=probtest, conlwr=conlwr, conupr=conupr )
        end if
!
        if ( present(edof) ) then
            edof(:nf) = edof2(:nf)
        end if
!
        if ( present(bandwidth) ) then
            bandwidth(:nf) = edof2(:nf)/( two*real( n, stnd ) )
        end if
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( edof2 )
!
    end if
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE power_spctrm2_rv
! __________________________________
!
    end subroutine power_spctrm2_rv
!
! =========================================================================================
!
    subroutine power_spctrm2_rm( mat, l, psmat, freq, edof, bandwidth, conlwr, conupr, initfft,      &
                                 overlap, normpsd, nsmooth, trend, trend2, win, taperp, l0, probtest )
!
! Purpose
! _______
!
!   Subroutine POWER_SPCTRM2 computes Fast Fourier Transform (FFT) estimates
!   of the power spectra of the multi-channel real time series MAT (e.g. each row
!   of MAT contains a time series).
!
!   The Power Spectral Density (PSD) estimates are returned in units which are
!   the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   MAT           (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                 On entry, the multi-channel real time series for which the power 
!                 spectra must be estimated. Each row of MAT is a real time series.
!                 If TREND=1, 2 or 3,  MAT is used as workspace and is transformed.
!
!                 Size(MAT,2) must be greater or equal to 4.
!
!   L             (INPUT) integer(i4b)
!                 On entry, an integer used to segment the time series. L is the
!                 length of the segments. L must be a positive even integer less
!                 or equal to size(MAT,2), but greater or equal to 4.
!
!                 Spectral computations are at (L/2)+1 frequencies if the optional
!                 argument L0 is absent and are at ((L+L0)/2)+1 frequencies if L0 is
!                 present (L0 is the number of zeros added to each segment).
!
!                 Suggested values for L+L0 are 16, 32, 64 or 128 (e.g. an integer power of two,
!                 in order to speed the computations).
!
!   PSMAT         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the Power Spectral Density (PSD) estimates of each row of MAT.
!
!                 The shape of PSMAT must verify:
!
!                 - size(PSMAT,1) =  size(MAT,1) ;
!                 - size(PSMAT,2) = ((L+L0)/2) + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/(L+L0)
!                 for i=1,2, ... , ((L+L0)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = (L+L0)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!                 EDOF must verify:  size(EDOF) = ((L+L0)/2) + 1 .
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the bandwidth of the power spectrum estimates.
!
!                 BANDWIDTH must verify:  size(BANDWIDTH) = ((L+L0)/2) + 1 .
!
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the
!                 PSMAT(:,:) argument) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!                 CONLWR must verify:  size(CONLWR) = ((L+L0)/2) + 1 .
!
!                 CONUPR must verify:  size(CONUPR) = ((L+L0)/2) + 1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine POWER_SPCTRM2 in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine POWER_SPCTRM2. This call to INITFFT must have the following form: 
!
!                      call init_fft( (/ size(MAT,1), (L+L0)/2 /), dim=2_i4b )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   POWER_SPCTRM2 and a call to END_FFT is also done before leaving
!                   subroutine POWER_SPCTRM2.
!
!                 The default is INITFFT=true .
!
!   OVERLAP       (INPUT, OPTIONAL) logical(lgl)
!                 If:
!
!                 - OVERLAP = false, the subroutine segments the data 
!                   without any overlapping.
!                 - OVERLAP=true, the subroutine overlaps the segments
!                   by one half of their length (which is equal to L).
!
!                 In both cases, zeros are eventually added to each segment (if argument L0 is present)
!                 and each segment will be FFT'd, and the resulting periodograms
!                 will averaged together to obtain a Power Spectrum Density estimate at the
!                 ((L+L0)/2)+1 frequencies.
!
!                 The default is OVERLAP=false .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the PSD estimates are normalized in such
!                   a way that the total area under the power spectrum is equal to the variance of
!                   the time series VEC.
!                 - NORMPSD is set to false, the sum of the PSD estimates
!                   (e.g. sum(PSMAT(:,2:),dim=2) ) is equal to the variance of the corresponding time series.
!
!                 The default is NORMPSD=true .
!
!   NSMOOTH       (INPUT, OPTIONAL) integer(i4b)
!                 If NSMOOTH is used, the PSD estimates are computed by smoothing 
!                 the periodogram with Daniell weights (e.g. a simple moving average).
!
!                 On entry, NSMOOTH gives the length of the Daniell filter to be applied.
!
!                 Setting NSMOOTH=0 on entry is equivalent to omit the optional argument
!                 NSMOOTH. Otherwise, NSMOOTH must be odd, greater than 2 and less or
!                 equal to (L+L0)/2+1 .
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The means of the time series are removed before computing the spectra
!                 - TREND=+2 The drifts from time series are removed before computing the spectra
!                 - TREND=+3 The least-squares lines from time series are removed before
!                   computing the spectra.
!
!                 For other values of TREND nothing is done before estimating the power and cross spectra.
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   TREND2        (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND2=+1 The mean of the time segment is removed before computing the spectrum
!                   on this segment.
!                 - TREND2=+2 The drift from the time segment is removed before computing the spectrum
!                   on this segment.
!                 - TREND2=+3 The least-squares line from the time segment is removed before
!                   computing the spectrum on this segment.
!
!                 For other values of TREND2 nothing is done before estimating the power spectrum
!                 on each segment.
!
!                 The default is TREND2=0, e.g. nothing is done before estimating the power spectrum
!                 on each segment.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power spectrum. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!
!   L0            (INPUT, OPTIONAL) integer(i4b)
!                 The number of zeros added to each time segment in order to obtain more finely
!                 spaced spectral estimates. L+L0 must be a positive even integer.
!
!                 The default is L0=0, e.g. no zeros are added to each time segment.
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the series
!   are padded with zero on the right such that the length of the resulting time series is 
!   evenly divisible by L (a positive even integer). The length, N, of this resulting time
!   series is the first integer greater than or equal to size(MAT,2) which is evenly divisible
!   by L. If size(MAT,2) is not evenly divisible by L, N is equal to size(MAT,2)+L-mod(size(MAT,2),L).
!   
!   Optionally, the mean or the trend may also be removed from each time segment (e.g. TREND2=1,2,3).
!   Optionally, zeros may be added to each time segment (e.g. the optional arguemnt L0) if more
!   finely spaced spectral esimates are desired.
!   
!   The stability of the PSD estimates depends on the averaging process. That is, the greater
!   the number of segments ( N/L if OVERLAP=false and (2N/L)-1 if OVERLAP=true), the more
!   stable the resulting PSD estimates.
!   
!   Optionally, theses PSD estimates may then be smoothed again in the frequency
!   domain by modified Daniell filters (e.g. if argument NSMOOTH is used).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, assert_eq, merror, arth
    use Reals_Constants,   only : zero, one, two, pi
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error15, tseries_error18,    &
                                  tseries_error19, tseries_error24, tseries_error25, tseries_error58,   &
                                  tseries_error59, tseries_error61, tseries_error67
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:,:) :: mat
    real(stnd), intent(out),    dimension(:,:) :: psmat
!
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out), dimension(:), optional :: freq, edof, bandwidth, conlwr, conupr
!
    integer(i4b),  intent(in)                         :: l
    integer(i4b),  intent(in),               optional :: trend, trend2, win, l0, nsmooth
!
    logical(lgl),  intent(in), optional :: initfft, overlap, normpsd
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, n, n2, nf, k, l2, ld2, nseg, win2, trendb, trend2b,   &
                    i, i1, i2, step, sl, sld2
    integer      :: iok
!
    real(stnd)                              :: c1, c2, sumw
    real(stnd), dimension(l)                :: wk
    real(stnd), dimension(:),   allocatable :: edof2
    real(stnd), dimension(:,:), allocatable :: seg
!
    complex(stnd), dimension(:,:), allocatable :: cwk
!
    logical(lgl)  :: normpsd2, initfft2, overlap2, smooth, out_freq, out_dof
!
#ifdef _OPENMP
!
    real(stnd), dimension(:,:), allocatable :: psmatb
!
    logical  :: test_par
!
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='power_spctrm2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m =  assert_eq( int(size(mat,1),i4b) ,          &
                    int(size(psmat,1),i4b) ,        &
                    name_proc )
    if ( m<=0_i4b ) return
!
    n =  size( mat, 2 )
!
    if ( l<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    if ( l>n )     &
    call merror( name_proc//tseries_error19 )
!
!   CHECK IF l IS A POSITIVE EVEN INTEGER.
!
    ld2 = l/2_i4b
!
    if ( l/=2_i4b*ld2 )     &
    call merror( name_proc//tseries_error18 )
!
!   DETERMINE THE LENGTH OF THE SEGMENTS.
!
    if ( present(l0) ) then
!
        if ( l0<0_i4b )     &
        call merror( name_proc//tseries_error24 )
!
        sl = l + l0
!
!       CHECK IF sl IS A POSITIVE EVEN INTEGER.
!
        sld2 = sl/2_i4b
!
        if ( sl/=2_i4b*sld2 )     &
        call merror( name_proc//tseries_error25 )
!
    else
!
        sl   = l
        sld2 = ld2
!
    end if
!
    nf  = sld2 + 1_i4b
    call assert( logical(int(size(psmat,2),i4b)==nf,lgl),    &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth )        &
              .or. present( conlwr ) .or. present( conupr ) 
!
    if ( out_dof ) then
!
        if ( present(probtest) ) then
            if ( zero>=probtest .or. probtest>=one ) then
                call merror( name_proc//tseries_error59 )
            end if
        end if
!
        if ( present( edof ) ) then
            call assert( logical(int(size(edof),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( bandwidth ) ) then
            call assert( logical(int(size(bandwidth),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conlwr ) ) then
            call assert( logical(int(size(conlwr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conupr ) ) then
            call assert( logical(int(size(conupr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
    end if
!
    if ( present(trend) ) then
        trendb = trend
    else
        trendb = 1_i4b
    end if
!
    if ( present(trend2) ) then
        trend2b = trend2
    else
        trend2b = 0_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
    else
        win2 = 3_i4b
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( (/ m, sld2/), dim=2_i4b )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    overlap2 = false
!
    if ( present(overlap) ) then
        overlap2 = overlap
    end if
!
    smooth = false
!
    if ( present(nsmooth) ) then
!
        if ( nsmooth/=0_i4b ) then
!
!           CHECK THE INPUT VALUE FOR THE LENGTH OF THE DANIELL FILTER.
!
            if ( nsmooth<3_i4b .or. nsmooth>nf )      &
            call merror( name_proc//tseries_error67 )
!
            if ( (nsmooth/2_i4b)*2_i4b==nsmooth )   &
            call merror( name_proc//tseries_error61 )
!
            smooth = true
!
        end if
!
    end if
!
!   FIND THE FIRST INTEGER GREATER THAN OR EQUAL TO n THAT IS EVENLY
!   DIVISIBLE BY l.
!
    i = mod( n, l )
!
    if ( i==0_i4b ) then
        n2 = n
    else
        n2 = n + l - i
    end if
!
!   DETERMINE THE NUMBER OF SEGMENTS.
!
    if ( overlap2 ) then
        nseg = (2*n2/l) - 1_i4b
        step = ld2
    else
        nseg = n2/l
        step = l
    end if
!
    if ( out_dof ) then
!
        if ( nseg/=1_i4b .and. overlap2 .and. win2/=1_i4b .and. win2/=3_i4b  )    &
        call merror( name_proc//tseries_error58 )
!
    end if
!
!   CALCULATE DATA WINDOW IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:l)  = data_window( l, win2, taperp=taperp )
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:l), wk(:l) )
!
    else
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:l)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( l, stnd )
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trendb>=1_i4b .and. trendb<=3_i4b ) then
        call detrend_rm( mat(:m,:n), trendb )
    end if
!
!   ZERO OUTPUT PSD ESTIMATES FOR SUMMING OVER THE m SEGMENTS.
!
    psmat(:m,:nf) = zero
!
!   COMPUTE PSD ESTIMATES BY SUMMING OVER THE m SEGMENTS.
!
#ifdef _OPENMP
    i1 = omp_get_num_procs()
    i2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               i1>1_i4b                     .and.      &
               i2>1_i4b                     .and.      &
               nseg>=i2
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(k,i,i1,i2,l2,iok,psmatb,seg,cwk)
!
!       ALLOCATE WORK MATRIX.
!
        allocate( psmatb(m,nf), seg(m,sl), cwk(m,nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ZERO PARTIAL OUTPUT PSD ESTIMATES FOR SUMMING OVER EACH THREAD.
!
        psmatb(:m,:nf) = zero
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(:m,l+1_i4b:sl) = zero
        end if    
!
!$OMP DO SCHEDULE(STATIC)
!
        do k=1_i4b, nseg
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg(:m,:l2) = mat(:m,i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENT IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rm( seg(:m,:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:m,:l2) = seg(:m,:l2)*spread( wk(:l2), dim=1, ncopies=m )
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(:m,l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT OF THE TIME SEGMENT.
!
            call real_fft( seg(:m,:sl), cwk(:m,:nf), true )
!
!           UPDATE POWER SPECTRUM ESIMATES.
!
            psmatb(:m,:sld2) = psmatb(:m,:sld2) + real(cwk(:m,:sld2),stnd)**2 + aimag(cwk(:m,:sld2))**2
            psmatb(:m,nf)    = psmatb(:m,nf)    + real(cwk(:m,nf),stnd)**2
!
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (updatepsd2)
        psmat(:m,:nf) = psmat(:m,:nf) + psmatb(:m,:nf)
!$OMP END CRITICAL (updatepsd2)
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( psmatb, seg, cwk )
!
!$OMP END PARALLEL
!
    else
!
#endif
!
!       ALLOCATE WORK MATRICES.
!
        allocate( seg(m,sl), cwk(m,nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(:m,l+1_i4b:sl) = zero
        end if    
!
        do k=1_i4b, nseg
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg(:m,:l2) = mat(:m,i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENT IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rm( seg(:m,:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:m,:l2) = seg(:m,:l2)*spread( wk(:l2), dim=1, ncopies=m )
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(:m,l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT FOR THE SERIES.
!
            call real_fft( seg(:m,:sl), cwk(:m,:nf), true )
!
!           UPDATE POWER SPECTRUM ESIMATES.
!
            psmat(:m,:sld2) = psmat(:m,:sld2) + real(cwk(:m,:sld2),stnd)**2 + aimag(cwk(:m,:sld2))**2
            psmat(:m,nf)    = psmat(:m,nf)    + real(cwk(:m,nf),stnd)**2
!
        end do
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( seg, cwk )
!
#ifdef _OPENMP
!
    end if
!
#endif
!
!   NORMALIZE THE POWER SPECTRA ESTIMATES.
!
    c1 = one/( sumw*real(l*nseg, stnd) )
!
    if ( normpsd2 ) then
!
        c2 = (real( l, stnd)/pi)*c1
        i1 = 1_i4b
        i2 = nf
!
    else
        c2 = two*c1
        i1 = 2_i4b
        i2 = sld2
!
        psmat(:m,1_i4b) = c1*psmat(:m,1_i4b)
        psmat(:m,nf)    = c1*psmat(:m,nf)
!
    end if
!
    psmat(:m,i1:i2) = c2*psmat(:m,i1:i2)
!
!   SMOOTH THE POWER SPECTRA ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
        call dan_filter_rm( psmat(:m,:nf), nsmooth, sym=one )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( sl, stnd )
        freq(:nf) = arth( zero, c1, nf )
    end if
!
!   OUTPUT DEGREES OF FREEDOM, CONFIDENCE LIMIT FACTORS
!   OR BANDWIDTH IF REQUIRED.
!
    if ( out_dof ) then
!
        allocate( edof2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
        i = sl - l
!
        edof2(:nf) = estim_dof2( wk(:l), i, win=win, nsmooth=nsmooth, nseg=nseg, overlap=overlap )
!
        if ( present(conlwr) .or. present(conupr) ) then
            call comp_conflim_rv( edof2(:nf), probtest=probtest, conlwr=conlwr, conupr=conupr )
        end if
!
        if ( present(edof) ) then
            edof(:nf) = edof2(:nf)
        end if
!
        if ( present(bandwidth) ) then
            bandwidth(:nf) = edof2(:nf)/( two*real( n, stnd ) )
        end if
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( edof2 )
!
    end if
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE power_spctrm2_rm
! __________________________________
!
    end subroutine power_spctrm2_rm
!
! =========================================================================================
!
    subroutine cross_spctrm2_rv( vec, vec2, l, psvec, psvec2, phase, coher, freq, edof, bandwidth,   &
                                 conlwr, conupr, testcoher, ampli, co_spect, quad_spect, prob_coher, &
                                 initfft, overlap, normpsd, nsmooth, trend, trend2, win, taperp, l0, &
                                 probtest )
!
! Purpose
! _______
!
!   Subroutine CROSS_SPCTRM2 computes Fast Fourier Transform (FFT) estimates
!   of the power and cross spectra of two real time series.
!
!   The Power Spectral Density (PSD) and Cross Spectral Density (CSD) estimates are returned
!   in units which are the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the first real time series for which the power and cross spectra
!                 must be estimated.
!                 If TREND=1, 2 or 3,  VEC is used as workspace and is transformed.
!
!                 Size(VEC) must be greater or equal to 4.
!
!   VEC2          (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the second real time series for which the power and cross spectra
!                 must be estimated.
!                 If TREND=1, 2 or 3,  VEC2 is used as workspace and is transformed.
!
!                 VEC2 must verify:  size(VEC2) = size(VEC).
!
!   L             (INPUT) integer(i4b)
!                 On entry, an integer used to segment the time series. L is the
!                 length of the segments. L must be a positive even integer, less
!                 or equal to size(VEC), but greater or equal to 4.
!
!                 Spectral computations are at (L/2)+1 frequencies if the optional
!                 argument L0 is absent and are at ((L+L0)/2)+1 frequencies if L0 is
!                 present (L0 is the number of zeros added to each segment).
!
!                 Suggested values for L+L0 are 16, 32, 64 or 128 (e.g. an integer power of two,
!                 in order to speed the computations).
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = ((L+L0)/2) + 1 .
!
!   PSVEC2        (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC2.  
!
!                 PSVEC2 must verify:  size(PSVEC2) = ((L+L0)/2) + 1 .
!
!   PHASE         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the phase of the cross spectrum, given in fractions
!                 of a circle (e.g. on the closed interval (0,1) ).
!
!                 PHASE must verify:  size(PHASE) = ((L+L0)/2) + 1 .
!
!   COHER         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the squared coherency  estimates for all frequencies.
!
!                 COHER must verify:  size(COHER) = ((L+L0)/2) + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/(L+L0)
!                 for i=1,2, ... , ((L+L0)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = (L+L0)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!                 EDOF must verify:  size(EDOF) = ((L+L0)/2) + 1 .
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the bandwidth of the power spectrum estimates.
!
!                 BANDWIDTH must verify:  size(BANDWIDTH) = ((L+L0)/2) + 1 .
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the
!                 PSVEC(:) and PSVEC2(:) arguments) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!                 CONLWR must verify:  size(CONLWR) = ((L+L0)/2) + 1 .
!
!                 CONUPR must verify:  size(CONUPR) = ((L+L0)/2) + 1 .
!   
!   TESTCOHER     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, this argument specifies the critical value for testing the null
!                 hypothesis that the squared coherency is zero at the PROBTEST * 100% significance
!                 level (e.g. elements of COHER(:) less than TESTCOHER(:) should be regarded as not
!                 significantly different from zero at the PROBTEST * 100% significance level).
!
!                 TESTCOHER must verify:  size(TESTCOHER) = ((L+L0)/2) + 1 .
!
!   AMPLI         (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the cross-amplitude spectrum.  
!
!                 AMPLI must verify:  size(AMPLI) = ((L+L0)/2) + 1 .
!
!   CO_SPECT      (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the co-spectrum (e.g. the real part of cross-spectrum).  
!
!                 CO_SPECT must verify:  size(CO_SPECT) = ((L+L0)/2) + 1 .
!
!   QUAD_SPECT    (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the quadrature spectrum (e.g. the imaginary part of
!                 cross-spectrum with a minus sign).  
!
!                 QUAD_SPECT must verify:  size(QUAD_SPECT) = ((L+L0)/2) + 1 .
!
!   PROB_COHER    (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the probabilities that the computed sample squared coherencies
!                 came from an ergodic stationary bivariate process with (corresponding)
!                 squared coherencies equal to zero.
!
!                 PROB_COHER  must verify:  size(PROB_COHER) = ((L+L0)/2)+1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine CROSS_SPCTRM2 in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine CROSS_SPCTRM2. This call to INITFFT must have the following form: 
!
!                       call init_fft( (L+L0)/2 )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   CROSS_SPCTRM2 and a call to END_FFT is also done before leaving
!                   subroutine CROSS_SPCTRM2.
!
!                 The default is INITFFT=true .
!
!   OVERLAP       (INPUT, OPTIONAL) logical(lgl)
!                 If:
!
!                 - OVERLAP = false, the subroutine segments the data 
!                   without any overlapping.
!                 - OVERLAP=true, the subroutine overlaps the segments
!                   by one half of their length (which is equal to L).
!
!                 In both cases, zeros are eventually added to each segment (if argument L0 is present)
!                 and each segment will be FFT'd, and the resulting periodograms
!                 will averaged together to obtain a Power Spectrum Density estimate at the
!                 ((L+L0)/2)+1 frequencies.
!
!                 The default is OVERLAP=false .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the power and cross spectra estimates are normalized
!                   in such a way that the total area under the power spectrum is equal to the variance
!                   of the time series VEC and VEC2.
!                 - NORMPSD = false, the sum of the PSD estimates
!                   (e.g. sum(PSVEC(2:)) and sum(PSVEC2(2:)) ) is equal to the variance of the corresponding
!                   time series.
!
!                 The default is NORMPSD=true .
!
!   NSMOOTH       (INPUT, OPTIONAL) integer(i4b)
!                 if NSMOOTH is used, the PSD and CSD estimates are computed by smoothing 
!                 the periodogram with Daniell weights (e.g. a simple moving average).
!
!                 On entry, NSMOOTH gives the length of the Daniell filter to be applied.
!
!                 Setting NSMOOTH=0 on entry is equivalent to omit the optional argument
!                 NSMOOTH. Otherwise, NSMOOTH must be odd, greater than 2 and less or
!                 equal to (L+L0)/2+1 .
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The mean of the two time series is removed before computing the spectra
!                 - TREND=+2 The drift from the two time series is removed before computing the spectra
!                 - TREND=+3 The least-squares line from the two time series is removed before
!                   computing the spectra.
!
!                 For other values of TREND nothing is done before estimating the power and cross spectra.
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   TREND2        (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND2=+1 The mean of the time segment is removed before computing the cross-spectrum
!                   on this segment.
!                 - TREND2=+2 The drift from the time segment is removed before computing the cross-spectrum
!                   on this segment.
!                 - TREND2=+3 The least-squares line from the time segment is removed before
!                   computing the cross-spectrum on this segment.
!
!                 For other values of TREND2 nothing is done before estimating the cross-spectrum
!                 on each segment.
!
!                 The default is TREND2=0, e.g. nothing is done before estimating the power spectrum
!                 on each segment.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power and cross spectra. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!
!   L0            (INPUT, OPTIONAL) integer(i4b)
!                 The number of zeros added to each time segment in order to obtain more finely
!                 spaced spectral estimates. L+L0 must be a positive even integer.
!
!                 The default is L0=0, e.g. no zeros are added to each time segment.
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ) and the critical value for
!                 testing the null hypothesis that the squared coherency is zero (e.g.
!                 the TESTCOHER optional argument).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the two time series (e.g. TREND=1,2,3), the series
!   are padded with zero on the right such that the length of the resulting two time series is 
!   evenly divisible by L (a positive even integer). The length, N, of these resulting time
!   series is the first integer greater than or equal to size(VEC) which is evenly divisible
!   by L. If size(VEC) is not evenly divisible by L, N is equal to size(VEC)+L-mod(size(VEC),L).
!   
!   Optionally, the mean or the trend may also be removed from each time segment (e.g. TREND2=1,2,3).
!   Optionally, zeros may be added to each time segment (e.g. the optional arguemnt L0) if more
!   finely spaced spectral esimates are desired.
!   
!   The stability of the power and cross spectra  estimates depends on the averaging process.
!   That is, the greater the number of segments ( N/L if OVERLAP=false and (2N/L)-1 if OVERLAP=true),
!   the more stable the resulting power and cross spectra estimates.
!   
!   Optionally, these power and cross spectra estimates may then be smoothed again
!   in the frequency domain by modified Daniell filters (e.g. if argument NSMOOTH is used).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, assert_eq, arth
    use Reals_Constants,   only : zero, one, two, pi, twopi
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error15, tseries_error17,    &
                                  tseries_error18, tseries_error24, tseries_error25, tseries_error58,   &
                                  tseries_error59, tseries_error61, tseries_error67
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : probf2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:) :: vec, vec2
    real(stnd), intent(out),    dimension(:) :: psvec, psvec2, phase, coher
!
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out), dimension(:), optional :: freq, ampli, co_spect, quad_spect, prob_coher, &
                                                       edof, bandwidth, conlwr, conupr, testcoher
!
    integer(i4b),  intent(in)                         :: l
    integer(i4b),  intent(in),               optional :: trend, trend2, win, l0, nsmooth
!
    logical(lgl),  intent(in), optional :: normpsd, initfft, overlap
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, n2, nf, k, l2, ld2, m, win2, trendb, trend2b,      &
                    i, i1, i2, step, sl, sld2
    integer      :: iok
!
    real(stnd)                            :: c1, c2, sumw, df, df2
    real(stnd), dimension(l)              :: wk
    real(stnd), dimension(:), allocatable :: seg1, seg2, edof2
!
    complex(stnd), dimension(:), allocatable  :: cwk1, cwk2
!
    logical(lgl)  :: normpsd2, initfft2, overlap2, smooth, out_freq, out_ampli,  &
                     out_co_spect, out_quad_spect, out_prob_coher, out_dof
!
#ifdef _OPENMP
!
    real(stnd), dimension(:), allocatable :: psvecb, psvec2b, rcospect, icospect
!
    logical  :: test_par
!
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='cross_spctrm2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n =  assert_eq( int(size(vec),i4b) ,          &
                    int(size(vec2),i4b) ,         &
                    name_proc )
!
    if ( l<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    if ( l>n )     &
    call merror( name_proc//tseries_error17 )
!
!   CHECK IF l IS A POSITIVE EVEN INTEGER.
!
    ld2 = l/2_i4b
!
    if ( l/=2_i4b*ld2 )     &
    call merror( name_proc//tseries_error18 )
!
!   DETERMINE THE LENGTH OF THE SEGMENTS.
!
    if ( present(l0) ) then
!
        if ( l0<0_i4b )     &
        call merror( name_proc//tseries_error24 )
!
        sl = l + l0
!
!       CHECK IF sl IS A POSITIVE EVEN INTEGER.
!
        sld2 = sl/2_i4b
!
        if ( sl/=2_i4b*sld2 )     &
        call merror( name_proc//tseries_error25 )
!
    else
!
        sl   = l
        sld2 = ld2
!
    end if
!
    nf  = sld2 + 1_i4b
    call assert( logical(int(size(psvec),i4b)==nf,lgl),   &
                 logical(int(size(psvec2),i4b)==nf,lgl),  &
                 logical(int(size(phase),i4b)==nf,lgl),   &
                 logical(int(size(coher),i4b)==nf,lgl),   &
                 name_proc )
!
    out_freq = present( freq )
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_ampli = present( ampli )
    if ( out_ampli ) then
        call assert( logical(int(size(ampli),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_co_spect = present( co_spect )
    if ( out_co_spect ) then
        call assert( logical(int(size(co_spect),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_quad_spect = present( quad_spect )
    if ( out_quad_spect ) then
        call assert( logical(int(size(quad_spect),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_prob_coher = present( prob_coher )
    if ( out_prob_coher ) then
        call assert( logical(int(size(prob_coher),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth ) .or. out_prob_coher   &
              .or. present( conlwr ) .or. present( conupr ) .or. present( testcoher ) 
!
    if ( out_dof ) then
!
        if ( present(probtest) ) then
            if ( zero>=probtest .or. probtest>=one ) then
                call merror( name_proc//tseries_error59 )
            end if
        end if
!
        if ( present( edof ) ) then
            call assert( logical(int(size(edof),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( bandwidth ) ) then
            call assert( logical(int(size(bandwidth),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conlwr ) ) then
            call assert( logical(int(size(conlwr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conupr ) ) then
            call assert( logical(int(size(conupr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
    end if
!
    if ( present(trend) ) then
        trendb = trend
    else
        trendb = 1_i4b
    end if
!
    if ( present(trend2) ) then
        trend2b = trend2
    else
        trend2b = 0_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
    else
        win2 = 3_i4b
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( sld2 )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    overlap2 = false
!
    if ( present(overlap) ) then
        overlap2 = overlap
    end if
!
    smooth = false
!
    if ( present(nsmooth) ) then
!
        if ( nsmooth/=0_i4b ) then
!
!           CHECK THE INPUT VALUE FOR THE LENGTH OF THE DANIELL FILTER.
!
            if ( nsmooth<3_i4b .or. nsmooth>nf )      &
            call merror( name_proc//tseries_error67 )
!
            if ( (nsmooth/2_i4b)*2_i4b==nsmooth )   &
            call merror( name_proc//tseries_error61 )
!
            smooth = true
!
        end if
!
    end if
!
!   FIND THE FIRST INTEGER GREATER THAN OR EQUAL TO n THAT IS EVENLY
!   DIVISIBLE BY l.
!
    i = mod( n, l )
!
    if ( i==0_i4b ) then
        n2 = n
    else
        n2 = n + l - i
    end if
!
!   DETERMINE THE NUMBER OF SEGMENTS.
!
    if ( overlap2 ) then
        m    = (2*n2/l) - 1_i4b
        step = ld2
    else
        m    = n2/l
        step = l
    end if
!
    if ( out_dof ) then
!
        if ( m/=1_i4b .and. overlap2 .and. win2/=1_i4b .and. win2/=3_i4b  )    &
        call merror( name_proc//tseries_error58 )
!
    end if
!
!   CALCULATE DATA WINDOW IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:l)  = data_window( l, win2, taperp=taperp )
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:l), wk(:l) )
!
    else
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:l)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( l, stnd )
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trendb>=1_i4b .and. trendb<=3_i4b ) then
!
        call detrend_rv( vec(:n), trendb )
        call detrend_rv( vec2(:n), trendb )
!
    end if
!
!   ZERO OUTPUT POWER AND CROSS SPECTRA ESTIMATES FOR SUMMING OVER THE m SEGMENTS.
!
    psvec(:nf)  = zero
    psvec2(:nf) = zero
!
!   USE THE ARRAYS coher AND phase TO STORE THE REAL AND IMAGINARY PARTS OF THE
!   CROSS SPECTRUM, RESPECTIVELY.
!
    coher(:nf)  = zero
    phase(:nf)  = zero
!
!   COMPUTE PSD AND CSD ESTIMATES BY SUMMING OVER THE m SEGMENTS.
!
#ifdef _OPENMP
    i1 = omp_get_num_procs()
    i2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               i1>1_i4b                     .and.      &
               i2>1_i4b                     .and.      &
               m>=i2
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(k,i,i1,i2,l2,iok,psvecb,psvec2b,seg1,seg2,cwk1,cwk2,rcospect,icospect)
!
!       ALLOCATE WORK VECTORS.
!
        allocate( psvecb(nf), psvec2b(nf), rcospect(nf), icospect(nf),     &
                  seg1(sl), seg2(sl), cwk1(nf), cwk2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ZERO PARTIAL OUTPUT PSD AND CSD ESTIMATES FOR SUMMING OVER EACH THREAD.
!
        psvecb(:nf)  = zero
        psvec2b(:nf) = zero
!
        rcospect(:nf) = zero
        icospect(:nf) = zero
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg1(l+1_i4b:sl) = zero
            seg2(l+1_i4b:sl) = zero
        end if    
!
!$OMP DO SCHEDULE(STATIC)
!
        do k=1_i4b, m
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg1(:l2) = vec(i1:i2)
            seg2(:l2) = vec2(i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENTS IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rv( seg1(:l2), trend2b )
                call detrend_rv( seg2(:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW TO THE SERIES.
!
            if ( win2/=2_i4b ) then
                seg1(:l2) = seg1(:l2)*wk(:l2)
                seg2(:l2) = seg2(:l2)*wk(:l2)
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg1(l2+1_i4b:l) = zero
                seg2(l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT FOR THE FIRST SERIES.
!
            call real_fft( seg1(:sl), cwk1(:nf), true)
!
!           COMPUTE FFT FOR THE SECOND SERIES.
!
            call real_fft( seg2(:sl), cwk2(:nf), true)
!
!           UPDATE SPECTRAL ESTIMATES.
!
            psvecb(:sld2) = psvecb(:sld2) + real(cwk1(:sld2),stnd)**2 + aimag(cwk1(:sld2))**2
            psvecb(nf)    = psvecb(nf)    + real(cwk1(nf),stnd)**2
!
            psvec2b(:sld2) = psvec2b(:sld2) + real(cwk2(:sld2),stnd)**2 + aimag(cwk2(:sld2))**2
            psvec2b(nf)    = psvec2b(nf)    + real(cwk2(nf),stnd)**2
!
!           UPDATE CO-SPECTRUM ESTIMATES.
!
            rcospect(:sld2)  = rcospect(:sld2)  + real(cwk1(:sld2),stnd)*real(cwk2(:sld2),stnd) +   &
                                                  aimag(cwk1(:sld2))*aimag(cwk2(:sld2))
            rcospect(nf)     = rcospect(nf)     + real(cwk1(nf),stnd)*real(cwk2(nf),stnd)
!
!           UPDATE QUADRATURE-SPECTRUM ESTIMATES.
!
            icospect(:sld2)  = icospect(:sld2)  + aimag(cwk1(:sld2))*real(cwk2(:sld2),stnd)     -   &
                                                  real(cwk1(:sld2),stnd)*aimag(cwk2(:sld2))
!
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (updatepscd)
!
        psvec(:nf)    = psvec(:nf) + psvecb(:nf)
        psvec2(:nf)   = psvec2(:nf) + psvec2b(:nf)
!
        coher(:nf)    = coher(:nf)   + rcospect(:nf) 
        phase(:sld2)  = phase(:sld2) + icospect(:sld2)
!
!$OMP END CRITICAL (updatepscd)
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( psvecb, psvec2b, seg1, seg2, cwk1, cwk2, rcospect, icospect )
!
!$OMP END PARALLEL
!
!       REALLOCATE WORK VECTORS.
!
        allocate( seg1(nf), seg2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
    else
!
#endif
!
!       ALLOCATE WORK VECTORS.
!
        allocate( seg1(sl), seg2(sl), cwk1(nf), cwk2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg1(l+1_i4b:sl) = zero
            seg2(l+1_i4b:sl) = zero
        end if    
!
        do k=1_i4b, m
!
            i = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg1(:l2) = vec(i1:i2)
            seg2(:l2) = vec2(i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENTS IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rv( seg1(:l2), trend2b )
                call detrend_rv( seg2(:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW TO THE SERIES.
!
            if ( win2/=2_i4b ) then
                seg1(:l2) = seg1(:l2)*wk(:l2)
                seg2(:l2) = seg2(:l2)*wk(:l2)
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg1(l2+1_i4b:l) = zero
                seg2(l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT FOR THE FIRST SERIES.
!
            call real_fft( seg1(:sl), cwk1(:nf), true)
!
!           COMPUTE FFT FOR THE SECOND SERIES.
!
            call real_fft( seg2(:sl), cwk2(:nf), true)
!
!           UPDATE POWER SPECTRUM ESIMATES.
!
            psvec(:sld2) = psvec(:sld2) + real(cwk1(:sld2),stnd)**2 + aimag(cwk1(:sld2))**2
            psvec(nf)    = psvec(nf)    + real(cwk1(nf),stnd)**2
!
            psvec2(:sld2) = psvec2(:sld2) + real(cwk2(:sld2),stnd)**2 + aimag(cwk2(:sld2))**2
            psvec2(nf)    = psvec2(nf)    + real(cwk2(nf),stnd)**2
!
!           UPDATE CO-SPECTRUM ESIMATES.
!
            coher(:sld2)  = coher(:sld2)  + real(cwk1(:sld2),stnd)*real(cwk2(:sld2),stnd) +   &
                                            aimag(cwk1(:sld2))*aimag(cwk2(:sld2))
            coher(nf)    = coher(nf)    + real(cwk1(nf),stnd)*real(cwk2(nf),stnd)
!
!           UPDATE QUADRATURE-SPECTRUM ESIMATES.
!
            phase(:sld2)  = phase(:sld2)  + aimag(cwk1(:sld2))*real(cwk2(:sld2),stnd)    -   &
                                            real(cwk1(:sld2),stnd)*aimag(cwk2(:sld2))
!
        end do
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( cwk1, cwk2 )
!        
#ifdef _OPENMP
!
    end if
!
#endif
!
!   NORMALIZE THE POWER AND CROSS SPECTRA ESTIMATES.
!
    c1 = one/( sumw*real(l*m, stnd) )
!
    if ( normpsd2 ) then
!
        c2 = (real( l, stnd)/pi)*c1
        i1 = 1_i4b
        i2 = nf
!
    else
!
        c2 = two*c1
        i1 = 2_i4b
        i2 = sld2
!
        psvec(1_i4b)  = c1*psvec(1_i4b)
        psvec2(1_i4b) = c1*psvec2(1_i4b)
!
        seg1(1_i4b) = c1*coher(1_i4b)
        seg2(1_i4b) = c1*phase(1_i4b)
!
        psvec(nf)  = c1*psvec(nf)
        psvec2(nf) = c1*psvec2(nf)
!
        seg1(nf)  = c1*coher(nf)
        seg2(nf)  = c1*phase(nf)
!
    end if
!
    psvec(i1:i2)  = c2*psvec(i1:i2)
    psvec2(i1:i2) = c2*psvec2(i1:i2)
!
    seg1(i1:i2) = c2*coher(i1:i2)
    seg2(i1:i2) = c2*phase(i1:i2)
!
!   SMOOTH THE POWER AND CROSS SPECTRA ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
!       SMOOTH POWER SPECTRUM OF THE FIRST SERIES.
!
        call dan_filter_rv( psvec(:nf), nsmooth, sym=one )
!
!       SMOOTH POWER SPECTRUM OF THE SECOND SERIES.
!
        call dan_filter_rv( psvec2(:nf), nsmooth, sym=one )
!
!       SMOOTH CO-SPECTRUM.
!
        call dan_filter_rv( seg1(:nf), nsmooth, sym=one )
!
!       SMOOTH QUADRATURE-SPECTRUM.
!
        call dan_filter_rv( seg2(:nf), nsmooth, sym=-one )
!
    end if
!
    if ( out_co_spect ) then
!
!       OUTPUT THE CO-SPECTRUM IF REQUIRED.
!
        co_spect(:nf) = seg1(:nf)
!
    end if
!
    if ( out_quad_spect ) then
!
!       OUTPUT THE QUADRATURE-SPECTRUM IF REQUIRED.
!
        quad_spect(:nf) = -seg2(:nf)
!
    end if
!
!   COMPUTE PHASE OF THE CROSS-SPECTRUM.
!
    where ( seg1(:nf)==zero .and. seg2(:nf)==zero )
        phase(:nf) = zero
    elsewhere
        coher(:nf)   = atan2( seg2(:nf), seg1(:nf) )
        phase(:nf) = (one/twopi)*coher(:nf) + merge( zero, one, coher(:nf)>=zero ) 
    end where
!
!   COMPUTE MAGNITUDE OF THE CROSS-SPECTRUM.
!
    coher(:nf) = seg1(:nf)**2 + seg2(:nf)**2
    seg2(:nf)  = coher(:nf)
!
!   COMPUTE THE SQUARED COHERENCE FOR ALL FREQUENCIES.
!
    where ( psvec(:nf)/=zero .and. psvec2(:nf)/=zero )
        seg1(:nf)  = psvec(:nf)*psvec2(:nf)
        coher(:nf) = min( seg2(:nf)/seg1(:nf), one )
    elsewhere
        coher(:nf) = zero
    end where 
!
!   COMPUTE AND OUTPUT AMPLITUDE OF THE CROSS-SPECTRUM IF REQUIRED.
!
    if ( out_ampli ) then
!
        ampli(:nf) = sqrt( seg2(:nf) )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( sl, stnd )
        freq(:nf) = arth(zero, c1, nf)
    end if
!
!   OUTPUT DEGREES OF FREEDOM, BANDWIDTH, CONFIDENCE LIMIT FACTORS 
!   OR PROBABILITIES FOR SQUARED COHERENCIES IF REQUIRED.
!
    if ( out_dof ) then
!
        allocate( edof2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
        i = sl - l
!
        edof2(:nf) = estim_dof2( wk(:l), i, win=win, nsmooth=nsmooth, nseg=m, overlap=overlap )
!
        if ( present(edof) ) then
            edof(:nf) = edof2(:nf)
        end if
!
        if ( present(bandwidth) ) then
            bandwidth(:nf) = edof2(:nf)/( two*real( n, stnd ) )
        end if
!
        if ( present(conlwr) .or. present(conupr) .or. present(testcoher) ) then
            call comp_conflim_rv( edof2(:nf), probtest=probtest, conlwr=conlwr,       &
                                  conupr=conupr, testcoher=testcoher )
        end if
!
!       COMPUTE AND OUTPUT SIGNIFICANCE PROBABILITY FOR SQUARED COHERENCIES.
!
        if ( out_prob_coher ) then
!
            sumw = nan()
!
            do i = 1_i4b, nf
!            
                df = edof2(i)  
!            
                if ( df>two ) then
!            
                    c1 = coher(i)
!            
                    if ( c1/=one ) then
!            
                        c2 = ( (df/two - one)*c1 )/( one - c1 )
                        df2 = df - two
!            
                        prob_coher(i) = probf2( c2, two, df2, true )
!            
                    else
                        prob_coher(i) = zero
                    end if 
!
                else
!
                    prob_coher(i)    = sumw
!
                end if
!
            end do
!
        end if
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( edof2 )
!
    end if
!
!   DEALLOCATE WORK ARRAYS.
!
    deallocate( seg1, seg2 )
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE cross_spctrm2_rv
! __________________________________
!
    end subroutine cross_spctrm2_rv
!
! =========================================================================================
!
    subroutine cross_spctrm2_rm( vec, mat, l, psvec, psmat, phase, coher, freq, edof, bandwidth,     &
                                 conlwr, conupr, testcoher, ampli, co_spect, quad_spect, prob_coher, &
                                 initfft, overlap, normpsd, nsmooth, trend, trend2, win, taperp, l0, &
                                 probtest )
!
! Purpose
! _______
!
!   Subroutine CROSS_SPCTRM2 computes Fast Fourier Transform (FFT) estimates
!   of the power and cross spectra of the real time series, VEC, and the multi-channel
!   real time series MAT.
!
!   The Power Spectral Density (PSD) and Cross Spectral Density (CSD) estimates are returned
!   in units which are the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the real time series for which the power and cross spectra
!                 must be estimated.
!                 If TREND=1, 2 or 3,  VEC is used as workspace and is transformed.
!
!                 Size(VEC) must be greater or equal to 4.
!
!   MAT           (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                 On entry, the multi-channel real time series for which the power 
!                 and cross spectra must be estimated. Each row of MAT is a real time series.
!                 If TREND=1, 2 or 3,  MAT is used as workspace and is transformed.
!
!                 The shape of MAT must verify:  size(MAT,2) = size(VEC).
!
!   L             (INPUT) integer(i4b)
!                 On entry, an integer used to segment the time series. L is the
!                 length of the segments. L must be a positive even integer, less
!                 or equal to size(VEC), but greater or equal to 4.
!
!                 Spectral computations are at (L/2)+1 frequencies if the optional
!                 argument L0 is absent and are at ((L+L0)/2)+1 frequencies if L0 is
!                 present (L0 is the number of zeros added to each segment).
!
!                 Suggested values for L+L0 are 16, 32, 64 or 128 (e.g. an integer power of two,
!                 in order to speed the computations).
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = ((L+L0)/2) + 1 .
!
!   PSMAT         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the Power Spectral Density (PSD) estimates of each row of MAT.
!
!                 The shape of PSMAT must verify:
!
!                 - size(PSMAT,1) =  size(MAT,1) ;
!                 - size(PSMAT,2) = ((L+L0)/2) + 1 .
!
!   PHASE         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the phase of the cross spectrum, given in fractions
!                 of a circle (e.g. on the closed interval (0,1) ).
!
!                 The shape of PHASE must verify:
!
!                 - size(PHASE,1) =  size(MAT,1) ;
!                 - size(PHASE,2) = ((L+L0)/2) + 1 .
!
!   COHER         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the squared coherency  estimates for all frequencies.
!
!                 The shape of COHER must verify:
!
!                 - size(COHER,1) =  size(MAT,1) ;
!                 - size(COHER,2) = ((L+L0)/2) + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/(L+L0)
!                 for i=1,2, ... , ((L+L0)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = (L+L0)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!                 EDOF must verify:  size(EDOF) = ((L+L0)/2) + 1 .
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, the bandwidth of the power spectrum estimates.
!
!                 BANDWIDTH must verify:  size(BANDWIDTH) = ((L+L0)/2) + 1 .
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the PSVEC(:) and
!                 PSMAT(:,:) arguments) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!                 CONLWR must verify:  size(CONLWR) = ((L+L0)/2) + 1 .
!
!                 CONUPR must verify:  size(CONUPR) = ((L+L0)/2) + 1 .
!   
!   TESTCOHER     (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On output, this argument specifies the critical value for testing the null
!                 hypothesis that the squared coherency is zero at the PROBTEST * 100% significance
!                 level (e.g. elements of COHER(:,:) less than TESTCOHER(:) should be regarded as not
!                 significantly different from zero at the PROBTEST * 100% significance level).
!
!                 TESTCOHER must verify:  size(TESTCOHER) = ((L+L0)/2) + 1 .
!
!   AMPLI         (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the cross-amplitude spectra.
!
!                 The shape of AMPLI must verify:
!
!                 - size(AMPLI,1) =  size(MAT,1) ;
!                 - size(AMPLI,2) = ((L+L0)/2) + 1 .
!
!   CO_SPECT      (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the co-spectra (e.g. the real part of cross-spectra).  
!
!                 The shape of CO_SPECT must verify:
!
!                 - size(CO_SPECT,1) =  size(MAT,1) ;
!                 - size(CO_SPECT,2) = ((L+L0)/2) + 1 .
!
!   QUAD_SPECT    (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the quadrature spectrum (e.g. the imaginary part of
!                 cross-spectrum with a minus sign).  
!
!                 The shape of QUAD_SPECT must verify:
!
!                 - size(QUAD_SPECT,1) =  size(MAT,1) ;
!                 - size(QUAD_SPECT,2) = ((L+L0)/2) + 1 .
!
!   PROB_COHER    (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the probabilities that the computed sample squared coherencies
!                 came from an ergodic stationary bivariate process with (corresponding)
!                 squared coherencies equal to zero.
!
!                 The shape of PROB_COHER must verify:
!
!                 size(PROB_COHER,1) =  size(MAT,1) ;
!                 size(PROB_COHER,2) =  ((L+L0)/2) + 1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine CROSS_SPCTRM2 in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine CROSS_SPCTRM2. This call to INITFFT must have the following form: 
!
!                        call init_fft( (/ size(MAT,1), (L+L0)/2 /), dim=2_i4b )
!
!                 - INITFFT is set to true, the call to INIT_FFT is done inside subroutine
!                   CROSS_SPCTRM2 and a call to END_FFT is also done before leaving
!                   subroutine CROSS_SPCTRM2.
!
!                 The default is INITFFT=true .
!
!   OVERLAP       (INPUT, OPTIONAL) logical(lgl)
!                 If:
!
!                 - OVERLAP = false, the subroutine segments the data 
!                   without any overlapping.
!                 - OVERLAP=true, the subroutine overlaps the segments
!                   by one half of their length (which is equal to L).
!
!                 In both cases, zeros are eventually added to each segment (if argument L0 is present)
!                 and each segment will be FFT'd, and the resulting periodograms
!                 will averaged together to obtain a Power Spectrum Density estimate at the
!                 ((L+L0)/2)+1 frequencies.
!
!                 The default is OVERLAP=false .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the power and cross spectra estimates are
!                   normalized in such a way that the total area under the power spectra is equal to
!                   the variance of the time series contained in VEC and in each row of MAT.
!                 - NORMPSD is set to false, the sum of the PSD estimates
!                   (e.g. sum(PSVEC(2:)) and sum(PSMAT(:,2:),dim=2) ) is equal to the variance of
!                   the corresponding time series.
!
!                 The default is NORMPSD=true .
!
!   NSMOOTH       (INPUT, OPTIONAL) integer(i4b)
!                 If NSMOOTH is used, the PSD estimates are computed by smoothing 
!                 the periodogram with Daniell weights (e.g. a simple moving average).
!
!                 On entry, NSMOOTH gives the length of the Daniell filter to be applied.
!
!                 Setting NSMOOTH=0 on entry is equivalent to omit the optional argument
!                 NSMOOTH. Otherwise, NSMOOTH must be odd, greater than 2 and less or
!                 equal to (L+L0)/2+1 .
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The means of the time series are removed before computing the spectra
!                 - TREND=+2 The drifts from time series are removed before computing the spectra
!                 - TREND=+3 The least-squares lines from time series are removed before
!                   computing the spectra.
!
!                 For other values of TREND nothing is done before estimating the power and cross spectra.
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   TREND2        (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND2=+1 The mean of the time segment is removed before computing the cross-spectrum
!                   on this segment.
!                 - TREND2=+2 The drift from the time segment is removed before computing the cross-spectrum
!                   on this segment.
!                 - TREND2=+3 The least-squares line from the time segment is removed before
!                   computing the cross-spectrum on this segment.
!
!                 For other values of TREND2 nothing is done before estimating the cross-spectrum
!                 on each segment.
!
!                 The default is TREND2=0, e.g. nothing is done before estimating the power spectrum
!                 on each segment.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power and cross spectra. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!
!   L0            (INPUT, OPTIONAL) integer(i4b)
!                 The number of zeros added to each time segment in order to obtain more finely
!                 spaced spectral estimates. L+L0 must be a positive even integer.
!
!                 The default is L0=0, e.g. no zeros are added to each time segment.
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ) and the critical value for
!                 testing the null hypothesis that the squared coherency is zero (e.g.
!                 the TESTCOHER optional argument).
!
!                 PROBTEST must verify   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the series
!   are padded with zero on the right such that the length of the resulting time series is 
!   evenly divisible by L (a positive even integer). The length, N, of these resulting time
!   series is the first integer greater than or equal to size(VEC) which is evenly divisible
!   by L. If size(VEC) is not evenly divisible by L, N is equal to size(VEC)+L-mod(size(VEC),L).
!   
!   Optionally, the mean or the trend may also be removed from each time segment (e.g. TREND2=1,2,3).
!   Optionally, zeros may be added to each time segment (e.g. the optional arguemnt L0) if more
!   finely spaced spectral esimates are desired.
!   
!   The stability of the power and cross spectra  estimates depends on the averaging process.
!   That is, the greater the number of segments ( N/L if OVERLAP=false and (2N/L)-1 if OVERLAP=true),
!   the more stable the resulting power and cross spectra estimates.
!   
!   Optionally, these power and cross spectra estimates may then be smoothed again
!   in the frequency domain by modified Daniell filters (e.g. if argument NSMOOTH is used).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, assert_eq, arth
    use Reals_Constants,   only : zero, one, two, pi, twopi
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error15, tseries_error17,   &
                                  tseries_error18, tseries_error21, tseries_error24, tseries_error25,  &
                                  tseries_error58, tseries_error59, tseries_error61, tseries_error67
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : probf2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:)   :: vec
    real(stnd), intent(inout),  dimension(:,:) :: mat
    real(stnd), intent(out),    dimension(:)   :: psvec
    real(stnd), intent(out),    dimension(:,:) :: psmat, phase, coher
!
    real(stnd), intent(in),                  optional :: taperp, probtest
    real(stnd), intent(out), dimension(:),   optional :: freq, edof, bandwidth, conlwr, conupr, testcoher
    real(stnd), intent(out), dimension(:,:), optional :: ampli, co_spect, quad_spect, prob_coher
!
    integer(i4b),  intent(in)                         :: l
    integer(i4b),  intent(in),               optional :: trend, trend2, win, l0, nsmooth
!
    logical(lgl),  intent(in), optional :: normpsd, initfft, overlap
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, mp1, n, n2, nf, nseg, k, l2, ld2, win2, trendb, trend2b,   &
                    i, i1, i2, step, sl, sld2
    integer      :: iok
!
    real(stnd)                                :: c1, c2, sumw, df, df2
    real(stnd), dimension(l)                  :: wk
    real(stnd), dimension(:),   allocatable   :: temp, edof2
    real(stnd), dimension(:,:), allocatable   :: seg, magni
!
    complex(stnd), dimension(:,:), allocatable :: cwk
!
    logical(lgl)  :: normpsd2, initfft2, overlap2, smooth, out_freq,   &
                     out_ampli, out_co_spect, out_quad_spect,          &
                     out_prob_coher, out_dof
!
#ifdef _OPENMP
!
    real(stnd), dimension(:,:), allocatable :: psmatb, rcospect, icospect
!
    logical  :: test_par
!
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='cross_spctrm2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m =  assert_eq( int(size(mat,1),i4b) ,          &
                    int(size(psmat,1),i4b) ,        &
                    int(size(phase,1),i4b) ,        &
                    int(size(coher,1),i4b) ,        &
                    name_proc )
!
    if ( m<=0_i4b )     &
    call merror( name_proc//tseries_error21 )
!
    n =  assert_eq( int(size(vec),i4b) ,          &
                    int(size(mat,2),i4b) ,        &
                    name_proc )
!
    if ( l<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    if ( l>n )     &
    call merror( name_proc//tseries_error17 )
!
!   CHECK IF l IS A POSITIVE EVEN INTEGER.
!
    ld2 = l/2_i4b
!
    if ( l/=2_i4b*ld2 )     &
    call merror( name_proc//tseries_error18 )
!
!   DETERMINE THE LENGTH OF THE SEGMENTS.
!
    if ( present(l0) ) then
!
        if ( l0<0_i4b )     &
        call merror( name_proc//tseries_error24 )
!
        sl = l + l0
!
!       CHECK IF sl IS A POSITIVE EVEN INTEGER.
!
        sld2 = sl/2_i4b
!
        if ( sl/=2_i4b*sld2 )     &
        call merror( name_proc//tseries_error25  )
!
    else
!
        sl   = l
        sld2 = ld2
!
    end if
!
    nf  = sld2 + 1_i4b
!
    call assert( logical(int(size(psvec),i4b)==nf,lgl),     &
                 logical(int(size(psmat,2),i4b)==nf,lgl),   &
                 logical(int(size(phase,2),i4b)==nf,lgl),   &
                 logical(int(size(coher,2),i4b)==nf,lgl),   &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_ampli = present( ampli )
!
    if ( out_ampli ) then
        call assert( logical(int(size(ampli,1),i4b)==m,lgl),   &
                     logical(int(size(ampli,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_co_spect = present( co_spect )
!
    if ( out_co_spect ) then
        call assert( logical(int(size(co_spect,1),i4b)==m,lgl),   &
                     logical(int(size(co_spect,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_quad_spect = present( quad_spect )
!
    if ( out_quad_spect ) then
        call assert( logical(int(size(quad_spect,1),i4b)==m,lgl),   &
                     logical(int(size(quad_spect,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_prob_coher = present( prob_coher )
!
    if ( out_prob_coher ) then
        call assert( logical(int(size(prob_coher,1),i4b)==m,lgl),   &
                     logical(int(size(prob_coher,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth ) .or. out_prob_coher   &
              .or. present( conlwr ) .or. present( conupr ) .or. present( testcoher ) 
!
    if ( out_dof ) then
!
        if ( present(probtest) ) then
            if ( zero>=probtest .or. probtest>=one ) then
                call merror( name_proc//tseries_error59 )
            end if
        end if
!
        if ( present( edof ) ) then
            call assert( logical(int(size(edof),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( bandwidth ) ) then
            call assert( logical(int(size(bandwidth),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conlwr ) ) then
            call assert( logical(int(size(conlwr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
        if ( present( conupr ) ) then
            call assert( logical(int(size(conupr),i4b)==nf,lgl),    &
                         name_proc )
        end if
!
    end if
!
    if ( present(trend) ) then
        trendb = trend
    else
        trendb = 1_i4b
    end if
!
    if ( present(trend2) ) then
        trend2b = trend2
    else
        trend2b = 0_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
    else
        win2 = 3_i4b
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( (/ m, sld2/), dim=2_i4b )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    overlap2 = false
!
    if ( present(overlap) ) then
        overlap2 = overlap
    end if
!
    smooth = false
!
    if ( present(nsmooth) ) then
!
        if ( nsmooth/=0_i4b ) then
!
!           CHECK THE INPUT VALUE FOR THE LENGTH OF THE DANIELL FILTER.
!
            if ( nsmooth<3_i4b .or. nsmooth>nf )      &
            call merror( name_proc//tseries_error67 )
!
            if ( (nsmooth/2_i4b)*2_i4b==nsmooth )   &
            call merror( name_proc//tseries_error61 )
!
            smooth = true
!
        end if
!
    end if
!
!   FIND THE FIRST INTEGER GREATER THAN OR EQUAL TO n THAT IS EVENLY
!   DIVISIBLE BY l.
!
    i = mod( n, l )
!
    if ( i==0_i4b ) then
        n2 = n
    else
        n2 = n + l - i
    end if
!
!   DETERMINE THE NUMBER OF SEGMENTS.
!
    if ( overlap2 ) then
        nseg = (2*n2/l) - 1_i4b
        step = ld2
    else
        nseg = n2/l
        step = l
    end if
!
    if ( out_dof ) then
!
        if ( nseg/=1_i4b .and. overlap2 .and. win2/=1_i4b .and. win2/=3_i4b  )    &
        call merror( name_proc//tseries_error58 )
!
    end if
!
!   CALCULATE DATA WINDOW IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:l)  = data_window( l, win2, taperp=taperp )
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:l), wk(:l) )
!
    else
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:l)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( l, stnd )
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trendb>=1_i4b .and. trendb<=3_i4b ) then
        call detrend_rv( vec(:n),    trendb )
        call detrend_rm( mat(:m,:n), trendb )
    end if
!
!   ZERO OUTPUT POWER AND CROSS SPECTRA ESTIMATES FOR SUMMING OVER THE nseg SEGMENTS.
!
    psvec(:nf)    = zero
    psmat(:m,:nf) = zero
!
!   USE THE ARRAYS coher AND phase TO STORE THE REAL AND IMAGINARY PARTS OF THE
!   CROSS SPECTRUM, RESPECTIVELY.
!
    coher(:m,:nf)  = zero
    phase(:m,:nf)  = zero
!
    mp1 = m + 1_i4b
!
!   COMPUTE PSD AND CSD ESTIMATES BY SUMMING OVER THE m SEGMENTS.
!
#ifdef _OPENMP
    i1 = omp_get_num_procs()
    i2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               i1>1_i4b                     .and.      &
               i2>1_i4b                     .and.      &
               nseg>=i2
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(k,i,i1,i2,l2,iok,psmatb,seg,cwk,rcospect,icospect)
!
!       ALLOCATE WORK MATRICES.
!
        allocate( seg(mp1,sl), cwk(mp1,nf), psmatb(mp1,nf), rcospect(m,nf), icospect(m,nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ZERO PARTIAL OUTPUT PSD AND CSD ESTIMATES FOR SUMMING OVER EACH THREAD.
!
        psmatb(:mp1,:nf) = zero
!
        rcospect(:m,:nf) = zero
        icospect(:m,:nf) = zero
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(:mp1,l+1_i4b:sl) = zero
        end if    
!
!$OMP DO SCHEDULE(STATIC)
!
        do k=1_i4b, nseg
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg(1_i4b,:l2)     = vec(i1:i2)
            seg(2_i4b:mp1,:l2) = mat(:m,i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENTS IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rm( seg(:mp1,:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:mp1,:l2) = seg(:mp1,:l2)*spread( wk(:l2), dim=1, ncopies=mp1 )
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(:mp1,l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT FOR THE SERIES.
!
            call real_fft( seg(:mp1,:sl), cwk(:mp1,:nf), true )
!
!           UPDATE POWER SPECTRUM ESIMATES.
!
            psmatb(:mp1,:sld2) = psmatb(:mp1,:sld2) + real(cwk(:mp1,:sld2),stnd)**2 + aimag(cwk(:mp1,:sld2))**2
            psmatb(:mp1,nf)    = psmatb(:mp1,nf)    + real(cwk(:mp1,nf),stnd)**2
!
!           UPDATE CO-SPECTRUM ESIMATES.
!
            rcospect(:m,:sld2)  = rcospect(:m,:sld2)                                                  +   &
                  spread(real(cwk(1_i4b,:sld2),stnd),dim=1,ncopies=m)*real(cwk(2_i4b:mp1,:sld2),stnd) +   &
                  spread(aimag(cwk(1_i4b,:sld2)),dim=1,ncopies=m)*aimag(cwk(2_i4b:mp1,:sld2))
            rcospect(:m,nf)    = rcospect(:m,nf) + real(cwk(1_i4b,nf),stnd)*real(cwk(2_i4b:mp1,nf),stnd)
!
!           UPDATE QUADRATURE-SPECTRUM ESIMATES.
!
            icospect(:m,:sld2)  = icospect(:m,:sld2)                                              +   &
                  spread(aimag(cwk(1_i4b,:sld2)),dim=1,ncopies=m)*real(cwk(2_i4b:mp1,:sld2),stnd) -   &
                  spread(real(cwk(1_i4b,:sld2),stnd),dim=1,ncopies=m)*aimag(cwk(2_i4b:mp1,:sld2))
!
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (updatepscd2)
!
        psvec(:nf)    = psvec(:nf)    + psmatb(1_i4b,:nf)
        psmat(:m,:nf) = psmat(:m,:nf) + psmatb(2_i4b:mp1,:nf)
!
        coher(:m,:nf)    = coher(:m,:nf)  + rcospect(:m,:nf) 
        phase(:m,:sld2)  = phase(:m,:sld2) + icospect(:m,:sld2)
!
!$OMP END CRITICAL (updatepscd2)
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( seg, cwk, psmatb, rcospect, icospect )
!
!$OMP END PARALLEL
!
    else
!
#endif
!
!       ALLOCATE WORK MATRICES.
!
        allocate( seg(mp1,sl), cwk(mp1,nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(:mp1,l+1_i4b:sl) = zero
        end if    
!
        do k=1_i4b, nseg
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg(1_i4b,:l2)     = vec(i1:i2)
            seg(2_i4b:mp1,:l2) = mat(:m,i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENTS IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rm( seg(:mp1,:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:mp1,:l2) = seg(:mp1,:l2)*spread( wk(:l2), dim=1, ncopies=mp1 )
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(:mp1,l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT FOR THE SERIES.
!
            call real_fft( seg(:mp1,:sl), cwk(:mp1,:nf), true )
!
!           UPDATE POWER SPECTRUM ESIMATES.
!
            psvec(:sld2) = psvec(:sld2)                    +   &
                           real(cwk(1_i4b,:sld2),stnd)**2  +   &
                           aimag(cwk(1_i4b,:sld2))**2
            psvec(nf)   = psvec(nf) + real(cwk(1_i4b,nf),stnd)**2
!
            psmat(:m,:sld2) = psmat(:m,:sld2)                     +   &
                              real(cwk(2_i4b:mp1,:sld2),stnd)**2  +   &
                              aimag(cwk(2_i4b:mp1,:sld2))**2
            psmat(:m,nf)   = psmat(:m,nf) + real(cwk(2_i4b:mp1,nf),stnd)**2
!
!           UPDATE CO-SPECTRUM ESIMATES.
!
            coher(:m,:sld2)  = coher(:m,:sld2)                                                            +   &
                      spread(real(cwk(1_i4b,:sld2),stnd),dim=1,ncopies=m)*real(cwk(2_i4b:mp1,:sld2),stnd) +   &
                      spread(aimag(cwk(1_i4b,:sld2)),dim=1,ncopies=m)*aimag(cwk(2_i4b:mp1,:sld2))
            coher(:m,nf)    = coher(:m,nf) + real(cwk(1_i4b,nf),stnd)*real(cwk(2_i4b:mp1,nf),stnd)
!
!           UPDATE QUADRATURE-SPECTRUM ESIMATES.
!
            phase(:m,:sld2)  = phase(:m,:sld2)                                                           +   &
                         spread(aimag(cwk(1_i4b,:sld2)),dim=1,ncopies=m)*real(cwk(2_i4b:mp1,:sld2),stnd) -   &
                         spread(real(cwk(1_i4b,:sld2),stnd),dim=1,ncopies=m)*aimag(cwk(2_i4b:mp1,:sld2))
!
        end do
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( seg, cwk )
!
#ifdef _OPENMP
!
    end if
!
#endif
!
!   NORMALIZE THE POWER AND CROSS SPECTRA ESTIMATES.
!
    c1 = one/( sumw*real(l*nseg, stnd) )
!
    if ( normpsd2 ) then
!
        c2 = (real( l, stnd)/pi)*c1
        i1 = 1_i4b
        i2 = nf
!
    else
        c2 = two*c1
        i1 = 2_i4b
        i2 = sld2
!
        psvec(1_i4b)    = c1*psvec(1_i4b)
        psmat(:m,1_i4b) = c1*psmat(:m,1_i4b)
!
        coher(:m,1_i4b)  = c1*coher(:m,1_i4b)
        phase(:m,1_i4b)  = c1*phase(:m,1_i4b)
!
        psvec(nf)    = c1*psvec(nf)
        psmat(:m,nf) = c1*psmat(:m,nf)
!
        coher(:m,nf) = c1*coher(:m,nf)
        phase(:m,nf)  = c1*phase(:m,nf)
!
    end if
!
    psvec(i1:i2)    = c2*psvec(i1:i2)
    psmat(:m,i1:i2) = c2*psmat(:m,i1:i2)
!
    coher(:m,i1:i2)   = c2*coher(:m,i1:i2)
    phase(:m,i1:i2)   = c2*phase(:m,i1:i2)
!
!   SMOOTH THE POWER AND CROSS SPECTRA ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
!       SMOOTH POWER SPECTRUM OF THE FIRST SERIES.
!
        call dan_filter_rv( psvec(:nf), nsmooth, sym=one )
!
!       SMOOTH POWER SPECTRUM OF THE OTHER SERIES.
!
        call dan_filter_rm( psmat(:m,:nf), nsmooth, sym=one )
!
!       SMOOTH CO-SPECTRUM.
!
        call dan_filter_rm( coher(:m,:nf), nsmooth, sym=one )
!
!       SMOOTH QUADRATURE-SPECTRUM.
!
        call dan_filter_rm( phase(:m,:nf), nsmooth, sym=-one )
!
    end if
!
    if ( out_co_spect ) then
!
!       OUTPUT THE CO-SPECTRUM IF REQUIRED.
!
        co_spect(:m,:nf) = coher(:m,:nf)
!
    end if
!
    if ( out_quad_spect ) then
!
!       OUTPUT THE QUADRATURE-SPECTRUM IF REQUIRED.
!
        quad_spect(:m,:nf) = -phase(:m,:nf)
!
    end if
!
!   ALLOCATE WORK MATRICES.
!
    allocate( temp(m), magni(m,nf), stat = iok )
!
    if ( iok/=0 ) call merror( name_proc//allocate_error )
!
    do k=1_i4b, nf
!
!       COMPUTE MAGNITUDE OF THE CROSS-SPECTRUM.
!
        magni(:m,k) = coher(:m,k)**2 + phase(:m,k)**2
!
!       COMPUTE PHASE OF THE CROSS-SPECTRUM.
!
        where ( coher(:m,k)/=zero .or. phase(:m,k)/=zero )
            temp(:m)   = atan2( phase(:m,k), coher(:m,k) )
            phase(:m,k) = (one/twopi)*temp(:m) + merge( zero, one, temp(:m)>=zero ) 
        end where
!
!       COMPUTE THE SQUARED COHERENCY.
!
        if ( psvec(k)/=zero ) then
!
            where ( psmat(:m,k)/=zero )
                coher(:m,k) = min( magni(:m,k)/(psvec(k)*psmat(:m,k)), one )
            elsewhere
                coher(:m,k) = zero
            end where
!
        else
!
            coher(:m,k) = zero     
!
        end if
!
    end do   
!
    if ( out_ampli ) then
!
!       COMPUTE AND OUTPUT AMPLITUDE OF THE CROSS-SPECTRUM IF REQUIRED.
!
        ampli(:m,:nf) = sqrt( magni(:m,:nf) )
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( sl, stnd )
        freq(:nf) = arth(zero, c1, nf)
    end if
!
!   OUTPUT DEGREES OF FREEDOM, BANDWIDTH, CONFIDENCE LIMIT FACTORS 
!   OR PROBABILITIES FOR SQUARED COHERENCIES IF REQUIRED.
!
    if ( out_dof ) then
!
        allocate( edof2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
        i = sl - l
!
        edof2(:nf) = estim_dof2( wk(:l), i, win=win, nsmooth=nsmooth, nseg=nseg, overlap=overlap )
!
        if ( present(edof) ) then
            edof(:nf) = edof2(:nf)
        end if
!
        if ( present(bandwidth) ) then
            bandwidth(:nf) = edof2(:nf)/( two*real( n, stnd ) )
        end if
!
        if ( present(conlwr) .or. present(conupr) .or. present(testcoher) ) then
            call comp_conflim_rv( edof2(:nf), probtest=probtest, conlwr=conlwr,       &
                                  conupr=conupr, testcoher=testcoher )
        end if
!
!       COMPUTE AND OUTPUT SIGNIFICANCE PROBABILITY FOR SQUARED COHERENCIES.
!
        if ( out_prob_coher ) then
!
            sumw = nan()
!
            do i = 1_i4b, nf
!            
                df = edof2(i)   
!            
                if ( df>two ) then
!            
                    where ( coher(:m,i)/=one )
                        temp(:m) = ( (df/two - one)*coher(:m,i) )/( one - coher(:m,i) )
                    elsewhere
                        temp(:m) = zero
                    end where
!            
                    df2 = df - two
!            
                    prob_coher(:m,i) = probf2( temp(:m), two, df2, true )
!            
                    where( coher(:m,i)==one ) prob_coher(:m,i) = zero
!
                else
!
                    prob_coher(:m,i) = sumw
!
                end if
!
            end do
!
        end if
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( edof2 )
!
    end if
!
!   DEALLOCATE WORK ARRAYS.
!
    deallocate( temp, magni )
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE cross_spctrm2_rm
! __________________________________
!
    end subroutine cross_spctrm2_rm
!
! =========================================================================================
!
    subroutine power_spectrum_rv( vec, psvec, freq, fftvec, edof, bandwidth, conlwr, conupr,   &
                                  initfft, normpsd, smooth_param, trend, win, taperp, probtest )
!
! Purpose
! _______
!
!   Subroutine POWER_SPECTRUM computes a Fast Fourier Transform (FFT) estimate
!   of the power spectrum of a real time series, VEC. The real valued sequence
!   VEC must be of even length.
!
!   The Power Spectral Density (PSD) estimates are returned in units which are
!   the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the real time series for which the power spectrum
!                 must be estimated.
!                 If WIN/=2 or TREND=1, 2 or 3,  VEC is used
!                 as workspace and is transformed.
!
!                 Size(VEC) must be an even (positive) integer greater or equal to 4.
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = size(VEC)/2 + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/size(VEC)
!                 for i=1,2, ... , (size(VEC)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = size(VEC)/2 + 1 .
!
!   FFTVEC        (OUTPUT, OPTIONAL) complex(stnd), dimension(:)
!                 On exit, a complex vector of length (size(VEC)/2)+1 containing
!                 the Fast Fourier Transform of the product of the (detrended, e.g. 
!                 the TREND argument) real time series VEC with the choosen window 
!                 function (e.g. The WIN argument).
!
!                 FFTVEC must verify:  size(FFTVEC) = size(VEC)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the bandwidth of the power spectrum estimates.
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the
!                 PSVEC(:) argument) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine POWER_SPECTRUM in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine POWER_SPECTRUM. This call to INITFFT must have the following form: 
!
!                       call init_fft( size(VEC)/2 )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   POWER_SPECTRUM and a call to END_FFT is also done before leaving
!                   subroutine POWER_SPECTRUM.
!
!                 The default is INITFFT=true .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the PSD estimates are normalized in such
!                   a way that the total area under the power spectrum is equal to the variance of
!                   the time series VEC.
!                 - NORMPSD is set to false, the sum of the PSD estimates
!                   (e.g. sum( PSVEC(2:) ) is equal to the variance of the time series.
!
!                 The default is NORMPSD=true .
!
!   SMOOTH_PARAM  (INPUT, OPTIONAL) integer(i4b), dimension(:)
!                 if SMOOTH_PARAM is used, the PSD estimates are computed by repeated
!                 smoothing of the periodogram with modified Daniell weights.
!
!                 On entry, SMOOTH_PARAM(:) gives the array of the half-lengths of the
!                 modified Daniell filters to be applied.
!
!                 All the values in SMOOTH_PARAM(:) must be greater than 0 and less
!                 than size(VEC)/2+1 .
!
!                 Size(SMOOTH_PARAM) must be greater or equal to 1.
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The mean of the time series is removed before computing the spectrum
!                 - TREND=+2 The drift from the time series is removed before computing the spectrum
!                   by using the formula: drift = (VEC(size(VEC)) - VEC(1))/(size(VEC) - 1)
!                 - TREND=+3 The least-squares line from the time series is removed before
!                   computing the spectrum.
!
!                 For other values of TREND nothing is done before estimating the power spectrum.
!
!                 The default is TREND=1, e.g. the mean of the time series is removed before the
!                 computations.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power spectrum. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the selected
!   data window (e.g. WIN=1,2,3,4,5,6) is applied to the time series and the PSD estimates are
!   computed by the FFT of this transformed time series. Optionally, theses PSD estimates
!   may then be smoothed in the frequency domain by modified Daniell filters (e.g. if
!   SMOOTH_PARAM is used).
!
!   The computed equivalent number of degrees of freedom and bandwidth must be divided by two
!   for the zero and Nyquist frequencies.
!
!   Furthermore, the computed equivalent number of degrees of freedom, bandwidth, lower and
!   upper (1-PROBTEST) * 100% confidence limit factors are not right near the zero and Nyquist frequencies
!   if the PSD estimates have been smoothed by modified Daniell filters. The reason is that
!   POWER_SPECTRUM assumes that smoothing involves averaging independent frequency
!   ordinates. This is true except near the zero and Nyquist frequencies where an average
!   may contain contributions from negative frequencies, which are identical to and hence not
!   independent of positive frequency spectral values. Thus, the number of degrees of freedom
!   in PSD estimates near the 0 and Nyquist frequencies are as little as half the number
!   of degrees of freedom of the spectral estimates away from these frequency extremes if
!   the optional argument SMOOTH_PARAM is used.
!      
!   If the optional argument SMOOTH_PARAM is used, the computed equivalent number of degrees
!   of freedom, bandwidth, lower and upper (1-PROBTEST) * 100% confidence limit factors are right
!   for PSD estimates  at frequencies
!      
!        (i-1)/Size(VEC)  for i= (nparam+1)/2 + 1   to  ( Size(VEC) - nparam + 1)/2
!      
!   where nparam = 2 * (2+sum(SMOOTH_PARAM(:)))- 1,
!   (e.g. for frequencies i/Size(VEC) for i = (nparam+1)/2, ... , ( Size(VEC) - nparam - 1)/2   ).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, arth
    use Reals_Constants,   only : zero, half, one, two, pi, c5_m2
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : tseries_error10, tseries_error14, tseries_error15,    &
                                  tseries_error50, tseries_error54, tseries_error59
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : pinvq2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:) :: vec
    real(stnd), intent(out),    dimension(:) :: psvec
!
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out),               optional :: edof, bandwidth, conlwr, conupr
    real(stnd), intent(out), dimension(:), optional :: freq
!
    complex(stnd), intent(out), dimension(:), optional :: fftvec
!
    integer(i4b),  intent(in),               optional :: trend, win
    integer(i4b),  intent(in), dimension(:), optional :: smooth_param
!
    logical(lgl),  intent(in), optional :: normpsd, initfft
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, nd2, nf,  win2, trend2, nparam
!
    real(stnd)                       :: c1, c2, sumw, edof2, probtest2
    real(stnd), dimension(size(vec)) :: wk
!
    complex(stnd), dimension(size(vec)/2+1)  :: cwk
!
    logical(lgl)  :: normpsd2, initfft2, smooth, out_freq, out_fftvec, out_dof
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='power_spectrum'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n   = size( vec )
!
    if ( n<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    nd2 = n/2_i4b
    if ( n/=2*nd2 )   &
    call merror( name_proc//tseries_error14 )
!
    nf  = nd2 + 1_i4b
!
    call assert( logical(int(size(psvec),i4b)==nf,lgl),     &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_fftvec = present( fftvec )
!
    if ( out_fftvec ) then
        call assert( logical(int(size(fftvec),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth )        &
              .or. present( conlwr ) .or. present( conupr ) 
!
    probtest2 = c5_m2
!
    if ( present(probtest) ) then
        if ( zero>=probtest .or. probtest>=one ) then
            call merror( name_proc//tseries_error59 )
        else
            probtest2 = probtest
        end if
    end if
!
    if ( present(trend) ) then
        trend2 = trend
    else
        trend2 = 1_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
!
    else
!
        win2 = 3_i4b
!
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( nd2 )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    smooth = false
!
    if ( present(smooth_param) ) then
!
        nparam = size( smooth_param )
!
        if ( nparam<=0_i4b  )    &
        call merror( name_proc//tseries_error50 )
!
!       CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIEL FILTERS.
!
        if ( any( smooth_param(:nparam)<=0_i4b .or. smooth_param(:nparam)>=nf  ) )   &
        call merror( name_proc//tseries_error54 )
!
        smooth = true
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trend2>=1_i4b .and. trend2<=3_i4b ) then
        call detrend_rv( vec(:n), trend2 )
    end if
!
!   TAPER THE TIME SERIES IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:n)  = data_window( n, win2, taperp=taperp )
!
!       APPLY DATA WINDOW.
!
        vec(:n) = vec(:n)*wk(:n)
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:n), wk(:n) )
!
    else    
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:n)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( n, stnd )
!
    end if
!
!   COMPUTE PSD ESTIMATES
!
    call real_fft( vec(:n), cwk(:nf), true)
!
    psvec(:nd2) = real(cwk(:nd2),stnd)**2 + aimag(cwk(:nd2))**2
    psvec(nf)   = real(cwk(nf),stnd)**2
!
!   NORMALIZE THE PSD ESTIMATES.
!
    c1 = one/(sumw*real( n, stnd))
!
    if ( normpsd2 ) then
!
!       CHATFIELD DEFINITION OF THE SPECTRAL DENSITY FUNCTION.
!
        c2 = real( n, stnd)/pi
        psvec(1_i4b:nf) = (c1*c2)*psvec(1_i4b:nf)
!
    else
!
!       POLLOCK DEFINITION OF THE PERIODOGRAM. IF THE ESTIMATES
!       ARE DIVIDED BY 4*PI, THIS IS THE POLLOCK DEFINITION OF
!       THE SPECTRAL DENSITY FUNCTION.
!
        psvec(1_i4b)     = c1*psvec(1_i4b)
        psvec(2_i4b:nd2) = (two*c1)*psvec(2_i4b:nd2)
        psvec(nf)        = c1*psvec(nf)
!
    end if
!
!   SMOOTH THE PSD ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
        call moddan_filter_rv( psvec(:nf), smooth_param(:nparam), sym=one )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( n, stnd )
        freq(:nf) = arth(zero, c1, nf)
    end if
!
!   OUTPUT FFT OF THE DATA IF REQUIRED.
!
    if ( out_fftvec ) then
        fftvec(:nf) = cwk(:nf)
    end if
!
!   OUTPUT DEGREES OF FREEDOM, BANDWIDTH OR CONFIDENCE LIMIT FACTORS IF REQUIRED.
!
    if ( out_dof ) then
!
        edof2 = estim_dof( wk(:n), smooth_param=smooth_param(:nparam) )
!
        if ( present(edof) ) then
            edof = edof2
        end if
!
        if ( present(bandwidth) ) then
            bandwidth = edof2/( two*real( n, stnd ) )
        end if
!
        c1 = probtest2*half
        c2 = one - c1
!
        if ( present(conlwr) ) then
            if ( edof2>=half ) then
                conlwr = edof2/pinvq2( c2, edof2 )
            else
                conlwr = nan()
            end if
        end if
!
        if ( present(conupr) ) then
            if ( edof2>=half ) then
                conupr = edof2/pinvq2( c1, edof2 )
            else
                conupr = nan()
            end if
        end if
!
    end if
!
!   DEALLOCATE WORK ARRAYS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE power_spectrum_rv
! ___________________________________
!
    end subroutine power_spectrum_rv
!
! =========================================================================================
!
    subroutine power_spectrum_rm( mat, psmat, freq, fftmat, edof, bandwidth, conlwr, conupr,   &
                                  initfft, normpsd, smooth_param, trend, win, taperp, probtest )
!
! Purpose
! _______
!
!   Subroutine POWER_SPECTRUM computes a Fast Fourier Transform (FFT) estimate
!   of the power spectra of the rows of the real matrix, MAT. size(MAT,2) must
!   be of even length.
!
!   The Power Spectral Density (PSD) estimates are returned in units which are
!   the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   MAT           (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                 On entry, the real time series for which power spectra
!                 must be estimated. Each row of MAT is a real time series.
!                 If WIN/=2 or TREND=1, 2 or 3,  MAT is used
!                 as workspace and is transformed.
!
!                 Size(MAT,2) must be an even (positive) integer greater or equal to 4.
!
!   PSMAT         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix containing the Power Spectral Density
!                 (PSD) estimates for each row of the real matrix MAT.  
!
!                 The shape of PSMAT must verify:
!
!                 - size(PSMAT,1) = size(MAT,1) ;
!                 - size(PSMAT,2) = size(MAT,2)/2 + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(MAT,2)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/size(VEC)
!                 for i=1,2, ... , (size(MAT,2)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = size(MAT,2)/2 + 1 .
!
!   FFTMAT        (OUTPUT, OPTIONAL) complex(stnd), dimension(:,:)
!                 On exit, a complex matrix containing the Fast Fourier Transform 
!                 of the product of the (detrended, e.g. the TREND argument) 
!                 real time series in each row of MAT with the choosen window 
!                 function (e.g. The WIN argument).
!
!                 The shape of FFTMAT must verify:
!
!                 - size(FFTMAT,1) = size(MAT,1) ;
!                 - size(FFTMAT,2) = size(MAT,2)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the bandwidth of the power spectrum estimates.
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100%
!                 confidence limit factors, respectively. Multiply the PSD estimates (e.g. the
!                 PSMAT(:,:) argument) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine POWER_SPECTRUM in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine POWER_SPECTRUM. This call to INITFFT must have the following form: 
!
!                       call init_fft( (/ size(MAT,1), size(MAT,2)/2 /), dim=2_i4b )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   POWER_SPECTRUM and a call to END_FFT is also done before leaving
!                   subroutine POWER_SPECTRUM
!
!                 The default is INITFFT=true .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the PSD estimates are normalized in such
!                   a way that the total area under the power spectrum is equal to the variance of
!                   the time series MAT.
!                 - NORMPSD = false, the sum of the PSD estimates
!                   for each row of MAT (e.g. sum( PSMAT(:,2:), dim=2 ) is equal to the variance of
!                   the corresponding time series.
!
!                 The default is NORMPSD=true .
!
!   SMOOTH_PARAM  (INPUT, OPTIONAL) integer(i4b), dimension(:)
!                 if SMOOTH_PARAM is used, the PSD estimates are computed by repeated
!                 smoothing of the periodogram with modified Daniell weights.
!
!                 On entry, SMOOTH_PARAM(:) gives the array of the half-lengths of the
!                 modified Daniell filters to be applied.
!
!                 All the values in SMOOTH_PARAM(:) must be greater than 0 and less
!                 than size(MAT,2)/2 + 1 .
!
!                 Size(SMOOTH_PARAM) must be greater or equal to 1.
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If
!
!                 - TREND=+1 The means of the time series are removed before computing the spectra
!                 - TREND=+2 The drifts from the time series are removed before computing the spectra
!                   by using the formula: drift(i) = (MAT(i,size(MAT,2)) - MAT(i,1))/(size(MAT,2) - 1)
!                 - TREND=+3 The least-squares lines from the time series are removed before
!                   computing the spectra.
!
!                 For other values of TREND nothing is done before estimating the power spectra.
!
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power spectrum. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the selected
!   data window (e.g. WIN=1,2,3,4,5,6) is applied to the time series and the PSD estimates are
!   computed by the FFT of these transformed time series. Optionally, theses PSD estimates
!   may then be smoothed in the frequency domain by modified Daniell filters (e.g. if 
!   SMOOTH_PARAM is used).
!
!   The computed equivalent number of degrees of freedom and bandwidth must be divided by two
!   for the zero and Nyquist frequencies.
!
!   Furthermore, the computed equivalent number of degrees of freedom, bandwidth, lower and
!   upper (1-PROBTEST) * 100% confidence limit factors are not right near the zero and Nyquist frequencies
!   if the PSD estimates have been smoothed by modified Daniell filters. The reason is that
!   POWER_SPECTRUM assumes that smoothing involves averaging independent frequency
!   ordinates. This is true except near the zero and Nyquist frequencies where an average
!   may contain contributions from negative frequencies, which are identical to and hence not
!   independent of positive frequency spectral values. Thus, the number of degrees of freedom
!   in PSD estimates near the 0 and Nyquist frequencies are as little as half the number
!   of degrees of freedom of the spectral estimates away from these frequency extremes if
!   the optional argument SMOOTH_PARAM is used.
!      
!   If the optional argument SMOOTH_PARAM is used, the computed equivalent number of degrees
!   of freedom, bandwidth, lower and upper (1-PROBTEST) * 100% confidence limit factors are right
!   for PSD estimates  at frequencies
!      
!        (i-1)/Size(MAT,2)  for i= (nparam+1)/2 + 1   to  ( Size(MAT,2) - nparam + 1)/2
!      
!   where nparam = 2 * (2+sum(SMOOTH_PARAM(:)))- 1,
!   (e.g. for frequencies i/Size(MAT,2) for i = (nparam+1)/2, ... , ( Size(MAT,2) - nparam - 1)/2   ).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, arth
    use Reals_Constants,   only : zero, half, one, two, pi, c5_m2
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : tseries_error10, tseries_error15, tseries_error16,    &
                                  tseries_error50, tseries_error55, tseries_error59
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : pinvq2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:,:) :: mat
!
    real(stnd), intent(out), dimension(:,:)         :: psmat
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out),               optional :: edof, bandwidth, conlwr, conupr
    real(stnd), intent(out), dimension(:), optional :: freq
!
    complex(stnd), intent(out), dimension(:,:), optional :: fftmat
!
    integer(i4b),  intent(in), dimension(:), optional :: smooth_param
    integer(i4b),  intent(in),               optional :: trend, win
!
    logical(lgl),  intent(in), optional :: normpsd, initfft
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, n, nd2, nf, k,  win2, trend2, nparam
!
    real(stnd)                               :: c1, c2, sumw, edof2, probtest2
    real(stnd), dimension(size(mat,2))       :: wk
!
    complex(stnd), dimension(size(mat,1),size(mat,2)/2+1)  :: cwk
!
    logical(lgl)  :: normpsd2, initfft2, smooth, out_freq, out_fftmat, out_dof
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='power_spectrum'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m   = size( mat, 1 )
    n   = size( mat, 2 )
!
    if ( m<=0_i4b ) return
!
    if ( n<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    nd2 = n/2_i4b
!
    if ( n/=2*nd2 )   &
    call merror( name_proc//tseries_error16 )
!
    nf  = nd2 + 1_i4b
!
    call assert( logical(int(size(psmat,1),i4b)==m,lgl),     &
                 logical(int(size(psmat,2),i4b)==nf,lgl),    &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_fftmat = present( fftmat )
!
    if ( out_fftmat ) then
        call assert( logical(int(size(fftmat,1),i4b)==m,lgl),     &
                     logical(int(size(fftmat,2),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth )        &
              .or. present( conlwr ) .or. present( conupr ) 
!
    probtest2 = c5_m2
!
    if ( present(probtest) ) then
        if ( zero>=probtest .or. probtest>=one ) then
            call merror( name_proc//tseries_error59 )
        else
            probtest2 = probtest
        end if
    end if
!
    if ( present(trend) ) then
        trend2 = trend
    else
        trend2 = 1_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
!
    else
!
        win2 = 3_i4b
!
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( (/ m, nd2/), dim=2_i4b )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    smooth = false
!
    if ( present(smooth_param) ) then
!
        nparam = size( smooth_param )
!
        if ( nparam<=0_i4b  )    &
        call merror( name_proc//tseries_error50 )
!
!       CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIEL FILTERS.
!
        if ( any( smooth_param(:nparam)<=0_i4b .or. smooth_param(:nparam)>=nf  ) )   &
        call merror( name_proc//tseries_error55 )
!
        smooth = true
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trend2>=1_i4b .and. trend2<=3_i4b ) then
        call detrend_rm( mat(:m,:n), trend2 )
    end if
!
!   TAPER THE TIME SERIES IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:n)  = data_window( n, win2, taperp=taperp )
!
!       APPLY DATA WINDOW.
!
        do k = 1_i4b, n
            mat(:m,k) = wk(k)*mat(:m,k)
        end do
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:n), wk(:n) )
!
    else    
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:n)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( n, stnd )
!
    end if
!
!   COMPUTE THE PSD ESTIMATES.
!
    call real_fft( mat(:m,:n), cwk(:m,:nf), true)
!
    psmat(:m,:nd2) = real(cwk(:m,:nd2),stnd)**2 + aimag(cwk(:m,:nd2))**2
    psmat(:m,nf)   = real(cwk(:m,nf),stnd)**2
!
!   NORMALIZE THE PSD ESTIMATES.
!
    c1 = one/(sumw*real( n, stnd))
!
    if ( normpsd2 ) then
!
!       CHATFIELD DEFINITION OF THE SPECTRAL DENSITY FUNCTION.
!
        c2 = real( n, stnd)/pi
        psmat(:m,:nf) = (c1*c2)*psmat(:m,:nf)
!
    else
!
!       POLLOCK DEFINITION OF THE PERIODOGRAM. IF THE ESTIMATES
!       ARE DIVIDED BY 4*PI, THIS IS THE POLLOCK DEFINITION OF
!       THE SPECTRAL DENSITY FUNCTION.
!
        psmat(:m,1_i4b)     = c1*psmat(:m,1_i4b)
        psmat(:m,2_i4b:nd2) = (two*c1)*psmat(:m,2_i4b:nd2)
        psmat(:m,nf)        = c1*psmat(:m,nf)
!
    end if
!
!   SMOOTH THE PSD ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
        call moddan_filter_rm( psmat(:m,:nf), smooth_param(:nparam), sym=one )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( n, stnd )
        freq(:nf) = arth(zero, c1, nf)
    end if
!
!   OUTPUT FFT OF THE DATA IF REQUIRED.
!
    if ( out_fftmat ) then
        fftmat(:m,:nf) = cwk(:m,:nf)
    end if
!
!   OUTPUT DEGREES OF FREEDOM, CONFIDENCE LIMIT FACTORS
!   OR BANDWIDTH IF REQUIRED.
!
    if ( out_dof ) then
!
        edof2 = estim_dof( wk(:n), smooth_param=smooth_param(:nparam) )
!
        if ( present(edof) ) then
            edof = edof2
        end if
!
        if ( present(bandwidth) ) then
            bandwidth = edof2/( two*real( n, stnd ) )
        end if
!
        c1 = probtest2*half
        c2 = one - c1
!
        if ( present(conlwr) ) then
            if ( edof2>=half ) then
                conlwr = edof2/pinvq2( c2, edof2 )
            else
                conlwr = nan()
            end if
        end if
!
        if ( present(conupr) ) then
            if ( edof2>=half ) then
                conupr = edof2/pinvq2( c1, edof2 )
            else
                conupr = nan()
            end if
        end if
!
    end if
!
!   DEALLOCATE WORK ARRAYS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE power_spectrum_rm
! ___________________________________
!
    end subroutine power_spectrum_rm
!
! =========================================================================================
!
    subroutine cross_spectrum_rv( vec, vec2, psvec, psvec2, phase, coher, freq, edof, bandwidth, &
                                  conlwr, conupr, testcoher,  ampli, co_spect, quad_spect,       &
                                  prob_coher, initfft, normpsd, smooth_param, trend, win,        &
                                  taperp, probtest )
!
! Purpose
! _______
!
!   Subroutine CROSS_SPECTRUM computes Fast Fourier Transform (FFT) estimates
!   of the power and cross spectra of two real time series, VEC and VEC2. The real valued
!   sequences VEC and VEC2  must be of even length.
!
!   The Power Spectral Density (PSD) and Cross Spectral Density (CSD) estimates are returned
!   in units which are the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the first real time series for which the power and cross spectra
!                 must be estimated.
!                 If WIN/=2 or TREND=1, 2 or 3,  VEC is used as workspace and is transformed.
!
!                 Size(VEC) must be an even (positive) integer greater or equal to 4.
!
!   VEC2          (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the second real time series for which the power and cross spectra
!                 must be estimated.
!                 If WIN/=2 or TREND=1, 2 or 3,  VEC2 is used as workspace and is transformed.
!
!                 VEC2 must verify:  size(VEC2) = size(VEC).
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = size(VEC)/2 + 1 .
!
!   PSVEC2        (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC2)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC2.  
!
!                 PSVEC2 must verify:  size(PSVEC2) = size(VEC)/2 + 1 .
!
!   PHASE         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the phase of the cross spectrum, given in fractions
!                 of a circle (e.g. on the closed interval (0,1) ).
!
!                 PHASE must verify:  size(PHASE) = size(VEC)/2 + 1 .
!
!   COHER         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the squared coherency  estimates for all frequencies.
!
!                 COHER must verify:  size(COHER) = size(VEC)/2 + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/size(VEC)
!                 for i=1,2, ... , (size(VEC)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = size(VEC)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 and cross spectrum estimates.
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the bandwidth of the power and cross spectrum estimates.
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the PSVEC(:) and 
!                 PSVEC2(:) arguments ) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!   
!   TESTCOHER     (OUTPUT, OPTIONAL) real(stnd)
!                 On output, this argument specifies the critical value for testing the null
!                 hypothesis that the squared coherency is zero at the PROBTEST * 100% significance
!                 level (e.g. elements of COHER(:) less than TESTCOHER should be regarded as not
!                 significantly different from zero at the PROBTEST * 100% significance level).
!
!   AMPLI         (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the cross-amplitude spectrum.  
!
!                 AMPLI must verify:  size(AMPLI) = (size(VEC)/2) + 1 .
!
!   CO_SPECT      (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1  containing
!                 the co-spectrum (e.g. the real part of cross-spectrum).  
!
!                 CO_SPECT must verify:  size(CO_SPECT) = (size(VEC)/2) + 1 .
!
!   QUAD_SPECT    (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the quadrature spectrum (e.g. the imaginary part of
!                 cross-spectrum with a minus sign).  
!
!                 QUAD_SPECT must verify:  size(QUAD_SPECT) = (size(VEC)/2) + 1 .
!
!   PROB_COHER    (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the probabilities that the computed sample squared coherencies
!                 came from an ergodic stationary bivariate process with (corresponding)
!                 squared coherencies equal to zero.
!
!                 PROB_COHER  must verify:  size(PROB_COHER) = (size(VEC)/2) + 1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine CROSS_SPECTRUM in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine CROSS_SPECTRUM. This call to INITFFT must have the following form: 
!
!                      call init_fft( size(VEC)/2 )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   CROSS_SPECTRUM and a call to END_FFT is also done before leaving
!                   subroutine CROSS_SPECTRUM.
!
!                 The default is INITFFT=true .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the power and cross spectra estimates are normalized
!                   in such a way that the total area under the power spectrum is equal to the variance
!                   of the time series VEC and VEC2.
!                 - NORMPSD = false, the sum of the PSD estimates
!                   (e.g. sum(PSVEC(2:)) and sum(PSVEC2(2:)) ) is equal to the variance of the corresponding
!                   time series.
!
!                 The default is NORMPSD=true .
!
!   SMOOTH_PARAM  (INPUT, OPTIONAL) integer(i4b), dimension(:)
!                 if SMOOTH_PARAM is used, the power and cross spectra estimates are computed by repeated
!                 smoothing of the periodograms and cross-periodogram with modified Daniell weights.
!
!                 On entry, SMOOTH_PARAM(:) gives the array of the half-lengths of the
!                 modified Daniell filters to be applied.
!
!                 All the values in SMOOTH_PARAM(:) must be greater than 0 and less
!                 than size(VEC)/2+1 .
!
!                 Size(SMOOTH_PARAM) must be greater or equal to 1.
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The mean of the two time series is removed before computing 
!                   the power and cross spectra.
!                 - TREND=+2 The drift from the two time series is removed before computing
!                   the power and cross spectra.
!                 - TREND=+3 The least-squares line from the two time series is removed before
!                   computing the power and cross spectra.
!
!                 For other values of TREND nothing is done before estimating the power and cross spectra.
!
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power and cross spectra. If
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ) and the critical value for
!                 testing the null hypothesis that the squared coherency is zero (e.g.
!                 the TESTCOHER optional argument).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the selected
!   data window (e.g. WIN=1,2,3,4,5,6) is applied to the time series and the PSD and CSD estimates are
!   computed by the FFT of these transformed time series. Optionally, theses PSD and CSD estimates
!   may then be smoothed in the frequency domain by modified Daniell filters (e.g. if argument 
!   SMOOTH_PARAM is used).
!
!   The computed equivalent number of degrees of freedom and bandwidth must be divided by two
!   for the zero and Nyquist frequencies.
!
!   Furthermore, the computed equivalent number of degrees of freedom, bandwidth, lower and
!   upper (1-PROBTEST) * 100% confidence limit factors and critical value for the squared coherency
!   (e.g. arguments EDOF, BANDWIDTH, CONLWR, CONUPR and TESTCOHER) are not right near the zero
!   and Nyquist frequencies if the PSD estimates have been smoothed by modified Daniell filters.
!   The reason is that CROSS_SPECTRUM assumes that smoothing involves averaging independent frequency
!   ordinates. This is true except near the zero and Nyquist frequencies where an average
!   may contain contributions from negative frequencies, which are identical to and hence not
!   independent of positive frequency spectral values. Thus, the number of degrees of freedom
!   in PSD estimates near the 0 and Nyquist frequencies are as little as half the number
!   of degrees of freedom of the spectral estimates away from these frequency extremes if
!   the optional argument SMOOTH_PARAM is used.
!      
!   If the optional argument SMOOTH_PARAM is used, the computed equivalent number of degrees
!   of freedom, bandwidth, lower and upper (1-PROBTEST) * 100% confidence limit factors and
!   critical value for the squared coherency are right for PSD estimates  at frequencies
!      
!        (i-1)/size(VEC)  for i= (nparam+1)/2 + 1   to ( size(VEC) - nparam + 1)/2
!      
!   where nparam = 2 * (2+sum(SMOOTH_PARAM(:)))- 1,
!   (e.g. for frequencies i/size(VEC) for i = (nparam+1)/2, ... , (size(VEC)-nparam-1)/2 ) .
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, assert_eq, arth
    use Reals_Constants,   only : zero, half, one, two, four, pi, twopi, c5_m2
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : tseries_error10, tseries_error14, tseries_error15,    &
                                  tseries_error50, tseries_error54, tseries_error59
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : probf2, pinvf2, pinvq2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:) :: vec, vec2
    real(stnd), intent(out),    dimension(:) :: psvec, psvec2, phase, coher
!
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out),               optional :: edof, bandwidth, conlwr, conupr, testcoher
    real(stnd), intent(out), dimension(:), optional :: freq, ampli, co_spect, quad_spect, prob_coher
!
    integer(i4b),  intent(in),               optional :: trend, win
    integer(i4b),  intent(in), dimension(:), optional :: smooth_param
!
    logical(lgl),  intent(in), optional :: normpsd, initfft
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, nd2, nf, k1, k2, win2, trend2, nparam
!
    real(stnd)                           :: c1, c2, sumw, edof2, probtest2, con
    real(stnd), dimension(size(vec))     :: wk
    real(stnd), dimension(size(vec)/2+1) :: rcospect, icospect
!
    complex(stnd), dimension(size(vec)/2+1)  :: cwk1, cwk2
!
    logical(lgl)  :: normpsd2, initfft2, smooth, out_freq, out_ampli,      &
                     out_co_spect, out_quad_spect, out_prob_coher, out_dof
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='cross_spectrum'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n =  assert_eq( int(size(vec),i4b) ,           &
                    int(size(vec2),i4b) ,          &
                    name_proc )
!
    if ( n<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    nd2 = n/2_i4b
!
    if ( n/=2*nd2 )   &
    call merror( name_proc//tseries_error14 )
!
    nf  = nd2 + 1_i4b
!
    call assert( logical(int(size(psvec),i4b)==nf,lgl),   &
                 logical(int(size(psvec2),i4b)==nf,lgl),  &
                 logical(int(size(phase),i4b)==nf,lgl),   &
                 logical(int(size(coher),i4b)==nf,lgl),   &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_ampli = present( ampli )
    if ( out_ampli ) then
        call assert( logical(int(size(ampli),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_co_spect = present( co_spect )
!
    if ( out_co_spect ) then
        call assert( logical(int(size(co_spect),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_quad_spect = present( quad_spect )
!
    if ( out_quad_spect ) then
        call assert( logical(int(size(quad_spect),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_prob_coher = present( prob_coher )
!
    if ( out_prob_coher ) then
        call assert( logical(int(size(prob_coher),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth ) .or. out_prob_coher   &
              .or. present( conlwr ) .or. present( conupr ) .or. present( testcoher ) 
!
    probtest2 = c5_m2
!
    if ( present(probtest) ) then
        if ( zero>=probtest .or. probtest>=one ) then
            call merror( name_proc//tseries_error59 )
        else
            probtest2 = probtest
        end if
    end if
!
    if ( present(trend) ) then
        trend2 = trend
    else
        trend2 = 1_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
!
    else
!
        win2 = 3_i4b
!
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( nd2 )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    smooth = false
!
    if ( present(smooth_param) ) then
!
        nparam = size( smooth_param )
!
        if ( nparam<=0_i4b  )    &
        call merror( name_proc//tseries_error50 )
!
!       CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIEL FILTERS.
!
        if ( any( smooth_param(:nparam)<=0_i4b .or. smooth_param(:nparam)>=nf ) )   &
        call merror( name_proc//tseries_error54 )
!
        smooth = true
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trend2>=1_i4b .and. trend2<=3_i4b ) then
!
        call detrend_rv( vec(:n),  trend2 )
        call detrend_rv( vec2(:n), trend2 )
!
    end if
!
!   TAPER THE TIME SERIES IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:n)  = data_window( n, win2, taperp=taperp )
!
!       APPLY DATA WINDOW.
!
        vec(:n) = vec(:n)*wk(:n)
!
        vec2(:n) = vec2(:n)*wk(:n)
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:n), wk(:n) )
!
    else    
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:n)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( n, stnd )
!
    end if
!
!   COMPUTE DEGREES OF FREEDOM IF REQUIRED.
!
    if ( out_dof ) then
!
        edof2 = estim_dof( wk(:n), smooth_param=smooth_param(:nparam) )
!
    end if
!
!   COMPUTE PSD AND CSD ESTIMATES .
!
!   COMPUTE FFT FOR THE FIRST SERIES.
!
    call real_fft( vec(:n), cwk1(:nf), true)
!
!   COMPUTE FFT FOR THE SECOND SERIES.
!
    call real_fft( vec2(:n), cwk2(:nf), true)
!
!   COMPUTE POWER SPECTRUM ESIMATES.
!
    psvec(:nd2) = real(cwk1(:nd2),stnd)**2 + aimag(cwk1(:nd2))**2
    psvec(nf)   = real(cwk1(nf),stnd)**2
!
    psvec2(:nd2) = real(cwk2(:nd2),stnd)**2 + aimag(cwk2(:nd2))**2
    psvec2(nf)   = real(cwk2(nf),stnd)**2
!
!   COMPUTE CO-SPECTRUM ESIMATES.
!
    rcospect(:nd2)  = real(cwk1(:nd2),stnd)*real(cwk2(:nd2),stnd) +   &
                      aimag(cwk1(:nd2))*aimag(cwk2(:nd2))
    rcospect(nf)    = real(cwk1(nf),stnd)*real(cwk2(nf),stnd)
!
!   COMPUTE QUADRATURE-SPECTRUM ESIMATES.
!
    icospect(:nd2)  = aimag(cwk1(:nd2))*real(cwk2(:nd2),stnd)    -   &
                      real(cwk1(:nd2),stnd)*aimag(cwk2(:nd2))
    icospect(nf)    = zero
!
!   NORMALIZE THE POWER AND CROSS SPECTRA ESTIMATES.
!
    c1 = one/(sumw*real( n, stnd))
!
    if ( normpsd2 ) then
!
!       CHATFIELD DEFINITION OF THE SPECTRAL DENSITY FUNCTION.
!
        c2 = (real( n, stnd)/pi)*c1
        k1 = 1_i4b
        k2 = nf
!
    else
!
!       POLLOCK DEFINITION OF THE PERIODOGRAM. IF THE ESTIMATES
!       ARE DIVIDED BY 4*PI, THIS IS THE POLLOCK DEFINITION OF
!       THE SPECTRAL DENSITY FUNCTION.
!
        c2 = two*c1
        k1 = 2_i4b
        k2 = nd2
!
        psvec(1_i4b)  = c1*psvec(1_i4b)
        psvec2(1_i4b) = c1*psvec2(1_i4b)
!
        rcospect(1_i4b) = c1*rcospect(1_i4b)
        icospect(1_i4b) = c1*icospect(1_i4b)
!
        psvec(nf)  = c1*psvec(nf)
        psvec2(nf) = c1*psvec2(nf)
!
        rcospect(nf)  = c1*rcospect(nf)
        icospect(nf)  = c1*icospect(nf)
!
    end if
!
    psvec(k1:k2)  = c2*psvec(k1:k2)
    psvec2(k1:k2) = c2*psvec2(k1:k2)
!
    rcospect(k1:k2) = c2*rcospect(k1:k2)
    icospect(k1:k2) = c2*icospect(k1:k2)
!
!   SMOOTH THE POWER AND CROSS SPECTRA ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
!       SMOOTH POWER SPECTRUM OF THE FIRST SERIES.
!
        call moddan_filter_rv( psvec(:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH POWER SPECTRUM OF THE SECOND SERIES.
!
        call moddan_filter_rv( psvec2(:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH CO-SPECTRUM.
!
        call moddan_filter_rv( rcospect(:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH QUADRATURE-SPECTRUM.
!
        call moddan_filter_rv( icospect(:nf), smooth_param(:nparam), sym=-one )
!
    end if
!
    if ( out_co_spect ) then
!
!       OUTPUT THE CO-SPECTRUM IF REQUIRED.
!
        co_spect(:nf) = rcospect(:nf)
!
    end if
!
    if ( out_quad_spect ) then
!
!       OUTPUT THE QUADRATURE-SPECTRUM IF REQUIRED.
!
        quad_spect(:nf) = -icospect(:nf)
!
    end if
!
!   COMPUTE PHASE OF THE CROSS-SPECTRUM.
!
    where ( rcospect(:nf)==zero .and. icospect(:nf)==zero )
        phase(:nf) = zero
    elsewhere
        wk(:nf)   = atan2( icospect(:nf), rcospect(:nf) )
        phase(:nf) = (one/twopi)*wk(:nf) + merge( zero, one, wk(:nf)>=zero ) 
    end where
!
!   COMPUTE MAGNITUDE OF THE CROSS-SPECTRUM.
!
    wk(:nf) = rcospect(:nf)**2 + icospect(:nf)**2
!
!   COMPUTE THE SQUARED COHERENCIES FOR ALL FREQUENCIES.
!
    where ( psvec(:nf)/=zero .and. psvec2(:nf)/=zero )
        rcospect(:nf)  = psvec(:nf)*psvec2(:nf)
        coher(:nf)     = min( wk(:nf)/rcospect(:nf), one )
    elsewhere
        coher(:nf) = zero
    end where 
!
    if ( out_ampli ) then
!
!       COMPUTE AND OUTPUT AMPLITUDE OF THE CROSS-SPECTRUM IF REQUIRED.
!
        ampli(:nf) = sqrt( wk(:nf) )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( n, stnd )
        freq(:nf) = arth(zero, c1, nf)
!
    end if
!
!   OUTPUT DEGREES OF FREEDOM, BANDWIDTH, CONFIDENCE LIMIT FACTORS 
!   OR PROBABILITIES FOR SQUARED COHERENCIES IF REQUIRED.
!
    if ( out_dof ) then
!
        if ( present(edof) ) then
            edof = edof2
        end if
!
        if ( present(bandwidth) ) then
            bandwidth = edof2/( two*real( n, stnd ) )
        end if
!
        c1 = probtest2*half
        c2 = one - c1
!
        sumw = nan()
!
        if ( present(conlwr) ) then
            if ( edof2>=half ) then
                conlwr = edof2/pinvq2( c2, edof2 )
            else
                conlwr = sumw
            end if
        end if
!
        if ( present(conupr) ) then
            if ( edof2>=half ) then
                conupr = edof2/pinvq2( c1, edof2 )
            else
                conupr = sumw
            end if
        end if
!
        if ( present(testcoher) ) then
            if ( edof2>two ) then
!
!                c1 = two/(edof2-two)
!                c2 = one - probtest2**c1
!                testcoher = max( min( one, c2 ), zero )
!
                 con = one - probtest2
                 c2  = edof2 - two
                 c1  = two*pinvf2( con, two, c2 )
                 testcoher = c1/( c2 + c1 )
            else
                testcoher = sumw
            end if
        end if
!
!       COMPUTE AND OUTPUT SIGNIFICANCE PROBABILITY FOR SQUARED COHERENCIES.
!
        if ( out_prob_coher ) then
!
            c2 = edof2 - two    
!
            if ( c2>zero ) then
!
                where( coher(2_i4b:nd2)/=one )
                    rcospect(2_i4b:nd2) = ( (edof2/two - one)*coher(2_i4b:nd2) )/  &
                                          ( one - coher(2_i4b:nd2) )
                elsewhere
                    rcospect(2_i4b:nd2) = zero
                end where 
!            
                prob_coher(2_i4b:nd2) = probf2( rcospect(2_i4b:nd2), two, c2, true )
!            
                where( coher(2_i4b:nd2)==one ) prob_coher(2_i4b:nd2) = zero
!
            else
!
                prob_coher(2_i4b:nd2) = sumw
!
            end if
!            
            c2 = edof2/two - two    
!
            if ( c2>zero ) then
!            
                do k1 = 1_i4b, nf, nd2
!            
                    c1 = coher(k1)
!            
                    if ( c1/=one ) then
!            
                        con = ( (edof2/four - one)*c1 )/( one - c1 )
!            
                        prob_coher(k1) = probf2( con, two, c2, true )
!            
                    else
                        prob_coher(k1) = zero
                    end if 
!
                end do
!
            else
!
                prob_coher(1_i4b) = sumw
                prob_coher(nf)    = sumw
!
            end if
!
        end if
!
    end if
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE cross_spectrum_rv
! ___________________________________
!
    end subroutine cross_spectrum_rv
!
! =========================================================================================
!
    subroutine cross_spectrum_rm( vec, mat, psvec, psmat, phase, coher, freq, edof, bandwidth, &
                                  conlwr, conupr, testcoher, ampli, co_spect, quad_spect,      &
                                  prob_coher, initfft, normpsd, smooth_param, trend, win,      &
                                  taperp, probtest )
!
! Purpose
! _______
!
!   Subroutine CROSS_SPECTRUM computes Fast Fourier Transform (FFT) estimates
!   of the power and cross spectra of the real time series, VEC, and the multi-channel
!   real time series MAT.
!
!   The Power Spectral Density (PSD) and Cross Spectral Density (CSD) estimates are returned
!   in units which are the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the real time series for which the power and cross spectra
!                 must be estimated.
!                 If WIN/=2 or TREND=1, 2 or 3,  VEC is used as workspace and is transformed.
!
!                 Size(VEC) must be an even (positive) integer greater or equal to 4.
!
!   MAT           (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                 On entry, the multi-channel real time series for which the power 
!                 and cross spectra must be estimated. Each row of MAT is a real time series.
!                 If WIN/=2 or TREND=1, 2 or 3,  MAR is used as workspace and is transformed.
!
!                 The shape of MAT must verify:  size(MAT,2) = size(VEC).
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = size(VEC)/2 + 1 .
!
!   PSMAT         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the Power Spectral Density (PSD) estimates of each row of MAT.
!
!                 The shape of PSMAT must verify:
!
!                 - size(PSMAT,1) =  size(MAT,1) ;
!                 - size(PSMAT,2) =  size(VEC)/2 + 1 .
!
!   PHASE         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the phase of the cross spectrum, given in fractions
!                 of a circle (e.g. on the closed interval (0,1) ).
!
!                 The shape of PHASE must verify:
!
!                 - size(PHASE,1) =  size(MAT,1) ;
!                 - size(PHASE,2) =  size(VEC)/2 + 1 .
!
!   COHER         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the squared coherency  estimates for all frequencies.
!
!                 The shape of COHER must verify:
!
!                 - size(COHER,1) =  size(MAT,1) ;
!                 - size(COHER,2) =  size(VEC)/2 + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length (size(VEC)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/size(VEC)
!                 for i=1,2, ... , (size(VEC)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = size(VEC)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 and cross spectrum estimates.
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the bandwidth of the power and cross spectrum estimates.
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the PSVEC(:) and 
!                 PSMAT(:,:) arguments ) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!   
!   TESTCOHER     (OUTPUT, OPTIONAL) real(stnd)
!                 On output, this argument specifies the critical value for testing the null
!                 hypothesis that the squared coherency is zero at the PROBTEST * 100% significance
!                 level (e.g. elements of COHER(:,:) less than TESTCOHER should be regarded as not
!                 significantly different from zero at the PROBTEST * 100% significance level).
!
!   AMPLI         (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the cross-amplitude spectra.
!
!                 The shape of AMPLI must verify:
!
!                 - size(AMPLI,1) =  size(MAT,1) ;
!                 - size(AMPLI,2) =  (size(VEC)/2) + 1 .
!
!   CO_SPECT      (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the co-spectra (e.g. the real part of cross-spectra).  
!
!                 The shape of CO_SPECT must verify:
!
!                 - size(CO_SPECT,1) =  size(MAT,1) ;
!                 - size(CO_SPECT,2) =  (size(VEC)/2) + 1 .
!
!   QUAD_SPECT    (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the quadrature spectrum (e.g. the imaginary part of
!                 cross-spectrum with a minus sign).  
!
!                 The shape of QUAD_SPECT must verify:
!
!                 - size(QUAD_SPECT,1) =  size(MAT,1) ;
!                 - size(QUAD_SPECT,2) =  (size(VEC)/2) + 1 .
!
!   PROB_COHER    (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and (size(VEC)/2)+1 columns
!                 containing the probabilities that the computed sample squared coherencies
!                 came from an ergodic stationary bivariate process with (corresponding)
!                 squared coherencies equal to zero.
!
!                 The shape of PROB_COHER must verify:
!
!                 - size(PROB_COHER,1) =  size(MAT,1) ;
!                 - size(PROB_COHER,2) =  (size(VEC)/2) + 1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine CROSS_SPECTRUM in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine CROSS_SPECTRUM. This call to INITFFT must have the following form: 
!
!                        call init_fft( (/ size(MAT,1), size(MAT,2)/2 /), dim=2_i4b )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   CROSS_SPECTRUM and a call to END_FFT is also done before leaving
!                   subroutine CROSS_SPECTRUM.
!
!                 The default is INITFFT=true .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the power and cross spectra estimates are
!                   normalized in such a way that the total area under the power spectra is equal to
!                   the variance of the time series contained in VEC and in each row of MAT.
!                 - NORMPSD = false, the sum of the PSD estimates
!                   (e.g. sum(PSVEC(2:)) and sum(PSMAT(:,2:),dim=2) ) is equal to the variance of
!                   the corresponding time series.
!
!                 The default is NORMPSD=true .
!
!   SMOOTH_PARAM  (INPUT, OPTIONAL) integer(i4b), dimension(:)
!                 if SMOOTH_PARAM is used, the power and cross spectra estimates are computed by repeated
!                 smoothing of the periodograms and cross-periodogram with modified Daniell weights.
!
!                 On entry, SMOOTH_PARAM(:) gives the array of the half-lengths of the
!                 modified Daniell filters to be applied.
!
!                 All the values in SMOOTH_PARAM(:) must be greater than 0 and less
!                 than (size(VEC)/2)+1 .
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The means of the time series are removed before computing the power and cross spectra
!                 - TREND=+2 The drifts from time series are removed before computing the power and cross spectra
!                 - TREND=+3 The least-squares lines from time series are removed before
!                   computing the power and cross spectra.
!
!                 For other values of TREND nothing is done before estimating the power and cross spectra.
!
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power and cross spectra. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ) and the critical value for
!                 testing the null hypothesis that the squared coherency is zero (e.g.
!                 the TESTCOHER optional argument).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the selected
!   data window (e.g. WIN=1,2,3,4,5,6) is applied to the time series and the PSD and CSD estimates are
!   computed by the FFT of these transformed time series. Optionally, theses PSD and CSD estimates
!   may then be smoothed in the frequency domain by modified Daniell filters (e.g. if argument
!   SMOOTH_PARAM is used).
!
!   The computed equivalent number of degrees of freedom and bandwidth must be divided by two
!   for the zero and Nyquist frequencies.
!
!   Furthermore, the computed equivalent number of degrees of freedom, bandwidth, lower and
!   upper (1-PROBTEST) * 100% confidence limit factors and critical value for the squared coherency
!   (e.g. arguments EDOF, BANDWIDTH, CONLWR, CONUPR and TESTCOHER) are not right near the zero
!   and Nyquist frequencies if the PSD estimates have been smoothed by modified Daniell filters.
!   The reason is that CROSS_SPECTRUM assumes that smoothing involves averaging independent frequency
!   ordinates. This is true except near the zero and Nyquist frequencies where an average
!   may contain contributions from negative frequencies, which are identical to and hence not
!   independent of positive frequency spectral values. Thus, the number of degrees of freedom
!   in PSD estimates near the 0 and Nyquist frequencies are as little as half the number
!   of degrees of freedom of the spectral estimates away from these frequency extremes if
!   the optional argument SMOOTH_PARAM is used.
!      
!   If the optional argument SMOOTH_PARAM is used, the computed equivalent number of degrees
!   of freedom, bandwidth, lower and upper (1-PROBTEST) * 100% confidence limit factors and
!   critical value for the squared coherency are right for PSD estimates  at frequencies
!      
!        (i-1)/size(VEC)  for i= (nparam+1)/2 + 1   to ( size(VEC) - nparam + 1)/2
!      
!   where nparam = 2 * (2+sum(SMOOTH_PARAM(:)))- 1,
!   (e.g. for frequencies i/size(VEC) for i = (nparam+1)/2, ... , (size(VEC)-nparam-1)/2 ) .
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!      
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, assert_eq, arth
    use Reals_Constants,   only : zero, half, one, two, four, pi, twopi, c5_m2
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error14, tseries_error15,    &
                                  tseries_error21, tseries_error50, tseries_error54, tseries_error59
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : probf2, pinvf2, pinvq2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:)   :: vec
    real(stnd), intent(inout),  dimension(:,:) :: mat
!
    real(stnd), intent(out), dimension(:)             :: psvec
    real(stnd), intent(out), dimension(:,:)           :: psmat, phase, coher
    real(stnd), intent(in),                  optional :: taperp, probtest
    real(stnd), intent(out),                 optional :: edof, bandwidth, conlwr, conupr, testcoher
    real(stnd), intent(out), dimension(:),   optional :: freq
    real(stnd), intent(out), dimension(:,:), optional :: ampli, co_spect, quad_spect, prob_coher
!
    integer(i4b),  intent(in),               optional :: trend, win
    integer(i4b),  intent(in), dimension(:), optional :: smooth_param
!
    logical(lgl),  intent(in), optional :: normpsd, initfft
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, n, nd2, nf, k, k1, k2, win2, trend2, nparam
    integer      :: iok
!
    real(stnd)                                :: c1, c2, sumw, edof2, probtest2, con
    real(stnd), dimension(size(vec))          :: wk
    real(stnd), dimension(:),   allocatable   :: temp
    real(stnd), dimension(:,:), allocatable   :: magni
!
    complex(stnd), dimension(:),   allocatable :: cwk1
    complex(stnd), dimension(:,:), allocatable :: cwk2
!
    logical(lgl)  :: normpsd2, initfft2, smooth, out_freq, out_ampli, out_co_spect,    &
                     out_quad_spect, out_prob_coher, out_dof
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='cross_spectrum'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m =  assert_eq( int(size(mat,1),i4b) ,          &
                    int(size(psmat,1),i4b) ,        &
                    int(size(phase,1),i4b) ,        &
                    int(size(coher,1),i4b) ,        &
                    name_proc )
!
    if ( m<=0_i4b )     &
    call merror( name_proc//tseries_error21 )
!
    n =  assert_eq( int(size(vec),i4b) ,          &
                    int(size(mat,2),i4b) ,        &
                    name_proc )
!
    if ( n<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    nd2 = n/2_i4b
!
    if ( n/=2*nd2 )   &
    call merror( name_proc//tseries_error14 )
!
    nf  = nd2 + 1_i4b
!
    call assert( logical(int(size(psvec),i4b)==nf,lgl),     &
                 logical(int(size(psmat,2),i4b)==nf,lgl),   &
                 logical(int(size(phase,2),i4b)==nf,lgl),   &
                 logical(int(size(coher,2),i4b)==nf,lgl),   &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_ampli = present( ampli )
!
    if ( out_ampli ) then
        call assert( logical(int(size(ampli,1),i4b)==m,lgl),   &
                     logical(int(size(ampli,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_co_spect = present( co_spect )
!
    if ( out_co_spect ) then
        call assert( logical(int(size(co_spect,1),i4b)==m,lgl),   &
                     logical(int(size(co_spect,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_quad_spect = present( quad_spect )
!
    if ( out_quad_spect ) then
        call assert( logical(int(size(quad_spect,1),i4b)==m,lgl),   &
                     logical(int(size(quad_spect,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_prob_coher = present( prob_coher )
!
    if ( out_prob_coher ) then
        call assert( logical(int(size(prob_coher,1),i4b)==m,lgl),   &
                     logical(int(size(prob_coher,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth ) .or. out_prob_coher   &
              .or. present( conlwr ) .or. present( conupr ) .or. present( testcoher ) 
!
    probtest2 = c5_m2
!
    if ( present(probtest) ) then
        if ( zero>=probtest .or. probtest>=one ) then
            call merror( name_proc//tseries_error59 )
        else
            probtest2 = probtest
        end if
    end if
!
    if ( present(trend) ) then
        trend2 = trend
    else
        trend2 = 1_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
    else
        win2 = 3_i4b
    end if
!
    initfft2 = true
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( (/ m, nd2/), dim=2_i4b )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    smooth = false
!
    if ( present(smooth_param) ) then
!
        nparam = size( smooth_param )
!
        if ( nparam<=0_i4b  )    &
        call merror( name_proc//tseries_error50 )
!
!       CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIEL FILTERS.
!
        if ( any( smooth_param(:nparam)<=0_i4b .or. smooth_param(:nparam)>=nf ) )   &
        call merror( name_proc//tseries_error54 )
!
        smooth = true
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trend2>=1_i4b .and. trend2<=3_i4b ) then
        call detrend_rv( vec(:n),    trend2 )
        call detrend_rm( mat(:m,:n), trend2 )
    end if
!
!   TAPER THE TIME SERIES IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:n)  = data_window( n, win2, taperp=taperp )
!
!       APPLY DATA WINDOW.
!
        vec(:n) = vec(:n)*wk(:n)
!
        mat(:m,:n) = mat(:m,:n)*spread( wk(:n), dim=1, ncopies=m )
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:n), wk(:n) )
!
    else    
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:n)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( n, stnd )
!
    end if
!
!   COMPUTE DEGREES OF FREEDOM IF REQUIRED.
!
    if ( out_dof ) then
!
        edof2 = estim_dof( wk(:n), smooth_param=smooth_param(:nparam) )
!
    end if
!
!   COMPUTE PSD AND CSD ESTIMATES.
!
!   ALLOCATE WORK MATRIX.
!
    allocate( cwk1(nf), cwk2(m,nf), stat = iok )
!
    if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!   COMPUTE FFT FOR THE SERIES.
!
    call real_fft( vec(:n), cwk1(:nf), true)
!
    call real_fft( mat(:m,:n), cwk2(:m,:nf), true )
!
!   COMPUTE POWER SPECTRUM ESIMATES.
!
    psvec(:nd2) = real(cwk1(:nd2),stnd)**2 + aimag(cwk1(:nd2))**2
    psvec(nf)   = real(cwk1(nf),stnd)**2
!
    psmat(:m,:nd2) = real(cwk2(:m,:nd2),stnd)**2  + aimag(cwk2(:m,:nd2))**2
    psmat(:m,nf)   = real(cwk2(:m,nf),stnd)**2
!
!   COMPUTE CO-SPECTRUM ESIMATES.
!
    coher(:m,:nd2)  = spread(real(cwk1(:nd2),stnd),dim=1,ncopies=m)*real(cwk2(:m,:nd2),stnd) +   &
                      spread(aimag(cwk1(:nd2)),    dim=1,ncopies=m)*aimag(cwk2(:m,:nd2))
    coher(:m,nf)    = real(cwk1(nf),stnd)*real(cwk2(:m,nf),stnd)
!
!   COMPUTE QUADRATURE-SPECTRUM ESIMATES.
!
    phase(:m,:nd2)  = spread(aimag(cwk1(:nd2)),    dim=1,ncopies=m)*real(cwk2(:m,:nd2),stnd) -   &
                      spread(real(cwk1(:nd2),stnd),dim=1,ncopies=m)*aimag(cwk2(:m,:nd2))
    phase(:m,nf)    = zero
!
!   DEALLOCATE WORK ARRAYS.
!
    deallocate( cwk1, cwk2 )
!
!   NORMALIZE THE POWER AND CROSS SPECTRA ESTIMATES.
!
    c1 = one/(sumw*real( n, stnd))
!
    if ( normpsd2 ) then
!
!       CHATFIELD DEFINITION OF THE SPECTRAL DENSITY FUNCTION.
!
        c2 = (real( n, stnd)/pi)*c1
        k1 = 1_i4b
        k2 = nf
!
    else
!
!       POLLOCK DEFINITION OF THE PERIODOGRAM. IF THE ESTIMATES
!       ARE DIVIDED BY 4*PI, THIS IS THE POLLOCK DEFINITION OF
!       THE SPECTRAL DENSITY FUNCTION.
!
        c2 = two*c1
        k1 = 2_i4b
        k2 = nd2
!
        psvec(1_i4b)    = c1*psvec(1_i4b)
        psmat(:m,1_i4b) = c1*psmat(:m,1_i4b)
!
        coher(:m,1_i4b)  = c1*coher(:m,1_i4b)
        phase(:m,1_i4b)  = c1*phase(:m,1_i4b)
!
        psvec(nf)    = c1*psvec(nf)
        psmat(:m,nf) = c1*psmat(:m,nf)
!
        coher(:m,nf)  = c1*coher(:m,nf)
        phase(:m,nf)  = c1*phase(:m,nf)
!
    end if
!
    psvec(k1:k2)    = c2*psvec(k1:k2)
    psmat(:m,k1:k2) = c2*psmat(:m,k1:k2)
!
    coher(:m,k1:k2)   = c2*coher(:m,k1:k2)
    phase(:m,k1:k2)   = c2*phase(:m,k1:k2)
!
!   SMOOTH THE POWER AND CROSS SPECTRA ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
!       SMOOTH POWER SPECTRUM OF THE FIRST SERIES.
!
        call moddan_filter_rv( psvec(:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH POWER SPECTRUM OF THE OTHER SERIES.
!
        call moddan_filter_rm( psmat(:m,:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH CO-SPECTRUM.
!
        call moddan_filter_rm( coher(:m,:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH QUADRATURE-SPECTRUM.
!
        call moddan_filter_rm( phase(:m,:nf), smooth_param(:nparam), sym=-one )
!
    end if
!
    if ( out_co_spect ) then
!
!       OUTPUT THE CO-SPECTRUM IF REQUIRED.
!
        co_spect(:m,:nf) = coher(:m,:nf)
!
    end if
!
    if ( out_quad_spect ) then
!
!       OUTPUT THE QUADRATURE-SPECTRUM IF REQUIRED.
!
        quad_spect(:m,:nf) = -phase(:m,:nf)
!
    end if
!
!   ALLOCATE WORK MATRICES.
!
    allocate( temp(m), magni(m,nf), stat = iok )
!
    if ( iok/=0 ) call merror( name_proc//allocate_error )
!
    do k = 1_i4b, nf
!
!       COMPUTE MAGNITUDE OF THE CROSS-SPECTRUM.
!
        magni(:m,k) = coher(:m,k)**2 + phase(:m,k)**2
!
!       COMPUTE PHASE OF THE CROSS-SPECTRUM.
!
        where ( coher(:m,k)/=zero .or. phase(:m,k)/=zero )
            temp(:m)   = atan2( phase(:m,k), coher(:m,k) )
            phase(:m,k) = (one/twopi)*temp(:m) + merge( zero, one, temp(:m)>=zero ) 
        end where
!
!       COMPUTE THE SQUARED COHERENCY.
!
        if ( psvec(k)/=zero ) then
!
            where( psmat(:m,k)/=zero )
                coher(:m,k) = min( magni(:m,k)/(psvec(k)*psmat(:m,k)), one )
            elsewhere
                coher(:m,k) = zero
            end where
!
        else
!
            coher(:m,k) = zero
!
        end if
!
    end do   
!
    if ( out_ampli ) then
!
!       OUTPUT AMPLITUDE OF THE CROSS-SPECTRUM IF REQUIRED.
!
        ampli(:m,:nf) = sqrt( magni(:m,:nf) )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
!
        c1 = one/real( n, stnd )
        freq(:nf) = arth(zero, c1, nf)
!
    end if
!
!   OUTPUT DEGREES OF FREEDOM, BANDWIDTH, CONFIDENCE LIMIT FACTORS 
!   OR PROBABILITIES FOR SQUARED COHERENCIES IF REQUIRED.
!
    if ( out_dof ) then
!
        if ( present(edof) ) then
            edof = edof2
        end if
!
        if ( present(bandwidth) ) then
            bandwidth = edof2/( two*real( n, stnd ) )
        end if
!
        c1 = probtest2*half
        c2 = one - c1
!
        sumw = nan()
!
        if ( present(conlwr) ) then
            if ( edof2>=half ) then
                conlwr = edof2/pinvq2( c2, edof2 )
            else
                conlwr = sumw
            end if
        end if
!
        if ( present(conupr) ) then
            if ( edof2>=half ) then
                conupr = edof2/pinvq2( c1, edof2 )
            else
                conupr = sumw
            end if
        end if
!
        if ( present(testcoher) ) then
!
            if ( edof2>two ) then
!
!                c1 = two/(edof2-two)
!                c2 = one - probtest2**c1
!                testcoher = max( min( one, c2 ), zero )
!
                 con = one - probtest2
                 c2  = edof2 - two
                 c1  = two*pinvf2( con, two, c2 )
                 testcoher = c1/( c2 + c1 )
            else
                testcoher = sumw
            end if
!
        end if
!
!       COMPUTE AND OUTPUT SIGNIFICANCE PROBABILITY FOR SQUARED COHERENCIES.
!
        if ( out_prob_coher ) then
!
            c2 = edof2 - two    
!
            if ( c2>zero ) then
!
                where( coher(:m,2_i4b:nd2)/=one )
                    magni(:m,2_i4b:nd2) = ( (edof2/two - one)*coher(:m,2_i4b:nd2) )/  &
                                          ( one - coher(:m,2_i4b:nd2) )
                elsewhere
                    magni(:m,2_i4b:nd2) = zero
                end where 
!            
                prob_coher(:m,2_i4b:nd2) = probf2( magni(:m,2_i4b:nd2), two, c2, true )
!            
                where( coher(:m,2_i4b:nd2)==one ) prob_coher(:m,2_i4b:nd2) = zero
!
            else
!
                prob_coher(:m,2_i4b:nd2) = sumw
!
            end if
!            
            c2 = edof2/two - two    
!
            if ( c2>zero ) then
!            
                do k1 = 1_i4b, nf, nd2
!            
                    where( coher(:m,k1)/=one )
                        temp(:m) = ( (edof2/four - one)*coher(:m,k1) )/( one - coher(:m,k1) )
                    elsewhere
                        temp(:m) = zero
                    end where
!            
                    prob_coher(:m,k1) = probf2( temp(:m), two, c2, true )
!            
                    where( coher(:m,k1)==one ) prob_coher(:m,k1) = zero
!
                end do
!
            else
!
                prob_coher(:m,1_i4b) = sumw
                prob_coher(:m,nf)    = sumw
!
            end if
!
        end if
!
    end if
!
!   DEALLOCATE WORK ARRAYS.
!
    deallocate( temp, magni )
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE cross_spectrum_rm
! ___________________________________
!
    end subroutine cross_spectrum_rm
!
! =========================================================================================
!
    subroutine power_spectrum2_rv( vec, l, psvec, freq, edof, bandwidth, conlwr, conupr, initfft,  &
                                   overlap, normpsd, smooth_param, trend, trend2, win, taperp, l0, &
                                   probtest )
!
! Purpose
! _______
!
!   Subroutine POWER_SPECTRUM2 computes a Fast Fourier Transform (FFT) estimate
!   of the power spectrum of a real time series.
!
!   The Power Spectral Density (PSD) estimates are returned in units which are
!   the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the real time series for which the power spectrum
!                 must be estimated.
!                 If TREND=1, 2 or 3,  VEC is used as workspace and is transformed.
!
!                 Size(VEC) must be greater or equal to 4.
!
!   L             (INPUT) integer(i4b)
!                 On entry, an integer used to segment the time series. L is the
!                 length of the segments. L must be a positive even integer, less
!                 or equal to size(VEC), but greater or equal to 4.
!
!                 Spectral computations are at (L/2)+1 frequencies if the optional
!                 argument L0 is absent and are at ((L+L0)/2)+1 frequencies if L0 is
!                 present (L0 is the number of zeros added to each segment).
!
!                 Suggested values for L+L0 are 16, 32, 64 or 128 (e.g. an integer power of two,
!                 in order to speed the computations).
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = ((L+L0)/2) + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/(L+L0)
!                 for i=1,2, ... , ((L+L0)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = (L+L0)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the bandwidth of the power spectrum estimates.
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the
!                 PSVEC(:) argument) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine POWER_SPECTRUM2 in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine POWER_SPECTRUM2. This call to INITFFT must have the following form: 
!
!                       call init_fft( (L+L0)/2 )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   POWER_SPECTRUM2 and a call to END_FFT is also done before leaving
!                   subroutine POWER_SPECTRUM2.
!
!                 The default is INITFFT=true .
!
!   OVERLAP       (INPUT, OPTIONAL) logical(lgl)
!                 If:
!
!                 - OVERLAP = false, the subroutine segments the data 
!                   without any overlapping.
!                 - OVERLAP = true, the subroutine overlaps the segments
!                   by one half of their length (which is equal to L).
!
!                 In both cases, zeros are eventually added to each segment (if argument L0 is present)
!                 and each segment will be FFT'd, and the resulting periodograms
!                 will averaged together to obtain a Power Spectrum Density estimate at the
!                 ((L+L0)/2)+1 frequencies.
!
!                 The default is OVERLAP=false .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the PSD estimates are normalized in such
!                   a way that the total area under the power spectrum is equal to the variance of
!                   the time series VEC.
!                 - NORMPSD = false, the sum of the PSD estimates
!                   (e.g. sum( PSVEC(2:) ) is equal to the variance of the time series.
!
!                 The default is NORMPSD=true .
!
!   SMOOTH_PARAM  (INPUT, OPTIONAL) integer(i4b), dimension(:)
!                 If SMOOTH_PARAM is used, the PSD estimates are computed by repeated
!                 smoothing of the periodogram with modified Daniell weights.
!
!                 On entry, SMOOTH_PARAM(:) gives the array of the half-lengths of the
!                 modified Daniell filters to be applied.
!
!                 All the values in SMOOTH_PARAM(:) must be greater than 0 and less
!                 than ((L+L0)/2) + 1 .
!
!                 Size(SMOOTH_PARAM) must be greater or equal to 1.
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The mean of the time series is removed before computing the spectrum
!                 - TREND=+2 The drift from the time series is removed before computing the spectrum
!                   by using the formula: drift = (VEC(size(VEC)) - VEC(1))/(size(VEC) - 1)
!                 - TREND=+3 The least-squares line from the time series is removed before
!                   computing the spectrum.
!
!                 For other values of TREND nothing is done before estimating the power spectrum.
!
!                 The default is TREND=1, e.g. the mean of the time series is removed before the
!                 computations.
!
!   TREND2        (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND2=+1 The mean of the time segment is removed before computing the spectrum
!                   on this segment.
!                 - TREND2=+2 The drift from the time segment is removed before computing the spectrum
!                   on this segment.
!                 - TREND2=+3 The least-squares line from the time segment is removed before
!                   computing the spectrum on this segment.
!
!                 For other values of TREND2 nothing is done before estimating the power spectrum
!                 on each segment.
!
!                 The default is TREND2=0, e.g. nothing is done before estimating the power spectrum
!                 on each segment.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power spectrum. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!
!   L0            (INPUT, OPTIONAL) integer(i4b)
!                 The number of zeros added to each time segment in order to obtain more finely
!                 spaced spectral estimates. L+L0 must be a positive even integer.
!
!                 The default is L0=0, e.g. no zeros are added to each time segment.
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the series
!   is padded with zero on the right such that the length of the resulting time series is 
!   evenly divisible by L (a positive even integer). The length, N, of this resulting time
!   series is the first integer greater than or equal to size(VEC) which is evenly divisible
!   by L. If size(VEC) is not evenly divisible by L, N is equal to size(VEC)+L-mod(size(VEC),L).
!   
!   Optionally, the mean or the trend may also be removed from each time segment (e.g. TREND2=1,2,3).
!   Optionally, zeros may be added to each time segment (e.g. the optional arguemnt L0) if more
!   finely spaced spectral esimates are desired.
!   
!   The stability of the PSD estimates depends on the averaging process. That is, the greater
!   the number of segments ( N/L if OVERLAP=false and (2N/L)-1 if OVERLAP=true), the more
!   stable the resulting PSD estimates.
!   
!   Optionally, theses PSD estimates may then be smoothed again in the frequency
!   domain by modified Daniell filters (e.g. if argument SMOOTH_PARAM is used).
!
!   The computed equivalent number of degrees of freedom and bandwidth must be divided by two
!   for the zero and Nyquist frequencies.
!
!   Furthermore, the computed equivalent number of degrees of freedom, bandwidth, lower and
!   upper (1-PROBTEST) * 100% confidence limit factors are not right near the zero and Nyquist frequencies
!   if the PSD estimates have been smoothed by modified Daniell filters. The reason is that
!   POWER_SPECTRUM2 assumes that smoothing involves averaging independent frequency
!   ordinates. This is true except near the zero and Nyquist frequencies where an average
!   may contain contributions from negative frequencies, which are identical to and hence not
!   independent of positive frequency spectral values. Thus, the number of degrees of freedom
!   in PSD estimates near the 0 and Nyquist frequencies are as little as half the number
!   of degrees of freedom of the spectral estimates away from these frequency extremes if
!   the optional argument SMOOTH_PARAM is used.
!      
!   If the optional argument SMOOTH_PARAM is used, the computed equivalent number of degrees
!   of freedom, bandwidth, lower and upper (1-PROBTEST) * 100% confidence limit factors are
!   right for PSD estimates at frequencies
!      
!              (i-1)/(L+L0)  for i= (nparam+1)/2 + 1  to  ( (L+L0) - nparam + 1)/2
!      
!   where nparam = 2 * (2+sum(SMOOTH_PARAM(:)))- 1,
!   (e.g. for frequencies i/(L+L0) for i = (nparam+1)/2, ... , ( (L+L0) - nparam - 1)/2   ).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!      
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, arth
    use Reals_Constants,   only : zero, half, one, two, pi, c5_m2
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error15, tseries_error17,    &
                                  tseries_error18, tseries_error24, tseries_error25, tseries_error50,   &
                                  tseries_error56, tseries_error58, tseries_error59
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : pinvq2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout), dimension(:) :: vec
    real(stnd), intent(out),   dimension(:) :: psvec
!
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out),               optional :: edof, bandwidth, conlwr, conupr
    real(stnd), intent(out), dimension(:), optional :: freq
!
    integer(i4b),  intent(in)                         :: l
    integer(i4b),  intent(in),               optional :: trend, trend2, win, l0
    integer(i4b),  intent(in), dimension(:), optional :: smooth_param
!
    logical(lgl),  intent(in), optional :: initfft, overlap, normpsd
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, n2, nf, k, l2, ld2, m, win2, trendb, trend2b,    &
                    i, i1, i2, step, sl, sld2, nparam
    integer      :: iok
!
    real(stnd)                            :: c1, c2, sumw, edof2, probtest2
    real(stnd), dimension(l)              :: wk
    real(stnd), dimension(:), allocatable :: seg
!
    complex(stnd), dimension(:), allocatable  :: cwk
!
    logical(lgl)  :: normpsd2, overlap2, initfft2, smooth, out_freq, out_dof
!
#ifdef _OPENMP
!
    real(stnd), dimension(:), allocatable :: psvecb
!
    logical  :: test_par
!
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='power_spectrum2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n   = size( vec )
!
    if ( l<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    if ( l>n )     &
    call merror( name_proc//tseries_error17 )
!
!   CHECK IF l IS A POSITIVE EVEN INTEGER.
!
    ld2 = l/2_i4b
!
    if ( l/=2_i4b*ld2 )     &
    call merror( name_proc//tseries_error18 )
!
!   DETERMINE THE LENGTH OF THE SEGMENTS.
!
    if ( present(l0) ) then
!
        if ( l0<0_i4b )     &
        call merror( name_proc//tseries_error24 )
!
        sl = l + l0
!
!       CHECK IF sl IS A POSITIVE EVEN INTEGER.
!
        sld2 = sl/2_i4b
!
        if ( sl/=2_i4b*sld2 )     &
        call merror( name_proc//tseries_error25  )
!
    else
!
        sl   = l
        sld2 = ld2
!
    end if
!
    nf  = sld2 + 1_i4b
    call assert( logical(int(size(psvec),i4b)==nf,lgl),    &
                 name_proc )
!
    out_freq = present( freq )
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth )        &
              .or. present( conlwr ) .or. present( conupr ) 
!
    probtest2 = c5_m2
    if ( present(probtest) ) then
        if ( zero>=probtest .or. probtest>=one ) then
            call merror( name_proc//tseries_error59 )
        else
            probtest2 = probtest
        end if
    end if
!
    if ( present(trend) ) then
        trendb = trend
    else
        trendb = 1_i4b
    end if
!
    if ( present(trend2) ) then
        trend2b = trend2
    else
        trend2b = 0_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
    else
        win2 = 3_i4b
    end if
!
    initfft2 = true
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( sld2 )
    end if
!
    normpsd2 = true
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    overlap2 = false
    if ( present(overlap) ) then
        overlap2 = overlap
    end if
!
    smooth = false
    if ( present(smooth_param) ) then
!
        nparam = size( smooth_param )
!
        if ( nparam<=0_i4b  )    &
        call merror( name_proc//tseries_error50 )
!
!       CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIEL FILTERS.
!
        if ( any( smooth_param(:nparam)<=0_i4b .or. smooth_param(:nparam)>=nf ) )   &
        call merror( name_proc//tseries_error56 )
!
        smooth = true
!
    end if
!
!   FIND THE FIRST INTEGER GREATER THAN OR EQUAL TO n THAT IS EVENLY
!   DIVISIBLE BY l.
!
    i = mod( n, l )
!
    if ( i==0_i4b ) then
        n2 = n
    else
        n2 = n + l - i
    end if
!
!   DETERMINE THE NUMBER OF SEGMENTS.
!
    if ( overlap2 ) then
        m    = (2*n2/l) - 1_i4b
        step = ld2
    else
        m    = n2/l
        step = l
    end if
!
!   CALCULATE DATA WINDOW IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:l)  = data_window( l, win2, taperp=taperp )
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:l), wk(:l) )
!
    else
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:l)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( l, stnd )
!
    end if
!
!   COMPUTE DEGREES OF FREEDOM IF REQUIRED.
!
    if ( out_dof ) then
!
        if ( m/=1_i4b .and. overlap2 .and. win2/=1_i4b .and. win2/=3_i4b  )    &
        call merror( name_proc//tseries_error58 )
!
        edof2 = estim_dof( wk(:l), win=win, smooth_param=smooth_param(:nparam), &
                           l0=l0, nseg=m, overlap=overlap )
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trendb>=1_i4b .and. trendb<=3_i4b ) then
        call detrend_rv( vec(:n), trendb )
    end if
!
!   ZERO OUTPUT PSD ESTIMATES FOR SUMMING OVER THE m SEGMENTS.
!
    psvec(:nf) = zero
!
!   COMPUTE PSD ESTIMATES BY SUMMING OVER THE m SEGMENTS.
!
#ifdef _OPENMP
    i1 = omp_get_num_procs()
    i2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               i1>1_i4b                     .and.      &
               i2>1_i4b                     .and.      &
               m>=i2
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(k,i,i1,i2,l2,iok,psvecb,seg,cwk)
!
!       ALLOCATE WORK VECTORS.
!
        allocate( psvecb(nf), seg(sl), cwk(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ZERO PARTIAL OUTPUT PSD ESTIMATES FOR SUMMING OVER EACH THREAD.
!
        psvecb(:nf) = zero
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(l+1_i4b:sl) = zero
        end if    
!
!$OMP DO SCHEDULE(STATIC)
!
        do k=1_i4b, m
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg(:l2) = vec(i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENT IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rv( seg(:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:l2) = seg(:l2)*wk(:l2)
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT OF THE TIME SEGMENT.
!
            call real_fft( seg(:sl), cwk(:nf), true)
!
!           UPDATE SPECTRAL ESTIMATES.
!
            psvecb(:sld2) = psvecb(:sld2) + real(cwk(:sld2),stnd)**2 + aimag(cwk(:sld2))**2
            psvecb(nf)    = psvecb(nf)    + real(cwk(nf),stnd)**2
!        
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (updatepsd)
        psvec(:nf) = psvec(:nf) + psvecb(:nf)
!$OMP END CRITICAL (updatepsd)
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( psvecb, seg, cwk )
!
!$OMP END PARALLEL
!
    else
!
#endif
!
!       ALLOCATE WORK VECTORS.
!
        allocate( seg(sl), cwk(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(l+1_i4b:sl) = zero
        end if    
!
        do k=1_i4b, m
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
            l2 = i2 - i1 + 1_i4b
!
            seg(:l2) = vec(i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENT IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rv( seg(:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:l2) = seg(:l2)*wk(:l2)
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT OF THE TIME SEGMENT.
!
            call real_fft( seg(:sl), cwk(:nf), true)
!
            psvec(:sld2) = psvec(:sld2) + real(cwk(:sld2),stnd)**2 + aimag(cwk(:sld2))**2
            psvec(nf)    = psvec(nf)    + real(cwk(nf),stnd)**2
!        
        end do
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( seg, cwk )
!        
#ifdef _OPENMP
!
    end if
!
#endif
!
!   NORMALIZE THE PSD ESTIMATES.
!
    c1 = one/( sumw*real(l*m, stnd) )
!
    if ( normpsd2 ) then
!
!       CHATFIELD DEFINITION OF THE SPECTRAL DENSITY FUNCTION.
!
        c2 = real( l, stnd)/pi
        psvec(1_i4b:nf) = (c1*c2)*psvec(1_i4b:nf)
!
    else
!
!       POLLOCK DEFINITION OF THE PERIODOGRAM. IF THE ESTIMATES
!       ARE DIVIDED BY 4*PI, THIS IS THE POLLOCK DEFINITION OF
!       THE SPECTRAL DENSITY FUNCTION.
!
        psvec(1_i4b)      = c1*psvec(1_i4b)
        psvec(2_i4b:sld2) = (two*c1)*psvec(2_i4b:sld2)
        psvec(nf)         = c1*psvec(nf)
!
    end if
!
!   SMOOTH THE PSD ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
        call moddan_filter_rv( psvec(:nf), smooth_param(:nparam), sym=one )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( sl, stnd )
        freq(:nf) = arth(zero, c1, nf)
    end if
!
!   OUTPUT DEGREES OF FREEDOM OR BANDWIDTH IF REQUIRED.
!
    if ( out_dof ) then
!
        if ( present(edof) ) then
            edof = edof2
        end if
!
        if ( present(bandwidth) ) then
            bandwidth = edof2/( two*real( n, stnd ) )
        end if
!
        c1 = probtest2*half
        c2 = one - c1
!
        if ( present(conlwr) ) then
            if ( edof2>=half ) then
                conlwr = edof2/pinvq2( c2, edof2 )
            else
                conlwr = nan()
            end if
        end if
!
        if ( present(conupr) ) then
            if ( edof2>=half ) then
                conupr = edof2/pinvq2( c1, edof2 )
            else
                conupr = nan()
            end if
        end if
!
    end if
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE power_spectrum2_rv
! ____________________________________
!
    end subroutine power_spectrum2_rv
!
! =========================================================================================
!
    subroutine power_spectrum2_rm( mat, l, psmat, freq, edof, bandwidth, conlwr, conupr, initfft,  &
                                   overlap, normpsd, smooth_param, trend, trend2, win, taperp, l0, &
                                   probtest )
!
! Purpose
! _______
!
!   Subroutine POWER_SPECTRUM2 computes Fast Fourier Transform (FFT) estimates
!   of the power spectra of the multi-channel real time series MAT (e.g. each row
!   of MAT contains a time series).
!
!   The Power Spectral Density (PSD) estimates are returned in units which are
!   the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   MAT           (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                 On entry, the multi-channel real time series for which the power 
!                 spectra must be estimated. Each row of MAT is a real time series.
!                 If TREND=1, 2 or 3,  MAT is used as workspace and is transformed.
!
!                 Size(MAT,2) must be greater or equal to 4.
!
!   L             (INPUT) integer(i4b)
!                 On entry, an integer used to segment the time series. L is the
!                 length of the segments. L must be a positive even integer less
!                 or equal to size(MAT,2), but greater or equal to 4.
!
!                 Spectral computations are at (L/2)+1 frequencies if the optional
!                 argument L0 is absent and are at ((L+L0)/2)+1 frequencies if L0 is
!                 present (L0 is the number of zeros added to each segment).
!
!                 Suggested values for L+L0 are 16, 32, 64 or 128 (e.g. an integer power of two,
!                 in order to speed the computations).
!
!   PSMAT         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the Power Spectral Density (PSD) estimates of each row of MAT.
!
!                 The shape of PSMAT must verify:
!
!                 - size(PSMAT,1) =  size(MAT,1) ;
!                 - size(PSMAT,2) = ((L+L0)/2) + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/(L+L0)
!                 for i=1,2, ... , ((L+L0)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = (L+L0)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 spectrum estimates.
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the bandwidth of the power spectrum estimates.
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the
!                 PSMAT(:,:) argument) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine POWER_SPECTRUM2 in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine POWER_SPECTRUM2. This call to INITFFT must have the following form: 
!
!                       call init_fft( (/ size(MAT,1), (L+L0)/2 /), dim=2_i4b )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   POWER_SPECTRUM2 and a call to END_FFT is also done before leaving
!                   subroutine POWER_SPECTRUM2.
!
!                 The default is INITFFT=true .
!
!   OVERLAP       (INPUT, OPTIONAL) logical(lgl)
!                 If:
!
!                 - OVERLAP = false, the subroutine segments the data 
!                   without any overlapping.
!                 - OVERLAP = true, the subroutine overlaps the segments
!                   by one half of their length (which is equal to L).
!
!                 In both cases, zeros are eventually added to each segment (if argument L0 is present)
!                 and each segment will be FFT'd, and the resulting periodograms
!                 will averaged together to obtain a Power Spectrum Density estimate at the
!                 ((L+L0)/2)+1 frequencies.
!
!                 The default is OVERLAP=false .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the PSD estimates are normalized in such
!                   a way that the total area under the power spectrum is equal to the variance of
!                   the corresponding time series in MAT.
!                 - NORMPSD = false, the sum of the PSD
!                   estimates (e.g. sum(PSMAT(:,2:),dim=2) ) is equal to the variance of the corresponding time series.
!
!                 The default is NORMPSD=true .
!
!   SMOOTH_PARAM  (INPUT, OPTIONAL) integer(i4b), dimension(:)
!                 If SMOOTH_PARAM is used, the PSD estimates are computed by repeated
!                 smoothing of the periodogram with modified Daniell weights.
!
!                 On entry, SMOOTH_PARAM(:) gives the array of the half-lengths of the
!                 modified Daniell filters to be applied.
!
!                 All the values in SMOOTH_PARAM(:) must be greater than 0 and less
!                 than ((L+L0)/2)+1 .
!
!                 Size(SMOOTH_PARAM) must be greater or equal to 1.
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The means of the time series are removed before computing the spectra
!                 - TREND=+2 The drifts from time series are removed before computing the spectra
!                 - TREND=+3 The least-squares lines from time series are removed before
!                   computing the spectra.
!
!                 For other values of TREND nothing is done before estimating the power and cross spectra.
!
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   TREND2        (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND2=+1 The mean of the time segment is removed before computing the spectrum
!                   on this segment.
!                 - TREND2=+2 The drift from the time segment is removed before computing the spectrum
!                   on this segment.
!                 - TREND2=+3 The least-squares line from the time segment is removed before
!                   computing the spectrum on this segment.
!
!                 For other values of TREND2 nothing is done before estimating the power spectrum
!                 on each segment.
!
!                 The default is TREND2=0, e.g. nothing is done before estimating the power spectrum
!                 on each segment.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power spectrum. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!
!   L0            (INPUT, OPTIONAL) integer(i4b)
!                 The number of zeros added to each time segment in order to obtain more finely
!                 spaced spectral estimates. L+L0 must be a positive even integer.
!
!                 The default is L0=0, e.g. no zeros are added to each time segment.
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the series
!   are padded with zero on the right such that the length of the resulting time series is 
!   evenly divisible by L (a positive even integer). The length, N, of this resulting time
!   series is the first integer greater than or equal to size(MAT,2) which is evenly divisible
!   by L. If size(MAT,2) is not evenly divisible by L, N is equal to size(MAT,2)+L-mod(size(MAT,2),L).
!   
!   Optionally, the mean or the trend may also be removed from each time segment (e.g. TREND2=1,2,3).
!   Optionally, zeros may be added to each time segment (e.g. the optional arguemnt L0) if more
!   finely spaced spectral esimates are desired.
!   
!   The stability of the PSD estimates depends on the averaging process. That is, the greater
!   the number of segments ( N/L if OVERLAP=false and (2N/L)-1 if OVERLAP=true), the more
!   stable the resulting PSD estimates.
!   
!   Optionally, theses PSD estimates may then be smoothed again in the frequency
!   domain by modified Daniell filters (e.g. if argument SMOOTH_PARAM is used).
!
!   The computed equivalent number of degrees of freedom and bandwidth must be divided by two
!   for the zero and Nyquist frequencies.
!
!   Furthermore, the computed equivalent number of degrees of freedom, bandwidth, lower and
!   upper (1-PROBTEST) * 100% confidence limit factors are not right near the zero and Nyquist frequencies
!   if the PSD estimates have been smoothed by modified Daniell filters. The reason is that
!   POWER_SPECTRUM2 assumes that smoothing involves averaging independent frequency
!   ordinates. This is true except near the zero and Nyquist frequencies where an average
!   may contain contributions from negative frequencies, which are identical to and hence not
!   independent of positive frequency spectral values. Thus, the number of degrees of freedom
!   in PSD estimates near the 0 and Nyquist frequencies are as little as half the number
!   of degrees of freedom of the spectral estimates away from these frequency extremes if
!   the optional argument SMOOTH_PARAM is used.
!      
!   If the optional argument SMOOTH_PARAM is used, the computed equivalent number of degrees
!   of freedom, bandwidth, lower and upper (1-PROBTEST) * 100% confidence limit factors are
!   right for PSD estimates at frequencies
!      
!              (i-1)/(L+L0)  for i= (nparam+1)/2 + 1  to  ( (L+L0) - nparam + 1)/2
!      
!   where nparam = 2 * (2+sum(SMOOTH_PARAM(:)))- 1,
!   (e.g. for frequencies i/(L+L0) for i = (nparam+1)/2, ... , ( (L+L0) - nparam - 1)/2   ).
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, assert_eq, merror, arth
    use Reals_Constants,   only : zero, half, one, two, pi, c5_m2
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error15, tseries_error18,    &
                                  tseries_error19, tseries_error24, tseries_error25, tseries_error50,   &
                                  tseries_error56, tseries_error58, tseries_error59
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : pinvq2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:,:) :: mat
    real(stnd), intent(out),    dimension(:,:) :: psmat
!
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out),               optional :: edof, bandwidth, conlwr, conupr
    real(stnd), intent(out), dimension(:), optional :: freq
!
    integer(i4b),  intent(in)                         :: l
    integer(i4b),  intent(in),               optional :: trend, trend2, win, l0
    integer(i4b),  intent(in), dimension(:), optional :: smooth_param
!
    logical(lgl),  intent(in), optional :: initfft, overlap, normpsd
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, n, n2, nf, k, l2, ld2, nseg, win2, trendb, trend2b,   &
                    i, i1, i2, step, sl, sld2, nparam
    integer      :: iok
!
    real(stnd)                              :: c1, c2, sumw, edof2, probtest2
    real(stnd), dimension(l)                :: wk
    real(stnd), dimension(:,:), allocatable :: seg
!
    complex(stnd), dimension(:,:), allocatable :: cwk
!
    logical(lgl)  :: normpsd2, initfft2, overlap2, smooth, out_freq, out_dof
!
#ifdef _OPENMP
!
    real(stnd), dimension(:,:), allocatable :: psmatb
!
    logical  :: test_par
!
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='power_spectrum2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m =  assert_eq( int(size(mat,1),i4b) ,          &
                    int(size(psmat,1),i4b) ,        &
                    name_proc )
    if ( m<=0_i4b ) return
!
    n =  size( mat, 2 )
!
    if ( l<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    if ( l>n )     &
    call merror( name_proc//tseries_error19 )
!
!   CHECK IF l IS A POSITIVE EVEN INTEGER.
!
    ld2 = l/2_i4b
!
    if ( l/=2_i4b*ld2 )     &
    call merror( name_proc//tseries_error18 )
!
!   DETERMINE THE LENGTH OF THE SEGMENTS.
!
    if ( present(l0) ) then
!
        if ( l0<0_i4b )     &
        call merror( name_proc//tseries_error24 )
!
        sl = l + l0
!
!       CHECK IF sl IS A POSITIVE EVEN INTEGER.
!
        sld2 = sl/2_i4b
!
        if ( sl/=2_i4b*sld2 )     &
        call merror( name_proc//tseries_error25 )
!
    else
!
        sl   = l
        sld2 = ld2
!
    end if
!
    nf  = sld2 + 1_i4b
    call assert( logical(int(size(psmat,2),i4b)==nf,lgl),    &
                 name_proc )
!
    out_freq = present( freq )
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth )        &
              .or. present( conlwr ) .or. present( conupr ) 
!
    probtest2 = c5_m2
    if ( present(probtest) ) then
        if ( zero>=probtest .or. probtest>=one ) then
            call merror( name_proc//tseries_error59 )
        else
            probtest2 = probtest
        end if
    end if
!
    if ( present(trend) ) then
        trendb = trend
    else
        trendb = 1_i4b
    end if
!
    if ( present(trend2) ) then
        trend2b = trend2
    else
        trend2b = 0_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
    else
        win2 = 3_i4b
    end if
!
    initfft2 = true
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( (/ m, sld2/), dim=2_i4b )
    end if
!
    normpsd2 = true
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    overlap2 = false
    if ( present(overlap) ) then
        overlap2 = overlap
    end if
!
    smooth = false
    if ( present(smooth_param) ) then
!
        nparam = size( smooth_param )
!
        if ( nparam<=0_i4b  )    &
        call merror( name_proc//tseries_error50 )
!
!       CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIEL FILTERS.
!
        if ( any( smooth_param(:nparam)<=0_i4b .or. smooth_param(:nparam)>=nf ) )   &
        call merror( name_proc//tseries_error56 )
!
        smooth = true
!
    end if
!
!   FIND THE FIRST INTEGER GREATER THAN OR EQUAL TO n THAT IS EVENLY
!   DIVISIBLE BY l.
!
    i = mod( n, l )
!
    if ( i==0_i4b ) then
        n2 = n
    else
        n2 = n + l - i
    end if
!
!   DETERMINE THE NUMBER OF SEGMENTS.
!
    if ( overlap2 ) then
        nseg = (2*n2/l) - 1_i4b
        step = ld2
    else
        nseg = n2/l
        step = l
    end if
!
!   CALCULATE DATA WINDOW IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:l)  = data_window( l, win2, taperp=taperp )
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:l), wk(:l) )
!
    else
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:l)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( l, stnd )
!
    end if
!
!   COMPUTE DEGREES OF FREEDOM IF REQUIRED.
!
    if ( out_dof ) then
!
        if ( nseg/=1_i4b .and. overlap2 .and. win2/=1_i4b .and. win2/=3_i4b  )    &
        call merror( name_proc//tseries_error58 )
!
        edof2 = estim_dof( wk(:l), win=win, smooth_param=smooth_param(:nparam), &
                           l0=l0, nseg=nseg, overlap=overlap )
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trendb>=1_i4b .and. trendb<=3_i4b ) then
        call detrend_rm( mat(:m,:n), trendb )
    end if
!
!   ZERO OUTPUT PSD ESTIMATES FOR SUMMING OVER THE m SEGMENTS.
!
    psmat(:m,:nf) = zero
!
!   COMPUTE PSD ESTIMATES BY SUMMING OVER THE m SEGMENTS.
!
#ifdef _OPENMP
    i1 = omp_get_num_procs()
    i2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               i1>1_i4b                     .and.      &
               i2>1_i4b                     .and.      &
               nseg>=i2
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(k,i,i1,i2,l2,iok,psmatb,seg,cwk)
!
!       ALLOCATE WORK MATRIX.
!
        allocate( psmatb(m,nf), seg(m,sl), cwk(m,nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ZERO PARTIAL OUTPUT PSD ESTIMATES FOR SUMMING OVER EACH THREAD.
!
        psmatb(:m,:nf) = zero
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(:m,l+1_i4b:sl) = zero
        end if    
!
!$OMP DO SCHEDULE(STATIC)
!
        do k=1_i4b, nseg
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg(:m,:l2) = mat(:m,i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENT IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rm( seg(:m,:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:m,:l2) = seg(:m,:l2)*spread( wk(:l2), dim=1, ncopies=m )
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(:m,l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT OF THE TIME SEGMENT.
!
            call real_fft( seg(:m,:sl), cwk(:m,:nf), true )
!
!           UPDATE POWER SPECTRUM ESIMATES.
!
            psmatb(:m,:sld2) = psmatb(:m,:sld2) + real(cwk(:m,:sld2),stnd)**2 + aimag(cwk(:m,:sld2))**2
            psmatb(:m,nf)    = psmatb(:m,nf)    + real(cwk(:m,nf),stnd)**2
!
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (updatepsd2)
        psmat(:m,:nf) = psmat(:m,:nf) + psmatb(:m,:nf)
!$OMP END CRITICAL (updatepsd2)
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( psmatb, seg, cwk )
!
!$OMP END PARALLEL
!
    else
!
#endif
!
!       ALLOCATE WORK MATRICES.
!
        allocate( seg(m,sl), cwk(m,nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(:m,l+1_i4b:sl) = zero
        end if    
!
        do k=1_i4b, nseg
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg(:m,:l2) = mat(:m,i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENT IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rm( seg(:m,:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:m,:l2) = seg(:m,:l2)*spread( wk(:l2), dim=1, ncopies=m )
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(:m,l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT FOR THE SERIES.
!
            call real_fft( seg(:m,:sl), cwk(:m,:nf), true )
!
!           UPDATE POWER SPECTRUM ESIMATES.
!
            psmat(:m,:sld2) = psmat(:m,:sld2) + real(cwk(:m,:sld2),stnd)**2 + aimag(cwk(:m,:sld2))**2
            psmat(:m,nf)    = psmat(:m,nf)    + real(cwk(:m,nf),stnd)**2
!
        end do
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( seg, cwk )
!
#ifdef _OPENMP
!
    end if
!
#endif
!
!   NORMALIZE THE POWER SPECTRA ESTIMATES.
!
    c1 = one/( sumw*real(l*nseg, stnd) )
!
    if ( normpsd2 ) then
!
        c2 = (real( l, stnd)/pi)*c1
        i1 = 1_i4b
        i2 = nf
!
    else
        c2 = two*c1
        i1 = 2_i4b
        i2 = sld2
!
        psmat(:m,1_i4b) = c1*psmat(:m,1_i4b)
        psmat(:m,nf)    = c1*psmat(:m,nf)
!
    end if
!
    psmat(:m,i1:i2) = c2*psmat(:m,i1:i2)
!
!   SMOOTH THE POWER SPECTRA ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
        call moddan_filter_rm( psmat(:m,:nf), smooth_param(:nparam), sym=one )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
        c1 = one/real( sl, stnd )
        freq(:nf) = arth( zero, c1, nf )
    end if
!
!   OUTPUT DEGREES OF FREEDOM, CONFIDENCE LIMIT FACTORS
!   OR BANDWIDTH IF REQUIRED.
!
    if ( out_dof ) then
!
        if ( present(edof) ) then
            edof = edof2
        end if
!
        if ( present(bandwidth) ) then
            bandwidth = edof2/( two*real( n, stnd ) )
        end if
!
        c1 = probtest2*half
        c2 = one - c1
!
        if ( present(conlwr) ) then
            if ( edof2>=half ) then
                conlwr = edof2/pinvq2( c2, edof2 )
            else
                conlwr = nan()
            end if
        end if
!
        if ( present(conupr) ) then
            if ( edof2>=half ) then
                conupr = edof2/pinvq2( c1, edof2 )
            else
                conupr = nan()
            end if
        end if
!
    end if
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE power_spectrum2_rm
! ____________________________________
!
    end subroutine power_spectrum2_rm
!
! =========================================================================================
!
    subroutine cross_spectrum2_rv( vec, vec2, l, psvec, psvec2, phase, coher, freq, edof, bandwidth,    &
                                   conlwr, conupr, testcoher, ampli, co_spect, quad_spect, prob_coher,  &
                                   initfft, overlap, normpsd, smooth_param, trend, trend2, win, taperp, &
                                   l0, probtest )
!
! Purpose
! _______
!
!   Subroutine CROSS_SPECTRUM2 computes Fast Fourier Transform (FFT) estimates
!   of the power and cross spectra of two real time series.
!
!   The Power Spectral Density (PSD) and Cross Spectral Density (CSD) estimates are returned
!   in units which are the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the first real time series for which the power and cross spectra
!                 must be estimated.
!                 If TREND=1, 2 or 3,  VEC is used as workspace and is transformed.
!
!                 Size(VEC) must be greater or equal to 4.
!
!   VEC2          (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the second real time series for which the power and cross spectra
!                 must be estimated.
!                 If TREND=1, 2 or 3,  VEC2 is used as workspace and is transformed.
!
!                 VEC2 must verify:  size(VEC2) = size(VEC).
!
!   L             (INPUT) integer(i4b)
!                 On entry, an integer used to segment the time series. L is the
!                 length of the segments. L must be a positive even integer, less
!                 or equal to size(VEC), but greater or equal to 4.
!
!                 Spectral computations are at (L/2)+1 frequencies if the optional
!                 argument L0 is absent and are at ((L+L0)/2)+1 frequencies if L0 is
!                 present (L0 is the number of zeros added to each segment).
!
!                 Suggested values for L+L0 are 16, 32, 64 or 128 (e.g. an integer power of two,
!                 in order to speed the computations).
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = ((L+L0)/2) + 1 .
!
!   PSVEC2        (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC2.  
!
!                 PSVEC2 must verify:  size(PSVEC2) = ((L+L0)/2) + 1 .
!
!   PHASE         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the phase of the cross spectrum, given in fractions
!                 of a circle (e.g. on the closed interval (0,1) ).
!
!                 PHASE must verify:  size(PHASE) = ((L+L0)/2) + 1 .
!
!   COHER         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the squared coherency  estimates for all frequencies.
!
!                 COHER must verify:  size(COHER) = ((L+L0)/2) + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/(L+L0)
!                 for i=1,2, ... , ((L+L0)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = (L+L0)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 and cross spectrum estimates.
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the bandwidth of the power and cross spectrum estimates.
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the PSVEC(:) and 
!                 PSVEC2(:) arguments ) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!   
!   TESTCOHER     (OUTPUT, OPTIONAL) real(stnd)
!                 On output, this argument specifies the critical value for testing the null
!                 hypothesis that the squared coherency is zero at the PROBTEST * 100% significance
!                 level (e.g. elements of COHER(:) less than TESTCOHER should be regarded as not
!                 significantly different from zero at the PROBTEST * 100% significance level).
!
!   AMPLI         (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the cross-amplitude spectrum.  
!
!                 AMPLI must verify:  size(AMPLI) = ((L+L0)/2) + 1 .
!
!   CO_SPECT      (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the co-spectrum (e.g. the real part of cross-spectrum).  
!
!                 CO_SPECT must verify:  size(CO_SPECT) = ((L+L0)/2) + 1 .
!
!   QUAD_SPECT    (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the quadrature spectrum (e.g. the imaginary part of
!                 cross-spectrum with a minus sign).  
!
!                 QUAD_SPECT must verify:  size(QUAD_SPECT) = ((L+L0)/2) + 1 .
!
!   PROB_COHER    (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the probabilities that the computed sample squared coherencies
!                 came from an ergodic stationary bivariate process with (corresponding)
!                 squared coherencies equal to zero.
!
!                 PROB_COHER  must verify:  size(PROB_COHER) = ((L+L0)/2)+1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if INITFFT is set to false, it is assumed that a call to subroutine
!                 INIT_FFT has been done before calling subroutine CROSS_SPECTRUM2 in order to 
!                 sets up constants and functions for use by subroutine FFT which is called inside
!                 subroutine CROSS_SPECTRUM2 (the call to INITFFT must have the following form: 
!
!                      call init_fft( (L+L0)/2 )
!
!                 If INITFFT is set to true, the call to INIT_FFT is done inside subroutine
!                 CROSS_SPECTRUM2 and a call to END_FFT is also done before leaving
!                 subroutine CROSS_SPECTRUM2.
!
!                 The default is INITFFT=true .
!
!   OVERLAP       (INPUT, OPTIONAL) logical(lgl)
!                 If:
!
!                 - OVERLAP = false, the subroutine segments the data 
!                   without any overlapping.
!                 - OVERLAP = true, the subroutine overlaps the segments
!                   by one half of their length (which is equal to L).
!
!                 In both cases, zeros are eventually added to each segment (if argument L0 is present)
!                 and each segment will be FFT'd, and the resulting periodograms
!                 will averaged together to obtain a Power Spectrum Density estimate at the
!                 ((L+L0)/2)+1 frequencies.
!
!                 The default is OVERLAP=false .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if NORMPSD is set to true, the power and cross spectra estimates are normalized
!                 in such a way that the total area under the power spectrum is equal to the variance
!                 of the time series VEC and VEC2.
!                 If NORMPSD is set to false, the sum of the PSD estimates
!                 (e.g. sum(PSVEC(2:)) and sum(PSVEC2(2:)) ) is equal to the variance of the corresponding
!                 time series.
!
!                 The default is NORMPSD=true .
!
!   SMOOTH_PARAM  (INPUT, OPTIONAL) integer(i4b), dimension(:)
!                 If SMOOTH_PARAM is used, the power and cross spectra estimates are computed by repeated
!                 smoothing of the periodograms and cross-periodogram with modified Daniell weights.
!
!                 On entry, SMOOTH_PARAM(:) gives the array of the half-lengths of the
!                 modified Daniell filters to be applied.
!
!                 All the values in SMOOTH_PARAM(:) must be greater than 0 and less
!                 than (L+L0)/2 + 1 .
!
!                 Size(SMOOTH_PARAM) must be greater or equal to 1.
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The mean of the two time series is removed before computing the spectra
!                 - TREND=+2 The drift from the two time series is removed before computing the spectra
!                 - TREND=+3 The least-squares line from the two time series is removed before
!                   computing the spectra.
!
!                 For other values of TREND nothing is done before estimating the power and cross spectra.
!
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   TREND2        (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND2=+1 The mean of the time segment is removed before computing the cross-spectrum
!                   on this segment.
!                 - TREND2=+2 The drift from the time segment is removed before computing the cross-spectrum
!                   on this segment.
!                 - TREND2=+3 The least-squares line from the time segment is removed before
!                   computing the cross-spectrum on this segment.
!
!                 For other values of TREND2 nothing is done before estimating the cross-spectrum
!                 on each segment.
!
!                 The default is TREND2=0, e.g. nothing is done before estimating the power spectrum
!                 on each segment.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power and cross spectra. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!
!   L0            (INPUT, OPTIONAL) integer(i4b)
!                 The number of zeros added to each time segment in order to obtain more finely
!                 spaced spectral estimates. L+L0 must be a positive even integer.
!
!                 The default is L0=0, e.g. no zeros are added to each time segment.
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ) and the critical value for
!                 testing the null hypothesis that the squared coherency is zero (e.g.
!                 the TESTCOHER optional argument).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the two time series (e.g. TREND=1,2,3), the series
!   are padded with zero on the right such that the length of the resulting two time series is 
!   evenly divisible by L (a positive even integer). The length, N, of these resulting time
!   series is the first integer greater than or equal to size(VEC) which is evenly divisible
!   by L. If size(VEC) is not evenly divisible by L, N is equal to size(VEC)+L-mod(size(VEC),L).
!   
!   Optionally, the mean or the trend may also be removed from each time segment (e.g. TREND2=1,2,3).
!   Optionally, zeros may be added to each time segment (e.g. the optional arguemnt L0) if more
!   finely spaced spectral esimates are desired.
!   
!   The stability of the power and cross spectra  estimates depends on the averaging process.
!   That is, the greater the number of segments ( N/L if OVERLAP=false and (2N/L)-1 if OVERLAP=true),
!   the more stable the resulting power and cross spectra estimates.
!   
!   Optionally, these power and cross spectra estimates may then be smoothed again
!   in the frequency domain by modified Daniell filters (e.g. if argument SMOOTH_PARAM is used).
!
!   The computed equivalent number of degrees of freedom and bandwidth must be divided by two
!   for the zero and Nyquist frequencies.
!
!   Furthermore, the computed equivalent number of degrees of freedom, bandwidth, lower and
!   upper (1-PROBTEST) * 100% confidence limit factors and critical value for the squared coherency
!   (e.g. arguments EDOF, BANDWIDTH, CONLWR, CONUPR and TESTCOHER) are not right near the zero
!   and Nyquist frequencies if the PSD estimates have been smoothed by modified Daniell filters.
!   The reason is that CROSS_SPECTRUM2 assumes that smoothing involves averaging independent frequency
!   ordinates. This is true except near the zero and Nyquist frequencies where an average
!   may contain contributions from negative frequencies, which are identical to and hence not
!   independent of positive frequency spectral values. Thus, the number of degrees of freedom
!   in PSD estimates near the 0 and Nyquist frequencies are as little as half the number
!   of degrees of freedom of the spectral estimates away from these frequency extremes if
!   the optional argument SMOOTH_PARAM is used.
!      
!   If the optional argument SMOOTH_PARAM is used, the computed equivalent number of degrees
!   of freedom, bandwidth, lower and upper (1-PROBTEST) * 100% confidence limit factors and
!   critical value for the squared coherency are right for PSD estimates  at frequencies
!      
!        (i-1)/(L+L0)  for i= (nparam+1)/2 + 1   to ( (L+L0) - nparam + 1)/2
!      
!   where nparam = 2 * (2+sum(SMOOTH_PARAM(:)))- 1,
!   (e.g. for frequencies i/(L+L0) for i = (nparam+1)/2, ... , ((L+L0)-nparam-1)/2 ) .
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, assert_eq, arth
    use Reals_Constants,   only : zero, half, one, two, four, pi, twopi, c5_m2
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error15, tseries_error17,    &
                                  tseries_error18, tseries_error24, tseries_error25, tseries_error50,   &
                                  tseries_error56, tseries_error58, tseries_error59
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : probf2, pinvf2, pinvq2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:) :: vec, vec2
    real(stnd), intent(out),    dimension(:) :: psvec, psvec2, phase, coher
!
    real(stnd), intent(in),                optional :: taperp, probtest
    real(stnd), intent(out),               optional :: edof, bandwidth, conlwr, conupr, testcoher
    real(stnd), intent(out), dimension(:), optional :: freq, ampli, co_spect, quad_spect, prob_coher
!
    integer(i4b),  intent(in)                         :: l
    integer(i4b),  intent(in),               optional :: trend, trend2, win, l0
    integer(i4b),  intent(in), dimension(:), optional :: smooth_param
!
    logical(lgl),  intent(in), optional :: normpsd, initfft, overlap
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, n2, nf, k, l2, ld2, m, win2, trendb, trend2b,      &
                    i, i1, i2, step, sl, sld2, nparam
    integer      :: iok
!
    real(stnd)                            :: c1, c2, sumw, edof2, probtest2, con
    real(stnd), dimension(l)              :: wk
    real(stnd), dimension(:), allocatable :: seg1, seg2
!
    complex(stnd), dimension(:), allocatable  :: cwk1, cwk2
!
    logical(lgl)  :: normpsd2, initfft2, overlap2, smooth, out_freq, out_ampli,  &
                     out_co_spect, out_quad_spect, out_prob_coher, out_dof
!
#ifdef _OPENMP
!
    real(stnd), dimension(:), allocatable :: psvecb, psvec2b, rcospect, icospect
!
    logical  :: test_par
!
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='cross_spectrum2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n =  assert_eq( int(size(vec),i4b) ,          &
                    int(size(vec2),i4b) ,         &
                    name_proc )
!
    if ( l<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    if ( l>n )     &
    call merror( name_proc//tseries_error17 )
!
!   CHECK IF l IS A POSITIVE EVEN INTEGER.
!
    ld2 = l/2_i4b
!
    if ( l/=2_i4b*ld2 )     &
    call merror( name_proc//tseries_error18 )
!
!   DETERMINE THE LENGTH OF THE SEGMENTS.
!
    if ( present(l0) ) then
!
        if ( l0<0_i4b )     &
        call merror( name_proc//tseries_error24 )
!
        sl = l + l0
!
!       CHECK IF sl IS A POSITIVE EVEN INTEGER.
!
        sld2 = sl/2_i4b
!
        if ( sl/=2_i4b*sld2 )     &
        call merror( name_proc//tseries_error25 )
!
    else
!
        sl   = l
        sld2 = ld2
!
    end if
!
    nf  = sld2 + 1_i4b
!
    call assert( logical(int(size(psvec),i4b)==nf,lgl),   &
                 logical(int(size(psvec2),i4b)==nf,lgl),  &
                 logical(int(size(phase),i4b)==nf,lgl),   &
                 logical(int(size(coher),i4b)==nf,lgl),   &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_ampli = present( ampli )
!
    if ( out_ampli ) then
        call assert( logical(int(size(ampli),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_co_spect = present( co_spect )
!
    if ( out_co_spect ) then
        call assert( logical(int(size(co_spect),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_quad_spect = present( quad_spect )
!
    if ( out_quad_spect ) then
        call assert( logical(int(size(quad_spect),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_prob_coher = present( prob_coher )
!
    if ( out_prob_coher ) then
        call assert( logical(int(size(prob_coher),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth ) .or. out_prob_coher   &
              .or. present( conlwr ) .or. present( conupr ) .or. present( testcoher ) 
!
    probtest2 = c5_m2
!
    if ( present(probtest) ) then
        if ( zero>=probtest .or. probtest>=one ) then
            call merror( name_proc//tseries_error59 )
        else
            probtest2 = probtest
        end if
    end if
!
    if ( present(trend) ) then
        trendb = trend
    else
        trendb = 1_i4b
    end if
!
    if ( present(trend2) ) then
        trend2b = trend2
    else
        trend2b = 0_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
!
    else
!
        win2 = 3_i4b
!
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( sld2 )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    overlap2 = false
!
    if ( present(overlap) ) then
        overlap2 = overlap
    end if
!
    smooth = false
!
    if ( present(smooth_param) ) then
!
        nparam = size( smooth_param )
!
        if ( nparam<=0_i4b  )    &
        call merror( name_proc//tseries_error50 )
!
!       CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIEL FILTERS.
!
        if ( any( smooth_param(:nparam)<=0_i4b .or. smooth_param(:nparam)>=nf ) )   &
        call merror( name_proc//tseries_error56 )
!
        smooth = true
!
    end if
!
!   FIND THE FIRST INTEGER GREATER THAN OR EQUAL TO n THAT IS EVENLY
!   DIVISIBLE BY l.
!
    i = mod( n, l )
!
    if ( i==0_i4b ) then
        n2 = n
    else
        n2 = n + l - i
    end if
!
!   DETERMINE THE NUMBER OF SEGMENTS.
!
    if ( overlap2 ) then
        m    = (2*n2/l) - 1_i4b
        step = ld2
    else
        m    = n2/l
        step = l
    end if
!
!   CALCULATE DATA WINDOW IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:l)  = data_window( l, win2, taperp=taperp )
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:l), wk(:l) )
!
    else
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:l)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( l, stnd )
!
    end if
!
!   COMPUTE DEGREES OF FREEDOM IF REQUIRED.
!
    if ( out_dof ) then
!
        if ( m/=1_i4b .and. overlap2 .and. win2/=1_i4b .and. win2/=3_i4b  )    &
        call merror( name_proc//tseries_error58 )
!
        edof2 = estim_dof( wk(:l), win=win, smooth_param=smooth_param(:nparam), &
                           l0=l0, nseg=m, overlap=overlap )
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trendb>=1_i4b .and. trendb<=3_i4b ) then
!
        call detrend_rv( vec(:n), trendb )
        call detrend_rv( vec2(:n), trendb )
!
    end if
!
!   ZERO OUTPUT POWER AND CROSS SPECTRA ESTIMATES FOR SUMMING OVER THE m SEGMENTS.
!
    psvec(:nf)  = zero
    psvec2(:nf) = zero
!
!   USE THE ARRAYS coher AND phase TO STORE THE REAL AND IMAGINARY PARTS OF THE
!   CROSS SPECTRUM, RESPECTIVELY.
!
    coher(:nf)  = zero
    phase(:nf)  = zero
!
!   COMPUTE PSD AND CSD ESTIMATES BY SUMMING OVER THE m SEGMENTS.
!
#ifdef _OPENMP
    i1 = omp_get_num_procs()
    i2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               i1>1_i4b                     .and.      &
               i2>1_i4b                     .and.      &
               m>=i2
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(k,i,i1,i2,l2,iok,psvecb,psvec2b,seg1,seg2,cwk1,cwk2,rcospect,icospect)
!
!       ALLOCATE WORK VECTORS.
!
        allocate( psvecb(nf), psvec2b(nf), rcospect(nf), icospect(nf),     &
                  seg1(sl), seg2(sl), cwk1(nf), cwk2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ZERO PARTIAL OUTPUT PSD AND CSD ESTIMATES FOR SUMMING OVER EACH THREAD.
!
        psvecb(:nf)  = zero
        psvec2b(:nf) = zero
!
        rcospect(:nf) = zero
        icospect(:nf) = zero
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg1(l+1_i4b:sl) = zero
            seg2(l+1_i4b:sl) = zero
        end if    
!
!$OMP DO SCHEDULE(STATIC)
!
        do k=1_i4b, m
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg1(:l2) = vec(i1:i2)
            seg2(:l2) = vec2(i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENTS IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rv( seg1(:l2), trend2b )
                call detrend_rv( seg2(:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW TO THE SERIES.
!
            if ( win2/=2_i4b ) then
                seg1(:l2) = seg1(:l2)*wk(:l2)
                seg2(:l2) = seg2(:l2)*wk(:l2)
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg1(l2+1_i4b:l) = zero
                seg2(l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT FOR THE FIRST SERIES.
!
            call real_fft( seg1(:sl), cwk1(:nf), true)
!
!           COMPUTE FFT FOR THE SECOND SERIES.
!
            call real_fft( seg2(:sl), cwk2(:nf), true)
!
!           UPDATE SPECTRAL ESTIMATES.
!
            psvecb(:sld2) = psvecb(:sld2) + real(cwk1(:sld2),stnd)**2 + aimag(cwk1(:sld2))**2
            psvecb(nf)    = psvecb(nf)    + real(cwk1(nf),stnd)**2
!
            psvec2b(:sld2) = psvec2b(:sld2) + real(cwk2(:sld2),stnd)**2 + aimag(cwk2(:sld2))**2
            psvec2b(nf)    = psvec2b(nf)    + real(cwk2(nf),stnd)**2
!
!           UPDATE CO-SPECTRUM ESTIMATES.
!
            rcospect(:sld2)  = rcospect(:sld2)  + real(cwk1(:sld2),stnd)*real(cwk2(:sld2),stnd) +   &
                                                  aimag(cwk1(:sld2))*aimag(cwk2(:sld2))
            rcospect(nf)     = rcospect(nf)     + real(cwk1(nf),stnd)*real(cwk2(nf),stnd)
!
!           UPDATE QUADRATURE-SPECTRUM ESTIMATES.
!
            icospect(:sld2)  = icospect(:sld2)  + aimag(cwk1(:sld2))*real(cwk2(:sld2),stnd)     -   &
                                                  real(cwk1(:sld2),stnd)*aimag(cwk2(:sld2))
!
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (updatepscd)
!
        psvec(:nf)    = psvec(:nf) + psvecb(:nf)
        psvec2(:nf)   = psvec2(:nf) + psvec2b(:nf)
!
        coher(:nf)    = coher(:nf)   + rcospect(:nf) 
        phase(:sld2)  = phase(:sld2) + icospect(:sld2)
!
!$OMP END CRITICAL (updatepscd)
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( psvecb, psvec2b, seg1, seg2, cwk1, cwk2, rcospect, icospect )
!
!$OMP END PARALLEL
!
!       REALLOCATE WORK VECTORS.
!
        allocate( seg1(nf), seg2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
    else
!
#endif
!
!       ALLOCATE WORK VECTORS.
!
        allocate( seg1(sl), seg2(sl), cwk1(nf), cwk2(nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg1(l+1_i4b:sl) = zero
            seg2(l+1_i4b:sl) = zero
        end if    
!
        do k=1_i4b, m
!
            i = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg1(:l2) = vec(i1:i2)
            seg2(:l2) = vec2(i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENTS IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rv( seg1(:l2), trend2b )
                call detrend_rv( seg2(:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW TO THE SERIES.
!
            if ( win2/=2_i4b ) then
                seg1(:l2) = seg1(:l2)*wk(:l2)
                seg2(:l2) = seg2(:l2)*wk(:l2)
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg1(l2+1_i4b:l) = zero
                seg2(l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT FOR THE FIRST SERIES.
!
            call real_fft( seg1(:sl), cwk1(:nf), true)
!
!           COMPUTE FFT FOR THE SECOND SERIES.
!
            call real_fft( seg2(:sl), cwk2(:nf), true)
!
!           UPDATE POWER SPECTRUM ESIMATES.
!
            psvec(:sld2) = psvec(:sld2) + real(cwk1(:sld2),stnd)**2 + aimag(cwk1(:sld2))**2
            psvec(nf)    = psvec(nf)    + real(cwk1(nf),stnd)**2
!
            psvec2(:sld2) = psvec2(:sld2) + real(cwk2(:sld2),stnd)**2 + aimag(cwk2(:sld2))**2
            psvec2(nf)    = psvec2(nf)    + real(cwk2(nf),stnd)**2
!
!           UPDATE CO-SPECTRUM ESIMATES.
!
            coher(:sld2)  = coher(:sld2)  + real(cwk1(:sld2),stnd)*real(cwk2(:sld2),stnd) +   &
                                            aimag(cwk1(:sld2))*aimag(cwk2(:sld2))
            coher(nf)    = coher(nf)    + real(cwk1(nf),stnd)*real(cwk2(nf),stnd)
!
!           UPDATE QUADRATURE-SPECTRUM ESIMATES.
!
            phase(:sld2)  = phase(:sld2)  + aimag(cwk1(:sld2))*real(cwk2(:sld2),stnd)    -   &
                                            real(cwk1(:sld2),stnd)*aimag(cwk2(:sld2))
!
        end do
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( cwk1, cwk2 )
!        
#ifdef _OPENMP
!
    end if
!
#endif
!
!   NORMALIZE THE POWER AND CROSS SPECTRA ESTIMATES.
!
    c1 = one/( sumw*real(l*m, stnd) )
!
    if ( normpsd2 ) then
!
        c2 = (real( l, stnd)/pi)*c1
        i1 = 1_i4b
        i2 = nf
!
    else
!
        c2 = two*c1
        i1 = 2_i4b
        i2 = sld2
!
        psvec(1_i4b)  = c1*psvec(1_i4b)
        psvec2(1_i4b) = c1*psvec2(1_i4b)
!
        seg1(1_i4b) = c1*coher(1_i4b)
        seg2(1_i4b) = c1*phase(1_i4b)
!
        psvec(nf)  = c1*psvec(nf)
        psvec2(nf) = c1*psvec2(nf)
!
        seg1(nf)  = c1*coher(nf)
        seg2(nf)  = c1*phase(nf)
!
    end if
!
    psvec(i1:i2)  = c2*psvec(i1:i2)
    psvec2(i1:i2) = c2*psvec2(i1:i2)
!
    seg1(i1:i2) = c2*coher(i1:i2)
    seg2(i1:i2) = c2*phase(i1:i2)
!
!   SMOOTH THE POWER AND CROSS SPECTRA ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
!       SMOOTH POWER SPECTRUM OF THE FIRST SERIES.
!
        call moddan_filter_rv( psvec(:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH POWER SPECTRUM OF THE SECOND SERIES.
!
        call moddan_filter_rv( psvec2(:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH CO-SPECTRUM.
!
        call moddan_filter_rv( seg1(:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH QUADRATURE-SPECTRUM.
!
        call moddan_filter_rv( seg2(:nf), smooth_param(:nparam), sym=-one )
!
    end if
!
    if ( out_co_spect ) then
!
!       OUTPUT THE CO-SPECTRUM IF REQUIRED.
!
        co_spect(:nf) = seg1(:nf)
!
    end if
!
    if ( out_quad_spect ) then
!
!       OUTPUT THE QUADRATURE-SPECTRUM IF REQUIRED.
!
        quad_spect(:nf) = -seg2(:nf)
!
    end if
!
!   COMPUTE PHASE OF THE CROSS-SPECTRUM.
!
    where ( seg1(:nf)==zero .and. seg2(:nf)==zero )
        phase(:nf) = zero
    elsewhere
        coher(:nf)   = atan2( seg2(:nf), seg1(:nf) )
        phase(:nf) = (one/twopi)*coher(:nf) + merge( zero, one, coher(:nf)>=zero ) 
    end where
!
!   COMPUTE MAGNITUDE OF THE CROSS-SPECTRUM.
!
    coher(:nf) = seg1(:nf)**2 + seg2(:nf)**2
    seg2(:nf)  = coher(:nf)
!
!   COMPUTE THE SQUARED COHERENCE FOR ALL FREQUENCIES.
!
    where ( psvec(:nf)/=zero .and. psvec2(:nf)/=zero )
        seg1(:nf)  = psvec(:nf)*psvec2(:nf)
        coher(:nf) = min( seg2(:nf)/seg1(:nf), one )
    elsewhere
        coher(:nf) = zero
    end where 
!
    if ( out_ampli ) then
!
!       COMPUTE AND OUTPUT AMPLITUDE OF THE CROSS-SPECTRUM IF REQUIRED.
!
        ampli(:nf) = sqrt( seg2(:nf) )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
!
        c1 = one/real( sl, stnd )
        freq(:nf) = arth(zero, c1, nf)
!
    end if
!
!   OUTPUT DEGREES OF FREEDOM, BANDWIDTH, CONFIDENCE LIMIT FACTORS 
!   OR PROBABILITIES FOR SQUARED COHERENCIES IF REQUIRED.
!
    if ( out_dof ) then
!
        if ( present(edof) ) then
            edof = edof2
        end if
!
        if ( present(bandwidth) ) then
            bandwidth = edof2/( two*real( n, stnd ) )
        end if
!
        c1 = probtest2*half
        c2 = one - c1
!
        sumw = nan()
!
        if ( present(conlwr) ) then
            if ( edof2>=half ) then
                conlwr = edof2/pinvq2( c2, edof2 )
            else
                conlwr = sumw
            end if
        end if
!
        if ( present(conupr) ) then
            if ( edof2>=half ) then
                conupr = edof2/pinvq2( c1, edof2 )
            else
                conupr = sumw
            end if
        end if
!
        if ( present(testcoher) ) then
            if ( edof2>two ) then
!
!                c1 = two/(edof2-two)
!                c2 = one - probtest2**c1
!                testcoher = max( min( one, c2 ), zero )
!
                 con = one - probtest2
                 c2  = edof2 - two
                 c1  = two*pinvf2( con, two, c2 )
                 testcoher = c1/( c2 + c1 )
            else
                testcoher = sumw
            end if
        end if
!
!       COMPUTE AND OUTPUT SIGNIFICANCE PROBABILITY FOR SQUARED COHERENCIES.
!
        if ( out_prob_coher ) then
!
            c2 = edof2 - two    
!
            if ( c2>zero ) then
!
                where( coher(2_i4b:sld2)/=one )
                    seg1(2_i4b:sld2) = ( (edof2/two - one)*coher(2_i4b:sld2) )/  &
                                      ( one - coher(2_i4b:sld2) )
                elsewhere
                    seg1(2_i4b:sld2) = zero
                end where 
!            
                prob_coher(2_i4b:sld2) = probf2( seg1(2_i4b:sld2), two, c2, true )
!            
                where( coher(2_i4b:sld2)==one ) prob_coher(2_i4b:sld2) = zero
!
            else
!
                prob_coher(2_i4b:sld2) = sumw
!
            end if
!            
            c2 = edof2/two - two    
!
            if ( c2>zero ) then
!            
                do i1 = 1_i4b, nf, sld2
!            
                    c1 = coher(i1)
!            
                    if ( c1/=one ) then
!            
                        con = ( (edof2/four - one)*c1 )/( one - c1 )
!            
                        prob_coher(i1) = probf2( con, two, c2, true )
!            
                    else
                        prob_coher(i1) = zero
                    end if 
!
                end do
!
            else
!
                prob_coher(1_i4b) = sumw
                prob_coher(nf)    = sumw
!
            end if
!
        end if
!
    end if
!
!   DEALLOCATE WORK ARRAYS.
!
    deallocate( seg1, seg2 )
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE cross_spectrum2_rv
! ____________________________________
!
    end subroutine cross_spectrum2_rv
!
! =========================================================================================
!
    subroutine cross_spectrum2_rm( vec, mat, l, psvec, psmat, phase, coher, freq, edof, bandwidth, &
                                   conlwr, conupr, testcoher, ampli, co_spect, quad_spect,         &
                                   prob_coher, initfft, overlap, normpsd, smooth_param, trend,     &
                                   trend2, win, taperp, l0, probtest )
!
! Purpose
! _______
!
!   Subroutine CROSS_SPECTRUM2 computes Fast Fourier Transform (FFT) estimates
!   of the power and cross spectra of the real time series, VEC, and the multi-channel
!   real time series MAT.
!
!   The Power Spectral Density (PSD) and Cross Spectral Density (CSD) estimates are returned
!   in units which are the square of the data (if NORMPSD=false) or in spectral density units (if
!   NORMPSD=true).
!
!
! Arguments
! _________
!
!   VEC           (INPUT/OUTPUT) real(stnd), dimension(:)
!                 On entry, the real time series for which the power and cross spectra
!                 must be estimated.
!                 If TREND=1, 2 or 3,  VEC is used as workspace and is transformed.
!
!                 Size(VEC) must be greater or equal to 4.
!
!   MAT           (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                 On entry, the multi-channel real time series for which the power 
!                 and cross spectra must be estimated. Each row of MAT is a real time series.
!                 If TREND=1, 2 or 3,  MAT is used as workspace and is transformed.
!
!                 The shape of MAT must verify:  size(MAT,2) = size(VEC).
!
!   L             (INPUT) integer(i4b)
!                 On entry, an integer used to segment the time series. L is the
!                 length of the segments. L must be a positive even integer, less
!                 or equal to size(VEC), but greater or equal to 4.
!
!                 Spectral computations are at (L/2)+1 frequencies if the optional
!                 argument L0 is absent and are at ((L+L0)/2)+1 frequencies if L0 is
!                 present (L0 is the number of zeros added to each segment).
!
!                 Suggested values for L+L0 are 16, 32, 64 or 128 (e.g. an integer power of two,
!                 in order to speed the computations).
!
!   PSVEC         (OUTPUT) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the Power Spectral Density (PSD) estimates of VEC.  
!
!                 PSVEC must verify:  size(PSVEC) = ((L+L0)/2) + 1 .
!
!   PSMAT         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the Power Spectral Density (PSD) estimates of each row of MAT.
!
!                 The shape of PSMAT must verify:
!
!                 - size(PSMAT,1) =  size(MAT,1) ;
!                 - size(PSMAT,2) = ((L+L0)/2) + 1 .
!
!   PHASE         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the phase of the cross spectrum, given in fractions
!                 of a circle (e.g. on the closed interval (0,1) ).
!
!                 The shape of PHASE must verify:
!
!                 - size(PHASE,1) =  size(MAT,1) ;
!                 - size(PHASE,2) = ((L+L0)/2) + 1 .
!
!   COHER         (OUTPUT) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the squared coherency  estimates for all frequencies.
!
!                 The shape of COHER must verify:
!
!                 - size(COHER,1) =  size(MAT,1) ;
!                 - size(COHER,2) = ((L+L0)/2) + 1 .
!
!   FREQ          (OUTPUT, OPTIONAL) real(stnd), dimension(:)
!                 On exit, a real vector of length ((L+L0)/2)+1 containing
!                 the frequencies at which the spectral quantities are calculated
!                 in cycles per unit of time.
!
!                 The spectral estimates are taken at frequencies (i-1)/(L+L0)
!                 for i=1,2, ... , ((L+L0)/2 + 1).
!
!                 FREQ must verify:  size(FREQ) = (L+L0)/2 + 1 .
!
!   EDOF          (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the equivalent number of degrees of freedom of the power
!                 and cross spectrum estimates.
!
!   BANDWIDTH     (OUTPUT, OPTIONAL) real(stnd)
!                 On exit, the bandwidth of the power and cross spectrum estimates.
!   
!   CONLWR        (OUTPUT, OPTIONAL) real(stnd)
!   
!   CONUPR        (OUTPUT, OPTIONAL) real(stnd)
!                 On output, these arguments specify the lower and upper (1-PROBTEST) * 100% confidence
!                 limit factors, respectively. Multiply the PSD estimates (e.g. the PSVEC(:) and 
!                 PSMAT(:,:) arguments ) by these constants to get the lower and upper
!                 limits of a (1-PROBTEST) * 100% confidence interval for the PSD estimates.
!   
!   TESTCOHER     (OUTPUT, OPTIONAL) real(stnd)
!                 On output, this argument specifies the critical value for testing the null
!                 hypothesis that the squared coherency is zero at the PROBTEST * 100% significance
!                 level (e.g. elements of COHER(:,:) less than TESTCOHER should be regarded as not
!                 significantly different from zero at the PROBTEST * 100% significance level).
!
!   AMPLI         (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the cross-amplitude spectra.
!
!                 The shape of AMPLI must verify:
!
!                 - size(AMPLI,1) =  size(MAT,1) ;
!                 - size(AMPLI,2) = ((L+L0)/2) + 1 .
!
!   CO_SPECT      (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the co-spectra (e.g. the real part of cross-spectra).  
!
!
!                 The shape of CO_SPECT must verify:
!
!                 - size(CO_SPECT,1) =  size(MAT,1) ;
!                 - size(CO_SPECT,2) = ((L+L0)/2) + 1 .
!
!   QUAD_SPECT    (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the quadrature spectrum (e.g. the imaginary part of
!                 cross-spectrum with a minus sign).  
!
!                 The shape of QUAD_SPECT must verify:
!
!                 - size(QUAD_SPECT,1) =  size(MAT,1) ;
!                 - size(QUAD_SPECT,2) = ((L+L0)/2) + 1 .
!
!   PROB_COHER    (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!                 On exit, a real matrix with size(MAT,1) rows and ((L+L0)/2) + 1 columns
!                 containing the probabilities that the computed sample squared coherencies
!                 came from an ergodic stationary bivariate process with (corresponding)
!                 squared coherencies equal to zero.
!
!                 The shape of PROB_COHER must verify:
!
!                 - size(PROB_COHER,1) =  size(MAT,1) ;
!                 - size(PROB_COHER,2) =  ((L+L0)/2) + 1 .
!
!   INITFFT       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - INITFFT = false, it is assumed that a call to subroutine
!                   INIT_FFT has been done before calling subroutine CROSS_SPECTRUM2 in order to 
!                   sets up constants and functions for use by subroutine FFT which is called inside
!                   subroutine CROSS_SPECTRUM2. This call to INITFFT must have the following form: 
!
!                      call init_fft( (/ size(MAT,1), (L+L0)/2 /), dim=2_i4b )
!
!                 - INITFFT = true, the call to INIT_FFT is done inside subroutine
!                   CROSS_SPECTRUM2 and a call to END_FFT is also done before leaving
!                   subroutine CROSS_SPECTRUM2.
!
!                 The default is INITFFT=true .
!
!   OVERLAP       (INPUT, OPTIONAL) logical(lgl)
!                 If:
!
!                 - OVERLAP = false, the subroutine segments the data 
!                   without any overlapping.
!                 - OVERLAP = true, the subroutine overlaps the segments
!                   by one half of their length (which is equal to L).
!
!                 In both cases, zeros are eventually added to each segment (if argument L0 is present)
!                 and each segment will be FFT'd, and the resulting periodograms
!                 will averaged together to obtain a Power Spectrum Density estimate at the
!                 ((L+L0)/2)+1 frequencies.
!
!                 The default is OVERLAP=false .
!
!   NORMPSD       (INPUT, OPTIONAL) logical(lgl)
!                 On entry, if:
!
!                 - NORMPSD = true, the power and cross spectra estimates are
!                   normalized in such a way that the total area under the power spectra is equal to
!                   the variance of the time series contained in VEC and in each row of MAT.
!                 - NORMPSD = false, the sum of the PSD estimates
!                   (e.g. sum(PSVEC(2:)) and sum(PSMAT(:,2:),dim=2) ) is equal to the variance of
!                   the corresponding time series.
!
!                 The default is NORMPSD=true .
!
!   SMOOTH_PARAM  (INPUT, OPTIONAL) integer(i4b), dimension(:)
!                 If SMOOTH_PARAM is used, the power and cross spectra estimates are computed by repeated
!                 smoothing of the periodograms and cross-periodogram with modified Daniell weights.
!
!                 On entry, SMOOTH_PARAM(:) gives the array of the half-lengths of the
!                 modified Daniell filters to be applied.
!
!                 All the values in SMOOTH_PARAM(:) must be greater than 0 and less
!                 than ((L+L0)/2)+1 .
!
!   TREND         (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND=+1 The means of the time series are removed before computing the spectra
!                 - TREND=+2 The drifts from time series are removed before computing the spectra
!                 - TREND=+3 The least-squares lines from time series are removed before
!                   computing the spectra.
!
!                 For other values of TREND nothing is done before estimating the power and cross spectra.
!
!                 The default is TREND=1, e.g. the means of the time series are removed before the
!                 computations.
!
!   TREND2        (INPUT, OPTIONAL) integer(i4b)
!                 If:
!
!                 - TREND2=+1 The mean of the time segment is removed before computing the cross-spectrum
!                   on this segment.
!                 - TREND2=+2 The drift from the time segment is removed before computing the cross-spectrum
!                   on this segment.
!                 - TREND2=+3 The least-squares line from the time segment is removed before
!                   computing the cross-spectrum on this segment.
!
!                 For other values of TREND2 nothing is done before estimating the cross-spectrum
!                 on each segment.
!
!                 The default is TREND2=0, e.g. nothing is done before estimating the power spectrum
!                 on each segment.
!
!   WIN           (INPUT, OPTIONAL) integer(i4b)
!                 On entry, this argument specify the data window used in the computations of the
!                 power and cross spectra. If:
!
!                 - WIN=+1 The Bartlett window is used
!                 - WIN=+2 The square window is used
!                 - WIN=+3 The Welch window is used
!                 - WIN=+4 The Hann window is used
!                 - WIN=+5 The Hamming window is used
!                 - WIN=+6 A split-cosine-bell window is used
!
!                 The default is WIN=3, e.g. the Welch window is used.
!   
!   TAPERP        (INPUT, OPTIONAL) real(stnd)
!                 The total percentage of the data to be tapered if WIN=6.
!                 TAPERP must be greater than zero and less or equal to one,
!                 otherwise the default value is used.
!
!                 The default is 0.2 .
!
!   L0            (INPUT, OPTIONAL) integer(i4b)
!                 The number of zeros added to each time segment in order to obtain more finely
!                 spaced spectral estimates. L+L0 must be a positive even integer.
!
!                 The default is L0=0, e.g. no zeros are added to each time segment.
!          
!   PROBTEST      (INPUT, OPTIONAL) real(stnd)
!                 On entry, a probability. PROBTEST is the critical probability which
!                 is used to determine the lower and upper confidence limit factors (e.g.
!                 the optional arguments CONLWR and CONUPR ) and the critical value for
!                 testing the null hypothesis that the squared coherency is zero (e.g.
!                 the TESTCOHER optional argument).
!
!                 PROBTEST must verify:   0. < P < 1.
!
!                 The default is 0.05 .
!
!
! Further Details
! _______________
!
!   After removing the mean or the trend from the time series (e.g. TREND=1,2,3), the series
!   are padded with zero on the right such that the length of the resulting time series is 
!   evenly divisible by L (a positive even integer). The length, N, of these resulting time
!   series is the first integer greater than or equal to size(VEC) which is evenly divisible
!   by L. If size(VEC) is not evenly divisible by L, N is equal to size(VEC)+L-mod(size(VEC),L).
!   
!   Optionally, the mean or the trend may also be removed from each time segment (e.g. TREND2=1,2,3).
!   Optionally, zeros may be added to each time segment (e.g. the optional arguemnt L0) if more
!   finely spaced spectral esimates are desired.
!   
!   The stability of the power and cross spectra  estimates depends on the averaging process.
!   That is, the greater the number of segments ( N/L if OVERLAP=false and (2N/L)-1 if OVERLAP=true),
!   the more stable the resulting power and cross spectra estimates.
!   
!   Optionally, these power and cross spectra estimates may then be smoothed again
!   in the frequency domain by modified Daniell filters (e.g. if argument SMOOTH_PARAM is used).
!
!   The computed equivalent number of degrees of freedom and bandwidth must be divided by two
!   for the zero and Nyquist frequencies.
!
!   Furthermore, the computed equivalent number of degrees of freedom, bandwidth, lower and
!   upper (1-PROBTEST) * 100% confidence limit factors and critical value for the squared coherency
!   (e.g. arguments EDOF, BANDWIDTH, CONLWR, CONUPR and TESTCOHER) are not right near the zero
!   and Nyquist frequencies if the PSD estimates have been smoothed by modified Daniell filters.
!   The reason is that CROSS_SPECTRUM2 assumes that smoothing involves averaging independent frequency
!   ordinates. This is true except near the zero and Nyquist frequencies where an average
!   may contain contributions from negative frequencies, which are identical to and hence not
!   independent of positive frequency spectral values. Thus, the number of degrees of freedom
!   in PSD estimates near the 0 and Nyquist frequencies are as little as half the number
!   of degrees of freedom of the spectral estimates away from these frequency extremes if
!   the optional argument SMOOTH_PARAM is used.
!      
!   If the optional argument SMOOTH_PARAM is used, the computed equivalent number of degrees
!   of freedom, bandwidth, lower and upper (1-PROBTEST) * 100% confidence limit factors and
!   critical value for the squared coherency are right for PSD estimates  at frequencies
!      
!        (i-1)/(L+L0)  for i= (nparam+1)/2 + 1   to ( (L+L0) - nparam + 1)/2
!      
!   where nparam = 2 * (2+sum(SMOOTH_PARAM(:)))- 1,
!   (e.g. for frequencies i/(L+L0) for i = (nparam+1)/2, ... , ((L+L0)-nparam-1)/2 ) .
!
!   For definitions, more details and algorithm, see:
!
!   (1) Bloomfield, P., 1976:
!            Fourier analysis of time series- An introduction.
!            John Wiley and Sons, New York.
!
!   (2) Welch, P.D., 1967:
!           The use of Fast Fourier Transform for the estimation of power
!           spectra: A method based on time averaging over short, modified periodograms.
!           IEEE trans. on audio and electroacoustics, Vol. Au-15, 2, 70-73.
!
!   (3) Diggle, P.J., 1990:
!           Time series: a biostatistical introduction.
!           Clarendon Press, Oxford.
!
!      
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,         only : assert, merror, assert_eq, arth
    use Reals_Constants,   only : zero, half, one, two, four, pi, twopi, c5_m2
    use Num_Constants,     only : nan
    use Logical_Constants, only : true, false
    use Char_Constants,    only : allocate_error, tseries_error10, tseries_error15, tseries_error17,   &
                                  tseries_error18, tseries_error21, tseries_error24, tseries_error25,  &
                                  tseries_error50, tseries_error56, tseries_error58, tseries_error59
    use FFT_Procedures,    only : real_fft, init_fft, end_fft
    use Prob_Procedures,   only : probf2, pinvf2, pinvq2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout),  dimension(:)   :: vec
    real(stnd), intent(inout),  dimension(:,:) :: mat
    real(stnd), intent(out),    dimension(:)   :: psvec
    real(stnd), intent(out),    dimension(:,:) :: psmat, phase, coher
!
    real(stnd), intent(in),                  optional :: taperp, probtest
    real(stnd), intent(out),                 optional :: edof, bandwidth, conlwr, conupr, testcoher
    real(stnd), intent(out), dimension(:),   optional :: freq
    real(stnd), intent(out), dimension(:,:), optional :: ampli, co_spect, quad_spect, prob_coher
!
    integer(i4b),  intent(in)                         :: l
    integer(i4b),  intent(in),               optional :: trend, trend2, win, l0
    integer(i4b),  intent(in), dimension(:), optional :: smooth_param
!
    logical(lgl),  intent(in), optional :: normpsd, initfft, overlap
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: m, mp1, n, n2, nf, nseg, k, l2, ld2, win2, trendb, trend2b,   &
                    i, i1, i2, step, sl, sld2, nparam
    integer      :: iok
!
    real(stnd)                                :: c1, c2, sumw, edof2, probtest2, con
    real(stnd), dimension(l)                  :: wk
    real(stnd), dimension(:),   allocatable   :: temp
    real(stnd), dimension(:,:), allocatable   :: seg, magni
!
    complex(stnd), dimension(:,:), allocatable :: cwk
!
    logical(lgl)  :: normpsd2, initfft2, overlap2, smooth, out_freq,   &
                     out_ampli, out_co_spect, out_quad_spect,          &
                     out_prob_coher, out_dof
!
#ifdef _OPENMP
!
    real(stnd), dimension(:,:), allocatable :: psmatb, rcospect, icospect
!
    logical  :: test_par
!
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='cross_spectrum2'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m =  assert_eq( int(size(mat,1),i4b) ,          &
                    int(size(psmat,1),i4b) ,        &
                    int(size(phase,1),i4b) ,        &
                    int(size(coher,1),i4b) ,        &
                    name_proc )
!
    if ( m<=0_i4b )     &
    call merror( name_proc//tseries_error21 )
!
    n =  assert_eq( int(size(vec),i4b) ,          &
                    int(size(mat,2),i4b) ,        &
                    name_proc )
!
    if ( l<4_i4b )        &
    call merror( name_proc//tseries_error10 )
!
    if ( l>n )     &
    call merror( name_proc//tseries_error17 )
!
!   CHECK IF l IS A POSITIVE EVEN INTEGER.
!
    ld2 = l/2_i4b
!
    if ( l/=2_i4b*ld2 )     &
    call merror( name_proc//tseries_error18 )
!
!   DETERMINE THE LENGTH OF THE SEGMENTS.
!
    if ( present(l0) ) then
!
        if ( l0<0_i4b )     &
        call merror( name_proc//tseries_error24 )
!
        sl = l + l0
!
!       CHECK IF sl IS A POSITIVE EVEN INTEGER.
!
        sld2 = sl/2_i4b
!
        if ( sl/=2_i4b*sld2 )     &
        call merror( name_proc//tseries_error25  )
!
    else
!
        sl   = l
        sld2 = ld2
!
    end if
!
    nf  = sld2 + 1_i4b
!
    call assert( logical(int(size(psvec),i4b)==nf,lgl),     &
                 logical(int(size(psmat,2),i4b)==nf,lgl),   &
                 logical(int(size(phase,2),i4b)==nf,lgl),   &
                 logical(int(size(coher,2),i4b)==nf,lgl),   &
                 name_proc )
!
    out_freq = present( freq )
!
    if ( out_freq ) then
        call assert( logical(int(size(freq),i4b)==nf,lgl),    &
                     name_proc )
    end if
!
    out_ampli = present( ampli )
!
    if ( out_ampli ) then
        call assert( logical(int(size(ampli,1),i4b)==m,lgl),   &
                     logical(int(size(ampli,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_co_spect = present( co_spect )
!
    if ( out_co_spect ) then
        call assert( logical(int(size(co_spect,1),i4b)==m,lgl),   &
                     logical(int(size(co_spect,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_quad_spect = present( quad_spect )
!
    if ( out_quad_spect ) then
        call assert( logical(int(size(quad_spect,1),i4b)==m,lgl),   &
                     logical(int(size(quad_spect,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_prob_coher = present( prob_coher )
!
    if ( out_prob_coher ) then
        call assert( logical(int(size(prob_coher,1),i4b)==m,lgl),   &
                     logical(int(size(prob_coher,2),i4b)==nf,lgl),  &
                     name_proc )
    end if
!
    out_dof = present( edof ) .or. present( bandwidth ) .or. out_prob_coher   &
              .or. present( conlwr ) .or. present( conupr ) .or. present( testcoher ) 
!
    probtest2 = c5_m2
!
    if ( present(probtest) ) then
        if ( zero>=probtest .or. probtest>=one ) then
            call merror( name_proc//tseries_error59 )
        else
            probtest2 = probtest
        end if
    end if
!
    if ( present(trend) ) then
        trendb = trend
    else
        trendb = 1_i4b
    end if
!
    if ( present(trend2) ) then
        trend2b = trend2
    else
        trend2b = 0_i4b
    end if
!
    if ( present(win) ) then
!
        if ( win<1_i4b .or. win>6_i4b  )     &
        call merror( name_proc//tseries_error15 )
!
        win2 = win
    else
        win2 = 3_i4b
    end if
!
    initfft2 = true
!
    if ( present(initfft) ) then
        initfft2 = initfft
    end if
!
    if ( initfft2 ) then
       call init_fft( (/ m, sld2/), dim=2_i4b )
    end if
!
    normpsd2 = true
!
    if ( present(normpsd) ) then
        normpsd2 = normpsd
    end if
!
    overlap2 = false
!
    if ( present(overlap) ) then
        overlap2 = overlap
    end if
!
    smooth = false
!
    if ( present(smooth_param) ) then
!
        nparam = size( smooth_param )
!
        if ( nparam<=0_i4b  )    &
        call merror( name_proc//tseries_error50 )
!
!       CHECK THE INPUT VALUES FOR THE HALF-LENGTHS OF THE DANIEL FILTERS.
!
        if ( any( smooth_param(:nparam)<=0_i4b .or. smooth_param(:nparam)>=nf ) )   &
        call merror( name_proc//tseries_error56 )
!
        smooth = true
!
    end if
!
!   FIND THE FIRST INTEGER GREATER THAN OR EQUAL TO n THAT IS EVENLY
!   DIVISIBLE BY l.
!
    i = mod( n, l )
!
    if ( i==0_i4b ) then
        n2 = n
    else
        n2 = n + l - i
    end if
!
!   DETERMINE THE NUMBER OF SEGMENTS.
!
    if ( overlap2 ) then
        nseg = (2*n2/l) - 1_i4b
        step = ld2
    else
        nseg = n2/l
        step = l
    end if
!
    mp1 = m + 1_i4b
!
!   CALCULATE DATA WINDOW IF REQUIRED.
!
    if ( win2/=2_i4b ) then
!
!       CALCULATE DATA WINDOW.
!
        wk(:l)  = data_window( l, win2, taperp=taperp )
!
!       COMPUTE SUM OF SQUARES OF DATA WINDOW.
!
        sumw = dot_product( wk(:l), wk(:l) )
!
    else
!
!       CALCULATE RECTANGULAR WINDOW.
!
        wk(:l)  = one
!
!       COMPUTE SUM OF SQUARES OF RECTANGULAR WINDOW.
!
        sumw = real( l, stnd )
!
    end if
!
!   COMPUTE DEGREES OF FREEDOM IF REQUIRED.
!
    if ( out_dof ) then
!
        if ( nseg/=1_i4b .and. overlap2 .and. win2/=1_i4b .and. win2/=3_i4b  )    &
        call merror( name_proc//tseries_error58 )
!
        edof2 = estim_dof( wk(:l), win=win, smooth_param=smooth_param(:nparam), &
                           l0=l0, nseg=nseg, overlap=overlap )
!
    end if
!
!   REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SERIES IF REQUIRED.
!
    if ( trendb>=1_i4b .and. trendb<=3_i4b ) then
        call detrend_rv( vec(:n),    trendb )
        call detrend_rm( mat(:m,:n), trendb )
    end if
!
!   ZERO OUTPUT POWER AND CROSS SPECTRA ESTIMATES FOR SUMMING OVER THE nseg SEGMENTS.
!
    psvec(:nf)    = zero
    psmat(:m,:nf) = zero
!
!   USE THE ARRAYS coher AND phase TO STORE THE REAL AND IMAGINARY PARTS OF THE
!   CROSS SPECTRUM, RESPECTIVELY.
!
    coher(:m,:nf)  = zero
    phase(:m,:nf)  = zero
!
!   COMPUTE PSD AND CSD ESTIMATES BY SUMMING OVER THE m SEGMENTS.
!
#ifdef _OPENMP
    i1 = omp_get_num_procs()
    i2 = omp_get_max_threads()
    test_par = .not.( omp_in_parallel() )   .and.      &
               i1>1_i4b                     .and.      &
               i2>1_i4b                     .and.      &
               nseg>=i2
!
    if ( test_par ) then
!
!$OMP PARALLEL PRIVATE(k,i,i1,i2,l2,iok,psmatb,seg,cwk,rcospect,icospect)
!
!       ALLOCATE WORK MATRICES.
!
        allocate( seg(mp1,sl), cwk(mp1,nf), psmatb(mp1,nf), rcospect(m,nf), icospect(m,nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ZERO PARTIAL OUTPUT PSD AND CSD ESTIMATES FOR SUMMING OVER EACH THREAD.
!
        psmatb(:mp1,:nf) = zero
!
        rcospect(:m,:nf) = zero
        icospect(:m,:nf) = zero
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(:mp1,l+1_i4b:sl) = zero
        end if    
!
!$OMP DO SCHEDULE(STATIC)
!
        do k=1_i4b, nseg
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg(1_i4b,:l2)     = vec(i1:i2)
            seg(2_i4b:mp1,:l2) = mat(:m,i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENTS IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rm( seg(:mp1,:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:mp1,:l2) = seg(:mp1,:l2)*spread( wk(:l2), dim=1, ncopies=mp1 )
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(:mp1,l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT FOR THE SERIES.
!
            call real_fft( seg(:mp1,:sl), cwk(:mp1,:nf), true )
!
!           UPDATE POWER SPECTRUM ESIMATES.
!
            psmatb(:mp1,:sld2) = psmatb(:mp1,:sld2) + real(cwk(:mp1,:sld2),stnd)**2 + aimag(cwk(:mp1,:sld2))**2
            psmatb(:mp1,nf)    = psmatb(:mp1,nf)    + real(cwk(:mp1,nf),stnd)**2
!
!           UPDATE CO-SPECTRUM ESIMATES.
!
            rcospect(:m,:sld2)  = rcospect(:m,:sld2)                                                  +   &
                  spread(real(cwk(1_i4b,:sld2),stnd),dim=1,ncopies=m)*real(cwk(2_i4b:mp1,:sld2),stnd) +   &
                  spread(aimag(cwk(1_i4b,:sld2)),dim=1,ncopies=m)*aimag(cwk(2_i4b:mp1,:sld2))
            rcospect(:m,nf)    = rcospect(:m,nf) + real(cwk(1_i4b,nf),stnd)*real(cwk(2_i4b:mp1,nf),stnd)
!
!           UPDATE QUADRATURE-SPECTRUM ESIMATES.
!
            icospect(:m,:sld2)  = icospect(:m,:sld2)                                              +   &
                  spread(aimag(cwk(1_i4b,:sld2)),dim=1,ncopies=m)*real(cwk(2_i4b:mp1,:sld2),stnd) -   &
                  spread(real(cwk(1_i4b,:sld2),stnd),dim=1,ncopies=m)*aimag(cwk(2_i4b:mp1,:sld2))
!
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (updatepscd2)
!
        psvec(:nf)    = psvec(:nf)    + psmatb(1_i4b,:nf)
        psmat(:m,:nf) = psmat(:m,:nf) + psmatb(2_i4b:mp1,:nf)
!
        coher(:m,:nf)    = coher(:m,:nf)  + rcospect(:m,:nf) 
        phase(:m,:sld2)  = phase(:m,:sld2) + icospect(:m,:sld2)
!
!$OMP END CRITICAL (updatepscd2)
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( seg, cwk, psmatb, rcospect, icospect )
!
!$OMP END PARALLEL
!
    else
!
#endif
!
!       ALLOCATE WORK MATRICES.
!
        allocate( seg(mp1,sl), cwk(mp1,nf), stat = iok )
!
        if ( iok/=0 ) call merror( name_proc//allocate_error )
!
!       ADD ZEROS IN ORDER TO OBTAIN FINELY SPACED ESTIMATES.
!
        if ( l/=sl ) then
            seg(:mp1,l+1_i4b:sl) = zero
        end if    
!
        do k=1_i4b, nseg
!
            i  = (k-1_i4b)*step
            i1 = i + 1_i4b
            i2 = min( i + l, n )
!
            l2 = i2 - i1 + 1_i4b
!
            seg(1_i4b,:l2)     = vec(i1:i2)
            seg(2_i4b:mp1,:l2) = mat(:m,i1:i2)
!
!           REMOVE MEAN, DRIFT OR LINEAR LEAST SQUARES LINE FROM THE SEGMENTS IF REQUIRED.
!
            if ( trend2b>=1_i4b .and. trend2b<=3_i4b ) then
                call detrend_rm( seg(:mp1,:l2), trend2b )
            end if
!
!           APPLY DATA WINDOW.
!
            if ( win2/=2_i4b ) then
                seg(:mp1,:l2) = seg(:mp1,:l2)*spread( wk(:l2), dim=1, ncopies=mp1 )
            end if
!
!           PADD WITH ZEROS THE RESULTANT TIME SERIES IF NECESSARY.
!
            if ( l2/=l ) then
                seg(:mp1,l2+1_i4b:l) = zero
            end if
!
!           COMPUTE FFT FOR THE SERIES.
!
            call real_fft( seg(:mp1,:sl), cwk(:mp1,:nf), true )
!
!           UPDATE POWER SPECTRUM ESIMATES.
!
            psvec(:sld2) = psvec(:sld2)                    +   &
                           real(cwk(1_i4b,:sld2),stnd)**2  +   &
                           aimag(cwk(1_i4b,:sld2))**2
            psvec(nf)   = psvec(nf) + real(cwk(1_i4b,nf),stnd)**2
!
            psmat(:m,:sld2) = psmat(:m,:sld2)                     +   &
                              real(cwk(2_i4b:mp1,:sld2),stnd)**2  +   &
                              aimag(cwk(2_i4b:mp1,:sld2))**2
            psmat(:m,nf)   = psmat(:m,nf) + real(cwk(2_i4b:mp1,nf),stnd)**2
!
!           UPDATE CO-SPECTRUM ESIMATES.
!
            coher(:m,:sld2)  = coher(:m,:sld2)                                                            +   &
                      spread(real(cwk(1_i4b,:sld2),stnd),dim=1,ncopies=m)*real(cwk(2_i4b:mp1,:sld2),stnd) +   &
                      spread(aimag(cwk(1_i4b,:sld2)),dim=1,ncopies=m)*aimag(cwk(2_i4b:mp1,:sld2))
            coher(:m,nf)    = coher(:m,nf) + real(cwk(1_i4b,nf),stnd)*real(cwk(2_i4b:mp1,nf),stnd)
!
!           UPDATE QUADRATURE-SPECTRUM ESIMATES.
!
            phase(:m,:sld2)  = phase(:m,:sld2)                                                           +   &
                         spread(aimag(cwk(1_i4b,:sld2)),dim=1,ncopies=m)*real(cwk(2_i4b:mp1,:sld2),stnd) -   &
                         spread(real(cwk(1_i4b,:sld2),stnd),dim=1,ncopies=m)*aimag(cwk(2_i4b:mp1,:sld2))
!
        end do
!
!       DEALLOCATE WORK ARRAYS.
!
        deallocate( seg, cwk )
!
#ifdef _OPENMP
!
    end if
!
#endif
!
!   NORMALIZE THE POWER AND CROSS SPECTRA ESTIMATES.
!
    c1 = one/( sumw*real(l*nseg, stnd) )
!
    if ( normpsd2 ) then
!
        c2 = (real( l, stnd)/pi)*c1
        i1 = 1_i4b
        i2 = nf
!
    else
        c2 = two*c1
        i1 = 2_i4b
        i2 = sld2
!
        psvec(1_i4b)    = c1*psvec(1_i4b)
        psmat(:m,1_i4b) = c1*psmat(:m,1_i4b)
!
        coher(:m,1_i4b)  = c1*coher(:m,1_i4b)
        phase(:m,1_i4b)  = c1*phase(:m,1_i4b)
!
        psvec(nf)    = c1*psvec(nf)
        psmat(:m,nf) = c1*psmat(:m,nf)
!
        coher(:m,nf) = c1*coher(:m,nf)
        phase(:m,nf)  = c1*phase(:m,nf)
!
    end if
!
    psvec(i1:i2)    = c2*psvec(i1:i2)
    psmat(:m,i1:i2) = c2*psmat(:m,i1:i2)
!
    coher(:m,i1:i2)   = c2*coher(:m,i1:i2)
    phase(:m,i1:i2)   = c2*phase(:m,i1:i2)
!
!   SMOOTH THE POWER AND CROSS SPECTRA ESTIMATES IF REQUIRED.
!
    if ( smooth ) then
!
!       SMOOTH POWER SPECTRUM OF THE FIRST SERIES.
!
        call moddan_filter_rv( psvec(:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH POWER SPECTRUM OF THE OTHER SERIES.
!
        call moddan_filter_rm( psmat(:m,:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH CO-SPECTRUM.
!
        call moddan_filter_rm( coher(:m,:nf), smooth_param(:nparam), sym=one )
!
!       SMOOTH QUADRATURE-SPECTRUM.
!
        call moddan_filter_rm( phase(:m,:nf), smooth_param(:nparam), sym=-one )
!
    end if
!
    if ( out_co_spect ) then
!
!       OUTPUT THE CO-SPECTRUM IF REQUIRED.
!
        co_spect(:m,:nf) = coher(:m,:nf)
!
    end if
!
    if ( out_quad_spect ) then
!
!       OUTPUT THE QUADRATURE-SPECTRUM IF REQUIRED.
!
        quad_spect(:m,:nf) = -phase(:m,:nf)
!
    end if
!
!   ALLOCATE WORK MATRICES.
!
    allocate( temp(m), magni(m,nf), stat = iok )
!
    if ( iok/=0 ) call merror( name_proc//allocate_error )
!
    do k=1_i4b, nf
!
!       COMPUTE MAGNITUDE OF THE CROSS-SPECTRUM.
!
        magni(:m,k) = coher(:m,k)**2 + phase(:m,k)**2
!
!       COMPUTE PHASE OF THE CROSS-SPECTRUM.
!
        where ( coher(:m,k)/=zero .or. phase(:m,k)/=zero )
            temp(:m)   = atan2( phase(:m,k), coher(:m,k) )
            phase(:m,k) = (one/twopi)*temp(:m) + merge( zero, one, temp(:m)>=zero ) 
        end where
!
!       COMPUTE THE SQUARED COHERENCY.
!
        if ( psvec(k)/=zero ) then
!
            where ( psmat(:m,k)/=zero )
                coher(:m,k) = min( magni(:m,k)/(psvec(k)*psmat(:m,k)), one )
            elsewhere
                coher(:m,k) = zero
            end where
!
        else
!
            coher(:m,k) = zero     
!
        end if
!
    end do   
!
    if ( out_ampli ) then
!
!       COMPUTE AND OUTPUT AMPLITUDE OF THE CROSS-SPECTRUM IF REQUIRED.
!
        ampli(:m,:nf) = sqrt( magni(:m,:nf) )
!
    end if
!
!   OUTPUT THE FREQUENCIES IF REQUIRED.
!
    if ( out_freq ) then
!
        c1 = one/real( sl, stnd )
        freq(:nf) = arth(zero, c1, nf)
!
    end if
!
!   OUTPUT DEGREES OF FREEDOM, BANDWIDTH, CONFIDENCE LIMIT FACTORS 
!   OR PROBABILITIES FOR SQUARED COHERENCIES IF REQUIRED.
!
    if ( out_dof ) then
!
        if ( present(edof) ) then
            edof = edof2
        end if
!
        if ( present(bandwidth) ) then
            bandwidth = edof2/( two*real( n, stnd ) )
        end if
!
        c1 = probtest2*half
        c2 = one - c1
!
        sumw = nan()
!
        if ( present(conlwr) ) then
            if ( edof2>=half ) then
                conlwr = edof2/pinvq2( c2, edof2 )
            else
                conlwr = sumw
            end if
        end if
!
        if ( present(conupr) ) then
            if ( edof2>=half ) then
                conupr = edof2/pinvq2( c1, edof2 )
            else
                conupr = sumw
            end if
        end if
!
        if ( present(testcoher) ) then
            if ( edof2>two ) then
!
!                c1 = two/(edof2-two)
!                c2 = one - probtest2**c1
!                testcoher = max( min( one, c2 ), zero )
!
                 con = one - probtest2
                 c2  = edof2 - two
                 c1  = two*pinvf2( con, two, c2 )
                 testcoher = c1/( c2 + c1 )
            else
                testcoher = sumw
            end if
        end if
!
!       COMPUTE AND OUTPUT SIGNIFICANCE PROBABILITY FOR SQUARED COHERENCIES.
!
        if ( out_prob_coher ) then
!
            c2 = edof2 - two    
!
            if ( c2>zero ) then
!
                where( coher(:m,2_i4b:sld2)/=one )
                    magni(:m,2_i4b:sld2) = ( (edof2/two - one)*coher(:m,2_i4b:sld2) )/  &
                                          ( one - coher(:m,2_i4b:sld2) )
                elsewhere
                    magni(:m,2_i4b:sld2) = zero
                end where 
!            
                prob_coher(:m,2_i4b:sld2) = probf2( magni(:m,2_i4b:sld2), two, c2, true )
!            
                where( coher(:m,2_i4b:sld2)==one ) prob_coher(:m,2_i4b:sld2) = zero
!
            else
!
                prob_coher(:m,2_i4b:sld2) = sumw
!
            end if
!            
            c2 = edof2/two - two    
!
            if ( c2>zero ) then
!            
                do i1 = 1_i4b, nf, sld2
!            
                    where( coher(:m,i1)/=one )
                        temp(:m) = ( (edof2/four - one)*coher(:m,i1) )/( one - coher(:m,i1) )
                    elsewhere
                        temp(:m)  = zero
                    end where
!
                    prob_coher(:m,i1) = probf2( temp(:m), two, c2, true )
!            
                    where( coher(:m,i1)==one ) prob_coher(:m,i1) = zero
!
                end do
!
            else
!
                prob_coher(:m,1_i4b) = sumw
                prob_coher(:m,nf)    = sumw
!
            end if
!
        end if
!
    end if
!
!   DEALLOCATE WORK ARRAYS.
!
    deallocate( temp, magni )
!
!   DEALLOCATE WORK ARRAYS USED IN THE FFT COMPUTATIONS IF REQUIRED.
!
    if ( initfft2 ) then
       call end_fft( )
    end if
!
!
! END OF SUBROUTINE cross_spectrum2_rm
! ____________________________________
!
    end subroutine cross_spectrum2_rm
!
!
! =========================================================================================
!                                 PRIVATE SUBROUTINES
! =========================================================================================
!
!
    subroutine ess( y, len, ideg, njump, userw, rw, ys )
!
!
! Purpose
! _______
!                                        
!   Undocumented internal subroutine for STL_RV, STL_RM, TREND_RV, TREND_RM, STLEZ_RV
!   and STLEZ_RM.
!
!
! Arguments
! _________
!                                                                              
!   Y      (INPUT) real(stnd), dimension(:)
!
!   LEN    (INPUT) integer(i4b)
!
!   IDEG   (INPUT) integer(i4b)
!
!   NJUMP  (INPUT) integer(i4b)
!
!   USERW  (INPUT) logical(lgl)
!
!   RW     (INPUT) real(stnd), dimension(size(y))
!
!   YS     (OUPUT) real(stnd), dimension(size(y))
!
!
! Further Details
! _______________
!            
!   ESS is a private subroutine.
!
!
! _________________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:),       intent(in)  :: y
    real(stnd), dimension(size(y)), intent(in)  :: rw
    real(stnd), dimension(size(y)), intent(out) :: ys
!
    integer(i4b), intent(in)  :: len, ideg, njump
!
    logical(lgl), intent(in)  :: userw
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd)   :: delta
!
    integer(i4b) :: n, newnj, nleft, nright, nsh, k, i, j
!
    logical(lgl)  :: ok
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( y )
!
    if ( n<2_i4b ) then
!
        if ( n==1_i4b ) ys(1_i4b) = y(1_i4b)
!
        return
!
    end if
!
    newnj  = min( njump, n-1_i4b )
!
    ys(1_i4b:n) = y(1_i4b:n)
!
    if ( len>=n ) then
!
        nleft  = 1_i4b
        nright = n
!
        do i = 1_i4b, n, newnj
!
            call est( nleft, nright, n, len, ideg, y(nleft:nright),     &
                      real(i,stnd), ys(i), userw, rw(nleft:nright), ok  )
!
        end do
!
    else
!
        nleft  = 1_i4b
        nright = len
        nsh    = ( len + 1_i4b )/2_i4b
!
        if ( newnj==1_i4b ) then
!
               do i = 1_i4b, n
!
                   if (  i>nsh  .and. nright/=n ) then
                       nleft  = nleft + 1_i4b
                       nright = nright + 1_i4b
                   end if
!
                   call est( nleft, nright, n, len, ideg, y(nleft:nright),    &
                             real(i,stnd), ys(i), userw, rw(nleft:nright), ok )
               end do
!
        else
!
            do i = 1_i4b, n, newnj
!
                if ( i>=nsh ) then
!
                    if ( i>=n-nsh+1_i4b ) then
                        nleft  = n - len + 1_i4b
                        nright = n
                    else
                        nleft  = i - nsh + 1_i4b
                        nright = len + i - nsh
                    end if
!
                end if
!
                call est( nleft, nright, n, len, ideg, y(nleft:nright),     &
                          real(i,stnd), ys(i),  userw, rw(nleft:nright), ok )
!
            end do
!
        end if
!
    end if
!
    if ( newnj/=1_i4b ) then
!
        do i = 1_i4b, n-newnj, newnj
!
            delta = ( ys(i+newnj) - ys(i) )/real( newnj, stnd )
!
            do j = i+1_i4b, i+newnj-1_i4b
                ys(j) = ys(i) + delta*real( j-i, stnd )
            end do
!
        end do
!
        k = ( (n-1_i4b)/newnj )*newnj + 1_i4b
!
        if ( k/=n ) then
!
            call est( nleft, nright, n, len, ideg, y(nleft:nright),     &
                      real(n,stnd), ys(n),  userw, rw(nleft:nright), ok )
!
            if ( k/=n-1_i4b ) then
!
                delta = ( ys(n) - ys(k) )/real( n-k, stnd )
!
                do j = k+1_i4b, n-1_i4b
                    ys(j) = ys(k) + delta*real( j-k , stnd )
                end do
!
            end if
!
        end if
!
    end if
!
!
! END OF SUBROUTINE ess
! _____________________
!
    end subroutine ess
!
! ===========================================================================================
!
    subroutine est( nleft, nright, n, len, ideg, y, xs, ys, userw, rw, ok )
!
! Purpose
! _______
!                                                                              
!   Undocumented internal subroutine for STL_RV, STL_RM, TREND_RV, TREND_RM, STLEZ_RV
!   and STLEZ_RM.
!
!
! Arguments
! _________
!                                                                              
!   NLEFT  (INPUT) integer(i4b)
!
!   NRIGHT (INPUT) integer(i4b)
!
!   N      (INPUT) integer(i4b)
!
!   LEN    (INPUT) integer(i4b)
!
!   IDEG   (INPUT) integer(i4b)
!
!   Y      (INPUT) real(stnd), dimension(nleft:nright)
!
!   XS     (INPUT) real(stnd)
!
!   YS     (INPUT/OUTPUT) real(stnd)
!
!   USERW  (INPUT) logical(lgl)
!
!   RW     (INPUT) real(stnd), dimension(nleft:nright)
!
!   OK     (OUTPUT) logical(lgl)
!
!
! Further Details
! _______________
!            
!   EST is a private subroutine.
!
!
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Reals_Constants,   only : zero, one, c0_999, c1_m3, c1_m4
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in)  :: nleft, nright, n, len, ideg
!
    real(stnd),                          intent(in)     :: xs
    real(stnd),                          intent(inout)  :: ys
    real(stnd), dimension(nleft:nright), intent(in)     :: y, rw
!
    logical(lgl), intent(in)  :: userw
    logical(lgl), intent(out) :: ok
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd)                          :: range, h, h1, h9, a, b, c, r, mx, mx2,   &
                                           covyx, covyx2, covxx2, varx, varx2
    real(stnd), dimension(nleft:nright) :: w, x
!
    integer(i4b) :: j, ngood
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    range = real( n, stnd ) - one
    h = max( xs-real(nleft,stnd), real(nright,stnd)-xs )
    if ( len>n ) h = h + real( (len-n)/2_i4b, stnd )
!
    h9 = c0_999*h
    h1 = c1_m3*h
!
!   COMPUTE WEIGHTS.
!
    do j = nleft, nright
!
        x(j) = real( j, stnd )
        r = abs( x(j) - xs )
!
        if ( r<=h9 ) then
!
            if ( r<=h1 ) then
                w(j) = one
            else
                w(j) = (one-(r/h)**3_i4b)**3_i4b
            end if
!
        else
!
            w(j) = zero
!
        end if
!
    end do
!
    if ( userw ) then
        w(nleft:nright) = rw(nleft:nright)*w(nleft:nright)
    end if
!
    a = sum( w(nleft:nright) )
    ok = a>zero
!
    if ( ok ) then
!
        ngood = count( w(nleft:nright)/=zero )
!
!       MAKE SUM OF WEIGHTS EQUAL TO ONE.
!
        w(nleft:nright) = w(nleft:nright)/a
!
!       COMPUTE WEIGHTED MEAN OF THE DATA.
!
        ys = dot_product( w(nleft:nright), y(nleft:nright) )
!
        if ( (h>zero) .and. (ideg>0_i4b) .and. (ngood>2_i4b) ) then
!
            mx = dot_product( w(nleft:nright), x(nleft:nright) )
!
            if ( ideg==1_i4b ) then
!
!               USE LINEAR FIT.
!
                varx  = zero
                covyx = zero
!
                do j = nleft, nright
!
                    a = x(j) - mx
!
                    varx  = varx  + w(j)*(a**2_i4b)
                    covyx = covyx + w(j)*( a*(y(j)-ys) )
!
                end do
!
                if ( sqrt(varx)>c1_m4*range ) then
!
!                   POINTS ARE SPREAD OUT ENOUGH TO COMPUTE SLOPE.
!
                    a  = covyx/varx
                    ys = a*(xs-mx) + ys
!
                end if
!
            else
!
!               USE QUADRATIC FIT.
!
                mx2 = dot_product( w(nleft:nright), x(nleft:nright)**2 )
!
                varx   = zero
                varx2  = zero
                covyx  = zero
                covyx2 = zero
                covxx2 = zero
!
                do j = nleft, nright
!
                    a     = x(j) - mx
                    varx  = varx + w(j)*(a*a)
                    b     = x(j)*x(j) - mx2
                    varx2 = varx2 + w(j)*(b*b)
                    c     = y(j) - ys
!
                    covyx  = covyx + w(j)*(a*c)
                    covyx2 = covyx2 + w(j)*(b*c)
                    covxx2 = covxx2 + w(j)*(a*b)
!
                end do
!
                r = varx*varx2 - covxx2**2
!
                if ( abs(r)>c1_m4 ) then
!
!                   QUADRATIC FIT IS POSSIBLE.
!
                    a  = ( (varx*covyx2) - (covyx*covxx2) )/r
                    b  = ( (varx2*covyx) - (covyx2*covxx2) )/r
                    c  = ys - b*mx - a*mx2
                    ys = (a*xs + b)*xs + c
!
                else if ( sqrt(varx)>c1_m4*range ) then
!
!                   USE LINEAR FIT OTHERWISE.
!
                    a  = covyx/varx
                    ys = a*(xs-mx) + ys
!
                end if
!
            end if
!
        end if
!
    end if
!
!
! END OF SUBROUTINE est
! _____________________
!
    end subroutine est
!
! ===========================================================================================
!
    subroutine rwts( y, fit, rw )
!
! Purpose
! _______
!                                          
!   Compute robustness weights for the stl or trend analysis of vector Y.
!
!
! Arguments
! _________
!                                                                              
!   Y    (INPUT) real(stnd), dimension(:)
!        On entry, the input time series.
!
!   FIT  (INPUT) real(stnd), dimension(size(y))
!        On entry, one component of the TREND or STL analysis of the vector Y.
!
!
!   RW   (OUTPUT) real(stnd), dimension(size(y))
!        On output, the robustness weights for each observations of Y.
!
!
! Further Details
! _______________
!            
!   RWTS is a private subroutine.
!
!                                                                                                 
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Reals_Constants,  only : zero, one, six, c0_999, c1_m3
    use Stat_Procedures,  only : valmed
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:),       intent(in)  :: y
    real(stnd), dimension(size(y)), intent(in)  :: fit
    real(stnd), dimension(size(y)), intent(out) :: rw
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd)   :: cmad, c9, c1
!
    integer(i4b) :: n
!
    logical(lgl) , dimension(size(y)) :: maskrw
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( y )
!
    rw(1_i4b:n) = abs( y(1_i4b:n) - fit(1_i4b:n) )
    cmad = six*valmed( rw(1_i4b:n) )
!
    c9 = c0_999*cmad
    c1 = c1_m3*cmad
!
    maskrw(1_i4b:n) = rw(1_i4b:n)>c1   
    where ( maskrw(1_i4b:n) .and. rw(1_i4b:n)<=c9 )
         rw(1_i4b:n) = ( one - (rw(1_i4b:n)/cmad)**2_i4b )**2_i4b
    elsewhere
         rw(1_i4b:n) = merge( zero, one, maskrw(1_i4b:n) )
    end where
!
!
! END OF SUBROUTINE rwts
! ______________________
!
    end subroutine rwts
!
! ===========================================================================================
!
    subroutine ss( y, np, ns, isdeg, nsjump, userw, rw, season )
!
! Purpose
! _______
!                                                                              
!   Undocumented internal subroutine for slt_rv, stl_rm, stlez_rv, stlez_rm.
!
!
! Arguments
! _________
!                                                                              
!   Y      (INPUT) real(stnd), dimension(:)
!
!   NP     (INPUT) integer(i4b)
!
!   NS     (INPUT) integer(i4b)
!
!   ISDEG  (INPUT) integer(i4b)
!
!   NSJUMP (INPUT) integer(i4b)
!
!   USERW  (INPUT) logical(lgl)
!
!   RW     (INPUT) real(stnd), dimension(size(y))
!
!   SEASON (OUTPUT) real(stnd), dimension(size(y)+2*np)
!
!
! Further Details
! _______________
!            
!   SS is a private subroutine.
!
!
! _________________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Reals_Constants,  only : zero
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in)  :: np, ns, isdeg, nsjump
!
    real(stnd), dimension(:),            intent(in)  :: y
    real(stnd), dimension(size(y)),      intent(in)  :: rw
    real(stnd), dimension(size(y)+2*np), intent(out) :: season
!
    logical(lgl), intent(in)  :: userw
!
!
! SPECIFICATIONS FOR VARIABLES
! ____________________________
!
    real(stnd)                     :: xs
    real(stnd), dimension(size(y)) :: work1, work2, work3
!
    integer(i4b) :: i, j, k, m, n, nright, nleft
!
    logical(lgl) :: ok
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( y )
!
    do j = 1_i4b, np
!
        k = (n-j)/np + 1_i4b
!
        do i = 1_i4b, k
            work1(i) = y((i-1_i4b)*np+j)
        end do
!
        if ( userw ) then
!
            do i = 1_i4b, k
                work3(i) = rw((i-1_i4b)*np+j)
            end do
!
        end if
!
        call ess( work1(1_i4b:k), ns, isdeg, nsjump, userw,     &
                  work3(1_i4b:k), work2(2_i4b:k+1_i4b) )
!
        xs = zero
        nright = min( ns, k )
        work2(1_i4b) = work2(2_i4b)
!
        call est( 1_i4b, nright, k, ns, isdeg, work1(1_i4b:nright), xs,     &
                  work2(1_i4b),  userw, work3(1_i4b:nright), ok )
!
        xs = real( k+1_i4b, stnd )
        nleft = max( 1_i4b, k-ns+1_i4b)
        work2(k+2_i4b) = work2(k+1_i4b)
!
        call est( nleft, k, k, ns, isdeg, work1(nleft:k), xs, work2(k+2_i4b), userw,    &
                  work3(nleft:k), ok )
!
        do m = 1_i4b, k+2_i4b
!
            season((m-1_i4b)*np+j) = work2(m)
!
        end do
!
    end do
!
!
! END OF SUBROUTINE ss
! ____________________
!
    end subroutine ss
!
!
! =========================================================================================
!
!
! **************************************
! END OF MODULE Time_Series_Procedures *
! **************************************
!

