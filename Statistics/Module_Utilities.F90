module Utilities
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
! MODULE EXPORTING GENERAL AND COMPUTING UTILITIES.                            *
!                                                                              * 
! MANY OF THESE ROUTINES ARE ADAPTED AND EXTENDED FROM PUBLIC DOMAIN ROUTINES  *
! FROM Numerical Recipes.                                                      *
!                                                                              *
! LATEST REVISION : 24/09/2020                                                 *
!                                                                              *
! ##############################################################################
!                                                                              *
! ******************************************************************************
!                                                                              *
! OPENMP DIRECTIVES IN THIS MODULE ARE ACTIVATED WITH THE C                    *
! PROCESSOR MACRO _OPENMP3 AND, BY DEFAULT, THESE DIRECTIVES                   *
! ARE ACTIVATED IF OPENMP IS USED.                                             *
!                                                                              *
! HOWEVER, SINCE PARALLELIZATION IS GENERALLY DONE AT HIGHER                   *
! LEVEL IN THE STATPACK SOFTWARE, THESE OPENMP DIRECTIVES MAY                  *
! BE DESACTIVATED WITH THE C PROCESSOR MACRO _NOOPENMP3 DURING                 *
! COMPILATION.                                                                 *
!                                                                              *
! MORE GENERALLY, THE C PROCESSOR MACROS USED IN THIS MODULE ARE:              *
!                                                                              *
!  _OPENMP3     FOR ACTIVATING OPENMP PARALLELIZATION                          *
!  _F2003       FOR ACTIVATING FORTRAN 2003 CONSTRUCTS                         *
!  _BLAS        FOR USING BLAS WHEN POSSIBLE                                   *
!  _DOT_PRODUCT FOR REPLACING THE dot_product INTRINSIC FUNCTION               *
!               WITH STATPACK dot_product2 FUNCTION OR OTHER                   *
!               CONSTRUCTS                                                     *
!                                                                              *
! ******************************************************************************
!
#ifdef _OPENMP
!
#ifndef _NOOPENMP3
#define _OPENMP3
#endif
!
#else
!
#undef _OPENMP3
!
#endif
!
! USED MODULES
! ============
!
    use Select_Parameters,  only : lgl, i2b, i4b, stnd,                    &
                                   defunit, blksz_util,                    &
                                   npar_arth, npar2_arth,                  &
                                   npar_geop, npar2_geop,                  &
                                   npar_cumsum, npar_cumprod,              &
                                   npar_poly, npar_polyterm
    use Logical_Constants
#ifdef _OPENMP3
    use Select_Parameters,  only : omp_limit, omp_limit2
    use omp_lib,            only : omp_get_max_threads, omp_in_parallel
#endif
!
!   use Reals_Constants
!   use Char_Constants
!   use Num_Constants
!   use BLAS_interfaces
!
! STRONG TYPING IMPOSED 
! =====================
!    
    implicit none
!
! PUBLIC ENTITIES 
! ===============
!    
! ALL SUBROUTINES, FUNCTIONS, VARIABLES AND PARAMETERS ARE PRIVATE BY DEFAULT.
!
    private
    public :: array_copy, swap,                                               &
              ifirstloc, imaxloc, iminloc,                                    &
              assert, assert_eq, merror,                                      &
              arth, geop, cumsum, cumprod,                                    &
              poly, poly_term, zroots_unity,                                  &
              update_rk1, update_rk2,                                         &
              outerprod, outerdiv,                                            &
              outersum, outerdiff,                                            &
              outerand, outeror, triangle,                                    &
              abse, lassq, norm,                                              &
              scatter_add, scatter_max,                                       &
              diagadd, diagmult, get_diag, put_diag,                          &
              unit_matrix,                                                    &
              lascl, norme, pythag,                                           &
              transpose2, dot_product2, mmproduct, matmul2
#ifdef _F2003
    public :: mvalloc
#endif
!
! GENERIC INTERFACES FOR ROUTINES WITH OVERLOADED VERSIONS
! ========================================================
!
    interface transpose2
        module procedure    transpose_rm, transpose_cm, transpose_im, transpose_lm
    end interface
!
    interface dot_product2
        module procedure    dot_product_rv, dot_product_cv, dot_product_iv, dot_product_lv
    end interface
!
    interface mmproduct
        module procedure    mmproduct_rv_rm,  mmproduct_rm_rv, mmproduct_rm_rm,  &
                            mmproduct_cv_cm,  mmproduct_cm_cv, mmproduct_cm_cm
    end interface
!
    interface matmul2
        module procedure    matmul_rv_rm,  matmul_rm_rv, matmul_rm_rm,            &
                            matmul_cv_cm,  matmul_cm_cv, matmul_cm_cm,            &
                            matmul_iv_im,  matmul_im_iv, matmul_im_im,            &
                            matmul_lv_lm,  matmul_lm_lv, matmul_lm_lm
    end interface
!
    interface array_copy
        module procedure    array_copy_iv, array_copy_rv, array_copy_cv
    end interface
!
    interface swap
        module procedure     swap_i, swap_iv, swap_im,                        &
                             swap_r, swap_rv, swap_rm,                        &
                             swap_c, swap_cv, swap_cm,                        &
                             masked_swap_i, masked_swap_iv, masked_swap_im,   &
                             masked_swap_r, masked_swap_rv, masked_swap_rm,   &
                             masked_swap_c, masked_swap_cv, masked_swap_cm
    end interface
!
#ifdef _F2003
    interface mvalloc
        module procedure     mvalloc_iv, mvalloc_im,   &
                             mvalloc_rv, mvalloc_rm,   &
                             mvalloc_cv, mvalloc_cm,   &
                             mvalloc_hv
    end interface
#endif
!
    interface imaxloc
        module procedure     imaxloc_i, masked_imaxloc_i,    &
                             imaxloc_r, masked_imaxloc_r
    end interface
!
    interface iminloc
        module procedure     iminloc_i, masked_iminloc_i,    &
                             iminloc_r, masked_iminloc_r
    end interface
!
    interface assert
        module procedure     assert1, assert2, assert3, assert4, assert_v
    end interface
!
    interface assert_eq
        module procedure     assert_eq2, assert_eq3, assert_eq4, assert_eqn
    end interface
!
    interface arth
        module procedure     arth_i, arth_iv, arth_r, arth_rv, arth_c, arth_cv
    end interface
!
    interface geop
        module procedure     geop_i, geop_iv, geop_r, geop_rv, geop_c, geop_cv       
    end interface
!
    interface cumsum
        module procedure     cumsum_i, cumsum_r, cumsum_c
    end interface
!
    interface cumprod
        module procedure     cumprod_i, cumprod_r, cumprod_c
    end interface
!
    interface poly
        module procedure     poly_rr, poly_rrv, poly_rc, poly_cc, poly_msk_rrv
    end interface
!
    interface poly_term
        module procedure     poly_term_rr, poly_term_cc
    end interface
!
    interface update_rk1
        module procedure     update_rk1_i, update_rk1_r, update_rk1_c
    end interface
!
    interface update_rk2
        module procedure     update_rk2_i, update_rk2_r, update_rk2_c
    end interface
!
    interface outerprod
        module procedure     outerprod_i, outerprod_r, outerprod_c
    end interface
!
    interface outerdiv
        module procedure     outerdiv_r, outerdiv_c
    end interface
!
    interface outersum
        module procedure     outersum_i, outersum_r, outersum_c
    end interface
!
    interface outerdiff
        module procedure     outerdiff_i, outerdiff_r, outerdiff_c
    end interface
!    
    interface abse
        module procedure     abse_rv, abse_rm, abse_dim_rm,          &
                             abse_cv, abse_cm, abse_dim_cm
    end interface
!    
    interface lassq
        module procedure     lassq_rv, lassq_rm,        &
                             lassq_cv, lassq_cm
    end interface
!    
    interface norm
        module procedure     norm_rv, norm_rm, norm_dim_rm,      &
                             norm_cv, norm_cm, norm_dim_cm
    end interface
!    
! SERIAL ALGORITHMS (BUT MAY BE PARALLELIZED WITH OpenMP)
!
    interface scatter_add
        module procedure     scatter_add_i, scatter_add_r, scatter_add_c
    end interface
!
    interface scatter_max
        module procedure     scatter_max_i, scatter_max_r
    end interface
!
    interface diagadd
        module procedure     diagadd_rv, diagadd_r, diagadd_cv, diagadd_c
    end interface
!    
    interface diagmult
        module procedure     diagmult_rv, diagmult_r, diagmult_cv, diagmult_c
    end interface
!    
    interface get_diag
        module procedure     get_diag_r, get_diag_c
    end interface
!    
    interface put_diag
        module procedure     put_diag_rv, put_diag_r, put_diag_cv, put_diag_c
    end interface
!    
    interface unit_matrix
        module procedure     unit_matrix_r, unit_matrix_c
    end interface
!    
!    
    interface lascl
        module procedure     lascl_r,  lascl_rv, lascl_rm, type_lascl_rm,           &
                             masked_lascl_r, masked_lascl_rv, masked_lascl_rm
    end interface
!    
    interface norme
        module procedure    norme_rv, norme_rm
    end interface
!
!
! MODULE VARIABLES
! ================
!    
!
!   
! MODULE PARAMETERS 
! =================
!
! BLOCK SIZE FOR THE OPTIMIZED MULTIPLICATION ALGORITHMS IN matmul_rm_rm, matmul_cm_cm,
! matmul_im_im AND matmul_lm_lm SUBROUTINES AND FOR THE OPTIMIZED TRANSPOSITION ALGORITHMS
! IN transpose_rm, transpose_cm, transpose_im AND transpose_lm FUNCTIONS.
!
    integer(i4b), parameter :: blksz = blksz_util
!
!
! =========================================================================================
!
                                  contains
!                                 ========
!
! =========================================================================================
!                 ROUTINES WHICH MAY BE USED TO REPLACE INTRINSIC PROCEDURES.
! =========================================================================================
!
    function transpose_rm( mat )
!
! Purpose
! _______
!
!   Transpose the real matrix MAT
!
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:,:), intent(in)         :: mat
    real(stnd), dimension(size(mat,2),size(mat,1)) :: transpose_rm
!
    integer(i4b) :: i, i2, j, j2, n, m, p, rem
#ifdef _OPENMP3
    logical      :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
    p = max( n, m )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( p<=blksz ) then
!
        transpose_rm(:m,:n)  = transpose( mat(:n,:m) )
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
!       PARALLEL AND BLOCKED TRANSPOSITION OF mat.
!
        if ( n/=m ) then
!
!           RECTANGULAR MATRIX.
!
!$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(i,i2,j,j2), COLLAPSE(2)
!
            do j = 1_i4b, m, blksz
!
                do i = 1_i4b, n, blksz
!
                    j2 = min( m, j+blksz-1_i4b )
                    i2 = min( n, i+blksz-1_i4b )
!
                    transpose_rm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
!
                end do
!
            end do
!
!$OMP END PARALLEL DO
!  
        else
!
!           SQUARE MATRIX.
!
!$OMP PARALLEL PRIVATE(i,i2,j,j2,p,rem,m)
!
            p   = n/blksz
            rem = n - p*blksz
            m   = blksz*( p - 1_i4b ) + 1_i4b
!
!$OMP DO SCHEDULE(STATIC)
!
            do j = 1_i4b, m, blksz
!
                j2 = j + blksz - 1_i4b
!
                transpose_rm(j:j2,j:j2) = transpose( mat(j:j2,j:j2) )
!
                do i = j+blksz, m, blksz
!
                    i2 = i + blksz - 1_i4b
!
                    transpose_rm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
                    transpose_rm(i:i2,j:j2) = transpose( mat(j:j2,i:i2) )
!
                end do
!
            end do
!
!$OMP END DO NOWAIT
!
            if ( rem>0_i4b ) then
!
                m = p*blksz + 1_i4b
!
!$OMP DO SCHEDULE(STATIC)
!
                do j = 1_i4b, n/rem-1_i4b
!
                    i  = (j-1_i4b)*rem + 1_i4b
                    i2 = j*rem
!
                    transpose_rm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                    transpose_rm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
                end do
!
!$OMP END DO NOWAIT
!
!$OMP SINGLE
!
                i  = ( n/rem - 2_i4b )*rem + 1_i4b
                i2 = m - 1_i4b
!
                transpose_rm(m:n,m:n)  = transpose( mat(m:n,m:n) )
                transpose_rm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                transpose_rm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
!$OMP END SINGLE
!
            end if
!
!$OMP END PARALLEL
!
        end if
!
    else
#endif
!
!       BLOCKED TRANSPOSITION OF mat.
!
        if ( n/=m ) then
!
!           RECTANGULAR MATRIX.
!
            do j = 1_i4b, m, blksz
!
                j2 = min( m, j+blksz-1_i4b )
!
                do i = 1_i4b, n, blksz
!
                    i2 = min( n, i+blksz-1_i4b )
!
                    transpose_rm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
!
                end do
!
            end do
!  
        else
!
!           SQUARE MATRIX.
!
            p   = n/blksz
            rem = n - p*blksz
            m   = blksz*( p - 1_i4b ) + 1_i4b
!
            do j = 1_i4b, m, blksz
!
                j2 = j + blksz - 1_i4b
!
                transpose_rm(j:j2,j:j2) = transpose( mat(j:j2,j:j2) )
!
                do i = j+blksz, m, blksz
!
                    i2 = i + blksz - 1_i4b
!
                    transpose_rm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
                    transpose_rm(i:i2,j:j2) = transpose( mat(j:j2,i:i2) )
!
                end do
!
            end do
!
            if ( rem>0_i4b ) then
!
                m = p*blksz + 1_i4b
!
                do j = 1_i4b, n/rem-1_i4b
!
                    i  = (j-1_i4b)*rem + 1_i4b
                    i2 = j*rem
!
                    transpose_rm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                    transpose_rm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
                end do
!
                i  = ( n/rem - 2_i4b )*rem + 1_i4b
                i2 = m - 1_i4b
!
                transpose_rm(m:n,m:n)  = transpose( mat(m:n,m:n) )
                transpose_rm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                transpose_rm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
            end if
!
        end if
!
#ifdef _OPENMP3
    end if
#endif
!  
    end function transpose_rm
!
! =========================================================================================
!
    function transpose_cm( mat )
!
! Purpose
! _______
!
!   Transpose the complex matrix MAT
!
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:,:), intent(in)         :: mat
    complex(stnd), dimension(size(mat,2),size(mat,1)) :: transpose_cm
!
    integer(i4b) :: i, i2, j, j2, n, m, p, rem
#ifdef _OPENMP3
    logical      :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
    p = max( n, m )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( p<=blksz ) then
!
        transpose_cm(:m,:n)  = transpose( mat(:n,:m) )
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
!       PARALLEL AND BLOCKED TRANSPOSITION OF mat.
!
        if ( n/=m ) then
!
!           RECTANGULAR MATRIX.
!
!$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(i,i2,j,j2), COLLAPSE(2)
!
            do j = 1_i4b, m, blksz
!
                do i = 1_i4b, n, blksz
!
                    j2 = min( m, j+blksz-1_i4b )
                    i2 = min( n, i+blksz-1_i4b )
!
                    transpose_cm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
!
                end do
!
            end do
!
!$OMP END PARALLEL DO
!  
        else
!
!           SQUARE MATRIX.
!
!$OMP PARALLEL PRIVATE(i,i2,j,j2,p,rem,m)
!
            p   = n/blksz
            rem = n - p*blksz
            m   = blksz*( p - 1_i4b ) + 1_i4b
!
!$OMP DO SCHEDULE(STATIC)
!
            do j = 1_i4b, m, blksz
!
                j2 = j + blksz - 1_i4b
!
                transpose_cm(j:j2,j:j2) = transpose( mat(j:j2,j:j2) )
!
                do i = j+blksz, m, blksz
!
                    i2 = i + blksz - 1_i4b
!
                    transpose_cm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
                    transpose_cm(i:i2,j:j2) = transpose( mat(j:j2,i:i2) )
!
                end do
!
            end do
!
!$OMP END DO NOWAIT
!
            if ( rem>0_i4b ) then
!
                m = p*blksz + 1_i4b
!
!$OMP DO SCHEDULE(STATIC)
!
                do j = 1_i4b, n/rem-1_i4b
!
                    i  = (j-1_i4b)*rem + 1_i4b
                    i2 = j*rem
!
                    transpose_cm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                    transpose_cm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
                end do
!
!$OMP END DO NOWAIT
!
!$OMP SINGLE
!
                i  = ( n/rem - 2_i4b )*rem + 1_i4b
                i2 = m - 1_i4b
!
                transpose_cm(m:n,m:n)  = transpose( mat(m:n,m:n) )
                transpose_cm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                transpose_cm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
!$OMP END SINGLE
!
            end if
!
!$OMP END PARALLEL
!
        end if
!
    else
#endif
!
!       BLOCKED TRANSPOSITION OF mat.
!
        if ( n/=m ) then
!
!           RECTANGULAR MATRIX.
!
            do j = 1_i4b, m, blksz
!
                j2 = min( m, j+blksz-1_i4b )
!
                do i = 1_i4b, n, blksz
!
                    i2 = min( n, i+blksz-1_i4b )
!
                    transpose_cm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
!
                end do
!
            end do
!  
        else
!
!           SQUARE MATRIX.
!
            p   = n/blksz
            rem = n - p*blksz
            m   = blksz*( p - 1_i4b ) + 1_i4b
!
            do j = 1_i4b, m, blksz
!
                j2 = j + blksz - 1_i4b
!
                transpose_cm(j:j2,j:j2) = transpose( mat(j:j2,j:j2) )
!
                do i = j+blksz, m, blksz
!
                    i2 = i + blksz - 1_i4b
!
                    transpose_cm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
                    transpose_cm(i:i2,j:j2) = transpose( mat(j:j2,i:i2) )
!
                end do
!
            end do
!
            if ( rem>0_i4b ) then
!
                m = p*blksz + 1_i4b
!
                do j = 1_i4b, n/rem-1_i4b
!
                    i  = (j-1_i4b)*rem + 1_i4b
                    i2 = j*rem
!
                    transpose_cm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                    transpose_cm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
                end do
!
                i  = ( n/rem - 2_i4b )*rem + 1_i4b
                i2 = m - 1_i4b
!
                transpose_cm(m:n,m:n)  = transpose( mat(m:n,m:n) )
                transpose_cm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                transpose_cm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
            end if
!
        end if
!
#ifdef _OPENMP3
    end if
#endif
!  
    end function transpose_cm
!
! =========================================================================================
!
    function transpose_im( mat )
!
! Purpose
! _______
!
!   Transpose the integer matrix MAT
!
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:,:), intent(in)         :: mat
    integer(i4b), dimension(size(mat,2),size(mat,1)) :: transpose_im
!
    integer(i4b) :: i, i2, j, j2, n, m, p, rem
#ifdef _OPENMP3
    logical      :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
    p = max( n, m )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( p<=blksz ) then
!
        transpose_im(:m,:n)  = transpose( mat(:n,:m) )
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
!       PARALLEL AND BLOCKED TRANSPOSITION OF mat.
!
        if ( n/=m ) then
!
!           RECTANGULAR MATRIX.
!
!$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(i,i2,j,j2), COLLAPSE(2)
!
            do j = 1_i4b, m, blksz
!
                do i = 1_i4b, n, blksz
!
                    j2 = min( m, j+blksz-1_i4b )
                    i2 = min( n, i+blksz-1_i4b )
!
                    transpose_im(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
!
                end do
!
            end do
!
!$OMP END PARALLEL DO
!  
        else
!
!           SQUARE MATRIX.
!
!$OMP PARALLEL PRIVATE(i,i2,j,j2,p,rem,m)
!
            p   = n/blksz
            rem = n - p*blksz
            m   = blksz*( p - 1_i4b ) + 1_i4b
!
!$OMP DO SCHEDULE(STATIC)
!
            do j = 1_i4b, m, blksz
!
                j2 = j + blksz - 1_i4b
!
                transpose_im(j:j2,j:j2) = transpose( mat(j:j2,j:j2) )
!
                do i = j+blksz, m, blksz
!
                    i2 = i + blksz - 1_i4b
!
                    transpose_im(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
                    transpose_im(i:i2,j:j2) = transpose( mat(j:j2,i:i2) )
!
                end do
!
            end do
!
!$OMP END DO NOWAIT
!
            if ( rem>0_i4b ) then
!
                m = p*blksz + 1_i4b
!
!$OMP DO SCHEDULE(STATIC)
!
                do j = 1_i4b, n/rem-1_i4b
!
                    i  = (j-1_i4b)*rem + 1_i4b
                    i2 = j*rem
!
                    transpose_im(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                    transpose_im(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
                end do
!
!$OMP END DO NOWAIT
!
!$OMP SINGLE
!
                i  = ( n/rem - 2_i4b )*rem + 1_i4b
                i2 = m - 1_i4b
!
                transpose_im(m:n,m:n)  = transpose( mat(m:n,m:n) )
                transpose_im(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                transpose_im(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
!$OMP END SINGLE
!
            end if
!
!$OMP END PARALLEL
!
        end if
!
    else
#endif
!
!       BLOCKED TRANSPOSITION OF mat.
!
        if ( n/=m ) then
!
!           RECTANGULAR MATRIX.
!
            do j = 1_i4b, m, blksz
!
                j2 = min( m, j+blksz-1_i4b )
!
                do i = 1_i4b, n, blksz
!
                    i2 = min( n, i+blksz-1_i4b )
!
                    transpose_im(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
!
                end do
!
            end do
!  
        else
!
!           SQUARE MATRIX.
!
            p   = n/blksz
            rem = n - p*blksz
            m   = blksz*( p - 1_i4b ) + 1_i4b
!
            do j = 1_i4b, m, blksz
!
                j2 = j + blksz - 1_i4b
!
                transpose_im(j:j2,j:j2) = transpose( mat(j:j2,j:j2) )
!
                do i = j+blksz, m, blksz
!
                    i2 = i + blksz - 1_i4b
!
                    transpose_im(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
                    transpose_im(i:i2,j:j2) = transpose( mat(j:j2,i:i2) )
!
                end do
!
            end do
!
            if ( rem>0_i4b ) then
!
                m = p*blksz + 1_i4b
!
                do j = 1_i4b, n/rem-1_i4b
!
                    i  = (j-1_i4b)*rem + 1_i4b
                    i2 = j*rem
!
                    transpose_im(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                    transpose_im(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
                end do
!
                i  = ( n/rem - 2_i4b )*rem + 1_i4b
                i2 = m - 1_i4b
!
                transpose_im(m:n,m:n)  = transpose( mat(m:n,m:n) )
                transpose_im(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                transpose_im(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
            end if
!
        end if
!
#ifdef _OPENMP3
    end if
#endif
!  
    end function transpose_im
!
! =========================================================================================
!
    function transpose_lm( mat )
!
! Purpose
! _______
!
!   Transpose the logical matrix MAT
!
!
! __________________________________________________________________________________________
!
    logical(lgl), dimension(:,:), intent(in)         :: mat
    logical(lgl), dimension(size(mat,2),size(mat,1)) :: transpose_lm
!
    integer(i4b) :: i, i2, j, j2, n, m, p, rem
#ifdef _OPENMP3
    logical      :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
    p = max( n, m )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( p<=blksz ) then
!
        transpose_lm(:m,:n)  = transpose( mat(:n,:m) )
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
!       PARALLEL AND BLOCKED TRANSPOSITION OF mat.
!
        if ( n/=m ) then
!
!           RECTANGULAR MATRIX.
!
!$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(i,i2,j,j2), COLLAPSE(2)
!
            do j = 1_i4b, m, blksz
!
                do i = 1_i4b, n, blksz
!
                    j2 = min( m, j+blksz-1_i4b )
                    i2 = min( n, i+blksz-1_i4b )
!
                    transpose_lm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
!
                end do
!
            end do
!
!$OMP END PARALLEL DO
!  
        else
!
!           SQUARE MATRIX.
!
!$OMP PARALLEL PRIVATE(i,i2,j,j2,p,rem,m)
!
            p   = n/blksz
            rem = n - p*blksz
            m   = blksz*( p - 1_i4b ) + 1_i4b
!
!$OMP DO SCHEDULE(STATIC)
!
            do j = 1_i4b, m, blksz
!
                j2 = j + blksz - 1_i4b
!
                transpose_lm(j:j2,j:j2) = transpose( mat(j:j2,j:j2) )
!
                do i = j+blksz, m, blksz
!
                    i2 = i + blksz - 1_i4b
!
                    transpose_lm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
                    transpose_lm(i:i2,j:j2) = transpose( mat(j:j2,i:i2) )
!
                end do
!
            end do
!
!$OMP END DO NOWAIT
!
            if ( rem>0_i4b ) then
!
                m = p*blksz + 1_i4b
!
!$OMP DO SCHEDULE(STATIC)
!
                do j = 1_i4b, n/rem-1_i4b
!
                    i  = (j-1_i4b)*rem + 1_i4b
                    i2 = j*rem
!
                    transpose_lm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                    transpose_lm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
                end do
!
!$OMP END DO NOWAIT
!
!$OMP SINGLE
!
                i  = ( n/rem - 2_i4b )*rem + 1_i4b
                i2 = m - 1_i4b
!
                transpose_lm(m:n,m:n)  = transpose( mat(m:n,m:n) )
                transpose_lm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                transpose_lm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
!$OMP END SINGLE
!
            end if
!
!$OMP END PARALLEL
!
        end if
!
    else
#endif
!
!       BLOCKED TRANSPOSITION OF mat.
!
        if ( n/=m ) then
!
!           RECTANGULAR MATRIX.
!
            do j = 1_i4b, m, blksz
!
                j2 = min( m, j+blksz-1_i4b )
!
                do i = 1_i4b, n, blksz
!
                    i2 = min( n, i+blksz-1_i4b )
!
                    transpose_lm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
!
                end do
!
            end do
!  
        else
!
!           SQUARE MATRIX.
!
            p   = n/blksz
            rem = n - p*blksz
            m   = blksz*( p - 1_i4b ) + 1_i4b
!
            do j = 1_i4b, m, blksz
!
                j2 = j + blksz - 1_i4b
!
                transpose_lm(j:j2,j:j2) = transpose( mat(j:j2,j:j2) )
!
                do i = j+blksz, m, blksz
!
                    i2 = i + blksz - 1_i4b
!
                    transpose_lm(j:j2,i:i2) = transpose( mat(i:i2,j:j2) )
                    transpose_lm(i:i2,j:j2) = transpose( mat(j:j2,i:i2) )
!
                end do
!
            end do
!
            if ( rem>0_i4b ) then
!
                m = p*blksz + 1_i4b
!
                do j = 1_i4b, n/rem-1_i4b
!
                    i  = (j-1_i4b)*rem + 1_i4b
                    i2 = j*rem
!
                    transpose_lm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                    transpose_lm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
                end do
!
                i  = ( n/rem - 2_i4b )*rem + 1_i4b
                i2 = m - 1_i4b
!
                transpose_lm(m:n,m:n)  = transpose( mat(m:n,m:n) )
                transpose_lm(m:n,i:i2) = transpose( mat(i:i2,m:n) )
                transpose_lm(i:i2,m:n) = transpose( mat(m:n,i:i2) )
!
            end if
!
        end if
!
#ifdef _OPENMP3
    end if
#endif
!  
    end function transpose_lm
!
! =========================================================================================
!
    function dot_product_rv( vecx, vecy )
!
! Purpose
! _______
!
!   Forms the dot product of two real vectors
!
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : dot
#else
    use Reals_Constants,  only : zero
#endif
!
    real(stnd), dimension(:), intent(in) :: vecx, vecy
    real(stnd)                           :: dot_product_rv
!
    real(stnd)   :: tmp
    integer(i4b) :: i, k, n
!
    i = size( vecx )
    k = size( vecy )
    n  = min( i, k )
!
#ifdef _BLAS
!
    tmp = dot( n, vecx, 1_i4b, vecy, 1_i4b )
!
#else
!
    tmp = zero
!
    do i = 1_i4b, n
        tmp = tmp + vecx(i)*vecy(i)
    end do
!
#endif
!
    dot_product_rv = tmp
! 
    end function dot_product_rv
!
! =========================================================================================
!
    function dot_product_cv( vecx, vecy )
!
! Purpose
! _______
!
!   Forms the dot product of two complex vectors, conjugating the first vector
!
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : dot
#else
    use Reals_Constants,  only : zero
#endif
!
    complex(stnd), dimension(:), intent(in) :: vecx, vecy
    complex(stnd)                           :: dot_product_cv
!
    complex(stnd) :: tmp
    integer(i4b)  :: i, k, n
!
    i = size( vecx )
    k = size( vecy )
    n  = min( i, k )
!
#ifdef _BLAS
!
    tmp = dot( n, vecx, 1_i4b, vecy, 1_i4b )
!
#else
!
    tmp = cmplx( zero, zero, kind=stnd )
!
    do i = 1_i4b, n
        tmp = tmp + conjg( vecx(i) )*vecy(i)
    end do
!
#endif
!
    dot_product_cv = tmp
! 
    end function dot_product_cv
!
! =========================================================================================
!
    function dot_product_iv( vecx, vecy )
!
! Purpose
! _______
!
!   Forms the dot product of two integer vectors
!
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: vecx, vecy
    integer(i4b)                           :: dot_product_iv
!
    integer(i4b) :: i, k, n, tmp
!
    i = size( vecx )
    k = size( vecy )
    n  = min( i, k )
!
    tmp = 0_i4b
!
    do i = 1_i4b, n
        tmp = tmp + vecx(i)*vecy(i)
    end do
!
    dot_product_iv = tmp
! 
    end function dot_product_iv
!
! =========================================================================================
!
    function dot_product_lv( vecx, vecy )
!
! Purpose
! _______
!
!   Forms the dot product of two logical vectors
! __________________________________________________________________________________________
!
    logical(lgl), dimension(:), intent(in) :: vecx, vecy
    logical(lgl) :: dot_product_lv
!
    integer(i4b) :: nx, ny, n
!
    nx = size( vecx )
    ny = size( vecy )
    n  = min( nx, ny )
!
    dot_product_lv = any( vecx(:n) .and. vecy(:n) )
! 
    end function dot_product_lv
!
! =========================================================================================
!
    function mmproduct_rv_rm( vec, mat )
!
! Purpose
! _______
!
!   Multiplies the real vector VEC by the real matrix MAT
!
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : zero
!
    real(stnd), dimension(:,:), intent(in) :: mat
    real(stnd), dimension(:),   intent(in) :: vec
    real(stnd), dimension(size(mat,2))     :: mmproduct_rv_rm
!
    real(stnd), dimension(size(mat,2)) :: tmpvec
    integer(i4b) :: i, n, m
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        mmproduct_rv_rm(1_i4b:m) = zero
!
        return
!
    end if
!
    do i = 1_i4b, m
        tmpvec(i) = dot_product( vec(1_i4b:n), mat(1_i4b:n,i) )
!       tmpvec(i) = sum( vec(1_i4b:n)*mat(1_i4b:n,i) )
    end do
!
    mmproduct_rv_rm(1_i4b:m) = tmpvec(1_i4b:m)
! 
    end function mmproduct_rv_rm
!
! =========================================================================================
!
    function mmproduct_rm_rv( mat, vec2 )
!
! Purpose
! _______
!
!   Multiplies the real matrix MAT by the real vector VEC2
!
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : zero
!
    real(stnd), dimension(:,:), intent(in) :: mat
    real(stnd), dimension(:),   intent(in) :: vec2
    real(stnd), dimension(size(mat,1))     :: mmproduct_rm_rv
!
    real(stnd), dimension(size(mat,1)) :: tmpvec
    integer(i4b)                       :: i, n, m
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        mmproduct_rm_rv(1_i4b:n) = zero
!
        return
!
    end if
!
    tmpvec(1_i4b:n) = mat(1_i4b:n,1_i4b)*vec2(1_i4b)
!
    do i = 2_i4b, m
        tmpvec(1_i4b:n) = tmpvec(1_i4b:n) + mat(1_i4b:n,i)*vec2(i)
    end do
!
    mmproduct_rm_rv(1_i4b:n) = tmpvec(1_i4b:n)
!
    end function mmproduct_rm_rv
!
! =========================================================================================
!
    function mmproduct_rm_rm( mat1, mat2 )
!
! Purpose
! _______
!
!   Multiplies the real matrix MAT1 by the real matrix MAT2
!
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : zero
!
    real(stnd), dimension(:,:), intent(in)           :: mat1, mat2
    real(stnd), dimension(size(mat1,1),size(mat2,2)) :: mmproduct_rm_rm
!
    real(stnd), dimension(size(mat1,1)) :: tmpvec
    integer(i4b)                        :: i, j, n, m, p
!
    n = size( mat1, 1 )
    m = size( mat2, 2 )
    p = min( size(mat1,2), size(mat2,1) )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( p<1_i4b ) then
!
        mmproduct_rm_rm(1_i4b:n,1_i4b:m) = zero
!
        return
!
    end if
!
    do i = 1_i4b, m
!
        tmpvec(1_i4b:n) = mat1(1_i4b:n,1_i4b)*mat2(1_i4b,i)
!
        do j = 2_i4b, p
            tmpvec(1_i4b:n) = tmpvec(1_i4b:n) + mat1(1_i4b:n,j)*mat2(j,i)
        end do
!
        mmproduct_rm_rm(1_i4b:n,i) = tmpvec(1_i4b:n)
!
    end do
!  
    end function mmproduct_rm_rm
!
! =========================================================================================
!
    function mmproduct_cv_cm( vec, mat )
!
! Purpose
! _______
!
!   Multiplies the complex vector VEC by the complex matrix MAT
!
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : zero
!
    complex(stnd), dimension(:,:), intent(in) :: mat
    complex(stnd), dimension(:),   intent(in) :: vec
    complex(stnd), dimension(size(mat,2))     :: mmproduct_cv_cm
!
    complex(stnd), dimension(size(mat,2)) :: tmpvec
    integer(i4b) :: i, n, m
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        mmproduct_cv_cm(1_i4b:m) = cmplx( zero, zero, kind=stnd )
!
        return
!
    end if
!
    do i = 1_i4b, m
        tmpvec(i) = sum( vec(1_i4b:n)*mat(1_i4b:n,i) )
    end do
!
    mmproduct_cv_cm(1_i4b:m) = tmpvec(1_i4b:m)
! 
    end function mmproduct_cv_cm
!
! =========================================================================================
!
    function mmproduct_cm_cv( mat, vec2 )
!
! Purpose
! _______
!
!   Multiplies the complex matrix MAT by the complex vector VEC2
!
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : zero
!
    complex(stnd), dimension(:,:), intent(in) :: mat
    complex(stnd), dimension(:),   intent(in) :: vec2
    complex(stnd), dimension(size(mat,1))     :: mmproduct_cm_cv
!
    complex(stnd), dimension(size(mat,1)) :: tmpvec
    integer(i4b)                          :: i, n, m
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        mmproduct_cm_cv(1_i4b:m) = cmplx( zero, zero, kind=stnd )
!
        return
!
    end if
!
    tmpvec(1_i4b:n) = mat(1_i4b:n,1_i4b)*vec2(1_i4b)
!
    do i = 2_i4b, m
        tmpvec(1_i4b:n) = tmpvec(1_i4b:n) + mat(1_i4b:n,i)*vec2(i)
    end do
!
    mmproduct_cm_cv(1_i4b:n) = tmpvec(1_i4b:n)
!
    end function mmproduct_cm_cv
!
! =========================================================================================
!
    function mmproduct_cm_cm( mat1, mat2 )
!
! Purpose
! _______
!
!   Multiplies the complex matrix MAT1 by the complex matrix MAT2
!
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : zero
!
    complex(stnd), dimension(:,:), intent(in)           :: mat1, mat2
    complex(stnd), dimension(size(mat1,1),size(mat2,2)) :: mmproduct_cm_cm
!
    complex(stnd), dimension(size(mat1,1)) :: tmpvec
    integer(i4b)                           :: i, j, n, m, p
!
    n = size( mat1, 1 )
    m = size( mat2, 2 )
    p = min( size(mat1,2), size(mat2,1) )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( p<1_i4b ) then
!
        mmproduct_cm_cm(1_i4b:n,1_i4b:m) = cmplx( zero, zero, kind=stnd )
!
        return
!
    end if
!
    do i = 1_i4b, m
!
        tmpvec(1_i4b:n) = mat1(1_i4b:n,1_i4b)*mat2(1_i4b,i)
!
        do j = 2_i4b, p
            tmpvec(1_i4b:n) = tmpvec(1_i4b:n) + mat1(1_i4b:n,j)*mat2(j,i)
        end do
!
         mmproduct_cm_cm(1_i4b:n,i) = tmpvec(1_i4b:n)
!
    end do
!  
    end function mmproduct_cm_cm
!
! =========================================================================================
!
    function matmul_rv_rm( vec, mat )
!
! Purpose
! _______
!
!   Multiplies the real vector VEC by the real matrix MAT
!
!
! Further Details
! _______________
!
!   This function will use the BLAS through the BLAS_interfaces module if the C processor
!   macro _BLAS is activated during compilation. Furthermore, if the _OPENMP3 macro is
!   activated during compilation, this function will be parallelized with OPENMP.
!
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : dot
#endif
    use Reals_Constants,  only : zero
!
    real(stnd), dimension(:,:), intent(in) :: mat
    real(stnd), dimension(:),   intent(in) :: vec
    real(stnd), dimension(size(mat,2))     :: matmul_rv_rm
!
    real(stnd), dimension(size(mat,2)) :: tmpvec
    integer(i4b)                       :: i, n, m
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        matmul_rv_rm(1_i4b:m) = zero
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  m>1_i4b           .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
!$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(i)
!
        do i = 1_i4b, m
#ifdef _BLAS
            tmpvec(i) = dot( n, vec, 1_i4b,  mat(1_i4b:n,i), 1_i4b )
#else
!
!BUG: For compiling STATPACK with the C processor macro _DOT_PRODUCT and some versions of
!     the ifort compiler, we cannot use the intrinsic dot_product function here.
!
#ifdef _DOT_PRODUCT
            tmpvec(i) = sum( vec(1_i4b:n)*mat(1_i4b:n,i) )
#else
            tmpvec(i) = dot_product( vec(1_i4b:n), mat(1_i4b:n,i) )
#endif
!
#endif
        end do
!
!$OMP END PARALLEL DO
!
    else
#endif
!
        do i = 1_i4b, m
#ifdef _BLAS
            tmpvec(i) = dot( n, vec, 1_i4b,  mat(1_i4b:n,i), 1_i4b )
#else
!
!BUG: For compiling STATPACK with the C processor macro _DOT_PRODUCT and some versions of
!     the ifort compiler, we cannot use the intrinsic dot_product function here.
!
#ifdef _DOT_PRODUCT
            tmpvec(i) = sum( vec(1_i4b:n)*mat(1_i4b:n,i) )
#else
            tmpvec(i) = dot_product( vec(1_i4b:n), mat(1_i4b:n,i) )
#endif
!
#endif
        end do
!
#ifdef _OPENMP3
    end if
#endif
!
    matmul_rv_rm(1_i4b:m) = tmpvec(1_i4b:m)
! 
    end function matmul_rv_rm
!
! =========================================================================================
!
    function matmul_rm_rv( mat, vec2 )
!
! Purpose
! _______
!
!   Multiplies the real matrix MAT by the real vector VEC2
!
!
! Further Details
! _______________
!
!   This function will use the BLAS through the BLAS_interfaces module if the C processor
!   macro _BLAS is activated during compilation. Furthermore, if the _OPENMP3 macro is
!   activated during compilation, this function will be parallelized with OPENMP.
!
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : axpy
#endif
    use Reals_Constants,  only : zero
!
    real(stnd), dimension(:,:), intent(in) :: mat
    real(stnd), dimension(:),   intent(in) :: vec2
    real(stnd), dimension(size(mat,1))     :: matmul_rm_rv
!
    real(stnd), dimension(size(mat,1)) :: tmpvec
    integer(i4b)                       :: i, n, m
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        matmul_rm_rv(1_i4b:n) = zero
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  m>1_i4b           .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
        matmul_rm_rv(1_i4b:n) = mat(1_i4b:n,1_i4b)*vec2(1_i4b)
!
!$OMP PARALLEL PRIVATE(i,tmpvec)
!
        tmpvec(1_i4b:n) = zero
!
!$OMP DO SCHEDULE(STATIC)
!
        do i = 2_i4b, m
#ifdef _BLAS
            call axpy( n, vec2(i), mat(1_i4b:n,i), 1_i4b, tmpvec, 1_i4b )
#else
            tmpvec(1_i4b:n) = tmpvec(1_i4b:n) + mat(1_i4b:n,i)*vec2(i)
#endif
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (mulrv)
!
        matmul_rm_rv(1_i4b:n) = matmul_rm_rv(1_i4b:n) + tmpvec(1_i4b:n)
!
!$OMP END CRITICAL (mulrv)
!
!$OMP END PARALLEL
!
    else
#endif
!
        tmpvec(1_i4b:n) = mat(1_i4b:n,1_i4b)*vec2(1_i4b)
!
        do i = 2_i4b, m
#ifdef _BLAS
            call axpy( n, vec2(i), mat(1_i4b:n,i), 1_i4b, tmpvec, 1_i4b )
#else
            tmpvec(1_i4b:n) = tmpvec(1_i4b:n) + mat(1_i4b:n,i)*vec2(i)
#endif
        end do
!
        matmul_rm_rv(1_i4b:n) = tmpvec(1_i4b:n)
!
#ifdef _OPENMP3
    end if
#endif
!
    end function matmul_rm_rv
!
! =========================================================================================
!
    function matmul_rm_rm( mat1, mat2 )
!
! Purpose
! _______
!
!   Multiplies the real matrix MAT1 by the real matrix MAT2
!
!
! Further Details
! _______________
!
!   This function will use the BLAS through the BLAS_interfaces module if the C processor
!   macro _BLAS is activated during compilation. On the other hand, if the _BLAS macro is
!   not activated and the _OPENMP3 macro is activated during compilation, this function
!   will be parallelized with OPENMP if the matrices are big enough.
!
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use Reals_Constants,  only : zero, one
    use BLAS_interfaces,  only : gemm
#else
    use Reals_Constants,  only : zero
#endif
!
    real(stnd), dimension(:,:), intent(in)           :: mat1, mat2
    real(stnd), dimension(size(mat1,1),size(mat2,2)) :: matmul_rm_rm
!
#ifdef _BLAS
    integer(i4b) :: n, m, p, q
#else
!
    integer(i4b) :: n, m, p
#ifdef _OPENMP3
!
    integer(i4b) :: i, i2, j, j2, k, k2, q
!
    logical :: test_par, ompparallel
#endif
#endif
!
    n = size( mat1, 1 )
    m = size( mat2, 2 )
    p = min( size(mat1,2), size(mat2,1) )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( p<1_i4b ) then
!
        matmul_rm_rm(1_i4b:n,1_i4b:m) = zero
!
        return
!
    end if
!
#ifdef _BLAS
!
    q = size( mat2, 1 )
!
    call gemm('N', 'N', n, m, p, one, mat1, n, mat2, q, zero, matmul_rm_rm, n )
!
#else
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel                 .and. &
                  (n+m)*p>=omp_limit               .and. &
                  min(m,n)>max(blksz,i*omp_limit2) .and. &
                  i>1_i4b
!
    if ( test_par ) then
!
        q = min( p, blksz )    
!
!$OMP PARALLEL DO SCHEDULE(STATIC),                    &
!$OMP             PRIVATE(i,i2,j,j2,k,k2), COLLAPSE(2)
!
        do i = 1_i4b, n, blksz
!
            do j = 1_i4b, m, blksz
!
                i2 = min( n, i+blksz-1_i4b )
                j2 = min( m, j+blksz-1_i4b )
!
                matmul_rm_rm(i:i2,j:j2) = matmul( mat1(i:i2,1_i4b:q), mat2(1_i4b:q,j:j2) )
!                matmul_rm_rm(i:i2,j:j2) = mmproduct_rm_rm( mat1(i:i2,1_i4b:q), mat2(1_i4b:q,j:j2) )
!
                do k = blksz+1_i4b, p, blksz
!
                    k2 = min( p, k+blksz-1_i4b )
!
                    matmul_rm_rm(i:i2,j:j2) = matmul_rm_rm(i:i2,j:j2) + matmul( mat1(i:i2,k:k2), mat2(k:k2,j:j2) )
!                    matmul_rm_rm(i:i2,j:j2) = matmul_rm_rm(i:i2,j:j2) + mmproduct_rm_rm( mat1(i:i2,k:k2), mat2(k:k2,j:j2) )
!
                end do
!
            end do
!
        end do
!
!$OMP END PARALLEL DO
!
    else
#endif
!
         matmul_rm_rm(1_i4b:n,1_i4b:m) = matmul( mat1(1_i4b:n,1_i4b:p), mat2(1_i4b:p,1_i4b:m) )
!
#ifdef _OPENMP3
    end if
#endif
!
#endif
!  
    end function matmul_rm_rm
!
! =========================================================================================
!
    function matmul_cv_cm( vec, mat )
!
! Purpose
! _______
!
!   Multiplies the complex vector VEC by the complex matrix MAT
!
!
! Further Details
! _______________
!
!   This function will use the BLAS through the BLAS_interfaces module if the C processor
!   macro _BLAS is activated during compilation. Furthermore, if the _OPENMP3 macro is
!   activated during compilation, this function will be parallelized with OPENMP.
!
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : dotu
#endif
    use Reals_Constants,  only : zero
!
    complex(stnd), dimension(:,:), intent(in) :: mat
    complex(stnd), dimension(:),   intent(in) :: vec
    complex(stnd), dimension(size(mat,2))     :: matmul_cv_cm
!
    complex(stnd), dimension(size(mat,2)) :: tmpvec
    integer(i4b)                          :: i, n, m
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        matmul_cv_cm(1_i4b:m) = cmplx( zero, zero, kind=stnd )
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  m>1_i4b           .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
!$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(i)
!
        do i = 1_i4b, m
#ifdef _BLAS
            tmpvec(i) = dotu( n, vec, 1_i4b,  mat(1_i4b:n,i), 1_i4b )
#else
            tmpvec(i) = sum( vec(1_i4b:n)*mat(1_i4b:n,i) )
#endif
        end do
!
!$OMP END PARALLEL DO
!
    else
#endif
!
        do i = 1_i4b, m
#ifdef _BLAS
            tmpvec(i) = dotu( n, vec, 1_i4b,  mat(1_i4b:n,i), 1_i4b )
#else
            tmpvec(i) = sum( vec(1_i4b:n)*mat(1_i4b:n,i) )
#endif
        end do
!
#ifdef _OPENMP3
    end if
#endif
!
    matmul_cv_cm(1_i4b:m) = tmpvec(1_i4b:m)
! 
    end function matmul_cv_cm
!
! =========================================================================================
!
    function matmul_cm_cv( mat, vec2 )
!
! Purpose
! _______
!
!   Multiplies the complex matrix MAT by the complex vector VEC2
!
!
! Further Details
! _______________
!
!   This function will use the BLAS through the BLAS_interfaces module if the C processor
!   macro _BLAS is activated during compilation. Furthermore, if the _OPENMP3 macro is
!   activated during compilation, this function will be parallelized with OPENMP.
!
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : axpy
#endif
    use Reals_Constants,  only : zero
!
    complex(stnd), dimension(:,:), intent(in) :: mat
    complex(stnd), dimension(:),   intent(in) :: vec2
    complex(stnd), dimension(size(mat,1))     :: matmul_cm_cv
!
    complex(stnd), dimension(size(mat,1)) :: tmpvec
    integer(i4b)                          :: i, n, m
!
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        matmul_cm_cv(1_i4b:m) = cmplx( zero, zero, kind=stnd )
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  m>1_i4b           .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
        matmul_cm_cv(1_i4b:n) = mat(1_i4b:n,1_i4b)*vec2(1_i4b)
!
!$OMP PARALLEL PRIVATE(i,tmpvec)
!
        tmpvec(1_i4b:n) = cmplx( zero, zero, kind=stnd )
!
!$OMP DO SCHEDULE(STATIC)
!
        do i = 2_i4b, m
#ifdef _BLAS
            call axpy( n, vec2(i), mat(1_i4b:n,i), 1_i4b, tmpvec, 1_i4b )
#else
            tmpvec(1_i4b:n) = tmpvec(1_i4b:n) + mat(1_i4b:n,i)*vec2(i)
#endif
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (mulcv)
!
        matmul_cm_cv(1_i4b:n) = matmul_cm_cv(1_i4b:n) + tmpvec(1_i4b:n)
!
!$OMP END CRITICAL (mulcv)
!
!$OMP END PARALLEL
!
    else
#endif
!
        tmpvec(1_i4b:n) = mat(1_i4b:n,1_i4b)*vec2(1_i4b)
!
        do i = 2_i4b, m
#ifdef _BLAS
            call axpy( n, vec2(i), mat(1_i4b:n,i), 1_i4b, tmpvec, 1_i4b )
#else
            tmpvec(1_i4b:n) = tmpvec(1_i4b:n) + mat(1_i4b:n,i)*vec2(i)
#endif
        end do
!
        matmul_cm_cv(1_i4b:n) = tmpvec(1_i4b:n)
!
#ifdef _OPENMP3
    end if
#endif
!
    end function matmul_cm_cv
!
! =========================================================================================
!
    function matmul_cm_cm( mat1, mat2 )
!
! Purpose
! _______
!
!   Multiplies the complex matrix MAT1 by the complex matrix MAT2
!
!
! Further Details
! _______________
!
!   This function will use the BLAS through the BLAS_interfaces module if the C processor
!   macro _BLAS is activated during compilation. On the other hand, if the _BLAS macro is
!   not activated and the _OPENMP3 macro is activated during compilation, this function
!   will be parallelized with OPENMP if the matrices are big enough.
!
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use Reals_Constants,  only : zero, one
    use BLAS_interfaces,  only : gemm
#else
    use Reals_Constants,  only : zero
#endif
!
    complex(stnd), dimension(:,:), intent(in)           :: mat1, mat2
    complex(stnd), dimension(size(mat1,1),size(mat2,2)) :: matmul_cm_cm
!
#ifdef _BLAS
    complex(stnd) :: alpha, beta
!
    integer(i4b)  :: n, m, p, q
#else
!
    integer(i4b)  :: n, m, p
#ifdef _OPENMP3
!
    integer(i4b)  :: i, i2, j, j2, k, k2, q
!
    logical :: test_par, ompparallel
#endif
#endif
!
    n = size( mat1, 1 )
    m = size( mat2, 2 )
    p = min( size(mat1,2), size(mat2,1) )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( p<1_i4b ) then
!
        matmul_cm_cm(1_i4b:n,1_i4b:m) = cmplx( zero, zero, kind=stnd )
!
        return
!
    end if
!
#ifdef _BLAS
!
    q = size( mat2, 1 )
!
    alpha = cmplx( one, zero, kind=stnd )
    beta  = cmplx( zero, zero, kind=stnd )
!
    call gemm('N', 'N', n, m, p, alpha, mat1, n, mat2, q, beta, matmul_cm_cm, n )
!
#else
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel                 .and. &
                  (n+m)*p>=omp_limit               .and. &
                  min(m,n)>max(blksz,i*omp_limit2) .and. &
                  i>1_i4b
!
    if ( test_par ) then
!
        q = min( p, blksz )    
!
!$OMP PARALLEL DO SCHEDULE(STATIC),   &
!$OMP             PRIVATE(i,i2,j,j2,k,k2), COLLAPSE(2)
!
        do i = 1_i4b, n, blksz
!
            do j = 1_i4b, m, blksz
!
                i2 = min( n, i+blksz-1_i4b )
                j2 = min( m, j+blksz-1_i4b )
!
                matmul_cm_cm(i:i2,j:j2) = matmul( mat1(i:i2,1_i4b:q), mat2(1_i4b:q,j:j2) )
!                matmul_cm_cm(i:i2,j:j2) = mmproduct_cm_cm( mat1(i:i2,1_i4b:q), mat2(1_i4b:q,j:j2) )
!
                do k = blksz+1_i4b, p, blksz
!
                    k2 = min( p, k+blksz-1_i4b )
!
                    matmul_cm_cm(i:i2,j:j2) = matmul_cm_cm(i:i2,j:j2) + matmul( mat1(i:i2,k:k2), mat2(k:k2,j:j2) )
!                    matmul_cm_cm(i:i2,j:j2) = matmul_cm_cm(i:i2,j:j2) + mmproduct_cm_cm( mat1(i:i2,k:k2), mat2(k:k2,j:j2) )
!
                end do
!
            end do
!
        end do
!
!$OMP END PARALLEL DO
!
    else
#endif
!
         matmul_cm_cm(1_i4b:n,1_i4b:m) = matmul( mat1(1_i4b:n,1_i4b:p), mat2(1_i4b:p,1_i4b:m) )
!
#ifdef _OPENMP3
    end if
#endif
!
#endif
!  
    end function matmul_cm_cm
!
! =========================================================================================
!
    function matmul_iv_im( vec, mat )
!
! Purpose
! _______
!
!   Multiplies the integer vector VEC by the integer matrix MAT
!
!
! Further Details
! _______________
!
!   If the _OPENMP3 macro is activated during compilation, this function
!   will be parallelized with OPENMP.
!
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:,:), intent(in) :: mat
    integer(i4b), dimension(:),   intent(in) :: vec
    integer(i4b), dimension(size(mat,2))     :: matmul_iv_im
!
    integer(i4b), dimension(size(mat,2)) :: tmpvec
    integer(i4b)                         :: i, n, m
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        matmul_iv_im(1_i4b:m) = 0_i4b
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  m>1_i4b           .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
!$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(i)
!
        do i = 1_i4b, m
!
!BUG: For compiling STATPACK with the C processor macro _DOT_PRODUCT and some versions of
!     the ifort compiler, we cannot use the intrinsic dot_product function here.
!
#ifdef _DOT_PRODUCT
            tmpvec(i) = sum( vec(1_i4b:n)*mat(1_i4b:n,i) )
#else
            tmpvec(i) = dot_product( vec(1_i4b:n), mat(1_i4b:n,i) )
#endif
!
        end do
!
!$OMP END PARALLEL DO
!
    else
#endif
!
        do i = 1_i4b, m
!
!BUG: For compiling STATPACK with the C processor macro _DOT_PRODUCT and some versions of
!     the ifort compiler, we cannot use the intrinsic dot_product function here.
!
#ifdef _DOT_PRODUCT
            tmpvec(i) = sum( vec(1_i4b:n)*mat(1_i4b:n,i) )
#else
            tmpvec(i) = dot_product( vec(1_i4b:n), mat(1_i4b:n,i) )
#endif
!
        end do
!
#ifdef _OPENMP3
    end if
#endif
!
    matmul_iv_im(1_i4b:m) = tmpvec(1_i4b:m)
! 
    end function matmul_iv_im
!
! =========================================================================================
!
    function matmul_im_iv( mat, vec2 )
!
! Purpose
! _______
!
!   Multiplies the integer matrix MAT by the integer vector VEC2
!
!
! Further Details
! _______________
!
!   If the _OPENMP3 macro is activated during compilation, this function
!   will be parallelized with OPENMP.
!
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:,:), intent(in) :: mat
    integer(i4b), dimension(:),   intent(in) :: vec2
    integer(i4b), dimension(size(mat,1))     :: matmul_im_iv
!
    integer(i4b), dimension(size(mat,1)) :: tmpvec
    integer(i4b)                         :: i, n, m
!
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        matmul_im_iv(1_i4b:n) = 0_i4b
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  m>1_i4b           .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
        matmul_im_iv(1_i4b:n) = mat(1_i4b:n,1_i4b)*vec2(1_i4b)
!
!$OMP PARALLEL PRIVATE(i,tmpvec)
!
        tmpvec(1_i4b:n) = 0_i4b
!
!$OMP DO SCHEDULE(STATIC)
!
        do i = 2_i4b, m
            tmpvec(1_i4b:n) = tmpvec(1_i4b:n) + mat(1_i4b:n,i)*vec2(i)
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (muliv)
!
        matmul_im_iv(1_i4b:n) = matmul_im_iv(1_i4b:n) + tmpvec(1_i4b:n)
!
!$OMP END CRITICAL (muliv)
!
!$OMP END PARALLEL
!
    else
#endif
!
        tmpvec(1_i4b:n) = mat(1_i4b:n,1_i4b)*vec2(1_i4b)
!
        do i = 2_i4b, m
            tmpvec(1_i4b:n) = tmpvec(1_i4b:n) + mat(1_i4b:n,i)*vec2(i)
        end do
!
        matmul_im_iv(1_i4b:n) = tmpvec(1_i4b:n)
!
#ifdef _OPENMP3
    end if
#endif
!
    end function matmul_im_iv
!
! =========================================================================================
!
    function matmul_im_im( mat1, mat2 )
!
! Purpose
! _______
!
!   Multiplies the integer matrix MAT1 by the integer matrix MAT2
!
!
! Further Details
! _______________
!
!   If the _OPENMP3 macro is activated during compilation, this function
!   will be parallelized with OPENMP if the matrices are big enough.
!
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:,:), intent(in)           :: mat1, mat2
    integer(i4b), dimension(size(mat1,1),size(mat2,2)) :: matmul_im_im
!
    integer(i4b) :: n, m, p
#ifdef _OPENMP3
!
    integer(i4b) :: i, i2, j, j2, k, k2, q
!
    logical :: test_par, ompparallel
#endif
!
    n = size( mat1, 1 )
    m = size( mat2, 2 )
    p = min( size(mat1,2), size(mat2,1) )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( p<1_i4b ) then
!
        matmul_im_im(1_i4b:n,1_i4b:m) = 0_i4b
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel                 .and. &
                  (n+m)*p>=omp_limit               .and. &
                  min(m,n)>max(blksz,i*omp_limit2) .and. &
                  i>1_i4b
!
    if ( test_par ) then
!
        q = min( p, blksz )    
!
!$OMP PARALLEL DO SCHEDULE(STATIC),    &
!$OMP             PRIVATE(i,i2,j,j2,k,k2), COLLAPSE(2)
!
        do i = 1_i4b, n, blksz
!
            do j = 1_i4b, m, blksz
!
                i2 = min( n, i+blksz-1_i4b )
                j2 = min( m, j+blksz-1_i4b )
!
                matmul_im_im(i:i2,j:j2) = matmul( mat1(i:i2,1_i4b:q), mat2(1_i4b:q,j:j2) )
!
                do k = blksz+1_i4b, p, blksz
!
                    k2 = min( p, k+blksz-1_i4b )
!
                    matmul_im_im(i:i2,j:j2) = matmul_im_im(i:i2,j:j2) + matmul( mat1(i:i2,k:k2), mat2(k:k2,j:j2) )
!
                end do
!
            end do
!
        end do
!
!$OMP END PARALLEL DO
!
    else
#endif
!
        matmul_im_im(1_i4b:n,1_i4b:m) = matmul( mat1(1_i4b:n,1_i4b:p), mat2(1_i4b:p,1_i4b:m) )
!
#ifdef _OPENMP3
    end if
#endif
!  
    end function matmul_im_im
!
! =========================================================================================
!
    function matmul_lv_lm( vec, mat )
!
! Purpose
! _______
!
!   Multiplies the logical vector VEC by the logical matrix MAT
!
!
! Further Details
! _______________
!
!   If the _OPENMP3 macro is activated during compilation, this function
!   will be parallelized with OPENMP.
!
!
! __________________________________________________________________________________________
!
    logical(lgl), dimension(:,:), intent(in) :: mat
    logical(lgl), dimension(:),   intent(in) :: vec
    logical(lgl), dimension(size(mat,2))     :: matmul_lv_lm
!
    logical(lgl), dimension(size(mat,2)) :: tmpvec
    integer(i4b)                         :: i, n, m
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        matmul_lv_lm(1_i4b:m) = false
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  m>1_i4b           .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
!$OMP PARALLEL DO SCHEDULE(STATIC), PRIVATE(i)
!
        do i = 1_i4b, m
!
!BUG: For compiling STATPACK with the C processor macro _DOT_PRODUCT and some versions of
!     the ifort compiler, we cannot use the intrinsic dot_product function here.
!
#ifdef _DOT_PRODUCT
            tmpvec(i) = any( vec(1_i4b:n) .and. mat(1_i4b:n,i) )
#else
            tmpvec(i) = dot_product( vec(1_i4b:n), mat(1_i4b:n,i) )
#endif
!
        end do
!
!$OMP END PARALLEL DO
!
    else
#endif
!
        do i = 1_i4b, m
!
!BUG: For compiling STATPACK with the C processor macro _DOT_PRODUCT and some versions of
!     the ifort compiler, we cannot use the intrinsic dot_product function here.
!
#ifdef _DOT_PRODUCT
            tmpvec(i) = any( vec(1_i4b:n) .and. mat(1_i4b:n,i) )
#else
            tmpvec(i) = dot_product( vec(1_i4b:n), mat(1_i4b:n,i) )
#endif
!
        end do
!
#ifdef _OPENMP3
    end if
#endif
!
    matmul_lv_lm(1_i4b:m) = tmpvec(1_i4b:m)
! 
    end function matmul_lv_lm
!
! =========================================================================================
!
    function matmul_lm_lv( mat, vec2 )
!
! Purpose
! _______
!
!   Multiplies the logical matrix MAT by the logical vector VEC2
!
!
! Further Details
! _______________
!
!   If the _OPENMP3 macro is activated during compilation, this function
!   will be parallelized with OPENMP.
!
!
! __________________________________________________________________________________________
!
    logical(lgl), dimension(:,:), intent(in) :: mat
    logical(lgl), dimension(:),   intent(in) :: vec2
    logical(lgl), dimension(size(mat,1))     :: matmul_lm_lv
!
    logical(lgl), dimension(size(mat,1)) :: tmpvec
    integer(i4b)                         :: i, n, m
!
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<1_i4b .or. m<1_i4b ) then
!
        matmul_lm_lv(1_i4b:n) = false
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  m>1_i4b           .and.      &
                  i>1_i4b
!
    if ( test_par ) then
!
        matmul_lm_lv(1_i4b:n) = mat(1_i4b:n,1_i4b) .and. vec2(1_i4b)
!
!$OMP PARALLEL PRIVATE(i,tmpvec)
!
        tmpvec(1_i4b:n) = false
!
!$OMP DO SCHEDULE(STATIC)
!
        do i = 2_i4b, m
            tmpvec(1_i4b:n) = tmpvec(1_i4b:n) .or. ( mat(1_i4b:n,i) .and. vec2(i) )
        end do
!
!$OMP END DO NOWAIT
!
!$OMP CRITICAL (mullv)
!
        matmul_lm_lv(1_i4b:n) = matmul_lm_lv(1_i4b:n) .or. tmpvec(1_i4b:n)
!
!$OMP END CRITICAL (mullv)
!
!$OMP END PARALLEL
!
    else
#endif
!
        matmul_lm_lv(1_i4b:n) = mat(1_i4b:n,1_i4b) .and. vec2(1_i4b)
!
        do i = 2_i4b, m
            tmpvec(1_i4b:n) = tmpvec(1_i4b:n) .or. ( mat(1_i4b:n,i) .and. vec2(i) )
        end do
!
        matmul_lm_lv(1_i4b:n) = tmpvec(1_i4b:n)
!
#ifdef _OPENMP3
    end if
#endif
!
    end function matmul_lm_lv
!
! =========================================================================================
!
    function matmul_lm_lm( mat1, mat2 )
!
! Purpose
! _______
!
!   Multiplies the logical matrix MAT1 by the logical matrix MAT2
!
!
! Further Details
! _______________
!
!   If the _OPENMP3 macro is activated during compilation, this function
!   will be parallelized with OPENMP if the matrices are big enough.
!
!
! __________________________________________________________________________________________
!
    logical(lgl), dimension(:,:), intent(in)           :: mat1, mat2
    logical(lgl), dimension(size(mat1,1),size(mat2,2)) :: matmul_lm_lm
!
    integer(i4b) :: n, m, p
#ifdef _OPENMP3
!
    integer(i4b) :: i, i2, j, j2, k, k2, q
!
    logical :: test_par, ompparallel
#endif
!
    n = size( mat1, 1 )
    m = size( mat2, 2 )
    p = min( size(mat1,2), size(mat2,1) )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( p<1_i4b ) then
!
        matmul_lm_lm(1_i4b:n,1_i4b:m) = false
!
        return
!
    end if
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel                 .and. &
                  (n+m)*p>=omp_limit               .and. &
                  min(m,n)>max(blksz,i*omp_limit2) .and. &
                  i>1_i4b
!
    if ( test_par ) then
!
        q = min( p, blksz )    
!
!$OMP PARALLEL DO SCHEDULE(STATIC),    &
!$OMP             PRIVATE(i,i2,j,j2,k,k2), COLLAPSE(2)
!
        do i = 1_i4b, n, blksz
!
            do j = 1_i4b, m, blksz
!
                i2 = min( n, i+blksz-1_i4b )
                j2 = min( m, j+blksz-1_i4b )
!
                matmul_lm_lm(i:i2,j:j2) = matmul( mat1(i:i2,1_i4b:q), mat2(1_i4b:q,j:j2) )
!
                do k = blksz+1_i4b, p, blksz
!
                    k2 = min( p, k+blksz-1_i4b )
!
                    matmul_lm_lm(i:i2,j:j2) = matmul_lm_lm(i:i2,j:j2) .or. matmul( mat1(i:i2,k:k2), mat2(k:k2,j:j2) )
!
                end do
!
            end do
!
        end do
!
!$OMP END PARALLEL DO
!
    else
#endif
!
        matmul_lm_lm(1_i4b:n,1_i4b:m) = matmul( mat1(1_i4b:n,1_i4b:p), mat2(1_i4b:p,1_i4b:m) )
!
#ifdef _OPENMP3
    end if
#endif
!  
    end function matmul_lm_lm
!
! =========================================================================================
!     ROUTINES THAT COPY ARRAY WHERE SIZE OF SOURCE IS NOT KNOWN IN ADVANCE.
!     ADAPTED AND EXTENDED FROM Numerical Recipes.
! =========================================================================================
!
!
    subroutine array_copy_iv( src, dest, n_copied, n_not_copied )
!
! Purpose
! _______
!
!   Copies to a destination integer array DEST the one-dimensional integer array SRC, 
!   or as much of SRC as will fit in DEST.
!
!   Returns the number of components copied as N_COPIED, and the number of components
!   not copied as N_NOT_COPIED.
!
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in)  :: src
    integer(i4b), dimension(:), intent(out) :: dest
    integer(i4b),               intent(out) :: n_copied, n_not_copied
!
    n_copied     = min( size(src), size(dest) )
    n_not_copied = size( src ) - n_copied
!
    dest(1_i4b:n_copied) = src(1_i4b:n_copied)
!
    end subroutine array_copy_iv
!
! =========================================================================================
!
    subroutine array_copy_rv( src, dest, n_copied, n_not_copied )
!
! Purpose
! _______
!
!   Copies to a destination real array DEST the one-dimensional real array SRC, 
!   or as much of SRC as will fit in DEST.
!
!   Returns the number of components copied as N_COPIED, and the number of components
!   not copied as N_NOT_COPIED.
!
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : copy
#endif
!
    real(stnd), dimension(:), intent(in)  :: src
    real(stnd), dimension(:), intent(out) :: dest
    integer(i4b),             intent(out) :: n_copied, n_not_copied
!
    n_copied     = min( size(src), size(dest) )
    n_not_copied = size( src ) - n_copied
!
#ifdef _BLAS
    call copy( n_copied, src(1_i4b:n_copied), 1_i4b, dest(1_i4b:n_copied), 1_i4b )
#else
    dest(1_i4b:n_copied) = src(1_i4b:n_copied)
#endif
!
    end subroutine array_copy_rv
!
! =========================================================================================
!
    subroutine array_copy_cv( src, dest, n_copied, n_not_copied )
!
! Purpose
! _______
!
!   Copies to a destination complex array DEST the one-dimensional complex array SRC, 
!   or as much of SRC as will fit in DEST.
!
!   Returns the number of components copied as N_COPIED, and the number of components
!   not copied as N_NOT_COPIED.
!
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : copy
#endif
!
    complex(stnd), dimension(:), intent(in)  :: src
    complex(stnd), dimension(:), intent(out) :: dest
    integer(i4b),                intent(out) :: n_copied, n_not_copied
!
    n_copied     = min( size(src), size(dest) )
    n_not_copied = size( src ) - n_copied
!
#ifdef _BLAS
    call copy( n_copied, src(1_i4b:n_copied), 1_i4b, dest(1_i4b:n_copied), 1_i4b )
#else
    dest(1_i4b:n_copied) = src(1_i4b:n_copied)
#endif
!
    end subroutine array_copy_cv
!
!
! =========================================================================================
!                       ROUTINES THAT SWAP THE CONTENTS OF a AND b.
!                       ADAPTED AND EXTENDED FROM Numerical Recipes.
! =========================================================================================
!
!
    subroutine swap_i( a, b )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the integers A and B
!
!
! __________________________________________________________________________________________
!
    integer(i4b), intent(inout) :: a, b
!
    integer(i4b) :: dum
!
    dum = a
    a   = b
    b   = dum
!
    end subroutine swap_i
!
! =========================================================================================
!
    subroutine swap_r( a, b )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the reals A and B
!    
!
! __________________________________________________________________________________________
!
    real(stnd), intent(inout) :: a, b
!
    real(stnd) :: dum
!
    dum = a
    a   = b
    b   = dum
!
    end subroutine swap_r
!
! =========================================================================================
!
    subroutine swap_c( a, b )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the complex A and B
!    
!
! __________________________________________________________________________________________
!
    complex(stnd), intent(inout) :: a, b
!
    complex(stnd) :: dum
!
    dum = a
    a   = b
    b   = dum
!
    end subroutine swap_c
!
! =========================================================================================
!
    subroutine swap_iv( a, b )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the one-dimensional integer arrays A and B
!    
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(inout) :: a, b
!
    integer(i4b)                     :: n, na, nb
    integer(i4b), dimension(size(a)) :: dum
!
    na = size( a )
    nb = size( b )
    n  = min( na, nb )
!
    dum(:n) = a(:n)
    a(:n)   = b(:n)
    b(:n)   = dum(:n)
!
    end subroutine swap_iv
!
! =========================================================================================
!
    subroutine swap_rv( a, b )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the one-dimensional real arrays A and B
!    
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : swap_blas=>swap
#endif
    real(stnd), dimension(:), intent(inout) :: a, b
!
    integer(i4b) :: n, na, nb
#ifndef _BLAS
    real(stnd), dimension(size(a)) :: dum
#endif
!
    na = size( a )
    nb = size( b )
    n  = min( na, nb )
!
#ifdef _BLAS
    call swap_blas( n, a(:n), 1_i4b, b(:n), 1_i4b )
#else
    dum(:n) = a(:n)
    a(:n)   = b(:n)
    b(:n)   = dum(:n)
#endif
!
    end subroutine swap_rv
!
! =========================================================================================
!
    subroutine swap_cv( a, b )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the one-dimensional complex arrays A and B
!    
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : swap_blas=>swap
#endif
!
    complex(stnd), dimension(:), intent(inout) :: a, b
!
    integer(i4b) :: n, na, nb
#ifndef _BLAS
    complex(stnd), dimension(size(a)) :: dum
#endif
!
    na = size( a )
    nb = size( b )
    n  = min( na, nb )
!
#ifdef _BLAS
    call swap_blas( n, a(:n), 1_i4b, b(:n), 1_i4b )
#else
    dum(:n) = a(:n)
    a(:n)   = b(:n)
    b(:n)   = dum(:n)
#endif
!
    end subroutine swap_cv
!
! =========================================================================================
!
    subroutine swap_im( a, b )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the two-dimensional integer arrays A and B
!    
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:,:), intent(inout) :: a, b
!
    integer(i4b) :: n, na, nb, m, i
    integer(i4b), dimension(size(a,1)) :: dum
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    na = size( a, 1 )
    nb = size( b, 1 )
    n  = min( na, nb )
!
    na = size( a, 2 )
    nb = size( b, 2 )
    m  = min( na, nb )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n*m>=omp_limit     .and.      &
                  m>1_i4b            .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC),    &
!$OMP             PRIVATE(i,dum)
#endif
!
    do i = 1_i4b, m
!
        dum(:n) = a(:n,i)
        a(:n,i) = b(:n,i)
        b(:n,i) = dum(:n)
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine swap_im
!
! =========================================================================================
!
    subroutine swap_rm( a, b )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the two-dimensional real arrays A and B
!    
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : swap_blas=>swap
#endif
!
    real(stnd), dimension(:,:), intent(inout) :: a, b
!
    integer(i4b)  :: n, na, nb, m, i
#ifndef _BLAS
    real(stnd), dimension(size(a,1)) :: dum
#endif
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    na = size( a, 1 )
    nb = size( b, 1 )
    n  = min( na, nb )
!
    na = size( a, 2 )
    nb = size( b, 2 )
    m  = min( na, nb )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n*m>=omp_limit     .and.      &
                  m>1_i4b            .and.      &
                  i>1_i4b
!
#ifdef _BLAS
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC),    &
!$OMP             PRIVATE(i)
#else
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC),    &
!$OMP             PRIVATE(i,dum)
#endif
!
#endif
!
    do i = 1_i4b, m
!
#ifdef _BLAS
        call swap_blas( n, a(:n,i), 1_i4b, b(:n,i), 1_i4b )
#else
        dum(:n) = a(:n,i)
        a(:n,i) = b(:n,i)
        b(:n,i) = dum(:n)
#endif
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine swap_rm
!
! =========================================================================================
!
    subroutine swap_cm( a, b )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the two-dimensional complex arrays A and B
!    
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : swap_blas=>swap
#endif
!
    complex(stnd), dimension(:,:), intent(inout) :: a, b
!
    integer(i4b)  :: n, na, nb, m, i
#ifndef _BLAS
    complex(stnd), dimension(size(a,1)) :: dum
#endif
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    na = size( a, 1 )
    nb = size( b, 1 )
    n  = min( na, nb )
!
    na = size( a, 2 )
    nb = size( b, 2 )
    m  = min( na, nb )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n*m>=omp_limit     .and.      &
                  m>1_i4b            .and.      &
                  i>1_i4b
!
#ifdef _BLAS
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC),    &
!$OMP             PRIVATE(i)
#else
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC),    &
!$OMP             PRIVATE(i,dum)
#endif
!
#endif
!
    do i = 1_i4b, m
!
#ifdef _BLAS
        call swap_blas( n, a(:n,i), 1_i4b, b(:n,i), 1_i4b )
#else
        dum(:n) = a(:n,i)
        a(:n,i) = b(:n,i)
        b(:n,i) = dum(:n)
#endif
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine swap_cm
!
! =========================================================================================
!
    subroutine masked_swap_i( a, b, mask )
!
! Purpose
! _______
!
!   Swap the integers A and B if MASK=true
!    
!
! __________________________________________________________________________________________
!
    integer(i4b), intent(inout) :: a, b
    logical(lgl), intent(in)    :: mask
!
    integer(i4b) :: swp
!
    if ( mask ) then
        swp = a
        a   = b
        b   = swp
    end if
!
    end subroutine masked_swap_i
!
! =========================================================================================
!
     subroutine masked_swap_r( a, b, mask )
!
! Purpose
! _______
!
!   Swap the reals A and B if MASK=true
!    
!
! __________________________________________________________________________________________
!
    real(stnd),   intent(inout) :: a, b
    logical(lgl), intent(in)    :: mask
!
    real(stnd) :: swp
!
    if ( mask ) then
        swp = a
        a   = b
        b   = swp
    end if
!
    end subroutine masked_swap_r
!
! =========================================================================================
!
    subroutine masked_swap_c( a, b, mask )
!
! Purpose
! _______
!
!   Swap the complex A and B if MASK=true
!    
!
! __________________________________________________________________________________________
!
    complex(stnd), intent(inout) :: a, b
    logical(lgl),   intent(in)   :: mask
!
    complex(stnd) :: swp
!
    if ( mask ) then
        swp = a
        a   = b
        b   = swp
    end if
!
    end subroutine masked_swap_c
!
! =========================================================================================
!
    subroutine masked_swap_iv( a, b, mask )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the one-dimensional integer arrays A and B,
!   if the corresponding element of the one-dimensional logical array MASK is true.
!    
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(inout) :: a, b
    logical(lgl), dimension(:), intent(in)    :: mask
!
    integer(i4b)                     :: n, na, nb, nmask
    integer(i4b), dimension(size(a)) :: swp
!
    na    = size( a )
    nb    = size( b )
    nmask = size( mask )
    n     = min( na, nb, nmask )
!
    where ( mask(:n) )
        swp(:n) = a(:n)
        a(:n)   = b(:n)
        b(:n)   = swp(:n)
    end where
!
    end subroutine masked_swap_iv
!
! =========================================================================================
!
    subroutine masked_swap_rv( a, b, mask )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the one-dimensional real arrays A and B,
!   if the corresponding element of the one-dimensional logical array MASK is true.
!    
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:),   intent(inout) :: a, b
    logical(lgl), dimension(:), intent(in)    :: mask
!
    integer(i4b)                   :: n, na, nb, nmask
    real(stnd), dimension(size(a)) :: swp
!
    na    = size( a )
    nb    = size( b )
    nmask = size( mask )
    n     = min( na, nb, nmask )
!
    where ( mask(:n) )
        swp(:n) = a(:n)
        a(:n)   = b(:n)
        b(:n)   = swp(:n)
    end where
!
    end subroutine masked_swap_rv
!
! =========================================================================================
!
    subroutine masked_swap_cv( a, b, mask )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the one-dimensional complex arrays A and B,
!   if the corresponding element of the one-dimensional logical array MASK is true.
!    
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:), intent(inout) :: a, b
    logical(lgl),   dimension(:), intent(in)   :: mask
!
    integer(i4b)                      :: n, na, nb, nmask
    complex(stnd), dimension(size(a)) :: swp
!
    na    = size( a )
    nb    = size( b )
    nmask = size( mask )
    n     = min( na, nb, nmask )
!
    where ( mask(:n) )
        swp(:n) = a(:n)
        a(:n)   = b(:n)
        b(:n)   = swp(:n)
    end where
!
    end subroutine masked_swap_cv
!
! =========================================================================================
!
    subroutine masked_swap_im( a, b, mask )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the two-dimensional integer arrays A and B,
!   if the corresponding element of the two-dimensional logical array MASK is true.
!    
!
! __________________________________________________________________________________________
!
    integer(i4b),   dimension(:,:), intent(inout) :: a, b
    logical(lgl), dimension(:,:),   intent(in)    :: mask
!
    integer(i4b) :: n, na, nb, nmask, m, i
    integer(i4b), dimension(size(a,1)) :: swp
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    na    = size( a, 1 )
    nb    = size( b, 1 )
    nmask = size( mask, 1 )
    n     = min( na, nb, nmask )
!
    na    = size( a, 2 )
    nb    = size( b, 2 )
    nmask = size( mask, 2 )
    m     = min( na, nb, nmask )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n*m>=omp_limit     .and.      &
                  m>1_i4b            .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC),    &
!$OMP             PRIVATE(i,swp)
#endif
!
    do i = 1_i4b, m
!
        where ( mask(:n,i) )
!
            swp(:n) = a(:n,i)
            a(:n,i) = b(:n,i)
            b(:n,i) = swp(:n)
!
        end where
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine masked_swap_im
!
! =========================================================================================
!
    subroutine masked_swap_rm( a, b, mask )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the two-dimensional real arrays A and B,
!   if the corresponding element of the two-dimensional logical array MASK is true.
!    
!
! __________________________________________________________________________________________
!
    real(stnd),   dimension(:,:), intent(inout) :: a, b
    logical(lgl), dimension(:,:), intent(in)    :: mask
!
    integer(i4b) :: n, na, nb, nmask, m, i
    real(stnd), dimension(size(a,1)) :: swp
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    na    = size( a, 1 )
    nb    = size( b, 1 )
    nmask = size( mask, 1 )
    n     = min( na, nb, nmask )
!
    na    = size( a, 2 )
    nb    = size( b, 2 )
    nmask = size( mask, 2 )
    m     = min( na, nb, nmask )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n*m>=omp_limit     .and.      &
                  m>1_i4b            .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC),    &
!$OMP             PRIVATE(i,swp)
#endif
!
    do i = 1_i4b, m
!
        where ( mask(:n,i) )
!
            swp(:n) = a(:n,i)
            a(:n,i) = b(:n,i)
            b(:n,i) = swp(:n)
!
        end where
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine masked_swap_rm
!
! =========================================================================================
!
    subroutine masked_swap_cm( a, b, mask )
!
! Purpose
! _______
!
!   Swap the corresponding elements of the two-dimensional complex arrays A and B,
!   if the corresponding element of the two-dimensional logical array MASK is true.
!    
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:,:), intent(inout) :: a, b
    logical(lgl),   dimension(:,:), intent(in)   :: mask
!
    integer(i4b) :: n, na, nb, nmask, m, i
    complex(stnd), dimension(size(a,1)) :: swp
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    na    = size( a, 1 )
    nb    = size( b, 1 )
    nmask = size( mask, 1 )
    n     = min( na, nb, nmask )
!
    na    = size( a, 2 )
    nb    = size( b, 2 )
    nmask = size( mask, 2 )
    m     = min( na, nb, nmask )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n*m>=omp_limit     .and.      &
                  m>1_i4b            .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC),    &
!$OMP             PRIVATE(i,swp)
#endif
!
    do i = 1_i4b, m
!
        where ( mask(:n,i) )
!
            swp(:n) = a(:n,i)
            a(:n,i) = b(:n,i)
            b(:n,i) = swp(:n)
!
        end where
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine masked_swap_cm
!
!
! =========================================================================================
!            ROUTINES THAT REALLOCATE AN ALLOCATABLE ARRAY WHILE PRESERVING ITS CONTENT.
!                       REQUIRED FORTRAN2003.
! =========================================================================================
!
!
#ifdef _F2003
    subroutine mvalloc_iv( p, n, ialloc )
!
! Purpose
! _______
!
!   Reallocates an allocatable array P to an integer one dimensional array with a new size N,
!   while preserving its contents.
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), allocatable, intent(inout) :: p
    integer(i4b), intent(in)  :: n
    integer,      intent(out) :: ialloc
!
    integer(i4b), dimension(:), allocatable :: temp
    integer(i4b) :: nold, i
!
    if ( .not.allocated( p ) ) then
!
        allocate( p(n), stat=ialloc )
!
    else
!
        allocate( temp(n), stat=ialloc )
!     
        if ( ialloc==0 )  then
!
            nold = size( p )
            i = min( nold, n )
            temp(1_i4b:i) = p(1_i4b:i)
            call move_alloc( temp, p )
!
        end if
!
    end if
!
    end subroutine mvalloc_iv
!
! =========================================================================================
!
    subroutine mvalloc_rv( p, n, ialloc )
!
! Purpose
! _______
!
!   Reallocates an allocatable array P to a real one dimensional array with a new size N,
!   while preserving its contents.
!
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:), allocatable, intent(inout) :: p
    integer(i4b), intent(in)  :: n
    integer,      intent(out) :: ialloc
!
    real(stnd), dimension(:), allocatable :: temp
    integer(i4b) :: nold, i
!
    if ( .not.allocated( p ) ) then
!
        allocate( p(n), stat=ialloc )
!
    else
!
        allocate( temp(n), stat=ialloc )
!     
        if ( ialloc==0 )  then
!
            nold = size( p )
            i = min( nold, n )
            temp(1_i4b:i) = p(1_i4b:i)
            call move_alloc( temp, p )
!
        end if
!
    end if
!
    end subroutine mvalloc_rv
!
! =========================================================================================
!
    subroutine mvalloc_cv( p, n, ialloc )
!
! Purpose
! _______
!
!   Reallocates an allocatable array P to a complex one dimensional array with a new size N,
!   while preserving its contents.
!
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:), allocatable, intent(inout) :: p
    integer(i4b), intent(in)  :: n
    integer,      intent(out) :: ialloc
!
    complex(stnd), dimension(:), allocatable :: temp
    integer(i4b) :: nold, i
!
    if ( .not.allocated( p ) ) then
!
        allocate( p(n), stat=ialloc )
!
    else
!
        allocate( temp(n), stat=ialloc )
!     
        if ( ialloc==0 )  then
!
            nold = size( p )
            i = min( nold, n )
            temp(1_i4b:i) = p(1_i4b:i)
            call move_alloc( temp, p )
!
        end if
!
    end if
!
    end subroutine mvalloc_cv
!
! =========================================================================================
!
    subroutine mvalloc_hv( p, n, ialloc )
!
! Purpose
! _______
!
!   Reallocates an allocatable array P to a character one dimensional array with a new size N,
!   while preserving its contents.
!
!
! __________________________________________________________________________________________
!
    character(1), dimension(:), allocatable, intent(inout) :: p
    integer(i4b), intent(in)  :: n
    integer,      intent(out) :: ialloc
!
    character(1), dimension(:), allocatable :: temp
    integer(i4b) :: nold, i
!
    if ( .not.allocated( p ) ) then
!
        allocate( p(n), stat=ialloc )
!
    else
!
        allocate( temp(n), stat=ialloc )
!     
        if ( ialloc==0 )  then
!
            nold = size( p )
            i = min( nold, n )
            temp(1_i4b:i) = p(1_i4b:i)
            call move_alloc( temp, p )
!
        end if
!
    end if
!
    end subroutine mvalloc_hv
!
! =========================================================================================
!
    subroutine mvalloc_im( p, n, m, ialloc )
!
! Purpose
! _______
!
!   Reallocates an allocatable array P to an integer two dimensional array with a new shape (N,M)
!   while preserving its contents.
!
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:,:), allocatable, intent(inout) :: p
    integer(i4b), intent(in)  :: n, m
    integer,      intent(out) :: ialloc
!
    integer(i4b), dimension(:,:), allocatable :: temp
    integer(i4b) :: old, i, j
!
    if ( .not.allocated( p ) ) then
!
        allocate( p(n,m), stat=ialloc )
!
    else
!
        allocate( temp(n,m), stat=ialloc )
!     
        if ( ialloc==0 )  then
!
            old = size( p, 1 )
            i = min( old, n )
            old = size( p, 2 )
            j = min( old, m )

            temp(1_i4b:i,1_i4b:j) = p(1_i4b:i,1_i4b:j)
            call move_alloc( temp, p )
!
        end if
!
    end if
!
    end subroutine mvalloc_im
!
! =========================================================================================
!
    subroutine mvalloc_rm( p, n, m, ialloc )
!
! Purpose
! _______
!
!   Reallocates an allocatable array P to a real two dimensional array with a new shape (N,M)
!   while preserving its contents.
!
!
! __________________________________________________________________________________________
!
    real(stnd),   dimension(:,:), allocatable, intent(inout) :: p
    integer(i4b), intent(in)  :: n, m
    integer,      intent(out) :: ialloc
!
    real(stnd), dimension(:,:), allocatable :: temp
    integer(i4b) :: old, i, j
!
    if ( .not.allocated( p ) ) then
!
        allocate( p(n,m), stat=ialloc )
!
    else
!
        allocate( temp(n,m), stat=ialloc )
!     
        if ( ialloc==0 )  then
!
            old = size( p, 1 )
            i = min( old, n )
            old = size( p, 2 )
            j = min( old, m )

            temp(1_i4b:i,1_i4b:j) = p(1_i4b:i,1_i4b:j)
            call move_alloc( temp, p )
!
        end if
!
    end if
!
    end subroutine mvalloc_rm
!
! =========================================================================================
!
    subroutine mvalloc_cm( p, n, m, ialloc )
!
! Purpose
! _______
!
!   Reallocates an allocatable array P to a complex two dimensional array with a new shape (N,M)
!   while preserving its contents.
!
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:,:), allocatable, intent(inout) :: p
    integer(i4b), intent(in)  :: n, m
    integer,      intent(out) :: ialloc
!
    complex(stnd), dimension(:,:), allocatable :: temp
    integer(i4b) :: old, i, j
!
    if ( .not.allocated( p ) ) then
!
        allocate( p(n,m), stat=ialloc )
!
    else
!
        allocate( temp(n,m), stat=ialloc )
!     
        if ( ialloc==0 )  then
!
            old = size( p, 1 )
            i = min( old, n )
            old = size( p, 2 )
            j = min( old, m )

            temp(1_i4b:i,1_i4b:j) = p(1_i4b:i,1_i4b:j)
            call move_alloc( temp, p )
!
        end if
!
    end if
!
    end subroutine mvalloc_cm
#endif
!
!
! =========================================================================================
!                   ROUTINES RETURNING A LOCATION AS AN INTEGER VALUE.
!                   ADAPTED AND EXTENDED FROM Numerical Recipes.
! =========================================================================================
!
!
    function ifirstloc( mask )
!
! Purpose
! _______
!
!   Returns the index of the first location, in a one-dimensional logical MASK,
!   that has the value true, or returns size(MASK)+1 if all components of MASK
!   are false .
!
!
! __________________________________________________________________________________________
!
    logical(lgl), dimension(:), intent(in) :: mask
!
    integer(i4b), dimension(1) :: loc
    integer(i4b) :: ifirstloc
    integer(i4b), dimension(size(mask)) :: tmp
!
    where ( mask )
        tmp = 1_i4b
    elsewhere
        tmp = 0_i4b
    end where
    loc = maxloc( tmp )
    ifirstloc = loc(1_i4b)
    if ( .not.mask(ifirstloc) ) ifirstloc = size( mask ) + 1
!
    end function ifirstloc
!
! =========================================================================================
!
    function imaxloc_i( arr )
!
! Purpose
! _______
!
!   Returns location of the one-dimensional integer array ARR maximum as an integer. 
!
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: arr
!
    integer(i4b), dimension(1) :: imax
    integer(i4b) :: imaxloc_i
!
    if ( size(arr)==0 ) then
        imaxloc_i = 1_i4b
    else
        imax= maxloc( arr )
        imaxloc_i = imax(1_i4b)
    end if
!
    end function imaxloc_i
!
! =========================================================================================
!
    function masked_imaxloc_i( arr, mask )
!
! Purpose
! _______
!
!   Returns location as an integer of the maximum in the elements of the one-dimensional 
!   integer array ARR under the control of the one-dimensional logical array MASK.
!   Returns size(MASK)+1 if all components of MASK are false .
!    
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: arr
    logical(lgl), dimension(:), intent(in) :: mask
!
    integer(i4b), dimension(1) :: imax
    integer(i4b) :: masked_imaxloc_i
!
    if ( any(mask) ) then
        imax = maxloc( arr, mask )
        masked_imaxloc_i = imax(1_i4b)
    else
        masked_imaxloc_i = size( arr ) + 1_i4b
    end if
!
    end function masked_imaxloc_i
!
! =========================================================================================
!
    function imaxloc_r( arr )
!
! Purpose
! _______
!
!   Returns location of the one-dimensional real array ARR maximum as an integer. 
!    
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:), intent(in) :: arr
!
    integer(i4b), dimension(1) :: imax
    integer(i4b) :: imaxloc_r
!
    if ( size(arr)==0 ) then
        imaxloc_r = 1_i4b
    else
        imax = maxloc( arr )
        imaxloc_r = imax(1_i4b)
    end if
!
    end function imaxloc_r
!
! =========================================================================================
!
    function masked_imaxloc_r( arr, mask )
!
! Purpose
! _______
!
!   Returns location as an integer of the maximum in the elements of the one-dimensional 
!   real array ARR under the control of the one-dimensional logical array MASK.
!   Returns size(MASK)+1 if all components of MASK are false .
!    
!
! __________________________________________________________________________________________
!
    real(stnd),   dimension(:), intent(in) :: arr
    logical(lgl), dimension(:), intent(in) :: mask
!
    integer(i4b), dimension(1) :: imax
    integer(i4b) :: masked_imaxloc_r
!
    if ( any(mask) ) then
        imax = maxloc( arr, mask )
        masked_imaxloc_r = imax(1_i4b)
    else
        masked_imaxloc_r = size( arr ) + 1_i4b
    end if
!
    end function masked_imaxloc_r
!
! =========================================================================================
!
    function iminloc_i( arr )
!
! Purpose
! _______
!
!   Returns location of the one-dimensional integer array ARR minimum as an integer. 
!    
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: arr
!
    integer(i4b), dimension(1) :: imin
    integer(i4b) :: iminloc_i
!
    if ( size(arr)==0 ) then
        iminloc_i = 1_i4b
    else
        imin = minloc( arr )
        iminloc_i = imin(1_i4b)
    end if
!
    end function iminloc_i
!
! =========================================================================================
!
    function masked_iminloc_i( arr, mask )
!
! Purpose
! _______
!
!   Returns location as an integer of the minimum in the elements of the one-dimensional 
!   integer array ARR under the control of the one-dimensional logical array MASK.
!   Returns size(MASK)+1 if all components of MASK are false .
!    
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: arr
    logical(lgl), dimension(:), intent(in) :: mask
!
    integer(i4b), dimension(1) :: imin
    integer(i4b) :: masked_iminloc_i
!
    if ( any(mask) ) then
        imin = minloc( arr, mask )
        masked_iminloc_i = imin(1_i4b)
    else
        masked_iminloc_i = size( arr ) + 1_i4b
    end if
!
    end function masked_iminloc_i
!
! =========================================================================================
!
    function iminloc_r( arr )
!
! Purpose
! _______
!
!   Returns location of the one-dimensional real array ARR minimum as an integer. 
!    
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:), intent(in) :: arr
!
    integer(i4b), dimension(1) :: imin
    integer(i4b) :: iminloc_r
!
    if ( size(arr)==0 ) then
       iminloc_r = 1_i4b
    else
       imin = minloc( arr )
       iminloc_r = imin(1_i4b)
    end if
!
    end function iminloc_r
!
! =========================================================================================
!
    function masked_iminloc_r( arr, mask )
!
! Purpose
! _______
!
!   Returns location as an integer of the minimum in the elements of the one-dimensional 
!   real array ARR under the control of the one-dimensional logical array MASK.
!   Returns size(MASK)+1 if all components of MASK are false .
!    
!
! __________________________________________________________________________________________
!
    real(stnd),   dimension(:), intent(in) :: arr
    logical(lgl), dimension(:), intent(in) :: mask
!
    integer(i4b), dimension(1) :: imin
    integer(i4b) :: masked_iminloc_r
!
    if ( any(mask) ) then
        imin = minloc( arr, mask )
        masked_iminloc_r = imin(1_i4b)
    else
        masked_iminloc_r = size( arr ) + 1_i4b
    end if
!
    end function masked_iminloc_r
!
!
! =========================================================================================
!                 ROUTINES FOR ARGUMENT CHECKING AND ERROR HANDLING.
!                 ADAPTED AND EXTENDED FROM Numerical Recipes.
! =========================================================================================
!
!
    subroutine assert1( n1, string )
!
! Purpose
! _______
!
!   Exit with error message STRING, if logical argument n1 is false .
!
!
! __________________________________________________________________________________________
!
    use Char_Constants,   only : utilities_error1
!
    character(len=*), intent(in) :: string
    logical(lgl),     intent(in) :: n1
!
    if ( .not.n1 )                      &
    call merror( utilities_error1//string )
!
    end subroutine assert1
!
! =========================================================================================
!
    subroutine assert2( n1, n2, string )
!
! Purpose
! _______
!
!   Exit with error message STRING, if any of the logical arguments n1, n2 are false .
!    
!
! __________________________________________________________________________________________
!
    use Char_Constants,   only : utilities_error1
!
    character(len=*), intent(in) :: string
    logical(lgl),     intent(in) :: n1, n2
!
    if ( .not.( n1 .and. n2 ) )           &
    call merror( utilities_error1//string )
!
    end subroutine assert2
!
! =========================================================================================
!
    subroutine assert3( n1, n2, n3, string )
!
! Purpose
! _______
!
!   Exit with error message STRING, if any of the logical arguments n1, n2, n3 are false .
!    
!
! __________________________________________________________________________________________
!
    use Char_Constants,   only : utilities_error1
!
    character(len=*), intent(in) :: string
    logical(lgl),     intent(in) :: n1, n2, n3
!
    if ( .not.( n1 .and. n2 .and. n3 ) )   &
    call merror( utilities_error1//string )
!
    end subroutine assert3
!
! =========================================================================================
!
    subroutine assert4( n1, n2, n3, n4, string )
!
! Purpose
! _______
!
!   Exit with error message STRING, if any of the logical arguments n1, n2, n3, n4
!   are false .
!    
!
! __________________________________________________________________________________________
!
    use Char_Constants,   only : utilities_error1
!
    character(len=*), intent(in) :: string
    logical(lgl),     intent(in) :: n1, n2, n3, n4
!
    if ( .not.( n1 .and. n2 .and. n3 .and. n4 ) )    &
    call merror( utilities_error1//string )
!
    end subroutine assert4
!
! =========================================================================================
!
    subroutine assert_v( n, string )
!
! Purpose
! _______
!
!   Exit with error message STRING, if any of the elements of the one-dimensional logical 
!   array N are false .
!    
!
! __________________________________________________________________________________________
!
    use Char_Constants,   only : utilities_error1
!
    character(len=*), intent(in) :: string
    logical(lgl),     intent(in) :: n(:)
!
    if ( .not.all( n ) )          &
    call merror( utilities_error1//string )
!
    end subroutine assert_v
!
! =========================================================================================
!
    function assert_eq2( n1, n2, string )
!
! Purpose
! _______
!
!   Exit with error message STRING, if the integer arguments n1, n2 are not equal.
!
!
! __________________________________________________________________________________________
!
    use Char_Constants,   only : utilities_error2
!
    character(len=*), intent(in) :: string
    integer(i4b),     intent(in) :: n1, n2
!
    integer(i4b) :: assert_eq2
!
    if ( n1==n2 ) then
        assert_eq2 = n1
    else
        call merror( utilities_error2//string )
    end if
!
    end function assert_eq2
!
! =========================================================================================
!
    function assert_eq3( n1, n2, n3, string )
!
! Purpose
! _______
!
!   Exit with error message STRING, if the integer arguments n1, n2, n3 
!   are not all equal.
!    
!
! __________________________________________________________________________________________
!
    use Char_Constants,   only : utilities_error2
!
    character(len=*), intent(in) :: string
    integer(i4b),     intent(in) :: n1, n2, n3
!
    integer(i4b) :: assert_eq3
!
    if ( n1==n2 .and. n2==n3 ) then
        assert_eq3 = n1
    else
        call merror( utilities_error2//string )
    end if
!
    end function assert_eq3
!
! =========================================================================================
!
    function assert_eq4( n1, n2, n3, n4, string )
!
! Purpose
! _______
!
!   Exit with error message STRING, if the integer arguments n1, n2, n3, n4
!   are not all equal.
!   
!
! __________________________________________________________________________________________
!
    use Char_Constants,   only : utilities_error2
!
    character(len=*), intent(in) :: string
    integer(i4b),     intent(in) :: n1, n2, n3, n4
!
    integer(i4b) :: assert_eq4
!
    if ( n1==n2 .and. n2==n3 .and. n3==n4 ) then
        assert_eq4 = n1
    else
        call merror( utilities_error2//string )
    end if
!
    end function assert_eq4
!
! =========================================================================================
!
    function assert_eqn( nn, string )
!
! Purpose
! _______
!
!   Exit with error message STRING, if the elements of the one-dimensional integer array NN
!   are not all equal.
!   
!
! __________________________________________________________________________________________
!
    use Char_Constants,   only : utilities_error2
!
    character(len=*), intent(in) :: string
    integer(i4b),     intent(in) :: nn(:)
!
    integer(i4b) :: assert_eqn
!
    if ( all( nn(2:)==nn(1) ) ) then
        assert_eqn = nn(1)
    else
        call merror( utilities_error2//string )
    end if
!
    end function assert_eqn
!
! =========================================================================================
!
    subroutine merror( string, ierror )
!
! Purpose
! _______
!
!   Report error message STRING and optional error number IERROR and stop.
!
!
! __________________________________________________________________________________________
!
    use Char_Constants,   only : utilities_error3
!
    character(len=*),       intent(in) :: string
    integer(i4b), optional, intent(in) :: ierror
!
    if ( present(ierror) ) then
         write( defunit, * ) 'MERROR: ', string, ', ERROR= ',ierror
    else
         write( defunit, * ) 'MERROR: ',string
    end if
!
    write( defunit, * ) utilities_error3
!
    stop
!
    end subroutine merror
!
!
! =========================================================================================
!                   ARRAY FUNCTIONS RETURNING AN ARITHMETIC PROGRESSION.
!                   ADAPTED AND EXTENDED FROM Numerical Recipes.
! =========================================================================================
!
!
    function arth_i( first, increment, n )
!
! Purpose
! _______
!
!   Returns an one-dimensional integer array of length N containing an arithmetic 
!   progression whose first value is FIRST and whose increment is INCREMENT.
!
!
! __________________________________________________________________________________________
!
    integer(i4b), intent(in) :: first, increment, n
!
    integer(i4b), dimension(n) :: arth_i
    integer(i4b) :: k, k2, temp
!
    if ( n>0_i4b ) arth_i(1_i4b) = first
!
    if ( n<=npar_arth ) then
!
        do k = 2_i4b, n
            arth_i(k) = arth_i(k-1_i4b) + increment
        end do
!
    else
!
        do k = 2_i4b, npar2_arth
            arth_i(k) = arth_i(k-1_i4b) + increment
        end do
!
        temp = increment*npar2_arth
        k = npar2_arth
!
        do
!
            if ( k>=n ) exit
!
            k2 = k + k
            arth_i(k+1_i4b:min(k2,n)) = temp + arth_i(1_i4b:min(k,n-k))
            temp = temp + temp
            k = k2
!
        end do
!
    end if
!
    end function arth_i
!
! =========================================================================================
!
    function arth_r( first, increment, n )
!
! Purpose
! _______
!
!   Returns an one-dimensional real array of length N containing an arithmetic 
!   progression whose first value is FIRST and whose increment is INCREMENT.
!   
!
! __________________________________________________________________________________________
!
    real(stnd),   intent(in) :: first, increment
    integer(i4b), intent(in) :: n
!
    real(stnd), dimension(n) :: arth_r
    real(stnd) :: temp
    integer(i4b) :: k, k2
!
    if ( n>0_i4b ) arth_r(1_i4b) = first
!
    if ( n<=npar_arth ) then
!
        do k = 2_i4b, n
            arth_r(k) = arth_r(k-1_i4b) + increment
        end do
!
    else
!
        do k = 2_i4b, npar2_arth
            arth_r(k) = arth_r(k-1_i4b) + increment
        end do
!
        temp = increment*npar2_arth
        k = npar2_arth
!
        do
!
            if ( k>=n ) exit
!
            k2= k + k
            arth_r(k+1_i4b:min(k2,n)) = temp + arth_r(1_i4b:min(k,n-k))
            temp = temp + temp
            k = k2
!
        end do
!
    end if
!
    end function arth_r
!
! =========================================================================================
!
    function arth_c( first, increment, n )
!
! Purpose
! _______
!
!   Returns an one-dimensional complex array of length N containing an arithmetic 
!   progression whose first value is FIRST and whose increment is INCREMENT.
!   
!
! __________________________________________________________________________________________
!
    complex(stnd),  intent(in) :: first, increment
    integer(i4b),    intent(in) :: n
!
    complex(stnd), dimension(n) :: arth_c
    complex(stnd)               :: temp
    integer(i4b)                 :: k, k2
!
    if ( n>0_i4b ) arth_c(1_i4b) = first
!
    if ( n<=npar_arth ) then
!
        do k = 2_i4b, n
            arth_c(k) = arth_c(k-1_i4b) + increment
        end do
!
    else
!
        do k = 2_i4b, npar2_arth
            arth_c(k) = arth_c(k-1_i4b) + increment
        end do
!
        temp = increment*npar2_arth
        k = npar2_arth
!
        do
!
            if ( k>=n ) exit
!
            k2= k + k
            arth_c(k+1_i4b:min(k2,n)) = temp + arth_c(1_i4b:min(k,n-k))
            temp = temp + temp
            k = k2
!
        end do
!
    end if
!
    end function arth_c
!
! =========================================================================================
!
    function arth_iv( first, increment, n )
!
! Purpose
! _______
!
!   Returns a two-dimensional integer array containing size(FIRST) = size(INCREMENT) 
!   arithmetic progressions of length N whose first values are FIRST(:) and whose 
!   increments are INCREMENT(:).
!
!   It is assumed that the vector arguments FIRST and INCREMENT have the same length.
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: first, increment
    integer(i4b),               intent(in) :: n
!
    integer(i4b), dimension(size(first),n) :: arth_iv
    integer(i4b)                           :: l, k, k2, temp(size(first))
!
    if ( n>0_i4b ) arth_iv(:,1_i4b) = first(:)
!
    if ( n<=npar_arth ) then
!
        do k = 2_i4b, n
            arth_iv(:,k) = arth_iv(:,k-1_i4b) + increment(:)
        end do
!
    else
!
        do k = 2_i4b, npar2_arth
            arth_iv(:,k) = arth_iv(:,k-1_i4b) + increment(:)
        end do
!
        temp = increment*npar2_arth
        k = npar2_arth
!
        do
!
            if ( k>=n ) exit
!
            l = min( k, n-k )
            k2 = k + k
            arth_iv(:,k+1_i4b:min(k2,n)) = arth_iv(:,1_i4b:l) + spread(temp,2,l)
            temp = temp + temp
            k = k2
!
        end do
!
    end if
!
    end function arth_iv
!
! =========================================================================================
!
    function arth_rv( first, increment, n )
!
! Purpose
! _______
!
!   Returns a two-dimensional real array containing size(FIRST) = size(INCREMENT) 
!   arithmetic progressions of length N whose first values are FIRST(:) and whose 
!   increments are INCREMENT(:).
!
!   It is assumed that the vector arguments FIRST and INCREMENT have the same length.
!
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:),  intent(in) :: first, increment
    integer(i4b),              intent(in) :: n
!
    real(stnd), dimension(size(first),n) :: arth_rv
    real(stnd), dimension(size(first))   :: temp
    integer(i4b)                         :: l, k, k2
!
    if ( n>0_i4b ) arth_rv(:,1_i4b) = first(:)
!
    if ( n<=npar_arth ) then
!
        do k = 2_i4b, n
            arth_rv(:,k) = arth_rv(:,k-1_i4b) + increment(:)
        end do
!
    else
!
        do k = 2_i4b, npar2_arth
            arth_rv(:,k) = arth_rv(:,k-1_i4b) + increment(:)
        end do
!
        temp = increment*npar2_arth
        k = npar2_arth
!
        do
!
            if ( k>=n ) exit
!
            l = min( k, n-k )
            k2 = k + k
            arth_rv(:,k+1_i4b:min(k2,n)) = arth_rv(:,1_i4b:l) + spread(temp,2,l)
            temp = temp + temp
            k = k2
!
        end do
!
    end if
!
    end function arth_rv
!
! =========================================================================================
!
    function arth_cv( first, increment, n )
!
! Purpose
! _______
!
!   Returns a two-dimensional complex array containing size(FIRST) = size(INCREMENT) 
!   arithmetic progressions of length N whose first values are FIRST(:) and whose 
!   increments are INCREMENT(:).
!
!   It is assumed that the vector arguments FIRST and INCREMENT have the same length.
!
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:),  intent(in) :: first, increment
    integer(i4b),                 intent(in) :: n
!
    complex(stnd), dimension(size(first),n) :: arth_cv
    complex(stnd), dimension(size(first))   :: temp
    integer(i4b)                         :: l, k, k2
!
    if ( n>0_i4b ) arth_cv(:,1_i4b) = first(:)
!
    if ( n<=npar_arth ) then
!
        do k = 2_i4b, n
            arth_cv(:,k) = arth_cv(:,k-1_i4b) + increment(:)
        end do
!
    else
!
        do k = 2_i4b, npar2_arth
            arth_cv(:,k) = arth_cv(:,k-1_i4b) + increment(:)
        end do
!
        temp = increment*npar2_arth
        k = npar2_arth
!
        do
!
            if ( k>=n ) exit
!
            l = min( k, n-k )
            k2 = k + k
            arth_cv(:,k+1_i4b:min(k2,n)) = arth_cv(:,1_i4b:l) + spread(temp,2,l)
            temp = temp + temp
            k = k2
!
        end do
!
    end if
!
    end function arth_cv
!
!
! =========================================================================================
!                 ARRAY FUNCTIONS RETURNING A GEOMETRIC PROGRESSION.
!                 ADAPTED FROM Numerical Recipes.
! =========================================================================================
!
!
    function geop_i( first, factor, n )
!
! Purpose
! _______
!
!   Returns an one-dimensional integer array of length N containing a geometric 
!   progression whose first value is FIRST and whose multiplier is FACTOR.
!
!
! __________________________________________________________________________________________
!
    integer(i4b), intent(in) :: first, factor, n
!
    integer(i4b), dimension(n) :: geop_i
    integer(i4b)               :: k, k2, temp
!
    if ( n>0_i4b ) geop_i(1_i4b) = first
!
    if ( n<=npar_geop ) then
!
        do k = 2_i4b, n
            geop_i(k) = geop_i(k-1_i4b)*factor
        end do
!
    else
!
        do k = 2_i4b, npar2_geop
            geop_i(k) = geop_i(k-1_i4b)*factor
        end do
!
        temp = factor**npar2_geop
        k = npar2_geop
!
        do
!
            if ( k>=n ) exit
!
                k2 = k + k
                geop_i(k+1_i4b:min(k2,n)) = temp*geop_i(1_i4b:min(k,n-k))
                temp = temp*temp
                k = k2
!
        end do
!
    end if
!
    end function geop_i
!
! =========================================================================================
!
    function geop_r( first, factor, n )
!
! Purpose
! _______
!
!   Returns an one-dimensional real array of length N containing a geometric 
!   progression whose first value is FIRST and whose multiplier is FACTOR.
!   
!
! __________________________________________________________________________________________
!
    real(stnd),   intent(in) :: first, factor
    integer(i4b), intent(in) :: n
!
    real(stnd), dimension(n) :: geop_r
    real(stnd)               :: temp
    integer(i4b)             :: k, k2
!
    if ( n>0_i4b ) geop_r(1_i4b) = first
!
    if ( n<=npar_geop ) then
!
        do k = 2_i4b, n
            geop_r(k) = geop_r(k-1_i4b)*factor
        end do
!
    else
!
        do k = 2_i4b, npar2_geop
            geop_r(k) = geop_r(k-1_i4b)*factor
        end do
!
        temp = factor**npar2_geop
        k = npar2_geop
!
        do
!
            if ( k>=n ) exit
!
            k2 = k + k
            geop_r(k+1_i4b:min(k2,n)) = temp*geop_r(1_i4b:min(k,n-k))
            temp = temp*temp
            k = k2
!
        end do
!
    end if
!
    end function geop_r
!
! =========================================================================================
!
    function geop_c( first, factor, n )
!
! Purpose
! _______
!
!   Returns an one-dimensional complex array of length N containing a geometric 
!   progression whose first value is FIRST and whose multiplier is FACTOR.
!   
!
! __________________________________________________________________________________________
!
    complex(stnd), intent(in) :: first, factor
    integer(i4b),  intent(in) :: n
!
    complex(stnd), dimension(n) :: geop_c
    complex(stnd)               :: temp
    integer(i4b)                :: k, k2
!
    if ( n>0_i4b ) geop_c(1_i4b) = first
!
    if ( n<=npar_geop ) then
!
        do k = 2_i4b, n
            geop_c(k) = geop_c(k-1_i4b)*factor
        end do
!
    else
!
        do k = 2_i4b, npar2_geop
            geop_c(k) = geop_c(k-1_i4b)*factor
        end do
!
        temp = factor**npar2_geop
        k = npar2_geop
!
        do
!
            if ( k>=n ) exit
!
            k2 = k + k
            geop_c(k+1_i4b:min(k2,n)) = temp*geop_c(1_i4b:min(k,n-k))
            temp = temp*temp
            k = k2
!
        end do
!
    end if
!
    end function geop_c
!
! =========================================================================================
!
    function geop_iv( first, factor, n )
!
! Purpose
! _______
!
!   Returns a two-dimensional integer array containing size(FIRST) = size(FACTOR) 
!   geometric progressions of length N whose first values are FIRST(:) and whose 
!   multipliers are FACTOR(:).
!
!   It is assumed that the vector arguments FIRST and FACTOR have the same length.
!   
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: first, factor
    integer(i4b),               intent(in) :: n
!
    integer(i4b), dimension(size(first),n) :: geop_iv
    integer(i4b)                           :: l, k, k2, temp(size(first))
!
    if ( n>0_i4b ) geop_iv(:,1_i4b) = first(:)
!
    if ( n<=npar_geop ) then
!
        do k = 2_i4b, n
            geop_iv(:,k) = geop_iv(:,k-1_i4b)*factor(:)
        end do
!
    else
!
        do k = 2_i4b, npar2_geop
            geop_iv(:,k) = geop_iv(:,k-1_i4b)*factor(:)
        end do
!
        temp = factor**npar2_geop
        k = npar2_geop
!
        do
!
            if ( k>= n ) exit
!
            l = min( k, n-k )
            k2 = k + k
            geop_iv(:,k+1_i4b:min(k2,n)) = geop_iv(:,1_i4b:l)*spread(temp,2,l)
            temp = temp*temp
            k = k2
!
        end do
!
    end if
!
    end function geop_iv
!
! =========================================================================================
!
    function geop_rv( first, factor, n )
!
! Purpose
! _______
!
!   Returns a two-dimensional real array containing size(FIRST) = size(FACTOR) 
!   geometric progressions of length N whose first values are FIRST(:) and whose 
!   multipliers are FACTOR(:).
!
!   It is assumed that the vector arguments FIRST and FACTOR have the same length.
!   
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:), intent(in) :: first, factor
    integer(i4b),             intent(in) :: n
!
    real(stnd), dimension(size(first),n) :: geop_rv
    real(stnd), dimension(size(first))   :: temp
    integer(i4b)                         :: l, k, k2
!
    if ( n>0_i4b ) geop_rv(:,1_i4b) = first(:)
!
    if ( n<=npar_geop ) then
!
        do k = 2_i4b, n
            geop_rv(:,k) = geop_rv(:,k-1_i4b)*factor(:)
        end do
!
    else
!
        do k = 2_i4b, npar2_geop
            geop_rv(:,k) = geop_rv(:,k-1_i4b)*factor(:)
        end do
!
        temp = factor**npar2_geop
        k = npar2_geop
!
        do
!
            if ( k>= n ) exit
!
            l = min( k, n-k )
            k2 = k+k
            geop_rv(:,k+1_i4b:min(k2,n)) = geop_rv(:,1_i4b:l)*spread(temp,2,l)
            temp = temp*temp
            k = k2
!
        end do
!
    end if
!
    end function geop_rv
!
! =========================================================================================
!
    function geop_cv( first, factor, n )
!
! Purpose
! _______
!
!   Returns a two-dimensional complex array containing size(FIRST) = size(FACTOR) 
!   geometric progressions of length N whose first values are FIRST(:) and whose 
!   multipliers are FACTOR(:).
!
!   It is assumed that the vector arguments FIRST and FACTOR have the same length.
!   
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:), intent(in) :: first, factor
    integer(i4b),                intent(in) :: n
!
    complex(stnd), dimension(size(first),n) :: geop_cv
    complex(stnd), dimension(size(first))   :: temp
    integer(i4b)                            :: l, k, k2
!
    if ( n>0_i4b ) geop_cv(:,1_i4b) = first(:)
!
    if ( n<=npar_geop ) then
!
        do k = 2_i4b, n
            geop_cv(:,k) = geop_cv(:,k-1_i4b)*factor(:)
        end do
!
    else
!
        do k = 2_i4b, npar2_geop
            geop_cv(:,k) = geop_cv(:,k-1_i4b)*factor(:)
        end do
!
        temp = factor**npar2_geop
        k = npar2_geop
!
        do
!
            if ( k>= n ) exit
!
            k2 = k+k
            l = min( k, n-k )
            geop_cv(:,k+1_i4b:min(k2,n)) = geop_cv(:,1_i4b:l)*spread(temp,2,l)
            temp = temp*temp
            k = k2
!
        end do
!
    end if
!
    end function geop_cv
!
!
! =========================================================================================
!             ARRAY FUNCTIONS RETURNING A CUMULATIVE SUM ON AN ARRAY.
!             ADAPTED FROM Numerical Recipes.
! =========================================================================================
!
!
    recursive function cumsum_i( arr, seed ) result( ans )
!
! Purpose
! _______
!
!   Returns a rank one integer array containing the cumulative sum 
!   of the rank one integer array ARR. If the optional argument SEED is present,
!   it is added to all components of the result.
!
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: arr
    integer(i4b), optional,     intent(in) :: seed
!
    integer(i4b), dimension(size(arr)) :: ans
    integer(i4b)                       :: n, j
!
    n = size( arr )
!
    if ( n==0_i4b ) return
!
    if ( present(seed) )  then
        ans(1_i4b) = arr(1_i4b) + seed
    else
        ans(1_i4b) = arr(1_i4b)
    end if
!
    if ( n<npar_cumsum ) then
!
        do j = 2_i4b, n
            ans(j) = ans(j-1_i4b) + arr(j)
        end do
!
    else
!
        ans(2_i4b:n:2_i4b) = cumsum_i(arr(2_i4b:n:2_i4b)+arr(1_i4b:n-1_i4b:2_i4b),seed)
        ans(3_i4b:n:2_i4b) = ans(2_i4b:n-1_i4b:2_i4b)+arr(3_i4b:n:2_i4b)
!
    end if
!
    end function cumsum_i
!
! =========================================================================================
!
    recursive function cumsum_r( arr, seed ) result( ans )
!
! Purpose
! _______
!
!   Returns a rank one real array containing the cumulative sum 
!   of the rank one real array ARR. If the optional argument SEED is present,
!   it is added to all components of the result.
!   
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:), intent(in) :: arr
    real(stnd), optional,     intent(in) :: seed
!
    real(stnd), dimension(size(arr)) :: ans
    integer(i4b)                     :: n, j
!
    n = size( arr )
!
    if ( n==0_i4b ) return
!
    if ( present(seed) )  then
        ans(1_i4b) = arr(1_i4b) + seed
    else
        ans(1_i4b) = arr(1_i4b)
    end if
!
    if ( n<npar_cumsum ) then
!
        do j = 2_i4b, n
            ans(j) = ans(j-1_i4b) + arr(j)
        end do
!
    else
!
        ans(2_i4b:n:2_i4b) = cumsum_r(arr(2_i4b:n:2_i4b)+arr(1_i4b:n-1_i4b:2_i4b),seed)
        ans(3_i4b:n:2_i4b) = ans(2_i4b:n-1_i4b:2_i4b)+arr(3_i4b:n:2_i4b)
!
    end if
!
    end function cumsum_r
!
! =========================================================================================
!
    recursive function cumsum_c( arr, seed ) result( ans )
!
! Purpose
! _______
!
!   Returns a rank one complex array containing the cumulative sum 
!   of the rank one complex array ARR. If the optional argument SEED is present,
!   it is added to all components of the result.
!   
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:), intent(in) :: arr
    complex(stnd), optional ,    intent(in) :: seed
!
    complex(stnd), dimension(size(arr)) :: ans
    integer(i4b)                        :: n, j
!
    n = size( arr )
!
    if ( n==0_i4b ) return
!
    if ( present(seed) )  then
        ans(1_i4b) = arr(1_i4b) + seed
    else
        ans(1_i4b) = arr(1_i4b)
    end if
!
    if ( n<npar_cumsum ) then
!
        do j = 2_i4b, n
            ans(j) = ans(j-1_i4b) + arr(j)
        end do
!
    else
!
        ans(2_i4b:n:2_i4b) = cumsum_c(arr(2_i4b:n:2_i4b)+arr(1_i4b:n-1_i4b:2_i4b),seed)
        ans(3_i4b:n:2_i4b) = ans(2_i4b:n-1_i4b:2_i4b)+arr(3_i4b:n:2_i4b)
!
    end if
!
    end function cumsum_c
!
!
! =========================================================================================
!                 ARRAY FUNCTIONS RETURNING A CUMULATIVE PRODUCT ON AN ARRAY.
!                 ADAPTED FROM Numerical Recipes.
! =========================================================================================
!
!
    recursive function cumprod_i( arr, seed ) result( ans )
!
! Purpose
! _______
!
!   Returns a rank one integer array containing the cumulative product 
!   of the rank one integer array ARR. If the optional argument SEED is present,
!   it is multiplied into all components of the result.
!   
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: arr
    integer(i4b), optional,     intent(in) :: seed
!
    integer(i4b), dimension(size(arr)) :: ans
    integer(i4b) :: n,j
!
    n = size( arr )
!
    if (n == 0_i4b) return
!
    if ( present(seed) )  then
        ans(1_i4b) = arr(1_i4b)*seed
    else
        ans(1_i4b) = arr(1_i4b)
    end if
!
    if ( n<npar_cumprod ) then
!
        do j = 2_i4b, n
            ans(j) = ans(j-1_i4b)*arr(j)
        end do
!
    else
!
        ans(2_i4b:n:2_i4b) = cumprod_i(arr(2_i4b:n:2_i4b)*arr(1_i4b:n-1_i4b:2_i4b),seed)
        ans(3_i4b:n:2_i4b) = ans(2_i4b:n-1_i4b:2_i4b)*arr(3_i4b:n:2_i4b)
!
    end if
!
    end function cumprod_i
!
! =========================================================================================
!
    recursive function cumprod_r( arr, seed ) result( ans )
!
! Purpose
! _______
!
!   Returns a rank one real array containing the cumulative product 
!   of the rank one real array ARR. If the optional argument SEED is present,
!   it is multiplied into all components of the result.
!   
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:), intent(in) :: arr
    real(stnd), optional,     intent(in) :: seed
!
    real(stnd), dimension(size(arr)) :: ans
    integer(i4b) :: n,j
!
    n = size( arr )
!
    if (n == 0_i4b) return
!
    if ( present(seed) )  then
        ans(1_i4b) = arr(1_i4b)*seed
    else
        ans(1_i4b) = arr(1_i4b)
    end if
!
    if ( n<npar_cumprod ) then
!
        do j = 2_i4b, n
            ans(j) = ans(j-1_i4b)*arr(j)
        end do
!
    else
!
        ans(2_i4b:n:2_i4b) = cumprod_r(arr(2_i4b:n:2_i4b)*arr(1_i4b:n-1_i4b:2_i4b),seed)
        ans(3_i4b:n:2_i4b) = ans(2_i4b:n-1_i4b:2_i4b)*arr(3_i4b:n:2_i4b)
!
    end if
!
    end function cumprod_r
!
! =========================================================================================
!
    recursive function cumprod_c( arr, seed ) result( ans )
!
! Purpose
! _______
!
!   Returns a rank one complex array containing the cumulative product 
!   of the rank one complex array ARR. If the optional argument SEED is present,
!   it is multiplied into all components of the result.
!   
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:),  intent(in) :: arr
    complex(stnd), optional ,     intent(in) :: seed
!
    complex(stnd), dimension(size(arr)) :: ans
    integer(i4b) :: n,j
!
    n = size( arr )
!
    if (n == 0_i4b) return
!
    if ( present(seed) )  then
        ans(1_i4b) = arr(1_i4b)*seed
    else
        ans(1_i4b) = arr(1_i4b)
    end if
!
    if ( n<npar_cumprod ) then
!
        do j = 2_i4b, n
            ans(j) = ans(j-1_i4b)*arr(j)
        end do
!
    else
!
        ans(2_i4b:n:2_i4b) = cumprod_c(arr(2_i4b:n:2_i4b)*arr(1_i4b:n-1_i4b:2_i4b),seed)
        ans(3_i4b:n:2_i4b) = ans(2_i4b:n-1_i4b:2_i4b)*arr(3_i4b:n:2_i4b)
!
    end if
!
    end function cumprod_c
!
!
! =========================================================================================
!                            FUNCTIONS FOR POLYNOMIAL EVALUATION.
!                            ADAPTED FROM Numerical Recipes.
! =========================================================================================
!
!
    function poly_rr( x, coeffs )
!
! Purpose
! _______
!
!   Returns a real scalar containing the result of evaluating the polynomial P(X) for X real
!   with one-dimensional real coefficient vector COEFFS
!
!               P(X) = COEFFS(1)  +  COEFFS(2) * X  +  COEFFS(3) * X**(2)  + ...
! 
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : zero
    use Char_Constants,   only : allocate_error
!
    real(stnd),               intent(in) :: x
    real(stnd), dimension(:), intent(in) :: coeffs
!
    real(stnd) :: poly_rr
    real(stnd) :: pow
    real(stnd), dimension(:), allocatable :: vec
    integer(i4b) :: i, n, nn
    integer      :: iok
!
    character(len=*),  parameter :: name_proc='poly'
!
    n = size( coeffs )
!
    if ( n<=0_i4b ) then
!
        poly_rr = zero
!
    else if ( n<npar_poly ) then
!
        poly_rr = coeffs(n)
!
        do i = n-1_i4b, 1_i4b, -1_i4b
            poly_rr = x*poly_rr + coeffs(i)
        end do
!
    else
!
        allocate( vec(n+1_i4b), stat = iok )
!
        if ( iok/=0 ) then
            call merror( name_proc//allocate_error )
        end if
!
        pow = x
        vec(1_i4b:n) = coeffs
!
        do
!
            vec(n+1_i4b) = zero
            nn = ishft(n+1_i4b,-1_i4b)
            vec(1_i4b:nn) = vec(1_i4b:n:2_i4b) + pow*vec(2_i4b:n+1_i4b:2_i4b)
!
            if ( nn==1_i4b ) exit
!
            pow = pow*pow
            n = nn
!
        end do
!
        poly_rr = vec(1_i4b)
!
        deallocate( vec )
!
    end if
!
    end function poly_rr
!
! =========================================================================================
!
    function poly_rc( x, coeffs )
!
! Purpose
! _______
!
!   Returns a complex scalar containing the result of evaluating the polynomial P(X) 
!   for X complex with one-dimensional real coefficient vector COEFFS
!
!               P(X) = COEFFS(1)  +  COEFFS(2) * X  +  COEFFS(3) * X**(2)  +  ...
!   
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : zero
    use Char_Constants,   only : allocate_error
!
    complex(stnd),            intent(in) :: x
    real(stnd), dimension(:), intent(in) :: coeffs
!
    complex(stnd) :: poly_rc
    complex(stnd) :: pow
    complex(stnd), dimension(:), allocatable :: vec
    integer(i4b) :: i, n, nn
    integer      :: iok
!
    character(len=*),  parameter :: name_proc='poly'
!
    n = size( coeffs )
!
    if ( n<=0_i4b ) then
!
        poly_rc = zero
!
    else if ( n<npar_poly ) then
!
        poly_rc = coeffs(n)
!
        do i = n-1_i4b, 1_i4b, -1_i4b
            poly_rc = x*poly_rc + coeffs(i)
        end do
!
    else
!
        allocate( vec(n+1_i4b), stat = iok )
!
        if ( iok/=0 ) then
            call merror( name_proc//allocate_error )
        end if
!
        pow = x
        vec(1_i4b:n) = coeffs
!
        do
!
            vec(n+1_i4b) = zero
            nn = ishft(n+1_i4b,-1_i4b)
            vec(1_i4b:nn) = vec(1_i4b:n:2_i4b) + pow*vec(2_i4b:n+1_i4b:2_i4b)
!
            if ( nn==1_i4b ) exit
!
            pow = pow*pow
            n = nn
!
        end do
!
        poly_rc = vec(1_i4b)
!
        deallocate( vec )
!
    end if
!
    end function poly_rc
!
! =========================================================================================
!
    function poly_cc( x, coeffs )
!
! Purpose
! _______
!
!   Returns a complex scalar containing the result of evaluating the polynomial P(X) 
!   for X complex with one-dimensional complex coefficient vector COEFFS
!
!               P(X) = COEFFS(1)  +  COEFFS(2) * X  +  COEFFS(3) * X**(2)  +  ...
!   
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : zero
    use Char_Constants,   only : allocate_error
!
    complex(stnd),               intent(in) :: x
    complex(stnd), dimension(:), intent(in) :: coeffs
!
    complex(stnd) :: poly_cc
    complex(stnd) :: pow
    complex(stnd), dimension(:), allocatable :: vec
    integer(i4b) :: i, n, nn
    integer      :: iok
!
    character(len=*),  parameter :: name_proc='poly'
!
    n = size( coeffs )
!
    if ( n<=0_i4b ) then
!
        poly_cc = zero
!
    else if ( n<npar_poly ) then
!
        poly_cc = coeffs(n)
!
        do i = n-1_i4b, 1_i4b, -1_i4b
            poly_cc = x*poly_cc+coeffs(i)
        end do
!
    else
!
        allocate( vec(n+1_i4b), stat = iok )
!
        if ( iok/=0 ) then
            call merror( name_proc//allocate_error )
        end if
!
        pow = x
        vec(1_i4b:n) = coeffs
!
        do
!
            vec(n+1_i4b) = zero
            nn = ishft(n+1_i4b,-1_i4b)
            vec(1_i4b:nn) = vec(1_i4b:n:2_i4b) + pow*vec(2_i4b:n+1_i4b:2_i4b)
!
            if ( nn==1_i4b ) exit
!
            pow = pow*pow
            n = nn
!
        end do
!
        poly_cc = vec(1_i4b)
!
        deallocate( vec )
!
    end if
!
    end function poly_cc
!
! =========================================================================================
!
    function poly_rrv( x, coeffs )
!
! Purpose
! _______
!
!   Returns a real vector containing the results of evaluating the polynomials P(X(:))
!   for X(:) real with one-dimensional real coefficient vector COEFFS
!
!               P(X(:)) = COEFFS(1)  +  COEFFS(2) * X(:)  +  COEFFS(3) * X(:)**(2)  +  ...
!   
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : zero
!
    real(stnd), dimension(:), intent(in) :: coeffs, x
!
    real(stnd), dimension(size(x)) :: poly_rrv
    integer(i4b) :: i, n, m
!
    m = size( coeffs )
    n = size( x )
!
    if ( m<=0_i4b ) then
!
        poly_rrv = zero
!
    else if ( m<n .or. m<npar_poly ) then
!
        poly_rrv = coeffs(m)
!
        do i = m-1_i4b, 1_i4b, -1_i4b
            poly_rrv = x*poly_rrv + coeffs(i)
        end do
!
    else
!
        do i = 1_i4b, n
            poly_rrv(i) = poly_rr( x(i), coeffs )
        end do
!
    end if
!
    end function poly_rrv
!
! =========================================================================================
!
    function poly_msk_rrv( x, coeffs, mask )
!
! Purpose
! _______
!
!   Returns a real vector containing the results of evaluating the polynomials P(X(:))
!   for X(:) real with one-dimensional real coefficient vector COEFFS
!
!               P(X(:)) = COEFFS(1)  +  COEFFS(2) * X(:)  +  COEFFS(3) * X(:)**(2)  +  ...
!   
!   under the control of the logical argument MASK. If MASK(i) = false, the polynomial is
!   not evaluated at X(i).
!   
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : zero
!
    real(stnd),   dimension(:), intent(in) :: coeffs, x
    logical(lgl), dimension(:), intent(in) :: mask
!
    real(stnd), dimension(size(x)) :: poly_msk_rrv
!
    poly_msk_rrv = unpack( poly_rrv(pack(x,mask),coeffs), mask, zero )
!
    end function poly_msk_rrv
!
! =========================================================================================
!
    recursive function poly_term_rr( coeffs, x ) result( u )
!
! Purpose
! _______
!
!   Returns a real array of size(COEFFS) containing the partial cumulants of the polynomial
!   with real coefficients COEFFS evaluated at the real scalar X.
!   On entry, the coefficients in COEFFS are arranged from highest order to lowest-order
!   coefficients.
!   
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:), intent(in) :: coeffs
    real(stnd),               intent(in) :: x
!
    real(stnd), dimension(size(coeffs)) :: u
    integer(i4b) :: n, j
!
    n = size( coeffs )
!
    if ( n<=0_i4b ) return
!
    u(1) = coeffs(1_i4b)
!
    if ( n<npar_polyterm ) then
!
        do j = 2_i4b, n
            u(j) = coeffs(j) + x*u(j-1_i4b)
        end do
!
    else
!
        u(2_i4b:n:2_i4b) = poly_term_rr( coeffs(2_i4b:n:2_i4b) + coeffs(1_i4b:n-1_i4b:2_i4b)*x, x*x )
        u(3_i4b:n:2_i4b) = coeffs(3_i4b:n:2_i4b) + x*u(2_i4b:n-1_i4b:2_i4b)
!
    end if
!
    end function poly_term_rr
!
! =========================================================================================
!
    recursive function poly_term_cc( coeffs, x ) result( u )
!
! Purpose
! _______
!
!   Returns a complex array of size(COEFFS) containing the partial cumulants of the polynomial
!   with complex coefficients COEFFS evaluated at the complex scalar X.
!   On entry, the coefficients in COEFFS are arranged from highest order to lowest-order
!   coefficients.
!   
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:), intent(in) :: coeffs
    complex(stnd),               intent(in) :: x
!
    complex(stnd), dimension(size(coeffs)) :: u
    integer(i4b) :: n, j
!
    n = size( coeffs )
!
    if ( n<=0_i4b ) return
!
    u(1) = coeffs(1_i4b)
!
    if ( n<npar_polyterm ) then
!
        do j = 2_i4b, n
            u(j) = coeffs(j) + x*u(j-1_i4b)
        end do
!
    else
!
        u(2_i4b:n:2_i4b) = poly_term_cc( coeffs(2_i4b:n:2_i4b) + coeffs(1_i4b:n-1_i4b:2_i4b)*x, x*x )
        u(3_i4b:n:2_i4b) = coeffs(3_i4b:n:2_i4b) + x*u(2_i4b:n-1_i4b:2_i4b)
!
    end if
!
    end function poly_term_cc
!
! =========================================================================================
!
    function zroots_unity( n, nn )
!
! Purpose
! _______
!
!   Complex function returning a complex array containing nn consecutive powers of 
!   the nth complex root of unity.
!
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : one, twopi
!
    integer(i4b), intent(in) :: n, nn
!
    complex(stnd), dimension(nn) :: zroots_unity
    integer(i4b) :: k
    real(stnd)   :: theta
!
    zroots_unity(1_i4b) = one
    theta = twopi/real( n, stnd )
    k = 1_i4b
!
    do
!
        if ( k>=nn ) exit
!
        zroots_unity(k+1_i4b) = cmplx( cos(k*theta), sin(k*theta), stnd )
        zroots_unity(k+2_i4b:min(2_i4b*k,nn)) = zroots_unity(k+1_i4b)*           &
                                                zroots_unity(2_i4b:min(k,nn-k))
        k = 2_i4b*k
!
    end do
!
    end function zroots_unity
!
!
! =========================================================================================
!       ROUTINES FOR OUTER OPERATIONS ON VECTORS (PARALLEL OR VECTORIZED VERSIONS)
!       ADAPTED AND EXTENDED FROM Numerical Recipes.
! =========================================================================================
!
!
    subroutine update_rk1_i( mat, u, v )
!
! Purpose
! _______
!
!   Updates the integer matrix MAT with the outer sum of the two integer vectors U and V :
!
!                        MAT = MAT + U * V'
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:,:), intent(inout) :: mat
    integer(i4b), dimension(:),   intent(in)    :: u, v
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( u )
    n = size( v )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
!
        mat(1_i4b:m,i) = mat(1_i4b:m,i) + u(1_i4b:m)*v(i)
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine update_rk1_i
!
! =========================================================================================
!
    subroutine update_rk1_r( mat, u, v )
!
! Purpose
! _______
!
!   Updates the real matrix MAT with the outer sum of the two reals vectors U and V :
!
!                        MAT = MAT + U * V'
!   
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : axpy
#endif
!
    real(stnd), dimension(:,:), intent(inout) :: mat
    real(stnd), dimension(:),   intent(in)    :: u, v
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( u )
    n = size( v )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
!
#ifdef _BLAS
        call axpy( m, v(i),  u(1_i4b:m),  1_i4b, mat(1_i4b:m,i), 1_i4b )
#else
        mat(1_i4b:m,i) = mat(1_i4b:m,i) + u(1_i4b:m)*v(i)
#endif
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine update_rk1_r
!
! =========================================================================================
!
    subroutine update_rk1_c( mat, u, v )
!
! Purpose
! _______
!
!   Updates the complex matrix MAT with the outer sum of the two complex vectors U and V :
!
!                        MAT = MAT + U * V'
!   
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : axpy
#endif
!
    complex(stnd), dimension(:,:), intent(inout) :: mat
    complex(stnd), dimension(:),   intent(in)    :: u, v
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( u )
    n = size( v )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
!
#ifdef _BLAS
        call axpy( m, v(i),  u(1_i4b:m),  1_i4b, mat(1_i4b:m,i), 1_i4b )
#else
        mat(1_i4b:m,i) = mat(1_i4b:m,i) + u(1_i4b:m)*v(i)
#endif
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine update_rk1_c
!
! =========================================================================================
!
    subroutine update_rk2_i( mat, u, v, u2, v2 )
!
! Purpose
! _______
!
!   Updates the integer matrix MAT with the outer sums of the integer vectors U, V, U2 and V2:
!
!                        MAT = MAT + U * V' + U2 * V2'
!
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:,:), intent(inout) :: mat
    integer(i4b), dimension(:),   intent(in)    :: u, v, u2, v2
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( u )
    n = size( v )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
!
        mat(1_i4b:m,i) = mat(1_i4b:m,i) + u(1_i4b:m)*v(i) + u2(1_i4b:m)*v2(i)
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine update_rk2_i
!
! =========================================================================================
!
    subroutine update_rk2_r( mat, u, v, u2, v2 )
!
! Purpose
! _______
!
!   Updates the real matrix MAT with the outer sums of the real vectors U, V, U2 and V2:
!
!                        MAT = MAT + U * V' + U2 * V2'
!   
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : axpy
#endif
!
    real(stnd), dimension(:,:), intent(inout) :: mat
    real(stnd), dimension(:),   intent(in)    :: u, v, u2, v2
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( u )
    n = size( v )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
!
#ifdef _BLAS
        call axpy( m, v(i),  u(1_i4b:m),  1_i4b, mat(1_i4b:m,i), 1_i4b )
        call axpy( m, v2(i), u2(1_i4b:m), 1_i4b, mat(1_i4b:m,i), 1_i4b )
#else
        mat(1_i4b:m,i) = mat(1_i4b:m,i) + u(1_i4b:m)*v(i) + u2(1_i4b:m)*v2(i)
#endif
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine update_rk2_r
!
! =========================================================================================
!
    subroutine update_rk2_c( mat, u, v, u2, v2 )
!
! Purpose
! _______
!
!   Updates the complex matrix MAT with the outer sums of the complex vectors U, V, U2 and V2:
!
!                        MAT = MAT + U * V' + U2 * V2'
!   
!
! __________________________________________________________________________________________
!
#ifdef _BLAS
    use BLAS_interfaces,  only : axpy
#endif
!
    complex(stnd), dimension(:,:), intent(inout) :: mat
    complex(stnd), dimension(:),   intent(in)    :: u, v, u2, v2
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( u )
    n = size( v )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
!
#ifdef _BLAS
        call axpy( m, v(i),  u(1_i4b:m),  1_i4b, mat(1_i4b:m,i), 1_i4b )
        call axpy( m, v2(i), u2(1_i4b:m), 1_i4b, mat(1_i4b:m,i), 1_i4b )
#else
        mat(1_i4b:m,i) = mat(1_i4b:m,i) + u(1_i4b:m)*v(i) + u2(1_i4b:m)*v2(i)
#endif
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine update_rk2_c
!
! =========================================================================================
!
    function outerprod_i( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer product of the two integer vectors A and B .
!
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: a, b
!
    integer(i4b), dimension(size(a),size(b)) :: outerprod_i
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
!
        outerprod_i(1_i4b:m,i) = a(1_i4b:m)*b(i)
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outerprod_i
!
! =========================================================================================
!
    function outerprod_r( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer product of the two real vectors A and B .
!   
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:), intent(in) :: a, b
!
    real(stnd), dimension(size(a),size(b)) :: outerprod_r
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
!
        outerprod_r(1_i4b:m,i) = a(1_i4b:m)*b(i)
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outerprod_r
!
! =========================================================================================
!
    function outerprod_c( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer product of the two complex vectors A and B .
!   
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:), intent(in) :: a, b
!
    complex(stnd), dimension(size(a),size(b)) :: outerprod_c
    integer(i4b) :: i, n, m
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
!
        outerprod_c(1_i4b:m,i) = a(1_i4b:m)*b(i)
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outerprod_c
!
! =========================================================================================
!
    function outerdiv_r( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer quotient of the two real vectors A and B .
!
!
! Further Details
! _______________
!
!   It is assumed that none of the elements of B is zero.
!
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : one
!
    real(stnd), dimension(:), intent(in) :: a, b
!
    real(stnd), dimension(size(a),size(b)) :: outerdiv_r
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
        outerdiv_r(1_i4b:m,i) =  a(1_i4b:m)*(one/b(i))
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outerdiv_r
!
! =========================================================================================
!
    function outerdiv_c( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer quotient of the two complex vectors A and B .
!
!
! Further Details
! _______________
!
!   It is assumed that none of the elements of B is zero.
!   
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : one
!
    complex(stnd), dimension(:), intent(in) :: a, b
!
    complex(stnd), dimension(size(a),size(b)) :: outerdiv_c
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
        outerdiv_c(1_i4b:m,i) =  a(1_i4b:m)*(one/b(i))
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outerdiv_c
!
! =========================================================================================
!
    function outersum_i( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer sum of the two integer vectors A and B .
!
!
! Further Details
! _______________
!
!   It is assumed that none of the elements of B is zero.
!   
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: a, b
!
    integer(i4b), dimension(size(a),size(b)) :: outersum_i
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
        outersum_i(1_i4b:m,i) =  a(1_i4b:m) + b(i)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outersum_i
!
! =========================================================================================
!
    function outersum_r( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer sum of the two real vectors A and B .
!   
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:), intent(in) :: a, b
!
    real(stnd), dimension(size(a),size(b)) :: outersum_r
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
        outersum_r(1_i4b:m,i) =  a(1_i4b:m) + b(i)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outersum_r
!
! =========================================================================================
!
    function outersum_c( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer sum of the two complex vectors A and B .
!   
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:), intent(in) :: a, b
!
    complex(stnd), dimension(size(a),size(b)) :: outersum_c
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
        outersum_c(1_i4b:m,i) =  a(1_i4b:m) + b(i)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outersum_c
!
! =========================================================================================
!
    function outerdiff_i( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer difference of the two integer vectors A and B .
!   
!
! __________________________________________________________________________________________
!
    integer(i4b), dimension(:), intent(in) :: a, b
!
    integer(i4b), dimension(size(a), size(b)) :: outerdiff_i
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
        outerdiff_i(1_i4b:m,i) =  a(1_i4b:m) - b(i)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outerdiff_i
!
! =========================================================================================
!
    function outerdiff_r( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer difference of the two real vectors A and B .
!   
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:), intent(in) :: a, b
!
    real(stnd), dimension(size(a), size(b)) :: outerdiff_r
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
        outerdiff_r(1_i4b:m,i) =  a(1_i4b:m) - b(i)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outerdiff_r
!
! =========================================================================================
!
    function outerdiff_c( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer difference of the two complex vectors A and B .
!  
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:), intent(in) :: a, b
!
    complex(stnd), dimension(size(a), size(b)) :: outerdiff_c
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
        outerdiff_c(1_i4b:m,i) =  a(1_i4b:m) - b(i)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outerdiff_c
!
! =========================================================================================
!
    function outerand( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer logical AND of two logical vectors A and B .
!   
!
! __________________________________________________________________________________________
!
    logical(lgl), dimension(:), intent(in) :: a, b
!
    logical(lgl), dimension(size(a),size(b)) :: outerand
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
        outerand(1_i4b:m,i) =  a(1_i4b:m) .and. b(i)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outerand
!
! =========================================================================================
!
    function outeror( a, b )
!
! Purpose
! _______
!
!   Returns a matrix that is the outer logical OR of two logical vectors A and B .
!   
!
! __________________________________________________________________________________________
!
    logical(lgl), dimension(:), intent(in) :: a, b
!
    logical(lgl), dimension(size(a),size(b)) :: outeror
!
    integer(i4b) :: i, m, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( a )
    n = size( b )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  n*m>=omp_limit    .and.      &
                  n>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
    do i = 1_i4b, n
        outeror(1_i4b:m,i) =  a(1_i4b:m) .or. b(i)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end function outeror
!
! =========================================================================================
!
    function triangle( upper, j, k, extra )
!
! Purpose
! _______
!
!   Return an upper (if UPPER=true) or lower (if UPPER=false) triangular logical mask.
!
! __________________________________________________________________________________________
!
    logical(lgl), intent(in) :: upper
!
    integer(i4b),           intent(in) :: j, k
    integer(i4b), optional, intent(in) :: extra
!
    logical(lgl), dimension(j,k) :: triangle
    integer(i4b) :: n, ind(j), i
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = 0_i4b
    if ( present(extra) ) n = extra
!
    ind = arth_i( 1_i4b, 1_i4b, j )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel  .and.      &
                  k*j>=omp_limit    .and.      &
                  k>1_i4b           .and.      &
                  i>1_i4b
#endif
!
    if ( upper ) then
!
#ifdef _OPENMP3
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
        do i = 1_i4b, k
            triangle(1_i4b:j,i) = ind - i < n
        end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    else
!
#ifdef _OPENMP3
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
        do i = 1_i4b, k
            triangle(1_i4b:j,i) = ind - i > -n
        end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end if
!
    end function triangle
!
!
! =========================================================================================
!        ROUTINES FOR COMPUTING NORM AND SUM OF SQUARES ON MATRICES AND VECTORS
!                        (PARALLEL OR VECTORIZED VERSIONS).
! =========================================================================================
!
!
    function abse_rv( vec )
!
! Purpose
! _______
!
!   Return the ordinary L2 norm of the real vector VEC.
!
!
! ___________________________________________________________________________________________
!
    real(stnd), intent(in) :: vec(:)
!
    real(stnd)   :: abse_rv
    integer(i4b) :: n
!
    n = size( vec )
    abse_rv = sqrt( sum( vec(:n)**2 ) )
!
    end function abse_rv
!
! =========================================================================================
!
    function abse_cv( vec )
!
! Purpose
! _______
!
!   Return the ordinary L2 norm of the complex vector VEC as a real scalar.
!
!
! ___________________________________________________________________________________________
!
    complex(stnd), intent(in) :: vec(:)
!
    real(stnd)   :: abse_cv
    integer(i4b) :: n
!
    n = size( vec )
    abse_cv = sqrt( sum( real(vec(:n))**2 + aimag(vec(:n))**2 ) )
!
    end function abse_cv
!
! =========================================================================================
!
    function abse_rm( mat )
!
! Purpose
! _______
!
!   Return the Froebenius norm of the real matrix MAT.
!
!
! ___________________________________________________________________________________________
!
    real(stnd), intent(in) :: mat(:,:)
!
    real(stnd)   :: abse_rm
    integer(i4b) :: n, m 
!
    n = size( mat, 1 )
    m = size( mat, 2 )
    abse_rm = sqrt( sum( mat(:n,:m)**2 ) )
!
    end function abse_rm
!
! =========================================================================================
!
    function abse_cm( mat )
!
! Purpose
! _______
!
!   Return the Froebenius norm of the complex matrix MAT as a real scalar.
!
!
! ___________________________________________________________________________________________
!
    complex(stnd), intent(in) :: mat(:,:)
!
    real(stnd)   :: abse_cm
    integer(i4b) :: n, m 
!
    n = size( mat, 1 )
    m = size( mat, 2 )
    abse_cm = sqrt( sum( real(mat(:n,:m))**2 + aimag(mat(:n,:m))**2 ) )
!
    end function abse_cm
!
! =========================================================================================
!
    function abse_dim_rm( mat, dim )
!
! Purpose
! _______
!
!   Return the ordinary L2 norm of the column vectors (DIM=2) or the row vectors (DIM=1) 
!   of the real matrix MAT as a real vector.
!
!
! ___________________________________________________________________________________________
!
    real(stnd),   intent(in) :: mat(:,:)
    integer(i4b), intent(in) :: dim
!
    real(stnd)   :: abse_dim_rm(size(mat,int(dim)))
    integer(i4b) :: n, m 
!
    n = size( mat, 1 )
    m = size( mat, 2 )
    abse_dim_rm(:) = sqrt( sum( mat(:n,:m)**2, dim=int(3_i4b-dim) ) )
!
    end function abse_dim_rm
!
! =========================================================================================
!
    function abse_dim_cm( mat, dim )
!
! Purpose
! _______
!
!   Return the ordinary L2 norm of the column vectors (DIM=2) or the row vectors (DIM=1) 
!   of the complex matrix MAT as a real vector.
!
!
! ___________________________________________________________________________________________
!
    complex(stnd), intent(in) :: mat(:,:)
    integer(i4b),  intent(in) :: dim
!
    real(stnd)   :: abse_dim_cm(size(mat,int(dim)))
    integer(i4b) :: n, m 
!
    n = size( mat, 1 )
    m = size( mat, 2 )
    abse_dim_cm(:) = sqrt( sum( real(mat(:n,:m))**2 + aimag(mat(:n,:m))**2, dim=int(3_i4b-dim )) )
!
    end function abse_dim_cm
!
! =========================================================================================
!
    subroutine lassq_rv( vec, scal, ssq )
!
! Purpose
! _______
!
!   LASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**(2) ) * smsq = sum( VEC**(2) ) + ( scale**(2) )*ssq,
!
!   The value of  ssq  is assumed to be non-negative 
!   and  scl  returns the value
!
!     scl = max( scale, maxval( abs( VEC ) ) ).
!
!   scale and ssq must be supplied in SCAL and SSQ and
!   scl and smsq are overwritten on SCAL and SSQ respectively.
!
!
! Arguments
! _________
!
!   VEC     (INPUT) real(stnd), dimension(:)
!           The real vector for which a scaled sum of squares is computed.
!
!   SCAL    (INPUT/OUTPUT) real(stnd)
!           On entry, the value  scale  in the equation above.
!           On exit, SCAL is overwritten with  scl , the scaling factor
!           for the sum of squares.
!
!   SSQ     (INPUT/OUTPUT) real(stnd)
!           On entry, the value  ssq  in the equation above.
!           On exit, SSQ is overwritten with  smsq , the basic sum of
!           squares from which  scl  has been factored out.
!
!
! ___________________________________________________________________________________________
!
    use Reals_Constants,   only : zero, one
    use Num_Constants,     only : lamch, machsmlnum
!
    real(stnd), dimension(:), intent(in)    :: vec
    real(stnd),               intent(inout) :: scal, ssq
!
    real(stnd)   :: hitest, smax, sqmax, cutlo, sfmin, sfmax
    integer(i4b) :: n
!
    n = size( vec )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<=0_i4b ) return
!
    cutlo  = sqrt( machsmlnum )
    hitest = one/( cutlo*real( n+1_i4b, stnd ) )
    if ( ssq==zero ) scal = one
!
!   FIND THE MAXIMUM VALUE IN vec .
!
    smax = maxval( abs(vec(:n)) )
!
    sqmax = max( scal*sqrt(ssq), smax )
!
    if ( scal==one .and. sqmax>cutlo .and. sqmax<hitest ) then
!
!       IF scal=1 AND sqmax IS GREATER THAN cutlo AND LESS THAN hitest, 
!       NO SCALING SHOULD BE NEEDED.
!
        ssq =  ssq + sum( vec(:n)**2 )
!
    else if ( smax>zero ) then
!
!       SCALE BY sqmax IF scal=1, OTHERWISE SCALE BY max(sqmax,scal) .
!
        sfmin = lamch( 's' )
        sfmax = one/sfmin
        sqmax = min(max( sfmin, sqmax ), sfmax )
!
        if ( scal==one .or. scal<sqmax ) then
            ssq  = ( ssq*( scal/sqmax ) )*( scal/sqmax )
            scal = sqmax
        end if
!
!       ADD THE SUM OF SQUARES OF VALUES OF vec SCALED BY scal .
!
        ssq =  ssq + sum( (vec(:n)/scal)**2 )
!
    end if
!
    end subroutine lassq_rv
!
! =========================================================================================
!
    subroutine lassq_cv( vec, scal, ssq )
!
! Purpose
! _______
!
!   LASSQ_CV  returns the values  scl  and  smsq  such that
!
!     ( scl**(2) ) * smsq = dot_product( VEC, VEC ) + ( scale**(2) ) * ssq,
!
!   The value of  ssq  is assumed to be non-negative 
!   and  scl  returns the value
!
!     scl = max( scale, maxval(abs(real(VEC))), maxval(abs(aimag(VEC))) ).
!
!   scale and ssq must be supplied in SCAL and SSQ and
!   scl and smsq are overwritten on SCAL and SSQ respectively.
!
!
! Arguments
! _________
!
!   VEC     (INPUT) complex(stnd), dimension(:)
!           The complex vector for which a scaled sum of squares is computed.
!
!   SCAL    (INPUT/OUTPUT) real(stnd)
!           On entry, the value  scale  in the equation above.
!           On exit, SCAL is overwritten with  scl , the scaling factor
!           for the sum of squares.
!
!   SSQ     (INPUT/OUTPUT) real(stnd)
!           On entry, the value  ssq  in the equation above.
!           On exit, SSQ is overwritten with  smsq , the basic sum of
!           squares from which  scl  has been factored out.
!
!
! ___________________________________________________________________________________________
!
    complex(stnd), dimension(:), intent(in)    :: vec
    real(stnd),                  intent(inout) :: scal, ssq
!
    integer(i4b) :: n
!
    n = size( vec )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<=0_i4b ) return
!
    call lassq_rv( real(vec(:n),stnd), scal, ssq )
    call lassq_rv( aimag(vec(:n)), scal, ssq )
!
    end subroutine lassq_cv
!
! =========================================================================================
!
    subroutine lassq_rm( mat, scal, ssq )
!
! Purpose
! _______
!
!   LASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**(2) ) * smsq = sum( MAT**(2) ) + ( scale**(2) ) * ssq,
!
!   The value of  ssq  is assumed to be non-negative 
!   and  scl  returns the value
!
!     scl = max( scale, maxval(abs(MAT)) ).
!
!   scale and ssq must be supplied in SCAL and SSQ and
!   scl and smsq are overwritten on SCAL and SSQ respectively.
!
!
! Arguments
! _________
!
!   MAT     (INPUT) real(stnd), dimension(:,:)
!           The matrix for which a scaled sum of squares is computed.
!
!   SCAL    (INPUT/OUTPUT) real(stnd)
!           On entry, the value  scale  in the equation above.
!           On exit, SCAL is overwritten with  scl , the scaling factor
!           for the sum of squares.
!
!   SSQ     (INPUT/OUTPUT) real(stnd)
!           On entry, the value  ssq  in the equation above.
!           On exit, SSQ is overwritten with  smsq , the basic sum of
!           squares from which  scl  has been factored out.
!
!
! ___________________________________________________________________________________________
!
    use Reals_Constants,   only : zero, one
    use Num_Constants,     only : lamch, machsmlnum
!
    real(stnd), dimension(:,:), intent(in)    :: mat
    real(stnd),                 intent(inout) :: scal, ssq
!
    real(stnd)   :: hitest, smax, sqmax, cutlo, sfmin, sfmax
    integer(i4b) :: n, m, nm
!
    n  = size( mat, 1 )
    m  = size( mat, 2 )
    nm = n*m
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( nm<=0_i4b ) return
!
    cutlo  = sqrt( machsmlnum )
    hitest = one / ( cutlo*real( nm+1_i4b, stnd ) )
    if ( ssq==zero ) scal = one
!
!   FIND THE MAXIMUM VALUE IN mat .
!
    smax = maxval( abs(mat(:n,:m)) )
!
    sqmax = max( scal*sqrt(ssq), smax )
!
    if ( scal==one .and. sqmax>cutlo .and. sqmax<hitest ) then
!
!       IF scal=1 AND sqmax IS GREATER THAN cutlo AND LESS THAN hitest, 
!       NO SCALING SHOULD BE NEEDED.
!
        ssq =  ssq + sum( mat(:n,:m)**2 )
!
    else if ( smax>zero ) then
!
!       SCALE BY sqmax IF scal=1, OTHERWISE SCALE BY max(sqmax,scal) .
!
        sfmin = lamch( 's' )
        sfmax = one/sfmin
        sqmax = min(max( sfmin, sqmax ), sfmax )
!
        if ( scal==one .or. scal<sqmax ) then
            ssq  = ( ssq*( scal/sqmax ) )*( scal/sqmax )
            scal = sqmax
        end if
!
!       ADD THE SUM OF SQUARES OF VALUES OF vec SCALED BY scal .
!
        ssq =  ssq + sum( (mat(:n,:m)/scal)**2 )
!
    end if
!
    end subroutine lassq_rm
!
! =========================================================================================
!
    subroutine lassq_cm( mat, scal, ssq )
!
! Purpose
! _______
!
!   LASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**(2) ) * smsq = sum( MAT * conjg(MAT) ) + ( scale**(2) ) * ssq,
!
!   The value of  ssq  is assumed to be non-negative 
!   and  scl  returns the value
!
!     scl = max( scale, maxval(abs(real(MAT))), maxval(abs(aimag(MAT))) ).
!
!   scale and ssq must be supplied in SCAL and SSQ and
!   scl and smsq are overwritten on SCAL and SSQ respectively.
!
!
! Arguments
! _________
!
!   MAT     (INPUT) complex(stnd), dimension(:,:)
!           The complex matrix for which a scaled sum of squares is computed.
!
!   SCAL    (INPUT/OUTPUT) real(stnd)
!           On entry, the value  scale  in the equation above.
!           On exit, SCAL is overwritten with  scl , the scaling factor
!           for the sum of squares.
!
!   SSQ     (INPUT/OUTPUT) real(stnd)
!           On entry, the value  ssq  in the equation above.
!           On exit, SSQ is overwritten with  smsq , the basic sum of
!           squares from which  scl  has been factored out.
!
!
! ___________________________________________________________________________________________
!
    complex(stnd), dimension(:,:), intent(in)    :: mat
    real(stnd),                    intent(inout) :: scal, ssq
!
    integer(i4b) :: n, m
!   
    n = size( mat, 1 )
    m = size( mat, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n*m<=0_i4b ) return
!
    call lassq_rm( real(mat(:n,:m),stnd), scal, ssq )
    call lassq_rm( aimag(mat(:n,:m)), scal, ssq )
!
    end subroutine lassq_cm
!
! =========================================================================================
!
    function norm_rv( vec )
!
! Purpose
! _______
!
!   Return the Euclidean norm of the real vector VEC via the function name, so that
!
!        norm_rv := sqrt( VEC * VEC )
!
!   This is done without destructive underflow or overflow.
!
!
! ___________________________________________________________________________________________
!
    use Reals_Constants,   only : zero, one
    use Num_Constants,     only : lamch, machsmlnum
!
    real(stnd), dimension(:), intent(in) :: vec
!
    real(stnd)   :: norm_rv
    real(stnd)   :: hitest, smax, cutlo, sfmin, sfmax
    integer(i4b) :: n
!
    n = size( vec )
!
    select case ( n )
!
        case( :0_i4b )
!
            norm_rv = zero
!
        case( 1_i4b )
!
            norm_rv = abs( vec(1_i4b) )
!
        case( 2_i4b: )
!
            cutlo  = sqrt( machsmlnum )
            hitest = one/( cutlo*real( n, stnd ) )
!
!           FIND THE MAXIMUM VALUE IN vec .
!
            smax = maxval( abs(vec(:n)) )
!
            if ( smax>cutlo .and. smax<hitest ) then
!
!               IF  smax IS GREATER THAN cutlo AND LESS THAN hitest, 
!               NO SCALING SHOULD BE NEEDED.
!
                norm_rv = sqrt( sum( vec(:n)**2 ) )
!
            else if ( smax>zero ) then
!
!               COMPUTE THE SUM OF SQUARES OF VALUES OF vec SCALED BY smax .
!
                sfmin   = lamch( 's' )
                sfmax   = one/sfmin
                smax    = min(max( sfmin, smax ), sfmax ) 
                norm_rv = sqrt( sum( (vec(:n)/smax)**2 ) )*smax
!
            else
!
                norm_rv = zero
!
            end if
!
    end select
!
    end function norm_rv
!
! =========================================================================================
!
    function norm_cv( vec )
!
! Purpose
! _______
!
!   Return the Euclidean norm of the complex vector VEC via the function name, so that
!
!        norm_cv := sqrt( dot_product(VEC,VEC) )
!
!   This is done without destructive underflow or overflow.
!
!
! ___________________________________________________________________________________________
!
    use Reals_Constants,  only : zero, one
!
    complex(stnd), dimension(:), intent(in) :: vec
!
    real(stnd) :: norm_cv
!
    real(stnd)   :: scal, ssq
    integer(i4b) :: n
!
    n = size( vec )
!
    select case ( n )
!
        case( :0_i4b )
!
            norm_cv = zero
!
        case( 1_i4b )
!
            norm_cv = abs( vec(1_i4b) )
!
        case( 2_i4b: )
!
            scal = one
            ssq = zero
            call lassq_rv( real(vec(:n),stnd), scal, ssq )
            call lassq_rv( aimag(vec(:n)), scal, ssq )
            norm_cv = scal*sqrt( ssq )
!
    end select
!
    end function norm_cv
!
! =========================================================================================
!
    function norm_rm( mat )
!
! Purpose
! _______
!
!   Return the Froebenius norm of the real matrix MAT via the function name, so that
!
!        norm_rm := sqrt( MAT * MAT )
!
!   This is done without destructive underflow or overflow.
!
!
! ___________________________________________________________________________________________
!
    use Reals_Constants,   only : zero, one
    use Num_Constants,     only : lamch, machsmlnum
!
    real(stnd), dimension(:,:), intent(in) :: mat
!
    real(stnd)   :: norm_rm
    real(stnd)   :: hitest, smax, cutlo, sfmin, sfmax
    integer(i4b) :: n, m, nm
!
    n  = size( mat, 1 )
    m  = size( mat, 2 )
    nm = n*m
!
    select case ( nm )
!
        case( :0_i4b )
!
            norm_rm = zero
!
        case( 1_i4b )
!
            norm_rm = abs( mat(1_i4b,1_i4b) )
!
        case( 2_i4b: )
!
            cutlo  = sqrt( machsmlnum )
            hitest = one/( cutlo*real( nm, stnd ) )
!
!           FIND THE MAXIMUM VALUE IN mat .
!
            smax = maxval( abs(mat(:n,:m)) )
!
            if ( smax>cutlo .and. smax<hitest ) then
!
!               IF  smax IS GREATER THAN cutlo AND LESS THAN hitest, 
!               NO SCALING SHOULD BE NEEDED.
!
                norm_rm = sqrt( sum( mat(:n,:m)**2 ) )
!
            else if ( smax>zero ) then
!
!               COMPUTE THE SUM OF SQUARES OF VALUES OF mat SCALED BY smax .
!
                sfmin   = lamch( 's' )
                sfmax   = one/sfmin
                smax    = min(max( sfmin, smax ), sfmax ) 
                norm_rm = sqrt( sum( (mat(:n,:m)/smax)**2 ) )*smax
!
            else
!
                norm_rm = zero
!
            end if
!
    end select
!
    end function norm_rm
!
! =========================================================================================
!
    function norm_cm( mat )
!
! Purpose
! _______
!
!   Return the Froebenius norm of the complex matrix MAT via the function name, so that
!
!        norm_cm := sqrt( MAT * conjg(MAT) )
!
!   This is done without destructive underflow or overflow.
!
!
! ___________________________________________________________________________________________
!
    use Reals_Constants,  only : zero, one
!
    complex(stnd), dimension(:,:), intent(in) :: mat
!
    real(stnd)   :: norm_cm
    real(stnd)   :: scal, ssq
    integer(i4b) :: n, m
!   
    n = size( mat, 1 )
    m = size( mat, 2 )
!
    select case ( n*m )
!
        case( :0_i4b )
!
            norm_cm = zero
!
        case( 1_i4b )
!
            norm_cm = abs( mat(1_i4b,1_i4b) )
!
        case( 2_i4b: )
!
            scal = one
            ssq = zero
            call lassq_rm( real(mat(:n,:m),stnd), scal, ssq )
            call lassq_rm( aimag(mat(:n,:m)), scal, ssq )
            norm_cm = scal*sqrt( ssq )
!
    end select
!
    end function norm_cm
!
! =========================================================================================
!
    function norm_dim_rm( mat, dim )
!
! Purpose
! _______
!
!   Return the Euclidean norms of the columns (DIM=2) or of the rows (DIM=1) of a real 
!   matrix MAT via the function name, so that
!
!        norm_dim_rm := sqrt( sum(MAT * MAT,dim=3-dim) )
!
!   This is done without destructive underflow or overflow.
!
!
! ___________________________________________________________________________________________
!
    use Reals_Constants,  only : zero, one
    use Num_Constants,     only : lamch, machsmlnum
!
    real(stnd), dimension(:,:), intent(in) :: mat
    integer(i4b),               intent(in) :: dim
!
    real(stnd), dimension(size(mat,int(dim))) :: norm_dim_rm
!
    real(stnd)   :: smax(size(mat,int(dim))), hitest, cutlo, sfmin, sfmax
    integer(i4b) :: n, m, dim2, i
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( mat, int(dim) )
!
    if ( m<=0_i4b ) return
!
    dim2 = 3_i4b-dim
    n    = size( mat, int(dim2) )
!
    select case( n )
!
        case( :0_i4b )
!
            smax(:m) = zero
!
        case( 1_i4b )
!
            if ( dim==1_i4b ) then
                smax(:m) = abs( mat(:,1_i4b) )
            else
                smax(:m) = abs( mat(1_i4b,:) )
            end if
!            smax(:m) = merge( abs(mat(1_i4b,:)), abs(mat(:,1_i4b)), dim==2_i4b)
!
        case( 2_i4b: )
!
            cutlo  = sqrt( machsmlnum )
            hitest = one/( cutlo*real( n, stnd ) )
!
            sfmin = lamch( 's' )
            sfmax = one/sfmin
!
#ifdef _OPENMP3
            i           = omp_get_max_threads()
            ompparallel = omp_in_parallel()
            test_par    = .not.ompparallel   .and.      &
                          n*m>=omp_limit     .and.      &
                          m>1_i4b            .and.      &
                          i>1_i4b
#endif
!
!           FIND THE NORM OF EACH ROW OR COLUMN IN mat .
!
            if ( dim==1_i4b ) then
!
#ifdef _OPENMP3
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
!
                do i = 1_i4b, m
!
!                   FIND THE MAXIMUM FOR EACH ROW IN mat .
!
                    smax(i) = maxval( abs(mat(i,:n)) )
!
                    if ( smax(i)>cutlo .and. smax(i)<hitest ) then
!
!                       IF  smax(i) IS GREATER THAN cutlo AND smax(i) IS LESS THAN hitest, 
!                       NO SCALING SHOULD BE NEEDED.
!
                        smax(i) = sqrt( sum( mat(i,:n)**2 ) )
!
                    else
!
!                       COMPUTE THE SUM OF SQUARES OF VALUES FOR EACH ROW
!                       OF mat SCALED BY smax .
!
                        smax(i) = min( max( sfmin, smax(i) ), sfmax ) 
                        smax(i) = sqrt( sum( (mat(i,:n)/smax(i))**2 ) )*smax(i)
!
                    end if
!
                end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!               
            else
!
#ifdef _OPENMP3
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
                do i = 1_i4b, m
!
!                   FIND THE MAXIMUM FOR EACH COLUMN IN mat .
!
                    smax(i) = maxval( abs(mat(:n,i)) )
!
                    if ( smax(i)>cutlo .and. smax(i)<hitest ) then
!
!                       IF  smax(i) IS GREATER THAN cutlo AND smax(i) IS LESS THAN hitest, 
!                       NO SCALING SHOULD BE NEEDED.
!
                        smax(i) = sqrt( sum( mat(:n,i)**2 ) )
!
                    else
!
!                       COMPUTE THE SUM OF SQUARES OF VALUES FOR EACH COLUMN
!                       OF mat SCALED BY smax .
!
                        smax(i) = min( max( sfmin, smax(i) ), sfmax ) 
                        smax(i) = sqrt( sum( (mat(:n,i)/smax(i))**2 ) )*smax(i)
!
                    end if
!
                end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!               
            end if
!
    end select
!
    norm_dim_rm(:m) = smax(:m)
!
    end function norm_dim_rm
!
! =========================================================================================
!
    function norm_dim_cm( mat, dim )
!
! Purpose
! _______
!
!   Return the Euclidean norms of the columns (DIM=2) or of the rows (DIM=1) of a complex 
!   matrix MAT via the function name, so that
!
!        norm_dim_cm := sqrt( sum(MAT * conjg(MAT),dim=3-dim) )
!
!   This is done without destructive underflow or overflow.
!
!
! ___________________________________________________________________________________________
!
    use Reals_Constants,  only : zero, one
    use Num_Constants,    only : lamch, machsmlnum
!
    complex(stnd), dimension(:,:), intent(in) :: mat
    integer(i4b),                  intent(in) :: dim
!
    real(stnd), dimension(size(mat,int(dim))) :: norm_dim_cm
!
    real(stnd)   :: smax(size(mat,int(dim))), hitest, cutlo, sfmin, sfmax
    integer(i4b) :: n, m, dim2, i
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( mat, int(dim) )
!
    if ( m<=0_i4b ) return
!
    dim2 = 3_i4b-dim
    n = size( mat, int(dim2) )
!
    select case( n )
!
        case( :0_i4b )
!
            smax(:m) = zero
!
        case( 1_i4b )
!
            if ( dim==1_i4b ) then
                smax(:m) = abs( mat(:,1_i4b) )
            else
                smax(:m) = abs( mat(1_i4b,:) )
            end if
!            smax(:m) = merge( abs(mat(1_i4b,:)), abs(mat(:,1_i4b)), dim==2_i4b)
!
        case( 2_i4b: )
!
            cutlo  = sqrt( machsmlnum )
            hitest = one/( cutlo*real( 2_i4b*n, stnd ) )
!
            sfmin = lamch( 's' )
            sfmax = one/sfmin
!
#ifdef _OPENMP3
            i           = omp_get_max_threads()
            ompparallel = omp_in_parallel()
            test_par    = .not.ompparallel   .and.      &
                          n*m>=omp_limit     .and.      &
                          m>1_i4b            .and.      &
                          i>1_i4b
#endif
!
!           FIND THE NORM OF EACH ROW OR COLUMN IN mat .
!
            if ( dim==1_i4b ) then
!
#ifdef _OPENMP3
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
                do i = 1_i4b, m
!
!                   FIND THE MAXIMUM FOR EACH ROW IN mat .
!
                    smax(i) = maxval( max( abs(real(mat(i,:n))), abs(aimag(mat(i,:n))) ) ) 
!
                    if ( smax(i)>cutlo .and. smax(i)<hitest ) then
!
!                       IF  smax(i) IS GREATER THAN cutlo AND smax(i) IS LESS THAN hitest, 
!                       NO SCALING SHOULD BE NEEDED.
!
                        smax(i) = sqrt( sum( real(mat(i,:n))**2+aimag(mat(i,:n))**2 ) )
!
                    else
!
!                       COMPUTE THE SUM OF SQUARES OF VALUES FOR EACH ROW
!                       OF mat SCALED BY smax .
!
                        smax(i) = min( max( sfmin, smax(i) ), sfmax ) 
                        smax(i) = sqrt( sum( (real( mat(i,:n))/smax(i))**2  +             &
                                             (aimag(mat(i,:n))/smax(i))**2  ) )*smax(i)
!
                    end if
!
                end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!               
            else
!
#ifdef _OPENMP3
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC), PRIVATE(i)
#endif
                do i = 1_i4b, m
!
!                   FIND THE MAXIMUM FOR EACH COLUMN IN mat .
!
                    smax(i) = maxval( max( abs(real(mat(:n,i))), abs(aimag(mat(:n,i))) ) ) 
!
                    if ( smax(i)>cutlo .and. smax(i)<hitest ) then
!
!                       IF  smax(i) IS GREATER THAN cutlo AND smax(i) IS LESS THAN hitest, 
!                       NO SCALING SHOULD BE NEEDED.
!
                        smax(i) = sqrt( sum( real(mat(:n,i))**2+aimag(mat(:n,i))**2 ) )
!
                    else
!
!                       COMPUTE THE SUM OF SQUARES OF VALUES FOR EACH COLUMN
!                       OF mat SCALED BY smax .
!
                        smax(i) = min( max( sfmin, smax(i) ), sfmax ) 
                        smax(i) = sqrt( sum( (real( mat(:n,i))/smax(i))**2  +             &
                                             (aimag(mat(:n,i))/smax(i))**2  ) )*smax(i)
!
                    end if
!
                end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!               
            end if
!
    end select
!
    norm_dim_cm(:m) = smax(:m)
!
    end function norm_dim_cm
!
!
! =========================================================================================
!                  ROUTINES FOR SCATTER-WITH-COMBINE (MAY BE PARALLELIZED WITH OpenMP).
!                              ADAPTED FROM Numerical Recipes.
! =========================================================================================
!
!
    subroutine scatter_add_i( dest, source, dest_index )
!
! Purpose
! _______
!
!   Adds each component of the integer vector SOURCE into a component of the integer
!   vector DEST specified by the index vector DEST_INDEX.
!
!
! __________________________________________________________________________________________
!
    integer(i4b),   dimension(:), intent(inout) :: dest
    integer(i4b),   dimension(:), intent(in)    :: source
    integer(i4b),   dimension(:), intent(in)    :: dest_index
!
    integer(i4b) :: m, n, j, i
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = assert_eq2( int( size(source), i4b),       &
                    int( size(dest_index), i4b),   &
                    'scatter_add' )
    m = size( dest )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par)    &
!$OMP         ,PRIVATE(i,j)       &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
!
        i = dest_index(j)
        if ( i>0_i4b .and. i<=m ) dest(i) = dest(i) + source(j)
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine scatter_add_i
!
! =========================================================================================
!
    subroutine scatter_add_r( dest, source, dest_index )
!
! Purpose
! _______
!
!   Adds each component of the real vector SOURCE into a component of the real
!   vector DEST specified by the index vector DEST_INDEX.
!   
!
! __________________________________________________________________________________________
!
    real(stnd),   dimension(:), intent(inout) :: dest
    real(stnd),   dimension(:), intent(in)    :: source
    integer(i4b), dimension(:), intent(in)    :: dest_index
!
    integer(i4b) :: m, n, j, i
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = assert_eq2( int( size(source), i4b),       &
                    int( size(dest_index), i4b),   &
                    'scatter_add' )
    m = size( dest )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par)     &
!$OMP         ,PRIVATE(i,j)        &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
!
        i = dest_index(j)
        if ( i>0_i4b .and. i<=m ) dest(i) = dest(i) + source(j)
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine scatter_add_r
!
! =========================================================================================
!
    subroutine scatter_add_c( dest, source, dest_index )
!
! Purpose
! _______
!
!   Adds each component of the complex vector SOURCE into a component of the complex
!   vector DEST specified by the index vector DEST_INDEX.
!   
!
! __________________________________________________________________________________________
!
    complex(stnd),   dimension(:), intent(inout) :: dest
    complex(stnd),   dimension(:), intent(in)    :: source
    integer(i4b),    dimension(:), intent(in)    :: dest_index
!
    integer(i4b) :: m, n, j, i
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = assert_eq2( int( size(source), i4b),       &
                    int( size(dest_index), i4b),   &
                    'scatter_add' )
    m = size( dest )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par)    &
!$OMP         ,PRIVATE(i,j)       &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
!
        i = dest_index(j)
        if ( i>0_i4b .and. i<=m ) dest(i) = dest(i) + source(j)
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine scatter_add_c
!
! =========================================================================================
!
    subroutine scatter_max_r( dest, source, dest_index )
!
! Purpose
! _______
!
!   Takes the max operation between each component of the real vector SOURCE and
!   a component of the real vector DEST specified by the index vector DEST_INDEX,
!   replacing the component of DEST with the value obtained.
!   
!
! __________________________________________________________________________________________
!
    real(stnd),   dimension(:), intent(inout) :: dest
    real(stnd),   dimension(:), intent(in)    :: source
    integer(i4b), dimension(:), intent(in)    :: dest_index
!
    integer(i4b) :: m, n, j, i
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = assert_eq2( int( size(source), i4b),       &
                    int( size(dest_index), i4b),   &
                    'scatter_max' )
    m = size( dest )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(i,j)           &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
!
        i = dest_index(j)
        if ( i>0_i4b .and. i<=m ) dest(i) = max( dest(i), source(j) )
!
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine scatter_max_r
!
! =========================================================================================
!
    subroutine scatter_max_i( dest, source, dest_index )
!
! Purpose
! _______
!
!   Takes the max operation between each component of the integer vector SOURCE and
!   a component of the integer vector DEST specified by the index vector DEST_INDEX,
!   replacing the component of DEST with the value obtained.
!   
!
! __________________________________________________________________________________________
!
    integer(i4b),   dimension(:), intent(inout) :: dest
    integer(i4b),   dimension(:), intent(in)    :: source
    integer(i4b),   dimension(:), intent(in)    :: dest_index
!
    integer(i4b) :: m, n, j, i
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = assert_eq2( int( size(source), i4b),       &
                    int( size(dest_index), i4b),   &
                    'scatter_max' )
    m = size( dest )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(i,j)           &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        i = dest_index(j)
        if ( i>0_i4b .and. i<=m ) dest(i) = max( dest(i), source(j) )
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine scatter_max_i
!
!
! =========================================================================================
!            ROUTINES FOR SKEW OPERATIONS ON MATRICES (PARALLELIZED WITH OpenMP).
!                      ADAPTED AND EXTENDED FROM Numerical Recipes.
! =========================================================================================
!
!
    subroutine diagadd_rv( mat, diag )
!
! Purpose
! _______
!
!   Adds real vector DIAG to the diagonal of real matrix MAT.
!
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:,:), intent(inout) :: mat
    real(stnd), dimension(:),   intent(in)    :: diag
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = assert_eq2( int( size(diag), i4b),                       &
                    int( min(size(mat,1), size(mat,2)), i4b),    &
                    'diagadd' )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = mat(j,j) + diag(j)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine diagadd_rv
!
! =========================================================================================
!
    subroutine diagadd_cv( mat, diag )
!
! Purpose
! _______
!
!   Adds complex vector DIAG to the diagonal of complex matrix MAT.
!
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:,:), intent(inout) :: mat
    complex(stnd), dimension(:),   intent(in)    :: diag
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = assert_eq2( int( size(diag), i4b),                       &
                    int( min(size(mat,1), size(mat,2)), i4b),    &
                    'diagadd' )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = mat(j,j) + diag(j)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine diagadd_cv
!
! =========================================================================================
!
    subroutine diagadd_r( mat, diag )
!
! Purpose
! _______
!
!   Adds real scalar DIAG to the diagonal of real matrix MAT.
!
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:,:), intent(inout) :: mat
    real(stnd),                 intent(in)    :: diag
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = min( size(mat,1), size(mat,2) )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = mat(j,j) + diag
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine diagadd_r
!
! =========================================================================================
!
    subroutine diagadd_c( mat, diag )
!
! Purpose
! _______
!
!   Adds complex scalar DIAG to the diagonal of complex matrix MAT.
!
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:,:), intent(inout) :: mat
    complex(stnd),                 intent(in)    :: diag
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = min( size(mat,1), size(mat,2) )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = mat(j,j) + diag
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine diagadd_c
!
! =========================================================================================
!
    subroutine diagmult_rv( mat, diag )
!
! Purpose
! _______
!
!   Multiplies real vector DIAG into the diagonal of real matrix MAT.
!
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:,:), intent(inout) :: mat
    real(stnd), dimension(:),   intent(in)    :: diag
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = assert_eq2( int( size(diag), i4b),                       &
                    int( min(size(mat,1), size(mat,2)), i4b),    &
                    'diagmult' )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = mat(j,j)*diag(j)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine diagmult_rv
!
! =========================================================================================
!
    subroutine diagmult_cv( mat, diag )
!
! Purpose
! _______
!
!   Multiplies complex vector DIAG into the diagonal of complex matrix MAT.
!
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:,:), intent(inout) :: mat
    complex(stnd), dimension(:),   intent(in)    :: diag
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = assert_eq2( int( size(diag), i4b),                       &
                    int( min(size(mat,1), size(mat,2)), i4b),    &
                    'diagmult' )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = mat(j,j)*diag(j)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine diagmult_cv
!
! =========================================================================================
!
    subroutine diagmult_r( mat, diag )
!
! Purpose
! _______
!
!   Multiplies real scalar DIAG into the diagonal of real matrix MAT.
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:,:), intent(inout) :: mat
    real(stnd),                 intent(in)    :: diag
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = min( size(mat,1), size(mat,2) )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = mat(j,j)*diag
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine diagmult_r
!
! =========================================================================================
!
    subroutine diagmult_c( mat, diag )
!
! Purpose
! _______
!
!   Multiplies complex scalar DIAG into the diagonal of complex matrix MAT.
!
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:,:), intent(inout) :: mat
    complex(stnd),                 intent(in)    :: diag
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = min( size(mat,1), size(mat,2) )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = mat(j,j)*diag
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine diagmult_c
!
! =========================================================================================
!
    function get_diag_r( mat )
!
! Purpose
! _______
!
!   Returns as a vector the diagonal of real matrix MAT.
!
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:,:), intent(in) :: mat
!
    real(stnd), dimension(min(size(mat,1),size(mat,2))) :: get_diag_r
    real(stnd), dimension(min(size(mat,1),size(mat,2))) :: tmp
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = min( size(mat,1), size(mat,2) )

#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        tmp(j) = mat(j,j)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    get_diag_r = tmp
!
    end function get_diag_r
!
! =========================================================================================
!
    function get_diag_c( mat )
!
! Purpose
! _______
!
!   Returns as a vector the diagonal of complex matrix MAT.
!
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:,:), intent(in) :: mat
!
    complex(stnd), dimension(min(size(mat,1),size(mat,2))) :: get_diag_c
    complex(stnd), dimension(min(size(mat,1),size(mat,2))) :: tmp
    integer(i4b)  :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = min( size(mat,1), size(mat,2) )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        tmp(j) = mat(j,j)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    get_diag_c = tmp
!
    end function get_diag_c
!
! =========================================================================================
!
    subroutine put_diag_rv( diag, mat )
!
! Purpose
! _______
!
!   Set the diagonal of real matrix MAT to the values of the real vector DIAG.
!
!
! __________________________________________________________________________________________
!
    real(stnd), dimension(:),   intent(in)    :: diag
    real(stnd), dimension(:,:), intent(inout) :: mat
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = assert_eq2( int( size(diag), i4b),                       &
                    int( min(size(mat,1), size(mat,2)), i4b),    &
                    'put_diag' )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = diag(j)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine put_diag_rv
!
! =========================================================================================
!
    subroutine put_diag_cv( diag, mat )
!
! Purpose
! _______
!
!   Set the diagonal of complex matrix MAT to the values of the complex vector DIAG.
!
!
! __________________________________________________________________________________________
!
    complex(stnd), dimension(:),   intent(in)    :: diag
    complex(stnd), dimension(:,:), intent(inout) :: mat
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = assert_eq2( int( size(diag), i4b),                       &
                    int( min(size(mat,1), size(mat,2)), i4b),    &
                    'put_diag' )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = diag(j)
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine put_diag_cv
!
! =========================================================================================
!
    subroutine put_diag_r( diag, mat )
!
! Purpose
! _______
!
!   Set the diagonal of real matrix MAT to the value of the real scalar DIAG.
!
!
! __________________________________________________________________________________________
!
    real(stnd),                 intent(in)    :: diag
    real(stnd), dimension(:,:), intent(inout) :: mat
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = min( size(mat,1), size(mat,2) )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = diag
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine put_diag_r
!
! =========================================================================================
!
    subroutine put_diag_c( diag, mat )
!
! Purpose
! _______
!
!   Set the diagonal of complex matrix MAT to the value of the complex scalar DIAG.
!
!
! __________________________________________________________________________________________
!
    complex(stnd),                 intent(in)    :: diag
    complex(stnd), dimension(:,:), intent(inout) :: mat
!
    integer(i4b) :: j, n
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    n = min( size(mat,1), size(mat,2) )
!
#ifdef _OPENMP3
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  n>=omp_limit       .and.      &
                  j>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(j)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do j = 1_i4b, n
        mat(j,j) = diag
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    end subroutine put_diag_c
!
! =========================================================================================
!
    subroutine unit_matrix_r( mat )
!
! Purpose
! _______
!
!   Set the real matrix MAT to be a unit real matrix (if it is square).
!
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : one, zero
!
    real(stnd), dimension(:,:), intent(out) :: mat
!
    integer(i4b) :: i, n, m, nm
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( mat, 1 )
    n = size( mat, 2 )
    nm = min( m, n )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  m*nm>=omp_limit    .and.      &
                  nm>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(i)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do i = 1_i4b, nm
        mat(:m,i) = zero
        mat(i,i) = one
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    mat(:m,nm+1_i4b:n) = zero
!
    end subroutine unit_matrix_r
!
! =========================================================================================
!
    subroutine unit_matrix_c( mat )
!
! Purpose
! _______
!
!   Set the complex matrix MAT to be a unit complex matrix (if it is square).
!
!
! __________________________________________________________________________________________
!
    use Reals_Constants,  only : one, zero
!
    complex(stnd), dimension(:,:), intent(out) :: mat
!
    integer(i4b) :: i, n, m, nm
#ifdef _OPENMP3
    logical :: test_par, ompparallel
#endif
!
    m = size( mat, 1 )
    n = size( mat, 2 )
    nm = min( m, n )
!
#ifdef _OPENMP3
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    test_par    = .not.ompparallel   .and.      &
                  m*nm>=omp_limit    .and.      &
                  nm>1_i4b           .and.      &
                  i>1_i4b
!
!$OMP PARALLEL DO IF(test_par)        &
!$OMP         ,PRIVATE(i)             &
!$OMP         ,SCHEDULE(STATIC)
#endif
!
    do i = 1_i4b, nm
        mat(:m,i) = zero
        mat(i,i) = one
    end do
!
#ifdef _OPENMP3
!$OMP END PARALLEL DO
#endif
!
    mat(:m,nm+1_i4b:n) = zero
!
    end subroutine unit_matrix_c
!
!
! =========================================================================================
!                     ROUTINES FOR SCALING MATRICES, VECTORS AND SCALARS.
! =========================================================================================
!
!
    subroutine lascl_r( x, cfrom, cto )
!
! Purpose
! _______
!
!   LASCL multiplies the real scalar  X  by the real scalar  CTO/CFROM .
!   This is done without over/underflow as long as the final result
!   CTO * X/CFROM  does not over/underflow. 
!   CFROM  must be nonzero.
!
!
! Arguments
! _________
!
!   X               (INPUT/OUTPUT) real(stnd)
!                   The real to be multiplied by CTO/CFROM.
!
!   CFROM, CTO      (INPUT) real(stnd)
!                   The real X is multiplied by CTO/CFROM.
!
!
! Further Details
! _______________
!
!   This subroutine is adapted from the routine DLASCL in LAPACK77 (version 3) with
!   improvements suggested by E. Anderson. See:
!
!   (1) Anderson, E., 2002:
!           LAPACK3E -- A Fortran90-enhanced version of LAPACK.
!           Lapack Working Note 158, University of Tennessee.
!
!
! ___________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Num_Constants,    only : lamch
    use Reals_Constants,  only : one, zero
    use Char_Constants,   only : utilities_error4
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout) :: x
    real(stnd), intent(in)    :: cfrom, cto
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: cfromc, ctoc, cfrom1, cto1, mul, smlnum, bignum
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE INPUT ARGUMENTS.
!
    if ( cfrom==zero ) call merror( 'lascl'//utilities_error4 )
!
!   GET MACHINE PARAMETERS.
!
    smlnum = sqrt( lamch( 's' ) )
    bignum = one/smlnum
!
    cfromc = cfrom
    ctoc   = cto
!
    cfrom1 = cfromc*smlnum
    cto1   = ctoc/bignum
!
    if ( abs(cfrom1)>abs(ctoc) .and. ctoc/=zero ) then
!
!       PRE-MULTIPLY x BY smlnum IF cfrom IS LARGE COMPARED TO cto .
!
        x      = x*smlnum
        cfromc = cfrom1
!
    else if ( abs(cto1)>abs(cfromc) ) then
!
!       PRE-MULTIPLY x BY bignum IF cfrom IS SMALL COMPARED TO cto .
!
        x    = x*bignum
        ctoc = cto1
!
    end if
!
!   MULTIPLY x BY cto/cfrom .
!
    mul = ctoc/cfromc
    x   = x*mul
!
!
! END OF SUBROUTINE lascl_r
! _________________________
!
    end subroutine lascl_r
!
! =========================================================================================
!
    subroutine lascl_rv( x, cfrom, cto )
!
! Purpose
! _______
!
!   LASCL multiplies the real vector  X  by the real scalar  CTO/CFROM .
!   This is done without over/underflow as long as the final result
!   CTO * X(i)/CFROM  does not over/underflow for i = 1 to size( X). 
!
!   CFROM  must be nonzero.
!
!
! Arguments
! _________
!
!   X               (INPUT/OUTPUT) real(stnd), dimension(:)
!                   The real vector to be multiplied by CTO/CFROM.
!
!   CFROM, CTO      (INPUT) real(stnd)
!                   The real vector X is multiplied by CTO/CFROM.
!
!
! Further Details
! _______________
!
!   This subroutine is adapted from the routine DLASCL in LAPACK77 (version 3) with
!   improvements suggested by E. Anderson. See:
!
!   (1) Anderson, E., 2002:
!           LAPACK3E -- A Fortran90-enhanced version of LAPACK.
!           Lapack Working Note 158, University of Tennessee.
!
!
! ___________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Num_Constants,    only : lamch
    use Reals_Constants,  only : one, zero
    use Char_Constants,   only : utilities_error4
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout) :: x(:)
    real(stnd), intent(in)    :: cfrom, cto
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: cfromc, ctoc, cfrom1, cto1, mul, smlnum, bignum
!
    integer(i4b) :: n
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE INPUT ARGUMENTS.
!
    if ( cfrom==zero ) call merror( 'lascl'//utilities_error4 )
!
!   QUICK RETURN IF POSSIBLE.
!
    n = size( x )
    if ( n==0_i4b ) return
!
!   GET MACHINE PARAMETERS.
!
    smlnum = sqrt( lamch( 's' ) )
    bignum = one/smlnum
!
    cfromc = cfrom
    ctoc   = cto
!
    cfrom1 = cfromc*smlnum
    cto1   = ctoc/bignum
!
    if ( abs(cfrom1)>abs(ctoc) .and. ctoc/=zero ) then
!
!       PRE-MULTIPLY x BY smlnum IF cfrom IS LARGE COMPARED TO cto .
!
        x(:n)  = x(:n)*smlnum
        cfromc = cfrom1
!
    else if ( abs(cto1)>abs(cfromc) ) then
!
!       PRE-MULTIPLY x BY bignum IF cfrom IS SMALL COMPARED TO cto .
!
        x(:n) = x(:n)*bignum
        ctoc  = cto1
!
    end if
!
!   MULTIPLY x BY cto/cfrom .
!
    mul   = ctoc/cfromc
    x(:n) = x(:n)*mul
!
!
! END OF SUBROUTINE lascl_rv
! __________________________
!
    end subroutine lascl_rv
!
! =========================================================================================
!
    subroutine lascl_rm( x, cfrom, cto )
!
! Purpose
! _______
!
!   LASCL multiplies the real matrix  X  by the real scalar  CTO/CFROM .
!   This is done without over/underflow as long as the final result
!   CTO * X(i,j)/CFROM  does not over/underflow for i = 1 to size( X, 1) and
!   j = 1 to size( X, 2). 
!
!   CFROM  must be nonzero.
!
!
! Arguments
! _________
!
!   X               (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                   The real matrix to be multiplied by CTO/CFROM.
!
!   CFROM, CTO      (INPUT) real(stnd)
!                   The real matrix X is multiplied by CTO/CFROM.
!
!
! Further Details
! _______________
!
!   This subroutine is adapted from the routine DLASCL in LAPACK77 (version 3) with
!   improvements suggested by E. Anderson. See:
!
!   (1) Anderson, E., 2002:
!           LAPACK3E -- A Fortran90-enhanced version of LAPACK.
!           Lapack Working Note 158, University of Tennessee.
!
!
! ___________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Num_Constants,    only : lamch
    use Reals_Constants,  only : one, zero
    use Char_Constants,   only : utilities_error4
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout) :: x(:,:)
    real(stnd), intent(in)    :: cfrom, cto
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: cfromc, ctoc, cfrom1, cto1, mul, smlnum, bignum
!
    integer(i4b) :: m, n
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE INPUT ARGUMENTS.
!
    if ( cfrom==zero ) call merror( 'lascl'//utilities_error4 )
    m = size( x, dim=1 )
    n = size( x, dim=2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<=0_i4b .or. m<=0_i4b ) return
!
!   GET MACHINE PARAMETERS.
!
    smlnum = sqrt( lamch( 's' ) )
    bignum = one/smlnum
!
    cfromc = cfrom
    ctoc   = cto
!
    cfrom1 = cfromc*smlnum
    cto1   = ctoc/bignum
!
    if ( abs(cfrom1)>abs(ctoc) .and. ctoc/=zero ) then
!
!       PRE-MULTIPLY x BY smlnum IF cfrom IS LARGE COMPARED TO cto .
!
        x(:m,:n)  = x(:m,:n)*smlnum
        cfromc    = cfrom1
!
    else if ( abs(cto1)>abs(cfromc) ) then
!
!       PRE-MULTIPLY x BY bignum IF cfrom IS SMALL COMPARED TO cto .
!
        x(:m,:n) = x(:m,:n)*bignum
        ctoc     = cto1
!
    end if
!
!   MULTIPLY x BY cto/cfrom .
!
    mul      = ctoc/cfromc
    x(:m,:n) = x(:m,:n)*mul
!
!
! END OF SUBROUTINE lascl_rm
! __________________________
!
    end subroutine lascl_rm
!
! =========================================================================================
!
    subroutine type_lascl_rm( x, cfrom, cto, type )
!
! Purpose
! _______
!
!   LASCL multiplies the real matrix  X  by the real scalar  CTO/CFROM .
!   This is done without over/underflow as long as the final result
!   CTO * X(i,j)/CFROM  does not over/underflow for i = 1 to size( X, 1) and
!   j = 1 to size( X, 2). 
!
!   CFROM  must be nonzero.
!
!   TYPE specifies that X may be full, upper triangular, lower triangular or
!   upper Hessenberg.
!
!
! Arguments
! _________
!
!   X               (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                   The real matrix to be multiplied by CTO/CFROM.
!
!   CFROM, CTO      (INPUT) real(stnd)
!                   The real matrix X is multiplied by CTO/CFROM.
!
!   TYPE            (INPUT) character*1
!                   TYPE indices the storage type of the input matrix.
!                   = 'L' or 'l': X is a lower triangular matrix.
!                   = 'U' or 'u': X is a upper triangular matrix.
!                   = 'H' or 'h': X is a upper Hessenberg matrix.
!                   = 'G' or 'g': X is a full matrix.
!                   = any other character: X is assumed to be a full matrix.
!
!
! Further Details
! _______________
!
!   This subroutine is adapted from the routine DLASCL in LAPACK77 (version 3) with
!   improvements suggested by E. Anderson. See:
!
!   (1) Anderson, E., 2002:
!           LAPACK3E -- A Fortran90-enhanced version of LAPACK.
!           Lapack Working Note 158, University of Tennessee.
!
!
! ___________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Num_Constants,    only : lamch
    use Reals_Constants,  only : one, zero
    use Char_Constants,   only : utilities_error4
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout) :: x(:,:)
    real(stnd), intent(in)    :: cfrom, cto
!
    character, intent(in) :: type
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: cfromc, ctoc, cfrom1, cto1, mul, smlnum, bignum
!
    integer(i4b) :: k, j, m, n
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE INPUT ARGUMENTS.
!
    if ( cfrom==zero ) call merror( 'lascl'//utilities_error4 )
!
    n = size( x, 2 )
    m = size( x, 1 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<=0_i4b .or. m<=0_i4b ) return
!
!   GET MACHINE PARAMETERS.
!
    smlnum = sqrt( lamch( 's' ) )
    bignum = one/smlnum
!
    cfromc = cfrom
    ctoc   = cto
!
    cfrom1 = cfromc*smlnum
    cto1   = ctoc/bignum
!
    if ( abs(cfrom1)>abs(ctoc) .and. ctoc/=zero ) then
!
!       PRE-MULTIPLY x BY smlnum IF cfrom IS LARGE COMPARED TO cto .
!
        select case ( type )
!
            case( achar(76), achar(108) )
!
!               LOWER CASE.
!
                do j = 1_i4b, min(n,m)
                    x(j:m,j) = smlnum*x(j:m,j)
                end do
!
            case( achar(85), achar(117) )
!
!               UPPER CASE.
!
                do j = 1_i4b, n
                    k = min( j, m )
                    x(:k,j) = smlnum*x(:k,j)
                end do
!
            case( achar(72), achar(104) )
!
!               HESSENBERG CASE.
!
                do j = 1_i4b, n
                    k = min( j+1_i4b, m )
                    x(:k,j) = smlnum*x(:k,j)
                end do
!
            case default
!
!               GENERAL CASE.
!
                x(:m,:n) = smlnum*x(:m,:n)
!
        end select
!
        cfromc = cfrom1
!
    else if ( abs(cto1)>abs(cfromc) ) then
!
!       PRE-MULTIPLY x BY bignum IF cfrom IS SMALL COMPARED TO cto .
!
        select case ( type )
!
            case( achar(76), achar(108) )
!
!               LOWER CASE.
!
                do j = 1_i4b, min(n,m)
                    x(j:m,j) = bignum*x(j:m,j)
                end do
!
            case( achar(85), achar(117) )
!
!               UPPER CASE.
!
                do j = 1_i4b, n
                    k = min( j, m )
                    x(:k,j) = bignum*x(:k,j)
                end do
!
            case( achar(72), achar(104) )
!
!               HESSENBERG CASE.
!
                do j = 1_i4b, n
                    k = min( j+1_i4b, m )
                    x(:k,j) = bignum*x(:k,j)
                end do
!
            case default
!
!               GENERAL CASE.
!
                x(:m,:n) = bignum*x(:m,:n)
!
        end select
!
        ctoc = cto1
!
    end if
!
!   MULTIPLY x BY cto/cfrom .
!
    mul = ctoc/cfromc
!
    select case ( type )
!
        case( achar(76), achar(108) )
!
!           LOWER CASE.
!
            do j = 1_i4b, min(n,m)
                x(j:m,j) = mul*x(j:m,j)
            end do
!
        case( achar(85), achar(117) )
!
!           UPPER CASE.
!
            do j = 1_i4b, n
                k = min( j, m )
                x(:k,j) = mul*x(:k,j)
            end do
!
        case( achar(72), achar(104) )
!
!           HESSENBERG CASE.
!
            do j = 1_i4b, n
                k = min( j+1_i4b, m )
                x(:k,j) = mul*x(:k,j)
            end do
!
        case default
!
!           GENERAL CASE.
!
            x(:m,:n) = mul*x(:m,:n)
!
    end select
!
!
! END OF SUBROUTINE type_lascl_rm
! _______________________________
!
    end subroutine type_lascl_rm
!
! =========================================================================================
!
    subroutine masked_lascl_r( x, cfrom, cto, mask )
!
! Purpose
! _______
!
!   LASCL multiplies the real scalar  X  by the real scalar  CTO/CFROM
!   under the control of the logical argument MASK .
!   This is done without over/underflow as long as the final result
!   CTO * X/CFROM  does not over/underflow. 
!
!   CFROM  must be nonzero.
!
!
! Arguments
! _________
!
!   X               (INPUT/OUTPUT) real(stnd)
!                   The real to be multiplied by CTO/CFROM.
!
!   CFROM, CTO      (INPUT) real(stnd)
!                   The real X is multiplied by CTO/CFROM
!                   if MASK=true.
!
!   MASK            (INPUT) logical(lgl)
!                   The logical mask : if MASK=true
!                   the multiplication is done, otherwise
!                   X is left unchanged.
!
!
! Further Details
! _______________
!
!   This subroutine is adapted from the routine DLASCL in LAPACK77 (version 3) with
!   improvements suggested by E. Anderson. See:
!
!   (1) Anderson, E., 2002:
!           LAPACK3E -- A Fortran90-enhanced version of LAPACK.
!           Lapack Working Note 158, University of Tennessee.
!
!
! ___________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Num_Constants,    only : lamch
    use Reals_Constants,  only : one, zero
    use Char_Constants,   only : utilities_error4
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout) :: x
    real(stnd),   intent(in)    :: cfrom, cto
!
    logical(lgl), intent(in)    :: mask
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: cfromc, ctoc, cfrom1, cto1, mul, smlnum, bignum
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE INPUT ARGUMENTS.
!
    if ( cfrom==zero ) call merror( 'lascl'//utilities_error4 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( .not.mask ) return
!
!   GET MACHINE PARAMETERS.
!
    smlnum = sqrt( lamch( 's' ) )
    bignum = one/smlnum
!
    cfromc = cfrom
    ctoc   = cto
!
    cfrom1 = cfromc*smlnum
    cto1   = ctoc/bignum
!
    if ( abs(cfrom1)>abs(ctoc) .and. ctoc/=zero ) then
!
!       PRE-MULTIPLY x BY smlnum IF cfrom IS LARGE COMPARED TO cto .
!
        x      = x*smlnum
        cfromc = cfrom1
!
    else if ( abs(cto1)>abs(cfromc) ) then
!
!       PRE-MULTIPLY x BY bignum IF cfrom IS SMALL COMPARED TO cto .
!
        x    = x*bignum
        ctoc = cto1
!
    end if
!
!   MULTIPLY x BY cto/cfrom .
!
    mul = ctoc/cfromc
    x   = x*mul
!
!
! END OF SUBROUTINE masked_lascl_r
! ________________________________
!
    end subroutine masked_lascl_r
!
! =========================================================================================
!
    subroutine masked_lascl_rv( x, cfrom, cto, mask )
!
! Purpose
! _______
!
!   LASCL multiplies the real vector  X  by the real scalar  CTO/CFROM
!   under the control of the logical argument MASK .
!   This is done without over/underflow as long as the final result
!   CTO * X(i)/CFROM  does not over/underflow for i = 1 to size( X). 
!
!   CFROM  must be nonzero.
!
!
! Arguments
! _________
!
!   X               (INPUT/OUTPUT) real(stnd), dimension(:)
!                   The real vector to be multiplied by CTO/CFROM.
!
!   CFROM, CTO      (INPUT) real(stnd)
!                   The real X(i) is multiplied by CTO/CFROM if MASK(i)=true.
!
!   MASK            (INPUT) logical(lgl), dimension(:)
!                   The logical mask : if MASK(i)=true
!                   the multiplication is done, otherwise
!                   X(i) is left unchanged.
!
!
! Further Details
! _______________
!
!   This subroutine is adapted from the routine DLASCL in LAPACK77 (version 3) with
!   improvements suggested in reference (1).
!
!   The sizes of X and MASK must match.
!
!   For further details, see:
!
!   (1) Anderson, E., 2002:
!           LAPACK3E -- A Fortran90-enhanced version of LAPACK.
!           Lapack Working Note 158, University of Tennessee.
!
!
! ___________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Num_Constants,    only : lamch
    use Reals_Constants,  only : one, zero
    use Char_Constants,   only : utilities_error4
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout) :: x(:)
    real(stnd),   intent(in)    :: cfrom, cto
!
    logical(lgl), intent(in)    :: mask(:)
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: cfromc, ctoc, cfrom1, cto1, mul, smlnum, bignum
!
    integer(i4b) :: n
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='lascl'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE INPUT ARGUMENTS.
!
    if ( cfrom==zero ) call merror( name_proc//utilities_error4 )
    n =  assert_eq2( int(size(x),i4b)    ,   &
                     int(size(mask),i4b) ,   &
                     name_proc  )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<=0_i4b .or. .not.any(mask) ) return
!
!   GET MACHINE PARAMETERS.
!
    smlnum = sqrt( lamch( 's' ) )
    bignum = one/smlnum
!
    cfromc = cfrom
    ctoc   = cto
!
    cfrom1 = cfromc*smlnum
    cto1   = ctoc/bignum
!
    if ( abs(cfrom1)>abs(ctoc) .and. ctoc/=zero ) then
!
!       PRE-MULTIPLY x BY smlnum IF cfrom IS LARGE COMPARED TO cto .
!
        where( mask(:n) )  x(:n)  = x(:n)*smlnum
        cfromc = cfrom1
!
    else if ( abs(cto1)>abs(cfromc) ) then
!
!       PRE-MULTIPLY x BY bignum IF cfrom IS SMALL COMPARED TO cto .
!
        where( mask(:n) )  x(:n) = x(:n)*bignum
        ctoc  = cto1
!
    end if
!
!   MULTIPLY x BY cto/cfrom .
!
    mul = ctoc/cfromc
    where( mask(:n) )  x(:n) = x(:n)*mul
!
!
! END OF SUBROUTINE masked_lascl_rv
! _________________________________
!
    end subroutine masked_lascl_rv
!
! =========================================================================================
!
    subroutine masked_lascl_rm( x, cfrom, cto, mask )
!
! Purpose
! _______
!
!   LASCL multiplies the real matrix  X  by the real scalar  CTO/CFROM 
!   under the control of the logical argument MASK .
!   This is done without over/underflow as long as the final result
!   CTO * X(i,j)/CFROM  does not over/underflow for i = 1 to size( X, 1) and
!   j = 1 to size( X, 2). 
!
!   CFROM  must be nonzero.
!
!
! Arguments
! _________
!
!   X               (INPUT/OUTPUT) real(stnd), dimension(:,:)
!                   The real matrix to be multiplied by CTO/CFROM.
!
!   CFROM, CTO      (INPUT) real(stnd)
!                   The real X(i,j) is multiplied by CTO/CFROM if MASK(i,j)=true.
!
!   MASK            (INPUT) logical(lgl), dimension(:,:)
!                   The logical mask : if MASK(i,j)=true
!                   the multiplication is done, otherwise
!                   X(i,j) is left unchanged.
!
!
! Further Details
! _______________
!
!   This subroutine is adapted from the routine DLASCL in LAPACK77 (version 3) with
!   improvements suggested in reference (1).
!
!   The sizes of X and MASK must match.
!
!   For further details, see:
!
!   (1) Anderson, E., 2002:
!           LAPACK3E -- A Fortran90-enhanced version of LAPACK.
!           Lapack Working Note 158, University of Tennessee.
!
!
! ___________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Num_Constants,    only : lamch
    use Reals_Constants,  only : one, zero
    use Char_Constants,   only : utilities_error4
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout) :: x(:,:)
    real(stnd),   intent(in)    :: cfrom, cto
!
    logical(lgl), intent(in)    :: mask(:,:)
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: cfromc, ctoc, cfrom1, cto1, mul, smlnum, bignum
!
    integer(i4b) :: m, n
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='lascl'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE INPUT ARGUMENTS.
!
    if ( cfrom==zero ) call merror( name_proc//utilities_error4 )
    m =  assert_eq2( int(size(x,1),i4b)    ,   &
                     int(size(mask,1),i4b) ,   &
                     name_proc  )
    n =  assert_eq2( int(size(x,2),i4b)    ,   &
                     int(size(mask,2),i4b) ,   &
                     name_proc  )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<=0_i4b .or. m<=0_i4b .or. .not.any(mask) ) return
!
!   GET MACHINE PARAMETERS.
!
    smlnum = sqrt( lamch( 's' ) )
    bignum = one/smlnum
!
    cfromc = cfrom
    ctoc   = cto
!
    cfrom1 = cfromc*smlnum
    cto1   = ctoc/bignum
!
    if ( abs(cfrom1)>abs(ctoc) .and. ctoc/=zero ) then
!
!       PRE-MULTIPLY x BY smlnum IF cfrom IS LARGE COMPARED TO cto .
!
        where( mask(:m,:n) )  x(:m,:n)  = x(:m,:n)*smlnum
        cfromc = cfrom1
!
    else if ( abs(cto1)>abs(cfromc) ) then
!
!       PRE-MULTIPLY x BY bignum IF cfrom IS SMALL COMPARED TO cto .
!
        where( mask(:m,:n) )  x(:m,:n) = x(:m,:n)*bignum
        ctoc = cto1
!
    end if
!
!   MULTIPLY x BY cto/cfrom .
!
    mul = ctoc/cfromc
    where( mask(:m,:n) )  x(:m,:n) = x(:m,:n)*mul
!
!
! END OF SUBROUTINE masked_lascl_rm
! _________________________________
!
    end subroutine masked_lascl_rm
!
!
! =========================================================================================
!    OTHER ROUTINES FOR COMPUTING THE 2-NORM WITHOUT DESTRUCTIVE UNDERFLOW OR OVERFLOW
! =========================================================================================
!
!
    function norme_rv( vec )
!
! Purpose
! _______
!
!   This function computes the 2-norm (i.e. the Euclidean norm) of the
!   vector VEC of length n, with due regard to avoiding overflow and underflow.
!
!
! Arguments
! _________
!
!   VEC   (INPUT) real(stnd), dimension(:)
!         On entry, the real vector VEC.
!
!
! Further Details
! _______________
!
!   The routine is based on snrm2 from the blas (in linpack), but this
!   version is written in Fortran 90. It is machine independent. The
!   algorithm is described in reference (1).
!
!   The machine constants MACHTINY (the smallest magnitude), MACHBASE(base
!   of the machine), and MACHEPS (epsilon) are used to calculate the
!   constants cutlo and cuthi:
!
!       cutlo = sqrt( machsmlnum ) = sqrt( MACHTINY/(MACHEPS * MACHBASE) )
!       cuthi = one/cutlo
!
!   Three different cases must be considered when calculating the norm:
!
!       (1)  All components of VEC are below cutlo.
!
!            To avoid underflow, each component is divided by sqrt(min)/n
!            and then the regular Euclidean norm of this modified vector
!            is calculated.  This result is then multiplied by
!            sqrt(min)/n  in order to get the correct value for the norm.
!
!       (2)  One or more components are greater than cuthi.
!
!            To avoid overflow, the same method as in case (1) is used
!            with a scaling factor of   sqrt(max) * n .
!
!       (3)  All components are less than cuthi, with at least one
!            component greater than cutlo.
!
!            The regular formula for the Euclidean norm is used.
!
!   For more details, see:
!
!   (1) Blue, J.L., 1978:
!          A portable Fortran program to find the Euclidean norm of a vector.
!          ACM Trans. Math. Soft., Vol. 4, No 1, 1978, pp.15-23.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Reals_Constants,  only : zero, one
    use Num_Constants,    only : machsmlnum
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:), intent(in) :: vec
    real(stnd)                           :: norme_rv
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: hitest, summ, xmax, vi, cutlo
!
    integer(i4b) :: i, n
    integer(i2b) :: rang
!
!
! PARAMETERS
! __________
!
    integer(i2b), parameter ::   null   = 0,     &
                                 small  = 1,     &
                                 normal = 2,     &
                                 large  = 3
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    n = size( vec )
!
    if ( n<=0_i4b ) then
        norme_rv = zero
        return
    end if
!
    cutlo  = sqrt( machsmlnum )
    hitest = one/( cutlo*real( n, stnd ) )
!
    summ = zero
    rang = null
!
!   EVALUATE THE NORM BY ACCUMULATING A SCALED SUM OF SQUARES AND
!   ADJUSTING THE SCALING AS NUMBERS OF INCREASINGLY LARGE MAGNITUDE
!   ARE FOUND.
!
    do  i = 1_i4b, n
!
        vi = abs( vec(i) )
!
        select case( rang )
!
            case( normal )
!
                if ( vi<hitest ) then
                    summ = summ + vi*vi
                else
                    rang = large
                    xmax = vi
                    summ = one + (summ/vi)/vi
                end if
!
            case( small )
!
                if ( vi<=cutlo ) then
                    if ( vi<=xmax ) then
                        summ = summ + (vi/xmax)**2
                    else
                        summ = one + (xmax/vi)**2
                        xmax = vi
                    end if
                else if ( vi>=hitest ) then
                    rang = large
                    xmax = vi
                    summ = one + (summ/vi)/vi
                else
                    rang = normal
                    summ = (summ*xmax)*xmax + vi*vi
                end if
!
            case( large )
!
                if ( vi<=xmax ) then
                    summ = summ + (vi/xmax)**2
                else
                    summ = one + summ*(xmax/vi)**2
                    xmax = vi
                end if
!
            case( null )
!
                if ( vi==zero  ) then
!
!                    JUST FALL THROUGH...
!
                else if ( vi<=cutlo ) then
                    rang = small
                    xmax = vi
                    summ = one
                else if ( vi>=hitest ) then
                    rang = large
                    xmax = vi
                    summ = one
                else
                    rang = normal
                    summ = vi*vi
                end if
!
        end select
!
    end do
!
!   COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
    select case( rang )
        case( normal, null ) ; norme_rv = sqrt( summ )
        case default         ; norme_rv = xmax*sqrt( summ )
    end select
!
!
! END OF FUNCTION norme_rv
! ________________________
!
    end function norme_rv
!
! =========================================================================================
!
    function norme_rm( mat )
!
! Purpose
! _______
!
!   This function computes the 2-norm (i.e. the Frobenius norm) of the matrix
!   MAT of size n, with due regard to avoiding overflow and underflow.
!
!
! Arguments
! _________
!
!   MAT    (INPUT) real(stnd), dimension(:,:)
!          On entry, the real matrix MAT.
!
!
! Further Details
! _______________
!
!   The routine is based on snrm2 from the blas (in linpack), but this
!   version is written in Fortran 90. It is machine independent. The
!   algorithm is described in reference (1).
!
!   The machine constants MACHTINY (the smallest magnitude), MACHBASE(base
!   of the machine), and MACHEPS (epsilon) are used to calculate the
!   constants cutlo and cuthi:
!
!       cutlo = sqrt( machsmlnum ) = sqrt( MACHTINY/(MACHEPS * MACHBASE) )
!       cuthi = one/cutlo
!
!   Three different cases must be considered when calculating the norm:
!
!       (1)  All components of MAT are below cutlo.
!
!            To avoid underflow, each component is divided by sqrt(min)/n
!            and then the regular Euclidean norm of this modified vector
!            is calculated.  This result is then multiplied by
!            sqrt(min)/n  in order to get the correct value for the norm.
!
!       (2)  One or more components are greater than cuthi.
!
!            To avoid overflow, the same method as in case (1) is used
!            with a scaling factor of   sqrt(max) * n .
!
!       (3)  All components are less than cuthi, with at least one
!            component greater than cutlo.
!
!            The regular formula for the Frobenius norm is used.
!
!   For more details, see:
!
!   (1) Blue, J.L., 1978:
!          A portable Fortran program to find the Euclidean norm of a vector.
!          ACM Trans. Math. Soft., Vol. 4, No 1, 1978, pp.15-23.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Reals_Constants,  only : zero, one
    use Num_Constants,    only : machsmlnum
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), dimension(:,:), intent(in) :: mat
    real(stnd)                             :: norme_rm
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)   :: hitest, summ, xmax, mij, cutlo
!
    integer(i4b) :: i, j, k, n
    integer(i2b) :: rang
!
!
! PARAMETERS
! __________
!
    integer(i2b), parameter ::   null   = 0,     &
                                 small  = 1,     &
                                 normal = 2,     &
                                 large  = 3
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    n = size( mat )
    k = size( mat, dim=1 )
!
    if ( n<=0_i4b ) then
        norme_rm = zero
        return
    end if
!
    cutlo  = sqrt( machsmlnum )
    hitest = one/( cutlo*real( n, stnd ) )
!
    summ  = zero
    rang  = null
!
!   EVALUATE THE NORM BY ACCUMULATING A SCALED SUM OF SQUARES AND
!   ADJUSTING THE SCALING AS NUMBERS OF INCREASINGLY LARGE MAGNITUDE
!   ARE FOUND.
!
    do  j = 1_i4b, size( mat, dim=2)
!
        do  i = 1_i4b, k
!
            mij = abs( mat(i,j) )
!
            select case( rang )
!
                case( normal )
!
                    if ( mij<hitest ) then
                        summ = summ + mij*mij
                    else
                        rang = large
                        xmax = mij
                        summ = one + (summ/mij)/mij
                    end if
!
                case( small )
!
                    if ( mij<=cutlo ) then
                        if ( mij<=xmax ) then
                            summ = summ + (mij/xmax)**2
                        else
                            summ = one + (xmax/mij)**2
                            xmax = mij
                        end if
                    else if ( mij>=hitest ) then
                        rang = large
                        xmax = mij
                        summ = one + (summ/mij)/mij
                    else
                        rang = normal
                        summ = (summ*xmax)*xmax + mij*mij
                    end if
!
                case( large )
!
                    if ( mij<=xmax ) then
                        summ = summ + (mij/xmax)**2
                    else
                        summ = one + summ*(xmax/mij)**2
                        xmax = mij
                    end if
!
                case( null )
!
                    if ( mij==zero ) then
!
!                       JUST FALL THROUGH...
!
                    else if ( mij<=cutlo ) then
                        rang = small
                        xmax = mij
                        summ = one
                    else if ( mij>=hitest ) then
                        rang = large
                        xmax = mij
                        summ = one
                    else
                        rang = normal
                        summ = mij*mij
                    end if
!
            end select
!
        end do
!
    end do
!
!   COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
    select case (rang)
        case( normal, null ) ; norme_rm = sqrt( summ )
        case default         ; norme_rm = xmax*sqrt( summ )
    end select
!
!
! END OF FUNCTION norme_rm
! ________________________
!
    end function norme_rm
!
! =========================================================================================
!
    function pythag( a, b )
!
! Purpose
! _______
!
!   Computes sqrt( a * a  +  b * b ) without destructive underflow or overflow.
!
!
! ___________________________________________________________________________________________
!
    use Reals_Constants,  only : one, zero
!
    real(stnd), intent(in) :: a, b
!
    real(stnd) :: pythag
!
    real(stnd) :: absa, absb
!
    absa = abs( a )
    absb = abs( b )
!
    if ( absa>absb ) then
!
        pythag = absa*sqrt( one+(absb/absa)**2 )
!
    else if ( absb==zero ) then
!
        pythag = zero
!
    else
!
        pythag = absb*sqrt( one+(absa/absb)**2 )
!
    end if
!
    end function pythag
!
! =========================================================================================
!
!
! *************************
! END OF MODULE Utilities *
! *************************
!    
end module Utilities
