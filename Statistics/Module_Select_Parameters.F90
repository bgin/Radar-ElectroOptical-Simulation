module Select_Parameters
! ####################################################################################
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
! ####################################################################################
!                                                                                    *
! ************************************************************************************
! THIS MODULE PROVIDES A CONVENIENT WAY OF SELECTING :                               * 
!                                                                                    *
! - THE PRECISION (KIND TYPES) REQUIRED FOR A COMPUTATION.                           *
!                                                                                    *
! - THE SIZE (KIND TYPES) OF INTEGER OR LOGICAL VARIABLES.                           *
!                                                                                    *
! - THE DEFAULT PRINTING UNIT.                                                       *
!                                                                                    *
! - THE DIFFERENT BLOCK SIZES FOR LINEAR ALGEBRA SUBROUTINES.                        *
!                                                                                    *
! - THE PARAMETERS FOR OpenMP COMPILATION.                                           *
!                                                                                    *
! - THE PARAMETERS FOR CROSSOVER FROM SERIAL TO PARALLEL ALGORITHMS.                 *
!                                                                                    *
! - THE PARAMETERS FOR THE STATPACK TESTING PROGRAMS.                                *
!                                                                                    *
! - THE LOCATION OF THE URANDOM DEVICE ON YOUR SYSTEM IF IT EXISTS.                  *
!                                                                                    * 
! IN ORDER TO CHANGE THE DEFAULT VALUES AND MAKE YOUR OWN CHOICE FOR THESE           *
! PARAMETERS, YOU MUST EDIT THE FILE Module_Select_Parameters.F90 AND FOLLOW THE     *
! INSTRUCTIONS IN THIS FILE.                                                         * 
!                                                                                    * 
!                                                                                    * 
! LATEST REVISION : 24/09/2020                                                       *
!                                                                                    *
! ####################################################################################
!                                                                                    *
! ************************************************************************************
!
!
! USED MODULES 
! ============
!
! ----------------------------------------------------------------------
! By simply ensuring that a leading '!' appears on all but exactly one !
! of the following use statements, and then recompiling all routines,  !
! the size of integer variables can be changed.                        !
! No harm will be done if short integers are made the same as i4b or   !
! i8b integers.                                                        !
! ----------------------------------------------------------------------
!
   use The_Kinds, only : i1b      , i2b      , i4b      , i8b
!   use The_Kinds, only : i1b=>i4b1, i2b=>i4b2, i4b      , i8b
!   use The_Kinds, only : i1b=>i4b1, i2b=>i4b2, i4b      , i8b=>i4b8
!   use The_Kinds, only : i1b=>i8b1, i2b=>i8b2, i4b=>i8b4, i8b
!   use The_Kinds, only : i1b      , i2b      , i4b=>i8b4, i8b
!
!
! -----------------------------------------------------------------------
! By simply ensuring that a leading '!' appears on all but exactly one  !
! of the following use statements, and then recompiling all routines,   !
! the precision of an entire real or complex computation can be altered.!
!                                                                       !
! A few computations are preferably done in higher precision 'extd'. So,!
! the kind type 'extd' should be such that the underlying hardware will !
! select a higher precision for kind 'extd' than for kind 'stnd', if    !
! this is feasible.  If a higher precision is not readily available,    !
! the same value may be used as for 'stnd'.                             !
! -----------------------------------------------------------------------
!
!   use The_Kinds, only : stnd=>sp,  extd=>dp
!
!   use The_Kinds, only : stnd=>dp,  extd=>qp
!
!   use The_Kinds, only : stnd=>sp,  extd=>sp2
!
   use The_Kinds, only : stnd=>dp,  extd=>dp2
!
!   use The_Kinds, only : stnd=>qp,  extd=>qp2
!
!   use The_Kinds, only : stnd=>low,     extd=>normal
!
!   use The_Kinds, only : stnd=>normal,  extd=>extended
!
!   use The_Kinds, only : stnd=>low,     extd=>low2
!
!   use The_Kinds, only : stnd=>normal,  extd=>normal2
!
! ----------------------------------------------------------------------
! By simply ensuring that a leading '!' appears on all but exactly one !
! of the following use statements, and then recompiling all routines,  !
! the size of logical variables can be changed.                        !
! ----------------------------------------------------------------------
!
   use The_Kinds, only : lgl=>logic
!   
!   use The_Kinds, only : lgl=>logic0
!   
!   use The_Kinds, only : lgl=>logic1
!   
!   use The_Kinds, only : lgl=>logic2
!   
!   use The_Kinds, only : lgl=>logic4
!   
!
! STRONG TYPING IMPOSED 
! =====================
!    
    implicit none
!   
!
! PUBLIC ENTITIES 
! ===============
!
    private
    public :: i1b, i2b, i4b, i8b, lgl, stnd, extd,        &
              defunit,                                    &
              blksz_util, blksz_lin, blksz_fft, blksz_qr, &
              blksz_eig, blksz2_eig, blksz2_svd,          &
              max_francis_steps_svd,                      &
              max_francis_steps_eig,                      &
              omp_limit, omp_limit2, omp_chunk,           &
              npar_arth, npar2_arth,                      &
              npar_geop, npar2_geop,                      &
              npar_cumsum, npar_cumprod,                  &
              npar_poly, npar_polyterm,                   &
              n1_def, n2_def, n3_def,                     &
              urandom_file
!
!
! MODULE PARAMETERS 
! =================
!
!
! DEFAULT PRINTING UNIT.
!
!
    integer(i4b), parameter :: defunit = 6
!
!
! BLOCK SIZES FOR LINEAR ALGEBRA SUBROUTINES.
!
!
! --------------------------------------------------------
! Block size for matrix multiplication and transposition !
! algorithms in module Utilities .                       !
! --------------------------------------------------------
!
    integer(i4b), parameter ::  blksz_util = 20_i4b
!
! -----------------------------------------------
! Block size for LU and Cholesky algorithms in  !
! module Lin_Procedures .                       !
! -----------------------------------------------
!
    integer(i4b), parameter ::  blksz_lin = 20_i4b
!
! ----------------------------------------
! Block size for parallel transposition  !
! algorithm in module FFT_Procedures .   !
! ----------------------------------------
!
    integer(i4b), parameter ::  blksz_fft = 20_i4b
!
! ----------------------------------------------
! Block size for QR and LQ algorithms in       !
! module QR_Procedures .                       !
! ----------------------------------------------
!
    integer(i4b), parameter ::  blksz_qr = 30_i4b
!
! ---------------------------------------------------
! Block sizes for reduction to tridiagonal form     !
! algorithms in module Eig_Procedures . blksz2_eig  !
! is only used by symtrid_cmp2 subroutine .         !
! ---------------------------------------------------
!
    integer(i4b), parameter ::  blksz_eig = 30_i4b
!
    integer(i4b), parameter ::  blksz2_eig = 30_i4b
!
! ---------------------------------------------------
! Block sizes for reduction to bidiagonal form      !
! algorithms in module SVD_Procedures . blksz2_svd  !
! is only used by bd_cmp2 and bd_cmp3 subroutines . !
! ---------------------------------------------------
!
    integer(i4b), parameter ::  blksz2_svd = 30_i4b
!
! -----------------------------------------------------
! Default block size for the wave front algorithm     !
! used in the bidiagonal QR SVD subroutines in module !
! SVD_Procedures .                                    !
! -----------------------------------------------------
!
    integer(i4b), parameter :: max_francis_steps_svd = 10_i4b
!
! ------------------------------------------------------
! Default block size for the wave front algorithm      !
! used in the tridiagonal QR EVD subroutines in module !
! Eig_Procedures .                                     !
! ------------------------------------------------------
!
    integer(i4b), parameter :: max_francis_steps_eig = 10_i4b
!
!
! PARAMETERS FOR OpenMP COMPILATION.
!
!
! -------------------------------------------------------
! Parallel regions working on a matrix are serialized   !
! if the size of the matrix is less then omp_limit .    !
! -------------------------------------------------------
!
    integer(i4b), parameter ::  omp_limit = 10000
!
! -------------------------------------------------------------
! Parallel regions working on a fixed matrix are serialized   !
! if the number of columns (or rows) of the matrix attributed !
! to each thread is less then omp_limit2 .                    !
! -------------------------------------------------------------
!
    integer(i4b), parameter :: omp_limit2 = 4
!
! --------------------------------------------------------------
! When do loops in a parallel region are executed in parallel  !
! with dynamic scheduling, the default chunk size is omp_chunk !
! --------------------------------------------------------------
!
    integer(i4b), parameter ::  omp_chunk = 16
!
!
! PARAMETERS FOR CROSSOVER FROM SERIAL TO PARALLEL ALGORITHMS,
! USED ONLY BY MODULE Utilities.
!
!
    integer(i4b), parameter ::  npar_arth     = 16,              &
                                npar2_arth    = 8,               &
                                npar_geop     = 4,               &
                                npar2_geop    = 2,               &
                                npar_cumsum   = 16,              &
                                npar_cumprod  = 8,               &
                                npar_poly     = 8,               &
                                npar_polyterm = 8
!
!
! PARAMETERS FOR THE statpack TESTING PROGRAMS.
!
!
! -------------------------------------------------------
! When testing statpack routines, these parameters are  !
! used to determine the size and number of matrices     !
! generated in the testing programs. Typically, a suite !
! of matrices of size n by n (or n by 2n), for n=n1_def !
! to n2_def by n3_def will be used in the testing       !
! programs. n1_def can be greater than n2_def, but in   !
! that case n3_def must be negative.                    !
! -------------------------------------------------------
!
    integer(i4b), parameter ::  n1_def    = 110,            &
                                n2_def    = 18,             &
                                n3_def    = -8
!
!
! PARAMETER FOR THE LOCATION OF THE URANDOM DEVICE ON YOUR SYSTEM
! IF IT EXISTS. FOR UNIX SYSTEM, NORMALLY IT IS /dev/urandom .
! THIS PARAMETER IS ONLY USED BY MODULE Random.
!
!
    character(len=*),  parameter :: urandom_file='/dev/urandom'
!
!
! **********************************
! END OF MODULE Select_Parameters  *
! **********************************
!
end module Select_Parameters
