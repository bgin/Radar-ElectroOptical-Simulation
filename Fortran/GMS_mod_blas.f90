

#include "Config.fpp"

module mod_blas


 !===================================================================================85
 !---------------------------- DESCRIPTION ------------------------------------------85
 !
 !
 !
 !          Module  name:
 !                         'mod_blas
 !          
 !          Purpose:
  !                      This module contains an explicitly vectorized complex
  !                      blas implementation, which is based on packed AoS
  !                      complex data type (AVX512c8f64_t).
 !
 !          History:
 !                        Date: 29-11-2019
 !                        Time: 10:41 GMT+2
 !                        
 !          Version:
 !
 !                      Major: 1
 !                      Minor: 0
 !                      Micro: 0
 !
 !          Author:  
 !                      Bernard Gingold
 !                      As per original authors request -- the name of the
 !                      modified subroutine/function was changed to include the 'gms' prefix.
 !                 
 !          References:
 !         
 !                            Based on LAPACK package
 !                            Copyright (c) 1992-2013 The University of Tennessee and The University
 !                            of Tennessee Research Foundation.  All rights
 !                            reserved.
 !                            Copyright (c) 2000-2013 The University of California Berkeley. All
 !                            rights reserved.
 !                            Copyright (c) 2006-2013 The University of Colorado Denver.  All rights
 !                            reserved.
 !
 !                            $COPYRIGHT$
 !
 !                            Additional copyrights may follow
 !         
 !          E-mail:
 !                  
 !                      beniekg@gmail.com
 !==================================================================================85
    ! Tab:5 col - Type and etc.. definitions
    ! Tab:10,11 col - Type , function and subroutine code blocks.

    use module_kinds, only : int4,dp
    use mod_avx512c8f64
    implicit none

     !=====================================================59
     !  File and module information:
     !  version,creation and build date, author,description
     !=====================================================59

    ! Major version
    integer(kind=int4), parameter, public :: MOD_BLAS_MAJOR = 1
    ! MInor version
    integer(kind=int4), parameter, public :: MOD_BLAS_MINOR = 0
    ! Micro version
    integer(kind=int4), parameter, public :: MOD_BLAS_MICRO = 0
    ! Module full version
    integer(kind=int4), parameter, public :: MOD_BLAS_FULLVER = &
         1000*MOD_BLAS_MAJOR+100*MOD_BLAS_MINOR+10*MOD_BLAS_MICRO
    !Module creation date
    character(*),       parameter, public :: MOD_BLAS_CREATION_DATE = "29-11-2019 10:55 +00200 (FRI 29 NOV 2019 GMT+2)"
    ! Module build date
    character(*),       parameter, public :: MOD_BLAS_BUILD_DATE    = __DATE__ " " __TIME__
    ! Module author info
    character(*)        parameter, public :: MOD_BLAS_AUTHOR = "LAPACK original authors[all rights reserved] -- This version was  modified by Bernard Gingold, contact: beniekg@gmail.com"
    ! Module short description
    character(*)        parameter, public :: MOD_BLAS_SYNOPSIS = "Explicitly vectorized complex  blas implementation, which is based on packed AoS complex data type (AVX512c8f64_t) "
                                                                  
    ! public

  contains

!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2017
!*
!*> \ingroup complex16_blas_level1
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>     jack dongarra, 3/11/78.
!*>     modified 12/3/93, array(1) declarations changed to array(*)
!       Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)
!*> \endverbatim
!*>

    subroutine gms_zaxpy(n,za,zx,incx,zy,incy)
#if defined __INTEL_COMPILER      
      !DIR$ ATTRIBUTES CODE_ALIGNED : 32 :: gms_zaxpy
      !DIR$ ATTRIBUTES VECTOR :: gms_zaxpy
#endif
      use mod_vecconsts, only : v8_n0
      integer(kind=int4),                intent(in),value    :: n
      type(AVX512c8f64_t),               intent(in)          :: za
      type(AVX512c8f64_t), dimension(*), intent(in)          :: zx
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED zx:64
#endif
      integer(kind=int4),                intent(in),value    :: incx
      type(AVX512c8f64_t), dimension(*), intent(inout)       :: zy
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED zy:64
#endif
      integer(kind=int4),                intent(in),value    :: incy
      ! Locals
      integer(kind=int4), automatic :: i,ix,iy
      ! EXec code .....
      if(n<=0) return
      if(all(cabs_zmm8c8(za) == v8_n0.v)) return
      if(incx==1 .and. incy==1) then
         ! *        code for both increments equal to 1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#endif
         do i = 1,n
            zy(i) = zy(i)+za*zx(i)
         end do

         !*        code for unequal increments or equal increments
         !*          not equal to 1
      else
         ix=1
         iy=1
         if(incx<0) ix=(-n+1)*incx+1
         if(incy<0) iy=(-n+1)*incy+1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#endif
         do i = 1,n
            zy(iy) = zy(iy)+za*zx(ix)
            ix = ix+incx
            iy = iy+incy
         end do
      end if
    end subroutine gms_zaxpy

!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2017
!*
!*> \ingroup complex16_blas_level1
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>     jack dongarra, linpack, 4/11/78.
!*>     modified 12/3/93, array(1) declarations changed to array(*)
!       Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)    
!*> \endverbatim
!*>

    subroutine gms_zcopy(n,zx,incx,zy,incy)
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zcopy
      !DIR$ ATTRIBUTES VECTOR :: gms_zcopy
#endif
      integer(kind=int4),                intent(in),value  :: n
      type(AVX512c8f64_t), dimension(*), intent(in)        :: zx
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED zx:64
#endif
      integer(kind=int4),                intent(in),value  :: incx
      type(AVX512c8f64_t), dimension(*), intent(out)       :: zy
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED zy:64
#endif
      integer(kind=dint4),               intent(in),value  :: incy
      ! LOcals
      integer(kind=int4), automatic :: i,ix,iy
      ! EXec code ...
      if(n<=0) return
      if(incx==1 .and. incy==1) then
         
         !  code for both increments equal to 1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(4)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(3)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(4)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(3)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
         do i = 1,n
            zy(i) = zx(i)
         end do

         !  code for unequal increments or equal increments
         !*          not equal to 1
       else
         ix=1
         iy=1
         if(incx<0) ix=(-n+1)*incx+1
         if(incy<0) iy=(-n+1)*incy+1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(4)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(3)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(4)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(3)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif        
         do i = 1,n
            zy(iy) = zx(ix)
            ix = ix+incx
            iy = iy+incy
         end do
      end if
    end subroutine gms_zcopy

!     Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2017
!*
!1*> \ingroup complex16_blas_level1
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>     jack dongarra, 3/11/78.
!*>     modified 12/3/93, array(1) declarations changed to array(*)
!       Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)     
!*> \endverbatim
    
    function gms_zdotc(n,zx,incx,zy,incy) result(dotc)
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zdotc
      !DIR$ ATTRIBUTES VECTOR :: gms_zdotc
#endif
      integer(kind=int4),                intent(in),value :: n
      type(AVX512c8f64_t), dimension(*), intent(in)       :: zx
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED zx:64
#endif
      integer(kind=int4),                intent(in),value :: incx
      type(AVX512c8f64_t), dimension(*), intent(in)       :: zy
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED zy:64
#endif
      integer(kind=int4),                intent(in),value :: incy
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES ALIGN : 64 :: dotc
#endif
      type(AVX512c8f64_t) :: dotc
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES ALIGN : 64 :: ztemp
#endif
      type(AVX512c8f64_t), automatic :: ztemp
      integer(kind=int4), automatic :: i,ix,iy
      ! EXec code ....
      !
      if(n<=0) return
      ztemp = default_init()
      if(incx==1 .and. incy==1) then
         !   code for both increments equal to 1
#if defined __INTEL_COMPILER         
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#endif
         do i = 1,n
            ztemp = ztemp+conjugate(zx(i))*zy(i)
         end do
      else
         !   code for unequal increments or equal increments
         !*          not equal to 1
         ix=1
         iy=1
         if(incx<0) ix=(-n+1)*incx+1
         if(incy<0) iy=(-n+1)*incy+1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#endif
         do i = 1,n
            ztemp = ztemp+conjugate(zx(ix))*zy(iy)
            ix = ix+incx
            iy = iy+incy
         end do
         zdotc = default_init()
      end if  
      zdotc = ztemp
    end function gms_zdotc

!     Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2017
!*
!1*> \ingroup complex16_blas_level1
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>     jack dongarra, 3/11/78.
!*>     modified 12/3/93, array(1) declarations changed to array(*)
!       Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)     
!*> \endverbatim

    function gms_zdotu(n,zx,incx,zy,incy) result(zdotu)
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zdotu
      !DIR$ ATTRIBUTES VECTOR :: gms_zdotu
#endif
      integer(kind=int4),                intent(in),value :: n
      type(AVX512c8f64_t), dimension(*), intent(in)       :: zx
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED zx:64
#endif
      integer(kind=int4),                intent(in),value :: incx
      type(AVX512c8f64_t), dimension(*), intent(in)       :: zy
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED zy:64
#endif
      integer(kind=int4),                intent(in),value :: incy
      ! LOcals
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES ALIGN : 64 :: zdotu
#endif
      type(AVX512c8f64_t) :: zdotu
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES ALIGN : 64 :: ztemp
#endif
      type(AVX512c8f64_t), automatic :: ztemp
      integer(kind=int4),  automatic :: i,ix,iy
      ! Exec code ....
      if(n<=0) return
      ztemp = default_init()
      if(incx==1 .and. incy==1) then

         !  code for both increments equal to 1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#endif        
         do i = 1,n
            ztemp = ztemp+zx(i)*zy(i)
         end do

         !  code for unequal increments or equal increments
         !*          not equal to 1
      else
         ix=1
         iy=1
         if(incx<0) ix=(-n+1)*incx+1
         if(incy<0) iy=(-n+1)*incy+1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#endif
         do i = 1,n
            ztemp = ztemp+zx(ix)*zy(iy)
            ix = ix+incx
            iy = iy+incy
         end do
      end if
      zdotu = default_init()
      zdotu = ztemp
    end function gms_zdotu

!    *
!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*  Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)     
!*> \date December 2016

    subroutine gms_zdrot(n,cx,incx,cy,incy,c,s)
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zdrot
      !DIR$ ATTRIBUTES VECTOR :: gms_zdrot
      use mod_vectypes, only : ZMM8r8_t
      integer(kind=int4),                intent(in),value    :: n
      type(AVX512c8f64_t), dimension(*), intent(inout)       :: cx
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED cx:64
#endif
      integer(kind=int4),                intent(in),value    :: incx
      type(AVX512c8f64_t), dimension(*), intent(inout)       :: cy
#if defined __INTEL_COMPILER      
      !DIR$ ASSUME_ALIGNED cy:64
#endif
      integer(kind=int4),                intent(in),value    :: incy
      type(ZMM8r8_t),                    intent(in)          :: c ! scalar extended to vector
     
      type(ZMM8r8_t),                    intent(in)          :: s ! scalar extended to vector
#if defined __INTEL_COMPILER     
      !DIR$ ATTRIBUTES ALIGN : 64 :: ztemp
#endif
      type(AVX512c8f64_t), automatic :: ztemp
      integer(kind=int4),  automatic :: i,ix,iy
      ! EXec code ...
      if(n<=0) return
      ctemp = default_init()
      if(incx==1 .and. incy==1) then
         !  code for both increments equal to 1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#endif
         do i = 1,n
            ctemp = c*cx(i)+s*cy(i)
            cy(i) = c*cy(i)-s*cx(i)
            cx(i) = ctemp
         end do

         !  code for unequal increments or equal increments not equal
         !*          to 1
      else
         ix=1
         iy=1
         if(incx<0) ix=(-n+1)*incx+1
         if(incy<0) iy=(-n+1)*incy+1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#endif
         do i = 1,n
            ctemp  = c*cx(ix)+s*cy(iy)
            cy(iy) = c*cy(iy)-s*cx(ix)
            cx(ix) = ctemp
            ix = ix+incx
            iy = iy+incy
         end do
       end if
    end subroutine gms_zdrot

!     Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2017
!1*
!*> \ingroup complex16_blas_level1
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>     jack dongarra, 3/11/78.
!*>     modified 3/93 to return if incx .le. 0.
!*>     modified 12/3/93, array(1) declarations changed to array(*)
!       Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)       
!*> \endverbatim

    subroutine gms_zdscal(n,da,zx,incx)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zdscal
      !DIR$ ATTRIBUTES VECTOR :: gms_zdscal
#endif
       integer(kind=int4),                intent(in),value    :: n
       type(AVX512c8f64_t),               intent(in)          :: da
       type(AVX512c8f64_t), dimension(*), intent(inout)       :: zx
#if defined __INTEL_COMPILER
       !DIR$ ASSUME_ALIGNED zx:64
#endif
       integer(kind=int4),                intent(in),value    :: incx
       ! LOcals
       integer(kind=int4), automatic :: i,nincx
       ! Exec code ....
       if(n<=0 .or. incx<0) return
       if(incx==1) then
          !  code for increment equal to 1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                  
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
          do i = 1,n
             zx(i) = da*zx(i)
          end do
       else
          !   *        code for increment not equal to 1
          nincx=n*incx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                  
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
          do i = 1,nincx,incx
             zx(i) = da*zx(i)
          end do
       end if
    end subroutine gms_zdscal

!    *> \date December 2016
!*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
!*>     Richard Hanson, Sandia National Labs.
!  Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
!*> \endverbatim

    subroutine gms_zgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zgbmv
#endif
        character(len=1),                      intent(in),value :: trans
        integer(kind=int4),                    intent(in),value :: m
        integer(kind=int4),                    intent(in),value :: n
        integer(kind=int4),                    intent(in),value :: kl
        integer(kind=int4),                    intent(in),value :: ku
        type(AVX512c8f64_t),                   intent(in)       :: alpha
        type(AVX512c8f64_t), dimension(lda,*), intent(in)       :: a
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED a:64
#endif
        integer(kind=int4),                    intent(in),value :: lda
        type(AVX512c8f64_t), dimension(*),     intent(in)       :: x
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED x:64
#endif
        integer(kind=int4),                    intent(in),value :: incx
        type(AVX512c8f64_t),                   intent(in)       :: beta
        type(AVX512c8f64_t), dimension(*),     intent(inout)    :: y
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED y:64
#endif
        integer(kind=int4),                    intent(in),value :: incy
        ! LOcals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
        !
        integer(kind=dp),    automatic :: i,info,ix,iy,j,jk,k,kup1,kx,ky,lenx,leny
        logical(kind=int4),  automatic :: noconj
        logical(kind=int1),  automatic :: beq0,aeq0,bneq0,aneq0,beq1,aeq1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ONE
#endif
        type(AVX512c8f64_t), parameter :: ONE  = AVX512c8f64_t([1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                                                                1.0_dp,1.0_dp,1.0_dp,1.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                 0.0_dp,0.0_dp,0.0_dp,0.0_dp])
 
#if defined __INTEL_COMPILER                                                               
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
            
             
        ! EXec code ....
        info = 0
        if(.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. &
           .not.lsame(trans,'C')) then
           info = 1
        else if(m<0) then
           info = 2
        else if(n<0) then
           info = 3
        else if(kl<0) then
           info = 4
        else if(ku<0) then
           info = 5
        else if(lda<(kl+ku+1)) then
           info = 8
        else if(incx==0) then
           info = 10
        else if(incy==0) then
           info = 13
        end if
        if(info/=0) then
           call xerbla('GMS_ZGBMV',info)
           return
        end if
        !  Quick return if possible.
        aeqz=all(alpha==ZERO)
        beq1=all(beta==ONE)
        if((m==0) .or. (n==0) .or. &
             ((aeqz) .and. (beq1))) return
        noconj = lsame(trans,'T')
        !  Set  LENX  and  LENY, the lengths of the vectors x and y, and set
        !*     up the start points in  X  and  Y.
        if(lsame(trans,'N')) then
           lenx = n
           leny = m
        else
           lenx = m
           leny = n
        end if
        if(incx>0) then
           kx = 1
        else
           kx = 1-(lenx-1)*incx
        end if
        if(incy>0) then
           ky = 1
        else
           ky = 1-(leny-1)*incy
        end if
        ! *     Start the operations. In this version the elements of A are
        ! *     accessed sequentially with one pass through the band part of A.
        ! *
        ! *     First form  y := beta*y.
        VCZERO = default_init()
        bneq1=all(beta/=ONE)
        beq0=all(beta==ZERO)
        if(bneq1) then
           if(incy==1) then
              if(beq0) then
#if defined __INTEL_COMPILER
                 !DIR$ VECTOR ALIGNED
                 !DIR$ VECTOR ALWAYS
#endif
                 do i=1,leny
                    y(i) = ZERO
                 end do
              else
#if defined __INTEL_COMPILER
                 !DIR$ VECTOR ALIGNED
                 !DIR$ VECTOR ALWAYS
#endif
                 do i=1,leny
                    y(i) = beta*y(i)
                 end do
              end if
           else
              iy=ky
              if(beq0) then
#if defined __INTEL_COMPILER
                 !DIR$ VECTOR ALIGNED
                 !DIR$ VECTOR ALWAYS
#endif
                 do i=1,leny
                    y(iy) = ZERO
                    iy = iy+incy
                 end do
              else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                  
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                 do i=1,leny
                    y(iy) = beta*y(iy)
                    iy = iy+incy
                 end do
              end if
           end if
        end if
        if(aeq0) return
        kup1 = ku+1
        if(lsame(trans,'N')) then
           !   Form  y := alpha*A*x + y.
           jx=kx
           if(incy==1) then
              do j=1,n
                 temp=alpha*x(jx)
                 k=kup1-j
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#endif
                 do i=max(1,j-ku),min(m,j+kl)
                    y(i) = y(i)+temp*a(k+1,j)
                 end do
                 jx = jx+incx
              end do
           else
              do j=1,n
                 temp=alpha*x(jx)
                 iy=ky
                 k=kup1-j
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#endif
                 do i=max(1,j-ku),min(m,j+kl)
                    y(iy) = y(iy)+temp*a(k+1,j)
                    iy=iy+incy
                 end do
                 jx=jx+incx
                 if(j>ku) ky=ky+incy
              end do
           end if
        else
           !  Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
           jy=ky
           if(incx==1) then
              do j=1,n
                 temp=ZERO
                 k=kup1-j
                 if(noconj) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#endif
                    do i=max(1,j-ku),min(m,j+kl)
                       temp = temp+a(k+i,j)*x(i)
                    end do
                 else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#endif
                    do i=,max(1,j-ku),min(m,j+kl)
                       temp = temp+conjugate(a(k+i,j))*x(i)
                    end do
                 end if
                 y(jy) = y(jy)+alpha*temp
                 jy=jy+incy
              end do
           else
              do j=1,n
                 temp=ZERO
                 ix=kx
                 k=kup1-j
                 if(noconj) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#endif
                    do i=max(1,j-ku),min(m,j+kl)
                       temp = temp+a(k+i,j)*x(ix)
                       ix = ix+incx
                    end do
                 else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#endif
                    do i=max(1,j-ku),min(m,j+kl)
                       temp = temp+conjugate(a(k+i,j)*x(ix)
                       ix = ix+incx
                    end do
                 end if
                 y(jy) = y(jy)+alpha*temp
                 jy = jy+incy
                 if(j>ku) kx = kx+incx
              end do
           end if
        end if
        
    end subroutine gms_zgbmv

!    *  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level3
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 3 Blas routine.
!*>
!1*>  -- Written on 8-February-1989.
!*>     Jack Dongarra, Argonne National Laboratory.
!*>     Iain Duff, AERE Harwell.
!*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
    !*>     Sven Hammarling, Numerical Algorithms Group Ltd.
!     Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
!1*> \endverbatim
    !*>

    subroutine gms_zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zgemm
#endif
       character(len=1),                      intent(in),value    :: transa
       character(len=1),                      intent(in),value    :: transb
       integer(kind=int4),                    intent(in),value    :: m
       integer(kind=int4),                    intent(in),value    :: n
       integer(kind=int4),                    intent(in),value    :: k
       type(AVX512c8f64_t),                   intent(in)          :: alpha
       type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
#if defined __INTEL_COMPILER
       !DIR$ ASSUME_ALIGNED a:64
#endif
       integer(kind=int4),                    intent(in),value    :: lda
       type(AVX512c8f64_t), dimension(ldb,*), intent(in)          :: b
#if defined __INTEL_COMPILER
       !DIR$ ASSUME_ALIGNED b:64
#endif
       integer(kind=int4),                    intent(in),value    :: ldb
       type(AVX512c8f64_t),                   intent(in)          :: beta
       type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
#if defined __INTEL_COMPILER
       !DIR$ ASSUME_ALIGNED c:64
#endif
       integer(kind=int4),                    intent(in),value    :: ldc
       ! Locals
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
       type(AVX512c8f64_t), automatic :: temp
       integer(kind=int4),  automatic :: i,info,j,l,ncola,nrowa,nrowb
       logical(kind=int4),  automatic :: conja,conjb,nota,notb
       logical(kind=int1),  automatic :: aeq0,beq0,beq1,bneq1
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES ALIGN : 64 :: ONE
#endif
       type(AVX512c8f64_t), parameter :: ONE  = AVX512c8f64_t([1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                                                                1.0_dp,1.0_dp,1.0_dp,1.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
#if defined __INTEL_COMPILER       
       !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
       type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
       ! EXec code ....
!      Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!*     B  respectively are to be  transposed but  not conjugated  and set
!*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
!1*    and the number of rows of  B  respectively.
        nota  = lsame(transa,'N')
        notb  = lsame(transb,'N')
        conja = lsame(transa,'C')
        conjb = lsame(transb,'C')
        if(nota) then
          nrowa = m
          ncola = k
        else
          nrowa = k
          ncola = m
        end if
        if(notb) then
          nrowb = k
        else
          nrowb = n
        end if
        !    Test the input parameters.
        info = 0
        if((.not.nota) .and. (.not.conja) .and. &
           (.not.lsame(transa,'T'))) then
           info = 1
        else if((.not.notb) .and. (.not.conjb) .and. &
           (.not.lsame(transb,'T'))) then
           info = 2
        else if(m<0) then
           info = 3
        else if(n<0) then
           info = 4
        else if(k<0) then
           info = 5
        else if(lda < max(1,nrowa)) then
           info = 8
        else if(ldb < max(1,nrowb)) then
           info = 10
        else if(ldc < max(1,m)) then
           info = 13
        end if
        if(info/=0) then
           call xerbla('GMS_ZGEMM',info)
        end if
        aeq0 = all(alpha==ZERO)
        beq1 = all(beta==ONE)
        ! Early exit
        if((m==0) .or. (n==0) .or. &
             (((aeq0) .or. (k==0)) .and. (beq1))) return
        !   And when  alpha.eq.zero.
        beq0 = all(beta==ZERO)
        if(aeq0) then
           if(beq0) then
              do j=1,n
#if defined __INTEL_COMPILER
                 !DIR$ VECTOR ALIGNED
                 !DIR$ VECTOR ALWAYS
#endif
                 do i=1,m
                    c(i,j) = ZERO
                 end do
              end do
           else
              do j=1,n
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                 do i=1,m
                    c(i,j) = beta*c(i,j)
                 end do
              end do
           end if
           return
        end if
        !  Start the operations.
        bneq1 = all(beta/=ONE)
        if(notb) then
           if(nota) then
              !  Form  C := alpha*A*B + beta*C.
              do j=1,n
                 if(beq0) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#endif
                    do i=1,m
                       c(i,j) = ZERO
                    end do
                 else if(bneq0) then
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=1,m
                       c(i,j) = beta*c(i,j)
                    end do
                 end if
                 do l=1,k
                    temp = alpha*b(l,j)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=1,m
                       c(i,j) = c(i,j)+temp*a(i,l)
                    end do
                 end do
              end do
           else if(conja) then
              !   Form  C := alpha*A**H*B + beta*C.
              do j=1,n
                 do i=1,m
                    temp = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do l=1,k
                       temp = temp+conjugate(a(l,i))*b(l,j)
                    end do
                    if(beq0) then
                       c(i,j) = alpha*temp
                    else
                       c(i,j) = alpha*temp+beta*c(i,j)
                    end if
                 end do
              end do
           else
              !  Form  C := alpha*A**T*B + beta*C
              do j=1,n
                 do i=1,m
                    temp = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do l=1,k
                       temp = temp+a(l,i)*b(l,j)
                    end do
                    if(beq0) then
                       c(i,j) = alpha*temp
                    else
                       c(i,j) = alpha*temp+beta*c(i,j)
                    end if
                 end do
              end do
           end if
        else if(nota) then
               if(conjb) then
                 !  Form  C := alpha*A*B**H + beta*C.
                  do j=1,n
                     if(beq0) then
#if defined __INTEL_COMPILER
                        !DIR$ VECTOR ALIGNED
                        !DIR$ VECTOR ALWAYS
#endif
                        do i=1,m
                           c(i,j) = ZERO
                        end do
                     else if(bneq1) then
#if defined __INTEL_COMPILER
                        !DIR$ VECTOR ALIGNED
                        !DIR$ VECTOR ALWAYS
#endif
                        do i=1,m
                           c(i,j) = beta*c(i,j)
                        end do
                     end if
                     do l=1,k
                        temp = alpha*conjugate(b(j,l))
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                        do i=1,m
                           c(i,j) = c(i,j)+temp*a(i,l)
                        end do
                     end do
                  end do
               else
                  !  Form  C := alpha*A*B**T + beta*C
                  do j=1,n
                     if(beq0) then
#if defined __INTEL_COMPILER
                        !DIR$ VECTOR ALIGNED
                        !DIR$ VECTOR ALWAYS
#endif
                        do i=1,m
                           c(i,j) = ZERO
                        end do
                     else if(bneq1) then
                        do i=1,m
                           c(i,j) = beta*c(i,j)
                        end do
                     end if
                     do l=1,k
                        temp = alpha*b(j,l)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                        do i=1,m
                           c(i,j) = c(i,j)+temp*a(i,l)
                        end do
                     end do
                  end do
               end if
            else if(conja) then
                  if(conjb) then
                     !   Form  C := alpha*A**H*B**H + beta*C.
                     do j=1,n
                        do i=1,m
                           temp = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                           do l=1,k
                              temp = temp+conjugate(a(l,i))*conjugate(b(j,l))
                           end do
                           if(beq0) then
                              c(i,j) = alpha*temp
                           else
                              c(i,j) = alpha*temp+beta*c(i,j)
                           end if
                        end do
                     end do
                  else
                     !   Form  C := alpha*A**H*B**T + beta*C
                     do j=1,n
                        do i=1,m
                           temp = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                           do l=1,k
                              temp = temp+conjugate(a(l,i))*b(j,l)
                           end do
                           if(beq0) then
                              c(i,j) = alpha*temp
                           else
                              c(i,j) = alpha*temp+beta*c(i,j)
                           end if
                        end do
                     end do
                  end if
               else
                  if(conjb) then
                     !  Form  C := alpha*A**T*B**H + beta*C
                     do j=1,n
                        do i=1,m
                           temp = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                           do l=1,k
                              temp = temp+a(l,i)*conjugate(b(j,l))
                           end do
                           if(beq0) then
                              c(i,j) = alpha*temp
                           else
                              c(i,j) = alpha*temp+beta*c(i,j)
                           end if
                        end do
                     end do
                  else
                     !  Form  C := alpha*A**T*B**T + beta*C
                     do j=1,n
                        do i=1,m
                           temp = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                           do l=1,k
                              temp = temp+a(l,i)*b(j,l)
                           end do
                           if(beq0) then
                              c(i,j) = alpha*temp
                           else
                              c(i,j) = alpha*temp+beta*c(i,j)
                           end if
                        end do
                     end do
                  end if
               end if
               ! End of GMS_ZGEMM
      end subroutine gms_zgemm

!       Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
      !*>     Richard Hanson, Sandia National Labs.
      !  Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
      !*> \endverbatim

      subroutine gms_zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zgmev
#endif
          character(len=1),                      intent(in),value :: trans
          integer(kind=int4),                    intent(in),value :: m
          integer(kind=int4),                    intent(in),value :: n
          type(AVX512c8f64_t),                   intent(in)       :: alpha
          type(AVX512c8f64_t), dimension(lda,*), intent(in)       :: a
#if defined __INTEL_COMPILER
          !DIR$ ASSUME_ALIGNED a:64
#endif
          integer(kind=int4),                    intent(in),value :: lda
          type(AVX512c8f64_t), dimension(*),     intent(in)       :: x
#if defined __INTEL_COMPILER
          !DIR$ ASSUME_ALIGNED x:64
#endif
          integer(kind=int4),                    intent(in),value    :: incx
          type(AVX512c8f64_t),                   intent(in)          :: beta
          type(AVX512c8f64_t), dimension(*),     intent(inout)       :: y
#if defined __INTEL_COMPILER
          !DIR$ ASSUME_ALIGNED y:64
#endif
          integer(kind=int4),                    intent(in),value    :: incy
          ! LOcals
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
          type(AVX512c8f64_t), automatic :: temp
          integer(kind=int4),  automatic :: i,info,ix,iy,j,jx.jy,kx,ky,lenx,leny
          logical(kind=int4),  automatic :: noconj
          logical(kind=int1),  automatic :: aeq0,beq1,bneq1,beq0
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: ONE
#endif
          type(AVX512c8f64_t), parameter :: ONE  = AVX512c8f64_t([1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                                                                1.0_dp,1.0_dp,1.0_dp,1.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
          type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
          ! EXec code .....
          info = 0
          if(.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. &
             .not.lsame(trans,'C')) then
             info = 1
          else if(m<0) then
             info = 2
          else if(n<0) then
             info = 3
          else if(lda < max(1,m)) then
             info = 6
          else if(incx==0) then
             info = 8
          else if(incy==0) then
             info = 11
          end if
          if(info/=0) then
             call xerbla('GMS_ZGEMV',info)
             return
          end if
          aeq0  = .false.
          beq1  = .false.
          bneq1 = .false.
          beq0  = .false.
          ! Quick return if possible.
          aeq0 = all(alpha==ZERO)
          beq1 = all(beta==ONE)
          if((m==0)  .or. (n==0) .or. &
               ((aeq0) .and. (beq1))) return
          noconj = lsame(trans,'T')
          !  Set  LENX  and  LENY, the lengths of the vectors x and y, and set
          !  *     up the start points in  X  and  Y.
          if (lsame(trans,'N')) then
             lenx = n
             leny = m
          else
             lenx = m
             leny = n
          end if
          if (incx>0) then
             kx = 1
          else
             kx = 1-(lenx-1)*incx
          end if
          if (incx>0) then
             ky = 1
          else
             ky = 1-(leny-1)*incy
          end if
          ! *     Start the operations. In this version the elements of A are
          ! *     accessed sequentially with one pass through A.
          ! *
          ! *     First form  y := beta*y.
          bneq1 = all(beta/=ONE)
          beq0  = all(beta==ZERO)
          if(bneq1) then
             if(incy==1) then
                if(beq0)  then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                   !DIR$ VECTOR ALWAYS
#endif
                    do i=1,leny
                       y(i) = ZERO
                    end do
                 else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#endif
                    do i=1,leny
                       y(i) = beta*leny(i)
                    end do
                 end if
              else
                 iy = ky
                 if(beq0) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#endif
                    do i=1,leny
                       y(iy) = ZERO
                       iy = iy+incy
                    end do
                 else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#endif
                    do i=1,leny
                       y(iy) = beta*y(iy)
                       iy = iy+incy
                    end do
                 end if
              end if
           end if
           if(aeq0) return
           if(lsame(trans,'N')) then
              ! Form  y := alpha*A*x + y.
              jx = kx
              if(incy==1) then
                 do j=1,n
                    temp = alpha*x(jx)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=1,m
                       y(i) = y(i)+temp*a(i,j)
                    end do
                    jx = jx+incx
                 end do
              else
                 do j=1,n
                    temp = alpha*x(jx)
                    iy = ky
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=1,m
                       y(iy) = y(iy)+temp*a(i,j)
                       iy = iy+incy
                    end do
                    jx = jx+incx
                 end do
              end if
           else
              !  Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
              jy = ky
              if(incx==1) then
                 do j=1,n
                    temp = zero
                    if(noconj) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=1,m
                          temp = temp+a(i,j)*x(i)
                       end do
                    else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=1,m
                          temp = temp+conjugate(a(i,j))*x(i)
                       end do
                    end if
                    y(jy) = y(jy)+alpha*temp
                    jy = jy+incy
                 end do
              else
                 do j=1,n
                    temp = zero
                    ix = ky
                    if (noconj) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=1,m
                          temp = temp+a(i,j)*x(ix)
                          ix = ix+incx
                       end do
                    else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=1,m
                          temp = temp+conjugate(a(i,j))*x(ix)
                          ix = ix+incx
                       end do
                    end if
                       y(jy) = y(jy)+alpha*temp
                       jy = jy+incy
                  end do
               end if
            end if
            ! End of ZGEMV
     end subroutine gms_zgemv

!      Authors:
!1*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!1*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>
!*>  -- Written on 22-October-1986.
!1*>     Jack Dongarra, Argonne National Lab.
!1*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
     !*>     Richard Hanson, Sandia National Labs.
     !   Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
!*> \endverbatim
     !*>

     subroutine gms_zgerc(m,n,alpha,x,incx,y,incy,a,lda)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zgerc
#endif
        integer(kind=int4),                    intent(in),value    :: m
        integer(kind=int4),                    intent(in),value    :: n
        type(AVX512c8f64_t),                   intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED x:64
#endif
        integer(kind=int4),                    intent(in),value    :: incx
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: y
#if defined __INTEL_COMPILER
        !DIR$ ASUME_ALIGNED y:64
#endif
        integer(kind=int4),                    intent(in),value    :: incy
        type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED a:64
#endif
        integer(kind=int4),                    intent(in),value    :: lda
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
        integer(kind=int4),  automatic :: i,info,ix,j,jy,kx
        logical(kind=int1),  automatic :: aeq0
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        ! EXec code .....
        info = 0
        if(m<0) then
           info = 1
        else if(n<0) then
           info = 2
        else if(incx==0) then
           info = 5
        else if(incy==0) then
           info = 7
        else if(lds < max(1,m)) then
           info = 9
        end if
        if(info/=0) then
           call xerbla('ZGERC',info)
           return
        end if
        !  Quick return if possible.
        aeq0 = .false.
        aeq0 = all(alpha==ZERO)
        if((m==0) .or. (n==0) .or. (aeq0)) return
        ! Start the operations. In this version the elements of A are
        ! *     accessed sequentially with one pass through A.
        if(incy>0) then
           jy = 1
        else
           jy = 1-(n-1)*incy
        end if
        if(incx==1) then
           do j=1,n
              if(all(y(jy)/=ZERO)) then
                 temp = alpha*conjugate(y(jy))
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                 do i=1,m
                    a(i,j) = a(i,j)+x(i)*temp
                 end do
              end if
              jy = jy+incy
           end do
        else
           if(incx>0) then
              kx = 1
           else
              kx = 1-(m-1)*incx
           end if
           do j=1,n
              if(all(y(jy)/=ZERO)) then
                 temp = alpha*conjugate(y(jy))
                 ix = kx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                 do i=1,m
                    a(i,j) = a(i,j)+x(ix)*temp
                    ix = ix+incx
                 end do
              end if
              jy = jy+incy
           end do
        end if
        ! End of zgerc
     end subroutine gms_zgerc

!      Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!1*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
     !*>     Richard Hanson, Sandia National Labs.
     !   Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
     !*> \endverbatim

     subroutine gms_zgeru(m,n,alpha,x,incx,y,incy,a,lda)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zgeru
#endif
        integer(kind=int4),                    intent(in),value    :: m
        integer(kind=int4),                    intent(in),value    :: n
        type(AVX512c8f64_t),                   intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED x:64
#endif
        integer(kind=int4),                    intent(in),value    :: incx
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: y
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED y:64
#endif
        integer(kind=int4),                    intent(in),value    :: incy
        type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
        integer(kind=int4),                    intent(in),value    :: lda
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
        integer(kind=int4),  automatic :: i,info,ix,j,jy,kx
        logical(kind=int1),  automatic :: aeq0
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        !  Test the input parameters.
        info = 0
        aeq0 = .false.
        if(m<0) then
           info = 1
        else if(n<0) then
           info = 2
        else if(incx==0) then
           info = 5
        else if(incy==0) then
           info = 7
        else if(lda < max(1,m)) then
           info = 9
        end if
        if(info/=0) then
           call xerbla('GMS_ZGERU',info)
           return
        end if
        !  Quick return if possible.
        aeq0 = all(alpha==ZERO)
        if((m==0) .or. (n==0) .or. (aeq0)) return
        !   Start the operations. In this version the elements of A are
        ! *     accessed sequentially with one pass through A.
        if(incy>0) then
           jy = 1
        else
           jy = 1-(n-1)*incy
        end if
        if(incx==1) then
           do j=1,n
              if(all(y(jy)/=ZERO)) then
                 temp = alpha*y(jy)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                 do i=1,m
                    a(i,j) = a(i,j)+x(i)*temp
                 end do
              end if
              jy = jy+incy
           end do
        else
           if(incx>0) then
              kx=1
           else
              kx = 1-(m-1)*incx
           end if
           do j=1,n
              if(all(y(jy)/=ZERO)) then
                 temp = alpha*y(jy)
                 ix = kx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                 do i=1,m
                    a(i,j) = a(i,j)+x(ix)*temp
                    ix = ix+incx
                 end do
              end if
              jy = jy+incy
           end do
        end if
        ! End of ZGERU
     end subroutine gms_zgeru

     subroutine gms_zhbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zhbmv
#endif
        character(len=1),                      intent(in),value    :: uplo
        integer(kind=int4),                    intent(in),value    :: n
        integer(kind=int4),                    intent(in),value    :: k
        type(AVX512c8f64_t),                   intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED a:64
#endif
        integer(kind=int4),                    intent(in),value    :: lda
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED x:64
#endif
        integer(kind=int4),                    intent(in),value    :: incx
        type(AVX512c8f64_t),                   intent(in)          :: beta
        type(AVX512c8f64_t), dimension(*),     intent(inout)       :: y
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED y:64
#endif
        integer(kind=int4),                    intent(in),value    :: incy
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
        type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
        type(AVX512c8f64_t), automatic :: temp2
        integer(kind=int4),  automatic :: i,info,ix,iy,j,jx,jy,kplus1,kx,ky,l
        logical(kind=int1),  automatic :: aeq0,beq1,beq0,bneq1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ONE
#endif
        type(AVX512c8f64_t), parameter :: ONE  = AVX512c8f64_t([1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                                                                1.0_dp,1.0_dp,1.0_dp,1.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        ! EXec code ....
        ! Test the input parameters.
        
        info = 1
        if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
           info = 1
        else if(n<0) then
           info = 2
        else if(k<0) then
           info = 3
        else if(lda<(k+1)) then
           info = 6
        else if(incx==0) then
           info = 8
        else if(incy==0) then
           info = 11
        end if
        if(info/=0) then
           call xerbla('GMS_ZHBMV',info)
           return
        end if
        ! Quick return if possible
        aeq0 = .false.
        beq1 = .false.
        aeq0 = all(alpha==ZERO)
        beq1 = all(beta==ONE)
        if((n==0) .or. ((aeq0) .and. (beq1))) return
        !  Set up the start points in  X  and  Y.
        if(incx>0) then
           kx = 1
        else
           kx = 1-(n-1)*incx
        end if
        if(incy>0) then
           ky = 1
        else
           ky = 1-(n-1)*incx
        end if
!      Start the operations. In this version the elements of the array A
!*     are accessed sequentially with one pass through A.
!*
        !1*     First form  y := beta*y.
        bneq1 = .false.
        beq0  = .false.
        bneq1 = all(beta/=ONE)
        beq0  = all(beta==ZERO)
        if(bneq1) then
           if(incy==1) then
              if(beq0) then
#if defined __INTEL_COMPILER
                   !DIR$ VECTOR ALIGNED
                 !DIR$ VECTOR ALWAYS
#endif
                   do i=1,n
                      y(i) = ZERO
                   end do
                else
#if defined __INTEL_COMPILER
                   !DIR$ VECTOR ALIGNED
                   !DIR$ VECTOR ALWAYS
#endif
                   do i=1,m
                      y(i) = beta*y(i)
                   end do
                end if
             else
                iy = ky
                if(beq0) then
#if defined __INTEL_COMPILER
                   !DIR$ VECTOR ALIGNED
                   !DIR$ VECTOR ALWAYS
#endif
                   do i=1,n
                      y(iy) = ZERO
                      iy = iy+incy
                   end do
                else
#if defined __INTEL_COMPILER
                   !DIR$ VECTOR ALIGNED
                   !DIR$ VECTOR ALWAYS
#endif
                   do i=1,n
                      y(iy) = beta*y(iy)
                      iy = iy+incy
                   end do
                end if
             end if
          end if
          if(aeq0) return
          if(lsame(uplo,'U')) then
             !  Form  y  when upper triangle of A is stored.
             kplus1 = k+1
             if((incx==1) .and. (incy==1)) then
                do j=1,n
                   temp1 = alpha*x(j)
                   temp2 = ZERO
                   l = kplus1-j
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do u=max(1,j-k),j-1
                      y(i) = y(i)+temp1*a(l+i,j)
                      temp2 = temp2+conjugate(a(l+i,j))*x(i)
                   end do
                   y(j) = y(j)+temp1*a(kplus1,j).re+alpha*temp2
                end do
             else
                jx = kx
                jy = ky
                do j=1,n
                   temp1 = alpha*x(jx)
                   temp2 = ZERO
                   ix = kx
                   iy = ky
                   l = kplus1-j
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do i=max(1,j-k),j-1
                      y(iy) = y(iy)+temp1*a(l+i,j)
                      temp2 = temp2+conjugate(a(l+i,j))*x(ix)
                      ix = ix+incx
                      iy = iy+incy
                   end if
                   y(jy) = y(jy)+temp1*a(kplus1,j).re+alpha*temp2
                   jx = jx+incx
                   jy = jy+incy
                   if(j>k) then
                      kx = kx+incx
                      ky = ky+incy
                   end if
                end do
             end if
          else
             !  Form  y  when lower triangle of A is stored.
             if((incx==1) .and. (incy==1)) then
                do j=1,n
                   temp1 = alpha*x(j)
                   temp2 = ZERO
                   y(j) = y(j)+temp1*a(1,j).re
                   l = 1-j
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do i=j+1,min(n,j+k)
                      y(i) = y(i)+temp1*a(l+i,j)
                      temp2 = temp2+conjugate(a(l+i,j))*x(i)
                   end do
                   y(j) = y(j)+alpha*temp2
                end do
             else
                jx = kx
                jy = ky
                do j=1,n
                   temp1 = alpha*x(jx)
                   temp2 = ZERO
                   y(jy) = y(jy)+temp1*a(1,j).re
                   l = 1-j
                   ix = jx
                   iy = jy
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do i=j+1,min(n,j+k)
                      ix = ix+incx
                      iy = iy+incy
                      y(iy) = y(iy) + temp1*a(l+i,j)
                      temp2 = temp2+conjugate(a(l+i,j))*x(ix)
                   end do
                   y(jy) = y(jy)+alpha*temp2
                   jx = jx+incx
                   jy = jy+incy
                end do
             end if
          end if
          ! End of zhbmv
      end subroutine gms_zhbmv

!        Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level3
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 3 Blas routine.
!*>
!*>  -- Written on 8-February-1989.
!*>     Jack Dongarra, Argonne National Laboratory.
!*>     Iain Duff, AERE Harwell.
!*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
      !*>     Sven Hammarling, Numerical Algorithms Group Ltd.
      !Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
!*> \endverbatim
!*>
     !*  =====================================================================

      subroutine gms_zhemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zhemm
#endif
          character(len=1),                      intent(in),value    :: side
          character(len=1),                      intent(in),value    :: uplo
          integer(kind=int4),                    intent(in),value    :: m
          integer(kind=int4),                    intent(in),value    :: n
          type(AVX512c8f64_t),                   intent(in)          :: alpha
          type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
#if defined __INTEL_COMPILER
          !DIR$ ASSUME_ALIGNED a:64
#endif
          integer(kind=int4),                    intent(in),value    :: lda
          type(AVX512c8f64_t), dimension(ldb,*), intent(in)          :: b
#if defined __INTEL_COMPILER
          !DIR$ ASSUME_ALIGNED b:64
#endif
          integer(kind=int4),                    intent(in),value    :: ldb
          type(AVX512c8f64_t),                   intent(in)          :: beta
          type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
#if defined __INTEL_COMPILER
          !DIR$ ASSUME_ALIGNED c:64
#endif
          integer(kind=int4),                    intent(in),value    :: ldc
          ! LOcals
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: temp1,temp2
#endif
          type(AVX512c8f64_t), automatic :: temp1,temp2
          integer(kind=int4),  automatic :: i,info,j,k,nrowa
          logical(kind=int4),  automatic :: upper
          logical(kind=int1),  automatic :: aeq0,beq1,beq0
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: ONE
#endif
          type(AVX512c8f64_t), parameter :: ONE  = AVX512c8f64_t([1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                                                                1.0_dp,1.0_dp,1.0_dp,1.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
          type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
          ! EXec code ....
          !  Set NROWA as the number of rows of A.
          if(lsame(side,'L')) then
             nrowa = m
          else
             nrowa = n
          end if
          upper = lsame(uplo,'U')
          !  Test the input parameters.
          aeq0 = .false.
          beq1 = .false.
          beq0 = .false.
          info = 0
          if((.not.lsame(side,'L')) .and. (.not.lsame(side,'R'))) then
             info = 1
          else if((.not.upper) .and. (not.lsame(uplo,'L'))) then
             info = 2
          else if(m<0) then
             info = 3
          else if(n<0) then
             info = 4
          else if(lda < max(1,nrowa)) then
             info = 7
          else if(ldb < max(1,m)) then
             info = 9
          else if(ldc < max(1,m)) then
             info = 13
          end if
          if(info/=0) then
             call xerbla('GMS_ZHEMM',info)
             return
          end if
          !  Quick return if possible.
          aeq0 = all(alpha==ZERO)
          beq1 = all(beta==ONE)
          if((m==0) .or. (n==0) .or. &
               ((aeq0) .and. (beq1))) return
          !  And when  alpha.eq.zero.
          beq0 = all(beta==ZERO)
          if(aeq0) then
             if(beq0) then
                do j=1,n
#if defined __INTEL_COMPILER
                   !DIR$ VECTOR ALIGNED
                   !DIR$ VECTOR ALWAYS
#endif
                   do i=1,m
                      c(i,j) = ZERO
                   end do
                end do
             else
                do j=1,n
#if defined __INTEL_COMPILER
                   !DIR$ VECTOR ALIGNED
                   !DIR$ VECTOR ALWAYS
#endif
                   do i=1,m
                      c(i,j) = beta*c(i,j)
                   end do
                end do
             end if
             return
          end if
          !  Start the operations.
          if(lsame(side,'L')) then
             !  Form  C := alpha*A*B + beta*C.
             if(upper) then
                do j=1,n
                   do i=1,m
                      temp1 = alpha*b(i,j)
                      temp2 = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                      do k=1,i-1
                         c(k,j) = c(k,j)+temp1*a(k,i)
                         temp2 = temp2+b(k,j)*conjugate(a(k,i))
                      end do
                      if(beq0) then
                         c(i,j) = temp1*a(i,i).re+alpha*temp2
                      else
                         c(i,j) = beta*c(i,j)+temp1*a(i,i).re+ &
                              alpha*temp2
                      end if
                   end do
                end do
             else
                do j=1,n
                   do i=m,1,-1
                      temp1 = alpha*b(i,j)
                      temp2 = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                      do k=i+1,m
                         c(k,j) = c(k,j)+temp1*a(k,i)
                         temp2 = temp2+b(k,j)*conjugate(a(k,i))
                      end do
                      if(beq0) then
                         c(i,j) = temp1*a(i,i).re+alpha*temp2
                      else
                         c(i,j) = beta*c(i,j)+temp1*a(i,i).re + &
                              alpha*temp2
                      end if
                   end do
                end do
             end if
          else
             !   Form  C := alpha*B*A + beta*C.
             do j=1,n
                temp1 = alpha*a(j,j).re
                if(beq0) then
                   !DIR$ VECTOR ALIGNED
                   !DIR$ VECTOR ALWAYS
                   do i=1,m
                      c(i,j) = temp*b(i,j)
                   end do
                else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do i=1,m
                      c(i,j) = beta*c(i,j)+temp1*b(i,j)
                   end do
                end if
                do k=1,j-1
                   if(upper) then
                      temp1 = alpha*a(k,j)
                   else
                      temp1 = alpha*conjugate(a(j,k))
                   end if
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do i=1,m
                      c(i,j) = c(i,j)+temp1*b(i,k)
                   end do
                end do
                do k=j+1,n
                   if(upper) then
                      temp1 = alpha*conjugate(a(j,k))
                   else
                      temp1 = alpha*a(k,j)
                   end if
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do i=1,m
                      c(i,j) = c(i,j)+temp1*b(i,k)
                   end do
                end do
             end do
          end if
          !End of ZHEMM
     end subroutine gms_zhemm

!      Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
     !*>     Richard Hanson, Sandia National Labs.
     !!Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
!*> \endverbatim
     !*>

     subroutine gms_zhemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zhemv
#endif
        character(len=1),                       intent(in),value    :: uplo
        integer(kind=int4),                     intent(in),value    :: n
        type(AVX512c8f64_t),                    intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(lda,*),  intent(in)          :: a
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED a:64
#endif
        integer(kind=int4),                     intent(in),value    :: lda
        type(AVX512c8f64_t), dimension(*),      intent(in)          :: x
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED x:64
#endif
        integer(kind=int4),                     intent(in),value    :: incx
        type(AVX512c8f64_t),                    intent(in)          :: beta
        type(AVX512c8f64_t), dimension(*),      intent(inout)       :: y
        integer(kind=int4),                     intent(in),value    :: incy
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
        type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
        type(AVX512c8f64_t), automatic :: temp2
        integer(kind=int4),  automatic :: i,info,ix,iy,j,jx,jy,kx,ky
        logical(kind=int1),  automatic :: aeq0,beq1,bneq1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ONE
#endif
        type(AVX512c8f64_t), parameter :: ONE  = AVX512c8f64_t([1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                                                                1.0_dp,1.0_dp,1.0_dp,1.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        !  Test the input parameters.
        info  = 0
        aeq0  = .false.
        beq1  = .false.
        bneq1 = .false.
        if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
           info = 1
        else if(n<0) then
           info = 2
        else if(lda < max(n,n)) then
           info = 5
        else if(incx==0) then
           info = 7
        else if(incy==0) then
           info = 10
        end if
        if(info/=0) then
           call xerbla('GMS_ZHEMV',info)
           return
        end if
        aeq0 = all(alpha==ZERO)
        beq1 = all(beta==ONE)
        !   Quick return if possible.
        if((n==0) .or. ((aeq0) .and. (beq1))) return
        !   Set up the start points in  X  and  Y.
        if(incx>0) then
           kx = 1
        else
           kx = 1-(n-1)*incx
        end if
        if(incy>0) then
           ky = 1
        else
           ky = 1-(n-1)*incy
        end if
        !  Start the operations. In this version the elements of A are
        !*     accessed sequentially with one pass through the triangular part
        !*     of A.
        !*
        bneq1 = all(beta/=ONE)
        !*     First form  y := beta*y.
        if(bneq1) then
           if(incy==1) then
              if(beq0) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                 !DIR$ VECTOR ALWAYS
#endif
                    do i=1,n
                       y(i) = ZERO
                    end do
                 else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#endif
                    do i=1,n
                       y(i) = beta*y(i)
                    end do
                 end if
              else
                 iy = ky
                 if(beq0) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#endif
                    do i=1,n
                       y(iy) = ZERO
                       iy = iy+incy
                    end do
                 else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#endif
                    do i=1,n
                       y(iy) = beta*y(iy)
                       iy = iy+incy
                    end do
                 end if
              end if
           end if
           if(aeq0) return
           if(lsame(uplo,'U')) then
              !Form  y  when A is stored in upper triangle.
              if((incx==1) .and. (incy==1)) then
                 do j=1,n
                    temp1 = alpha*x(j)
                    temp2 = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=1,j-1
                       y(i) = y(i)+temp1*a(i,j)
                       temp2 = temp2+conjugate(a(i,j))*x(i)
                    end do
                    y(j) = y(j)+temp1*a(j,j).re+alpha*temp2
                 end do
              else
                 jx = kx
                 jy = ky
                 do j=1,n
                    temp1 = alpha*x(jx)
                    temp2 = ZERO
                    ix = kx
                    iy = ky
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=1,j-1
                       y(iy) = y(iy)+temp1*a(i,j)
                       temp2 = temp2+conjugate(a(i,j))*x(ix)
                       ix = ix+incx
                       iy = iy+incy
                    end do
                    y(jy) = y(jy)+temp1*a(j,j).re+alpha*temp2
                    jx = jx+incx
                    jy = jy+incy
                 end do
              end if
           else
              !  Form  y  when A is stored in lower triangle.
              if((incx==1) .and. (incy==1)) then
                 do j=1,n
                    temp1 = alpha*x(j)
                    temp2 = ZERO
                    y(j) = y(j)+temp1*a(j,j).re
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=j+1,n
                       y(i) = y(i)+temp1*a(i,j)
                       temp2 = temp2+conjugate(a(i,j))*x(i)
                    end do
                    y(j) = y(j)+alpha*temp2
                 end do
              else
                 jx = kx
                 jy = ky
                 do j=1,n
                    temp1 = alpha*x(jx)
                    temp2 = ZERO
                    y(jy) = y(jy)+temp1*a(j,j).re
                    ix = jx
                    iy = jy
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=j+1,n
                       ix = ix+incx
                       iy = iy+incy
                       y(iy) = y(iy)+temp1*a(i,j)
                       temp2 = temp2+conjugate(a(i,j))*x(ix)
                    end do
                    y(jy) = y(jy)+alpha*temp2
                    jx = jx+incx
                    jy = jy+incy
                 end do
              end if
           end if
           !End of ZHEMV
     end subroutine gms_zhemv

!      Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>
!1*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
     !*>     Richard Hanson, Sandia National Labs.
     ! !!Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
!*> \endverbatim
!*>
     subroutine gms_zher(uplo,n,alpha,x,incx,a,lda)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zher
#endif
        character(len=1),                      intent(in),value    :: uplo
        integer(kind=int4),                    intent(in),value    :: n
        type(ZMM8r4_t),                        intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED x:64
#endif
        integer(kind=int4),                    intent(in),value    :: incx
        type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED a:64
#endif
        integer(kind=int4),                    intent(in),value    :: lda
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
        integer(kind=int4),  automatic :: i,info,ix,j,jx,kx
        logical(kind=int1),  automatic :: aeq0
        !
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        aeq0 = .false.
        info = 0
        if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
           info = 1
        else if(n<0) then
           info = 2
        else if(incx==0) then
           info = 5
        else if(lda < max(1,n)) then
           info = 7
        end if
        if(info/=0) then
           call xerbla('GMS_ZHER',info)
           return
        end if
        !  Quick return if possible.
        aeq0 = all(alpha==ZERO.re)
        if((n==0) .or. (aeq0)) return
        !  Set the start point in X if the increment is not unity.
        if(incx<=0) then
           kx = 1-(n-1)*incx
        else
           kx = 1
        end if
        !   Start the operations. In this version the elements of A are
        !*     accessed sequentially with one pass through the triangular part
        !*     of A.
        if(lsame(uplo,'U')) then
           !   Form  A  when A is stored in upper triangle.
           if(incx==1) then
              do j=1,n
                 if(all(x(j)/=ZERO)) then
                    temp = alpha*conjugate(x(j))
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=1,j-1
                       a(i,j) = a(i,j)+x(i)*temp
                    end do
                    a(j,j) = a(j,j).re+x(j).re*temp.re
                 else
                    a(j,j) = a(j,j).re
                 end if
              end do
           else
              jx = kx
              do j=1,n
                 if(all(x(jx)/=ZERO)) then
                    temp = alpha*conjugate(x(jx))
                    ix = kx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=1,j-1
                       a(i,j) = a(i,j)+x(ix)*temp
                       ix = ix+incx
                    end do
                    a(j,j) = a(j,j).re+x(j).re*temp.re
                 else
                    a(j,j) = a(j,j).re
                 end if
              end do
           else
              jx = kx
              do j=1,n
                 if(all(x(jx)/=ZERO)) then
                    temp = alpha*conjugate(x(jx))
                    ix = kx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=1,j-1
                       a(i,j) = a(i,j)+x(ix)*temp
                       ix = ix+incx
                    end do
                    a(j,j) = a(j,j).re+x(jx).re*temp.re
                 else
                    a(j,j) = a(j,j).re
                 end if
                 jx = jx+incx
              end do
           end if
        else
           !  Form  A  when A is stored in lower triangle.
           if(incx==1) then
              do j=1,n
                 if(all(x(j)/=ZERO)) then
                    temp = alpha*conjugate(x(j))
                    a(j,j) = a(j,j).re+temp.re*x(j).re
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                    do i=j+1,n
                       a(i,j) = a(i,j)+x(i)*temp
                    end do
                 else
                    a(j,j) = a(j,j).re
                 end if
              end do
           else
              jx = kx
              do j=1,n
                 if(all(x(jx)/=ZERO)) then
                    temp = alpha*conjugate(x(jx))
                    a(j,j) = a(j,j).re+temp.re*x(jx).re
                    ix = jx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                  
                    do i=j+1,n
                       ix = ix+incx
                       a(i,j) = a(i,j)+x(ix)*temp
                    end do
                 else
                    a(j,j) = a(j,j).re
                 end if
                 jx = jx+incx
              end do
           end if
        end if
        ! End of ZHER
    end subroutine gms_zher

!    Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!1*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2017
!*
!*> \ingroup complex16_blas_level1
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!1*>
!*>     jack dongarra, 3/11/78.
!*>     modified 3/93 to return if incx .le. 0.
    !*>     modified 12/3/93, array(1) declarations changed to array(*)
    !  !!Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
    !*> \endverbatim

    subroutine gms_zscal(n,za,zx,incx)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zscal
      !DIR$ ATTIRBUTES VECTOR  :: gms_zscal
#endif
       integer(kind=int4),                     intent(in),value    :: n
       type(AVX512c8f64_t),                    intent(in)          :: za
       type(AVX512c8f64_t), dimension(*),      intent(inout)       :: zx
#if defined __INTEL_COMPILER
       !DIR$ ASSUME_ALIGNED zx:64
#endif
       integer(kind=int4),                     intent(in),value    :: incx
       ! Locals
       integer(kind=int4), automatic :: i,nincx
       ! Exec code ....
       if(n<=0 .or. incx<0) return
       if(incx==1) then
          !   code for increment equal to 1
#if defined __INTEL_COMPILER
          !DIR$ VECTOR ALIGNED
          !DIR$ VECTOR ALWAYS
#endif
         
          !GCC$ UNROLL(4)
          do i=1,n
             zx(i) = za*zx(i)
          end do
       else
          !  code for increment not equal to 1
          nincx = n*incx
#if defined __INTEL_COMPILER
          !DIR$ VECTOR ALIGNED
          !DIR$ VECTOR ALWAYS
          do i=1,nincx,incx
             zx(i) = za*zx(i)
          end do
       end if 
     end subroutine gms_zscal

!      Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!1*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!1*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
     !1*>     Richard Hanson, Sandia National Labs.
     ! !!Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
     !*> \endverbatim

     subroutine gms_zher2(uplo,n,alpha,x,incx,y,incy,a,lda)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 64 :: gms_zher2
#endif
         character(len=1),                      intent(in),value    :: uplo
         integer(kind=int4),                    intent(in),value    :: n
         type(AVX512c8f64_t),                   intent(in)          :: alpha
         type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED x:64
#endif
         integer(kind=int4),                    intent(in),value    :: incx
         type(AVX512c8f64_t), dimension(*),     intent(in)          :: y
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED y:64
#endif
         integer(kind=int4),                    intent(in),value    :: incy
         type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED a:64
#endif
         integer(kind=int4),                    intent(in),value    :: lda
         ! LOcals ....
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
         type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
         type(AVX512c8f64_t), automatic :: temp2
         integer(kind=int4),  automatic :: i,info,ix,iy,j,jx,jy,kx,ky
         logical(kind=int1),  automatic :: aeq0
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
         type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
         ! EXec code ....
         aeq0 = .false.
         info = 0
         if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
            info = 1
         else if(n<0) then
            info = 2
         else if(incx==0) then
            info = 5
         else if(incy==0) then
            info = 7
         else if(lda < max(1,n)) then
            info = 9
         end if
         if(info/=0) then
            call xerbla('GMS_ZHER2',info)
            return
         end if
         !  Quick return if possible.
         aeq0 = all(alpha==ZERO)
         if((n==0) .or. (aeq0)) return
         !  Set up the start points in X and Y if the increments are not both
         !*     unity.
         if((incx/=1) .or. (incy/=1)) then
            if(incx>0) then
               kx = 1
            else
               kx = 1-(n-1)*incx
            end if
            if(incy.0) then
               ky = 1
            else
               ky = 1-(n-1)*incy
            end if
            jx = kx
            jy = ky
         end if
         !  Start the operations. In this version the elements of A are
         !*     accessed sequentially with one pass through the triangular part
         !*     of A.
         if(lsame(uplo,'U')) then
            !  Form  A  when A is stored in the upper triangle.
            if((incx==1) .and. (incy==1)) then
               do j=1,n
                  if((all(x(j)/=ZERO)) .or. (all(y(j)/=ZERO))) then
                     temp1 = alpha*conjugate(y(j))
                     temp2 = conjugate(alpha*x(j))
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                     do i=1,j-1
                        a(i,j) = a(i,j)+x(i)*temp1+y(i)*temp2
                     end do
                     a(j,j) = a(j,j).re+x(j).re*temp1.re+x(j).re*temp2.re
                  else
                     a(j,j) = a(j,j).re
                  end if
               end do
            else
               do j=1,n
                  if((all(x(jx)/=ZERO)) .or. (all(y(jy)/=ZERO))) then
                     temp1 = alpha*conjugate(y(jy))
                     temp2 = conjugate(alpha*x(jx))
                     ix = kx
                     iy = ky
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                    
                     do i=1,j-1
                        a(i,j) = a(i,j)+x(ix)*temp1+y(iy)*temp2
                        ix = ix+incx
                        iy = iy+incy
                     end do
                     a(j,j) = a(j,j).re+x(jx).re*temp1.re+y(jy).re*temp2.re
                  else
                     a(j,j) = a(j,j).re
                  end if
                  jx = jx+incx
                  jy = jy+incy
               end do
            end if
         else
            !  Form  A  when A is stored in the lower triangle.
            if((incx==1) .and. (incy==1)) then
               do j=1,n
                  if((all(x(j)/=ZERO)) .or. (all(y(j)/=ZERO))) then
                     temp1 = alpha*conjugate(y(j))
                     temp2 = conjugate(alpha*x(j))
                     a(j,j) = a(j,j).re+ &
                          x(j).re*temp1.re+y(j).re*temp2.re
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                 
                     do i=j+1,n
                        a(i,j) = a(i,j)+x(i)*temp1+y(i)*temp2
                     end do
                  else
                     a(j,j) = a(j,j).re
                  end if
               end do
            else
               do j=1,n
                  if((all(x(jx)/=ZERO)) .or. (all(y(jx)/=ZERO))) then
                     temp1 = alpha*conjugate(y(jy))
                     temp2 = conjugate(alpha*x(jx))
                     a(j,j) = a(j,j).re + &
                          x(jx).re*temp1.re+y(jy).re*temp2.re
                     ix = jx
                     iy = jy
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                 
                     do i=j+1,n
                        ix = ix+incx
                        iy = iy+incy
                        a(i,j) = a(i,j)+x(ix)*temp1+y(iy)*temp2
                     end do
                  else
                     a(j,j) = a(j,j).re
                  end if
                  jx = jx+incx
                  jy = jy+incy
               end do
            end if
         end if
         ! End of ZHER2
       end subroutine gms_zher2

!        Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level3
!*
!*> \par Further Details:
!*  =====================
!1*>
!*> \verbatim
!*>
!*>  Level 3 Blas routine.
!*>
!*>  -- Written on 8-February-1989.
!*>     Jack Dongarra, Argonne National Laboratory.
!*>     Iain Duff, AERE Harwell.
!*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!*>     Sven Hammarling, Numerical Algorithms Group Ltd.
!*>
!*>  -- Modified 8-Nov-93 to set C(J,J) to DBLE( C(J,J) ) when BETA = 1.
       !*>     Ed Anderson, Cray Research Inc.
       !!Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
!*> \endverbatim
       !*>

       subroutine gms_zher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zher2k
#endif
          use mod_vecconst, only : v8_n1
          character(len=1),                      intent(in),value    :: uplo
          character(len=1),                      intent(in),value    :: trans
          integer(kind=int4),                    intent(in),value    :: n
          integer(kind=int4),                    intent(in),value    :: k
          type(AVX512c8f64_t),                   intent(in)          :: alpha
          type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
#if defined __INTEL_COMPILER
          !DIR$ ASSUME_ALIGNED a:64
#endif
          integer(kind=int4),                    intent(in),value    :: lda
          type(AVX512c8f64_t), dimension(ldb,*), intent(in)          :: b
#if defined __INTEL_COMPILER
          !DIR$ ASSUME_ALIGNED b:64
#endif
          integer(kind=int4),                    intent(in),value    :: ldb
          type(ZMM8r8_t),                        intent(in)          :: beta
          type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
#if defined __INTEL_COMPILER
          !DIR$ ASSUME_ALIGNED c:64
#endif
          integer(kind=int4),                    intent(in),value    :: ldc
          ! Locals
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
          type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
          type(AVX512c8f64_t), automatic :: temp2
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: temp3
#endif
          type(AVX512c8f64_t), automatic :: temp3
          integer(kind=int4),  automatic :: i,info,j,l,nrowa
          logical(kind=int4),  automatic :: upper
          logical(kind=int1),  automatic :: aeq0,beq1,beq0,bneq1
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
          type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
          ! EXec code ....
          ! Test the input parameters.
          
          aeq0  = .false.
          beq1  = .false.
          beq0  = .false.
          bneq1 = .false.
          if(lsame(trans,'N')) then
             nrowa = n
          else
             nrowa = k
          end if
          upper = lsame(uplo,'U')
          info = 0
          if((.not.upper) .and. (.not.lsame(uplo,'L'))) then
             info = 1
          else if((.not.lsame(trans,'N')) .and. &
               (.not.lsame(trans,'C'))) then
             info = 2
          else if(n<0) then
             info = 3
          else if(k<0) then
             info = 4
          else if(lda < max(1,nrowa)) then
             info = 7
          else if(ldb < max(1,nrowa)) then
             info = 9
          else if(ldc < max(1,n)) then
             info = 12
          end if
          if(info/=0) then
             call xerbla('GMS_ZHER2K',info)
             return
          end if
          !  Quick return if possible.
          aeq0 = all(alpha==ZERO)
          beq1 = all(beta.v==v8_n1.v)
          if((n==0) .or. ((aeq0) .or. &
               (k==0) .and. (beq1))) return
          !  And when  alpha.eq.zero.
          beq0 = all(beta.v==ZERO.re)
          if(aeq0) then
             if(upper) then
                if(beq0)  then
                   do j=1,n
                      !DIR$ VECTOR ALIGNED
                      !DIR$ VECTOR ALWAYS
                      do i=1,j
                         c(i,j) = ZERO
                      end do
                   end do
                else
                   do j=1,n
                      !DIR$ VECTOR ALIGNED
                      !DIR$ VECTOR ALWAYS
                      do i=1,j-1
                         c(i,j) = beta*c(i,j)
                      end do
                      c(j,j) = beta.v*c(i,j).re
                   end do
                end if
             else
                if(beq0) then
                   do j=1,n
                      !DIR$ VECTOR ALIGNED
                      !DIR$ VECTOR ALWAYS
                      do i=j,n
                         c(i,j) = ZERO
                      end do
                   end do
                else
                   do j=1,n
                      c(j,j) = beta.v*c(j,j).re
                      !DIR$ VECTOR ALIGNED
                      !DIR$ VECTOR ALWAYS
                      do i=j+1,n
                         c(i,j) = beta*c(i,j)
                      end do
                   end do
                end if
             end if
             return
          end if
          !  Start the operations.
          bneq1 = all(beta.v/=v8_n1.v)
          if(lsame(trans,'N')) then
             !   Form  C := alpha*A*B**H + conjg( alpha )*B*A**H +
             !*                   C.
             if(upper) then
                do j=1,n
                   if(beq0) then
                      !DIR$ VECTOR ALIGNED
                      !DIR$ VECTOR ALWAYS
                      do i=1,j
                         c(i,j) = ZERO
                      end do
                   else if(bneq1) then
                      !DIR$ VECTOR ALIGNED
                      !DIR$ VECTOR ALWAYS
                      do i=1,j-1
                         c(i,j) = beta*c(i,j)
                      end do
                      c(j,j) = beta.v*c(j,j).re
                   else
                      c(j,j) = c(j,j).re
                   end if
                   do l=1,k
                      if((all(a(j,l)/=ZERO)) .or. (all(b(j,l)/=ZERO))) then
                         temp1 = alpha*conjugate(b(j,l))
                         temp2 = conjugate(alpha*a(j,l))
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                       
                         do i=1,j-1
                            c(i,j) = c(i,j)+a(i,l)*temp1+ &
                                 b(i,l)*temp2
                         end do
                         c(j,j) = c(j,j).re+a(j,l).re*temp1.re+&
                              b(j,l).re*temp2.re
                      end if
                   end do
                end do
             else
                do j=1,n
                   if(beq0) then
                      !DIR$ VECTOR ALIGNED
                      !DIR$ VECTOR ALWAYS
                      do i=j,n
                         c(i,j) = ZERO
                      end do
                   else if(bneq1) then
                      !DIR$ VECTOR ALIGNED
                      !DIR$ VECTOR ALWAYS
                      do i=j+1,n
                         c(i,j) = beta*c(i,j)
                      end do
                      c(j,j) = beta.v*c(j,j).re
                   else
                      c(j,j) = c(j,j).re
                   end if
                   do l=1,k
                      if((all(a(j,l)/=ZERO)) .or. (all(b(j,l)/=ZERO))) then
                         temp1 = alpha*conjugate(b(j,l))
                         temp2 = conjugate(alpha*a(j,l))
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                      
                         do i=j+1,n
                            c(i,j) = c(i,j)+a(i,l)*temp1+&
                                 b(i,l)*temp2
                         end do
                         c(j,j) = c(j,j).re+a(j,l).re*temp1.re+&
                              b(j,l).re*temp2.re
                      end if
                   end do
                end do
             end if
          else
             !   Form  C := alpha*A**H*B + conjg( alpha )*B**H*A +
             if(upper) then
                temp3 = default_init()
                do j=1,n
                   do i=1,j
                      temp1 = ZERO
                      temp2 = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                   
                      do l=1,k
                         temp1 = temp1+conjugate(a(l,i))*b(l,j)
                         temp2 = temp2+conjugate(b(l,i))*a(l,j)
                      end do
                      if(i==j) then
                         if(beq0) then
                            temp3=conjugate(alpha)
                            c(j,j) = alpha.re*temp1.re+temp3.re* &
                                     temp2.re
                         else
                            temp3=conjugate(alpha)
                            c(j,j) = beta.v*c(j,j).re+alpha.re* &
                                 temp1.re+temp3.re*temp2.re
                         end if
                      else
                         if(beq0) then
                            c(i,j) = alpha*temp1+conjugate(alpha)*temp2
                         else
                            c(i,j) = beta*c(i,j)+alpha*temp1+ &
                                 conjugate(alpha)*temp2
                         end if
                      end if
                   end do
                end do
             else
                do j=1,n
                   do i=1,n
                      temp1 = ZERO
                      temp2 = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                   
                      do l=1,k
                         temp1 = temp1+conjugate(a(l,i))*b(l,j)
                         temp2 = temp2+conjugate(b(l,i))*a(l,j)
                      end do   
                      if(i==j) then
                         if(beq0) then
                            temp3=conjugate(alpha)
                            c(j,j) = alpha.re+temp1.re+temp3.re* &
                                 temp2.re
                         else
                            temp3=conjugate(alpha)
                            c(j,j) = beta.v*c(j,j).re+alpha.re*temp1.re+ &
                                 temp3.re*temp2.re
                         end if
                      else
                         if(beq0) then
                            c(i,j) = alpha*temp1+conjugate(alpha)*temp2
                         else
                            c(i,j) = beta*c(i,j)+alpha*temp1+ &
                                 conjugate(alpha)*temp2
                         end if
                      end if
                   end do
                end do
             end if
          end if
          ! End of ZHER2K
     end subroutine gms_zher2k

!      Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level3
!*
!*> \par Further Details:
!*  =====================
!1*>
!*> \verbatim
!*>
!*>  Level 3 Blas routine.
!*>
!*>  -- Written on 8-February-1989.
!*>     Jack Dongarra, Argonne National Laboratory.
!*>     Iain Duff, AERE Harwell.
!*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!*>     Sven Hammarling, Numerical Algorithms Group Ltd.
!*>
!*>  -- Modified 8-Nov-93 to set C(J,J) to DBLE( C(J,J) ) when BETA = 1.
     !*>     Ed Anderson, Cray Research Inc.
     !Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
     !*> \endverbatim

     subroutine gms_zherk(uplo,trans,n,alpha,a,lda,beta,c,ldc)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zherk
#endif
        use mod_vecconst, only : v8_n1,v8_n0
        character(len=1),                      intent(in),value    :: uplo
        character(len=1),                      intent(in),value    :: trans
        integer(kind=int4),                    intent(in),value    :: n
        type(ZMM8r8_t),                        intent(in)          :: alpha
        type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED a:64
#endif
        integer(kind=int4),                    intent(in),value    :: lda
        type(ZMM8r8_t),                        intent(in)          :: beta
        type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
#if defined __INTEL_COMPILER
        !DIR$ ASUME_ALIGNED c:64
#endif
        integer(kind=int4),                    intent(in),value    :: ldc
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: rtemp
#endif
        type(ZMM8r8_t),      automatic :: rtemp
        integer(kind=int4),  automatic :: i,info,j,l,nrowa
        logical(kind=int4),  automatic :: upper
        logical(kind=int1),  automatic :: aeq0,beq1,beq0,bneq1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        ! EXec code....
        if(lsame(trans,'N')) then
           nrowa = n
        else
           nrowa = k
        end if
        upper = lsame(uplo,'U')
        info = 0
        if((.not.upper) .and. (.not.lsame(uplo,'L'))) then
           info = 1
        else if((.not.lsame(trans,'N')) .and. &
             (.not.lsame(trans,'C'))) then
           info = 2
        else if(n<0) then
           info = 3
        else if(k<0) then
           info = 4
        else if(lda < max(1,nrowa)) then
           info = 7
        else if(ldc < max(1,n)) then
           info = 10
        endif
        if(info/=0) then
           call xerbla('GMS_ZHERK',info)
           return
        end if
        aeq0 = .false.
        beq1 = .false.
        aeq0 = all(alpha.v==v8_n0.v)
        beq1 = all(beta.v==v8_n1.v)
        if((n==0) .or. (aeq0) .or. &
             ((k==0) .and. (beq1))) return
        !  And when  alpha.eq.zero.
        beq0 = .false.
        beq0 = all(beta.v==v8_n0.v)
        if(aeq0) then
           if(upper) then
              if(beq0)  then
                 do j=1,n
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    do i=1,j
                       c(i,j) = ZERO
                    end do
                 end do
              else
                 do j=1,n
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    do i=1,j-1
                       c(i,j) = beta*c(i,j)
                    end do
                    c(j,j) = beta.v*c(j,j).re
                 end do
              end if
           else
              if(beq0) then
                 do j=1,n
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    do i=j,n
                       c(i,j) = ZERO
                    end do
                 end do
              else
                 do j=1,n
                    c(j,j) = beta.v*c(j,j).re
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    do i=j+1,n
                       c(i,j) = beta*c(i,j)
                    end do
                 end do
              end if
           end if
           return
        end if
        !  Start the operations.
        if(lsame(trans,'N')) then
           bneq1 = .false.
           bneq1 = all(beta.v/=v8_n1.v)
           !   Form  C := alpha*A*A**H + beta*C.
           if(upper) then
              do j=1,n
                 if(beq0) then
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    do i=1,j
                       c(i,j) = ZERO
                    end do
                 else if(bneq1) then
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    do i=1,j-1
                       c(i,j) = beta*c(i,j)
                    end do
                    c(j,j) = beta.v*c(j,j).re
                 else
                    c(j,j) = c(j,j).re
                 end if
                 do l=1,k
                    if((all(a(j,l)/=ZERO))) then
                       temp = alpha*conjugate(a(j,l))
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=1,j-1
                          c(i,j) = c(i,j)+temp*a(i,l)
                       end do
                       c(j,j) = c(j,j).re+temp.re*a(i,l).re
                    end if
                 end do
              end do
           else
              do j=1,n
                 if(beq0) then
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    do i=j,n
                       c(i,j) = ZERO
                    end do
                 else if(bneq1) then
                    c(j,j) = beta.v*c(j,j).re
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    do i=j+1,n
                       c(i,j) = beta*c(i,j)
                    end do
                 else
                    c(j,j) = c(j,j).re
                 end if
                 do l=1,k
                    if(all(a(j,l)/=ZERO)) then
                       temp = alpha*conjugate(a(j,l))
                       c(j,j) = c(j,j).re+temp.re*a(j,l).re
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=j+1,n
                          c(i,j) = c(i,j)+temp*a(i,l)
                       end do
                    end if
                 end do
               end do
            end if
         else
            ! Form  C := alpha*A**H*A + beta*C.
            if(upper) then
               do j=1,n
                  do i=1,j-1
                     temp = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                     do l=1,k
                        temp = temp+conjugate(a(l,i))*a(l,j)
                     end do
                     if(beq0) then
                        c(i,j) = alpha*temp
                     else
                        c(i,j) = alpha*temp+beta*c(i,j)
                     end if
                  end do
                  rtemp = v8_n0
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#endif
                  do l=1,k
                     rtemp = rtemp+conjugate(a(l,j))*a(l,j)
                  end do
                  if(beq0) then
                     c(j,j) = alpha*rtemp
                  else
                     c(j,j) = alpha*rtemp+beta.v*c(j,j).re
                  end if
                  do i=j+1,n
                     temp = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                     do l=1,k
                        temp = temp+conjugate(a(l,i))*a(l,j)
                     end do
                     if(beq0) then
                        c(i,j) = alpha*temp
                     else
                        c(i,j) = alpha*temp+beta*c(i,j)
                     end if
                  end do
               end do
            end if
         end if
         !End of ZHERK
     end subroutine gms_zherk

!      Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
     !*>     Richard Hanson, Sandia National Labs.
     !Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features) 
!*> \endverbatim
     !*>

     subroutine gms_zhpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zhpmv
#endif
         character(len=1),                     intent(in),value :: uplo
         integer(kind=int4),                   intent(in),value :: n
         type(AVX512c8f64_t),                  intent(in)       :: alpha
         type(AVX512c8f64_t), dimension(*),    intent(in)       :: ap
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED ap:64
#endif
         type(AVX512c8f64_t), dimension(*),    intent(in)       :: x
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED x:64
#endif
         integer(kind=int4),                   intent(in),value :: incx
         type(AVX512c8f64_t),                  intent(in)       :: beta
         type(AVX512c8f64_t), dimension(*),    intent(inout)    :: y
         integer(kind=int4),                   intent(in),value :: incy
         ! Locals
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
         type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
         type(AVX512c8f64_t), automatic :: temp2
         integer(kind=int4),  automatic :: i,info,ix,iy,j,jx,jy,k,kk,kx,ky
         logical(kind=int1),  automatic :: aeq0,beq1,bneq1,beq0
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: ONE
#endif
         type(AVX512c8f64_t), parameter :: ONE  = AVX512c8f64_t([1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                                                                1.0_dp,1.0_dp,1.0_dp,1.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
         type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
         ! Exec code .....
         info = 0
         if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
            info = 1
         else if(n<0) then
            info = 2
         else if(incx==0) then
            info = 6
         else if(incy==0) then
            info = 9
         end if
         if(info/=0) then
            call xerbla('GMS_ZHPMV',info)
            return
         end if
         aeq0 = .false.
         beq1 = .false.
         aeq0 = all(alpha==ZERO)
         beq1 = all(beta==ONE)
         !   Quick return if possible.
         if((n==0) .or. (aeq0 .and. beq1)) return
         !  Set up the start points in  X  and  Y.
         if(incx>0) then
            kx = 1
         else
            kx = 1-(n-1)*incx
         end if
         if(incy>0) then
            ky = 1
         else
            ky = 1-(n-1)*incy
         end if
         !      Start the operations. In this version the elements of the array AP
         !*     are accessed sequentially with one pass through AP.
         !*
         !*     First form  y := beta*y.
         bneq1 = .false.
         bneq1 = all(beta/=ONE)
         beq0  = .false.
         beq0  = all(beta==ZERO)
         if(bneq1) then
            if(incy==1) then
                if(beq0)  then
                  
                   !DIR$ VECTOR ALIGNED
                   !DIR$ VECTOR ALWAYS
                   !DIR$ UNROLL(4)
                   do i=1,n
                      y(i) = ZERO
                   end do
                else
                   !DIR$ VECTOR ALIGNED
                   !DIR$ VECTOR ALWAYS
                   !DIR$ UNROLL(10)
                   do i=1,n
                      y(i) = beta*y(i)
                   end do
                end if
             else
                iy = ky
                if(beq0) then
                   !DIR$ VECTOR ALIGNED
                   !DIR$ VECTOR ALWAYS
                   !DIR$ UNROLL(4)
                   do i=1,n
                      y(iy) = ZERO
                      iy = iy+incy
                   end do
                else
                   !DIR$ VECTOR ALIGNED
                   !DIR$ VECTOR ALWAYS
                   !DIR$ UNROLL(10)
                   do i=1,n
                      y(iy) = beta*y(iy)
                      iy = iy+incy
                   end do
                end if
             end if
          end if
          if(aeq0) return
          kk = 1
          if(lsame(uplo,'U')) then
             !  Form  y  when AP contains the upper triangle.
             if((incx==1) .and. (incy==1)) then
               
                do j=1,n
                   temp1 = alpha*x(j)
                   temp2 = ZERO
                   k = kk
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do i=1,j-1
                      y(i) = y(i)+temp1*ap(k)
                      temp2 = temp2+conjugate(ap(k))*x(i)
                      k = k+1
                   end do
                   y(j) = y(j)+temp1*ap(kk+j-1).re+alpha*temp2
                   kk = kk+j
                end do
             else
                jx = kx
                jy = ky
                do j=1,n
                   temp1 = alpha*x(jx)
                   temp2 = ZERO
                   ix = kx
                   iy = ky
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do k=kk,kk+j-2
                      y(iy) = y(iy)+temp1*ap(k)
                      temp2 = temp2+conjugate(ap(k))*x(ix)
                      ix = ix+incx
                      iy = iy+incy
                   end do
                   y(jy) = y(jy)+temp1*ap(kk+j-1).re+alpha*temp2
                   jx = jx+incx
                   jy = jy+incy
                   kk = kk+j
                end do
             end if
          else
             !   Form  y  when AP contains the lower triangle.
             if((incx==1) .and. (incy==1)) then
               
                do j=1,n
                   temp1 = alpha*x(j)
                   temp2 = ZERO
                   y(j) = y(j)+temp1*ap(kk).re
                   k = kk+1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do i=j+1,n
                      y(i) = y(i)+temp1*ap(k)
                      temp2 = temp2+conjugate(ap(k))*x(i)
                      k = k+1
                   end do
                   y(j) = y(j)+alpha*temp2
                   kk = kk+(n-j+1)
                end do
             else
                jx = kx
                jy = ky
                do j=1,n
                   temp1 = alpha*x(jx)
                   temp2 = ZERO
                   y(jy) = y(jy)+temp1*ap(kk).re
                   ix = jx
                   iy = jy
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do k=kk+1,kk+n-j
                      ix = ix+incx
                      iy = iy+incy
                      y(iy) = y(iy)+temp1*ap(k)
                      temp2 = temp2+conjugate(ap(k))*x(ix)
                   end do
                   y(jy) = y(jy)+alpha*temp2
                   jx = jx+incx
                   jy = jy+incy
                   kk = kk+(n-j+1)
                end do
             end if
          end if
          !End of ZHPMV
     end subroutine gms_zhpmv

!     *  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
     !*>     Richard Hanson, Sandia National Labs.
     !Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)
     !*> \endverbatim

     subroutine gms_zhpr(uplo,n,alpha,x,incx,ap)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zhpr
#endif
         use mod_vecconsts, only : v8_n0
         character(len=1),                  intent(in),value :: uplo
         integer(kind=int4),                intent(in),value :: n
         type(ZMM8r8_t),                    intent(in)       :: alpha
         type(AVX512c8f64_t), dimension(*), intent(in)       :: x
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED x:64
#endif
         integer(kind=int4),                intent(in),value :: incx
         type(AVX512c8f64_t), dimension(*), intent(inout)    :: ap
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED ap:64
#endif
         ! Locals
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
         type(AVX512c8f64_t), automatic :: temp
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: vtemp
#endif
         type(ZMM8r8_t),      automatic :: vtemp
         integer(kind=int4),  automatic :: i,info,ix,j,jx,k,kk,kx
         integer(kind=int1),  automatic :: aeq0
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
         type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
         ! Test the input parameters.
         info = 0
         if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
            info = 1
         else if(n<0) then
            info = 2
         else if(incx==0) then
            info = 5
         end if
         if(info/=0) then
            call xerbla('GMS_ZHPR',info)
            return
         end if
         !  Quick return if possible.
         aeq0 = .false.
         aeq0 = all(alpha.v==ZERO.re)
         if((n==0) .or. (aeq0)) return
         ! Set the start point in X if the increment is not unity.
         if(incx<=0) then
            kx = 1-(n-1)*incx
         else if(incx/=1) then
            kx = 1
         end if
         !   Start the operations. In this version the elements of the array AP
         !*     are accessed sequentially with one pass through AP.
         kk = 1
         vtemp = v8_n0
         if(lsame(uplo,'U')) then
            !  Form  A  when upper triangle is stored in AP.
            if(incx==1) then
               do j=1,n
                  if(all(x(j)/=ZERO)) then
                     vtemp = conjugate(x(j))
                     temp = zmm8r81x_init(alpha.v*vtemp.v)
                     k = kk
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                     do i=1,j-1
                        ap(k) = ap(k)+x(i)*temp
                        k = k+1
                     end do
                     ap(kk+j-1).re = ap(kk+j-1).re+x(j).re*temp.v
                  else
                     ap(kk+j-1).re = ap(kk+j-1).re
                  end if
                  kk = kk+j
               end do
            else
               jx = kx
               do j=1,n
                  if(all(x(jx)/=ZERO)) then
                     vtemp = conjugate(x(j))
                     ix = kx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                     do k=kk,kk+j-2
                        ap(k) = ap(k)+x(ix)*temp
                        ix = ix+incx
                     end do
                     ap(kk+j-1).re = ap(kk+j-1).re+x(jx).re*temp.v
                  else
                     ap(kk+j-1).re = ap(kk+j-1).re
                  end if
                  jx = jx+incx
                  kk = kk+j
               end do
            end if
         else
            ! Form  A  when lower triangle is stored in AP.
            if(incx==1) then
               do j=1,n
                  if(all(x(j)/=ZERO)) then
                     vtemp = conjugate(x(j))
                     temp = zmm8r81x_init(alpha.v*vtemp.v)
                     ap(kk).re = ap(kk).re+temp.v*x(j).re
                     k = kk+1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                     do i=j+1,n
                        ap(k) = ap(k)+x(i)*temp
                        k = k+1
                     end do
                  else
                     ap(kk).re = ap(kk).re
                  end if
                  kk = kk+n-j+1
               end do
            else
               jx = kx
               do j=1,n
                  if(all(x(jx)/=ZERO)) then
                     vtemp = conjugate(x(j))
                     temp = zmm8r81x_init(alpha.v*vtemp.v)
                     ap(kk).re = app(kk).re+temp.v*x(jx).re
                     ix = jx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                     do k=kk+1,kk+n-j
                        ix = ix+incx
                        ap(k) = ap(k)+x(ix)*temp
                     end do
                  else
                     ap(kk).re = ap(kk).re
                  end if
                  jx = jx+incx
                  kk = kk+n-j+1
               end do
            end if
         end if
         ! End of ZHPR
     end subroutine gms_zhpr

!      Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level3
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 3 Blas routine.
!*>
!*>  -- Written on 8-February-1989.
!*>     Jack Dongarra, Argonne National Laboratory.
!*>     Iain Duff, AERE Harwell.
!*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
     !*>     Sven Hammarling, Numerical Algorithms Group Ltd.
     !Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)
     !*> \endverbatim

     subroutine gms_zsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zsymm
#endif
         character(len=1),                       intent(in),value :: side
         character(len=1),                       intent(in),value :: uplo
         integer(kind=int4),                     intent(in),value :: m
         integer(kind=int4),                     intent(in),value :: n
         type(AVX512c8f64_t),                    intent(in)       :: alpha
         type(AVX512c8f64_t), dimension(lda,*),  intent(in)       :: a
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED a:64
#endif
         integer(kind=int4),                     intent(in),value :: lda
         type(AVX512c8f64_t), dimension(ldb,*),  intent(in)       :: b
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED b:64
#endif
         integer(kind=int4),                     intent(in),value :: ldb
         type(AVX512c8f64_t),                    intent(in)       :: beta
         type(AVX512c8f64_t), dimension(ldc,*),  intent(inout)    :: c
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED c:64
#endif
         integer(kind=int4),                     intent(in)       :: ldc
         ! Locals
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
         type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
         type(AVX512c8f64_t), automatic :: temp2
         integer(kind=int4),  automatic :: i,info,j,k,nrowa
         logical(kind=int4),  automatic :: upper
         logical(kind=int1),  automatic :: aeq0,beq1,beq0
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: ONE
#endif
         type(AVX512c8f64_t), parameter :: ONE  = AVX512c8f64_t([1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                                                                1.0_dp,1.0_dp,1.0_dp,1.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
         type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
         ! EXec code ...
         !  Set NROWA as the number of rows of A.
         if(lsame(side,'L')) then
            nrowa = m
         else
            nrowa = n
         end if
         upper = lsame(uplo,'U')
         !   Test the input parameters.
         info = 0
         if((.not.lsame(side,'L')) .and. (.not.lsame(side,'R'))) then
            info = 1
         else if((.not.upper) .and. (.not.lsame(uplo,'L'))) then
            info = 2
         else if(m<0) then
            info = 3
         else if(n<0) then
            info = 4
         else if(lda < max(1,nrowa)) then
            info = 7
         else if(ldb < max(1,m)) then
            info = 9
         else if(ldc < max(1,m)) then
            info = 12
         end if
         if(info/=0) then
            call xerbla('GMS_ZSYMM',info)
            return
         end if
         aeq0 = .false.
         beq1 = .false.
         aeq0 = all(alpha==ZERO)
         beq1 = all(beta==ONE)
         ! Quick return if possible.
         if((m==0) .or. (n==0) .or. &
            ((aeq0) .and. (beq1))) return
         !  And when  alpha.eq.zero.
         if(aeq0) then
            if(beq0) then
               do j=1,n
                  !DIR$ VECTOR ALIGNED
                  !DIR$ VECTOR ALWAYS
                  !DIR$ UNROLL(6)
                  do i=1,n
                     c(i,j) = ZERO
                  end do
               end do
            else
               do j=1,n
                  !DIR$ VECTOR ALIGNED
                  !DIR$ VECTOR ALWAYS
                  !DIR$ UNROLL(8)
                  do i=1,m
                     c(i,j) = beta*c(i,j)
                  end do
               end do
            end if
            return
         end if
         !  Start the operations.
         if(lsame(side,'L')) then
            ! Form  C := alpha*A*B + beta*C.
            if(upper) then
               do j=1,n
                  do i=1,m
                     temp1 = alpha*b(i,j)
                     temp2 = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#endif
                     do k=1,i-1
                        c(k,j) = c(k,j)+temp1*a(k,i)
                        temp2 = temp2+b(k,j)*a(k,i)
                     end do
                     if(beq0) then
                        c(i,j) = temp1*a(i,i)+alpha*temp2
                     else
                        c(i,j) = beta*c(i,j)+temp1*a(i,i) + &
                             alpha*temp2
                     end if
                  end do
               end do
            else
               do j=1,n
                  do i=m,1,-1
                     temp1 = alpha*b(i,j)
                     temp2 = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                     do k=i+1,m
                        c(k,j) = c(k,j)+temp1*a(k,i)
                        temp2 = temp2+b(k,j)*a(k,i)
                     end do
                     if(beq0) then
                        c(i,j) = temp1*a(i,i)+alpha*temp2
                     else
                        c(i,j) = beta*c(i,j)+temp1*a(i,i) + &
                             alpha*temp2
                     end if
                  end do
               end do
            end if
         else
            !  Form  C := alpha*B*A + beta*C.
            do j=1,n
               temp1 = alpha*a(j,j)
               if(beq0) then
                      !DIR$ VECTOR ALIGNED
                      !DIR$ VECTOR ALWAYS
                      !DIR$ UNROLL(8)
                  do i=1,m
                     c(i,j) = temp1*b(i,j)
                  end do
               else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                  do i=1,m
                     c(i,j) = beta*c(i,j)+temp1*b(i,j)
                  end do
               end if
               do k=1,j-1
                  if(upper) then
                     temp1 = alpha*a(k,j)
                  else
                     temp1 = alpha*a(j,k)
                  end if
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do i=1,m
                      c(i,j) = c(i,j)+temp1*b(i,k)
                   end do
                end do
                do k=j+1,n
                   if(upper) then
                      temp1 = alpha*a(j,k)
                   else
                      temp1 = alpha*a(k,j)
                   end if
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do i=1,m
                      c(i,j) = c(i,j)+temp1*b(i,k)
                   end do
                end do
             end do
          end if
          !End of ZSYMM
     end subroutine gms_zsymm

!*  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!1*> \ingroup complex16_blas_level3
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 3 Blas routine.
!*>
!*>  -- Written on 8-February-1989.
!*>     Jack Dongarra, Argonne National Laboratory.
!*>     Iain Duff, AERE Harwell.
!*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
     !*>     Sven Hammarling, Numerical Algorithms Group Ltd.
      !Modified by Bernard Gingold on 29-11-2019 (removing built-in complex*16 data type,using modern Fortran features)
!*> \endverbatim     
     subroutine gms_zsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zsyr2k
#endif
        character(len=1),                      intent(in),value :: uplo
        character(len=1),                      intent(in),value :: trans
        integer(kind=int4),                    intent(in),value :: n
        integer(kind=int4),                    intent(in),value :: k
        type(AVX512c8f64_t),                   intent(in)       :: alpha
        type(AVX512c8f64_t), dimension(lda,*), intent(in)       :: a
#if defined __INTEL_COMPILER
        !DIR$ ASSUMED_ALIGNED a:64
#endif
        integer(kind=int4),                    intent(in),value :: lda
        type(AVX512c8f64_t), dimension(ldb,*), intent(in)       :: b
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED b:64
#endif
        integer(kind=int4),                    intent(in),value :: ldb
        type(AVX512c8f64_t),                   intent(in)       :: beta
        type(AVX512c8f64_t), dimension(ldc,*), intent(inout)    :: c
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED c:64
#endif
        integer(kind=int4),                    intent(in)       :: ldc
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
        type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
        type(AVX512c8f64_t), automatic :: temp2
        integer(kind=int4),  automatic :: i,info,j,l,nrowa
        logical(kind=int4),  automatic :: upper
        logical(kind=int1),  automatic :: aeq0,beq1,beq0,bneq1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ONE
#endif
        type(AVX512c8f64_t), parameter :: ONE  = AVX512c8f64_t([1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                                                                1.0_dp,1.0_dp,1.0_dp,1.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        ! Test the input parameters.
        if(lsame(trans,'N')) then
           nrowa = n
        else
           nrowa = k
        end if
        upper = lsame(uplo,'U')
        info  = 0
        if((.not.upper) .and. (.not.lsame(uplo,'L'))) then
           info = 1
        else if((.not.lsame(trans,'N')) .and. &
             (.not.lsame(trans,'T'))) then
           info = 2
        else if(n<0) then
           info = 3
        else if(k<0) then
           info = 4
        else if(lda < max(1,nrowa)) then
           info = 7
        else if(ldb < max(1,nrowa)) then
           info = 9
        else if(ldc < max(1,n)) then
           info = 12
        end if
        if(info/=0) then
           call xerbla('GMS_ZSYRK2',info)
           return
        end if
        !   Quick return if possible.
        aeq0 = .false.
        beq1 = .false.
        aeq0 = all(alpha==ZERO)
        beq1 = all(beta==ONE)
        if((n==0) .or. (((aeq0) .or. &
             (k==0)) .and. (beq1))) return
        !  And when  alpha.eq.zero.
        if(aeq0) then
           if(upper) then
              if(beq0) then
                 do j=1,n
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ UNROLL(6)
                    do i=1,j
                       c(i,j) = ZERO
                    end do
                 end do
              else
                 do j=1,n
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
                    do i=1,j
                       c(i,j) = beta*c(i,j)
                    end do
                 end do
              end if
           else
              if(beq0) then
                 do j=1,n
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ UNROLL(6)
                    do i=j,n
                       c(i,j) = ZERO
                    end do
                 end do
              else
                 do j=1,n
#if defined __INTEL_COMPILER                    
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#endif
                    do i=j,n
                       c(i,j) = beta*c(i,j)
                    end do
                 end do
              end if
           end if
           return
        end if
        !  Start the operations.
        if(upper) then
           do j=1,n
              if(beq0) then
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ UNROLL(6)
                 do i=1,j
                    c(i,j) = ZERO
                 end do
              else if(bneq1) then
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif 
                 do i=1,j
                    c(i,j) = beta*c(i,j)
                 end do
              end if
              do l=1,k
                 if(all(a(j,l)/=ZERO) .or. all(b(j,l)/=ZERO)) then
                    temp1 = alpha*b(j,l)
                    temp2 = alpha*a(j,l)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                 
                    do i=1,j
                      
                       c(i,j) = c(i,j)+a(i,l)*temp1 +
                       b(i,l)*temp2
                    end do
                 end if
              end do
           end do
        else
           do j=1,n
              if(beq0) then
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ UNROLL(6)
                 do i=j,n
                    c(i,j) = ZERO
                 end do
              else if(bneq1) then
                  !DIR$ VECTOR ALIGNED
                  !DIR$ VECTOR ALWAYS
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif  
                 do i=j,n
                    c(i,j) = beta*c(i,j)
                 end do
              end if
              do l=1,k
                 if(all(a(j,l)/=ZERO) .or. all(b(j,l)/=ZERO)) then
                    temp1 = alpha*b(j,l)
                    temp2 = alpha*a(j,l)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                   
                    do i=j,n
                       c(i,j) = c(i,j)+a(i,l)*temp1 + &
                            b(i,l)*temp2
                    end do
                 end if
              end do
           end do
        end if
     else
        !  Form  C := alpha*A**T*B + alpha*B**T*A + C.
        if(upper) then
           do j=1,n
              do i=1,j
                 temp1 = ZERO
                 temp2 = ZERO
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                 do l=1,k
                    temp1 = temp1+a(l,i)*b(l,j)
                    temp2 = temp2+b(l,i)*a(l,j)
                 end do
                 if(beq0) then
                    c(i,j) = alpha*temp1+alpha*temp2
                 else
                    c(i,j) = beta*c(i,j)+alpha*temp1 + &
                         alpha*temp2
                 end if
              end do
           end do
        else
           do j=1,n
              do i=j,n
                 temp1 = ZERO
                 temp = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                 do l=1,k
                    temp1 = temp1+a(l,i)*b(l,j)
                    temp2 = temp2+b(l,i)*a(l,j)
                 end do
                 if(beq0) then
                    c(i,j) = alpha*temp1+alpha*temp2
                 else
                    c(i,j) = beta*c(i,j)+alpha*temp1 + &
                         alpha*temp2
                 end if
              end do
           end do
        end if
     end if
      ! End of ZSYR2K
   end subroutine gms_zsyr2k

!     Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level3
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 3 Blas routine.
!*>
!*>  -- Written on 8-February-1989.
!*>     Jack Dongarra, Argonne National Laboratory.
!*>     Iain Duff, AERE Harwell.
!*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
   !*>     Sven Hammarling, Numerical Algorithms Group Ltd.
   ! !Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)
     !*> \endverbatim

     subroutine gms_zsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
#if defined __INTEL_COMPILER     
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zsyrk
#endif
        character(len=1),                      intent(in),value  :: uplo
        character(len=1),                      intent(in),value  :: trans
        integer(kind=int4),                    intent(in),value  :: n
        integer(kind=int4),                    intent(in),value  :: k
        type(AVX512c8f64_t),                   intent(in)        :: alpha
        type(AVX512c8f64_t), dimension(lda,*), intent(in)        :: a
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED a:64
#endif
        integer(kind=int4),                    intent(in),value  :: lda
        type(AVX512c8f64_t),                   intent(in)        :: beta
        type(AVX512c8f64_t), dimension(ldc,*), intent(inout)     :: c
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED c:64
#endif
        integer(kind=int4),                    intent(in),value  :: ldc
        ! LOcals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
        integer(kind=int4),  automatic :: i,info,j,l,nrowa
        logical(kind=int4),  automatic :: upper
        logical(kind=int1),  automatic :: aeq0,beq1,beq0,bneq1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ONE
#endif
        type(AVX512c8f64_t), parameter :: ONE  = AVX512c8f64_t([1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                                                                1.0_dp,1.0_dp,1.0_dp,1.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        ! Test the input parameters.
        if(lsame(trans,'N')) then
           nrowa = n
        else
           nrowa = k
        end if
        upper = lsame(uplo,'U')
        info = 0
        if((.not.upper) .and. (.not.lsame(uplo,'L'))) then
           info = 1
        else if((.not.lsame(trans,'N')) .and. &
             (.not.lsame(trans,'T'))) then
           info = 2
        else if(n<0) then
           info = 3
        else if(k<0) then
           info = 4
        else if(lda < max(1,nrowa)) then
           info = 7
        else if(ldc < max(1,n)) then
           info = 10
        end if
        if(info/=0) then
           call xerbla('GMS_ZSYRK',info)
           return
        end do
        ! Quick return if possible.
        aeq0 = .false.
        beq1 = .false.
        aeq0 = all(alpha==ZERO)
        beq1 = all(beta==ONE)
        if((n==0) .or. (((aeq0) .or. &
             (k==0)) .and.(beq1))) return
        !  And when  alpha.eq.zero.
        if(aeq0) then
           if(upper) then
              if(beq0) then
                 do j=1,n
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ UNROLL(6)
                    do i=1,j
                       c(i,j) = ZERO
                    end do
                 end do
              else
                 do j=1,n
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#endif
                    do i=1,j
                       c(i,j) = beta*c(i,j)
                    end do
                 end do
              end if
           else
              if(beq0) then
                 do j=1,n
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ UNROLL(6)
                    do i=j,n
                       c(i,j) = ZERO
                    end do
                 end do
              else
                 do j=1,n
#if defined __INTEL_COMPILER                    
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#endif
                    do i=j,n
                       c(i,j) = beta*c(i,j)
                    end do
                 end do
              end if
           end if
           return
        end if
        !  Start the operations.
        bneq1 = .false.
        beq0  = .false.
        bneq1 = all(beta/=ONE)
        beq0  = all(beta==ZERO)
        if(lsame(trans,'N')) then
           !  Form  C := alpha*A*A**T + beta*C.
           if(upper) then
              do j=1,n
                 if(beq0) then
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ UNROLL(6)
                    do i=1,j
                       c(i,j) = ZERO
                    end do
                 else if(bneq1) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#endif
                    do i=1,j
                       c(i,j) = beta*c(i,j)
                    end do
                 end if
                 do l=1,k
                    if(all(a(j,l)/=ZERO)) then
                       temp = alpha*a(j,l)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=1,j
                          c(i,j) = c(i,j)+temp*a(i,l)
                       end do
                    end if
                 end do
              end do
           else
              do j=1,n
                 if(beq0) then
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ UNROLL(6)
                    do i=j,n
                       c(i,j) = ZERO
                    end do
                 else if(bneq1) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(4)
#endif
#endif
                    do i=j,n  
                      c(i,j) = beta*c(i,j)
                   end do
                end if
                do l=1,k
                   if(all(a(j,l)/=ZERO)) then
                      temp = alpha*a(j,l)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                      do i=j,n
                         c(i,j) = c(i,j)+temp*a(i,l)
                      end do
                   end if
                end do
             end do
          end if
       else
          !   Form  C := alpha*A**T*A + beta*C.
          if(upper) then
             do j=1,n
                do i=1,j
                   temp = ZERO
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do l=1,k
                      temp = temp+a(l,i)*a(l,j)
                   end do
                   if(beq0) then
                      c(i,j) = alpha*temp
                   else
                      c(i,j) = alpha*temp+beta*c(i,j)
                   end if
                end do
             end do
          else
             do j=1,n
                do i=1,n
                   temp = ZERO
#if defined __INTE_COMPILER
                     !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                   do l=1,k
                      temp = temp+a(l,i)*a(l,j)
                   end do
                   if(beq0) then
                      c(i,j) = alpha*temp
                   else
                      c(i,j) = alpha*temp+beta*c(i,j)
                   end if
                end do
             end do
          end if
       end if
       ! End of ZSYRK
     end subroutine gms_zsyrk

!  Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
     !*>     Richard Hanson, Sandia National Labs.
      !Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)
     !*> \endverbatim

     subroutine gms_ztbsv(uplo,trans,diag,n,k,a,lda,x,incx)
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_ztbsv
#endif
        character(len=1),                         intent(in),value  :: uplo
        character(len=1),                         intent(in),value  :: trans
        character(len=1),                         intent(in),value  :: diag
        integer(kind=int4),                       intent(in),value  :: n
        integer(kind=int4),                       intent(in),value  :: k
        type(AVX512c8f64_t), dimension(lda,*),    intent(in)        :: a
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED a:64
#endif
        integer(kind=int4),                       intent(in),value  :: lda
        type(AVX512c8f64_t), dimension(*),        intent(inout)     :: x
#if defined __INTEL_COMPILER
        !DIR$ ASSUME_ALIGNED x:64
#endif
        integer(kind=int4),                       intent(in),value  :: incx
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
        integer(kind=int4),  automatic :: i,info,ix,j,jx,kplus1,kx,l
        logical(kind=int4),  automatic :: noconj,nounit
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        ! Test the input parameters.
        info = 0
        if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
           info = 1
        else if(.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. &
             .not.lsame(trans,'C')) then
           info = 2
        else if(.not.lsame(diag,'U') .and. .not.lsame(diag,'N')) then
           info = 3
        else if(n<0) then
           info = 4
        else if(k<0) then
           info = 5
        else if(lda<(k+1)) then
           info = 7
        else if(incx==0) then
           info = 9
        end if
        if(info/=0) then
           call xerbla('GMS_ZTBSV',info)
           return
        end if
        !  Quick return if possible.
        if(n==0) return
        noconj = lsame(trans,'T')
        nounit = lsame(diag,'N')
        !   Set up the start point in X if the increment is not unity. This
        !   *     will be  ( N - 1 )*INCX  too small for descending loops.
        if(incx<0) then
           kx = 1-(n-1)*incx
        else if(incx/=1) then
           kx = 1
        end if
        !   Start the operations. In this version the elements of A are
        !*     accessed by sequentially with one pass through A.
        if(lsame(trans,'N')) then
           !   Form  x := inv( A )*x.
           if(lsame(uplo,'U')) then
              kplus1 = k+1
              if(incx==1) then
                 do j=n,1,-1
                    if(all(x(j)/=ZERO)) then
                       l = kplus1-j
                       if(nounit) x(j)/a(kplus1,j)
                       temp = x(j)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=j-1,max(1,j-k),-1
                          x(i) = x(i)-temp*a(l+i,j)
                       end do
                    end if
                 end do
              else
                 kx = kx+(n-1)*incx
                 jx = kx
                 do j=n,1,-1
                    kx = kx-incx
                    if(all(x(jx)/=ZERO)) then
                       ix = kx
                       l = kplus1-j
                       if(nounit) x(jx) = x(jx)/a(kplus1,j)
                       temp = x(jx)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=j-1,max(1,j-k),-1
                          x(ix) = x(ix)-temp*a(l+i,j)
                       end do
                    end if
                    jx = kx-incx
                 end do
              end if
           else
              if(incx==1) then
                 do j=1,n
                    if(all(x(j)/=ZERO)) then
                       l = 1-j
                       if(nounit) x(j) = x(j)/a(1,j)
                       temp = x(j)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=j+1,min(n,j+k)
                          x(i) = x(i)-temp*a(l+i,j)
                       end do
                    end if
                 end do
              else
                 jx = kx
                 do j=1,n
                    kx = kx+incx
                    if(all(x(jx)/=ZERO)) then
                       ix = kx
                       l = 1-j
                       if(nounit) x(jx) = x(jx)/a(1,j)
                       temp = x(jx)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                  
                       do i=j+1,min(n,j+k)
                          x(ix) = x(ix)-temp*a(l+i,j)
                          ix = ix+incx
                       end do
                    end if
                    jx = jx+incx
                 end do
              end if
           end if
        else
           !  Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
           if(lsame(uplo,'U')) then
              kplus1 = k+1
              if(incx==1) then
                 do j=1,n
                    temp = x(j)
                    l = kplus1-j
                    if(noconj) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                  
                       do i=max(1,j-k),j-1
                          temp = temp-a(l+i,j)*x(i)
                       end do
                       if(nounit) temp = temp/a(kplus1,j)
                    else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                 
                       do i=max(1,j-k),j-1
                          temp = temp-conjugate(a(l+i,j))*x(i)
                       end do
                       if(nounit) temp = temp/conjugate(a(kplus1,j))
                    end if
                    x(j) = temp
                 end do
              else
                 jx = kx
                 do j=1,n
                    temp = x(jx)
                    ix = kx
                    l = kplus1-j
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                 
                       do i=max(1,j-k),j-1
                          temp = temp-a(l+i,j)*x(ix)
                          ix = ix+incx
                       end do
                       if(nounit) temp = temp/a(kplus1,j)
                    else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                 
                       do i=max(1,j-k),j-1
                          temp = temp-conjugate(a(l+i,j))*x(ix)
                          ix = ix+incx
                       end do
                       if(nounit) temp = temp/conjugate(a(kplus1,j))
                    end if
                    x(jx) = temp
                    jx = jx+incx
                    if(j>k) kx = kx+incx
                 end do
              end if
           else
              if(incx==1) then
                 do j=n,1,-1
                    temp = x(j)
                    l =1-j
                    if(noconj) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                   
                       do i=min(n,j+k),j+1,-1
                          temp = temp-a(l+i,j)*x(i)
                       end do
                       if(nounit) temp = temp/a(1,j)
                    else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=min(n,j+k),j+1,-1
                          temp = temp-conjugate(a(l+i,j))*x(i)
                       end do
                       if(nounit) temp = temp/conjugate(a(1,j))
                    end if
                    x(j) = temp
                 end do
              else
                 kx = kx+(n-1)*incx
                 jx = kx
                 do j=n,1,-1
                    temp = x(jx)
                    ix = kx
                    l = 1-j
                    if(noconj) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                       do i=min(n,j+k),j+1,-1
                          temp = temp-a(l+i,j)*x(i)
                          ix = ix+incx
                       end do
                       if(nounit) temp = temp/a(1,j)
                    else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                      
                       do i=min(n,j+k),j+1,-1
                          temp = temp-conjugate(a(l+i,j))*x(ix)
                          ix = ix+incx
                       end do
                       if(nounit) temp = temp/conjugate(a(1,j))
                    end if
                    x(jx) = temp
                    jx = jx-incx
                    if((n-j)>=k) kx = kx-incx
                 end do
              end if
           end if
        end if
        !End of ZTBSV
     end subroutine gms_ztbsv

! Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date December 2016
!*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!*>
!*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
     !*>     Richard Hanson, Sandia National Labs.
     !Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)
     !*> \endverbatim

     subroutine gms_ztpmv(uplo,trans,diag,n,ap,x,incx)
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_ztpmv
#endif
         character(len=1),                         intent(in),value :: uplo
         character(len=1),                         intent(in),value :: trans
         character(len=1),                         intent(in),value :: diag
         integer(kind=int4),                       intent(in),value :: n
         type(AVX512c8f64_t), dimension(*),        intent(in)       :: ap
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED ap:64
#endif
         type(AVX512c8f64_t), dimension(*),        intent(inout)    :: x
#if defined __INTEL_COMPILER
         !DIR$ ASSUME_ALIGNED x:64
#endif
         integer(kind=int4),                       intent(in),value :: incx
         ! LOcals
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
         type(AVX512c8f64_t), automatic :: temp
         integer(kind=int4),  automatic :: i,info,ix,j,jx,k,kk,kx
         logical(kind=int4),  automatic :: noconj,nounit
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        ! Exec code ....
        info = 0
        if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
           info = 1
        else if(.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. &
             .not.lsame(trans,'C')) then
           info = 2
        else if(.not.lsame(diag,'U') .and. .not.lsame(diag,'N')) then
           info = 3
        else if(n<0) then
           info = 4
        else if(incx==0) then
           info = 7
        end if
        if(info/=0) then
           call xerbla('GMS_ZTPMV',info)
           return
        end if
        !  Quick return if possible.
        if(n==0) return
        noconj = lsame(trans,'T')
        nounit = lsame(diag,'N')
        ! Set up the start point in X if the increment is not unity. This
        !*     will be  ( N - 1 )*INCX  too small for descending loops.
        if(incx<=0) then
           kx = 1-(n-1)*incx
        else if(incx/=1) then
           kx = 1
        end if
        !   Start the operations. In this version the elements of AP are
        !*     accessed sequentially with one pass through AP.
        if(lsame(trans,'N')) then
           !    Form  x:= A*x.
           if(lsame(uplo,'U')) then
              kk = 1
              if(incx==1) then
                 do j=1,n
                    if(all(x(j)/=ZERO)) then
                       temp = x(j)
                       k = kk
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                        
                       do i=1,j-1
                          x(i) = x(i)+temp*ap(k)
                          k = k+1
                       end do
                       if(nounit) x(j) = x(j)*ap(kk+j-1)
                    end if
                    kk = kk+j
                 end do
              else
                 jx = kx
                 do j=1,n
                    if(all(x(jx)/=ZERO)) then
                       temp = x(jx)
                       ix = kx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                        
                       do k=kk,kk+j-2
                          x(ix) = x(ix)+temp*ap(k)
                          ix = ix+incx
                       end do
                       if(nounit) x(jx) = x(jx)*ap(kk+j-1)
                    end if
                    jx = jx+incx
                    kk = kk+j
                 end do
              end if
           else
              kk = (n*(n+1))/2
              if(incx==1) then
                 do j=n,1,-1
                    if(all(x(j)/=ZERO)) then
                       temp = x(j)
                       k = kk
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif  
                       do i=n,j+1,-1
                          x(i) = x(i)+temp*ap(k)
                          k = k-1
                       end do
                       if(nounity) x(j) = x(j)*ap(kk-n+j)
                    end if
                    kk = kk-(n-j+1)
                 end do
              else
                 kx = kx+(n-1)*incx
                 jx = kx
                 do j=n,1,-1
                    if(all(x(jx)/=ZERO)) then
                       temp = x(jx)
                       ix = kx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif   
                       do k=kk,kk-(n-(j+1)),-1
                          x(ix) = x(ix)+temp*ap(k)
                          ix = ix+incx
                       end do
                       if(nounit) x(jx) = x(jx)*ap(kk-n+j)
                    end if
                    jx = jx-incx
                    kk = kk-(n-j+1)
                 end do
              end if
            end if
         else
            !   Form  x := A**T*x  or  x := A**H*x.
            if(lsame(uplo,'U')) then
               kk = (n*(n+1))/2
               if(incx==1) then
                  do j=n,1,-1
                     temp = x(j)
                     k = kk-1
                     if(noconj) then
                        if(nounit) temp = temp*ap(kk)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                        
                        do i=j-1,1,-1
                           temp = temp+ap(k)*x(i)
                           k = k-1
                        end do
                     else
                        if(nounit) temp = temp*conjugate(ap(kk))
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                        
                        do i=j-1,1,-1
                           temp = temp+conjugate(ap(k))*x(i)
                           k = k-1
                        end do
                     end if
                     x(j) = temp
                     kk = kk-j
                  end do
               else
                  jx = kx+(n-1)*incx
                  do j=n,1,-1
                     temp = x(jx)
                     ix = jx
                     if(noconj) then
                        if(nounit) temp = temp*ap(kk)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                         
                        do k=kk-1,kk-j+1,-1
                           ix = ix-incx
                           temp = temp+ap(k)*x(ix)
                        end do
                     else
                        if(nounit) temp = temp*conjugate(ap(kk))
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                        
                        do k=kk-1,kk-j+1,-1
                           ix = ix-incx
                           temp = temp+conjugate(ap(k))*x(ix)
                        end do
                     end if
                     x(jx) = temp
                     jx = jx-incx
                     kk = kk-j
                  end do
               end if
            else
               kk = 1
               if(incx==1) then
                  do j=1,n
                     temp = x(j)
                     k = kk+1
                     if(noconj) then
                        if(nounit) temp = temp*ap(kk)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                         
                        do i=j+1,n
                           temp = temp+ap(k)*x(i)
                           k = k+1
                        end do
                     else
                        if(nounit) temp = temp+conjugate(ap(kk))
                        do i=j+1,n
                           temp = temp+conjugate(ap(k))*x(i)
                           k = k+1
                        end do
                     end if
                     x(j) = temp
                     kk = kk+(n-j+1)
                  end do
               else
                  jx = kx
                  do j=1,n
                     temp = x(jx)
                     ix = jx
                     if(noconj) then
                        if(nounit) temp = temp*ap(kk)
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                         
                        do k=kk+1,kk+n-j
                           ix = ix+incx
                           temp = temp+ap(k)*x(ix)
                        end do
                     else
                        if(nounit) temp = temp*conjugate(ap(k))
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                         
                       do k=kk+1,kk+n-j
                          ix = ix+incx
                          temp = temp+conjugate(ap(k))*x(ix)
                       end do
                    end if
                    x(jx) = temp
                    jx = jx+incx
                    kk = kk+(n-j+1)
                 end do
              end if
            end if
         end if
         ! End of ZTPMV
      end subroutine gms_ztpmv

! Authors:
!*  ========
!*
!*> \author Univ. of Tennessee
!1*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!1*> \date December 2016
!*
!*> \ingroup complex16_blas_level2
!*
!*> \par Further Details:
!*  =====================
!*>
!*> \verbatim
!*>
!*>  Level 2 Blas routine.
!*>
!1*>  -- Written on 22-October-1986.
!*>     Jack Dongarra, Argonne National Lab.
!*>     Jeremy Du Croz, Nag Central Office.
!*>     Sven Hammarling, Nag Central Office.
      !*>     Richard Hanson, Sandia National Labs.
      !Modified by Bernard Gingold on 29-11-2019 (removing build-in complex*16 data type,using modern Fortran features)
!*> \endverbatim
      !*>

      subroutine gms_ztpsv(uplo,trans,diag,n,ap,x,incx)
#if defined __INTEL_COMPILER
           !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_ztpsv
#endif
           character(len=1),                    intent(in),value :: uplo
           character(len=1),                    intent(in),value :: trans
           character(len=1),                    intent(in),value :: diag
           integer(kind=int4),                  intent(in),value :: n
           type(AVX512c8f64_t), dimension(*),   intent(in)       :: ap
#if defined __INTEL_COMPILER
           !DIR$ ASSUME_ALIGNED ap:64
#endif
           type(AVX512c8f64_t), dimension(*),   intent(inout)    :: x
#if defined __INTEL_COMPILER
           !DIR$ ASSUME_ALIGNED x:64
#endif
           integer(kind=int4),                  intent(in),value :: incx
           ! LOcals
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
           type(AVX512c8f64_t), automatic :: temp
           integer(kind=int4),  automatic :: i,info,ix,j,jx,k,kk,kx
           logical(kind=int4),  automatic :: noconj,nounit           
         
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
           type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
           ! Exec code ....
           info = 0
           if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
              info = 1
           else if(.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. &
              .not.lsame(trans,'C')) then
              info = 2
           else if(.not.lsame(diag,'U') .and. .not.lsame(diag,'N')) then
              info = 3
           else if(n<0) then
              info = 4
           else if(incx==0) then
              info = 7
           end if
           if(info/=0) then
              call xerbla('GMS_ZTPSV',info)
              return
           end if
           ! Quick return if possible.
           if(n==0) return
           noconj = lsame(trans,'T')
           nounit = lsame(diag,'N')
           !  Set up the start point in X if the increment is not unity. This
           !  *     will be  ( N - 1 )*INCX  too small for descending loops.
           if(incx<=0) then
              kx = 1-(n-1)*incx
           else if(incx/=1) then
              kx = 1
           end if
           !  Start the operations. In this version the elements of AP are
           !*     accessed sequentially with one pass through AP.
           if(lsame(trans,'N')) then
              !  Form  x := inv( A )*x.
              if(lsame(uplo,'U')) then
                 kk = (n*(n+1))/2
                 if(incx==1) then
                    do j=n,1,-1
                       if(all(x(j)/=ZERO)) then
                          if(nounit) x(j) = x(j)/ap(kk)
                          temp = x(j)
                          k = kk-1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                          
                          do i=j-1,1,-1
                             x(i) = x(i)-temp*ap(k)
                             k = k-1
                          end do
                       end if
                       kk = kk-j
                    end do
                 else
                    jx = kx+(n-1)*incx
                    do j=n,1,-1
                       if(all(x(jx)/=ZERO)) then
                          if(nounit) x(jx) = x(jx)/ap(kk)
                          temp = x(jx)
                          ix = jx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                          
                          do k=kk-1,kk-j+1,-1
                             ix = ix-incx
                             x(ix) = x(ix)-temp*ap(k)
                          end do
                       end if
                       jx = jx-incx
                       kk = kk-j
                    end do
                 end if
              else
                 kk = 1
                 if(incx==1) then
                    do j=1,n
                       if(all(x(j)/=ZERO)) then
                          if(nounit) x(j) = x(j)/ap(kk)
                          temp = x(j)
                          k = kk+1
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                          
                          do i=j+1,n
                             x(i) = x(i)-temp*ap(k)
                             k = k+1
                          end do
                       end if
                       kk = kk+(n-j+1)
                    end do
                 else
                    jx = kx
                    do j=1,n
                       if(all((jx)/=ZERO)) then
                          if(nounit) x(jx) = x(jx)/ap(kk)
                          temp = x(jx)
                          ix = jx
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                           
                          do k=kk+1,kk+n-j
                             ix = ix+incx
                             x(ix) = x(ix)-temp*ap(k)
                          end do
                       end if
                       jx = jx+incx
                       kk = kk+(n-j+1)
                    end do
                 end if
              end if
           else
              !  Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
              if(lsame(uplo,'U')) then
                 kk = 1
                 if(incx==1) then
                    do j=1,n
                       temp = x(j)
                       k = kk
                       if(noconj) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                          
                          do i=1,j-1
                             temp = temp-ap(k)*x(i)
                             k = k+1
                          end do
                          if(nounit) temp = temp/ap(kk+j-1)
                       else
                          do i=1,j-1
                             temp = temp-conjugate(ap(k))*x(i)
                             k = k+1
                          end do
                          if(nounit) temp = temp/conjugate(ap(kk+j-1))
                       end if
                       x(j) = temp
                       kk = kk+j
                    end do
                 else
                    jx = kx
                    do j=1,n
                       temp = x(jx)
                       ix = kx
                       if(noconj) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                          
                          do k=kk,kk+j-2
                             temp = temp-ap(k)*x(ix)
                             ix = ix+incx
                          end do
                          if(nounit) temp = temp/ap(kk+j-1)
                       else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                          
                          do k=kk,kk+j-2
                             temp = temp-conjugate(ap(k))*x(ix)
                             ix = ix+incx
                          end do
                          if(nounit) temp = temp/conjugate(ap(kk+j-1))
                       end if
                       x(jx) = temp
                       jx = jx+incx
                       kk = kk+j
                    end do
                 end if
              else
                 kk = (n*(n+1))/2
                 if(incx==1) then
                    do j=n,1,-1
                       temp = x(j)
                       k = kk
                       if(noconj) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                           
                          do i=n,j+1,-1
                             temp = temp-ap(k)*x(i0
                             k = k-1
                          end do
                          if(nounit) temp = temp/ap(kk-n+j)
                       else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                          
                          do i=n,j+1,-1
                             temp = temp/conjugate(ap(kk-n+j))
                             k = k-1
                          end do
                          if(nounit) temp = temp/conjugate(ap(kk-n+j))
                       end if
                       x(j) = temp
                       kk = kk-(n-j+1)
                    end do
                 else
                    kx = kx+(n-1)*incx
                    jx = kx
                    do j=n,1,-1
                       temp = x(jx)
                       ix = kx
                       if(noconj) then
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif
                          
                          do k=kk,kk-(n-(j+1)),-1
                             temp = temp-ap(k)*x(ix)
                             ix = ix-incx
                          end do
                          if(nounit) temp = temp/ap(kk-n+j)
                       else
#if defined __INTEL_COMPILER
                    !DIR$ VECTOR ALIGNED
                    !DIR$ VECTOR ALWAYS
                    !DIR$ FMA
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !DIR$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !DIR$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !DIR$ UNROLL(5)
#endif
#elif defined __GFORTRAN__
#if (CURRENT_PROCESSOR_ARCH_NAME) == 1
                    !GCC$ UNROLL(10)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 2
                    !GCC$ UNROLL(8)
#elif (CURRENT_PROCESSOR_ARCH_NAME) == 3
                    !GCC$ UNROLL(4)
#endif
#endif                           
                          do k=kk,kk-(n-(j+1)),-1
                             temp = temp-conjugate(ap(k))*x(ix)
                             ix = ix-incx
                          end do
                          if(nounit) temp = temp/conjugate(ap(kk-n+j))
                       end if
                       x(jx) = temp
                       jx = jx-incx
                       kk = kk-(n-j+1)
                    end do
                 end if
              end if
           end if
           ! End of ZTPSV
      end subroutine gms_ztpsv
        
    
end module mod_blas
