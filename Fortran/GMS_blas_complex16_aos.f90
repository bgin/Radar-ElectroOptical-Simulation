






 
    !use module_kinds, only : i4,dp
    !use mod_avx512c8f64
    !implicit none

     !=====================================================59
     !  File and module information:
     !  version,creation and build date, author,description
     !=====================================================59

    ! Major version
    !integer(kind=i4), parameter, public :: MOD_BLAS_MAJOR = 1
    ! MInor version
    !integer(kind=i4), parameter, public :: MOD_BLAS_MINOR = 0
    ! Micro version
    !integer(kind=i4), parameter, public :: MOD_BLAS_MICRO = 0
    ! Module full version
    !integer(kind=i4), parameter, public :: MOD_BLAS_FULLVER = &
    !     1000*MOD_BLAS_MAJOR+100*MOD_BLAS_MINOR+10*MOD_BLAS_MICRO
    !Module creation date
    !character(*),       parameter, public :: MOD_BLAS_CREATION_DATE = "29-11-2019 10:55 +00200 (FRI 29 NOV 2019 GMT+2)"
    ! Module build date
    !character(*),       parameter, public :: MOD_BLAS_BUILD_DATE    = __DATE__ " " __TIME__
    ! Module author info
   ! character(*)        parameter, public :: MOD_BLAS_AUTHOR = "LAPACK original authors[all rights reserved] -- This version was  modified by Bernard Gingold, contact: beniekg@gmail.com"
    ! Module short description
    !character(*)        parameter, public :: MOD_BLAS_SYNOPSIS = "Explicitly vectorized complex  blas implementation, which is based on packed AoS complex data type (AVX512c8f64_t) "
                                                                  
    ! public

 

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
subroutine zaxpy(n,za,zx,incx,zy,incy) !GCC$ ATTRIBUTES inline :: zaxpy !GCC$ ATTRIBUTES aligned(32) :: zaxpy
#if defined(__ICC) || defined(__INTEL_COMPILER)
  subroutine zaxpy(n,za,zx,incx,zy,incy)
      !DIR$ ATTRIBUTES FORCEINLINE :: gms_zaxpy
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zaxpy
       !DIR$ OPTIMIZE : 3
    !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zaxpy
#endif
      use module_kinds, only : i4
      use mod_avx512c8f64
      use omp_lib
      use mod_vecconsts, only : v8_n0
      implicit none
      integer(kind=i4),                             intent(in),value    :: n
      type(AVX512c8f64_t),                            intent(in)          :: za
      type(AVX512c8f64_t), dimension(:), allocatable, intent(in)          :: zx
      integer(kind=i4),                intent(in),value    :: incx
      type(AVX512c8f64_t), dimension(*), intent(inout)       :: zy
      integer(kind=i4),                intent(in),value    :: incy
      ! Locals
      integer(kind=i4), automatic :: i,ix,iy
      ! EXec code .....
      if(n<=0) return
      if(all(cabs_zmm8c8(za) == v8_n0.v)) return
      if(incx==1 .and. incy==1) then
         ! *        code for both increments equal to 1
         !$OMP SIMD ALIGNED(zy:64,zy) LINEAR(i:1) UNROLL PARTIAL(4)
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
          !$OMP SIMD ALIGNED(zy:64,zy) LINEAR(i:1) UNROLL PARTIAL(4)
         do i = 1,n
            zy(iy) = zy(iy)+za*zx(ix)
            ix = ix+incx
            iy = iy+incy
         end do
      end if
    end subroutine zaxpy

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
    subroutine zcopy(n,zx,incx,zy,incy) !GCC$ ATTRIBUTES inline :: zcopy !GCC$ ATTRIBUTES aligned(32) :: zcopy
#if defined(__INTEL_COMPILER) || defined(__ICC)
     subroutine zcopy(n,zx,incx,zy,incy)  
        !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: gms_zcopy
       !DIR$ ATTRIBUTES FORCEINLINE :: gms_copy
        !DIR$ OPTIMIZE : 3
        !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: gms_zcopy
#endif
        use module_kinds, only : i4
      use mod_avx512c8f64
      use omp_lib
      implicit none
      integer(kind=i4),                                intent(in),value  :: n
      type(AVX512c8f64_t), dimension(:), allocatable,  intent(in)        :: zx
      integer(kind=i4),                                intent(in),value  :: incx
      type(AVX512c8f64_t), dimension(:), allocatable,  intent(out)       :: zy
      integer(kind=i4),                                intent(in),value  :: incy
      ! LOcals
      integer(kind=i4), automatic :: i,ix,iy
      ! EXec code ...
      if(n<=0) return
      if(incx==1 .and. incy==1) then
         
         !  code for both increments equal to 1
      
         !$OMP SIMD ALIGNED(zy:64,zx) LINEAR(I:1) UNROLL PARTIAL(4)
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
           !$OMP SIMD ALIGNED(zy:64,zx) LINEAR(I:1) UNROLL PARTIAL(4)
         do i = 1,n
            zy(iy) = zx(ix)
            ix = ix+incx
            iy = iy+incy
         end do
      end if
    end subroutine zcopy

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))    
    function zdotc(n,zx,incx,zy,incy) result(dotc) !GCC$ ATTRIBUTES inline :: zdotc !GCC$ ATTRIBUTES aligned(32) :: zdotc
#if defined __INTEL_COMPILER
    function zdotc(n,zx,incx,zy,incy) result(dotc) 
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zdotc
       !DIR$ ATTRIBUTES FORCEINLINE :: zdotc
        !DIR$ OPTIMIZE : 3
        !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zdotc
#endif
         use module_kinds, only : i4
      use mod_avx512c8f64
      use omp_lib
      implicit none
      integer(kind=i4),                              intent(in),value :: n
      type(AVX512c8f64_t), dimension(:), allocatable,  intent(in)       :: zx
      integer(kind=i4),                              intent(in),value :: incx
      type(AVX512c8f64_t), dimension(:), allocatable,  intent(in)       :: zy
      integer(kind=i4),                              intent(in),value :: incy
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES ALIGN : 64 :: dotc
#endif
      type(AVX512c8f64_t) :: dotc
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES ALIGN : 64 :: ztemp
#endif
      type(AVX512c8f64_t), automatic :: ztemp
      integer(kind=i4), automatic :: i,ix,iy
      ! EXec code ....
      !
      if(n<=0) return
      ztemp = default_init()
      if(incx==1 .and. incy==1) then
         !   code for both increments equal to 1
         !$OMP SIMD ALIGNED(zx:64,zy) REDUCTION(+:ztemp) UNROLL PARTIAL(4)
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
         !$OMP SIMD ALIGNED(zx:64,zy) REDUCTION(+:ztemp) UNROLL PARTIAL(4)
         do i = 1,n
            ztemp = ztemp+conjugate(zx(ix))*zy(iy)
            ix = ix+incx
            iy = iy+incy
         end do
         zdotc = default_init()
      end if  
      zdotc = ztemp
    end function zdotc

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))   
    function zdotu(n,zx,incx,zy,incy) result(zdotu) !GCC$ ATTRIBUTES inline :: zdotu !GCC$ ATTRIBUTES aligned(32) :: zdotu
#if defined(__INTEL_COMPILER) || defined(__ICC)
     function zdotu(n,zx,incx,zy,incy) result(zdotu)  
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zdotu
       !DIR$ ATTRIBUTES FORCEINLINE :: zdotu
        !DIR$ OPTIMIZE : 3
        !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zdotu
#endif
      use module_kinds, only : i4
      use mod_avx512c8f64
      use omp_lib
      implicit none
      integer(kind=i4),                                intent(in),value :: n
      type(AVX512c8f64_t), dimension(:), allocatable,  intent(in)       :: zx
      integer(kind=i4),                                intent(in),value :: incx
      type(AVX512c8f64_t), dimension(:), allocatable,  intent(in)       :: zy
      integer(kind=i4),                                intent(in),value :: incy
      ! LOcals
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES ALIGN : 64 :: zdotu
#endif
      type(AVX512c8f64_t) :: zdotu
#if defined __INTEL_COMPILER
      !DIR$ ATTRIBUTES ALIGN : 64 :: ztemp
#endif
      type(AVX512c8f64_t), automatic :: ztemp
      integer(kind=i4),  automatic :: i,ix,iy
      ! Exec code ....
      if(n<=0) return
      ztemp = default_init()
      if(incx==1 .and. incy==1) then

         !  code for both increments equal to 1
          !$OMP SIMD ALIGNED(zx:64,zy) LINEAR(I:1) REDUCTION(+:ztemp) UNROLL PARTIAL(4)
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
         !$OMP SIMD ALIGNED(zx:64,zy) LINEAR(I:1) REDUCTION(+:ztemp) UNROLL PARTIAL(4)
         do i = 1,n
            ztemp = ztemp+zx(ix)*zy(iy)
            ix = ix+incx
            iy = iy+incy
         end do
      end if
      zdotu = default_init()
      zdotu = ztemp
    end function zdotu

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))   
    subroutine zdrot(n,cx,incx,cy,incy,c,s) !GCC$ ATTRIBUTES inline :: zdrot !GCC$ ATTRIBUTES aligned(32) :: zdrot
#elif defined(__ICC) || defined(__INTEL_COMPILER)
    subroutine zdrot(n,cx,incx,cy,incy,c,s)  
      !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zdrot
      !DIR$ ATTRIBUTES FORCEINLINE :: zdrot
        !DIR$ OPTIMIZE : 3
        !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zdrot
#endif
      use mod_vectypes, only : ZMM8r8_t
        use module_kinds, only : i4
      use mod_avx512c8f64
      use omp_lib
      implicit none
      integer(kind=i4),                                 intent(in),value    :: n
      type(AVX512c8f64_t), dimension(:), allocatable,   intent(inout)       :: cx
      integer(kind=i4),                                 intent(in),value    :: incx
      type(AVX512c8f64_t), dimension(:),allocatable,    intent(inout)       :: cy
      integer(kind=i4),                                 intent(in),value    :: incy
      type(ZMM8r8_t),                                   intent(in)          :: c ! scalar extended to vector
      type(ZMM8r8_t),                                   intent(in)          :: s ! scalar extended to vector
#if defined __INTEL_COMPILER     
      !DIR$ ATTRIBUTES ALIGN : 64 :: ztemp
#endif
      type(AVX512c8f64_t), automatic :: ztemp
      integer(kind=i4),  automatic :: i,ix,iy
      ! EXec code ...
      
      ctemp = default_init()
      if(incx==1 .and. incy==1) then
         !  code for both increments equal to 1
         !$OMP SIMD ALIGNED(cx:64,cy) LINEAR(i:1) UNROLL PARTIAL(4)
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
           !$OMP SIMD ALIGNED(cx:64,cy) LINEAR(i:1) UNROLL PARTIAL(4)
         do i = 1,n
            ctemp  = c*cx(ix)+s*cy(iy)
            cy(iy) = c*cy(iy)-s*cx(ix)
            cx(ix) = ctemp
            ix = ix+incx
            iy = iy+incy
         end do
       end if
    end subroutine zdrot

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))   
    subroutine zdscal(n,da,zx,incx) !GCC$ ATTRIBUTES inline :: zdscal !GCC$ ATTRIBUTES aligned(32) :: zdscal
#elif defined(__INTEL_COMPILER) || defined(__ICC)
    subroutine zdscal(n,da,zx,incx)  
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zdscal
      !DIR$ ATTRIBUTES FORCEINLINE :: zdscal
       !DIR$ OPTIMIZE : 3
        !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zdscal
#endif
        use module_kinds, only : i4
      use mod_avx512c8f64
      use omp_lib
      implicit none
       integer(kind=i4),                               intent(in),value    :: n
       type(AVX512c8f64_t),                            intent(in)          :: da
       type(AVX512c8f64_t), dimension(:), allocatable, intent(inout)       :: zx
       integer(kind=i4),                               intent(in),value    :: incx
       ! LOcals
       integer(kind=i4), automatic :: i,nincx
       ! Exec code ....
       
       if(incx==1) then
          !  code for increment equal to 1
          !$OMP SIMD ALIGNED(zx:64) LINEAR(i:1) UNROLL PARTIAL(4)
          do i = 1,n
             zx(i) = da*zx(i)
          end do
       else
          !   *        code for increment not equal to 1
          nincx=n*incx
          !$OMP SIMD ALIGNED(zx:64) LINEAR(i:1) UNROLL PARTIAL(4)
          do i = 1,nincx,incx
             zx(i) = da*zx(i)
          end do
       end if
    end subroutine zdscal

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))   
subroutine zgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy) !GCC$ ATTRIBUTES hot :: zgbmv !GCC$ ATTRIBUTES aligned(32) :: zgbmv !GCC$ ATTRIBUTES no_stack_protector :: zgbmv 
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zgbmv
  !DIR$ OPTIMIZE : 3
   !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zgbmv
#endif
          use module_kinds, only : i4
          use mod_avx512c8f64
          use omp_lib
          implicit none
        character(len=1),                                    intent(in),value :: trans
        integer(kind=i4),                                    intent(in),value :: m
        integer(kind=i4),                                    intent(in),value :: n
        integer(kind=i4),                                    intent(in),value :: kl
        integer(kind=i4),                                    intent(in),value :: ku
        type(AVX512c8f64_t),                                 intent(in)       :: alpha
        !type(AVX512c8f64_t), dimension(lda,*), intent(in)       :: a
        type(AVX512c8f64_t), dimension(:,:), allocatable,    intent(in)       :: a
        integer(kind=i4),                                    intent(in),value :: lda
       ! type(AVX512c8f64_t), dimension(*),     intent(in)       :: x
         type(AVX512c8f64_t), dimension(:), allocatable,     intent(in)       :: x
        integer(kind=i4),                                    intent(in),value :: incx
        type(AVX512c8f64_t),                                 intent(in)       :: beta
        !type(AVX512c8f64_t), dimension(*),     intent(inout)    :: y
        type(AVX512c8f64_t), dimension(:), allocatable,      intent(in)       :: x
        integer(kind=i4),                                    intent(in),value :: incy
        ! LOcals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
        !
        integer(kind=dp),    automatic :: i,info,ix,iy,j,jk,k,kup1,kx,ky,lenx,leny
        logical(kind=i4),  automatic :: noconj
        logical(kind=i1),  automatic :: beq0,aeq0,bneq0,aneq0,beq1,aeq1
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
       ! info = 0
       ! if(.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. &
       !    .not.lsame(trans,'C')) then
       !    info = 1
       ! else if(m<0) then
       !    info = 2
       ! else if(n<0) then
       !    info = 3
       ! else if(kl<0) then
       !    info = 4
       ! else if(ku<0) then
       !    info = 5
       ! else if(lda<(kl+ku+1)) then
       !    info = 8
       ! else if(incx==0) then
       !    info = 10
       ! else if(incy==0) then
       !    info = 13
       ! end if
       ! if(info/=0) then
       !    call xerbla('GMS_ZGBMV',info)
       !    return
       ! end if
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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))                 
                 !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(4)
#elif defined(__INTEL_COMPILER) || defined(__ICC)
                 !DIR$ VECTOR ALIGNED
                 !DIR$ VECTOR NONTEMPORAL
#endif
                 do i=1,leny
                    y(i) = ZERO
                 end do
              else
                 !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(4)
                 do i=1,leny
                    y(i) = beta*y(i)
                 end do
              end if
           else
              iy=ky
              if(beq0) then
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))                 
                 !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(4)
#elif defined(__INTEL_COMPILER) || defined(__ICC)
                 !DIR$ VECTOR ALIGNED
                 !DIR$ VECTOR NONTEMPORAL
#endif
                 do i=1,leny
                    y(iy) = ZERO
                    iy = iy+incy
                 end do
              else
                  !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(4)
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
              !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,k,jx) IF(n>=200)
              do j=1,n
                 temp=alpha*x(jx)
                 k=kup1-j
                 do i=max(1,j-ku),min(m,j+kl)
                    y(i) = y(i)+temp*a(k+1,j)
                 end do
                 jx = jx+incx
              end do
              !$OMP END PARALLEL DO
           else
              !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,iy,k,jx,ky) 
              do j=1,n
                 temp=alpha*x(jx)
                 iy=ky
                 k=kup1-j
                 do i=max(1,j-ku),min(m,j+kl)
                    y(iy) = y(iy)+temp*a(k+1,j)
                    iy=iy+incy
                 end do
                 jx=jx+incx
                 if(j>ku) ky=ky+incy
              end do
              !$OMP END PARALLEL DO
           end if
        else
           !  Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
           jy=ky
           if(incx==1) then
              !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,k,jy) 
              do j=1,n
                 temp=ZERO
                 k=kup1-j
                 if(noconj) then
                    do i=max(1,j-ku),min(m,j+kl)
                       temp = temp+a(k+i,j)*x(i)
                    end do
                 else
                    do i=max(1,j-ku),min(m,j+kl)
                       temp = temp+conjugate(a(k+i,j))*x(i)
                    end do
                 end if
                 y(jy) = y(jy)+alpha*temp
                 jy=jy+incy
              end do
              !$OMP END PARALLEL DO
           else
              !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,k,jy,kx) 
              do j=1,n
                 temp=ZERO
                 ix=kx
                 k=kup1-j
                 if(noconj) then
                    do i=max(1,j-ku),min(m,j+kl)
                       temp = temp+a(k+i,j)*x(ix)
                       ix = ix+incx
                    end do
                 else
                    do i=max(1,j-ku),min(m,j+kl)
                       temp = temp+conjugate(a(k+i,j)*x(ix)
                       ix = ix+incx
                    end do
                 end if
                 y(jy) = y(jy)+alpha*temp
                 jy = jy+incy
                 if(j>ku) kx = kx+incx
              end do
              !$OMP END PARALLEL DO
           end if
        end if
        
    end subroutine zgbmv

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))   
subroutine zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) !GCC$ ATTRIBUTES hot :: zgemm !GCC$ ATTRIBUTES aligned(32) :: zgemm !GCC$ ATTRIBUTES no_stack_protector :: zgemm
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zgemm
  !DIR$ OPTIMIZE : 3
   !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zgemm
#endif
          use module_kinds, only : i4,i1
          use mod_avx512c8f64
          use omp_lib
          implicit none
       character(len=1),                                 intent(in),value    :: transa
       character(len=1),                                 intent(in),value    :: transb
       integer(kind=i4),                                 intent(in),value    :: m
       integer(kind=i4),                                 intent(in),value    :: n
       integer(kind=i4),                                 intent(in),value    :: k
       type(AVX512c8f64_t),                              intent(in)          :: alpha
       !type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
       type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: a
       integer(kind=i4),                                 intent(in),value    :: lda
       !type(AVX512c8f64_t), dimension(ldb,*), intent(in)          :: b
       type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: a
       integer(kind=i4),                                 intent(in),value    :: ldb
       type(AVX512c8f64_t),                              intent(in)          :: beta
       !type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
       type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: a
       integer(kind=i4),                                 intent(in),value    :: ldc
       ! Locals
#if defined __INTEL_COMPILER
       !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
       type(AVX512c8f64_t), automatic :: temp
       integer(kind=i4),  automatic :: i,info,j,l,ncola,nrowa,nrowb
       logical(kind=i4),  automatic :: conja,conjb,nota,notb
       logical(kind=i1),  automatic :: aeq0,beq0,beq1,bneq1
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
     !   info = 0
     !   if((.not.nota) .and. (.not.conja) .and. &
      !     (.not.lsame(transa,'T'))) then
     !      info = 1
     !   else if((.not.notb) .and. (.not.conjb) .and. &
     !      (.not.lsame(transb,'T'))) then
     !      info = 2
     !   else if(m<0) then
     !      info = 3
     !   else if(n<0) then
    !       info = 4
     !   else if(k<0) then
     !      info = 5
     !   else if(lda < max(1,nrowa)) then
     !      info = 8
     !   else if(ldb < max(1,nrowb)) then
     !      info = 10
     !   else if(ldc < max(1,m)) then
    ! !      info = 13
    !    end if
    !    if(info/=0) then
    !       call xerbla('GMS_ZGEMM',info)
    !    end if
        aeq0 = all(alpha==ZERO)
        beq1 = all(beta==ONE)
        ! Early exit
        if((m==0) .or. (n==0) .or. &
             (((aeq0) .or. (k==0)) .and. (beq1))) return
        !   And when  alpha.eq.zero.
        beq0 = all(beta==ZERO)
        if(aeq0) then
           if(beq0) then
              !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
              do j=1,n
                 !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                 do i=1,m
                    c(i,j) = ZERO
                 end do
              end do
              !$OMP END PARALLEL DO
           else
               !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
              do j=1,n
                  !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                 do i=1,m
                    c(i,j) = beta*c(i,j)
                 end do
              end do
              !$OMP END PARALLEL DO
           end if
           return
        end if
        !  Start the operations.
        bneq1 = all(beta/=ONE)
        if(notb) then
           if(nota) then
              !  Form  C := alpha*A*B + beta*C.
              !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
              do j=1,n
                 if(beq0) then
                    !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=1,m
                       c(i,j) = ZERO
                    end do
                 else if(bneq0) then
                    !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=1,m
                       c(i,j) = beta*c(i,j)
                    end do
                 end if
                 do l=1,k
                    temp = alpha*b(l,j)
                    !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=1,m
                       c(i,j) = c(i,j)+temp*a(i,l)
                    end do
                 end do
              end do
              !$OMP END PARALLEL DO
           else if(conja) then
              !   Form  C := alpha*A**H*B + beta*C.
               !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
              do j=1,n
                 do i=1,m
                    temp = ZERO
                    !$OMP SIMD ALIGNED(a:64,b) LINEAR(l:1) REDUCTION(+:temp)
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
              !$OMP END PARALLEL DO
           else
              !  Form  C := alpha*A**T*B + beta*C
               !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
              do j=1,n
                 do i=1,m
                    temp = ZERO
                      !$OMP SIMD ALIGNED(a:64,b) LINEAR(l:1) REDUCTION(+:temp)
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
              !$OMP END PARALLEL DO
           end if
        else if(nota) then
               if(conjb) then
                  !  Form  C := alpha*A*B**H + beta*C.
                   !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                  do j=1,n
                     if(beq0) then
                        !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                        do i=1,m
                           c(i,j) = ZERO
                        end do
                     else if(bneq1) then
                        !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                        do i=1,m
                           c(i,j) = beta*c(i,j)
                        end do
                     end if
                     do l=1,k
                        temp = alpha*conjugate(b(j,l))
                        !$OMP SIMD ALIGNED(c:64,a) LINEAR(i:1) UNROLL PARTIAL(6)
                        do i=1,m
                           c(i,j) = c(i,j)+temp*a(i,l)
                        end do
                     end do
                  end do
                  !$OMP END PARALLEL DO
               else
                  !  Form  C := alpha*A*B**T + beta*C
                   !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                  do j=1,n
                     if(beq0) then
                          !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                        do i=1,m
                           c(i,j) = ZERO
                        end do
                     else if(bneq1) then
                          !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                        do i=1,m
                           c(i,j) = beta*c(i,j)
                        end do
                     end if
                     do l=1,k
                        temp = alpha*b(j,l)
                         !$OMP SIMD ALIGNED(c:64,a) LINEAR(i:1) UNROLL PARTIAL(6)
                        do i=1,m
                           c(i,j) = c(i,j)+temp*a(i,l)
                        end do
                     end do
                  end do
                  !$OMP END PARALLEL DO
               end if
            else if(conja) then
                  if(conjb) then
                     !   Form  C := alpha*A**H*B**H + beta*C.
                       !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                     do j=1,n
                        do i=1,m
                           temp = ZERO
                            !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
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
                     !$OMP END PARALLEL DO
                  else
                     !   Form  C := alpha*A**H*B**T + beta*C
                     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                     do j=1,n
                        do i=1,m
                           temp = ZERO
                            !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
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
                     !$OMP END PARALLEL DO
                  end if
               else
                  if(conjb) then
                     !  Form  C := alpha*A**T*B**H + beta*C
                      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                     do j=1,n
                        do i=1,m
                           temp = ZERO
                            !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
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
                     !$OMP END PARALLEL DO
                  else
                     !  Form  C := alpha*A**T*B**T + beta*C
                     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                     do j=1,n
                        do i=1,m
                           temp = ZERO
                              !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
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
                     !$OMP END PARALLEL DO
                  end if
               end if
               ! End of GMS_ZGEMM
      end subroutine zgemm

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))   
subroutine zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy) !GCC$ ATTRIBUTES hot :: zgemv !GCC$ ATTRIBUTES aligned(32) :: zgemv !GCC$ ATTRIBUTES no_stack_protector :: zgemv
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zgmev
  !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zgemv
#endif
          use module_kinds, only : i4,i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

          character(len=1),                      intent(in),value :: trans
          integer(kind=i4),                    intent(in),value :: m
          integer(kind=i4),                    intent(in),value :: n
          type(AVX512c8f64_t),                   intent(in)       :: alpha
          !type(AVX512c8f64_t), dimension(lda,*), intent(in)       :: a
          type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)       :: a
          integer(kind=i4),                    intent(in),value :: lda
          !type(AVX512c8f64_t), dimension(*),     intent(in)       :: x
          type(AVX512c8f64_t), dimension(:), allocatable,intent(in)       :: a
          integer(kind=i4),                    intent(in),value    :: incx
          type(AVX512c8f64_t),                   intent(in)          :: beta
          !type(AVX512c8f64_t), dimension(*),     intent(inout)       :: y
          type(AVX512c8f64_t), dimension(:), allocatable,intent(in)       :: a
          integer(kind=i4),                    intent(in),value    :: incy
          ! LOcals
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
          type(AVX512c8f64_t), automatic :: temp
          integer(kind=i4),  automatic :: i,info,ix,iy,j,jx.jy,kx,ky,lenx,leny
          logical(kind=i4),  automatic :: noconj
          logical(kind=i1),  automatic :: aeq0,beq1,bneq1,beq0
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
  !        info = 0
  !        if(.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. &
   !          .not.lsame(trans,'C')) then
  !!           info = 1
   !       else if(m<0) then
   !          info = 2
   !       else if(n<0) then
   !          info = 3
   !       else if(lda < max(1,m)) then
   !          info = 6
   !       else if(incx==0) then
   !          info = 8
    !      else if(incy==0) then
    !         info = 11
    !      end if
    !      if(info/=0) then
   ! !         call xerbla('GMS_ZGEMV',info)
   !          return
   !       end if
   !       aeq0  = .false.
   !       beq1  = .false.
   !       bneq1 = .false.
   !       beq0  = .false.
          ! Quick return if possible.
   !       aeq0 = all(alpha==ZERO)
   !       beq1 = all(beta==ONE)
   !       if((m==0)  .or. (n==0) .or. &
   !            ((aeq0) .and. (beq1))) return
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
                   !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=1,leny
                       y(i) = ZERO
                    end do
                 else
                    !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=1,leny
                       y(i) = beta*leny(i)
                    end do
                 end if
              else
                 iy = ky
                 if(beq0) then
                      !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=1,leny
                       y(iy) = ZERO
                       iy = iy+incy
                    end do
                 else
                      !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(6)
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
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,jx) 
                 do j=1,n
                    temp = alpha*x(jx)
                    !$OMP SIMD ALIGNED(y:64,a) LINEAR(i:1) UNROLL PARTIAL(10)
                    do i=1,m
                       y(i) = y(i)+temp*a(i,j)
                    end do
                    jx = jx+incx
                 end do
                 !$OMP END PARALLEL DO
              else
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,iy,jx) 
                 do j=1,n
                    temp = alpha*x(jx)
                    iy = ky
                     !$OMP SIMD ALIGNED(y:64,a) LINEAR(i:1) UNROLL PARTIAL(10)
                    do i=1,m
                       y(iy) = y(iy)+temp*a(i,j)
                       iy = iy+incy
                    end do
                    jx = jx+incx
                 end do
                 !$OMP END PARALLEL DO
              end if
           else
              !  Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
              jy = ky
              if(incx==1) then
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,jy)
                 do j=1,n
                    temp = zero
                    if(noconj) then
                        !$OMP SIMD ALIGNED(y:64,a) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                       do i=1,m
                          temp = temp+a(i,j)*x(i)
                       end do
                    else
                        !$OMP SIMD ALIGNED(y:64,a) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                       do i=1,m
                          temp = temp+conjugate(a(i,j))*x(i)
                       end do
                    end if
                    y(jy) = y(jy)+alpha*temp
                    jy = jy+incy
                 end do
                 !$OMP END PARALLEL DO
              else
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jy) 
                 do j=1,n
                    temp = zero
                    ix = ky
                    if (noconj) then
                        !$OMP SIMD ALIGNED(y:64,a) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                       do i=1,m
                          temp = temp+a(i,j)*x(ix)
                          ix = ix+incx
                       end do
                    else
                        !$OMP SIMD ALIGNED(y:64,a) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                       do i=1,m
                          temp = temp+conjugate(a(i,j))*x(ix)
                          ix = ix+incx
                       end do
                    end if
                       y(jy) = y(jy)+alpha*temp
                       jy = jy+incy
                    end do
                    !$OMP END PARALLEL DO
               end if
            end if
            ! End of ZGEMV
     end subroutine zgemv

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zgerc(m,n,alpha,x,incx,y,incy,a,lda) !GCC$ ATTRIBUTES hot :: zgerc !GCC$ ATTRIBUTES aligned(32) :: zgerc !GCC$ ATTRIBUTES no_stack_protector :: zgerc
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zgerc(m,n,alpha,x,incx,y,incy,a,lda)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zgerc
   !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zgerc
#endif
          use module_kinds, only : i4
          use mod_avx512c8f64
          use omp_lib
          implicit none

        integer(kind=int4),                    intent(in),value    :: m
        integer(kind=int4),                    intent(in),value    :: n
        type(AVX512c8f64_t),                   intent(in)          :: alpha
        !type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
        type(AVX512c8f64_t), dimension(:), allocatable,     intent(in)          :: x
        integer(kind=int4),                    intent(in),value    :: incx
        !type(AVX512c8f64_t), dimension(*),     intent(in)          :: y
        type(AVX512c8f64_t), dimension(:), allocatable,    intent(in)          :: y
        integer(kind=int4),                    intent(in),value    :: incy
        !type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
        type(AVX512c8f64_t), dimension(:,:), allocatable,    intent(in)          :: x
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
      !  info = 0
      !  if(m<0) then
       !    info = 1
       ! else if(n<0) then
       !    info = 2
       ! else if(incx==0) then
      !     info = 5
       ! else if(incy==0) then
       !    info = 7
      !!  else if(lds < max(1,m)) then
      !     info = 9
      !  end if
      !  if(info/=0) then
     !!      call xerbla('ZGERC',info)
      !     return
     !   end if
        !  Quick return if possible.
      !  aeq0 = .false.
      !  aeq0 = all(alpha==ZERO)
      !  if((m==0) .or. (n==0) .or. (aeq0)) return
        ! Start the operations. In this version the elements of A are
        ! *     accessed sequentially with one pass through A.
        if(incy>0) then
           jy = 1
        else
           jy = 1-(n-1)*incy
        end if
        if(incx==1) then
           !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,jy) 
           do j=1,n
              if(all(y(jy)/=ZERO)) then
                 temp = alpha*conjugate(y(jy))
                 !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
                 do i=1,m
                    a(i,j) = a(i,j)+x(i)*temp
                 end do
              end if
              jy = jy+incy
           end do
           !$OMP END PARALLEL DO
        else
           if(incx>0) then
              kx = 1
           else
              kx = 1-(m-1)*incx
           end if
           !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,jy)
           do j=1,n
              if(all(y(jy)/=ZERO)) then
                 temp = alpha*conjugate(y(jy))
                 ix = kx
                  !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
                 do i=1,m
                    a(i,j) = a(i,j)+x(ix)*temp
                    ix = ix+incx
                 end do
              end if
              jy = jy+incy
           end do
        end if
        ! End of zgerc
     end subroutine zgerc

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zgeru(m,n,alpha,x,incx,y,incy,a,lda) !GCC$ ATTRIBUTES hot :: zgeru !GCC$ ATTRIBUTES aligned(32) :: zgeru !GCC$ ATTRIBUTES no_stack_protector :: zgeru
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zgeru(m,n,alpha,x,incx,y,incy,a,lda)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zgeru
    !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zgeru
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

        integer(kind=i4),                    intent(in),value    :: m
        integer(kind=i4),                    intent(in),value    :: n
        type(AVX512c8f64_t),                   intent(in)          :: alpha
        !type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
        type(AVX512c8f64_t), dimension(:), allocatable,    intent(in)          :: x
        integer(kind=i4),                    intent(in),value    :: incx
        !type(AVX512c8f64_t), dimension(*),     intent(in)          :: y
        type(AVX512c8f64_t), dimension(:), allocatable,     intent(in)          :: y
        integer(kind=i4),                    intent(in),value    :: incy
        !type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
        type(AVX512c8f64_t), dimension(:,:), allocatable,     intent(in)          :: x
        integer(kind=i4),                    intent(in),value    :: lda
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
        integer(kind=i4),  automatic :: i,info,ix,j,jy,kx
        logical(kind=i1),  automatic :: aeq0
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        !  Test the input parameters.
!        info = 0
!        aeq0 = .false.
 !       if(m<0) then
!           info = 1
!        else if(n<0) then
 !          info = 2
 !       else if(incx==0) then
 !          info = 5
 !       else if(incy==0) then
 !          info = 7
 !       else if(lda < max(1,m)) then
 !          info = 9
 !       end if
 !       if(info/=0) then
 !          call xerbla('GMS_ZGERU',info)
 !!          return
 !       end if
        !  Quick return if possible.
 !       aeq0 = all(alpha==ZERO)
 !       if((m==0) .or. (n==0) .or. (aeq0)) return
        !   Start the operations. In this version the elements of A are
        ! *     accessed sequentially with one pass through A.
        if(incy>0) then
           jy = 1
        else
           jy = 1-(n-1)*incy
        end if
        if(incx==1) then
           !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,jy) 
           do j=1,n
              if(all(y(jy)/=ZERO)) then
                 temp = alpha*y(jy)
                 !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
                 do i=1,m
                    a(i,j) = a(i,j)+x(i)*temp
                 end do
              end if
              jy = jy+incy
           end do
           !$OMP END PARALLEL DO
        else
           if(incx>0) then
              kx=1
           else
              kx = 1-(m-1)*incx
           end if
            !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jy)
           do j=1,n
              if(all(y(jy)/=ZERO)) then
                 temp = alpha*y(jy)
                 ix = kx
                   !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
                 do i=1,m
                    a(i,j) = a(i,j)+x(ix)*temp
                    ix = ix+incx
                 end do
              end if
              jy = jy+incy
           end do
           !$OMP END PARALLEL DO
        end if
        ! End of ZGERU
     end subroutine zgeru

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))       
subroutine zhbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy) !GCC$ ATTRIBUTES hot :: zhbmv !GCC$ ATTRIBUTES aligned(32) :: zhbmv !GCC$ ATTRIBUTES no_stack_protector :: zhbmv
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zhbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy) 
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zhbmv
   !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zhbmv
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

        character(len=1),                      intent(in),value    :: uplo
        integer(kind=i4),                    intent(in),value    :: n
        integer(kind=i4),                    intent(in),value    :: k
        type(AVX512c8f64_t),                   intent(in)          :: alpha
        !type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
         type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: a
        integer(kind=i4),                    intent(in),value    :: lda
        !type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
        type(AVX512c8f64_t), dimension(:), allocatable,    intent(in)          :: x
        integer(kind=i4),                    intent(in),value    :: incx
        type(AVX512c8f64_t),                   intent(in)          :: beta
        !type(AVX512c8f64_t), dimension(*),     intent(inout)       :: y
        type(AVX512c8f64_t), dimension(:), allocatable,    intent(in)          :: y
        integer(kind=i4),                    intent(in),value    :: incy
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
        type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
        type(AVX512c8f64_t), automatic :: temp2
        integer(kind=i4),  automatic :: i,info,ix,iy,j,jx,jy,kplus1,kx,ky,l
        logical(kind=i1),  automatic :: aeq0,beq1,beq0,bneq1
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
        
        !info = 1
        !if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
        !   info = 1
        !else if(n<0) then
        !   info = 2
        !else if(k<0) then
        !   info = 3
        !else if(lda<(k+1)) then
        !   info = 6
       ! else if(incx==0) then
        !   info = 8
       ! else if(incy==0) then
        !   info = 11
       ! end if
        !if(info/=0) then
        !   call xerbla('GMS_ZHBMV',info)
       !    return
        !end if
        ! Quick return if possible
       ! aeq0 = .false.
       ! beq1 = .false.
       ! aeq0 = all(alpha==ZERO)
       ! beq1 = all(beta==ONE)
       ! if((n==0) .or. ((aeq0) .and. (beq1))) return
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
                   !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(8)
                   do i=1,n
                      y(i) = ZERO
                   end do
                else
                    !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(6)
                   do i=1,m
                      y(i) = beta*y(i)
                   end do
                end if
             else
                iy = ky
                if(beq0) then
                    !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(8)
                   do i=1,n
                      y(iy) = ZERO
                      iy = iy+incy
                   end do
                else
                    !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(6)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,l) 
                do j=1,n
                   temp1 = alpha*x(j)
                   temp2 = ZERO
                   l = kplus1-j
                   !$OMP SIMD ALIGNED(y:64,a) LINEAR(i:1) UNROLL PARTIAL(10) REDUCTION(+:temp2)
                   do i=max(1,j-k),j-1
                      y(i) = y(i)+temp1*a(l+i,j)
                      temp2 = temp2+conjugate(a(l+i,j))*x(i)
                   end do
                   y(j) = y(j)+temp1*a(kplus1,j).re+alpha*temp2
                end do
                !$OMP END PARALLEL DO
             else
                jx = kx
                jy = ky
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,ix,iy,jx,jy,kx,ky,l) 
                do j=1,n
                   temp1 = alpha*x(jx)
                   temp2 = ZERO
                   ix = kx
                   iy = ky
                   l = kplus1-j
                   !$OMP  SIMD ALIGNED(y:64,a) LINEAR(i:1) UNROLL PARTIAL(10) REDUCTION(+:temp2)
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
                !$OMP END PARALLEL DO
             end if
          else
             !  Form  y  when lower triangle of A is stored.
             if((incx==1) .and. (incy==1)) then
                   !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,l) 
                do j=1,n
                   temp1 = alpha*x(j)
                   temp2 = ZERO
                   y(j) = y(j)+temp1*a(1,j).re
                   l = 1-j
                    !$OMP  SIMD ALIGNED(y:64,a) LINEAR(i:1) UNROLL PARTIAL(10) REDUCTION(+:temp2)
                   do i=j+1,min(n,j+k)
                      y(i) = y(i)+temp1*a(l+i,j)
                      temp2 = temp2+conjugate(a(l+i,j))*x(i)
                   end do
                   y(j) = y(j)+alpha*temp2
                end do
                !$OMP END PARALLEL DO
             else
                jx = kx
                jy = ky
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,ix,iy,jx,jy,l) 
                do j=1,n
                   temp1 = alpha*x(jx)
                   temp2 = ZERO
                   y(jy) = y(jy)+temp1*a(1,j).re
                   l = 1-j
                   ix = jx
                   iy = jy
                   !$OMP  SIMD ALIGNED(y:64,a) LINEAR(i:1) UNROLL PARTIAL(10) REDUCTION(+:temp2)
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
                !$OMP END PARALLEL DO
             end if
          end if
          ! End of zhbmv
      end subroutine zhbmv

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zhemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc) !GCC$ ATTRIBUTES hot :: zhemm !GCC$ ATTRIBUTES aligned(32) :: zhemm !GCC$ ATTRIBUTES no_stack_protector :: zhemm
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zhemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zhemm
    !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zhemm
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

          character(len=1),                      intent(in),value    :: side
          character(len=1),                      intent(in),value    :: uplo
          integer(kind=i4),                    intent(in),value    :: m
          integer(kind=i4),                    intent(in),value    :: n
          type(AVX512c8f64_t),                   intent(in)          :: alpha
          !type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
          type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: a
          integer(kind=i4),                    intent(in),value    :: lda
          !type(AVX512c8f64_t), dimension(ldb,*), intent(in)          :: b
          type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: a
          integer(kind=i4),                    intent(in),value    :: ldb
          type(AVX512c8f64_t),                   intent(in)          :: beta
          !type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
          type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: a
          integer(kind=i4),                    intent(in),value    :: ldc
          ! LOcals
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: temp1,temp2
#endif
          type(AVX512c8f64_t), automatic :: temp1,temp2
          integer(kind=i4),  automatic :: i,info,j,k,nrowa
          logical(kind=i4),  automatic :: upper
          logical(kind=i1),  automatic :: aeq0,beq1,beq0
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
          !aeq0 = .false.
          !beq1 = .false.
          beq0 = .false.
          info = 0
         ! if((.not.lsame(side,'L')) .and. (.not.lsame(side,'R'))) then
          !   info = 1
         ! else if((.not.upper) .and. (not.lsame(uplo,'L'))) then
          !   info = 2
          !else if(m<0) then
          !   info = 3
          !else if(n<0) then
          !   info = 4
          !else if(lda < max(1,nrowa)) then
          !   info = 7
          !else if(ldb < max(1,m)) then
          !   info = 9
          !else if(ldc < max(1,m)) then
          !   info = 13
          !end if
          !if(info/=0) then
          !   call xerbla('GMS_ZHEMM',info)
          !   return
          !end if
          !  Quick return if possible.
          aeq0 = all(alpha==ZERO)
          beq1 = all(beta==ONE)
          if((m==0) .or. (n==0) .or. &
               ((aeq0) .and. (beq1))) return
          !  And when  alpha.eq.zero.
          beq0 = all(beta==ZERO)
          if(aeq0) then
             if(beq0) then
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                do j=1,n
                   !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                   do i=1,m
                      c(i,j) = ZERO
                   end do
                end do
                !$OMP END PARALLEL DO
             else
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                do j=1,n
                    !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                   do i=1,m
                      c(i,j) = beta*c(i,j)
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
             return
          end if
          !  Start the operations.
          if(lsame(side,'L')) then
             !  Form  C := alpha*A*B + beta*C.
             if(upper) then
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2) 
                do j=1,n
                   do i=1,m
                      temp1 = alpha*b(i,j)
                      temp2 = ZERO
                      !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1) REDUCTION(+:temp2) UNROLL PARTIAL(10)
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
                !$OMP END PARALLEL DO
             else
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2) 
                do j=1,n
                   do i=m,1,-1
                      temp1 = alpha*b(i,j)
                      temp2 = ZERO
                       !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1) REDUCTION(+:temp2) UNROLL PARTIAL(10)
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
                !$OMP END PARALLEL DO
             end if
          else
             !   Form  C := alpha*B*A + beta*C.
             !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1)
             do j=1,n
                temp1 = alpha*a(j,j).re
                if(beq0) then
                     !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1)  UNROLL PARTIAL(6)
                   do i=1,m
                      c(i,j) = temp*b(i,j)
                   end do
                else
                   !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1)  UNROLL PARTIAL(6)
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
                   !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1)  UNROLL PARTIAL(6)
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
                   !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1)  UNROLL PARTIAL(6)
                   do i=1,m
                      c(i,j) = c(i,j)+temp1*b(i,k)
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          end if
          !End of ZHEMM
     end subroutine zhemm

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zhemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy) !GCC$ ATTRIBUTES hot :: zhemv !GCC$ ATTRIBUTES aligned(32) :: zhemv !GCC$ ATTRIBUTES no_stack_protector :: zhemv
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zhemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zhemv
    !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zhemv
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

        character(len=1),                       intent(in),value    :: uplo
        integer(kind=i4),                     intent(in),value    :: n
        type(AVX512c8f64_t),                    intent(in)          :: alpha
        !type(AVX512c8f64_t), dimension(lda,*),  intent(in)          :: a
        type(AVX512c8f64_t), dimension(:,:), allocatable,  intent(in)          :: a
        integer(kind=i4),                     intent(in),value    :: lda
        !type(AVX512c8f64_t), dimension(*),      intent(in)          :: x
        type(AVX512c8f64_t), dimension(:),  allocatable,  intent(in)          :: a
        integer(kind=i4),                     intent(in),value    :: incx
        type(AVX512c8f64_t),                    intent(in)          :: beta
        !type(AVX512c8f64_t), dimension(*),      intent(inout)       :: y
        type(AVX512c8f64_t), dimension(:),  allocatable, intent(in)          :: a
        integer(kind=i4),                     intent(in),value    :: incy
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
        type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
        type(AVX512c8f64_t), automatic :: temp2
        integer(kind=i4),  automatic :: i,info,ix,iy,j,jx,jy,kx,ky
        logical(kind=i1),  automatic :: aeq0,beq1,bneq1
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
       ! aeq0  = .false.
        !beq1  = .false.
        bneq1 = .false.
        !if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
        !   info = 1
        !else if(n<0) then
        !   info = 2
        !else if(lda < max(n,n)) then
        !   info = 5
        !else if(incx==0) then
        !   info = 7
        !else if(incy==0) then
        !   info = 10
        !end if
        !if(info/=0) then
        !   call xerbla('GMS_ZHEMV',info)
        !   return
        !end if
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
                    !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=1,n
                       y(i) = ZERO
                    end do
                 else
                    !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=1,n
                       y(i) = beta*y(i)
                    end do
                 end if
              else
                 iy = ky
                 if(beq0) then
                 !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=1,n
                       y(iy) = ZERO
                       iy = iy+incy
                    end do
                 else
                    !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(6)
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
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2) 
                 do j=1,n
                    temp1 = alpha*x(j)
                    temp2 = ZERO
                     !$OMP SIMD ALIGNED(y:64,a,x) LINEAR(i:1) REDIUCTION(+:temp2) UNROLL PARTIAL(6)
                    do i=1,j-1
                       y(i) = y(i)+temp1*a(i,j)
                       temp2 = temp2+conjugate(a(i,j))*x(i)
                    end do
                    y(j) = y(j)+temp1*a(j,j).re+alpha*temp2
                 end do
                 !$OMP END PARALLEL DO
              else
                 jx = kx
                 jy = ky
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,ix,iy,jx,jy) 
                 do j=1,n
                    temp1 = alpha*x(jx)
                    temp2 = ZERO
                    ix = kx
                    iy = ky
                    !$OMP SIMD ALIGNED(y:64,a,x) LINEAR(i:1) REDIUCTION(+:temp2) UNROLL PARTIAL(6)
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
                 !$OMP END PARALLEL DO
              end if
           else
              !  Form  y  when A is stored in lower triangle.
              if((incx==1) .and. (incy==1)) then
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2)
                 do j=1,n
                    temp1 = alpha*x(j)
                    temp2 = ZERO
                    y(j) = y(j)+temp1*a(j,j).re
                     !$OMP SIMD ALIGNED(y:64,a,x) LINEAR(i:1) REDIUCTION(+:temp2) UNROLL PARTIAL(6)
                    do i=j+1,n
                       y(i) = y(i)+temp1*a(i,j)
                       temp2 = temp2+conjugate(a(i,j))*x(i)
                    end do
                    y(j) = y(j)+alpha*temp2
                 end do
                 !OMP END PARALLEL DO
              else
                 jx = kx
                 jy = ky
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,ix,iy,jx,jy) 
                 do j=1,n
                    temp1 = alpha*x(jx)
                    temp2 = ZERO
                    y(jy) = y(jy)+temp1*a(j,j).re
                    ix = jx
                    iy = jy
                      !$OMP SIMD ALIGNED(y:64,a,x) LINEAR(i:1) REDIUCTION(+:temp2) UNROLL PARTIAL(6)
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
                 !$OMP END PARALLEL DO
              end if
           end if
           !End of ZHEMV
     end subroutine zhemv

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zher(uplo,n,alpha,x,incx,a,lda) !GCC$ ATTRIBUTES hot :: zher !GCC$ ATTRIBUTES aligned(32) :: zher !GCC$ ATTRIBUTES no_stack_protector :: zher
#if defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zher(uplo,n,alpha,x,incx,a,lda)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zher
   !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zhemr
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

        character(len=1),                      intent(in),value    :: uplo
        integer(kind=i4),                    intent(in),value    :: n
        type(ZMM8r4_t),                        intent(in)          :: alpha
        !type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
        type(AVX512c8f64_t), dimension(:), allocatable,    intent(in)          :: x
        integer(kind=i4),                    intent(in),value    :: incx
        !type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
        type(AVX512c8f64_t), dimension(:,:),  allocatable,   intent(in)          :: a
        integer(kind=i4),                    intent(in),value    :: lda
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
        integer(kind=i4),  automatic :: i,info,ix,j,jx,kx
        logical(kind=i1),  automatic :: aeq0
        !
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
       ! aeq0 = .false.
        info = 0
       ! if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
       !    info = 1
       ! else if(n<0) then
       !    info = 2
       ! else if(incx==0) then
       !    info = 5
       ! else if(lda < max(1,n)) then
       !    info = 7
       ! end if
       ! if(info/=0) then
       !    call xerbla('GMS_ZHER',info)
        !   return
       ! end if
        !  Quick return if possible.
       ! aeq0 = all(alpha==ZERO.re)
       ! if((n==0) .or. (aeq0)) return
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
              !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
              do j=1,n
                 if(all(x(j)/=ZERO)) then
                    temp = alpha*conjugate(x(j))
                    !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
                    do i=1,j-1
                       a(i,j) = a(i,j)+x(i)*temp
                    end do
                    a(j,j) = a(j,j).re+x(j).re*temp.re
                 else
                    a(j,j) = a(j,j).re
                 end if
              end do
              !$OMP END PARALLEL DO
           else
              jx = kx
              !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix) 
              do j=1,n
                 if(all(x(jx)/=ZERO)) then
                    temp = alpha*conjugate(x(jx))
                    ix = kx
                     !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
                    do i=1,j-1
                       a(i,j) = a(i,j)+x(ix)*temp
                       ix = ix+incx
                    end do
                    a(j,j) = a(j,j).re+x(j).re*temp.re
                 else
                    a(j,j) = a(j,j).re
                 end if
              end do
              !$OMP END PARALLEL DO
           else
              jx = kx
               !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jx) 
              do j=1,n
                 if(all(x(jx)/=ZERO)) then
                    temp = alpha*conjugate(x(jx))
                    ix = kx
                    !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
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
              !$OMP END PARALLEL DO
           end if
        else
           !  Form  A  when A is stored in lower triangle.
           if(incx==1) then
               !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp) 
              do j=1,n
                 if(all(x(j)/=ZERO)) then
                    temp = alpha*conjugate(x(j))
                    a(j,j) = a(j,j).re+temp.re*x(j).re
                      !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
                    do i=j+1,n
                       a(i,j) = a(i,j)+x(i)*temp
                    end do
                 else
                    a(j,j) = a(j,j).re
                 end if
              end do
              !$OMP END PARALLEL DO
           else
              jx = kx
               !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jx) 
              do j=1,n
                 if(all(x(jx)/=ZERO)) then
                    temp = alpha*conjugate(x(jx))
                    a(j,j) = a(j,j).re+temp.re*x(jx).re
                    ix = jx
                      !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
                    do i=j+1,n
                       ix = ix+incx
                       a(i,j) = a(i,j)+x(ix)*temp
                    end do
                 else
                    a(j,j) = a(j,j).re
                 end if
                 jx = jx+incx
              end do
              !$OMP END PARALLEL DO
           end if
        end if
        ! End of ZHER
    end subroutine zher

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zscal(n,za,zx,incx) !GCC$ ATTRIBUTES inline :: zscal !GCC$ ATTRIBUTES aligned(32) :: zscal
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zscal(n,za,zx,incx)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zscal
  !DIR$ ATTRIBUTES FORCEINLINE :: zscal
      !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zscal
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none
          
       integer(kind=i4),                     intent(in),value    :: n
       type(AVX512c8f64_t),                    intent(in)          :: za
       !type(AVX512c8f64_t), dimension(*),      intent(inout)       :: zx
       type(AVX512c8f64_t), dimension(:), allocatable,    intent(inout)       :: zx
       integer(kind=i4),                     intent(in),value    :: incx
       ! Locals
       integer(kind=i4), automatic :: i,nincx
       ! Exec code ....
       
       if(incx==1) then
          !$OMP SIMD ALIGNED(zx:64) LINEAR(i:1) UNROLL PARTIAL(6)
          do i=1,n
             zx(i) = za*zx(i)
          end do
       else
          !  code for increment not equal to 1
          nincx = n*incx
#if defined __INTEL_COMPILER
           !$OMP SIMD ALIGNED(zx:64)  UNROLL PARTIAL(6)
          do i=1,nincx,incx
             zx(i) = za*zx(i)
          end do
       end if 
     end subroutine zscal

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zher2(uplo,n,alpha,x,incx,y,incy,a,lda) !GCC$ ATTRIBUTES hot :: zher2 !GCC$ ATTRIBUTES aligned(32) :: zher2 !GCC$ ATTRIBUTES no_stack_protector :: zher2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zher2(uplo,n,alpha,x,incx,y,incy,a,lda)
  !DIR$ ATTRIBUTES CODE_ALIGN : 64 :: zher2
     !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zher2
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

         character(len=i1),                      intent(in),value    :: uplo
         integer(kind=i4),                    intent(in),value    :: n
         type(AVX512c8f64_t),                   intent(in)          :: alpha
         !type(AVX512c8f64_t), dimension(*),     intent(in)          :: x
          type(AVX512c8f64_t), dimension(:), allocatable,     intent(in)          :: x
         integer(kind=i4),                    intent(in),value    :: incx
         !type(AVX512c8f64_t), dimension(*),     intent(in)          :: y
         type(AVX512c8f64_t), dimension(:), allocatable,    intent(in)          :: y
         integer(kind=i4),                    intent(in),value    :: incy
         !type(AVX512c8f64_t), dimension(lda,*), intent(inout)       :: a
         type(AVX512c8f64_t), dimension(:,:), allocatable,     intent(in)          :: x 
         integer(kind=i4),                    intent(in),value    :: lda
         ! LOcals ....
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
         type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
         type(AVX512c8f64_t), automatic :: temp2
         integer(kind=i4),  automatic :: i,info,ix,iy,j,jx,jy,kx,ky
         logical(kind=i1),  automatic :: aeq0
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
         type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
         ! EXec code ....
       !  aeq0 = .false.
       !  info = 0
       !  if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
       !     info = 1
       !  else if(n<0) then
       !     info = 2
       !  else if(incx==0) then
       !     info = 5
       !  else if(incy==0) then
       !     info = 7
       !  else if(lda < max(1,n)) then
       !     info = 9
       !  end if
       !  if(info/=0) then
       !     call xerbla('GMS_ZHER2',info)
       !     return
       !  end if
         !  Quick return if possible.
       !  aeq0 = all(alpha==ZERO)
       !  if((n==0) .or. (aeq0)) return
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
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2)
               do j=1,n
                  if((all(x(j)/=ZERO)) .or. (all(y(j)/=ZERO))) then
                     temp1 = alpha*conjugate(y(j))
                     temp2 = conjugate(alpha*x(j))
                     !$OMP SIMD ALIGNED(a:64,x,y) LINEAR(i:1) UNROLL PARTIAL(10)
                     do i=1,j-1
                        a(i,j) = a(i,j)+x(i)*temp1+y(i)*temp2
                     end do
                     a(j,j) = a(j,j).re+x(j).re*temp1.re+x(j).re*temp2.re
                  else
                     a(j,j) = a(j,j).re
                  end if
               end do
               !$OMP END PARALLEL DO
            else
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,ix,iy,jx,jy) 
               do j=1,n
                  if((all(x(jx)/=ZERO)) .or. (all(y(jy)/=ZERO))) then
                     temp1 = alpha*conjugate(y(jy))
                     temp2 = conjugate(alpha*x(jx))
                     ix = kx
                     iy = ky
                       !$OMP SIMD ALIGNED(a:64,x,y) LINEAR(i:1) 
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
               !$OMP END PARALLEL DO
            end if
         else
            !  Form  A  when A is stored in the lower triangle.
            if((incx==1) .and. (incy==1)) then
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2) 
               do j=1,n
                  if((all(x(j)/=ZERO)) .or. (all(y(j)/=ZERO))) then
                     temp1 = alpha*conjugate(y(j))
                     temp2 = conjugate(alpha*x(j))
                     a(j,j) = a(j,j).re+ &
                          x(j).re*temp1.re+y(j).re*temp2.re
                      !$OMP SIMD ALIGNED(a:64,x,y) LINEAR(i:1) UNROLL PARTIAL(10)
                     do i=j+1,n
                        a(i,j) = a(i,j)+x(i)*temp1+y(i)*temp2
                     end do
                  else
                     a(j,j) = a(j,j).re
                  end if
               end do
               !$OMP END PARALLEL DO
            else
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,ix,iy,jx,jy) 
               do j=1,
                  if((all(x(jx)/=ZERO)) .or. (all(y(jx)/=ZERO))) then
                     temp1 = alpha*conjugate(y(jy))
                     temp2 = conjugate(alpha*x(jx))
                     a(j,j) = a(j,j).re + &
                          x(jx).re*temp1.re+y(jy).re*temp2.re
                     ix = jx
                     iy = jy
                      !$OMP SIMD ALIGNED(a:64,x,y) LINEAR(i:1)
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
               !$OMP END PARALLEL DO
            end if
         end if
         ! End of ZHER2
       end subroutine zher2

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc) !GCC$ ATTRIBUTES hot :: zherk2 !GCC$ ATTRIBUTES aligned(32) :: zherk2 !GCC$ ATTRIBUTES no_stack_protector :: zherk2
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zher2k
    !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zher2k
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

          use mod_vecconst, only : v8_n1
          character(len=1),                      intent(in),value    :: uplo
          character(len=1),                      intent(in),value    :: trans
          integer(kind=i4),                    intent(in),value    :: n
          integer(kind=i4),                    intent(in),value    :: k
          type(AVX512c8f64_t),                   intent(in)          :: alpha
          !type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
          type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: a
          integer(kind=i4),                    intent(in),value    :: lda
          !type(AVX512c8f64_t), dimension(ldb,*), intent(in)          :: b
           type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: b
          integer(kind=i4),                    intent(in),value    :: ldb
          type(ZMM8r8_t),                        intent(in)          :: beta
          !type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
          type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: c
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
          integer(kind=i4),  automatic :: i,info,j,l,nrowa
          logical(kind=i4),  automatic :: upper
          logical(kind=i1),  automatic :: aeq0,beq1,beq0,bneq1
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
          type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
          ! EXec code ....
          ! Test the input parameters.
          
          !aeq0  = .false.
          !beq1  = .false.
          beq0  = .false.
          bneq1 = .false.
          if(lsame(trans,'N')) then
             nrowa = n
          else
             nrowa = k
          end if
          upper = lsame(uplo,'U')
         ! info = 0
         ! if((.not.upper) .and. (.not.lsame(uplo,'L'))) then
         !    info = 1
         ! else if((.not.lsame(trans,'N')) .and. &
          !     (.not.lsame(trans,'C'))) then
          !   info = 2
         ! else if(n<0) then
         !    info = 3
         ! else if(k<0) then
         !    info = 4
         !! else if(lda < max(1,nrowa)) then
         !    info = 7
         ! else if(ldb < max(1,nrowa)) then
         !    info = 9
        !  else if(ldc < max(1,n)) then
        !     info = 12
        !  end if
        !  if(info/=0) then
       !!      call xerbla('GMS_ZHER2K',info)
        !     return
       !   end if
          !  Quick return if possible.
         ! aeq0 = all(alpha==ZERO)
         ! beq1 = all(beta.v==v8_n1.v)
         ! if((n==0) .or. ((aeq0) .or. &
         !      (k==0) .and. (beq1))) return
          !  And when  alpha.eq.zero.
          beq0 = all(beta.v==ZERO.re)
          if(aeq0) then
             if(upper) then
                if(beq0)  then
                   !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                   do j=1,n
                      !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                      do i=1,j
                         c(i,j) = ZERO
                      end do
                   end do
                   !$OMP END PARALLEL DO
                else
                     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                   do j=1,n
                       !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                      do i=1,j-1
                         c(i,j) = beta*c(i,j)
                      end do
                      c(j,j) = beta.v*c(i,j).re
                   end do
                   !$OMP END PARALLEL DO
                end if
             else
                if(beq0) then
                   !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                   do j=1,n
                       !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                      do i=j,n
                         c(i,j) = ZERO
                      end do
                   end do
                   !$OMP END PARALLEL DO
                else
                    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                   do j=1,n
                      c(j,j) = beta.v*c(j,j).re
                       !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                      do i=j+1,n
                         c(i,j) = beta*c(i,j)
                      end do
                   end do
                   !$OMP END PARALLEL DO
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
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                do j=1,n
                   if(beq0) then
                        !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                      do i=1,j
                         c(i,j) = ZERO
                      end do
                   else if(bneq1) then
                        !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
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
                          !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1) UNROLL PARTIAL(6)
                         do i=1,j-1
                            c(i,j) = c(i,j)+a(i,l)*temp1+ &
                                 b(i,l)*temp2
                         end do
                         c(j,j) = c(j,j).re+a(j,l).re*temp1.re+&
                              b(j,l).re*temp2.re
                      end if
                   end do
                end do
                !$OMP END PARALLEL DO
             else
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                do j=1,n
                   if(beq0) then
                         !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                      do i=j,n
                         c(i,j) = ZERO
                      end do
                   else if(bneq1) then
                         !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
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
                            !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1) UNROLL PARTIAL(6)
                         do i=j+1,n
                            c(i,j) = c(i,j)+a(i,l)*temp1+&
                                 b(i,l)*temp2
                         end do
                         c(j,j) = c(j,j).re+a(j,l).re*temp1.re+&
                              b(j,l).re*temp2.re
                      end if
                   end do
                end do
                !$OMP END PARALLEL DO
             end if
          else
             !   Form  C := alpha*A**H*B + conjg( alpha )*B**H*A +
             if(upper) then
                temp3 = default_init()
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                do j=1,n
                   do i=1,j
                      temp1 = ZERO
                      temp2 = ZERO
                        !$OMP SIMD ALIGNED(a:64,b) REDUCTION(+:temp1,temp2) LINEAR(i:1) UNROLL PARTIAL(10)
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
                !$OMP END PARALLEL DO
             else
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                do j=1,n
                   do i=1,n
                      temp1 = ZERO
                      temp2 = ZERO
                       !$OMP SIMD ALIGNED(a:64,b) REDUCTION(+:temp1,temp2) LINEAR(i:1) UNROLL PARTIAL(10)
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
                !$OMP END PARALLEL DO
             end if
          end if
          ! End of ZHER2K
     end subroutine zher2k

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zherk(uplo,trans,n,alpha,a,lda,beta,c,ldc) !GCC$ ATTRIBUTES hot :: zherk !GCC$ ATTRIBUTES aligned(32) :: zherk !GCC$ ATTRIBUTES no_stack_protector :: zherk
#elif defined(__INTEL_COMPILER) || defined(__ICC)
 subroutine zherk(uplo,trans,n,alpha,a,lda,beta,c,ldc) 
   !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zherk
     !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zherk
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

        use mod_vecconst, only : v8_n1,v8_n0
        character(len=1),                      intent(in),value    :: uplo
        character(len=1),                      intent(in),value    :: trans
        integer(kind=i4),                    intent(in),value    :: n
        type(ZMM8r8_t),                        intent(in)          :: alpha
        !type(AVX512c8f64_t), dimension(lda,*), intent(in)          :: a
         type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: a
        integer(kind=i4),                    intent(in),value    :: lda
        type(ZMM8r8_t),                        intent(in)          :: beta
        !type(AVX512c8f64_t), dimension(ldc,*), intent(inout)       :: c
         type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)          :: c
        integer(kind=i4),                    intent(in),value    :: ldc
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: rtemp
#endif
        type(ZMM8r8_t),      automatic :: rtemp
        integer(kind=i4),  automatic :: i,info,j,l,nrowa
        logical(kind=i4),  automatic :: upper
        logical(kind=i1),  automatic :: aeq0,beq1,beq0,bneq1
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
        !info = 0
        !if((.not.upper) .and. (.not.lsame(uplo,'L'))) then
        !   info = 1
        !else if((.not.lsame(trans,'N')) .and. &
        !     (.not.lsame(trans,'C'))) then
        !   info = 2
        !else if(n<0) then
        !   info = 3
        !else if(k<0) then
        !   info = 4
        !else if(lda < max(1,nrowa)) then
        !!   info = 7
        !else if(ldc < max(1,n)) then
        !   info = 10
        !endif
        !if(info/=0) then
        !!   call xerbla('GMS_ZHERK',info)
        !   return
        !end if
        !aeq0 = .false.
        !beq1 = .false.
        !aeq0 = all(alpha.v==v8_n0.v)
        !beq1 = all(beta.v==v8_n1.v)
        !if((n==0) .or. (aeq0) .or. &
        !     ((k==0) .and. (beq1))) return
        !  And when  alpha.eq.zero.
        beq0 = .false.
        beq0 = all(beta.v==v8_n0.v)
        if(aeq0) then
           if(upper) then
              if(beq0)  then
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                 do j=1,n
                    !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=1,j
                       c(i,j) = ZERO
                    end do
                 end do
              else
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                 do j=1,n
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=1,j-1
                       c(i,j) = beta*c(i,j)
                    end do
                    c(j,j) = beta.v*c(j,j).re
                 end do
                 !$OMP END PARALLEL DO
              end if
           else
              if(beq0) then
                   !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                 do j=1,n
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=j,n
                       c(i,j) = ZERO
                    end do
                 end do
                 !$OMP END PARALLEL DO
              else
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
                 do j=1,n
                    c(j,j) = beta.v*c(j,j).re
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=j+1,n
                       c(i,j) = beta*c(i,j)
                    end do
                 end do
                 !$OMP END PARALLEL DO
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
               !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
              do j=1,n
                 if(beq0) then
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=1,j
                       c(i,j) = ZERO
                    end do
                 else if(bneq1) then
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
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
                        !$OMP SIMD ALIGNED(c:64,a) LINEAR(i:1) UNROLL PARTIAL(10)
                       do i=1,j-1
                          c(i,j) = c(i,j)+temp*a(i,l)
                       end do
                       c(j,j) = c(j,j).re+temp.re*a(i,l).re
                    end if
                 end do
              end do
              !$OMP END PARALLEL DO
           else
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
              do j=1,n
                 if(beq0) then
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(10)
                    do i=j,n
                       c(i,j) = ZERO
                    end do
                 else if(bneq1) then
                    c(j,j) = beta.v*c(j,j).re
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
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
                        !$OMP SIMD ALIGNED(c:64,a) LINEAR(i:1) UNROLL PARTIAL(10)
                       do i=j+1,n
                          c(i,j) = c(i,j)+temp*a(i,l)
                       end do
                    end if
                 end do
              end do
              !$OMP END PARALLEL DO
            end if
         else
            ! Form  C := alpha*A**H*A + beta*C.
            if(upper) then
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
               do j=1,n
                  do i=1,j-1
                     temp = ZERO
                       !$OMP SIMD ALIGNED(c:64,a) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
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
                   !$OMP SIMD ALIGNED(c:64,a) LINEAR(i:1) REDUCTION(+:rtemp) UNROLL PARTIAL(10)
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
                      !$OMP SIMD ALIGNED(c:64,a) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
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
               !$OMP END PARALLEL DO
            end if
         end if
         !End of ZHERK
     end subroutine zherk

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zhpmv(uplo,n,alpha,ap,x,incx,beta,y,incy) !GCC$ ATTRIBUTES hot :: zhpmv !GCC$ ATTRIBUTES aligned(32) :: zhpmv !GCC$ ATTRIBUTES no_stack_protector :: zhpmv
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zhpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zhpmv
    !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zhpmv

          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none
#endif
         character(len=1),                     intent(in),value :: uplo
         integer(kind=i4),                   intent(in),value :: n
         type(AVX512c8f64_t),                  intent(in)       :: alpha
         !type(AVX512c8f64_t), dimension(*),    intent(in)       :: ap
         type(AVX512c8f64_t), dimension(:), allocatable,    intent(in)       :: ap
         !type(AVX512c8f64_t), dimension(*),    intent(in)       :: x
         type(AVX512c8f64_t), dimension(:), allocatable,    intent(in)       :: x
         integer(kind=int4),                   intent(in),value :: incx
         type(AVX512c8f64_t),                  intent(in)       :: beta
         !type(AVX512c8f64_t), dimension(*),    intent(inout)    :: y
         type(AVX512c8f64_t), dimension(:), allocatable,    intent(in)       :: y
         integer(kind=i4),                   intent(in),value :: incy
         ! Locals
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
         type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
         type(AVX512c8f64_t), automatic :: temp2
         integer(kind=i4),  automatic :: i,info,ix,iy,j,jx,jy,k,kk,kx,ky
         logical(kind=i1),  automatic :: aeq0,beq1,bneq1,beq0
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
         !info = 0
         !if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
         !   info = 1
         !else if(n<0) then
         !   info = 2
        ! else if(incx==0) then
         !   info = 6
        ! else if(incy==0) then
        !    info = 9
        ! end if
        ! if(info/=0) then
        !    call xerbla('GMS_ZHPMV',info)
        !    return
        ! end if
         !aeq0 = .false.
         !beq1 = .false.
        ! aeq0 = all(alpha==ZERO)
        ! beq1 = all(beta==ONE)
         !   Quick return if possible.
        ! if((n==0) .or. (aeq0 .and. beq1)) return
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
                   !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(8)
                   do i=1,n
                      y(i) = ZERO
                   end do
                else
                    !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(6)
                   do i=1,n
                      y(i) = beta*y(i)
                   end do
                end if
             else
                iy = ky
                if(beq0) then
                     !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(8)
                   do i=1,n
                      y(iy) = ZERO
                      iy = iy+incy
                   end do
                else
                     !$OMP SIMD ALIGNED(y:64) LINEAR(i:1) UNROLL PARTIAL(6)
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
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,k,kk) 
                do j=1,n
                   temp1 = alpha*x(j)
                   temp2 = ZERO
                   k = kk
                    !$OMP SIMD ALIGNED(y:64,ap,x) LINEAR(i:1) REDUCTION(+:temp2) UNROLL PARTIAL(10)
                   do i=1,j-1
                      y(i) = y(i)+temp1*ap(k)
                      temp2 = temp2+conjugate(ap(k))*x(i)
                      k = k+1
                   end do
                   y(j) = y(j)+temp1*ap(kk+j-1).re+alpha*temp2
                   kk = kk+j
                end do
                !$OMP END PARALLEL DO
             else
                jx = kx
                jy = ky
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,ix,iy,jx,jy,kk) 
                do j=1,n
                   temp1 = alpha*x(jx)
                   temp2 = ZERO
                   ix = kx
                   iy = ky
                    !$OMP SIMD ALIGNED(y:64,ap,x) LINEAR(i:1) REDUCTION(+:temp2) UNROLL PARTIAL(10)
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
                !$OMP END PARALLEL DO
             end if
          else
             !   Form  y  when AP contains the lower triangle.
             if((incx==1) .and. (incy==1)) then
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,k,kk) 
                do j=1,n
                   temp1 = alpha*x(j)
                   temp2 = ZERO
                   y(j) = y(j)+temp1*ap(kk).re
                   k = kk+1
                    !$OMP SIMD ALIGNED(y:64,ap,x) LINEAR(i:1) REDUCTION(+:temp2) UNROLL PARTIAL(10)
                   do i=j+1,n
                      y(i) = y(i)+temp1*ap(k)
                      temp2 = temp2+conjugate(ap(k))*x(i)
                      k = k+1
                   end do
                   y(j) = y(j)+alpha*temp2
                   kk = kk+(n-j+1)
                end do
                !$OMP END PARALLEL DO
             else
                jx = kx
                jy = ky
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1,temp2,ix,iy,jx,jy,kk) 
                do j=1,n
                   temp1 = alpha*x(jx)
                   temp2 = ZERO
                   y(jy) = y(jy)+temp1*ap(kk).re
                   ix = jx
                   iy = jy
                    !$OMP SIMD ALIGNED(y:64,ap,x) LINEAR(i:1) REDUCTION(+:temp2) UNROLL PARTIAL(10)
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
                !$OMP END PARALLEL DO
             end if
          end if
          !End of ZHPMV
     end subroutine zhpmv

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zhpr(uplo,n,alpha,x,incx,ap) !GCC$ ATTRIBUTES hot :: zhpr !GCC$ ATTRIBUTES aligned(32) :: zhpr !GCC$ ATTRIBUTES no_stack_protector :: zhpr
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zhpr(uplo,n,alpha,x,incx,ap)
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zhpr
   !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zhpr
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none
       

         use mod_vecconsts, only : v8_n0
         character(len=1),                  intent(in),value :: uplo
         integer(kind=i4),                intent(in),value :: n
         type(ZMM8r8_t),                    intent(in)       :: alpha
         !type(AVX512c8f64_t), dimension(*), intent(in)       :: x
          type(AVX512c8f64_t), dimension(:),allocatable,  intent(in)       :: x
         integer(kind=int4),                intent(in),value :: incx
         !type(AVX512c8f64_t), dimension(*), intent(inout)    :: ap
         type(AVX512c8f64_t), dimension(:), allocatable, intent(in)       :: x
         ! Locals
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
         type(AVX512c8f64_t), automatic :: temp
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: vtemp
#endif
         type(ZMM8r8_t),      automatic :: vtemp
         integer(kind=i4),  automatic :: i,info,ix,j,jx,k,kk,kx
         integer(kind=i1),  automatic :: aeq0
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
         type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                                [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp])
         ! Test the input parameters.
         !info = 0
         !if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
         !   info = 1
        ! else if(n<0) then
         !   info = 2
         !else if(incx==0) then
         !   info = 5
         !end if
         if(info/=0) then
         !   call xerbla('GMS_ZHPR',info)
         !   return
         !end if
         !  Quick return if possible.
        ! aeq0 = .false.
        ! aeq0 = all(alpha.v==ZERO.re)
         !if((n==0) .or. (aeq0)) return
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
               !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,vtemp,temp,k,kk) IF(n>=10)
               do j=1,n
                  if(all(x(j)/=ZERO)) then
                     vtemp = conjugate(x(j))
                     temp = zmm8r81x_init(alpha.v*vtemp.v)
                     k = kk
                     !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
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
               !$OMP END PARALLEL DO
            else
               jx = kx
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,vtemp,ix,jx,kk) 
               do j=1,n
                  if(all(x(jx)/=ZERO)) then
                     vtemp = conjugate(x(j))
                     ix = kx
                      !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
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
               !$OMP END PARALLEL DO
            end if
         else
            ! Form  A  when lower triangle is stored in AP.
            if(incx==1) then
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,vtemp,temp,k,kk) 
               do j=1,n
                  if(all(x(j)/=ZERO)) then
                     vtemp = conjugate(x(j))
                     temp = zmm8r81x_init(alpha.v*vtemp.v)
                     ap(kk).re = ap(kk).re+temp.v*x(j).re
                     k = kk+1
                       !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
                     do i=j+1,n
                        ap(k) = ap(k)+x(i)*temp
                        k = k+1
                     end do
                  else
                     ap(kk).re = ap(kk).re
                  end if
                  kk = kk+n-j+1
               end do
               !$OMP END PARALLEL DO
            else
               jx = kx
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,vtemp,temp,ix,jx,kk) 
               do j=1,n
                  if(all(x(jx)/=ZERO)) then
                     vtemp = conjugate(x(j))
                     temp = zmm8r81x_init(alpha.v*vtemp.v)
                     ap(kk).re = app(kk).re+temp.v*x(jx).re
                     ix = jx
                       !$OMP SIMD ALIGNED(a:64,x) LINEAR(i:1) UNROLL PARTIAL(10)
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
               !$OMP END PARALLEL DO
            end if
         end if
         ! End of ZHPR
     end subroutine zhpr

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

#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))       
subroutine zsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc) !GCC$ ATTRIBUTES hot :: zsymm !GCC$ ATTRIBUTES aligned(32) :: zsymm !GCC$ ATTRIBUTES no_stack_protector :: zsymm
#elif defined __INTEL_COMPILER
subroutine zsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zsymm
   !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zsymm
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

         character(len=1),                       intent(in),value :: side
         character(len=1),                       intent(in),value :: uplo
         integer(kind=i4),                     intent(in),value :: m
         integer(kind=i4),                     intent(in),value :: n
         type(AVX512c8f64_t),                    intent(in)       :: alpha
         !type(AVX512c8f64_t), dimension(lda,*),  intent(in)       :: a
          type(AVX512c8f64_t), dimension(:,:),  allocatable, intent(in)       :: a
         integer(kind=i4),                     intent(in),value :: lda
         !type(AVX512c8f64_t), dimension(ldb,*),  intent(in)       :: b
         type(AVX512c8f64_t), dimension(:,:),  allocatable, intent(in)       :: b
         integer(kind=i4),                     intent(in),value :: ldb
         type(AVX512c8f64_t),                    intent(in)       :: beta
         !type(AVX512c8f64_t), dimension(ldc,*),  intent(inout)    :: c
         type(AVX512c8f64_t), dimension(:,:),  allocatable, intent(in)       :: c
         integer(kind=i4),                     intent(in)       :: ldc
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
         !info = 0
         !if((.not.lsame(side,'L')) .and. (.not.lsame(side,'R'))) then
         !   info = 1
         !else if((.not.upper) .and. (.not.lsame(uplo,'L'))) then
         !   info = 2
         !else if(m<0) then
         !   info = 3
         !else if(n<0) then
         !   info = 4
         !else if(lda < max(1,nrowa)) then
         !   info = 7
         !else if(ldb < max(1,m)) then
         !   info = 9
         !else if(ldc < max(1,m)) then
         !   info = 12
         !end if
         !if(info/=0) then
         !   call xerbla('GMS_ZSYMM',info)
         !   return
         !end if
         !aeq0 = .false.
         !beq1 = .false.
         !aeq0 = all(alpha==ZERO)
         !beq1 = all(beta==ONE)
         ! Quick return if possible.
         !if((m==0) .or. (n==0) .or. &
        !    ((aeq0) .and. (beq1))) return
         !  And when  alpha.eq.zero.
         if(aeq0) then
            if(beq0) then
               !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
               do j=1,n
                  !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                  do i=1,n
                     c(i,j) = ZERO
                  end do
               end do
               !$OMP END PARALLEL DO
            else
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
               do j=1,n
                    !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
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
               !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
               do j=1,n
                  do i=1,m
                     temp1 = alpha*b(i,j)
                     temp2 = ZERO
                      !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1) REDUCTION(+:temp2) UNROLL PARTIAL(6)
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
               !$OMP END PARALLEL DO
            else
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j) 
               do j=1,n
                  do i=m,1,-1
                     temp1 = alpha*b(i,j)
                     temp2 = ZERO
                       !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1) REDUCTION(+:temp2) UNROLL PARTIAL(6)
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
               !$OMP END PARALLEL DO
            end if
         else
            !  Form  C := alpha*B*A + beta*C.
             !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp1) 
            do j=1,n
               temp1 = alpha*a(j,j)
               if(beq0) then
                    !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1) UNROLL PARTIAL(6)  
                  do i=1,m
                     c(i,j) = temp1*b(i,j)
                  end do
               else
                   !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1) REDUCTION(+:temp2) UNROLL PARTIAL(10)
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
                    !$OMP SIMD ALIGNED(c:64,b) LINEAR(i:1) UNROLL PARTIAL(6)
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
                    !$OMP SIMD ALIGNED(c:64,b) LINEAR(i:1) REDUCTION(+:temp2) UNROLL PARTIAL(6)
                   do i=1,m
                      c(i,j) = c(i,j)+temp1*b(i,k)
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          end if
          !End of ZSYMM
     end subroutine zsymm

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))        
subroutine zsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc) !GCC$ ATTRIBUTES hot :: zsyr2k !GCC$ ATTRIBUTES aligned(32) :: zsyr2k !GCC$ ATTRIBUTES no_stack_protector :: zsyr2k
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zsyr2k
    !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zsyr2k
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

        character(len=1),                      intent(in),value :: uplo
        character(len=1),                      intent(in),value :: trans
        integer(kind=i4),                    intent(in),value :: n
        integer(kind=i4),                    intent(in),value :: k
        type(AVX512c8f64_t),                   intent(in)       :: alpha
        !type(AVX512c8f64_t), dimension(lda,*), intent(in)       :: a
        type(AVX512c8f64_t), dimension(:,:), allocatable,  intent(in)       :: a
        integer(kind=i4),                    intent(in),value :: lda
        !type(AVX512c8f64_t), dimension(ldb,*), intent(in)       :: b
        type(AVX512c8f64_t), dimension(:,:), allocatable,  intent(in)       :: b
        integer(kind=i4),                    intent(in),value :: ldb
        type(AVX512c8f64_t),                   intent(in)       :: beta
        !type(AVX512c8f64_t), dimension(ldc,*), intent(inout)    :: c
        type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)       :: c 
        integer(kind=i4),                    intent(in)       :: ldc
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp1
#endif
        type(AVX512c8f64_t), automatic :: temp1
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp2
#endif
        type(AVX512c8f64_t), automatic :: temp2
        integer(kind=i4),  automatic :: i,info,j,l,nrowa
        logical(kind=i4),  automatic :: upper
        logical(kind=i1),  automatic :: aeq0,beq1,beq0,bneq1
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
        !info  = 0
        !if((.not.upper) .and. (.not.lsame(uplo,'L'))) then
        !   info = 1
        !e!lse if((.not.lsame(trans,'N')) .and. &
        !     (.not.lsame(trans,'T'))) then
        !   info = 2
        !else if(n<0) then
        !   info = 3
        !else if(k<0) then
        !   info = 4
        !else if(lda < max(1,nrowa)) then
        !   info = 7
        !else if(ldb < max(1,nrowa)) then
        !   info = 9
        !else if(ldc < max(1,n)) then
        !   info = 12
        !end if
        !if(info/=0) then
        !   call xerbla('GMS_ZSYRK2',info)
        !   return
        !end if
        !   Quick return if possible.
        aeq0 = .false.
        beq1 = .false.
        aeq0 = all(alpha==ZERO)
        beq1 = all(beta==ONE)
        !if((n==0) .or. (((aeq0) .or. &
        !     (k==0)) .and. (beq1))) return
        !  And when  alpha.eq.zero.
        if(aeq0) then
           if(upper) then
              if(beq0) then
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEAFULT(SHARED) PRIVATE(j)
                 do j=1,n
                    !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=1,j
                       c(i,j) = ZERO
                    end do
                 end do
                 !$OMP END PARALLEL DO
              else
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEAFULT(SHARED) PRIVATE(j)
                 do j=1,n
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=1,j
                       c(i,j) = beta*c(i,j)
                    end do
                 end do
              end if
           else
              if(beq0) then
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEAFULT(SHARED) PRIVATE(j) 
                 do j=1,n
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=j,n
                       c(i,j) = ZERO
                    end do
                 end do
                 !$OMP END PARALLEL DO
              else
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEAFULT(SHARED) PRIVATE(j)
                 do j=1,n
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=j,n
                       c(i,j) = beta*c(i,j)
                    end do
                 end do
                 !$OMP END PARALLEL DO
              end if
           end if
           return
        end if
        !  Start the operations.
        if(upper) then
           !$OMP PARALLEL DO SCHEDULE(STATIC) DEAFULT(SHARED) PRIVATE(j) 
           do j=1,n
              if(beq0) then
                   !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8) 
                 do i=1,j
                    c(i,j) = ZERO
                 end do
              else if(bneq1) then
                    !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                 do i=1,j
                    c(i,j) = beta*c(i,j)
                 end do
              end if
              do l=1,k
                 if(all(a(j,l)/=ZERO) .or. all(b(j,l)/=ZERO)) then
                    temp1 = alpha*b(j,l)
                    temp2 = alpha*a(j,l)
                     !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=1,j
                       c(i,j) = c(i,j)+a(i,l)*temp1 + &
                       b(i,l)*temp2
                    end do
                 end if
              end do
           end do
           !$OMP END PARALLEL DO
        else
            !$OMP PARALLEL DO SCHEDULE(STATIC) DEAFULT(SHARED) PRIVATE(j)  
           do j=1,n
              if(beq0) then
                  !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)   
                 do i=j,n
                    c(i,j) = ZERO
                 end do
              else if(bneq1) then
                 !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                 do i=j,n
                    c(i,j) = beta*c(i,j)
                 end do
              end if
              do l=1,k
                 if(all(a(j,l)/=ZERO) .or. all(b(j,l)/=ZERO)) then
                     !$OMP SIMD ALIGNED(c:64,a,b) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=j,n
                       c(i,j) = c(i,j)+a(i,l)*temp1 + &
                            b(i,l)*temp2
                    end do
                 end if
              end do
           end do
           !$OMP END PARALLEL DO
        end if
     else
        !  Form  C := alpha*A**T*B + alpha*B**T*A + C.
        if(upper) then
           !$OMP PARALLEL DO SCHEDULE(STATIC) DEAFULT(SHARED) PRIVATE(j)  
           do j=1,n
              do i=1,j
                 temp1 = ZERO
                 temp2 = ZERO
                   !$OMP SIMD ALIGNED(a:64,b) LINEAR(i:1) REDUCTION(+:temp1,temp2) UNROLL PARTIAL(10) 
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
           !$OMP END PARALLEL DO
        else
            !$OMP PARALLEL DO SCHEDULE(STATIC) DEAFULT(SHARED) PRIVATE(j)  
           do j=1,n
              do i=j,n
                 temp1 = ZERO
                 temp = ZERO
                   !$OMP SIMD ALIGNED(a:64,b) LINEAR(i:1) REDUCTION(+:temp1,temp2) UNROLL PARTIAL(10) 
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
           !$OMP END PARALLEL DO
        end if
     end if
      ! End of ZSYR2K
   end subroutine zsyr2k

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine zsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc) !GCC$ ATTRIBUTES hot :: zsyrk !GCC$ ATTRIBUTES aligned(32) :: zsyrk !GCC$ ATTRIBUTES no_stack_protector :: zsyrk
#if defined(__INTEL_COMPILER) || defined(__ICC)
subroutine zsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)  
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: zsyrk
   !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: zsyrk
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

        character(len=1),                      intent(in),value  :: uplo
        character(len=1),                      intent(in),value  :: trans
        integer(kind=i4),                    intent(in),value  :: n
        integer(kind=i4),                    intent(in),value  :: k
        type(AVX512c8f64_t),                   intent(in)        :: alpha
        !type(AVX512c8f64_t), dimension(lda,*), intent(in)        :: a
         type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)        :: a
        integer(kind=i4),                    intent(in),value  :: lda
        type(AVX512c8f64_t),                   intent(in)        :: beta
        !type(AVX512c8f64_t), dimension(ldc,*), intent(inout)     :: c
         type(AVX512c8f64_t), dimension(:,:), allocatable, intent(in)        :: c
        integer(kind=i4),                    intent(in),value  :: ldc
        ! LOcals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
        integer(kind=i4),  automatic :: i,info,j,l,nrowa
        logical(kind=i4),  automatic :: upper
        logical(kind=i1),  automatic :: aeq0,beq1,beq0,bneq1
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
        !info = 0
        !if((.not.upper) .and. (.not.lsame(uplo,'L'))) then
        !   info = 1
        !else if((.not.lsame(trans,'N')) .and. &
       !      (.not.lsame(trans,'T'))) then
        !   info = 2
        !else if(n<0) then
        !   info = 3
        !else if(k<0) then
        !   info = 4
        !else if(lda < max(1,nrowa)) then
        !   info = 7
        !else if(ldc < max(1,n)) then
        !   info = 10
        !end if
        !if(info/=0) then
        !   call xerbla('GMS_ZSYRK',info)
        !   return
        !end do
        ! Quick return if possible.
        aeq0 = .false.
        beq1 = .false.
        aeq0 = all(alpha==ZERO)
        beq1 = all(beta==ONE)
        !if((n==0) .or. (((aeq0) .or. &
        !     (k==0)) .and.(beq1))) return
        !  And when  alpha.eq.zero.
        if(aeq0) then
           if(upper) then
              if(beq0) then
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                 do j=1,n
                    !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=1,j
                       c(i,j) = ZERO
                    end do
                 end do
                 !$OMP END PARALLEL DO
              else
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                 do j=1,n
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=1,j
                       c(i,j) = beta*c(i,j)
                    end do
                 end do
              end if
           else
              if(beq0) then
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                 do j=1,n
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=j,n
                       c(i,j) = ZERO
                    end do
                 end do
              else
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
                 do j=1,n
                      !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=j,n
                       c(i,j) = beta*c(i,j)
                    end do
                 end do
                 !$OMP END PARALLEL DO
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
               !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
              do j=1,n
                 if(beq0) then
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=1,j
                       c(i,j) = ZERO
                    end do
                 else if(bneq1) then
                     !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=1,j
                       c(i,j) = beta*c(i,j)
                    end do
                 end if
                 do l=1,k
                    if(all(a(j,l)/=ZERO)) then
                       temp = alpha*a(j,l)
                        !$OMP SIMD ALIGNED(c:64,a) LINEAR(i:1) UNROLL PARTIAL(10)
                       do i=1,j
                          c(i,j) = c(i,j)+temp*a(i,l)
                       end do
                    end if
                 end do
              end do
              !$OMP END PARALLEL DO
           else
                !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j)
              do j=1,n
                 if(beq0) then
                      !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(8)
                    do i=j,n
                       c(i,j) = ZERO
                    end do
                 else if(bneq1) then
                       !$OMP SIMD ALIGNED(c:64) LINEAR(i:1) UNROLL PARTIAL(6)
                    do i=j,n  
                      c(i,j) = beta*c(i,j)
                   end do
                end if
                do l=1,k
                   if(all(a(j,l)/=ZERO)) then
                      temp = alpha*a(j,l)
                        !$OMP SIMD ALIGNED(c:64,a) LINEAR(i:1) UNROLL PARTIAL(10)
                      do i=j,n
                         c(i,j) = c(i,j)+temp*a(i,l)
                      end do
                   end if
                end do
             end do
             !$OMP END PARALLEL DO
          end if
       else
          !   Form  C := alpha*A**T*A + beta*C.
          if(upper) then
             !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,i) COLLAPSE(2)
             do j=1,n
                do i=1,j
                   temp = ZERO
                    !$OMP SIMD ALIGNED(c:64,a) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
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
             !$OMP END PARALLEL DO
          else
             !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,i) COLLAPSE(2)
             do j=1,n
                do i=1,n
                   temp = ZERO
                    !$OMP SIMD ALIGNED(c:64,a) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
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
             !$OMP END PARALLEL DO
          end if
       end if
       ! End of ZSYRK
     end subroutine zsyrk

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine ztbsv(uplo,trans,diag,n,k,a,lda,x,incx) !GCC$ ATTRIBUTES hot :: ztbsv !GCC$ ATTRIBUTES aligned(32) :: ztbsv !GCC$ ATTRIBUTES no_stack_protector :: ztbsv
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine ztbsv(uplo,trans,diag,n,k,a,lda,x,incx)  
       !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ztbsv
        !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ztbsv
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none
        character(len=1),                         intent(in),value  :: uplo
        character(len=1),                         intent(in),value  :: trans
        character(len=1),                         intent(in),value  :: diag
        integer(kind=i4),                       intent(in),value  :: n
        integer(kind=i4),                       intent(in),value  :: k
        !type(AVX512c8f64_t), dimension(lda,*),    intent(in)        :: a
        type(AVX512c8f64_t), dimension(:,:), allocatable,    intent(in)        :: a
        integer(kind=i4),                       intent(in),value  :: lda
        !type(AVX512c8f64_t), dimension(*),        intent(inout)     :: x
        type(AVX512c8f64_t), dimension(:), allocatable,    intent(in)        :: a
        integer(kind=i4),                       intent(in),value  :: incx
        ! Locals
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
        type(AVX512c8f64_t), automatic :: temp
        integer(kind=i4),  automatic :: i,info,ix,j,jx,kplus1,kx,l
        logical(kind=i4),  automatic :: noconj,nounit
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        ! Test the input parameters.
        !info = 0
        !if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
         !  info = 1
        !else if(.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. &
        !     .not.lsame(trans,'C')) then
        !   info = 2
        !else if(.not.lsame(diag,'U') .and. .not.lsame(diag,'N')) then
        !   info = 3
        !else if(n<0) then
        !   info = 4
        !else if(k<0) then
        !   info = 5
        !else if(lda<(k+1)) then
        !   info = 7
        !else if(incx==0) then
        !   info = 9
        !end if
        !if(info/=0) then
         !  call xerbla('GMS_ZTBSV',info)
        !   return
        !end if
        !  Quick return if possible.
        !if(n==0) return
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
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,l,temp) 
                 do j=n,1,-1
                    if(all(x(j)/=ZERO)) then
                       l = kplus1-j
                       if(nounit) x(j) = x(j)/a(kplus1,j)
                       temp = x(j)
                       !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) UNROLL PARTIAL(6)
                       do i=j-1,max(1,j-k),-1
                          x(i) = x(i)-temp*a(l+i,j)
                       end do
                    end if
                 end do
                 !$OMP END PARALLEL DO
              else
                 kx = kx+(n-1)*incx
                 jx = kx
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,kx,ix,jx,l,temp) 
                 do j=n,1,-1
                    kx = kx-incx
                    if(all(x(jx)/=ZERO)) then
                       ix = kx
                       l = kplus1-j
                       if(nounit) x(jx) = x(jx)/a(kplus1,j)
                       temp = x(jx)
                         !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) UNROLL PARTIAL(6)
                       do i=j-1,max(1,j-k),-1
                          x(ix) = x(ix)-temp*a(l+i,j)
                       end do
                    end if
                    jx = kx-incx
                 end do
                 !$OMP END PARALLEL DO
              end if
           else
              if(incx==1) then
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,l,temp) 
                 do j=1,n
                    if(all(x(j)/=ZERO)) then
                       l = 1-j
                       if(nounit) x(j) = x(j)/a(1,j)
                       temp = x(j)
                         !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) UNROLL PARTIAL(6)
                       do i=j+1,min(n,j+k)
                          x(i) = x(i)-temp*a(l+i,j)
                       end do
                    end if
                 end do
                 !$OMP END PARALLEL DO
              else
                 jx = kx
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,kx,ix,jx,l,temp) 
                 do j=1,n
                    kx = kx+incx
                    if(all(x(jx)/=ZERO)) then
                       ix = kx
                       l = 1-j
                       if(nounit) x(jx) = x(jx)/a(1,j)
                       temp = x(jx)
                         !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) UNROLL PARTIAL(6)
                       do i=j+1,min(n,j+k)
                          x(ix) = x(ix)-temp*a(l+i,j)
                          ix = ix+incx
                       end do
                    end if
                    jx = jx+incx
                 end do
                 !$OMP END PARALLEL DO
              end if
           end if
        else
           !  Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
           if(lsame(uplo,'U')) then
              kplus1 = k+1
              if(incx==1) then
                   !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,l,temp) 
                 do j=1,n
                    temp = x(j)
                    l = kplus1-j
                    if(noconj) then
                       !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) REDUCTION(-:temp) UNROLL PARTIAL(6)
                       do i=max(1,j-k),j-1
                          temp = temp-a(l+i,j)*x(i)
                       end do
                       if(nounit) temp = temp/a(kplus1,j)
                    else
                         !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) REDUCTION(-:temp) UNROLL PARTIAL(6)
                       do i=max(1,j-k),j-1
                          temp = temp-conjugate(a(l+i,j))*x(i)
                       end do
                       if(nounit) temp = temp/conjugate(a(kplus1,j))
                    end if
                    x(j) = temp
                 end do
                 !$OMP END PARALLEL DO
              else
                 jx = kx
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,l,ix,temp,jx,kx) 
                 do j=1,n
                    temp = x(jx)
                    ix = kx
                    l = kplus1-j
                          !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) REDUCTION(-:temp) UNROLL PARTIAL(6)
                       do i=max(1,j-k),j-1
                          temp = temp-a(l+i,j)*x(ix)
                          ix = ix+incx
                       end do
                       if(nounit) temp = temp/a(kplus1,j)
                    else
                         !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) REDUCTION(-:temp) UNROLL PARTIAL(6)
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
                 !$OMP END PARALLEL DO
              end if
           else
              if(incx==1) then
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,l,temp) 
                 do j=n,1,-1
                    temp = x(j)
                    l =1-j
                    if(noconj) then
                        !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) REDUCTION(-:temp) UNROLL PARTIAL(6)
                       do i=min(n,j+k),j+1,-1
                          temp = temp-a(l+i,j)*x(i)
                       end do
                       if(nounit) temp = temp/a(1,j)
                    else
                         !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) REDUCTION(-:temp) UNROLL PARTIAL(6)
                       do i=min(n,j+k),j+1,-1
                          temp = temp-conjugate(a(l+i,j))*x(i)
                       end do
                       if(nounit) temp = temp/conjugate(a(1,j))
                    end if
                    x(j) = temp
                 end do
                 !$OMP END PARALLEL DO
              else
                 kx = kx+(n-1)*incx
                 jx = kx
                     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,l,ix,temp,jx,kx) 
                 do j=n,1,-1
                    temp = x(jx)
                    ix = kx
                    l = 1-j
                    if(noconj) then
                         !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) REDUCTION(-:temp) UNROLL PARTIAL(6)
                       do i=min(n,j+k),j+1,-1
                          temp = temp-a(l+i,j)*x(i)
                          ix = ix+incx
                       end do
                       if(nounit) temp = temp/a(1,j)
                    else
                          !$OMP SIMD ALIGNED(x:64,a) LINEAR(i:1) REDUCTION(-:temp) UNROLL PARTIAL(6)
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
                 !$OMP END PARALLEL DO
              end if
           end if
        end if
        !End of ZTBSV
     end subroutine ztbsv

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine ztpmv(uplo,trans,diag,n,ap,x,incx) !GCC$ ATTRIBUTES Hot :: ztpmv !GCC$ ATTRIBUTES aligned(32) :: ztpmv !GCC$ ATTRIBUTES no_stack_protector :: ztpmv
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine ztpmv(uplo,trans,diag,n,ap,x,incx)
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ztpmv
      !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ztpmv
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

   
         character(len=1),                         intent(in),value :: uplo
         character(len=1),                         intent(in),value :: trans
         character(len=1),                         intent(in),value :: diag
         integer(kind=i4),                       intent(in),value :: n
         !type(AVX512c8f64_t), dimension(*),        intent(in)       :: ap
         type(AVX512c8f64_t), dimension(:), allocatable,       intent(in)       :: ap
         !type(AVX512c8f64_t), dimension(*),        intent(inout)    :: x
         type(AVX512c8f64_t), dimension(:), allocatable,        intent(in)       :: x
         integer(kind=i4),                       intent(in),value :: incx
         ! LOcals
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
         type(AVX512c8f64_t), automatic :: temp
         integer(kind=i4),  automatic :: i,info,ix,j,jx,k,kk,kx
         logical(kind=i4),  automatic :: noconj,nounit
#if defined __INTEL_COMPILER
        !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
        type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        ! Exec code ....
        !info = 0
       ! if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
        !   info = 1
        !else if(.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. &
        !     .not.lsame(trans,'C')) then
        !   info = 2
        !else if(.not.lsame(diag,'U') .and. .not.lsame(diag,'N')) then
        !   info = 3
        !else if(n<0) then
        !   info = 4
        !else if(incx==0) then
        !   info = 7
        !end if
        !if(info/=0) then
        !   call xerbla('GMS_ZTPMV',info)
        !   return
        !end if
        !  Quick return if possible.
        !if(n==0) return
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
                 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,k,kk)
                 do j=1,n
                    if(all(x(j)/=ZERO)) then
                       temp = x(j)
                       k = kk
                        !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) UNROLL PARTIAL(10)
                       do i=1,j-1
                          x(i) = x(i)+temp*ap(k)
                          k = k+1
                       end do
                       if(nounit) x(j) = x(j)*ap(kk+j-1)
                    end if
                    kk = kk+j
                 end do
                 !$OMP END PARALLEL DO
              else
                 jx = kx
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jx,kk)
                 do j=1,n
                    if(all(x(jx)/=ZERO)) then
                       temp = x(jx)
                       ix = kx
                         !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) UNROLL PARTIAL(10)
                       do k=kk,kk+j-2
                          x(ix) = x(ix)+temp*ap(k)
                          ix = ix+incx
                       end do
                       if(nounit) x(jx) = x(jx)*ap(kk+j-1)
                    end if
                    jx = jx+incx
                    kk = kk+j
                 end do
                 !$OMP END PARALLEL DO
              end if
           else
              kk = (n*(n+1))/2
              if(incx==1) then
                  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,k,kk)
                 do j=n,1,-1
                    if(all(x(j)/=ZERO)) then
                       temp = x(j)
                       k = kk
                        !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) UNROLL PARTIAL(10)
                       do i=n,j+1,-1
                          x(i) = x(i)+temp*ap(k)
                          k = k-1
                       end do
                       if(nounity) x(j) = x(j)*ap(kk-n+j)
                    end if
                    kk = kk-(n-j+1)
                 end do
                 !$OMP END PARALLEL DO
              else
                 kx = kx+(n-1)*incx
                 jx = kx
                    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jx,kk)
                 do j=n,1,-1
                    if(all(x(jx)/=ZERO)) then
                       temp = x(jx)
                       ix = kx
                        !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) UNROLL PARTIAL(10)
                       do k=kk,kk-(n-(j+1)),-1
                          x(ix) = x(ix)+temp*ap(k)
                          ix = ix+incx
                       end do
                       if(nounit) x(jx) = x(jx)*ap(kk-n+j)
                    end if
                    jx = jx-incx
                    kk = kk-(n-j+1)
                 end do
                 !$OMP END PARALLEL DO
              end if
            end if
         else
            !   Form  x := A**T*x  or  x := A**H*x.
            if(lsame(uplo,'U')) then
               kk = (n*(n+1))/2
               if(incx==1) then
                   !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,k,kk)
                  do j=n,1,-1
                     temp = x(j)
                     k = kk-1
                     if(noconj) then
                        if(nounit) temp = temp*ap(kk)
                         !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                        do i=j-1,1,-1
                           temp = temp+ap(k)*x(i)
                           k = k-1
                        end do
                     else
                        if(nounit) temp = temp*conjugate(ap(kk))
                        !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                        do i=j-1,1,-1
                           temp = temp+conjugate(ap(k))*x(i)
                           k = k-1
                        end do
                     end if
                     x(j) = temp
                     kk = kk-j
                  end do
                  !$OMP END PARALLEL DO
               else
                  jx = kx+(n-1)*incx
                    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jx,kk)
                  do j=n,1,-1
                     temp = x(jx)
                     ix = jx
                     if(noconj) then
                        if(nounit) temp = temp*ap(kk)
                         !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                        do k=kk-1,kk-j+1,-1
                           ix = ix-incx
                           temp = temp+ap(k)*x(ix)
                        end do
                     else
                        if(nounit) temp = temp*conjugate(ap(kk))
                        !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                        do k=kk-1,kk-j+1,-1
                           ix = ix-incx
                           temp = temp+conjugate(ap(k))*x(ix)
                        end do
                     end if
                     x(jx) = temp
                     jx = jx-incx
                     kk = kk-j
                  end do
                  !$OMP END PARALLEL DO
               end if
            else
               kk = 1
               if(incx==1) then
                    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,k,kk)
                  do j=1,n
                     temp = x(j)
                     k = kk+1
                     if(noconj) then
                        if(nounit) temp = temp*ap(kk)
                        !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                        do i=j+1,n
                           temp = temp+ap(k)*x(i)
                           k = k+1
                        end do
                     else
                        if(nounit) temp = temp+conjugate(ap(kk))
                         !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                        do i=j+1,n
                           temp = temp+conjugate(ap(k))*x(i)
                           k = k+1
                        end do
                     end if
                     x(j) = temp
                     kk = kk+(n-j+1)
                  end do
                  !$OMP END PARALLEL DO
               else
                  jx = kx
                    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jx,kk)
                  do j=1,n
                     temp = x(jx)
                     ix = jx
                     if(noconj) then
                        if(nounit) temp = temp*ap(kk)
                          !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                        do k=kk+1,kk+n-j
                           ix = ix+incx
                           temp = temp+ap(k)*x(ix)
                        end do
                     else
                        if(nounit) temp = temp*conjugate(ap(k))
                          !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(+:temp) UNROLL PARTIAL(10)
                       do k=kk+1,kk+n-j
                          ix = ix+incx
                          temp = temp+conjugate(ap(k))*x(ix)
                       end do
                    end if
                    x(jx) = temp
                    jx = jx+incx
                    kk = kk+(n-j+1)
                 end do
                 !$OMP END PARALLEL DO
              end if
            end if
         end if
         ! End of ZTPMV
      end subroutine ztpmv

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
#if defined(__GFORTRAN__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))  
subroutine ztpsv(uplo,trans,diag,n,ap,x,incx) !GCC$ ATTRIBUTES hot :: ztpsv !GCC$ ATTRIBUTES aligned(32) :: ztpsv !GCC$ ATTRIBUTES no_stack_protector :: ztpsv
#elif defined(__INTEL_COMPILER) || defined(__ICC)
subroutine ztpsv(uplo,trans,diag,n,ap,x,incx)
  !DIR$ ATTRIBUTES CODE_ALIGN : 32 :: ztpsv
     !DIR$ OPTIMIZE : 3
  !DIR$ ATTRIBUTES OPTIMIZATION_PARAMETER: TARGET_ARCH=skylake_avx512 :: ztpmv
#endif
          use module_kinds, only : i4, i1
          use mod_avx512c8f64
          use omp_lib
          implicit none

           character(len=1),                    intent(in),value :: uplo
           character(len=1),                    intent(in),value :: trans
           character(len=1),                    intent(in),value :: diag
           integer(kind=i4),                  intent(in),value :: n
           !type(AVX512c8f64_t), dimension(*),   intent(in)       :: ap
           type(AVX512c8f64_t), dimension(:), allocatable,   intent(in)       :: ap
           !type(AVX512c8f64_t), dimension(*),   intent(inout)    :: x
           type(AVX512c8f64_t), dimension(:), allocatable,  intent(in)       :: x
           integer(kind=i4),                  intent(in),value :: incx
           ! LOcals
#if defined __INTEL_COMPILER
         !DIR$ ATTRIBUTES ALIGN : 64 :: temp
#endif
           type(AVX512c8f64_t), automatic :: temp
           integer(kind=i4),  automatic :: i,info,ix,j,jx,k,kk,kx
           logical(kind=i4),  automatic :: noconj,nounit           
         
#if defined __INTEL_COMPILER
          !DIR$ ATTRIBUTES ALIGN : 64 :: ZERO
#endif
           type(AVX512c8f64_t), parameter :: ZERO = AVX512c8f64_t([0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                                0.0_dp,0.0_dp,0.0_dp,0.0_dp],&
                                                               [0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                                               0.0_dp,0.0_dp,0.0_dp,0.0_dp])
           ! Exec code ....
          ! info = 0
           !if(.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then
          !    info = 1
          ! else if(.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. &
          !    .not.lsame(trans,'C')) then
          !    info = 2
          ! else if(.not.lsame(diag,'U') .and. .not.lsame(diag,'N')) then
          !    info = 3
          ! else if(n<0) then
          !    info = 4
          ! else if(incx==0) then
          !    info = 7
          ! end if
          ! if(info/=0) then
          !    call xerbla('GMS_ZTPSV',info)
          !    return
          ! end if
           ! Quick return if possible.
           !if(n==0) return
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
                    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,k,kk)
                    do j=n,1,-1
                       if(all(x(j)/=ZERO)) then
                          if(nounit) x(j) = x(j)/ap(kk)
                          temp = x(j)
                          k = kk-1
                          !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) UNROLL(6)
                          do i=j-1,1,-1
                             x(i) = x(i)-temp*ap(k)
                             k = k-1
                          end do
                       end if
                       kk = kk-j
                    end do
                    !$OMP END PARALLEL DO
                 else
                    jx = kx+(n-1)*incx
                     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jx,kk)
                    do j=n,1,-1
                       if(all(x(jx)/=ZERO)) then
                          if(nounit) x(jx) = x(jx)/ap(kk)
                          temp = x(jx)
                          ix = jx
                            !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) UNROLL(6)
                          do k=kk-1,kk-j+1,-1
                             ix = ix-incx
                             x(ix) = x(ix)-temp*ap(k)
                          end do
                       end if
                       jx = jx-incx
                       kk = kk-j
                    end do
                    !$OMP END PARALLEL DO
                 end if
              else
                 kk = 1
                 if(incx==1) then
                     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,k,kk)
                    do j=1,n
                       if(all(x(j)/=ZERO)) then
                          if(nounit) x(j) = x(j)/ap(kk)
                          temp = x(j)
                          k = kk+1
                            !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) UNROLL(6)
                          do i=j+1,n
                             x(i) = x(i)-temp*ap(k)
                             k = k+1
                          end do
                       end if
                       kk = kk+(n-j+1)
                    end do
                    !$OMP END PARALLEL DO
                 else
                    jx = kx
                     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jx,kk)
                    do j=1,n
                       if(all((jx)/=ZERO)) then
                          if(nounit) x(jx) = x(jx)/ap(kk)
                          temp = x(jx)
                          ix = jx
                           !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) UNROLL(6)
                          do k=kk+1,kk+n-j
                             ix = ix+incx
                             x(ix) = x(ix)-temp*ap(k)
                          end do
                       end if
                       jx = jx+incx
                       kk = kk+(n-j+1)
                    end do
                    !$OMP END PARALLEL DO
                 end if
              end if
           else
              !  Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
              if(lsame(uplo,'U')) then
                 kk = 1
                 if(incx==1) then
                      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,k,kk)
                    do j=1,n
                       temp = x(j)
                       k = kk
                       if(noconj) then
                            !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) UNROLL(6)
                          do i=1,j-1
                             temp = temp-ap(k)*x(i)
                             k = k+1
                          end do
                          if(nounit) temp = temp/ap(kk+j-1)
                       else
                            !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) UNROLL(6)
                          do i=1,j-1
                             temp = temp-conjugate(ap(k))*x(i)
                             k = k+1
                          end do
                          if(nounit) temp = temp/conjugate(ap(kk+j-1))
                       end if
                       x(j) = temp
                       kk = kk+j
                    end do
                    !$OMP END PARALLEL DO
                 else
                    jx = kx
                      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jx,kk)
                    do j=1,n
                       temp = x(jx)
                       ix = kx
                       if(noconj) then
                           !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(-:temp) UNROLL(6)
                          do k=kk,kk+j-2
                             temp = temp-ap(k)*x(ix)
                             ix = ix+incx
                          end do
                          if(nounit) temp = temp/ap(kk+j-1)
                       else
                           !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(-:temp) UNROLL(6)
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
                    !$OMP END PARALLEL DO
                 end if
              else
                 kk = (n*(n+1))/2
                 if(incx==1) then
                     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,k,kk)
                    do j=n,1,-1
                       temp = x(j)
                       k = kk
                       if(noconj) then
                           !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(-:temp) UNROLL(6)
                          do i=n,j+1,-1
                             temp = temp-ap(k)*x(i)
                             k = k-1
                          end do
                          if(nounit) temp = temp/ap(kk-n+j)
                       else
                           !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(-:temp) UNROLL(6)
                          do i=n,j+1,-1
                             temp = temp/conjugate(ap(kk-n+j))
                             k = k-1
                          end do
                          if(nounit) temp = temp/conjugate(ap(kk-n+j))
                       end if
                       x(j) = temp
                       kk = kk-(n-j+1)
                    end do
                    !$OMP END PARALLEL DO
                 else
                    kx = kx+(n-1)*incx
                    jx = kx
                       !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j,temp,ix,jx,kk)
                    do j=n,1,-1
                       temp = x(jx)
                       ix = kx
                       if(noconj) then
                           !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(-:temp) UNROLL(6)
                          do k=kk,kk-(n-(j+1)),-1
                             temp = temp-ap(k)*x(ix)
                             ix = ix-incx
                          end do
                          if(nounit) temp = temp/ap(kk-n+j)
                       else
                           !$OMP SIMD ALIGNED(x:64,ap) LINEAR(i:1) REDUCTION(-:temp) UNROLL(6)
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
                    !$OMP END PARALLEL DO
                 end if
              end if
           end if
           ! End of ZTPSV
      end subroutine ztpsv
        
    

