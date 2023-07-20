



!/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!*/

module rcs_cylinder_zmm16r4


!===================================================================================85
 !---------------------------- DESCRIPTION ------------------------------------------85
 !
 !
 !
 !          Module  name:
 !                         rcs_cylinder_zmm16r4
 !          
 !          Purpose:
 !                        Various characteristics of analytically derived Radar
 !                        Cross Section of cylindrical objects  
 !                        Based  on George T. Ruck, Donald E. Barrick , William D. Stuart , 
 !                        - "Radar Cross Section Handbook 1 and 2" (1970, Kluwer Academic Plenum Publishers) 
 !                        This module contains only explicitly vectorized (SIMD)
 !                        
 !          History:
 !                        Date: 08-07-2023
 !                        Time: 14:21 GMT+2
 !                        
 !          Version:
 !
 !                      Major: 1
 !                      Minor: 0
 !                      Micro: 0
 !
 !          Author:  
 !                      Bernard Gingold
 !          
 !                 
 !          References:
 !         
 !                      George T. Ruck, Donald E. Barrick , William D. Stuart
 !                      Radar Cross Section Handbook 1 and 2" (1970, Kluwer Academic Plenum Publishers)     
 !         
 !          E-mail:
 !                  
 !                      beniekg@gmail.com
!==================================================================================85
    ! Tab:5 col - Type and etc.. definitions
    ! Tab:10,11 col - Type , function and subroutine code blocks.

    use mod_kinds,    only : i4,sp
    use mod_vectypes, only : ZMM16r4_t
    use avx512_cvec16_v2
    
    public
    implicit none
    
     ! Major version
     integer(kind=i4),  parameter :: RCS_CYLINDER_ZMM16R4_MAJOR = 1
     ! Minor version
     integer(kind=i4),  parameter :: RCS_CYLINDER_ZMM16R4_MINOR = 0
     ! Micro version
     integer(kind=i4),  parameter :: RCS_CYLINDER_ZMM16R4_MICRO = 0
     ! Full version
     integer(kind=i4),  parameter :: RCS_CYLINDER_ZMM16R4_FULLVER =   &
            1000*RCS_CYLINDER_ZMM16R4_MAJOR+100*RCS_CYLINDER_ZMM16R4_MINOR+10*RCS_CYLINDER_ZMM16R4_MICRO
     ! Module creation date
     character(*),        parameter :: RCS_CYLINDER_ZMM16R4_CREATE_DATE = "18-07-2022 14:30 +00200 (TUE 18 JUL 2023 GMT+2)"
     ! Module build date
     character(*),        parameter :: RCS_CYLINDER_ZMM16R4_BUILD_DATE  = __DATE__ " " __TIME__
     ! Module author info
     character(*),        parameter :: RCS_CYLINDER_ZMM16R4_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com"
     ! Short description
     character(*),        parameter :: RCS_CYLINDER_ZMM16R4_SYNOPSIS    = "Analytical Cylindrical objects RCS characteristics and models explicitly vectorized (SIMD)."
    
#ifndef __RCS_CYLINDER_PF_CACHE_HINT__
#define __RCS_CYLINDER_PF_CACHE_HINT__ 1
#endif 
    
     contains
     
     
                   !/* 
                   !      Low frequency scattering widths (k0a << 1).
                   !      Backscatter scattering width for E-field 
                   !      cylinder-parallel,formula 4.1-19
                   ! */
                   
                   
             pure function rcs_f419_zmm16r4(a,k0a) result(rcs)
                  
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f419_zmm16r4
                   !dir$ attributes forceinline :: rcs_f419_zmm16r4
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f419_zmm16r4
                   type(ZMM16r4_t),  intent(in) :: a
                   type(ZMM16r4_t),  intent(in) :: k0a
                   type(ZMM16r4_t) :: rcs
                   ! locals
                   type(ZMM16r4_t), parameter :: C9869604401089358618834490999876 =             &
                                                 ZMM16r4_t(9.869604401089358618834490999876_sp)
                   type(ZMM16r4_t), parameter :: C2467401100272339654708622749969 =             &
                                                 ZMM16r4_t(2.467401100272339654708622749969_sp)
                   type(ZMM16r4_t), parameter :: C08905  = ZMM16r4_t(0.8905_sp)
                   ZMM16r4_t, automatic :: num,arg,ln,ln2,den
                   num.v = a.v*C9869604401089358618834490999876.v
                   arg.v = k0a.v*C08905.v
                   ln.v  = log(arg.v)
                   ln2.v = ln.v*ln.v
                   den.v = k0a.v*ln2.v+C2467401100272339654708622749969.v
                   rcs.v = num.v/den.v
             end function rcs_f419_zmm16r4
             
             
             subroutine rcs_f419_zmm16r4_unroll16x(pa,pk0a,prcs,n,PF_DIST)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f419_zmm16r4_unroll16x
                   !dir$ attributes forceinline :: rcs_f419_zmm16r4_unroll16x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f419_zmm16r4_unroll16x
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out)   :: prcs 
                   integer(kind=i4)               , intent(in)    :: n
                   integer(kind=i4)               , intent(in)    :: PF_DIST
                   ! Locals
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,16)
                   if(m/=0) then
                      do i=1,m
                         a0.v       = pa(i).v
                         k0a0.v     = pk0a(i).v
                         rcs0       = rcs_f419_zmm16r4(a0,k0a0)
                         prcs(i).v  = rcs0.v
                      end do
                      if(n<16) return
                   end if
                   m1=m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                   do i=m1,n,16
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif       
                        a0.v         = pa(i+0).v
                        k0a0.v       = pk0a(i+0).v
                        rcs0         = rcs_f419_zmm16r4(a0,k0a0)
                        prcs(i+0).v  = rcs0.v  
                        a1.v         = pa(i+1).v
                        k0a1.v       = pk0a(i+1).v
                        rcs1         = rcs_f419_zmm16r4(a1,k0a1)
                        prcs(i+1).v  = rcs1.v 
                        a2.v         = pa(i+2).v
                        k0a2.v       = pk0a(i+2).v
                        rcs2         = rcs_f419_zmm16r4(a2,k0a2)
                        prcs(i+2).v  = rcs2.v  
                        a3.v         = pa(i+3).v
                        k0a3.v       = pk0a(i+3).v
                        rcs3         = rcs_f419_zmm16r4(a3,k0a3)
                        prcs(i+3).v  = rcs3.v  
                        a0.v         = pa(i+4).v
                        k0a0.v       = pk0a(i+4).v
                        rcs0         = rcs_f419_zmm16r4(a0,k0a0)
                        prcs(i+4).v  = rcs0.v  
                        a1.v         = pa(i+5).v
                        k0a1.v       = pk0a(i+5).v
                        rcs1         = rcs_f419_zmm16r4(a1,k0a1)
                        prcs(i+5).v  = rcs1.v 
                        a2.v         = pa(i+6).v
                        k0a2.v       = pk0a(i+6).v
                        rcs2         = rcs_f419_zmm16r4(a2,k0a2)
                        prcs(i+6).v  = rcs2.v
                        a3.v         = pa(i+7).v
                        k0a3.v       = pk0a(i+7).v
                        rcs3         = rcs_f419_zmm16r4(a3,k0a3)
                        prcs(i+7).v  = rcs3.v  
                        a0.v         = pa(i+8).v
                        k0a0.v       = pk0a(i+8).v
                        rcs0         = rcs_f419_zmm16r4(a0,k0a0)
                        prcs(i+8).v  = rcs0.v  
                        a1.v         = pa(i+9).v
                        k0a1.v       = pk0a(i+9).v
                        rcs1         = rcs_f419_zmm16r4(a1,k0a1)
                        prcs(i+9).v  = rcs1.v 
                        a2.v         = pa(i+10).v
                        k0a2.v       = pk0a(i+10).v
                        rcs2         = rcs_f419_zmm16r4(a2,k0a2)
                        prcs(i+10).v = rcs2.v
                        a3.v         = pa(i+11).v
                        k0a3.v       = pk0a(i+11).v
                        rcs3         = rcs_f419_zmm16r4(a3,k0a3)
                        prcs(i+11).v = rcs3.v  
                        a0.v         = pa(i+12).v
                        k0a0.v       = pk0a(i+12).v
                        rcs0         = rcs_f419_zmm16r4(a0,k0a0)
                        prcs(i+12).v = rcs0.v  
                        a1.v         = pa(i+13).v
                        k0a1.v       = pk0a(i+13).v
                        rcs1         = rcs_f419_zmm16r4(a1,k0a1)
                        prcs(i+13).v = rcs1.v 
                        a2.v         = pa(i+14).v
                        k0a2.v       = pk0a(i+14).v
                        rcs2         = rcs_f419_zmm16r4(a2,k0a2)
                        prcs(i+14).v = rcs2.v
                        a3.v         = pa(i+15).v
                        k0a3.v       = pk0a(i+15).v
                        rcs3         = rcs_f419_zmm16r4(a3,k0a3)
                        prcs(i+15).v = rcs3.v                      
                   end do
             end subroutine rcs_f419_zmm16r4_unroll16x


             subroutine rcs_f419_zmm16r4_unroll12x(pa,pk0a,prcs,n,PF_DIST)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f419_zmm16r4_unroll12x
                   !dir$ attributes forceinline :: rcs_f419_zmm16r4_unroll12x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f419_zmm16r4_unroll12x
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out)   :: prcs 
                   integer(kind=i4)               , intent(in)    :: n
                   integer(kind=i4)               , intent(in)    :: PF_DIST
                   ! Locals
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,12)
                   if(m/=0) then
                      do i=1,m
                         a0.v       = pa(i).v
                         k0a0.v     = pk0a(i).v
                         rcs0       = rcs_f419_zmm16r4(a0,k0a0)
                         prcs(i).v  = rcs0.v
                      end do
                      if(n<12) return
                   end if
                   m1=m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                   do i=m1,n,12
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif       
                        a0.v         = pa(i+0).v
                        k0a0.v       = pk0a(i+0).v
                        rcs0         = rcs_f419_zmm16r4(a0,k0a0)
                        prcs(i+0).v  = rcs0.v  
                        a1.v         = pa(i+1).v
                        k0a1.v       = pk0a(i+1).v
                        rcs1         = rcs_f419_zmm16r4(a1,k0a1)
                        prcs(i+1).v  = rcs1.v 
                        a2.v         = pa(i+2).v
                        k0a2.v       = pk0a(i+2).v
                        rcs2         = rcs_f419_zmm16r4(a2,k0a2)
                        prcs(i+2).v  = rcs2.v  
                        a3.v         = pa(i+3).v
                        k0a3.v       = pk0a(i+3).v
                        rcs3         = rcs_f419_zmm16r4(a3,k0a3)
                        prcs(i+3).v  = rcs3.v  
                        a0.v         = pa(i+4).v
                        k0a0.v       = pk0a(i+4).v
                        rcs0         = rcs_f419_zmm16r4(a0,k0a0)
                        prcs(i+4).v  = rcs0.v  
                        a1.v         = pa(i+5).v
                        k0a1.v       = pk0a(i+5).v
                        rcs1         = rcs_f419_zmm16r4(a1,k0a1)
                        prcs(i+5).v  = rcs1.v 
                        a2.v         = pa(i+6).v
                        k0a2.v       = pk0a(i+6).v
                        rcs2         = rcs_f419_zmm16r4(a2,k0a2)
                        prcs(i+6).v  = rcs2.v
                        a3.v         = pa(i+7).v
                        k0a3.v       = pk0a(i+7).v
                        rcs3         = rcs_f419_zmm16r4(a3,k0a3)
                        prcs(i+7).v  = rcs3.v  
                        a0.v         = pa(i+8).v
                        k0a0.v       = pk0a(i+8).v
                        rcs0         = rcs_f419_zmm16r4(a0,k0a0)
                        prcs(i+8).v  = rcs0.v  
                        a1.v         = pa(i+9).v
                        k0a1.v       = pk0a(i+9).v
                        rcs1         = rcs_f419_zmm16r4(a1,k0a1)
                        prcs(i+9).v  = rcs1.v 
                        a2.v         = pa(i+10).v
                        k0a2.v       = pk0a(i+10).v
                        rcs2         = rcs_f419_zmm16r4(a2,k0a2)
                        prcs(i+10).v = rcs2.v
                        a3.v         = pa(i+11).v
                        k0a3.v       = pk0a(i+11).v
                        rcs3         = rcs_f419_zmm16r4(a3,k0a3)
                        prcs(i+11).v = rcs3.v  
                    end do
             end subroutine rcs_f419_zmm16r4_unroll12x
             
             
             subroutine rcs_f419_zmm16r4_unroll8x(pa,pk0a,prcs,n,PF_DIST)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f419_zmm16r4_unroll8x
                   !dir$ attributes forceinline :: rcs_f419_zmm16r4_unroll8x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f419_zmm16r4_unroll8x
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out)   :: prcs 
                   integer(kind=i4)               , intent(in)    :: n
                   integer(kind=i4)               , intent(in)    :: PF_DIST
                   ! Locals
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,8)
                   if(m/=0) then
                      do i=1,m
                         a0.v       = pa(i).v
                         k0a0.v     = pk0a(i).v
                         rcs0       = rcs_f419_zmm16r4(a0,k0a0)
                         prcs(i).v  = rcs0.v
                      end do
                      if(n<8) return
                   end if
                   m1=m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                   do i=m1,n,8
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif       
                        a0.v         = pa(i+0).v
                        k0a0.v       = pk0a(i+0).v
                        rcs0         = rcs_f419_zmm16r4(a0,k0a0)
                        prcs(i+0).v  = rcs0.v  
                        a1.v         = pa(i+1).v
                        k0a1.v       = pk0a(i+1).v
                        rcs1         = rcs_f419_zmm16r4(a1,k0a1)
                        prcs(i+1).v  = rcs1.v 
                        a2.v         = pa(i+2).v
                        k0a2.v       = pk0a(i+2).v
                        rcs2         = rcs_f419_zmm16r4(a2,k0a2)
                        prcs(i+2).v  = rcs2.v  
                        a3.v         = pa(i+3).v
                        k0a3.v       = pk0a(i+3).v
                        rcs3         = rcs_f419_zmm16r4(a3,k0a3)
                        prcs(i+3).v  = rcs3.v  
                        a0.v         = pa(i+4).v
                        k0a0.v       = pk0a(i+4).v
                        rcs0         = rcs_f419_zmm16r4(a0,k0a0)
                        prcs(i+4).v  = rcs0.v  
                        a1.v         = pa(i+5).v
                        k0a1.v       = pk0a(i+5).v
                        rcs1         = rcs_f419_zmm16r4(a1,k0a1)
                        prcs(i+5).v  = rcs1.v 
                        a2.v         = pa(i+6).v
                        k0a2.v       = pk0a(i+6).v
                        rcs2         = rcs_f419_zmm16r4(a2,k0a2)
                        prcs(i+6).v  = rcs2.v
                        a3.v         = pa(i+7).v
                        k0a3.v       = pk0a(i+7).v
                        rcs3         = rcs_f419_zmm16r4(a3,k0a3)
                        prcs(i+7).v  = rcs3.v  
                    end do
             end subroutine rcs_f419_zmm16r4_unroll8x
             
             
             subroutine rcs_f419_zmm16r4_unroll4x(pa,pk0a,prcs,n,PF_DIST)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f419_zmm16r4_unroll4x
                   !dir$ attributes forceinline :: rcs_f419_zmm16r4_unroll4x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f419_zmm16r4_unroll4x
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out)   :: prcs 
                   integer(kind=i4)               , intent(in)    :: n
                   integer(kind=i4)               , intent(in)    :: PF_DIST
                   ! Locals
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,4)
                   if(m/=0) then
                      do i=1,m
                         a0.v       = pa(i).v
                         k0a0.v     = pk0a(i).v
                         rcs0       = rcs_f419_zmm16r4(a0,k0a0)
                         prcs(i).v  = rcs0.v
                      end do
                      if(n<4) return
                   end if
                   m1=m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                   do i=m1,n,4
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif       
                        a0.v         = pa(i+0).v
                        k0a0.v       = pk0a(i+0).v
                        rcs0         = rcs_f419_zmm16r4(a0,k0a0)
                        prcs(i+0).v  = rcs0.v  
                        a1.v         = pa(i+1).v
                        k0a1.v       = pk0a(i+1).v
                        rcs1         = rcs_f419_zmm16r4(a1,k0a1)
                        prcs(i+1).v  = rcs1.v 
                        a2.v         = pa(i+2).v
                        k0a2.v       = pk0a(i+2).v
                        rcs2         = rcs_f419_zmm16r4(a2,k0a2)
                        prcs(i+2).v  = rcs2.v  
                        a3.v         = pa(i+3).v
                        k0a3.v       = pk0a(i+3).v
                        rcs3         = rcs_f419_zmm16r4(a3,k0a3)
                        prcs(i+3).v  = rcs3.v  
                    end do
             end subroutine rcs_f419_zmm16r4_unroll4x
             
             
             subroutine rcs_f419_zmm16r4_rolled(pa,pk0a,prcs,n,PF_DIST)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f419_zmm16r4_unroll4x
                   !dir$ attributes forceinline :: rcs_f419_zmm16r4_unroll4x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f419_zmm16r4_unroll4x
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out)   :: prcs 
                   integer(kind=i4)               , intent(in)    :: n
                   integer(kind=i4)               , intent(in)    :: PF_DIST
                   ! Locals
                   type(ZMM16r4_t), automatic :: a0
                   type(ZMM16r4_t), automatic :: k0a0
                   type(ZMM16r4_t), automatic :: rcs0
                   integer(kind=i4) :: i
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                   do i=1,n
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif       
                        a0.v         = pa(i+0).v
                        k0a0.v       = pk0a(i+0).v
                        rcs0         = rcs_f419_zmm16r4(a0,k0a0)
                        prcs(i+0).v  = rcs0.v  
                     end do
             end subroutine rcs_f419_zmm16r4_rolled
             
             
             
             ! /* 
             !           Low frequency scattering widths (k0a << 1).
             !            Backscatter scattering width for H-field 
             !            cylinder-parallel,formula 4.1-20
             !       */
             
             
             pure function rcs_f4120_zmm16r4(a,k0a) result(rcs)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4120_zmm16r4
                   !dir$ attributes forceinline :: rcs_f4120_zmm16r4
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4120_zmm16r4
                   type(ZMM16r4_t),  intent(in) :: a
                   type(ZMM16r4_t),  intent(in) :: k0a
                   type(ZMM16r4_t) :: rcs
                   ! Locals
                   type(ZMM16r4_t), parameter :: C9869604401089358618834490999876 = &
                                             ZMM16r4_t(9.869604401089358618834490999876_sp)
                   type(ZMM16r4_t), parameter :: C225 = ZMM16r4_t(2.25_sp)
                   type(ZMM16r4_t), automatic :: pi2a,k0a3,t0
                   k0a3.v = k0a.v*k0a.v*k0a.v
                   t0.v   = C225.v*k0a3.v
                   pi2a.v = a.v*C9869604401089358618834490999876.v
                   rcs.v  = pi2a.v*t0.v
             end function rcs_f4120_zmm16r4
             
             
             subroutine rcs_f4120_zmm16r4_unroll16x(pa,pk0a,prcs,n,PF_DIST)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4120_zmm16r4_unroll16x
                   !dir$ attributes forceinline :: rcs_f4120_zmm16r4_unroll16x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4120_zmm16r4_unroll16x
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out):: prcs
                   integer(kind=i4),                intent(in) :: n
                   integer(kind=i4),                intent(in) :: PF_DIST
                   !Locals
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,16)
                   if(m/=0) then
                      do i=1,m
                         a0.v      = pa(i).v
                         k0a.v     = pk0a(i).v
                         rcs.v     = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i).v = rcs.v
                      end do
                      if(n<16) return
                   end if
                   m1 = m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                   do i=m1,n,16
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif                    
                         a0.v        = pa(i+0).v
                         k0a.v       = pk0a(i+0).v
                         rcs.v       = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i+0).v = rcs.v
                         a1.v        = pa(i+1).v
                         k0a1.v      = pk0a(i+1).v
                         rcs1.v      = rcs_f4120_zmm16r4(a1,k0a1)
                         prcs(i+1).v = rcs1.v
                         a2.v        = pa(i+2).v
                         k0a2.v      = pk0a(i+2).v
                         rcs2.v      = rcs_f4120_zmm16r4(a2,k0a2)
                         prcs(i+2).v = rcs2.v
                         a3.v        = pa(i+3).v
                         k0a3.v      = pk0a(i+3).v
                         rcs3.v      = rcs_f4120_zmm16r4(a03,k0a3)
                         prcs(i+3).v = rcs3.v
                         a0.v        = pa(i+4).v
                         k0a.v       = pk0a(i+4).v
                         rcs.v       = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i+4).v = rcs.v
                         a1.v        = pa(i+5).v
                         k0a1.v      = pk0a(i+5).v
                         rcs1.v      = rcs_f4120_zmm16r4(a1,k0a1)
                         prcs(i+5).v = rcs1.v
                         a2.v        = pa(i+6).v
                         k0a2.v      = pk0a(i+6).v
                         rcs2.v      = rcs_f4120_zmm16r4(a2,k0a2)
                         prcs(i+6).v = rcs2.v
                         a3.v        = pa(i+7).v
                         k0a3.v      = pk0a(i+7).v
                         rcs3.v      = rcs_f4120_zmm16r4(a03,k0a3)
                         prcs(i+7).v = rcs3.v
                         a0.v        = pa(i+8).v
                         k0a.v       = pk0a(i+8).v
                         rcs.v       = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i+8).v = rcs.v
                         a1.v        = pa(i+9).v
                         k0a1.v      = pk0a(i+9).v
                         rcs1.v      = rcs_f4120_zmm16r4(a1,k0a1)
                         prcs(i+9).v = rcs1.v
                         a2.v        = pa(i+10).v
                         k0a2.v      = pk0a(i+10).v
                         rcs2.v      = rcs_f4120_zmm16r4(a2,k0a2)
                         prcs(i+10).v = rcs2.v
                         a3.v        = pa(i+11).v
                         k0a3.v      = pk0a(i+11).v
                         rcs3.v      = rcs_f4120_zmm16r4(a03,k0a3)
                         prcs(i+11).v = rcs3.v
                         a0.v        = pa(i+12).v
                         k0a.v       = pk0a(i+12).v
                         rcs.v       = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i+12).v = rcs.v
                         a1.v        = pa(i+13).v
                         k0a1.v      = pk0a(i+13).v
                         rcs1.v      = rcs_f4120_zmm16r4(a1,k0a1)
                         prcs(i+13).v = rcs1.v
                         a2.v        = pa(i+14).v
                         k0a2.v      = pk0a(i+14).v
                         rcs2.v      = rcs_f4120_zmm16r4(a2,k0a2)
                         prcs(i+14).v = rcs2.v
                         a3.v        = pa(i+15).v
                         k0a3.v      = pk0a(i+15).v
                         rcs3.v      = rcs_f4120_zmm16r4(a03,k0a3)
                         prcs(i+15).v = rcs3.v
                   end do
             end subroutine rcs_f4120_zmm16r4_unroll16x
             
             
             subroutine rcs_f4120_zmm16r4_unroll12x(pa,pk0a,prcs,n,PF_DIST)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4120_zmm16r4_unroll12x
                   !dir$ attributes forceinline :: rcs_f4120_zmm16r4_unroll12x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4120_zmm16r4_unroll12x
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out):: prcs
                   integer(kind=i4),                intent(in) :: n
                   integer(kind=i4),                intent(in) :: PF_DIST
                   !Locals
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,12)
                   if(m/=0) then
                      do i=1,m
                         a0.v      = pa(i).v
                         k0a.v     = pk0a(i).v
                         rcs.v     = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i).v = rcs.v
                      end do
                      if(n<12) return
                   end if
                   m1 = m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                   do i=m1,n,12
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif                    
                         a0.v        = pa(i+0).v
                         k0a.v       = pk0a(i+0).v
                         rcs.v       = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i+0).v = rcs.v
                         a1.v        = pa(i+1).v
                         k0a1.v      = pk0a(i+1).v
                         rcs1.v      = rcs_f4120_zmm16r4(a1,k0a1)
                         prcs(i+1).v = rcs1.v
                         a2.v        = pa(i+2).v
                         k0a2.v      = pk0a(i+2).v
                         rcs2.v      = rcs_f4120_zmm16r4(a2,k0a2)
                         prcs(i+2).v = rcs2.v
                         a3.v        = pa(i+3).v
                         k0a3.v      = pk0a(i+3).v
                         rcs3.v      = rcs_f4120_zmm16r4(a03,k0a3)
                         prcs(i+3).v = rcs3.v
                         a0.v        = pa(i+4).v
                         k0a.v       = pk0a(i+4).v
                         rcs.v       = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i+4).v = rcs.v
                         a1.v        = pa(i+5).v
                         k0a1.v      = pk0a(i+5).v
                         rcs1.v      = rcs_f4120_zmm16r4(a1,k0a1)
                         prcs(i+5).v = rcs1.v
                         a2.v        = pa(i+6).v
                         k0a2.v      = pk0a(i+6).v
                         rcs2.v      = rcs_f4120_zmm16r4(a2,k0a2)
                         prcs(i+6).v = rcs2.v
                         a3.v        = pa(i+7).v
                         k0a3.v      = pk0a(i+7).v
                         rcs3.v      = rcs_f4120_zmm16r4(a03,k0a3)
                         prcs(i+7).v = rcs3.v
                         a0.v        = pa(i+8).v
                         k0a.v       = pk0a(i+8).v
                         rcs.v       = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i+8).v = rcs.v
                         a1.v        = pa(i+9).v
                         k0a1.v      = pk0a(i+9).v
                         rcs1.v      = rcs_f4120_zmm16r4(a1,k0a1)
                         prcs(i+9).v = rcs1.v
                         a2.v        = pa(i+10).v
                         k0a2.v      = pk0a(i+10).v
                         rcs2.v      = rcs_f4120_zmm16r4(a2,k0a2)
                         prcs(i+10).v = rcs2.v
                         a3.v        = pa(i+11).v
                         k0a3.v      = pk0a(i+11).v
                         rcs3.v      = rcs_f4120_zmm16r4(a03,k0a3)
                         prcs(i+11).v = rcs3.v
                     end do
             end subroutine rcs_f4120_zmm16r4_unroll12x
             

             subroutine rcs_f4120_zmm16r4_unroll8x(pa,pk0a,prcs,n,PF_DIST)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4120_zmm16r4_unroll8x
                   !dir$ attributes forceinline :: rcs_f4120_zmm16r4_unroll8x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4120_zmm16r4_unroll8x
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out):: prcs
                   integer(kind=i4),                intent(in) :: n
                   integer(kind=i4),                intent(in) :: PF_DIST
                   !Locals
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,8)
                   if(m/=0) then
                      do i=1,m
                         a0.v      = pa(i).v
                         k0a.v     = pk0a(i).v
                         rcs.v     = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i).v = rcs.v
                      end do
                      if(n<8) return
                   end if
                   m1 = m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                   do i=m1,n,8
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif                    
                         a0.v        = pa(i+0).v
                         k0a.v       = pk0a(i+0).v
                         rcs.v       = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i+0).v = rcs.v
                         a1.v        = pa(i+1).v
                         k0a1.v      = pk0a(i+1).v
                         rcs1.v      = rcs_f4120_zmm16r4(a1,k0a1)
                         prcs(i+1).v = rcs1.v
                         a2.v        = pa(i+2).v
                         k0a2.v      = pk0a(i+2).v
                         rcs2.v      = rcs_f4120_zmm16r4(a2,k0a2)
                         prcs(i+2).v = rcs2.v
                         a3.v        = pa(i+3).v
                         k0a3.v      = pk0a(i+3).v
                         rcs3.v      = rcs_f4120_zmm16r4(a03,k0a3)
                         prcs(i+3).v = rcs3.v
                         a0.v        = pa(i+4).v
                         k0a.v       = pk0a(i+4).v
                         rcs.v       = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i+4).v = rcs.v
                         a1.v        = pa(i+5).v
                         k0a1.v      = pk0a(i+5).v
                         rcs1.v      = rcs_f4120_zmm16r4(a1,k0a1)
                         prcs(i+5).v = rcs1.v
                         a2.v        = pa(i+6).v
                         k0a2.v      = pk0a(i+6).v
                         rcs2.v      = rcs_f4120_zmm16r4(a2,k0a2)
                         prcs(i+6).v = rcs2.v
                         a3.v        = pa(i+7).v
                         k0a3.v      = pk0a(i+7).v
                         rcs3.v      = rcs_f4120_zmm16r4(a03,k0a3)
                         prcs(i+7).v = rcs3.v
                     end do
             end subroutine rcs_f4120_zmm16r4_unroll8x
             

             subroutine rcs_f4120_zmm16r4_unroll4x(pa,pk0a,prcs,n,PF_DIST)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4120_zmm16r4_unroll4x
                   !dir$ attributes forceinline :: rcs_f4120_zmm16r4_unroll4x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4120_zmm16r4_unroll4x
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out):: prcs
                   integer(kind=i4),                intent(in) :: n
                   integer(kind=i4),                intent(in) :: PF_DIST
                   !Locals
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,4)
                   if(m/=0) then
                      do i=1,m
                         a0.v      = pa(i).v
                         k0a.v     = pk0a(i).v
                         rcs.v     = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i).v = rcs.v
                      end do
                      if(n<4) return
                   end if
                   m1 = m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                   do i=m1,n,4
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif                    
                         a0.v        = pa(i+0).v
                         k0a.v       = pk0a(i+0).v
                         rcs.v       = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i+0).v = rcs.v
                         a1.v        = pa(i+1).v
                         k0a1.v      = pk0a(i+1).v
                         rcs1.v      = rcs_f4120_zmm16r4(a1,k0a1)
                         prcs(i+1).v = rcs1.v
                         a2.v        = pa(i+2).v
                         k0a2.v      = pk0a(i+2).v
                         rcs2.v      = rcs_f4120_zmm16r4(a2,k0a2)
                         prcs(i+2).v = rcs2.v
                         a3.v        = pa(i+3).v
                         k0a3.v      = pk0a(i+3).v
                         rcs3.v      = rcs_f4120_zmm16r4(a03,k0a3)
                         prcs(i+3).v = rcs3.v
                      end do
             end subroutine rcs_f4120_zmm16r4_unroll4x
             
         
             subroutine rcs_f4120_zmm16r4_rolled(pa,pk0a,prcs,n,PF_DIST)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4120_zmm16r4_rolled
                   !dir$ attributes forceinline :: rcs_f4120_zmm16r4_rolled
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4120_zmm16r4_rolled
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out):: prcs
                   integer(kind=i4),                intent(in) :: n
                   integer(kind=i4),                intent(in) :: PF_DIST
                   !Locals
                   type(ZMM16r4_t), automatic :: a0
                   type(ZMM16r4_t), automatic :: k0a0
                   type(ZMM16r4_t), automatic :: rcs0
                   integer(kind=i4) :: i
                  
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                   do i=1,n
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif                    
                         a0.v        = pa(i+0).v
                         k0a.v       = pk0a(i+0).v
                         rcs.v       = rcs_f4120_zmm16r4(a0,k0a)
                         prcs(i+0).v = rcs.v
                     end do
             end subroutine rcs_f4120_zmm16r4_rolled
             
             
             ! /*
             !           Bistatic scattering widths, E-field cylinder axis-parallel
             !           Formula 4.1-21
             !      */
             
             pure function rcs_f4121_zmm16r4(a,k0a) result(rcs)
                  
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4121_zmm16r4
                   !dir$ attributes forceinline :: rcs_f4121_zmm16r4
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4121_zmm16r4
                   type(ZMM16r4_t), intent(in) :: a
                   type(ZMM16r4_t), intent(in) :: k0a
                   type(ZMM16r4_t) :: rcs
                   rcs = rcs_f4120_zmm16r4(a,k0a)
             end function rcs_f4121_zmm16r4
             
             
              !/*
              !          Bistatic scattering widths, H-field cylinder axis-parallel
              !          Formula 4.1-22
              !     */ 
              
              
              pure function rcs_f4122_zmm16r4(phi,a,k0a) result(rcs)
              
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4122_zmm16r4
                   !dir$ attributes forceinline :: rcs_f4122_zmm16r4
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4122_zmm16r4 
                   type(ZMM16r4_t),   intent(in) :: phi
                   type(ZMM16r4_t),   intent(in) :: a
                   type(ZMM16r4_t),   intent(in) :: k0a
                   type(ZMM16r4_t) :: rcs
                   ! Locals
                   type(ZMM16r4_t), parameter :: C9869604401089358618834490999876 = &
                                              ZMM16r4_t(9.869604401089358618834490999876_sp)
                   type(ZMM16r4_t), parameter :: C05 = ZMM16r4_t(0.5_sp)
                   type(ZMM16r4_t), automatic :: pi2a,k0a3,cosp,frac,sqr
                   pi2a.v = a.v*C9869604401089358618834490999876.v
                   k0a3.v = k0a.v*k0a.v*k0a.v
                   cosp.v = cos(phi.v)
                   frac.v = C05.v*cosp.v
                   sqr.v  = frac.v*frac.v
                   rcs.v  = pi2a.v*k0a3.v*sqr.v
              end function rcs_f4122_zmm16r4
              
              
              subroutine rcs_f4122_zmm16r4_unroll16x(pphi,pa,pk0a,prcs,n,PF_DIST)
              
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4122_zmm16r4_unroll16x
                   !dir$ attributes forceinline :: rcs_f4122_zmm16r4_unroll16x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4122_zmm16r4_unroll16x
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pphi
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out):: prcs
                   integer(kind=i4),                intent(in) :: n
                   integer(kind=i4),                intent(in) :: PF_DIST
                   ! Locals
                   type(ZMM16r4_t), automatic :: phi0,phi1,phi2,phi3
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,16)
                   if(m/=0) then
                      do i=1,m
                         phi0.v    = pphi(i).v
                         a0.v      = pa(i).v
                         k0a0.v    = pk0a(i).v
                         rcs0.v    = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i).v = rcs0.v
                      end do
                      if(n<16) return
                   end if
                   m1 = m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ assume_aligned pphi:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                  do i=m1,n,16
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif                  
                         phi0.v      = pphi(i+0).v
                         a0.v        = pa(i+0).v
                         k0a0.v      = pk0a(i+0).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+0).v = rcs0.v
                         phi1.v      = pphi(i+1).v
                         a1.v        = pa(i+1).v
                         k0a1.v      = pk0a(i+1).v
                         rcs1.v      = rcs_f4122_zmm16r4(phi1,a1,k0a1)
                         prcs(i+1).v = rcs1.v
                         phi2.v      = pphi(i+2).v
                         a2.v        = pa(i+2).v
                         k0a2.v      = pk0a(i+2).v
                         rcs2.v      = rcs_f4122_zmm16r4(phi2,a2,k0a2)
                         prcs(i+2).v = rcs2.v   
                         phi3.v      = pphi(i+3).v
                         a3.v        = pa(i+3).v
                         k0a3.v      = pk0a(i+3).v
                         rcs3.v      = rcs_f4122_zmm16r4(phi3,a3,k0a3)
                         prcs(i+3).v = rcs3.v
                         phi0.v      = pphi(i+4).v
                         a0.v        = pa(i+4).v
                         k0a0.v      = pk0a(i+4).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+4).v = rcs0.v
                         phi1.v      = pphi(i+5).v
                         a1.v        = pa(i+5).v
                         k0a1.v      = pk0a(i+5).v
                         rcs1.v      = rcs_f4122_zmm16r4(phi1,a1,k0a1)
                         prcs(i+5).v = rcs1.v
                         phi2.v      = pphi(i+6).v
                         a2.v        = pa(i+6).v
                         k0a2.v      = pk0a(i+6).v
                         rcs2.v      = rcs_f4122_zmm16r4(phi2,a2,k0a2)
                         prcs(i+6).v = rcs2.v   
                         phi3.v      = pphi(i+7).v
                         a3.v        = pa(i+7).v
                         k0a3.v      = pk0a(i+7).v
                         rcs3.v      = rcs_f4122_zmm16r4(phi3,a3,k0a3)
                         prcs(i+7).v = rcs3.v
                         phi0.v      = pphi(i+8).v
                         a0.v        = pa(i+8).v
                         k0a0.v      = pk0a(i+8).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+8).v = rcs0.v
                         phi1.v      = pphi(i+9).v
                         a1.v        = pa(i+9).v
                         k0a1.v      = pk0a(i+9).v
                         rcs1.v      = rcs_f4122_zmm16r4(phi1,a1,k0a1)
                         prcs(i+9).v = rcs1.v
                         phi2.v      = pphi(i+10).v
                         a2.v        = pa(i+10).v
                         k0a2.v      = pk0a(i+10).v
                         rcs2.v      = rcs_f4122_zmm16r4(phi2,a2,k0a2)
                         prcs(i+10).v= rcs2.v   
                         phi3.v      = pphi(i+11).v
                         a3.v        = pa(i+11).v
                         k0a3.v      = pk0a(i+11).v
                         rcs3.v      = rcs_f4122_zmm16r4(phi3,a3,k0a3)
                         prcs(i+11).v = rcs3.v
                         phi0.v      = pphi(i+12).v
                         a0.v        = pa(i+12).v
                         k0a0.v      = pk0a(i+12).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+12).v = rcs0.v
                         phi1.v      = pphi(i+13).v
                         a1.v        = pa(i+13).v
                         k0a1.v      = pk0a(i+13).v
                         rcs1.v      = rcs_f4122_zmm16r4(phi1,a1,k0a1)
                         prcs(i+13).v = rcs1.v
                         phi2.v      = pphi(i+14).v
                         a2.v        = pa(i+14).v
                         k0a2.v      = pk0a(i+14).v
                         rcs2.v      = rcs_f4122_zmm16r4(phi2,a2,k0a2)
                         prcs(i+14).v = rcs2.v   
                         phi3.v      = pphi(i+15).v
                         a3.v        = pa(i+15).v
                         k0a3.v      = pk0a(i+15).v
                         rcs3.v      = rcs_f4122_zmm16r4(phi3,a3,k0a3)
                         prcs(i+15).v = rcs3.v
                  end do
              end subroutine rcs_f4122_zmm16r4_unroll16x
              
             
              subroutine rcs_f4122_zmm16r4_unroll12x(pphi,pa,pk0a,prcs,n,PF_DIST)
              
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4122_zmm16r4_unroll12x
                   !dir$ attributes forceinline :: rcs_f4122_zmm16r4_unroll12x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4122_zmm16r4_unroll12x
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pphi
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out):: prcs
                   integer(kind=i4),                intent(in) :: n
                   integer(kind=i4),                intent(in) :: PF_DIST
                   ! Locals
                   type(ZMM16r4_t), automatic :: phi0,phi1,phi2,phi3
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,12)
                   if(m/=0) then
                      do i=1,m
                         phi0.v    = pphi(i).v
                         a0.v      = pa(i).v
                         k0a0.v    = pk0a(i).v
                         rcs0.v    = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i).v = rcs0.v
                      end do
                      if(n<12) return
                   end if
                   m1 = m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ assume_aligned pphi:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                  do i=m1,n,12
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif                  
                         phi0.v      = pphi(i+0).v
                         a0.v        = pa(i+0).v
                         k0a0.v      = pk0a(i+0).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+0).v = rcs0.v
                         phi1.v      = pphi(i+1).v
                         a1.v        = pa(i+1).v
                         k0a1.v      = pk0a(i+1).v
                         rcs1.v      = rcs_f4122_zmm16r4(phi1,a1,k0a1)
                         prcs(i+1).v = rcs1.v
                         phi2.v      = pphi(i+2).v
                         a2.v        = pa(i+2).v
                         k0a2.v      = pk0a(i+2).v
                         rcs2.v      = rcs_f4122_zmm16r4(phi2,a2,k0a2)
                         prcs(i+2).v = rcs2.v   
                         phi3.v      = pphi(i+3).v
                         a3.v        = pa(i+3).v
                         k0a3.v      = pk0a(i+3).v
                         rcs3.v      = rcs_f4122_zmm16r4(phi3,a3,k0a3)
                         prcs(i+3).v = rcs3.v
                         phi0.v      = pphi(i+4).v
                         a0.v        = pa(i+4).v
                         k0a0.v      = pk0a(i+4).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+4).v = rcs0.v
                         phi1.v      = pphi(i+5).v
                         a1.v        = pa(i+5).v
                         k0a1.v      = pk0a(i+5).v
                         rcs1.v      = rcs_f4122_zmm16r4(phi1,a1,k0a1)
                         prcs(i+5).v = rcs1.v
                         phi2.v      = pphi(i+6).v
                         a2.v        = pa(i+6).v
                         k0a2.v      = pk0a(i+6).v
                         rcs2.v      = rcs_f4122_zmm16r4(phi2,a2,k0a2)
                         prcs(i+6).v = rcs2.v   
                         phi3.v      = pphi(i+7).v
                         a3.v        = pa(i+7).v
                         k0a3.v      = pk0a(i+7).v
                         rcs3.v      = rcs_f4122_zmm16r4(phi3,a3,k0a3)
                         prcs(i+7).v = rcs3.v
                         phi0.v      = pphi(i+8).v
                         a0.v        = pa(i+8).v
                         k0a0.v      = pk0a(i+8).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+8).v = rcs0.v
                         phi1.v      = pphi(i+9).v
                         a1.v        = pa(i+9).v
                         k0a1.v      = pk0a(i+9).v
                         rcs1.v      = rcs_f4122_zmm16r4(phi1,a1,k0a1)
                         prcs(i+9).v = rcs1.v
                         phi2.v      = pphi(i+10).v
                         a2.v        = pa(i+10).v
                         k0a2.v      = pk0a(i+10).v
                         rcs2.v      = rcs_f4122_zmm16r4(phi2,a2,k0a2)
                         prcs(i+10).v= rcs2.v   
                         phi3.v      = pphi(i+11).v
                         a3.v        = pa(i+11).v
                         k0a3.v      = pk0a(i+11).v
                         rcs3.v      = rcs_f4122_zmm16r4(phi3,a3,k0a3)
                         prcs(i+11).v = rcs3.v
                      end do
              end subroutine rcs_f4122_zmm16r4_unroll12x
              

              subroutine rcs_f4122_zmm16r4_unroll8x(pphi,pa,pk0a,prcs,n,PF_DIST)
              
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4122_zmm16r4_unroll8x
                   !dir$ attributes forceinline :: rcs_f4122_zmm16r4_unroll8x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4122_zmm16r4_unroll8x
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pphi
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out):: prcs
                   integer(kind=i4),                intent(in) :: n
                   integer(kind=i4),                intent(in) :: PF_DIST
                   ! Locals
                   type(ZMM16r4_t), automatic :: phi0,phi1,phi2,phi3
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,8)
                   if(m/=0) then
                      do i=1,m
                         phi0.v    = pphi(i).v
                         a0.v      = pa(i).v
                         k0a0.v    = pk0a(i).v
                         rcs0.v    = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i).v = rcs0.v
                      end do
                      if(n<8) return
                   end if
                   m1 = m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ assume_aligned pphi:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                  do i=m1,n,8
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif                  
                         phi0.v      = pphi(i+0).v
                         a0.v        = pa(i+0).v
                         k0a0.v      = pk0a(i+0).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+0).v = rcs0.v
                         phi1.v      = pphi(i+1).v
                         a1.v        = pa(i+1).v
                         k0a1.v      = pk0a(i+1).v
                         rcs1.v      = rcs_f4122_zmm16r4(phi1,a1,k0a1)
                         prcs(i+1).v = rcs1.v
                         phi2.v      = pphi(i+2).v
                         a2.v        = pa(i+2).v
                         k0a2.v      = pk0a(i+2).v
                         rcs2.v      = rcs_f4122_zmm16r4(phi2,a2,k0a2)
                         prcs(i+2).v = rcs2.v   
                         phi3.v      = pphi(i+3).v
                         a3.v        = pa(i+3).v
                         k0a3.v      = pk0a(i+3).v
                         rcs3.v      = rcs_f4122_zmm16r4(phi3,a3,k0a3)
                         prcs(i+3).v = rcs3.v
                         phi0.v      = pphi(i+4).v
                         a0.v        = pa(i+4).v
                         k0a0.v      = pk0a(i+4).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+4).v = rcs0.v
                         phi1.v      = pphi(i+5).v
                         a1.v        = pa(i+5).v
                         k0a1.v      = pk0a(i+5).v
                         rcs1.v      = rcs_f4122_zmm16r4(phi1,a1,k0a1)
                         prcs(i+5).v = rcs1.v
                         phi2.v      = pphi(i+6).v
                         a2.v        = pa(i+6).v
                         k0a2.v      = pk0a(i+6).v
                         rcs2.v      = rcs_f4122_zmm16r4(phi2,a2,k0a2)
                         prcs(i+6).v = rcs2.v   
                         phi3.v      = pphi(i+7).v
                         a3.v        = pa(i+7).v
                         k0a3.v      = pk0a(i+7).v
                         rcs3.v      = rcs_f4122_zmm16r4(phi3,a3,k0a3)
                         prcs(i+7).v = rcs3.v
                     end do
              end subroutine rcs_f4122_zmm16r4_unroll8x
              

              subroutine rcs_f4122_zmm16r4_unroll4x(pphi,pa,pk0a,prcs,n,PF_DIST)
              
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4122_zmm16r4_unroll4x
                   !dir$ attributes forceinline :: rcs_f4122_zmm16r4_unroll4x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4122_zmm16r4_unroll4x
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pphi
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out):: prcs
                   integer(kind=i4),                intent(in) :: n
                   integer(kind=i4),                intent(in) :: PF_DIST
                   ! Locals
                   type(ZMM16r4_t), automatic :: phi0,phi1,phi2,phi3
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,4)
                   if(m/=0) then
                      do i=1,m
                         phi0.v    = pphi(i).v
                         a0.v      = pa(i).v
                         k0a0.v    = pk0a(i).v
                         rcs0.v    = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i).v = rcs0.v
                      end do
                      if(n<4) return
                   end if
                   m1 = m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ assume_aligned pphi:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                  do i=m1,n,4
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif                  
                         phi0.v      = pphi(i+0).v
                         a0.v        = pa(i+0).v
                         k0a0.v      = pk0a(i+0).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+0).v = rcs0.v
                         phi1.v      = pphi(i+1).v
                         a1.v        = pa(i+1).v
                         k0a1.v      = pk0a(i+1).v
                         rcs1.v      = rcs_f4122_zmm16r4(phi1,a1,k0a1)
                         prcs(i+1).v = rcs1.v
                         phi2.v      = pphi(i+2).v
                         a2.v        = pa(i+2).v
                         k0a2.v      = pk0a(i+2).v
                         rcs2.v      = rcs_f4122_zmm16r4(phi2,a2,k0a2)
                         prcs(i+2).v = rcs2.v   
                         phi3.v      = pphi(i+3).v
                         a3.v        = pa(i+3).v
                         k0a3.v      = pk0a(i+3).v
                         rcs3.v      = rcs_f4122_zmm16r4(phi3,a3,k0a3)
                         prcs(i+3).v = rcs3.v
                         phi0.v      = pphi(i+4).v
                         a0.v        = pa(i+4).v
                         k0a0.v      = pk0a(i+4).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+4).v = rcs0.v
                     end do
              end subroutine rcs_f4122_zmm16r4_unroll4x
              

              subroutine rcs_f4122_zmm16r4_rolled(pphi,pa,pk0a,prcs,n,PF_DIST)
              
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4122_zmm16r4_rolled
                   !dir$ attributes forceinline :: rcs_f4122_zmm16r4_rolled
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4122_zmm16r4_rolled
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pphi
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in) :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out):: prcs
                   integer(kind=i4),                intent(in) :: n
                   integer(kind=i4),                intent(in) :: PF_DIST
                   ! Locals
                   type(ZMM16r4_t), automatic :: phi0
                   type(ZMM16r4_t), automatic :: a0,
                   type(ZMM16r4_t), automatic :: k0a0
                   type(ZMM16r4_t), automatic :: rcs0
                   integer(kind=i4) :: i
                  
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ assume_aligned pphi:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                  do i=1,n
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pphi(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif                  
                         phi0.v      = pphi(i+0).v
                         a0.v        = pa(i+0).v
                         k0a0.v      = pk0a(i+0).v
                         rcs0.v      = rcs_f4122_zmm16r4(phi0,a0,k0a0)
                         prcs(i+0).v = rcs0.v
                     end do
              end subroutine rcs_f4122_zmm16r4_rolled
              
              
               ! /*
               !        Forward scattering widths, E-field.
               !        Formula 4.1-23
               !    */
               
               pure function rcs_f4123_zmm16r4(a,k0a) result(rcs)
               
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4123_zmm16r4
                   !dir$ attributes forceinline :: rcs_f4123_zmm16r4
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4123_zmm16r4
                   type(ZMM16r4_t), intent(in) :: a
                   type(ZMM16r4_t), intent(in) :: k0a
                   type(ZMM16r4_t) :: rcs
                   rcs = rcs_f4120_zmm16r4(a,k0a)
               end function rcs_f4123_zmm16r4
               
               
               !/*
               !        Forward scattering widths, H-field.
               !        Formula 4.1-24
               !    */
               
               pure function rcs_f4124_zmm16r4(a,k0a) result(rcs)
               
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f4124_zmm16r4
                   !dir$ attributes forceinline :: rcs_f4124_zmm16r4
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f4124_zmm16r4
                   type(ZMM16r4_t), intent(in) :: a
                   type(ZMM16r4_t), intent(in) :: k0a
                   type(ZMM16r4_t) :: rcs
                   ! Locals
                   type(ZMM16r4_t), parameter :: C9869604401089358618834490999876 = &
                                                   ZMM16r4_t(9.869604401089358618834490999876_sp)
                   type(ZMM16r4_t), parameter :: C025 = ZMM16r4_t(0.25_sp)
                   type(ZMM16r4_t), automatic :: pi2a,k0a3
                   k0a3.v = k0a.v*k0a.v*k0a.v
                   pi2a.v = C9869604401089358618834490999876.v*a.v
                   rcs.v  = pi2a.v*k0a3.v*C025.v
               end function rcs_f4124_zmm16r4
               
               
               subroutine rcs_f424_zmm16r4_unroll16x(pa,pk0a,prcs,n,PF_DIST)
             
                   !dir$ optimize:3
                   !dir$ attributes code_align : 32 :: rcs_f424_zmm16r4_unroll16x
                   !dir$ attributes forceinline :: rcs_f424_zmm16r4_unroll16x
                   !dir$ attributes optimization_parameter:"target_arch=skylake-avx512" :: rcs_f424_zmm16r4_unroll16x
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pa
                   type(ZMM16r4_t), dimension(1:n), intent(in)    :: pk0a
                   type(ZMM16r4_t), dimension(1:n), intent(out)   :: prcs 
                   integer(kind=i4)               , intent(in)    :: n
                   integer(kind=i4)               , intent(in)    :: PF_DIST
                   ! Locals
                   type(ZMM16r4_t), automatic :: a0,a1,a2,a3
                   type(ZMM16r4_t), automatic :: k0a0,k0a1,k0a2,k0a3
                   type(ZMM16r4_t), automatic :: rcs0,rcs1,rcs2,rcs3
                   integer(kind=i4) :: i,m,m1
                   m = mod(n,16)
                   if(m/=0) then
                      do i=1,m
                         a0.v       = pa(i).v
                         k0a0.v     = pk0a(i).v
                         rcs0       = rcs_f424_zmm16r4(a0,k0a0)
                         prcs(i).v  = rcs0.v
                      end do
                      if(n<16) return
                   end if
                   m1=m+1
                    !dir$ assume_aligned pa:64
                    !dir$ assume_aligned pk0a:64
                    !dir$ assume_aligned prcs:64
                    !dir$ vector aligned
                    !dir$ ivdep
                    !dir$ vector vectorlength(4)
                    !dir$ vector multiple_gather_scatter_by_shuffles 
                    !dir$ vector always
                   do i=m1,n,16
#if (__RCS_CYLINDER_PF_CACHE_HINT__) == 1
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T0,.true.,.false.)
#elif (__RCS_CYLINDER_PF_CACHE_HINT__) == 2   
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T1,.true.,.false.)   
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 3
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_T2,.true.,.false.) 
#elif (_RCS_CYLINDER_PF_CACHE_HINT__) == 4
                       call mm_prefetch(pa(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.)
                       call mm_prefetch(pk0a(i+PF_DIST),FOR_K_PREFETCH_NTA,.true.,.false.) 
#endif       
                        a0.v         = pa(i+0).v
                        k0a0.v       = pk0a(i+0).v
                        rcs0         = rcs_f424_zmm16r4(a0,k0a0)
                        prcs(i+0).v  = rcs0.v  
                        a1.v         = pa(i+1).v
                        k0a1.v       = pk0a(i+1).v
                        rcs1         = rcs_f424_zmm16r4(a1,k0a1)
                        prcs(i+1).v  = rcs1.v 
                        a2.v         = pa(i+2).v
                        k0a2.v       = pk0a(i+2).v
                        rcs2         = rcs_f424_zmm16r4(a2,k0a2)
                        prcs(i+2).v  = rcs2.v  
                        a3.v         = pa(i+3).v
                        k0a3.v       = pk0a(i+3).v
                        rcs3         = rcs_f424_zmm16r4(a3,k0a3)
                        prcs(i+3).v  = rcs3.v  
                        a0.v         = pa(i+4).v
                        k0a0.v       = pk0a(i+4).v
                        rcs0         = rcs_f424_zmm16r4(a0,k0a0)
                        prcs(i+4).v  = rcs0.v  
                        a1.v         = pa(i+5).v
                        k0a1.v       = pk0a(i+5).v
                        rcs1         = rcs_f424_zmm16r4(a1,k0a1)
                        prcs(i+5).v  = rcs1.v 
                        a2.v         = pa(i+6).v
                        k0a2.v       = pk0a(i+6).v
                        rcs2         = rcs_f424_zmm16r4(a2,k0a2)
                        prcs(i+6).v  = rcs2.v
                        a3.v         = pa(i+7).v
                        k0a3.v       = pk0a(i+7).v
                        rcs3         = rcs_f424_zmm16r4(a3,k0a3)
                        prcs(i+7).v  = rcs3.v  
                        a0.v         = pa(i+8).v
                        k0a0.v       = pk0a(i+8).v
                        rcs0         = rcs_f424_zmm16r4(a0,k0a0)
                        prcs(i+8).v  = rcs0.v  
                        a1.v         = pa(i+9).v
                        k0a1.v       = pk0a(i+9).v
                        rcs1         = rcs_f424_zmm16r4(a1,k0a1)
                        prcs(i+9).v  = rcs1.v 
                        a2.v         = pa(i+10).v
                        k0a2.v       = pk0a(i+10).v
                        rcs2         = rcs_f424_zmm16r4(a2,k0a2)
                        prcs(i+10).v = rcs2.v
                        a3.v         = pa(i+11).v
                        k0a3.v       = pk0a(i+11).v
                        rcs3         = rcs_f424_zmm16r4(a3,k0a3)
                        prcs(i+11).v = rcs3.v  
                        a0.v         = pa(i+12).v
                        k0a0.v       = pk0a(i+12).v
                        rcs0         = rcs_f424_zmm16r4(a0,k0a0)
                        prcs(i+12).v = rcs0.v  
                        a1.v         = pa(i+13).v
                        k0a1.v       = pk0a(i+13).v
                        rcs1         = rcs_f424_zmm16r4(a1,k0a1)
                        prcs(i+13).v = rcs1.v 
                        a2.v         = pa(i+14).v
                        k0a2.v       = pk0a(i+14).v
                        rcs2         = rcs_f424_zmm16r4(a2,k0a2)
                        prcs(i+14).v = rcs2.v
                        a3.v         = pa(i+15).v
                        k0a3.v       = pk0a(i+15).v
                        rcs3         = rcs_f424_zmm16r4(a3,k0a3)
                        prcs(i+15).v = rcs3.v                      
                   end do
             end subroutine rcs_f424_zmm16r4_unroll16x

               
               
               
              










end module rcs_cylinder_zmm16r4
