!***************************************************************************************************
!>
!  A suite of Fortran routines for interpolating values between
!  one-, two-, three-, and four-dimensional arrays defined on uniform or nonuniform
!  orthogonal grids. This operation is commonly referred to as "regridding." Linear
!  or cubic interpolation can be selected independently in each dimension.
!  Extrapolation is not allowed. The subroutines in REGRIDPACK cannot be used to
!  transfer values on nonorthogonal (randomly scattered) data grids.
!
!### History
!  * John C. Adams (NCAR 1997) : original REGRIDPACK
!  * Jacob Williams, Oct 2019, modernized and refactored
!  * Bernard Gingold, Jun 2022, minor adaptations and optimizations
!                               to better suite my simulation projects
!                               Added OpenMP parallelization and 
!                               optimization ifort directives.

    module regridpack

    
    use mod_kinds, only : i4, dp
    implicit none

    private

   ! interface regrid

        !low level routines:
   !    module procedure :: rgrd1, rgrd2, rgrd3 !, rgrd4
   !    module procedure :: rgrd1u, rgrd2u, rgrd3u! , rgrd4u

   !   module procedure :: rgrd1_wrapper, rgrd2_wrapper, rgrd3_wrapper !, rgrd4_wrapper

   !end interface
   ! public :: regrid
    public :: rgrd1_wrapper, rgrd2_wrapper, rgrd3_wrapper 
    public :: rgrd1, rgrd2, rgrd3
    public :: rgrd1u, rgrd2u, rgrd3u

#ifndef (AUTO_VECTORIZE)
#define AUTO_VECTORIZE 0
#endif

! Forced Compiler auto-vectorization is possible, although
! highly counterproductive, which is caused by the gather
! indexing operations.
!=========== An example ====================!

#if 0

    subroutine lint1(nx,p,mx,q,ix,dx)
       !dir$ attributes forceinline :: lint1
       !dir$ attributes code_align : 32 :: lint1
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: lint1
    implicit none

    integer(kind=4) :: mx,ix(mx),nx,ii,i
    real(kind=8) :: p(nx),q(mx),dx(mx)
    real(kind=8), automatic :: pi,pi1
    !dir$ assume_aligned q:64
    !dir$ assume_aligned p:64
    !dir$ assume_aligned ix:64
    !dir$ assume_aligned dx:64
    !dir$ code_align(32)
   !dir$ ivdep
   !dir$ vector aligned
   !dir$ vector always
    !dir$ unroll(4)
    !!dir$ novector
    do ii=1,mx
        i = ix(ii)
        pi = p(i)
        pi1 =p(i+1)
        !q(ii) = p(i)+dx(ii)*(p(i+1)-p(i))
        q(ii) = pi+dx(ii)*(pi1-pi)
    end do

    end subroutine lint1  

    !Corresponding dissassembly i.e. loop only.

 vmovq  xmm4,QWORD PTR [r8+rdx*4]
 vmovd  r11d,xmm4
 vpaddd xmm7,xmm4,xmm0
 vpshufd xmm5,xmm4,0x39
 vpshufd xmm8,xmm7,0x39
 movsxd r11,r11d
 vmovq  xmm13,QWORD PTR [r8+rdx*4+0x8]
 vpshufd xmm14,xmm13,0x39
 vpaddd xmm16,xmm13,xmm0
 vpshufd xmm17,xmm16,0x39
 vmovsd xmm6,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm5
 vmovq  xmm22,QWORD PTR [r8+rdx*4+0x10]
 vpshufd xmm23,xmm22,0x39
 vpaddd xmm25,xmm22,xmm0
 movsxd r11,r11d
 vpshufd xmm26,xmm25,0x39
 vmovq  xmm31,QWORD PTR [r8+rdx*4+0x18]
 vpshufd xmm3,xmm31,0x39
 vpaddd xmm2,xmm31,xmm0
 vmovhpd xmm12,xmm6,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm7
 vpshufd xmm1,xmm2,0x39
 movsxd r11,r11d
 vmovsd xmm9,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm8
 movsxd r11,r11d
 vmovhpd xmm10,xmm9,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm13
 vsubpd xmm11,xmm10,xmm12
 movsxd r11,r11d
 vfmadd231pd xmm12,xmm11,XMMWORD PTR [r9+rdx*8]
 vmovsd xmm15,QWORD PTR [rsi+r11*8-0x8]
 vmovupd XMMWORD PTR [rcx+rdx*8],xmm12
 vmovd  r11d,xmm14
 movsxd r11,r11d
 vmovhpd xmm21,xmm15,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm16
 movsxd r11,r11d
 vmovsd xmm18,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm17
 movsxd r11,r11d
 vmovhpd xmm19,xmm18,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm22
 vsubpd xmm20,xmm19,xmm21
 movsxd r11,r11d
 nop
 vfmadd231pd xmm21,xmm20,XMMWORD PTR [r9+rdx*8+0x10]
 vmovsd xmm24,QWORD PTR [rsi+r11*8-0x8]
 vmovupd XMMWORD PTR [rcx+rdx*8+0x10],xmm21
 vmovd  r11d,xmm23
 movsxd r11,r11d
 vmovhpd xmm30,xmm24,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm25
 movsxd r11,r11d
 vmovsd xmm27,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm26
 movsxd r11,r11d
 vmovhpd xmm28,xmm27,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm31
 vsubpd xmm29,xmm28,xmm30
 movsxd r11,r11d
 vfmadd231pd xmm30,xmm29,XMMWORD PTR [r9+rdx*8+0x20]
 vmovsd xmm4,QWORD PTR [rsi+r11*8-0x8]
 vmovupd XMMWORD PTR [rcx+rdx*8+0x20],xmm30
 vmovd  r11d,xmm3
 movsxd r11,r11d
 vmovhpd xmm5,xmm4,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm2
 movsxd r11,r11d
 vmovsd xmm2,QWORD PTR [rsi+r11*8-0x8]
 vmovd  r11d,xmm1
 movsxd r11,r11d
 vmovhpd xmm1,xmm2,QWORD PTR [rsi+r11*8-0x8]
 vsubpd xmm3,xmm1,xmm5
 vfmadd231pd xmm5,xmm3,XMMWORD PTR [r9+rdx*8+0x30]
 vmovupd XMMWORD PTR [rcx+rdx*8+0x30],xmm5
 add    rdx,0x8
 cmp    rdx,rdi
 jb     4038c0 <lint1_+0x40>

  !========== Forced scalar version ==================!
  Loop only

 movsxd r11,DWORD PTR [rbp+r8*1+0x0]
 inc    rax
 vmovsd xmm1,QWORD PTR [r10+r11*8-0x8]
 vmovsd xmm0,QWORD PTR [r10+r11*8]
 movsxd r11,DWORD PTR [rbp+r8*1+0x4]
 vsubsd xmm2,xmm0,xmm1
 vmovsd xmm4,QWORD PTR [r10+r11*8-0x8]
 vmovsd xmm3,QWORD PTR [r10+r11*8]
 movsxd r11,DWORD PTR [rbp+r8*1+0x8]
 vsubsd xmm5,xmm3,xmm4
 vfmadd132sd xmm2,xmm1,QWORD PTR [rdx+r9*1]
 vmovsd xmm7,QWORD PTR [r10+r11*8-0x8]
 vmovsd xmm6,QWORD PTR [r10+r11*8]
 movsxd r11,DWORD PTR [rbp+r8*1+0xc]
 add    rbp,0x10
 vfmadd132sd xmm5,xmm4,QWORD PTR [rdx+r9*1+0x8]
 vsubsd xmm8,xmm6,xmm7
 vmovsd xmm10,QWORD PTR [r10+r11*8-0x8]
 vmovsd xmm9,QWORD PTR [r10+r11*8]
 vfmadd132sd xmm8,xmm7,QWORD PTR [rdx+r9*1+0x10]
 vsubsd xmm11,xmm9,xmm10
 vfmadd132sd xmm11,xmm10,QWORD PTR [rdx+r9*1+0x18]
 vmovsd QWORD PTR [rdx+rcx*1],xmm2
 vmovsd QWORD PTR [rdx+rcx*1+0x8],xmm5
 vmovsd QWORD PTR [rdx+rcx*1+0x10],xmm8
 vmovsd QWORD PTR [rdx+rcx*1+0x18],xmm11
 add    rdx,0x20
 cmp    rax,rsi
 jb     4038c0 <lint1_+0x40>

 !========= using qopt_zmm_usage=high ==================!

 
 vmovdqu ymm0,YMMWORD PTR [r11+rdi*4]
 vmovdqu ymm4,YMMWORD PTR [r11+rdi*4+0x20]
 vmovdqu ymm8,YMMWORD PTR [r11+rdi*4+0x40]
 vmovdqu ymm12,YMMWORD PTR [r11+rdi*4+0x60]
 kxnorw k1,k1,k1
 vpcmpeqb k2,xmm0,xmm0
 vpcmpeqb k3,xmm0,xmm0
 vpcmpeqb k4,xmm0,xmm0
 vpcmpeqb k5,xmm0,xmm0
 vpcmpeqb k6,xmm0,xmm0
 vpcmpeqb k7,xmm0,xmm0
 vpxord zmm1,zmm1,zmm1
 vgatherdpd zmm1{k1},QWORD PTR [r14+ymm0*8]
 vpcmpeqb k1,xmm0,xmm0
 vpxord zmm2,zmm2,zmm2
 vpxord zmm5,zmm5,zmm5
 vpxord zmm6,zmm6,zmm6
 vpxord zmm9,zmm9,zmm9
 vpxord zmm10,zmm10,zmm10
 vpxord zmm13,zmm13,zmm13
 vpxord zmm14,zmm14,zmm14
 vgatherdpd zmm2{k2},QWORD PTR [r14+ymm0*8-0x8]
 vsubpd zmm3,zmm1,zmm2
 vfmadd132pd zmm3,zmm2,ZMMWORD PTR [r10+rdi*8]
 vmovupd ZMMWORD PTR [r15+rdi*8],zmm3
 vgatherdpd zmm5{k3},QWORD PTR [r14+ymm4*8]
 vgatherdpd zmm6{k4},QWORD PTR [r14+ymm4*8-0x8]
 vsubpd zmm7,zmm5,zmm6
 vfmadd132pd zmm7,zmm6,ZMMWORD PTR [r10+rdi*8+0x40]
 vmovupd ZMMWORD PTR [r15+rdi*8+0x40],zmm7
 vgatherdpd zmm9{k5},QWORD PTR [r14+ymm8*8]
 vgatherdpd zmm10{k6},QWORD PTR [r14+ymm8*8-0x8]
 vsubpd zmm11,zmm9,zmm10
 vfmadd132pd zmm11,zmm10,ZMMWORD PTR [r10+rdi*8+0x80]
 vmovupd ZMMWORD PTR [r15+rdi*8+0x80],zmm11
 vgatherdpd zmm13{k7},QWORD PTR [r14+ymm12*8]
 vgatherdpd zmm14{k1},QWORD PTR [r14+ymm12*8-0x8]
 vsubpd zmm15,zmm13,zmm14
 vfmadd132pd zmm15,zmm14,ZMMWORD PTR [r10+rdi*8+0xc0]
 vmovupd ZMMWORD PTR [r15+rdi*8+0xc0],zmm15
 add    rdi,0x20
 cmp    rdi,rbx
 jb     4038c0 <lint1_+0x40>


!==================================================================!
!   More complex loop example

 do ii=1,mx
        i = ix(ii)
       dxm(ii) = (xx(ii)-x(i))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2)) /((x(i-1)-x(i))*(x(i-1)-x(i+1))*(x(i-1)-x(i+2)))
        dx(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2))/((x(i)-x(i-1))*(x(i)-x(i+1))*(x(i)-x(i+2)))
       dxp(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+2)) /((x(i+1)-x(i-1))*(x(i+1)-x(i))*(x(i+1)-x(i+2)))
        dxpp(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+1))/((x(i+2)-x(i-1))*(x(i+2)-x(i))*(x(i+2)-x(i+1)))
       !t0         = xx(ii)
       !t1         = x(i)
       !t2         = x(i-1)
       !t3         = x(i+1)
      ! t4         = x(i+2)
      ! c0         = t0-t1
      ! c1         = t0-t2
      ! c2         = t0-t3
     !  c3         = t0-t4
       !c4         = t2-t1
      ! dxm(ii)     = (c0*c1*c3)/((t2-t1))*(t2-t3)*(t2-t4)
      ! dx(ii)      = (c1*c2*c3)/((t1-t2))*(t1-t3)*(t1-t4)
      ! dxp(ii)     = (c1*c0*c3)/((t3-t2))*(t2-t1)*(t3-t4)
      ! dxpp(ii)    = (c1*c0*c2)/((t4-t2))*(t4-t1)*(t4-t3)   
      ! dxm(ii)    = ((t0-t1)*(t0-t2)*(t0-t4))/((t2-t1))*(t2-t3)*(t2-t4)
      ! dx(ii)     = ((t0-t2)*(t0-t3)*(t0-t4))/((t1-t2))*(t1-t3)*(t1-t4)
      ! dxp(ii)    = ((t0-t2)*(t0-t1)*(t0-t4))/((t3-t2))*(t2-t1)*(t3-t4)
      ! dxpp(ii)   = ((t0-t2)*(t0-t1)*(t0-t3))/((t4-t2))*(t4-t1)*(t4-t3)   
    end do

 vmovdqu ymm2,YMMWORD PTR [r9+rdx*4]
 vmovups zmm1,ZMMWORD PTR [rcx+rdx*8]
 kxnorw k2,k2,k2
 kxnorw k3,k3,k3
 kxnorw k4,k4,k4
 kxnorw k7,k7,k7
 kxnorw k6,k6,k6
 vpaddd ymm0,ymm2,YMMWORD PTR [rip+0x96697]        # 49a040 <__NLITPACK_0.0.2+0x20>
 vpxord zmm8,zmm8,zmm8
 vpxord zmm7,zmm7,zmm7
 vpxord zmm11,zmm11,zmm11
 vpcmpeqb k1,xmm0,xmm0
 vgatherdpd zmm8{k2},QWORD PTR [r11+ymm0*8+0x8]
 vgatherdpd zmm7{k3},QWORD PTR [r11+ymm0*8]
 vgatherdpd zmm11{k4},QWORD PTR [r11+ymm0*8-0x8]
 vsubpd zmm4,zmm1,zmm8
 vsubpd zmm3,zmm1,zmm7
 vsubpd zmm9,zmm11,zmm7
 vsubpd zmm10,zmm11,zmm8
 vmulpd zmm5,zmm3,zmm4
 vmulpd zmm13,zmm9,zmm10
 kxnorw k2,k2,k2
 kxnorw k4,k4,k4
 vpxord zmm12,zmm12,zmm12
 vpxord zmm25,zmm25,zmm25
 vpxord zmm28,zmm28,zmm28
 vpxord zmm24,zmm24,zmm24
 vpxord zmm29,zmm29,zmm29
 vpxord zmm10,zmm10,zmm10
 vpxord zmm9,zmm9,zmm9
 vgatherdpd zmm12{k1},QWORD PTR [r11+ymm0*8+0x10]
 kxnorw k1,k1,k1
 vsubpd zmm14,zmm11,zmm12
 vsubpd zmm6,zmm1,zmm12
 vmulpd zmm15,zmm13,zmm14
 vmulpd zmm17,zmm5,zmm6
 vrcp14pd zmm18,zmm15
 vfnmadd213pd zmm15,zmm18,QWORD PTR [rip+0x9661e]{1to8}        # 49a080 <__NLITPACK_0.0.2+0x60>
 vfpclasspd k0,zmm18,0x1e
 vmulpd zmm16,zmm15,zmm15
 knotw  k5,k0
 vfmadd213pd zmm18{k5},zmm15,zmm18
 vfmadd213pd zmm18{k5},zmm16,zmm18
 kxnorw k5,k5,k5
 vmulpd zmm19,zmm17,zmm18
 vmovupd ZMMWORD PTR [r10+rdx*8],zmm19
 vgatherdpd zmm29{k6},QWORD PTR [r11+ymm0*8+0x10]
 vgatherdpd zmm25{k7},QWORD PTR [r11+ymm0*8+0x8]
 vgatherdpd zmm28{k1},QWORD PTR [r11+ymm0*8]
 vgatherdpd zmm24{k2},QWORD PTR [r11+ymm0*8-0x8]
 vpcmpeqb k2,xmm0,xmm0
 vsubpd zmm27,zmm28,zmm25
 vsubpd zmm26,zmm28,zmm24
 vsubpd zmm31,zmm28,zmm29
 vsubpd zmm20,zmm1,zmm24
 vsubpd zmm21,zmm1,zmm25
 vsubpd zmm23,zmm1,zmm29
 vmulpd zmm30,zmm26,zmm27
 vmulpd zmm22,zmm20,zmm21
 kxnorw k6,k6,k6
 kxnorw k1,k1,k1
 vmulpd zmm24,zmm30,zmm31
 vmulpd zmm2,zmm22,zmm23
 vrcp14pd zmm3,zmm24
 vfnmadd213pd zmm24,zmm3,QWORD PTR [rip+0x96577]{1to8}        # 49a080 <__NLITPACK_0.0.2+0x60>
 vfpclasspd k0,zmm3,0x1e
 vmulpd zmm25,zmm24,zmm24
 knotw  k3,k0
 vfmadd213pd zmm3{k3},zmm24,zmm3
 vfmadd213pd zmm3{k3},zmm25,zmm3
 vpcmpeqb k3,xmm0,xmm0
 vmulpd zmm4,zmm2,zmm3
 vmovupd ZMMWORD PTR [r13+rdx*8+0x0],zmm4
 vgatherdpd zmm10{k6},QWORD PTR [r11+ymm0*8]
 vgatherdpd zmm9{k1},QWORD PTR [r11+ymm0*8-0x8]
 vsubpd zmm6,zmm1,zmm10
 vsubpd zmm5,zmm1,zmm9
 vmulpd zmm7,zmm5,zmm6
 vpxord zmm13,zmm13,zmm13
 vpxord zmm14,zmm14,zmm14
 vpxord zmm30,zmm30,zmm30
 vpxord zmm27,zmm27,zmm27
 vpxord zmm26,zmm26,zmm26
 vpxord zmm2,zmm2,zmm2
 vgatherdpd zmm13{k5},QWORD PTR [r11+ymm0*8+0x8]
 vgatherdpd zmm14{k4},QWORD PTR [r11+ymm0*8+0x10]
 vpcmpeqb k4,xmm0,xmm0
 vpcmpeqb k5,xmm0,xmm0
 vsubpd zmm11,zmm13,zmm9
 vsubpd zmm12,zmm13,zmm10
 vsubpd zmm16,zmm13,zmm14
 vsubpd zmm8,zmm1,zmm14
 vmulpd zmm15,zmm11,zmm12
 vmulpd zmm19,zmm7,zmm8
 vmulpd zmm17,zmm15,zmm16
 vrcp14pd zmm20,zmm17
 vfnmadd213pd zmm17,zmm20,QWORD PTR [rip+0x964ab]{1to8}        # 49a080 <__NLITPACK_0.0.2+0x60>
 vfpclasspd k7,zmm20,0x1e
 vmulpd zmm18,zmm17,zmm17
 knotw  k1,k7
 vfmadd213pd zmm20{k1},zmm17,zmm20
 vfmadd213pd zmm20{k1},zmm18,zmm20
 vmulpd zmm21,zmm19,zmm20
 vmovupd ZMMWORD PTR [r12+rdx*8],zmm21
 nop
 vgatherdpd zmm30{k2},QWORD PTR [r11+ymm0*8+0x10]
 vgatherdpd zmm2{k3},QWORD PTR [r11+ymm0*8+0x8]
 vgatherdpd zmm27{k4},QWORD PTR [r11+ymm0*8]
 vgatherdpd zmm26{k5},QWORD PTR [r11+ymm0*8-0x8]
 vsubpd zmm29,zmm30,zmm27
 vsubpd zmm28,zmm30,zmm26
 vsubpd zmm0,zmm1,zmm26
 vsubpd zmm22,zmm1,zmm27
 vsubpd zmm1,zmm1,zmm2
 vsubpd zmm2,zmm30,zmm2
 vmulpd zmm31,zmm28,zmm29
 vmulpd zmm23,zmm0,zmm22
 vmulpd zmm26,zmm31,zmm2
 vmulpd zmm0,zmm23,zmm1
 vrcp14pd zmm1,zmm26
 vfnmadd213pd zmm26,zmm1,QWORD PTR [rip+0x96415]{1to8}        # 49a080 <__NLITPACK_0.0.2+0x60>
 vfpclasspd k0,zmm1,0x1e
 vmulpd zmm27,zmm26,zmm26
 knotw  k6,k0
 vfmadd213pd zmm1{k6},zmm26,zmm1
 vfmadd213pd zmm1{k6},zmm27,zmm1
 vmulpd zmm3,zmm0,zmm1
 vmovupd ZMMWORD PTR [rdi+rdx*8],zmm3
 add    rdx,0x8
 cmp    rdx,rax
 jb     403980 <cubnmx_+0x100>


!======================================================================================!
  do ii=1,mx
        i = ix(ii)
       !dxm(ii) = (xx(ii)-x(i))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2)) /((x(i-1)-x(i))*(x(i-1)-x(i+1))*(x(i-1)-x(i+2)))
       ! dx(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2))/((x(i)-x(i-1))*(x(i)-x(i+1))*(x(i)-x(i+2)))
       !dxp(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+2)) /((x(i+1)-x(i-1))*(x(i+1)-x(i))*(x(i+1)-x(i+2)))
      !  dxpp(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+1))/((x(i+2)-x(i-1))*(x(i+2)-x(i))*(x(i+2)-x(i+1)))
       t0         = xx(ii)
       t1         = x(i)
       t2         = x(i-1)
       t3         = x(i+1)
       t4         = x(i+2)
      ! c0         = t0-t1
      ! c1         = t0-t2
      ! c2         = t0-t3
     !  c3         = t0-t4
       !c4         = t2-t1
      ! dxm(ii)     = (c0*c1*c3)/((t2-t1))*(t2-t3)*(t2-t4)
      ! dx(ii)      = (c1*c2*c3)/((t1-t2))*(t1-t3)*(t1-t4)
      ! dxp(ii)     = (c1*c0*c3)/((t3-t2))*(t2-t1)*(t3-t4)
      ! dxpp(ii)    = (c1*c0*c2)/((t4-t2))*(t4-t1)*(t4-t3)   
       dxm(ii)    = ((t0-t1)*(t0-t2)*(t0-t4))/((t2-t1))*(t2-t3)*(t2-t4)
       dx(ii)     = ((t0-t2)*(t0-t3)*(t0-t4))/((t1-t2))*(t1-t3)*(t1-t4)
       dxp(ii)    = ((t0-t2)*(t0-t1)*(t0-t4))/((t3-t2))*(t2-t1)*(t3-t4)
       dxpp(ii)   = ((t0-t2)*(t0-t1)*(t0-t3))/((t4-t2))*(t4-t1)*(t4-t3)   
    end do


 vmovdqu ymm4,YMMWORD PTR [r9+rdx*4]
 vmovups zmm14,ZMMWORD PTR [rcx+rdx*8]
 vpcmpeqb k4,xmm0,xmm0
 vpcmpeqb k2,xmm0,xmm0
 vpaddd ymm5,ymm4,YMMWORD PTR [rip+0x9569f]        # 499040 <__NLITPACK_0.0.2+0x20>
 kxnorw k3,k3,k3
 kxnorw k1,k1,k1
 vpxord zmm24,zmm24,zmm24
 vpxord zmm11,zmm11,zmm11
 vpxord zmm3,zmm3,zmm3
 vpxord zmm2,zmm2,zmm2
 vgatherdpd zmm24{k3},QWORD PTR [r11+ymm5*8]
 vgatherdpd zmm11{k4},QWORD PTR [r11+ymm5*8-0x8]
 vgatherdpd zmm3{k1},QWORD PTR [r11+ymm5*8+0x10]
 vgatherdpd zmm2{k2},QWORD PTR [r11+ymm5*8+0x8]
 vsubpd zmm0,zmm11,zmm24
 vsubpd zmm6,zmm14,zmm24
 vsubpd zmm15,zmm14,zmm11
 vsubpd zmm17,zmm14,zmm3
 vsubpd zmm27,zmm11,zmm2
 vxorpd zmm18,zmm0,QWORD PTR [rip+0x95678]{1to8}        # 499080 <__NLITPACK_0.0.2+0x60>
 vsubpd zmm4,zmm11,zmm3
 vsubpd zmm22,zmm24,zmm2
 vsubpd zmm5,zmm24,zmm3
 vxorpd zmm28,zmm27,QWORD PTR [rip+0x9565c]{1to8}        # 499080 <__NLITPACK_0.0.2+0x60>
 vrcp14pd zmm9,zmm0
 vrcp14pd zmm21,zmm18
 vrcp14pd zmm31,zmm28
 vmulpd zmm1,zmm6,zmm15
 vfnmadd213pd zmm18,zmm21,QWORD PTR [rip+0x95642]{1to8}        # 499088 <__NLITPACK_0.0.2+0x68>
 vfnmadd213pd zmm28,zmm31,QWORD PTR [rip+0x95638]{1to8}        # 499088 <__NLITPACK_0.0.2+0x68>
 vfpclasspd k0,zmm9,0x1e
 vfpclasspd k6,zmm21,0x1e
 vmulpd zmm30,zmm1,zmm17
 vmulpd zmm19,zmm18,zmm18
 vmulpd zmm29,zmm28,zmm28
 knotw  k5,k0
 knotw  k7,k6
 vfpclasspd k0,zmm31,0x1e
 vfmadd213pd zmm21{k7},zmm18,zmm21
 knotw  k1,k0
 vfmadd213pd zmm21{k7},zmm19,zmm21
 vfmadd213pd zmm31{k1},zmm28,zmm31
 vmovaps zmm7,zmm0
 vfnmadd213pd zmm7,zmm9,QWORD PTR [rip+0x955e3]{1to8}        # 499088 <__NLITPACK_0.0.2+0x68>
 vfmadd213pd zmm31{k1},zmm29,zmm31
 vmulpd zmm8,zmm7,zmm7
 vfmadd213pd zmm9{k5},zmm7,zmm9
 vxorpd zmm7,zmm5,QWORD PTR [rip+0x955bf]{1to8}        # 499080 <__NLITPACK_0.0.2+0x60>
 vfmadd213pd zmm9{k5},zmm8,zmm9
 vsubpd zmm8,zmm2,zmm3
 vmulpd zmm10,zmm30,zmm9
 vmulpd zmm12,zmm27,zmm10
 vxorpd zmm10,zmm8,QWORD PTR [rip+0x9559d]{1to8}        # 499080 <__NLITPACK_0.0.2+0x60>
 vmulpd zmm13,zmm12,zmm4
 vsubpd zmm12,zmm14,zmm2
 vmovupd ZMMWORD PTR [r10+rdx*8],zmm13
 vmulpd zmm16,zmm15,zmm12
 vmulpd zmm20,zmm16,zmm17
 vmulpd zmm23,zmm20,zmm21
 vmulpd zmm25,zmm22,zmm23
 vmulpd zmm26,zmm25,zmm5
 vmulpd zmm25,zmm30,zmm31
 vmovupd ZMMWORD PTR [r13+rdx*8+0x0],zmm26
 vmulpd zmm26,zmm1,zmm12
 vmulpd zmm0,zmm0,zmm25
 vxorpd zmm1,zmm4,QWORD PTR [rip+0x95548]{1to8}        # 499080 <__NLITPACK_0.0.2+0x60>
 vmulpd zmm2,zmm0,zmm8
 vrcp14pd zmm27,zmm1
 vmovupd ZMMWORD PTR [r12+rdx*8],zmm2
 vfnmadd213pd zmm1,zmm27,QWORD PTR [rip+0x95533]{1to8}        # 499088 <__NLITPACK_0.0.2+0x68>
 vfpclasspd k2,zmm27,0x1e
 vmulpd zmm3,zmm1,zmm1
 knotw  k3,k2
 vfmadd213pd zmm27{k3},zmm1,zmm27
 vfmadd213pd zmm27{k3},zmm3,zmm27
 vmulpd zmm6,zmm26,zmm27
 vmulpd zmm9,zmm6,zmm7
 vmulpd zmm11,zmm9,zmm10
 vmovupd ZMMWORD PTR [rdi+rdx*8],zmm11
 add    rdx,0x8
 cmp    rdx,rax
 jb     403980 <cubnmx_+0x100>



#endif

    contains
!***************************************************************************************************

!**************************************************************************
!>
!  Wrapper to rgrd1.  Allocates the work arrays internally.

    subroutine rgrd1_wrapper(x,p,xx,q,intpol,ier)
        !dir$ attributes code_align : 32 :: rgrd1_wrapper
        !dir$ optimize : 3
    implicit none

    real(kind=dp),dimension(:),intent(in)     :: x            !! original x
    real(kind=dp),dimension(:),intent(in)     :: p            !! original p(x)
    real(kind=dp),dimension(:),intent(in)     :: xx           !! regridded xx
    real(kind=dp),dimension(:),intent(out)    :: q            !! regridded q(xx)
    integer(kind=i4),intent(in)                   :: intpol
    integer(kind=i4),intent(out)                  :: ier          !! status code:
                                                         !!
                                                         !! * 0    : no errors
                                                         !! * 1-6 : error [see original code]
                                                         !! * 10  : input vectors are the wrong size
                                                         !! * 100 : out of memory

    integer(kind=i4) :: lw, liw
    integer(kind=i4) :: nx, mx
    integer(kind=i4) :: np, nq
    real(kind=dp),dimension(:),allocatable :: w
    integer(kind=i4),dimension(:),allocatable :: iw
    integer(kind=i4) :: ierr1, ierr2

    !get array sizes:

    nx = size(x)
    np = size(p)

    mx = size(xx)
    nq = size(q)

    if (nx/=np .or. mx/=nq) then
        !Error: vectors are the wrong size
        ier = 10
        return
    end if

    !allocate work matrices:

    select case(intpol)
    case(1)
        lw = mx
    case(3)
        lw = 4*mx
    case default
        ier = 6     !Error: invalid intpol value
        return
    end select

    liw = mx

    allocate(w(lw),   stat=ierr1)
    allocate(iw(liw), stat=ierr2)

    if (ierr1==0 .and. ierr2==0) then
        !call the main routine:
        call rgrd1(nx,x,p,mx,xx,q,intpol,w,lw,iw,liw,ier)
    else
        !error: out of memory
        ier = 100
    end if

    !clean up:
    if (allocated(w)) deallocate(w)
    if (allocated(iw)) deallocate(iw)

    end subroutine rgrd1_wrapper
!**************************************************************************

!**************************************************************************
!>
!  Wrapper to rgrd2.  Allocates the work arrays internally.

    subroutine rgrd2_wrapper(x,y,p,xx,yy,q,intpol,ier)
        !dir$ attributes code_align : 32 :: rgrd2_wrapper
        !dir$ optimize : 3
 
    implicit none

    real(kind=dp),dimension(:),intent(in)     :: x              !! original x
    real(kind=dp),dimension(:),intent(in)     :: y              !! original y
    real(kind=dp),dimension(:,:),intent(in)   :: p              !! original p(x,y)
    real(kind=dp),dimension(:),intent(in)     :: xx             !! regridded xx
    real(kind=dp),dimension(:),intent(in)     :: yy             !! regridded yy
    real(kind=dp),dimension(:,:),intent(out)  :: q              !! regridded q(xx,yy)
    integer(kind=i4),dimension(2),intent(in)      :: intpol
    integer(kind=i4),intent(out)                  :: ier            !! * 0    : no errors
                                                           !! * 1-6 : error [see original code]
                                                           !! * 10  : input vectors are the wrong size
                                                           !! * 100 : out of memory

    integer(kind=i4) :: lw, liw
    integer(kind=i4) :: nx, ny, mx, my
    integer(kind=i4),dimension(2) :: np, nq
    integer(kind=i4) :: lwx, lwy
    real(kind=dp),dimension(:),allocatable :: w
    integer(kind=i4),dimension(:),allocatable :: iw
    integer(kind=i4) :: ierr1, ierr2

    !get array sizes:

    nx = size(x)
    ny = size(y)
    np(1) = size(p,1)
    np(2) = size(p,2)

    mx = size(xx)
    my = size(yy)
    nq(1) = size(q,1)
    nq(2) = size(q,2)

    if (nx/=np(1) .or. ny/=np(2) .or. mx/=nq(1) .or. my/=nq(2)) then

        !Error: vectors are the wrong size
        ier = 10
        return

    end if

    !allocate work matrices:

    select case(intpol(1))
    case(1)
        lwx = mx
    case(3)
        lwx = 4*mx
        case default
        ier = 6     !Error: invalid intpol value
        return
    end select

    select case(intpol(2))
    case(1)
        lwy = my+2*mx
    case(3)
        lwy = 4*(mx+my)
    end select

    lw  = lwx + lwy
    liw = mx + my

    allocate(w(lw),   stat=ierr1)
    allocate(iw(liw), stat=ierr2)

    if (ierr1==0 .and. ierr2==0) then

        !call the main routine:
        call rgrd2(nx,ny,x,y,p,mx,my,xx,yy,q,intpol,w,lw,iw,liw,ier)

    else

        !error: out of memory
        ier = 100

    end if

    !clean up:

    deallocate(w)
    deallocate(iw)

    end subroutine rgrd2_wrapper
!**************************************************************************

!**************************************************************************
!>
!  Wrapper to rgrd3.  Allocates the work arrays internally.

    subroutine rgrd3_wrapper(x,y,z,p,xx,yy,zz,q,intpol,ier)
        !dir$ attributes code_align : 32 :: rgrd3_wrapper
        !dir$ optimize : 3
    implicit none

    real(kind=dp),dimension(:),intent(in)              :: x            !! original x
    real(kind=dp),dimension(:),intent(in)              :: y            !! original y
    real(kind=dp),dimension(:),intent(in)              :: z            !! original z
    real(kind=dp),dimension(:,:,:),intent(in)          :: p            !! original p(x,y,z)
    real(kind=dp),dimension(:),intent(in)              :: xx           !! regridded xx
    real(kind=dp),dimension(:),intent(in)              :: yy           !! regridded yy
    real(kind=dp),dimension(:),intent(in)              :: zz           !! regridded zz
    real(kind=dp),dimension(:,:,:),intent(out)         :: q            !! regridded q(xx,yy,zz)
    integer(kind=i4),dimension(3),intent(in)               :: intpol
    integer(kind=i4),intent(out)                           :: ier          !! * 0   : no errors
                                                                  !! * 1-6 : error [see original code]
                                                                  !! * 10  : input vectors are the wrong size
                                                                  !! * 100 : out of memory

    integer(kind=i4) :: nx, ny, nz, mx, my, mz
    integer(kind=i4),dimension(3) :: np, nq
    integer(kind=i4) :: lw, liw
    integer(kind=i4) :: lwx, lwy, lwz
    real(kind=dp),dimension(:),allocatable :: w
    integer(kind=i4),dimension(:),allocatable :: iw
    integer(kind=i4) :: ierr1, ierr2

    !get array sizes:

    nx = size(x)
    ny = size(y)
    nz = size(z)
    np(1) = size(p,1)
    np(2) = size(p,2)
    np(3) = size(p,3)

    mx = size(xx)
    my = size(yy)
    mz = size(zz)
    nq(1) = size(q,1)
    nq(2) = size(q,2)
    nq(3) = size(q,3)

    if (nx/=np(1) .or. ny/=np(2) .or. nz/=np(3) .or. mx/=nq(1) .or. my/=nq(2) .or. mz/=nq(3)) then

        !Error: vectors are the wrong size
        ier = 10
        return

    end if

    !allocate work matrices:

    select case(intpol(1))
    case(1)
        lwx = mx
    case(3)
        lwx = 4*mx
        case default
        ier = 6     !Error: invalid intpol value
        return
    end select

    select case(intpol(2))
    case(1)
        lwy = my+2*mx
    case(3)
        lwy = 4*(mx+my)
    end select
    select case(intpol(3))
    case(1)
        lwz = 2*mx*my+mz
    case(3)
        lwz = 4*(mx*my+mz)
    end select

    lw  = lwx + lwy + lwz
    liw = mx + my + mz

    allocate(w(lw),   stat=ierr1)
    allocate(iw(liw), stat=ierr2)

    if (ierr1==0 .and. ierr2==0) then

        !call the main routine:
        call rgrd3(nx,ny,nz,x,y,z,p,mx,my,mz,xx,yy,zz,q,intpol,w,lw,iw,liw,ier)

    else

        !error: out of memory
        ier = 100

    end if

    !clean up:

    deallocate(w)
    deallocate(iw)

    end subroutine rgrd3_wrapper
!**************************************************************************

!**************************************************************************
!>
!  Wrapper to rgrd4.  Allocates the work arrays internally.

    subroutine rgrd4_wrapper(x,y,z,t,p,xx,yy,zz,tt,q,intpol,ier)
        !dir$ attributes code_align : 32 :: rgrd4_wrapper
        !dir$ optimize : 3
    implicit none

    real(kind=dp),dimension(:),intent(in)          :: x            !! original x
    real(kind=dp),dimension(:),intent(in)          :: y            !! original y
    real(kind=dp),dimension(:),intent(in)          :: z            !! original z
    real(kind=dp),dimension(:),intent(in)          :: t            !! original t
    real(kind=dp),dimension(:,:,:,:),intent(in)    :: p            !! original p(x,y,z,t)
    real(kind=dp),dimension(:),intent(in)          :: xx           !! regridded xx
    real(kind=dp),dimension(:),intent(in)          :: yy           !! regridded yy
    real(kind=dp),dimension(:),intent(in)          :: zz           !! regridded zz
    real(kind=dp),dimension(:),intent(in)          :: tt           !! regridded tt
    real(kind=dp),dimension(:,:,:,:),intent(out)   :: q            !! regridded q(xx,yy,zz,tt)
    integer(kind=i4),dimension(4),intent(in)           :: intpol
    integer(kind=i4),intent(out)                       :: ier          !! * 0    : no errors
                                                              !! * 1-6 : error [see original code]
                                                              !! * 10  : input vectors are the wrong size
                                                              !! * 100 : out of memory

    integer(kind=i4) :: nx, ny, nz, nt, mx, my, mz, mt
    integer(kind=i4),dimension(4) :: np, nq
    integer(kind=i4) :: lw, liw
    integer(kind=i4) :: lwx, lwy, lwz, lwt
    real(kind=dp),dimension(:),allocatable :: w
    integer(kind=i4),dimension(:),allocatable :: iw
    integer(kind=i4) :: ierr1, ierr2

    !get array sizes:

    nx = size(x)
    ny = size(y)
    nz = size(z)
    nt = size(t)

    np(1) = size(p,1)
    np(2) = size(p,2)
    np(3) = size(p,3)
    np(4) = size(p,4)

    mx = size(xx)
    my = size(yy)
    mz = size(zz)
    mt = size(tt)

    nq(1) = size(q,1)
    nq(2) = size(q,2)
    nq(3) = size(q,3)
    nq(4) = size(q,4)

    if (nx/=np(1) .or. ny/=np(2) .or. nz/=np(3) .or. nt/=np(4) .or. &
        mx/=nq(1).or. my/=nq(2) .or. mz/=nq(3) .or. mt/=nq(4)) then

        !Error: vectors are the wrong size
        ier = 10
        return

    end if

    !allocate work matrices:

    select case(intpol(1))
    case(1)
        lwx = mx
    case(3)
        lwx = 4*mx
        case default
        ier = 6     !Error: invalid intpol value
        return
    end select

    select case(intpol(2))
    case(1)
        lwy = my+2*mx
    case(3)
        lwy = 4*(mx+my)
    end select

    select case(intpol(3))
    case(1)
        lwz = 2*mx*my+mz
    case(3)
        lwz = 4*(mx*my+mz)
    end select

    select case(intpol(4))
    case(1)
        lwt = 2*mx*my*mz+mt
    case(3)
        lwt = 4*(mx*my*mz+mt)
    end select

    lw  = lwx + lwy + lwz + lwt
    liw = mx + my + mz + mt

    allocate(w(lw),   stat=ierr1)
    allocate(iw(liw), stat=ierr2)

    if (ierr1==0 .and. ierr2==0) then

        !call the main routine:
        call rgrd4(nx,ny,nz,nt,x,y,z,t,p,mx,my,mz,mt,xx,yy,zz,tt,q,intpol,w,lw,iw,liw,ier)

    else

        !error: out of memory
        ier = 100

    end if

    !clean up:

    deallocate(w)
    deallocate(iw)

    end subroutine rgrd4_wrapper
!**************************************************************************

!**************************************************************************
!>
!  subroutine rgrd1 interpolates the values p(i) on the grid x(i)
!  for i=1,...,nx onto q(ii) on the grid xx(ii),ii=1,...,mx.
!
!### requirements
!
!  x must be a strictly increasing grid and xx must be an increasing
!  grid (see ier = 4).  in addition the interval
!
!    [xx(1),xx(mx)]
!
!  must lie within the interval
!
!    [x(1),x(nx)].
!
!  extrapolation is not allowed (see ier=3).  if these intervals
!  are identical and the x and xx grids are UNIFORM then subroutine
!  rgrd1u should be used in place of rgrd1.

    subroutine rgrd1(nx,x,p,mx,xx,q,intpol,w,lw,iw,liw,ier)
        !dir$ attributes code_align : 32 :: rgrd1
        !dir$ optimize : 3
    implicit none

    integer(kind=i4),intent(in) :: nx                       !! the integer(kind=i4) dimension of the grid
                                                   !! vector x and the dimension of p.
                                                   !! nx > 1 if intpol = 1 or nx > 3 if
                                                   !! intpol = 3 is required.
    real(kind=dp),dimension(nx),intent(in) :: x         !! a real(kind=dp) nx vector of strictly
                                                   !! increasing values which defines the x
                                                   !! grid on which p is given.
    real(kind=dp),dimension(nx),intent(in) :: p         !! a real(kind=dp) nx vector of values given on the x grid
    integer(kind=i4),intent(in) :: mx                       !! the integer(kind=i4) dimension of the grid vector
                                                   !! xx and the dimension of q.
                                                   !! mx > 0 is required.
    real(kind=dp),dimension(mx),intent(in)  :: xx       !! a real(kind=dp) mx vector of increasing values which defines the
                                                   !! grid on which q is defined.  xx(1) < x(1) or xx(mx) > x(nx)
                                                   !! is not allowed (see ier = 3)
    integer(kind=i4),intent(in) :: intpol                   !! an integer(kind=i4) which sets linear or cubic
                                                   !! interpolation as follows:
                                                   !!
                                                   !! * intpol = 1 sets linear interpolation
                                                   !! * intpol = 3 sets cubic interpolation
                                                   !!
                                                   !! values other than 1 or 3 in intpol are not allowed (ier = 6).
    real(kind=dp),intent(out) :: q(mx)                  !! a real(kind=dp) mx vector of values on the xx grid which are
                                                   !! interpolated from p on the x grid
    integer(kind=i4),intent(in) :: lw                       !! the integer(kind=i4) length of the real(kind=dp) work space w.  let
                                                   !!
                                                   !!  * lwmin = mx     if intpol(1) = 1
                                                   !!  * lwmin = 4*mx   if intpol(1) = 3
                                                   !!
                                                   !! then lw must be greater than or equal to lwmin
    real(kind=dp),dimension(lw),intent(inout) :: w      !! a real(kind=dp) work space of length at least
                                                   !! lw which must be provided in the
                                                   !! routine calling rgrd1
    integer(kind=i4),intent(in) :: liw                      !! the length of the integer(kind=i4) work space iw.
                                                   !! liw must be greater than or equal to mx.
    integer(kind=i4),dimension(*),intent(inout) :: iw       !! an integer(kind=i4) work space of length at least
                                                   !! liw which must be provided in the
                                                   !! routine calling rgrd1
    integer(kind=i4),intent(out) :: ier                     !! an integer(kind=i4) error flag set as follows:
                                                   !!
                                                   !! * ier = 0 if no errors in input arguments are detected
                                                   !! * ier = 1 if mx < 1
                                                   !! * ier = 2 if nx < 2 when intpol=1 or nx < 4 when intpol=3
                                                   !! * ier = 3 if xx(1) < x(1) or x(nx) < xx(mx)
                                                   !!    to avoid this flag when end points are intended to be the
                                                   !!    same but may differ slightly due to roundoff error, they
                                                   !!    should be set exactly in the calling routine (e.g., if both
                                                   !!    grids have the same x boundaries then xx(1)=x(1) and xx(mx)=x(nx)
                                                   !!    should be set before calling rgrd1)
                                                   !! * ier = 4 if the x grid is not strictly monotonically increasing
                                                   !!    or if the xx grid is not montonically increasing.  more
                                                   !!    precisely if:
                                                   !!    x(i+1) <= x(i) for some i such that 1 <= i < nx (or)
                                                   !!    xx(ii+1) < xx(ii) for some ii such that 1 <= ii < mx
                                                   !! * ier = 5 if lw or liw is too small (insufficient work space)
                                                   !! * ier = 6 if intpol is not equal to 1 or 3

    integer(kind=i4) :: i,ii,i1,i2,i3,i4

    ! check arguments for errors

    ! check xx grid resolution
    ier = 1
    if (mx < 1) return

    ! check intpol
    ier = 6
    if (intpol/=1 .and. intpol/=3) return

    ! check x grid resolution
    ier = 2
    if (intpol==1 .and. nx<2) return
    if (intpol==3 .and. nx<4) return

    ! check xx grid contained in x grid
    ier = 3
    if (xx(1)<x(1) .or. xx(mx)>x(nx)) return

    ! check montonicity of grids
    do i=2,nx
        if (x(i-1)>=x(i)) then
            ier = 4
            return
        end if
    end do
    do ii=2,mx
        if (xx(ii-1)>xx(ii)) then
            ier = 4
            return
        end if
    end do

    ! check minimum work space lengths
    ier = 5
    if (intpol==1) then
        if (lw < mx) return
    else
        if (lw < 4*mx) return
    end if
    if (liw < mx) return

    ! arguments o.k.

    ier = 0

    if (intpol==1) then
        ! linear interpolation in x
        call linmx(nx,x,mx,xx,iw,w)
        call lint1(nx,p,mx,q,iw,w)
    else
        ! cubic interpolation in x
        i1 = 1
        i2 = i1+mx
        i3 = i2+mx
        i4 = i3+mx
        call cubnmx(nx,x,mx,xx,iw,w(i1),w(i2),w(i3),w(i4))
        call cubt1(nx,p,mx,q,iw,w(i1),w(i2),w(i3),w(i4))
    end if

    end subroutine rgrd1
!**************************************************************************

!**************************************************************************
!>
!  linearly interpolate p on x onto q on xx

    subroutine lint1(nx,p,mx,q,ix,dx)
       !dir$ attributes forceinline :: lint1
       !dir$ attributes code_align : 32 :: lint1
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: lint1
    implicit none

    integer(kind=i4) :: mx,ix(mx),nx,ii,i
    real(kind=dp) :: p(nx),q(mx),dx(mx)
    real(kind=dp), automatic :: t0,t1
    !dir$ assume_aligned q:64
    !dir$ assume_aligned p:64
    !dir$ assume_aligned ix:64
    !dir$ assume_aligned dx:64
#if (AUTO_VECTORIZE) == 1
    !dir$ ivdep
    !dir$ vector aligned
    !dir$ vector always
    !dir$ vector vectorlength(8)
    
#else
   
    !dir$ novector
#endif
    do ii=1,mx
        i  = ix(ii)
        t0 = p(i)
        t1 = p(i+1)
        q(ii)  = t0+dx(ii)*(t1-t0)
        !q(ii) = p(i)+dx(ii)*(p(i+1)-p(i))
    end do

    end subroutine lint1
!**************************************************************************

!**************************************************************************
!>
! cubically interpolate p on x to q on xx

    subroutine cubt1(nx,p,mx,q,ix,dxm,dx,dxp,dxpp)
       !dir$ attributes forceinline :: cubt1
       !dir$ attributes code_align : 32 :: cubt1
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: cubt1
    implicit none

    integer(kind=i4) :: mx,ix(mx),nx,i,ii
    real(kind=dp) :: p(nx),q(mx),dxm(mx),dx(mx),dxp(mx),dxpp(mx)
    real(kind=dp), automatic :: t0,t1,t2
    !dir$ assume_aligned q:64
    !dir$ assume_aligned p:64
    !dir$ assume_aligned ix:64
    !dir$ assume_aligned dx:64
    !dir$ assume_aligned dxm:64
    !dir$ assume_aligned dxp:64
    !dir$ assume_aligned dxpp:64
#if (AUTO_VECTORIZE) == 1
    !dir$ ivdep
    !dir$ vector aligned
    !dir$ vector always
    !dir$ vector vectorlength(8)
   
#else
   
    !dir$ novector
#endif
    do ii=1,mx
        i  = ix(ii)
        t0 = p(i-1)
        t1 = p(i)
        t2 = p(i+1)
        t3 = p(i+2)
        q(ii)  = dxm(ii)-t0+dx(ii)*t1+dxp(ii)*t2+dxpp(ii)*t3
        !q(ii) = dxm(ii)*p(i-1)+dx(ii)*p(i)+dxp(ii)*p(i+1)+dxpp(ii)*p(i+2)
    end do

    end subroutine cubt1
!**************************************************************************

!**************************************************************************
!>
!  set x grid pointers for xx grid and interpolation scale terms

    subroutine linmx(nx,x,mx,xx,ix,dx)
       !dir$ attributes forceinline :: linmx
       !dir$ attributes code_align : 32 :: linmx
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: linmx
    implicit none

    real(kind=dp) :: x(*),xx(*),dx(*)
    integer(kind=i4) :: ix(*),isrt,ii,i,nx,mx
    real(kind=dp), automatic :: t0,t1,t2
    isrt = 1

    do ii=1,mx
        ! find x(i) s.t. x(i) < xx(ii) <= x(i+1)
        do i=isrt,nx-1
            if (x(i+1) >= xx(ii)) then
                isrt = i
                ix(ii) = i
                exit
            end if
        end do
    end do

    ! set linear scale term
   
    !dir$ assume_aligned ix:64
    !dir$ assume_aligned dx:64
    !dir$ assume_aligned xx:64
    !dir$ assume_aligned x:64
#if (AUTO_VECTORIZE) == 1
    !dir$ ivdep
    !dir$ vector aligned
    !dir$ vector always
    !dir$ vector vectorlength(8)
   
#else
   
    !dir$ novector
#endif
    do ii=1,mx
        i  = ix(ii)
        t0 = xx(ii)
        t1 = x(i)
        t2 = x(i+1)
        !dx(ii) = (xx(ii)-x(i))/(x(i+1)-x(i))
        dx(ii) = (t0-t1)/(t2-t1)
    end do

    end subroutine linmx
!**************************************************************************

!**************************************************************************
!>
!
    subroutine cubnmx(nx,x,mx,xx,ix,dxm,dx,dxp,dxpp)
       !dir$ attributes forceinline :: cubnmx
       !dir$ attributes code_align : 32 :: cubnmx
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: cubnmx
    implicit none

    real(kind=dp) :: x(*),xx(*),dxm(*),dx(*),dxp(*),dxpp(*)
    integer(kind=i4) :: ix(*),mx,nx,i,ii,isrt
    real(kind=dp), automatic :: t0,t1,t2,t3,t4
    real(kind=dp), automatic :: c0,c1,c2,c3,c4
    isrt = 1

    do ii=1,mx
        ! set i in [2,nx-2] closest s.t.
        ! x(i-1),x(i),x(i+1),x(i+2) can interpolate xx(ii)
        do i=isrt,nx-1
            if (x(i+1) >= xx(ii)) then
                ix(ii) = min(nx-2,max(2,i))
                isrt = ix(ii)
                exit
            end if
        end do
    end do

    ! set cubic scale terms
    !dir$ assume_aligned ix:64
    !dir$ assume_aligned dx:64
    !dir$ assume_aligned xx:64
    !dir$ assume_aligned x:64
    !dir$ assume_aligned dxm:64
    !dir$ assume_aligned dxp:64
    !dir$ assume_aligned dxpp:64
#if (AUTO_VECTORIZE) == 1
    !dir$ ivdep
    !dir$ vector aligned
    !dir$ vector always
    !dir$ vector vectorlength(8)
    !dir$ unroll(8)
#else
    !dir$ unroll(8)
    !dir$ novector
#endif
    do ii=1,mx
        i = ix(ii)
       ! dxm(ii) = (xx(ii)-x(i))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2)) /((x(i-1)-x(i))*(x(i-1)-x(i+1))*(x(i-1)-x(i+2)))
       ! dx(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2))/((x(i)-x(i-1))*(x(i)-x(i+1))*(x(i)-x(i+2)))
       ! dxp(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+2)) /((x(i+1)-x(i-1))*(x(i+1)-x(i))*(x(i+1)-x(i+2)))
       ! dxpp(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+1))/((x(i+2)-x(i-1))*(x(i+2)-x(i))*(x(i+2)-x(i+1)))
       t0         = xx(ii)
       t1         = x(i)
       t2         = x(i-1)
       t3         = x(i+1)
       t4         = x(i+2)
       dxm(ii)    = ((t0-t1)*(t0-t2)*(t0-t4))/((t2-t1))*(t2-t3)*(t2-t4)
       dx(ii)     = ((t0-t2)*(t0-t3)*(t0-t4))/((t1-t2))*(t1-t3)*(t1-t4)
       dxp(ii)    = ((t0-t2)*(t0-t1)*(t0-t4))/((t3-t2))*(t2-t1)*(t3-t4)
       dxpp(ii)   = ((t0-t2)*(t0-t1)*(t0-t3))/((t4-t2))*(t4-t1)*(t4-t3)   
    end do

    end subroutine cubnmx
!**************************************************************************

!**************************************************************************
!>
!  subroutine rgrd1u interpolates the nx vector p onto
!  the mx vector q. it is assumed that p and q are
!  values on uniform nx and mx grids which subdivide
!  the same interval (INCLUDING END POINTS).  if p and
!  q are values on nonuniform grids and/or if q is defined
!  on a grid which lies within the p grid then subroutine
!  rgrd1 should be used.
!
!### method
!
!  linear or cubic interpolation (see intpol) is used when the
!  mx uniform grid is not a subgrid of the nx uniform grid (i.e.,
!  whenever mx-1 does not divide nx-1).  q is set directly from
!  p in the subgrid case.

    subroutine rgrd1u(nx,p,mx,q,intpol,w,lw,iw,liw,ier)
       !dir$ attributes code_align : 32 :: rgrd1u
       !dir$ optimize : 3
    implicit none

    integer(kind=i4),intent(in)  :: intpol       !! an integer(kind=i4) which sets linear or cubic interpolation as follows:
                                        !!
                                        !! * intpol = 1 sets linear interpolation
                                        !! * intpol = 3 sets cubic interpolation
                                        !!
                                        !! values other than 1 or 3 in intpol are not allowed (ier = 4).
    integer(kind=i4),intent(in)  :: liw          !! the integer(kind=i4) length of the integer(kind=i4) work space iw in the routine calling rgrd1u.
                                        !! liw must be greater than or equal to mx.
    integer(kind=i4),intent(inout)  :: iw(liw)   !! an integer(kind=i4) work space of length liw
    integer(kind=i4),intent(in)  :: lw           !! the integer(kind=i4) length of the work space w in the routine calling rgrd1u.
                                        !! if mx-1 divides nx-1 then the mx uniform grid is a subgrid of
                                        !! the nx uniform grid.  in this case let lwmin = 1.  otherwise
                                        !! let lwmin = mx if intpol = 1 or lwmin = mx if intpol = 3.
                                        !! then lw must be greater than or equal to lwmin (see ier=4).
    integer(kind=i4),intent(in) :: nx            !! the integer(kind=i4) dimension of p.  nx > 1 if intpol = 1 or
                                        !! nx > 3 if intpol = 3 is required (see ier = 2).
    integer(kind=i4),intent(in)  :: mx           !! the integer(kind=i4) dimension of q.  mx > 1 is required (see ier = 1)
    real(kind=dp),intent(in) :: p(nx)        !! a real(kind=dp) nx dimensioned vector of given values
    real(kind=dp),intent(inout) :: w(lw)     !! a real(kind=dp) work space of length lw.
    real(kind=dp),dimension(mx),intent(out) :: q  !! a real(kind=dp) mx dimensioned vector of values which are interpolated from p.
    integer(kind=i4),intent(out) :: ier  !! an integer(kind=i4) error flag set as follows:
                                !!
                                !! * ier = 0 if no errors in input arguments are detected
                                !! * ier = 1 if  mx < 2
                                !! * ier = 2 if nx < 2 when intpol=1 or nx < 4 when intpol=3.
                                !! * ier = 3 if intpol is not equal to 1 or 3
                                !! * ier = 4 if lw or liw is too small (insufficient work space)

    integer(kind=i4)  :: inmx,isubx,i2,i3,i4,i5,lwmin

    ! check input arguments

    ! check mx
    ier = 1
    if (mx < 2) return

    ! check intpol
    ier = 3
    if (intpol/=1 .and. intpol/=3) return

    ! check nx
    ier = 2
    if (intpol==1 .and. nx<2) return
    if (intpol==3 .and. nx<4) return

    ! set subgrid integer(kind=i4) indicator
    inmx = (nx-1)/(mx-1)
    isubx = nx - inmx*(mx-1)

    ! set minimum and check work space
    ier = 4
    if (isubx/=1) then
        if (intpol==1) lwmin = mx
        if (intpol==3) lwmin = 4*mx
    else
        lwmin = 1
    end if
    if (lw < lwmin) return
    if (liw < mx) return

    ! input arguments o.k.

    ier = 0

    ! preset pointers

    i2 = 1
    i3 = 1
    i4 = 1
    i5 = 1
    if (intpol == 1) then
        ! linear interpolation in x
        if (isubx /= 1) then
            call linmxu(nx,mx,iw,w)
        end if
        call lint1u(nx,p,mx,q,iw,w,inmx,isubx)
    else
        ! cubic interpolation in x
        if (isubx /= 1) then
            i2 = 1
            i3 = i2+mx
            i4 = i3+mx
            i5 = i4+mx
            call cubnmxu(nx,mx,iw,w(i2),w(i3),w(i4),w(i5))
        end if
        call cubt1u(nx,p,mx,q,iw,w(i2),w(i3),w(i4),w(i5),inmx,isubx)
    end if

    end subroutine rgrd1u
!**************************************************************************

!**************************************************************************
!>
    subroutine lint1u(nx,p,mx,q,ix,dx,inmx,isubx)
       !dir$ attributes forceinline :: lint1u
       !dir$ attributes code_align : 32 :: lint1u
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: lint1u
    implicit none

    integer(kind=i4) :: nx,mx,ix(mx),inmx,isubx,i,ii
    real(kind=dp) :: p(nx),q(mx),dx(mx)
    real(kind=dp), automatic :: t0,t1
    if (isubx == 1) then
        ! mx grid is subset of nx grid so q can be set directly
    !dir$ assume_aligned q:64
    !dir$ assume_aligned p:64
#if (AUTO_VECTORIZE) == 1
    !dir$ ivdep
    !dir$ vector aligned
    !dir$ vector always
    !dir$ vector vectorlength(8)
#else
    !dir$ novector
#endif
        do ii=1,mx
            i = inmx*(ii-1)+1
            q(ii) = p(i)
        end do
    else
! linearly interpolate
    !dir$ assume_aligned ix:64
    !dir$ assume_aligned dx:64
    !dir$ assume_aligned q:64
    !dir$ assume_aligned p:64
#if (AUTO_VECTORIZE) == 1
    !dir$ ivdep
    !dir$ vector aligned
    !dir$ vector always
    !dir$ vector vectorlength(8)
#else
    !dir$ novector
#endif
        do ii=1,mx
            i  = ix(ii)
            t0 = p(i)
            t1 = p(i+1) 
            q(ii) = t0+dx(ii)*(t1-t0)
            !q(ii) = p(i)+dx(ii)*(p(i+1)-p(i))
        end do
    end if

    end subroutine lint1u
!**************************************************************************

!**************************************************************************
!>
    subroutine cubt1u(nx,p,mx,q,ix,dxm,dx,dxp,dxpp,inmx,isubx)
       !dir$ attributes forceinline :: cubt1u
       !dir$ attributes code_align : 32 :: cubt1u
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: cubt1u
    implicit none

    integer(kind=i4)  :: nx,mx,ix(mx),inmx,isubx,i,ii
    real(kind=dp) :: p(nx),q(mx)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)
    real(kind=dp), automatic :: t0,t1,t2,t3
    if (isubx == 1) then
        ! mx grid is subset of nx grid so q can be set directly
    !dir$ assume_aligned q:64
    !dir$ assume_aligned p:64
#if (AUTO_VECTORIZE) == 1
    !dir$ ivdep
    !dir$ vector aligned
    !dir$ vector always
    !dir$ vector vectorlength(8)
#else
    !dir$ novector
#endif
        do ii=1,mx
            i = inmx*(ii-1)+1
            q(ii) = p(i)
        end do
    else
        ! cubically interpolate on uniform grid
    !dir$ assume_aligned dxm:64
    !dir$ assume_aligned dxp:64
    !dir$ assume_aligned dxpp:64
    !dir$ assume_aligned ix:64
    !dir$ assume_aligned dx:64
    !dir$ assume_aligned q:64
    !dir$ assume_aligned p:64
#if (AUTO_VECTORIZE) == 1
    !dir$ ivdep
    !dir$ vector aligned
    !dir$ vector always
    !dir$ vector vectorlength(8)
#else
    !dir$ novector
#endif
        do ii=1,mx
            i  = ix(ii)
            t0 = p(i-1)
            t1 = p(i)
            t2 = p(i+1)
            t3 = p(i+2)
            q(ii) = (dxm(ii)*t0+dx(ii)*t1+dxp(ii)*t2+dxpp(ii)*t3)
            !q(ii)=(dxm(ii)*p(i-1)+dx(ii)*p(i)+dxp(ii)*p(i+1)+dxpp(ii)*p(i+2))
        end do
    end if

    end subroutine cubt1u
!**************************************************************************

!**************************************************************************
!>
!  set linear interpolation terms

    subroutine linmxu(nx,mx,ix,dx)
       !dir$ attributes forceinline :: linmxu
       !dir$ attributes code_align : 32 :: linmxu
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: linmxu
    implicit none

    integer(kind=i4) :: nx,mx,ix(mx),i,ii
    real(kind=dp) :: dx(mx),dnx,dmx,x,xx

    ! set "virtual" uniform increments
    dnx = 1.0_dp/(nx-1)
    dmx = 1.0_dp/(mx-1)

    ! set ix(ii) = i  s.t. i,i+1 can interpolate for ii
  
    !dir$ assume_aligned ix:64
    !dir$ assume_aligned dx:64
#if (AUTO_VECTORIZE) == 1
    !dir$ ivdep
    !dir$ vector aligned
    !dir$ vector always
    !dir$ vector vectorlength(8)
#else
    !dir$ novector
#endif
    do ii=1,mx
        xx = (ii-1)*dmx
        ix(ii) = min(int(xx/dnx)+1,nx-1)
        ! set scale term for linear
        i = ix(ii)
        x = (i-1)*dnx
        dx(ii) = (xx-x)/dnx
    end do

    end subroutine linmxu
!**************************************************************************

!**************************************************************************
!>
! set cubic interpolation terms
    !dimension length "mx" must by multiplicty of 64 bytes, i.e. 8 doubles!!
    subroutine cubnmxu(nx,mx,ix,dxm,dx,dxp,dxpp)
       !dir$ attributes forceinline :: cubnmxu
       !dir$ attributes code_align : 32 :: cubnmxu
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: cubnmxu
    implicit none

    integer(kind=i4) :: nx,mx,ix(mx),i,ii
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx),dnx,dmx,odnx3
    real(kind=dp) :: xx,xim,xi,xip,xipp

    ! set "virtual" uniform increments
    dnx = 1.0_dp/(nx-1)
    dmx = 1.0_dp/(mx-1)
    odnx3 = 1.0_dp/(6.0_dp*dnx*dnx*dnx)

    ! set i=ix(ii) in [2,nx-2] such that
    ! i-1,i,i+1,i+2 can be used to interpolate at ii
    !dir$ assume_aligned dxp:64
    !dir$ assume_aligned dxpp:64
    !dir$ assume_aligned dxm:64
    !dir$ assume_aligned ix:64
    !dir$ assume_aligned dx:64
    
#if (AUTO_VECTORIZE) == 1
    !dir$ assume (mod(mx,8).eq.0)
    !dir$ ivdep
    !dir$ vector aligned
    !dir$ vector always
    !dir$ vector vectorlength(8)
#else
    !dir$ novector
#endif
    do ii=1,mx
        xx = (ii-1)*dmx
        ix(ii) = min(max(int(xx/dnx)+1,2),nx-2)
        i = ix(ii)
        ! set scale terms for cubic
        xi = (i-1)*dnx
        xim = xi-dnx
        xip = xi+dnx
        xipp = xip+dnx
        dxm(ii) = -(xx-xi)*(xx-xip)*(xx-xipp)*odnx3
        dx(ii) = 3.0_dp*(xx-xim)*(xx-xip)*(xx-xipp)*odnx3
        dxp(ii) = -3.0_dp*(xx-xim)*(xx-xi)*(xx-xipp)*odnx3
        dxpp(ii) = (xx-xim)*(xx-xi)*(xx-xip)*odnx3
    end do

    end subroutine cubnmxu
!**************************************************************************

!**************************************************************************
!>
!  subroutine rgrd2 interpolates the values p(i,j) on the orthogonal
!  grid (x(i),y(j)) for i=1,...,nx and j=1,...,ny onto q(ii,jj) on the
!  orthogonal grid (xx(ii),yy(jj)) for ii=1,...,mx and jj=1,...,my.
!
!### method
!
!  linear or cubic interpolation is used (independently) in
!  each direction (see argument intpol).
!
!### requirements
!
!  each of the x,y grids must be strictly montonically increasing
!  and each of the xx,yy grids must be montonically increasing (see
!  ier = 4).  in addition the (X,Y) region
!
!    [xx(1),xx(mx)] X [yy(1),yy(my)]
!
!  must lie within the (X,Y) region
!
!    [x(1),x(nx)] X [y(1),y(ny)].
!
!  extrapolation is not allowed (see ier=3).  if these (X,Y)
!  regions are identical and the orthogonal grids are UNIFORM
!  in each direction then subroutine rgrd2u
!  should be used instead of rgrd2.

    subroutine rgrd2(nx,ny,x,y,p,mx,my,xx,yy,q,intpol,w,lw,iw,liw,ier)
          !dir$ attributes code_align : 32 :: rgrd2
          !dir$ optimize : 3
    implicit none

    integer(kind=i4),intent(in)  :: nx      !! the integer(kind=i4) dimension of the grid vector x and the first dimension
                                   !! of p.  nx > 1 if intpol(1) = 1 or nx > 3 if intpol(1) = 3 is required.
    integer(kind=i4),intent(in)  :: ny      !! the integer(kind=i4) dimension of the grid vector y and the second dimension
                                   !! of p.  ny > 1 if intpol(2) = 1 or ny > 3 if intpol(2) = 3 is required.
    integer(kind=i4),intent(in)  :: mx      !! the integer(kind=i4) dimension of the grid vector xx and the first dimension
                                   !! of q.  mx > 0 is required.
    integer(kind=i4),intent(in)  :: my      !! the integer(kind=i4) dimension of the grid vector yy and the second dimension
                                   !! of q.  my > 0 is required.
    integer(kind=i4),intent(in)  :: lw      !! the integer(kind=i4) length of the real(kind=dp) work space w.  let
                                   !!
                                   !! * lwx = mx                if intpol(1) = 1
                                   !! * lwx = 4*mx              if intpol(1) = 3
                                   !! * lwy = my+2*mx           if intpol(2) = 1
                                   !! * lwy = 4*(mx+my)         if intpol(2) = 3
                                   !!
                                   !! then lw must be greater than or equal to lwx+lwy
    integer(kind=i4),intent(in)  :: liw  !! the integer(kind=i4) length of the integer(kind=i4) work space iw.  liw must be at least mx+my
    integer(kind=i4),intent(out)  :: ier  !! an integer(kind=i4) error flag set as follows:
                                 !!
                                 !! * ier = 0 if no errors in input arguments are detected
                                 !! * ier = 1 if  min(mx,my) < 1
                                 !! * ier = 2 if nx < 2 when intpol(1)=1 or nx < 4 when intpol(1)=3 (or)
                                 !!   ny < 2 when intpol(2)=1 or ny < 4 when intpol(2)=3
                                 !! * ier = 3 if xx(1) < x(1) or x(nx) < xx(mx) (or)
                                 !!   yy(1) < y(1) or y(ny) < yy(my) (or)
                                 !!   to avoid this flag when end points are intended to be the
                                 !!   same but may differ slightly due to roundoff error, they
                                 !!   should be set exactly in the calling routine (e.g., if both
                                 !!   grids have the same y boundaries then yy(1)=y(1) and yy(my)=y(ny)
                                 !!   should be set before calling rgrd2)
                                 !! * ier = 4 if one of the grids x,y is not strictly monotonically
                                 !!   increasing or if one of the grids xx,yy is not
                                 !!   montonically increasing.  more precisely if:
                                 !!
                                 !!    * x(i+1) <= x(i) for some i such that 1 <= i < nx (or)
                                 !!    * y(j+1) <= y(j) for some j such that 1 <= j < ny (or)
                                 !!    * xx(ii+1) < xx(ii) for some ii such that 1 <= ii < mx (or)
                                 !!    * yy(jj+1) < yy(jj) for some jj such that 1 <= jj < my
                                 !! * ier = 5 if lw or liw is to small (insufficient work space)
                                 !! * ier = 6 if intpol(1) or intpol(2) is not equal to 1 or 3
    integer(kind=i4),intent(in)  :: intpol(2)    !! an integer(kind=i4) vector of dimension 2 which sets linear or cubic
                                        !! interpolation in the x,y directions as follows:
                                        !!
                                        !! * intpol(1) = 1 sets linear interpolation in the x direction
                                        !! * intpol(1) = 3 sets cubic interpolation in the x direction.
                                        !! * intpol(2) = 1 sets linear interpolation in the y direction
                                        !! * intpol(2) = 3 sets cubic interpolation in the y direction.
                                        !!
                                        !! values other than 1 or 3 in intpol are not allowed (ier = 5).
    integer(kind=i4),intent(inout)  :: iw(liw)   !! an integer(kind=i4) work space of length at least
                                        !! liw which must be provided in the
                                        !! routine calling rgrd2
    real(kind=dp),intent(in) :: x(nx)    !! a real(kind=dp) nx vector of strictly increasing values which defines the x
                                    !! portion of the orthogonal grid on which p is given
    real(kind=dp),intent(in) :: y(ny)    !! a real(kind=dp) ny vector of strictly increasing values which defines the y
                                    !! portion of the orthogonal grid on which p is given
    real(kind=dp),intent(in) :: p(nx,ny) !! a real(kind=dp) nx by ny array of values given on the orthogonal (x,y) grid
    real(kind=dp),intent(in) :: xx(mx)   !! a real(kind=dp) mx vector of increasing values which defines the x portion of the
                                    !! orthogonal grid on which q is defined.  xx(1) < x(1) or xx(mx) > x(nx)
                                    !! is not allowed (see ier = 3)
    real(kind=dp),intent(in) :: yy(my)   !! a real(kind=dp) my vector of increasing values which defines the y portion of the
                                    !! orthogonal grid on which q is defined.  yy(1) < y(1) or yy(my) > y(ny)
                                    !! is not allowed (see ier = 3)
    real(kind=dp),intent(out) :: q(mx,my)    !! a real(kind=dp) mx by my array of values on the (xx,yy) grid which are
                                        !! interpolated from p on the (x,y) grid
    real(kind=dp),intent(inout) :: w(lw)     !! a real(kind=dp) work space of length at least
                                        !! lw which must be provided in the
                                        !! routine calling rgrd2

    integer(kind=i4)  :: i,ii,j,jj,j2,j3,j4,j5,j6,j7,j8,j9,i2,i3,i4,i5
    integer(kind=i4)  :: jy,lwx,lwy

    ! check input arguments

    ! check (xx,yy) grid resolution
    ier = 1
    if (min(mx,my) < 1) return

    ! check intpol
    ier = 6
    if (intpol(1)/=1 .and. intpol(1)/=3) return
    if (intpol(2)/=1 .and. intpol(2)/=3) return

    ! check (x,y) grid resolution
    ier = 2
    if (intpol(1)==1 .and. nx<2) return
    if (intpol(1)==3 .and. nx<4) return
    if (intpol(2)==1 .and. ny<2) return
    if (intpol(2)==3 .and. ny<4) return

    ! check work space lengths
    ier = 5
    if (intpol(1)==1) then
        lwx = mx
    else
        lwx = 4*mx
    end if
    if (intpol(2)==1) then
        lwy = my+2*mx
    else
        lwy = 4*(mx+my)
    end if
    if (lw < lwx+lwy) return
    if (liw < mx+my) return

    ! check (xx,yy) grid contained in (x,y) grid
    ier = 3
    if (xx(1)<x(1) .or. xx(mx)>x(nx)) return
    if (yy(1)<y(1) .or. yy(my)>y(ny)) return

    ! check montonicity of grids
    ier = 4
    do i=2,nx
        if (x(i-1)>=x(i)) return
    end do
    do j=2,ny
        if (y(j-1)>=y(j)) return
    end do
    do ii=2,mx
        if (xx(ii-1)>xx(ii)) return
    end do
    do jj=2,my
        if (yy(jj-1)>yy(jj)) return
    end do

    ! arguments o.k.

    ier = 0

    ! set pointer in integer(kind=i4) work space

    jy = mx+1
    if (intpol(2) ==1) then

        ! linearly interpolate in y

        j2 = 1
        j3 = j2
        j4 = j3+my
        j5 = j4
        j6 = j5
        j7 = j6
        j8 = j7+mx
        j9 = j8+mx

        ! set y interpolation indices and scales and linearly interpolate

        call linmx(ny,y,my,yy,iw(jy),w(j3))
        i2 = j9

        ! set work space portion and indices which depend on x interpolation

        if (intpol(1) == 1) then
            i3 = i2
            i4 = i3
            i5 = i4
            call linmx(nx,x,mx,xx,iw,w(i3))
        else
            i3 = i2+mx
            i4 = i3+mx
            i5 = i4+mx
            call cubnmx(nx,x,mx,xx,iw,w(i2),w(i3),w(i4),w(i5))
        end if
        call lint2(nx,ny,p,mx,my,q,intpol,iw(jy),w(j3),w(j7),w(j8),iw,w(i2),w(i3),w(i4),w(i5))

    else

        ! cubically interpolate in y, set indice pointers

        j2 = 1
        j3 = j2+my
        j4 = j3+my
        j5 = j4+my
        j6 = j5+my
        j7 = j6+mx
        j8 = j7+mx
        j9 = j8+mx
        call cubnmx(ny,y,my,yy,iw(jy),w(j2),w(j3),w(j4),w(j5))
        i2 =  j9+mx

        ! set work space portion and indices which depend on x interpolation

        if (intpol(1) == 1) then
            i3 = i2
            i4 = i3
            i5 = i4
            call linmx(nx,x,mx,xx,iw,w(i3))
        else
            i3 = i2+mx
            i4 = i3+mx
            i5 = i4+mx
            call cubnmx(nx,x,mx,xx,iw,w(i2),w(i3),w(i4),w(i5))
        end if
        call cubt2(nx,ny,p,mx,my,q,intpol,iw(jy),w(j2),w(j3),&
                   w(j4),w(j5),w(j6),w(j7),w(j8),w(j9),iw,w(i2),w(i3),w(i4),w(i5))

    end if

    end subroutine rgrd2
!**************************************************************************

!**************************************************************************
!>
!  linearly interpolate in y

    subroutine lint2(nx,ny,p,mx,my,q,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
       !dir$ attributes code_align : 32 :: lint2
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: lint2
       
    implicit none

    integer(kind=i4)  :: nx,ny,mx,my,intpol(2),jy(my),ix(mx)
    integer(kind=i4)  :: jsave,j,jj,ii
    real(kind=dp) :: p(nx,ny),q(mx,my)
    real(kind=dp) :: pj(mx),pjp(mx),dy(my)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)

    if (intpol(1)==1) then

        ! linear in x

        jsave = -1

        do jj=1,my
            j = jy(jj)
            if (j==jsave) then
                ! j pointer has not moved since last pass (no updates or interpolation)
            else if (j==jsave+1) then
                ! update j and interpolate j+1
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ vector aligned
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pj(ii) = pjp(ii)
                end do
                call lint1(nx,p(1,j+1),mx,pjp,ix,dx)
            else
                ! interpolate j,j+1in pj,pjp on xx mesh
                call lint1(nx,p(1,j),mx,pj,ix,dx)
                call lint1(nx,p(1,j+1),mx,pjp,ix,dx)
            end if

            ! save j pointer for next pass

            jsave = j

            ! linearly interpolate q(ii,jj) from pjp,pj in y direction
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ assume_aligned dy:64
                !dir$ vector aligned
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
            do ii=1,mx
                q(ii,jj) = pj(ii)+dy(jj)*(pjp(ii)-pj(ii))
            end do
        end do

    else

        ! cubic in x

        jsave = -1
        do jj=1,my
            j = jy(jj)
            if (j==jsave) then
                ! j pointer has not moved since last pass (no updates or interpolation)
            else if (j==jsave+1) then
                ! update j and interpolate j+1
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pj(ii) = pjp(ii)
                end do
                call cubt1(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp)
            else
                ! interpolate j,j+1 in pj,pjp on xx mesh
                call cubt1(nx,p(1,j),mx,pj,ix,dxm,dx,dxp,dxpp)
                call cubt1(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp)
            end if

            ! save j pointer for next pass

            jsave = j

            ! linearly interpolate q(ii,jj) from pjp,pj in y direction
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ assume_aligned dy:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
            do ii=1,mx
                q(ii,jj) = pj(ii)+dy(jj)*(pjp(ii)-pj(ii))
            end do
        end do

    end if

    end subroutine lint2
!**************************************************************************

!**************************************************************************
!>
    subroutine cubt2(nx,ny,p,mx,my,q,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
       !dir$ attributes code_align : 32 :: cubt2
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: cubt2
    implicit none

    integer(kind=i4)  :: nx,ny,mx,my,intpol(2),jy(my),ix(mx)
    integer(kind=i4)  :: jsave,j,jj,ii
    real(kind=dp) :: p(nx,ny),q(mx,my)
    real(kind=dp) :: pjm(mx),pj(mx),pjp(mx),pjpp(mx)
    real(kind=dp) :: dym(my),dy(my),dyp(my),dypp(my)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)

    if (intpol(1)==1) then

        ! linear in x

        jsave = -3
        do jj=1,my

            ! load closest four j lines containing interpolate on xx mesh
            ! for j-1,j,j+1,j+2 in pjm,pj,pjp,pjpp

            j = jy(jj)
            if (j==jsave) then
                ! j pointer has not moved since last pass (no updates or interpolation)
            else if (j==jsave+1) then
                ! update j-1,j,j+1 and interpolate j+2
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pjm(ii) = pj(ii)
                    pj(ii) = pjp(ii)
                    pjp(ii) = pjpp(ii)
                end do
                call lint1(nx,p(1,j+2),mx,pjpp,ix,dx)
            else if (j==jsave+2) then
                ! update j-1,j and interpolate j+1,j+2
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pjm(ii) = pjp(ii)
                    pj(ii) = pjpp(ii)
                end do
                call lint1(nx,p(1,j+1),mx,pjp,ix,dx)
                call lint1(nx,p(1,j+2),mx,pjpp,ix,dx)
            else if (j==jsave+3) then
                ! update j-1 and interpolate j,j+1,j+2
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pjm(ii) = pjpp(ii)
                end do
                call lint1(nx,p(1,j),mx,pj,ix,dx)
                call lint1(nx,p(1,j+1),mx,pjp,ix,dx)
                call lint1(nx,p(1,j+2),mx,pjpp,ix,dx)
            else
                ! interpolate all four j-1,j,j+1,j+2
                call lint1(nx,p(1,j-1),mx,pjm,ix,dx)
                call lint1(nx,p(1,j),mx,pj,ix,dx)
                call lint1(nx,p(1,j+1),mx,pjp,ix,dx)
                call lint1(nx,p(1,j+2),mx,pjpp,ix,dx)
            end if

            ! save j pointer for next pass

            jsave = j

            ! cubically interpolate q(ii,jj) from pjm,pj,pjp,pjpp in y direction
                !dir$ assume_aligned q:64
                !dir$ assume_aligned dym:64
                !dir$ assume_aligned dy:64
                !dir$ assume_aligned dyp:64
                !dir$ assume_aligned dypp:64
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                !dir$ fma
            do ii=1,mx
                q(ii,jj) = dym(jj)*pjm(ii)+dy(jj)*pj(ii)+dyp(jj)*pjp(ii)+dypp(jj)*pjpp(ii)
            end do
        end do
        return

    else

        ! cubic in x

        jsave = -3
        do jj=1,my

            ! load closest four j lines containing interpolate on xx mesh
            ! for j-1,j,j+1,j+2 in pjm,pj,pjp,pjpp

            j = jy(jj)
            if (j==jsave) then
                ! j pointer has not moved since last pass (no updates or interpolation)
            else if (j==jsave+1) then
                ! update j-1,j,j+1 and interpolate j+2
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pjm(ii) = pj(ii)
                    pj(ii) = pjp(ii)
                    pjp(ii) = pjpp(ii)
                end do
                call cubt1(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp)
            else if (j==jsave+2) then
                ! update j-1,j and interpolate j+1,j+2
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pjm(ii) = pjp(ii)
                    pj(ii) = pjpp(ii)
                end do
                call cubt1(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp)
                call cubt1(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp)
            else if (j==jsave+3) then
                ! update j-1 and interpolate j,j+1,j+2
                
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pjm(ii) = pjpp(ii)
                end do
                call cubt1(nx,p(1,j),mx,pj,ix,dxm,dx,dxp,dxpp)
                call cubt1(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp)
                call cubt1(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp)
            else
                ! interpolate all four j-1,j,j+1,j+2
                call cubt1(nx,p(1,j-1),mx,pjm,ix,dxm,dx,dxp,dxpp)
                call cubt1(nx,p(1,j),mx,pj,ix,dxm,dx,dxp,dxpp)
                call cubt1(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp)
                call cubt1(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp)
            end if

            ! save j pointer for next pass

            jsave = j

            ! cubically interpolate q(ii,jj) from pjm,pj,pjp,pjpp in y direction
                !dir$ assume_aligned q:64
                !dir$ assume_aligned dym:64
                !dir$ assume_aligned dy:64
                !dir$ assume_aligned dyp:64
                !dir$ aassume_aligned dypp:64
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
            do ii=1,mx
                q(ii,jj) = dym(jj)*pjm(ii)+dy(jj)*pj(ii)+dyp(jj)*pjp(ii)+dypp(jj)*pjpp(ii)
            end do
        end do

    end if

    end subroutine cubt2
!**************************************************************************

!**************************************************************************
!>
!  subroutine rgrd2u interpolates the nx by ny array p onto
!  the mx by my array q.  linear or cubic interpolation is
!  used in each direction (see argument intpol).  it is assumed
!  that p and q are values on uniform nx by ny and mx by my grids
!  superimposed on the same rectangle (INCLUDING BOUNDARIES).
!  if p and q are values on nonuniform orthogonal grids and/or
!  if the grid on which q is defined lies within the p grid
!  then subroutine rgrd2 should be used.
!
!### method
!
!  linear or cubic interpolation (see intpol) is used in each
!  direction for which the q grid is not a subgrid of the p grid.
!  [the mx (my) uniform grid is a subgrid of the nx (ny) uniform
!  grid if and only if mx-1 (my-1) divides nx-1 (ny-1)].
!  values are set directly without (the need for) interpolation
!  in subgrid directions.

    subroutine rgrd2u(nx,ny,p,mx,my,q,intpol,w,lw,iw,liw,ier)
           !dir$ attributes code_align : 32 :: rgrd2u
           !dir$ optimize : 3
    implicit none

    integer(kind=i4),intent(in)  :: nx   !! the integer(kind=i4) first dimension of p.  nx > 1 if intpol(1) = 1 or
                                !! nx > 3 if intpol(1) = 3 is required (see ier = 2).
    integer(kind=i4),intent(in)  :: ny   !! the integer(kind=i4) second dimension of p.  ny > 1 if intpol(2) = 1 or
                                !! ny > 3 if intpol(2) = 3 is required (see ier = 2).
    integer(kind=i4),intent(in)  :: mx   !! the integer(kind=i4) first dimension of q.  mx > 1 is required (see ier = 1)
    integer(kind=i4),intent(in)  :: my   !! the integer(kind=i4) second dimension of q. my > 1 is required (see ier = 1)
    integer(kind=i4),intent(in)  :: intpol(2)    !! an integer(kind=i4) vector of dimension 2 which sets linear or cubic
                                        !! interpolation in each of the x,y directions as follows:
                                        !!
                                        !! * intpol(1) = 1 sets linear interpolation in the x direction
                                        !! * intpol(1) = 3 sets cubic interpolation in the x direction.
                                        !! * intpol(2) = 1 sets linear interpolation in the y direction
                                        !! * intpol(2) = 3 sets cubic interpolation in the y direction.
                                        !!
                                        !! values other than 1 or 3 in intpol are not allowed (ier = 3).
    integer(kind=i4),intent(in)  :: lw   !! the integer(kind=i4) length of the work space w.
                                !!
                                !! * let lwx = 1 if mx-1 divides nx-1; otherwise
                                !!   let lwx = mx if intpol(1) = 1 or
                                !!   let lwx = 4*mx if intpol(1) = 3
                                !! * let lwy = 0 if my-1 divides ny-1; otherwise
                                !!   let lwy = 2*mx+my if intpol(2) = 1 or
                                !!   let lwy = 4*(mx+my)  if intpol(2) = 3
                                !!
                                !! then lw must be greater than or equal to lwx+lwy
    integer(kind=i4),intent(in)  :: liw  !! the integer(kind=i4) length of the integer(kind=i4) work space iw.
                                !! liw must be greater than or equal to mx+my.
    integer(kind=i4),intent(out)  :: ier !! an integer(kind=i4) error flag set as follows:
                                !!
                                !! * ier = 0 if no errors in input arguments are detected
                                !! * ier = 1 if  min(mx,my) < 2
                                !! * ier = 2 if nx < 2 when intpol(1)=1 or nx < 4 when intpol(1)=3 (or)
                                !!   ny < 2 when intpol(2)=1 or ny < 4 when intpol(2)=3.
                                !! * ier = 3 if intpol(1) or intpol(2) is not equal to 1 or 3
                                !! * ier = 4 if lw or liw is to small (insufficient work space)
    real(kind=dp),intent(in) :: p(nx,ny)   !! a real(kind=dp) nx by ny array of given values
    real(kind=dp),intent(out) :: q(mx,my)  !! a real(kind=dp) mx by my array of values which are interpolated from p.
    real(kind=dp),intent(inout) :: w(lw)     !! a real(kind=dp) work space of length at
                                        !! least lw which must be provided in the
                                        !! routine calling rgrd2u
    integer(kind=i4),intent(inout)  :: iw(liw)   !! an integer(kind=i4) work space of length at least
                                        !! liw which must be provided in the
                                        !! routine calling rgrd2u

    integer(kind=i4)  :: inmx,jnmy,isubx,jsuby,lwx,lwy,jy
    integer(kind=i4)  :: j2,j3,j4,j5,j6,j7,j8,j9,i2,i3,i4,i5

    ! check input aarguments

    ! check mx,my
    ier = 1
    if (min(mx,my) < 2) return

    ! check intpol
    ier = 3
    if (intpol(1)/=1 .and. intpol(1)/=3) return
    if (intpol(2)/=1 .and. intpol(2)/=3) return

    ! check nx,ny
    ier = 2
    if (intpol(1)==1 .and. nx<2) return
    if (intpol(1)==3 .and. nx<4) return
    if (intpol(2)==1 .and. ny<2) return
    if (intpol(2)==3 .and. ny<4) return

    ! set subgrid indicators
    inmx = (nx-1)/(mx-1)
    jnmy = (ny-1)/(my-1)
    isubx = nx - inmx*(mx-1)
    jsuby = ny - jnmy*(my-1)

    ! check work space length input
    ier = 4
    lwx = 1
    lwy = 0
    if (isubx/=1) then
        if (intpol(1)==1) then
            lwx = mx
        else
            lwx = mx
        end if
    end if
    if (jsuby/=1) then
        if (intpol(2)==1) then
            lwy = (my+2*mx)
        else
            lwy = 4*(mx+my)
        end if
    end if
    if (lw < lwx+lwy) return
    if (liw < mx+my) return

    ! input arguments o.k.

    ier = 0
    jy = mx+1

    ! preset work space pointers

    j2 = 1
    j3 = j2
    j4 = j3
    j5 = j4
    j6 = j5
    j7 = j6
    j8 = j7
    j9 = j8
    i2 = j9
    i3 = i2
    i4 = i3
    i5 = i4

    if (intpol(2) ==1) then

        ! linearly interpolate in y

        if (jsuby/=1) then
            j2 = 1
            j3 = j2
            j4 = j3+my
            j5 = j4
            j6 = j5
            j7 = j6
            j8 = j7+mx
            j9 = j8+mx
            ! set y interpolation indices and scales and linearly interpolate
            call linmxu(ny,my,iw(jy),w(j3))
            i2 = j9
        end if

        ! set work space portion and indices which depend on x interpolation

        if (isubx/=1) then
            if (intpol(1) == 1) then
                i3 = i2
                i4 = i3
                i5 = i4
                call linmxu(nx,mx,iw,w(i3))
            else
                i3 = i2+mx
                i4 = i3+mx
                i5 = i4+mx
                call cubnmxu(nx,mx,iw,w(i2),w(i3),w(i4),w(i5))
            end if
        end if
        call lint2u(nx,ny,p,mx,my,q,intpol,iw(jy),w(j3),w(j7),w(j8),iw,&
                    w(i2),w(i3),w(i4),w(i5),inmx,jnmy,isubx,jsuby)
        return

    else

        ! cubically interpolate in y, set indice pointers

        if (jsuby/=1) then
            j2 = 1
            j3 = j2+my
            j4 = j3+my
            j5 = j4+my
            j6 = j5+my
            j7 = j6+mx
            j8 = j7+mx
            j9 = j8+mx
            ! set y interpolation indices and scales and cubically interpolate in y
            call cubnmxu(ny,my,iw(jy),w(j2),w(j3),w(j4),w(j5))
            i2 =  j9+mx
        end if

        ! set work space portion and indices which depend on x interpolation

        if (isubx/=1) then
            if (intpol(1) == 1) then
                i3 = i2
                i4 = i3
                i5 = i4
                call linmxu(nx,mx,iw,w(i3))
            else
                i3 = i2+mx
                i4 = i3+mx
                i5 = i4+mx
                call cubnmxu(nx,mx,iw,w(i2),w(i3),w(i4),w(i5))
            end if
        end if
        call cubt2u(nx,ny,p,mx,my,q,intpol,iw(jy),w(j2),w(j3),w(j4),&
                    w(j5),w(j6),w(j7),w(j8),w(j9),iw,w(i2),w(i3),w(i4),w(i5),&
                    inmx,jnmy,isubx,jsuby)
        return
    end if

    end subroutine
!**************************************************************************

!**************************************************************************
!>
!  linearly interpolate p onto q in y

    subroutine lint2u(nx,ny,p,mx,my,q,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
       !dir$ attributes forceinline :: lint2u
       !dir$ attributes code_align : 32 :: lint2u
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: lint2u
    implicit none

    integer(kind=i4) :: nx,ny,mx,my,intpol(2),jy(my),ix(mx),inmx,jnmy,isubx,jsuby
    integer(kind=i4) :: j,jj,ii,jsave
    real(kind=dp) :: p(nx,ny),q(mx,my)
    real(kind=dp) :: dy(my),pj(mx),pjp(mx)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)

    if (intpol(1) == 1) then

        ! linear in x

        if (jsuby == 1) then
            ! my grid is subset of ny grid
            do jj=1,my
                j = jnmy*(jj-1)+1
                call lint1u(nx,p(1,j),mx,q(1,jj),ix,dx,inmx,isubx)
            end do
            return
        end if

        jsave = -1
        do jj=1,my
            j = jy(jj)
            if (j == jsave) then
                ! pointer has not moved, no interpolation in pj,pjp necessary
            else if (j == jsave+1) then
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
#if (AUTO_VECTORIZE) == 1
                !!dir$ ivdep
                !dir$ vector aligned
                !dir$ vector always
                !dir$ vector vectorlength(8)
                !dir$ assume (mod(mx,8) .eq. 0)
#else
   
    !dir$ novector
#endif                
                do ii=1,mx
                    pj(ii) = pjp(ii)
                end do
                call lint1u(nx,p(1,j+1),mx,pjp,ix,dx,inmx,isubx)
            else
                call lint1u(nx,p(1,j),mx,pj,ix,dx,inmx,isubx)
                call lint1u(nx,p(1,j+1),mx,pjp,ix,dx,inmx,isubx)
            end if
            ! update pointer
            jsave = j
                !dir$ assume_aligned dy:64
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
#if (AUTO_VECTORIZE) == 1
                !!dir$ ivdep
                !dir$ vector aligned
                !dir$ vector always
                !dir$ vector vectorlength(8)
                !dir$ assume (mod(mx,8) .eq. 0)
#else
   
    !dir$ novector
#endif                         
            do ii=1,mx
                q(ii,jj) = pj(ii)+dy(jj)*(pjp(ii)-pj(ii))
            end do
        end do

    else

        ! cubic in x

        if (jsuby == 1) then
            ! my grid is subset of ny grid
            do jj=1,my
                j = jnmy*(jj-1)+1
                call cubt1u(nx,p(1,j),mx,q(1,jj),ix,dxm,dx,dxp,dxpp,inmx,isubx)
            end do
            return
        end if

        jsave = -1
        do jj=1,my
            j = jy(jj)
            if (j == jsave) then
                ! no interpolation in pj,pjp necessary
            else if (j == jsave+1) then
                do ii=1,mx
                    pj(ii) = pjp(ii)
                end do
                call cubt1u(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp,inmx,isubx)
            else
                call cubt1u(nx,p(1,j),mx,pj,ix,dxm,dx,dxp,dxpp,inmx,isubx)
                call cubt1u(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp,inmx,isubx)
            end if
            ! update pointer
            jsave = j
                !dir$ assume_aligned dy:64
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
#if (AUTO_VECTORIZE) == 1
                !!dir$ ivdep
                !dir$ vector aligned
                !dir$ vector always
                !dir$ vector vectorlength(8)
                !dir$ assume (mod(mx,8) .eq. 0)
#else
   
                !dir$ novector
#endif               
            do ii=1,mx
                q(ii,jj) = pj(ii)+dy(jj)*(pjp(ii)-pj(ii))
            end do
        end do

    end if

    end subroutine lint2u
!**************************************************************************

!**************************************************************************
!>
!  cubically interpolate p onto q in y

    subroutine cubt2u(nx,ny,p,mx,my,q,intpol,jy,dym,dy,dyp,dypp,&
                      pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,&
                      isubx,jsuby)
       !dir$ attributes forceinline :: cubt2u
       !dir$ attributes code_align : 32 :: cubt2u
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: cubt2u
    implicit none

    integer(kind=i4)  :: nx,ny,mx,my,intpol(2),jy(my),ix(mx),inmx,jnmy,isubx,jsuby
    integer(kind=i4)  :: j,jj,ii,jsave
    real(kind=dp) :: p(nx,ny),q(mx,my)
    real(kind=dp) :: dym(my),dy(my),dyp(my),dypp(my)
    real(kind=dp) :: pjm(mx),pj(mx),pjp(mx),pjpp(mx)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)

    if (intpol(1) == 1) then

        ! linear in x

        if (jsuby == 1) then
            ! my grid is subset of ny grid
            do jj=1,my
                j = jnmy*(jj-1)+1
                call lint1u(nx,p(1,j),mx,q(1,jj),ix,dx,inmx,isubx)
            end do
            return
        end if

        jsave = -3
        do jj=1,my
            j = jy(jj)

            ! load pjm,pj,pjp,pjpp

            if (j==jsave) then
                ! no updates or x interpolation necessary
            else if (j==jsave+1) then
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ assume_aligned pjp:64
                !dir$ vector aligned
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pjm(ii) = pj(ii)
                    pj(ii) = pjp(ii)
                    pjp(ii) = pjpp(ii)
                end do
                call lint1u(nx,p(1,j+2),mx,pjpp,ix,dx,inmx,isubx)
            else if (j==jsave+2) then
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ assume_aligned pjp:64
                !dir$ vector aligned
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pjm(ii) = pjp(ii)
                    pj(ii) = pjpp(ii)
                end do
                call lint1u(nx,p(1,j+1),mx,pjp,ix,dx,inmx,isubx)
                call lint1u(nx,p(1,j+2),mx,pjpp,ix,dx,inmx,isubx)
            else if (j==jsave+3) then
                do ii=1,mx
                    pjm(ii) = pjpp(ii)
                end do
                call lint1u(nx,p(1,j),mx,pj,ix,dx,inmx,isubx)
                call lint1u(nx,p(1,j+1),mx,pjp,ix,dx,inmx,isubx)
                call lint1u(nx,p(1,j+2),mx,pjpp,ix,dx,inmx,isubx)
            else
                ! load all four (no updates)
                call lint1u(nx,p(1,j-1),mx,pjm,ix,dx,inmx,isubx)
                call lint1u(nx,p(1,j),mx,pj,ix,dx,inmx,isubx)
                call lint1u(nx,p(1,j+1),mx,pjp,ix,dx,inmx,isubx)
                call lint1u(nx,p(1,j+2),mx,pjpp,ix,dx,inmx,isubx)
            end if
            ! update pointer
            jsave = j
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ assume_aligned dy:64
                !dir$ assume_aligned dym:64
                !dir$ assume_aligned dyp:64
                !dir$ assume_aligned dypp:64
                !dir$ assume_aligned pjpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
            do ii=1,mx
                q(ii,jj) = dym(jj)*pjm(ii) + dy(jj)*pj(ii) + dyp(jj)*pjp(ii) + dypp(jj)*pjpp(ii)
            end do
        end do

    else

        ! cubic in x

        if (jsuby == 1) then
            ! my grid is subset of ny grid
            do jj=1,my
                j = jnmy*(jj-1)+1
                call cubt1u(nx,p(1,j),mx,q(1,jj),ix,dxm,dx,dxp,dxpp,inmx,isubx)
            end do
            return
        end if

        jsave = -3
        do jj=1,my
            j = jy(jj)

            ! load pjm,pj,pjp,pjpp

            if (j==jsave) then
                ! no updates or x interpolation necessary
            else if (j==jsave+1) then
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ assume_aligned pjp:64
                !dir$ vector aligned
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pjm(ii) = pj(ii)
                    pj(ii) = pjp(ii)
                    pjp(ii) = pjpp(ii)
                end do
                call cubt1u(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp,inmx,isubx)
            else if (j==jsave+2) then
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ assume_aligned pjp:64
                !dir$ vector aligned
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pjm(ii) = pjp(ii)
                    pj(ii) = pjpp(ii)
                end do
                call cubt1u(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp,inmx,isubx)
                call cubt1u(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp,inmx,isubx)
            else if (j==jsave+3) then
                 
                !dir$ assume_aligned pjm:64
                !dir$ assume_aligned pjpp:64
                !dir$ vector aligned
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do ii=1,mx
                    pjm(ii) = pjpp(ii)
                end do
                call cubt1u(nx,p(1,j),mx,pj,ix,dxm,dx,dxp,dxpp,inmx,isubx)
                call cubt1u(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp,inmx,isubx)
                call cubt1u(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp,inmx,isubx)
            else
                ! load all four (no updates)
                call cubt1u(nx,p(1,j-1),mx,pjm,ix,dxm,dx,dxp,dxpp,inmx,isubx)
                call cubt1u(nx,p(1,j),mx,pj,ix,dxm,dx,dxp,dxpp,inmx,isubx)
                call cubt1u(nx,p(1,j+1),mx,pjp,ix,dxm,dx,dxp,dxpp,inmx,isubx)
                call cubt1u(nx,p(1,j+2),mx,pjpp,ix,dxm,dx,dxp,dxpp,inmx,isubx)
            end if
            ! update pointer
            jsave = j
                !dir$ assume_aligned pj:64
                !dir$ assume_aligned pjp:64
                !dir$ assume_aligned dy:64
                !dir$ assume_aligned dym:64
                !dir$ assume_aligned dyp:64
                !dir$ assume_aligned dypp:64
                !dir$ assume_aligned pjpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
            do ii=1,mx
                q(ii,jj) = dym(jj)*pjm(ii) + dy(jj)*pj(ii) + dyp(jj)*pjp(ii) + dypp(jj)*pjpp(ii)
            end do
        end do
        return

    end if

    end subroutine cubt2u
!**************************************************************************

!**************************************************************************
!>
!  subroutine rgrd3 interpolates the values p(i,j,k) on the orthogonal
!  grid (x(i),y(j),z(k)) for i=1,...,nx; j=1,...,ny; k=1,...,nz
!  onto q(ii,jj,kk) on the orthogonal grid (xx(ii),yy(jj),zz(kk)) for
!  ii=1,...,mx; jj=1,...,my; kk=1,...,mz.
!
!### method
!
!  linear or cubic interpolation is used (independently) in
!  each direction (see argument intpol).
!
!### requirements
!
!  each of the x,y,z grids must be strictly montonically increasing
!  and each of the xx,yy,zz grids must be montonically increasing
!  (see ier = 4).  in addition the (X,Y,Z) region
!
!    [xx(1),xx(mx)] X [yy(1),yy(my)] X [zz(1),zz(mz)]
!
!  must lie within the (X,Y,Z) region
!
!    [x(1),x(nx)] X [y(1),y(ny)] X [z(1),z(nz)].
!
!  extrapolation is not allowed (see ier=3).  if these (X,Y,Z)
!  regions are identical and the orthogonal grids are UNIFORM
!  in each direction then subroutine rgrd3u
!  should be used instead of rgrd3.

    subroutine rgrd3(nx,ny,nz,x,y,z,p,mx,my,mz,xx,yy,zz,q,intpol,w,lw,iw,liw,ier)
       !dir$ attributes code_align : 32 :: rgrd3
       !dir$ optimize : 3
    implicit none

    integer(kind=i4),intent(in)      :: nx   !! the integer(kind=i4) dimension of the grid vector x and the first dimension of p.
                                    !! nx > 1 if intpol(1) = 1 or nx > 3 if intpol(1) = 3 is required.
    integer(kind=i4),intent(in)      :: ny   !! the integer(kind=i4) dimension of the grid vector y and the second dimension of p.
                                    !! ny > 1 if intpol(2) = 1 or ny > 3 if intpol(2) = 3 is required.
    integer(kind=i4),intent(in)      :: nz   !! the integer(kind=i4) dimension of the grid vector z and the third dimension of p.
                                    !! nz > 1 if intpol(3) = 1 or nz > 3 if intpol(3) = 3 is required.
    integer(kind=i4),intent(in)      :: mx   !! the integer(kind=i4) dimension of the grid vector xx and the first dimension of q.
                                    !! mx > 0 is required.
    integer(kind=i4),intent(in)      :: my   !! the integer(kind=i4) dimension of the grid vector yy and the second dimension of q.
                                    !! my > 0 is required.
    integer(kind=i4),intent(in)      :: mz   !! the integer(kind=i4) dimension of the grid vector zz and the third dimension of q.
                                    !! mz > 0 is required.
    integer(kind=i4),intent(in)      :: lw   !! the integer(kind=i4) length of the real(kind=dp) work space w.  let
                                    !!
                                    !! * lwx = mx              if intpol(1) = 1
                                    !! * lwx = 4*mx            if intpol(1) = 3
                                    !! * lwy = my+2*mx         if intpol(2) = 1
                                    !! * lwy = 4*(mx+my)       if intpol(2) = 3
                                    !! * lwz = 2*mx*my+mz      if intpol(3) = 1
                                    !! * lwz = 4*(mx*my+mz)    if intpol(3) = 3
                                    !!
                                    !! then lw must be greater than or equal to lwx+lwy+lwz
    integer(kind=i4),intent(in)      :: liw  !! the integer(kind=i4) length of the integer(kind=i4) work space iw.  liw must be at least mx+my+mz
    integer(kind=i4),intent(out)     :: ier  !! an integer(kind=i4) error flag set as follows:
                                    !!
                                    !! * ier = 0 if no errors in input arguments are detected
                                    !! * ier = 1 if  min(mx,my,mz) < 1
                                    !! * ier = 2 if nx < 2 when intpol(1)=1 or nx < 4 when intpol(1)=3 (or)
                                    !!                     ny < 2 when intpol(2)=1 or ny < 4 when intpol(2)=3 (or)
                                    !!                     nz < 2 when intpol(3)=1 or nz < 4 when intpol(3)=3.
                                    !! * ier = 3 if xx(1) < x(1) or x(nx) < xx(mx) (or)
                                    !!                     yy(1) < y(1) or y(ny) < yy(my) (or)
                                    !!                     zz(1) < z(1) or z(nz) < zz(mz)
                                    !!    to avoid this flag when end points are intended to be the
                                    !!    same but may differ slightly due to roundoff error, they
                                    !!    should be set exactly in the calling routine (e.g., if both
                                    !!    grids have the same y boundaries then yy(1)=y(1) and yy(my)=y(ny)
                                    !!    should be set before calling rgrd3)
                                    !! * ier = 4 if one of the grids x,y,z is not strictly monotonically
                                    !!   increasing or if one of the grids xx,yy,zz is not
                                    !!   montonically increasing.  more precisely if:
                                    !!
                                    !!    * x(i+1) <= x(i) for some i such that 1 <= i < nx (or)
                                    !!    * y(j+1) <= y(j) for some j such that 1 <= j < ny (or)
                                    !!    * z(k+1) <= z(k) for some k such that 1 <= k < nz (or)
                                    !!    * xx(ii+1) < xx(ii) for some ii such that 1 <= ii < mx (or)
                                    !!    * yy(jj+1) < yy(jj) for some jj such that 1 <= jj < my (or)
                                    !!    * zz(kk+1) < zz(k)  for some kk such that 1 <= kk < mz
                                    !! * ier = 5 if lw or liw is too small (insufficient work space)
                                    !! * ier = 6 if any of intpol(1),intpol(2),intpol(3) is not equal to 1 or 3
    real(kind=dp),intent(in)     :: x(nx)  !! a real(kind=dp) nx vector of strictly increasing values which defines the x
                                      !! portion of the orthogonal grid on which p is given
    real(kind=dp),intent(in)     :: y(ny)  !! a real(kind=dp) ny vector of strictly increasing values which defines the y
                                      !! portion of the orthogonal grid on which p is given
    real(kind=dp),intent(in)     :: z(nz)  !! a real(kind=dp) nz vector of strictly increasing values which defines the z
                                      !! portion of the orthogonal grid on which p is given
    real(kind=dp),intent(in)     :: p(nx,ny,nz)  !! a real(kind=dp) nx by ny by nz array of values given on the (x,y,z) grid
    real(kind=dp),intent(in)     :: xx(mx)   !! a real(kind=dp) mx vector of increasing values which defines the x portion of the
                                        !! orthogonal grid on which q is defined.  xx(1) < x(1) or xx(mx) > x(nx)
                                        !! is not allowed (see ier = 3)
    real(kind=dp),intent(in)     :: yy(my)   !! a real(kind=dp) my vector of increasing values which defines the y portion of the
                                        !! orthogonal grid on which q is defined.  yy(1) < y(1) or yy(my) > y(ny)
                                        !! is not allowed (see ier = 3)
    real(kind=dp),intent(in)     :: zz(mz)   !! a real(kind=dp) mz vector of increasing values which defines the z portion of the
                                        !! orthogonal grid on which q is defined.  zz(1) < z(1) or zz(mz) > z(nz)
                                        !! is not allowed (see ier = 3)
    real(kind=dp),intent(out)    :: q(mx,my,mz)  !! a real(kind=dp) mx by my by mz array of values
                                            !! on the (xx,yy,zz) grid which are
                                            !! interpolated from p on the (x,y,z) grid
    real(kind=dp),intent(inout)  :: w(lw)    !! a real(kind=dp) work space of length at least
                                        !! lw which must be provided in the
                                        !! routine calling rgrd3
    integer(kind=i4),intent(in)      :: intpol(3)    !! an integer(kind=i4) vector of dimension 3 which sets linear or cubic
                                            !! interpolation in each of the x,y,z directions as follows:
                                            !!
                                            !! * intpol(1) = 1 sets linear interpolation in the x direction
                                            !! * intpol(1) = 3 sets cubic interpolation in the x direction.
                                            !! * intpol(2) = 1 sets linear interpolation in the y direction
                                            !! * intpol(2) = 3 sets cubic interpolation in the y direction.
                                            !! * intpol(3) = 1 sets linear interpolation in the z direction
                                            !! * intpol(3) = 3 sets cubic interpolation in the z direction.
                                            !!
                                            !! values other than 1 or 3 in intpol are not allowed (ier = 5).
    integer(kind=i4),intent(inout)   :: iw(liw)  !! an integer(kind=i4) work space of length at least
                                        !! liw which must be provided in the
                                        !! routine calling rgrd3

    integer(kind=i4)  :: i,ii,j,jj,k,kk
    integer(kind=i4)  :: i2,i3,i4,i5
    integer(kind=i4)  :: j2,j3,j4,j5,j6,j7,j8,j9
    integer(kind=i4)  :: k2,k3,k4,k5,k6,k7,k8,k9
    integer(kind=i4)  :: lwx,lwy,lwz,jy,kz,mxmy

    ! check input arguments

    ! check (xx,yy,zz) grid resolution
    ier = 1
    if (min(mx,my,mz) < 1) return

    ! check intpol
    ier = 6
    if (intpol(1)/=1 .and. intpol(1)/=3) return
    if (intpol(2)/=1 .and. intpol(2)/=3) return
    if (intpol(3)/=1 .and. intpol(3)/=3) return

    ! check (x,y,z) grid resolution
    ier = 2
    if (intpol(1)==1 .and. nx<2) return
    if (intpol(1)==3 .and. nx<4) return
    if (intpol(2)==1 .and. ny<2) return
    if (intpol(2)==3 .and. ny<4) return
    if (intpol(3)==1 .and. nz<2) return
    if (intpol(3)==3 .and. nz<4) return

    ! check work space length input and set minimum
    ier = 5
    mxmy = mx*my
    if (intpol(1)==1) then
        lwx = mx
    else
        lwx = 4*mx
    end if
    if (intpol(2)==1) then
        lwy = (my+2*mx)
    else
        lwy = 4*(my+mx)
    end if
    if (intpol(3)==1) then
        lwz = (2*mxmy+mz)
    else
        lwz = 4*(mxmy+mz)
    end if
    if (lw < lwx+lwy+lwz) return
    if (liw < mx+my+mz) return

    ! check (xx,yy,zz) grid contained in (x,y,z) grid
    ier = 3
    if (xx(1)<x(1) .or. xx(mx)>x(nx)) return
    if (yy(1)<y(1) .or. yy(my)>y(ny)) return
    if (zz(1)<z(1) .or. zz(mz)>z(nz)) return

    ! check montonicity of grids
    ier = 4
    do i=2,nx
        if (x(i-1)>=x(i)) return
    end do
    do j=2,ny
        if (y(j-1)>=y(j)) return
    end do
    do k=2,nz
        if (z(k-1)>=z(k)) return
    end do
    do ii=2,mx
        if (xx(ii-1)>xx(ii)) return
    end do
    do jj=2,my
        if (yy(jj-1)>yy(jj)) return
    end do
    do kk=2,mz
        if (zz(kk-1)>zz(kk)) return
    end do

    ! arguments o.k.

    ier = 0
    jy = mx+1
    kz = mx+my+1
    if (intpol(3)==1) then

        ! linearly interpolate in nz, set work space pointers and scales

        k2 = 1
        k3 = k2
        k4 = k3+mz
        k5 = k4
        k6 = k5
        k7 = k6
        k8 = k7+mxmy
        k9 = k8+mxmy
        call linmx(nz,z,mz,zz,iw(kz),w(k3))
        j2 = k9

        ! set indices and scales which depend on y interpolation

        if (intpol(2) == 1) then
            ! linear in y
            j3 = j2
            j4 = j3+my
            j5 = j4
            j6 = j5
            j7 = j6
            j8 = j7+mx
            j9 = j8+mx
            call linmx(ny,y,my,yy,iw(jy),w(j3))
            i2 = j9
        else
            ! cubic in y
            j3 = j2+my
            j4 = j3+my
            j5 = j4+my
            j6 = j5+my
            j7 = j6+mx
            j8 = j7+mx
            j9 = j8+mx
            call cubnmx(ny,y,my,yy,iw(jy),w(j2),w(j3),w(j4),w(j5))
            i2 = j9+mx
        end if

        ! set indices and scales which depend on x interpolation

        if (intpol(1) == 1) then
            ! linear in x
            i3 = i2
            i4 = i3
            i5 = i4
            call linmx(nx,x,mx,xx,iw,w(i3))
        else
            ! cubic in x
            i3 = i2+mx
            i4 = i3+mx
            i5 = i4+mx
            call cubnmx(nx,x,mx,xx,iw,w(i2),w(i3),w(i4),w(i5))
        end if
        call lint3(nx,ny,nz,p,mx,my,mxmy,mz,q,intpol,iw(kz),&
                    w(k3),w(k7),w(k8),iw(jy),w(j2),w(j3),w(j4),w(j5),w(j6),&
                    w(j7),w(j8),w(j9),iw,w(i2),w(i3),w(i4),w(i5))

    else

        ! cubically interpolate in z

        k2 = 1
        k3 = k2+mz
        k4 = k3+mz
        k5 = k4+mz
        k6 = k5+mz
        k7 = k6+mxmy
        k8 = k7+mxmy
        k9 = k8+mxmy
        call cubnmx(nz,z,mz,zz,iw(kz),w(k2),w(k3),w(k4),w(k5))
        j2 = k9+mxmy

        ! set indices which depend on y interpolation

        if (intpol(2) == 1) then
            j3 = j2
            j4 = j3+my
            j5 = j4
            j6 = j5
            j7 = j6
            j8 = j7+mx
            j9 = j8+mx
            call linmx(ny,y,my,yy,iw(jy),w(j3))
            i2 = j9
        else
            j3 = j2+my
            j4 = j3+my
            j5 = j4+my
            j6 = j5+my
            j7 = j6+mx
            j8 = j7+mx
            j9 = j8+mx
            call cubnmx(ny,y,my,yy,iw(jy),w(j2),w(j3),w(j4),w(j5))
            i2 = j9+mx
        end if

        ! set work space portion and indices which depend on x interpolation

        if (intpol(1) == 1) then
            i3 = i2
            i4 = i3
            i5 = i4
            call linmx(nx,x,mx,xx,iw,w(i3))
        else
            i3 = i2+mx
            i4 = i3+mx
            i5 = i4+mx
            call cubnmx(nx,x,mx,xx,iw,w(i2),w(i3),w(i4),w(i5))
        end if
        call cubt3(nx,ny,nz,p,mx,my,mxmy,mz,q,intpol,&
                    iw(kz),w(k2),w(k3),w(k4),w(k5),w(k6),w(k7),w(k8),w(k9),&
                    iw(jy),w(j2),w(j3),w(j4),w(j5),w(j6),w(j7),w(j8),w(j9),&
                    iw,w(i2),w(i3),w(i4),w(i5))

    end if

    end subroutine rgrd3
!**************************************************************************

!**************************************************************************
!>
!  linearly interpolate in z direction

    subroutine lint3(nx,ny,nz,p,mx,my,mxmy,mz,q,intpol,kz,&
                     dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,&
                     pjpp,ix,dxm,dx,dxp,dxpp)
       !dir$ attributes code_align : 32 :: lint3
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: lint3
    implicit none

    integer(kind=i4)  :: nx,ny,nz,mx,my,mz,mxmy
    real(kind=dp) :: p(nx,ny,nz),q(mxmy,mz)
    real(kind=dp) :: dz(mz),pk(mxmy),pkp(mxmy)
    real(kind=dp) :: dym(my),dy(my),dyp(my),dypp(my)
    real(kind=dp) :: pjm(mx),pj(mx),pjp(mx),pjpp(mx)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)
    integer(kind=i4)  :: intpol(3),kz(mz),jy(my),ix(mx)
    integer(kind=i4)  :: k,kk,iijj,ksave

    if (intpol(2) == 1) then

        ! linear in y

        ksave = -1
        do kk=1,mz
            k = kz(kk)
            if (k==ksave) then
                ! k pointer has not moved since last pass (no updates or interpolation)
            else if (k==ksave+1) then
                ! update k and interpolate k+1
                 !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pk(iijj) = pkp(iijj)
                end do
                call lint2(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
            else
                ! interpolate k,k+1 in pk,pkp on xx,yy mesh
                call lint2(nx,ny,p(1,1,k),mx,my,pk,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
                call lint2(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
            end if

            ! save k pointer for next pass

            ksave = k

            ! linearly interpolate q(ii,jj,k) from pk,pkp in z direction
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned dz:64
                !dir$ assume_aligned q:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
            do iijj=1,mxmy
                q(iijj,kk) = pk(iijj)+dz(kk)*(pkp(iijj)-pk(iijj))
            end do
        end do

    else

        ! cubic in y

        ksave = -1
        do kk=1,mz
            k = kz(kk)
            if (k==ksave) then
                ! k pointer has not moved since last pass (no updates or interpolation)
            else if (k==ksave+1) then
                ! update k and interpolate k+1
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pk(iijj) = pkp(iijj)
                end do
                call cubt2(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            else
                ! interpolate k,k+1 in pk,pkp on xx,yy mesh
                call cubt2(nx,ny,p(1,1,k),mx,my,pk,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call cubt2(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            end if

            ! save k pointer for next pass

            ksave = k

            ! linearly interpolate q(ii,jj,k) from pk,pkp in z direction
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned dz:64
                !dir$ assume_aligned q:64
                !dir$ vector aligned
                !dir$ assume (mod(mx,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
            do iijj=1,mxmy
                q(iijj,kk) = pk(iijj)+dz(kk)*(pkp(iijj)-pk(iijj))
            end do
        end do

    end if

    end subroutine lint3
!**************************************************************************

!**************************************************************************
!>
!  cubically interpolate in z

    subroutine cubt3(nx,ny,nz,p,mx,my,mxmy,mz,q,intpol,&
                     kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,&
                     jy,dym,dy,dyp,dypp,pjm,pj,&
                     pjp,pjpp,ix,dxm,dx,dxp,dxpp)
       !dir$ attributes code_align : 32 :: cubt3
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: cubt3
    implicit none

    integer(kind=i4)  :: nx,ny,nz,mx,my,mxmy,mz,k,kk,ksave,iijj
    real(kind=dp) :: p(nx,ny,nz),q(mxmy,mz)
    real(kind=dp) :: pkm(mxmy),pk(mxmy),pkp(mxmy),pkpp(mxmy)
    real(kind=dp) :: pjm(mx),pj(mx),pjp(mx),pjpp(mx)
    real(kind=dp) :: dzm(mz),dz(mz),dzp(mz),dzpp(mz)
    real(kind=dp) :: dym(my),dy(my),dyp(my),dypp(my)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)
    integer(kind=i4)  :: intpol(3),kz(mz),jy(my),ix(mx)

    if (intpol(2) == 1) then

        ! linear in y

        ksave = -3
        do kk=1,mz
            k = kz(kk)
            if (k==ksave) then
                ! k pointer has not moved since last pass (no updates or interpolation)
            else if (k==ksave+1) then
                ! update k-1,k,k+1 and interpolate k+2
                !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pk(iijj)
                    pk(iijj) = pkp(iijj)
                    pkp(iijj) = pkpp(iijj)
                end do
                call lint2(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
            else if (k==ksave+2) then
                ! update k-1,k and interpolate k+1,k+2
                 !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pkp(iijj)
                    pk(iijj) = pkpp(iijj)
                end do
                call lint2(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
                call lint2(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
            else if (k==ksave+3) then
                ! update k-1 and interpolate k,k+1,k+2
                 !dir$ assume_aligned pkm:64
                 !dir$ assume_aligned pkpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pkpp(iijj)
                end do
                call lint2(nx,ny,p(1,1,k),mx,my,pk,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
                call lint2(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
                call lint2(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
            else
                ! interpolate all four k-1,k,k+1,k+2
                call lint2(nx,ny,p(1,1,k-1),mx,my,pkm,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
                call lint2(nx,ny,p(1,1,k),mx,my,pk,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
                call lint2(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
                call lint2(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp)
            end if

            ! save k pointer for next pass

            ksave = k

            ! cubically interpolate q(ii,jj,kk) from pkm,pk,pkp,pkpp in z direction
                !dir$ assume_aligned q:64
                !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ assume_aligned dzm:64
                !dir$ assume_aligned dz:64
                !dir$ assume_aligned dzp:64
                !dir$ assume_aligned dzpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ assume (mod(mz,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
            do iijj=1,mxmy
                q(iijj,kk) = dzm(kk)*pkm(iijj) + dz(kk)*pk(iijj) + dzp(kk)*pkp(iijj) + dzpp(kk)*pkpp(iijj)
            end do
        end do

    else

        ! cubic in y

        ksave = -3
        do kk=1,mz
            k = kz(kk)
            if (k==ksave) then
                ! k pointer has not moved since last pass (no updates or interpolation)
            else if (k==ksave+1) then
                ! update k-1,k,k+1 and interpolate k+2
                !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pk(iijj)
                    pk(iijj) = pkp(iijj)
                    pkp(iijj) = pkpp(iijj)
                end do
                call cubt2(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            else if (k==ksave+2) then
                ! update k-1,k and interpolate k+1,k+2
                !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pkp(iijj)
                    pk(iijj) = pkpp(iijj)
                end do
                call cubt2(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call cubt2(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            else if (k==ksave+3) then
                ! update k-1 and interpolate k,k+1,k+2
                 !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pkpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pkpp(iijj)
                end do
                call cubt2(nx,ny,p(1,1,k),mx,my,pk,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call cubt2(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call cubt2(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            else
                ! interpolate all four k-1,k,k+1,k+2
                call cubt2(nx,ny,p(1,1,k-1),mx,my,pkm,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call cubt2(nx,ny,p(1,1,k),mx,my,pk,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call cubt2(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call cubt2(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            end if

            ! save k pointer for next pass

            ksave = k

            ! cubically interpolate q(ii,jj,kk) from pkm,pk,pkp,pkpp in z direction

               !dir$ assume_aligned q:64
                !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ assume_aligned dzm:64
                !dir$ assume_aligned dz:64
                !dir$ assume_aligned dzp:64
                !dir$ assume_aligned dzpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ assume (mod(mz,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8) 
            do iijj=1,mxmy
                q(iijj,kk) = dzm(kk)*pkm(iijj) + dz(kk)*pk(iijj) + dzp(kk)*pkp(iijj) + dzpp(kk)*pkpp(iijj)
            end do
        end do

    end if

    end subroutine cubt3
!**************************************************************************

!**************************************************************************
!>
!  subroutine rgrd3u interpolates the nx by ny by nz array p onto
!  the mx by my by mz array q.  it is assumed that p and q are
!  values on uniform nx by ny by nz and mx by my by mz grids which
!  are superimposed on the same box region (INCLUDING BOUNDARIES).
!  if p and q are values on nonuniform orthogonal grids and/or
!  if the grid on which q is defined lies within the p grid then
!  subroutine rgrd3 should be used.
!
!### method
!
!  linear or cubic interpolation (see intpol) is used in each
!  direction for which the q grid is not a subgrid of the p grid.
!  [the mx (my,mz) uniform grid is a subgrid of the nx (ny,nz) uniform
!  grid if and only if mx-1 (my-1,nz-1) divides nx-1 (ny-1,nz-1)].
!  Values are set directly without (the need for) interpolation
!  in subgrid directions.

    subroutine rgrd3u(nx,ny,nz,p,mx,my,mz,q,intpol,w,lw,iw,liw,ier)
       !dir$ attributes code_align : 32 :: rgrd3u
       !dir$ optimize : 3
    implicit none

    integer(kind=i4),intent(in)     :: nx    !! the integer(kind=i4) first dimension of p.  nx > 1 if intpol(1) = 1 or
                                    !! nx > 3 if intpol(1) = 3 is required (see ier = 2).
    integer(kind=i4),intent(in)     :: ny    !! the integer(kind=i4) second dimension of p.  ny > 1 if intpol(2) = 1 or
                                    !! ny > 3 if intpol(2) = 3 is required (see ier = 2).
    integer(kind=i4),intent(in)     :: nz    !! the integer(kind=i4) third dimension of p.  nz > 1 if intpol(3) = 1 or
                                    !! nz > 3 if intpol(3) = 3 is required (see ier = 2)
    integer(kind=i4),intent(in)     :: mx    !! the integer(kind=i4) first dimension of q.  mx > 1 is required (see ier = 1)
    integer(kind=i4),intent(in)     :: my    !! the integer(kind=i4) second dimension of q. my > 1 is required (see ier = 1)
    integer(kind=i4),intent(in)     :: mz    !! the integer(kind=i4) third dimension of q. mz > 1 is required (see ier = 1)
    integer(kind=i4),intent(in)     :: intpol(3) !! an integer(kind=i4) vector of dimension 3 which sets linear or cubic
                                        !! interpolation in each of the x,y,z directions as follows:
                                        !!
                                        !! * intpol(1) = 1 sets linear interpolation in the x direction
                                        !! * intpol(1) = 3 sets cubic interpolation in the x direction.
                                        !! * intpol(2) = 1 sets linear interpolation in the y direction
                                        !! * intpol(2) = 3 sets cubic interpolation in the y direction.
                                        !! * intpol(3) = 1 sets linear interpolation in the z direction
                                        !! * intpol(3) = 3 sets cubic interpolation in the z direction.
                                        !!
                                        !! values other than 1 or 3 in intpol are not allowed (ier = 3).
    integer(kind=i4),intent(in)     :: lw        !! the integer(kind=i4) length of the real(kind=dp) work space w.
                                        !!
                                        !! * let lwx = 1 if mx-1 divides nx-1; otherwise
                                        !! * let lwx = mx if intpol(1) = 1 or
                                        !! * let lwx = 4*mx if intpol(1) = 3
                                        !! * let lwy = 0 if my-1 divides ny-1; otherwise
                                        !! * let lwy = my+2*mx if intpol(2) = 1 or
                                        !! * let lwy = 4*(mx+my) if intpol(2) = 3
                                        !! * let lwz = 0 if mz-1 divides nz-1; otherwise
                                        !! * let lwz = 2*mx*my+mz if intpol(3) = 1 or
                                        !! * let lwz = 4*(mx*my+mz) if intpol(3) = 3
                                        !!
                                        !! then lw must be greater than or equal to lwx+lwy+lwz
    integer(kind=i4),intent(in)     :: liw       !! the integer(kind=i4) length of the integer(kind=i4) work space iw.
                                        !! liw must be greater than or equal to mx+my+mz
    integer(kind=i4),intent(inout)  :: iw(liw)   !! an integer(kind=i4) work space of length at least liw
                                        !! which must be provided in the
                                        !! routine calling rgrd3u
    integer(kind=i4),intent(out)    :: ier       !! an integer(kind=i4) error flag set as follows:
                                        !!
                                        !! * ier = 0 if no errors in input arguments are detected
                                        !! * ier = 1 if  min(mx,my,mz) < 2
                                        !! * ier = 2 if nx < 2 when intpol(1)=1 or nx < 4 when intpol(1)=3 (or)
                                        !!   ny < 2 when intpol(2)=1 or ny < 4 when intpol(2)=3 (or)
                                        !!   nz < 2 when intpol(3)=1 or nz < 4 when intpol(3)=3.
                                        !! * ier = 3 if any of intpol(1),intpol(2),intpol(3) is not equal to 1 or 3
                                        !! * ier = 4 if lw or liw is too small (insufficient work space)
    real(kind=dp),intent(in)    :: p(nx,ny,nz)   !! a real(kind=dp) nx by ny by nz array of given values
    real(kind=dp),intent(out)   :: q(mx,my,mz)   !! a real(kind=dp) mx by my by mz array of values which are interpolated from p.
    real(kind=dp),intent(inout) :: w(lw) !! a real(kind=dp) work space of length at least lw
                                    !! which must be provided in the
                                    !! routine calling rgrd3u

    integer(kind=i4)  :: inmx,jnmy,knmz,isubx,jsuby,ksubz
    integer(kind=i4)  :: lwx,lwy,lwz,mxmy,jy,kz
    integer(kind=i4)  :: i2,i3,i4,i5
    integer(kind=i4)  :: j2,j3,j4,j5,j6,j7,j8,j9
    integer(kind=i4)  :: k2,k3,k4,k5,k6,k7,k8,k9

    ! check input arguments

    ! check mx,my,mz
    ier = 1
    if (min(mx,my,mz) < 1) return

    ! check intpol
    ier = 3
    if (intpol(1)/=1 .and. intpol(1)/=3) return
    if (intpol(2)/=1 .and. intpol(2)/=3) return
    if (intpol(3)/=1 .and. intpol(3)/=3) return

    ! check nx,ny,nz
    ier = 2
    if (intpol(1)==1 .and. nx<2) return
    if (intpol(1)==3 .and. nx<4) return
    if (intpol(2)==1 .and. ny<2) return
    if (intpol(2)==3 .and. ny<4) return
    if (intpol(3)==1 .and. nz<2) return
    if (intpol(3)==3 .and. nz<4) return

    ! set subgrid indicators
    inmx = (nx-1)/(mx-1)
    jnmy = (ny-1)/(my-1)
    knmz = (nz-1)/(mz-1)
    isubx = nx - inmx*(mx-1)
    jsuby = ny - jnmy*(my-1)
    ksubz = nz - knmz*(mz-1)

    ! check work space lengths
    ier = 4
    mxmy = mx*my
    lwx = 1
    if (isubx/=1) then
        if (intpol(1)==1) then
            lwx = mx
        else
            lwx = 4*mx
        end if
    end if
    lwy = 0
    if (jsuby/=1) then
        if (intpol(2)==1) then
            lwy = (2*mx+my)
        else
            lwy = 4*(mx+my)
        end if
    end if
    lwz = 0
    if (ksubz/=1) then
        if (intpol(3)==1) then
            lwz = (2*mxmy+mz)
        else
            lwz = 4*(mxmy+mz)
        end if
    end if
    if (lw < lwx+lwy+lwz) return
    if (liw < mx+my+mz) return

    ! arguments o.k.

    ier = 0
    jy = mx+1
    kz = mx+my+1

    ! preset work space pointers

    k2 = 1
    k3 = 1
    k4 = 1
    k5 = 1
    k6 = 1
    k7 = 1
    k8 = 1
    k9 = 1
    j2 = 1
    j3 = 1
    j4 = 1
    j5 = 1
    j6 = 1
    j7 = 1
    j8 = 1
    j9 = 1
    i2 = 1
    i3 = 1
    i4 = 1
    i5 = 1

    if (intpol(3)==1) then
        if (ksubz/=1) then
            ! linearly interpolate in nz, set work space pointers
            k2 = 1
            k3 = k2
            k4 = k3+mz
            k5 = k4
            k6 = k5
            k7 = k6
            k8 = k7+mxmy
            k9 = k8+mxmy
            ! set z interpolation indices and scales
            call linmxu(nz,mz,iw(kz),w(k3))
            j2 = k9
            i2 = k9
        end if

        if (jsuby/=1) then
            if (intpol(2) == 1) then
                ! linear in y
                j3 = j2
                j4 = j3+my
                j5 = j4
                j6 = j5
                j7 = j6
                j8 = j7+mx
                j9 = j8+mx
                call linmxu(ny,my,iw(jy),w(j3))
                i2 = j9
            else
                ! cubic in y
                j3 = j2+my
                j4 = j3+my
                j5 = j4+my
                j6 = j5+my
                j7 = j6+mx
                j8 = j7+mx
                j9 = j8+mx
                call cubnmxu(ny,my,iw(jy),w(j2),w(j3),w(j4),w(j5))
                i2 = j9+mx
            end if
        end if

        if (isubx/=1) then
            if (intpol(1) == 1) then
                ! linear in x
                i3 = i2
                i4 = i3
                i5 = i4
                call linmxu(nx,mx,iw,w(i3))
            else
                ! cubic in x
                i3 = i2+mx
                i4 = i3+mx
                i5 = i4+mx
                call cubnmxu(nx,mx,iw,w(i2),w(i3),w(i4),w(i5))
            end if
        end if

        ! linearly interpolate p onto q in z

        call lint3u(nx,ny,nz,p,mx,my,mxmy,mz,q,intpol,iw(kz),w(k3),&
                    w(k7),w(k8),iw(jy),w(j2),w(j3),w(j4),w(j5),w(j6),w(j7),&
                    w(j8),w(j9),iw,w(i2),w(i3),w(i4),w(i5),&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)

    else

        ! cubically interpolate in z

        if (ksubz/=1) then
            k2 = 1
            k3 = k2+mz
            k4 = k3+mz
            k5 = k4+mz
            k6 = k5+mz
            k7 = k6+mxmy
            k8 = k7+mxmy
            k9 = k8+mxmy
            call cubnmxu(nz,mz,iw(kz),w(k2),w(k3),w(k4),w(k5))
            j2 = k9+mxmy
            i2 = j2
        end if

        if (jsuby/=1) then
            if (intpol(2) == 1) then
                j3 = j2
                j4 = j3+my
                j5 = j4
                j6 = j5
                j7 = j6
                j8 = j7+mx
                j9 = j8+mx
                call linmxu(ny,my,iw(jy),w(j3))
                i2 = j9
            else
                j3 = j2+my
                j4 = j3+my
                j5 = j4+my
                j6 = j5+my
                j7 = j6+mx
                j8 = j7+mx
                j9 = j8+mx
                call cubnmxu(ny,my,iw(jy),w(j2),w(j3),w(j4),w(j5))
                i2 = j9+mx
            end if
        end if

        if (isubx/=1) then
            if (intpol(1) == 1) then
                i3 = i2
                i4 = i3
                i5 = i4
                call linmxu(nx,mx,iw,w(i3))
            else
                i3 = i2+mx
                i4 = i3+mx
                i5 = i4+mx
                call cubnmxu(nx,mx,iw,w(i2),w(i3),w(i4),w(i5))
            end if
        end if

        ! cubically interpolate p onto q in z

        call cubt3u(nx,ny,nz,p,mx,my,mxmy,mz,q,intpol,&
                    iw(kz),w(k2),w(k3),w(k4),w(k5),w(k6),w(k7),w(k8),w(k9),&
                    iw(jy),w(j2),w(j3),w(j4),w(j5),w(j6),w(j7),w(j8),w(j9),&
                    iw,w(i2),w(i3),w(i4),w(i5),&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)

    end if

    end subroutine rgrd3u
!**************************************************************************

!**************************************************************************
!>
!  linearly interpolate in z direction

    subroutine lint3u(nx,ny,nz,p,mx,my,mxmy,mz,q,intpol,kz,&
                      dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,&
                      pjpp,ix,dxm,dx,dxp,dxpp,&
                      inmx,jnmy,knmz,isubx,jsuby,ksubz)
       !dir$ attributes code_align : 32 :: lint3u
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter: "TARGET_ARCH=skylake_avx512" :: lint3u
    implicit none

    integer(kind=i4)  :: nx,ny,nz,mx,my,mz,mxmy,intpol(3),kz(mz),jy(my),ix(mx)
    integer(kind=i4)  :: inmx,jnmy,knmz,isubx,jsuby,ksubz
    integer(kind=i4)  :: kk,k,iijj,ksave
    real(kind=dp) :: p(nx,ny,nz),q(mxmy,mz)
    real(kind=dp) :: dz(mz),pk(mxmy),pkp(mxmy)
    real(kind=dp) :: dym(my),dy(my),dyp(my),dypp(my)
    real(kind=dp) :: pjm(mx),pj(mx),pjp(mx),pjpp(mx)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)

    if (intpol(2) == 1) then

        ! linear in y

        if (ksubz == 1) then
            ! mz grid is subset of nz grid
            do kk=1,mz
                k = knmz*(kk-1)+1
                call lint2u(nx,ny,p(1,1,k),mx,my,q(1,kk),intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
            end do
            return
        end if

        ksave = -1
        do kk=1,mz
            k = kz(kk)
            if (k==ksave) then
                ! k pointer has not moved since last pass (no updates or interpolation)
            else if (k==ksave+1) then
                ! update k and interpolate k+1
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pk(iijj) = pkp(iijj)
                end do
                call lint2u(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
            else
                ! interpolate k,k+1 in pk,pkp
                call lint2u(nx,ny,p(1,1,k),mx,my,pk,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
                call lint2u(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
            end if

            ! save k pointer for next pass

            ksave = k

            ! linearly interpolate q(ii,jj,k) from pk,pkp in z direction
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned q:64
                !dir$ assume_aligned dz:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ assume (mod(mz,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
            do iijj=1,mxmy
                q(iijj,kk) = pk(iijj)+dz(kk)*(pkp(iijj)-pk(iijj))
            end do
        end do

    else

        ! cubic in y

        if (ksubz == 1) then
            ! mz grid is subset of nz grid
            do kk=1,mz
                k = knmz*(kk-1)+1
                call cubt2u(nx,ny,p(1,1,k),mx,my,q(1,kk),intpol,&
                    jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,isubx,jsuby)
            end do
            return
        end if

        ksave = -1
        do kk=1,mz
            k = kz(kk)
            if (k==ksave) then
                ! k pointer has not moved since last pass (no updates or interpolation)
            else if (k==ksave+1) then
                ! update k and interpolate k+1
                 !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pk(iijj) = pkp(iijj)
                end do
                call cubt2u(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,&
                    jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,isubx,jsuby)
            else
                ! interpolate k,k+1 in pk,pkp
                call cubt2u(nx,ny,p(1,1,k),mx,my,pk,intpol,&
                    jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,isubx,jsuby)
                call cubt2u(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,&
                    jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,isubx,jsuby)
            end if

            ! save k pointer for next pass

            ksave = k

            ! linearly interpolate q(ii,jj,k) from pk,pkp in z direction
               !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned q:64
                !dir$ assume_aligned dz:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ assume (mod(mz,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
            do iijj=1,mxmy
                q(iijj,kk) = pk(iijj)+dz(kk)*(pkp(iijj)-pk(iijj))
            end do
        end do

    end if

    end subroutine lint3u
!**************************************************************************

!**************************************************************************
!>
!  cubically interpolate in z

    subroutine cubt3u(nx,ny,nz,p,mx,my,mxmy,mz,q,intpol,&
                      kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,&
                      pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                      inmx,jnmy,knmz,isubx,jsuby,ksubz)
       !dir$ attributes code_align : 32 :: cubt3u
       !dir$ optimize : 3
       !dir$ attributes optimization_parameter:"TARGET_ARCH=skylake_avx512" :: cubt3u
    implicit none

    integer(kind=i4)  :: nx,ny,nz,mx,my,mz,mxmy,intpol(3),kz(mz),jy(my),ix(mx)
    integer(kind=i4)  :: inmx,jnmy,knmz,isubx,jsuby,ksubz
    integer(kind=i4)  :: kk,k,iijj,ksave
    real(kind=dp) :: p(nx,ny,nz),q(mxmy,mz)
    real(kind=dp) :: dzm(mz),dz(mz),dzp(mz),dzpp(mz)
    real(kind=dp) :: pkm(mxmy),pk(mxmy),pkp(mxmy),pkpp(mxmy)
    real(kind=dp) :: dym(my),dy(my),dyp(my),dypp(my)
    real(kind=dp) :: pjm(mx),pj(mx),pjp(mx),pjpp(mx)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)

    if (intpol(2) == 1) then

        ! linear in y
        if (ksubz == 1) then
            ! mz grid is subset of nz grid
            do kk=1,mz
                k = knmz*(kk-1)+1
                call lint2u(nx,ny,p(1,1,k),mx,my,q(1,kk),intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
            end do
            return
        end if

        ! mz not a subgrid of nz

        ksave = -3
        do kk=1,mz
            k = kz(kk)
            if (k==ksave) then
                ! k pointer has not moved since last pass (no updates or interpolation)
            else if (k==ksave+1) then
                ! update k-1,k,k+1 and interpolate k+2
                 !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pk(iijj)
                    pk(iijj) = pkp(iijj)
                    pkp(iijj) = pkpp(iijj)
                end do
                call lint2u(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
            else if (k==ksave+2) then
                ! update k-1,k and interpolate k+1,k+2
                 !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pkp(iijj)
                    pk(iijj) = pkpp(iijj)
                end do
                call lint2u(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
                call lint2u(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
            else if (k==ksave+3) then
                ! update k-1 and interpolate k,k+1,k+2
                 !dir$ assume_aligned pkm:64
                 !dir$ assume_aligned pkp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pkpp(iijj)
                end do
                call lint2u(nx,ny,p(1,1,k),mx,my,pk,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
                call lint2u(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
                call lint2u(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
            else
                ! interpolate all four k-1,k,k+1,k+2
                call lint2u(nx,ny,p(1,1,k-1),mx,my,pkm,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
                call lint2u(nx,ny,p(1,1,k),mx,my,pk,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
                call lint2u(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
                call lint2u(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,jy,dy,pj,pjp,ix,dxm,dx,dxp,dxpp,inmx,jnmy,isubx,jsuby)
            end if

            ! save k pointer for next pass

            ksave = k

            ! cubically interpolate q(ii,jj,kk) from pkm,pk,pkp,pkpp in z direction
                !dir$ assume_aligned q:64
                !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ assume_aligned dzm:64
                !dir$ assume_aligned dz:64
                !dir$ assume_aligned dzp:64
                !dir$ assume_aligned dzpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ assume (mod(mz,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8) 
            do iijj=1,mxmy
                q(iijj,kk) = dzm(kk)*pkm(iijj) + dz(kk)*pk(iijj) + dzp(kk)*pkp(iijj) + dzpp(kk)*pkpp(iijj)
            end do
        end do

    else

        ! cubic in y

        if (ksubz == 1) then
            ! mz grid is subset of nz grid
            do kk=1,mz
                k = knmz*(kk-1)+1
                call cubt2u(nx,ny,p(1,1,k),mx,my,q(1,kk),intpol,jy,dym,dy,&
                    dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,isubx,jsuby)

            end do
            return
        end if

        ksave = -3
        do kk=1,mz
            k = kz(kk)
            if (k==ksave) then
                ! k pointer has not moved since last pass (no updates or interpolation)
            else if (k==ksave+1) then
                ! update k-1,k,k+1 and interpolate k+2
                !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pk(iijj)
                    pk(iijj) = pkp(iijj)
                    pkp(iijj) = pkpp(iijj)
                end do
                call cubt2u(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,isubx,jsuby)
            else if (k==ksave+2) then
                ! update k-1,k and interpolate k+1,k+2
               !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pkp(iijj)
                    pk(iijj) = pkpp(iijj)
                end do
                call cubt2u(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,isubx,jsuby)
                call cubt2u(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,isubx,jsuby)
            else if (k==ksave+3) then
                ! update k-1 and interpolate k,k+1,k+2
                !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pkpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8)
                do iijj=1,mxmy
                    pkm(iijj) = pkpp(iijj)
                end do
                call cubt2u(nx,ny,p(1,1,k),mx,my,pk,intpol,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,isubx,jsuby)
                call cubt2u(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,isubx,jsuby)
                call cubt2u(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,isubx,jsuby)
            else
                ! interpolate all four k-1,k,k+1,k+2
                call cubt2u(nx,ny,p(1,1,k-1),mx,my,pkm,intpol,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,isubx,jsuby)
                call cubt2u(nx,ny,p(1,1,k),mx,my,pk,intpol,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,isubx,jsuby)
                call cubt2u(nx,ny,p(1,1,k+1),mx,my,pkp,intpol,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,isubx,jsuby)
                call cubt2u(nx,ny,p(1,1,k+2),mx,my,pkpp,intpol,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,isubx,jsuby)
            end if

            ! save k pointer for next pass

            ksave = k

            ! cubically interpolate q(ii,jj,kk) from pkm,pk,pkp,pkpp in z direction
               !dir$ assume_aligned q:64
                !dir$ assume_aligned pkm:64
                !dir$ assume_aligned pk:64
                !dir$ assume_aligned pkp:64
                !dir$ assume_aligned pkpp:64
                !dir$ assume_aligned dzm:64
                !dir$ assume_aligned dz:64
                !dir$ assume_aligned dzp:64
                !dir$ assume_aligned dzpp:64
                !dir$ vector aligned
                !dir$ assume (mod(mxmy,8) .eq. 0)
                !dir$ assume (mod(mz,8) .eq. 0)
                !dir$ vector vectorlength(8)
                !dir$ vector always
                !dir$ unroll(8) 
            do iijj=1,mxmy
                q(iijj,kk) = dzm(kk)*pkm(iijj) + dz(kk)*pk(iijj) + dzp(kk)*pkp(iijj) + dzpp(kk)*pkpp(iijj)
            end do
        end do

    end if

    end subroutine cubt3u
!**************************************************************************

!**************************************************************************
!>
!  subroutine rgrd4 interpolates the values p(i,j,k,l) on the orthogonal
!  grid (x(i),y(j),z(k),t(l)) for i=1,...,nx;j=1,...,ny;k=1,...,nz;l=1,...,nt
!  onto q(ii,jj,kk,ll) on the orthogonal grid (xx(ii),yy(jj),zz(kk),tt(ll))
!  for ii=1,...,mx;jj=1,...,my;kk=1,...,mz;ll=1,...,mt
!
!### method
!
!  linear or cubic interpolation is used (independently) in
!  each direction (see argument intpol).
!
!### requirements
!
!  each of the x,y,z,t grids must be strictly montonically increasing
!  and each of the xx,yy,zz,tt grids must be montonically increasing
!  (see ier = 4).  in addition the (X,Y,Z,T) region of the q grid
!
!    [xx(1),xx(mx)] X [yy(1),yy(my)] X [zz(1),zz(mz)] X [tt(1),tt(my)]
!
!  must lie within the (X,Y,Z,T) region of the p grid
!
!    [x(1),x(nx)] X [y(1),y(ny)] X [z(1),z(nz)] X [t(1),t(nt)].
!
!  extrapolation is not allowed (see ier=3).  if these (X,Y,Z,T)
!  regions are identical and the orthogonal grids are UNIFORM
!  in each direction then subroutine rgrd4u
!  should be used instead of rgrd4.
#if 0
    subroutine rgrd4(nx,ny,nz,nt,x,y,z,t,p,mx,my,mz,mt,xx,yy,zz,tt,q,intpol,w,lw,iw,liw,ier)

    implicit none

    integer(kind=i4),intent(in)  :: nx   !! the integer(kind=i4) dimension of the grid vector x and the first dimension of p.
                                !! nx > 1 if intpol(1) = 1 or nx > 3 if intpol(1) = 3 is required.
    integer(kind=i4),intent(in)  :: ny   !! the integer(kind=i4) dimension of the grid vector y and the second dimension of p.
                                !! ny > 1 if intpol(2) = 1 or ny > 3 if intpol(2) = 3 is required.
    integer(kind=i4),intent(in)  :: nz   !! the integer(kind=i4) dimension of the grid vector z and the third dimension of p.
                                !! nz > 1 if intpol(3) = 1 or nz > 3 if intpol(3) = 3 is required.
    integer(kind=i4),intent(in)  :: nt   !! the integer(kind=i4) dimension of the grid vector t and the fourth dimension of p.
                                !! nt > 1 if intpol(4) = 1 or nt > 3 if intpol(4) = 3 is required.
    integer(kind=i4),intent(in)  :: mx   !! the integer(kind=i4) dimension of the grid vector xx and the first dimension
                                !! of q.  mx > 0 is required.
    integer(kind=i4),intent(in)  :: my   !! the integer(kind=i4) dimension of the grid vector yy and the second dimension
                                !! of q.  my > 0 is required.
    integer(kind=i4),intent(in)  :: mz   !! the integer(kind=i4) dimension of the grid vector zz and the third dimension
                                !! of q.  mz > 0 is required.
    integer(kind=i4),intent(in)  :: mt   !! the integer(kind=i4) dimension of the grid vector tt and the fourth dimension
                                !! of q.  mt > 0 is required.
    integer(kind=i4),intent(in)  :: lw   !! the integer(kind=i4) length of the real(kind=dp) work space w.  let
                                !!
                                !! * lwx = mx                  if intpol(1) = 1
                                !! * lwx = 4*mx                if intpol(1) = 3
                                !! * lwy = my+2*mx             if intpol(2) = 1
                                !! * lwy = 4*(my+mx)           if intpol(2) = 3
                                !! * lwz = 2*mx*my+mz          if intpol(3) = 1
                                !! * lwz = 4*(mx*my+mz)        if intpol(3) = 3
                                !! * lwt = 2*mx*my*mz+mt       if intpol(4) = 1
                                !! * lwt = 4*(mx*my*mz+mt)     if intpol(4) = 3
                                !!
                                !! then lw must be greater than or equal to lwx+lwy+lwz+lwt
    integer(kind=i4),intent(in)  :: liw  !! the integer(kind=i4) length of the integer(kind=i4) work space iw.  liw must be at least mx+my+mz+mt
    integer(kind=i4),intent(out)  :: ier     !! an integer(kind=i4) error flag set as follows:
                                    !!
                                    !! * ier = 0 if no errors in input arguments are detected
                                    !! * ier = 1 if  min(mx,my,mz,mt) < 1
                                    !! * ier = 2 if nx < 2 when intpol(1)=1 or nx < 4 when intpol(1)=3 (or)
                                    !!    * ny < 2 when intpol(2)=1 or ny < 4 when intpol(2)=3 (or)
                                    !!    * nz < 2 when intpol(3)=1 or nz < 4 when intpol(3)=3 (or)
                                    !!    * nt < 2 when intpol(4)=1 or nt < 4 when intpol(4)=3
                                    !! * ier = 3 if xx(1) < x(1) or x(nx) < xx(mx) (or)
                                    !!    * yy(1) < y(1) or y(ny) < yy(my) (or)
                                    !!    * zz(1) < z(1) or z(nz) < zz(mz) (or)
                                    !!    * tt(1) < t(1) or t(nt) < tt(mt)
                                    !!   to avoid this flag when end points are intended to be the
                                    !!   same but may differ slightly due to roundoff error, they
                                    !!   should be set exactly in the calling routine (e.g., if both
                                    !!   grids have the same y boundaries then yy(1)=y(1) and yy(my)=y(ny)
                                    !!   should be set before calling rgrd4)
                                    !!
                                    !! * ier = 4 if one of the grids x,y,z,t is not strictly monotonically
                                    !!   increasing or if one of the grids xx,yy,zz,tt is not
                                    !!   montonically increasing.  more precisely if:
                                    !!
                                    !!    * x(i+1) <= x(i) for some i such that 1 <= i < nx (or)
                                    !!    * y(j+1) <= y(j) for some j such that 1 <= j < ny (or)
                                    !!    * z(k+1) <= z(k) for some k such that 1 <= k < nz (or)
                                    !!    * t(l+1) <= t(l) for some l such that 1 <= l < nt (or)
                                    !!    * xx(ii+1) < xx(ii) for some ii such that 1 <= ii < mx (or)
                                    !!    * yy(jj+1) < yy(jj) for some jj such that 1 <= jj < my (or)
                                    !!    * zz(kk+1) < zz(k)  for some kk such that 1 <= kk < mz (or)
                                    !!    * tt(ll+1) < tt(l)  for some ll such that 1 <= ll < mt
                                    !! * ier = 5 if lw or liw is too small (insufficient work space)
                                    !! * ier = 6 if any of intpol(1),intpol(2),intpol(3),intpol(4)
                                    !!                 is not equal to 1 or 3

    integer(kind=i4),intent(inout)  :: iw(liw)   !! an integer(kind=i4) work space of length at least liw
                                        !! which must be provided in the routine calling rgrd4
    integer(kind=i4),intent(in)  :: intpol(4)    !! an integer(kind=i4) vector of dimension 4 which sets linear or cubic
                                        !! interpolation in each of the x,y,z,t directions as follows:
                                        !!
                                        !! * intpol(1) = 1 sets linear interpolation in the x direction
                                        !! * intpol(1) = 3 sets cubic interpolation in the x direction.
                                        !! * intpol(2) = 1 sets linear interpolation in the y direction
                                        !! * intpol(2) = 3 sets cubic interpolation in the y direction.
                                        !! * intpol(3) = 1 sets linear interpolation in the z direction
                                        !! * intpol(3) = 3 sets cubic interpolation in the z direction.
                                        !! * intpol(4) = 1 sets linear interpolation in the t direction
                                        !! * intpol(4) = 3 sets cubic interpolation in the t direction.
                                        !!
                                        !! values other than 1 or 3 in intpol are not allowed (ier = 6).
    real(kind=dp),intent(in) :: x(nx)    !! a real(kind=dp) nx vector of strictly increasing values which defines the x
                                    !! portion of the orthogonal grid on which p is given
    real(kind=dp),intent(in) :: y(ny)    !! a real(kind=dp) ny vector of strictly increasing values which defines the y
                                    !! portion of the orthogonal grid on which p is given
    real(kind=dp),intent(in) :: z(nz)    !! a real(kind=dp) nz vector of strictly increasing values which defines the z
                                    !! portion of the orthogonal grid on which p is given
    real(kind=dp),intent(in) :: t(nt)    !! a real(kind=dp) nt vector of strictly increasing values which defines the t
                                    !! portion of the orthogonal grid on which p is given
    real(kind=dp),intent(in) :: p(nx,ny,nz,nt)   !! a real(kind=dp) nx by ny by nz by nt array of values given on the (x,y,z,t) grid
    real(kind=dp),intent(inout) :: w(lw) !! a real(kind=dp) work space of length at least lw
                                    !! which must be provided in the routine calling rgrd4
    real(kind=dp),intent(in) :: xx(mx)   !! a real(kind=dp) mx vector of increasing values which defines the x portion of the
                                    !! orthogonal grid on which q is defined.  xx(1) < x(1) or xx(mx) > x(nx)
                                    !! is not allowed (see ier = 3)
    real(kind=dp),intent(in) :: yy(my)   !! a real(kind=dp) my vector of increasing values which defines the y portion of the
                                    !! orthogonal grid on which q is defined.  yy(1) < y(1) or yy(my) > y(ny)
                                    !! is not allowed (see ier = 3)
    real(kind=dp),intent(in) :: zz(mz)   !! a real(kind=dp) mz vector of increasing values which defines the z portion of the
                                    !! orthogonal grid on which q is defined.  zz(1) < z(1) or zz(mz) > z(nz)
                                    !! is not allowed (see ier = 3)
    real(kind=dp),intent(in) :: tt(mt)   !! a real(kind=dp) mt vector of increasing values which defines the t portion of the
                                    !! orthogonal grid on which q is defined.  tt(1) < t(1) or tt(mt) > t(nt)
                                    !! is not allowed (see ier = 3)
    real(kind=dp),intent(out) :: q(mx,my,mz,mt)  !! a real(kind=dp) mx by my by mz by mt array of values on the (xx,yy,zz,tt) grid
                                            !! which are interpolated from p on the (x,y,z,t) grid

    integer(kind=i4)  :: l2,l3,l4,l5,l6,l7,l8,l9
    integer(kind=i4)  :: k2,k3,k4,k5,k6,k7,k8,k9
    integer(kind=i4)  :: j2,j3,j4,j5,j6,j7,j8,j9
    integer(kind=i4)  :: i2,i3,i4,i5
    integer(kind=i4)  :: lwx,lwy,lwz,lwt,mxmy,mxmymz
    integer(kind=i4)  :: ii,jj,kk,ll,i,j,k,l
    integer(kind=i4)  :: jy,kz,lt

    ! check input arguments

    ! check (xx,yy,zz,tt) grid resolution
    ier = 1
    if (min(mx,my,mz,mt) < 1) return

    ! check intpol
    ier = 6
    if (intpol(1)/=1 .and. intpol(1)/=3) return
    if (intpol(2)/=1 .and. intpol(2)/=3) return
    if (intpol(3)/=1 .and. intpol(3)/=3) return
    if (intpol(4)/=1 .and. intpol(4)/=3) return

    ! check (x,y,z,t) grid resolution
    ier = 2
    if (intpol(1)==1 .and. nx<2) return
    if (intpol(1)==3 .and. nx<4) return
    if (intpol(2)==1 .and. ny<2) return
    if (intpol(2)==3 .and. ny<4) return
    if (intpol(3)==1 .and. nz<2) return
    if (intpol(3)==3 .and. nz<4) return
    if (intpol(4)==1 .and. nt<2) return
    if (intpol(4)==3 .and. nt<4) return

    ! check work space length input and set minimum
    ier = 5
    mxmy = mx*my
    mxmymz = mxmy*mz
    if (intpol(1)==1) then
        lwx = mx
    else
        lwx = 4*mx
    end if
    if (intpol(2)==1) then
        lwy = (my+2*mx)
    else
        lwy = 4*(mx+my)
    end if
    if (intpol(3)==1) then
        lwz = (2*mxmy+mz)
    else
        lwz = 4*(mxmy+mz)
    end if
    if (intpol(4)==1) then
        lwt = (2*mxmymz+mt)
    else
        lwt = 4*(mxmymz+mt)
    end if
    if (lw < lwx+lwy+lwz+lwt) return
    if (liw < mx+my+mz+mt) return

    ! check (xx,yy,zz,tt) grid contained in (x,y,z,t) grid
    ier = 3
    if (xx(1)<x(1) .or. xx(mx)>x(nx)) return
    if (yy(1)<y(1) .or. yy(my)>y(ny)) return
    if (zz(1)<z(1) .or. zz(mz)>z(nz)) return
    if (tt(1)<t(1) .or. tt(mt)>t(nt)) return

    ! check montonicity of grids
    ier = 4
    do i=2,nx
        if (x(i-1)>=x(i)) return
    end do
    do j=2,ny
        if (y(j-1)>=y(j)) return
    end do
    do k=2,nz
        if (z(k-1)>=z(k)) return
    end do
    do l=2,nt
        if (t(l-1)>=t(l)) return
    end do
    do ii=2,mx
        if (xx(ii-1)>xx(ii)) return
    end do
    do jj=2,my
        if (yy(jj-1)>yy(jj)) return
    end do
    do kk=2,mz
        if (zz(kk-1)>zz(kk)) return
    end do
    do ll=2,mt
        if (tt(ll-1)>tt(ll)) return
    end do

    ! arguments o.k.

    ier = 0

    ! set pointers for integer(kind=i4) work space iw

    jy = mx+1
    kz = mx+my+1
    lt = mx+my+mz+1

    if (intpol(4)==1) then

        ! linearly interpolate in nt, set work space pointers and scales

        l2 = 1
        l3 = l2
        l4 = l3+mt
        l5 = l4
        l6 = l5
        l7 = l6
        l8 = l7+mxmymz
        l9 = l8+mxmymz
        call linmx(nt,t,mt,tt,iw(lt),w(l3))
        k2 = l9

        if (intpol(3)==1) then
            ! linear in z
            k3 = k2
            k4 = k3+mz
            k5 = k4
            k6 = k5
            k7 = k6
            k8 = k7+mxmy
            k9 = k8+mxmy
            call linmx(nz,z,mz,zz,iw(kz),w(k3))
            j2 = k9
        else
            ! cubic in z
            k3 = k2+mz
            k4 = k3+mz
            k5 = k4+mz
            k6 = k5+mz
            k7 = k6+mxmy
            k8 = k7+mxmy
            k9 = k8+mxmy
            call cubnmx(nz,z,mz,zz,iw(kz),w(k2),w(k3),w(k4),w(k5))
            j2 = k9+mxmy
        end if

        if (intpol(2) == 1) then
            ! linear in y
            j3 = j2
            j4 = j3+my
            j5 = j4
            j6 = j5
            j7 = j6
            j8 = j7+mx
            j9 = j8+mx
            call linmx(ny,y,my,yy,iw(jy),w(j3))
            i2 = j9
        else
            ! cubic in y
            j3 = j2+my
            j4 = j3+my
            j5 = j4+my
            j6 = j5+my
            j7 = j6+mx
            j8 = j7+mx
            j9 = j8+mx
            call cubnmx(ny,y,my,yy,iw(jy),w(j2),w(j3),w(j4),w(j5))
            i2 = j9+mx
        end if

        if (intpol(1) == 1) then
            ! linear in x
            i3 = i2
            i4 = i3
            i5 = i4
            call linmx(nx,x,mx,xx,iw,w(i3))
        else
            ! cubic in x
            i3 = i2+mx
            i4 = i3+mx
            i5 = i4+mx
            call cubnmx(nx,x,mx,xx,iw,w(i2),w(i3),w(i4),w(i5))
        end if

        ! linearly interpolate in t

        call lint4(nx,ny,nz,nt,p,mx,my,mz,mt,mxmy,mxmymz,q,intpol,&
                    iw(lt),w(l3),w(l7),w(l8),&
                    iw(kz),w(k2),w(k3),w(k4),w(k5),w(k6),w(k7),w(k8),w(k9),&
                    iw(jy),w(j2),w(j3),w(j4),w(j5),w(j6),w(j7),w(j8),w(j9),&
                    iw,w(i2),w(i3),w(i4),w(i5))

    else

        ! cubically interpolate in t

        l2 = 1
        l3 = l2+mt
        l4 = l3+mt
        l5 = l4+mt
        l6 = l5+mt
        l7 = l6+mxmymz
        l8 = l7+mxmymz
        l9 = l8+mxmymz
        call cubnmx(nt,t,mt,tt,iw(lt),w(l2),w(l3),w(l4),w(l5))
        k2 = l9+mxmymz

        if (intpol(3)==1) then
            ! linear in z
            k3 = k2
            k4 = k3+mz
            k5 = k4
            k6 = k5
            k7 = k6
            k8 = k7+mxmy
            k9 = k8+mxmy
            call linmx(nz,z,mz,zz,iw(kz),w(k3))
            j2 = k9
        else
            ! cubic in z
            k3 = k2+mz
            k4 = k3+mz
            k5 = k4+mz
            k6 = k5+mz
            k7 = k6+mxmy
            k8 = k7+mxmy
            k9 = k8+mxmy
            call cubnmx(nz,z,mz,zz,iw(kz),w(k2),w(k3),w(k4),w(k5))
            j2 = k9+mxmy
        end if

        if (intpol(2) == 1) then
            j3 = j2
            j4 = j3+my
            j5 = j4
            j6 = j5
            j7 = j6
            j8 = j7+mx
            j9 = j8+mx
            call linmx(ny,y,my,yy,iw(jy),w(j3))
            i2 = j9
        else
            j3 = j2+my
            j4 = j3+my
            j5 = j4+my
            j6 = j5+my
            j7 = j6+mx
            j8 = j7+mx
            j9 = j8+mx
            call cubnmx(ny,y,my,yy,iw(jy),w(j2),w(j3),w(j4),w(j5))
            i2 = j9+mx
        end if

        ! set work space portion and indices which depend on x interpolation

        if (intpol(1) == 1) then
            i3 = i2
            i4 = i3
            i5 = i4
            call linmx(nx,x,mx,xx,iw,w(i3))
        else
            i3 = i2+mx
            i4 = i3+mx
            i5 = i4+mx
            call cubnmx(nx,x,mx,xx,iw,w(i2),w(i3),w(i4),w(i5))
        end if

        ! cubically interpolate in t

        call cubt4(nx,ny,nz,nt,p,mx,my,mz,mt,mxmy,mxmymz,q,intpol,&
                    iw(lt),w(l2),w(l3),w(l4),w(l5),w(l6),w(l7),w(l8),w(l9),&
                    iw(kz),w(k2),w(k3),w(k4),w(k5),w(k6),w(k7),w(k8),w(k9),&
                    iw(jy),w(j2),w(j3),w(j4),w(j5),w(j6),w(j7),w(j8),w(j9),&
                    iw,w(i2),w(i3),w(i4),w(i5))

    end if

    end subroutine rgrd4
#endif
!**************************************************************************

!**************************************************************************
!>
!  linearly interpolate in t direction
#if 0
    subroutine lint4(nx,ny,nz,nt,p,mx,my,mz,mt,mxmy,mxmymz,q,intpol,&
                     lt,dt,pt,ptp,kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,&
                     jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)

    implicit none

    integer(kind=i4)  :: nx,ny,nz,nt,mx,my,mz,mt,mxmy,mxmymz,lsave,ll,l,iijjkk
    integer(kind=i4)  :: lt(mt),kz(mz),jy(my),ix(mx),intpol(4)
    real(kind=dp) :: p(nx,ny,nz,nt),q(mxmymz,mt)
    real(kind=dp) :: dt(mt),pt(mxmymz),ptp(mxmymz)
    real(kind=dp) :: dzm(mz),dz(mz),dzp(mz),dzpp(mz)
    real(kind=dp) :: pkm(mxmy),pk(mxmy),pkp(mxmy),pkpp(mxmy)
    real(kind=dp) :: dym(my),dy(my),dyp(my),dypp(my)
    real(kind=dp) :: pjm(mx),pj(mx),pjp(mx),pjpp(mx)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)

    if (intpol(3) == 1) then

        ! linear in z

        lsave = -1
        do ll=1,mt
            l = lt(ll)
            if (l==lsave) then
                ! l pointer has not moved since last pass (no updates or interpolation)
            else if (l==lsave+1) then
                ! update l and interpolate l+1
                do iijjkk=1,mxmymz
                    pt(iijjkk) = ptp(iijjkk)
                end do
                call lint3(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,&
                            kz,dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,&
                            dxm,dx,dxp,dxpp)
            else
                ! interpolate l,l+1 in pt,ptp on xx,yy,zz mesh
                call lint3(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call lint3(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            end if

            ! save l pointer for next pass

            lsave = l

            ! linearly interpolate q(ii,jj,,kk,ll) from pt,ptp in t direction

            do iijjkk=1,mxmymz
                q(iijjkk,ll) = pt(iijjkk)+dt(ll)*(ptp(iijjkk)-pt(iijjkk))
            end do
        end do

    else

        ! cubic in z

        lsave = -1
        do ll=1,mt
            l = lt(ll)
            if (l==lsave) then
                ! l pointer has not moved since last pass (no updates or interpolation)
            else if (l==lsave+1) then
                ! update l and interpolate l+1
                do iijjkk=1,mxmymz
                    pt(iijjkk) = ptp(iijjkk)
                end do
                call cubt3(nx,ny,nt,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,&
                            kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            else
                ! interpolate l,l+1 in pt,ptp on xx,yy,zz mesh
                call cubt3(nx,ny,nt,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,&
                            kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call cubt3(nx,ny,nt,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,&
                            kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,&
                            jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            end if

            ! save l pointer for next pass

            lsave = l

            ! linearly interpolate q(ii,jj,kk,ll) from pt,ptp in t direction

            do iijjkk=1,mxmymz
                q(iijjkk,ll) = pt(iijjkk)+dt(ll)*(ptp(iijjkk)-pt(iijjkk))
            end do

        end do

    end if

    end subroutine lint4
#endif
!**************************************************************************

!**************************************************************************
!>
!  cubically interpolate in t
#if 0
    subroutine cubt4(nx,ny,nz,nt,p,mx,my,mz,mt,mxmy,mxmymz,q,intpol,&
                     lt,dtm,dt,dtp,dtpp,ptm,pt,ptp,ptpp,&
                     kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,&
                     jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                     ix,dxm,dx,dxp,dxpp)

    implicit none

    integer(kind=i4)  :: nx,ny,nz,nt,mx,my,mz,mt,mxmy,mxmymz,lsave,ll,l,iijjkk
    integer(kind=i4)  :: lt(mt),kz(mz),jy(my),ix(mx),intpol(4)
    real(kind=dp) :: p(nx,ny,nz,nt),q(mxmymz,mt)
    real(kind=dp) :: dtm(mt),dt(mt),dtp(mt),dtpp(mt)
    real(kind=dp) :: ptm(mxmymz),pt(mxmymz),ptp(mxmymz),ptpp(mxmymz)
    real(kind=dp) :: dzm(mz),dz(mz),dzp(mz),dzpp(mz)
    real(kind=dp) :: pkm(mxmy),pk(mxmy),pkp(mxmy),pkpp(mxmy)
    real(kind=dp) :: dym(my),dy(my),dyp(my),dypp(my)
    real(kind=dp) :: pjm(mx),pj(mx),pjp(mx),pjpp(mx)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)

    if (intpol(3) == 1) then

        ! linear in z

        lsave = -3
        do ll=1,mt
            l = lt(ll)
            if (l==lsave) then
                ! l pointer has not moved since last pass (no updates or interpolation)
            else if (l==lsave+1) then
                ! update l-1,l,l+1 and interpolate l+2
                do iijjkk=1,mxmymz
                    ptm(iijjkk) = pt(iijjkk)
                    pt(iijjkk) = ptp(iijjkk)
                    ptp(iijjkk) = ptpp(iijjkk)
                end do
                call lint3(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            else if (l==lsave+2) then
                ! update l-1,l and interpolate l+1,l+2
                do iijjkk=1,mxmymz
                    ptm(iijjkk) = ptp(iijjkk)
                    pt(iijjkk) = ptpp(iijjkk)
                end do
                call lint3(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call lint3(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            else if (l==lsave+3) then
                ! update l-1 and interpolate l,l+1,l+2

                do iijjkk=1,mxmymz
                    ptm(iijjkk) = ptpp(iijjkk)
                end do
                call lint3(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call lint3(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call lint3(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            else
                ! interpolate all four l-1,l,l+1,l+2
                call lint3(nx,ny,nz,p(1,1,1,l-1),mx,my,mxmy,mz,ptm,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call lint3(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call lint3(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
                call lint3(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp)
            end if

            ! save l pointer for next pass

            lsave = l

            ! cubically interpolate q(ii,jj,kk,ll) from ptm,pt,ptp,ptpp in t direction

            do iijjkk=1,mxmymz
                q(iijjkk,ll) = dtm(ll)*ptm(iijjkk) + dt(ll)*pt(iijjkk) + dtp(ll)*ptp(iijjkk) + dtpp(ll)*ptpp(iijjkk)
            end do
        end do

    else

        ! cubic in z

        lsave = -3
        do ll=1,mt
            l = lt(ll)
            if (l==lsave) then
                ! l pointer has not moved since last pass (no updates or interpolation)
            else if (l==lsave+1) then
                ! update l-1,l,l+1 and interpolate l+2
                do iijjkk=1,mxmymz
                    ptm(iijjkk) = pt(iijjkk)
                    pt(iijjkk) = ptp(iijjkk)
                    ptp(iijjkk) = ptpp(iijjkk)
                end do
                call cubt3(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,dzm,&
                            dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                            ix,dxm,dx,dxp,dxpp)
            else if (l==lsave+2) then
                ! update l-1,l and interpolate l+1,l+2
                do iijjkk=1,mxmymz
                    ptm(iijjkk) = ptp(iijjkk)
                    pt(iijjkk) = ptpp(iijjkk)
                end do
                call cubt3(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,dzm,&
                            dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                            ix,dxm,dx,dxp,dxpp)
                call cubt3(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,dzm,&
                            dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                            ix,dxm,dx,dxp,dxpp)
            else if (l==lsave+3) then
                ! update l-1 and interpolate l,l+1,l+2
                do iijjkk=1,mxmymz
                    ptm(iijjkk) = ptpp(iijjkk)
                end do
                call cubt3(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,kz,dzm,&
                        dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                        ix,dxm,dx,dxp,dxpp)
                call cubt3(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,dzm,&
                        dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                        ix,dxm,dx,dxp,dxpp)
                call cubt3(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,dzm,&
                        dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                        ix,dxm,dx,dxp,dxpp)
            else
                ! interpolate all four l-1,l,l+1,l+2
                call cubt3(nx,ny,nz,p(1,1,1,l-1),mx,my,mxmy,mz,ptm,intpol,kz,dzm,&
                            dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                            ix,dxm,dx,dxp,dxpp)
                call cubt3(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,kz,dzm,&
                            dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                            ix,dxm,dx,dxp,dxpp)
                call cubt3(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,dzm,&
                            dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                            ix,dxm,dx,dxp,dxpp)
                call cubt3(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,dzm,&
                            dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                            ix,dxm,dx,dxp,dxpp)
            end if

            ! save l pointer for next pass

            lsave = l

            ! cubically interpolate q(ii,jj,kk,ll) from ptm,pt,ptp,ptpp in t direction

            do iijjkk=1,mxmymz
                q(iijjkk,ll) = dtm(ll)*ptm(iijjkk) + dt(ll)*pt(iijjkk) + dtp(ll)*ptp(iijjkk) + dtpp(ll)*ptpp(iijjkk)
            end do
        end do

    end if

    end subroutine cubt4
#endif
!**************************************************************************

!**************************************************************************
!>
!  subroutine rgrd4u interpolates the nx by ny by nz by nt array p onto
!  the mx by my by mz by mt array q.  it is assumed that p and q are
!  values on uniform nx by ny by nz by nt and mx by my by mz by mt grids
!  which are superimposed on the same box region (INCLUDING BOUNDARIES).
!  if p and q are values on nonuniform orthogonal grids and/or if the grid
!  on which q is defined lies within the p grid then subroutine rgrd4
!  should be used.
!
!### method
!
!  linear or cubic interpolation (see intpol) is used in each
!  direction for which the q grid is not a subgrid of the p grid.
!  [the mx (my,mz,mt) uniform grid is a subgrid of the nx (ny,nz,nt)
!  uniform grid if and only if mx-1 (my-1,nz-1,nt-1) divides nx-1
!  (ny-1,nz-1,nt-1)].  Values are set directly without (the need for)
!  interpolation in subgrid directions.
#if 0
    subroutine rgrd4u(nx,ny,nz,nt,p,mx,my,mz,mt,q,intpol,w,lw,iw,liw,ier)

    implicit none

    integer(kind=i4),intent(in)      :: nx   !! the integer(kind=i4) first dimension of p.  nx > 1 if intpol(1) = 1 or
                                    !! nx > 3 if intpol(1) = 3 is required (see ier = 2).
    integer(kind=i4),intent(in)      :: ny   !! the integer(kind=i4) second dimension of p.  ny > 1 if intpol(2) = 1 or
                                    !! ny > 3 if intpol(2) = 3 is required (see ier = 2).
    integer(kind=i4),intent(in)      :: nz   !! the integer(kind=i4) third dimension of p.  nz > 1 if intpol(3) = 1 or
                                    !! nz > 3 if intpol(3) = 3 is required (see ier = 2)
    integer(kind=i4),intent(in)      :: nt   !! the integer(kind=i4) fourth dimension of p.  nt > 1 if intpol(4) = 1 or
                                    !! nt > 3 if intpol(4) = 3 is required (see ier=2)
    integer(kind=i4),intent(in)      :: mx   !! the integer(kind=i4) first dimension of q.  mx > 1 is required (see ier = 1)
    integer(kind=i4),intent(in)      :: my   !! the integer(kind=i4) second dimension of q. my > 1 is required (see ier = 1)
    integer(kind=i4),intent(in)      :: mz   !! the integer(kind=i4) third dimension of q. mz > 1 is required (see ier = 1)
    integer(kind=i4),intent(in)      :: mt   !! the integer(kind=i4) fourth dimension of q. mt > 1 is required (see ier = 1)
    integer(kind=i4),intent(in)      :: intpol(4)    !! an integer(kind=i4) vector of dimension 4 which sets linear or cubic
                                            !! interpolation in each of the x,y,z,t directions as follows:
                                            !!
                                            !! * intpol(1) = 1 sets linear interpolation in the x direction
                                            !! * intpol(1) = 3 sets cubic interpolation in the x direction.
                                            !! * intpol(2) = 1 sets linear interpolation in the y direction
                                            !! * intpol(2) = 3 sets cubic interpolation in the y direction.
                                            !! * intpol(3) = 1 sets linear interpolation in the z direction
                                            !! * intpol(3) = 3 sets cubic interpolation in the z direction.
                                            !! * intpol(4) = 1 sets linear interpolation in the t direction
                                            !! * intpol(4) = 3 sets cubic interpolation in the t direction.
                                            !!
                                            !! values other than 1 or 3 in intpol are not allowed (ier = 3).
    integer(kind=i4),intent(in)      :: liw          !! the integer(kind=i4) length of the integer(kind=i4) work space iw.
                                            !! liw must be at least mx+my+mz+mt
    integer(kind=i4),intent(inout)   :: iw(liw)      !! an integer(kind=i4) work space of length at least liw
                                            !! which must be provided in the routine calling rgrd4u
    integer(kind=i4),intent(in)      :: lw           !! the integer(kind=i4) length of the work space w.
                                            !!
                                            !! * let lwx = 1 if mx-1 divides nx-1; otherwise
                                            !!   let lwx = mx if intpol(1) = 1 or
                                            !!   let lwx = 4*mx if intpol(1) = 3
                                            !! * let lwy = 0 if my-1 divides ny-1; otherwise
                                            !!   let lwy = my+2*mx if intpol(2) = 1 or
                                            !!   let lwy = 4*(mx+my) if intpol(2) = 3
                                            !! * let lwz = 0 if mz-1 divides nz-1; otherwise
                                            !!   let lwz = 2*mx*my+mz if intpol(3) = 1 or
                                            !!   let lwz = 4*(mx*my+mz) if intpol(3) = 3
                                            !! * let lwt = 0 if mt-1 divides nt-1; otherwise
                                            !!   let lwt = 2*mx*my*mz+mt if intpol(4) = 1 or
                                            !!   let lwt = 4*(mx*my*mz+mt) if intpol(4) = 3
                                            !!
                                            !! then lw must be greater than or equal to lwx+lwy+lwz+lwt
    integer(kind=i4),intent(out)     :: ier          !! an integer(kind=i4) error flag set as follows:
                                            !!
                                            !! * ier = 0 if no errors in input arguments are detected
                                            !! * ier = 1 if min(mx,my,mz,mt) < 2
                                            !! * ier = 2 if nx < 2 when intpol(1)=1 or nx < 4 when intpol(1)=3 (or)
                                            !!                     ny < 2 when intpol(2)=1 or ny < 4 when intpol(2)=3 (or)
                                            !!                     nz < 2 when intpol(3)=1 or nz < 4 when intpol(3)=3 (or)
                                            !!                     nt < 2 when intpol(4)=1 or nt < 4 when intpol(4)=3.
                                            !! * ier = 3 if any of intpol(1),intpol(2),intpol(3),intpol(4)  is not
                                            !!                 equal to 1 or 3.
                                            !! * ier = 4 if lw or liw is too small (insufficient work space)
    real(kind=dp),intent(in)     :: p(nx,ny,nz,nt)   !! a real(kind=dp) nx by ny by nz by nt array of given values
    real(kind=dp),intent(out)    :: q(mx,my,mz,mt)   !! a real(kind=dp) mx by my by mz by mt array of values which are interpolated from p.
    real(kind=dp),intent(inout)  :: w(lw)    !! a real(kind=dp) work space of length at least lw
                                        !! which must be provided in the routine calling rgrd4u


    integer(kind=i4) :: inmx,jnmy,knmz,lnmt,isubx,jsuby,ksubz,lsubt
    integer(kind=i4) :: mxmy,mxmymz,lwx,lwy,lwz,lwt,jy,kz,lt
    integer(kind=i4) :: i2,i3,i4,i5
    integer(kind=i4) :: j2,j3,j4,j5,j6,j7,j8,j9
    integer(kind=i4) :: k2,k3,k4,k5,k6,k7,k8,k9
    integer(kind=i4) :: l2,l3,l4,l5,l6,l7,l8,l9

    ! check input arguments

    ! check mx,my,mz,mt
    ier = 1
    if (min(mx,my,mz,mt) < 1) return

    ! check intpol
    ier = 3
    if (intpol(1)/=1 .and. intpol(1)/=3) return
    if (intpol(2)/=1 .and. intpol(2)/=3) return
    if (intpol(3)/=1 .and. intpol(3)/=3) return
    if (intpol(4)/=1 .and. intpol(4)/=3) return

    ! check nx,ny,nz,nt
    ier = 2
    if (intpol(1)==1 .and. nx<2) return
    if (intpol(1)==3 .and. nx<4) return
    if (intpol(2)==1 .and. ny<2) return
    if (intpol(2)==3 .and. ny<4) return
    if (intpol(3)==1 .and. nz<2) return
    if (intpol(3)==3 .and. nz<4) return
    if (intpol(4)==1 .and. nt<2) return
    if (intpol(4)==3 .and. nt<4) return

    ! set subgrid indicators
    inmx = (nx-1)/(mx-1)
    jnmy = (ny-1)/(my-1)
    knmz = (nz-1)/(mz-1)
    lnmt = (nt-1)/(mt-1)
    isubx = nx - inmx*(mx-1)
    jsuby = ny - jnmy*(my-1)
    ksubz = nz - knmz*(mz-1)
    lsubt = nt - lnmt*(mt-1)

    ! check work space length input
    ier = 4
    mxmy = mx*my
    mxmymz = mxmy*mz
    lwx = 1
    if (isubx/=1) then
        if (intpol(1)==1) then
            lwx = mx
        else
            lwx = 4*mx
        end if
    end if
    lwy = 0
    if (jsuby/=1) then
        if (intpol(2)==1) then
            lwy = (2*mx+my)
        else
            lwy = 4*my+4*mx
        end if
    end if
    lwz = 0
    if (ksubz/=1) then
        if (intpol(3)==1) then
            lwz = (2*mxmy+mz)
        else
            lwz = 4*mxmy+4*mz
        end if
    end if
    lwt = 0
    if (lsubt/=1) then
        if (intpol(4)==1) then
            lwt = (2*mxmymz+mt)
        else
            lwt = 4*mxmymz+4*mt
        end if
    end if

    if (lw < lwx+lwy+lwz+lwt) return
    if (liw < mx+my+mz+mt) return

    ! arguments o.k.

    ier = 0
    jy = mx+1
    kz = mx+my+1
    lt = mx+my+mz+1

    if (intpol(4)==1) then

        ! linearly interpolate in nt, set work space pointers and scales

        l2 = 1
        l3 = l2
        l4 = l3+mt
        l5 = l4
        l6 = l5
        l7 = l6
        l8 = l7+mxmymz
        l9 = l8+mxmymz
        call linmxu(nt,mt,iw(lt),w(l3))
        k2 = l9
        if (intpol(3)==1) then
            ! linear in z
            k3 = k2
            k4 = k3+mz
            k5 = k4
            k6 = k5
            k7 = k6
            k8 = k7+mxmy
            k9 = k8+mxmy
            call linmxu(nz,mz,iw(kz),w(k3))
            j2 = k9
        else
            ! cubic in z
            k3 = k2+mz
            k4 = k3+mz
            k5 = k4+mz
            k6 = k5+mz
            k7 = k6+mxmy
            k8 = k7+mxmy
            k9 = k8+mxmy
            call cubnmxu(nz,mz,iw(kz),w(k2),w(k3),w(k4),w(k5))
            j2 = k9+mxmy
        end if

        if (intpol(2) == 1) then
            ! linear in y
            j3 = j2
            j4 = j3+my
            j5 = j4
            j6 = j5
            j7 = j6
            j8 = j7+mx
            j9 = j8+mx
            call linmxu(ny,my,iw(jy),w(j3))
            i2 = j9
        else
            ! cubic in y
            j3 = j2+my
            j4 = j3+my
            j5 = j4+my
            j6 = j5+my
            j7 = j6+mx
            j8 = j7+mx
            j9 = j8+mx
            call cubnmxu(ny,my,iw(jy),w(j2),w(j3),w(j4),w(j5))
            i2 = j9+mx
        end if

        if (intpol(1) == 1) then
            ! linear in x
            i3 = i2
            i4 = i3
            i5 = i4
            call linmxu(nx,mx,iw,w(i3))
        else
            ! cubic in x
            i3 = i2+mx
            i4 = i3+mx
            i5 = i4+mx
            call cubnmxu(nx,mx,iw,w(i2),w(i3),w(i4),w(i5))
        end if

        ! linearly interpolate in t

        call lint4u(nx,ny,nz,nt,p,mx,my,mz,mt,mxmy,mxmymz,q,intpol,&
                    iw(lt),w(l3),w(l7),w(l8),&
                    iw(kz),w(k2),w(k3),w(k4),w(k5),w(k6),w(k7),w(k8),w(k9),&
                    iw(jy),w(j2),w(j3),w(j4),w(j5),w(j6),w(j7),w(j8),w(j9),&
                    iw,w(i2),w(i3),w(i4),w(i5),&
                    inmx,jnmy,knmz,lnmt,isubx,jsuby,ksubz,lsubt)

    else

        ! cubically interpolate in t

        l2 = 1
        l3 = l2+mt
        l4 = l3+mt
        l5 = l4+mt
        l6 = l5+mt
        l7 = l6+mxmymz
        l8 = l7+mxmymz
        l9 = l8+mxmymz
        call cubnmxu(nt,mt,iw(lt),w(l2),w(l3),w(l4),w(l5))
        k2 = l9+mxmymz

        if (intpol(3)==1) then
            ! linear in z
            k3 = k2
            k4 = k3+mz
            k5 = k4
            k6 = k5
            k7 = k6
            k8 = k7+mxmy
            k9 = k8+mxmy
            call linmxu(nz,mz,iw(kz),w(k3))
            j2 = k9
        else
            ! cubic in z
            k3 = k2+mz
            k4 = k3+mz
            k5 = k4+mz
            k6 = k5+mz
            k7 = k6+mxmy
            k8 = k7+mxmy
            k9 = k8+mxmy
            call cubnmxu(nz,mz,iw(kz),w(k2),w(k3),w(k4),w(k5))
            j2 = k9+mxmy
        end if

        if (intpol(2) == 1) then
            j3 = j2
            j4 = j3+my
            j5 = j4
            j6 = j5
            j7 = j6
            j8 = j7+mx
            j9 = j8+mx
            call linmxu(ny,my,iw(jy),w(j3))
            i2 = j9
        else
            j3 = j2+my
            j4 = j3+my
            j5 = j4+my
            j6 = j5+my
            j7 = j6+mx
            j8 = j7+mx
            j9 = j8+mx
            call cubnmxu(ny,my,iw(jy),w(j2),w(j3),w(j4),w(j5))
            i2 = j9+mx
        end if

        ! set work space portion and indices which depend on x interpolation

        if (intpol(1) == 1) then
            i3 = i2
            i4 = i3
            i5 = i4
            call linmxu(nx,mx,iw,w(i3))
        else
            i3 = i2+mx
            i4 = i3+mx
            i5 = i4+mx
            call cubnmxu(nx,mx,iw,w(i2),w(i3),w(i4),w(i5))
        end if

        ! cubically interpolate in t

        call cubt4u(nx,ny,nz,nt,p,mx,my,mz,mt,mxmy,mxmymz,q,intpol,&
                    iw(lt),w(l2),w(l3),w(l4),w(l5),w(l6),w(l7),w(l8),w(l9),&
                    iw(kz),w(k2),w(k3),w(k4),w(k5),w(k6),w(k7),w(k8),w(k9),&
                    iw(jy),w(j2),w(j3),w(j4),w(j5),w(j6),w(j7),w(j8),w(j9),&
                    iw,w(i2),w(i3),w(i4),w(i5),&
                    inmx,jnmy,knmz,lnmt,isubx,jsuby,ksubz,lsubt)

    end if

    end subroutine rgrd4u
#endif
!**************************************************************************

!**************************************************************************
!>
!  linearly interpolate in t direction
#if 0
    subroutine lint4u(nx,ny,nz,nt,p,mx,my,mz,mt,mxmy,mxmymz,q,intpol,&
                      lt,dt,pt,ptp,kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,&
                      jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                      inmx,jnmy,knmz,lnmt,isubx,jsuby,ksubz,lsubt)

    implicit none

    integer(kind=i4)  :: nx,ny,nz,nt,mx,my,mz,mt
    integer(kind=i4)  :: mxmy,mxmymz
    real(kind=dp) :: p(nx,ny,nz,nt),q(mxmymz,mt)
    integer(kind=i4)  :: inmx,jnmy,knmz,lnmt,isubx,jsuby,ksubz,lsubt
    real(kind=dp) :: dt(mt),pt(mxmymz),ptp(mxmymz)
    real(kind=dp) :: dzm(mz),dz(mz),dzp(mz),dzpp(mz)
    real(kind=dp) :: pkm(mxmy),pk(mxmy),pkp(mxmy),pkpp(mxmy)
    real(kind=dp) :: dym(my),dy(my),dyp(my),dypp(my)
    real(kind=dp) :: pjm(mx),pj(mx),pjp(mx),pjpp(mx)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)
    integer(kind=i4)  :: lt(mt),kz(mz),jy(my),ix(mx),intpol(4)
    integer(kind=i4)  :: l,ll,lsave,iijjkk

    if (intpol(3) == 1) then
        ! linear in z
        if (lsubt == 1) then
            ! mt grid is subset of nt grid
            do ll=1,mt
                l = lnmt*(ll-1)+1
                call lint3u(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,q(1,ll),intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
            end do
            return
        end if

        lsave = -1
        do ll=1,mt
            l = lt(ll)
            if (l==lsave) then
                ! l pointer has not moved since last pass (no updates or interpolation)
            else if (l==lsave+1) then
                ! update l and interpolate l+1
                do iijjkk=1,mxmymz
                    pt(iijjkk) = ptp(iijjkk)
                end do
                call lint3u(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
            else
                ! interpolate l,l+1 in pt,ptp on xx,yy,zz mesh
                call lint3u(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call lint3u(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
            end if

            ! save l pointer for next pass

            lsave = l

            ! linearly interpolate q(ii,jj,,kk,ll) from pt,ptp in t direction

            do iijjkk=1,mxmymz
                q(iijjkk,ll) = pt(iijjkk)+dt(ll)*(ptp(iijjkk)-pt(iijjkk))
            end do
        end do

    else

        ! cubic in z

        if (lsubt == 1) then
            ! mt grid is subset of nt grid
            do ll=1,mt
                l = lnmt*(ll-1)+1
                call cubt3u(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,q(1,ll),intpol,&
                    kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,&
                    pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
            end do
            return
        end if

        lsave = -1
        do ll=1,mt
            l = lt(ll)
            if (l==lsave) then
                ! l pointer has not moved since last pass (no updates or interpolation)
            else if (l==lsave+1) then
                ! update l and interpolate l+1
                do iijjkk=1,mxmymz
                    pt(iijjkk) = ptp(iijjkk)
                end do
                call cubt3u(nx,ny,nt,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,&
                    kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,&
                    jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
            else
                ! interpolate l,l+1 in pt,ptp on xx,yy,zz mesh
                call cubt3u(nx,ny,nt,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,&
                    kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,&
                    jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call cubt3u(nx,ny,nt,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,&
                    kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,&
                    jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
            end if

            ! save l pointer for next pass
            lsave = l

            ! linearly interpolate q(ii,jj,kk,ll) from pt,ptp in t direction
            do iijjkk=1,mxmymz
                q(iijjkk,ll) = pt(iijjkk)+dt(ll)*(ptp(iijjkk)-pt(iijjkk))
            end do

        end do

    end if

    end subroutine lint4u
#endif
!**************************************************************************

!**************************************************************************
!>
!  cubically interpolate in t
#if 0
    subroutine cubt4u(nx,ny,nz,nt,p,mx,my,mz,mt,mxmy,mxmymz,q,intpol,&
                        lt,dtm,dt,dtp,dtpp,ptm,pt,ptp,ptpp,&
                        kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,&
                        jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,&
                        ix,dxm,dx,dxp,dxpp,&
                        inmx,jnmy,knmz,lnmt,isubx,jsuby,ksubz,lsubt)

    implicit none

    integer(kind=i4)  :: nx,ny,nz,nt,mx,my,mz,mt
    integer(kind=i4)  :: inmx,jnmy,knmz,lnmt,isubx,jsuby,ksubz,lsubt
    integer(kind=i4)  :: mxmy,mxmymz
    real(kind=dp) :: p(nx,ny,nz,nt),q(mxmymz,mt)
    real(kind=dp) :: ptm(mxmymz),pt(mxmymz),ptp(mxmymz),ptpp(mxmymz)
    real(kind=dp) :: dtm(mt),dt(mt),dtp(mt),dtpp(mt)
    real(kind=dp) :: dzm(mz),dz(mz),dzp(mz),dzpp(mz)
    real(kind=dp) :: pkm(mxmy),pk(mxmy),pkp(mxmy),pkpp(mxmy)
    real(kind=dp) :: dym(my),dy(my),dyp(my),dypp(my)
    real(kind=dp) :: pjm(mx),pj(mx),pjp(mx),pjpp(mx)
    real(kind=dp) :: dxm(mx),dx(mx),dxp(mx),dxpp(mx)
    integer(kind=i4)  :: lt(mt),kz(mz),jy(my),ix(mx),intpol(4)
    integer(kind=i4)  :: l,ll,iijjkk,lsave

    if (intpol(3) == 1) then

        ! linear in z

        if (lsubt == 1) then
            ! mt grid is subset of nt grid
            do ll=1,mt
                l = lnmt*(ll-1)+1
                call lint3u(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,q(1,ll),intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
            end do
            return
        end if
        lsave = -3
        do ll=1,mt
            l = lt(ll)
            if (l==lsave) then
                ! l pointer has not moved since last pass (no updates or interpolation)
            else if (l==lsave+1) then
                ! update l-1,l,l+1 and interpolate l+2
                do iijjkk=1,mxmymz
                    ptm(iijjkk) = pt(iijjkk)
                    pt(iijjkk) = ptp(iijjkk)
                    ptp(iijjkk) = ptpp(iijjkk)
                end do
                call lint3u(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
            else if (l==lsave+2) then
                ! update l-1,l and interpolate l+1,l+2
                do iijjkk=1,mxmymz
                    ptm(iijjkk) = ptp(iijjkk)
                    pt(iijjkk) = ptpp(iijjkk)
                end do
                call lint3u(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call lint3u(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
            else if (l==lsave+3) then
                ! update l-1 and interpolate l,l+1,l+2
                do iijjkk=1,mxmymz
                    ptm(iijjkk) = ptpp(iijjkk)
                end do
                call lint3u(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,kz,dz,&
                            pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call lint3u(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call lint3u(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
            else
                ! interpolate all four l-1,l,l+1,l+2
                call lint3u(nx,ny,nz,p(1,1,1,l-1),mx,my,mxmy,mz,ptm,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call lint3u(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call lint3u(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call lint3u(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,&
                            dz,pk,pkp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                            inmx,jnmy,knmz,isubx,jsuby,ksubz)
            end if

            ! save l pointer for next pass

            lsave = l

            ! cubically interpolate q(ii,jj,kk,ll) from ptm,pt,ptp,ptpp in t direction

            do iijjkk=1,mxmymz
                q(iijjkk,ll) = dtm(ll)*ptm(iijjkk) + dt(ll)*pt(iijjkk) + dtp(ll)*ptp(iijjkk) + dtpp(ll)*ptpp(iijjkk)
            end do
        end do

    else

        ! cubic in z

        if (lsubt == 1) then

            ! mt grid is subset of nt grid

            do ll=1,mt
                l = lnmt*(ll-1)+1
                call cubt3u(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,q(1,ll),intpol,&
                    kz,dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,&
                    pjm,pj,pjp,pjpp,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
            end do
            return
        end if
        lsave = -3
        do ll=1,mt
            l = lt(ll)
            if (l==lsave) then

                ! l pointer has not moved since last pass (no updates or interpolation)

            else if (l==lsave+1) then

                ! update l-1,l,l+1 and interpolate l+2

                do iijjkk=1,mxmymz
                    ptm(iijjkk) = pt(iijjkk)
                    pt(iijjkk) = ptp(iijjkk)
                    ptp(iijjkk) = ptpp(iijjkk)
                end do
                call cubt3u(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,&
                    dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp&
                    ,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
            else if (l==lsave+2) then

                ! update l-1,l and interpolate l+1,l+2

                do iijjkk=1,mxmymz
                    ptm(iijjkk) = ptp(iijjkk)
                    pt(iijjkk) = ptpp(iijjkk)
                end do
                call cubt3u(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,&
                    dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp&
                    ,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call cubt3u(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,&
                    dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp&
                    ,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
            else if (l==lsave+3) then

                ! update l-1 and interpolate l,l+1,l+2

                do iijjkk=1,mxmymz
                    ptm(iijjkk) = ptpp(iijjkk)
                end do
                call cubt3u(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,kz,&
                    dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp&
                    ,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call cubt3u(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,&
                    dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp&
                    ,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call cubt3u(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,&
                    dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp&
                    ,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
            else

                ! interpolate all four l-1,l,l+1,l+2

                call cubt3u(nx,ny,nz,p(1,1,1,l-1),mx,my,mxmy,mz,ptm,intpol,kz,&
                    dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp&
                    ,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call cubt3u(nx,ny,nz,p(1,1,1,l),mx,my,mxmy,mz,pt,intpol,kz,&
                    dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp&
                    ,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call cubt3u(nx,ny,nz,p(1,1,1,l+1),mx,my,mxmy,mz,ptp,intpol,kz,&
                    dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp&
                    ,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
                call cubt3u(nx,ny,nz,p(1,1,1,l+2),mx,my,mxmy,mz,ptpp,intpol,kz,&
                    dzm,dz,dzp,dzpp,pkm,pk,pkp,pkpp,jy,dym,dy,dyp,dypp,pjm,pj,pjp,pjpp&
                    ,ix,dxm,dx,dxp,dxpp,&
                    inmx,jnmy,knmz,isubx,jsuby,ksubz)
            end if

            ! save l pointer for next pass
            lsave = l

            ! cubically interpolate q(ii,jj,kk,ll) from ptm,pt,ptp,ptpp in t direction
            do iijjkk=1,mxmymz
                q(iijjkk,ll) = dtm(ll)*ptm(iijjkk) + dt(ll)*pt(iijjkk) + dtp(ll)*ptp(iijjkk) + dtpp(ll)*ptpp(iijjkk)
            end do
        end do

    end if

    end subroutine cubt4u
#endif
!**************************************************************************

!***************************************************************************************************
    end module regridpack
!***************************************************************************************************
