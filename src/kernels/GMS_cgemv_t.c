



//=========================================================================
// Modified and optimized version of OpenBLAS zscal and kernels zscal_kernel_x
// Programmer: Bernard Gingold, contact: beniekg@gmail.com
// 20-02-2021 14:58  +00200
// Original copyright below
//========================================================================
/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include <cstdint>
#include "GMS_cgemv_t.h" 

 // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)

void cgemv_kernel_4x4(const int32_t n,
		      float ** __restrict ap,
		      float * __restrict  x,
		      float * __restrict  y,
		      float * __restrict  alpha) {
      int32_t i = 0;
      bool isAligned32 = false;
      isAligned32 = ((uintptr_t)x & 0x1F) != 0ULL && ((uintptr_t)y & 0x1F) != 0ULL;
      if(!isAligned32) {

          __asm__ __volatile__ (

           "vzeroupper			 \n\t"

	"vxorps		%%ymm8 , %%ymm8 , %%ymm8 	\n\t" // temp
	"vxorps		%%ymm9 , %%ymm9 , %%ymm9 	\n\t" // temp
	"vxorps		%%ymm10, %%ymm10, %%ymm10	\n\t" // temp
	"vxorps		%%ymm11, %%ymm11, %%ymm11	\n\t" // temp
	"vxorps		%%ymm12, %%ymm12, %%ymm12	\n\t" // temp
	"vxorps		%%ymm13, %%ymm13, %%ymm13	\n\t"
	"vxorps		%%ymm14, %%ymm14, %%ymm14	\n\t"
	"vxorps		%%ymm15, %%ymm15, %%ymm15	\n\t"

        "testq          $0x04, %1                      \n\t"
        "jz             2f                      \n\t"

	"vmovups	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
	"vmovups	(%5,%0,4), %%ymm5               \n\t" // 4 complex values from a1

	"vmovups	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	
	"vmovups	(%6,%0,4), %%ymm6	        \n\t" // 4 complex values from a2
	"vmovups	(%7,%0,4), %%ymm7               \n\t" // 4 complex values from a3

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm6 , %%ymm0, %%ymm12      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm6 , %%ymm1, %%ymm13      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm7 , %%ymm0, %%ymm14      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm7 , %%ymm1, %%ymm15      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	
        "addq		$8  , %0	  	 	        \n\t"
	"subq	        $4  , %1			        \n\t"		

        "2:                                  \n\t"
	"cmpq           $0, %1                         \n\t"
        "je             3f                      \n\t"

		".align 16				        \n\t"
	"1:				        \n\t"
        "prefetcht0      256(%4,%0,4)                   \n\t"
	"vmovups	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
        "prefetcht0      256(%5,%0,4)                   \n\t"
	"vmovups	(%5,%0,4), %%ymm5               \n\t" // 4 complex values from a1

        "prefetcht0      256(%2,%0,4)                   \n\t"
	"vmovups	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	
        "prefetcht0      256(%6,%0,4)                   \n\t"
	"vmovups	(%6,%0,4), %%ymm6	        \n\t" // 4 complex values from a2
        "prefetcht0      256(%7,%0,4)                   \n\t"
	"vmovups	(%7,%0,4), %%ymm7               \n\t" // 4 complex values from a3

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm6 , %%ymm0, %%ymm12      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm6 , %%ymm1, %%ymm13      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm7 , %%ymm0, %%ymm14      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm7 , %%ymm1, %%ymm15      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

	"vmovups       32(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
	"vmovups       32(%5,%0,4), %%ymm5              \n\t" // 4 complex values from a1

	"vmovups	  32(%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts

	"vmovups       32(%6,%0,4), %%ymm6	        \n\t" // 4 complex values from a2
	"vmovups       32(%7,%0,4), %%ymm7              \n\t" // 4 complex values from a3

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm6 , %%ymm0, %%ymm12      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm6 , %%ymm1, %%ymm13      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm7 , %%ymm0, %%ymm14      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm7 , %%ymm1, %%ymm15      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

        "addq		$16 , %0	  	 	        \n\t"
	"subq	        $8  , %1			        \n\t"		
	"jnz		1b		        \n\t"

        "3:                                   \n\t"

        "vbroadcastss    (%8)  , %%xmm0                \n\t"  // value from alpha
        "vbroadcastss   4(%8)  , %%xmm1                \n\t"  // value from alpha


#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vpermilps      $0xb1 , %%ymm9 , %%ymm9                \n\t"
        "vpermilps      $0xb1 , %%ymm11, %%ymm11               \n\t"
        "vpermilps      $0xb1 , %%ymm13, %%ymm13               \n\t"
        "vpermilps      $0xb1 , %%ymm15, %%ymm15               \n\t"
        "vaddsubps      %%ymm9 , %%ymm8, %%ymm8                \n\t" 
        "vaddsubps      %%ymm11, %%ymm10, %%ymm10              \n\t"
        "vaddsubps      %%ymm13, %%ymm12, %%ymm12              \n\t"
        "vaddsubps      %%ymm15, %%ymm14, %%ymm14              \n\t"
#else
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
        "vpermilps      $0xb1 , %%ymm10, %%ymm10               \n\t"
        "vpermilps      $0xb1 , %%ymm12, %%ymm12               \n\t"
        "vpermilps      $0xb1 , %%ymm14, %%ymm14               \n\t"
        "vaddsubps      %%ymm8 , %%ymm9 , %%ymm8               \n\t"
        "vaddsubps      %%ymm10, %%ymm11, %%ymm10              \n\t"
        "vaddsubps      %%ymm12, %%ymm13, %%ymm12              \n\t"
        "vaddsubps      %%ymm14, %%ymm15, %%ymm14              \n\t"
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
        "vpermilps      $0xb1 , %%ymm10, %%ymm10               \n\t"
        "vpermilps      $0xb1 , %%ymm12, %%ymm12               \n\t"
        "vpermilps      $0xb1 , %%ymm14, %%ymm14               \n\t"
#endif

	"vmovsd         (%3), %%xmm4			\n\t" // read y
	"vmovsd        8(%3), %%xmm5			\n\t"
	"vmovsd       16(%3), %%xmm6			\n\t"
	"vmovsd       24(%3), %%xmm7			\n\t"

	"vextractf128   $1, %%ymm8 , %%xmm9		      \n\t"
	"vextractf128   $1, %%ymm10, %%xmm11	      	      \n\t"
	"vextractf128   $1, %%ymm12, %%xmm13		      \n\t"
	"vextractf128   $1, %%ymm14, %%xmm15		      \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm11, %%xmm10      \n\t"
	"vaddps		%%xmm12, %%xmm13, %%xmm12      \n\t"
	"vaddps		%%xmm14, %%xmm15, %%xmm14      \n\t"

	"vshufpd        $0x1, %%xmm8 , %%xmm8 , %%xmm9   \n\t"
	"vshufpd        $0x1, %%xmm10, %%xmm10, %%xmm11  \n\t"
	"vshufpd        $0x1, %%xmm12, %%xmm12, %%xmm13  \n\t"
	"vshufpd        $0x1, %%xmm14, %%xmm14, %%xmm15  \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm11, %%xmm10      \n\t"
	"vaddps		%%xmm12, %%xmm13, %%xmm12      \n\t"
	"vaddps		%%xmm14, %%xmm15, %%xmm14      \n\t"


        "vmulps         %%xmm8 , %%xmm1 , %%xmm9              \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm8 , %%xmm0 , %%xmm8              \n\t"  // t_r * alpha_r , t_i * alpha_r
        "vmulps         %%xmm10, %%xmm1 , %%xmm11             \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm10, %%xmm0 , %%xmm10             \n\t"  // t_r * alpha_r , t_i * alpha_r
        "vmulps         %%xmm12, %%xmm1 , %%xmm13             \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm12, %%xmm0 , %%xmm12             \n\t"  // t_r * alpha_r , t_i * alpha_r
        "vmulps         %%xmm14, %%xmm1 , %%xmm15             \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm14, %%xmm0 , %%xmm14             \n\t"  // t_r * alpha_r , t_i * alpha_r

#if !defined(XCONJ)
        "vpermilps      $0xb1 , %%xmm9 , %%xmm9                \n\t"
        "vpermilps      $0xb1 , %%xmm11, %%xmm11               \n\t"
        "vpermilps      $0xb1 , %%xmm13, %%xmm13               \n\t"
        "vpermilps      $0xb1 , %%xmm15, %%xmm15               \n\t"
        "vaddsubps      %%xmm9 , %%xmm8, %%xmm8               \n\t"
        "vaddsubps      %%xmm11, %%xmm10, %%xmm10             \n\t"
        "vaddsubps      %%xmm13, %%xmm12, %%xmm12             \n\t"
        "vaddsubps      %%xmm15, %%xmm14, %%xmm14             \n\t"
#else
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
        "vpermilps      $0xb1 , %%xmm10, %%xmm10               \n\t"
        "vpermilps      $0xb1 , %%xmm12, %%xmm12               \n\t"
        "vpermilps      $0xb1 , %%xmm14, %%xmm14               \n\t"
        "vaddsubps      %%xmm8 , %%xmm9 , %%xmm8              \n\t"
        "vaddsubps      %%xmm10, %%xmm11, %%xmm10             \n\t"
        "vaddsubps      %%xmm12, %%xmm13, %%xmm12             \n\t"
        "vaddsubps      %%xmm14, %%xmm15, %%xmm14             \n\t"
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
        "vpermilps      $0xb1 , %%xmm10, %%xmm10               \n\t"
        "vpermilps      $0xb1 , %%xmm12, %%xmm12               \n\t"
        "vpermilps      $0xb1 , %%xmm14, %%xmm14               \n\t"
#endif


	"vaddps		%%xmm8 , %%xmm4 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm5 , %%xmm10      \n\t"
	"vaddps		%%xmm12, %%xmm6 , %%xmm12      \n\t"
	"vaddps		%%xmm14, %%xmm7 , %%xmm14      \n\t"

	"vmovsd	%%xmm8 ,   (%3)			\n\t"
	"vmovsd	%%xmm10,  8(%3)			\n\t"
	"vmovsd	%%xmm12, 16(%3)			\n\t"
	"vmovsd	%%xmm14, 24(%3)			\n\t"

	"vzeroupper			 \n\t"

	:
          "+r" (i),	// 0	
	  "+r" (n)  	// 1
	:
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (ap[0]),  // 4
          "r" (ap[1]),  // 5
          "r" (ap[2]),  // 6
          "r" (ap[3]),  // 7
          "r" (alpha)   // 8
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	  );
        
      } else {

              __asm__ __volatile__ (

           "vzeroupper			 \n\t"

	"vxorps		%%ymm8 , %%ymm8 , %%ymm8 	\n\t" // temp
	"vxorps		%%ymm9 , %%ymm9 , %%ymm9 	\n\t" // temp
	"vxorps		%%ymm10, %%ymm10, %%ymm10	\n\t" // temp
	"vxorps		%%ymm11, %%ymm11, %%ymm11	\n\t" // temp
	"vxorps		%%ymm12, %%ymm12, %%ymm12	\n\t" // temp
	"vxorps		%%ymm13, %%ymm13, %%ymm13	\n\t"
	"vxorps		%%ymm14, %%ymm14, %%ymm14	\n\t"
	"vxorps		%%ymm15, %%ymm15, %%ymm15	\n\t"

        "testq          $0x04, %1                      \n\t"
        "jz             2f                      \n\t"

	"vmovaps	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
	"vmovaps	(%5,%0,4), %%ymm5               \n\t" // 4 complex values from a1

	"vmovups	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	
	"vmovaps	(%6,%0,4), %%ymm6	        \n\t" // 4 complex values from a2
	"vmovaps	(%7,%0,4), %%ymm7               \n\t" // 4 complex values from a3

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm6 , %%ymm0, %%ymm12      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm6 , %%ymm1, %%ymm13      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm7 , %%ymm0, %%ymm14      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm7 , %%ymm1, %%ymm15      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	
        "addq		$8  , %0	  	 	        \n\t"
	"subq	        $4  , %1			        \n\t"		

        "2:                                  \n\t"
	"cmpq           $0, %1                         \n\t"
        "je             3f                      \n\t"

		".align 16				        \n\t"
	"1:				        \n\t"
        "prefetcht0      256(%4,%0,4)                   \n\t"
	"vmovaps	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
        "prefetcht0      256(%5,%0,4)                   \n\t"
	"vmovaps	(%5,%0,4), %%ymm5               \n\t" // 4 complex values from a1

        "prefetcht0      256(%2,%0,4)                   \n\t"
	"vmovaps	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	
        "prefetcht0      256(%6,%0,4)                   \n\t"
	"vmovaps	(%6,%0,4), %%ymm6	        \n\t" // 4 complex values from a2
        "prefetcht0      256(%7,%0,4)                   \n\t"
	"vmovaps	(%7,%0,4), %%ymm7               \n\t" // 4 complex values from a3

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm6 , %%ymm0, %%ymm12      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm6 , %%ymm1, %%ymm13      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm7 , %%ymm0, %%ymm14      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm7 , %%ymm1, %%ymm15      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

	"vmovaps       32(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
	"vmovaps       32(%5,%0,4), %%ymm5              \n\t" // 4 complex values from a1

	"vmovaps	  32(%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts

	"vmovaps       32(%6,%0,4), %%ymm6	        \n\t" // 4 complex values from a2
	"vmovaps       32(%7,%0,4), %%ymm7              \n\t" // 4 complex values from a3

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm6 , %%ymm0, %%ymm12      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm6 , %%ymm1, %%ymm13      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm7 , %%ymm0, %%ymm14      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm7 , %%ymm1, %%ymm15      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

        "addq		$16 , %0	  	 	        \n\t"
	"subq	        $8  , %1			        \n\t"		
	"jnz		1b		        \n\t"

        "3:                                   \n\t"

        "vbroadcastss    (%8)  , %%xmm0                \n\t"  // value from alpha
        "vbroadcastss   4(%8)  , %%xmm1                \n\t"  // value from alpha


#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vpermilps      $0xb1 , %%ymm9 , %%ymm9                \n\t"
        "vpermilps      $0xb1 , %%ymm11, %%ymm11               \n\t"
        "vpermilps      $0xb1 , %%ymm13, %%ymm13               \n\t"
        "vpermilps      $0xb1 , %%ymm15, %%ymm15               \n\t"
        "vaddsubps      %%ymm9 , %%ymm8, %%ymm8                \n\t" 
        "vaddsubps      %%ymm11, %%ymm10, %%ymm10              \n\t"
        "vaddsubps      %%ymm13, %%ymm12, %%ymm12              \n\t"
        "vaddsubps      %%ymm15, %%ymm14, %%ymm14              \n\t"
#else
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
        "vpermilps      $0xb1 , %%ymm10, %%ymm10               \n\t"
        "vpermilps      $0xb1 , %%ymm12, %%ymm12               \n\t"
        "vpermilps      $0xb1 , %%ymm14, %%ymm14               \n\t"
        "vaddsubps      %%ymm8 , %%ymm9 , %%ymm8               \n\t"
        "vaddsubps      %%ymm10, %%ymm11, %%ymm10              \n\t"
        "vaddsubps      %%ymm12, %%ymm13, %%ymm12              \n\t"
        "vaddsubps      %%ymm14, %%ymm15, %%ymm14              \n\t"
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
        "vpermilps      $0xb1 , %%ymm10, %%ymm10               \n\t"
        "vpermilps      $0xb1 , %%ymm12, %%ymm12               \n\t"
        "vpermilps      $0xb1 , %%ymm14, %%ymm14               \n\t"
#endif

	"vmovsd         (%3), %%xmm4			\n\t" // read y
	"vmovsd        8(%3), %%xmm5			\n\t"
	"vmovsd       16(%3), %%xmm6			\n\t"
	"vmovsd       24(%3), %%xmm7			\n\t"

	"vextractf128   $1, %%ymm8 , %%xmm9		      \n\t"
	"vextractf128   $1, %%ymm10, %%xmm11	      	      \n\t"
	"vextractf128   $1, %%ymm12, %%xmm13		      \n\t"
	"vextractf128   $1, %%ymm14, %%xmm15		      \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm11, %%xmm10      \n\t"
	"vaddps		%%xmm12, %%xmm13, %%xmm12      \n\t"
	"vaddps		%%xmm14, %%xmm15, %%xmm14      \n\t"

	"vshufpd        $0x1, %%xmm8 , %%xmm8 , %%xmm9   \n\t"
	"vshufpd        $0x1, %%xmm10, %%xmm10, %%xmm11  \n\t"
	"vshufpd        $0x1, %%xmm12, %%xmm12, %%xmm13  \n\t"
	"vshufpd        $0x1, %%xmm14, %%xmm14, %%xmm15  \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm11, %%xmm10      \n\t"
	"vaddps		%%xmm12, %%xmm13, %%xmm12      \n\t"
	"vaddps		%%xmm14, %%xmm15, %%xmm14      \n\t"


        "vmulps         %%xmm8 , %%xmm1 , %%xmm9              \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm8 , %%xmm0 , %%xmm8              \n\t"  // t_r * alpha_r , t_i * alpha_r
        "vmulps         %%xmm10, %%xmm1 , %%xmm11             \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm10, %%xmm0 , %%xmm10             \n\t"  // t_r * alpha_r , t_i * alpha_r
        "vmulps         %%xmm12, %%xmm1 , %%xmm13             \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm12, %%xmm0 , %%xmm12             \n\t"  // t_r * alpha_r , t_i * alpha_r
        "vmulps         %%xmm14, %%xmm1 , %%xmm15             \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm14, %%xmm0 , %%xmm14             \n\t"  // t_r * alpha_r , t_i * alpha_r

#if !defined(XCONJ)
        "vpermilps      $0xb1 , %%xmm9 , %%xmm9                \n\t"
        "vpermilps      $0xb1 , %%xmm11, %%xmm11               \n\t"
        "vpermilps      $0xb1 , %%xmm13, %%xmm13               \n\t"
        "vpermilps      $0xb1 , %%xmm15, %%xmm15               \n\t"
        "vaddsubps      %%xmm9 , %%xmm8, %%xmm8               \n\t"
        "vaddsubps      %%xmm11, %%xmm10, %%xmm10             \n\t"
        "vaddsubps      %%xmm13, %%xmm12, %%xmm12             \n\t"
        "vaddsubps      %%xmm15, %%xmm14, %%xmm14             \n\t"
#else
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
        "vpermilps      $0xb1 , %%xmm10, %%xmm10               \n\t"
        "vpermilps      $0xb1 , %%xmm12, %%xmm12               \n\t"
        "vpermilps      $0xb1 , %%xmm14, %%xmm14               \n\t"
        "vaddsubps      %%xmm8 , %%xmm9 , %%xmm8              \n\t"
        "vaddsubps      %%xmm10, %%xmm11, %%xmm10             \n\t"
        "vaddsubps      %%xmm12, %%xmm13, %%xmm12             \n\t"
        "vaddsubps      %%xmm14, %%xmm15, %%xmm14             \n\t"
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
        "vpermilps      $0xb1 , %%xmm10, %%xmm10               \n\t"
        "vpermilps      $0xb1 , %%xmm12, %%xmm12               \n\t"
        "vpermilps      $0xb1 , %%xmm14, %%xmm14               \n\t"
#endif


	"vaddps		%%xmm8 , %%xmm4 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm5 , %%xmm10      \n\t"
	"vaddps		%%xmm12, %%xmm6 , %%xmm12      \n\t"
	"vaddps		%%xmm14, %%xmm7 , %%xmm14      \n\t"

	"vmovsd	%%xmm8 ,   (%3)			\n\t"
	"vmovsd	%%xmm10,  8(%3)			\n\t"
	"vmovsd	%%xmm12, 16(%3)			\n\t"
	"vmovsd	%%xmm14, 24(%3)			\n\t"

	"vzeroupper			 \n\t"

	:
          "+r" (i),	// 0	
	  "+r" (n)  	// 1
	:
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (ap[0]),  // 4
          "r" (ap[1]),  // 5
          "r" (ap[2]),  // 6
          "r" (ap[3]),  // 7
          "r" (alpha)   // 8
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	  );
      }
}

 // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)

void cgemv_kernel_4x2(const int32_t n,
                      float ** __restrict ap,
		      float *  __restrict x,
		      float *  __restrict y,
		      float *  __restrict alpha) {

       int32_t i = 0;
       bool isAligned32 = false;
       isAligned32 = ((uintptr_t)x & 0x1F) != 0ULL && ((uintptr_t)y & 0x1F) != 0ULL;
       if(!isAligned32) {
        
            __asm__ __volatile__ (

            "vzeroupper			 \n\t"

	"vxorps		%%ymm8 , %%ymm8 , %%ymm8 	\n\t" // temp
	"vxorps		%%ymm9 , %%ymm9 , %%ymm9 	\n\t" // temp
	"vxorps		%%ymm10, %%ymm10, %%ymm10	\n\t" // temp
	"vxorps		%%ymm11, %%ymm11, %%ymm11	\n\t" // temp

        "testq          $0x04, %1                      \n\t"
        "jz             2f                    \n\t"

	"vmovups	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
	"vmovups	(%5,%0,4), %%ymm5               \n\t" // 4 complex values from a1

	"vmovups	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	
        "addq		$8  , %0	  	 	        \n\t"
	"subq	        $4  , %1			        \n\t"		

        "2:                                  \n\t"
	"cmpq           $0, %1                         \n\t"
        "je             3f                      \n\t"

		".align 16				        \n\t"
	"1:				        \n\t"
        "prefetcht0      256(%4,%0,4)                   \n\t"
	"vmovups	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
        "prefetcht0      256(%5,%0,4)                   \n\t"
	"vmovups	(%5,%0,4), %%ymm5               \n\t" // 4 complex values from a1

        "prefetcht0      256(%2,%0,4)                   \n\t"
	"vmovups	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	
	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

	"vmovups       32(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
	"vmovups       32(%5,%0,4), %%ymm5              \n\t" // 4 complex values from a1

	"vmovups	  32(%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

        "addq		$16 , %0	  	 	        \n\t"
	"subq	        $8  , %1			        \n\t"		
	"jnz		1b		        \n\t"

        "3:                                   \n\t"

        "vbroadcastss    (%6)  , %%xmm0                \n\t"  // value from alpha
        "vbroadcastss   4(%6)  , %%xmm1                \n\t"  // value from alpha


#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vpermilps      $0xb1 , %%ymm9 , %%ymm9                \n\t"
        "vpermilps      $0xb1 , %%ymm11, %%ymm11               \n\t"
        "vaddsubps      %%ymm9 , %%ymm8, %%ymm8                \n\t" 
        "vaddsubps      %%ymm11, %%ymm10, %%ymm10              \n\t"
#else
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
        "vpermilps      $0xb1 , %%ymm10, %%ymm10               \n\t"
        "vaddsubps      %%ymm8 , %%ymm9 , %%ymm8               \n\t"
        "vaddsubps      %%ymm10, %%ymm11, %%ymm10              \n\t"
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
        "vpermilps      $0xb1 , %%ymm10, %%ymm10               \n\t"
#endif

	"vmovsd         (%3), %%xmm4			\n\t" // read y
	"vmovsd        8(%3), %%xmm5			\n\t"

	"vextractf128   $1, %%ymm8 , %%xmm9		      \n\t"
	"vextractf128   $1, %%ymm10, %%xmm11	      	      \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm11, %%xmm10      \n\t"

	"vshufpd        $0x1, %%xmm8 , %%xmm8 , %%xmm9   \n\t"
	"vshufpd        $0x1, %%xmm10, %%xmm10, %%xmm11  \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm11, %%xmm10      \n\t"

        "vmulps         %%xmm8 , %%xmm1 , %%xmm9              \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm8 , %%xmm0 , %%xmm8              \n\t"  // t_r * alpha_r , t_i * alpha_r
        "vmulps         %%xmm10, %%xmm1 , %%xmm11             \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm10, %%xmm0 , %%xmm10             \n\t"  // t_r * alpha_r , t_i * alpha_r

#if !defined(XCONJ)
        "vpermilps      $0xb1 , %%xmm9 , %%xmm9                \n\t"
        "vpermilps      $0xb1 , %%xmm11, %%xmm11               \n\t"
        "vaddsubps      %%xmm9 , %%xmm8, %%xmm8               \n\t"
        "vaddsubps      %%xmm11, %%xmm10, %%xmm10             \n\t"
#else
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
        "vpermilps      $0xb1 , %%xmm10, %%xmm10               \n\t"
        "vaddsubps      %%xmm8 , %%xmm9 , %%xmm8              \n\t"
        "vaddsubps      %%xmm10, %%xmm11, %%xmm10             \n\t"
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
        "vpermilps      $0xb1 , %%xmm10, %%xmm10               \n\t"
#endif


	"vaddps		%%xmm8 , %%xmm4 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm5 , %%xmm10      \n\t"

	"vmovsd	%%xmm8 ,   (%3)			\n\t"
	"vmovsd	%%xmm10,  8(%3)			\n\t"

	"vzeroupper			 \n\t"

	:
          "+r" (i),	// 0	
	  "+r" (n)  	// 1
	:
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (ap[0]),  // 4
          "r" (ap[1]),  // 5
          "r" (alpha)   // 6
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	 );

	    
     }
      else {

               __asm__ __volatile__ (

            "vzeroupper			 \n\t"

	"vxorps		%%ymm8 , %%ymm8 , %%ymm8 	\n\t" // temp
	"vxorps		%%ymm9 , %%ymm9 , %%ymm9 	\n\t" // temp
	"vxorps		%%ymm10, %%ymm10, %%ymm10	\n\t" // temp
	"vxorps		%%ymm11, %%ymm11, %%ymm11	\n\t" // temp

        "testq          $0x04, %1                      \n\t"
        "jz             2f                    \n\t"

	"vmovaps	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
	"vmovaps	(%5,%0,4), %%ymm5               \n\t" // 4 complex values from a1

	"vmovaps	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	
        "addq		$8  , %0	  	 	        \n\t"
	"subq	        $4  , %1			        \n\t"		

        "2:                                  \n\t"
	"cmpq           $0, %1                         \n\t"
        "je             3f                      \n\t"

		".align 16				        \n\t"
	"1:				        \n\t"
        "prefetcht0      256(%4,%0,4)                   \n\t"
	"vmovaps	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
        "prefetcht0      256(%5,%0,4)                   \n\t"
	"vmovaps	(%5,%0,4), %%ymm5               \n\t" // 4 complex values from a1

        "prefetcht0      256(%2,%0,4)                   \n\t"
	"vmovaps	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	
	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

	"vmovaps       32(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0
	"vmovaps       32(%5,%0,4), %%ymm5              \n\t" // 4 complex values from a1

	"vmovaps	  32(%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	"vfmadd231ps      %%ymm5 , %%ymm0, %%ymm10      \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm5 , %%ymm1, %%ymm11      \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

        "addq		$16 , %0	  	 	        \n\t"
	"subq	        $8  , %1			        \n\t"		
	"jnz		1b		        \n\t"

        "3:                                   \n\t"

        "vbroadcastss    (%6)  , %%xmm0                \n\t"  // value from alpha
        "vbroadcastss   4(%6)  , %%xmm1                \n\t"  // value from alpha


#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vpermilps      $0xb1 , %%ymm9 , %%ymm9                \n\t"
        "vpermilps      $0xb1 , %%ymm11, %%ymm11               \n\t"
        "vaddsubps      %%ymm9 , %%ymm8, %%ymm8                \n\t" 
        "vaddsubps      %%ymm11, %%ymm10, %%ymm10              \n\t"
#else
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
        "vpermilps      $0xb1 , %%ymm10, %%ymm10               \n\t"
        "vaddsubps      %%ymm8 , %%ymm9 , %%ymm8               \n\t"
        "vaddsubps      %%ymm10, %%ymm11, %%ymm10              \n\t"
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
        "vpermilps      $0xb1 , %%ymm10, %%ymm10               \n\t"
#endif

	"vmovsd         (%3), %%xmm4			\n\t" // read y
	"vmovsd        8(%3), %%xmm5			\n\t"

	"vextractf128   $1, %%ymm8 , %%xmm9		      \n\t"
	"vextractf128   $1, %%ymm10, %%xmm11	      	      \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm11, %%xmm10      \n\t"

	"vshufpd        $0x1, %%xmm8 , %%xmm8 , %%xmm9   \n\t"
	"vshufpd        $0x1, %%xmm10, %%xmm10, %%xmm11  \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm11, %%xmm10      \n\t"

        "vmulps         %%xmm8 , %%xmm1 , %%xmm9              \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm8 , %%xmm0 , %%xmm8              \n\t"  // t_r * alpha_r , t_i * alpha_r
        "vmulps         %%xmm10, %%xmm1 , %%xmm11             \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm10, %%xmm0 , %%xmm10             \n\t"  // t_r * alpha_r , t_i * alpha_r

#if !defined(XCONJ)
        "vpermilps      $0xb1 , %%xmm9 , %%xmm9                \n\t"
        "vpermilps      $0xb1 , %%xmm11, %%xmm11               \n\t"
        "vaddsubps      %%xmm9 , %%xmm8, %%xmm8               \n\t"
        "vaddsubps      %%xmm11, %%xmm10, %%xmm10             \n\t"
#else
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
        "vpermilps      $0xb1 , %%xmm10, %%xmm10               \n\t"
        "vaddsubps      %%xmm8 , %%xmm9 , %%xmm8              \n\t"
        "vaddsubps      %%xmm10, %%xmm11, %%xmm10             \n\t"
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
        "vpermilps      $0xb1 , %%xmm10, %%xmm10               \n\t"
#endif


	"vaddps		%%xmm8 , %%xmm4 , %%xmm8       \n\t"
	"vaddps		%%xmm10, %%xmm5 , %%xmm10      \n\t"

	"vmovsd	%%xmm8 ,   (%3)			\n\t"
	"vmovsd	%%xmm10,  8(%3)			\n\t"

	"vzeroupper			 \n\t"

	:
          "+r" (i),	// 0	
	  "+r" (n)  	// 1
	:
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (ap[0]),  // 4
          "r" (ap[1]),  // 5
          "r" (alpha)   // 6
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	 ); 
      }
}


 // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
void cgemv_kernel_4x1(const int32_t n,
                      float * __restrict ap,
		      float * __restrict x,
		      float * __restrict y,
		      float * __restrict alpha) {
          int32_t i = 0;
          bool isAligned32 = ((uintptr_t)x & 0x1F) != 0ULL && ((uintptr_t)y & 0x1F) != 0ULL;
	  if(!isAligned32) {

               __asm__ __volatile__ (

               "vzeroupper			 \n\t"

	"vxorps		%%ymm8 , %%ymm8 , %%ymm8 	\n\t" // temp
	"vxorps		%%ymm9 , %%ymm9 , %%ymm9 	\n\t" // temp

        "testq          $0x04, %1                      \n\t"
        "jz             2f                    \n\t"

	"vmovups	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0

	"vmovups	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	
        "addq		$8  , %0	  	 	        \n\t"
	"subq	        $4  , %1			        \n\t"		

        "2:                                  \n\t"
	"cmpq           $0, %1                         \n\t"
        "je             3f                      \n\t"

		".align 16				        \n\t"
	"1:				        \n\t"
        "prefetcht0      256(%4,%0,4)                   \n\t"
	"vmovups	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0

        "prefetcht0      256(%2,%0,4)                   \n\t"
	"vmovups	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	
	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

	"vmovups       32(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0

	"vmovups	  32(%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

        "addq		$16 , %0	  	 	        \n\t"
	"subq	        $8  , %1			        \n\t"		
	"jnz		1b		        \n\t"

        "3:                                   \n\t"

        "vbroadcastss    (%5)  , %%xmm0                \n\t"  // value from alpha
        "vbroadcastss   4(%5)  , %%xmm1                \n\t"  // value from alpha


#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vpermilps      $0xb1 , %%ymm9 , %%ymm9                \n\t"
        "vaddsubps      %%ymm9 , %%ymm8, %%ymm8                \n\t" 
#else
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
        "vaddsubps      %%ymm8 , %%ymm9 , %%ymm8               \n\t"
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
#endif

	"vmovsd         (%3), %%xmm4			\n\t" // read y

	"vextractf128   $1, %%ymm8 , %%xmm9		      \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"

	"vshufpd        $0x1, %%xmm8 , %%xmm8 , %%xmm9   \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"

        "vmulps         %%xmm8 , %%xmm1 , %%xmm9              \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm8 , %%xmm0 , %%xmm8              \n\t"  // t_r * alpha_r , t_i * alpha_r

#if !defined(XCONJ)
        "vpermilps      $0xb1 , %%xmm9 , %%xmm9                \n\t"
        "vaddsubps      %%xmm9 , %%xmm8, %%xmm8               \n\t"
#else
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
        "vaddsubps      %%xmm8 , %%xmm9 , %%xmm8              \n\t"
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
#endif


	"vaddps		%%xmm8 , %%xmm4 , %%xmm8       \n\t"

	"vmovsd	%%xmm8 ,   (%3)			\n\t"

	"vzeroupper			 \n\t"

	:
          "+r" (i),	// 0	
	  "+r" (n)  	// 1
	:
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (ap),     // 4
          "r" (alpha)   // 5
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	       );
	  }
	  else {

                 __asm__ __volatile__ (

               "vzeroupper			 \n\t"

	"vxorps		%%ymm8 , %%ymm8 , %%ymm8 	\n\t" // temp
	"vxorps		%%ymm9 , %%ymm9 , %%ymm9 	\n\t" // temp

        "testq          $0x04, %1                      \n\t"
        "jz             2f                    \n\t"

	"vmovaps	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0

	"vmovaps	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 
	
        "addq		$8  , %0	  	 	        \n\t"
	"subq	        $4  , %1			        \n\t"		

        "2:                                  \n\t"
	"cmpq           $0, %1                         \n\t"
        "je             3f                      \n\t"

		".align 16				        \n\t"
	"1:				        \n\t"
        "prefetcht0      256(%4,%0,4)                   \n\t"
	"vmovaps	(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0

        "prefetcht0      256(%2,%0,4)                   \n\t"
	"vmovaps	    (%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts
	
	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

	"vmovaps       32(%4,%0,4), %%ymm4	        \n\t" // 4 complex values from a0

	"vmovaps	  32(%2,%0,4)  , %%ymm6		\n\t" // 4 complex values from x
	"vpermilps        $0xb1, %%ymm6, %%ymm7		\n\t" // exchange real and imap parts
	"vblendps $0x55, %%ymm6, %%ymm7, %%ymm0         \n\t" // only the real parts
	"vblendps $0x55, %%ymm7, %%ymm6, %%ymm1         \n\t" // only the imag parts

	"vfmadd231ps      %%ymm4 , %%ymm0, %%ymm8       \n\t" // ar0*xr0,al0*xr0,ar1*xr1,al1*xr1 
	"vfmadd231ps      %%ymm4 , %%ymm1, %%ymm9       \n\t" // ar0*xl0,al0*xl0,ar1*xl1,al1*xl1 

        "addq		$16 , %0	  	 	        \n\t"
	"subq	        $8  , %1			        \n\t"		
	"jnz		1b		        \n\t"

        "3:                                   \n\t"

        "vbroadcastss    (%5)  , %%xmm0                \n\t"  // value from alpha
        "vbroadcastss   4(%5)  , %%xmm1                \n\t"  // value from alpha


#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vpermilps      $0xb1 , %%ymm9 , %%ymm9                \n\t"
        "vaddsubps      %%ymm9 , %%ymm8, %%ymm8                \n\t" 
#else
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
        "vaddsubps      %%ymm8 , %%ymm9 , %%ymm8               \n\t"
        "vpermilps      $0xb1 , %%ymm8 , %%ymm8                \n\t"
#endif

	"vmovsd         (%3), %%xmm4			\n\t" // read y

	"vextractf128   $1, %%ymm8 , %%xmm9		      \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"

	"vshufpd        $0x1, %%xmm8 , %%xmm8 , %%xmm9   \n\t"

	"vaddps		%%xmm8 , %%xmm9 , %%xmm8       \n\t"

        "vmulps         %%xmm8 , %%xmm1 , %%xmm9              \n\t"  // t_r * alpha_i , t_i * alpha_i
        "vmulps         %%xmm8 , %%xmm0 , %%xmm8              \n\t"  // t_r * alpha_r , t_i * alpha_r

#if !defined(XCONJ)
        "vpermilps      $0xb1 , %%xmm9 , %%xmm9                \n\t"
        "vaddsubps      %%xmm9 , %%xmm8, %%xmm8               \n\t"
#else
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
        "vaddsubps      %%xmm8 , %%xmm9 , %%xmm8              \n\t"
        "vpermilps      $0xb1 , %%xmm8 , %%xmm8                \n\t"
#endif


	"vaddps		%%xmm8 , %%xmm4 , %%xmm8       \n\t"

	"vmovsd	%%xmm8 ,   (%3)			\n\t"

	"vzeroupper			 \n\t"

	:
          "+r" (i),	// 0	
	  "+r" (n)  	// 1
	:
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (ap),     // 4
          "r" (alpha)   // 5
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	       );
	  }
}



#include <cstring> // memset
#define CGEMV_T_NBMAX 2048

 // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)

int32_t cgemv_t(const int32_t m,
             const int32_t n,
	     const float alpha_r,
	     const float alpha_i,
	     float * __restrict a,
	     const int32_t lda,
	     float * __restrict x,
	     const int32_t inc_x,
	     float * __restrict y,
	     const int32_t inc_y,
	     float * __restrict buffer) {
  
    __ATTR_ALIGN__(32) float *ap[8] = {};
    float * __restrict a_ptr = NULL;
    float * __restrict x_ptr = NULL;
    float * __restrict y_ptr = NULL;
    float * __restrict xbuffer = NULL;
    __ATTR_ALIGN__(32) float ybuffer[8] = {};
    float alpha[2] = {};
    int32_t i,j;
    int32_t n1,m1,m2,m3,n2,lda4;
      
      if ( m < 1 || n < 1) return(0);
      
        inc_x *= 2;
        inc_y *= 2;
        lda   *= 2;
	lda4    = lda << 2;

	xbuffer = buffer;
	
	n1 = n  >> 2 ;
	n2 = n  &  3 ;
	
	m3 = m & 3 ;
	m1 = m - m3;
	m2 = (m & (NBMAX-1)) - m3 ;
	
	alpha[0] = alpha_r;
	alpha[1] = alpha_i;

	int32_t NB = NBMAX;

	while ( NB == NBMAX )
	{
		
		m1 -= NB;
		if ( m1 < 0)
		{
			if ( m2 == 0 ) break;	
			NB = m2;
		}
		
		y_ptr = y;
		a_ptr = a;
		x_ptr = x;
		ap[0] = a_ptr;
		ap[1] = a_ptr + lda;
		ap[2] = ap[1] + lda;
		ap[3] = ap[2] + lda;
		if ( inc_x != 2 )
			copy_x(NB,x_ptr,xbuffer,inc_x);
		else
			xbuffer = x_ptr;
		
		if ( inc_y == 2 )
		{

			for( i = 0; i < n1 ; i++)
			{
				cgemv_kernel_4x4(NB,ap,xbuffer,y_ptr,alpha);
				ap[0] += lda4;
				ap[1] += lda4;
				ap[2] += lda4;
				ap[3] += lda4;
				a_ptr += lda4;
				y_ptr += 8;
				
			}

			if ( n2 & 2 )
			{
				cgemv_kernel_4x2(NB,ap,xbuffer,y_ptr,alpha);
				a_ptr += lda * 2;
				y_ptr += 4;

			}

			if ( n2 & 1 )
			{
				cgemv_kernel_4x1(NB,a_ptr,xbuffer,y_ptr,alpha);
				/* a_ptr += lda;
				y_ptr += 2; */

			}

		}
		else
		{

			for( i = 0; i < n1 ; i++)
			{
				memset(ybuffer,0,32);
				cgemv_kernel_4x4(NB,ap,xbuffer,ybuffer,alpha);
				ap[0] += lda4;
				ap[1] += lda4;
				ap[2] += lda4;
				ap[3] += lda4;
				a_ptr += lda4;

				y_ptr[0] += ybuffer[0];
				y_ptr[1] += ybuffer[1];
				y_ptr  += inc_y;
				y_ptr[0] += ybuffer[2];
				y_ptr[1] += ybuffer[3];
				y_ptr  += inc_y;
				y_ptr[0] += ybuffer[4];
				y_ptr[1] += ybuffer[5];
				y_ptr  += inc_y;
				y_ptr[0] += ybuffer[6];
				y_ptr[1] += ybuffer[7];
				y_ptr  += inc_y;

			}

			for( i = 0; i < n2 ; i++)
			{
				memset(ybuffer,0,32);
				cgemv_kernel_4x1(NB,a_ptr,xbuffer,ybuffer,alpha);
				a_ptr += lda;
				y_ptr[0] += ybuffer[0];
				y_ptr[1] += ybuffer[1];
				y_ptr  += inc_y;

			}

		}
		a += 2 * NB;
		x += NB * inc_x;	
	}



	if ( m3 == 0 ) return(0);

        x_ptr = x;
        j=0;
        a_ptr = a;
        y_ptr = y;

	if ( m3 == 3 )
	{

                float temp_r ;
                float  temp_i ;
		float  x0 = x_ptr[0];
		float  x1 = x_ptr[1];
		x_ptr += inc_x;
		float x2 = x_ptr[0];
		float x3 = x_ptr[1];
		x_ptr += inc_x;
		float x4 = x_ptr[0];
		float x5 = x_ptr[1];
	        while ( j < n)
        	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       	temp_r  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 + a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 - a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 + a_ptr[3] * x2; 
                       	temp_r += a_ptr[4] * x4 - a_ptr[5] * x5;
                       	temp_i += a_ptr[4] * x5 + a_ptr[5] * x4;
#else

                       	temp_r  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 - a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 + a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 - a_ptr[3] * x2; 
                       	temp_r += a_ptr[4] * x4 + a_ptr[5] * x5;
                       	temp_i += a_ptr[4] * x5 - a_ptr[5] * x4;
#endif

#if !defined(XCONJ) 
                	y_ptr[0] += alpha_r * temp_r - alpha_i * temp_i;
                	y_ptr[1] += alpha_r * temp_i + alpha_i * temp_r;
#else
                	y_ptr[0] += alpha_r * temp_r + alpha_i * temp_i;
                	y_ptr[1] -= alpha_r * temp_i - alpha_i * temp_r;
#endif

                	a_ptr += lda;
                	y_ptr += inc_y;
                	j++;
        	}
        	return(0);
	}


	if ( m3 == 2 )
	{

                float temp_r ;
                float temp_i ;
                float temp_r1 ;
                float temp_i1 ;
		float x0 = x_ptr[0];
		float x1 = x_ptr[1];
		x_ptr += inc_x;
		float x2 = x_ptr[0];
		float x3 = x_ptr[1];
		float ar = alpha[0];
		float ai = alpha[1];

	        while ( j < ( n & -2 ))
        	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       	temp_r  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 + a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 - a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 + a_ptr[3] * x2; 
                	a_ptr += lda;
                       	temp_r1  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i1  = a_ptr[0] * x1 + a_ptr[1] * x0; 
                       	temp_r1 += a_ptr[2] * x2 - a_ptr[3] * x3; 
                       	temp_i1 += a_ptr[2] * x3 + a_ptr[3] * x2; 
#else

                       	temp_r  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 - a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 + a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 - a_ptr[3] * x2; 
                	a_ptr += lda;
                       	temp_r1  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i1  = a_ptr[0] * x1 - a_ptr[1] * x0; 
                       	temp_r1 += a_ptr[2] * x2 + a_ptr[3] * x3; 
                       	temp_i1 += a_ptr[2] * x3 - a_ptr[3] * x2; 
#endif

#if !defined(XCONJ) 
                	y_ptr[0] += ar * temp_r - ai * temp_i;
                	y_ptr[1] += ar * temp_i + ai * temp_r;
                	y_ptr += inc_y;
                	y_ptr[0] += ar * temp_r1 - ai * temp_i1;
                	y_ptr[1] += ar * temp_i1 + ai * temp_r1;
#else
                	y_ptr[0] += ar * temp_r + ai * temp_i;
                	y_ptr[1] -= ar * temp_i - ai * temp_r;
                	y_ptr += inc_y;
                	y_ptr[0] += ar * temp_r1 + ai * temp_i1;
                	y_ptr[1] -= ar * temp_i1 - ai * temp_r1;
#endif

                	a_ptr += lda;
                	y_ptr += inc_y;
                	j+=2;
        	}


	        while ( j < n)
        	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       	temp_r  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 + a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 - a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 + a_ptr[3] * x2; 
#else

                       	temp_r  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 - a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 + a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 - a_ptr[3] * x2; 
#endif

#if !defined(XCONJ) 
                	y_ptr[0] += ar * temp_r - ai * temp_i;
                	y_ptr[1] += ar * temp_i + ai * temp_r;
#else
                	y_ptr[0] += ar * temp_r + ai * temp_i;
                	y_ptr[1] -= ar * temp_i - ai * temp_r;
#endif

                	a_ptr += lda;
                	y_ptr += inc_y;
                	j++;
        	}

        	return(0);
	}


	if ( m3 == 1 )
	{

                float temp_r ;
                float temp_i ;
                float temp_r1 ;
                float temp_i1 ;
		float x0 = x_ptr[0];
		float x1 = x_ptr[1];
		float ar = alpha[0];
		float ai = alpha[1];

	        while ( j < ( n & -2 ))
        	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       	temp_r  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 + a_ptr[1] * x0; 
                	a_ptr += lda;
                       	temp_r1  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i1  = a_ptr[0] * x1 + a_ptr[1] * x0; 
#else

                       	temp_r  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 - a_ptr[1] * x0; 
                	a_ptr += lda;
                       	temp_r1  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i1  = a_ptr[0] * x1 - a_ptr[1] * x0; 
#endif

#if !defined(XCONJ) 
                	y_ptr[0] += ar * temp_r - ai * temp_i;
                	y_ptr[1] += ar * temp_i + ai * temp_r;
                	y_ptr += inc_y;
                	y_ptr[0] += ar * temp_r1 - ai * temp_i1;
                	y_ptr[1] += ar * temp_i1 + ai * temp_r1;
#else
                	y_ptr[0] += ar * temp_r + ai * temp_i;
                	y_ptr[1] -= ar * temp_i - ai * temp_r;
                	y_ptr += inc_y;
                	y_ptr[0] += ar * temp_r1 + ai * temp_i1;
                	y_ptr[1] -= ar * temp_i1 - ai * temp_r1;
#endif

                	a_ptr += lda;
                	y_ptr += inc_y;
                	j+=2;
        	}

	        while ( j < n)
        	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       	temp_r  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 + a_ptr[1] * x0; 
#else

                       	temp_r  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 - a_ptr[1] * x0; 
#endif

#if !defined(XCONJ) 
                	y_ptr[0] += ar * temp_r - ai * temp_i;
                	y_ptr[1] += ar * temp_i + ai * temp_r;
#else
                	y_ptr[0] += ar * temp_r + ai * temp_i;
                	y_ptr[1] -= ar * temp_i - ai * temp_r;
#endif

                	a_ptr += lda;
                	y_ptr += inc_y;
                	j++;
        	}
        	return(0);
	}

	return(0);

}


