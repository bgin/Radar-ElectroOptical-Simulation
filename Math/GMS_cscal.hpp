

#ifndef __GMS_CSCAL_HPP__
#define __GMS_CSCAL_HPP__


//=========================================================================
// Modified and optimized version of OpenBLAS cscal and kernels cscal_kernel_x
// Programmer: Bernard Gingold, contact: beniekg@gmail.com
// 27-02-2021 09:58  +00200
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
#include "GMS_config.h"

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
__attribute__((no_stack_protector))
static inline
void cscal_kernel_16(const int32_t n,
                     float * __restrict alpha,
		     float * __restrict x) {

      const bool isAligned32 = ((uintptr_t)alpha & 0x1F) != 0ULL && ((uintptr_t)y & 0x1F) != 0ULL;
      if(!isAligned32) {

          __asm__ __volatile__ (

               	"vbroadcastss		(%2), %%ymm0		    \n\t"  // da_r	
	"vbroadcastss          4(%2), %%ymm1		    \n\t"  // da_i 	

	"addq	$128, %1				    \n\t"

	"vmovups	-128(%1), %%ymm4		    \n\t"
	"vmovups	 -96(%1), %%ymm5		    \n\t"
	"vmovups	 -64(%1), %%ymm6		    \n\t"
	"vmovups	 -32(%1), %%ymm7		    \n\t"

	"vpermilps	$0xb1 , %%ymm4, %%ymm12		    \n\t"
	"vpermilps	$0xb1 , %%ymm5, %%ymm13		    \n\t"
	"vpermilps	$0xb1 , %%ymm6, %%ymm14 	    \n\t"
	"vpermilps	$0xb1 , %%ymm7, %%ymm15		    \n\t"

	"subq	        $16, %0			            \n\t"		
	"jz	2f					    \n\t"

	".p2align 4				            \n\t"
	"1:				            	    \n\t"
#if (GMS_MAN_PREFETCH) == 1
	"prefetcht0     256(%1)				    \n\t"
#endif
	// ".align 2				            \n\t"

	"vmulps		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmovups	   0(%1), %%ymm4		    \n\t"
	"vmulps		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmovups	  32(%1), %%ymm5		    \n\t"
	"vmulps		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	"vmovups	  64(%1), %%ymm6		    \n\t"
	"vmulps		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 
	"vmovups	  96(%1), %%ymm7		    \n\t"

	"vmulps		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	"vaddsubps	%%ymm12 , %%ymm8 , %%ymm8	    \n\t"
	"vmulps		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	"vaddsubps	%%ymm13 , %%ymm9 , %%ymm9	    \n\t"
	"vmulps		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	"vaddsubps	%%ymm14 , %%ymm10, %%ymm10	    \n\t"
	"vmulps		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	"vaddsubps	%%ymm15 , %%ymm11, %%ymm11	    \n\t"

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm4, %%ymm12		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm5, %%ymm13		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm6, %%ymm14 	    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm7, %%ymm15		    \n\t"

	"addq		$128 ,%1  	 	            \n\t"
	"subq	        $16, %0			            \n\t"		
	"jnz		1b		             	    \n\t"

	"2:				            	    \n\t"


	"vmulps		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmulps		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmulps		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	"vmulps		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 

	"vmulps		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	"vaddsubps	%%ymm12 , %%ymm8 , %%ymm8	    \n\t"
	"vmulps		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	"vaddsubps	%%ymm13 , %%ymm9 , %%ymm9	    \n\t"
	"vmulps		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	"vaddsubps	%%ymm14 , %%ymm10, %%ymm10	    \n\t"
	"vmulps		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	"vaddsubps	%%ymm15 , %%ymm11, %%ymm11	    \n\t"

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"

	"vzeroupper					    \n\t"

	:
	  "+r" (n),  	// 0
          "+r" (x)      // 1
        :
          "r" (alpha)   // 2
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

               	"vbroadcastss		(%2), %%ymm0		    \n\t"  // da_r	
	"vbroadcastss          4(%2), %%ymm1		    \n\t"  // da_i 	

	"addq	$128, %1				    \n\t"

	"vmovaps	-128(%1), %%ymm4		    \n\t"
	"vmovaps	 -96(%1), %%ymm5		    \n\t"
	"vmovaps	 -64(%1), %%ymm6		    \n\t"
	"vmovaps	 -32(%1), %%ymm7		    \n\t"

	"vpermilps	$0xb1 , %%ymm4, %%ymm12		    \n\t"
	"vpermilps	$0xb1 , %%ymm5, %%ymm13		    \n\t"
	"vpermilps	$0xb1 , %%ymm6, %%ymm14 	    \n\t"
	"vpermilps	$0xb1 , %%ymm7, %%ymm15		    \n\t"

	"subq	        $16, %0			            \n\t"		
	"jz	2f					    \n\t"

	".p2align 4				            \n\t"
	"1:				            	    \n\t"
#if (GMS_MAN_PREFETCH) == 1
	"prefetcht0     256(%1)				    \n\t"
#endif
	// ".align 2				            \n\t"

	"vmulps		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmovaps	   0(%1), %%ymm4		    \n\t"
	"vmulps		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmovaps	  32(%1), %%ymm5		    \n\t"
	"vmulps		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	"vmovaps	  64(%1), %%ymm6		    \n\t"
	"vmulps		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 
	"vmovaps	  96(%1), %%ymm7		    \n\t"

	"vmulps		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	"vaddsubps	%%ymm12 , %%ymm8 , %%ymm8	    \n\t"
	"vmulps		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	"vaddsubps	%%ymm13 , %%ymm9 , %%ymm9	    \n\t"
	"vmulps		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	"vaddsubps	%%ymm14 , %%ymm10, %%ymm10	    \n\t"
	"vmulps		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	"vaddsubps	%%ymm15 , %%ymm11, %%ymm11	    \n\t"

	"vmovaps	%%ymm8 , -128(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm4, %%ymm12		    \n\t"
	"vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm5, %%ymm13		    \n\t"
	"vmovaps	%%ymm10,  -64(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm6, %%ymm14 	    \n\t"
	"vmovaps	%%ymm11,  -32(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm7, %%ymm15		    \n\t"

	"addq		$128 ,%1  	 	            \n\t"
	"subq	        $16, %0			            \n\t"		
	"jnz		1b		             	    \n\t"

	"2:				            	    \n\t"


	"vmulps		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmulps		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmulps		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	"vmulps		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 

	"vmulps		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	"vaddsubps	%%ymm12 , %%ymm8 , %%ymm8	    \n\t"
	"vmulps		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	"vaddsubps	%%ymm13 , %%ymm9 , %%ymm9	    \n\t"
	"vmulps		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	"vaddsubps	%%ymm14 , %%ymm10, %%ymm10	    \n\t"
	"vmulps		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	"vaddsubps	%%ymm15 , %%ymm11, %%ymm11	    \n\t"

	"vmovaps	%%ymm8 , -128(%1)		    \n\t"
	"vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovaps	%%ymm10,  -64(%1)		    \n\t"
	"vmovaps	%%ymm11,  -32(%1)		    \n\t"

	"vzeroupper					    \n\t"

	:
	  "+r" (n),  	// 0
          "+r" (x)      // 1
        :
          "r" (alpha)   // 2
	: "cc",
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	  );
      }
}


__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
__attribute__((no_stack_protector))
static inline
void cscal_kernel_16_zero_r(const int32_t n,
                            float * __restrict alpha,
			    float * __restrict x) {
      const bool isALigned32 = ((uintptr_t)alpha & 0x1F) != 0ULL && ((uintptr_t)y & 0x1F) != 0ULL;
      if(!isAligned32) {

           __asm__ __volatile__ (

         "vxorps	           %%ymm0, %%ymm0, %%ymm0	    \n\t"
	"vbroadcastss          4(%2), %%ymm1		    \n\t"  // da_i 	

	"addq	$128, %1				    \n\t"

	"vmovups	-128(%1), %%ymm4		    \n\t"
	"vmovups	 -96(%1), %%ymm5		    \n\t"
	"vmovups	 -64(%1), %%ymm6		    \n\t"
	"vmovups	 -32(%1), %%ymm7		    \n\t"

	"vpermilps	$0xb1 , %%ymm4, %%ymm12		    \n\t"
	"vpermilps	$0xb1 , %%ymm5, %%ymm13		    \n\t"
	"vpermilps	$0xb1 , %%ymm6, %%ymm14 	    \n\t"
	"vpermilps	$0xb1 , %%ymm7, %%ymm15		    \n\t"

	"subq	        $16, %0			            \n\t"		
	"jz	2f					    \n\t"

	".p2align 4				            \n\t"
	"1:				            	    \n\t"
#if (GMS_MAN_PREFETCH) == 1
	"prefetcht0     128(%1)				    \n\t"
#endif
	// ".align 2				            \n\t"

	"vmovups	   0(%1), %%ymm4		    \n\t"
	"vmovups	  32(%1), %%ymm5		    \n\t"
	"vmovups	  64(%1), %%ymm6		    \n\t"
	"vmovups	  96(%1), %%ymm7		    \n\t"

	"vmulps		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	"vaddsubps	%%ymm12 , %%ymm0 , %%ymm8	    \n\t"
	"vmulps		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	"vaddsubps	%%ymm13 , %%ymm0 , %%ymm9	    \n\t"
	"vmulps		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	"vaddsubps	%%ymm14 , %%ymm0 , %%ymm10	    \n\t"
	"vmulps		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	"vaddsubps	%%ymm15 , %%ymm0 , %%ymm11	    \n\t"

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm4, %%ymm12		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm5, %%ymm13		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm6, %%ymm14 	    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm7, %%ymm15		    \n\t"

	"addq		$128 ,%1  	 	            \n\t"
	"subq	        $16, %0			            \n\t"		
	"jnz		1b		             	    \n\t"

	"2:				            	    \n\t"

	"vmulps		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	"vaddsubps	%%ymm12 , %%ymm0 , %%ymm8	    \n\t"
	"vmulps		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	"vaddsubps	%%ymm13 , %%ymm0 , %%ymm9	    \n\t"
	"vmulps		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	"vaddsubps	%%ymm14 , %%ymm0 , %%ymm10	    \n\t"
	"vmulps		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	"vaddsubps	%%ymm15 , %%ymm0 , %%ymm11	    \n\t"

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"

	"vzeroupper					    \n\t"

	:
	  "+r" (n),  	// 0
          "+r" (x)      // 1
        :
          "r" (alpha)   // 2
	: "cc", // "0", "1",
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	   );
      }
      else {

          __asm__ __volatile__ (

         "vxorps	           %%ymm0, %%ymm0, %%ymm0	    \n\t"
	"vbroadcastss          4(%2), %%ymm1		    \n\t"  // da_i 	

	"addq	$128, %1				    \n\t"

	"vmovaps	-128(%1), %%ymm4		    \n\t"
	"vmovaps	 -96(%1), %%ymm5		    \n\t"
	"vmovaps	 -64(%1), %%ymm6		    \n\t"
	"vmovaps	 -32(%1), %%ymm7		    \n\t"

	"vpermilps	$0xb1 , %%ymm4, %%ymm12		    \n\t"
	"vpermilps	$0xb1 , %%ymm5, %%ymm13		    \n\t"
	"vpermilps	$0xb1 , %%ymm6, %%ymm14 	    \n\t"
	"vpermilps	$0xb1 , %%ymm7, %%ymm15		    \n\t"

	"subq	        $16, %0			            \n\t"		
	"jz	2f					    \n\t"

	".p2align 4				            \n\t"
	"1:				            	    \n\t"
#if (GMS_MAN_PREFETCH) == 1
	"prefetcht0     128(%1)				    \n\t"
#endif
	// ".align 2				            \n\t"

	"vmovaps	   0(%1), %%ymm4		    \n\t"
	"vmovaps	  32(%1), %%ymm5		    \n\t"
	"vmovaps	  64(%1), %%ymm6		    \n\t"
	"vmovaps	  96(%1), %%ymm7		    \n\t"

	"vmulps		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	"vaddsubps	%%ymm12 , %%ymm0 , %%ymm8	    \n\t"
	"vmulps		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	"vaddsubps	%%ymm13 , %%ymm0 , %%ymm9	    \n\t"
	"vmulps		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	"vaddsubps	%%ymm14 , %%ymm0 , %%ymm10	    \n\t"
	"vmulps		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	"vaddsubps	%%ymm15 , %%ymm0 , %%ymm11	    \n\t"

	"vmovaps	%%ymm8 , -128(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm4, %%ymm12		    \n\t"
	"vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm5, %%ymm13		    \n\t"
	"vmovaps	%%ymm10,  -64(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm6, %%ymm14 	    \n\t"
	"vmovaps	%%ymm11,  -32(%1)		    \n\t"
	"vpermilps	$0xb1 , %%ymm7, %%ymm15		    \n\t"

	"addq		$128 ,%1  	 	            \n\t"
	"subq	        $16, %0			            \n\t"		
	"jnz		1b		             	    \n\t"

	"2:				            	    \n\t"

	"vmulps		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	"vaddsubps	%%ymm12 , %%ymm0 , %%ymm8	    \n\t"
	"vmulps		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	"vaddsubps	%%ymm13 , %%ymm0 , %%ymm9	    \n\t"
	"vmulps		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	"vaddsubps	%%ymm14 , %%ymm0 , %%ymm10	    \n\t"
	"vmulps		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	"vaddsubps	%%ymm15 , %%ymm0 , %%ymm11	    \n\t"

	"vmovaps	%%ymm8 , -128(%1)		    \n\t"
	"vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovaps	%%ymm10,  -64(%1)		    \n\t"
	"vmovaps	%%ymm11,  -32(%1)		    \n\t"

	"vzeroupper					    \n\t"

	:
	  "+r" (n),  	// 0
          "+r" (x)      // 1
        :
          "r" (alpha)   // 2
	: "cc", // "0", "1",
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	   );
      }
}


__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
__attribute__((no_stack_protector))
static inline
void cscal_kernel_16_zero_i(const int32_t n,
                            float * __restrict alpha,
                            float * __restrict x) {

           const bool isAligned32 =
	              ((uintptr_t)alpha & 0x1F) != 0ULL && ((uintptr_t)y & 0x1F) != 0ULL;
		      
	   if(!isAligned32) {

           __asm__  __volatile__ (
	
	"vbroadcastss		(%2), %%ymm0		    \n\t"  // da_r	

	"addq	$128, %1				    \n\t"

	"vmovups	-128(%1), %%ymm4		    \n\t"
	"vmovups	 -96(%1), %%ymm5		    \n\t"
	"vmovups	 -64(%1), %%ymm6		    \n\t"
	"vmovups	 -32(%1), %%ymm7		    \n\t"


	"subq	        $16, %0			            \n\t"		
	"jz	2f					    \n\t"

	".p2align 4				            \n\t"
	"1:				            	    \n\t"
#if (GMS_MAN_PREFETCH) == 1
	"prefetcht0     256(%1)				    \n\t"
#endif
	// ".align 2				            \n\t"

	"vmulps		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmovups	   0(%1), %%ymm4		    \n\t"
	"vmulps		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmovups	  32(%1), %%ymm5		    \n\t"
	"vmulps		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	"vmovups	  64(%1), %%ymm6		    \n\t"
	"vmulps		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 
	"vmovups	  96(%1), %%ymm7		    \n\t"

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"

	"addq		$128 ,%1  	 	            \n\t"
	"subq	        $16, %0			            \n\t"		
	"jnz		1b		             	    \n\t"

	"2:				            	    \n\t"


	"vmulps		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmulps		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmulps		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	"vmulps		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"

	"vzeroupper					    \n\t"

	:
	  "+r" (n),  	// 0
          "+r" (x)      // 1
        :
          "r" (alpha)   // 2
	: "cc", //"%0", "%1",
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

     }
     else {

             __asm__  __volatile__ (
	
	"vbroadcastss		(%2), %%ymm0		    \n\t"  // da_r	

	"addq	$128, %1				    \n\t"

	"vmovaps	-128(%1), %%ymm4		    \n\t"
	"vmovaps	 -96(%1), %%ymm5		    \n\t"
	"vmovaps	 -64(%1), %%ymm6		    \n\t"
	"vmovaps	 -32(%1), %%ymm7		    \n\t"


	"subq	        $16, %0			            \n\t"		
	"jz	2f					    \n\t"

	".p2align 4				            \n\t"
	"1:				            	    \n\t"
#if (GMS_MAN_PREFETCH) == 1
	"prefetcht0     256(%1)				    \n\t"
#endif
	// ".align 2				            \n\t"

	"vmulps		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmovaps	   0(%1), %%ymm4		    \n\t"
	"vmulps		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmovaps	  32(%1), %%ymm5		    \n\t"
	"vmulps		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	"vmovaps	  64(%1), %%ymm6		    \n\t"
	"vmulps		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 
	"vmovaps	  96(%1), %%ymm7		    \n\t"

	"vmovaps	%%ymm8 , -128(%1)		    \n\t"
	"vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovaps	%%ymm10,  -64(%1)		    \n\t"
	"vmovaps	%%ymm11,  -32(%1)		    \n\t"

	"addq		$128 ,%1  	 	            \n\t"
	"subq	        $16, %0			            \n\t"		
	"jnz		1b		             	    \n\t"

	"2:				            	    \n\t"


	"vmulps		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmulps		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmulps		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	"vmulps		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 

	"vmovaps	%%ymm8 , -128(%1)		    \n\t"
	"vmovaps	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovaps	%%ymm10,  -64(%1)		    \n\t"
	"vmovaps	%%ymm11,  -32(%1)		    \n\t"

	"vzeroupper					    \n\t"

	:
	  "+r" (n),  	// 0
          "+r" (x)      // 1
        :
          "r" (alpha)   // 2
	: "cc", //"%0", "%1",
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

     }
}

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
__attribute__((no_stack_protector))
static inline
void cscal_kernel_16_zero(const int32_t n,
                          float * __restrict alpha,
			  float * __restrict x) {

         const bool isAligned32 =
	          ((uintptr_t)alpha & 0x1F) != 0ULL && ((uintptr_t)y & 0x1F) != 0ULL;

	 if(!isAligned32) {

	       __asm__ __volatile__ (

                    	"vxorps	           %%ymm0, %%ymm0, %%ymm0	    \n\t"

	                "addq	$128, %1				    \n\t"

	                ".p2align 4				            \n\t"
	                "1:				            	    \n\t"

	              //  "prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"

	                "vmovups	%%ymm0 , -128(%1)		    \n\t"
	                "vmovups	%%ymm0 ,  -96(%1)		    \n\t"
	                "vmovups	%%ymm0 ,  -64(%1)		    \n\t"
	                "vmovups	%%ymm0 ,  -32(%1)		    \n\t"

	                "addq		$128 ,%1  	 	            \n\t"
	                "subq	        $16, %0			            \n\t"		
	                "jnz		1b		             	    \n\t"

	                "vzeroupper					    \n\t"

	                : 
	                "+r" (n),  	// 0
                        "+r" (x)      // 1
                        :
                        "r" (alpha)   // 2
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

                    	"vxorps	           %%ymm0, %%ymm0, %%ymm0	    \n\t"

	                "addq	$128, %1				    \n\t"

	                ".p2align 4				            \n\t"
	                "1:				            	    \n\t"

	              //  "prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"

	                "vmovaps	%%ymm0 , -128(%1)		    \n\t"
	                "vmovaps	%%ymm0 ,  -96(%1)		    \n\t"
	                "vmovaps	%%ymm0 ,  -64(%1)		    \n\t"
	                "vmovaps	%%ymm0 ,  -32(%1)		    \n\t"

	                "addq		$128 ,%1  	 	            \n\t"
	                "subq	        $16, %0			            \n\t"		
	                "jnz		1b		             	    \n\t"

	                "vzeroupper					    \n\t"

	                : 
	                "+r" (n),  	// 0
                        "+r" (x)      // 1
                        :
                        "r" (alpha)   // 2
	                : "cc",
	                "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	                "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	                "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	                "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	                "memory"
	       ); 
	 }
}

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
__attribute__((no_stack_protector))
static inline
void cscal_kernel_inc_8(const int32_t n,
                        float * __restrict alpha,
			float * __restrict x,
			const int32_t inc_x) {
      float da_r = alpha[0];
      float da_i = alpha[1];
      float t0,t1,t2,t3;
      const int32_t inc_x2 = inc_x+inc_x;
      const int32_t inc_x3 = inc_x2+inc_x;
      int32_t i;
      for( i=0; i<n; i+=4 ){
	
		t0 = da_r * x[0]         - da_i *x[1];	
		t1 = da_r * x[inc_x]     - da_i *x[inc_x  + 1];	
		t2 = da_r * x[inc_x2]    - da_i *x[inc_x2 + 1];	
		t3 = da_r * x[inc_x3]    - da_i *x[inc_x3 + 1];	

		x[1]               = da_i * x[0]       + da_r * x[1];
		x[inc_x  +1]       = da_i * x[inc_x]   + da_r * x[inc_x  +1];
		x[inc_x2 +1]       = da_i * x[inc_x2]  + da_r * x[inc_x2 +1];
		x[inc_x3 +1]       = da_i * x[inc_x3]  + da_r * x[inc_x3 +1];

		x[0]        = t0;
		x[inc_x]    = t1;
		x[inc_x2]   = t2;
		x[inc_x3]   = t3;

		x+=4*inc_x;

	}
}


__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
__attribute__((no_stack_protector))
static inline
int32_t  cscal(const int32_t n,
               float da_r,
	       float da_i,
	       float * __restrict x,
	       const int32_t inc_x,
               float * __restrict y,
               const int32_t inc_y) {

      float alpha[2] = {};
      float temp0,temp1;
      int32_t i = 0,j = 0;
      if ( inc_x != 1 )
	{
		inc_x <<= 1;

		if ( da_r == 0.0f )
		{

			const int32_t n1 = n & -2;

			if ( da_i == 0.0f )
			{

				while(j < n1)
				{
			
					x[i]=0.0;
					x[i+1]=0.0;
					x[i+inc_x]=0.0;
					x[i+1+inc_x]=0.0;
					i += 2*inc_x ;
					j+=2;

				}

				while(j < n)
				{
			
					x[i]=0.0;
					x[i+1]=0.0;
					i += inc_x ;
					j++;

				}

			}
			else
			{

				while(j < n1)
				{
			
					temp0        = -da_i * x[i+1];
					x[i+1]       =  da_i * x[i];
					x[i]         =  temp0;
					temp1        = -da_i * x[i+1+inc_x];
					x[i+1+inc_x] =  da_i * x[i+inc_x];
					x[i+inc_x]   =  temp1;
					i += 2*inc_x ;
					j+=2;

				}

				while(j < n)
				{
			
					temp0        = -da_i * x[i+1];
					x[i+1]       =  da_i * x[i];
					x[i]         =  temp0;
					i += inc_x ;
					j++;

				}



			}

		}
		else
		{


			if ( da_i == 0.0f )
			{
				const int32_t n1 = n & -2;

				while(j < n1)
				{
			
					temp0        =  da_r * x[i];
					x[i+1]       =  da_r * x[i+1];
					x[i]         =  temp0;
					temp1        =  da_r * x[i+inc_x];
					x[i+1+inc_x] =  da_r * x[i+1+inc_x];
					x[i+inc_x]   =  temp1;
					i += 2*inc_x ;
					j+=2;

				}

				while(j < n)
				{
			
					temp0        =  da_r * x[i];
					x[i+1]       =  da_r * x[i+1];
					x[i]         =  temp0;
					i += inc_x ;
					j++;

				}

			}
			else
			{

				const int32_t n1 = n & -8;
				if ( n1 > 0 )
				{
					alpha[0] = da_r;
			                alpha[1] = da_i;
					cscal_kernel_inc_8(n1, alpha, x, inc_x);
					j = n1 ;
					i = n1 * inc_x;
				}

				while(j < n)
				{

					temp0        =  da_r * x[i]   - da_i * x[i+1];
					x[i+1]       =  da_r * x[i+1] + da_i * x[i];
					x[i]         =  temp0;
					i += inc_x ;
					j++;

				}

			}

		}

		return(0);
	}


	const int32_t n1 = n & -16;
	if ( n1 > 0 )
	{

		alpha[0] = da_r;
		alpha[1] = da_i;
	
		if ( da_r == 0.0f )
			if ( da_i == 0.0f )
				cscal_kernel_16_zero(n1 , alpha , x);
			else
				cscal_kernel_16_zero_r(n1 , alpha , x);
		else
			if ( da_i == 0 )
				cscal_kernel_16_zero_i(n1 , alpha , x);
			else
				cscal_kernel_16(n1 , alpha , x);

		i = n1 << 1;
		j = n1;
	}


	if ( da_r == 0.0f )
	{

		if ( da_i == 0.0f )
		{

			while(j < n)
			{
		
					x[i]=0.0f;
					x[i+1]=0.0f;
					i += 2 ;
					j++;

			}

		}
		else
		{

			while(j < n)
			{
			
				temp0        = -da_i * x[i+1];
				x[i+1]       =  da_i * x[i];
				x[i]         =  temp0;
				i += 2 ;
				j++;

			}

		}

	}
	else
	{

		if ( da_i == 0.0f )
		{

			while(j < n)
			{
			
					temp0        =  da_r * x[i];
					x[i+1]       =  da_r * x[i+1];
					x[i]         =  temp0;
					i += 2 ;
					j++;

			}

		}
		else
		{

			const int32_t n2 = n & -2;

			while(j < n2)
			{

				temp0        =  da_r * x[i]   - da_i * x[i+1];
				temp1        =  da_r * x[i+2] - da_i * x[i+3];
				x[i+1]       =  da_r * x[i+1] + da_i * x[i];
				x[i+3]       =  da_r * x[i+3] + da_i * x[i+2];
				x[i]         =  temp0;
				x[i+2]       =  temp1;
				i += 4 ;
				j+=2;

			}

			while(j < n)
			{

				temp0        =  da_r * x[i]   - da_i * x[i+1];
				x[i+1]       =  da_r * x[i+1] + da_i * x[i];
				x[i]         =  temp0;
				i += 2 ;
				j++;

			}

		}		

	}

	return(0);
}


	       

#endif /*__GMS_CSCAL_HPP__*/
