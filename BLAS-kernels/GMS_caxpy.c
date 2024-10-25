



//=========================================================================
// Modified and optimized version of OpenBLAS caxpy and kernel caxpy_kernel_8
// Programmer: Bernard Gingold, contact: beniekg@gmail.com
// 13-09-2020 11:14 AM +00200
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


#include "GMS_caxpy.h" // must define 'CONJ'

// TODO:
// Count the number of uops and implement version for SKX,CLK CPUs.
__ATTR_HOT__

void caxpy_kernel_8(const int32_t n,
		    float * __restrict __ATTR_ALIGN__(32) x,
		    const  float * __restrict __ATTR_ALIGN__(32) y,
		    float * __restrict  __ATTR_ALIGN__(8) alpha) {
#if (CONJ) == 1
  __attribute__((aligned(32)))float mvec[8] =
                       {-1.0f,1.0f,-1.0f,1.0f,-1.0f,1.0f,-1.0f,1.0f};
#else
  __attribute__((aligned(32)))float mvec[8] = 
                       {1.0f,-1.0f,1.0f,-1.0f,1.0f,-1.0f,1.0f,-1.0f};
#endif
      int32_t i = 0;
      if(((uintptr_t)x & 0x1F) != 0ULL &&
              ((uintptr_t)y & 0x1F) != 0ULL ) {

	   __asm__ __volatile__ (

              "vzeroupper					    \n\t"
	      "vbroadcastss    (%4), %%ymm0		    \n\t"  // real part of alpha
	      "vbroadcastss    4(%4), %%ymm1		    \n\t"  // imag part of alpha
#if (CONJ) == 1
   	     "vmulps		(%5), %%ymm1 , %%ymm1		    \n\t"
#else
	     "vmulps		(%5), %%ymm0 , %%ymm0		    \n\t"
#endif

	     ".p2align 4				            \n\t"
	     "1:				            \n\t"

	     "vmovups        (%2,%0,4), %%ymm5                   \n\t" // 4 complex values from x
	     ".p2align 1					    \n\t"
	     "vmovups      32(%2,%0,4), %%ymm7                   \n\t" // 4 complex values from x
	     "vmovups      64(%2,%0,4), %%ymm9                   \n\t" // 4 complex values from x
	     "vmovups      96(%2,%0,4), %%ymm11                  \n\t" // 4 complex values from x

	     "vmovups     128(%2,%0,4), %%ymm12                  \n\t" // 4 complex values from x
	     "vmovups     160(%2,%0,4), %%ymm13                  \n\t" // 4 complex values from x
	     "vmovups     192(%2,%0,4), %%ymm14                  \n\t" // 4 complex values from x
	     "vmovups     224(%2,%0,4), %%ymm15                  \n\t" // 4 complex values from x

	     "vpermilps	$0xb1 , %%ymm5 , %%ymm4 	    \n\t"  // exchange real and imag part
	     "vpermilps	$0xb1 , %%ymm7 , %%ymm6 	    \n\t"  // exchange real and imag part
	     "vpermilps	$0xb1 , %%ymm9 , %%ymm8 	    \n\t"  // exchange real and imag part
	     "vpermilps	$0xb1 , %%ymm11, %%ymm10 	    \n\t"  // exchange real and imag part

	     "vfmadd213ps    (%3,%0,4), %%ymm0 , %%ymm5          \n\t"
             ".p2align 1					    \n\t"
	     "vfmadd213ps  32(%3,%0,4), %%ymm0 , %%ymm7          \n\t"
	     "vfmadd213ps  64(%3,%0,4), %%ymm0 , %%ymm9          \n\t"
	     "vfmadd213ps  96(%3,%0,4), %%ymm0 , %%ymm11         \n\t"

	     "vfmadd231ps	%%ymm1 , %%ymm4 , %%ymm5   \n\t"
	     "vfmadd231ps	%%ymm1 , %%ymm6 , %%ymm7   \n\t"
	     "vfmadd231ps	%%ymm1 , %%ymm8 , %%ymm9   \n\t"
	     "vfmadd231ps	%%ymm1 , %%ymm10, %%ymm11  \n\t"

	     "vpermilps	$0xb1 , %%ymm12, %%ymm4 	    \n\t"  // exchange real and imag part
	     "vpermilps	$0xb1 , %%ymm13, %%ymm6 	    \n\t"  // exchange real and imag part
	     "vpermilps	$0xb1 , %%ymm14, %%ymm8 	    \n\t"  // exchange real and imag part
	     "vpermilps	$0xb1 , %%ymm15, %%ymm10 	    \n\t"  // exchange real and imag part

	     "vfmadd213ps 128(%3,%0,4), %%ymm0 , %%ymm12         \n\t"
	     "vfmadd213ps 160(%3,%0,4), %%ymm0 , %%ymm13         \n\t"
	     "vfmadd213ps 192(%3,%0,4), %%ymm0 , %%ymm14         \n\t"
	     "vfmadd213ps 224(%3,%0,4), %%ymm0 , %%ymm15         \n\t"

	     "vfmadd231ps	%%ymm1 , %%ymm4 , %%ymm12  \n\t"
	     "vfmadd231ps	%%ymm1 , %%ymm6 , %%ymm13  \n\t"
	     "vfmadd231ps	%%ymm1 , %%ymm8 , %%ymm14  \n\t"
	     "vfmadd231ps	%%ymm1 , %%ymm10, %%ymm15  \n\t"

	     "vmovups	%%ymm5 ,   (%3,%0,4)		    \n\t"
	     ".p2align 1					    \n\t"
	     "vmovups	%%ymm7 , 32(%3,%0,4)		    \n\t"
	     "vmovups	%%ymm9 , 64(%3,%0,4)		    \n\t"
	     "vmovups	%%ymm11, 96(%3,%0,4)		    \n\t"

	     "vmovups	%%ymm12,128(%3,%0,4)		    \n\t"
	     "vmovups	%%ymm13,160(%3,%0,4)		    \n\t"
	     "vmovups	%%ymm14,192(%3,%0,4)		    \n\t"
	     "vmovups	%%ymm15,224(%3,%0,4)		    \n\t"

	     "addq		$64, %0	  	 	             \n\t"
	     "subq	        $32, %1			             \n\t"		
	     "jnz		1b		             \n\t"
	     "vzeroupper					    \n\t"

	:
          "+r" (i),	// 0	
	  "+r" (n)  	// 1
        :
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (alpha),  // 4
          "r" (mvec)    // 5
	: "cc", 
	  "%xmm0", "%xmm1",
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15", 
	  "memory"
	   );

      }
       else {

            __asm__ __volatile__ (

                  "vzeroupper					    \n\t"
	          "vbroadcastss    (%4), %%ymm0		    \n\t"  // real part of alpha
	          "vbroadcastss    4(%4), %%ymm1		    \n\t"  // imag part of alpha
#if (CONJ) == 1
   	          "vmulps		(%5), %%ymm1 , %%ymm1		    \n\t"
#else
	          "vmulps		(%5), %%ymm0 , %%ymm0		    \n\t"
#endif

	          ".p2align 4				            \n\t"
	          "1:				            \n\t"

	          "vmovaps        (%2,%0,4), %%ymm5                   \n\t" // 4 complex values from x
	          ".p2align 1					    \n\t"
	          "vmovaps      32(%2,%0,4), %%ymm7                   \n\t" // 4 complex values from x
	          "vmovaps      64(%2,%0,4), %%ymm9                   \n\t" // 4 complex values from x
	          "vmovaps      96(%2,%0,4), %%ymm11                  \n\t" // 4 complex values from x

	          "vmovaps     128(%2,%0,4), %%ymm12                  \n\t" // 4 complex values from x
	          "vmovaps     160(%2,%0,4), %%ymm13                  \n\t" // 4 complex values from x
	          "vmovaps     192(%2,%0,4), %%ymm14                  \n\t" // 4 complex values from x
	          "vmovaps     224(%2,%0,4), %%ymm15                  \n\t" // 4 complex values from x

	          "vpermilps	$0xb1 , %%ymm5 , %%ymm4 	    \n\t"  // exchange real and imag part
	          "vpermilps	$0xb1 , %%ymm7 , %%ymm6 	    \n\t"  // exchange real and imag part
	          "vpermilps	$0xb1 , %%ymm9 , %%ymm8 	    \n\t"  // exchange real and imag part
	          "vpermilps	$0xb1 , %%ymm11, %%ymm10 	    \n\t"  // exchange real and imag part

	          "vfmadd213ps    (%3,%0,4), %%ymm0 , %%ymm5          \n\t"
                  ".p2align 1					    \n\t"
	          "vfmadd213ps  32(%3,%0,4), %%ymm0 , %%ymm7          \n\t"
	          "vfmadd213ps  64(%3,%0,4), %%ymm0 , %%ymm9          \n\t"
	          "vfmadd213ps  96(%3,%0,4), %%ymm0 , %%ymm11         \n\t"

	          "vfmadd231ps	%%ymm1 , %%ymm4 , %%ymm5   \n\t"
	          "vfmadd231ps	%%ymm1 , %%ymm6 , %%ymm7   \n\t"
	          "vfmadd231ps	%%ymm1 , %%ymm8 , %%ymm9   \n\t"
	          "vfmadd231ps	%%ymm1 , %%ymm10, %%ymm11  \n\t"

	          "vpermilps	$0xb1 , %%ymm12, %%ymm4 	    \n\t"  // exchange real and imag part
	          "vpermilps	$0xb1 , %%ymm13, %%ymm6 	    \n\t"  // exchange real and imag part
	          "vpermilps	$0xb1 , %%ymm14, %%ymm8 	    \n\t"  // exchange real and imag part
	          "vpermilps	$0xb1 , %%ymm15, %%ymm10 	    \n\t"  // exchange real and imag part

	          "vfmadd213ps 128(%3,%0,4), %%ymm0 , %%ymm12         \n\t"
	          "vfmadd213ps 160(%3,%0,4), %%ymm0 , %%ymm13         \n\t"
	          "vfmadd213ps 192(%3,%0,4), %%ymm0 , %%ymm14         \n\t"
	          "vfmadd213ps 224(%3,%0,4), %%ymm0 , %%ymm15         \n\t"

	          "vfmadd231ps	%%ymm1 , %%ymm4 , %%ymm12  \n\t"
	          "vfmadd231ps	%%ymm1 , %%ymm6 , %%ymm13  \n\t"
	          "vfmadd231ps	%%ymm1 , %%ymm8 , %%ymm14  \n\t"
	          "vfmadd231ps	%%ymm1 , %%ymm10, %%ymm15  \n\t"

	          "vmovaps	%%ymm5 ,   (%3,%0,4)		    \n\t"
	          ".p2align 1					    \n\t"
	          "vmovaps	%%ymm7 , 32(%3,%0,4)		    \n\t"
	          "vmovaps	%%ymm9 , 64(%3,%0,4)		    \n\t"
	          "vmovaps	%%ymm11, 96(%3,%0,4)		    \n\t"

	          "vmovaps	%%ymm12,128(%3,%0,4)		    \n\t"
	          "vmovaps	%%ymm13,160(%3,%0,4)		    \n\t"
	          "vmovaps	%%ymm14,192(%3,%0,4)		    \n\t"
	          "vmovaps	%%ymm15,224(%3,%0,4)		    \n\t"

	          "addq		$64, %0	  	 	             \n\t"
	          "subq	        $32, %1			             \n\t"		
	          "jnz		1b		             \n\t"
	          "vzeroupper					    \n\t"

	          :
                    "+r" (i),	// 0	
	            "+r" (n)  	// 1
                  :
                    "r" (x),      // 2
                    "r" (y),      // 3
                    "r" (alpha),  // 4
                    "r" (mvec)    // 5
	          : "cc", 
	            "%xmm0", "%xmm1",
	            "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	            "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	            "%xmm12", "%xmm13", "%xmm14", "%xmm15", 
	            "memory"


	    );
       }
	   
}



int caxpy(const int32_t n,
          const float da_r,
	  const float da_i,
	  const float * __restrict __ATTR_ALIGN__(32) x,
	  const int32_t incx,
	  float * __restrict __ATTR_ALIGN__(32) y,
	  const int32_t incy) {

            __attribute__((aligned(8)))float da[2] = {};
	    int32_t i  = 0;
            int32_t ix = 0;
	    int32_t iy = 0;
	    if(n <= 0) {
               return (0);
	    }

	    if((incx == 1) && (incy == 1)) {
	    
                int32_t n1 = n & -32;
		if(n1) {
                   da[0] = da_r;
		   da[1] = da_i;
		   caxpy_kernel_8(n1,x,y,da);
		}
		 i = n1;
#if defined __GNUC__ && !defined __INTEL_COMPILER
            x = (float*)__builtin_assume_aligned(x,32);
	    y = (float*)__builtin_assume_aligned(y,32);
#elif defined __ICC || defined __INTEL_COMPILER
            __assume_aligned(x,32); // Probably without effect in case of inline assembly
	    __assume_aligned(y,32);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),max(16),avg(8)
                for(; i < n; ++i) {
#if !defined (CONJ)
                    y[ix]   += (da_r * x[ix]   - da_i * x[ix+1]);
		    y[ix+1] += (da_r * x[ix+1] + da_i * x[ix]);
#else
                    y[ix]   += (da_r * x[ix]   + da_i * x[ix+1]);
		    y[ix+1] -= (da_r * x[ix+1] - da_i * x[ix]);
#endif
                     ix += 2;
		}
		return (0);
	    }
	    incx *= 2;
	    incy *= 2;
#if defined __GNUC__ && !defined __INTEL_COMPILER
            x = (float*)__builtin_assume_aligned(x,32);
	    y = (float*)__builtin_assume_aligned(y,32);
#elif defined __ICC || defined __INTEL_COMPILER
            __assume_aligned(x,32); // Probably without effect in case of inline assembly
	    __assume_aligned(y,32);
#endif
#if defined __ICC || defined __INTEL_COMPILER
#pragma ivdep
#pragma simd vectorlengthfor(4)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma GCC ivdep
#endif
            for(; i != n; ++i) {
#if !defined (CONJ)
                  y[iy]   +=  (da_r * x[ix]   -  da_i * x[ix+1]);
		  y[iy+1] +=  (da_r * x[ix+1] +  da_i * x[ix]);
#else
                  y[iy]   +=  (da_r * x[ix]   -  da_i * x[ix+1]);
		  y[iy+1] -=  (da_r * x[ix+1] +  da_i * x[ix]);
#endif
                  ix += incx;
		  iy += incy;
	    }
	    return (0);
}










