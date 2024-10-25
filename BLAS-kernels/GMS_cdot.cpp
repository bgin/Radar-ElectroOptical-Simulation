



//=========================================================================
// Modified and optimized version of OpenBLAS zscal and kernels zscal_kernel_x
// Programmer: Bernard Gingold, contact: beniekg@gmail.com
// 20-02-2021 12:27  +00200
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


#include "GMS_cdot.h"

 

void cdot_kernel_16(const int32_t n,
		    float * __restrict x,
		    float * __restrict y,
		    float * __restrict dot) {
     bool isAligned32 =  ((uintptr_t)x & 0x1F) != 0ULL && ((uintptr_t)y & 0x1F) != 0ULL;
     int32_t i = 0;
     if(!isAligned32) {

         __asm__ __volatile__ (

	          "vzeroupper					     \n\t"
	          "vxorps		%%ymm0, %%ymm0, %%ymm0	             \n\t"
	          "vxorps		%%ymm1, %%ymm1, %%ymm1	             \n\t"
	          "vxorps		%%ymm2, %%ymm2, %%ymm2	             \n\t"
	          "vxorps		%%ymm3, %%ymm3, %%ymm3	             \n\t"
	          "vxorps		%%ymm4, %%ymm4, %%ymm4	             \n\t"
	          "vxorps		%%ymm5, %%ymm5, %%ymm5	             \n\t"
	          "vxorps		%%ymm6, %%ymm6, %%ymm6	             \n\t"
	          "vxorps		%%ymm7, %%ymm7, %%ymm7	             \n\t"

	          ".p2align 4			             \n\t"
	          "1:				             \n\t"
                  "vmovups                  (%2,%0,4), %%ymm8          \n\t"  // 2 * x
                  "vmovups                32(%2,%0,4), %%ymm9          \n\t"  // 2 * x

                  "vmovups                  (%3,%0,4), %%ymm12         \n\t"  // 2 * y
                  "vmovups                32(%3,%0,4), %%ymm13         \n\t"  // 2 * y

                  "vmovups                64(%2,%0,4), %%ymm10         \n\t"  // 2 * x
                  "vmovups                96(%2,%0,4), %%ymm11         \n\t"  // 2 * x

                  "vmovups                64(%3,%0,4), %%ymm14         \n\t"  // 2 * y
                  "vmovups                96(%3,%0,4), %%ymm15         \n\t"  // 2 * y

	          "vfmadd231ps       %%ymm8 , %%ymm12, %%ymm0     \n\t"  // x_r * y_r, x_i * y_i
	          "vfmadd231ps       %%ymm9 , %%ymm13, %%ymm1     \n\t"  // x_r * y_r, x_i * y_i
	          "vpermilps      $0xb1 , %%ymm12, %%ymm12               \n\t"
	          "vpermilps      $0xb1 , %%ymm13, %%ymm13               \n\t"

	          "vfmadd231ps       %%ymm10, %%ymm14, %%ymm2     \n\t"  // x_r * y_r, x_i * y_i
	          "vfmadd231ps       %%ymm11, %%ymm15, %%ymm3     \n\t"  // x_r * y_r, x_i * y_i
	          "vpermilps      $0xb1 , %%ymm14, %%ymm14               \n\t"
	          "vpermilps      $0xb1 , %%ymm15, %%ymm15               \n\t"

	          "vfmadd231ps       %%ymm8 , %%ymm12, %%ymm4     \n\t"  // x_r * y_i, x_i * y_r
	          "addq		$32 , %0	  	 	             \n\t"
	          "vfmadd231ps       %%ymm9 , %%ymm13, %%ymm5     \n\t"  // x_r * y_i, x_i * y_r
	          "vfmadd231ps       %%ymm10, %%ymm14, %%ymm6     \n\t"  // x_r * y_i, x_i * y_r
	          "subq	        $16 , %1			             \n\t"		
	          "vfmadd231ps       %%ymm11, %%ymm15, %%ymm7     \n\t"  // x_r * y_i, x_i * y_r

	          "jnz		1b		             \n\t"

	          "vaddps        %%ymm0, %%ymm1, %%ymm0	\n\t"
	          "vaddps        %%ymm2, %%ymm3, %%ymm2	\n\t"
	          "vaddps        %%ymm0, %%ymm2, %%ymm0	\n\t"

	          "vaddps        %%ymm4, %%ymm5, %%ymm4	\n\t"
	          "vaddps        %%ymm6, %%ymm7, %%ymm6	\n\t"
	          "vaddps        %%ymm4, %%ymm6, %%ymm4	\n\t"

	          "vextractf128 $1 , %%ymm0 , %%xmm1	\n\t"
	          "vextractf128 $1 , %%ymm4 , %%xmm5	\n\t"

	         "vaddps        %%xmm0, %%xmm1, %%xmm0	\n\t"
	         "vaddps        %%xmm4, %%xmm5, %%xmm4	\n\t"

	         "vmovups       %%xmm0,    (%4)		\n\t"
	         "vmovups       %%xmm4,  16(%4)		\n\t"
	         "vzeroupper					     \n\t"

	         :
                 "+r" (i),	// 0	
	         "+r" (n)  	// 1
                 :
                 "r" (x),      // 2
                 "r" (y),      // 3
                 "r" (dot)     // 4
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

	          "vzeroupper					     \n\t"
	          "vxorps		%%ymm0, %%ymm0, %%ymm0	             \n\t"
	          "vxorps		%%ymm1, %%ymm1, %%ymm1	             \n\t"
	          "vxorps		%%ymm2, %%ymm2, %%ymm2	             \n\t"
	          "vxorps		%%ymm3, %%ymm3, %%ymm3	             \n\t"
	          "vxorps		%%ymm4, %%ymm4, %%ymm4	             \n\t"
	          "vxorps		%%ymm5, %%ymm5, %%ymm5	             \n\t"
	          "vxorps		%%ymm6, %%ymm6, %%ymm6	             \n\t"
	          "vxorps		%%ymm7, %%ymm7, %%ymm7	             \n\t"

	          ".p2align 4			             \n\t"
	          "1:				             \n\t"
                  "vmovaps                  (%2,%0,4), %%ymm8          \n\t"  // 2 * x
                  "vmovaps                32(%2,%0,4), %%ymm9          \n\t"  // 2 * x

                  "vmovaps                  (%3,%0,4), %%ymm12         \n\t"  // 2 * y
                  "vmovaps                32(%3,%0,4), %%ymm13         \n\t"  // 2 * y

                  "vmovaps                64(%2,%0,4), %%ymm10         \n\t"  // 2 * x
                  "vmovaps                96(%2,%0,4), %%ymm11         \n\t"  // 2 * x

                  "vmovaps                64(%3,%0,4), %%ymm14         \n\t"  // 2 * y
                  "vmovaps                96(%3,%0,4), %%ymm15         \n\t"  // 2 * y

	          "vfmadd231ps       %%ymm8 , %%ymm12, %%ymm0     \n\t"  // x_r * y_r, x_i * y_i
	          "vfmadd231ps       %%ymm9 , %%ymm13, %%ymm1     \n\t"  // x_r * y_r, x_i * y_i
	          "vpermilps      $0xb1 , %%ymm12, %%ymm12               \n\t"
	          "vpermilps      $0xb1 , %%ymm13, %%ymm13               \n\t"

	          "vfmadd231ps       %%ymm10, %%ymm14, %%ymm2     \n\t"  // x_r * y_r, x_i * y_i
	          "vfmadd231ps       %%ymm11, %%ymm15, %%ymm3     \n\t"  // x_r * y_r, x_i * y_i
	          "vpermilps      $0xb1 , %%ymm14, %%ymm14               \n\t"
	          "vpermilps      $0xb1 , %%ymm15, %%ymm15               \n\t"

	          "vfmadd231ps       %%ymm8 , %%ymm12, %%ymm4     \n\t"  // x_r * y_i, x_i * y_r
	          "addq		$32 , %0	  	 	             \n\t"
	          "vfmadd231ps       %%ymm9 , %%ymm13, %%ymm5     \n\t"  // x_r * y_i, x_i * y_r
	          "vfmadd231ps       %%ymm10, %%ymm14, %%ymm6     \n\t"  // x_r * y_i, x_i * y_r
	          "subq	        $16 , %1			             \n\t"		
	          "vfmadd231ps       %%ymm11, %%ymm15, %%ymm7     \n\t"  // x_r * y_i, x_i * y_r

	          "jnz		1b		             \n\t"

	          "vaddps        %%ymm0, %%ymm1, %%ymm0	\n\t"
	          "vaddps        %%ymm2, %%ymm3, %%ymm2	\n\t"
	          "vaddps        %%ymm0, %%ymm2, %%ymm0	\n\t"

	          "vaddps        %%ymm4, %%ymm5, %%ymm4	\n\t"
	          "vaddps        %%ymm6, %%ymm7, %%ymm6	\n\t"
	          "vaddps        %%ymm4, %%ymm6, %%ymm4	\n\t"

	          "vextractf128 $1 , %%ymm0 , %%xmm1	\n\t"
	          "vextractf128 $1 , %%ymm4 , %%xmm5	\n\t"

	         "vaddps        %%xmm0, %%xmm1, %%xmm0	\n\t"
	         "vaddps        %%xmm4, %%xmm5, %%xmm4	\n\t"

	         "vmovaps       %%xmm0,    (%4)		\n\t"
	         "vmovaps       %%xmm4,  16(%4)		\n\t"
	         "vzeroupper					     \n\t"

	         :
                 "+r" (i),	// 0	
	         "+r" (n)  	// 1
                 :
                 "r" (x),      // 2
                 "r" (y),      // 3
                 "r" (dot)     // 4
	         : "cc", 
	         "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	         "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	         "%xmm8", "%xmm9", "%xmm10", "%xmm11",
	         "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	          "memory"
	 );
     }
}    
		      



void std::complex<float> cdot(const int32_t n,
                              float * __restrict x,
			      const int32_t inc_x,
			      float * __restrict y,
			      const int32_t inc_y) {
        __ATTR_ALIGN__(32) float dot[8] = {0.0f,0.0f,0.0f,0.0f,
	                                   0.0f,0.0f,0.0f,0.0f};
	int i,ix,iy;
        if((inc_x==1) && (inc_y == 1)) {
            const int32_t n1 = n & -16;
	    if(n1) {
               cdot_kernel_16(n1,x,y,dot)l
	       	dot[0] += dot[2];
		dot[1] += dot[3];
		dot[4] += dot[6];
		dot[5] += dot[7];
	   }
	   i = n1;
	   int32_t j = i+i;
	   while(i < n) {
                 dot[0] += x[j]   * y[j]   ;
		 dot[1] += x[j+1] * y[j+1] ;
		 dot[4] += x[j]   * y[j+1] ;
		 dot[5] += x[j+1] * y[j]   ;

		 j+=2;
		 i++ ;
	   }
	}
	else {

	      	i=0;
		ix=0;
		iy=0;
		inc_x *= 2;
		inc_y *= 2;
		while(i < n) {
		
  		  dot[0] += x[ix]   * y[iy]   ;
		  dot[1] += x[ix+1] * y[iy+1] ;
		  dot[4] += x[ix]   * y[iy+1] ;
		  dot[5] += x[ix+1] * y[iy]   ;
		  ix  += inc_x ;
		  iy  += inc_y ;
		  i++ ;

		}
	}
#if defined(CONJ)
                   std::complex<float> res = (dot[0]-dot[1],dot[4]+dot[5]);
#else
                   std::complex<float> res = (dot[0]+dot[1],dot[4]-dot[5]);
#endif
       return (res);		   
}








