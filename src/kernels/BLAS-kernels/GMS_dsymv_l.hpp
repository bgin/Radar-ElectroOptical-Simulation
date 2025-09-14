

#ifndef __GMS_DSYMV_L_HPP__
#define __GMS_DSYMV_L_HPP__

//=========================================================================
// Modified and optimized version of OpenBLAS dsymv_l and kernels kernel_dsymv_l
// Programmer: Bernard Gingold, contact: beniekg@gmail.com
// 28-02-2021 15:58  +00200
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
#if defined __GNUC__ && !defined __INTEL_COMPILER
#include <omp.h>
#endif
#include "GMS_config.h" 

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
__attribute__((no_stack_protector))
static inline
void dsymv_kernel_4x4(const int32_t from,
			const int32_t to,
			double ** __restrict a,
			double *  __restrict x,
			double *  __restrict y,
			double *  __restrict temp1,
			double *  __restrict temp2) {
      	__asm__  __volatile__
	(
	"vzeroupper				     \n\t"
	"vxorpd		%%ymm0 , %%ymm0 , %%ymm0     \n\t"	// temp2[0]
	"vxorpd		%%ymm1 , %%ymm1 , %%ymm1     \n\t"	// temp2[1]
	"vxorpd		%%ymm2 , %%ymm2 , %%ymm2     \n\t"	// temp2[2]
	"vxorpd		%%ymm3 , %%ymm3 , %%ymm3     \n\t"	// temp2[3]
	"vbroadcastsd   (%8),    %%ymm4	             \n\t"	// temp1[0]
	"vbroadcastsd  8(%8),    %%ymm5	             \n\t"	// temp1[1]
	"vbroadcastsd 16(%8),    %%ymm6	             \n\t"	// temp1[1]
	"vbroadcastsd 24(%8),    %%ymm7	             \n\t"	// temp1[1]

	".p2align 4				     \n\t"
	"1:				     \n\t"

	"vmovups	(%3,%0,8), %%ymm9	           \n\t"  // 2 * y
	"vmovups	(%2,%0,8), %%ymm8	           \n\t"  // 2 * x

	"vmovups	(%4,%0,8), %%ymm12	           \n\t"  // 2 * a
	"vmovups	(%5,%0,8), %%ymm13	           \n\t"  // 2 * a
	"vmovups	(%6,%0,8), %%ymm14	           \n\t"  // 2 * a
	"vmovups	(%7,%0,8), %%ymm15	           \n\t"  // 2 * a

	"vfmadd231pd	%%ymm4, %%ymm12 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231pd	%%ymm8, %%ymm12 , %%ymm0  \n\t"  // temp2 += x * a

	"vfmadd231pd	%%ymm5, %%ymm13 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231pd	%%ymm8, %%ymm13 , %%ymm1  \n\t"  // temp2 += x * a

	"vfmadd231pd	%%ymm6, %%ymm14 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231pd	%%ymm8, %%ymm14 , %%ymm2  \n\t"  // temp2 += x * a

	"vfmadd231pd	%%ymm7, %%ymm15 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231pd	%%ymm8, %%ymm15 , %%ymm3  \n\t"  // temp2 += x * a
	"addq		$4 , %0	  	 	      \n\t"

	"vmovups	%%ymm9 ,  -32(%3,%0,8)		   \n\t"

	"cmpq		%0 , %1			      \n\t"
	"jnz		1b		      \n\t"

	"vmovsd		  (%9), %%xmm4		      \n\t"
	"vmovsd		 8(%9), %%xmm5		      \n\t"
	"vmovsd		16(%9), %%xmm6		      \n\t"
	"vmovsd		24(%9), %%xmm7		      \n\t"

	"vextractf128 $0x01, %%ymm0 , %%xmm12	      \n\t"
	"vextractf128 $0x01, %%ymm1 , %%xmm13	      \n\t"
	"vextractf128 $0x01, %%ymm2 , %%xmm14	      \n\t"
	"vextractf128 $0x01, %%ymm3 , %%xmm15	      \n\t"

	"vaddpd	        %%xmm0, %%xmm12, %%xmm0	      \n\t"
	"vaddpd	        %%xmm1, %%xmm13, %%xmm1	      \n\t"
	"vaddpd	        %%xmm2, %%xmm14, %%xmm2	      \n\t"
	"vaddpd	        %%xmm3, %%xmm15, %%xmm3	      \n\t"

	"vhaddpd        %%xmm0, %%xmm0, %%xmm0  \n\t"
	"vhaddpd        %%xmm1, %%xmm1, %%xmm1  \n\t"
	"vhaddpd        %%xmm2, %%xmm2, %%xmm2  \n\t"
	"vhaddpd        %%xmm3, %%xmm3, %%xmm3  \n\t"

	"vaddsd		%%xmm4, %%xmm0, %%xmm0  \n\t"
	"vaddsd		%%xmm5, %%xmm1, %%xmm1  \n\t"
	"vaddsd		%%xmm6, %%xmm2, %%xmm2  \n\t"
	"vaddsd		%%xmm7, %%xmm3, %%xmm3  \n\t"

	"vmovsd         %%xmm0 ,  (%9)		\n\t"	// save temp2
	"vmovsd         %%xmm1 , 8(%9)		\n\t"	// save temp2
	"vmovsd         %%xmm2 ,16(%9)		\n\t"	// save temp2
	"vmovsd         %%xmm3 ,24(%9)		\n\t"	// save temp2
	"vzeroupper				     \n\t"

	:
          "+r" (from)	// 0	
        :
	  "r" (to),  	// 1
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (a[0]),	// 4
          "r" (a[1]),	// 5
          "r" (a[2]),	// 6
          "r" (a[3]),	// 8
          "r" (temp1),  // 8
          "r" (temp2)   // 9
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

}


__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
__attribute__((no_stack_protector))
static inline
int32_t dsymv_l(const int32_t m,
                const int32_t offset,
		const double alpha,
		double * __restrict a,
		const int32_t lda,
		double * __restrict x,
		const int32_t inc_x,
		double * __restrict y,
		const int32_t inc_y,
		double * __restrict buffer) {

      double tmp1[4] = {};
      double tmp2[4] = {};
      double temp1,temp2;
      double *ap[4] = {};
      int32_t i,ix,iy,j,jx,jy;
#if defined __GNUC__ && !defined __INTEL_COMPILER
      a = (float)__builtin_assume_aligned(a,64);
      x = (float)__builtin_assume_aligned(x,64);
      y = (float)__builtin_assume_aligned(y,64);
#elif defined __ICC || defined __INTEL_COMPILER
      __assume_aligned(a,64);
      __assume_aligned(x,64);
      __assume_aligned(y,64);
#endif
      if ( (inc_x != 1) || (inc_y != 1) )
	{

		jx = 0;
		jy = 0;

		for (j=0; j<offset; j++)
		{
			temp1 = alpha * x[jx];
			temp2 = 0.0f;
			y[jy] += temp1 * a[j*lda+j];
			iy = jy;
			ix = jx;
#if defined __ICC || defined __INTEL_COMPILER
#pragma vector aligned
#pragma simd reduction(+:y,temp2)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:y,temp2) aligned(a:64,y:64) 
#endif
			for (i=j+1; i<m; i++)
			{
				ix += inc_x;
				iy += inc_y;
				y[iy] += temp1 * a[j*lda+i];
				temp2 += a[j*lda+i] * x[ix];
			
			}
			y[jy] += alpha * temp2;
			jx    += inc_x;
			jy    += inc_y;
		}
		return(0);
	}

	const int32_t offset1 = (offset/4)*4;

	for (j=0; j<offset1; j+=4)
	{
		tmp1[0] = alpha * x[j];
		tmp1[1] = alpha * x[j+1];
		tmp1[2] = alpha * x[j+2];
		tmp1[3] = alpha * x[j+3];
		tmp2[0] = 0.0;
		tmp2[1] = 0.0;
		tmp2[2] = 0.0;
		tmp2[3] = 0.0;
		ap[0]   = &a[j*lda];
		ap[1]   = ap[0] + lda;
		ap[2]   = ap[1] + lda;
		ap[3]   = ap[2] + lda;
		y[j]   += tmp1[0] * ap[0][j];
		y[j+1] += tmp1[1] * ap[1][j+1];
		y[j+2] += tmp1[2] * ap[2][j+2];
		y[j+3] += tmp1[3] * ap[3][j+3];
		const int32_t from = j+1;
		if ( m - from >=12 )
		{
			const int32_t m2 = (m/4)*4;
			for (i=j+1; i<j+4; i++)
			{
				y[i] += tmp1[0] * ap[0][i];
				tmp2[0] += ap[0][i] * x[i];
			}

			for (i=j+2; i<j+4; i++)
			{
				y[i] += tmp1[1] * ap[1][i];
				tmp2[1] += ap[1][i] * x[i];
			}

			for (i=j+3; i<j+4; i++)
			{
				y[i] += tmp1[2] * ap[2][i];
				tmp2[2] += ap[2][i] * x[i];
			}

			if ( m2 > j+4 )
				dsymv_kernel_4x4(j+4,m2,ap,x,y,tmp1,tmp2);


			for (i=m2; i<m; i++)
			{
				y[i] += tmp1[0] * ap[0][i];
				tmp2[0] += ap[0][i] * x[i];

				y[i] += tmp1[1] * ap[1][i];
				tmp2[1] += ap[1][i] * x[i];

				y[i] += tmp1[2] * ap[2][i];
				tmp2[2] += ap[2][i] * x[i];

				y[i] += tmp1[3] * ap[3][i];
				tmp2[3] += ap[3][i] * x[i];

			}


		}
		else
		{

			for (i=j+1; i<j+4; i++)
			{
				y[i] += tmp1[0] * ap[0][i];
				tmp2[0] += ap[0][i] * x[i];
			}

			for (i=j+2; i<j+4; i++)
			{
				y[i] += tmp1[1] * ap[1][i];
				tmp2[1] += ap[1][i] * x[i];
			}

			for (i=j+3; i<j+4; i++)
			{
				y[i] += tmp1[2] * ap[2][i];
				tmp2[2] += ap[2][i] * x[i];
			}

			for (i=j+4; i<m; i++)
			{
				y[i] += tmp1[0] * ap[0][i];
				tmp2[0] += ap[0][i] * x[i];

				y[i] += tmp1[1] * ap[1][i];
				tmp2[1] += ap[1][i] * x[i];

				y[i] += tmp1[2] * ap[2][i];
				tmp2[2] += ap[2][i] * x[i];

				y[i] += tmp1[3] * ap[3][i];
				tmp2[3] += ap[3][i] * x[i];

			}

		}
		y[j]   += alpha * tmp2[0];
		y[j+1] += alpha * tmp2[1];
		y[j+2] += alpha * tmp2[2];
		y[j+3] += alpha * tmp2[3];
	}


	for (j=offset1; j<offset; j++)
	{
		temp1 = alpha * x[j];
		temp2 = 0.0;
		y[j] += temp1 * a[j*lda+j];
                const int32_t  from = j+1;
		if ( m - from >=8 )
		{
			int32_t j1 = ((from + 4)/4)*4;
			int32_t j2 = (m/4)*4;
#if defined __ICC || defined __INTEL_COMPILER
#pragma vector aligned
#pragma simd reduction(+:y,temp2)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:y,temp2) aligned(a:64,y:64) 
#endif
			for (i=from; i<j1; i++)
			{
				y[i] += temp1 * a[j*lda+i];
				temp2 += a[j*lda+i] * x[i];
			
			}
#if defined __ICC || defined __INTEL_COMPILER
#pragma vector aligned
#pragma simd reduction(+:y,temp2)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:y,temp2) aligned(a:64,y:64) 
#endif
			for (i=j1; i<j2; i++)
			{
				y[i] += temp1 * a[j*lda+i];
				temp2 += a[j*lda+i] * x[i];
			
			}
#if defined __ICC || defined __INTEL_COMPILER
#pragma vector aligned
#pragma simd reduction(+:y,temp2)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:y,temp2) aligned(a:64,y:64) 
#endif
			for (i=j2; i<m; i++)
			{
				y[i] += temp1 * a[j*lda+i];
				temp2 += a[j*lda+i] * x[i];
			
			}

		}
		else
		{
	#if defined __ICC || defined __INTEL_COMPILER
#pragma vector aligned
#pragma simd reduction(+:y,temp2)
#elif defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd reduction(+:y,temp2) aligned(a:64,y:64) 
#endif	
			for (i=from; i<m; i++)
			{
				y[i] += temp1 * a[j*lda+i];
				temp2 += a[j*lda+i] * x[i];
			
			}

		}
		y[j] += alpha * temp2;
	}
	return(0);
}



#endif /*__GMS_DSYMV_L_HPP__*/
