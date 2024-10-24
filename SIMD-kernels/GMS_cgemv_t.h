

#ifndef __GMS_CGEMV_T_H__
#define __GMS_CGEMV_T_H__

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
#include "GMS_config.h" 

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
void cgemv_kernel_4x4(const int32_t n,
		      float ** __restrict ap,
		      float * __restrict  x,
		      float * __restrict  y,
		      float * __restrict  alpha); 

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
void cgemv_kernel_4x2(const int32_t n,
                      float ** __restrict ap,
		      float *  __restrict x,
		      float *  __restrict y,
		      float *  __restrict alpha); 

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
void cgemv_kernel_4x1(const int32_t n,
                      float * __restrict ap,
		      float * __restrict x,
		      float * __restrict y,
		      float * __restrict alpha); 
		      
		      
// default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__ATTR_ALWAYS_INLINE__
static inline
void copy_x(const int32_t n,
	    float * __restrict src,
	    float * __restrict dest,
	    const int32_t inc_src) {

        int32_t i;
        for ( i=0; i<n; i++ )
        {
                *dest     = *src;
                *(dest+1) = *(src+1);
                dest+=2;
                src += inc_src;
        }
}

#include <cstring> // memset
#define CGEMV_T_NBMAX 2048

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
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
	     float * __restrict buffer); 
	     
	     
	     
#endif /*__GMS_CGEMV_T_H__*/
