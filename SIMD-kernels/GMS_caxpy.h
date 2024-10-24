

#ifndef __GMS_CAXPY_H__
#define __GMS_CAXPY_H__

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

#include <cstdint>
#include "GMS_config.h" // must define 'CONJ'

// TODO:
// Count the number of uops and implement version for SKX,CLK CPUs.
__ATTR_HOT__
__ATTR_ALIGN__(32)
void caxpy_kernel_8(const int32_t n,
		    float * __restrict __ATTR_ALIGN__(32) x,
		    const  float * __restrict __ATTR_ALIGN__(32) y,
		    float * __restrict  __ATTR_ALIGN__(8) alpha); 

__ATTR_HOT__
__ATTR_ALIGN__(16)
int caxpy(const int32_t n,
          const float da_r,
	  const float da_i,
	  const float * __restrict __ATTR_ALIGN__(32) x,
	  const int32_t incx,
	  float * __restrict __ATTR_ALIGN__(32) y,
	  const int32_t incy); 








#endif /*__GMS_CAXPY_H__*/
