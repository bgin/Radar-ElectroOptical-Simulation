

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
___attribute__((no_stack_protector))
void cscal_kernel_16(const int32_t n,
                     float * __restrict alpha,
		     float * __restrict x); 

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__attribute__((no_stack_protector))
void cscal_kernel_16_zero_r(const int32_t n,
                            float * __restrict alpha,
			    float * __restrict x); 

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__attribute__((no_stack_protector))
void cscal_kernel_16_zero_i(const int32_t n,
                            float * __restrict alpha,
                            float * __restrict x); 
                            
__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__attribute__((no_stack_protector))
void cscal_kernel_16_zero(const int32_t n,
                          float * __restrict alpha,
			  float * __restrict x); 
			  
__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__attribute__((no_stack_protector))
void cscal_kernel_inc_8(const int32_t n,
                        float * __restrict alpha,
			float * __restrict x,
			const int32_t inc_x);
			

__ATTR_HOT__
__ATTR_ALIGN__(32) // default alignment on boundary of 32-bytes (should be tested against ICC or GCC alignment logic)
__attribute__((no_stack_protector))
int32_t  cscal(const int32_t n,
               float da_r,
	       float da_i,
	       float * __restrict x,
	       const int32_t inc_x,
               float * __restrict y,
               const int32_t inc_y); 

	       

#endif /*__GMS_CSCAL_H__*/
