

#ifndef __GMS_SWAP_AVX_UNROLL8X_H__
#define __GMS_SWAP_AVX_UNROLL8X_H__

/*
     This version is based on BLIS framework original implementation.
     Few changes were introduced in order to more easily adapt it 
     to 'Guided Missile Simulation' project.
     
   The original authors license
   BLIS
   An object-based framework for developing high-performance BLAS-like
   libraries.
   Copyright (C) 2016 - 2019, Advanced Micro Devices, Inc.
   Copyright (C) 2018 - 2020, The University of Texas at Austin. All rights reserved.
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name(s) of the copyright holder(s) nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
     
*/

namespace file_version 
{

    const unsigned int gGMS_SWAP_AVX_UNROLL8X_MAJOR = 1U;
    const unsigned int gGMS_SWAP_AVX_UNROLL8X_MINOR = 0U;
    const unsigned int gGMS_SWAP_AVX_UNROLL8X_MICRO = 0U;
    const unsigned int gGMS_SWAP_AVX_UNROLL8X_FULLVER =
      1000U*gGMS_SWAP_AVX_UNROLL8X_MAJOR+
      100U*gGMS_SWAP_AVX_UNROLL8X_MINOR+
      10U*gGMS_SWAP_AVX_UNROLL8X_MICRO;
    const char * const pgGMS_SWAP_AVX_UNROLL8X_CREATION_DATE = "24-08-2021 14:26 AM +00200 (TUE 24 AUG 2021 GMT+2)";
    const char * const pgGMS_SWAP_AVX_UNROLL8X_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_SWAP_AVX_UNROLL8X_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_SWAP_AVX_UNROLL8X_DESCRIPTION   = "AVX optimized SWAP kernels.";

}

#include <cstdint>
#include "GMS_config.h"

#if !defined(SWAP_AVX_UNROLL8X_SOFT_PREFETCH)
#define SWAP_AVX_UNROLL8X_SOFT_PREFETCH 1 
#endif 

#if !defined(SWAP_AVX_UNROLL8X_INTERLEAVE_SIMD_OPS)
#define SWAP_AVX_UNROLL8X_INTERLEAVE_SIMD_OPS 1
#endif 


namespace gms {

         namespace math {


                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void sswap_u_ymm8r4_unroll8x(const std::size_t n,
		                                      float * __restrict x,
					                                const std::size_t incx,
					                                float * __restrict y,
					                                const std::size_t incy); 




		  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
            void sswap_a_ymm8r4_unroll8x(const std::size_t n,
		                                     float * __restrict __ATTR_ALIGN__(32) x,
					                               const std::size_t incx,
					                               float * __restrict __ATTR_ALIGN__(32) y,
					                               const std::size_t incy); 
					  




		 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
             void dswap_u_ymm4r8_unroll8x(const std::size_t n,
		                                      double * __restrict x,
						                              const std::size_t incx,
						                              double * __restrict y,
						                              const std::size_t incy); 
						

		 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
             void dswap_a_ymm4r8_unroll8x(const std::size_t n,
		                                      double * __restrict __ATTR_ALIGN__(32) x,
						                              const std::size_t incx,
						                              double * __restrict __ATTR_ALIGN__(32) y,
						                              const std::size_t incy); 
						



		 


	 

    } // math

} // gms
































#endif /*__GMS_SWAP_AVX_UNROLL8X_H__*/
