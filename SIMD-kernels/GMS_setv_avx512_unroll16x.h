

#ifndef __GMS_SETV_AVX512_UNROLL16X_H__
#define __GMS_SETV_AVX512_UNROLL16X_H__  190920210954

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

    const unsigned int GMS_SETV_AVX512_UNROLL16X_MAJOR = 1U;
    const unsigned int GMS_SETV_AVX512_UNROLL16X_MINOR = 0U;
    const unsigned int GMS_SETV_AVX512_UNROLL16X_MICRO = 0U;
    const unsigned int GMS_SETV_AVX512_UNROLL16X_FULLVER =
      1000U*GMS_SETV_AVX512_UNROLL16X_MAJOR+
      100U*GMS_SETV_AVX512_UNROLL16X_MINOR+
      10U*GMS_SETV_AVX512_UNROLL16X_MICRO;
    const char * const GMS_SETV512_AVX_UNROLL16X_CREATION_DATE = "19-09-2021 09:54 AM +00200 (SUN 19 SEP 2021 GMT+2)";
    const char * const GMS_SETV512_AVX_UNROLL16X_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_SETV512_AVX_UNROLL16X_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_SETV512_AVX_UNROLL16X_DESCRIPTION   = "AVX512 optimized SETV kernels.";

}

#include <cstdint>
#include "GMS_config.h"


namespace  gms 
{

           namespace math 
           {



	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           void ssetv_u_zmm16r4_unroll16x(const int32_t n,
		                                        const float alpha,
						                                float * __restrict x,
						                                const int32_t incx); 


		
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	         
		         void ssetv_a_zmm16r4_unroll16x(const int32_t n,
		                                  const float alpha,
						                          float * __restrict __ATTR_ALIGN__(64) x,
						                          const int32_t incx); 
						  


		   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	         
		        void dsetv_u_zmm8r8_unroll16x(const int32_t n,
		                                 const double alpha,
						                         double * __restrict x,
						                         const int32_t incx); 
						 


	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	         
		        void dsetv_a_zmm8r8_unroll16x(const int32_t n,
		                                 const double alpha,
						                         double * __restrict __ATTR_ALIGN__(64) x,
						                         const int32_t incx); 
      }

}









#endif /*__GMS_SETV_AVX512_UNROLL16X_H__*/
