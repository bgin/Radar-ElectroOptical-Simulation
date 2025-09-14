
#ifndef __GMS_SCALV_AVX512_UNROLL8X_HPP__
#define __GMS_SCALV_AVX512_UNROLL8X_HPP__


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

namespace file_version {

    const unsigned int gGMS_SCALV_AVX512_UNROLL8X_MAJOR = 1U;
    const unsigned int gGMS_SCALV_AVX512_UNROLL8X_MINOR = 0U;
    const unsigned int gGMS_SCALV_AVX512_UNROLL8X_MICRO = 0U;
    const unsigned int gGMS_SCALV_AVX512_UNROLL8X_FULLVER =
      1000U*gGMS_SCALV_AVX512_UNROLL8X_MAJOR+
      100U*gGMS_SCALV_AVX512_UNROLL8X_MINOR+
      10U*gGMS_SCALV_AVX512_UNROLL8X_MICRO;
    const char * const pgGMS_SCALV_AVX512_UNROLL8X_CREATION_DATE = "12-09-2021 11:07 AM +00200 (SAT 11 SEP 2021 GMT+2)";
    const char * const pgGMS_SCALV_AVX512_UNROLL8X_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_SCALV_AVX512_UNROLL8X_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_SCALV_AVX512_UNROLL8X_DESCRIPTION   = "AVX512 optimized SCALV kernels."

}

#include <cstdint>
#include <immintrin.h>
#include <omp.h>
#include "GMS_config.h"
#include "GMS_setv_avx512_unroll16x.hpp"

namespace  gms {

        namespace  math {

                namespace  {

		    typedef union {

		        __m512 v;
			__ATTR_ALIGN__(64) float f[16];
		    } zmm16r4_t;

		    typedef union {
                        __m512d v;
			__ATTR_ALIGN__(64) double f[8];
		    } zmm8r8_t;

	     }


	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void sscalv_u_zmm16r4_unroll8x(const int32_t n,
		                                  const float alpha,
						  float * __restrict x,
						  const int32_t incx) {
                        if(__builtin_expect(0==n,0) ||
			   __builtin_expect(1.0f==alpha,0)) {
                           return;
			}
			if(__builtin_expect(0.0f==alpha,0)) {
                           // call setv here!!
			  ssetv_u_zmm16r4_unroll16x(n,
			                            alpha,
						    x,incx);
			  return;
			}
			__ATTR_ALIGN__(64) zmm16r4_t xv[8];
			zmm16r4_t alphav;
			int32_t i;

			if(__builtin_expect(1==incx,1)) {

                             alphav.v = _mm512_broadcast_ss(alpha);
			     for(i = 0; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[0].v = _mm512_loadu_ps(&x[i+0]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  _mm512_storeu_ps(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_loadu_ps(&x[i+16]);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  _mm512_storeu_ps(&x[i+16],xv[1].v);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv[2].v = _mm512_loadu_ps(&x[i+32]);
				  xv[2].v = _mm512_mul_ps(alphav.v,xv[2].v);
				  _mm512_storeu_ps(&x[i+32],xv[2].v);
				  xv[3].v = _mm512_loadu_ps(&x[i+48]);
				  xv[3].v = _mm512_mul_ps(alphav.v,xv[3].v);
				  _mm512_storeu_ps(&x[i+48],xv[3].v);
				  _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				  xv[4].v = _mm512_loadu_ps(&x[i+64]);
				  xv[4].v = _mm512_mul_ps(alphav.v,xv[4].v);
				  _mm512_storeu_ps(&x[i+64],xv[4].v);
				  xv[5].v = _mm512_loadu_ps(&x[i+80]);
				  xv[5].v = _mm512_mul_ps(alphav.v,xv[5].v);
				  _mm512_storeu_ps(&x[i+80],xv[5].v);
				  _mm_prefetch((const char*)&x[i+128],_MM_HINT_T0);
				  xv[6].v = _mm512_loadu_ps(&x[i+96]);
				  xv[6].v = _mm512_mul_ps(alphav.v,xv[6].v);
				  _mm512_storeu_ps(&x[i+96],xv[6].v);
				  xv[7].v = _mm512_loadu_ps(&x[i+112]);
				  xv[7].v = _mm512_mul_ps(alphav.v,xv[7].v);
				  _mm512_storeu_ps(&x[i+112],xv[7].v);
#else
                                  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+128],_MM_HINT_T0);
				  xv[0].v = _mm512_loadu_ps(&x[i+0]);
				  xv[1].v = _mm512_loadu_ps(&x[i+16]);
				  xv[2].v = _mm512_loadu_ps(&x[i+32]);
				  xv[3].v = _mm512_loadu_ps(&x[i+48]);
				  xv[4].v = _mm512_loadu_ps(&x[i+64]);
				  xv[5].v = _mm512_loadu_ps(&x[i+80]);
				  xv[6].v = _mm512_loadu_ps(&x[i+96]);
				  xv[7].v = _mm512_loadu_ps(&x[i+112]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  xv[2].v = _mm512_mul_ps(alphav.v,xv[2].v);
				  xv[3].v = _mm512_mul_ps(alphav.v,xv[3].v);
				  xv[4].v = _mm512_mul_ps(alphav.v,xv[4].v);
				  xv[5].v = _mm512_mul_ps(alphav.v,xv[5].v);
				  xv[6].v = _mm512_mul_ps(alphav.v,xv[6].v);
				  xv[7].v = _mm512_mul_ps(alphav.v,xv[7].v);
				  _mm512_storeu_ps(&x[i+0],xv[0].v);
				  _mm512_storeu_ps(&x[i+16],xv[1].v);
				  _mm512_storeu_ps(&x[i+32],xv[2].v);
				  _mm512_storeu_ps(&x[i+48],xv[3].v);
				  _mm512_storeu_ps(&x[i+64],xv[4].v);
				  _mm512_storeu_ps(&x[i+80],xv[5].v);
				  _mm512_storeu_ps(&x[i+96],xv[6].v);
				  _mm512_storeu_ps(&x[i+112],xv[7].v);
				  
#endif
			     }

			     for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_loadu_ps(&x[i+0]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  _mm512_storeu_ps(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_loadu_ps(&x[i+16]);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  _mm512_storeu_ps(&x[i+16],xv[1].v);
				  xv[2].v = _mm512_loadu_ps(&x[i+32]);
				  xv[2].v = _mm512_mul_ps(alphav.v,xv[2].v);
				  _mm512_storeu_ps(&x[i+32],xv[2].v);
				  xv[3].v = _mm512_loadu_ps(&x[i+48]);
				  xv[3].v = _mm512_mul_ps(alphav.v,xv[3].v);
#else
                                  xv[0].v = _mm512_loadu_ps(&x[i+0]);
				  xv[1].v = _mm512_loadu_ps(&x[i+16]);
				  xv[2].v = _mm512_loadu_ps(&x[i+32]);
				  xv[3].v = _mm512_loadu_ps(&x[i+48]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  xv[2].v = _mm512_mul_ps(alphav.v,xv[2].v);
				  xv[3].v = _mm512_mul_ps(alphav.v,xv[3].v);
				  _mm512_storeu_ps(&x[i+0],xv[0].v);
				  _mm512_storeu_ps(&x[i+16],xv[1].v);
				  _mm512_storeu_ps(&x[i+32],xv[2].v);
				  _mm512_storeu_ps(&x[i+48],xv[3].v);
#endif
			     }

			     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_loadu_ps(&x[i+0]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  _mm512_storeu_ps(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_loadu_ps(&x[i+16]);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  _mm512_storeu_ps(&x[i+16],xv[1].v);
#else
                                  xv[0].v = _mm512_loadu_ps(&x[i+0]);
				  xv[1].v = _mm512_loadu_ps(&x[i+16]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  _mm512_storeu_ps(&x[i+0],xv[0].v);
				  _mm512_storeu_ps(&x[i+16],xv[1].v);
#endif
			     }

			     for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_loadu_ps(&x[i+0]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  _mm512_storeu_ps(&x[i+0],xv[0].v);
#else
                                   xv[0].v = _mm512_loadu_ps(&x[i+0]);
				   xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				   _mm512_storeu_ps(&x[i+0],xv[0].v); 
#endif
			     }

			     for(; (i+0) < n; i += 1) {
                                 x[i] *= alpha;
			      }
			}
			 else {
                                   for(i = 0; i != n; ++i) {
                                        *x *= alpha;
				         x += incx;
			           }

			 }

		  }



		   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void sscalv_a_zmm16r4_unroll8x(const int32_t n,
		                                  const float alpha,
						  float * __restrict __ATTR_ALIGN__(64) x,
						  const int32_t incx) {
                        if(__builtin_expect(0==n,0) ||
			   __builtin_expect(1.0f==alpha,0)) {
                           return;
			}
			if(__builtin_expect(0.0f==alpha,0)) {
                           // call setv here!!
			    ssetv_a_zmm16r4_unroll16x(n,
			                            alpha,
						    x,incx);
			  return;
			}
			__ATTR_ALIGN__(64) zmm16r4_t xv[8];
			zmm16r4_t alphav;
			int32_t i;

			if(__builtin_expect(1==incx,1)) {

                             alphav.v = _mm512_broadcast_ss(alpha);
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (float*)__builtin_assume_aligned(x,64);
#endif
			     for(i = 0; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[0].v = _mm512_load_ps(&x[i+0]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  _mm512_store_ps(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_load_ps(&x[i+16]);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  _mm512_store_ps(&x[i+16],xv[1].v);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv[2].v = _mm512_load_ps(&x[i+32]);
				  xv[2].v = _mm512_mul_ps(alphav.v,xv[2].v);
				  _mm512_store_ps(&x[i+32],xv[2].v);
				  xv[3].v = _mm512_load_ps(&x[i+48]);
				  xv[3].v = _mm512_mul_ps(alphav.v,xv[3].v);
				  _mm512_store_ps(&x[i+48],xv[3].v);
				  _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				  xv[4].v = _mm512_load_ps(&x[i+64]);
				  xv[4].v = _mm512_mul_ps(alphav.v,xv[4].v);
				  _mm512_store_ps(&x[i+64],xv[4].v);
				  xv[5].v = _mm512_load_ps(&x[i+80]);
				  xv[5].v = _mm512_mul_ps(alphav.v,xv[5].v);
				  _mm512_store_ps(&x[i+80],xv[5].v);
				  _mm_prefetch((const char*)&x[i+128],_MM_HINT_T0);
				  xv[6].v = _mm512_load_ps(&x[i+96]);
				  xv[6].v = _mm512_mul_ps(alphav.v,xv[6].v);
				  _mm512_store_ps(&x[i+96],xv[6].v);
				  xv[7].v = _mm512_load_ps(&x[i+112]);
				  xv[7].v = _mm512_mul_ps(alphav.v,xv[7].v);
				  _mm512_store_ps(&x[i+112],xv[7].v);
#else
                                  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+128],_MM_HINT_T0);
				  xv[0].v = _mm512_load_ps(&x[i+0]);
				  xv[1].v = _mm512_load_ps(&x[i+16]);
				  xv[2].v = _mm512_load_ps(&x[i+32]);
				  xv[3].v = _mm512_load_ps(&x[i+48]);
				  xv[4].v = _mm512_load_ps(&x[i+64]);
				  xv[5].v = _mm512_load_ps(&x[i+80]);
				  xv[6].v = _mm512_load_ps(&x[i+96]);
				  xv[7].v = _mm512_load_ps(&x[i+112]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  xv[2].v = _mm512_mul_ps(alphav.v,xv[2].v);
				  xv[3].v = _mm512_mul_ps(alphav.v,xv[3].v);
				  xv[4].v = _mm512_mul_ps(alphav.v,xv[4].v);
				  xv[5].v = _mm512_mul_ps(alphav.v,xv[5].v);
				  xv[6].v = _mm512_mul_ps(alphav.v,xv[6].v);
				  xv[7].v = _mm512_mul_ps(alphav.v,xv[7].v);
				  _mm512_store_ps(&x[i+0],xv[0].v);
				  _mm512_store_ps(&x[i+16],xv[1].v);
				  _mm512_store_ps(&x[i+32],xv[2].v);
				  _mm512_store_ps(&x[i+48],xv[3].v);
				  _mm512_store_ps(&x[i+64],xv[4].v);
				  _mm512_store_ps(&x[i+80],xv[5].v);
				  _mm512_store_ps(&x[i+96],xv[6].v);
				  _mm512_store_ps(&x[i+112],xv[7].v);
				  
#endif
			     }

			     for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_load_ps(&x[i+0]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  _mm512_store_ps(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_load_ps(&x[i+16]);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  _mm512_store_ps(&x[i+16],xv[1].v);
				  xv[2].v = _mm512_load_ps(&x[i+32]);
				  xv[2].v = _mm512_mul_ps(alphav.v,xv[2].v);
				  _mm512_store_ps(&x[i+32],xv[2].v);
				  xv[3].v = _mm512_load_ps(&x[i+48]);
				  xv[3].v = _mm512_mul_ps(alphav.v,xv[3].v);
				  _mm512_store_ps(&x[i+48],xv[3].v);
#else
                                  xv[0].v = _mm512_load_ps(&x[i+0]);
				  xv[1].v = _mm512_load_ps(&x[i+16]);
				  xv[2].v = _mm512_load_ps(&x[i+32]);
				  xv[3].v = _mm512_load_ps(&x[i+48]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  xv[2].v = _mm512_mul_ps(alphav.v,xv[2].v);
				  xv[3].v = _mm512_mul_ps(alphav.v,xv[3].v);
				  _mm512_store_ps(&x[i+0],xv[0].v);
				  _mm512_store_ps(&x[i+16],xv[1].v);
				  _mm512_store_ps(&x[i+32],xv[2].v);
				  _mm512_store_ps(&x[i+48],xv[3].v);
#endif
			     }

			     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_load_ps(&x[i+0]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  _mm512_store_ps(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_load_ps(&x[i+16]);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  _mm512_store_ps(&x[i+16],xv[1].v);
#else
                                  xv[0].v = _mm512_load_ps(&x[i+0]);
				  xv[1].v = _mm512_load_ps(&x[i+16]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_ps(alphav.v,xv[1].v);
				  _mm512_store_ps(&x[i+0],xv[0].v);
				  _mm512_store_ps(&x[i+16],xv[1].v);
#endif
			     }

			     for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_load_ps(&x[i+0]);
				  xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				  _mm512_store_ps(&x[i+0],xv[0].v);
#else
                                   xv[0].v = _mm512_load_ps(&x[i+0]);
				   xv[0].v = _mm512_mul_ps(alphav.v,xv[0].v);
				   _mm512_store_ps(&x[i+0],xv[0].v); 
#endif
			     }

			     for(; (i+0) < n; i += 1) {
                                 x[i] *= alpha;
			        }
			}
			 else {
                                   for(i = 0; i != n; ++i) {
                                        *x *= alpha;
				         x += incx;
			           }

			 }

		  }



		   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void sscalv_a_zmm16r4_unroll8x_omp(const int32_t n,
		                                      const float alpha,
						      float * __restrict __ATTR_ALIGN__(64) x,
						      const int32_t incx) {
                        if(__builtin_expect(0==n,0) ||
			   __builtin_expect(1.0f==alpha,0)) {
                           return;
			}
			if(__builtin_expect(0.0f==alpha,0)) {
                           // call setv here!!
			    ssetv_a_zmm16r4_unroll16x(n,
			                            alpha,
						    x,incx);
			  return;
			}
		        zmm16r4_t xv0;
			zmm16r4_t xv1;
			zmm16r4_t xv2;
			zmm16r4_t xv3;
			zmm16r4_t xv4;
			zmm16r4_t xv5;
			zmm16r4_t xv6;
			zmm16r4_t xv7;
			zmm16r4_t alphav;
			int32_t i;
                        int32_t ii;
			ii = 0;
			if(__builtin_expect(1==incx,1)) {

                             alphav.v = _mm512_broadcast_ss(alpha);
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (float*)__builtin_assume_aligned(x,64);
#endif
#pragma omp parallel for schedule(static,128) default(none) private(i)  \
              shared(x,n) lastprivate(ii,xv0,xv1,xv2,xv3,xv4,xv5,xv6,xv7)
			     for(i = 0; (i+127) < n; i += 128) {
                                   ii = i;
                                  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv0.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[i+0]),
						    alphav.v);
				  _mm512_store_ps(&x[i+0],xv0.v);
				  xv1.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[i+16]),
						    alphav.v);
				  _mm512_store_ps(&x[i+16],xv1.v);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv2.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[i+32]),
						    alphav.v);
				  _mm512_store_ps(&x[i+32],xv2.v);
				  xv3.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[i+48]),
						    alphav.v);
				  _mm512_store_ps(&x[i+48],xv3.v);
				  _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				  xv4.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[i+64]),
						    alphav.v);
				  _mm512_store_ps(&x[i+64],xv4.v);
				  xv5.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[i+80]),
						    alphav.v);
				  _mm512_store_ps(&x[i+80],xv5.v);
				  _mm_prefetch((const char*)&x[i+128],_MM_HINT_T0);
				  xv6.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[i+96]),
						    alphav.v);
				  _mm512_store_ps(&x[i+96],xv6.v);
				  xv7.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[i+112]),
						    alphav.v);
				  _mm512_store_ps(&x[i+112],xv7.v);
				 
			     }

			     for(; (ii+63) < n; ii += 64) {

			         xv0.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[ii+0]),
						    alphav.v);
				  _mm512_store_ps(&x[ii+0],xv0.v);
				  xv1.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[ii+16]),
						    alphav.v);
				  _mm512_store_ps(&x[ii+16],xv1.v);
				 
				  xv2.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[ii+32]),
						    alphav.v);
				  _mm512_store_ps(&x[ii+32],xv2.v);
				  xv3.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[ii+48]),
						    alphav.v);
				  _mm512_store_ps(&x[ii+48],xv3.v);
                                 

			     }

			     for(; (ii+31) < n; ii += 32) {

			         xv0.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[ii+0]),
						    alphav.v);
				  _mm512_store_ps(&x[ii+0],xv0.v);
				  xv1.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[ii+16]),
						    alphav.v);
				  _mm512_store_ps(&x[ii+16],xv1.v);
                                 

			     }

			     for(; (ii+15) < n; ii += 16) {

                                  xv0.v = _mm512_mul_ps(
				                    _mm512_load_ps(&x[ii+0]),
						    alphav.v);
				  _mm512_store_ps(&x[ii+0],xv0.v);

			     }

			     for(; (ii+0) < n; ii += 1) {
                                 x[ii] *= alpha;
			      }
			}
			 else {
                                   for(i = 0; i != n; ++i) {
                                        *x *= alpha;
				         x += incx;
			           }

			 }

		  }



		   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void dscalv_u_zmm8r8_unroll8x(const int32_t n,
		                                 const double alpha,
						 double * __restrict x,
						 const int32_t incx) {
                        if(__builtin_expect(0==n,0) ||
			   __builtin_expect(1.0==alpha,0)) {
                           return;
			}
			if(__builtin_expect(0.0==alpha,0)) {
                           // call setv here!!
			    dsetv_u_zmm8r8_unroll16x(n,
			                            alpha,
						    x,incx);
			    return;
			}
			__ATTR_ALIGN__(64) zmm16r4_t xv[8];
			zmm8r8_t alphav;
			int32_t i;
			if(__builtin_expect(1==incx,1)) {
                             alphav.v = _mm512_broadcast_sd(alpha);
                             for(i = 0; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				  xv[0].v = _mm512_loadu_pd(&x[i+0]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  _mm512_storeu_pd(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_loadu_pd(&x[i+8]);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  _mm512_storeu_pd(&x[i+8],xv[1].v);
				  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[2].v = _mm512_loadu_pd(&x[i+16]);
				  xv[2].v = _mm512_mul_pd(alphav.v,xv[2].v);
				  _mm512_storeu_pd(&x[i+16],xv[2].v);
				  xv[3].v = _mm512_loadu_pd(&x[i+24]);
				  xv[3].v = _mm512_mul_pd(alphav.v,xv[3].v);
				  _mm512_storeu_pd(&x[i+24],xv[3].v);
				  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				  xv[4].v = _mm512_loadu_pd(&x[i+32]);
				  xv[4].v = _mm512_mul_pd(alphav.v,xv[4].v);
				  _mm512_storeu_pd(&x[i+32],xv[4].v);
				  xv[5].v = _mm512_loadu_pd(&x[i+40]);
				  xv[5].v = _mm512_mul_pd(alphav.v,xv[5].v);
				  _mm512_storeu_pd(&x[i+40],xv[5].v);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv[6].v = _mm512_loadu_pd(&x[i+48]);
				  xv[6].v = _mm512_mul_pd(alphav.v,xv[6].v);
				  _mm512_storeu_pd(&x[i+48],xv[6].v);
				  xv[7].v = _mm512_loadu_pd(&x[i+56]);
				  xv[7].v = _mm512_mul_pd(alphav.v,xv[7].v);
				  _mm512_storeu_pd(&x[i+56],xv[7].v);
#else
                                  _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv[0].v = _mm512_loadu_pd(&x[i+0]);
				  xv[1].v = _mm512_loadu_pd(&x[i+8]);
				  xv[2].v = _mm512_loadu_pd(&x[i+16]);
				  xv[3].v = _mm512_loadu_pd(&x[i+24]);
				  xv[4].v = _mm512_loadu_pd(&x[i+32]);
				  xv[5].v = _mm512_loadu_pd(&x[i+40]);
				  xv[6].v = _mm512_loadu_pd(&x[i+48]);
				  xv[7].v = _mm512_loadu_pd(&x[i+56]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  xv[2].v = _mm512_mul_pd(alphav.v,xv[2].v);
				  xv[3].v = _mm512_mul_pd(alphav.v,xv[3].v);
				  xv[4].v = _mm512_mul_pd(alphav.v,xv[4].v);
				  xv[5].v = _mm512_mul_pd(alphav.v,xv[5].v);
				  xv[6].v = _mm512_mul_pd(alphav.v,xv[6].v);
				  xv[7].v = _mm512_mul_pd(alphav.v,xv[7].v);
				  _mm512_storeu_pd(&x[i+0],xv[0].v);
				  _mm512_storeu_pd(&x[i+8],xv[1].v);
				  _mm512_storeu_pd(&x[i+16],xv[2].v);
				  _mm512_storeu_pd(&x[i+24],xv[3].v);
				  _mm512_storeu_pd(&x[i+32],xv[4].v);
				  _mm512_storeu_pd(&x[i+40],xv[5].v);
				  _mm512_storeu_pd(&x[i+48],xv[6].v);
				  _mm512_storeu_pd(&x[i+56],xv[7].v);
#endif
			     }

			     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_loadu_pd(&x[i+0]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  _mm512_storeu_pd(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_loadu_pd(&x[i+8]);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  _mm512_storeu_pd(&x[i+8],xv[1].v);
				  //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[2].v = _mm512_loadu_pd(&x[i+16]);
				  xv[2].v = _mm512_mul_pd(alphav.v,xv[2].v);
				  _mm512_storeu_pd(&x[i+16],xv[2].v);
				  xv[3].v = _mm512_loadu_pd(&x[i+24]);
				  xv[3].v = _mm512_mul_pd(alphav.v,xv[3].v);
				  _mm512_storeu_pd(&x[i+24],xv[3].v);
#else
                                  xv[0].v = _mm512_loadu_pd(&x[i+0]);
				  xv[1].v = _mm512_loadu_pd(&x[i+8]);
				  xv[2].v = _mm512_loadu_pd(&x[i+16]);
				  xv[3].v = _mm512_loadu_pd(&x[i+24]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  xv[2].v = _mm512_mul_pd(alphav.v,xv[2].v);
				  xv[3].v = _mm512_mul_pd(alphav.v,xv[3].v);
				  _mm512_storeu_pd(&x[i+0],xv[0].v);
				  _mm512_storeu_pd(&x[i+8],xv[1].v);
				  _mm512_storeu_pd(&x[i+16],xv[2].v);
				  _mm512_storeu_pd(&x[i+24],xv[3].v);
#endif
			     }

			     for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_loadu_pd(&x[i+0]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  _mm512_storeu_pd(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_loadu_pd(&x[i+8]);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  _mm512_storeu_pd(&x[i+8],xv[1].v);
#else
                                  xv[0].v = _mm512_loadu_pd(&x[i+0]);
				  xv[1].v = _mm512_loadu_pd(&x[i+8]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  _mm512_storeu_pd(&x[i+0],xv[0].v);
				  _mm512_storeu_pd(&x[i+8],xv[1].v);
#endif

			     }

			     for(; (i+7) < n; i += 8) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_loadu_pd(&x[i+0]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  _mm512_storeu_pd(&x[i+0],xv[0].v);
#else
                                  xv[0].v = _mm512_loadu_pd(&x[i+0]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  _mm512_storeu_pd(&x[i+0],xv[0].v);
#endif

			     }

			     for(; (i+0) < n; i += 1) {
                                   x[i] *= alpha;
			     }
			}
			else {

			      for(i = 0; i != n; ++i) {
                                   *x *= alpha;
				    x += incx;
			       }
			}

		  }



		   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void dscalv_a_zmm8r8_unroll8x(const int32_t n,
		                                 const double alpha,
						 double * __restrict __ATTR_ALIGN__(64) x,
						 const int32_t incx) {
                        if(__builtin_expect(0==n,0) ||
			   __builtin_expect(1.0==alpha,0)) {
                           return;
			}
			if(__builtin_expect(0.0==alpha,0)) {
                           // call setv here!!
			     dsetv_a_zmm8r8_unroll16x(n,
			                            alpha,
						    x,incx);
			     return;
			}
			__ATTR_ALIGN__(64) zmm16r4_t xv[8];
			zmm8r8_t alphav;
			int32_t i;
			if(__builtin_expect(1==incx,1)) {
                             alphav.v = _mm512_broadcast_sd(alpha);
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (double*)__builtin_assume_aligned(x,64);
#endif
                             for(i = 0; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				  xv[0].v = _mm512_load_pd(&x[i+0]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  _mm512_store_pd(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_load_pd(&x[i+8]);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  _mm512_store_pd(&x[i+8],xv[1].v);
				  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[2].v = _mm512_load_pd(&x[i+16]);
				  xv[2].v = _mm512_mul_pd(alphav.v,xv[2].v);
				  _mm512_store_pd(&x[i+16],xv[2].v);
				  xv[3].v = _mm512_load_pd(&x[i+24]);
				  xv[3].v = _mm512_mul_pd(alphav.v,xv[3].v);
				  _mm512_store_pd(&x[i+24],xv[3].v);
				  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				  xv[4].v = _mm512_load_pd(&x[i+32]);
				  xv[4].v = _mm512_mul_pd(alphav.v,xv[4].v);
				  _mm512_store_pd(&x[i+32],xv[4].v);
				  xv[5].v = _mm512_load_pd(&x[i+40]);
				  xv[5].v = _mm512_mul_pd(alphav.v,xv[5].v);
				  _mm512_store_pd(&x[i+40],xv[5].v);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv[6].v = _mm512_load_pd(&x[i+48]);
				  xv[6].v = _mm512_mul_pd(alphav.v,xv[6].v);
				  _mm512_store_pd(&x[i+48],xv[6].v);
				  xv[7].v = _mm512_load_pd(&x[i+56]);
				  xv[7].v = _mm512_mul_pd(alphav.v,xv[7].v);
				  _mm512_store_pd(&x[i+56],xv[7].v);
#else
                                  _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv[0].v = _mm512_load_pd(&x[i+0]);
				  xv[1].v = _mm512_load_pd(&x[i+8]);
				  xv[2].v = _mm512_load_pd(&x[i+16]);
				  xv[3].v = _mm512_load_pd(&x[i+24]);
				  xv[4].v = _mm512_load_pd(&x[i+32]);
				  xv[5].v = _mm512_load_pd(&x[i+40]);
				  xv[6].v = _mm512_load_pd(&x[i+48]);
				  xv[7].v = _mm512_load_pd(&x[i+56]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  xv[2].v = _mm512_mul_pd(alphav.v,xv[2].v);
				  xv[3].v = _mm512_mul_pd(alphav.v,xv[3].v);
				  xv[4].v = _mm512_mul_pd(alphav.v,xv[4].v);
				  xv[5].v = _mm512_mul_pd(alphav.v,xv[5].v);
				  xv[6].v = _mm512_mul_pd(alphav.v,xv[6].v);
				  xv[7].v = _mm512_mul_pd(alphav.v,xv[7].v);
				  _mm512_store_pd(&x[i+0],xv[0].v);
				  _mm512_store_pd(&x[i+8],xv[1].v);
				  _mm512_store_pd(&x[i+16],xv[2].v);
				  _mm512_store_pd(&x[i+24],xv[3].v);
				  _mm512_store_pd(&x[i+32],xv[4].v);
				  _mm512_store_pd(&x[i+40],xv[5].v);
				  _mm512_store_pd(&x[i+48],xv[6].v);
				  _mm512_store_pd(&x[i+56],xv[7].v);
#endif
			     }

			     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_load_pd(&x[i+0]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  _mm512_store_pd(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_load_pd(&x[i+8]);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  _mm512_store_pd(&x[i+8],xv[1].v);
				  //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv[2].v = _mm512_load_pd(&x[i+16]);
				  xv[2].v = _mm512_mul_pd(alphav.v,xv[2].v);
				  _mm512_store_pd(&x[i+16],xv[2].v);
				  xv[3].v = _mm512_load_pd(&x[i+24]);
				  xv[3].v = _mm512_mul_pd(alphav.v,xv[3].v);
				  _mm512_store_pd(&x[i+24],xv[3].v);
#else
                                  xv[0].v = _mm512_load_pd(&x[i+0]);
				  xv[1].v = _mm512_load_pd(&x[i+8]);
				  xv[2].v = _mm512_load_pd(&x[i+16]);
				  xv[3].v = _mm512_load_pd(&x[i+24]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  xv[2].v = _mm512_mul_pd(alphav.v,xv[2].v);
				  xv[3].v = _mm512_mul_pd(alphav.v,xv[3].v);
				  _mm512_store_pd(&x[i+0],xv[0].v);
				  _mm512_store_pd(&x[i+8],xv[1].v);
				  _mm512_store_pd(&x[i+16],xv[2].v);
				  _mm512_store_pd(&x[i+24],xv[3].v);
#endif
			     }

			     for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_load_pd(&x[i+0]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  _mm512_store_pd(&x[i+0],xv[0].v);
				  xv[1].v = _mm512_load_pd(&x[i+8]);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  _mm512_store_pd(&x[i+8],xv[1].v);
#else
                                  xv[0].v = _mm512_load_pd(&x[i+0]);
				  xv[1].v = _mm512_load_pd(&x[i+8]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  xv[1].v = _mm512_mul_pd(alphav.v,xv[1].v);
				  _mm512_store_pd(&x[i+0],xv[0].v);
				  _mm512_store_pd(&x[i+8],xv[1].v);
#endif

			     }

			     for(; (i+7) < n; i += 8) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  xv[0].v = _mm512_load_pd(&x[i+0]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  _mm512_store_pd(&x[i+0],xv[0].v);
#else
                                  xv[0].v = _mm512_load_pd(&x[i+0]);
				  xv[0].v = _mm512_mul_pd(alphav.v,xv[0].v);
				  _mm512_store_pd(&x[i+0],xv[0].v);
#endif

			     }

			     for(; (i+0) < n; i += 1) {
                                   x[i] *= alpha;
			     }
			}
			else {

			      for(i = 0; i != n; ++i) {
                                   *x *= alpha;
				    x += incx;
			       }
			}

		  }



		   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void dscalv_a_zmm16r4_unroll8x_omp(const int32_t n,
		                                      const double alpha,
						      double * __restrict __ATTR_ALIGN__(64) x,
						      const int32_t incx) {
                        if(__builtin_expect(0==n,0) ||
			   __builtin_expect(1.0==alpha,0)) {
                           return;
			}
			if(__builtin_expect(0.0==alpha,0)) {
                           // call setv here!!
			    dsetv_a_zmm8r8_unroll16x(n,
			                            alpha,
						    x,incx);
			    return;
			}
		        zmm8r8_t xv0;
			zmm8r8_t xv1;
			zmm8r8_t xv2;
			zmm8r8_t xv3;
			zmm8r8_t xv4;
			zmm8r8_t xv5;
			zmm8r8_t xv6;
			zmm8r8_t xv7;
			zmm8r8_t alphav;
			int32_t i;
                        int32_t ii;
			ii = 0;
			if(__builtin_expect(1==incx,1)) {

                             alphav.v = _mm512_broadcast_sd(alpha);
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (dobule*)__builtin_assume_aligned(x,64);
#endif
#pragma omp parallel for schedule(static,64) default(none) private(i)  \
              shared(x,n) lastprivate(ii,xv0,xv1,xv2,xv3,xv4,xv5,xv6,xv7)
			     for(i = 0; (i+63) < n; i += 64) {
                                   ii = i;
                                  _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				  xv0.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[i+0]),
						    alphav.v);
				  _mm512_store_pd(&x[i+0],xv0.v);
				  xv1.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[i+8]),
						    alphav.v);
				  _mm512_store_pd(&x[i+8],xv1.v);
				  _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv2.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[i+16]),
						    alphav.v);
				  _mm512_store_pd(&x[i+16],xv2.v);
				  xv3.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[i+24]),
						    alphav.v);
				  _mm512_store_pd(&x[i+24],xv3.v);
				  _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				  xv4.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[i+32]),
						    alphav.v);
				  _mm512_store_pd(&x[i+32],xv4.v);
				  xv5.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[i+40]),
						    alphav.v);
				  _mm512_store_pd(&x[i+40],xv5.v);
				  _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				  xv6.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[i+48]),
						    alphav.v);
				  _mm512_store_pd(&x[i+48],xv6.v);
				  xv7.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[i+56]),
						    alphav.v);
				  _mm512_store_ps(&x[i+56],xv7.v);
				 
			     }

			     for(; (ii+31) < n; ii += 32) {

			          xv0.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[ii+0]),
						    alphav.v);
				  _mm512_store_pd(&x[ii+0],xv0.v);
				  xv1.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[ii+8]),
						    alphav.v);
				  _mm512_store_pd(&x[ii+8],xv1.v);
				  //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				  xv2.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[ii+16]),
						    alphav.v);
				  _mm512_store_pd(&x[ii+16],xv2.v);
				  xv3.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[ii+24]),
						    alphav.v);
				  _mm512_store_pd(&x[ii+24],xv3.v);
                                 

			     }

			     for(; (ii+15) < n; ii += 16) {

			          xv0.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[ii+0]),
						    alphav.v);
				  _mm512_store_pd(&x[ii+0],xv0.v);
				  xv1.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[ii+8]),
						    alphav.v);
				  _mm512_store_pd(&x[ii+8],xv1.v);
                                 

			     }

			     for(; (ii+7) < n; ii += 8) {

                                  xv0.v = _mm512_mul_pd(
				                    _mm512_load_pd(&x[ii+0]),
						    alphav.v);
				  _mm512_store_pd(&x[ii+0],xv0.v);

			     }

			     for(; (ii+0) < n; ii += 1) {
                                 x[ii] *= alpha;
			      }
			}
			 else {
                                   for(i = 0; i != n; ++i) {
                                        *x *= alpha;
				         x += incx;
			           }

			 }

		  }


		   
      }


}




#endif /*__GMS_SCALV_AVX512_UNROLL8X_HPP__*/
