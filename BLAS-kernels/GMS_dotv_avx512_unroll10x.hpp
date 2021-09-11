
#ifndef __GMS_DOTV_AVX512_UNROLL10X_HPP__
#define __GMS_DOTV_AVX512_UNROLL10X_HPP__


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

    const unsigned int gGMS_DOTV_AVX512_UNROLL10X_MAJOR = 1U;
    const unsigned int gGMS_DOTV_AVX512_UNROLL10X_MINOR = 0U;
    const unsigned int gGMS_DOTV_AVX512_UNROLL10X_MICRO = 0U;
    const unsigned int gGMS_DOTV_AVX512_UNROLL10X_FULLVER =
      1000U*gGMS_DOTV_AVX512_UNROLL10X_MAJOR+
      100U*gGMS_DOTV_AVX512_UNROLL10X_MINOR+
      10U*gGMS_DOTV_AVX512_UNROLL10X_MICRO;
    const char * const pgGMS_DOTV_AVX512_UNROLL10X_CREATION_DATE = "28-08-2021 11:00 AM +00200 (SAT 28 AUG 2021 GMT+2)";
    const char * const pgGMS_DOTV_AVX512_UNROLL10X_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_DOTV_AVX512_UNROLL10X_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_DOTV_AVX512_UNROLL10X_DESCRIPTION   = "AVX512 optimized DOTV kernels."

}

#include <cstdint>
#include <immintrin.h>
#include <omp.h>
#include "GMS_config.h"


namespace  gms {

          namespace math {

                   /*
                        Helper unions
                     */
	    namespace {

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
		   void sdotv_u_zmm16r4_unroll10x(const int32_t n,
		                                  float * __restrict x,
						  const int32_t incx,
						  float * __restrict y,
						  const int32_t incy,
						  float * __restrict rho) {
                              
                           if(__builtin_expect(0==n,0)) { return;}
			   
			   __ATTR_ALIGN__(64) __m512    xv[10];
			   __ATTR_ALIGN__(64) __m512    yv[10];
			   __ATTR_ALIGN__(64) zmm16r4_t rhov[10];
			   const __m512 vz = _mm512_setzero_ps();
			   float rho0;
			   int32_t i;
			   rho0 = 0.0f;
			   if(__builtin_expect(1==incx,1) &&
			      __builtin_expect(1==incy,1)) {

			        rhov[0].v = vz;
			        rhov[1].v = vz;
			        rhov[2].v = vz;
			        rhov[3].v = vz;
			        rhov[4].v = vz;
			        rhov[5].v = vz;
			        rhov[6].v = vz;
			        rhov[7].v = vz;
			        rhov[8].v = vz;
			        rhov[9].v = vz;

				for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                     _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				     xv[0]     = _mm512_loadu_ps(&x[i+0]);
				     yv[0]     = _mm512_loadu_ps(&y[i+0]);
				     rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				     xv[1]     = _mm512_loadu_ps(&x[i+16]);
				     yv[1]     = _mm512_loadu_ps(&y[i+16]);
				     rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v);
				     _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+64],_MM_HINT_T0);
				     xv[2]     = _mm512_loadu_ps(&x[i+32]);
				     yv[2]     = _mm512_loadu_ps(&y[i+32]);
				     rhov[2].v = _mm512_fmadd_ps(xv[2],yv[2],rhov[2].v);
				     xv[3]     = _mm512_loadu_ps(&x[i+48]);
				     yv[3]     = _mm512_loadu_ps(&y[i+48]);
				     rhov[3].v = _mm512_fmadd_ps(xv[3],yv[3],rhov[3].v);
				     _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+96],_MM_HINT_T0);
				     xv[4]     = _mm512_loadu_ps(&x[i+64]);
				     yv[4]     = _mm512_loadu_ps(&y[i+64]);
				     rhov[4].v = _mm512_fmadd_ps(xv[4],yv[4],rhov[4].v);
				     xv[5]     = _mm512_loadu_ps(&x[i+80]);
				     yv[5]     = _mm512_loadu_ps(&y[i+80]);
				     rhov[5].v = _mm512_fmadd_ps(xv[5],yv[5],rhov[5].v);
				     _mm_prefetch((const char*)&x[i+128],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+128],_MM_HINT_T0);
				     xv[6]     = _mm512_loadu_ps(&x[i+96]);
				     yv[6]     = _mm512_loadu_ps(&y[i+96]);
				     rhov[6].v = _mm512_fmadd_ps(xv[6],yv[6],rhov[6].v);
				     xv[7]     = _mm512_loadu_ps(&x[i+112]);
				     yv[7]     = _mm512_loadu_ps(&y[i+112]);
				     rhov[7].v = _mm512_fmadd_ps(xv[7],yv[7],rhov[7].v);
				     _mm_prefetch((const char*)&x[i+160],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+160],_MM_HINT_T0);
				     xv[8]     = _mm512_loadu_ps(&x[i+128]);
				     yv[8]     = _mm512_loadu_ps(&y[i+128]);
				     rhov[8].v = _mm512_fmadd_ps(xv[8],yv[8],rhov[8].v);
				     xv[9]     = _mm512_loadu_ps(&x[i+144]);
				     yv[9]     = _mm512_loadu_ps(&y[i+144]);
				     rhov[9].v = _mm512_fmadd_ps(xv[9],yv[9],rhov[9].v);
				     
#else
                                     _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+96],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+128],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+128],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+160],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+160],_MM_HINT_T0);
				     xv[0]     = _mm512_loadu_ps(&x[i+0]);
				     xv[1]     = _mm512_loadu_ps(&x[i+16]);
				     xv[2]     = _mm512_loadu_ps(&x[i+32]);
				     xv[3]     = _mm512_loadu_ps(&x[i+48]);
				     xv[4]     = _mm512_loadu_ps(&x[i+64]);
				     xv[5]     = _mm512_loadu_ps(&x[i+80]);
				     xv[6]     = _mm512_loadu_ps(&x[i+96]);
				     xv[7]     = _mm512_loadu_ps(&x[i+112]);
				     xv[8]     = _mm512_loadu_ps(&x[i+128]);
				     xv[9]     = _mm512_loadu_ps(&x[i+144]);
				     yv[0]     = _mm512_loadu_ps(&y[i+0]);
				     yv[1]     = _mm512_loadu_ps(&y[i+16]);
				     yv[2]     = _mm512_loadu_ps(&y[i+32]);
				     yv[3]     = _mm512_loadu_ps(&y[i+48]);
				     yv[4]     = _mm512_loadu_ps(&y[i+64]);
				     yv[5]     = _mm512_loadu_ps(&y[i+80]);
				     yv[6]     = _mm512_loadu_ps(&y[i+96]);
				     yv[7]     = _mm512_loadu_ps(&y[i+112]);
				     yv[8]     = _mm512_loadu_ps(&y[i+128]);
				     yv[9]     = _mm512_loadu_ps(&y[i+144]);
				     rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				     rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v);
				     rhov[2].v = _mm512_fmadd_ps(xv[2],yv[2],rhov[2].v);
				     rhov[3].v = _mm512_fmadd_ps(xv[3],yv[3],rhov[3].v);
				     rhov[4].v = _mm512_fmadd_ps(xv[4],yv[4],rhov[4].v);
				     rhov[5].v = _mm512_fmadd_ps(xv[5],yv[5],rhov[5].v);
				     rhov[6].v = _mm512_fmadd_ps(xv[6],yv[6],rhov[6].v);
				     rhov[7].v = _mm512_fmadd_ps(xv[7],yv[7],rhov[7].v);
				     rhov[8].v = _mm512_fmadd_ps(xv[8],yv[8],rhov[8].v);
				     rhov[9].v = _mm512_fmadd_ps(xv[9],yv[9],rhov[9].v);
				     
#endif

				}

				      rhov[0].v += rhov[5].v;
		                      rhov[1].v += rhov[6].v;
		                      rhov[2].v += rhov[7].v;
		                      rhov[3].v += rhov[8].v;
		                      rhov[4].v += rhov[9].v;

				for(; (i+79) < n; i += 80) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                       xv[0]     = _mm512_loadu_ps(&x[i+0]);
				       yv[0]     = _mm512_loadu_ps(&y[i+0]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				       xv[1]     = _mm512_loadu_ps(&x[i+16]);
				       yv[1]     = _mm512_loadu_ps(&y[i+16]);
				       rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v);
				       //_mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				       //_mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				       xv[2]     = _mm512_loadu_ps(&x[i+32]);
				       yv[2]     = _mm512_loadu_ps(&y[i+32]);
				       rhov[2].v = _mm512_fmadd_ps(xv[2],yv[2],rhov[2].v);
				       xv[3]     = _mm512_loadu_ps(&x[i+48]);
				       yv[3]     = _mm512_loadu_ps(&y[i+48]);
				       rhov[3].v = _mm512_fmadd_ps(xv[3],yv[3],rhov[3].v);
				       //_mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				       //_mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				       xv[4]     = _mm512_loadu_ps(&x[i+64]);
				       yv[4]     = _mm512_loadu_ps(&y[i+64]);
				       rhov[4].v = _mm512_fmadd_ps(xv[4],yv[4],rhov[4].v);
#else 
                                       xv[0]     = _mm512_loadu_ps(&x[i+0]);
				       xv[1]     = _mm512_loadu_ps(&x[i+16]);
				       xv[2]     = _mm512_loadu_ps(&x[i+32]);
				       xv[3]     = _mm512_loadu_ps(&x[i+48]);
				       xv[4]     = _mm512_loadu_ps(&x[i+64]);
				       yv[0]     = _mm512_loadu_ps(&y[i+0]);
				       yv[1]     = _mm512_loadu_ps(&y[i+16]);
				       yv[2]     = _mm512_loadu_ps(&y[i+32]);
				       yv[3]     = _mm512_loadu_ps(&y[i+48]);
				       yv[4]     = _mm512_loadu_ps(&y[i+64]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				       rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v);
				       rhov[2].v = _mm512_fmadd_ps(xv[2],yv[2],rhov[2].v);
				       rhov[3].v = _mm512_fmadd_ps(xv[3],yv[3],rhov[3].v);
				       rhov[4].v = _mm512_fmadd_ps(xv[4],yv[4],rhov[4].v);
#endif

				}

				     rhov[0].v += rhov[2].v;
		                     rhov[1].v += rhov[3].v;
		                     rhov[0].v += rhov[4].v;

				for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                       xv[0]     = _mm512_loadu_ps(&x[i+0]);
				       yv[0]     = _mm512_loadu_ps(&y[i+0]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				       xv[1]     = _mm512_loadu_ps(&x[i+16]);
				       yv[1]     = _mm512_loadu_ps(&y[i+16]);
				       rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v); 
#else
                                       xv[0]     = _mm512_loadu_ps(&x[i+0]);
				       xv[1]     = _mm512_loadu_ps(&x[i+16]);
				       yv[0]     = _mm512_loadu_ps(&y[i+0]);
				       yv[1]     = _mm512_loadu_ps(&y[i+16]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				       rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v);
#endif

				}

				    	rhov[0].v += rhov[1].v;

			        for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                       xv[0]     = _mm512_loadu_ps(&x[i+0]);
				       yv[0]     = _mm512_loadu_ps(&y[i+0]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
#else
                                       xv[0]     = _mm512_loadu_ps(&x[i+0]);
				       yv[0]     = _mm512_loadu_ps(&y[i+0]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
#endif
				}

				for(; (i+0) < n; i += 1) {
                                      rho0 += x[i] * y[i];
				}

				rho0 += _mm512_reduce_add_ps(rhov[0]);

			 }
			 else {

			         for (i = 0; i < n; ++i ) {
					const float x0c = *x0;
			                const float y0c = *y0;
			                rho0 += x0c * y0c;
			                x0 += incx;
			                y0 += incy;
		                    }
			 }
			 *rho = rho0;

		   }


		   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void sdotv_a_zmm16r4_unroll10x(const int32_t n,
		                                  float * __restrict __ATTR_ALIGN__(64) x,
						  const int32_t incx,
						  float * __restrict __ATTR_ALIGN__(64) y,
						  const int32_t incy,
						  float * __restrict rho) {
                              
                           if(__builtin_expect(0==n,0)) { return;}
			   
			   __ATTR_ALIGN__(64) __m512    xv[10];
			   __ATTR_ALIGN__(64) __m512    yv[10];
			   __ATTR_ALIGN__(64) zmm16r4_t rhov[10];
			   const __m512 vz = _mm512_setzero_ps();
			   float rho0;
			   int32_t i;
			   rho0 = 0.0f;
			   if(__builtin_expect(1==incx,1) &&
			      __builtin_expect(1==incy,1)) {

			        rhov[0].v = vz;
			        rhov[1].v = vz;
			        rhov[2].v = vz;
			        rhov[3].v = vz;
			        rhov[4].v = vz;
			        rhov[5].v = vz;
			        rhov[6].v = vz;
			        rhov[7].v = vz;
			        rhov[8].v = vz;
			        rhov[9].v = vz;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,32);
			     __assume_aligned(y,32);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (float*)__builtin_assume_aligned(x,32);
			     y = (float*)__builtin_assume_aligned(y,32);
#endif
				for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                     _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				     xv[0]     = _mm512_load_ps(&x[i+0]);
				     yv[0]     = _mm512_load_ps(&y[i+0]);
				     rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				     xv[1]     = _mm512_load_ps(&x[i+16]);
				     yv[1]     = _mm512_load_ps(&y[i+16]);
				     rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v);
				     _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+64],_MM_HINT_T0);
				     xv[2]     = _mm512_load_ps(&x[i+32]);
				     yv[2]     = _mm512_load_ps(&y[i+32]);
				     rhov[2].v = _mm512_fmadd_ps(xv[2],yv[2],rhov[2].v);
				     xv[3]     = _mm512_load_ps(&x[i+48]);
				     yv[3]     = _mm512_load_ps(&y[i+48]);
				     rhov[3].v = _mm512_fmadd_ps(xv[3],yv[3],rhov[3].v);
				     _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+96],_MM_HINT_T0);
				     xv[4]     = _mm512_load_ps(&x[i+64]);
				     yv[4]     = _mm512_load_ps(&y[i+64]);
				     rhov[4].v = _mm512_fmadd_ps(xv[4],yv[4],rhov[4].v);
				     xv[5]     = _mm512_load_ps(&x[i+80]);
				     yv[5]     = _mm512_load_ps(&y[i+80]);
				     rhov[5].v = _mm512_fmadd_ps(xv[5],yv[5],rhov[5].v);
				     _mm_prefetch((const char*)&x[i+128],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+128],_MM_HINT_T0);
				     xv[6]     = _mm512_load_ps(&x[i+96]);
				     yv[6]     = _mm512_load_ps(&y[i+96]);
				     rhov[6].v = _mm512_fmadd_ps(xv[6],yv[6],rhov[6].v);
				     xv[7]     = _mm512_load_ps(&x[i+112]);
				     yv[7]     = _mm512_load_ps(&y[i+112]);
				     rhov[7].v = _mm512_fmadd_ps(xv[7],yv[7],rhov[7].v);
				     _mm_prefetch((const char*)&x[i+160],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+160],_MM_HINT_T0);
				     xv[8]     = _mm512_load_ps(&x[i+128]);
				     yv[8]     = _mm512_load_ps(&y[i+128]);
				     rhov[8].v = _mm512_fmadd_ps(xv[8],yv[8],rhov[8].v);
				     xv[9]     = _mm512_load_ps(&x[i+144]);
				     yv[9]     = _mm512_load_ps(&y[i+144]);
				     rhov[9].v = _mm512_fmadd_ps(xv[9],yv[9],rhov[9].v);
				     
#else
                                     _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+96],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+128],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+128],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+160],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+160],_MM_HINT_T0);
				     xv[0]     = _mm512_load_ps(&x[i+0]);
				     xv[1]     = _mm512_load_ps(&x[i+16]);
				     xv[2]     = _mm512_load_ps(&x[i+32]);
				     xv[3]     = _mm512_load_ps(&x[i+48]);
				     xv[4]     = _mm512_load_ps(&x[i+64]);
				     xv[5]     = _mm512_load_ps(&x[i+80]);
				     xv[6]     = _mm512_load_ps(&x[i+96]);
				     xv[7]     = _mm512_load_ps(&x[i+112]);
				     xv[8]     = _mm512_load_ps(&x[i+128]);
				     xv[9]     = _mm512_load_ps(&x[i+144]);
				     yv[0]     = _mm512_load_ps(&y[i+0]);
				     yv[1]     = _mm512_load_ps(&y[i+16]);
				     yv[2]     = _mm512_load_ps(&y[i+32]);
				     yv[3]     = _mm512_load_ps(&y[i+48]);
				     yv[4]     = _mm512_load_ps(&y[i+64]);
				     yv[5]     = _mm512_load_ps(&y[i+80]);
				     yv[6]     = _mm512_load_ps(&y[i+96]);
				     yv[7]     = _mm512_load_ps(&y[i+112]);
				     yv[8]     = _mm512_load_ps(&y[i+128]);
				     yv[9]     = _mm512_load_ps(&y[i+144]);
				     rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				     rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v);
				     rhov[2].v = _mm512_fmadd_ps(xv[2],yv[2],rhov[2].v);
				     rhov[3].v = _mm512_fmadd_ps(xv[3],yv[3],rhov[3].v);
				     rhov[4].v = _mm512_fmadd_ps(xv[4],yv[4],rhov[4].v);
				     rhov[5].v = _mm512_fmadd_ps(xv[5],yv[5],rhov[5].v);
				     rhov[6].v = _mm512_fmadd_ps(xv[6],yv[6],rhov[6].v);
				     rhov[7].v = _mm512_fmadd_ps(xv[7],yv[7],rhov[7].v);
				     rhov[8].v = _mm512_fmadd_ps(xv[8],yv[8],rhov[8].v);
				     rhov[9].v = _mm512_fmadd_ps(xv[9],yv[9],rhov[9].v);
				     
#endif

				}

				      rhov[0].v += rhov[5].v;
		                      rhov[1].v += rhov[6].v;
		                      rhov[2].v += rhov[7].v;
		                      rhov[3].v += rhov[8].v;
		                      rhov[4].v += rhov[9].v;

				for(; (i+79) < n; i += 80) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                       xv[0]     = _mm512_load_ps(&x[i+0]);
				       yv[0]     = _mm512_load_ps(&y[i+0]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				       xv[1]     = _mm512_load_ps(&x[i+16]);
				       yv[1]     = _mm512_load_ps(&y[i+16]);
				       rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v);
				       //_mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				       //_mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				       xv[2]     = _mm512_load_ps(&x[i+32]);
				       yv[2]     = _mm512_load_ps(&y[i+32]);
				       rhov[2].v = _mm512_fmadd_ps(xv[2],yv[2],rhov[2].v);
				       xv[3]     = _mm512_load_ps(&x[i+48]);
				       yv[3]     = _mm512_load_ps(&y[i+48]);
				       rhov[3].v = _mm512_fmadd_ps(xv[3],yv[3],rhov[3].v);
				       //_mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				       //_mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				       xv[4]     = _mm512_load_ps(&x[i+64]);
				       yv[4]     = _mm512_load_ps(&y[i+64]);
				       rhov[4].v = _mm512_fmadd_ps(xv[4],yv[4],rhov[4].v);
#else 
                                       xv[0]     = _mm512_load_ps(&x[i+0]);
				       xv[1]     = _mm512_load_ps(&x[i+16]);
				       xv[2]     = _mm512_load_ps(&x[i+32]);
				       xv[3]     = _mm512_load_ps(&x[i+48]);
				       xv[4]     = _mm512_load_ps(&x[i+64]);
				       yv[0]     = _mm512_load_ps(&y[i+0]);
				       yv[1]     = _mm512_load_ps(&y[i+16]);
				       yv[2]     = _mm512_load_ps(&y[i+32]);
				       yv[3]     = _mm512_load_ps(&y[i+48]);
				       yv[4]     = _mm512_load_ps(&y[i+64]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				       rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v);
				       rhov[2].v = _mm512_fmadd_ps(xv[2],yv[2],rhov[2].v);
				       rhov[3].v = _mm512_fmadd_ps(xv[3],yv[3],rhov[3].v);
				       rhov[4].v = _mm512_fmadd_ps(xv[4],yv[4],rhov[4].v);
#endif

				}

				     rhov[0].v += rhov[2].v;
		                     rhov[1].v += rhov[3].v;
		                     rhov[0].v += rhov[4].v;

				for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                       xv[0]     = _mm512_load_ps(&x[i+0]);
				       yv[0]     = _mm512_load_ps(&y[i+0]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				       xv[1]     = _mm512_load_ps(&x[i+16]);
				       yv[1]     = _mm512_load_ps(&y[i+16]);
				       rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v); 
#else
                                       xv[0]     = _mm512_load_ps(&x[i+0]);
				       xv[1]     = _mm512_load_ps(&x[i+16]);
				       yv[0]     = _mm512_load_ps(&y[i+0]);
				       yv[1]     = _mm512_load_ps(&y[i+16]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
				       rhov[1].v = _mm512_fmadd_ps(xv[1],yv[1],rhov[1].v);
#endif

				}

				    	rhov[0].v += rhov[1].v;

			        for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                       xv[0]     = _mm512_load_ps(&x[i+0]);
				       yv[0]     = _mm512_load_ps(&y[i+0]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
#else
                                       xv[0]     = _mm512_load_ps(&x[i+0]);
				       yv[0]     = _mm512_load_ps(&y[i+0]);
				       rhov[0].v = _mm512_fmadd_ps(xv[0],yv[0],rhov[0].v);
#endif
				}

				for(; (i+0) < n; i += 1) {
                                      rho0 += x[i] * y[i];
				}

				rho0 += _mm512_reduce_add_ps(rhov[0]);

			 }
			 else {

			         for (i = 0; i < n; ++i ) {
					const float x0c = *x0;
			                const float y0c = *y0;
			                rho0 += x0c * y0c;
			                x0 += incx;
			                y0 += incy;
		                    }
			 }
			 *rho = rho0;

		   }



		   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void sdotv_a_zmm16r4_unroll10x_omp(const int32_t n,
		                                      float * __restrict __ATTR_ALIGN__(64) x,
						      const int32_t incx,
						      float * __restrict __ATTR_ALIGN__(64) y,
						      const int32_t incy,
						      float * __restrict rho) {
                              
                           if(__builtin_expect(0==n,0)) { return;}
			   
			   __ATTR_ALIGN__(64) __m512    xv[10];
			   __ATTR_ALIGN__(64) __m512    yv[10];
			   //__ATTR_ALIGN__(64) zmm16r4_t rhov[10];
			   zmm16r4_t rhov0;
			   zmm16r4_t rhov1;
			   zmm16r4_t rhov2;
			   zmm16r4_t rhov3;
			   zmm16r4_t rhov4;
			   zmm16r4_t rhov5;
			   zmm16r4_t rhov6;
			   zmm16r4_t rhov7;
			   zmm16r4_t rhov8;
			   zmm16r4_t rhov9;
			   const __m512 vz = _mm512_setzero_ps();
			   float rho0;
			   int32_t i;
			   int32_t ii;
			   ii = 0;
			   rho0 = 0.0f;
			   if(__builtin_expect(1==incx,1) &&
			      __builtin_expect(1==incy,1)) {

			       rhov0 = vz;
			       rhov1 = vz;
			       rhov2 = vz;
			       rhov3 = vz;
			       rhov4 = vz;
			       rhov5 = vz;
			       rhov6 = vz;
			       rhov7 = vz;
			       rhov8 = vz;
			       rhov9 = vz;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,32);
			     __assume_aligned(y,32);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (float*)__builtin_assume_aligned(x,32);
			     y = (float*)__builtin_assume_aligned(y,32);
#endif
#pragma omp parallel for schedule(static,160) default(none) private(i) shared(x,y,n) \
            lastprivate(ii,rhov0,rhov1,rhov2,rhov3,rhov4,rhov5,rhov6,rhov7,rhov8,rhov9)
				for(i = 0; (i+159) < n; i += 160) {
                                      ii = i;
                                     _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				     rhov0.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+0]),
				                              _mm512_load_ps(&y[i+0]),
							      rhov0.v);               
				     rhov1.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+16]),
				                              _mm512_load_ps(&y[i+16]),
							      rhov1.v);   
				     _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+64],_MM_HINT_T0);
				     rhov2.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+32]),
				                              _mm512_load_ps(&y[i+32]),
							      rhov2.v);
				     rhov3.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+48]),
				                              _mm512_load_ps(&y[i+48]),
							      rhov3.v);   
				     _mm_prefetch((const char*)&x[i+96],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+96],_MM_HINT_T0);
				     rhov4.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+64]),
				                              _mm512_load_ps(&y[i+64]),
							      rhov4.v);
				     rhov5.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+80]),
				                              _mm512_load_ps(&y[i+80]),
							      rhov5.v);   
				     _mm_prefetch((const char*)&x[i+128],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+128],_MM_HINT_T0);
				     rhov6.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+96]),
				                              _mm512_load_ps(&y[i+96]),
							      rhov6.v);   
				     rhov7.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+112]),
				                              _mm512_load_ps(&y[i+112]),
							      rhov7.v);   
				     _mm_prefetch((const char*)&x[i+160],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+160],_MM_HINT_T0);
				     rhov8.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+128]),
				                              _mm512_load_ps(&y[i+128]),
							      rhov8.v);   
				     rhov9.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+144]),
				                              _mm512_load_ps(&y[i+144]),
							      rhov9.v);   
				    				                                     
				}

				      rhov0.v += rhov.v;
		                      rhov1.v += rhov6.v;
		                      rhov2.v += rhov7.v;
		                      rhov3.v += rhov8.v;
		                      rhov4.v += rhov9.v;

				for(; (i+79) < n; i += 80) {

                                     rhov0.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+0]),
				                              _mm512_load_ps(&y[i+0]),
							      rhov0.v);               
				     rhov1.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+16]),
				                              _mm512_load_ps(&y[i+16]),
							      rhov1.v);   
				    
				     rhov2.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+32]),
				                              _mm512_load_ps(&y[i+32]),
							      rhov2.v);
				     rhov3.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+48]),
				                              _mm512_load_ps(&y[i+48]),
							      rhov3.v);   
				      rhov4.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+64]),
				                              _mm512_load_ps(&y[i+64]),
							      rhov4.v);
				
                                     
				}

				     rhov0.v += rhov2.v;
		                     rhov1.v += rhov3.v;
		                     rhov0.v += rhov4.v;

				for(; (i+31) < n; i += 32) {

                                     rhov0.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+0]),
				                              _mm512_load_ps(&y[i+0]),
							      rhov0.v);               
				     rhov1.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+16]),
				                              _mm512_load_ps(&y[i+16]),
							      rhov1.v);   

				}

				    	rhov0.v += rhov1.v;

			        for(; (i+15) < n; i += 16) {

                                      rhov0.v =
				              _mm512_fmadd_ps(_mm512_load_ps(&x[i+0]),
				                              _mm512_load_ps(&y[i+0]),
							      rhov0.v);   

				}

				for(; (i+0) < n; i += 1) {
                                      rho0 += x[i] * y[i];
				}

				rho0 += _mm512_reduce_add_ps(rhov0.v);

			 }
			 else {

			         for (i = 0; i < n; ++i ) {
					const float x0c = *x0;
			                const float y0c = *y0;
			                rho0 += x0c * y0c;
			                x0 += incx;
			                y0 += incy;
		                    }
			 }
			 *rho = rho0;

		   }



		   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void ddotv_u_zmm8r8_unroll10x(const int32_t,
		                                 double * __restrict x,
						 const int32_t incx,
						 double * __restrict y,
						 const int32_t incy,
						 double * __restrict rho) {

                             if(__builtin_expect(0==n,0)) { return;}
			   
			   __ATTR_ALIGN__(64) __m512d    xv[10];
			   __ATTR_ALIGN__(64) __m512d    yv[10];
			   __ATTR_ALIGN__(64) zmm8r8_t rhov[10];
			   const __m512d vz = _mm512_setzero_pd();
			   double rho0;
			   int32_t i;
			   rho0 = 0.0;
			   if(__builtin_expect(1==incx,1) &&
			      __builtin_expect(1==incy,1)) {

			        rhov[0].v = vz;
			        rhov[1].v = vz;
			        rhov[2].v = vz;
			        rhov[3].v = vz;
			        rhov[4].v = vz;
			        rhov[5].v = vz;
			        rhov[6].v = vz;
			        rhov[7].v = vz;
			        rhov[8].v = vz;
			        rhov[9].v = vz;

				for(i = 0; (i+79) < n; i += 80) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                     _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+16],_MM_HINT_T0);
				     xv[0] = _mm512_loadu_pd(&x[i+0]);
				     yv[0] = _mm512_loadu_pd(&y[i+0]);
				     rhov[0].v = _mm512_fmadd_pd(xv[0],yv[0],rhov[0].v);
				     xv[1] = _mm512_loadu_pd(&x[i+8]);
				     yv[1] = _mm512_loadu_pd(&y[i+8]);
				     rhov[1].v = _mm512_fmadd_pd(xv[1],yv[1],rhov[1].v);
				     _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				     xv[2] = _mm512_loadu_pd(&x[i+16]);
				     yv[2] = _mm512_loadu_pd(&y[i+16]);
				     rhov[2].v = _mm512_fmadd_pd(xv[2],yv[2],rhov[2].v);
				     xv[3] = _mm512_loadu_pd(&x[i+24]);
				     yv[3] = _mm512_loadu_pd(&y[i+24]);
				     rhov[3].v = _mm512_fmadd_pd(xv[3],yv[3],rhov[3].v);
				     _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+48],_MM_HINT_T0);
				     xv[4] = _mm512_loadu_pd(&x[i+32]);
				     yv[4] = _mm512_loadu_pd(&y[i+32]);
				     rhov[4].v = _mm512_fmadd_pd(xv[4],yv[4],rhov[4].v);
				     xv[5] = _mm512_loadu_pd(&x[i+40]);
				     yv[5] = _mm512_loadu_pd(&y[i+40]);
				     rhov[5].v = _mm512_fmadd_pd(xv[5],yv[5],rhov[5].v);
				     _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+64],_MM_HINT_T0);
				     xv[6] = _mm512_loadu_pd(&x[i+48]);
				     yv[6] = _mm512_loadu_pd(&y[i+48]);
				     rhov[6].v = _mm512_fmadd_pd(xv[6],yv[6],rhov[6].v);
				     xv[7] = _mm512_loadu_pd(&x[i+56]);
				     yv[7] = _mm512_loadu_pd(&y[i+56]);
				     rhov[7].v = _mm512_fmadd_pd(xv[7],yv[7],rhov[7].v);
				     _mm_prefetch((const char*)&x[i+80],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+80],_MM_HINT_T0);
				     xv[8] = _mm512_loadu_pd(&x[i+64]);
				     yv[8] = _mm512_loadu_pd(&y[i+64]);
				     rhov[8].v = _mm512_fmadd_pd(xv[8],yv[8],rhov[8].v);
				     xv[9] = _mm512_loadu_pd(&x[i+72]);
				     yv[9] = _mm512_loadu_pd(&y[i+72]);
				     rhov[9].v = _mm512_fmadd_pd(xv[9],yv[9],rhov[9].v);
#else
                                     _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+16],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+48],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+80],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+80],_MM_HINT_T0);
				     xv[0] = _mm512_loadu_pd(&x[i+0]);
				     xv[1] = _mm512_loadu_pd(&x[i+8]);
				     xv[2] = _mm512_loadu_pd(&x[i+16]);
				     xv[3] = _mm512_loadu_pd(&x[i+24]);
				     xv[4] = _mm512_loadu_pd(&x[i+32]);
				     xv[5] = _mm512_loadu_pd(&x[i+40]);
				     xv[6] = _mm512_loadu_pd(&x[i+48]);
				     xv[7] = _mm512_loadu_pd(&x[i+56]);
				     xv[8] = _mm512_loadu_pd(&x[i+64]);
				     xv[9] = _mm512_loadu_pd(&x[i+72]);
				     yv[0] = _mm512_loadu_pd(&y[i+0]);
				     yv[1] = _mm512_loadu_pd(&y[i+8]);
				     yv[2] = _mm512_loadu_pd(&y[i+16]);
				     yv[3] = _mm512_loadu_pd(&y[i+24]);
				     yv[4] = _mm512_loadu_pd(&y[i+32]);
				     yv[5] = _mm512_loadu_pd(&y[i+40]);
				     yv[6] = _mm512_loadu_pd(&y[i+48]);
				     yv[7] = _mm512_loadu_pd(&y[i+56]);
				     yv[8] = _mm512_loadu_pd(&y[i+64]);
				     yv[9] = _mm512_loadu_pd(&y[i+72]);
				     rhov[0].v = _mm256_fmadd_pd(xv[0],yv[0],rhov[0].v);
			             rhov[1].v = _mm256_fmadd_pd(xv[1],yv[1],rhov[1].v);
			             rhov[2].v = _mm256_fmadd_pd(xv[2],yv[2],rhov[2].v);
			             rhov[3].v = _mm256_fmadd_pd(xv[3],yv[3],rhov[3].v);
			             rhov[4].v = _mm256_fmadd_pd(xv[4],yv[4],rhov[4].v);
			             rhov[5].v = _mm256_fmadd_pd(xv[5],yv[5],rhov[5].v);
			             rhov[6].v = _mm256_fmadd_pd(xv[6],yv[6],rhov[6].v);
			             rhov[7].v = _mm256_fmadd_pd(xv[7],yv[7],rhov[7].v);
			             rhov[8].v = _mm256_fmadd_pd(xv[8],yv[8] rhov[8].v);
			             rhov[9].v = _mm256_fmadd_pd(xv[9],yv[9],rhov[9].v);
#endif

				}

				   	rhov[0].v += rhov[5].v;
		                        rhov[1].v += rhov[6].v;
		                        rhov[2].v += rhov[7].v;
		                        rhov[3].v += rhov[8].v;
		                        rhov[4].v += rhov[9].v;

			       for(; (i+39) < n; i += 40) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                     xv[0] = _mm512_loadu_pd(&x[i+0]);
				     yv[0] = _mm512_loadu_pd(&y[i+0]);
				     rhov[0].v = _mm512_fmadd_pd(xv[0],yv[0],rhov[0].v);
				     xv[1] = _mm512_loadu_pd(&x[i+8]);
				     yv[1] = _mm512_loadu_pd(&y[i+8]);
				     rhov[1].v = _mm512_fmadd_pd(xv[1],yv[1],rhov[1].v);
				     //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     //_mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				     xv[2] = _mm512_loadu_pd(&x[i+16]);
				     yv[2] = _mm512_loadu_pd(&y[i+16]);
				     rhov[2].v = _mm512_fmadd_pd(xv[2],yv[2],rhov[2].v);
				     xv[3] = _mm512_loadu_pd(&x[i+24]);
				     yv[3] = _mm512_loadu_pd(&y[i+24]);
				     rhov[3].v = _mm512_fmadd_pd(xv[3],yv[3],rhov[3].v);
				     //_mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				     // _mm_prefetch((const char*)&y[i+48],_MM_HINT_T0);
				     xv[4] = _mm512_loadu_pd(&x[i+32]);
				     yv[4] = _mm512_loadu_pd(&y[i+32]);
				     rhov[4].v = _mm512_fmadd_pd(xv[4],yv[4],rhov[4].v);
#else
                                     xv[0] = _mm512_loadu_pd(&x[i+0]);
				     xv[1] = _mm512_loadu_pd(&x[i+8]);
				     xv[2] = _mm512_loadu_pd(&x[i+16]);
				     xv[3] = _mm512_loadu_pd(&x[i+24]);
				     xv[4] = _mm512_loadu_pd(&x[i+32]);
				     yv[0] = _mm512_loadu_pd(&y[i+0]);
				     yv[1] = _mm512_loadu_pd(&y[i+8]);
				     yv[2] = _mm512_loadu_pd(&y[i+16]);
				     yv[3] = _mm512_loadu_pd(&y[i+24]);
				     yv[4] = _mm512_loadu_pd(&y[i+32]);
				     rhov[0].v = _mm256_fmadd_pd(xv[0],yv[0],rhov[0].v);
			             rhov[1].v = _mm256_fmadd_pd(xv[1],yv[1],rhov[1].v);
			             rhov[2].v = _mm256_fmadd_pd(xv[2],yv[2],rhov[2].v);
			             rhov[3].v = _mm256_fmadd_pd(xv[3],yv[3],rhov[3].v);
			             rhov[4].v = _mm256_fmadd_pd(xv[4],yv[4],rhov[4].v);
#endif

				 }

				     rhov[0].v += rhov[2].v;
				     rhov[1].v += rhov[3].v;
				     rhov[0].v += rhov[4].v;

				 for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                        xv[0] = _mm512_loadu_pd(&x[i+0]);
				        yv[0] = _mm512_loadu_pd(&y[i+0]);
				        rhov[0].v = _mm512_fmadd_pd(xv[0],yv[0],rhov[0].v);
				        xv[1] = _mm512_loadu_pd(&x[i+8]);
				        yv[1] = _mm512_loadu_pd(&y[i+8]);
				        rhov[1].v = _mm512_fmadd_pd(xv[1],yv[1],rhov[1].v);
				     //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     //_mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				        xv[2] = _mm512_loadu_pd(&x[i+16]);
				        yv[2] = _mm512_loadu_pd(&y[i+16]);
				        rhov[2].v = _mm512_fmadd_pd(xv[2],yv[2],rhov[2].v);
				        xv[3] = _mm512_loadu_pd(&x[i+24]);
				        yv[3] = _mm512_loadu_pd(&y[i+24]);
				        rhov[3].v = _mm512_fmadd_pd(xv[3],yv[3],rhov[3].v);
#else
                                        xv[0] = _mm512_loadu_pd(&x[i+0]);
				        xv[1] = _mm512_loadu_pd(&x[i+8]);
				        xv[2] = _mm512_loadu_pd(&x[i+16]);
				        xv[3] = _mm512_loadu_pd(&x[i+24]);
					yv[0] = _mm512_loadu_pd(&y[i+0]);
				        yv[1] = _mm512_loadu_pd(&y[i+8]);
				        yv[2] = _mm512_loadu_pd(&y[i+16]);
				        yv[3] = _mm512_loadu_pd(&y[i+24]);
					rhov[0].v = _mm256_fmadd_pd(xv[0],yv[0],rhov[0].v);
			                rhov[1].v = _mm256_fmadd_pd(xv[1],yv[1],rhov[1].v);
			                rhov[2].v = _mm256_fmadd_pd(xv[2],yv[2],rhov[2].v);
			                rhov[3].v = _mm256_fmadd_pd(xv[3],yv[3],rhov[3].v);
#endif

				 }

				     rhov[0].v += rhov[2].v;
				     rhov[1].v += rhov[3].v;

				 for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                        xv[0] = _mm512_loadu_pd(&x[i+0]);
				        yv[0] = _mm512_loadu_pd(&y[i+0]);
				        rhov[0].v = _mm512_fmadd_pd(xv[0],yv[0],rhov[0].v);
				        xv[1] = _mm512_loadu_pd(&x[i+8]);
				        yv[1] = _mm512_loadu_pd(&y[i+8]);
				        rhov[1].v = _mm512_fmadd_pd(xv[1],yv[1],rhov[1].v);
#else
                                        xv[0] = _mm512_loadu_pd(&x[i+0]);
				        xv[1] = _mm512_loadu_pd(&x[i+8]);
					yv[0] = _mm512_loadu_pd(&y[i+0]);
				        yv[1] = _mm512_loadu_pd(&y[i+8]);
					rhov[0].v = _mm256_fmadd_pd(xv[0],yv[0],rhov[0].v);
			                rhov[1].v = _mm256_fmadd_pd(xv[1],yv[1],rhov[1].v);
#endif
				 }
				     
                                       rhov[0].v += rhov[1].v;
			  

			        for(; (i+7) < n; i += 8) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                        xv[0] = _mm512_loadu_pd(&x[i+0]);
				        yv[0] = _mm512_loadu_pd(&y[i+0]);
				        rhov[0].v = _mm512_fmadd_pd(xv[0],yv[0],rhov[0].v);
#else
                                         xv[0] = _mm512_loadu_pd(&x[i+0]);
					 yv[0] = _mm512_loadu_pd(&y[i+0]);
					 rhov[0].v = _mm256_fmadd_pd(xv[0],yv[0],rhov[0].v);	
#endif
				}

				for(; (i+0) < n; i += 1) {
                                      rho0 += x[i] * y[i];
				}

				  rho0 += _mm512_reduce_add_pd(rhov[0].v);

			  }
			  else {

			           	for (i = 0; i < n; ++i) {
					     const double x0c = *x0;
			                     const double y0c = *y0;
			                     rho0 += x0c * y0c;
			                     x0 += incx;
			                     y0 += incy;
				      }
			  }

			  *rho = rho0;

		   }




		   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void ddotv_a_zmm8r8_unroll10x(const int32_t,
		                                 double * __restrict __ATTR_ALIGN__(64) x,
						 const int32_t incx,
						 double * __restrict __ATTR_ALIGN__(64)y,
						 const int32_t incy,
						 double * __restrict rho) {

                             if(__builtin_expect(0==n,0)) { return;}
			   
			   __ATTR_ALIGN__(64) __m512d    xv[10];
			   __ATTR_ALIGN__(64) __m512d    yv[10];
			   __ATTR_ALIGN__(64) zmm8r8_t rhov[10];
			   const __m512d vz = _mm512_setzero_pd();
			   double rho0;
			   int32_t i;
			   rho0 = 0.0;
			   if(__builtin_expect(1==incx,1) &&
			      __builtin_expect(1==incy,1)) {

			        rhov[0].v = vz;
			        rhov[1].v = vz;
			        rhov[2].v = vz;
			        rhov[3].v = vz;
			        rhov[4].v = vz;
			        rhov[5].v = vz;
			        rhov[6].v = vz;
			        rhov[7].v = vz;
			        rhov[8].v = vz;
			        rhov[9].v = vz;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,64);
			     __assume_aligned(y,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (double*)__builtin_assume_aligned(x,64);
			     y = (double*)__builtin_assume_aligned(y,64);
#endif
				for(i = 0; (i+79) < n; i += 80) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                     _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+16],_MM_HINT_T0);
				     xv[0] = _mm512_load_pd(&x[i+0]);
				     yv[0] = _mm512_load_pd(&y[i+0]);
				     rhov[0].v = _mm512_fmadd_pd(xv[0],yv[0],rhov[0].v);
				     xv[1] = _mm512_load_pd(&x[i+8]);
				     yv[1] = _mm512_load_pd(&y[i+8]);
				     rhov[1].v = _mm512_fmadd_pd(xv[1],yv[1],rhov[1].v);
				     _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				     xv[2] = _mm512_load_pd(&x[i+16]);
				     yv[2] = _mm512_load_pd(&y[i+16]);
				     rhov[2].v = _mm512_fmadd_pd(xv[2],yv[2],rhov[2].v);
				     xv[3] = _mm512_load_pd(&x[i+24]);
				     yv[3] = _mm512_load_pd(&y[i+24]);
				     rhov[3].v = _mm512_fmadd_pd(xv[3],yv[3],rhov[3].v);
				     _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+48],_MM_HINT_T0);
				     xv[4] = _mm512_load_pd(&x[i+32]);
				     yv[4] = _mm512_load_pd(&y[i+32]);
				     rhov[4].v = _mm512_fmadd_pd(xv[4],yv[4],rhov[4].v);
				     xv[5] = _mm512_load_pd(&x[i+40]);
				     yv[5] = _mm512_load_pd(&y[i+40]);
				     rhov[5].v = _mm512_fmadd_pd(xv[5],yv[5],rhov[5].v);
				     _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+64],_MM_HINT_T0);
				     xv[6] = _mm512_load_pd(&x[i+48]);
				     yv[6] = _mm512_load_pd(&y[i+48]);
				     rhov[6].v = _mm512_fmadd_pd(xv[6],yv[6],rhov[6].v);
				     xv[7] = _mm512_load_pd(&x[i+56]);
				     yv[7] = _mm512_load_pd(&y[i+56]);
				     rhov[7].v = _mm512_fmadd_pd(xv[7],yv[7],rhov[7].v);
				     _mm_prefetch((const char*)&x[i+80],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+80],_MM_HINT_T0);
				     xv[8] = _mm512_load_pd(&x[i+64]);
				     yv[8] = _mm512_load_pd(&y[i+64]);
				     rhov[8].v = _mm512_fmadd_pd(xv[8],yv[8],rhov[8].v);
				     xv[9] = _mm512_load_pd(&x[i+72]);
				     yv[9] = _mm512_load_pd(&y[i+72]);
				     rhov[9].v = _mm512_fmadd_pd(xv[9],yv[9],rhov[9].v);
#else
                                     _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+16],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+48],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&x[i+80],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+80],_MM_HINT_T0);
				     xv[0] = _mm512_load_pd(&x[i+0]);
				     xv[1] = _mm512_load_pd(&x[i+8]);
				     xv[2] = _mm512_load_pd(&x[i+16]);
				     xv[3] = _mm512_load_pd(&x[i+24]);
				     xv[4] = _mm512_load_pd(&x[i+32]);
				     xv[5] = _mm512_load_pd(&x[i+40]);
				     xv[6] = _mm512_load_pd(&x[i+48]);
				     xv[7] = _mm512_load_pd(&x[i+56]);
				     xv[8] = _mm512_load_pd(&x[i+64]);
				     xv[9] = _mm512_load_pd(&x[i+72]);
				     yv[0] = _mm512_load_pd(&y[i+0]);
				     yv[1] = _mm512_load_pd(&y[i+8]);
				     yv[2] = _mm512_load_pd(&y[i+16]);
				     yv[3] = _mm512_load_pd(&y[i+24]);
				     yv[4] = _mm512_load_pd(&y[i+32]);
				     yv[5] = _mm512_load_pd(&y[i+40]);
				     yv[6] = _mm512_load_pd(&y[i+48]);
				     yv[7] = _mm512_load_pd(&y[i+56]);
				     yv[8] = _mm512_load_pd(&y[i+64]);
				     yv[9] = _mm512_load_pd(&y[i+72]);
				     rhov[0].v = _mm256_fmadd_pd(xv[0],yv[0],rhov[0].v);
			             rhov[1].v = _mm256_fmadd_pd(xv[1],yv[1],rhov[1].v);
			             rhov[2].v = _mm256_fmadd_pd(xv[2],yv[2],rhov[2].v);
			             rhov[3].v = _mm256_fmadd_pd(xv[3],yv[3],rhov[3].v);
			             rhov[4].v = _mm256_fmadd_pd(xv[4],yv[4],rhov[4].v);
			             rhov[5].v = _mm256_fmadd_pd(xv[5],yv[5],rhov[5].v);
			             rhov[6].v = _mm256_fmadd_pd(xv[6],yv[6],rhov[6].v);
			             rhov[7].v = _mm256_fmadd_pd(xv[7],yv[7],rhov[7].v);
			             rhov[8].v = _mm256_fmadd_pd(xv[8],yv[8] rhov[8].v);
			             rhov[9].v = _mm256_fmadd_pd(xv[9],yv[9],rhov[9].v);
#endif

				}

				   	rhov[0].v += rhov[5].v;
		                        rhov[1].v += rhov[6].v;
		                        rhov[2].v += rhov[7].v;
		                        rhov[3].v += rhov[8].v;
		                        rhov[4].v += rhov[9].v;

			       for(; (i+39) < n; i += 40) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                     xv[0] = _mm512_load_pd(&x[i+0]);
				     yv[0] = _mm512_load_pd(&y[i+0]);
				     rhov[0].v = _mm512_fmadd_pd(xv[0],yv[0],rhov[0].v);
				     xv[1] = _mm512_load_pd(&x[i+8]);
				     yv[1] = _mm512_load_pd(&y[i+8]);
				     rhov[1].v = _mm512_fmadd_pd(xv[1],yv[1],rhov[1].v);
				     //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     //_mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				     xv[2] = _mm512_load_pd(&x[i+16]);
				     yv[2] = _mm512_load_pd(&y[i+16]);
				     rhov[2].v = _mm512_fmadd_pd(xv[2],yv[2],rhov[2].v);
				     xv[3] = _mm512_load_pd(&x[i+24]);
				     yv[3] = _mm512_load_pd(&y[i+24]);
				     rhov[3].v = _mm512_fmadd_pd(xv[3],yv[3],rhov[3].v);
				     //_mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				     // _mm_prefetch((const char*)&y[i+48],_MM_HINT_T0);
				     xv[4] = _mm512_load_pd(&x[i+32]);
				     yv[4] = _mm512_load_pd(&y[i+32]);
				     rhov[4].v = _mm512_fmadd_pd(xv[4],yv[4],rhov[4].v);
#else
                                     xv[0] = _mm512_load_pd(&x[i+0]);
				     xv[1] = _mm512_load_pd(&x[i+8]);
				     xv[2] = _mm512_load_pd(&x[i+16]);
				     xv[3] = _mm512_load_pd(&x[i+24]);
				     xv[4] = _mm512_load_pd(&x[i+32]);
				     yv[0] = _mm512_load_pd(&y[i+0]);
				     yv[1] = _mm512_load_pd(&y[i+8]);
				     yv[2] = _mm512_load_pd(&y[i+16]);
				     yv[3] = _mm512_load_pd(&y[i+24]);
				     yv[4] = _mm512_load_pd(&y[i+32]);
				     rhov[0].v = _mm256_fmadd_pd(xv[0],yv[0],rhov[0].v);
			             rhov[1].v = _mm256_fmadd_pd(xv[1],yv[1],rhov[1].v);
			             rhov[2].v = _mm256_fmadd_pd(xv[2],yv[2],rhov[2].v);
			             rhov[3].v = _mm256_fmadd_pd(xv[3],yv[3],rhov[3].v);
			             rhov[4].v = _mm256_fmadd_pd(xv[4],yv[4],rhov[4].v);
#endif

				 }

				     rhov[0].v += rhov[2].v;
				     rhov[1].v += rhov[3].v;
				     rhov[0].v += rhov[4].v;

				 for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                        xv[0] = _mm512_load_pd(&x[i+0]);
				        yv[0] = _mm512_load_pd(&y[i+0]);
				        rhov[0].v = _mm512_fmadd_pd(xv[0],yv[0],rhov[0].v);
				        xv[1] = _mm512_load_pd(&x[i+8]);
				        yv[1] = _mm512_load_pd(&y[i+8]);
				        rhov[1].v = _mm512_fmadd_pd(xv[1],yv[1],rhov[1].v);
				     //_mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     //_mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				        xv[2] = _mm512_load_pd(&x[i+16]);
				        yv[2] = _mm512_load_pd(&y[i+16]);
				        rhov[2].v = _mm512_fmadd_pd(xv[2],yv[2],rhov[2].v);
				        xv[3] = _mm512_load_pd(&x[i+24]);
				        yv[3] = _mm512_load_pd(&y[i+24]);
				        rhov[3].v = _mm512_fmadd_pd(xv[3],yv[3],rhov[3].v);
#else
                                        xv[0] = _mm512_load_pd(&x[i+0]);
				        xv[1] = _mm512_load_pd(&x[i+8]);
				        xv[2] = _mm512_load_pd(&x[i+16]);
				        xv[3] = _mm512_load_pd(&x[i+24]);
					yv[0] = _mm512_load_pd(&y[i+0]);
				        yv[1] = _mm512_load_pd(&y[i+8]);
				        yv[2] = _mm512_load_pd(&y[i+16]);
				        yv[3] = _mm512_load_pd(&y[i+24]);
					rhov[0].v = _mm256_fmadd_pd(xv[0],yv[0],rhov[0].v);
			                rhov[1].v = _mm256_fmadd_pd(xv[1],yv[1],rhov[1].v);
			                rhov[2].v = _mm256_fmadd_pd(xv[2],yv[2],rhov[2].v);
			                rhov[3].v = _mm256_fmadd_pd(xv[3],yv[3],rhov[3].v);
#endif

				 }

				     rhov[0].v += rhov[2].v;
				     rhov[1].v += rhov[3].v;

				 for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                        xv[0] = _mm512_load_pd(&x[i+0]);
				        yv[0] = _mm512_load_pd(&y[i+0]);
				        rhov[0].v = _mm512_fmadd_pd(xv[0],yv[0],rhov[0].v);
				        xv[1] = _mm512_load_pd(&x[i+8]);
				        yv[1] = _mm512_load_pd(&y[i+8]);
				        rhov[1].v = _mm512_fmadd_pd(xv[1],yv[1],rhov[1].v);
#else
                                        xv[0] = _mm512_load_pd(&x[i+0]);
				        xv[1] = _mm512_load_pd(&x[i+8]);
					yv[0] = _mm512_load_pd(&y[i+0]);
				        yv[1] = _mm512_load_pd(&y[i+8]);
					rhov[0].v = _mm256_fmadd_pd(xv[0],yv[0],rhov[0].v);
			                rhov[1].v = _mm256_fmadd_pd(xv[1],yv[1],rhov[1].v);
#endif
				 }
				     
                                       rhov[0].v += rhov[1].v;
			  

			        for(; (i+7) < n; i += 8) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                        xv[0] = _mm512_load_pd(&x[i+0]);
				        yv[0] = _mm512_load_pd(&y[i+0]);
				        rhov[0].v = _mm512_fmadd_pd(xv[0],yv[0],rhov[0].v);
#else
                                         xv[0] = _mm512_load_pd(&x[i+0]);
					 yv[0] = _mm512_load_pd(&y[i+0]);
					 rhov[0].v = _mm256_fmadd_pd(xv[0],yv[0],rhov[0].v);	
#endif
				}

				for(; (i+0) < n; i += 1) {
                                      rho0 += x[i] * y[i];
				}

				  rho0 += _mm512_reduce_add_pd(rhov[0].v);

			  }
			  else {

			           	for (i = 0; i < n; ++i) {
					     const double x0c = *x0;
			                     const double y0c = *y0;
			                     rho0 += x0c * y0c;
			                     x0 += incx;
			                     y0 += incy;
				      }
			  }

			  *rho = rho0;

		   }



		   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
		   void ddotv_a_zmm8r8_unroll10x_omp(const int32_t,
		                                     double * __restrict __ATTR_ALIGN__(64) x,
						     const int32_t incx,
						     double * __restrict __ATTR_ALIGN__(64)y,
						     const int32_t incy,
						     double * __restrict rho) {

                             if(__builtin_expect(0==n,0)) { return;}
			   
			   __ATTR_ALIGN__(64) __m512d    xv[10];
			   __ATTR_ALIGN__(64) __m512d    yv[10];
			  // __ATTR_ALIGN__(64) zmm8r8_t rhov[10];
			   zmm8r8_t rhov0;
			   zmm8r8_t rhov1;
			   zmm8r8_t rhov2;
			   zmm8r8_t rhov3;
			   zmm8r8_t rhov4;
			   zmm8r8_t rhov5;
			   zmm8r8_t rhov6;
			   zmm8r8_t rhov7;
			   zmm8r8_t rhov8;
			   zmm8r8_t rhov9;
			   const __m512d vz = _mm512_setzero_pd();
			   double rho0;
			   int32_t i;

			   int32_t ii;
			   rho0 = 0.0;
			   ii = 0;
			   if(__builtin_expect(1==incx,1) &&
			      __builtin_expect(1==incy,1)) {

			        rhov0.v = vz;
				rhov1.v = vz;
				rhov2.v = vz;
				rhov3.v = vz;
				rhov4.v = vz;
				rhov5.v = vz;
				rhov6.v = vz;
				rhov7.v = vz;
				rhov8.v = vz;
				rhov9.v = vz;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                             __assume_aligned(x,64);
			     __assume_aligned(y,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
                             x = (double*)__builtin_assume_aligned(x,64);
			     y = (double*)__builtin_assume_aligned(y,64);
#endif
#pragma omp parallel for schedule(static,80) default(none) private(i) shared(x,y,n) \
            lastprivate(ii,rhov0,rhov1,rhov2,rhov3,rhov4,rhov5,rhov6,rhov7,rhov8,rhov9) 
				for(i = 0; (i+79) < n; i += 80) {
				      ii = i;

                                     _mm_prefetch((const char*)&x[i+16],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+16],_MM_HINT_T0);
				     rhov0.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+0]),
						                _mm512_load_pd(&y[i+0]),
								rhov0.v);
				      rhov1.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+8]),
						                _mm512_load_pd(&y[i+8]),
								rhov1.v);
				     _mm_prefetch((const char*)&x[i+32],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+32],_MM_HINT_T0);
				      rhov2.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+16]),
						                _mm512_load_pd(&y[i+16]),
								rhov2.v);
				      rhov3.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+24]),
						                _mm512_load_pd(&y[i+24]),
								rhov3.v);
				     _mm_prefetch((const char*)&x[i+48],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+48],_MM_HINT_T0);
				      rhov4.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+32]),
						                _mm512_load_pd(&y[i+32]),
								rhov4.v);
				      rhov5.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+40]),
						                _mm512_load_pd(&y[i+40]),
								rhov5.v);
				     _mm_prefetch((const char*)&x[i+64],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+64],_MM_HINT_T0);
				      rhov6.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+48]),
						                _mm512_load_pd(&y[i+48]),
								rhov6.v);
				      rhov7.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+56]),
						                _mm512_load_pd(&y[i+56]),
								rhov7.v);
				     _mm_prefetch((const char*)&x[i+80],_MM_HINT_T0);
				     _mm_prefetch((const char*)&y[i+80],_MM_HINT_T0);
				     rhov8.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+64]),
						                _mm512_load_pd(&y[i+64]),
								rhov8.v);
				     rhov9.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+72]),
						                _mm512_load_pd(&y[i+72]),
								rhov9.v);
				   

			       }

				   	rhov[0].v += rhov[5].v;
		                        rhov[1].v += rhov[6].v;
		                        rhov[2].v += rhov[7].v;
		                        rhov[3].v += rhov[8].v;
		                        rhov[4].v += rhov[9].v;

			       for(; (ii+39) < n; ii += 40) {
                                      rhov0.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+0]),
						                _mm512_load_pd(&y[i+0]),
								rhov0.v);
				      rhov1.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+8]),
						                _mm512_load_pd(&y[i+8]),
								rhov1.v);
				    
				      rhov2.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+16]),
						                _mm512_load_pd(&y[i+16]),
								rhov2.v);
				      rhov3.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+24]),
						                _mm512_load_pd(&y[i+24]),
								rhov3.v);
				    
				      rhov4.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+32]),
						                _mm512_load_pd(&y[i+32]),
								rhov4.v);
                                  

				 }

				     rhov[0].v += rhov[2].v;
				     rhov[1].v += rhov[3].v;
				     rhov[0].v += rhov[4].v;

				 for(; (ii+31) < n; ii += 32) {

				      rhov0.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+0]),
						                _mm512_load_pd(&y[i+0]),
								rhov0.v);
				      rhov1.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+8]),
						                _mm512_load_pd(&y[i+8]),
								rhov1.v);
				    
				      rhov2.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+16]),
						                _mm512_load_pd(&y[i+16]),
								rhov2.v);
				      rhov3.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+24]),
						                _mm512_load_pd(&y[i+24]),
								rhov3.v);
                                       

				 }

				     rhov[0].v += rhov[2].v;
				     rhov[1].v += rhov[3].v;

				 for(; (ii+15) < n; ii += 16) {
                                      rhov0.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+0]),
						                _mm512_load_pd(&y[i+0]),
								rhov0.v);
				      rhov1.v =
				                _mm512_fmadd_pd(_mm512_load_pd(&x[i+8]),
						                _mm512_load_pd(&y[i+8]),
								rhov1.v);
				    
                                     

				 }
				     
                                       rhov[0].v += rhov[1].v;
			  

			        for(; (ii+7) < n; ii += 8) {

                                    rhov0.v =
				              _mm512_fmadd_pd(_mm512_load_pd(&x[i+0]),
						                _mm512_load_pd(&y[i+0]),
								rhov0.v);

				}

				for(; (ii+0) < n; ii += 1) {
                                      rho0 += x[ii] * y[ii];
				}

				  rho0 += _mm512_reduce_add_pd(rhov[0].v);

			  }
			  else {

			           	for (i = 0; i < n; ++i) {
					     const double x0c = *x0;
			                     const double y0c = *y0;
			                     rho0 += x0c * y0c;
			                     x0 += incx;
			                     y0 += incy;
				      }
			  }

			  *rho = rho0;

		   }







		   










		   
    } // math



} // gms


























#endif /*__GMS_DOTV_AVX512_UNROLL10X_HPP__*/
