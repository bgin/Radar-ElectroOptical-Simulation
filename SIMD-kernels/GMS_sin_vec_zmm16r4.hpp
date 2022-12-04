

#ifndef __GMS_SIN_VEC_ZMM16R4_HPP__
#define __GMS_SIN_VEC_ZMM16R4_HPP__ 041220220149

/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

namespace file_version {

    const unsigned int GMS_SIN_VEC_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_SIN_VEC_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_SIN_VEC_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_SIN_VEC_ZMM16R4_FULLVER =
      1000U*GMS_SIN_VEC_ZMM16R4_MAJOR+
      100U*GMS_SIN_VEC_ZMM16R4_MINOR+
      10U*GMS_SIN_VEC_ZMM16R4_MICRO;
    const char * const GMS_SIN_VEC_ZMM16R4_CREATION_DATE = "04-12-2022 10:49 AM +00200 (SUN 04 DEC 2022 GMT+2)";
    const char * const GMS_SIN_VEC_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_SIN_VEC_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_SIN_VEC_ZMM16R4_DESCRIPTION   = "AVX512 optimized vector of sin values."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_cephes.h"

namespace  gms {


         namespace  math {

/*
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_sin_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void sinv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                                   float * __restrict __ATTR_ALIGN__(64) y,
                                                   const __m512 a,
                                                   const __m512 b,
                                                   const __m512 c,
                                                   const int32_t n) {

                       if(__builtin_expect(0==n,0)) {return;}
                       volatile __m512 d;
                       int32_t i;
                       // Start the preload phase.
                       _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                       volatile __m512 first = _mm512_sin_ps(_mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

                       // Main processing loop starts.
                       for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           register const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                           register const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                           register const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                           register const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                           register const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                           register const __m512 zmm9      = _mm512_sin_ps(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           register const __m512 zmm10      = _mm512_load_ps(&x[i+80]); 
                           register const __m512 zmm11      = _mm512_sin_ps(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
                           register const __m512 zmm12      = _mm512_load_ps(&x[i+96]); 
                           register const __m512 zmm13      = _mm512_sin_ps(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
                           register const __m512 zmm14      = _mm512_load_ps(&x[i+112]); 
                           register const __m512 zmm15      = _mm512_sin_ps(zmm14);
                           _mm512_store_ps(&y[i+112], zmm15);
                           register const __m512 zmm16      = _mm512_load_ps(&x[i+128]); 
                           register const __m512 zmm17      = _mm512_sin_ps(zmm16);
                           _mm512_store_ps(&y[i+128], zmm17);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           register const __m512 zmm18      = _mm512_load_ps(&x[i+144]); 
                           register const __m512 zmm19      = _mm512_sin_ps(zmm18);
                           _mm512_store_ps(&y[i+144], zmm19);
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           register const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                           register const __m512 zmm10     = _mm512_load_ps(&x[i+80]);
                           register const __m512 zmm12     = _mm512_load_ps(&x[i+96]); 
                           register const __m512 zmm14     = _mm512_load_ps(&x[i+112]);  
                           register const __m512 zmm16     = _mm512_load_ps(&x[i+128]); 
                           register const __m512 zmm18     = _mm512_load_ps(&x[i+144]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           register const __m512 zmm9      = _mm512_sin_ps(zmm8);
                           register const __m512 zmm11     = _mm512_sin_ps(zmm10);
                           register const __m512 zmm13     = _mm512_sin_ps(zmm12);
                           register const __m512 zmm15     = _mm512_sin_ps(zmm14);
                           register const __m512 zmm17     = _mm512_sin_ps(zmm16);
                           register const __m512 zmm19     = _mm512_sin_ps(zmm18);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
                           _mm512_store_ps(&y[i+96], zmm13);
                           _mm512_store_ps(&y[i+112], zmm15);
                           _mm512_store_ps(&y[i+128], zmm17);
                           _mm512_store_ps(&y[i+144], zmm19);
#endif
                       }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                           register const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                           register const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                           register const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                           register const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                           register const __m512 zmm9      = _mm512_sin_ps(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           register const __m512 zmm10      = _mm512_load_ps(&x[i+80]); 
                           register const __m512 zmm11      = _mm512_sin_ps(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
                           register const __m512 zmm12      = _mm512_load_ps(&x[i+96]); 
                           register const __m512 zmm13      = _mm512_sin_ps(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
                           register const __m512 zmm14      = _mm512_load_ps(&x[i+112]); 
                           register const __m512 zmm15      = _mm512_sin_ps(zmm14);
                           _mm512_store_ps(&y[i+112], zmm15);
#else
                           register const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                           register const __m512 zmm10     = _mm512_load_ps(&x[i+80]);
                           register const __m512 zmm12     = _mm512_load_ps(&x[i+96]); 
                           register const __m512 zmm14     = _mm512_load_ps(&x[i+112]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           register const __m512 zmm9      = _mm512_sin_ps(zmm8);
                           register const __m512 zmm11     = _mm512_sin_ps(zmm10);
                           register const __m512 zmm13     = _mm512_sin_ps(zmm12);
                           register const __m512 zmm15     = _mm512_sin_ps(zmm14);
                           register const __m512 zmm17     = _mm512_sin_ps(zmm16);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
                           _mm512_store_ps(&y[i+96], zmm13);
                           _mm512_store_ps(&y[i+112], zmm15);
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                           register const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                           register const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                           register const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#else
                           register const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                           register const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#else
                           register const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
#else
                           register const __m512 zmm0      = _mm512_load_ps(&x[i+0]);
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0); 
                           _mm512_store_ps(&y[i+0], zmm1);
#endif
                     } 

                     for(; (i+0) < n; i += 1) {
                           y[i] = ceph_sinf(x[i]);
                     }
               }


              /*
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_sin_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void sinv_zmm16r4_unroll_10x_u(const float * __restrict  x,
                                                  float * __restrict  y,
                                                  const __m512 a,
                                                  const __m512 b,
                                                  const __m512 c,
                                                  const int32_t n) {

                       if(__builtin_expect(0==n,0)) {return;}
                       volatile __m512 d;
                       int32_t i;
                       // Start the preload phase.
                       _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                       volatile __m512 first = _mm512_sin_ps(_mm512_loadu_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

                       // Main processing loop starts.
                       for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           register const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           register const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           register const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           register const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           register const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                           register const __m512 zmm9      = _mm512_sin_ps(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           register const __m512 zmm10      = _mm512_loadu_ps(&x[i+80]); 
                           register const __m512 zmm11      = _mm512_sin_ps(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           register const __m512 zmm12      = _mm512_loadu_ps(&x[i+96]); 
                           register const __m512 zmm13      = _mm512_sin_ps(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           register const __m512 zmm14      = _mm512_loadu_ps(&x[i+112]); 
                           register const __m512 zmm15      = _mm512_sin_ps(zmm14);
                           _mm512_storeu_ps(&y[i+112], zmm15);
                           register const __m512 zmm16      = _mm512_loadu_ps(&x[i+128]); 
                           register const __m512 zmm17      = _mm512_sin_ps(zmm16);
                           _mm512_storeu_ps(&y[i+128], zmm17);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           register const __m512 zmm18      = _mm512_loadu_ps(&x[i+144]); 
                           register const __m512 zmm19      = _mm512_sin_ps(zmm18);
                           _mm512_storeu_ps(&y[i+144], zmm19);
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           register const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                           register const __m512 zmm10     = _mm512_loadu_ps(&x[i+80]);
                           register const __m512 zmm12     = _mm512_loadu_ps(&x[i+96]); 
                           register const __m512 zmm14     = _mm512_loadu_ps(&x[i+112]);  
                           register const __m512 zmm16     = _mm512_loadu_ps(&x[i+128]); 
                           register const __m512 zmm18     = _mm512_loadu_ps(&x[i+144]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           register const __m512 zmm9      = _mm512_sin_ps(zmm8);
                           register const __m512 zmm11     = _mm512_sin_ps(zmm10);
                           register const __m512 zmm13     = _mm512_sin_ps(zmm12);
                           register const __m512 zmm15     = _mm512_sin_ps(zmm14);
                           register const __m512 zmm17     = _mm512_sin_ps(zmm16);
                           register const __m512 zmm19     = _mm512_sin_ps(zmm18);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           _mm512_storeu_ps(&y[i+112], zmm15);
                           _mm512_storeu_ps(&y[i+128], zmm17);
                           _mm512_storeu_ps(&y[i+144], zmm19);
#endif
                       }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           register const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           register const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           register const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           register const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                           register const __m512 zmm9      = _mm512_sin_ps(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           register const __m512 zmm10      = _mm512_loadu_ps(&x[i+80]); 
                           register const __m512 zmm11      = _mm512_sin_ps(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           register const __m512 zmm12      = _mm512_loadu_ps(&x[i+96]); 
                           register const __m512 zmm13      = _mm512_sin_ps(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           register const __m512 zmm14      = _mm512_loadu_ps(&x[i+112]); 
                           register const __m512 zmm15      = _mm512_sin_ps(zmm14);
                           _mm512_storeu_ps(&y[i+112], zmm15);
#else
                           register const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                           register const __m512 zmm10     = _mm512_loadu_ps(&x[i+80]);
                           register const __m512 zmm12     = _mm512_loadu_ps(&x[i+96]); 
                           register const __m512 zmm14     = _mm512_loadu_ps(&x[i+112]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           register const __m512 zmm9      = _mm512_sin_ps(zmm8);
                           register const __m512 zmm11     = _mm512_sin_ps(zmm10);
                           register const __m512 zmm13     = _mm512_sin_ps(zmm12);
                           register const __m512 zmm15     = _mm512_sin_ps(zmm14);
                           register const __m512 zmm17     = _mm512_sin_ps(zmm16);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           _mm512_storeu_ps(&y[i+112], zmm15);
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           register const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           register const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           register const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#else
                           register const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           register const __m512 zmm5      = _mm512_sin_ps(zmm4);
                           register const __m512 zmm7      = _mm512_sin_ps(zmm6);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           register const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#else
                           register const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           register const __m512 zmm3      = _mm512_sin_ps(zmm2);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
#else
                           register const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]);
                           register const __m512 zmm1      = _mm512_sin_ps(zmm0); 
                           _mm512_storeu_ps(&y[i+0], zmm1);
#endif
                     } 

                     for(; (i+0) < n; i += 1) {
                           y[i] = ceph_sinf(x[i]);
                     }
               }

               
               /*
                    Calls non-SVML implementation of sine function
                    SLEEF version is inlined.
                */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void sinv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                                  float * __restrict __ATTR_ALIGN__(64) y,
                                                  const int32_t n) {

                       if(__builtin_expect(0==n,0)) {return;}
                       int32_t i;
                      
                       for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           register const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                           register const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const vfloat zmm3      = xsinf(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                           register const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const vfloat zmm5      = xsinf(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                           register const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const vfloat zmm7      = xsinf(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                           register const vfloat zmm8      = _mm512_load_ps(&x[i+64]); 
                           register const vfloat zmm9      = xsinf(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           register const vfloat zmm10      = _mm512_load_ps(&x[i+80]); 
                           register const vfloat zmm11      = xsinf(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
                           register const vfloat zmm12      = _mm512_load_ps(&x[i+96]); 
                           register const vfloat zmm13      = xsinf(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
                           register const vfloat zmm14      = _mm512_load_ps(&x[i+112]); 
                           register const vfloat zmm15      = xsinf(zmm14);
                           _mm512_store_ps(&y[i+112], zmm15);
                           register const vfloat zmm16      = _mm512_load_ps(&x[i+128]); 
                           register const vfloat zmm17      = xsinf(zmm16);
                           _mm512_store_ps(&y[i+128], zmm17);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           register const vfloat zmm18      = _mm512_load_ps(&x[i+144]); 
                           register const vfloat zmm19      = xsinf(zmm18);
                           _mm512_store_ps(&y[i+144], zmm19);
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           register const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const vfloat zmm8      = _mm512_load_ps(&x[i+64]); 
                           register const vfloat zmm10     = _mm512_load_ps(&x[i+80]);
                           register const vfloat zmm12     = _mm512_load_ps(&x[i+96]); 
                           register const vfloat zmm14     = _mm512_load_ps(&x[i+112]);  
                           register const vfloat zmm16     = _mm512_load_ps(&x[i+128]); 
                           register const vfloat zmm18     = _mm512_load_ps(&x[i+144]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           register const vfloat zmm3      = xsinf(zmm2);
                           register const vfloat zmm5      = xsinf(zmm4);
                           register const vfloat zmm7      = xsinf(zmm6);
                           register const vfloat zmm9      = xsinf(zmm8);
                           register const vfloat zmm11     = xsinf(zmm10);
                           register const vfloat zmm13     = xsinf(zmm12);
                           register const vfloat zmm15     = xsinf(zmm14);
                           register const vfloat zmm17     = xsinf(zmm16);
                           register const vfloat zmm19     = xsinf(zmm18);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
                           _mm512_store_ps(&y[i+96], zmm13);
                           _mm512_store_ps(&y[i+112], zmm15);
                           _mm512_store_ps(&y[i+128], zmm17);
                           _mm512_store_ps(&y[i+144], zmm19);
#endif
                       }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                           register const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const vfloat zmm3      = xsinf(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                           register const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const vfloat zmm5      = xsinf(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                           register const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const vfloat zmm7      = xsinf(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                           register const vfloat zmm8      = _mm512_load_ps(&x[i+64]); 
                           register const vfloat zmm9      = xsinf(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           register const vfloat zmm10      = _mm512_load_ps(&x[i+80]); 
                           register const vfloat zmm11      = xsinf(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
                           register const vfloat zmm12      = _mm512_load_ps(&x[i+96]); 
                           register const vfloat zmm13      = xsinf(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
                           register const vfloat zmm14      = _mm512_load_ps(&x[i+112]); 
                           register const vfloat zmm15      = xsinf(zmm14);
                           _mm512_store_ps(&y[i+112], zmm15);
#else
                           register const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const vfloat zmm8      = _mm512_load_ps(&x[i+64]); 
                           register const vfloat zmm10     = _mm512_load_ps(&x[i+80]);
                           register const vfloat zmm12     = _mm512_load_ps(&x[i+96]); 
                           register const vfloat zmm14     = _mm512_load_ps(&x[i+112]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           register const vfloat zmm3      = xsinf(zmm2);
                           register const vfloat zmm5      = xsinf(zmm4);
                           register const vfloat zmm7      = xsinf(zmm6);
                           register const vfloat zmm9      = xsinf(zmm8);
                           register const vfloat zmm11     = xsinf(zmm10);
                           register const vfloat zmm13     = xsinf(zmm12);
                           register const vfloat zmm15     = xsinf(zmm14);
                           register const vfloat zmm17     = xsinf(zmm16);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
                           _mm512_store_ps(&y[i+96], zmm13);
                           _mm512_store_ps(&y[i+112], zmm15);
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                           register const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const vfloat zmm3      = xsinf(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                           register const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const vfloat zmm5      = xsinf(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                           register const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const vfloat zmm7      = xsinf(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#else
                           register const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const vfloat zmm4      = _mm512_load_ps(&x[i+32]); 
                           register const vfloat zmm6      = _mm512_load_ps(&x[i+48]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           register const vfloat zmm3      = xsinf(zmm2);
                           register const vfloat zmm5      = xsinf(zmm4);
                           register const vfloat zmm7      = xsinf(zmm6);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                           register const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const vfloat zmm3      = xsinf(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#else
                           register const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const vfloat zmm2      = _mm512_load_ps(&x[i+16]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           register const vfloat zmm3      = xsinf(zmm2);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const vfloat zmm0      = _mm512_load_ps(&x[i+0]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
#else
                           register const vfloat zmm0      = _mm512_load_ps(&x[i+0]);
                           register const vfloat zmm1      = xsinf(zmm0); 
                           _mm512_store_ps(&y[i+0], zmm1);
#endif
                     } 

                     for(; (i+0) < n; i += 1) {
                           y[i] = ceph_sinf(x[i]);
                     }
               } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
	           static inline
                   void sinv_zmm16r4_unroll_10x_u(const float * __restrict  x,
                                                  float * __restrict  y,
                                                  const int32_t n) {

                       if(__builtin_expect(0==n,0)) {return;}
                       int32_t i;
                      
                       for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           register const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           register const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const vfloat zmm3      = xsinf(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           register const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const vfloat zmm5      = xsinf(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           register const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const vfloat zmm7      = xsinf(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           register const vfloat zmm8      = _mm512_loadu_ps(&x[i+64]); 
                           register const vfloat zmm9      = xsinf(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           register const vfloat zmm10      = _mm512_loadu_ps(&x[i+80]); 
                           register const vfloat zmm11      = xsinf(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           register const vfloat zmm12      = _mm512_loadu_ps(&x[i+96]); 
                           register const vfloat zmm13      = xsinf(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           register const vfloat zmm14      = _mm512_loadu_ps(&x[i+112]); 
                           register const vfloat zmm15      = xsinf(zmm14);
                           _mm512_storeu_ps(&y[i+112], zmm15);
                           register const vfloat zmm16      = _mm512_loadu_ps(&x[i+128]); 
                           register const vfloat zmm17      = xsinf(zmm16);
                           _mm512_storeu_ps(&y[i+128], zmm17);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           register const vfloat zmm18      = _mm512_loadu_ps(&x[i+144]); 
                           register const vfloat zmm19      = xsinf(zmm18);
                           _mm512_storeu_ps(&y[i+144], zmm19);
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           register const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const vfloat zmm8      = _mm512_loadu_ps(&x[i+64]); 
                           register const vfloat zmm10     = _mm512_loadu_ps(&x[i+80]);
                           register const vfloat zmm12     = _mm512_loadu_ps(&x[i+96]); 
                           register const vfloat zmm14     = _mm512_loadu_ps(&x[i+112]);  
                           register const vfloat zmm16     = _mm512_loadu_ps(&x[i+128]); 
                           register const vfloat zmm18     = _mm512_loadu_ps(&x[i+144]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           register const vfloat zmm3      = xsinf(zmm2);
                           register const vfloat zmm5      = xsinf(zmm4);
                           register const vfloat zmm7      = xsinf(zmm6);
                           register const vfloat zmm9      = xsinf(zmm8);
                           register const vfloat zmm11     = xsinf(zmm10);
                           register const vfloat zmm13     = xsinf(zmm12);
                           register const vfloat zmm15     = xsinf(zmm14);
                           register const vfloat zmm17     = xsinf(zmm16);
                           register const vfloat zmm19     = xsinf(zmm18);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           _mm512_storeu_ps(&y[i+112], zmm15);
                           _mm512_storeu_ps(&y[i+128], zmm17);
                           _mm512_storeu_ps(&y[i+144], zmm19);
#endif
                       }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           register const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const vfloat zmm3      = xsinf(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           register const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const vfloat zmm5      = xsinf(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           register const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const vfloat zmm7      = xsinf(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           register const vfloat zmm8      = _mm512_loadu_ps(&x[i+64]); 
                           register const vfloat zmm9      = xsinf(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           register const vfloat zmm10      = _mm512_load_ps(&x[i+80]); 
                           register const vfloat zmm11      = xsinf(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           register const vfloat zmm12      = _mm512_loadu_ps(&x[i+96]); 
                           register const vfloat zmm13      = xsinf(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           register const vfloat zmm14      = _mm512_loadu_ps(&x[i+112]); 
                           register const vfloat zmm15      = xsinf(zmm14);
                           _mm512_storeu_ps(&y[i+112], zmm15);
#else
                           register const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const vfloat zmm8      = _mm512_loadu_ps(&x[i+64]); 
                           register const vfloat zmm10     = _mm512_loadu_ps(&x[i+80]);
                           register const vfloat zmm12     = _mm512_loadu_ps(&x[i+96]); 
                           register const vfloat zmm14     = _mm512_loadu_ps(&x[i+112]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           register const vfloat zmm3      = xsinf(zmm2);
                           register const vfloat zmm5      = xsinf(zmm4);
                           register const vfloat zmm7      = xsinf(zmm6);
                           register const vfloat zmm9      = xsinf(zmm8);
                           register const vfloat zmm11     = xsinf(zmm10);
                           register const vfloat zmm13     = xsinf(zmm12);
                           register const vfloat zmm15     = xsinf(zmm14);
                           register const vfloat zmm17     = xsinf(zmm16);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                           _mm512_storeu_ps(&y[i+112], zmm15);
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           register const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const vfloat zmm3      = xsinf(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           register const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const vfloat zmm5      = xsinf(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           register const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const vfloat zmm7      = xsinf(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#else
                           register const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const vfloat zmm4      = _mm512_loadu_ps(&x[i+32]); 
                           register const vfloat zmm6      = _mm512_loadu_ps(&x[i+48]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           register const vfloat zmm3      = xsinf(zmm2);
                           register const vfloat zmm5      = xsinf(zmm4);
                           register const vfloat zmm7      = xsinf(zmm6);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           register const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const vfloat zmm3      = xsinf(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#else
                           register const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const vfloat zmm2      = _mm512_loadu_ps(&x[i+16]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           register const vfloat zmm3      = xsinf(zmm2);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           register const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]); 
                           register const vfloat zmm1      = xsinf(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
#else
                           register const vfloat zmm0      = _mm512_loadu_ps(&x[i+0]);
                           register const vfloat zmm1      = xsinf(zmm0); 
                           _mm512_storeu_ps(&y[i+0], zmm1);
#endif
                     } 

                     for(; (i+0) < n; i += 1) {
                           y[i] = ceph_sinf(x[i]);
                     }
               } 

                    

        } // math 
 

} // gms

















#endif /*__GMS_SIN_VEC_ZMM16R4_HPP__*/
