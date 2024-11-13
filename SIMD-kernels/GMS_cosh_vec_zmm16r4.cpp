


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




#include <immintrin.h>
#include "GMS_cosh_vec_zmm16r4.h"
#include "GMS_cephes.h"


/*
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cosh_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/
                 
                   void gms::math::coshv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) x,
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
                       volatile __m512 first = _mm512_cosh_ps(_mm512_load_ps(&x[0])); //L1I cache miss.
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
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const __m512 zmm10      = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cosh_ps(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
                            const __m512 zmm12      = _mm512_load_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_cosh_ps(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
                            const __m512 zmm14      = _mm512_load_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_cosh_ps(zmm14);
                           _mm512_store_ps(&y[i+112], zmm15);
                            const __m512 zmm16      = _mm512_load_ps(&x[i+128]); 
                            const __m512 zmm17      = _mm512_cosh_ps(zmm16);
                           _mm512_store_ps(&y[i+128], zmm17);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                            const __m512 zmm18      = _mm512_load_ps(&x[i+144]); 
                            const __m512 zmm19      = _mm512_cosh_ps(zmm18);
                           _mm512_store_ps(&y[i+144], zmm19);
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm10     = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm12     = _mm512_load_ps(&x[i+96]); 
                            const __m512 zmm14     = _mm512_load_ps(&x[i+112]);  
                            const __m512 zmm16     = _mm512_load_ps(&x[i+128]); 
                            const __m512 zmm18     = _mm512_load_ps(&x[i+144]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                            const __m512 zmm11     = _mm512_cosh_ps(zmm10);
                            const __m512 zmm13     = _mm512_cosh_ps(zmm12);
                            const __m512 zmm15     = _mm512_cosh_ps(zmm14);
                            const __m512 zmm17     = _mm512_cosh_ps(zmm16);
                            const __m512 zmm19     = _mm512_cosh_ps(zmm18);
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
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const __m512 zmm10      = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cosh_ps(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
                            const __m512 zmm12      = _mm512_load_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_cosh_ps(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
                            const __m512 zmm14      = _mm512_load_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_cosh_ps(zmm14);
                           _mm512_store_ps(&y[i+112], zmm15);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm10     = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm12     = _mm512_load_ps(&x[i+96]); 
                            const __m512 zmm14     = _mm512_load_ps(&x[i+112]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                            const __m512 zmm11     = _mm512_cosh_ps(zmm10);
                            const __m512 zmm13     = _mm512_cosh_ps(zmm12);
                            const __m512 zmm15     = _mm512_cosh_ps(zmm14);
                            const __m512 zmm17     = _mm512_cosh_ps(zmm16);
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
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0); 
                           _mm512_store_ps(&y[i+0], zmm1);
#endif
                     } 

                     for(; (i+0) < n; i += 1) {
                           y[i] = cephes_coshf(x[i]);
                     }
               }


                  
                   void gms::math::coshv_mask_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) x, 
                                                       const float * __restrict __ATTR_ALIGN__(64) z,
                                                       float * __restrict __ATTR_ALIGN__(64) y,
                                                       const __mmask16 * __restrict __ATTR_ALIGN__(64) m,
                                                       const __m512 a,
                                                       const __m512 b,
                                                       const __m512 c,
                                                       const int32_t n,
                                                       const bool additive) {

                     if(__builtin_expect(0==n,0)) {return;}
                     volatile __m512 d;
                     const __
                     int32_t i,j;
                       // Start the preload phase.
                     _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                     _mm_prefetch((const char *)&z[0],_MM_HINT_T0);
                     _mm_prefetch((const char *)&m[0],_MM_HINT_T0);
                     volatile __m512 first = _mm512_mask_cosh_ps(_mm512_load_ps(&z[0],m[0],
                                                                        _mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                     for(int32_t j=0; j != 14; ++j) {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                     }
                     j = 0;
                     if(additive) {

                          // Main processing loop starts.
                        for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+64],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                           _mm512_store_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(_mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_store_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                           _mm512_store_ps(&y[i+64],zmm14);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+128],_MM_HINT_T0);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                           _mm512_store_ps(&y[i+80],zmm17);
                            const __m512 zmm18  = _mm512_load_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_load_ps(&z[i+96]);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm18,zmm19),m[j+6],zmm18);
                           _mm512_store_ps(&y[i+96],zmm20);
                            const __m512 zmm21  = _mm512_load_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_load_ps(&z[i+112]);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm21,zmm22),m[j+7],zmm21);
                           _mm512_store_ps(&y[i+112],zmm23);
                            const __m512 zmm24  = _mm512_load_ps(&x[i+128]);
                            const __m512 zmm25  = _mm512_load_ps(&z[i+128]);
                            const __m512 zmm26  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm24,zmm25),m[j+8],zmm24);
                           _mm512_store_ps(&y[i+128],zmm26;
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+144],_MM_HINT_T0); 
                            const __m512 zmm27  = _mm512_load_ps(&x[i+144]);
                            const __m512 zmm28  = _mm512_load_ps(&z[i+144]);
                            const __m512 zmm29  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm27,zmm28),m[j+9],zmm27);
                           _mm512_store_ps(&y[i+144],zmm29;
                           j += 16;
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+144],_MM_HINT_T0); 
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm18  = _mm512_load_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_load_ps(&z[i+96]);
                            const __m512 zmm21  = _mm512_load_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_load_ps(&z[i+112]);
                            const __m512 zmm24  = _mm512_load_ps(&x[i+128]);
                            const __m512 zmm25  = _mm512_load_ps(&z[i+128]);
                            const __m512 zmm27  = _mm512_load_ps(&x[i+144]);
                            const __m512 zmm28  = _mm512_load_ps(&z[i+144]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm0,zmm1),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm3,zmm4),m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                            const __m512 zmm14 = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm18,zmm19),m[j+6],zmm18);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm21,zmm22),m[j+7],zmm21);
                            const __m512 zmm26  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm24,zmm25),m[j+8],zmm24);
                            const __m512 zmm29  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm27,zmm28),m[j+9],zmm27);
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm512_store_ps(&y[i+48],zmm11);
                           _mm512_store_ps(&y[i+64],zmm14);
                           _mm512_store_ps(&y[i+80],zmm17); 
                           _mm512_store_ps(&y[i+96],zmm20);
                           _mm512_store_ps(&y[i+112],zmm23);
                           _mm512_store_ps(&y[i+128],zmm26);
                           _mm512_store_ps(&y[i+144],zmm29);
                           j += 16;
#endif
                       }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                           _mm512_store_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(_mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_store_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                           _mm512_store_ps(&y[i+64],zmm14);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                           _mm512_store_ps(&y[i+80],zmm17);
                            const __m512 zmm18  = _mm512_load_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_load_ps(&z[i+96]);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm18,zmm19),m[j+6],zmm18);
                           _mm512_store_ps(&y[i+96],zmm20);
                            const __m512 zmm21  = _mm512_load_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_load_ps(&z[i+112]);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm21,zmm22),m[j+7],zmm21);
                           _mm512_store_ps(&y[i+112],zmm23);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm18  = _mm512_load_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_load_ps(&z[i+96]);
                            const __m512 zmm21  = _mm512_load_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_load_ps(&z[i+112]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(
                                                                _mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(
                                                                _mm512_add_ps(zmm18,zmm19),m[j+6],zmm18);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(
                                                                _mm512_add_ps(zmm21,zmm22),m[j+7],zmm21);
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm512_store_ps(&y[i+48],zmm11);
                           _mm512_store_ps(&y[i+64],zmm14);
                           _mm512_store_ps(&y[i+80],zmm17); 
                           _mm512_store_ps(&y[i+96],zmm20);
                           _mm512_store_ps(&y[i+112],zmm23);
                           j += 16;
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                           _mm512_store_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(_mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_store_ps(&y[i+48],zmm11);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm512_store_ps(&y[i+48],zmm11);
                           j += 16;
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           j += 16;
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                           j += 16;
#endif
                     } 

                     for(; (i+0) < n; i += 1) {
                           if(m[i])
                              y[i] = cephes_coshf(x[i]+z[i]);
                           else
                              y[i] = cephes_coshf(x[i]);
                     }


                       }
                        else {
                       // Main processing loop starts.
                          for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+64],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                           _mm512_store_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_store_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                           _mm512_store_ps(&y[i+64],zmm14);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+128],_MM_HINT_T0);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(zmm16,m[j+5],zmm15);
                           _mm512_store_ps(&y[i+80],zmm17);
                            const __m512 zmm18  = _mm512_load_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_load_ps(&z[i+96]);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(zmm19,m[j+6],zmm18);
                           _mm512_store_ps(&y[i+96],zmm20);
                            const __m512 zmm21  = _mm512_load_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_load_ps(&z[i+112]);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(zmm22,m[j+7],zmm21);
                           _mm512_store_ps(&y[i+112],zmm23);
                            const __m512 zmm24  = _mm512_load_ps(&x[i+128]);
                            const __m512 zmm25  = _mm512_load_ps(&z[i+128]);
                            const __m512 zmm26  = _mm512_mask_cosh_ps(zmm25,m[j+8],zmm24);
                           _mm512_store_ps(&y[i+128],zmm26;
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+144],_MM_HINT_T0); 
                            const __m512 zmm27  = _mm512_load_ps(&x[i+144]);
                            const __m512 zmm28  = _mm512_load_ps(&z[i+144]);
                            const __m512 zmm29  = _mm512_mask_cosh_ps(zmm28,m[j+9],zmm27);
                           _mm512_store_ps(&y[i+144],zmm29;
                           j += 16;
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+144],_MM_HINT_T0); 
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm18  = _mm512_load_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_load_ps(&z[i+96]);
                            const __m512 zmm21  = _mm512_load_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_load_ps(&z[i+112]);
                            const __m512 zmm24  = _mm512_load_ps(&x[i+128]);
                            const __m512 zmm25  = _mm512_load_ps(&z[i+128]);
                            const __m512 zmm27  = _mm512_load_ps(&x[i+144]);
                            const __m512 zmm28  = _mm512_load_ps(&z[i+144]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                            const __m512 zmm14 = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(zmm16,m[j+5],zmm15);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(zmm19,m[j+6],zmm18);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(zmm22,m[j+7],zmm21);
                            const __m512 zmm26  = _mm512_mask_cosh_ps(zmm25,m[j+8],zmm24);
                            const __m512 zmm29  = _mm512_mask_cosh_ps(zmm28,m[j+9],zmm27);
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm512_store_ps(&y[i+48],zmm11);
                           _mm512_store_ps(&y[i+64],zmm14);
                           _mm512_store_ps(&y[i+80],zmm17); 
                           _mm512_store_ps(&y[i+96],zmm20);
                           _mm512_store_ps(&y[i+112],zmm23);
                           _mm512_store_ps(&y[i+128],zmm26);
                           _mm512_store_ps(&y[i+144],zmm29);
                           j += 16;
#endif
                       }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                           _mm512_store_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_store_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                           _mm512_store_ps(&y[i+64],zmm14);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(zmm16,m[j+5],zmm15);
                           _mm512_store_ps(&y[i+80],zmm17);
                            const __m512 zmm18  = _mm512_load_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_load_ps(&z[i+96]);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(zmm19,m[j+6],zmm18);
                           _mm512_store_ps(&y[i+96],zmm20);
                            const __m512 zmm21  = _mm512_load_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_load_ps(&z[i+112]);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(zmm22,m[j+7],zmm21);
                           _mm512_store_ps(&y[i+112],zmm23);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm18  = _mm512_load_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_load_ps(&z[i+96]);
                            const __m512 zmm21  = _mm512_load_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_load_ps(&z[i+112]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                            const __m512 zmm14 = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(zmm16,m[j+5],zmm15);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(zmm19,m[j+6],zmm18);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(zmm22,m[j+7],zmm21);
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm512_store_ps(&y[i+48],zmm11);
                           _mm512_store_ps(&y[i+64],zmm14);
                           _mm512_store_ps(&y[i+80],zmm17); 
                           _mm512_store_ps(&y[i+96],zmm20);
                           _mm512_store_ps(&y[i+112],zmm23);
                           j += 16;
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                           _mm512_store_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_store_ps(&y[i+48],zmm11);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm512_store_ps(&y[i+48],zmm11);
                           j += 16;
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
#else
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           j += 16;
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                           j += 16;

                     } 

                     for(; (i+0) < n; i += 1) {
                           if(m[i])
                              y[i] = cephes_coshf(x[i]+z[i]);
                           else
                              y[i] = cephes_coshf(x[i]);
                     }

                 }
             }


                 
                   void gms::math::coshv_mask_zmm16r4_unroll_10x_u(const float * __restrict  x, 
                                                       const float * __restrict  z,
                                                       float * __restrict  y,
                                                       const __mmask16 *  m,
                                                       const __m512 a,
                                                       const __m512 b,
                                                       const __m512 c,
                                                       const int32_t n,
                                                       const bool additive) {

                     if(__builtin_expect(0==n,0)) {return;}
                     volatile __m512 d;
                     const __
                     int32_t i,j;
                       // Start the preload phase.
                     _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                     _mm_prefetch((const char *)&z[0],_MM_HINT_T0);
                     _mm_prefetch((const char *)&m[0],_MM_HINT_T0);
                     volatile __m512 first = _mm512_mask_cosh_ps(_mm512_load_ps(&z[0],m[0],
                                                                        _mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                     for(int32_t j=0; j != 14; ++j) {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }
                     j = 0;
                     if(additive) {

                          // Main processing loop starts.
                        for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+64],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(_mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+128],_MM_HINT_T0);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                           _mm512_storeu_ps(&y[i+80],zmm17);
                            const __m512 zmm18  = _mm512_loadu_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_loadu_ps(&z[i+96]);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm18,zmm19),m[j+6],zmm18);
                           _mm512_storeu_ps(&y[i+96],zmm20);
                            const __m512 zmm21  = _mm512_loadu_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_loadu_ps(&z[i+112]);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm21,zmm22),m[j+7],zmm21);
                           _mm512_storeu_ps(&y[i+112],zmm23);
                            const __m512 zmm24  = _mm512_loadu_ps(&x[i+128]);
                            const __m512 zmm25  = _mm512_loadu_ps(&z[i+128]);
                            const __m512 zmm26  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm24,zmm25),m[j+8],zmm24);
                           _mm512_storeu_ps(&y[i+128],zmm26;
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+144],_MM_HINT_T0); 
                            const __m512 zmm27  = _mm512_loadu_ps(&x[i+144]);
                            const __m512 zmm28  = _mm512_loadu_ps(&z[i+144]);
                            const __m512 zmm29  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm27,zmm28),m[j+9],zmm27);
                           _mm512_storeu_ps(&y[i+144],zmm29;
                           j += 16;
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+144],_MM_HINT_T0); 
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                            const __m512 zmm18  = _mm512_loadu_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_loadu_ps(&z[i+96]);
                            const __m512 zmm21  = _mm512_loadu_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_loadu_ps(&z[i+112]);
                            const __m512 zmm24  = _mm512_loadu_ps(&x[i+128]);
                            const __m512 zmm25  = _mm512_loadu_ps(&z[i+128]);
                            const __m512 zmm27  = _mm512_loadu_ps(&x[i+144]);
                            const __m512 zmm28  = _mm512_loadu_ps(&z[i+144]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm0,zmm1),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm3,zmm4),m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                            const __m512 zmm14 = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm18,zmm19),m[j+6],zmm18);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm21,zmm22),m[j+7],zmm21);
                            const __m512 zmm26  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm24,zmm25),m[j+8],zmm24);
                            const __m512 zmm29  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm27,zmm28),m[j+9],zmm27);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                           _mm512_storeu_ps(&y[i+80],zmm17); 
                           _mm512_storeu_ps(&y[i+96],zmm20);
                           _mm512_storeu_ps(&y[i+112],zmm23);
                           _mm512_storeu_ps(&y[i+128],zmm26);
                           _mm512_storeu_ps(&y[i+144],zmm29);
                           j += 16;
#endif
                       }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(_mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                           _mm512_storeu_ps(&y[i+80],zmm17);
                            const __m512 zmm18  = _mm512_loadu_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_loadu_ps(&z[i+96]);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm18,zmm19),m[j+6],zmm18);
                           _mm512_storeu_ps(&y[i+96],zmm20);
                            const __m512 zmm21  = _mm512_loadu_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_loadu_ps(&z[i+112]);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm21,zmm22),m[j+7],zmm21);
                           _mm512_storeu_ps(&y[i+112],zmm23);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                            const __m512 zmm18  = _mm512_loadu_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_loadu_ps(&z[i+96]);
                            const __m512 zmm21  = _mm512_loadu_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_loadu_ps(&z[i+112]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(
                                                                _mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(
                                                                _mm512_add_ps(zmm18,zmm19),m[j+6],zmm18);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(
                                                                _mm512_add_ps(zmm21,zmm22),m[j+7],zmm21);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                           _mm512_storeu_ps(&y[i+80],zmm17); 
                           _mm512_storeu_ps(&y[i+96],zmm20);
                           _mm512_storeu_ps(&y[i+112],zmm23);
                           j += 16;
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(_mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           j += 16;
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           j += 16;
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           j += 16;
                     } 

                     for(; (i+0) < n; i += 1) {
                           if(m[i])
                              y[i] = cephes_coshf(x[i]+z[i]);
                           else
                              y[i] = cephes_coshf(x[i]);
                     }


                       }
                        else {
                       // Main processing loop starts.
                          for(i = 0; (i+159) < n; i += 160) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+64],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+128],_MM_HINT_T0);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(zmm16,m[j+5],zmm15);
                           _mm512_storeu_ps(&y[i+80],zmm17);
                            const __m512 zmm18  = _mm512_loadu_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_loadu_ps(&z[i+96]);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(zmm19,m[j+6],zmm18);
                           _mm512_storeu_ps(&y[i+96],zmm20);
                            const __m512 zmm21  = _mm512_loadu_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_loadu_ps(&z[i+112]);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(zmm22,m[j+7],zmm21);
                           _mm512_storeu_ps(&y[i+112],zmm23);
                            const __m512 zmm24  = _mm512_loadu_ps(&x[i+128]);
                            const __m512 zmm25  = _mm512_loadu_ps(&z[i+128]);
                            const __m512 zmm26  = _mm512_mask_cosh_ps(zmm25,m[j+8],zmm24);
                           _mm512_storeu_ps(&y[i+128],zmm26;
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+144],_MM_HINT_T0); 
                            const __m512 zmm27  = _mm512_loadu_ps(&x[i+144]);
                            const __m512 zmm28  = _mm512_loadu_ps(&z[i+144]);
                            const __m512 zmm29  = _mm512_mask_cosh_ps(zmm28,m[j+9],zmm27);
                           _mm512_storeu_ps(&y[i+144],zmm29;
                           j += 16;
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+144],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+144],_MM_HINT_T0); 
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                            const __m512 zmm18  = _mm512_loadu_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_loadu_ps(&z[i+96]);
                            const __m512 zmm21  = _mm512_loadu_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_loadu_ps(&z[i+112]);
                            const __m512 zmm24  = _mm512_loadu_ps(&x[i+128]);
                            const __m512 zmm25  = _mm512_loadu_ps(&z[i+128]);
                            const __m512 zmm27  = _mm512_loadu_ps(&x[i+144]);
                            const __m512 zmm28  = _mm512_loadu_ps(&z[i+144]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                            const __m512 zmm14 = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(zmm16,m[j+5],zmm15);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(zmm19,m[j+6],zmm18);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(zmm22,m[j+7],zmm21);
                            const __m512 zmm26  = _mm512_mask_cosh_ps(zmm25,m[j+8],zmm24);
                            const __m512 zmm29  = _mm512_mask_cosh_ps(zmm28,m[j+9],zmm27);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                           _mm512_storeu_ps(&y[i+80],zmm17); 
                           _mm512_storeu_ps(&y[i+96],zmm20);
                           _mm512_storeu_ps(&y[i+112],zmm23);
                           _mm512_storeu_ps(&y[i+128],zmm26);
                           _mm512_storeu_ps(&y[i+144],zmm29);
                           j += 16;
#endif
                       }

                       for(; (i+127) < n; i += 128) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_masku_sin_ps(zmm16,m[j+5],zmm15);
                           _mm512_storeu_ps(&y[i+80],zmm17);
                            const __m512 zmm18  = _mm512_loadu_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_loadu_ps(&z[i+96]);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(zmm19,m[j+6],zmm18);
                           _mm512_storeu_ps(&y[i+96],zmm20);
                            const __m512 zmm21  = _mm512_loadu_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_loadu_ps(&z[i+112]);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(zmm22,m[j+7],zmm21);
                           _mm512_storeu_ps(&y[i+112],zmm23);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                            const __m512 zmm18  = _mm512_loadu_ps(&x[i+96]);
                            const __m512 zmm19  = _mm512_loadu_ps(&z[i+96]);
                            const __m512 zmm21  = _mm512_loadu_ps(&x[i+112]);
                            const __m512 zmm22  = _mm512_loadu_ps(&z[i+112]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                            const __m512 zmm14 = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(zmm16,m[j+5],zmm15);
                            const __m512 zmm20  = _mm512_mask_cosh_ps(zmm19,m[j+6],zmm18);
                            const __m512 zmm23  = _mm512_mask_cosh_ps(zmm22,m[j+7],zmm21);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                           _mm512_storeu_ps(&y[i+80],zmm17); 
                           _mm512_storeu_ps(&y[i+96],zmm20);
                           _mm512_storeu_ps(&y[i+112],zmm23);
                           j += 16;
#endif
 
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           j += 16;
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           j += 16;
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           j += 16;
                     } 

                     for(; (i+0) < n; i += 1) {
                           if(m[i])
                              y[i] = cephes_coshf(x[i]+z[i]);
                           else
                              y[i] = cephes_coshf(x[i]);
                     }

                 }
             }




              /*
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cosh_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/
                 
                   void gms::math::coshv_zmm16r4_unroll_10x_u(const float * __restrict  x,
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
                       volatile __m512 first = _mm512_cosh_ps(_mm512_loadu_ps(&x[0])); //L1I cache miss.
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
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const __m512 zmm10      = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cosh_ps(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                            const __m512 zmm12      = _mm512_loadu_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_cosh_ps(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                            const __m512 zmm14      = _mm512_loadu_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_cosh_ps(zmm14);
                           _mm512_storeu_ps(&y[i+112], zmm15);
                            const __m512 zmm16      = _mm512_loadu_ps(&x[i+128]); 
                            const __m512 zmm17      = _mm512_cosh_ps(zmm16);
                           _mm512_storeu_ps(&y[i+128], zmm17);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                            const __m512 zmm18      = _mm512_loadu_ps(&x[i+144]); 
                            const __m512 zmm19      = _mm512_cosh_ps(zmm18);
                           _mm512_storeu_ps(&y[i+144], zmm19);
#else
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm10     = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm12     = _mm512_loadu_ps(&x[i+96]); 
                            const __m512 zmm14     = _mm512_loadu_ps(&x[i+112]);  
                            const __m512 zmm16     = _mm512_loadu_ps(&x[i+128]); 
                            const __m512 zmm18     = _mm512_loadu_ps(&x[i+144]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                            const __m512 zmm11     = _mm512_cosh_ps(zmm10);
                            const __m512 zmm13     = _mm512_cosh_ps(zmm12);
                            const __m512 zmm15     = _mm512_cosh_ps(zmm14);
                            const __m512 zmm17     = _mm512_cosh_ps(zmm16);
                            const __m512 zmm19     = _mm512_cosh_ps(zmm18);
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
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const __m512 zmm10      = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cosh_ps(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                            const __m512 zmm12      = _mm512_loadu_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_cosh_ps(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                            const __m512 zmm14      = _mm512_loadu_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_cosh_ps(zmm14);
                           _mm512_storeu_ps(&y[i+112], zmm15);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm10     = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm12     = _mm512_loadu_ps(&x[i+96]); 
                            const __m512 zmm14     = _mm512_loadu_ps(&x[i+112]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                            const __m512 zmm11     = _mm512_cosh_ps(zmm10);
                            const __m512 zmm13     = _mm512_cosh_ps(zmm12);
                            const __m512 zmm15     = _mm512_cosh_ps(zmm14);
                            const __m512 zmm17     = _mm512_cosh_ps(zmm16);
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
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#endif
                      }

                     for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#endif
                     }

                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0); 
                           _mm512_storeu_ps(&y[i+0], zmm1);
#endif
                     } 

                     for(; (i+0) < n; i += 1) {
                           y[i] = cephes_coshf(x[i]);
                     }
               }

               
              

                 

 /*             
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cosh_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/

                 
                   void gms::math::coshv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) x,
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
                       volatile __m512 first = _mm512_cosh_ps(_mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

                       // Main processing loop starts here
                       for(i = 0; (i+95) < n; i += 96) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                            const __m512 zmm10      = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cosh_ps(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
#else
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm10     = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                            const __m512 zmm11     = _mm512_cosh_ps(zmm10);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
#endif
                     }

                     for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
#endif
                     }

                    for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3); 
#endif

                    }

                    for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
#endif
                    }

                    for(; (i+0) < n; i += 1) {
                           y[i] = cephes_coshf(x[i]);
                     }
               }


               /*             
      These sine vector kernels call Intel SVML library
      sine function implementation which is not inlined
      by the ICC/ICPC compilers, hence a pre-load call
      to _mm512_cosh_ps and warmup loop are inserted in
      order to mitigate as far as it is possible the 
      issue of impossibility of caching ahead of time. 
*/

                  
                   void gms::math::coshv_zmm16r4_unroll_6x_u(const float * __restrict  x,
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
                       volatile __m512 first = _mm512_cosh_ps(_mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

                       // Main processing loop starts here
                       for(i = 0; (i+95) < n; i += 96) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                            const __m512 zmm10      = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_cosh_ps(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
#else
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm10     = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                            const __m512 zmm9      = _mm512_cosh_ps(zmm8);
                            const __m512 zmm11     = _mm512_cosh_ps(zmm10);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm512_storeu_ps(&y[i+80], zmm11);
#endif
                     }

                     for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                            const __m512 zmm5      = _mm512_cosh_ps(zmm4);
                            const __m512 zmm7      = _mm512_cosh_ps(zmm6);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#endif
                     }

                    for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                            const __m512 zmm3      = _mm512_cosh_ps(zmm2);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3); 
#endif

                    }

                    for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_cosh_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
#endif
                    }

                    for(; (i+0) < n; i += 1) {
                           y[i] = cephes_coshf(x[i]);
                     }
               }


                

                   
                   void gms::math::coshv_mask_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) x, 
                                                       const float * __restrict __ATTR_ALIGN__(64) z,
                                                       float * __restrict __ATTR_ALIGN__(64) y,
                                                       const __mmask16 * __restrict __ATTR_ALIGN__(64) m,
                                                       const __m512 a,
                                                       const __m512 b,
                                                       const __m512 c,
                                                       const int32_t n,
                                                       const bool additive) {

                     if(__builtin_expect(0==n,0)) {return;}
                     volatile __m512 d;
                     const __
                     int32_t i,j;
                       // Start the preload phase.
                     _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                     _mm_prefetch((const char *)&z[0],_MM_HINT_T0);
                     _mm_prefetch((const char *)&m[0],_MM_HINT_T0);
                     volatile __m512 first = _mm512_mask_cosh_ps(_mm512_load_ps(&z[0],m[0],
                                                                        _mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                     for(int32_t j=0; j != 14; ++j) {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }
                     j = 0;
                     if(additive) {

                          // Main processing loop starts.
                        for(i = 0; (i+95) < n; i += 96) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+48],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+96],_MM_HINT_T0);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(_mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_store_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                           _mm512_store_ps(&y[i+64],zmm14);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                           _mm512_store_ps(&y[i+80],zmm17);
                           j += 16;
#else
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+96],_MM_HINT_T0);
                          
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                          
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm0,zmm1),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm3,zmm4),m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                            const __m512 zmm14 = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                          
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm512_store_ps(&y[i+48],zmm11);
                           _mm512_store_ps(&y[i+64],zmm14);
                           _mm512_store_ps(&y[i+80],zmm17); 
                           j += 16;
                          
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                           _mm512_store_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(_mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_store_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                           _mm512_store_ps(&y[i+64],zmm14);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                           _mm512_store_ps(&y[i+80],zmm17);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                          
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(
                                                                _mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                          
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm512_store_ps(&y[i+48],zmm11);
                           _mm512_store_ps(&y[i+64],zmm14);
                           _mm512_store_ps(&y[i+80],zmm17); 
                           j += 16;
#endif
                       }

                       for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                          
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                          
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           j += 16;
#endif
                      }

                     for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                           j += 16;

                     }

                    
                     for(; (i+0) < n; i += 1) {
                           if(m[i])
                              y[i] = cephes_coshf(x[i]+z[i]);
                           else
                              y[i] = cephes_coshf(x[i]);
                     }


                       }
                        else {
                       // Main processing loop starts.
                          for(i = 0; (i+95) < n; i += 96) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+48],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+96],_MM_HINT_T0);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_store_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                           _mm512_store_ps(&y[i+64],zmm14);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(zmm16,m[j+5],zmm15);
                           _mm512_store_ps(&y[i+80],zmm17);
                           j += 16;
#else
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+96],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_load_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_load_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_load_ps(&z[i+80]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                            const __m512 zmm14 = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(zmm16,m[j+5],zmm15);
                         
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm512_store_ps(&y[i+48],zmm11);
                           _mm512_store_ps(&y[i+64],zmm14);
                           _mm512_store_ps(&y[i+80],zmm17); 
                           j += 16;
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                           _mm512_store_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_store_ps(&y[i+48],zmm11);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_load_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_load_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_load_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_load_ps(&z[i+48]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm512_store_ps(&y[i+48],zmm11);
                           j += 16;
#endif
                       }

                       for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_store_ps(&y[i+16],zmm5);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_load_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_load_ps(&z[i+16]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_store_ps(&y[i+0],zmm2);
                           _mm512_store_ps(&y[i+16],zmm5);
                           j += 16;
#endif
                      }

                  
                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_load_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                           j += 16;
                     } 

                     for(; (i+0) < n; i += 1) {
                           if(m[i])
                              y[i] = cephes_coshf(x[i]+z[i]);
                           else
                              y[i] = cephes_coshf(x[i]);
                     }

                 }
             }


                  
                   void gms::math::coshv_mask_zmm16r4_unroll_6x_u(const float * __restrict  x, 
                                                       const float * __restrict  z,
                                                       float * __restrict  y,
                                                       const __mmask16 *  m,
                                                       const __m512 a,
                                                       const __m512 b,
                                                       const __m512 c,
                                                       const int32_t n,
                                                       const bool additive) {

                     if(__builtin_expect(0==n,0)) {return;}
                     volatile __m512 d;
                     const __
                     int32_t i,j;
                       // Start the preload phase.
                     _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                     _mm_prefetch((const char *)&z[0],_MM_HINT_T0);
                     _mm_prefetch((const char *)&m[0],_MM_HINT_T0);
                     volatile __m512 first = _mm512_mask_cosh_ps(_mm512_loadu_ps(&z[0],m[0],
                                                                        _mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                     for(int32_t j=0; j != 14; ++j) {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

                     if(additive) {

                          // Main processing loop starts.
                        for(i = 0; (i+95) < n; i += 96) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+48],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+96],_MM_HINT_T0);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(_mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                           _mm512_storeu_ps(&y[i+80],zmm17);
                           j += 16;
#else
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+96],_MM_HINT_T0);
                          
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                          
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm0,zmm1),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm3,zmm4),m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                            const __m512 zmm14 = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(
                                                                   _mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                          
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                           _mm512_storeu_ps(&y[i+80],zmm17); 
                           j += 16;
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(_mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                           _mm512_storeu_ps(&y[i+80],zmm17);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                          
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm6,zmm7),m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm9,zmm10),m[j+3],zmm9);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm12,zmm13),m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(
                                                                _mm512_add_ps(zmm15,zmm16),m[j+5],zmm15);
                          
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                           _mm512_storeu_ps(&y[i+80],zmm17); 
                           j += 16;
#endif
                       }

                       for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           j += 16;
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                          
                            const __m512 zmm2  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(
                                                               _mm512_add_ps(zmm4,zmm3),m[j+1],zmm3);
                          
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           j += 16;
#endif
                      }

                     for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(_mm512_add_ps(zmm1,zmm0),m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           j += 16;

                     }

                    
                     for(; (i+0) < n; i += 1) {
                           if(m[i])
                              y[i] = cephes_coshf(x[i]+z[i]);
                           else
                              y[i] = cephes_coshf(x[i]);
                     }


                       }
                        else {
                       // Main processing loop starts.
                          for(i = 0; (i+95) < n; i += 96) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+48],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                           _mm512_store_ps(&y[i+32],zmm8);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+96],_MM_HINT_T0);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm14  = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(zmm16,m[j+5],zmm15);
                           _mm512_storeu_ps(&y[i+80],zmm17);
                           j += 16;
#else
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&y[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&m[i+96],_MM_HINT_T0);
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm12  = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm13  = _mm512_loadu_ps(&z[i+64]);
                            const __m512 zmm15  = _mm512_loadu_ps(&x[i+80]);
                            const __m512 zmm16  = _mm512_loadu_ps(&z[i+80]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                            const __m512 zmm14 = _mm512_mask_cosh_ps(zmm13,m[j+4],zmm12);
                            const __m512 zmm17  = _mm512_mask_cosh_ps(zmm16,m[j+5],zmm15);
                         
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           _mm512_storeu_ps(&y[i+64],zmm14);
                           _mm512_storeu_ps(&y[i+80],zmm17); 
                           j += 16;
#endif
                       }

                       for(; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           j += 16;
                          
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm6  = _mm512_loadu_ps(&x[i+32]);
                            const __m512 zmm7  = _mm512_loadu_ps(&z[i+32]);
                            const __m512 zmm9  = _mm512_loadu_ps(&x[i+48]);
                            const __m512 zmm10 = _mm512_loadu_ps(&z[i+48]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                            const __m512 zmm8  = _mm512_mask_cosh_ps(zmm7,m[j+2],zmm6);
                            const __m512 zmm11 = _mm512_mask_cosh_ps(zmm10,m[j+3],zmm9);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           _mm512_storeu_ps(&y[i+32],zmm8);
                           _mm512_storeu_ps(&y[i+48],zmm11);
                           j += 16;
#endif
                       }

                       for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_store_ps(&y[i+0],zmm2);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           j += 16;
                         
#else
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm3  = _mm512_loadu_ps(&x[i+16]);
                            const __m512 zmm4  = _mm512_loadu_ps(&z[i+16]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                            const __m512 zmm5  = _mm512_mask_cosh_ps(zmm4,m[j+1],zmm3);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           _mm512_storeu_ps(&y[i+16],zmm5);
                           j += 16;
#endif
                      }

                  
                     for(; (i+15) < n; i += 15) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                            const __m512 zmm0  = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1  = _mm512_loadu_ps(&z[i+0]);
                            const __m512 zmm2  = _mm512_mask_cosh_ps(zmm1,m[j+0],zmm0);
                           _mm512_storeu_ps(&y[i+0],zmm2);
                           j += 16;
                     } 

                     for(; (i+0) < n; i += 1) {
                           if(m[i])
                              y[i] = cephes_coshf(x[i]+z[i]);
                           else
                              y[i] = cephes_coshf(x[i]);
                     }

                 }
             }







              

                    

   
