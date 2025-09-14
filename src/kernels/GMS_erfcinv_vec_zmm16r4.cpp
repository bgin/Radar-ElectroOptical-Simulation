



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


#include <cmath>
#include "GMS_erfcinv_vec_zmm16r4.h"




                  
void gms::math::erfcinvv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                              float * __restrict __ATTR_ALIGN__(64) y,
                                              const __m512 a,
                                              const __m512 b,
                                              const __m512 c,
                                              const int32_t n) 
{

                       if(__builtin_expect(0==n,0)) {return;}
                       volatile __m512 d;
                       int32_t i;
                       // Start the preload phase.
                       _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                       volatile __m512 first = _mm512_erfcinv_ps(_mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) 
                       {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

                       // Main processing loop starts.
                       for(i = 0; (i+159) < n; i += 160) 
                       {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
#endif
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
#endif
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
#endif
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
#endif
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
#endif                          
                            const __m512 zmm10      = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_erfcinv_ps(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
#endif
                            const __m512 zmm12      = _mm512_load_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_erfcinv_ps(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+112],_MM_HINT_T0);
#endif
                            const __m512 zmm14      = _mm512_load_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_erfcinv_ps(zmm14);
                           _mm512_store_ps(&y[i+112], zmm15);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
#endif
                            const __m512 zmm16      = _mm512_load_ps(&x[i+128]); 
                            const __m512 zmm17      = _mm512_erfcinv_ps(zmm16);
                           _mm512_store_ps(&y[i+128], zmm17);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
#endif                          
                            const __m512 zmm18      = _mm512_load_ps(&x[i+144]); 
                            const __m512 zmm19      = _mm512_erfcinv_ps(zmm18);
                           _mm512_store_ps(&y[i+144], zmm19);
#else
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+112],_MM_HINT_T0);
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
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                            const __m512 zmm11     = _mm512_erfcinv_ps(zmm10);
                            const __m512 zmm13     = _mm512_erfcinv_ps(zmm12);
                            const __m512 zmm15     = _mm512_erfcinv_ps(zmm14);
                            const __m512 zmm17     = _mm512_erfcinv_ps(zmm16);
                            const __m512 zmm19     = _mm512_erfcinv_ps(zmm18);
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

                       for(; (i+127) < n; i += 128) 
                       {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const __m512 zmm10      = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_erfcinv_ps(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
                            const __m512 zmm12      = _mm512_load_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_erfcinv_ps(zmm12);
                           _mm512_store_ps(&y[i+96], zmm13);
                            const __m512 zmm14      = _mm512_load_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_erfcinv_ps(zmm14);
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
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                            const __m512 zmm11     = _mm512_erfcinv_ps(zmm10);
                            const __m512 zmm13     = _mm512_erfcinv_ps(zmm12);
                            const __m512 zmm15     = _mm512_erfcinv_ps(zmm14);
                            const __m512 zmm17     = _mm512_erfcinv_ps(zmm16);
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

                       for(; (i+63) < n; i += 64) 
                       {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
#endif
                      }

                     for(; (i+31) < n; i += 32) 
                     {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
#endif
                     }

                     for(; (i+15) < n; i += 15) 
                     {

                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]);
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0); 
                           _mm512_store_ps(&y[i+0], zmm1);

                     } 

                    for(; (i+0) < n; i += 1) 
                    {
                          y[i] = 1.0f/std::erf(1.0f-x[i]); //std::erfc(1.0f-x[i]);
                    }
}

                  
                
void gms::math::erfcinvv_zmm16r4_unroll_10x_u(const float * __restrict  x,
                                              float * __restrict  y,
                                              const __m512 a,
                                              const __m512 b,
                                              const __m512 c,
                                              int32_t n) 
{

                       if(__builtin_expect(0==n,0)) {return;}
                       volatile __m512 d;
                       int32_t i;
                       // Start the preload phase.
                       _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                       volatile __m512 first = _mm512_erfcinv_ps(_mm512_loadu_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) 
                       {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }
                       
#if (ERFCINV_VEC_ZMM16R4_PEEL_LOOP) == 1
                        
                        while(((uintptr_t)&y & 63) && n) 
                        {
                              const float xt  = *x; 
                              const float arg = 1.0f-xt;
                              *y              = 1.0f/std::erf(arg);
                              x++;
                              y++;
                              n--; 
                        }
#endif 
                       // Main processing loop starts.
                       for(i = 0; (i+159) < n; i += 160) 
                       {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
#endif
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
#endif
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
#endif
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
#endif
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
#endif                         
                            const __m512 zmm10      = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_erfcinv_ps(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
#endif
                            const __m512 zmm12      = _mm512_loadu_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_erfcinv_ps(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+112],_MM_HINT_T0);
#endif
                            const __m512 zmm14      = _mm512_loadu_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_erfcinv_ps(zmm14);
                           _mm512_storeu_ps(&y[i+112], zmm15);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
#endif
                            const __m512 zmm16      = _mm512_loadu_ps(&x[i+128]); 
                            const __m512 zmm17      = _mm512_erfcinv_ps(zmm16);
                           _mm512_storeu_ps(&y[i+128], zmm17);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+144],_MM_HINT_T0);
#endif                          
                            const __m512 zmm18      = _mm512_loadu_ps(&x[i+144]); 
                            const __m512 zmm19      = _mm512_erfcinv_ps(zmm18);
                           _mm512_storeu_ps(&y[i+144], zmm19);
#else
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+96],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+112],_MM_HINT_T0);
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
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                            const __m512 zmm11     = _mm512_erfcinv_ps(zmm10);
                            const __m512 zmm13     = _mm512_erfcinv_ps(zmm12);
                            const __m512 zmm15     = _mm512_erfcinv_ps(zmm14);
                            const __m512 zmm17     = _mm512_erfcinv_ps(zmm16);
                            const __m512 zmm19     = _mm512_erfcinv_ps(zmm18);
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

                       for(; (i+127) < n; i += 128) 
                       {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm_prefetch((const char *)&x[i+128],_MM_HINT_T0);
                            const __m512 zmm10      = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_erfcinv_ps(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
                            const __m512 zmm12      = _mm512_loadu_ps(&x[i+96]); 
                            const __m512 zmm13      = _mm512_erfcinv_ps(zmm12);
                           _mm512_storeu_ps(&y[i+96], zmm13);
                            const __m512 zmm14      = _mm512_loadu_ps(&x[i+112]); 
                            const __m512 zmm15      = _mm512_erfcinv_ps(zmm14);
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
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                            const __m512 zmm11     = _mm512_erfcinv_ps(zmm10);
                            const __m512 zmm13     = _mm512_erfcinv_ps(zmm12);
                            const __m512 zmm15     = _mm512_erfcinv_ps(zmm14);
                            const __m512 zmm17     = _mm512_erfcinv_ps(zmm16);
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

                       for(; (i+63) < n; i += 64) 
                       {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#endif
                      }

                     for(; (i+31) < n; i += 32) 
                     {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#endif
                     }

                     for(; (i+15) < n; i += 15) 
                     {

                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]);
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0); 
                           _mm512_storeu_ps(&y[i+0], zmm1);

                     } 

                     for(; (i+0) < n; i += 1) 
                     {
                          y[i] = 1.0f/std::erf(1.0f-x[i]); //std::erfc(1.0f-x[i]);
                     }
}

               
              

                 

 

                  
void gms::math::erfcinvv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) x,
                                             float * __restrict __ATTR_ALIGN__(64) y,
                                             const __m512 a,
                                             const __m512 b,
                                             const __m512 c,
                                             const int32_t n) 
{

                       if(__builtin_expect(0==n,0)) {return;}
                       volatile __m512 d;
                       int32_t i;
                       // Start the preload phase.
                       _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                       volatile __m512 first = _mm512_erfcinv_ps(_mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) 
                       {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

                       // Main processing loop starts here
                       for(i = 0; (i+95) < n; i += 96) 
                       {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
#endif 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
#endif 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
#endif 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
#endif                          
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                           _mm512_store_ps(&y[i+64], zmm9);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
#endif 
                            const __m512 zmm10      = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_erfcinv_ps(zmm10);
                           _mm512_store_ps(&y[i+80], zmm11);
#else
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_load_ps(&x[i+64]);
                            const __m512 zmm10     = _mm512_load_ps(&x[i+80]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                            const __m512 zmm11     = _mm512_erfcinv_ps(zmm10);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
                           _mm512_store_ps(&y[i+64], zmm9);
                           _mm512_store_ps(&y[i+80], zmm11);
#endif
                     }

                     for(; (i+63) < n; i += 64) 
                     {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                           _mm512_store_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_store_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_load_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_load_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3);
                           _mm512_store_ps(&y[i+32], zmm5);
                           _mm512_store_ps(&y[i+48], zmm7);
#endif
                     }

                    for(; (i+31) < n; i += 32) 
                    {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_store_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_load_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_store_ps(&y[i+0], zmm1);
                           _mm512_store_ps(&y[i+16], zmm3); 
#endif

                    }

                    for(; (i+15) < n; i += 16) 
                    {

                            const __m512 zmm0      = _mm512_load_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_store_ps(&y[i+0], zmm1);

                    }

                    for(; (i+0) < n; i += 1) 
                    {
                          y[i] = 1.0f/std::erf(1.0f-x[i]); //std::erfc(1.0f-x[i]);
                    }
}


  

void gms::math::erfcinvv_zmm16r4_unroll_6x_u(const float * __restrict  x,
                                             float * __restrict  y,
                                             const __m512 a,
                                             const __m512 b,
                                             const __m512 c,
                                             int32_t n) 
{

                       if(__builtin_expect(0==n,0)) {return;}
                       volatile __m512 d;
                       int32_t i;
                       // Start the preload phase.
                       _mm_prefetch((const char *)&x[0],_MM_HINT_T0);
                       volatile __m512 first = _mm512_erfcinv_ps(_mm512_load_ps(&x[0])); //L1I cache miss.
                       // Warmup loop of ~100 cycles (worst case scenario) of memory-fetch machine 
                       // code instructions, needed to keep core busy while waiting on instructions
                       // arrival, this is done to prevent the logic from progressing towards main 
                       // loop.
                       for(int32_t j=0; j != 14; ++j) 
                       {
                           d = _mm512_add_ps(d,_mm512_fmadd_ps(a,b,c));
                       }

#if (ERFCINV_VEC_ZMM16R4_PEEL_LOOP) == 1
                       
                        while(((uintptr_t)&y & 63) && n) 
                        {
                              const float xt  = *x; 
                              const float arg = 1.0f-xt;
                              *y              = 1.0f/std::erf(arg);
                              x++;
                              y++;
                              n--; 
                        }
#endif 
                       // Main processing loop starts here
                       for(i = 0; (i+95) < n; i += 96) {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
#endif 
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
#endif 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
#endif 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
#endif 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
#endif                           
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]); 
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                           _mm512_storeu_ps(&y[i+64], zmm9);
#if (ERFCINV_VEC_ZMM16R4_SOFT_PREFETCH) == 1
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
#endif 
                            const __m512 zmm10      = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm11      = _mm512_erfcinv_ps(zmm10);
                           _mm512_storeu_ps(&y[i+80], zmm11);
#else
                           _mm_prefetch((const char *)&x[i+0],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+16],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+32],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+48],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+64],_MM_HINT_T0);
                           _mm_prefetch((const char *)&x[i+80],_MM_HINT_T0);
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm8      = _mm512_loadu_ps(&x[i+64]);
                            const __m512 zmm10     = _mm512_loadu_ps(&x[i+80]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                            const __m512 zmm9      = _mm512_erfcinv_ps(zmm8);
                            const __m512 zmm11     = _mm512_erfcinv_ps(zmm10);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
                           _mm512_storeu_ps(&y[i+64], zmm9);
                           _mm512_storeu_ps(&y[i+80], zmm11);
#endif
                     }

                     for(; (i+63) < n; i += 64) 
                     {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm4      = _mm512_loadu_ps(&x[i+32]); 
                            const __m512 zmm6      = _mm512_loadu_ps(&x[i+48]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                            const __m512 zmm5      = _mm512_erfcinv_ps(zmm4);
                            const __m512 zmm7      = _mm512_erfcinv_ps(zmm6);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3);
                           _mm512_storeu_ps(&y[i+32], zmm5);
                           _mm512_storeu_ps(&y[i+48], zmm7);
#endif
                     }

                    for(; (i+31) < n; i += 32) 
                    {
#if (ERFCINV_VEC_ZMM16R4_INTERLEAVE_SIMD_OPS) == 1
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_storeu_ps(&y[i+16], zmm3);
#else
                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm2      = _mm512_loadu_ps(&x[i+16]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                            const __m512 zmm3      = _mm512_erfcinv_ps(zmm2);
                           _mm512_storeu_ps(&y[i+0], zmm1);
                           _mm512_storeu_ps(&y[i+16], zmm3); 
#endif

                    }

                    for(; (i+15) < n; i += 16) 
                    {

                            const __m512 zmm0      = _mm512_loadu_ps(&x[i+0]); 
                            const __m512 zmm1      = _mm512_erfcinv_ps(zmm0);
                           _mm512_storeu_ps(&y[i+0], zmm1);

                    }

                     for(; (i+0) < n; i += 1) 
                     {
                          y[i] = 1.0f/std::erf(1.0f-x[i]); //std::erfc(1.0f-x[i]);
                     }
}


                

  


              

                    

