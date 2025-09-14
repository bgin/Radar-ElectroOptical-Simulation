



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
#include "GMS_cmean_cprod_zmm16r4.h"




                   void gms::math::cmean_prod_u10x_zmm16r4_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict yre,
                                                  const float * __restrict yim,
                                                  float * __restrict mre,
                                                  float * __restrict mim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) {return;}
                          __m512 zmm0,zmm1,zmm2,zmm3;
                          __m512 zmm4,zmm5,zmm6,zmm7;
                          __m512 zmm8,zmm9,zmm10,zmm11;
                          __m512 zmm12,zmm13,zmm14,zmm15;
                          __m512 zmm16,zmm17,zmm18,zmm19;
                          __m512 zmm20,zmm21,zmm22,zmm23;
                          __m512 zmm24,zmm25,zmm26,zmm27;
                          __m512 zmm28,zmm29,zmm30,zmm31;
                         __m512 redr[10] = {_mm512_setzero_ps()};
                         __m512 redi[10] = {_mm512_setzero_ps()};
                          float re,im;
                         constexpr float invN = 1.0f/static_cast<float>(n);
                         int32_t i; 
                         re = 0.0f;
                         im = 0.0f;
                         for(i = 0; (i+159) < n; i += 160) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             redr[0]  = _mm512_add_ps(redr[0],_mm512_fmsub_ps(zmm0,zmm1,
                                                              _mm512_mul_ps(zmm2,zmm3))); // rep
                             redi[0]  = _mm512_add_ps(redi[0],_mm512_fmadd_ps(zmm2,zmm1,
                                                              _mm512_mul_ps(zmm0,zmm3))); // imp
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             redr[1] = _mm512_add_ps(redr[1],_mm512_fmsub_ps(zmm6,zmm7,
                                                             _mm512_mul_ps(zmm8,zmm9))); // rep
                             redi[1] = _mm512_add_ps(redi[1],_mm512_fmadd_ps(zmm8,zmm7,
                                                     _mm512_mul_ps(zmm6,zmm9))); // imp
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                             zmm12  = _mm512_loadu_ps(&xre[i+32]);
                             zmm13  = _mm512_loadu_ps(&yre[i+32]);
                             zmm14  = _mm512_loadu_ps(&xim[i+32]);
                             zmm15  = _mm512_loadu_ps(&yim[i+32]);
                             redr[2] = _mm512_add_ps(redr[2],_mm512_fmsub_ps(zmm12,zmm13,
                                                     _mm512_mul_ps(zmm14,zmm15))); // rep
                             redi[2] = _mm512_add_ps(redi[2],_mm512_fmadd_ps(zmm14,zmm13,
                                                     _mm512_mul_ps(zmm12,zmm15))); // imp
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             redr[3] = _mm512_add_ps(redr[3],_mm512_fmsub_ps(zmm18,zmm19,
                                                     _mm512_mul_ps(zmm20,zmm21))); // rep
                             redi[3] = _mm512_add_ps(redi[3],_mm512_fmadd_ps(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21))); // imp
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0);
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             redr[4] = _mm512_add_ps(redr[3],_mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27))); // rep
                             redi[4] = _mm512_add_ps(redi[3],_mm512_fmadd_ps(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27))); // imp
                             zmm30 = _mm512_loadu_ps(&xre[i+80]);
                             zmm31 = _mm512_loadu_ps(&yre[i+80]);
                             zmm0  = _mm512_loadu_ps(&xim[i+80]);
                             zmm1  = _mm512_loadu_ps(&yim[i+80]);
                             redr[5]  = _mm512_add_ps(redr[5],_mm512_fmsub_ps(zmm30,zmm31,
                                               _mm512_mul_ps(zmm0,zmm1))); // rep
                             redi[5]  = _mm512_add_ps(redi[5],_mm512_fmadd_ps(zmm0,zmm31,
                                               _mm512_mul_ps(zmm30,zmm1))); // imp
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                             zmm4  = _mm512_loadu_ps(&xre[i+96]);
                             zmm5  = _mm512_loadu_ps(&yre[i+96]);
                             zmm6  = _mm512_loadu_ps(&xim[i+96]);
                             zmm7  = _mm512_loadu_ps(&yim[i+96]);
                             redr[6]  = _mm512_add_ps(redr[6],_mm512_fmsub_ps(zmm4,zmm5,
                                               _mm512_mul_ps(zmm6,zmm7))); // rep
                             redi[6]  = _mm512_add_ps(redi[6],_mm512_fmadd_ps(zmm6,zmm5,
                                               _mm512_mul_ps(zmm4,zmm7))); // imp
                             zmm10 = _mm512_loadu_ps(&xre[i+112]);
                             zmm11 = _mm512_loadu_ps(&yre[i+112]);
                             zmm12 = _mm512_loadu_ps(&xim[i+112]);
                             zmm13 = _mm512_loadu_ps(&yim[i+112]);
                             redr[7] = _mm512_add_ps(redr[7],_mm512_fmsub_ps(zmm10,zmm11,
                                               _mm512_mul_ps(zmm12,zmm13))); // rep
                             redi[7] = _mm512_add_ps(redi[7],_mm512_fmadd_ps(zmm13,zmm12,
                                               _mm512_mul_ps(zmm10,zmm14))); // imp
                             _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+160],_MM_HINT_T0);
                             zmm16 = _mm512_loadu_ps(&xre[i+128]);
                             zmm17 = _mm512_loadu_ps(&yre[i+128]);
                             zmm18 = _mm512_loadu_ps(&xim[i+128]);
                             zmm19 = _mm512_loadu_ps(&yim[i+128]);
                             redr[8] = _mm512_add_ps(redr[8],_mm512_fmsub_ps(zmm16,zmm17,
                                               _mm512_mul_ps(zmm18,zmm19))); // rep
                             redi[8] = _mm512_add_ps(redi[8],_mm512_fmadd_ps(zmm18,zmm17,
                                               _mm512_mul_ps(zmm16,zmm19))); // imp
                             zmm22 = _mm512_loadu_ps(&xre[i+144]);
                             zmm23 = _mm512_loadu_ps(&yre[i+144]);
                             zmm24 = _mm512_loadu_ps(&xim[i+144]);
                             zmm25 = _mm512_loadu_ps(&yim[i+144]);
                             redr[9] = _mm512_add_ps(redr[9],_mm512_fmsub_ps(zmm22,zmm23,
                                               _mm512_mul_ps(zmm24,zmm25))); // rep
                             redi[9] = _mm512_add_ps(redi[9],_mm512_fmadd_ps(zmm24,zmm23,
                                               _mm512_mul_ps(zmm22,zmm25))); // imp
                           
                        }

                             redr[0] = _mm512_add_ps(redr[0],redr[5]);
                             redi[0] = _mm512_add_ps(redi[0],redi[5]);
                             redr[1] = _mm512_add_ps(redr[1],redr[6]);
                             redi[1] = _mm512_add_ps(redi[1],redi[6]);
                             redr[2] = _mm512_add_ps(redr[2],redr[7]);
                             redi[2] = _mm512_add_ps(redi[2],redi[7]);
                             redr[3] = _mm512_add_ps(redr[3],redr[8]);
                             redi[3] = _mm512_add_ps(redi[3],redi[8]);
                             redr[4] = _mm512_add_ps(redr[4],redr[9]);
                             redi[4] = _mm512_add_ps(redi[4],redi[9]); 

                         for(; (i+79) < n; i += 80) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             redr[0]  = _mm512_add_ps(redr[0],_mm512_fmsub_ps(zmm0,zmm1,
                                                               _mm512_mul_ps(zmm2,zmm3))); // rep
                             redi[0]  = _mm512_add_ps(redi[0],_mm512_fmadd_ps(zmm2,zmm1,
                                                              _mm512_mul_ps(zmm0,zmm3))); // imp
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             redr[1] = _mm512_add_ps(redr[1],_mm512_fmsub_ps(zmm6,zmm7,
                                                             _mm512_mul_ps(zmm8,zmm9))); // rep
                             redi[1] = _mm512_add_ps(redi[1],_mm512_fmadd_ps(zmm8,zmm7,
                                                     _mm512_mul_ps(zmm6,zmm9))); // imp
                             zmm12  = _mm512_loadu_ps(&xre[i+32]);
                             zmm13  = _mm512_loadu_ps(&yre[i+32]);
                             zmm14  = _mm512_loadu_ps(&xim[i+32]);
                             zmm15  = _mm512_loadu_ps(&yim[i+32]);
                             redr[2] = _mm512_add_ps(redr[2],_mm512_fmsub_ps(zmm12,zmm13,
                                                     _mm512_mul_ps(zmm14,zmm15))); // rep
                             redi[2] = _mm512_add_ps(redi[2],_mm512_fmadd_ps(zmm14,zmm13,
                                                     _mm512_mul_ps(zmm12,zmm15))); // imp
                             zmm18 = _mm512_loadu_ps(&xre[i+48]);
                             zmm19 = _mm512_loadu_ps(&yre[i+48]);
                             zmm20 = _mm512_loadu_ps(&xim[i+48]);
                             zmm21 = _mm512_loadu_ps(&yim[i+48]);
                             redr[3] = _mm512_add_ps(redr[3],_mm512_fmsub_ps(zmm18,zmm19,
                                                     _mm512_mul_ps(zmm20,zmm21))); // rep
                             redi[3] = _mm512_add_ps(redi[3],_mm512_fmadd_ps(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21))); // imp
                             zmm24 = _mm512_loadu_ps(&xre[i+64]);
                             zmm25 = _mm512_loadu_ps(&yre[i+64]);
                             zmm26 = _mm512_loadu_ps(&xim[i+64]);
                             zmm27 = _mm512_loadu_ps(&yim[i+64]);
                             redr[4] = _mm512_add_ps(redr[4],_mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27))); // rep
                             redi[4] = _mm512_add_ps(redi[4],_mm512_fmadd_ps(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27))); // imp
                        }

                          redr[0] = _mm512_add_ps(redr[0],redr[2]);
                          redi[0] = _mm512_add_ps(redi[0],redi[2]);
                          redr[1] = _mm512_add_ps(redr[1],redr[3]);
                          redi[1] = _mm512_add_ps(redi[1],redi[3]);
                          redr[0] = _mm512_add_ps(redr[0],redr[4]);
                          redi[0] = _mm512_add_ps(redi[0],redi[4]);

                          for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             redr[0]  = _mm512_add_ps(redr[0],_mm512_fmsub_ps(zmm0,zmm1,
                                                               _mm512_mul_ps(zmm2,zmm3))); // rep
                             redi[0]  = _mm512_add_ps(redi[0],_mm512_fmadd_ps(zmm2,zmm1,
                                                              _mm512_mul_ps(zmm0,zmm3))); // imp
                             zmm6  = _mm512_loadu_ps(&xre[i+16]);
                             zmm7  = _mm512_loadu_ps(&yre[i+16]);
                             zmm8  = _mm512_loadu_ps(&xim[i+16]);
                             zmm9  = _mm512_loadu_ps(&yim[i+16]);
                             redr[1] = _mm512_add_ps(redr[1],_mm512_fmsub_ps(zmm6,zmm7,
                                                             _mm512_mul_ps(zmm8,zmm9))); // rep
                             redi[1] = _mm512_add_ps(redi[1],_mm512_fmadd_ps(zmm8,zmm7,
                                                     _mm512_mul_ps(zmm6,zmm9))); // imp
                        }

                          redr[0] = _mm512_add_ps(redr[0],redr[1]);
                          redi[0] = _mm512_add_ps(redi[0],redi[1]);

                          for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&yre[i+0]);
                             zmm2  = _mm512_loadu_ps(&xim[i+0]);
                             zmm3  = _mm512_loadu_ps(&yim[i+0]);
                             redr[0]  = _mm512_add_ps(redr[0],_mm512_fmsub_ps(zmm0,zmm1,
                                                               _mm512_mul_ps(zmm2,zmm3))); // rep
                             redi[0]  = _mm512_add_ps(redi[0],_mm512_fmadd_ps(zmm2,zmm1,
                                                              _mm512_mul_ps(zmm0,zmm3))); // imp
                        }

                          for(; (i+0) < n; i += 1) {
                              const float xr = xre[i];
                              const float yr = yre[i];
                              const float xi = xim[i];
                              const float yi = yim[i];
                              re += (xr*yr)-(xi*yi);
                              im += (xi*yr)+(xr*yi);
                               
                        }
                          re += _mm512_reduce_add_ps(redr[0]);
                          *mre = re*invN;
                          im += _mm512_reduce_add_ps(redi[0]);
                          *mim = im*invN;
                } 


                
                   void gms::cmath::cmean_prod_u10x_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) yre,
                                                  const float * __restrict __ATTR_ALIGN__(64) yim,
                                                  float * __restrict mre,
                                                  float * __restrict mim,
                                                  const int32_t n) {

                         if(__builtin_expect(0==n,0)) {return;}
                          __m512 zmm0,zmm1,zmm2,zmm3;
                          __m512 zmm4,zmm5,zmm6,zmm7;
                          __m512 zmm8,zmm9,zmm10,zmm11;
                          __m512 zmm12,zmm13,zmm14,zmm15;
                          __m512 zmm16,zmm17,zmm18,zmm19;
                          __m512 zmm20,zmm21,zmm22,zmm23;
                          __m512 zmm24,zmm25,zmm26,zmm27;
                          __m512 zmm28,zmm29,zmm30,zmm31;
                         __ATTR_ALIGN__(64) __m512 redr[10] = {_mm512_setzero_ps()};
                         __ATTR_ALIGN__(64) __m512 redi[10] = {_mm512_setzero_ps()};
                          float re,im;
                         constexpr float invN = 1.0f/static_cast<float>(n);
                         int32_t i; 
                         re = 0.0f;
                         im = 0.0f;
                         for(i = 0; (i+159) < n; i += 160) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+32],_MM_HINT_T0);
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             redr[0]  = _mm512_add_ps(redr[0],_mm512_fmsub_ps(zmm0,zmm1,
                                                              _mm512_mul_ps(zmm2,zmm3))); // rep
                             redi[0]  = _mm512_add_ps(redi[0],_mm512_fmadd_ps(zmm2,zmm1,
                                                              _mm512_mul_ps(zmm0,zmm3))); // imp
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             redr[1] = _mm512_add_ps(redr[1],_mm512_fmsub_ps(zmm6,zmm7,
                                                             _mm512_mul_ps(zmm8,zmm9))); // rep
                             redi[1] = _mm512_add_ps(redi[1],_mm512_fmadd_ps(zmm8,zmm7,
                                                     _mm512_mul_ps(zmm6,zmm9))); // imp
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+64],_MM_HINT_T0);
                             zmm12  = _mm512_load_ps(&xre[i+32]);
                             zmm13  = _mm512_load_ps(&yre[i+32]);
                             zmm14  = _mm512_load_ps(&xim[i+32]);
                             zmm15  = _mm512_load_ps(&yim[i+32]);
                             redr[2] = _mm512_add_ps(redr[2],_mm512_fmsub_ps(zmm12,zmm13,
                                                     _mm512_mul_ps(zmm14,zmm15))); // rep
                             redi[2] = _mm512_add_ps(redi[2],_mm512_fmadd_ps(zmm14,zmm13,
                                                     _mm512_mul_ps(zmm12,zmm15))); // imp
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             redr[3] = _mm512_add_ps(redr[3],_mm512_fmsub_ps(zmm18,zmm19,
                                                     _mm512_mul_ps(zmm20,zmm21))); // rep
                             redi[3] = _mm512_add_ps(redi[3],_mm512_fmadd_ps(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21))); // imp
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+96],_MM_HINT_T0);
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             redr[4] = _mm512_add_ps(redr[3],_mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27))); // rep
                             redi[4] = _mm512_add_ps(redi[3],_mm512_fmadd_ps(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27))); // imp
                             zmm30 = _mm512_load_ps(&xre[i+80]);
                             zmm31 = _mm512_load_ps(&yre[i+80]);
                             zmm0  = _mm512_load_ps(&xim[i+80]);
                             zmm1  = _mm512_load_ps(&yim[i+80]);
                             redr[5]  = _mm512_add_ps(redr[5],_mm512_fmsub_ps(zmm30,zmm31,
                                               _mm512_mul_ps(zmm0,zmm1))); // rep
                             redi[5]  = _mm512_add_ps(redi[5],_mm512_fmadd_ps(zmm0,zmm31,
                                               _mm512_mul_ps(zmm30,zmm1))); // imp
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+128],_MM_HINT_T0);
                             zmm4  = _mm512_load_ps(&xre[i+96]);
                             zmm5  = _mm512_load_ps(&yre[i+96]);
                             zmm6  = _mm512_load_ps(&xim[i+96]);
                             zmm7  = _mm512_load_ps(&yim[i+96]);
                             redr[6]  = _mm512_add_ps(redr[6],_mm512_fmsub_ps(zmm4,zmm5,
                                               _mm512_mul_ps(zmm6,zmm7))); // rep
                             redi[6]  = _mm512_add_ps(redi[6],_mm512_fmadd_ps(zmm6,zmm5,
                                               _mm512_mul_ps(zmm4,zmm7))); // imp
                             zmm10 = _mm512_load_ps(&xre[i+112]);
                             zmm11 = _mm512_load_ps(&yre[i+112]);
                             zmm12 = _mm512_load_ps(&xim[i+112]);
                             zmm13 = _mm512_load_ps(&yim[i+112]);
                             redr[7] = _mm512_add_ps(redr[7],_mm512_fmsub_ps(zmm10,zmm11,
                                               _mm512_mul_ps(zmm12,zmm13))); // rep
                             redi[7] = _mm512_add_ps(redi[7],_mm512_fmadd_ps(zmm13,zmm12,
                                               _mm512_mul_ps(zmm10,zmm14))); // imp
                             _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&yim[i+160],_MM_HINT_T0);
                             zmm16 = _mm512_load_ps(&xre[i+128]);
                             zmm17 = _mm512_load_ps(&yre[i+128]);
                             zmm18 = _mm512_load_ps(&xim[i+128]);
                             zmm19 = _mm512_load_ps(&yim[i+128]);
                             redr[8] = _mm512_add_ps(redr[8],_mm512_fmsub_ps(zmm16,zmm17,
                                               _mm512_mul_ps(zmm18,zmm19))); // rep
                             redi[8] = _mm512_add_ps(redi[8],_mm512_fmadd_ps(zmm18,zmm17,
                                               _mm512_mul_ps(zmm16,zmm19))); // imp
                             zmm22 = _mm512_load_ps(&xre[i+144]);
                             zmm23 = _mm512_load_ps(&yre[i+144]);
                             zmm24 = _mm512_load_ps(&xim[i+144]);
                             zmm25 = _mm512_load_ps(&yim[i+144]);
                             redr[9] = _mm512_add_ps(redr[9],_mm512_fmsub_ps(zmm22,zmm23,
                                               _mm512_mul_ps(zmm24,zmm25))); // rep
                             redi[9] = _mm512_add_ps(redi[9],_mm512_fmadd_ps(zmm24,zmm23,
                                               _mm512_mul_ps(zmm22,zmm25))); // imp
                           
                        }

                             redr[0] = _mm512_add_ps(redr[0],redr[5]);
                             redi[0] = _mm512_add_ps(redi[0],redi[5]);
                             redr[1] = _mm512_add_ps(redr[1],redr[6]);
                             redi[1] = _mm512_add_ps(redi[1],redi[6]);
                             redr[2] = _mm512_add_ps(redr[2],redr[7]);
                             redi[2] = _mm512_add_ps(redi[2],redi[7]);
                             redr[3] = _mm512_add_ps(redr[3],redr[8]);
                             redi[3] = _mm512_add_ps(redi[3],redi[8]);
                             redr[4] = _mm512_add_ps(redr[4],redr[9]);
                             redi[4] = _mm512_add_ps(redi[4],redi[9]); 

                         for(; (i+79) < n; i += 80) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             redr[0]  = _mm512_add_ps(redr[0],_mm512_fmsub_ps(zmm0,zmm1,
                                                               _mm512_mul_ps(zmm2,zmm3))); // rep
                             redi[0]  = _mm512_add_ps(redi[0],_mm512_fmadd_ps(zmm2,zmm1,
                                                              _mm512_mul_ps(zmm0,zmm3))); // imp
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             redr[1] = _mm512_add_ps(redr[1],_mm512_fmsub_ps(zmm6,zmm7,
                                                             _mm512_mul_ps(zmm8,zmm9))); // rep
                             redi[1] = _mm512_add_ps(redi[1],_mm512_fmadd_ps(zmm8,zmm7,
                                                     _mm512_mul_ps(zmm6,zmm9))); // imp
                             zmm12  = _mm512_load_ps(&xre[i+32]);
                             zmm13  = _mm512_load_ps(&yre[i+32]);
                             zmm14  = _mm512_load_ps(&xim[i+32]);
                             zmm15  = _mm512_load_ps(&yim[i+32]);
                             redr[2] = _mm512_add_ps(redr[2],_mm512_fmsub_ps(zmm12,zmm13,
                                                     _mm512_mul_ps(zmm14,zmm15))); // rep
                             redi[2] = _mm512_add_ps(redi[2],_mm512_fmadd_ps(zmm14,zmm13,
                                                     _mm512_mul_ps(zmm12,zmm15))); // imp
                             zmm18 = _mm512_load_ps(&xre[i+48]);
                             zmm19 = _mm512_load_ps(&yre[i+48]);
                             zmm20 = _mm512_load_ps(&xim[i+48]);
                             zmm21 = _mm512_load_ps(&yim[i+48]);
                             redr[3] = _mm512_add_ps(redr[3],_mm512_fmsub_ps(zmm18,zmm19,
                                                     _mm512_mul_ps(zmm20,zmm21))); // rep
                             redi[3] = _mm512_add_ps(redi[3],_mm512_fmadd_ps(zmm20,zmm19,
                                               _mm512_mul_ps(zmm18,zmm21))); // imp
                             zmm24 = _mm512_load_ps(&xre[i+64]);
                             zmm25 = _mm512_load_ps(&yre[i+64]);
                             zmm26 = _mm512_load_ps(&xim[i+64]);
                             zmm27 = _mm512_load_ps(&yim[i+64]);
                             redr[4] = _mm512_add_ps(redr[4],_mm512_fmsub_ps(zmm24,zmm25,
                                               _mm512_mul_ps(zmm26,zmm27))); // rep
                             redi[4] = _mm512_add_ps(redi[4],_mm512_fmadd_ps(zmm26,zmm25,
                                               _mm512_mul_ps(zmm24,zmm27))); // imp
                        }

                          redr[0] = _mm512_add_ps(redr[0],redr[2]);
                          redi[0] = _mm512_add_ps(redi[0],redi[2]);
                          redr[1] = _mm512_add_ps(redr[1],redr[3]);
                          redi[1] = _mm512_add_ps(redi[1],redi[3]);
                          redr[0] = _mm512_add_ps(redr[0],redr[4]);
                          redi[0] = _mm512_add_ps(redi[0],redi[4]);

                          for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             redr[0]  = _mm512_add_ps(redr[0],_mm512_fmsub_ps(zmm0,zmm1,
                                                               _mm512_mul_ps(zmm2,zmm3))); // rep
                             redi[0]  = _mm512_add_ps(redi[0],_mm512_fmadd_ps(zmm2,zmm1,
                                                              _mm512_mul_ps(zmm0,zmm3))); // imp
                             zmm6  = _mm512_load_ps(&xre[i+16]);
                             zmm7  = _mm512_load_ps(&yre[i+16]);
                             zmm8  = _mm512_load_ps(&xim[i+16]);
                             zmm9  = _mm512_load_ps(&yim[i+16]);
                             redr[1] = _mm512_add_ps(redr[1],_mm512_fmsub_ps(zmm6,zmm7,
                                                             _mm512_mul_ps(zmm8,zmm9))); // rep
                             redi[1] = _mm512_add_ps(redi[1],_mm512_fmadd_ps(zmm8,zmm7,
                                                     _mm512_mul_ps(zmm6,zmm9))); // imp
                        }

                          redr[0] = _mm512_add_ps(redr[0],redr[1]);
                          redi[0] = _mm512_add_ps(redi[0],redi[1]);

                          for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&yre[i+0]);
                             zmm2  = _mm512_load_ps(&xim[i+0]);
                             zmm3  = _mm512_load_ps(&yim[i+0]);
                             redr[0]  = _mm512_add_ps(redr[0],_mm512_fmsub_ps(zmm0,zmm1,
                                                               _mm512_mul_ps(zmm2,zmm3))); // rep
                             redi[0]  = _mm512_add_ps(redi[0],_mm512_fmadd_ps(zmm2,zmm1,
                                                              _mm512_mul_ps(zmm0,zmm3))); // imp
                        }

                          for(; (i+0) < n; i += 1) {
                              const float xr = xre[i];
                              const float yr = yre[i];
                              const float xi = xim[i];
                              const float yi = yim[i];
                              re += (xr*yr)-(xi*yi);
                              im += (xi*yr)+(xr*yi);
                               
                        }
                          re += _mm512_reduce_add_ps(redr[0]);
                          *mre = re*invN;
                          im += _mm512_reduce_add_ps(redi[0]);
                          *mim = im*invN;
                } 
  
   
