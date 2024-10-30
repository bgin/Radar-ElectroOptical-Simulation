



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
#include "GMS_carg_vec_zmm16r4.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_cephes.h"




                   void gms::math::cargv_zmm16r4_unroll_16x_a(const float * __restrict __ATTR_ALIGN__(64) re,
                                                   const float * __restrict __ATTR_ALIGN__(64) im,
                                                   float * __restrict __ATTR_ALIGN__(64) carg,
                                                   const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           vfloat zmm0,zmm1,zmm2,zmm3;
                           vfloat zmm4,zmm5,zmm6,zmm7;
                           vfloat zmm8,zmm9,zmm10,zmm11;
                           vfloat zmm12,zmm13,zmm14,zmm15;
                           vfloat zmm16,zmm17,zmm18,zmm19;
                           vfloat zmm20,zmm21,zmm22,zmm23;
                           vfloat zmm24,zmm25,zmm26,zmm27;
                           vfloat zmm28,zmm29,zmm30,zmm31;
                          int32_t i;

                          for(i = 0; (i+255) < n; i += 256) {
                              _mm_prefetch((const char*)&re[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm4  = _mm512_load_ps(&re[i+32]);
                              zmm5  = _mm512_load_ps(&im[i+32]);
                              _mm512_store_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_load_ps(&re[i+48]);
                              zmm7  = _mm512_load_ps(&im[i+48]);
                              _mm512_store_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm8  = _mm512_load_ps(&re[i+64]);
                              zmm9  = _mm512_load_ps(&im[i+64]);
                              _mm512_store_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_load_ps(&re[i+80]);
                              zmm11 = _mm512_load_ps(&im[i+80]);
                              _mm512_store_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                              _mm_prefetch((const char*)&re[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+128],_MM_HINT_T0);
                              zmm12 = _mm512_load_ps(&re[i+96]);
                              zmm13 = _mm512_load_ps(&im[i+96]);
                              _mm512_store_ps(&carg[i+96], xatan2f(zmm12,zmm13));
                              zmm14 = _mm512_load_ps(&re[i+112]);
                              zmm15 = _mm512_load_ps(&im[i+112]);
                              _mm512_store_ps(&carg[i+112], xatan2f(zmm14,zmm15));
                              _mm_prefetch((const char*)&re[i+160],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+160],_MM_HINT_T0);
                              zmm16 = _mm512_load_ps(&re[i+128]);
                              zmm17 = _mm512_load_ps(&im[i+128]);
                              _mm512_store_ps(&carg[i+128], xatan2f(zmm16,zmm17));
                              zmm18 = _mm512_load_ps(&re[i+144]);
                              zmm19 = _mm512_load_ps(&im[i+144]);
                              _mm512_store_ps(&carg[i+144], xatan2f(zmm18,zmm19));
                              _mm_prefetch((const char*)&re[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+192],_MM_HINT_T0);
                              zmm20 = _mm512_load_ps(&re[i+160]);
                              zmm21 = _mm512_load_ps(&im[i+160]);
                              _mm512_store_ps(&carg[i+144], xatan2f(zmm20,zmm21));
                              zmm22 = _mm512_load_ps(&re[i+176]);
                              zmm23 = _mm512_load_ps(&im[i+176]);
                              _mm512_store_ps(&carg[i+176], xatan2f(zmm22,zmm23));
                              _mm_prefetch((const char*)&re[i+224],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+224],_MM_HINT_T0);
                              zmm24 = _mm512_load_ps(&re[i+192]);
                              zmm25 = _mm512_load_ps(&im[i+192]);
                              _mm512_store_ps(&carg[i+192], xatan2f(zmm24,zmm25));
                              zmm26 = _mm512_load_ps(&re[i+208]);
                              zmm27 = _mm512_load_ps(&im[i+208]);
                              _mm512_store_ps(&carg[i+208], xatan2f(zmm26,zmm27));
                              _mm_prefetch((const char*)&re[i+256],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+256],_MM_HINT_T0);
                              zmm28 = _mm512_load_ps(&re[i+224]);
                              zmm29 = _mm512_load_ps(&im[i+224]);
                              _mm512_store_ps(&carg[i+224], xatan2f(zmm28,zmm29));
                              zmm30 = _mm512_load_ps(&re[i+240]);
                              zmm31 = _mm512_load_ps(&im[i+240]);
                              _mm512_store_ps(&carg[i+240], xatan2f(zmm30,zmm31));
                        }

                         for(; (i+191) < n; i += 192) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_load_ps(&re[i+32]);
                              zmm5  = _mm512_load_ps(&im[i+32]);
                              _mm512_store_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_load_ps(&re[i+48]);
                              zmm7  = _mm512_load_ps(&im[i+48]);
                              _mm512_store_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              zmm8  = _mm512_load_ps(&re[i+64]);
                              zmm9  = _mm512_load_ps(&im[i+64]);
                              _mm512_store_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_load_ps(&re[i+80]);
                              zmm11 = _mm512_load_ps(&im[i+80]);
                              _mm512_store_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                              zmm12 = _mm512_load_ps(&re[i+96]);
                              zmm13 = _mm512_load_ps(&im[i+96]);
                              _mm512_store_ps(&carg[i+96], xatan2f(zmm12,zmm13));
                              zmm14 = _mm512_load_ps(&re[i+112]);
                              zmm15 = _mm512_load_ps(&im[i+112]);
                              _mm512_store_ps(&carg[i+112], xatan2f(zmm14,zmm15));
                              zmm16 = _mm512_load_ps(&re[i+128]);
                              zmm17 = _mm512_load_ps(&im[i+128]);
                              _mm512_store_ps(&carg[i+128], xatan2f(zmm16,zmm17));
                              zmm18 = _mm512_load_ps(&re[i+144]);
                              zmm19 = _mm512_load_ps(&im[i+144]);
                              _mm512_store_ps(&carg[i+144], xatan2f(zmm18,zmm19));
                              zmm20 = _mm512_load_ps(&re[i+160]);
                              zmm21 = _mm512_load_ps(&im[i+160]);
                              _mm512_store_ps(&carg[i+144], xatan2f(zmm20,zmm21));
                              zmm22 = _mm512_load_ps(&re[i+176]);
                              zmm23 = _mm512_load_ps(&im[i+176]);
                              _mm512_store_ps(&carg[i+176], xatan2f(zmm22,zmm23));
                        }

                       for(; (i+127) < n; i += 128) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_load_ps(&re[i+32]);
                              zmm5  = _mm512_load_ps(&im[i+32]);
                              _mm512_store_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_load_ps(&re[i+48]);
                              zmm7  = _mm512_load_ps(&im[i+48]);
                              _mm512_store_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              zmm8  = _mm512_load_ps(&re[i+64]);
                              zmm9  = _mm512_load_ps(&im[i+64]);
                              _mm512_store_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_load_ps(&re[i+80]);
                              zmm11 = _mm512_load_ps(&im[i+80]);
                              _mm512_store_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                              zmm12 = _mm512_load_ps(&re[i+96]);
                              zmm13 = _mm512_load_ps(&im[i+96]);
                              _mm512_store_ps(&carg[i+96], xatan2f(zmm12,zmm13));
                              zmm14 = _mm512_load_ps(&re[i+112]);
                              zmm15 = _mm512_load_ps(&im[i+112]);
                              _mm512_store_ps(&carg[i+112], xatan2f(zmm14,zmm15));
                       }

                       for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_load_ps(&re[i+32]);
                              zmm5  = _mm512_load_ps(&im[i+32]);
                              _mm512_store_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_load_ps(&re[i+48]);
                              zmm7  = _mm512_load_ps(&im[i+48]);
                              _mm512_store_ps(&carg[i+48], xatan2f(zmm6,zmm7));   
                       }

                       for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                      }

                       for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                      }

                       for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float xim = im[i];
                              carg[i]         = ceph_atan2f(xre,xim);
                      }
               }


                 
                   void gms::math::cargv_zmm16r4_unroll_16x_u(const float * __restrict  re,
                                                   const float * __restrict  im,
                                                   float * __restrict  carg,
                                                   const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           vfloat zmm0,zmm1,zmm2,zmm3;
                           vfloat zmm4,zmm5,zmm6,zmm7;
                           vfloat zmm8,zmm9,zmm10,zmm11;
                           vfloat zmm12,zmm13,zmm14,zmm15;
                           vfloat zmm16,zmm17,zmm18,zmm19;
                           vfloat zmm20,zmm21,zmm22,zmm23;
                           vfloat zmm24,zmm25,zmm26,zmm27;
                           vfloat zmm28,zmm29,zmm30,zmm31;
                          int32_t i;

                          for(i = 0; (i+255) < n; i += 256) {
                              _mm_prefetch((const char*)&re[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm4  = _mm512_loadu_ps(&re[i+32]);
                              zmm5  = _mm512_loadu_ps(&im[i+32]);
                              _mm512_storeu_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_loadu_ps(&re[i+48]);
                              zmm7  = _mm512_loadu_ps(&im[i+48]);
                              _mm512_storeu_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm8  = _mm512_loadu_ps(&re[i+64]);
                              zmm9  = _mm512_loadu_ps(&im[i+64]);
                              _mm512_storeu_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_loadu_ps(&re[i+80]);
                              zmm11 = _mm512_loadu_ps(&im[i+80]);
                              _mm512_storeu_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                              _mm_prefetch((const char*)&re[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+128],_MM_HINT_T0);
                              zmm12 = _mm512_loadu_ps(&re[i+96]);
                              zmm13 = _mm512_loadu_ps(&im[i+96]);
                              _mm512_storeu_ps(&carg[i+96], xatan2f(zmm12,zmm13));
                              zmm14 = _mm512_loadu_ps(&re[i+112]);
                              zmm15 = _mm512_loadu_ps(&im[i+112]);
                              _mm512_storeu_ps(&carg[i+112], xatan2f(zmm14,zmm15));
                              _mm_prefetch((const char*)&re[i+160],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+160],_MM_HINT_T0);
                              zmm16 = _mm512_loadu_ps(&re[i+128]);
                              zmm17 = _mm512_loadu_ps(&im[i+128]);
                              _mm512_storeu_ps(&carg[i+128], xatan2f(zmm16,zmm17));
                              zmm18 = _mm512_loadu_ps(&re[i+144]);
                              zmm19 = _mm512_loadu_ps(&im[i+144]);
                              _mm512_storeu_ps(&carg[i+144], xatan2f(zmm18,zmm19));
                              _mm_prefetch((const char*)&re[i+192],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+192],_MM_HINT_T0);
                              zmm20 = _mm512_loadu_ps(&re[i+160]);
                              zmm21 = _mm512_loadu_ps(&im[i+160]);
                              _mm512_storeu_ps(&carg[i+144], xatan2f(zmm20,zmm21));
                              zmm22 = _mm512_loadu_ps(&re[i+176]);
                              zmm23 = _mm512_loadu_ps(&im[i+176]);
                              _mm512_storeu_ps(&carg[i+176], xatan2f(zmm22,zmm23));
                              _mm_prefetch((const char*)&re[i+224],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+224],_MM_HINT_T0);
                              zmm24 = _mm512_loadu_ps(&re[i+192]);
                              zmm25 = _mm512_loadu_ps(&im[i+192]);
                              _mm512_storeu_ps(&carg[i+192], xatan2f(zmm24,zmm25));
                              zmm26 = _mm512_loadu_ps(&re[i+208]);
                              zmm27 = _mm512_loadu_ps(&im[i+208]);
                              _mm512_storeu_ps(&carg[i+208], xatan2f(zmm26,zmm27));
                              _mm_prefetch((const char*)&re[i+256],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+256],_MM_HINT_T0);
                              zmm28 = _mm512_loadu_ps(&re[i+224]);
                              zmm29 = _mm512_loadu_ps(&im[i+224]);
                              _mm512_storeu_ps(&carg[i+224], xatan2f(zmm28,zmm29));
                              zmm30 = _mm512_loadu_ps(&re[i+240]);
                              zmm31 = _mm512_loadu_ps(&im[i+240]);
                              _mm512_storeu_ps(&carg[i+240], xatan2f(zmm30,zmm31));
                        }

                         for(; (i+191) < n; i += 192) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_loadu_ps(&re[i+32]);
                              zmm5  = _mm512_loadu_ps(&im[i+32]);
                              _mm512_storeu_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_loadu_ps(&re[i+48]);
                              zmm7  = _mm512_loadu_ps(&im[i+48]);
                              _mm512_storeu_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              zmm8  = _mm512_loadu_ps(&re[i+64]);
                              zmm9  = _mm512_loadu_ps(&im[i+64]);
                              _mm512_storeu_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_loadu_ps(&re[i+80]);
                              zmm11 = _mm512_loadu_ps(&im[i+80]);
                              _mm512_storeu_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                              zmm12 = _mm512_loadu_ps(&re[i+96]);
                              zmm13 = _mm512_loadu_ps(&im[i+96]);
                              _mm512_storeu_ps(&carg[i+96], xatan2f(zmm12,zmm13));
                              zmm14 = _mm512_loadu_ps(&re[i+112]);
                              zmm15 = _mm512_loadu_ps(&im[i+112]);
                              _mm512_storeu_ps(&carg[i+112], xatan2f(zmm14,zmm15));
                              zmm16 = _mm512_loadu_ps(&re[i+128]);
                              zmm17 = _mm512_loadu_ps(&im[i+128]);
                              _mm512_storeu_ps(&carg[i+128], xatan2f(zmm16,zmm17));
                              zmm18 = _mm512_loadu_ps(&re[i+144]);
                              zmm19 = _mm512_loadu_ps(&im[i+144]);
                              _mm512_storeu_ps(&carg[i+144], xatan2f(zmm18,zmm19));
                              zmm20 = _mm512_loadu_ps(&re[i+160]);
                              zmm21 = _mm512_loadu_ps(&im[i+160]);
                              _mm512_storeu_ps(&carg[i+144], xatan2f(zmm20,zmm21));
                              zmm22 = _mm512_loadu_ps(&re[i+176]);
                              zmm23 = _mm512_loadu_ps(&im[i+176]);
                              _mm512_storeu_ps(&carg[i+176], xatan2f(zmm22,zmm23));
                        }

                       for(; (i+127) < n; i += 128) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_loadu_ps(&re[i+32]);
                              zmm5  = _mm512_loadu_ps(&im[i+32]);
                              _mm512_storeu_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_loadu_ps(&re[i+48]);
                              zmm7  = _mm512_loadu_ps(&im[i+48]);
                              _mm512_storeu_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              zmm8  = _mm512_loadu_ps(&re[i+64]);
                              zmm9  = _mm512_loadu_ps(&im[i+64]);
                              _mm512_storeu_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_loadu_ps(&re[i+80]);
                              zmm11 = _mm512_loadu_ps(&im[i+80]);
                              _mm512_storeu_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                              zmm12 = _mm512_loadu_ps(&re[i+96]);
                              zmm13 = _mm512_loadu_ps(&im[i+96]);
                              _mm512_storeu_ps(&carg[i+96], xatan2f(zmm12,zmm13));
                              zmm14 = _mm512_loadu_ps(&re[i+112]);
                              zmm15 = _mm512_loadu_ps(&im[i+112]);
                              _mm512_storeu_ps(&carg[i+112], xatan2f(zmm14,zmm15));
                       }

                       for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_loadu_ps(&re[i+32]);
                              zmm5  = _mm512_loadu_ps(&im[i+32]);
                              _mm512_storeu_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_loadu_ps(&re[i+48]);
                              zmm7  = _mm512_loadu_ps(&im[i+48]);
                              _mm512_storeu_ps(&carg[i+48], xatan2f(zmm6,zmm7));   
                       }

                       for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                      }

                       for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                      }

                       for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float xim = im[i];
                              carg[i]         = ceph_atan2f(xre,xim);
                      }
               }


                
                   void gms::math::cargv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) re,
                                                   const float * __restrict __ATTR_ALIGN__(64) im,
                                                   float * __restrict __ATTR_ALIGN__(64) carg,
                                                   const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           vfloat zmm0,zmm1,zmm2,zmm3;
                           vfloat zmm4,zmm5,zmm6,zmm7;
                           vfloat zmm8,zmm9,zmm10,zmm11;
                           vfloat zmm12,zmm13,zmm14,zmm15;
                           vfloat zmm16,zmm17,zmm18,zmm19;
                          int32_t i;

                          for(i = 0; (i+159) < n; i += 160) {
                              _mm_prefetch((const char*)&re[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm4  = _mm512_load_ps(&re[i+32]);
                              zmm5  = _mm512_load_ps(&im[i+32]);
                              _mm512_store_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_load_ps(&re[i+48]);
                              zmm7  = _mm512_load_ps(&im[i+48]);
                              _mm512_store_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm8  = _mm512_load_ps(&re[i+64]);
                              zmm9  = _mm512_load_ps(&im[i+64]);
                              _mm512_store_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_load_ps(&re[i+80]);
                              zmm11 = _mm512_load_ps(&im[i+80]);
                              _mm512_store_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                              _mm_prefetch((const char*)&re[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+128],_MM_HINT_T0);
                              zmm12 = _mm512_load_ps(&re[i+96]);
                              zmm13 = _mm512_load_ps(&im[i+96]);
                              _mm512_store_ps(&carg[i+96], xatan2f(zmm12,zmm13));
                              zmm14 = _mm512_load_ps(&re[i+112]);
                              zmm15 = _mm512_load_ps(&im[i+112]);
                              _mm512_store_ps(&carg[i+112], xatan2f(zmm14,zmm15));
                              _mm_prefetch((const char*)&re[i+160],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+160],_MM_HINT_T0);
                              zmm16 = _mm512_load_ps(&re[i+128]);
                              zmm17 = _mm512_load_ps(&im[i+128]);
                              _mm512_store_ps(&carg[i+128], xatan2f(zmm16,zmm17));
                              zmm18 = _mm512_load_ps(&re[i+144]);
                              zmm19 = _mm512_load_ps(&im[i+144]);
                              _mm512_store_ps(&carg[i+144], xatan2f(zmm18,zmm19));
                         }

                       

                       for(; (i+127) < n; i += 128) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_load_ps(&re[i+32]);
                              zmm5  = _mm512_load_ps(&im[i+32]);
                              _mm512_store_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_load_ps(&re[i+48]);
                              zmm7  = _mm512_load_ps(&im[i+48]);
                              _mm512_store_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              zmm8  = _mm512_load_ps(&re[i+64]);
                              zmm9  = _mm512_load_ps(&im[i+64]);
                              _mm512_store_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_load_ps(&re[i+80]);
                              zmm11 = _mm512_load_ps(&im[i+80]);
                              _mm512_store_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                              zmm12 = _mm512_load_ps(&re[i+96]);
                              zmm13 = _mm512_load_ps(&im[i+96]);
                              _mm512_store_ps(&carg[i+96], xatan2f(zmm12,zmm13));
                              zmm14 = _mm512_load_ps(&re[i+112]);
                              zmm15 = _mm512_load_ps(&im[i+112]);
                              _mm512_store_ps(&carg[i+112], xatan2f(zmm14,zmm15));
                       }

                       for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_load_ps(&re[i+32]);
                              zmm5  = _mm512_load_ps(&im[i+32]);
                              _mm512_store_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_load_ps(&re[i+48]);
                              zmm7  = _mm512_load_ps(&im[i+48]);
                              _mm512_store_ps(&carg[i+48], xatan2f(zmm6,zmm7));   
                       }

                       for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                      }

                       for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                      }

                       for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float xim = im[i];
                              carg[i]         = ceph_atan2f(xre,xim);
                      }
               }


                 
                   void gms::math::cargv_zmm16r4_unroll_10x_u(const float * __restrict  re,
                                                   const float * __restrict  im,
                                                   float * __restrict  carg,
                                                   const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           vfloat zmm0,zmm1,zmm2,zmm3;
                           vfloat zmm4,zmm5,zmm6,zmm7;
                           vfloat zmm8,zmm9,zmm10,zmm11;
                           vfloat zmm12,zmm13,zmm14,zmm15;
                           vfloat zmm16,zmm17,zmm18,zmm19;
                          int32_t i;

                          for(i = 0; (i+159) < n; i += 160) {
                              _mm_prefetch((const char*)&re[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm4  = _mm512_loadu_ps(&re[i+32]);
                              zmm5  = _mm512_loadu_ps(&im[i+32]);
                              _mm512_storeu_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_loadu_ps(&re[i+48]);
                              zmm7  = _mm512_loadu_ps(&im[i+48]);
                              _mm512_storeu_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm8  = _mm512_loadu_ps(&re[i+64]);
                              zmm9  = _mm512_loadu_ps(&im[i+64]);
                              _mm512_storeu_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_loadu_ps(&re[i+80]);
                              zmm11 = _mm512_loadu_ps(&im[i+80]);
                              _mm512_storeu_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                              _mm_prefetch((const char*)&re[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+128],_MM_HINT_T0);
                              zmm12 = _mm512_loadu_ps(&re[i+96]);
                              zmm13 = _mm512_loadu_ps(&im[i+96]);
                              _mm512_storeu_ps(&carg[i+96], xatan2f(zmm12,zmm13));
                              zmm14 = _mm512_loadu_ps(&re[i+112]);
                              zmm15 = _mm512_loadu_ps(&im[i+112]);
                              _mm512_storeu_ps(&carg[i+112], xatan2f(zmm14,zmm15));
                              _mm_prefetch((const char*)&re[i+160],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+160],_MM_HINT_T0);
                              zmm16 = _mm512_loadu_ps(&re[i+128]);
                              zmm17 = _mm512_loadu_ps(&im[i+128]);
                              _mm512_storeu_ps(&carg[i+128], xatan2f(zmm16,zmm17));
                              zmm18 = _mm512_loadu_ps(&re[i+144]);
                              zmm19 = _mm512_loadu_ps(&im[i+144]);
                              _mm512_storeu_ps(&carg[i+144], xatan2f(zmm18,zmm19));
                             
                        }

                       

                       for(; (i+127) < n; i += 128) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_loadu_ps(&re[i+32]);
                              zmm5  = _mm512_loadu_ps(&im[i+32]);
                              _mm512_storeu_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_loadu_ps(&re[i+48]);
                              zmm7  = _mm512_loadu_ps(&im[i+48]);
                              _mm512_storeu_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              zmm8  = _mm512_loadu_ps(&re[i+64]);
                              zmm9  = _mm512_loadu_ps(&im[i+64]);
                              _mm512_storeu_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_loadu_ps(&re[i+80]);
                              zmm11 = _mm512_loadu_ps(&im[i+80]);
                              _mm512_storeu_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                              zmm12 = _mm512_loadu_ps(&re[i+96]);
                              zmm13 = _mm512_loadu_ps(&im[i+96]);
                              _mm512_storeu_ps(&carg[i+96], xatan2f(zmm12,zmm13));
                              zmm14 = _mm512_loadu_ps(&re[i+112]);
                              zmm15 = _mm512_loadu_ps(&im[i+112]);
                              _mm512_storeu_ps(&carg[i+112], xatan2f(zmm14,zmm15));
                       }

                       for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_loadu_ps(&re[i+32]);
                              zmm5  = _mm512_loadu_ps(&im[i+32]);
                              _mm512_storeu_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_loadu_ps(&re[i+48]);
                              zmm7  = _mm512_loadu_ps(&im[i+48]);
                              _mm512_storeu_ps(&carg[i+48], xatan2f(zmm6,zmm7));   
                       }

                       for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                      }

                       for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                      }

                       for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float xim = im[i];
                              carg[i]         = ceph_atan2f(xre,xim);
                      }
               }


                
                   void gms::math::cargv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) re,
                                                   const float * __restrict __ATTR_ALIGN__(64) im,
                                                   float * __restrict __ATTR_ALIGN__(64) carg,
                                                   const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           vfloat zmm0,zmm1,zmm2,zmm3;
                           vfloat zmm4,zmm5,zmm6,zmm7;
                           vfloat zmm8,zmm9,zmm10,zmm11;
                          int32_t i;

                          for(i = 0; (i+95) < n; i += 96) {
                              _mm_prefetch((const char*)&re[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm4  = _mm512_load_ps(&re[i+32]);
                              zmm5  = _mm512_load_ps(&im[i+32]);
                              _mm512_store_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_load_ps(&re[i+48]);
                              zmm7  = _mm512_load_ps(&im[i+48]);
                              _mm512_store_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm8  = _mm512_load_ps(&re[i+64]);
                              zmm9  = _mm512_load_ps(&im[i+64]);
                              _mm512_store_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_load_ps(&re[i+80]);
                              zmm11 = _mm512_load_ps(&im[i+80]);
                              _mm512_store_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                              
                         }

                      for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_load_ps(&re[i+32]);
                              zmm5  = _mm512_load_ps(&im[i+32]);
                              _mm512_store_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_load_ps(&re[i+48]);
                              zmm7  = _mm512_load_ps(&im[i+48]);
                              _mm512_store_ps(&carg[i+48], xatan2f(zmm6,zmm7));   
                       }

                       for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_load_ps(&re[i+16]);
                              zmm3  = _mm512_load_ps(&im[i+16]);
                              _mm512_store_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                      }

                       for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_load_ps(&re[i+0]);
                              zmm1  = _mm512_load_ps(&im[i+0]);
                              _mm512_store_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                      }

                       for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float xim = im[i];
                              carg[i]         = ceph_atan2f(xre,xim);
                      }
               }


                  
                   void gms::math::cargv_zmm16r4_unroll_6x_u(const float * __restrict  re,
                                                   const float * __restrict  im,
                                                   float * __restrict  carg,
                                                   const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           vfloat zmm0,zmm1,zmm2,zmm3;
                           vfloat zmm4,zmm5,zmm6,zmm7;
                           vfloat zmm8,zmm9,zmm10,zmm11;
                          int32_t i;

                          for(i = 0; (i+95) < n; i += 96) {
                              _mm_prefetch((const char*)&re[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+32],_MM_HINT_T0);
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              _mm_prefetch((const char*)&re[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+64],_MM_HINT_T0);
                              zmm4  = _mm512_loadu_ps(&re[i+32]);
                              zmm5  = _mm512_loadu_ps(&im[i+32]);
                              _mm512_storeu_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_loadu_ps(&re[i+48]);
                              zmm7  = _mm512_loadu_ps(&im[i+48]);
                              _mm512_storeu_ps(&carg[i+48], xatan2f(zmm6,zmm7));
                              _mm_prefetch((const char*)&re[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&im[i+96],_MM_HINT_T0);
                              zmm8  = _mm512_loadu_ps(&re[i+64]);
                              zmm9  = _mm512_loadu_ps(&im[i+64]);
                              _mm512_storeu_ps(&carg[i+64], xatan2f(zmm8,zmm9));
                              zmm10 = _mm512_loadu_ps(&re[i+80]);
                              zmm11 = _mm512_loadu_ps(&im[i+80]);
                              _mm512_storeu_ps(&carg[i+80], xatan2f(zmm10,zmm11));
                                                           
                        }

                       

                      for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                              zmm4  = _mm512_loadu_ps(&re[i+32]);
                              zmm5  = _mm512_loadu_ps(&im[i+32]);
                              _mm512_storeu_ps(&carg[i+32], xatan2f(zmm4,zmm5));
                              zmm6  = _mm512_loadu_ps(&re[i+48]);
                              zmm7  = _mm512_loadu_ps(&im[i+48]);
                              _mm512_storeu_ps(&carg[i+48], xatan2f(zmm6,zmm7));   
                       }

                       for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                              zmm2  = _mm512_loadu_ps(&re[i+16]);
                              zmm3  = _mm512_loadu_ps(&im[i+16]);
                              _mm512_storeu_ps(&carg[i+16], xatan2f(zmm2,zmm3));
                      }

                       for(; (i+15) < n; i += 16) {
                              zmm0  = _mm512_loadu_ps(&re[i+0]);
                              zmm1  = _mm512_loadu_ps(&im[i+0]);
                              _mm512_storeu_ps(&carg[i+0], xatan2f(zmm0,zmm1));
                      }

                       for(; (i+0) < n; i += 1) {
                              const float xre = re[i];
                              const float xim = im[i];
                              carg[i]         = ceph_atan2f(xre,xim);
                      }
               }



