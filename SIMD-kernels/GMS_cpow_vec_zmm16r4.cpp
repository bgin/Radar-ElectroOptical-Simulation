


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
#include "GMS_cpow_vec_zmm16r4.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_cephes.h"





                
                   void gms::math::cpowv_zmm16r4_unroll_8x_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict vn,
                                                  float * __restrict cpowr,
                                                  float * __restrict cpowi,
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
                           vfloat zmmx;
                          int32_t i;

                          for(i = 0; (i+127) < n; i += 128) {
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+32], _MM_HINT_T0);
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_loadu_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_loadu_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_loadu_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_storeu_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_storeu_ps(&cpowi[i+16],zmm19);
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+64], _MM_HINT_T0);
                              zmm20 = _mm512_loadu_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_loadu_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_loadu_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_storeu_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_storeu_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_loadu_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_loadu_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_loadu_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_storeu_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_storeu_ps(&cpowi[i+48],zmm7);
                              _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+96], _MM_HINT_T0);
                              zmm8  = _mm512_loadu_ps(&xre[i+64]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_loadu_ps(&xim[i+64]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm9,zmm11));
                              zmmx  = _mm512_loadu_ps(&vn[i+64]);
                              zmm13 = xatanf(_mm512_div_ps(zmm10,zmm8));
                              zmm14 = xpowf(zmm12,zmmx);
                              zmm15 = _mm512_mul_ps(zmmx,zmm13);
                              zmm16 = _mm512_mul_ps(zmm14,xcosf(zmm15));
                              _mm512_storeu_ps(&cpowr[i+64],zmm16);
                              zmm17 = _mm512_mul_ps(zmm14,xsinf(zmm15));
                              _mm512_storeu_ps(&cpowi[i+64],zmm17);
                              zmm18  = _mm512_loadu_ps(&xre[i+80]);
                              zmm19  = _mm512_mul_ps(zmm18,zmm18);
                              zmm20 = _mm512_loadu_ps(&xim[i+80]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_sqrt_ps(_mm512_add_ps(zmm19,zmm21));
                              zmmx  = _mm512_loadu_ps(&vn[i+80]);
                              zmm23 = xatanf(_mm512_div_ps(zmm20,zmm18));
                              zmm24 = xpowf(zmm22,zmmx);
                              zmm25 = _mm512_mul_ps(zmmx,zmm23);
                              zmm26 = _mm512_mul_ps(zmm24,xcosf(zmm25));
                              _mm512_storeu_ps(&cpowr[i+80],zmm26);
                              zmm27 = _mm512_mul_ps(zmm24,xsinf(zmm25));
                              _mm512_storeu_ps(&cpowi[i+80],zmm27);
                              _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+128], _MM_HINT_T0);
                              zmm28 = _mm512_loadu_ps(&xre[i+96]);
                              zmm29 = _mm512_mul_ps(zmm28,zmm28);
                              zmm30 = _mm512_loadu_ps(&xim[i+96]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_sqrt_ps(_mm512_add_ps(zmm29,zmm31));
                              zmmx  = _mm512_loadu_ps(&vn[i+96]);
                              zmm1  = xatanf(_mm512_div_ps(zmm30,zmm28));
                              zmm2  = xpowf(zmm0,zmmx);
                              zmm3  = _mm512_mul_ps(zmmx,zmm1);
                              zmm4  = _mm512_mul_ps(zmm2,xcosf(zmm3));
                              _mm512_storeu_ps(&cpowr[i+96],zmm4);
                              zmm5  = _mm512_mul_ps(zmm2,xsinf(zmm3));
                              _mm512_storeu_ps(&cpowi[i+96],zmm5);
                              zmm6  = _mm512_loadu_ps(&xre[i+112]);
                              zmm7  = _mm512_mul_ps(zmm6,zmm6);
                              zmm8  = _mm512_loadu_ps(&xim[i+112]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_sqrt_ps(_mm512_add_ps(zmm7,zmm9));
                              zmmx  = _mm512_loadu_ps(&vn[i+112]);
                              zmm11 = xatanf(_mm512_div_ps(zmm8,zmm6));
                              zmm12 = xpowf(zmm10,zmmx);
                              zmm13 = _mm512_mul_ps(zmmx,zmm11);
                              zmm14 = _mm512_mul_ps(zmm12,xcosf(zmm13));
                              _mm512_storeu_ps(&cpowr[i+112],zmm14);
                              zmm15 = _mm512_mul_ps(zmm12,xsinf(zmm13));
                              _mm512_storeu_ps(&cpowi[i+112],zmm15);
                         }

                          for(; (i+95) < n; i += 96) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_loadu_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_loadu_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_loadu_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_storeu_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_storeu_ps(&cpowi[i+16],zmm19);
                              zmm20 = _mm512_loadu_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_loadu_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_loadu_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_storeu_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_storeu_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_loadu_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_loadu_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_loadu_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_storeu_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_storeu_ps(&cpowi[i+48],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+64]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_loadu_ps(&xim[i+64]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm9,zmm11));
                              zmmx  = _mm512_loadu_ps(&vn[i+64]);
                              zmm13 = xatanf(_mm512_div_ps(zmm10,zmm8));
                              zmm14 = xpowf(zmm12,zmmx);
                              zmm15 = _mm512_mul_ps(zmmx,zmm13);
                              zmm16 = _mm512_mul_ps(zmm14,xcosf(zmm15));
                              _mm512_storeu_ps(&cpowr[i+64],zmm16);
                              zmm17 = _mm512_mul_ps(zmm14,xsinf(zmm15));
                              _mm512_storeu_ps(&cpowi[i+64],zmm17);
                              zmm18  = _mm512_loadu_ps(&xre[i+80]);
                              zmm19  = _mm512_mul_ps(zmm18,zmm18);
                              zmm20 = _mm512_loadu_ps(&xim[i+80]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_sqrt_ps(_mm512_add_ps(zmm19,zmm21));
                              zmmx  = _mm512_loadu_ps(&vn[i+80]);
                              zmm23 = xatanf(_mm512_div_ps(zmm20,zmm18));
                              zmm24 = xpowf(zmm22,zmmx);
                              zmm25 = _mm512_mul_ps(zmmx,zmm23);
                              zmm26 = _mm512_mul_ps(zmm24,xcosf(zmm25));
                              _mm512_storeu_ps(&cpowr[i+80],zmm26);
                              zmm27 = _mm512_mul_ps(zmm24,xsinf(zmm25));
                              _mm512_storeu_ps(&cpowi[i+80],zmm27);
                         }

                          for(; (i+79) < n; i += 80) {
                               zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_loadu_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_loadu_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_loadu_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_storeu_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_storeu_ps(&cpowi[i+16],zmm19);
                              zmm20 = _mm512_loadu_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_loadu_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_loadu_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_storeu_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_storeu_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_loadu_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_loadu_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_loadu_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_storeu_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_storeu_ps(&cpowi[i+48],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+64]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_loadu_ps(&xim[i+64]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm9,zmm11));
                              zmmx  = _mm512_loadu_ps(&vn[i+64]);
                              zmm13 = xatanf(_mm512_div_ps(zmm10,zmm8));
                              zmm14 = xpowf(zmm12,zmmx);
                              zmm15 = _mm512_mul_ps(zmmx,zmm13);
                              zmm16 = _mm512_mul_ps(zmm14,xcosf(zmm15));
                              _mm512_storeu_ps(&cpowr[i+64],zmm16);
                              zmm17 = _mm512_mul_ps(zmm14,xsinf(zmm15));
                              _mm512_storeu_ps(&cpowi[i+64],zmm17);
                         }

                          for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_loadu_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_loadu_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_loadu_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_storeu_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_storeu_ps(&cpowi[i+16],zmm19);
                              zmm20 = _mm512_loadu_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_loadu_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_loadu_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_storeu_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_storeu_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_loadu_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_loadu_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_loadu_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_storeu_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_storeu_ps(&cpowi[i+48],zmm7);
                         }

                          for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_loadu_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_loadu_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_loadu_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_storeu_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_storeu_ps(&cpowi[i+16],zmm19);
                         }

                          for(; (i+15) < n; i += 15) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                         }

                           for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = z0*z0;
                              const float z2 = xim[i];
                              const float z3 = z2*z2;
                              const float z4 = ceph_sqrtf(z1+z3);
                              const float zx = vn[i];
                              const float z5 = ceph_atanf(z2/z0);
                              const float z6 = ceph_powf(z4,zx);
                              const float z7 = zx*z5;
                              const float z8 = z6*ceph_cosf(z7);
                              cpowr[i]       = z8;
                              const float z9 = z6*ceph_sinf(z7);
                              cpowi[i]       = z9;
                         }
                 }


               
                   void gms::math::cpowv_zmm16r4_unroll_8x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) vn,
                                                  float * __restrict __ATTR_ALIGN__(64) cpowr,
                                                  float * __restrict __ATTR_ALIGN__(64) cpowi,
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
                           vfloat zmmx;
                          int32_t i;

                          for(i = 0; (i+127) < n; i += 128) {
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+32], _MM_HINT_T0);
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_load_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_load_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_load_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&cpowi[i+16],zmm19);
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+64], _MM_HINT_T0);
                              zmm20 = _mm512_load_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_load_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_load_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_store_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_store_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_load_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_load_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_load_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&cpowi[i+48],zmm7);
                              _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+96], _MM_HINT_T0);
                              zmm8  = _mm512_load_ps(&xre[i+64]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_load_ps(&xim[i+64]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm9,zmm11));
                              zmmx  = _mm512_load_ps(&vn[i+64]);
                              zmm13 = xatanf(_mm512_div_ps(zmm10,zmm8));
                              zmm14 = xpowf(zmm12,zmmx);
                              zmm15 = _mm512_mul_ps(zmmx,zmm13);
                              zmm16 = _mm512_mul_ps(zmm14,xcosf(zmm15));
                              _mm512_store_ps(&cpowr[i+64],zmm16);
                              zmm17 = _mm512_mul_ps(zmm14,xsinf(zmm15));
                              _mm512_store_ps(&cpowi[i+64],zmm17);
                              zmm18  = _mm512_load_ps(&xre[i+80]);
                              zmm19  = _mm512_mul_ps(zmm18,zmm18);
                              zmm20 = _mm512_load_ps(&xim[i+80]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_sqrt_ps(_mm512_add_ps(zmm19,zmm21));
                              zmmx  = _mm512_load_ps(&vn[i+80]);
                              zmm23 = xatanf(_mm512_div_ps(zmm20,zmm18));
                              zmm24 = xpowf(zmm22,zmmx);
                              zmm25 = _mm512_mul_ps(zmmx,zmm23);
                              zmm26 = _mm512_mul_ps(zmm24,xcosf(zmm25));
                              _mm512_store_ps(&cpowr[i+80],zmm26);
                              zmm27 = _mm512_mul_ps(zmm24,xsinf(zmm25));
                              _mm512_store_ps(&cpowi[i+80],zmm27);
                              _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+128], _MM_HINT_T0);
                              zmm28 = _mm512_load_ps(&xre[i+96]);
                              zmm29 = _mm512_mul_ps(zmm28,zmm28);
                              zmm30 = _mm512_load_ps(&xim[i+96]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_sqrt_ps(_mm512_add_ps(zmm29,zmm31));
                              zmmx  = _mm512_load_ps(&vn[i+96]);
                              zmm1  = xatanf(_mm512_div_ps(zmm30,zmm28));
                              zmm2  = xpowf(zmm0,zmmx);
                              zmm3  = _mm512_mul_ps(zmmx,zmm1);
                              zmm4  = _mm512_mul_ps(zmm2,xcosf(zmm3));
                              _mm512_store_ps(&cpowr[i+96],zmm4);
                              zmm5  = _mm512_mul_ps(zmm2,xsinf(zmm3));
                              _mm512_store_ps(&cpowi[i+96],zmm5);
                              zmm6  = _mm512_load_ps(&xre[i+112]);
                              zmm7  = _mm512_mul_ps(zmm6,zmm6);
                              zmm8  = _mm512_load_ps(&xim[i+112]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_sqrt_ps(_mm512_add_ps(zmm7,zmm9));
                              zmmx  = _mm512_load_ps(&vn[i+112]);
                              zmm11 = xatanf(_mm512_div_ps(zmm8,zmm6));
                              zmm12 = xpowf(zmm10,zmmx);
                              zmm13 = _mm512_mul_ps(zmmx,zmm11);
                              zmm14 = _mm512_mul_ps(zmm12,xcosf(zmm13));
                              _mm512_store_ps(&cpowr[i+112],zmm14);
                              zmm15 = _mm512_mul_ps(zmm12,xsinf(zmm13));
                              _mm512_store_ps(&cpowi[i+112],zmm15);
                         }

                          for(; (i+95) < n; i += 96) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_load_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_load_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_load_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&cpowi[i+16],zmm19);
                              zmm20 = _mm512_load_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_load_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_load_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_store_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_store_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_load_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_load_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_load_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&cpowi[i+48],zmm7);
                              zmm8  = _mm512_load_ps(&xre[i+64]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_load_ps(&xim[i+64]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm9,zmm11));
                              zmmx  = _mm512_load_ps(&vn[i+64]);
                              zmm13 = xatanf(_mm512_div_ps(zmm10,zmm8));
                              zmm14 = xpowf(zmm12,zmmx);
                              zmm15 = _mm512_mul_ps(zmmx,zmm13);
                              zmm16 = _mm512_mul_ps(zmm14,xcosf(zmm15));
                              _mm512_store_ps(&cpowr[i+64],zmm16);
                              zmm17 = _mm512_mul_ps(zmm14,xsinf(zmm15));
                              _mm512_store_ps(&cpowi[i+64],zmm17);
                              zmm18  = _mm512_load_ps(&xre[i+80]);
                              zmm19  = _mm512_mul_ps(zmm18,zmm18);
                              zmm20 = _mm512_load_ps(&xim[i+80]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_sqrt_ps(_mm512_add_ps(zmm19,zmm21));
                              zmmx  = _mm512_load_ps(&vn[i+80]);
                              zmm23 = xatanf(_mm512_div_ps(zmm20,zmm18));
                              zmm24 = xpowf(zmm22,zmmx);
                              zmm25 = _mm512_mul_ps(zmmx,zmm23);
                              zmm26 = _mm512_mul_ps(zmm24,xcosf(zmm25));
                              _mm512_store_ps(&cpowr[i+80],zmm26);
                              zmm27 = _mm512_mul_ps(zmm24,xsinf(zmm25));
                              _mm512_store_ps(&cpowi[i+80],zmm27);
                         }

                          for(; (i+79) < n; i += 80) {
                               zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_load_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_load_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_load_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&cpowi[i+16],zmm19);
                              zmm20 = _mm512_load_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_load_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_load_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_store_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_store_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_load_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_load_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_load_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&cpowi[i+48],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+64]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_load_ps(&xim[i+64]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm9,zmm11));
                              zmmx  = _mm512_load_ps(&vn[i+64]);
                              zmm13 = xatanf(_mm512_div_ps(zmm10,zmm8));
                              zmm14 = xpowf(zmm12,zmmx);
                              zmm15 = _mm512_mul_ps(zmmx,zmm13);
                              zmm16 = _mm512_mul_ps(zmm14,xcosf(zmm15));
                              _mm512_store_ps(&cpowr[i+64],zmm16);
                              zmm17 = _mm512_mul_ps(zmm14,xsinf(zmm15));
                              _mm512_store_ps(&cpowi[i+64],zmm17);
                         }

                          for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_load_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_load_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_load_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&cpowi[i+16],zmm19);
                              zmm20 = _mm512_load_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_load_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_load_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_store_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_store_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_load_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_load_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_load_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&cpowi[i+48],zmm7);
                         }

                          for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_load_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_load_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_load_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&cpowi[i+16],zmm19);
                         }

                          for(; (i+15) < n; i += 15) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                         }

                           for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = z0*z0;
                              const float z2 = xim[i];
                              const float z3 = z2*z2;
                              const float z4 = ceph_sqrtf(z1+z3);
                              const float zx = vn[i];
                              const float z5 = ceph_atanf(z2/z0);
                              const float z6 = ceph_powf(z4,zx);
                              const float z7 = zx*z5;
                              const float z8 = z6*ceph_cosf(z7);
                              cpowr[i]       = z8;
                              const float z9 = z6*ceph_sinf(z7);
                              cpowi[i]       = z9;
                         }
                 }


                 
                   void gms::math::cpowv_zmm16r4_unroll_6x_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict vn,
                                                  float * __restrict cpowr,
                                                  float * __restrict cpowi,
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
                           vfloat zmmx;
                          int32_t i;

                          for(i = 0; (i+95) < n; i += 96) {
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+32], _MM_HINT_T0);
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_loadu_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_loadu_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_loadu_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_storeu_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_storeu_ps(&cpowi[i+16],zmm19);
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+64], _MM_HINT_T0);
                              zmm20 = _mm512_loadu_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_loadu_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_loadu_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_storeu_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_storeu_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_loadu_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_loadu_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_loadu_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_storeu_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_storeu_ps(&cpowi[i+48],zmm7);
                              _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+96], _MM_HINT_T0);
                              zmm8  = _mm512_loadu_ps(&xre[i+64]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_loadu_ps(&xim[i+64]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm9,zmm11));
                              zmmx  = _mm512_loadu_ps(&vn[i+64]);
                              zmm13 = xatanf(_mm512_div_ps(zmm10,zmm8));
                              zmm14 = xpowf(zmm12,zmmx);
                              zmm15 = _mm512_mul_ps(zmmx,zmm13);
                              zmm16 = _mm512_mul_ps(zmm14,xcosf(zmm15));
                              _mm512_storeu_ps(&cpowr[i+64],zmm16);
                              zmm17 = _mm512_mul_ps(zmm14,xsinf(zmm15));
                              _mm512_storeu_ps(&cpowi[i+64],zmm17);
                              zmm18  = _mm512_loadu_ps(&xre[i+80]);
                              zmm19  = _mm512_mul_ps(zmm18,zmm18);
                              zmm20 = _mm512_loadu_ps(&xim[i+80]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_sqrt_ps(_mm512_add_ps(zmm19,zmm21));
                              zmmx  = _mm512_loadu_ps(&vn[i+80]);
                              zmm23 = xatanf(_mm512_div_ps(zmm20,zmm18));
                              zmm24 = xpowf(zmm22,zmmx);
                              zmm25 = _mm512_mul_ps(zmmx,zmm23);
                              zmm26 = _mm512_mul_ps(zmm24,xcosf(zmm25));
                              _mm512_storeu_ps(&cpowr[i+80],zmm26);
                              zmm27 = _mm512_mul_ps(zmm24,xsinf(zmm25));
                              _mm512_storeu_ps(&cpowi[i+80],zmm27);
                             
                         }

                        
                          for(; (i+79) < n; i += 80) {
                               zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_loadu_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_loadu_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_loadu_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_storeu_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_storeu_ps(&cpowi[i+16],zmm19);
                              zmm20 = _mm512_loadu_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_loadu_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_loadu_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_storeu_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_storeu_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_loadu_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_loadu_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_loadu_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_storeu_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_storeu_ps(&cpowi[i+48],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+64]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_loadu_ps(&xim[i+64]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm9,zmm11));
                              zmmx  = _mm512_loadu_ps(&vn[i+64]);
                              zmm13 = xatanf(_mm512_div_ps(zmm10,zmm8));
                              zmm14 = xpowf(zmm12,zmmx);
                              zmm15 = _mm512_mul_ps(zmmx,zmm13);
                              zmm16 = _mm512_mul_ps(zmm14,xcosf(zmm15));
                              _mm512_storeu_ps(&cpowr[i+64],zmm16);
                              zmm17 = _mm512_mul_ps(zmm14,xsinf(zmm15));
                              _mm512_storeu_ps(&cpowi[i+64],zmm17);
                         }

                          for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_loadu_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_loadu_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_loadu_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_storeu_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_storeu_ps(&cpowi[i+16],zmm19);
                              zmm20 = _mm512_loadu_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_loadu_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_loadu_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_storeu_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_storeu_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_loadu_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_loadu_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_loadu_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_storeu_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_storeu_ps(&cpowi[i+48],zmm7);
                         }

                          for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_loadu_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_loadu_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_loadu_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_storeu_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_storeu_ps(&cpowi[i+16],zmm19);
                         }

                          for(; (i+15) < n; i += 15) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                         }

                           for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = z0*z0;
                              const float z2 = xim[i];
                              const float z3 = z2*z2;
                              const float z4 = ceph_sqrtf(z1+z3);
                              const float zx = vn[i];
                              const float z5 = ceph_atanf(z2/z0);
                              const float z6 = ceph_powf(z4,zx);
                              const float z7 = zx*z5;
                              const float z8 = z6*ceph_cosf(z7);
                              cpowr[i]       = z8;
                              const float z9 = z6*ceph_sinf(z7);
                              cpowi[i]       = z9;
                         }
                 }


             
                   void gms::math::cpowv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) vn,
                                                  float * __restrict __ATTR_ALIGN__(64) cpowr,
                                                  float * __restrict __ATTR_ALIGN__(64) cpowi,
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
                           vfloat zmmx;
                          int32_t i;

                          for(i = 0; (i+95) < n; i += 96) {
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+32], _MM_HINT_T0);
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_load_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_load_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_load_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&cpowi[i+16],zmm19);
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+64], _MM_HINT_T0);
                              zmm20 = _mm512_load_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_load_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_load_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_store_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_store_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_load_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_load_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_load_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&cpowi[i+48],zmm7);
                              _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+96], _MM_HINT_T0);
                              zmm8  = _mm512_load_ps(&xre[i+64]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_load_ps(&xim[i+64]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm9,zmm11));
                              zmmx  = _mm512_load_ps(&vn[i+64]);
                              zmm13 = xatanf(_mm512_div_ps(zmm10,zmm8));
                              zmm14 = xpowf(zmm12,zmmx);
                              zmm15 = _mm512_mul_ps(zmmx,zmm13);
                              zmm16 = _mm512_mul_ps(zmm14,xcosf(zmm15));
                              _mm512_store_ps(&cpowr[i+64],zmm16);
                              zmm17 = _mm512_mul_ps(zmm14,xsinf(zmm15));
                              _mm512_store_ps(&cpowi[i+64],zmm17);
                              zmm18  = _mm512_load_ps(&xre[i+80]);
                              zmm19  = _mm512_mul_ps(zmm18,zmm18);
                              zmm20 = _mm512_load_ps(&xim[i+80]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_sqrt_ps(_mm512_add_ps(zmm19,zmm21));
                              zmmx  = _mm512_load_ps(&vn[i+80]);
                              zmm23 = xatanf(_mm512_div_ps(zmm20,zmm18));
                              zmm24 = xpowf(zmm22,zmmx);
                              zmm25 = _mm512_mul_ps(zmmx,zmm23);
                              zmm26 = _mm512_mul_ps(zmm24,xcosf(zmm25));
                              _mm512_store_ps(&cpowr[i+80],zmm26);
                              zmm27 = _mm512_mul_ps(zmm24,xsinf(zmm25));
                              _mm512_store_ps(&cpowi[i+80],zmm27);
                             
                         }


                          for(; (i+79) < n; i += 80) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_load_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_load_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_load_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&cpowi[i+16],zmm19);
                              zmm20 = _mm512_load_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_load_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_load_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_store_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_store_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_load_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_load_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_load_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&cpowi[i+48],zmm7);
                              zmm8  = _mm512_loadu_ps(&xre[i+64]);
                              zmm9  = _mm512_mul_ps(zmm8,zmm8);
                              zmm10 = _mm512_load_ps(&xim[i+64]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_sqrt_ps(_mm512_add_ps(zmm9,zmm11));
                              zmmx  = _mm512_load_ps(&vn[i+64]);
                              zmm13 = xatanf(_mm512_div_ps(zmm10,zmm8));
                              zmm14 = xpowf(zmm12,zmmx);
                              zmm15 = _mm512_mul_ps(zmmx,zmm13);
                              zmm16 = _mm512_mul_ps(zmm14,xcosf(zmm15));
                              _mm512_store_ps(&cpowr[i+64],zmm16);
                              zmm17 = _mm512_mul_ps(zmm14,xsinf(zmm15));
                              _mm512_store_ps(&cpowi[i+64],zmm17);
                         }

                          for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_load_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_load_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_load_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&cpowi[i+16],zmm19);
                              zmm20 = _mm512_load_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_load_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_load_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_store_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_store_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_load_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_load_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_load_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&cpowi[i+48],zmm7);
                         }

                          for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_load_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_load_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_load_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&cpowi[i+16],zmm19);
                         }

                          for(; (i+15) < n; i += 15) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                         }

                           for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = z0*z0;
                              const float z2 = xim[i];
                              const float z3 = z2*z2;
                              const float z4 = ceph_sqrtf(z1+z3);
                              const float zx = vn[i];
                              const float z5 = ceph_atanf(z2/z0);
                              const float z6 = ceph_powf(z4,zx);
                              const float z7 = zx*z5;
                              const float z8 = z6*ceph_cosf(z7);
                              cpowr[i]       = z8;
                              const float z9 = z6*ceph_sinf(z7);
                              cpowi[i]       = z9;
                         }
                 }


                  
                   void gms::math::cpowv_zmm16r4_unroll_4x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                  const float * __restrict __ATTR_ALIGN__(64) xim,
                                                  const float * __restrict __ATTR_ALIGN__(64) vn,
                                                  float * __restrict __ATTR_ALIGN__(64) cpowr,
                                                  float * __restrict __ATTR_ALIGN__(64) cpowi,
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
                           vfloat zmmx;
                          int32_t i;

                          for(i = 0; (i+63) < n; i += 64) {
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+32], _MM_HINT_T0);
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_load_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_load_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_load_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&cpowi[i+16],zmm19);
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+64], _MM_HINT_T0);
                              zmm20 = _mm512_load_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_load_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_load_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_store_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_store_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_load_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_load_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_load_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_store_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_store_ps(&cpowi[i+48],zmm7);
                                                          
                         }

                           for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_load_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_load_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_load_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_store_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_store_ps(&cpowi[i+16],zmm19);
                         }

                          for(; (i+15) < n; i += 15) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_load_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_load_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_store_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_store_ps(&cpowi[i+0],zmm9);
                         }

                           for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = z0*z0;
                              const float z2 = xim[i];
                              const float z3 = z2*z2;
                              const float z4 = ceph_sqrtf(z1+z3);
                              const float zx = vn[i];
                              const float z5 = ceph_atanf(z2/z0);
                              const float z6 = ceph_powf(z4,zx);
                              const float z7 = zx*z5;
                              const float z8 = z6*ceph_cosf(z7);
                              cpowr[i]       = z8;
                              const float z9 = z6*ceph_sinf(z7);
                              cpowi[i]       = z9;
                         }
                 }


                  
                   void gms::math::cpowv_zmm16r4_unroll_4x_u(const float * __restrict xre,
                                                  const float * __restrict xim,
                                                  const float * __restrict vn,
                                                  float * __restrict cpowr,
                                                  float * __restrict cpowi,
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
                           vfloat zmmx;
                          int32_t i;

                          for(i = 0; (i+63) < n; i += 64) {
                              _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+32], _MM_HINT_T0);
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_loadu_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_loadu_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_loadu_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_storeu_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_storeu_ps(&cpowi[i+16],zmm19);
                              _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                              _mm_prefetch((const char*)&vn[i+64], _MM_HINT_T0);
                              zmm20 = _mm512_loadu_ps(&xre[i+32]);
                              zmm21 = _mm512_mul_ps(zmm20,zmm20);
                              zmm22 = _mm512_loadu_ps(&xim[i+32]);
                              zmm23 = _mm512_mul_ps(zmm22,zmm22);
                              zmm24 = _mm512_sqrt_ps(_mm512_add_ps(zmm21,zmm23));
                              zmmx  = _mm512_loadu_ps(&vn[i+32]);
                              zmm25 = xatanf(_mm512_div_ps(zmm22,zmm20));
                              zmm26 = xpowf(zmm24,zmmx);
                              zmm27 = _mm512_mul_ps(zmmx,zmm25);
                              zmm28 = _mm512_mul_ps(zmm26,xcosf(zmm27));
                              _mm512_storeu_ps(&cpowr[i+32],zmm28);
                              zmm29 = _mm512_mul_ps(zmm26,xsinf(zmm27));
                              _mm512_storeu_ps(&cpowi[i+32],zmm29);
                              zmm30 = _mm512_loadu_ps(&xre[i+48]);
                              zmm31 = _mm512_mul_ps(zmm30,zmm30);
                              zmm0  = _mm512_loadu_ps(&xim[i+48]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_sqrt_ps(_mm512_add_ps(zmm31,zmm1));
                              zmmx  = _mm512_loadu_ps(&vn[i+48]);
                              zmm3  = xatanf(_mm512_div_ps(zmm0,zmm30));
                              zmm4  = xpowf(zmm2,zmmx);
                              zmm5  = _mm512_mul_ps(zmmx,zmm3);
                              zmm6  = _mm512_mul_ps(zmm4,xcosf(zmm5));
                              _mm512_storeu_ps(&cpowr[i+48],zmm6);
                              zmm7  = _mm512_mul_ps(zmm4,xsinf(zmm5));
                              _mm512_storeu_ps(&cpowi[i+48],zmm7);
                              
                         }
                                                 
                          for(; (i+31) < n; i += 32) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                              zmm10 = _mm512_loadu_ps(&xre[i+16]);
                              zmm11 = _mm512_mul_ps(zmm10,zmm10);
                              zmm12 = _mm512_loadu_ps(&xim[i+16]);
                              zmm13 = _mm512_mul_ps(zmm12,zmm12);
                              zmm14 = _mm512_sqrt_ps(_mm512_add_ps(zmm11,zmm13));
                              zmmx  = _mm512_loadu_ps(&vn[i+16]);
                              zmm15 = xatanf(_mm512_div_ps(zmm12,zmm10));
                              zmm16 = xpowf(zmm14,zmmx);
                              zmm17 = _mm512_mul_ps(zmmx,zmm15);
                              zmm18 = _mm512_mul_ps(zmm16,xcosf(zmm17));
                              _mm512_storeu_ps(&cpowr[i+16],zmm18);
                              zmm19 = _mm512_mul_ps(zmm16,xsinf(zmm17));
                              _mm512_storeu_ps(&cpowi[i+16],zmm19);
                         }

                          for(; (i+15) < n; i += 15) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                              zmm1  = _mm512_mul_ps(zmm0,zmm0);
                              zmm2  = _mm512_loadu_ps(&xim[i+0]);
                              zmm3  = _mm512_mul_ps(zmm2,zmm2);
                              zmm4  = _mm512_sqrt_ps(_mm512_add_ps(zmm1,zmm3));
                              zmmx  = _mm512_loadu_ps(&vn[i+0]);
                              zmm5  = xatanf(_mm512_div_ps(zmm2,zmm0));
                              zmm6  = xpowf(zmm4,zmmx);
                              zmm7  = _mm512_mul_ps(zmmx,zmm5);
                              zmm8  = _mm512_mul_ps(zmm6,xcosf(zmm7));
                              _mm512_storeu_ps(&cpowr[i+0],zmm8);
                              zmm9  = _mm512_mul_ps(zmm6,xsinf(zmm7));
                              _mm512_storeu_ps(&cpowi[i+0],zmm9);
                         }

                           for(; (i+0) < n; i += 1) {
                              const float z0 = xre[i];
                              const float z1 = z0*z0;
                              const float z2 = xim[i];
                              const float z3 = z2*z2;
                              const float z4 = ceph_sqrtf(z1+z3);
                              const float zx = vn[i];
                              const float z5 = ceph_atanf(z2/z0);
                              const float z6 = ceph_powf(z4,zx);
                              const float z7 = zx*z5;
                              const float z8 = z6*ceph_cosf(z7);
                              cpowr[i]       = z8;
                              const float z9 = z6*ceph_sinf(z7);
                              cpowi[i]       = z9;
                         }
                 }



    
