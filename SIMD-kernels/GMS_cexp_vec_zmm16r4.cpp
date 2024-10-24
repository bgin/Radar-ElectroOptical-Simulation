



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
#include "GMS_cexp_vec_zmm16r4.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_cephes.h"




                   void gms::math::cexpv_zmm16r4_unroll_16x_u(const float * __restrict xre,
                                                   const float * __restrict xim,
                                                   float * __restrict cexpr,
                                                   float * __restrict cexpi,
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
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19);
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             zmm20  = _mm512_loadu_ps(&xre[i+64]);
                             zmm21  = _mm512_loadu_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_storeu_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_storeu_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_loadu_ps(&xre[i+80]);
                             zmm26  = _mm512_loadu_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_storeu_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_storeu_ps(&cexpi[i+80],zmm29);
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             zmm30  = _mm512_loadu_ps(&xre[i+96]);
                             zmm31  = _mm512_loadu_ps(&xim[i+96]);
                             zmm0  = xexpf(zmm30);
                             zmm1  = _mm512_mul_ps(zmm0,xcosf(zmm31));
                             _mm512_storeu_ps(&cexpr[i+96],zmm1);
                             zmm2  = _mm512_mul_ps(zmm0,xsinf(zmm31));
                             _mm512_storeu_ps(&cexpi[i+96],zmm2);
                             zmm3  = _mm512_loadu_ps(&xre[i+112]);
                             zmm4  = _mm512_loadu_ps(&xim[i+112]);
                             zmm5  = xexpf(zmm3);
                             zmm6  = _mm512_mul_ps(zmm5,xcosf(zmm4));
                             _mm512_storeu_ps(&cexpr[i+112],zmm6);
                             zmm7  = _mm512_mul_ps(zmm5,xsinf(zmm4));
                             _mm512_storeu_ps(&cexpi[i+112],zmm7);
                             _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                             zmm8  = _mm512_loadu_ps(&xre[i+128]);
                             zmm9  = _mm512_loadu_ps(&xim[i+128]);
                             zmm10 = xexpf(zmm8);
                             zmm11 = _mm512_mul_ps(zmm10,xcosf(zmm9));
                             _mm512_storeu_ps(&cexpr[i+128],zmm11);
                             zmm12 = _mm512_mul_ps(zmm10,xsinf(zmm9));
                             _mm512_storeu_ps(&cexpi[i+128],zmm12);
                             zmm13 = _mm512_loadu_ps(&xre[i+144]);
                             zmm14 = _mm512_loadu_ps(&xim[i+144]);
                             zmm15 = xexpf(zmm13);
                             zmm16 = _mm512_mul_ps(zmm15,xcosf(zmm14));
                             _mm512_storeu_ps(&cexpr[i+144],zmm26);
                             zmm17 = _mm512_mul_ps(zmm15,xsinf(zmm14));
                             _mm512_storeu_ps(&cexpi[i+144],zmm17);
                             _mm_prefetch((const char*)&xre[i+192],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+192],_MM_HINT_T0);
                             zmm18 = _mm512_loadu_ps(&xre[i+160]);
                             zmm19 = _mm512_loadu_ps(&xim[i+160]);
                             zmm20 = xexpf(zmm18);
                             zmm21 = _mm512_mul_ps(zmm18,xcosf(zmm19));
                             _mm512_storeu_ps(&cexpr[i+160],zmm21);
                             zmm22 = _mm512_mul_ps(zmm18,xsinf(zmm19));
                             _mm512_storeu_ps(&cexpi[i+160],zmm22);
                             zmm23 = _mm512_loadu_ps(&xre[i+176]);
                             zmm24 = _mm512_loadu_ps(&xim[i+176]);
                             zmm25 = xexpf(zmm23);
                             zmm26 = _mm512_mul_ps(zmm25,xcosf(zmm24));
                             _mm512_storeu_ps(&cexpr[i+176],zmm26);
                             zmm27 = _mm512_mul_ps(zmm25,xsinf(zmm24));
                             _mm512_storeu_ps(&cexpi[i+176],zmm27);
                             _mm_prefetch((const char*)&xre[i+224],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+224],_MM_HINT_T0);
                             zmm28 = _mm512_loadu_ps(&xre[i+192]);
                             zmm29 = _mm512_loadu_ps(&xim[i+192]);
                             zmm30 = xexpf(zmm28);
                             zmm31 = _mm512_mul_ps(zmm28,xcosf(zmm29));
                             _mm512_storeu_ps(&cexpr[i+192],zmm31);
                             zmm0 = _mm512_mul_ps(zmm28,xsinf(zmm29));
                             _mm512_storeu_ps(&cexpi[i+192],zmm0);
                             zmm1 = _mm512_loadu_ps(&xre[i+208]);
                             zmm2 = _mm512_loadu_ps(&xim[i+208]);
                             zmm3 = xexpf(zmm1);
                             zmm4 = _mm512_mul_ps(zmm3,xcosf(zmm2));
                             _mm512_storeu_ps(&cexpr[i+208],zmm4);
                             zmm5 = _mm512_mul_ps(zmm3,xsinf(zmm2));
                             _mm512_storeu_ps(&cexpi[i+208],zmm5);
                             _mm_prefetch((const char*)&xre[i+256],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+256],_MM_HINT_T0);
                             zmm6 = _mm512_loadu_ps(&xre[i+224]);
                             zmm7 = _mm512_loadu_ps(&xim[i+224]);
                             zmm8 = xexpf(zmm6);
                             zmm9 = _mm512_mul_ps(zmm8,xcosf(zmm7));
                             _mm512_storeu_ps(&cexpr[i+224],zmm9);
                             zmm10= _mm512_mul_ps(zmm8,xsinf(zmm7));
                             _mm512_storeu_ps(&cexpi[i+224],zmm10);
                             zmm11 = _mm512_loadu_ps(&xre[i+240]);
                             zmm12 = _mm512_loadu_ps(&xim[i+240]);
                             zmm13 = xexpf(zmm11);
                             zmm14 = _mm512_mul_ps(zmm13,xcosf(zmm12));
                             _mm512_storeu_ps(&cexpr[i+240],zmm14);
                             zmm15 = _mm512_mul_ps(zmm13,xsinf(zmm12));
                             _mm512_storeu_ps(&cexpi[i+240],zmm15);
                        }

                         for(; (i+191) < n; i += 192) {
                                zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_loadu_ps(&xre[i+64]);
                             zmm21  = _mm512_loadu_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_storeu_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_storeu_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_loadu_ps(&xre[i+80]);
                             zmm26  = _mm512_loadu_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_storeu_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_storeu_ps(&cexpi[i+80],zmm29);
                             zmm30  = _mm512_loadu_ps(&xre[i+96]);
                             zmm31  = _mm512_loadu_ps(&xim[i+96]);
                             zmm0  = xexpf(zmm30);
                             zmm1  = _mm512_mul_ps(zmm0,xcosf(zmm31));
                             _mm512_storeu_ps(&cexpr[i+96],zmm1);
                             zmm2  = _mm512_mul_ps(zmm0,xsinf(zmm31));
                             _mm512_storeu_ps(&cexpi[i+96],zmm2);
                             zmm3  = _mm512_loadu_ps(&xre[i+112]);
                             zmm4  = _mm512_loadu_ps(&xim[i+112]);
                             zmm5  = xexpf(zmm3);
                             zmm6  = _mm512_mul_ps(zmm5,xcosf(zmm4));
                             _mm512_storeu_ps(&cexpr[i+112],zmm6);
                             zmm7  = _mm512_mul_ps(zmm5,xsinf(zmm4));
                             _mm512_storeu_ps(&cexpi[i+112],zmm7);
                             zmm8  = _mm512_loadu_ps(&xre[i+128]);
                             zmm9  = _mm512_loadu_ps(&xim[i+128]);
                             zmm10 = xexpf(zmm8);
                             zmm11 = _mm512_mul_ps(zmm10,xcosf(zmm9));
                             _mm512_storeu_ps(&cexpr[i+128],zmm11);
                             zmm12 = _mm512_mul_ps(zmm10,xsinf(zmm9));
                             _mm512_storeu_ps(&cexpi[i+128],zmm12);
                             zmm13 = _mm512_loadu_ps(&xre[i+144]);
                             zmm14 = _mm512_loadu_ps(&xim[i+144]);
                             zmm15 = xexpf(zmm13);
                             zmm16 = _mm512_mul_ps(zmm15,xcosf(zmm14));
                             _mm512_storeu_ps(&cexpr[i+144],zmm26);
                             zmm17 = _mm512_mul_ps(zmm15,xsinf(zmm14));
                             _mm512_storeu_ps(&cexpi[i+144],zmm17);
                             zmm18 = _mm512_loadu_ps(&xre[i+160]);
                             zmm19 = _mm512_loadu_ps(&xim[i+160]);
                             zmm20 = xexpf(zmm18);
                             zmm21 = _mm512_mul_ps(zmm18,xcosf(zmm19));
                             _mm512_storeu_ps(&cexpr[i+160],zmm21);
                             zmm22 = _mm512_mul_ps(zmm18,xsinf(zmm19));
                             _mm512_storeu_ps(&cexpi[i+160],zmm22);
                             zmm23 = _mm512_loadu_ps(&xre[i+176]);
                             zmm24 = _mm512_loadu_ps(&xim[i+176]);
                             zmm25 = xexpf(zmm23);
                             zmm26 = _mm512_mul_ps(zmm25,xcosf(zmm24));
                             _mm512_storeu_ps(&cexpr[i+176],zmm26);
                             zmm27 = _mm512_mul_ps(zmm25,xsinf(zmm24));
                             _mm512_storeu_ps(&cexpi[i+176],zmm27);
                        }

                         for(; (i+127) < n; i += 128) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_loadu_ps(&xre[i+64]);
                             zmm21  = _mm512_loadu_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_storeu_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_storeu_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_loadu_ps(&xre[i+80]);
                             zmm26  = _mm512_loadu_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_storeu_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_storeu_ps(&cexpi[i+80],zmm29);
                             zmm30  = _mm512_loadu_ps(&xre[i+96]);
                             zmm31  = _mm512_loadu_ps(&xim[i+96]);
                             zmm0  = xexpf(zmm30);
                             zmm1  = _mm512_mul_ps(zmm0,xcosf(zmm31));
                             _mm512_storeu_ps(&cexpr[i+96],zmm1);
                             zmm2  = _mm512_mul_ps(zmm0,xsinf(zmm31));
                             _mm512_storeu_ps(&cexpi[i+96],zmm2);
                             zmm3  = _mm512_loadu_ps(&xre[i+112]);
                             zmm4  = _mm512_loadu_ps(&xim[i+112]);
                             zmm5  = xexpf(zmm3);
                             zmm6  = _mm512_mul_ps(zmm5,xcosf(zmm4));
                             _mm512_storeu_ps(&cexpr[i+112],zmm6);
                             zmm7  = _mm512_mul_ps(zmm5,xsinf(zmm4));
                             _mm512_storeu_ps(&cexpi[i+112],zmm7);
                        }

                         for(; (i+79) < n; i += 80) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_loadu_ps(&xre[i+64]);
                             zmm21  = _mm512_loadu_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_storeu_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_storeu_ps(&cexpi[i+64],zmm24);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19); 
                        }

                         for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                             const float re = xre[i];
                             const float im = xim[i];
                             const float x  = ceph_expf(re);
                             const float y  = x * ceph_sinf(im);
                             cexpr[i]       = y;
                             const float z  = x * ceph_cosf(im);
                             cexpi[i]       = z; 
                        }
                }


                  
                   void gms::math::cexpv_zmm16r4_unroll_16x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                   const float * __restrict __ATTR_ALIGN__(64) xim,
                                                   float * __restrict __ATTR_ALIGN__(64) cexpr,
                                                   float * __restrict __ATTR_ALIGN__(64) cexpi,
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
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19);
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             zmm20  = _mm512_load_ps(&xre[i+64]);
                             zmm21  = _mm512_load_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_store_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_store_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_load_ps(&xre[i+80]);
                             zmm26  = _mm512_load_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_store_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_store_ps(&cexpi[i+80],zmm29);
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             zmm30  = _mm512_load_ps(&xre[i+96]);
                             zmm31  = _mm512_load_ps(&xim[i+96]);
                             zmm0  = xexpf(zmm30);
                             zmm1  = _mm512_mul_ps(zmm0,xcosf(zmm31));
                             _mm512_store_ps(&cexpr[i+96],zmm1);
                             zmm2  = _mm512_mul_ps(zmm0,xsinf(zmm31));
                             _mm512_store_ps(&cexpi[i+96],zmm2);
                             zmm3  = _mm512_load_ps(&xre[i+112]);
                             zmm4  = _mm512_load_ps(&xim[i+112]);
                             zmm5  = xexpf(zmm3);
                             zmm6  = _mm512_mul_ps(zmm5,xcosf(zmm4));
                             _mm512_store_ps(&cexpr[i+112],zmm6);
                             zmm7  = _mm512_mul_ps(zmm5,xsinf(zmm4));
                             _mm512_store_ps(&cexpi[i+112],zmm7);
                             _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                             zmm8  = _mm512_load_ps(&xre[i+128]);
                             zmm9  = _mm512_load_ps(&xim[i+128]);
                             zmm10 = xexpf(zmm8);
                             zmm11 = _mm512_mul_ps(zmm10,xcosf(zmm9));
                             _mm512_store_ps(&cexpr[i+128],zmm11);
                             zmm12 = _mm512_mul_ps(zmm10,xsinf(zmm9));
                             _mm512_store_ps(&cexpi[i+128],zmm12);
                             zmm13 = _mm512_load_ps(&xre[i+144]);
                             zmm14 = _mm512_load_ps(&xim[i+144]);
                             zmm15 = xexpf(zmm13);
                             zmm16 = _mm512_mul_ps(zmm15,xcosf(zmm14));
                             _mm512_store_ps(&cexpr[i+144],zmm26);
                             zmm17 = _mm512_mul_ps(zmm15,xsinf(zmm14));
                             _mm512_store_ps(&cexpi[i+144],zmm17);
                             _mm_prefetch((const char*)&xre[i+192],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+192],_MM_HINT_T0);
                             zmm18 = _mm512_load_ps(&xre[i+160]);
                             zmm19 = _mm512_load_ps(&xim[i+160]);
                             zmm20 = xexpf(zmm18);
                             zmm21 = _mm512_mul_ps(zmm18,xcosf(zmm19));
                             _mm512_store_ps(&cexpr[i+160],zmm21);
                             zmm22 = _mm512_mul_ps(zmm18,xsinf(zmm19));
                             _mm512_store_ps(&cexpi[i+160],zmm22);
                             zmm23 = _mm512_load_ps(&xre[i+176]);
                             zmm24 = _mm512_load_ps(&xim[i+176]);
                             zmm25 = xexpf(zmm23);
                             zmm26 = _mm512_mul_ps(zmm25,xcosf(zmm24));
                             _mm512_store_ps(&cexpr[i+176],zmm26);
                             zmm27 = _mm512_mul_ps(zmm25,xsinf(zmm24));
                             _mm512_store_ps(&cexpi[i+176],zmm27);
                             _mm_prefetch((const char*)&xre[i+224],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+224],_MM_HINT_T0);
                             zmm28 = _mm512_load_ps(&xre[i+192]);
                             zmm29 = _mm512_load_ps(&xim[i+192]);
                             zmm30 = xexpf(zmm28);
                             zmm31 = _mm512_mul_ps(zmm28,xcosf(zmm29));
                             _mm512_store_ps(&cexpr[i+192],zmm31);
                             zmm0 = _mm512_mul_ps(zmm28,xsinf(zmm29));
                             _mm512_store_ps(&cexpi[i+192],zmm0);
                             zmm1 = _mm512_load_ps(&xre[i+208]);
                             zmm2 = _mm512_load_ps(&xim[i+208]);
                             zmm3 = xexpf(zmm1);
                             zmm4 = _mm512_mul_ps(zmm3,xcosf(zmm2));
                             _mm512_store_ps(&cexpr[i+208],zmm4);
                             zmm5 = _mm512_mul_ps(zmm3,xsinf(zmm2));
                             _mm512_store_ps(&cexpi[i+208],zmm5);
                             _mm_prefetch((const char*)&xre[i+256],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+256],_MM_HINT_T0);
                             zmm6 = _mm512_load_ps(&xre[i+224]);
                             zmm7 = _mm512_load_ps(&xim[i+224]);
                             zmm8 = xexpf(zmm6);
                             zmm9 = _mm512_mul_ps(zmm8,xcosf(zmm7));
                             _mm512_store_ps(&cexpr[i+224],zmm9);
                             zmm10= _mm512_mul_ps(zmm8,xsinf(zmm7));
                             _mm512_store_ps(&cexpi[i+224],zmm10);
                             zmm11 = _mm512_load_ps(&xre[i+240]);
                             zmm12 = _mm512_load_ps(&xim[i+240]);
                             zmm13 = xexpf(zmm11);
                             zmm14 = _mm512_mul_ps(zmm13,xcosf(zmm12));
                             _mm512_store_ps(&cexpr[i+240],zmm14);
                             zmm15 = _mm512_mul_ps(zmm13,xsinf(zmm12));
                             _mm512_store_ps(&cexpi[i+240],zmm15);
                        }

                         for(; (i+191) < n; i += 192) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_load_ps(&xre[i+64]);
                             zmm21  = _mm512_load_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_store_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_store_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_load_ps(&xre[i+80]);
                             zmm26  = _mm512_load_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_store_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_store_ps(&cexpi[i+80],zmm29);
                             zmm30  = _mm512_load_ps(&xre[i+96]);
                             zmm31  = _mm512_load_ps(&xim[i+96]);
                             zmm0  = xexpf(zmm30);
                             zmm1  = _mm512_mul_ps(zmm0,xcosf(zmm31));
                             _mm512_store_ps(&cexpr[i+96],zmm1);
                             zmm2  = _mm512_mul_ps(zmm0,xsinf(zmm31));
                             _mm512_store_ps(&cexpi[i+96],zmm2);
                             zmm3  = _mm512_load_ps(&xre[i+112]);
                             zmm4  = _mm512_load_ps(&xim[i+112]);
                             zmm5  = xexpf(zmm3);
                             zmm6  = _mm512_mul_ps(zmm5,xcosf(zmm4));
                             _mm512_store_ps(&cexpr[i+112],zmm6);
                             zmm7  = _mm512_mul_ps(zmm5,xsinf(zmm4));
                             _mm512_store_ps(&cexpi[i+112],zmm7);
                             zmm8  = _mm512_load_ps(&xre[i+128]);
                             zmm9  = _mm512_load_ps(&xim[i+128]);
                             zmm10 = xexpf(zmm8);
                             zmm11 = _mm512_mul_ps(zmm10,xcosf(zmm9));
                             _mm512_store_ps(&cexpr[i+128],zmm11);
                             zmm12 = _mm512_mul_ps(zmm10,xsinf(zmm9));
                             _mm512_store_ps(&cexpi[i+128],zmm12);
                             zmm13 = _mm512_load_ps(&xre[i+144]);
                             zmm14 = _mm512_load_ps(&xim[i+144]);
                             zmm15 = xexpf(zmm13);
                             zmm16 = _mm512_mul_ps(zmm15,xcosf(zmm14));
                             _mm512_store_ps(&cexpr[i+144],zmm26);
                             zmm17 = _mm512_mul_ps(zmm15,xsinf(zmm14));
                             _mm512_store_ps(&cexpi[i+144],zmm17);
                             zmm18 = _mm512_load_ps(&xre[i+160]);
                             zmm19 = _mm512_load_ps(&xim[i+160]);
                             zmm20 = xexpf(zmm18);
                             zmm21 = _mm512_mul_ps(zmm18,xcosf(zmm19));
                             _mm512_store_ps(&cexpr[i+160],zmm21);
                             zmm22 = _mm512_mul_ps(zmm18,xsinf(zmm19));
                             _mm512_store_ps(&cexpi[i+160],zmm22);
                             zmm23 = _mm512_load_ps(&xre[i+176]);
                             zmm24 = _mm512_load_ps(&xim[i+176]);
                             zmm25 = xexpf(zmm23);
                             zmm26 = _mm512_mul_ps(zmm25,xcosf(zmm24));
                             _mm512_store_ps(&cexpr[i+176],zmm26);
                             zmm27 = _mm512_mul_ps(zmm25,xsinf(zmm24));
                             _mm512_store_ps(&cexpi[i+176],zmm27);
                        }

                         for(; (i+127) < n; i += 128) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_load_ps(&xre[i+64]);
                             zmm21  = _mm512_load_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_store_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_store_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_load_ps(&xre[i+80]);
                             zmm26  = _mm512_load_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_store_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_store_ps(&cexpi[i+80],zmm29);
                             zmm30  = _mm512_load_ps(&xre[i+96]);
                             zmm31  = _mm512_load_ps(&xim[i+96]);
                             zmm0  = xexpf(zmm30);
                             zmm1  = _mm512_mul_ps(zmm0,xcosf(zmm31));
                             _mm512_store_ps(&cexpr[i+96],zmm1);
                             zmm2  = _mm512_mul_ps(zmm0,xsinf(zmm31));
                             _mm512_store_ps(&cexpi[i+96],zmm2);
                             zmm3  = _mm512_load_ps(&xre[i+112]);
                             zmm4  = _mm512_load_ps(&xim[i+112]);
                             zmm5  = xexpf(zmm3);
                             zmm6  = _mm512_mul_ps(zmm5,xcosf(zmm4));
                             _mm512_store_ps(&cexpr[i+112],zmm6);
                             zmm7  = _mm512_mul_ps(zmm5,xsinf(zmm4));
                             _mm512_store_ps(&cexpi[i+112],zmm7);
                        }

                         for(; (i+79) < n; i += 80) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_load_ps(&xre[i+64]);
                             zmm21  = _mm512_load_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mu_ps(zmm22,xcosf(zmm21));
                             _mm512_store_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_store_ps(&cexpi[i+64],zmm24);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19); 
                        }

                         for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                             const float re = xre[i];
                             const float im = xim[i];
                             const float x  = ceph_expf(re);
                             const float y  = x * ceph_sinf(im);
                             cexpr[i]       = y;
                             const float z  = x * ceph_cosf(im);
                             cexpi[i]       = z; 
                        }
                }


                 
                 
                   void gms::math::cexpv_zmm16r4_unroll_10x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                   const float * __restrict __ATTR_ALIGN__(64) xim,
                                                   float * __restrict __ATTR_ALIGN__(64) cexpr,
                                                   float * __restrict __ATTR_ALIGN__(64) cexpi,
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

                         for(i = 0; (i+159) < n; i += 160) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19);
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             zmm20  = _mm512_load_ps(&xre[i+64]);
                             zmm21  = _mm512_load_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_store_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_store_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_load_ps(&xre[i+80]);
                             zmm26  = _mm512_load_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_store_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_store_ps(&cexpi[i+80],zmm29);
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             zmm30  = _mm512_load_ps(&xre[i+96]);
                             zmm31  = _mm512_load_ps(&xim[i+96]);
                             zmm0  = xexpf(zmm30);
                             zmm1  = _mm512_mul_ps(zmm0,xcosf(zmm31));
                             _mm512_store_ps(&cexpr[i+96],zmm1);
                             zmm2  = _mm512_mul_ps(zmm0,xsinf(zmm31));
                             _mm512_store_ps(&cexpi[i+96],zmm2);
                             zmm3  = _mm512_load_ps(&xre[i+112]);
                             zmm4  = _mm512_load_ps(&xim[i+112]);
                             zmm5  = xexpf(zmm3);
                             zmm6  = _mm512_mul_ps(zmm5,xcosf(zmm4));
                             _mm512_store_ps(&cexpr[i+112],zmm6);
                             zmm7  = _mm512_mul_ps(zmm5,xsinf(zmm4));
                             _mm512_store_ps(&cexpi[i+112],zmm7);
                             _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                             zmm8  = _mm512_load_ps(&xre[i+128]);
                             zmm9  = _mm512_load_ps(&xim[i+128]);
                             zmm10 = xexpf(zmm8);
                             zmm11 = _mm512_mul_ps(zmm10,xcosf(zmm9));
                             _mm512_store_ps(&cexpr[i+128],zmm11);
                             zmm12 = _mm512_mul_ps(zmm10,xsinf(zmm9));
                             _mm512_store_ps(&cexpi[i+128],zmm12);
                             zmm13 = _mm512_load_ps(&xre[i+144]);
                             zmm14 = _mm512_load_ps(&xim[i+144]);
                             zmm15 = xexpf(zmm13);
                             zmm16 = _mm512_mul_ps(zmm15,xcosf(zmm14));
                             _mm512_store_ps(&cexpr[i+144],zmm26);
                             zmm17 = _mm512_mul_ps(zmm15,xsinf(zmm14));
                             _mm512_store_ps(&cexpi[i+144],zmm17);
                         }

                       
                         for(; (i+127) < n; i += 128) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_load_ps(&xre[i+64]);
                             zmm21  = _mm512_load_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_store_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_store_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_load_ps(&xre[i+80]);
                             zmm26  = _mm512_load_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_store_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_store_ps(&cexpi[i+80],zmm29);
                             zmm30  = _mm512_load_ps(&xre[i+96]);
                             zmm31  = _mm512_load_ps(&xim[i+96]);
                             zmm0  = xexpf(zmm30);
                             zmm1  = _mm512_mul_ps(zmm0,xcosf(zmm31));
                             _mm512_store_ps(&cexpr[i+96],zmm1);
                             zmm2  = _mm512_mul_ps(zmm0,xsinf(zmm31));
                             _mm512_store_ps(&cexpi[i+96],zmm2);
                             zmm3  = _mm512_load_ps(&xre[i+112]);
                             zmm4  = _mm512_load_ps(&xim[i+112]);
                             zmm5  = xexpf(zmm3);
                             zmm6  = _mm512_mul_ps(zmm5,xcosf(zmm4));
                             _mm512_store_ps(&cexpr[i+112],zmm6);
                             zmm7  = _mm512_mul_ps(zmm5,xsinf(zmm4));
                             _mm512_store_ps(&cexpi[i+112],zmm7);
                        }

                         for(; (i+79) < n; i += 80) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_load_ps(&xre[i+64]);
                             zmm21  = _mm512_load_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mu_ps(zmm22,xcosf(zmm21));
                             _mm512_store_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_store_ps(&cexpi[i+64],zmm24);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19); 
                        }

                         for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                             const float re = xre[i];
                             const float im = xim[i];
                             const float x  = ceph_expf(re);
                             const float y  = x * ceph_sinf(im);
                             cexpr[i]       = y;
                             const float z  = x * ceph_cosf(im);
                             cexpi[i]       = z; 
                        }
                }


                 
                   void gms::math::cexpv_zmm16r4_unroll_10x_u(const float * __restrict xre,
                                                   const float * __restrict xim,
                                                   float * __restrict cexpr,
                                                   float * __restrict cexpi,
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

                         for(i = 0; (i+159) < n; i += 160) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19);
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             zmm20  = _mm512_loadu_ps(&xre[i+64]);
                             zmm21  = _mm512_loadu_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_storeu_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_storeu_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_loadu_ps(&xre[i+80]);
                             zmm26  = _mm512_loadu_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_storeu_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_storeu_ps(&cexpi[i+80],zmm29);
                             _mm_prefetch((const char*)&xre[i+128],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+128],_MM_HINT_T0);
                             zmm30  = _mm512_loadu_ps(&xre[i+96]);
                             zmm31  = _mm512_loadu_ps(&xim[i+96]);
                             zmm0  = xexpf(zmm30);
                             zmm1  = _mm512_mul_ps(zmm0,xcosf(zmm31));
                             _mm512_storeu_ps(&cexpr[i+96],zmm1);
                             zmm2  = _mm512_mul_ps(zmm0,xsinf(zmm31));
                             _mm512_storeu_ps(&cexpi[i+96],zmm2);
                             zmm3  = _mm512_loadu_ps(&xre[i+112]);
                             zmm4  = _mm512_loadu_ps(&xim[i+112]);
                             zmm5  = xexpf(zmm3);
                             zmm6  = _mm512_mul_ps(zmm5,xcosf(zmm4));
                             _mm512_storeu_ps(&cexpr[i+112],zmm6);
                             zmm7  = _mm512_mul_ps(zmm5,xsinf(zmm4));
                             _mm512_storeu_ps(&cexpi[i+112],zmm7);
                             _mm_prefetch((const char*)&xre[i+160],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+160],_MM_HINT_T0);
                             zmm8  = _mm512_loadu_ps(&xre[i+128]);
                             zmm9  = _mm512_loadu_ps(&xim[i+128]);
                             zmm10 = xexpf(zmm8);
                             zmm11 = _mm512_mul_ps(zmm10,xcosf(zmm9));
                             _mm512_storeu_ps(&cexpr[i+128],zmm11);
                             zmm12 = _mm512_mul_ps(zmm10,xsinf(zmm9));
                             _mm512_storeu_ps(&cexpi[i+128],zmm12);
                             zmm13 = _mm512_loadu_ps(&xre[i+144]);
                             zmm14 = _mm512_loadu_ps(&xim[i+144]);
                             zmm15 = xexpf(zmm13);
                             zmm16 = _mm512_mul_ps(zmm15,xcosf(zmm14));
                             _mm512_storeu_ps(&cexpr[i+144],zmm26);
                             zmm17 = _mm512_mul_ps(zmm15,xsinf(zmm14));
                             _mm512_storeu_ps(&cexpi[i+144],zmm17);
                         }

                        
                         for(; (i+127) < n; i += 128) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_loadu_ps(&xre[i+64]);
                             zmm21  = _mm512_loadu_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_storeu_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_storeu_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_loadu_ps(&xre[i+80]);
                             zmm26  = _mm512_loadu_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_storeu_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_storeu_ps(&cexpi[i+80],zmm29);
                             zmm30  = _mm512_loadu_ps(&xre[i+96]);
                             zmm31  = _mm512_loadu_ps(&xim[i+96]);
                             zmm0  = xexpf(zmm30);
                             zmm1  = _mm512_mul_ps(zmm0,xcosf(zmm31));
                             _mm512_storeu_ps(&cexpr[i+96],zmm1);
                             zmm2  = _mm512_mul_ps(zmm0,xsinf(zmm31));
                             _mm512_storeu_ps(&cexpi[i+96],zmm2);
                             zmm3  = _mm512_loadu_ps(&xre[i+112]);
                             zmm4  = _mm512_loadu_ps(&xim[i+112]);
                             zmm5  = xexpf(zmm3);
                             zmm6  = _mm512_mul_ps(zmm5,xcosf(zmm4));
                             _mm512_storeu_ps(&cexpr[i+112],zmm6);
                             zmm7  = _mm512_mul_ps(zmm5,xsinf(zmm4));
                             _mm512_storeu_ps(&cexpi[i+112],zmm7);
                        }

                         for(; (i+79) < n; i += 80) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_loadu_ps(&xre[i+64]);
                             zmm21  = _mm512_loadu_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_storeu_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_storeu_ps(&cexpi[i+64],zmm24);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19); 
                        }

                         for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                             const float re = xre[i];
                             const float im = xim[i];
                             const float x  = ceph_expf(re);
                             const float y  = x * ceph_sinf(im);
                             cexpr[i]       = y;
                             const float z  = x * ceph_cosf(im);
                             cexpi[i]       = z; 
                        }
                }


                   
                   void gms::math::cexpv_zmm16r4_unroll_6x_u(const float * __restrict xre,
                                                   const float * __restrict xim,
                                                   float * __restrict cexpr,
                                                   float * __restrict cexpi,
                                                   const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           vfloat zmm0,zmm1,zmm2,zmm3;
                           vfloat zmm4,zmm5,zmm6,zmm7;
                           vfloat zmm8,zmm9,zmm10,zmm11;
                           vfloat zmm12,zmm13,zmm14,zmm15;
                           vfloat zmm16,zmm17,zmm18,zmm19;
                           vfloat zmm20,zmm21,zmm22,zmm23;
                           vfloat zmm24,zmm25,zmm26,zmm27;
                           vfloat zmm28,zmm29;
                          int32_t i;

                         for(i = 0; (i+95) < n; i += 96) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19);
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             zmm20  = _mm512_loadu_ps(&xre[i+64]);
                             zmm21  = _mm512_loadu_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_storeu_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_storeu_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_loadu_ps(&xre[i+80]);
                             zmm26  = _mm512_loadu_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_storeu_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_storeu_ps(&cexpi[i+80],zmm29);
                             
                         }

                        
                        
                         for(; (i+79) < n; i += 80) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_loadu_ps(&xre[i+64]);
                             zmm21  = _mm512_loadu_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_storeu_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_storeu_ps(&cexpi[i+64],zmm24);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_loadu_ps(&xre[i+32]);
                             zmm11  = _mm512_loadu_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_storeu_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_storeu_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_loadu_ps(&xre[i+48]);
                             zmm16  = _mm512_loadu_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_storeu_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_storeu_ps(&cexpi[i+48],zmm19); 
                        }

                         for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_loadu_ps(&xre[i+16]);
                             zmm6  = _mm512_loadu_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_storeu_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_storeu_ps(&cexpi[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_loadu_ps(&xre[i+0]);
                             zmm1  = _mm512_loadu_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_storeu_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_storeu_ps(&cexpi[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                             const float re = xre[i];
                             const float im = xim[i];
                             const float x  = ceph_expf(re);
                             const float y  = x * ceph_sinf(im);
                             cexpr[i]       = y;
                             const float z  = x * ceph_cosf(im);
                             cexpi[i]       = z; 
                        }
                }



                   void gms::math::cexpv_zmm16r4_unroll_6x_a(const float * __restrict __ATTR_ALIGN__(64) xre,
                                                   const float * __restrict __ATTR_ALIGN__(64) xim,
                                                   float * __restrict __ATTR_ALIGN__(64) cexpr,
                                                   float * __restrict __ATTR_ALIGN__(64) cexpi,
                                                   const int32_t n) {

                         if(__builtin_expect(0==n,0)) { return;}
                           vfloat zmm0,zmm1,zmm2,zmm3;
                           vfloat zmm4,zmm5,zmm6,zmm7;
                           vfloat zmm8,zmm9,zmm10,zmm11;
                           vfloat zmm12,zmm13,zmm14,zmm15;
                           vfloat zmm16,zmm17,zmm18,zmm19;
                           vfloat zmm20,zmm21,zmm22,zmm23;
                           vfloat zmm24,zmm25,zmm26,zmm27;
                           vfloat zmm28,zmm29;
                          int32_t i;

                         for(i = 0; (i+95) < n; i += 96) {
                             _mm_prefetch((const char*)&xre[i+32],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+32],_MM_HINT_T0);
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             _mm_prefetch((const char*)&xre[i+64],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+64],_MM_HINT_T0);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19);
                             _mm_prefetch((const char*)&xre[i+96],_MM_HINT_T0);
                             _mm_prefetch((const char*)&xim[i+96],_MM_HINT_T0);
                             zmm20  = _mm512_load_ps(&xre[i+64]);
                             zmm21  = _mm512_load_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mul_ps(zmm22,xcosf(zmm21));
                             _mm512_store_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_store_ps(&cexpi[i+64],zmm24);
                             zmm25  = _mm512_load_ps(&xre[i+80]);
                             zmm26  = _mm512_load_ps(&xim[i+80]);
                             zmm27  = xexpf(zmm25);
                             zmm28  = _mm512_mul_ps(zmm27,xcosf(zmm26));
                             _mm512_store_ps(&cexpr[i+80],zmm28);
                             zmm29  = _mm512_mul_ps(zmm27,xsinf(zmm26));
                             _mm512_store_ps(&cexpi[i+80],zmm29);
                             
                         }

                       
                        
                         for(; (i+79) < n; i += 80) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19);
                             zmm20  = _mm512_load_ps(&xre[i+64]);
                             zmm21  = _mm512_load_ps(&xim[i+64]);
                             zmm22  = xexpf(zmm20);
                             zmm23  = _mm512_mu_ps(zmm22,xcosf(zmm21));
                             _mm512_store_ps(&cexpr[i+64],zmm23);
                             zmm24  = _mm512_mul_ps(zmm22,xsinf(zmm21));
                             _mm512_store_ps(&cexpi[i+64],zmm24);
                        }

                         for(; (i+63) < n; i += 64) {
                              zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                             zmm10  = _mm512_load_ps(&xre[i+32]);
                             zmm11  = _mm512_load_ps(&xim[i+32]);
                             zmm12  = xexpf(zmm10);
                             zmm13  = _mm512_mul_ps(zmm12,xcosf(zmm11));
                             _mm512_store_ps(&cexpr[i+32],zmm13);
                             zmm14  = _mm512_mul_ps(zmm12,xsinf(zmm11));
                             _mm512_store_ps(&cexpi[i+32],zmm14);
                             zmm15  = _mm512_load_ps(&xre[i+48]);
                             zmm16  = _mm512_load_ps(&xim[i+48]);
                             zmm17  = xexpf(zmm15);
                             zmm18  = _mm512_mul_ps(zmm17,xcosf(zmm16));
                             _mm512_store_ps(&cexpr[i+48],zmm18);
                             zmm19  = _mm512_mul_ps(zmm17,xsinf(zmm16));
                             _mm512_store_ps(&cexpi[i+48],zmm19); 
                        }

                         for(; (i+31) < n; i += 32) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                             zmm5  = _mm512_load_ps(&xre[i+16]);
                             zmm6  = _mm512_load_ps(&xim[i+16]);
                             zmm7  = xexpf(zmm5);
                             zmm8  = _mm512_mul_ps(zmm7,xcosf(zmm6));
                             _mm512_store_ps(&cexpr[i+16],zmm8);
                             zmm9  = _mm512_mul_ps(zmm7,xsinf(zmm6));
                             _mm512_store_ps(&cexpi[i+16],zmm9);
                        }

                         for(; (i+15) < n; i += 16) {
                             zmm0  = _mm512_load_ps(&xre[i+0]);
                             zmm1  = _mm512_load_ps(&xim[i+0]);
                             zmm2  = xexpf(zmm0);
                             zmm3  = _mm512_mul_ps(zmm2,xcosf(zmm1));
                             _mm512_store_ps(&cexpr[i+0],zmm3);
                             zmm4  = _mm512_mul_ps(zmm2,xsinf(zmm1));
                             _mm512_store_ps(&cexpi[i+0],zmm4);
                        }

                         for(; (i+0) < n; i += 1) {
                             const float re = xre[i];
                             const float im = xim[i];
                             const float x  = ceph_expf(re);
                             const float y  = x * ceph_sinf(im);
                             cexpr[i]       = y;
                             const float z  = x * ceph_cosf(im);
                             cexpi[i]       = z; 
                        }
                }


                  


